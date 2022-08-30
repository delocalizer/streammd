"""
Read a SAM file from STDIN, mark duplicates in a single pass and stream
processed records to STDOUT.

Input must begin with a valid SAM header, followed by qname-grouped records.

Default log level is 'INFO' — set to something else with the LOG_LEVEL
environment variable.
"""
from importlib import metadata
from math import ceil
from multiprocessing import Process, Queue
from multiprocessing.managers import SharedMemoryManager
import argparse
import json
import logging
import os
import re
import sys
from .bloomfilter import BloomFilter

LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(os.environ.get('LOG_LEVEL', 'INFO'))
LOGGER.addHandler(logging.StreamHandler())

DEFAULT_FPRATE = 1e-6
DEFAULT_LOGINTERVAL = 1000000
DEFAULT_METRICS = 'streammd-metrics.json'
DEFAULT_NITEMS = 1e9
DEFAULT_NWORKERS = 4

ALIGNMENTS = 'ALIGNMENTS'
ALIGNMENTS_MARKED_DUPLICATE = 'ALIGNMENTS_MARKED_DUPLICATE'
READ_PAIRS = 'READ_PAIRS'
READ_PAIRS_MARKED_DUPLICATE = 'READ_PAIRS_MARKED_DUPLICATE'
READ_PAIR_DUPLICATE_FRACTION = 'READ_PAIR_DUPLICATE_FRACTION'
UNIQUE_ITEMS_APPROXIMATE = 'UNIQUE_ITEMS_APPROXIMATE'

MSG_ALIGNMENTS = 'alignments seen: %s'
MSG_ALIGNMENTS_MARKED_DUPLICATE = 'alignments marked duplicate: %s'
MSG_BATCHSIZE = 'running with batchsize=%s'
MSG_NOHEADER = 'no header lines detected'
MSG_QNAMEGRP = 'singleton %s: input is not paired reads or not qname grouped'
MSG_QSIZE = 'after %s qnames seen, approx input queue size is %s'
MSG_READ_PAIRS = 'read pairs seen: %s'
MSG_READ_PAIRS_MARKED_DUPLICATE = 'read pairs marked duplicate: %s'
MSG_READ_PAIR_DUPLICATE_FRACTION = 'read pair duplicate fraction: %0.4f'
MSG_UNIQUE_ITEMS_APPROXIMATE = 'approximate number of unique items: %s'
MSG_VERSION = 'streammd version %s'

CIGAR_CONSUMES_REF = {'M', 'D', 'N', '=', 'X'}
FLAG_UNMAPPED = 4
FLAG_REVERSE = 16
FLAG_SECONDARY = 256
FLAG_DUPLICATE = 1024
FLAG_SUPPLEMENTARY = 2048
PGID = f'{__name__.partition(".")[0]}'
PGTAG = 'PG:Z'
RE_CIGAR = re.compile(r'(?:(\d+)([MIDNSHPX=]))').findall
RE_LEADING_S = re.compile(r'^(\d+)S').search
RE_TRAILING_S = re.compile(r'(\d+)S$').search
SAM_OPTS_IDX = 11
SENTINEL = 'STOP'
# DEL sorts last in ascii
DEL = b'\x7F'.decode('ascii')
UNMAPPED = (DEL, -1, '')
VERSION = metadata.version(PGID)


def input_alnfile(infd, outfd, inq, nconsumers, batchsize=None):
    """
    Read records from a qname-grouped SAM file input stream and enqueue them
    in batches.

    Header lines are written directly to the output stream.

    Args:
        infd: Input stream file name or descriptor.
        outfd: Output stream file descriptor.
        inq: multiprocessing.Queue to put SAM records.
        nconsumers: Number of consumer processes.
        batchsize: Number of qnames per batch put into inq. The reader operates
            at fixed speed so bigger batch size => lower rate at which batches
            are added to the queue. Up to a point this reduces queue overhead
            (fewer put/get ops) but beyond that the rate of supply to the queue
            falls below the number required to keep all consumers sufficiently
            fed, and they end up blocking on get calls. If a batchsize value is
            not specified, the heuristic 400/nconsumers is used. Empirically
            this works quite well at least for 1 <= nconsumers <= 8.

    Returns:
        None
    """
    n_qname = 0
    batch = []
    qname_group = []
    qname_last = None
    header_done = None
    headlines = []
    batchsize = batchsize or ceil(400/nconsumers)
    LOGGER.info(MSG_BATCHSIZE, batchsize)
    with open(infd) as infh:
        for line in infh:
            if header_done:
                qname = line.partition('\t')[0]
                if qname == qname_last:
                    qname_group.append(line)
                else:
                    if len(qname_group) < 2:
                        raise ValueError(MSG_QNAMEGRP % qname)
                    n_qname += 1
                    if n_qname % DEFAULT_LOGINTERVAL == 0:
                        LOGGER.debug(MSG_QSIZE, n_qname, inq.qsize())
                    batch.append(qname_group)
                    if len(batch) == batchsize:
                        inq.put(batch)
                        batch = []
                    qname_last = qname
                    qname_group = [line]
            elif line.startswith('@'):
                headlines.append(line)
            else:
                if not headlines:
                    raise ValueError(MSG_NOHEADER)
                headlines.append(pgline(headlines[-1]))
                header_done = ''.join(headlines)
                os.write(outfd, header_done.encode('ascii'))
                qname_last = line.partition('\t')[0]
                qname_group = [line]
    batch.append(qname_group)
    inq.put(batch)
    for _ in range(nconsumers):
        inq.put(SENTINEL)


def markdup(record):
    """
    Mark the record as duplicate by updating in-place the FLAG, adding or
    updating the PG:Z: tag, and returning 1.

    If the record is unmapped, no update is done and 0 is returned.

    Args:
        record: list of SAM record str tokens.

    Retuns:
        number of duplicates marked: 1 if updated or 0 if not (unmapped read).
    """
    flag = int(record[1])
    # Replicate Picard MarkDuplicates behaviour: only an aligned read is
    # marked as duplicate.
    if flag & FLAG_UNMAPPED:
        return 0
    record[1] = str(flag + FLAG_DUPLICATE)
    pg_old, pg_new = None, f'{PGTAG}:{PGID}'
    for i, opt in enumerate(record[SAM_OPTS_IDX:], SAM_OPTS_IDX):
        if opt.startswith(PGTAG):
            pg_old = opt
            record[i] = pg_new
    if not pg_old:
        record.append(pg_new)
    return 1


def markdups(bfconfig, inq, outq, outfd):
    """
    Process SAM file records.

    Args:
        bfconfig: Bloom filter configuration dict.
        inq: multiprocessing.Queue to get batches of qname grouped SAM
             records.
        outq: multiprocessing.Queue to put results.
        outfd: output stream file descriptor.

    Results are added to the queue as:

        {
            READ_PAIRS: n_qname,                    # read pairs (qnames) seen
            READ_PAIRS_MARKED_DUPLICATE: n_rp_dup,  # read pairs marked dup
            ALIGNMENTS: n_align,                    # alignments seen
            ALIGNMENTS_MARKED_DUPLICATE: n_aln_dup  # alignments marked dup
        }
    """
    bf = BloomFilter.copy(bfconfig)
    n_qname, n_rp_dup, n_align, n_aln_dup = 0, 0, 0, 0
    while True:
        batch = inq.get()
        outlines = []
        if batch == SENTINEL:
            break
        for group in batch:
            n_qname += 1
            n_align += len(group)
            qnamegrp = [r.strip().split('\t') for r in group]
            if (ends := readends(qnamegrp)):
                ends_str = [f'{end[0]}_{end[1]}{end[2]}' for end in ends]
                if (
                        # one end mapped and is dupe
                        (ends[1] == UNMAPPED and not bf.add(ends_str[0])) or
                        # both ends mapped and frag is dupe
                        not bf.add(''.join(ends_str))
                ):
                    n_rp_dup += 1
                    for read in qnamegrp:
                        n_aln_dup += markdup(read)
            for read in qnamegrp:
                outlines.append('\t'.join(read))
        # Write the batch as a batch. In contrast to sys.stdout.write,
        # os.write is atomic so we don't have to care about locking or
        # using an output queue.
        out = '\n'.join(outlines) + '\n'
        os.write(outfd, out.encode('ascii'))
    outq.put({
        READ_PAIRS:                  n_qname,
        READ_PAIRS_MARKED_DUPLICATE: n_rp_dup,
        ALIGNMENTS:                  n_align,
        ALIGNMENTS_MARKED_DUPLICATE: n_aln_dup,
    })


def mem_calc(n, p):
    """
    Returns approximate memory requirement in GB for n items and target maximum
    false positive rate p.
    """
    m, _ = BloomFilter.optimal_m_k(n, p)
    return m / 8 / 1024 ** 3


def output_metrics(metrics, metfh):
    """
    Output metrics.

    Args:
        metrics: Dict of metrics.
        metfh: Open file handle for writing.
    """
    LOGGER.info(MSG_UNIQUE_ITEMS_APPROXIMATE,
                metrics[UNIQUE_ITEMS_APPROXIMATE])
    LOGGER.info(MSG_ALIGNMENTS, metrics[ALIGNMENTS])
    LOGGER.info(MSG_ALIGNMENTS_MARKED_DUPLICATE,
                metrics[ALIGNMENTS_MARKED_DUPLICATE])
    LOGGER.info(MSG_READ_PAIRS, metrics[READ_PAIRS])
    LOGGER.info(MSG_READ_PAIRS_MARKED_DUPLICATE,
                metrics[READ_PAIRS_MARKED_DUPLICATE])
    LOGGER.info(MSG_READ_PAIR_DUPLICATE_FRACTION,
                metrics[READ_PAIR_DUPLICATE_FRACTION])
    metfh.write(
        # kludge to output rounded floats
        # https://stackoverflow.com/a/29066406/6705037
        json.dumps(
            json.loads(
                json.dumps(metrics),
                parse_float=lambda x: round(float(x), 4)),
            indent=2,
            sort_keys=True))


def parse_cmdargs(args):
    """
    Returns: Parsed arguments
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',
                        default=0,
                        help='Input file (default=STDIN).')
    parser.add_argument('--output',
                        default=1,
                        help='Output file (default=STDOUT).')
    parser.add_argument('--metrics',
                        default=DEFAULT_METRICS,
                        help=('Output metrics file '
                              f'(default={DEFAULT_METRICS}).'))
    parser.add_argument('-n', '--n-items',
                        type=int,
                        default=DEFAULT_NITEMS,
                        help=('Expected maximum number of read pairs n '
                              f'(default={DEFAULT_NITEMS:.2E}).'))
    parser.add_argument('-p', '--fp-rate',
                        type=float,
                        default=DEFAULT_FPRATE,
                        help=('Target maximum false positive rate when n items '
                              f'are stored (default={DEFAULT_FPRATE:.2E}).'))
    parser.add_argument('--consumer-processes',
                        type=int,
                        default=DEFAULT_NWORKERS,
                        help=('Number of hashing processes '
                              f'(default={DEFAULT_NWORKERS}).'))
    parser.add_argument('--input-batch-size',
                        type=int,
                        help=('Specify the number of SAM records in each '
                              'batch for the input queue. If not specified, '
                              'the heuristic 400/nconsumers is used.'))
    parser.add_argument('--mem-calc',
                        type=float,
                        nargs=2,
                        metavar=('N_ITEMS', 'FP_RATE'),
                        help=('Print approximate memory requirement in GB '
                              'for n items and target maximum false positive '
                              'rate p.'))
    parser.add_argument('--version',
                        action='version',
                        version=VERSION)
    return parser.parse_args(args)


def pgline(last):
    """
    Return the @PG header line for data processed by this tool.

    Args:
        last: the last line of the header as read before processing.
    """
    tags = [
        f'ID:{PGID}',
        f'PN:{PGID}',
        f'CL:{" ".join(sys.argv)}',
        f'VN:{VERSION}'
    ]
    PP = None
    tkns = last.strip().split('\t')
    if tkns[0] == '@PG':
        prev = {tag: value for tag, value in
                (tkn.split(':', 1) for tkn in tkns[1:])}
        PP = prev.get('ID')
    if PP:
        tags.insert(2, f'PP:{PP}')
    return '\t'.join(['@PG'] + tags) + '\n'


def readends(qnamegrp):
    """
    Calculate ends of the fragment, accounting for soft-clipped bases.

    Args:
        qnamegrp: QNAME group of SAM records, each record supplied as a list
            of str — i.e. the result of calling .split(TAB) on a SAM text line.

    Returns:
        None if there are no aligned primary reads, otherwise a coord-sorted
        pair of ends:

            [(left_rname, left_ref_start, left_orientation),
                (right_rname, right_ref_end, right_orientation)]

        a single unmapped end always appears last with the value UNMAPPED.
    """
    ends = [UNMAPPED, UNMAPPED]
    idx = 0
    for read in qnamegrp:
        flag      = int(read[1])
        rname     =     read[2]
        ref_start = int(read[3])
        cigar     =     read[5]
        # use only the primary alignments for end calculation
        if flag & FLAG_SECONDARY or flag & FLAG_SUPPLEMENTARY:
            continue
        # unmapped
        if flag & FLAG_UNMAPPED:
            pass
        # forward
        elif not flag & FLAG_REVERSE:
            # Leading soft clips.
            leading_s = int(match[1]) if (match := RE_LEADING_S(cigar)) else 0
            ends[idx] = rname, ref_start - leading_s, 'F'
        # reverse
        else:
            # Trailing soft clips
            trailing_s = int(match[1]) if (match := RE_TRAILING_S(cigar)) else 0
            ref_end = ref_start
            for num, op in RE_CIGAR(cigar):
                ref_end += (int(num) if op in CIGAR_CONSUMES_REF else 0)
            ends[idx] = rname, ref_end + trailing_s, 'R'
        idx += 1
    assert idx == 2

    # Canonical ordering: l < r and UNMAPPED is always last by construction.
    ends.sort()
    return None if ends == [UNMAPPED, UNMAPPED] else ends


def main():
    """
    Run as CLI script.
    """
    args = parse_cmdargs(sys.argv[1:])
    if args.mem_calc:
        print(f'{mem_calc(*args.mem_calc):0.3f}GB')
        sys.exit(0)
    LOGGER.info(MSG_VERSION, VERSION)
    LOGGER.info(' '.join(sys.argv))
    nconsumers = args.consumer_processes
    inputbatchsize = args.input_batch_size
    inq = Queue(20 * nconsumers)
    outq = Queue(nconsumers)
    with (SharedMemoryManager() as smm,
          open(args.input) as infh,
          open(args.output, 'w') as outfh,
          open(args.metrics, 'w') as metfh):
        infd, outfd = infh.fileno(), outfh.fileno()
        reader = Process(target=input_alnfile,
                         args=(infd, outfd, inq, nconsumers, inputbatchsize))
        reader.start()
        bf = BloomFilter(smm, args.n_items, args.fp_rate)
        consumers = [
            Process(target=markdups, args=(bf.config, inq, outq, outfd))
            for _ in range(nconsumers)
        ]
        list(map(lambda x: x.start(), consumers))
        list(map(lambda x: x.join(), consumers))
        reader.join()
        counts = [outq.get() for _ in range(nconsumers)]
        metrics = {k: sum(d[k] for d in counts) for k in counts[0]}
        metrics[UNIQUE_ITEMS_APPROXIMATE] = bf.count()
        metrics[READ_PAIR_DUPLICATE_FRACTION] = (
                metrics[READ_PAIRS_MARKED_DUPLICATE]/metrics[READ_PAIRS])
        output_metrics(metrics, metfh)


if __name__ == '__main__':
    main()
