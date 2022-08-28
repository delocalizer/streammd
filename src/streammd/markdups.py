"""
Read a SAM file from STDIN, mark duplicates in a single pass and stream
processed records to STDOUT.

Input must begin with a valid SAM header, followed by qname-grouped records.

Default log level is 'INFO' — set to something else with the LOG_LEVEL
environment variable.
"""
from importlib import metadata
from itertools import repeat
from math import ceil
from multiprocessing import Pool, Process, Queue
from multiprocessing.managers import SharedMemoryManager
import argparse
import json
import logging
import os
import sys
from pysam import AlignmentHeader, AlignedSegment
from .bloomfilter import BloomFilter

LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(os.environ.get('LOG_LEVEL', 'INFO'))
LOGGER.addHandler(logging.StreamHandler())

DEFAULT_FPRATE = 1e-6
DEFAULT_NITEMS = 1e9
DEFAULT_NWORKERS = 4
DEFAULT_INQSIZE = 100
DEFAULT_METRICS = 'streammd-metrics.json'

ALIGNMENTS = 'ALIGNMENTS'
ALIGNMENTS_MARKED_DUPLICATE = 'ALIGNMENTS_MARKED_DUPLICATE'
READ_PAIRS = 'READ_PAIRS'
READ_PAIRS_MARKED_DUPLICATE = 'READ_PAIRS_MARKED_DUPLICATE'
READ_PAIR_DUPLICATE_FRACTION = 'READ_PAIR_DUPLICATE_FRACTION'
UNIQUE_ITEMS_APPROXIMATE = 'UNIQUE_ITEMS_APPROXIMATE'

MSG_BATCHSIZE = 'running with batchsize=%s'
MSG_ALIGNMENTS = 'alignments seen: %s'
MSG_ALIGNMENTS_MARKED_DUPLICATE = 'alignments marked duplicate: %s'
MSG_READ_PAIRS = 'read pairs seen: %s'
MSG_READ_PAIRS_MARKED_DUPLICATE = 'read pairs marked duplicate: %s'
MSG_READ_PAIR_DUPLICATE_FRACTION = 'read pair duplicate fraction: %0.4f'
MSG_UNIQUE_ITEMS_APPROXIMATE = 'approximate number of unique items: %s'

MSG_NOHEADER = 'no header lines detected'
MSG_QNAMEGRP = 'singleton %s: input is not paired reads or not qname grouped'
MSG_VERSION = 'streammd version %s'

PGID = f'{__name__.partition(".")[0]}'
SENTINEL = 'STOP'
# refID in SAM spec is int32 so first element is > any legal value.
UNMAPPED = (2**31, -1, '')
VERSION = metadata.version(PGID)


def input_alnfile(infd, outfd, headerq, inq, nconsumers, batchsize=None):
    """
    Read records from a qname-grouped SAM file input stream and enqueue them
    in batches.

    Header lines are written directly to the output stream and also to the
    header Value.

    Args:
        infd: Input stream file name or descriptor.
        outfd: Output stream file descriptor.
        headerq: multiprocessing.Queue to put header text.
        inq: multiprocessing.Queue to put SAM records.
        nconsumers: Number of consumer processes.
        batchsize: Number of qnames per batch put to inq. If not specified, the
                   heuristic 250/nconsumers is used.

    Returns:
        None
    """
    batch = []
    group = []
    groupid = None
    header_txt = None
    headlines = []
    batchsize = batchsize or ceil(250/nconsumers)
    LOGGER.info(MSG_BATCHSIZE, batchsize)
    with open(infd) as infh:
        for line in infh:
            if line.startswith('@'):
                headlines.append(line)
            else:
                if not header_txt:
                    if not headlines:
                        raise ValueError(MSG_NOHEADER)
                    headlines.append(pgline(headlines[-1]))
                    header_txt = ''.join(headlines)
                    os.write(outfd, header_txt.encode('ascii'))
                    headerq.put(header_txt)
                record = line.strip()
                qname = record.partition('\t')[0]
                if qname == groupid:
                    group.append(record)
                else:
                    if group:
                        if not len(group) > 1:
                            raise ValueError(MSG_QNAMEGRP % qname)
                        batch.append(group)
                        if len(batch) == batchsize:
                            inq.put(batch)
                            batch = []
                    groupid = qname
                    group = [record]
    batch.append(group)
    inq.put(batch)
    for _ in range(nconsumers):
        inq.put(SENTINEL)


def markdups(bfconfig, header, inq, outq, outfd):
    """
    Process SAM file records.

    Args:
        bfconfig: Bloom filter configuration dict.
        header: SAM file header as string.
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
    ah = AlignmentHeader.from_text(header)
    n_qname, n_rp_dup, n_align, n_aln_dup = 0, 0, 0, 0
    while True:
        batch = inq.get()
        if batch == SENTINEL:
            break
        for group in batch:
            n_qname += 1
            n_align += len(group)
            alignments = [AlignedSegment.fromstring(r, ah) for r in group]
            if not (ends := readends(alignments)):
                continue
            ends_str = [f'{end[0]}_{end[1]}{end[2]}' for end in ends]
            if ends[1] == UNMAPPED and not bf.add(ends_str[0]):
                n_rp_dup += 1
                for a in alignments:
                    # Replicate Picard MarkDuplicates behaviour: only the
                    # aligned read is marked as duplicate.
                    if a.is_mapped:
                        n_aln_dup += 1
                        a.flag += 1024
                        a.set_tag('PG', PGID, 'Z')
            elif not bf.add(''.join(ends_str)):
                n_rp_dup += 1
                for a in alignments:
                    n_aln_dup += 1
                    a.flag += 1024
                    a.set_tag('PG', PGID, 'Z')
            # Write the group as a group. In contrast to sys.stdout.write,
            # os.write is atomic so we don't have to care about locking or
            # using an output queue.
            out = '\n'.join(a.to_string() for a in alignments) + '\n'
            os.write(outfd, out.encode('ascii'))
    outq.put({
        READ_PAIRS: n_qname,
        READ_PAIRS_MARKED_DUPLICATE: n_rp_dup,
        ALIGNMENTS: n_align,
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
                              'batch for the input queue. If not specified '
                              'the heuristic 250/nconsumers is used.'))
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


def readends(alignments):
    """
    Calculate ends of the fragment, accounting for soft-clipped bases.

    Args:
        alignments: Qname group tuple of AlignedSegment instances.

    Returns:
        None if there are no aligned reads, otherwise a coordinate-sorted pair
        of ends:

            [(left_refid, left_pos, left_orientation),
                (right_refid, right_pos, right_orientation)]

        a single unmapped end always appears last with the value UNMAPPED.
    """
    r12 = [None, None]
    ends = [UNMAPPED, UNMAPPED]

    # Pick the primary alignments.
    for alignment in alignments:
        if not (alignment.is_secondary or alignment.is_supplementary):
            if alignment.is_read1:
                r12[0] = alignment
            elif alignment.is_read2:
                r12[1] = alignment

    # Bail if neither aligns.
    if all(r.is_unmapped for r in r12):
        return None

    # Calculate the ends.
    for i, r in enumerate(r12):
        if r.is_unmapped:
            pass
        elif r.is_forward:
            # Leading soft clips.
            front_s = r.cigar[0][1] if r.cigar[0][0] == 4 else 0
            ends[i] = r.reference_id, r.reference_start - front_s, 'F'
        elif r.is_reverse:
            # Trailing soft clips.
            back_s = r.cigar[-1][1] if r.cigar[-1][0] == 4 else 0
            ends[i] = r.reference_id, r.reference_end + back_s, 'R'

    # Canonical ordering: l < r and UNMAPPED is always last by construction.
    ends.sort()
    return ends


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
    headerq = Queue(1)
    inq = Queue(DEFAULT_INQSIZE)
    outq = Queue(nconsumers)
    with (SharedMemoryManager() as smm,
          open(args.input) as infh,
          open(args.output, 'w') as outfh,
          open(args.metrics, 'w') as metfh):
        infd, outfd = infh.fileno(), outfh.fileno()
        reader = Process(target=input_alnfile,
                         args=(infd, outfd, headerq, inq, nconsumers,
                               inputbatchsize))
        reader.start()
        header = headerq.get()
        bf = BloomFilter(smm, args.n_items, args.fp_rate)
        consumers = [                                                       
            Process(target=markdups, args=(bf.config, header, inq, outq, outfd))
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
