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
from typing import List, Literal, Optional, TextIO, Tuple, TypedDict
import argparse
import json
import logging
import os
import re
import sys
from .bloomfilter import BloomFilter, BloomFilterConfig

LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(os.environ.get('LOG_LEVEL', 'INFO'))
LOGGER.addHandler(logging.StreamHandler())

DEFAULT_FPRATE = 1e-6
DEFAULT_LOGINTERVAL = 1000000
DEFAULT_METRICS = 'streammd-metrics.json'
DEFAULT_NITEMS = 1e9
DEFAULT_NWORKERS = 4

MSG_ALIGNMENTS = 'alignments seen: %s'
MSG_ALIGNMENTS_MARKED_DUPLICATE = 'alignments marked duplicate: %s'
MSG_BATCHSIZE = 'running with batchsize=%s'
MSG_NOHEADER = 'no header lines detected'
MSG_NOTSINGLE = ('%s: expected 1 primary alignment, got %s. Input is not '
                 'single-end reads?')
MSG_NOTPAIRED = ('%s: expected 2 primary alignments, got %s. Input is not '
                 'paired-end reads or not qname grouped?')
MSG_QSIZE = 'after %s qnames seen, approx input queue size is %s'
MSG_TEMPLATES = 'templates seen: %s'
MSG_TEMPLATES_MARKED_DUPLICATE = 'templates marked duplicate: %s'
MSG_TEMPLATE_DUPLICATE_FRACTION = 'template duplicate fraction: %0.4f'
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
UNMAPPED = (DEL, 0, '')
VERSION = metadata.version(PGID)


class Metrics(TypedDict, total=False):
    """
    Metrics from duplicate marking processes.
    """
    ALIGNMENTS: int
    ALIGNMENTS_MARKED_DUPLICATE: int
    TEMPLATES: int
    TEMPLATES_MARKED_DUPLICATE: int
    TEMPLATE_DUPLICATE_FRACTION: float


def input_alnfile(infd: int|str,
                  outfd: int,
                  inq: Queue,
                  nconsumers: int,
                  batchsize:Optional[int]=None) -> None:
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
    # the single producer process is the bottleneck so we have some fussy
    # optimizations in here at the cost of readability — e.g. increment
    # batch_sz rather than calling len(batch).
    n_qname = 0
    batch = []
    batch_sz = 0
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
                    n_qname += 1
                    if n_qname % DEFAULT_LOGINTERVAL == 0:
                        LOGGER.debug(MSG_QSIZE, n_qname, inq.qsize())
                    batch.append(qname_group)
                    batch_sz += 1
                    if batch_sz == batchsize:
                        inq.put(batch)
                        batch = []
                        batch_sz = 0
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


def markdup(record: List[str]) -> int:
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
    record[1] = str(flag | FLAG_DUPLICATE)
    pg_old, pg_new = None, f'{PGTAG}:{PGID}'
    for i, opt in enumerate(record[SAM_OPTS_IDX:], SAM_OPTS_IDX):
        if opt.startswith(PGTAG):
            pg_old = opt
            record[i] = pg_new
    if not pg_old:
        record.append(pg_new)
    return 1


def markdups(bfconfig: BloomFilterConfig,
             inq: Queue,
             outq: Queue,
             outfd: int,
             reads_per_template: Literal[1, 2],
             strip_previous: bool=False) -> None:
    """
    Process SAM file records.

    Args:
        bfconfig: Bloom filter configuration dict.
        inq: multiprocessing.Queue to get batches of qname grouped SAM
             records.
        outq: multiprocessing.Queue to put results.
        outfd: output stream file descriptor.
        reads_per_template: 1 or 2.
        strip_previous: unset duplicate flag bit for any reads that have it
            set and are no longer considered duplicate (default=False). Not
            necessary unless records have previously been through a duplicate
            marking step, in which case it is strongly recommended for
            sensible results.

    Results are added to the queue as a Metrics instance:

        {
            TEMPLATES:                   n_tpl,     # templates (qnames) seen
            TEMPLATES_MARKED_DUPLICATE:  n_tpl_dup, # templates marked dup
            ALIGNMENTS:                  n_aln,     # alignments seen
            ALIGNMENTS_MARKED_DUPLICATE: n_aln_dup  # alignments marked dup
        }
    """
    bf = BloomFilter.copy(bfconfig)
    n_tpl, n_tpl_dup, n_aln, n_aln_dup = 0, 0, 0, 0
    while True:
        batch = inq.get()
        outlines = []
        if batch == SENTINEL:
            break
        for group in batch:
            n_tpl += 1
            n_aln += len(group)
            qnamegrp = [r.strip().split('\t') for r in group]
            ends = readends(qnamegrp)
            if len(ends) != reads_per_template:
                errmsg = (MSG_NOTSINGLE if reads_per_template == 1 else
                          MSG_NOTPAIRED)
                raise ValueError(errmsg % (qnamegrp[0][0], len(ends)))
            ends_str = ''.join([f'{end[0]}_{end[1]}{end[2]}' for end in ends])
            # sort order => if 1st is unmapped, all are
            if ends[0] == UNMAPPED:
                pass
            # ends already seen => dupe
            elif not bf.add(ends_str):
                n_tpl_dup += 1
                for read in qnamegrp:
                    n_aln_dup += markdup(read)
            elif strip_previous:
                for read in qnamegrp:
                    unmarkdup(read)
            for read in qnamegrp:
                outlines.append('\t'.join(read))
        # Write the batch as a batch. In contrast to sys.stdout.write,
        # os.write is atomic so we don't have to care about locking or
        # using an output queue.
        out = '\n'.join(outlines) + '\n'
        os.write(outfd, out.encode('ascii'))
    outq.put(Metrics({
        'TEMPLATES':                   n_tpl,
        'TEMPLATES_MARKED_DUPLICATE':  n_tpl_dup,
        'ALIGNMENTS':                  n_aln,
        'ALIGNMENTS_MARKED_DUPLICATE': n_aln_dup,
    }))


def output_metrics(metrics: List[Metrics], metfh: TextIO) -> None:
    """
    Output aggregate metrics.

    Args:
        metrics: list of duplicate marking metrics from output queue.
        metfh: Open file handle for writing.
    """
    agg = Metrics(
        ALIGNMENTS=sum(m['ALIGNMENTS'] for m in metrics),
        ALIGNMENTS_MARKED_DUPLICATE=sum(
            m['ALIGNMENTS_MARKED_DUPLICATE'] for m in metrics),
        TEMPLATES=sum(m['TEMPLATES'] for m in metrics),
        TEMPLATES_MARKED_DUPLICATE=sum(
            m['TEMPLATES_MARKED_DUPLICATE'] for m in metrics))
    agg['TEMPLATE_DUPLICATE_FRACTION'] = (
            agg['TEMPLATES_MARKED_DUPLICATE']/agg['TEMPLATES'])

    LOGGER.info(MSG_ALIGNMENTS, agg['ALIGNMENTS'])
    LOGGER.info(MSG_ALIGNMENTS_MARKED_DUPLICATE,
                agg['ALIGNMENTS_MARKED_DUPLICATE'])
    LOGGER.info(MSG_TEMPLATES, agg['TEMPLATES'])
    LOGGER.info(MSG_TEMPLATES_MARKED_DUPLICATE,
                agg['TEMPLATES_MARKED_DUPLICATE'])
    LOGGER.info(MSG_TEMPLATE_DUPLICATE_FRACTION,
                agg['TEMPLATE_DUPLICATE_FRACTION'])
    metfh.write(
        # kludge to output rounded floats
        # https://stackoverflow.com/a/29066406/6705037
        json.dumps(
            json.loads(
                json.dumps(agg),
                parse_float=lambda x: round(float(x), 4)),
            indent=2,
            sort_keys=True))


def parse_cmdargs(args: List[str]) -> argparse.Namespace:
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
    templatereads = parser.add_mutually_exclusive_group(required=True)
    templatereads.add_argument('--single',
                               dest='reads_per_template',
                               action='store_const',
                               const=1)
    templatereads.add_argument('--paired',
                               dest='reads_per_template',
                               action='store_const',
                               const=2)
    parser.add_argument('-n', '--n-items',
                        type=int,
                        default=DEFAULT_NITEMS,
                        help=('Expected maximum number of templates n '
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
    parser.add_argument('--strip-previous',
                        action='store_true',
                        help=('Unset duplicate flag bit for any reads that '
                              'have it set and are no longer considered '
                              'duplicate (default=False). Not required unless '
                              'records have previously been through a '
                              'duplicate marking step, in which case it is '
                              'strongly recommended for sensible results.'))
    parser.add_argument('--version',
                        action='version',
                        version=VERSION)
    return parser.parse_args(args)


def pgline(last: str) -> str:
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


def readends(qnamegrp: List[List[str]]) -> List[Tuple[str, int, str]]:
    """
    Calculate ends of the fragment, accounting for soft-clipped bases.

    Args:
        qnamegrp: QNAME group of SAM records, each record supplied as a list
            of str — i.e. the result of calling .split(TAB) on a SAM text line.

    Returns:
        For single-end reads, a one-element list:

            [(rname, ref_start, orientation)]

        For paired-end reads, a two-element coordinate-sorted list:

            [(left_rname, left_ref_start, left_orientation),
                (right_rname, right_ref_end, right_orientation)]

        Unmapped reads have the special end value `UNMAPPED` that sorts after
        any mapped end. Hence pairs with no aligned primary reads return:

            [UNMAPPED, UNMAPPED]

        and pairs with one aligned primary read return:

            [(rname, ref_start, orientation), UNMAPPED]
    """
    ends = []
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
            ends.append(UNMAPPED)
        # forward
        elif not flag & FLAG_REVERSE:
            # Leading soft clips.
            leading_s = int(match[1]) if (match := RE_LEADING_S(cigar)) else 0
            ends.append((rname, ref_start - leading_s, 'F'))
        # reverse
        else:
            # Trailing soft clips
            trailing_s = int(match[1]) if (match := RE_TRAILING_S(cigar)) else 0
            ref_end = ref_start
            for num, op in RE_CIGAR(cigar):
                ref_end += (int(num) if op in CIGAR_CONSUMES_REF else 0)
            ends.append((rname, ref_end + trailing_s, 'R'))

    # Canonical ordering: l < r and UNMAPPED is always last by construction.
    ends.sort()
    return ends


def unmarkdup(record: List[str]) -> None:
    """
    If the record is currently marked as a duplicate, update it in place to
    remove the duplicate flag and add or update the PG:Z: tag. This handles
    the case where a record was marked duplicate by an earlier PG but is not
    considered as duplicate by streammd.

    Args:
        record: list of SAM record str tokens.
    """
    flag = int(record[1])
    if flag | FLAG_DUPLICATE:
        record[1] = str(flag ^ FLAG_DUPLICATE)
        pg_old, pg_new = None, f'{PGTAG}:{PGID}'
        for i, opt in enumerate(record[SAM_OPTS_IDX:], SAM_OPTS_IDX):
            if opt.startswith(PGTAG):
                pg_old = opt
                record[i] = pg_new
        if not pg_old:
            record.append(pg_new)


def main() -> None:
    """
    Run as CLI script.
    """
    args = parse_cmdargs(sys.argv[1:])
    LOGGER.info(MSG_VERSION, VERSION)
    LOGGER.info(' '.join(sys.argv))
    nconsumers = args.consumer_processes
    reads = args.reads_per_template
    inq: 'Queue[str]' = Queue(20 * nconsumers)
    outq: 'Queue[Metrics]' = Queue(nconsumers)
    with (SharedMemoryManager() as smm,
          open(args.input) as infh,
          open(args.output, 'w') as outfh,
          open(args.metrics, 'w') as metfh):
        infd, outfd = infh.fileno(), outfh.fileno()
        reader_args = (infd, outfd, inq, nconsumers, args.input_batch_size)
        reader = Process(target=input_alnfile, args=reader_args)
        reader.start()
        bf = BloomFilter(smm, args.n_items, args.fp_rate)
        markdups_args = (bf.config, inq, outq, outfd, args.reads_per_template,
                         args.strip_previous)
        consumers = [
            Process(target=markdups, args=markdups_args)
            for _ in range(nconsumers)
        ]
        list(map(lambda x: x.start(), consumers))
        list(map(lambda x: x.join(), consumers))
        reader.join()
        LOGGER.info(MSG_UNIQUE_ITEMS_APPROXIMATE, bf.count())
        output_metrics([outq.get() for _ in range(nconsumers)], metfh)


if __name__ == '__main__':
    main()
