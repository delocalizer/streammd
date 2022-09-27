"""
Read a SAM file from STDIN, mark duplicates in a single pass and stream
processed records to STDOUT.

Input must begin with a valid SAM header, followed by qname-grouped records.

Default log level is 'INFO' — set to something else with the LOG_LEVEL
environment variable.
"""
from humanfriendly import format_size
from importlib import metadata
from multiprocessing import Process, SimpleQueue
from multiprocessing.managers import SharedMemoryManager
from typing import List, Literal, Tuple, TypedDict
import argparse
import json
import logging
import os
import re
import sys
from .bloomfilter import BloomFilter, BloomFilterConfig, NoMemorySolution

ROOTLOG = logging.getLogger()
ROOTLOG.setLevel(os.environ.get('LOG_LEVEL', 'INFO'))
ROOTLOG.addHandler(logging.StreamHandler())

LOGGER = logging.getLogger(__name__)

DEFAULT_FPRATE = 1e-6
DEFAULT_LOGINTERVAL = 1000000
DEFAULT_MEM = '4GiB'
DEFAULT_METRICS = 'streammd-metrics.json'
DEFAULT_NITEMS = 1e9
DEFAULT_NWORKERS = 2
DEFAULT_WORKQBATCHSIZE = 500

MSG_ALIGNMENTS = 'alignments seen: %s'
MSG_ALIGNMENTS_MARKED_DUPLICATE = 'alignments marked duplicate: %s'
MSG_MEM = 'Running with minimum required memory %s'
MSG_PARAMS = 'mem=%s; n=%.2E; p=%.2E; workers=%s; batchsize=%s;'
MSG_NOHEADER = 'no header lines detected'
MSG_NOTSINGLE = ('%s: expected 1 primary alignment, got %s. Input is not '
                 'single-end reads?')
MSG_NOTPAIRED = ('%s: expected 2 primary alignments, got %s. Input is not '
                 'paired-end reads or not qname grouped?')
MSG_NQNAMES = 'qnames read: %s'
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


def markdup(record: List[str]) -> None:
    """
    Mark the record as duplicate by updating in-place the FLAG and adding or
    updating the PG:Z: tag.

    Args:
        record: list of SAM record str tokens.

    Returns:
        None
    """
    flag = int(record[1])
    record[1] = str(flag | FLAG_DUPLICATE)
    pg_old, pg_new = None, f'{PGTAG}:{PGID}'
    for i, opt in enumerate(record[SAM_OPTS_IDX:], SAM_OPTS_IDX):
        if opt.startswith(PGTAG):
            pg_old = opt
            record[i] = pg_new
    if not pg_old:
        record.append(pg_new)


def markdups(bfconfig: BloomFilterConfig,
             workq: SimpleQueue,
             outq: SimpleQueue,
             resultq: SimpleQueue,
             reads_per_template: Literal[1, 2],
             strip_previous: bool=False) -> None:
    """
    Process SAM file records.

    Args:
        bfconfig: Bloom filter configuration dict.
        workq: SimpleQueue to get batches of qname grouped SAM records.
        outq: SimpleQueue to put processed SAM records.
        resultq: SimpleQueue to put metrics.
        reads_per_template: 1 or 2.
        strip_previous: unset duplicate flag bit for any reads that have it
            set and are no longer considered duplicate (default=False). Not
            necessary unless records have previously been through a duplicate
            marking step, in which case it is strongly recommended for
            sensible results.

    Results are added to the result queue as a Metrics instance:

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
        batch = workq.get()
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
                    markdup(read)
                    n_aln_dup += 1
            elif strip_previous:
                for read in qnamegrp:
                    unmarkdup(read)
            for read in qnamegrp:
                outlines.append('\t'.join(read))
        outq.put('\n'.join(outlines) + '\n')
    resultq.put(Metrics({
        'TEMPLATES':                   n_tpl,
        'TEMPLATES_MARKED_DUPLICATE':  n_tpl_dup,
        'ALIGNMENTS':                  n_aln,
        'ALIGNMENTS_MARKED_DUPLICATE': n_aln_dup,
    }))


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
    parser.add_argument('--single',
                        dest='reads_per_template',
                        action='store_const',
                        help='Accept single-ended reads as input (default is '
                        'paired-end).',
                        const=1,
                        default=2)
    parser.add_argument('-m', '--mem',
                        default=DEFAULT_MEM,
                        help='Human-friendly byte size for the Bloom filter '
                        f'(default="{DEFAULT_MEM}"). Higher mem => lower k => '
                        'faster hashing. A value of mem that is a power of 2 '
                        'also optimizes filter update performance.')
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
    parser.add_argument('-w', '--workers',
                        type=int,
                        default=DEFAULT_NWORKERS,
                        help=('Number of worker processes '
                              f'(default={DEFAULT_NWORKERS}).'))
    parser.add_argument('--strip-previous',
                        action='store_true',
                        help=('Unset duplicate flag bit for any reads that '
                              'have it set and are no longer considered '
                              'duplicate (default=False). Not required unless '
                              'records have previously been through a '
                              'duplicate marking step, in which case it is '
                              'strongly recommended for sensible results.'))
    parser.add_argument('--workq-batch-size',
                        type=int,
                        default=DEFAULT_WORKQBATCHSIZE,
                        help=('Specify the number of SAM records in each '
                              'batch for the work queue '
                              f'(default={DEFAULT_WORKQBATCHSIZE}). '
                              'The default value performs well in most cases; '
                              'adjustments for SAM record size and number of '
                              'worker processes may yield small performance '
                              'improvements.'))
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
    prev = None
    tkns = last.strip().split('\t')
    if tkns[0] == '@PG':
        prev = dict(tkn.split(':', 1) for tkn in tkns[1:])['ID']
    if prev:
        tags.insert(2, f'PP:{prev}')
    return '\t'.join(['@PG'] + tags) + '\n'


def read_input(inf: int|str,
               outq: SimpleQueue,
               workq: SimpleQueue,
               nworkers: int,
               batchsize: int=500) -> None:
    """
    Read records from a qname-grouped SAM file input stream and enqueue them
    for processing in batches.

    Header lines are written directly to the output queue.

    Args:
        inf: Input stream file name or descriptor.
        outq: multiprocesing.Queue to put header.
        workq: SimpleQueue to put SAM records.
        nworkers: Number of workers.
        batchsize: Number of qnames per batch put into the workq (default=500).
            The single producer runs at fixed speed so bigger batch size =>
            lower rate at which batches are added to the queue. Up to a point
            this reduces queue overhead (fewer put/get ops) but beyond that
            the rate of supply falls below the value required to keep all
            consumers continually fed and they waste time blocking on get
            calls.

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
    header = None
    headlines = []
    with open(inf, 'rt', encoding='utf-8') as infh:
        for line in infh:
            if header:
                qname = line.partition('\t')[0]
                if qname == qname_last:
                    qname_group.append(line)
                else:
                    n_qname += 1
                    if n_qname % DEFAULT_LOGINTERVAL == 0:
                        LOGGER.debug(MSG_NQNAMES, n_qname)
                    batch.append(qname_group)
                    batch_sz += 1
                    if batch_sz == batchsize:
                        workq.put(batch)
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
                header = ''.join(headlines)
                outq.put(header)
                qname_last = line.partition('\t')[0]
                qname_group = [line]
        batch.append(qname_group)
        workq.put(batch)
    for _ in range(nworkers):
        workq.put(SENTINEL)


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
    if flag & FLAG_DUPLICATE:
        record[1] = str(flag ^ FLAG_DUPLICATE)
        pg_old, pg_new = None, f'{PGTAG}:{PGID}'
        for i, opt in enumerate(record[SAM_OPTS_IDX:], SAM_OPTS_IDX):
            if opt.startswith(PGTAG):
                pg_old = opt
                record[i] = pg_new
        if not pg_old:
            record.append(pg_new)


def write_metrics(metrics: List[Metrics], metf: int|str) -> None:
    """
    Output aggregate metrics.

    Args:
        metrics: list of duplicate marking metrics from output queue.
        metf: Metrics file name or descriptor.
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
    with open(metf, 'wt', encoding='utf-8') as metfh:
        metfh.write(
            # kludge to output rounded floats
            # https://stackoverflow.com/a/29066406/6705037
            json.dumps(
                json.loads(
                    json.dumps(agg),
                    parse_float=lambda x: round(float(x), 4)),
                indent=2,
                sort_keys=True))


def write_output(outq: SimpleQueue, out: int|str):
    """
    Consume text items from outq and write them to out.

    Args:
        outq: SimpleQueue to get header and processed qname groups of SAM
            records.
        out: Output stream file name or descriptor.

    """
    with open(out, 'wt', encoding='utf-8') as outfh:
        while True:
            item = outq.get()
            if item == SENTINEL:
                outfh.flush()
                break
            outfh.write(item)


def main() -> None:
    """
    Run as CLI script.
    """
    args = parse_cmdargs(sys.argv[1:])
    inf = args.input
    out = args.output
    metf = args.metrics
    reads = args.reads_per_template
    mem = args.mem
    nitems = args.n_items
    fprate = args.fp_rate
    nworkers = args.workers
    strip_previous = args.strip_previous
    batchsize = args.workq_batch_size
    workq: 'SimpleQueue[str]' = SimpleQueue()
    outq: 'SimpleQueue[str]' = SimpleQueue()
    resultq: 'SimpleQueue[Metrics]' = SimpleQueue()

    LOGGER.info(MSG_VERSION, VERSION)
    LOGGER.info(' '.join(sys.argv))
    LOGGER.info(MSG_PARAMS, mem, nitems, fprate, nworkers, batchsize)

    with SharedMemoryManager() as smm:
        reader = Process(
            target=read_input, args=(inf, outq, workq, nworkers, batchsize))
        writer = Process(target=write_output, args=(outq, out))
        reader.start()
        writer.start()
        try:
            bf = BloomFilter(smm, nitems, fprate, mem)
        except NoMemorySolution as nms:
            LOGGER.warning(nms)
            bf = BloomFilter(smm, nitems, fprate)
            LOGGER.warning(MSG_MEM, format_size(bf.m // 8, binary=bf.mpow2))
        workers = [
            Process(
                target=markdups,
                args=(bf.config, workq, outq, resultq, reads, strip_previous))
            for _ in range(nworkers)
        ]
        list(map(lambda x: x.start(), workers))
        list(map(lambda x: x.join(), workers))
        outq.put(SENTINEL)
        reader.join()
        writer.join()
        LOGGER.info(MSG_UNIQUE_ITEMS_APPROXIMATE, bf.count())
        write_metrics([resultq.get() for _ in range(nworkers)], metf)


if __name__ == '__main__':
    main()
