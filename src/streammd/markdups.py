"""
Read a SAM file from STDIN, mark duplicates in a single pass and stream
processed records to STDOUT.

Input must begin with a valid SAM header, followed by qname-grouped records.
Currently only paired reads are handled.

Default log level is 'INFO' — set to something else with the LOG_LEVEL
environment variable.
"""
from importlib import metadata
from itertools import repeat
from multiprocessing import Manager, Pool, Process
from multiprocessing.managers import SharedMemoryManager
import argparse
import logging
import os
import sys
from pysam import AlignmentHeader, AlignedSegment
from .bloomfilter import BloomFilter

DEFAULT_FPRATE = 1e-6
DEFAULT_NITEMS = int(1e9)
DEFAULT_NWORKERS = 8
DEFAULT_SAMQSIZE = 1000

LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(os.environ.get('LOG_LEVEL', 'INFO'))
LOGGER.addHandler(logging.StreamHandler())

MSG_NQNAME = 'qnames seen: %s'
MSG_NUNIQUE = 'approximate count of unique pairs: %s'
MSG_NDUP = 'approximate count of duplicates: %s'
MSG_DUPFRAC = 'approximate duplicate fraction: %0.4f'
MSG_VERSION = 'streammd version %s'
SENTINEL = 'STOP'


def samrecords(headerq, samq, nconsumers, batchsize=50, infd=0, outfd=1):
    """
    Read records from a qname-grouped SAM file input stream and enqueue them
    in batches.

    Header lines are written directly to the output stream and also to the
    header queue.

    Args:
        headerq: multiprocessing.Queue to put header.
        samq: multiprocessing.Queue to put SAM records.
        nconsumers: number of consumer processes.
        batchsize: number of lines per batch in samq (default=50).
        infd: input stream file descriptor (default=0).
        outfd: output stream file descriptor (default=1).

    Returns:
        None
    """
    group = []
    groupid = None
    header = None
    headlines = []
    batch = []
    for line in os.fdopen(infd):
        if line.startswith('@'):
            headlines.append(line)
            os.write(outfd, line.encode('ascii'))
        else:
            if not header:
                header = ''.join(headlines)
                for _ in range(nconsumers):
                    headerq.put(header)
            record = line.strip()
            qname = record.partition('\t')[0]
            if qname == groupid:
                group.append(record)
            else:
                if group:
                    batch.append(group)
                    if len(batch) == batchsize:
                        samq.put(batch)
                        batch = []
                groupid = qname
                group = [record]
    batch.append(group)
    samq.put(batch)
    for _ in range(nconsumers):
        samq.put(SENTINEL)


def markdups(bfconfig, headerq, samq, outfd=1):
    """
    Process SAM file records.

    Args:
        bfconfig: Bloom filter configuration dict.
        headerq: multiprocessing.Queue to get header.
        samq: multiprocessing.Queue to get batches of paired SAM records.
        outfd: output stream file descriptor (default=1).

    Returns:
        (n_qname, n_dupe): number of qnames and number of duplicates seen.
    """
    bf = BloomFilter.copy(bfconfig)
    header = AlignmentHeader.from_text(headerq.get())
    n_qname, n_dup = 0, 0
    while True:
        batch = samq.get()
        if batch == SENTINEL:
            break
        for group in batch:
            n_qname += 1
            alignments = [AlignedSegment.fromstring(r, header) for r in group]
            ends = readends(alignments)
            if ends and bf.add(ends):
                n_dup += 1
                for a in alignments:
                    a.flag += 1024
            # Write the group as a group. In contrast to sys.stdout.write,
            # os.write is atomic so we don't have to care about locking or
            # using an output queue.
            out = '\n'.join(a.to_string() for a in alignments) + '\n'
            os.write(outfd, out.encode('ascii'))
    return (n_qname, n_dup)


def readends(alignments):
    """
    Calculate ends of the fragment, accounting for soft-clipped bases.

    Args:
        alignments: qname group tuple of AlignedSegment instances.
    """
    r1, r2 = None, None
    # pick the primary alignments
    for alignment in alignments:
        if not (alignment.is_secondary or alignment.is_supplementary):
            if alignment.is_read1:
                r1 = alignment
            elif alignment.is_read2:
                r2 = alignment
    ends = [None, None]
    # TODO: allow single ends to be mapped
    if r1.is_unmapped or r2.is_unmapped:
        return None
    for i, r in enumerate((r1, r2)):
        if r.is_forward:
            front_s = r.cigar[0][1] if r.cigar[0][0] == 4 else 0
                                                          # leading soft clips
            ends[i] = r.reference_name, r.reference_start - front_s
        else:
            back_s = r.cigar[-1][1] if r.cigar[-1][0] == 4 else 0
                                                        # trailing soft clips
            ends[i] = r.reference_name, r.reference_end + back_s
    orientation = '^' if r1.is_reverse ^ r2.is_reverse else '='
    # canonical ordering: l < r
    ends.sort()
    (l_rname, l_pos), (r_rname, r_pos) = ends
    return f'{orientation}{l_rname}{l_pos}{r_rname}{r_pos}'


def parse_cmdargs(args):
    """
    Returns: Parsed arguments
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-n', '--n-items',
                        type=int,
                        default=DEFAULT_NITEMS,
                        help=('expected maximum number of read pairs n '
                              f'(default={DEFAULT_NITEMS}).'))
    parser.add_argument('-p', '--fp-rate',
                        type=float,
                        default=DEFAULT_FPRATE,
                        help=('target maximum false positive rate when n '
                              f'items are stored (default={DEFAULT_FPRATE}).'))
    parser.add_argument('--consumer-processes',
                        type=int,
                        default=DEFAULT_NWORKERS,
                        help=('Number of hashing processes '
                              f'(default={DEFAULT_NWORKERS}).'))
    parser.add_argument('--queue-size',
                        type=int,
                        default=DEFAULT_SAMQSIZE,
                        help=('size of the SAM record queue '
                              f'(default={DEFAULT_SAMQSIZE})'))
    parser.add_argument('--version',
                        action='version',
                        version=metadata.version('streammd'))
    return parser.parse_args(args)


def main():
    """
    Run as CLI script
    """
    args = parse_cmdargs(sys.argv[1:])
    LOGGER.info(MSG_VERSION, metadata.version('streammd'))
    LOGGER.info(' '.join(sys.argv))
    manager = Manager()
    headerq = manager.Queue(args.consumer_processes)
    samq = manager.Queue(args.queue_size)
    nconsumers = args.consumer_processes
    producer = Process(target=samrecords, args=(headerq, samq, nconsumers))
    producer.start()
    with SharedMemoryManager() as smm, Pool(nconsumers) as pool:
        bf = BloomFilter(smm, args.n_items, args.fp_rate)
        mdargs = repeat((bf.config, headerq, samq), nconsumers)
        counts = pool.starmap(markdups, mdargs)
        n_qname, n_dup = [sum(col) for col in zip(*counts)]
        n_unique = bf.count()
        producer.join()
    LOGGER.info(MSG_NQNAME, n_qname)
    LOGGER.info(MSG_NUNIQUE, n_unique)
    LOGGER.info(MSG_NDUP, n_dup)
    LOGGER.info(MSG_DUPFRAC, (n_dup/n_qname))


if __name__ == '__main__':
    main()
