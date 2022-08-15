"""
Mark duplicates on SAM file input stream.
"""
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import Queue, Process
from multiprocessing.managers import SharedMemoryManager
import logging
import os
from bloomfilter import BloomFilter
from pysam import AlignmentHeader, AlignedSegment

DEFAULT_FPRATE = 1e-6
DEFAULT_NITEMS = 1e9
DEFAULT_NWORKERS = 10
LOGGER = logging.getLogger(__name__)
SENTINEL = 'STOP'


def markdups(bfconfig, headerq, samq, outfd=1):
    """
    Process SAM file records.

    Args:
        bfconfig: Bloom filter configuration dict.
        headerq: multiprocessing.Queue to get header.
        samq: multiprocessing.Queue to get batches of paired SAM records.
        outfd: output stream file descriptor (default=1).

    Returns:
        None
    """
    bf = BloomFilter.copy(bfconfig)
    header = AlignmentHeader.from_text(headerq.get())
    while True:
        batch = samq.get()
        if batch == SENTINEL:
            break
        for group in batch:
            alignments = [AlignedSegment.fromstring(r, header) for r in group]
            ends = readends(alignments)
            if ends and bf.add(ends):
                for a in alignments:
                    a.flag += 1024
            # Write the group as a group. In contrast to sys.stdout.write,
            # os.write is atomic so we don't have to care about locking or
            # using an output queue.
            out = '\n'.join(a.to_string() for a in alignments) + '\n'
            os.write(outfd, out.encode('utf-8'))


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
    # TODO: allow single ends to be mapped
    if r1.is_unmapped or r2.is_unmapped:
        return None
    ends = [None, None]
    for i, r in enumerate((r1, r2)):
        if r.is_forward:
            leading_S = r.cigar[0][1] if r.cigar[0][0] == 4 else 0
                                                          # leading soft clips
            ends[i] = r.reference_name, r.reference_start - leading_S
        else:
            trailing_S = r.cigar[-1][1] if r.cigar[-1][0] == 4 else 0
                                                        # trailing soft clips
            ends[i] = r.reference_name, r.reference_end + trailing_S
    orientation = '^' if r1.is_reverse ^ r2.is_reverse else '='
    # define canonical ordering: l < r
    ends.sort()
    (l_rname, l_pos), (r_rname, r_pos) = ends
    return f'{orientation}{l_rname}{l_pos}{r_rname}{r_pos}'


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


def main():

    nconsumers = DEFAULT_NWORKERS
    headerq = Queue(nconsumers)
    samq = Queue(1000)
    producer = Process(target=samrecords, args=(headerq, samq, nconsumers))
    producer.start()
    with SharedMemoryManager() as smm:
        nitems = DEFAULT_NITEMS
        fprate = DEFAULT_FPRATE
        bf = BloomFilter(smm, nitems, fprate)
        consumers = [
            Process(target=markdups, args=(bf.config, headerq, samq))
            for _ in range(nconsumers)
        ]
        list(map(lambda x: x.start(), consumers))
        list(map(lambda x: x.join(), consumers))
        producer.join()


if __name__ == '__main__':
    main()
