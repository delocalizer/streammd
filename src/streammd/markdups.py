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
DEFAULT_NITEMS = 1e6
DEFAULT_NWORKERS = 4
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
        for r1, r2 in batch:
            pair = [AlignedSegment.fromstring(r, header) for r in (r1, r2)]
            ends = readends(pair)
            if ends and bf.add(ends):
                for alignment in pair:
                    alignment.flag += 1024
            # Write the pair as a pair. In contrast to sys.stdout.write,
            # os.write is atomic so we don't have to care about locking or
            # using an output queue.
            out = pair[0].to_string() + '\n' + pair[1].to_string() + '\n'
            os.write(outfd, out.encode('utf-8'))


def readends(pair):
    """
    Generate a signature for a read, symmetric under flip of orientation.
    """
    r1, r2 = pair
    assert r1.qname == r2.qname, f'{r1.qname} != {r2.qname}'
    # to start let's just consider mapped pairs
    # TODO: allow single ends to be mapped
    if r1.is_unmapped or r2.is_unmapped:
        return None
    # TODO: handle ff and rr pairs
    if not (r1.is_reverse ^ r2.is_reverse):
        return None
    ends = list(sorted((
        (r1.reference_name, r1.pos - r1.qstart + 1),
        (r2.reference_name, r2.pos - r2.qstart + 1))))
    with open('/tmp/readends', 'a') as fh:
        fh.write(str(ends))
    return f'{ends[0][0]}{ends[0][1]}{ends[1][0]}{ends[1][1]}'


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
    # FIXME this doesn't handle reads with secondary and supplementary
    # alignments
    def batch(lines):
        """
        For a list of reads [r1.1, r1.2, r2.1, r2.2, ...] return the list of
        pairs [(r1.1, r1.2), (r2.1, r2.2), ...]
        """
        return [(r1, r2) for (r1, r2) in zip(lines[::2], lines[1::2])]

    header = None
    headlines = []
    samlines = []
    for line in os.fdopen(infd):
        if line.startswith('@'):
            headlines.append(line)
            os.write(outfd, line.encode('ascii'))
        else:
            if not header:
                header = ''.join(headlines)
                for _ in range(nconsumers):
                    headerq.put(header)
            samlines.append(line.strip())
            if len(samlines) == 2 * batchsize:
                samq.put(batch(samlines))
                samlines = []
    samq.put(batch(samlines))
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

#from random import shuffle
#
#def add_batch(bf_config, items):
#    bf = BloomFilter.copy(bf_config)
#    dupes = 0
#    for item in items:
#        present = bf.add(item)
#        if present:
#            dupes += 1
#    return dupes
#
#
#def main():
#    """
#    Test it
#    """
#    nconsumers = 1 
#    with SharedMemoryManager() as smm,\
#            ProcessPoolExecutor(max_workers=nconsumers) as ppe:
#
#        target_size = int(1e7)
#        bf = BloomFilter(smm, target_size, 1e-7)
#
#        # Here we assume we know the list of items in advance, so we can
#        # construct N chunks of items in advance. When we're reading through a
#        # bam in real time that won't be the case but we can do something like
#        # create 1 reader per contig, so 'items' is the list of contigs and
#        # the mapped func takes a contig name as an argument. TODO What to do
#        # for coordinate unsorted bams though?
#
#        # this creates 4 duplicates of every unique value
#        items = list(str(i) for i in range(int(target_size/5))) * 5
#        shuffle(items)
#
#
#        def chunker(l, n):
#            """
#            yield striped chunks
#            """
#            for i in range(0, n):
#                yield l[i::n]
#
#        chunks = chunker(items, nconsumers)
#
#        batch_args = ((bf.config, chunk) for chunk in chunks)
#        dupes = ppe.map(add_batch, *zip(*batch_args))
#        ppe.shutdown()
#
#        print(f'{len(items)} total items added (true)')
#        print(f'{bf.count()} unique items added (approx)')
#        print(f'{sum(dupes)} duplicates')
#        check = [0, 1, 10, 100, 1000, 10000, 100000, 1000000, 2000000, 10000000]
#        for j in check:
#            in_filt = str(j) in bf
#            print(f'{j} {"is" if in_filt else "is NOT"} present')
#    #for k in range(target_size, 2*target_size):
#    #    in_filt = str(k) in bf
#    #    if in_filt:
#    #        print(f'FP {k}')
#
#if __name__ == '__main__':
#    main()
