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
        samq: multiprocessing.Queue to get batches of SAM records.
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
        for line in batch:
            alignment = AlignedSegment.fromstring(line, header)
            ends = readends(alignment)
            if ends and bf.add(ends):
                alignment.flag += 1024
            # in contrast to sys.stdout.write, os.write is atomic so we get
            # whole lines in the output and we don't have to care about
            # locking or using an output queue
            os.write(outfd, (alignment.to_string()+'\n').encode('utf-8'))


def readends(alignment):
    """
    Generate a signature for a read, symmetric under flip of orientation.
    """
    flag = alignment.flag
    read_unmapped = bool(flag & 4)
    mate_unmapped = bool(flag & 8)
    read_reverse = bool(flag & 16)
    mate_reverse = bool(flag & 32)
    first_in_pair = bool(flag & 64)
    second_in_pair = bool(flag & 128)
    # to start let's just consider mapped pairs
    # TODO: allow single ends to be mapped
    if read_unmapped or mate_unmapped:
        return None
    # TODO: handle ff and rr pairs
    if not (read_reverse ^ mate_reverse):
        return None
    ends = (
        (alignment.reference_name,
            alignment.reference_start - alignment.query_alignment_start + 1),
        (alignment.next_reference_name,
            alignment.next_reference_start))
    return f'{ends[0][0]}{ends[0][1]}{ends[1][0]}{ends[1][1]}'


def samrecords(headerq, samq, nconsumers, batchsize=50, infd=0, outfd=1):
    """
    Read records from a SAM file input stream and enqueue them in batches.

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
    samlines = []
    headlines = []
    header = None
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
            if len(samlines) == batchsize:
                samq.put(samlines)
                samlines = []
    samq.put(samlines)
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
