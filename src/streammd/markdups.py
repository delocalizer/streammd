#!/usr/bin/env
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import Manager, Pool, Process
from multiprocessing.managers import SharedMemoryManager
from random import shuffle
from sys import stdin, stdout
import os
from bloomfilter import BloomFilter
from pysam import AlignedSegment


def samrecords(infd, queue, outfd):
    """
    Read records from a SAM file input stream and add them to queue.
    Header lines are written directly to the output stream.
    """
    for line in os.fdopen(infd):
        if line.startswith('@'):
            os.write(outfd, line.encode('utf-8'))
        else:
            queue.put(line)


def markdups(queue, outfd):
    """
    Process SAM record.
    """
    while True:
        raw = queue.get()
        os.write(outfd, (raw[:-1] + f'\tprocessed by {os.getpid()}\n').encode('utf-8'))
    

def main():
    manager = Manager()
    samq = manager.Queue(1000)
    producer = Process(target=samrecords, args=(0, samq, 1))
    producer.start()
    nworkers = 5
    pool = Pool(nworkers)
    results = []
    for _ in range(nworkers):
        results.append(pool.apply_async(markdups, args=(samq, 1)))
    producer.join()
    producer.close()


if __name__ == '__main__':
    main()

# def read_bam(
#def add_batch(bf_vars, bf_bits, items):
#    bf = BloomFilter.copy(bf_vars, bf_bits)
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
#    workers = 10 
#    with SharedMemoryManager() as smm,\
#            ProcessPoolExecutor(max_workers=workers) as ppe:
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
#        chunks = chunker(items, workers)
#
#        batch_args = ((bf.shl_vars.shm.name, bf.shm_bits.name, chunk)
#                      for chunk in chunks)
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


if __name__ == '__main__':
    main()
