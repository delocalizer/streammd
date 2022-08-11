"""
Mark duplicates on SAM file input stream.
"""
from multiprocessing import Manager, Pool, Process
from multiprocessing.managers import SharedMemoryManager
import os
from bloomfilter import BloomFilter
from pysam import AlignmentHeader, AlignedSegment


MAX_HEADER_SIZE = 1024 ** 2  # 1M
DEFAULT_WORKERS = 8

def samrecords(queue, sharedv, hlock, infd=0, outfd=1):
    """
    Read records from a SAM file input stream and add them to queue.

    Header lines are written directly to the output stream and also to the
    dict key sharedv['header']. When all header lines are written, hlock is
    released.
    
    Args:
        queue: multiprocessing.Queue to put SAM records
        sharedv: multiprocessing.Manager.dict; header text is written to the
            value of the 'header' key.
        hlock: released once sharedv['header'] contains the header text.
        infd: input stream file descriptor (default=0)
        outfd: output stream file descriptor (default=1)

    Returns:
        None
    """
    headlines = []
    header_complete = False
    for line in os.fdopen(infd):
        if line.startswith('@'):
            headlines.append(line)
            # os.write is atomic (unlike sys.stdout.write)
            os.write(outfd, line.encode('ascii'))
        else:
            if not header_complete:
                sharedv['header'] = ''.join(headlines)
                header_complete = True
                hlock.release()
            queue.put(line)


def markdups(queue, sharedv, hlock, outfd=1):
    """
    Process SAM record.

    Args:
        queue: multiprocessing.Queue to put SAM records
        sharedv: multiprocessing.Manager.dict; 'header' key is header text. 
        hlist: ShareableList; the 0th element stores the header text.
        outfd: output stream file descriptor (default=1)
    
    Returns:
        None
    """
    hlock.acquire()
    header = AlignmentHeader.from_text(sharedv['header'])
    hlock.release()
    while True:
        raw = queue.get()
        alignment = AlignedSegment.fromstring(raw, header)
        # os.write is atomic (unlike sys.stdout.write)
        os.write(outfd, (alignment.to_string() + '\n').encode('ascii'))
    

def main():
    manager = Manager()
    samq = manager.Queue(10000)
    sharedv = manager.dict()
    hlock = manager.Lock()
    hlock.acquire()
    producer = Process(target=samrecords, args=(samq, sharedv, hlock))
    producer.start()
    nworkers = DEFAULT_WORKERS
    pool = Pool(nworkers)
    for _ in range(nworkers):
        pool.apply_async(markdups, args=(samq, sharedv, hlock))
    # producer is done once stdin is exhausted
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
