"""
Mark duplicates on SAM file input stream.
"""
from multiprocessing import JoinableQueue, Process
from multiprocessing.managers import SharedMemoryManager
import os
from bloomfilter import BloomFilter
from pysam import AlignmentHeader, AlignedSegment

DEFAULT_WORKERS = 8
SENTINEL = 'STOP'


def samrecords(samq, headerq, workers, batchsize=50, infd=0, outfd=1):
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
    batch = []
    headlines = []
    header_done = False
    for line in os.fdopen(infd):
        if line.startswith('@'):
            headlines.append(line)
            # os.write is atomic (unlike sys.stdout.write)
            os.write(outfd, line.encode('ascii'))
        else:
            if not header_done:
                header_done = True
                header = ''.join(headlines)
                for _ in range(workers):
                    headerq.put(header)
            batch.append(line)
            if len(batch) == batchsize:
                samq.put(batch)
                batch = []
    samq.put(batch)
    for _ in range(workers):
        samq.put(SENTINEL)


def markdups(samq, headerq, outfd=1):
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
    header = AlignmentHeader.from_text(headerq.get(timeout=1))
    headerq.task_done()
    while True:
        lines = samq.get()
        if lines == SENTINEL:
            samq.task_done()
            break
        for line in lines:
            alignment = AlignedSegment.fromstring(line, header)
            # os.write is atomic (unlike sys.stdout.write)
            os.write(outfd, (alignment.to_string() + '\tdone\n').encode('ascii'))
        samq.task_done()


def main():
    workers = DEFAULT_WORKERS
    samq = JoinableQueue(10000)
    headerq = JoinableQueue(10)
    producer = Process(target=samrecords, args=(samq, headerq, workers))
    producer.start()
    consumers = [
        Process(target=markdups, args=(samq, headerq))
        for _ in range(workers)
    ]
    for c in consumers:
        c.start()
    for c in consumers:
        c.join()
        c.close()
    producer.join()
    producer.close()
    headerq.join()
    headerq.close()
    samq.join()
    samq.close()


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

