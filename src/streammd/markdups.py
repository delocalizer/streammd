from concurrent.futures import ProcessPoolExecutor
from multiprocessing.managers import SharedMemoryManager
from bloomfilter import BloomFilter


def add_batch(bf_vars, bf_bits, items):
    bf = BloomFilter.copy(bf_vars, bf_bits)
    for item in items:
        bf.add(item)


def main():
    """
    Test it
    """
    workers = 6 
    with SharedMemoryManager() as smm,\
            ProcessPoolExecutor(max_workers=workers) as ppe:

        bf = BloomFilter(smm, int(1e6), 1e-1, hard_limit=False)

        # Here we assume we know the list of items in advance, so we can
        # construct N chunks of items in advance. When we're reading through a
        # bam in real time that won't be the case but we can do something like
        # create 1 reader per contig, so 'items' is the list of contigs and
        # the mapped func takes a contig name as an argument.
        items = list(str(i) for i in range(1000001))

        def chunker(l, n):
            """
            yield striped chunks
            """
            for i in range(0, n):
                yield l[i::n]

        chunks = chunker(items, workers)

        batch_args = ((bf.shl_vars.shm.name, bf.shm_bits.name, chunk)
                      for chunk in chunks)
        ppe.map(add_batch, *zip(*batch_args)) 

    check = [0, 1, 10, 100, 1000, 10000, 100000, 1000000, 2000000, 10000000]
    for j in check:
        in_filt = str(j) in bf
        print(f'{j} {"is" if in_filt else "is NOT"} present')


if __name__ == '__main__':
    main()
