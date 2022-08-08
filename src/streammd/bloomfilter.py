"""
A Bloom filter implementation.
"""

from math import ceil, log
from multiprocessing.shared_memory import SharedMemory
from sys import getsizeof, stderr
from bitarray import bitarray
from xxhash import xxh3_64_intdigest


PRIMES = [
    2, 3, 5, 7, 521, 11, 523, 13, 17, 19, 23, 29, 541, 31, 37, 41, 43, 47, 53,
    59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137,
    139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223,
    227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307,
    311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397,
    401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487,
    491, 499, 503, 509
]



class BloomFilter:
    """
    A Bloom filter implementation.
    """
    default_n = int(1e9)
    default_p = 1e-9

    def __init__(self, n=default_n, p=default_p):
        """
        Args:
            n: number of items to store (default=1e9).
            p: false positive rate when n items are stored (default=1e-9).
            TODO: hard_limit to toggle behaviour when stored > n
            TODO: make self.nstored a SharedMemory instance
            TODO: accept (bits, nstored) to allow initialization from
                  existing
        """
        self.n = n
        self.p = p
        self.m, self.k = self.bf_m_k(self.n, self.p)
        self.hash = self.get_hasher(self.m, self.k)

        # tmp create to find correct size
        _ba = bitarray(self.m)
        ba_size = getsizeof(_ba)
        del _ba
        # initialize
        self.shm = SharedMemory(size=ba_size, create=True)
        self.bits = bitarray(buffer=self.shm.buf)
        self.bits[:] = 0

    def __del__(self):
        del self.bits
        self.shm.close()
        self.shm.unlink()

    @classmethod
    def bf_m_k(cls, n, p):
        """
        Return the optimal number of bits m and optimal number of hash
        functions k for a Bloom filter to store n items with with max false
        positive rate of p.

        Args:
            n: number of items to store.
            p: desired false positive rate when n items are stored.

        Returns:
            (m, k)
        """
        m = ceil(-n * log(p) / (log(2) ** 2))
        k = ceil(log(2) * m/n)
        return (m, k)

    @classmethod
    def get_hasher(cls, m, k, seeds=None):
        """
        Return a function that hashes to k independent integer outputs of
        size m.

        Args:
            m: size of the output hash space.
            k: number of hash functions.
            seeds: (optional) an iterable of k integer seeds. If not supplied
            the first k primes will be used.
        """
        seeds = list(seeds) if seeds else sorted(list(PRIMES[:k]))
        assert len(seeds) == k
        def hasher(item):
            return [xxh3_64_intdigest(item, seed) % m for seed in seeds]
        return hasher

    def get(self, item):
        """
        Test if an item is present.

        Args:
            item: item to check for.

        Returns True if item is present, False otherwise.
        """
        for pos in self.hash(item):
            if self.bits[pos] == 0:
                return False
        return True

    def set(self, item):
        """
        Add an item.

        Args:
            item: item to add.

        Returns True if item is already present, False otherwise.
        """
        present = True
        for pos in self.hash(item):
            if self.bits[pos] == 0:
                present = False
                self.bits[pos] = 1
#        if not present:
#            self.stored += 1
#        if self.stored > self.n:
#            msg = f'{self.stored} stored items exceeds target: {self.n}'
#            if self.hard_limit:
#                raise ValueError(msg)
#            else:
#                print(msg, file=stderr)
        return present


def main():
    bf = BloomFilter(int(1e6), 1e-6)

    for i in range(1000001):
        bf.set(str(i))

    check = [0, 1, 10, 100, 1000, 10000000, 2000000]
    for j in check:
        in_filt = bf.get(str(j))
        print(f'{j} {"is" if in_filt else "is not"} present')


if __name__ == '__main__':
    main()
