"""
A Bloom filter implementation.
"""
import logging
from math import ceil, log
from multiprocessing.shared_memory import ShareableList, SharedMemory
from sys import getsizeof
from bitarray import bitarray
from xxhash import xxh3_64_intdigest


LOGGER = logging.getLogger(__name__)
PRIMES = [
    2, 3, 5, 7, 521, 11, 523, 13, 17, 19, 23, 29, 541, 31, 37, 41, 43, 47, 53,
    59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137,
    139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223,
    227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307,
    311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397,
    401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487,
    491, 499, 503, 509
]
MSG_NADDED_GT_N = 'approx number of added items %s now exceeds target: %s'


class BloomFilter:
    """
    A Bloom filter implementation.
    """
    default_n = int(1e9)
    default_p = 1e-9

    def __init__(self, smm, n=default_n, p=default_p, hard_limit=True):
        """
        Args:
            smm: SharedMemoryManager instance.
            n: number of items to add (default=1e9).
            p: false positive rate when n items are added (default=1e-9).
            hard_limit: raise ValueError when nadded > n; if False just log
                        a warning.
        """
        m, k = self.optimal_m_k(n, p)
        # set up memory
        self.shl_vars = smm.ShareableList([None] * 6)
        self.shm_bits = smm.SharedMemory(getsizeof(bitarray(m)))

        # initialize
        self.n = self.shl_vars[0] = n
        self.p = self.shl_vars[1] = p
        self.m = self.shl_vars[2] = m
        self.k = self.shl_vars[3] = k
        self.nadded = self.shl_vars[4] = 0  # see add(self, item)
        self.hard_limit = self.shl_vars[5] = hard_limit
        self.hash = self.hasher(self.m, self.k)
        self.bits = bitarray(buffer=self.shm_bits.buf)
        self.bits[:] = 0

    def __contains__(self, item):
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

    def __del__(self):
        del self.bits

    def add(self, item):
        """
        Add an item.

        The value of self.nadded is incremented iff the filter misses so it is
        only accurate to within the false positive rate p and is always <=
        the true number of unique items added.

        Args:
            item: item to add.

        Returns True if item is already present, False otherwise.
        """
        present = True
        for pos in self.hash(item):
            if self.bits[pos] == 0:
                present = False
                self.bits[pos] = 1
        if not present:
            self.nadded += 1
        if self.nadded > self.n:
            if self.hard_limit:
                raise ValueError(MSG_NADDED_GT_N % (self.nadded, self.n))
            LOGGER.warning(MSG_NADDED_GT_N, self.nadded, self.n)
        return present

    @classmethod
    def copy(cls, shl_vars, shm_bits):
        """
        Copy state from another BloomFilter instance.

        Args:
            shl_vars: name of ShareableList instance.
            shm_bits: name of SharedMemory instance.
        """
        instance = object.__new__(cls)
        instance.shl_vars = ShareableList(name=shl_vars)
        instance.shm_bits = SharedMemory(name=shm_bits)
        instance.n = instance.shl_vars[0]
        instance.p = instance.shl_vars[1]
        instance.m = instance.shl_vars[2]
        instance.k = instance.shl_vars[3]
        instance.nadded = instance.shl_vars[4]
        instance.hard_limit = instance.shl_vars[5]
        instance.hash = instance.hasher(instance.m, instance.k)
        instance.bits = bitarray(buffer=instance.shm_bits.buf)
        return instance

    @classmethod
    def hasher(cls, m, k, seeds=None):
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
        def _hasher(item):
            return [xxh3_64_intdigest(item, seed) % m for seed in seeds]
        return _hasher

    @classmethod
    def optimal_m_k(cls, n, p):
        """
        Return the optimal number of bits m and optimal number of hash
        functions k for a Bloom filter containing n items with a false
        positive rate of p.

        Args:
            n: number of items to add.
            p: desired false positive rate after n items are added.

        Returns:
            (m, k)
        """
        m = ceil(-n * log(p) / (log(2) ** 2))
        k = ceil(log(2) * m/n)
        return (m, k)
