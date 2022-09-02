"""
A Bloom filter implementation.
"""
import logging
from dataclasses import dataclass
from math import ceil, log
from multiprocessing.shared_memory import SharedMemory
from multiprocessing.managers import SharedMemoryManager
from sys import getsizeof
from typing import Callable, List, Optional, Tuple, TypedDict
from bitarray import bitarray
# tests as essentially the same speed as xxhash but much better distribution
from farmhash import hash64withseed


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
MSG_INIT = 'BloomFilter initialized with n=%s, p=%s'
MSG_NADDED_GT_N = 'approx number of added items %s now exceeds target: %s'


@dataclass
class BloomFilterConfig:
    """
    Configuration for a BloomFilter.
    """
    shm_bits: str
    n: int
    p: float
    m: int
    k: int


class BloomFilter:
    """
    A Bloom filter implementation.
    """
    default_n: int = int(1e9)
    default_p: float = 1e-9

    def __init__(self,
                 smm: SharedMemoryManager,
                 n: int=default_n,
                 p: float=default_p) -> None:
        """
        Args:
            smm: SharedMemoryManager instance.
            n: number of items to add (default=1e9).
            p: false positive rate when n items are added (default=1e-9).
        """
        self.n, self.p = n, p
        self.m, self.k = self.optimal_m_k(n, p)
        self.hash = self.hasher(self.m, self.k)

        # set up memory
        self.shm_bits = smm.SharedMemory( getsizeof(bitarray(self.m)))

        # initialize
        self.bits = bitarray(buffer=self.shm_bits.buf)
        self.bits[:] = 0

    def __contains__(self, item: bytes|str) -> bool:
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

    def __del__(self) -> None:
        """
        Cleanup on aisle 3.
        """
        del self.bits

    def add(self, item: bytes|str) -> bool:
        """
        Add an item.

        Args:
            item: item to add.

        Returns True if item was added, False if it was already present.

        Note:
        For the sake of speed, this implementation does not use locking and
        provides no guarantee of consistency across multiple threads or
        processes accessing the shared memory. This potentially matters in
        cases where the filter is being updated at the same time it is being
        used to check for membership, but even then the specific requirements
        of the application determine if it's a real problem. A couple of
        practical matters to consider:
        1. Much more time is spent in the hash and modulus functions than in
           getting/setting bits.
        2. Fewer threads = less of an issue.
        3. Randomly and infrequently distributed duplicates in the inputs =
           less of an issue.
        If in doubt tests should be run with real data to check if the
        behaviour of this implementation is acceptable. My testing showed that
        using multiprocessing.Lock slowed performance by about 20x.
        """
        added = False
        for pos in self.hash(item):
            if self.bits[pos] == 0:
                added = True
                self.bits[pos] = 1
        return added

    @property
    def config(self) -> BloomFilterConfig:
        """
        Returns configuration for this instance.
        """
        return BloomFilterConfig(
                self.shm_bits.name, self.n, self.p, self.m, self.k)

    def count(self) -> int:
        """
        Return the approximate number of items stored.

        Ref: Swamidass & Baldi (2007)
        """
        # Implementation note:
        # Using xxhash.xxh3_64_intdigest instead of farmhash in hasher(...)
        # produces a poor result from this approximation: sometimes ~ 20% away
        # from the real value. This suggests it isn't well distributed.
        return ceil((-self.m/self.k) * log(1 - (self.bits.count(1)/self.m)))

    @classmethod
    def copy(cls, config: BloomFilterConfig) -> 'BloomFilter':
        """
        Copy state from another BloomFilter instance, referencing the same
        shared memory.

        Args:
            config: configuration dict of the template instance.

        Returns:
            BloomFilter
        """
        instance = object.__new__(cls)
        instance.shm_bits = SharedMemory(name=config.shm_bits)
        instance.n = config.n
        instance.p = config.p
        instance.m = config.m
        instance.k = config.k
        instance.hash = instance.hasher(instance.m, instance.k)
        instance.bits = bitarray(buffer=instance.shm_bits.buf)
        return instance

    @classmethod
    def hasher(cls,
               m: int,
               k: int,
               seeds: Optional[List[int]]=None
               ) -> Callable[[bytes|str], List[int]]:
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
            return [hash64withseed(item, seed) % m for seed in seeds]
        return _hasher

    @classmethod
    def optimal_m_k(cls, n: int, p: float) -> Tuple[int, int]:
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
