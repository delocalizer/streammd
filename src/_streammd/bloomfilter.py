"""
A Bloom filter implementation.

Or should I say, "yet another...". I justify this one with the fact that it's
pretty fast, and usable with shared memory in multiprocessing applications
thanks to bitarray's memoryview interface.
"""
import logging
from dataclasses import dataclass
from math import ceil, log
from multiprocessing.shared_memory import SharedMemory
from multiprocessing.managers import SharedMemoryManager
from sys import getsizeof
from typing import Callable, List, Optional, Tuple
from bitarray import bitarray
# tests as essentially the same speed as xxhash but better distribution
from farmhash import hash64withseed
from humanfriendly import parse_size

KMAX = 100
LOGGER = logging.getLogger(__name__)
MSG_INIT = 'BloomFilter initialized with n=%s, p=%s, m=%s, k=%s'
MSG_NOMEM = 'No solution for mem=%s with n=%s p=%s k<=%s'
# any 2 primes should do
SEED1 = 43
SEED2 = 9967


class NoMemorySolution(ValueError):
    """
    Raise when target n, p cannot bet met with the memory specification.
    """


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

    def __init__(self,
                 smm: SharedMemoryManager,
                 n: int,
                 p: float,
                 mem: Optional[str]=None) -> None:
        """
        Args:
            smm: SharedMemoryManager instance.
            n: Number of items to add.
            p: False positive rate when n items are added.
            mem: Optional human-friendly byte size for the bitarray e.g.
                "4GiB". A value of mem that is a power of 2 optimizes bitarray
                update performance. If mem is not supplied then the optimal
                (minimum) value for m that satisfies n and p will be used,
                usually at the cost of higher k (lower performance).
                NoMemorySolution is raised if the target n and p cannot be
                achieved with the requested mem.
        """
        self.n, self.p = int(n), p
        self.m, self.k = (self.m_k_mem(self.n, self.p, mem) if mem else
                          self.m_k_min(self.n, self.p))
        self.hash = self._hasher()

        # set up memory
        self.shm_bits = smm.SharedMemory(getsizeof(bitarray(self.m)))

        # initialize
        self.bits = bitarray(buffer=self.shm_bits.buf)
        self.bits[:] = 0
        LOGGER.info(MSG_INIT, self.n, self.p, self.m, self.k)

    @classmethod
    def copy(cls, config: BloomFilterConfig) -> 'BloomFilter':
        """
        Copy state from another BloomFilter instance, referencing the same
        shared memory.

        Args:
            config: Configuration of the template instance.

        Returns:
            BloomFilter
        """
        instance = object.__new__(cls)
        instance.shm_bits = SharedMemory(name=config.shm_bits)
        instance.n = config.n
        instance.p = config.p
        instance.m = config.m
        instance.k = config.k
        instance.hash = instance._hasher()
        instance.bits = bitarray(buffer=instance.shm_bits.buf)
        return instance

    @classmethod
    def m_k_mem(cls, n: int, p: float, mem: str) -> Tuple[int, int]:
        """
        Return the number of bits m and number of hash functions k for a Bloom
        filter containing n items with a false positive rate of p, where m
        will occupy (approximately) mem bytes. For fixed n and p, k is highly
        sensitive to m around the memory-optimal value, and even a slightly
        higher value for m can result in a significantly lower value for k and
        thus higher performance.

        Args:
            n: Number of items to add.
            p: Desired false positive rate after n items are added.
            mem: Human-friendly byte size for the bitarray e.g. "4GiB". A
                value that is a power of 2 optimizes bitarray update
                performance. NoMemorySolution is raised if the target n and p
                cannot be achieved with the requested mem.

        Returns:
            (m, k)
        """
        m = parse_size(mem) * 8
        # no closed form estimate for k in this case, solve by evaluation
        for k in range(1, KMAX+1):
            if pow(1 - pow(1-1/m, k*n), k) < p:
                return (m, k)
        raise NoMemorySolution(MSG_NOMEM % (mem, n, p, KMAX))

    @classmethod
    def m_k_min(cls, n: int, p: float) -> Tuple[int, int]:
        """
        Return the memory-optimal number of bits m and number of hash
        functions k for a Bloom filter containing n items with a false
        positive rate of p.

        Args:
            n: Number of items to add.
            p: Desired false positive rate after n items are added.

        Returns:
            (m, k)
        """
        m = ceil(-n * log(p) / (log(2) ** 2))
        k = ceil(log(2) * m/n)
        return (m, k)

    def __contains__(self, item: bytes|str) -> bool:
        """
        Test if an item is present.

        Args:
            item: Item to check for.

        Returns True if item is present, False otherwise.
        """
        return all(self.bits[pos] for pos in self.hash(item))

    def __del__(self) -> None:
        """
        Cleanup on aisle 3.
        """
        try:
            del self.bits
        except AttributeError:
            pass

    @property
    def config(self) -> BloomFilterConfig:
        """
        Returns configuration for this instance.
        """
        return BloomFilterConfig(
                self.shm_bits.name, self.n, self.p, self.m, self.k)

    @property
    def mpow2(self) -> bool:
        """
        Returns True if self.m is a power of 2, False otherwise.
        """
        return self.m & (self.m - 1) == 0

    def add(self, item: bytes|str) -> bool:
        """
        Add an item.

        Args:
            item: Item to add.

        Returns True if item was added, False if it was already present.

        Note:
        For the sake of speed, this implementation does not use locking and
        provides no guarantee of consistency across multiple threads or
        processes accessing the shared memory. This -potentially- matters in
        cases where the filter is being updated at the same time it is being
        used to check for membership, but even then the specific requirements
        of the application determine if it's a real problem. A couple of
        practical matters to consider:
        1. Fewer threads/processes = less of an issue.
        2. Randomly and infrequently distributed duplicates in the inputs =
           less of an issue.
        If in doubt tests should be run using real data for your application
        to check if the behaviour of this implementation is acceptable.
        """
        added = False
        for pos in self.hash(item):
            if not self.bits[pos]:
                self.bits[pos] = added = True
        return added

    def count(self) -> int:
        """
        Return the approximate number of items stored.

        Ref: Swamidass & Baldi (2007) https://doi.org/10.1021/ci600358f
        """
        return ceil((-self.m/self.k) * log(1 - (self.bits.count(1)/self.m)))

    def _hasher(self) -> Callable[[bytes|str], List[int]]:
        """
        Return a function that hashes to k independent integer outputs of
        size m.

        k linear combinations of just 2 independent hashes has the same
        asymptotic behaviour as k independent hashes.

        Ref: Kirsch & Mitzenmacher (2006) https://doi.org/10.1007/11841036_42
        """
        m, ks = self.m, tuple(range(self.k))
        mask = m - 1 if self.mpow2 else None

        def _hasher_mod(u):
            h1 = hash64withseed(u, SEED1)
            h2 = hash64withseed(u, SEED2)
            return [(h1 + i*h2) % m for i in ks]

        def _hasher_msk(u):
            h1 = hash64withseed(u, SEED1)
            h2 = hash64withseed(u, SEED2)
            return [(h1 + i*h2) & mask for i in ks]

        if self.mpow2:
            return _hasher_msk
        return _hasher_mod
