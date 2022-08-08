from bitarray import bitarray
from functools import partial
from math import ceil, log
from multiprocessing.shared_memory import SharedMemory
from sys import getsizeof
from xxhash import xxh3_64_intdigest


def bf_m_k(n, p):
    """
    Return the optimal number of bits m and optimal number of hash functions k
    for a Bloom filter to store n items with with max false positive rate of p.

    Args:
        n: number of elements to store.
        p: desired false positive rate when n elements are stored.

    Returns:
        (m, k)
    """
    m = ceil(-n * log(p) / (log(2) ** 2))
    k = ceil(log(2) * m/n)
    return (m, k)


def get_hasher(m, k, seeds=None):
    """
    Return a function that hashes to k independent integer outputs of size m.

    Args:
        m: size of the output hash space.
        k: number of hash functions.
        seeds: (optional) an iterable of k integer seeds. If not supplied, the
               first m primes will be used.
    """
    seeds = list(seeds) if seeds else sorted(list(primes(k)))
    assert len(seeds) == k
    def hasher(element):
        return [xxh3_64_intdigest(element, seed) % m for seed in seeds]
    return hasher


def primes(n):
    """
    Return first n primes.
    """
    # naïve but doesn't need to be fast
    def prime(i, primes):
        for prime in primes:
            if not (i == prime or i % prime):
                return False
        primes.add(i)
        return i
    primes = {2}
    i, p = 2, 0
    while True:
        if prime(i, primes):
            p += 1
            if p == n:
                return primes
        i += 1


DEFAULT_N = 1e9
DEFAULT_P = 1e-9

(m, k) = bf_m_k(DEFAULT_N, DEFAULT_P)

# create just to find size
_ba = bitarray(m)
ba_size = getsizeof(_ba)
del _ba

# initialize Bloom filter
sm = SharedMemory(size=ba_size, create=True)
ba = bitarray(buffer=sm.buf)
ba[:] = 0
hasher = get_hasher(m, k)


def get(ba, hasher, element):
    """
    Test if an element is present in a Bloom filter.

    Args:
        ba: Bloom filter bit array.
        hasher: function to hash into the array.
        element: element to test for.

    Returns True if element is present, False otherwise.
    """
    indexes = hasher(element)
    present = True
    for idx in indexes:
        if ba[idx] == 0:
            return False
    return present


def set(ba, hasher, element):
    """
    Add an element to Bloom filter.

    Args:
        ba: Bloom filter bit array.
        hasher: function to hash into the array.
        element: element to add.

    Returns True if element is already present, False otherwise.
    """
    indexes = hasher(element)
    present = True
    for idx in indexes:
        if ba[idx] == 0:
            present = False
            ba[idx] = 1
    return present


for i in range(1000000):
    set(ba, hasher, str(i))

check = [0, 1, 10, 100, 1000, 10000000, 2000000]
for j in check:
    present = get(ba, hasher, str(j))
    print(f'{check} {"is" if present else "is not"} present')
   

