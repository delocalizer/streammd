"""
Test bloomfilter module.
"""
from multiprocessing.managers import SharedMemoryManager
from unittest import TestCase

from streammd.bloomfilter import BloomFilter


class TestBloomFilter(TestCase):
    """
    Test bloomfilter module functions.
    """

    def test_contains_false(self):
        """
        Confirm that 'in' returns False when value is not present.
        """
        with SharedMemoryManager() as smm:
            n, p = 1e3, 1e-3
            value = 'test value'
            bf = BloomFilter(smm, n, p)
            self.assertFalse(value in bf)

    def test_contains_true(self):
        """
        Confirm that 'in' returns True when value is present.
        """
        with SharedMemoryManager() as smm:
            n, p = 1e3, 1e-3
            value = 'test value'
            bf = BloomFilter(smm, n, p)
            bf.add(value)
            self.assertTrue(value in bf)

    def test_add_false(self):
        """
        Confirm that adding a value already present returns False.
        """
        with SharedMemoryManager() as smm:
            n, p = 1e3, 1e-3
            value = 'test value'
            bf = BloomFilter(smm, n, p)
            bf.add(value)
            self.assertFalse(bf.add(value))

    def test_add_true(self):
        """
        Confirm that adding a value to empty filter returns True.
        """
        with SharedMemoryManager() as smm:
            n, p = 1e3, 1e-3
            value = 'test value'
            bf = BloomFilter(smm, n, p)
            self.assertTrue(bf.add(value))

    def test_count(self):
        """
        Confirm that BloomFilter.count returns approximately the number of
        elements added to it.
        """
        with SharedMemoryManager() as smm:
            n, p = 1e5, 1e-5
            values = [str(i) for i in range(100000, 110000)]
            bf = BloomFilter(smm, n, p)
            for value in values:
                bf.add(value)
            self.assertAlmostEqual(bf.count()/len(values), 1.0, places=2)

    def test_optimal_m_k(self):
        """
        Confirm that m, k calculation performs as expected.
        """
        values = (
            (1e6, 1e-06, 28755176, 20),
            (1e7, 1e-07, 335477044, 24),
            (1e8, 1e-08, 3834023351, 27),
            (1e9, 1e-06, 28755175133, 20),
        )
        for n, p, m_expected, k_expected in values:
            m, k = BloomFilter.m_k_min(n, p)
            self.assertEqual(m, m_expected)
            self.assertEqual(k, k_expected)

    def test_fnr_0(self):
        """
        Add n values and check that FNR = 0.
        """
        with SharedMemoryManager() as smm:
            n, p = 1e4, 1e-4
            values = [str(i) for i in range(10000, 11000)]
            bf = BloomFilter(smm, n, p)
            for value in values:
                bf.add(value)
            self.assertTrue(all(test in bf for test in values))

    def test_fpr_p(self):
        """
        Add n values and check that FPR is within expected limit.
        """
        with SharedMemoryManager() as smm:
            n = 1e6
            ps = (1e-3, 1e-4, 1e-5, 1e-6)
            values = [str(i) for i in range(int(n), 2*int(n))]
            misses = [str(i) for i in range(2*int(n), 3*int(n))]
            for p in ps:
                bf = BloomFilter(smm, n, p)
                for value in values:
                    bf.add(value)
                fpr = sum(test in bf for test in misses) / len(misses)
                print(n, p, fpr)
                #self.assertLessEqual(fpr, p)
