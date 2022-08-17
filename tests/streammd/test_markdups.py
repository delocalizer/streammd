"""""
Test markdups module.
"""
from multiprocessing.managers import SharedMemoryManager
from unittest import TestCase

from streammd.bloomfilter import BloomFilter
from streammd.markdups import *


class TestMarkDups(TestCase):
    """
    Test markdups module functions.
    """

    def test_samrecords_header(self):
        """
        Confirm that ValueError is raised if SAM file input lacks header.
        """
        pass

    def test_samrecords_qnamegrouped(self):
        """
        Confirm that ValueError is raised if SAM file records are not
        grouped by qname.
        """
        pass

    def test_samrecords_operation(self):
        """
        Confirm that SAM file header is written to outfd and headerq, and
        SAM file records are written to samq.
        """
        pass

    def test_readends_1(self):
        """
        Confirm that ends of a duplicate pair in same orientation are the
        same.
        """
        pass

    def test_readends_2(self):
        """
        Confirm that ends of a duplicate pair in opposite orientation are the
        same.
        """
        pass

    def test_readends_3(self):
        """
        Confirm that soft clipping at fwd and reverse ends is handled
        as expected.
        """
        pass

    def test_readends_4(self):
        """
        Confirm that one end unmapped is handled as expected.
        """
        pass

    def test_readends_5(self):
        """
        Confirm that both ends unmapped is handled as expected.
        """
        pass

    def test_markdups_1(self):
        """
        Confirm that duplicates are marked as expected.
        """
        pass

    def test_main(self):
        """
        Confirm that main() operates as expected.
        """
        pass
