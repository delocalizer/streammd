"""""
Test markdups module.
"""
from multiprocessing.managers import SharedMemoryManager
from multiprocessing import Queue
from importlib.resources import files
from unittest import TestCase

from streammd.bloomfilter import BloomFilter
from streammd.markdups import *

RESOURCES = files('tests.streammd.resources')


class TestMarkDups(TestCase):
    """
    Test markdups module functions.
    """

    def test_samrecords_header(self):
        """
        Confirm that ValueError is raised if SAM file input lacks header.
        """
        samq = Queue(1000)
        headerq = Queue(1000)
        with RESOURCES.joinpath('no_header.sam') as infile:
            with self.assertRaises(ValueError, msg=MSG_NOHEADER):
                samrecords(headerq, samq, 1, 50, infile)

    def test_samrecords_qnamegrouped(self):
        """
        Confirm that ValueError is raised if SAM file records are not
        grouped by qname.
        """
        samq = Queue(1000)
        headerq = Queue(1000)
        with RESOURCES.joinpath('not_qnamegrouped.sam') as infile:
            with self.assertRaises(ValueError, msg=MSG_QNAMEGRP):
                samrecords(headerq, samq, 1, 50, infile)

    def test_samrecords_read_and_enqueue(self):
        """
        Confirm that SAM file header is written headerq and SAM file records
        are written to samq in batched groups as expected.
        """
        samq = Queue(1000)
        headerq = Queue(1000)
        nconsumers = 1
        with RESOURCES.joinpath('6_good_records.sam') as infile:
            samrecords(headerq, samq, nconsumers, 50, infile)
        # one header per consumer
        self.assertEqual(headerq.qsize(), nconsumers)
        # one batch (3 QNAME groups < batchsize) + one sentinel per consumer
        self.assertEqual(samq.qsize(), 1 + nconsumers)
        # 3 QNAME groups in the batch, each with one pair
        batch = samq.get()
        self.assertEqual(len(batch), 3)
        for group in batch:
            self.assertEqual(len(group), 2)

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
