"""""
Test markdups module.
"""
from multiprocessing.managers import SharedMemoryManager
from multiprocessing import Queue
from importlib.resources import files
from tempfile import NamedTemporaryFile
from unittest import TestCase

from pysam import AlignmentFile

from streammd.bloomfilter import BloomFilter
from streammd.markdups import (DEFAULT_FPRATE,
                               DEFAULT_NITEMS,
                               MSG_NOHEADER,
                               MSG_QNAMEGRP,
                               UNMAPPED,
                               markdups,
                               readends,
                               samrecords)

RESOURCES = files('tests.streammd.resources')


class TestMarkDups(TestCase):
    """
    Test markdups module functions.
    """

    def test_samrecords_header(self):
        """
        Confirm that ValueError is raised if SAM file input lacks header.
        """
        samq = Queue(100)
        headerq = Queue(10)
        nconsumers = 1
        with (RESOURCES.joinpath('no_header.sam') as inf,
                NamedTemporaryFile() as out):
            with self.assertRaises(ValueError, msg=MSG_NOHEADER):
                samrecords(headerq, samq, nconsumers, infd=inf,
                                                      outfd=out.fileno())

    def test_samrecords_qnamegrouped(self):
        """
        Confirm that ValueError is raised if SAM file records are not
        grouped by qname.
        """
        samq = Queue(100)
        headerq = Queue(10)
        nconsumers = 1
        with (RESOURCES.joinpath('not_qnamegrouped.sam') as inf,
                NamedTemporaryFile() as out):
            with self.assertRaises(ValueError, msg=MSG_QNAMEGRP):
                samrecords(headerq, samq, nconsumers, infd=inf,
                                                      outfd=out.fileno())

    def test_samrecords_batch_and_enqueue(self):
        """
        Confirm that SAM file header is written headerq and SAM file records
        are written to samq in batched groups as expected.
        """
        samq = Queue(1000)
        headerq = Queue(1000)
        nconsumers = 2
        with (RESOURCES.joinpath('6_good_records.sam') as inf,
                NamedTemporaryFile() as out):
            samrecords(headerq, samq, nconsumers, infd=inf,
                                                  outfd=out.fileno())
        # One header per consumer.
        self.assertEqual(headerq.qsize(), nconsumers)
        # One batch (3 QNAME groups < batchsize) + one sentinel per consumer.
        self.assertEqual(samq.qsize(), 1 + nconsumers)
        # 3 QNAME groups in the batch, each with one pair.
        batch = samq.get()
        self.assertEqual(len(batch), 3)
        for group in batch:
            self.assertEqual(len(group), 2)

    def test_readends_1(self):
        """
        Confirm that calculated ends of a duplicate pair in same orientation
        are the same.
        """
        sam = AlignmentFile(RESOURCES.joinpath('2_pairs_same_orientation.sam'))
        records = list(iter(sam))
        pair_1 = (records[0], records[1])
        pair_2 = (records[2], records[3])
        self.assertTrue(pair_1[0].is_read1 and pair_1[0].is_forward)
        self.assertTrue(pair_1[1].is_read2 and pair_1[1].is_reverse)
        self.assertTrue(pair_2[0].is_read1 and pair_2[0].is_forward)
        self.assertTrue(pair_2[1].is_read2 and pair_2[1].is_reverse)
        ends_1 = readends(pair_1)
        ends_2 = readends(pair_2)
        self.assertEqual(ends_1, ends_2)

    def test_readends_2(self):
        """
        Confirm that calculated ends of a duplicate pair in opposite
        orientation are the same.
        """
        sam = AlignmentFile(RESOURCES.joinpath(
            '2_pairs_opposite_orientation.sam'))
        records = list(iter(sam))
        pair_1 = (records[0], records[1])
        pair_2 = (records[2], records[3])
        self.assertTrue(pair_1[0].is_read1 and pair_1[0].is_forward)
        self.assertTrue(pair_1[1].is_read2 and pair_1[1].is_reverse)
        self.assertTrue(pair_2[0].is_read1 and pair_2[0].is_reverse)
        self.assertTrue(pair_2[1].is_read2 and pair_2[1].is_forward)
        ends_1 = readends(pair_1)
        ends_2 = readends(pair_2)
        self.assertEqual(ends_1, ends_2)

    def test_readends_3(self):
        """
        Confirm that soft clipping at fwd and reverse ends is handled
        as expected.
        """
        sam = AlignmentFile(RESOURCES.joinpath('2_pairs_soft_clipping.sam'))
        records = list(iter(sam))
        pair_1 = (records[0], records[1])
        pair_2 = (records[2], records[3])
        self.assertTrue(pair_1[0].is_read1 and pair_1[0].is_reverse)
        self.assertTrue(pair_1[1].is_read2 and pair_1[1].is_forward)
        self.assertTrue(pair_2[0].is_read1 and pair_2[0].is_forward)
        self.assertTrue(pair_2[1].is_read2 and pair_2[1].is_reverse)
        # Alignments all have different pos.
        self.assertNotEqual(pair_1[0].pos, pair_2[0].pos)
        self.assertNotEqual(pair_1[0].pos, pair_2[1].pos)
        self.assertNotEqual(pair_1[1].pos, pair_2[0].pos)
        self.assertNotEqual(pair_1[1].pos, pair_2[1].pos)
        # But accounting for soft clipping the template ends are the same.
        ends_1 = readends(pair_1)
        ends_2 = readends(pair_2)
        self.assertEqual(ends_1, ends_2)

    def test_readends_4(self):
        """
        Confirm that one end unmapped is handled as expected.
        """
        sam = AlignmentFile(RESOURCES.joinpath('1_pair_one_unmapped_end.sam'))
        records = list(iter(sam))
        pair_1 = (records[0], records[1])
        self.assertTrue(pair_1[0].is_read1 and pair_1[0].is_unmapped)
        self.assertTrue(pair_1[1].is_read2 and pair_1[1].is_mapped)
        ends_1 = readends(pair_1)
        # Unmapped read end is sorted to last.
        self.assertEqual(ends_1[-1], UNMAPPED)

    def test_readends_5(self):
        """
        Confirm that both ends unmapped is handled as expected.
        """
        sam = AlignmentFile(RESOURCES.joinpath('1_pair_two_unmapped_ends.sam'))
        records = list(iter(sam))
        pair_1 = (records[0], records[1])
        self.assertTrue(pair_1[0].is_read1 and pair_1[0].is_unmapped)
        self.assertTrue(pair_1[1].is_read2 and pair_1[1].is_unmapped)
        ends_1 = readends(pair_1)
        # If both reads are unmapped readends returns None.
        self.assertIsNone(ends_1)

    def test_markdups_1(self):
        """
        Confirm that duplicates are marked as expected.
        """
        samq = Queue(1000)
        headerq = Queue(1000)
        nconsumers = 1
        expected = [
            (alignment.qname, alignment.flag) for alignment in
            AlignmentFile(RESOURCES.joinpath('test.qname.streammd.sam'))]
        with (SharedMemoryManager() as smm,
                RESOURCES.joinpath('test.qname.sam') as inf,
                NamedTemporaryFile() as out):
            samrecords(headerq, samq, nconsumers, infd=inf, outfd=out.fileno())
            bf = BloomFilter(smm, DEFAULT_NITEMS, DEFAULT_FPRATE)
            counts = markdups(bf.config, headerq, samq, outfd=out.fileno())
            n_qname, n_align, n_dup = counts
            self.assertEqual(n_qname, 2027)
            self.assertEqual(n_align, 4058)
            self.assertEqual(n_dup, 2037)
            result = [
                (alignment.qname, alignment.flag) for alignment in
                AlignmentFile(out.name)]
            self.assertEqual(result, expected)
            

    def test_main(self):
        """
        Confirm that main() operates as expected.
        """
