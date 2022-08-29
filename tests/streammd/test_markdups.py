"""""
Test markdups module.
"""
from multiprocessing.managers import SharedMemoryManager
from multiprocessing import Manager, Queue
from importlib.resources import files
from tempfile import NamedTemporaryFile
from unittest import TestCase
from unittest.mock import patch
import contextlib
import io
import json
import sys

from pysam import AlignmentFile

from streammd.bloomfilter import BloomFilter
from streammd.markdups import (DEFAULT_FPRATE,
                               DEFAULT_NITEMS,
                               ALIGNMENTS,
                               ALIGNMENTS_MARKED_DUPLICATE,
                               READ_PAIRS,
                               READ_PAIRS_MARKED_DUPLICATE,
                               MSG_NOHEADER,
                               MSG_QNAMEGRP,
                               PGID,
                               UNMAPPED,
                               VERSION,
                               input_alnfile,
                               main,
                               markdups,
                               pgline,
                               readends)

RESOURCES = files('tests.streammd.resources')


class TestMarkDups(TestCase):
    """
    Test markdups module functions.
    """

    def test_input_alnfile_header(self):
        """
        Confirm that ValueError is raised if SAM file input lacks header.
        """
        nconsumers = 1
        headerq = Queue(1)
        inq = Queue(100)
        with (RESOURCES.joinpath('no_header.sam') as inf,
                NamedTemporaryFile() as out):
            with self.assertRaises(ValueError, msg=MSG_NOHEADER):
                input_alnfile(inf, out.fileno(), headerq, inq, nconsumers)

    def test_input_alnfile_qnamegrouped(self):
        """
        Confirm that ValueError is raised if SAM file records are not
        grouped by qname.
        """
        nconsumers = 1
        headerq = Queue(1)
        inq = Queue(100)
        with (RESOURCES.joinpath('not_qnamegrouped.sam') as inf,
                NamedTemporaryFile() as out):
            with self.assertRaises(ValueError, msg=MSG_QNAMEGRP):
                input_alnfile(inf, out.fileno(), headerq, inq, nconsumers)

    def test_input_alnfile_batch_and_enqueue(self):
        """
        Confirm that SAM file header is written headerq and SAM file records
        are written to samq in batched groups as expected.
        """
        nconsumers = 2
        headerq = Queue(1)
        inq = Queue(100)
        with (RESOURCES.joinpath('6_good_records.sam') as inf,
                NamedTemporaryFile() as out):
            input_alnfile(inf, out.fileno(), headerq, inq, nconsumers)
        # Header is set.
        self.assertTrue(headerq.qsize(), 1)
        # One batch (3 QNAME groups < batchsize) + one sentinel per consumer.
        self.assertEqual(inq.qsize(), 1 + nconsumers)
        # 3 QNAME groups in the batch, each with one pair.
        batch = inq.get()
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

    def test_markdups(self):
        """
        Confirm that duplicates are marked as expected.
        """
        nconsumers = 1
        headerq = Queue(1)
        inq = Queue(100)
        outq = Queue(nconsumers)
        expected = [
            (alignment.qname, alignment.flag) for alignment in
            AlignmentFile(RESOURCES.joinpath('test.qname.streammd.sam'))]
        with (SharedMemoryManager() as smm,
                RESOURCES.joinpath('test.qname.sam') as inf,
                NamedTemporaryFile() as out):
            outfd=out.fileno()
            input_alnfile(inf, outfd, headerq, inq, nconsumers)
            header = headerq.get()
            bf = BloomFilter(smm, DEFAULT_NITEMS, DEFAULT_FPRATE)
            markdups(bf.config, header, inq, outq, outfd)
            counts = outq.get()
            self.assertEqual(counts[READ_PAIRS], 2027)
            self.assertEqual(counts[READ_PAIRS_MARKED_DUPLICATE], 1018)
            self.assertEqual(counts[ALIGNMENTS], 4058)
            self.assertEqual(counts[ALIGNMENTS_MARKED_DUPLICATE], 2037)
            result = [
                (alignment.qname, alignment.flag) for alignment in
                AlignmentFile(out.name)]
            self.assertEqual(result, expected)

    def test_main_markdups_1(self):
        """
        Confirm that main() markdups operates as expected.
        """
        expected = [
            (alignment.qname, alignment.flag) for alignment in
            AlignmentFile(RESOURCES.joinpath('test.qname.streammd.sam'))]
        with (RESOURCES.joinpath('test.qname.sam') as inf,
                NamedTemporaryFile() as out):
            testargs = list(
                map(str, ('streammd', '--consumer-processes', 1, '--input',
                          inf, '--output', out.name, '--metrics',
                          '/dev/null')))
            with patch.object(sys, 'argv', testargs):
                main()
            result = [
                (alignment.qname, alignment.flag) for alignment in
                AlignmentFile(out.name)]
            self.assertEqual(result, expected)

    def test_main_markdups_2(self):
        """
        Confirm that main() markdups --metrics operates as expected.
        """
        expected = {
            'ALIGNMENTS': 4058,
            'ALIGNMENTS_MARKED_DUPLICATE': 2037,
            'READ_PAIRS': 2027,
            'READ_PAIRS_MARKED_DUPLICATE': 1018,
            'READ_PAIR_DUPLICATE_FRACTION': 0.5022,
            'UNIQUE_ITEMS_APPROXIMATE': 1011
        }
        with (RESOURCES.joinpath('test.qname.sam') as inf,
                NamedTemporaryFile() as out,
                NamedTemporaryFile() as met):
            testargs = list(
                map(str, ('streammd', '--consumer-processes', 1, '--input',
                          inf, '--output', out.name, '--metrics', met.name)))
            outstr = io.StringIO()
            with (patch.object(sys, 'argv', testargs),
                    self.assertLogs('streammd.markdups', level='INFO') as log):
                main()
                self.assertEqual(expected, json.load(met))

    def test_main_memcalc(self):
        """
        Confirm that main() --mem-calc operates as expected.
        """
        expected = '0.003GB\n'
        testargs = list(
            map(str, ('streammd', '--mem-calc', '1000000', '0.000001')))
        with patch.object(sys, 'argv', testargs):
            outstr = io.StringIO()
            with contextlib.redirect_stdout(outstr):
                with self.assertRaises(SystemExit):
                    main()
                self.assertEqual(outstr.getvalue(), expected)

    def test_pgline_1(self):
        """
        Confirm that pgline operates as expected when last line of the header
        is @SQ.
        """
        last = '@SQ\tSN:chrMT\tLN:16569'
        result = pgline(last)
        tkns = result.strip().split('\t')
        self.assertEqual(tkns[0], '@PG')
        tags = {tag: value for tag, value in
                (tkn.split(':', 1) for tkn in tkns[1:])}
        self.assertEqual(tags.get('ID'), PGID)
        self.assertEqual(tags.get('PN'), PGID)
        self.assertFalse(tags.get('PP'))
        self.assertEqual(tags.get('VN'), VERSION)

    def test_pgline_2(self):
        """
        Confirm that pgline operates as expected when last line of the header
        is @PG.
        """
        last = ('@PG\tID:bwa\tPN:bwa\tVN:0.7.15-r1140\t'
                'CL:bwa mem -t 6 -p -K 100000000 '
                '/reference/genomes/GRCh37_ICGC_standard_v2/indexes/'
                'BWAKIT_0.7.12/GRCh37_ICGC_standard_v2.fa -\n')
        result = pgline(last)
        tkns = result.strip().split('\t')
        self.assertEqual(tkns[0], '@PG')
        tags = {tag: value for tag, value in
                (tkn.split(':', 1) for tkn in tkns[1:])}
        self.assertEqual(tags.get('ID'), PGID)
        self.assertEqual(tags.get('PN'), PGID)
        self.assertEqual(tags.get('PP'), 'bwa')
        self.assertEqual(tags.get('VN'), VERSION)

    def test_unmapped(self):
        """
        Confirm that records with both read and mate unmapped appear in the
        output.
        """
        nconsumers = 1
        headerq = Queue(1)
        inq = Queue(100)
        outq = Queue(nconsumers)
        expected = [
            ('HWI-ST1213:151:C1DTBACXX:2:2207:13476:31678', 77),
            ('HWI-ST1213:151:C1DTBACXX:2:2207:13476:31678', 141)]
        with (SharedMemoryManager() as smm,
                RESOURCES.joinpath('test.unmapped.sam') as inf,
                NamedTemporaryFile() as out):
            outfd=out.fileno()
            input_alnfile(inf, outfd, headerq, inq, nconsumers)
            header = headerq.get()
            bf = BloomFilter(smm, DEFAULT_NITEMS, DEFAULT_FPRATE)
            markdups(bf.config, header, inq, outq, outfd)
            counts = outq.get()
            self.assertEqual(counts[READ_PAIRS], 1)
            self.assertEqual(counts[READ_PAIRS_MARKED_DUPLICATE], 0)
            self.assertEqual(counts[ALIGNMENTS], 2)
            self.assertEqual(counts[ALIGNMENTS_MARKED_DUPLICATE], 0)
            result = [
                (alignment.qname, alignment.flag) for alignment in
                AlignmentFile(out.name)]
            self.assertEqual(result, expected)
