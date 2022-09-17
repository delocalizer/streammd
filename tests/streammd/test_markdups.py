"""""
Test markdups module.
"""
from multiprocessing.managers import SharedMemoryManager
from multiprocessing import SimpleQueue
from importlib.resources import files
from tempfile import NamedTemporaryFile
from unittest import TestCase
from unittest.mock import patch
import io
import json
import sys

from pysam import AlignmentFile

from streammd.bloomfilter import BloomFilter
from streammd.markdups import (DEFAULT_FPRATE,
                               MSG_NOHEADER,
                               MSG_NOTSINGLE,
                               MSG_NOTPAIRED,
                               PGID,
                               SENTINEL,
                               UNMAPPED,
                               VERSION,
                               read_input,
                               main,
                               markdups,
                               pgline,
                               readends)

RESOURCES = files('tests.streammd.resources')


class TestMarkDups(TestCase):
    """
    Test markdups module functions.
    """

    def test_read_input_header(self):
        """
        Confirm that ValueError is raised if SAM file input lacks header.
        """
        nworkers = 1
        workq = SimpleQueue()
        inf = RESOURCES.joinpath('no_header.sam')
        with (self.assertRaises(ValueError, msg=MSG_NOHEADER),
                NamedTemporaryFile('wt') as out):
            read_input(inf, out.name, workq, nworkers)

    def test_read_input_batch_and_enqueue(self):
        """
        Confirm that SAM file records are written to workq in batched groups as
        expected.
        """
        nworkers = 2
        workq = SimpleQueue()
        inf = RESOURCES.joinpath('6_good_records.sam')
        with NamedTemporaryFile('wt') as out:
            read_input(inf, out.name, workq, nworkers)
        # One batch (3 QNAME groups < batchsize) + one sentinel per consumer.
        batch = workq.get()
        for _ in range(nworkers):
            self.assertEqual(workq.get(), SENTINEL)
        self.assertTrue(workq.empty())
        # 3 QNAME groups in the batch, each with one pair.
        self.assertEqual(len(batch), 3)
        for group in batch:
            self.assertEqual(len(group), 2)

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

    def test_readends_1(self):
        """
        Confirm that calculated ends of aligned single-end reads with same
        rname, pos and orientation are the same.
        """
        sam = AlignmentFile(RESOURCES.joinpath('2_reads_same_orientation.sam'))
        records = list(iter(sam))
        read_1, read_2 = records[0], records[1]
        self.assertEqual(read_1.rname, read_2.rname)
        self.assertEqual(read_1.pos, read_2.pos)
        self.assertEqual(read_1.is_forward,  read_2.is_forward)
        read_1s, read_2s = [r.to_string().split('\t')
                            for r in (read_1, read_2)]
        ends_1 = readends([read_1s])
        ends_2 = readends([read_2s])
        self.assertEqual(ends_1, ends_2)

    def test_readends_2(self):
        """
        Confirm that calculated ends of aligned paired-end reads with same
        rname and pos but opposite orientation are different.
        """
        sam = AlignmentFile(RESOURCES.joinpath('2_reads_diff_orientation.sam'))
        records = list(iter(sam))
        read_1, read_2 = records[0], records[1]
        self.assertEqual(read_1.rname, read_2.rname)
        self.assertEqual(read_1.pos, read_2.pos)
        self.assertNotEqual(read_1.is_forward,  read_2.is_forward)
        read_1s, read_2s = [r.to_string().split('\t')
                            for r in (read_1, read_2)]
        ends_1 = readends([read_1s])
        ends_2 = readends([read_2s])
        self.assertNotEqual(ends_1, ends_2)

    def test_readends_3(self):
        """
        Confirm that calculated ends of a duplicate pair in same orientation
        are the same.
        """
        sam = AlignmentFile(RESOURCES.joinpath('2_pairs_same_orientation.sam'))
        records = list(iter(sam))
        pair_1, pair_2 = (records[0], records[1]), (records[2], records[3])
        self.assertTrue(pair_1[0].is_read1 and pair_1[0].is_forward)
        self.assertTrue(pair_1[1].is_read2 and pair_1[1].is_reverse)
        self.assertTrue(pair_2[0].is_read1 and pair_2[0].is_forward)
        self.assertTrue(pair_2[1].is_read2 and pair_2[1].is_reverse)
        pair_1s, pair_2s = [(r1.to_string().split('\t'),
                             r2.to_string().split('\t'))
                            for (r1, r2) in (pair_1, pair_2)]
        ends_1 = readends(pair_1s)
        ends_2 = readends(pair_2s)
        self.assertEqual(ends_1, ends_2)

    def test_readends_4(self):
        """
        Confirm that calculated ends of a duplicate pair in opposite
        orientation are the same.
        """
        sam = AlignmentFile(RESOURCES.joinpath(
            '2_pairs_opposite_orientation.sam'))
        records = list(iter(sam))
        pair_1, pair_2 = (records[0], records[1]), (records[2], records[3])
        self.assertTrue(pair_1[0].is_read1 and pair_1[0].is_forward)
        self.assertTrue(pair_1[1].is_read2 and pair_1[1].is_reverse)
        self.assertTrue(pair_2[0].is_read1 and pair_2[0].is_reverse)
        self.assertTrue(pair_2[1].is_read2 and pair_2[1].is_forward)
        pair_1s, pair_2s = [(r1.to_string().split('\t'),
                             r2.to_string().split('\t'))
                            for (r1, r2) in (pair_1, pair_2)]
        ends_1 = readends(pair_1s)
        ends_2 = readends(pair_2s)
        self.assertEqual(ends_1, ends_2)

    def test_readends_5(self):
        """
        Confirm that soft clipping at fwd and reverse ends is handled
        as expected.
        """
        sam = AlignmentFile(RESOURCES.joinpath('2_pairs_soft_clipping.sam'))
        records = list(iter(sam))
        pair_1, pair_2 = (records[0], records[1]), (records[2], records[3])
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
        pair_1s, pair_2s = [(r1.to_string().split('\t'),
                             r2.to_string().split('\t'))
                            for (r1, r2) in (pair_1, pair_2)]
        ends_1 = readends(pair_1s)
        ends_2 = readends(pair_2s)
        self.assertEqual(ends_1, ends_2)

    def test_readends_6(self):
        """
        Confirm that pair with one end unmapped is handled as expected.
        """
        sam = AlignmentFile(RESOURCES.joinpath('1_pair_one_unmapped_end.sam'))
        records = list(iter(sam))
        pair = (records[0], records[1])
        self.assertTrue(pair[0].is_read1 and pair[0].is_unmapped)
        self.assertTrue(pair[1].is_read2 and pair[1].is_mapped)
        pair_s = [r.to_string().split('\t') for r in pair]
        ends = readends(pair_s)
        # Unmapped read was first but unmapped end is sorted to last.
        self.assertEqual(ends[-1], UNMAPPED)

    def test_readends_7(self):
        """
        Confirm that pair with both ends unmapped is handled as expected.
        """
        sam = AlignmentFile(RESOURCES.joinpath('1_pair_two_unmapped_ends.sam'))
        records = list(iter(sam))
        pair = (records[0], records[1])
        self.assertTrue(pair[0].is_read1 and pair[0].is_unmapped)
        self.assertTrue(pair[1].is_read2 and pair[1].is_unmapped)
        pair_s = [r.to_string().split('\t') for r in pair]
        ends = readends(pair_s)
        self.assertEqual(ends, [UNMAPPED, UNMAPPED])

    def test_markdups_notpaired(self):
        """
        Confirm that ValueError is raised if reads_per_template == 2 and
        SAM file records are not grouped by qname.
        """
        nworkers = 1
        workq = SimpleQueue()
        resultq = SimpleQueue()
        inf = RESOURCES.joinpath('not_paired.sam')
        with (SharedMemoryManager() as smm,
                NamedTemporaryFile('wt') as out):
            read_input(inf, out.name, workq, nworkers)
            bf = BloomFilter(smm, 100, DEFAULT_FPRATE)
            with self.assertRaises(ValueError, msg=MSG_NOTPAIRED):
                markdups(bf.config, workq, out.name, resultq, 2)

    def test_markdups_notsingle(self):
        """
        Confirm that ValueError is raised if reads_per_template == 1 and
        SAM file records contain multiple primary alignments per template.
        """
        nworkers = 1
        workq = SimpleQueue()
        resultq = SimpleQueue()
        inf = RESOURCES.joinpath('not_single.sam')
        with (SharedMemoryManager() as smm,
                NamedTemporaryFile('wt') as out):
            read_input(inf, out.name, workq, nworkers)
            bf = BloomFilter(smm, 100, DEFAULT_FPRATE)
            with self.assertRaises(ValueError, msg=MSG_NOTSINGLE):
                markdups(bf.config, workq, out.name, resultq, 1)

    def test_markdups_3(self):
        """
        Confirm that duplicates are marked as expected on single-end reads.
        """
        nworkers = 1
        workq = SimpleQueue()
        resultq = SimpleQueue()
        expected = [
            (alignment.qname, alignment.flag) for alignment in
            AlignmentFile(RESOURCES.joinpath('test.single.streammd.sam'))]
        inf = RESOURCES.joinpath('test.single.sam')
        with (SharedMemoryManager() as smm,
                NamedTemporaryFile('wt') as out):
            read_input(inf, out.name, workq, nworkers)
            bf = BloomFilter(smm, 100, DEFAULT_FPRATE)
            markdups(bf.config, workq, out.name, resultq, 1)
            counts = resultq.get()
            self.assertEqual(counts['TEMPLATES'], 4)
            self.assertEqual(counts['TEMPLATES_MARKED_DUPLICATE'], 1)
            self.assertEqual(counts['ALIGNMENTS'], 4)
            self.assertEqual(counts['ALIGNMENTS_MARKED_DUPLICATE'], 1)
            result = [(alignment.qname, alignment.flag) for alignment in
                    AlignmentFile(out.name)]
            self.assertEqual(result, expected)

    def test_markdups_4(self):
        """
        Confirm that duplicates are marked as expected on paired-end reads.
        """
        nworkers = 1
        workq = SimpleQueue()
        resultq = SimpleQueue()
        expected = [
            (alignment.qname, alignment.flag) for alignment in
            AlignmentFile(RESOURCES.joinpath('test.paired.streammd.sam'))]
        inf = RESOURCES.joinpath('test.paired.sam')
        with (SharedMemoryManager() as smm,
                NamedTemporaryFile('wt') as out):
            read_input(inf, out.name, workq, nworkers)
            bf = BloomFilter(smm, 100, DEFAULT_FPRATE)
            markdups(bf.config, workq, out.name, resultq, 2)
            counts = resultq.get()
            self.assertEqual(counts['TEMPLATES'], 2)
            self.assertEqual(counts['TEMPLATES_MARKED_DUPLICATE'], 1)
            self.assertEqual(counts['ALIGNMENTS'], 4)
            self.assertEqual(counts['ALIGNMENTS_MARKED_DUPLICATE'], 2)
            result = [
                (alignment.qname, alignment.flag) for alignment in
                AlignmentFile(out.name)]
            self.assertEqual(result, expected)

    def test_markdups_5(self):
        """
        Confirm that paired records with both read and mate unmapped still
        appear in the output.
        """
        nworkers = 1
        workq = SimpleQueue()
        resultq = SimpleQueue()
        expected = [
            ('HWI-ST1213:151:C1DTBACXX:2:2207:13476:31678', 77),
            ('HWI-ST1213:151:C1DTBACXX:2:2207:13476:31678', 141)]
        inf = RESOURCES.joinpath('test.unmapped.sam')
        with (SharedMemoryManager() as smm,
                NamedTemporaryFile('wt') as out):
            read_input(inf, out.name, workq, nworkers)
            bf = BloomFilter(smm, 100, DEFAULT_FPRATE)
            markdups(bf.config, workq, out.name, resultq, 2)
            counts = resultq.get()
            self.assertEqual(counts['TEMPLATES'], 1)
            self.assertEqual(counts['TEMPLATES_MARKED_DUPLICATE'], 0)
            self.assertEqual(counts['ALIGNMENTS'], 2)
            self.assertEqual(counts['ALIGNMENTS_MARKED_DUPLICATE'], 0)
            result = [
                (alignment.qname, alignment.flag) for alignment in
                AlignmentFile(out.name)]
            self.assertEqual(result, expected)

    def test_markdups_6(self):
        """
        Confirm that when strip_previous is False (default) that previously
        marked duplicate flag remains on read no longer considered duplicate.
        """
        nworkers = 1
        workq = SimpleQueue()
        resultq = SimpleQueue()
        inf = RESOURCES.joinpath('test.previousdup.sam')
        with (SharedMemoryManager() as smm,
                NamedTemporaryFile('wt') as out):
            read_input(inf, out.name, workq, nworkers)
            bf = BloomFilter(smm, 100, DEFAULT_FPRATE)
            markdups(bf.config, workq, out.name, resultq, 1)
            result = [(alignment.flag, alignment.get_tag('PG'))
                      for alignment in AlignmentFile(out.name)]
            self.assertEqual(result[0], (1024, 'MarkDuplicates'))  # previous
            self.assertEqual(result[1], (1024, PGID))              # new

    def test_markdups_7(self):
        """
        Confirm that when strip_previous is True that previously marked
        duplicate flag is removed on read no longer considered duplicate.
        """
        nworkers = 1
        workq = SimpleQueue()
        resultq = SimpleQueue()
        inf = RESOURCES.joinpath('test.previousdup.sam')
        with (SharedMemoryManager() as smm,
                NamedTemporaryFile('wt') as out):
            read_input(inf, out.name, workq, nworkers)
            bf = BloomFilter(smm, 100, DEFAULT_FPRATE)
            markdups(bf.config, workq, out.name, resultq, 1, True)
            result = [(alignment.flag, alignment.get_tag('PG'))
                      for alignment in AlignmentFile(out.name)]
            self.assertEqual(result[0], (0, PGID))     # updated
            self.assertEqual(result[1], (1024, PGID))  # new

    def test_main_markdups_1(self):
        """
        Confirm that markdups.main() operates as expected on a larger input
        file: 4058 alignments from 2027 templates, with a variety of different
        flags, orientations, and complex CIGAR strings.

        The test reference output was generated from the test input by Picard
        MarkDuplicates (2.23.8).

        The input SAM records in 'test.paired_full.sam' have been ordered with
        higher quality reads in a duplicate set occuring first so that in this
        case we expect the output of streammd to EXACTLY match that from
        Picard MarkDuplicates (assuming no false positives from streammd,
        which happens to be the case for this data with default --fp-rate and
        seeds).

        In the general case we expect only that the duplicate counts should
        match, since Picard MarkDuplicates picks the highest quality read in a
        set as the original, where streammd must pick the first.
        """
        expected = [
            (alignment.qname, alignment.flag) for alignment in
            AlignmentFile(RESOURCES.joinpath('test.paired_full.picardmd.sam'))]
        with (RESOURCES.joinpath('test.paired_full.sam') as inf,
                NamedTemporaryFile() as out):
            testargs = list(
                map(str, ('streammd', '--workers', 1, '--input', inf,
                    '--output', out.name, '--metrics', '/dev/null',
                    '-n', 1000000, '-m', '4MiB')))
            with patch.object(sys, 'argv', testargs):
                main()
            result = [
                (alignment.qname, alignment.flag) for alignment in
                AlignmentFile(out.name)]
            self.assertEqual(result, expected)

    def test_main_markdups_2(self):
        """
        Confirm that markdups.main() --metrics operates as expected.
        """
        expected = {
            'ALIGNMENTS': 4058,
            'ALIGNMENTS_MARKED_DUPLICATE': 2039,
            'TEMPLATES': 2027,
            'TEMPLATES_MARKED_DUPLICATE': 1018,
            'TEMPLATE_DUPLICATE_FRACTION': 0.5022
        }
        with (RESOURCES.joinpath('test.paired_full.sam') as inf,
                NamedTemporaryFile() as out,
                NamedTemporaryFile() as met):
            testargs = list(
                map(str, ('streammd', '--workers', 1, '--input', inf,
                    '--output', out.name, '--metrics', met.name,
                    '-n', 1000000)))
            outstr = io.StringIO()
            with (patch.object(sys, 'argv', testargs),
                    self.assertLogs('streammd.markdups', level='INFO') as log):
                main()
                self.assertEqual(expected, json.load(met))
