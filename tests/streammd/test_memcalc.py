"""""
Test memcalc module.
"""
from unittest import TestCase
from unittest.mock import patch
import contextlib
import io
import sys

from streammd.memcalc import main
from streammd.bloomfilter import KMAX, MSG_NOMEM

class TestMemCalc(TestCase):
    """
    Test memcalc module functions.
    """

    def test_main_memcalc_min(self):
        """
        Confirm that memcalc.main() prints mem; k when --mem is not specified.
        """
        expected = 'mem=3.59 MB; k=20\n'
        testargs = list(
            map(str, ('memcalc', '1000000', '0.000001')))
        with patch.object(sys, 'argv', testargs):
            outstr = io.StringIO()
            with contextlib.redirect_stdout(outstr):
                main()
                self.assertEqual(outstr.getvalue(), expected)

    def test_main_memcalc_mem(self):
        """
        Confirm that memcalc.main() prints mem; k when --mem is specified.
        """
        expected = 'mem=10MB; k=5\n'
        testargs = list(
            map(str, ('memcalc', '1000000', '0.000001', '10MB')))
        with patch.object(sys, 'argv', testargs):
            outstr = io.StringIO()
            with contextlib.redirect_stdout(outstr):
                main()
                self.assertEqual(outstr.getvalue(), expected)

    def test_main_memcalc_nosoln(self):
        """
        Confirm that memcalc.main() prints warning to stderr when --mem is
        specified but is too low for a solution.
        """
        expected = MSG_NOMEM % ('1MB', 1000000, 0.000001, KMAX)
        testargs = list(
            map(str, ('memcalc', '1000000', '0.000001', '1MB')))
        with patch.object(sys, 'argv', testargs):
            outstr = io.StringIO()
            with (self.assertRaises(SystemExit),
                    contextlib.redirect_stderr(outstr)):
                main()
                self.assertEqual(outstr.getvalue(), expected)
