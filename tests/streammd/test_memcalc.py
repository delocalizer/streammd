"""""
Test memcalc module.
"""
from unittest import TestCase
from unittest.mock import patch
import contextlib
import io
import sys

from streammd.memcalc import main

class TestMemCalc(TestCase):
    """
    Test memcalc module functions.
    """

    def test_main_memcalc(self):
        """
        Confirm that memcalc.main() operates as expected.
        """
        expected = '0.003GB\n'
        testargs = list(
            map(str, ('memcalc', '1000000', '0.000001')))
        with patch.object(sys, 'argv', testargs):
            outstr = io.StringIO()
            with contextlib.redirect_stdout(outstr):
                main()
                self.assertEqual(outstr.getvalue(), expected)
