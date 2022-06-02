from __future__ import print_function
from __future__ import division

from sphire.bin_py3 import sp_relion2sphire as oldfu
from ..sphire.bin import sp_relion2sphire as fu
import unittest

try:
    # python 3.4+ should use builtin unittest.mock not mock package
    from unittest.mock import patch
except ImportError:
    from mock import patch

try:
    from StringIO import StringIO  # python2 case
except ImportError:
    # python3 case. You will get an error because 'sys.stdout.write(msg)' presents in the library not in the test!!
    from io import StringIO
import sys


"""
WHAT IS MISSING:
1) 'prepare_outdir_log' is just crating the output directory for log stuff. I 'm not going to test it
2) cannot run from pyvharm because it cannot find 'e2proc2d.py'

RESULT AND KNOWN ISSUES
Some compatibility tests for the following functions fail!!!


In these tests there is a bug --> syntax error:


In these tests there is a strange behavior:

"""


class Test_helperFunctions(unittest.TestCase):
    def test_get_cmd_line(self):
        cmdLine=["this", "is", "a", "test"]
        with patch.object(sys, 'argv', cmdLine):
            return_new = fu.get_cmd_line()
            return_old = oldfu.get_cmd_line()
            self.assertEqual(return_new,return_old)
            self.assertEqual(return_new, 'Shell line command: this  is  a  test  ')

    def test_makerealpath(self):
        return_new = fu.makerelpath(p1="aa/vv",p2="aa2/vv2")
        return_old = oldfu.makerelpath(p1="aa/vv",p2="aa2/vv2")
        self.assertEqual(return_new, return_old)
        self.assertEqual(return_new, '../../aa2/vv2')

    def test_estimate_angle(self):
        return_new = fu.estimate_angle(coords_a=[15,32], coords_b=[18,38])
        return_old = oldfu.estimate_angle(coords_a=[15,32], coords_b=[18,38])
        self.assertEqual(return_new, return_old)
        self.assertEqual(return_new, -63.43494882292201)

    def test_estimate_angle_same_coords(self):
        return_new = fu.estimate_angle(coords_a=[15,32], coords_b=[15,32])
        return_old = oldfu.estimate_angle(coords_a=[15,32], coords_b=[15,32])
        self.assertEqual(return_new, return_old)
        self.assertEqual(return_new, -0.0)

    def test_estimate_angle_index_error_coords_a(self):
        with self.assertRaises(IndexError) as cm_new:
            fu.estimate_angle(coords_a=[], coords_b=[18,38])
        with self.assertRaises(IndexError) as cm_old:
            oldfu.estimate_angle(coords_a=[], coords_b=[18,38])
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_estimate_angle_index_error_coords_b(self):
        with self.assertRaises(IndexError) as cm_new:
            fu.estimate_angle(coords_b=[], coords_a=[18,38])
        with self.assertRaises(IndexError) as cm_old:
            oldfu.estimate_angle(coords_b=[], coords_a=[18,38])
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))
