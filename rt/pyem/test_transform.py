#!/bin/env python

from EMAN2 import *
import unittest
from test import test_support
import testlib

class TestTransform(unittest.TestCase):
        
    def test_get_rotation(self):
        alt = 1.45232928554
        az = -0.60170830102
        phi = 0

        t = Transform3D(az, alt, phi)
        rot = t.get_rotation(EULER_EMAN)

        testlib.assertfloat(self, az, float(rot["az"]))
        testlib.assertfloat(self, alt, float(rot["alt"]))
        testlib.assertfloat(self, phi, float(rot["phi"]))


def test_main():
    test_support.run_unittest(TestTransform)

if __name__ == '__main__':
    test_main()



