#!/bin/env python

from EMAN2 import *
import unittest
from test import test_support
import math

class TestTransform(unittest.TestCase):

	def setUp(self):
		self.delta =  0.0001

	def assertfloat(self, f1, f2):
        diff = f1 - f2
        if math.fabs(diff) > self.delta:
            self.assertEqual(f1, f2)
		
	def test_get_rotation(self):
		az = -0.60170830102
		alt = 1.45232928554
		phi = 0

		t = Transform(EULER_EMAN, alt, az, phi)
		rot = t.get_rotation(EULER_EMAN)

		self.assertfloat(az, float(rot["az"]))
		self.assertfloat(alt, float(rot["alt"]))
		self.assertfloat(phi, float(rot["phi"]))


def test_main():
    test_support.run_unittest(TestTransform)

if __name__ == '__main__':
    test_main()



