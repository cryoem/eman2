#!/bin/env python

from EMAN2 import *
import unittest
from test import test_support

class TestVec(unittest.TestCase):
	
	def test_as_list(self):
		v = Vec3i(1,2,3)
		vlist = v.as_list()
		self.assertEqual(vlist, [1,2,3])

	def test_cplusplus_vec_to_python_list(self):
		e = EMData()
		e.set_size(32,32,1)
		intp = e.calc_min_location()
		self.assertEqual(intp, (0,0,0))
		self.assertEqual(type(intp), type((1,2)))
		
	def test_python_list_to_cplusplus_vec(self):
		e = EMData()
		e.set_size(32,32,1)
		
		e2 = EMData()
		e2.set_size(12,12,1)
		e.insert_clip(e2, [1,1,0])
		e.insert_clip(e2, (1,1,0))

		e2.translate(Vec3f(1,2,3))
		e2.translate((1,2,3))
		e2.translate([1,2,3])

	def test_vec_transform_op(self):
		m1 = Transform()
		v1 = Vec3f(1.5, 2.5, 3.5)
		v2 = v1 * m1
		self.assertEqual(v2, v1)

	def test_vec_ops(self):		
		v1 = Vec3f(1,2,3)
		v2 = Vec3f(2,4,6)

		v1 = v1 + v2
		v1list = v1.as_list()
		assert(v1list == [3, 6, 9])




def test_main():
    test_support.run_unittest(TestVec)

if __name__ == '__main__':
    test_main()



