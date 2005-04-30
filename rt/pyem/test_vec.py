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
        m1 = Transform3D()
        v1 = Vec3f(1.5, 2.5, 3.5)
        v2 = v1 * m1
        self.assertEqual(v2, v1)


    def test_vec_funcs(self):
        a1 = Vec3i()
        self.assertEqual(a1, Vec3i(0,0,0))

        a1 = Vec3i(1, 2, 4)
        self.assertEqual(a1, Vec3i(1, 2, 4))

        a2 = Vec3i(3,4,0)
        self.assertEqual(a2.length(), 5)

        n2 = a2.normalize()
        self.assertEqual(a2.length(), 0)
        self.assertEqual(int(n2), 5) 

        a3 = Vec3i(1, 2, 3)
        self.assertEqual(a3, Vec3i(1,2,3))

        a4 = Vec3f(3.0, 4.0, 0)
        a4.normalize()
        self.assertEqual(a4.length(), 1)

        a5 = Vec3i(1,2,3)
        a6 = Vec3i(2,4,5)
        dot_result = a5.dot(a6)
        self.assertEqual(dot_result, 25)

        a7 = Vec3f(a6)
        a8 = a7.cross(Vec3f(a5))
        a9 = a7.cross(Vec3f(a5))

        self.assertEqual(a8, a9)
        self.assertEqual(a9, Vec3f(2,-1,0))

        a8.set_value(1.1, 2.2, 3)

        self.assertEqual(a8, Vec3f(1.1, 2.2, 3))

        #a8[2] = 12
        #a8[1] += 2.2
        #self.assertEqual(a8, Vec3f(1.1, 4.4, 12))


    def test_vec_ops(self):
        v1 = Vec3f(1, 2, 3)
        v2 = Vec3f(2, 4, 6)
        v3 = Vec3f(v1)
       
        v1 += v2
        self.assertEqual(v1.as_list(), [3,6,9])

        v4 = v2 + v3
        self.assertEqual(v4, v1)
        
        v3 -= v2
        v5 = v3

        self.assertEqual(v3.as_list(), [-1,-2,-3])
        v6 = v5 - v2

        self.assertEqual(v3, Vec3f(-1,-2,-3))

        v3 *= 3
        self.assertEqual(v6, v3)
        self.assert_(v3 != v2)

        v3 *= 0
        v7 = v2 * 1
        self.assertEqual(v3, Vec3f())
        self.assertEqual(v7, v2)

        v8 = Vec3i(1,2,3)
        v8 *= 3
        self.assertEqual(v8, Vec3i(3,6,9))

        v9 = v8 * 3
        self.assertEqual(v9, Vec3i(9, 18, 27))

        v9 /= 3
        self.assertEqual(v9, v8)
        v9 = v8 / 3

        self.assertEqual(v9, Vec3i(1,2,3))



def test_main():
    test_support.run_unittest(TestVec)

if __name__ == '__main__':
    test_main()



