#!/bin/env python

from EMAN2 import *
import unittest
from test import test_support
from pyemtbx.exceptions import *
import testlib
import sys

class TestPixel(unittest.TestCase):
    
    def test_pixel(self):
        x = 1
        y = 2
        z = 3
        v = 2.5

        p1 = Pixel(x, y, z, v)
        p2 = Pixel(x, y, z, v)
        self.assertEqual(p1, p2)

        p3 = Pixel(1, 2, 3, 2.3)
        p4 = Pixel(11, 2, 3, 2.3)
        self.assert_(p3 != p4)

        self.assertEqual(p1.x, x)
        self.assertEqual(p1.y, y)
        self.assertEqual(p1.z, z)
        self.assertEqual(p1.value, v)



class TestBoost(unittest.TestCase):
    
    def test_overloads(self):
        e = EMData()
        e.set_size(10,10,1)
        # The following core dump
        #e.make_rotational_footprint()

    def test_defaultargs(self):
        image1 = EMData()
        image1.read_image(TestUtil.get_debug_image("samesize1.mrc"))

        image2 = EMData()
        image2.read_image(TestUtil.get_debug_image("samesize2.mrc"))

        #image3 = image1.calc_ccf(image2)
        #testlib.check_emdata(image3, sys.argv[0])




class TestException(unittest.TestCase):
    
    def test_FileAccessException(self):
        e = EMData()
        e.set_size(10, 10, 1)

        try:
            e.read_image("__notexistingfile__.mrc")
        except RuntimeError, runtime_err:
            err_type = exception_type(runtime_err)
            self.assertEqual(err_type, "FileAccessException")
            

    def test_NotExistingObjectException(self):
        e = EMData()
        e.set_size(100, 100, 1)

        try:
            e.filter("NotExistintFilter_kfjda")
        except RuntimeError, runtime_err:
            err_type = exception_type(runtime_err)
            self.assertEqual(err_type, "NotExistingObjectException")

    def test_ImageFormatException(self):
        fake_img = TestUtil.get_debug_image("fake.mrc")
        e = EMData()
        try:
            e.read_image(fake_img)
        except RuntimeError, runtime_err:
            err_type = exception_type(runtime_err)
            self.assertEqual(err_type, "ImageFormatException")

"""
class TestEMObject(unittest.TestCase):

    def test_basic_functionality(self):
        n = 2
        f = 3.5

        e1 = EMObject(n)
        e2 = EMObject(n)

        self.assertEqual(e1, e2)

        e3 = EMObject(f)
        e4 = EMObject(f)

        self.assertEqual(e3, e4)

        self.assertEqual(int(e1), n)
        self.assertEqual(float(e3), f)

        n2 = 2.0
        self.assertEqual((float(EMObject(float(n2)))), n2)

        str1 = "hello"
        self.assertEqual(str1, EMObject(str1).to_str())

        nx = 10
        ny = 12
        nz = 2
        emdata = EMData()
        emdata.set_size(nx, ny, nz)
        emobj = EMObject(emdata)
        emdata2 = emobj.to_EMAN_EMData()
        self.assertEqual(emdata2.get_xsize(), nx)
        self.assertEqual(emdata2.get_ysize(), ny)
        self.assertEqual(emdata2.get_zsize(), nz)
"""        


def test_main():
    test_support.run_unittest(TestPixel, TestBoost, TestException)

if __name__ == '__main__':
    test_main()


