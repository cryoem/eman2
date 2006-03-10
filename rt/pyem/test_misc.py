#!/usr/bin/env python

from EMAN2 import *
import unittest
from test import test_support
from pyemtbx.exceptions import *
import testlib
import sys

class TestPixel(unittest.TestCase):
    """miscellaneous tests"""
    
    def test_pixel(self):
        """test Pixel class ................................."""
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
    """some boost related test .............................."""
    
    def test_defaultargs(self):
        """test default arguments ..........................."""
        imgfile1 = "test_defaultargs_1.mrc"
        imgfile2 = "test_defaultargs_2.mrc"
        TestUtil.make_image_file(imgfile1, MRC)
        TestUtil.make_image_file(imgfile2, MRC)
        
        image1 = EMData()
        image1.read_image(imgfile1)

        image2 = EMData()
        image2.read_image(imgfile2)

        image3 = image1.calc_ccf(image2)
        testlib.check_emdata(image3, sys.argv[0])

        os.unlink(imgfile1)
        os.unlink(imgfile2)

class TestException(unittest.TestCase):
    """exception related tests"""
    
    def test_FileAccessException(self):
        """test file access exception ......................."""
        e = EMData()
        e.set_size(10, 10, 1)

        try:
            e.read_image("__notexistingfile__.mrc")
        except RuntimeError, runtime_err:
            err_type = exception_type(runtime_err)
            self.assertEqual(err_type, "FileAccessException")

    def test_NotExistingObjectException(self):
        """test not existing object exception ..............."""
        e = EMData()
        e.set_size(100, 100, 1)

        try:
            e.process_inplace("eman1.NotExistintFilter_kfjda")
        except RuntimeError, runtime_err:
            err_type = exception_type(runtime_err)
            self.assertEqual(err_type, "NotExistingObjectException")

    def test_ImageFormatException(self):
        """test image format exception ......................"""
        fake_img = TestUtil.get_debug_image("fake.mrc")
        e = EMData()
        try:
            e.read_image(fake_img)
        except RuntimeError, runtime_err:
            err_type = exception_type(runtime_err)
            self.assertEqual(err_type, "FileAccessException")
            
class TestRegion(unittest.TestCase):
    """tests for class Region"""
    
    def test_region(self):
        """test region construction ........................."""
        e = Region()
        self.assertEqual(e.get_ndim(), 0)
        
        e2 = Region(5, 32)
        self.assertEqual(e2.get_ndim(), 1)
        
        e3 = Region(5, 6, 32, 32)
        self.assertEqual(e3.get_ndim(), 2)
        
        e4 = Region(5, 6, 7, 32, 32, 32)
        self.assertEqual(e4.get_ndim(), 3)
        
        e5 = Region(5.1, 32.2)
        self.assertEqual(e5.get_ndim(), 1)
        
        e6 = Region(5.2, 6.4, 32.5, 32.6)
        self.assertEqual(e6.get_ndim(), 2)
        
        e7 = Region(5.1, 6.2, 7.3, 32.4, 32.5, 32.6)
        self.assertEqual(e7.get_ndim(), 3)
        
        e8 = Region(e7)
        self.assertEqual(e8.get_ndim(), 3)
        

def test_main():
    test_support.run_unittest(TestPixel, TestBoost, TestException, TestRegion)

if __name__ == '__main__':
    test_main()


