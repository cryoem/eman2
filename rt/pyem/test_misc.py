#!/usr/bin/env python

#
# Author: Liwei Peng, 01/30/2005 (sludtke@bcm.edu)
# Copyright (c) 2000-2006 Baylor College of Medicine
#
# This software is issued under a joint BSD/GNU license. You may use the
# source code in this file under either license. However, note that the
# complete EMAN2 and SPARX software packages have some GPL dependencies,
# so you are responsible for compliance with the licenses of these packages
# if you opt to use BSD licensing. The warranty disclaimer below holds
# in either instance.
#
# This complete copyright notice must be included in any revised version of the
# source code. Additional authorship citations may be added, but existing
# author citations must be preserved.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  2111-1307 USA
#
#

from EMAN2 import *
import unittest
from pyemtbx.exceptions import *
import testlib
import sys
from optparse import OptionParser

IS_TEST_EXCEPTION = False

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
        TestUtil.make_image_file(imgfile1, IMAGE_MRC)
        TestUtil.make_image_file(imgfile2, IMAGE_MRC)
        
        image1 = EMData()
        image1.read_image(imgfile1)

        image2 = EMData()
        image2.read_image(imgfile2)

        image3 = image1.calc_ccf(image2)
        testlib.check_emdata(image3, sys.argv[0])

        testlib.safe_unlink(imgfile1)
        testlib.safe_unlink(imgfile2)

class TestException(unittest.TestCase):
    """exception related tests"""
    
    def test_FileAccessException(self):
        """test file access exception ......................."""
        e = EMData()
        e.set_size(10, 10, 1)
        
        if(IS_TEST_EXCEPTION):
            try:
                e.read_image("__notexistingfile__.mrc")
            except RuntimeError, runtime_err:
                err_type = exception_type(runtime_err)
                self.assertEqual(err_type, "FileAccessException")

    def test_NotExistingObjectException(self):
        """test not existing object exception ..............."""
        e = EMData()
        e.set_size(100, 100, 1)
        
        if(IS_TEST_EXCEPTION):
            try:
                e.process_inplace("NotExistintFilter_kfjda")
            except RuntimeError, runtime_err:
                err_type = exception_type(runtime_err)
                self.assertEqual(err_type, "NotExistingObjectException")

    def test_ImageFormatException(self):
        """test image format exception ......................"""
        fake_img = TestUtil.get_debug_image("fake.mrc")
        e = EMData()
        
        if(IS_TEST_EXCEPTION):
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
    p = OptionParser()
    p.add_option('--t', action='store_true', help='test exception', default=False )
    global IS_TEST_EXCEPTION
    opt, args = p.parse_args()
    if opt.t:
        IS_TEST_EXCEPTION = True
    Log.logger().set_level(-1)  #perfect solution for quenching the Log error information, thank Liwei
    suite1 = unittest.TestLoader().loadTestsFromTestCase(TestPixel)
    suite2 = unittest.TestLoader().loadTestsFromTestCase(TestBoost)
    suite3 = unittest.TestLoader().loadTestsFromTestCase(TestException)
    suite4 = unittest.TestLoader().loadTestsFromTestCase(TestRegion)
    unittest.TextTestRunner(verbosity=2).run(suite1)
    unittest.TextTestRunner(verbosity=2).run(suite2)
    unittest.TextTestRunner(verbosity=2).run(suite3)
    unittest.TextTestRunner(verbosity=2).run(suite4)

if __name__ == '__main__':
    test_main()
