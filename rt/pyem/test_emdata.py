#!/usr/bin/env python

import EMAN2
from EMAN2 import *
from pyemtbx.exceptions import *
import unittest
from test import test_support
import testlib
import sys
import math
import os

class TestEMData(unittest.TestCase):
    """this is the unit test for EMData class"""

    def test_default_args(self):
        """test default constructor of EMData ..............."""
        e = EMData()
        self.assertEqual(e.get_attr("apix_x"), 1.0)
        self.assertEqual(e.get_attr("apix_y"), 1.0)
        self.assertEqual(e.get_attr("apix_z"), 1.0)
        self.assertEqual(e.get_attr("is_complex"), 0)
        self.assertEqual(e.get_attr("is_complex_ri"), 1)
        self.assertEqual(e.get_xsize(), 0)
        self.assertEqual(e.get_ysize(), 0)
        self.assertEqual(e.get_zsize(), 0) 

    def test_copy(self):
        """test copy() function ............................."""
        e = EMData()
        e.set_size(32,32,32)
        e.to_zero()
        e.process_inplace("testimage.noise.uniform.rand")
        e2 = e.copy()
        
        self.assertEqual(e.get_attr_dict(), e2.get_attr_dict())
        for k in range(32):
            for j in range(32):
                for i in range(32):
                    self.assertEqual(e.get_3dview()[i][j][k], e2.get_3dview()[i][j][k])

    def test_copy_head(self):
        """test copy_head() function ........................"""
        e = EMData()
        e.set_size(32,32,32)
        e.to_zero()
        e.process_inplace("testimage.noise.uniform.rand")
        e2 = e.copy_head()
        
        dict1 = e.get_attr_dict()
        dict2 = e2.get_attr_dict()
        self.assertEqual(dict1['nx'], dict2['nx'])
        self.assertEqual(dict1['ny'], dict2['ny'])
        self.assertEqual(dict1['nz'], dict2['nz'])
        self.assertEqual(dict1.keys(), dict2.keys())
    
    def test_get_clip(self):
        """test get_clip() function ........................."""
        e = EMData()
        e.set_size(32,32,32)
        e.to_zero()
        e.process_inplace("testimage.noise.uniform.rand")
        
        #test inclusive clip
        e2 = e.get_clip(Region(10,10,10, 20,20,20))
        d1 = e.get_3dview()
        d2 = e2.get_3dview()
        for k in range(10):
            for j in range(10):
                for i in range(10):
                    self.assertEqual(d1[i+10][j+10][k+10], d2[i][j][k])
        
        #test padding zero for clip larger than original image
        e3 = e.get_clip(Region(30,30,30, 35,35,35))
        d3 = e3.get_3d_view()
        for k in range(2):
            for j in range(2):
                for i in range(2):
                    self.assertEqual(d1[i+30][j+30][k+30], d3[i][j][k])
        for k in range(3):
            for j in range(3):
                for i in range(3):
                    self.assertEqual(d3[i+2][j+2][k+2], 0)

    def test_insert_clip(self):
        """test insert_clip() function ......................"""
        e = EMData()
        e.set_size(32,32,32)
        e.to_zero()
        e.process_inplace("testimage.noise.uniform.rand")
        
        e2 = EMData()
        e2.set_size(10,10,10)
        e2.to_zero()
        e2.process_inplace("testimage.noise.uniform.rand")
        
        #test exception if the clip is out side of the image
        self.assertRaises( RuntimeError, e.insert_clip, e2, (30,30,30) )
        try:
            e.insert_clip( e2,(30,30,30) )
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "ImageFormatException")
        
        e.insert_clip(e2, (16,16,16))
        d1 = e.get_3dview()
        d2 = e2.get_3dview()
        for k in range(10):
            for j in range(10):
                for i in range(10):
                    self.assertEqual(d1[i+16][j+16][k+16], d2[i][j][k])

    def test_get_top_half(self):
        """test get_top_half() function ....................."""
        e = EMData()
        e.set_size(32,32,32)
        e.to_zero()
        e.process_inplace("testimage.noise.uniform.rand")
        e2 = e.get_top_half()
        d = e.get_3dview()
        d2 = e2.get_3dview()
        for k in range(16):
            for j in range(32):
                for i in range(32):
                    self.assertEqual(d[k+16][j][i], d2[k][j][i]) #(nz/2, nz] is the top half
        
        e3 = EMData()
        e3.set_size(32,32,1)
        e3.to_one()
        #get_top_half() should raise exception when apply to 2D image
        self.assertRaises( RuntimeError, e3.get_top_half, )
        try:
            e3.get_top_half()
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "ImageDimensionException")

    def test_insert_scaled_sum(self):
        """test insert_scaled_sum() function ................"""
        e = EMData()
        e.set_size(64,64,64)
        e.process_inplace("testimage.noise.uniform.rand")
        e2 = EMData()
        e2.set_size(32,32,32)
        e2.to_one()
        e.insert_scaled_sum(e2, (32,32,32))
        
        e3 = EMData()
        e3.set_size(32,1,1)
        e3.to_zero()
        e4 = EMData()
        e4.set_size(12,1,1)
        e4.to_one()
        #insert_scaled_sum() will raise exception for 1D image
        Log.logger().set_level(-1)    #perfect solution for quenching the Log error information, thank Liwei
        self.assertRaises( RuntimeError, e3.insert_scaled_sum, e4, (0,0,0))
        try:
            e3.insert_scaled_sum(e4, (0,0,0))
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "ImageDimensionException")
            
    def test_window_center(self):
        """test window_center() function ...................."""
        e = EMData()
        e.set_size(64,64,64)
        e.process_inplace("testimage.noise.uniform.rand")
        e2 = e.window_center(32)
        
        #window_padded() only work for real data, raise exception for complex
        e.set_complex(True)
        self.assertRaises( RuntimeError, e.window_center, 32 )
        try:
            e2 = e.window_center(32)
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "ImageFormatException")
            
        #window_padded() apply to cubic real space image only
        e3 = EMData()
        e3.set_size(64,64,62)
        e3.process_inplace("testimage.noise.uniform.rand")
        self.assertRaises( RuntimeError, e3.window_center, 32 )
        try:
            e4 = e3.window_center(32)
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "ImageFormatException")
        
        #test 2-D support
        e7 = EMData()
        e7.set_size(64,64,1)
        e7.process_inplace("testimage.noise.uniform.rand")
        e8 = e7.window_center(32)
            
        #test 1-D support 
        e5 = EMData()
        e5.set_size(64,1,1)
        e5.process_inplace("testimage.noise.uniform.rand")
        e6 = e5.window_center(32)   
            
    def test_center_origin(self):
        """test center_origin() function ...................."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace("testimage.noise.uniform.rand")
        e.center_origin()
        
        #center_origin() only apply to real image
        e.set_complex(True)
        Log.logger().set_level(-1)
        self.assertRaises( RuntimeError, e.center_origin, )
        try:
            e.center_origin()
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "ImageFormatException")
    
    def test_center_origin_fft(self):
        """test center_origin_fft() function ................"""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace("testimage.noise.uniform.rand")
        e.set_complex(True)
        e.center_origin_fft()
        
        #center_origin_fft() apply to complex image only
        e.set_complex(False)
        self.assertRaises( RuntimeError, e.center_origin_fft, )
        try:
            e.center_origin_fft()
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "ImageFormatException")
        
    def test_zeropad_ntimes(self):
        """test zeropad_ntimes() function ..................."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace("testimage.noise.uniform.rand")
        
        #test default argument npad=4
        e2 = e.zeropad_ntimes()
        d = e.get_3dview()
        d2 = e2.get_3dview()
        for z in range(32):
            for y in range(32):
                for x in range(32):
                    self.assertEqual( d[z][y][x], d2[z+48][y+48][x+48])
        
        #test argument npad=3
        e3 = e.zeropad_ntimes(3)
        d3 = e3.get_3dview()
        for z in range(32):
            for y in range(32):
                for x in range(32):
                    self.assertEqual( d[z][y][x], d3[z+32][y+32][x+32])
                    
    def test_pad_fft(self):
        """test pad_fft() function .........................."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace("testimage.noise.uniform.rand")
        
        #test default argument
        e2 = e.pad_fft()
        
        #test arbitrary argument
        e3 = e.pad_fft(3)
        
    def test_postift_depad_corner_inplace(self):
        """test postift_depad_corner_inplace() .............."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace("testimage.noise.uniform.rand")
        e2 = e.do_fft()
        e3 = e2.do_ift()
        e3.postift_depad_corner_inplace()        
        #test the correctness for this function, something I don't quite understand
        #d = e.get_3dview()
        #d3 = e3.get_3dview()
        #import math
        #for z in range(31):
        #    for y in range(31):
        #        for x in range(31):
        #            self.assertEqual( math.ceil((d[z][y][x]-d3[z][y][x])*1000), 0 )

    def test_real2FH(self):
        """test real2FH() function .........................."""
        e = EMData()
        e.set_size(31,31,1)
        e.process_inplace("testimage.noise.uniform.rand")
        #Log.logger().set_level(-1)
        #import sys
        #outfile = "out.txt" 
        #sys.stdout = open(outfile,"w")
        e3 = e.real2FH(1.0)
        
        #real2FH apply to 2D/Square/Real/odd image
        e2 = EMData()
        e2.set_size(31,31,31)
        self.assertRaises( RuntimeError, e2.real2FH, 1.0)
        try:
            e3 = e2.real2FH(1.0)
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "ImageFormatException")
        
        e2.set_size(31,21,1)
        self.assertRaises( RuntimeError, e2.real2FH, 1.0)
        try:
            e3 = e2.real2FH(1.0)
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "ImageFormatException")
            
        e2.set_size(32,32,1)
        self.assertRaises( RuntimeError, e2.real2FH, 1.0)
        try:
            e3 = e2.real2FH(1.0)
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "ImageFormatException")
        
        e2.set_size(31,31,1)
        e2.set_complex(True)
        self.assertRaises( RuntimeError, e2.real2FH, 1.0)
        try:
            e2.real2FH(1.0)
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "ImageFormatException")
        
    def test_FH2F(self):
        """test FH2F() function ............................."""
        e = EMData()
        e.set_size(31,31,1)
        e.process_inplace("testimage.noise.uniform.rand")
        e2 = e.real2FH(1.0)
        e3 = e2.FH2F(31, 1.0)
        
        #for image not FH, should raise exception
        Log.logger().set_level(-1)
        self.assertRaises( RuntimeError, e.FH2F, 31, 1.0)
        try:
            e.FH2F(31, 1.0)
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "ImageFormatException")

    def test_do_fft(self):
        """test do_fft()/do_ift() function .................."""
        #test 3D image
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace("testimage.noise.uniform.rand")
        e2 = e.do_fft()
        e3 = e2.do_ift()
        
        for i in range(32):
            for j in range(32):
                for k in range(32):
                    self.assertAlmostEqual(e.get_value_at(i,j,k), e3.get_value_at(i,j,k), 3)
        
        #test 1D image
        e4 = EMData()
        e4.set_size(32,1,1)
        e.process_inplace("testimage.noise.uniform.rand")
        e5 = e4.do_fft()
        e6 = e5.do_ift()
        
        for l in range(32):
                self.assertAlmostEqual(e4.get_value_at(l, 0, 0), e6.get_value_at(l, 0, 0), 3)
        
        #test 2D image
        e7 = EMData()
        e7.set_size(32,32,1)
        e7.process_inplace("testimage.noise.uniform.rand")
        e8 = e7.do_fft()
        e9 = e8.do_ift()
        
        for m in range(32):
            for n in range(32):
                self.assertAlmostEqual(e7.get_value_at(m,n,0), e9.get_value_at(m,n,0), 3)
        
        #do_fft() only apply to real image
        self.assertRaises( RuntimeError, e2.do_fft, )
        try:
            e2.do_fft()
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "ImageFormatException")
        
        #do_ift() only apply to complex image
        self.assertRaises( RuntimeError, e3.do_ift, )
        try:
            e3.do_ift()
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "ImageFormatException")
    
    def test_do_fft_complex_value(self):
        """test the complex image values after FFT .........."""
        #test 1D image(5)
        e = EMData()
        e.set_size(5,1,1)
        for i in range(5):
            e.set_value_at(i,0,0, i+1)
        f = e.do_fft()
        self.assertAlmostEqual(f.get_value_at(0,0,0), 15, 2)
        self.assertAlmostEqual(f.get_value_at(1,0,0), 0, 2)
        self.assertAlmostEqual(f.get_value_at(2,0,0), -2.5, 2)
        self.assertAlmostEqual(f.get_value_at(3,0,0), 3.44, 2)
        self.assertAlmostEqual(f.get_value_at(4,0,0), -2.5, 2)
        self.assertAlmostEqual(f.get_value_at(5,0,0), 0.81, 2)
        
        #test 1D image(4)
        e2 = EMData()
        e2.set_size(4,1,1)
        for i in range(4):
            e2.set_value_at(i,0,0, i+1)
        f2 = e2.do_fft()
        self.assertAlmostEqual(f2.get_value_at(0,0,0), 10, 2)
        self.assertAlmostEqual(f2.get_value_at(1,0,0), 0, 2)
        self.assertAlmostEqual(f2.get_value_at(2,0,0), -2, 2)
        self.assertAlmostEqual(f2.get_value_at(3,0,0), 2, 2)
        self.assertAlmostEqual(f2.get_value_at(4,0,0), -2, 2)
        self.assertAlmostEqual(f2.get_value_at(5,0,0), 0, 2)
        
        #test 2D image(2x2)
        e3 = EMData()
        e3.set_size(2,2,1)
        value = 0
        for i in range(2):
            for j in range(2):
                value += 1
                e3.set_value_at(j, i, 0, value)
        #print e3.get_2dview() 
        f3 = e3.do_fft()
        self.assertAlmostEqual(f3.get_value_at(0,0,0), 10, 2)
        self.assertAlmostEqual(f3.get_value_at(1,0,0), 0, 2)
        self.assertAlmostEqual(f3.get_value_at(2,0,0), -2, 2)
        self.assertAlmostEqual(f3.get_value_at(3,0,0), 0, 2)
        self.assertAlmostEqual(f3.get_value_at(0,1,0), -4, 2)
        self.assertAlmostEqual(f3.get_value_at(1,1,0), 0, 2)
        self.assertAlmostEqual(f3.get_value_at(2,1,0), 0, 2)
        self.assertAlmostEqual(f3.get_value_at(3,1,0), 0, 2)
        
        #test 2D image(3x2)
        e4 = EMData()
        e4.set_size(3,2,1)
        value = 0
        for i in range(2):
            for j in range(3):
                value += 1
                e4.set_value_at(j, i, 0, value)
        #print '\n', e4.get_2dview() 
        f4 = e4.do_fft()
        self.assertAlmostEqual(f4.get_value_at(0,0,0), 21, 2)
        self.assertAlmostEqual(f4.get_value_at(1,0,0), 0, 2)
        self.assertAlmostEqual(f4.get_value_at(2,0,0), -3, 2)
        self.assertAlmostEqual(f4.get_value_at(3,0,0), 1.732, 2)
        self.assertAlmostEqual(f4.get_value_at(0,1,0), -9, 2)
        self.assertAlmostEqual(f4.get_value_at(1,1,0), 0, 2)
        self.assertAlmostEqual(f4.get_value_at(2,1,0), 0, 2)
        self.assertAlmostEqual(f4.get_value_at(3,1,0), 0, 2)
        
        #test 3D image(2x2x2)
        e5 = EMData()
        e5.set_size(2,2,2)
        value = 0
        for i in range(2):
            for j in range(2):
                for k in range(2):
                    value += 1
                    e5.set_value_at(k, j, i, value)
        #print '\n', e5.get_3dview()
        f5 = e5.do_fft()
        self.assertAlmostEqual(f5.get_value_at(0,0,0), 36, 2)
        self.assertAlmostEqual(f5.get_value_at(1,0,0), 0, 2)
        self.assertAlmostEqual(f5.get_value_at(2,0,0), -4, 2)
        self.assertAlmostEqual(f5.get_value_at(3,0,0), 0, 2)
        self.assertAlmostEqual(f5.get_value_at(0,1,0), -8, 2)
        self.assertAlmostEqual(f5.get_value_at(1,1,0), 0, 2)
        self.assertAlmostEqual(f5.get_value_at(2,1,0), 0, 2)
        self.assertAlmostEqual(f5.get_value_at(3,1,0), 0, 2)
        self.assertAlmostEqual(f5.get_value_at(0,0,1), -16, 2)
        self.assertAlmostEqual(f5.get_value_at(1,0,1), 0, 2)
        self.assertAlmostEqual(f5.get_value_at(2,0,1), 0, 2)
        self.assertAlmostEqual(f5.get_value_at(3,0,1), 0, 2)
        self.assertAlmostEqual(f5.get_value_at(0,1,1), 0, 2)
        self.assertAlmostEqual(f5.get_value_at(1,1,1), 0, 2)
        self.assertAlmostEqual(f5.get_value_at(2,1,1), 0, 2)
        self.assertAlmostEqual(f5.get_value_at(3,1,1), 0, 2)
        
        #test 3D image(3x2x2)
        e6 = EMData()
        e6.set_size(3,2,2)
        value = 0
        for i in range(2):
            for j in range(2):
                for k in range(3):
                    value += 1
                    e6.set_value_at(k, j, i, value)
        #print '\n', e6.get_3dview()
        f6 = e6.do_fft()
        #print '\n', f6.get_3dview()
        self.assertAlmostEqual(f6.get_value_at(0,0,0), 78, 2)
        self.assertAlmostEqual(f6.get_value_at(1,0,0), 0, 2)
        self.assertAlmostEqual(f6.get_value_at(2,0,0), -6, 2)
        self.assertAlmostEqual(f6.get_value_at(3,0,0), 3.464, 2)
        self.assertAlmostEqual(f6.get_value_at(0,1,0), -18, 2)
        self.assertAlmostEqual(f6.get_value_at(1,1,0), 0, 2)
        self.assertAlmostEqual(f6.get_value_at(2,1,0), 0, 2)
        self.assertAlmostEqual(f6.get_value_at(3,1,0), 0, 2)
        self.assertAlmostEqual(f6.get_value_at(0,0,1), -36, 2)
        self.assertAlmostEqual(f6.get_value_at(1,0,1), 0, 2)
        self.assertAlmostEqual(f6.get_value_at(2,0,1), 0, 2)
        self.assertAlmostEqual(f6.get_value_at(3,0,1), 0, 2)
        self.assertAlmostEqual(f6.get_value_at(0,1,1), 0, 2)
        self.assertAlmostEqual(f6.get_value_at(1,1,1), 0, 2)
        self.assertAlmostEqual(f6.get_value_at(2,1,1), 0, 2)
        self.assertAlmostEqual(f6.get_value_at(3,1,1), 0, 2)
    
    def test_do_fft_inplace_value(self):
        """test inplace fft/ift numerical correctness ......."""
        #test 1D image
        e = EMData()
        e.set_size(4,1,1)
        for i in range(4):
            e.set_value_at(i,0,0, i+1)
        e.do_fft_inplace()
        self.assertEqual(e.get_xsize(), 6)
        self.assertEqual(e.is_complex(), True)
        self.assertEqual(e.is_complex_x(), True)
        self.assertAlmostEqual(e.get_value_at(0,0,0), 10, 2)
        self.assertAlmostEqual(e.get_value_at(1,0,0), 0, 2)
        self.assertAlmostEqual(e.get_value_at(2,0,0), -2, 2)
        self.assertAlmostEqual(e.get_value_at(3,0,0), 2, 2)
        self.assertAlmostEqual(e.get_value_at(4,0,0), -2, 2)
        self.assertAlmostEqual(e.get_value_at(5,0,0), 0, 2)
        
        e.do_ift_inplace()
        self.assertEqual(e.get_xsize(), 4)
        self.assertEqual(e.is_complex(), False)
        self.assertAlmostEqual(e.get_value_at(0,0,0), 1, 2)
        self.assertAlmostEqual(e.get_value_at(1,0,0), 2, 2)
        self.assertAlmostEqual(e.get_value_at(2,0,0), 3, 2)
        self.assertAlmostEqual(e.get_value_at(3,0,0), 4, 2)
        
        e2 = EMData()
        e2.set_size(5,1,1)
        for i in range(5):
            e2.set_value_at(i,0,0, i+1)
        e2.do_fft_inplace()
        self.assertEqual(e2.get_xsize(), 6)
        self.assertEqual(e2.is_complex(), True)
        self.assertAlmostEqual(e2.get_value_at(0,0,0), 15, 2)
        self.assertAlmostEqual(e2.get_value_at(1,0,0), 0, 2)
        self.assertAlmostEqual(e2.get_value_at(2,0,0), -2.5, 2)
        self.assertAlmostEqual(e2.get_value_at(3,0,0), 3.44, 2)
        self.assertAlmostEqual(e2.get_value_at(4,0,0), -2.5, 2)
        self.assertAlmostEqual(e2.get_value_at(5,0,0), 0.81, 2)
        
        e2.do_ift_inplace()
        self.assertEqual(e2.get_xsize(), 5)
        self.assertEqual(e2.is_complex(), False)
        self.assertAlmostEqual(e2.get_value_at(0,0,0), 1, 2)
        self.assertAlmostEqual(e2.get_value_at(1,0,0), 2, 2)
        self.assertAlmostEqual(e2.get_value_at(2,0,0), 3, 2)
        self.assertAlmostEqual(e2.get_value_at(3,0,0), 4, 2)
        self.assertAlmostEqual(e2.get_value_at(4,0,0), 5, 2)
        
        #test 2D image(2x2)
        e3 = EMData()
        e3.set_size(2,2,1)
        value = 0
        for i in range(2):
            for j in range(2):
                value += 1
                e3.set_value_at(j, i, 0, value)
        #print e3.get_2dview() 
        e3.do_fft_inplace()
        self.assertEqual(e3.get_xsize(), 4)
        self.assertEqual(e3.get_ysize(), 2)
        self.assertEqual(e3.is_complex(), True)
        self.assertAlmostEqual(e3.get_value_at(0,0,0), 10, 2)
        self.assertAlmostEqual(e3.get_value_at(1,0,0), 0, 2)
        self.assertAlmostEqual(e3.get_value_at(2,0,0), -2, 2)
        self.assertAlmostEqual(e3.get_value_at(3,0,0), 0, 2)
        self.assertAlmostEqual(e3.get_value_at(0,1,0), -4, 2)
        self.assertAlmostEqual(e3.get_value_at(1,1,0), 0, 2)
        self.assertAlmostEqual(e3.get_value_at(2,1,0), 0, 2)
        self.assertAlmostEqual(e3.get_value_at(3,1,0), 0, 2)
        
        e3.do_ift_inplace()
        #print '\n', e3.get_2dview() 
        self.assertEqual(e3.get_xsize(), 2)
        self.assertEqual(e3.get_ysize(), 2)
        self.assertEqual(e3.is_complex(), False)
        self.assertAlmostEqual(e3.get_value_at(0,0,0), 1, 2)
        self.assertAlmostEqual(e3.get_value_at(1,0,0), 2, 2)
        self.assertAlmostEqual(e3.get_value_at(0,1,0), 3, 2)
        self.assertAlmostEqual(e3.get_value_at(1,1,0), 4, 2)
        
        #test 2D image(3x2)
        e4 = EMData()
        e4.set_size(3,2,1)
        value = 0
        for i in range(2):
            for j in range(3):
                value += 1
                e4.set_value_at(j, i, 0, value)
        #print '\n', e4.get_2dview() 
        e4.do_fft_inplace()
        self.assertEqual(e4.get_xsize(), 4)
        self.assertEqual(e4.get_ysize(), 2)
        self.assertEqual(e4.is_complex(), True)
        self.assertAlmostEqual(e4.get_value_at(0,0,0), 21, 2)
        self.assertAlmostEqual(e4.get_value_at(1,0,0), 0, 2)
        self.assertAlmostEqual(e4.get_value_at(2,0,0), -3, 2)
        self.assertAlmostEqual(e4.get_value_at(3,0,0), 1.732, 2)
        self.assertAlmostEqual(e4.get_value_at(0,1,0), -9, 2)
        self.assertAlmostEqual(e4.get_value_at(1,1,0), 0, 2)
        self.assertAlmostEqual(e4.get_value_at(2,1,0), 0, 2)
        self.assertAlmostEqual(e4.get_value_at(3,1,0), 0, 2)
        
        e4.do_ift_inplace()
        self.assertEqual(e4.get_xsize(), 3)
        self.assertEqual(e4.get_ysize(), 2)
        self.assertEqual(e4.is_complex(), False)
        self.assertAlmostEqual(e4.get_value_at(0,0,0), 1, 2)
        self.assertAlmostEqual(e4.get_value_at(1,0,0), 2, 2)
        self.assertAlmostEqual(e4.get_value_at(2,0,0), 3, 2)
        self.assertAlmostEqual(e4.get_value_at(0,1,0), 4, 2)
        self.assertAlmostEqual(e4.get_value_at(1,1,0), 5, 2)
        self.assertAlmostEqual(e4.get_value_at(2,1,0), 6, 2)
        
        #test 3D image(2x2x2)
        e5 = EMData()
        e5.set_size(2,2,2)
        value = 0
        for i in range(2):
            for j in range(2):
                for k in range(2):
                    value += 1
                    e5.set_value_at(k, j, i, value)
        #print '\n', e5.get_3dview()
        e5.do_fft_inplace()
        self.assertEqual(e5.get_xsize(), 4)
        self.assertEqual(e5.get_ysize(), 2)
        self.assertEqual(e5.get_zsize(), 2)
        self.assertEqual(e5.is_complex(), True)
        self.assertAlmostEqual(e5.get_value_at(0,0,0), 36, 2)
        self.assertAlmostEqual(e5.get_value_at(1,0,0), 0, 2)
        self.assertAlmostEqual(e5.get_value_at(2,0,0), -4, 2)
        self.assertAlmostEqual(e5.get_value_at(3,0,0), 0, 2)
        self.assertAlmostEqual(e5.get_value_at(0,1,0), -8, 2)
        self.assertAlmostEqual(e5.get_value_at(1,1,0), 0, 2)
        self.assertAlmostEqual(e5.get_value_at(2,1,0), 0, 2)
        self.assertAlmostEqual(e5.get_value_at(3,1,0), 0, 2)
        self.assertAlmostEqual(e5.get_value_at(0,0,1), -16, 2)
        self.assertAlmostEqual(e5.get_value_at(1,0,1), 0, 2)
        self.assertAlmostEqual(e5.get_value_at(2,0,1), 0, 2)
        self.assertAlmostEqual(e5.get_value_at(3,0,1), 0, 2)
        self.assertAlmostEqual(e5.get_value_at(0,1,1), 0, 2)
        self.assertAlmostEqual(e5.get_value_at(1,1,1), 0, 2)
        self.assertAlmostEqual(e5.get_value_at(2,1,1), 0, 2)
        self.assertAlmostEqual(e5.get_value_at(3,1,1), 0, 2)
        
        e5.do_ift_inplace()
        self.assertEqual(e5.get_xsize(), 2)
        self.assertEqual(e5.get_ysize(), 2)
        self.assertEqual(e5.get_zsize(), 2)
        self.assertEqual(e5.is_complex(), False)
        self.assertAlmostEqual(e5.get_value_at(0,0,0), 1, 2)
        self.assertAlmostEqual(e5.get_value_at(1,0,0), 2, 2)
        self.assertAlmostEqual(e5.get_value_at(0,1,0), 3, 2)
        self.assertAlmostEqual(e5.get_value_at(1,1,0), 4, 2)
        self.assertAlmostEqual(e5.get_value_at(0,0,1), 5, 2)
        self.assertAlmostEqual(e5.get_value_at(1,0,1), 6, 2)
        self.assertAlmostEqual(e5.get_value_at(0,1,1), 7, 2)
        self.assertAlmostEqual(e5.get_value_at(1,1,1), 8, 2)
        
        #test 3D image(3x2x2)
        e6 = EMData()
        e6.set_size(3,2,2)
        value = 0
        for i in range(2):
            for j in range(2):
                for k in range(3):
                    value += 1
                    e6.set_value_at(k, j, i, value)
        e6.do_fft_inplace()
        self.assertEqual(e6.is_complex(), True)
        self.assertEqual(e6.get_xsize(), 4)
        self.assertEqual(e6.get_ysize(), 2)
        self.assertEqual(e6.get_zsize(), 2)
        self.assertAlmostEqual(e6.get_value_at(0,0,0), 78, 2)
        self.assertAlmostEqual(e6.get_value_at(1,0,0), 0, 2)
        self.assertAlmostEqual(e6.get_value_at(2,0,0), -6, 2)
        self.assertAlmostEqual(e6.get_value_at(3,0,0), 3.464, 2)
        self.assertAlmostEqual(e6.get_value_at(0,1,0), -18, 2)
        self.assertAlmostEqual(e6.get_value_at(1,1,0), 0, 2)
        self.assertAlmostEqual(e6.get_value_at(2,1,0), 0, 2)
        self.assertAlmostEqual(e6.get_value_at(3,1,0), 0, 2)
        self.assertAlmostEqual(e6.get_value_at(0,0,1), -36, 2)
        self.assertAlmostEqual(e6.get_value_at(1,0,1), 0, 2)
        self.assertAlmostEqual(e6.get_value_at(2,0,1), 0, 2)
        self.assertAlmostEqual(e6.get_value_at(3,0,1), 0, 2)
        self.assertAlmostEqual(e6.get_value_at(0,1,1), 0, 2)
        self.assertAlmostEqual(e6.get_value_at(1,1,1), 0, 2)
        self.assertAlmostEqual(e6.get_value_at(2,1,1), 0, 2)
        self.assertAlmostEqual(e6.get_value_at(3,1,1), 0, 2)
        
        e6.do_ift_inplace()
        self.assertEqual(e6.is_complex(), False)
        self.assertEqual(e6.get_xsize(), 3)
        self.assertEqual(e6.get_ysize(), 2)
        self.assertEqual(e6.get_zsize(), 2)
        self.assertAlmostEqual(e6.get_value_at(0,0,0), 1, 2)
        self.assertAlmostEqual(e6.get_value_at(1,0,0), 2, 2)
        self.assertAlmostEqual(e6.get_value_at(2,0,0), 3, 2)
        self.assertAlmostEqual(e6.get_value_at(0,1,0), 4, 2)
        self.assertAlmostEqual(e6.get_value_at(1,1,0), 5, 2)
        self.assertAlmostEqual(e6.get_value_at(2,1,0), 6, 2)
        self.assertAlmostEqual(e6.get_value_at(0,0,1), 7, 2)
        self.assertAlmostEqual(e6.get_value_at(1,0,1), 8, 2)
        self.assertAlmostEqual(e6.get_value_at(2,0,1), 9, 2)
        self.assertAlmostEqual(e6.get_value_at(0,1,1), 10, 2)
        self.assertAlmostEqual(e6.get_value_at(1,1,1), 11, 2)
        self.assertAlmostEqual(e6.get_value_at(2,1,1), 12, 2)
        
    def test_do_fft_inplace(self):
        """test do_fft_inplace()/do_ift_inplace ............."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace("testimage.noise.uniform.rand")
        
        #test unpadded real image
        e2 = e.do_fft_inplace()    
        
        #test padded real image
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace("testimage.noise.uniform.rand")
        e2 = e.pad_fft()
        e3 = e2.do_fft_inplace()
        
        e.set_complex(True)
        e.do_ift_inplace()
        
        #do_ift_inplace() only apply to complex image
        e4 = EMData()
        e4.set_size(32,32,32)
        e4.process_inplace("testimage.noise.uniform.rand")
        self.assertRaises( RuntimeError, e4.do_ift_inplace, )
        try:
            e4.do_ift_inplace()
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "ImageFormatException")
        
        #do_fft_inplace() only apply to real image
        e5 = EMData()
        e5.set_size(32,32,32)
        e5.process_inplace("testimage.noise.uniform.rand")
        e6 = e5.do_fft()
        self.assertRaises( RuntimeError, e6.do_fft_inplace, )
        try:
            e6.do_fft_inplace()
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "ImageFormatException")

    def test_get_fft_amplitude(self):
        """test get_fft_amplitude() function ................"""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace("testimage.noise.uniform.rand")
        e2 = e.do_fft()
        e3 = e2.get_fft_amplitude()
        self.assertEqual( e3.is_complex(), False)
        
        #this function only apply to complex image
        self.assertRaises( RuntimeError, e.get_fft_amplitude, )
        try:
            e.get_fft_amplitude()
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "ImageFormatException")

    def test_get_fft_phase(self):
        """test get_fft_phase() function ...................."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace("testimage.noise.uniform.rand")
        e2 = e.do_fft()
        e3 = e2.get_fft_phase()
        self.assertEqual( e3.is_complex(), False)
        
        #this function only apply to complex image
        self.assertRaises( RuntimeError, e.get_fft_phase, )
        try:
            e.get_fft_phase()
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "ImageFormatException")
            
    def test_get_fft_amplitude2D(self):
        """test get_fft_amplitude2D() function .............."""
        e = EMData()
        e.set_size(32,32,1)
        e.process_inplace("testimage.noise.uniform.rand")
        e2 = e.do_fft()
        e3 = e2.get_fft_amplitude2D()
        self.assertEqual( e3.is_complex(), False)
        
        #this function only apply to 2D image
        e4 = EMData()
        e4.set_size(32,32,32)
        self.assertRaises( RuntimeError, e4.get_fft_amplitude2D, )
        try:
            e4.get_fft_amplitude2D()
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "ImageFormatException")
        
        #this function only apply to complex image
        self.assertRaises( RuntimeError, e.get_fft_amplitude2D, )
        try:
            e.get_fft_amplitude2D()
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "ImageFormatException")
    
    def test_render_amp8(self):
        """test render_amp8() function ......................"""
        e = EMData()
        e.set_size(32,32,1)
        e.process_inplace("testimage.noise.uniform.rand")
        str = e.render_amp8(0, 0, 32, 32, 96, 1.2, 1, 254, 100, 200, 3)
        
        #only apply to 2D image
        e2 = EMData()
        e2.set_size(32,32,32)
        e2.process_inplace("testimage.noise.uniform.rand")
        self.assertRaises( RuntimeError, e2.render_amp8, 0, 0, 32, 32, 96, 1.2, 1, 254, 100, 200, 3)
        try:
            str = e2.render_amp8(0, 0, 32, 32, 96, 1.2, 1, 254, 100, 200, 3)
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "ImageDimensionException")
            
    def test_ri2ap_ap2ri(self):
        """test ri2ap()/ap2ri() function ...................."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace("testimage.noise.uniform.rand")
        e2 = e.do_fft()
        e2.ri2ap()
        e2.ap2ri()
         
    def test_scale(self):
        """test scale() function ............................"""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace("testimage.noise.uniform.rand")
        d = e.get_3dview()
        
        e2 = e.scale(2.0)
        self.assertEqual(e2, None) #this function return None(void in C++)    

    def test_make_rotational_footprint(self):
        """test make_rotational_footprint() function ........"""
        e = EMData()
        e.set_size(64,64,64)
        e.to_one()
        e.make_rotational_footprint()
        
        #test for bad input, only even sized image accepted, ImageFormatException raised 
        e2 = EMData()
        e.set_size(31,31,31)
        e.to_one()
        try:
            e.make_rotational_footprint()
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "ImageFormatException")

    def test_to_zero(self):
        """test to_zero() function .........................."""
        from random import randint
        nx = randint(1,100)
        ny = randint(1,100)
        nz = randint(1,100)
        e = EMData()
        e.set_size(nx,ny,nz)
        e.to_zero()
        for k in range(nz):
            for j in range(ny):
                for i in range(nx):
                    self.assertEqual(e.get_value_at(i, j, k), 0.0)
    
    def test_to_one(self):
        """test to_one() function ..........................."""
        from random import randint
        nx = randint(1,100)
        ny = randint(1,100)
        nz = randint(1,100)
        e = EMData()
        e.set_size(nx,ny,nz)
        e.to_one()
        for k in range(nz):
            for j in range(ny):
                for i in range(nx):
                    self.assertEqual(e.get_value_at(i, j, k), 1.0)
           
    def test_real_operator_unary(self):
        """test real unary operator of EMData ..............."""
        e = EMData()
        e.set_size(32,32,32)
        e.to_one()
        
        e += 0.5
        for i in range(32):
            for j in range(32):
                for k in range(32):
                    self.assertEqual(e.get_value_at(i, j, k), 1.5)
        
        e -= 0.5
        for i in range(32):
            for j in range(32):
                for k in range(32):
                    self.assertEqual(e.get_value_at(i, j, k), 1.0)
        
        e *= 2.0
        for i in range(32):
            for j in range(32):
                for k in range(32):
                    self.assertEqual(e.get_value_at(i, j, k), 2.0)
        
        e /= 2.0
        for i in range(32):
            for j in range(32):
                for k in range(32):
                    self.assertEqual(e.get_value_at(i, j, k), 1.0)
 
        e2 = EMData()
        e2.set_size(32,32,32)
        e2.to_one()

        e += e2
        for i in range(32):
            for j in range(32):
                for k in range(32):
                    self.assertEqual(e.get_value_at(i, j, k), 2.0)
        
        e -= e2
        for i in range(32):
            for j in range(32):
                for k in range(32):
                    self.assertEqual(e.get_value_at(i, j, k), 1.0)
        
        e *= 2.0
        e2 *= 3.0
        e *= e2
        for i in range(32):
            for j in range(32):
                for k in range(32):
                    self.assertEqual(e.get_value_at(i, j, k), 6.0)
        
        e /= e2
        for i in range(32):
            for j in range(32):
                for k in range(32):
                    self.assertEqual(e.get_value_at(i, j, k), 2.0)
  
    def test_multi_array_2d(self):
        """test multi_array_2d real.........................."""
        nx = 16
        ny = 32

        e = EMData()
        e.set_size(nx, ny)
        array = e.get_2dview()
        self.assertEqual(type(array).__name__, 'array')
        self.assertEqual(array.typecode(), "f")
        self.assertEqual(array.shape, (ny, nx))

        for i in range(ny):
            for j in range(nx):
                self.assertEqual(e.get_value_at(j,i), array[i][j])

        e *= 2
        for i in range(ny):
            for j in range(nx):
                self.assertEqual(e.get_value_at(j,i), array[i][j])

        array.savespace(1)
        array *= 2
        
        for i in range(ny):
            for j in range(nx):
                self.assertEqual(e.get_value_at(j,i), array[i][j])

    def test_multi_array_3d(self):
        """test multi_array_3d real.........................."""
        nx = 8
        ny = 16
        nz = 4

        e = EMData()
        e.set_size(nx, ny, nz)
        array = e.get_3dview()
        self.assertEqual(type(array).__name__, 'array')
        self.assertEqual(array.typecode(), "f")
        self.assertEqual(array.shape, (nz, ny, nx))

        for i in range(nz):
            for j in range(ny):
                for k in range(nx):
                    self.assertEqual(e.get_value_at(k,j,i), array[i][j][k])

        e *= 2
        for i in range(nz):
            for j in range(ny):
                for k in range(nx):
                    self.assertEqual(e.get_value_at(k,j,i), array[i][j][k])

        array.savespace(1)
        array *= 2
        for i in range(nz):
            for j in range(ny):
                for k in range(nx):
                    self.assertEqual(e.get_value_at(k,j,i), array[i][j][k])
      
    def test_multi_array_c2d(self):
        """test multi_array 2d complex ......................"""
        nx = 16
        ny = 16
        infile = "test_multi_array_c2d.mrc"
        TestUtil.make_image_file(infile, MRC, EM_FLOAT, nx, ny)

        e = EMData()
        e.read_image(infile)
        fft = e.do_fft()
        nx2 = fft.get_xsize()
        array = fft.get_2dcview()
        
        self.assertEqual(type(array).__name__, 'array')
        self.assertEqual(array.typecode(), "F")
        self.assertEqual(array.shape, (ny, nx2/2))

        for i in range(ny):
            for j in range(0, nx2, 2):
                c1 = array[i][j/2]
                testlib.assertfloat(self, fft.get_value_at(j,i), c1.real)
                testlib.assertfloat(self, fft.get_value_at(j+1,i), c1.imag)


        e *= 2
        for i in range(ny):
            for j in range(0, nx2, 2):
                c1 = array[i][j/2]
                testlib.assertfloat(self, fft.get_value_at(j,i), c1.real)
                testlib.assertfloat(self, fft.get_value_at(j+1,i), c1.imag)

        array.savespace(1)
        array *= 2
        for i in range(ny):
            for j in range(0, nx2, 2):
                c1 = array[i][j/2]
                testlib.assertfloat(self, fft.get_value_at(j,i), c1.real)
                testlib.assertfloat(self, fft.get_value_at(j+1,i), c1.imag)

        os.unlink(infile)

    def test_multi_array_c3d(self):
        """test multi_array 3d complex ......................"""
        nx = 8
        ny = 16
        nz = 4
        infile = "test_multi_array_c3d.mrc"
        TestUtil.make_image_file(infile, MRC, EM_FLOAT, nx, ny, nz)

        e = EMData()
        e.read_image(infile)
        fft = e.do_fft()

        nx2 = fft.get_xsize()
        
        array = fft.get_3dcview()
        self.assertEqual(type(array).__name__, 'array')
        self.assertEqual(array.typecode(), "F")
        self.assertEqual(array.shape, (nz, ny, nx2/2))

        for i in range(nz):
            for j in range(ny):
                for k in range(0,nx2,2):
                    c1 = array[i][j][k/2]
                    testlib.assertfloat(self, fft.get_value_at(k,j,i), c1.real)
                    testlib.assertfloat(self, fft.get_value_at(k+1,j,i), c1.imag)

        e *= 2
        for i in range(nz):
            for j in range(ny):
                for k in range(0,nx2,2):
                    c1 = array[i][j][k/2]
                    testlib.assertfloat(self, fft.get_value_at(k,j,i), c1.real)
                    testlib.assertfloat(self, fft.get_value_at(k+1,j,i), c1.imag)


        array.savespace(1)
        array *= 2
        for i in range(nz):
            for j in range(ny):
                for k in range(0,nx2,2):
                    c1 = array[i][j][k/2]
                    testlib.assertfloat(self, fft.get_value_at(k,j,i), c1.real)
                    testlib.assertfloat(self, fft.get_value_at(k+1,j,i), c1.imag)

        os.unlink(infile)
        
    def test_translate(self):
        """test translate() function ........................"""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace("testimage.noise.uniform.rand")
        
        #integer translation, no interpolation involved
        e.translate(2,2,2)
        
        #float translation, interpolation involved
        e.translate(1.2, 2.1, 3.3)
        
        #Vec3i is tuple 3 of int in Python
        e.translate((1,2,3))
        
        #Vec3f is tuple 3 of float in Python
        #e.translate((1.1,2.1,3.1))     #problem here, Vec3f not reconized from flaot tuple 3

    def test_rotate_translate(self):
        """test rotate_translate ............................"""
        infile = "test_rotate_translate.mrc"
        TestUtil.make_image_file(infile, MRC, EM_FLOAT, 16,16,16)
        
        x=EMData()
        x.read_image(infile)

        alt = 1.0329837512591338
        az = 3.7260642381912579
        phi = 5.7671541529246966

        x.rotate_translate(az, alt, phi, 2, 2, 2)
        testlib.check_emdata(x, sys.argv[0])

        x.rotate_translate(az, alt, phi, 6, 6, 6)
        testlib.check_emdata(x, sys.argv[0])

        os.unlink(infile)
                
    def test_rotate_2d(self):
        """test rotate_2d ..................................."""
        infile = "test_rotate_2d.mrc"
        TestUtil.make_image_file(infile, MRC, EM_FLOAT, 24, 32)

        outfile1 = "test_rotate_2d_out_1.mrc"
        outfile2 = "test_rotate_2d_out_2.mrc"
        
        a = EMData()
        a.read_image(infile)
        b=a.copy()
        b.rotate(0, 0, math.pi/4)
        b.write_image(outfile1)
        
        # verify b
        
        b=a.copy()
        b.rotate(0, 0, math.pi/2)
        b.write_image(outfile2)

        # verify b

        os.unlink(infile)
        os.unlink(outfile1)
        os.unlink(outfile2)
        
    def test_rotate_3d(self):
        """test rotate_3d ..................................."""
        a = EMData()
        a.set_size(72,72,72)
        a.to_one()
        
        b=a.copy()
        b.rotate(0,0,math.pi/4)
        testlib.check_emdata(b, sys.argv[0])
    
        b=a.copy()
        b.rotate(0,0,math.pi/2)
        testlib.check_emdata(b, sys.argv[0])

    def test_rotate_x(self):
        """test rotate_x() function ........................."""
        e = EMData()
        e.set_size(32,32,1)
        e.to_one()
        
        e.rotate_x(10)
        
        #apply only to 2D image
        e2 = EMData()
        e2.set_size(32,32,32)
        e2.to_one()
        self.assertRaises( RuntimeError, e2.rotate_x, 10)
        try:
            e2.rotate_x(10)
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "ImageDimensionException")
            
    def test_rotate_180(self):
        """test rotate_180() function ......................."""
        e = EMData()
        e.set_size(32,32,1)
        e.to_one()
        
        e.rotate_180()
        
        #apply only to square 2D image
        e2 = EMData()
        e2.set_size(32,24,1)
        e2.to_one()
        self.assertRaises( RuntimeError, e2.rotate_180, )
        try:
            e2.rotate_180()
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "ImageFormatException")
            
        #exception to 3D image
        e3 = EMData()
        e3.set_size(32,32,32)
        e3.to_one()
        self.assertRaises( RuntimeError, e3.rotate_180, )
        try:
            e3.rotate_180()
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "ImageDimensionException")

    def test_dot_rotate_translate(self):
        """test dot_rotate_translate() functon .............."""
        e =EMData()
        e.set_size(32,32,1)
        e.process_inplace("testimage.noise.uniform.rand")
        e2 = EMData()
        e2.set_size(32,32,1)
        e2.process_inplace("testimage.noise.uniform.rand")
        e.dot_rotate_translate(e2, 2.0, 3.0, 1.0)
        
        #two image must be the same size
        e3 =EMData()
        e3.set_size(32,32,1)
        e3.process_inplace("testimage.noise.uniform.rand")
        e4 = EMData()
        e4.set_size(24,24,1)
        e4.process_inplace("testimage.noise.uniform.rand")
        self.assertRaises( RuntimeError, e3.dot_rotate_translate, e4, 2.0, 3.0, 1.0)
        try:
            e3.dot_rotate_translate(e4, 2.0, 3.0, 1.0)
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "ImageFormatException")

        #two image must be 2D
        e5 =EMData()
        e5.set_size(32,32,32)
        e5.process_inplace("testimage.noise.uniform.rand")
        e6 = EMData()
        e6.set_size(32,32,32)
        e6.process_inplace("testimage.noise.uniform.rand")
        self.assertRaises( RuntimeError, e5.dot_rotate_translate, e6, 2.0, 3.0, 1.0)
        try:
            e5.dot_rotate_translate(e6, 2.0, 3.0, 1.0)
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "ImageDimensionException")
        
    def test_little_big_dot(self):
        """test little_big_dot() function ..................."""
        big = EMData()
        big.set_size(256,256,1)
        big.process_inplace("testimage.noise.uniform.rand")
        small = EMData()
        small.set_size(32,32,1)
        small.process_inplace("testimage.noise.uniform.rand")
        
        e = big.little_big_dot(small)
        self.assertNotEqual(e, None)
        e2 = big.little_big_dot(small, True)
        self.assertNotEqual(e2, None)
        
        #this image only can be 1D/2D
        e3 = EMData()
        e3.set_size(256,256,256)
        e3.to_one()
        self.assertRaises( RuntimeError, e3.little_big_dot, small)
        try:
            e3.little_big_dot(small)
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "ImageDimensionException")

    def test_do_radon(self):
        """test do_radon() function ........................."""
        e = EMData()
        e.set_size(32,32,1)
        e.process_inplace("testimage.noise.uniform.rand")
        
        e2 = e.do_radon()
        
        #this function only apply to square 2D image
        e3 = EMData()
        e3.set_size(32,24,1)
        e3.process_inplace("testimage.noise.uniform.rand")
        self.assertRaises( RuntimeError, e3.do_radon,)
        try:
            e3.do_radon()
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "ImageFormatException")
        
        e4 = EMData()
        e4.set_size(32,32,32)
        e4.process_inplace("testimage.noise.uniform.rand")
        self.assertRaises( RuntimeError, e4.do_radon,)
        try:
            e4.do_radon()
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "ImageDimensionException")

    def test_calc_ccf(self):
        """test calc_ccf() function ........................."""
        #for two 1 D images
        e = EMData()
        e.set_size(24,1,1)
        e.process_inplace("testimage.noise.uniform.rand")
        e2 = EMData()
        e2.set_size(24,1,1)
        e2.process_inplace("testimage.noise.uniform.rand")
        re = e.calc_ccf(e2)
        self.assertEqual(re.is_complex(), False)
        
        #for two 2 D images
        e3 = EMData()
        e3.set_size(24,24,1)
        e3.process_inplace("testimage.noise.uniform.rand")
        e4 = EMData()
        e4.set_size(24,24,1)
        e4.process_inplace("testimage.noise.uniform.rand")
        re = e3.calc_ccf(e4)
        self.assertEqual(re.is_complex(), False)
        
        #for two 3 D images
        e5 = EMData()
        e5.set_size(24,24,24)
        e5.process_inplace("testimage.noise.uniform.rand")
        e6 = EMData()
        e6.set_size(24,24,24)
        e6.process_inplace("testimage.noise.uniform.rand")
        re = e5.calc_ccf(e6)
        self.assertEqual(re.is_complex(), False)
        
        #try different fpflag other than default CIRCULANT
        #import EMAN
        #from EMAN2 import fundamentals
        re = e5.calc_ccf(e6, fp_flag.CIRCULANT_NORMALIZED)
        re = e5.calc_ccf(e6, fp_flag.PADDED)    
        re = e5.calc_ccf(e6, fp_flag.PADDED_NORMALIZED)    
        re = e5.calc_ccf(e6, fp_flag.PADDED_LAG)    
        re = e5.calc_ccf(e6, fp_flag.PADDED_NORMALIZED_LAG) 
        
    def no_test_calc_ccfx(self):
        """test calc_ccfx() function ........................"""
        e = EMData()
        e.set_size(24,24,1)
        e.process_inplace("testimage.noise.uniform.rand")
        e2 = EMData()
        e2.set_size(24,24,1)
        e2.process_inplace("testimage.noise.uniform.rand")
        re = e.calc_ccfx(e2)
        self.assertEqual(re.is_complex(), False)
        self.assertEqual(re.get_xsize(), 24)
        self.assertEqual(re.get_ysize(), 1)
        self.assertEqual(re.get_zsize(), 1)
        
        Log.logger().set_level(-1)
        
        #test non-default parameter
        re = e.calc_ccfx(e2, 3, 20, True)
        self.assertEqual(re.is_complex(), False)
        self.assertEqual(re.get_xsize(), 24)
        #self.assertEqual(re.get_ysize(), 1)    #fail here, need investigate further
        self.assertEqual(re.get_zsize(), 1)
        
        #the input image can not be None
        e3 = None
        self.assertRaises( RuntimeError, e.calc_ccfx, e3)
        try:
            e.calc_ccfx(e3)
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "NullPointerException")
        
        #two image must be same size
        e4 = EMData()
        e4.set_size(32,32,1)
        self.assertRaises( RuntimeError, e.calc_ccfx, e4)
        try:
            e.calc_ccfx(e4)
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "ImageFormatException")
            
        #this funtion can not be used by 3D image
        e5 = EMData()
        e5.set_size(24,24,24)
        e6 = EMData()
        e6.set_size(24,24,24)
        self.assertRaises( RuntimeError, e5.calc_ccfx, e6)
        try:
            e5.calc_ccfx(e6)
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "ImageDimensionException")
            
    def test_calc_mutual_correlation(self):
        """test calc_mutual_correlation() function .........."""
        e = EMData()
        e.set_size(24,24,24)
        e.process_inplace("testimage.noise.uniform.rand")
        e2 = EMData()
        e2.set_size(24,24,24)
        e2.process_inplace("testimage.noise.uniform.rand")
        
        #test for default argument
        e3 = e.calc_mutual_correlation(e2)
        
        #test for non-default argument
        e4 = e2.do_fft()
        e4.ap2ri()
        e.calc_mutual_correlation(e2, True, e4)
        
        #two image must be in same size
        e5 = EMData()
        e5.set_size(24,24,24)
        e5.process_inplace("testimage.noise.uniform.rand")
        e6 = EMData()
        e6.set_size(32,32,32)
        e6.process_inplace("testimage.noise.uniform.rand")
        self.assertRaises( RuntimeError, e5.calc_mutual_correlation, e6)
        try:
            e5.calc_mutual_correlation(e6)
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "ImageFormatException")
            
    def test_unwrap(self):
        """test unwrap() function ..........................."""
        e = EMData()
        e.set_size(32,32,1)
        e.process_inplace("testimage.noise.uniform.rand")
        
        #test for default argument
        e2 = e.unwrap()
        self.assertNotEqual(e2, None) 
        
        #test non-default arguments
        e3 = e.unwrap(1,1,0,1,1,True)
        self.assertNotEqual(e3, None)
        
        #this function only apply to 2D image
        e4 = EMData()
        e4.set_size(24,24,24)
        self.assertRaises( RuntimeError, e4.unwrap,)
        try:
            e4.unwrap()
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "ImageDimensionException")
            
    def test_mean_shrink(self):
        """test mean_shrink() function ......................"""
        e = EMData()
        e.set_size(64,64,64)
        e.process_inplace("testimage.noise.uniform.rand")
        
        e.mean_shrink(4)
        self.assertEqual(e.get_xsize(), 16)
        self.assertEqual(e.get_xsize(), 16)
        self.assertEqual(e.get_xsize(), 16)
        
        e2 = EMData()
        e2.set_size(30,30,1)
        e2.process_inplace("testimage.noise.uniform.rand")
        e2.mean_shrink(1.5)    
        
        #shrink factor 1.5 only support 2D image
        self.assertRaises( RuntimeError, e.mean_shrink, 1.5)
        try:
            e.mean_shrink(1.5)
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "InvalidValueException")
            
        #shrink factor must be >1
        self.assertRaises( RuntimeError, e.mean_shrink, 0)
        try:
            e.mean_shrink(0)
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "InvalidValueException")

    def test_median_shrink(self):
        """test median_shrink() function ...................."""
        e = EMData()
        e.set_size(64,64,64)
        e.process_inplace("testimage.noise.uniform.rand")
        
        e.median_shrink(4)
        self.assertEqual(e.get_xsize(), 16)
        self.assertEqual(e.get_xsize(), 16)
        self.assertEqual(e.get_xsize(), 16)
        
        #shrink factor must be >1 
        self.assertRaises( RuntimeError, e.median_shrink, 0)
        try:
            e.median_shrink(0)
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "InvalidValueException")
        
        #image size must be divisible by shrink factor
        self.assertRaises( RuntimeError, e.median_shrink, 5)
        try:
            e.median_shrink(5)
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "InvalidValueException")

    def test_apply_radial_func(self):
        """test apply_radial_func() functon ................."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace("testimage.noise.uniform.rand")
        self.assertEqual(e.is_complex(), False)
        e.do_fft_inplace()
        self.assertEqual(e.is_complex(), True)
        e.apply_radial_func(0, 2, (1.0,2.0,3.0))
        e.apply_radial_func(0, 2, (1.0,2.0,3.0), False)
        
    def test_add_incoherent(self):
        """test add_incoherent() function ..................."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace("testimage.noise.uniform.rand")
        e2 = EMData()
        e2.set_size(32,32,32)
        e2.process_inplace("testimage.noise.uniform.rand")
        e3 = e.do_fft()
        e4 = e2.do_fft()
        e3.add_incoherent(e4)
        
        #two images must be the same size
        e5 = EMData()
        e5.set_size(24,24,24)
        e5.process_inplace("testimage.noise.uniform.rand")
        e6 = e5.do_fft()
        self.assertRaises( RuntimeError, e3.add_incoherent, e6)
        try:
            e3.add_incoherent(e6)
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "ImageFormatException")
        
        #two image must be complex
        self.assertRaises( RuntimeError, e.add_incoherent, e2)
        try:
            e.add_incoherent(e2)
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "ImageFormatException")
        
    def test_calc_fourier_shell_correlation(self):
        """test calc_fourier_shell_correlation() ............"""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace("testimage.noise.uniform.rand")
        e2 = EMData()
        e2.set_size(32,32,32)
        e2.process_inplace("testimage.noise.uniform.rand")
        e.calc_fourier_shell_correlation(e2)
        
        #input image can not be null
        e3 = None
        self.assertRaises( RuntimeError, e.calc_fourier_shell_correlation, e3)
        try:
            e.calc_fourier_shell_correlation(e3)
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "NullPointerException")
            
        #two image must be in same size
        e4 = EMData()
        e4.set_size(24,24,24)
        self.assertRaises( RuntimeError, e.calc_fourier_shell_correlation, e4)
        try:
            e.calc_fourier_shell_correlation(e4)
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "ImageDimensionException")
            
    def test_calc_hist(self):
        """test calc_hist() function ........................"""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace("testimage.noise.uniform.rand")
        hist = e.calc_hist()
        
    def test_calc_az_dist(self):
        """test calc_az_dist() function ....................."""
        e = EMData()
        e.set_size(32,32,1)    #this function apply to 2D only
        e.process_inplace("testimage.noise.uniform.rand")
        data = e.calc_az_dist(32, 0.0, 0.1, 0.0, 2.0)
        
        #this function not apply to 3D image
        e2 = EMData()
        e2.set_size(32,32,32)
        e2.process_inplace("testimage.noise.uniform.rand")
        self.assertRaises( RuntimeError, e2.calc_az_dist, 32, 0.0, 0.1, 0.0, 2.0)
        try:
            e2.calc_az_dist(32, 0.0, 0.1, 0.0, 2.0)
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "ImageDimensionException")
        
    def test_calc_dist(self):
        """test calc_dist() function ........................"""
        e = EMData()
        e.set_size(32,1,1)
        e.process_inplace("testimage.noise.uniform.rand")
        e2 = EMData()
        e2.set_size(32,32,32)     #second image must be the same x size as this image
        e2.process_inplace("testimage.noise.uniform.rand")
        dist = e.calc_dist(e2)
        #test for non-default argument
        dist2 = e.calc_dist(e2, 10)
        
        e3 = EMData()
        e3.set_size(32,32,1)
        self.assertRaises( RuntimeError, e3.calc_dist, e2)
        try:
            e3.calc_dist(e2)
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "ImageDimensionException")
            
        e4 = EMData()
        e4.set_size(24,1,1)
        self.assertRaises( RuntimeError, e4.calc_dist, e2)
        try:
            e4.calc_dist(e2)
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "ImageFormatException")
            
    def no_test_calc_flcf(self):
        """test calc_flcf() function ........................"""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace("testimage.noise.uniform.rand")
        e2 = EMData()
        e2.set_size(32,32,32)     
        e2.process_inplace("testimage.noise.uniform.rand")
        e3 = e.calc_flcf(e2)
        e3 = e.calc_flcf(e2, 60)
        
        e4 = None
        self.assertRaises( RuntimeError, e.calc_flcf, e4)
        try:
            e.calc_flcf(e4)
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "NullPointerException")
            
    def test_convolute(self):
        """test convolute() function ........................"""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace("testimage.noise.uniform.rand")
        e2 = EMData()
        e2.set_size(32,32,32)     
        e2.process_inplace("testimage.noise.uniform.rand")
        e3 = e.convolute(e2)        
    
    def test_has_ctff(self):
        """test has_ctff() function ........................."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace("testimage.noise.uniform.rand")
        result = e.has_ctff()
        self.assertEqual(result, False)
        
    def test_dot(self):
        """test dot() function .............................."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace("testimage.noise.uniform.rand")
        e2 = EMData()
        e2.set_size(32,32,32)     
        e2.process_inplace("testimage.noise.uniform.rand")
        f = e.dot(e2)
        
    def no_test_common_lines(self):
        """test common_lines() function ....................."""
        e = EMData()
        e.set_size(32,32,1)
        e.process_inplace("testimage.noise.uniform.rand")
        e2 = EMData()
        e2.set_size(32,32,1)     
        e2.process_inplace("testimage.noise.uniform.rand")
        e3 = e2.do_fft()
        e4 = EMData()
        e4.set_size(32,32,1)     
        e4.process_inplace("testimage.noise.uniform.rand")
        e5 = e4.do_fft()
        e.common_lines(e3, e5)
        
    def no_test_common_lines_real(self):
        """test common_lines_real() function ................"""
        e = EMData()
        e.set_size(32,32,1)
        e.process_inplace("testimage.noise.uniform.rand")
        e2 = EMData()
        e2.set_size(32,32,1)     
        e2.process_inplace("testimage.noise.uniform.rand")
        e3 = EMData()
        e3.set_size(32,32,1)     
        e3.process_inplace("testimage.noise.uniform.rand")
        e.common_lines_real(e2,e3)
        
    def test_cut_slice(self):
        """test cut_slice() function ........................"""
        e = EMData()
        e.set_size(32,32,1)
        e2 = EMData()
        e2.set_size(32,32,32)     
        e2.process_inplace("testimage.noise.uniform.rand")
        e.cut_slice(e2, 10.5)    #default argument
        e.cut_slice(e2, 11.6, None, False, 1, 1)
        
        e3 = None
        self.assertRaises( RuntimeError, e.cut_slice, e3, 5.1)
        try:
            e.cut_slice(e3, 5.1)
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "NullPointerException")
            
    def test_uncut_slice(self):
        """test uncut_slice() function ......................"""
        e = EMData()
        e.set_size(32,32,1)
        e.process_inplace("testimage.noise.uniform.rand")
        e2 = EMData()
        e2.set_size(32,32,32)     
        e2.process_inplace("testimage.noise.uniform.rand")
        e.uncut_slice(e2, 10.0)
        e.uncut_slice(e2, 10.0, None, 2, 3) #non-default argument
        
    def test_calc_center_density(self):
        """test calc_center_density() function .............."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace("testimage.noise.uniform.rand")
        d = e.calc_center_density()
        self.assertNotEqual(d, None)
        
        e2 = EMData()
        e2.set_size(24,24,24)
        d2 = e2.calc_center_density()
        self.assertEqual(d2, 0)
        
    def test_calc_sigma_diff(self):
        """test calc_sigma_diff() function .................."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace("testimage.noise.uniform.rand")
        s = e.calc_sigma_diff()
        self.assertNotEqual(s, None)
        
    def test_calc_min_location(self):
        """test calc_min/max_location() function ............"""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace("testimage.noise.uniform.rand")
        loc = e.calc_min_location()
        self.assertNotEqual(loc, None)
        loc2 = e.calc_max_location()
        self.assertNotEqual(loc2, None)
        index = e.calc_min_index()
        self.assertEqual(index, loc[0]+loc[1]*32+loc[2]*32*32)
        index2 = e.calc_max_index()
        self.assertEqual(index2, loc2[0]+loc2[1]*32+loc2[2]*32*32)
        
    def test_get_edge_mean(self):
        """test get_edge/circle_mean() function ............."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace("testimage.noise.uniform.rand")
        em = e.get_edge_mean()
        cm = e.get_circle_mean()
    
    def test_project(self):
        """test image projection ............................"""
        n = 20
        infile = "test_project.mrc"
        TestUtil.make_image_file(infile, MRC, EM_FLOAT, n, n, n)
        volume = EMData()
        volume.read_image(infile)
        pi = math.pi
        proj = volume.project("standard", { "alt" : pi/3, "az" : pi/5, "phi" : 1})
        self.assertEqual(proj.get_xsize(), n)
        self.assertEqual(proj.get_ysize(), n)
        self.assertEqual(proj.get_zsize(), 1)
        testlib.check_emdata(proj, sys.argv[0])
        os.unlink(infile)

    def test_calc_highest_locations(self):
        """test calculation of highest location ............."""
        infile = "test_calc_highest_locations.mrc"
        TestUtil.make_image_file(infile, MRC, EM_FLOAT, 40, 60)
        
        e = EMData()
        e.read_image(infile)
        pixels = e.calc_highest_locations(400)
        self.assertEqual(len(pixels),  1144)
        self.assertEqual(pixels[0], Pixel(19, 10, 0, 401.0))

        os.unlink(infile)

    def test_get_attr_dict(self):
        """test get_attr_dict() function ...................."""
        imgfile = "tablet.mrc"
        TestUtil.make_image_file(imgfile, MRC, EM_FLOAT, 40, 60)
        
        e = EMData()
        e.read_image(imgfile)
        d0 = e.get_attr_dict()
        e += 10
        d1 = e.get_attr_dict()        
        e += 100
        d2 = e.get_attr_dict()

        cur_attrlist = []        
        for mydict in (d0, d1, d2):
            for mykey in mydict.keys():
                if (not ("MRC" in mykey)):
                    cur_attrlist.append(mykey + "=" + str(mydict[mykey])+"\n")        
        
        os.unlink(imgfile)        

    def test_set_attr_dict(self):
        """test set/del_attr_dict() function ................"""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace("testimage.noise.uniform.rand")
        d = e.get_attr_dict()
        self.assertEqual(d['is_complex'], False)
        d['is_complex'] = True
        e.set_attr_dict(d)
        d2 = e.get_attr_dict()
        self.assertEqual(d2['is_complex'], True)
        
        e.set_attr('name1', 'Grant')
        e.set_attr('int1', 1000)
        d3 = e.get_attr_dict()
        self.assertEqual(d3['name1'], 'Grant')
        self.assertEqual(d3['int1'], 1000)
        e.del_attr_dict(['name1', 'int1', 'is_complex'])
        d4 = e.get_attr_dict()
        self.assertEqual(d4.has_key('name1'), False)
        self.assertEqual(d4.has_key('int1'), False)
        self.assertEqual(d4.has_key('is_complex'), False)
        self.assertEqual(d4.has_key('mean'), True)
        
    def test_get_clip(self):
        """test get_clip() function ........................."""
        nx = 32
        ny = 48
        filebase = "test_get_clip_" + str(os.getpid())
        infile = filebase + ".mrc"
        TestUtil.make_image_file(infile, MRC, EM_FLOAT, nx, ny)

        e = EMData()
        e.read_image(infile)

        region1 = Region(nx/4, ny/4, nx/2, ny/2)
        outfile1 = filebase + "_out1.mrc"
        e2 = e.get_clip(region1)
        e2.write_image(outfile1)

        self.assertEqual(e2.get_xsize(), nx/2)
        self.assertEqual(e2.get_ysize(), ny/2)
        
        region2 = Region(-nx/4, -ny/4, 2*nx, 2*ny)
        outfile2 = filebase + "_out2.mrc"
        e3 = e.get_clip(region2)
        e3.write_image(outfile2)
        
        self.assertEqual(e3.get_xsize(), nx*2)
        self.assertEqual(e3.get_ysize(), ny*2)

        os.unlink(infile)
        os.unlink(outfile1)
        os.unlink(outfile2)    
            
    def test_get_rotated_clip(self):
        """test get_rotated_clip() function ................."""
        a=EMData()
        a.set_size(100,100,100)
        a.to_one()
        b=a.get_rotated_clip(Transform3D([24,24,24], 0,0,0),[32,32,32],1.0)
        
    def test_complex_image(self):
        """test complex image ..............................."""
        nx = 16
        ny = 32
        infile = "test_complex_image_1.mrc"
        TestUtil.make_image_file(infile, MRC, EM_FLOAT_COMPLEX, nx, ny)

        e = EMData()
        e.read_image(infile)
        attrd = e.get_attr_dict()

        self.assertEqual(attrd["nx"], nx+2)
        self.assertEqual(attrd["ny"], ny)
        self.assertEqual(attrd["nz"], 1)
        self.assertEqual(attrd["datatype"], EM_FLOAT_COMPLEX)

        os.unlink(infile)        
        
    def test_set_value_at(self):
        """test set_value_at() .............................."""
        nx = 10
        ny = 20
        nz = 2
        infile = "test_set_value_at_in.mrc"
        TestUtil.make_image_file(infile, MRC, EM_FLOAT, nx, ny, nz)
        
        e = EMData()
        e.read_image(infile)

        for i in range(nz):
            for j in range(ny):
                for k in range(nx):
                    e.set_value_at(k,j,i, 1)

        os.unlink(infile)
        narray = EMNumPy.em2numpy(e)

        for i in range(nz):
            for j in range(ny):
                for k in range(nx):
                    self.assertEqual(e.get_value_at(k,j,i), 1)
                    self.assertEqual(e.sget_value_at(k,j,i), 1)
                    self.assertEqual(narray[i][j][k], 1)      
                    
        #test exception, if set value out of range
        e2 = EMData()
        e2.set_size(32,32,32)
        e2.process_inplace("testimage.noise.uniform.rand")
        self.assertRaises( RuntimeError, e2.set_value_at, 32, 1, 1, 1.0)
        try:
            e2.set_value_at(32, 1, 1, 1.0)
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "OutofRangeException")
        self.assertRaises( RuntimeError, e2.set_value_at, 31, 45, 1, 1.0)
        try:
            e2.set_value_at(31, 45, 1, 1.0)
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "OutofRangeException")
        self.assertRaises( RuntimeError, e2.set_value_at, 31, 1, 1000, 1.0)
        try:
            e2.set_value_at(31, 1, 1000, 1.0)
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "OutofRangeException")
            
    def test_sget_value_at_interp(self):
        """test sget_value_at_interp() function ............."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace("testimage.noise.uniform.rand")
        d = e.sget_value_at_interp(1.5, 1.5, 1.5)
        
        e2 = EMData()
        e2.set_size(32,32,32)
        e2.process_inplace("testimage.noise.uniform.rand")
        d2 = e2.sget_value_at_interp(1.5, 1.5)
                                  
    def test_calc_radial_dist(self):
        """test calc_radial_dist() function.................."""
        file1 = "test_calc_radial_dist.mrc"
        TestUtil.make_image_file(file1, MRC)
        
        e1 = EMData()
        e1.read_image(file1)
        ny = e1.get_ysize()
        dlist = e1.calc_radial_dist(ny, 0, 0.5)

        dlist2 = [0.0, 0.5732, 1.1464, 2.0, 4.5135, 5.7996, 9.2307, 12.1255, 16.4826, 19.1460, 25.6796, 30.2887, 36.4654, 41.4641, 49.9667, 1538.3833]

        for i in range(len(dlist)):
            testlib.assertfloat(self, dlist[i], dlist2[i])

        os.unlink(file1)

    def test_rotavg(self):
        """test rotavg() function ..........................."""
        e = EMData()
        e.set_size(32,32,1)
        e.process_inplace("testimage.noise.uniform.rand")
        
        re = e.rotavg()
        self.assertEqual(re.get_ysize(), 1)
        self.assertEqual(re.get_zsize(), 1)
        
        #this function not apply to 1D image
        e2 = EMData()
        e2.set_size(32,1,1)
        self.assertRaises( RuntimeError, e2.rotavg, )
        try:
            e2.rotavg()
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "ImageDimensionException")
            
    def test_image_arithmetic(self):
        """test image arithmetic operation .................."""
        #test add
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace("testimage.noise.uniform.rand")
        e2 = e + 2.0
        d = e.get_3dview()
        d2 = e2.get_3dview()
        for x in range(32):
            for y in range(32):
                for z in range(32):
                    self.assertAlmostEqual(d[z][y][x]+2.0, d2[z][y][x], 3)
        e3 = EMData()
        e3.set_size(32,32,32)
        e3.to_one()
        e4 = e2 + e3
        d4 = e4.get_3dview()
        for x in range(32):
            for y in range(32):
                for z in range(32):
                    self.assertAlmostEqual(d2[z][y][x]+1.0, d4[z][y][x], 3)
        #image must be the same size
        e5 = EMData()
        e5.set_size(24,24,24)
        self.assertRaises( RuntimeError, e.add, e5)
        try:
            e6 = e + e5
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "ImageFormatException")
        #image must be both real or both complex
        e6 = EMData()
        e6.set_size(32,32,32)
        e6.set_complex(True)
        self.assertRaises( RuntimeError, e.add, e6)
        try:
            e += e6
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "ImageFormatException")
        
        #test substract
        e7 = e - 3.0
        d7 = e7.get_3dview()
        for x in range(32):
            for y in range(32):
                for z in range(32):
                    self.assertAlmostEqual(d[z][y][x]-3.0, d7[z][y][x], 3)
        #image must be the same size
        self.assertRaises( RuntimeError, e.sub, e5)
        try:
            e -= e5
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "ImageFormatException")
        #image must be both real or both complex   
        self.assertRaises( RuntimeError, e.sub, e6) 
        try:
            e -= e6
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "ImageFormatException")
            
        #test mutiply
        e8 = e * 2
        d8 = e8.get_3dview()
        for x in range(32):
            for y in range(32):
                for z in range(32):
                    self.assertAlmostEqual(d[z][y][x] * 2.0, d8[z][y][x], 3)
        #image must be the same size
        self.assertRaises( RuntimeError, e.mult, e5)
        try:
            e *= e5
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "ImageFormatException")
        #image must be both real or both complex   
        self.assertRaises( RuntimeError, e.mult, e6)
        try:
            e *= e6
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "ImageFormatException") 
            
        #test division
        e9 = e / 5.0
        d9 = e9.get_3dview()
        for x in range(32):
            for y in range(32):
                for z in range(32):
                    self.assertAlmostEqual(d[z][y][x]/5.0, d9[z][y][x], 3)
        #image must be the same size
        self.assertRaises( RuntimeError, e.div, e5)
        try:
            e /= e5
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "ImageFormatException")
        #image must be both real or both complex   
        self.assertRaises( RuntimeError, e.div, e6)
        try:
            e /= e6
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "ImageFormatException")
        #division by zero exception
        self.assertRaises( RuntimeError, e.div, 0.0 )
        try:
            e/= 0.0
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "InvalidValueException")
        e7 = EMData()
        e7.set_size(32,32,32)
        e7.to_zero()
        self.assertRaises( RuntimeError, e.div, e7 )
        try:
            e/=e7
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "InvalidValueException")

    def test_stat_locations(self):
        """test locational stats ............................"""
        nx = 16
        ny = 24
        file1 = "test_stat_locations.mrc"
        TestUtil.make_image_file(file1, MRC, EM_FLOAT, nx, ny)
        
        e1 = EMData()
        e1.read_image(file1)

        max_index = e1.calc_max_index()
        self.assertEqual(max_index, 0)

        max_location = e1.calc_max_location()
        self.assertEqual(max_location, (0,0,0))

        max_index2 = max_location[2] * nx * ny + max_location[1] * nx + max_location[0]
        self.assertEqual(max_index, max_index2)
        
        attrdict = e1.get_attr_dict()
        max = attrdict["maximum"]
        self.assertEqual(max, 208)

        max2 = e1.get_value_at(max_location[2], max_location[1], max_location[0])
        self.assertEqual(max, max2)
        
        os.unlink(file1)

    def test_image_overwrite(self):
        """test overwriting a image ........................."""
        file1 = "test_image_overwrite.mrc"
        nx = 16
        ny = 18
        nz = 2
        TestUtil.make_image_file(file1, MRC, EM_FLOAT, nx, ny, nz)
        e = EMData()
        e.read_image(file1)
        
        nx2 = nx*2
        ny2 = ny*2
        nz2 = nz
        e.set_size(nx2, ny2, nz2)
        e.write_image(file1)
        
        e2 = EMData()
        e2.read_image(file1)
        self.assertEqual(nx2, e2.get_xsize())
        self.assertEqual(ny2, e2.get_ysize())
        self.assertEqual(nz2, e2.get_zsize())

        os.unlink(file1)
        
    def test_ctf(self):
        """test ctf ........................................."""
        infile = "test_ctf_in.mrc"
        TestUtil.make_image_file(infile, MRC)
        ctf = SimpleCtf()
        d = {"defocus":1, "bfactor":2, "amplitude":3, "ampcont":4, "noise1":5, "noise2":6, "noise3":7, "noise4":8, "voltage":9, "cs":10,"apix":11}
        
        ctf.from_dict(d)

        e = EMData()
        e.read_image(infile)
        e.set_ctf(ctf)

        outfile = "test_ctf_out.mrc"
        e.write_image(outfile)

        e2 = EMData()
        e2.read_image(outfile)
        ctf2 = e2.get_ctf()

        ctfstr = ctf2.to_string()
        self.assertEqual(ctfstr, "1 2 3 4 5 6 7 8 9 10 11")
        os.unlink(infile)
        os.unlink(outfile)
        
    def test_get_translation(self):
        """test get/set_translation() function .............."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace("testimage.noise.uniform.rand")
        e.set_translation(1.0,2.0,3.0)
        t = e.get_translation()
        self.assertEqual(t.at(0), 1.0)
        self.assertEqual(t.at(1), 2.0)
        self.assertEqual(t.at(2), 3.0)
        
        e.set_translation(Vec3f(2.0,3.0,4.0))
        t1 = e.get_translation()
        self.assertEqual(t1.at(0), 2.0)
        self.assertEqual(t1.at(1), 3.0)
        self.assertEqual(t1.at(2), 4.0)
    
    def test_get_transform(self):
        """test get_transform() function ...................."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace("testimage.noise.uniform.rand")
        t = e.get_transform()    
    
    def test_set_rotation(self):
        """test set_rotation() function ....................."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace("testimage.noise.uniform.rand")
        e.set_rotation(0.23, 0.45, 0.67)
    
    def test_set_size(self):
        """test set_size() function ........................."""
        e = EMData()
        e.set_size(32)    #test default argument
        self.assertEqual(e.get_xsize(), 32)
        self.assertEqual(e.get_ysize(), 1)
        self.assertEqual(e.get_zsize(), 1)
        e.set_size(23,34,56)    #test non-default argument
        self.assertEqual(e.get_xsize(), 23)
        self.assertEqual(e.get_ysize(), 34)
        self.assertEqual(e.get_zsize(), 56)
        
        #test for exception
        self.assertRaises( RuntimeError, e.set_size, 0, 32, 32)
        try:
            e.set_size(0, 32, 32)
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "InvalidValueException")
            
        self.assertRaises( RuntimeError, e.set_size, 32, -32, 32)
        try:
            e.set_size(32, -32, 32)
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "InvalidValueException")
            
        self.assertRaises( RuntimeError, e.set_size, 32, 32, -32)
        try:
            e.set_size(32, -32, 32)
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "InvalidValueException")
            
        #test set_complex_size()
        e2 = EMData()
        e2.set_size(32)
        e2.process_inplace("testimage.noise.uniform.rand")
        e3 = e2.do_fft()
        e3.set_complex_size(32)
        self.assertEqual(e3.get_xsize(), 64)
        self.assertEqual(e3.get_ysize(), 1)
        self.assertEqual(e3.get_zsize(), 1)
        
    def test_set_path(self):
        """test set_path() function ........................."""
        e = EMData()
        e.set_size(32,32,32)
        e.set_path("new_path")
        e.set_pathnum(10) 
        
    def test_get_row(self):
        """test get_row() function .........................."""
        e = EMData()
        e.set_size(32, 32)
        e.process_inplace("testimage.noise.uniform.rand")
        e2 = e.get_row(1)
        self.assertEqual(e2.get_xsize(), 32)
        self.assertEqual(e2.get_ysize(), 1)
        self.assertEqual(e2.get_zsize(), 1)
        d = e.get_2dview()
        d2 = e2.get_3dview()
        for x in range(32):
            self.assertEqual(d[1][x], d2[0][0][x])
            
        #this function only apply to 1D/2D image
        e3 = EMData()
        e3.set_size(32,32,32)
        self.assertRaises( RuntimeError, e3.get_row, 1) 
        try:
            e3.get_row(1)
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "ImageDimensionException")
            
    def test_set_row(self):
        """test set_row() function .........................."""
        e = EMData()
        e.set_size(32, 32)
        e.process_inplace("testimage.noise.uniform.rand")
        e2 = EMData()
        e2.set_size(32,1)
        e.set_row(e2, 1)
        self.assertEqual(e.get_xsize(), 32)
        self.assertEqual(e.get_ysize(), 32)
        
        d = e.get_2dview()
        d2 = e2.get_3dview()
        for x in range(32):
            self.assertEqual(d[1][x], d2[0][0][x])
            
        #this function only apply to 1D/2D image
        e3 = EMData()
        e3.set_size(32,32,32)
        self.assertRaises( RuntimeError, e3.set_row, e2, 1)
        try:
            e3.set_row(e2, 1)
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "ImageDimensionException")
            
    def test_get_col(self):
        """test get_col() function .........................."""
        e = EMData()
        e.set_size(32, 32)
        e.process_inplace("testimage.noise.uniform.rand")
        e2 = e.get_col(1)
        self.assertEqual(e2.get_xsize(), 32)
        self.assertEqual(e2.get_ysize(), 1)
        self.assertEqual(e2.get_zsize(), 1)
        d = e.get_2dview()
        d2 = e2.get_3dview()
        for y in range(32):
            self.assertEqual(d[y][1], d2[0][0][y])
            
        #this function only apply to 2D image
        e3 = EMData()
        e3.set_size(32,32,32)
        self.assertRaises( RuntimeError, e3.get_col, 1) 
        try:
            e3.get_col(1)
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "ImageDimensionException")
            
    def test_set_col(self):
        """test set_col() function .........................."""
        e = EMData()
        e.set_size(32, 32)
        e.process_inplace("testimage.noise.uniform.rand")
        e2 = EMData()
        e2.set_size(32,1)
        e.set_col(e2, 1)
        
        d = e.get_2dview()
        d2 = e2.get_3dview()
        for y in range(32):
            self.assertEqual(d[y][1], d2[0][0][y])
            
        #this function only apply to 2D image
        e3 = EMData()
        e3.set_size(32,32,32)
        self.assertRaises( RuntimeError, e3.set_col, e2, 1) 
        try:
            e3.set_col(e2, 1)
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "ImageDimensionException")
    
    def test_set_attr(self):
        """test set/get/del_attr() function ................."""
        e = EMData()
        e.set_size(32, 32, 32)
        e.process_inplace("testimage.noise.uniform.rand")
        self.assertEqual(e.get_attr("is_complex"), False)
        e.set_attr("is_complex", True)
        self.assertEqual(e.get_attr("is_complex"), True)
        
        e.del_attr('is_complex')
        d = e.get_attr_dict()
        self.assertEqual(d.has_key('is_complex'), False)
        
    def test_boolean_check(self):
        """test some boolean check in EMData ................"""
        e = EMData()
        e.set_size(32, 32, 32)
        e.process_inplace("testimage.noise.uniform.rand")

        self.assertEqual(e.is_shuffled(), False)
        e.set_shuffled(True)
        self.assertEqual(e.is_shuffled(), True)
        
        self.assertEqual(e.is_FH(), False)
        e.set_FH(True)
        self.assertEqual(e.is_FH(), True)
        
        self.assertEqual(e.is_complex(), False)
        e.set_complex(True)
        self.assertEqual(e.is_complex(), True)
        
        self.assertEqual(e.is_real(), False)
        e.set_complex(False)
        self.assertEqual(e.is_real(), True)
        
        self.assertEqual(e.is_complex_x(), False)
        e2 = EMData()
        e2.set_size(32, 1, 1)
        e2.process_inplace("testimage.noise.uniform.rand")
        self.assertEqual(e2.is_complex_x(), False)
        e3 = e2.do_fft()
        self.assertEqual(e3.is_complex_x(), True)
        
        self.assertEqual(e.is_flipped(), False)
        e.set_flipped(True)
        self.assertEqual(e.is_flipped(), True)
        
        self.assertEqual(e.is_ri(), True)
        e4 = e.do_fft()
        self.assertEqual(e4.is_ri(), True)
        e4.set_ri(False)
        self.assertEqual(e4.is_ri(), False)
        
        self.assertEqual(e.is_fftpadded(), False)
        e.set_fftpad(True)
        self.assertEqual(e.is_fftpadded(), True)
        
        e5 = EMData()
        e5.set_size(31,31,31)
        e6 = e5.do_fft()
        self.assertEqual(e6.is_fftodd(), True)
        e7 = EMData()
        e7.set_size(32,32,32)
        e8 = e7.do_fft()
        self.assertEqual(e8.is_fftodd(), False)
        
    
    def no_test_statistics(self):
        """test statistics of EMData ........................"""
        e = EMData()
        e.set_size(10,10,1)
        e.process_inplace("testimage.circlesphere",{"radius":4})
        descriptive_statistics(e)
        f = norm_pad_ft(e, False, True, 2)
        descriptive_statistics(f)
        f.do_ift_inplace()
        descriptive_statistics(f)
        g = f*10
        descriptive_statistic   
        
    def test_more_image_arithmetic(self):
        """test some more image arithmetics ................."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        d = e.get_3dview() 
        e2 = e.power(2)
        d2 = e2.get_3dview()
        dd = e.get_3dview()
        for x in range(32):
            for y in range(32):
                for z in range(32):
                    self.assertEqual( d[z][y][x], dd[z][y][x] )
                    self.assertAlmostEqual( dd[z][y][x]**2, d2[z][y][x], 3 ) 
        
        c = e.do_fft()
        c2 = c.real()
        c3 = c.imag()
        self.assertEqual( c2.get_xsize()*2, c.get_xsize() )
        self.assertEqual( c3.get_xsize()*2, c.get_xsize() )
        d = c.get_3dview()
        d2 = c2.get_3dview()
        d3 = c3.get_3dview()
        for x in range(16):
            for y in range(32):
                for z in range(32):
                    self.assertEqual( d2[z][y][x], d[z][y][x*2] )
                    self.assertEqual( d3[z][y][x], d[z][y][x*2+1] )
        self.assertRaises( RuntimeError, e.imag, )
        try:
            ee = e.imag()
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "Invalid")
                    
        e5 = e.real2complex()    #test default argument
        e6 = e.real2complex(1.0)    #test non-default argument
        self.assertEqual(e5.get_xsize(), e.get_xsize()*2)
        self.assertEqual(e6.get_xsize(), e.get_xsize()*2)
        d = e.get_3dview()
        d5 = e5.get_3dview()
        d6 = e6.get_3dview()
        for x in range(32):
            for y in range(32):
                for z in range(32):
                    self.assertEqual( d5[z][y][x*2+1], 0 )
                    self.assertEqual( d6[z][y][x*2+1], 1.0 )
                    self.assertEqual( d5[z][y][x*2], d[z][y][x] )
                    self.assertEqual( d5[z][y][x*2], d[z][y][x] )

    def test_image_amplitude_phase(self):
        """test amplitude()/phase() function ................"""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        e.do_fft_inplace()
        e.ri2ap()
        e2 = e.amplitude()
        e3 = e.phase()
        d = e.get_3dview()
        d2 = e2.get_3dview()
        d3 = e3.get_3dview()
        for x in range(16):
            for y in range(32):
                for z in range(32):
                    self.assertEqual( d2[z][y][x], d[z][y][x*2] )
                    self.assertEqual( d3[z][y][x], d[z][y][x*2+1] )
        
        #amplitude()/phase() only apply to complex image in amplitude/phase format
        f = EMData()
        f.set_size(32,32)
        f.process_inplace('testimage.noise.uniform.rand')
        self.assertRaises( RuntimeError, f.amplitude, )
        self.assertRaises( RuntimeError, f.phase, )
        
        f.do_fft_inplace()
        self.assertRaises( RuntimeError, f.amplitude, )
        self.assertRaises( RuntimeError, f.phase, )
        
    def test_read_images(self):
        """test read_images() function ......................"""
        e = EMData()
        e.set_size(32,32,32)
        
        e.process_inplace('testimage.noise.uniform.rand')
        file1 = 'uniformrand.mrc'
        e.write_image(file1)
        im = EMData.read_images(file1)
        
        os.unlink(file1)
        
    def test_rot_trans2D(self):
        """test rot_trans2D() function ......................"""
        e = EMData()
        e.set_size(32,32,1)
        e.process_inplace('testimage.noise.uniform.rand')
        e2 = e.rot_trans2D(0.12)    #default argument
        e3 = e.rot_trans2D(0.12, 0.2, 0.3)    #non-default argument
        
        #this functin only apply to 2D image
        e4 = EMData()
        e4.set_size(32,1,1)
        self.assertRaises( RuntimeError, e4.rot_trans2D, 0.12)
        try:
            e4.rot_trans2D(0.12)
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "ImageDimensionException")
            
        e5 = EMData()
        e5.set_size(32,32,32)
        self.assertRaises( RuntimeError, e5.rot_trans2D, 0.12)
        try:
            e5.rot_trans2D(0.12)
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "ImageDimensionException")
    
    def test_rot_scale_trans2D(self):
        """test rot_scale_trans2D() function ................"""
        e = EMData()
        e.set_size(32,32,1)
        e.process_inplace('testimage.noise.uniform.rand')
        e2 = e.rot_scale_trans2D(0.12)    #test default argument
        e3 = e.rot_scale_trans2D(0.12, 2.0, 0.2, 0.3)    #test non-default argument
        
        #this function not apply to 1D image
        e4 = EMData()
        e4.set_size(32,1,1)
        self.assertRaises( RuntimeError, e4.rot_scale_trans2D, 0.12)
        try:
            e4.rot_scale_trans2D(0.12)
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "ImageDimensionException")
            
        #this function can not apply to 3D image
        e5 = EMData()
        e5.set_size(32,32,32)
        self.assertRaises( RuntimeError, e5.rot_scale_trans2D, 0.12)
        try:
            e5.rot_scale_trans2D(0.12)
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "ImageDimensionException")
        
            
    def no_test_rotconvtrunc2d_kbi0(self):
        """test rotconvtrunc2d_kbi0() function .............."""
        e = EMData()
        e.set_size(32,32,1)
        e.process_inplace('testimage.noise.uniform.rand')
        e2 = e.rotconvtrunc2d_kbi0(0.1,0.2,10)
        
        #this function only apply to 2D image
        e3 = EMData()
        e3.set_size(32,1,1)
        self.assertRaises( RuntimeError, e3.rotconvtrunc2d_kbi0, 0.1, 0.2, 10)
        try:
            e3.rotconvtrunc2d_kbi0(0.1, 0.2, 10)
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "ImageDimensionException")
            
        e4 = EMData()
        e4.set_size(32,32,32)
        self.assertRaises( RuntimeError, e4.rotconvtrunc2d_kbi0, 0.1, 0.2, 10)
        try:
            e4.rotconvtrunc2d_kbi0(0.1, 0.2, 10)
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "ImageDimensionException")
            
    def no_test_gridrot2d_kbi0(self):
        """test gridrot2d_kbi0() function ..................."""
        e = EMData()
        e.set_size(32,32,1)
        e.process_inplace('testimage.noise.uniform.rand')
        e2 = e.gridrot2d_kbi0(0.1)    #default argument
        e3 = e.gridrot2d_kbi0(0.1, 1.5, 10)
        
    def test_symvol(self):
        """test symvol() function ..........................."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        e2 = e.symvol('CSYM')
        e3 = e.symvol('DSYM')
        #e4 = e.symvol('ICOS_SYM')    #todo problem here, need investigate
        #e5 = e.symvol('OCT_SYM')
        #e6 = e.symvol('ISYM')
    
    def test_fft_shuffle_is_shuffled(self):
        """test fft_shuffle correctly sets is_shuffled......."""
        e = EMData()
        e.set_size(4,4,1)
        e.to_zero()
        e[3,3] = 1.
        eft = e.do_fft()
        self.assertEqual(eft.is_shuffled(), False)
        eft.fft_shuffle()
        self.assertEqual(eft.is_shuffled(), True)
        eft.fft_shuffle()
        self.assertEqual(eft.is_shuffled(), False)

    def test_fft_shuffle_exact_toggle_odd(self):
        """test fft_shuffle() toggles for odd images........."""
        e_odd = EMData()
        e_odd.set_size(3,3,1)
        e_odd.to_zero()
        e_odd[1,1] = 1.
        e_oddft = e_odd.do_fft()
        e_oddft_copy = e_oddft.copy()
        e_oddft_copy.fft_shuffle()
        e_oddft_copy.fft_shuffle()
        diff = e_oddft_copy - e_oddft
        diff2 = diff*diff
        error = False
        if diff2.get_attr("maximum") > 1.e-10:
            error = True
        self.assertEqual(error, False)

    def test_fft_shuffle_exact_toggle_even(self):
        """test fft_shuffle() toggles for even images........"""
        e_even = EMData()
        e_even.set_size(4,4,1)
        e_even.to_zero()
        e_even[1,1] = 1.
        e_evenft = e_even.do_fft()
        e_evenft_copy = e_evenft.copy()
        e_evenft_copy.fft_shuffle()
        e_evenft_copy.fft_shuffle()
        diff = e_evenft_copy - e_evenft
        diff2 = diff*diff
        error = False
        if diff2.get_attr("maximum") > 1.e-10:
            error = True
        self.assertEqual(error, False)

    def test_getitem_real1d(self):
        """Test __getitem__ on a real 1-D image.............."""
        e = EMData()
        e.set_size(4,1,1)
        e.to_zero()
        e += 1.
        val = e[1]
        self.assertEqual(val, 1.)
    def test_getitem_real2d(self):
        """Test __getitem__ on a real 2-D image.............."""
        e = EMData()
        e.set_size(4,4,1)
        e.to_zero()
        e += 1.
        val = e[1,1]
        self.assertEqual(val, 1.)
    def test_getitem_real3d(self):
        """Test __getitem__ on a real 2-D image.............."""
        e = EMData()
        e.set_size(4,4,4)
        e.to_zero()
        e += 1.
        val = e[1,1,1]
        self.assertEqual(val, 1.)
    def test_setitem_real1d(self):
        """Test __setitem__ on a real 1-D image.............."""
        e = EMData()
        e.set_size(4,1,1)
        e.to_zero()
        e[1] = 3.
        val = e[1]
        self.assertEqual(val, 3.)
    def test_setitem_real2d(self):
        """Test __setitem__ on a real 2-D image.............."""
        e = EMData()
        e.set_size(4,4,1)
        e.to_zero()
        e[1,1] = 3.
        val = e[1,1]
        self.assertEqual(val, 3.)
    def test_setitem_real3d(self):
        """Test __setitem__ on a real 3-D image.............."""
        e = EMData()
        e.set_size(4,4,4)
        e.to_zero()
        e[1,1,1] = 3.
        val = e[1,1,1]
        self.assertEqual(val, 3.)
    def test_getitem_complex1d(self):
        """Test __getitem__ on a complex 1-D image..........."""
        e = EMData()
        e.set_size(4,1,1)
        e[3] = 1.
        e.set_complex(True)
        val = e[1]
        self.assertEqual(val, 1j)
    def test_getitem_complex2d(self):
        """Test __getitem__ on a complex 2-D image..........."""
        e = EMData()
        e.set_size(4,4,1)
        e[3,0] = 1.
        e.set_complex(True)
        val = e[1,0]
        self.assertEqual(val, 1j)
    def test_getitem_complex3d(self):
        """Test __getitem__ on a complex 3-D image..........."""
        e = EMData()
        e.set_size(4,4,4)
        e[3,0,0] = 1.
        e.set_complex(True)
        val = e[1,0,0]
        self.assertEqual(val, 1j)
    def test_setitem_complex1d(self):
        """Test __setitem__ on a complex 1-D image..........."""
        e = EMData()
        e.set_size(4,1,1)
        e.set_complex(True)
        e[1] = complex(1.,1.)
        val = e[1]
        self.assertEqual(val, 1+1j)
    def test_setitem_complex2d(self):
        """Test __setitem__ on a complex 2-D image..........."""
        e = EMData()
        e.set_size(4,4,1)
        e.set_complex(True)
        e[1,0] = complex(1.,1.)
        val = e[1,0]
        self.assertEqual(val, 1+1j)
    def test_setitem_complex3d(self):
        """Test __setitem__ on a complex 3-D image..........."""
        e = EMData()
        e.set_size(4,4,4)
        e.set_complex(True)
        e[1,0,0] = complex(1.,1.)
        val = e[1,0,0]
        self.assertEqual(val, 1+1j)

    def test_log(self):
        """test log() arithmatica function for image........."""
        from math import log
        e = EMData()
        e.set_size(2,2,2)
        e.to_one()
        coefficient = 6.0
        e = e*coefficient
        e1 = e.log()
        d = e1.get_3dview()
        for z in range(2):
            for y in range(2):
                for x in range(2):
                    self.assertAlmostEqual( d[z][y][x], log(coefficient))
                    
    def test_log10(self):
        """test log10() arithmatica function for image......."""
        from math import log10
        e = EMData()
        e.set_size(2,2,2)
        e.to_one()
        coefficient = 6.0
        e = e*coefficient
        e1 = e.log10()
        d = e1.get_3dview()
        for z in range(2):
            for y in range(2):
                for x in range(2):
                    self.assertAlmostEqual( d[z][y][x], log10(coefficient))
    
    def test_process(self):
        """test process(Processor)..........................."""
        p = Processors.get('testimage.circlesphere', {'radius':200})
        
        e = EMData()
        e.set_size(800,800)
        e.process_inplace(p)
        e2 = e.process(p)
        
        d1 = e.get_2dview()
        d2 = e2.get_2dview()
        for y in range(800):
            for x in range(800):
                self.assertEqual(d1[y][x], d2[y][x])

def test_main():
    test_support.run_unittest(TestEMData)

if __name__ == '__main__':
    test_main()



