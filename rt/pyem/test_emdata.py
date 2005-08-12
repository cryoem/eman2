#!/bin/env python

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

    def test_default_args(self):
        """test default constructor of EMData ..............."""
        e = EMData()
        self.assertEqual(e.get_attr("apix_x"), 1.0)
        self.assertEqual(e.get_attr("apix_y"), 1.0)
        self.assertEqual(e.get_attr("apix_z"), 1.0)
        self.assertEqual(e.get_attr("is_complex"), 0)
        self.assertEqual(e.get_attr("is_ri"), 0)
        self.assertEqual(e.get_xsize(), 0)
        self.assertEqual(e.get_ysize(), 0)
        self.assertEqual(e.get_zsize(), 0) 

    def test_copy(self):
        """test copy() function ............................."""
        e = EMData()
        e.set_size(32,32,32)
        e.to_zero()
        e.process("testimage.noise.uniform.rand")
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
        e.process("testimage.noise.uniform.rand")
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
        e.process("testimage.noise.uniform.rand")
        
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
        e.process("testimage.noise.uniform.rand")
        
        e2 = EMData()
        e2.set_size(10,10,10)
        e2.to_zero()
        e2.process("testimage.noise.uniform.rand")
        
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
        e.process("testimage.noise.uniform.rand")
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
        e.process("testimage.noise.uniform.rand")
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
            
    def test_window_padded(self):
        """test window_padded() function ...................."""
        e = EMData()
        e.set_size(64,64,64)
        e.process("testimage.noise.uniform.rand")
        e2 = e.window_padded(32)
        
        #window_padded() only work for real data, raise exception for complex
        e.set_complex(True)
        self.assertRaises( RuntimeError, e.window_padded, 32 )
        try:
            e2 = e.window_padded(32)
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "ImageFormatException")
            
        #window_padded() apply to cubic real space image only
        e3 = EMData()
        e3.set_size(64,64,62)
        e3.process("testimage.noise.uniform.rand")
        self.assertRaises( RuntimeError, e3.window_padded, 32 )
        try:
            e4 = e3.window_padded(32)
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "ImageFormatException")
            
    def test_center_origin(self):
        """test center_origin() function ...................."""
        e = EMData()
        e.set_size(32,32,32)
        e.process("testimage.noise.uniform.rand")
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
        e.process("testimage.noise.uniform.rand")
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
        e.process("testimage.noise.uniform.rand")
        
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
        e.process("testimage.noise.uniform.rand")
        
        #test default argument
        e2 = e.pad_fft()
        
        #test arbitrary argument
        e3 = e.pad_fft(3)
        
    def test_postift_depad_corner_inplace(self):
        """test postift_depad_corner_inplace() .............."""
        e = EMData()
        e.set_size(32,32,32)
        e.process("testimage.noise.uniform.rand")
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
        e.process("testimage.noise.uniform.rand")
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
        #os.unlink(outfile)
        
    def test_FH2F(self):
        """test FH2F() function ............................."""
        e = EMData()
        e.set_size(31,31,1)
        e.process("testimage.noise.uniform.rand")
        e2 = e.real2FH(1.0)
        e3 = e2.FH2F(31, 1.0)
        
        #for image not FH, should raise exception
        self.assertRaises( RuntimeError, e.FH2F, 31, 1.0)
        try:
            e.FH2F(31, 1.0)
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "ImageFormatException")

    def test_do_fft(self):
        """test do_fft()/do_ift() function .................."""
        e = EMData()
        e.set_size(32,32,32)
        e.process("testimage.noise.uniform.rand")
        e2 = e.do_fft()
        e3 = e2.do_ift()
        
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
            
    def test_do_fft_inplace(self):
        """test do_fft_inplace()/do_ift_inplace ............."""
        e = EMData()
        e.set_size(32,32,32)
        e.process("testimage.noise.uniform.rand")
        
        #test unpadded real image
        e2 = e.do_fft_inplace()    
        
        #test padded real image
        e = EMData()
        e.set_size(32,32,32)
        e.process("testimage.noise.uniform.rand")
        e2 = e.pad_fft()
        e3 = e2.do_fft_inplace()
        
        e.set_complex(True)
        e.do_ift_inplace()
        
        #do_ift_inplace() only apply to complex image
        e4 = EMData()
        e4.set_size(32,32,32)
        e4.process("testimage.noise.uniform.rand")
        self.assertRaises( RuntimeError, e4.do_ift_inplace, )
        try:
            e4.do_ift_inplace()
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "ImageFormatException")
        
        #do_fft_inplace() only apply to real image
        e5 = EMData()
        e5.set_size(32,32,32)
        e5.process("testimage.noise.uniform.rand")
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
        e.process("testimage.noise.uniform.rand")
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
        e.process("testimage.noise.uniform.rand")
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
        e.process("testimage.noise.uniform.rand")
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
        e.process("testimage.noise.uniform.rand")
        str = e.render_amp8(0, 0, 32, 32, 96, 1.2, 1, 254, 100, 200, 3)
        
        #only apply to 2D image
        e2 = EMData()
        e2.set_size(32,32,32)
        e2.process("testimage.noise.uniform.rand")
        self.assertRaises( RuntimeError, e2.render_amp8, 0, 0, 32, 32, 96, 1.2, 1, 254, 100, 200, 3)
        try:
            str = e2.render_amp8(0, 0, 32, 32, 96, 1.2, 1, 254, 100, 200, 3)
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "ImageDimensionException")
            
    def test_ri2ap_ap2ri(self):
        """test ri2ap()/ap2ri() function ...................."""
        e = EMData()
        e.set_size(32,32,32)
        e.process("testimage.noise.uniform.rand")
        e2 = e.do_fft()
        e2.ri2ap()
        e2.ap2ri()
         
    def test_scale(self):
        """test scale() function ............................"""
        e = EMData()
        e.set_size(32,32,32)
        e.process("testimage.noise.uniform.rand")
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
        e.process("testimage.noise.uniform.rand")
        
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
        """test rotate_2d ..................................."""
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
        e.process("testimage.noise.uniform.rand")
        e2 = EMData()
        e2.set_size(32,32,1)
        e2.process("testimage.noise.uniform.rand")
        e.dot_rotate_translate(e2, 2.0, 3.0, 1.0)
        
        #two image must be the same size
        e3 =EMData()
        e3.set_size(32,32,1)
        e3.process("testimage.noise.uniform.rand")
        e4 = EMData()
        e4.set_size(24,24,1)
        e4.process("testimage.noise.uniform.rand")
        self.assertRaises( RuntimeError, e3.dot_rotate_translate, e4, 2.0, 3.0, 1.0)
        try:
            e3.dot_rotate_translate(e4, 2.0, 3.0, 1.0)
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "ImageFormatException")

        #two image must be 2D
        e5 =EMData()
        e5.set_size(32,32,32)
        e5.process("testimage.noise.uniform.rand")
        e6 = EMData()
        e6.set_size(32,32,32)
        e6.process("testimage.noise.uniform.rand")
        self.assertRaises( RuntimeError, e5.dot_rotate_translate, e6, 2.0, 3.0, 1.0)
        try:
            e5.dot_rotate_translate(e6, 2.0, 3.0, 1.0)
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "ImageDimensionException")
        
    def test_little_big_dot(self):
        """test little_big_dot() function ..................."""
        big = EMData()
        big.set_size(256,256,1)
        big.process("testimage.noise.uniform.rand")
        small = EMData()
        small.set_size(32,32,1)
        small.process("testimage.noise.uniform.rand")
        
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
        e.process("testimage.noise.uniform.rand")
        
        e2 = e.do_radon()
        
        #this function only apply to square 2D image
        e3 = EMData()
        e3.set_size(32,24,1)
        e3.process("testimage.noise.uniform.rand")
        self.assertRaises( RuntimeError, e3.do_radon,)
        try:
            e3.do_radon()
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "ImageFormatException")
        
        e4 = EMData()
        e4.set_size(32,32,32)
        e4.process("testimage.noise.uniform.rand")
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
        e.process("testimage.noise.uniform.rand")
        e2 = EMData()
        e2.set_size(24,1,1)
        e2.process("testimage.noise.uniform.rand")
        re = e.calc_ccf(e2)
        self.assertEqual(re.is_complex(), False)
        
        #for two 2 D images
        e3 = EMData()
        e3.set_size(24,24,1)
        e3.process("testimage.noise.uniform.rand")
        e4 = EMData()
        e4.set_size(24,24,1)
        e4.process("testimage.noise.uniform.rand")
        re = e3.calc_ccf(e4)
        self.assertEqual(re.is_complex(), False)
        
        #for two 3 D images
        e5 = EMData()
        e5.set_size(24,24,24)
        e5.process("testimage.noise.uniform.rand")
        e6 = EMData()
        e6.set_size(24,24,24)
        e6.process("testimage.noise.uniform.rand")
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
        
    def test_calc_ccfx(self):
        """test calc_ccfx() function ........................"""
        e = EMData()
        e.set_size(24,24,1)
        e.process("testimage.noise.uniform.rand")
        e2 = EMData()
        e2.set_size(24,24,1)
        e2.process("testimage.noise.uniform.rand")
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
                    self.assertEqual(narray[i][j][k], 1)
                    
    def test_calc_radial_dist(self):
        """test calc_radial_dist() .........................."""
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
        
    def no_test_statistics(self):
        """test statistics of EMData ........................"""
        e = EMData()
        e.set_size(10,10,1)
        e.process("testimage.circlesphere",{"radius":4})
        descriptive_statistics(e)
        f = norm_pad_ft(e, False, True, 2)
        descriptive_statistics(f)
        f.do_ift_inplace()
        descriptive_statistics(f)
        g = f*10
        descriptive_statistic   
             
def test_main():
    test_support.run_unittest(TestEMData)

if __name__ == '__main__':
    test_main()



