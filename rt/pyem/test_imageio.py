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
import os
import sys
import testlib
import os
from pyemtbx.exceptions import *
from optparse import OptionParser
from EMAN2 import remove_file

IS_TEST_EXCEPTION = False


class ImageIOTester(unittest.TestCase):
	from EMAN2 import get_supported_3d_formats,get_supported_2d_stack_formats,get_supported_3d_stack_formats
	threed_formats = get_supported_3d_formats()
	twod_stack_formats = get_supported_2d_stack_formats()
	threed_stack_formats = get_supported_3d_stack_formats()
	def do_test_read_write(self,ext_type):
		"""
		test write-read of an image type. Automatically tests read/write fidelity
		in 2D, and also 3D, 2D stack, and 3D stack formats (where applicable)
		"""
		filename = 'testimage.'+ext_type
		
		# stage 1 - 2D and 3D (if applicable) testing
		images = []
		images.append(EMData(32,32,1))
		if ext_type in ImageIOTester.threed_formats:
			images.append(EMData(8,8,8))
		for e in images: # vtk works for 2D and 3D
			e.process_inplace("testimage.noise.uniform.rand")
			e.write_image(filename,0)
			f = EMData()
			if(e.get_zsize()>1):
				f.read_image(filename,0,False,None,True)
			else:
				f.read_image(filename,0)
			try:
				self.assertEqual(e==f,True)
			finally:
				remove_file(filename)
		
		images = []
		
		# stage 2 - 2D and 3D stack testing (if applicable)
		if ext_type in ImageIOTester.twod_stack_formats:
			images.append(EMData(32,32,1))
		if ext_type in ImageIOTester.threed_stack_formats:
			images.append(EMData(8,8,8))
			
		for e in images: # vtk works for 2D and 3D
			try:
				for i in [0,1]:
					e.process_inplace("testimage.noise.uniform.rand")
					e.write_image(filename,i)
					f = EMData()
					f.read_image(filename,i)
					self.assertEqual(e==f,True)
			finally:
				remove_file(filename)	
				
class TestEMIO(ImageIOTester):
	"""EM file IO test"""
	
	def test_negative_image_index_em(self):
		"""test ignore negative index for em ................"""
		filename = 'testimage.em'
		e = EMData()
		e.set_size(32,32)
		e.to_zero()
		try:
			e.write_image(filename, -1)
		except RuntimeError, runtime_err:
			err_type = exception_type(runtime_err)
			self.assertEqual(err_type, "ImageWriteException")
		finally:
			os.unlink(filename)
		
	def test_read_write_em(self):
		"""test write-read em  .............................."""
		self.do_test_read_write("em")
		
	def test_emio_write(self):
		"""test em image file write ........................."""
		infile = "test_emio_write_in.mrc"
		TestUtil.make_image_file(infile, IMAGE_MRC)
		e = EMData()
		e.read_image(infile)
		
		outfile = "test_emio_write_out.em"
		e.write_image(outfile, 0, IMAGE_EM)
		TestUtil.check_image(outfile)
		
		os.unlink(infile)
		os.unlink(outfile)
		testlib.unlink_data_header_files(outfile)

class TestPifIO(ImageIOTester):
	"""PIF file IO test"""
	def no_test_read_write_pif(self):
		"""test write-read pif .............................."""
		self.do_test_read_write("pif")

class TestFitsIO(ImageIOTester):
	"""fits file IO test"""
	''' Cant write to FITs '''
	def no_test_read_write_fits(self):
		"""test write-read fits .............................."""
		self.do_test_read_write("fits")

class TestIcosIO(ImageIOTester):
	"""ICOS file IO test"""
	
	def test_write_icos(self):
		"""test ignore negative index for icos .............."""
		filename = 'testimage.icos'
		e = EMData()
		e.set_size(32,32)
		e.to_zero()
		try: 
			e.write_image(filename, -1)
		except RuntimeError, runtime_err:
			err_type = exception_type(runtime_err)
			self.assertEqual(err_type, "ImageWriteException")
		finally:
			os.unlink(filename)
		
	def test_read_write_icos(self):
		"""test write-read icos.............................."""
		self.do_test_read_write("icos")

class TestPNGIO(ImageIOTester):
	"""PNG file IO test"""
	
	def test_write_png(self):
		"""test ignore negative index for png ..............."""
		filename = 'testimage.png'
		e = EMData()
		e.set_size(32,32)
		e.to_zero()
		try:
			e.write_image(filename, -1)	
		except RuntimeError, runtime_err:
			err_type = exception_type(runtime_err)
			self.assertEqual(err_type, "ImageWriteException")
		finally:
			os.unlink(filename)
		
	def no_test_read_write_png(self):
		"""test write-read png .............................."""
		self.do_test_read_write("png")  
		
class TestVTKIO(ImageIOTester):
	"""VTK file IO test"""

	def test_write_vtk(self):
		"""test ignore negative index for vtk ..............."""
		filename = 'testimage.vtk'
		e = EMData()
		e.set_size(32,32)
		e.to_zero()
		try:
			e.write_image(filename, -1)
		except RuntimeError, runtime_err:
			err_type = exception_type(runtime_err)
			self.assertEqual(err_type, "ImageWriteException")	
		finally:
			os.unlink(filename)
	
	def test_read_write_vtk(self):
		"""test write-read vtk .............................."""
		self.do_test_read_write("vtk")

class TestXPLORIO(ImageIOTester):
	"""XPLOR file IO test"""
	
	def test_write_xplor(self):
		"""test ignore negative index for xplor ............."""
		filename = 'testimage.xplor'
		e = EMData()
		e.set_size(32,32)
		e.to_zero()
		try: 
			e.write_image(filename, -1)
		except RuntimeError, runtime_err:
			err_type = exception_type(runtime_err)
			self.assertEqual(err_type, "ImageWriteException")
		finally:
			os.unlink(filename)

	def no_test_read_write_xplor(self):
		"""test write-read xplor  ..........................."""
		self.do_test_read_write("xplor")

class TestPGMIO(unittest.TestCase):
	"""PGM file IO test"""
	
	def test_read_write_pgm(self):
		"""test pgm file read/write ........................."""
		try:
			e = EMData()
			e.set_size(64,64)
			e.process_inplace('testimage.noise.uniform.rand')
			e.write_image('test_image.pgm')
			
			f = EMData()
			f.read_image('test_image.pgm')
			
			#self.assertEqual(f==e,True)
			
		finally:
			testlib.safe_unlink('test_image.pgm')
	
	def no_test_pgm_region_io(self):
		"""test pgm region io ..............................."""
		try:
			e = EMData()
			e.set_size(1024,1024)
			e.process_inplace('testimage.circlesphere', {'radius':300})
			e.write_image('test_circle.pgm')
			
			f = EMData()
			#f.read_image('test_circle.pgm', 0, False, Region(300,300,200,200))
			f.read_image('test_circle.pgm')
		finally:
			testlib.safe_unlink('test_circle.pgm')

class TestSpiderIO(ImageIOTester):
	"""spider file IO test"""

	def test_make_spider(self):
		"""test make spider image file ......................"""
		try:
			file1 = "test_make_spider_1.spi"
			nx1 = 100
			ny1 = 24
			TestUtil.make_image_file(file1, IMAGE_SPIDER, EM_FLOAT, nx1, nx1)
			err = TestUtil.verify_image_file(file1, IMAGE_SPIDER, EM_FLOAT, nx1, nx1)
			self.assertEqual(err, 0)
		finally:
			testlib.safe_unlink(file1)

	def test_overwrite_spider(self):
		"""test overwrite spider image file ................."""
		try:
			file1 = "test_overwrite_spider.spi"
			nx1 = 24
			ny1 = 32
			TestUtil.make_image_file(file1, IMAGE_SINGLE_SPIDER, EM_FLOAT, nx1, ny1)
			TestUtil.make_image_file(file1, IMAGE_SINGLE_SPIDER, EM_FLOAT, nx1*2, ny1*2)
		finally:
			testlib.safe_unlink(file1)  
		
	def test_write_spider_stack(self):
		"""test to write a spider stack image file  ........."""
		try:
			if(os.path.isfile('test.spi')):
				testlib.safe_unlink('test.spi')
			e = EMData() 
			e.set_size(100,100)
			e.process_inplace('testimage.squarecube', {'edge_length':20})
			e.write_image('test.spi', 0)
				
			e.process_inplace('testimage.circlesphere', {'radius':20})
			e.write_image('test.spi', 1)
			
			e.process_inplace('testimage.gaussian', {'sigma':20})
			e.write_image('test.spi', 2)
				
			e.process_inplace('testimage.sinewave', {'wavelength':20})
			e.set_attr('SPIDER.title', 'The fourth image in the stack')
			e.write_image('test.spi', 3)
			
			f = EMData()
			#read the overall herder
			f.read_image('test.spi', -1, True)
			d = f.get_attr_dict()
			img_num = d['SPIDER.maxim']
			
			#read the individual image from a stack
			for i in range(img_num):
				f.read_image('test.spi', i)
				self.assertEqual(f.is_complex(), False)
				self.assertEqual(f.get_xsize(), 100)
				self.assertEqual(f.get_ysize(), 100)
				self.assertEqual(f.get_zsize(), 1)
		finally:
			testlib.safe_unlink('test.spi')
		
	def test_write_transform_spider(self):
		"""test write spi header info from Transform object ."""
		
		filename = 'test_write_transform.'
		img = EMData(32,32)
		img.process_inplace('testimage.noise.uniform.rand')
		t3d = Transform()
		t3d.set_rotation({'type':'spider', 'phi':1.56, 'theta':2.56, 'psi':3.56})
		t3d.set_trans(10.78, 20.78, 30.78)
		t3d.set_scale(4.0)
		img.set_attr('xform.projection', t3d)
		img.write_image(filename+'spi')
		del img
		
		img2 = EMData(filename+'spi')
		self.assertAlmostEqual(img2.get_attr('SPIDER.phi'), 1.56, 3)
		self.assertAlmostEqual(img2.get_attr('SPIDER.theta'), 2.56, 3)
		self.assertAlmostEqual(img2.get_attr('SPIDER.gamma'), 3.56, 3)
		self.assertAlmostEqual(img2.get_attr('SPIDER.dx'), 10.78, 3)
		self.assertAlmostEqual(img2.get_attr('SPIDER.dy'), 20.78, 3)
		self.assertAlmostEqual(img2.get_attr('SPIDER.dz'), 30.78, 3)
		self.assertAlmostEqual(img2.get_attr('SPIDER.scale'), 4.0, 3)
		trans = img2.get_attr('xform.projection')
		d = trans.get_params('spider')
		self.assertAlmostEqual(d['phi'], 1.56, 3)
		self.assertAlmostEqual(d['theta'], 2.56, 3)
		self.assertAlmostEqual(d['psi'], 3.56, 3)
		self.assertAlmostEqual(d['tx'], 10.78, 3)
		self.assertAlmostEqual(d['ty'], 20.78, 3)
		self.assertAlmostEqual(d['tz'], 30.78, 3)
		self.assertAlmostEqual(d['scale'], 4.0, 3)
		del img2
		testlib.safe_unlink(filename+'spi')
		
		filename2 = 'test_write_transform2.'
		img = EMData(32,32,32)
		img.process_inplace('testimage.noise.uniform.rand')
		t3d = Transform()
		t3d.set_rotation({'type':'spider', 'phi':1.56, 'theta':2.56, 'psi':3.56})
		t3d.set_trans(10.78, 20.78, 30.78)
		t3d.set_scale(4.0)
		img.set_attr('xform.align3d', t3d)
		img.write_image(filename2+'spi')
		
		img2 = EMData(filename2+'spi')
		self.assertAlmostEqual(img2.get_attr('SPIDER.phi'), 1.56, 3)
		self.assertAlmostEqual(img2.get_attr('SPIDER.theta'), 2.56, 3)
		self.assertAlmostEqual(img2.get_attr('SPIDER.gamma'), 3.56, 3)
		self.assertAlmostEqual(img2.get_attr('SPIDER.dx'), 10.78, 3)
		self.assertAlmostEqual(img2.get_attr('SPIDER.dy'), 20.78, 3)
		self.assertAlmostEqual(img2.get_attr('SPIDER.dz'), 30.78, 3)
		self.assertAlmostEqual(img2.get_attr('SPIDER.scale'), 4.0, 3)
		trans2 = img2.get_attr('xform.align3d')
		d2 = trans2.get_params('spider')
		self.assertAlmostEqual(d2['phi'], 1.56, 3)
		self.assertAlmostEqual(d2['theta'], 2.56, 3)
		self.assertAlmostEqual(d2['psi'], 3.56, 3)
		self.assertAlmostEqual(d2['tx'], 10.78, 3)
		self.assertAlmostEqual(d2['ty'], 20.78, 3)
		self.assertAlmostEqual(d2['tz'], 30.78, 3)
		self.assertAlmostEqual(d2['scale'], 4.0, 3)
		testlib.safe_unlink(filename2+'spi')
		
	def test_read_write_spi(self):
		"""test write-read spi .............................."""
		self.do_test_read_write("spi")

class TestDF3IO(ImageIOTester):
	"""test read_write DF3 file ............................."""
	imgfile = 'test.df3'
	
	def test_32bit_image(self):
		"""test unsigned int df3 image I/O .................."""
		e = test_image()
		e.write_image(self.imgfile, -1, EMUtil.ImageType.IMAGE_DF3, False, None, EMUtil.EMDataType.EM_UINT)
		f = EMData()
		f.read_image(self.imgfile)
		self.assertEqual(f.get_attr('nx'), 128)
		self.assertEqual(f.get_attr('ny'), 128)
		self.assertEqual(f.get_attr('nz'), 1)
		testlib.safe_unlink(self.imgfile)
		
	def test_16bit_image(self):
		"""test unsigned short df3 image I/O ................"""
		e = test_image()
		e.write_image(self.imgfile, -1, EMUtil.ImageType.IMAGE_DF3, False, None, EMUtil.EMDataType.EM_USHORT)
		f = EMData()
		f.read_image(self.imgfile)
		self.assertEqual(f.get_attr('nx'), 128)
		self.assertEqual(f.get_attr('ny'), 128)
		self.assertEqual(f.get_attr('nz'), 1)
		testlib.safe_unlink(self.imgfile)
		
	def test_8bit_image(self):
		"""test unsigned char df3 image I/O ................."""
		e = test_image()
		e.write_image(self.imgfile, -1, EMUtil.ImageType.IMAGE_DF3, False, None, EMUtil.EMDataType.EM_UCHAR)
		f = EMData()
		f.read_image(self.imgfile)
		self.assertEqual(f.get_attr('nx'), 128)
		self.assertEqual(f.get_attr('ny'), 128)
		self.assertEqual(f.get_attr('nz'), 1)
		testlib.safe_unlink(self.imgfile)

class TestHdfIO(ImageIOTester):
	"""hdf file IO test ....................................."""
	
	def test_make_image(self):
		"""test make hdf image file ........................."""
		imgfile = "test_make_image_1.h5"
		TestUtil.make_image_file(imgfile, IMAGE_HDF, EM_FLOAT)
		err = TestUtil.verify_image_file(imgfile, IMAGE_HDF, EM_FLOAT)
		self.assertEqual(err, 0)
		testlib.safe_unlink(imgfile)
		
	def test_read_image(self):
		"""test read hdf image file ........................."""
		nx = 20
		ny = 30
		nz = 2
		imgfile = "test_read_image.h5"
		TestUtil.make_image_file(imgfile, IMAGE_HDF, EM_FLOAT, nx, ny, nz)

		e = EMData()
		e.read_image(imgfile)
		attrdict = e.get_attr_dict()
		
		self.assertEqual(attrdict["nx"], nx)
		self.assertEqual(attrdict["ny"], ny)
		self.assertEqual(attrdict["nz"], nz)
		self.assertEqual(attrdict["ImageEndian"], "big")
		self.assertEqual(attrdict["datatype"], EM_FLOAT)
		self.assertEqual(attrdict["is_complex"], 0)
		#self.assertEqual(attrdict["maximum"], 325.0)
		#self.assertEqual(attrdict["minimum"], 0.0)

		testlib.safe_unlink(imgfile)  

	def test_hdf_attr_transform(self):
		"""test Transform object as image attibute .........."""
		t = Transform()
		t.to_identity()
		e=EMData(2,2)
		e.set_attr('tr', t)
		testimage = 'testimage.hdf'
		e.write_image(testimage)
		g = EMData()
		g.read_image(testimage)
		tt = g.get_attr('tr')
		for i in range(3):
			for j in range(4):
				if i==j:
					self.assertAlmostEqual(tt.at(i,j), 1.0, 3)
				else:
					self.assertAlmostEqual(tt.at(i,j), 0.0, 3) 
		os.unlink(testimage)
		
	def test_hdf_attr_int_array(self):
		"""test int array as attribute ......................"""
		hdffile = 'testfile.hdf'
		e = EMData(64, 64)
		e.process_inplace('testimage.noise.uniform.rand')
		e.set_attr('int_array1', [1, 2, 3])
		e.set_attr('int_array2', (5, 6, 7))
		e.write_image(hdffile)
		
		e2 = EMData()
		e2.read_image(hdffile)
		d1 = e2.get_attr('int_array1')
		d2 = e2.get_attr('int_array2')
		self.assertEqual(d1, [1,2,3])
		self.assertEqual(d2, [5,6,7])
		os.unlink(hdffile)	  

	def test_hdf_attr_float_array(self):
		"""test float array as attribute ...................."""
		hdffile = 'testfile.hdf'
		e = EMData(64, 64)
		e.process_inplace('testimage.noise.uniform.rand')
		e.set_attr('float_array1', (1.1, 2.2, 3.3))
		e.write_image(hdffile)
		
		e2 = EMData()
		e2.read_image(hdffile)
		d1 = e2.get_attr('float_array1')
		import math
		for i in (1.1,2.2,3.3):
			self.assertAlmostEqual(i, d1[int(math.floor(i)-1)])
		#self.assertAlmostEquals(d1, (1.1,2.2,3.3), 3)
		#self.assertAlmostEquals(d2, (5.5,6.6,7.7), 3)
		os.unlink(hdffile)	

	def test_hdf_attr_boolean(self):
		"""test hdf file boolean attribute .................."""
		hdffile = 'testfile.hdf'
		e = test_image()
		e.set_attr('trueattr', True)
		e.set_attr('falseattr', False)
		e.write_image(hdffile)
		
		e2 = EMData(hdffile)
		self.assertTrue(e2.get_attr('trueattr'))
		self.assertFalse(e2.get_attr('falseattr'))
		os.unlink(hdffile)	
		
	def test_hdf_attr_int(self):
		"""test hdf file integer attribute .................."""
		hdffile = 'testfile.hdf'
		e = test_image()
		e.set_attr('intattr', 12345)
		e.write_image(hdffile)
		
		e2 = EMData(hdffile)
		self.assertEqual(12345, e2.get_attr('intattr'))
		os.unlink(hdffile)	
		
	def test_hdf_attr_float(self):
		"""test hdf file float attribute ...................."""
		hdffile = 'testfile.hdf'
		e = test_image()
		e.set_attr('floatattr', 12.345)
		e.write_image(hdffile)
		
		e2 = EMData(hdffile)
		self.assertAlmostEqual(12.345, e2.get_attr('floatattr'), 3)
		os.unlink(hdffile)	

	def test_hdf_attr_string(self):
		"""test hdf file string attribute ..................."""
		hdffile = 'testfile.hdf'
		e = test_image()
		e.set_attr('stringattr', 'Grant Tang')
		e.write_image(hdffile)
		
		e2 = EMData(hdffile)
		self.assertEqual('Grant Tang', e2.get_attr('stringattr'))
		os.unlink(hdffile)	

	def test_write_image(self):
		"""test write hdf image file ........................"""
		imgfile1 = "test_write_image_1.mrc"
		imgfile2 = "test_write_image_2.mrc"
		imgfile3 = "test_write_image_3.mrc"

		TestUtil.make_image_file(imgfile1, IMAGE_MRC)
		TestUtil.make_image_file(imgfile2, IMAGE_MRC)
		TestUtil.make_image_file(imgfile3, IMAGE_MRC)

		e1 = EMData()
		e1.read_image(imgfile1)

		e2 = EMData()
		e2.read_image(imgfile2)

		e3 = EMData()
		e3.read_image(imgfile3)

		outfile = "test_write_image_out.h5"
		e1.write_image(outfile, 0)
		e2.write_image(outfile, 1)
		e3.write_image(outfile, 2)
		
		os.unlink(imgfile1)
		os.unlink(imgfile2)
		os.unlink(imgfile3)
		os.unlink(outfile)
		
	def test_delete_attribute(self):
		"""test add and delete attribute from hdf file ......"""
		e1 = EMData()
		e1.set_size(32,32)
		e1.process_inplace('testimage.noise.uniform.rand')
		e1.set_attr('Grant', 10000)
		testimage = 'testimage.hdf'
		e1.write_image(testimage)
		del e1
		
		#test whether the attribute 'Grant' has been written to file
		e2 = EMData()
		e2.read_image(testimage)
		self.assertEqual(e2.get_attr('Grant'), 10000)
		
		if(IS_TEST_EXCEPTION):		
			#testwhether the attribute 'Grant' can be removed from file
			e2.del_attr('Grant')
			e2.write_image(testimage)
			del e2
			e3 = EMData()
			e3.read_image(testimage)
			try:
				no_such_attr = e3.get_attr('Grant')
			except RuntimeError, runtime_err:
				err_type = exception_type(runtime_err)
				self.assertEqual(err_type, "NotExistingObjectException")		   
		
		os.unlink(testimage)
		
	def test_read_write_hdf(self):
		"""test write-read hdf .............................."""
		self.do_test_read_write("hdf")

	def test_hdf_16bit_IO(self):
		"""test read/write 16 bit HDF5......................."""
		file = '16bit.hdf'
		e = test_image()
		e.write_image(file, -1, EMUtil.ImageType.IMAGE_HDF, False, None, EMUtil.EMDataType.EM_USHORT)
		e2 = EMData(file)
		self.assertEqual(e2.get_attr('datatype'), EMUtil.EMDataType.EM_USHORT)
		testlib.safe_unlink(file)
	
	def test_hdf_8bit_IO(self):
		"""test read/write 8 bit HDF5........................"""
		file = '8bit.hdf'
		e = test_image()
		e.write_image(file, -1, EMUtil.ImageType.IMAGE_HDF, False, None, EMUtil.EMDataType.EM_UCHAR)
		e2 = EMData(file)
		self.assertEqual(e2.get_attr('datatype'), EMUtil.EMDataType.EM_UCHAR)
		testlib.safe_unlink(file)
		
	def test_hdf_overwrite(self):
		"""test the overwriting of HDF5......................"""
		file = 'overwrite.hdf'
		e = EMData(2,2)
		e.to_zero()
		e.write_image(file)
		e.to_one()
		e.write_image(file)
		f = EMData(file)
		for i in range(2*2):
			self.assertEqual(1, f[i])
		testlib.safe_unlink(file)

class TestMrcIO(ImageIOTester):
	"""mrc file IO test"""
	def test_negative_image_index(self):
		"""test ignore negative image index ................."""
		filename = 'testimage.mrc'
		e = EMData()
		e.set_size(32,32)
		e.to_zero()
		try: 
			e.write_image(filename, -1)
		except RuntimeError, runtime_err:
				err_type = exception_type(runtime_err)
				self.assertEqual(err_type, "ImageWriteException")	
		finally:
			os.unlink(filename)
	
	def test_overwrite(self):
		"""test overwrite mrc image file ...................."""
		base = "overwrite_" + str(os.getpid())
		imgfile1 = base + "_1.mrc"
		TestUtil.make_image_file(imgfile1, IMAGE_MRC, EM_FLOAT, 10, 20, 1)
		TestUtil.make_image_file(imgfile1, IMAGE_MRC, EM_FLOAT, 30, 40, 1)
		e = EMData()
		e.read_image(imgfile1)
		self.assertEqual(TestUtil.verify_image_file(imgfile1, IMAGE_MRC, EM_FLOAT, 30, 40, 1), 0)
		
		imgfile2 = base + "_2.mrc"
		TestUtil.make_image_file(imgfile2, IMAGE_MRC, EM_FLOAT, 30, 40, 1)
		TestUtil.make_image_file(imgfile2, IMAGE_MRC, EM_FLOAT, 10, 20, 1)
		self.assertEqual(TestUtil.verify_image_file(imgfile2, IMAGE_MRC, EM_FLOAT, 10, 20, 1), 0)
		
		imgfile3 = base + "_3.mrc"
		TestUtil.make_image_file(imgfile3, IMAGE_MRC, EM_FLOAT, 30, 40, 1)
		TestUtil.make_image_file2(imgfile3, IMAGE_MRC, EM_FLOAT, 30, 40, 1)
		self.assertEqual(TestUtil.verify_image_file2(imgfile3, IMAGE_MRC, EM_FLOAT, 30, 40, 1), 0)

		os.unlink(imgfile1)
		os.unlink(imgfile2)
		os.unlink(imgfile3)

	def test_make_image_file(self):
		"""test make mrc image file ........................."""
		base = "test_make_image_file"
		img1 = base + "_c.mrc"
		img2 = base + "_s.mrc"
		img3 = base + "_f.mrc"
		img4 = base + "_sc.mrc"
		img5 = base + "_fc.mrc"
		
		TestUtil.make_image_file(img1, IMAGE_MRC, EMUtil.EMDataType.EM_UCHAR)
		self.assertEqual(TestUtil.verify_image_file(img1, IMAGE_MRC, EMUtil.EMDataType.EM_UCHAR), 0)
  
		TestUtil.make_image_file(img2, IMAGE_MRC, EM_USHORT)
		self.assertEqual(TestUtil.verify_image_file(img2, IMAGE_MRC, EM_USHORT), 0)
		
		TestUtil.make_image_file(img3, IMAGE_MRC, EM_FLOAT,64,64)
		self.assertEqual(TestUtil.verify_image_file(img3, IMAGE_MRC, EM_FLOAT, 64,64), 0)
		
		TestUtil.make_image_file(img4, IMAGE_MRC, EM_SHORT_COMPLEX)
		self.assertEqual(TestUtil.verify_image_file(img4, IMAGE_MRC, EM_SHORT_COMPLEX), 0)
		
		TestUtil.make_image_file(img5, IMAGE_MRC, EM_FLOAT_COMPLEX)
		self.assertEqual(TestUtil.verify_image_file(img5, IMAGE_MRC, EM_FLOAT_COMPLEX), 0)

		os.unlink(img1)
		os.unlink(img2)
		os.unlink(img3)
		os.unlink(img4)
		os.unlink(img5)

		img6 = base + "_3d.mrc"
		TestUtil.make_image_file(img6, IMAGE_MRC, EM_FLOAT, 16,16,10)
		os.unlink(img6)

	def test_complex_image(self):
		"""test complex mrc image file ......................"""
		imgfile1 = "test_complex_image.mrc"
		TestUtil.make_image_file(imgfile1, IMAGE_MRC, EM_FLOAT)

		imgfile2 = "test_complex_image_fft.mrc"

		e = EMData()
		e.read_image(imgfile1)
		fft = e.do_fft()
		fft.write_image(imgfile2)

		self.assertEqual(TestUtil.verify_image_file(imgfile2, IMAGE_MRC, EM_FLOAT_COMPLEX,
													e.get_xsize(), e.get_ysize(),
													e.get_zsize()), 0)
		
		os.unlink(imgfile1)
		os.unlink(imgfile2)
   
	def test_mrcio_label(self):
		"""test mrc file label .............................."""
		pid = str(os.getpid())
		infile = "test_mrcio_label_in_" + pid + ".mrc"
		TestUtil.make_image_file(infile, IMAGE_MRC)

		e = EMData()
		e.read_image(infile)
		label = "ByLiweiPeng"
		label_i = 3
		labelname = "MRC.label" + str(label_i)
		
		e.set_attr(labelname, label)

		outfile="test_mrcio_label_out_" + pid + ".mrc"		
		e.write_image(outfile, 0, IMAGE_MRC)
		
		e2 = EMData()
		e2.read_image(outfile)
		d = e2.get_attr_dict()
		nlabels = int(d["MRC.nlabels"])

		os.unlink(outfile)
		os.unlink(infile)
		
		self.assert_(nlabels > label_i)
		self.assertEqual(d[labelname], label)
		
	def no_test_write_transform_mrc(self):
		"""test write mrc header info from Transform object ."""
		filename = 'test_write_transform.'
		img = EMData(32,32)
		img.process_inplace('testimage.noise.uniform.rand')
		t3d = Transform()
		t3d.set_rotation({'type':'imagic', 'alpha':1.56, 'beta':2.56, 'gamma':3.56})
		t3d.set_trans(10.78, 20.78, 30.78)
		img.set_attr('xform.projection', t3d)
		img.write_image(filename+'mrc')
		del img
		
		img2 = EMData(filename+'mrc')
		self.assertAlmostEqual(img2.get_attr('MRC.alpha'), 1.56, 3)
		self.assertAlmostEqual(img2.get_attr('MRC.beta'), 2.56, 3)
		self.assertAlmostEqual(img2.get_attr('MRC.gamma'), 3.56, 3)
		self.assertAlmostEqual(img2.get_attr('origin_x'), 10.78, 3)
		self.assertAlmostEqual(img2.get_attr('origin_y'), 20.78, 3)
		self.assertAlmostEqual(img2.get_attr('origin_z'), 30.78, 3)
		trans = img2.get_attr('xform.projection')
		d = trans.get_params('imagic')
		self.assertAlmostEqual(d['alpha'], 0.0, 3)
		self.assertAlmostEqual(d['beta'], 0.0, 3)
		self.assertAlmostEqual(d['gamma'], 0.0, 3)
		self.assertAlmostEqual(d['MRC.nxstart'], 10.78, 3)
		self.assertAlmostEqual(d['MRC.nxstart'], 20.78, 3)
		self.assertAlmostEqual(d['MRC.nxstart'], 30.78, 3)
		del img2
		testlib.safe_unlink(filename+'mrc')
		
		filename2 = 'test_write_transform2.'
		img = EMData(32,32,32)
		img.process_inplace('testimage.noise.uniform.rand')
		t3d = Transform()
		t3d.set_rotation({'type':'imagic', 'alpha':1.56, 'beta':2.56, 'gamma':3.56})
		t3d.set_trans(10.78, 20.78, 30.78)
		img.set_attr('xform.align3d', t3d)
		img.write_image(filename2+'mrc')
		
		img2 = EMData(filename2+'mrc')
		self.assertAlmostEqual(img2.get_attr('MRC.alpha'), 1.56, 3)
		self.assertAlmostEqual(img2.get_attr('MRC.beta'), 2.56, 3)
		self.assertAlmostEqual(img2.get_attr('MRC.gamma'), 3.56, 3)
		self.assertAlmostEqual(img2.get_attr('origin_x'), 10.78, 3)
		self.assertAlmostEqual(img2.get_attr('origin_y'), 20.78, 3)
		self.assertAlmostEqual(img2.get_attr('origin_z'), 30.78, 3)
		trans2 = img2.get_attr('xform.align3d')
		d2 = trans2.get_params('imagic')
		self.assertAlmostEqual(d2['alpha'], 1.56, 3)
		self.assertAlmostEqual(d2['beta'], 2.56, 3)
		self.assertAlmostEqual(d2['gamma'], 3.56, 3)
		self.assertAlmostEqual(d2['tx'], 10.78, 3)
		self.assertAlmostEqual(d2['ty'], 20.78, 3)
		self.assertAlmostEqual(d2['tz'], 30.78, 3)
		testlib.safe_unlink(filename2+'mrc')

	def test_read_write_mrc(self):
		"""test write-read mrc .............................."""
		self.do_test_read_write("mrc")


class TestImagicIO(ImageIOTester):
	"""imagic file IO test"""
	
	def test_no_ext_filename(self):  
		"""test no extention file name ......................"""	  
		infile = "test_no_ext_filename.mrc"
		TestUtil.make_image_file(infile, IMAGE_MRC)

		outfile = "test_no_ext_filename_out"
		
		img = EMData()
		img.read_image(infile)
		img.write_image(outfile, 0, IMAGE_IMAGIC)

		(hedfile, imgfile) = testlib.get_imagic_filename_pair(outfile)
		os.unlink(infile)
		os.unlink(hedfile)
		os.unlink(imgfile)

	def test_append_image(self):
		"""test image appending ............................."""
		infile = "test_append_image.hed"
		TestUtil.make_image_file(infile, IMAGE_IMAGIC, EM_FLOAT, 16, 16, 10)
			   
		e = EMData()
		all_imgs = e.read_images(infile)

		outfile1 = "test_append_image_out_" + str(os.getpid()) + ".hed"
		
		for img in all_imgs:
			img.append_image(outfile1)
		
		(infilehed, infileimg) = testlib.get_imagic_filename_pair(infile)
		(outfilehed, outfileimg) = testlib.get_imagic_filename_pair(outfile1)

		os.unlink(infilehed)
		os.unlink(infileimg)
		os.unlink(outfilehed)
		os.unlink(outfileimg)

	def test_append_to_newfile(self):
		"""test append image to new file ...................."""
		infile = "test_append_to_newfile_in.mrc"
		outfile = "test_append_to_newfile_in.img"
		
		TestUtil.make_image_file(infile, IMAGE_MRC)
		e = EMData()
		e.read_image(infile)
		e.append_image(outfile, IMAGE_IMAGIC)

		# check e
		(outhed, outimg) = testlib.get_imagic_filename_pair(outfile)

		os.unlink(infile)
		os.unlink(outhed)
		os.unlink(outimg)

	def test_append_to_existing_file(self):
		"""test append image to existing file ..............."""
		img1 = "test_append_to_existing_file_1.hed"
		TestUtil.make_image_file(img1, IMAGE_IMAGIC)
		e = EMData()
		e.read_image(img1, 0, False, None, True)
		e.append_image(img1, IMAGE_IMAGIC)

		# verify here
		(hedfile, imgfile) = testlib.get_imagic_filename_pair(img1)
		os.unlink(hedfile)
		os.unlink(imgfile)
	
	#new Imagic4D does not take a 3D image as a stack of 2D
	def no_test_insert_to_newfile(self):
		"""test insert image to new file ...................."""
		img1 = "test_insert_to_newfile_in.hed"
		TestUtil.make_image_file(img1, IMAGE_IMAGIC)
		e = EMData()
		e.read_image(img1)
		outfile = "test_insert_to_newfile_out.hed"
		nimg = 4
		e.write_image(outfile, nimg-1, IMAGE_IMAGIC)

		nimg2 = EMUtil.get_image_count(outfile)
		self.assertEqual(nimg2, nimg)
		
		(hedfile1, imgfile1) = testlib.get_imagic_filename_pair(img1)
		(hedfile2, imgfile2) = testlib.get_imagic_filename_pair(outfile)
		
		os.unlink(hedfile1)
		os.unlink(imgfile1)
		os.unlink(hedfile2)
		os.unlink(imgfile2)

	def test_insert_beyond_existing_file(self):
		"""test insert image beyond existing file ..........."""
		infile = "insert_beyond_existing_in.hed"
		TestUtil.make_image_file(infile, IMAGE_IMAGIC)
		e = EMData()
		e.read_image(infile, 0, False, None, True)

		nimg1 = EMUtil.get_image_count(infile)
		self.assertEqual(nimg1, 1)

		n2 = 9
		e.write_image(infile, n2, IMAGE_IMAGIC)
		nimg2 = EMUtil.get_image_count(infile)
		self.assertEqual(nimg2, n2+1)

		# todo: verify images
		
		n3 = 14
		e.write_image(infile, n3, IMAGE_IMAGIC)
		nimg3 = EMUtil.get_image_count(infile)
		self.assertEqual(nimg3, n3+1)

		# todo: verify images
	
		(hedfile, imgfile) = testlib.get_imagic_filename_pair(infile)
		os.unlink(hedfile)
		os.unlink(imgfile)

	#New Imagic4D format does not take a 3D image as a stack of 2D image
	def no_test_insert_inside_existing_file(self):
		"""test insert image in existing file ..............."""
		infile = "test_insert_inside_existing_file_1.img"
		TestUtil.make_image_file(infile, IMAGE_IMAGIC, EM_FLOAT, 20, 30, 20)
		
		insertfile = "test_insert_inside_existing_file_2.mrc"
		TestUtil.make_image_file(insertfile, IMAGE_MRC, EM_FLOAT, 20, 30)
		e = EMData()
		e.read_image(insertfile)
		e.write_image(infile, 2, IMAGE_IMAGIC)

		# verify result
		
		(hedfile1, imgfile1) = testlib.get_imagic_filename_pair(infile)
		os.unlink(hedfile1)
		os.unlink(imgfile1)
		os.unlink(insertfile)
	
	def test_write_transform_to_euler(self):
		"""test write Transform as euler angles ............."""
		filename = 'test_write_transform.'
		img = EMData(32,32)
		img.process_inplace('testimage.noise.uniform.rand')
		t3d = Transform()
		t3d.set_rotation({'type':'imagic', 'alpha':1.56, 'beta':2.56, 'gamma':3.56})
		img.set_attr('xform.projection', t3d)
		img.write_image(filename+'img')
		del img
		
		img2 = EMData(filename+'img')
		self.assertAlmostEqual(img2.get_attr('euler_alpha'), 1.56, 3)
		self.assertAlmostEqual(img2.get_attr('euler_beta'), 2.56, 3)
		self.assertAlmostEqual(img2.get_attr('euler_gamma'), 3.56, 3)
		xform1 = img2.get_attr('xform.projection')
		d = xform1.get_rotation('imagic')
		self.assertAlmostEqual(d['alpha'], 1.56, 3)
		self.assertAlmostEqual(d['beta'], 2.56, 3)
		self.assertAlmostEqual(d['gamma'], 3.56, 3)
		del img2
		testlib.safe_unlink(filename+'img')
		testlib.safe_unlink(filename+'hed')
		
		filename2 = 'test_write_transform2.'
		img = EMData(32,32,32)
		img.process_inplace('testimage.noise.uniform.rand')
		t3d = Transform()
		t3d.set_rotation({'type':'imagic', 'alpha':1.56, 'beta':2.56, 'gamma':3.56})
		img.set_attr('xform.align3d', t3d)
		img.write_image(filename2+'img')
		
		img2 = EMData(filename2+'img')
		self.assertAlmostEqual(img2.get_attr('euler_alpha'), 1.56, 3)
		self.assertAlmostEqual(img2.get_attr('euler_beta'), 2.56, 3)
		self.assertAlmostEqual(img2.get_attr('euler_gamma'), 3.56, 3)
		xform2 = img2.get_attr('xform.align3d')
		d2 = xform2.get_rotation('imagic')
		self.assertAlmostEqual(d2['alpha'], 1.56, 3)
		self.assertAlmostEqual(d2['beta'], 2.56, 3)
		self.assertAlmostEqual(d2['gamma'], 3.56, 3)
		testlib.safe_unlink(filename2+'hed')
		testlib.safe_unlink(filename2+'img')
		
	def test_eman1ctf_io(self):
		"""test EMAN1Ctf object I/O ........................."""
		filename = 'test_imagic_ctf.'
		img = EMData(32,32)
		img.process_inplace('testimage.noise.uniform.rand')
		ctf1 = EMAN1Ctf()
		ctf1.from_vector((1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0))
		img.set_attr('ctf', ctf1)
		img.write_image(filename+'img')
		
		img2 = EMData(filename + 'img')
		ctf2 = img2.get_attr('ctf')
		self.assertEqual(ctf1.to_string(), ctf2.to_string())
		del ctf1, ctf2
		testlib.safe_unlink(filename+'hed')
		testlib.safe_unlink(filename+'img')
	
	def test_eman2ctf_io(self):
		"""test EMAN2Ctf object I/O ........................."""
		filename = 'test_imagic_ctf2.hdf'
		img = EMData(32,32)
		img.process_inplace('testimage.noise.uniform.rand')
		ctf1 = EMAN2Ctf()
		ctf1.from_vector((1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0, 0, 0))
		img.set_attr('ctf', ctf1)
		img.write_image(filename)
		
		img2 = EMData(filename)
		ctf2 = img2.get_attr('ctf')
		self.assertEqual(ctf1.to_string(), ctf2.to_string())
		del ctf1, ctf2
		testlib.safe_unlink(filename)
		
	def test_read_write_img(self):
		"""test write-read img .............................."""
		self.do_test_read_write("img")

#	def 

class TestImageIO(unittest.TestCase):
	"""image data IO test"""
		
	def test_emdata_overwriting(self):
		"""test emdata overwriting .........................."""
		imgfile1 = "test_emdata_overwriting_1.mrc"
		imgfile2 = "test_emdata_overwriting_2.mrc"
		TestUtil.make_image_file(imgfile1, IMAGE_MRC, EM_FLOAT, 16,32,2)
		TestUtil.make_image_file(imgfile2, IMAGE_MRC, EM_FLOAT, 32, 24)

		a=EMData()
		b=EMData()
		
		a.read_image(imgfile1, 0)
		b.read_image(imgfile2, 0)
		b.read_image(imgfile1, 0)

		os.unlink(imgfile1)
		os.unlink(imgfile2)

	def test_image_overwriting(self):
		"""test image overwriting ..........................."""
		e = EMData()
		e.set_size(20, 20, 1)
		outfile = "test_image_overwriting_" + str(os.getpid()) + ".mrc"
		e.write_image(outfile)
		e.write_image(outfile)
		os.unlink(outfile)
		
	def create_dummy_region(self,e):
		"""test region creation ............................."""
		nx = e.get_xsize()
		ny = e.get_ysize()
		nz = e.get_zsize()

		x0 = nx/4
		y0 = ny/4
		z0 = nz/4
		
		xsize = nx/2
		ysize = ny/2
		zsize = nz/2
		
		if zsize == 0:
			zsize = 1
			
		ndims = e.get_ndim()

		if ndims == 2:
			region_2d = Region(x0, y0, xsize, ysize)
			region_3d = None
		elif ndims == 3:
			region_2d = Region(x0, y0, z0, xsize, ysize, 1)
			region_3d = Region(x0, y0, z0, xsize, ysize, zsize)

		return (region_2d, region_3d)

	def region_read_test(self, imgtype, imgfile, outtype = None):	
		"""test region read ................................."""
		if not outtype:
			outtype = imgtype
			
		is_3d = False
		if imgtype == IMAGE_IMAGIC:
		   is_3d = True

		imgbase = Util.remove_filename_ext(imgfile)
		ext = Util.get_filename_ext(imgfile)

		readfile_2d = imgbase + "_read_region_2d." + ext
		readfile_3d = imgbase + "_read_region_3d." + ext

		e = EMData()
		e.read_image(imgfile, 0, False, None, is_3d)
	 
		TestUtil.check_image(imgfile, e)
		 
		#testlib.unlink_data_header_files(imgfile)	
		 
		(region_2d, region_3d) = self.create_dummy_region(e)

		e2 = EMData()
		e2.read_image(imgfile, 0, False, region_2d, is_3d)
		 
		TestUtil.check_image(readfile_2d, e2)
		e2.write_image(readfile_2d, 0, outtype)
		 
		if region_3d:
			e4 = EMData()
			e4.read_image(imgfile, 0, False, region_3d, is_3d)
			e4.write_image(readfile_3d, 0, outtype)
			TestUtil.check_image(readfile_3d, e4)

			if outtype == IMAGE_IMAGIC:
				(hed3d, img3d) = testlib.get_imagic_filename_pair(readfile_3d)
				os.unlink(hed3d)
				os.unlink(img3d)
			else:
				os.unlink(readfile_3d)
			 
		if outtype == IMAGE_IMAGIC:
			(hed2d, img2d) = testlib.get_imagic_filename_pair(readfile_2d)
			os.unlink(hed2d)
			os.unlink(img2d)
		else:
			os.unlink(readfile_2d)

	def region_write_test(self, imgtype, imgfile, outtype = None):
		"""test region writing .............................."""
		if not outtype:
			outtype = imgtype
		
		is_3d = False
		if imgtype == IMAGE_IMAGIC:
			is_3d = True

		imgbase = Util.remove_filename_ext(imgfile)
		ext = Util.get_filename_ext(imgfile)
		writefile_2d = imgbase + "_write_region_2d." + ext
		writefile_3d = imgbase + "_write_region_3d." + ext

		e = EMData()
		e.read_image(imgfile, 0, False, None, is_3d)

		TestUtil.check_image(imgfile, e)
		testlib.unlink_data_header_files(imgfile)
		 
		ndims = e.get_ndim()
		e.write_image(writefile_2d, 0, outtype)
		
		if ndims == 3:
			e.write_image(writefile_3d, 0, outtype)

		(region_2d, region_3d) = self.create_dummy_region(e)

		zsize = e.get_zsize()/2
		if zsize == 0:
			zsize = 1
			
		e3 = EMData()
		e3.set_size(e.get_xsize()/2, e.get_ysize()/2, zsize)
		e3.to_zero()

		image_index = 0
		if outtype == IMAGE_SPIDER:
			image_index = e.get_zsize()/2
		
		e3.write_image(writefile_2d, image_index, outtype, False, region_2d)
		TestUtil.check_image(writefile_2d)

		if ndims == 3:
			e3.write_image(writefile_3d, image_index, outtype, False, region_3d)
			TestUtil.check_image(writefile_3d)

			if outtype  == IMAGE_IMAGIC:
				(hed3d, img3d) = testlib.get_imagic_filename_pair(writefile_3d)
				os.unlink(hed3d)
				os.unlink(img3d)
			else:
				os.unlink(writefile_3d)
				
		if outtype == IMAGE_IMAGIC:
			(hed2d, img2d) = testlib.get_imagic_filename_pair(writefile_2d)
			os.unlink(hed2d)
			os.unlink(img2d)
		else:
			os.unlink(writefile_2d)

	def region_read_write_test(self, imgtype, imgfile, outtype = None):
		"""test region read and write ......................."""
		self.region_read_test(imgtype, imgfile, outtype)
		self.region_write_test(imgtype, imgfile, outtype)
		
		
	def test_region_equiv_to_clip(self):
		"""test read region is identical to get clip ........"""
		
		# note support for going beyond image boundaries, which was the main purpose of this function, is tested here,
		# but that getting the entire region of the image as a clip is also tested.
		
		from EMAN2 import get_supported_3d_formats
		
		#test for 2D image	
		size = 8
		e = EMData()
		e.set_size(size,size)
		e.process_inplace("testimage.noise.uniform.rand")
		
		fmts = ["mrc","spi","em","img", "hdf"]
		
		for fmt in fmts:
			name = "testregionimage."+fmt
			e.write_image(name)
			try:
				for i in range(-2,3):
					for j in range(-2,3):
							for l in range(-2,3):
								for m in range(-2,3):
									region = Region(i,j,size+l,size+m)
									f = EMData()
									f.read_image(name,0,False,region)
									g = e.get_clip(region)
									self.assertEqual(f==g,True)
			finally:
				remove_file(name)
		
		
		#test for 3D image
		size = 8
		e = EMData(size,size,size)
		e.process_inplace('testimage.axes')
		
		fmts = get_supported_3d_formats()
		fmts =  ["spi","mrc","icos","em", "hdf"]
		unsupported = ["img","xplor","pif","emim","vtk"] # so many that we can't use :(
		for fmt in fmts:
			name = "testregionimage."+fmt
			e.write_image(name)
			try:
				for i in range(-2,3):
					for j in range(-2,3):
						for k in range(-2,3):
							for l in range(-2,3):
								for m in range(-2,3):
									for n in range(-2,3):
										region = Region(i,j,k,size+l,size+m,size+n)
										f = EMData()
										f.read_image(name,0,False,region)
										g = e.get_clip(region)
										self.assertEqual(f==g,True)
			finally:
				remove_file(name)
		
	def test_mrcio_region(self):
		"""test mrc io region ..............................."""
		mrc2d = "test_mrcio_region_2d.mrc"
		mrc3d = "test_mrcio_region_3d.mrc"
		TestUtil.make_image_file(mrc2d, IMAGE_MRC, EM_FLOAT, 32,64)
		TestUtil.make_image_file(mrc3d, IMAGE_MRC, EM_FLOAT, 32,32,12)		
		self.region_read_write_test(IMAGE_MRC, mrc2d)
		self.region_read_write_test(IMAGE_MRC, mrc3d)
		os.unlink(mrc2d)
		os.unlink(mrc3d)

	#regional I/O should be work for Imagic4D, need to test this back later.
	def no_test_imagicio_region(self):
		"""test imagic io region ............................"""
		infile = "test_imagicio_region_11.hed"
		TestUtil.make_image_file(infile, IMAGE_IMAGIC, EM_FLOAT, 32,32,64)
		self.region_read_write_test(IMAGE_IMAGIC, infile)

		(hedfile1, imgfile1) = testlib.get_imagic_filename_pair(infile)
		os.unlink(hedfile1)
		os.unlink(imgfile1)

	def test_hdfio_region(self):
		"""test hdf io region ..............................."""
#		file1 = "test_hdfio_region_1.h5"
#		nx = 48
#		ny = 64
#		nz1 = 1
#		TestUtil.make_image_file(file1, IMAGE_HDF, EM_FLOAT, nx, ny, nz1)
#		self.region_read_write_test(IMAGE_HDF, file1)
#
#		file2 = "test_hdfio_region_2.h5"
#		nz2 = 12
#		TestUtil.make_image_file(file2, IMAGE_HDF, EM_FLOAT, nx, ny, nz2)
#		self.region_read_write_test(IMAGE_HDF, file2)
#
#		os.unlink(file1)
#		os.unlink(file2)
		
		file3 = "test2D.hdf"
		img = EMData(4,4)
		for i in range(4):
			for j in range(4):
				img.set_value_at(j, i, i*4+j)
		img.write_image(file3)
		
		f = EMData()
		ff = EMData(file3)
		fclip = ff.get_clip(Region(0,0,2,2))
		f.read_image(file3, 0, False, Region(0,0,2,2))
		self.assertEqual(f.get_value_at(0,0), fclip.get_value_at(0,0))
		self.assertEqual(f.get_value_at(0,1), fclip.get_value_at(0,1))
		self.assertEqual(f.get_value_at(1,0), fclip.get_value_at(1,0))
		self.assertEqual(f.get_value_at(1,1), fclip.get_value_at(1,1))
		
		g = EMData()
		g.read_image(file3, 0, False, Region(2,2,2,2))
		gclip = ff.get_clip(Region(2,2,2,2))
		self.assertEqual(g.get_value_at(0,0), gclip.get_value_at(0,0))
		self.assertEqual(g.get_value_at(0,1), gclip.get_value_at(0,1))
		self.assertEqual(g.get_value_at(1,0), gclip.get_value_at(1,0))
		self.assertEqual(g.get_value_at(1,1), gclip.get_value_at(1,1))

		os.unlink(file3)
		
		file4 = "test3D.hdf"
		img2 = EMData(4,4,4)
		for i in range(4):
			for j in range(4):
				for k in range(4):
					img2.set_value_at(k, j, i, i*4*4+j*4+k)
		img2.write_image(file4)
		
		f2 = EMData()
		ff2 = EMData(file4)
		f2.read_image(file4, 0, False, Region(0,0,0,2,2,2))
		f2clip = ff2.get_clip(Region(0,0,0,2,2,2))
		self.assertEqual(f2.get_value_at(0,0,0), f2clip.get_value_at(0,0,0))
		self.assertEqual(f2.get_value_at(0,1,0), f2clip.get_value_at(0,1,0))
		self.assertEqual(f2.get_value_at(1,0,0), f2clip.get_value_at(1,0,0))
		self.assertEqual(f2.get_value_at(1,1,0), f2clip.get_value_at(1,1,0))
		self.assertEqual(f2.get_value_at(0,0,1), f2clip.get_value_at(0,0,1))
		self.assertEqual(f2.get_value_at(0,1,1), f2clip.get_value_at(0,1,1))
		self.assertEqual(f2.get_value_at(1,0,1), f2clip.get_value_at(1,0,1))
		self.assertEqual(f2.get_value_at(1,1,1), f2clip.get_value_at(1,1,1))
		
		os.unlink(file4)

"""
	def  test_spiderio_region(self):
		file1 = "test_spiderio_region_1.h5"
		
  
	def test_spiderio_region(self):
		self.region_read_write_test(IMAGE_SINGLE_SPIDER, "spider-single.spi")
		self.region_read_write_test(IMAGE_SPIDER, "spider-stack.spi")

	def test_emio_region(self):
		self.region_read_write_test(EM, "20s2d.em")
		self.region_read_write_test(EM, "stack3d.em")

	def test_pifio_region(self):
		self.region_read_write_test(PIF, "sv-3d.pif")

	def test_xplorio_region(self):
		self.region_read_write_test(XPLOR, "2f.xplor")
"""

def test_main():
	platform = get_platform()

	p = OptionParser()
	p.add_option('--t', action='store_true', help='test exception', default=False )
	global IS_TEST_EXCEPTION
	opt, args = p.parse_args()
	if opt.t:
		IS_TEST_EXCEPTION = True
	Log.logger().set_level(-1)  #perfect solution for quenching the Log error information, thank Liwei
	suite1 = unittest.TestLoader().loadTestsFromTestCase(TestEMIO)
	unittest.TextTestRunner(verbosity=2).run(suite1)
	
	suite2 = unittest.TestLoader().loadTestsFromTestCase(TestIcosIO)
	unittest.TextTestRunner(verbosity=2).run(suite2)
	
	if platform!='Windows' and platform!='win32':
		suite3 = unittest.TestLoader().loadTestsFromTestCase(TestPNGIO)
		unittest.TextTestRunner(verbosity=2).run(suite3)
	
	suite4 = unittest.TestLoader().loadTestsFromTestCase(TestVTKIO)
	unittest.TextTestRunner(verbosity=2).run(suite4)
	
	suite5 = unittest.TestLoader().loadTestsFromTestCase(TestXPLORIO)
	unittest.TextTestRunner(verbosity=2).run(suite5)
	
	suite6 = unittest.TestLoader().loadTestsFromTestCase(TestPGMIO)
	unittest.TextTestRunner(verbosity=2).run(suite6)
	
	suite7 = unittest.TestLoader().loadTestsFromTestCase(TestSpiderIO)
	unittest.TextTestRunner(verbosity=2).run(suite7)
	
	suite8 = unittest.TestLoader().loadTestsFromTestCase(TestImageIO)
	unittest.TextTestRunner(verbosity=2).run(suite8)
	
	suite9 = unittest.TestLoader().loadTestsFromTestCase(TestHdfIO)
	unittest.TextTestRunner(verbosity=2).run(suite9)
	
	suite10 = unittest.TestLoader().loadTestsFromTestCase(TestMrcIO)
	unittest.TextTestRunner(verbosity=2).run(suite10)
	
	suite11 = unittest.TestLoader().loadTestsFromTestCase(TestImagicIO)
	unittest.TextTestRunner(verbosity=2).run(suite11)
	
	suite12 = unittest.TestLoader().loadTestsFromTestCase(TestPifIO)
	unittest.TextTestRunner(verbosity=2).run(suite12)
	
	#suite13 = unittest.TestLoader().loadTestsFromTestCase(TestFitsIO)
	#inittest.TextTestRunner(verbosity=2).run(suite13)
	
	suite14 = unittest.TestLoader().loadTestsFromTestCase(TestDF3IO)	
	unittest.TextTestRunner(verbosity=2).run(suite14) 
	
if __name__ == '__main__':
	test_main()
