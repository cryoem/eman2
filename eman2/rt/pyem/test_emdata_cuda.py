#!/usr/bin/env python

#
# Author: David Woolford March 2nd 2009
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

import EMAN2
from EMAN2 import *
from pyemtbx.exceptions import *
import unittest
import testlib
import sys
import math
import os
from optparse import OptionParser

IS_TEST_EXCEPTION = False

class TestEMDataCuda(unittest.TestCase):
	"""this is the unit test that verifies the CUDA functionality of the EMData class"""
	
	def test_cuda_phase_comparator(self):
		"""test cuda phase comparator ......................."""
		for z in [1,15,16]:
			for y in [15,16]:
				for x in [15,16]:
					e = EMData()
					e.set_size(x,y,1)
					e.process_inplace("testimage.noise.uniform.rand")
					f = e.copy()
					c = e.cmp("phase",f)
					self.assertAlmostEqual(c,0,4)
	
	def test_cuda_calc_ccfx(self):
		"""test cuda calc_ccfx .............................."""
		for y in [15,16]:
			for x in [15,16]:
				e = EMData()
				e.set_size(x,y,1)
				e.process_inplace("testimage.noise.uniform.rand")
				f = e.copy()
				c = e.calc_ccfx(e,0,-1,True)
				f.set_gpu_rw_current()
				d = f.calc_ccfx(e,0,-1,True)
				self.__test_image_almost_equality(c,d,3)
	
	def test_cuda_get_clip(self):
		"""test cuda get_clip ..............................."""
		for z in [1,15,16]:
			for y in [15,16]:
				for x in [15,16]:
					e = EMData()
			        e.set_size(x,y,z)
			        e.process_inplace("testimage.noise.uniform.rand")
        
			        # Just lose the edge pixels
			        if z == 1:
			        	clip_region_1 = Region(1,1,x-2,y-2)
			        else:
			        	clip_region_1 = Region(1,1,1,x-2,y-2,z-2)
			        
			        r = e.get_clip(clip_region_1)
			        e.set_gpu_rw_current()
			        r2 = e.get_clip(clip_region_1)
			        self.__test_image_equality(r,r2)
			        
			        e = EMData()
			        e.set_size(x,y,z)
			        e.process_inplace("testimage.noise.uniform.rand")
			        if z == 1:
			        	clip_region_1 = Region(-1,-1,x+2,y+2)
			        else:
			        	clip_region_1 = Region(-1,-1,-1,x+2,y+2,z+2)
			        r = e.get_clip(clip_region_1)
			        e.set_gpu_rw_current()
			        r2 = e.get_clip(clip_region_1)
			        self.__test_image_equality(r,r2)
			        
	def __test_image_almost_equality(self,a,b,precision):
		for k in range(a.get_zsize()):
			for j in range(a.get_ysize()):
				for i in range(a.get_xsize()):
					self.assertAlmostEqual(a.get_value_at(i,j,k), b.get_value_at(i,j,k), precision)
	
	def __test_image_equality(self,a,b):
		for k in range(a.get_zsize()):
			for j in range(a.get_ysize()):
				for i in range(a.get_xsize()):
					self.assertEqual(a.get_value_at(i,j,k), b.get_value_at(i,j,k))
	
	def test_cuda_unwrap(self):
		"""test cuda unwrap ................................."""
		for x in [32,33]:
			a = test_image(0,size=(x,x))
			b = a.unwrap()
			a.set_gpu_rw_current()
			c = a.unwrap()
			
			self.__test_image_almost_equality(b,c,1)
		
	def test_cuda_ft_fidelity(self):
		"""test cuda fft/ift equals input ..................."""
		
		for z in [1,15,16,27]:
			for y in [15,16,27]:
				for x in [15,16,27]:
					#print x,y,z
					a = EMData(x,y,z)
					a.process_inplace('testimage.noise.uniform.rand')
					b = a.do_fft_cuda()
					c = b.do_ift_cuda()
					
					attrs = ["get_xsize","get_ysize","get_zsize"]
					for attr in attrs:
						self.assertEqual(getattr(c,attr)(),getattr(a,attr)())
					
					self.__test_image_almost_equality(a,c,3)
	def test_cuda_ccf(self):
		"""test cuda ccf equals cpu ccf ....................."""
		for y in [15,16]:
			for x in [15,16]:
				a = test_image(0,size=(x,y))
				b = a.calc_ccf(a)
				# Not textured
				c = a.calc_ccf_cuda(a,False,False)
				self.__test_image_almost_equality(b,c,1)
				
				# Textured
				c = a.calc_ccf_cuda(a,True,False)
				self.__test_image_almost_equality(b,c,1)
		for z in [15,16]:	
			for y in [15,16]:
				for x in [15,16]:
					#print x,y,z
					a = test_image_3d(0,size=(x,y,z))
					b = a.calc_ccf(a)
					# Not textured
					c = a.calc_ccf_cuda(a,False,False)
					self.__test_image_almost_equality(b,c,1)
										
					#Textured there is a problem with textured at the moment, something to do with pitched pointers on the CUDA side
#					c = a.calc_ccf_cuda(a,True)
#					for k in range(c.get_zsize()):
#						for j in range(c.get_ysize()):
#							for i in range(c.get_xsize()):
#								self.assertAlmostEqual(c.get_value_at(i,j,k), b.get_value_at(i,j,k), 1)
							
		
					
	def test_cuda_2d_fft(self):
		"""test cuda 2D fft equals cpu fft .................."""
		for y in [15,16,27]:
			for x in [15,16,27]:
				a = test_image(0,size=(x,y))
				b = a.do_fft()
				c = a.do_fft_cuda()
				self.__test_image_almost_equality(b,c,3)						
	def test_cuda_3d_fft(self):
		"""test cuda 3D fft equals cpu fft .................."""
		for z in [15,16,27]:
			for y in [15,16,27]:
				for x in [15,16,27]:
					a = test_image_3d(0,size=(x,x,x))
					b = a.do_fft()
					c = a.do_fft_cuda()
					self.__test_image_almost_equality(b,c,3)
								
	def test_cuda_basic_mult(self):
		"""test cuda basic multiplication ..................."""
		self.__test_cuda_basic_operation("mult",2.0)
	
	def test_cuda_basic_division(self):
		"""test cuda basic division ........................."""
		self.__test_cuda_basic_operation("div",2.0)
	
	def test_cuda_basic_addition(self):
		"""test cuda basic addition ........................."""
		self.__test_cuda_basic_operation("add",2.0)
		
	def test_cuda_basic_subtraction(self):
		"""test cuda basic subtraction ......................"""
		self.__test_cuda_basic_operation("sub",2.0)
		
	def __test_cuda_basic_operation(self,op,arg):
		"""test cuda basic multiplication ..................."""
		test_suite = [test_image(1,size=(32,32)), test_image(1,size=(33,33)), test_image(1,size=(32,33)),test_image(1,size=(33,32))]
		test_suite.extend([test_image_3d(0,size=(32,32,32)), test_image_3d(0,size=(33,33,33))])
		for a in test_suite:
			#a = EMData(x,x)
			#a.process_inplace('testimage.noise.uniform.rand')
			b = a.copy()
			b.set_gpu_rw_current()
			getattr(b,op)(arg)
			getattr(a,op)(arg)
			self.__test_image_almost_equality(a,b,8)
					
	def test_cuda_standard_projector(self):
		"""test cuda basic projection ......................."""
		#print ""
		#print "The problem with projection is to do with the CPU version, not the GPU version"
		for x in [15,16]:
			a = EMData(x,x,x)
			a.process_inplace('testimage.noise.uniform.rand')
			b = a.project("cuda_standard",Transform())
			c = a.project("standard",Transform())
			self.__test_image_almost_equality(b,c,3)
			
	def test_dt_cpu_gpurw_cpu(self):
		"""test data transfer cpu->gpurw->cpu................"""
		test_suite = [test_image(1,size=(32,32)), test_image(1,size=(33,33)), test_image(1,size=(32,33)), test_image(1,size=(33,32)),test_image(1,size=(600,800))]
		test_suite.extend([test_image_3d(0,size=(32,32,32)), test_image_3d(0,size=(33,33,33))])
		for a in test_suite:
			b = a.copy()
			a._copy_cpu_to_gpu_rw()
			a._copy_gpu_rw_to_cpu()
			self.assertEqual(a==b,True)
	def test_dt_cpu_gpuro_cpu(self):
		"""test data transfer cpu->gpuro->cpu................"""
		test_suite = [test_image(1,size=(32,32)), test_image(1,size=(33,33)), test_image(1,size=(32,33)), test_image(1,size=(33,32)),test_image(1,size=(600,800))]
		test_suite.extend([test_image_3d(0,size=(32,32,32)), test_image_3d(0,size=(33,33,33))])
		for a in test_suite:
			b = a.copy()
			a._copy_cpu_to_gpu_ro()
			a._copy_gpu_ro_to_cpu()
			self.assertEqual(a==b,True)
		
	def test_dt_cpu_gpurw_gpuro_gpurw_cpu(self):
		"""test data transfer cpu->gpurw->gpuro->gpurw->cpu.."""
		# There is a problem here testing odd sized images. This is because gmem requires pitched pointers to deal with odd dimensions (or so I think)
		test_suite = [test_image(1,size=(32,32)),test_image(1,size=(48,48))]
		test_suite.extend([test_image_3d(0,size=(32,32,32)),test_image_3d(0,size=(48,48,48))])
		for a in test_suite:
			b = a.copy()
			b._copy_cpu_to_gpu_rw()
			b._copy_gpu_rw_to_gpu_ro()
			b._copy_gpu_ro_to_gpu_rw()
			b._copy_gpu_rw_to_cpu()
			self.assertEqual(a==b,True)
	
	def test_dt_cpu_gpuro_gpurw_cpu(self):
		# There is a problem here testing odd sized images. This is because gmem requires pitched pointers to deal with odd dimensions (or so I think)
		"""test data transfer cpu->gpuro->gpurw->cpu ........"""
		test_suite = [test_image(1,size=(32,32)),test_image(1,size=(48,48))]
		test_suite.extend([test_image_3d(0,size=(32,32,32)),test_image_3d(0,size=(48,48,48))])
		for a in test_suite:
			b = a.copy()
			b._copy_cpu_to_gpu_ro()
			b._copy_gpu_ro_to_gpu_rw()
			b._copy_gpu_rw_to_cpu()
			self.assertEqual(a==b,True)
	
def test_main():
	p = OptionParser()
	p.add_option('--t', action='store_true', help='test exception', default=False )
	global IS_TEST_EXCEPTION
	opt, args = p.parse_args()
	if opt.t:
		IS_TEST_EXCEPTION = True
	Log.logger().set_level(-1)  #perfect solution for quenching the Log error information, thank Liwei
	
	
	suite = unittest.TestLoader().loadTestsFromTestCase(TestEMDataCuda)
	unittest.TextTestRunner(verbosity=2).run(suite)
	
	testlib.safe_unlink('mydb2')
	testlib.safe_unlink('mydb1')
	testlib.safe_unlink('mydb')

if __name__ == '__main__':
	if EMUtil.cuda_available(): test_main()
