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
	
	# Known issues
	# cuda irregular sized FFTs are not invertible (test_image(0,size=(15,16)))
	# cuda driver has bugs (March 2009) relating to non power of two sized data
	
	def test_cuda_ft_fidelity(self):
		"""test cuda fft/ift equals input ..........................."""
		
		for z in [1]:
			for y in [15,16]:
				for x in [15,16]:
					#print x,y,z
					a = EMData(x,y,z)
					a.process_inplace('testimage.noise.uniform.rand')
					b = a.do_fft_cuda()
					c = b.do_ift_cuda()
					
					attrs = ["get_xsize","get_ysize","get_zsize"]
					for attr in attrs:
						self.assertEqual(getattr(c,attr)(),getattr(a,attr)())
					
					for k in range(c.get_zsize()):
						for j in range(c.get_ysize()):
							for i in range(c.get_xsize()):
								self.assertAlmostEqual(c.get_value_at(i,j,k), a.get_value_at(i,j,k), 3)
	
	def test_cuda_ccf(self):
		"""test cuda ccf equals cpu ccf ............................."""
		a = test_image(0,size=(32,32))
		b = a.calc_ccf(a)
		c = a.calc_ccf_cuda(a)
		c.process_inplace("xform.flip",{"axis":"x"})
		c.process_inplace("xform.flip",{"axis":"y"})
		#b.process_inplace("normalize")
		#c.process_inplace("normalize")
		for k in range(c.get_zsize()):
			for j in range(c.get_ysize()):
				for i in range(c.get_xsize()):
					self.assertAlmostEqual(c.get_value_at(i,j,k), b.get_value_at(i,j,k), 1)
					
	def test_cuda_2d_square_fft(self):
		"""test cuda 2D square fft equals cpu fft ..................."""
		for x in [15,16]:
			a = test_image(0,size=(x,x))
			b = a.do_fft()
			c = a.do_fft_cuda()
			for k in range(c.get_zsize()):
				for j in range(c.get_ysize()):
					for i in range(c.get_xsize()):
						self.assertAlmostEqual(c.get_value_at(i,j,k), b.get_value_at(i,j,k), 3)
						
	def test_cuda_3d_square_fft(self):
		"""test cuda 3D square fft equals cpu fft ..................."""
		for x in [15,16]:
			a = test_image_3d(0,size=(x,x,x))
			b = a.do_fft()
			c = a.do_fft_cuda()
			for k in range(c.get_zsize()):
				for j in range(c.get_ysize()):
					for i in range(c.get_xsize()):
						self.assertAlmostEqual(c.get_value_at(i,j,k), b.get_value_at(i,j,k), 3)
						
	def test_cuda_basic_mult(self):
		"""test cuda basic multiplication ..........................."""
		for x in [15,16]:
			a = EMData(x,x)
			a.process_inplace('testimage.noise.uniform.rand')
			b = a.copy()
			a.process_inplace("cuda.math.mult",{"scale":2.0})
			a.process_inplace("cuda.math.mult",{"scale":1.0/2.0})
			for k in range(a.get_zsize()):
				for j in range(a.get_ysize()):
					for i in range(a.get_xsize()):
						self.assertAlmostEqual(a.get_value_at(i,j,k), b.get_value_at(i,j,k), 8)
		
	
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
