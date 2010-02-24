#!/usr/bin/env python

#
# Author: David Woolford, 09/06/2007 (woolford@bcm.edu)
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
import unittest,os,sys
import testlib
from pyemtbx.exceptions import *
import numpy
from optparse import OptionParser

FFT_RECON_ZERO_TOLERANCE = 1.e-3
IS_TEST_EXCEPTION = False

# This testing code will perform many tests than it currently does.

class TestReconstructor(unittest.TestCase):
	"""Reconstructor test"""
	
	def notest_insert_remove(self):
		"""test insert and remove .........................."""
		n = 16
		recon_name = "fourier"
		params = {}
		params["x_in"] = n
		params["y_in"] = n
		params["mode"] = 2
		params["quiet"] = True
		#recon=Reconstructors.get(recon_name, params)
		#recon.setup()
		e = EMData()
		e.set_size(n,n,1)
		e.process_inplace('testimage.noise.uniform.rand')
		
		ee = e.copy()
		for mode in range(1,8):
			params["mode"] = mode
			for az in range(0,180,5):
				for alt in range(0,90,5):
					recon=Reconstructors.get(recon_name, params)
					recon.setup()
					transform = Transform3D(EULER_EMAN,az,alt,0)
					p = {}
					p["weight"] = 1.0
					recon.insert_params(p)
					recon.insert_slice(e,transform)
					#d = recon.get_emdata()
					#dd = d.get_fft_amplitude();
					#dd.write_image("resulta.mrc")
					
					
					p["weight"] = -1.0
					recon.insert_params(p)
					recon.insert_slice(ee,transform)
					
					d = recon.get_emdata()
					#dd = d.get_fft_amplitude();
					#dd.write_image("resultb.mrc")
					
					for k in range(d.get_zsize()):
						for j in range(d.get_ysize()):
							for i in range(d.get_xsize()):
								#print d.get_value_at(i,j,k)
								assert abs(d.get_value_at(i,j,k)) < FFT_RECON_ZERO_TOLERANCE
 
	def no_test_FourierReconstructor(self):
		"""test FourierReconstructor ........................"""
		n = 32
		e1 = EMData()
		e1.set_size(n,n,1)
		e1.process_inplace('testimage.noise.uniform.rand')
		#e1.do_fft_inplace()
		e2 = EMData()
		e2.set_size(n,n,1)
		e2.process_inplace('testimage.noise.uniform.rand')
		#e2.do_fft_inplace()
		e3 = EMData()
		e3.set_size(n,n,1)
		e3.process_inplace('testimage.noise.uniform.rand')
		#e3.do_fft_inplace()

		r = Reconstructors.get('fourier', {'size':(n,n), 'mode':'gauss_2', 'sym':'c1'})
		r.setup()
		r.insert_slice(e1, Transform({'type':'eman', 'alt':1.56, 'az':2.56, 'phi':3.56}))
		r.insert_slice(e2, Transform({'type':'eman', 'alt':2.56, 'az':3.56, 'phi':4.56}))
		r.insert_slice(e3, Transform({'type':'eman', 'alt':3.56, 'az':4.56, 'phi':5.56}))
		result = r.finish()
		
		testlib.safe_unlink('density.mrc')
	
	def no_test_WienerFourierReconstructor(self):
		"""test WienerFourierReconstructor .................."""
		a = 1
		e1 = test_image()
		e1.do_fft_inplace()
		e2 = test_image()
		e2.do_fft_inplace()
		e3 = test_image()
		e3.do_fft_inplace()
		
		r = Reconstructors.get('wiener_fourier', {'size':10, 'mode':1, 'padratio':0.5, 'snr':(0.2, 3.4, 5.6)})
		r.setup()
		r.insert_slice(e1, Transform({'type':'eman', 'alt':1.56, 'az':2.56, 'phi':3.56}))
		r.insert_slice(e2, Transform({'type':'eman', 'alt':1.56, 'az':2.56, 'phi':3.56}))
		r.insert_slice(e3, Transform({'type':'eman', 'alt':1.56, 'az':2.56, 'phi':3.56}))
		result = r.finish()
		
	def test_BackProjectionReconstructor(self):
		"""test BackProjectionReconstructor ................."""
		e1 = EMData()
		e1.set_size(32,32,1)
		e1.process_inplace('testimage.noise.uniform.rand')
		e2 = EMData()
		e2.set_size(32,32,1)
		e2.process_inplace('testimage.noise.uniform.rand')
		e3 = EMData()
		e3.set_size(32,32,1)
		e3.process_inplace('testimage.noise.uniform.rand')
		
		r = Reconstructors.get('back_projection', {'size':32, 'weight':0.8, 'sym':'c3'})
		r.setup()
		r.insert_slice(e1, Transform({'type':'eman', 'alt':1.56, 'az':2.56, 'phi':3.56}))
		r.insert_slice(e2, Transform({'type':'eman', 'alt':1.56, 'az':2.56, 'phi':3.56}))
		r.insert_slice(e3, Transform({'type':'eman', 'alt':1.56, 'az':2.56, 'phi':3.56}))
		result = r.finish(True)
		
	def no_test_nn4Reconstructor(self):
		"""test nn4Reconstructor ............................"""
		e1 = EMData()
		e1.set_size(34,32,32)
		e1.process_inplace('testimage.noise.uniform.rand')
		e2 = EMData()
		e2.set_size(34,32,32)
		e2.process_inplace('testimage.noise.uniform.rand')
		e3 = EMData()
		e3.set_size(34,32,32)
		e3.process_inplace('testimage.noise.uniform.rand')
		
		Log.logger().set_level(-1)    #no log message printed out
		r = Reconstructors.get('nn4', {'size':32, 'npad':1, 'symmetry':'CSYM'})
		r.setup()
		r.insert_slice(e1, Transform3D(EULER_EMAN, 0,0,0))
		r.insert_slice(e2, Transform3D(EULER_EMAN, 0,0,0))
		r.insert_slice(e3, Transform3D(EULER_EMAN, 0,0,0))
		result = r.finish()
	
	def no_test_ReverseGriddingReconstructor(self):
		"""test ReverseGriddingReconstructor ................"""
		e1 = EMData()
		e1.set_size(32,32,1)
		e1.process_inplace('testimage.noise.uniform.rand')
		e2 = EMData()
		e2.set_size(32,32,1)
		e2.process_inplace('testimage.noise.uniform.rand')
		e3 = EMData()
		e3.set_size(32,32,1)
		e3.process_inplace('testimage.noise.uniform.rand')
		
		r = Reconstructors.get('reverse_gridding', {'size':32, 'weight':0.8, 'npad':1})
		r.setup()
		r.insert_slice(e1, Transform3D(EULER_EMAN, 0,0,0))
		r.insert_slice(e2, Transform3D(EULER_EMAN, 0,0,0))
		r.insert_slice(e3, Transform3D(EULER_EMAN, 0,0,0))
		result = r.finish()

def test_main():
	p = OptionParser()
	p.add_option('--t', action='store_true', help='test exception', default=False )
	global IS_TEST_EXCEPTION
	opt, args = p.parse_args()
	if opt.t:
		IS_TEST_EXCEPTION = True
	suite = unittest.TestLoader().loadTestsFromTestCase(TestReconstructor)
	unittest.TextTestRunner(verbosity=2).run(suite)

if __name__ == '__main__':
    test_main()

