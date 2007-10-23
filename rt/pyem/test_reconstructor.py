#!/usr/bin/env python

#
# Author: Steven Ludtke, 04/10/2003 (sludtke@bcm.edu)
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
from test import test_support
import testlib
from pyemtbx.exceptions import *
import numpy

FFT_RECON_ZERO_TOLERANCE = 1.e-3

# This testing code will perform many tests than it currently does.

class TestReconstructor(unittest.TestCase):
	"""Reconstructor test"""
	
	def test_insert_remove(self):
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
 
def test_main():
    test_support.run_unittest(TestReconstructor)

if __name__ == '__main__':
    test_main()
