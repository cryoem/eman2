#!/usr/bin/env    python

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

class TestAverager(unittest.TestCase):
    """averager test"""
    
    def test_ImageAverager(self):
        """test ImageAverager ..............................."""
        
class TestConstructor(unittest.TestCase):
    """constructor test"""
    
    def test_FourierReconstructor(self):
        """test FourierReconstructor ........................"""
        e1 = EMData()
        e1.set_size(32,32,1)
        e1.process_inplace('testimage.noise.uniform.rand')
        #e1.do_fft_inplace()
        e2 = EMData()
        e2.set_size(32,32,1)
        e2.process_inplace('testimage.noise.uniform.rand')
        #e2.do_fft_inplace()
        e3 = EMData()
        e3.set_size(32,32,1)
        e3.process_inplace('testimage.noise.uniform.rand')
        #e3.do_fft_inplace()
        
        r = Reconstructors.get('fourier', {'size':10, 'mode':1, 'weight':0.5, 'dlog':2, 'sym':'DSYM'})
        r.setup()
        r.insert_slice(e1, Transform3D(EULER_EMAN, 0,0,0))
        r.insert_slice(e2, Transform3D(EULER_EMAN, 0,0,0))
        r.insert_slice(e3, Transform3D(EULER_EMAN, 0,0,0))
        result = r.finish()
    
    def test_WienerFourierReconstructor(self):
        """test WienerFourierReconstructor .................."""
        e1 = EMData()
        e1.set_size(32,32,1)
        e1.process_inplace('testimage.noise.uniform.rand')
        e1.do_fft_inplace()
        e2 = EMData()
        e2.set_size(32,32,1)
        e2.process_inplace('testimage.noise.uniform.rand')
        e2.do_fft_inplace()
        e3 = EMData()
        e3.set_size(32,32,1)
        e3.process_inplace('testimage.noise.uniform.rand')
        e3.do_fft_inplace()
        
        r = Reconstructors.get('wiener_fourier', {'size':10, 'mode':1, 'padratio':0.5, 'snr':(0.2, 3.4, 5.6)})
        r.setup()
        r.insert_slice(e1, Transform3D(EULER_EMAN, 0,0,0))
        r.insert_slice(e2, Transform3D(EULER_EMAN, 0,0,0))
        r.insert_slice(e3, Transform3D(EULER_EMAN, 0,0,0))
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
        
        r = Reconstructors.get('back_projection', {'size':32, 'weight':0.8})
        r.setup()
        r.insert_slice(e1, Transform3D(EULER_EMAN, 0,0,0))
        r.insert_slice(e2, Transform3D(EULER_EMAN, 0,0,0))
        r.insert_slice(e3, Transform3D(EULER_EMAN, 0,0,0))
        result = r.finish()
        
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

class TestProjector(unittest.TestCase):
    """test Projector"""
    
    def test_GaussFFTProjector(self):
        """test GaussFFTProjector ..........................."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        
        e.project('gauss_fft', {'alt':1.234, 'az':1.345, 'phi':1.54, 'mode':1})

        
def test_main():
    test_support.run_unittest(TestAverager, TestConstructor, TestProjector)

if __name__ == '__main__':
    test_main()