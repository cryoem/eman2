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

class TestAligner(unittest.TestCase):
    """aligner test"""
    
    def test_TranslationalAligner(self):
        """test TranslationalAligner ........................"""
        e = EMData()
        e.set_size(32,32,1)
        e.process_inplace('testimage.noise.uniform.rand')
        
        e2 = EMData()
        e2.set_size(32,32,1)
        e2.process_inplace('testimage.noise.uniform.rand')
   
        e.align('translational', e2, {})
        e.align('translational', e2, {'intonly':1, 'maxshift':2})
        
        #translational not support 3D image
        e3 = EMData()
        e3.set_size(32,32,32)
        e4 = EMData()
        e4.set_size(32,32,32)
        #self.assertRaises( RuntimeError, e3.align, 'translational', e4)
        
    def test_Translational3DAligner(self):
        """test Translational3DAligner ......................"""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        
        e2 = EMData()
        e2.set_size(32,32,32)
        e2.process_inplace('testimage.noise.uniform.rand')
   
        e.align('translational3d', e2, {"intonly":1})
   
    def test_RotationalAligner(self):
        """test RotationalAligner ..........................."""
        e = EMData()
        e.set_size(32,32,1)
        e.process_inplace('testimage.noise.uniform.rand')
        
        e2 = EMData()
        e2.set_size(32,32,1)
        e2.process_inplace('testimage.noise.uniform.rand')
        
        e.align('rotational', e2)

    def no_test_RotatePrecenterAligner(self):
        """test RotatePrecenterAligner ......................"""
        e = EMData()
        e.set_size(32,32,1)
        e.process_inplace('testimage.noise.uniform.rand')
        
        e2 = EMData()
        e2.set_size(32,32,1)
        e2.process_inplace('testimage.noise.uniform.rand')
        
        e.align('rotate_precenter', e2)
        
    def test_RotateCHAligner(self):
        """test RotateCHAligner ............................."""
        e = EMData()
        e.set_size(32,32,1)
        e.process_inplace('testimage.noise.uniform.rand')
        
        e2 = EMData()
        e2.set_size(32,32,1)
        e2.process_inplace('testimage.noise.uniform.rand')
        
        e.align('rotate_ch', e2, {'irad':1, 'orad':2})
        
    def test_RotateTranslateAligner(self):
        """test RotateTranslateAligner ......................"""
        e = EMData()
        e.set_size(32,32,1)
        e.process_inplace('testimage.noise.uniform.rand')
        
        e2 = EMData()
        e2.set_size(32,32,1)
        e2.process_inplace('testimage.noise.uniform.rand')
        
        e.align('rotate_translate', e2, {'maxshift':1})
        
    def no_test_RotateTranslateBestAligner(self):
        """test RotateTranslateBestAligner .................."""
        e = EMData()
        e.set_size(32,32,1)
        e.process_inplace('testimage.noise.uniform.rand')
        
        e2 = EMData()
        e2.set_size(32,32,1)
        e2.process_inplace('testimage.noise.uniform.rand')
        
        img = e.align('rotate_translate_best', e2, {'maxshift':1, 'snr':(1.0, 2.0, 3.0)})
        
    def test_RotateTranslateRadonAligner(self):
        """test RotateTranslateRadonAligner ................."""
        Log.logger().set_level(-1)    #no log message printed out
        e = EMData()
        e.set_size(32,32,1)
        e.process_inplace('testimage.noise.uniform.rand')
        
        e2 = EMData()
        e2.set_size(32,32,1)
        e2.process_inplace('testimage.noise.uniform.rand')
        
        e3 = EMData()
        e3.set_size(32,32,1)
        e3.process_inplace('testimage.noise.uniform.rand')
        
        e4 = EMData()
        e4.set_size(32,32,1)
        e4.process_inplace('testimage.noise.uniform.rand')
        
        img = e.align('rotate_translate_radon', e2, {'maxshift':2, 'radonwith':e3, 'radonthis':e4})
        
        import os
        testlib.safe_unlink('radon.hed')
        testlib.safe_unlink('radon.img')
        testlib.safe_unlink('racf.hed')
        testlib.safe_unlink('racf.img')
   
    def test_RotateFlipAligner(self):
        """test RotateFlipAligner ..........................."""
        e = EMData()
        e.set_size(32,32,1)
        e.process_inplace('testimage.noise.uniform.rand')
        
        e2 = EMData()
        e2.set_size(32,32,1)
        e2.process_inplace('testimage.noise.uniform.rand')
        
        e3 = EMData()
        e3.set_size(32,32,1)
        e3.process_inplace('testimage.noise.uniform.rand')
        
        e.align('rotate_flip', e2, {'flip':e3, 'imask':2})
        
    def test_RotateTranslateFlipAligner(self):
        """test RotateTranslateFlipAligner .................."""
        e = EMData()
        e.set_size(32,32,1)
        e.process_inplace('testimage.noise.uniform.rand')
        
        e2 = EMData()
        e2.set_size(32,32,1)
        e2.process_inplace('testimage.noise.uniform.rand')
        
        e3 = EMData()
        e3.set_size(32,32,1)
        e3.process_inplace('testimage.noise.uniform.rand')
        
        e.align('rotate_translate_flip', e2, {'flip':e3, 'usedot':1, 'maxshift':2})
        
    def no_test_RTFSlowAligner(self):
        """test RTFSlowAligner .............................."""
        e = EMData()
        e.set_size(32,32,1)
        e.process_inplace('testimage.noise.uniform.rand')
        
        e2 = EMData()
        e2.set_size(32,32,1)
        e2.process_inplace('testimage.noise.uniform.rand')
        
        e3 = EMData()
        e3.set_size(32,32,1)
        e3.process_inplace('testimage.noise.uniform.rand')
        
        e.align('rtf_slow', e2, {'flip':e3, 'maxshift':2})
        
        #RTFSlowestAligner eliminated
    def no_test_RTFSlowestAligner(self):
        """test RTFSlowestAligner ..........................."""
        e = EMData()
        e.set_size(32,32,1)
        e.process_inplace('testimage.noise.uniform.rand')
        
        e2 = EMData()
        e2.set_size(32,32,1)
        e2.process_inplace('testimage.noise.uniform.rand')
        
        e3 = EMData()
        e3.set_size(32,32,1)
        e3.process_inplace('testimage.noise.uniform.rand')
        
        e.align('rtf_slowest', e2, {'flip':e3, 'maxshift':2})
        
    def no_test_RTFBestAligner(self):
        """test RTFBestAligner .............................."""
        e = EMData()
        e.set_size(32,32,1)
        e.process_inplace('testimage.noise.uniform.rand')
        
        e2 = EMData()
        e2.set_size(32,32,1)
        e2.process_inplace('testimage.noise.uniform.rand')
        
        e3 = EMData()
        e3.set_size(32,32,1)
        e3.process_inplace('testimage.noise.uniform.rand')
        
        e.align('rtf_best', e2, {'flip':e3, 'maxshift':2, 'snr':(1.0, 2.0, 3.0)})
        
    def no_test_RTFRadonAligner(self):
        """test RTFRadonAligner ............................."""
        e = EMData()
        e.set_size(32,32,1)
        e.process_inplace('testimage.noise.uniform.rand')
        
        e2 = EMData()
        e2.set_size(32,32,1)
        e2.process_inplace('testimage.noise.uniform.rand')
        
        e3 = EMData()
        e3.set_size(32,32,1)
        e3.process_inplace('testimage.noise.uniform.rand')
        
        e4 = EMData()
        e4.set_size(32,32,1)
        e4.process_inplace('testimage.noise.uniform.rand')
        
        e5 = EMData()
        e5.set_size(32,32,1)
        e5.process_inplace('testimage.noise.uniform.rand')
        
        e6 = EMData()
        e6.set_size(32,32,1)
        e6.process_inplace('testimage.noise.uniform.rand')
   
        e.align('rtf_radon', e2, {'maxshift':2, 'thisf':e3, 'radonwith':e4, \
                'radonthis':e5, 'radonthisf':e6})
                
    def test_RefineAligner(self):
        """test RefineAligner ..............................."""
        e = EMData()
        e.set_size(32,32,1)
        e.process_inplace('testimage.noise.uniform.rand')
        
        e2 = EMData()
        e2.set_size(32,32,1)
        e2.process_inplace('testimage.noise.uniform.rand')
        
        e.align('refine', e2)
        e.align('refine', e2, {'mode':1, 'az':1.2, 'dx':2, 'dy':3.4})
   
def test_main():
    test_support.run_unittest(TestAligner)

if __name__ == '__main__':
    test_main()
    