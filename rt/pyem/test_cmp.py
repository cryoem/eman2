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

class TestCmp(unittest.TestCase):
    """cmp test"""
    
    def test_variance(self):
        """test variance ...................................."""
        imgfile1 = "test_variance_1.mrc"
        TestUtil.make_image_file(imgfile1, IMAGE_MRC)
        e1 = EMData()
        e1.read_image(imgfile1)

        e2 = e1.copy()
        score = e2.cmp("SqEuclidean", e1)
        self.assertEqual(score, 0)
        testlib.safe_unlink(imgfile1)
        
    def test_basic_cmp(self):
        """test basic cmp ..................................."""
        imgfile1 = "test_basic_cmp_1.hed"
        TestUtil.make_image_file(imgfile1, IMAGE_IMAGIC, EM_FLOAT, 16,16,4)

        e1 = EMData()
        e1.read_image(imgfile1, 1)

        e2 = EMData()
        e2.read_image(imgfile1, 2)

        #e1.write_image("test_basic_cmp_out_1.mrc")
        #e2.write_image("test_basic_cmp_out_2.mrc")

        dot_score = e2.cmp("dot", e1, {"negative":0, "normalize":1})
#        self.assertEqual(dot_score, 19944.0)    #todo: dot score not match, anything wrong?

        variance_score = e2.cmp("SqEuclidean", e1)
#        self.assertEqual(variance_score, 0)    #todo: score not match, anything wrong?
        
        phase_score = e2.cmp("phase", e1, {})
#        testlib.assertfloat(self, phase_score, 1.6488)    #todo: score not match, anything wrong?
        
        frc_score = e2.cmp("frc", e1, {})
#        testlib.assertfloat(self, frc_score, -0.4011)    #todo: score not match, anything wrong?

        (hed1,img1) = testlib.get_imagic_filename_pair(imgfile1)
        testlib.safe_unlink(hed1)
        testlib.safe_unlink(img1)
    
    def test_QuadMinDotCmp(self):
        """test QuadMinDotCmp ..............................."""
        e = EMData()
        e.set_size(64,64,1)
        e.process_inplace('testimage.noise.uniform.rand')
        
        e2 = EMData()
        e2.set_size(64,64,1)
        e2.process_inplace('testimage.noise.uniform.rand')
        
        score  = e.cmp('quadmindot', e2, {})    #default argument
        score2 = e.cmp('quadmindot', e2, {'negative':0, 'normalize':1})
        
    def test_OptVarianceCmp(self):
        """test OptVarianceCmp .............................."""
        e = EMData()
        e.set_size(64,64,1)
        e.process_inplace('testimage.noise.uniform.rand')
        
        e2 = EMData()
        e2.set_size(64,64,1)
        e2.process_inplace('testimage.noise.uniform.rand')
        
        score  = e.cmp('optvariance', e2, {})    #default argument
        score2 = e.cmp('optvariance', e2, {'invert':1, 'keepzero':1, 'matchfilt':2, 'radweight':2, 'debug':1})
        
        testlib.safe_unlink('a.hdf')
        testlib.safe_unlink('dbug.optvar.txt')
        
    def test_FRCCmp(self):
        """test FRCCmp ......................................"""
        e = EMData()
        e.set_size(64,64,1)
        e.process_inplace('testimage.noise.uniform.rand')
        
        e2 = EMData()
        e2.set_size(64,64,1)
        e2.process_inplace('testimage.noise.uniform.rand')
        
        score  = e.cmp('frc', e2, {})
        #score2 = e.cmp('frc', e2, {'snr':(1.0, 2.0)})    #todo: segmentation fault
        
    def test_PhaseCmp(self):
        """test PhaseCmp ...................................."""
        e = EMData()
        e.set_size(64,64,1)
        e.process_inplace('testimage.noise.uniform.rand')
        
        e2 = EMData()
        e2.set_size(64,64,1)
        e2.process_inplace('testimage.noise.uniform.rand')
        
        score  = e.cmp('phase', e2, {})
        
    def test_SqEuclideanCmp(self):
        """test SqEuclideanCmp .............................."""
        e = EMData()
        e.set_size(64,64,1)
        e.process_inplace('testimage.noise.uniform.rand')
        
        e2 = EMData()
        e2.set_size(64,64,1)
        e2.process_inplace('testimage.noise.uniform.rand')
        
        score  = e.cmp('SqEuclidean', e2, {})
        
    def test_DotCmp(self):
        """test DotCmp ......................................"""
        e = EMData()
        e.set_size(64,64,1)
        e.process_inplace('testimage.noise.uniform.rand')
        
        e2 = EMData()
        e2.set_size(64,64,1)
        e2.process_inplace('testimage.noise.uniform.rand')
        
        score  = e.cmp('dot', e2, {})
        score  = e.cmp('dot', e2, {'negative':1, 'normalize':1})
def test_main():
    test_support.run_unittest(TestCmp)

if __name__ == '__main__':
    test_main()