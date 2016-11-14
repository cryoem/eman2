#!/usr/bin/env python

#
# Author: Grant Tang, 09/03/2005 (gtang@bcm.edu)
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
from optparse import OptionParser

IS_TEST_EXCEPTION = False

class TestCmp(unittest.TestCase):
    """cmp test"""
    
    def test_variance(self):
        """test variance ...................................."""
        imgfile1 = "test_variance_1.mrc"
        TestUtil.make_image_file(imgfile1, IMAGE_MRC)
        e1 = EMData()
        e1.read_image(imgfile1)

        e2 = e1.copy()
        score = e2.cmp("sqeuclidean", e1)
        self.assertEqual(score, 0)
        testlib.safe_unlink(imgfile1)
        
    def test_basic_cmp(self):
        """test basic cmp ..................................."""
        e1 = EMData(16,16)
        e1.process_inplace('testimage.noise.uniform.rand')

        e2 = EMData(16,16)
        e2.process_inplace('testimage.noise.uniform.rand')

        dot_score = e2.cmp("dot", e1, {"negative":0, "normalize":1})

        variance_score = e2.cmp("sqeuclidean", e1)
        
        phase_score = e2.cmp("phase", e1, {})
        
        e1.do_fft_inplace()
        e2.do_fft_inplace()
        frc_score = e2.cmp("frc", e1, {})
    
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
        e.do_fft_inplace()
        
        e2 = EMData()
        e2.set_size(64,64,1)
        e2.process_inplace('testimage.noise.uniform.rand')
        e2.do_fft_inplace()
        
        score  = e.cmp('frc', e2, {})
        #score2 = e.cmp('frc', e2, {'snr':(1.0, 2.0)})    #todo: segmentation fault
        
        # TODO - THE FRC SEEMS BROKEN, THIS TEST CODE FAILS - d.woolford
        # the frc (integral) of an image compared to itself should always be (in EMAN2 negative) 1
        # Here this assertion is tested for all combinations of all even odd combinations
        # of 2D and 3D images - the image tested against is random noise
        #n = 16
        #for i in range(n-1,n+1):
		#	for j in range(n-1,n+1):
		#		for k in [1,n-1,n]:
		#			e3 = EMData()
		#			e3.set_size(i,j,k)
		#			e3.process_inplace('testimage.noise.uniform.rand')
        #           e3.do_fft_inplace()
        #		    neg_one  = e3.cmp('frc', e3.copy(), {})
		#		    self.assertAlmostEqual(neg_one,-1, places=6)
        
    def test_PhaseCmp(self):
        """test PhaseCmp ...................................."""
        #THIS TEST WILL BE FIXED SOON BY DAVE
        e = EMData()
        e.set_size(64,64,1)
        e.process_inplace('testimage.noise.uniform.rand')
        
        e2 = EMData()
        e2.set_size(64,64,1)
        e2.process_inplace('testimage.noise.uniform.rand')
        
        score  = e.cmp('phase', e2, {})
        # the phase residual of an image compared to itself should always be zero
        # Here this assertion is tested for all combinations of all even odd combinations
        # of 2D and 3D images - the image tested against is random noise
        n = 16
        for i in range(n-1,n+1):
			for j in range(n-1,n+1):
				for k in [1,n-1,n]:
					e3 = EMData()
					e3.set_size(i,j,k)
					#print "%d %d %d" %(i,j,k)
					e3.process_inplace('testimage.noise.uniform.rand')
					f = e3.copy()
					zero  = e3.cmp('phase', f, {})
					#print "%d %d %d %f" %(i,j,k,zero)
					self.assertAlmostEqual(zero,0, places=2)

        
    def test_SqEuclideanCmp(self):
        """test SqEuclideanCmp .............................."""
        e = EMData()
        e.set_size(64,64,1)
        e.process_inplace('testimage.noise.uniform.rand')
        
        e2 = EMData()
        e2.set_size(64,64,1)
        e2.process_inplace('testimage.noise.uniform.rand')
        
        score  = e.cmp('sqeuclidean', e2, {})
        # the square euclidiean distance difference of an image and itself should always be zero
        zero  = e2.cmp('sqeuclidean', e2, {})

        # the square euclidiean distance difference of an image and itself should always be zero
        # Here this assertion is tested for all combinations of all even odd combinations
        # of 2D and 3D images - the image tested against is random noise
        n = 16
        for i in range(n-1,n+1):
			for j in range(n-1,n+1):
				for k in [1,n-1,n]:
					e3 = EMData()
					e3.set_size(i,j,k)
					e3.process_inplace('testimage.noise.uniform.rand')
					zero  = e3.cmp('sqeuclidean', e3.copy(), {})
					self.assertAlmostEqual(zero,0, places=6)
        
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
        n = 16
        # the normalized dot product of an image and itself should always be -1 (by default the dot product is negated)
        # Here this assertion is tested for all combinations of all even odd combinations
        # of 2D and 3D images - the image tested against is random noise
        n = 16
        for i in range(n-1,n+1):
			for j in range(n-1,n+1):
				for k in [1,n-1,n]:
					e3 = EMData()
					e3.set_size(i,j,k)
					e3.process_inplace('testimage.noise.uniform.rand')
					neg_one  = e3.cmp('dot', e3.copy(), {"normalize":1})
					self.assertAlmostEqual(neg_one,-1, places=6)
        

def test_main():
    p = OptionParser()
    p.add_option('--t', action='store_true', help='test exception', default=False )
    global IS_TEST_EXCEPTION
    opt, args = p.parse_args()
    if opt.t:
        IS_TEST_EXCEPTION = True
    Log.logger().set_level(-1)  #perfect solution for quenching the Log error information, thank Liwei
    suite = unittest.TestLoader().loadTestsFromTestCase(TestCmp)
    unittest.TextTestRunner(verbosity=2).run(suite)

if __name__ == '__main__':
    test_main()
