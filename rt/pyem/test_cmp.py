#!/usr/bin/env    python

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
        TestUtil.make_image_file(imgfile1, MRC)
        e1 = EMData()
        e1.read_image(imgfile1)

        e2 = e1.copy()
        score = e2.cmp("SqEuclidean", e1)
        self.assertEqual(score, 0)
        testlib.safe_unlink(imgfile1)
        
    def test_basic_cmp(self):
        """test basic cmp ..................................."""
        imgfile1 = "test_basic_cmp_1.hed"
        TestUtil.make_image_file(imgfile1, IMAGIC, EM_FLOAT, 16,16,4)

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