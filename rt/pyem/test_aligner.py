#!/usr/bin/env    python

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
        e.process('testimage.noise.uniform.rand')
        
        e2 = EMData()
        e2.set_size(32,32,1)
        e2.process('testimage.noise.uniform.rand')
   
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
        e.process('testimage.noise.uniform.rand')
        
        e2 = EMData()
        e2.set_size(32,32,32)
        e2.process('testimage.noise.uniform.rand')
   
        e.align('translational3d', e2)
   
    def test_RotationalAligner(self):
        """test RotationalAligner ..........................."""
        e = EMData()
        e.set_size(32,32,1)
        e.process('testimage.noise.uniform.rand')
        
        e2 = EMData()
        e2.set_size(32,32,1)
        e2.process('testimage.noise.uniform.rand')
        
        e.align('rotational', e2)

    def test_RotatePrecenterAligner(self):
        """test RotatePrecenterAligner ......................"""
        e = EMData()
        e.set_size(32,32,1)
        e.process('testimage.noise.uniform.rand')
        
        e2 = EMData()
        e2.set_size(32,32,1)
        e2.process('testimage.noise.uniform.rand')
        
        e.align('rotate_precenter', e2)
        
    def test_RotateCHAligner(self):
        """test RotateCHAligner ............................."""
        e = EMData()
        e.set_size(32,32,1)
        e.process('testimage.noise.uniform.rand')
        
        e2 = EMData()
        e2.set_size(32,32,1)
        e2.process('testimage.noise.uniform.rand')
        
        e.align('rotate_ch', e2, {'irad':1, 'orad':2})
        
    def test_RotateTranslateAligner(self):
        """test RotateTranslateAligner ......................"""
        e = EMData()
        e.set_size(32,32,1)
        e.process('testimage.noise.uniform.rand')
        
        e2 = EMData()
        e2.set_size(32,32,1)
        e2.process('testimage.noise.uniform.rand')
        
        e.align('rotate_translate', e2, {'maxshift':1})
        
    def test_RotateTranslateBestAligner(self):
        """test RotateTranslateBestAligner .................."""
        e = EMData()
        e.set_size(32,32,1)
        e.process('testimage.noise.uniform.rand')
        
        e2 = EMData()
        e2.set_size(32,32,1)
        e2.process('testimage.noise.uniform.rand')
        
        img = e.align('rotate_translate_best', e2, {'maxshift':1, 'snr':(1.0, 2.0, 3.0)})
        
    def no_test_RotateTranslateRadonAligner(self):
        """test RotateTranslateRadonAligner ................."""
        e = EMData()
        e.set_size(32,32,1)
        e.process('testimage.noise.uniform.rand')
        
        e2 = EMData()
        e2.set_size(32,32,1)
        e2.process('testimage.noise.uniform.rand')
        
        e3 = EMData()
        e3.set_size(32,32,1)
        e3.process('testimage.noise.uniform.rand')
        
        e4 = EMData()
        e4.set_size(32,32,1)
        e4.process('testimage.noise.uniform.rand')
        
        img = e.align('rotate_translate_radon', e2, {'maxshift':2, 'radonwith':e3, 'radonthis':e4})
   
    def test_RotateFlipAligner(self):
        """test RotateFlipAligner ..........................."""
        e = EMData()
        e.set_size(32,32,1)
        e.process('testimage.noise.uniform.rand')
        
        e2 = EMData()
        e2.set_size(32,32,1)
        e2.process('testimage.noise.uniform.rand')
        
        e3 = EMData()
        e3.set_size(32,32,1)
        e3.process('testimage.noise.uniform.rand')
        
        e.align('rotate_flip', e2, {'flip':e3, 'imask':2})
        
    def test_RotateTranslateFlipAligner(self):
        """test RotateTranslateFlipAligner .................."""
        e = EMData()
        e.set_size(32,32,1)
        e.process('testimage.noise.uniform.rand')
        
        e2 = EMData()
        e2.set_size(32,32,1)
        e2.process('testimage.noise.uniform.rand')
        
        e3 = EMData()
        e3.set_size(32,32,1)
        e3.process('testimage.noise.uniform.rand')
        
        e.align('rotate_translate_flip', e2, {'flip':e3, 'usedot':1, 'maxshift':2})
        
    def test_RTFSlowAligner(self):
        """test RTFSlowAligner .............................."""
        e = EMData()
        e.set_size(32,32,1)
        e.process('testimage.noise.uniform.rand')
        
        e2 = EMData()
        e2.set_size(32,32,1)
        e2.process('testimage.noise.uniform.rand')
        
        e3 = EMData()
        e3.set_size(32,32,1)
        e3.process('testimage.noise.uniform.rand')
        
        e.align('rtf_slow', e2, {'flip':e3, 'maxshift':2})
        
    def no_test_RTFSlowestAligner(self):
        """test RTFSlowestAligner ..........................."""
        e = EMData()
        e.set_size(32,32,1)
        e.process('testimage.noise.uniform.rand')
        
        e2 = EMData()
        e2.set_size(32,32,1)
        e2.process('testimage.noise.uniform.rand')
        
        e3 = EMData()
        e3.set_size(32,32,1)
        e3.process('testimage.noise.uniform.rand')
        
        e.align('rtf_slowest', e2, {'flip':e3, 'maxshift':2})
        
    def test_RTFBestAligner(self):
        """test RTFBestAligner .............................."""
        e = EMData()
        e.set_size(32,32,1)
        e.process('testimage.noise.uniform.rand')
        
        e2 = EMData()
        e2.set_size(32,32,1)
        e2.process('testimage.noise.uniform.rand')
        
        e3 = EMData()
        e3.set_size(32,32,1)
        e3.process('testimage.noise.uniform.rand')
        
        e.align('rtf_best', e2, {'flip':e3, 'maxshift':2, 'snr':(1.0, 2.0, 3.0)})
        
    def no_test_RTFRadonAligner(self):
        """test RTFRadonAligner ............................."""
        e = EMData()
        e.set_size(32,32,1)
        e.process('testimage.noise.uniform.rand')
        
        e2 = EMData()
        e2.set_size(32,32,1)
        e2.process('testimage.noise.uniform.rand')
        
        e3 = EMData()
        e3.set_size(32,32,1)
        e3.process('testimage.noise.uniform.rand')
        
        e4 = EMData()
        e4.set_size(32,32,1)
        e4.process('testimage.noise.uniform.rand')
        
        e5 = EMData()
        e5.set_size(32,32,1)
        e5.process('testimage.noise.uniform.rand')
        
        e6 = EMData()
        e6.set_size(32,32,1)
        e6.process('testimage.noise.uniform.rand')
   
        e.align('rtf_radon', e2, {'maxshift':2, 'thisf':e3, 'radonwith':e4, \
                'radonthis':e5, 'radonthisf':e6})
                
    def no_test_RefineAligner(self):
        """test RefineAligner ..............................."""
        e = EMData()
        e.set_size(32,32,1)
        e.process('testimage.noise.uniform.rand')
        
        e2 = EMData()
        e2.set_size(32,32,1)
        e2.process('testimage.noise.uniform.rand')
        
        e.align('refine', e2, {'mode':1, 'snr':(1.0, 2.0, 3.0), 'alot':1.2, \
                'az':1.2, 'phi':2.3, 'dx':2, 'dy':3.4, 'dz':2.6})
   
def test_main():
    test_support.run_unittest(TestAligner)

if __name__ == '__main__':
    test_main()
    