#!/usr/bin/env python

from EMAN2 import *
import unittest,os,sys
from test import test_support
import testlib
from pyemtbx.exceptions import *
import Numeric

class TestProcessor(unittest.TestCase):
    """Processor test"""
    
    def test_get_processor_list(self):
        """test get processor list .........................."""
        processor_names = Processors.get_list()
        self.assertEqual(len(processor_names), 116)

        try:
            f2 = Processors.get("_nosuchfilter___")
        except RuntimeError, runtime_err:
            err_type = exception_type(runtime_err)
            self.assertEqual(err_type, "NotExistingObjectException")
            
    def test_filter_lowpass_tophat(self):
        """test filter.lowpass.tophat processor ............."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        e.do_fft_inplace()
        self.assertEqual(e.is_complex(), True)
        
        e.process_inplace('filter.highpass.tophat', {'cutoff_frequency':0.3})
        
    def test_filter_highpass_tophat(self):
        """test filter.lowpass.tophat processor ............."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        e.do_fft_inplace()
        self.assertEqual(e.is_complex(), True)
        
        e.process_inplace('filter.lowpass.tophat', {'cutoff_frequency':0.3})
        
    def test_filter_bandpass_tophat(self):
        """test filter.bandpass.tophat processor ............"""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        e.do_fft_inplace()
        self.assertEqual(e.is_complex(), True)
        
        e.process_inplace('filter.bandpass.tophat', \
            {'low_cutoff_frequency':0.1, 'high_cutoff_frequency':0.45})
            
    def test_filter_homomorphic_tophat(self):
        """test filter.homomorphic.tophat processor ........."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        e.do_fft_inplace()
        self.assertEqual(e.is_complex(), True)
        
        e.process_inplace('filter.homomorphic.tophat', \
            {'low_cutoff_frequency':0.05, 'high_cutoff_frequency':0.45, 'value_at_zero_frequency':1.0})

    def test_filter_lowpass_gauss(self):
        """test filter.lowpass.gauss processor .............."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        e.do_fft_inplace()
        self.assertEqual(e.is_complex(), True)
        
        e.process_inplace('filter.lowpass.gauss', {'sigma':0.5})
        
    def test_filter_highpass_gauss(self):
        """test filter.highpass.gauss processor ............."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        e.do_fft_inplace()
        self.assertEqual(e.is_complex(), True)
        
        e.process_inplace('filter.highpass.gauss', {'sigma':0.5})
        
    def test_filter_bandpass_gauss(self):
        """test filter.bandpass.gauss processor ............."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        e.do_fft_inplace()
        self.assertEqual(e.is_complex(), True)
        
        e.process_inplace('filter.bandpass.gauss', {'sigma':0.5, 'center':0.1})
        
    def test_filter_homomorphic_gauss(self):
        """test filter.homomorphic.gauss processor .........."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        e.do_fft_inplace()
        self.assertEqual(e.is_complex(), True)
        
        e.process_inplace('filter.homomorphic.gauss', {'sigma':0.5, 'value_at_zero_frequency':1.0})
        
    def test_filter_gaussinverse(self):
        """test filter.gaussinverse processor ..............."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        e.do_fft_inplace()
        self.assertEqual(e.is_complex(), True)
        
        e.process_inplace('filter.gaussinverse', {'sigma':0.5})
        
    #gsl: bessel_I0.c:216: ERROR: overflow
    #Default GSL error handler invoked. Aborted
    def no_test_filter_kaiserI0inverse(self):
        """test filter.kaiser_io_inverse processor .........."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        e.do_fft_inplace()
        self.assertEqual(e.is_complex(), True)
        
        e.process_inplace('filter.kaiser_io_inverse')
        
    def test_filter_kaisersinhinverse(self):
        """test filter.kaisersinhinverse processor .........."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        e.do_fft_inplace()
        self.assertEqual(e.is_complex(), True)
        
        e.process_inplace('filter.kaisersinhinverse')
        
    def test_filter_radialtable(self):
        """test filter.radialtable processor ................"""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        e.do_fft_inplace()
        self.assertEqual(e.is_complex(), True)
        
        e.process_inplace('filter.radialtable', {'table':(0.2,0.2,0.3)})
        
    def test_filter_lowpass_butterworth(self):
        """test filter.lowpass.butterworth processor ........"""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        e.do_fft_inplace()
        self.assertEqual(e.is_complex(), True)
        
        e.process_inplace('filter.lowpass.butterworth', \
            {'low_cutoff_frequency':0.1, 'high_cutoff_frequency':0.45})
            
    def test_filter_highpass_butterworth(self):
        """test filter.highpass.butterworth processor ......."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        e.do_fft_inplace()
        self.assertEqual(e.is_complex(), True)
        
        e.process_inplace('filter.highpass.butterworth', \
            {'low_cutoff_frequency':0.1, 'high_cutoff_frequency':0.45})
            
    def test_filter_homomorphic_butterworth(self):
        """test filter.homomorphic.butterworth processor ...."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        e.do_fft_inplace()
        self.assertEqual(e.is_complex(), True)
        
        e.process_inplace('filter.homomorphic.butterworth', \
            {'low_cutoff_frequency':0.1, 'high_cutoff_frequency':0.45, 'value_at_zero_frequency':1.0})
            
    def test_filter_lowpass_tanh(self):
        """test filter.lowpass.tanh processor ..............."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        e.do_fft_inplace()
        self.assertEqual(e.is_complex(), True)
        
        e.process_inplace('filter.lowpass.tanh', {'cutoff_frequency':0.25, 'fall_off':0.1})
        
    def test_filter_highpass_tanh(self):
        """test filter.highpass.tanh processor .............."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        e.do_fft_inplace()
        self.assertEqual(e.is_complex(), True)
        
        e.process_inplace('filter.highpass.tanh', {'cutoff_frequency':0.25, 'fall_off':0.1})
        
    def test_filter_homomorphic_tanh(self):
        """test filter.homomorphic.tanh processor ..........."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        e.do_fft_inplace()
        self.assertEqual(e.is_complex(), True)
        
        e.process_inplace('filter.homomorphic.tanh', \
            {'cutoff_frequency':0.25, 'fall_off':0.1, 'value_at_zero_frequency':1.0})
            
    def test_filter_bandpass_tanh(self):
        """test filter.bandpass.tanh processor .............."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        e.do_fft_inplace()
        self.assertEqual(e.is_complex(), True)
        
        e.process_inplace('filter.bandpass.tanh', \
            {'low_cutoff_frequency':0.1, 'Low_fall_off':0.15, 'high_cutoff_frequency':0.45, \
             'high_fall_off':0.15, 'fall_off':1.25})
             
    def test_eman1_filter_lowpass_sharp(self):
        """test eman1.filter.lowpass.sharp processor ........"""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        e.do_fft_inplace()
        self.assertEqual(e.is_complex(), True)
        
        e.process_inplace('eman1.filter.lowpass.sharp', {'lowpass':0.1})
        
    def test_eman1_filter_highpass_sharp(self):
        """test eman1.filter.highpass.sharp processor ......."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        e.do_fft_inplace()
        self.assertEqual(e.is_complex(), True)
        
        e.process_inplace('eman1.filter.highpass.sharp', {'highpass':0.45})
        
    def test_eman1_filter_lowpass_gaussian(self):
        """test eman1.filter.lowpass.gaussian processor ....."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        e.do_fft_inplace()
        self.assertEqual(e.is_complex(), True)
        
        e.process_inplace('eman1.filter.lowpass.gaussian', {'lowpass':0.1})
        
    def test_eman1_filter_highpass_gaussian(self):
        """test eman1.filter.highpass.gaussian processor ...."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        e.do_fft_inplace()
        self.assertEqual(e.is_complex(), True)
        
        e.process_inplace('eman1.filter.highpass.gaussian', {'highpass':0.45})
        
    def test_eman1_filter_lowpass_tanh(self):
        """test eman1.filter.lowpass.tanh processor ........."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        e.do_fft_inplace()
        self.assertEqual(e.is_complex(), True)
        
        e.process_inplace('eman1.filter.lowpass.tanh', {'lowpass':0.1})
        
    def test_eman1_filter_highpass_tanh(self):
        """test eman1.filter.highpass.tanh processor ........"""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        e.do_fft_inplace()
        self.assertEqual(e.is_complex(), True)
        
        e.process_inplace('eman1.filter.highpass.tanh', {'highpass':0.45})
        
    def test_eman1_filter_highpass_butterworth(self):
        """test eman1.filter.highpass.butterworth processor ."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        e.do_fft_inplace()
        self.assertEqual(e.is_complex(), True)
        
        e.process_inplace('eman1.filter.highpass.butterworth', {'highpass':0.45})
        
    def test_eman1_filter_ramp(self):
        """test eman1.filter.ramp processor ................."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        e.do_fft_inplace()
        self.assertEqual(e.is_complex(), True)
        
        e.process_inplace('eman1.filter.ramp', {'intercept':0.25, 'slope':0.3})
        
    def test_eman1_math_absvalue(self):
        """test eman1.math.absvalue processor ..............."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        d = e.get_3dview()
        e2 = e * -1
        d2 = e2.get_3dview()
        for x in range(32):
            for y in range(32):
                for z in range(32):
                    self.assertAlmostEqual(d[z][y][x], -d2[z][y][x], 3)
        
        e2.process_inplace('eman1.math.absvalue')
        d2 = e2.get_3dview()
        for x in range(32):
            for y in range(32):
                for z in range(32):
                    self.assertAlmostEqual(d[z][y][x], d2[z][y][x], 3)
                    
    def test_eman1_threshold_notzero(self):
        """test eman1.threshold.notzero processor ..........."""
        e = EMData()
        e.set_size(32,32,32)
        
        e.to_zero()
        e.process_inplace('eman1.threshold.notzero')
        d = e.get_3dview()
        for x in range(32):
            for y in range(32):
                for z in range(32):
                    self.assertEqual(d[z][y][x], 0)
        
        e.to_one()
        e *= 0.5
        e.process_inplace('eman1.threshold.notzero')
        d = e.get_3dview()
        for x in range(32):
            for y in range(32):
                for z in range(32):
                    self.assertEqual(d[z][y][x], 1)
   
    def test_eman1_math_squared(self):
        """test eman1.math.squared processor ................"""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        d = e.get_3dview()
        e2 = e.copy()
        e2.process_inplace('eman1.math.squared')
        d2 = e2.get_3dview()
        for x in range(32):
            for y in range(32):
                for z in range(32):
                    self.assertAlmostEqual(d[z][y][x]**2, d2[z][y][x], 3)
                    
    def test_eman1_math_sqrt(self):
        """test eman1.math.sqrt processor ..................."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        d = e.get_3dview()
        e2 = e.copy()
        e2.process_inplace('eman1.math.sqrt')
        d2 = e2.get_3dview()
        from math import sqrt
        for x in range(32):
            for y in range(32):
                for z in range(32):
                    self.assertAlmostEqual(sqrt(d[z][y][x]), d2[z][y][x], 3) 
                    
    def test_eman1_threshold_belowtozero(self):
        """test eman1.threshold.belowtozero processor ......."""
        e = EMData()
        e.set_size(32,32,32)
        e.to_one()
        e *= 0.4
        
        e.process_inplace('eman1.threshold.belowtozero', {'minval':0.5})
        d = e.get_3dview()
        for x in range(32):
            for y in range(32):
                for z in range(32):
                    self.assertAlmostEqual(d[z][y][x], 0, 3)
                            
    def test_BinarizeProcessor(self):
        """test binary processor ............................"""
        imgfile1 = "test_BinarizeFilter.mrc"
        TestUtil.make_image_file(imgfile1, MRC)
        e = EMData()
        e.read_image(imgfile1)
        fnum = 1000
        f1 = Processors.get("eman1.threshold.binary", {'value': fnum})
        new_params = f1.get_params()
        self.assertEqual(float(new_params["value"]), fnum)
        f1.process_inplace(e)
        testlib.check_emdata(e, sys.argv[0])
        os.unlink(imgfile1)

    def test_eman1_threshold_compress(self):
        """test eman1.threshold.compress processor .........."""
        e = EMData()
        e.set_size(32,32,32)
        e.to_one()
        e *= 0.45
        e.process_inplace('eman1.threshold.compress', {'range':0.2, 'value':0.5})
        d = e.get_3dview()
        for x in range(32):
            for y in range(32):
                for z in range(32):
                    self.assertAlmostEqual(d[z][y][x], 0.5, 3)
        
        e2 = EMData()
        e2.set_size(32,32,32)
        e2.to_one()
        e2 *= 0.05
        e2.process_inplace('eman1.threshold.compress', {'range':0.2, 'value':0.5})
        d2 = e2.get_3dview()
        for x in range(32):
            for y in range(32):
                for z in range(32):
                    self.assertAlmostEqual(d2[z][y][x], 0.2+0.05, 3)
        
    def test_eman1_math_linear(self):
        """test eman1.math.linear processor ................."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        d = e.get_3dview()
        
        e2 = e.copy()
        e2.process_inplace('eman1.math.linear', {'scale':3.23, 'shift':2.56})
        d2 = e2.get_3dview()
        for x in range(32):
            for y in range(32):
                for z in range(32):
                    self.assertAlmostEqual(d2[z][y][x], d[z][y][x]*3.23 + 2.56, 3)
        
    def test_eman1_math_exp(self):
        """test eman1.math.exp processor ...................."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        d = e.get_3dview()
        
        e2 = e.copy()
        e2.process_inplace('eman1.math.exp', {'low':2.56, 'high':3.23})
        d2 = e2.get_3dview()
        from math import exp
        for x in range(32):
            for y in range(32):
                for z in range(32):
                    self.assertAlmostEqual(d2[z][y][x], exp(d[z][y][x]/2.56 - 3.23), 3)
    
    def test_RangeThreshold(self):
        """test range threshhold processor .................."""
        imgfile1 = "test_RangeThreshold.mrc"
        TestUtil.make_image_file(imgfile1, MRC)
        e = EMData()
        e.read_image(imgfile1)
        low = 10
        high = 20
        f1 = Processors.get("eman1.threshold.binaryrange", {"low":low, "high":high});
        d1 = f1.get_params()
        self.assertEqual(d1["low"], low)
        self.assertEqual(d1["high"], high)

        #xc = 12
        rw = 12.5
        #f2 = Processors.get("eman1.mask.ringmean", {"xc":xc, "ring_width":rw})
        f2 = Processors.get("eman1.mask.ringmean", {"ring_width":rw})
        d2 = f2.get_params()
        #self.assertEqual(d2["xc"], xc)
        self.assertEqual(d2["ring_width"], rw)

        outfile1 = "test_RangeThreshold_out.mrc"
        
        e.process_inplace("eman1.threshold.binary", {"value": 200})
        e.write_image(outfile1)

        os.unlink(imgfile1)
        os.unlink(outfile1)

    def test_eman1_math_sigma(self):
        """test eman1.math.sigma processor .................."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        d = e.get_3dview()
        
        e2 = e.copy()
        mean = e2.get_attr('mean')
        sigma = e2.get_attr('sigma')
        v1 = 1
        v2 = 1
        e2.process_inplace('eman1.math.sigma', {'value1':v1, 'value2':v2})
        d2 = e2.get_3dview()
        for x in range(32):
            for y in range(32):
                for z in range(32):
                    if(d[z][y][x] < (mean - v2 * sigma) or d[z][y][x] > (mean + v1 * sigma)):
                        self.assertAlmostEqual(d2[z][y][x], mean, 3)
                    else:
                        self.assertAlmostEqual(d2[z][y][x], d[z][y][x], 3)
                        
    def test_eman1_math_log(self):
        """test eman1.math.log processor ...................."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        e -= 0.5
        d = e.get_3dview()
        
        e2 = e.copy()
        max = e2.get_attr('maximum')
        e2.process_inplace('eman1.math.log')
        d2 = e2.get_3dview()
        from math import log10
        for x in range(32):
            for y in range(32):
                for z in range(32):
                    if( d[z][y][x] > 0 ):
                        self.assertAlmostEqual(d2[z][y][x], log10(d[z][y][x]/max), 3)
                    else:
                        self.assertAlmostEqual(d2[z][y][x], 0.0, 3)
    
    def test_eman1_mask_sharp(self):
        """test eman1.mask.sharp processor .................."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e.process_inplace('eman1.mask.sharp', {'value':0.5})
        
    def test_eman1_mask_ringmean(self):
        """test eman1.mask.ringmean processor ..............."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e.process_inplace('eman1.mask.ringmean', {'ring_width':10})
        
    def test_eman1_mask_noise(self):
        """test eman1.mask.noise processor .................."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e.process_inplace('eman1.mask.noise')
        
    def test_eman1_mask_gaussian(self):
        """test eman1.mask.gaussian processor ..............."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e.process_inplace('eman1.mask.gaussian')
        
    def test_eman1_math_gausskernelfix(self):
        """test eman1.math.gausskernelfix processor ........."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e.process_inplace('eman1.math.gausskernelfix', {'gauss_width':10.0})
        
    def test_eman1_math_toradiussqr(self):
        """test eman1.math.toradiussqr processor ............"""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e.process_inplace('eman1.math.toradiussqr')
        
    def test_eman1_math_toradius(self):
        """test eman1.math.toradius processor ..............."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e.process_inplace('eman1.math.toradius')
        
    #need fix, this processor has problem    
    def no_test_eman1_complex_normpixels(self):
        """test eman1.complex.normpixels processor .........."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        e.do_fft_inplace()
        self.assertEqual(e.is_complex(), True)
        
        e0 = e.get_fft_phase()    #phase before nomalization
        d0 = e0.get_3dview()
        
        e.process_inplace('eman1.complex.normpixels')
        
        e2 = e.get_fft_amplitude()    #amplitude to 1.0
        d2 = e2.get_3dview()
        e3 = e.get_fft_phase()    #phase after normalization, unchanged
        d3 = e3.get_3dview()
        for x in range(32):
            for y in range(32):
                for z in range(32):
                    self.assertAlmostEqual(d2[z][y][x], 1.0, 3)
                    self.assertAlmostEqual(d0[z][y][x]. d3[z][y][x], 3)
                    
    def no_test_eman1_math_laplacian(self):
        """test eman1.math.laplacian processor .............."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e.process_inplace('eman1.math.laplacian')
        
    def test_eman1_mask_contract(self):
        """test eman1.mask.contract processor ..............."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e.process_inplace('eman1.mask.contract')
        
        e2 = EMData()
        e2.set_size(32,32,32)
        e2.to_zero()
        e2.process_inplace('eman1.mask.contract')
        d2 = e2.get_3dview()
        for x in range(32):
            for y in range(32):
                for z in range(32):
                    self.assertAlmostEqual(d2[z][y][x], 0, 3)
                    
    def test_eman1_filter_median(self):
        """test eman1.filter.median processor ..............."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e.process_inplace('eman1.filter.median')
        
    def test_eman1_math_localsigma(self):
        """test eman1.math.localsigma processor ............."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e.process_inplace('eman1.math.localsigma')
        
    def test_eman1_math_localmax(self):
        """test eman1.math.localmax processor ..............."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e.process_inplace('eman1.math.localmax')
        
    def test_eman1_math_submax(self):
        """test eman1.math.submax processor ................."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e.process_inplace('eman1.math.submax')
        
    def test_eman1_mask_onlypeaks(self):
        """test eman1.mask.onlypeaks processor .............."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e.process_inplace('eman1.mask.onlypeaks', {'npeaks':2})
        
    def no_test_eman1_filter_blockrange(self):
        """test eman1.filter.blockrange processor ..........."""
        e = EMData()
        e.set_size(32,32,1)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e.process_inplace('eman1.filter.blockrange', {'cal_half_width':0.2, 'fill_half_width':0.3})
        
        #3D image not supported by this processor
        e2 = EMData()
        e2.set_size(32,32,32)
        Log.logger().set_level(-1)    #no log message printed out
        self.assertRaises( RuntimeError, e2.process, 'eman1.filter.blockrange', {'cal_half_width':0.2, 'fill_half_width':0.3} )
        try:
            e2.process_inplace('eman1.filter.blockrange', {'value1':0.2, 'value2':0.3})
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "ImageDimensionException")
        
    def no_test_eman1_filter_blockcutoff(self):
        """test eman1.filter.blockcutoff processor .........."""
        e = EMData()
        e.set_size(32,32,1)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e.process_inplace('eman1.filter.blockcutoff', {'value1':0.2, 'value2':0.3})
        
        #3D image not supported by this processor
        e2 = EMData()
        e2.set_size(32,32,32)
        Log.logger().set_level(-1)    #no log message printed out
        self.assertRaises( RuntimeError, e2.process, 'eman1.filter.blockcutoff', {'value1':0.2, 'value2':0.3} )
        try:
            e2.process_inplace('eman1.filter.blockcutoff', {'value1':0.2, 'value2':0.3})
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "ImageDimensionException")
        
    def no_test_eman1_math_lineargradientfix(self):
        """test eman1.math.lineargradientfix processor ......"""
        e = EMData()
        e.set_size(32,32, 1)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e.process_inplace('eman1.math.lineargradientfix')
        
        #3D image not supported by this processor
        e2 = EMData()
        e2.set_size(32,32,32)
        Log.logger().set_level(-1)    #no log message printed out
        self.assertRaises( RuntimeError, e2.process, 'eman1.math.lineargradientfix' )
        try:
            e2.process_inplace('eman1.math.lineargradientfix')
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "ImageDimensionException")
        
    def test_filter_ramp(self):
        """test filter.ramp processor ......................."""
        e = EMData()
        e.set_size(32,32, 1)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e.process_inplace('filter.ramp')
        
        #3D image not supported by this processor
        e2 = EMData()
        e2.set_size(32,32,32)
        Log.logger().set_level(-1)    #no log message printed out
        self.assertRaises( RuntimeError, e2.process_inplace, 'filter.ramp' )
        try:
            e2.process_inplace('filter.ramp')
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "ImageDimensionException")
        
    def test_eman1_math_verticalstripefix(self):
        """test eman1.math.verticalstripefix processor ......"""
        e = EMData()
        e.set_size(32,32, 32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e.process_inplace('eman1.math.verticalstripefix')
        
    def test_eman1_math_realtofft(self):
        """test eman1.math.realtofft processor .............."""
        e = EMData()
        e.set_size(32,32,1)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e.process_inplace('eman1.math.realtofft')
        
        #3D image not supported by this processor
        e2 = EMData()
        e2.set_size(32,32,32)
        Log.logger().set_level(-1)    #no log message printed out
        self.assertRaises( RuntimeError, e2.process_inplace, 'eman1.math.realtofft' )
        try:
            e2.process_inplace('eman1.math.realtofft')
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "ImageDimensionException")
            
        #apply to real image only
        e3 = EMData()
        e3.set_size(32,32,1)
        e3.process_inplace('testimage.noise.uniform.rand')
        e3.do_fft_inplace()
        self.assertRaises( RuntimeError, e3.process_inplace, 'eman1.math.realtofft' )
        try:
            e3.process_inplace('eman1.math.realtofft')
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "ImageFormatException")
            
    def test_eman1_mask_zeroedgefill(self):
        """test eman1.mask.zeroedgefill processor ..........."""
        e = EMData()
        e.set_size(32,32,1)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e.process_inplace('eman1.mask.zeroedgefill')
        
        #3D image not supported by this processor
        e2 = EMData()
        e2.set_size(32,32,32)
        Log.logger().set_level(-1)    #no log message printed out
        self.assertRaises( RuntimeError, e2.process_inplace, 'eman1.mask.zeroedgefill' )
        try:
            e2.process_inplace('eman1.mask.zeroedgefill')
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "ImageDimensionException")
            
    def test_eman1_mask_beamstop(self):
        """test eman1.mask.beamstop processor ..............."""
        e = EMData()
        e.set_size(32,32,1)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e.process_inplace('eman1.mask.beamstop', {'value1':-1.2, 'value2':10, 'value3':8})
        
        #3D image not supported by this processor
        e2 = EMData()
        e2.set_size(32,32,32)
        Log.logger().set_level(-1)    #no log message printed out
        self.assertRaises( RuntimeError, e2.process_inplace, 'eman1.mask.beamstop', \
                            {'value1':-1.2, 'value2':10, 'value3':8} )
        try:
            e2.process_inplace('eman1.mask.beamstop', {'value1':-1.2, 'value2':10, 'value3':8})
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "ImageDimensionException")
            
    def test_eman1_mask_dampedzeroedgefill(self):
        """test eman1.mask.dampedzeroedgefill processor ....."""
        e = EMData()
        e.set_size(32,32,1)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e.process_inplace('eman1.mask.dampedzeroedgefill')
        
        #3D image not supported by this processor
        e2 = EMData()
        e2.set_size(32,32,32)
        Log.logger().set_level(-1)    #no log message printed out
        self.assertRaises( RuntimeError, e2.process_inplace, 'eman1.mask.dampedzeroedgefill' )
        try:
            e2.process_inplace('eman1.mask.dampedzeroedgefill')
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "ImageDimensionException")
            
    def test_eman1_math_averageovery(self):
        """test eman1.math.averageovery processor ..........."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e.process_inplace('eman1.math.averageovery')
        
    def test_eman1_mask_zeroedge2d(self):
        """test eman1.mask.zeroedge2d processor ............."""
        e = EMData()
        e.set_size(32,32,1)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e.process_inplace('eman1.mask.zeroedge2d', {'x0':1, 'y0':1, 'x1':25, 'y1':30})
        
        #3D image not supported by this processor
        e2 = EMData()
        e2.set_size(32,32,32)
        Log.logger().set_level(-1)    #no log message printed out
        self.assertRaises( RuntimeError, e2.process_inplace, 'eman1.mask.zeroedge2d', {'x0':1, 'y0':1, 'x1':25, 'y1':30})
        try:
            e2.process_inplace('eman1.mask.zeroedge2d', {'x0':1, 'y0':1, 'x1':25, 'y1':30})
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "ImageDimensionException")
            
    def test_eman1_mask_zeroedge3d(self):
        """test eman1.mask.zeroedge3d processor ............."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e.process_inplace('eman1.mask.zeroedge3d', {'x0':2, 'y0':3, 'z0':4, 'x1':20, 'y1':22, 'z1':26})
        
        #3D image only
        e2 = EMData()
        e2.set_size(32,32,1)
        Log.logger().set_level(-1)    #no log message printed out
        self.assertRaises( RuntimeError, e2.process_inplace, 'eman1.mask.zeroedge3d', {'x0':2, 'y0':3, 'z0':4, 'x1':20, 'y1':22, 'z1':26})
        try:
            e2.process_inplace('eman1.mask.zeroedge3d', {'x0':2, 'y0':3, 'z0':4, 'x1':20, 'y1':22, 'z1':26})
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "ImageDimensionException")
            
    def test_eman1_bilateral(self):
        """test eman1.bilateral processor ..................."""
        e = EMData()
        e.set_size(16,16,16)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e.process_inplace('eman1.bilateral', {'distance_sigma':0.3, 'value_sigma':0.4, 'niter':2, 'half_width':5})
        
        e2 = EMData()
        e2.set_size(32,32)
        e2.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e2.is_complex(), False)
        e.process_inplace('eman1.bilateral', {'distance_sigma':0.3, 'value_sigma':0.4, 'niter':2, 'half_width':5})
            
        
    def test_eman1_normalize_unitlen(self):
        """test eman1.normalize.unitlen processor ..........."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e.process_inplace('eman1.normalize.unitlen')
        
    def test_eman1_normalize_unitsum(self):
        """test eman1.normalize.unitsum processor ..........."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e.process_inplace('eman1.normalize.unitsum')
        
    def test_eman1_normalize(self):
        """test eman1.normalize processor ..................."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e.process_inplace('eman1.normalize')
        
    def test_eman1_normalize_mask(self):
        """test eman1.normalize.mask processor .............."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e2 = EMData()
        e2.set_size(10,10,10)
        e2.process_inplace('testimage.noise.uniform.rand')
        
        e.process_inplace('eman1.normalize.mask', {'mask':e2, 'no_sigma':3})
        
    def test_eman1_normalize_edgemean(self):
        """test eman1.normalize.edgemean processor .........."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e.process_inplace('eman1.normalize.edgemean')
        
    def test_eman1_normalize_circlemean(self):
        """test eman1.normalize.circlemean processor ........"""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e.process_inplace('eman1.normalize.circlemean')
        
    def test_eman1_normalize_lredge(self):
        """test eman1.normalize.lredge processor ............"""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e.process_inplace('eman1.normalize.lredge')
        
    def test_eman1_normalize_maxmin(self):
        """test eman1.normalize.maxmin processor ............"""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e.process_inplace('eman1.normalize.maxmin')
        
    def test_eman1_normalize_rows(self):
        """test eman1.normalize.rows processor .............."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e.process_inplace('eman1.normalize.rows')
        
    def test_eman1_normalize_toimage(self):
        """test eman1.normalize.toimage processor ..........."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e2 = EMData()
        e2.set_size(32,32,32)
        e2.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e2.is_complex(), False)
        
        e.process_inplace('eman1.normalize.toimage', {'noisy':e2, 'keepzero':2, 'invert':1, \
                                                'mult':2.3, 'add':0.5})
                                                
    def test_eman1_normalize_tofile(self):
        """test eman1.normalize.tofile processor ............"""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e2 = EMData()
        e2.set_size(32,32,32)
        e2.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e2.is_complex(), False)
        filename = 'noise.mrc'
        e2.write_image(filename)
        
        e.process_inplace('eman1.normalize.tofile', {'noisyfile':filename, 'keepzero':2, \
                   'invert':1, 'mult':2.3, 'add':0.5 })
                   
        os.unlink(filename)
        
    def test_eman1_normalize_toimage_lsq(self):
        """test eman1.normalize.toimage.lsq processor ......."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e2 = EMData()
        e2.set_size(32,32,32)
        e2.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e2.is_complex(), False)
        
        e.process_inplace('eman1.normalize.toimage.lsq', {'to':e2, 'low_threshold':0.2, 'high_threshold':0.8})
        
    def test_eman1_math_radialaverage(self):
        """test eman1.math.radialaverage processor .........."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e.process_inplace('eman1.math.radialaverage')
        
    def test_eman1_math_radialsubtract(self):
        """test eman1.math.radialsubtract processor ........."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e.process_inplace('eman1.math.radialsubtract')
        
    def test_eman1_xform_flip(self):
        """test eman1.xform.flip processor .................."""
        e = EMData()
        e.set_size(2,2,2)
        self.assertEqual(e.is_complex(), False)
        e.set_value_at(0,0,0, 1)
        e.set_value_at(0,0,1, 2)
        e.set_value_at(0,1,0, 3)
        e.set_value_at(0,1,1, 4)
        e.set_value_at(1,0,0, 5)
        e.set_value_at(1,0,1, 6)
        e.set_value_at(1,1,0, 7)
        e.set_value_at(1,1,1, 8)
        
        e.process_inplace('eman1.xform.flip', {'axis':'x'})
        self.assertEqual(e.get_value_at(0, 0, 0), 5)
        self.assertEqual(e.get_value_at(0, 1, 0), 7)
        self.assertEqual(e.get_value_at(1, 0, 0), 1)
        self.assertEqual(e.get_value_at(1, 1, 0), 3)
        
        e.process_inplace('eman1.xform.flip', {'axis':'y'})
        self.assertEqual(e.get_value_at(0, 0, 0), 7)
        self.assertEqual(e.get_value_at(0, 1, 0), 5)
        self.assertEqual(e.get_value_at(1, 0, 0), 3)
        self.assertEqual(e.get_value_at(1, 1, 0), 1)
        
        e.process_inplace('eman1.xform.flip', {'axis':'z'})
        self.assertEqual(e.get_value_at(0, 0, 0), 8)
        self.assertEqual(e.get_value_at(0, 1, 0), 6)
        self.assertEqual(e.get_value_at(1, 0, 0), 4)
        self.assertEqual(e.get_value_at(1, 1, 0), 2)
    
    def test_eman1_math_addnoise(self):
        """test eman1.math.addnoise processor ..............."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e.process_inplace('eman1.math.addnoise', {'noise':0.89})
        
    def test_eman1_math_addsignoise(self):
        """test eman1.math.addsignoise processor ............"""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e.process_inplace('eman1.math.addsignoise')
        
    def test_eman1_addspectralnoise(self):
        """test eman1.addspectralnoise processor ............"""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        e.do_fft_inplace()
        self.assertEqual(e.is_complex(), True)
        
        e.process_inplace('eman1.addspectralnoise', {'n':2, 'x0':9.8, 'dx':1.2, \
                                    'y':(1.3,2.4), 'interpolation':1})
        
        #complex image only
        e2 = EMData()
        e2.set_size(32,32,32)
        e2.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e2.is_complex(), False)
        Log.logger().set_level(-1)    #no log message printed out
        self.assertRaises( RuntimeError, e2.process_inplace, 'eman1.addspectralnoise', \
                            {'n':2, 'x0':9.8, 'dx':1.2, 'y':(1.3,2.4), 'interpolation':1})
        try:
            e2.process_inplace('eman1.addspectralnoise', {'n':2, 'x0':9.8, 'dx':1.2, \
                                    'y':(1.3,2.4), 'interpolation':1})
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "ImageFormatException")
    
    def test_eman1_xform_fourierorigin(self):
        """test eman1.xform.fourierorigin processor ........."""
        e = EMData()
        e.set_size(32,32,1)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        e.do_fft_inplace()
        self.assertEqual(e.is_complex(), True)
        
        e.process_inplace('eman1.xform.fourierorigin')
        
    def test_eman1_xform_phaseorigin(self):
        """test eman1.xform.phaseorigin processor ..........."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e.process_inplace('eman1.xform.phaseorigin')
        
    def test_eman1_mask_auto2d(self):
        """test eman1.mask.auto2d processor ................."""
        e = EMData()
        e.set_size(32,32,1)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e.process_inplace('eman1.mask.auto2d', {'threshold':0.5, 'filter':0.1})
        
        #2D image only
        e2 = EMData()
        e2.set_size(32,32,32)
        e2.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        self.assertRaises( RuntimeError, e2.process_inplace, \
            'eman1.mask.auto2d', {'threshold':0.5, 'filter':0.1})
        try:
            e2.process_inplace('eman1.mask.auto2d', {'threshold':0.5, 'filter':0.1})
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "ImageDimensionException")
    
    def no_test_eman1_mask_auto3d_thresh(self):
        """test eman1.mask.auto3d.thresh processor .........."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e.process_inplace('eman1.mask.auto3d.thresh', {'threshold1':0.23, 'threshold2':0.86})
        
    def test_eman1_mask_auto3d(self):
        """test eman1.mask.auto3d processor ................."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e.process_inplace('eman1.mask.auto3d', {'radius':16, 'threshold':0.5, 'nshells':3})
        
    def test_eman1_mask_addshells(self):
        """test eman1.mask.addshells processor .............."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e.process_inplace('eman1.mask.addshells', {'nshells':3})
        
    def test_eman1_xform_centerofmass(self):
        """test eman1.xform.centerofmass processor .........."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e.process_inplace('eman1.xform.centerofmass', {'int_shift_only':2})
        
    def test_eman1_xform_centeracf(self):
        """test eman1.xform.centeracf processor ............."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e.process_inplace('eman1.xform.centeracf', {'is3d':1})
        
    def no_test_eman1_filter_snr(self):
        """test eman1.filter.snr processor .................."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        #need supply this snr file, check it back later
        snr_file = 'snrfile.txt'
        e.process_inplace('eman1.filter.snr', {'wiener':3, 'snrfile':snr_file})
        
    def no_test_eman1_filter_byfile(self):
        """test eman1.filter.byfile processor ..............."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        filter_file = 'filter.txt'
        e.process_inplace('eman1.filter.byfile', {'filename':filter_file})
        
    def no_test_eman1_misc_symsearch(self):
        """test eman1.misc.symsearch processor .............."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e2 = EMData()
        e2.set_size(32,32,32)
        e2.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e2.is_complex(), False)
        
        e.process_inplace('eman1.misc.symsearch', {'sym':['CSYM'], 'thresh':0.5, 'output_symlabel':1, 'symlabel_map':e2})
        
    def test_eman1_misc_localnorm(self):
        """test eman1.misc.localnorm processor .............."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e.process_inplace('eman1.misc.localnorm', {'threshold':0.4, 'radius':16, 'apix':0.8})
        f = e.process('eman1.misc.localnorm', {'threshold':0.4, 'radius':16, 'apix':0.8})
        
        os.unlink('norm.mrc')
        
    def test_eman1_mask_fromfile(self):
        """test eman1.mask.fromfile processor ..............."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e2 = EMData()
        e2.set_size(32,32,32)
        e2.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e2.is_complex(), False)
        filename = 'maskfile.mrc'
        e2.write_image(filename)
        
        e.process_inplace('eman1.mask.fromfile', {'filename':filename, 'ismaskset':1})
        
        os.unlink(filename)
        
    def test_eman1_mask_fromfile_sizediff(self):
        """test eman1.mask.fromfile.sizediff processor ......"""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e2 = EMData()
        e2.set_size(64,48,32)
        e2.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e2.is_complex(), False)
        filename = 'maskfile.mrc'
        e2.write_image(filename)
        
        e.process_inplace('eman1.mask.fromfile.sizediff', {'filename':filename})
        
        os.unlink(filename)
        
    def no_test_eman1_misc_setpowspec(self):
        """test eman1.misc.setpowspec processor ............."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        filename = 'powerspec.txt'
        
        e.process_inplace('eman1.misc.setpowspec', {'filename':filename})
        
    def test_eman1_mask_smart(self):
        """test eman1.mask.smart processor .................."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e.process_inplace('eman1.mask.smart', {'mask':1.1})
        
    def test_eman1_mask_addshells_gauss(self):
        """test eman1.mask.addshells.gauss processor ........"""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e.process_inplace('eman1.mask.addshells.gauss')
        
    def test_testimage_puregaussian(self):
        """test testimage.puregaussian processor ............"""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.puregaussian', {'sigma':12})
        
    def test_testimage_gaussian(self):
        """test testimage.gaussian processor ................"""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.puregaussian', {'sigma':12})
        
    def test_testimage_scurve(self):
        """test testimage.scurve processor .................."""
        e = EMData()
        e.set_size(32,32,1)
        e.process_inplace('testimage.scurve')
        
        #only for 2D image
        e2 = EMData()
        e2.set_size(32,32,32)
        self.assertRaises( RuntimeError, e2.process_inplace, 'testimage.scurve')
        try:
            e2.process_inplace('testimage.scurve')
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "ImageDimensionException")
    
    def test_testimage_sinewave(self):
        """test testimage.sinewave processor ................"""
        e = EMData()
        e.set_size(32,1,1)
        e.process_inplace('testimage.sinewave', {'wave_length':10})
        e.process_inplace('testimage.sinewave', {'wave_length':10, 'phase':90})
        
        e2 = EMData()
        e2.set_size(32,32,1)
        e2.process_inplace('testimage.sinewave', {'wave_length':10})
        e2.process_inplace('testimage.sinewave', {'wave_length':10, 'axis':'y', 'phase':90, 'az':45})
        
        e3 = EMData()
        e3.set_size(32,32,32)
        e3.process_inplace('testimage.sinewave', {'wave_length':10})
        e3.process_inplace('testimage.sinewave', {'wave_length':10, 'axis':'z', 'phase':90, 'az':30, 'alt':45, 'phi':60})
            
    def test_testimage_sinewavecircular(self):
        """test testimage.sinewave.circular processor ......."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.sinewave.circular', {'wave_length':10})
        e.process_inplace('testimage.sinewave.circular', {'wave_length':10, 'axis':'z', 'c':7, 'phase':3})
        
    def test_testimage_squarecube(self):
        """test testimage.squarecube processor .............."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.squarecube', {'edge_length':22})
        e.process_inplace('testimage.squarecube', {'edge_length':22, 'axis':'y', 'odd_edge':18, 'fill':'yes'})
        
    def test_testimage_circlesphere(self):
        """test testimage.circlesphere processor ............"""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.circlesphere', {'radius':20})
        e.process_inplace('testimage.circlesphere', {'radius':20, 'axis':'z', 'c':15, 'fill':'yes'})
        
    def test_testimage_noise_uniform_rand(self):
        """test testimage.noise.uniform.rand processor ......"""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        
    def test_testimage_noise_gauss(self):
        """test testimage.noise.gauss processor ............."""
        e1 = EMData()
        e1.set_size(8,8,8)
        e1.process_inplace('testimage.noise.gauss', {'mean':0.3, 'sigma':1.1, 'seed':6666})
        d1 = e1.get_3dview()
        
        e2 = EMData()
        e2.set_size(8, 8, 8)
        e2.process_inplace('testimage.noise.gauss', {'mean':0.3, 'sigma':1.1})
        d2 = e2.get_3dview()
        
        e3 = EMData()
        e3.set_size(8,8,8)
        e3.process_inplace('testimage.noise.gauss', {'mean':0.3, 'sigma':1.1, 'seed':6666})
        d3 = e1.get_3dview()
        
        e4 = EMData()
        e4.set_size(8, 8, 8)
        e4.process_inplace('testimage.noise.gauss', {'mean':0.3, 'sigma':1.1})
        d4 = e2.get_3dview()
        
        for k in range(8):
            for j in range(8):
                for i in range(8):
                    self.assertAlmostEqual(d1[i][j][k], d3[i][j][k])
                    self.assertAlmostEqual(d2[i][j][k], d4[i][j][k])
        
    #this filter.integercyclicshift2d processor is removed by Phani at 5/18/2006    
    def no_test_IntegerCyclicShift2DProcessor(self):
        """test filter.integercyclicshift2d processor........"""
        e = EMData()
        e.set_size(50,50,1)
        e.process_inplace('testimage.noise.uniform.rand')
        
        d1 = e.get_2dview()
        d3 = Numeric.array(d1)    #make a copy of d1, since d1 and d2 share the same memory
        e.process_inplace('filter.integercyclicshift2d', {'dx':10, 'dy':20})
        d2 = e.get_2dview()
        for x in range(50):
            if x+10 > 49:
                x2 = x+10-50
            else:
                x2 = x+10
            for y in range(50):
                if y+20 > 49:
                    y2 = y+20-50
                else:
                    y2 = y+20
                self.assertEqual(d3[y][x], d2[y2][x2])
                
        #this filter apply to 2D real image only
        e2 = e.do_fft()
        self.assertEqual(e2.is_complex(), True)
        self.assertRaises( RuntimeError, e2.process_inplace, 'filter.integercyclicshift2d', {'dx':10, 'dy':20})
        try:
            e2.process_inplace('filter.integercyclicshift2d', {'dx':10, 'dy':20})
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "ImageFormatException")
        
        e3 = EMData()
        e3.set_size(1,50,50)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertRaises( RuntimeError, e3.process_inplace, 'filter.integercyclicshift2d', {'dx':10, 'dy':20})
        try:
            e3.process_inplace('filter.integercyclicshift2d', {'dx':10, 'dy':20})
        except RuntimeError, runtime_err:
            self.assertEqual(exception_type(runtime_err), "ImageFormatException")
        
def test_main():
    test_support.run_unittest(TestProcessor)

if __name__ == '__main__':
    test_main()


