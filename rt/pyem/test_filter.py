#!/bin/env python

from EMAN2 import *
import unittest,os,sys
from test import test_support
import testlib
from pyemtbx.exceptions import *

class TestProcessor(unittest.TestCase):
    """Processor test"""
    
    def test_get_processor_list(self):
        """test get processor list .........................."""
        processor_names = Processors.get_list()
        self.assertEqual(len(processor_names), 107)

        try:
            f2 = Processors.get("_nosuchfilter___")
        except RuntimeError, runtime_err:
            err_type = exception_type(runtime_err)
            self.assertEqual(err_type, "NotExistingObjectException")
            
    def test_filter_lowpass_tophat(self):
        """test filter.lowpass.tophat processor ............."""
        e = EMData()
        e.set_size(32,32,32)
        e.process('testimage.noise.uniform.rand')
        e.do_fft_inplace()
        self.assertEqual(e.is_complex(), True)
        
        e.process('filter.highpass.tophat', {'Cutoff_frequency':0.3})
        
    def test_filter_highpass_tophat(self):
        """test filter.lowpass.tophat processor ............."""
        e = EMData()
        e.set_size(32,32,32)
        e.process('testimage.noise.uniform.rand')
        e.do_fft_inplace()
        self.assertEqual(e.is_complex(), True)
        
        e.process('filter.lowpass.tophat', {'Cutoff_frequency':0.3})
        
    def test_filter_bandpass_tophat(self):
        """test filter.bandpass.tophat processor ............"""
        e = EMData()
        e.set_size(32,32,32)
        e.process('testimage.noise.uniform.rand')
        e.do_fft_inplace()
        self.assertEqual(e.is_complex(), True)
        
        e.process('filter.bandpass.tophat', \
            {'Low_cutoff_frequency':0.1, 'High_cutoff_frequency':0.45})
            
    def test_filter_homomorphic_tophat(self):
        """test filter.homomorphic.tophat processor ........."""
        e = EMData()
        e.set_size(32,32,32)
        e.process('testimage.noise.uniform.rand')
        e.do_fft_inplace()
        self.assertEqual(e.is_complex(), True)
        
        e.process('filter.homomorphic.tophat', \
            {'Low_cutoff_frequency':0.05, 'High_cutoff_frequency':0.45, 'Value_at_zero_frequency':1.0})

    def test_filter_lowpass_gauss(self):
        """test filter.lowpass.gauss processor .............."""
        e = EMData()
        e.set_size(32,32,32)
        e.process('testimage.noise.uniform.rand')
        e.do_fft_inplace()
        self.assertEqual(e.is_complex(), True)
        
        e.process('filter.lowpass.gauss', {'Sigma':0.5})
        
    def test_filter_highpass_gauss(self):
        """test filter.highpass.gauss processor ............."""
        e = EMData()
        e.set_size(32,32,32)
        e.process('testimage.noise.uniform.rand')
        e.do_fft_inplace()
        self.assertEqual(e.is_complex(), True)
        
        e.process('filter.highpass.gauss', {'Sigma':0.5})
        
    def test_filter_bandpass_gauss(self):
        """test filter.bandpass.gauss processor ............."""
        e = EMData()
        e.set_size(32,32,32)
        e.process('testimage.noise.uniform.rand')
        e.do_fft_inplace()
        self.assertEqual(e.is_complex(), True)
        
        e.process('filter.bandpass.gauss', {'Sigma':0.5, 'Center':0.1})
        
    def test_filter_homomorphic_gauss(self):
        """test filter.homomorphic.gauss processor .........."""
        e = EMData()
        e.set_size(32,32,32)
        e.process('testimage.noise.uniform.rand')
        e.do_fft_inplace()
        self.assertEqual(e.is_complex(), True)
        
        e.process('filter.homomorphic.gauss', {'Sigma':0.5, 'Value_at_zero_frequency':1.0})
        
    def test_filter_gaussinverse(self):
        """test filter.gaussinverse processor ..............."""
        e = EMData()
        e.set_size(32,32,32)
        e.process('testimage.noise.uniform.rand')
        e.do_fft_inplace()
        self.assertEqual(e.is_complex(), True)
        
        e.process('filter.gaussinverse', {'Sigma':0.5})
        
    def test_filter_kaiserI0inverse(self):
        """test filter.kaiserI0inverse processor ............"""
        e = EMData()
        e.set_size(32,32,32)
        e.process('testimage.noise.uniform.rand')
        e.do_fft_inplace()
        self.assertEqual(e.is_complex(), True)
        
        e.process('filter.kaiserI0inverse')
        
    def test_filter_kaisersinhinverse(self):
        """test filter.kaisersinhinverse processor .........."""
        e = EMData()
        e.set_size(32,32,32)
        e.process('testimage.noise.uniform.rand')
        e.do_fft_inplace()
        self.assertEqual(e.is_complex(), True)
        
        e.process('filter.kaisersinhinverse')
        
    def test_filter_radialtable(self):
        """test filter.radialtable processor ................"""
        e = EMData()
        e.set_size(32,32,32)
        e.process('testimage.noise.uniform.rand')
        e.do_fft_inplace()
        self.assertEqual(e.is_complex(), True)
        
        e.process('filter.radialtable', {'Table':(0.2,0.2,0.3)})
        
    def test_filter_lowpass_butterworth(self):
        """test filter.lowpass.butterworth processor ........"""
        e = EMData()
        e.set_size(32,32,32)
        e.process('testimage.noise.uniform.rand')
        e.do_fft_inplace()
        self.assertEqual(e.is_complex(), True)
        
        e.process('filter.lowpass.butterworth', \
            {'Low_cutoff_frequency':0.1, 'High_cutoff_frequency':0.45})
            
    def test_filter_highpass_butterworth(self):
        """test filter.highpass.butterworth processor ......."""
        e = EMData()
        e.set_size(32,32,32)
        e.process('testimage.noise.uniform.rand')
        e.do_fft_inplace()
        self.assertEqual(e.is_complex(), True)
        
        e.process('filter.highpass.butterworth', \
            {'Low_cutoff_frequency':0.1, 'High_cutoff_frequency':0.45})
            
    def test_filter_homomorphic_butterworth(self):
        """test filter.homomorphic.butterworth processor ...."""
        e = EMData()
        e.set_size(32,32,32)
        e.process('testimage.noise.uniform.rand')
        e.do_fft_inplace()
        self.assertEqual(e.is_complex(), True)
        
        e.process('filter.homomorphic.butterworth', \
            {'Low_cutoff_frequency':0.1, 'High_cutoff_frequency':0.45, 'Value_at_zero_frequency':1.0})
            
    def test_filter_lowpass_tanh(self):
        """test filter.lowpass.tanh processor ..............."""
        e = EMData()
        e.set_size(32,32,32)
        e.process('testimage.noise.uniform.rand')
        e.do_fft_inplace()
        self.assertEqual(e.is_complex(), True)
        
        e.process('filter.lowpass.tanh', {'Cutoff_frequency':0.25, 'Fall_off':0.1})
        
    def test_filter_highpass_tanh(self):
        """test filter.highpass.tanh processor .............."""
        e = EMData()
        e.set_size(32,32,32)
        e.process('testimage.noise.uniform.rand')
        e.do_fft_inplace()
        self.assertEqual(e.is_complex(), True)
        
        e.process('filter.highpass.tanh', {'Cutoff_frequency':0.25, 'Fall_off':0.1})
        
    def test_filter_homomorphic_tanh(self):
        """test filter.homomorphic.tanh processor ..........."""
        e = EMData()
        e.set_size(32,32,32)
        e.process('testimage.noise.uniform.rand')
        e.do_fft_inplace()
        self.assertEqual(e.is_complex(), True)
        
        e.process('filter.homomorphic.tanh', \
            {'Cutoff_frequency':0.25, 'Fall_off':0.1, 'Value_at_zero_frequency':1.0})
            
    def test_filter_bandpass_tanh(self):
        """test filter.bandpass.tanh processor .............."""
        e = EMData()
        e.set_size(32,32,32)
        e.process('testimage.noise.uniform.rand')
        e.do_fft_inplace()
        self.assertEqual(e.is_complex(), True)
        
        e.process('filter.bandpass.tanh', \
            {'Low_cutoff_frequency':0.1, 'Low_fall_off':0.15, 'High_cutoff_frequency':0.45, \
             'High_fall_off':0.15, 'Fall_off':1.25})
             
    def test_eman1_filter_lowpass_sharp(self):
        """test eman1.filter.lowpass.sharp processor ........"""
        e = EMData()
        e.set_size(32,32,32)
        e.process('testimage.noise.uniform.rand')
        e.do_fft_inplace()
        self.assertEqual(e.is_complex(), True)
        
        e.process('eman1.filter.lowpass.sharp', {'lowpass':0.1})
        
    def test_eman1_filter_highpass_sharp(self):
        """test eman1.filter.highpass.sharp processor ......."""
        e = EMData()
        e.set_size(32,32,32)
        e.process('testimage.noise.uniform.rand')
        e.do_fft_inplace()
        self.assertEqual(e.is_complex(), True)
        
        e.process('eman1.filter.highpass.sharp', {'highpass':0.45})
        
    def test_eman1_filter_lowpass_gaussian(self):
        """test eman1.filter.lowpass.gaussian processor ....."""
        e = EMData()
        e.set_size(32,32,32)
        e.process('testimage.noise.uniform.rand')
        e.do_fft_inplace()
        self.assertEqual(e.is_complex(), True)
        
        e.process('eman1.filter.lowpass.gaussian', {'lowpass':0.1})
        
    def test_eman1_filter_highpass_gaussian(self):
        """test eman1.filter.highpass.gaussian processor ...."""
        e = EMData()
        e.set_size(32,32,32)
        e.process('testimage.noise.uniform.rand')
        e.do_fft_inplace()
        self.assertEqual(e.is_complex(), True)
        
        e.process('eman1.filter.highpass.gaussian', {'highpass':0.45})
        
    def test_eman1_filter_lowpass_tanh(self):
        """test eman1.filter.lowpass.tanh processor ........."""
        e = EMData()
        e.set_size(32,32,32)
        e.process('testimage.noise.uniform.rand')
        e.do_fft_inplace()
        self.assertEqual(e.is_complex(), True)
        
        e.process('eman1.filter.lowpass.tanh', {'lowpass':0.1})
        
    def test_eman1_filter_highpass_tanh(self):
        """test eman1.filter.highpass.tanh processor ........"""
        e = EMData()
        e.set_size(32,32,32)
        e.process('testimage.noise.uniform.rand')
        e.do_fft_inplace()
        self.assertEqual(e.is_complex(), True)
        
        e.process('eman1.filter.highpass.tanh', {'highpass':0.45})
        
    def test_eman1_filter_highpass_butterworth(self):
        """test eman1.filter.highpass.butterworth processor ."""
        e = EMData()
        e.set_size(32,32,32)
        e.process('testimage.noise.uniform.rand')
        e.do_fft_inplace()
        self.assertEqual(e.is_complex(), True)
        
        e.process('eman1.filter.highpass.butterworth', {'highpass':0.45})
        
    def test_eman1_filter_ramp(self):
        """test eman1.filter.ramp processor ................."""
        e = EMData()
        e.set_size(32,32,32)
        e.process('testimage.noise.uniform.rand')
        e.do_fft_inplace()
        self.assertEqual(e.is_complex(), True)
        
        e.process('eman1.filter.ramp', {'intercept':0.25, 'slope':0.3})
        
    def test_eman1_math_absvalue(self):
        """test eman1.math.absvalue processor ..............."""
        e = EMData()
        e.set_size(32,32,32)
        e.process('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        d = e.get_3dview()
        e2 = e * -1
        d2 = e2.get_3dview()
        for x in range(32):
            for y in range(32):
                for z in range(32):
                    self.assertAlmostEqual(d[z][y][x], -d2[z][y][x], 3)
        
        e2.process('eman1.math.absvalue')
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
        e.process('eman1.threshold.notzero')
        d = e.get_3dview()
        for x in range(32):
            for y in range(32):
                for z in range(32):
                    self.assertEqual(d[z][y][x], 0)
        
        e.to_one()
        e *= 0.5
        e.process('eman1.threshold.notzero')
        d = e.get_3dview()
        for x in range(32):
            for y in range(32):
                for z in range(32):
                    self.assertEqual(d[z][y][x], 1)
   
    def test_eman1_math_squared(self):
        """test eman1.math.squared processor ................"""
        e = EMData()
        e.set_size(32,32,32)
        e.process('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        d = e.get_3dview()
        e2 = e.copy()
        e2.process('eman1.math.squared')
        d2 = e2.get_3dview()
        for x in range(32):
            for y in range(32):
                for z in range(32):
                    self.assertAlmostEqual(d[z][y][x]**2, d2[z][y][x], 3)
                    
    def test_eman1_math_sqrt(self):
        """test eman1.math.sqrt processor ..................."""
        e = EMData()
        e.set_size(32,32,32)
        e.process('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        d = e.get_3dview()
        e2 = e.copy()
        e2.process('eman1.math.sqrt')
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
        
        e.process('eman1.threshold.belowtozero', {'minval':0.5})
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
        f1.process(e)
        testlib.check_emdata(e, sys.argv[0])
        os.unlink(imgfile1)

    def test_eman1_threshold_compress(self):
        """test eman1.threshold.compress processor .........."""
        e = EMData()
        e.set_size(32,32,32)
        e.to_one()
        e *= 0.45
        e.process('eman1.threshold.compress', {'range':0.2, 'value':0.5})
        d = e.get_3dview()
        for x in range(32):
            for y in range(32):
                for z in range(32):
                    self.assertAlmostEqual(d[z][y][x], 0.5, 3)
        
        e2 = EMData()
        e2.set_size(32,32,32)
        e2.to_one()
        e2 *= 0.05
        e2.process('eman1.threshold.compress', {'range':0.2, 'value':0.5})
        d2 = e2.get_3dview()
        for x in range(32):
            for y in range(32):
                for z in range(32):
                    self.assertAlmostEqual(d2[z][y][x], 0.2+0.05, 3)
        
    def test_eman1_math_linear(self):
        """test eman1.math.linear processor ................."""
        e = EMData()
        e.set_size(32,32,32)
        e.process('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        d = e.get_3dview()
        
        e2 = e.copy()
        e2.process('eman1.math.linear', {'scale':3.23, 'shift':2.56})
        d2 = e2.get_3dview()
        for x in range(32):
            for y in range(32):
                for z in range(32):
                    self.assertAlmostEqual(d2[z][y][x], d[z][y][x]*3.23 + 2.56, 3)
        
    def test_eman1_math_exp(self):
        """test eman1.math.exp processor ...................."""
        e = EMData()
        e.set_size(32,32,32)
        e.process('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        d = e.get_3dview()
        
        e2 = e.copy()
        e2.process('eman1.math.exp', {'low':2.56, 'high':3.23})
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

        xc = 12
        rw = 12.5
        f2 = Processors.get("eman1.mask.ringmean", {"xc":xc, "ring_width":rw})
        d2 = f2.get_params()
        self.assertEqual(d2["xc"], xc)
        self.assertEqual(d2["ring_width"], rw)

        outfile1 = "test_RangeThreshold_out.mrc"
        
        e.process("eman1.threshold.binary", {"value": 200})
        e.write_image(outfile1)

        os.unlink(imgfile1)
        os.unlink(outfile1)


        
class TestCmp(unittest.TestCase):
    """cmp test"""
    
    def test_variance(self):
        """test variance ...................................."""
        imgfile1 = "test_variance_1.mrc"
        TestUtil.make_image_file(imgfile1, MRC)
        e1 = EMData()
        e1.read_image(imgfile1)

        e2 = e1.copy()
        score = e2.cmp("variance", e1, {"keepzero": 0})
        self.assertEqual(score, 0)
        os.unlink(imgfile1)
        
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

        dot_score = e2.cmp("dot", e1, {"evenonly":0})
#        self.assertEqual(dot_score, 19944.0)    #todo: dot score not match, anything wrong?

        variance_score = e2.cmp("variance", e1, {"keepzero":1})
#        self.assertEqual(variance_score, 0)    #todo: score not match, anything wrong?
        
        phase_score = e2.cmp("phase", e1, {})
#        testlib.assertfloat(self, phase_score, 1.6488)    #todo: score not match, anything wrong?
        
        frc_score = e2.cmp("frc", e1, {})
#        testlib.assertfloat(self, frc_score, -0.4011)    #todo: score not match, anything wrong?

        (hed1,img1) = testlib.get_imagic_filename_pair(imgfile1)
        os.unlink(hed1)
        os.unlink(img1)

def test_main():
    test_support.run_unittest(TestProcessor, TestCmp)

if __name__ == '__main__':
    test_main()


