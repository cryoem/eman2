#!/usr/bin/env python

#
# Author: Grant Tang, 08/28/2003 (gtang@bcm.edu)
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

IS_TEST_EXCEPTION = False

class TestProcessor(unittest.TestCase):
    """Processor test"""
    
    def test_get_processor_list(self):
        """test get processor list .........................."""
        processor_names = Processors.get_list()
        self.assertEqual(len(processor_names), 170)
        
        if(IS_TEST_EXCEPTION):
            try:
                f2 = Processors.get("_nosuchfilter___")
            except RuntimeError, runtime_err:
                err_type = exception_type(runtime_err)
                self.assertEqual(err_type, "NotExistingObjectException")
                
    def test_transpose(self):
        """test xform.transpose processor ..................."""
        
        for y in [15,16]:
        	for x in [15,16]:
        		e = EMData(x,y,1)
                f = e.process('xform.transpose')
                f.process_inplace('xform.transpose')
                self.assertEqual(f==e,True)
                
                g = e.copy()
                g.process_inplace('xform.transpose')
                g.process_inplace('xform.transpose')
                self.assertEqual(g==e,True)
                
        if(IS_TEST_EXCEPTION):
        	# Check to make we throw for 3D and complex images
			e1 = EMData(2,2,2)
			e2 = EMData(2,2)
			e2.set_complex(True)
			for e in [e1,e2]:
				try:
				    e.process_inplace('xform.transpose')
				except RuntimeError, runtime_err:
					self.assertEqual(exception_type(runtime_err), "UnexpectedBehaviorException")
					
				try:
				    f =  e.process('xform.transpose')
				except RuntimeError, runtime_err:
					self.assertEqual(exception_type(runtime_err), "UnexpectedBehaviorException")
					
    def test_threshold_binary_fourier(self):
		"""test threshold.binary.fourier  ..................."""
		a = [test_image(0,(16,16)),test_image_3d(0,(16,16,16))]
		for e in a:
			af = e.do_fft()
			b = af.process("threshold.binary.fourier",{"value":0})
			self.assertAlmostEqual(b["mean"], 0.5, 8)
				
		if(IS_TEST_EXCEPTION):
			try:
				b = a[0].process("threshold.binary.fourier",{"value":0})
			except RuntimeError, runtime_err: # doesn't work on real images
				self.assertEqual(exception_type(runtime_err), "ImageFormatException")
				
			try:
				b = a[0].process("threshold.binary.fourier",{"value":-1})
			except RuntimeError, runtime_err: # doesn't work for negative thresholds
				self.assertEqual(exception_type(runtime_err), "InvalidParameterException")

    def test_flattenbackground(self):
        """test filter.flattenbackground processor .........."""
        e = EMData(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        e.process_inplace('filter.flattenbackground',{"radius":8})
		
        mask = EMData(16,16,16)
        mask.to_one()
        e.process_inplace('filter.flattenbackground',{"mask":mask})
		
        if(IS_TEST_EXCEPTION):
            # Check to make sure that when both the parameters are specified we throw
            # an InvalidParameterException
            try:
                e.process_inplace('filter.flattenbackground', {"radius":8,"mask":mask})
            except RuntimeError, runtime_err:
    			self.assertEqual(exception_type(runtime_err), "InvalidParameterException")
    		
    		# Check to make sure that atleast one parameter is specified
            try:
    			e.process_inplace('filter.flattenbackground', {})
            except RuntimeError, runtime_err:
    			self.assertEqual(exception_type(runtime_err), "InvalidParameterException")
    			
    		# make sure that we throw if the mask is to big
            mask = EMData(33,32,32)
            mask.to_one()
            try:
    			e.process_inplace('filter.flattenbackground',{"mask":mask})
            except RuntimeError, runtime_err:
    			self.assertEqual(exception_type(runtime_err), "ImageDimensionException")
		
    def test_filter_lowpass_tophat(self):
        """test filter.lowpass.tophat processor ............."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        e.do_fft_inplace()
        self.assertEqual(e.is_complex(), True)
        
        e.process_inplace('filter.highpass.tophat', {'cutoff_abs':0.3})
        
    def test_filter_highpass_tophat(self):
        """test filter.lowpass.tophat processor ............."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        e.do_fft_inplace()
        self.assertEqual(e.is_complex(), True)
        
        e.process_inplace('filter.lowpass.tophat', {'cutoff_abs':0.3})
        
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
        
        e.process_inplace('filter.lowpass.gauss', {'cutoff_abs':0.5})
        
    def test_filter_highpass_gauss(self):
        """test filter.highpass.gauss processor ............."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        e.do_fft_inplace()
        self.assertEqual(e.is_complex(), True)
        
        e.process_inplace('filter.highpass.gauss', {'cutoff_abs':0.5})
        
    def test_filter_bandpass_gauss(self):
        """test filter.bandpass.gauss processor ............."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        e.do_fft_inplace()
        self.assertEqual(e.is_complex(), True)
        
        e.process_inplace('filter.bandpass.gauss', {'cutoff_abs':0.5, 'center':0.1})
        
    def test_filter_homomorphic_gauss(self):
        """test filter.homomorphic.gauss processor .........."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        e.do_fft_inplace()
        self.assertEqual(e.is_complex(), True)
        
        e.process_inplace('filter.homomorphic.gauss', {'cutoff_abs':0.5, 'value_at_zero_frequency':1.0})
        
    def test_filter_gaussinverse(self):
        """test filter.gaussinverse processor ..............."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        e.do_fft_inplace()
        self.assertEqual(e.is_complex(), True)
        
        e.process_inplace('filter.gaussinverse', {'cutoff_abs':0.5})
        
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
        
    def no_test_filter_kaisersinhinverse(self):
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
        self.assertEqual(e.is_complex(), False)
        
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
        
        e.process_inplace('filter.lowpass.tanh', {'cutoff_abs':0.25, 'fall_off':0.1})
        
    def test_filter_highpass_tanh(self):
        """test filter.highpass.tanh processor .............."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        e.do_fft_inplace()
        self.assertEqual(e.is_complex(), True)
        
        e.process_inplace('filter.highpass.tanh', {'cutoff_abs':0.25, 'fall_off':0.1})
        
    def test_filter_homomorphic_tanh(self):
        """test filter.homomorphic.tanh processor ..........."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        e.do_fft_inplace()
        self.assertEqual(e.is_complex(), True)
        
        e.process_inplace('filter.homomorphic.tanh', \
            {'cutoff_abs':0.25, 'fall_off':0.1, 'value_at_zero_frequency':1.0})
            
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
    
    def test_mean_shrink(self):
        """test math.meanshrink processor ..................."""
        e = EMData()
        e.set_size(64,64,64)
        e.process_inplace("testimage.noise.uniform.rand")
        
        e.process_inplace("math.meanshrink",{"n":4})
        self.assertEqual(e.get_xsize(), 16)
        self.assertEqual(e.get_ysize(), 16)
        self.assertEqual(e.get_zsize(), 16)
        
        e2 = EMData()
        e2.set_size(30,30,1)
        e2.process_inplace("testimage.noise.uniform.rand")
        e2.process_inplace("math.meanshrink",{"n":1.5}) 
        
        if(IS_TEST_EXCEPTION):
            ##shrink factor 1.5 only support 2D image
            #self.assertRaises( InvalidValueException, e.process_inplace, 'math.meanshrink', {"n":1.5} )
            try:
               e.process_inplace("math.meanshrink",{"n":1.5})
            except RuntimeError, runtime_err:
                self.assertEqual(exception_type(runtime_err), "InvalidValueException")
                
            #shrink factor must be >1
            try:
                e.process_inplace("math.meanshrink",{"n":0})
            except RuntimeError, runtime_err:
                self.assertEqual(exception_type(runtime_err), "InvalidValueException")

    def test_median_shrink(self):
        """test math.medianshrink processor ................."""
        e = EMData()
        e.set_size(64,64,64)
        e.process_inplace("testimage.noise.uniform.rand")
        
        e.process_inplace("math.medianshrink",{"n":4})
        self.assertEqual(e.get_xsize(), 16)
        self.assertEqual(e.get_ysize(), 16)
        self.assertEqual(e.get_zsize(), 16)
        
        if(IS_TEST_EXCEPTION):
            #shrink factor must be >1 
            try:
                e.process_inplace("math.medianshrink",{"n":0})
            except RuntimeError, runtime_err:
                self.assertEqual(exception_type(runtime_err), "InvalidValueException")
            
            #image size must be divisible by shrink factor
            try:
                e.process_inplace("math.medianshrink",{"n":5})
            except RuntimeError, runtime_err:
                self.assertEqual(exception_type(runtime_err), "InvalidValueException")

    def test_eman1_math_absvalue(self):
        """test math.absvalue processor ....................."""
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
        
        e2.process_inplace('math.absvalue')
        d2 = e2.get_3dview()
        for x in range(32):
            for y in range(32):
                for z in range(32):
                    self.assertAlmostEqual(d[z][y][x], d2[z][y][x], 3)

    def test_threshold_clampminmax(self):
        """test threshold.clampminmax processor ............."""
        n = 32
        
        for i in [1,n]:
			e = EMData()
			e.set_size(n,n,i)
			e.process_inplace("testimage.noise.uniform.rand")
			
			cmax = e.get_attr("maximum")
			cmin = e.get_attr("minimum")
			
			nmax = cmax/2.0
			nmin = cmin/2.0
			
			a = {}
			a["minval"] = nmin
			a["maxval"] = nmax
			
			e.process_inplace("threshold.clampminmax", a)
			
			d = e.get_3dview()
			for z in range(i):
				for y in range(n):
					for x in range(n):
						assert d[z][y][x] >= nmin
						assert d[z][y][x] <= nmax
						
		
    def test_histogram_bin(self):
        """test historgram.bin processor ...................."""
        n = 32
        
        for i in [1,n]:
			for nbin in [1,128,256,512]:
				e = EMData()
				e.set_size(n,n,i)
				e.process_inplace("testimage.noise.uniform.rand")
				
				cmax = e.get_attr("mean")
				cmin = e.get_attr("mean")
				
				a = {}
				a["nbins"] = nbin
				#a["debug"] = 1
				e.process_inplace("histogram.bin", a)
				
				d = e.get_3dview()
				
				bins = []
				for z in range(i):
					for y in range(n):
						for x in range(n):
							streuth = False
							val = d[z][y][x] 
							for blah in bins:
								if ( val == blah ):
									streuth = True
									break
							
							if ( streuth == False ):
								bins.append(val)
				assert len(bins) <= nbin
						
    def test_threshold_clampminmax_nsigma(self):
        """test threshold.clampminmax.nsigma processor ......"""
        n = 32
        
        for i in [1,n]:
			for nsigma in [0.5,1,2]:
				e = EMData()
				e.set_size(n,n,i)
				e.process_inplace("testimage.noise.uniform.rand")
				
				cmax = e.get_attr("mean") + nsigma*e.get_attr("sigma")
				cmin = e.get_attr("mean") - nsigma*e.get_attr("sigma")
				
				a = {}
				a["nsigma"] = nsigma
				e.process_inplace("threshold.clampminmax.nsigma", a)
				
				d = e.get_3dview()
				for z in range(i):
					for y in range(n):
						for x in range(n):
							# the multiplication above in calculating cmax and cmin
							# sometimes causes some floating point precision issues
							# so the following approach overlooks cases where
							# the clamped values are slightly beyond the intended
							# clamping region
							if ( d[z][y][x] < cmin ):
								self.assertAlmostEqual(d[z][y][x] - cmin, 0, 6)
							if ( d[z][y][x] > cmax ):
								self.assertAlmostEqual(d[z][y][x] - cmax, 0, 6)
							
							# the above two if statements implicitly (in a fuzzy way) enforce
							# these two assert statements
							#assert d[z][y][x] >= cmin
							#assert d[z][y][x] <= cmax

    def test_threshold_notzero(self):
        """test threshold.notzero processor ................."""
        e = EMData()
        e.set_size(32,32,32)
        
        e.to_zero()
        e.process_inplace('threshold.notzero')
        d = e.get_3dview()
        for x in range(32):
            for y in range(32):
                for z in range(32):
                    self.assertEqual(d[z][y][x], 0)
        
        e.to_one()
        e *= 0.5
        e.process_inplace('threshold.notzero')
        d = e.get_3dview()
        for x in range(32):
            for y in range(32):
                for z in range(32):
                    self.assertEqual(d[z][y][x], 1)
   
    def test_math_squared(self):
        """test math.squared processor ......................"""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        d = e.get_3dview()
        e2 = e.copy()
        e2.process_inplace('math.squared')
        d2 = e2.get_3dview()
        for x in range(32):
            for y in range(32):
                for z in range(32):
                    self.assertAlmostEqual(d[z][y][x]**2, d2[z][y][x], 3)
                    
    def test_math_sqrt(self):
        """test math.sqrt processor ........................."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        d = e.get_3dview()
        e2 = e.copy()
        e2.process_inplace('math.sqrt')
        d2 = e2.get_3dview()
        from math import sqrt
        for x in range(32):
            for y in range(32):
                for z in range(32):
                    self.assertAlmostEqual(sqrt(d[z][y][x]), d2[z][y][x], 3) 
                    
    def test_threshold_belowtozero(self):
        """test threshold.belowtozero processor ............."""
        e = EMData()
        e.set_size(32,32,32)
        e.to_one()
        e *= 0.4
        
        e.process_inplace('threshold.belowtozero', {'minval':0.5})
        d = e.get_3dview()
        for x in range(32):
            for y in range(32):
                for z in range(32):
                    self.assertAlmostEqual(d[z][y][x], 0, 3)
                    
    def test_threshold_belowtozero_cut(self):
        """test threshold.belowtozero_cut processor ........."""
        e = EMData()
        e.set_size(32,32,32)
        e.to_one()
        e*=2
        
        d = e.get_3dview()
        for x in range(16):
            for y in range(32):
                for z in range(32):
                    d[x][y][z] -= 1.5
        
        e.process_inplace('threshold.belowtozero_cut', {'minval':1})
        
        for x in range(16):
            for y in range(32):
                for z in range(32):
                    self.assertAlmostEqual(d[x][y][z], 0.0, 3)
                    self.assertAlmostEqual(d[x+16][y][z], 1.0, 3)
                            
    def test_BinarizeProcessor(self):
        """test binary processor ............................"""
        imgfile1 = "test_BinarizeFilter.mrc"
        TestUtil.make_image_file(imgfile1, IMAGE_MRC)
        e = EMData()
        e.read_image(imgfile1)
        fnum = 1000
        f1 = Processors.get("threshold.binary", {'value': fnum})
        new_params = f1.get_params()
        self.assertEqual(float(new_params["value"]), fnum)
        f1.process_inplace(e)
        testlib.check_emdata(e, sys.argv[0])
        testlib.safe_unlink(imgfile1)

    def test_threshold_compress(self):
        """test threshold.compress processor ................"""
        e = EMData()
        e.set_size(32,32,32)
        e.to_one()
        e *= 0.45
        e.process_inplace('threshold.compress', {'range':0.2, 'value':0.5})
        d = e.get_3dview()
        for x in range(32):
            for y in range(32):
                for z in range(32):
                    self.assertAlmostEqual(d[z][y][x], 0.5, 3)
        
        e2 = EMData()
        e2.set_size(32,32,32)
        e2.to_one()
        e2 *= 0.05
        e2.process_inplace('threshold.compress', {'range':0.2, 'value':0.5})
        d2 = e2.get_3dview()
        for x in range(32):
            for y in range(32):
                for z in range(32):
                    self.assertAlmostEqual(d2[z][y][x], 0.2+0.05, 3)
        
    def test_math_linear(self):
        """test math.linear processor ......................."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        d = e.get_3dview()
        
        e2 = e.copy()
        e2.process_inplace('math.linear', {'scale':3.23, 'shift':2.56})
        d2 = e2.get_3dview()
        for x in range(32):
            for y in range(32):
                for z in range(32):
                    self.assertAlmostEqual(d2[z][y][x], d[z][y][x]*3.23 + 2.56, 3)
        
    def test_math_exp(self):
        """test math.exp processor .........................."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        d = e.get_3dview()
        
        e2 = e.copy()
        e2.process_inplace('math.exp', {'low':2.56, 'high':3.23})
        d2 = e2.get_3dview()
        from math import exp
        for x in range(32):
            for y in range(32):
                for z in range(32):
                    self.assertAlmostEqual(d2[z][y][x], exp(d[z][y][x]/2.56 - 3.23), 3)
    
    def test_RangeThreshold(self):
        """test range threshhold processor .................."""
        imgfile1 = "test_RangeThreshold.mrc"
        TestUtil.make_image_file(imgfile1, IMAGE_MRC)
        e = EMData()
        e.read_image(imgfile1)
        low = 10
        high = 20
        f1 = Processors.get("threshold.binaryrange", {"low":low, "high":high});
        d1 = f1.get_params()
        self.assertEqual(d1["low"], low)
        self.assertEqual(d1["high"], high)

        #xc = 12
        rw = 12.5
        #f2 = Processors.get("mask.ringmean", {"xc":xc, "ring_width":rw})
        f2 = Processors.get("mask.ringmean", {"ring_width":rw})
        d2 = f2.get_params()
        #self.assertEqual(d2["xc"], xc)
        self.assertEqual(d2["ring_width"], rw)

        outfile1 = "test_RangeThreshold_out.mrc"
        
        e.process_inplace("threshold.binary", {"value": 200})
        e.write_image(outfile1)

        testlib.safe_unlink(imgfile1)
        testlib.safe_unlink(outfile1)

    def test_math_sigma(self):
        """test math.sigma processor ........................"""
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
        e2.process_inplace('math.sigma', {'value1':v1, 'value2':v2})
        d2 = e2.get_3dview()
        for x in range(32):
            for y in range(32):
                for z in range(32):
                    if(d[z][y][x] < (mean - v2 * sigma) or d[z][y][x] > (mean + v1 * sigma)):
                        self.assertAlmostEqual(d2[z][y][x], mean, 3)
                    else:
                        self.assertAlmostEqual(d2[z][y][x], d[z][y][x], 3)
                        
    def test_math_log(self):
        """test math.log processor .........................."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        e -= 0.5
        d = e.get_3dview()
        
        e2 = e.copy()
        max = e2.get_attr('maximum')
        e2.process_inplace('math.log')
        d2 = e2.get_3dview()
        from math import log10
        for x in range(32):
            for y in range(32):
                for z in range(32):
                    if( d[z][y][x] > 0 ):
                        self.assertAlmostEqual(d2[z][y][x], log10(d[z][y][x]), 3)
                    else:
                        self.assertAlmostEqual(d2[z][y][x], 0.0, 3)
    
    def test_mask_sharp(self):
        """test mask.sharp processor ........................"""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e.process_inplace('mask.sharp', {'value':0.5})
        
    def test_mask_ringmean(self):
        """test mask.ringmean processor ....................."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e.process_inplace('mask.ringmean', {'ring_width':10})
        
    def test_mask_noise(self):
        """test mask.noise processor ........................"""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e.process_inplace('mask.noise')
        
    def test_mask_gaussian(self):
        """test mask.gaussian processor ....................."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e.process_inplace('mask.gaussian')
        
    def test_math_gausskernelfix(self):
        """test math.gausskernelfix processor ..............."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e.process_inplace('math.gausskernelfix', {'gauss_width':10.0})
        
    def test_math_toradiussqr(self):
        """test math.toradiussqr processor .................."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e.process_inplace('math.toradiussqr')
        
    def test_math_toradius(self):
        """test math.toradius processor ....................."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e.process_inplace('math.toradius')
        
    #need fix, this processor has problem    
    def no_test_complex_normpixels(self):
        """test complex.normpixels processor ................"""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        e.do_fft_inplace()
        self.assertEqual(e.is_complex(), True)
        
        e0 = e.get_fft_phase()    #phase before nomalization
        d0 = e0.get_3dview()
        
        e.process_inplace('complex.normpixels')
        
        e2 = e.get_fft_amplitude()    #amplitude to 1.0
        d2 = e2.get_3dview()
        e3 = e.get_fft_phase()    #phase after normalization, unchanged
        d3 = e3.get_3dview()
        for x in range(32):
            for y in range(32):
                for z in range(32):
                    self.assertAlmostEqual(d2[z][y][x], 1.0, 3)
                    self.assertAlmostEqual(d0[z][y][x]. d3[z][y][x], 3)
                    
    def no_test_math_laplacian(self):
        """test math.laplacian processor ...................."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e.process_inplace('math.laplacian')
        
    def test_mask_contract(self):
        """test mask.contract processor ....................."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e.process_inplace('mask.contract')
        
        e2 = EMData()
        e2.set_size(32,32,32)
        e2.to_zero()
        e2.process_inplace('mask.contract')
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
        
    def test_math_localsigma(self):
        """test math.localsigma processor ..................."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e.process_inplace('math.localsigma')
        
    def test_math_localmax(self):
        """test math.localmax processor ....................."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e.process_inplace('math.localmax')
        
    def test_math_submax(self):
        """test math.submax processor ......................."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e.process_inplace('math.submax')
        
    def test_mask_onlypeaks(self):
        """test mask.onlypeaks processor ...................."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e.process_inplace('mask.onlypeaks', {'npeaks':2})
        
    def no_test_eman1_filter_blockrange(self):
        """test eman1.filter.blockrange processor ..........."""
        e = EMData()
        e.set_size(32,32,1)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e.process_inplace('eman1.filter.blockrange', {'cal_half_width':0.2, 'fill_half_width':0.3})
        
        if(IS_TEST_EXCEPTION):
            #3D image not supported by this processor
            e2 = EMData()
            e2.set_size(32,32,32)
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
        
        if(IS_TEST_EXCEPTION):
            #3D image not supported by this processor
            e2 = EMData()
            e2.set_size(32,32,32)
            self.assertRaises( RuntimeError, e2.process, 'eman1.filter.blockcutoff', {'value1':0.2, 'value2':0.3} )
            try:
                e2.process_inplace('eman1.filter.blockcutoff', {'value1':0.2, 'value2':0.3})
            except RuntimeError, runtime_err:
                self.assertEqual(exception_type(runtime_err), "ImageDimensionException")
        
    def no_test_math_lineargradientfix(self):
        """test math.lineargradientfix processor ............"""
        e = EMData()
        e.set_size(32,32, 1)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e.process_inplace('math.lineargradientfix')
        
        if(IS_TEST_EXCEPTION):
            #3D image not supported by this processor
            e2 = EMData()
            e2.set_size(32,32,32)
            self.assertRaises( RuntimeError, e2.process, 'math.lineargradientfix' )
            try:
                e2.process_inplace('math.lineargradientfix')
            except RuntimeError, runtime_err:
                self.assertEqual(exception_type(runtime_err), "ImageDimensionException")
        
    def test_filter_ramp(self):
        """test filter.ramp processor ......................."""
        e = EMData()
        e.set_size(32,32, 1)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e.process_inplace('filter.ramp')
        
        if(IS_TEST_EXCEPTION):
            #3D image not supported by this processor
            e2 = EMData()
            e2.set_size(32,32,32)
            self.assertRaises( RuntimeError, e2.process_inplace, 'filter.ramp' )
            try:
                e2.process_inplace('filter.ramp')
            except RuntimeError, runtime_err:
                self.assertEqual(exception_type(runtime_err), "ImageDimensionException")
        
    def test_math_verticalstripefix(self):
        """test math.verticalstripefix processor ............"""
        e = EMData()
        e.set_size(32,32, 32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e.process_inplace('math.verticalstripefix')
        
    def test_math_realtofft(self):
        """test math.realtofft processor ...................."""
        e = EMData()
        e.set_size(32,32,1)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e.process_inplace('math.realtofft')
        
        if(IS_TEST_EXCEPTION):
            #3D image not supported by this processor
            e2 = EMData()
            e2.set_size(32,32,32)
            self.assertRaises( RuntimeError, e2.process_inplace, 'math.realtofft' )
            try:
                e2.process_inplace('math.realtofft')
            except RuntimeError, runtime_err:
                self.assertEqual(exception_type(runtime_err), "ImageDimensionException")
                
            #apply to real image only
            e3 = EMData()
            e3.set_size(32,32,1)
            e3.process_inplace('testimage.noise.uniform.rand')
            e3.do_fft_inplace()
            self.assertRaises( RuntimeError, e3.process_inplace, 'math.realtofft' )
            try:
                e3.process_inplace('math.realtofft')
            except RuntimeError, runtime_err:
                self.assertEqual(exception_type(runtime_err), "ImageFormatException")
            
    def test_mask_zeroedgefill(self):
        """test mask.zeroedgefill processor ................."""
        e = EMData()
        e.set_size(32,32,1)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e.process_inplace('mask.zeroedgefill')
        
        if(IS_TEST_EXCEPTION):
            #3D image not supported by this processor
            e2 = EMData()
            e2.set_size(32,32,32)
            self.assertRaises( RuntimeError, e2.process_inplace, 'mask.zeroedgefill' )
            try:
                e2.process_inplace('mask.zeroedgefill')
            except RuntimeError, runtime_err:
                self.assertEqual(exception_type(runtime_err), "ImageDimensionException")
            
    def test_mask_beamstop(self):
        """test mask.beamstop processor ....................."""
        e = EMData()
        e.set_size(32,32,1)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e.process_inplace('mask.beamstop', {'value1':-1.2, 'value2':10, 'value3':8})
        
        if(IS_TEST_EXCEPTION):
            #3D image not supported by this processor
            e2 = EMData()
            e2.set_size(32,32,32)
            self.assertRaises( RuntimeError, e2.process_inplace, 'mask.beamstop', \
                                {'value1':-1.2, 'value2':10, 'value3':8} )
            try:
                e2.process_inplace('mask.beamstop', {'value1':-1.2, 'value2':10, 'value3':8})
            except RuntimeError, runtime_err:
                self.assertEqual(exception_type(runtime_err), "ImageDimensionException")
            
    def test_mask_dampedzeroedgefill(self):
        """test mask.dampedzeroedgefill processor ..........."""
        e = EMData()
        e.set_size(32,32,1)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e.process_inplace('mask.dampedzeroedgefill')
        
        if(IS_TEST_EXCEPTION):
            #3D image not supported by this processor
            e2 = EMData()
            e2.set_size(32,32,32)
            self.assertRaises( RuntimeError, e2.process_inplace, 'mask.dampedzeroedgefill' )
            try:
                e2.process_inplace('mask.dampedzeroedgefill')
            except RuntimeError, runtime_err:
                self.assertEqual(exception_type(runtime_err), "ImageDimensionException")
            
    def test_math_averageovery(self):
        """test math.averageovery processor ................."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e.process_inplace('math.averageovery')
        
    def test_mask_zeroedge2d(self):
        """test mask.zeroedge2d processor ..................."""
        e = EMData()
        e.set_size(32,32,1)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e.process_inplace('mask.zeroedge2d', {'x0':1, 'y0':1, 'x1':25, 'y1':30})
        
        if(IS_TEST_EXCEPTION):
            #3D image not supported by this processor
            e2 = EMData()
            e2.set_size(32,32,32)
            self.assertRaises( RuntimeError, e2.process_inplace, 'mask.zeroedge2d', {'x0':1, 'y0':1, 'x1':25, 'y1':30})
            try:
                e2.process_inplace('mask.zeroedge2d', {'x0':1, 'y0':1, 'x1':25, 'y1':30})
            except RuntimeError, runtime_err:
                self.assertEqual(exception_type(runtime_err), "ImageDimensionException")
            
    def test_mask_zeroedge3d(self):
        """test mask.zeroedge3d processor ..................."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e.process_inplace('mask.zeroedge3d', {'x0':2, 'y0':3, 'z0':4, 'x1':20, 'y1':22, 'z1':26})
        
        if(IS_TEST_EXCEPTION):
            #3D image only
            e2 = EMData()
            e2.set_size(32,32,1)
            self.assertRaises( RuntimeError, e2.process_inplace, 'mask.zeroedge3d', {'x0':2, 'y0':3, 'z0':4, 'x1':20, 'y1':22, 'z1':26})
            try:
                e2.process_inplace('mask.zeroedge3d', {'x0':2, 'y0':3, 'z0':4, 'x1':20, 'y1':22, 'z1':26})
            except RuntimeError, runtime_err:
                self.assertEqual(exception_type(runtime_err), "ImageDimensionException")
            
    def test_bilateral(self):
        """test bilateral processor ........................."""
        e = EMData()
        e.set_size(16,16,16)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e.process_inplace('bilateral', {'distance_sigma':0.3, 'value_sigma':0.4, 'niter':2, 'half_width':5})
        
        e2 = EMData()
        e2.set_size(32,32)
        e2.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e2.is_complex(), False)
        e.process_inplace('bilateral', {'distance_sigma':0.3, 'value_sigma':0.4, 'niter':2, 'half_width':5})
            
        
    def test_normalize_unitlen(self):
        """test normalize.unitlen processor ................."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e.process_inplace('normalize.unitlen')
        
    def test_normalize_unitsum(self):
        """test normalize.unitsum processor ................."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e.process_inplace('normalize.unitsum')
        
    def test_normalize(self):
        """test normalize processor ........................."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
                
        e.process_inplace('normalize')
        sigma = e.get_attr('sigma')
        self.assertAlmostEqual(sigma, 1, 3)
        
        e2 = EMData(32,32)
        e2.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e2.is_complex(), False)
        e2.process_inplace('normalize')
        sigma = e2.get_attr('sigma')
        self.assertAlmostEqual(sigma, 1, 3)
        
        
    def test_normalize_mask(self):
        """test normalize.mask processor ...................."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e2 = EMData()
        e2.set_size(32,32,32)
        e2.process_inplace('testimage.noise.uniform.rand')
        
        e3 = e.process('normalize.mask', {'mask':e2, 'no_sigma':1})
        
        e4 = EMData()
        e4.set_size(16,16,16)
        e4.process_inplace('testimage.noise.uniform.rand')
        
        if(IS_TEST_EXCEPTION):
            try:
                e5 = e.process('normalize.mask', {'mask':e4, 'no_sigma':1})
            except RuntimeError, runtime_err:
                self.assertEqual(exception_type(runtime_err), "ImageDimensionException")
        
        
    def test_normalize_edgemean(self):
        """test normalize.edgemean processor ................"""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e.process_inplace('normalize.edgemean')
        
    def test_normalize_circlemean(self):
        """test normalize.circlemean processor .............."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e.process_inplace('normalize.circlemean')
        
    def test_normalize_lredge(self):
        """test normalize.lredge processor .................."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e.process_inplace('normalize.lredge')
        
    def test_normalize_maxmin(self):
        """test normalize.maxmin processor .................."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e.process_inplace('normalize.maxmin')
        
    def test_normalize_rows(self):
        """test normalize.rows processor ...................."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e.process_inplace('normalize.rows')
        
    def test_normalize_toimage(self):
        """test normalize.toimage processor ................."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e2 = EMData()
        e2.set_size(32,32,32)
        e2.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e2.is_complex(), False)
        
        e.process_inplace('normalize.toimage', {'noisy':e2, 'keepzero':2, 'invert':1, \
                                                'mult':2.3, 'add':0.5})
                                                
    def test_normalize_tofile(self):
        """test normalize.tofile processor .................."""
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
        
        e.process_inplace('normalize.tofile', {'noisyfile':filename, 'keepzero':2, \
                   'invert':1, 'mult':2.3, 'add':0.5 })
                   
        testlib.safe_unlink(filename)
        
    def test_normalize_toimage_lsq(self):
        """test normalize.toimage.lsq processor ............."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e2 = EMData()
        e2.set_size(32,32,32)
        e2.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e2.is_complex(), False)
        
        e.process_inplace('normalize.toimage.lsq', {'to':e2, 'low_threshold':0.2, 'high_threshold':0.8})
        
    def test_math_radialaverage(self):
		"""test math.radialaverage processor ................"""
		e = EMData()
		e.set_size(32,32,32)
		e.process_inplace('testimage.noise.uniform.rand')
		self.assertEqual(e.is_complex(), False)
		
		e.process_inplace('math.radialaverage')
		
		# test to make sure radial average is centered correctly for 2D
		# the 4 pixels that are exactly one pixel away from the center should be equal...
		# odd images
		e = EMData(5,5)
		e.set(2,3,1)
		e.process_inplace("math.radialaverage")
		self.assertEqual(e.get(2,3), e.get(2,1))
		self.assertEqual(e.get(2,3), e.get(3,2))
		self.assertEqual(e.get(2,3), e.get(1,2))
		
		# even images
		e = EMData(4,4)
		e.set(2,3,1)
		e.process_inplace("math.radialaverage")
		self.assertEqual(e.get(2,3), e.get(2,1))
		self.assertEqual(e.get(2,3), e.get(3,2))
		self.assertEqual(e.get(2,3), e.get(1,2))
		
		# test to make sure radial average is centered correctly for 3D
		# the 6 pixels that are exactly one pixel away from the center should be equal...
		e = EMData(5,5,5)
		e.set(2,2,3,1)
		e.process_inplace("math.radialaverage")
		self.assertEqual(e.get(2,2,3), e.get(2,2,1))
		self.assertEqual(e.get(2,2,3), e.get(2,3,2))
		self.assertEqual(e.get(2,2,3), e.get(2,1,2))
		self.assertEqual(e.get(2,2,3), e.get(3,2,2))
		self.assertEqual(e.get(2,2,3), e.get(1,2,2))
		
		e = EMData(4,4,4)
		e.set(2,2,3,1)
		e.process_inplace("math.radialaverage")
		self.assertEqual(e.get(2,2,3), e.get(2,2,1))
		self.assertEqual(e.get(2,2,3), e.get(2,3,2))
		self.assertEqual(e.get(2,2,3), e.get(2,1,2))
		self.assertEqual(e.get(2,2,3), e.get(3,2,2))
		self.assertEqual(e.get(2,2,3), e.get(1,2,2))
		
		
    def test_math_radialsubtract(self):
        """test math.radialsubtract processor ..............."""	
        e = EMData()
        e.set_size(32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e.process_inplace('math.radialsubtract')
        
        if(IS_TEST_EXCEPTION):
            e.set_size(32,32,32)
            e.process_inplace('testimage.noise.uniform.rand')
    		# It only currently works for 3D - if that ever changes this test should be updated
            try: e.process_inplace('testimage.noise.uniform.rand')
            except RuntimeError, runtime_err:
    			self.assertEqual(exception_type(runtime_err), "ImageDimensionException")
		
    def test_xform_flip(self):
		"""test xform.flip processor ........................"""
		
		# First test is just make sure flipped coordinates as are we expect.
		# Strategy - make some pixels non zero, flip the image, make sure the non
		# zero pixels end up where we expect them to be
		# X FLIPPING
		ims = [EMData(32,1,1),EMData(32,32,1),EMData(32,32,32),EMData(33,1,1),EMData(33,33,1),EMData(33,33,33)]
		for a in ims:
			
			nx = a.get_xsize()
			offset = nx%2==0
			for i in xrange(0,5):
				a.to_zero()
				lst = [Util.get_irand(offset,nx/2-offset) for i in range(3)] # go from 1 because even dimension flip 0s the 0 pixel
				for j in lst:
					for y in range(a.get_ysize()):
						for z in range(a.get_zsize()):
							a.set(j,y,z,1)
				a.process_inplace("xform.flip",{"axis":"x"})
				for j in lst:
					for y in range(a.get_ysize()):
						for z in range(a.get_zsize()):
							self.failIf(a.get(nx-1-j+offset,y,z) != 1)
		# Y FLIPPING		
		ims = [EMData(32,32,1),EMData(32,32,32),EMData(33,33,1),EMData(33,33,33)]
		for a in ims:
			ny = a.get_ysize()
			offset = ny%2==0
			for i in xrange(0,5):
				a.to_zero()
				lst = [Util.get_irand(offset,ny/2-offset) for i in range(3)] # go from 1 because even dimension flip 0s the 0 pixel
				for j in lst:
					for x in range(a.get_xsize()):
						for z in range(a.get_zsize()):
							a.set(x,j,z,1)
				a.process_inplace("xform.flip",{"axis":"y"})
				for j in lst:
					for x in range(a.get_xsize()):
						for z in range(a.get_zsize()):
							self.failIf(a.get(x,ny-1-j+offset,z) != 1)
		
		# Z FLIPPING				
		ims = [EMData(32,32,32),EMData(33,33,33)]
		for a in ims:
			nz = a.get_zsize()
			offset = nz%2==0
			for i in xrange(0,5):
				a.to_zero()
				lst = [Util.get_irand(offset,nz/2-offset) for i in range(3)] # go from 1 because even dimension flip 0s the 0 pixel
				for j in lst:
					for x in range(a.get_xsize()):
						for y in range(a.get_ysize()):
							a.set(x,y,j,1)
				a.process_inplace("xform.flip",{"axis":"z"})
				for j in lst:
					for x in range(a.get_xsize()):
						for y in range(a.get_ysize()):
							self.failIf(a.get(x,y,nz-1-j+offset) != 1)				
		
		# MIRROR TRANSFORMS IS THE SAME AS HORIZONTAL FLIP
		t = Transform()
		t.set_mirror(1)
		# This test should be True, It is an important part of ensuring the aligners are functioning accurately
		for n in [8,9]:
			a = test_image(1,size=(n,n))
			b = a.copy()
			a.process_inplace("xform.flip",{"axis":"x"})
			b.process_inplace("math.transform",{"transform":t})
			if n % 2 == 0:
				r = Region(1,0,n,n)
				aa = a.get_clip(r)
				bb = b.get_clip(r)
				self.assertEqual(aa==bb,True)
			self.assertEqual(a==b, True)
				
		# The 3D test is not as important as the 2D case (above), but it being true means we're in good shape - future
		# developers won't inadvertently make mistakes by interchanging the xform.flip and math.transform processors
		for n in [8,9]:
			a = test_image_3d(6,size=(n,n,n))
			b = a.copy()
			a.process_inplace("xform.flip",{"axis":"x"})
			b.process_inplace("math.transform",{"transform":t})
			if n % 2 == 0:
				r = Region(1,0,0,n,n,n)
				aa = a.get_clip(r)
				bb = b.get_clip(r)
				self.assertEqual(aa==bb,True)
			self.assertEqual(a==b, True)
	
		# ODD INVERTIBILITY
		a = EMData(33,1,1)
		a.process_inplace("testimage.noise.gauss")
		b = a.copy()
		b.process_inplace("xform.flip",{"axis":"x"})
		b.process_inplace("xform.flip",{"axis":"x"})
		self.assertEqual(a==b, True)
		  
		a = test_image(1,size=(33,33))
		b = a.copy()
		b.process_inplace("xform.flip",{"axis":"x"})
		b.process_inplace("xform.flip",{"axis":"x"})
		self.assertEqual(a==b, True)
		b.process_inplace("xform.flip",{"axis":"y"})
		b.process_inplace("xform.flip",{"axis":"y"})
		self.assertEqual(a==b, True)
		
		a = test_image_3d(6,size=(33,33,33))
		b = a.copy()
		b.process_inplace("xform.flip",{"axis":"x"})
		b.process_inplace("xform.flip",{"axis":"x"})
		self.assertEqual(a==b, True)
		b.process_inplace("xform.flip",{"axis":"y"})
		b.process_inplace("xform.flip",{"axis":"y"})
		self.assertEqual(a==b, True)
		b.process_inplace("xform.flip",{"axis":"z"})
		b.process_inplace("xform.flip",{"axis":"z"})
		self.assertEqual(a==b, True)
    
    	# EVEN INVERTIBILITY
    	# We have to do clipping in cases where the dimension is even, because 
    	# the equivalence of the clipped region is all we can guarantee
		a = EMData(32,1,1)
		a.process_inplace("testimage.noise.gauss")
		b = a.copy()
		b.process_inplace("xform.flip",{"axis":"x"})
		b.process_inplace("xform.flip",{"axis":"x"})
		r = Region(1,31)
		aa = a.get_clip(r)
		bb = b.get_clip(r)
		self.assertEqual(aa==bb, True)
		
		
		a = test_image(1,size=(32,32))
		b = a.copy()
		b.process_inplace("xform.flip",{"axis":"x"})
		b.process_inplace("xform.flip",{"axis":"x"})
		r = Region(1,0,31,32)
		aa = a.get_clip(r)
		bb = b.get_clip(r)
		self.assertEqual(aa==bb, True)
		a = test_image(1,size=(32,32))
		b = a.copy()
		b.process_inplace("xform.flip",{"axis":"y"})
		b.process_inplace("xform.flip",{"axis":"y"})
		r = Region(0,1,32,31)
		aa = a.get_clip(r)
		bb = b.get_clip(r)
		self.assertEqual(aa==bb, True)
		a = test_image_3d(6,size=(32,32,32))
		b = a.copy()
		b.process_inplace("xform.flip",{"axis":"x"})
		b.process_inplace("xform.flip",{"axis":"x"})
		r = Region(1,0,0,31,32,32)
		aa = a.get_clip(r)
		bb = b.get_clip(r)
		self.assertEqual(aa==bb, True)
		a = test_image_3d(6,size=(32,32,32))
		b = a.copy()
		b.process_inplace("xform.flip",{"axis":"y"})
		b.process_inplace("xform.flip",{"axis":"y"})
		r = Region(0,1,0,32,31,32)
		aa = a.get_clip(r)
		bb = b.get_clip(r)
		self.assertEqual(aa==bb, True)
		a = test_image_3d(6,size=(32,32,32))
		b = a.copy()
		b.process_inplace("xform.flip",{"axis":"z"})
		b.process_inplace("xform.flip",{"axis":"z"})
		r = Region(0,0,1,32,32,31)
		aa = a.get_clip(r)
		bb = b.get_clip(r)
		self.assertEqual(aa==bb, True)
    
    def test_math_addnoise(self):
        """test math.addnoise processor ....................."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e.process_inplace('math.addnoise', {'noise':0.89})
        
    def test_math_addsignoise(self):
        """test math.addsignoise processor .................."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e.process_inplace('math.addsignoise')
        
    def test_addspectralnoise(self):
        """test addspectralnoise processor .................."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        e.do_fft_inplace()
        self.assertEqual(e.is_complex(), True)
        
        e.process_inplace('addspectralnoise', {'n':2, 'x0':9.8, 'dx':1.2, \
                                    'y':(1.3,2.4), 'interpolation':1})
        
        if(IS_TEST_EXCEPTION):
            #complex image only
            e2 = EMData()
            e2.set_size(32,32,32)
            e2.process_inplace('testimage.noise.uniform.rand')
            self.assertEqual(e2.is_complex(), False)
            self.assertRaises( RuntimeError, e2.process_inplace, 'addspectralnoise', \
                                {'n':2, 'x0':9.8, 'dx':1.2, 'y':(1.3,2.4), 'interpolation':1})
            try:
                e2.process_inplace('addspectralnoise', {'n':2, 'x0':9.8, 'dx':1.2, \
                                        'y':(1.3,2.4), 'interpolation':1})
            except RuntimeError, runtime_err:
                self.assertEqual(exception_type(runtime_err), "ImageFormatException")
    
    def test_xform_fourierorigin(self):
		"""test xform.fourierorigin processor ..............."""
		e = EMData()
		e.set_size(32,32,1)
		e.process_inplace('testimage.noise.uniform.rand')
		self.assertEqual(e.is_complex(), False)
		e.do_fft_inplace()
		self.assertEqual(e.is_complex(), True)

		n = 16
		# test that 2D works
		for i in [n,n+1,n+2,n+3]:
			for j in [n,n+1,n+2,n+3]:
				for k in [1,n,n+1,n+2,n+3]:
					e = EMData()
					e.set_size(i,j,k)
				
					e.process_inplace("testimage.noise.uniform.rand")
					e.do_fft_inplace()
			
					d = e.copy()
	
					e.process_inplace("xform.fourierorigin.tocenter")
					e.process_inplace("xform.fourierorigin.tocorner")
					for kk in range(e.get_zsize()):
						for jj in range(e.get_ysize()):
							for ii in range(e.get_xsize()):
								self.assertEqual(e.get_value_at(ii,jj,kk), d.get_value_at(ii,jj,kk))

    def test_xform_phaseorigin_twostage(self):
		"""test xform.phaseorigin processor ................."""
		e = EMData()
		e.set_size(32,32,32)
		e.process_inplace('testimage.noise.uniform.rand')
		self.assertEqual(e.is_complex(), False)

		e.process_inplace('xform.phaseorigin.tocorner')
		e.process_inplace('xform.phaseorigin.tocenter')
		
		n = 8
		# first test that 1D works
		for i in [n,n+1]:
			e.set_size(i,1,1)
			e.process_inplace('testimage.noise.uniform.rand')
			f = e.copy()
			e.process_inplace('xform.phaseorigin.tocorner')
			e.process_inplace('xform.phaseorigin.tocenter')
			
			for ii in range(i):
				self.assertEqual(e.get_value_at(ii), f.get_value_at(ii))
		
		# now test that 2D works
		for i in [n,n+1,n+2,n+3]:
			for j in [n,n+1,n+2,n+3]:
				e.set_size(i,j,1)
				e.process_inplace('testimage.noise.uniform.rand')
				f = e.copy()
				e.process_inplace('xform.phaseorigin.tocorner')
				e.process_inplace('xform.phaseorigin.tocenter')
				
				for ii in range(i):
					for jj in range(j):
						self.assertEqual(e.get_value_at(ii,jj), f.get_value_at(ii,jj))
						
		# now test that 3D works
		for k in [n,n+1,n+2,n+3]:
			for j in [n,n+1,n+2,n+3]:
				for i in [n,n+1,n+2,n+3]:
					e.set_size(i,j,k)
					e.process_inplace('testimage.noise.uniform.rand')
					f = e.copy()
					e.process_inplace('xform.phaseorigin.tocorner')
					e.process_inplace('xform.phaseorigin.tocenter')

					for kk in range(k):
						for jj in range(j):
							for ii in range(i):
								self.assertEqual(e.get_value_at(ii,jj,kk), f.get_value_at(ii,jj,kk))

    def test_mask_auto2d(self):
        """test mask.auto2d processor ......................."""
        e = EMData()
        e.set_size(32,32,1)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e.process_inplace('mask.auto2d', {'threshold':0.5, 'filter':0.1})
        
        if(IS_TEST_EXCEPTION):
            #2D image only
            e2 = EMData()
            e2.set_size(32,32,32)
            e2.process_inplace('testimage.noise.uniform.rand')
            self.assertEqual(e.is_complex(), False)
            self.assertRaises( RuntimeError, e2.process_inplace, \
                'mask.auto2d', {'threshold':0.5, 'filter':0.1})
            try:
                e2.process_inplace('mask.auto2d', {'threshold':0.5, 'filter':0.1})
            except RuntimeError, runtime_err:
                self.assertEqual(exception_type(runtime_err), "ImageDimensionException")
    
    def no_test_mask_auto3d_thresh(self):
        """test mask.auto3d.thresh processor ................"""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e.process_inplace('mask.auto3d.thresh', {'threshold1':0.23, 'threshold2':0.86})
        
    def test_mask_auto3d(self):
        """test mask.auto3d processor ......................."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e.process_inplace('mask.auto3d', {'radius':16, 'threshold':0.5, 'nshells':3})
        
        testlib.safe_unlink('mask.mrc')
        
    def test_mask_addshells(self):
        """test mask.addshells processor ...................."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e.process_inplace('mask.addshells', {'nshells':3})
        
    def test_xform_phasecenterofmass(self):
		"""test xform.phasecenterofmass processor ..........."""
		self.centring_test("xform.phasecenterofmass",1)
    def test_xform_centerofmass(self):
		"""test xform.centerofmass processor ................"""
		self.centring_test("xform.centerofmass",1)
    def centring_test(self,processor_string,val=1):
		# 2D and 3D alignment test using a pixel value (val arg) offset one positive pixel from the origin in all directions
		# process_string should be either 
		if val == 0:
			print  "error, you can't use 0 as the value, the function is not capable of handling it"
			return
		for z in 1,8,9:
			for y in 8,9:
				for x in 8,9:
					
					e = EMData(x,y,z)
					if z != 1:
						e.set(x/2+1,y/2+1,z/2+1,val)
					else:
						e.set(x/2+1,y/2+1,val)
				
					e.process_inplace(processor_string)
					if val > 0:
						mx = e.calc_max_location()
					elif val < 0:
						mx = e.calc_min_location()
					
					self.failIf(mx[0] != x/2)
					self.failIf(mx[1] != y/2)
					
					if z != 1:
						self.failIf(mx[2] != z/2)

    def test_xform_centeracf(self):
		"""test xform.centeracf processor ..................."""
		self.centring_test("xform.centeracf",1)
		self.centring_test("xform.centeracf",-1)
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
        
    def est_misc_symsearch(self):
        """test misc.symsearch processor ...................."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e2 = EMData()
        e2.set_size(32,32,32)
        e2.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e2.is_complex(), False)
        
        e.process_inplace('misc.symsearch', {'sym':['CSYM'], 'thresh':0.5, 'output_symlabel':1, 'symlabel_map':e2})
        
    def test_misc_localnorm(self):
        """test misc.localnorm processor ...................."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        
        e.process_inplace('misc.localnorm', {'threshold':0.4, 'radius':16, 'apix':0.8})
        f = e.process('misc.localnorm', {'threshold':0.4, 'radius':16, 'apix':0.8})
        
        testlib.safe_unlink('norm.mrc')
        
    def test_mask_fromfile(self):
        """test mask.fromfile processor ....................."""
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
        
        e.process_inplace('mask.fromfile', {'filename':filename, 'ismaskset':1})
        
        testlib.safe_unlink(filename)
        
    def test_mask_fromfile_sizediff(self):
        """test mask.fromfile.sizediff processor ............"""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        e.set_attr('origin_row', 32)
        e.set_attr('origin_col', 32)
        e.set_attr('origin_sec', 32)
        
        e2 = EMData()
        e2.set_size(64,48,32)
        e2.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e2.is_complex(), False)
        filename = 'maskfile.mrc'
        e2.write_image(filename)
        
        e.process_inplace('mask.fromfile.sizediff', {'filename':filename})
        
        testlib.safe_unlink(filename)
        
    def no_test_misc_setpowspec(self):
        """test misc.setpowspec processor ..................."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        filename = 'powerspec.txt'
        
        e.process_inplace('misc.setpowspec', {'filename':filename})
        
    def test_mask_smart(self):
        """test mask.smart processor ........................"""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e.process_inplace('mask.smart', {'mask':1.1})
        
    def test_mask_addshells_gauss(self):
        """test mask.addshells.gauss processor .............."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        self.assertEqual(e.is_complex(), False)
        
        e.process_inplace('mask.addshells.gauss')
        
    def test_testimage_puregaussian(self):
        """test testimage.puregaussian processor ............"""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.gaussian', {'sigma':12})
        
    def test_testimage_gaussian(self):
        """test testimage.gaussian processor ................"""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.puregaussian', {'x_sigma':12, 'y_sigma':12, 'z_sigma':12})
        
    def test_testimage_scurve(self):
        """test testimage.scurve processor .................."""
        e = EMData()
        e.set_size(32,32,1)
        e.process_inplace('testimage.scurve')
        
        if(IS_TEST_EXCEPTION):
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
        e.process_inplace('testimage.sinewave', {'wavelength':10})
        e.process_inplace('testimage.sinewave', {'wavelength':10, 'phase':90})
        
        e2 = EMData()
        e2.set_size(32,32,1)
        e2.process_inplace('testimage.sinewave', {'wavelength':10})
        e2.process_inplace('testimage.sinewave', {'wavelength':10, 'axis':'y', 'phase':90, 'az':45})
        
        e3 = EMData()
        e3.set_size(32,32,32)
        e3.process_inplace('testimage.sinewave', {'wavelength':10})
        e3.process_inplace('testimage.sinewave', {'wavelength':10, 'axis':'z', 'phase':90, 'az':30, 'alt':45, 'phi':60})
            
    def test_testimage_sinewavecircular(self):
        """test testimage.sinewave.circular processor ......."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.sinewave.circular', {'wavelength':10})
        e.process_inplace('testimage.sinewave.circular', {'wavelength':10, 'axis':'z', 'c':7, 'phase':3})
        
    def test_testimage_squarecube(self):
        """test testimage.squarecube processor .............."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.squarecube', {'edge_length':22})
        e.process_inplace('testimage.squarecube', {'edge_length':22, 'axis':'y', 'odd_edge':18, 'fill':0})
        
    def test_testimage_circlesphere(self):
        """test testimage.circlesphere processor ............"""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.circlesphere', {'radius':20})
        e.process_inplace('testimage.circlesphere', {'radius':20, 'axis':'z', 'c':15, 'fill':1})
        
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
    
    def test_basis_wavelet(self):
        """test basis.wavelet processor ....................."""
        e1 = EMData()
        e1.set_size(64,64)
        e1.process_inplace('testimage.noise.uniform.rand')
        d1 = e1.get_2dview()
        
        e2 = e1.process('basis.wavelet', {'type':'daub', 'dir':1, 'ord':20})
        d2 = e2.get_2dview()
        
        e3 = e2.process('basis.wavelet', {'type':'daub', 'dir':-1, 'ord':20})        
        d3 = e3.get_2dview()
        
        for j in range(64):
            for i in range(64):
                self.assertAlmostEqual(d1[i][j], d3[i][j], 3)
    
    def test_basis_fft(self):
        """test basis.fft processor ........................."""
        e1 = EMData()
        e1.set_size(8,8,8)
        e1.process_inplace('testimage.noise.uniform.rand')
        d1 = e1.get_3dview()
        
        e2 = e1.process('basis.fft', {'dir':1})
        d2 = e2.get_3dview()
        
        e3 = e2.process('basis.fft', {'dir':-1})
        d3 = e3.get_3dview()
        
        e4 = e1.do_fft()
        d4 = e4.get_3dview()
        
        e5 = e2.do_ift()
        d5 = e5.get_3dview()
        
        for k in range(8):
            for j in range(8):
                for i in range(8):
                    self.assertAlmostEqual(d1[i][j][k], d3[i][j][k], 3)
                    self.assertAlmostEqual(d2[i][j][k], d4[i][j][k], 3)
                    self.assertAlmostEqual(d1[i][j][k], d5[i][j][k], 3)
        
    #this filter.integercyclicshift2d processor is removed by Phani at 5/18/2006    
    def no_test_IntegerCyclicShift2DProcessor(self):
        """test filter.integercyclicshift2d processor........"""
        e = EMData()
        e.set_size(50,50,1)
        e.process_inplace('testimage.noise.uniform.rand')
        
        d1 = e.get_2dview()
        d3 = numpy.array(d1)    #make a copy of d1, since d1 and d2 share the same memory
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
        
        if(IS_TEST_EXCEPTION):        
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
    p = OptionParser()
    p.add_option('--t', action='store_true', help='test exception', default=False )
    global IS_TEST_EXCEPTION
    opt, args = p.parse_args()
    if opt.t:
        IS_TEST_EXCEPTION = True
    Log.logger().set_level(-1)  #perfect solution for quenching the Log error information, thank Liwei
    suite = unittest.TestLoader().loadTestsFromTestCase(TestProcessor)
    unittest.TextTestRunner(verbosity=2).run(suite)

if __name__ == '__main__':
    test_main()
