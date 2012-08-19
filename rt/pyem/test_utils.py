#!/usr/bin/env python

#
# Author: Liwei Peng, 01/30/2005 (sludtke@bcm.edu)
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
import unittest
import os
import random
import math
import testlib
from optparse import OptionParser

IS_TEST_EXCEPTION = False

class TestUtils(unittest.TestCase):
    """test Util class"""
    
    def test_is_file_exist(self):
        """test is_file_exist() function ...................."""
        imgfile1 = "test_is_file_exist.mrc"
        TestUtil.make_image_file(imgfile1, IMAGE_MRC)
        result1 = Util.is_file_exist(imgfile1)
        self.assertEqual(result1, True)
        
        result2 = Util.is_file_exist("__nosuchafile__.dm3")
        self.assertEqual(result2, False)
        
        result3 = Util.is_file_exist("")
        self.assertEqual(result3, False)

        testlib.safe_unlink(imgfile1)

    def test_sstrcmp(self):
        """test sstrcmp() function .........................."""
        e1 = Util.sstrncmp("helloworld", "hello");
        self.assertEqual(e1, True);

        e2 = Util.sstrncmp("foobar", "bar");
        self.assertEqual(e2, False);

        e3 = Util.sstrncmp("", "bar");
        self.assertEqual(e3, False);

        e4 = Util.sstrncmp("bar", "");
        self.assertEqual(e4, True);

        e5 = Util.sstrncmp("", "");
        self.assertEqual(e5, True);

    def test_int2str(self):
        """test int2str() function .........................."""
        s1 = Util.int2str(123)
        self.assertEqual(s1, "123")

        s2 = Util.int2str(-1)
        self.assertEqual(s2, "-1")

    def test_change_filename_ext(self):
        """test change_filename_ext() function .............."""
        file1 = Util.change_filename_ext("hello.c", "cpp")
        file2 = Util.change_filename_ext("hello.", "mrc")
        file3 = Util.change_filename_ext("hello", "cpp")

        file4 = Util.change_filename_ext("hello.c", "")
        file5 = Util.change_filename_ext("hello.", "")
        file6 = Util.change_filename_ext("hello", "")

        self.assertEqual(file1, "hello.cpp")
        self.assertEqual(file2, "hello.mrc")
        self.assertEqual(file3, "hello.cpp")
        self.assertEqual(file4, "hello.c")
        self.assertEqual(file5, "hello.")
        self.assertEqual(file6, "hello")

    def test_remove_filename_ext(self):
        """test remove_filename_ext() function .............."""
        s1 = Util.remove_filename_ext("hello.cpp")
        self.assertEqual(s1, "hello")

        s1 = Util.remove_filename_ext("hello.")
        self.assertEqual(s1, "hello")

        s1 = Util.remove_filename_ext("hello")
        self.assertEqual(s1, "hello")

        s1 = Util.remove_filename_ext("")
        self.assertEqual(s1, "")
    
    def test_get_filename_ext(self):
        """test get_filename_ext() function ................."""
        s1 = Util.get_filename_ext("hello.cpp")
        self.assertEqual(s1, "cpp")
        
        s1 = Util.get_filename_ext("hello.")
        self.assertEqual(s1, "")
        
        s1 = Util.get_filename_ext("hello")
        self.assertEqual(s1, "")
        
        s1 = Util.get_filename_ext("")
        self.assertEqual(s1, "")    

    def test_sbasename(self):
        """test sbasename() function ........................"""
        b1 = Util.sbasename("hello.c")
        self.assertEqual(b1, "hello.c")
        
        import os
        if os.sys.platform == 'win32':
            b2 = Util.sbasename("\tmp\hello.mrc")
            b3 = Util.sbasename(".\test\hello")
        else:    
            b2 = Util.sbasename("/tmp/hello.mrc")
            b3 = Util.sbasename("./test/hello")
        self.assertEqual(b2, "hello.mrc")
        self.assertEqual(b3, "hello")

    def test_get_frand(self):
        """test get_frand() function ........................"""
        n1 = 10
        n2 = 20
        f1 = Util.get_frand(n1, n2)
        self.assert_(n1 <= f1)
        self.assert_(f1 <= n2)
        
    def test_get_gauss_rand(self):
        """test get_gauss_rand() function ..................."""
        mean = 0.0
        sigma = 2.71828
        gr = Util.get_gauss_rand(mean, sigma)
        
    def test_bilinear_interpolate(self):
        """test bilinear_interpolate() function ............."""
        i1 = Util.bilinear_interpolate(1.0, 1.0, 1.0, 1.0, 0.5, 0.5)
        self.assertAlmostEqual(i1, 1.0, 3)
        
        i2 = Util.bilinear_interpolate(1.0, 2.0, 3.0, 4.0, 0.5, 0.5)
        i3 = Util.bilinear_interpolate(1.0, 4.0, 2.0, 3.0, 0.5, 0.5)
        self.assertAlmostEqual(i2, i3, 3)

    def test_trilinear_interpolate(self):
        """test trilinear_interpolate() function ............"""
        i1 = Util.trilinear_interpolate(1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.5,0.5,0.5)
        self.assertAlmostEqual(i1, 1.0, 3)
        
        i2 = Util.trilinear_interpolate(1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,0.5,0.5,0.5)
        i3 = Util.trilinear_interpolate(8.0,4.0,2.0,7.0,5.0,3.0,6.0,1.0,0.5,0.5,0.5)
        self.assertAlmostEqual(i2, i3, 3)
        
    def test_round(self):
        """test round() function ............................"""
        self.assertEqual(Util.round(1), 1)
        self.assertEqual(Util.round(1.2), 1)
        self.assertEqual(Util.round(-2.4), -2)
        self.assertEqual(Util.round(1.7), 2)
        self.assertEqual(Util.round(-2.8), -3)

    def test_square(self):
        """test square() function ..........................."""
        for i in range(10):
            r1 = random.randint(1,100)
            self.assertEqual(Util.square(r1), r1*r1)

            r2 = r1 + 0.5
            self.assertEqual(Util.square(r2), r2*r2)

    def test_square_sum(self):
        """test square_sum() function ......................."""
        for i in range(10):
            r1 = random.randint(1,100) + 0.5
            r2 = random.randint(1,100) + 0.5
            s1 = r1*r1 + r2*r2
            s2 = Util.square_sum(r1,r2)            
            self.assertEqual(s1, s2)

    def test_hypot3(self):
        """test hypot3() function ..........................."""
        for i in range(10):
            r1 = random.randint(1,100) + 0.5
            r2 = random.randint(1,100) + 0.5
            r3 = random.randint(1,100) + 0.5
            h1 = math.sqrt(r1*r1+r2*r2+r3*r3)
            h2 = Util.hypot3(r1,r2,r3)
            testlib.assertfloat(self, h1, h2)
        
    def test_min(self):
        """test min() function .............................."""
        for i in range(10):
            r1 = random.randint(1,100)
            r2 = random.randint(1,100)
            r3 = random.randint(1,100)
            self.assertEqual(Util.get_min(r1, r2), min(r1, r2))
            self.assertEqual(Util.get_min(r1, r2, r3), min(r1, r2, r3))

    def test_max(self):
        """test max() function .............................."""
        for i in range(10):
            r1 = random.randint(1,100)
            r2 = random.randint(1,100)
            r3 = random.randint(1,100)
            self.assertEqual(Util.get_max(r1, r2), max(r1, r2))
            self.assertEqual(Util.get_max(r1, r2, r3), max(r1, r2, r3))

    def test_angle_sub_2pi(self):
        """test angle_sub_2pi() function ...................."""
        self.assertEqual(Util.angle_sub_2pi(3,1), 2)
        testlib.assertfloat(self, Util.angle_sub_2pi(6,1), 1.28318536);
        
    def test_angle_sub_pi(self):
        """test angle_sub_pi() function ....................."""
        self.assertEqual(Util.angle_sub_pi(2,1), 1)
        testlib.assertfloat(self, Util.angle_sub_pi(4,1), 0.14159)
        
    def test_get_time_label(self):
        """test get_time_label() function ..................."""
        s = Util.get_time_label()
        #ss = time.localtime()
        #sss=str(ss[1])+'/'+str(ss[2])+'/'+str(ss[0])+' '+str(ss[3])+':'+str(ss[4])
        #self.assertEqual(s[0:12],sss[0:12])
        #self.assertEqual(len(s), 15)
    
    def test_calc_best_fft_size(self):
        """test calc_best_fft_size() funciton ..............."""
        b = Util.calc_best_fft_size(1023)
        self.assertEqual(b, 1024)
        
    def test_fast_floor(self):
        """test fast_floor() function ......................."""
        f = 2.4
        self.assertEqual(Util.fast_floor(f), 2)
        
        f1 = -2.4
        self.assertEqual(Util.fast_floor(f1), -3)
        
    def test_agauss(self):
        """test agauss() function ..........................."""
        g = Util.agauss(1.2, 2.3, 3.4, 4.5, 5.6)
        
    def test_eman_copysign(self):
        """test eman_copysign() function ...................."""
        a = 12.3
        b1 = -5
        b2 = 5
        c1 = Util.eman_copysign(a, b1)
        c2 = Util.eman_copysign(a, b2)
        
        self.assertAlmostEqual(c1, -a, 3)
        self.assertAlmostEqual(c2, a, 3)
        
    def test_eman_erfc(self):
        """test eman_erfc() function ........................"""
        a = 123.456
        e = Util.eman_erfc(a)
    
    def test_set_randnum_seed(self):
        """test the set_randnum_seed() function ............."""
        a = []
        b = []
        SEED = 123456
        Util.set_randnum_seed(SEED)
        for i in xrange(10):
            a.append(Util.get_irand(1,100))
        
        Util.set_randnum_seed(SEED)
        for i in xrange(10):
            b.append(Util.get_irand(1,100))
        
        for i in xrange(10):
            self.assertEqual(a[i], b[i])
            
        seed = Util.get_randnum_seed()
        self.assertEqual(seed, SEED)
   
    # no more voea() functions
    def no_test_voea(self):
       """test voea() function ............................."""
       v = Util.voea(1.2)    #test default argument
       v2 = Util.voea(1.2, 20.0, 3.0, 401.0)
        
    #quadri() function changed, return image as EMData* in argument
    def no_test_quadri(self):
        """test quadri() function ..........................."""
        e = EMData()
        e.set_size(32,32,32)
        e.process_inplace('testimage.noise.uniform.rand')
        
        #f = Util.quadri(e, 2.0, 3.0)    #test default argument
        f2 = Util.quadri(e, 2.3, 3.4, 2)    #test non-default argument
        
class TestEMUtils(unittest.TestCase):
    """test EMUtil class"""
    
    def test_get_all_attibutes(self):
        """test get_all_attributes function ................."""
        e1 = test_image()
        e1.set_attr('name', 'Tom')
        e1.write_image('test.hdf', 0)
        
        e2 = test_image()
        e2.set_attr('name', 'Sam')
        e2.write_image('test.hdf', 1)
        
        e3 = test_image()
        e3.write_image('test.hdf', 2)
        
        l = EMUtil.get_all_attributes('test.hdf', 'name')
        self.assertEqual(l, ['Tom', 'Sam', None])
        
        testlib.safe_unlink('test.hdf')
    
    def test_is_same_size(self):
        """test is_same_size function ......................."""
        e1 = test_image()
        e2 = test_image()
        self.assertEqual(EMUtil.is_same_size(e1, e2), True)
        
        e3 = EMData()
        e3.set_size(32,32,32)
        e4 = EMData()
        e4.set_size(32,32,31)
        self.assertEqual(EMUtil.is_same_size(e3, e4), False)
        
        e5 = EMData()
        e5.set_size(32,32,32)
        e6 = EMData()
        e6.set_size(2,30,108)
        self.assertEqual(EMUtil.is_same_size(e5, e6), False)
        
    def test_is_complex_type(self):
        """test is_complex_type function ...................."""
        self.assertEqual(EMUtil.is_complex_type(EMUtil.EMDataType.EM_USHORT_COMPLEX), True)
        self.assertEqual(EMUtil.is_complex_type(EMUtil.EMDataType.EM_FLOAT_COMPLEX), True)
        self.assertEqual(EMUtil.is_complex_type(EMUtil.EMDataType.EM_SHORT_COMPLEX), True)
        
        self.assertEqual(EMUtil.is_complex_type(EMUtil.EMDataType.EM_CHAR), False)
        self.assertEqual(EMUtil.is_complex_type(EMUtil.EMDataType.EM_DOUBLE), False)
        self.assertEqual(EMUtil.is_complex_type(EMUtil.EMDataType.EM_FLOAT), False)
        self.assertEqual(EMUtil.is_complex_type(EMUtil.EMDataType.EM_INT), False)
        self.assertEqual(EMUtil.is_complex_type(EMUtil.EMDataType.EM_SHORT), False)
        self.assertEqual(EMUtil.is_complex_type(EMUtil.EMDataType.EM_UCHAR), False)
        self.assertEqual(EMUtil.is_complex_type(EMUtil.EMDataType.EM_UINT), False)
        self.assertEqual(EMUtil.is_complex_type(EMUtil.EMDataType.EM_UNKNOWN), False)
        self.assertEqual(EMUtil.is_complex_type(EMUtil.EMDataType.EM_USHORT), False)
        
    def test_get_euler_names(self):
        """test get_euler_names ............................."""
        l1 = EMUtil.get_euler_names('EMAN')
        self.assertEqual(l1, ['euler_alt', 'euler_az', 'euler_phi'])
        
        l2 = EMUtil.get_euler_names('MRC')
        self.assertEqual(l2, ['euler_theta', 'euler_phi', 'euler_omega'])
        
        l3 = EMUtil.get_euler_names('IMAGIC')
        self.assertEqual(l3, ['euler_alpha', 'euler_beta', 'euler_gamma'])
        
        l4 = EMUtil.get_euler_names('SPIDER')
        self.assertEqual(l4, ['euler_phi', 'euler_theta', 'euler_gamma'])
        
        l5 = EMUtil.get_euler_names('SPIN')
        self.assertEqual(l5, ['euler_q', 'euler_n1', 'euler_n2', 'euler_n3'])
        
        l6 = EMUtil.get_euler_names('SGIROT')
        self.assertEqual(l6, ['euler_q', 'euler_n1', 'euler_n2', 'euler_n3'])
        
        l7 = EMUtil.get_euler_names('QUATERNION')
        self.assertEqual(l7, ['euler_e0', 'euler_e1', 'euler_e2', 'euler_e3',])
        
    def test_is_same_ctf(self):
        """test is_same_ctf function ........................"""
        e1 = EMData()
        e1.set_size(32,32,1)
        e1.set_attr('ctf', [0.,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.])
        
        e2 = EMData()
        e2.set_size(64,64,1)
        e2.set_attr('ctf', [0.,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.])
        
        self.assertEqual(EMUtil.is_same_ctf(e1, e2), True)
        
        e1.set_attr('ctf', [0.,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.])
        e2.set_attr('ctf', [0.,1.,2.,3.,4.,5.,6.,7.,8.,9.,99.])
        self.assertEqual(EMUtil.is_same_ctf(e1, e2), False)
        
    def test_get_image_ext_type(self):
        """test get_image_ext_type function ................."""
        self.assertEqual(EMUtil.get_image_ext_type('mrc'), EMUtil.ImageType.IMAGE_MRC)
        self.assertEqual(EMUtil.get_image_ext_type('MRC'), EMUtil.ImageType.IMAGE_MRC)
        self.assertEqual(EMUtil.get_image_ext_type('tnf'), EMUtil.ImageType.IMAGE_MRC)
        self.assertEqual(EMUtil.get_image_ext_type('TNF'), EMUtil.ImageType.IMAGE_MRC)
        
        self.assertEqual(EMUtil.get_image_ext_type('dm3'), EMUtil.ImageType.IMAGE_DM3)
        self.assertEqual(EMUtil.get_image_ext_type('DM3'), EMUtil.ImageType.IMAGE_DM3)
        
        self.assertEqual(EMUtil.get_image_ext_type('dm4'), EMUtil.ImageType.IMAGE_DM4)
        self.assertEqual(EMUtil.get_image_ext_type('DM4'), EMUtil.ImageType.IMAGE_DM4)
        
        self.assertEqual(EMUtil.get_image_ext_type('spi'), EMUtil.ImageType.IMAGE_SPIDER)
        self.assertEqual(EMUtil.get_image_ext_type('SPI'), EMUtil.ImageType.IMAGE_SPIDER)
        self.assertEqual(EMUtil.get_image_ext_type('spider'), EMUtil.ImageType.IMAGE_SPIDER)
        self.assertEqual(EMUtil.get_image_ext_type('SPIDER'), EMUtil.ImageType.IMAGE_SPIDER)
        
        self.assertEqual(EMUtil.get_image_ext_type('spidersingle'), EMUtil.ImageType.IMAGE_SINGLE_SPIDER)
        self.assertEqual(EMUtil.get_image_ext_type('SPIDERSINGLE'), EMUtil.ImageType.IMAGE_SINGLE_SPIDER)
        self.assertEqual(EMUtil.get_image_ext_type('singlespider'), EMUtil.ImageType.IMAGE_SINGLE_SPIDER)
        self.assertEqual(EMUtil.get_image_ext_type('SINGLESPIDER'), EMUtil.ImageType.IMAGE_SINGLE_SPIDER)
        
        self.assertEqual(EMUtil.get_image_ext_type('img'), EMUtil.ImageType.IMAGE_IMAGIC)
        self.assertEqual(EMUtil.get_image_ext_type('IMG'), EMUtil.ImageType.IMAGE_IMAGIC)
        self.assertEqual(EMUtil.get_image_ext_type('hed'), EMUtil.ImageType.IMAGE_IMAGIC)
        self.assertEqual(EMUtil.get_image_ext_type('HED'), EMUtil.ImageType.IMAGE_IMAGIC)
        self.assertEqual(EMUtil.get_image_ext_type('imagic'), EMUtil.ImageType.IMAGE_IMAGIC)
        self.assertEqual(EMUtil.get_image_ext_type('IMAGIC'), EMUtil.ImageType.IMAGE_IMAGIC)
        
        self.assertEqual(EMUtil.get_image_ext_type('pgm'), EMUtil.ImageType.IMAGE_PGM)
        self.assertEqual(EMUtil.get_image_ext_type('PGM'), EMUtil.ImageType.IMAGE_PGM)
        
        self.assertEqual(EMUtil.get_image_ext_type('lst'), EMUtil.ImageType.IMAGE_LST)
        self.assertEqual(EMUtil.get_image_ext_type('LST'), EMUtil.ImageType.IMAGE_LST)
        
        self.assertEqual(EMUtil.get_image_ext_type('pif'), EMUtil.ImageType.IMAGE_PIF)
        self.assertEqual(EMUtil.get_image_ext_type('PIF'), EMUtil.ImageType.IMAGE_PIF)
        
        self.assertEqual(EMUtil.get_image_ext_type('png'), EMUtil.ImageType.IMAGE_PNG)
        self.assertEqual(EMUtil.get_image_ext_type('PNG'), EMUtil.ImageType.IMAGE_PNG)
        
        self.assertEqual(EMUtil.get_image_ext_type('h5'), EMUtil.ImageType.IMAGE_HDF)
        self.assertEqual(EMUtil.get_image_ext_type('H5'), EMUtil.ImageType.IMAGE_HDF)
        self.assertEqual(EMUtil.get_image_ext_type('hdf'), EMUtil.ImageType.IMAGE_HDF)
        self.assertEqual(EMUtil.get_image_ext_type('HDF'), EMUtil.ImageType.IMAGE_HDF)
        self.assertEqual(EMUtil.get_image_ext_type('hd5'), EMUtil.ImageType.IMAGE_HDF)
        self.assertEqual(EMUtil.get_image_ext_type('HD5'), EMUtil.ImageType.IMAGE_HDF)
        
        self.assertEqual(EMUtil.get_image_ext_type('tif'), EMUtil.ImageType.IMAGE_TIFF)
        self.assertEqual(EMUtil.get_image_ext_type('TIF'), EMUtil.ImageType.IMAGE_TIFF)
        self.assertEqual(EMUtil.get_image_ext_type('tiff'), EMUtil.ImageType.IMAGE_TIFF)
        self.assertEqual(EMUtil.get_image_ext_type('TIFF'), EMUtil.ImageType.IMAGE_TIFF)
        
        self.assertEqual(EMUtil.get_image_ext_type('vtk'), EMUtil.ImageType.IMAGE_VTK)
        self.assertEqual(EMUtil.get_image_ext_type('VTK'), EMUtil.ImageType.IMAGE_VTK)
        
        self.assertEqual(EMUtil.get_image_ext_type('hdr'), EMUtil.ImageType.IMAGE_SAL)
        self.assertEqual(EMUtil.get_image_ext_type('HDR'), EMUtil.ImageType.IMAGE_SAL)
        self.assertEqual(EMUtil.get_image_ext_type('sal'), EMUtil.ImageType.IMAGE_SAL)
        self.assertEqual(EMUtil.get_image_ext_type('SAL'), EMUtil.ImageType.IMAGE_SAL)
        
        self.assertEqual(EMUtil.get_image_ext_type('map'), EMUtil.ImageType.IMAGE_ICOS)
        self.assertEqual(EMUtil.get_image_ext_type('MAP'), EMUtil.ImageType.IMAGE_ICOS)
        self.assertEqual(EMUtil.get_image_ext_type('icos'), EMUtil.ImageType.IMAGE_ICOS)
        self.assertEqual(EMUtil.get_image_ext_type('ICOS'), EMUtil.ImageType.IMAGE_ICOS)
        
        self.assertEqual(EMUtil.get_image_ext_type('am'), EMUtil.ImageType.IMAGE_AMIRA)
        self.assertEqual(EMUtil.get_image_ext_type('AM'), EMUtil.ImageType.IMAGE_AMIRA)
        self.assertEqual(EMUtil.get_image_ext_type('amira'), EMUtil.ImageType.IMAGE_AMIRA)
        self.assertEqual(EMUtil.get_image_ext_type('AMIRA'), EMUtil.ImageType.IMAGE_AMIRA)
        
        self.assertEqual(EMUtil.get_image_ext_type('emim'), EMUtil.ImageType.IMAGE_EMIM)
        self.assertEqual(EMUtil.get_image_ext_type('EMIM'), EMUtil.ImageType.IMAGE_EMIM)
        
        self.assertEqual(EMUtil.get_image_ext_type('xplor'), EMUtil.ImageType.IMAGE_XPLOR)
        self.assertEqual(EMUtil.get_image_ext_type('XPLOR'), EMUtil.ImageType.IMAGE_XPLOR)
        
        self.assertEqual(EMUtil.get_image_ext_type('em'), EMUtil.ImageType.IMAGE_EM)
        self.assertEqual(EMUtil.get_image_ext_type('EM'), EMUtil.ImageType.IMAGE_EM)
        
        self.assertEqual(EMUtil.get_image_ext_type('dm2'), EMUtil.ImageType.IMAGE_GATAN2)
        self.assertEqual(EMUtil.get_image_ext_type('DM2'), EMUtil.ImageType.IMAGE_GATAN2)
        
        self.assertEqual(EMUtil.get_image_ext_type('v4l'), EMUtil.ImageType.IMAGE_V4L)
        self.assertEqual(EMUtil.get_image_ext_type('V4L'), EMUtil.ImageType.IMAGE_V4L)
        
        self.assertEqual(EMUtil.get_image_ext_type('jpg'), EMUtil.ImageType.IMAGE_JPEG)
        self.assertEqual(EMUtil.get_image_ext_type('JPG'), EMUtil.ImageType.IMAGE_JPEG)
        self.assertEqual(EMUtil.get_image_ext_type('jpeg'), EMUtil.ImageType.IMAGE_JPEG)
        self.assertEqual(EMUtil.get_image_ext_type('JPEG'), EMUtil.ImageType.IMAGE_JPEG)
        
        self.assertEqual(EMUtil.get_image_ext_type('xyz'), EMUtil.ImageType.IMAGE_UNKNOWN)
        
    def test_get_image_type(self):
        """test get_image_type function ....................."""
        e = test_image()
        e.write_image('mrcfile', 0, EMUtil.ImageType.IMAGE_MRC)
        self.assertEqual(EMUtil.get_image_type('mrcfile'), EMUtil.ImageType.IMAGE_MRC)
        testlib.safe_unlink('mrcfile')
        
        e2 = test_image()
        e2.write_image('hdffile', 0, EMUtil.ImageType.IMAGE_HDF)
        self.assertEqual(EMUtil.get_image_type('hdffile'), EMUtil.ImageType.IMAGE_HDF)
        testlib.safe_unlink('hdffile')
        
        #e.write_image('lstfile', 0, EMUtil.ImageType.IMAGE_LST)
        #self.assertEqual(EMUtil.get_image_type('lstfile'), EMUtil.ImageType.IMAGE_LST)
        #testlib.safe_unlink('lstfile')
        
        e3 = test_image()
        e3.write_image('spiderfile', 0, EMUtil.ImageType.IMAGE_SPIDER)
        self.assertEqual(EMUtil.get_image_type('spiderfile'), EMUtil.ImageType.IMAGE_SPIDER)
        testlib.safe_unlink('spiderfile')
        
        e3.write_image('sspiderfile', 0, EMUtil.ImageType.IMAGE_SINGLE_SPIDER)
        self.assertEqual(EMUtil.get_image_type('sspiderfile'), EMUtil.ImageType.IMAGE_SINGLE_SPIDER)
        testlib.safe_unlink('sspiderfile')
        
        #e.write_image('piffile', 0, EMUtil.ImageType.IMAGE_PIF)
        #self.assertEqual(EMUtil.get_image_type('piffile'), EMUtil.ImageType.IMAGE_PIF)
        #testlib.safe_unlink('piffile')
        
        platform = get_platform()
        if platform!='Windows' and platform!='win32':
            e4 = test_image()
            e4.write_image('pngfile', 0, EMUtil.ImageType.IMAGE_PNG)
            self.assertEqual(EMUtil.get_image_type('pngfile'), EMUtil.ImageType.IMAGE_PNG)
            testlib.safe_unlink('pngfile')
        
        e5 = test_image()
        e5.write_image('vtkfile', 0, EMUtil.ImageType.IMAGE_VTK)
        self.assertEqual(EMUtil.get_image_type('vtkfile'), EMUtil.ImageType.IMAGE_VTK)
        testlib.safe_unlink('vtkfile')
        
        e6 = test_image()
        e6.write_image('pgmfile', 0, EMUtil.ImageType.IMAGE_PGM)
        self.assertEqual(EMUtil.get_image_type('pgmfile'), EMUtil.ImageType.IMAGE_PGM)
        testlib.safe_unlink('pgmfile')
        
        e7 = test_image()
        e7.write_image('icosfile', 0, EMUtil.ImageType.IMAGE_ICOS)
        self.assertEqual(EMUtil.get_image_type('icosfile'), EMUtil.ImageType.IMAGE_ICOS)
        testlib.safe_unlink('icosfile')
        
        e8 = test_image()
        e8.write_image('xplorfile', 0, EMUtil.ImageType.IMAGE_XPLOR)
        self.assertEqual(EMUtil.get_image_type('xplorfile'), EMUtil.ImageType.IMAGE_XPLOR)
        testlib.safe_unlink('xplorfile')
        
        e9 = test_image()
        e9.write_image('emfile', 0, EMUtil.ImageType.IMAGE_EM)
        self.assertEqual(EMUtil.get_image_type('emfile'), EMUtil.ImageType.IMAGE_EM)
        testlib.safe_unlink('emfile')
        
        e10 = test_image()
        e10.write_image('imagicfile', 0, EMUtil.ImageType.IMAGE_IMAGIC)
        self.assertEqual(EMUtil.get_image_type('imagicfile.hed'), EMUtil.ImageType.IMAGE_IMAGIC)
        testlib.safe_unlink('imagicfile.hed')
        testlib.safe_unlink('imagicfile.img')
        
    def test_get_image_count(self):
        """test get_image_count function ...................."""
        e = test_image()
        e.write_image('test_image_count.mrc')
        self.assertEqual(EMUtil.get_image_count('test_image_count.mrc'), 1)
        testlib.safe_unlink('test_image_count.mrc')
        
        f = test_image()
        f.write_image('test_image_count.hdf', 0)
        f.write_image('test_image_count.hdf', 1)
        f.write_image('test_image_count.hdf', 2)
        self.assertEqual(EMUtil.get_image_count('test_image_count.hdf'), 3)
        testlib.safe_unlink('test_image_count.hdf')
        
        g = test_image()
        g.write_image('test_image_count.spi', 0, EMUtil.ImageType.IMAGE_SPIDER)
        g.write_image('test_image_count.spi', 1, EMUtil.ImageType.IMAGE_SPIDER)
        g.write_image('test_image_count.spi', 2, EMUtil.ImageType.IMAGE_SPIDER)
        self.assertEqual(EMUtil.get_image_count('test_image_count.spi'), 3)
        testlib.safe_unlink('test_image_count.spi')
        
    def test_get_imagetype_name(self):
        """test get_imagetype_name function ................."""
        self.assertEqual(EMUtil.get_imagetype_name(EMUtil.ImageType.IMAGE_V4L), 'V4L2')
        self.assertEqual(EMUtil.get_imagetype_name(EMUtil.ImageType.IMAGE_MRC), 'MRC')
        self.assertEqual(EMUtil.get_imagetype_name(EMUtil.ImageType.IMAGE_SPIDER), 'SPIDER')
        self.assertEqual(EMUtil.get_imagetype_name(EMUtil.ImageType.IMAGE_SINGLE_SPIDER), 'Single-SPIDER')
        self.assertEqual(EMUtil.get_imagetype_name(EMUtil.ImageType.IMAGE_IMAGIC), 'IMAGIC')
        self.assertEqual(EMUtil.get_imagetype_name(EMUtil.ImageType.IMAGE_PGM), 'PGM')
        self.assertEqual(EMUtil.get_imagetype_name(EMUtil.ImageType.IMAGE_LST), 'LST')
        self.assertEqual(EMUtil.get_imagetype_name(EMUtil.ImageType.IMAGE_PIF), 'PIF')
        self.assertEqual(EMUtil.get_imagetype_name(EMUtil.ImageType.IMAGE_PNG), 'PNG')
        self.assertEqual(EMUtil.get_imagetype_name(EMUtil.ImageType.IMAGE_HDF), 'HDF5')
        self.assertEqual(EMUtil.get_imagetype_name(EMUtil.ImageType.IMAGE_DM3), 'GatanDM3')
        self.assertEqual(EMUtil.get_imagetype_name(EMUtil.ImageType.IMAGE_DM3), 'GatanDM4')
        self.assertEqual(EMUtil.get_imagetype_name(EMUtil.ImageType.IMAGE_TIFF), 'TIFF')
        self.assertEqual(EMUtil.get_imagetype_name(EMUtil.ImageType.IMAGE_VTK), 'VTK')
        self.assertEqual(EMUtil.get_imagetype_name(EMUtil.ImageType.IMAGE_SAL), 'HDR')
        self.assertEqual(EMUtil.get_imagetype_name(EMUtil.ImageType.IMAGE_ICOS), 'ICOS_MAP')
        self.assertEqual(EMUtil.get_imagetype_name(EMUtil.ImageType.IMAGE_EMIM), 'EMIM')
        self.assertEqual(EMUtil.get_imagetype_name(EMUtil.ImageType.IMAGE_GATAN2), 'GatanDM2')
        self.assertEqual(EMUtil.get_imagetype_name(EMUtil.ImageType.IMAGE_JPEG), 'JPEG')
        self.assertEqual(EMUtil.get_imagetype_name(EMUtil.ImageType.IMAGE_AMIRA), 'AmiraMesh')
        self.assertEqual(EMUtil.get_imagetype_name(EMUtil.ImageType.IMAGE_XPLOR), 'XPLOR')
        self.assertEqual(EMUtil.get_imagetype_name(EMUtil.ImageType.IMAGE_EM), 'EM')
        self.assertEqual(EMUtil.get_imagetype_name(EMUtil.ImageType.IMAGE_UNKNOWN), 'unknown')
        
    def test_get_datatype_string(self):
        """test get_datatype_string function ................"""
        self.assertEqual(EMUtil.get_datatype_string(EMUtil.EMDataType.EM_CHAR), 'CHAR')
        self.assertEqual(EMUtil.get_datatype_string(EMUtil.EMDataType.EM_UCHAR), 'UNSIGNED CHAR')
        self.assertEqual(EMUtil.get_datatype_string(EMUtil.EMDataType.EM_SHORT), 'SHORT')
        self.assertEqual(EMUtil.get_datatype_string(EMUtil.EMDataType.EM_USHORT), 'UNSIGNED SHORT')
        self.assertEqual(EMUtil.get_datatype_string(EMUtil.EMDataType.EM_INT), 'INT')
        self.assertEqual(EMUtil.get_datatype_string(EMUtil.EMDataType.EM_UINT), 'UNSIGNED INT')
        self.assertEqual(EMUtil.get_datatype_string(EMUtil.EMDataType.EM_FLOAT), 'FLOAT')
        self.assertEqual(EMUtil.get_datatype_string(EMUtil.EMDataType.EM_DOUBLE), 'DOUBLE')
        self.assertEqual(EMUtil.get_datatype_string(EMUtil.EMDataType.EM_SHORT_COMPLEX), 'SHORT_COMPLEX')
        self.assertEqual(EMUtil.get_datatype_string(EMUtil.EMDataType.EM_USHORT_COMPLEX), 'USHORT_COMPLEX')
        self.assertEqual(EMUtil.get_datatype_string(EMUtil.EMDataType.EM_FLOAT_COMPLEX), 'FLOAT_COMPLEX')
        self.assertEqual(EMUtil.get_datatype_string(EMUtil.EMDataType.EM_UNKNOWN), 'UNKNOWN')
        
    def test_read_write_HDF_attribute(self):
        """test read/write single attribute from/to HDF5 file"""
        file = 'test.hdf'
        img = test_image()
        img.write_image(file, -1)
        img.write_image(file, -1)
        
        EMUtil.write_hdf_attribute(file, 'count', 100)
        EMUtil.write_hdf_attribute(file, 'count', 1000, 1)
        c100 = EMUtil.read_hdf_attribute(file, 'count')
        c1000 = EMUtil.read_hdf_attribute(file, 'count', 1)
        
        self.assertEqual(EMUtil.read_hdf_attribute(file, 'count'), 100)
        self.assertEqual(EMUtil.read_hdf_attribute(file, 'count', 1), 1000)
        
        EMUtil.delete_hdf_attribute(file, 'count')
        d = img.get_attr_dict()
        self.assertEqual(d.has_key('count'), False)
        
        testlib.safe_unlink(file)
        
def test_main():
    p = OptionParser()
    p.add_option('--t', action='store_true', help='test exception', default=False )
    global IS_TEST_EXCEPTION
    opt, args = p.parse_args()
    if opt.t:
        IS_TEST_EXCEPTION = True
    Log.logger().set_level(-1)  #perfect solution for quenching the Log error information, thank Liwei
    suite1 = unittest.TestLoader().loadTestsFromTestCase(TestUtils)
    suite2 = unittest.TestLoader().loadTestsFromTestCase(TestEMUtils)
    unittest.TextTestRunner(verbosity=2).run(suite1)
    unittest.TextTestRunner(verbosity=2).run(suite2)

if __name__ == '__main__':
    test_main()
