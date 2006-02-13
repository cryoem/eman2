#!/bin/env python

from EMAN2 import *
import unittest
from test import test_support
import os
import random
import math
import testlib

class TestUtils(unittest.TestCase):
    """test Util class"""
    
    def test_is_file_exist(self):
        """test is_file_exist() function ...................."""
        imgfile1 = "test_is_file_exist.mrc"
        TestUtil.make_image_file(imgfile1, MRC)
        result1 = Util.is_file_exist(imgfile1)
        self.assertEqual(result1, True)
        
        result2 = Util.is_file_exist("__nosuchafile__.dm3")
        self.assertEqual(result2, False)
        
        result3 = Util.is_file_exist("")
        self.assertEqual(result3, False)

        os.unlink(imgfile1)

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

        b2 = Util.sbasename("/tmp/hello.mrc")
        self.assertEqual(b2, "hello.mrc")

        b3 = Util.sbasename("./test/hello")
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
   
    def test_voea(self):
       """test voea() function ............................."""
       v = Util.voea(1.2)    #test default argument
       v2 = Util.voea(1.2, 20.0, 3.0, 401.0)
        
    #quadri() function changed, return image as EMData* in argument
    def no_test_quadri(self):
        """test quadri() function ..........................."""
        e = EMData()
        e.set_size(32,32,32)
        e.process('testimage.noise.uniform.rand')
        
        #f = Util.quadri(e, 2.0, 3.0)    #test default argument
        f2 = Util.quadri(e, 2.3, 3.4, 2)    #test non-default argument
        
def test_main():
    test_support.run_unittest(TestUtils)

if __name__ == '__main__':
    test_main()


