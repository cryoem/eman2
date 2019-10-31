"""
#!/usr/bin/env python
from __future__ import print_function

# Author: Adnan Ali 2019 (adnan.ali@mpi-dortmund.mpg.de)
# Author: Luca Lusnig 2019 (luca.lusnig@mpi-dortmund.mpg.de)
# Author: Markus Stabrin 2019 (markus.stabrin@mpi-dortmund.mpg.de)
# Author: Fabian Schoenfeld 2019 (fabian.schoenfeld@mpi-dortmund.mpg.de)
# Author: Thorsten Wagner 2019 (thorsten.wagner@mpi-dortmund.mpg.de)
# Author: Tapu Shaikh 2019 (tapu.shaikh@mpi-dortmund.mpg.de)
# Copyright (c) 2019 Max Planck Institute of Molecular Physiology
# Copyright (c) 2000-2006 The University of Texas - Houston Medical School
#
# This software is issued under a joint BSD/GNU license. You may use the
# source code in this file under either license. However, note that the
# complete EMAN2 and SPARX software packages have some GPL dependencies,
# so you are responsible for compliance with the licenses of these packages
# if you opt to use BSD licensing. The warranty disclaimer below holfds
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
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
#

from builtins import range
import unittest
from optparse import OptionParser
from EMAN2 import *

IS_TEST_EXCEPTION = False

# ====================================================================================================================
class TestMultiAlignStability(unittest.TestCase):
    #this is unit test for multi_align_stability(...) from pixel_error.py

    # ======================= TEST FOR 2 GROUPS OF TRANSFORMATIONS (L=2, N>=1)
    
    # a0, a1 - list of parameters (alpha0, sx0, sy0, mirror0, aplha1, sx1, sy1, mirror1, alpha2, ...)
    # d - diameter
    def internal_test_conf(self, a1, a2, d):
        from pixel_error import align_diff_params, multi_align_stability
        
        # number of images
        N = len(a1) / 4
        
        # ideal G matrices (we assumed that the last one always equals 0)
        alpha, sx, sy, mirror = align_diff_params(a1, a2)  # I am not sure that it is a global solution for L=2...
        G1 = Transform({"type":"2D","alpha":alpha,"tx":sx ,"ty":sy ,"mirror":mirror,"scale":1.0})
        G2 = Transform({"type":"2D","alpha":0.0  ,"tx":0.0,"ty":0.0,"mirror":0     ,"scale":1.0})
        
        #print "G1=", G1
        
        # ideal H matrices and pixel errors
        H = []
        pixel_error = []
        for i in range(N):
            T1 = Transform({"type":"2D","alpha":a1[4*i+0],"tx":a1[4*i+1],"ty":a1[4*i+2],"mirror":a1[4*i+3],"scale":1.0})
            T2 = Transform({"type":"2D","alpha":a2[4*i+0],"tx":a2[4*i+1],"ty":a2[4*i+2],"mirror":a2[4*i+3],"scale":1.0})
            GT1 = G1 * T1
            GT1_alpha = GT1.get_params("2D")["alpha"]
            GT1_tx    = GT1.get_params("2D")["tx"   ]
            GT1_ty    = GT1.get_params("2D")["ty"   ]
            GT2 = G2 * T2
            GT2_alpha = GT2.get_params("2D")["alpha"]
            GT2_tx    = GT2.get_params("2D")["tx"   ]
            GT2_ty    = GT2.get_params("2D")["ty"   ]
            # fit period
            while GT1_alpha < GT2_alpha - 180.0:
                GT1_alpha += 360.0
            while GT1_alpha > GT2_alpha + 180.0:
                GT1_alpha -= 360.0
            # H matrix
            H_alpha  = (GT1_alpha + GT2_alpha) / 2
            H_tx     = (GT1_tx    + GT2_tx   ) / 2
            H_ty     = (GT1_ty    + GT2_ty   ) / 2
            H_mirror = GT1.get_params("2D")["mirror"]
            self.assertEqual( H_mirror, GT2.get_params("2D")["mirror"] )
            H.append( Transform({"type":"2D","alpha":H_alpha,"tx":H_tx,"ty":H_ty,"mirror":H_mirror,"scale":1.0}) )
            #pixel error
            sum_sin = sin( GT1_alpha * pi / 180.0 ) + sin( GT2_alpha * pi / 180.0 )
            sum_cos = cos( GT1_alpha * pi / 180.0 ) + cos( GT2_alpha * pi / 180.0 )
            var_sx = (GT1_tx - GT2_tx)**2 / 2
            var_sy = (GT1_ty - GT2_ty)**2 / 2
            squared_pixel_error = (d/2)**2 * (1 - sqrt(sum_sin**2 + sum_cos**2) / 2) + var_sx + var_sy
            pixel_error.append( sqrt(squared_pixel_error) )
        
        # function being tested
        stable_set, mirror_consistent_rate, pix_err = multi_align_stability([a1,a2], err_thld=99999.0, print_individual=False, d=d)
        
        # verification of H matrices and pixel errors
        self.assertEqual( len(stable_set), N )
        for sse in stable_set:
            pixerr = sse[0]
            i      = sse[1]
            alpha  = sse[2][0]
            sx     = sse[2][1]
            sy     = sse[2][2]
            mirror = sse[2][3]
            # fit alpha into proper period
            while alpha < H[i].get_params("2D")["alpha"] - 180.0:
                alpha += 360.0
            while alpha > H[i].get_params("2D")["alpha"] + 180.0:
                alpha -= 360.0
            # allowed errors
            allowed_mismatch_angle  = 5.0
            allowed_mismatch_shift  = 0.5
            allowed_mismatch_pixerr = max( 0.1, (0.05*pixel_error[i]) )  # error <= 5%
            # validation
            self.assertAlmostEqual( alpha , H[i].get_params("2D")["alpha" ], delta=allowed_mismatch_angle )
            self.assertAlmostEqual( sx    , H[i].get_params("2D")["tx"    ], delta=allowed_mismatch_shift )
            self.assertAlmostEqual( sy    , H[i].get_params("2D")["ty"    ], delta=allowed_mismatch_shift )
            self.assertEqual      ( mirror, H[i].get_params("2D")["mirror"])
            #self.assertAlmostEqual( pixerr, pixel_error[i]                 , delta=allowed_mismatch_pixerr )  


    def internal_test_conf_without_mirror(self, T1, T2):
        diameters = [64]#, 75, 100, 128]
        N = len(T1) / 3
        a1 = []
        a2 = []
        for i in range(N):
            a1.append( T1[3*i+0] )
            a1.append( T1[3*i+1] )
            a1.append( T1[3*i+2] )
            a1.append( 0 )
            a2.append( T2[3*i+0] )
            a2.append( T2[3*i+1] )
            a2.append( T2[3*i+2] )
            a2.append( 0 )
        for d in diameters:
            self.internal_test_conf(a1, a2, d)

    def internal_test_conf_with_mirror(self, T1, T2):
        diameters = [64, 75, 100, 128]
        N = len(T1) / 3
        a1 = []
        a2 = []
        a1m = []
        a2m = []
        for i in range(N):
            a1.append( T1[3*i+0] )
            a1.append( T1[3*i+1] )
            a1.append( T1[3*i+2] )
            a1.append( 0 )
            a1m.append( T1[3*i+0] )
            a1m.append( T1[3*i+1] )
            a1m.append( T1[3*i+2] )
            a1m.append( 1 )    
            a2.append( T2[3*i+0] )
            a2.append( T2[3*i+1] )
            a2.append( T2[3*i+2] )
            a2.append( 0 )
            a2m.append( T2[3*i+0] )
            a2m.append( T2[3*i+1] )
            a2m.append( T2[3*i+2] )
            a2m.append( 1 )
        for d in diameters:
            self.internal_test_conf(a1m, a2, d)
            self.internal_test_conf(a1, a2m, d)
            self.internal_test_conf(a1m, a2m, d)


    # ======================= TESTS FOR 2 TRANSFORMATIONS (L=2, N=2)

    def test_two_images_without_mirrors(self):
        #test: L=2,  N=2,  without mirror ..........................
        a1 = [23.1, 0.0, 0.0, 23.1, 0.0, 0.0]
        a2 = [72.2, 0.0, 0.0, 72.2, 0.0, 0.0]
        self.internal_test_conf_without_mirror(a1, a2)
        a1 = [23.1, 2.0, -3.0, 23.1, 2.0, -3.0]
        a2 = [72.2, 1.0,  4.0, 72.2, 1.0,  4.0]
        self.internal_test_conf_without_mirror(a1, a2)
        a1 = [193.1,  0.0, -2.0,  23.1, 2.0, -3.0]
        a2 = [ 72.2, -1.0,  4.0, 144.2, 3.0, -1.0]
        self.internal_test_conf_without_mirror(a1, a2)

    def test_two_images_with_mirrors(self):
        #test: L=2,  N=2,  with mirror ..........................
        a1 = [23.1, 0.0, 0.0, 23.1, 0.0, 0.0]
        a2 = [72.2, 0.0, 0.0, 72.2, 0.0, 0.0]
        self.internal_test_conf_with_mirror(a1, a2)
        a1 = [23.1, 2.0, -3.0, 23.1, 2.0, -3.0]
        a2 = [72.2, 1.0,  4.0, 72.2, 1.0,  4.0]
        self.internal_test_conf_with_mirror(a1, a2)
        a1 = [193.1,  0.0, -2.0,  23.1, 2.0, -3.0]
        a2 = [ 72.2, -1.0,  4.0, 144.2, 3.0, -1.0]
        self.internal_test_conf_with_mirror(a1, a2)

    # ======================= TESTS FOR 4 TRANSFORMATIONS (L=2, N=4)
    
    def test_four_images_without_mirrors(self):
        # test: L=2,  N=4,  without mirror .......................... 
        a1 = [193.0, 0.0, 0.0,   122.0, 0.0, 0.0,   223.0, 0.0, 0.0,   213.0, 0.0, 0.0]
        a2 = [ 72.0, 0.0, 0.0,   164.0, 0.0, 0.0,   323.0, 0.0, 0.0,   117.0, 0.0, 0.0]
        self.internal_test_conf_without_mirror(a1, a2)
        a1 = [193.1,  0.0, -2.0,   122.1, 3.0, -3.0,   76.1, 5.0, -1.0,   291.8, 0.2, 0.0]
        a2 = [ 72.2, -1.0,  4.0,   164.2, 3.0, -1.0,   56.1, 1.0, -7.0,   193.1, 0.3, 1.0]
        self.internal_test_conf_without_mirror(a1, a2)

    def test_four_images_with_mirrors(self):
        # test: L=2,  N=4,  with mirror .......................... 
        a1 = [193.0, 0.0, 0.0,   122.0, 0.0, 0.0,   223.0, 0.0, 0.0,   213.0, 0.0, 0.0]
        a2 = [ 72.0, 0.0, 0.0,   164.0, 0.0, 0.0,   323.0, 0.0, 0.0,   117.0, 0.0, 0.0]
        self.internal_test_conf_with_mirror(a1, a2)
        a1 = [193.1,  0.0, -2.0,   122.1, 3.0, -3.0,   76.1, 5.0, -1.0,   291.8, 0.2, 0.0]
        a2 = [ 72.2, -1.0,  4.0,   164.2, 3.0, -1.0,   56.1, 1.0, -7.0,   193.1, 0.3, 1.0]
        self.internal_test_conf_with_mirror(a1, a2)

    # ======================= TESTS FOR 10 TRANSFORMATIONS (L=2, N=10)
    
    def test_ten_images_without_mirrors(self):
        # test: L=2,  N=10,  without mirror .......................... 
        from random import random
        a1 = [193.1,  0.0, -2.0,   122.1, 3.0, -3.0,    44.6, 6.0, -3.0,    23.1, 2.0,  1.0,    93.1, -2.0,  4.2,   223.1, 2.0,  2.1,    52.1,  5.1, -3.0,   213.3, 2.0, -1.2,    76.1, 5.0, -1.0,   291.8, 0.2, 0.0]
        a2 = [ 72.2, -1.0,  4.0,   164.2, 3.0, -1.0,   196.1, 2.0, -3.0,   239.1, 2.0, -3.0,   356.2,  3.1,  2.1,   323.1, 0.2, -0.8,   303.1, -2.1,  2.0,   117.6, 2.0,  2.0,    56.1, 1.0, -7.0,   193.1, 0.3, 1.0]
        self.internal_test_conf_without_mirror(a1, a2)

    def test_ten_images_with_mirrors(self):
        # test: L=2,  N=10,  with mirror .......................... 
        from random import random
        a1 = [193.1,  0.0, -2.0,   122.1, 3.0, -3.0,    44.6, 6.0, -3.0,    23.1, 2.0,  1.0,    93.1, -2.0,  4.2,   223.1, 2.0,  2.1,    52.1,  5.1, -3.0,   213.3, 2.0, -1.2,    76.1, 5.0, -1.0,   291.8, 0.2, 0.0]
        a2 = [ 72.2, -1.0,  4.0,   164.2, 3.0, -1.0,   196.1, 2.0, -3.0,   239.1, 2.0, -3.0,   356.2,  3.1,  2.1,   323.1, 0.2, -0.8,   303.1, -2.1,  2.0,   117.6, 2.0,  2.0,    56.1, 1.0, -7.0,   193.1, 0.3, 1.0]
        self.internal_test_conf_with_mirror(a1, a2)


def test_main():
    from EMAN2 import Log
    p = OptionParser()
    p.add_option('--t', action='store_true', help='test exception', default=False )
    global IS_TEST_EXCEPTION
    opt, args = p.parse_args()
    if opt.t:
        IS_TEST_EXCEPTION = True
    Log.logger().set_level(-1)  #perfect solution for quenching the Log error information, thank Liwei
    suite = unittest.TestLoader().loadTestsFromTestCase(TestMultiAlignStability)
    unittest.TextTestRunner(verbosity=2).run(suite)

if __name__ == '__main__':
    test_main()
"""

from __future__ import print_function
from __future__ import division

import unittest
from os import path
from mpi import *

from ..libpy import sp_global_def

mpi_init(0, [])
sp_global_def.BATCH = True
sp_global_def.MPI = True

ABSOLUTE_PATH = path.dirname(path.realpath(__file__))

from EMAN2 import Transform
from numpy import array_equal as numpy_array_equal
from test_module import   IMAGE_2D,IMAGE_3D,IMAGE_2D_REFERENCE

from sphire.libpy import sp_pixel_error as oldfu
from ..utils.SPHIRE.libpy import sp_pixel_error as fu

"""
pickle files stored under smb://billy.storage.mpi-dortmund.mpg.de/abt3/group/agraunser/transfer/Adnan/pickle files
"""

"""
WHAT IS MISSING:
0) in all the cases where the input file is an image. I did not test the case with a complex image. I was not able to generate it 


RESULT AND KNOWN ISSUES
Some compatibility tests for the following functions fail!!!
1) 

In these tests there is a bug --> syntax error:
1) 

In these tests there is a strange behavior:
1) 
"""






#   THESE FUNCTIONS ARE COMMENTED BECAUSE NOT INVOLVED IN THE PYTHON3 CONVERSION. THEY HAVE TO BE TESTED ANYWAY
""" start: new in sphire 1.3
class Test_angle_error(unittest.TestCase):

    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.angle_error()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.angle_error()
        self.assertEqual(str(cm_new.exception), "angle_error() takes at least 2 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_empty_ang1_returns_IndexError_list_index_out_of_range(self):
        with self.assertRaises(IndexError) as cm_new:
            fu.angle_error(ang1=[], ang2=[0.6804220676422119, 0.6526213884353638, 0.3333333432674408], delta_ang=0.1)
        with self.assertRaises(IndexError) as cm_old:
            oldfu.angle_error(ang1=[], ang2=[0.6804220676422119, 0.6526213884353638, 0.3333333432674408], delta_ang=0.1)
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_empty_ang2_returns_IndexError_list_index_out_of_range(self):
        with self.assertRaises(IndexError) as cm_new:
            fu.angle_error(ang2=[], ang1=[0.6804220676422119, 0.6526213884353638, 0.3333333432674408], delta_ang=0.1)
        with self.assertRaises(IndexError) as cm_old:
            oldfu.angle_error(ang2=[], ang1=[0.6804220676422119, 0.6526213884353638, 0.3333333432674408], delta_ang=0.1)
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_angle_error(self):
        return_new =fu.angle_error(ang1=[0.6804220676422119, 0.6526213884353638, 0.3333333432674408], ang2=[0.6804220676422119, 0.65262184353638, 0.3333333432674408], delta_ang=0.1)
        return_old =oldfu.angle_error(ang1=[0.6804220676422119, 0.6526213884353638, 0.3333333432674408], ang2=[0.6804220676422119, 0.65262184353638, 0.3333333432674408], delta_ang=0.1)
        self.assertEqual(return_new,return_old)
        self.assertEqual(return_new, 2.99999543075)


#todo: find a new image -->RuntimeError: NotExistingObjectException at /home/lusnig/src_sphire_1_3/eman2/libEM/emdata_metadata.cpp:1203: error with 'xform.align2d': 'The requested key does not exist' caught
class Test_align_diff(unittest.TestCase):
    data = [IMAGE_2D, IMAGE_2D_REFERENCE]
    data2 = [IMAGE_3D, IMAGE_3D]

    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.align_diff()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.align_diff()
        self.assertEqual(str(cm_new.exception), "align_diff() takes at least 1 argument (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_empty_data_returns_IndexError_list_index_out_of_range(self):
        with self.assertRaises(IndexError) as cm_new:
            fu.align_diff(data1=[], data2=None, suffix="_ideal")
        with self.assertRaises(IndexError) as cm_old:
            oldfu.align_diff(data1=[], data2=None, suffix="_ideal")
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_empty_data2_returns_IndexError_list_index_out_of_range(self):
        with self.assertRaises(IndexError) as cm_new:
            fu.align_diff(data1=self.data, data2=[], suffix="_ideal")
        with self.assertRaises(IndexError) as cm_old:
            oldfu.align_diff(data1=self.data, data2=[], suffix="_ideal")
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_align_diff(self):
        return_new = fu.align_diff(data1=self.data, data2=None, suffix="_ideal")
        return_old = oldfu.align_diff(data1=self.data, data2=None, suffix="_ideal")
        self.assertTrue(numpy_array_equal(return_new,return_old))

    def test_align_diff_3Dimgs_error(self):
        return_new = fu.align_diff(data1=self.data2, data2=None, suffix="_ideal")
        return_old = oldfu.align_diff(data1=self.data2, data2=None, suffix="_ideal")
        self.assertTrue(numpy_array_equal(return_new,return_old))

    def test_align_diff_with_data2(self):
        return_new = fu.align_diff(data1=self.data, data2= [IMAGE_2D, IMAGE_2D], suffix="_ideal")
        return_old = oldfu.align_diff(data1=self.data, data2=[IMAGE_2D, IMAGE_2D], suffix="_ideal")
        self.assertTrue(numpy_array_equal(return_new,return_old))


    def test_data_different_length(self):
        return_new = fu.align_diff(data1=self.data, data2=[IMAGE_2D_REFERENCE], suffix="_ideal")
        return_old = oldfu.align_diff(data1=self.data, data2=[IMAGE_2D_REFERENCE], suffix="_ideal")
        self.assertTrue(numpy_array_equal(return_new,return_old))
        self.assertTrue(numpy_array_equal(return_new, [0.0, 0.0, 0.0, 0]))


#todo: need text file
class Test_align_diff_textfile(unittest.TestCase):

    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.align_diff_textfile()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.align_diff_textfile()
        self.assertEqual(str(cm_new.exception), "align_diff_textfile() takes exactly 2 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_align_diff_textfile(self):
        oldv = oldfu.align_diff_textfile(textfile1=0, textfile2=0)
        v = fu.align_diff_textfile(textfile1=0, textfile2=0)
        pass


#todo: find a new image -->RuntimeError: NotExistingObjectException at /home/lusnig/src_sphire_1_3/eman2/libEM/emdata_metadata.cpp:1203: error with 'xform.align2d': 'The requested key does not exist' caught
class Test_ave_ali_err(unittest.TestCase):
    data = [IMAGE_2D, IMAGE_2D_REFERENCE]
    data2 = [IMAGE_3D, IMAGE_3D]
    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.ave_ali_err()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.ave_ali_err()
        self.assertEqual(str(cm_new.exception), "ave_ali_err() takes at least 1 argument (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_empty_data_returns_IndexError_list_index_out_of_range(self):
        with self.assertRaises(IndexError) as cm_new:
            fu.ave_ali_err(data1=[], data2=None, suffix="_ideal", r=25)
        with self.assertRaises(IndexError) as cm_old:
            oldfu.ave_ali_err(data1=[], data2=None, suffix="_ideal", r=25)
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_empty_data2_returns_IndexError_list_index_out_of_range(self):
        with self.assertRaises(IndexError) as cm_new:
            fu.ave_ali_err(data1=self.data, data2=[], suffix="_ideal", r=25)
        with self.assertRaises(IndexError) as cm_old:
            oldfu.ave_ali_err(data1=self.data, data2=[], suffix="_ideal", r=25)
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_ave_ali_err(self):
        return_new = fu.ave_ali_err(data1=self.data, data2=None, suffix="_ideal", r=25)
        return_old = oldfu.ave_ali_err(data1=self.data, data2=None, suffix="_ideal", r=25)
        self.assertTrue(numpy_array_equal(return_new, return_old))

    def test_ave_ali_err_3Dimgs_error(self):
        return_new = fu.ave_ali_err(data1=self.data2, data2=None, suffix="_ideal", r=25)
        return_old = oldfu.ave_ali_err(data1=self.data2, data2=None, suffix="_ideal", r=25)
        self.assertTrue(numpy_array_equal(return_new, return_old))

    def test_ave_ali_err_with_data2(self):
        return_new = fu.ave_ali_err(data1=self.data, data2=[IMAGE_2D, IMAGE_2D], suffix="_ideal", r=25)
        return_old = oldfu.ave_ali_err(data1=self.data, data2=[IMAGE_2D, IMAGE_2D], suffix="_ideal", r=25)
        self.assertTrue(numpy_array_equal(return_new, return_old))

    def test_data_different_length(self):
        return_new = fu.ave_ali_err(data1=self.data, data2=[IMAGE_2D_REFERENCE], suffix="_ideal", r=25)
        return_old = oldfu.ave_ali_err(data1=self.data, data2=[IMAGE_2D_REFERENCE], suffix="_ideal", r=25)
        self.assertTrue(numpy_array_equal(return_new, return_old))
        self.assertTrue(numpy_array_equal(return_new, [0.0, 0.0, 0.0, 0]))



#todo: need text file
class Test_ave_ali_err_textfile(unittest.TestCase):
    def test_ave_ali_err_textfile(self):
        oldv = oldfu.ave_ali_err_textfile(textfile1=0, textfile2=0,r=25)
        v = fu.ave_ali_err_textfile(textfile1=0, textfile2=0, r=25)
        pass



class Test_multi_align_diff_params(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.multi_align_diff_params()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.multi_align_diff_params()
        self.assertEqual(str(cm_new.exception), "multi_align_diff_params() takes at least 1 argument (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_empty_ali_params_returns_IndexError_list_index_out_of_range(self):
        with self.assertRaises(IndexError) as cm_new:
            fu.multi_align_diff_params(ali_params=[], verbose=0)
        with self.assertRaises(IndexError) as cm_old:
            oldfu.multi_align_diff_params(ali_params=[], verbose=0)
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_empty_ali_params0_returns_IndexError_list_index_out_of_range(self):
        with self.assertRaises(IndexError) as cm_new:
            fu.multi_align_diff_params(ali_params=[[],[1,2,3]], verbose=0)
        with self.assertRaises(IndexError) as cm_old:
            oldfu.multi_align_diff_params(ali_params=[[],[1,2,3]], verbose=0)
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_multi_align_diff_params(self):
        return_new = fu.multi_align_diff_params(ali_params=[[13, 23, 31,0,1, 23, 31,0,13, 3, 31,0], [13, 23, 31,0,0,1, 2, 3,13, 23, 31,0]], verbose=0)
        return_old = oldfu.multi_align_diff_params(ali_params=[[13, 23, 31,0,1, 23, 31,0,13, 3, 31,0], [13, 23, 31,0,0,1, 2, 3,13, 23, 31,0]], verbose=0)
        self.assertTrue(numpy_array_equal(return_new, return_old))
        self.assertTrue(numpy_array_equal(return_new, [[100.00000000000001, 0.6666666666666666, 0, 1, 0.0, 10.0, 0.0, 0]]))




class Test_calc_connect_list(unittest.TestCase):
    multi_align_results_new=fu.multi_align_diff_params(ali_params=[[13, 23, 31, 0, 1, 23, 31, 0, 13, 3, 31, 0], [13, 23, 31, 0, 0, 1, 2, 3, 13, 23, 31, 0]], verbose=0)
    multi_align_results_old = oldfu.multi_align_diff_params(ali_params=[[13, 23, 31, 0, 1, 23, 31, 0, 13, 3, 31, 0], [13, 23, 31, 0, 0, 1, 2, 3, 13, 23, 31, 0]], verbose=0)

    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.calc_connect_list()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.calc_connect_list()
        self.assertEqual(str(cm_new.exception), "calc_connect_list() takes at least 1 argument (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_calc_connect_list(self):
        return_old = oldfu.calc_connect_list(multi_align_results=self.multi_align_results_old, pixel_error_threshold = 5.0, mirror_consistency_threshold = 0.8)
        return_new = fu.calc_connect_list(multi_align_results=self.multi_align_results_new, pixel_error_threshold = 5.0, mirror_consistency_threshold = 0.8)
        self.assertTrue(numpy_array_equal(return_old[0], return_new[0]))
        self.assertEqual(return_old[1], return_new[1])
        self.assertTrue(numpy_array_equal([], return_new[0]))
        self.assertEqual(0, return_new[1])

    def test_empty_multi_align_results_returns_IndexError_list_index_out_of_range(self):
        with self.assertRaises(IndexError) as cm_new:
            fu.calc_connect_list(multi_align_results=[], pixel_error_threshold = 5.0, mirror_consistency_threshold = 0.8)
        with self.assertRaises(IndexError) as cm_old:
            oldfu.calc_connect_list(multi_align_results=[], pixel_error_threshold = 5.0, mirror_consistency_threshold = 0.8)
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_empty_multi_align_results0_returns_IndexError_list_index_out_of_range(self):
        with self.assertRaises(IndexError) as cm_new:
            fu.calc_connect_list(multi_align_results=[[], [1, 2, 3]], pixel_error_threshold = 5.0, mirror_consistency_threshold = 0.8)
        with self.assertRaises(IndexError) as cm_old:
            oldfu.calc_connect_list(multi_align_results=[[], [1, 2, 3]], pixel_error_threshold = 5.0, mirror_consistency_threshold = 0.8)
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))



#i do not test the consistency of the ali_params1 and ali_params2 params because they are used in the 'align_diff_params' call ... I tested them in 'Test_align_diff_params
class Test_ali_stable_list(unittest.TestCase):
    ali_params1 = [358.1422162327943, 1.5746676921844482, 2.359189510345459, 1, 335.90588778212316, 0.22259289026260376, 1.6781892776489258, 0, 81.05934436345304, 2.0682716369628906, 2.3159945011138916, 1, 88.94114972876307, 2.3903369903564453, 2.760911226272583, 0, 230.3400198158369, -1.2752182483673096, -1.0786179304122925, 0, 179.6066808748772, 0.3892550766468048, -0.7324743270874023, 0, 179.05318018749753, -0.694373607635498, -0.10197801142930984, 1, 267.2248149400371, 0.5244513750076294, 0.7750062346458435, 1, 168.26126969747668, 1.208488941192627, 1.3360750675201416, 0, 177.5181948630324, 1.0620733499526978, -1.5996335744857788, 0, 90.44012904554053, 0.8572239875793457, -1.572252869606018, 0]
    ali_params2 = [358.142216327943, 1.5746676921844482, 2.359189510345459, 1, 335.90588778212316, 0.22259289026260376, 1.6781892776489258, 0, 81.05934436345304, 2.0682716369628906, 2.3159945011138916, 1, 88.94114972876307, 2.3903369903564453, 2.760911226272583, 0, 230.3400198158369, -1.2752182483673096, -1.0786179304122925, 0, 179.6066808748772, 0.3892550766468048, -0.7324743270874023, 0, 179.05318018749753, -0.694373607635498, -0.10197801142930984, 1, 267.2248149400371, 0.5244513750076294, 0.7750062346458435, 1, 168.26126969747668, 1.208488941192627, 1.3360750675201416, 0, 177.5181948630324, 1.0620733499526978, -1.5996335744857788, 0, 90.44012904554053, 0.8572239875793457, -1.572252869606018, 0]

    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.ali_stable_list()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.ali_stable_list()
        self.assertEqual(str(cm_new.exception), "ali_stable_list() takes at least 3 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_ali_stable_list(self):
        return_new = oldfu.ali_stable_list(ali_params1=self.ali_params1, ali_params2=self.ali_params2, pixel_error_threshold=0.3, r=25)
        return_old = fu.ali_stable_list(ali_params1=self.ali_params1, ali_params2=self.ali_params2, pixel_error_threshold=0.3, r=25)
        self.assertTrue(numpy_array_equal(return_new, return_old))
        self.assertTrue(numpy_array_equal(return_new, [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]))


#todo: look the params
class Test_ave2dtransform(unittest.TestCase):
    ali_params = [358.142216327943, 1.5746676921844482, 2.359189510345459, 1, 335.90588778212316, 0.22259289026260376, 1.6781892776489258, 0, 81.05934436345304, 2.0682716369628906, 2.3159945011138916, 1, 88.94114972876307, 2.3903369903564453, 2.760911226272583, 0, 230.3400198158369, -1.2752182483673096, -1.0786179304122925, 0, 179.6066808748772, 0.3892550766468048, -0.7324743270874023, 0, 179.05318018749753, -0.694373607635498, -0.10197801142930984, 1, 267.2248149400371, 0.5244513750076294, 0.7750062346458435, 1, 168.26126969747668, 1.208488941192627, 1.3360750675201416, 0, 177.5181948630324, 1.0620733499526978, -1.5996335744857788, 0, 90.44012904554053, 0.8572239875793457, -1.572252869606018, 0]
    data=[ali_params,3]
    args=[ali_params, [x * 2 for x in ali_params] , [x /2 for x in ali_params] ]
    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.ave2dtransform()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.ave2dtransform()
        self.assertEqual(str(cm_new.exception), "ave2dtransform() takes at least 2 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_empty_args_returns_IndexError_list_index_out_of_range(self):
        with self.assertRaises(IndexError) as cm_new:
            fu.ave2dtransform(args=[], data=[1,2,3], return_avg_pixel_error = False)
        with self.assertRaises(IndexError) as cm_old:
            oldfu.ave2dtransform(args=[], data=[1,2,3], return_avg_pixel_error = False)
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_empty_args0_returns_IndexError_list_index_out_of_range(self):
        with self.assertRaises(IndexError) as cm_new:
            fu.ave2dtransform(args=[[],self.ali_params], data=[1,2,3], return_avg_pixel_error = False)
        with self.assertRaises(IndexError) as cm_old:
            oldfu.ave2dtransform(args=[], data=[1,2,3], return_avg_pixel_error = False)
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_empty_data_returns_IndexError_list_index_out_of_range(self):
        with self.assertRaises(IndexError) as cm_new:
            fu.ave2dtransform(args=[1, 2, 3], data=[], return_avg_pixel_error = False)
        with self.assertRaises(IndexError) as cm_old:
            oldfu.ave2dtransform(args= [1, 2, 3], data=[], return_avg_pixel_error = False)
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_return_NOavg_pixel_error(self):
        return_new = fu.ave2dtransform(args=self.args, data=self.data, return_avg_pixel_error = False)
        return_old = oldfu.ave2dtransform(args=self.args, data=self.data, return_avg_pixel_error = False)
        self.assertEqual(return_new[0],return_old[0])
        self.assertTrue(numpy_array_equal(return_new[1], return_old[1]))

    def test_return_avg_pixel_error(self):
        return_new = fu.ave2dtransform(args=self.args, data=self.data, return_avg_pixel_error = True)
        return_old = oldfu.ave2dtransform(args=self.args, data=self.data, return_avg_pixel_error = True)
        self.assertEqual(return_new,return_old)



#i do not test the consistency of the agls1 and agls2 params because they are used in the 'sp_utilities.rotation_between_anglesets' call ... I tested them in 'Test_rotation_between_anglesets
class Test_rotate_angleset_to_match(unittest.TestCase):
    agls1 =  [[0.0, 0.0, 1.0], [0.6804220676422119, 0.6526213884353638, 0.3333333432674408], [-0.4104178845882416, 0.8487909436225891, 0.3333333432674408], [-0.9340742230415344, -0.12803982198238373, 0.3333333432674408], [-0.16687190532684326, -0.927923858165741, 0.3333333432674408], [0.8309417366981506, -0.4454488158226013, 0.3333333432674408], [8.742277657347586e-08, 7.64274186065882e-15, -1.0], [0.9340742230415344, 0.12803970277309418, -0.3333333134651184], [0.16687177121639252, 0.927923858165741, -0.3333333134651184], [-0.8309418559074402, 0.44544869661331177, -0.3333333134651184], [-0.6804221272468567, -0.652621328830719, -0.3333333134651184], [0.41041797399520874, -0.8487908840179443, -0.3333333134651184]]
    agls2 = [[0.0, 0.0, 0.66], [0.44907856464385987, 0.4307301163673401, 0.22000000655651095], [-0.27087580382823945, 0.5602020227909088, 0.22000000655651095], [-0.6164889872074127, -0.08450628250837326, 0.22000000655651095], [-0.11013545751571656, -0.6124297463893891, 0.22000000655651095], [0.5484215462207794, -0.2939962184429169, 0.22000000655651095], [5.7699032538494066e-08, 5.044209628034821e-15, -0.66], [0.6164889872074127, 0.08450620383024215, -0.21999998688697817], [0.11013536900281906, 0.6124297463893891, -0.21999998688697817], [-0.5484216248989106, 0.2939961397647858, -0.21999998688697817], [-0.44907860398292543, -0.43073007702827454, -0.21999998688697817], [0.2708758628368378, -0.5602019834518432, -0.21999998688697817]]

    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.rotate_angleset_to_match()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.rotate_angleset_to_match()
        self.assertEqual(str(cm_new.exception), "rotate_angleset_to_match() takes exactly 2 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_return_avg_pixel_error(self):
        return_new = fu.rotate_angleset_to_match(agls1=self.agls1, agls2=self.agls2)
        return_old = oldfu.rotate_angleset_to_match(agls1=self.agls1, agls2=self.agls2)
        self.assertTrue(numpy_array_equal(return_new, return_old))
        self.assertTrue(numpy_array_equal( return_old, [[1.000000002215, 0.0, 0.0], [0.680422069857, 0.652621388435, 0.333333343267], [359.589582117627, 0.848790943622, 0.333333343267], [179.065925779173, 0.128039821983, 180.333333343267], [179.833128096888, 0.927923858166, 180.333333343267], [180.830941738913, 0.445448815823, 180.333333343267], [359.000000089638, 0.0, 0.0], [0.934074225257, 0.128039702772, 359.666666686535], [0.166871773431, 0.927923858166, 359.666666686535], [359.169058146308, 0.445448696613, 359.666666686535], [179.319577874968, 0.652621328831, 179.666666686535], [180.41041797621, 0.848790884018, 179.666666686535]]))



#todo: need stack data
class Test_ordersegments(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.ordersegments()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.ordersegments()
        self.assertEqual(str(cm_new.exception), "ordersegments() takes exactly 2 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_ordersegments(self):
        oldv = oldfu.ordersegments(infilaments=0, ptclcoords=0)
        v = fu.ordersegments(infilaments=0, ptclcoords=0)
        pass



class Test_mapcoords(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.mapcoords()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.mapcoords()
        self.assertEqual(str(cm_new.exception), "ordersegments() takes exactly 4 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_r_lowe_than1(self):
        return_old = oldfu.mapcoords(x=1, y=1, r=0.1, nx=5, ny=5)
        return_new = fu.mapcoords(x=1, y=1, r=0.1, nx=5, ny=5)
        self.assertTrue(numpy_array_equal(return_new, return_old))
        self.assertTrue(numpy_array_equal((0, 0), return_old))

    def test_r_higher_than1(self):
        return_old = oldfu.mapcoords(x=1, y=1, r=1.1, nx=5, ny=5)
        return_new = fu.mapcoords(x=1, y=1, r=1.1, nx=5, ny=5)
        self.assertTrue(numpy_array_equal(return_new, return_old))
        self.assertTrue(numpy_array_equal((0, 0), return_old))


#todo: need stack data
class Test_consistency_params(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.consistency_params()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.consistency_params()
        self.assertEqual(str(cm_new.exception), "ordersegments() takes at least 4 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_consistency_params(self):
        oldv = oldfu.consistency_params(stack="", dphi="", dp="", pixel_size="", phithr=2.5, ythr=1.5, THR=3)
        v = fu.consistency_params(stack="", dphi="", dp="", pixel_size="", phithr=2.5, ythr=1.5, THR=3)
        pass



#todo: need stack data
class Test_getnewhelixcoords(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.getnewhelixcoords()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.getnewhelixcoords()
        self.assertEqual(str(cm_new.exception), "ordersegments() takes at least 5 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_(self):
        oldv = oldfu.getnewhelixcoords(hcoordsname="", outdir="", ratio="",nx="",ny="", newpref="resampled_", boxsize=-1)
        v = fu.getnewhelixcoords(hcoordsname="", outdir="", ratio="",nx="",ny="", newpref="resampled_", boxsize=-1)
        pass


#todo: look into the params
class Test_helical_params_err(unittest.TestCase):
    r=29
    t1 = Transform({"az": 123.155, "alt": 119.102, "phi": 318.812, "tx": -0.00, "ty": -0.00, "tz": -0.00, "mirror": 0,"scale": 1.0000, "type":"eman"})
    t2 = Transform({"az": 146.608, "alt": 143.005, "phi": 146.656, "tx": -0.00, "ty": -0.00, "tz": -0.00, "mirror": 0,"scale": 1.0000, "type":"eman"})
    params2=[t2,t2]
    params1 = [t1,t1]
    fillist = [1,2,3]
    def test_wrong_number_params_returns_TypeError_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.helical_params_err()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.helical_params_err()
        self.assertEqual(str(cm_new.exception), "helical_params_err() takes exactly 3 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))


    def test_params1_emptyList_returns_IndexError(self):
        with self.assertRaises(IndexError) as cm_new:
            fu.helical_params_err(params1= [], params2= self.params2, fil_list=self.fillist)
        with self.assertRaises(IndexError) as cm_old:
            oldfu.helical_params_err(params1= [], params2= self.params2, fil_list=self.fillist)
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_params2_emptyList_returns_IndexError(self):
        with self.assertRaises(IndexError) as cm_new:
            fu.helical_params_err(params1= self.params1, params2= [], fil_list=self.fillist)
        with self.assertRaises(IndexError) as cm_old:
            oldfu.helical_params_err(params1= self.params1, params2= [], fil_list=self.fillist)
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_params1_0_emptyList_returns_IndexError(self):
        with self.assertRaises(IndexError) as cm_new:
            fu.helical_params_err(params1= [[],self.t2], params2= self.params2, fil_list=self.fillist)
        with self.assertRaises(IndexError) as cm_old:
            oldfu.helical_params_err(params1= [[],self.t2], params2= self.params2, fil_list=self.fillist)
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_params2_0_emptyList_returns_IndexError(self):
        with self.assertRaises(IndexError) as cm_new:
            fu.helical_params_err(params1= self.params1, params2= [[],self.t2], fil_list=self.fillist)
        with self.assertRaises(IndexError) as cm_old:
            oldfu.helical_params_err(params1= self.params1, params2= [[],self.t2], fil_list=self.fillist)
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))


    def test_fil_list_emptyList_returns_IndexError(self):
        with self.assertRaises(IndexError) as cm_new:
            fu.helical_params_err(params1= self.params1, params2= self.params2, fil_list=[])
        with self.assertRaises(IndexError) as cm_old:
            oldfu.helical_params_err(params1= self.params1, params2= self.params2, fil_list=[])
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_helical_params_err(self):
        return_new = fu.helical_params_err(params1= self.params1, params2= self.params2, fil_list=self.fillist)
        return_old = oldfu.helical_params_err(params1=self.params1, params2=self.params2, fil_list=self.fillist)
        for i,j in zip(return_new,return_old):
            self.assertTrue(numpy_array_equal(i, j))




 end: new in sphire 1.3"""



class Test_pixel_error_2D(unittest.TestCase):
    ali_params1 = [0.0, 0.0, 0.0]
    ali_params2 = [270.0, -1.0, -1.0]
    r=29

    def test_wrong_number_params_returns_TypeError_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.pixel_error_2D()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.pixel_error_2D()
        self.assertEqual(str(cm_new.exception), "pixel_error_2D() takes at least 2 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_pickle_file(self):
        return_new = fu.pixel_error_2D(ali_params1= self.ali_params1, ali_params2= self.ali_params2, r= self.r)
        return_old = oldfu.pixel_error_2D(ali_params1= self.ali_params1, ali_params2= self.ali_params2, r= self.r)
        self.assertEqual(return_new, return_old)
        self.assertEqual(return_new, 872.25000000000011)

    def test_r_is0(self):
        return_new = fu.pixel_error_2D(ali_params1= self.ali_params1, ali_params2= self.ali_params2, r= 0)
        return_old = oldfu.pixel_error_2D(ali_params1= self.ali_params1, ali_params2= self.ali_params2, r= 0)
        self.assertEqual(return_new, return_old)
        self.assertEqual(return_new, 2.25)


    def test_ali_params1_emptyList_returns_IndexError(self):
        with self.assertRaises(IndexError) as cm_new:
            fu.pixel_error_2D(ali_params1= [], ali_params2= self.ali_params2, r= self.r)
        with self.assertRaises(IndexError) as cm_old:
            oldfu.pixel_error_2D(ali_params1= [], ali_params2= self.ali_params2, r= self.r)
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_ali_params2_emptyList_returns_IndexError(self):
        with self.assertRaises(IndexError) as cm_new:
            fu.pixel_error_2D(ali_params1= self.ali_params1, ali_params2= [], r= self.r)
        with self.assertRaises(IndexError) as cm_old:
            oldfu.pixel_error_2D(ali_params1= self.ali_params1, ali_params2= [], r= self.r)
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_ali_params_both_emptyList_returns_IndexError(self):
        with self.assertRaises(IndexError) as cm_new:
            fu.pixel_error_2D(ali_params1= [], ali_params2= [], r= self.r)
        with self.assertRaises(IndexError) as cm_old:
            oldfu.pixel_error_2D(ali_params1= [], ali_params2= [], r= self.r)
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))


class Test_max_3D_pixel_error(unittest.TestCase):
    r=29
    t1 = Transform({"az": 123.155, "alt": 119.102, "phi": 318.812, "tx": -0.00, "ty": -0.00, "tz": -0.00, "mirror": 0,"scale": 1.0000, "type":"eman"})
    t2 = Transform({"az": 146.608, "alt": 143.005, "phi": 146.656, "tx": -0.00, "ty": -0.00, "tz": -0.00, "mirror": 0,"scale": 1.0000, "type":"eman"})
    def test_wrong_number_params_returns_TypeError_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.max_3D_pixel_error()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.max_3D_pixel_error()
        self.assertEqual(str(cm_new.exception), "max_3D_pixel_error() takes at least 2 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_pickle_file(self):
        return_new = fu.max_3D_pixel_error(t1=self.t1, t2=self.t2, r=self.r)
        return_old = oldfu.max_3D_pixel_error(t1=self.t1, t2=self.t2, r=self.r)
        self.assertEqual(return_new, return_old)
        self.assertEqual(return_new, 57.86697769165039)

    def test_r_is0(self):
        return_new = fu.max_3D_pixel_error(t1=self.t1, t2=self.t2, r=0)
        return_old = oldfu.max_3D_pixel_error(t1=self.t1, t2=self.t2, r=0)
        self.assertEqual(return_new, return_old)
        self.assertEqual(return_new, 0.0)

    def test_t1_emptyList_returns_IndexError(self):
        with self.assertRaises(IndexError) as cm_new:
            fu.max_3D_pixel_error(t1=[], t2=self.t2, r=self.r)
        with self.assertRaises(IndexError) as cm_old:
            oldfu.max_3D_pixel_error(t1=[], t2=self.t2, r=self.r)
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_t2_emptyList_returns_IndexError(self):
        with self.assertRaises(IndexError) as cm_new:
            fu.max_3D_pixel_error(t1=self.t1, t2=[], r=self.r)
        with self.assertRaises(IndexError) as cm_old:
            oldfu.max_3D_pixel_error(t1=self.t1, t2=[], r=self.r)
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_ali_params_both_emptyList_returns_IndexError(self):
        with self.assertRaises(IndexError) as cm_new:
            fu.max_3D_pixel_error(t1=[], t2=[], r=self.r)
        with self.assertRaises(IndexError) as cm_old:
            oldfu.max_3D_pixel_error(t1=[], t2=[], r=self.r)
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))



class Test_angle_ave(unittest.TestCase):
    angle1 = [345.74270725250244, 345.59336185455322, 345.83632707595825, 345.99127292633057, 346.02255821228027, 346.00519895553589, 347.74018049240112, 345.92735052108765, 348.68605613708496, 346.62910223007202, 359.24490451812744, 345.78871250152588, 346.73074722290039, 358.99380683898926, 345.88786840438843]

    def test_wrong_number_params_returns_TypeError_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.angle_ave()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.angle_ave()
        self.assertEqual(str(cm_new.exception), "angle_ave() takes exactly 1 argument (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_pickle_value(self):
        return_new = fu.angle_ave(angle1=self.angle1)
        return_old = oldfu.angle_ave(angle1=self.angle1)
        self.assertTrue(numpy_array_equal(return_new, return_old))
        self.assertTrue(numpy_array_equal(return_new, (348.04583441529024, 4.4154015011463477)))

    def test_emptyList_returns_ZeroDivisionError(self):
        with self.assertRaises(ZeroDivisionError) as cm_new:
            fu.angle_ave(angle1=[])
        with self.assertRaises(ZeroDivisionError) as cm_old:
            oldfu.angle_ave(angle1=[])
        self.assertEqual(str(cm_new.exception), "float division by zero")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))



class Test_angle_diff(unittest.TestCase):
    angle1 = [-358.1422162327943, 335.90588778212316, -81.05934436345304, 88.94114972876307, 230.3400198158369, 179.6066808748772, -179.05318018749753, -267.2248149400371, 168.26126969747668, 177.5181948630324, 90.44012904554053]
    angle2 = [-270.32559871309775, 359.28047486975726, -271.30870341528345, 65.26503368866545, 359.68416319428593, 182.75795469476327, -86.79404676631873, -357.78252122233323, 220.67100196439736, 347.68044245436516, 229.87355240518391]

    def test_wrong_number_params_returns_TypeError_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.angle_diff()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.angle_diff()
        self.assertEqual(str(cm_new.exception), "angle_diff() takes exactly 2 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_pickle_value(self):
        return_new = fu.angle_diff(angle1 = self.angle1,angle2= self.angle2)
        return_old = oldfu.angle_diff(angle1 = self.angle1,angle2= self.angle2)
        self.assertEqual(return_new, 88.918153887518429)
        self.assertEqual(return_new, return_old)

    def test_angle1_bigger_than_angle2(self):
        with self.assertRaises(IndexError) as cm_new:
            fu.angle_diff(angle1 = self.angle1+[2,2],angle2= self.angle2)
        with self.assertRaises(IndexError) as cm_old:
            oldfu.angle_diff(angle1 = self.angle1+[2,2],angle2= self.angle2)
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_angle2_bigger_than_angle1_printout_errorMsg(self):
        return_new = fu.angle_diff(angle1 = self.angle1,angle2= self.angle2+[2,2])
        return_old = oldfu.angle_diff(angle1 = self.angle1,angle2= self.angle2+[2,2])
        self.assertEqual(return_new, 88.918153887518429)
        self.assertEqual(return_new, return_old)

    def test_angle1_emptyList(self):
        return_new = fu.angle_diff(angle1 = [],angle2= self.angle2)
        return_old = oldfu.angle_diff(angle1 = [],angle2= self.angle2)
        self.assertEqual(return_new, 0.0)
        self.assertEqual(return_new, return_old)

    def test_angle2_emptyList_returns_IndexError(self):
        with self.assertRaises(IndexError) as cm_new:
            fu.angle_diff(angle1 = self.angle1,angle2= [])
        with self.assertRaises(IndexError) as cm_old:
            oldfu.angle_diff(angle1 = self.angle1,angle2= [])
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_angle_both_emptyList(self):
        return_new = fu.angle_diff(angle1 = [],angle2= [])
        return_old = oldfu.angle_diff(angle1 = [],angle2= [])
        self.assertEqual(return_new, 0.0)
        self.assertEqual(return_new, return_old)




class Test_angle_diff_sym(unittest.TestCase):
    angle1 = [[5.16760248821511, 102.82730387712962, 56.24999974223351, -0.0, -0.0], [5.394529010967219, 108.65382266507675, 57.65624899025124, -0.0, -0.0], [2.3675495391350125, 114.39230625756204, 54.84375196092924, -0.0, -0.0], [11.868847924292282, 120.4678965884309, 57.65624964621793, -0.0, -0.0], [6.624446177521037, 126.1456590935264, 54.84375044653484, -0.0, -0.0], [9.85254975964493, 139.8070854957106, 56.24999964448676, -0.0, -0.0], [10.284527175032494, 149.50936934645668, 57.65624974825016, -0.0, -0.0], [58.633634895542144, 152.74703226226322, 39.374999323555016, -0.0, -0.0], [4.697009771612514, 143.54534575495455, 50.62500035351184, -0.0, -0.0], [8.739223964927803, 153.34156227954617, 53.528990504798344, -0.0, -0.0], [13.215579639909052, 153.46194700482943, 56.24999989276381, -0.0, -0.0], [35.125741719584425, 155.98690869749788, 74.53125022294324, -0.0, -0.0], [24.241351890014528, 157.62937775317255, 46.40625150340372, -0.0, -0.0], [18.161124683110756, 81.83346531088333, 333.5816655275555, -0.0, -0.0], [14.119870407166076, 81.72454607009684, 332.9447033765355, -0.0, -0.0], [9.106476144519206, 79.64886950327454, 336.093749983792, -0.0, -0.0], [9.175438738145829, 89.35318537114819, 330.4687496386661, -0.0, -0.0], [7.175311305543758, 89.29928152249676, 333.2812488125304, -0.0, -0.0], [17.18271530409055, 91.50937668932787, 333.281249739689, -0.0, -0.0], [13.148079549129832, 87.51982195330109, 333.47210714792664, -0.0, -0.0], [9.00367345558179, 85.46766826184691, 331.87499935742085, -0.0, -0.0], [11.180777220614587, 91.34762616568598, 336.0937499007258, -0.0, -0.0], [5.179245407263252, 91.18588626611778, 331.8749991482967, -0.0, -0.0], [11.450337146268097, 95.23559046496233, 333.2812498303779, -0.0, -0.0], [7.433751874506456, 95.12734785957713, 336.0937502957362, -0.0, -0.0], [3.8208447469418445, 96.97028993475566, 334.68750081043464, -0.0, -0.0], [5.140075844399703, 87.30400152156635, 336.09374960687114, -0.0, -0.0], [4.643656184704199, 83.40966858984036, 333.2812496372485, -0.0, -0.0], [71.40261964863882, 94.91091772300047, 334.6875001281891, -0.0, -0.0], [10.603349995951802, 99.09333923199034, 334.68749995393057, -0.0, -0.0], [14.655192397786806, 99.20252805650983, 335.94796750543776, -0.0, -0.0], [13.538593406105079, 104.99304573002755, 336.0937493549564, -0.0, -0.0], [8.577894245539255, 99.03875710583709, 336.09374924950487, -0.0, -0.0], [10.182790116760444, 106.84275119823799, 334.6875002146449, -0.0, -0.0], [7.851333113661184, 97.07890586721807, 333.2812509877207, -0.0, -0.0], [9.866933990624162, 97.1332228238438, 331.87500111632653, -0.0, -0.0], [7.666315462821387, 100.95443140310614, 336.0937505837185, -0.0, -0.0], [11.244603045120286, 93.28971168469701, 336.0796199084946, -0.0, -0.0], [11.882774845480753, 97.187546373306, 336.12268113025533, -0.0, -0.0], [2.503367034950003, 98.8750603072273, 333.28125131573734, -0.0, -0.0], [9.398161895807235, 104.8814741425398, 334.6875000369289, -0.0, -0.0], [11.882774845480753, 97.187546373306, 334.6875011177296, -0.0, -0.0], [9.703812619306333, 101.00933703189095, 333.28124946500117, -0.0, -0.0], [4.527904954402473, 98.92961813936152, 336.09375126858004, -0.0, -0.0], [7.851333113661184, 97.07890586721807, 334.68750098436266, -0.0, -0.0], [7.23810893112136, 93.18173962163453, 333.28124941496867, -0.0, -0.0], [6.552746655416698, 98.98418321741717, 336.0937498838442, -0.0, -0.0], [5.235021200492042, 93.1277575502507, 334.6875004175527, -0.0, -0.0], [9.27134046402358, 102.93788788413879, 334.68750020714737, -0.0, -0.0], [5.835970636488611, 97.0245942754992, 336.0937490511087, -0.0, -0.0], [15.182021374401842, 91.45545851067713, 334.6875003445952, -0.0, -0.0], [9.241302521565729, 93.23572413398958, 336.0937499729066, -0.0, -0.0], [13.181375953848303, 91.40154176200272, 333.2812509830407, -0.0, -0.0], [12.629115548522591, 99.14792940839799, 334.6874997591077, -0.0, -0.0], [13.24801457530063, 93.3437024073173, 337.500000172, -0.0, -0.0], [15.467620867444467, 95.3438516377721, 336.09375026702287, -0.0, -0.0], [19.897058653402155, 101.28401964529849, 333.2812490442801, -0.0, -0.0], [15.915184504101504, 97.29621414748426, 337.49999906519645, -0.0, -0.0], [9.441957601783272, 95.18146658061052, 336.0937496887225, -0.0, -0.0], [13.779953550943318, 101.11917895309425, 334.68749940411817, -0.0, -0.0], [17.857634121295803, 101.22906160063684, 335.4797666079291, -0.0, -0.0], [17.484307516595152, 103.15920267814805, 337.5000003338388, -0.0, -0.0], [11.323892282464229, 102.99319712700886, 333.28125023694145, -0.0, -0.0], [15.818601298408666, 101.17411462087938, 337.5000007170344, -0.0, -0.0], [16.681583692662457, 99.25713457197845, 334.687499708489, -0.0, -0.0], [19.94857671233467, 97.40490735450291, 337.079010675876, -0.0, -0.0], [15.60962045500591, 105.04885124512585, 333.28124908958665, -0.0, -0.0], [21.593572086778565, 103.26993428093563, 336.0937488884011, -0.0, -0.0], [30.100090288944543, 101.55896424127138, 336.09375095478686, -0.0, -0.0], [28.018370162576304, 97.62237454203566, 334.687499771866, -0.0, -0.0], [25.513988462247426, 95.61458957120372, 334.6874998082104, -0.0, -0.0], [30.87527364673717, 99.63962468908254, 338.9062495754473, -0.0, -0.0], [31.874990040884242, 103.54698728725273, 337.50000039797715, -0.0, -0.0], [35.56504650091372, 95.88545259294263, 336.0937509638818, -0.0, -0.0], [39.303126016451245, 94.04586705852455, 337.9848935904527, -0.0, -0.0], [27.27510248619842, 93.72172158401631, 336.0937494424427, -0.0, -0.0], [18.708290336629204, 99.31175080426632, 334.6874998527977, -0.0, -0.0], [22.762658459045085, 99.42100732732668, 337.50000029396926, -0.0, -0.0], [13.89885770762399, 97.24187694122969, 338.90624954762507, -0.0, -0.0], [19.258920853617752, 93.50569178776301, 337.84103489688675, -0.0, -0.0], [11.146204697214003, 87.46587032396305, 337.499999237331, -0.0, -0.0], [15.467620867444467, 95.3438516377721, 334.68749980197634, -0.0, -0.0], [17.476531499835744, 95.39798951684963, 338.906250910386, -0.0, -0.0], [17.255170878212667, 93.45169230106706, 338.90624945820224, -0.0, -0.0], [17.931756284985397, 97.35055668796024, 337.4999999042536, -0.0, -0.0], [9.144246594507479, 87.41191644949221, 337.50000060410133, -0.0, -0.0], [16.455452947996136, 107.01177363163896, 340.31249921588466, -0.0, -0.0], [11.741692084016222, 101.06425227230656, 337.5000002600115, -0.0, -0.0], [13.458891659961509, 95.28971892217496, 337.49999922757036, -0.0, -0.0], [14.363937850820761, 106.95541497524451, 337.5000009635654, -0.0, -0.0], [17.68118909392578, 105.10467265590191, 339.62268155714503, -0.0, -0.0], [13.376903307961527, 103.04852041074584, 336.0937504239598, -0.0, -0.0], [7.219244624233198, 102.88258904815184, 334.68749892391844, -0.0, -0.0], [18.07484495396018, 108.99549899529111, 338.90624932076173, -0.0, -0.0], [7.506126483824204, 108.71072136056677, 336.0937507242717, -0.0, -0.0]]
    angle2 = [[21.297641186004853, 79.97745095305105, 57.65625182077656, -0.0, -0.0], [26.30217997306748, 76.23113845875804, 57.65624932460355, -0.0, -0.0], [29.358511191890443, 66.60987872663543, 55.30169707305913, -0.0, -0.0], [15.149302483452558, 58.46307341830345, 54.84375106444094, -0.0, -0.0], [18.528332563589473, 50.78936611013466, 57.656250796099016, -0.0, -0.0], [25.412645189110563, 41.26681410524694, 53.43749870198468, -0.0, -0.0], [26.118494919067558, 29.629787901538126, 46.40624897147546, -0.0, -0.0], [35.935343132914326, 27.950969903405397, 45.0, -0.0, -0.0], [30.64670316285857, 33.638480224545106, 47.812499223576594, -0.0, -0.0], [14.43598132692972, 27.370441262118504, 59.062499397206295, -0.0, -0.0], [19.254389939080298, 25.556052556681134, 54.8437502307637, -0.0, -0.0], [67.8574933356951, 24.924028412874847, 75.93749977799234, -0.0, -0.0], [13.345214600946619, 21.504981391057456, 42.18749921675186, -0.0, -0.0], [17.931756284985397, 97.35055668796024, 333.4150698401827, -0.0, -0.0], [11.180777220614587, 91.34762616568598, 331.87499958876185, -0.0, -0.0], [15.915184504101504, 97.29621414748426, 331.8749997188293, -0.0, -0.0], [21.175837019175745, 89.67659784741858, 333.2812502203766, -0.0, -0.0], [23.185108296154738, 91.67113929444578, 333.3251947339376, -0.0, -0.0], [13.17563428590529, 89.46099128532461, 331.87500003850175, -0.0, -0.0], [23.156279155997538, 87.78954872942698, 333.6575017496408, -0.0, -0.0], [21.154790888773945, 87.73560729140638, 331.8749996259746, -0.0, -0.0], [19.175806461820684, 89.62269673951845, 334.6875006879113, -0.0, -0.0], [25.175874555580577, 89.78439919513313, 331.87499998165066, -0.0, -0.0], [19.258920853617752, 93.50569178776301, 333.28125066659646, -0.0, -0.0], [27.053643998645086, 85.95413292369146, 334.68749900528707, -0.0, -0.0], [24.76493831041941, 83.95196603718402, 331.8750011628807, -0.0, -0.0], [22.753743514460965, 83.89776163095178, 334.6875003950914, -0.0, -0.0], [18.73074516637881, 83.78933637646092, 333.4909359753444, -0.0, -0.0], [26.240354037888906, 82.05121594323842, 334.6874995231353, -0.0, -0.0], [24.038603453923244, 78.11074374424224, 336.0937486480204, -0.0, -0.0], [20.123189745364826, 76.06459447131637, 334.6874994866322, -0.0, -0.0], [17.90582058621864, 77.94544781700415, 336.09375079384176, -0.0, -0.0], [21.76846376912684, 74.16828765962819, 331.87499946276114, -0.0, -0.0], [23.32829413725878, 80.03218081927426, 331.8750004859039, -0.0, -0.0], [20.742346015222637, 83.84355177937768, 333.28124901455203, -0.0, -0.0], [17.235307961838686, 79.86796188462372, 333.28125051705865, -0.0, -0.0], [20.18134074871284, 81.88791400100129, 331.87500053202626, -0.0, -0.0], [21.9947582534207, 78.05565667174946, 334.4925843920765, -0.0, -0.0], [16.140637032013714, 81.77900965371339, 337.2022401937973, -0.0, -0.0], [16.001394119760576, 75.95349894118092, 337.5000005736099, -0.0, -0.0], [13.939749465729264, 75.8979300077477, 336.09375024789136, -0.0, -0.0], [19.153227784186655, 87.68166417916409, 334.6874998254447, -0.0, -0.0], [21.038248582193006, 85.79201211540865, 333.28124917530556, -0.0, -0.0], [19.26664703037237, 79.92271029470429, 333.28125031416545, -0.0, -0.0], [20.18134074871284, 81.88791400100129, 335.30651879879395, -0.0, -0.0], [24.220952114037843, 81.99678891824239, 334.687499658415, -0.0, -0.0], [17.15159050224196, 87.62771895518732, 334.68750071708075, -0.0, -0.0], [22.201282396357783, 81.94235469274858, 334.68750097386874, -0.0, -0.0], [16.718937873778543, 83.73511534520249, 333.2812496237367, -0.0, -0.0], [23.043514819351657, 85.84605595290141, 334.68749989361834, -0.0, -0.0], [17.0272955263693, 85.68391315629866, 336.0937497603034, -0.0, -0.0], [15.14987443758487, 87.57377167219413, 336.093749390464, -0.0, -0.0], [17.18271530409055, 91.50937668932787, 333.281249739689, -0.0, -0.0], [15.203624421267165, 79.81320273908536, 336.093750362957, -0.0, -0.0], [17.17576286685079, 89.56879534280377, 337.499999943684, -0.0, -0.0], [14.706921707167226, 83.6808882266122, 336.0937497425133, -0.0, -0.0], [12.098826473487563, 81.67007466816304, 336.09374961329587, -0.0, -0.0], [15.02160701045571, 85.62985784549772, 336.0937503132991, -0.0, -0.0], [19.032841383475727, 85.73796466549993, 337.5000001820951, -0.0, -0.0], [15.860724873649104, 77.89032718690052, 334.6874994880993, -0.0, -0.0], [10.077502471194038, 81.61559643428924, 335.23565613413626, -0.0, -0.0], [15.530157726646536, 74.00014022346902, 335.5092780271573, -0.0, -0.0], [11.13920956205078, 79.70365726945819, 333.2812502412992, -0.0, -0.0], [12.694696318733676, 83.6266557893642, 336.58859302604475, -0.0, -0.0], [13.171592561362061, 79.75843412355391, 334.6874995705287, -0.0, -0.0], [9.106476144519206, 79.64886950327454, 337.49999888322805, -0.0, -0.0], [9.814953001434901, 75.78675294655795, 336.09375024137444, -0.0, -0.0], [3.81698347037117, 69.80298583131767, 336.09374987682844, -0.0, -0.0], [7.073387604522281, 79.59407341671378, 337.5000006466972, -0.0, -0.0], [69.43564848971883, 77.39371866612434, 337.0431507802318, -0.0, -0.0], [68.90256224029883, 79.31994488855202, 337.4538583929192, -0.0, -0.0], [70.93745260278942, 79.37479008466487, 338.90625004784096, -0.0, -0.0], [67.38624238486705, 77.33848081740646, 338.9062502696945, -0.0, -0.0], [60.93392144872331, 84.92676649799077, 336.0937500798579, -0.0, -0.0], [65.17721109912068, 90.86243373540891, 337.4999998828992, -0.0, -0.0], [3.1378599278328494, 87.25004057768219, 338.9062499383664, -0.0, -0.0], [4.643656184704199, 83.40966858984036, 337.49999944801027, -0.0, -0.0], [13.815208229199186, 77.8351942829876, 336.36907836422296, -0.0, -0.0], [13.015775850438331, 85.57579881267367, 338.90625022818614, -0.0, -0.0], [13.148079549129832, 87.51982195330109, 336.11328022053294, -0.0, -0.0], [17.255170878212667, 93.45169230106706, 338.90624945820224, -0.0, -0.0], [15.175706134571186, 89.5148935342978, 336.09374988299925, -0.0, -0.0], [10.68225901059121, 83.57241764071574, 338.90625002518703, -0.0, -0.0], [11.877602499687342, 75.84234920940848, 338.9062493960632, -0.0, -0.0], [16.140637032013714, 81.77900965371339, 337.4999999090718, -0.0, -0.0], [18.06254142500468, 76.009052783409, 338.6346440529167, -0.0, -0.0], [18.161124683110756, 81.83346531088333, 338.9062489929637, -0.0, -0.0], [13.449560754983963, 73.94406117392357, 337.61929356282656, -0.0, -0.0], [15.860724873649104, 77.89032718690052, 337.5000003023005, -0.0, -0.0], [19.950498765577194, 78.00055803807608, 338.9062492669124, -0.0, -0.0], [9.722901027948168, 77.72489512822938, 340.3125001332472, -0.0, -0.0], [17.61017292978873, 74.05620533778338, 337.5000008929808, -0.0, -0.0], [16.446428988822134, 72.0841448424426, 334.6874995786454, -0.0, -0.0], [14.460455758454842, 70.0898797223358, 337.5000007581206, -0.0, -0.0], [13.404246631068986, 64.23894693338787, 336.0937509878121, -0.0, -0.0]]
    def test_wrong_number_params_returns_TypeError_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.angle_diff_sym()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.angle_diff_sym()
        self.assertEqual(str(cm_new.exception), "angle_diff_sym() takes at least 2 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_pickle_value(self):
        return_new = fu.angle_diff_sym(angle1=self.angle1, angle2=self.angle2, simi=1)
        return_old = oldfu.angle_diff_sym(angle1=self.angle1, angle2=self.angle2, simi=1)
        self.assertEqual(return_old,return_new)
        self.assertEqual(12.485941961880933, return_new)

    def test_sim_is0_returns_ZeroDivisionError(self):
        with self.assertRaises(ZeroDivisionError) as cm_new:
            fu.angle_diff_sym(angle1=self.angle1, angle2=self.angle2, simi=0)
        with self.assertRaises(ZeroDivisionError) as cm_old:
            oldfu.angle_diff_sym(angle1=self.angle1, angle2=self.angle2, simi=0)
        self.assertEqual(str(cm_new.exception), "float division by zero")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_angle1_bigger_than_angle2(self):
        with self.assertRaises(IndexError) as cm_new:
            fu.angle_diff_sym(angle1=self.angle1+[2], angle2=self.angle2, simi=1)
        with self.assertRaises(IndexError) as cm_old:
            oldfu.angle_diff_sym(angle1=self.angle1+[2], angle2=self.angle2, simi=1)
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_angle2_bigger_than_angle1_printout_errorMsg(self):
        return_new = fu.angle_diff_sym(angle1=self.angle1, angle2=self.angle2+[2], simi=1)
        return_old = oldfu.angle_diff_sym(angle1=self.angle1, angle2=self.angle2+[2], simi=1)
        self.assertEqual(return_old,return_new)
        self.assertEqual(12.485941961880933, return_new)

    def test_angle1_emptylist_printout_errorMsg(self):
        return_new = fu.angle_diff_sym(angle1=[], angle2=self.angle2, simi=1)
        return_old = oldfu.angle_diff_sym(angle1=[], angle2=self.angle2, simi=1)
        self.assertEqual(return_old,return_new)
        self.assertEqual(0.0, return_new)

    def test_angle_both_emptylist(self):
        return_new = fu.angle_diff_sym(angle1=[], angle2=[], simi=1)
        return_old = oldfu.angle_diff_sym(angle1=[], angle2=[], simi=1)
        self.assertEqual(return_old,return_new)
        self.assertEqual(0.0, return_new)

    def test_angle_both_emptyList(self):
        return_new = fu.angle_diff_sym(angle1=[], angle2=[], simi=1)
        return_old = oldfu.angle_diff_sym(angle1=[], angle2=[], simi=1)
        self.assertEqual(return_old,return_new)
        self.assertEqual(0.0, return_new)




class Test_align_diff_params(unittest.TestCase):
    ali_params1 = [358.1422162327943, 1.5746676921844482, 2.359189510345459, 1, 335.90588778212316, 0.22259289026260376, 1.6781892776489258, 0, 81.05934436345304, 2.0682716369628906, 2.3159945011138916, 1, 88.94114972876307, 2.3903369903564453, 2.760911226272583, 0, 230.3400198158369, -1.2752182483673096, -1.0786179304122925, 0, 179.6066808748772, 0.3892550766468048, -0.7324743270874023, 0, 179.05318018749753, -0.694373607635498, -0.10197801142930984, 1, 267.2248149400371, 0.5244513750076294, 0.7750062346458435, 1, 168.26126969747668, 1.208488941192627, 1.3360750675201416, 0, 177.5181948630324, 1.0620733499526978, -1.5996335744857788, 0, 90.44012904554053, 0.8572239875793457, -1.572252869606018, 0]
    ali_params2 = [358.142216327943, 1.5746676921844482, 2.359189510345459, 1, 335.90588778212316, 0.22259289026260376, 1.6781892776489258, 0, 81.05934436345304, 2.0682716369628906, 2.3159945011138916, 1, 88.94114972876307, 2.3903369903564453, 2.760911226272583, 0, 230.3400198158369, -1.2752182483673096, -1.0786179304122925, 0, 179.6066808748772, 0.3892550766468048, -0.7324743270874023, 0, 179.05318018749753, -0.694373607635498, -0.10197801142930984, 1, 267.2248149400371, 0.5244513750076294, 0.7750062346458435, 1, 168.26126969747668, 1.208488941192627, 1.3360750675201416, 0, 177.5181948630324, 1.0620733499526978, -1.5996335744857788, 0, 90.44012904554053, 0.8572239875793457, -1.572252869606018, 0]

    def test_wrong_number_params_returns_TypeError_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.align_diff_params()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.align_diff_params()
        self.assertEqual(str(cm_new.exception), "align_diff_params() takes exactly 2 arguments (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_pickle_file(self):
        return_new = fu.align_diff_params(ali_params1= self.ali_params1, ali_params2= self.ali_params2)
        return_old = oldfu.align_diff_params(ali_params1= self.ali_params1, ali_params2= self.ali_params2)
        self.assertTrue(numpy_array_equal(return_new, return_old))
        self.assertTrue(numpy_array_equal(return_new, (359.99999999135014, 0.0, 0.0, 0)))

    def test_ali_params1_emptyList_printout_ErrorMsg(self):
        return_new = fu.align_diff_params(ali_params1= [], ali_params2= self.ali_params2)
        return_old = oldfu.align_diff_params(ali_params1= [], ali_params2= self.ali_params2)
        self.assertTrue(numpy_array_equal(return_new, return_old))
        self.assertTrue(numpy_array_equal(return_new, (0.0, 0.0, 0.0, 0)))

    def test_ali_params2_emptyList_printout_ErrorMsg(self):
        return_new = fu.align_diff_params(ali_params1= self.ali_params1, ali_params2= [])
        return_old = oldfu.align_diff_params(ali_params1= self.ali_params1, ali_params2= [])
        self.assertTrue(numpy_array_equal(return_new, return_old))
        self.assertTrue(numpy_array_equal(return_new, (0.0, 0.0, 0.0, 0)))

    def test_ali_params_both_emptyList_returns_ZeroDivisionError(self):
        with self.assertRaises(ZeroDivisionError) as cm_new:
            fu.align_diff_params(ali_params1= [], ali_params2= [])
        with self.assertRaises(ZeroDivisionError) as cm_old:
            oldfu.align_diff_params(ali_params1= [], ali_params2= [])
        self.assertEqual(str(cm_new.exception), "float division by zero")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))



class Test_multi_align_stability(unittest.TestCase):
    ali_params =  [[181.31269459154402, -2.5839810371398926, -2.0082650184631348, 1, 90.25404811898599, -1.171455979347229, 1.9825596809387207, 1, 158.28384987667707, 0.10442383587360382, 0.43271294236183167, 1, 177.09224922620902, 1.1460936069488525, 2.6825344562530518, 0, 224.5203330012341, 0.255629301071167, -0.563730001449585, 0, 91.16107461457933, -4.165151119232178, 1.8753869533538818, 0, 91.30484106770238, 1.4268516302108765, -2.4773027896881104, 0, 123.85811587900938, -5.611675262451172, -1.9485784769058228, 0, 172.59703867632064, -0.15322698652744293, -0.9023720026016235, 1, 211.98109486449883, -1.4229055643081665, -0.0031406283378601074, 0, 336.0930857469997, -1.6030802726745605, 2.498553991317749, 0, 181.61649229022615, 1.449699878692627, 1.3747519254684448, 1, 239.52141516359143, -2.080139636993408, 2.982375383377075, 1, 57.36939212806606, 2.8285839557647705, 1.273552417755127, 0, 270.6765818845557, 0.755837619304657, -0.12974362075328827, 0, 42.25069201252337, 2.3062591552734375, 1.6736345291137695, 1, 270.10040045188396, -0.8220255970954895, 2.439899444580078, 0, 31.10853891661594, -0.04069456458091736, -0.8928937315940857, 0, 89.40835835581117, 1.7607072591781616, 0.42965540289878845, 0, 179.72021342271486, 0.4846632182598114, -0.7986509799957275, 1, 178.18965916962765, 4.561519622802734, -1.4123663902282715, 0, 272.87569293897354, -0.7022020220756531, 0.4636093080043793, 1, 269.10935639342426, 0.7251760363578796, 3.3150687217712402, 0, 268.5329174739993, 1.863369107246399, 3.4191946983337402, 1, 338.5374143678607, 0.4387291669845581, 4.636652946472168, 1, 182.89386744019887, 3.6257309913635254, 0.6001610159873962, 1, 150.28444979455818, -4.206900596618652, -0.6949625015258789, 0, 358.1422162327943, 1.5746676921844482, 2.359189510345459, 1, 233.54686846204189, 0.5799525380134583, -0.08851173520088196, 0, 335.90588778212316, 0.22259289026260376, 1.6781892776489258, 0, 270.37806272671065, 0.8993915319442749, -0.4666183888912201, 0, 178.45095994531277, -0.6545171737670898, -7.941856384277344, 1, 81.05934436345304, 2.0682716369628906, 2.3159945011138916, 1, 88.94114972876307, 2.3903369903564453, 2.760911226272583, 0, 71.17695387089937, 4.591431140899658, 0.20732438564300537, 1, 301.18870540563, -1.9143840074539185, 0.5434432029724121, 0, 178.3262479018072, -4.607028961181641, -1.3345184326171875, 0, 180.6416916974331, 0.025997497141361237, 0.5040155649185181, 1, 230.3400198158369, -1.2752182483673096, -1.0786179304122925, 0, 247.94419511173476, 0.7020229697227478, -2.228996992111206, 1, 90.36840677809415, -1.8641890287399292, -3.533489465713501, 0, 338.96053751136895, -3.208756923675537, 0.3189549446105957, 0, 331.883712061379, -2.2741057872772217, -0.028077542781829834, 0, 179.6066808748772, 0.3892550766468048, -0.7324743270874023, 0, 1.023211479940733, 2.0373544692993164, 0.46514835953712463, 1, 180.83373783035353, -1.1226463317871094, 0.45663633942604065, 0, 263.9159129584964, -0.9649151563644409, -3.129274845123291, 1, 231.60703053363372, 1.204798936843872, 1.0629167556762695, 0, 335.35343857446617, -0.2024703025817871, 1.640763282775879, 0, 161.5998308962837, 2.3083269596099854, 1.563036322593689, 0, 274.17946114304044, 0.20257164537906647, 1.9253363609313965, 0, 333.9763718846628, 3.883352756500244, -2.4557318687438965, 0, 1.0879650258769762, -2.565249443054199, -3.21978759765625, 0, 2.7111044348811846, -0.17390495538711548, 2.368929862976074, 1, 218.0498613780321, 2.747969150543213, 0.44726407527923584, 0, 179.05318018749753, -0.694373607635498, -0.10197801142930984, 1, 2.2084069326124345, -1.9855765104293823, 0.6704190373420715, 1, 179.56593277553742, -4.927476406097412, -2.0986437797546387, 0, 180.13681411653684, -0.68706876039505, 3.342438220977783, 1, 192.62803291748114, 0.00803869217634201, -0.5794680714607239, 1, 54.85455243874288, -1.677180290222168, 1.4227328300476074, 0, 211.97135165250157, 0.7141838073730469, -1.2560396194458008, 1, 260.37483947003966, -5.181369781494141, -4.346599578857422, 0, 95.48371216913883, 1.2659498453140259, -1.6158771514892578, 0, 268.06889535064465, -1.5196741819381714, 1.1935900449752808, 1, 267.2248149400371, 0.5244513750076294, 0.7750062346458435, 1, 91.8551016618861, 0.5114099979400635, -2.622548818588257, 1, 359.9192382674421, 4.417739391326904, -3.760923147201538, 1, 167.97475546197109, -1.3366694450378418, 3.706094741821289, 1, 181.04376560956464, 4.5168962478637695, 2.8973283767700195, 1, 180.0648021683185, -0.09634430706501007, -1.843454122543335, 1, 89.7879016373518, -0.18748335540294647, 0.3684837222099304, 1, 266.93874127689224, -0.7439479231834412, -3.3092563152313232, 1, 325.9200222274459, -5.604404449462891, 2.824554204940796, 1, 241.9676970804552, -1.5311527252197266, 0.5358741879463196, 1, 3.3323932269925707, -0.24386408925056458, -0.08598934859037399, 1, 189.8085811586751, 0.26547837257385254, -1.8841137886047363, 1, 156.34082431095638, -0.5988807678222656, -0.368464857339859, 0, 5.861479067579755, -1.952226161956787, 2.6622698307037354, 1, 183.1619381154372, -0.054596513509750366, -2.522308111190796, 1, 191.06640793262434, 0.2508540153503418, -0.8341284394264221, 0, 90.00627636894232, 0.28989264369010925, -4.166843414306641, 1, 168.26126969747668, 1.208488941192627, 1.3360750675201416, 0, 115.85709459617763, 1.3527047634124756, -2.044482707977295, 0, 90.95989819269074, 0.7148043513298035, -0.2923847734928131, 0, 177.5181948630324, 1.0620733499526978, -1.5996335744857788, 0, 104.50916313104166, -0.9761509299278259, -0.08903886377811432, 0, 270.37962914134505, 0.1329900026321411, 1.5415745973587036, 0, 87.46173443135923, -0.7523272037506104, 1.6808485984802246, 0, 30.04302416999718, -0.43151432275772095, -1.7011327743530273, 1, 173.1377185987646, -0.2973138988018036, 2.1424126625061035, 1, 213.3448374804534, -0.43501055240631104, -3.688663959503174, 1, 12.31113174798294, -4.407229423522949, 2.5972700119018555, 0, 359.99374148443155, -0.6971377730369568, 0.7342761754989624, 1, 265.93550920260486, -0.24208670854568481, -0.08130913227796555, 0, 358.89853480964916, -0.39681124687194824, 1.262823224067688, 0, 90.44012904554053, 0.8572239875793457, -1.572252869606018, 0, 300.1424519922064, -0.33729061484336853, -0.7680656313896179, 1, 177.1296823645375, 2.6842546463012695, -2.505058765411377, 1, 353.2451033862821, -2.9232842922210693, -5.195979118347168, 1], [247.071662782531, -1.3337912559509277, -2.152374744415283, 0, 181.06554509219995, 0.9991251230239868, 1.1618562936782837, 0, 270.01553535419134, 0.36221086978912354, -4.2572174072265625, 1, 88.39909136932091, -0.6564781665802002, 1.262850284576416, 1, 87.14155914318746, -1.8537553548812866, 1.4869365692138672, 1, 90.76525567671831, -3.1914520263671875, 0.8333014249801636, 0, 88.98941874352857, 1.5610458850860596, -4.417323589324951, 1, 305.8383890799799, 3.8291423320770264, 4.320676326751709, 0, 316.6905004224449, 0.28083959221839905, 0.09310165792703629, 1, 271.8210268884862, -3.747871160507202, 0.3245396912097931, 0, 179.72049236325367, 2.462404489517212, -1.9218522310256958, 0, 269.2794084723601, 0.432804137468338, -1.4050098657608032, 1, 343.196904881169, 1.475010633468628, 0.7379008531570435, 1, 358.7027120345546, -0.6392879486083984, 4.055423736572266, 0, 175.51812515336346, -0.9356815814971924, 0.6862978339195251, 0, 69.15436053483751, -0.7013748288154602, 1.2504613399505615, 1, 134.08616056109602, -3.9310696125030518, -2.371527671813965, 1, 88.27803262932673, -1.7117913961410522, -2.4790396690368652, 0, 89.93222594366095, 0.7633799314498901, 1.4123553037643433, 1, 270.80361127505114, -3.821402072906494, -1.4273020029067993, 1, 1.6281673297288235, -2.4409921169281006, 2.6261820793151855, 1, 262.94523231666903, 3.1980841159820557, 0.8268386125564575, 0, 123.60915397190736, -3.3081367015838623, -1.7677745819091797, 0, 88.3964360195968, -0.023905247449874878, 2.6299943923950195, 1, 358.2700181093825, 1.0996569395065308, 0.18788741528987885, 0, 267.5573015798942, -1.1060279607772827, -2.6403470039367676, 1, 268.21566586684673, -2.4858200550079346, -1.079756498336792, 0, 1.1175450715780042, 4.713932991027832, 3.21557354927063, 0, 76.37942609575155, -2.443915843963623, -0.16436150670051575, 0, 344.73157870479184, 2.406846284866333, 2.150850296020508, 1, 114.84708244058277, -1.1125738620758057, 2.124868869781494, 0, 270.84238529271886, -6.907739162445068, 0.9706274271011353, 0, 134.80240560588413, 0.9718772768974304, -1.012520432472229, 0, 306.10948676108484, -1.9570173025131226, -1.934586524963379, 1, 176.3296760343071, -0.0661172866821289, -3.4240264892578125, 0, 1.3395791053350092, 0.42474156618118286, -2.0913829803466797, 0, 91.45084980280613, -0.8926230072975159, -3.622603178024292, 0, 269.6529006804276, 0.5043893456459045, -0.01729593053460121, 1, 0.6297005074123022, -0.009120743721723557, 0.6702684164047241, 1, 265.7841897157143, 0.9091724753379822, -1.26608145236969, 0, 88.87517808265684, 2.188490152359009, -1.5027315616607666, 0, 302.6200393838466, -5.8395490646362305, -0.8695988655090332, 0, 216.65174292324988, 0.19253510236740112, -1.4482128620147705, 1, 1.5793017164748646, -3.3075497150421143, 2.827359914779663, 1, 90.4515123557062, -0.5302762985229492, -4.024670124053955, 0, 356.39736436704436, 0.1565800905227661, -0.43126562237739563, 0, 0.6615014887419176, -0.9944058656692505, 1.3027127981185913, 0, 133.07753086113533, -1.1824122667312622, -0.3795945346355438, 0, 310.50633657533274, -4.564042091369629, 3.2415976524353027, 1, 180.34091948500665, 0.6881660223007202, 0.7504172325134277, 0, 265.38505565729025, 4.889618396759033, 2.3359663486480713, 1, 194.05589114795802, -0.6630386114120483, 3.166377544403076, 0, 179.9042129514465, 2.4931700229644775, 6.272096157073975, 0, 91.38321523622426, 1.3645554780960083, 0.25289493799209595, 1, 279.39112201608793, -1.9096702337265015, 0.42907121777534485, 0, 271.70564882315676, -2.0688493251800537, 0.7578786015510559, 0, 119.6775376156264, 1.5107234716415405, 1.4524766206741333, 0, 357.4554848730085, 0.939584493637085, 0.10306869447231293, 0, 359.27012066373584, 0.7885023355484009, -7.331338405609131, 1, 356.12824190825285, -3.1654365062713623, 0.3507530093193054, 1, 88.19117195848206, -1.7450981140136719, 6.076779365539551, 0, 266.335995087882, 1.3911774158477783, -1.1844356060028076, 0, 88.53166226459989, 1.7213499546051025, 4.464653015136719, 0, 358.14815997789685, 1.4087002277374268, 2.46138072013855, 1, 0.8299068252524734, 0.2510281205177307, 0.475003719329834, 1, 358.7989711850176, -0.2185131013393402, -1.566283941268921, 0, 178.71661665032414, -2.6130130290985107, 0.3455228805541992, 1, 270.95705504466895, 5.806636810302734, 2.3157684803009033, 0, 179.60964204615198, 1.458612322807312, 0.9132119417190552, 0, 268.2878494440546, 0.08261682093143463, -3.4624459743499756, 1, 268.6997723391242, -3.8447105884552, 0.007019639015197754, 0, 180.22045611604148, -3.633929491043091, -0.7999054193496704, 0, 12.394103409710988, -1.2525818347930908, 0.19295716285705566, 1, 247.60050771065096, -6.130857944488525, -4.754083156585693, 1, 334.3032144370768, 1.0314335823059082, 0.6069421768188477, 0, 88.85487913527035, -2.0843801498413086, -0.8033629655838013, 0, 179.35090537948523, -0.37422609329223633, -3.815830945968628, 1, 1.3577482616423755, -0.2775101065635681, 1.1040524244308472, 0, 233.20163767822714, 0.7648578882217407, -3.44146728515625, 0, 179.94714975338385, 1.0897048711776733, -5.5204758644104, 1, 88.56860164311836, -0.26459044218063354, 1.40012526512146, 0, 202.72219384398366, -3.03318190574646, 0.9554083347320557, 0, 272.69643552813363, -1.14628005027771, -4.405988693237305, 1, 0.09803658369342699, -1.7499648332595825, 0.11192885041236877, 0, 295.0379067058912, -2.7685956954956055, 0.4989778399467468, 0, 91.4607524466645, 3.6427314281463623, -0.10108429193496704, 1, 91.01588132475773, 4.07078218460083, -0.40299543738365173, 0, 177.4094366842803, 1.450539231300354, 0.18852785229682922, 0, 357.0172333153138, -3.776323080062866, 1.1278741359710693, 0, 119.22890019197273, -0.8344163298606873, -0.0809997022151947, 0, 89.6638119139233, -1.156471610069275, -1.0460052490234375, 1, 319.58795203231205, -4.8294291496276855, 1.5625371932983398, 1, 182.82267097219102, 0.7306990623474121, -2.637016534805298, 0, 90.84104417820315, 1.7443972826004028, 0.6715249419212341, 0, 183.83918512813423, -0.7496392726898193, 2.8092567920684814, 1, 49.28124021939194, 2.130016803741455, 1.0053449869155884, 0, 134.592079949396, 1.6461796760559082, -2.4171810150146484, 1, 183.06911230259388, 0.8909564018249512, 1.0478039979934692, 0, 205.0047543064103, 1.6242637634277344, -2.5631320476531982, 0, 175.41306859999446, 3.437636137008667, 1.0944865942001343, 1], [105.29634933002345, 0.8873910903930664, -0.7997310161590576, 1, 272.65270469975064, 1.272582769393921, 1.9658623933792114, 0, 358.6524153101049, -2.224431276321411, -1.4160702228546143, 1, 24.12981149441228, -0.0061986446380615234, -1.7733607292175293, 1, 251.98529823984256, 0.2760835289955139, -1.5690689086914062, 0, 359.8511270917563, -1.774480938911438, -6.206929683685303, 1, 180.12297391777867, -2.449492931365967, -2.4776430130004883, 0, 88.82242441619591, -0.47678864002227783, -4.7554192543029785, 1, 354.4975235448929, 0.27895745635032654, -0.09859801828861237, 0, 3.348234094394551, -0.6907907724380493, 1.7997089624404907, 0, 272.7278494176578, -0.003068000078201294, -1.4545012712478638, 0, 177.2854221722577, -0.37185272574424744, -4.607682228088379, 1, 230.73660112951828, -1.2454447746276855, 1.0812126398086548, 0, 268.6968756122271, -4.032615661621094, -1.639441967010498, 1, 124.48487190598046, 0.2684389352798462, -0.03782404959201813, 1, 178.66318702859311, -1.4929673671722412, -3.825385332107544, 0, 49.75136315104843, 4.027347087860107, -5.026365756988525, 1, 294.8594221954073, 4.03443717956543, -0.29709696769714355, 1, 0.43925280704196806, -2.382516622543335, 3.783426523208618, 0, 185.38433102718898, -0.8736227750778198, -3.7355124950408936, 1, 267.80226222596684, -3.380021810531616, -4.647661209106445, 1, 0.1657409059351419, -1.5816720724105835, -3.270862340927124, 0, 211.33354526107647, -1.563599944114685, 1.8611624240875244, 0, 88.96835802577546, 2.002023220062256, 2.666111469268799, 0, 12.321440208029154, 0.5621973276138306, 2.0826048851013184, 1, 176.28715753020273, -0.32948869466781616, -1.2414969205856323, 0, 244.73558702326608, -3.130171298980713, 0.7322671413421631, 0, 93.01874055952247, 3.0574140548706055, -4.818017482757568, 0, 165.11773002763832, -0.6709222197532654, 0.24995937943458557, 0, 268.35466860617083, -3.4364829063415527, -0.21045765280723572, 1, 1.1107602788936086, -0.5362369418144226, -3.8927876949310303, 0, 268.4332895534161, -6.942432880401611, 0.6794074177742004, 0, 288.2039803126255, -0.3724100887775421, 0.15781790018081665, 0, 218.06389293478043, 1.6959348917007446, -0.6179931163787842, 1, 267.95634977366205, -3.3851089477539062, -0.8360757231712341, 0, 359.7573521952283, 1.4780824184417725, -1.0746318101882935, 0, 268.7979197776867, 1.7867072820663452, 0.6816751956939697, 1, 180.9238600460254, 0.012355613522231579, -0.49598827958106995, 1, 226.1777044971778, -0.4720751941204071, -0.47590744495391846, 1, 274.1290781472676, 0.5717906951904297, -3.379434823989868, 0, 87.34395862716482, 0.23001253604888916, -1.5363950729370117, 1, 94.46491258451158, 4.816094875335693, 4.519120216369629, 1, 225.41233206865326, 0.6615402698516846, -3.5869951248168945, 0, 1.0037803256500655, -4.335628986358643, 2.8115110397338867, 1, 269.4466066262449, -0.48239246010780334, -1.9760277271270752, 1, 181.152899193386, 1.8591870069503784, -0.597276508808136, 0, 1.6713189602035101, -1.00045907497406, 0.32046088576316833, 1, 86.86244849508203, -1.6520953178405762, 0.825973391532898, 1, 86.46168341262324, 5.410813331604004, 2.839141368865967, 0, 270.5460119386429, 0.7288897633552551, -2.6907572746276855, 0, 92.74146443535706, -3.0551018714904785, 0.21138572692871094, 1, 344.1328238554575, 5.171586036682129, -5.7143144607543945, 1, 177.0094335304731, -0.7705074548721313, 5.234849452972412, 1, 178.78899341224215, 0.3566676080226898, -3.351263999938965, 0, 141.48433457459416, 0.5067294239997864, -0.8160825967788696, 0, 359.85240423499727, 1.6957247257232666, 0.09485982358455658, 0, 358.5533881375062, -1.9737849235534668, -1.4568852186203003, 1, 356.60131933563247, 1.6990517377853394, 4.169312000274658, 1, 271.6303325187775, 7.35759162902832, 0.4859154224395752, 1, 226.51907632477483, 5.263437747955322, -3.4628939628601074, 1, 181.67972795736773, 5.142825126647949, 0.4018394947052002, 0, 269.71768141367636, 2.3287312984466553, -3.2594826221466064, 1, 177.13293301677945, 3.5566134452819824, -2.6606216430664062, 0, 299.96690504393905, 0.017043232917785645, 2.861538887023926, 1, 268.92280343643756, -1.4453186988830566, -1.783364176750183, 1, 1.927669022606575, -1.2021944522857666, 1.4799168109893799, 0, 358.8272953478933, 2.6328060626983643, -1.3503600358963013, 1, 92.68048166199868, -7.871467590332031, -2.0465545654296875, 0, 181.3553952675547, -2.4421634674072266, 3.962128162384033, 1, 271.77873376402147, 4.931306838989258, -1.6172147989273071, 1, 358.75669955895245, 1.0540004968643188, 1.8668705224990845, 0, 0.9484712367727945, 2.5942771434783936, -2.22936749458313, 0, 88.92034645144541, 1.8954397439956665, 1.300750970840454, 0, 16.945480797387006, 0.00109100341796875, 3.637495756149292, 0, 58.253340738642976, 2.4130353927612305, 0.09060752391815186, 0, 89.70801473668509, -0.09613664448261261, -0.7620460987091064, 0, 91.53129463576676, 4.745028018951416, -2.5451412200927734, 1, 88.85259389934205, 1.1108918190002441, -0.6742971539497375, 0, 268.3445406111423, -0.35594838857650757, -3.225529432296753, 1, 91.39120809841494, -0.40348950028419495, 3.095320701599121, 0, 240.67139582084212, 0.8889416456222534, -1.1136155128479004, 1, 339.0454884693399, 1.5620453357696533, 1.979921579360962, 1, 178.91265389439366, 3.4722094535827637, -0.8720005750656128, 1, 266.57563446847354, -0.06393247842788696, -0.7553210258483887, 1, 26.690203971279782, 0.12943387031555176, 1.8596086502075195, 0, 358.9065599392202, -0.023185912519693375, 1.643980860710144, 1, 268.4969758978855, -4.084559917449951, 0.2236994057893753, 1, 268.8265228151947, 1.2138270139694214, -4.433650016784668, 0, 1.0119770607796166, -1.706548810005188, 0.35304442048072815, 1, 120.18035415704036, -1.3383698463439941, -0.9315799474716187, 0, 2.7495326524153265, 1.8851451873779297, -3.2567687034606934, 1, 1.4600172578312607, -3.527421236038208, 5.41217565536499, 0, 269.111609440622, 0.38443732261657715, 1.1464262008666992, 0, 88.3226001115354, 1.7132048606872559, 0.7475265264511108, 0, 219.7000628520578, 0.5300384163856506, 0.5383257269859314, 0, 211.9680672317204, -0.356535404920578, -1.2747799158096313, 1, 0.3302795001459203, -1.4235786199569702, 1.8775123357772827, 1, 91.37212162455089, 2.021193504333496, 4.786551475524902, 1, 121.68462299513658, 2.535135269165039, -0.9120575785636902, 1, 358.3957386130166, -4.461540222167969, -1.9417345523834229, 0], [157.35198495205555, 0.3000476658344269, -2.11444354057312, 1, 358.6251210971613, -3.875901222229004, -4.2744646072387695, 0, 2.978662892165292, -0.3276173174381256, -1.3481621742248535, 0, 0.9825075981488491, -0.29130375385284424, -0.6160370111465454, 1, 156.21827885632732, 2.1538023948669434, 3.9841244220733643, 0, 91.39302131010875, -2.1338067054748535, 2.8433163166046143, 0, 93.60010401183916, 1.3892855644226074, -1.5344338417053223, 1, 325.77147156469397, 5.37094783782959, 0.5402934551239014, 1, 87.36825937672452, -1.1113907098770142, -0.31958580017089844, 0, 14.822944861528379, 3.8037590980529785, 1.8445496559143066, 0, 280.779154471958, -1.0020865201950073, -0.2703554630279541, 0, 12.467543724686887, 2.1102166175842285, 2.1807572841644287, 0, 179.7427225067583, -0.6240844130516052, -0.2827993333339691, 1, 123.6777305477976, 1.4706132411956787, -1.5239405632019043, 1, 62.36985175447399, -1.6129392385482788, -0.7079445719718933, 1, 181.2901639773627, -1.6667282581329346, -3.7529380321502686, 0, 14.093247507533691, 2.566265106201172, 0.20766335725784302, 1, 212.73129702850366, -0.5350686311721802, 3.655214309692383, 0, 90.18403172448804, 2.7695703506469727, 1.4025629758834839, 0, 301.1912795009019, -3.8380751609802246, 2.9306087493896484, 0, 276.0631013725601, -2.8070719242095947, -2.230607748031616, 1, 1.537710776182354, -1.6058628559112549, -1.2327741384506226, 0, 359.87127542782963, 3.311872720718384, -3.7691867351531982, 1, 358.62198834093397, -4.677587032318115, 1.9377728700637817, 0, 104.0816614923145, 2.7909979820251465, -1.8391573429107666, 0, 290.5553761107699, -1.932698130607605, 1.1617510318756104, 0, 246.53966220296832, -2.9844584465026855, 3.0632338523864746, 0, 272.37426991137136, -4.193409442901611, 2.826467275619507, 1, 181.42348052184224, -0.5240727663040161, 2.4276583194732666, 0, 220.37606411364231, -1.4962594509124756, -1.932018756866455, 0, 0.9798189605281387, -3.5269007682800293, -3.842702865600586, 0, 0.10056791799508939, 0.8776018023490906, 4.919729232788086, 0, 43.77123213381575, 1.6866464614868164, 1.712049126625061, 1, 72.17318008264591, -1.7034744024276733, 0.20458677411079407, 0, 42.41764141259288, 1.16348397731781, 0.8506182432174683, 0, 97.79687018091577, -3.523622989654541, -3.023449420928955, 0, 286.01321665407835, -0.013983815908432007, 0.6740539073944092, 1, 178.97859333851687, 1.011199712753296, 0.5223840475082397, 1, 82.26451487467894, 1.0656101703643799, -2.8661632537841797, 0, 97.23690060966965, 0.47832855582237244, 2.287458896636963, 0, 275.4515147940725, 0.08397042751312256, 2.5489561557769775, 0, 5.853676305378977, -2.41149640083313, 4.7202301025390625, 1, 89.84426020912757, 1.0441622734069824, 1.0218191146850586, 1, 20.89024698522807, -2.1864662170410156, 3.7620737552642822, 0, 213.87790818489665, -1.13359797000885, 0.15711438655853271, 0, 294.9214332521224, -0.8745731115341187, -0.5538784265518188, 1, 270.170288084823, -0.30008819699287415, -3.0085010528564453, 1, 268.3589815955007, 0.6017334461212158, 0.10211408883333206, 1, 359.62951853364456, -2.522864818572998, 3.5594797134399414, 1, 19.496182342429215, -1.5052179098129272, 0.7933467626571655, 0, 88.76843686524727, -2.084130048751831, 1.0204615592956543, 1, 359.5860600376021, 3.4548099040985107, -5.885735988616943, 1, 297.647718094174, 4.43270206451416, -2.890636682510376, 0, 39.08589710720983, 0.4887070655822754, 2.6409764289855957, 0, 179.74244355036893, 0.9105468988418579, -4.950046062469482, 1, 150.21547530278778, 1.7990429401397705, -2.5312092304229736, 1, 88.5305411412174, 0.6186021566390991, -0.9745205044746399, 1, 269.68485115746347, -2.0719118118286133, 1.9318674802780151, 1, 26.15497691700884, -0.8487429022789001, -3.3050856590270996, 0, 271.20465514303646, 1.3280426263809204, -5.163588523864746, 1, 306.92163838477, -3.235578775405884, -0.23159730434417725, 1, 357.6763916184892, -1.216653823852539, -1.3630906343460083, 0, 1.0704709485224466, -3.309469699859619, 5.898046970367432, 1, 298.02331187939006, -1.3187944889068604, 6.9182329177856445, 1, 268.8113736826525, 0.5162340998649597, 0.25488579273223877, 1, 178.67224209365747, 1.3074662685394287, -2.4089906215667725, 1, 347.1589350132631, 1.1237711906433105, -2.2094640731811523, 1, 150.13772267257883, -5.461596965789795, 4.667445182800293, 1, 269.773685920877, 2.8892664909362793, 3.5466394424438477, 0, 271.9420480358905, 4.8589019775390625, -3.6301159858703613, 0, 180.08072376634797, -1.0940383672714233, 0.15798038244247437, 0, 2.6562623109910604, -0.4700917601585388, -2.166661024093628, 0, 286.1518335853213, 1.2422616481781006, 0.4056287407875061, 1, 1.5541207772235457, 0.08949600160121918, 5.479471683502197, 1, 0.21405503695065925, -0.8032843470573425, 1.1026870012283325, 0, 358.7340402241108, -1.2580351829528809, 0.8723999857902527, 1, 134.9109939198651, 0.2820870876312256, -2.282655715942383, 1, 183.110911977537, -1.8077186346054077, -2.0020148754119873, 0, 312.9184553020664, -3.2496516704559326, -1.3669207096099854, 1, 358.5183692104714, 0.8757420182228088, 1.544638991355896, 0, 344.8676844252849, -1.2976150512695312, -0.5886860489845276, 1, 36.77916332504795, 1.504914402961731, 1.7345775365829468, 1, 91.27360946599258, -0.07436316460371017, -0.5434443950653076, 0, 3.796676975958829, -0.5426462888717651, 3.1517841815948486, 0, 91.15640754858747, 2.713388681411743, -0.335198312997818, 0, 180.53028346154952, 0.9672955274581909, -2.6531903743743896, 0, 21.221584052365444, 3.8023219108581543, 6.115654945373535, 1, 146.31826763017494, -0.24190106987953186, 1.188057541847229, 1, 89.99192118693621, 0.32246649265289307, 2.71256422996521, 0, 27.181239788049954, -1.2356106042861938, -1.2637487649917603, 1, 87.74804945402647, -2.0807907581329346, -2.122603416442871, 0, 51.945823063888014, 4.122188091278076, 6.610276222229004, 0, 59.18343535709449, 0.07302430272102356, 1.0445618629455566, 1, 178.54761604008934, -0.25888603925704956, -1.74147367477417, 0, 270.23682832792474, 0.7566215395927429, 0.9329351186752319, 0, 90.17639280062528, 1.277515172958374, 2.368542432785034, 1, 227.37639509685874, -1.0812950134277344, -1.5854573249816895, 0, 267.0276618283013, 1.3950247764587402, -5.76900053024292, 0, 64.92490797573352, 1.304852843284607, 3.4318957328796387, 1, 202.63411009591488, 1.7115936279296875, 0.17050188779830933, 1], [152.87877922058289, 0.008170455694198608, -1.1945573091506958, 1, 41.200266438845254, -4.806244373321533, -4.006312370300293, 1, 167.63629391945358, 1.4846200942993164, 3.7686469554901123, 0, 270.66778417778744, 1.605905532836914, -1.2995022535324097, 1, 345.271046526974, -2.5494539737701416, -0.4399448037147522, 0, 267.3617362363516, 1.2834454774856567, -1.7333629131317139, 0, 269.8763823517168, -1.4881651401519775, 2.4409637451171875, 0, 245.98002414568404, -0.49324214458465576, 4.9746222496032715, 0, 270.60057881753863, 1.1381239891052246, 1.2563565969467163, 1, 178.77965924442654, -1.1675994396209717, -1.7815622091293335, 0, 87.89801186200006, 3.8412225246429443, 2.5956356525421143, 0, 357.5783729364584, 0.606671154499054, -0.3898153305053711, 0, 269.1442251644139, -0.30423808097839355, 1.6209845542907715, 1, 156.35342990671714, -0.500913143157959, -2.4778151512145996, 1, 350.87222829235935, 1.14683198928833, -1.5955928564071655, 0, 323.04514569827194, -0.2150554656982422, 4.580379486083984, 0, 269.9268937092279, 0.1693125218153, 3.4386746883392334, 1, 91.46715763345252, -3.7951905727386475, -0.32941821217536926, 0, 90.18156408506252, 3.7631676197052, -0.6004769802093506, 1, 0.017265211284209218, 3.520387649536133, 3.7999472618103027, 1, 189.66039403215564, 2.868098735809326, -4.0951247215271, 1, 1.1876080929920205, 0.42199647426605225, -0.2842307984828949, 1, 272.29026326948906, 0.9080077409744263, 3.2697227001190186, 1, 127.40248249200337, 6.449554443359375, 3.41408371925354, 1, 183.99466750584023, -3.386709690093994, -3.9282166957855225, 0, 169.98132712981203, 0.3307725787162781, -4.224429607391357, 0, 344.2014891060639, 1.927056074142456, 4.201601028442383, 0, 270.32559871309775, -5.280389785766602, 4.6804046630859375, 0, 2.072466809688861, 1.5320508480072021, -1.4702565670013428, 1, 359.28047486975726, 0.8826894164085388, 0.4522290825843811, 1, 92.38603949905445, -2.714226722717285, 4.57764196395874, 0, 354.52559023418553, 1.3909473419189453, 5.077116012573242, 1, 271.30870341528345, -1.379101037979126, 1.0653531551361084, 0, 65.26503368866545, -7.51224422454834, 0.6308622360229492, 1, 71.3263318230623, 3.337596893310547, -2.339717388153076, 0, 89.71373320805488, -0.07851585745811462, -0.47391659021377563, 1, 178.7901092302685, 0.31785881519317627, 1.8077296018600464, 1, 179.2081475578872, -1.0003442764282227, 1.490593433380127, 1, 359.68416319428593, 1.9797884225845337, 0.6810514330863953, 1, 94.0316069531936, -0.6478455662727356, 2.380932569503784, 0, 180.42001248309518, -4.546445369720459, -0.12524443864822388, 0, 87.36131781393448, 4.31234884262085, 3.082127332687378, 1, 251.437354273979, 0.27988773584365845, -0.9809080362319946, 1, 182.75795469476327, 3.248690128326416, -2.894796848297119, 1, 95.85780480830253, 1.5927436351776123, 0.8129660487174988, 0, 269.6959018806715, -0.5498109459877014, -1.8737729787826538, 0, 269.05013322896576, -0.2743726670742035, -1.0140670537948608, 0, 358.9406346951239, 2.07716965675354, 0.4340667724609375, 0, 0.08376673082048992, -0.4916444718837738, 5.576444625854492, 0, 180.33737898649773, 2.693974018096924, 1.7386653423309326, 0, 91.62459260517123, -4.058239936828613, 0.1801467090845108, 0, 6.076177274182228, 2.0905637741088867, -3.149524688720703, 0, 272.69021269450394, 0.24398383498191833, -0.5156718492507935, 1, 4.6441604639481495, -1.1715809106826782, 1.457697868347168, 0, 351.51436996396586, -0.47298017144203186, 3.9272878170013428, 0, 86.79404676631873, 0.12927137315273285, -0.689812958240509, 0, 270.3702092028966, -0.5676143765449524, 3.9938137531280518, 1, 295.30900369307176, -1.5091888904571533, 3.9694161415100098, 0, 270.4379940096651, 6.345916748046875, 0.6465567350387573, 0, 299.8138263614649, -2.174241542816162, -4.671873092651367, 0, 233.99842612126378, -0.728750467300415, -3.160938262939453, 0, 269.85295057706594, 0.31333673000335693, -0.2701147794723511, 1, 1.031914878880541, -4.367305278778076, 2.9143152236938477, 0, 2.6699374496239483, 5.826936721801758, 7.15090274810791, 0, 104.55796899820501, -0.8173916339874268, -1.0731167793273926, 1, 357.78252122233323, -1.3060228824615479, 1.3891466856002808, 0, 358.8018036050349, 2.6334080696105957, -1.3491898775100708, 0, 46.13713422779104, -1.8215866088867188, -7.870553016662598, 1, 19.093591051217214, 3.372358560562134, -1.0649809837341309, 0, 180.91523408840015, 2.4786956310272217, 0.9396582245826721, 0, 359.7457313694598, -0.9050368666648865, -0.16045686602592468, 0, 142.66811290906614, -0.803165078163147, -1.6361286640167236, 1, 180.03583431240642, 3.2623534202575684, -3.9216527938842773, 0, 165.6984187284102, 2.132969379425049, -4.079410552978516, 1, 178.69831316076335, 1.8546302318572998, -2.058077812194824, 1, 283.29058215023656, 0.042521581053733826, -0.2550606429576874, 1, 1.6754425983172043, 1.4697848558425903, 1.7691303491592407, 1, 273.45009317849264, -2.0514583587646484, 0.8213512897491455, 0, 269.79235411291796, -1.4445736408233643, -1.2191438674926758, 0, 268.47343674451946, -2.411769151687622, -4.150338172912598, 0, 88.49557584562315, -0.2138664275407791, -0.599524736404419, 1, 270.86319925265013, -1.2424945831298828, 3.185892343521118, 1, 220.67100196439736, -0.37275558710098267, 0.4023868441581726, 1, 127.16424799902528, 0.7327858805656433, -0.8690142631530762, 0, 196.4740856162765, -0.47292473912239075, -0.6105522513389587, 0, 347.68044245436516, -1.7464954853057861, 3.348593235015869, 1, 271.6399025354337, -4.152002334594727, -2.551466941833496, 1, 196.63197013929494, 0.660304844379425, 4.105541229248047, 1, 134.58486905179132, -2.118746042251587, 3.088695764541626, 1, 268.113892220159, -0.29403215646743774, -0.5318340063095093, 1, 270.0504255274134, 4.162621021270752, 0.03553861752152443, 1, 89.66485263279762, 5.316567897796631, 0.6952872276306152, 1, 263.7920166910981, 2.2647669315338135, 1.3934054374694824, 0, 91.69854642657745, 0.7250503897666931, -0.3245762586593628, 1, 90.17160772969731, -0.7435786128044128, 3.06618595123291, 0, 120.46092521948727, 2.652642250061035, 0.032141804695129395, 1, 229.87355240518391, -0.3847227394580841, -0.8923634886741638, 1, 270.6464982199856, 1.0506783723831177, -3.8456761837005615, 1, 0.34836713818129444, -3.551510810852051, 0.6579376459121704, 0, 359.1189637491582, -1.4860488176345825, -1.8391355276107788, 1]]
    mir_stab_thld = 0.0
    grp_err_thld = 10000.0
    err_thld = 0.7
    d = 58
    def test_wrong_number_params_returns_TypeError_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.multi_align_stability()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.multi_align_stability()
        self.assertEqual(str(cm_new.exception), "multi_align_stability() takes at least 1 argument (0 given)")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_pickle_value(self):
        return_new = fu.multi_align_stability(ali_params=self.ali_params, mir_stab_thld=self.mir_stab_thld, grp_err_thld=self.grp_err_thld, err_thld=self.err_thld, print_individual=False, d=self.d)
        return_old = oldfu.multi_align_stability(ali_params=self.ali_params, mir_stab_thld=self.mir_stab_thld, grp_err_thld=self.grp_err_thld, err_thld=self.err_thld, print_individual=False, d=self.d)
        self.assertTrue(numpy_array_equal(return_new, ([], 0.11, 21.663624836876036)))
        self.assertTrue(numpy_array_equal(return_new, return_old))

    def test_ali_params_emptyList_returns_IndexError(self):
        with self.assertRaises(IndexError) as cm_new:
            fu.multi_align_stability(ali_params=[], mir_stab_thld=self.mir_stab_thld, grp_err_thld=self.grp_err_thld, err_thld=self.err_thld, print_individual=False, d=self.d)
        with self.assertRaises(IndexError) as cm_old:
            oldfu.multi_align_stability(ali_params=[], mir_stab_thld=self.mir_stab_thld, grp_err_thld=self.grp_err_thld, err_thld=self.err_thld, print_individual=False, d=self.d)
        self.assertEqual(str(cm_new.exception), "list index out of range")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_grp_err_thld_is0(self):
        return_new = fu.multi_align_stability(ali_params=self.ali_params, mir_stab_thld=self.mir_stab_thld, grp_err_thld=self.grp_err_thld, err_thld=0, print_individual=False, d=self.d)
        return_old = oldfu.multi_align_stability(ali_params=self.ali_params, mir_stab_thld=self.mir_stab_thld, grp_err_thld=self.grp_err_thld, err_thld=0, print_individual=False, d=self.d)
        self.assertTrue(numpy_array_equal(return_new, ([], 0.11, 21.663624836876036)))
        self.assertTrue(numpy_array_equal(return_new, return_old))

    def test_d_is0(self):
        return_new = fu.multi_align_stability(ali_params=self.ali_params, mir_stab_thld=self.mir_stab_thld, grp_err_thld=self.grp_err_thld, err_thld=self.err_thld, print_individual=False, d=0)
        return_old = oldfu.multi_align_stability(ali_params=self.ali_params, mir_stab_thld=self.mir_stab_thld, grp_err_thld=self.grp_err_thld, err_thld=self.err_thld, print_individual=False, d=0)
        self.assertTrue(numpy_array_equal(return_new, ([], 0.11, 2.5553015852192766)))
        self.assertTrue(numpy_array_equal(return_new, return_old))




"""
@unittest.skip("adnan's tests")
class Test_lib_pixel_error_compare(unittest.TestCase):

    def test_pixel_error_2D_true_should_return_equal_objects(self):
        print("Hello testing")
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/pixel_error/pixel_error.pixel_error_2D")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum[0])

        (ali_params1, ali_params2, r) = argum[0]

        return_new = fu.pixel_error_2D(ali_params1, ali_params2, r)
        return_old = oldfu.pixel_error_2D(ali_params1, ali_params2, r)


        self.assertEqual(return_new, return_old)


    def test_max_3D_pixel_error_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/pixel_error/pixel_error.max_3D_pixel_error")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum[0])

        (t1, t2, r) = argum[0]

        return_new = fu.max_3D_pixel_error(t1, t2, r)
        return_old = oldfu.max_3D_pixel_error(t1, t2, r)

        self.assertEqual(return_new, return_old)


    def test_angle_ave_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/pixel_error/pixel_error.angle_ave")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum[0][0])

        (angle1) = argum[0][0]

        return_new = fu.angle_ave(angle1)
        return_old = oldfu.angle_ave(angle1)

        self.assertEqual(return_new, return_old)


    def test_angle_diff_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/pixel_error/pixel_error.angle_diff")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum[0])

        (angle1,angle2) = argum[0]

        return_new = fu.angle_diff(angle1,angle2)
        return_old = oldfu.angle_diff(angle1,angle2)

        self.assertEqual(return_new, return_old)


    def test_angle_diff_sym_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/pixel_error/pixel_error.angle_diff_sym")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum[0])

        (angle1,angle2, simi) = argum[0]

        return_new = fu.angle_diff_sym(angle1,angle2,simi)
        return_old = oldfu.angle_diff_sym(angle1,angle2,simi)

        self.assertEqual(return_new, return_old)


    def test_align_diff_params_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/pixel_error/pixel_error.align_diff_params")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum[0])

        (ali_params1, ali_params2) = argum[0]

        return_new = fu.align_diff_params(ali_params1, ali_params2)
        return_old = oldfu.align_diff_params(ali_params1, ali_params2)

        self.assertEqual(return_new, return_old)



    def test_multi_align_stability_true_should_return_equal_objects(self):
        filepath = os.path.join(ABSOLUTE_PATH, "pickle files/pixel_error/pixel_error.multi_align_stability")
        with open(filepath, 'rb') as rb:
            argum = pickle.load(rb)

        print(argum[0])

        (ali_params, mir_stab_thld, grp_err_thld, err_thld,print_individual, d) = argum[0]

        return_new = fu.multi_align_stability(ali_params, mir_stab_thld, grp_err_thld, err_thld, print_individual, d)
        return_old = oldfu.multi_align_stability(ali_params, mir_stab_thld, grp_err_thld, err_thld, print_individual, d)

        self.assertEqual(return_new, return_old)
"""





if __name__ == '__main__':
    unittest.main()
