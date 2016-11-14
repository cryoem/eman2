#!/usr/bin/env python

#
# Author: Piotr Pawliczek, 09/06/2012 
# Copyright (c) 2000-2006 The University of Texas - Houston Medical School
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

import unittest
from optparse import OptionParser

IS_TEST_EXCEPTION = False

# ====================================================================================================================
class TestMultiAlignStability(unittest.TestCase):
    """this is unit test for multi_align_stability(...) from pixel_error.py"""

    # ======================= TEST FOR 2 GROUPS OF TRANSFORMATIONS (L=2, N>=1)
    
    # a0, a1 - list of parameters (alpha0, sx0, sy0, mirror0, aplha1, sx1, sy1, mirror1, alpha2, ...)
    # d - diameter
    def internal_test_conf(self, a1, a2, d):
        from EMAN2 import *
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
        for i in xrange(N):
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
        """test: L=2,  N=2,  without mirror .........................."""
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
        """test: L=2,  N=2,  with mirror .........................."""
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
        """test: L=2,  N=4,  without mirror .........................."""
        a1 = [193.0, 0.0, 0.0,   122.0, 0.0, 0.0,   223.0, 0.0, 0.0,   213.0, 0.0, 0.0]
        a2 = [ 72.0, 0.0, 0.0,   164.0, 0.0, 0.0,   323.0, 0.0, 0.0,   117.0, 0.0, 0.0]
        self.internal_test_conf_without_mirror(a1, a2)
        a1 = [193.1,  0.0, -2.0,   122.1, 3.0, -3.0,   76.1, 5.0, -1.0,   291.8, 0.2, 0.0]
        a2 = [ 72.2, -1.0,  4.0,   164.2, 3.0, -1.0,   56.1, 1.0, -7.0,   193.1, 0.3, 1.0]
        self.internal_test_conf_without_mirror(a1, a2)

    def test_four_images_with_mirrors(self):
        """test: L=2,  N=4,  with mirror .........................."""
        a1 = [193.0, 0.0, 0.0,   122.0, 0.0, 0.0,   223.0, 0.0, 0.0,   213.0, 0.0, 0.0]
        a2 = [ 72.0, 0.0, 0.0,   164.0, 0.0, 0.0,   323.0, 0.0, 0.0,   117.0, 0.0, 0.0]
        self.internal_test_conf_with_mirror(a1, a2)
        a1 = [193.1,  0.0, -2.0,   122.1, 3.0, -3.0,   76.1, 5.0, -1.0,   291.8, 0.2, 0.0]
        a2 = [ 72.2, -1.0,  4.0,   164.2, 3.0, -1.0,   56.1, 1.0, -7.0,   193.1, 0.3, 1.0]
        self.internal_test_conf_with_mirror(a1, a2)

    # ======================= TESTS FOR 10 TRANSFORMATIONS (L=2, N=10)
    
    def test_ten_images_without_mirrors(self):
        """test: L=2,  N=10,  without mirror .........................."""
        from random import random
        a1 = [193.1,  0.0, -2.0,   122.1, 3.0, -3.0,    44.6, 6.0, -3.0,    23.1, 2.0,  1.0,    93.1, -2.0,  4.2,   223.1, 2.0,  2.1,    52.1,  5.1, -3.0,   213.3, 2.0, -1.2,    76.1, 5.0, -1.0,   291.8, 0.2, 0.0]
        a2 = [ 72.2, -1.0,  4.0,   164.2, 3.0, -1.0,   196.1, 2.0, -3.0,   239.1, 2.0, -3.0,   356.2,  3.1,  2.1,   323.1, 0.2, -0.8,   303.1, -2.1,  2.0,   117.6, 2.0,  2.0,    56.1, 1.0, -7.0,   193.1, 0.3, 1.0]
        self.internal_test_conf_without_mirror(a1, a2)

    def test_ten_images_with_mirrors(self):
        """test: L=2,  N=10,  with mirror .........................."""
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
