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
class TestCorrelationFunctions(unittest.TestCase):
    """this is unit test for [acs]cf*(...) from fundamentals.py"""

    def internal_correlation(self, A, B, circulant, center): # A, B - images, circulant - bool (False - zero padded), center - bool
        from EMAN2 import EMData
        from utilities import model_blank
        from fundamentals import cyclic_shift
        anx = A.get_xsize()
        any = A.get_ysize()
        anz = A.get_zsize()
        self.assertEqual(anx, B.get_xsize())
        self.assertEqual(any, B.get_ysize())
        self.assertEqual(anz, B.get_zsize())
        snx = 2*anx
        sny = 2*any
        snz = 2*anz
        S = model_blank(snx, sny, snz)
        if circulant:
            tx = snx
            ty = sny
            tz = snz
        else:
            tx = anx
            ty = any
            tz = anz
        for x in xrange(tx):
            for y in xrange(ty):
                for z in xrange(tz):
                    S.set_value_at(x, y, z, A.get_value_at( (x)%anx, (y)%any, (z)%anz ))
        if center:
            S = cyclic_shift(S, anx/2, any/2, anz/2)
        R = model_blank(anx, any, anz)
        for x in xrange(anx):
            for y in xrange(any):
                for z in xrange(anz):
                    s = 0.0
                    for x2 in xrange(anx):
                        for y2 in xrange(any):
                            for z2 in xrange(anz):
                                s += S.get_value_at(x+x2, y+y2, z+z2) * B.get_value_at(x2, y2, z2)
                    R.set_value_at(x, y, z, s)
        return R

    def internal_assert_almostEquals(self, A, B):
        from EMAN2 import EMData
        nx = A.get_xsize()
        ny = A.get_ysize()
        nz = A.get_zsize()
        self.assertEqual(nx, B.get_xsize())
        self.assertEqual(ny, B.get_ysize())
        self.assertEqual(nz, B.get_zsize())
        for x in xrange(nx):
            for y in xrange(ny):
                for z in xrange(nz):
                    delta = abs(A.get_value_at(x,y,z)) / 100.0 # allowed error: 1% of value
                    if delta < 0.001:
                        delta = 0.001
                    self.assertAlmostEqual(A.get_value_at(x,y,z), B.get_value_at(x,y,z), delta=delta)

    def internal_check_ccf(self, A, B, AB_circ, AB_zero, AB_circ_center, AB_zero_center):
        from EMAN2 import EMData
        from global_def import Util
        
        R_circ = self.internal_correlation(A, B, True , False)
        R_zero = self.internal_correlation(A, B, False, False)
        R_circ_center = self.internal_correlation(A, B, True , True)
        R_zero_center = self.internal_correlation(A, B, False, True)
        
        self.internal_assert_almostEquals( R_circ, AB_circ )
        self.internal_assert_almostEquals( R_zero, AB_zero )
        self.internal_assert_almostEquals( R_circ_center, AB_circ_center )
        self.internal_assert_almostEquals( R_zero_center, AB_zero_center )

    def internal_test_image(self, nx, ny=1, nz=1):
        from EMAN2 import EMData
        e = EMData()
        e.set_size(nx, ny, nz)
        e.process_inplace("testimage.tomo.objects")
        return  e

    def internal_test_image2(self, nx, ny=1, nz=1):
        from EMAN2 import EMData, display
        from fundamentals import cyclic_shift, mirror
        e = EMData()
        e.set_size(nx, ny, nz)
        e.process_inplace("testimage.tomo.objects")
        e = cyclic_shift(e, nx/2, ny/3, nz/5)
        e = mirror(e)
        return  e

    # ======================= TESTS FOR acf* functions

    def test_acf_circle_2D_20x30(self):
        """test acf*: circle 2D, 20x30.........................."""
        from utilities import model_circle
        from fundamentals import acf, acfp
        A = model_circle(7, 20, 30)
        self.internal_check_ccf(A, A, acf(A,False), acfp(A,False), acf(A), acfp(A))

    def test_acf_circle_2D_21x31(self):
        """test acf*: circle 2D, 21x31.........................."""
        from utilities import model_circle
        from fundamentals import acf, acfp
        A = model_circle(7, 21, 31)
        self.internal_check_ccf(A, A, acf(A,False), acfp(A,False), acf(A), acfp(A))

    def test_acf_circle_2D_31x20(self):
        """test acf*: circle 2D, 31x20.........................."""
        from utilities import model_circle
        from fundamentals import acf, acfp
        A = model_circle(7, 31, 20)
        self.internal_check_ccf(A, A, acf(A,False), acfp(A,False), acf(A), acfp(A))

    def test_acf_objects_2D_20x30(self):
        """test acf*: objects 2D, 20x30.........................."""
        from fundamentals import acf, acfp
        A = self.internal_test_image(20, 30)
        self.internal_check_ccf(A, A, acf(A,False), acfp(A,False), acf(A), acfp(A))

    def test_acf_objects_2D_21x31(self):
        """test acf*: objects 2D, 21x31.........................."""
        from fundamentals import acf, acfp
        A = self.internal_test_image(21, 31)
        self.internal_check_ccf(A, A, acf(A,False), acfp(A,False), acf(A), acfp(A))

    def test_acf_objects_2D_31x20(self):
        """test acf*: objects 2D, 31x20.........................."""
        from fundamentals import acf, acfp
        A = self.internal_test_image(31, 20)
        self.internal_check_ccf(A, A, acf(A,False), acfp(A,False), acf(A), acfp(A))

    # ======================= TESTS FOR ccf* functions

    def test_ccf_circle_2D_20x30(self):
        """test ccf*: circle 2D, 20x30.........................."""
        from utilities import model_circle
        from fundamentals import ccf, ccfp
        A = model_circle(7, 20, 30)
        B = model_circle(4, 20, 30)
        self.internal_check_ccf(A, B, ccf(A,B,False), ccfp(A,B,False), ccf(A,B), ccfp(A,B))

    def test_ccf_circle_2D_21x31(self):
        """test ccf*: circle 2D, 21x31.........................."""
        from utilities import model_circle
        from fundamentals import ccf, ccfp
        A = model_circle(7, 21, 31)
        B = model_circle(4, 21, 31)
        self.internal_check_ccf(A, B, ccf(A,B,False), ccfp(A,B,False), ccf(A,B), ccfp(A,B))

    def test_ccf_circle_2D_31x20(self):
        """test ccf*: circle 2D, 31x20.........................."""
        from utilities import model_circle
        from fundamentals import ccf, ccfp
        A = model_circle(7, 31, 20)
        B = model_circle(4, 31, 20)
        self.internal_check_ccf(A, B, ccf(A,B,False), ccfp(A,B,False), ccf(A,B), ccfp(A,B))

    def test_ccf_objects_2D_20x30(self):
        """test ccf*: objects 2D, 20x30.........................."""
        from fundamentals import ccf, ccfp
        A = self.internal_test_image(20, 30)
        B = self.internal_test_image2(20, 30)
        self.internal_check_ccf(A, B, ccf(A,B,False), ccfp(A,B,False), ccf(A,B), ccfp(A,B))

    def test_ccf_objects_2D_21x31(self):
        """test ccf*: objects 2D, 21x31.........................."""
        from fundamentals import ccf, ccfp
        A = self.internal_test_image(21, 31)
        B = self.internal_test_image2(21, 31)
        self.internal_check_ccf(A, B, ccf(A,B,False), ccfp(A,B,False), ccf(A,B), ccfp(A,B))

    def test_ccf_objects_2D_31x20(self):
        """test ccf*: objects 2D, 31x20.........................."""
        from fundamentals import ccf, ccfp
        A = self.internal_test_image(31, 20)
        B = self.internal_test_image2(31, 20)
        self.internal_check_ccf(A, B, ccf(A,B,False), ccfp(A,B,False), ccf(A,B), ccfp(A,B))

    # ======================= TESTS FOR cnv* functions

    def test_cnv_circle_2D_20x30(self):
        """test cnv*: circle 2D, 20x30.........................."""
        from utilities import model_circle
        from fundamentals import cnv, cnvp, mirror
        A = model_circle(7, 20, 30)
        B = model_circle(4, 20, 30)
        self.internal_check_ccf(A, mirror(mirror(B,'x'),'y'), cnv(A,B,False), cnvp(A,B,False), cnv(A,B), cnvp(A,B))

    def test_cnv_circle_2D_21x31(self):
        """test cnv*: circle 2D, 21x31.........................."""
        from utilities import model_circle
        from fundamentals import cnv, cnvp, mirror
        A = model_circle(7, 21, 31)
        B = model_circle(4, 21, 31)
        self.internal_check_ccf(A, mirror(mirror(B,'x'),'y'), cnv(A,B,False), cnvp(A,B,False), cnv(A,B), cnvp(A,B))

    def test_cnv_circle_2D_31x20(self):
        """test cnv*: circle 2D, 31x20.........................."""
        from utilities import model_circle
        from fundamentals import cnv, cnvp, mirror
        A = model_circle(7, 31, 20)
        B = model_circle(4, 31, 20)
        self.internal_check_ccf(A, mirror(mirror(B,'x'),'y'), cnv(A,B,False), cnvp(A,B,False), cnv(A,B), cnvp(A,B))

    def test_cnv_objects_2D_20x30(self):
        """test cnv*: objects 2D, 20x30.........................."""
        from fundamentals import cnv, cnvp, mirror
        A = self.internal_test_image(20, 30)
        B = self.internal_test_image2(20, 30)
        self.internal_check_ccf(A, mirror(mirror(B,'x'),'y'), cnv(A,B,False), cnvp(A,B,False), cnv(A,B), cnvp(A,B))

    def test_cnv_objects_2D_21x31(self):
        """test cnv*: objects 2D, 21x31.........................."""
        from fundamentals import cnv, cnvp, mirror
        A = self.internal_test_image(21, 31)
        B = self.internal_test_image2(21, 31)
        self.internal_check_ccf(A, mirror(mirror(B,'x'),'y'), cnv(A,B,False), cnvp(A,B,False), cnv(A,B), cnvp(A,B))

    def test_cnv_objects_2D_31x20(self):
        """test cnv*: objects 2D, 31x20.........................."""
        from fundamentals import cnv, cnvp, mirror
        A = self.internal_test_image(31, 20)
        B = self.internal_test_image2(31, 20)
        self.internal_check_ccf(A, mirror(mirror(B,'x'),'y'), cnv(A,B,False), cnvp(A,B,False), cnv(A,B), cnvp(A,B))


def test_main():
    from EMAN2 import Log
    p = OptionParser()
    p.add_option('--t', action='store_true', help='test exception', default=False )
    global IS_TEST_EXCEPTION
    opt, args = p.parse_args()
    if opt.t:
        IS_TEST_EXCEPTION = True
    Log.logger().set_level(-1)  #perfect solution for quenching the Log error information, thank Liwei
    suite = unittest.TestLoader().loadTestsFromTestCase(TestCorrelationFunctions)
    unittest.TextTestRunner(verbosity=2).run(suite)

if __name__ == '__main__':
    test_main()
