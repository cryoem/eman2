#!/usr/bin/env python

#
# Author: Piotr Pawliczek, 10/25/2012 
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

	def internal_correlation(self, A, B, center, circulant, normalized, lag_normalization): # A, B - images, circulant - bool (False - zero padded), center - bool, normalized - bool
		from EMAN2 import EMData
		from utilities import model_blank
		from fundamentals import cyclic_shift
		from math import sqrt
		anx = A.get_xsize()
		any = A.get_ysize()
		anz = A.get_zsize()
		self.assertEqual(anx, B.get_xsize())
		self.assertEqual(any, B.get_ysize())
		self.assertEqual(anz, B.get_zsize())
		snx = 2*anx
		sny = 2*any
		snz = 2*anz
		if normalized:
			A = A.copy()
			B = B.copy()
			A.sub(A.get_attr("mean"))
			B.sub(B.get_attr("mean"))
			A.div(A.get_attr("sigma") * sqrt(anz) * sqrt(any) * sqrt(anx))
			B.div(B.get_attr("sigma") * sqrt(anz) * sqrt(any) * sqrt(anx))
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
		if lag_normalization:
			cx = anx/2
			cy = any/2
			cz = anz/2
			for x in xrange(anx):
				x_center = abs(x-cx)
				x_lag = 1 + (x_center * 1.0) / (anx - x_center)
				for y in xrange(any):
					y_center = abs(y-cy)
					y_lag = 1 + (y_center * 1.0) / (any - y_center)
					for z in xrange(anz):
						z_center = abs(z-cz)
						z_lag = 1 + (z_center * 1.0) / (anz - z_center)
						R.set_value_at(x, y, z, R.get_value_at(x,y,z) * x_lag * y_lag * z_lag )
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

	def internal_check_ccf_center(self, A, B, AB_circ, AB_circ_norm, AB_zero, AB_zero_norm, AB_lag, AB_lag_norm, center):
		from EMAN2 import EMData
		from global_def import Util
		
		R_circ          = self.internal_correlation(A, B, center, True , False, False)
		R_zero          = self.internal_correlation(A, B, center, False, False, False)
		R_circ_norm     = self.internal_correlation(A, B, center, True , True , False)
		R_zero_norm     = self.internal_correlation(A, B, center, False, True , False)
		R_zero_lag      = self.internal_correlation(A, B, center, False, False, True )
		R_zero_lag_norm = self.internal_correlation(A, B, center, False, True , True )
		
		self.internal_assert_almostEquals( R_circ         , AB_circ      )
		self.internal_assert_almostEquals( R_zero         , AB_zero      )
		self.internal_assert_almostEquals( R_circ_norm    , AB_circ_norm )
		self.internal_assert_almostEquals( R_zero_norm    , AB_zero_norm )
		self.internal_assert_almostEquals( R_zero_lag     , AB_lag       )
		self.internal_assert_almostEquals( R_zero_lag_norm, AB_lag_norm  )

	def internal_check_ccf(self, A, B, AB_circ, AB_circ_norm, AB_zero, AB_zero_norm, AB_lag, AB_lag_norm, cent_AB_circ, cent_AB_circ_norm, cent_AB_zero, cent_AB_zero_norm, cent_AB_lag, cent_AB_lag_norm):
		self.internal_check_ccf_center( A, B,      AB_circ,      AB_circ_norm,      AB_zero,      AB_zero_norm,      AB_lag,      AB_lag_norm, False )
		self.internal_check_ccf_center( A, B, cent_AB_circ, cent_AB_circ_norm, cent_AB_zero, cent_AB_zero_norm, cent_AB_lag, cent_AB_lag_norm, True  )

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
		from fundamentals import acf, acfn, acfnp, acfnpl, acfp, acfpl
		A = model_circle(7, 20, 30)
		self.internal_check_ccf(A, A, acf(A,False), acfn(A,False), acfp(A,False), acfnp(A,False), acfpl(A, False), acfnpl(A, False)
									, acf(A,True ), acfn(A,True ), acfp(A,True ), acfnp(A,True ), acfpl(A, True ), acfnpl(A, True ) )

	def test_acf_circle_2D_21x31(self):
		"""test acf*: circle 2D, 21x31.........................."""
		from utilities import model_circle
		from fundamentals import acf, acfn, acfnp, acfnpl, acfp, acfpl
		A = model_circle(7, 21, 31)
		self.internal_check_ccf(A, A, acf(A,False), acfn(A,False), acfp(A,False), acfnp(A,False), acfpl(A, False), acfnpl(A, False)
									, acf(A,True ), acfn(A,True ), acfp(A,True ), acfnp(A,True ), acfpl(A, True ), acfnpl(A, True ) )

	def test_acf_circle_2D_31x20(self):
		"""test acf*: circle 2D, 31x20.........................."""
		from utilities import model_circle
		from fundamentals import acf, acfn, acfnp, acfnpl, acfp, acfpl
		A = model_circle(7, 31, 20)
		self.internal_check_ccf(A, A, acf(A,False), acfn(A,False), acfp(A,False), acfnp(A,False), acfpl(A, False), acfnpl(A, False)
									, acf(A,True ), acfn(A,True ), acfp(A,True ), acfnp(A,True ), acfpl(A, True ), acfnpl(A, True ) )

	def test_acf_objects_2D_20x30(self):
		"""test acf*: objects 2D, 20x30.........................."""
		from fundamentals import acf, acfn, acfnp, acfnpl, acfp, acfpl
		A = self.internal_test_image(20, 30)
		self.internal_check_ccf(A, A, acf(A,False), acfn(A,False), acfp(A,False), acfnp(A,False), acfpl(A, False), acfnpl(A, False)
									, acf(A,True ), acfn(A,True ), acfp(A,True ), acfnp(A,True ), acfpl(A, True ), acfnpl(A, True ) )

	def test_acf_objects_2D_21x31(self):
		"""test acf*: objects 2D, 21x31.........................."""
		from fundamentals import acf, acfn, acfnp, acfnpl, acfp, acfpl
		A = self.internal_test_image(21, 31)
		self.internal_check_ccf(A, A, acf(A,False), acfn(A,False), acfp(A,False), acfnp(A,False), acfpl(A, False), acfnpl(A, False)
									, acf(A,True ), acfn(A,True ), acfp(A,True ), acfnp(A,True ), acfpl(A, True ), acfnpl(A, True ) )

	def test_acf_objects_2D_31x20(self):
		"""test acf*: objects 2D, 31x20.........................."""
		from fundamentals import acf, acfn, acfnp, acfnpl, acfp, acfpl
		A = self.internal_test_image(31, 20)
		self.internal_check_ccf(A, A, acf(A,False), acfn(A,False), acfp(A,False), acfnp(A,False), acfpl(A, False), acfnpl(A, False)
									, acf(A,True ), acfn(A,True ), acfp(A,True ), acfnp(A,True ), acfpl(A, True ), acfnpl(A, True ) )

	# ======================= TESTS FOR ccf* functions
	
	def test_ccf_circle_2D_20x30(self):
		"""test ccf*: circle 2D, 20x30.........................."""
		from utilities import model_circle
		from fundamentals import ccf, ccfn, ccfnp, ccfnpl, ccfp, ccfpl
		A = model_circle(7, 20, 30)
		B = model_circle(4, 20, 30)
		self.internal_check_ccf(A, B, ccf(A,B,False), ccfn(A,B,False), ccfp(A,B,False), ccfnp(A,B,False), ccfpl(A,B,False), ccfnpl(A,B,False)
									, ccf(A,B,True ), ccfn(A,B,True ), ccfp(A,B,True ), ccfnp(A,B,True ), ccfpl(A,B,True ), ccfnpl(A,B,True ) )
	
	def test_ccf_circle_2D_21x31(self):
		"""test ccf*: circle 2D, 21x31.........................."""
		from utilities import model_circle
		from fundamentals import ccf, ccfn, ccfnp, ccfnpl, ccfp, ccfpl
		A = model_circle(7, 21, 31)
		B = model_circle(4, 21, 31)
		self.internal_check_ccf(A, B, ccf(A,B,False), ccfn(A,B,False), ccfp(A,B,False), ccfnp(A,B,False), ccfpl(A,B,False), ccfnpl(A,B,False)
									, ccf(A,B,True ), ccfn(A,B,True ), ccfp(A,B,True ), ccfnp(A,B,True ), ccfpl(A,B,True ), ccfnpl(A,B,True ) )

	def test_ccf_circle_2D_31x20(self):
		"""test ccf*: circle 2D, 31x20.........................."""
		from utilities import model_circle
		from fundamentals import ccf, ccfn, ccfnp, ccfnpl, ccfp, ccfpl
		A = model_circle(7, 31, 20)
		B = model_circle(4, 31, 20)
		self.internal_check_ccf(A, B, ccf(A,B,False), ccfn(A,B,False), ccfp(A,B,False), ccfnp(A,B,False), ccfpl(A,B,False), ccfnpl(A,B,False)
									, ccf(A,B,True ), ccfn(A,B,True ), ccfp(A,B,True ), ccfnp(A,B,True ), ccfpl(A,B,True ), ccfnpl(A,B,True ) )

	def test_ccf_objects_2D_20x30(self):
		"""test ccf*: objects 2D, 20x30.........................."""
		from fundamentals import ccf, ccfn, ccfnp, ccfnpl, ccfp, ccfpl
		A = self.internal_test_image(20, 30)
		B = self.internal_test_image2(20, 30)
		self.internal_check_ccf(A, B, ccf(A,B,False), ccfn(A,B,False), ccfp(A,B,False), ccfnp(A,B,False), ccfpl(A,B,False), ccfnpl(A,B,False)
									, ccf(A,B,True ), ccfn(A,B,True ), ccfp(A,B,True ), ccfnp(A,B,True ), ccfpl(A,B,True ), ccfnpl(A,B,True ) )

	def test_ccf_objects_2D_21x31(self):
		"""test ccf*: objects 2D, 21x31.........................."""
		from fundamentals import ccf, ccfn, ccfnp, ccfnpl, ccfp, ccfpl
		A = self.internal_test_image(21, 31)
		B = self.internal_test_image2(21, 31)
		self.internal_check_ccf(A, B, ccf(A,B,False), ccfn(A,B,False), ccfp(A,B,False), ccfnp(A,B,False), ccfpl(A,B,False), ccfnpl(A,B,False)
									, ccf(A,B,True ), ccfn(A,B,True ), ccfp(A,B,True ), ccfnp(A,B,True ), ccfpl(A,B,True ), ccfnpl(A,B,True ) )

	def test_ccf_objects_2D_31x20(self):
		"""test ccf*: objects 2D, 31x20.........................."""
		from fundamentals import ccf, ccfn, ccfnp, ccfnpl, ccfp, ccfpl
		A = self.internal_test_image(31, 20)
		B = self.internal_test_image2(31, 20)
		self.internal_check_ccf(A, B, ccf(A,B,False), ccfn(A,B,False), ccfp(A,B,False), ccfnp(A,B,False), ccfpl(A,B,False), ccfnpl(A,B,False)
									, ccf(A,B,True ), ccfn(A,B,True ), ccfp(A,B,True ), ccfnp(A,B,True ), ccfpl(A,B,True ), ccfnpl(A,B,True ) )
	
	# ======================= TESTS FOR cnv* functions
'''
	def test_cnv_circle_2D_20x30(self):
		"""test cnv*: circle 2D, 20x30.........................."""
		from utilities import model_circle
		from fundamentals import cnv, cnvn, cnvnp, cnvnpl, cnvp, cnvpl, mirror
		A = model_circle(7, 20, 30)
		B = model_circle(4, 20, 30)
		C = mirror(mirror(B,'x'),'y')
		self.internal_check_ccf(A, C, cnv(A,B,False), cnvn(A,B,False), cnvp(A,B,False), cnvnp(A,B,False), cnvpl(A,B,False), cnvnpl(A,B,False)
									, cnv(A,B,True ), cnvn(A,B,True ), cnvp(A,B,True ), cnvnp(A,B,True ), cnvpl(A,B,True ), cnvnpl(A,B,True ) )

	def test_cnv_circle_2D_21x31(self):
		"""test cnv*: circle 2D, 21x31.........................."""
		from utilities import model_circle
		from fundamentals import cnv, cnvn, cnvnp, cnvnpl, cnvp, cnvpl, mirror
		A = model_circle(7, 21, 31)
		B = model_circle(4, 21, 31)
		C = mirror(mirror(B,'x'),'y')
		self.internal_check_ccf(A, C, cnv(A,B,False), cnvn(A,B,False), cnvp(A,B,False), cnvnp(A,B,False), cnvpl(A,B,False), cnvnpl(A,B,False)
									, cnv(A,B,True ), cnvn(A,B,True ), cnvp(A,B,True ), cnvnp(A,B,True ), cnvpl(A,B,True ), cnvnpl(A,B,True ) )

	def test_cnv_circle_2D_31x20(self):
		"""test cnv*: circle 2D, 31x20.........................."""
		from utilities import model_circle
		from fundamentals import cnv, cnvn, cnvnp, cnvnpl, cnvp, cnvpl, mirror
		A = model_circle(7, 31, 20)
		B = model_circle(4, 31, 20)
		C = mirror(mirror(B,'x'),'y')
		self.internal_check_ccf(A, C, cnv(A,B,False), cnvn(A,B,False), cnvp(A,B,False), cnvnp(A,B,False), cnvpl(A,B,False), cnvnpl(A,B,False)
									, cnv(A,B,True ), cnvn(A,B,True ), cnvp(A,B,True ), cnvnp(A,B,True ), cnvpl(A,B,True ), cnvnpl(A,B,True ) )

	def test_cnv_objects_2D_20x30(self):
		"""test cnv*: objects 2D, 20x30.........................."""
		from fundamentals import cnv, cnvn, cnvnp, cnvnpl, cnvp, cnvpl, mirror
		A = self.internal_test_image(20, 30)
		B = self.internal_test_image2(20, 30)
		C = mirror(mirror(B,'x'),'y')
		self.internal_check_ccf(A, C, cnv(A,B,False), cnvn(A,B,False), cnvp(A,B,False), cnvnp(A,B,False), cnvpl(A,B,False), cnvnpl(A,B,False)
									, cnv(A,B,True ), cnvn(A,B,True ), cnvp(A,B,True ), cnvnp(A,B,True ), cnvpl(A,B,True ), cnvnpl(A,B,True ) )

	def test_cnv_objects_2D_21x31(self):
		"""test cnv*: objects 2D, 21x31.........................."""
		from fundamentals import cnv, cnvn, cnvnp, cnvnpl, cnvp, cnvpl, mirror
		A = self.internal_test_image(21, 31)
		B = self.internal_test_image2(21, 31)
		C = mirror(mirror(B,'x'),'y')
		self.internal_check_ccf(A, C, cnv(A,B,False), cnvn(A,B,False), cnvp(A,B,False), cnvnp(A,B,False), cnvpl(A,B,False), cnvnpl(A,B,False)
									, cnv(A,B,True ), cnvn(A,B,True ), cnvp(A,B,True ), cnvnp(A,B,True ), cnvpl(A,B,True ), cnvnpl(A,B,True ) )

	def test_cnv_objects_2D_31x20(self):
		"""test cnv*: objects 2D, 31x20.........................."""
		from fundamentals import cnv, cnvn, cnvnp, cnvnpl, cnvp, cnvpl, mirror
		A = self.internal_test_image(31, 20)
		B = self.internal_test_image2(31, 20)
		C = mirror(mirror(B,'x'),'y')
		self.internal_check_ccf(A, C, cnv(A,B,False), cnvn(A,B,False), cnvp(A,B,False), cnvnp(A,B,False), cnvpl(A,B,False), cnvnpl(A,B,False)
									, cnv(A,B,True ), cnvn(A,B,True ), cnvp(A,B,True ), cnvnp(A,B,True ), cnvpl(A,B,True ), cnvnpl(A,B,True ) )
'''

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
