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
from test import test_support
import testlib

class TestTransform(unittest.TestCase):
	"""this is the unit test for Transform3D class"""
		
	def test_get_rotation(self):
		"""test rotation ...................................."""
		alt = 1.45232928554
		az = -0.60170830102
		phi = 0

		t = Transform3D(az, alt, phi)
		rot = t.get_rotation(EULER_EMAN)

		#testlib.assertfloat(self, az, float(rot["az"]))
		self.assertAlmostEqual(az, rot["az"], 3)
		testlib.assertfloat(self, alt, float(rot["alt"]))
		testlib.assertfloat(self, phi, float(rot["phi"]))
		
	def test_trans_after_rotation(self):
		"""test translation after rotation .................."""
		alt = 1.45232928554
		az = -0.60170830102
		phi = 0
		
		t = Transform3D((0.0,0.0,0.0), az, alt, phi,(1.0,2.0,3.0),)
		rot = t.get_rotation(EULER_EMAN)
		tran = t.get_posttrans()
		
		#testlib.assertfloat(self, az, float(rot["az"]))
		self.assertAlmostEqual(az, rot["az"], 3)
		testlib.assertfloat(self, alt, float(rot["alt"]))
		testlib.assertfloat(self, phi, float(rot["phi"]))
		testlib.assertfloat(self, tran.at(0), 1.0)
		testlib.assertfloat(self, tran.at(1), 2.0)
		testlib.assertfloat(self, tran.at(2), 3.0)
		
	def no_test_trans_before_rotation(self):
		"""test translation before rotation ................."""
		t = Transform3D(Transform3D.EulerType.EMAN, {'az':-0.60170830102, 'alt':1.45232928554,'phi':0})
		
	def test_pre_post_trans_rotation(self):
		"""test pretrans,rotation,post trans constructor ...."""
		alt = 1.45232928554
		az = -0.60170830102
		phi = 0
		t = Transform3D((1.0,2.0,3.0), az, alt, phi, (4.0,5.0,6.0))
		
		rot = t.get_rotation()    #default argument is EULER_EMAN
		#testlib.assertfloat(self, az, float(rot["az"]))
		self.assertAlmostEqual(az, rot["az"], 3)
		testlib.assertfloat(self, alt, float(rot["alt"]))
		testlib.assertfloat(self, phi, float(rot["phi"]))
		
		
		tran = t.get_posttrans()
		self.assertAlmostEqual(tran.at(0), 4.0, 3)
		self.assertAlmostEqual(tran.at(1), 5.0, 3)
		self.assertAlmostEqual(tran.at(2), 6.0, 3)
		
		tran = t.get_pretrans()
		self.assertAlmostEqual(tran.at(0), 1.0, 3)
		self.assertAlmostEqual(tran.at(1), 2.0, 3)
		self.assertAlmostEqual(tran.at(2), 3.0, 3)
		
	def test_set_posttrans(self):
		"""test set/get_posttrans() function ................"""
		t = Transform3D()
		t.set_posttrans((1.0,2.0,3.0))
		
		tran = t.get_posttrans()
		self.assertAlmostEqual(tran.at(0), 1.0, 3)
		self.assertAlmostEqual(tran.at(1), 2.0, 3)
		self.assertAlmostEqual(tran.at(2), 3.0, 3)
		
	def test_apply_scale(self):
		"""test apply_scale() function ......................"""
		alt = 1.45232928554
		az = -0.60170830102
		phi = 0
		t = Transform3D((1.0,2.0,3.0), az, alt, phi, (4.0,5.0,6.0))
		
		t.apply_scale(0.4)
		self.assertAlmostEqual(t.get_scale(), 0.4, 3)
		
		#tran = t.get_posttrans()
		#self.assertAlmostEqual(tran.at(0), 4.0*1.5, 3)
		#self.assertAlmostEqual(tran.at(1), 5.0*1.5, 3)
		#self.assertAlmostEqual(tran.at(2), 6.0*1.5, 3)

	def test_set_scale(self):
		"""test set_scale() function ........................"""
		alt = 1.45232928554
		az = -0.60170830102
		phi = 0
		t = Transform3D((1.0,2.0,3.0), az, alt, phi, (4.0,5.0,6.0))
		
		t.set_scale(2.5)
		self.assertAlmostEqual(t.get_scale(), 2.5, 3)
		
	def test_orthogonalize(self):
		"""test orthogonalize() function ...................."""
		alt = 1.45232928554
		az = -0.60170830102
		phi = 0
		t = Transform3D((1.0,2.0,3.0), az, alt, phi, (4.0,5.0,6.0))
		
		t.set_scale(2.3)
		self.assertAlmostEqual(t.get_scale(), 2.3, 3)
		t.orthogonalize()
		self.assertAlmostEqual(t.get_scale(), 1.0, 3)

	def test_get_sym(self):
		"""test get_sym() function .........................."""
		alt = 1.45232928554
		az = -0.60170830102
		phi = 0
		t = Transform3D((1.0,2.0,3.0), az, alt, phi, (4.0,5.0,6.0))
		
		tt = t.get_sym('c1', 1)
	
	def test_set_center(self):
		"""test set_center/pretrans() function .............."""
		alt = 1.45232928554
		az = -0.60170830102
		phi = 0
		t = Transform3D((1.0,2.0,3.0), az, alt, phi, (4.0,5.0,6.0))
		
		t.set_center((2.0,3.0,4.0))
		t.set_pretrans((1.1,2.2,3.3))
	
	def test_to_identity(self):
		"""test to_identity() function ......................"""
		alt = 1.45232928554
		az = -0.60170830102
		phi = 0
		t = Transform3D((1.0,2.0,3.0), az, alt, phi, (4.0,5.0,6.0))
		
		t.to_identity()
		self.assertEqual(t.is_identity(), True)
		for i in range(3):
			col = t.get_matrix3_col(i)
			for j in range(3):
				if j==i:
					self.assertAlmostEqual(col.at(j), 1.0, 3)
	def test_symc_reduce(self):
		"""test csym reduce ................................."""
		syms = []
			
		syms.append(Symmetries.get("c",{"nsym":2}))
		syms.append(Symmetries.get("c",{"nsym":3}))
		syms.append(Symmetries.get("c",{"nsym":4}))
		syms.append(Symmetries.get("c",{"nsym":4}))
		syms.append(Symmetries.get("c",{"nsym":5}))
		syms.append(Symmetries.get("c",{"nsym":6}))
		syms.append(Symmetries.get("c",{"nsym":7}))
		syms.append(Symmetries.get("c",{"nsym":8}))
		syms.append(Symmetries.get("c",{"nsym":9}))
		syms.append(Symmetries.get("c",{"nsym":10}))
		syms.append(Symmetries.get("c",{"nsym":11}))
		syms.append(Symmetries.get("c",{"nsym":12}))
		for sym in syms:
			n = sym.get_nsym()
			azmax = 360.0/n
			eulers = sym.gen_orientations("eman",{"delta":8})
			for euler in eulers:
				az = euler.get_rotation()["az"]
				alt = euler.get_rotation()["alt"]
				#print az,alt,n
				for i in range(1,n):
					#print i
					T = Transform3D(az,alt,0)
					T1 = sym.get_sym(i)
					T2 = T*T1
					#print i,i
					A = sym.reduce(T2,0)
					result = A.get_rotation()
					azsoln = result["az"] %azmax
					if ( azsoln < azmax and (azmax - azsoln) < 0.001 ): azsoln = 0.0
					#print i,i,i
					#print result["az"],az,result["az"] % (720.0/n),az% (720.0/n),720.0/n
					self.assertAlmostEqual(azsoln,az%azmax, 3)
					self.assertAlmostEqual(result["alt"],alt, 3)
	def test_symd_reduce(self):
		"""test dsym reduce ................................."""
		syms = []
		syms.append(Symmetries.get("d",{"nsym":1}))
		syms.append(Symmetries.get("d",{"nsym":2}))
		syms.append(Symmetries.get("d",{"nsym":3}))
		syms.append(Symmetries.get("d",{"nsym":4}))
		syms.append(Symmetries.get("d",{"nsym":4}))
		syms.append(Symmetries.get("d",{"nsym":5}))
		syms.append(Symmetries.get("d",{"nsym":6}))
		syms.append(Symmetries.get("d",{"nsym":7}))
		syms.append(Symmetries.get("d",{"nsym":8}))
		syms.append(Symmetries.get("d",{"nsym":9}))
		syms.append(Symmetries.get("d",{"nsym":10}))
		syms.append(Symmetries.get("d",{"nsym":11}))
		syms.append(Symmetries.get("d",{"nsym":12}))
		for sym in syms:
			n = sym.get_nsym()
			azmax = 720.0/n
			
			eulers = sym.gen_orientations("eman",{"delta":8})
			for euler in eulers:
				az = euler.get_rotation()["az"]
				alt = euler.get_rotation()["alt"]
				#print az,alt,n
				for i in range(1,n):
					#print i
					T = Transform3D(az,alt,0)
					T1 = sym.get_sym(i)
					T2 = T*T1
					#print i,i
					A = sym.reduce(T2,0)
					result = A.get_rotation()
					azsoln = result["az"] %azmax
					if ( azsoln < azmax and (azmax - azsoln) < 0.001 ): azsoln = 0.0
					#print i,i,i
					#print result["az"],az,result["az"] % (720.0/n),az% (720.0/n),720.0/n
					self.assertAlmostEqual(azsoln,az%azmax, 3)
					self.assertAlmostEqual(result["alt"],alt, 3)
	def test_symtet_reduce(self):
		"""test tetsym reduce ..............................."""
		syms = []
		syms.append(Symmetries.get("tet",{}))
		for sym in syms:
			n = sym.get_nsym()
			azmax = 120.0
			
			eulers = sym.gen_orientations("eman",{"delta":4})
			for euler in eulers:
				az = euler.get_rotation()["az"]
				alt = euler.get_rotation()["alt"]
				#print az,alt,n
				for i in range(1,n):
					#print i
					T = Transform3D(az,alt,0)
					T1 = sym.get_sym(i)
					T2 = T*T1
					#print i,i
					A = sym.reduce(T2,0)
					result = A.get_rotation()
					azsoln = result["az"] %azmax
					if ( azsoln < azmax and (azmax - azsoln) < 0.001 ): azsoln = 0.0
					#print i,i,i
					#print result["az"],az,result["az"] % (720.0/n),az% (720.0/n),720.0/n
					self.assertAlmostEqual(azsoln,az%azmax, 3)
					self.assertAlmostEqual(result["alt"],alt, 3)
	def test_symoct_reduce(self):
		"""test octsym reduce ..............................."""
		syms = []
		syms.append(Symmetries.get("oct",{}))
		for sym in syms:
			n = sym.get_nsym()
			azmax = 90.0
			
			eulers = sym.gen_orientations("eman",{"delta":4})
			for euler in eulers:
				az = euler.get_rotation()["az"]
				alt = euler.get_rotation()["alt"]
				#print az,alt,n
				for i in range(1,n):
					#print i
					T = Transform3D(az,alt,0)
					T1 = sym.get_sym(i)
					T2 = T*T1
					#print i,i
					A = sym.reduce(T2,0)
					result = A.get_rotation()
					azsoln = result["az"] %azmax
					if ( azsoln < azmax and (azmax - azsoln) < 0.001 ): azsoln = 0.0
					#print i,i,i
					#print result["az"],az,result["az"] % (720.0/n),az% (720.0/n),720.0/n
					self.assertAlmostEqual(azsoln,az%azmax, 3)
					self.assertAlmostEqual(result["alt"],alt, 3)
	def test_symicos_reduce(self):
		"""test icossym reduce .............................."""
		syms = []
			
		syms.append(Symmetries.get("icos",{}))
		for sym in syms:
			n = sym.get_nsym()
			azmax = 72.0
			
			eulers = sym.gen_orientations("eman",{"delta":4})
			for euler in eulers:
				az = euler.get_rotation()["az"]
				alt = euler.get_rotation()["alt"]
				#print az,alt,n
				for i in range(1,n):
					#print i
					T = Transform3D(az,alt,0)
					T1 = sym.get_sym(i)
					T2 = T*T1
					#print i,i
					A = sym.reduce(T2,0)
					result = A.get_rotation()
					azsoln = result["az"] %azmax
					if ( azsoln < azmax and (azmax - azsoln) < 0.001 ): azsoln = 0.0
					#print i,i,i
					#print result["az"],az,result["az"] % (720.0/n),az% (720.0/n),720.0/n
					self.assertAlmostEqual(azsoln,az%azmax, 3)
					self.assertAlmostEqual(result["alt"],alt, 3)

	def no_test_angles2tfvec(self):    #no callable in Python
		"""test angles2tfvec() function ....................."""
		t = Transform3D.angles2tfvec(Transform3D.EulerType.EMAN, (1.45232928554, -0.60170830102, 0))

def test_main():
	test_support.run_unittest(TestTransform)

if __name__ == '__main__':
	test_main()
