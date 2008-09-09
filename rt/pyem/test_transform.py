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
import math


class TestTransform(unittest.TestCase):
	"""this is the unit test for Transform class"""
	#transforms = []
	#d = {"type":"2d","alpha":self.get_angle_rand()}
	#t = Transform(d)
	#transforms.append(t)
	#d = {"type":"eman","az":self.get_angle_rand(),"alt":self.get_angle_rand(0,179),"phi":self.get_angle_rand()}
	#t = Transform(d)
	#transforms.append(t)
	#d = {"dx":0.2,"dy":3.0}
	#t = Transform(d)
	#transforms.append(t)

	def get_angle_rand(self,lower=-359,upper=359):
		return Util.get_frand(lower,upper)
	
	def test_get_trans(self):
		"""test get trans ................................."""
		for scale in [1.0,2.0]:
			for mirror in [True,False]:
				dx = Util.get_frand(-30.0,30.0)
				dy = Util.get_frand(-30.0,30.0)
				dz = Util.get_frand(-30.0,30.0)
				d = {"tx":dx,"mirror":mirror,"scale":scale}
				t = Transform(d)
				v = t.get_trans()
				self.assertAlmostEqual(v[0], dx, 5)
			
				d = {"tx":dx,"ty":dy,"mirror":mirror,"scale":scale}
				t = Transform(d)
				v = t.get_trans()
				self.assertAlmostEqual(v[0], dx, 5)
				self.assertAlmostEqual(v[1], dy, 5)
				
				d = {"tx":dx,"ty":dy,"tz":dz,"mirror":mirror,"scale":scale}
				t = Transform(d)
				v = t.get_trans()
				self.assertAlmostEqual(v[0], dx, 5)
				self.assertAlmostEqual(v[1], dy, 5)
				self.assertAlmostEqual(v[2], dz, 5)
	
	def test_get_trans_2d(self):
		"""test get trans 2d .............................."""
		for scale in [1.0,2.0]:
			for mirror in [True,False]:
				dx = Util.get_frand(-30.0,30.0)
				dy = Util.get_frand(-30.0,30.0)
				dz = Util.get_frand(-30.0,30.0)
				d = {"tx":dx,"mirror":mirror,"scale":scale,"type":"2d"}
				t = Transform(d)
				v = t.get_trans_2d()
				self.assertAlmostEqual(v[0], dx, 5)
			
				d = {"tx":dx,"ty":dy,"mirror":mirror,"scale":scale,"type":"2d"}
				t = Transform(d)
				v = t.get_trans_2d()
				self.assertAlmostEqual(v[0], dx, 5)
				self.assertAlmostEqual(v[1], dy, 5)
				
				# should even work in this case, dz is ignored
				#d = {"tx":dx,"ty":dy,"tz":dz,"mirror":mirror,"scale":scale,"type":"2d"}
				#t = Transform(d)
				#v = t.get_trans_2d()
				#self.assertAlmostEqual(v[0], dx, 5)
				#self.assertAlmostEqual(v[1], dy, 5)
	
	def test_get_rotation(self):
		"""test get rotation .............................."""	
		
		
		for scale in [1.0,2.0]:
			for mirror in [False,True]:
				#2D convention
				d = {"type":"2d","alpha":self.get_angle_rand(),"mirror":mirror,"scale":scale}
				t = Transform(d)
				rot = t.get_rotation("2d")
				self.assertEqual("2d",rot["type"])
		
				self.assertAlmostEqual(d["alpha"]%360.0, rot["alpha"]%360.0, 3)
				
				# EMAN convention
				d = {"type":"eman","az":self.get_angle_rand(),"alt":self.get_angle_rand(0,179),"phi":self.get_angle_rand(),"mirror":mirror,"scale":scale}
				t = Transform(d)
				rot = t.get_rotation("eman")
				self.assertEqual("eman",rot["type"])
		
				self.assertAlmostEqual(d["az"]%360.0, rot["az"]%360.0, 3)
				self.assertAlmostEqual(d["alt"]%180.0, rot["alt"]%360.0, 3)
				self.assertAlmostEqual(d["phi"]%360.0, rot["phi"]%360.0, 3)
				
				# SPIDER convention
				d = {"type":"spider","phi":self.get_angle_rand(),"theta":self.get_angle_rand(0,179),"psi":self.get_angle_rand(),"mirror":mirror,"scale":scale}
				t = Transform(d)
				rot = t.get_rotation("spider")
				self.assertEqual("spider",rot["type"])
		
				self.assertAlmostEqual(d["phi"]%360.0, rot["phi"]%360.0, 3)
				self.assertAlmostEqual(d["theta"]%360.0, rot["theta"]%360.0, 3)
				self.assertAlmostEqual(d["psi"]%360.0, rot["psi"]%360.0, 3)
			
				# IMAGIC convention
				d = {"type":"imagic","alpha":self.get_angle_rand(),"beta":self.get_angle_rand(0,179),"gamma":self.get_angle_rand(),"mirror":mirror,"scale":scale}
				t = Transform(d)
				rot = t.get_rotation("imagic")
				self.assertEqual("imagic",rot["type"])
		
				self.assertAlmostEqual(d["alpha"]%360.0, rot["alpha"]%360.0, 3)
				self.assertAlmostEqual(d["beta"]%360.0, rot["beta"]%360.0, 3)
				self.assertAlmostEqual(d["gamma"]%360.0, rot["gamma"]%360.0, 3)
				
				# MRC convention
				d = {"type":"mrc","phi":self.get_angle_rand(),"theta":self.get_angle_rand(0,179),"omega":self.get_angle_rand(),"mirror":mirror,"scale":scale}
				t = Transform(d)
				rot = t.get_rotation("mrc")
				self.assertEqual("mrc",rot["type"])
		
				self.assertAlmostEqual(d["phi"]%360.0, rot["phi"]%360.0, 3)
				self.assertAlmostEqual(d["theta"]%360.0, rot["theta"]%360.0, 3)
				self.assertAlmostEqual(d["omega"]%360.0, rot["omega"]%360.0, 3)
				
		# XTILT convention
		d = {"type":"xyz","xtilt":self.get_angle_rand(),"ytilt":self.get_angle_rand(),"ztilt":self.get_angle_rand()}
		t = Transform(d)
		rot = t.get_rotation("xyz")
		self.assertEqual("xyz",rot["type"])
		t1 = Transform(rot)
		self.assert_matrix_equality(t,t1)
		
		# SPIN convention
		for e0 in [-10,10]:
			# solving a rotation matrix to deduce a quaternion style rotation
			# has more than one solution.
			n = Vec3f(1,-1,-.5)
			norm = n.normalize()
			d = {"type":"spin","Omega":e0,"n1":n[0],"n2":n[1],"n3":n[2]}
			t = Transform(d)
			rot = t.get_rotation("spin")
			self.assertEqual("spin",rot["type"])
			t1 = Transform(rot)
			
			# check to make sure the rotation  matrix has exactly the same form
			self.assert_matrix_equality(t,t1)
		
		# QUATERNION convention
		for alpha in [-10,10]:
			# solving a rotation matrix to deduce a quaternion style rotation
			# has more than one solution.
			
			e0 = math.cos(alpha*math.pi/180.0)
			n = Vec3f(1,-1,-.5)
			norm = n.normalize()
			sin_alpha = math.sin(alpha*math.pi/180.0)
			e1 = sin_alpha*n[0]
			e2 = sin_alpha*n[1]
			e3 = sin_alpha*n[2]
			d = {"type":"quaternion","e0":e0,"e1":e1,"e2":e2,"e3":e3}
			t = Transform(d)
			rot = t.get_rotation("quaternion")
			self.assertEqual("quaternion",rot["type"])
			t1 = Transform(rot)
			
			# check to make sure the rotation  matrix has exactly the same form
			self.assert_matrix_equality(t,t1)
					
		# SGIROT convention
		for e0 in [-10,10]:
			# solving a rotation matrix to deduce a quaternion style rotation
			# has more than one solution.
			n = Vec3f(1,-1,-.5)
			norm = n.normalize()
			d = {"type":"sgirot","q":e0,"n1":n[0],"n2":n[1],"n3":n[2]}
			t = Transform(d)
			rot = t.get_rotation("sgirot")
			self.assertEqual("sgirot",rot["type"])
			t1 = Transform(rot)
			
			# check to make sure the rotation  matrix has exactly the same form
			self.assert_matrix_equality(t,t1)
			
		#MATRIX convention
		#note this is bad because we are not constructing a matrix
		#that reflects a true rotation
		d = {"type":"eman","az":3,"alt":5,"phi":-1}
		t = Transform(d)
		m11, m12, m13 = t.at(0,0),t.at(0,1),t.at(0,2)
		m21, m22, m23 = t.at(1,0),t.at(1,1),t.at(1,2)
		m31, m32, m33 = t.at(2,0),t.at(2,1),t.at(2,2)
		d = {"type":"matrix","m11":m11,"m12":m12,"m13":m13,"m21":m21,"m22":m22,"m23":m23,"m31":m31,"m32":m32,"m33":m33}
		t = Transform(d)
		rot = t.get_rotation("matrix")
		self.assertEqual("matrix",rot["type"])
		self.assertAlmostEqual(d["m11"],rot["m11"], 3)
		self.assertAlmostEqual(m12,rot["m12"], 3)
		self.assertAlmostEqual(m13,rot["m13"], 3)
		self.assertAlmostEqual(m21,rot["m21"], 3)
		self.assertAlmostEqual(m22,rot["m22"], 3)
		self.assertAlmostEqual(m23,rot["m23"], 3)
		self.assertAlmostEqual(m31,rot["m31"], 3)
		self.assertAlmostEqual(m32,rot["m32"], 3)
		self.assertAlmostEqual(m33,rot["m33"], 3)
		
	def test_set_get_scale(self):
		"""test set/get scale ............................."""
		for scale in [1.0,2.0]:
			for mirror in [True, False]:
				d = {"type":"2d","alpha":self.get_angle_rand(),"mirror":mirror,"scale":scale}
				t = Transform(d)
				self.assertAlmostEqual(scale,t.get_scale(), 5)
				
				d = {"type":"eman","az":self.get_angle_rand(),"alt":self.get_angle_rand(0,179),"phi":self.get_angle_rand(),"mirror":mirror,"scale":scale}
				t = Transform(d)
				self.assertAlmostEqual(scale,t.get_scale(), 5)
				
				d = {"mirror":mirror,"scale":scale}
				t = Transform(d)
				self.assertAlmostEqual(scale,t.get_scale(), 5)
					
	def test_set_get_mirror(self):
		"""test set/get mirror............................."""
		for scale in [1.0,2.0]:
			for mirror in [True, False]:
				d = {"type":"2d","alpha":self.get_angle_rand(),"mirror":mirror,"scale":scale}
				t = Transform(d)
				self.assertEqual(mirror,t.get_mirror())
				
				d = {"type":"eman","az":self.get_angle_rand(),"alt":self.get_angle_rand(0,179),"phi":self.get_angle_rand(),"mirror":mirror,"scale":scale}
				t = Transform(d)
				self.assertEqual(mirror,t.get_mirror())
				
				d = {"mirror":mirror,"scale":scale}
				t = Transform(d)
				self.assertEqual(mirror,t.get_mirror())
	def test_get_set_params(self):
		"""test set/get params............................."""
		t = Transform()
		t.set_params({"type":"eman","az":10,"alt":150,"scale":2.0,"mirror":True,"tx":3.4})
		d = t.get_params("eman")
		s = Transform(d)
		self.assert_matrix_equality(s,t)
		d = t.get_params("spider")
		s = Transform(d)
		self.assert_matrix_equality(s,t)
		d = t.get_params("mrc")
		s = Transform(d)
		self.assert_matrix_equality(s,t)
		d = t.get_params("imagic")
		s = Transform(d)
		self.assert_matrix_equality(s,t)
		d = t.get_params("quaternion")
		s = Transform(d)
		self.assert_matrix_equality(s,t)
		d = t.get_params("sgirot")
		s = Transform(d)
		self.assert_matrix_equality(s,t)
		d = t.get_params("spin")
		s = Transform(d)
		self.assert_matrix_equality(s,t)
		d = t.get_params("matrix")
		s = Transform(d)
		self.assert_matrix_equality(s,t)
		d = t.get_params("xyz")
		s = Transform(d)
		self.assert_matrix_equality(s,t)
	def test_get_set_params_2d(self):
		"""test set/get params 2d.........................."""
		t = Transform()
		t.set_params({"type":"2d","alpha":10,"scale":2.0,"mirror":True,"tx":3.4,"ty":0.0})
		d = t.get_params("2d") # no euler type required because there is only one ("2d")
		s = Transform(d) # s is the same as t
		self.assert_matrix_equality(s,t)
	def test_multiplication(self):
		"""test multiplication............................."""
		for v in [Vec3f(1,1,1),[1,1,1]]:
			t = Transform()
			v = Vec3f(1,1,1)
			v_d = t*v
			self.assertAlmostEqual(v[0],v_d[0], 5)
			self.assertAlmostEqual(v[1],v_d[1], 5)
			self.assertAlmostEqual(v[2],v_d[2], 5)
			
			t = Transform()
			scale = 2.3
			t.set_scale(scale)
			v = Vec3f(1,1,1)
			v_d = t*v
			self.assertAlmostEqual(scale*v[0],v_d[0], 5)
			self.assertAlmostEqual(scale*v[1],v_d[1], 5)
			self.assertAlmostEqual(scale*v[2],v_d[2], 5)
			
			t = Transform()
			t.set_mirror(True)
			v = Vec3f(1,1,1)
			v_d = t*v
			self.assertAlmostEqual(-v[0],v_d[0], 5)
			self.assertAlmostEqual(v[1],v_d[1], 5)
			self.assertAlmostEqual(v[2],v_d[2], 5)
			
			t = Transform()
			dx = 2
			dy = -1
			dz = .23232234
			t.set_trans(dx,dy,dz)
			v = Vec3f(1,1,1)
			v_d = t*v
			self.assertAlmostEqual(v[0]+dx,v_d[0], 5)
			self.assertAlmostEqual(v[1]+dy,v_d[1], 5)
			self.assertAlmostEqual(v[2]+dz,v_d[2], 5)
	
	def test_transform(self):
		"""test transform.................................."""
		t = Transform()
		v = Vec3f(1,1,1)
		v_d = t.transform(v)
		self.assertAlmostEqual(v[0],v_d[0], 5)
		self.assertAlmostEqual(v[1],v_d[1], 5)
		self.assertAlmostEqual(v[2],v_d[2], 5)
		
		t = Transform()
		scale = 2.3
		t.set_scale(scale)
		v = Vec3f(1,1,1)
		v_d = t.transform(v)
		self.assertAlmostEqual(scale*v[0],v_d[0], 5)
		self.assertAlmostEqual(scale*v[1],v_d[1], 5)
		self.assertAlmostEqual(scale*v[2],v_d[2], 5)
		
		t = Transform()
		t.set_mirror(True)
		v = Vec3f(1,1,1)
		v_d = t.transform(v)
		self.assertAlmostEqual(-v[0],v_d[0], 5)
		self.assertAlmostEqual(v[1],v_d[1], 5)
		self.assertAlmostEqual(v[2],v_d[2], 5)
		
		t = Transform()
		dx = 2
		dy = -1
		dz = .23232234
		t.set_trans(dx,dy,dz)
		v = Vec3f(1,1,1)
		v_d = t.transform(v)
		self.assertAlmostEqual(v[0]+dx,v_d[0], 5)
		self.assertAlmostEqual(v[1]+dy,v_d[1], 5)
		self.assertAlmostEqual(v[2]+dz,v_d[2], 5)
		
	def test_multiplication_2d(self):
		"""test multiplication 2d.........................."""
		t = Transform()
		v = Vec2f(1,1)
		v_d = t*v
		self.assertAlmostEqual(v[0],v_d[0], 5)
		self.assertAlmostEqual(v[1],v_d[1], 5)
		
		t = Transform()
		scale = 2.3
		t.set_scale(scale)
		v = Vec2f(1,1)
		v_d = t*v
		self.assertAlmostEqual(scale*v[0],v_d[0], 5)
		self.assertAlmostEqual(scale*v[1],v_d[1], 5)
		
		t = Transform()
		t.set_mirror(True)
		v = Vec2f(1,1)
		v_d = t*v
		self.assertAlmostEqual(-v[0],v_d[0], 5)
		self.assertAlmostEqual(v[1],v_d[1], 5)
		
		t = Transform()
		dx = 2
		dy = -1.0032023
		t.set_trans(dx,dy)
		v = Vec2f(1,1)
		v_d = t*v
		self.assertAlmostEqual(v[0]+dx,v_d[0], 5)
		self.assertAlmostEqual(v[1]+dy,v_d[1], 5)
	
	def test_transform_2d(self):
		"""test transform 2d..............................."""
		t = Transform()
		v = Vec2f(1,1)
		v_d = t.transform(v)
		self.assertAlmostEqual(v[0],v_d[0], 5)
		self.assertAlmostEqual(v[1],v_d[1], 5)
		
		t = Transform()
		scale = 2.3
		t.set_scale(scale)
		v = Vec2f(1,1)
		v_d = t.transform(v)
		self.assertAlmostEqual(scale*v[0],v_d[0], 5)
		self.assertAlmostEqual(scale*v[1],v_d[1], 5)
		
		t = Transform()
		t.set_mirror(True)
		v = Vec2f(1,1)
		v_d = t.transform(v)
		self.assertAlmostEqual(-v[0],v_d[0], 5)
		self.assertAlmostEqual(v[1],v_d[1], 5)
		
		t = Transform()
		dx = 2
		dy = -1.0032023
		t.set_trans(dx,dy)
		v = Vec2f(1,1)
		v_d = t.transform(v)
		self.assertAlmostEqual(v[0]+dx,v_d[0], 5)
		self.assertAlmostEqual(v[1]+dy,v_d[1], 5)
		
#	def test_set_get_post_x_mirror(self):
#		"""test set/get post_x_mirror ....................."""
#		for t in TestTransform.transforms:
#			t.set_post_x_mirror(False)
#			self.assertEqual(False,t.get_post_x_mirror())
#		
#		for t in TestTransform.transforms:
#			t.set_post_x_mirror(True)
#			self.assertEqual(True,t.get_post_x_mirror())
#	
#	def test_set_get_trans(self):
#		"""test set/get trans ........................."""
#		x = 1.3
#		y = 3.3
#		z = -9.0
#		post_trans = Vec3f(x,y,z)
#		for t in TestTransform.transforms:
#			t.set_trans(post_trans)
#			v = t.get_trans()
#			self.assertAlmostEqual(x,v[0], 5)
#			self.assertAlmostEqual(y,v[1], 5)
#			self.assertAlmostEqual(z,v[2], 5)
#			
	def test_inverse_invert(self):
		"""test inverse/invert ............................"""
		no_trans = {}
		two_trans = {"tx":1.023,"ty":-1.002}
		three_trans = {"tx":.023,"ty":431.22002,"tz":120.02}
		for scale in [1.0,2.0]:
			for mirror in [True, False]:
				for i,trans in enumerate([no_trans,two_trans,three_trans]):
					if i < 2:
						d = {"type":"2d","alpha":self.get_angle_rand(),"mirror":mirror,"scale":scale}
						t = Transform(d)
						trans1 = {}
						for c in trans: trans1[c] = trans[c]
						trans1["type"] = "2d"
						t.set_params(trans1)
						s = t.inverse()
						self.assert_identity(s*t)
						self.assert_identity(t*s)
						s = Transform(t)
						s.invert()
						self.assert_identity(s*t)
						self.assert_identity(t*s)
					
					d = {"type":"eman","az":self.get_angle_rand(),"alt":self.get_angle_rand(0,179),"phi":self.get_angle_rand(),"mirror":mirror,"scale":scale}
					t = Transform(d)
					t.set_params(trans)
					s = t.inverse()
					self.assert_identity(s*t)
					self.assert_identity(t*s)
					s = Transform(t)
					s.invert()
					self.assert_identity(s*t)
					self.assert_identity(t*s)
		
					d = {"mirror":mirror,"scale":scale}
					t = Transform(d)
					t.set_params(trans)
					s = t.inverse()
					self.assert_identity(s*t)
					self.assert_identity(t*s)
					s = Transform(t)
					s.invert()
					self.assert_identity(s*t)
					self.assert_identity(t*s)
	def test_copy_construction(self):
		"""test copy construction.........................."""
		three_trans = {"tx":.023,"ty":431.22002,"tz":120.02}
		d = {"type":"eman","az":self.get_angle_rand(),"alt":self.get_angle_rand(0,179),"phi":self.get_angle_rand()}
		t = Transform(d)
		t.set_params(three_trans) # I.E. now we have all matrix elements filled
		t1 = Transform(t)
		self.assert_matrix_equality(t,t1)
#	
	def test_transpose(self):
		"""test transpose and transpose inplace............"""
		
		for mirror in [True, False]:
			d = {"type":"2d","alpha":self.get_angle_rand(),"mirror":mirror}
			t = Transform(d)
			t1 = t.transpose()
			self.assert_identity(t1*t,3)
			self.assert_identity(t*t1,3)
			
			d = {"type":"eman","az":self.get_angle_rand(),"alt":self.get_angle_rand(0,179),"phi":self.get_angle_rand()}
			t = Transform(d)
			t = Transform(d)
			t1 = t.transpose()
			self.assert_identity(t1*t,3)
			self.assert_identity(t*t1,3)

			t = Transform()
			t = Transform(d)
			t1 = t.transpose()
			self.assert_identity(t1*t,3)
			self.assert_identity(t*t1,3)
		
	def test_get_pre_trans(self):
		"""test get pre trans.............................."""
		dx = .023
		dy = 431.220002
		dz = 120.02
		three_trans = {"tx":dx,"ty":dy,"tz":dz}
		t = Transform(three_trans)
		v = t.get_pre_trans()
		
		self.assertAlmostEqual(v[0],dx, 5)
		self.assertAlmostEqual(v[1],dy, 5)
		self.assertAlmostEqual(v[2],dz, 5)
		
		scale = Util.get_frand(1.00001,100.0)
		t.set_scale(scale)
		v = t.get_pre_trans()
		
		self.assertAlmostEqual(v[0]*scale,dx, 3)
		self.assertAlmostEqual(v[1]*scale,dy, 3)
		self.assertAlmostEqual(v[2]*scale,dz, 3)
		
		t.set_mirror(True)
		v = t.get_pre_trans()
		self.assertAlmostEqual(v[0]*scale,dx, 3)
		self.assertAlmostEqual(v[1]*scale,dy, 3)
		self.assertAlmostEqual(v[2]*scale,dz, 3)
		
		# finally perhaps the most important test
		d = {"type":"eman","az":self.get_angle_rand(),"alt":self.get_angle_rand(0,179),"phi":self.get_angle_rand()}
		t.set_params(d)
		v = t.get_pre_trans()
		pre_trans = Transform()
		pre_trans.set_trans(v)
		
		rot = t.get_rotation("eman")
		scale = t.get_scale()
		mirror = t.get_mirror()
		
		without_trans = Transform(rot)
		without_trans.set_scale(scale)
		without_trans.set_mirror(mirror)
		
		#without_trans.printme()
		#pre_trans.printme()
		#(without_trans*pre_trans).printme()
		#(pre_trans*without_trans).printme()
		self.assert_matrix_equality(without_trans*pre_trans,t)
	
	def assert_identity(self,t2,n=4):
		for j in range(n):
			for i in range(n):
				if i == j:
					self.assertAlmostEqual(t2.at(i,j),1, 3)
				else:
					self.assertAlmostEqual(t2.at(i,j),0, 3)
	
	def assert_matrix_equality(self,t,t1):
		for j in range(4):
			for i in range(4):
				if t1.at(i,j) != 0:
					self.assertAlmostEqual(t.at(i,j),t1.at(i,j), 3)
				else:
					self.assertAlmostEqual(t.at(i,j),t1.at(i,j), 3)

class TestSymmetry3D(unittest.TestCase):
	def assert_reduction_works(self,i,az,alt,azmax,sym):
		T = Transform({"type":"eman","az":az,"alt":alt,"phi":0})
		T1 = sym.get_sym(i)
		T2 = T*T1
		A = sym.reduce(T2,0)
		#print i,i
		result = A.get_rotation("eman")
		azsoln = result["az"] %azmax
		if ( azsoln < azmax and (azmax - azsoln) < 0.001 ): azsoln = 0.0
		#print i,i,i
		#print result["az"],az,result["az"] % (720.0/n),az% (720.0/n),720.0/n
		
		if (az%azmax) == 0 or azsoln == 0:
			self.assertAlmostEqual(azsoln,(az%azmax), 3)
		else:
			self.assertAlmostEqual(azsoln/(az%azmax),1.0, 3)
		if alt == 0 or result["alt"] == 0:
			self.assertAlmostEqual(result["alt"],alt, 3)
		else:
			self.assertAlmostEqual(result["alt"]/alt,1.0, 3)
	
	def test_symc_reduce(self):
		"""test csym reduce ..............................."""
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
			eulers = sym.gen_orientations("eman",{"delta":12})
			for euler in eulers:
				rot = euler.get_rotation("eman")
				az = rot["az"]
				alt = rot["alt"]
				#print az,alt,n
				for i in range(1,n):
					self.assert_reduction_works(i,az,alt,azmax,sym)
	def test_symd_reduce(self):
		"""test dsym reduce ..............................."""
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
			
			eulers = sym.gen_orientations("eman",{"delta":12})
			for euler in eulers:
				rot = euler.get_rotation("eman")
				az = rot["az"]
				alt = rot["alt"]
				#print az,alt,n
				for i in range(1,n):
					self.assert_reduction_works(i,az,alt,azmax,sym)
	def test_symtet_reduce(self):
		"""test tetsym reduce ............................."""
		syms = []
		syms.append(Symmetries.get("tet",{}))
		for sym in syms:
			n = sym.get_nsym()
			azmax = 120.0
			
			eulers = sym.gen_orientations("eman",{"delta":12})
			for euler in eulers:
				rot = euler.get_rotation("eman")
				az = rot["az"]
				alt = rot["alt"]
				#print az,alt,n
				for i in range(1,n):
					self.assert_reduction_works(i,az,alt,azmax,sym)
		
	def test_symoct_reduce(self):
		"""test octsym reduce ............................."""
		syms = []
		syms.append(Symmetries.get("oct",{}))
		for sym in syms:
			n = sym.get_nsym()
			azmax = 90.0
			
			eulers = sym.gen_orientations("eman",{"delta":12})
			for euler in eulers:
				rot = euler.get_rotation("eman")
				az = rot["az"]
				alt = rot["alt"]
				#print az,alt,n
				for i in range(1,n):
					self.assert_reduction_works(i,az,alt,azmax,sym)

	def test_symicos_reduce(self):
		"""test icossym reduce ............................"""
		syms = []
			
		syms.append(Symmetries.get("icos",{}))
		for sym in syms:
			n = sym.get_nsym()
			azmax = 72.0
			
			eulers = sym.gen_orientations("eman",{"delta":12})
			for euler in eulers:
				rot = euler.get_rotation("eman")
				az = rot["az"]
				alt = rot["alt"]
				#print az,alt,n
				for i in range(1,n):
					self.assert_reduction_works(i,az,alt,azmax,sym)

class TestTransform3D(unittest.TestCase):
	"""this is the unit test for Transform3D class"""
	
	#def test_get_set_params(self):
		#"""test get/set params .............................."""
		#az = 
		#t = Transform3D(az, alt, phi)
		#rot = t.get_rotation(EULER_EMAN)

		##testlib.assertfloat(self, az, float(rot["az"]))
		#self.assertAlmostEqual(az, rot["az"], 3)
		#testlib.assertfloat(self, alt, float(rot["alt"]))
		#testlib.assertfloat(self, phi, float(rot["phi"]))
		
	def test_get_rotation(self):
		"""test get rotation ................................"""
		alt = 1.45232928554
		az = -0.60170830102
		phi = 10.0
	
		# EMAN convention
		t = Transform3D(az, alt, phi)
		rot = t.get_rotation(EULER_EMAN)

		self.assertAlmostEqual(az%360.0, rot["az"]%360.0, 3)
		self.assertAlmostEqual(alt%360.0, rot["alt"]%360.0, 3)
		self.assertAlmostEqual(phi%360.0, rot["phi"]%360.0, 3)
		
		# SPIDER convention
		phi = az
		theta = alt
		psi = phi
		t = Transform3D(EULER_SPIDER,phi,theta,psi)
		rot = t.get_rotation(EULER_SPIDER)

		self.assertAlmostEqual(phi%360.0, rot["phi"]%360.0, 3)
		self.assertAlmostEqual(theta%360.0, rot["theta"]%360.0, 3)
		self.assertAlmostEqual(psi%360.0, rot["psi"]%360.0, 3)
	
		# IMAGIC convention
		alpha = az
		beta = alt
		gamma = phi
		t = Transform3D(EULER_IMAGIC,alpha,beta,gamma)
		rot = t.get_rotation(EULER_IMAGIC)

		self.assertAlmostEqual(alpha%360.0, rot["alpha"]%360.0, 3)
		self.assertAlmostEqual(beta%360.0, rot["beta"]%360.0, 3)
		self.assertAlmostEqual(gamma%360.0, rot["gamma"]%360.0, 3)
		
		# MRC convention
		phi = az
		theta = alt
		omega = phi
		t = Transform3D(EULER_MRC,phi,theta,omega)
		rot = t.get_rotation(EULER_MRC)

		self.assertAlmostEqual(phi%360.0, rot["phi"]%360.0, 3)
		self.assertAlmostEqual(theta%360.0, rot["theta"]%360.0, 3)
		self.assertAlmostEqual(omega%360.0, rot["omega"]%360.0, 3)
		
		# XTILT convention
		xtilt = az
		ytilt = alt
		ztilt = phi
		t = Transform3D(EULER_XYZ,xtilt,ytilt,ztilt)
		rot = t.get_rotation(EULER_XYZ)

		self.assertAlmostEqual(xtilt%360.0, rot["xtilt"]%360.0, 3)
		self.assertAlmostEqual(ytilt%360.0, rot["ytilt"]%360.0, 3)
		self.assertAlmostEqual(ztilt%360.0, rot["ztilt"]%360.0, 3)
	
		# SPIN convention
		for e0 in [-10,10]:
			# solving a rotation matrix to deduce a quaternion style rotation
			# has more than one solution.
			n = Vec3f(1,-1,-.5)
			norm = n.normalize()
			t = Transform3D(EULER_SPIN,e0,n[0],n[1],n[2])
			rot = t.get_rotation(EULER_SPIN)
			t1 = Transform3D(EULER_SPIN,rot["Omega"],rot["n1"],rot["n2"],rot["n3"])
			
			# check to make sure the rotation  matrix has exactly the same form
			for j in range(0,2):
				for i in range(0,2):
					self.assertAlmostEqual(t.at(i,j), t1.at(i,j), 3)
		
		# QUATERNION convention
		for alpha in [-10,10]:
			# solving a rotation matrix to deduce a quaternion style rotation
			# has more than one solution.
			
			e0 = math.cos(alpha*math.pi/180.0)
			n = Vec3f(1,-1,-.5)
			norm = n.normalize()
			sin_alpha = math.sin(alpha*math.pi/180.0)
			e1 = sin_alpha*n[0]
			e2 = sin_alpha*n[1]
			e3 = sin_alpha*n[2]
			t = Transform3D(EULER_QUATERNION,e0,e1,e2,e3)
			rot = t.get_rotation(EULER_QUATERNION)
			t1 = Transform3D(EULER_QUATERNION,rot["e0"],rot["e1"],rot["e2"],rot["e3"])
			
			# check to make sure the rotation  matrix has exactly the same form
			for j in range(0,2):
				for i in range(0,2):
					self.assertAlmostEqual(t.at(i,j), t1.at(i,j), 3)
					
		# SGIROT convention
		for e0 in [-10,10]:
			# solving a rotation matrix to deduce a quaternion style rotation
			# has more than one solution.
			n = Vec3f(1,-1,-.5)
			norm = n.normalize()
			t = Transform3D(EULER_SGIROT,e0,n[0],n[1],n[2])
			rot = t.get_rotation(EULER_SGIROT)
			t1 = Transform3D(EULER_SGIROT,rot["q"],rot["n1"],rot["n2"],rot["n3"])
			
			# check to make sure the rotation  matrix has exactly the same form
			for j in range(0,2):
				for i in range(0,2):
					self.assertAlmostEqual(t.at(i,j), t1.at(i,j), 3)	
	
		# retrospective note = we don't support retrieving EULER_MATRIX style
		# rotations from python even though it's possible in C++
		# MATRIX convention
		# note this is bad because we are not constructing a matrix
		# that reflects a true rotation
		#m11, m12, m13 = 1,2,3
		#m21, m22, m23 = 4,5,6
		#m31, m32, m33 = 7,8,9
		#t = Transform3D(m11,m12,m13,m21,m22,m23,m31,m32,m33)
		#rot = t.get_rotation(EULER_MATRIX)
		#self.assertAlmostEqual(m11,rot["m11"], 3)
		#self.assertAlmostEqual(m12,rot["m12"], 3)
		#self.assertAlmostEqual(m13,rot["m13"], 3)
		#self.assertAlmostEqual(m21,rot["m21"], 3)
		#self.assertAlmostEqual(m22,rot["m22"], 3)
		#self.assertAlmostEqual(m23,rot["m23"], 3)
		#self.assertAlmostEqual(m31,rot["m31"], 3)
		#self.assertAlmostEqual(m32,rot["m32"], 3)
		#self.assertAlmostEqual(m33,rot["m33"], 3)	
		
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

	def no_test_angles2tfvec(self):    #no callable in Python
		"""test angles2tfvec() function ....................."""
		t = Transform3D.angles2tfvec(Transform3D.EulerType.EMAN, (1.45232928554, -0.60170830102, 0))	
		
def test_main():
	print 'Testing Transform'
	test_support.run_unittest(TestTransform)
	
	print '\nTesting Symmetry3D'
	test_support.run_unittest(TestSymmetry3D)
	
	print "\nTesting Transform3D"
	test_support.run_unittest(TestTransform3D)


if __name__ == '__main__':
	test_main()
