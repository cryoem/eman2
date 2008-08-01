#!/usr/bin/env python

#
# Author: Grant Tang, 09/01/2005 (gtang@bcm.edu)
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
from EMAN2db import *
import unittest,os,sys
from test import test_support
import testlib
from pyemtbx.exceptions import *
from math import *

class TestAligner(unittest.TestCase):
	"""aligner test"""
	
	def test_TranslationalAligner(self):
		"""test TranslationalAligner ........................"""
		e = EMData()
		e.set_size(32,32,1)
		e.process_inplace('testimage.noise.uniform.rand')
		
		e2 = EMData()
		e2.set_size(32,32,1)
		e2.process_inplace('testimage.noise.uniform.rand')
		
		e.align('translational', e2, {})
		e.align('translational', e2, {'intonly':1, 'maxshift':2})
		
		# Test 2D behavior
		for z in (32,33):
			for y in (32,33):
				for x in (32,33): # does not work yet for odd x - rotational footprint fails
					size = (x,y,z)
					ref = test_image_3d(0,size)
					for i in range(0,20):
						e =test_image_3d(0,size)
						dx = Util.get_frand(-3,3)
						dy = Util.get_frand(-3,3)
						dz = Util.get_frand(-3,3)
						e.translate(dx,dy,dz)
						g = e.align("translational",ref,{},"dot",{})
						self.failIf(fabs(g.get_attr("align.dx") + dx) > 1)
						self.failIf(fabs(g.get_attr("align.dy") + dy) > 1)
						self.failIf(fabs(g.get_attr("align.dz") + dz) > 1)
						
								
		# Test 2D behavior
		for y in (32,33):
			for x in (32,33): # does not work yet for odd x - rotational footprint fails
				size = (x,y)
				images = []
				ref = test_image(0,size)
				for i in range(0,20):
					e =test_image(0,size)
					dx = Util.get_frand(-3,3)
					dy = Util.get_frand(-3,3)
					e.translate(dx,dy,0)
					g = e.align("translational",ref,{},"dot",{})
					self.failIf(fabs(g.get_attr("align.dx") + dx) > 1)
					self.failIf(fabs(g.get_attr("align.dy") + dy) > 1)

				
		# Test 1D behavior
		#for i in [n-1]: #FIXME - do for i in [n-1,n] for some reason n even fails
			#e = EMData()
			#e.set_size(i,1,1);
			#e.to_zero()
			#e.process_inplace("testimage.x")
			
			#for dx in [-1,0,1]:
				#f = e.copy()
				#f.translate(dx,0,0)
				
				#g = f.align('translational', e, {'maxshift':3})

				#sdx = g.get_attr("align.dx")
					
				##print "%d %d" %(dx,sdx)
				
				#self.assertEqual(-sdx,dx)

	def test_RotationalAligner(self):
		"""test RotationalAligner ..........................."""
		e = EMData()
		e.set_size(32,32,1)
		e.process_inplace('testimage.noise.uniform.rand')
		
		e2 = EMData()
		e2.set_size(32,32,1)
		e2.process_inplace('testimage.noise.uniform.rand')
		
		e.align('rotational', e2)
		
		for y in [32]:
			for x in [32]: # does not work yet for odd x - rotational footprint fails
				size = (x,y)
				ref = test_image(0,size)
				for i in range(0,20):
					e =test_image(0,size)
					az = Util.get_frand(0,360)
					e.rotate(az,0,0)
					g = e.align("rotational",ref,{},"dot",{})
					
					result = fabs(g.get_attr("align.az") + az)
					#print g.get_attr("align.az"), az
					if result > 180 and result < 360:
						result = 360-result
					if result > 360: result = result-360
					self.failIf( result > 3 ) # 3 seems accurate enough
				
	def no_test_RotatePrecenterAligner(self):
		"""test RotatePrecenterAligner ......................"""
		e = EMData()
		#e.set_size(32,32,1)
		#e.process_inplace('testimage.noise.uniform.rand')
		
		e2 = EMData()
		e2.set_size(32,32,1)
		e2.process_inplace('testimage.noise.uniform.rand')
		
		e.align('rotate_precenter', e2)
		
	def test_RotateCHAligner(self):
		"""test RotateCHAligner ............................."""
		e = EMData()
		e.set_size(32,32,1)
		e.process_inplace('testimage.noise.uniform.rand')
		
		e2 = EMData()
		e2.set_size(32,32,1)
		e2.process_inplace('testimage.noise.uniform.rand')
		
		e.align('rotate_ch', e2, {'irad':1, 'orad':2})
		
	def test_RotateTranslateAligner(self):
		"""test RotateTranslateAligner ......................"""
		e = EMData()
		e.set_size(32,32,1)
		e.process_inplace('testimage.noise.uniform.rand')
		
		e2 = EMData()
		e2.set_size(32,32,1)
		e2.process_inplace('testimage.noise.uniform.rand')
		
		e.align('rotate_translate', e2, {'maxshift':1})
		
		for y in [32]:
			for x in [32]:# does not work yet for odd x - rotational footprint fails
				size = (x,y)
				ref = test_image(7,size)
				for i in range(0,20):
					e = ref.copy()
					dx = Util.get_frand(-3,3)
					dy = Util.get_frand(-3,3)
					az = Util.get_frand(0,360)
					t = Transform3D(az,0,0)
					t.set_pretrans(dx,dy,0)
					e.rotate_translate(t)
					g = e.align("rotate_translate",ref,{},"dot")
					t1 = Transform3D(g.get_attr("align.az"),0,0)
					t1.set_posttrans(g.get_attr("align.dx"),g.get_attr("align.dy"),0)
					v = Vec3f(0,1,0)
					vd = v*t
					vd = vd*t1
					#print vd[0],vd[1],vd[2],az,g.get_attr("align.az"),dx,dy
					self.failIf(fabs(vd[0]) > 0.15)
					self.failIf(fabs(vd[1]-1) > 0.15)

		
	def no_test_RotateTranslateBestAligner(self):
		"""test RotateTranslateBestAligner .................."""
		e = EMData()
		e.set_size(32,32,1)
		e.process_inplace('testimage.noise.uniform.rand')
		
		e2 = EMData()
		e2.set_size(32,32,1)
		e2.process_inplace('testimage.noise.uniform.rand')
		
		img = e.align('rotate_translate_best', e2, {'maxshift':1, 'snr':(1.0, 2.0, 3.0)})
		
	def test_RotateTranslateRadonAligner(self):
		"""test RotateTranslateRadonAligner ................."""
		Log.logger().set_level(-1)    #no log message printed out
		e = EMData()
		e.set_size(32,32,1)
		e.process_inplace('testimage.noise.uniform.rand')
		
		e2 = EMData()
		e2.set_size(32,32,1)
		e2.process_inplace('testimage.noise.uniform.rand')
		
		e3 = EMData()
		e3.set_size(32,32,1)
		e3.process_inplace('testimage.noise.uniform.rand')
		
		e4 = EMData()
		e4.set_size(32,32,1)
		e4.process_inplace('testimage.noise.uniform.rand')
		
		img = e.align('rotate_translate_radon', e2, {'maxshift':2, 'radonwith':e3, 'radonthis':e4})
		
		import os
		testlib.safe_unlink('radon.hed')
		testlib.safe_unlink('radon.img')
		testlib.safe_unlink('racf.hed')
		testlib.safe_unlink('racf.img')

	def test_RotateFlipAligner(self):
		"""test RotateFlipAligner ..........................."""
		e = EMData()
		e.set_size(32,32,1)
		e.process_inplace('testimage.noise.uniform.rand')
		
		e2 = EMData()
		e2.set_size(32,32,1)
		e2.process_inplace('testimage.noise.uniform.rand')
		
		e3 = EMData()
		e3.set_size(32,32,1)
		e3.process_inplace('testimage.noise.uniform.rand')
		
		#e.align('rotate_flip', e2, {'flip':e3, 'imask':2})
		e.align('rotate_flip', e2, {'imask':2})
		
	def test_RotateTranslateFlipAligner(self):
		"""test RotateTranslateFlip Aligner ................."""
		e = EMData()
		e.set_size(32,32,1)
		e.process_inplace('testimage.noise.uniform.rand')
		
		e2 = EMData()
		e2.set_size(32,32,1)
		e2.process_inplace('testimage.noise.uniform.rand')
		
		e.align('rotate_translate_flip', e2, {'maxshift':1})
		
		for y in [32]:
			for x in [32]:# does not work yet for odd x - rotational footprint fails
				size = (x,y)
				ref = test_image(7,size)
				for i in range(0,20):
					e = ref.copy()
					flip = Util.get_irand(0,1)
					if flip: e.process_inplace("xform.flip",{"axis":"x"})
					dx = Util.get_frand(-3,3)
					dy = Util.get_frand(-3,3)
					az = Util.get_frand(0,360)
					
					t = Transform3D(az,0,0)
					t.set_pretrans(dx,dy,0)
					e.rotate_translate(t)
					g = e.align("rotate_translate_flip",ref,{"rfp_mode":2},"phase")
					t1 = Transform3D(g.get_attr("align.az"),0,0)
					t1.set_posttrans(g.get_attr("align.dx"),g.get_attr("align.dy"),0)
					
					if g.get_attr("align.flip"):
						t = Transform3D(360-az,0,0)
						t.set_pretrans(dx,dy,0)
					v = Vec3f(0,1,0)
					vd = v*t
					vd = vd*t1
					print vd[0],vd[1],vd[2],az,g.get_attr("align.az"),dx,dy,flip,g.get_attr("align.flip")
					self.failIf(fabs(vd[0]) > 0.2)
					self.failIf(fabs(vd[1]-1) > 0.15)
					self.failIf(flip != g.get_attr("align.flip"))
		
	def no_test_RTFSlowAligner(self):
		"""test RTFSlowAligner .............................."""
		e = EMData()
		e.set_size(32,32,1)
		e.process_inplace('testimage.noise.uniform.rand')
		
		e2 = EMData()
		e2.set_size(32,32,1)
		e2.process_inplace('testimage.noise.uniform.rand')
		
		e3 = EMData()
		e3.set_size(32,32,1)
		e3.process_inplace('testimage.noise.uniform.rand')
		
		e.align('rtf_slow', e2, {'flip':e3, 'maxshift':2})
		
		#RTFSlowestAligner eliminated
	def no_test_RTFSlowestAligner(self):
		"""test RTFSlowestAligner ..........................."""
		e = EMData()
		e.set_size(32,32,1)
		e.process_inplace('testimage.noise.uniform.rand')
		
		e2 = EMData()
		e2.set_size(32,32,1)
		e2.process_inplace('testimage.noise.uniform.rand')
		
		e3 = EMData()
		e3.set_size(32,32,1)
		e3.process_inplace('testimage.noise.uniform.rand')
		
		e.align('rtf_slowest', e2, {'flip':e3, 'maxshift':2})
		
	def no_test_RTFBestAligner(self):
		"""test RTFBestAligner .............................."""
		e = EMData()
		e.set_size(32,32,1)
		e.process_inplace('testimage.noise.uniform.rand')
		
		e2 = EMData()
		e2.set_size(32,32,1)
		e2.process_inplace('testimage.noise.uniform.rand')
		
		e3 = EMData()
		e3.set_size(32,32,1)
		e3.process_inplace('testimage.noise.uniform.rand')
		
		e.align('rtf_best', e2, {'flip':e3, 'maxshift':2, 'snr':(1.0, 2.0, 3.0)})
		
	def no_test_RTFRadonAligner(self):
		"""test RTFRadonAligner ............................."""
		e = EMData()
		e.set_size(32,32,1)
		e.process_inplace('testimage.noise.uniform.rand')
		
		e2 = EMData()
		e2.set_size(32,32,1)
		e2.process_inplace('testimage.noise.uniform.rand')
		
		e3 = EMData()
		e3.set_size(32,32,1)
		e3.process_inplace('testimage.noise.uniform.rand')
		
		e4 = EMData()
		e4.set_size(32,32,1)
		e4.process_inplace('testimage.noise.uniform.rand')
		
		e5 = EMData()
		e5.set_size(32,32,1)
		e5.process_inplace('testimage.noise.uniform.rand')
		
		e6 = EMData()
		e6.set_size(32,32,1)
		e6.process_inplace('testimage.noise.uniform.rand')

		e.align('rtf_radon', e2, {'maxshift':2, 'thisf':e3, 'radonwith':e4, \
				'radonthis':e5, 'radonthisf':e6})
				
	def test_RefineAligner(self):
		"""test RefineAligner ..............................."""
		e = EMData()
		e.set_size(32,32,1)
		e.process_inplace('testimage.noise.uniform.rand')
		
		e2 = EMData()
		e2.set_size(32,32,1)
		e2.process_inplace('testimage.noise.uniform.rand')
		
		e.align('refine', e2)
		e.align('refine', e2, {'mode':1, 'az':1.2, 'dx':2, 'dy':3.4})
		
		
		for y in [32]:
			for x in [32]:# does not work yet for odd x - rotational footprint fails
				size = (x,y)
				ref = test_image(8,size)
				for i in range(0,20):
					e = ref.copy()
					dx = Util.get_frand(-3,3)
					dy = Util.get_frand(-3,3)
					az = Util.get_frand(-5,5)
					t = Transform3D(az,0,0)
					t.set_pretrans(dx,dy,0)
					e.rotate_translate(t)
					g = e.align("refine",ref,{},"dot",{"normalize":1})
					t1 = Transform3D(g.get_attr("align.az"),0,0)
					t1.set_posttrans(g.get_attr("align.dx"),g.get_attr("align.dy"),0)
					v = Vec3f(0,1,0)
					vd = v*t
					vd = vd*t1
					#print vd[0],vd[1],vd[2],az,g.get_attr("align.az"),dx,dy, g.get_attr("align.dx"),g.get_attr("align.dy")
					self.failIf(fabs(vd[0]) > 0.1)
					self.failIf(fabs(vd[1]-1) > 0.02)
		
		#n = 128
		##dx = 0.4
		##dy = 0.9
		#div = 4.0
		##scale=1.0
		#inc = [ i for i in range(-int(div),int(div)+1)]
		#for i in range(-int(div),int(div)+1):
			#inc[i] = inc[i]/div

		#e7 = EMData()
		#e8 = EMData()
		#e7.set_size(32,32,1)
		#e7.process_inplace("testimage.noise.uniform.rand")
		#e8.set_size(32,32,1)
		#e8.process_inplace("testimage.noise.uniform.rand")
		#result2 = e8.cmp("dot",e8, {"normalize":1})
		#result1 = e8.dot_rotate_translate(e8, 0.0, 0.0, 0.0)
		#result3 = e8.dot_rotate_translate(e8, 1.0, 2.0, 3.0)
		#e7.rotate_translate(1.0,0,0,0,0,0,2.0,3.0,0.0)
		#print "Results %f %f %f" %(result1,result2, result3)

		#print ""
		##for i in range(0,2):
		#i = 0
		#for precision in [0.01]:
			#for maxiter in [100]:
				#result = EMData()
				#result.set_size(len(inc), len(inc),len(inc))
				#for dx in inc:
					#for dy in inc:
						#for az in inc:
							#t3d = Transform3D(EULER_EMAN,az,0,0)
							#t3d.set_posttrans(dx,dy)
							#e3 = EMData()
							#e3.set_size(n,n,1)
							#e3.process_inplace('testimage.squarecube', {'axis':'x', 'edge_length':10, 'fill':1})
							#e3.process_inplace('filter.lowpass.gauss', {'cutoff_pixels':n/8, 'sigma':0.5})
							
							##etmp = EMData()
							##etmp.set_size(n,n,1)
							##etmp.process_inplace('testimage.squarecube', {'axis':'x', 'edge_length':10, 'fill':1})
							##e3 = e3 + etmp
							
							#if ( i == 0 ):
								#e3.write_image("test_square.img")
							
							#e4 = e3.copy()
							
							#e3.rotate_translate(t3d)
											
							#e5 = e4.align('refine', e3, {'mode':2, 'dx':-1, 'dy':-1, 'az':-1, 'stepx':2, 'stepy':2, 'stepaz':2, 'maxiter':maxiter, 'precision':precision}, 'dot', {'normalize':1} )
							
							#if ( i == 0 ): (e5-e3).write_image("subtraction.img")
							
							#soln_dx = e5.get_attr("align.dx")
							#soln_dy = e5.get_attr("align.dy")
							#soln_az = e5.get_attr("align.az")
							
						
		
							#dx_diff = dx - soln_dx
							#dx_diff *= dx_diff
							#dy_diff = dy - soln_dy
							#dy_diff *= dy_diff
							#daz_diff = az - soln_az;
							#daz_diff *= daz_diff
		
							#x = int((dx + 1)*div)
							#y = int((dy + 1)*div)
							#z = int((az + 1)*div)
							
							#intensity = sqrt(dx_diff+dy_diff+daz_diff)
							
							#if (intensity > 1):
								#print "##############################################"
								#print "intensity diff is %f" %intensity
								#print "iter and precision is %d %f" %(maxiter, precision)
								#print "x align %f %f" %(dx, soln_dx)
								#print "y align %f %f" %(dy, soln_dy)
								#print "az align %f %f" %(az, soln_az)
								
							#result.set_value_at(x,y,z,intensity)
				
				#print "iter and precision are %d %f" %(maxiter, precision)
				#print "Stats are %f %f" %(result.get_attr("mean"), result.get_attr("sigma"))
				
		#result.write_image("result.mrc")

	
def test_main():
	test_support.run_unittest(TestAligner)

if __name__ == '__main__':
	test_main()
