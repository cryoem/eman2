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
import testlib
from pyemtbx.exceptions import *
from math import *
from optparse import OptionParser

IS_TEST_EXCEPTION = False

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
				for x in (32,33):
					size = (x,y,z)
					ref = test_image_3d(0,size)
					for i in range(0,20):
						e = ref.copy()
						dx = Util.get_frand(-3,3)
						dy = Util.get_frand(-3,3)
						dz = Util.get_frand(-3,3)
						e.translate(dx,dy,dz)
						g = e.align("translational",ref,{},"dot",{})
						t =  g.get_attr("xform.align3d")
						params = t.get_params("eman")
						self.failIf(fabs(params["tx"]+ dx) > 1)
						self.failIf(fabs(params["ty"] + dy) > 1)
						self.failIf(fabs(params["tz"] + dz) > 1)
						
						f = e.process("xform",{"transform":t})
						self.assertEqual(f==g,True)
						
								
		# Test 2D behavior
		for y in (32,33):
			for x in (32,33): # does not work yet for odd x - rotational footprint fails
				size = (x,y)
				images = []
				ref = test_image(0,size)
				for i in range(0,20):
					e = ref.copy()
					dx = Util.get_frand(-3,3)
					dy = Util.get_frand(-3,3)
					e.translate(dx,dy,0)
					g = e.align("translational",ref,{},"dot",{})
					t =  g.get_attr("xform.align2d")
					params = t.get_params("2d")
					self.failIf(fabs(params["tx"]+ dx) > 1)
					self.failIf(fabs(params["ty"] + dy) > 1)
					
					f = e.process("xform",{"transform":t})
					self.assertEqual(f==g,True)


	def test_RotationalAligner(self):
		"""test RotationalAligner ..........................."""
		e = EMData()
		e.set_size(32,32,1)
		e.process_inplace('testimage.noise.uniform.rand')
		
		e2 = EMData()
		e2.set_size(32,32,1)
		e2.process_inplace('testimage.noise.uniform.rand')
		
		e.align('rotational', e2)
		
		for y in [32,33]:
			for x in [32,33]:
				size = (x,y)
				ref = test_image(5,size)
				ref.translate(4,5,0) # give it handedness
				for i in range(0,20):
					e = ref.copy()
					az = Util.get_frand(0,360)
					e.transform(Transform({"type":"2d","alpha":az}))
					g = ref.align("rotational",e,{},"dot",{})
					t =  g.get_attr("xform.align2d")
					params = t.get_params("2d")
					result = fabs(params["alpha"] - az)
					#print g.get_attr("align.az"), az
					if result > 180 and result < 360:
						result = 360-result
					if result > 360: result = result-360
					self.failIf( result > 3 ) # 3 seems accurate enough
						
					# we have to do it this way because of double to float conversions
					f = ref.process("xform",{"transform":t})
					dif = f-g
					dif.process_inplace("math.absvalue")
					self.failIf(dif["mean"] > 0.01)
				
	def no_test_RotatePrecenterAligner(self):
		"""test RotatePrecenterAligner ......................"""

		e = EMData()
		e.set_size(32,32,1)
		e.process_inplace('testimage.noise.uniform.rand')
		
		e2 = EMData()
		e2.set_size(32,32,1)
		e2.process_inplace('testimage.noise.uniform.rand')
		
		e.align('rotate_precenter', e2)
		
		for y in [32,33]:
			for x in [32,33]:
				size = (x,y)
				ref = test_image(0,size)
				for i in range(0,20):
					e = ref.copy()
					az = Util.get_frand(0,360)
					e.transform(Transform({"type":"2d","alpha":az}))
					g = e.align("rotate_precenter",ref,{},"dot",{})
					t =  g.get_attr("xform.align2d")
					params = t.get_params("2d")
					result = fabs(params["alpha"] - az)
					print params["alpha"],az
					#print g.get_attr("align.az"), az
					if result > 180 and result < 360:
						result = 360-result
					if result > 360: result = result-360
					
					self.failIf( result > 3 ) # 3 seems accurate enough
		

	def test_RotateTranslateAligner(self):
		"""test RotateTranslateAligner ......................"""
		e = EMData()
		e.set_size(32,32,1)
		e.process_inplace('testimage.noise.uniform.rand')
		
		e2 = EMData()
		e2.set_size(32,32,1)
		e2.process_inplace('testimage.noise.uniform.rand')
		
		e.align('rotate_translate', e2, {'maxshift':1})
		
		for y in [32,33]:
			for x in [32,33]:
				size = (x,y)
				ref = test_image(0,size)
				ref.translate(4,5,0) # give it handedness
				for i in range(0,20):
					e = ref.copy()
					dx = Util.get_frand(-3,3)
					dy = Util.get_frand(-3,3)
					az = Util.get_frand(0,360)
					t = Transform({"type":"2d","alpha":az})
					t.set_trans(Vec2f(dx,dy))
					e.transform(t)
					g = ref.align("rotate_translate",e,{},"phase")
					t =  g.get_attr("xform.align2d")
#					params = t.get_params("2d")
#					self.failIf(fabs(params["tx"] - dx) > 1)
#					self.failIf(fabs(params["ty"] - dy) > 1)
#					
#					result = fabs(params["alpha"] - az)
#					#print g.get_attr("align.az"), az
#					if result > 180 and result < 360:
#						result = 360-result
#					if result > 360: result = result-360
#					self.failIf( result > 3 ) # 3 seems accurate enough
					
					# we have to do it this way because of double to float conversions
					f = ref.process("xform",{"transform":t})
					dif = f-g
					dif.process_inplace("math.absvalue")
					self.failIf(dif["mean"] > 0.01)


	def test_RotateFlipAligner(self):
		"""test RotateFlipAligner ..........................."""
		e = EMData()
		e.set_size(32,32,1)
		e.process_inplace('testimage.noise.uniform.rand')
		
		e2 = EMData()
		e2.set_size(32,32,1)
		e2.process_inplace('testimage.noise.uniform.rand')
		
		e.align('rotate_flip', e2, {'imask':2})
		
		for y in [64,65]:
			for x in [64,65]:
				size = (x,y)
				ref = test_image(0,size)
				ref.translate(4,5,0) # give it handedness
				#ref = test_image(0,size)
				for i in range(0,20):
					e = ref.copy()
					az = Util.get_frand(0,360)
					mirror = Util.get_irand(0,1)
					if mirror:e.process_inplace("xform.flip",{"axis":"x"})
					e.transform(Transform({"type":"2d","alpha":az}))
					g = ref.align("rotate_flip",e,{},"dot",{})
					t1 =  g.get_attr("xform.align2d")
					t1.invert()
					params = t1.get_params("2d")
					#(t1*t).printme()
					#print params
					#print az,mirror
					result = fabs(params["alpha"]+ az)
					if result > 180 and result < 360:
						result = 360-result
					if result > 360: result = result-360
					self.failIf( result > 3 ) # 3 seems accurate enough
					self.failIf( t1.get_mirror() != mirror)
					
#					f = e.process("xform",{"transform":t1})
#					dif = f-g
#					dif.process_inplace("math.absvalue")
#					print dif["mean"]
#					self.failIf(dif["mean"] > 0.01)
#	
	def run_rtf_aligner_test(self,aligner_name,aligner_params={},debug=False):
		
		for y in [64,65]:
			for x in [64,65]:
				size = (x,y)
				ref = test_image(0,size)
				ref = ref + test_image(5,size)
				ref.translate(4,5,0) # give it handedness
				for i in range(0,20):
					e = ref.copy()
					dx = Util.get_frand(-3,3)
					dy = Util.get_frand(-3,3)
					az = Util.get_frand(0,360)
					
					mirror = Util.get_irand(0,1)
					if mirror:e.process_inplace("xform.flip",{"axis":"x"})
					
					t = Transform({"type":"2d","alpha":az})
					t.set_pre_trans(Vec2f(dx,dy))
					e.transform(t)
					g = e.align(aligner_name,ref,aligner_params,"dot")
					t =  g.get_attr("xform.align2d")
					params = t.get_params("2d")
					if debug:
						print params
						print az,dx,dy,mirror
					self.failIf(fabs(params["tx"] + dx) > 2)
					self.failIf(fabs(params["ty"] + dy) > 2)
					
					result = fabs( (params["alpha"] + az) %360 )
					
					if result > 180 and result < 360:
						result = 360-result
					if result > 360: result = result-360
					self.failIf( result > 5 ) # 5 seems accurate enough
					self.failIf( t.get_mirror() != mirror)
	
	def test_RotateTranslateFlipAligner(self):
		"""test RotateTranslateFlip Aligner ................."""
		e = EMData()
		e.set_size(32,32,1)
		e.process_inplace('testimage.noise.uniform.rand')
		
		e2 = EMData()
		e2.set_size(32,32,1)
		e2.process_inplace('testimage.noise.uniform.rand')
		
		e.align('rotate_translate_flip', e2, {'maxshift':1})
		
		for y in [64,65]:
			for x in [64,65]:
				size = (x,y)
				ref = test_image(0,size)
				ref.translate(4,5,0) # give it handedness
				for i in range(0,20):
					mirror = Util.get_irand(0,1)
					
					e = ref.copy()
					dx = Util.get_frand(-3,3)
					dy = Util.get_frand(-3,3)
					az = Util.get_frand(0,360)
					t = Transform({"type":"2d","alpha":az})
					t.set_trans(Vec2f(dx,dy))
					t.set_mirror(mirror)
					e.transform(t)
					e.write_image("e.hdf",-1)
					g = ref.align("rotate_translate_flip",e,{},"phase")
					g.write_image("e.hdf",-1)
					t1 =  g.get_attr("xform.align2d")
#					params = t1.get_params("2d")
#					params2 = t.get_params("2d")
#					t.invert()
#					print params2
#					print params
#					print fabs(params["tx"] - params2["tx"])
#					print fabs(params["ty"] - params2["ty"])
#							
#					self.failIf(fabs(params["tx"] - params2["tx"]) > 1)
#					self.failIf(fabs(params["ty"] - params2["ty"]) > 1)
#					result = fabs(params["alpha"] - params2["alpha"])
#					#print g.get_attr("align.az"), az
#					if result > 180 and result < 360:
#						result = 360-result
#					if result > 360: result = result-360
#					self.failIf( result > 3 ) # 3 seems accurate enough
#					
#					
#					self.failIf( params2["mirror"] != params["mirror"] ) # 3 seems accurate enough
#					
					# we have to do it this way because of double to float conversions
					f = ref.process("xform",{"transform":t1})
					dif = f-g
					dif.process_inplace("math.absvalue")
					self.failIf(dif["mean"] > 0.01)
		testlib.safe_unlink('e.hdf')

	def test_RTF_slow_exhaustive_aligner(self):
		"""test RTFSlowExhaustiveAligner Aligner ............"""
		e = EMData()
		e.set_size(32,32,1)
		e.process_inplace('testimage.noise.uniform.rand')
		
		e2 = EMData()
		e2.set_size(32,32,1)
		e2.process_inplace('testimage.noise.uniform.rand')
		
		e.align('rtf_slow_exhaustive', e2, {'maxshift':10})
		
		#self.run_rtf_aligner_test("rtf_slow_exhaustive",debug=True)
		
	def test_RTF_exhaustive_aligner(self):
		"""test RTFExhaustiveAligner Aligner ................"""
		e = test_image()
		
		e2 = test_image()
		
		e.align('rtf_exhaustive', e2)
		#self.run_rtf_aligner_test("rtf_exhaustive")
		
	def test_RefineAligner(self):
		"""test RefineAligner ..............................."""
		e = EMData()
		e.set_size(32,32,1)
		e.process_inplace('testimage.noise.uniform.rand')
		
		e2 = EMData()
		e2.set_size(32,32,1)
		e2.process_inplace('testimage.noise.uniform.rand')
		
		e.align('refine', e2)
		t = Transform()
		e.align('refine', e2, {'mode':1, "xform.align2d":t})
		
		
		#for y in [32]:
			#for x in [32]:# does not work yet for odd x - rotational footprint fails
				#size = (x,y)
				#ref = test_image(0,size)
				#ref.translate(2,2,0)
				#for i in range(0,20):
					#e = ref.copy()
					#dx = Util.get_frand(-3,3)
					#dy = Util.get_frand(-3,3)
					#az = Util.get_frand(-5,5)
					#t = Transform({"type":"2d","alpha":az})
					#t.set_pre_trans(Vec2f(dx,dy))
					#e.transform(t)
					#g = e.align("refine",ref,{},"dot")
					
					#t =  g.get_attr("xform.align2d")
					#params = t.get_params("2d")
					##print params
					##print dx,dy,az
					#self.failIf(fabs(params["tx"] + dx) > 1)
					#self.failIf(fabs(params["ty"] + dy) > 1)
					
					#result = fabs(params["alpha"] + az)
					#if result > 180 and result < 360:
						#result = 360-result
					#if result > 360: result = result-360
					#self.failIf( result > 3 ) # 3 seems accurate enough
					
		
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
							#e3.process_inplace('filter.lowpass.gauss', {'cutoff_pixels':n/8, 'cutoff_abs':0.5})
							
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
				
	
def test_main():
    p = OptionParser()
    p.add_option('--t', action='store_true', help='test exception', default=False )
    global IS_TEST_EXCEPTION
    opt, args = p.parse_args()
    if opt.t:
        IS_TEST_EXCEPTION = True
    Log.logger().set_level(-1)  #perfect solution for quenching the Log error information, thank Liwei
    suite = unittest.TestLoader().loadTestsFromTestCase(TestAligner)
    unittest.TextTestRunner(verbosity=2).run(suite)

if __name__ == '__main__':
	test_main()
