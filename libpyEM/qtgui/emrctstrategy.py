#!/usr/bin/env python
#
# Author: John Flanagan, 04/08/2011 (jfflanag@bcm.edu)
# Edited by: Stephen Murray (scmurray@bcm.edu) May 2014
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
# Foundation, Inc., 59 Temple Place, Suite 330, Boston MA 02111-1307 USA
#
#
import numpy, math

class Strategy:
	''' This is a base class for the strategy to use for pcik event hadeling'''
	def __init__ (self, mediator):
		self.mediator = mediator
	
	# Run this function to do initial calculation, when widget is loaded
	def initial_calculations(self):
		raise NotImplementedError("Subclass must implement abstract method")
	
	# Run this function to respond to user input from the GUI(caller)
	def configure_strategy(self, caller):
		raise NotImplementedError("Subclass must implement abstract method")
	
	# Respond to signals form the GUI
	def handle_strategy_signal(self, signal):
		raise NotImplementedError("Subclass must implement abstract method")
	
	# Handle pick events
	def pickevent(self, caller, x, y):
		raise NotImplementedError("Subclass must implement abstract method")
	
	# Handle unpick events
	def unpickevent(self, box_num):
		raise NotImplementedError("Subclass must implement abstract method")
	
	# Handle move events
	def moveevent(self):
		raise NotImplementedError("Subclass must implement abstract method")

	def imagesaveevent(self, image):
		raise NotImplementedError("Subclass must implement abstract method")
		
class Strategy2IMGMan(Strategy):
	''' This is a derived class for the strategy to use for pcik event hadeling, more classes can be added'''
	def __init__ (self, mediator):
		Strategy.__init__(self, mediator)
	
	def initial_calculations(self):
		pass
	
	def configure_strategy(self, caller):
		pass
	
	def handle_strategy_signal(self, signal):
		pass
	
	def imagesaveevent(self, image):
		pass
		
	def pickevent(self, caller, x, y):
		if caller == self.mediator.untilt_win:
			if self.mediator.tilt_win.boxes.boxpopulation < self.mediator.untilt_win.boxes.boxpopulation:
				print "Error, you need to selct an untilted partilce pair, before you select a new tilted one"
				return False
		if caller == self.mediator.tilt_win:
			if (self.mediator.tilt_win.boxes.boxpopulation == 0 and self.mediator.untilt_win.boxes.boxpopulation == 0):
				print "Error, you first need to pick an untilted particle"
				return False
			if self.mediator.untilt_win.boxes.boxpopulation < self.mediator.tilt_win.boxes.boxpopulation:
				print "Error, you need to selct an untilted partilce pair, before you select a new tilted one"
				return False
		return True
		
	def unpickevent(self, caller, box_num):
		if caller == self.mediator.untilt_win:
			if len(self.mediator.tilt_win.boxes.boxlist)-1 >= box_num:
				self.mediator.tilt_win.boxes.remove_box(box_num, self.mediator.boxsize)
				self.mediator.tilt_win.update_mainwin()
		if caller == self.mediator.tilt_win:
			if len(self.mediator.untilt_win.boxes.boxlist)-1 >= box_num:
				self.mediator.untilt_win.boxes.remove_box(box_num, self.mediator.boxsize)
				self.mediator.untilt_win.update_mainwin()
		return True
		
	def moveevent(self):
		return True
		
class Strategy2IMGPair(Strategy):
	''' This is a derived class for the strategy to use for pcik event hadeling, more classes can be added'''
	def __init__ (self, mediator):
		Strategy.__init__(self, mediator)
		self.A = None
		self.invA = None
		self.tiltangle = None
		self.minpp_for_xform = 3
		self.cont_update_boxes = False
		self.centertilts = False
		
	# This function is called after boxes are loaded from the DB (to intialize the strategy
	def initial_calculations(self):
		if self.mediator.untilt_win.boxes.boxpopulation == self.mediator.tilt_win.boxes.boxpopulation:
				if self.mediator.untilt_win.boxes.boxpopulation > self.minpp_for_xform:
					self.compute_transform()
					self.compute_tilt_angle()
	
	# This function is called to configure or reconfigure the strategy alogorithm (sort of sets the state)
	def configure_strategy(self, caller):
		self.minpp_for_xform = caller.minpp_for_xform
		self.cont_update_boxes =  caller.updateboxes
		self.centertilts = caller.centertilts
	
	def handle_strategy_signal(self, signal):
		if signal == "updateboxes":
			self.on_update_boxes()
		if signal == "centerboxes":
			self.on_center_boxes()
		
	def pickevent(self, caller, x, y):
		# Pick tilted particle
		if caller == self.mediator.untilt_win:
			if self.mediator.untilt_win.boxes.boxpopulation == self.mediator.tilt_win.boxes.boxpopulation:
				if self.mediator.untilt_win.boxes.boxpopulation >= self.minpp_for_xform:

					# Compute transform
					self.compute_transform()
					
					# Compute tilt angle
					self.compute_tilt_angle()
					
					# Use the transfomration matrix to compute the tilted angle
					# I could just use the affine matrix, but better to use just the rotational part to reduce error
					currX = [x,y,1]
					currY = numpy.dot(self.A,currX)
					self.mediator.tilt_win.boxes.append_box(currY[0],currY[1])
				
					# Center boxes
					if self.centertilts:
						self.centerboxes(self.mediator.tilt_win)
						
					self.mediator.tilt_win.update_mainwin()
					self.mediator.tilt_win.update_particles()
					
					if self.cont_update_boxes:
						self.update_boxes()
					
		#pick untilted particle
		if caller == self.mediator.tilt_win:
			if (self.mediator.tilt_win.boxes.boxpopulation == 0 and self.mediator.untilt_win.boxes.boxpopulation == 0):
				print "Error, you first need to pick an untilted particle"
				return False
			if self.mediator.untilt_win.boxes.boxpopulation == self.mediator.tilt_win.boxes.boxpopulation:
				if self.mediator.tilt_win.boxes.boxpopulation >= self.minpp_for_xform:
					# Compute transform
					self.compute_transform()
					# Compute tilt angle
					self.compute_tilt_angle()
					
					# Compute untilted box coords
					currX = [x,y,1]
					currY = numpy.dot(self.invA,currX)	# Inverse of A
					self.mediator.untilt_win.boxes.append_box(currY[0],currY[1])
					
					# Center boxes
					if self.centertilts:
						self.centerboxes(self.mediator.untilt_win)
						
					self.mediator.untilt_win.update_mainwin()
					self.mediator.untilt_win.update_particles()
					
					if self.cont_update_boxes:
						self.update_boxes()	
		return True
	
	def compute_transform(self):
		#print "calc matrix"
		# Find the transformation matrix
		Xrow1 = []
		Xrow2 = []
		Xrow3 = []
		Yrow1 = []
		Yrow2 = []
		Yrow3 = []
		for i,box in enumerate(self.mediator.untilt_win.boxes.boxlist):
			Xrow1.append(box.x)
			Xrow2.append(box.y)
			Xrow3.append(1.0)
			Yrow1.append(self.mediator.tilt_win.boxes.boxlist[i].x)
			Yrow2.append(self.mediator.tilt_win.boxes.boxlist[i].y)
			Yrow3.append(1.0)
		X = numpy.array([Xrow1,Xrow2,Xrow3])
		Y = numpy.array([Yrow1,Yrow2,Yrow3])
		pinvX = numpy.linalg.pinv(X)	# Use pseduoinverse to find the best transformation matrix, A, in a least squares sense
		self.A = numpy.dot(Y,pinvX)
		self.invA = numpy.linalg.inv(self.A)
		
		# Enable buttons
		self.set_gui_buttons(True)
		self.compute_mask()

	def compute_tiltaxis(self):
		""" Must have already computed tilt angle for this to work!"""
		if self.A != None and self.tiltangle != None:
			rotA = numpy.array([[self.A[0,0],self.A[0,1]],[self.A[1,0],self.A[1,1]]])
			tan_phi = (rotA[0,0] - rotA[1,1]*math.cos(math.radians(self.tiltangle)))/(rotA[1,0]*math.cos(math.radians(self.tiltangle)) + rotA[0,1])
			phi = math.atan2((rotA[0,0] - rotA[1,1]*math.cos(math.radians(self.tiltangle))),(rotA[1,0]*math.cos(math.radians(self.tiltangle)) + rotA[0,1]))
			self.dphi = math.degrees(phi)
			
			sin_gamma = rotA[0,0]*math.sin(phi) + rotA[0,1]*math.cos(phi)
			try:
				gamma = math.asin(sin_gamma)
			except:
				gamma = math.pi/2
				
			self.dgamma = math.degrees(gamma) 
			#print rotA
			
			self.mediator.control_window.pair_picker_tool.tiltaxis.setText(("%3.2f"%self.dphi)+u'\u00B0')
			self.mediator.control_window.pair_picker_tool.gamma.setText(("%3.2f"%self.dgamma)+u'\u00B0')
			# Save tilt data
			self.mediator.tilt_win.boxes.save_tiltdata_to_db([self.tiltangle, self.dphi, self.dgamma])
			self.mediator.untilt_win.boxes.save_tiltdata_to_db([self.tiltangle, self.dphi, self.dgamma])
		
	def compute_mask(self):
		if self.A != None:
			v1 = numpy.dot(self.A,[0,0,1])
			v2 = numpy.dot(self.A,[self.mediator.untilt_win.win_xsize,0,1])
			v3 = numpy.dot(self.A,[self.mediator.untilt_win.win_xsize,self.mediator.untilt_win.win_ysize,1])
			v4 = numpy.dot(self.A,[0,self.mediator.untilt_win.win_ysize,1])
			self.mediator.tilt_win.paint_mask(v1[0],v1[1], v2[0],v2[1],v3[0],v3[1],v4[0],v4[1])
			self.mediator.tilt_win.update_mainwin()
		
			vinv1 = numpy.dot(self.invA,[0,0,1])
			vinv2 = numpy.dot(self.invA,[self.mediator.tilt_win.win_xsize,0,1])
			vinv3 = numpy.dot(self.invA,[self.mediator.tilt_win.win_xsize,self.mediator.tilt_win.win_ysize,1])
			vinv4 = numpy.dot(self.invA,[0,self.mediator.tilt_win.win_ysize,1])
			self.mediator.untilt_win.paint_mask(vinv1[0],vinv1[1], vinv2[0],vinv2[1],vinv3[0],vinv3[1],vinv4[0],vinv4[1])
			self.mediator.untilt_win.update_mainwin()
	
	def compute_tilt_angle(self):
		#self.compute_tilt_angle_phil()
		
		if self.A != None:
			# Use the transformation matrix to compute the tilt angle
			rotA = numpy.array([[self.A[0,0],self.A[0,1]],[self.A[1,0],self.A[1,1]]])
			detA = numpy.linalg.det(self.A)	# The determinate is is COS of the tilt angle
			try:
				self.tiltangle = math.degrees(math.acos(detA))
				self.mediator.control_window.pair_picker_tool.tiltangle.setText(("%3.2f"%self.tiltangle)+u'\u00B0')
			except:
				self.mediator.control_window.pair_picker_tool.tiltangle.setText("Det(A) > 1")
			self.compute_tiltaxis()
		
	def compute_tilt_angle_phil(self):
		if self.A != None:
			# SVD of A
			rotA = numpy.array([[self.A[0,0],self.A[0,1]],[self.A[1,0],self.A[1,1]]])
			U, D, V = numpy.linalg.svd(rotA)
			# single values are ranked by numpy
			self.tiltangle = math.degrees(math.acos(D[1]))
			self.mediator.control_window.pair_picker_tool.tiltangle.setText(("%3.2f"%self.tiltangle)+u'\u00B0')
			# compute tilt axis
			self.dphi = math.degrees(math.atan2(U[0][1],U[0][0]))
			self.dgamma = math.degrees(math.atan2(V[0][1],V[0][0]))
			# save and display data
			self.mediator.control_window.pair_picker_tool.tiltaxis.setText(("%3.2f"%self.dphi)+u'\u00B0')
			self.mediator.control_window.pair_picker_tool.gamma.setText(("%3.2f"%self.dgamma)+u'\u00B0')
			# Save tilt data
			self.mediator.tilt_win.boxes.save_tiltdata_to_db([self.tiltangle, self.dphi, self.dgamma])
			self.mediator.untilt_win.boxes.save_tiltdata_to_db([self.tiltangle, self.dphi, self.dgamma])
			
			print D
			
	def update_boxes(self):
		for i,box in enumerate(self.mediator.untilt_win.boxes.boxlist):
			# Compute tilted box
			boxX = [self.mediator.untilt_win.boxes.boxlist[i].x,self.mediator.untilt_win.boxes.boxlist[i].y,1.0]
			boxY = numpy.dot(self.A,boxX)
			# Set tilted box
			self.mediator.tilt_win.boxes.boxlist[i].x = boxY[0]
			self.mediator.tilt_win.boxes.boxlist[i].y = boxY[1]
			self.mediator.tilt_win.boxes.shapelist[i] = None
			self.mediator.tilt_win.boxes.labellist[i] = None
			
		# Update tilted window
		self.mediator.tilt_win.boxes.save_boxes_to_db()
		self.mediator.tilt_win.window.update_shapes(self.mediator.tilt_win.boxes.get_box_shapes(self.mediator.boxsize))
		self.mediator.tilt_win.window.updateGL()
	
	def centerboxes(self, window):
		image = window.boxes.boxlist[len(window.boxes.boxlist)-1].get_image(window.filename,self.mediator.boxsize,"normalize.edgemean")
		centerimg = image.process("xform.centeracf")
		tv = centerimg.get_attr("xform.align2d").get_trans()
		window.boxes.boxlist[len(window.boxes.boxlist)-1].move(-tv[0],-tv[1])
		window.boxes.shapelist[len(window.boxes.boxlist)-1] = None
		window.boxes.labellist[len(window.boxes.boxlist)-1] = None
	
	def imagesaveevent(self, image):
		if self.A != None:
			if image.has_attr("tiltaxis"): image.set_attr("tiltaxis", self.dphi)
			if image.has_attr("tiltgamma"): image.set_attr("tiltgamma", self.dgamma)
			if image.has_attr("tiltangle"): image.set_attr("tiltangle", self.tiltangle)
	
	def unpickevent(self, caller, box_num):
		if caller == self.mediator.untilt_win:
			if len(self.mediator.tilt_win.boxes.boxlist)-1 >= box_num:
				self.mediator.tilt_win.boxes.remove_box(box_num, self.mediator.boxsize)
				if self.mediator.untilt_win.boxes.boxpopulation <= self.minpp_for_xform or self.mediator.tilt_win.boxes.boxpopulation <= self.minpp_for_xform:
					self.set_gui_buttons(False)
					self.A = None
					self.invA = None
				self.mediator.tilt_win.update_mainwin()
		if caller == self.mediator.tilt_win:
			if len(self.mediator.untilt_win.boxes.boxlist)-1 >= box_num:
				self.mediator.untilt_win.boxes.remove_box(box_num, self.mediator.boxsize)
				if self.mediator.untilt_win.boxes.boxpopulation <= self.minpp_for_xform or self.mediator.tilt_win.boxes.boxpopulation <= self.minpp_for_xform:
					self.set_gui_buttons(False)
					self.A = None
					self.invA = None
				self.mediator.untilt_win.update_mainwin()
			
		return True
	
	# This toggle the contol buttons on and off (also controls the mask)
	def set_gui_buttons(self, toggle):
		self.mediator.control_window.pair_picker_tool.mask_combobox.setEnabled(toggle)
		self.mediator.control_window.pair_picker_tool.upboxes_but.setEnabled(toggle)
		self.mediator.control_window.pair_picker_tool.centerboxes_but.setEnabled(toggle)
		if toggle == False:
			self.mediator.tilt_win.boxes.add_mask(None)
			self.mediator.untilt_win.boxes.add_mask(None)
		
	def moveevent(self):
		return True

	def on_update_boxes(self):
		self.compute_transform()
		self.compute_tilt_angle()
		self.update_boxes()
		
	def on_center_boxes(self):
		for i,image in enumerate(self.mediator.tilt_win.boxes.get_particle_images(self.mediator.tilt_win.filename,self.mediator.boxsize)):
			centerimg = image.process("xform.centeracf")
			tv = centerimg.get_attr("xform.align2d").get_trans()
			self.mediator.tilt_win.boxes.boxlist[i].move(-tv[0],-tv[1])
			self.mediator.tilt_win.boxes.shapelist[i] = None
			self.mediator.tilt_win.boxes.labellist[i] = None
		# Update tilted window
		self.mediator.tilt_win.boxes.save_boxes_to_db()
		self.mediator.tilt_win.window.update_shapes(self.mediator.tilt_win.boxes.get_box_shapes(self.mediator.boxsize))
		self.mediator.tilt_win.window.updateGL()