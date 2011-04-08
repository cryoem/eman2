#!/usr/bin/env python
#
# Author: Steven Ludtke, 11/01/2007 (sludtke@bcm.edu)
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
		self.minpp_for_xfrom = 3
		self.cont_update_boxes = False
	
	def initial_calculations(self):
		if self.mediator.untilt_win.boxes.boxpopulation == self.mediator.tilt_win.boxes.boxpopulation:
				if self.mediator.untilt_win.boxes.boxpopulation >= self.minpp_for_xfrom:
					self.compute_transform()
					self.compute_tilt_angle()
					self.mediator.control_window.pair_picker_tool.upboxes_but.setEnabled(True)
					
	def configure_strategy(self, caller):
		self.minpp_for_xfrom = caller.minpp_for_xfrom
		self.cont_update_boxes =  caller.updateboxes
	
	def handle_strategy_signal(self, signal):
		if signal == "updateboxes":
			self.on_update_boxes()
		
	def pickevent(self, caller, x, y):
		# Pick tilted particle
		if caller == self.mediator.untilt_win:
			if self.mediator.untilt_win.boxes.boxpopulation == self.mediator.tilt_win.boxes.boxpopulation:
				if self.mediator.untilt_win.boxes.boxpopulation >= self.minpp_for_xfrom:

					# Compute transform
					self.compute_transform()
					
					# Compute tilt angle
					self.compute_tilt_angle()
					
					# Use the transfomration matrix to compute the tilted angle
					# I could just use the affine matrix, but better to use just the rotational part to reduce error
					currX = [x,y,1]
					currY = numpy.dot(self.A,currX)
					self.mediator.tilt_win.boxes.append_box(currY[0],currY[1])
					self.mediator.tilt_win.update_mainwin()	
					
					#Talk back to GUI:
					self.mediator.control_window.pair_picker_tool.upboxes_but.setEnabled(True)
					
					if self.cont_update_boxes:
						self.update_boxes()
		#pick untilted particle
		if caller == self.mediator.tilt_win:
			if (self.mediator.tilt_win.boxes.boxpopulation == 0 and self.mediator.untilt_win.boxes.boxpopulation == 0):
				print "Error, you first need to pick an untilted particle"
				return False
			if self.mediator.untilt_win.boxes.boxpopulation == self.mediator.tilt_win.boxes.boxpopulation:
				if self.mediator.tilt_win.boxes.boxpopulation >= self.minpp_for_xfrom:
					# Compute transform
					self.compute_transform()
					# Compute tilt angle
					self.compute_tilt_angle()
					
					currX = [x,y,1]
					currY = numpy.dot(numpy.linalg.inv(self.A),currX)	# Inverse of A
					self.mediator.untilt_win.boxes.append_box(currY[0],currY[1])
					self.mediator.untilt_win.update_mainwin()
					
					#Talk back to GUI:
					self.mediator.control_window.pair_picker_tool.upboxes_but.setEnabled(True)
					
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
	
	def compute_tilt_angle(self):
		# Use the transformation matrix to compute the tilt angle
		rotA = numpy.array([[self.A[0,0],self.A[0,1]],[self.A[1,0],self.A[1,1]]])
		detA = numpy.linalg.det(self.A)	# The determinate is is COS of the tilt angle
		try:
			tiltangle = math.degrees(math.acos(detA))
			self.mediator.control_window.pair_picker_tool.tiltangle.setText(("%3.2f"%tiltangle)+u'\u00B0')
		except:
			self.mediator.control_window.pair_picker_tool.tiltangle.setText("Det(A) > 1")	
		
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

	def on_update_boxes(self):
		self.compute_transform()
		self.compute_tilt_angle()
		self.update_boxes()