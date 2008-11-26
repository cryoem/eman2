#!/usr/bin/env python

#
# Author: David Woolford 11/25/08 (woolford@bcm.edu
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

from emapplication import EMStandAloneApplication
from emimage3dsym import EM3DSymViewerModule,EMSymInspector
from emglobjects import EMImage3DGUIModule
from PyQt4 import QtGui,QtCore
from OpenGL import GL,GLU,GLUT
from OpenGL.GL import *
from OpenGL.GLU import *
import weakref
from optparse import OptionParser
from EMAN2 import Util, E2init, E2end,EMANVERSION,is_2d_image_mx, EMUtil, db_open_dict, EMData
from emimagemx import EMImageMXModule
import os
import sys

def get_eulers_from(filename):
	eulers = []
	n = EMUtil.get_image_count(filename)
	for i in range(n):
		h = get_header(filename,i)
		try: p = h["xform.projection"]
		except:
			print "image",i,"doesn't have the xform.projection attribute"
			return None
		
		eulers.append(p)
		
	return eulers


def get_ptcl_from(filename):
	ptcl = []
	n = EMUtil.get_image_count(filename)
	for i in range(n):
		h = get_header(filename,i)
		try: p = h["ptcl_repr"]
		except:
			print "image",i,"doesn't have the ptcl_repr attribute"
			return None
		ptcl.append(p)
#		
#	norm_ptcl = normalize_ptcl(ptcl)
	return ptcl

def check_projections_match_averages(projection_file, average_file):
	fine, message = is_2d_image_mx(projection_file)
	fine2, message2 = is_2d_image_mx(average_file)
	
	if not fine: print message # just print the messages first
	if not fine2: print message2
	
	if not fine or not fine2: return None
	
	if EMUtil.get_image_count(projection_file) != EMUtil.get_image_count(average_file):
		print "image count for projection and averages files don't match", EMUtil.get_image_count(projection_file), EMUtil.get_image_count(average_file)
		return None
	
	eulers = []
	ptcl = []
	for i in range(EMUtil.get_image_count(projection_file)):
		h1 = get_header(projection_file)
		p1 = h1["xform.projection"] # projection should definitely have this attribute
		eulers.append(p1)
		
		h2 = get_header(average_file)
		p2 = h2["ptcl_repr"] # projection should definitely have this attribute
		ptcl_append(p2)
		
	return eulers, ptcl

def normalize_ptcl(ptcl):
	mn = min(ptcl)
	mx = max(ptcl)
	diff = float(mx-mn)
	norm = [ (val-mn)/(diff) for val in ptcl ]
	return norm		
	
def get_header(filename,i):
	if filename[0:4] == "bdb:":
		db = db_open_dict(filename)
		return db.get_header(i)
	else:
		read_header_only = True
		e = EMData()
		e.read_image(filename,i,read_header_only)
		return e.get_attr_dict()

def get_normalize_colors(ptcls):
	mn = min(ptcls)
	mx = max(ptcls)
	df = float(mx-mn)
	colors = []
	for val in ptcls:
		val = (val-mn)/(df)
		if val < 0.5:
			frac = val/0.5
			colors.append((1.0,frac,frac,1.0))
		elif val > 0.5:
			frac = (-val+1.0)/0.5
			colors.append((frac,frac,1.0,1.0))
		else:
			colors.append((1.0,1.0,1.0,1.0))
			
	return colors
	
def main():
	progname = os.path.basename(sys.argv[0])
	usage = """%prog [options] <class image>
	
Asymmetric unit viewer for EMAN2."""

	parser = OptionParser(usage=usage,version=EMANVERSION)
	
	parser.add_option("--sym",type="string",help="Guiding symmetry",default="icos")
	parser.add_option("--projections",type="string",help="File containing projections",default=None)
	parser.add_option("--averages",type="string",help="File containing averages",default=None)
	
	(options, args) = parser.parse_args()
	
	eulers = None
	ptcls = None
	if options.projections != None and options.averages != None:
		eulers = check_projections_match_averages(options.projections,options.averages)
		if eulers == None:
			print "failed to get eulers from",options.projections,options.averages
			sys.exit(1)

	elif options.projections != None:
		fine, message = is_2d_image_mx(options.projections)
		if not fine:
			print message
			sys.exit(1)
		else: 
			eulers = get_eulers_from(options.projections)
			if eulers == None:
				print "failed to get eulers from",options.projections
				sys.exit(1)
	elif options.averages != None:
		fine, message = is_2d_image_mx(options.averages)
		if not fine:
			print message
			sys.exit(1)
		else: 
			eulers = get_eulers_from(options.averages)
			if eulers == None:
				print "failed to get eulers from",options.averages
				sys.exit(1)
			ptcls = get_ptcl_from(options.averages)
			if ptcls == None:
				print "failed to get ptcls from",options.averages
				sys.exit(1)
	elif len(args) > 1:
		print "not supported yet"
		sys.exit(1)
	
	
	logid=E2init(sys.argv)
	
	em_app = EMStandAloneApplication()
	window = EMAsymmetricUnitViewer(application=em_app)
	if eulers != None:
		window.specify_eulers(eulers)
	else:
		print "no eulers"
	if ptcls != None:
		window.specify_colors(get_normalize_colors(ptcls))
		
	window.projection_file = options.projections
	window.average_file = options.averages
	window.set_sym(options.sym)
	em_app.show()
	em_app.execute()
	
	E2end(logid)

class InputEventsHandler:
	'''
	Perhaps the final installation in terms of what I think is the best design for the mouse events
	handling class. Others exist in emimage2d, emimagemx, and e2boxer. Eventually they should all use the same approach, and 
	I vote for using this one
	'''
	def __init__(self,parent):
		self.parent = weakref.ref(parent)
		
	def mousePressEvent(self,event):
		pass
	
	def mouseReleaseEvent(self,event):
		pass
	
	def mouseMoveEvent(self,event):
		pass
	
	def mouseDoubleClickEvent(self,event):
		pass
	
	def keyPressEvent(self,event):
		pass

	def wheelEvent(self,event):
		pass
	
class InputEventsManager(InputEventsHandler):
	def __init__(self):
		InputEventsHandler.__init__(self,self)
		self.current_events_handler = None
		
	def mousePressEvent(self,event):
		if self.current_events_handler != None:
			self.current_events_handler.mousePressEvent(event)
	
	def mouseReleaseEvent(self,event):
		if self.current_events_handler != None:
			self.current_events_handler.mouseReleaseEvent(event)
	
	def mouseMoveEvent(self,event):
		if self.current_events_handler != None:
			self.current_events_handler.mouseMoveEvent(event)
	
	def mouseDoubleClickEvent(self,event):
		if self.current_events_handler != None:
			self.current_events_handler.mouseDoubleClickEvent(event)
	
	def keyPressEvent(self,event):
		if self.current_events_handler != None:
			self.current_events_handler.keyPressEvent(event)
	
	def wheelEvent(self,event):
		if self.current_events_handler != None:
			self.current_events_handler.wheelEvent(event)
	
class ClassOrientationEvents(InputEventsHandler,QtCore.QObject): 
	def __init__(self,parent):
		InputEventsHandler.__init__(self,parent)
		QtCore.QObject.__init__(self)
		self.old_intersection = -1
		self.old_color = None
		
	def mouseReleaseEvent(self,event):
		m,p,v = self.parent().model_matrix.tolist(),self.parent().vdtools.wproj.tolist(),self.parent().vdtools.wview.tolist()
		
		model_matrix = []
		proj_matrix = []
		view_matrix = []
		
		
		for val in m:
			if isinstance(val,list): model_matrix.extend(val)
			else: modul_matrix.append(val)
			
		for val in p:
			if isinstance(val,list): proj_matrix.extend(val)
			else: proj_matrix.append(val)
			
		for val in v:
			if isinstance(val,list): view_matrix.extend(val)
			else: view_matrix.append(val)
			
		points = self.parent().points
		mouse_x = event.x()
		mouse_y = view_matrix[-1]-event.y()
		intersection = Util.nearest_projected_points(model_matrix,proj_matrix,view_matrix,points,float(mouse_x),float(mouse_y),6.0)
		
		new_colors = {}
		
		if intersection >= 0:
			new_colors[intersection] = (1.0,1.0,0,1)
			
			if self.old_intersection >= 0:
				new_colors[self.old_intersection] = self.old_color
			
			self.old_intersection = intersection
			self.old_color = self.parent().point_colors[intersection]			
		
		else:
			if self.old_intersection >= 0:
				new_colors[self.old_intersection] = self.old_color
				self.old_intersection = -1
			else: return
		
		if len(new_colors) > 0:
			self.parent().set_point_colors(new_colors)
			self.parent().updateGL()
			if intersection >= 0:self.emit(QtCore.SIGNAL("point_selected"),intersection)	
		
class NavigationEvents(InputEventsHandler):
	def __init__(self,parent):
		InputEventsHandler.__init__(self,parent)
		
	def mousePressEvent(self,event):
		EMImage3DGUIModule.mousePressEvent(self.parent(),event)
	
	def mouseReleaseEvent(self,event):
		EMImage3DGUIModule.mouseReleaseEvent(self.parent(),event)
	
	def mouseMoveEvent(self,event):
		EMImage3DGUIModule.mouseMoveEvent(self.parent(),event)
	
	def mouseDoubleClickEvent(self,event):
		EMImage3DGUIModule.mouseDoubleClickEvent(self.parent(),event)
	
	def keyPressEvent(self,event):
		EMImage3DGUIModule.keyPressEvent(self.parent(),event)
		
	def wheelEvent(self,event):
		EMImage3DGUIModule.wheelEvent(self.parent(),event)
		


class EMAsymmetricUnitViewer(InputEventsManager,EM3DSymViewerModule):
	def __init__(self,application):
		
		EM3DSymViewerModule.__init__(self,application)
		InputEventsManager.__init__(self)
		
		self.__init_events_handlers()
		self.projection_file = None
		self.average_file = None
		self.mx_viewer = None
		#self.mousePressEvent = InputEventsManager.mousePressEvent
		
	def set_projection_file(self,projection_file): self.projection_file = projection_file
	def get_inspector(self):
		if not self.inspector : self.inspector=EMAsymmetricUnitInspector(self)
		return self.inspector
	
	def __init_events_handlers(self):
		self.events_mode = "navigate"
		self.events_handlers = {}
		self.events_handlers["navigate"] = NavigationEvents(self)
		self.events_handlers["inspect"] = ClassOrientationEvents(self)
		self.connect(self.events_handlers["inspect"],QtCore.SIGNAL("point_selected"), self.au_point_selected)
		self.current_events_handler = self.events_handlers["navigate"]
	
	def au_point_selected(self,i):
		try:
			a = EMData(self.projection_file,i)
		except: a = None
		
		try:
			b = EMData(self.average_file,i)
			if a == None:
				try:
					a = EMData(b.get_attr("projection_image"),b.get_attr("projection_image_idx"))
				except:
					a = None
		except: b = None
		if a == None and b == None: return
		
		first = False
		if self.mx_viewer == None:
			first = True
			self.mx_viewer = EMImageMXModule(data=None,application=self.application())
			self.application().show_specific(self.mx_viewer)
			
		
		disp = []
		if a != None: disp.append(a)
		if b != None: disp.append(b)

		self.mx_viewer.set_data(disp)
		if first: self.mx_viewer.optimally_resize()
		self.mx_viewer.updateGL()
	def set_events_mode(self,mode):
		if not self.events_handlers.has_key(mode):
			print "error, unknown events mode", mode
			return
		
		else:
			self.current_events_handler = self.events_handlers[mode]
			
		
	
#	def get_inspector(self):
#		pass
	
class EMAsymmetricUnitInspector(EMSymInspector):
	def __init__(self,target) :
		EMSymInspector.__init__(self,target)
		hbl = QtGui.QHBoxLayout()
		self.mouse_modes = ["navigate","inspect"]
		buttons = []
		check_first = True
		self.mouse_mode_buttons = []
		for choice in self.mouse_modes:
			button = QtGui.QRadioButton(str(choice))
			self.mouse_mode_buttons.append(button) # for determining which button is on
 			if check_first:
				button.setChecked(True)
				check_first = False
			hbl.addWidget( button)
			buttons.append(button)
			self.connect(button,QtCore.SIGNAL("clicked(bool)"),self.on_mouse_mode_clicked)
		
		
		groupbox = QtGui.QGroupBox("Mouse mode")
		groupbox.setToolTip("Set the current mouse mode")
		groupbox.setLayout(hbl)
		
		self.vbl2.addWidget(groupbox)

	def on_mouse_mode_clicked(self,bool):
		for button in self.mouse_mode_buttons:
			if button.isChecked():
				self.target().set_events_mode(str(button.text()))

if __name__ == '__main__':
	main()