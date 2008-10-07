#!/usr/bin/env python

#
# Author: David Woolford 10/01/2008 (woolford@bcm.edu)
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

from PyQt4 import QtGui
import sys

class EMGUIModule:
	FTGL = "ftgl"
	GLUT = "glut"
	def __init__(self,application=None): 
		self.application = application
		if application != None: application.attach_child(self)
		
	def set_app(self,application): self.application = application
	def get_app(self): return self.application
	def get_qt_widget(self): raise
	
	#def load_font_renderer(self):
		#try:
			#self.font_render_mode = EMGUIModule.FTGL
			#self.font_renderer = get_3d_font_renderer()
			#self.font_renderer.set_face_size(16)
			#self.font_renderer.set_font_mode(FTGLFontMode.TEXTURE)
		#except:
			#self.font_render_mode = EMGUIModule.GLUT
	

class EMApplication:
	def __init__(self,qt_application_control=True):
		if qt_application_control:
			self.app = QtGui.QApplication(sys.argv)
		else: self.app = None
	def attach_child(self,child):
		raise
	
	def detach_child(self,child):
		raise
	
	def show(self):
		raise
	
	def show_specific(self,object):
		raise
	
	def execute(self):
		if self.app != None:
			return sys.exit(self.app.exec_())
		else: return
		
	def quit(self):
		if self.app != None:
			self.app.quit()
			
	def setOverrideCursor(self,cursor_type):
		if self.app != None:
			self.app.setOverrideCursor(cursor_type)

class EMStandAloneApplication(EMApplication):
	
	def __init__(self,qt_application_control=True):
		EMApplication.__init__(self,qt_application_control)
		
		self.children = []
		
	def attach_child(self,child):
		for i in self.children:
			if i == child:
				print "error, can't attach the same child twice",child
				return
			
		self.children.append(child)
	
	def show(self):
		for child in self.children:
			widget = child.get_qt_widget()
			if widget.isVisible() == False:
				widget.show()
	
	def show_specific(self,child):
		widget = child.get_qt_widget()
		if widget.isVisible() == False:
			widget.show()
		

	def close_child(self,child):
		for child_ in self.children:
			if child == child_:
				widget = child.get_qt_widget()
				widget.close()
				return
			
		print "error, attempt to close a child that did not belong to this application"
		

class EMQtWidgetModule(EMGUIModule):
	def __init__(self,qt_widget,application):
		EMGUIModule.__init__(self,application)
		self.qt_widget = qt_widget
	
	def get_qt_widget(self):
		return self.qt_widget

