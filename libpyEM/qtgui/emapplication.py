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

from PyQt4 import QtGui, QtCore, QtOpenGL
from PyQt4.QtCore import Qt
import sys
from emimageutil import EMParentWin
from EMAN2 import remove_directories_from_name, get_image_directory,get_3d_font_renderer, E2end,get_platform
import EMAN2db
import weakref
from libpyGLUtils2 import *

try: from PyQt4 import QtWebKit
except: pass

class ModuleEventsManager: 
	'''
	Coordinates events of the various modules.
	To begin with this is only the close event, then I added the idle event
	'''
	def __init__(self,target,module): #TODO: determine if this class should be deprecated and its functionality integrated into self.target
		self.target = weakref.ref(target)
		self.module = weakref.ref(module)
		try:
			emitter = self.module().emitter()
			print "It has an emitter() function! Check whether it should be converted to a QWidget subclass."
		except AttributeError:
			emitter = self.module() #Ross's hack to get this to work with QWidget's as well
			
		QtCore.QObject.connect(emitter, QtCore.SIGNAL("module_closed"), self.module_closed)
		QtCore.QObject.connect(emitter, QtCore.SIGNAL("module_idle"), self.module_idle)
	
		QtCore.QObject.connect(emitter, QtCore.SIGNAL("ok"), self.module_ok) # yes, redundant, but time is short
		QtCore.QObject.connect(emitter, QtCore.SIGNAL("cancel"), self.module_cancel)# yes, redundant, but time is short
		
	
	def module_closed(self):
		self.disconnect_object()
		self.target().module_closed(self.module())
		
	def module_idle(self):
		self.disconnect_object()
		self.target().module_idle(self.module())
		
	def module_ok(self,*args,**kargs):
		self.disconnect_object()
		self.module().close()
		
	def module_cancel(self):
		self.disconnect_object()
		self.module().close()
	
	
	def disconnect_object(self):
		try:
			emitter = self.module().emitter()
		except AttributeError:
			emitter = self.module() #Ross's hack to get this to work with QWidget's as well
			
		QtCore.QObject.disconnect(emitter, QtCore.SIGNAL("module_closed"), self.module_closed)
		QtCore.QObject.disconnect(emitter, QtCore.SIGNAL("module_idle"), self.module_idle)
	
		QtCore.QObject.disconnect(emitter, QtCore.SIGNAL("ok"), self.module_ok) # yes, redundant, but time is short
		QtCore.QObject.disconnect(emitter, QtCore.SIGNAL("cancel"), self.module_cancel)# yes, redundant, but time is short

class EMGLWidget(QtOpenGL.QGLWidget):
	"""
	This class encapsulates the use of the EMParentWin to provide a status bar with a size grip on Mac.
	It also handles much of the inspector behavior, displays help in a web browser, and provides 
	a self.busy attribute to prevent updateGL() from redrawing before all changes to display parameters are in place. 
	"""
	
	def hide(self):
		if self.qt_parent:
			self.qt_parent.hide()
			
	def resize(self, w, h):
		if self.qt_parent:
			QtOpenGL.QGLWidget.resize(self, w, h)
			if get_platform()=="Darwin" : self.qt_parent.resize(w, h+22)
			else : self.qt_parent.resize(w+4, h+4)
			
	def show(self):
		if self.qt_parent:
			self.qt_parent.show()
			
	def setWindowTitle(self, title):
		if self.qt_parent:
			self.qt_parent.setWindowTitle(title)
		else:
			self.setWindowTitle(title)
			
	def __init__(self, parent=None,enable_timer=False, application_control=True,winid=None):
		QtOpenGL.QGLWidget.__init__(self,parent)
		self.qt_parent = EMParentWin(self, enable_timer)
		
		self.inspector = None # a Qt Widget for changing display parameters, setting the data, accessing metadata, etc.
		self.winid=winid # a 'unique' identifier for the window used to restore locations on the screen
		
		self.image_change_count =  0# this is important when the user has more than one display instance of the same image, for instance in e2.py if 
		app = get_application()
		if app != None and application_control:
			app.attach_child(self)
		
		self.application_control = application_control
		self.file_name = ""
		self.disable_inspector = False
		
		self.makeCurrent()
		self.font_renderer = get_3d_font_renderer()
		if self.font_renderer == None : 
			print "FTGL not initialized. Please make sure you have FTGL installed, and configured for compilation in CMake. If using a binary, please report a problem to sludtke@bcm.edu."
			sys.exit(1)
		self.font_renderer.set_face_size(16)
		self.font_renderer.set_font_mode(FTGLFontMode.TEXTURE)
			
		self.busy = False #updateGL() does nothing when self.busy == True
	
	def closeEvent(self, event):
		if self.inspector:
			self.inspector.close()
		QtOpenGL.QGLWidget.closeEvent(self, event)
		self.qt_parent.close()
		self.emit(QtCore.SIGNAL("module_closed")) # this could be a useful signal, especially for something like the selector module, which can potentially show a lot of images but might want to close them all when it is closed
		event.accept()
		
	def display_web_help(self,url="http://blake.bcm.edu/emanwiki/e2display"):
		try:
			try:
				import webbrowser
				webbrowser.open(url)
			except:
				try:
					test = self.browser
				except: 
					self.browser = QtWebKit.QWebView()
					self.browser.load(QtCore.QUrl())
					self.browser.resize(800,800)
				
				if not self.browser.isVisible(): self.browser.show(url)
		except:
			pass

	def enable_inspector(self,val=True): 
		self.disable_inspector = not val
		
			
	def show_inspector(self,force=0):
		if self.disable_inspector: 
			return
		
		self.emit(QtCore.SIGNAL("inspector_shown")) # debug only
		app = get_application()
		if app == None:
			print "can't show an inspector with having an associated application"
		
		if not force and self.inspector==None : 
			return
		if not self.inspector : 
			self.inspector = self.get_inspector()
			if self.inspector == None: 
				return # sometimes this happens
		if not app.child_is_attached(self.inspector):
			app.attach_child(self.inspector)
		app.show_specific(self.inspector)
		
	def update_inspector_texture(self):
		if self.inspector != None:
			self.inspector.update()
			
	def updateGL(self):
		if self.busy:
			return
		else:
			QtOpenGL.QGLWidget.updateGL(self)

class EMInstance:
	'''
	Holds a reference to an instance, supports a static interface
	'''
	instance = None
	
	def __init__(self):pass
	
	def get_instance(): return EMInstance.instance
	get_instance = staticmethod(get_instance)
	def set_instance(instance): EMInstance.instance = instance
	set_instance = staticmethod(set_instance)
	
em_app_instance = EMInstance()

get_application = em_app_instance.get_instance
#def get_application() : return QtGui.qApp


class EMApp(QtGui.QApplication):
	def __init__(self):
		self.children = []
		
		# Stuff for display synchronization in e2.py
		self.timer_function = None
		self.tmr = None
		
		QtGui.QApplication.__init__(self, sys.argv)
		
		style=QtGui.QStyleFactory.create("Cleanlooks")
		
		if style==None:
			print "Note: standard Cleanlooks style not available, controls may be distorted. Using ",
			
			# the first one should work, but we have the loop, just in case
			for s in list(QtGui.QStyleFactory.keys()):
				style=QtGui.QStyleFactory.create(s)
				if style!=None: 
					print s
					break

		if style!=None: self.setStyle(style)
		
		if em_app_instance.get_instance() == None:
			em_app_instance.set_instance(self)
	

	def child_is_attached(self,query_child):
		return (query_child in self.children)

	def detach_child(self,child):
		if child not in self.children:
			print "EMApp.detach_child() error, EMApp instance doesn't contain this child:", child
			return
		
		while child in self.children:
			self.children.remove(child)
			
	def attach_child(self,child):
		if child in self.children:
			print "error, can't attach the same child twice:", child
			return
			
		self.children.append(child)
		
	def isVisible(self,child):
		if child != None:
			return child.isVisible()
		else: return False
	
	def show(self):
		for child in self.children:
			if child.isVisible() == False:
				child.show()
				
	def close_specific(self,child,inspector_too=True):
		for i,child_ in enumerate(self.children):
			if child == child_:
				self.children.pop(i) # need to double check that this is the correct behavior
				if child != None: 
					child.close()
#					widget.deleteLater() #TODO: see if this is still needed
				if inspector_too and child.inspector != None:
					inspector = child.get_inspector()
					inspector.close()
				return
			
		#print "couldn't close",child
		
	def execute(self, logid=None):
		self.exec_()
		print logid
		if logid: E2end(logid) # We need to log the end of the process, don't we....
		return sys.exit()
		
	def hide_specific(self,child,inspector_too=True):
		for child_ in self.children:
			if child == child_:
				child.hide()
				inspector = child.get_inspector()
				if inspector != None:
					inspector.hide()
				return
			
		print "couldn't hide",child
		

	def show_specific(self,child):
		for child_ in self.children:
			if child == child_:
#				print "show",child
				if child.isVisible() == False:
					child.show()
					child.setFocus()
				child.raise_()
				return
	
		# if we make it here than we automatically attach the child
		self.attach_child(child)
		if child.isVisible() == False:
			child.show()
			child.setFocus()

	def start_timer(self,interval,function):
		print "START APP TIMER"
	
		if self.tmr != None:
			print "can't start a timer, already have one running. Call stop_timer first"
			#FIXME, add support for mutliple timers
			return
	
		self.tmr=QtCore.QTimer()
		self.tmr.setInterval(interval)
		QtCore.QObject.connect(self.tmr,QtCore.SIGNAL("timeout()"), function)
		self.tmr.start()
		
		self.timer_function = function
		
	
	def stop_timer(self):
		print "STOP APP TIMER"
		if self.tmr != None:
			QtCore.QObject.disconnect(self.tmr, QtCore.SIGNAL("timeout()"), self.timer_function)
			self.tmr = None
			self.timer_function = None
		else:
			print "warning, can't stop a timer when there is none"

		
	
class EMProgressDialog(QtGui.QProgressDialog):
	def __init__(self,label_text,cancel_button_text, minimum, maximum, parent = None):
		QtGui.QProgressDialog.__init__(self,label_text,cancel_button_text, minimum, maximum, parent)
		self.setWindowIcon(QtGui.QIcon(get_image_directory() + "/eman.png"))


def error(msg,title="Almost"):
	EMErrorMessageDisplay.run(msg,title)
	
class EMErrorMessageDisplay:
	'''
	Has a static error display function which is very useful
	'''
	def __init__(self): pass
	def run(error_message,title="Almost"):
		'''
		error_message is a list of error messages
		'''
		msg = QtGui.QMessageBox()
		msg.setWindowTitle(title)
		msg.setWindowIcon(QtGui.QIcon(get_image_directory() + "/eman.png"))
		mes = ""
		if isinstance(error_message,tuple): error_message=list(error_message)
		if isinstance(error_message,list):
			for error in error_message:
				mes += error
				
				if len(error) > 0 and error[-1] != '.':
					# correct my own inconsistencies....AWESOME
					mes += '.'
				if error != error_message[-1]: mes += "\n"
		else:
			mes = error_message
			if len(error_message) > 0 and error_message[-1] != '.':
				# correct my own inconsistencies....AWESOME
				mes += '.'
		msg.setText(mes)
		msg.exec_()
	
	run = staticmethod(run)
			
