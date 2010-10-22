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
from EMAN2 import remove_directories_from_name, get_image_directory
import EMAN2db
import weakref
import warnings

try: from PyQt4 import QtWebKit
except: pass

class ModuleEventsManager:
	'''
	Coordinates events of the various modules.
	To begin with this is only the close event, then I added the idle event
	'''
	def __init__(self,target,module):
		warnings.warn("Developers: ModuleEventsManager class should probably be deprecated and eventually removed.", DeprecationWarning)
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
	FTGL = "ftgl"
	GLUT = "glut"
	def hide(self):
		if self.qt_parent:
			self.qt_parent.hide()
	def resize(self, w, h):
		if self.qt_parent:
			QtOpenGL.QGLWidget.resize(self, w, h)
			self.qt_parent.resize(w, h)
	def show(self):
		if self.qt_parent:
			self.qt_parent.show()
	def setWindowTitle(self, title):
		if self.qt_parent:
			self.qt_parent.setWindowTitle(title)
	def __init__(self, parent=None,enable_timer=False, application_control=True,winid=None):
		QtOpenGL.QGLWidget.__init__(self,parent)
		self.qt_parent = EMParentWin(self, enable_timer)
		
		self.core_object =  QtCore.QObject()
		self.under_qt_control = True #TODO: figure out which value works better -- eventually eliminate this Qt control vs GL control was from 3D desktop days
		self.suppress_inspector = False # turn on to suppress showing the inspector
		self.inspector = None # this should be a qt widget, otherwise referred to as an inspector in eman
		self.winid=winid # a 'unique' identifier for the window used to restore locations on the screen
		
		self.image_change_count =  0# this is important when the user has more than one display instance of the same image, for instance in e2.py if 
		app = get_application()
		if app != None and application_control:
			app.attach_child(self)
		
		self.application_control = application_control
		self.file_name = ""
		self.disable_inspector = False
	
	def closeEvent(self, event):
		QtOpenGL.QGLWidget.closeEvent(self, event)
		self.qt_parent.close()
#		print 'signal: "module_closed"'
		self.emit(QtCore.SIGNAL("module_closed")) # this could be a useful signal, especially for something like the selector module, which can potentially show a lot of images but might want to close them all when it is closed
		
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

	def enable_inspector(self,val=True): self.disable_inspector = not val
	def load_font_renderer(self):
		try:
			self.font_render_mode = EMGLWidget.FTGL
			self.font_renderer = get_3d_font_renderer()
			self.font_renderer.set_face_size(16)
			self.font_renderer.set_font_mode(FTGLFontMode.TEXTURE)
		except:
			self.font_render_mode = EMGLWidget.GLUT
	def show_inspector(self,force=0):
		if self.disable_inspector: return
		self.emit(QtCore.SIGNAL("inspector_shown")) # debug only
		app = get_application()
		if app == None:
			print "can't show an inspector with having an associated application"
		
		if self.suppress_inspector: return
		if not force and self.inspector==None : return
		if not self.inspector : 
			self.inspector = self.get_inspector()
			if self.inspector == None: return # sometimes this happens
		if not app.child_is_attached(self.inspector):
			app.attach_child(self.inspector)
		app.show_specific(self.inspector)
	def update_inspector_texture(self):
		if self.inspector != None:
			self.inspector.update()

class EMGUIModule:
	FTGL = "ftgl"
	GLUT = "glut"
	
	def __init__(self,ensure_gl_context=False,application_control=True,winid=None):
		warnings.warn("EMGUIModule.__init__()--use EMGLWidget instead", DeprecationWarning)
		self.core_object =  QtCore.QObject()
		self.under_qt_control = False
		self.em_qt_inspector_widget = None # should be = EMQtWidgetModule(application) somewhere 
		self.suppress_inspector = False # turn on to suppress showing the inspector
		self.inspector = None # this should be a qt widget, otherwise referred to as an inspector in eman
		self.winid=None # a 'unique' identifier for the window used to restore locations on the screen
		
		self.parent = None # EMParentWin
		
		self.gl_parent = None # EMImage3DWidget, EMImage2DWidget, EMImageMXWidget etc
		
		self.gl_widget = None # the thing that nows the dimensions of the object (width,height), could be EMImage3DWidget or EM3DGLVolume etc
		
		self.gl_context_parent = None # should be the thing that owns the associated context
		self.qt_context_parent = None # should the the thing that will emit signals
		
		self.image_change_count =  0# this is important when the user has more than one display instance of the same image, for instance in e2.py if 
		self.application_control = application_control
		self.file_name = ""
		app = get_application()
		if app != None and application_control:
			app.attach_child(self)

		if ensure_gl_context and app != None:
			app.ensure_gl_context(self)
		
		self.emit_events = False
		self.disable_inspector = False
		
		self.winid=winid

		
	def enable_emit_events(self,val=True):
		#print "set emit events to",val
		self.emit_events = val
	
	def is_emitting(self): return self.emit_events
	
	def get_emit_signals_and_connections(self): return {}
	
	
	
	def setWindowTitle(self,title):
		pass
		
	def enable_inspector(self,val=True): self.disable_inspector = not val
	def emitter(self):
		return self.core_object
	
	def show(self):
		try:
			vis=self.qt_context_parent.isVisible()
		except : vis=False
		get_application().show_specific(self)
		
		if self.under_qt_control and self.winid!=None and not vis  :
			hdb=EMAN2db.EMAN2DB.open_db()
			hdb.open_dict("window_placement")
			db=hdb.window_placement

			if db[self.winid]!=None :
				self.qt_context_parent.restoreGeometry(db[self.winid])
	
	def emit(self,*args,**kargs):
		self.core_object.emit(*args,**kargs)
		
	def __del__(self):
		self.core_object.deleteLater()
		
	def updateGL(self):
		if self.gl_widget != None and self.under_qt_control:
			self.gl_widget.updateGL()
	
	def is_visible(self):
		return self.qt_context_parent.isVisible()
	
	def set_gl_context_parent(self,parent): self.gl_context_parent = parent
	def get_gl_context_parent(self): return self.gl_context_parent
	
	def set_qt_context_parent(self,parent): self.qt_context_parent = parent
	def get_qt_context_parent(self): return self.qt_context_parent
	
	def set_parent(self,parent): self.parent = parent
	def get_parent(self): return self.parent()
	
	def set_gl_parent(self,parent): self.gl_parent = parent
	def get_gl_parent(self): return self.gl_parent
	
	def set_gl_widget(self,gl_widget): self.gl_widget = gl_widget
	def get_gl_widget(self): return self.gl_widget
	
	def get_inspector(self): raise # this need to be supplied
	
	def get_last_render_image_display_count(self): return self.image_change_count
	
	def get_em_inspector(self):
		return self.em_qt_inspector_widget
	
	def show_inspector(self,force=0):
		if self.disable_inspector: return
		self.emit(QtCore.SIGNAL("inspector_shown")) # debug only
		app = get_application()
		if app == None:
			print "can't show an inspector with having an associated application"
		
		if self.suppress_inspector: return
		if not force and self.inspector==None : return
		if not self.inspector : 
			self.inspector = self.get_inspector()
			if self.inspector == None: return # sometimes this happens
		if not self.em_qt_inspector_widget:
			self.em_qt_inspector_widget = EMQtWidgetModule(self.inspector)
			self.em_qt_inspector_widget.setWindowTitle(remove_directories_from_name(self.file_name))
			if self.gl_widget != None:
				try:
					self.em_qt_inspector_widget.set_selected(self.gl_widget.decoration.is_selected())
				except: pass

		if not app.child_is_attached(self.em_qt_inspector_widget):
			app.attach_child(self.em_qt_inspector_widget)
			
		app.show_specific(self.em_qt_inspector_widget)
	
	def update_inspector_texture(self):
		if self.em_qt_inspector_widget != None:
			self.em_qt_inspector_widget.force_texture_update()
	
	def set_winid(self,winid):
		"""This is a unique name for this window for window positioning purposes"""
		self.winid=winid

		hdb=EMAN2db.EMAN2DB.open_db()
		hdb.open_dict("window_placement")
		db=hdb.window_placement

		if db[self.winid]!=None :
			self.qt_context_parent.restoreGeometry(db[self.winid])
		else :
			db[self.winid]=self.qt_context_parent.saveGeometry()

	
	def closeEvent(self,event):
		if self.under_qt_control and self.winid!=None and self.qt_context_parent.isVisible():
			hdb=EMAN2db.EMAN2DB.open_db()
			hdb.open_dict("window_placement")
			db=hdb.window_placement

			db[self.winid]=self.qt_context_parent.saveGeometry()

		
		if get_application() != None:
			if self.em_qt_inspector_widget != None: 
				self.em_qt_inspector_widget.closeEvent(event)
			get_application().close_specific(self)
		self.emit(QtCore.SIGNAL("module_closed")) # this could be a useful signal, especially for something like the selector module, which can potentially show a lot of images but might want to close them all when it is closed
		
	def load_font_renderer(self):
		try:
			self.font_render_mode = EMGUIModule.FTGL
			self.font_renderer = get_3d_font_renderer()
			self.font_renderer.set_face_size(16)
			self.font_renderer.set_font_mode(FTGLFontMode.TEXTURE)
		except:
			self.font_render_mode = EMGUIModule.GLUT
			
	def mouseDoubleClickEvent(self,event):
		pass
	
#	def emit(self,*args,**kargs):
#		self.application.get_qt_emitter(self).emit(*args,**kargs)
		
	def make_connections(self,calling_object):
		'''
		Should overwrite this function in the subclass 
		'''
		pass
	
	def set_inspector_selected(self,bool):
		if self.em_qt_inspector_widget != None: self.em_qt_inspector_widget.set_selected(bool)

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
			#self.browser2 = QtGui.QTextBrowser()
			##url = QtCore.QUrl("http://blake.bcm.edu/emanwiki/e2display")
			#url = QtCore.QUrl("http://www.google.com")
			#url.setPort(80)
			##print url.port()
			#self.browser2.setSource(url)
			##print browser2.port()
			#self.browser2.show()
			#self.browser2.resize(800,800)

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

class EMApplication:
	def __init__(self,qt_application_control=True):
		warnings.warn("EMApplication.__init__()--use EMApp instead", DeprecationWarning)
		self.children = []
		if qt_application_control:
			self.app = QtGui.QApplication(sys.argv)
		else: self.app = None
		
		self.qt_emission_registry = {}
		
		if em_app_instance.get_instance() == None:
			em_app_instance.set_instance(self)
	def child_is_attached(self,query_child):
		if query_child in self.children: return True
		else: return False
	
	def get_app(self): return self.app
	
	def set_app(self,app): self.app = app # hack
	
	def attach_child(self,child):
		raise
	
	def detach_child(self,child):
		raise
	
	def show(self):
		raise
	
	def isVisible(self,child): raise
	
	def show_specific(self,object):
		raise
	
	def hide_specific(self,child,inspector_too=True):
		raise
		
	def close_specific(self,child,inspector_too=True):
		raise
	
	def execute(self):
		if self.app != None:
			return sys.exit(self.app.exec_())
		else: return
		
	def quit(self):
		if self.app != None:
			self.app.quit()
	
	def ensure_gl_context(self,child): raise
			
	def setOverrideCursor(self,cursor_type):
		if self.app != None:
			self.app.setOverrideCursor(cursor_type)
			
	def exec_(self):
		self.app.exec_()
#	def __getattr__( self, name ):
#		try: return getattr(self,name)
#		except:	return getattr( self.app, name )

	def get_qt_emitter(self,child):
		raise
	
	def get_qt_gl_updategl_target(self,child):
		raise

	def register_qt_emitter(self,child,emitter):
	
		self.qt_emission_registry[child] = emitter
	
	
	def deregister_qt_emitter(self,child):
		self.qt_emission_registry.pop(child)

	#def starting_up(self):
		#if self.app != None:
			#if QtGui.QApplication.startingUp():
			#a = QtGui.QApplication.QApplication(sys.argv)

	def processEvents(self):
		self.app.processEvents()
	

class EMStandAloneApplication(EMApplication):
	def __init__(self,qt_application_control=True):
		warnings.warn("EMStandAloneApplication.__init__()--use EMApp instead", DeprecationWarning)
		
		# Stuff for display synchronization in e2.py
		self.timer_function = None
		self.tmr = None
		
		EMApplication.__init__(self,qt_application_control)
	
#	def __del__(self):
#		print "stand alone application death"
		
	def detach_child(self,child):
		for i,child_ in enumerate(self.children):
			if child_ == child:
				self.children.pop(i)
				return
	
		print "error, can't detach a child that doesn't belong to this",child,child.get_child()
	
	def attach_child(self,child):
		for i in self.children:
			if i == child:
				print "error, can't attach the same child twice",child
				return
			
		self.children.append(child)
		
	def ensure_gl_context(self,child):
		try: child.get_qt_widget().initGL()
		except: child.get_qt_widget().glInit()
		
	def processEvents(self):
		self.app.processEvents()
	
	def isVisible(self,child):
		if child.gl_widget != None:
			return child.gl_widget.isVisible()
		else: return False
	
	def show(self):
		for child in self.children:
#			print "showing all",child
			widget = child.get_qt_widget()
			if widget.isVisible() == False:
				widget.show()
				
	def close_specific(self,child,inspector_too=True):
		for i,child_ in enumerate(self.children):
			if child == child_:
				widget = child.get_qt_widget()
				self.children.pop(i) # need to double check that this is the correct behavior
#				print "popped",widget
				if widget != None: 
					widget.close() # widget was already closed ?
					widget.deleteLater()
				if inspector_too and child.inspector != None:
					inspector = child.get_inspector()
					inspector.close()
				return
			
		#print "couldn't close",child
	
	def hide_specific(self,child,inspector_too=True):
		for child_ in self.children:
			if child == child_:
				widget = child.get_qt_widget()
				widget.hide()
				inspector = child.get_inspector()
				if inspector != None:
					inspector.hide()
				return
			
		print "couldn't hide",child
		
	def get_qt_emitter(self,child):
		if isinstance(child,QtCore.QObject):
			return child
	
		# WARNING, THIS FUNCTIONALITY COULD CURRENTLY BE BROKEN DUE TO THE LINE I ADDED ABOVE
		for child_ in self.children:
			if child == child_:
				return child.get_qt_widget()
		
		try:
			return self.qt_emission_registry[child]
		except:
			print "couldn't get emitter", child
			return None

	
	def get_qt_gl_updategl_target(self,child):
		return self.get_qt_emitter(child) # it's the same thing for now

	
	def show_specific(self,child):
		for child_ in self.children:
			if child == child_:
				widget = child.get_qt_widget()
#				print "show",child
				if widget.isVisible() == False:
					widget.show()
					widget.setFocus()
				widget.raise_()
				return
	
		# if we make it here than we automatically attach the child
		
		self.attach_child(child)
		widget = child.get_qt_widget()
		if widget.isVisible() == False:
			widget.show()
			widget.setFocus()
		

		
	def __call__( *args, **kwargs ):
		return QtGui.qApp

	def exec_loop( *args, **kwargs ):
		pass

	def start_timer(self,interval,function):
	
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
		if self.tmr != None:
			QtCore.QObject.disconnect(self.tmr, QtCore.SIGNAL("timeout()"), self.timer_function)
			self.tmr = None
			self.timer_function = None
		else:
			print "warning, can't stop a timer when there is none"



class EMApp(QtGui.QApplication):
	def __init__(self,qt_application_control=True):
		self.children = []
		
		# Stuff for display synchronization in e2.py
		self.timer_function = None
		self.tmr = None
		
		QtGui.QApplication.__init__(self, sys.argv)
		
		if em_app_instance.get_instance() == None:
			em_app_instance.set_instance(self)
	
#	def __del__(self):
#		print "stand alone application death"
	def child_is_attached(self,query_child):
		if query_child in self.children: return True
		else: return False
	def detach_child(self,child):
		for i,child_ in enumerate(self.children):
			if child_ == child:
				self.children.pop(i)
				return
	
		print "error, can't detach a child that doesn't belong to this",child,child.get_child()
	
	def attach_child(self,child):
		for i in self.children:
			if i == child:
				print "error, can't attach the same child twice",child
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
		
	def execute(self):
		return sys.exit(self.exec_())
		
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
		if self.tmr != None:
			QtCore.QObject.disconnect(self.tmr, QtCore.SIGNAL("timeout()"), self.timer_function)
			self.tmr = None
			self.timer_function = None
		else:
			print "warning, can't stop a timer when there is none"



class EMQtWidgetModule(EMGUIModule):
	def __init__(self,qt_widget,winid=None):
		warnings.warn("EMQtWidgetModule.__init__()--use EMGLWidget, QWidget, or QDialog instead", DeprecationWarning)
		self.qt_widget = qt_widget
		self.gl_widget = None
		EMGUIModule.__init__(self,winid=winid)
		
		self.selected = False
	
#	def __del__(self): 
#		print "qt widget module death"
#		import sys
#		print sys.getrefcount(self.qt_widget),"for self.widget"
#		
#		self.qt_widget.deleteLater()
		
	#def __del__(self):
		#print "in qt delete"
		#if self.qt_widget != None:
			#print "yo"
			#self.qt_widget.clear_gl_memory(self)
	
	def set_selected(self,bool):
		self.selected = bool
		if self.gl_widget != None:self.gl_widget.set_selected(self.selected)
			
	def get_qt_widget(self):
		self.under_qt_control = True
#		print "returning self.qt_widget"
		return self.qt_widget
	
	def get_desktop_hint(self):
		return self.qt_widget.get_desktop_hint()

	def keyPressEvent(self,event):
		self.qt_widget.keyPressEvent(event)

	def get_child(self):
		return self.qt_widget

	def lock_texture(self):
		self.gl_widget.lock_texture()
	
	def unlock_texture(self):
		self.gl_widget.unlock_texture()
		
	def width(self): return self.qt_widget.width()
	
	def height(self): return self.qt_widget.height()
	
	def closeEvent(self,event) :
		if self.qt_widget != None:
			self.qt_widget.close()
				
		app = get_application()
		if app != None:
			app.close_specific(self)
		
		self.gl_widget = None
		self.qt_widget = None
		self.emit(QtCore.SIGNAL("module_closed")) # this could be a useful signal, especially for something like the selector module, which can potentially show a lot of images but might want to close them all when it is closed
		
	def connect_qt_pop_up_application_event(self,signal):
		get_application().connect_qt_pop_up_application_event(signal,self)
	
	
	def force_texture_update(self,val=True):
		if not self.under_qt_control and val == True:
			if get_application().isVisible(self):
				#self.gl_widget.drawable.gen_texture = True
				self.gl_widget.drawable.set_refresh_dl()
				
	def get_inspector(self):
		'''

		'''
		return None
	
	def setWindowTitle(self,title):
		if self.qt_widget != None:
			self.qt_widget.setWindowTitle(title)
			
	
class EMProgressDialog(QtGui.QProgressDialog):
	def __init__(self,label_text,cancel_button_text, minimum, maximum, parent = None):
		QtGui.QProgressDialog.__init__(self,label_text,cancel_button_text, minimum, maximum, parent)
		self.setWindowIcon(QtGui.QIcon(get_image_directory() + "/eman.png"))
	def get_desktop_hint(self): #TODO: get_desktop_hint() functions should probably be removed everywhere--remnant of 3D desktop days --Ross
		return "dialog"
			