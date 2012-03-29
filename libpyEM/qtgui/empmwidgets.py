#!/usr/bin/env python
#
# Author: John Flanagan Oct 20th 2011 (jfflanag@bcm.edu)
# Copyright (c) 2000-2011 Baylor College of Medicine
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

# These are widgets that e2projectmanger instatiates to make a GUI interface for the e2 programs.There should be enough widgets to represent
# just about any e2program, but if you feel the desire to make a new one, just subclass PMBaseWidget, and implemnt getValue and setValue.
# You may also need to reimplemnt getArgument (which returns the argument used in calling the e2program), if the default will not work for you.
# In addition, you'll need to add a line in the class PMGUIWidget (e2projectmanager) to instatiate the widget based on the value of 'guitype' 

from EMAN2db import db_check_dict
import sys, math, weakref
from PyQt4 import QtCore, QtGui
from PyQt4.QtCore import Qt
from emselector import EMSelectorDialog	# This will be replaced by something more sensible in the future
import re, os, glob
from embrowser import EMBrowserWidget
from empmtabwidgets import *

class PMComboBox(QtGui.QComboBox):
	""" Reimplment the QComboBox to remove wheel widget activation """
	def __init__(self):
		QtGui.QComboBox.__init__(self)
		
	def wheelEvent(self, event):
		""" Ignore wheelevents is not poped up """
		event.ignore()

			
class PMBaseWidget(QtGui.QWidget):
	""" A base widget upon which all the other PM widgets are derived """
	def __init__(self, name, mode="",returnNone=False):
		QtGui.QWidget.__init__(self)
		
		self.postional = False
		self.name = name
		self.setMode(mode)
		self.errormessage = None
		self.noarg = False
		self.returnNone = returnNone
		
	def getValue(self):
		raise NotImplementedError("Sub class must reimplemnt 'getValue' function")
	
	def setValue(self, value, quiet=False):
		raise NotImplementedError("Sub class must reimplemnt 'setValue' function")
		 
	def getName(self):
		return self.name
		
	def setPositional(self, position):
		self.postional = position
		
	def getPositional(self):
		return self.postional
	
	def setMode(self, mode):
		self.mode = mode
	
	def getMode(self):
		return self.mode
		
	def getArgument(self):
		# If the value is None or blank then do not yeild an option. Obviously this will nver happen for Int bool or float
		if str(self.getValue()) == "" or self.noarg:
			return ""
		elif  str(self.getValue()).upper() == "NONE" and not self.returnNone:
			return ""
		else:
			""" There are two sorts of arguments: Posional and optional """
			if self.getPositional():
				return str(self.getValue())
			else:
				return "--%s=%s"%(self.getName(), str(self.getValue()))
	
	def setErrorMessage(self, message):
		self.errormessage = message
		
	def getErrorMessage(self):
		return self.errormessage
	
class PMIntEntryWidget(PMBaseWidget):
	""" A Widget for geting Int values. Type and range is checked """
		
	@staticmethod 
	def copyWidget(widget):
		""" Basically a copy constructor to get around QT and python limitations """
		return PMIntEntryWidget(widget.getName(), widget.getValue(), widget.getMode(), widget.lrange, widget.urange, widget.postional, widget.initdefault)
		
	def __init__(self, name, value, mode, lrange=None, urange=None, postional=False, initdefault=None):
		PMBaseWidget.__init__(self, name, mode) 
		gridbox = QtGui.QGridLayout()
		label = QtGui.QLabel(name)
		self.intbox = QtGui.QLineEdit()
		gridbox.addWidget(label, 0, 0)
		gridbox.addWidget(self.intbox, 0, 1)
		self.setLayout(gridbox)

		self.value = None
		self.lrange = lrange
		self.urange = urange
		self.initdefault = initdefault
		self.setValue(value)
		self.setPositional(postional)

		QtCore.QObject.connect(self.intbox,QtCore.SIGNAL("editingFinished()"),self._on_intchanged)
		
	def _on_intchanged(self, quiet=False):
		if str(self.intbox.text()).upper() == "NONE":
			self.value = None
			self.setErrorMessage(None)
			return
		try:
			self.value = int(self.intbox.text())
			self.setErrorMessage(None)
			self._confirm_bounds()
		except ValueError:
			self.intbox.setText("") 
			self.setErrorMessage("Invalid type, Int neeeded in %s"%self.getName())
			if self.isVisible() and not quiet: self.emit(QtCore.SIGNAL("pmmessage(QString)"),"Invalid type, Int neeeded in %s"%self.getName())
			
	def _confirm_bounds(self):
		if self.lrange != None and (self.value < self.lrange):
			self.intbox.setText(str(self.lrange))
			if self.isVisible(): self.emit(QtCore.SIGNAL("pmmessage(QString)"),"Value too low for '%s', clipping to '%d'"%(self.name,self.lrange))
		if self.urange != None and (self.value > self.urange):
			self.intbox.setText(str(self.urange))
			if self.isVisible(): self.emit(QtCore.SIGNAL("pmmessage(QString)"),"Value too high for '%s', clipping to '%d'"%(self.name,self.urange))
			
	def getValue(self):
		return self.value
		
	def setValue(self, value, quiet=False):
		self.intbox.setText(str(value))
		self._on_intchanged(quiet=quiet)
		
	def setEnabled(self, state):
		self.intbox.setEnabled(state)

class PMShrinkEntryWidget(PMIntEntryWidget):
	""" A widget for shink options. If this entry is set to <= 1 then no argument is returned """
	
	@staticmethod 
	def copyWidget(widget):
		""" Basically a copy constructor to get around QT and python limitations """
		return PMShrinkEntryWidget(widget.getName(), widget.getValue(), widget.getMode(), widget.lrange, widget.postional, widget.initdefault)
		
	def __init__(self, name, value, mode, lrange=None, postional=False, initdefault=None):
		PMIntEntryWidget.__init__(self, name, value, mode, lrange=lrange, urange=None, postional=postional, initdefault=initdefault)
		
	def _on_intchanged(self, quiet=False):
		if str(self.intbox.text()).upper() == "NONE" or str(self.intbox.text()) == "":
			self.value = self.lrange - 1
			self.intbox.setText(str(self.lrange-1))
			self.setErrorMessage(None)
			return
		try:
			self.value = int(self.intbox.text())
			self.setErrorMessage(None)
			self._confirm_bounds()
		except ValueError:
			self.value = self.lrange - 1
			self.intbox.setText(str(self.lrange-1)) 
			self.setErrorMessage("Invalid type, Int neeeded in %s"%self.getName())
			if self.isVisible() and not quiet: self.emit(QtCore.SIGNAL("pmmessage(QString)"),"Invalid type, Int neeeded in %s"%self.getName())
		
	def _confirm_bounds(self):
		if self.lrange != None and (self.value < self.lrange):
			self.value = self.lrange - 1
			self.intbox.setText(str(self.lrange-1))
			self.noarg = True
			return
		self.noarg = False
		
class PMFloatEntryWidget(PMBaseWidget):
	""" A Widget for geting Float values. Type and range is checked """
	
	@staticmethod 
	def copyWidget(widget):
		""" Basically a copy constructor to get around QT and python limitations """
		return PMFloatEntryWidget(widget.getName(), widget.getValue(), widget.getMode(), widget.lrange, widget.urange, widget.postional, widget.initdefault)
		
	def __init__(self, name, value, mode, lrange=None, urange=None, postional=False, initdefault=None):
		PMBaseWidget.__init__(self, name, mode) 
		gridbox = QtGui.QGridLayout()
		label = QtGui.QLabel(name)
		self.floatbox = QtGui.QLineEdit()
		gridbox.addWidget(label, 0, 0)
		gridbox.addWidget(self.floatbox, 0, 1)
		self.setLayout(gridbox)
		
		self.value = None
		self.lrange = lrange
		self.urange = urange
		self.initdefault = initdefault
		self.setValue(value)
		self.setPositional(postional)
		
		QtCore.QObject.connect(self.floatbox,QtCore.SIGNAL("editingFinished()"),self._on_floatchanged)
		
	def _on_floatchanged(self, quiet=False):
		if str(self.floatbox.text()).upper() == "NONE":
			self.value = None
			self.setErrorMessage(None)
			return
		try:
			self.value = float(self.floatbox.text())
			self.setErrorMessage(None)
			self._confirm_bounds()
		except ValueError:
			self.floatbox.setText("") 
			self.setErrorMessage("Invalid type, float needed in '%s'"%self.getName())
			if self.isVisible() and not quiet: self.emit(QtCore.SIGNAL("pmmessage(QString)"),"Invalid type, float needed in '%s'"%self.getName())
			
	def _confirm_bounds(self):
		if self.lrange and (self.value < self.lrange):
			self.floatbox.setText(str(self.lrange))
			if self.isVisible(): self.emit(QtCore.SIGNAL("pmmessage(QString)"),"Value too low for '%s', clipping to '%f'"%(self.name,self.lrange))
		if self.urange and (self.value > self.urange):
			self.floatbox.setText(str(self.urange))
			if self.isVisible(): self.emit(QtCore.SIGNAL("pmmessage(QString)"),"Value too high for '%s', clipping to '%f'"%(self.name,self.urange))
			
	def getValue(self):
		return self.value
		
	def setValue(self, value, quiet=False):
		self.floatbox.setText(str(value))
		self._on_floatchanged(quiet=quiet)
	
	def setEnabled(self, state):
		self.floatbox.setEnabled(state)
		
class PMStringEntryWidget(PMBaseWidget):
	""" A Widget for geting String values. Type is checked """
	
	@staticmethod 
	def copyWidget(widget):
		""" Basically a copy constructor to get around QT and python limitations """
		return PMStringEntryWidget(widget.getName(), widget.getValue(), widget.getMode(), widget.postional, widget.initdefault, widget.returnNone)
		
	def __init__(self, name, string, mode, postional=False, initdefault=None, returnNone=False):
		PMBaseWidget.__init__(self, name, mode, returnNone=returnNone) 
		gridbox = QtGui.QGridLayout()
		label = QtGui.QLabel(name)
		self.stringbox = QtGui.QLineEdit()
		gridbox.addWidget(label, 0, 0)
		gridbox.addWidget(self.stringbox, 0, 1)
		self.setLayout(gridbox)
		
		self.initdefault = initdefault
		self.setValue(string)
		self.setPositional(postional)
		
		QtCore.QObject.connect(self.stringbox,QtCore.SIGNAL("editingFinished()"),self._on_stringchanged)
	
	def _on_stringchanged(self):
		self.string = str(self.stringbox.text())
		
	def getValue(self):
		# What to do with None tpye values? For strings, just set None to "". This should be equivilent
		return self.string
		
	def setValue(self, string, quiet=False):
		if str(string).upper() == "NONE" and not self.returnNone: string = ""	# If none is passed set to blank
		self.stringbox.setText(str(string))
		self._on_stringchanged()
		self.setErrorMessage(None)
	
	def setEnabled(self, state):
		self.stringbox.setEnabled(state)
		
class PMHeaderWidget(PMBaseWidget):
	""" A widget to add a header """
	def __init__(self, name, header):
		PMBaseWidget.__init__(self, name)
		gridbox = QtGui.QGridLayout()
		self.header = QtGui.QLabel()
		font = QtGui.QFont()
		font.setBold(True)
		self.header.setFont(font)
		gridbox.addWidget(self.header)
		self.setLayout(gridbox)

		self.setValue(header)
		
	def getValue(self):
		return str(self.header.text())
		
	def setValue(self, header):
		self.header.setText(header)
		self.setErrorMessage(None)
		
	def getArgument(self):
		""" Obvioulsy the hear does give an argument """
		return None

class PMBoolWidget(PMBaseWidget):
	""" A Widget for getting Bool values. Type is checked """
	
	@staticmethod 
	def copyWidget(widget):
		""" Basically a copy constructor to get around QT and python limitations """
		return PMBoolWidget(widget.getName(), widget.getValue(), widget.getMode(), widget.initdefault)
		
	def __init__(self, name, boolvalue, mode, initdefault=None):
		PMBaseWidget.__init__(self, name, mode)
		gridbox = QtGui.QGridLayout()
		self.boolbox = QtGui.QCheckBox(name)
		gridbox.addWidget(self.boolbox, 0, 0)
		self.setLayout(gridbox)
		

		self.boolvalue = boolvalue
		self.initdefault = initdefault
		self.setValue(self.boolvalue)
		
		QtCore.QObject.connect(self.boolbox,QtCore.SIGNAL("stateChanged(int)"),self._on_boolchanged)
	
	def _on_boolchanged(self):
		self.boolvalue = self.boolbox.isChecked()
		
	def getValue(self):
		return self.boolvalue
		
	def setValue(self, boolvalue, quiet=False):
		self.boolbox.setChecked(boolvalue)
		self. _on_boolchanged()
		self.setErrorMessage(None)
		
	def getArgument(self):
		""" Only return argument if set to true"""
		if self.getValue():
			return "--%s"%(self.getName())
		else:
			return ""
		
class PMFileNameWidget(PMBaseWidget):
	""" A Widget for geting filenames. Type is checked """
	@staticmethod 
	def copyWidget(widget):
		""" Basically a copy constructor to get around QT and python limitations """
		return PMFileNameWidget(widget.getName(), widget.filename, widget.mode, widget.browser, widget.postional, widget.initdefault, widget.checkfileexist, infolabels=False) 
		
	def __init__(self, name, filename, mode, browser, postional=True, initdefault=None, checkfileexist=True, infolabels=True):
		PMBaseWidget.__init__(self, name, mode) 
		gridbox = QtGui.QGridLayout()
		label = QtGui.QLabel(name)
		self.browser = browser
		self.filenamebox = QtGui.QLineEdit()
		self.browsebutton = QtGui.QPushButton("Browse")
		self.infolabel = QtGui.QLabel("Num Images: None")
		gridbox.addWidget(label, 0, 0)
		gridbox.addWidget(self.filenamebox, 0, 1)
		gridbox.addWidget(self.browsebutton, 0, 2)
		if infolabels: gridbox.addWidget(self.infolabel, 1, 1, 1, 2)
		self.setLayout(gridbox)

		self.initdefault = initdefault
		self.checkfileexist= checkfileexist
		self.setValue(filename)
		self.setPositional(postional)
		
		QtCore.QObject.connect(self.filenamebox,QtCore.SIGNAL("editingFinished()"),self._on_filenamechanged)
		QtCore.QObject.connect(self.browsebutton,QtCore.SIGNAL("clicked()"),self._on_clicked)
	
	def _on_filenamechanged(self):
		self.setValue(str(self.filenamebox.text()))
		
	def _on_cancel(self):
		self.window.close()
		
	def _on_ok(self):
		filename = ""
		for f in self.window.getResult():
			filename += (" "+f)
		self.setValue(filename[1:])
		self.window.close()
		
	def _on_clicked(self):
		self.window = eval(self.browser)
		QtCore.QObject.connect(self.window, QtCore.SIGNAL("ok"),self._on_ok)
		QtCore.QObject.connect(self.window, QtCore.SIGNAL("cancel"),self._on_cancel)
		self.window.show()
		
	def getValue(self):
		return self.filename
		
	def setValue(self, filename, quiet=False):
		# Check to see if the file exists
		filename = str(filename)
		if filename.upper() == "NONE": filename=""	# If none is passed set to blank
		# In some cases a file is optional
		if self.checkfileexist:
			# We need to check if the field is blank
			if filename == "": 
				self._onBadFile(filename, quiet)
				return
			# In not blank then check to ensure each file is 'ok'. Not that a list of files is accepted
			if not self._checkfiles(filename): return
		
		# If all is well, then  we are happy
		self.filename = filename
		self.filenamebox.setText(filename)
		self.setErrorMessage(None)
		self.emit(QtCore.SIGNAL("pmfilename(QString)"),self.getValue())
	
	def _checkfiles(self, filename):
		files = filename.split()
		numimages = 0
		for f in files:
			if not os.access(f, os.F_OK) and not db_check_dict(f):
				self._onBadFile(f)
				# Display the rubbish file to the user
				self.filenamebox.setText(filename)
				return False
			numimages += EMUtil.get_image_count(f)
		self.infolabel.setText("Num Images: %d"%numimages)
		return True
		
	def _onBadFile(self, filename, quiet=False):
		self.filename = None
		self.setErrorMessage("File '%s' from field '%s' does not exist"%(filename,self.getName()))
		if self.isVisible() and not quiet: self.emit(QtCore.SIGNAL("pmmessage(QString)"),"File '%s' from field '%s' does not exist"%(filename,self.getName()))

class PMDirectoryWidget(PMBaseWidget):
	""" A Widget for display dircories of a certian type """
	
	@staticmethod 
	def copyWidget(widget):
		""" Basically a copy constructor to get around QT and python limitations """
		return PMDirectoryWidget(widget.getName(), widget.dirbasename, widget.getValue(), widget.getMode(), widget.postional, widget.initdefault)
		
	def __init__(self, name, dirbasename, default, mode, postional=False, initdefault=None):
		PMBaseWidget.__init__(self, name, mode) 
		self.dirbasename = dirbasename
		gridbox = QtGui.QGridLayout()
		label = QtGui.QLabel(name)
		self.combobox = PMComboBox()
		gridbox.addWidget(label, 0, 0)
		gridbox.addWidget(self.combobox, 0, 1)
		self.setLayout(gridbox)
		
		self.initdefault = initdefault
		self.setValue(default)
		self.setPositional(postional)
		
		self.connect(self.combobox, QtCore.SIGNAL("activated(const QString &)"), self.setValue)
				
	def updateDirs(self):
		for idx in xrange(self.combobox.count()):
			self.combobox.removeItem(self.combobox.count()-1)
		# This extra code allows use to have more than one type of directory
		patterns = self.dirbasename.split("|")
		dirs = []
		for pattern in patterns:
			dirs.extend(glob.glob("%s*"%pattern))
		for directory in sorted(dirs):
			self.combobox.addItem(str(directory))
			
	def getValue(self):
		return str(self.combobox.currentText())
		
	def setValue(self, value, quiet=False):
		self.updateDirs()
		idx = self.combobox.findText(str(value))
		if idx > -1:
			self.combobox.setCurrentIndex(idx)
		self.setErrorMessage(None)
		
class PMComboWidget(PMBaseWidget):
	""" A Widget for combo boxes. Type is checked """
	
	@staticmethod 
	def copyWidget(widget):
		""" Basically a copy constructor to get around QT and python limitations """
		return PMComboWidget(widget.getName(), widget.choices, widget.getValue(), widget.getMode(), widget.postional, widget.datatype, widget.initdefault, widget.returnNone)	

	def __init__(self, name, choices, default, mode, postional=False,  datatype=str, initdefault=None, returnNone=False):
		PMBaseWidget.__init__(self, name, mode, returnNone=returnNone) 
		gridbox = QtGui.QGridLayout()
		label = QtGui.QLabel(name)
		self.combobox = PMComboBox()
		gridbox.addWidget(label, 0, 0)
		gridbox.addWidget(self.combobox, 0, 1)
		self.setLayout(gridbox)

		self.initdefault = initdefault
		self.datatype=datatype	# Must be int, float or str
		
		self.choices = sorted(choices)
		for choice in self.choices:
			self.combobox.addItem(str(choice))
		self.setValue(default)
		self.setPositional(postional)
		
		self.connect(self.combobox, QtCore.SIGNAL("activated(const QString &)"), self.setValue)
				
	def getValue(self):
		if str(self.combobox.currentText()) == "None":
			# In the e2 programs we actually need to specify None otherwise default behaviour will be implmented
			return "None"
		return self.datatype(self.combobox.currentText())
		
	def setValue(self, value, quiet=False):
		if value == '': value = "None"
		idx = self.combobox.findText(str(value))
		if idx > -1:
			self.combobox.setCurrentIndex(idx)
			self.setErrorMessage(None)
			return
		else:
			self.setErrorMessage("Value '%s' not found in combobox '%s'"%(value,self.getName()))
			if not quiet: self.emit(QtCore.SIGNAL("pmmessage(QString)"),"Value '%s' not found in combobox '%s'"%(value,self.getName()))
			return
		
			
class PMComboParamsWidget(PMBaseWidget):
	""" A Widget for combo boxes. Type is checked. For the combobox with params the datatype is always str """
	
	@staticmethod 
	def copyWidget(widget):
		""" Basically a copy constructor to get around QT and python limitations """
		return PMComboParamsWidget(widget.getName(), widget.choices, widget.getValue(), widget.getMode(), widget.postional, widget.initdefault, widget.returnNone)	
		
	def __init__(self, name, choices, default, mode, postional=False, initdefault=None, returnNone=False):
		PMBaseWidget.__init__(self, name, mode, returnNone=returnNone) 
		gridbox = QtGui.QGridLayout()
		label = QtGui.QLabel(name)
		self.combobox = PMComboBox()
		plabel = QtGui.QLabel("params:")
		self.params = QtGui.QLineEdit()
		gridbox.addWidget(label, 0, 0)
		gridbox.addWidget(self.combobox, 0, 1)
		gridbox.addWidget(plabel, 0, 2)
		gridbox.addWidget(self.params, 0, 3)
		self.setLayout(gridbox)

		self.choices = sorted(choices)
		for choice in self.choices:
			self.combobox.addItem(str(choice))
		self.combobox.addItem('None')
		
		self.initdefault = initdefault
		self.setValue(default)
		self.setPositional(postional)
		
		self.connect(self.combobox, QtCore.SIGNAL("activated(const QString &)"), self.setValue)
				
	def getValue(self):
		""" Return the value. Concatenate if necessary """
		if str(self.combobox.currentText()) == "None":
			# In the e2 programs we actually need to specify None otherwise default behaviour will be implmented
			return "None"
		if self.params.text() == "":
			return str(self.combobox.currentText())
		else:
			return str(self.combobox.currentText())+":"+self.params.text()
		
	def setValue(self, value, quiet=False):
		# First parse the value, as it may contain both a options and params
		if value == '': value = "None"
		values = self._parsedefault(str(value))
		# Next process the parsed value
		idx = self.combobox.findText(str(values[0]))
		if idx > -1:
			self.combobox.setCurrentIndex(idx)
		else:
			self.setErrorMessage("Value '%s' not found in combobox '%s'"%(values[0],self.getName()))
			if not quiet: self.emit(QtCore.SIGNAL("pmmessage(QString)"),"Value '%s' not found in combobox '%s'"%(values[0],self.getName()))
			return
		if len(values) == 2: self.params.setText(values[1])
		self.setErrorMessage(None)
			
	def _parsedefault(self, default):
		default=str(default)
		if default.find(":") != -1:
			return [default[:default.find(":")], default[default.find(":")+1:]]
		else:
			return [default]
		
class PMSymWidget(PMBaseWidget):
	""" A widget for getting/setting symmetry input """
	
	@staticmethod 
	def copyWidget(widget):
		""" Basically a copy constructor to get around QT and python limitations """
		return PMSymWidget(widget.getName(), widget.getValue(), widget.getMode(), widget.initdefault)	
		
	def __init__(self, name, default, mode, initdefault=None):
		PMBaseWidget.__init__(self, name, mode)
		gridbox = QtGui.QGridLayout()
		label = QtGui.QLabel(name)
		label.setAlignment(QtCore.Qt.AlignVCenter | QtCore.Qt.AlignRight)
		self.combobox = PMComboBox()
		self.symnumbox = PMIntEntryWidget("Symmetry Number", 0, mode, lrange=0)
		gridbox.addWidget(label, 0, 0)
		gridbox.addWidget(self.combobox, 0, 1)
		gridbox.addWidget(self.symnumbox, 0, 2, 1, 2)
		self.setLayout(gridbox)
		
		for i in ['icos','oct','tet','c','d','h']: self.combobox.addItem(i)
		
		self.connect(self.symnumbox,QtCore.SIGNAL("pmmessage(QString)"),self._on_message)
		
		self.initdefault = initdefault
		self.setValue(default)
	
	def _on_message(self, message):
		self.emit(QtCore.SIGNAL("pmmessage(QString)"),message)
		
	def getValue(self):
		""" Return the symmetry value """
		if str(self.combobox.currentText()) in ['icos','oct','tet']:
			return str(self.combobox.currentText())
		else:
			return str(self.combobox.currentText())+str(self.symnumbox.getValue())
		
	def setValue(self, value, quiet=False):
		""" Set the value. For example c1 is split into c and 1 and then set """
		if not value: return
		defsym = re.findall('[a-zA-Z]+', value)
		defsymnum = re.findall('[0-9]+', value)
		defsym = reduce(lambda x, y: x+y, defsym)
		# Deal with icos, tet, oct cases
		if defsymnum: 
			defsymnum = reduce(lambda x, y: x+y, defsymnum)
		else:
			defsymnum = 0
		
		idx = self.combobox.findText(defsym)
		if idx > -1:
			self.combobox.setCurrentIndex(idx)
		else:
			self.setErrorMessage("'%s' not a valid symmetry!!!"%value)
			if not quiet: self.emit(QtCore.SIGNAL("pmmessage(QString)"),"'%s' not a valid symmetry!!!"%value)
			return
		self.symnumbox.setValue(defsymnum)
		self.setErrorMessage(None)
		
	def getErrorMessage(self):
		if self.errormessage: return self.errormessage
		if self.symnumbox.getErrorMessage(): return self.symnumbox.getErrorMessage()

class PMMultiSymWidget(PMBaseWidget):
	""" A widget for getting/setting symmetry input from multiple models and this widget is runtime dynamic """
	
	@staticmethod 
	def copyWidget(widget):
		""" Basically a copy constructor to get around QT and python limitations """
		return PMMultiSymWidget(widget.getName(), widget.getMode(), widget.initdefault)
		
	def __init__(self, name, mode, initdefault=None):
		PMBaseWidget.__init__(self, name, mode)
		self.gridbox = QtGui.QVBoxLayout()
		self.gridbox.setContentsMargins(0,0,0,0)
		self.stackedwidget = QtGui.QStackedWidget()
		self.multisymwidgetlist = []
		self.gridbox.addWidget(self.stackedwidget)
		self.setLayout(self.gridbox)
		self.lastsymvalue = None
		self.initdefault = initdefault
	
	def setValue(self, value, quiet=False):
		if not value: return
		values = value.split(",")
		for i, widget in enumerate(self.multisymwidgetlist):
			if i > len(values) - 1: break	# Can add more defaults than exist!!!
			widget.setValue(values[i])
		self.lastsymvalue = value
		
	def getValue(self):
		multisym = ""
		for widget in self.multisymwidgetlist:
			multisym += ","+widget.getValue()
		multisym = multisym[1:]
		return multisym
		
	def update(self, files): 
		fileslist = str(files).split()
		amount = 0
		self.multisymwidgetlist = []
		# First remove the old widget
		if self.stackedwidget.currentIndex() == -1:
			widget = self.stackedwidget.widget(0)
			self.stackedwidget.removeWidget(widget)
			del(widget)
		# Then add new one
		multisym = QtGui.QWidget()
		vbox = QtGui.QVBoxLayout()
		for i,f in enumerate(fileslist):
			widget = PMSymWidget("Model%d"%i,"c1",self.getMode())
			self.multisymwidgetlist.append(widget)
			vbox.addWidget(widget)
			amount += 60	# A complete HACK
		if self.lastsymvalue: self.setValue(self.lastsymvalue)
		vbox.setContentsMargins(0,0,0,0)
		multisym.setLayout(vbox)
		self.stackedwidget.insertWidget(0, multisym)
		self.stackedwidget.setCurrentIndex(0)
		
		# Then resize the GUI widget, this was just hacked together
		self.setMinimumHeight(amount)
		self.setMaximumHeight(amount)
		
	def getErrorMessage(self):
		if self.errormessage: return self.errormessage
		for widget in self.multisymwidgetlist:
			if widget.getErrorMessage(): return widget.getErrorMessage()
			
class PMAutoMask3DWidget(PMBaseWidget):
	""" A Widget for getting automask 3D input """
	
	@staticmethod 
	def copyWidget(widget):
		""" Basically a copy constructor to get around QT and python limitations """
		return PMAutoMask3DWidget(widget.getName(), widget.getValue(), widget.getMode(), widget.initdefault)
		
	def __init__(self, name, default, mode, initdefault=None):
		PMBaseWidget.__init__(self, name, mode)
		gridbox = QtGui.QGridLayout()
		self.automask3dbool = QtGui.QCheckBox("Auto Mask 3D")
		self.params = []
		self.params.append(PMFloatEntryWidget("Threshold", 0.8, mode))
		self.params.append(PMIntEntryWidget("NMax", 30, mode))
		self.params.append(PMIntEntryWidget("Radius", 30, mode))
		self.params.append(PMIntEntryWidget("Mask Dilations", 5, mode))
		self.params.append(PMIntEntryWidget("Post Gaussian Dialations", 5, mode))
		gridbox.addWidget(self.automask3dbool, 0, 0)
		gridbox.addWidget(self.params[0], 1, 0)
		gridbox.addWidget(self.params[1], 1, 1)
		gridbox.addWidget(self.params[2], 1, 2)
		gridbox.addWidget(self.params[3], 2, 0)
		gridbox.addWidget(self.params[4], 2, 1, 1, 2)
		self.setLayout(gridbox)
		self.setValue(default)
		self.initdefault = initdefault
		
		QtCore.QObject.connect(self.automask3dbool,QtCore.SIGNAL("stateChanged(int)"),self._on_boolchanged)
		self.connect(self.params[0],QtCore.SIGNAL("pmmessage(QString)"),self._on_message)
		self.connect(self.params[1],QtCore.SIGNAL("pmmessage(QString)"),self._on_message)
		self.connect(self.params[2],QtCore.SIGNAL("pmmessage(QString)"),self._on_message)
		self.connect(self.params[3],QtCore.SIGNAL("pmmessage(QString)"),self._on_message)
		self.connect(self.params[4],QtCore.SIGNAL("pmmessage(QString)"),self._on_message)
		
	def _on_boolchanged(self):
		for widget in self.params:
			widget.setEnabled(self.automask3dbool.isChecked())
	
	def _on_message(self, message):
		self.emit(QtCore.SIGNAL("pmmessage(QString)"),message)
		
	def setValue(self, value, quiet=False):
		# if value is "" of None, set bool to false
		if not value:
			self.automask3dbool.setChecked(False)
			self._on_boolchanged()
			return
		# Otherwise parse input and set
		self.automask3dbool.setChecked(True)
		automaskval = value.split(',')
		for i, param in enumerate(self.params):
			param.setValue(automaskval[i])
			
	def getValue(self):
		if not self.automask3dbool.isChecked(): return ""
		value = ""
		# concatenate things
		for i in xrange(len(self.params)):
			value = value+","+str(self.params[i].getValue())
		value = value[1:]
		return value
		
	def getErrorMessage(self):
		if self.errormessage: return self.errormessage
		if self.params[0].getErrorMessage(): return self.params[0].getErrorMessage()
		if self.params[1].getErrorMessage(): return self.params[1].getErrorMessage()
		if self.params[2].getErrorMessage(): return self.params[2].getErrorMessage()
		if self.params[3].getErrorMessage(): return self.params[3].getErrorMessage()
		if self.params[4].getErrorMessage(): return self.params[4].getErrorMessage()
		