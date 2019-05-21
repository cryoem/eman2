#!/usr/bin/env python
from __future__ import print_function
from __future__ import absolute_import
from __future__ import division
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

# These are widgets that e2projectmanager instantiates to make a GUI interface for the e2 programs.There should be enough widgets to represent
# just about any e2program, but if you feel the desire to make a new one, just subclass PMBaseWidget, and implement getValue and setValue.
# You may also need to reimplement getArgument (which returns the argument used in calling the e2program), if the default will not work for you.
# In addition, you'll need to add a line in the class PMGUIWidget (e2projectmanager) to instantiate the widget based on the value of 'guitype'

from past.utils import old_div
from builtins import range
from EMAN2db import db_check_dict
import sys, math, weakref
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtCore import Qt
from .emselector import EMSelectorDialog	# This will be replaced by something more sensible in the future
import re, os, glob
from .embrowser import EMBrowserWidget
from .empmtabwidgets import *
from functools import reduce

class PMComboBox(QtWidgets.QComboBox):
	""" Reimplement the QComboBox to remove wheel widget activation """
	def __init__(self):
		QtWidgets.QComboBox.__init__(self)

	def wheelEvent(self, event):
		""" Ignore wheelevents is not poped up """
		event.ignore()


class PMBaseWidget(QtWidgets.QWidget):
	""" A base widget upon which all the other PM widgets are derived """
	pmmessage = QtCore.pyqtSignal(str)

	def __init__(self, name, mode="",returnNone=False):
		QtWidgets.QWidget.__init__(self)

		self.postional = False
		self.name = name
		self.setMode(mode)
		self.errormessage = None
		self.noarg = False
		self.returnNone = returnNone

	def getValue(self):
		""" Get the value for this widget """
		raise NotImplementedError("Sub class must reimplement 'getValue' function")

	def setValue(self, value, quiet=False):
		""" Set the value for this widget """
		raise NotImplementedError("Sub class must reimplement 'setValue' function")

	def getName(self):
		""" Return the name of this widget """
		return self.name

	def setPositional(self, position):
		""" Set whether or not this is a positional argument """
		self.postional = position

	def getPositional(self):
		""" Return true if this is a positional argument rather than an option """
		return self.postional

	def setMode(self, mode):
		""" Set the DB sharing mode """
		self.mode = mode

	def getMode(self):
		""" Get the DB sharing mode """
		return self.mode

	def getArgument(self):
		# If the value is None or blank then do not yield an option. Obviously this will never happen for int, bool or float
		if str(self.getValue()) == "" or self.noarg:
			return ""
		elif  str(self.getValue()).upper() == "NONE" and not self.returnNone:
			return ""
		else:
			""" There are two sorts of arguments: Positional and optional """
			if self.getPositional():
				return str(self.getValue())
			else:
				return "--%s=%s"%(self.getName(), str(self.getValue()))

	def setErrorMessage(self, message):
		""" Set the error message """
		self.errormessage = message

	def getErrorMessage(self):
		""" Return the error message if none then this returns blank """
		return self.errormessage

class PMIntEntryWidget(PMBaseWidget):
	""" A widget for getting int values. Type and range are checked """

	@staticmethod
	def copyWidget(widget):
		""" Basically a copy constructor to get around QT and python limitations """
		return PMIntEntryWidget(widget.getName(), widget.getValue(), widget.getMode(), widget.lrange, widget.urange, widget.postional, widget.initdefault)

	def __init__(self, name, value, mode, lrange=None, urange=None, postional=False, initdefault=None):
		PMBaseWidget.__init__(self, name, mode)
		self.value = None
		self.lrange = lrange
		self.urange = urange
		self.initdefault = initdefault
		self.setPositional(postional)

		gridbox = QtWidgets.QGridLayout()
		label = QtWidgets.QLabel(name)
		self.intbox = QtWidgets.QLineEdit()
		gridbox.addWidget(label, 0, 0)
		gridbox.addWidget(self.intbox, 0, 1)
		self.setLayout(gridbox)

		self.intbox.editingFinished.connect(self._on_intchanged)

		self.setValue(value)

	def _on_intchanged(self, quiet=False):
		if str(self.intbox.text()).upper() == "NONE" or str(self.intbox.text()) == "":
			self.value = None
			self.setErrorMessage(None)
			return
		try:
			if len(str(self.intbox.text()).strip())==0 : self.value=0
			else: self.value = int(self.intbox.text())
			self.setErrorMessage(None)
			self._confirm_bounds()
		except ValueError:
			self.intbox.setText("")
			self.setErrorMessage("Invalid type, int needed in %s"%self.getName())
			if self.isVisible() and not quiet: self.pmmessage.emit("Invalid type, int needed in %s"%self.getName())

	def _confirm_bounds(self):
		if self.lrange != None and (self.value < self.lrange):
			self.intbox.setText(str(self.lrange))
			if self.isVisible(): self.pmmessage.emit("Value too low for '%s', clipping to '%d'"%(self.name,self.lrange))
		if self.urange != None and (self.value > self.urange):
			self.intbox.setText(str(self.urange))
			if self.isVisible(): self.pmmessage.emit("Value too high for '%s', clipping to '%d'"%(self.name,self.urange))

	def getValue(self):
		return self.value

	def setValue(self, value, quiet=False):
		self.intbox.setText(str(value))
		self._on_intchanged(quiet=quiet)

	def setEnabled(self, state):
		self.intbox.setEnabled(state)

class PMShrinkEntryWidget(PMIntEntryWidget):
	""" A widget for shrink options. If this entry is set to <= 1 then no argument is returned """

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
			self.setErrorMessage("Invalid type, int needed in %s"%self.getName())
			if self.isVisible() and not quiet: self.pmmessage.emit("Invalid type, int needed in %s"%self.getName())

	def _confirm_bounds(self):
		if self.lrange != None and (self.value < self.lrange):
			self.value = self.lrange - 1
			self.intbox.setText(str(self.lrange-1))
			self.noarg = True
			return
		self.noarg = False

class PMFloatEntryWidget(PMBaseWidget):
	""" A widget for getting float values. Type and range are checked """

	@staticmethod
	def copyWidget(widget):
		""" Basically a copy constructor to get around QT and python limitations """
		return PMFloatEntryWidget(widget.getName(), widget.getValue(), widget.getMode(), widget.lrange, widget.urange, widget.postional, widget.initdefault)

	def __init__(self, name, value, mode, lrange=None, urange=None, postional=False, initdefault=None):
		PMBaseWidget.__init__(self, name, mode)
		self.value = None
		self.lrange = lrange
		self.urange = urange
		self.initdefault = initdefault
		self.setPositional(postional)

		gridbox = QtWidgets.QGridLayout()
		label = QtWidgets.QLabel(name)
		self.floatbox = QtWidgets.QLineEdit()
		gridbox.addWidget(label, 0, 0)
		gridbox.addWidget(self.floatbox, 0, 1)
		self.setLayout(gridbox)

		self.floatbox.editingFinished.connect(self._on_floatchanged)

		self.setValue(value)

	def _on_floatchanged(self, quiet=False):
		if str(self.floatbox.text()).upper() == "NONE" or str(self.floatbox.text()) == "":
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
			if self.isVisible() and not quiet: self.pmmessage.emit("Invalid type, float needed in '%s'"%self.getName())

	def _confirm_bounds(self):
		if self.lrange and (self.value < self.lrange):
			self.floatbox.setText(str(self.lrange))
			if self.isVisible(): self.pmmessage.emit("Value too low for '%s', clipping to '%f'"%(self.name,self.lrange))
		if self.urange and (self.value > self.urange):
			self.floatbox.setText(str(self.urange))
			if self.isVisible(): self.pmmessage.emit("Value too high for '%s', clipping to '%f'"%(self.name,self.urange))

	def getValue(self):
		return self.value

	def setValue(self, value, quiet=False):
		self.floatbox.setText(str(value))
		self._on_floatchanged(quiet=quiet)

	def setEnabled(self, state):
		self.floatbox.setEnabled(state)

class PMStringEntryWidget(PMBaseWidget):
	""" A widget for getting string values. Type is checked """

	@staticmethod
	def copyWidget(widget):
		""" Basically a copy constructor to get around QT and python limitations """
		return PMStringEntryWidget(widget.getName(), widget.getValue(), widget.getMode(), widget.postional, widget.initdefault, widget.returnNone)

	def __init__(self, name, string, mode, postional=False, initdefault=None, returnNone=False):
		PMBaseWidget.__init__(self, name, mode, returnNone=returnNone)
		self.initdefault = initdefault
		self.setPositional(postional)

		gridbox = QtWidgets.QGridLayout()
		label = QtWidgets.QLabel(name)
		self.stringbox = QtWidgets.QLineEdit()
		gridbox.addWidget(label, 0, 0)
		gridbox.addWidget(self.stringbox, 0, 1)
		self.setLayout(gridbox)

		self.stringbox.editingFinished.connect(self._on_stringchanged)

		self.setValue(string)

	def _on_stringchanged(self):
		self.string = str(self.stringbox.text())

	def getValue(self):
		# What to do with None type values? For strings, just set None to "". This should be equivalent
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

		gridbox = QtWidgets.QGridLayout()
		self.header = QtWidgets.QLabel()
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
		""" Obviously the header does give an argument """
		return None

class PMBoolWidget(PMBaseWidget):
	""" A widget for getting bool values. Type is checked """

	@staticmethod
	def copyWidget(widget):
		""" Basically a copy constructor to get around QT and python limitations """
		return PMBoolWidget(widget.getName(), widget.getValue(), widget.getMode(), widget.initdefault)

	def __init__(self, name, boolvalue, mode, initdefault=None):
		PMBaseWidget.__init__(self, name, mode)
		self.boolvalue = boolvalue
		self.initdefault = initdefault

		gridbox = QtWidgets.QGridLayout()
		self.boolbox = QtWidgets.QCheckBox(name)
		gridbox.addWidget(self.boolbox, 0, 0)
		self.setLayout(gridbox)

		self.boolbox.stateChanged[int].connect(self._on_boolchanged)

		self.setValue(self.boolvalue)

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
	""" A widget for getting filenames. Type is checked """
	pmfilename = QtCore.pyqtSignal(str)
	@staticmethod
	def copyWidget(widget):
		""" Basically a copy constructor to get around QT and python limitations """
		return PMFileNameWidget(widget.getName(), widget.filename, widget.mode, widget.browser, widget.postional, widget.initdefault, widget.checkfileexist, infolabels=False)

	def __init__(self, name, filename, mode, browser, postional=True, initdefault=None, checkfileexist=True, infolabels=True):
		PMBaseWidget.__init__(self, name, mode)
		self.initdefault = initdefault
		self.checkfileexist= checkfileexist
		self.setPositional(postional)

		gridbox = QtWidgets.QGridLayout()
		label = QtWidgets.QLabel(name)
		self.browser = browser
		self.filenamebox = QtWidgets.QLineEdit()
		self.browsebutton = QtWidgets.QPushButton("Browse")
		self.infolabel = QtWidgets.QLabel("Num Images: None")
		gridbox.addWidget(label, 0, 0)
		gridbox.addWidget(self.filenamebox, 0, 1)
		gridbox.addWidget(self.browsebutton, 0, 2)
		if infolabels: gridbox.addWidget(self.infolabel, 1, 1, 1, 2)
		self.setLayout(gridbox)

		self.filenamebox.editingFinished.connect(self._on_filenamechanged)
		self.browsebutton.clicked.connect(self._on_clicked)

		self.setValue(filename)

	def _on_filenamechanged(self):
		self.setValue(str(self.filenamebox.text()))

	def _on_cancel(self):
		self.window.close()
		self.window=None

	def _on_ok(self):
		filename = ""
		for f in self.window.getResult():
			filename += (" "+f)
		self.setValue(filename[1:])
		self.window=None

	def _on_clicked(self):
		self.window = eval(self.browser)
		self.window.ok.connect(self._on_ok)
		self.window.cancel.connect(self._on_cancel)
		self.window.setAttribute(QtCore.Qt.WA_DeleteOnClose)
		self.window.show()

	def getValue(self):
		return self.filename

	def setValue(self, filename, quiet=False):
		# Check to see if the file exists
		filename = str(filename)
		if filename.upper() == "NONE": filename=""	# If none is passed set to blank
		# Positional arguments must be space delimited for multiple files, whereas options must be comma delimited
		if not self.getPositional():
			filename = filename.replace(" ",",")
		# In some cases a file is optional
		#if self.checkfileexist:
			## We need to check if the field is blank
			#if filename == "":
				#self._onBadFile(filename, quiet)
				#return
			## In not blank then check to ensure each file is 'ok'. Not that a list of files is accepted
			#if not self._checkfiles(filename): return

		# If all is well, then we are happy
		self.filename = filename
		self.filenamebox.setText(filename)
		self.setErrorMessage(None)
		self.pmfilename.emit(self.getValue())

	def _checkfiles(self, filename):
		# Positional arguments must be space delimited for multiple files, whereas options must be comma delimited
		if self.getPositional():
			files = filename.split()
		else:
			files = filename.split(",")
		
		# If we have too many files, the user will have to wait a LOOONG time for the image check, so we just skip it
		if len(files)>20 :
			try:
				tst=EMData(files[0],0)
				nx,ny,nz=tst["nx"],tst["ny"],tst["nz"]
			except: nx,ny,nz=0,0,0

			if nx>0: self.infolabel.setText("Files: %d   %dx%dx%d"%(len(files),nx,ny,nz))
			else : self.infolabel.setText("Files: %d"%(len(files)))

			return True
			
		# Check each file
		numimages = 0
		nx,ny,nz=0,0,0
		for f in files:
			if not os.access(f, os.F_OK) and not db_check_dict(f):
				self._onBadFile(f)
				# Display the rubbish file to the user
				self.filenamebox.setText(filename)
				return False
			try:
				numimages += EMUtil.get_image_count(f)
				tst=EMData(files[0],0)
				nx,ny,nz=tst["nx"],tst["ny"],tst["nz"]
			except:
				nx,ny,nz=0,0,0
		if nx>0: self.infolabel.setText("Files: %d Images: %d  %dx%dx%d"%(len(files),numimages,nx,ny,nz))
		else : self.infolabel.setText("Files: %d Images: %d"%(len(files),numimages))

		return True

	def _onBadFile(self, filename, quiet=False):
		self.filename = None
		self.setErrorMessage("File '%s' from field '%s' does not exist"%(filename,self.getName()))
		if self.isVisible() and not quiet: self.pmmessage.emit("File '%s' from field '%s' does not exist"%(filename,self.getName()))

class PMDirectoryWidget(PMBaseWidget):
	""" A widget for display directories of a certain type """

	@staticmethod
	def copyWidget(widget):
		""" Basically a copy constructor to get around QT and python limitations """
		return PMDirectoryWidget(widget.getName(), widget.dirbasename, widget.getValue(), widget.getMode(), widget.postional, widget.initdefault)

	def __init__(self, name, dirbasename, default, mode, postional=False, initdefault=None):
		PMBaseWidget.__init__(self, name, mode)
		self.dirbasename = dirbasename
		self.initdefault = initdefault
		self.setPositional(postional)

		gridbox = QtWidgets.QGridLayout()
		label = QtWidgets.QLabel(name)
		self.combobox = PMComboBox()
		gridbox.addWidget(label, 0, 0)
		gridbox.addWidget(self.combobox, 0, 1)
		self.setLayout(gridbox)

		self.combobox.activated[str].connect(self.setValue)

		self.setValue(default)

	def updateDirs(self):
		for idx in range(self.combobox.count()):
			self.combobox.removeItem(self.combobox.count()-1)
		# This extra code allows us to have more than one type of directory
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
	""" A widget for combo boxes. Type is checked """

	@staticmethod
	def copyWidget(widget):
		""" Basically a copy constructor to get around QT and python limitations """
		return PMComboWidget(widget.getName(), widget.choices, widget.getValue(), widget.getMode(), widget.postional, widget.datatype, widget.initdefault, widget.returnNone)

	def __init__(self, name, choices, default, mode, postional=False,  datatype=str, initdefault=None, returnNone=False):
		PMBaseWidget.__init__(self, name, mode, returnNone=returnNone)
		self.initdefault = initdefault
		self.datatype=datatype	# Must be int, float or str
		self.setPositional(postional)

		gridbox = QtWidgets.QGridLayout()
		label = QtWidgets.QLabel(name)
		self.combobox = PMComboBox()
		gridbox.addWidget(label, 0, 0)
		gridbox.addWidget(self.combobox, 0, 1)
		self.setLayout(gridbox)
		# Load combo box
		self.choices = sorted(choices)
		for choice in self.choices:
			self.combobox.addItem(str(choice))

		self.combobox.activated[str].connect(self.setValue)

		self.setValue(default)

	def getValue(self):
		if str(self.combobox.currentText()) == "None":
			# In the e2 programs we actually need to specify None otherwise default behaviour will be implemented
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
			if not quiet: self.pmmessage.emit("Value '%s' not found in combobox '%s'"%(value,self.getName()))
			return


class PMComboParamsWidget(PMBaseWidget):
	""" A widget for combo boxes. Type is checked. For the combobox with params the datatype is always str """

	@staticmethod
	def copyWidget(widget):
		""" Basically a copy constructor to get around QT and python limitations """
		return PMComboParamsWidget(widget.getName(), widget.choices, widget.getValue(), widget.getMode(), widget.postional, widget.initdefault, widget.returnNone)

	def __init__(self, name, choices, default, mode, postional=False, initdefault=None, returnNone=False):
		PMBaseWidget.__init__(self, name, mode, returnNone=returnNone)
		self.initdefault = initdefault
		self.setPositional(postional)

		gridbox = QtWidgets.QGridLayout()
		label = QtWidgets.QLabel(name)
		self.combobox = PMComboBox()
		plabel = QtWidgets.QLabel("params:")
		self.params = QtWidgets.QLineEdit()
		gridbox.addWidget(label, 0, 0)
		gridbox.addWidget(self.combobox, 0, 1)
		gridbox.addWidget(plabel, 0, 2)
		gridbox.addWidget(self.params, 0, 3)
		self.setLayout(gridbox)
		# Load choices
		self.choices = sorted(choices)
		for choice in self.choices:
			self.combobox.addItem(str(choice))
		self.combobox.addItem('None')

		self.combobox.activated[str].connect(self.setValue)

		self.setValue(default)

	def getValue(self):
		""" Return the value. Concatenate if necessary """
		if str(self.combobox.currentText()) == "None":
			# In the e2 programs we actually need to specify None otherwise default behaviour will be implemented
			return "None"
		if self.params.text() == "":
			return str(self.combobox.currentText())
		else:
			return str(self.combobox.currentText())+":"+self.params.text()

	def setValue(self, value, quiet=False):
		# First parse the value, as it may contain both options and params
		if value == '': value = "None"
		values = self._parsedefault(str(value))
		# Next process the parsed value
		idx = self.combobox.findText(str(values[0]))
		if idx > -1:
			self.combobox.setCurrentIndex(idx)
		else:
			self.setErrorMessage("Value '%s' not found in combobox '%s'"%(values[0],self.getName()))
			if not quiet: self.pmmessage.emit("Value '%s' not found in combobox '%s'"%(values[0],self.getName()))
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
		self.initdefault = initdefault

		gridbox = QtWidgets.QGridLayout()
		label = QtWidgets.QLabel(name)
		label.setAlignment(QtCore.Qt.AlignVCenter | QtCore.Qt.AlignRight)
		self.combobox = PMComboBox()
		self.symnumbox = PMIntEntryWidget("Sym Number", 0, mode, lrange=0)
		gridbox.addWidget(label, 0, 0)
		gridbox.addWidget(self.combobox, 0, 1)
		gridbox.addWidget(self.symnumbox, 0, 2, 1, 2)
		self.setLayout(gridbox)

		for i in ['icos','oct','tet','c','d','h']: self.combobox.addItem(i)

		self.symnumbox.pmmessage[str].connect(self._on_message)

		self.setValue(default)

	def _on_message(self, message):
		self.pmmessage.emit(message)

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
			if not quiet: self.pmmessage.emit("'%s' not a valid symmetry!!!"%value)
			return
		self.symnumbox.setValue(defsymnum)
		self.setErrorMessage(None)

	def getErrorMessage(self):
		if self.errormessage: return self.errormessage
		if self.symnumbox.getErrorMessage(): return self.symnumbox.getErrorMessage()

class PMAutoMask3DWidget(PMBaseWidget):
	""" A widget for getting automask 3D input """

	@staticmethod
	def copyWidget(widget):
		""" Basically a copy constructor to get around QT and python limitations """
		return PMAutoMask3DWidget(widget.getName(), widget.getValue(), widget.getMode(), widget.initdefault)

	def __init__(self, name, default, mode, initdefault=None):
		PMBaseWidget.__init__(self, name, mode)
		self.initdefault = initdefault

		gridbox = QtWidgets.QGridLayout()
		self.automask3dbool = QtWidgets.QCheckBox("Auto Mask 3D")
		self.params = []
		self.params.append(PMFloatEntryWidget("Threshold", 0.8, mode))
		self.params.append(PMIntEntryWidget("Radius", 30, mode))
		self.params.append(PMIntEntryWidget("Mask Dilations", 5, mode))
		self.params.append(PMIntEntryWidget("Gaussian Dilations", 5, mode))
		self.params.append(PMIntEntryWidget("NMax", 30, mode))
		gridbox.addWidget(self.automask3dbool, 0, 0)
		gridbox.addWidget(self.params[0], 1, 0)
		gridbox.addWidget(self.params[1], 1, 1)
		gridbox.addWidget(self.params[2], 1, 2)
		gridbox.addWidget(self.params[3], 2, 0)
		gridbox.addWidget(self.params[4], 2, 1)
		self.setLayout(gridbox)

		self.automask3dbool.stateChanged[int].connect(self._on_boolchanged)
		self.params[0].pmmessage[str].connect(self._on_message)
		self.params[1].pmmessage[str].connect(self._on_message)
		self.params[2].pmmessage[str].connect(self._on_message)
		self.params[3].pmmessage[str].connect(self._on_message)
		self.params[4].pmmessage[str].connect(self._on_message)

		self.setValue(default)

	def _on_boolchanged(self):
		for widget in self.params:
			widget.setEnabled(self.automask3dbool.isChecked())

	def _on_message(self, message):
		self.pmmessage.emit(message)

	def setValue(self, value, quiet=False):
		# If value is "" of None, set bool to false
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
		# Concatenate things
		for i in range(len(self.params)):
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

class PMTableBase(PMBaseWidget):
	""" A base widget for making tables. To make a table class subclass this base and add a line to PMGUIWidget in e2projectmanager.py
	in order to instantiate it using directions from an e2program. See Wiki for more info """
	def __init__(self, name, mode, postional=False, initdefault=None):
		PMBaseWidget.__init__(self, name, mode)
		self.setPositional(postional)
		self.initdefault = initdefault

		gridbox = QtWidgets.QGridLayout()
		self.tablewidget = QtWidgets.QTableWidget()
		gridbox.addWidget(self.tablewidget, 0, 0)
		self.tablewidget.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)	# Readonly table
		self.setLayout(gridbox)

	def updateTable(self):
		""" Update FSC table. You must implement this function to build the table"""
		raise NotImplementedError("Sub class must reimplement 'getValue' function")

	def getValue(self):
		""" Must implement this to return a value from the table """
		raise NotImplementedError("Sub class must reimplement 'getValue' function")

	def setValue(self):
		""" Must implement this function to set a value in the table """
		raise NotImplementedError("Sub class must reimplement 'getValue' function")

class PMFSCTableWidget(PMTableBase):
	""" A widget for generating FSC tables"""
	pmmessage = QtCore.pyqtSignal(str)

	@staticmethod
	def copyWidget(widget):
		""" Basically a copy constructor to get around QT and python limitations """
		return PMFSCTableWidget(widget.getName(), widget.getValue(), widget.getMode(), widget.postional, widget.initdefault)

	def __init__(self, name, default, mode, postional=False, initdefault=None, resize=False):
		PMTableBase.__init__(self, name, mode, postional, initdefault)

		# Table stuff
		self.tablewidget.setColumnCount(4)
		self.tablewidget.setHorizontalHeaderLabels(["Refine", "# Iter", "Masked .143", "Unmasked .143"])
		self.tablewidget.horizontalHeader().setSectionResizeMode(QtWidgets.QHeaderView.ResizeToContents)
		self.tablewidget.horizontalHeader().setHighlightSections(False)
		self.tablewidget.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)	# select rows
		self.tablewidget.setSelectionMode(QtWidgets.QAbstractItemView.SingleSelection)	# single selection

		self.tablewidget.setRowCount(0)
		self.patterns = ["refine","frealign","multi"]
		self.tablewidget.cellDoubleClicked[int, int].connect(self.loadFSC)

		# Now update table
		self.setValue(default)
		if resize: self.resize(550,300)

	def getValue(self):
		value = self.tablewidget.item(self.tablewidget.currentRow(),0)
		if value:
			return value.text()
		else:
			self.pmmessage.emit("Warning NOTHING IS SELECTED IN THE TABLE!!!")
			return ""

	def setValue(self, value, quiet=False):
		self.updateTable()
		str(value)
		wlist = self.tablewidget.findItems(str(value), QtCore.Qt.MatchExactly)
		# Only set if something was found
		if wlist:
			self.tablewidget.setCurrentItem(wlist[0])

	def loadFSC(self, row, col):
		""" Display the FSC curve. This is a callback for double clicking"""
		if not self.tablewidget.item(row, 1):
			msg = "Rubbish!!! No FSC curves to plot."
			print(msg)
			self.pmmessage.emit("Rubbish!!! No FSC curves to plot.")
			return
		
		fsccmd=["e2display.py --plot"]
		path=str(self.tablewidget.item(row, 0).text())
		fls=["{}/{}".format(path,i) for i in os.listdir(path) if i[:11]=="fsc_masked_"]
		fsccmd.extend(sorted(fls))
		fsccmd=" ".join(fsccmd)
		
		# Now load the FSC curves
		msg = "Loading FSC curves, please wait..."
		print(msg)
		self.pmmessage.emit(msg)
		subprocess.Popen(fsccmd, shell=True)

	def updateTable(self):
		""" Update FSC table"""
		dirs = []
		for pattern in self.patterns:
			dirs.extend(glob.glob("%s*"%pattern))
		
		dirs=[i for i in dirs if os.path.isdir(i)]
		self.tablewidget.setRowCount(len(dirs))

		for i, directory in enumerate(sorted(dirs)):
			# Load each directory
			qwi_dirname = QtWidgets.QTableWidgetItem(str(directory))
			self.tablewidget.setItem(i, 0, qwi_dirname)

			fscs=sorted([ii for ii in os.listdir(directory) if ii[:11]=="fsc_masked_"])
			niter=len(fscs)
			if "fsc_masked_00.txt" in fscs : niter-=1
			
			self.tablewidget.setItem(i, 1, QtWidgets.QTableWidgetItem(str(niter)))
			
			try:
				# We use a running average of 5 points to compute the threshold
				xyd=XYData()
				xyd.read_file("{}/{}".format(directory,fscs[-1]))
				for ii in range(2,xyd.get_size()-2):
					v=old_div((xyd.get_y(ii-2)+xyd.get_y(ii-1)+xyd.get_y(ii)+xyd.get_y(ii+1)+xyd.get_y(ii+2)),5.0)
					if v<0.143 : break
				
				self.tablewidget.setItem(i,2,QtWidgets.QTableWidgetItem("{:1.1f}".format(old_div(1.0,xyd.get_x(ii-1)))))
			except:
				self.tablewidget.setItem(i,2,QtWidgets.QTableWidgetItem("?"))

			try:
				# We use a running average of 5 points to compute the threshold
				xyd=XYData()
				xyd.read_file("{}/fsc_un{}".format(directory,fscs[-1][4:]))
				for ii in range(2,xyd.get_size()-2):
					v=old_div((xyd.get_y(ii-2)+xyd.get_y(ii-1)+xyd.get_y(ii)+xyd.get_y(ii+1)+xyd.get_y(ii+2)),5.0)
					if v<0.143 : break
				
				self.tablewidget.setItem(i,3,QtWidgets.QTableWidgetItem("{:1.1f}".format(old_div(1.0,xyd.get_x(ii-1)))))
			except:
				self.tablewidget.setItem(i,3,QtWidgets.QTableWidgetItem("?"))

			

	# I lifted this code from Daivid's SPR workflow module
	def find_first_point_5_crossing(self,xaxis,yaxis,thr=0.5):
		'''
		Find the first 0.5 crossing in the FSC - interpolate and try to return an accurate estimate
		Disregards the first five entries
		Assumes the Nyquist frequency is correct
		'''
		idx = 0
		if len(yaxis) > 5:
			idx = 6

		soln = -1
		while (idx < len(yaxis)-1 ):
			if yaxis[idx] >= thr and yaxis[idx+1] <= thr:
				v1 = yaxis[idx]
				v2 = yaxis[idx+1]
				if v1 != v2:
					d = v1-v2
					offset = v1-thr
					interp = old_div(offset,d)
					soln = idx+interp
				else: soln = idx
				break
			else:
				idx += 1

		try:
			if soln == -1:
				return "???"
			elif int(soln) == soln:
				return "%.1f" %(old_div(1.0,xaxis(soln)))
			else:
				# Interpolated frequency
				return "%.1f" %(old_div(1.0,(old_div(soln,len(yaxis))*xaxis[-1])))
		except:
			return "invalid"
