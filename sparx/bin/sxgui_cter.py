#!/usr/bin/env python
from __future__ import print_function
#
# Author: Toshio Moriya, 12/21/2015 (toshio.moriya@mpi-dortmund.mpg.de)
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
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
#
#

from EMAN2 import *
from EMAN2db import db_open_dict, db_close_dict, db_check_dict, db_list_dicts
from eman2_gui.emshape import *
from eman2_gui.valslider import *
from eman2_gui.emplot2d import EMPlot2DWidget
from eman2_gui.emapplication import EMApp
from eman2_gui.emimage2d import EMImage2DWidget

from OpenGL import GL,GLUT
from math import *
import os
import sys
# from numpy import array,arange
import numpy as np
import traceback

try:
	from PyQt4 import QtCore, QtGui, QtOpenGL
	from PyQt4.QtCore import Qt
	from PyQt4.QtCore import QTimer
except:
	print("Warning: PyQt4 must be installed")
	sys.exit(1)

from sparx import *
from optparse import OptionParser
from statistics import hist_list

'''
Scipy now calls numpy 1.15, which generates numerous warnings of the form 
"RuntimeWarning: numpy.dtype size changed, may indicate binary incompatibility. Expected 96, got 88".
Filterwarnings suppreses this message.
'''
import warnings
warnings.filterwarnings("ignore", message="numpy.dtype")
from scipy import interpolate, optimize

from morphology import ampcont2angle

def main():
	progname = os.path.basename(sys.argv[0])
	usage = progname + """  cter_ctf_file 
	This GUI application is designed for the evaluation of micrographs using the parameters outputed by CTER.
	"""
	parser = OptionParser(usage, version=SPARXVERSION)
	# No options!!! Does not need to call parser.add_option()
	
	(options, args) = parser.parse_args(sys.argv[1:])
	
	if len(args) > 2:
		print("see usage " + usage)
		sys.exit()
	
	app=EMApp()
	
	cter_ctf_file = None
	if len(args) == 1 and args[0] != "":
		cter_ctf_file = args[0]
	# else: # Do nothing
	
	# Make sure main window is shown, raised, and activated upon startup.
	# gui=SXGuiCter(cter_ctf_file)
	gui = SXGuiCter()
	gui.show()
	gui.raise_()
	gui.activateWindow()
	
	# read CTER partres file if necessary
	if cter_ctf_file != None:
		gui.readCterPartresFile(os.path.relpath(cter_ctf_file))
	
#	try:
#		# gui.wimgmicthumb.qt_parent.raise_() # wimgmicthumb should not be visible upon start-up
#		# gui.wfft.qt_parent.raise_()
#		gui.wplotrotavgcoarse.qt_parent.raise_()
#		gui.wplotrotavgfine.qt_parent.raise_()
#		gui.whistparam.qt_parent.raise_()
#		gui.wscatterparam.qt_parent.raise_()
#		gui.raise_()
#		gui.activateWindow()
#	except: 
#		print "Received unexpected exception in main(): ", sys.exc_info()[0]
#		exc_type, exc_value, exc_traceback = sys.exc_info()
#		traceback.print_exception(exc_type, exc_value, exc_traceback)
#		# MRK_NOTE: 2015/12/17 Toshio Moriya
#		# Another way to print out exception info...
#		# lines = traceback.format_exception(exc_type, exc_value, exc_traceback)
#		# print ''.join('!! ' + line for line in lines)
#		pass
	
	# NOTE: 2016/03/08 Toshio Moriya
	# Unfortunately, EMApp.execute() will print logid (e.g. "None") upon the exit...
	app.execute()

class SXListWidget(QtGui.QListWidget):
	"""Exactly like a normal list widget but intercepts a few keyboard events"""
	def keyPressEvent(self,event):
		if event.key() in (Qt.Key_Up,Qt.Key_Down) :
			QtGui.QListWidget.keyPressEvent(self,event)
			return
		
		self.emit(QtCore.SIGNAL("keypress"),event)

class SXPlot2DWidget(EMPlot2DWidget):
	
	mouseup = QtCore.pyqtSignal(QtGui.QMouseEvent)

	def full_refresh(self):
		'''
		This function is called from resizeGL and from the inspector when somebody toggles the display of a line
		'''
		self.needupd=1
		self.del_shapes(("xcross","ycross","lcross","Circle")) # for EMPlot2DInspector
		self.del_shapes(("error_astig","error_def","error_ctf")) # for SXGuiCter.wplotrotavgcoarse & SXGuiCter.wplotrotavgfine
		# self.del_shapes(("error_astig","error_def")) # for SXGuiCter.wplotrotavgcoarse & SXGuiCter.wplotrotavgfine
		#self.del_shapes(("hist_param_shape_value","hist_param_shape_unapply_threshold_lower","hist_param_shape_apply_threshold_lower",
				   #"hist_param_shape_unapply_threshold_upper", "hist_param_shape_apply_threshold_upper", "hist_param_shape_label")) # for SXGuiCter.whistparam
		# self.del_shapes(("hist_param_shape_label")) # for SXGuiCter.whistparam
		# self.del_shapes(("plot_param_shape_label")) # for SXGuiCter.wscatterparam
	
	def mouseReleaseEvent(self, event):
		EMPlot2DWidget.mouseReleaseEvent(self,event)
		if event.button()==Qt.LeftButton:
			self.mouseup.emit(event)  #self.emit(QtCore.SIGNAL("mouseup"),event)

class SXLogoButton(QtGui.QLabel):
	def __init__(self, imagename, logo_width, parent = None):
		super(SXLogoButton, self).__init__(parent)

		# Width of logo image
		logo_file_path = '{0}{1}'.format(get_image_directory(), imagename)

		# Style of widget
		self.setFixedSize(logo_width, logo_width)
		self.customButtonStyle = """
			SXLogoButton {{background-color: rgba(0, 0, 0, 0); border: 0px solid black; border-radius: 0px; image: url("{0}");}}
			""".format(logo_file_path)

		# Set style and add click event
		self.setStyleSheet(self.customButtonStyle)
		
	def add_sxmenu_item_btn_widget(self, sxmenu_item_btn_subarea_widget):
		sxmenu_item_btn_subarea_widget.addWidget(self)
		
class SXGuiCter(QtGui.QWidget):
# 	def __init__(self, cter_ctf_file = None):
	def __init__(self):
		"""Implements the CTF fitting dialog using various EMImage and EMPlot2D widgets
		'data' is a list of (filename,ctf,im_1d,bg_1d,quality)
		'parms' is [box size,ctf,box coord,set of excluded boxnums,quality,oversampling]
		"""
		#try:
			#from emimage2d import EMImage2DWidget
		#except:
			#print("Cannot import EMAN image GUI objects (EMImage2DWidget)")
			#sys.exit(1)
		
		super(SXGuiCter, self).__init__(None)
		
###		# MRK_TEST: Flags to test experimental functions
###		self.is_enable_max_power = False # MRK_TEST: Can this be an option in future?
		
		# NOTE: 2016/03/15 Toshio Moriya
		# Due to the representation error of float number, 
		# default thresholds using max/min value of each parameter can be different after
		# saving/loading the threshold values from the file.
		# This problem is solved by round() function for thresholding parameter values
		# 
		# About Precision: 
		# Double precision numbers have 53 bits (16 digits ~~ 15.95 digits = 53*log10(2)) of precision and 
		# regular floats have 24 bits (8 digits ~~ 7.22 digits = 24*log10(2)) of precision. 
		# The floating point in python uses double precision to store the values, and 
		# usually taken to be 15 for practical purposes
		# 
		self.round_ndigits = 15
		
#		self.setWindowIcon(QtGui.QIcon(get_image_directory() + "ctf.png"))
		self.setWindowIcon(QtGui.QIcon(get_image_directory()+"sparxicon.png"))

#		# NOTE: 2016/03/08 Toshio Moriya
#		# Checked the following window flags and found out ...
#		# (1) I could not find a way to disabled a close button at all.
#		# (2) Qt.CustomizeWindowHint is best to disabled a minimize AND maximize button.
#		# (3) With Qt.WindowTitleHint and Qt.WindowSystemMenuHint, I could not find a way to disable maximize button.
#		# (4) It seems to be impossible to turn off Qt.Window flag here...
#		print "MRK_DEBUG: Qt.Widget Window Flags: 0x%08x " % (Qt.Widget)
#		print "MRK_DEBUG: Qt.Window Window Flags: 0x%08x " % (Qt.Window)
#		print "MRK_DEBUG: Qt.Dialog Window Flags: 0x%08x " % (Qt.Dialog)
#		print "MRK_DEBUG: Qt.Popup Window Flags: 0x%08x " % (Qt.Popup)
#		print "MRK_DEBUG: Qt.Tool Window Flags: 0x%08x " % (Qt.Tool)
#		print "MRK_DEBUG: Qt.SubWindow Window Flags: 0x%08x " % (Qt.SubWindow)
#		
#		print "MRK_DEBUG: Original Window Flags: 0x%08x " % (self.windowFlags())
#		
#		self.setWindowFlags((self.windowFlags()| Qt.CustomizeWindowHint)) # Turns off the default window title hints.
#		self.setWindowFlags((self.windowFlags()| Qt.WindowTitleHint)) # Gives the window a title bar.
#		self.setWindowFlags((self.windowFlags()| Qt.WindowSystemMenuHint)) # Disabled minimize icon button in window title bar
		# self.setWindowFlags((self.windowFlags()| Qt.CustomizeWindowHint) & ~Qt.WindowMinimizeButtonHint) # OK
		# self.setWindowFlags((self.windowFlags()| Qt.CustomizeWindowHint) & ~Qt.WindowMaximizeButtonHint) # OK
#		self.setWindowFlags((self.windowFlags()| Qt.CustomizeWindowHint) & ~Qt.WindowCloseButtonHint) # NG
#		self.setWindowFlags((self.windowFlags()| Qt.WindowTitleHint) & ~Qt.WindowMinimizeButtonHint) # OK
#		self.setWindowFlags((self.windowFlags()| Qt.WindowTitleHint) & ~Qt.WindowMaximizeButtonHint) # NG
#		self.setWindowFlags((self.windowFlags()| Qt.WindowTitleHint) & ~Qt.WindowCloseButtonHint) # NG
#		self.setWindowFlags((self.windowFlags()| Qt.WindowSystemMenuHint) & ~Qt.WindowMinimizeButtonHint) # OK
#		self.setWindowFlags((self.windowFlags()| Qt.WindowSystemMenuHint) & ~Qt.WindowMaximizeButtonHint) # NG
#		self.setWindowFlags((self.windowFlags()| Qt.WindowSystemMenuHint) & ~Qt.WindowCloseButtonHint) # NG
#		self.setWindowFlags(Qt.Widget| Qt.CustomizeWindowHint) # It looks like, we always need Qt.CustomizeWindowHint to override window flags
#		self.setWindowFlags(Qt.Widget| Qt.CustomizeWindowHint | Qt.WindowTitleHint | Qt.WindowSystemMenuHint) # OK, but it is still Qt.Window ...
#		self.setWindowFlags(Qt.Dialog| Qt.CustomizeWindowHint) # It looks like, we always need Qt.CustomizeWindowHint to override window flags
#		self.setWindowFlags(Qt.Dialog| Qt.CustomizeWindowHint | Qt.WindowTitleHint | Qt.WindowSystemMenuHint) # OK, but it is still Qt.Window ...
#		self.setWindowFlags(Qt.SubWindow| Qt.CustomizeWindowHint) # It looks like, we always need Qt.CustomizeWindowHint to override window flags
#		self.setWindowFlags(Qt.SubWindow| Qt.CustomizeWindowHint | Qt.WindowTitleHint | Qt.WindowSystemMenuHint) # OK, but it is still Qt.Window ...
#		self.setWindowFlags(Qt.Widget & ~Qt.WindowSystemMenuHint) # 
#		
#		print "MRK_DEBUG: Edited Window Flags: 0x%08x " % (self.windowFlags())
#		
#		window_flags = self.windowFlags()
#		if Qt.WindowCloseButtonHint == (window_flags & Qt.WindowCloseButtonHint):
#			print "MRK_DEBUG: Window Close Button 0x%08x should be ON: 0x%08x" % (Qt.WindowCloseButtonHint, window_flags & Qt.WindowCloseButtonHint)
#		else:
#			print "MRK_DEBUG: Window Close Button 0x%08x should be OFF: 0x%08x" % (Qt.WindowCloseButtonHint, window_flags & Qt.WindowCloseButtonHint)
#			
		self.installEventFilter(self) # Necessary for self.eventFilter()
		
		# Enumerate indices
		self.enumerate_params_old()
		self.enumerate_params_new()
		self.map_value_items()
		self.enumerate_1drot_indices()
		self.enumerate_sort_indices()
		self.enumerate_map_list_indices()
		self.map_sorted_items()
		self.enumerate_histograms()
		self.enumerate_histogram_pulldowns()
		self.map_histogram_items()
		self.enumerate_threshold_controls()
		self.enumerate_threshold_map_indices()
		self.map_threshold_control()
		self.enumerate_1drot_graphs()
		self.enumerate_graph_items()
		self.map_display_checkboxes()
		self.enumerate_threshold_status()
		self.enumerate_threshold_entries()
		self.map_threshold_status()
		
		self.initialize_paths()
		self.initialize_sorting()
		
		self.child_status_list = []  # NOT USED?
		
		self.set_default_display()
		
		self.initialize_popup_windows()
		
		# Emit signals
#		self.wfft.connect(self.wfft,QtCore.SIGNAL("mousedown"),self.fftmousedown)
#		self.wfft.connect(self.wfft,QtCore.SIGNAL("mousedrag"),self.fftmousedrag)
#		self.wfft.connect(self.wfft,QtCore.SIGNAL("mouseup")  ,self.fftmouseup)
#		self.wimgmicthumb.connect(self.wimgmicthumb,QtCore.SIGNAL("mousedown"),self.imgmicthumbmousedown)
#		self.wimgmicthumb.connect(self.wimgmicthumb,QtCore.SIGNAL("mousedrag"),self.imgmicthumbmousedrag)
#		self.wimgmicthumb.connect(self.wimgmicthumb,QtCore.SIGNAL("mouseup")  ,self.imgmicthumbmouseup)
#		self.wplotrotavgcoarse.connect(self.wplotrotavgcoarse,QtCore.SIGNAL("mousedown"),self.plotmousedown)
#		self.wplotrotavgfine.connect(self.wplotrotavgfine,QtCore.SIGNAL("mousedown"),self.plotmousedown)
		self.whistparam.mouseup.connect(self.histparammouseup)  #self.whistparam.connect(self.whistparam,QtCore.SIGNAL("mouseup"),self.histparammouseup)
		self.wscatterparam.mouseup.connect(self.plotparammouseup)  #self.wscatterparam.connect(self.wscatterparam,QtCore.SIGNAL("mouseup"),self.plotparammouseup)
		
		self.draw_main_window()
		self.signal_handler()
		self.set_size_popup_windows()
		
		# Try to recover sizes & positions of windows of the previous GUI session
		E2loadappwin("sxgui_cter","main",self)
		E2loadappwin("sxgui_cter","fft",self.wfft.qt_parent)
		E2loadappwin("sxgui_cter","imgmicthumb",self.wimgmicthumb.qt_parent)
		E2loadappwin("sxgui_cter","plotcoarse",self.wplotrotavgcoarse.qt_parent)
		E2loadappwin("sxgui_cter","plotfine",self.wplotrotavgfine.qt_parent)
		E2loadappwin("sxgui_cter","histparam",self.whistparam.qt_parent)
		E2loadappwin("sxgui_cter","plotparam",self.wscatterparam.qt_parent)
		
#		if self.cter_entry_list:
# #			self.wfft.show()
#			self.whistparam.show()
#			self.wplotrotavgcoarse.show()
		
		### This section is responsible for background updates
		self.busy = False
#		self.needupdate = True
		self.needredisp = False
#		self.procthread = None
#		self.errors = None # used to communicate errors back from the reprocessing thread
		
		self.timer = QTimer()
		self.timer.timeout.connect(self.timeOut)
		self.timer.start(100)
		
#		# Finally, read CTER partres file if necessary
#		if cter_ctf_file != None:
#			self.readCterPartresFile(os.path.relpath(cter_ctf_file))
		
	def enumerate_params_old(self):
		#
		# NOTE: 2017/11/22 Toshio Moriya
		# The following code is to support the old format of CTER partres file. It should be removed near future
		# 
		# Define enumerators for the old format of CTER partres file.
		i_enum = -1
		i_enum += 1; self.idx_old_cter_def          = i_enum # defocus [um]
		i_enum += 1; self.idx_old_cter_cs           = i_enum # Cs [mm]
		i_enum += 1; self.idx_old_cter_vol          = i_enum # voltage[kV]
		i_enum += 1; self.idx_old_cter_apix         = i_enum # pixel size [A]
		i_enum += 1; self.idx_old_cter_bfactor      = i_enum # B-factor [A^2]
		i_enum += 1; self.idx_old_cter_ac           = i_enum # amplitude contrast [%]
		i_enum += 1; self.idx_old_cter_astig_amp    = i_enum # astigmatism amplitude [um]
		i_enum += 1; self.idx_old_cter_astig_ang    = i_enum # astigmatism angle [degree]
		i_enum += 1; self.idx_old_cter_sd_def       = i_enum # std dev of defocus [um]
		i_enum += 1; self.idx_old_cter_sd_ac        = i_enum # std dev of amplitude contrast [%]
		i_enum += 1; self.idx_old_cter_sd_astig_amp = i_enum # std dev of astigmatism amp [A]
		i_enum += 1; self.idx_old_cter_sd_astig_ang = i_enum # std dev of astigmatism angle [degree]
		i_enum += 1; self.idx_old_cter_cv_def       = i_enum # coefficient of variation of defocus [%]
		i_enum += 1; self.idx_old_cter_cv_astig_amp = i_enum # coefficient of variation of astigmatism amp [%]
		i_enum += 1; self.idx_old_cter_spectra_diff = i_enum # average of differences between with- and without-astig. experimental 1D spectra at extrema
		i_enum += 1; self.idx_old_cter_error_def    = i_enum # frequency at which signal drops by 50% due to estimated error of defocus alone [1/A]
		i_enum += 1; self.idx_old_cter_error_astig  = i_enum # frequency at which signal drops by 50% due to estimated error of defocus and astigmatism [1/A]
		i_enum += 1; self.idx_old_cter_error_ctf    = i_enum # limit frequency by CTF error [1/A] 
		i_enum += 1; self.idx_old_cter_mic_name     = i_enum # micrograph name
		i_enum += 1; self.n_idx_old_cter            = i_enum

	def enumerate_params_new(self):
		# Define enumerators for mapping of parameter value items (line edit widgets)
		i_enum = -1
		i_enum += 1; self.idx_cter_id               = i_enum # <extra> entry id
		i_enum += 1; self.idx_cter_select           = i_enum # <extra> selected state
		i_enum += 1; self.idx_cter_def              = i_enum # defocus [um]
		i_enum += 1; self.idx_cter_cs               = i_enum # Cs [mm]
		i_enum += 1; self.idx_cter_vol              = i_enum # voltage[kV]
		i_enum += 1; self.idx_cter_apix             = i_enum # pixel size [A]
		i_enum += 1; self.idx_cter_bfactor          = i_enum # B-factor [A^2]
		i_enum += 1; self.idx_cter_total_ac         = i_enum # total amplitude contrast [%]
		i_enum += 1; self.idx_cter_astig_amp        = i_enum # astigmatism amplitude [um]
		i_enum += 1; self.idx_cter_astig_ang        = i_enum # astigmatism angle [degree]
		i_enum += 1; self.idx_cter_sd_def           = i_enum # std dev of defocus [um]
		i_enum += 1; self.idx_cter_sd_total_ac      = i_enum # std dev of total amplitude contrast [%]
		i_enum += 1; self.idx_cter_sd_astig_amp     = i_enum # std dev of astigmatism amp [A]
		i_enum += 1; self.idx_cter_sd_astig_ang     = i_enum # std dev of astigmatism angle [degree]
		i_enum += 1; self.idx_cter_cv_def           = i_enum # coefficient of variation of defocus [%]
		i_enum += 1; self.idx_cter_cv_astig_amp     = i_enum # coefficient of variation of astigmatism amp [%]
		i_enum += 1; self.idx_cter_error_def        = i_enum # frequency at which signal drops by 50% due to estimated error of defocus alone [1/A]
		i_enum += 1; self.idx_cter_error_astig      = i_enum # frequency at which signal drops by 50% due to estimated error of defocus and astigmatism [1/A]
		i_enum += 1; self.idx_cter_error_ctf        = i_enum # limit frequency by CTF error [1/A] 
		i_enum += 1; self.idx_cter_max_freq         = i_enum # visual-impression-based maximum frequency limit [A] (e.g. max frequency of relion; CCC between neighbour zero-crossing pair)
		i_enum += 1; self.idx_cter_reserved         = i_enum # reserved spot for maximum frequency limit or error criterion. possibly originated from external program (e.g. CTF figure of merit of RELION)
		i_enum += 1; self.idx_cter_const_ac         = i_enum # constant amplitude contrast [%]
		i_enum += 1; self.idx_cter_phase_shift      = i_enum # phase shift [degrees]
		i_enum += 1; self.idx_cter_mic_name         = i_enum # micrograph name
###		if self.is_enable_max_power == True: 
###			i_enum += 1; self.idx_cter_max_power    = i_enum # MRK_TEST: <extra> maximum power in experimental rotational average (with astigmatism)
		i_enum += 1; self.n_idx_cter                = i_enum
		self.n_idx_cter_extra = 2
###		if self.is_enable_max_power == True: 
###			self.n_idx_cter_extra += 1 # MRK_TEST:
		
		# Define enumerators for items of each entry in value_map_list 
		i_enum = -1
		i_enum += 1; self.idx_cter_item_label   =  i_enum
		i_enum += 1; self.idx_cter_item_widget  =  i_enum
		i_enum += 1; self.n_idx_cter_item       =  i_enum
		
	def map_value_items(self):
		# Map parameter value items (line edit widgets)
		self.value_map_list = [None] * self.n_idx_cter
		self.value_map_list[self.idx_cter_id]           = ["CTER ID", None]
		self.value_map_list[self.idx_cter_select]       = ["Select", None]
		self.value_map_list[self.idx_cter_def]          = ["Defocus [um]", None]
		self.value_map_list[self.idx_cter_cs]           = ["Cs [mm]", None]
		self.value_map_list[self.idx_cter_vol]          = ["Voltage [kV]", None]
		self.value_map_list[self.idx_cter_apix]         = ["Pixel Size [A]", None]
		self.value_map_list[self.idx_cter_bfactor]      = ["B-factor [A^2]", None]
		self.value_map_list[self.idx_cter_total_ac]     = ["Total Amp. Contrast [%]", None]
		self.value_map_list[self.idx_cter_astig_amp]    = ["Astig. Amp.[um]", None]
		self.value_map_list[self.idx_cter_astig_ang]    = ["Astig. Ang.[deg]", None]
		self.value_map_list[self.idx_cter_sd_def]       = ["Defocus SD [um]", None]
		self.value_map_list[self.idx_cter_sd_total_ac]  = ["Total Amp. Contrast SD [%]", None]
		self.value_map_list[self.idx_cter_sd_astig_amp] = ["Astig. Amp. SD [um]", None]
		self.value_map_list[self.idx_cter_sd_astig_ang] = ["Astig. Ang. SD [deg]", None]
		self.value_map_list[self.idx_cter_cv_def]       = ["Defocus CV [%]", None]
		self.value_map_list[self.idx_cter_cv_astig_amp] = ["Astig. Amp. CV [%]", None]
		self.value_map_list[self.idx_cter_error_def]    = ["Defocus Freq. Limit [1/A]", None]
		self.value_map_list[self.idx_cter_error_astig]  = ["Astig. Freq. Limit [1/A]", None]
		self.value_map_list[self.idx_cter_error_ctf]    = ["CTF Freq. Limit [1/A]", None]
		self.value_map_list[self.idx_cter_max_freq]     = ["Max Freq. [A]", None]
		self.value_map_list[self.idx_cter_reserved]     = ["Reserved", None]
		self.value_map_list[self.idx_cter_const_ac]     = ["Const. Amp. Contrast [%]", None]
		self.value_map_list[self.idx_cter_phase_shift]  = ["Phase Shift [deg]", None]
		self.value_map_list[self.idx_cter_mic_name]     = ["Micrograph", None]
#		self.value_map_list[self.idx_cter_pwrot_name]   = ["PW. Rot. File", None]
###		if self.is_enable_max_power == True: 
###			self.value_map_list[self.idx_cter_max_power] = ["Max Power", None] # MRK_TEST:
		
	def enumerate_1drot_indices(self):
		# Define enumerators for curves of 1D power spectrum & CTF fitting plot
		i_enum = -1
		i_enum += 1; self.idx_rotinf_cter_id        = i_enum # line number == cter id
		i_enum += 1; self.idx_rotinf_freq           = i_enum # spatial frequency (1/A)
		i_enum += 1; self.idx_rotinf_exp_no_astig   = i_enum # experimental rotational average (no astigmatism)
		i_enum += 1; self.idx_rotinf_fit_no_astig   = i_enum # fitted rotational average (no astigmatism)
		i_enum += 1; self.idx_rotinf_exp_with_astig = i_enum # experimental rotational average (with astigmatism)
		i_enum += 1; self.idx_rotinf_fit_with_astig = i_enum # fitted rotational average (with astigmatism)
		i_enum += 1; self.idx_rotinf_exp_background = i_enum # experimental rotational average, background-subtracted
		i_enum += 1; self.idx_rotinf_fit_envelope   = i_enum # fitted rotational average, with envelope applied
		i_enum += 1; self.n_idx_rotinf              = i_enum
		
	def enumerate_sort_indices(self):
		# Define enumerators for mapping of sorting items (combo box widget)
		i_enum = -1
		i_enum += 1; self.idx_sort_id           = i_enum
		i_enum += 1; self.idx_sort_mic_names    = i_enum
		i_enum += 1; self.idx_sort_def          = i_enum
		i_enum += 1; self.idx_sort_total_ac     = i_enum
		i_enum += 1; self.idx_sort_astig_amp    = i_enum
		i_enum += 1; self.idx_sort_astig_ang    = i_enum
###		i_enum += 1; self.idx_sort_sd_def       = i_enum
		i_enum += 1; self.idx_sort_cv_def       = i_enum
		i_enum += 1; self.idx_sort_sd_total_ac  = i_enum
###		i_enum += 1; self.idx_sort_sd_astig_amp = i_enum
		i_enum += 1; self.idx_sort_cv_astig_amp = i_enum
		i_enum += 1; self.idx_sort_sd_astig_ang = i_enum
		i_enum += 1; self.idx_sort_error_def    = i_enum
		i_enum += 1; self.idx_sort_error_astig  = i_enum
		i_enum += 1; self.idx_sort_error_ctf    = i_enum
		i_enum += 1; self.idx_sort_max_freq     = i_enum
		i_enum += 1; self.idx_sort_reserved     = i_enum
		i_enum += 1; self.idx_sort_phase_shift  = i_enum
###		if self.is_enable_max_power == True: 
###			i_enum += 1; self.idx_sort_max_power = i_enum # MRK_TEST:
		i_enum += 1; self.n_idx_sort            = i_enum
		
	def enumerate_map_list_indices(self):
		# Define enumerators for items of each entry in sort_map_list
		i_enum = -1
		i_enum += 1; self.idx_sort_item_idx_cter =  i_enum
		i_enum += 1; self.n_idx_sort_item        =  i_enum
		
	def map_sorted_items(self):
		# Map sorting items (combo box widget)
		# Includes mapping from idx_sort to idx_cter
		self.sort_map_list = [None] * self.n_idx_sort
		self.sort_map_list[self.idx_sort_id]           = [self.idx_cter_id]
		self.sort_map_list[self.idx_sort_mic_names]    = [self.idx_cter_mic_name]
		self.sort_map_list[self.idx_sort_def]          = [self.idx_cter_def]
		self.sort_map_list[self.idx_sort_total_ac]     = [self.idx_cter_total_ac]
		self.sort_map_list[self.idx_sort_astig_amp]    = [self.idx_cter_astig_amp]
		self.sort_map_list[self.idx_sort_astig_ang]    = [self.idx_cter_astig_ang]
###		self.sort_map_list[self.idx_sort_sd_def]       = [self.idx_cter_sd_def]
		self.sort_map_list[self.idx_sort_cv_def]       = [self.idx_cter_cv_def]
		self.sort_map_list[self.idx_sort_sd_total_ac]  = [self.idx_cter_sd_total_ac]
###		self.sort_map_list[self.idx_sort_sd_astig_amp] = [self.idx_cter_sd_astig_amp]
		self.sort_map_list[self.idx_sort_cv_astig_amp] = [self.idx_cter_cv_astig_amp]
		self.sort_map_list[self.idx_sort_sd_astig_ang] = [self.idx_cter_sd_astig_ang]
		self.sort_map_list[self.idx_sort_error_def]    = [self.idx_cter_error_def]
		self.sort_map_list[self.idx_sort_error_astig]  = [self.idx_cter_error_astig]
		self.sort_map_list[self.idx_sort_error_ctf]    = [self.idx_cter_error_ctf]
		self.sort_map_list[self.idx_sort_max_freq]     = [self.idx_cter_max_freq]
		self.sort_map_list[self.idx_sort_reserved]     = [self.idx_cter_reserved]
		self.sort_map_list[self.idx_sort_phase_shift]  = [self.idx_cter_phase_shift]
###		if self.is_enable_max_power == True: 
###			self.sort_map_list[self.idx_sort_max_power] = [self.idx_cter_max_power] # MRK_TEST:
		
	def enumerate_histograms(self):
		# Define enumerators for mapping of histogram items (combo box widget) and threshold setting (line edit widgets)
		i_enum = -1
		i_enum += 1; self.idx_hist_def          = i_enum
		i_enum += 1; self.idx_hist_total_ac     = i_enum
		i_enum += 1; self.idx_hist_astig_amp    = i_enum
		i_enum += 1; self.idx_hist_astig_ang    = i_enum
###		i_enum += 1; self.idx_hist_sd_def       = i_enum
		i_enum += 1; self.idx_hist_cv_def       = i_enum
		i_enum += 1; self.idx_hist_sd_total_ac  = i_enum
###		i_enum += 1; self.idx_hist_sd_astig_amp = i_enum
		i_enum += 1; self.idx_hist_cv_astig_amp = i_enum
		i_enum += 1; self.idx_hist_sd_astig_ang = i_enum
		i_enum += 1; self.idx_hist_error_def    = i_enum
		i_enum += 1; self.idx_hist_error_astig  = i_enum
		i_enum += 1; self.idx_hist_error_ctf    = i_enum
		i_enum += 1; self.idx_hist_max_freq     = i_enum
		i_enum += 1; self.idx_hist_reserved     = i_enum
		i_enum += 1; self.idx_hist_phase_shift  = i_enum
###		if self.is_enable_max_power == True: 
###			i_enum += 1; self.idx_hist_max_power = i_enum  # MRK_TEST:
		i_enum += 1; self.n_idx_hist            = i_enum
		
	def enumerate_histogram_pulldowns(self):
		# Define enumerators for items of each entry in hist_map_list
		i_enum = -1
		i_enum += 1; self.idx_hist_item_idx_cter                = i_enum
		i_enum += 1; self.idx_hist_item_idx_sort                = i_enum
		i_enum += 1; self.idx_hist_item_val_min                 = i_enum
		i_enum += 1; self.idx_hist_item_val_max                 = i_enum
		i_enum += 1; self.idx_hist_item_unapply_threshold_lower = i_enum
		i_enum += 1; self.idx_hist_item_unapply_threshold_upper = i_enum
		i_enum += 1; self.idx_hist_item_unapply_widget_lower    = i_enum
		i_enum += 1; self.idx_hist_item_unapply_widget_upper    = i_enum
		i_enum += 1; self.idx_hist_item_apply_threshold_lower   = i_enum
		i_enum += 1; self.idx_hist_item_apply_threshold_upper   = i_enum
		i_enum += 1; self.idx_hist_item_apply_widget_lower      = i_enum
		i_enum += 1; self.idx_hist_item_apply_widget_upper      = i_enum
		i_enum += 1; self.n_idx_hist_item                       = i_enum
		
	def map_histogram_items(self):
		# Map histogram items (combo box widget) and threshold setting (line edit widgets)
		# Includes mapping from idx_hist to idx_cter and idx_sort
		self.hist_map_list = [None] * self.n_idx_hist
		self.hist_map_list[self.idx_hist_def]          = [self.idx_cter_def, self.idx_sort_def, 0, 5, 0, 5, None, None, 0, 5, None, None]
		self.hist_map_list[self.idx_hist_total_ac]     = [self.idx_cter_total_ac, self.idx_sort_total_ac, 0, 100, 0, 100, None, None, 0, 100, None, None]
		self.hist_map_list[self.idx_hist_astig_amp]    = [self.idx_cter_astig_amp, self.idx_sort_astig_amp, 0, 1, 0, 1, None, None, 0, 1, None, None]
		self.hist_map_list[self.idx_hist_astig_ang]    = [self.idx_cter_astig_ang, self.idx_sort_astig_ang, 0, 180, 0, 180, None, None, 0, 180, None, None]
###		self.hist_map_list[self.idx_hist_sd_def]       = [self.idx_cter_sd_def, self.idx_sort_sd_def, 0, 5, 0, 5, None, None, 0, 5, None, None]
		self.hist_map_list[self.idx_hist_cv_def]       = [self.idx_cter_cv_def, self.idx_sort_cv_def, 0, 5, 0, 5, None, None, 0, 5, None, None]
		self.hist_map_list[self.idx_hist_sd_total_ac]  = [self.idx_cter_sd_total_ac, self.idx_sort_sd_total_ac, 0, 100, 0, 100, None, None, 0, 100, None, None]
###		self.hist_map_list[self.idx_hist_sd_astig_amp] = [self.idx_cter_sd_astig_amp, self.idx_sort_sd_astig_amp, 0, 1, 0, 1, None, None, 0, 1, None, None]
		self.hist_map_list[self.idx_hist_cv_astig_amp] = [self.idx_cter_cv_astig_amp, self.idx_sort_cv_astig_amp, 0, 1, 0, 1, None, None, 0, 1, None, None]
		self.hist_map_list[self.idx_hist_sd_astig_ang] = [self.idx_cter_sd_astig_ang, self.idx_sort_sd_astig_ang, 0, 180, 0, 180, None, None, 0, 180, None, None]
		self.hist_map_list[self.idx_hist_error_def]    = [self.idx_cter_error_def, self.idx_sort_error_def, 0, 10, 0, 10, None, None, 0, 10, None, None]
		self.hist_map_list[self.idx_hist_error_astig]  = [self.idx_cter_error_astig, self.idx_sort_error_astig, 0, 10, 0, 10, None, None, 0, 10, None, None]
		self.hist_map_list[self.idx_hist_error_ctf]    = [self.idx_cter_error_ctf, self.idx_sort_error_ctf, 0, 10, 0, 10, None, None, 0, 10, None, None]
		self.hist_map_list[self.idx_hist_max_freq]     = [self.idx_cter_max_freq, self.idx_sort_max_freq, 0, 10, 0, 10, None, None, 0, 10, None, None]
		self.hist_map_list[self.idx_hist_reserved]     = [self.idx_cter_reserved, self.idx_sort_reserved, 0, 10, 0, 10, None, None, 0, 10, None, None]
		self.hist_map_list[self.idx_hist_phase_shift]  = [self.idx_cter_phase_shift, self.idx_sort_phase_shift, 0, 180, 0, 180, None, None, 0, 180, None, None]
###		if self.is_enable_max_power == True: 
###			sself.hist_map_list[self.idx_hist_max_power] = [self.idx_cter_max_power, self.idx_sort_max_power, 0, 99999, 0, 99999, None, None, 0, 99999, None, None] # MRK_TEST:
		
	def enumerate_threshold_controls(self):
		# Define enumerators for threshold control selection
		i_enum = -1
		i_enum += 1; self.idx_threshold_control_lower     = i_enum
		i_enum += 1; self.idx_threshold_control_upper     = i_enum
		i_enum += 1; self.idx_threshold_control_edit_only = i_enum
		i_enum += 1; self.n_idx_threshold_control         = i_enum
		
	def enumerate_threshold_map_indices(self):
		# Define enumerators for items of each entry in threshold_control_map_list
		i_enum = -1
		i_enum += 1; self.idx_threshold_control_item_label = i_enum
		i_enum += 1; self.idx_threshold_control_item_color = i_enum
		i_enum += 1; self.n_idx_threshold_control_item     = i_enum
		
	def map_threshold_control(self):
		# Mapping for threshold control (combo box widget)
		self.threshold_control_map_list = [None] * self.n_idx_threshold_control
		self.threshold_control_map_list[self.idx_threshold_control_lower]     = ["Lower (blue)", "blue"]
		self.threshold_control_map_list[self.idx_threshold_control_upper]     = ["Upper (red)", "red"]
		self.threshold_control_map_list[self.idx_threshold_control_edit_only] = ["Edit Only", "black"]
		
	def enumerate_1drot_graphs(self):
		# Define enumerators for display curve selection of 1D power spectrum & CTF fitting plot
		i_enum = -1
		i_enum += 1; self.idx_graph_exp_no_astig   = i_enum
		i_enum += 1; self.idx_graph_fit_no_astig   = i_enum
		i_enum += 1; self.idx_graph_exp_with_astig = i_enum
		i_enum += 1; self.idx_graph_fit_with_astig = i_enum
		i_enum += 1; self.idx_graph_exp_background = i_enum
		i_enum += 1; self.idx_graph_fit_envelope   = i_enum
		i_enum += 1; self.n_idx_graph              = i_enum
		
	def enumerate_graph_items(self):
		# Define enumerators for items of each entry in graph_map_list
		i_enum = -1
		i_enum += 1; self.idx_graph_item_name   = i_enum
		i_enum += 1; self.idx_graph_item_label  = i_enum
		i_enum += 1; self.idx_graph_idx_rotinf  = i_enum
		i_enum += 1; self.idx_graph_item_widget = i_enum
		i_enum += 1; self.n_idx_graph_item      = i_enum
		
	def map_display_checkboxes(self):
		# Map for graph display setting (check box widgets)
		self.graph_map_list = [None] * self.n_idx_graph
		self.graph_map_list[self.idx_graph_exp_no_astig]   = ["exp_no_astig",   "Exp. No Astig (Black)",   self.idx_rotinf_exp_no_astig,   None]
		self.graph_map_list[self.idx_graph_fit_no_astig]   = ["fit_no_astig",   "Fit. No Astig (Blue)",    self.idx_rotinf_fit_no_astig,   None]
		self.graph_map_list[self.idx_graph_exp_with_astig] = ["exp_with_astig", "Exp. with Astig (Red)",   self.idx_rotinf_exp_with_astig, None]
		self.graph_map_list[self.idx_graph_fit_with_astig] = ["fit_with_astig", "Fit. with Astig (Green)", self.idx_rotinf_fit_with_astig, None]
		self.graph_map_list[self.idx_graph_exp_background] = ["exp_background", "Exp. No Backg. (Olive)",  self.idx_rotinf_exp_background, None]
		self.graph_map_list[self.idx_graph_fit_envelope]   = ["fit_envelope",   "Fit. Envelope (Cyan)",    self.idx_rotinf_fit_envelope,   None]
		
	def enumerate_threshold_status(self):
		# Define enumerators for threshold apply status
		i_enum = -1
		i_enum += 1; self.idx_thresholdset_unapplied = i_enum
		i_enum += 1; self.idx_thresholdset_applied   = i_enum
		i_enum += 1; self.n_idx_thresholdset         = i_enum
		
	def enumerate_threshold_entries(self):
		# Define enumerators for items of each entry in thresholdset_map_list
		i_enum = -1
		i_enum += 1; self.idx_thresholdset_item_label  = i_enum
		i_enum += 1; self.n_idx_thresholdset_item      = i_enum
		
	def map_threshold_status(self):
		# Map for threshold set (combo box widget)
		self.thresholdset_map_list = [None] * self.n_idx_thresholdset
		self.thresholdset_map_list[self.idx_thresholdset_unapplied] = ["Unapplied"]
		self.thresholdset_map_list[self.idx_thresholdset_applied]   = ["Applied"]
		
	def initialize_paths(self):
		self.cter_partres_file_path  = None
		self.cter_entry_list         = None
		self.cter_mic_file_path      = None
		self.cter_micthumb_file_path = None
		self.cter_pwrot_file_path    = None
		self.cter_fft_file_path      = None
		
	def initialize_sorting(self):
		self.curentry = None
		self.cursortidx = 0
		self.cursortorder = False
		self.cursortselect = False
		self.curhistidx = 0
		self.curthresholdcontrol = 0
		self.curentryperbin = 10
		self.cursyncsort = False
		self.curthresholdset = 0
	
	def set_default_display(self):
		self.curplotrotavgdisplay = False  # True  # (now off by default)
		self.curplotrotzoomdisplay = True
		self.curimgmicthumbdisplay = False  # True  # (now off by default)
		#self.curhistdisable = False
		self.curhistogramdisplay = True
		self.curscatterdisplay = True
		self.curplotfixscale = 1.1  # 5  (applied envelope and subtracted background -- can still override from GUI)
		self.curfftdisplay = False
		
	def initialize_popup_windows(self):
		# NOTE: 2016/03/09 Toshio Moriya
		# To set window flags of EMGLWidget (SXPlot2DWidget and EMImage2DWidget) window,
		# we have to go through its qt_parent attribute to call setWindowTitle()...
		# 
		self.wfft=EMImage2DWidget()
		self.wfft.setWindowTitle("sxgui_cter - 2D FFT")
		self.wfft.mmode="app"  # NOT USED?
		#self.wfft.qt_parent.setWindowFlags((self.qt_parent.wfft.windowFlags()| Qt.CustomizeWindowHint) & ~Qt.WindowMinimizeButtonHint) # Disabled minimize icon button in window title bar
		self.wfft.qt_parent.setWindowFlags((self.wfft.qt_parent.windowFlags()| Qt.CustomizeWindowHint) & ~Qt.WindowMinimizeButtonHint) # Disabled minimize icon button in window title bar
		self.is_wfft_minimized = False
		
		self.wimgmicthumb=EMImage2DWidget()
		self.wimgmicthumb.setWindowTitle("sxgui_cter - Micrograph Thumbnail")
		self.wimgmicthumb.mmode="app"  # NOT USED?
		self.wimgmicthumb.qt_parent.setWindowFlags((self.wimgmicthumb.qt_parent.windowFlags()| Qt.CustomizeWindowHint) & ~Qt.WindowMinimizeButtonHint) # Disabled minimize icon button in window title bar
		self.is_wimgmicthumb_minimized = False
		
		self.wplotrotavgcoarse=SXPlot2DWidget()
		self.wplotrotavgcoarse.setWindowTitle("sxgui_cter - Plot")
		self.wplotrotavgcoarse.qt_parent.setWindowFlags((self.wplotrotavgcoarse.qt_parent.windowFlags()| Qt.CustomizeWindowHint) & ~Qt.WindowMinimizeButtonHint) # Disabled minimize icon button in window title bar
		self.is_wplotrotavgcoarse_minimized = False
		
		self.wplotrotavgfine=SXPlot2DWidget()
		self.wplotrotavgfine.setWindowTitle("sxgui_cter - Plot Zoom")
		self.wplotrotavgfine.qt_parent.setWindowFlags((self.wplotrotavgfine.qt_parent.windowFlags()| Qt.CustomizeWindowHint) & ~Qt.WindowMinimizeButtonHint) # Disabled minimize icon button in window title bar
		self.is_wplotrotavgfine_minimized = False
		
		self.whistparam=SXPlot2DWidget()
		self.whistparam.setWindowTitle("sxgui_cter - Histogram")
		self.whistparam.qt_parent.setWindowFlags((self.whistparam.qt_parent.windowFlags()| Qt.CustomizeWindowHint) & ~Qt.WindowMinimizeButtonHint) # Disabled minimize icon button in window title bar
		self.is_whistparam_minimized = False
		
		self.wscatterparam=SXPlot2DWidget()
		self.wscatterparam.setWindowTitle("sxgui_cter - Sort Plot")
		self.wscatterparam.qt_parent.setWindowFlags((self.wscatterparam.qt_parent.windowFlags()| Qt.CustomizeWindowHint) & ~Qt.WindowMinimizeButtonHint) # Disabled minimize icon button in window title bar
		self.is_wscatterparam_minimized = False
	
	def draw_main_window(self):
		# Place layout inside QWidget
		templayout = QtGui.QHBoxLayout(self)
		templayout.setContentsMargins(0,0,0,0)
		mainwidget = QtGui.QWidget(self)
		mainwidget.setObjectName("MainWidgetObject")
		
		# Color scheme
		background_image_file_path = '{0}sxgui_background.png'.format(get_image_directory())
		mainwidget.setStyleSheet("QWidget#MainWidgetObject {{background-image: url('{0}')}}".format(background_image_file_path))
		mainlayout = QtGui.QHBoxLayout(mainwidget)
		mainlayout.setContentsMargins(12,12,12,12)
		templayout.addWidget(mainwidget)
		
		# --------------------------------------------------------------------------------
		# Columns 1-3
		# --------------------------------------------------------------------------------
		
		leftwidget = QtGui.QWidget(self)
		leftwidget.setObjectName("LeftWidgetObject")
		leftwidget.setStyleSheet("QWidget#LeftWidgetObject {background-color: rgba(229, 229, 229, 208); border-radius: 15px;}")
		
		leftcolumn = QtGui.QVBoxLayout(leftwidget)
		leftcolumn.setContentsMargins(0,10,0,10)
		mainlayout.addWidget(leftwidget)
		
		labelwidth = 90
		editwidth = 100
		sublabelwidth = 140
		
		self.pbopencter = QtGui.QPushButton("Open CTER partres file")
		#self.pbopencter.setStyleSheet("QPushButton {color:gray; }")  # (This doesn't do anything.)
		leftcolumn.addWidget(self.pbopencter)
		
		self.add_centered_label("<b>Selection Summary:</b>", leftcolumn)
		
		self.vbnentry = ValBox(self,(0,10000),None,0)
		self.add_label_with_value("Num. of entries", self.vbnentry, leftcolumn)
		
		self.vbuncheckcounts = ValBox(self,(0,1000000),None,0)
		self.add_label_with_value("Unchecked", self.vbuncheckcounts, leftcolumn, style_sheet="color: rgb(0,0,0);")
		
		self.vbuncheckratio = ValBox(self,(0,1.0),None,0)
		self.add_label_with_value("Ratio", self.vbuncheckratio, leftcolumn, style_sheet="color: rgb(0,0,0);", intonly=False)

		self.add_centered_label("", leftcolumn)  # spacer
		self.add_centered_label("<b>Electron Microscopy:</b>", leftcolumn)
		
		# Voltage
		self.add_label_with_cter_param(self.idx_cter_vol, leftcolumn, 0, 500, style_sheet = "color: rgb(127,127,127);", labelwidth=labelwidth)
		
		# Spherical aberration
		self.add_label_with_cter_param(self.idx_cter_cs, leftcolumn, 0, 5, style_sheet = "color: rgb(127,127,127);",labelwidth=labelwidth)
		
		# Pixel size
		self.add_label_with_cter_param(self.idx_cter_apix, leftcolumn, 0,500, style_sheet = "color: rgb(127,127,127);",labelwidth=labelwidth)
		
		self.add_centered_label("", leftcolumn)  # spacer
		self.add_centered_label("<b>Display Windows:</b>", leftcolumn)
		
		self.cbrotavgdisplay = CheckBox(None,None,self.curplotrotavgdisplay)
		self.add_label_with_checkbox("Rot. Avg. Plot", self.cbrotavgdisplay, leftcolumn, labelwidth=sublabelwidth)
		
		self.cbrotzoomdisplay = CheckBox(None,None,self.curplotrotzoomdisplay)
		self.add_label_with_checkbox("Rot. Avg. Plot Zoom", self.cbrotzoomdisplay, leftcolumn, labelwidth=sublabelwidth)
		
		self.cbhistogramdisplay = CheckBox(None,None,self.curhistogramdisplay)
		self.add_label_with_checkbox("Histogram", self.cbhistogramdisplay, leftcolumn, labelwidth=sublabelwidth)
		
		self.cbscatterdisplay = CheckBox(None,None,self.curscatterdisplay)
		self.add_label_with_checkbox("Sort Plot", self.cbscatterdisplay, leftcolumn, labelwidth=sublabelwidth)
		
		self.cbmicthumbdisplay = CheckBox(None,None,self.curimgmicthumbdisplay)
		self.add_label_with_checkbox("Micrograph Thumbnail", self.cbmicthumbdisplay, leftcolumn, labelwidth=sublabelwidth)
		
		self.cbfftdisplay = CheckBox(None,None,self.curfftdisplay)
		self.add_label_with_checkbox("2D Power Spectrum", self.cbfftdisplay, leftcolumn, labelwidth=sublabelwidth)
		
		self.add_centered_label("", leftcolumn)  # spacer
		self.add_centered_label("<b>Display Curves:</b>", leftcolumn)
		
		for idx_graph in range(self.n_idx_graph):
			self.graph_map_list[idx_graph][self.idx_graph_item_widget] = CheckBox(None,None,True)
			self.add_label_with_checkbox(self.graph_map_list[idx_graph][self.idx_graph_item_label], 
					self.graph_map_list[idx_graph][self.idx_graph_item_widget], leftcolumn, labelwidth=sublabelwidth)
		
		self.vbplotfixscale = ValBox(self,(0,99999),None,self.curplotfixscale)  # default <- self.curplotfixscale
		self.add_label_with_value("Plot Fix Scale", self.vbplotfixscale, leftcolumn, labelwidth=sublabelwidth, 
							enabled=True, intonly=False, style_sheet="color: rgb(0,0,0);")
		
		self.add_centered_label("", leftcolumn)  # spacer
		
		self.pbrefreshgraphs = QtGui.QPushButton("Refresh Graphs")
		self.pbrefreshgraphs.setEnabled(False)
		leftcolumn.addWidget(self.pbrefreshgraphs)
		
		leftcolumn.addStretch(1)
		
		# --------------------------------------------------------------------------------
		# 2nd column
		# --------------------------------------------------------------------------------
		
		secondcolumn = QtGui.QVBoxLayout()
		secondcolumn.setContentsMargins(0,0,0,0)
		mainlayout.addLayout(secondcolumn)
		
		# plot list and plot mode combobox
		row_span_entry_list = 27  # length of file list (QListWidget)
		self.lbentry = SXListWidget(self)
		self.lbentry.setSizePolicy(QtGui.QSizePolicy.Preferred,QtGui.QSizePolicy.Expanding)
		self.lbentry.setMinimumWidth(220)
		secondcolumn.addWidget(self.lbentry)
		
		labelwidth = 180
		editwidth = 100
		sublabelwidth = editwidth
		sidemargin = 10
		topmargin = 3
		bottommargin = 3
		borderwidth = 1
		
		# --------------------------------------------------------------------------------
		# The main big layout (VBox), I'll split up into 4 HBoxLayouts: upper, threshold, bottom, menu
		# --------------------------------------------------------------------------------

		biglayout = QtGui.QVBoxLayout()
		biglayout.setContentsMargins(0,0,0,0)
		mainlayout.addLayout(biglayout)
		
		cterwidget = QtGui.QWidget(self)
		cterwidget.setObjectName("CterWidgetObject")
		cterwidget.setStyleSheet("QWidget#CterWidgetObject {background-color: rgba(229, 229, 229, 208); border-radius: 15px;}")
		biglayout.addWidget(cterwidget)
		
		chartlayout = QtGui.QVBoxLayout(cterwidget)
		chartlayout.setContentsMargins(0,10,10,10)
		
		upperlayout = QtGui.QHBoxLayout()
		upperlayout.setContentsMargins(0,0,0,0)
		chartlayout.addLayout(upperlayout)
		
		first2rowslayout = QtGui.QVBoxLayout()
		first2rowslayout.setContentsMargins(0,0,0,0)
		upperlayout.addLayout(first2rowslayout)
		
		self.add_centered_label("<b>Current Entry Info:</b>", first2rowslayout, labelwidth=350)  # hardwired labelwidth
		
		self.ssortedid = ValBox(self,(0,10000),None,0)
		self.add_label_with_value("Sorted ID", self.ssortedid, first2rowslayout, style_sheet="color: rgb(0,0,0);", labelwidth=labelwidth, editwidth=editwidth)
		
		self.add_label_with_cter_param(self.idx_cter_id, first2rowslayout, 0,10000, intonly=True, labelwidth=labelwidth, editwidth=editwidth)
		
		# Add image
		logolayout = QtGui.QVBoxLayout()
		logolayout.setContentsMargins(0,0,0,0)
		upperlayout.addLayout(logolayout)
		
		logo = SXLogoButton("sxgui_pictograph_cter.png", 64, parent=self)
		logo.add_sxmenu_item_btn_widget(logolayout)
		
		# Layout for thresholds and column labels
		threshlayout = QtGui.QHBoxLayout()
		threshlayout.maximumSize()
		chartlayout.addLayout(threshlayout)
		
		# Draw borderless box to preserve spacing
		cterFrame = QtGui.QWidget(self)
		cterFrame.setContentsMargins(0,0,0,0)
		cterFrame.setStyleSheet("border: %spx solid transparent;" % borderwidth)
		threshlayout.addWidget(cterFrame)
		
		# Layout for CTER columns: select through phase shift
		cterlayout = QtGui.QVBoxLayout(cterFrame)
		cterlayout.setContentsMargins(0,0,0,bottommargin)
		
		# Selection flag
		self.add_label_with_cter_param(self.idx_cter_select, cterlayout, 0, 1, intonly=True, labelwidth=labelwidth, editwidth=editwidth)

		# Draw boxes around limits
		unappliedLimitFrame = QtGui.QWidget(self)  # QGridLayout doesn't have borders, so I'm enclosing it inside a QFrame, which can
		unappliedLimitFrame.setContentsMargins(0,0,0,0)
		unappliedLimitFrame.setStyleSheet("border: %spx solid rgb(0,0,0);" % borderwidth)
		threshlayout.addWidget(unappliedLimitFrame)
		
		# Layout for unapplied threshholds
		unappliedlayout = QtGui.QVBoxLayout(unappliedLimitFrame)
		unappliedlayout.setObjectName("unappliedlayout")
		unappliedlayout.setContentsMargins(sidemargin,topmargin,sidemargin,bottommargin)
		
		self.add_centered_label("<b>Unapplied Thresholds:</b>", unappliedlayout, style_sheet="border: 0px;")
		
		# Applied limits
		appliedLimitFrame = QtGui.QWidget(self)  # QGridLayout doesn't have borders, so I'm enclosing it inside a QFrame, which can
		appliedLimitFrame.setStyleSheet("border: %spx solid rgb(0,0,0);" % borderwidth)
		threshlayout.addWidget(appliedLimitFrame)
		
		# Layout for applied threshholds
		appliedlayout = QtGui.QVBoxLayout(appliedLimitFrame)
		appliedlayout.setObjectName("appliedlayout")
		appliedlayout.setContentsMargins(sidemargin,topmargin,sidemargin,bottommargin)
		
		self.add_centered_label("<b>Applied Thresholds:</b>", appliedlayout, style_sheet="border: 0px;")
		
		# Write table of thresholds 
		for idx_hist in range(self.n_idx_hist):
			#print("MRK_DEBUG", idx_hist, self.hist_map_list[idx_hist][self.idx_hist_item_idx_cter], self.value_map_list[idx_cter][0], self.hist_map_list[idx_hist][self.idx_hist_item_apply_threshold_lower], self.hist_map_list[idx_hist][self.idx_hist_item_apply_threshold_upper])
			self.add_value_widget_with_label(idx_hist, cterlayout, labelwidth=labelwidth)#, editwidth=5)
			self.add_widgets_unapplied_threshold(idx_hist, unappliedlayout, labelwidth=labelwidth, editwidth=editwidth)
			self.add_widgets_applied_threshold(  idx_hist, appliedlayout, labelwidth=labelwidth, editwidth=editwidth)
			
		# Last two CTER rows
		last2rowFrame = QtGui.QWidget(self)
		last2rowFrame.setContentsMargins(0,0,0,0)
		last2rowFrame.setStyleSheet("border: %spx solid transparent;" % borderwidth)
		chartlayout.addWidget(last2rowFrame)
		last2rowlayout = QtGui.QVBoxLayout(last2rowFrame)
		last2rowlayout.setContentsMargins(0,0,0,0)
		
		# Amplitude contrast and B-factor (These used to be grayed out.)
		self.add_label_with_cter_param(self.idx_cter_const_ac, last2rowlayout, 0,1600, 
								 labelwidth=labelwidth, editwidth=editwidth, enabled=False, style_sheet = "color: rgb(127,127,127);")
		self.add_label_with_cter_param(self.idx_cter_bfactor, last2rowlayout, 0,1600, 
								 labelwidth=labelwidth, editwidth=editwidth, enabled=False, style_sheet = "color: rgb(127,127,127);")
		
		# ---------------------------------
		# Settings layout
		# ---------------------------------
		
		settingswidget = QtGui.QWidget(self)
		settingswidget.setObjectName("SettingsWidgetObject")
		settingswidget.setStyleSheet("QWidget#SettingsWidgetObject {background-color: rgba(229, 229, 229, 208); border-radius: 15px;}")
		
		settingslayout = QtGui.QHBoxLayout(settingswidget)
		#biglayout.addLayout(settingslayout)
		biglayout.addWidget(settingswidget)
		
		# There will be three VBox layouts: sort, histogram/plot, and save/load
		sortlayout = QtGui.QVBoxLayout()
		sortlayout.setContentsMargins(0,0,0,0)
		settingslayout.addLayout(sortlayout)
		
		self.add_centered_label("<b>Sort CTER Partres Entries:</b>", sortlayout)
		
		# ---------------------------------
		# Pulldown menu options for sorting
		# ---------------------------------
		
		self.ssort=QtGui.QComboBox(self)
		for map_entry in self.sort_map_list:
			idx_cter = map_entry[self.idx_sort_item_idx_cter]
			self.ssort.addItem(self.value_map_list[idx_cter][self.idx_cter_item_label])
		self.ssort.setCurrentIndex(self.cursortidx)
		sortlayout.addWidget(self.ssort)
		
		# Checkbox to reverse order
		self.cbsortorder = CheckBox(None,None,self.cursortorder)
		self.add_label_with_checkbox("Descending", self.cbsortorder, sortlayout, labelwidth=sublabelwidth)
		
		# Checkbox to list deselected files first
		self.cbsortselect = CheckBox(None,None,self.cursortselect)
		self.add_label_with_checkbox("Sort Select", self.cbsortselect, sortlayout, labelwidth=sublabelwidth)
		
		sortlayout.addStretch(1)
		#self.add_centered_label("", sortlayout)  # spacer
		
		self.pbreapplysort = QtGui.QPushButton("Reapply Sort")
		self.pbreapplysort.setEnabled(False)
		sortlayout.addWidget(self.pbreapplysort)
		
		# ---------------------------------
		# Histogram & Plot Settings
		# ---------------------------------

		histplotlayout = QtGui.QVBoxLayout()
		#histplotlayout.setContentsMargins(0,3,0,0)  # I don't know why these are needed to align the layouts in the HBox
		histplotlayout.setContentsMargins(0,0,0,0)  # I don't know why these are needed to align the layouts in the HBox
		settingslayout.addLayout(histplotlayout)
		
		self.add_centered_label("<b>Histogram & Plot Settings:</b>", histplotlayout)
		
		# Pulldown menu options for histogram
		self.shist = QtGui.QComboBox(self)
		self.shist.setMaximumWidth(250)
		for map_entry in self.hist_map_list:
			idx_cter = map_entry[self.idx_hist_item_idx_cter]
			self.shist.addItem(self.value_map_list[idx_cter][self.idx_cter_item_label])
		self.shist.setCurrentIndex(self.curhistidx)
		histplotlayout.addWidget(self.shist)
		
		# Pulldown menu to move lower/upper threshold
		self.add_label_with_pulldown("Move Threshold", histplotlayout, labelwidth=100, menuwidth=140)
		
		# Checkbox to sort according to histogrammed/plotted values
		self.cbsyncsort = CheckBox(None,None,self.cursyncsort)
		self.add_label_with_checkbox("Sync. w/Sort", self.cbsyncsort, histplotlayout, labelwidth=95)
		
		self.vsentryperbin = self.add_label_with_param("counts/bin", histplotlayout, 0,10000, labelwidth=95, intonly=True)
		
		histplotlayout.addStretch(1)
		
		self.pbapplyallthreshold = QtGui.QPushButton("Apply All Thresholds")
		self.pbapplyallthreshold.setMaximumWidth(250)
		self.pbapplyallthreshold.setEnabled(False)
		histplotlayout.addWidget(self.pbapplyallthreshold)
		
		# ---------------------------------
		# Save/Load thresholds
		# ---------------------------------

		saveloadlayout = QtGui.QVBoxLayout()
		saveloadlayout.setContentsMargins(0,0,0,0)
		saveloadlayout.setSpacing(8.5)
		settingslayout.addLayout(saveloadlayout)
		
		self.add_centered_label("<b>Save/Load Thresholds:</b>", saveloadlayout)
		
		# Pulldown menu to save unapplied/applied thresholds
		self.sthresholdset=QtGui.QComboBox(self)
		for map_entry in self.thresholdset_map_list:
			self.sthresholdset.addItem(map_entry[self.idx_thresholdset_item_label])
		self.sthresholdset.setCurrentIndex(self.curthresholdset)
		saveloadlayout.addWidget(self.sthresholdset)
		
		# Save/Load threshold buttons
		self.pbsavethresholdset = QtGui.QPushButton("Save")
		self.pbloadthresholdset = QtGui.QPushButton("Load")
		self.add_two_buttons(self.pbsavethresholdset, self.pbloadthresholdset, saveloadlayout)
		
		self.add_centered_label("<b>Save Selection:</b>", saveloadlayout)
		
		# Prefix for output files
		self.vfilesuffix = StringBox(self,None,"Trial00")
		self.add_label_with_textbox("File Suffix", self.vfilesuffix, saveloadlayout)
		
		saveloadlayout.addStretch(1)
		
		self.pbsaveselection = QtGui.QPushButton("Save Selection")
		self.pbsaveselection.setEnabled(False)
		saveloadlayout.addWidget(self.pbsaveselection)
		
		# Pads at the bottom
		biglayout.addStretch(1)
		
		self.setWindowTitle("sxgui_cter - Control Panel")
		
	def signal_handler(self):
		self.pbopencter.clicked[bool].connect(self.openCterPartres)
		
		self.cbrotavgdisplay.valueChanged.connect(self.newRotAvgDisplay)
		self.cbrotzoomdisplay.valueChanged.connect(self.newRotZoomDisplay)
		self.cbmicthumbdisplay.valueChanged.connect(self.newMicThumbDisplay)
		self.cbhistogramdisplay.valueChanged.connect(self.newHistogramDisplay)
		self.cbscatterdisplay.valueChanged.connect(self.newScatterDisplay)
		self.cbfftdisplay.valueChanged.connect(self.newFFTDisplay)
		
		for idx_graph in range(self.n_idx_graph):
			self.graph_map_list[idx_graph][self.idx_graph_item_widget].valueChanged.connect(self.updatePlotVisibility)
		self.vbplotfixscale.valueChanged.connect(self.newPlotFixScale)
		self.pbrefreshgraphs.clicked[bool].connect(self.refreshGraphs)
		
		self.lbentry.currentRowChanged[int].connect(self.newEntry)
#		QtCore.QObject.connect(self.lbentry,QtCore.SIGNAL("keypress"),self.entryKey)
		self.lbentry.itemChanged.connect(self.updateEntrySelect)
		
		self.ssort.currentIndexChanged[int].connect(self.newSort)
		self.cbsortorder.valueChanged.connect(self.newSortOrder)
		self.cbsortselect.valueChanged.connect(self.newSortSelect)
		self.pbreapplysort.clicked[bool].connect(self.reapplySort)
		
		#for idx_hist in range(self.n_idx_hist):
			#self.hist_map_list[idx_hist][self.idx_hist_item_unapply_widget_lower].valueChanged.connect(self.newThresholdLower)
			#self.hist_map_list[idx_hist][self.idx_hist_item_unapply_widget_upper].valueChanged.connect(self.newThresholdUpper)
		
		self.shist.currentIndexChanged[int].connect(self.newHistogramRow)
		self.sthresholdcontrol.currentIndexChanged[int].connect(self.newThresholdControl)
		self.cbsyncsort.valueChanged.connect(self.newSyncSort)
		self.vsentryperbin.valueChanged.connect(self.newEntryPerBin)
		self.pbapplyallthreshold.clicked[bool].connect(self.applyAllThresholds)
		
		self.sthresholdset.currentIndexChanged[int].connect(self.newThresholdSet)
		self.pbsavethresholdset.clicked[bool].connect(self.saveThresholdSet)
		self.pbloadthresholdset.clicked[bool].connect(self.loadThresholdSet)
		
		self.pbsaveselection.clicked[bool].connect(self.saveSelection)
		
	def set_size_popup_windows(self):
		# Set default sizes & positions of windows in case this is the first time to run in this project directory
		# (I figured these values out by printing the width and height in resize event)
		main_win_width = 1166
		main_win_height = 726
		child_win_width = 1166
		child_win_height = 726
		img_size = 512
		
		win_left = 0
		win_top = 0
		win_left_shift = 30
		win_top_shift = 30
		
		self.resize(main_win_width, main_win_height)
		self.move(win_left, win_top); win_left += win_left_shift; win_top += win_top_shift
		self.wscatterparam.qt_parent.resize(child_win_width,child_win_height)
		self.wscatterparam.qt_parent.move(win_left, win_top); win_left += win_left_shift; win_top += win_top_shift
		self.whistparam.qt_parent.resize(child_win_width,child_win_height)
		self.whistparam.qt_parent.move(win_left, win_top); win_left += win_left_shift; win_top += win_top_shift
		self.wplotrotavgfine.qt_parent.resize(child_win_width,child_win_height)
		self.wplotrotavgfine.qt_parent.move(win_left, win_top); win_left += win_left_shift; win_top += win_top_shift
		self.wplotrotavgcoarse.qt_parent.resize(child_win_width,child_win_height)
		self.wplotrotavgcoarse.qt_parent.move(win_left, win_top); win_left += win_left_shift; win_top += win_top_shift
		self.wimgmicthumb.set_data(model_blank(img_size,img_size, bckg=1.0)) # resize does not work if no image is set
		self.wimgmicthumb.qt_parent.resize(child_win_width,child_win_height)
		self.wimgmicthumb.qt_parent.move(win_left, win_top); win_left += win_left_shift; win_top += win_top_shift
		self.wimgmicthumb.scroll_to(-1 * img_size,-1 * img_size)
		self.wfft.set_data(model_blank(img_size,img_size, bckg=1.0)) # resize does not work if no image is set
		self.wfft.qt_parent.resize(child_win_width,child_win_height)
		self.wfft.qt_parent.move(win_left, win_top); win_left += win_left_shift; win_top += win_top_shift
		self.wfft.scroll_to(-1 * img_size,-1 * img_size)
		
#		# The following are obsolete after 2016/07/04
#		win_height = 512  # Let use the same height for all windows
#		win_height_margin = 46
#		main_win_width = 1200
#		graph_win_width = 980
#		img_win_width = win_height
#		# Top Left
#		win_top = 0
#		win_left = 0
#		win_width = graph_win_width
#		self.whistparam.qt_parent.resize(win_width,win_height)
#		self.whistparam.qt_parent.move(win_left,win_top)
#		self.wscatterparam.qt_parent.resize(win_width,win_height)
#		self.wscatterparam.qt_parent.move(win_left,win_top)
#		# Top Right
#		win_left = graph_win_width
#		win_width = main_win_width;
#		self.resize(win_width,win_height)
#		self.move(win_left,win_top)
#		# Bottom Left
#		win_top = win_height + win_height_margin; 
#		win_left = 0
#		win_width = graph_win_width
#		self.wplotrotavgcoarse.qt_parent.resize(win_width,win_height)
#		self.wplotrotavgcoarse.qt_parent.move(win_left,win_top)
#		self.wplotrotavgfine.qt_parent.resize(win_width,win_height)
#		self.wplotrotavgfine.qt_parent.move(win_left,win_top)
#		# Bottom Right
#		# Set the image window
#		win_left = graph_win_width
#		win_width = img_win_width
#		img_size = 512
#		# scale_factor = float(win_width)/img_size
#		self.wimgmicthumb.set_data(model_blank(img_size,img_size, bckg=1.0)) # resize does not work if no image is set
#		self.wimgmicthumb.qt_parent.resize(win_width,win_height)
#		self.wimgmicthumb.qt_parent.move(win_left,win_top)
#		self.wimgmicthumb.scroll_to(-1 * img_size,-1 * img_size)
#		# self.wimgmicthumb.set_scale(scale_factor)
		
	
	def add_centered_label(self, labeltext, target, labelwidth=None, style_sheet=None):
		temp_label = QtGui.QLabel(labeltext,self)
		temp_label.setAlignment(QtCore.Qt.AlignHCenter|QtCore.Qt.AlignVCenter)
		if labelwidth: temp_label.setMaximumWidth(labelwidth)
		if style_sheet: temp_label.setStyleSheet("border: 0px;")
		target.addWidget(temp_label)
		
	def add_label_with_value(self, labeltext, valueentry, target, 
						  labelwidth=90, enabled=False, intonly=True, style_sheet="color: rgb(127,127,127);", editwidth=80, maxwidth=100):
		# Label
		temp_hbox = QtGui.QHBoxLayout()
		temp_hbox.setContentsMargins(0,0,0,0)
		temp_label = QtGui.QLabel(labeltext, self)
		temp_label.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignVCenter)
		#temp_label.setMinimumSize(QtCore.QSize(labelwidth,20))
		temp_label.setFixedSize(QtCore.QSize(labelwidth,20))
		#temp_label.setMaximumWidth(labelwidth)
		temp_hbox.addWidget(temp_label)
		
		# Value
		valueentry.setEnabled(enabled)
		valueentry.intonly = intonly
		valueentry.text.setStyleSheet(style_sheet)
		valueentry.text.setMinimumSize(QtCore.QSize(editwidth,0))
		valueentry.text.setMaximumWidth(maxwidth)
		temp_hbox.addWidget(valueentry)
		temp_hbox.addStretch(1)
		target.addLayout(temp_hbox, 0)
		
	def add_label_with_param(self, labeltext, target, val_min, val_max, 
						  intonly = False, style_sheet = "color: rgb(0,0,0);", labelwidth=80, editwidth=80, maxwidth=100, enabled=False):
		# Label
		temp_hbox = QtGui.QHBoxLayout()
		temp_hbox.setContentsMargins(0,0,0,0)
		temp_label = QtGui.QLabel(labeltext,self)
		temp_label.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignVCenter)
		temp_label.setFixedSize(QtCore.QSize(labelwidth,20))
		temp_hbox.addWidget(temp_label)
		
		# Value
		val_default = val_min
		valueentry = ValBox(self,(val_min,val_max),None,val_default)
		valueentry.setEnabled(enabled)
		valueentry.intonly = intonly
		valueentry.text.setStyleSheet(style_sheet)
		valueentry.text.setMinimumSize(QtCore.QSize(editwidth,0))
		valueentry.text.setMaximumWidth(maxwidth)
		temp_hbox.addWidget(valueentry)
		temp_hbox.addStretch(1)
		target.addLayout(temp_hbox, 0)

		return valueentry
	
	def add_label_with_cter_param(self, idx_cter, target, val_min, val_max, 
						  intonly = False, style_sheet = "color: rgb(0,0,0);", labelwidth=80, editwidth=80, maxwidth=100, enabled=False):
		
		# Specific for self.value_map_list, calls add_label_with_param()
		labeltext = self.value_map_list[idx_cter][self.idx_cter_item_label]
		self.value_map_list[idx_cter][self.idx_cter_item_widget] = self.add_label_with_param(labeltext, target, val_min, val_max, 
							intonly=intonly, style_sheet=style_sheet, labelwidth=labelwidth, editwidth=editwidth, maxwidth=maxwidth, enabled=enabled)
	
	def add_value_widget_with_label(self, idx_hist, target, intonly = False, labelwidth=80, editwidth=80):
		# Specific for CTER parameters select to phase shift, calls add_label_with_cter_param()
		val_min = self.hist_map_list[idx_hist][self.idx_hist_item_val_min]
		val_max = self.hist_map_list[idx_hist][self.idx_hist_item_val_max]
		
		# Add widget for parameter value
		self.add_label_with_cter_param(self.hist_map_list[idx_hist][self.idx_hist_item_idx_cter], target, val_min, val_max, 
							intonly=intonly, labelwidth=labelwidth, editwidth=editwidth)
		
	def add_label_with_checkbox(self, labeltext, checkbutton, target, labelwidth=140):
		temp_hbox = QtGui.QHBoxLayout()
		temp_hbox.setContentsMargins(0,0,0,0)
		temp_label = QtGui.QLabel(labeltext, self)
		temp_label.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignVCenter)
		temp_label.setFixedSize(QtCore.QSize(labelwidth,20))
		temp_hbox.addWidget(temp_label)
		
		temp_hbox.addWidget(checkbutton)
		target.addLayout(temp_hbox, 0)
	
	def add_widgets_unapplied_threshold(self, idx_hist, target, labelwidth = 80, editwidth = 80):
		val_min = self.hist_map_list[idx_hist][self.idx_hist_item_val_min]
		val_max = self.hist_map_list[idx_hist][self.idx_hist_item_val_max]
		
		# Label
		temp_hbox = QtGui.QHBoxLayout()
		temp_hbox.setContentsMargins(0,0,0,0)
		
		# Add widget for unapplied thresholds
		self.hist_map_list[idx_hist][self.idx_hist_item_unapply_widget_lower] = ValBox(self,(val_min,val_max),None,val_min,labelwidth)
		self.hist_map_list[idx_hist][self.idx_hist_item_unapply_widget_lower].setEnabled(False)
		temp_hbox.addWidget(self.hist_map_list[idx_hist][self.idx_hist_item_unapply_widget_lower])
		
		self.hist_map_list[idx_hist][self.idx_hist_item_unapply_widget_lower].text.setMinimumSize(QtCore.QSize(editwidth,0))
		self.hist_map_list[idx_hist][self.idx_hist_item_unapply_widget_upper] = ValBox(self,(val_min,val_max),None,val_max,labelwidth)
		self.hist_map_list[idx_hist][self.idx_hist_item_unapply_widget_upper].setEnabled(False)
		self.hist_map_list[idx_hist][self.idx_hist_item_unapply_widget_upper].text.setStyleSheet("color: rgb(255,0,0);")
		self.hist_map_list[idx_hist][self.idx_hist_item_unapply_widget_upper].text.setMinimumSize(QtCore.QSize(editwidth,0))
		temp_hbox.addWidget(self.hist_map_list[idx_hist][self.idx_hist_item_unapply_widget_upper])
		target.addLayout(temp_hbox, 0)
		
	def add_widgets_applied_threshold(self, idx_hist, target, labelwidth = 80, editwidth = 80):
		val_min = self.hist_map_list[idx_hist][self.idx_hist_item_val_min]
		val_max = self.hist_map_list[idx_hist][self.idx_hist_item_val_max]
		
		# Label
		temp_hbox = QtGui.QHBoxLayout()
		temp_hbox.setContentsMargins(0,0,0,0)
		
		# Add widget for applied thresholds
		self.hist_map_list[idx_hist][self.idx_hist_item_apply_widget_lower] = ValBox(self,(val_min,val_max),None,val_min,labelwidth)
		self.hist_map_list[idx_hist][self.idx_hist_item_apply_widget_lower].setEnabled(False)
		self.hist_map_list[idx_hist][self.idx_hist_item_apply_widget_lower].text.setStyleSheet("color: rgb(0,0,255);")
		self.hist_map_list[idx_hist][self.idx_hist_item_apply_widget_lower].text.setMinimumSize(QtCore.QSize(editwidth,0))
		temp_hbox.addWidget(self.hist_map_list[idx_hist][self.idx_hist_item_apply_widget_lower])
		
		self.hist_map_list[idx_hist][self.idx_hist_item_apply_widget_upper] = ValBox(self,(val_min,val_max),None,val_max,labelwidth)
		self.hist_map_list[idx_hist][self.idx_hist_item_apply_widget_upper].setEnabled(False)
		self.hist_map_list[idx_hist][self.idx_hist_item_apply_widget_upper].text.setStyleSheet("color: rgb(255,0,0);")
		self.hist_map_list[idx_hist][self.idx_hist_item_apply_widget_upper].text.setMinimumSize(QtCore.QSize(editwidth,0))
		temp_hbox.addWidget(self.hist_map_list[idx_hist][self.idx_hist_item_apply_widget_upper])
		target.addLayout(temp_hbox, 0)
		
	def add_label_with_pulldown(self, labeltext, target, 
						  labelwidth=90, menuwidth=100):
		# Label
		temp_hbox = QtGui.QHBoxLayout()
		temp_hbox.setContentsMargins(0,0,0,0)
		temp_label = QtGui.QLabel(labeltext, self)
		temp_label.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignVCenter)
		temp_label.setFixedSize(QtCore.QSize(labelwidth,20))
		temp_hbox.addWidget(temp_label)
		
		# Pulldown/ComboBox (may want to generalize this someday)
		self.sthresholdcontrol = QtGui.QComboBox(self)
		self.sthresholdcontrol.setMaximumWidth(menuwidth)
		#self.sthresholdcontrol.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignVCenter)
		for idx_threshold_control in range(self.n_idx_threshold_control):
			map_entry = self.threshold_control_map_list[idx_threshold_control]
			self.sthresholdcontrol.addItem(map_entry[self.idx_threshold_control_item_label])
			self.sthresholdcontrol.setItemData(idx_threshold_control, QtGui.QColor(map_entry[self.idx_threshold_control_item_color]), Qt.TextColorRole);
		self.sthresholdcontrol.setCurrentIndex(self.curthresholdcontrol)
		temp_hbox.addWidget(self.sthresholdcontrol)
			
		target.addLayout(temp_hbox, 0)
	
	def add_two_buttons(self, button1, button2, target):
		temp_hbox = QtGui.QHBoxLayout()
		temp_hbox.setContentsMargins(0,0,0,0)
		button1.setEnabled(False)
		temp_hbox.addWidget(button1)
		button2.setEnabled(False)
		temp_hbox.addWidget(button2)
		target.addLayout(temp_hbox, 0)
	
	def add_label_with_textbox(self, labeltext, stringbox, target, labelwidth=140):#, valueentry, target, labelwidth=90):
		# Label
		temp_hbox = QtGui.QHBoxLayout()
		temp_hbox.setContentsMargins(0,0,0,0)
		temp_label = QtGui.QLabel(labeltext, self)
		temp_label.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignVCenter)
		#temp_label.setMinimumSize(QtCore.QSize(labelwidth,20))
		temp_hbox.addWidget(temp_label)
		temp_hbox.addWidget(stringbox)
		target.addLayout(temp_hbox, 0)
	
	def readCterPartresFile(self, file_path):
		"""Read all entries from a CTER partres file into the list box"""
		
		if not os.path.exists(file_path):
			QtGui.QMessageBox.warning(None,"Warning","Cannot find CTER partres file (%s). Please check the file path." % (file_path))
			return
		
		if os.path.basename(file_path).find("partres") == -1:
			QtGui.QMessageBox.warning(None,"Warning","Invalid file name for CTER partres file (%s). The file name must contain \"partres\"." % (file_path))
			return
		
		if file_path[-1*len(".txt"):] != ".txt":
			QtGui.QMessageBox.warning(None,"Warning","Invalid file extension for CTER partres file (%s). The file extension must be \".txt\"." % (file_path))
			return
		
		new_entry_list = read_text_row(file_path)
		if len(new_entry_list) == 0:
			QtGui.QMessageBox.warning(self, "Warning", "Specified CTER partres file (%s) does not contain any entry. Please check the file." % (file_path))
			return
		assert len(new_entry_list) > 0, "MRK_DEBUG"
		
		# NOTE: 2017/03/20 Toshio Moriya
		# The following code is to support the old format of CTER partres file. It should be removed near future
		
		if len(new_entry_list[0]) == self.n_idx_old_cter:
			QtGui.QMessageBox.warning(None,"Warning",
							 "The format of CTER partres file (%s) might be old. We will stop supporting this format near future. \n\nPlease consider rerun CTER. Alternatively, please save selection now to generate *_partres_select.txt with the new format then open it. (Note: when you save selecction, *_partres_select.txt and *_partres_discard.txt will be saved with the new CTER partres format.)" 
							 % (file_path))
			# Continue processing for now (2017/03/20 Toshio Moriya)
		elif len(new_entry_list[0]) != self.n_idx_cter - self.n_idx_cter_extra:
			QtGui.QMessageBox.warning(None,"Warning","The number of columns (%d) has to be %d in %s." % (len(new_entry_list[0]), self.n_idx_cter - self.n_idx_cter_extra, file_path))
			return
		
		self.busy = True
		
		# print("MRK_DEBUG: ")
		# print("MRK_DEBUG: Detected %s entries in %s" % (len(new_entry_list), file_path))
		# print("MRK_DEBUG: Num. of Columns is %d in %s" % (len(new_entry_list[0]), file_path))
		# print("MRK_DEBUG: ")
		for cter_id in range(len(new_entry_list)):
			#
			# NOTE: 2017/11/22 Toshio Moriya
			# This CTEF file format is the latest
			#
			if len(new_entry_list[cter_id]) == self.n_idx_cter - self.n_idx_cter_extra:
				# Add extra items first to make sure indices match
				extended_entry = []
				extended_entry = extended_entry + [cter_id]                              # self.idx_cter_id , <extra> entry id
				extended_entry = extended_entry + [1]                                    # self.idx_cter_select  <extra> selected state
				
				extended_entry = extended_entry + new_entry_list[cter_id]                # original entry
				
###				if self.is_enable_max_power == True: 
###					extended_entry = extended_entry + [0.0] # MRK_TEST: self.idx_cter_max_power, <extra> maximum power in experimental rotational average (with astigmatism)
				
				# Store the extended entry to entry list
				new_entry_list[cter_id] = extended_entry
				
###				# MRK_TEST: Set max value of pwrot related to this micrograph
###				if self.is_enable_max_power == True: # MRK_TEST: 
###					new_cter_mic_file_path = new_entry_list[cter_id][self.idx_cter_mic_name]
###					# print "MRK_DEBUG: new_cter_mic_file_path := ", new_cter_mic_file_path
###					mic_basename_root = os.path.splitext(os.path.basename(new_cter_mic_file_path))
###					# print "MRK_DEBUG: mic_basename_root := ", mic_basename_root
###					new_cter_pwrot_file_path = os.path.join(os.path.dirname(file_path), "pwrot", "%s_rotinf.txt" % (mic_basename_root))
###					# print "MRK_DEBUG: new_cter_pwrot_file_path := ", new_cter_pwrot_file_path
###					new_rotinf_table = read_text_file(new_cter_pwrot_file_path, ncol=-1) # MRK_TEST: 
###					new_entry_list[cter_id][self.idx_cter_max_power] = max(new_rotinf_table[self.idx_rotinf_exp_with_astig]) # MRK_TEST: 
				
				# Always set selection state to 1 (selected)
				new_entry_list[cter_id][self.idx_cter_select] = 1
			#
			# NOTE: 2017/11/22 Toshio Moriya
			# The following code is to support the old format of CTER partres. It should be removed near future
			# 
			elif len(new_entry_list[cter_id]) == self.n_idx_old_cter:
				# assume amplitude amplitude contrast is total amplitude constrast in [%], estimated as variable Volta phase shift. Conver it to [deg].
				# Also, assuming constant amplitude contrast is zero since there is no information available in the old format.
				#from morphology import ampcont2angle
				total_phase_shift = ampcont2angle(new_entry_list[cter_id][self.idx_old_cter_ac])
				
				# Add extra items first to make sure indices match
				extended_entry = []
				extended_entry = extended_entry + [cter_id]                                # self.idx_cter_id , <extra> entry id
				extended_entry = extended_entry + [1]                                      # self.idx_cter_select  <extra> selected state
				
				extended_entry = extended_entry + [new_entry_list[cter_id][self.idx_old_cter_def]]          # self.idx_cter_def
				extended_entry = extended_entry + [new_entry_list[cter_id][self.idx_old_cter_cs]]           # self.idx_cter_cs
				extended_entry = extended_entry + [new_entry_list[cter_id][self.idx_old_cter_vol]]          # self.idx_cter_vol
				extended_entry = extended_entry + [new_entry_list[cter_id][self.idx_old_cter_apix]]         # self.idx_cter_apix
				extended_entry = extended_entry + [new_entry_list[cter_id][self.idx_old_cter_bfactor]]      # self.idx_cter_bfactor
				extended_entry = extended_entry + [new_entry_list[cter_id][self.idx_old_cter_ac]]           # self.idx_cter_total_ac
				extended_entry = extended_entry + [new_entry_list[cter_id][self.idx_old_cter_astig_amp]]    # self.idx_cter_astig_amp
				extended_entry = extended_entry + [new_entry_list[cter_id][self.idx_old_cter_astig_ang]]    # self.idx_cter_astig_ang
				extended_entry = extended_entry + [new_entry_list[cter_id][self.idx_old_cter_sd_def]]       # self.idx_cter_sd_def
				extended_entry = extended_entry + [new_entry_list[cter_id][self.idx_old_cter_sd_ac]]        # self.idx_cter_sd_total_ac
				extended_entry = extended_entry + [new_entry_list[cter_id][self.idx_old_cter_sd_astig_amp]] # self.idx_cter_sd_astig_amp
				extended_entry = extended_entry + [new_entry_list[cter_id][self.idx_old_cter_sd_astig_amp]] # self.idx_cter_sd_astig_ang
				extended_entry = extended_entry + [new_entry_list[cter_id][self.idx_old_cter_cv_def]]       # self.idx_cter_cv_def 
				extended_entry = extended_entry + [new_entry_list[cter_id][self.idx_old_cter_cv_astig_amp]] # self.idx_cter_cv_astig_amp 
				# Ignore self.idx_old_cter_spectra_diff
				extended_entry = extended_entry + [new_entry_list[cter_id][self.idx_old_cter_error_def]]    # self.idx_cter_error_def 
				extended_entry = extended_entry + [new_entry_list[cter_id][self.idx_old_cter_error_astig]]  # self.idx_cter_error_astig 
				extended_entry = extended_entry + [new_entry_list[cter_id][self.idx_old_cter_error_ctf]]    # self.idx_cter_error_ctf 
				extended_entry = extended_entry + [0.5/new_entry_list[cter_id][self.idx_old_cter_apix]]     # self.idx_cter_max_freq. Set to Nyquist frequency.
				extended_entry = extended_entry + [0.0]                                                     # self.idx_cter_reserved
				extended_entry = extended_entry + [0.0]                                                     # self.idx_cter_const_ac 
				extended_entry = extended_entry + [total_phase_shift]                                       # self.idx_cter_phase_shift 
				extended_entry = extended_entry + [new_entry_list[cter_id][self.idx_old_cter_mic_name ]]    # self.idx_cter_mic_name 
				
				# Store the extended entry to entry list
				new_entry_list[cter_id] = extended_entry
				
###				# MRK_TEST: Set max value of pwrot related to this micrograph
###				if self.is_enable_max_power == True: # MRK_TEST: 
###					new_cter_mic_file_path = new_entry_list[cter_id][self.idx_cter_mic_name]
###					# print "MRK_DEBUG: new_cter_mic_file_path := ", new_cter_mic_file_path
###					mic_basename_root = os.path.splitext(os.path.basename(new_cter_mic_file_path))
###					# print "MRK_DEBUG: mic_basename_root := ", mic_basename_root
###					new_cter_pwrot_file_path = os.path.join(os.path.dirname(file_path), "pwrot", "%s_rotinf.txt" % (mic_basename_root))
###					# print "MRK_DEBUG: new_cter_pwrot_file_path := ", new_cter_pwrot_file_path
###					new_rotinf_table = read_text_file(new_cter_pwrot_file_path, ncol=-1) # MRK_TEST: 
###					new_entry_list[cter_id][self.idx_cter_max_power] = max(new_rotinf_table[self.idx_rotinf_exp_with_astig]) # MRK_TEST: 
				
				# Always set selection state to 1 (selected)
				new_entry_list[cter_id][self.idx_cter_select] = 1
			#
			# NOTE: 2017/11/22 Toshio Moriya
			# Removed following code because it is too old at this point.
			#
			# NOTE: 2017/03/20 Toshio Moriya
			# The following code is to support the old format of CTER partres. It should be removed near future
			#
			# elif len(new_entry_list[cter_id]) == self.n_idx_cter - self.n_idx_cter_extra - 1:
			# 	# Add extra items first to make sure indices match
			# 	extended_entry = []
			# 	extended_entry = extended_entry + [cter_id]                              # self.idx_cter_id , <extra> entry id
			# 	extended_entry = extended_entry + [1]                                    # self.idx_cter_select  <extra> selected state
			# 	
			# 	# original entry
			# 	for original_entry_id in range(0, self.idx_cter_sd_total_ac): # From self.idx_cter_def to self.idx_cter_sd_def
			# 		extended_entry = extended_entry + [new_entry_list[cter_id][original_entry_id]]
			# 	extended_entry = extended_entry + [0.0] # self.idx_cter_sd_total_ac
			# 	for original_entry_id in range(self.idx_cter_sd_total_ac, len(new_entry_list[0])): # From self.idx_cter_sd_astig_amp to self.idx_cter_mic_name 
			# 		extended_entry = extended_entry + [new_entry_list[cter_id][original_entry_id]]
			# 	
#			# 	extended_entry = extended_entry + [""]                                   # self.idx_cter_pwrot_name, <extra> CTER power spectrum rotational average file name
#			# 	extended_entry = extended_entry + [0.5]                                  # self.idx_cter_error_ctf, <extra> limit frequency by CTF error 
			# 	if self.is_enable_max_power == True: extended_entry = extended_entry + [0.0] # MRK_TEST: self.idx_cter_max_power, <extra> maximum power in experimental rotational average (with astigmatism)
			# 	
			# 	# Store the extended entry to entry list
			# 	new_entry_list[cter_id] = extended_entry
			# 	
			# 	# MRK_TEST: Set max value of pwrot related to this micrograph
			# 	if self.is_enable_max_power == True: # MRK_TEST: 
			# 		new_cter_mic_file_path = new_entry_list[cter_id][self.idx_cter_mic_name]
#			# 		print "MRK_DEBUG: new_cter_mic_file_path := ", new_cter_mic_file_path
			# 		mic_basename_root = os.path.splitext(os.path.basename(new_cter_mic_file_path))
#			# 		print "MRK_DEBUG: mic_basename_root := ", mic_basename_root
			# 		new_cter_pwrot_file_path = os.path.join(os.path.dirname(file_path), "pwrot", "%s_rotinf.txt" % (mic_basename_root))
#			# 		print "MRK_DEBUG: new_cter_pwrot_file_path := ", new_cter_pwrot_file_path
			# 		new_rotinf_table = read_text_file(new_cter_pwrot_file_path, ncol=-1) # MRK_TEST: 
			# 		new_entry_list[cter_id][self.idx_cter_max_power] = max(new_rotinf_table[self.idx_rotinf_exp_with_astig]) # MRK_TEST: 
			# 	
			# 	# Always set selection state to 1 (selected)
			# 	new_entry_list[cter_id][self.idx_cter_select] = 1
			else: 
				assert False, "MRK_DEBUG: Unreachable code! Found Invalid number of columns (%d) in %s" % (len(new_entry_list[0]), file_path)
		
		# now set the new status
		
		self.cter_partres_file_path = file_path
		self.cter_entry_list = new_entry_list
		
		# Set the values and ranges of thresholds
		for idx_hist in range(self.n_idx_hist):
			idx_cter = self.hist_map_list[idx_hist][self.idx_hist_item_idx_cter]
			val_min = min(self.cter_entry_list, key=lambda x:x[idx_cter])[idx_cter]
			val_min = round(val_min, self.round_ndigits)
			val_max = max(self.cter_entry_list, key=lambda x:x[idx_cter])[idx_cter]
			val_max = round(val_max, self.round_ndigits)
			self.hist_map_list[idx_hist][self.idx_hist_item_val_min] = val_min
			self.hist_map_list[idx_hist][self.idx_hist_item_val_max] = val_max
			self.hist_map_list[idx_hist][self.idx_hist_item_unapply_threshold_lower] = val_min
			self.hist_map_list[idx_hist][self.idx_hist_item_unapply_threshold_upper] = val_max
			# self.hist_map_list[idx_hist][self.idx_hist_item_unapply_widget_lower].setRange(val_min, val_max)
			self.hist_map_list[idx_hist][self.idx_hist_item_unapply_widget_lower].setValue(val_min)
			# self.hist_map_list[idx_hist][self.idx_hist_item_unapply_widget_upper].setRange(val_min, val_max)
			self.hist_map_list[idx_hist][self.idx_hist_item_unapply_widget_upper].setValue(val_max)
			self.hist_map_list[idx_hist][self.idx_hist_item_apply_threshold_lower] = val_min
			self.hist_map_list[idx_hist][self.idx_hist_item_apply_threshold_upper] = val_max
			# self.hist_map_list[idx_hist][self.idx_hist_item_apply_widget_lower].setRange(val_min, val_max)
			self.hist_map_list[idx_hist][self.idx_hist_item_apply_widget_lower].setValue(val_min)
			# self.hist_map_list[idx_hist][self.idx_hist_item_apply_widget_upper].setRange(val_min, val_max)
			self.hist_map_list[idx_hist][self.idx_hist_item_apply_widget_upper].setValue(val_max)
			
			idx_sort = self.hist_map_list[idx_hist][self.idx_hist_item_idx_sort]
			if val_min == val_max:
				self.shist.model().item(idx_hist).setEnabled(False)
				self.ssort.model().item(idx_sort).setEnabled(False)
				self.hist_map_list[idx_hist][self.idx_hist_item_unapply_widget_lower].text.setStyleSheet("color: rgb(127,127,127);")
				self.hist_map_list[idx_hist][self.idx_hist_item_unapply_widget_upper].text.setStyleSheet("color: rgb(127,127,127);")
				self.hist_map_list[idx_hist][self.idx_hist_item_apply_widget_lower].text.setStyleSheet("color: rgb(127,127,127);")
				self.hist_map_list[idx_hist][self.idx_hist_item_apply_widget_upper].text.setStyleSheet("color: rgb(127,127,127);")
			else:
				assert val_min < val_max, "MRK_DEBUG"
				self.shist.model().item(idx_hist).setEnabled(True)
				self.ssort.model().item(idx_sort).setEnabled(True)
				self.hist_map_list[idx_hist][self.idx_hist_item_unapply_widget_lower].text.setStyleSheet("color: rgb(0,0,255);")
				self.hist_map_list[idx_hist][self.idx_hist_item_unapply_widget_upper].text.setStyleSheet("color: rgb(255,0,0);")
				self.hist_map_list[idx_hist][self.idx_hist_item_apply_widget_lower].text.setStyleSheet("color: rgb(0,0,255);")
				self.hist_map_list[idx_hist][self.idx_hist_item_apply_widget_upper].text.setStyleSheet("color: rgb(255,0,0);")
		
		# Set disable status of histogram
		if self.hist_map_list[self.curhistidx][self.idx_hist_item_val_min] == self.hist_map_list[self.curhistidx][self.idx_hist_item_val_max]:
			idx_cter = self.hist_map_list[self.curhistidx][self.idx_hist_item_idx_cter]
			#self.curhistdisable=True
			self.curhistogramdisplay = False
			if self.whistparam.isVisible():
				self.whistparam.hide()
			if self.wscatterparam.isVisible():
				self.wscatterparam.hide()
			# Error message of this condition should be displayed at the end of this function for smooth visual presentation
		
		self.updateEntryList()
		
		# Set the number of entries
		self.vbnentry.setValue(len(self.cter_entry_list))
		
		# Set the range of histogram bin
		# self.vsentryperbin.setRange(1,len(self.cter_entry_list))
		self.vsentryperbin.setValue(self.curentryperbin)
		
		self.updateUncheckCounts()
		
		# Enable buttons
		self.pbrefreshgraphs.setEnabled(True)
		self.pbreapplysort.setEnabled(True)
		self.pbapplyallthreshold.setEnabled(True)
		self.pbsavethresholdset.setEnabled(True)
		self.pbloadthresholdset.setEnabled(True)
		self.pbsaveselection.setEnabled(True)
		
#		# Disable rotational average plot display at the beginning of loading new dataset
#		if self.wplotrotavgcoarse.isVisible():
#			self.wplotrotavgcoarse.hide()
#		if self.wplotrotavgfine.isVisible():
#			self.wplotrotavgfine.hide()
#		self.cbrotavgdisplay.setValue(False)
#		self.curplotrotavgdisplay = False
		
#		# Disable micrograph display at the beginning of loading new dataset
#		if self.wimgmicthumb.isVisible():
#			self.wimgmicthumb.hide()
#		self.cbmicthumbdisplay.setValue(False)
#		self.curimgmicthumbdisplay = False
		
		cter_pwrot_dir = os.path.join(os.path.dirname(self.cter_partres_file_path), "pwrot")
		# print "MRK_DEBUG: cter_pwrot_dir = \"%s\" in readCterPartresFile() "% (cter_pwrot_dir)
		if os.path.exists(cter_pwrot_dir):
			# if not self.cbrotavgdisplay.getEnabled(): # MRK_NOTE: 2017/11/22 Toshio Moriya: This method does not work as I expected
			self.cbrotavgdisplay.setEnabled(True)
			for idx_graph in range(self.n_idx_graph):
				self.graph_map_list[idx_graph][self.idx_graph_item_widget].setEnabled(True)
			self.vbplotfixscale.setEnabled(True)
		else:
			# if self.cbrotavgdisplay.getEnabled(): # MRK_NOTE: 2017/11/22 Toshio Moriya: This method does not work as I expected
			self.cbrotavgdisplay.setEnabled(False)
			for idx_graph in range(self.n_idx_graph):
				self.graph_map_list[idx_graph][self.idx_graph_item_widget].setEnabled(False)
			self.vbplotfixscale.setEnabled(False)
			# Error message of this condition should be displayed at the end of this function for smooth visual presentation
			# QtGui.QMessageBox.warning(None,"Warning","Cannot find \"%s\" sub-directory associated with specified CTER partres file (%s). Please check your project directory. \n\nPower spectrum rotational average plots display option is disabled for this session." % (cter_pwrot_dir, self.cter_partres_file_path))
		
		cter_micthumb_dir = os.path.join(os.path.dirname(self.cter_partres_file_path), "micthumb")
		# print "MRK_DEBUG: cter_micthumb_dir = \"%s\" in readCterPartresFile() "% (cter_micthumb_dir)
		if os.path.exists(cter_micthumb_dir):
			# if not self.cbmicthumbdisplay.getEnabled(): # MRK_NOTE: 2016/03/22 Toshio Moriya: This method does not work as I expected
			self.cbmicthumbdisplay.setEnabled(True)
		else:
			# if self.cbmicthumbdisplay.getEnabled() != self.curimgmicthumbdisplay: # MRK_NOTE: 2016/03/22 Toshio Moriya: This method does not work as I expected
			self.cbmicthumbdisplay.setEnabled(False)
			# Error message of this condition should be displayed at the end of this function for smooth visual presentation
			# QtGui.QMessageBox.warning(None,"Warning","Cannot find \"%s\" sub-directory associated with specified CTER partres file (%s). Please check your project directory. \n\nMicrograph thumbnail display option is disabled for this session." % (cter_micthumb_dir, self.cter_partres_file_path))
			
		# NOTE: 2016/01/03 Toshio Moriya
		# Force update related plots to hide too much scaling delay...
		self.updateImgMicThumb(False)
		self.updateHist()
		self.updatePlotParam()
		self.updateFFT()
		####self.updatePlotCurves()  # REDUNDANT?
		
###		print("MRK_DEBUG: ")
###		print("MRK_DEBUG: readCterPartresFile(): self.curplotrotavgdisplay = ", self.curplotrotavgdisplay)
###		print("MRK_DEBUG: readCterPartresFile(): self.wplotrotavgcoarse.isVisible() = ", self.wplotrotavgcoarse.isVisible())
###		print("MRK_DEBUG: readCterPartresFile(): self.wplotrotavgcoarse.isVisible() = ", self.wplotrotavgcoarse.isVisible())
		
		# Coarse plot
		if self.curplotrotavgdisplay:
			if not self.wplotrotavgcoarse.isVisible():
				self.wplotrotavgcoarse.show()
		else:
			assert not self.curplotrotavgdisplay, "MRK_DEBUG"
			if self.wplotrotavgcoarse.isVisible():
				self.wplotrotavgcoarse.hide()
		
		# Zoomed plot
		if self.curplotrotzoomdisplay:
			if not self.wplotrotavgfine.isVisible():
				self.wplotrotavgfine.show()
		else:
			assert not self.curplotrotavgdisplay, "MRK_DEBUG"
			if self.wplotrotavgfine.isVisible():
				self.wplotrotavgfine.hide()
###		print("MRK_DEBUG: readCterPartresFile(): self.wplotrotavgcoarse.isVisible() = ", self.wplotrotavgcoarse.isVisible())
###		print("MRK_DEBUG: readCterPartresFile(): self.wplotrotavgcoarse.isVisible() = ", self.wplotrotavgcoarse.isVisible())
		
###		print("MRK_DEBUG: ")
###		print("MRK_DEBUG: readCterPartresFile(): self.curimgmicthumbdisplay = ", self.curimgmicthumbdisplay)
###		print("MRK_DEBUG: readCterPartresFile(): self.wimgmicthumb.isVisible() = ", self.wimgmicthumb.isVisible())
		if self.curimgmicthumbdisplay:
			if not self.wimgmicthumb.isVisible():
				self.wimgmicthumb.show()
		else:
			assert not self.curimgmicthumbdisplay, "MRK_DEBUG"
			if self.wimgmicthumb.isVisible():
				self.wimgmicthumb.hide()
###		print("MRK_DEBUG: readCterPartresFile(): self.wimgmicthumb.isVisible() = ", self.wimgmicthumb.isVisible())
		
		if self.curfftdisplay:
			if not self.wfft.isVisible():
				self.wfft.show()
		else:
			assert not self.curfftdisplay, "MRK_DEBUG"
			if self.wfft.isVisible():
				self.wfft.hide()
		
		self.busy = False
		self.needredisp = True
		
		if self.hist_map_list[self.curhistidx][self.idx_hist_item_val_min] == self.hist_map_list[self.curhistidx][self.idx_hist_item_val_max]:
			param_label = self.value_map_list[idx_cter][self.idx_cter_item_label]
			QtGui.QMessageBox.information(self, "Information",
								 "All entries have the same selected parameter values (%s). \n\nParameter Histogram & Plot will not be shown." 
								 % (param_label))
		
		if not os.path.exists(cter_pwrot_dir):
			QtGui.QMessageBox.warning(None,"Warning",
							 "Cannot find \"%s\" sub-directory associated with specified CTER partres file (%s). Please check your project directory. \n\nPower spectrum rotational average plots display option is disabled for this session." 
							 % (cter_pwrot_dir, self.cter_partres_file_path))

		if not os.path.exists(cter_micthumb_dir):
			QtGui.QMessageBox.warning(None,"Warning",
							 "Cannot find \"%s\" sub-directory associated with specified CTER partres file (%s). Please check your project directory. \n\nMicrograph thumbnail display option is disabled for this session." 
							 % (cter_micthumb_dir, self.cter_partres_file_path))

#		assert self.isVisible(), "MRK_DEBUG"
#		self.raise_()
#		self.activateWindow()
	
	def openCterPartres(self,val=None):
		"""Open CTER partres file"""
		
		file_path = str(QtGui.QFileDialog.getOpenFileName(self, "Open CTER partres File", options = QtGui.QFileDialog.DontUseNativeDialog))
		if file_path == "": return
		
		self.readCterPartresFile(os.path.relpath(file_path))
	
	def updateHist(self, error_display = True):
		if self.whistparam == None: return # it's closed/not visible
		if self.cter_partres_file_path == None: return # no cter ctf file is selected
		if self.cter_entry_list == None: return # no cter ctf file is selected
		#if self.curhistdisable == True: return # do nothing while it is hidden
		if not self.curhistogramdisplay:
			#print('self.curhistogramdisplay',self.curhistogramdisplay)
			return
	
		val_list = []
		
		# Create Histogram for selected parameter
		idx_cter = self.hist_map_list[self.curhistidx][self.idx_hist_item_idx_cter]
		for cter_entry in self.cter_entry_list:
			val_list.append(cter_entry[idx_cter])
		
		n_bin = 1
		if len(self.cter_entry_list) < self.curentryperbin:
			self.curentryperbin = len(self.cter_entry_list)
			self.vsentryperbin.setValue(self.curentryperbin)
			# self.vsentryperbin.setRange(1,len(self.cter_entry_list))
		elif self.curentryperbin < 1:
			self.curentryperbin = 1
			self.vsentryperbin.setValue(self.curentryperbin)
			# self.vsentryperbin.setRange(1,len(self.cter_entry_list))
		n_bin = len(self.cter_entry_list)/ self.curentryperbin
		assert len(val_list) >= n_bin, "MRK_DEBUG"
		assert n_bin > 0, "MRK_DEBUG"
		#from statistics import hist_list
		hist_x_list, hist_y_list = hist_list(val_list, n_bin)
		
		# Pad with zero for better visual impression...
		hist_x_list += [max(val_list)]
		hist_y_list += [0]
		self.whistparam.set_data((hist_x_list,hist_y_list),"hist_param",quiet=False,color=0,linetype=0,symtype=0)
		
		# MRK_NOTE: 2015/12/17 Toshio Moriya
		# This may NOT be good place to update the following information...
		idx_cter = self.hist_map_list[self.curhistidx][self.idx_hist_item_idx_cter]
		param_label = self.value_map_list[idx_cter][self.idx_cter_item_label]
		
		# Get value for current micrograph
		param_val = self.cter_entry_list[self.curentry][idx_cter]
		
		# shape_name = "hist_param_shape_value"
		# scr_x, scr_y = self.whistparam.plot2scr(param_val, 0.0)
		# self.whistparam.add_shape(shape_name,EMShape(("scrline",0,1,0,scr_x,self.whistparam.scrlim[1],scr_x,self.whistparam.scrlim[1]+self.whistparam.scrlim[3],3)))
		
		# Get ordinate range in order to overlay current micrograph's value on histogram
		val_min = round(min(hist_y_list), self.round_ndigits)
		val_max = round(max(hist_y_list), self.round_ndigits)
		
		# threshold_lower_label = "Lower(Blue)"
		unapply_threshold_lower_val = round(self.hist_map_list[self.curhistidx][self.idx_hist_item_unapply_threshold_lower], self.round_ndigits)
		apply_threshold_lower_val = round(self.hist_map_list[self.curhistidx][self.idx_hist_item_apply_threshold_lower], self.round_ndigits)
		# threshold_upper_label = "Upper(Red)"
		unapply_threshold_upper_val = round(self.hist_map_list[self.curhistidx][self.idx_hist_item_unapply_threshold_upper], self.round_ndigits)
		apply_threshold_upper_val = round(self.hist_map_list[self.curhistidx][self.idx_hist_item_apply_threshold_upper], self.round_ndigits)
		
		#print('updateHist', param_val, val_min, val_max)
		#print('updateHist', [param_val, param_val], [val_min, val_max])
		self.whistparam.set_data(([param_val, param_val], [val_min, val_max]),"selected_val",quiet=False,color=3, linetype=1)  # default linetype was hard to see
		self.whistparam.set_data(([unapply_threshold_lower_val, unapply_threshold_lower_val], [val_min, val_max]),"unapply_threshold_lower_val",quiet=False,color=1,linetype=1)
		self.whistparam.set_data(([apply_threshold_lower_val, apply_threshold_lower_val], [val_min, val_max]),"apply_threshold_lower_val",quiet=False,color=1)
		self.whistparam.set_data(([unapply_threshold_upper_val, unapply_threshold_upper_val], [val_min, val_max]),"unapply_threshold_upper_val",quiet=False,color=2,linetype=1)
		self.whistparam.set_data(([apply_threshold_upper_val, apply_threshold_upper_val], [val_min, val_max]),"apply_threshold_upper_val",quiet=False,color=2)
		
		# shape_name = "hist_param_shape_unapply_threshold_lower"
		# scr_x, scr_y = self.whistparam.plot2scr(unapply_threshold_lower_val, 0.0)
		#self.whistparam.add_shape(shape_name,EMShape(("scrline",0,0,0.5,scr_x,self.whistparam.scrlim[1],scr_x,self.whistparam.scrlim[1]+self.whistparam.scrlim[3],1)))
		
		# shape_name = "hist_param_shape_apply_threshold_lower"
		# scr_x, scr_y = self.whistparam.plot2scr(apply_threshold_lower_val, 0.0)
		# self.whistparam.add_shape(shape_name,EMShape(("scrline",0,0,1,scr_x,self.whistparam.scrlim[1],scr_x,self.whistparam.scrlim[1]+self.whistparam.scrlim[3],1)))
		
		# shape_name = "hist_param_shape_unapply_threshold_upper"
		# scr_x, scr_y = self.whistparam.plot2scr(unapply_threshold_upper_val, 0.0)
		# self.whistparam.add_shape(shape_name,EMShape(("scrline",0.5,0,0,scr_x,self.whistparam.scrlim[1],scr_x,self.whistparam.scrlim[1]+self.whistparam.scrlim[3],1)))
		
		# shape_name = "hist_param_shape_apply_threshold_upper"
		# scr_x, scr_y = self.whistparam.plot2scr(apply_threshold_upper_val, 0.0)
		# self.whistparam.add_shape(shape_name,EMShape(("scrline",1,0,0,scr_x,self.whistparam.scrlim[1],scr_x,self.whistparam.scrlim[1]+self.whistparam.scrlim[3],1)))
		
		# shape_name = "hist_param_shape_label"
		# if self.curhistidx in [self.idx_hist_error_def, self.idx_hist_error_astig, self.idx_hist_error_ctf] and param_val > 0.0:
			#self.whistparam.add_shape(shape_name,EMShape(("scrlabel",0,0,0,self.whistparam.scrlim[0]+30,self.whistparam.scrlim[1]+self.whistparam.scrlim[3]-18,"%s(Green) %1.3g (%1.3g), %s %1.3g, %s %1.3g"%(param_label,param_val,
																												#1/param_val,threshold_lower_label,apply_threshold_lower_val,threshold_upper_label,apply_threshold_upper_val),120.0,-1)))
		# 	# self.whistparam.add_shape(shape_name,EMShape(("scrlabel",0,0,0,self.whistparam.scrlim[0]+30,self.whistparam.scrlim[1]+self.whistparam.scrlim[3]-18,"%s %1.5g (%1.5g)"%(param_label,param_val,1/param_val),120.0,-1)))
		# else:
			#self.whistparam.add_shape(shape_name,EMShape(("scrlabel",0,0,0,self.whistparam.scrlim[0]+30,self.whistparam.scrlim[1]+self.whistparam.scrlim[3]-18,"%s(Green) %1.3g, %s %1.3g, %s %1.3g"%(param_label,param_val,
																															#threshold_lower_label,apply_threshold_lower_val,threshold_upper_label,apply_threshold_upper_val),120.0,-1)))
		# 	# self.whistparam.add_shape(shape_name,EMShape(("scrlabel",0,0,0,self.whistparam.scrlim[0]+30,self.whistparam.scrlim[1]+self.whistparam.scrlim[3]-18,"%s %1.5g"%(param_label,param_val),120.0,-1)))
		
		self.whistparam.setAxisParms(param_label,"Image Counts")
		# x_margin = (hist_x_list[-1] - hist_x_list[0]) * 0.05
		# NOTE: 2016/01/02 Toshio Moriya
		# Disable manual rescale for now and use autoscale
		# self.whistparam.rescale(min(val_list),max(val_list),0,max(hist_y_list) * 1.05)
		self.whistparam.autoscale(True)
	
	def updatePlotParam(self, error_display = True):
		if self.wscatterparam == None: return # it's closed/not visible
		if self.cter_partres_file_path == None: return # no cter ctf file is selected
		if self.cter_entry_list == None: return # no cter ctf file is selected
		if not self.curscatterdisplay and not self.curhistogramdisplay: return
		
		x_list = []
		y_list = []
		
		# Create graph for selected parameter
		idx_cter = self.hist_map_list[self.curhistidx][self.idx_hist_item_idx_cter]
		for cter_id in range(len(self.cter_entry_list)):
			x_list.append(cter_id)
			y_list.append(self.cter_entry_list[cter_id][idx_cter])
		# self.wscatterparam.set_data((x_list,y_list),"plot_param",quiet=False,color=0)
		self.wscatterparam.set_data((x_list,y_list),"plot_param",quiet=False,color=0,linetype=0,symtype=0)
		
		# Create graph for single parameter value of selected entry
		# MRK_NOTE: 2015/12/17 Toshio Moriya
		# This may NOT be good place to update the following information...
		idx_cter = self.hist_map_list[self.curhistidx][self.idx_hist_item_idx_cter]
		param_label = self.value_map_list[idx_cter][self.idx_cter_item_label]
		param_val = round(self.cter_entry_list[self.curentry][idx_cter], self.round_ndigits)
		
		# threshold_lower_label = "Lower(Blue)"
		unapply_threshold_lower_val = round(self.hist_map_list[self.curhistidx][self.idx_hist_item_unapply_threshold_lower], self.round_ndigits)
		apply_threshold_lower_val = round(self.hist_map_list[self.curhistidx][self.idx_hist_item_apply_threshold_lower], self.round_ndigits)
		# threshold_upper_label = "Upper(Red)"
		unapply_threshold_upper_val = round(self.hist_map_list[self.curhistidx][self.idx_hist_item_unapply_threshold_upper], self.round_ndigits)
		apply_threshold_upper_val = round(self.hist_map_list[self.curhistidx][self.idx_hist_item_apply_threshold_upper], self.round_ndigits)
		
		y_list = [param_val]*len(x_list)
		self.wscatterparam.set_data((x_list,y_list),"selected_val",quiet=False,color=3, linetype=1)  # default linetype was hard to see
		y_list = [unapply_threshold_lower_val]*len(x_list)
		self.wscatterparam.set_data((x_list,y_list),"unapply_threshold_lower_val",quiet=False,color=1,linetype=1)
		y_list = [apply_threshold_lower_val]*len(x_list)
		self.wscatterparam.set_data((x_list,y_list),"apply_threshold_lower_val",quiet=False,color=1)
		y_list = [unapply_threshold_upper_val]*len(x_list)
		self.wscatterparam.set_data((x_list,y_list),"unapply_threshold_upper_val",quiet=False,color=2,linetype=1)
		y_list = [apply_threshold_upper_val]*len(x_list)
		self.wscatterparam.set_data((x_list,y_list),"apply_threshold_upper_val",quiet=False,color=2)
		
		# shape_name = "plot_param_shape_label"
		# if self.curhistidx in [self.idx_hist_error_def, self.idx_hist_error_astig, self.idx_hist_error_ctf] and param_val > 0.0:
		# 	self.wscatterparam.add_shape(shape_name,EMShape(("scrlabel",0,0,0,self.wscatterparam.scrlim[0]+30,self.wscatterparam.scrlim[1]+self.wscatterparam.scrlim[3]-18,"%s(Green) %1.3g (%1.3g), %s %1.3g, %s %1.3g"%(param_label,param_val,1/param_val,threshold_lower_label,apply_threshold_lower_val,threshold_upper_label,apply_threshold_upper_val),120.0,-1)))
		# 	# self.wscatterparam.add_shape(shape_name,EMShape(("scrlabel",0,0,0,self.wscatterparam.scrlim[0]+30,self.wscatterparam.scrlim[1]+self.wscatterparam.scrlim[3]-18,"%s(Green) %1.5g (%1.5g)"%(param_label,param_val,1/param_val),120.0,-1)))
		# else:
		# 	self.wscatterparam.add_shape(shape_name,EMShape(("scrlabel",0,0,0,self.wscatterparam.scrlim[0]+30,self.wscatterparam.scrlim[1]+self.wscatterparam.scrlim[3]-18,"%s(Green) %1.3g, %s %1.3g, %s %1.3g"%(param_label,param_val,threshold_lower_label,apply_threshold_lower_val,threshold_upper_label,apply_threshold_upper_val),120.0,-1)))
		# 	# self.wscatterparam.add_shape(shape_name,EMShape(("scrlabel",0,0,0,self.wscatterparam.scrlim[0]+30,self.wscatterparam.scrlim[1]+self.wscatterparam.scrlim[3]-18,"%s(Green) %1.5g"%(param_label,param_val),120.0,-1)))
		
		self.wscatterparam.setAxisParms("Sorted Image ID", param_label)
		# NOTE: 2016/01/02 Toshio Moriya
		# Use autoscale for now
		self.wscatterparam.autoscale(True)
	
	def updatePlotCurves(self, error_display = True):
		if not self.curplotrotavgdisplay and not self.curplotrotzoomdisplay:
			return # Both power-spectrum rotational-average plot displays are unchecked
		if self.wplotrotavgcoarse == None and self.wplotrotavgfine == None: return # closed/not visible
		if self.cter_pwrot_file_path == None: return # no cter entry is selected
		
		# print "MRK_DEBUG: self.cter_pwrot_file_path =\"%s\" in updatePlotCurves() "% (self.cter_pwrot_file_path)
		if not os.path.exists(self.cter_pwrot_file_path):
			if self.wplotrotavgcoarse.isVisible():
				self.wplotrotavgcoarse.hide()
			if self.wplotrotavgfine.isVisible():
				self.wplotrotavgfine.hide()
			if self.wfft.isVisible():
				self.wfft.hide()
			if error_display and os.path.exists(os.path.dirname(self.cter_pwrot_file_path)):
				QtGui.QMessageBox.warning(None,"Warning","Cannot find file cter_pwrot_file_path (%s). Please check the contents of pwrot directory. \n\nPlots will not be shown." 
							  % (self.cter_pwrot_file_path))
			return
		assert os.path.exists(self.cter_pwrot_file_path), "MRK_DEBUG"
		
		# Now update the plots
		self.rotinf_table = read_text_file(self.cter_pwrot_file_path, ncol=-1)
		
		# print "MRK_DEBUG: Last entry of the 1st colum should be a micrograph name %s which is same as " % os.path.basename(self.rotinf_table[0][-1])
		
#		mic_basename_rotinf = os.path.basename(self.rotinf_table[0][-1]) # last entry of 1st colum should be associated micrograph
#		mic_basename_partres = os.path.basename(self.cter_entry_list[self.curentry][self.idx_cter_mic_name])
#		
#		if mic_basename_rotinf != mic_basename_partres:
#			QtGui.QMessageBox.warning(None,"Warning","Micrograph name (%s) in %s is not same as name (%s) in %s." % (mic_basename_rotinf, os.path.basename(self.cter_pwrot_file_path), mic_basename_partres, os.path.basename(self.cter_partres_file_path)))
#			return
		
		# global_min = float("inf")
		# global_max = float("-inf")
		
		# Spatial frequency
		spFreqList = self.rotinf_table[self.idx_rotinf_freq]
		
		# Loop through power-spectrum profiles
		for idx_graph in range(self.n_idx_graph):
			# List of plots, containing: label, title, index, checkboxValue
			plotList = self.graph_map_list[idx_graph]
			
			# Partres file has 4 columns
			columnValues = self.rotinf_table[plotList[self.idx_graph_idx_rotinf]]
		
			self.wplotrotavgcoarse.set_data((spFreqList,columnValues), plotList[self.idx_graph_item_name], quiet=False, color=idx_graph, linetype=0)
			self.wplotrotavgfine.set_data((spFreqList, columnValues),  plotList[self.idx_graph_item_name], quiet=False, color=idx_graph, linetype=0)
			# val_min = min(columnValues)
			# val_max = max(columnValues)
			# if global_min > val_min:
			# 	global_min = val_min
			# if global_max < val_max:
			# 	global_max = val_max
		
		
		# NOTE: 2016/01/02 Toshio Moriya
		# Disable manual rescale for now and use autoscale
		# self.wplotrotavgcoarse.rescale(spFreqList[0],spFreqList[-1],0.0,1.0)
		#self.wplotrotavgcoarse.autoscale(True)
		self.wplotrotavgfine.rescale(spFreqList[0],spFreqList[-1],0.0,self.curplotfixscale)
		
		nyquist_freq = spFreqList[-1]
		# print "MRK_DEBUG: nyquist_freq = %1.5g" % nyquist_freq
		
		fineLimits   = self.wplotrotavgfine.scrlim
		#print ('self.wplotrotavgfine.scrlim',type(self.wplotrotavgfine.scrlim), self.wplotrotavgfine.scrlim)
		coarseLimits = self.wplotrotavgcoarse.scrlim
		
		error_name = "error_astig"
		error_label = "Astig. Limit"
		error_freq = self.cter_entry_list[self.curentry][self.idx_cter_error_astig]
		# print "MRK_DEBUG: %s= %1.5g" % (error_name, error_freq)
		if error_freq > 0.0 and error_freq <= nyquist_freq:
			error_scr_x, error_scr_y = self.wplotrotavgcoarse.plot2scr(error_freq, 0.0)
			self.wplotrotavgcoarse.add_shape(error_name,EMShape(("scrline",0,0,0.5,error_scr_x,coarseLimits[1],error_scr_x,coarseLimits[1]+coarseLimits[3],1)))
			# self.wplotrotavgcoarse.set_data(([error_freq, error_freq], [global_min, global_max]),"astig_error_freq_limit",quiet=False,color=0,linetype=0)
			self.wplotrotavgcoarse.add_shape("%s_freq"%(error_name),
								  EMShape(("scrlabel",0,0,0,error_scr_x-260,coarseLimits[1]+coarseLimits[3]-18,"%s %1.5g (%1.5g)"%(error_label,error_freq,1.0/error_freq),120.0,-1)))
			error_scr_x, error_scr_y = self.wplotrotavgfine.plot2scr(error_freq, 0.0)
			self.wplotrotavgfine.add_shape(error_name,EMShape(("scrline",0,0,0.5,error_scr_x,fineLimits[1],error_scr_x,fineLimits[1]+fineLimits[3],1)))
			# self.wplotrotavgfine.set_data(([error_freq, error_freq], [global_min, global_max]),"astig_error_freq_limit",quiet=False,color=0,linetype=0)
			self.wplotrotavgfine.add_shape("%s_freq"%(error_name),
								  EMShape(("scrlabel",0,0,0,error_scr_x-260,fineLimits[1]+fineLimits[3]-18,"%s %1.5g (%1.5g)"%(error_label,error_freq,1.0/error_freq),120.0,-1)))
		
		error_name = "error_def"
		error_label = "Defocus Limit"
		error_freq = self.cter_entry_list[self.curentry][self.idx_cter_error_def]
		# print "MRK_DEBUG: %s= %1.5g" % (error_name, error_freq)
		if error_freq > 0.0 and error_freq <= nyquist_freq:
			error_scr_x, error_scr_y = self.wplotrotavgcoarse.plot2scr(error_freq, 0.0)
			self.wplotrotavgcoarse.add_shape(error_name,EMShape(("scrline",0.5,0,0,error_scr_x,coarseLimits[1],error_scr_x,coarseLimits[1]+coarseLimits[3],1)))
			# self.wplotrotavgcoarse.set_data(([error_freq, error_freq], [global_min, global_max]),"defocus_error_freq_limit",quiet=False,color=0,linetype=0)
			self.wplotrotavgcoarse.add_shape("%s_freq"%(error_name),
									EMShape(("scrlabel",0,0,0,error_scr_x-260,coarseLimits[1]+coarseLimits[3]-36,"%s %1.5g (%1.5g)"%(error_label,error_freq,1.0/error_freq),120.0,-1)))
			error_scr_x, error_scr_y = self.wplotrotavgfine.plot2scr(error_freq, 0.0)
			self.wplotrotavgfine.add_shape(error_name,EMShape(("scrline",0.5,0,0,error_scr_x,fineLimits[1],error_scr_x,fineLimits[1]+fineLimits[3],1)))
			# self.wplotrotavgfine.set_data(([error_freq, error_freq], [global_min, global_max]),"defocus_error_freq_limit",quiet=False,color=0,linetype=0)
			self.wplotrotavgfine.add_shape("%s_freq"%(error_name),
								  EMShape(("scrlabel",0,0,0,error_scr_x-260,fineLimits[1]+fineLimits[3]-36,"%s %1.5g (%1.5g)"%(error_label,error_freq,1.0/error_freq),120.0,-1)))
		
		error_name = "error_ctf"
		error_label = "CTF Limit"
		error_freq = self.cter_entry_list[self.curentry][self.idx_cter_error_ctf]
		# print "MRK_DEBUG: %s= %1.5g" % (error_name, error_freq)
		if error_freq > 0.0 and error_freq <= nyquist_freq:
			error_scr_x, error_scr_y = self.wplotrotavgcoarse.plot2scr(error_freq, 0.0)
			self.wplotrotavgcoarse.add_shape(error_name,EMShape(("scrline",0,0.5,0,error_scr_x,coarseLimits[1],error_scr_x,coarseLimits[1]+coarseLimits[3],1)))
			# self.wplotrotavgcoarse.set_data(([error_freq, error_freq], [global_min, global_max]),"ctf_freq_limit")
			self.wplotrotavgcoarse.add_shape("%s_freq"%(error_name),
									EMShape(("scrlabel",0,0,0,error_scr_x-260,coarseLimits[1]+coarseLimits[3]-54,"%s %1.5g (%1.5g)"%(error_label,error_freq,1.0/error_freq),120.0,-1)))
			error_scr_x, error_scr_y = self.wplotrotavgfine.plot2scr(error_freq, 0.0)
			self.wplotrotavgfine.add_shape(error_name,
								  EMShape(("scrline",0,0.5,0,error_scr_x,fineLimits[1],error_scr_x,fineLimits[1]+fineLimits[3],1)))
			# self.wplotrotavgfine.set_data(([error_freq, error_freq], [global_min, global_max]),"ctf_freq_limit")
			self.wplotrotavgfine.add_shape("%s_freq"%(error_name),
								  EMShape(("scrlabel",0,0,0,error_scr_x-260,fineLimits[1]+fineLimits[3]-54,"%s %1.5g (%1.5g)"%(error_label,error_freq,1.0/error_freq),120.0,-1)))
		
		self.wplotrotavgcoarse.setAxisParms("frequency (1/"+ "$\AA$" +")","power spectrum")
		self.wplotrotavgfine.setAxisParms("frequency (1/"+ "$\AA$" +")","power spectrum")
###		self.wplotrotavgcoarse.setAxisParms("frequency (1/"+ "$\AA$" +")","amplitude")
###		self.wplotrotavgfine.setAxisParms("frequency (1/"+ "$\AA$" +")","amplitude")
		
		# Update plot
		self.updatePlotVisibility()
	
		
	def newRotAvgDisplay(self,val=None):
		"""Change rotational average plot display status."""
		# assert self.cbrotavgdisplay.getEnabled() == True, "MRK_DEBUG: 2017/11/22 Toshio Moriya: This method does not work as I expected"
		if self.curplotrotavgdisplay != val:
			# now set the new display status
			self.curplotrotavgdisplay = val
		else:
			assert self.curplotrotavgdisplay == val, "MRK_DEBUG"
			# The status did not change, there is nothing to do nothing
			return
		
		if self.cter_pwrot_file_path == None: return # no cter entry is selected
		
		if not os.path.exists(self.cter_pwrot_file_path):
			assert not self.wplotrotavgcoarse.isVisible(), "MRK_DEBUG"
			#assert not self.wplotrotavgfine.isVisible(), "MRK_DEBUG"
			return
		assert os.path.exists(self.cter_pwrot_file_path), "MRK_DEBUG"
		
		if self.curplotrotavgdisplay and not self.wplotrotavgcoarse.isVisible():
			self.wplotrotavgcoarse.show()
			self.needredisp = True
		elif not self.curplotrotavgdisplay and self.wplotrotavgcoarse.isVisible():
			self.wplotrotavgcoarse.hide()

	def newRotZoomDisplay(self,val=None):
		"""Change rotational average plot display status."""
		
		if self.curplotrotzoomdisplay != val:
			# now set the new display status
			self.curplotrotzoomdisplay = val
		else:
			assert self.curplotrotzoomdisplay == val, "MRK_DEBUG"
			# The status did not change, there is nothing to do
			return
		
		if self.cter_pwrot_file_path == None: return # no cter entry is selected
		
		if not os.path.exists(self.cter_pwrot_file_path):
			assert not self.wplotrotavgfine.isVisible(), "MRK_DEBUG"
			return
		assert os.path.exists(self.cter_pwrot_file_path), "MRK_DEBUG"
		
		if self.curplotrotzoomdisplay and not self.wplotrotavgfine.isVisible():
			self.wplotrotavgfine.show()
			self.needredisp = True
		elif not self.curplotrotzoomdisplay and self.wplotrotavgfine.isVisible():
			self.wplotrotavgfine.hide()
		
	def newHistogramDisplay(self,val=None):
		"""Change histogram display status."""
		
		if self.curhistogramdisplay != val:
			# now set the new display status
			self.curhistogramdisplay = val
		else:
			assert self.curhistogramdisplay == val, "MRK_DEBUG"
			# The status did not change, there is nothing to do
			return
		
		if self.cter_pwrot_file_path == None: return # no cter entry is selected
		
		if not os.path.exists(self.cter_pwrot_file_path):
			assert not self.wplotrotavgfine.isVisible(), "MRK_DEBUG"
			return
		assert os.path.exists(self.cter_pwrot_file_path), "MRK_DEBUG"
		
		if self.curhistogramdisplay and not self.whistparam.isVisible():
			self.whistparam.show()
			self.needredisp = True
		elif not self.curhistogramdisplay and self.whistparam.isVisible():
			self.whistparam.hide()
		
	def newScatterDisplay(self,val=None):
		"""Change sort plot display status."""
		
		#print('newScatterDisplay', val)
		if self.curscatterdisplay != val:
			# now set the new display status
			self.curscatterdisplay = val
		else:
			assert self.curscatterdisplay == val, "MRK_DEBUG"
			# The status did not change, there is nothing to do
			return
		
		if self.cter_pwrot_file_path == None: return # no cter entry is selected
		
		if not os.path.exists(self.cter_pwrot_file_path):
			assert not self.wplotrotavgfine.isVisible(), "MRK_DEBUG"
			return
		assert os.path.exists(self.cter_pwrot_file_path), "MRK_DEBUG"
		
		if self.curscatterdisplay and not self.wscatterparam.isVisible():
			self.wscatterparam.show()
			self.needredisp = True
		elif not self.curscatterdisplay and self.wscatterparam.isVisible():
			self.wscatterparam.hide()
		
	def updateEntryList(self):
		"""Updated entry list box after sorting of CTER partres entries based on current setting."""
		
		# sort CTER partres entry list
		assert self.cter_entry_list != None, "MRK_DEBUG"
		if self.cursortidx != self.idx_sort_mic_names:
			self.cter_entry_list = sorted(self.cter_entry_list, key=lambda x: x[self.sort_map_list[self.cursortidx][self.idx_sort_item_idx_cter]], reverse=self.cursortorder)
		else:
			assert (self.cursortidx == self.idx_sort_mic_names)
			self.cter_entry_list = sorted(self.cter_entry_list, key=lambda x: os.path.basename(x[self.sort_map_list[self.cursortidx][self.idx_sort_item_idx_cter]]), reverse=self.cursortorder)
		
		if self.cursortselect: 
			# Additionaly, sort cter entry list by select state
			self.cter_entry_list = sorted(self.cter_entry_list, key=lambda x: x[self.idx_cter_select])
		# else: # Do nothing
		
		# Refresh entry list box
		self.lbentry.clear()
		newItemflags = Qt.ItemFlags(Qt.ItemIsSelectable)|Qt.ItemFlags(Qt.ItemIsEnabled)|Qt.ItemFlags(Qt.ItemIsUserCheckable)
		for cter_entry in self.cter_entry_list:
			newItem = QtGui.QListWidgetItem(os.path.basename(cter_entry[self.idx_cter_mic_name]))
			newItem.setFlags(newItemflags)
			if cter_entry[self.idx_cter_select] == 1:
				newItem.setCheckState(Qt.Checked)
			else: 
				assert cter_entry[self.idx_cter_select] == 0, "MRK_DEBUG"
				newItem.setCheckState(Qt.Unchecked)
			self.lbentry.addItem(newItem)
			# self.lbentry.addItem(os.path.basename(cter_entry[self.idx_cter_mic_name]))
			
		self.newEntry(0)
		self.lbentry.setCurrentRow(0)

	def updateImgMicThumb(self, error_display = True):
		if not self.curimgmicthumbdisplay: return # Micrograph thumbnail display is unchecked
		if self.wimgmicthumb == None: return # it's closed/not visible
		if self.cter_micthumb_file_path == None: return # no cter entry is selected
		
		# print "MRK_DEBUG: self.cter_micthumb_file_path =\"%s\" in updateImgMicThumb() "% (self.cter_micthumb_file_path)
		if not os.path.exists(self.cter_micthumb_file_path):
			if self.wimgmicthumb.isVisible():
				self.wimgmicthumb.hide()
			if error_display and os.path.exists(os.path.dirname(self.cter_micthumb_file_path)):
				QtGui.QMessageBox.warning(None,"Warning","Cannot find micrograph thumbnail (%s). Please check your micrograph thumbnail directory. \n\nMicrograph thumbnail will not be shown." % (self.cter_micthumb_file_path))
			return
		assert os.path.exists(self.cter_micthumb_file_path), "MRK_DEBUG"
		
		# Now update the image
		micthumb_img = EMData(self.cter_micthumb_file_path) # read the image from disk
		self.wimgmicthumb.set_data(micthumb_img)
		self.wimgmicthumb.setWindowTitle("sxgui_cter - Micrograph Thumbnail- %s, %s" 
				% (os.path.basename(self.cter_entry_list[self.curentry][self.idx_cter_mic_name]), 
				os.path.basename(self.cter_micthumb_file_path)))
		
#		# print "MRK_DEBUG: self.cter_mic_file_path =\"%s\" in updateImgMicThumb() "% (self.cter_mic_file_path)
#		if not os.path.exists(self.cter_mic_file_path):
#			QtGui.QMessageBox.warning(None,"Warning","Cannot find micrograph (%s). Please check your micrograph directory. \n\nA blank image is shown." % (self.cter_mic_file_path))
#			mic_img = EMData() # Set empty image...
#			img_size = 4096 # Use the most typical image size?!!!
#			mic_img = model_blank(img_size,img_size, bckg=1.0)
#		else:
#			mic_img = EMData(self.cter_mic_file_path) # read the image from disk
#		self.wimage.set_data(mic_img)
#		self.wimage.setWindowTitle("sxgui_cter - Micrograph - %s, %s" % (os.path.basename(self.cter_entry_list[self.curentry][self.idx_cter_mic_name]), os.path.basename(self.cter_mic_file_path)))
		
	def newMicThumbDisplay(self,val=None):
		"""Change micrograph thumbnail display status."""
		# assert self.cbmicthumbdisplay.getEnabled() == True, "MRK_DEBUG: 2016/03/22 Toshio Moriya: This method does not work as I expected"
		if self.curimgmicthumbdisplay != val:
			# now set the new display status
			self.curimgmicthumbdisplay = val
		else: 
			assert self.curimgmicthumbdisplay == val, "MRK_DEBUG"
			# The status did not change, there is nothing to do nothing
			return
		
		if self.cter_micthumb_file_path == None: return # no cter entry is selected
		
		# print "MRK_DEBUG: self.cter_micthumb_file_path =\"%s\" in newMicThumbDisplay() "% (self.cter_micthumb_file_path)
		if not os.path.exists(self.cter_micthumb_file_path):
			assert not self.wimgmicthumb.isVisible(), "MRK_DEBUG"
			return
		assert os.path.exists(self.cter_micthumb_file_path), "MRK_DEBUG"
		
		if self.curimgmicthumbdisplay and not self.wimgmicthumb.isVisible():
			self.wimgmicthumb.show()
			self.needredisp = True
		elif not self.curimgmicthumbdisplay and self.wimgmicthumb.isVisible():
			self.wimgmicthumb.hide()
		
	def updateFFT(self, error_display = True):
		if not self.curfftdisplay: return # FFT display is unchecked
		if self.wfft == None: return # it's closed/not visible
		if self.cter_fft_file_path == None: 
			# Try directory of partres file
			pwsdir = os.path.join(os.path.dirname(self.cter_partres_file_path), "power2d")
			#if os.path.exists(pwsdir): self.cter_fft_file_path = pwsdir
			print("MRK_DEBUG: pwsdir", pwsdir)
			return
		
		if not os.path.exists(self.cter_fft_file_path):
			if self.wfft.isVisible():
				self.wfft.hide()
			if error_display and os.path.exists(os.path.dirname(self.cter_fft_file_path)):
				QtGui.QMessageBox.warning(None,"Warning","Cannot find power spectrum (%s). Please check your power spectrum directory. \n\nPower spectrum will not be shown." % (self.cter_micthumb_file_path))
			return
		assert os.path.exists(self.cter_fft_file_path), "MRK_DEBUG"
		
		# Now update the image
		fft_img = EMData(self.cter_fft_file_path) # read the image from disk
		self.wfft.set_data(fft_img)
		self.wfft.setWindowTitle("sxgui_cter - 2D FFT %s, %s" 
				% (os.path.basename(self.cter_entry_list[self.curentry][self.idx_cter_mic_name]), 
				os.path.basename(self.cter_fft_file_path)))
		
	def newFFTDisplay(self,val=None):
		"""Change FFT display status."""
		if self.curfftdisplay != val:
			# now set the new display status
			self.curfftdisplay = val
		else: 
			assert self.curfftdisplay == val, "MRK_DEBUG"
			# The status did not change, there is nothing to do 
			return
		
		if self.cter_fft_file_path == None: 
			return # no cter entry is selected
		
		if not os.path.exists(self.cter_fft_file_path):
			assert not self.wfft.isVisible(), "MRK_DEBUG"
			pwsdir = os.path.join(os.path.dirname(self.cter_partres_file_path), "power2d")
			#print("cter_fft_file_path exists: ", os.path.exists(self.cter_fft_file_path))
			QtGui.QMessageBox.warning(None,"Warning",
							 "Cannot find \"%s\" sub-directory associated with specified CTER partres file (%s). Please check your project directory. \n\nPower-spectrum display option is disabled for this session." 
							 % (pwsdir, self.cter_partres_file_path))
			#print("curfftdisplay", self.curfftdisplay)
			self.curfftdisplay = False
			self.cbfftdisplay.setValue(False)
			self.cbfftdisplay.setEnabled(False)
			return
		assert os.path.exists(self.cter_fft_file_path), "MRK_DEBUG"
		
		if self.curfftdisplay and not self.wfft.isVisible():
			self.wfft.show()
			self.needredisp = True
		elif not self.curfftdisplay and self.wfft.isVisible():
			self.wfft.hide()
		
	def newEntry(self,currow):
		"""called when a new data set is selected from the CTER partres entry list box."""
		assert self.cter_partres_file_path != None, "MRK_DEBUG"
		assert self.cter_entry_list != None, "MRK_DEBUG"
		
		# always update the current row of cter entry list 
		# to get associated micrograph path and pwrot file path 
		self.curentry = currow # row can be the same even after resorting of the cter entry list
		
		# Get associated pwrot file path of current entry
#		new_cter_pwrot_file_path = self.cter_entry_list[self.curentry][self.idx_cter_pwrot_name]
		
		# Get associated micrograph path of current entry
		new_cter_mic_file_path = self.cter_entry_list[self.curentry][self.idx_cter_mic_name]
		
		# Generate associated micthumb & pwrot file path of current entry
		mic_basename_root = os.path.splitext(os.path.basename(new_cter_mic_file_path))[0]
		new_cter_micthumb_file_path = os.path.join(os.path.dirname(self.cter_partres_file_path), "micthumb", "%s_thumb.hdf" % (mic_basename_root))
		new_cter_pwrot_file_path = os.path.join(os.path.dirname(self.cter_partres_file_path), "pwrot", "%s_rotinf.txt" % (mic_basename_root))
		new_cter_fft_file_path = os.path.join(os.path.dirname(self.cter_partres_file_path), "power2d", "%s_pws.hdf" % (mic_basename_root))
			
		# Changing row does not always change the pwrot file path after resorting of the cter entry list
		# If same, skip the following processes
		if self.cter_pwrot_file_path == new_cter_pwrot_file_path: 
			assert self.cter_micthumb_file_path == new_cter_micthumb_file_path, "MRK_DEBUG"
			assert self.cter_mic_file_path == new_cter_mic_file_path, "MRK_DEBUG"
			assert self.cter_fft_file_path == new_cter_fft_file_path, "MRK_DEBUG"
			return
		
		# now set the new item
		assert self.cter_pwrot_file_path != new_cter_pwrot_file_path, "MRK_DEBUG"
		self.cter_pwrot_file_path = new_cter_pwrot_file_path
		assert self.cter_micthumb_file_path != new_cter_micthumb_file_path, "MRK_DEBUG"
		self.cter_micthumb_file_path = new_cter_micthumb_file_path
		assert self.cter_mic_file_path != new_cter_mic_file_path , "MRK_DEBUG"
		self.cter_mic_file_path = new_cter_mic_file_path
		assert self.cter_fft_file_path != new_cter_fft_file_path , "MRK_DEBUG"
		self.cter_fft_file_path = new_cter_fft_file_path
		
		# print "MRK_DEBUG: Row No. %d (CTER partres entry No. %d) is selected from cter entry list box" % (self.curentry, self.cter_entry_list[self.curentry][self.idx_cter_id])
		
		self.ssortedid.setValue(self.curentry,True)
		
		for idx_cter in range(self.n_idx_cter):
			if self.value_map_list[idx_cter][self.idx_cter_item_widget] is not None:
				self.value_map_list[idx_cter][self.idx_cter_item_widget].setValue(self.cter_entry_list[self.curentry][idx_cter],True)
			# else:
			# 	assert self.value_map_list[idx_cter][self.idx_cter_item_widget] is None, "MRK_DEBUG"
			# 	# Do nothing
			# 	print("MRK_DEBUG: ")
			# 	print("MRK_DEBUG: No widget is assigned to this item! Ignoring this item...")
			# 	print("MRK_DEBUG: idx_cter = ", idx_cter)
			# 	print("MRK_DEBUG: self.idx_cter_item_widget = ", self.idx_cter_item_widget)
			# 	print("MRK_DEBUG: self.value_map_list[idx_cter][self.idx_cter_item_widget] = ", self.value_map_list[idx_cter][self.idx_cter_item_widget])
			# 	print("MRK_DEBUG: ")
		
		# Use blue (lower) & red (higher) fonts to indicate the value is not between applied threshold ranage
		for idx_hist in range(self.n_idx_hist):
			lower_threshold = round(self.hist_map_list[idx_hist][self.idx_hist_item_apply_threshold_lower], self.round_ndigits)
			upper_threshold = round(self.hist_map_list[idx_hist][self.idx_hist_item_apply_threshold_upper], self.round_ndigits)
			idx_cter = self.hist_map_list[idx_hist][self.idx_hist_item_idx_cter]
			#print("MRK_DEBUG", idx_hist, self.value_map_list[idx_cter][0], self.hist_map_list[idx_hist][self.idx_hist_item_apply_threshold_lower], self.hist_map_list[idx_hist][self.idx_hist_item_apply_threshold_upper])
			param_val = round(self.cter_entry_list[self.curentry][idx_cter], self.round_ndigits)
			assert self.value_map_list[idx_cter][self.idx_cter_item_widget] is not None, "MRK_DEBUG: A widget must be assigned to this item since it is a histogram item"
			if lower_threshold < upper_threshold:
				if lower_threshold <= param_val and param_val <= upper_threshold:
						self.value_map_list[idx_cter][self.idx_cter_item_widget].text.setStyleSheet("color: rgb(0,0,0);")
				elif param_val < lower_threshold:
					self.value_map_list[idx_cter][self.idx_cter_item_widget].text.setStyleSheet("color: rgb(0,0,255);")
				else:
					assert upper_threshold < param_val, "MRK_DEBUG"
					self.value_map_list[idx_cter][self.idx_cter_item_widget].text.setStyleSheet("color: rgb(255,0,0);")
			else:
				assert lower_threshold == upper_threshold, "MRK_DEBUG"
				assert lower_threshold == param_val and param_val == upper_threshold, "MRK_DEBUG"
				##print(self.value_map_list[idx_cter][self.idx_cter_item_widget].text.text())
				#print("MRK_DEBUG", idx_hist, self.value_map_list[idx_cter][0], self.hist_map_list[idx_hist][self.idx_hist_item_apply_threshold_lower], self.hist_map_list[idx_hist][self.idx_hist_item_apply_threshold_upper])
				self.value_map_list[idx_cter][self.idx_cter_item_widget].text.setStyleSheet("color: rgb(127,127,127);")
		
#		self.wfft.setWindowTitle("sxgui_cter - 2D FFT - "+fsp.split("/")[-1])
		self.wplotrotavgcoarse.setWindowTitle("sxgui_cter - Plot - %s, %s" % (os.path.basename(self.cter_entry_list[self.curentry][self.idx_cter_mic_name]), os.path.basename(self.cter_pwrot_file_path)))
		self.wplotrotavgfine.setWindowTitle("sxgui_cter - Plot Zoom- %s, %s" % (os.path.basename(self.cter_entry_list[self.curentry][self.idx_cter_mic_name]), os.path.basename(self.cter_pwrot_file_path)))
		
#		self.cter_mic_file_path = new_cter_mic_file_path
	
		self.needredisp = True
	
	def updateEntrySelect(self, entry):
		"""called when check status of an cter entry in list box is changed."""
		assert self.cter_partres_file_path != None, "MRK_DEBUG"
		assert self.cter_entry_list != None, "MRK_DEBUG"
		
		newSelect = 1
		if entry.checkState() == Qt.Unchecked:
			newSelect = 0
		entry_row = self.lbentry.row(entry)
		self.cter_entry_list[entry_row][self.idx_cter_select] = newSelect
		
		if self.curentry == entry_row:
			assert self.value_map_list[self.idx_cter_select][self.idx_cter_item_widget] is not None, "MRK_DEBUG: It should be impossible to select item without any associated widgets!!!"
			self.value_map_list[self.idx_cter_select][self.idx_cter_item_widget].setValue(self.cter_entry_list[self.curentry][self.idx_cter_select],True)
			
		self.updateUncheckCounts()
	
	def updateUncheckCounts(self):
		"""called whenever checked status of cter entries change."""
		assert self.cter_partres_file_path != None, "MRK_DEBUG"
		assert self.cter_entry_list != None, "MRK_DEBUG"
		
		assert len(self.cter_entry_list) > 0, "MRK_DEBUG"
		n_entry = len(self.cter_entry_list)
		uncheck_counts  = n_entry
		for cter_entry in self.cter_entry_list:
			uncheck_counts -= cter_entry[self.idx_cter_select]
		assert uncheck_counts >= 0 and uncheck_counts <= n_entry, "MRK_DEBUG"
		
		self.vbuncheckcounts.setValue(uncheck_counts,True)
		self.vbuncheckratio.setValue(float(uncheck_counts)/n_entry,True)
	
	def reapplySort(self,item = None):
		"""Called when reapply button is clicked."""
		if self.cter_partres_file_path == None: return # no cter ctf file is selected
		if self.cter_entry_list == None: return # no cter ctf file is selected
		
		self.updateEntryList()
	
	def newSort(self,cursortidx):
		"""Sort CTER partres entries by selected parameter values."""
		if self.cursortidx == cursortidx: return
		
		# now set the new item
		self.cursortidx = cursortidx
		
		if self.cter_partres_file_path == None: return # no cter ctf file is selected
		if self.cter_entry_list == None: return # no cter ctf file is selected
		
		self.updateEntryList()
	
	def newSortOrder(self, sortorder):
		"""Change sorting order of CTER partres entries."""
		if self.cursortorder == sortorder: return
		
		# now set the new status
		self.cursortorder = sortorder
		
		if self.cter_partres_file_path == None: return # no cter ctf file is selected
		if self.cter_entry_list == None: return # no cter ctf file is selected
		
		self.updateEntryList()
	
	def newSortSelect(self, sortselect):
		"""Change sort select status of CTER partres entries."""
		if self.cursortselect == sortselect: return
		
		# now set the new status
		self.cursortselect = sortselect
		
		if self.cter_partres_file_path == None: return # no cter ctf file is selected
		if self.cter_entry_list == None: return # no cter ctf file is selected
		
		self.updateEntryList()
	
	def newThresholdLower(self):
		threshold_lower = round(self.hist_map_list[self.curhistidx][self.idx_hist_item_unapply_widget_lower].getValue(), self.round_ndigits)
		if threshold_lower < round(self.hist_map_list[self.curhistidx][self.idx_hist_item_val_min], self.round_ndigits):
			threshold_lower = round(self.hist_map_list[self.curhistidx][self.idx_hist_item_val_min], self.round_ndigits)
			self.hist_map_list[self.curhistidx][self.idx_hist_item_unapply_widget_lower].setValue(threshold_lower)
		elif threshold_lower > round(self.hist_map_list[self.curhistidx][self.idx_hist_item_val_max], self.round_ndigits):
			threshold_lower = round(self.hist_map_list[self.curhistidx][self.idx_hist_item_val_max], self.round_ndigits)
			self.hist_map_list[self.curhistidx][self.idx_hist_item_unapply_widget_lower].setValue(threshold_lower)
		# else: # Do nothing
		
		# now set the new threshold
		self.hist_map_list[self.curhistidx][self.idx_hist_item_unapply_threshold_lower] = threshold_lower
		
		self.needredisp=True
	
	def newThresholdUpper(self):
		threshold_upper = round(self.hist_map_list[self.curhistidx][self.idx_hist_item_unapply_widget_upper].getValue(), self.round_ndigits)
		if threshold_upper < round(self.hist_map_list[self.curhistidx][self.idx_hist_item_val_min], self.round_ndigits):
			threshold_upper = round(self.hist_map_list[self.curhistidx][self.idx_hist_item_val_min], self.round_ndigits)
			self.hist_map_list[self.curhistidx][self.idx_hist_item_unapply_widget_upper].setValue(threshold_upper)
		elif threshold_upper > round(self.hist_map_list[self.curhistidx][self.idx_hist_item_val_max], self.round_ndigits):
			threshold_upper = round(self.hist_map_list[self.curhistidx][self.idx_hist_item_val_max], self.round_ndigits)
			self.hist_map_list[self.curhistidx][self.idx_hist_item_unapply_widget_upper].setValue(threshold_upper)
		# else: # Do nothing
		
		# now set the new threshold
		self.hist_map_list[self.curhistidx][self.idx_hist_item_unapply_threshold_upper] = threshold_upper
		
		self.needredisp=True
	
	def newHistogramRow(self,currow):
		"called when a new row is selected from the Histogram list box"
		
		if self.curhistidx == currow: return
		
		#print('newHistogramRow: old %s %s, new %s %s' % (self.curhistidx, self.value_map_list[self.curhistidx][self.idx_cter_item_label],
			#currow, self.value_map_list[currow][self.idx_cter_item_label]))
		#print('newHistogramRow: self.whistparam.isVisible', self.whistparam.isVisible())
		#print('newHistogramRow: self.wscatterparam.isVisible', self.wscatterparam.isVisible())
			
		# Disable old item
		if self.curthresholdcontrol == self.idx_threshold_control_lower:
			self.hist_map_list[self.curhistidx][self.idx_hist_item_unapply_widget_lower].setEnabled(False)
		elif self.curthresholdcontrol == self.idx_threshold_control_upper:
			self.hist_map_list[self.curhistidx][self.idx_hist_item_unapply_widget_upper].setEnabled(False)
		else:
			assert self.curthresholdcontrol == self.idx_threshold_control_edit_only, "MRK_DEBUG"
			self.hist_map_list[self.curhistidx][self.idx_hist_item_unapply_widget_lower].setEnabled(False)
			self.hist_map_list[self.curhistidx][self.idx_hist_item_unapply_widget_upper].setEnabled(False)
		
		# now set the new item and enable it
		self.curhistidx=currow
		# print "MRK_DEBUG: Row No. %d is selected from histogram list box" % (self.curhistidx)
		
		# Check if the all selected parameter values are same 
		if self.hist_map_list[self.curhistidx][self.idx_hist_item_val_min] == self.hist_map_list[self.curhistidx][self.idx_hist_item_val_max]:
			idx_cter = self.hist_map_list[self.curhistidx][self.idx_hist_item_idx_cter]
			param_label = self.value_map_list[idx_cter][self.idx_cter_item_label]
			#self.curhistdisable=True
			self.curhistogramdisplay = False
			if self.whistparam.isVisible():
				self.whistparam.hide()
			if self.wscatterparam.isVisible():
				self.wscatterparam.hide()
			QtGui.QMessageBox.information(self, 
								 "Information","All entries have the same selected parameter values (%s). \n\nParameter Histogram & Plot will not be shown." 
								 % (param_label))
		else:
			if self.curthresholdcontrol == self.idx_threshold_control_lower:
				self.hist_map_list[self.curhistidx][self.idx_hist_item_unapply_widget_lower].setEnabled(True)
			elif self.curthresholdcontrol == self.idx_threshold_control_upper:
				self.hist_map_list[self.curhistidx][self.idx_hist_item_unapply_widget_upper].setEnabled(True)
			else:
				assert self.curthresholdcontrol == self.idx_threshold_control_edit_only, "MRK_DEBUG"
				self.hist_map_list[self.curhistidx][self.idx_hist_item_unapply_widget_lower].setEnabled(True)
				self.hist_map_list[self.curhistidx][self.idx_hist_item_unapply_widget_upper].setEnabled(True)
			
			idx_cter = self.hist_map_list[self.curhistidx][self.idx_hist_item_idx_cter]
			param_label = self.value_map_list[idx_cter][self.idx_cter_item_label]
			self.whistparam.setWindowTitle("sxgui_cter - %s Histogram" % (param_label))
			
			if self.cursyncsort == True:
				idx_sort = self.hist_map_list[self.curhistidx][self.idx_hist_item_idx_sort]
				if (idx_sort != self.cursortidx):
					self.newSort(idx_sort)
					self.ssort.setCurrentIndex(idx_sort)
				# else: assert idx_sort == self.cursortidx, "MRK_DEBUG" # Do nothing
			# else: assert self.cursyncsort == False, "MRK_DEBUG" # Do nothing
			
			if self.cter_partres_file_path == None: return # no cter ctf file is selected
			if self.cter_entry_list == None: return # no cter ctf file is selected
			
			#self.curhistdisable=False
			
			#### (These statements will open the corresponding windows if they are not already open.)
			#if not self.whistparam.isVisible():
				#self.whistparam.show()
			#if not self.wscatterparam.isVisible():
				#self.wscatterparam.show()
			
			# NOTE: 2016/01/03 Toshio Moriya
			# Force update related plots for scaling delay...
			self.updateHist()
			self.updatePlotParam()
			if self.cursyncsort == True:
				self.updateEntryList()
			
			self.needredisp = True
	
	def newThresholdControl(self, currow):
		"called when a new row is selected from the Threshold Control list box"
		
		if self.curthresholdcontrol == currow: return
		
		# Disable old item
		if self.curthresholdcontrol == self.idx_threshold_control_lower:
			self.hist_map_list[self.curhistidx][self.idx_hist_item_unapply_widget_lower].setEnabled(False)
		elif self.curthresholdcontrol == self.idx_threshold_control_upper:
			self.hist_map_list[self.curhistidx][self.idx_hist_item_unapply_widget_upper].setEnabled(False)
		else:
			assert self.curthresholdcontrol == self.idx_threshold_control_edit_only, "MRK_DEBUG"
			self.hist_map_list[self.curhistidx][self.idx_hist_item_unapply_widget_lower].setEnabled(False)
			self.hist_map_list[self.curhistidx][self.idx_hist_item_unapply_widget_upper].setEnabled(False)
		
		# now set the new item and enalble it
		self.curthresholdcontrol=currow
		# print "MRK_DEBUG: Row No. %d is selected from histogram list box" % (self.curhistidx)
		
		if self.curthresholdcontrol == self.idx_threshold_control_lower:
			self.hist_map_list[self.curhistidx][self.idx_hist_item_unapply_widget_lower].setEnabled(True)
		elif self.curthresholdcontrol == self.idx_threshold_control_upper:
			self.hist_map_list[self.curhistidx][self.idx_hist_item_unapply_widget_upper].setEnabled(True)
		else:
			assert self.curthresholdcontrol == self.idx_threshold_control_edit_only, "MRK_DEBUG"
			self.hist_map_list[self.curhistidx][self.idx_hist_item_unapply_widget_lower].setEnabled(True)
			self.hist_map_list[self.curhistidx][self.idx_hist_item_unapply_widget_upper].setEnabled(True)
	
	def newSyncSort(self, syncsort):
		"""Change sync sort enable state."""
		if self.cursyncsort == syncsort: return
		
		# now set the new status
		self.cursyncsort = syncsort
		
		self.ssort.setEnabled(not self.cursyncsort)
		
		if self.cursyncsort == True:
			idx_sort = self.hist_map_list[self.curhistidx][self.idx_hist_item_idx_sort]
			if (idx_sort != self.cursortidx):
				self.newSort(idx_sort)
				self.ssort.setCurrentIndex(idx_sort)
			# else: assert idx_sort == self.cursortidx, "MRK_DEBUG" # Do nothing
		# else: assert self.cursyncsort == False, "MRK_DEBUG" # Do nothing
	
	def newEntryPerBin(self,curentryperbin):
		if self.curentryperbin == curentryperbin: return
		
		# now set the new entry per bin
		self.curentryperbin = curentryperbin
		
		if self.cter_partres_file_path == None: return # no cter ctf file is selected
		if self.cter_entry_list == None: return # no cter ctf file is selected
		
		# NOTE: 2016/01/03 Toshio Moriya
		# Force update related plots for scaling delay...
		self.updateHist()
		self.updatePlotParam()
		
		self.needredisp = True
	
	def applyAllThresholds(self,item = None):
		if self.cter_partres_file_path == None: return # no cter ctf file is selected
		if self.cter_entry_list == None: return # no cter ctf file is selected
		
		reply = QtGui.QMessageBox.question(self, "Warning", "Applying all threshold setting will wipe the previous selection states including manual setting. Do you really want to continue?", QtGui.QMessageBox.Yes | QtGui.QMessageBox.No, QtGui.QMessageBox.No)
		if reply == QtGui.QMessageBox.No:
			return
		
		# Set set the select status of all cter entries based on the threshold values
		for cter_entry in self.cter_entry_list:
			new_select_state = 1
			for idx_hist in range(self.n_idx_hist):
				threshold_lower = round(self.hist_map_list[idx_hist][self.idx_hist_item_unapply_threshold_lower], self.round_ndigits)
				threshold_upper = round(self.hist_map_list[idx_hist][self.idx_hist_item_unapply_threshold_upper], self.round_ndigits)
				self.hist_map_list[idx_hist][self.idx_hist_item_apply_threshold_lower] = threshold_lower
				self.hist_map_list[idx_hist][self.idx_hist_item_apply_threshold_upper] = threshold_upper
				self.hist_map_list[idx_hist][self.idx_hist_item_apply_widget_lower].setValue(threshold_lower)
				self.hist_map_list[idx_hist][self.idx_hist_item_apply_widget_upper].setValue(threshold_upper)
				idx_cter = self.hist_map_list[idx_hist][self.idx_hist_item_idx_cter]
				param_val = round(cter_entry[idx_cter], self.round_ndigits)
				if param_val < threshold_lower or threshold_upper < param_val:
					# print "MRK_DEBUG: Param #%d diselected entry #%04d with (param_val, threshold_lower, threshold_upper) = (%1.15g, %1.15g, %1.15g)" % (idx_hist, idx_cter, param_val, threshold_lower, threshold_upper)
					new_select_state = 0
				# else: # Do nothing
			cter_entry[self.idx_cter_select] = new_select_state
			
		self.updateEntryList()
		self.updateUncheckCounts()
	
	def newThresholdSet(self, currow):
		"called when a new row is selected from the Threshold Set list box"
		
		if self.curthresholdset == currow: return
		# now set the new item and enalble it
		self.curthresholdset=currow
	
	def writeThresholdSet(self, file_path_out, idx_thresholdset):
		assert self.cter_partres_file_path != None, "MRK_DEBUG"
		assert self.cter_entry_list != None, "MRK_DEBUG"
		
		file_out = open(file_path_out,"w")
		
		# Write lines to check consistency upon loading
		file_out.write("# @@@@@ gui_cter thresholds - ")
		# file_out.write(EMANVERSION + " (CVS" + CVSDATESTAMP[6:-2] +")")
		file_out.write(EMANVERSION + " (GITHUB: " + DATESTAMP +")" )
		file_out.write(" @@@@@ \n")
		file_out.write("# Associated CTER Partres File == %s\n" % (self.cter_partres_file_path))
		file_out.write("# Saved Threshold Set == %s\n" % (self.thresholdset_map_list[idx_thresholdset][self.idx_thresholdset_item_label]))
		file_out.write("# [Parameter Id] [Parameter Name] [Lower Threshold] [Upper Threshold]\n")
		
		# Assigne the index of target threshold values
		idx_threshold_lower = self.idx_hist_item_unapply_threshold_lower
		idx_threshold_upper = self.idx_hist_item_unapply_threshold_upper
		if idx_thresholdset == self.idx_thresholdset_applied:
			idx_threshold_lower = self.idx_hist_item_apply_threshold_lower
			idx_threshold_upper = self.idx_hist_item_apply_threshold_upper
		
		for idx_hist in range(self.n_idx_hist):
			map_entry = self.hist_map_list[idx_hist]
			idx_cter = map_entry[self.idx_hist_item_idx_cter]
			param_label = self.value_map_list[idx_cter][self.idx_cter_item_label]
			# NOTE: 2016/01/26 Toshio Moriya
			# Use the precision for double to minimise precision loss by save & load operations 
			file_out.write("%2d %s == %1.15g %1.15g \n" % (idx_hist, param_label, round(map_entry[idx_threshold_lower], self.round_ndigits), round(map_entry[idx_threshold_upper], self.round_ndigits)))
		
		file_out.close()
	
	def readThresholdSet(self, file_path_in, idx_thresholdset):
		assert self.cter_partres_file_path != None, "MRK_DEBUG"
		assert self.cter_entry_list != None, "MRK_DEBUG"
		
		file_in = open(file_path_in,"r")
		
		# Check if this parameter file is threshold
		line_in = file_in.readline()
		if line_in.find("@@@@@ gui_cter thresholds") != -1:
			# loop through the rest of lines
			for line_in in file_in:
				if line_in[0] == "#":
					continue
					
				tokens_in = line_in.split("==")
				assert len(tokens_in) == 2, "MRK_DEBUG"
				tokens_label = tokens_in[0].split()
				assert len(tokens_label) >= 2, "MRK_DEBUG"
				idx_hist = int(tokens_label[0]) # index of hist_map_list
				map_entry = self.hist_map_list[idx_hist]
				tokens_val = tokens_in[1].split()
				assert len(tokens_val) == 2, "MRK_DEBUG"
				threshold_lower = round(float(tokens_val[0]), self.round_ndigits)
				threshold_upper = round(float(tokens_val[1]), self.round_ndigits)
				map_entry[self.idx_hist_item_unapply_threshold_lower] = threshold_lower
				map_entry[self.idx_hist_item_unapply_threshold_upper] = threshold_upper
				map_entry[self.idx_hist_item_unapply_widget_lower].setValue(threshold_lower)
				map_entry[self.idx_hist_item_unapply_widget_upper].setValue(threshold_upper)
				self.newThresholdLower()
				self.newThresholdUpper()
			
			if idx_thresholdset == self.idx_thresholdset_applied:
				self.applyAllThresholds()
		else:
			QtGui.QMessageBox.warning(self, "Warning", "The specified file is not threshold file.")
		
		file_in.close()
	
	def saveThresholdSet(self,item = None):
		if self.cter_partres_file_path == None: return # no cter ctf file is selected
		if self.cter_entry_list == None: return # no cter ctf file is selected
		
		title_string = "Save %s Thresholds" % self.thresholdset_map_list[self.curthresholdset][self.idx_thresholdset_item_label]
		file_path_out = str(QtGui.QFileDialog.getSaveFileName(self, title_string, options = QtGui.QFileDialog.DontUseNativeDialog))
		if file_path_out == "": return
		
		self.writeThresholdSet(os.path.relpath(file_path_out), self.curthresholdset)
	
	def loadThresholdSet(self,item = None):
		if self.cter_partres_file_path == None: return # no cter ctf file is selected
		if self.cter_entry_list == None: return # no cter ctf file is selected
		
		reply = QtGui.QMessageBox.question(self, "Warning", "Loading thresholds will wipe the previous threshold setting. Do you really want to continue?", QtGui.QMessageBox.Yes | QtGui.QMessageBox.No, QtGui.QMessageBox.No)
		if reply == QtGui.QMessageBox.No:
			return
		
		title_string = "Load %s Thresholds" % self.thresholdset_map_list[self.curthresholdset][self.idx_thresholdset_item_label]
		file_path_in = str(QtGui.QFileDialog.getOpenFileName(self, title_string, options = QtGui.QFileDialog.DontUseNativeDialog))
		if file_path_in == "": return
		
		self.readThresholdSet(os.path.relpath(file_path_in), self.curthresholdset)
	
	def saveSelection(self,item = None):
		if self.cter_partres_file_path == None: return # no cter ctf file is selected
		if self.cter_entry_list == None: return # no cter ctf file is selected
		
		assert os.path.basename(self.cter_partres_file_path).find("partres") != -1, "MRK_DEBUG"
		assert self.cter_partres_file_path[-1*len(".txt"):] == ".txt", "MRK_DEBUG"
#		assert os.path.dirname(self.cter_partres_file_path)[-1*len("partres"):] == "partres", "MRK_DEBUG"
		
		file_suffix = self.vfilesuffix.getValue()
		file_path_out_select = os.path.join(os.path.dirname(self.cter_partres_file_path), "%s_partres_select.txt" % (file_suffix))
		file_path_out_discard = os.path.join(os.path.dirname(self.cter_partres_file_path), "%s_partres_discard.txt" % (file_suffix))
		file_path_out_mic_select = os.path.join(os.path.dirname(self.cter_partres_file_path), "%s_micrographs_select.txt" % (file_suffix))
		file_path_out_mic_discard = os.path.join(os.path.dirname(self.cter_partres_file_path), "%s_micrographs_discard.txt" % (file_suffix))
		file_path_out_thresholds = os.path.join(os.path.dirname(self.cter_partres_file_path), "%s_thresholds.txt" % (file_suffix))
		
		existing_file_path = None
		if os.path.exists(file_path_out_select):
			existing_file_path = file_path_out_select
		elif os.path.exists(file_path_out_discard):
			existing_file_path = file_path_out_discard
		elif os.path.exists(file_path_out_mic_select):
			existing_file_path = file_path_out_mic_select
		elif os.path.exists(file_path_out_mic_discard):
			existing_file_path = file_path_out_mic_discard
		elif os.path.exists(file_path_out_thresholds):
			existing_file_path = file_path_out_thresholds
		# else: # Do nothing
		
		if existing_file_path != None:
			reply = QtGui.QMessageBox.question(self, "Warning", "The file (%s) already exists. Do you want to overwrite the file?" % (existing_file_path), QtGui.QMessageBox.Yes | QtGui.QMessageBox.No, QtGui.QMessageBox.No)
			if reply == QtGui.QMessageBox.No:
				return
		
		# Save selection in CTER Format
		file_out_select = open(file_path_out_select,"w")
		file_out_discard = open(file_path_out_discard,"w")
		
		save_cter_entry_list = sorted(self.cter_entry_list, key=lambda x: x[self.idx_cter_id])
# 		idx_cter_ignore_list = [self.idx_cter_id, self.idx_cter_pwrot_name, self.idx_cter_error_ctf, self.idx_cter_select]
		idx_cter_ignore_list = [self.idx_cter_id, self.idx_cter_select]
###		if self.is_enable_max_power == True: 
###			idx_cter_ignore_list.append(self.idx_cter_max_power)
		for cter_entry in save_cter_entry_list:
			file_out = file_out_select
			if cter_entry[self.idx_cter_select] == 0:
				file_out = file_out_discard
			# else: assert cter_entry[self.idx_cter_select] == 1, "MRK_DEBUG" # do nothing
			
			for idx_cter in range(self.n_idx_cter):
				if idx_cter in idx_cter_ignore_list:
					# Do nothing
					continue
				elif idx_cter in [self.idx_cter_mic_name]:
					file_out.write("  %s" % cter_entry[idx_cter])
				else:
					file_out.write("  %12.5g" % cter_entry[idx_cter])
			file_out.write("\n")
		
		file_out_select.close()
		file_out_discard.close()
		
		# Save selection in Micrograph List Format
		file_out_mic_select = open(file_path_out_mic_select,"w")
		file_out_mic_discard = open(file_path_out_mic_discard,"w")
		
		for cter_entry in save_cter_entry_list:
			file_out = file_out_mic_select
			if cter_entry[self.idx_cter_select] == 0:
				file_out = file_out_mic_discard
			# else: assert cter_entry[self.idx_cter_select] == 1, "MRK_DEBUG" # do nothing
			file_out.write("  %s\n" % os.path.basename(cter_entry[self.idx_cter_mic_name]))
		
		file_out_mic_select.close()
		file_out_mic_discard.close()
		
		# Save the associated applied threshold 
		self.writeThresholdSet(file_path_out_thresholds, self.idx_thresholdset_applied) 
		
		QtGui.QMessageBox.information(self, "Information","The following files are saved in %s:\n\nCTER Partres File - Selected: %s\n\nCTER Partres FIle - Discarded: %s\n\nMicrograph - Selected: %s\n\nMicrograph - Discarded: %s\n\nApplied Threshold Set: %s" 
								% (os.path.dirname(self.cter_partres_file_path), os.path.basename(file_path_out_select), os.path.basename(file_path_out_discard), os.path.basename(file_path_out_mic_select), os.path.basename(file_path_out_mic_discard), os.path.basename(file_path_out_thresholds)))
	
	def timeOut(self):
		if self.busy: return
		
		# Redisplay before spawning thread for more interactive display
		if self.needredisp:
			try:
				self.redisplay()
			except:
				print("Received unexpected exception from redisplay() in timeOut(): ")
				exc_type, exc_value, exc_traceback = sys.exc_info()
				traceback.print_exception(exc_type, exc_value, exc_traceback)
				# MRK_NOTE: 2015/12/17 Toshio Moriya
				# Another way to print out exception info...
				# lines = traceback.format_exception(exc_type, exc_value, exc_traceback)
				# print ''.join('!! ' + line for line in lines)
				pass
	
	def redisplay(self):
		self.needredisp=False
		self.busy=True
		
		if self.cter_entry_list != None:
			is_child_shown = False
			if not self.wfft.isVisible() and self.is_wfft_minimized:
				self.wfft.show()
				is_child_shown = True
			if not self.wimgmicthumb.isVisible() and self.curimgmicthumbdisplay and not self.is_wimgmicthumb_minimized and os.path.exists(self.cter_micthumb_file_path):
				self.wimgmicthumb.show()
				is_child_shown = True
			if not self.wplotrotavgcoarse.isVisible() and self.curplotrotavgdisplay and not self.is_wplotrotavgcoarse_minimized and os.path.exists(self.cter_pwrot_file_path):
				self.wplotrotavgcoarse.show()
				is_child_shown = True
			if not self.wplotrotavgfine.isVisible() and self.curplotrotzoomdisplay and not self.is_wplotrotavgfine_minimized and os.path.exists(self.cter_pwrot_file_path):
				self.wplotrotavgfine.show()
				is_child_shown = True
			#if not self.whistparam.isVisible() and not self.curhistdisable and not self.is_whistparam_minimized:
			if not self.whistparam.isVisible() and self.curhistogramdisplay and not self.is_whistparam_minimized:
				self.whistparam.show()
				is_child_shown = True
			#if not self.wscatterparam.isVisible() and not self.curhistdisable and not self.is_wscatterparam_minimized:
			if not self.wscatterparam.isVisible() and self.curscatterdisplay and not self.is_wscatterparam_minimized:
				self.wscatterparam.show()
				is_child_shown = True
			if is_child_shown == True:
				self.raise_()
				self.activateWindow()
				
		self.updateImgMicThumb()
		self.updateHist()
		self.updatePlotParam()
		self.updatePlotCurves()
		self.updateFFT()
		
		self.busy = False
		
#	def entryKey(self,event):
#		if event.key()>=Qt.Key_0 and event.key()<=Qt.Key_9 :
#			q=int(event.key())-Qt.Key_0
#			self.squality.setValue(q)
#		elif event.key() == Qt.Key_Left:
#			self.sdef.setValue(self.sdef.getValue()-0.03)
#		elif event.key() == Qt.Key_Right:
#			self.sdef.setValue(self.sdef.getValue()+0.03)
#		elif event.key()==Qt.Key_I :
#			print "MRK_DEBUG: Not used now"
#			self.doImport()
#		elif event.key()==Qt.Key_U :
#			print "MRK_DEBUG: Not used now"
#			self.unImport()
#		elif event.key()==Qt.Key_F :
#			print "MRK_DEBUG: Not used now"
#			self.doRefit()
	
	def eventFilter(self, source, event):
		if event.type() == QtCore.QEvent.WindowStateChange:
			# print "MRK_DEBUG: Hello QtCore.QEvent.WindowStateChange"
			if self.windowState() & QtCore.Qt.WindowMinimized:
				# print "MRK_DEBUG: sxgui main window has minimized"
				assert self.isMinimized() == True, "MRK_DEBUG"
				#
				# NOTE: 2016/03/09 Toshio Moriya
				# Minimize icon button of child window should be disabled
				#
				if self.cter_entry_list != None:
					if self.wfft.isVisible() and not self.is_wfft_minimized:
						self.wfft.hide()
						self.is_wfft_minimized = True
					if self.wimgmicthumb.isVisible() and not self.is_wimgmicthumb_minimized:
						assert self.curimgmicthumbdisplay, "MRK_DEBUG"
						self.wimgmicthumb.hide()
						self.is_wimgmicthumb_minimized = True
					if self.wplotrotavgcoarse.isVisible() == True and not self.is_wplotrotavgcoarse_minimized:
						assert self.curplotrotavgdisplay, "MRK_DEBUG"
						self.wplotrotavgcoarse.hide()
						self.is_wplotrotavgcoarse_minimized = True
					if self.wplotrotavgfine.isVisible() and not self.is_wplotrotavgfine_minimized:
						assert self.curplotrotzoomdisplay, "MRK_DEBUG"
						self.wplotrotavgfine.hide()
						self.is_wplotrotavgfine_minimized = True
					if self.whistparam.isVisible() == True and self.is_whistparam_minimized == False:
						#assert self.curhistdisable == False, "MRK_DEBUG"
						assert self.curhistogramdisplay == True, "MRK_DEBUG"
						self.whistparam.hide()
						self.is_whistparam_minimized = True
					if self.wscatterparam.isVisible() == True and self.is_wscatterparam_minimized == False:
						#assert self.curhistdisable == False, "MRK_DEBUG"
						assert self.curhistogramdisplay == True, "MRK_DEBUG"
						self.wscatterparam.hide()
						self.is_wscatterparam_minimized = True
			else:
				# print "MRK_DEBUG: sxgui main window has not minimized"
				assert self.isMinimized() == False, "MRK_DEBUG"
				#
				# NOTE: 2016/03/09 Toshio Moriya
				# Minimize icon button of child window should be disabled
				#
				if self.cter_entry_list != None:
					if self.is_wfft_minimized == True:
						assert self.wfft.isVisible() == False and self.curfftdisplay == True, "MRK_DEBUG"
						self.wfft.show()
						self.is_wfft_minimized = False
					if self.is_wimgmicthumb_minimized:
						assert not self.wimgmicthumb.isVisible() and self.curimgmicthumbdisplay, "MRK_DEBUG"
						self.wimgmicthumb.show()
						self.is_wimgmicthumb_minimized = False
					if self.is_wplotrotavgcoarse_minimized:
						assert not self.wplotrotavgcoarse.isVisible() and self.curplotrotavgdisplay, "MRK_DEBUG"
						self.wplotrotavgcoarse.show()
						self.is_wplotrotavgcoarse_minimized = False
					if self.is_wplotrotavgfine_minimized:
						assert not self.wplotrotavgfine.isVisible() and self.curplotrotzoomdisplay, "MRK_DEBUG"
						self.wplotrotavgfine.show()
						self.is_wplotrotavgfine_minimized = False
					if self.is_whistparam_minimized:
						#assert not self.whistparam.isVisible() and not self.curhistdisable, "MRK_DEBUG"
						assert not self.whistparam.isVisible() and self.curhistogramdisplay, "MRK_DEBUG"
						self.whistparam.show()
						self.is_whistparam_minimized = False
					if self.is_wscatterparam_minimized:
						#assert not self.wscatterparam.isVisible() and not self.curhistdisable, "MRK_DEBUG"
						assert not self.wscatterparam.isVisible() and self.curhistogramdisplay, "MRK_DEBUG"
						self.wscatterparam.show()
						self.is_wscatterparam_minimized = False
				if self.isVisible() : # Depends on timing at startup, this can happen?!!
					self.raise_()
					self.activateWindow()
		elif event.type() == QtCore.QEvent.WindowActivate:
			# print "MRK_DEBUG: sxgui main window has gained focus (become active)"
			# 
			# NOTE: 2016/03/08 Toshio Moriya
			# To raise EMGLWidget (SXPlot2DWidget and EMImage2DWidget) window,
			# we have to go through its qt_parent attribute to call raise_()...
			# 
			if self.cter_entry_list != None:
				if self.wfft.isVisible() == True:
					self.wfft.qt_parent.raise_()
				if self.wimgmicthumb.isVisible() == True:
					self.wimgmicthumb.qt_parent.raise_()
				if self.wplotrotavgcoarse.isVisible() == True:
					self.wplotrotavgcoarse.qt_parent.raise_()
				if self.wplotrotavgfine.isVisible() == True:
					self.wplotrotavgfine.qt_parent.raise_()
				if self.whistparam.isVisible() == True:
					self.whistparam.qt_parent.raise_()
				if self.wscatterparam.isVisible() == True:
					self.wscatterparam.qt_parent.raise_()
				assert self.isVisible(), "MRK_DEBUG"
				self.raise_()
		# elif event.type()== QtCore.QEvent.WindowDeactivate:
		# 	print "MRK_DEBUG: sxgui main window has lost focus (beome deactive)"
		# elif event.type()== QtCore.QEvent.FocusIn:
		# 	print "MRK_DEBUG: sxgui main has gained keyboard focus"
		# elif event.type()== QtCore.QEvent.FocusOut:
		# 	print "MRK_DEBUG: sxgui main has lost keyboard focus"
		
		return super(SXGuiCter, self).eventFilter(source, event)
	
	def closeEvent(self,event):
		# Save current window layout
		E2saveappwin("sxgui_cter","main",self)
		if self.cter_entry_list != None:
			E2saveappwin("sxgui_cter","fft",self.wfft.qt_parent)
			E2saveappwin("sxgui_cter","imgmicthumb",self.wimgmicthumb.qt_parent)
			E2saveappwin("sxgui_cter","plotcoarse",self.wplotrotavgcoarse.qt_parent)
			E2saveappwin("sxgui_cter","plotfine",self.wplotrotavgfine.qt_parent)
			E2saveappwin("sxgui_cter","plotparam",self.wscatterparam.qt_parent)
			E2saveappwin("sxgui_cter","histparam",self.whistparam.qt_parent)
		
		# close all child windows
		if self.wfft: self.wfft.close()
		if self.wimgmicthumb: self.wimgmicthumb.close()
		if self.wplotrotavgcoarse: self.wplotrotavgcoarse.close()
		if self.wplotrotavgfine: self.wplotrotavgfine.close()
		if self.whistparam: self.whistparam.close()
		if self.wscatterparam: self.wscatterparam.close()
		
		event.accept()
		# NOTE: Toshio Moriya 2017/11/23
		# The following code seems to be printing "None" to standard output
		# upon exiting the application. But why???
		# print ("MRK_DEBUG: calling QtGui.qApp.exit(0)")
		QtGui.qApp.exit(0)
		
	def updatePlotVisibility(self,val=None):
		if self.wplotrotavgcoarse == None: return # it's closed/not visible
		if self.wplotrotavgfine == None: return # it's closed/not visible
		if self.cter_pwrot_file_path == None: return # no cter entry is selected
		
		# REDUNDANT?
		#for idx_graph in range(self.n_idx_graph):
			#item_widget = self.graph_map_list[idx_graph][self.idx_graph_item_widget]
			#name = self.graph_map_list[idx_graph][self.idx_graph_item_name]
		
		for idx_graph in range(self.n_idx_graph):
			item_widget = self.graph_map_list[idx_graph][self.idx_graph_item_widget]
			name = self.graph_map_list[idx_graph][self.idx_graph_item_name]
			
			#print('updatePlotVisibility', name, idx_graph, self.wplotrotavgfine.visibility[name])
			if self.wplotrotavgfine.visibility[name] != item_widget.getValue():
				self.wplotrotavgfine.visibility[name] = item_widget.getValue()
				self.wplotrotavgfine.full_refresh()
				self.wplotrotavgfine.updateGL()
				
			if self.wplotrotavgcoarse.visibility[name] != item_widget.getValue():
				self.wplotrotavgcoarse.visibility[name] = item_widget.getValue()
				self.wplotrotavgcoarse.full_refresh()
				self.wplotrotavgcoarse.updateGL()
	
	def newPlotFixScale(self,curplotfixscale):
		if self.curplotfixscale == curplotfixscale: return
		
		# now set the new entry per bin
		self.curplotfixscale = curplotfixscale
		
		if self.cter_partres_file_path == None: return # no cter ctf file is selected
		if self.cter_entry_list == None: return # no cter ctf file is selected
		
		# NOTE: 2016/01/03 Toshio Moriya
		# Force update related plots for scaling delay...
		# self.updatePlotParam()
		
		self.needredisp = True
	
	def refreshGraphs(self,item):
		self.needredisp = True
	
	def plotparammouseup(self,event):
		if self.curthresholdcontrol == self.idx_threshold_control_edit_only:
			return
		
		# Swap control if shift button is pressed
		is_not_reverse_control = True
		modifiers = event.modifiers()
		if modifiers & QtCore.Qt.ShiftModifier:
			is_not_reverse_control = False
		
		plot_x,plot_y=self.wscatterparam.scr2plot(event.x(),event.y())
		if self.curthresholdcontrol == self.idx_threshold_control_lower:
			if is_not_reverse_control:
				self.hist_map_list[self.curhistidx][self.idx_hist_item_unapply_widget_lower].setValue(plot_y)
				self.newThresholdLower()
			else:
				self.hist_map_list[self.curhistidx][self.idx_hist_item_unapply_widget_upper].setValue(plot_y)
				self.newThresholdUpper()
		else:
			assert self.curthresholdcontrol == self.idx_threshold_control_upper, "MRK_DEBUG"
			if is_not_reverse_control:
				self.hist_map_list[self.curhistidx][self.idx_hist_item_unapply_widget_upper].setValue(plot_y)
				self.newThresholdUpper()
			else:
				self.hist_map_list[self.curhistidx][self.idx_hist_item_unapply_widget_lower].setValue(plot_y)
				self.newThresholdLower()
	
	def histparammouseup(self,event):
		if self.curthresholdcontrol == self.idx_threshold_control_edit_only:
			return
		
		# Swap control if shift button is pressed
		is_not_reverse_control = True
		modifiers = event.modifiers()
		if modifiers & QtCore.Qt.ShiftModifier:
			is_not_reverse_control = False
		
		hist_x,hist_y=self.whistparam.scr2plot(event.x(),event.y())
		if self.curthresholdcontrol == self.idx_threshold_control_lower:
			if is_not_reverse_control:
				self.hist_map_list[self.curhistidx][self.idx_hist_item_unapply_widget_lower].setValue(hist_x)
				self.newThresholdLower()
			else:
				self.hist_map_list[self.curhistidx][self.idx_hist_item_unapply_widget_upper].setValue(hist_x)
				self.newThresholdUpper()
				
		else:
			assert self.curthresholdcontrol == self.idx_threshold_control_upper, "MRK_DEBUG"
			if is_not_reverse_control:
				self.hist_map_list[self.curhistidx][self.idx_hist_item_unapply_widget_upper].setValue(hist_x)
				self.newThresholdUpper()
			else:
				self.hist_map_list[self.curhistidx][self.idx_hist_item_unapply_widget_lower].setValue(hist_x)
				self.newThresholdLower()

if __name__ == "__main__":
	main()
