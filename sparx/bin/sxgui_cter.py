#!/usr/bin/env python
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
from OpenGL import GL,GLUT
from math import *
import os
import sys
# from numpy import array,arange
import traceback

try:
	from PyQt4 import QtCore, QtGui, QtOpenGL
	from PyQt4.QtCore import Qt
	from PyQt4.QtCore import QTimer
	from emshape import *
	from valslider import *
	from emplot2d import EMPlot2DWidget
except:
	print "Warning: PyQt4 must be installed"
	sys.exit(1)

from sparx import *
from optparse import OptionParser

def main():
	progname = os.path.basename(sys.argv[0])
	usage = progname + """  cter_ctf_file 
	This GUI application is designed for the evaluation of micrographs using the parameters outputed by CTER.
	"""
	parser = OptionParser(usage, version=SPARXVERSION)
	# No options!!! Does not need to call parser.add_option()
	
	(options, args) = parser.parse_args(sys.argv[1:])
	
	if len(args) > 2:
		print "see usage " + usage
		sys.exit()
	
	from emapplication import EMApp
	app=EMApp()
	
	cter_ctf_file = None
	if len(args) == 1 and args[0] != "":
		cter_ctf_file = args[0]
	# else: # Do nothig
	
	# Make sure main window is shown, raised, and activated upon startup.
	# gui=SXGuiCter(cter_ctf_file)
	gui=SXGuiCter()
	gui.show()
	gui.raise_()
	gui.activateWindow()
	
	# read CTER CTF file if necessary
	if cter_ctf_file != None:
		gui.readCterCtfFile(os.path.relpath(cter_ctf_file))
	
#	try:
#		# gui.wimgmicthumb.qt_parent.raise_() # wimgmicthumb should not be visible upon start-up
#		# gui.wfft.qt_parent.raise_()
#		gui.wplotrotavgcoarse.qt_parent.raise_()
#		gui.wplotrotavgfine.qt_parent.raise_()
#		gui.whistparam.qt_parent.raise_()
#		gui.wplotparam.qt_parent.raise_()
#		gui.raise_()
#		gui.activateWindow()
#	except: 
#		print "Recieved unexpected exception in main(): ", sys.exc_info()[0]
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
	
	def full_refresh(self):
		'''
		This function is called from resizeGL and from the inspector when somebody toggles the display of a line
		'''
		self.needupd=1
		self.del_shapes(("xcross","ycross","lcross","Circle")) # for EMPlot2DInspector
		self.del_shapes(("error_astig","error_def","error_ctf")) # for SXGuiCter.wplotrotavgcoarse & SXGuiCter.wplotrotavgfine
		# self.del_shapes(("error_astig","error_def")) # for SXGuiCter.wplotrotavgcoarse & SXGuiCter.wplotrotavgfine
		# self.del_shapes(("hist_param_shape_value","hist_param_shape_unapply_threshold_lower","hist_param_shape_apply_threshold_lower","hist_param_shape_unapply_threshold_upper", "hist_param_shape_apply_threshold_upper", "hist_param_shape_label")) # for SXGuiCter.whistparam
		# self.del_shapes(("hist_param_shape_label")) # for SXGuiCter.whistparam
		# self.del_shapes(("plot_param_shape_label")) # for SXGuiCter.wplotparam
	
	def mouseReleaseEvent(self, event):
		EMPlot2DWidget.mouseReleaseEvent(self,event)
		if event.button()==Qt.LeftButton:
			self.emit(QtCore.SIGNAL("mouseup"),event)

class SXGuiCter(QtGui.QWidget):
# 	def __init__(self, cter_ctf_file = None):
	def __init__(self):
		"""Implements the CTF fitting dialog using various EMImage and EMPlot2D widgets
		'data' is a list of (filename,ctf,im_1d,bg_1d,quality)
		'parms' is [box size,ctf,box coord,set of excluded boxnums,quality,oversampling]
		"""
		try:
			from emimage2d import EMImage2DWidget
		except:
			print "Cannot import EMAN image GUI objects (EMImage2DWidget)"
			sys.exit(1)
		
# 		QtGui.QWidget.__init__(self,None)
		super(SXGuiCter, self).__init__(None)
		
		# MRK_TEST: Flags to test experimental functions
		self.is_enable_max_power = False # MRK_TEST: Can this be an option in future?
		
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
		
		i_enum = -1
		i_enum += 1; self.idx_cter_id           = i_enum # <extra> entry id
		i_enum += 1; self.idx_cter_select       = i_enum # <extra> selected state
		i_enum += 1; self.idx_cter_def          = i_enum # defocus (ym)
		i_enum += 1; self.idx_cter_cs           = i_enum # Cs (mm)
		i_enum += 1; self.idx_cter_vol          = i_enum # voltage(kV)
		i_enum += 1; self.idx_cter_apix         = i_enum # pixel size (A)
		i_enum += 1; self.idx_cter_bfactor      = i_enum # B-factor (A^2)
		i_enum += 1; self.idx_cter_ac           = i_enum # amp contrast (%)
		i_enum += 1; self.idx_cter_astig_amp    = i_enum # astigmatism amplitude (um)
		i_enum += 1; self.idx_cter_astig_ang    = i_enum # astigmatism angle
		i_enum += 1; self.idx_cter_sd_def       = i_enum # std dev of defocus (um)
		i_enum += 1; self.idx_cter_sd_astig_amp = i_enum # std dev of ast amp (A)
		i_enum += 1; self.idx_cter_sd_astig_ang = i_enum # std dev of ast angle
		i_enum += 1; self.idx_cter_cv_def       = i_enum # coefficient of variation of defocus (%)
		i_enum += 1; self.idx_cter_cv_astig_amp = i_enum # coefficient of variation of ast amp (%)
		i_enum += 1; self.idx_cter_error_def    = i_enum # frequency at which signal drops by 50% due to estimated error of defocus alone (1/A)
		i_enum += 1; self.idx_cter_error_astig  = i_enum # frequency at which signal drops by 50% due to estimated error of defocus and astigmatism (1/A)
		i_enum += 1; self.idx_cter_error_ctf    = i_enum # limit frequency by CTF error 
		i_enum += 1; self.idx_cter_mic_name     = i_enum # Micrograph name
#		i_enum += 1; self.idx_cter_pwrot_name   = i_enum # <extra> CTER power spectrum rotational average file name
		if self.is_enable_max_power == True: i_enum += 1; self.idx_cter_max_power = i_enum # MRK_TEST: <extra> maximum power in experimental rotational average (with astigmatism)
		i_enum += 1; self.n_idx_cter            = i_enum
		self.n_idx_cter_extra = 2
		if self.is_enable_max_power == True: self.n_idx_cter_extra += 1 # MRK_TEST:
		
		i_enum = -1
		i_enum += 1; self.idx_cter_item_label   =  i_enum
		i_enum += 1; self.idx_cter_item_widget  =  i_enum
		i_enum += 1; self.n_idx_cter_item       =  i_enum
		
		# Mapping for parameter value items (line edit widgets)
		self.value_map_list = [None] * self.n_idx_cter
		self.value_map_list[self.idx_cter_id]           = ["CTER ID", None]
		self.value_map_list[self.idx_cter_select]       = ["Select", None]
		self.value_map_list[self.idx_cter_def]          = ["Defocus [um]", None]
		self.value_map_list[self.idx_cter_cs]           = ["Cs [mm]", None]
		self.value_map_list[self.idx_cter_vol]          = ["Voltage [kV]", None]
		self.value_map_list[self.idx_cter_apix]         = ["Pixel Size [A]", None]
		self.value_map_list[self.idx_cter_bfactor]      = ["B-factor [A^2]", None]
		self.value_map_list[self.idx_cter_ac]           = ["Amp. Contrast [%]", None]
		self.value_map_list[self.idx_cter_astig_amp]    = ["Astig. Amp.[um]", None]
		self.value_map_list[self.idx_cter_astig_ang]    = ["Astig. Ang.[deg]", None]
		self.value_map_list[self.idx_cter_sd_def]       = ["Defocus SD [um]", None]
		self.value_map_list[self.idx_cter_sd_astig_amp] = ["Astig. Amp. SD [um]", None]
		self.value_map_list[self.idx_cter_sd_astig_ang] = ["Astig. Ang. SD [deg]", None]
		self.value_map_list[self.idx_cter_cv_def]       = ["Defocus CV [%]", None]
		self.value_map_list[self.idx_cter_cv_astig_amp] = ["Astig. Amp. CV [%]", None]
		self.value_map_list[self.idx_cter_error_def]    = ["Defocus Freq. Limit [1/A]", None]
		self.value_map_list[self.idx_cter_error_astig]  = ["Astig. Freq. Limit [1/A]", None]
		self.value_map_list[self.idx_cter_error_ctf]    = ["CTF Freq. Limit [1/A]", None]
		self.value_map_list[self.idx_cter_mic_name]     = ["Micrograph", None]
#		self.value_map_list[self.idx_cter_pwrot_name]   = ["PW. Rot. File", None]
		if self.is_enable_max_power == True: self.value_map_list[self.idx_cter_max_power] = ["Max Power", None] # MRK_TEST:
		
		i_enum = -1
		i_enum += 1; self.idx_rotinf_cter_id        = i_enum # line number == cter id
		i_enum += 1; self.idx_rotinf_freq           = i_enum # spatial frequency (1/A)
		i_enum += 1; self.idx_rotinf_exp_no_astig   = i_enum # experimental rotational average (no astigmatism)
		i_enum += 1; self.idx_rotinf_fit_no_astig   = i_enum # fitted rotational average (no astigmatism)
		i_enum += 1; self.idx_rotinf_exp_with_astig = i_enum # experimental rotational average (with astigmatism)
		i_enum += 1; self.idx_rotinf_fit_with_astig = i_enum # fitted rotational average (with astigmatism)
		i_enum += 1; self.n_idx_rotinf              = i_enum
		
		i_enum = -1
		i_enum += 1; self.idx_sort_id           = i_enum
		i_enum += 1; self.idx_sort_def          = i_enum
		i_enum += 1; self.idx_sort_astig_amp    = i_enum
		i_enum += 1; self.idx_sort_astig_ang    = i_enum
		i_enum += 1; self.idx_sort_sd_def       = i_enum
		i_enum += 1; self.idx_sort_sd_astig_amp = i_enum
		i_enum += 1; self.idx_sort_sd_astig_ang = i_enum
		i_enum += 1; self.idx_sort_cv_def       = i_enum
		i_enum += 1; self.idx_sort_cv_astig_amp = i_enum
		i_enum += 1; self.idx_sort_error_def    = i_enum
		i_enum += 1; self.idx_sort_error_astig  = i_enum
		i_enum += 1; self.idx_sort_error_ctf    = i_enum
		if self.is_enable_max_power == True: i_enum += 1; self.idx_sort_max_power = i_enum # MRK_TEST:
		i_enum += 1; self.n_idx_sort            = i_enum
		
		i_enum = -1
		i_enum += 1; self.idx_sort_item_idx_cter =  i_enum
		i_enum += 1; self.n_idx_sort_item        =  i_enum
		
		# Mapping for sorting item (combo box widget)
		# Includes mapping from idx_sort to idx_cter
		self.sort_map_list = [None] * self.n_idx_sort
		self.sort_map_list[self.idx_sort_id]           = [self.idx_cter_id]
		self.sort_map_list[self.idx_sort_def]          = [self.idx_cter_def]
		self.sort_map_list[self.idx_sort_astig_amp]    = [self.idx_cter_astig_amp]
		self.sort_map_list[self.idx_sort_astig_ang]    = [self.idx_cter_astig_ang]
		self.sort_map_list[self.idx_sort_sd_def]       = [self.idx_cter_sd_def]
		self.sort_map_list[self.idx_sort_sd_astig_amp] = [self.idx_cter_sd_astig_amp]
		self.sort_map_list[self.idx_sort_sd_astig_ang] = [self.idx_cter_sd_astig_ang]
		self.sort_map_list[self.idx_sort_cv_def]       = [self.idx_cter_cv_def]
		self.sort_map_list[self.idx_sort_cv_astig_amp] = [self.idx_cter_cv_astig_amp]
		self.sort_map_list[self.idx_sort_error_def]    = [self.idx_cter_error_def]
		self.sort_map_list[self.idx_sort_error_astig]  = [self.idx_cter_error_astig]
		self.sort_map_list[self.idx_sort_error_ctf]    = [self.idx_cter_error_ctf]
		if self.is_enable_max_power == True: self.sort_map_list[self.idx_sort_max_power] = [self.idx_cter_max_power] # MRK_TEST:
		
		i_enum = -1
		i_enum += 1; self.idx_hist_def          = i_enum
		i_enum += 1; self.idx_hist_astig_amp    = i_enum
		i_enum += 1; self.idx_hist_astig_ang    = i_enum
		i_enum += 1; self.idx_hist_sd_def       = i_enum
		i_enum += 1; self.idx_hist_sd_astig_amp = i_enum
		i_enum += 1; self.idx_hist_sd_astig_ang = i_enum
		i_enum += 1; self.idx_hist_cv_def       = i_enum
		i_enum += 1; self.idx_hist_cv_astig_amp = i_enum
		i_enum += 1; self.idx_hist_error_def    = i_enum
		i_enum += 1; self.idx_hist_error_astig  = i_enum
		i_enum += 1; self.idx_hist_error_ctf    = i_enum
		if self.is_enable_max_power == True: i_enum += 1; self.idx_hist_max_power = i_enum  # MRK_TEST:
		i_enum += 1; self.n_idx_hist            = i_enum
		
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
		
		# Mapping for histogram items (combo box widget) and threshold setting (line edit widgets)
		# Includes mapping from idx_hist to idx_cter and idx_sort
		self.hist_map_list = [None] * self.n_idx_hist
		self.hist_map_list[self.idx_hist_def]          = [self.idx_cter_def, self.idx_sort_def, 0, 5, 0, 5, None, None, 0, 5, None, None]
		self.hist_map_list[self.idx_hist_astig_amp]    = [self.idx_cter_astig_amp, self.idx_sort_astig_amp, 0, 1, 0, 1, None, None, 0, 1, None, None]
		self.hist_map_list[self.idx_hist_astig_ang]    = [self.idx_cter_astig_ang, self.idx_sort_astig_ang, 0, 180, 0, 180, None, None, 0, 180, None, None]
		self.hist_map_list[self.idx_hist_sd_def]       = [self.idx_cter_sd_def, self.idx_sort_sd_def, 0, 5, 0, 5, None, None, 0, 5, None, None]
		self.hist_map_list[self.idx_hist_sd_astig_amp] = [self.idx_cter_sd_astig_amp, self.idx_sort_sd_astig_amp, 0, 1, 0, 1, None, None, 0, 1, None, None]
		self.hist_map_list[self.idx_hist_sd_astig_ang] = [self.idx_cter_sd_astig_ang, self.idx_sort_sd_astig_ang, 0, 180, 0, 180, None, None, 0, 180, None, None]
		self.hist_map_list[self.idx_hist_cv_def]       = [self.idx_cter_cv_def, self.idx_sort_cv_def, 0, 5, 0, 5, None, None, 0, 5, None, None]
		self.hist_map_list[self.idx_hist_cv_astig_amp] = [self.idx_cter_cv_astig_amp, self.idx_sort_cv_astig_amp, 0, 1, 0, 1, None, None, 0, 1, None, None]
		self.hist_map_list[self.idx_hist_error_def]    = [self.idx_cter_error_def, self.idx_sort_error_def, 0, 10, 0, 10, None, None, 0, 10, None, None]
		self.hist_map_list[self.idx_hist_error_astig]  = [self.idx_cter_error_astig, self.idx_sort_error_astig, 0, 10, 0, 10, None, None, 0, 10, None, None]
		self.hist_map_list[self.idx_hist_error_ctf]    = [self.idx_cter_error_ctf, self.idx_sort_error_ctf, 0, 10, 0, 10, None, None, 0, 10, None, None]
		if self.is_enable_max_power == True: self.hist_map_list[self.idx_hist_max_power] = [self.idx_cter_max_power, self.idx_sort_max_power, 0, 99999, 0, 99999, None, None, 0, 99999, None, None] # MRK_TEST:
		
		i_enum = -1
		i_enum += 1; self.idx_threshold_control_lower     = i_enum
		i_enum += 1; self.idx_threshold_control_upper     = i_enum
		i_enum += 1; self.idx_threshold_control_edit_only = i_enum
		i_enum += 1; self.n_idx_threshold_control         = i_enum
		
		i_enum = -1
		i_enum += 1; self.idx_threshold_control_item_label = i_enum
		i_enum += 1; self.idx_threshold_control_item_color = i_enum
		i_enum += 1; self.n_idx_threshold_control_item     = i_enum
		
		# Mapping for threshold control (combo box widget)
		self.threshold_control_map_list = [None] * self.n_idx_threshold_control
		self.threshold_control_map_list[self.idx_threshold_control_lower]     = ["Lower (blue)", "blue"]
		self.threshold_control_map_list[self.idx_threshold_control_upper]     = ["Upper (red)", "red"]
		self.threshold_control_map_list[self.idx_threshold_control_edit_only] = ["Edit Only", "black"]
		
		i_enum = -1
		i_enum += 1; self.idx_graph_exp_no_astig   = i_enum
		i_enum += 1; self.idx_graph_fit_no_astig   = i_enum
		i_enum += 1; self.idx_graph_exp_with_astig = i_enum
		i_enum += 1; self.idx_graph_fit_with_astig = i_enum
		i_enum += 1; self.n_idx_graph              = i_enum
		
		i_enum = -1
		i_enum += 1; self.idx_graph_item_name   = i_enum
		i_enum += 1; self.idx_graph_item_label  = i_enum
		i_enum += 1; self.idx_graph_idx_rotinf  = i_enum
		i_enum += 1; self.idx_graph_item_widget = i_enum
		i_enum += 1; self.n_idx_graph_item      = i_enum
		
		# Mapping for graph display setting (check box widgets)
		self.graph_map_list = [None] * self.n_idx_graph
		self.graph_map_list[self.idx_graph_exp_no_astig]   = ["exp_no_astig", "Exp. No Astig (Black)", self.idx_rotinf_exp_no_astig, None]
		self.graph_map_list[self.idx_graph_fit_no_astig]   = ["fit_no_astig", "Fit. No Astig (Blue)", self.idx_rotinf_fit_no_astig, None]
		self.graph_map_list[self.idx_graph_exp_with_astig] = ["exp_with_astig", "Exp. with Astig (Red)", self.idx_rotinf_exp_with_astig, None]
		self.graph_map_list[self.idx_graph_fit_with_astig] = ["fit_with_astig", "Fit. with Astig (Greed)", self.idx_rotinf_fit_with_astig, None]
		
		i_enum = -1
		i_enum += 1; self.idx_thresholdset_unapplied = i_enum
		i_enum += 1; self.idx_thresholdset_applied   = i_enum
		i_enum += 1; self.n_idx_thresholdset         = i_enum
		
		i_enum = -1
		i_enum += 1; self.idx_thresholdset_item_label  = i_enum
		i_enum += 1; self.n_idx_thresholdset_item      = i_enum
		
		# Mapping for threshold set (combo box widget)
		self.thresholdset_map_list = [None] * self.n_idx_thresholdset
		self.thresholdset_map_list[self.idx_thresholdset_unapplied] = ["Unapplied"]
		self.thresholdset_map_list[self.idx_thresholdset_applied]   = ["Applied"]
		
		self.cter_partres_file_path  = None
		self.cter_entry_list         = None
		self.cter_mic_file_path      = None
		self.cter_micthumb_file_path = None
		self.cter_pwrot_file_path    = None
		
		self.curimgmicthumbdisplay=True
		self.curplotfixscale=5
		self.curentry=None
		self.cursort=0
		self.cursortoder=False
		self.cursortselect = False
		self.curhist=0
		self.curhistdisable=False
		self.curthresholdcontrol = 0
		self.curentryperbin=10
		self.cursyncsort = False
		self.curthresholdset=0
		
		self.child_status_list = [] 
		
		#
		# NOTE: 2016/03/09 Toshio Moriya
		# To set window flags of EMGLWidget (SXPlot2DWidget and EMImage2DWidget) window,
		# we have to go through its qt_parent attribute to call setWindowTitle()...
		# 
#		self.wfft=EMImage2DWidget()
#		self.wfft.setWindowTitle("sxgui_cter - 2D FFT")
#		self.wfft.mmode="app"
#		self.wfft.qt_parent.setWindowFlags((self.qt_parent.wfft.windowFlags()| Qt.CustomizeWindowHint) & ~Qt.WindowMinimizeButtonHint) # Disabled minimize icon button in window title bar
#		self.is_wfft_minimized = False
		
		self.wimgmicthumb=EMImage2DWidget()
		self.wimgmicthumb.setWindowTitle("sxgui_cter - Micrograph Thumbnail")
		self.wimgmicthumb.mmode="app"
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
		
		self.wplotparam=SXPlot2DWidget()
		self.wplotparam.setWindowTitle("sxgui_cter - Sort Plot")
		self.wplotparam.qt_parent.setWindowFlags((self.wplotparam.qt_parent.windowFlags()| Qt.CustomizeWindowHint) & ~Qt.WindowMinimizeButtonHint) # Disabled minimize icon button in window title bar
		self.is_wplotparam_minimized = False
		
#		self.wfft.connect(self.wfft,QtCore.SIGNAL("mousedown"),self.fftmousedown)
#		self.wfft.connect(self.wfft,QtCore.SIGNAL("mousedrag"),self.fftmousedrag)
#		self.wfft.connect(self.wfft,QtCore.SIGNAL("mouseup")  ,self.fftmouseup)
#		self.wimgmicthumb.connect(self.wimgmicthumb,QtCore.SIGNAL("mousedown"),self.imgmicthumbmousedown)
#		self.wimgmicthumb.connect(self.wimgmicthumb,QtCore.SIGNAL("mousedrag"),self.imgmicthumbmousedrag)
#		self.wimgmicthumb.connect(self.wimgmicthumb,QtCore.SIGNAL("mouseup")  ,self.imgmicthumbmouseup)
#		self.wplotrotavgcoarse.connect(self.wplotrotavgcoarse,QtCore.SIGNAL("mousedown"),self.plotmousedown)
#		self.wplotrotavgfine.connect(self.wplotrotavgfine,QtCore.SIGNAL("mousedown"),self.plotmousedown)
		self.whistparam.connect(self.whistparam,QtCore.SIGNAL("mouseup"),self.histparammouseup)
		self.wplotparam.connect(self.wplotparam,QtCore.SIGNAL("mouseup"),self.plotparammouseup)
		
		# This object is itself a widget we need to set up
		self.gbl = QtGui.QGridLayout(self)
		self.gbl.setMargin(8)
		self.gbl.setSpacing(6)
		
		# --------------------------------------------------------------------------------
		# Columns 1-3
		# --------------------------------------------------------------------------------
		grid_col = 0;
		col_span_label = 1
		col_span_edit = 2
		col_span_sublabel = 2
		col_span_subedit = 1
		assert(col_span_label + col_span_edit == col_span_sublabel + col_span_subedit) # MRK_DEBUG
		col_span = col_span_label + col_span_edit
		grid_row = 0
		
		labelwidth=70
		sublabelwidth=140
		
		self.pbopencter=QtGui.QPushButton("Open CTER CTF file")
		self.gbl.addWidget(self.pbopencter,grid_row,grid_col,1,col_span)
		grid_row += 1
		
		# Make space
		grid_row+=1
		
		temp_label=QtGui.QLabel("<b>Selection Summary:</b>",self)
		temp_label.setAlignment(QtCore.Qt.AlignHCenter|QtCore.Qt.AlignVCenter)
		self.gbl.addWidget(temp_label,grid_row,grid_col,1,col_span)
		grid_row += 1
		
		temp_label=QtGui.QLabel("Num. of entries",self)
		temp_label.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignVCenter)
		temp_label.setMinimumSize(QtCore.QSize(labelwidth,20))
		self.gbl.addWidget(temp_label,grid_row,grid_col,1,col_span_label)
		self.vbnentry=ValBox(self,(0,10000),None,0)
		self.vbnentry.setEnabled(False)
		self.vbnentry.intonly=True
		self.gbl.addWidget(self.vbnentry,grid_row,grid_col+col_span_label,1,col_span_edit)
		grid_row+=1
		
		temp_label=QtGui.QLabel("Unchecked",self)
		temp_label.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignVCenter)
		temp_label.setMinimumSize(QtCore.QSize(labelwidth,20))
		self.gbl.addWidget(temp_label,grid_row,grid_col,1,col_span_label)
		self.vbuncheckcounts=ValBox(self,(0,1000000),None,0)
		self.vbuncheckcounts.setEnabled(False)
		self.vbuncheckcounts.intonly=True
		self.gbl.addWidget(self.vbuncheckcounts,grid_row,grid_col+col_span_label,1,col_span_edit)
		grid_row+=1
		
		temp_label=QtGui.QLabel("Ratio",self)
		temp_label.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignVCenter)
		temp_label.setMinimumSize(QtCore.QSize(labelwidth,20))
		self.gbl.addWidget(temp_label,grid_row,grid_col,1,col_span_label)
		self.vbuncheckratio=ValBox(self,(0,1.0),None,0)
		self.vbuncheckratio.setEnabled(False)
		self.gbl.addWidget(self.vbuncheckratio,grid_row,grid_col+col_span_label,1,col_span_edit)
		grid_row+=1
		
		# Make space
		grid_row+=1
		
		temp_label=QtGui.QLabel("<b>Electron Microscopy:</b>",self)
		temp_label.setAlignment(QtCore.Qt.AlignHCenter|QtCore.Qt.AlignVCenter)
		self.gbl.addWidget(temp_label,grid_row,grid_col,1,col_span)
		grid_row += 1
		
		self.add_value_widget(self.idx_cter_vol,0,500,grid_row,grid_col,col_span_label,col_span_edit,labelwidth=labelwidth)
		grid_row+=1
		
		self.add_value_widget(self.idx_cter_cs,0,5,grid_row,grid_col,col_span_label,col_span_edit,labelwidth=labelwidth)
		grid_row+=1
		
		self.add_value_widget(self.idx_cter_ac,0,100,grid_row,grid_col,col_span_label,col_span_edit,labelwidth=labelwidth)
		grid_row+=1
		
		self.add_value_widget(self.idx_cter_apix,0,500,grid_row,grid_col,col_span_label,col_span_edit,labelwidth=labelwidth)
		grid_row+=1
		
		# Make space
		grid_row+=1
		
		temp_label=QtGui.QLabel("<b>Display Micrograph:</b>",self)
		temp_label.setAlignment(QtCore.Qt.AlignHCenter|QtCore.Qt.AlignVCenter)
		self.gbl.addWidget(temp_label,grid_row,grid_col,1,col_span)
		grid_row += 1
		
		temp_label=QtGui.QLabel("Open Window",self)
		temp_label.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignVCenter)
		temp_label.setMinimumSize(QtCore.QSize(sublabelwidth,20))
		self.gbl.addWidget(temp_label,grid_row,grid_col, 1, col_span_sublabel)
		self.cbmicthumbdisplay=CheckBox(None,None,self.curimgmicthumbdisplay)
#		self.cbmicthumbdisplay.setEnabled(False)
		self.gbl.addWidget(self.cbmicthumbdisplay,grid_row,grid_col+col_span_sublabel,1,col_span_subedit)
		grid_row+=1
		
		# Make space
		grid_row+=1
		
		temp_label=QtGui.QLabel("<b>Display Curves:</b>",self)
		temp_label.setAlignment(QtCore.Qt.AlignHCenter|QtCore.Qt.AlignVCenter)
		self.gbl.addWidget(temp_label,grid_row,grid_col,1,col_span)
		grid_row += 1
		
		for idx_graph in xrange(self.n_idx_graph):
			temp_label=QtGui.QLabel(self.graph_map_list[idx_graph][self.idx_graph_item_label],self)
			temp_label.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignVCenter)
			temp_label.setMinimumSize(QtCore.QSize(sublabelwidth,20))
			self.gbl.addWidget(temp_label,grid_row,grid_col,1,col_span_sublabel)
			self.graph_map_list[idx_graph][self.idx_graph_item_widget]=CheckBox(None,None,True)
			self.gbl.addWidget(self.graph_map_list[idx_graph][self.idx_graph_item_widget],grid_row,grid_col+col_span_sublabel,1,col_span_subedit)
			grid_row += 1
		
		temp_label=QtGui.QLabel("Plot Fix Scale",self)
		temp_label.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignVCenter)
		temp_label.setMinimumSize(QtCore.QSize(sublabelwidth,20))
		self.gbl.addWidget(temp_label,grid_row,grid_col, 1, col_span_sublabel)
		self.vbplotfixscale=ValBox(self,(0,99999),None,self.curplotfixscale)
		self.gbl.addWidget(self.vbplotfixscale,grid_row,grid_col+col_span_sublabel,1,col_span_subedit)
		grid_row+=1
		
		self.pbrefreshgraphs=QtGui.QPushButton("Refresh Graphs")
		self.pbrefreshgraphs.setEnabled(False)
		self.gbl.addWidget(self.pbrefreshgraphs,grid_row,grid_col,1,col_span)
		grid_row += 1
		
		# --------------------------------------------------------------------------------
		# Columns 4
		# --------------------------------------------------------------------------------
		grid_col += col_span
		col_span = 1
		grid_row = 0
		
		# plot list and plot mode combobox
		row_span_entry_list = 24
		self.lbentry=SXListWidget(self) # self.lbentry=e2ctf.MyListWidget(self)
		self.lbentry.setSizePolicy(QtGui.QSizePolicy.Preferred,QtGui.QSizePolicy.Expanding)
		self.lbentry.setMinimumWidth(220)
		self.gbl.addWidget(self.lbentry,grid_row,grid_col,row_span_entry_list,col_span)
		grid_row += row_span_entry_list
		
		grid_row_entry_list = grid_row - 1
		grid_col_entry_list = grid_col
		
		# --------------------------------------------------------------------------------
		# Columns 5-7 (for Micrograph/CTER Entry), 7-8 (for Unapplied Threshold), 9-10 (for Applied Threshold)
		# --------------------------------------------------------------------------------
		grid_col += col_span
		col_span_1st_label = 2
		col_span_1st_edit = 1
		col_span_1st_subspace = 1
		col_span_1st_sublabel = 1
		col_span_1st_subedit = 1
		col_span_2nd_sublabel = col_span_1st_sublabel
		col_span_2nd_subedit = col_span_1st_subedit
		col_span_3rd_sublabel = col_span_1st_sublabel
		col_span_3rd_subedit = col_span_1st_subedit
		assert(col_span_1st_label + col_span_1st_edit == col_span_1st_subspace + col_span_1st_sublabel + col_span_1st_subedit)
		col_span_1st = col_span_1st_label + col_span_1st_edit
		col_span_1st_sub = col_span_1st_sublabel + col_span_1st_subedit
		col_span_2nd = 2
		col_span_3rd = 2
		col_span = col_span_1st + col_span_2nd + col_span_3rd
		grid_row = 0
		
		grid_col_1st = grid_col
		grid_col_1st_sub = grid_col + col_span_1st_subspace
		grid_col_2nd = grid_col_1st + col_span_1st
		grid_col_3rd = grid_col_2nd + col_span_2nd
		
		labelwidth=160
		editwidth=100
		sublabelwidth=editwidth
		
		temp_label=QtGui.QLabel("<b>Current Entry Info:</b>",self)
		temp_label.setAlignment(QtCore.Qt.AlignHCenter|QtCore.Qt.AlignVCenter)
		self.gbl.addWidget(temp_label,grid_row,grid_col_1st_sub,1,col_span_1st_sub)
		temp_label=QtGui.QLabel("<b>Unapplied Thresholds:</b>",self)
		temp_label.setAlignment(QtCore.Qt.AlignHCenter|QtCore.Qt.AlignVCenter)
		self.gbl.addWidget(temp_label,grid_row,grid_col_2nd,1,col_span_2nd)
		temp_label=QtGui.QLabel("<b>Applied Thresholds:</b>",self)
		temp_label.setAlignment(QtCore.Qt.AlignHCenter|QtCore.Qt.AlignVCenter)
		self.gbl.addWidget(temp_label,grid_row,grid_col_3rd,1,col_span_3rd)
		grid_row += 1
		
		temp_label=QtGui.QLabel("Sorted ID",self)
		temp_label.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignVCenter)
		temp_label.setMinimumSize(QtCore.QSize(labelwidth,20))
		self.gbl.addWidget(temp_label,grid_row,grid_col_1st,1,col_span_1st_label)
		self.ssortedid=ValBox(self,(0,10000),None,0)
		self.ssortedid.setEnabled(False)
		self.ssortedid.intonly=True
		self.gbl.addWidget(self.ssortedid,grid_row,grid_col_1st+col_span_1st_label,1,col_span_1st_edit)
		grid_row+=1
		
		self.add_value_widget(self.idx_cter_id,0,10000,grid_row,grid_col_1st,col_span_1st_label,col_span_1st_edit,intonly=True,labelwidth=labelwidth,editwidth=editwidth)
		grid_row+=1
		
		self.add_value_widget(self.idx_cter_select,0,1,grid_row,grid_col_1st,col_span_1st_label,col_span_1st_edit,intonly=True,labelwidth=labelwidth,editwidth=editwidth)
		grid_row+=1
		
		for idx_hist in xrange(self.n_idx_hist):
			self.add_value_widget_with_threshold(idx_hist,grid_row,grid_col_1st,grid_col_2nd,grid_col_3rd,col_span_1st_label,col_span_1st_edit,labelwidth=labelwidth,editwidth=editwidth)
			grid_row+=1
		self.hist_map_list[0][self.idx_hist_item_unapply_widget_lower].setEnabled(True)
		# self.hist_map_list[0][self.idx_hist_item_unapply_widget_upper].setEnabled(True)
		
		self.add_value_widget(self.idx_cter_bfactor,0,1600,grid_row,grid_col_1st,col_span_1st_label,col_span_1st_edit,labelwidth=labelwidth,editwidth=editwidth)
		grid_row+=1
		
		# make space
		grid_row += 1
		
		temp_label=QtGui.QLabel("<b>Sort CTER Entries:</b>",self)
		temp_label.setAlignment(QtCore.Qt.AlignHCenter|QtCore.Qt.AlignVCenter)
		self.gbl.addWidget(temp_label,grid_row,grid_col_1st_sub,1,col_span_1st_sub)
		
		temp_label=QtGui.QLabel("<b>Histogram & Plot Settings:</b>",self)
		temp_label.setAlignment(QtCore.Qt.AlignHCenter|QtCore.Qt.AlignVCenter)
		self.gbl.addWidget(temp_label,grid_row,grid_col_2nd,1,col_span_2nd)
		
		temp_label=QtGui.QLabel("<b>Save/Load Thresholds:</b>",self)
		temp_label.setAlignment(QtCore.Qt.AlignHCenter|QtCore.Qt.AlignVCenter)
		self.gbl.addWidget(temp_label,grid_row,grid_col_3rd,1,col_span_3rd)
		grid_row += 1
		
		self.ssort=QtGui.QComboBox(self)
		for map_entry in self.sort_map_list:
			idx_cter = map_entry[self.idx_sort_item_idx_cter]
			self.ssort.addItem(self.value_map_list[idx_cter][self.idx_cter_item_label])
		self.ssort.setCurrentIndex(self.cursort)
		self.gbl.addWidget(self.ssort,grid_row,grid_col_1st_sub,1,col_span_1st_sub)
		
		self.shist=QtGui.QComboBox(self)
		for map_entry in self.hist_map_list:
			idx_cter = map_entry[self.idx_hist_item_idx_cter]
			self.shist.addItem(self.value_map_list[idx_cter][self.idx_cter_item_label])
		self.shist.setCurrentIndex(self.curhist)
		self.gbl.addWidget(self.shist,grid_row,grid_col_2nd,1,col_span_2nd)
		
		self.sthresholdset=QtGui.QComboBox(self)
		for map_entry in self.thresholdset_map_list:
			self.sthresholdset.addItem(map_entry[self.idx_thresholdset_item_label])
		self.sthresholdset.setCurrentIndex(self.curthresholdset)
		self.gbl.addWidget(self.sthresholdset,grid_row,grid_col_3rd,1,col_span_3rd)
		grid_row += 1
		
		temp_label=QtGui.QLabel("Decending",self)
		temp_label.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignVCenter)
		temp_label.setMinimumSize(QtCore.QSize(sublabelwidth,20))
		self.gbl.addWidget(temp_label,grid_row,grid_col_1st_sub,1,col_span_1st_sublabel)
		self.cbsortoder=CheckBox(None,None,self.cursortoder)
		self.gbl.addWidget(self.cbsortoder,grid_row,grid_col_1st_sub+col_span_1st_sublabel,1,col_span_1st_subedit)
		
		temp_label=QtGui.QLabel("Move Threshold",self)
		temp_label.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignVCenter)
		temp_label.setMinimumSize(QtCore.QSize(sublabelwidth,20))
		self.gbl.addWidget(temp_label,grid_row,grid_col_2nd,1,col_span_2nd_sublabel)
		self.sthresholdcontrol=QtGui.QComboBox(self)
		# self.sthresholdcontrol.setStyleSheet("color: rgb(255,0,0);") # NOTE: Toshio Moriya 2016/01/22 Unfortunately, this will over write the individual item color...
		for idx_threshold_control in xrange(self.n_idx_threshold_control):
			map_entry = self.threshold_control_map_list[idx_threshold_control]
			self.sthresholdcontrol.addItem(map_entry[self.idx_threshold_control_item_label])
			self.sthresholdcontrol.setItemData(idx_threshold_control, QtGui.QColor(map_entry[self.idx_threshold_control_item_color]), Qt.TextColorRole);
		self.sthresholdcontrol.setCurrentIndex(self.curthresholdcontrol)
		self.gbl.addWidget(self.sthresholdcontrol,grid_row,grid_col_2nd+col_span_2nd_sublabel,1,col_span_2nd_subedit)
		
		self.pbsavethresholdset=QtGui.QPushButton("Save")
		self.pbsavethresholdset.setEnabled(False)
		self.gbl.addWidget(self.pbsavethresholdset,grid_row,grid_col_3rd,1,col_span_3rd_sublabel)
		self.pbloadthresholdset=QtGui.QPushButton("Load")
		self.pbloadthresholdset.setEnabled(False)
		self.gbl.addWidget(self.pbloadthresholdset,grid_row,grid_col_3rd+col_span_3rd_sublabel,1,col_span_3rd_subedit)
		grid_row += 1
		
		temp_label=QtGui.QLabel("Sort Select",self)
		temp_label.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignVCenter)
		temp_label.setMinimumSize(QtCore.QSize(sublabelwidth,20))
		self.gbl.addWidget(temp_label,grid_row,grid_col_1st_sub,1,col_span_1st_sublabel)
		self.cbsortselect=CheckBox(None,None,self.cursortselect)
		self.gbl.addWidget(self.cbsortselect,grid_row,grid_col_1st_sub+col_span_1st_sublabel,1,col_span_1st_subedit)
		
		temp_label=QtGui.QLabel("Sync. Sort",self)
		temp_label.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignVCenter)
		temp_label.setMinimumSize(QtCore.QSize(sublabelwidth,20))
		self.gbl.addWidget(temp_label,grid_row,grid_col_2nd,1,col_span_2nd_sublabel)
		self.cbsyncsort=CheckBox(None,None,self.cursyncsort)
		self.gbl.addWidget(self.cbsyncsort,grid_row,grid_col_2nd+col_span_2nd_sublabel,1,col_span_2nd_subedit)
		
		temp_label=QtGui.QLabel("<b>Save Selection:</b>",self)
		temp_label.setAlignment(QtCore.Qt.AlignHCenter|QtCore.Qt.AlignVCenter)
		self.gbl.addWidget(temp_label,grid_row,grid_col_3rd,1,col_span_3rd)
		grid_row += 1
		
		temp_label=QtGui.QLabel("counts/bin",self)
		temp_label.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignVCenter)
		temp_label.setMinimumSize(QtCore.QSize(sublabelwidth,20))
		self.gbl.addWidget(temp_label,grid_row,grid_col_2nd, 1, col_span_2nd_sublabel)
		self.vsentryperbin=ValBox(self,(0,10000),None,self.curentryperbin)
		self.vsentryperbin.setIntonly(True)
		self.gbl.addWidget(self.vsentryperbin,grid_row,grid_col_2nd+col_span_2nd_sublabel,1,col_span_2nd_subedit)
		
		temp_label=QtGui.QLabel("File Suffix",self)
		temp_label.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignVCenter)
		temp_label.setMinimumSize(QtCore.QSize(sublabelwidth,20))
		self.gbl.addWidget(temp_label,grid_row,grid_col_3rd,1,col_span_3rd_sublabel)
		self.vfilesuffix=StringBox(self,None,"Trial00")
		self.gbl.addWidget(self.vfilesuffix,grid_row,grid_col_3rd+col_span_3rd_sublabel,1,col_span_3rd_subedit)
		grid_row += 1
		
		self.pbreaplysort=QtGui.QPushButton("Reapply Sort")
		self.pbreaplysort.setEnabled(False)
		self.gbl.addWidget(self.pbreaplysort,grid_row,grid_col_1st_sub,1,col_span_1st_sub)
		
		self.pbapplyallthreshold=QtGui.QPushButton("Apply All Thresholds")
		self.pbapplyallthreshold.setEnabled(False)
		self.gbl.addWidget(self.pbapplyallthreshold,grid_row,grid_col_2nd,1,col_span_2nd)
		
		self.pbsaveselection=QtGui.QPushButton("Save Selection")
		self.pbsaveselection.setEnabled(False)
		self.gbl.addWidget(self.pbsaveselection,grid_row,grid_col_3rd,1,col_span_3rd)
		grid_row += 1
		
		# --------------------------------------------------------------------------------
		# Set spacer for global grid layout
		# --------------------------------------------------------------------------------
		self.gbl.setRowStretch(grid_row_entry_list, self.gbl.rowStretch(grid_row_entry_list) + 1) # Give a priority to the last empty row of the entry box list for stretching
		self.gbl.setColumnStretch(grid_col_entry_list, self.gbl.columnStretch(grid_col_entry_list) + 1) # Give a priority to the entry list box column for stretching
		
		# --------------------------------------------------------------------------------
		# Set signal handler
		# --------------------------------------------------------------------------------
		QtCore.QObject.connect(self.pbopencter, QtCore.SIGNAL("clicked(bool)"),self.openCter)
		
		QtCore.QObject.connect(self.cbmicthumbdisplay, QtCore.SIGNAL("valueChanged"),self.newMicThumbDisplay)
		
		for idx_graph in xrange(self.n_idx_graph):
			QtCore.QObject.connect(self.graph_map_list[idx_graph][self.idx_graph_item_widget], QtCore.SIGNAL("valueChanged"),self.updatePlotVisibility)
		QtCore.QObject.connect(self.vbplotfixscale, QtCore.SIGNAL("valueChanged"),self.newPlotFixScale)
		QtCore.QObject.connect(self.pbrefreshgraphs, QtCore.SIGNAL("clicked(bool)"),self.refreshGraphs)
		
		QtCore.QObject.connect(self.lbentry,QtCore.SIGNAL("currentRowChanged(int)"),self.newEntry)
#		QtCore.QObject.connect(self.lbentry,QtCore.SIGNAL("keypress"),self.entryKey)
		QtCore.QObject.connect(self.lbentry,QtCore.SIGNAL("itemChanged(QListWidgetItem*)"),self.updateEntrySelect)
		
		QtCore.QObject.connect(self.ssort,QtCore.SIGNAL("currentIndexChanged(int)"),self.newSort)
		QtCore.QObject.connect(self.cbsortoder, QtCore.SIGNAL("valueChanged"),self.newSortOrder)
		QtCore.QObject.connect(self.cbsortselect, QtCore.SIGNAL("valueChanged"),self.newSortSelect)
		QtCore.QObject.connect(self.pbreaplysort, QtCore.SIGNAL("clicked(bool)"),self.reapplySort)
		
		for idx_hist in xrange(self.n_idx_hist):
			QtCore.QObject.connect(self.hist_map_list[idx_hist][self.idx_hist_item_unapply_widget_lower],QtCore.SIGNAL("valueChanged"),self.newThresholdLower)
			QtCore.QObject.connect(self.hist_map_list[idx_hist][self.idx_hist_item_unapply_widget_upper],QtCore.SIGNAL("valueChanged"),self.newThresholdUpper)
			# QtCore.QObject.connect(self.hist_map_list[idx_hist][self.idx_hist_item_unapply_widget_lower],QtCore.SIGNAL("valueChanged"),self.updateHist)
			# QtCore.QObject.connect(self.hist_map_list[idx_hist][self.idx_hist_item_unapply_widget_lower],QtCore.SIGNAL("valueChanged"),self.updatePlotParam)
			# QtCore.QObject.connect(self.hist_map_list[idx_hist][self.idx_hist_item_unapply_widget_upper],QtCore.SIGNAL("valueChanged"),self.updateHist)
			# QtCore.QObject.connect(self.hist_map_list[idx_hist][self.idx_hist_item_unapply_widget_upper],QtCore.SIGNAL("valueChanged"),self.updatePlotParam)
		
		QtCore.QObject.connect(self.shist,QtCore.SIGNAL("currentIndexChanged(int)"),self.newHist)
		QtCore.QObject.connect(self.sthresholdcontrol,QtCore.SIGNAL("currentIndexChanged(int)"),self.newThresholdControl)
		QtCore.QObject.connect(self.cbsyncsort, QtCore.SIGNAL("valueChanged"),self.newSyncSort)
		QtCore.QObject.connect(self.vsentryperbin, QtCore.SIGNAL("valueChanged"),self.newEntryPerBin)
		QtCore.QObject.connect(self.pbapplyallthreshold, QtCore.SIGNAL("clicked(bool)"),self.applyAllThresholds)
		
		QtCore.QObject.connect(self.sthresholdset,QtCore.SIGNAL("currentIndexChanged(int)"),self.newThresholdSet)
		QtCore.QObject.connect(self.pbsavethresholdset, QtCore.SIGNAL("clicked(bool)"),self.saveThresholdSet)
		QtCore.QObject.connect(self.pbloadthresholdset, QtCore.SIGNAL("clicked(bool)"),self.loadThresholdSet)
		
		QtCore.QObject.connect(self.pbsaveselection, QtCore.SIGNAL("clicked(bool)"),self.saveSelection)
		
		self.setWindowTitle("sxgui_cter - Control Panel")
		
		# Set default sizes & positions of windows in case this is the first time to run in this project directory
		# figured these values out by printing the width and height in resize event
		win_height = 512  # Let use the same height for all windows
		win_height_margin = 46
		main_win_width = 1200
		graph_win_width = 980
		img_win_width = win_height
		# Top Left
		win_top = 0
		win_left = 0
		win_width = graph_win_width
		self.whistparam.qt_parent.resize(win_width,win_height)
		self.whistparam.qt_parent.move(win_left,win_top)
		self.wplotparam.qt_parent.resize(win_width,win_height)
		self.wplotparam.qt_parent.move(win_left,win_top)
		# Top Right
		win_left = graph_win_width
		win_width = main_win_width;
		self.resize(win_width,win_height)
		self.move(win_left,win_top)
		# Bottom Left
		win_top = win_height + win_height_margin; 
		win_left = 0
		win_width = graph_win_width
		self.wplotrotavgcoarse.qt_parent.resize(win_width,win_height)
		self.wplotrotavgcoarse.qt_parent.move(win_left,win_top)
		self.wplotrotavgfine.qt_parent.resize(win_width,win_height)
		self.wplotrotavgfine.qt_parent.move(win_left,win_top)
		# Bottom Right
		# Set the image window
		win_left = graph_win_width
		win_width = img_win_width
		img_size = 512
		# scale_factor = float(win_width)/img_size
		self.wimgmicthumb.set_data(model_blank(img_size,img_size, bckg=1.0)) # resize does not work if no image is set
		self.wimgmicthumb.qt_parent.resize(win_width,win_height)
		self.wimgmicthumb.qt_parent.move(win_left,win_top)
		self.wimgmicthumb.scroll_to(-1 * img_size,-1 * img_size)
		# self.wimgmicthumb.set_scale(scale_factor)
		
		# Try to recover sizes & positions of windows of the previous GUI session
		E2loadappwin("sxgui_cter","main",self)
#		E2loadappwin("sxgui_cter","fft",self.wfft.qt_parent)
		E2loadappwin("sxgui_cter","imgmicthumb",self.wimgmicthumb.qt_parent)
		E2loadappwin("sxgui_cter","plotcoarse",self.wplotrotavgcoarse.qt_parent)
		E2loadappwin("sxgui_cter","plotfine",self.wplotrotavgfine.qt_parent)
		E2loadappwin("sxgui_cter","histparam",self.whistparam.qt_parent)
		E2loadappwin("sxgui_cter","plotparam",self.wplotparam.qt_parent)
		
#		if self.cter_entry_list:
# #			self.wfft.show()
#			self.whistparam.show()
#			self.wplotrotavgcoarse.show()
		
		### This section is responsible for background updates
		self.busy=False
#		self.needupdate=True
 		self.needredisp=False
#		self.procthread=None
#		self.errors=None # used to communicate errors back from the reprocessing thread
		
		self.timer=QTimer()
		QtCore.QObject.connect(self.timer, QtCore.SIGNAL("timeout()"), self.timeOut)
		self.timer.start(100)
		
#		# Finally, read CTER CTF file if necessary
#		if cter_ctf_file != None:
#			self.readCterCtfFile(os.path.relpath(cter_ctf_file))
		
	def add_value_widget(self, idx_cter, val_min, val_max, grid_row, grid_col, col_span_label, col_span_edit, intonly = False, labelwidth = 80, editwidth = 80):
		param_label = self.value_map_list[idx_cter][self.idx_cter_item_label]
		temp_label=QtGui.QLabel(param_label,self)
		temp_label.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignVCenter)
		temp_label.setMinimumSize(QtCore.QSize(labelwidth,20))
		self.gbl.addWidget(temp_label,grid_row,grid_col,1,col_span_label)
		val_default = val_min
		val_widget = ValBox(self,(val_min,val_max),None,val_default)
		val_widget.setEnabled(False)
		val_widget.intonly=intonly
		val_widget.text.setMinimumSize(QtCore.QSize(editwidth,0))
		self.gbl.addWidget(val_widget,grid_row,grid_col+col_span_label,1,col_span_edit)
		self.value_map_list[idx_cter][self.idx_cter_item_widget] = val_widget
	
	def add_value_widget_with_threshold(self, idx_hist, grid_row, grid_col_1st, grid_col_2nd, grid_col_3rd, col_span_1st_label, col_span_1st_edit, intonly = False, labelwidth = 80, editwidth = 80):
		val_min = self.hist_map_list[idx_hist][self.idx_hist_item_val_min]
		val_max = self.hist_map_list[idx_hist][self.idx_hist_item_val_max]
		# Add widget for parameter value
		self.add_value_widget(self.hist_map_list[idx_hist][self.idx_hist_item_idx_cter],val_min,val_max,grid_row,grid_col_1st,col_span_1st_label,col_span_1st_edit,intonly=intonly,labelwidth=labelwidth)
		# Add widget for unapplied thresholds
		self.hist_map_list[idx_hist][self.idx_hist_item_unapply_widget_lower]=ValBox(self,(val_min,val_max),None,val_min,labelwidth)
		self.hist_map_list[idx_hist][self.idx_hist_item_unapply_widget_lower].setEnabled(False)
		self.hist_map_list[idx_hist][self.idx_hist_item_unapply_widget_lower].text.setStyleSheet("color: rgb(0,0,255);")
		self.hist_map_list[idx_hist][self.idx_hist_item_unapply_widget_lower].text.setMinimumSize(QtCore.QSize(editwidth,0))
		self.gbl.addWidget(self.hist_map_list[idx_hist][self.idx_hist_item_unapply_widget_lower],grid_row,grid_col_2nd)
		self.hist_map_list[idx_hist][self.idx_hist_item_unapply_widget_upper]=ValBox(self,(val_min,val_max),None,val_max,labelwidth)
		self.hist_map_list[idx_hist][self.idx_hist_item_unapply_widget_upper].setEnabled(False)
		self.hist_map_list[idx_hist][self.idx_hist_item_unapply_widget_upper].text.setStyleSheet("color: rgb(255,0,0);")
		self.hist_map_list[idx_hist][self.idx_hist_item_unapply_widget_upper].text.setMinimumSize(QtCore.QSize(editwidth,0))
		self.gbl.addWidget(self.hist_map_list[idx_hist][self.idx_hist_item_unapply_widget_upper],grid_row,grid_col_2nd+1)
		# Add widget for applied thresholds
		self.hist_map_list[idx_hist][self.idx_hist_item_apply_widget_lower]=ValBox(self,(val_min,val_max),None,val_min,labelwidth)
		self.hist_map_list[idx_hist][self.idx_hist_item_apply_widget_lower].setEnabled(False)
		self.hist_map_list[idx_hist][self.idx_hist_item_apply_widget_lower].text.setStyleSheet("color: rgb(0,0,255);")
		self.hist_map_list[idx_hist][self.idx_hist_item_apply_widget_lower].text.setMinimumSize(QtCore.QSize(editwidth,0))
		self.gbl.addWidget(self.hist_map_list[idx_hist][self.idx_hist_item_apply_widget_lower],grid_row,grid_col_3rd)
		self.hist_map_list[idx_hist][self.idx_hist_item_apply_widget_upper]=ValBox(self,(val_min,val_max),None,val_max,labelwidth)
		self.hist_map_list[idx_hist][self.idx_hist_item_apply_widget_upper].setEnabled(False)
		self.hist_map_list[idx_hist][self.idx_hist_item_apply_widget_upper].text.setStyleSheet("color: rgb(255,0,0);")
		self.hist_map_list[idx_hist][self.idx_hist_item_apply_widget_upper].text.setMinimumSize(QtCore.QSize(editwidth,0))
		self.gbl.addWidget(self.hist_map_list[idx_hist][self.idx_hist_item_apply_widget_upper],grid_row,grid_col_3rd+1)
	
	def readCterCtfFile(self, file_path):
		"""Read all entries from a CTER CTF file into the list box"""
		
		if not os.path.exists(file_path):
			QtGui.QMessageBox.warning(None,"Warning","Can not find CTER CTF File (%s). Please check the file path." % (file_path))
			return
		
		if os.path.basename(file_path).find("partres") == -1:
			QtGui.QMessageBox.warning(None,"Warning","Invalid file name for CTER CTF File (%s). The file name must contain \"partres\"." % (file_path))
			return
		
		if file_path[-1*len(".txt"):] != ".txt":
			QtGui.QMessageBox.warning(None,"Warning","Invalid file extension for CTER CTF File (%s). The file extension must be \".txt\"." % (file_path))
			return
		
#		if os.path.dirname(file_path)[-1*len("partres"):] != "partres":
#			QtGui.QMessageBox.warning(None,"Warning","Invalid file path for CTER CTF File (%s). The file must be in \"partres\" directory." % (file_path))
#			return
		
		new_entry_list = read_text_row(file_path)
		if len(new_entry_list) == 0:
			QtGui.QMessageBox.warning(self, "Warning", "Specified CTER CTF file (%s) does not contain any entry. Please check the file." % (file_path))
			return
		
		cter_pwrot_dir = os.path.join(os.path.dirname(file_path), "pwrot")
		if not os.path.exists(cter_pwrot_dir):
			QtGui.QMessageBox.warning(self, "Warning", "Can not find \"%s\" sub-directory associated with specified CTER CTF file (%s). Please check your project directory." % (cter_pwrot_dir, file_path))
			return
		
		# print "MRK_DEBUG: Detected %s entries in %s" % (len(new_entry_list), file_path)
		# print "MRK_DEBUG: Num. of Columns is %d in %s" % (len(new_entry_list[0]), file_path)
		if len(new_entry_list[0]) == self.n_idx_cter - self.n_idx_cter_extra:
			# This CTEF file format is original one (before around 2016/01/29)
			for cter_id in xrange(len(new_entry_list)):
				# Add extra items first to make sure indices match
				new_entry_list[cter_id] = [cter_id] +  new_entry_list[cter_id]           # self.idx_cter_id , <extra> entry id
				new_entry_list[cter_id] = [1] + new_entry_list[cter_id]                  # self.idx_cter_select  <extra> selected state
#				new_entry_list[cter_id] = new_entry_list[cter_id] + [""]                 # self.idx_cter_pwrot_name, <extra> CTER power spectrum rotational average file name
#				new_entry_list[cter_id] = new_entry_list[cter_id] + [0.5]                # self.idx_cter_error_ctf, <extra> limit frequency by CTF error 
				if self.is_enable_max_power == True: new_entry_list[cter_id] = new_entry_list[cter_id] + [0.0] # MRK_TEST: self.idx_cter_max_power, <extra> maximum power in experimental rotational average (with astigmatism)
				
#				# Cut off frequency components higher than CTF limit 
#				cter_box_size = 512 # NOTE: Toshio Moriya 2016/03/15: This is temporary. Remove this after adding CTF limit to cter output file
#				cter_def = new_entry_list[cter_id][self.idx_cter_def]
#				cter_cs = new_entry_list[cter_id][self.idx_cter_cs]
#				cter_vol = new_entry_list[cter_id][self.idx_cter_vol]
#				cter_apix = new_entry_list[cter_id][self.idx_cter_apix]
#				cter_limit_ab_freq, cter_limit_freq = ctflimit(cter_box_size, cter_def, cter_cs, cter_vol, cter_apix)
#				# NOTE: 2015/12/16 Toshio Moriya
#				# Limiting_frequency[cycle/A]: xr[cycle/A]. Max is Nyquist frequency = 1.0[cycle]/(2*apix[A/pixel]). <UNIT: [cycle/(A/pixel)/[pixel])] => [(cycle*pixel)/(A*pixel] => [cycle/A]>
#				# limiting_period(Angstrom resolution)[A/cycle]: 1.0/xr[cycle/A]. Min is Nyquist period = (2*apix[A/pixel]). <UNIT: [1/(cycle/A)] = [A/cycle]>
#				# Width of Fourier pixel [pixel/A]: fwpix = 1.0[pixel]/(2*apix[A/pixel])/box_half[pixel] = 1[pixel]/fullsize[A]) <UNIT: [pixel/(A/pixel)/(pixel)] = [pixel*(pixel/A)*(1/pixel) = [pixel/A]>
#				# Limiting_absolute_frequency [cycle/pixel] int(xr/fwpix+0.5) = <Unit:[(cycle/A)/(pixel/A)] = [(cycle/A)*(A/pixel)] = [cycle/pixel]>
#				# return  Limiting_abs_frequency [cycle/pixel]Limiting_frequency[cycle/A] <= int(xr/fwpix+0.5),xr
#				new_entry_list[cter_id][self.idx_cter_error_ctf] = cter_limit_freq
				
				# Set associated pwrot file path 
#				assert os.path.dirname(file_path).find("partres") != -1 # MRK_DEBUG
#				cter_pwrot_dir = os.path.dirname(file_path).replace("partres", "pwrot")
#				cter_pwrot_dir = os.path.join(os.path.dirname(file_path), "pwrot")
#				new_cter_pwrot_file_path = os.path.join(cter_pwrot_dir, "rotinf%04d.txt" % new_entry_list[cter_id][self.idx_cter_id])
#				new_entry_list[cter_id][self.idx_cter_pwrot_name] = new_cter_pwrot_file_path
				
				# MRK_TEST: Set max value of pwrot related to this micrograph
				if self.is_enable_max_power == True: # MRK_TEST: 
					new_cter_mic_file_path = new_entry_list[cter_id][self.idx_cter_mic_name]
#					print "MRK_DEBUG: new_cter_mic_file_path := ", new_cter_mic_file_path
					mic_basename_root = os.path.splitext(os.path.basename(new_cter_mic_file_path))
#					print "MRK_DEBUG: mic_basename_root := ", mic_basename_root
					new_cter_pwrot_file_path = os.path.join(os.path.dirname(file_path), "pwrot", "%s_rotinf.txt" % (mic_basename_root))
#					print "MRK_DEBUG: new_cter_pwrot_file_path := ", new_cter_pwrot_file_path
					new_rotinf_table = read_text_file(new_cter_pwrot_file_path, ncol=-1) # MRK_TEST: 
					new_entry_list[cter_id][self.idx_cter_max_power] = max(new_rotinf_table[self.idx_rotinf_exp_with_astig]) # MRK_TEST: 
				
				# Always set selection state to 1 (selected)
				new_entry_list[cter_id][self.idx_cter_select] = 1
		else: 
			# This CTEF file format must be current (after around 2016/01/29) or output of this script
			assert len(new_entry_list[0]) == self.n_idx_cter, "MRK_DEBUG: The number of columns (%d) have to be %d or %d in %s" % (len(new_entry_list[0]), self.n_idx_cter - self.n_idx_cter_extra, self.n_idx_cter, file_path)
		
		# now set the new status
		self.cter_partres_file_path = file_path
		self.cter_entry_list = new_entry_list
		
		# Set the values and ranges of thresholds
		for idx_hist in xrange(self.n_idx_hist):
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
		
		# Set disable status of histogram
		if self.hist_map_list[self.curhist][self.idx_hist_item_val_min] == self.hist_map_list[self.curhist][self.idx_hist_item_val_max]:
			idx_cter = self.hist_map_list[self.curhist][self.idx_hist_item_idx_cter]
			self.curhistdisable=True
			if self.whistparam.isVisible():
				self.whistparam.hide()
			if self.wplotparam.isVisible():
				self.wplotparam.hide()
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
		self.pbreaplysort.setEnabled(True)
		self.pbapplyallthreshold.setEnabled(True)
		self.pbsavethresholdset.setEnabled(True)
		self.pbloadthresholdset.setEnabled(True)
		self.pbsaveselection.setEnabled(True)
		
#		# Always disable micrograph display at the beginning of loading new dataset
#		if self.wimgmicthumb.isVisible() == True:
#			self.wimgmicthumb.hide()
#		self.cbmicthumbdisplay.setValue(False)
#		self.curimgmicthumbdisplay = False
		
		cter_micthumb_dir = os.path.join(os.path.dirname(self.cter_partres_file_path), "micthumb")
		# print "MRK_DEBUG: cter_micthumb_dir = \"%s\" in readCterCtfFile() "% (cter_micthumb_dir)
		if os.path.exists(cter_micthumb_dir):
			# if self.cbmicthumbdisplay.getEnabled() == False: # MRK_NOTE: 2016/03/22 Toshio Moriya: This method does not work as I expected
			self.cbmicthumbdisplay.setEnabled(True)
		else:
			# if self.cbmicthumbdisplay.getEnabled() == True: # MRK_NOTE: 2016/03/22 Toshio Moriya: This method does not work as I expected
			self.cbmicthumbdisplay.setEnabled(False)
			if self.curimgmicthumbdisplay == True:
				self.curimgmicthumbdisplay = False
				self.cbmicthumbdisplay.setValue(self.curimgmicthumbdisplay)
			if self.wimgmicthumb.isVisible():
				self.wimgmicthumb.hide()
			# Error message of this condition should be displayed at the end of this function for smooth visual presentation
			
		# NOTE: 2016/01/03 Toshio Moriya
		# Force update related plots to hide too much scaling delay...
		self.updateImgMicThumb(False)
		self.updateHist()
		self.updatePlotParam()
		
		self.needredisp = True
		
		if self.hist_map_list[self.curhist][self.idx_hist_item_val_min] == self.hist_map_list[self.curhist][self.idx_hist_item_val_max]:
			param_label = self.value_map_list[idx_cter][self.idx_cter_item_label]
			QtGui.QMessageBox.information(self, "Information","All entries have the same selected paramter values (%s). \n\nParameter Histogram & Plot will not be shown" % (param_label))
		
		if not os.path.exists(cter_micthumb_dir):
			QtGui.QMessageBox.warning(None,"Warning","Can not find \"%s\" sub-directory associated with specified CTER CTF file (%s). Please check your project directory. \n\nMicrograph thumbnail display option is disabled for this session." % (cter_micthumb_dir, self.cter_partres_file_path))
		
#		assert(self.isVisible()) 
#		self.raise_()
#		self.activateWindow()
	
	def openCter(self,val=None):
		"""Open CTER CTF file"""
		
		file_path = str(QtGui.QFileDialog.getOpenFileName(self, "Open CTER CTF File", options = QtGui.QFileDialog.DontUseNativeDialog))
		if file_path == "": return
		
		self.readCterCtfFile(os.path.relpath(file_path))
	
	def updateHist(self, error_display = True):
		if self.whistparam == None: return # it's closed/not visible
		if self.cter_partres_file_path == None: return # no cter ctf file is selected
		if self.cter_entry_list == None: return # no cter ctf file is selected
		if self.curhistdisable == True: return # do nothing while it is hidden
		
		val_list = []
		
		# Create Histogram for selected paramter
		idx_cter = self.hist_map_list[self.curhist][self.idx_hist_item_idx_cter]
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
		assert len(val_list) >= n_bin
		assert n_bin > 0
		from statistics import hist_list
		hist_x_list, hist_y_list = hist_list(val_list, n_bin)
		
		# Pad with zero for better visual impression...
		hist_x_list += [max(val_list)]
		hist_y_list += [0]
		self.whistparam.set_data((hist_x_list,hist_y_list),"hist_param",quiet=False,color=0,linetype=0,symtype=0)
		
		# MRK_NOTE: 2015/12/17 Toshio Moriya
		# This may NOT be good place to update the following information...
		idx_cter = self.hist_map_list[self.curhist][self.idx_hist_item_idx_cter]
		param_label = self.value_map_list[idx_cter][self.idx_cter_item_label]
		param_val = self.cter_entry_list[self.curentry][idx_cter]
		
		# shape_name = "hist_param_shape_value"
		# scr_x, scr_y = self.whistparam.plot2scr(param_val, 0.0)
		# self.whistparam.add_shape(shape_name,EMShape(("scrline",0,1,0,scr_x,self.whistparam.scrlim[1],scr_x,self.whistparam.scrlim[1]+self.whistparam.scrlim[3],3)))
		
		val_min = round(min(hist_y_list), self.round_ndigits)
		val_max = round(max(hist_y_list), self.round_ndigits)
		# threshold_lower_label = "Lower(Blue)"
		unapply_threshold_lower_val = round(self.hist_map_list[self.curhist][self.idx_hist_item_unapply_threshold_lower], self.round_ndigits)
		apply_threshold_lower_val = round(self.hist_map_list[self.curhist][self.idx_hist_item_apply_threshold_lower], self.round_ndigits)
		# threshold_upper_label = "Upper(Red)"
		unapply_threshold_upper_val = round(self.hist_map_list[self.curhist][self.idx_hist_item_unapply_threshold_upper], self.round_ndigits)
		apply_threshold_upper_val = round(self.hist_map_list[self.curhist][self.idx_hist_item_apply_threshold_upper], self.round_ndigits)
		
		self.whistparam.set_data(([param_val, param_val], [val_min, val_max]),"selected_val",quiet=False,color=3)
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
		# if self.curhist in [self.idx_hist_error_def, self.idx_hist_error_astig, self.idx_hist_error_ctf] and param_val > 0.0:
		# 	self.whistparam.add_shape(shape_name,EMShape(("scrlabel",0,0,0,self.whistparam.scrlim[0]+30,self.whistparam.scrlim[1]+self.whistparam.scrlim[3]-18,"%s(Green) %1.3g (%1.3g), %s %1.3g, %s %1.3g"%(param_label,param_val,1/param_val,threshold_lower_label,apply_threshold_lower_val,threshold_upper_label,apply_threshold_upper_val),120.0,-1)))
		# 	# self.whistparam.add_shape(shape_name,EMShape(("scrlabel",0,0,0,self.whistparam.scrlim[0]+30,self.whistparam.scrlim[1]+self.whistparam.scrlim[3]-18,"%s %1.5g (%1.5g)"%(param_label,param_val,1/param_val),120.0,-1)))
		# else:
		# 	self.whistparam.add_shape(shape_name,EMShape(("scrlabel",0,0,0,self.whistparam.scrlim[0]+30,self.whistparam.scrlim[1]+self.whistparam.scrlim[3]-18,"%s(Green) %1.3g, %s %1.3g, %s %1.3g"%(param_label,param_val,threshold_lower_label,apply_threshold_lower_val,threshold_upper_label,apply_threshold_upper_val),120.0,-1)))
		# 	# self.whistparam.add_shape(shape_name,EMShape(("scrlabel",0,0,0,self.whistparam.scrlim[0]+30,self.whistparam.scrlim[1]+self.whistparam.scrlim[3]-18,"%s %1.5g"%(param_label,param_val),120.0,-1)))
		
		self.whistparam.setAxisParms(param_label,"Image Counts")
		# x_margin = (hist_x_list[-1] - hist_x_list[0]) * 0.05
		# NOTE: 2016/01/02 Toshio Moriya
		# Disable manual rescale for now and use autoscale
		# self.whistparam.rescale(min(val_list),max(val_list),0,max(hist_y_list) * 1.05)
		self.whistparam.autoscale(True)
	
	def updatePlotParam(self, error_display = True):
		if self.wplotparam == None: return # it's closed/not visible
		if self.cter_partres_file_path == None: return # no cter ctf file is selected
		if self.cter_entry_list == None: return # no cter ctf file is selected
		if self.curhistdisable == True: return # do nothing while it is hidden
		
		x_list = []
		y_list = []
		
		# Create graph for selected paramter
		idx_cter = self.hist_map_list[self.curhist][self.idx_hist_item_idx_cter]
		for cter_id in xrange(len(self.cter_entry_list)):
			x_list.append(cter_id)
			y_list.append(self.cter_entry_list[cter_id][idx_cter])
		# self.wplotparam.set_data((x_list,y_list),"plot_param",quiet=False,color=0)
		self.wplotparam.set_data((x_list,y_list),"plot_param",quiet=False,color=0,linetype=0,symtype=0)
		
		# Create graph for single paramter value of selected entry
		# MRK_NOTE: 2015/12/17 Toshio Moriya
		# This may NOT be good place to update the following information...
		idx_cter = self.hist_map_list[self.curhist][self.idx_hist_item_idx_cter]
		param_label = self.value_map_list[idx_cter][self.idx_cter_item_label]
		param_val = round(self.cter_entry_list[self.curentry][idx_cter], self.round_ndigits)
		
		# threshold_lower_label = "Lower(Blue)"
		unapply_threshold_lower_val = round(self.hist_map_list[self.curhist][self.idx_hist_item_unapply_threshold_lower], self.round_ndigits)
		apply_threshold_lower_val = round(self.hist_map_list[self.curhist][self.idx_hist_item_apply_threshold_lower], self.round_ndigits)
		# threshold_upper_label = "Upper(Red)"
		unapply_threshold_upper_val = round(self.hist_map_list[self.curhist][self.idx_hist_item_unapply_threshold_upper], self.round_ndigits)
		apply_threshold_upper_val = round(self.hist_map_list[self.curhist][self.idx_hist_item_apply_threshold_upper], self.round_ndigits)
		
		y_list = [param_val]*len(x_list)
		self.wplotparam.set_data((x_list,y_list),"selected_val",quiet=False,color=3)
		y_list = [unapply_threshold_lower_val]*len(x_list)
		self.wplotparam.set_data((x_list,y_list),"unapply_threshold_lower_val",quiet=False,color=1,linetype=1)
		y_list = [apply_threshold_lower_val]*len(x_list)
		self.wplotparam.set_data((x_list,y_list),"apply_threshold_lower_val",quiet=False,color=1)
		y_list = [unapply_threshold_upper_val]*len(x_list)
		self.wplotparam.set_data((x_list,y_list),"unapply_threshold_upper_val",quiet=False,color=2,linetype=1)
		y_list = [apply_threshold_upper_val]*len(x_list)
		self.wplotparam.set_data((x_list,y_list),"apply_threshold_upper_val",quiet=False,color=2)
		
		# shape_name = "plot_param_shape_label"
		# if self.curhist in [self.idx_hist_error_def, self.idx_hist_error_astig, self.idx_hist_error_ctf] and param_val > 0.0:
		# 	self.wplotparam.add_shape(shape_name,EMShape(("scrlabel",0,0,0,self.wplotparam.scrlim[0]+30,self.wplotparam.scrlim[1]+self.wplotparam.scrlim[3]-18,"%s(Green) %1.3g (%1.3g), %s %1.3g, %s %1.3g"%(param_label,param_val,1/param_val,threshold_lower_label,apply_threshold_lower_val,threshold_upper_label,apply_threshold_upper_val),120.0,-1)))
		# 	# self.wplotparam.add_shape(shape_name,EMShape(("scrlabel",0,0,0,self.wplotparam.scrlim[0]+30,self.wplotparam.scrlim[1]+self.wplotparam.scrlim[3]-18,"%s(Green) %1.5g (%1.5g)"%(param_label,param_val,1/param_val),120.0,-1)))
		# else:
		# 	self.wplotparam.add_shape(shape_name,EMShape(("scrlabel",0,0,0,self.wplotparam.scrlim[0]+30,self.wplotparam.scrlim[1]+self.wplotparam.scrlim[3]-18,"%s(Green) %1.3g, %s %1.3g, %s %1.3g"%(param_label,param_val,threshold_lower_label,apply_threshold_lower_val,threshold_upper_label,apply_threshold_upper_val),120.0,-1)))
		# 	# self.wplotparam.add_shape(shape_name,EMShape(("scrlabel",0,0,0,self.wplotparam.scrlim[0]+30,self.wplotparam.scrlim[1]+self.wplotparam.scrlim[3]-18,"%s(Green) %1.5g"%(param_label,param_val),120.0,-1)))
		
		self.wplotparam.setAxisParms("Sorted Image ID", param_label)
		# NOTE: 2016/01/02 Toshio Moriya
		# Use autoscale for now
		self.wplotparam.autoscale(True)
	
	def updatePlot(self, error_display = True):
		if self.wplotrotavgcoarse == None: return # it's closed/not visible
		if self.wplotrotavgfine == None: return # it's closed/not visible
		if self.cter_pwrot_file_path == None: return # no cter entry is selected
		
		# Now update the plots
		if not os.path.exists(self.cter_pwrot_file_path):
			if self.wplotrotavgcoarse.isVisible():
				self.wplotrotavgcoarse.hide()
			if self.wplotrotavgfine.isVisible():
				self.wplotrotavgfine.hide()
			if error_display:
				QtGui.QMessageBox.warning(None,"Warning","Can not find file cter_pwrot_file_path (%s). Please check the contents of pwrot directory. \n\nPlots will not be shown." % (self.cter_pwrot_file_path))
			return
			
		self.rotinf_table = read_text_file(self.cter_pwrot_file_path, ncol=-1)
		
		# print "MRK_DEBUG: Last entry of the 1st colum should be a micrograph name %s which is same as " % os.path.basename(self.rotinf_table[0][-1])
		
#		mic_basename_rotinf = os.path.basename(self.rotinf_table[0][-1]) # last entry of 1st colum should be associated micrograph
#		mic_basename_partres = os.path.basename(self.cter_entry_list[self.curentry][self.idx_cter_mic_name])
#		
#		if mic_basename_rotinf != mic_basename_partres:
#			QtGui.QMessageBox.warning(None,"Warning","Micrograph name (%s) in %s is not same as name (%s) in %s " % (mic_basename_rotinf, os.path.basename(self.cter_pwrot_file_path), mic_basename_partres, os.path.basename(self.cter_partres_file_path)))
#			return
		
		# global_min = float("inf")
		# global_max = float("-inf")
		for idx_graph in xrange(self.n_idx_graph):
			self.wplotrotavgcoarse.set_data((self.rotinf_table[self.idx_rotinf_freq],self.rotinf_table[self.graph_map_list[idx_graph][self.idx_graph_idx_rotinf]]),self.graph_map_list[idx_graph][self.idx_graph_item_name],quiet=False,color=idx_graph)
			self.wplotrotavgfine.set_data((self.rotinf_table[self.idx_rotinf_freq],self.rotinf_table[self.graph_map_list[idx_graph][self.idx_graph_idx_rotinf]]),self.graph_map_list[idx_graph][self.idx_graph_item_name],quiet=False,color=idx_graph)
			# val_min = min(self.rotinf_table[self.graph_map_list[idx_graph][self.idx_graph_idx_rotinf]])
			# val_max = max(self.rotinf_table[self.graph_map_list[idx_graph][self.idx_graph_idx_rotinf]])
			# if global_min > val_min:
			# 	global_min = val_min
			# if global_max < val_max:
			# 	global_max = val_max
		
		# NOTE: 2016/01/02 Toshio Moriya
		# Disable manual rescale for now and use autoscale
		# self.wplotrotavgcoarse.rescale(self.rotinf_table[self.idx_rotinf_freq][0],self.rotinf_table[self.idx_rotinf_freq][-1],0.0,1.0)
		self.wplotrotavgcoarse.autoscale(True)
		self.wplotrotavgfine.rescale(self.rotinf_table[self.idx_rotinf_freq][0],self.rotinf_table[self.idx_rotinf_freq][-1],0.0,self.curplotfixscale)
		
		nyquist_freq = self.rotinf_table[self.idx_rotinf_freq][-1]
		# print "MRK_DEBUG: nyquist_freq = %1.5g" % nyquist_freq
		
		error_name = "error_astig"
		error_label = "Astig. Limit"
		error_freq = self.cter_entry_list[self.curentry][self.idx_cter_error_astig]
		# print "MRK_DEBUG: %s= %1.5g" % (error_name, error_freq)
		if error_freq > 0.0 and error_freq <= nyquist_freq:
			error_scr_x, error_scr_y = self.wplotrotavgcoarse.plot2scr(error_freq, 0.0)
			self.wplotrotavgcoarse.add_shape(error_name,EMShape(("scrline",0,0,0.5,error_scr_x,self.wplotrotavgcoarse.scrlim[1],error_scr_x,self.wplotrotavgcoarse.scrlim[1]+self.wplotrotavgcoarse.scrlim[3],1)))
			# self.wplotrotavgcoarse.set_data(([error_freq, error_freq], [global_min, global_max]),"astig_error_freq_limit",quiet=False,color=0,linetype=0)
			self.wplotrotavgcoarse.add_shape("%s_freq"%(error_name),EMShape(("scrlabel",0,0,0,error_scr_x-260,self.wplotrotavgcoarse.scrlim[1]+self.wplotrotavgcoarse.scrlim[3]-18,"%s %1.5g (%1.5g)"%(error_label,error_freq,1.0/error_freq),120.0,-1)))
			error_scr_x, error_scr_y = self.wplotrotavgfine.plot2scr(error_freq, 0.0)
			self.wplotrotavgfine.add_shape(error_name,EMShape(("scrline",0,0,0.5,error_scr_x,self.wplotrotavgfine.scrlim[1],error_scr_x,self.wplotrotavgfine.scrlim[1]+self.wplotrotavgfine.scrlim[3],1)))
			# self.wplotrotavgfine.set_data(([error_freq, error_freq], [global_min, global_max]),"astig_error_freq_limit",quiet=False,color=0,linetype=0)
			self.wplotrotavgfine.add_shape("%s_freq"%(error_name),EMShape(("scrlabel",0,0,0,error_scr_x-260,self.wplotrotavgfine.scrlim[1]+self.wplotrotavgfine.scrlim[3]-18,"%s %1.5g (%1.5g)"%(error_label,error_freq,1.0/error_freq),120.0,-1)))
		
		error_name = "error_def"
		error_label = "Defocus Limit"
		error_freq = self.cter_entry_list[self.curentry][self.idx_cter_error_def]
		# print "MRK_DEBUG: %s= %1.5g" % (error_name, error_freq)
		if error_freq > 0.0 and error_freq <= nyquist_freq:
			error_scr_x, error_scr_y = self.wplotrotavgcoarse.plot2scr(error_freq, 0.0)
			self.wplotrotavgcoarse.add_shape(error_name,EMShape(("scrline",0.5,0,0,error_scr_x,self.wplotrotavgcoarse.scrlim[1],error_scr_x,self.wplotrotavgcoarse.scrlim[1]+self.wplotrotavgcoarse.scrlim[3],1)))
			# self.wplotrotavgcoarse.set_data(([error_freq, error_freq], [global_min, global_max]),"defocus_error_freq_limit",quiet=False,color=0,linetype=0)
			self.wplotrotavgcoarse.add_shape("%s_freq"%(error_name),EMShape(("scrlabel",0,0,0,error_scr_x-260,self.wplotrotavgcoarse.scrlim[1]+self.wplotrotavgcoarse.scrlim[3]-36,"%s %1.5g (%1.5g)"%(error_label,error_freq,1.0/error_freq),120.0,-1)))
			error_scr_x, error_scr_y = self.wplotrotavgfine.plot2scr(error_freq, 0.0)
			self.wplotrotavgfine.add_shape(error_name,EMShape(("scrline",0.5,0,0,error_scr_x,self.wplotrotavgfine.scrlim[1],error_scr_x,self.wplotrotavgfine.scrlim[1]+self.wplotrotavgfine.scrlim[3],1)))
			# self.wplotrotavgfine.set_data(([error_freq, error_freq], [global_min, global_max]),"defocus_error_freq_limit",quiet=False,color=0,linetype=0)
			self.wplotrotavgfine.add_shape("%s_freq"%(error_name),EMShape(("scrlabel",0,0,0,error_scr_x-260,self.wplotrotavgfine.scrlim[1]+self.wplotrotavgfine.scrlim[3]-36,"%s %1.5g (%1.5g)"%(error_label,error_freq,1.0/error_freq),120.0,-1)))
		
		error_name = "error_ctf"
		error_label = "CTF Limit"
		error_freq = self.cter_entry_list[self.curentry][self.idx_cter_error_ctf]
		# print "MRK_DEBUG: %s= %1.5g" % (error_name, error_freq)
		if error_freq > 0.0 and error_freq <= nyquist_freq:
			error_scr_x, error_scr_y = self.wplotrotavgcoarse.plot2scr(error_freq, 0.0)
			self.wplotrotavgcoarse.add_shape(error_name,EMShape(("scrline",0,0.5,0,error_scr_x,self.wplotrotavgcoarse.scrlim[1],error_scr_x,self.wplotrotavgcoarse.scrlim[1]+self.wplotrotavgcoarse.scrlim[3],1)))
			# self.wplotrotavgcoarse.set_data(([error_freq, error_freq], [global_min, global_max]),"ctf_freq_limit")
			self.wplotrotavgcoarse.add_shape("%s_freq"%(error_name),EMShape(("scrlabel",0,0,0,error_scr_x-260,self.wplotrotavgcoarse.scrlim[1]+self.wplotrotavgcoarse.scrlim[3]-54,"%s %1.5g (%1.5g)"%(error_label,error_freq,1.0/error_freq),120.0,-1)))
			error_scr_x, error_scr_y = self.wplotrotavgfine.plot2scr(error_freq, 0.0)
			self.wplotrotavgfine.add_shape(error_name,EMShape(("scrline",0,0.5,0,error_scr_x,self.wplotrotavgfine.scrlim[1],error_scr_x,self.wplotrotavgfine.scrlim[1]+self.wplotrotavgfine.scrlim[3],1)))
			# self.wplotrotavgfine.set_data(([error_freq, error_freq], [global_min, global_max]),"ctf_freq_limit")
			self.wplotrotavgfine.add_shape("%s_freq"%(error_name),EMShape(("scrlabel",0,0,0,error_scr_x-260,self.wplotrotavgfine.scrlim[1]+self.wplotrotavgfine.scrlim[3]-54,"%s %1.5g (%1.5g)"%(error_label,error_freq,1.0/error_freq),120.0,-1)))
		
		self.wplotrotavgcoarse.setAxisParms("frequency (1/"+ "$\AA$" +")","power spectrum")
		self.wplotrotavgfine.setAxisParms("frequency (1/"+ "$\AA$" +")","power spectrum")
		
		self.updatePlotVisibility()
	
	def updateEntryList(self):
		"""Updated entry list box after sorting of CTER entries based on current setting."""
		
		# sort CTER entry list
		assert (self.cter_entry_list != None)
		self.cter_entry_list = sorted(self.cter_entry_list, key=lambda x: x[self.sort_map_list[self.cursort][self.idx_sort_item_idx_cter]], reverse=self.cursortoder)
		
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
				assert(cter_entry[self.idx_cter_select] == 0)
				newItem.setCheckState(Qt.Unchecked)
			self.lbentry.addItem(newItem)
			# self.lbentry.addItem(os.path.basename(cter_entry[self.idx_cter_mic_name]))
			
		self.newEntry(0)
		self.lbentry.setCurrentRow(0)

	def updateImgMicThumb(self, error_display = True):
		if not self.curimgmicthumbdisplay: return # Micrograph thumbnail display is disabled
		if self.wimgmicthumb == None: return # it's closed/not visible
		if self.cter_micthumb_file_path == None: return # no cter entry is selected
		
		# print "MRK_DEBUG: self.cter_micthumb_file_path =\"%s\" in updateImgMicThumb() "% (self.cter_micthumb_file_path)
		if not os.path.exists(self.cter_micthumb_file_path):
			# QtGui.QMessageBox.warning(None,"Warning","Can not find micrograph thumbnail (%s). Please check your micrograph thumbnail directory. \n\nA blank image is shown." % (self.cter_micthumb_file_path))
			if self.wimgmicthumb.isVisible():
				self.wimgmicthumb.hide()
			if error_display:
				QtGui.QMessageBox.warning(None,"Warning","Can not find micrograph thumbnail (%s). Please check your micrograph thumbnail directory. \n\nMicrograph thumbnail will not be shown." % (self.cter_micthumb_file_path))
			return
		
		micthumb_img = EMData(self.cter_micthumb_file_path) # read the image from disk
		self.wimgmicthumb.set_data(micthumb_img)
		self.wimgmicthumb.setWindowTitle("sxgui_cter - Micrograph Thumbnail- %s, %s" % (os.path.basename(self.cter_entry_list[self.curentry][self.idx_cter_mic_name]), os.path.basename(self.cter_pwrot_file_path)))
		
#		# print "MRK_DEBUG: self.cter_mic_file_path =\"%s\" in updateImgMicThumb() "% (self.cter_mic_file_path)
#		if not os.path.exists(self.cter_mic_file_path):
#			QtGui.QMessageBox.warning(None,"Warning","Can not find micrograph (%s). Please check your micrograph directory. \n\nA blank image is shown." % (self.cter_mic_file_path))
#			mic_img = EMData() # Set empty image...
#			img_size = 4096 # Use the most typical image size?!!!
#			mic_img = model_blank(img_size,img_size, bckg=1.0)
#		else:
#			mic_img = EMData(self.cter_mic_file_path) # read the image from disk
#		self.wimage.set_data(mic_img)
#		self.wimage.setWindowTitle("sxgui_cter - Micrograph - %s, %s" % (os.path.basename(self.cter_entry_list[self.curentry][self.idx_cter_mic_name]), os.path.basename(self.cter_pwrot_file_path)))
		
	def newMicThumbDisplay(self,val=None):
		"""Change micrograph thumbnail display status."""
		# assert(self.cbmicthumbdisplay.getEnabled() == True)  # MRK_NOTE: 2016/03/22 Toshio Moriya: This method does not work as I expected
		
		if self.curimgmicthumbdisplay == val: return
		
		# now set the new status
		self.curimgmicthumbdisplay = val
		
		if self.curimgmicthumbdisplay and not self.wimgmicthumb.isVisible():
			self.wimgmicthumb.show()
			self.needredisp = True
		elif not self.curimgmicthumbdisplay and self.wimgmicthumb.isVisible():
			self.wimgmicthumb.hide()
		
	def newEntry(self,currow):
		"""called when a new data set is selected from the CTER Entry list box."""
		assert(self.cter_partres_file_path != None)
		assert(self.cter_entry_list != None)
		
		# always update the current row of cter entry list 
		# to get associated micrograph path and pwrot file path 
		self.curentry = currow # row can be the same even after resorting of the cter entry list
		
		# Get associated pwrot file path of current entry
#		new_cter_pwrot_file_path = self.cter_entry_list[self.curentry][self.idx_cter_pwrot_name]
		
		# Get associated micrograph path of current entry
		new_cter_mic_file_path = self.cter_entry_list[self.curentry][self.idx_cter_mic_name]
		
		# Generate associated micthumb & pwrot file path of current entry
		mic_basename_root, mic_extension = os.path.splitext(os.path.basename(new_cter_mic_file_path))
		new_cter_micthumb_file_path = os.path.join(os.path.dirname(self.cter_partres_file_path), "micthumb", "%s_thumb%s" % (mic_basename_root, mic_extension))
		new_cter_pwrot_file_path = os.path.join(os.path.dirname(self.cter_partres_file_path), "pwrot", "%s_rotinf.txt" % (mic_basename_root))
		
		# Changing row does not always change the pwrot file path after resorting of the cter entry list
		# If same, skip the following processes
		if self.cter_pwrot_file_path == new_cter_pwrot_file_path: 
			assert(self.cter_micthumb_file_path == new_cter_micthumb_file_path)
			assert(self.cter_mic_file_path == new_cter_mic_file_path)
			return
		
		# now set the new item
		assert(self.cter_pwrot_file_path != new_cter_pwrot_file_path)
		self.cter_pwrot_file_path = new_cter_pwrot_file_path
		assert(self.cter_micthumb_file_path != new_cter_micthumb_file_path)
		self.cter_micthumb_file_path = new_cter_micthumb_file_path
		assert(self.cter_mic_file_path != new_cter_mic_file_path)
		self.cter_mic_file_path = new_cter_mic_file_path
		
#		# Now update the image
#		self.updateImgMicThumb()
		
		# print "MRK_DEBUG: Row No. %d (CTER Entry No. %d) is selected from cter entry list box" % (self.curentry, self.cter_entry_list[self.curentry][self.idx_cter_id])
		
		self.ssortedid.setValue(self.curentry,True)
		
		for idx_cter in xrange(self.n_idx_cter):
#			if idx_cter not in [self.idx_cter_mic_name, self.idx_cter_pwrot_name]:
			if idx_cter not in [self.idx_cter_mic_name]:
				self.value_map_list[idx_cter][self.idx_cter_item_widget].setValue(self.cter_entry_list[self.curentry][idx_cter],True)
		
		# Use red font to indicate the value is not between applied threshold ranage
		for idx_hist in xrange(self.n_idx_hist):
			lower_threshold = round(self.hist_map_list[idx_hist][self.idx_hist_item_apply_threshold_lower], self.round_ndigits)
			upper_threshold = round(self.hist_map_list[idx_hist][self.idx_hist_item_apply_threshold_upper], self.round_ndigits)
			idx_cter = self.hist_map_list[idx_hist][self.idx_hist_item_idx_cter]
			param_val = round(self.cter_entry_list[self.curentry][idx_cter], self.round_ndigits)
			if lower_threshold <= param_val and param_val <= upper_threshold:
				self.value_map_list[idx_cter][self.idx_cter_item_widget].text.setStyleSheet("color: rgb(0,0,0);")
			elif param_val < lower_threshold:
				self.value_map_list[idx_cter][self.idx_cter_item_widget].text.setStyleSheet("color: rgb(0,0,255);")
			else:
				assert(upper_threshold < param_val)
				self.value_map_list[idx_cter][self.idx_cter_item_widget].text.setStyleSheet("color: rgb(255,0,0);")
		
#		self.wfft.setWindowTitle("sxgui_cter - 2D FFT - "+fsp.split("/")[-1])
		self.wplotrotavgcoarse.setWindowTitle("sxgui_cter - Plot - %s, %s" % (os.path.basename(self.cter_entry_list[self.curentry][self.idx_cter_mic_name]), os.path.basename(self.cter_pwrot_file_path)))
		self.wplotrotavgfine.setWindowTitle("sxgui_cter - Plot Zoom- %s, %s" % (os.path.basename(self.cter_entry_list[self.curentry][self.idx_cter_mic_name]), os.path.basename(self.cter_pwrot_file_path)))
		
#		self.cter_mic_file_path = new_cter_mic_file_path
	
		self.needredisp = True
	
	def updateEntrySelect(self, entry):
		"""called when check status of an cter entry in list box is changed."""
		assert(self.cter_partres_file_path != None)
		assert(self.cter_entry_list != None)
		
		newSelect = 1
		if entry.checkState() == Qt.Unchecked:
			newSelect = 0
		entry_row = self.lbentry.row(entry)
		self.cter_entry_list[entry_row][self.idx_cter_select] = newSelect
		
		if self.curentry == entry_row:
			self.value_map_list[self.idx_cter_select][self.idx_cter_item_widget].setValue(self.cter_entry_list[self.curentry][self.idx_cter_select],True)
			
		self.updateUncheckCounts()
	
	def updateUncheckCounts(self):
		"""called whenever checked status of cter entries change."""
		assert(self.cter_partres_file_path != None)
		assert(self.cter_entry_list != None)
		
		assert(len(self.cter_entry_list) > 0)
		n_entry = len(self.cter_entry_list)
		uncheck_counts  = n_entry
		for cter_entry in self.cter_entry_list:
			uncheck_counts -= cter_entry[self.idx_cter_select]
		assert(uncheck_counts >= 0 and uncheck_counts <= n_entry)
		
		self.vbuncheckcounts.setValue(uncheck_counts,True)
		self.vbuncheckratio.setValue(float(uncheck_counts)/n_entry,True)
	
	def reapplySort(self,item = None):
		"""Called when reapply button is clicked."""
		if self.cter_partres_file_path == None: return # no cter ctf file is selected
		if self.cter_entry_list == None: return # no cter ctf file is selected
		
		self.updateEntryList()
	
	def newSort(self,cursort):
		"""Sort CTER entries by selected parameter values."""
		if self.cursort == cursort: return
		
		# now set the new item
		self.cursort = cursort
		
		if self.cter_partres_file_path == None: return # no cter ctf file is selected
		if self.cter_entry_list == None: return # no cter ctf file is selected
		
		self.updateEntryList()
	
	def newSortOrder(self, sortoder):
		"""Change sorting order of CTER entries."""
		if self.cursortoder == sortoder: return
		
		# now set the new status
		self.cursortoder = sortoder
		
		if self.cter_partres_file_path == None: return # no cter ctf file is selected
		if self.cter_entry_list == None: return # no cter ctf file is selected
		
		self.updateEntryList()
	
	def newSortSelect(self, sortselect):
		"""Change sort select status of CTER entries."""
		if self.cursortselect == sortselect: return
		
		# now set the new status
		self.cursortselect = sortselect
		
		if self.cter_partres_file_path == None: return # no cter ctf file is selected
		if self.cter_entry_list == None: return # no cter ctf file is selected
		
		self.updateEntryList()
	
	def newThresholdLower(self):
		threshold_lower = round(self.hist_map_list[self.curhist][self.idx_hist_item_unapply_widget_lower].getValue(), self.round_ndigits)
		if threshold_lower < round(self.hist_map_list[self.curhist][self.idx_hist_item_val_min], self.round_ndigits):
			threshold_lower = round(self.hist_map_list[self.curhist][self.idx_hist_item_val_min], self.round_ndigits)
			self.hist_map_list[self.curhist][self.idx_hist_item_unapply_widget_lower].setValue(threshold_lower)
		elif threshold_lower > round(self.hist_map_list[self.curhist][self.idx_hist_item_val_max], self.round_ndigits):
			threshold_lower = round(self.hist_map_list[self.curhist][self.idx_hist_item_val_max], self.round_ndigits)
			self.hist_map_list[self.curhist][self.idx_hist_item_unapply_widget_lower].setValue(threshold_lower)
		# else: # Do nothing
		
		# now set the new threshold
		self.hist_map_list[self.curhist][self.idx_hist_item_unapply_threshold_lower] = threshold_lower
		
		self.needredisp=True
	
	def newThresholdUpper(self):
		threshold_upper = round(self.hist_map_list[self.curhist][self.idx_hist_item_unapply_widget_upper].getValue(), self.round_ndigits)
		if threshold_upper < round(self.hist_map_list[self.curhist][self.idx_hist_item_val_min], self.round_ndigits):
			threshold_upper = round(self.hist_map_list[self.curhist][self.idx_hist_item_val_min], self.round_ndigits)
			self.hist_map_list[self.curhist][self.idx_hist_item_unapply_widget_upper].setValue(threshold_upper)
		elif threshold_upper > round(self.hist_map_list[self.curhist][self.idx_hist_item_val_max], self.round_ndigits):
			threshold_upper = round(self.hist_map_list[self.curhist][self.idx_hist_item_val_max], self.round_ndigits)
			self.hist_map_list[self.curhist][self.idx_hist_item_unapply_widget_upper].setValue(threshold_upper)
		# else: # Do nothing
		
		# now set the new threshold
		self.hist_map_list[self.curhist][self.idx_hist_item_unapply_threshold_upper] = threshold_upper
		
		self.needredisp=True
	
	def newHist(self,currow):
		"called when a new row is selected from the Histogram list box"
		
		if self.curhist == currow: return
		
		# Disable old item
		if self.curthresholdcontrol == self.idx_threshold_control_lower:
			self.hist_map_list[self.curhist][self.idx_hist_item_unapply_widget_lower].setEnabled(False)
		elif self.curthresholdcontrol == self.idx_threshold_control_upper:
			self.hist_map_list[self.curhist][self.idx_hist_item_unapply_widget_upper].setEnabled(False)
		else:
			assert(self.curthresholdcontrol == self.idx_threshold_control_edit_only)
			self.hist_map_list[self.curhist][self.idx_hist_item_unapply_widget_lower].setEnabled(False)
			self.hist_map_list[self.curhist][self.idx_hist_item_unapply_widget_upper].setEnabled(False)
		
		# now set the new item and enalble it
		self.curhist=currow
		# print "MRK_DEBUG: Row No. %d is selected from histogram list box" % (self.curhist)
		
		# Check if the all selected parameter values are same 
		if self.hist_map_list[self.curhist][self.idx_hist_item_val_min] == self.hist_map_list[self.curhist][self.idx_hist_item_val_max]:
			idx_cter = self.hist_map_list[self.curhist][self.idx_hist_item_idx_cter]
			param_label = self.value_map_list[idx_cter][self.idx_cter_item_label]
			self.curhistdisable=True
			if self.whistparam.isVisible():
				self.whistparam.hide()
			if self.wplotparam.isVisible():
				self.wplotparam.hide()
			QtGui.QMessageBox.information(self, "Information","All entries have the same selected paramter values (%s). \n\nParameter Histogram & Plot will not be shown" % (param_label))
		else:
			if self.curthresholdcontrol == self.idx_threshold_control_lower:
				self.hist_map_list[self.curhist][self.idx_hist_item_unapply_widget_lower].setEnabled(True)
			elif self.curthresholdcontrol == self.idx_threshold_control_upper:
				self.hist_map_list[self.curhist][self.idx_hist_item_unapply_widget_upper].setEnabled(True)
			else:
				assert(self.curthresholdcontrol == self.idx_threshold_control_edit_only)
				self.hist_map_list[self.curhist][self.idx_hist_item_unapply_widget_lower].setEnabled(True)
				self.hist_map_list[self.curhist][self.idx_hist_item_unapply_widget_upper].setEnabled(True)
			
			idx_cter = self.hist_map_list[self.curhist][self.idx_hist_item_idx_cter]
			param_label = self.value_map_list[idx_cter][self.idx_cter_item_label]
			self.whistparam.setWindowTitle("sxgui_cter - %s Histogram" % (param_label))
			self.wplotparam.setWindowTitle("sxgui_cter - %s Sort Plot" % (param_label))
			
			if self.cursyncsort == True:
				idx_sort = self.hist_map_list[self.curhist][self.idx_hist_item_idx_sort]
				if (idx_sort != self.cursort):
					self.newSort(idx_sort)
					self.ssort.setCurrentIndex(idx_sort)
				# else: assert(idx_sort == self.cursort) # Do nothing
			# else: assert(self.cursyncsort == False) # Do nothing
			
			if self.cter_partres_file_path == None: return # no cter ctf file is selected
			if self.cter_entry_list == None: return # no cter ctf file is selected
			
			self.curhistdisable=False
			if not self.whistparam.isVisible():
				self.whistparam.show()
			if not self.wplotparam.isVisible():
				self.wplotparam.show()
			
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
			self.hist_map_list[self.curhist][self.idx_hist_item_unapply_widget_lower].setEnabled(False)
		elif self.curthresholdcontrol == self.idx_threshold_control_upper:
			self.hist_map_list[self.curhist][self.idx_hist_item_unapply_widget_upper].setEnabled(False)
		else:
			assert(self.curthresholdcontrol == self.idx_threshold_control_edit_only)
			self.hist_map_list[self.curhist][self.idx_hist_item_unapply_widget_lower].setEnabled(False)
			self.hist_map_list[self.curhist][self.idx_hist_item_unapply_widget_upper].setEnabled(False)
		
		# now set the new item and enalble it
		self.curthresholdcontrol=currow
		# print "MRK_DEBUG: Row No. %d is selected from histogram list box" % (self.curhist)
		
		if self.curthresholdcontrol == self.idx_threshold_control_lower:
			self.hist_map_list[self.curhist][self.idx_hist_item_unapply_widget_lower].setEnabled(True)
		elif self.curthresholdcontrol == self.idx_threshold_control_upper:
			self.hist_map_list[self.curhist][self.idx_hist_item_unapply_widget_upper].setEnabled(True)
		else:
			assert(self.curthresholdcontrol == self.idx_threshold_control_edit_only)
			self.hist_map_list[self.curhist][self.idx_hist_item_unapply_widget_lower].setEnabled(True)
			self.hist_map_list[self.curhist][self.idx_hist_item_unapply_widget_upper].setEnabled(True)
	
	def newSyncSort(self, syncsort):
		"""Change scyn sort enable state."""
		if self.cursyncsort == syncsort: return
		
		# now set the new status
		self.cursyncsort = syncsort
		
		self.ssort.setEnabled(not self.cursyncsort)
		
		if self.cursyncsort == True:
			idx_sort = self.hist_map_list[self.curhist][self.idx_hist_item_idx_sort]
			if (idx_sort != self.cursort):
				self.newSort(idx_sort)
				self.ssort.setCurrentIndex(idx_sort)
			# else: assert(idx_sort == self.cursort) # Do nothing
		# else: assert(self.cursyncsort == False) # Do nothing
	
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
			for idx_hist in xrange(self.n_idx_hist):
				threshold_lower = round(self.hist_map_list[idx_hist][self.idx_hist_item_unapply_threshold_lower], self.round_ndigits)
				threshold_upper = round(self.hist_map_list[idx_hist][self.idx_hist_item_unapply_threshold_upper], self.round_ndigits)
				self.hist_map_list[idx_hist][self.idx_hist_item_apply_threshold_lower] = threshold_lower
				self.hist_map_list[idx_hist][self.idx_hist_item_apply_threshold_upper] = threshold_upper
				self.hist_map_list[idx_hist][self.idx_hist_item_apply_widget_lower].setValue(threshold_lower)
				self.hist_map_list[idx_hist][self.idx_hist_item_apply_widget_upper].setValue(threshold_upper)
				idx_cter = self.hist_map_list[idx_hist][self.idx_hist_item_idx_cter]
				param_val = round(cter_entry[idx_cter], self.round_ndigits)
				if param_val < threshold_lower or threshold_upper < param_val:
					print "MRK_DEBUG: Param #%d diselected entry #%04d with (param_val, threshold_lower, threshold_upper) = (%1.15g, %1.15g, %1.15g)" % (idx_hist, idx_cter, param_val, threshold_lower, threshold_upper)
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
		assert(self.cter_partres_file_path != None)
		assert(self.cter_entry_list != None)
		
		file_out = open(file_path_out,"w")
		
		# Write lines to check consistency upon loading
		file_out.write("# @@@@@ gui_cter thresholds - ")
		# file_out.write(EMANVERSION + " (CVS" + CVSDATESTAMP[6:-2] +")")
		file_out.write(EMANVERSION + " (GITHUB: " + DATESTAMP +")" )
		file_out.write(" @@@@@ \n")
		file_out.write("# Associated CTER CTF File == %s\n" % (self.cter_partres_file_path))
		file_out.write("# Saved Threshold Set == %s\n" % (self.thresholdset_map_list[idx_thresholdset][self.idx_thresholdset_item_label]))
		file_out.write("# [Paramter Id] [Paramter Name] [Lower Threshold] [Upper Threshold]\n")
		
		# Assigne the index of target threshold values
		idx_threshold_lower = self.idx_hist_item_unapply_threshold_lower
		idx_threshold_upper = self.idx_hist_item_unapply_threshold_upper
		if idx_thresholdset == self.idx_thresholdset_applied:
			idx_threshold_lower = self.idx_hist_item_apply_threshold_lower
			idx_threshold_upper = self.idx_hist_item_apply_threshold_upper
		
		for idx_hist in xrange(self.n_idx_hist):
			map_entry = self.hist_map_list[idx_hist]
			idx_cter = map_entry[self.idx_hist_item_idx_cter]
			param_label = self.value_map_list[idx_cter][self.idx_cter_item_label]
			# NOTE: 2016/01/26 Toshio Moriya
			# Use the precision for double to minimise precision loss by save & load operations 
			file_out.write("%2d %s == %1.15g %1.15g \n" % (idx_hist, param_label, round(map_entry[idx_threshold_lower], self.round_ndigits), round(map_entry[idx_threshold_upper], self.round_ndigits)))
		
		file_out.close()
	
	def readThresholdSet(self, file_path_in, idx_thresholdset):
		assert(self.cter_partres_file_path != None)
		assert(self.cter_entry_list != None)
		
		file_in = open(file_path_in,"r")
		
		# Check if this parameter file is threshold
		line_in = file_in.readline()
		if line_in.find("@@@@@ gui_cter thresholds") != -1:
			# loop through the rest of lines
			for line_in in file_in:
				if line_in[0] == "#":
					continue
					
				tokens_in = line_in.split("==")
				assert(len(tokens_in) == 2)
				tokens_label = tokens_in[0].split()
				assert(len(tokens_label) >= 2)
				idx_hist = int(tokens_label[0]) # index of hist_map_list
				map_entry = self.hist_map_list[idx_hist]
				tokens_val = tokens_in[1].split()
				assert(len(tokens_val) == 2)
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
		
		assert(os.path.basename(self.cter_partres_file_path).find("partres") != -1)
		assert(self.cter_partres_file_path[-1*len(".txt"):] == ".txt")
#		assert(os.path.dirname(self.cter_partres_file_path)[-1*len("partres"):] == "partres")
		
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
		idx_cter_ignore_list = [self.idx_cter_id, self.idx_cter_select, self.idx_cter_error_ctf]
		if self.is_enable_max_power == True: idx_cter_ignore_list.append(self.idx_cter_max_power)
		for cter_entry in save_cter_entry_list:
			file_out = file_out_select
			if cter_entry[self.idx_cter_select] == 0:
				file_out = file_out_discard
			# else: assert(cter_entry[self.idx_cter_select] == 1) # do nothing
			
			for idx_cter in xrange(self.n_idx_cter):
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
			# else: assert(cter_entry[self.idx_cter_select] == 1) # do nothing
			file_out.write("  %s\n" % cter_entry[self.idx_cter_mic_name])
		
		file_out_mic_select.close()
		file_out_mic_discard.close()
		
		# Save the associated applied threshold 
		self.writeThresholdSet(file_path_out_thresholds, self.idx_thresholdset_applied) 
		
		QtGui.QMessageBox.information(self, "Information","The following files are saved in %s:\n\nCTER CTF List - Selected: %s\n\nCTER CTF List - Discarded: %s\n\nMicrograph - Selected: %s\n\nMicrograph - Discarded: %s\n\nApplied Threshold Set: %s" % (os.path.dirname(self.cter_partres_file_path), file_path_out_select, file_path_out_discard, file_path_out_mic_select, file_path_out_mic_discard, file_path_out_thresholds))
	
	def timeOut(self):
		if self.busy: return
		
		# Redisplay before spawning thread for more interactive display
		if self.needredisp:
			try:
				self.redisplay()
			except:
				print "Recieved unexpected exception from redisplay() in timeOut(): "
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
#			if not self.wfft.isVisible() and self.is_wfft_minimized:
#				self.wfft.show()
#				is_child_shown = True
			if not self.wimgmicthumb.isVisible() and self.curimgmicthumbdisplay and not self.is_wimgmicthumb_minimized:
				self.wimgmicthumb.show()
				is_child_shown = True
			if not self.wplotrotavgcoarse.isVisible() and not self.is_wplotrotavgcoarse_minimized:
				self.wplotrotavgcoarse.show()
				is_child_shown = True
			if not self.wplotrotavgfine.isVisible() and not self.is_wplotrotavgfine_minimized:
				self.wplotrotavgfine.show()
				is_child_shown = True
			if not self.whistparam.isVisible() and not self.curhistdisable and not self.is_whistparam_minimized:
				self.whistparam.show()
				is_child_shown = True
			if not self.wplotparam.isVisible() and not self.curhistdisable and not self.is_wplotparam_minimized:
				self.wplotparam.show()
				is_child_shown = True
			if is_child_shown == True:
				self.raise_()
				self.activateWindow()
				
		self.updateImgMicThumb()
		self.updateHist()
		self.updatePlotParam()
		self.updatePlot()
		
		self.busy=False
		
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
				assert(self.isMinimized() == True)
				#
				# NOTE: 2016/03/09 Toshio Moriya
				# Minimize icon button of child window should be disabled
				#
				if self.cter_entry_list != None:
					# if self.wfft.isVisible() and not self.is_wfft_minimized:
					# 	self.wfft.hide()
					# 	self.is_wfft_minimized = True
					if self.wimgmicthumb.isVisible() == True and self.is_wimgmicthumb_minimized == False:
						assert(self.curimgmicthumbdisplay == True)
						self.wimgmicthumb.hide()
						self.is_wimgmicthumb_minimized = True
					if self.wplotrotavgcoarse.isVisible() == True and self.is_wplotrotavgcoarse_minimized == False:
						self.wplotrotavgcoarse.hide()
						self.is_wplotrotavgcoarse_minimized = True
					if self.wplotrotavgfine.isVisible() == True and self.is_wplotrotavgfine_minimized == False:
						self.wplotrotavgfine.hide()
						self.is_wplotrotavgfine_minimized = True
					if self.whistparam.isVisible() == True and self.is_whistparam_minimized == False:
						assert(self.curhistdisable == False)
						self.whistparam.hide()
						self.is_whistparam_minimized = True
					if self.wplotparam.isVisible() == True and self.is_wplotparam_minimized == False:
						assert(self.curhistdisable == False )
						self.wplotparam.hide()
						self.is_wplotparam_minimized = True
			else:
				# print "MRK_DEBUG: sxgui main window has not minimized"
				assert(self.isMinimized() == False)
				#
				# NOTE: 2016/03/09 Toshio Moriya
				# Minimize icon button of child window should be disabled
				#
				if self.cter_entry_list != None:
					# if self.is_wfft_minimized == True:
					# 	assert(self.wfft.isVisible() == False and self.curfftdisplay == True)
					# 	self.wfft.show()
					# 	self.is_wfft_minimized = False
					if self.is_wimgmicthumb_minimized == True:
						assert(self.wimgmicthumb.isVisible() == False and self.curimgmicthumbdisplay == True)
						self.wimgmicthumb.show()
						self.is_wimgmicthumb_minimized = False
					if self.is_wplotrotavgcoarse_minimized == True:
						assert(self.wplotrotavgcoarse.isVisible() == False)
						self.wplotrotavgcoarse.show()
						self.is_wplotrotavgcoarse_minimized = False
					if self.is_wplotrotavgfine_minimized == True:
						assert(self.wplotrotavgfine.isVisible() == False)
						self.wplotrotavgfine.show()
						self.is_wplotrotavgfine_minimized = False
					if self.is_whistparam_minimized == True:
						assert(self.whistparam.isVisible() == False and self.curhistdisable == False)
						self.whistparam.show()
						self.is_whistparam_minimized = False
					if self.is_wplotparam_minimized == True:
						assert(self.wplotparam.isVisible() == False and self.curhistdisable == False )
						self.wplotparam.show()
						self.is_wplotparam_minimized = False
				if self.isVisible() == True: # Depends on timing at startup, this can happen?!!
					self.raise_()
					self.activateWindow()
		elif event.type() == QtCore.QEvent.WindowActivate:
			# print "MRK_DEBUG: sxgui main window has gained focus (beome active)"
			# 
			# NOTE: 2016/03/08 Toshio Moriya
			# To raise EMGLWidget (SXPlot2DWidget and EMImage2DWidget) window,
			# we have to go through its qt_parent attribute to call raise_()...
			# 
			if self.cter_entry_list != None:
				# if self.wfft.isVisible() == True:
				#	self.wfft.qt_parent.raise_()
				if self.wimgmicthumb.isVisible() == True:
					self.wimgmicthumb.qt_parent.raise_()
				if self.wplotrotavgcoarse.isVisible() == True:
					self.wplotrotavgcoarse.qt_parent.raise_()
				if self.wplotrotavgfine.isVisible() == True:
					self.wplotrotavgfine.qt_parent.raise_()
				if self.whistparam.isVisible() == True:
					self.whistparam.qt_parent.raise_()
				if self.wplotparam.isVisible() == True:
					self.wplotparam.qt_parent.raise_()
				assert(self.isVisible()) 
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
#			E2saveappwin("sxgui_cter","fft",self.wfft.qt_parent)
			E2saveappwin("sxgui_cter","imgmicthumb",self.wimgmicthumb.qt_parent)
			E2saveappwin("sxgui_cter","plotcoarse",self.wplotrotavgcoarse.qt_parent)
			E2saveappwin("sxgui_cter","plotfine",self.wplotrotavgfine.qt_parent)
			E2saveappwin("sxgui_cter","plotparam",self.wplotparam.qt_parent)
			E2saveappwin("sxgui_cter","histparam",self.whistparam.qt_parent)
		
		# close all child windows
		# if self.wfft:
		# 	self.wfft.close()
		if self.wimgmicthumb: self.wimgmicthumb.close()
		if self.wplotrotavgcoarse: self.wplotrotavgcoarse.close()
		if self.wplotrotavgfine: self.wplotrotavgfine.close()
		if self.whistparam: self.whistparam.close()
		if self.wplotparam: self.wplotparam.close()
		
		event.accept()
		QtGui.qApp.exit(0)
		
	def updatePlotVisibility(self,val=None):
		if self.wplotrotavgcoarse == None: return # it's closed/not visible
		if self.wplotrotavgfine == None: return # it's closed/not visible
		if self.cter_pwrot_file_path == None: return # no cter entry is selected
		
		for idx_graph in xrange(self.n_idx_graph):
			item_widget = self.graph_map_list[idx_graph][self.idx_graph_item_widget]
			name = self.graph_map_list[idx_graph][self.idx_graph_item_name]
			if self.wplotrotavgcoarse.visibility[name] != item_widget.getValue():
				self.wplotrotavgcoarse.visibility[name] = item_widget.getValue()
				self.wplotrotavgcoarse.full_refresh()
				self.wplotrotavgcoarse.updateGL()
		
		for idx_graph in xrange(self.n_idx_graph):
			item_widget = self.graph_map_list[idx_graph][self.idx_graph_item_widget]
			name = self.graph_map_list[idx_graph][self.idx_graph_item_name]
			if self.wplotrotavgfine.visibility[name] != item_widget.getValue():
				self.wplotrotavgfine.visibility[name] = item_widget.getValue()
				self.wplotrotavgfine.full_refresh()
				self.wplotrotavgfine.updateGL()
	
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
		
		plot_x,plot_y=self.wplotparam.scr2plot(event.x(),event.y())
		if self.curthresholdcontrol == self.idx_threshold_control_lower:
			if is_not_reverse_control:
				self.hist_map_list[self.curhist][self.idx_hist_item_unapply_widget_lower].setValue(plot_y)
				self.newThresholdLower()
			else:
				self.hist_map_list[self.curhist][self.idx_hist_item_unapply_widget_upper].setValue(plot_y)
				self.newThresholdUpper()
		else:
			assert(self.curthresholdcontrol == self.idx_threshold_control_upper)
			if is_not_reverse_control:
				self.hist_map_list[self.curhist][self.idx_hist_item_unapply_widget_upper].setValue(plot_y)
				self.newThresholdUpper()
			else:
				self.hist_map_list[self.curhist][self.idx_hist_item_unapply_widget_lower].setValue(plot_y)
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
				self.hist_map_list[self.curhist][self.idx_hist_item_unapply_widget_lower].setValue(hist_x)
				self.newThresholdLower()
			else:
				self.hist_map_list[self.curhist][self.idx_hist_item_unapply_widget_upper].setValue(hist_x)
				self.newThresholdUpper()
				
		else:
			assert(self.curthresholdcontrol == self.idx_threshold_control_upper)
			if is_not_reverse_control:
				self.hist_map_list[self.curhist][self.idx_hist_item_unapply_widget_upper].setValue(hist_x)
				self.newThresholdUpper()
			else:
				self.hist_map_list[self.curhist][self.idx_hist_item_unapply_widget_lower].setValue(hist_x)
				self.newThresholdLower()

if __name__ == "__main__":
	main()
