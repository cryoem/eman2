#!/usr/bin/env python
#
# Author: Toshio Moriya, 2015/12/21 (toshio.moriya@mpi-dortmund.mpg.de)
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

def main():
	from emapplication import EMApp
	app=EMApp()

	gui=SXGuiCter()
	gui.show()

	try:
		gui.wimage.raise_()
#		gui.wfft.raise_()
		gui.whist.raise_()
		gui.wplot.raise_()
		gui.raise_()
	except: 
		print "Recieved unexpected exception in main(): ", sys.exc_info()[0]
		exc_type, exc_value, exc_traceback = sys.exc_info()
		traceback.print_exception(exc_type, exc_value, exc_traceback)
		# MRK_NOTE: 2015/12/17 Toshio Moriya
		# Another way to print out exception info...
		# lines = traceback.format_exception(exc_type, exc_value, exc_traceback)
		# print ''.join('!! ' + line for line in lines)
		pass
	
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
		self.del_shapes(("error_astig","error_astig_freq","error_def","error_def_freq","error_ctf","error_ctf_freq")) # for SXGuiCter.wplot
		self.del_shapes(("param_value","param_value_label")) # for SXGuiCter.whist
		
class SXGuiCter(QtGui.QWidget):
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
		
		QtGui.QWidget.__init__(self,None)
		self.setWindowIcon(QtGui.QIcon(get_image_directory() + "ctf.png"))
		
		self.idx_cter_id           =  0 # <extra> entry id
		self.idx_cter_def          =  1 # defocus (ym)
		self.idx_cter_cs           =  2 # Cs (mm)
		self.idx_cter_vol          =  3 # voltage(kV)
		self.idx_cter_apix         =  4 # pixel size (A)
		self.idx_cter_bfactor      =  5 # B-factor (A^2)
		self.idx_cter_ac           =  6 # amp contrast (%)
		self.idx_cter_astig_amp    =  7 # astigmatism amplitude (um)
		self.idx_cter_astig_ang    =  8 # astigmatism angle
		self.idx_cter_sd_def       =  9 # std dev of defocus (um)
		self.idx_cter_sd_astig_amp = 10 # std dev of ast amp (A)
		self.idx_cter_sd_astig_ang = 11 # std dev of ast angle
		self.idx_cter_error_def    = 12 # frequency at which signal drops by 50% due to estimated error of defocus alone (1/A)
		self.idx_cter_error_astig  = 13 # frequency at which signal drops by 50% due to estimated error of defocus and astigmatism (1/A)
		self.idx_cter_mic_name     = 14 # Micrograph name
		self.idx_cter_box_size     = 15 # <extra> frequency at which CTF oscillation exceeds Nyquist of micrograph (1/A)
		self.idx_cter_error_ctf    = 16 # <extra> frequency at which CTF oscillation exceeds Nyquist of micrograph (1/A)
		self.n_idx_cter            = 17
		self.n_idx_cter_extra      = 3
		
		self.cter_box_size = 512

		# Hold the widget for parameter values of current entry
		self.value_map_list = [None] * self.n_idx_cter
		
		self.idx_rotinf_id             = 0 # line number
		self.idx_rotinf_freq           = 1 # spatial frequency (1/A)
		self.idx_rotinf_exp_no_astig   = 2 # experimental rotational average (no astigmatism)
		self.idx_rotinf_fit_no_astig   = 3 # fitted rotational average (no astigmatism)
		self.idx_rotinf_exp_with_astig = 4 # experimental rotational average (with astigmatism)
		self.idx_rotinf_fit_with_astig = 5 # fitted rotational average (with astigmatism)
		self.n_idx_rotinf              = 6
		
		self.idx_sort_id           =  0
		self.idx_sort_def          =  1
		self.idx_sort_astig_amp    =  2
		self.idx_sort_astig_ang    =  3
		self.idx_sort_sd_def       =  4
		self.idx_sort_sd_astig_amp =  5
		self.idx_sort_sd_astig_ang =  6
		self.idx_sort_error_def    =  7
		self.idx_sort_error_astig  =  8
		self.idx_sort_error_ctf    =  9
		self.n_idx_sort            = 10
		
		self.idx_sort_item_param_label =  0 
		self.idx_sort_item_cter_id     =  1 
		self.n_idx_sort_item           =  2 
		
		# Map idx_sort to idx_cter for sorting item
		self.sort_map_list = [None] * self.n_idx_sort
		self.sort_map_list[self.idx_sort_id]           = ["CTER ID", self.idx_cter_id]
		self.sort_map_list[self.idx_sort_def]          = ["Defocus", self.idx_cter_def]
		self.sort_map_list[self.idx_sort_astig_amp]    = ["Astig. Amp.", self.idx_cter_astig_amp]
		self.sort_map_list[self.idx_sort_astig_ang]    = ["Astig. Ang.", self.idx_cter_astig_ang]
		self.sort_map_list[self.idx_sort_sd_def]       = ["Defocus SD", self.idx_cter_sd_def]
		self.sort_map_list[self.idx_sort_sd_astig_amp] = ["Astig. Amp. SD", self.idx_cter_sd_astig_amp]
		self.sort_map_list[self.idx_sort_sd_astig_ang] = ["Astig. Ang. SD", self.idx_cter_sd_astig_ang]
		self.sort_map_list[self.idx_sort_error_def]    = ["Defocus Error", self.idx_cter_error_def]
		self.sort_map_list[self.idx_sort_error_astig]  = ["Astig. Error", self.idx_cter_error_astig]
		self.sort_map_list[self.idx_sort_error_ctf]    = ["CTF Error", self.idx_cter_error_ctf]
		
		self.idx_hist_def          =  0
		self.idx_hist_astig_amp    =  1
		self.idx_hist_astig_ang    =  2
		self.idx_hist_sd_def       =  3
		self.idx_hist_sd_astig_amp =  4
		self.idx_hist_sd_astig_ang =  5
		self.idx_hist_error_def    =  6
		self.idx_hist_error_astig  =  7
		self.idx_hist_error_ctf    =  8
		self.n_idx_hist            =  9
		
		self.idx_hist_item_param_label =  0
		self.idx_hist_item_idx_cter     =  1 
		self.idx_hist_item_widget      =  2
		self.n_idx_hist_item           =  3
		
		# Map idx_hist to idx_cter for histogram item and slider widget for threshold setting
		self.hist_map_list = [None] * self.n_idx_hist
		self.hist_map_list[self.idx_hist_def]          = ["Defocus", self.idx_cter_def, None]
		self.hist_map_list[self.idx_hist_astig_amp]    = ["Astig. Amp.", self.idx_cter_astig_amp, None]
		self.hist_map_list[self.idx_hist_astig_ang]    = ["Astig. Ang.", self.idx_cter_astig_ang, None]
		self.hist_map_list[self.idx_hist_sd_def]       = ["Defocus SD", self.idx_cter_sd_def, None]
		self.hist_map_list[self.idx_hist_sd_astig_amp] = ["Astig. Amp. SD", self.idx_cter_sd_astig_amp, None]
		self.hist_map_list[self.idx_hist_sd_astig_ang] = ["Astig. Ang. SD", self.idx_cter_sd_astig_ang, None]
		self.hist_map_list[self.idx_hist_error_def]    = ["Defocus Error", self.idx_cter_error_def, None]
		self.hist_map_list[self.idx_hist_error_astig]  = ["Astig. Error", self.idx_cter_error_astig, None]
		self.hist_map_list[self.idx_hist_error_ctf]    = ["CTF Error", self.idx_cter_error_ctf, None]
		
		self.idx_graph_exp_no_astig   =  0
		self.idx_graph_fit_no_astig   =  1
		self.idx_graph_exp_with_astig =  2
		self.idx_graph_fit_with_astig =  3
		self.n_idx_graph              =  4

		self.idx_graph_item_name   =  0
		self.idx_graph_item_label  =  1
		self.idx_graph_idx_rotinf  =  2
		self.idx_graph_item_widget =  3
		self.n_idx_graph_item      =  4
		
		# Map check box widget for graph display setting
		self.graph_list = [None] * self.n_idx_graph
		self.graph_list[self.idx_graph_exp_no_astig]   = ["exp_no_astig", "Exp. No Astig (Black)", self.idx_rotinf_exp_no_astig, None]
		self.graph_list[self.idx_graph_fit_no_astig]   = ["fit_no_astig", "Fit. No Astig (Blue)", self.idx_rotinf_fit_no_astig, None]
		self.graph_list[self.idx_graph_exp_with_astig] = ["exp_with_astig", "Exp. with Astig (Red)", self.idx_rotinf_exp_with_astig, None]
		self.graph_list[self.idx_graph_fit_with_astig] = ["fit_with_astig", "Fit. with Astig (Greed)", self.idx_rotinf_fit_with_astig, None]
		
		self.cter_partres_file_path = None
		self.cter_entry_list        = None
		self.cter_mic_file_path     = None
		self.mic_img                = None
		self.cter_pwrot_file_path   = None
		
		self.curhist=0
		self.curentryperbin=10
		self.curentry=None
		self.cursort=0
		self.cursortoder=False
		
		self.wimage=EMImage2DWidget()
		self.wimage.setWindowTitle("sxgui_cter - Micrograph")

#		self.wfft=EMImage2DWidget()
#		self.wfft.setWindowTitle("sxgui_cter - 2D FFT")

		self.whist=SXPlot2DWidget()
		self.whist.setWindowTitle("sxgui_cter - Histogram")

		# self.wplot=EMPlot2DWidget()
		self.wplot=SXPlot2DWidget()
		self.wplot.setWindowTitle("sxgui_cter - Plot")

#		self.wimage.connect(self.wimage,QtCore.SIGNAL("mousedown"),self.imgmousedown)
#		self.wimage.connect(self.wimage,QtCore.SIGNAL("mousedrag"),self.imgmousedrag)
#		self.wimage.connect(self.wimage,QtCore.SIGNAL("mouseup")  ,self.imgmouseup)
#		self.wfft.connect(self.wfft,QtCore.SIGNAL("mousedown"),self.fftmousedown)
#		self.wfft.connect(self.wfft,QtCore.SIGNAL("mousedrag"),self.fftmousedrag)
#		self.wfft.connect(self.wfft,QtCore.SIGNAL("mouseup")  ,self.fftmouseup)
#		self.wplot.connect(self.wplot,QtCore.SIGNAL("mousedown"),self.plotmousedown)

		self.wimage.mmode="app"
#		self.wfft.mmode="app"

		# This object is itself a widget we need to set up
		self.gbl = QtGui.QGridLayout(self)
		self.gbl.setMargin(8)
		self.gbl.setSpacing(6)

		# --------------------------------------------------------------------------------
		# Columns 1-2
		# --------------------------------------------------------------------------------
		grid_col=0
		grid_row=0

		checkbox_label_width=160

		self.pbopencter=QtGui.QPushButton("Open CTER CTF file")
		self.gbl.addWidget(self.pbopencter,grid_row,grid_col,1,1)
		grid_row += 1
		
		# Make space
		grid_row+=1

		self.vbnentry=ValBox(self,(0,10000),"Num. of entries:",0,90)
		self.vbnentry.setEnabled(False)
		self.vbnentry.intonly=True
		self.gbl.addWidget(self.vbnentry,grid_row,grid_col)
		grid_row+=1

		self.value_map_list[self.idx_cter_vol]=ValBox(self,(0,500),"Voltage (kV):",0.0,90)
		self.value_map_list[self.idx_cter_vol].setEnabled(False)
		self.gbl.addWidget(self.value_map_list[self.idx_cter_vol],grid_row,grid_col)
		grid_row+=1

		self.value_map_list[self.idx_cter_cs]=ValBox(self,(0,5),"Cs (mm):",0.0,90)
		self.value_map_list[self.idx_cter_cs].setEnabled(False)
		self.gbl.addWidget(self.value_map_list[self.idx_cter_cs],grid_row,grid_col)
		grid_row+=1

		self.value_map_list[self.idx_cter_ac]=ValBox(self,(0,100),"Amp. Contrast:",0.0,90)
		self.value_map_list[self.idx_cter_ac].setEnabled(False)
		self.gbl.addWidget(self.value_map_list[self.idx_cter_ac],grid_row,grid_col)
		grid_row+=1

		self.value_map_list[self.idx_cter_apix]=ValBox(self,(0,500),"A/pix:",0.0,90)
		self.value_map_list[self.idx_cter_apix].setEnabled(False)
		self.gbl.addWidget(self.value_map_list[self.idx_cter_apix],grid_row,grid_col)
		grid_row+=1

		self.value_map_list[self.idx_cter_box_size]=ValBox(self,(128,1024),"Window Size:",0,90)
		self.value_map_list[self.idx_cter_box_size].setEnabled(False)
		self.value_map_list[self.idx_cter_box_size].intonly=True
		self.gbl.addWidget(self.value_map_list[self.idx_cter_box_size],grid_row,grid_col)
		grid_row+=1
		
		# Make space
		grid_row+=2

		self.lsort=QtGui.QLabel("Select Curves:",self)
		self.gbl.addWidget(self.lsort,grid_row,grid_col)
		grid_row += 1

		for index_graph in xrange(self.n_idx_graph):
			self.graph_list[index_graph][self.idx_graph_item_widget]=CheckBox(None,self.graph_list[index_graph][self.idx_graph_item_label],True,checkbox_label_width)
			self.gbl.addWidget(self.graph_list[index_graph][self.idx_graph_item_widget],grid_row,grid_col)
			grid_row += 1
		
		# Make space
		grid_row+=1
		
		self.pbrefreshgraphs=QtGui.QPushButton("Refresh Graphs")
		self.gbl.addWidget(self.pbrefreshgraphs,grid_row,grid_col,1,1)
		grid_row += 1

		# --------------------------------------------------------------------------------
		# Columns 3-4
		# --------------------------------------------------------------------------------
		grid_col+=2
		grid_row=0
		
		# plot list and plot mode combobox
		row_span = 17
		col_span = 2
		self.lbentry=SXListWidget(self) # self.lbentry=e2ctf.MyListWidget(self)
		self.lbentry.setSizePolicy(QtGui.QSizePolicy.Preferred,QtGui.QSizePolicy.Expanding)
		self.gbl.addWidget(self.lbentry,grid_row,grid_col,row_span,col_span)
		grid_row += row_span
		
		# --------------------------------------------------------------------------------
		# Columns 5-6 (for Micrograph/CTER Entry) and 7-8 (for Histogram)
		# --------------------------------------------------------------------------------
		grid_col+=col_span
		grid_row=0

		grid_col_next = grid_col + 2
		
		self.ssortedid=ValBox(self,(0,10000),"Sorted ID:",0,90)
		self.ssortedid.setEnabled(False)
		self.ssortedid.intonly=True
		self.gbl.addWidget(self.ssortedid,grid_row,grid_col)
		grid_row+=1
		
		self.value_map_list[self.idx_cter_id]=ValBox(self,(0,10000),"CTER ID:",0,90)
		self.value_map_list[self.idx_cter_id].setEnabled(False)
		self.value_map_list[self.idx_cter_id].intonly=True
		self.gbl.addWidget(self.value_map_list[self.idx_cter_id],grid_row,grid_col)
		grid_row+=1

		self.value_map_list[self.idx_cter_def]=ValBox(self,(0,5),"Defocus:",0.0,90)
		self.value_map_list[self.idx_cter_def].setEnabled(False)
		self.gbl.addWidget(self.value_map_list[self.idx_cter_def],grid_row,grid_col)
		self.hist_map_list[self.idx_hist_def][self.idx_hist_item_widget]=ValSlider(self,(0,5),None,0.0,90)
		self.hist_map_list[self.idx_hist_def][self.idx_hist_item_widget].setEnabled(True)
		self.gbl.addWidget(self.hist_map_list[self.idx_hist_def][self.idx_hist_item_widget],grid_row,grid_col_next,1,2)
		grid_row+=1

		self.value_map_list[self.idx_cter_astig_amp]=ValBox(self,(0,1),"Astig. Amp.:",0.0,90)
		self.value_map_list[self.idx_cter_astig_amp].setEnabled(False)
		self.gbl.addWidget(self.value_map_list[self.idx_cter_astig_amp],grid_row,grid_col)
		self.hist_map_list[self.idx_hist_astig_amp][self.idx_hist_item_widget]=ValSlider(self,(0,1),None,0.0,90)
		self.hist_map_list[self.idx_hist_astig_amp][self.idx_hist_item_widget].setEnabled(False)
		self.gbl.addWidget(self.hist_map_list[self.idx_hist_astig_amp][self.idx_hist_item_widget],grid_row,grid_col_next,1,2)
		grid_row+=1

		self.value_map_list[self.idx_cter_astig_ang]=ValBox(self,(0,180),"Astig. Ang.:",0.0,90)
		self.value_map_list[self.idx_cter_astig_ang].setEnabled(False)
		self.gbl.addWidget(self.value_map_list[self.idx_cter_astig_ang],grid_row,grid_col)
		self.hist_map_list[self.idx_hist_astig_ang][self.idx_hist_item_widget]=ValSlider(self,(0,180),None,0.0,90)
		self.hist_map_list[self.idx_hist_astig_ang][self.idx_hist_item_widget].setEnabled(False)
		self.gbl.addWidget(self.hist_map_list[self.idx_hist_astig_ang][self.idx_hist_item_widget],grid_row,grid_col_next,1,2)
		grid_row+=1

		self.value_map_list[self.idx_cter_bfactor]=ValBox(self,(0,1600),"B factor:",0.0,90)
		self.value_map_list[self.idx_cter_bfactor].setEnabled(False)
		self.gbl.addWidget(self.value_map_list[self.idx_cter_bfactor],grid_row,grid_col)
		grid_row+=1

		self.value_map_list[self.idx_cter_sd_def]=ValBox(self,(0,5),"Defocus SD:",0.0,90)
		self.value_map_list[self.idx_cter_sd_def].setEnabled(False)
		self.gbl.addWidget(self.value_map_list[self.idx_cter_sd_def],grid_row,grid_col)
		self.hist_map_list[self.idx_hist_sd_def][self.idx_hist_item_widget]=ValSlider(self,(0,5),None,0.0,90)
		self.hist_map_list[self.idx_hist_sd_def][self.idx_hist_item_widget].setEnabled(False)
		self.gbl.addWidget(self.hist_map_list[self.idx_hist_sd_def][self.idx_hist_item_widget],grid_row,grid_col_next,1,2)
		grid_row+=1

		self.value_map_list[self.idx_cter_sd_astig_amp]=ValBox(self,(0,1),"Astig. Amp. SD:",0.0,90)
		self.value_map_list[self.idx_cter_sd_astig_amp].setEnabled(False)
		self.gbl.addWidget(self.value_map_list[self.idx_cter_sd_astig_amp],grid_row,grid_col)
		self.hist_map_list[self.idx_hist_sd_astig_amp][self.idx_hist_item_widget]=ValSlider(self,(0,1),None,0.0,90)
		self.hist_map_list[self.idx_hist_sd_astig_amp][self.idx_hist_item_widget].setEnabled(False)
		self.gbl.addWidget(self.hist_map_list[self.idx_hist_sd_astig_amp][self.idx_hist_item_widget],grid_row,grid_col_next,1,2)
		grid_row+=1

		self.value_map_list[self.idx_cter_sd_astig_ang]=ValBox(self,(0,180),"Astig. Ang. SD:",0.0,90)
		self.value_map_list[self.idx_cter_sd_astig_ang].setEnabled(False)
		self.gbl.addWidget(self.value_map_list[self.idx_cter_sd_astig_ang],grid_row,grid_col)
		self.hist_map_list[self.idx_hist_sd_astig_ang][self.idx_hist_item_widget]=ValSlider(self,(0,180),None,0.0,90)
		self.hist_map_list[self.idx_hist_sd_astig_ang][self.idx_hist_item_widget].setEnabled(False)
		self.gbl.addWidget(self.hist_map_list[self.idx_hist_sd_astig_ang][self.idx_hist_item_widget],grid_row,grid_col_next,1,2)
		grid_row+=1

		self.value_map_list[self.idx_cter_error_def]=ValBox(self,(0,10),"Defocus Error:",0,90)
		self.value_map_list[self.idx_cter_error_def].setEnabled(False)
		self.gbl.addWidget(self.value_map_list[self.idx_cter_error_def],grid_row,grid_col)
		self.hist_map_list[self.idx_hist_error_def][self.idx_hist_item_widget]=ValSlider(self,(0,10),None,0,90)
		self.hist_map_list[self.idx_hist_error_def][self.idx_hist_item_widget].setEnabled(False)
		self.gbl.addWidget(self.hist_map_list[self.idx_hist_error_def][self.idx_hist_item_widget],grid_row,grid_col_next,1,2)
		grid_row+=1

		self.value_map_list[self.idx_cter_error_astig]=ValBox(self,(0,10),"Astig. Error:",0,90)
		self.value_map_list[self.idx_cter_error_astig].setEnabled(False)
		self.gbl.addWidget(self.value_map_list[self.idx_cter_error_astig],grid_row,grid_col)
		self.hist_map_list[self.idx_hist_error_astig][self.idx_hist_item_widget]=ValSlider(self,(0,10),None,0,90)
		self.hist_map_list[self.idx_hist_error_astig][self.idx_hist_item_widget].setEnabled(False)
		self.gbl.addWidget(self.hist_map_list[self.idx_hist_error_astig][self.idx_hist_item_widget],grid_row,grid_col_next,1,2)
		grid_row+=1

		self.value_map_list[self.idx_cter_error_ctf]=ValBox(self,(0,10),"CTF Error:",0,90)
		self.value_map_list[self.idx_cter_error_ctf].setEnabled(False)
		self.gbl.addWidget(self.value_map_list[self.idx_cter_error_ctf],grid_row,grid_col)
		self.hist_map_list[self.idx_hist_error_ctf][self.idx_hist_item_widget]=ValSlider(self,(0,10),None,0,90)
		self.hist_map_list[self.idx_hist_error_ctf][self.idx_hist_item_widget].setEnabled(False)
		self.gbl.addWidget(self.hist_map_list[self.idx_hist_error_ctf][self.idx_hist_item_widget],grid_row,grid_col_next,1,2)
		grid_row+=1
		
		# make space
		grid_row += 2

		self.lsort=QtGui.QLabel("Sort CTER Entries:",self)
		self.gbl.addWidget(self.lsort,grid_row,grid_col)
		# grid_row += 1
		
		self.lhist=QtGui.QLabel("Select Histogram:",self)
		self.gbl.addWidget(self.lhist,grid_row,grid_col_next)
		grid_row += 1
		
		self.ssort=QtGui.QComboBox(self)
		for map_entry in self.sort_map_list:
			self.ssort.addItem(map_entry[self.idx_sort_item_param_label])
		self.ssort.setCurrentIndex(self.cursort)
		self.gbl.addWidget(self.ssort,grid_row,grid_col,1,2)
		# grid_row += 1
		
		self.shist=QtGui.QComboBox(self)
		for map_entry in self.hist_map_list:
			self.shist.addItem(map_entry[self.idx_hist_item_param_label])
		self.shist.setCurrentIndex(self.curhist)
		self.gbl.addWidget(self.shist,grid_row,grid_col_next,1,2)
		grid_row += 1

		self.cbsortoder=CheckBox(None,"Decending",self.cursortoder)
		self.gbl.addWidget(self.cbsortoder,grid_row,grid_col)

		self.vsentryperbin=ValSlider(self,(0,10000),"Avg. counts per bin",self.curentryperbin,90)
		self.vsentryperbin.setIntonly(True)
		self.gbl.addWidget(self.vsentryperbin,grid_row,grid_col_next,1,1)
		grid_row += 1

		# this is just a spacer
		self.gbl.setRowStretch(grid_row,2)
		self.gbl.setColumnStretch(3,2)

		# --------------------------------------------------------------------------------
		# Set signal handler
		# --------------------------------------------------------------------------------
		QtCore.QObject.connect(self.pbopencter, QtCore.SIGNAL("clicked(bool)"),self.openCter)

		QtCore.QObject.connect(self.pbrefreshgraphs, QtCore.SIGNAL("clicked(bool)"),self.refreshGraphs)

		QtCore.QObject.connect(self.lbentry,QtCore.SIGNAL("currentRowChanged(int)"),self.newEntry)
#		QtCore.QObject.connect(self.lbentry,QtCore.SIGNAL("keypress"),self.entryKey)
		
		QtCore.QObject.connect(self.ssort,QtCore.SIGNAL("currentIndexChanged(int)"),self.newSort)
		QtCore.QObject.connect(self.cbsortoder, QtCore.SIGNAL("valueChanged"),self.newSortOrder)

		QtCore.QObject.connect(self.shist,QtCore.SIGNAL("currentIndexChanged(int)"),self.newHist)
		
		for idx_graph in xrange(self.n_idx_graph):
			QtCore.QObject.connect(self.graph_list[idx_graph][self.idx_graph_item_widget], QtCore.SIGNAL("valueChanged"),self.updatePlotVisibility)
		
		for idx_hist_def in xrange(self.n_idx_hist):
			QtCore.QObject.connect(self.hist_map_list[idx_hist_def][self.idx_hist_item_widget],QtCore.SIGNAL("valueChanged"),self.updateHist)

		QtCore.QObject.connect(self.vsentryperbin, QtCore.SIGNAL("valueChanged"),self.updateEntryPerBin)
		
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
		self.whist.qt_parent.resize(win_width,win_height)
		self.whist.qt_parent.move(win_left,win_top)
		# Top Right
		win_left = graph_win_width
		win_width = main_win_width;
		self.resize(win_width,win_height)
		self.move(win_left,win_top)
		# Bottom Left
		win_top = win_height + win_height_margin; 
		win_left = 0
		win_width = graph_win_width
		self.wplot.qt_parent.resize(win_width,win_height)
		self.wplot.qt_parent.move(win_left,win_top)
		# Bottom Right
		# Set the image window
		win_left = graph_win_width
		win_width = img_win_width
		img_size = 4096
		scale_factor = float(win_width)/img_size
		self.wimage.set_data(model_blank(img_size,img_size, bckg=1.0)) # resize does not work if no image is set
		self.wimage.qt_parent.resize(win_width,win_height)
		self.wimage.qt_parent.move(win_left,win_top)
		self.wimage.scroll_to(-1 * img_size,-1 * img_size)
		self.wimage.set_scale(scale_factor)
		
		# Try to recover sizes & positions of windows of the previous GUI session
		E2loadappwin("sxgui_cter","main",self)
		E2loadappwin("sxgui_cter","image",self.wimage.qt_parent)
#		E2loadappwin("sxgui_cter","fft",self.wfft.qt_parent)
		E2loadappwin("sxgui_cter","hist",self.whist.qt_parent)
		E2loadappwin("sxgui_cter","plot",self.wplot.qt_parent)


#		if self.cter_entry_list:
# #			self.wfft.show()
#			self.whist.show()
#			self.wplot.show()

		### This section is responsible for background updates
		self.busy=False
#		self.needupdate=True
		self.needredisp=False
#		self.procthread=None
#		self.errors=None # used to communicate errors back from the reprocessing thread

		self.timer=QTimer()
		QtCore.QObject.connect(self.timer, QtCore.SIGNAL("timeout()"), self.timeOut)
		self.timer.start(100)



# 	def entryKey(self,event):
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

	def closeEvent(self,event):
		E2saveappwin("sxgui_cter","main",self)
		if self.cter_entry_list != None:
			E2saveappwin("sxgui_cter","image",self.wimage.qt_parent)
#			E2saveappwin("sxgui_cter","fft",self.wfft.qt_parent)
			E2saveappwin("sxgui_cter","hist",self.whist.qt_parent)
			E2saveappwin("sxgui_cter","plot",self.wplot.qt_parent)
		
		event.accept()
		QtGui.qApp.exit(0)

	def updatePlotVisibility(self,val=None):
		if self.wplot == None: return # it's closed/not visible
		if self.cter_pwrot_file_path == None: return # no cter entry is selected
		
		for idx_graph  in xrange(self.n_idx_graph):
			item_widget = self.graph_list[idx_graph][self.idx_graph_item_widget]
			name = self.graph_list[idx_graph][self.idx_graph_item_name]
			if self.wplot.visibility[name] != item_widget.getValue():
				self.wplot.visibility[name] = item_widget.getValue()
				self.wplot.full_refresh()
				self.wplot.updateGL()

	def updateEntryPerBin(self,curentryperbin):
		if self.whist == None: return # it's closed/not visible
		if self.cter_partres_file_path == None: return # no cter ctf file is selected
		if self.cter_entry_list == None: return # no cter ctf file is selected
		
		if self.curentryperbin == curentryperbin:
			return
		assert self.curentryperbin != curentryperbin
		
		# now set the new entry per bin
		self.curentryperbin = curentryperbin
		
		self.needredisp = True
	
	def updateHist(self):
		if self.whist == None: return # it's closed/not visible
		if self.cter_partres_file_path == None: return # no cter ctf file is selected
		if self.cter_entry_list == None: return # no cter ctf file is selected
		
		val_list = []
		idx_cter = self.hist_map_list[self.curhist][self.idx_hist_item_idx_cter]
		for cter_entry in self.cter_entry_list:
			val_list.append(cter_entry[idx_cter])
		
		n_bin = 10
		assert self.vsentryperbin.getValue() != 0
		n_bin = len(self.cter_entry_list)/ self.curentryperbin
		assert len(val_list) >= n_bin
		from statistics import hist_list
		hist_x_list, hist_y_list = hist_list(val_list, n_bin)
		
		# Pad with zero for better visual impression...
		hist_x_list += [max(val_list)]
		hist_y_list += [0]
		self.whist.set_data((hist_x_list,hist_y_list),"hist_param",quiet=False,color=0,linetype=0,symtype=0)
		x_margin = (hist_x_list[-1] - hist_x_list[0]) * 0.05
		self.whist.rescale(min(val_list),max(val_list),0,max(hist_y_list) * 1.05)
		
		# MRK_NOTE: 2015/12/17 Toshio Moriya
		# This is not good place to update the following information...
		cter_id = self.hist_map_list[self.curhist][self.idx_hist_item_idx_cter]
		param_label = self.hist_map_list[self.curhist][self.idx_hist_item_param_label]
		param_val = self.cter_entry_list[self.curentry][cter_id]
		shape_name = "param_value"
		
		# self.whist.del_shapes((shape_name))
		
		scr_x, scr_y = self.whist.plot2scr(param_val, 0.0)
		self.whist.add_shape(shape_name,EMShape(("scrline",0,1,0,scr_x,self.whist.scrlim[1],scr_x,self.whist.scrlim[1]+self.whist.scrlim[3],1)))
		self.whist.add_shape("%s_val"%(shape_name),EMShape(("scrlabel",0,0,0,self.whist.scrlim[0]+30,self.whist.scrlim[1]+self.whist.scrlim[3]-18,"%s %1.5g"%(param_label,param_val),120.0,-1)))
		
		shape_name = "threshold"
		threshold_label = "Threshold"
		threshold_val = self.hist_map_list[self.curhist][self.idx_hist_item_widget].getValue()
		scr_x, scr_y = self.whist.plot2scr(threshold_val, 0.0)
		self.whist.add_shape(shape_name,EMShape(("scrline",1,0,0,scr_x,self.whist.scrlim[1],scr_x,self.whist.scrlim[1]+self.whist.scrlim[3],1)))
		self.whist.add_shape("%s_val"%(shape_name),EMShape(("scrlabel",0,0,0,self.whist.scrlim[0]+30,self.whist.scrlim[1]+self.whist.scrlim[3]-36,"%s %1.5g"%(threshold_label,threshold_val),120.0,-1)))
		
		self.whist.setAxisParms(param_label,"Counts")
	
	def updatePlot(self):
		if self.wplot == None: return # it's closed/not visible
		if self.cter_pwrot_file_path == None: return # no cter entry is selected
		
		# Now update the plots
		if not os.path.exists(self.cter_pwrot_file_path):
			QtGui.QMessageBox.warning(None,"Warning","Can not find file cter_pwrot_file_path (%s). Please check the contents of pwrot directory" % (self.cter_pwrot_file_path))
			return
			
		self.rotinf_table = read_text_file(self.cter_pwrot_file_path, ncol=-1)
		
		# print "MRK_DEBUG: Last entry of the 1st colum should be a micrograph name %s which is same as " % os.path.basename(self.rotinf_table[0][-1])
		
		mic_basename_rotinf = os.path.basename(self.rotinf_table[0][-1]) # last entry of 1st colum should be associated micrograph
		mic_basename_partres = os.path.basename(self.cter_entry_list[self.curentry][self.idx_cter_mic_name])
		
		if mic_basename_rotinf != mic_basename_partres:
			QtGui.QMessageBox.warning(None,"Warning","Micrograph name (%s) in %s is not same as name (%s) in %s " % (mic_basename_rotinf, os.path.basename(self.cter_pwrot_file_path), mic_basename_partres, os.path.basename(self.cter_partres_file_path)))
			return
		
		for idx_graph in xrange(self.n_idx_graph):
			self.wplot.set_data((self.rotinf_table[self.idx_rotinf_freq],self.rotinf_table[self.graph_list[idx_graph][self.idx_graph_idx_rotinf]]),self.graph_list[idx_graph][self.idx_graph_item_name],quiet=False,color=idx_graph)
		
		self.wplot.rescale(self.rotinf_table[self.idx_rotinf_freq][0],self.rotinf_table[self.idx_rotinf_freq][-1],0.0,1.0)
		
		nyquist_freq = self.rotinf_table[self.idx_rotinf_freq][-1]
		# print "MRK_DEBUG: nyquist_freq = %1.5g" % nyquist_freq
		
		error_name = "error_astig"
		error_label = "Astig. Limit"
		error_freq = self.cter_entry_list[self.curentry][self.idx_cter_error_astig]
		# print "MRK_DEBUG: %s= %1.5g" % (error_name, error_freq)
		if error_freq > 0.0 and error_freq <= nyquist_freq:
			error_scr_x, error_scr_y = self.wplot.plot2scr(error_freq, 0.0)
			self.wplot.add_shape(error_name,EMShape(("scrline",0,0,0.5,error_scr_x,self.wplot.scrlim[1],error_scr_x,self.wplot.scrlim[1]+self.wplot.scrlim[3],1)))
			self.wplot.add_shape("%s_freq"%(error_name),EMShape(("scrlabel",0,0,0,error_scr_x-260,self.wplot.scrlim[1]+self.wplot.scrlim[3]-18,"%s %1.5g (%1.5g)"%(error_label,error_freq,1.0/error_freq),120.0,-1)))
		
		error_name = "error_def"
		error_label = "Defocus Limit"
		error_freq = self.cter_entry_list[self.curentry][self.idx_cter_error_def]
		# print "MRK_DEBUG: %s= %1.5g" % (error_name, error_freq)
		if error_freq > 0.0 and error_freq <= nyquist_freq:
			error_scr_x, error_scr_y = self.wplot.plot2scr(error_freq, 0.0)
			self.wplot.add_shape(error_name,EMShape(("scrline",0.5,0,0,error_scr_x,self.wplot.scrlim[1],error_scr_x,self.wplot.scrlim[1]+self.wplot.scrlim[3],1)))
			self.wplot.add_shape("%s_freq"%(error_name),EMShape(("scrlabel",0,0,0,error_scr_x-260,self.wplot.scrlim[1]+self.wplot.scrlim[3]-36,"%s %1.5g (%1.5g)"%(error_label,error_freq,1.0/error_freq),120.0,-1)))
		
		error_name = "error_ctf"
		error_label = "CTF Limit"
		error_freq = self.cter_entry_list[self.curentry][self.idx_cter_error_ctf]
		# print "MRK_DEBUG: %s= %1.5g" % (error_name, error_freq)
		if error_freq > 0.0 and error_freq <= nyquist_freq:
			error_scr_x, error_scr_y = self.wplot.plot2scr(error_freq, 0.0)
			self.wplot.add_shape(error_name,EMShape(("scrline",0,0.5,0,error_scr_x,self.wplot.scrlim[1],error_scr_x,self.wplot.scrlim[1]+self.wplot.scrlim[3],1)))
			self.wplot.add_shape("%s_freq"%(error_name),EMShape(("scrlabel",0,0,0,error_scr_x-260,self.wplot.scrlim[1]+self.wplot.scrlim[3]-54,"%s %1.5g (%1.5g)"%(error_label,error_freq,1.0/error_freq),120.0,-1)))

		self.wplot.setAxisParms("frequency (1/"+ "$\AA$" +")","power spectrum")
		# self.wplot.setAxisParms("frequency (1/"+ "$\AA$" +")","power spectrum", "linear", "log")
		
		self.updatePlotVisibility()
	
	def timeOut(self):
		if self.busy : return
		
		# Redisplay before spawning thread for more interactive display
		if self.needredisp :
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
	
	def readCterCtfFile(self, file_path):
		"""Read all entries from a CTER CTF file into the list box"""
		
		if os.path.exists(file_path):
			new_entry_list = read_text_row(file_path)
			if len(new_entry_list) == 0:
				QMessageBox.warning(self, "No entry is found", "The specified CTER CTF file (%S) does not contain any entry. Please try it again." % new_cter_entry_list)
				return
				
			# print "MRK_DEBUG: Detected %s entries in %s" % (len(new_entry_list), file_path)
			# print "MRK_DEBUG: Num. of Columns is %d in %s" % (len(new_entry_list[0]), file_path)
			assert len(new_entry_list[0]) == self.n_idx_cter - self.n_idx_cter_extra, "MRK_DEBUG: The number of columns have to be %d in %s" % (self.n_idx_cter - self.n_idx_cter_extra, file_path)
			
			for i in xrange(len(new_entry_list)):
				new_entry_list[i] = [i] +  new_entry_list[i] + [self.cter_box_size, 0.5] # Add extra items first to make sure indices match
				# assert (options.import_ctf)
				# Cut off frequency components higher than CTF limit 
				cter_box_size = new_entry_list[i][self.idx_cter_box_size]
				cter_def = new_entry_list[i][self.idx_cter_def]
				cter_cs = new_entry_list[i][self.idx_cter_cs]
				cter_vol = new_entry_list[i][self.idx_cter_vol]
				cter_apix = new_entry_list[i][self.idx_cter_apix]
				cter_limit_ab_freq, cter_limit_freq = ctflimit(cter_box_size, cter_def, cter_cs, cter_vol, cter_apix)
				# Limiting_frequency[cycle/A]: xr[cycle/A]. Max is Nyquist frequency = 1.0[cycle]/(2*apix[A/pixel]). <UNIT: [cycle/(A/pixel)/[pixel])] => [(cycle*pixel)/(A*pixel] => [cycle/A]>
				# limiting_period(Angstrom resolution)[A/cycle]: 1.0/xr[cycle/A]. Min is Nyquist period = (2*apix[A/pixel]). <UNIT: [1/(cycle/A)] = [A/cycle]>
				# Width of Fourier pixel [pixel/A]: fwpix = 1.0[pixel]/(2*apix[A/pixel])/box_half[pixel] = 1[pixel]/fullsize[A]) <UNIT: [pixel/(A/pixel)/(pixel)] = [pixel*(pixel/A)*(1/pixel) = [pixel/A]>
				# Limiting_absolute_frequency [cycle/pixel] int(xr/fwpix+0.5) = <Unit:[(cycle/A)/(pixel/A)] = [(cycle/A)*(A/pixel)] = [cycle/pixel]>
				# return  Limiting_abs_frequency [cycle/pixel]Limiting_frequency[cycle/A] <= int(xr/fwpix+0.5),xr
				new_entry_list[i][self.idx_cter_error_ctf] = cter_limit_freq
			
			self.cter_partres_file_path = file_path
			self.cter_entry_list = new_entry_list
			self.sortByParam()
			self.lbentry.clear()
			for cter_entry in self.cter_entry_list:
				self.lbentry.addItem(os.path.basename(cter_entry[self.idx_cter_mic_name]))
			
			self.curentry = None # Always force to execute newEntry
			self.lbentry.setCurrentRow(0) # Should emit signal to call newEntry
			
			# Set the number of entries
			self.vbnentry.setValue(len(self.cter_entry_list))
			
			# Set the range of threshold slider
			for idx_hist in xrange(self.n_idx_hist):
				item_widget = self.hist_map_list[idx_hist][self.idx_hist_item_widget]
				idx_cter = self.hist_map_list[idx_hist][self.idx_hist_item_idx_cter]
				min_val = min(self.cter_entry_list, key=lambda x:x[idx_cter])[idx_cter]
				max_val = max(self.cter_entry_list, key=lambda x:x[idx_cter])[idx_cter]
				item_widget.setRange(min_val,max_val)
				item_widget.setValue(min_val)
				
			# Set the range of histogram bin
			self.vsentryperbin.setRange(1,len(self.cter_entry_list))
			self.vsentryperbin.setValue(self.curentryperbin)
			
	def openCter(self,val=None):
		"""Open CTER CTF file"""
		
		file_path = str(QtGui.QFileDialog.getOpenFileName(self, "Open CTER CTF File", options = QtGui.QFileDialog.DontUseNativeDialog))
		if file_path != "":
			self.readCterCtfFile(file_path)
		# else: # Do nothing
	
	def sortByParam(self):
		"""Do actual sorting of CTER entries based on current setting."""
		
		# sort CTER entry list
		assert (self.cter_entry_list != None)
		self.cter_entry_list = sorted(self.cter_entry_list, key=lambda x: x[self.sort_map_list[self.cursort][self.idx_sort_item_cter_id]], reverse=self.cursortoder)
		self.lbentry.clear()
		for cter_entry in self.cter_entry_list:
			self.lbentry.addItem(os.path.basename(cter_entry[self.idx_cter_mic_name]))
			
		self.curentry = None # Always force to execute newEntry
		self.lbentry.setCurrentRow(0) # Should emit signal to call newEntry
		
	def newSortOrder(self, sortoder):
		"""Change sorting order of CTER entries."""
		if self.cursortoder == sortoder:
			return

		# now set the new item
		self.cursortoder = sortoder
		# print "MRK_DEBUG: Sort oder (%s) is changed when Param entry %d (%s, %d). " % (self.cursortoder, self.cursort, self.sort_map_list[self.cursort][self.idx_sort_item_param_label], self.sort_map_list[self.cursort][self.idx_sort_item_cter_id])
		
		if self.cter_entry_list == None:
			return
		
		self.sortByParam()
		
	def newSort(self,cursort):
		"""Sort CTER entries by selected parameter values."""
		if self.cursort == cursort:
			return

		# now set the new item
		self.cursort = cursort
		# print "MRK_DEBUG: Param entry %d (%s, %d) is selected from combo box with sort oder (%s)" % (self.cursort, self.sort_map_list[self.cursort][self.idx_sort_item_param_label], self.sort_map_list[self.cursort][self.idx_sort_item_cter_id], self.cursortoder)
		
		if self.cter_entry_list == None:
			return
		
		self.sortByParam()
	
	def newHist(self,currow):
		"called when a new row is selected from the Histogram list box"
		
		if self.curhist == currow:
			return
		
		# Disable old item
		self.hist_map_list[self.curhist][self.idx_hist_item_widget].setEnabled(False)
		
		# now set the new item and enalble it
		self.curhist=currow
		# print "MRK_DEBUG: Row No. %d is selected from histogram list box" % (self.curhist)
		
		self.hist_map_list[self.curhist][self.idx_hist_item_widget].setEnabled(True)
		self.whist.setWindowTitle("sxgui_cter - %s Histogram" % (self.hist_map_list[self.curhist][self.idx_hist_item_param_label]))
	
		self.needredisp = True
	
	def refreshGraphs(self,item):
		self.needredisp = True
	
	def newEntry(self,currow):
		"called when a new data set is selected from the CTER Entry list box"
		
		# If list selection does not change, do nothing
		if self.curentry == currow:
			return
		assert self.curentry != currow
		
		new_cter_mic_file_path = self.cter_entry_list[currow][self.idx_cter_mic_name]
		
		if self.cter_mic_file_path == new_cter_mic_file_path:
			return
		assert self.cter_mic_file_path != new_cter_mic_file_path
		
		cter_pwrot_dir = os.path.dirname(self.cter_partres_file_path).replace("partresdir", "pwrot")
		new_cter_pwrot_file_path = os.path.join(cter_pwrot_dir, "rotinf%04d.txt" % self.cter_entry_list[currow][self.idx_cter_id])
		
		if self.cter_pwrot_file_path == new_cter_pwrot_file_path:
			return
		assert self.cter_pwrot_file_path != new_cter_pwrot_file_path
		
		
		# now set the new item
		self.curentry=currow
		self.cter_mic_file_path = new_cter_mic_file_path
		self.cter_pwrot_file_path = new_cter_pwrot_file_path
		# print "MRK_DEBUG: Row No. %d (CTER Entry No. %d) is selected from cter entry list box" % (self.curentry, self.cter_entry_list[self.curentry][self.idx_cter_id])

		# Now update the image
		if not os.path.exists(self.cter_mic_file_path):
			QtGui.QMessageBox.warning(None,"Warning","Can not find micrograph (%s). Please make sure the micrograph names in %s are correct." % (self.cter_mic_file_path, os.path.basename(self.cter_partres_file_path)))
			self.mic_img = EMData() # Set empty image...
		else:
			self.mic_img = EMData(self.cter_mic_file_path) # read the image from disk

		self.wimage.set_data(self.mic_img)
		
		self.ssortedid.setValue(self.curentry,True)
		
		for idx_cter in xrange(self.n_idx_cter):
			if idx_cter != self.idx_cter_mic_name:
				self.value_map_list[idx_cter].setValue(self.cter_entry_list[self.curentry][idx_cter],True)
		
		self.wimage.setWindowTitle("sxgui_cter - Micrograph - %s, %s" % (os.path.basename(self.cter_entry_list[self.curentry][self.idx_cter_mic_name]), os.path.basename(self.cter_pwrot_file_path)))
#		self.wfft.setWindowTitle("sxgui_cter - 2D FFT - "+fsp.split("/")[-1])
		self.wplot.setWindowTitle("sxgui_cter - Plot - %s, %s" % (os.path.basename(self.cter_entry_list[self.curentry][self.idx_cter_mic_name]), os.path.basename(self.cter_pwrot_file_path)))
		
		self.needredisp = True
	
	def redisplay(self):
		self.needredisp=False
		self.busy=True
		
		if self.cter_entry_list != None:
			
			if not self.wimage.isVisible():
				self.wimage.show()
#			if not self.wfft.isVisible():
#				self.wfft.show()
			if not self.wplot.isVisible():
				self.wplot.show()
			if not self.whist.isVisible():
				self.whist.show()
		
#		try:
#			if self.cter_entry_list != None:
#				self.wimage.raise_()
#				self.wfft.raise_()
#				self.wplot.raise_()
#				self.whist.raise_()
#		except: pass

		try:
			self.updatePlot()
		except:
			print "Recieved unexpected exception from updatePlot() in redisplay(): ", sys.exc_info()[0]
			raise
			
		try:
			self.updateHist()
		except:
			print "Recieved exception from updatePlot() in updateHist(): ", sys.exc_info()[0]
			raise
		
		self.busy=False
	
if __name__ == "__main__":
	main()
