#!/usr/bin/env python
#
# Author: Lan Dang, 03/17/2022 (dlan@bcm.edu)
import sys
from past.utils import old_div
import OpenGL
OpenGL.ERROR_CHECKING = False
from OpenGL.GL import *
from OpenGL.GLU import *
from PyQt5 import QtGui, QtWidgets, QtCore, QtOpenGL
from PyQt5.QtWidgets import QSplitter, QHBoxLayout
from PyQt5.QtCore import Qt
from EMAN2 import *
from EMAN2_utils import interp_points
from eman2_gui.emapplication import get_application, EMApp
from eman2_gui.emimage import EMImageWidget
from eman2_gui.emimage2d import EMImage2DWidget
from eman2_gui.emannotate2d import EMAnnotate2DWidget,EMSegTab
from eman2_gui.emimagemx import EMImageMXWidget
from eman2_gui.emscene3d import EMScene3D
from eman2_gui.emdataitem3d import EMDataItem3D, EMIsosurface
from eman2_gui.embrowser import EMBrowserWidget
from eman2_gui.emshape import EMShape
from eman2_gui.valslider import ValSlider,ValBox,StringBox,EMSpinWidget

import scipy.spatial.distance as scipydist
import scipy.ndimage as ndi
import tensorflow as tf
from tensorflow.keras.layers import Input, Conv2D, MaxPooling2D,Flatten,Dense,Dropout,BatchNormalization, concatenate, Conv2DTranspose
from tensorflow.keras.models import Model
from tensorflow.keras.optimizers import Adam, SGD
from tensorflow.keras.callbacks import ModelCheckpoint
from tensorflow.keras import backend as K
from tensorflow.keras.callbacks import History

import numpy as np
from eman2_gui.valslider import ValSlider
import weakref

from matplotlib.patches import Circle
import matplotlib.path as mplPath
import matplotlib.pyplot as plt

def main():
	usage="""annotate tomograms in a folder. annotated images will be saved in ./segs folder.
	 still developing

	[prog] --folder <path_to_tomogram_folder>

	"""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	#parser.add_pos_argument(name="tomogram",help="Specify a tomogram from which you want to extract particles.", default="", guitype='filebox', browser="EMTomoBoxesTable(withmodal=True,multiselect=False)", row=0, col=0,rowspan=1, colspan=2, mode="box3d,box2d")
	parser.add_argument("--folder",type=str, help="List the folder contain all tomograms to process", default="./tomograms/")

	parser.add_argument("--seg_folder",type=str, help="List the folder contain all annotation file", default="./segs/")
	parser.add_argument("--region_sz",type=int, help="Region size for Region I/O. -1 reads whole tomogram", default=680)
	#parser.add_argument("--alltomograms",default=False,help="Process all tomograms from tomograms folder")
	#parser.add_argument("--boxsize","-b",type=int,help="Box size in pixels",default=-1)

	(options, args) = parser.parse_args()

	#
	# if len(args) == 0:
	# 	print("INPUT ERROR: You must specify an image to display.")
	# 	sys.exit(1)


	#img = args[0]

	#imghdr = EMData(img,0,True)
	#options.apix = imghdr['apix_x']
	app = EMApp()
	awin = EMAnnotateWindow(app,options)




	awin.resize(1120,720)
	awin.show()


	x=app.exec_()
	#E2end(logid)
	sys.exit(0)


class EMAnnotateWindow(QtWidgets.QMainWindow):
	keypress = QtCore.pyqtSignal(QtGui.QKeyEvent)
	def __init__(self, application,options,data=None,annotate=None):
		super().__init__()
		self.app = weakref.ref(application)
		self.setMinimumSize(1120, 720)
		self.options=options
		self.setWindowTitle("Image Viewer")
		self.tom_folder = options.folder
		self.seg_folder = options.seg_folder

		try:
			os.mkdir(self.seg_folder)
			print("Directory", self.seg_folder, "is created" )
		except OSError as error:
			print("Directory",self.seg_folder,"already existed. Continue")
			pass

		self.tomogram_list = QtWidgets.QListWidget()
		self.tom_file_list = []
		for file_name in os.listdir(self.tom_folder):
			if file_name.endswith(".hdf"):
				self.tomogram_list.addItem(file_name)
				self.tom_file_list.append(file_name)
			# info=js_open_dict("info/annotate_"+file_name[:-4]+".json")
			# info["class"] = []
			# info["boxes"] = []
			# info["seg"] = ""
			# info.close()

		self.tomogram_list.setCurrentRow(0)
		self.data_file = str(os.path.join(self.tom_folder,self.tomogram_list.item(0).text()))
		hdr=EMData(self.data_file, 0,True)
		self.nx = hdr["nx"]
		self.ny = hdr["ny"]
		self.nz=hdr["nz"]

		#TODO
		#Need to copy header to annotation file
		#Need to call garbage collector to remove trash from buffer memory

		self.seg_path = os.path.join(self.seg_folder,self.tomogram_list.item(0).text()[0:-4]+"_seg.hdf")
		if not os.path.isfile(self.seg_path):
			seg_out = EMData(self.nx,self.ny,self.nz)
			#self.write_header(seg_out)
			seg_out.write_image(self.seg_path)
			del seg_out
		else:
			print("Seg file for the first one already exists. Continue ")
			pass

		#print("Nz ,iz",hdr["nz"],iz)
		# info=js_open_dict("info/annotate_"+self.tomogram_list.item(0).text()+".json")
		# try:
		# 	self.boxes=info["boxes"]
		# except:
		# 	self.boxes = []
		# info.close()
		#self.ann_file = "./segs/"+self.data_file[0:-4]+"_seg.hdf"
		print("File_name", self.data_file)
		self.data = EMData(self.data_file)
		self.apix=self.data['apix_x']

		try:
			self.annotation = EMData(annotate)
		except:
			self.annotation = None
			pass

		self.tomogram_list.currentRowChanged[int].connect(self.tomolist_current_change)


		self.centralWidget = QtWidgets.QWidget()
		self.setCentralWidget(self.centralWidget)
		self.gbl = QtWidgets.QGridLayout(self.centralWidget)



		self.gbl.setColumnStretch(0,1)
		self.gbl.setColumnStretch(1,3)
		self.gbl.setContentsMargins(8, 8, 8, 8)
		self.gbl.setSpacing(10)
		if options.region_sz == -1:
			self.img_view_region_size = self.nx
			print("Reading full images of size", self.nx)
		elif options.region_sz == 0:
			print("Region size needs to be greater than 0. Set region size to default value of 680.")
			self.img_view_region_size = 680
		else:
			self.img_view_region_size = options.region_sz
		self.thumbnail_size=220


		self.img_view = EMAnnotate2DWidget(sizehint=(680,680))
		#self.img_view.setSizePolicy(QtWidgets.QSizePolicy.Fixed,QtWidgets.QSizePolicy.Fixed)
		self.img_view.setMinimumSize(680, 680)

		#self.img_view.show()



		self.boxes = []
		self.unet = None

		#Thumbnail
		#Attribute error occurs when Thumbnail was initialized before the EMAnnotateWindow


		self.thumbnail = Thumbnail(current_file=self.data_file,target=self.img_view,app_target=self,tn_size=self.thumbnail_size)
		self.thumbnail.setSizePolicy(QtWidgets.QSizePolicy.Fixed,QtWidgets.QSizePolicy.Fixed)
		self.thumbnail.setMinimumSize(self.thumbnail_size, self.thumbnail_size)
		try:
			self.thumbnail.set_im()
		except:
			pass

		self.img_view_inspector = self.img_view.get_inspector()
		self.tomo_list_panel=QtWidgets.QWidget()
		tomo_vbl = QtWidgets.QGridLayout()
		tomo_vbl.addWidget(QtWidgets.QLabel("Tomograms"),0,0)
		tomo_vbl.addWidget(self.tomogram_list,1,0)
		tomo_vbl.setRowStretch(1,2)
		tomo_vbl.setRowStretch(2,1)
		tomo_vbl.setRowStretch(3,1)
		tomo_vbl.addWidget(self.thumbnail,2,0)
		zt_hbl = QtWidgets.QHBoxLayout()
		zt_hbl.addWidget(QtWidgets.QLabel("Z-thickness"))
		self.zt_spinbox = QtWidgets.QSpinBox(self)
		self.zt_spinbox.setValue(0)
		self.zt_spinbox.setMinimum(-1)
		self.zt_spinbox.setMaximum(self.nz//2)
		zt_hbl.addWidget(self.zt_spinbox)
		tomo_vbl.addLayout(zt_hbl,3,0)
		self.tomo_list_panel.setLayout(tomo_vbl)
		self.tomo_list_panel.setWindowTitle("Tomograms")
		self.tomo_list_panel.show()


		filter_vbl = QtWidgets.QVBoxLayout()
		self.ft_vbl=QtWidgets.QVBoxLayout()
		self.ft_vbl.addWidget(QtWidgets.QLabel("Filters"))

		filter_vbl.addLayout(self.ft_vbl)
		self.procbox1=StringBox(label="Process1:",value="filter.lowpass.gauss:cutoff_abs=0.125",showenable=0)
		self.procbox1.setFixedWidth(350)
		self.ft_vbl.addWidget(self.procbox1)

		self.procbox2=StringBox(label="Process2:",value="filter.highpass.gauss:cutoff_pixels=3",showenable=0)
		self.procbox2.setFixedWidth(350)
		self.ft_vbl.addWidget(self.procbox2)

		self.procbox3=StringBox(label="Process3:",value="math.linear:scale=5:shift=0",showenable=0)
		self.procbox3.setFixedWidth(350)
		self.ft_vbl.addWidget(self.procbox3)

		self.proclbl1=QtWidgets.QLabel("\t\tImage unchanged, display only!")
		self.ft_vbl.addWidget(self.proclbl1)

		self.procbox1.enableChanged.connect(self.do_filters)
		self.procbox1.textChanged.connect(self.do_filters)
		self.procbox2.enableChanged.connect(self.do_filters)
		self.procbox2.textChanged.connect(self.do_filters)
		self.procbox3.enableChanged.connect(self.do_filters)
		self.procbox3.textChanged.connect(self.do_filters)


		#self.thumbnail.setMinimumSize(300,300)


		#Basic tools
		basic_label = QtWidgets.QLabel("Basic tools")
		basic_vbl = QtWidgets.QVBoxLayout()
		basic_vbl.addWidget(basic_label)
		basic_button_l = QtWidgets.QVBoxLayout()

		self.brush_tab = QtWidgets.QWidget()
		self.btlay = QtWidgets.QVBoxLayout(self.brush_tab)
		self.btlay.addWidget(QtWidgets.QLabel("Use paint brush to annotate on tomogram.\nUse manual annotate tools panel"))
		self.basic_tab_num = 0
		self.linear_tab = QtWidgets.QWidget()
		self.ltlay = QtWidgets.QGridLayout(self.linear_tab)
		self.ltlay.setColumnStretch(70,70)

		self.points_label = QtWidgets.QLabel("Choose anchor points. \nCtrl+Click: start a contour.\nShift+Click: delete a point.")
		self.clear_button = QtWidgets.QPushButton("Clear points")
		self.interp_button = QtWidgets.QPushButton("Interpolate")
		self.tx_line_width=QtWidgets.QSpinBox(self)
		self.tx_line_width.setValue(5)
		self.tx_ann_class=QtWidgets.QSpinBox(self)
		self.tx_ann_class.setValue(1)

		self.ltlay.addWidget(self.points_label,0,0,2,2)
		self.ltlay.addWidget(QtWidgets.QLabel("Line Width"),2,2,1,2)
		self.ltlay.addWidget(self.tx_line_width,2,3,1,1)
		self.ltlay.addWidget(QtWidgets.QLabel("Ann Class"),3,2,1,2)
		self.ltlay.addWidget(self.tx_ann_class,3,3,1,1)
		self.tx_interp=QtWidgets.QSpinBox(self)
		self.tx_interp.setValue(20)
		self.ltlay.addWidget(self.interp_button,0,2,1,1)
		self.ltlay.addWidget(self.tx_interp,0,3,1,1)
		self.ltlay.addWidget(self.clear_button,1,2,1,1)
		self.ltlay.addWidget(QtWidgets.QLabel("Draw line to annotate file"),2,0,1,2)

		self.contour_tab = QtWidgets.QWidget()
		self.ctlay = QtWidgets.QGridLayout(self.contour_tab)
		self.ctlay.setColumnStretch(70,70)
		self.clear_contour_button = QtWidgets.QPushButton("Clear points")
		self.draw_contour_label = QtWidgets.QLabel("Choose anchor points. \nCtrl+Click: start a contour.\nShift+Click: delete a point.")
		self.fill_contour_checkbox = QtWidgets.QCheckBox("Fill contour")
		self.fill_contour_checkbox.setChecked(True)
		self.ct_line_width=QtWidgets.QSpinBox(self)
		self.ct_line_width.setValue(5)
		self.ct_ann_class=QtWidgets.QSpinBox(self)
		self.ct_ann_class.setValue(1)


		self.ctlay.addWidget(self.draw_contour_label,0,0,2,2)
		self.ctlay.addWidget(QtWidgets.QLabel("Draw line to annotate file"),2,0,1,2)
		self.ctlay.addWidget(QtWidgets.QLabel("Line Width"),2,2,1,1)
		self.ctlay.addWidget(self.ct_line_width,2,3,1,1)
		self.ctlay.addWidget(QtWidgets.QLabel("Ann Class"),3,2,1,1)
		self.ctlay.addWidget(self.ct_ann_class,3,3,1,1)
		self.ctlay.addWidget(self.clear_contour_button,0,2,1,1)
		self.ctlay.addWidget(self.fill_contour_checkbox,4,2,1,1)

		self.boxer_tab = QtWidgets.QWidget()
		self.bxlay = QtWidgets.QGridLayout(self.boxer_tab)
		self.bxlay.setColumnStretch(70,70)
		self.random_bx_bt = QtWidgets.QPushButton("Create Random Boxes")
		self.clear_bx_bt = QtWidgets.QPushButton("Clear Boxes")
		self.extract_bt = QtWidgets.QPushButton("Extract Boxes")


		self.random_bx_sb = StringBox(label="No Box:",value="20",showenable=-1)
		self.bsz_vs = ValSlider(self,value=64,rng=(0,300),rounding=0,label="Box Size")
		self.bsz_vs.setIntonly(1)


		self.bxlay.addWidget(QtWidgets.QLabel("Select regions \nfor training nnet"),0,0,1,1)
		self.bxlay.addWidget(self.bsz_vs,1,0,1,2)
		self.bxlay.addWidget(self.random_bx_bt,2,0,1,1)
		self.bxlay.addWidget(self.random_bx_sb,2,1,1,1)
		self.bxlay.addWidget(self.clear_bx_bt,3,0,1,1)
		self.bxlay.addWidget(self.extract_bt,3,1,1,1)


		self.basic_tab = QtWidgets.QTabWidget()
		# self.test_seg_tab = EMSegTab(target=self.img_view)
		# self.basic_tab.addTab(self.test_seg_tab,"Seg")
		self.basic_tab.addTab(self.brush_tab,"Brush")
		self.basic_tab.addTab(self.linear_tab,"Linear")
		self.basic_tab.addTab(self.contour_tab,"Contour")
		self.basic_tab.addTab(self.boxer_tab,"Boxer")

		basic_button_l.addWidget(self.basic_tab)
		basic_vbl.addLayout(basic_button_l)

		self.basic_tab.tabBarClicked[int].connect(self.basic_tab_change)
		self.basic_tab.currentChanged[int].connect(self.basic_tab_change)
		self.interp_button.clicked[bool].connect(self.interp_bt_clicked)
		self.clear_button.clicked[bool].connect(self.clear_points)
		self.clear_contour_button.clicked[bool].connect(self.clear_points)
		self.fill_contour_checkbox.stateChanged[int].connect(self.fill_contour_checkbox_changed)
		self.random_bx_bt.clicked[bool].connect(self.random_bx_bt_clicked)
		self.clear_bx_bt.clicked[bool].connect(self.clear_bx_bt_clicked)
		self.extract_bt.clicked[bool].connect(self.extract_bt_clicked)
		self.bsz_vs.valueChanged.connect(self.bsz_vs_value_changed)



		self.assisted_tab = QtWidgets.QTabWidget()
		self.nn_tab = NNet_Tab(target=self)
		self.morp_tab = Morp_Tab(target=self)
		self.binary_tab = Binary_Tab(target=self)
		self.spec_tab = Specific_Tab(target=self)
		self.assisted_tab.addTab(self.nn_tab,"NeuralNetwork")
		self.assisted_tab.addTab(self.binary_tab,"AutoDetect")
		self.assisted_tab.addTab(self.morp_tab,"Morphological")
		self.assisted_tab.addTab(self.spec_tab,"Specific")
		self.assisted_tab.currentChanged[int].connect(self.assisted_tab_changed)

		#assisted tab setup + function
		assisted_vbl = QtWidgets.QVBoxLayout()
		assisted_vbl.addWidget(QtWidgets.QLabel("Automated tools"))
		assisted_vbl.addWidget(self.assisted_tab)

		#Right Panel
		self.button_gbl = QtWidgets.QGridLayout()
		self.button_gbl.setColumnStretch(170,170)
		self.button_gbl.addLayout(filter_vbl,0,0,2,1)
		self.button_gbl.addLayout(basic_vbl,2,0,2,1)
		self.button_gbl.addLayout(assisted_vbl,4,0,2,1)



		self.test_button = QtWidgets.QPushButton("Test Button")
		self.button_gbl.addWidget(self.test_button,6,0,1,1)

		inspector_vbl = QtWidgets.QVBoxLayout()
		inspector_vbl.addWidget(QtWidgets.QLabel("Manual Annotate Tools"))
		inspector_vbl.addWidget(self.img_view_inspector)

		#self.gbl.addLayout(tomo_vbl,1,0,1,1)
		self.gbl.addWidget(self.img_view,0,0,1,1)
		self.gbl.addLayout(inspector_vbl,0,2,1,1)
		self.centralWidget.setLayout(self.gbl)

		self.control_panel = QtWidgets.QWidget()
		self.control_panel.setLayout(self.button_gbl)
		self.control_panel.setWindowTitle("Control Panel")
		self.control_panel.show()

		self.test_button.clicked[bool].connect(self.test_drawing_function)
		#TO FIX LATER
		#self.cb_group.buttonClicked[QtWidgets.QAbstractButton].connect(self.on_check_cb_group)
		# self.table_tom.cellClicked[int,int].connect(self.on_table_tom)
		# self.table_tom.currentCellChanged.connect(self.table_tom_cell_changed)

		#Need to fix
		self.lb_lines=QtWidgets.QLabel("")
		self.lb_lines.setWordWrap(True)

		pts=[]
		self.curve=Curve(img=self.img_view, points=pts )
		self.curve_shape_index = 0
		self.img_view.shapes[0]=self.curve

		self.contour=Contour(img=self.img_view, points=pts)
		self.contour_shape_index = 1
		self.img_view.shapes[1]=self.contour
		self.img_view.mouseup.connect(self.img_view_mouse_up)
		self.img_view.keypress.connect(self.key_press)
		#self.img_view.mousewheel.connect(self.img_view_wheel_event)
		self.img_view.mousedrag.connect(self.img_view_mouse_drag)
		self.zt_spinbox.valueChanged.connect(self.zt_spinbox_changed)

		E2loadappwin("e2annotate","main",self)
		E2loadappwin("e2annotate","controlpanel",self.control_panel)
		E2loadappwin("e2annotate","tomograms",self.tomo_list_panel)
		#self.update_label()

		glEnable(GL_POINT_SMOOTH)
		glEnable( GL_LINE_SMOOTH );
		glEnable( GL_POLYGON_SMOOTH );
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	def img_view_wheel_event(self, event):
		return

	def img_view_mouse_drag(self, event):
		if event.buttons()&Qt.RightButton:
			return
		else:
			return

	def get_data(self):
		return self.img_view.get_full_data()
	def get_annotation(self):
		return self.img_view.get_full_annotation()
	def get_inspector(self):
		return self.img_view_inspector

	def zchange(self,value):
		print(value)


	def do_filters(self):
		print("do filter")
		ret=[]
		if self.procbox1.getEnabled():
			try:
				nm,op=parsemodopt(self.procbox1.getValue())
				ret.append(Processors.get(nm,op))
			except:
				print("Error with processor: ",self.procbox1.getValue())

		if self.procbox2.getEnabled():
			try:
				nm,op=parsemodopt(self.procbox2.getValue())
				ret.append(Processors.get(nm,op))
			except:
				print("Error with processor2: ",self.procbox2.getValue())

		if self.procbox3.getEnabled():
			try:
				nm,op=parsemodopt(self.procbox3.getValue())
				ret.append(Processors.get(nm,op))
			except:
				print("Error with processor: ",self.procbox3.getValue())
		#print(ret)
		self.img_view.set_disp_proc(ret)

	def assisted_tab_changed(self):
		self.reset_morp_params()

	def reset_morp_params(self,reset_vs=False):
		# self.closing_thresh =32
		# self.opening_thresh =32
		self.binary_tab.closing_n_iters =1
		self.binary_tab.opening_n_iters =1
		self.morp_tab.morp_n_iters_sp.setValue(1)
		if reset_vs:
			self.binary_tab.bin_low_pass_vs.setValue(1)
			self.binary_tab.bin_threshold_vs.setValue(0.001)


	def tomolist_current_change(self, int):
		# info=js_open_dict("info/annotate_"+self.data_file[:-4]+".json")
		# #info["class"] = self.classes
		# info["boxes"] = self.boxes
		# info.close()
		print("Current item", self.tomogram_list.item(int).text())
		print(str(self.tomogram_list.item(int).text()))
		self.data_file = str(os.path.join(self.tom_folder,self.tomogram_list.item(int).text()))
		hdr=EMData(self.data_file, 0,True)
		self.nx=hdr["nx"]
		self.ny=hdr["ny"]
		self.nz=hdr["nz"]

		row_count = self.get_inspector().seg_tab.table_set.rowCount()
		for i in range(row_count):
			self.get_inspector().seg_tab.table_set.removeRow(row_count - i - 1)

		seg_path = os.path.join(self.seg_folder,self.tomogram_list.item(int).text()[0:-4]+"_seg.hdf")
		if not os.path.isfile(seg_path):
			seg_out = EMData(self.nx,self.ny,self.nz)
			#self.write_header(seg_out)
			seg_out.write_image(seg_path)
			del seg_out
		else:
			print("Seg file for tomogram {} already exists. Continue ".format(self.tomogram_list.item(int).text()))
			pass
		self.zt_spinbox.setMaximum(self.nz//2)
		#self.data = EMData(self.data_file)


		self.thumbnail.get_im(self.data_file)
		self.thumbnail.set_im()
		self.clear_shapes()

		#TODO - Fix zthick and hard code part

		self.data_xy = self.thumbnail.get_xy()
		print("data xy", self.data_xy)

		if self.get_annotation():
			print("Print annotation to file", self.seg_path)
			#self.write_header(self.get_annotation())
			self.write_out(self.get_annotation(), self.seg_path, self.cur_region)

			#self.get_annotation().write_image(self.seg_path, 0, IMAGE_HDF, False, self.cur_region)
			#self.img_view.inspector.seg_tab.save_all(outfile=self.seg_path, region=self.cur_region)
		else:
			print("Annotation is none.")
			pass
		self.set_imgview_data(round(self.data_xy[0]),round(self.data_xy[1]),self.img_view_region_size)
		self.seg_path = seg_path
		self.reset_morp_params(reset_vs=True)


		# info=js_open_dict("info/annotate_"+self.tomogram_list.item(int).text()+".json")
		# #self.classes = info["class"]
		# try:
		# 	self.boxes = info["boxes"]
		# except:
		# 	self.boxes = []
		#self.add_boxes()

	#need to write region out before setting new data

	def write_out(self,file_out, out_name, region):
		file_out["ann_name"] = 0.111
		file_out["ann_num"] = 2
		file_out.write_image(out_name, 0, IMAGE_HDF, False, region)


	def set_imgview_data(self,x,y,sz):
		iz = self.nz//2
		print("Nz ,iz",self.nz,iz)
		print("Img x,y,sz",x,y,sz)
		if self.nz == 1:
			self.zthick = 0
			self.cur_region = Region(x-sz//2,y-sz//2,sz,sz)
			#self.data = EMData(self.data_file, 0, False, Region(x-old_div(sz,2),y-old_div(sz,2),sz,sz))
		else:
			self.zthick = int(self.zt_spinbox.value())
			print("zthick",self.zthick)
			if self.zthick == -1:
				print(self.nz)
				self.cur_region = Region(x-old_div(sz,2),y-old_div(sz,2),0, sz, sz, self.nz)
				#self.data = EMData(self.data_file, 0, False, Region(x-old_div(sz,2),y-old_div(sz,2),0, sz, sz, self.nz))
			else:
				self.cur_region = Region(x-old_div(sz,2),y-old_div(sz,2),iz-self.zthick, sz, sz,self.zthick*2+1)
				print(self.cur_region)
				#self.data = EMData(self.data_file, 0, False, Region(x-old_div(sz,2),y-old_div(sz,2),iz-self.zthick, sz, sz,self.zthick*2+1))
		self.data = EMData(self.data_file, 0, False, self.cur_region)
		seg_path = os.path.join(self.seg_folder,self.tomogram_list.currentItem().text()[0:-4]+"_seg.hdf")


		self.annotate = EMData(seg_path, 0, False, self.cur_region)
		#self.read_metadata(seg_path)

		self.img_view.set_data(self.data, self.annotate)
		self.img_view.set_scale(1)
		self.img_view.set_origin(0,0)
		#print("Imgview, inspector, segtab",self.img_view,self.img_view.get_inspector(),self.img_view.inspector.seg_tab)
		#self.img_view.get_inspector().seg_tab.read_header(seg_path)
		self.img_view.get_inspector().seg_tab.update_sets()




	def zt_spinbox_changed(self,event):
		self.data_xy = self.thumbnail.get_xy()
		try:
			self.write_out(self.get_annotation(), self.seg_path, self.cur_region)
			#self.get_annotation().write_image(self.seg_path, 0, IMAGE_HDF, False, self.cur_region)
		except:#when annotation files is None
			pass
		self.set_imgview_data(self.data_xy[0],self.data_xy[1],self.img_view_region_size)
		self.reset_morp_params(reset_vs=True)


	def test_drawing_function(self):
		self.create_organelles_training_set()
		return

	def annotate_from_curve(self):
		self.x = self.img_view.get_data().get_xsize()
		self.y = self.img_view.get_data().get_ysize()
		mask = np.zeros((self.x,self.y))

		def trapez(y,y0,w):
			return np.clip(np.minimum(y+1+w/2-y0, -y+1+w/2+y0),0,1)
		def weighted_line(r0, c0, r1, c1, w, rmin=0, rmax=np.inf):
			# The algorithm below works fine if c1 >= c0 and c1-c0 >= abs(r1-r0).
			# If either of these cases are violated, do some switches.
			if c0 == c1 and r0 == r1:
				return 0,0,0
			if abs(c1-c0) < abs(r1-r0):
				# Switch x and y, and switch again when returning.
				xx, yy, val = weighted_line(c0, r0, c1, r1, w, rmin=rmin, rmax=rmax)
				return (yy, xx, val)
			# At this point we know that the distance in columns (x) is greater
			# than that in rows (y). Possibly one more switch if c0 > c1.
			if c0 > c1:
				return weighted_line(r1, c1, r0, c0, w, rmin=rmin, rmax=rmax)
			# The following is now always < 1 in abs
			slope = (r1-r0) / (c1-c0)
			# Adjust weight by the slope
			w *= np.sqrt(1+np.abs(slope)) / 2
			# We write y as a function of x, because the slope is always <= 1
			# (in absolute value)
			x = np.arange(c0, c1+1, dtype=float)
			y = x * slope + (c1*r0-c0*r1) / (c1-c0)
			# Now instead of 2 values for y, we have 2*np.ceil(w/2).
			# All values are 1 except the upmost and bottommost.
			thickness = np.ceil(w/2)
			yy = (np.floor(y).reshape(-1,1) + np.arange(-thickness-1,thickness+2).reshape(1,-1))
			xx = np.repeat(x, yy.shape[1])
			vals = trapez(yy, y.reshape(-1,1), w).flatten()
			yy = yy.flatten()
			# Exclude useless parts and those outside of the interval
			# to avoid parts outside of the picture
			mask = np.logical_and.reduce((yy >= rmin, yy < rmax, vals > 0))
			return (yy[mask].astype(int), xx[mask].astype(int), vals[mask])
		def is_point_in_polygon(point,vertices):
			points = np.array([[x[0],x[1]] for x in vertices])
			poly_path = mplPath.Path(points)
			return poly_path.contains_point(point)
		def polygon(points):
			min_r = int(min([x[0] for x in points]))
			max_r = int(max([x[0] for x in points]))
			min_c = int(min([x[1] for x in points]))
			max_c = int(max([x[1] for x in points]))
			rr = []
			cc = []

			for r in range(min_r,max_r+1):
				for c in range(min_c,max_c+1):
					if is_point_in_polygon([r,c],points):
						rr.append(r)
						cc.append(c)
			return np.array(rr, dtype=np.intp), np.array(cc, dtype=np.intp)

		if self.fill_type == "curve":
			if len(self.curve.points) < 2:
				print("Need at least 2 points to draw a line")
				return
			r = np.array([int(item[0]) for item in self.curve.points])
			c = np.array([int(item[1]) for item in self.curve.points])
			for i in range(len(r)-1):
				try:
					rr, cc, val = weighted_line(r[i], c[i], r[i+1],c[i+1], int(self.tx_line_width.text()))
					mask[rr, cc] = int(self.tx_ann_class.text())
				except:
					continue
			mask[r,c]=int(self.tx_ann_class.text())
		elif self.fill_type == "contour_fill":
			if len(self.contour.points) < 3:
				print("Need at least 3 points to draw a polygon")
				return
			rr, cc = polygon(self.contour.points)
			mask[rr,cc]=int(self.ct_ann_class.text())
		elif self.fill_type == "contour_nofill":
			if len(self.contour.points) < 3:
				print("Need at least 3 points to draw a polygon")
				return
			r_l = [int(item[0]) for item in self.contour.points]
			r_l.append(int(self.contour.points[0][0]))
			r = np.array(r_l)
			c_l = [int(item[1]) for item in self.contour.points]
			c_l.append(int(self.contour.points[0][1]))
			c = np.array(c_l)
			for i in range(len(r)-1):
				try:
					rr, cc, val = weighted_line(r[i], c[i], r[i+1],c[i+1], int(self.ct_line_width.text()))
					mask[rr, cc] = int(self.ct_ann_class.text())
				except:
					continue
			mask[r,c]=int(self.ct_ann_class.text())
		mask =  np.rot90(np.fliplr(mask))
		ann = self.img_view.get_annotation().numpy()
		ann += mask
		self.img_view.force_display_update(set_clip=1)
		self.img_view.updateGL()


	#####ALL BASIC TAB LINEAR AND CONTOUR METHODS
	def basic_tab_change(self, tab_num):

		self.basic_tab_num = tab_num
		# print(self.basic_tab_num)
		if tab_num==3:
			try:
				self.annotate_from_curve()
				self.curve.points = []
				self.contour.points =[]
				self.do_update()
			except:
				pass
			self.img_view.mouse_mode = 1
			return

		elif tab_num==0:
			try:
				self.annotate_from_curve()
				self.curve.points = []
				self.contour.points =[]
				self.do_update()
			except:
				pass
			self.extract_bt_clicked()
			self.img_view.mouse_mode = 6
			#print("Mouse mode is:",self.img_view.mouse_mode_dict[self.img_view.mouse_mode])
			self.img_view.show_inspector(6)
			return

		elif tab_num==1:
			try:
				self.annotate_from_curve()
				self.contour.points =[]
				self.do_update()
			except:
				pass
			self.fill_type = "curve"
			self.curve.points =[]
			self.img_view.mouse_mode = 1
			print("Mouse mode is:",self.img_view.mouse_mode_dict[self.img_view.mouse_mode])

			#Setting EMAnnotate2DWidget mouse mode to emit mode
			return
		else:
			try:
				self.annotate_from_curve()
				self.curve.points =[]
				self.do_update()
			except:
				pass
			self.contour.points =[]
			self.img_view.mouse_mode = 1
			if self.fill_contour_checkbox.isChecked():
				self.fill_type = "contour_fill"
			else:
				self.fill_type = "contour_nofill"
			print("Mouse mode is:",self.img_view.mouse_mode_dict[self.img_view.mouse_mode])
			#Setting EMAnnotate2DWidget mouse mode to emit mode
			return

	def fill_contour_checkbox_changed(self,int):

		if self.fill_contour_checkbox.isChecked():
			self.fill_type = "contour_fill"
		else:
			self.fill_type = "contour_nofill"

	# def update_label(self):
	# 	pts=np.array([p for p in self.curve.points if p[4]==self.curve.classid])
	# 	if len(pts)==0:
	# 		return
	# 	lb=np.unique(pts[:,3])
	# 	txt="{:d} curves\n".format(len(lb))
	# 	txt+="{:d} points\n".format(len(pts))
	# 	txt+='   '+','.join([str(np.sum(pts[:,3]==i)) for i in lb])
	# 	self.lb_lines.setText(txt)

	def do_update(self):
		self.img_view.shapechange=1
		#self.update_label()
		self.img_view.force_display_update(set_clip=0)
		self.img_view.updateGL()


	def interp_bt_clicked(self):

		nppts=np.array(self.curve.points)
		print(nppts)
		if len(self.curve.points) == 0:
			print("No anchor points to interpolate. Return")
			return
		sel=nppts[:,4]==self.curve.classid
		pts=nppts[sel].copy()
		otherpts=nppts[sel==False].copy()

		if len(pts)==0:
			return

		try:
			density=float(self.tx_interp.text())
			print("Interpolate to one points per {:.1f} pixel...".format(density))
		except:
			return

		pts_intp=[]
		kk=0
		for li in np.unique(pts[:,3]):
			pt=pts[pts[:,3]==li][:,:3].copy()
			pt=pt[np.append(True, np.linalg.norm(np.diff(pt, axis=0), axis=1)>0.1)]
			ln=np.linalg.norm(np.diff(pt, axis=0), axis=1)
			ln=np.sum(ln)
			#ln=np.linalg.norm(pt[-1]-pt[0])
			#     print np.round(ln)//2
			if len(pt)<2: continue
			#print len(pt), ln
			ipt=interp_points(pt, npt=np.round(ln)//density)

			pts_intp.append(np.hstack([ipt, kk+np.zeros((len(ipt), 1)), self.curve.classid+np.zeros((len(ipt), 1))]))
			kk+=1

		pts_intp=np.vstack(pts_intp)
		self.curve.points=otherpts.tolist()
		self.curve.points.extend(pts_intp.tolist())
		#print("Int points:",self.curve.points)
		self.do_update()
		#self.save_points()

	def clear_points(self,ask=True):
		if ask:
			choice = QtWidgets.QMessageBox.question(self, 'Clear curve points',
				'Clear all curve anchor points in the tomogram?', QtWidgets.QMessageBox.Yes, QtWidgets.QMessageBox.No)
			if choice == QtWidgets.QMessageBox.Yes:
				self.curve.points=[]
				self.contour.points=[]
		else:
			self.curve.points=[]
			self.contour.points=[]
		self.do_update()


		#self.save_points()
		return

	def get_zpos(self):

		try:
			zpos = self.img_view.inspector.ns.value
			return zpos
		except:
			#TODO
			#self.img_view.show_inspector(6)
			print("Except in get zpos")
			return 0

	def get_boxsize(self):
		return int(self.bsz_vs.value)
	def inside_box(self,i,x=-1,y=-1,z=0):
		#Checks to see if a point in image coordinates is inside box number n. If any value is negative, it will not be checked.
		box = self.boxes[i]
		bs=self.get_boxsize()/2

		rr=(x>=0)*((box[0]-x)**2) + (y>=0)*((box[1]-y) **2)
		#print("rr",rr,"xy",x,y,"box_xy",self.boxes[i])
		#+ (z>=0)*((box[2]-z)**2)
		# else:
		# 	rr=(x>=0)*((box[0]-x)**2) + (y>=0)*((box[1]-y) **2) + (z>=0)*(box[2]!=z)*(1e3*bs**2)
		if rr<=bs**2: print(i)
		return rr<=bs**2



	def key_press(self, event):
		if event.key() == Qt.Key_Up or event.key() == Qt.Key_Down:
			if self.img_view.data:
				print("Current zpos", self.get_zpos())
				self.update_img_view_box()
		else:
			return

	def update_img_view_box(self):
		current_z = self.get_zpos()
		self.clear_shapes(reset_boxes_list=False)
		for box in self.boxes:
			if box[3] != -1:
				if box[2] == current_z:
					box[3]=1
				else:
					box[3]=0
			else:
				continue
		self.add_boxes(size=int(self.bsz_vs.value))


	def img_view_mouse_up(self, event):
		#Boxer tab
		if self.basic_tab_num == 3:
			x,y=self.img_view.scr_to_img((event.x(),event.y()))
			z =self.img_view.zpos
			print(x,y,self.img_view.zpos)
			if not event.button()&Qt.LeftButton:
				return
			if event.modifiers()&Qt.ShiftModifier:
				for i in range(len(self.boxes)):
					if self.inside_box(i,x,y,0):
						self.boxes[i][3]=-1
						self.img_view.del_shape("box{}".format(i))
						self.img_view.updateGL()
			# 			#### remove point
			# 			self.curve.points=[p for i,p in enumerate(self.curve.points) if i!=ii]
			else:
				self.boxes.append([x, y, z,1])
				self.add_boxes(size=int(self.bsz_vs.value))
				print(self.boxes)
		#Contour tab
		elif self.basic_tab_num == 1:
			#print("Mouse is to draw contour")
			x,y=self.img_view.scr_to_img((event.x(),event.y()))
			if not event.button()&Qt.LeftButton:
				return
			if event.modifiers()&Qt.ControlModifier or event.modifiers()&Qt.ShiftModifier:
				#### check if clicking on an existing point
				pts=np.array(self.curve.points)

				if len(pts)<1:
					return
				ofpln=abs(pts[:,2]-self.img_view.list_idx)>3

				res=np.sqrt(np.sum((pts[:,:2]-np.array([x,y]))**2, axis=1))
				res[ofpln]+=1e5
				ii=np.argmin(res)

				if event.modifiers()&Qt.ShiftModifier:
					if res[ii]<20:
						#### remove point
						self.curve.points=[p for i,p in enumerate(self.curve.points) if i!=ii]

				elif event.modifiers()&Qt.ControlModifier:
					if res[ii]<20:
						#### select contour
						ci=self.curve.points[ii][3]
						print('select curve {:d}'.format(int(ci)))
						self.curve.select=ci

					else:
						print('new curve')
						self.annotate_from_curve()
						self.curve.points = []
						self.curve.add_point([x, y], True)

			else:
				#### add point
				self.curve.add_point([x, y]) #, self.img_view.list_idx
			self.do_update()
			print(self.curve.points)

		#Curve tab
		elif self.basic_tab_num == 2:
			x,y=self.img_view.scr_to_img((event.x(),event.y()))
			#if event.button()&Qt.LeftButton:
			if event.modifiers()&Qt.ControlModifier:
				#### interpolate curve from previous slices
				print('new contour')
				self.annotate_from_curve()
				self.contour.points = []
				self.contour.add_point([x, y], True)
			elif event.modifiers()&Qt.ShiftModifier:
				#### remove point
				pts=np.array(self.contour.points)
				pts=pts[pts[:,2]==self.img_view.list_idx,:]
				if len(pts)<1:
					return
				res=np.sqrt(np.sum((pts[:,:2]-np.array([x,y]))**2, axis=1))
				ii=np.argmin(res)
				if res[ii]<20:
					pts=[[p[0], p[1], self.img_view.list_idx, p[3]] for i,p in enumerate(pts) if i!=ii]
					self.contour.points=[p for p in self.contour.points if p[2]!=self.img_view.list_idx]
					self.contour.points.extend(pts)
			else:
				#### add point
				self.contour.add_point([x, y]) #, self.img_view.list_idx
			self.img_view.shapechange=1
			print(self.contour.points)
			self.img_view.updateGL()

		#Brush tab
		else:
			print("I don't know yet")


	#TODO
	def random_bx_bt_clicked(self):
		try:
			no_box = int(self.random_bx_sb.getValue())
		except:
			print("Number of boxes needed to be integer")
			return
		xs = np.random.randint(self.bsz_vs.value//2,self.img_view.get_data_dims()[0]-self.bsz_vs.value//2,no_box)
		ys = np.random.randint(self.bsz_vs.value//2,self.img_view.get_data_dims()[1]-self.bsz_vs.value//2,no_box)
		for i in range(no_box):
			self.boxes.append([xs[i],ys[i],self.get_zpos(),1])
		self.add_boxes(size=int(self.bsz_vs.value))
		#self.img_view.updateGL()


	def clear_bx_bt_clicked(self):
		print("Delete all boxes")
		self.clear_shapes()

	def extract_bt_clicked(self):
		print("Extract image patches of size",self.bsz_vs.value,"at positions:",self.boxes)
		outfile = "./particles/bg_temp.hdf"
		try:
			os.remove(outfile)
		except:
			pass
		for i in range(len(self.boxes)):
			x,y,z = int(self.boxes[i][0]),int(self.boxes[i][1]),int(self.boxes[i][2])
			bs = self.bsz_vs.value
			if self.boxes[i][3] == 1:
				r = self.data.get_clip(Region(x-bs//2,y-bs//2,z,bs,bs,1))
				r.write_image(outfile,-1)
				#r=self.data.get_clip(Region(x-bs//2,y-bs//2,z-bz//2,bs,bs,bz))
		self.clear_shapes()



	def bsz_vs_value_changed(self,event):
		print("New region size:",self.bsz_vs.value)
		#self.clear_shapes()
		self.img_view.updateGL()
		self.add_boxes(size=int(self.bsz_vs.value))

	def on_list_tomo_clicked(self):
		return

	def clear_shapes(self,reset_boxes_list=True):
		for i in range(len(self.boxes)):
			#x,y = int(self.boxes[i][0]),int(self.boxes[i][1])
			self.img_view.del_shape("box{}".format(i))
		#self.img_view.updateGL()
		if reset_boxes_list:
			self.boxes = []
		self.clear_points(ask=False)
		self.do_update()
		self.img_view.updateGL()

	def on_table_tom(self,row,col):
		print("Row", row, "Col", col, "Tomogram", self.table_tom.itemAt(row,col))
		print(self.table_tom.currentItem().text())

	def table_tom_cell_changed(self):
		#print("New Cell", self.table_tom.currentItem().text())
		self.img_view.set_data(EMData(self.table_tom.currentItem().text()),None)
	def set_data(self):
		return

	def update_sets(self):
		#Set the colors and flags of table set items
		for i in range(self.table_tom.rowCount()):
			key = int(self.table_tom.item(i,0).text())
			self.table_tom.item(i,0).setFlags(self.indexflags)
			self.table_tom.item(i,1).setFlags(self.itemflags)
			self.table_tom.item(i,0).setForeground(self.colors[key])
			self.table_tom.item(i,1).setForeground(self.colors[key])



	def add_boxes(self, size = 64):
		sz = size
		color = [0.1,0.1,0.3]
		for i in range(len(self.boxes)):

			x,y = int(self.boxes[i][0]),int(self.boxes[i][1])

			if self.boxes[i][3] == 1:
				self.img_view.add_shape("box{}".format(i),EMShape(("rect",color[0],color[1],color[2],x-old_div(sz,2),y-old_div(sz,2),x+old_div((sz+1),2),y+old_div((sz+1),2),2)))
		self.img_view.updateGL()

	def SaveJson(self):
		info=js_open_dict("./info/annotate_project.json")
		info["tom_list"]=self.tom_file_list
		info.close()


	def write_metadata(self, seg_path):
		file_out = EMData(seg_path,0,True)
		row_count = self.get_inspector().seg_tab.table_set.rowCount()
		print("Row Count",int(row_count))
		if row_count == 0:
			print("Table set is empty, no classes was saved in metadata")
		else:
			print("Writing classes to metadata")
			nums = [int(self.get_inspector().seg_tab.table_set.item(row,0).text()) for row in range(row_count)]
			names = [str(self.get_inspector().seg_tab.table_set.item(row,1).text()) for row in range(row_count)]
			serialize_name = json.dumps(names, default=lambda a: "[%s,%s]" % (str(type(a)), a.pk))
			file_out["ann_name"] = serialize_name
			file_out["ann_num"] = nums
		return

	def read_metadata(self,seg_path):
		file_in = EMData(seg_path,0,True)
		#try:
		keys = file_in["ann_num"]
		values = json.loads(file_in["ann_name"])
		#print(values)
		if type(keys) != list: keys = [keys]
		item_dict=(dict(zip(keys, values)))
		print(item_dict)
		if ret:
			return item_dict
		else:
			for key, value in item_dict.items():
				self.get_inspector().seg_tab.add_new_row(key,value)
			self.get_inspector().seg_tab.update_sets()
		# except:
		# 	print("Trouble reading headers information. Continue...")
		# 	pass

		return

	def closeEvent(self,event):
		"""Close everything when close the ImageViewer"""
		print("Exiting")
		#self.SaveJson()
		E2saveappwin("e2annotate","main",self)
		E2saveappwin("e2annotate","controlpanel",self.control_panel)
		E2saveappwin("e2annotate","tomograms",self.tomo_list_panel)
		#self.write_header(self.get_annotation())

		self.write_metadata(self.seg_path)
		self.write_out(self.get_annotation(), self.seg_path, self.cur_region)

		#self.get_annotation().write_image(self.seg_path, 0, IMAGE_HDF, False, self.cur_region)
		self.close()
		self.control_panel.close()
		self.tomo_list_panel.close()

class Thumbnail(EMImage2DWidget):
	def __init__(self,current_file=None,target=None,app_target = None,tn_size=220):
		super().__init__()
		#self.full_im = EMData(current_file,0)
		#self.thumbnail_image=self.target.get_data().process("math.meanshrink",{"n":(ceil(self.scale_fac))})

		if current_file:
			self.current_file = current_file
			self.get_im(self.current_file)

		self.scale_fac = round(self.img["nx"]/tn_size)
		self.thumbnail_image=self.img.process("math.meanshrink",{"n":(self.scale_fac)})
		self.size = self.thumbnail_image.get_sizes()[0]

		if target:
			self.target = target
		if app_target:
			self.app_target = app_target

		self.im_xsize= self.img["nx"]

		#WinView Box
		self.box_size=self.get_box_size()
		#print("box_size", self.box_size)
		self.x = self.box_size[0]/2
		self.y = self.box_size[1]/2
		self.set_im()
		self.app_target.set_imgview_data(self.get_xy()[0],self.get_xy()[1],self.app_target.img_view_region_size)


		#self.glass_scale
	def get_xy(self):
		return [self.x*self.get_scale_fac(), self.y*self.get_scale_fac()]

	def get_im(self, data_file):
		print("Current_file", data_file)
		hdr=EMData(data_file, 0,True)
		iz=hdr["nz"]//2
		print(hdr["nz"],iz)
		if hdr["nz"] == 1:
			self.img = EMData(data_file, 0, False)
		else:
			self.img = EMData(data_file, 0, False, Region(0,0,iz, hdr["nx"], hdr["ny"],1))

	def set_im(self):
		#print("Imxsize,scale,size",self.im_xsize,self.target.scale,self.size)
		#self.thumbnail_image=self.target.get_data().process("math.meanshrink",{"n":(ceil(self.scale_fac))})
		self.thumbnail_image=self.img.process("math.meanshrink",{"n":(round(self.scale_fac))})
		print("tn_image_size",self.thumbnail_image.get_sizes())
		self.size = self.thumbnail_image.get_sizes()[0]

		# tmp = self.thumbnail_image.copy().process_inplace("mult",{"value":-1})
		# print("Temp", tmp)
		self.set_data(self.thumbnail_image)
		self.add_box(self.box_size[0]/2+0,self.box_size[1]/2+0,self.box_size)
		#print("Boxsize", self.box_size, "self.iv_fixed_size", self.iv_size)

	def get_box_size(self):
		#self.im_scalesz = self.im_xsize * self.target.scale
		self.im_scalesz = self.im_xsize
		#print("Hahah",self.im_xsize,self.target.scale)
		#self.iv_size = [self.target.size().width(), self.target.size().height()]
		self.iv_size = [self.app_target.img_view_region_size, self.app_target.img_view_region_size]
		box_size = [self.size*(self.iv_size[0]/self.im_scalesz),self.size*(self.iv_size[1]/self.im_scalesz)]
		#print("IV size", self.iv_size, "Box_size", box_size)
		return box_size

	def get_scale_fac(self):

		return self.im_xsize*self.target.scale/self.size

	def mousePressEvent(self, event):
		return

	def mouseMoveEvent(self, event):
		lc=self.scr_to_img(event.x(),event.y())
		if event.buttons()&Qt.LeftButton:
			self.box_size = self.get_box_size()
			self.scale_fac = self.get_scale_fac()
			self.add_box((lc[0]),(lc[1]), (self.box_size))
		elif event.buttons()&Qt.RightButton:
			return

		else:
			return

	def mouseReleaseEvent(self, event):
		get_application().setOverrideCursor(Qt.ArrowCursor)
		print("Mouse release")
		lc=self.scr_to_img(event.x(),event.y())
		if event.button()==Qt.LeftButton:
			lc=self.scr_to_img(event.x(),event.y())
			xy = [lc[0],lc[1]]
			print("Box x,y,bound", lc[0], lc[1], self.box_size[0]/2, self.size-self.box_size[0]/2)
			if lc[0] <= self.box_size[0]/2:
				print("x is out of bound, move position to x =",self.box_size[0]/2)
				xy[0]=self.box_size[0]/2
			elif lc[0] >= self.size-self.box_size[0]/2:
				print("x is out of bound, move position to x =",self.size-self.box_size[0]/2)
				xy[0]=self.size-self.box_size[0]/2
			if lc[1] <= self.box_size[1]/2:
				print("y is out of bound, move position to y =",self.box_size[1]/2)
				xy[1]=self.box_size[1]/2
			elif lc[1] >= (self.size-self.box_size[1]/2):
				print("y is out of bound, move position to y =",self.size-self.box_size[1]/2)
				xy[1]=self.size-self.box_size[1]/2
			else:
				pass

			print("Mouse release at", lc, "mouse position set to", xy)
			self.box_size = self.get_box_size()
			self.scale_fac = self.get_scale_fac()
			print("Bzsz", self.box_size, "scale fac", self.scale_fac)
			self.add_box(xy[0],xy[1], self.box_size)
			if self.app_target.get_annotation():
				print("Print annotation to file", self.app_target.seg_path)
				try:
					#self.write_header(self.app_target.get_annotation())
					self.app_target.write_out(self.app_target.get_annotation(), self.app_target.seg_path, self.app_target.cur_region)

					#self.app_target.get_annotation().write_image(self.app_target.seg_path, 0, IMAGE_HDF, False, self.app_target.cur_region)
				#self.img_view.inspector.seg_tab.save_all(outfile=self.seg_path, region=self.cur_region)
				except:
					print("Cannot write to region out of bound. Continue")
					pass
			else:
				print("Annotation is none.")
				pass
			self.app_target.set_imgview_data((xy[0])*self.scale_fac,(xy[1])*self.scale_fac,self.app_target.img_view_region_size)
		else:
			print("or else")
			return

	#Handle winview on thumbnail
	def target_rmousedrag(self):
		self.box_size = self.get_box_size()
		self.scale_fac = self.get_scale_fac()
		self.x = (self.target.get_origin()[0])/self.scale_fac+self.box_size[0]/2
		self.y = (self.target.get_origin()[1])/self.scale_fac+self.box_size[1]/2
		#print("x,y", self.x, self.y)
	def update_box(self):
		self.box_size = self.get_box_size()
		#self.scale_fac = self.get_scale_fac()
		self.add_box(self.x,self.y, self.box_size)
		#print("x,y", self.x, self.y)

	def add_box(self, x, y, box_sz):
		#bound=[int(self.box_size//2),int(self.size-self.box_size//2)]
		x,y = x,y
		sz = (box_sz)
		print("box x,y,sz", x,y,sz)
		self.del_shape("WinView")
		self.add_shape("WinView",EMShape(("rectpoint",.5,.5,.1,x-old_div(sz[0],2),y-old_div(sz[1],2),x+old_div((sz[0]+1),2),y+old_div((sz[1]+1),2),2)))
		#self.add_shape("ClipView",EMShape(("rectpoint",1,.5,.1,x-old_div(sz[0],2),y-old_div(sz[1],2),x+old_div((sz[0]+31),2),y+old_div((sz[1]+31),2),2)))
		self.x = x
		self.y = y
		self.updateGL()

	def wheelEvent(self, event):
		return


class Curve(EMShape):
	def __init__(self, img, points):
		self.points=points
		print(points)
		self.isanimated = False
		self.shape=["scr_contour",0]
		self.image=img
		self.triangles=[]
		self.select=0
		self.classid=0

	def add_point(self, newpt=[], newcontour=False, optimize=True):
		#zpos=self.image.list_idx
		zpos = 0
		if len(self.points)==0 and len(newpt)>0:
			#self.points.append([newpt[0], newpt[1], zpos, self.select, self.classid])
			self.points.append([newpt[0], newpt[1], zpos, self.select, self.classid])

		if newcontour==False:

			ci=self.select
			#### separate the points on the current contour and everything else
			nppts=np.array(self.points)
			sel=np.logical_and(nppts[:,3]==ci, nppts[:,4]==self.classid)
			pts=nppts[sel,:3].copy()
			otherpts=nppts[sel==False].copy()

			if len(pts)<3 or optimize==False:
				if len(newpt)>0:
					self.points.append([newpt[0], newpt[1], zpos, ci, self.classid])
					return
			else:
				thr=1000.
				newpt=np.array([newpt[0], newpt[1], zpos])


			if len(newpt)>0:
				#### add a point
				dst=np.array(scipydist.cdist([newpt], pts)[0])
				rg=np.arange(len(dst), dtype=int)+1
				#print rg
				rg[-1]=rg[-2]
				d0=dst[0]
				d1=dst[-1]
				dst+=dst[rg]
				dst=np.append(d0*3, dst)
				dst[-1]=d1*3

				mi=np.argmin(dst)
				#print(mi, len(pts))
				#ci=pts[mi, 3]
				if mi==0:
					pts=np.vstack([ newpt, pts])
				if mi==len(pts):
					pts=np.vstack([pts, newpt])
				pts=np.insert(pts, mi+0, newpt, axis=0)


			allpts=[]
			#### a simple tsp solver...
			#idx=np.where(pts[:,3]==ci)[0]
			pp=pts.copy()
			path=np.arange(len(pp), dtype=int)

			dmat=scipydist.squareform(scipydist.pdist(pp))
			dmat+=np.eye(len(dmat))*thr
			if len(pp)>=3:
				niter=2000
			else:
				niter=0

			calc_dist=lambda path: np.sum([dmat[path[i], path[i+1]] for i in range(len(path)-1)])
			dst=calc_dist(path)
			nochange=0
			for k in range(niter):
				p0=path.copy()
				i0,i1= np.sort(np.random.randint(0,len(pp)+1, 2))
				if abs(i0-i1)%len(path)<2: continue
				path= np.hstack([path[:i0], path[i0:i1][::-1], path[i1:]])
				d=calc_dist(path)
				#print(i0, i1, d, dst)
				#print(path)
				if d>=dst:
					path=p0
					nochange+=1
					if nochange>200: break

				else:
					dst=d
					nochange=0


			allpts=[[pp[i][0], pp[i][1], pp[i][2], ci, self.classid] for i in path]

			self.points=otherpts.tolist()
			self.points.extend(allpts)
			#self.points=allpts#.tolist()

		else:
			#### start a new contour
			pts=np.array(self.points)
			ci=np.max(pts[:,3])+1
			self.select=ci
			self.points.append([newpt[0], newpt[1], zpos, ci, self.classid])


	def draw(self,d2s=None,col=None):
		zpos=self.image.list_idx
		curpts=[p for p in self.points if p[4]==self.classid]

		#print np.array(self.points)
		#print "#########"

		cid=np.unique([p[3] for p in curpts])
		for ci in cid:
			#### draw lines
			pts=np.array([[p[0], p[1], p[2]] for p in curpts if ci==p[3]])
			dzs=pts[:,2]-zpos
			lns=[]
			for i in range(len(pts)-1):
				if dzs[i]*dzs[i+1]<=0: ### line cross current plane
					p0=pts[i]
					p1=pts[i+1]
					if (p1[2]==p0[2]):
						q0=p0
						q1=p1
					else:
						dz=abs(old_div(float(zpos-p0[2]),(p1[2]-p0[2])))
						dp=old_div(1.,abs((p1[2]-p0[2])))
						#print dp, dz
						d0=max(0, (dz-dp))
						d1=min(1, (dz+dp))
						q0=p1*d0+p0*(1-d0)
						q1=p1*d1+p0*(1-d1)

					lns.append([q0[0], q0[1]])
					lns.append([q1[0], q1[1]])

			if ci==self.select:
				clr0=.3
			else:
				clr0=.7

			glColor3f( 1., clr0, clr0 );
			glLineWidth(3.)
			glEnableClientState(GL_VERTEX_ARRAY)
			glVertexPointerf(lns)
			glDrawArrays(GL_LINES, 0, len(lns))


		for p in curpts:
			#### draw nodes on the plane
			#pts=[[p[0], p[1]] for p in self.points if p[2]==zpos]
			s=14.-abs(p[2]-zpos)
			if s<=0: continue
			if p[3]==self.select:
				clr0=.3
			else:
				clr0=.7
			glColor3f( clr0, clr0, 1 );
			glPointSize(s)
			glBegin(GL_POINTS)
			glVertex(p[0], p[1], 0)

			glEnd()
			#glEnableClientState(GL_VERTEX_ARRAY)
			#glVertexPointerf(pts)
			#glDrawArrays(GL_POINTS, 0, len(pts))

		return

class Contour(EMShape):
	def __init__(self, img=None, points=[]):
		self.points=points
		self.isanimated = False
		self.shape=["scr_contour",0]
		self.image=img
		self.triangles=[]
		#self.lines=[]
#		if len(points)>3:
#		self.make_triangle()

	def add_point(self, newpt=[], newcontour=False):
		zpos=self.image.list_idx
		pts=np.array(self.points)
		if len(pts)>=3:
			pts=pts[pts[:,2]==zpos,:]

		if len(pts)<3:
			if len(newpt)>0:
				self.points.append([newpt[0], newpt[1], zpos, 0])
		else:
			thr=1000.
			if newcontour==False:
				cid=np.unique(pts[:,3])
				if len(newpt)>0:
					#### add a point
					dst=np.array(scipydist.cdist([newpt], pts[:,:2])[0])
					rg=np.arange(len(dst), dtype=int)+1
					for ci in cid:
						idx=np.where(pts[:,3]==ci)[0]
						rg[idx[-1]]=idx[0]
					#print rg
					dst+=dst[rg]

					mi=np.argmin(dst)
					ci=pts[mi, 3]
					pts=np.insert(pts, mi+1, np.append(newpt, [zpos, ci]), axis=0)
				allpts=[]
				for ci in cid:
					#### a simple tsp solver...
					idx=np.where(pts[:,3]==ci)[0]
					pp=pts[idx].copy()
					path=np.arange(len(pp), dtype=int)

					dmat=scipydist.squareform(scipydist.pdist(pp[:,:2]))
					dmat+=np.eye(len(dmat))*thr
					if len(pp)>3:
						niter=2000
					else:
						niter=0

					calc_dist=lambda path: np.sum([dmat[path[i], path[(i+1)%len(path)]] for i in range(len(path))])
					dst=calc_dist(path)
					nochange=0
					for k in range(niter):
						p0=path.copy()
						i0,i1= np.sort(np.random.randint(0,len(pp), 2))
						if abs(i0-i1)%len(path)<2: continue
						path= np.hstack([path[:i0+1], path[i0+1:i1+1][::-1], path[i1+1:]])
						d=calc_dist(path)
						if d>=dst:
							path=p0
							nochange+=1
							if nochange>200: break

						else:
							dst=d
							nochange=0


					allpts.extend([[pp[i][0], pp[i][1], zpos, ci] for i in path])

				self.points=[p for p in self.points if p[2]!=zpos]
				self.points.extend(allpts)

			else:
					#### start a new contour
					ci=np.max(pts[:,3])+1
					self.points.append([newpt[0], newpt[1], zpos, ci])

		np.savetxt("pts.save", self.points, fmt='%d')

	def next_slice(self):
		pts=np.array(self.points)
		mi=self.image.list_idx
		ii=np.argmin(abs(pts[:,2]-mi))
		if pts[ii,2]==mi: return
		last=pts[ii,2]
		pts=pts[pts[:,2]==last]
		#print(mi, last, pts.shape)
		img=self.image.data.numpy()

		vec=[]
		cid=np.unique(pts[:,3])
		rg0=np.arange(len(pts), dtype=int)-1
		rg1=np.arange(len(pts), dtype=int)+1
		for ci in cid:
			idx=np.where(pts[:,3]==ci)[0]
			rg1[idx[-1]]=idx[0]
			rg0[idx[0]]=idx[-1]

		#for i in range(len(pts)):
		vec= pts[rg0]-pts[rg1]
		vec=np.vstack([vec[:,1], -vec[:,0]]).T
		vec/=np.linalg.norm(vec, axis=1)[:, None]
		#vec.append(v)
		#vec=np.array(vec)
		pval=[]
		rg=np.arange(-2,2.1)
		for i in rg:
			p=(pts[:,:2]-vec*i).astype(int)
			out=(p[:,0]<0)+(p[:,1]<0)+(p[:,0]>=img.shape[0])+(p[:,1]>=img.shape[1])
			out=out>0
			p[out,:]=0
			pval.append(img[p[:,1], p[:,0]])
			pval[-1][out]=1000

		pval=np.array(pval)

		vp=rg[np.argmin(pval, 0)]
		p1=pts.copy()
		p1[:,:2]-=vp[:,None]*vec
		p1[:,:2]=np.round(p1[:,:2])
#		p1[-1]=p1[0]
		p1[:,2]=mi
		self.points.extend(p1.tolist())
	def draw(self,d2s=None,col=None):

		zpos=self.image.list_idx
		allpts=[[p[0], p[1], p[3]] for p in self.points if p[2]==zpos]
		#print(np.array(self.points))
		#print("#########")
		cid=np.unique([p[2] for p in allpts])
		for ci in cid:
			pts=[[p[0], p[1]] for p in allpts if p[2]==ci]
			if len(pts)>2:
				pts.append(pts[0])

			area=0.
			for i in range(len(pts)):
				p0=pts[i]
				p1=pts[(i+1)%len(pts)]
				area+=p0[0]*p1[1]-p0[1]*p1[0]
			area=abs(area/2.)
			print("Contour {:d}, area {:.1f} px^2".format(int(ci), area))

			glColor3f( 1, .3, .3 );
			glLineWidth(3.)
			glEnableClientState(GL_VERTEX_ARRAY)
			glVertexPointerf(pts)
			glDrawArrays(GL_LINE_STRIP, 0, len(pts))


			glColor3f( .3, .3, 1 );
			glPointSize(7.)
			glEnableClientState(GL_VERTEX_ARRAY)
			glVertexPointerf(pts)
			glDrawArrays(GL_POINTS, 0, len(pts))

class UNet():
	def __init__(self, infile=None, data=None, label=None, batchsz=50 ):

		#self.model=self.get_tiny_unet(tile_sz,tile_sz)
		self.model=self.get_unet(64,64)
		print("Getting a new unet model")

		self.datas=None
		self.labels=None
		# self.data_img=data
		# self.label_img=label
		#self.data_img.write_image("./neural_nets/test_data.hdf")
		if infile:
			self.datas, self.labels = self.load_particles_from_file(ptcls=infile)
		#self.create_edge_training_set()
			print("Done creating training set")
			print("Datas: ",self.datas.shape,"Labels: ",self.labels.shape)


	def create_edge_training_set(self,tile_sz=6,disk_sz=5):
		coor = np.where(self.label_img.numpy()>0)
		data_l=[]
		label_l=[]
		for i in range(15,len(coor[0])-15):
			x = int(coor[1][i])
			y = int(coor[0][i])
			d = self.data_img.get_clip(Region(x,y,tile_sz,tile_sz))
			l = self.label_img.get_clip(Region(x,y,tile_sz,tile_sz))
			# d.write_image("./neural_nets/data_stack_train.hdf",-1)
			# l.write_image("./neural_nets/label_stack_train.hdf",-1)
			data_l.append(d.numpy())
			label_l.append(l.numpy())
		# print(data_l[15])
		self.datas = np.asarray(data_l,dtype=np.float32)
		self.datas.reshape(-1,tile_sz,tile_sz,1)
		self.datas /=3.0
		self.labels=np.asarray(label_l,dtype=np.float32)

	def get_tiny_unet(self,inp_x=None,inp_y=None):
		inputs = Input((inp_x, inp_y, 1))
		conv1 = Conv2D(32, (5, 5), activation='relu', padding='same')(inputs)
		pool1 = MaxPooling2D(pool_size=(2, 2))(conv1)
		conv2 = Conv2D(64, (5, 5), activation='relu', padding='same')(pool1)
		up1 = concatenate([Conv2DTranspose(64, (2, 2), strides=(2, 2), padding='same')(conv2), conv1], axis=3)
		conv3 = Conv2D(64, (5, 5), activation='relu', padding='same')(up1)
		conv4 = Conv2D(1, (1, 1), activation='sigmoid')(conv3)
		model = Model(inputs=[inputs], outputs=[conv4])
		return model

	def load_particles_from_file(self, ptcls=None, ncopy=1, rng=None):
		is3d=False
		e=EMData(ptcls,0, True)
		tsz=max(e["nx"],e["ny"])
		nframe=EMUtil.get_image_count(ptcls)
		if nframe==1:
			nframe=e["nz"]
			if nframe>1:
				is3d=True
		num = nframe//2
		data=[]
		label=[]
		ntrain=-1
		for i in range(num):
			for nc in range(ncopy):
				if is3d:
					ptl=EMData(ptcls,0, False,Region(0,0,i*2,tsz,tsz,1))
				else:
					ptl=EMData(ptcls,i*2, False, Region(0,0,tsz,tsz))

				if ntrain<0 and ptl.get_attr_default("valid_set", 0)==1:
					ntrain=len(data)
				#ptl.process_inplace("threshold.belowtozero")

				ar=ptl.numpy().copy()
				#shp=np.shape(ar)
				data.append(ar)

				if is3d:
					ptl=EMData(ptcls,0, False,Region(0,0,i*2+1,tsz,tsz,1))
				else:
					ptl=EMData(ptcls,i*2+1, False, Region(0,0,tsz,tsz))
					#ptl.process_inplace("threshold.belowtozero")

				ar=ptl.numpy().copy()
				#shp=np.shape(ar)
				label.append(ar)

		if ntrain<0: ntrain=len(data)
		print("{:d} particles loaded, {:d} in training set, {:d} in validation set".format(len(data), ntrain, len(data)-ntrain))
		data=np.asarray(data,dtype=np.float32)

		print(data.shape)
		print("Std of particles: ",np.std(data))
		data/=3.

		label=np.asarray(label,dtype=np.float32)
		label/=(np.max(np.abs(label)))

		header=EMData(ptcls,0,True)
		shape=[header["nx"],header["ny"],header["nz"]]
		return self.normalize(data), label


	def get_unet(self,inp_x=None,inp_y=None):
		inputs = Input((inp_x, inp_y, 1))
		conv1 = Conv2D(32, (3, 3), activation='relu', padding='same')(inputs)
		conv1 = Conv2D(32, (3, 3), activation='relu', padding='same')(conv1)
		pool1 = MaxPooling2D(pool_size=(2, 2))(conv1)

		conv2 = Conv2D(64, (3, 3), activation='relu', padding='same')(pool1)
		conv2 = Conv2D(64, (3, 3), activation='relu', padding='same')(conv2)
		pool2 = MaxPooling2D(pool_size=(2, 2))(conv2)

		conv3 = Conv2D(128, (3, 3), activation='relu', padding='same')(pool2)
		conv3 = Conv2D(128, (3, 3), activation='relu', padding='same')(conv3)
		pool3 = MaxPooling2D(pool_size=(2, 2))(conv3)

		conv4 = Conv2D(256, (3, 3), activation='relu', padding='same')(pool3)
		conv4 = Conv2D(256, (3, 3), activation='relu', padding='same')(conv4)
		#pool4 = MaxPooling2D(pool_size=(2, 2))(conv4)

		#conv5 = Conv2D(512, (3, 3), activation='relu', padding='same')(pool4)
		#conv5 = Conv2D(512, (3, 3), activation='relu', padding='same')(conv5)

		#up6 = concatenate([Conv2DTranspose(256, (2, 2), strides=(2, 2), padding='same')(conv5), conv4], axis=3)
		#conv6 = Conv2D(256, (3, 3), activation='relu', padding='same')(up6)
		#conv6 = Conv2D(256, (3, 3), activation='relu', padding='same')(conv6)

		#up7 = concatenate([Conv2DTranspose(128, (2, 2), strides=(2, 2), padding='same')(conv6), conv3], axis=3)
		up7 = concatenate([Conv2DTranspose(128, (2, 2), strides=(2, 2), padding='same')(conv4), conv3], axis=3)
		conv7 = Conv2D(128, (3, 3), activation='relu', padding='same')(up7)
		conv7 = Conv2D(128, (3, 3), activation='relu', padding='same')(conv7)

		up8 = concatenate([Conv2DTranspose(64, (2, 2), strides=(2, 2), padding='same')(conv7), conv2], axis=3)
		conv8 = Conv2D(64, (3, 3), activation='relu', padding='same')(up8)
		conv8 = Conv2D(64, (3, 3), activation='relu', padding='same')(conv8)

		up9 = concatenate([Conv2DTranspose(32, (2, 2), strides=(2, 2), padding='same')(conv8), conv1], axis=3)
		conv9 = Conv2D(32, (3, 3), activation='relu', padding='same')(up9)
		conv9 = Conv2D(32, (3, 3), activation='relu', padding='same')(conv9)

		conv10 = Conv2D(1, (1, 1), activation='sigmoid')(conv9)

		model = Model(inputs=[inputs], outputs=[conv10])

		return model


	def train_unet(self,weights_out='./neural_nets/weights_temp.h5',no_epoch = 30, batch_sz= 50,val_split=0.2,learnrate=3e-4):
		def dice_coef(y_true, y_pred):
			y_true_f = K.flatten(y_true)
			y_pred_f = K.flatten(y_pred)
			intersection = K.sum(y_true_f * y_pred_f)
			return (2. * intersection ) / (K.sum(y_true_f) + K.sum(y_pred_f) )

		def dice_coef_loss(y_true, y_pred):
			return -dice_coef(y_true, y_pred)

		self.weights_out = weights_out
		print('-'*30)
		print('Loading and preprocessing train data...')
		print('-'*30)

		try:
			os.remove(self.weights_out)
			#os.remove('./neural_nets/weights_temp.h5')
		except:
			pass
		self.model.compile(optimizer=Adam(learning_rate=learnrate), loss=dice_coef_loss, metrics=[dice_coef])
		model_checkpoint = ModelCheckpoint(self.weights_out, monitor='val_loss', save_best_only=True)

		print('-'*30)
		print('Fitting model...')
		print('-'*30)
		self.model.fit(self.datas, self.labels, batch_size=batch_sz, epochs=no_epoch, verbose=1, shuffle=True,validation_split=val_split,callbacks=[model_checkpoint])
		#return model,history
		print("Done training model. Network saved to",weights_out)

	def load_model(self,weights_in, tiny=False):
		if tiny:
			model = self.get_tiny_unet()
		else:
			model = self.get_unet()
		model.load_weights(weights_in)
		return model

	# def apply_unet(self,tomogram):
	# 	m=tomogram.numpy()
	# 	print(m.shape)
	# 	model = self.load_model('./neural_nets/weights_temp.h5')
	# 	p=model.predict(m[None, :, :, None]/3.,verbose=1)
	# 	cout=from_numpy(p[0,:,:,0])
	# 	return cout

	def apply_unet(self,weights_in,tomogram=None,outfile=None):
		if tomogram==None:
			print("Need to specify tomogram to apply U-net")
			return

		# if weights_in==None:
		# 	print("Need to provide weights of U-net")
		# 	return
		if not os.path.exists("./neural_nets/weights_temp.h5"):
			print("There's no neural networks trained for this tomogram. Train a neural net first")
			return
		else:
			pass
		self.model = self.load_model('./neural_nets/weights_temp.h5')
		# nframe=EMUtil.get_image_count(tomogram)
		# is3d=False
		# ### deal with 3D volume or image stack
		# e=EMData(tomogram, 0, True)
		#
		# apix=e["apix_x"]
		# if nframe==1:
		# 	nframe=e["nz"]
		# 	if nframe>1:
		# 	#### input data is 3D volume
		# 		is3d=True

		is3d=True
		enx,eny=tomogram["nx"], tomogram["ny"]
		nframe=tomogram["nz"]
		tsz=max(enx,eny)
		output=EMData(enx, eny, nframe)


		print("Loading tomogram...")
		tomo_in=[]
		for nf in range(nframe):
			e0=tomogram.get_clip(Region((enx-tsz)//2,(eny-tsz)//2,nf,tsz,tsz,1))
			tomo_in.append(e0)
		print(len(tomo_in))
		for idx, img in enumerate(tomo_in):
			# if idx == 123:
			# 	plt.imshow(tomo_in[idx].numpy())
			m=img.numpy()
			p=self.model.predict(m[None, :, :, None]/3.,verbose=1)
			#p[p<0]=0
			cout=from_numpy(p[0,:,:,0])
			cout=cout.get_clip(Region((cout["nx"]-enx)//2,(cout["ny"]-eny)//2 ,enx, eny))
			#cout.scale(int(options.labelshrink))
			output.insert_clip(cout, [0,0,idx])

			sys.stdout.write("\r  {}/{} finished.".format(idx+1, len(tomo_in)))
			sys.stdout.flush()
		if outfile:
			output.write_image(outfile)
		return output

	def normalize(self,data):
		mean = np.mean(data)
		std = np.std(data)
		data = data - mean
		data = data/std
		return data
	def stack_outfile(self,*dats, outfile):
		im_x,im_y = dats[0].shape[1:]
		out_array = np.array([*zip(*dats)]).reshape(-1,im_x,im_y)
		#print(out_array.shape)
		out = from_numpy(out_array)
		out.write_image(outfile)

class NNet_Tab(QtWidgets.QWidget):
	def __init__(self,target):
		QtWidgets.QWidget.__init__(self,None)
		#self.target=weakref.ref(target)
		self.target = target

		self.bg_button = QtWidgets.QPushButton("Background")
		self.ann_button = QtWidgets.QPushButton("Classes")
		self.bg_button.setCheckable(True)
		#self.bg_button.setFixedWidth(100)
		self.ann_button.setCheckable(True)
		#self.bg_button.setFixedWidth(100)

		self.nn_cb_group = QtWidgets.QButtonGroup()
		self.nn_cb_group.addButton(self.bg_button,0)
		self.nn_cb_group.addButton(self.ann_button,1)



		self.train_class_button=QtWidgets.QPushButton("Train NNet")
		#self.train_class_button.setFixedWidth(70)
		self.train_no_iters_sb = StringBox(label="No iters",value="8",showenable=-1)
		#self.train_no_iters_sb.text.setFixedWidth(60)
		self.train_lr_sb = StringBox(label="Learnrate",value="3e-4",showenable=-1)
		#self.train_lr_sb.text.setFixedWidth(60)
		self.build_ts_button = QtWidgets.QPushButton("Build Trainset")
		#self.build_ts_button.setFixedWidth(120)
		self.build_class_sb=StringBox(label="Class ",value="1",showenable=-1)
		#self.build_class_sb.text.setFixedWidth(60)
		self.build_norep_sb=StringBox(label="No reps",value="5",showenable=-1)
		#self.build_norep_sb.text.setFixedWidth(60)
		#self.train_all_button=QtWidgets.QPushButton("Train NNet for all ")
		self.apply_button=QtWidgets.QPushButton("Apply NNet")
		self.apply_all_button=QtWidgets.QPushButton("Apply All")


		nnet_gbl=QtWidgets.QGridLayout(self)
		nnet_gbl.setColumnStretch(70,70)
		nnet_gbl.addWidget(self.bg_button,0,0,1,1)
		nnet_gbl.addWidget(self.ann_button,0,1,1,1)
		nnet_gbl.addWidget(self.build_ts_button,1,0,1,1)
		nnet_gbl.addWidget(self.build_class_sb,2,0,1,1)
		nnet_gbl.addWidget(self.build_norep_sb,3,0,1,1)

		nnet_gbl.addWidget(self.train_class_button,1,1,1,1)
		nnet_gbl.addWidget(self.train_no_iters_sb,2,1,1,1)
		nnet_gbl.addWidget(self.train_lr_sb,3,1,1,1)

		#nnet_gbl.addWidget(QtWidgets.QLabel("Apply to tomogram"),3,0,1,1)
		nnet_gbl.addWidget(self.apply_button,4,0,1,1)
		nnet_gbl.addWidget(self.apply_all_button,4,1,1,1)

		self.bg_button.clicked[bool].connect(self.bg_bt_clicked)
		self.ann_button.clicked[bool].connect(self.ann_bt_clicked)
		self.nn_cb_group.buttonClicked[QtWidgets.QAbstractButton].connect(self.on_check_nn_cb_group)
		self.train_class_button.clicked[bool].connect(self.train_class_bt_clicked)
				#self.train_all_button.clicked[bool].connect(self.train_all_bt_clicked)
		self.apply_button.clicked[bool].connect(self.apply_bt_clicked)
		self.build_ts_button.clicked[bool].connect(self.build_trainset)


	# def extract_bt_clicked(self):
	# 	print("Extract image patches of size",self.bsz_vs.value,"at positions:",self.boxes)
	# 	outfile = "./particles/bg_temp.hdf"
	# 	try:
	# 		os.remove(outfile)
	# 	except:
	# 		pass
	# 	for i in range(len(self.boxes)):
	# 		x,y,z = int(self.boxes[i][0]),int(self.boxes[i][1]),int(self.boxes[i][2])
	# 		bs = self.bsz_vs.value
	# 		if self.boxes[i][3] == 1:
	# 			r = self.target.data.get_clip(Region(x-bs//2,y-bs//2,z,bs,bs,1))
	# 			r.write_image(outfile,-1)
	# 			#r=self.data.get_clip(Region(x-bs//2,y-bs//2,z-bz//2,bs,bs,bz))
	# 	self.clear_shapes()

	def extract_region(self, iter=3,thresh=0.5):
		#enumerate label_stack by objects then find bounding box for each object
		reg_list = []
		#datas = to_numpy(self.get_data())
		labels = to_numpy(self.target.get_annotation())

		if len(labels.shape) == 2:#single 2D image
			#datas = np.expand_dims(datas, axis=0)
			labels = np.expand_dims(labels, axis=0)
		#print("Data shape, label shape", datas.shape,labels.shape)

		for i in range(len(labels)):
			open_lab=ndi.binary_opening(labels[i],iterations=iter)
			labeled,num = ndi.label(open_lab>thresh)
			cent_mass = ndi.center_of_mass(open_lab,labeled,[i+1 for i in range(num)])
			# open_lab = morphology.opening(labels[i],disk(3))
			# labeled= morphology.label(open_lab>thres,connectivity=2,return_num=False)
			# temp_im = datas[i]
			# regions = regionprops(labeled)
			# for region in regions:
			# 	#minr, minc, maxr, maxc = region.bbox
			# 	y,x = region.centroid
			# reg_list.append([x,y,i for (x,y) in cent_mass])
			print("mass at slice", i,cent_mass)
			for pair in cent_mass:
				reg_list.append([pair[1],pair[0],i])

		print(reg_list)
		return reg_list

	def create_organelles_training_set(self):
		self.reg_list = self.extract_region()
		datas = self.target.get_data()
		labels = self.target.get_annotation()

		d_outfile = "./particles/org_temp.hdf"
		l_outfile = "./particles/org_temp_seg.hdf"
		#print("Create image patches of size",self.bsz_vs.value,"at positions:",self.reg_list)
		try:
			os.remove(d_outfile)
			os.remove(l_outfile)
		except:
			pass
		bs = self.target.bsz_vs.value
		for i in range(len(self.reg_list)):
			x = int(self.reg_list[i][0])
			y = int(self.reg_list[i][1])
			z = int(self.reg_list[i][2])
			d = datas.get_clip(Region(x-bs//2,y-bs//2,z,bs,bs,1))
			l = labels.get_clip(Region(x-bs//2,y-bs//2,z,bs,bs,1))
			d.write_image(d_outfile,-1)
			l.write_image(l_outfile,-1)

	def build_trainset(self):
		print("Building trainset")
		self.create_organelles_training_set()
		os.system("e2tomoseg_buildtrainset.py --buildset --particles_raw=./particles/org_temp.hdf --particles_label=./particles/org_temp_seg.hdf --boxes_negative=./particles/bg_temp.hdf --ncopy={} --trainset_output=./particles/org_temp_trainset.hdf --validset=0.0".format(int(self.build_norep_sb.getValue())))


	def bg_bt_clicked(self):
		self.target.basic_tab.setCurrentIndex(3)
		print("Select background patches for training")
		return

	def ann_bt_clicked(self):
		self.target.basic_tab.setCurrentIndex(0)
		print("Annotate objects for training")
		self.target.img_view.show_inspector(6)
		return


	def train_class_bt_clicked(self):
		annotation_out=self.target.img_view.get_full_annotation().copy_head()
		annotation_out.to_zero()
		trainset_path = './particles/org_temp_trainset.hdf'
		if not os.path.exists(trainset_path):
			print("No trainset file detected in the ./particles folder. Abort")
		else:
			pass
		self.unet = UNet(infile="./particles/org_temp_trainset.hdf")
		print("Initialize Unet")
		try:
			os.listdir('./neural_nets')
		except:
			os.mkdir('./neural_nets')
		unet_wout = os.path.join("./neural_nets",self.target.tomogram_list.currentItem().text()[:-4]+"_nnet.h5")

		self.unet.train_unet(weights_out=unet_wout,batch_sz=50,val_split=0.2,no_epoch = int(self.train_no_iters_sb.getValue()),learnrate=float(self.train_lr_sb.getValue()))
		#print("Done training Unet. Network weights saved to", self.unet_wout)


	def train_all_bt_clicked(self):

		return


	def apply_bt_clicked(self):

		if not self.unet:
			self.unet = UNet()
		unet_win = os.path.join("./neural_nets",self.target.tomogram_list.currentItem().text()[:-4]+"_nnet.h5")
		if not os.path.exists(unet_win):
			print("No neural networks saved for this tomogram. Train a Unet first.")
			return
		pred_map = self.unet.apply_unet(weights_in=unet_win, tomogram=self.target.img_view.get_full_data())
		print("Done applying unet")
		self.target.img_view.full_annotation =pred_map
		self.target.img_view.force_display_update(set_clip=0)
		self.target.img_view.updateGL()
		return

	def on_check_nn_cb_group(self,cb):
		print(cb.text()+" is selected")
		if cb.text() == "Background":
			self.target.img_view.mouse_mode=1
		elif cb.text() =="Classes":
			self.target.img_view.mouse_mode=6
			self.target.img_view.show_inspector(6)
		else:
			return

class Morp_Tab(QtWidgets.QWidget):
	def __init__(self,target) :
		QtWidgets.QWidget.__init__(self,None)
		self.target = target
		self.morp_n_iters_sp = QtWidgets.QSpinBox()
		self.morp_n_iters_sp.setValue(3)
		self.morp_close_bt = QtWidgets.QPushButton("Closing")
		self.morp_open_bt = QtWidgets.QPushButton("Opening")
		self.morp_erode_bt = QtWidgets.QPushButton("Erosion")
		self.morp_dilate_bt = QtWidgets.QPushButton("Dilation")
		self.morp_label_bt = QtWidgets.QPushButton("Numbering")

		morp_gbl=QtWidgets.QGridLayout(self)
		#morp_gbl.setColumnStretch(70,70)
		morp_gbl.addWidget(QtWidgets.QLabel("No iters"), 0,0,1,1)
		morp_gbl.addWidget(self.morp_n_iters_sp, 0,1,1,1)
		morp_gbl.addWidget(self.morp_close_bt, 1,0,1,1)
		morp_gbl.addWidget(self.morp_open_bt, 1,1,1,1)
		morp_gbl.addWidget(self.morp_erode_bt, 2,0,1,1)
		morp_gbl.addWidget(self.morp_dilate_bt, 2,1,1,1)
		morp_gbl.addWidget(self.morp_label_bt, 3,0,1,1)

		self.morp_close_bt.clicked[bool].connect(self.do_morp_close)
		self.morp_open_bt.clicked[bool].connect(self.do_morp_open)
		self.morp_erode_bt.clicked[bool].connect(self.do_morp_erode)
		self.morp_dilate_bt.clicked[bool].connect(self.do_morp_dilate)
		self.morp_label_bt.clicked[bool].connect(self.do_morp_label)

	def do_morp_close(self):

		n_iters = int(self.morp_n_iters_sp.value())
		mask = self.target.get_annotation()
		self.target.annotate = from_numpy(ndi.binary_closing(to_numpy(mask),iterations=n_iters))
		self.target.img_view.set_data(self.target.data, self.target.annotate)
		del mask

	def do_morp_open(self):
		n_iters = int(self.morp_n_iters_sp.value())
		mask = self.target.get_annotation()
		self.target.annotate = from_numpy(ndi.binary_opening(to_numpy(mask),iterations=n_iters))
		self.target.img_view.set_data(self.target.data, self.target.annotate)
		del mask

	def do_morp_dilate(self):
		n_iters = int(self.morp_n_iters_sp.value())
		mask = self.target.get_annotation()
		self.target.annotate = from_numpy(ndi.binary_dilation(to_numpy(mask),iterations=n_iters))
		self.target.img_view.set_data(self.target.data, self.target.annotate)
		del mask



	def do_morp_erode(self):
		n_iters = int(self.morp_n_iters_sp.value())
		mask = self.target.get_annotation()
		self.target.annotate = from_numpy(ndi.binary_erosion(to_numpy(mask),iterations=n_iters))
		self.target.img_view.set_data(self.target.data, self.target.annotate)
		del mask

	def do_morp_label(self):
		n_iters = int(self.morp_n_iters_sp.value())
		mask, num = ndi.label(to_numpy(self.target.get_annotation()))
		self.target.annotate = from_numpy(mask)
		print("number of object detected:", num)
		self.target.img_view.set_data(self.target.data, self.target.annotate)
		del mask

class Binary_Tab(QtWidgets.QWidget):
	def __init__(self,target) :
		QtWidgets.QWidget.__init__(self,None)
		self.target = target

		#Set up Binary Tab and Function
		bin_gbl = QtWidgets.QGridLayout(self)
		self.bin_invert_cb = QtWidgets.QCheckBox("Dark Feature")
		self.bin_invert_cb.setChecked(True)
		#valslider set to 0 will paint all the annotate file
		self.bin_detect_bt = QtWidgets.QPushButton("Detect Feature")
		self.bin_fill_bt = QtWidgets.QPushButton("Fill")
		self.bin_trim_bt = QtWidgets.QPushButton("Trim")
		self.bin_low_pass_vs = ValSlider(value=1,rng=(0.001,1),rounding=2,label= "Cut-off Abs")
		self.bin_threshold_vs = ValSlider(value=0.001,rng=(0.001,1),rounding=2,label="Threshold  ")
		self.closing_n_iters =1
		self.opening_n_iters =1

		bin_gbl.addWidget(self.bin_invert_cb,0,0,1,1)
		bin_gbl.addWidget(self.bin_detect_bt,0,1,1,1)
		bin_gbl.addWidget(self.bin_low_pass_vs,1,0,1,2)
		bin_gbl.addWidget(self.bin_threshold_vs,2,0,1,2)
		bin_gbl.addWidget(self.bin_fill_bt,3,0,1,1)
		bin_gbl.addWidget(self.bin_trim_bt,3,1,1,1)

		#self.setLayout(bin_gbl)
		self.bin_detect_bt.clicked[bool].connect(self.bin_detect_bt_click)
		self.bin_fill_bt.clicked[bool].connect(self.do_area_closing)
		self.bin_trim_bt.clicked[bool].connect(self.do_area_opening)
		self.bin_low_pass_vs.valueChanged.connect(self.update_mask_from_vs)
		self.bin_threshold_vs.valueChanged.connect(self.update_mask_from_vs)


	def do_area_opening(self):
		self.target.morp_tab.morp_n_iters_sp.setValue(self.opening_n_iters)
		self.target.morp_tab.do_morp_open()
		self.opening_n_iters +=1

	def do_area_closing(self):
		self.target.morp_tab.morp_n_iters_sp.setValue(self.closing_n_iters)
		self.target.morp_tab.do_morp_close()
		self.closing_n_iters +=1


	def bin_detect_bt_click(self):
		if self.bin_low_pass_vs.value == 0.1:
			self.bin_low_pass_vs.setValue(0.101)
			self.bin_threshold_vs.setValue(0.601)
		self.bin_low_pass_vs.setValue(0.1)
		self.bin_threshold_vs.setValue(0.6)
		self.target.reset_morp_params()

	def update_mask_from_vs(self):
		if self.bin_invert_cb.isChecked():
			mult = -1
		else:
			mult= 1
		lp = self.bin_low_pass_vs.value
		thres = self.bin_threshold_vs.value
		#self.mask = self.get_annotation().process("threshold.binary",{"value":0.3})
		self.mask = self.target.get_annotation()
		self.lp_masked = mult*self.mask*self.target.get_data().process("filter.lowpass.gauss",{"cutoff_abs":lp})
		self.thres_mask = self.lp_masked.process("threshold.binary",{"value":thres})
		self.target.img_view.set_data(self.target.get_data(),self.thres_mask*self.mask)
		return


class Specific_Tab(QtWidgets.QWidget):
	def __init__(self,target) :
		QtWidgets.QWidget.__init__(self,None)
		self.target = target

		cell_vbl = QtWidgets.QVBoxLayout(self)
		cell_button_l = QtWidgets.QVBoxLayout()

		self.memb_tab = QtWidgets.QWidget()
		self.mtlay = QtWidgets.QVBoxLayout(self.memb_tab)
		self.mem_detect_button = QtWidgets.QPushButton("Auto detect in brushed region")
		self.bezier_button = QtWidgets.QPushButton("Draw bezier curve")
		self.double_mem_button = QtWidgets.QPushButton("Double membrane")
		self.mtlay.addWidget(self.mem_detect_button)
		self.mtlay.addWidget(self.bezier_button)
		self.mtlay.addWidget(self.double_mem_button)

		self.actin_tab = QtWidgets.QWidget()
		self.atlay = QtWidgets.QGridLayout(self.actin_tab)
		self.actin_button = QtWidgets.QPushButton("Actin")
		self.atlay.addWidget(self.actin_button)


		self.gra_tab = QtWidgets.QWidget()
		self.gtlay = QtWidgets.QVBoxLayout(self.gra_tab)
		self.dark_gra_button = QtWidgets.QPushButton("Dark granule")
		self.alpha_gra_button = QtWidgets.QPushButton("Alpha granule")
		self.ves_button = QtWidgets.QPushButton("Vesicle")
		self.gtlay.addWidget(self.dark_gra_button)
		self.gtlay.addWidget(self.alpha_gra_button)
		self.gtlay.addWidget(self.ves_button)

		self.cell_tab = QtWidgets.QTabWidget()
		self.cell_tab.addTab(self.memb_tab,"Membrane")
		self.cell_tab.addTab(self.actin_tab,"Actin/Microtubules")
		self.cell_tab.addTab(self.gra_tab,"Granule")
		cell_button_l.addWidget(self.cell_tab)
		cell_vbl.addLayout(cell_button_l)


if __name__ == '__main__':
	main()
