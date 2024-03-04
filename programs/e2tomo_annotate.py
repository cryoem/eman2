#!/usr/bin/env python
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
from scipy.spatial import ConvexHull
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
	parser.add_argument("--tomos",type=str,help="The tomograms to annotate",default="")
	parser.add_argument("--folder",type=str, help="List the folder contain all tomograms to process", default="./tomograms")

	parser.add_argument("--seg_folder",type=str, help="List the folder contain all annotation file", default="./segs/")
	parser.add_argument("--region_sz",type=int, help="Region size for Region I/O. -1 reads whole tomogram", default=-1)
	parser.add_argument("--zthick",type=int, help="Zthickness for Region I/O. -1 reads full z. Default=0 reading a single central slice.", default=0)
	parser.add_argument("--no_tmp",action="store_true", help="DO NOT saving temporary file for undo function. Make the program faster when handle big volume. Default=False.", default=False)
	parser.add_argument("--load_data_czi",action="store_true", help="Loading data from CZI cryoet_data_portal to open in an e2tomo_annotate project", default=False)

	(options, args) = parser.parse_args()

	app = EMApp()
	if options.load_data_czi:
		os.system("python ./e2import_from_czi.py")
		return
	else:
		pass

	#
	# if len(args) == 0:
	# 	print("INPUT ERROR: You must specify an image to display.")
	# 	sys.exit(1)

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
		#System Menu
		self.mfile=self.menuBar().addMenu("File")
		self.mfile_new_full=self.mfile.addAction("New Full Annotation")
		self.mfile_save_full=self.mfile.addAction("Save Full Annotation")
		self.mfile_del_all=self.mfile.addAction("Delete Full Annotation")

		self.mfile_quit=self.mfile.addAction("Close All Windows")

		self.medit=self.menuBar().addMenu("Edit")
		self.medit_undo=self.medit.addAction("Undo")
		self.medit_undo.setEnabled(False)

		self.mfile_del_all.triggered[bool].connect(self.menu_file_del_all)
		self.mfile_quit.triggered[bool].connect(self.menu_file_quit)
		self.mfile_new_full.triggered[bool].connect(self.menu_file_new_seg)
		self.medit_undo.triggered[bool].connect(self.menu_edit_undo)


		self.tom_folder = options.folder
		self.seg_folder = options.seg_folder
		self.no_tmp = options.no_tmp

		try:
			os.mkdir(self.seg_folder)
			print("Directory", self.seg_folder, "is created" )
		except OSError as error:
			print("Directory",self.seg_folder,"already existed. Continue")
			pass

		self.tomogram_list = QtWidgets.QListWidget()
		self.tom_file_list = []
		self.seg_grps = []

		self.tomo_tree = QtWidgets.QTreeWidget()

		if len(self.tom_folder)>0:
			tom_ls = sorted(os.listdir(self.tom_folder))
			# seg_ls = sorted(os.listdir(self.seg_folder))
			for file_name in tom_ls:
				if file_name.endswith(".hdf"):
					self.tom_file_list.append(file_name)

		elif len(options.tomos) == 0:
			print("Specify tomogram or tomograms folder for start annotation")
			sys.exit(0)
			return

		else:
			tom_ls = options.tomos.split(',')
			for file_name in tom_ls:
				if file_name.endswith(".hdf"):
					self.tom_file_list.append(file_name)
				else:
					print(file_name,"is not a valid file_name. Pass")


		print("Self.tom_file_list",self.tom_file_list)
		for t in range(len(self.tom_file_list)):
			file_name = self.tom_file_list[t]
			#Populate the Tree
			tomo_tree_item = QtWidgets.QTreeWidgetItem([file_name])
			self.tomo_tree.addTopLevelItem(tomo_tree_item)
			self.seg_grps.append(QtWidgets.QButtonGroup())
			i = -1
			seg_ls = sorted(os.listdir(self.seg_folder))
			for seg_f in seg_ls:
				# if seg_f.startswith(file_name[0:-4]) and seg_f.endswith(".hdf"):
				data_file = str(os.path.join(self.tom_folder,file_name))
				data_hdr = EMData(data_file, 0,True)
				if seg_f.startswith(base_name(data_file)) and seg_f.endswith(".hdf"):
					i +=1
					hdr = EMData(os.path.join(self.seg_folder,seg_f), 0,True)
					if [hdr["nx"],hdr["ny"],hdr["nz"]]!= [data_hdr["nx"],data_hdr["ny"],data_hdr["nz"]]:
						print("Mismatch dimension between full data and annotation map ",seg_f,". Skip")
						continue
					else:
						pass
					seg_item = QtWidgets.QTreeWidgetItem()
					tomo_tree_item.addChild(seg_item)
					seg_item.setFlags(Qt.ItemFlags(Qt.ItemIsEnabled)|Qt.ItemFlags(Qt.ItemIsUserCheckable))
					#seg_item.setCheckState(0,0)
					seg_button = QtWidgets.QRadioButton(seg_f)
					self.seg_grps[t].addButton(seg_button,i)
					self.tomo_tree.setItemWidget(seg_item,0,seg_button)

			try:
				self.seg_grps[t].button(0).setChecked(True)
			except:#No seg file available for some tomograms
				# seg_path = file_name[0:-4]+"_seg.hdf"
				print("No seg files available for tomogram", file_name)
				pass

		file_name = self.tom_file_list[0]
		self.data_file = str(os.path.join(self.tom_folder,file_name))
		hdr=EMData(self.data_file, 0,True)

		self.nx = hdr["nx"]
		self.ny = hdr["ny"]
		self.nz=hdr["nz"]
		self.apix = hdr["apix_x"]

		if self.tomo_tree.topLevelItem(0).childCount() > 0:
			seg_path = os.path.join(self.seg_folder,self.seg_grps[0].checkedButton().text())
			print("Seg file for current tomogram already exists at",seg_path)

		else:
			seg_path = os.path.join(self.seg_folder,base_name(self.data_file)+"_seg.hdf")
			print("Create seg_file",seg_path)
			seg_out = EMData(self.nx,self.ny,self.nz)
			seg_out.write_image(seg_path)
			del seg_out
			seg_item = QtWidgets.QTreeWidgetItem()
			try:
				self.tomo_tree.currentItem().addChild(seg_item)
				seg_item.setFlags(Qt.ItemFlags(Qt.ItemIsEnabled))

				#seg_item.setCheckState(0,0)
				seg_button = QtWidgets.QRadioButton(os.path.basename(seg_path))
				self.seg_grps[self.tomo_tree.indexOfTopLevelItem(self.tomo_tree.currentItem())].addButton(seg_button,0)
				seg_button.setChecked(True)
				self.tomo_tree.setItemWidget(seg_item,0,seg_button)
			except:
				pass

		self.seg_path = seg_path
		self.seg_info_path = self.seg_path[0:-4]+".json"
		try:
			f = open(self.seg_info_path, 'x')
		except:
			print("Info file for the first one already exists. Continue ")
			# if t == 0: #FIRST ITEM OF THE TREE
			# 	self.data_file = str(os.path.join(self.tom_folder,file_name))
			# 	hdr=EMData(self.data_file, 0,True)
			#
			# 	self.nx = hdr["nx"]
			# 	self.ny = hdr["ny"]
			# 	self.nz=hdr["nz"]
			# 	self.apix = hdr["apix_x"]
			#
			# 	self.seg_path = os.path.join(self.seg_folder,base_name(self.data_file)+"_seg.hdf")
			# 	print("First seg path", self.seg_path)
			#
			# 	#if not os.path.isfile(self.seg_path):
			# 	if self.tomo_tree.topLevelItem(0):
			# 		seg_out = EMData(self.nx,self.ny,self.nz)
			# 		#self.write_header(seg_out)
			# 		seg_out.write_image(self.seg_path)
			# 		del seg_out
			# 	else:
			# 		print("Seg file for the first one already exists. Continue ")
			# 		pass
			#
			# 	self.seg_info_path = self.seg_path[0:-4]+".json"
			# 	try:
			# 		f = open(self.seg_info_path, 'x')
			# 	except:
			# 		print("Info file for the first one already exists. Continue ")



		self.tomo_tree.setCurrentItem(self.tomo_tree.topLevelItem(0))
		self.tomo_tree.setHeaderHidden(True)
		self.popwidgets = []



		#TODO
		#Need to copy header to annotation file
		#Need to call garbage collector to remove trash from buffer memory

		print("File_name", self.data_file)
		self.data = EMData(self.data_file)/3.0
		self.apix=self.data['apix_x']

		try:
			self.annotation = EMData(annotate)
		except:
			self.annotation = None
			pass

		self.tomo_tree.currentItemChanged[QtWidgets.QTreeWidgetItem,QtWidgets.QTreeWidgetItem].connect(self.tomotree_current_item_change)
		for seg_grp in self.seg_grps:
			seg_grp.buttonClicked[QtWidgets.QAbstractButton].connect(self.seg_path_change)

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
		self.img_view.setSizePolicy(QtWidgets.QSizePolicy.Preferred,QtWidgets.QSizePolicy.Expanding)
		self.img_view.setMinimumSize(680, 680)

		#self.img_view.show()


		#boxes for Boxer tab
		self.boxes = []
		#2 anchor points for autodetect plane
		self.p1p2 = []
		self.fill_type = None
		self.unet = None
		self.setcolors = ['blue','green','red','cyan','purple','orange','yellow','hotpink','gold']

		self.zt_spinbox = QtWidgets.QSpinBox()

		self.zt_spinbox.setMinimum(-1)
		self.zt_spinbox.setMaximum(self.get_nz()//2)
		self.zt_spinbox.setValue(self.options.zthick)
		self.zc_spinbox = QtWidgets.QSpinBox()
		self.zc_spinbox.setMaximum(self.get_nz())
		self.zc_spinbox.setValue(self.get_nz()//2)
		print("ZCS",self.zc_spinbox.value())
		self.zc_spinbox.setMinimum(0)



		self.thumbnail = Thumbnail(current_file=self.data_file,target=self.img_view,app_target=self,tn_size=self.thumbnail_size)
		self.thumbnail.setSizePolicy(QtWidgets.QSizePolicy.Preferred,QtWidgets.QSizePolicy.Expanding)

		self.thumbnail.setMinimumSize(self.thumbnail_size, self.thumbnail_size)
		try:
			self.thumbnail.set_im()
		except:
			pass

		self.img_view_inspector = self.img_view.get_inspector()
		self.img_view_inspector.setSizePolicy(QtWidgets.QSizePolicy.Preferred,QtWidgets.QSizePolicy.Expanding)
		self.update_tree_set()

		self.tomo_list_panel=QtWidgets.QWidget()
		tomo_vbl = QtWidgets.QGridLayout()
		tomo_vbl.addWidget(QtWidgets.QLabel("Tomograms"),0,0)
		#tomo_vbl.addWidget(self.tomogram_list,1,0)
		tomo_vbl.addWidget(self.tomo_tree,1,0)
		#tomo_vbl.addWidget(self.tomo_tree,4,0)
		tomo_vbl.setRowStretch(1,2)
		tomo_vbl.setRowStretch(2,1)
		tomo_vbl.setRowStretch(3,1)
		tomo_vbl.addWidget(self.thumbnail,2,0)
		zt_hbl = QtWidgets.QHBoxLayout()



		zt_hbl.addWidget(QtWidgets.QLabel("cen"))
		zt_hbl.addWidget(self.zc_spinbox)
		zt_hbl.addWidget(QtWidgets.QLabel("z-thick"))
		zt_hbl.addWidget(self.zt_spinbox)




		tomo_vbl.addLayout(zt_hbl,3,0)
		self.tomo_list_panel.setLayout(tomo_vbl)
		self.tomo_list_panel.setWindowTitle("Tomograms List")
		self.tomo_list_panel.show()


		filter_vbl = QtWidgets.QVBoxLayout()
		self.ft_vbl=QtWidgets.QVBoxLayout()
		self.ft_vbl.addWidget(QtWidgets.QLabel("Filters"))

		filter_vbl.addLayout(self.ft_vbl)
		self.procbox1=StringBox(label="Process1:",value="filter.lowpass.gauss:cutoff_abs=0.125",showenable=0)
		self.procbox1.setMinimumWidth(350)
		self.ft_vbl.addWidget(self.procbox1)

		self.procbox2=StringBox(label="Process2:",value="filter.highpass.gauss:cutoff_pixels=3",showenable=0)
		self.procbox2.setMinimumWidth(350)
		self.ft_vbl.addWidget(self.procbox2)

		self.procbox3=StringBox(label="Process3:",value="math.linear:scale=5:shift=0",showenable=0)
		self.procbox3.setMinimumWidth(350)
		self.ft_vbl.addWidget(self.procbox3)

		self.proclbl1=QtWidgets.QLabel("\t\tImage unchanged, display only!")
		self.ft_vbl.addWidget(self.proclbl1)

		self.procbox1.enableChanged.connect(self.do_filters)
		self.procbox1.textChanged.connect(self.do_filters)
		self.procbox2.enableChanged.connect(self.do_filters)
		self.procbox2.textChanged.connect(self.do_filters)
		self.procbox3.enableChanged.connect(self.do_filters)
		self.procbox3.textChanged.connect(self.do_filters)


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

		self.points_label = QtWidgets.QLabel("Choose anchor points. \nCtrl+Click: start a curve.\nShift+Click: delete a point.")
		self.clear_button = QtWidgets.QPushButton("Clear points")
		self.interp_button = QtWidgets.QPushButton("Interpolate")
		self.tx_line_width=QtWidgets.QSpinBox(self)
		self.tx_line_width.setValue(5)

		self.ltlay.addWidget(self.points_label,0,0,2,2)
		self.ltlay.addWidget(QtWidgets.QLabel("Line Width"),2,2,1,2)
		self.ltlay.addWidget(self.tx_line_width,2,3,1,1)
		#self.ltlay.addWidget(QtWidgets.QLabel("Ann Class"),3,2,1,2)
		#self.ltlay.addWidget(self.tx_ann_class,3,3,1,1)
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

		self.ctlay.addWidget(self.draw_contour_label,0,0,2,2)
		self.ctlay.addWidget(QtWidgets.QLabel("Draw line to annotate file"),2,0,1,2)
		self.ctlay.addWidget(QtWidgets.QLabel("Line Width"),2,2,1,1)
		self.ctlay.addWidget(self.ct_line_width,2,3,1,1)
		# self.ctlay.addWidget(QtWidgets.QLabel("Ann Class"),3,2,1,1)
		# self.ctlay.addWidget(self.ct_ann_class,3,3,1,1)
		self.ctlay.addWidget(self.clear_contour_button,0,2,1,1)
		self.ctlay.addWidget(self.fill_contour_checkbox,4,2,1,1)

		self.boxer_tab = QtWidgets.QWidget()
		self.bxlay = QtWidgets.QGridLayout(self.boxer_tab)
		self.bxlay.setColumnStretch(70,70)


		self.random_bx_bt = QtWidgets.QPushButton("Create Random Boxes")
		self.random_sparse_checkbox = QtWidgets.QCheckBox("Sparse")
		self.random_sparse_checkbox.setChecked(True)

		self.clear_bx_bt = QtWidgets.QPushButton("Clear Boxes")
		self.extract_bt = QtWidgets.QPushButton("Extract Boxes")
		self.current_box_index_sp = QtWidgets.QSpinBox()
		self.current_box_index_sp.setValue(0)
		self.current_box_index_sp.setMinimum(0)



		self.random_bx_sb = StringBox(label="No Box:",value="20",showenable=-1)
		self.bsz_vs = ValSlider(self,value=64,rng=(0,300),rounding=0,label="Box Size")
		self.bsz_vs.setIntonly(1)


		self.bxlay.addWidget(QtWidgets.QLabel("Select regions \nfor training nnet"),0,0,1,1)
		hlay = QtWidgets.QHBoxLayout()
		hlay.addWidget(QtWidgets.QLabel("Current index"))
		hlay.addWidget(self.current_box_index_sp)
		self.bxlay.addLayout(hlay,0,1,1,1)
		self.bxlay.addWidget(self.bsz_vs,1,0,1,2)
		self.bxlay.addWidget(self.random_bx_bt,2,0,1,1)
		self.bxlay.addWidget(self.random_bx_sb,2,1,1,1)
		self.bxlay.addWidget(self.random_sparse_checkbox,2,2,1,1)
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

		#self.basic_tab.tabBarClicked[int].connect(self.basic_tab_change)
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
		self.ext_tab = Extrapolate_Tab(target=self)
		self.nn_tab = NNet_Tab(target=self)
		self.snn_tab = Simple_NNet_Tab(target=self)
		self.morp_tab = Morp_Tab(target=self)
		self.binary_tab = Binary_Tab(target=self)
		self.spec_tab = Specific_Tab(target=self)
		self.templ_tab = Templ_Match_Tab(target=self)
		self.stat_tab = Statistics_Tab(target=self)
		self.subtom_tab = Subtom_Tab(target=self)


		self.assisted_tab.addTab(self.ext_tab,"Extrapolate")
		self.assisted_tab.addTab(self.binary_tab,"AutoDetect")
		self.assisted_tab.addTab(self.nn_tab,"UNet")
		self.assisted_tab.addTab(self.snn_tab,"EMAN2NNet")
		self.assisted_tab.addTab(self.templ_tab,"Template")
		self.assisted_tab.addTab(self.stat_tab,"Statistics")
		self.assisted_tab.addTab(self.morp_tab,"Morphological")
		self.assisted_tab.addTab(self.spec_tab,"Specific")
		self.assisted_tab.addTab(self.subtom_tab,"SubTomogram")

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
		#self.test_button.setCheckable(True)
		self.button_gbl.addWidget(self.test_button,6,0,1,1)




		inspector_vbl = QtWidgets.QVBoxLayout()
		inspector_vbl.addWidget(QtWidgets.QLabel("Manual Annotate Tools"))
		inspector_vbl.addWidget(self.img_view_inspector)

		#self.gbl.addLayout(tomo_vbl,1,0,1,1)
		self.gbl.addWidget(self.img_view,0,0,1,2)
		self.gbl.addLayout(inspector_vbl,0,2,1,1)
		self.centralWidget.setLayout(self.gbl)

		self.control_panel = QtWidgets.QWidget()
		self.control_panel.setLayout(self.button_gbl)
		self.control_panel.setWindowTitle("Control Panel")
		self.control_panel.show()

		self.test_button.clicked[bool].connect(self.test_drawing_function)
		#self.undo_button.clicked[bool].connect(self.reverse_to_saved_state)

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
		self.img_view.mousedown.connect(self.img_view_mouse_down)
		self.img_view.keypress.connect(self.key_press)
		self.img_view.mousedrag.connect(self.img_view_mouse_drag)
		self.zt_spinbox.valueChanged.connect(self.z_spinbox_changed)
		self.zc_spinbox.valueChanged.connect(self.z_spinbox_changed)
		self.previous_z = self.get_zpos()

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

	def get_full_data_from_file(self):
		return EMData(self.data_file)
	def get_full_annotation_from_file(self):
		return EMData(self.seg_path)
	def get_inspector(self):
		return self.img_view_inspector

	def get_segtab(self):
		return self.img_view_inspector.seg_tab

	def get_treeset(self):
		return self.get_segtab().tree_set

	def get_current_class(self):
		return self.get_segtab().get_current_class()


	def zchange(self,value):
		print(value)


	def do_filters(self):
		#print("do filter")
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
		self.save_current_state()
		self.binary_tab.saved_ori = False
		self.reset_morp_params()


	def reset_morp_params(self,reset_vs=False):
		self.binary_tab.closing_n_iters =1
		self.binary_tab.opening_n_iters =1
		self.morp_tab.morp_n_iters_sp.setValue(1)
		if reset_vs:
			self.binary_tab.bin_low_pass_vs.setValue(1)
			self.binary_tab.bin_threshold_vs.setValue(0.001)




	def seg_path_change(self,seg_bt):
		#print(grp.text()+" is selected")
		self.seg_path = os.path.join(self.seg_folder,seg_bt.text())
		self.seg_info_path = self.seg_path[0:-4]+".json"
		self.reset_morp_params(reset_vs=False)
		self.set_imgview_data(round(self.data_xy[0]),round(self.data_xy[1]),self.img_view_region_size)
		self.update_tree_set()

		#print("Seg file for current tomogram already exists at",seg_path)

	def tomotree_item_change(self,item,int):
		if (item.parent()):

			print(item.text(int), "is changed")
		else:
			pass
		return

	def tomotree_current_item_change(self,current,previous):
		print(current.text(0))
		if (current.parent()):
			print("Current Tomogram",current.parent().text(0))
			data_file = current.parent().text(0)
			self.tomo_tree.setCurrentItem(current.parent(),0)
			return
		else:
			print("Current Tomogram", current.text(0))
			data_file = current.text(0)

		self.data_file = str(os.path.join(self.tom_folder,data_file))
		hdr=EMData(self.data_file, 0,True)
		self.nx=hdr["nx"]
		self.ny=hdr["ny"]
		self.nz=hdr["nz"]
		self.apix = hdr["apix_x"]
		#print('Nz',self.nz)
		self.get_inspector().seg_tab.write_treeset_json(self.seg_info_path)

		#seg_path = os.path.join(self.seg_folder,self.tomogram_list.item(int).text()[0:-4]+"_seg.hdf")
		if current.childCount() > 0:
			seg_path = os.path.join(self.seg_folder,self.seg_grps[self.tomo_tree.indexOfTopLevelItem(self.tomo_tree.currentItem())].checkedButton().text())
			print("Seg file for current tomogram already exists at",seg_path)
		else:
			seg_path = os.path.join(self.seg_folder,base_name(self.data_file)+"_seg.hdf")
			print("Create seg_file",seg_path)
			seg_out = EMData(self.nx,self.ny,self.nz)
			seg_out.write_image(seg_path)
			del seg_out
			seg_item = QtWidgets.QTreeWidgetItem()
			self.tomo_tree.currentItem().addChild(seg_item)
			seg_item.setFlags(Qt.ItemFlags(Qt.ItemIsEnabled))
			#seg_item.setCheckState(0,0)
			seg_button = QtWidgets.QRadioButton(os.path.basename(seg_path))
			self.seg_grps[self.tomo_tree.indexOfTopLevelItem(self.tomo_tree.currentItem())].addButton(seg_button,0)
			seg_button.setChecked(True)
			self.tomo_tree.setItemWidget(seg_item,0,seg_button)



		seg_info_path = os.path.join(seg_path[0:-4]+".json")
		#print("This is seg_info_path",seg_info_path)
		try:
			fp = open(seg_info_path, 'x')
		except:
			print("Info file for the current one already exists. Continue ")

		self.zt_spinbox.setMaximum(self.nz//2)
		#self.data = EMData(self.data_file)
		self.thumbnail.get_im(self.data_file)
		self.thumbnail.set_im()
		self.clear_shapes()

		#TODO - Fix zthick and hard code part

		self.data_xy = self.thumbnail.get_xy()
		#print("data xy", self.data_xy)

		if self.get_annotation():
			print("Print annotation to file", self.seg_path)

			#self.write_header(self.get_annotation())
			self.write_out(self.get_annotation(), self.seg_path, self.cur_region)

		else:
			print("Annotation is none.")
			pass

		self.zc_spinbox.setValue(self.nz//2)
		self.zc_spinbox.setMaximum(self.nz)


		#print("SegPath,newSegPath", self.seg_path,seg_path)
		self.seg_path = seg_path
		self.seg_info_path = seg_info_path
		self.reset_morp_params(reset_vs=False)
		self.set_imgview_data(round(self.data_xy[0]),round(self.data_xy[1]),self.img_view_region_size)
		self.update_tree_set()



	#need to write region out before setting new data
	def write_out(self,file_out, out_name, region):
		#self.set_scale=1
		if file_out:
			try:
				file_out.write_image(out_name, 0, IMAGE_HDF, False, region)
			except Exception as e:
				print("Error writing file out due to",e)
				return
		else:

			return


	def set_imgview_data(self,x,y,sz):
		try:
			cen = int(self.zc_spinbox.value())
		except:
			cen = self.nz//2

		print("Nz ,iz",self.nz,cen)
		print("Img x,y,sz",x,y,sz)
		if self.nz == 1:
			self.zthick = 0
			#self.cur_region = Region(x-sz//2,y-sz//2,sz,sz)
			self.cur_region = Region(max(0,floor(x-sz//2)),max(0,floor(y-sz//2)),sz,sz)

		else:
			try:
				self.zthick = int(self.zt_spinbox.value())
				print("zthick",self.zthick)
			except:
				self.zthick = 0
			if self.zthick == -1:
				#print(self.nz)
				self.zc_spinbox.setValue(self.nz//2)
				self.zc_spinbox.setMaximum(self.nz)
				self.cur_region = Region(max(0,floor(x-sz//2)),max(0,floor(y-sz//2)),0,sz,sz, self.nz)
				#self.data = EMData(self.data_file, 0, False, Region(x-old_div(sz,2),y-old_div(sz,2),0, sz, sz, self.nz))
			else:
				#self.cur_region = Region(x-sz//2,y-sz//2,cen-self.zthick, sz, sz,self.zthick*2+1)
				self.cur_region = Region(max(0,floor(x-sz//2)),max(0,floor(y-sz//2)),cen-self.zthick, sz, sz,self.zthick*2+1)
				print(self.cur_region)
				#self.data = EMData(self.data_file, 0, False, Region(x-old_div(sz,2),y-old_div(sz,2),iz-self.zthick, sz, sz,self.zthick*2+1))
		self.data = EMData(self.data_file, 0, False, self.cur_region)
		# seg_path = os.path.join(self.seg_folder,self.tomogram_list.currentItem().text()[0:-4]+"_seg.hdf")


		self.annotate = EMData(self.seg_path, 0, False, self.cur_region)
		#self.read_metadata(seg_path)

		self.img_view.set_data(self.data, self.annotate)
		self.img_view.set_scale(1)
		self.img_view.set_origin(0,0)
		#print("Imgview, inspector, segtab",self.img_view,self.img_view.get_inspector(),self.img_view.inspector.seg_tab)
		#self.img_view.get_inspector().seg_tab.read_header(seg_path)

	def update_tree_set(self,seg_info_path = ""):
		#self.get_inspector().seg_tab.tree_set.clear()
		if len(seg_info_path) == 0:
			#seg_info_path = os.path.join(self.seg_folder,'info',self.tomogram_list.currentItem().text()[0:-4]+"_seg_info.json")
			#seg_info_path = os.path.join(self.seg_folder,self.tomogram_list.currentItem().text()[0:-4]+"_seg.json")
			seg_info_path = self.seg_info_path
		#if os.path.isdir(seg_info_path):
		try:
			#self.get_treeset().clear()
			self.get_segtab().read_json_treeset(seg_info_path)
		except:
			print("No info file available. Using class info of project")
			return
		self.get_segtab().update_sets()
		# try:
		tree_topit = self.get_segtab().tree_set.topLevelItem(0)
		self.get_segtab().tree_set.setCurrentItem(tree_topit,1)
		# except:
		# 	print("Except")
		# 	pass





	def z_spinbox_changed(self,event):
		self.data_xy = self.thumbnail.get_xy()
		try:
			self.write_out(self.get_annotation(), self.seg_path, self.cur_region)
			#self.get_annotation().write_image(self.seg_path, 0, IMAGE_HDF, False, self.cur_region)
		except:
			print("Zt_change annotation",self.get_annotation())
			print("image cannot be write to disk")#when annotation files is None
			pass
		self.set_imgview_data(self.data_xy[0],self.data_xy[1],self.img_view_region_size)
		self.reset_morp_params(reset_vs=False)


	def test_drawing_function(self):
		self.save_current_state()
		return

	def save_current_state(self):
		if self.no_tmp:
			return
		temp_path = os.path.join(self.seg_folder,"current_temp.hdf")
		temp_info_path =  os.path.join(self.seg_folder,"info","current_set_info.json")
		if not os.path.isfile(temp_path):
			seg_out = EMData(self.nx,self.ny,self.nz)
			#self.write_header(seg_out)
			seg_out.write_image(temp_path)
			del seg_out
		else:
			#print("Seg file for the first one already exists. Continue ")
			pass

		try:
			self.write_out(self.get_annotation(), temp_path, self.cur_region)

			print("Saved current state to disk")
			#self.get_annotation().write_image(self.seg_path, 0, IMAGE_HDF, False, self.cur_region)
		except Exception as e:
			print("image cannot be write to disk due to", e)#when annotation files is None
			pass

		try:

			self.get_inspector().seg_tab.write_treeset_json(temp_info_path)
		except Exception as e:
			print("Cannot print metadata to info file due to", e)
			pass

		#self.undo_button.setEnabled(True)

		self.medit_undo.setEnabled(True)

		return


	def reverse_to_saved_state(self):
		if self.no_tmp:
			print("no_tmp mode does not allow undo function. Re-initiate the program to enable the undo function")
			return
		print("Back to previous state")
		temp_path = os.path.join(self.seg_folder,"current_temp.hdf")
		temp_info_path =  os.path.join(self.seg_folder,"info","current_set_info.json")

		self.annotate = EMData(temp_path, 0, False, self.cur_region)
		#self.read_metadata(seg_path)
		self.img_view.set_data(self.data, self.annotate)
		self.update_tree_set(seg_info_path = temp_info_path)
		#self.undo_button.setEnabled(False)
		self.medit_undo.setEnabled(False)
		return




	def create_trainset_from_boxes(self):
		print("TESTBUTTON CLICKED")
		print("Extract image and label patches of size",self.bsz_vs.value,"at positions:",self.boxes)
		outfile = "training_temp.hdf"
		try:
			os.remove(outfile)
		except:
			# print("File not Exiting")
			pass

		for i in range(len(self.boxes)):
			nz = self.data.get_sizes()[2]//2
			x,y,z = int(self.boxes[i][0]),int(self.boxes[i][1]),nz+int(self.boxes[i][2])

			bs = self.bsz_vs.value
			#if self.boxes[i][3] == 1:
			if self.boxes[i][3] != -1:
				print(self.boxes[i])
				r = self.data.get_clip(Region(x-bs//2,y-bs//2,z,bs,bs,1))
				l = self.get_annotation().get_clip(Region(x-bs//2,y-bs//2,z,bs,bs,1))
				r.write_image(outfile,-1)
				l.write_image(outfile,-1)
				#r=self.data.get_clip(Region(x-bs//2,y-bs//2,z-bz//2,bs,bs,bz))
		self.clear_shapes()
		return

	def annotate_from_curve(self,insert=None):
		try:
			if self.img_view.brush_mode == 0:
				self.annotate_from_curve_2D(insert=insert)
			else:
				self.annotate_from_curve_3D()
		except Exception as e:
			print(e)
			pass
		self.clear_shapes(reset_boxes_list=False)
		return

	def annotate_from_curve_2D(self, insert=None):
		self.save_current_state()
		#TODO: Clear points messing with annotate mechanism
		self.x = self.img_view.get_data().get_xsize()
		self.y = self.img_view.get_data().get_ysize()

		tsz=max(self.x,self.y)
		mask = np.zeros((self.x,self.y))
		#output=EMData(self.x, self.y, self.z)
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
			#mask = np.zeros((self.x,self.y,self.z))
			if len(self.curve.points) < 2:
				print("Need at least 2 points to draw a line")
				return
			#print(self.curve.points)
			r = np.array([int(item[0]) for item in self.curve.points])
			c = np.array([int(item[1]) for item in self.curve.points])
			z = np.array([int(item[2]) for item in self.curve.points])

			pt_dict = {}
			for i,pt in enumerate(self.curve.points):
				if pt[2] not in pt_dict.keys():
					pt_dict[pt[2]] = [i]
				else:
					pt_dict[pt[2]].append(i)

			print(pt_dict)
			if insert is None:
				ins = self.get_zpos()
			else:
				ins = insert

			for i in pt_dict[ins][:-1]:
				try:
					if self.curve.points[i][3] != -1:
						rr, cc, val = weighted_line(r[i], c[i], r[i+1],c[i+1], int(self.tx_line_width.text()))
						#mask[rr, cc] = int(self.tx_ann_class.text())
						#print("Len RR",len(rr))
						mask[rr, cc] = int(self.get_current_class())
						self.curve.points[i][3] = -1
					else:
						continue
				except:
					continue
				self.curve.points[pt_dict[ins][-1]][3] = -1
			#mask[r,c]=int(self.tx_ann_class.text())
			mask[r,c] = int(self.get_current_class())
			# 	mask[rr, cc] = int(self.tx_ann_class.text())
			# mask[r,c]=int(self.tx_ann_class.text())
		elif self.fill_type == "contour_fill":
			#mask = np.zeros((self.x,self.y))
			if len(self.contour.points) < 3:
				print("Need at least 3 points to draw a polygon")
				return
			rr, cc = polygon(self.contour.points)
			#print(self.contour.points)
			mask[rr,cc] = int(self.get_current_class())

		elif self.fill_type == "contour_nofill":
			#mask = np.zeros((self.x,self.y))
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
					#mask[rr, cc] = int(self.ct_ann_class.text())
					mask[rr,cc] = int(self.get_current_class())

				except:
					continue
			mask[r,c] = int(self.get_current_class())

		mask =  np.rot90(np.fliplr(mask))

		if insert is None:
			print("No insertion")
			ins_img = to_numpy(self.img_view.get_annotation())
			self.img_view.annotation = from_numpy(np.where(mask>0,mask,ins_img))
			self.img_view.force_display_update(set_clip=1)
			self.img_view.updateGL()
		else:
			center = self.zt_spinbox.value()
			if center == -1:
				center = self.nz//2
			# ins_img = to_numpy(self.img_view.get_full_annotation().get_clip(Region(0,0,center+int(insert),self.x,self.y,1)))
			# ann = from_numpy(np.where(mask>0,mask,ins_img))
			xform=Transform({"type":"eman","alt":self.img_view.alt,"az":self.img_view.az,"tx":self.img_view.nx//2,"ty":self.img_view.ny//2,"tz":self.img_view.nz//2+int(insert)})
			ins_img=to_numpy(self.img_view.get_full_annotation().get_rotated_clip(xform,(self.img_view.nx,self.img_view.ny,1),0))
			ann = from_numpy(np.where(mask>0,mask,ins_img))
			self.img_view.get_full_annotation().set_rotated_clip(xform,ann)
			self.img_view.force_display_update()
			self.img_view.updateGL()



	def annotate_from_curve_3D(self):
		self.save_current_state()
		def calc_middle_pts(start,end,interval):
			dist = np.linalg.norm(np.array(end)-np.array(start))
			if interval <1:
				interval=1
			no_pts = floor(dist//interval)
			if no_pts == 0:
				return []
			n_vec = (np.array(end)-np.array(start))/no_pts
			pts_l = []
			for i in range(no_pts):
				pts_l.append(start+n_vec*i)
			pts_l.append(end)
			return pts_l
		def is_point_in_polygon(point,vertices):
			points = np.array([[x[0],x[1]] for x in vertices])
			poly_path = mplPath.Path(points)
			return poly_path.contains_point(point)


		def find_bounding_square_box(vertices):
		# Initialize min and max coordinates for x and y
			min_x, min_y = float('inf'), float('inf')
			max_x, max_y = float('-inf'), float('-inf')

			# Iterate through each vertex to find min and max coordinates for x and y
			for x, y,_,_ in vertices:
				min_x, max_x = min(min_x, x), max(max_x, x)
				min_y, max_y = min(min_y, y), max(max_y, y)

			# Determine the length of the side of the square
			side_length = max(max_x - min_x, max_y - min_y) + 2
			if side_length % 2 == 1:
				side_length = side_length+1

			# Define the bottom-left and top-right corners of the square
			bottom_left = (min_x, min_y)
			top_right = (min_x + side_length, min_y + side_length)
			mid = (int(min_x) + side_length//2, int(min_y) + side_length//2)
			return mid,int(side_length)


		#
		# def find_bounding_box(vertices):
		# # Initialize min and max coordinates for x and y
		# 	min_x, min_y = float('inf'), float('inf')
		# 	max_x, max_y = float('-inf'), float('-inf')
		#
		# 	# Iterate through each vertex to find min and max coordinates for x and y
		# 	for x, y,_,_ in vertices:
		# 		min_x, max_x = min(min_x, x), max(max_x, x)
		# 		min_y, max_y = min(min_y, y), max(max_y, y)
		#
		# 	# Determine the length of the side of the square
		# 	side_length = max(max_x - min_x, max_y - min_y)
		# 	if side_length % 2 == 1:
		# 		side_length = side_length+1
		#
		# 	# Define the bottom-left and top-right corners of the square
		# 	bottom_left = (min_x, min_y)
		# 	top_right = (max_x, max_y)
		# 	mid = ((max_x+min_x)//2, (max_y+min_y)//2)
		# 	return mid, max_x-min_x,max_y-min_y

		def polygon(points):
			center, side_length = find_bounding_square_box(points)
			min_r = int(min([x[0] for x in points]))
			max_r = int(max([x[0] for x in points]))
			min_c = int(min([x[1] for x in points]))
			max_c = int(max([x[1] for x in points]))
			rr = []
			cc = []

			for r in range(min_r,max_r+1):
				for c in range(min_c,max_c+1):
					if is_point_in_polygon([r,c],points):
						rr.append(r-min_r)
						cc.append(c-min_c)
			return np.array(rr, dtype=np.intp), np.array(cc, dtype=np.intp)
			#mask =  np.rot90(np.fliplr(mask))

		def define_eroded_volume(middle_image, z = None):

			x, y = middle_image.shape[0],middle_image.shape[0]
			if z is None:
				z=middle_image.shape[0]
			vol = np.zeros((z,x,y))
			center_slice = z // 2
			for i in range(z):
				vol[i,:, :] = ndi.binary_erosion(middle_image,iterations=abs(z//2-i))
			vol[z//2,:, :] = middle_image
			return vol
		value = int(self.get_current_class())
		self.w_line=15
		if self.fill_type == "curve":
			if len(self.curve.points) < 2:
				print("Need at least 2 points to draw a line")
				return
			self.w_line = int(self.tx_line_width.text())//2
			pts_l = []
			for i in range(len(self.curve.points)-1):
				for p in calc_middle_pts(self.curve.points[i][0:3],self.curve.points[i+1][0:3],self.w_line/10):
					pts_l.append(p)

		elif self.fill_type == "contour_fill":
			if len(self.contour.points) < 3:
				print("Need at least 3 points for a contour")
				return

			#w_line = int(self.tx_line_width.text())
			# pts_l = sparse_polygon(self.contour.points)
			center, side_length = find_bounding_square_box(self.contour.points)
			print("Center,side_length", center, side_length)

			a = np.zeros((side_length,side_length))
			rr,cc = polygon(self.contour.points)
			a[rr,cc] = 1
			a = np.rot90(np.fliplr(a))
			fill_vol = define_eroded_volume(a,self.img_view.nz)
			# mask = np.zeros((self.img_view.nx,self.img_view.ny))
			# mask[rr,cc] = 1
			# fill_vol = define_eroded_volume(mask,self.img_view.nz)

		# 	mask[rr,cc] = int(self.get_current_class())
			#from_numpy(fill_vol).write_image("fillvol.hdf")

			#print("Fillvol shape",fill_vol.shape)


			xrot, yrot, zrot = self.img_view.reverse_transform_point([center[0],center[1],self.img_view.nz//2+self.img_view.zpos], self.img_view.get_xform())
			xform=Transform({"type":"eman","alt":self.img_view.alt,"az":self.img_view.az,"tx":self.img_view.nx//2+xrot,"ty":self.img_view.ny//2+yrot,"tz":self.img_view.nz//2+zrot})
			subvol=self.img_view.get_full_annotation().get_rotated_clip(xform,(side_length,side_length,self.img_view.nz),0)


			#subvol = from_numpy(np.ones((2*self.w_line,2*self.w_line,2*self.w_line)))
			#subvol.process_inplace("mask.paint",{"x":self.w_line,"y":self.w_line,"z":self.w_line,"r1":self.w_line,"v1":value,"r2":self.w_line,"v2":0})
			self.img_view.full_annotation.set_rotated_clip(xform,from_numpy(np.where(fill_vol>0,value,to_numpy(subvol))))
			#self.img_view.full_annotation.set_rotated_clip(xform,subvol*2)
			# full_vol=from_numpy(np.where(fill_vol>0,value,to_numpy(self.img_view.full_annotation)))
			# self.img_view.full_annotation = full_vol.process("xform",{"transform":self.img_view.xform})
			self.img_view.force_display_update()
			self.img_view.updateGL()
			return



		else:
			anchor_l = self.contour.points
			anchor_l.append(self.contour.points[0])
			self.w_line = int(self.ct_line_width.text())//2
			if self.w_line <5:
				interval = self.w_line
			elif self.w_line <20:
				interval = self.w_line/2
			else:
				interval = self.w_line/4
			pts_l = []

			for i in range(len(anchor_l)-1):

				for p in calc_middle_pts(anchor_l[i][0:3],anchor_l[i+1][0:3],interval):
					pts_l.append(p)


		for pts in pts_l:
			xrot, yrot, zrot = self.img_view.reverse_transform_point([pts[0],pts[1],self.img_view.nz//2+pts[2]], self.img_view.get_xform())
			xform=Transform({"type":"eman","alt":self.img_view.alt,"az":self.img_view.az,"tx":self.img_view.nx//2+xrot,"ty":self.img_view.ny//2+yrot,"tz":self.img_view.nz//2+zrot})
			subvol=self.img_view.get_full_annotation().get_rotated_clip(xform,(2*self.w_line,2*self.w_line,2*self.w_line),0)
			#subvol = from_numpy(np.ones((2*self.w_line,2*self.w_line,2*self.w_line)))
			fill_vol = subvol.process("mask.paint",{"x":self.w_line,"y":self.w_line,"z":self.w_line,"r1":self.w_line,"v1":value,"r2":self.w_line,"v2":0})
			#self.img_view.full_annotation.set_rotated_clip(xform,value*subvol.process("threshold.binaryrange",{"high":value+0.1,"low":value-0.1}))
			self.img_view.full_annotation.set_rotated_clip(xform,fill_vol)
		self.img_view.force_display_update()
		self.img_view.updateGL()




	#####ALL BASIC TAB LINEAR AND CONTOUR METHODS
	def basic_tab_change(self, tab_num):
		self.save_current_state()
		self.basic_tab_num = tab_num
		# print(self.basic_tab_num)
		if tab_num==3:

			try:
				self.annotate_from_curve(insert=self.get_zpos())
				self.curve.points = []
				self.contour.points =[]
				self.do_update()
			except:
				pass
			self.fill_type = None
			self.img_view.mouse_mode = 1
			return

		elif tab_num==0:
			#print("Zpos:",self.get_zpos())
			try:
				self.annotate_from_curve(insert=self.get_zpos())
				self.curve.points = []
				self.contour.points =[]
				self.do_update()
			except Exception as e:
				print("Can't paint due to", e)

				pass
			# if self.extract_bt:
			# self.extract_bt_clicked()
			self.img_view.mouse_mode = 5
			#print("Mouse mode is:",self.img_view.mouse_mode_dict[self.img_view.mouse_mode])
			#self.img_view.show_inspector(6)
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
			#print("Li",li)
			pt=pts[pts[:,3]==li][:,:3].copy()
			pt=pt[np.append(True, np.linalg.norm(np.diff(pt, axis=0), axis=1)>0.1)]
			ln=np.linalg.norm(np.diff(pt, axis=0), axis=1)
			ln=np.sum(ln)
			if len(pt)<2: continue
			#print len(pt), ln
			ipt=interp_points(pt, npt=np.round(ln)//density)
			pts_intp.append(np.hstack([ipt, li+np.zeros((len(ipt), 1)), self.curve.classid+np.zeros((len(ipt), 1))]))
			#print("Np_hstack",pts_intp)
			kk+=1

			# else:
			# 	continue

		pts_intp=np.vstack(pts_intp)
		self.curve.points=otherpts.tolist()
		self.curve.points.extend(pts_intp.tolist())
		self.do_update()

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

	def get_nz(self):
		hdr=EMData(self.data_file, 0,True)

		self.nz=hdr["nz"]
		#print("get_nz", self.nz)
		return self.nz

	def get_boxsize(self):
		return int(self.bsz_vs.value)
	def inside_box(self,i,x=-1,y=-1,z=0):
		#Checks to see if a point in image coordinates is inside box number n. If any value is negative, it will not be checked.
		box = self.boxes[i]
		bs=self.get_boxsize()/2

		rr=(x>=0)*((box[0]-x)**2) + (y>=0)*((box[1]-y) **2)
		if rr<=bs**2: print(i)
		return rr<=bs**2



	def key_press(self, event):
		if event.modifiers()&Qt.ControlModifier and event.key()==Qt.Key_Z :


			self.reverse_to_saved_state()
		elif event.key() == Qt.Key_Up or event.key() == Qt.Key_Down:
			if event.key() == Qt.Key_Up:
				self.previous_z = self.get_zpos()-1
			if event.key() == Qt.Key_Down:
				self.previous_z = self.get_zpos()+1
			if self.img_view.data:
				if self.basic_tab_num == 3:
					self.update_img_view_box()
				elif self.basic_tab_num == 2:
					if self.img_view.brush_mode == 0:
						self.annotate_from_curve(insert=self.previous_z)
						self.clear_shapes(reset_boxes_list=False)
					else:
						pass
				elif self.basic_tab_num == 1:
					try:
						if self.img_view.brush_mode == 0:
							self.annotate_from_curve(insert=self.previous_z)
						else:
							pass
					except:
						pass
		else:
			return

	def update_img_view_box(self):
		current_z = self.get_zpos()
		#print("Current z", current_z)
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


	def img_view_mouse_down(self, event):
		if self.basic_tab_num == 0:
			self.save_current_state()
		else:
			return


	def img_view_mouse_up(self, event):
		if self.basic_tab_num == 3:
			x,y=self.img_view.scr_to_img((event.x(),event.y()))
			z =self.img_view.zpos
			print(x,y,self.img_view.zpos)
			if not event.button()&Qt.LeftButton:
				return
			if event.modifiers()&Qt.ControlModifier:
				color = [0,0,255]
				if len(self.p1p2)>=1:
					if z != self.p1p2[0][2]:
						color = [128,128,255]
				self.p1p2.insert(0,[x,y,z])
				self.img_view.add_shape("p2",EMShape(("point",0,0,255,x,y,15)))
				try:
					self.img_view.add_shape("p1",EMShape(("point",color[0],color[1],color[2],self.p1p2[1][0],self.p1p2[1][1],15)))
				except:
					pass
				self.img_view.updateGL()

				if len(self.p1p2) > 2:
					self.p1p2.pop()
					#self.img_view.add_shape("p1",EMShape(("point",0,0,255,self.p1p2[0][0],self.p1p2[0][1],15)))
				if len(self.p1p2) == 2:
					print("Anchor point 1:{}, Anchor point 2:{}".format(self.p1p2[0],self.p1p2[1]))
				#color=QtGui.QColor.fromRgb(0,0,255)

					#self.img_view.add_shape("box{}".format(i),EMShape(("rect",color[0],color[1],color[2],x-sz//2,y-sz//2,x+(sz+1)//2,y+(sz+1)//2),2))



			elif event.modifiers()&Qt.ShiftModifier:
				for i in range(len(self.boxes)):
					if self.inside_box(i,x,y,0):
						self.boxes[i][3]=-1
						self.img_view.del_shape("box{}".format(i))
						self.img_view.updateGL()
			# 			#### remove point
			# 			self.curve.points=[p for i,p in enumerate(self.curve.points) if i!=ii]
			else:
				box_index = self.current_box_index_sp.value()
				self.boxes.append([x, y, z,1,box_index])
				self.add_boxes(size=int(self.bsz_vs.value))
				#print(self.boxes)

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
						pt_dict = {}
						for i,pt in enumerate(self.curve.points):
							if pt[2] not in pt_dict.keys():
								pt_dict[pt[2]] = [i]
							else:
								pt_dict[pt[2]].append(i)
						for ins in pt_dict.keys():
							self.annotate_from_curve(insert=ins)
						self.curve.points = []
						self.curve.add_point([x, y], True)

			else:
				#### add point
				self.curve.add_point([x, y]) #, self.img_view.list_idx
			self.do_update()
			#print(self.curve.points)

		#Curve tab
		elif self.basic_tab_num == 2:
			# self.save_current_state()
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
			#print(self.contour.points)
			self.img_view.updateGL()

		#Brush tab
		else:


			#self.save_current_state()
			# get_application().setOverrideCursor(Qt.ArrowCursor)
			# lc=self.img_view.scr_to_img(event.x(),event.y())
			# current_shapes = self.img_view.get_shapes()
			# if self.img_view.rmousedrag:
			# 	self.rmousedrag=None
			# if self.img_view.mouse_mode_dict[self.img_view.mouse_mode] == "emit":
			# 	#lc=self.scr_to_img(event.x(),event.y())
			# 	self.img_view.mouseup.emit(event, lc)
			# elif self.img_view.mouse_mode_dict[self.img_view.mouse_mode] == "measure":
			# 	if event.buttons()&Qt.LeftButton:
			# 		self.img_view.add_shape("MEAS",EMShape(("line",.5,.1,.5,current_shapes["MEAS"].shape[4],current_shapes["MEAS"].shape[5],lc[0],lc[1],2)))
			# # elif self.img_view.mouse_mode_dict[self.mouse_mode] == "draw":
			# # 	if event.button()==Qt.LeftButton:
			# # 		self.img_view.force_display_update()
			# # 		self.img_view.updateGL()
			# elif self.img_view.mouse_mode_dict[self.img_view.mouse_mode] == "seg":
			# 	if event.button()==Qt.LeftButton:
			# 		self.img_view.force_display_update(set_clip=1)
			#self.img_view.updateGL()
			#return
			# if event.button()==Qt.LeftButton:
				# self.img_view.force_display_update()
				# self.img_view.updateGL()
			#print("I don't know")
			return


	#TODO
	def do_morp_dilate(self):
		n_iters = int(self.morp_n_iters_sp.value())
		sels = self.get_target_selected()
		self.target.save_current_state()
		self.target.activateWindow()
		try:
			for sel in sels:
				val=int(sel.text(0))
				mask=self.target.get_annotation().process("threshold.binaryrange",{"high":val+0.1,"low":val-0.1})
				if mask:
					self.target.get_annotation().process_inplace("threshold.rangetozero",{"maxval":(val+0.1),"minval":(val-0.1)})
					self.target.annotate += val*from_numpy(ndi.binary_dilation(to_numpy(mask),iterations=n_iters))

					del mask
			self.target.img_view.set_data(self.target.data, self.target.annotate)
		except Exception as e:
			print(e)


	def random_bx_bt_clicked(self):

		try:
			no_box = int(self.random_bx_sb.getValue())
		except:
			print("Number of boxes needed to be integer")
			return
		pos_annotate = self.get_segtab().get_selected_annotate().process("threshold.binary",{"value":0.1})
		#Dilating selected so that the background is not overlapped.
		pos_np_dilated = ndi.binary_dilation(to_numpy(pos_annotate),iterations=self.bsz_vs.value//2+25)
		#self.img_view.set_data(self.get_annotation(),from_numpy(pos_np_dilated))
		possible_boxes = list(zip(*np.where(pos_np_dilated==0)))

		#ind = np.random.randint(0, len(possible_boxes), no_box)
		if self.random_sparse_checkbox.isChecked():
			ind = np.random.randint(0, len(possible_boxes)//100, no_box)*100
		else:
			ind = np.random.randint(0, len(possible_boxes), no_box)

		# xs = np.random.randint(self.bsz_vs.value//2,self.img_view.get_data_dims()[0]-self.bsz_vs.value//2,no_box)
		# ys = np.random.randint(self.bsz_vs.value//2,self.img_view.get_data_dims()[1]-self.bsz_vs.value//2,no_box)
		box_index = self.current_box_index_sp.value()
		for i in ind:
			self.boxes.append([possible_boxes[i][1],possible_boxes[i][0],self.get_zpos(),1,box_index])
			#self.boxes.append([xs[i],ys[i],self.get_zpos(),1,box_index])
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
			nz = self.data.get_sizes()[2]//2
			x,y,z = int(self.boxes[i][0]),int(self.boxes[i][1]),nz+int(self.boxes[i][2])
			print("X,Y,Z", x,y,z)
			bs = self.bsz_vs.value
			if self.boxes[i][3] != -1:
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

	# def on_table_tom(self,row,col):
	# 	print("Row", row, "Col", col, "Tomogram", self.table_tom.itemAt(row,col))
	# 	print(self.table_tom.currentItem().text())

	# def table_tom_cell_changed(self):
	# 	#print("New Cell", self.table_tom.currentItem().text())
	# 	self.img_view.set_data(EMData(self.table_tom.currentItem().text()),None)
	def set_data(self):
		return

	# def update_sets(self):
	# 	#Set the colors and flags of table set items
	# 	for i in range(self.table_tom.rowCount()):
	# 		key = int(self.table_tom.item(i,0).text())
	# 		self.table_tom.item(i,0).setFlags(self.indexflags)
	# 		self.table_tom.item(i,1).setFlags(self.itemflags)
	# 		self.table_tom.item(i,0).setForeground(self.colors[key])
	# 		self.table_tom.item(i,1).setForeground(self.colors[key])



	def add_boxes(self, size = 64):
		self.boxes = [box for box in self.boxes if box[3] != -1]
		sz = size
		#color = [0.1,0.1,0.3]
		#print("add boxes", self.boxes)
		for i in range(len(self.boxes)):
			#print(self.boxes[i])
			color=QtGui.QColor(self.setcolors[self.boxes[i][4]%len(self.setcolors)]).getRgbF()
			x,y = int(self.boxes[i][0]),int(self.boxes[i][1])
			if self.boxes[i][3] == 1:
				#self.img_view.add_shape("box{}".format(i),EMShape(("rect",color[0],color[1],color[2],x-sz//2,y-sz//2,x+(sz+1)//2,y+(sz+1)//2),2))
				self.img_view.add_shape("box{}".format(i),EMShape(("rect",color[0],color[1],color[2],x-old_div(sz,2),y-old_div(sz,2),x+old_div((sz+1),2),y+old_div((sz+1),2),2)))
		self.img_view.updateGL()

	def SaveJson(self):
		info=js_open_dict("./info/annotate_project.json")
		info["tom_list"]=self.tom_file_list
		info.close()

	#Function from system menu
	def menu_file_del_all(self):
		self.img_view.full_annotation = EMData(self.nx,self.ny,self.nz)
		self.img_view.force_display_update()
		self.img_view.updateGL()

	def menu_file_new_seg(self):

		hdr=EMData(self.data_file, 0,True)
		self.nx=hdr["nx"]
		self.ny=hdr["ny"]
		self.nz=hdr["nz"]
		self.apix = hdr["apix_x"]
		current_child_count = self.tomo_tree.currentItem().childCount()
		# 	seg_path = os.path.join(self.seg_folder,self.seg_grps[self.tomo_tree.indexOfTopLevelItem(self.tomo_tree.currentItem())].checkedButton().text())
		# 	print("Seg file for current tomogram already exists at",seg_path)
		# else:
		seg_path = os.path.join(self.seg_folder,"{}_seg_{}.hdf".format(base_name(self.data_file),str(current_child_count+1)))
		print("Create seg_file",seg_path)
		seg_out = EMData(self.nx,self.ny,self.nz)
		seg_out.write_image(seg_path)
		del seg_out
		seg_item = QtWidgets.QTreeWidgetItem()
		self.tomo_tree.currentItem().addChild(seg_item)
		seg_item.setFlags(Qt.ItemFlags(Qt.ItemIsEnabled))
		#seg_item.setCheckState(0,0)
		seg_button = QtWidgets.QRadioButton(os.path.basename(seg_path))
		self.seg_grps[self.tomo_tree.indexOfTopLevelItem(self.tomo_tree.currentItem())].addButton(seg_button,0)
		seg_button.setChecked(True)
		self.tomo_tree.setItemWidget(seg_item,0,seg_button)
		seg_info_path = os.path.join(seg_path[0:-4]+".json")

		self.seg_path = seg_path
		self.seg_info_path = seg_info_path

		self.data_xy = self.thumbnail.get_xy()
		self.reset_morp_params(reset_vs=False)
		self.set_imgview_data(round(self.data_xy[0]),round(self.data_xy[1]),self.img_view_region_size)
		self.update_tree_set()

	def menu_file_quit(self):
		self.close()
	def menu_edit_undo(self):
		self.reverse_to_saved_state()



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
		try:
			E2saveappwin("e2annotate","main",self)
			E2saveappwin("e2annotate","controlpanel",self.control_panel)
			E2saveappwin("e2annotate","tomograms",self.tomo_list_panel)
		except:
			pass
		#self.write_header(self.get_annotation())

		#self.write_metadata(self.seg_path)
		print(self.get_annotation(), self.seg_path, self.cur_region)
		try:
			self.get_inspector().seg_tab.write_treeset_json(self.seg_info_path)
		except:
			print("Cannot print class metadata to info file.")
			pass
		#self.get_treeset().clear()

		try:
			print(self.get_annotation().numpy().shape, self.cur_region)
			self.write_out(self.get_annotation(), self.seg_path, self.cur_region)
		except:
			print("Cannot write annotation to segs file.")
			pass

		self.control_panel.close()
		self.tomo_list_panel.close()

		for widget in self.popwidgets:
			print("Closing",widget)
			widget.close()

		#self.get_annotation().write_image(self.seg_path, 0, IMAGE_HDF, False, self.cur_region)
		self.close()






class Thumbnail(EMImage2DWidget):
	def __init__(self,current_file=None,target=None,app_target = None,tn_size=220):
		super().__init__()
		if current_file:
			self.current_file = current_file
			self.get_im(self.current_file)

		self.scale_fac = round(self.img["nx"]/tn_size)
		print(self.scale_fac)
		try:
			self.thumbnail_image=self.img.process("math.meanshrink",{"n":(self.scale_fac)})
		except:
			temp=to_numpy(self.img)
			self.scale_fac = temp.shape[0]
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
		self.iv_size = [self.app_target.img_view_region_size, self.app_target.img_view_region_size]
		box_size = [self.size*(self.iv_size[0]/self.im_scalesz),self.size*(self.iv_size[1]/self.im_scalesz)]
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
					self.app_target.write_out(self.app_target.get_annotation(), self.app_target.seg_path, self.app_target.cur_region)
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
		self.add_shape("WinView",EMShape(("rectpoint",.9,0,0,x-old_div(sz[0],2),y-old_div(sz[1],2),x+old_div((sz[0]+1),2),y+old_div((sz[1]+1),2),2)))
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
		czpos = self.image.zpos
		if len(self.points)==0 and len(newpt)>0:
			#self.points.append([newpt[0], newpt[1], zpos, self.select, self.classid])
			self.points.append([newpt[0], newpt[1], czpos, self.select, self.classid])

		if newcontour==False:
			ci=self.select
			#### separate the points on the current contour and everything else
			nppts=np.array(self.points)
			sel=np.logical_and(nppts[:,3]==ci, nppts[:,4]==self.classid)
			pts=nppts[sel,:3].copy()
			otherpts=nppts[sel==False].copy()

			if len(pts)<3 or optimize==False:
				if len(newpt)>0:
					self.points.append([newpt[0], newpt[1], czpos, ci, self.classid])
					return
			else:
				thr=1000.
				newpt=np.array([newpt[0], newpt[1], czpos])


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

		else:
			#### start a new contour
			pts=np.array(self.points)
			ci=np.max(pts[:,3])+1
			self.select=ci
			#self.current_zpos = self.image.zpos
			self.points.append([newpt[0], newpt[1], czpos, ci, self.classid])


	def draw(self,d2s=None,col=None):
		zpos=self.image.zpos
		curpts=[p for p in self.points if p[4]==self.classid]
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
		return

class Contour(EMShape):
	def __init__(self, img=None, points=[]):
		self.points=points
		self.isanimated = False
		self.shape=["scr_contour",0]
		self.image=img
		self.triangles=[]


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

		#self.model=self.get_tiny_unet()
		#self.model=self.get_unet()
		self.model = self.get_unet_23()
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

	def get_unet_23(inp_x=None,inp_y=None):
		inputs = Input((None, None, 1))
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
		pool4 = MaxPooling2D(pool_size=(2, 2))(conv4)

		conv5 = Conv2D(512, (3, 3), activation='relu', padding='same')(pool4)
		conv5 = Conv2D(512, (3, 3), activation='relu', padding='same')(conv5)

		up6 = concatenate([Conv2DTranspose(256, (2, 2), strides=(2, 2), padding='same')(conv5), conv4], axis=3)
		conv6 = Conv2D(256, (3, 3), activation='relu', padding='same')(up6)
		conv6 = Conv2D(256, (3, 3), activation='relu', padding='same')(conv6)

		up7 = concatenate([Conv2DTranspose(128, (2, 2), strides=(2, 2), padding='same')(conv6), conv3], axis=3)
		#up7 = concatenate([Conv2DTranspose(128, (2, 2), strides=(2, 2), padding='same')(conv4), conv3], axis=3)
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

	def train_unet(self,weights_out='./neural_nets/weights_temp.h5',no_epoch = 30, batch_sz= 128,val_split=0.2,learnrate=3e-4):
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
		self.model.compile(optimizer=Adam(learning_rate=learnrate,clipnorm=1), loss=dice_coef_loss, metrics=[dice_coef])
		model_checkpoint = ModelCheckpoint(self.weights_out, monitor='val_loss', save_best_only=True)

		print('-'*30)
		print('Fitting model...')
		print('-'*30)
		self.model.fit(self.datas, self.labels, batch_size=batch_sz, epochs=no_epoch,verbose=1, shuffle=True,validation_split=val_split,callbacks=[model_checkpoint])
		#return model,history
		print("Done training model. Network saved to",weights_out)

	def load_model(self,weights_in, tiny=False):
		if tiny:
			model = self.get_tiny_unet()
		else:
			#model = self.get_unet()
			#model = self.get_tiny_unet()
			model = self.get_unet_23()
		model.load_weights(weights_in)
		return model


	def apply_unet(self,weights_in,tomogram=None,outfile=None):
		if tomogram==None:
			print("Need to specify tomogram to apply U-net")
			return
		if not os.path.exists("./neural_nets/weights_temp.h5"):
			print("There's no neural networks trained for this tomogram. Train a neural net first")
			return
		else:
			pass
		is3d=True
		enx,eny=tomogram["nx"], tomogram["ny"]
		nframe=tomogram["nz"]
		tsz=max(enx,eny)
		output=EMData(enx, eny, nframe)


		print("Loading tomogram...")
		tomo_in=[]
		for nf in range(nframe):
			if enx %2 ==1:

				e0=tomogram.get_clip(Region((enx-1-tsz)//2,(eny-1-tsz)//2,nf,tsz,tsz,1))
			else:
				e0=tomogram.get_clip(Region((enx-tsz)//2,(eny-tsz)//2,nf,tsz,tsz,1))
			tomo_in.append(e0)
		print(len(tomo_in))
		for idx, img in enumerate(tomo_in):


			# if idx == 123:
			# 	plt.imshow(tomo_in[idx].numpy())
			im = img.get_clip(Region(0,0,960,960))
			#m=img.numpy()
			m=im.numpy()


			#m = self.normalize(m)

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
		return


class Simple_NNet_Tab(QtWidgets.QWidget):
	def __init__(self,target):
		QtWidgets.QWidget.__init__(self,None)
		#self.target=weakref.ref(target)
		self.target = target

		self.snnet_tab = QtWidgets.QTabWidget()
		self.snn_data = QtWidgets.QWidget()
		self.snn_train = QtWidgets.QWidget()
		self.snn_apply = QtWidgets.QWidget()
		self.snnet_tab.addTab(self.snn_data,"Create Train Data")
		self.snnet_tab.addTab(self.snn_train,"Train NNet")
		self.snnet_tab.addTab(self.snn_apply,"Apply NNet")

		self.map_file_lst = []
		self.no_maps = 5
		self.box_training_bt = QtWidgets.QPushButton("Create Training Data")

		self.train_nnet_sb_list = []

		self.trainset_sb = StringBox(label="trainset      ")
		self.trainset_browser_bt = QtWidgets.QPushButton("Browse")



		self.start_nnet_sb = StringBox(label="from_trained")
		self.start_nnet_browser_bt = QtWidgets.QPushButton("Browse")


		self.train_nnet_sb_list.append(StringBox(label="learnrate",value="0.0003",showenable=-1))
		self.train_nnet_sb_list.append(StringBox(label="niter  ",value="20",showenable=-1))
		self.train_nnet_sb_list.append(StringBox(label="ncopy    ",value="3",showenable=-1))
		self.train_nnet_sb_list.append(StringBox(label="batch ",value="20",showenable=-1))
		self.train_nnet_sb_list.append(StringBox(label="nkernel  ",value="40,40,1",showenable=-1))
		self.train_nnet_sb_list.append(StringBox(label="ksize  ",value="15,15,15",showenable=-1))
		self.train_nnet_sb_list.append(StringBox(label="poolsz   ",value="2,1,1",showenable=-1))
		self.train_nnet_sb_list.append(StringBox(label="device",value="cpu",showenable=-1))
		self.train_nnet_bt = QtWidgets.QPushButton("Train NNet")

		self.train_nnet_launcher = QtWidgets.QWidget()
		self.train_nnet_launcher.setWindowTitle("Command for train nnet")
		self.train_nnet_launcher.setMinimumWidth(500)
		self.train_nnet_cmd = QtWidgets.QTextEdit()
		self.train_nnet_launch_bt = QtWidgets.QPushButton("Launch")

		tnl_layout = QtWidgets.QVBoxLayout()
		tnl_layout.addWidget(self.train_nnet_cmd)
		tnl_layout.addWidget(self.train_nnet_launch_bt)
		self.train_nnet_launcher.setLayout(tnl_layout)



		self.trained_nnet_sb = StringBox(label="nnet")
		self.trained_nnet_browser_bt = QtWidgets.QPushButton("Browse")
		self.apply_nnet_bt = QtWidgets.QPushButton("Apply NNet")
		self.apply_eman2_thread_sb = StringBox(label="threads")
		self.apply_eman2_thread_sb.setValue("12")
		self.apply_outtag_sb = StringBox(label="outtag")
		self.apply_outtag_sb.setValue("eman2_seg")

		self.apply_nnet_launcher = QtWidgets.QWidget()
		self.apply_nnet_launcher.setWindowTitle("Command for apply nnet")
		self.apply_nnet_launcher.setMinimumWidth(500)
		self.apply_nnet_cmd = QtWidgets.QTextEdit()
		self.apply_nnet_launch_bt = QtWidgets.QPushButton("Launch")

		anl_layout = QtWidgets.QVBoxLayout()
		anl_layout.addWidget(self.apply_nnet_cmd)
		anl_layout.addWidget(self.apply_nnet_launch_bt)
		self.apply_nnet_launcher.setLayout(anl_layout)


		self.prob_maps_sb = StringBox(label="Prob Maps")
		self.prob_maps_browser_bt = QtWidgets.QPushButton("Browse")
		self.merger_tool_bt = QtWidgets.QPushButton("Merger Tool")


		gbl = QtWidgets.QGridLayout(self)
		# gbl.addWidget(self.segment_box_bt,0,0,1,4)
		# gbl.addWidget(self.prob_maps_sb,1,0,1,3)
		# gbl.addWidget(self.prob_maps_browser_bt,1,3,1,1)
		# gbl.addWidget(self.merger_tool_bt,2,0,1,4)
		gbl.addWidget(self.snnet_tab,0,0,1,4)


		data_gbl = QtWidgets.QGridLayout(self.snn_data)
		data_gbl.addWidget(self.box_training_bt,0,0,1,3)
		train_gbl = QtWidgets.QGridLayout(self.snn_train)
		train_gbl.addWidget(self.trainset_sb,0,0,1,3)
		train_gbl.addWidget(self.trainset_browser_bt,0,3,1,1)
		train_gbl.addWidget(self.start_nnet_sb,1,0,1,3)
		train_gbl.addWidget(self.start_nnet_browser_bt,1,3,1,1)
		for i in range(0,len(self.train_nnet_sb_list),2):
			train_gbl.addWidget(self.train_nnet_sb_list[i],2+i//2,0,1,2)
			train_gbl.addWidget(self.train_nnet_sb_list[i+1],2+i//2,2,1,2)


		train_gbl.addWidget(self.train_nnet_bt,6,0,1,4)


		apply_gbl = QtWidgets.QGridLayout(self.snn_apply)
		apply_gbl.addWidget(self.trained_nnet_sb,0,0,1,3)
		apply_gbl.addWidget(self.trained_nnet_browser_bt,0,3,1,1)
		apply_gbl.addWidget(self.apply_eman2_thread_sb,1,0,1,1)
		apply_gbl.addWidget(self.apply_nnet_bt,1,3,1,1)
		apply_gbl.addWidget(self.apply_outtag_sb,1,1,1,2)


		apply_gbl.addWidget(self.prob_maps_sb,2,0,1,3)
		apply_gbl.addWidget(self.prob_maps_browser_bt,2,3,1,1)
		apply_gbl.addWidget(self.merger_tool_bt,3,0,1,4)




		self.trainset_browser_bt.clicked[bool].connect(self.load_trainset)
		self.start_nnet_browser_bt.clicked[bool].connect(self.load_start_nnet)

		self.box_training_bt.clicked[bool].connect(self.box_training_bt_clicked)
		self.prob_maps_browser_bt.clicked[bool].connect(self.load_masks)
		self.merger_tool_bt.clicked[bool].connect(self.merger_tool_bt_clicked)
		self.train_nnet_bt.clicked[bool].connect(self.train_nnet_bt_clicked)
		self.train_nnet_launch_bt.clicked[bool].connect(self.train_nnet_launched)
		self.trained_nnet_browser_bt.clicked[bool].connect(self.load_trained_nnet)
		self.apply_nnet_bt.clicked[bool].connect(self.apply_nnet_bt_clicked)
		self.apply_nnet_launch_bt.clicked[bool].connect(self.apply_nnet_launched)


	def load_trainset(self):
		self.trainset_browser = EMBrowserWidget(withmodal=True,multiselect=False,startpath='./particles',dirregex='*_trainset.hdf')
		self.trainset_browser.ok.connect(self.load_trainset_browser_ok)
		self.trainset_browser.show()
		return

	def load_trainset_browser_ok(self):
		trainset = self.trainset_browser.getResult()
		self.trainset_sb.setValue(trainset[0])

	def load_start_nnet(self):
		self.start_nnet_browser = EMBrowserWidget(withmodal=True,multiselect=False,startpath='./neuralnets',dirregex='nnet_save*')
		self.start_nnet_browser.ok.connect(self.load_start_nnet_browser_ok)
		self.start_nnet_browser.show()
		return

	def load_start_nnet_browser_ok(self):
		start_nnet = self.start_nnet_browser.getResult()
		self.start_nnet_sb.setValue(start_nnet[0])
		return

	def load_trained_nnet(self):
		self.trained_nnet_browser = EMBrowserWidget(withmodal=True,multiselect=False,startpath='./neuralnets',dirregex='nnet_save*')
		self.trained_nnet_browser.ok.connect(self.load_trained_nnet_browser_ok)
		self.trained_nnet_browser.show()
		return

	def load_trained_nnet_browser_ok(self):
		trained_nnet = self.trained_nnet_browser.getResult()
		self.trained_nnet_sb.setValue(trained_nnet[0])
		return


	def visualize_seg_browser_ok(self):
		self.vis_seg_browser_ret = self.visualize_seg_browser.getResult()
		self.seg_file_le.setText(self.vis_seg_browser_ret[0])


	def load_masks(self):
		#Load probability maps to merge in EMAnnotateWidget
		self.openbrowser = EMBrowserWidget(withmodal=True,multiselect=True,startpath='./segmentations')
		self.openbrowser.ok.connect(self.load_masks_browser_ok)
		self.openbrowser.show()
		return

	def load_masks_browser_ok(self):
		self.browser_ret = (self.openbrowser.getResult())
		self.prob_maps_sb.setValue(",".join(self.browser_ret))



	def merger_tool_bt_clicked(self):
		map_file_lst = self.prob_maps_sb.getValue().split(",")
		print(map_file_lst)
		if len(map_file_lst) == 0:
			print("Specify probability maps for merging")
			return
		else:
			self.merger_tool = Merger_Tool(target=self.target, maps_lst=map_file_lst)
			self.target.popwidgets.append(self.merger_tool)
			return

	def box_training_bt_clicked(self):
		self.boxer_widget = Boxer_Widget(target=self.target)
		self.target.popwidgets.append(self.boxer_widget)

	def box_training_bt_clicked(self):
		self.boxer_widget = Boxer_Widget(target=self.target)
		self.target.popwidgets.append(self.boxer_widget)

	def train_nnet_bt_clicked(self):

		cmd_l = [self.trainset_sb.getValue()]
		for sb in self.train_nnet_sb_list:
			cmd_l.append(sb.getValue())
		self.tnn_cmd = "e2tomoseg_convnet.py --trainset={} --learnrate={} --niter={} --ncopy={} --batch={} --nkernel={} --ksize={} --poolsz={} --trainout --training --device={}".format(*cmd_l)
		if len(self.start_nnet_sb.getValue()) > 0:
			self.tnn_cmd += " --from_trained={}".format(self.start_nnet_sb.getValue())

		self.train_nnet_cmd.setText(self.tnn_cmd)
		self.train_nnet_launcher.show()
		return

	def train_nnet_launched(self):
		self.train_nnet_launcher.close()
		try:
			os.system(self.train_nnet_cmd.toPlainText())
		except Exception as e:
			print(e)
			pass
		#tplt_match_map = EMData("tmp_ccc_inv.hdf")
		return

	def apply_nnet_bt_clicked(self):

		self.ann_cmd = "e2tomoseg_convnet.py --nnet={} --tomograms={} --applying --outtag={} --threads={} --device=cpu".format(self.trained_nnet_sb.getValue(),self.target.data_file,self.apply_outtag_sb.getValue(),self.apply_eman2_thread_sb.getValue())
		#print(self.ann_cmd)
		self.apply_nnet_cmd.setText(self.ann_cmd)
		self.apply_nnet_launcher.show()
		return

	def apply_nnet_launched(self):
		self.apply_nnet_launcher.close()
		print("Apply nnet to current tomogram")
		try:
			os.system(self.apply_nnet_cmd.toPlainText())
		except Exception as e:
			print(e)
			pass
		return



class Boxer_Widget(QtWidgets.QWidget):
	def __init__(self,target):
		QtWidgets.QWidget.__init__(self,None)
		self.target = target
		self.jsonfile = info_name(self.target.data_file)
		print(self.jsonfile)
		if os.path.isfile(self.jsonfile):
			info = js_open_dict(self.jsonfile)
		else:
			print("No info file was detected for the current datafile. Please provide the correct info name in the box below")
			return

		self.sets={}
		self.boxsize={}
		self.all_boxes=[]
		if "class_list" in info:
			clslst=info["class_list"]
			for k in sorted(clslst.keys()):
				if type(clslst[k])==dict:
					self.sets[int(k)]=str(clslst[k]["name"])
					self.boxsize[int(k)]=int(clslst[k]["boxsize"])
				else:
					self.sets[int(k)]=str(clslst[k])
					self.boxsize[int(k)]=64
		#print(self.sets)
		self.set_table = QtWidgets.QTableWidget(len(self.sets),4)
		self.set_table.setHorizontalHeaderLabels(["Index","Class Name","Box Size","Stored Boxes"])
		for i in range(len(self.sets)):
			key = list(self.sets.keys())[i]
			self.set_table.setItem(i,0,QtWidgets.QTableWidgetItem(str(key)))
			self.set_table.setItem(i,1,QtWidgets.QTableWidgetItem(str(self.sets[key])))
			self.set_table.setItem(i,2,QtWidgets.QTableWidgetItem(str(self.boxsize[key])))
			#self.set_table.setItem(i,3,QtWidgets.QTableWidgetItem(str(self.boxsize[key])))

		self.extract_train_bt = QtWidgets.QPushButton("Extract Train Data")
		self.extract_train_le = QtWidgets.QLineEdit()
		self.extract_train_ncopy_sb = QtWidgets.QSpinBox()
		self.extract_train_ncopy_sb.setValue(1)
		self.extract_train_ncopy_sb.setMinimum(1)

		self.show_train_bt = QtWidgets.QPushButton("Show Train Data")



		self.visualize_seg_bt = QtWidgets.QPushButton("Visualize Seg")
		self.visualize_seg_browser = EMBrowserWidget(withmodal=True,multiselect=False,startpath='./particles')
		self.seg_file_le = QtWidgets.QLineEdit()
		self.browse_seg_file_bt = QtWidgets.QPushButton("Browse")
		self.visualize_seg_browser.ok.connect(self.visualize_seg_browser_ok)

		#self.extract_train_le = QtWidgets.QLineEdit()
		#self.extract_train_le.setText(os.path.join('./trainset/',base_name(self.target.data_file)[0:-4]+"_trainset.hdf"))

		if "boxes_3d" in info:
			box=info["boxes_3d"]

			stored_l = [i[5] for i in box]
			for i in range(len(self.sets)):
				key = list(self.sets.keys())[i]
				self.set_table.setItem(i,3,QtWidgets.QTableWidgetItem(str(stored_l.count(key))))
			for i,b in enumerate(box):

				#### X-center,Y-center,Z-center,method,[score,[class #]]
				bdf=[0,0,0,"manual",0.0, 0, 0]
				for j,bi in enumerate(b):  bdf[j]=bi
				if bdf[5] not in list(self.sets.keys()):
					clsi=int(bdf[5])
					self.sets[clsi]="particles_{:02d}".format(clsi)
					self.boxsize[clsi]=64
				self.all_boxes.append(bdf)

		et_hlay = QtWidgets.QHBoxLayout()
		et_hlay.addWidget(QtWidgets.QLabel("Ncopy:"))
		et_hlay.addWidget(self.extract_train_ncopy_sb)

		bw_gbl=QtWidgets.QGridLayout(self)
		bw_gbl.addWidget(self.set_table,0,0,2,3)
		bw_gbl.addWidget(self.visualize_seg_bt,2,0,1,1)
		bw_gbl.addWidget(self.seg_file_le,2,1,1,1)
		bw_gbl.addWidget(self.browse_seg_file_bt,2,2,1,1)
		bw_gbl.addWidget(self.extract_train_bt,3,0,1,1)
		bw_gbl.addWidget(self.extract_train_le,3,1,1,1)
		bw_gbl.addWidget(self.show_train_bt,4,0,1,1)

		bw_gbl.addLayout(et_hlay,3,2,1,1)


		self.set_table.currentItemChanged[QtWidgets.QTableWidgetItem,QtWidgets.QTableWidgetItem].connect(self.set_table_item_changed)
		self.set_table.itemClicked[QtWidgets.QTableWidgetItem].connect(self.on_item_clicked)
		self.extract_train_bt.clicked[bool].connect(self.extract_boxes)
		self.show_train_bt.clicked[bool].connect(self.show_train_data)
		self.browse_seg_file_bt.clicked[bool].connect(self.browse_seg_bt_clicked)
		self.visualize_seg_bt.clicked[bool].connect(self.visualize_seg_bt_clicked)
		self.setWindowTitle("Boxer Widget")
		self.show()


	def on_item_clicked(self, item):
		set_index = [int(self.set_table.item(self.set_table.row(item),0).text()) for item in self.set_table.selectedItems()]
		print("set_index", set_index)
		#set_index = int(self.set_table.item(self.set_table.row(item),0).text())
		sz = int(self.set_table.item(self.set_table.row(item),2).text())
		self.target.img_view.del_shapes()
		nz = self.target.get_nz()//2

		# for other_box in self.target.boxes:
		# 	other_box[3] = 0
		self.target.boxes = []
		for box in self.all_boxes:
			if box[5] in set_index:
				if box[2]-nz == 0:
					self.target.boxes.append([box[0],box[1],box[2]-nz,1,box[5]])
				else:
					self.target.boxes.append([box[0],box[1],box[2]-nz,0,box[5]])


		zs = [abs(box[2]) for box in self.target.boxes]
		print("Change z-thickness to", max(zs),"to visualize all boxes from file")
		self.target.zt_spinbox.setValue(max(zs))
		self.target.current_box_index_sp.setValue(set_index[0])
		self.target.basic_tab.setCurrentIndex(3)
		bs = self.target.bsz_vs.setValue(sz)
		self.target.add_boxes(sz)

		tentative_seg = self.get_seg_file_name()
		self.extract_train_le.setText(tentative_seg[0:-7]+"_trainset.hdf")
		if os.path.isfile(tentative_seg):
			self.seg_file_le.setText(tentative_seg)
		else:
			self.seg_file_le.setText("No seg file detected for current selections")
		return


	def get_seg_file_name(self):
		#print("Current Item", self.set_table.selectedItems())
		name_list = [self.set_table.item(self.set_table.row(item),1).text() for item in self.set_table.selectedItems()]
		p0=self.target.data_file.find('__')
		if p0>0:
			p1=self.target.data_file.rfind('.')
			filetag=self.target.data_file[p0:p1]
			if filetag[-1]!='_':
				filetag+='_'
		else:
			filetag="__"
		#name= self.set_table.item(self.set_table.row(item),1).text()+"_seg.hdf"
		name = "_".join(name_list)+"_seg.hdf"
		name=filetag+str(name)
		dr = "./particles/"
		return os.path.join(dr,base_name(self.target.data_file))+name



	def get_select_indexes(self):
		return [int(self.set_table.item(self.set_table.row(item),0).text()) for item in self.set_table.selectedItems()]

	def extract_boxes(self):
		try:
			print("Saving annotation at", self.target.seg_path)
			self.target.write_out(self.target.get_annotation(), self.target.seg_path, self.target.cur_region)
		except Exception as e:
			print("Cannot extract box due to error", e, "while writing annotation to segs file. Abort")
			return

		# if not os.path.exists("./trainset"):
		# 	os.mkdir("./trainset")

		print("Extract image and label patches of size",self.target.bsz_vs.value) #"at positions:",self.target.boxes)
		outfile_trainset = self.extract_train_le.text()
		tentative_seg = self.get_seg_file_name()
		outfile_raw = tentative_seg[0:-7]+"_raw.hdf"
		outfile_seg = tentative_seg[0:-7]+"_lab.hdf"

		try:
			os.remove(outfile_raw)
			os.remove(outfile_seg)
		except:
			pass

		print("Extract Boxes")
		for i in range(len(self.target.boxes)):
			nz = self.target.get_nz()//2
			x,y,z = int(self.target.boxes[i][0]),int(self.target.boxes[i][1]),int(self.target.boxes[i][2]+nz)
			bs = self.target.bsz_vs.value
			target_data =  self.target.get_full_data_from_file()
			target_annotation = self.target.get_full_annotation_from_file()
			#if self.boxes[i][3] == 1:
			# if self.target.boxes[i][3] != 1:
			# 	pass
			# else:
			box_region = Region(x-bs//2,y-bs//2,z,bs,bs,1)
			r = target_data.get_clip(Region(x-bs//2,y-bs//2,z,bs,bs,1))
			l = target_annotation.get_clip(Region(x-bs//2,y-bs//2,z,bs,bs,1))

			try:
				r.write_image(outfile_raw,-1)
				l.write_image(outfile_seg,-1)

			except Exception as e:
				print("Trainset file is not correctly formatted. Abort")
				print(e)
				return
			sys.stdout.write("\r{}/{} boxes written".format(i+1,len(self.target.boxes)))
			sys.stdout.flush()

		print("\nBuilding trainset")
		os.system("e2tomoseg_buildtrainset.py --buildset --particles_raw={} --particles_label={} --ncopy={} --trainset_output={} --validset=0.0".format(outfile_raw,outfile_seg,str(self.extract_train_ncopy_sb.value()),outfile_trainset))
		try:
			del target_data, target_annotation
		except:
			pass
		self.target.clear_shapes()
		return

	def show_train_data(self):

		trainset_path = self.extract_train_le.text()
		if not os.path.exists(trainset_path):
			print("No trainset file detected in the ./particles folder. Create one first")
			return
		else:
			os.system("e2display.py "+trainset_path)
			return


	def browse_seg_bt_clicked(self):
		self.visualize_seg_browser.show()
		return

	#self.visualize_seg_browser = EMBrowserWidget(withmodal=True,multiselect=False,startpath='./particles')

	def visualize_seg_browser_ok(self):
		self.vis_seg_browser_ret = self.visualize_seg_browser.getResult()
		self.seg_file_le.setText(self.vis_seg_browser_ret[0])

	def visualize_seg_bt_clicked(self):
		#print(self.all_boxes)
		seg_l = EMData.read_images(self.seg_file_le.text())
		#boxes = [box for box in self.all_boxes if box[5] == self.get_current_set_index()]
		boxes = [box for box in self.all_boxes if box[5] in self.get_select_indexes()]
		if (len(seg_l) != len(boxes)):
			print("Number of seg images mismatches with number of boxes loaded. Cannot load segs for visualize")
		else:
			sz = int(self.target.bsz_vs.value)
			for i in range(len(seg_l)):
				#region = Region(self.all_boxes[i][0]-sz//2,self.all_boxes[i][1]-sz//2,self.all_boxes[i][2],sz,sz,1)
				region = Region(boxes[i][0]-sz//2,boxes[i][1]-sz//2,boxes[i][2],sz,sz,1)
				#print("Region small",region)
				self.target.write_out(seg_l[i], self.target.seg_path, region)
				annotate = EMData(self.target.seg_path, 0, False, self.target.cur_region)
				#print("Cur Region",self.target.cur_region)
				self.target.img_view.set_data(self.target.get_data(),annotate)
		return


	def set_table_item_changed(self, current, previous):
		return
	def show_boxes(self, set_index):
		return

	def closeEvent(self,event):

		self.target.clear_shapes()
		self.close()

	# def add_boxes(self, size = 64):
	# 	sz = size
	# 	for i in range(len(self.boxes)):
	# 		color=QtGui.QColor(self.setcolors[self.boxes[i][5]%len(self.setcolors)]).getRgbF()
	# 		x,y = int(self.boxes[i][0]),int(self.boxes[i][1])
	# 		if self.boxes[i][3] == 1:
	# 			self.img_view.add_shape("box{}".format(i),EMShape(("rect",color[0],color[1],color[2],x-sz//2,y-sz//2,x+(sz+1)//2,y+(sz+1)//2),2))
	# 	self.img_view.updateGL()

class Merger_Tool(QtWidgets.QWidget):
	def __init__(self,target,maps_lst):
		QtWidgets.QWidget.__init__(self,None)
		self.target = target
		self.maps_lst = maps_lst
		self.maps = []
		for map_file in self.maps_lst:
			if len(map_file) == 0:
				print("Please specify probability map for merging. Abort")
				return
			hdr = EMData(map_file, 0,True)
			if [hdr["nx"],hdr["ny"],hdr["nz"]]!= [self.target.nx,self.target.ny,self.target.nz]:
				print("Mismatch dimension between full data and probablity map ",map_file,". Abort")
				return
			else:
				try:
					self.maps.append(EMData(map_file,0, False, self.target.cur_region))
				except Exception as e:
					print(e)
		self.vss = []
		self.bts = []
		self.wss = []
		self.bts_group = QtWidgets.QButtonGroup()
		self.create_class_bt = QtWidgets.QPushButton("Create Classes")
		self.class_root_le = QtWidgets.QLineEdit("Class")
		self.merger_full_bt = QtWidgets.QPushButton("Merge Full Annotation")
		self.merger_full_le = QtWidgets.QLineEdit()
		self.merger_full_le.setText(self.target.seg_path)
		mt_gbl=QtWidgets.QGridLayout(self)
		for i in range(len(self.maps)):
			self.bts.append(QtWidgets.QPushButton("Show Map "+str(i+1)))
			self.vss.append(ValSlider(value=0.5,rng=(0.001,1.2),rounding=2,label= "Thres "+str(i+1)))
			self.wss.append(ValSlider(value=1,rng=(0.000,1),rounding=2,label= "Weight "+str(i+1)))
			self.bts_group.addButton(self.bts[i],i)
			mt_gbl.addWidget(self.bts[i],i,0,1,1)
			mt_gbl.addWidget(self.vss[i],i,1,1,1)
			mt_gbl.addWidget(self.wss[i],i,2,1,1)
			self.vss[i].valueChanged.connect(self.update_from_slider)
			self.wss[i].valueChanged.connect(self.update_from_slider)
		mt_gbl.addWidget(self.create_class_bt,len(self.maps),0,1,1)
		mt_gbl.addWidget(self.class_root_le,len(self.maps),1,1,1)
		mt_gbl.addWidget(self.merger_full_bt,len(self.maps)+1,0,1,1)
		mt_gbl.addWidget(self.merger_full_le,len(self.maps)+1,1,1,1)
		self.bts_group.buttonClicked[QtWidgets.QAbstractButton].connect(self.on_check_bts_group)
		self.create_class_bt.clicked[bool].connect(self.create_class_from_merge)
		self.merger_full_bt.clicked[bool].connect(self.merge_full)
		self.setWindowTitle("Merger Tool Widget")
		self.show()
		#self.setCentralWidget(QtWidgets.QWidget)

	###TODO: NEED FURTHER TESTING ON UPDATE FROM SLIDER SAVED ORI WHEN TAB CHANGE AND SEL CHANGED
	def update_from_slider(self):
		zero_map = self.target.get_annotation().copy_head()
		zero_map.to_zero()
		self.norm_maps = [to_numpy(zero_map)]
		for i in range(len(self.vss)):
			map = self.maps[i].process("threshold.binary",{"value":self.vss[i].value})
			self.norm_maps.append(self.wss[i].value * to_numpy(map))
		stack_map = np.stack(self.norm_maps)
		self.target.annotate= from_numpy(np.argmax(stack_map,axis=0))
		self.target.img_view.set_data(self.target.data, self.target.annotate, keepcontrast=True)

	def on_check_bts_group(self,bt):
		id = self.bts_group.id(bt)
		map = self.maps[id].process("threshold.binary",{"value":self.vss[id].value})
		self.target.img_view.set_data(self.target.data, (id+1)*map, keepcontrast=True)

	def create_class_from_merge(self):
		name = self.class_root_le.text()
		self.target.get_segtab().tree_set.clear()
		for i in range(len(self.maps)):
			self.target.get_segtab().add_child(parent_item="root",child_l=[str(i+1),name+"_"+str(i+1),"-1"])
			self.target.get_segtab().update_sets()
		#self.close()

	def merge_full(self):
		print ("Start thresholding probability maps.")
		try:
			maps = [EMData(file_name) for file_name in self.maps_lst]
		except Exception as e:
			print(e)
			return
		zero_map = maps[0].copy_head()
		zero_map.to_zero()
		norm_maps = [to_numpy(zero_map)]
		for i in range(len(maps)):
			print ("***")
			map = maps[i].process("threshold.binary",{"value":self.vss[i].value})
			norm_maps.append(self.wss[i].value * to_numpy(map))

		print("Done processing probablity map. Calculating final map.")
		merge_out = from_numpy(np.argmax(np.stack(norm_maps),axis=0))
		try:
			merge_out.write_image(self.merger_full_le.text())
			print("Final map written to ",self.merger_full_le.text(),". Close Merger Tool.")
			del maps, norm_maps, merge_out
		except Exception as e:
			print(e)
		self.close()

class NNet_Tab(QtWidgets.QWidget):
	def __init__(self,target):
		QtWidgets.QWidget.__init__(self,None)
		#self.target=weakref.ref(target)
		self.target = target

		self.bg_button = QtWidgets.QPushButton("Background")
		self.ann_button = QtWidgets.QPushButton("Classes")
		self.bg_button.setCheckable(True)
		self.ann_button.setCheckable(True)

		self.nn_cb_group = QtWidgets.QButtonGroup()
		self.nn_cb_group.addButton(self.bg_button,0)
		self.nn_cb_group.addButton(self.ann_button,1)



		self.train_class_button=QtWidgets.QPushButton("Train NNet")
		self.train_no_iters_sb = StringBox(label="No iters",value="8",showenable=-1)
		self.train_lr_sb = StringBox(label="Learnrate",value="3e-4",showenable=-1)
		self.build_ts_button = QtWidgets.QPushButton("Build Trainset")
		self.vis_ts_button = QtWidgets.QPushButton("Visualize Trainset")

		#self.build_class_sb=StringBox(label="Class ",value="1",showenable=-1)
		self.build_norep_sb=StringBox(label="No reps",value="5",showenable=-1)
		self.apply_button=QtWidgets.QPushButton("Apply NNet")

		self.apply_all_button=QtWidgets.QPushButton("Apply Full")
		self.unet = None


		nnet_gbl=QtWidgets.QGridLayout(self)
		nnet_gbl.setColumnStretch(70,70)
		nnet_gbl.addWidget(self.bg_button,0,0,1,1)
		nnet_gbl.addWidget(self.ann_button,0,1,1,1)
		nnet_gbl.addWidget(self.build_ts_button,1,0,1,1)
		nnet_gbl.addWidget(self.build_norep_sb,1,1,1,1)
		nnet_gbl.addWidget(self.vis_ts_button,2,0,1,1)
		# nnet_gbl.addWidget(self.build_class_sb,3,0,1,1)
		nnet_gbl.addWidget(self.train_class_button,2,1,1,1)

		nnet_gbl.addWidget(self.train_no_iters_sb,3,0,1,1)
		nnet_gbl.addWidget(self.train_lr_sb,4,0,1,1)
		nnet_gbl.addWidget(self.apply_button,5,0,1,1)
		nnet_gbl.addWidget(self.apply_all_button,5,1,1,1)

		self.bg_button.clicked[bool].connect(self.bg_bt_clicked)
		self.ann_button.clicked[bool].connect(self.ann_bt_clicked)

		self.nn_cb_group.buttonClicked[QtWidgets.QAbstractButton].connect(self.on_check_nn_cb_group)
		self.train_class_button.clicked[bool].connect(self.train_class_bt_clicked)
		self.apply_button.clicked[bool].connect(self.apply_bt_clicked)
		self.build_ts_button.clicked[bool].connect(self.build_trainset)
		self.vis_ts_button.clicked[bool].connect(self.visualize_trainset)

	def extract_region(self, iter=3,thresh=0.5):
		reg_list = []
		#datas = to_numpy(self.get_data())
		labels = to_numpy(self.target.get_annotation())

		if len(labels.shape) == 2:#single 2D image
			#datas = np.expand_dims(datas, axis=0)
			labels = np.expand_dims(labels, axis=0)

		for i in range(len(labels)):
			#open_lab=ndi.binary_opening(labels[i],iterations=iter)
			open_lab = labels[i]
			labeled,num = ndi.label(open_lab>thresh)
			cent_mass = ndi.center_of_mass(open_lab,labeled,[i+1 for i in range(num)])
			print("mass at slice", i,cent_mass)
			for pair in cent_mass:
				reg_list.append([pair[1],pair[0],i])

		print(reg_list)
		return reg_list

	def get_selected_item(self):
		sels = self.target.get_segtab().tree_set.selectedItems()
		if len(sels) == 0 or len(sels) >1:
			print("Select a single class or group before quantification")
			return None
		return sels[0]

	def create_organelles_training_set(self):

		sel = self.get_selected_item()
		if sel is None:
			print('Must select class used for training the neural network')
			return
		self.reg_list = self.extract_region()
		datas = self.target.get_data()
		labels = self.target.get_segtab().get_whole_annotate(sel)
		#labels = self.target.get_annotation()

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
		current_class = self.target.get_treeset().currentItem().text(1)
		trainset_out = os.path.join("./particles",base_name(self.target.data_file)+"_"+current_class+"_trainset.hdf")
		os.system("e2tomoseg_buildtrainset.py --buildset --particles_raw=./particles/org_temp.hdf --particles_label=./particles/org_temp_seg.hdf --boxes_negative=./particles/bg_temp.hdf --ncopy={} --trainset_output={} --validset=0.0".format(int(self.build_norep_sb.getValue()),trainset_out))


	def visualize_trainset(self):
		current_class = self.target.get_treeset().currentItem().text(1)
		trainset_path = os.path.join("./particles",base_name(self.target.data_file)+"_"+current_class+"_trainset.hdf")
		if not os.path.exists(trainset_path):
			print("No trainset file detected in the ./particles folder. Create one first")
			return
		else:
			os.system("e2display.py "+trainset_path)
			return
	def bg_bt_clicked(self):
		self.target.basic_tab.setCurrentIndex(3)
		print("Select background patches for training")
		return

	def ann_bt_clicked(self):
		self.target.basic_tab.setCurrentIndex(0)
		print("Annotate objects for training")
		self.target.img_view.show_inspector(5)
		return


	def train_class_bt_clicked(self):
		annotation_out=self.target.img_view.get_full_annotation().copy_head()
		annotation_out.to_zero()
		#trainset_path = './particles/org_temp_trainset.hdf'
		current_class = self.target.get_treeset().currentItem().text(1)
		trainset_path = os.path.join("./particles",base_name(self.target.data_file)+"_"+current_class+"_trainset.hdf")
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
		unet_wout = os.path.join("./neural_nets",base_name(self.target.data_file)+"_nnet.h5")

		self.unet.train_unet(weights_out=unet_wout,batch_sz=50,val_split=0.2,no_epoch = int(self.train_no_iters_sb.getValue()),learnrate=float(self.train_lr_sb.getValue()))
		#print("Done training Unet. Network weights saved to", self.unet_wout)


	def train_all_bt_clicked(self):

		return


	def apply_bt_clicked(self):
		try:
			sel_val = int(self.get_selected_item().text(0))
		except:
			print("Must select a class of features")
			sel_val = 0

		if not self.unet:
			self.unet = UNet()

		#unet_win = os.path.join("./neural_nets",self.target.tomogram_list.currentItem().text()[:-4]+"_nnet.h5")
		unet_win = os.path.join("./neural_nets",base_name(self.target.data_file)+"_nnet.h5")
		if not os.path.exists(unet_win):
			print("No neural networks saved for this tomogram. Train a Unet first.")
			return
		bin_map = self.unet.apply_unet(weights_in=unet_win, tomogram=self.target.img_view.get_full_data())
		bin_map.write_image("./temp_unet.hdf")
		pred_map = sel_val*bin_map.process("threshold.binary",{"value":0.5})
		print("Done applying unet")
		self.target.save_current_state()
		self.target.activateWindow()
		self.target.img_view.full_annotation = pred_map
		self.target.img_view.force_display_update(set_clip=0)
		self.target.img_view.updateGL()
		del bin_map, pred_map
		return

	def on_check_nn_cb_group(self,cb):
		print(cb.text()+" is selected")
		if cb.text() == "Background":
			self.target.img_view.mouse_mode=1
		elif cb.text() =="Classes":
			self.target.img_view.mouse_mode=5
			self.target.img_view.show_inspector(5)
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
		self.morp_convex_bt = QtWidgets.QPushButton("ConvexHull")

		morp_gbl=QtWidgets.QGridLayout(self)
		#morp_gbl.setColumnStretch(70,70)
		morp_gbl.addWidget(QtWidgets.QLabel("No iters"), 0,0,1,1)
		morp_gbl.addWidget(self.morp_n_iters_sp, 0,1,1,1)
		morp_gbl.addWidget(self.morp_close_bt, 1,0,1,1)
		morp_gbl.addWidget(self.morp_open_bt, 1,1,1,1)
		morp_gbl.addWidget(self.morp_erode_bt, 2,0,1,1)
		morp_gbl.addWidget(self.morp_dilate_bt, 2,1,1,1)
		morp_gbl.addWidget(self.morp_convex_bt, 3,0,1,1)

		self.morp_close_bt.clicked[bool].connect(self.do_morp_close)
		self.morp_open_bt.clicked[bool].connect(self.do_morp_open)
		self.morp_erode_bt.clicked[bool].connect(self.do_morp_erode)
		self.morp_dilate_bt.clicked[bool].connect(self.do_morp_dilate)
		self.morp_convex_bt.clicked[bool].connect(self.calc_convex_hull)
		#self.morp_label_bt.clicked[bool].connect(self.do_morp_label)

	def get_target_selected(self):
		sels = self.target.get_segtab().tree_set.selectedItems()
		if len(sels) == 0 or len(sels) >1:
			print("Select a single class or grouping before morphological operation")
			return 0,None
		return self.target.get_segtab().get_whole_branch(sels[0])


	# def calc_convex_hull(self):
	# 	def is_point_in_polygon(point,vertices):
	# 		points = np.array([[x[0],x[1]] for x in vertices])
	# 		poly_path = mplPath.Path(points)
	# 		return poly_path.contains_point(point)
	# 	def polygon(points):
	# 		min_r = int(min([x[0] for x in points]))
	# 		max_r = int(max([x[0] for x in points]))
	# 		min_c = int(min([x[1] for x in points]))
	# 		max_c = int(max([x[1] for x in points]))
	# 		rr = []
	# 		cc = []
	# 		for r in range(min_r,max_r+1):
	# 			for c in range(min_c,max_c+1):
	# 				if is_point_in_polygon([r,c],points):
	# 					rr.append(r)
	# 					cc.append(c)
	# 		return np.array(rr), np.array(cc)
	# 	sels = self.get_target_selected()
	# 	if len(sels) > 1:
	# 		sels = sels[1:]
	# 	self.target.save_current_state()
	# 	self.target.activateWindow()
	# 	for i in range(len(sels)):
	# 		sel = sels[-i]
	# 		val=int(sel.text(0))
	# 		ann=self.target.get_annotation().process("threshold.binaryrange",{"high":val+0.1,"low":val-0.1})
	#
	# 		if ann:
	# 			mask = to_numpy(ann)
	#
	# 			self.target.get_annotation().process_inplace("threshold.rangetozero",{"maxval":(val+0.1),"minval":(val-0.1)})
	# 			p1 = np.where(mask>0)
	#
	# 			points = np.stack((p1[0], p1[1]), axis=-1)
	# 			try:
	# 				hull = ConvexHull(points=points,qhull_options='QG4')
	# 			except Exception as e:
	# 				print(val)
	# 				print("Cannot calc convex hull due to exception", e)
	# 				continue
	# 			pts = []
	# 			for simplex in hull.vertices:
	# 				pts.append([points[simplex,0],points[simplex,1]])
	# 			rr,cc = polygon(pts)
	# 			mask[rr,cc]=val
	# 			self.target.annotate += from_numpy(mask)
	# 			del mask, ann
	# 	self.target.img_view.set_data(self.target.data, self.target.annotate)

	def calc_convex_hull(self):
		def find_bounding_square_box(vertices):
		# Initialize min and max coordinates for x and y
			min_x, min_y = float('inf'), float('inf')
			max_x, max_y = float('-inf'), float('-inf')

			# Iterate through each vertex to find min and max coordinates for x and y
			for x, y,_,_ in vertices:
				min_x, max_x = min(min_x, x), max(max_x, x)
				min_y, max_y = min(min_y, y), max(max_y, y)

			# Determine the length of the side of the square
			side_length = max(max_x - min_x, max_y - min_y) + 2
			if side_length % 2 == 1:
				side_length = side_length+1

			# Define the bottom-left and top-right corners of the square
			bottom_left = (min_x, min_y)
			top_right = (min_x + side_length, min_y + side_length)
			mid = (int(min_x) + side_length//2, int(min_y) + side_length//2)
			return mid,int(side_length)


		def point_in_hull(point, hull, tolerance=1e-12):
			return all((np.dot(eq[:-1], point) + eq[-1] <= tolerance) for eq in hull.equations)

		def generate_square_convex_hull_matrix(points):
		# Calculate the convex hull of the input points
			print("P0",points[0])
			hull = ConvexHull(points)

			# Determine the bounding box of the convex hull
			min_coords = np.floor(np.min(hull.points, axis=0)).astype(int)
			max_coords = np.ceil(np.max(hull.points, axis=0)).astype(int)


			# Find the maximum span in all dimensions to make the matrix square
			span = np.max(max_coords - min_coords)


			# Initialize a square 3D matrix
			matrix_size = span + 1  # Adding 1 to include both ends
			center = min_coords + matrix_size//2
			matrix = np.zeros((matrix_size, matrix_size, matrix_size), dtype=int)

			# Adjust min and max to center the hull in the matrix
			offset = (np.array(matrix.shape) - (max_coords - min_coords)) // 2 - min_coords

			# Fill the matrix
			for point in np.ndindex(matrix.shape):
				adjusted_point = np.array(point) - offset
				if point_in_hull(adjusted_point, hull):
					matrix[point] = 1
			matrix_rearranged = np.transpose(matrix, (2, 0, 1))
			return matrix_rearranged, center
		def set_border_voxels_to_one(matrix):
			# Create a new matrix filled with zeros of the same shape as the input matrix
			border_matrix = np.zeros_like(matrix)
			# Set the border voxels to 1
			# For x dimension
			border_matrix[0, :, :] = 1
			border_matrix[-1, :, :] = 1
			# For y dimension
			border_matrix[:, 0, :] = 1
			border_matrix[:, -1, :] = 1
			# For z dimension
			border_matrix[:, :, 0] = 1
			border_matrix[:, :, -1] = 1
			return border_matrix

		sels = self.get_target_selected()
		if len(sels) > 1:
			sels = sels[1:]


		self.target.save_current_state()
		self.target.activateWindow()

		for i in range(len(sels)):
			sel = sels[-i]
			val=int(sel.text(0))
			border_matrix = set_border_voxels_to_one(np.where(to_numpy(self.target.img_view.get_full_annotation())==val,1,0))
			points = list(zip(*np.where(border_matrix==1)))
			print("Ps Len", len(points))
			fill_vol, center = generate_square_convex_hull_matrix(points)
			side_length = fill_vol.shape[0]


			xrot, yrot, zrot = self.target.img_view.reverse_transform_point([center[0],center[1],self.target.img_view.nz//2+center[2]], self.target.img_view.get_xform())
			xform=Transform({"type":"eman","alt":self.target.img_view.alt,"az":self.target.img_view.az,"tx":self.target.img_view.nx//2+xrot,"ty":self.target.img_view.ny//2+yrot,"tz":self.target.img_view.nz//2+zrot})
			subvol=self.target.img_view.get_full_annotation().get_rotated_clip(xform,(side_length,side_length,side_length),0)


			#subvol = from_numpy(np.ones((2*self.w_line,2*self.w_line,2*self.w_line)))
			#subvol.process_inplace("mask.paint",{"x":self.w_line,"y":self.w_line,"z":self.w_line,"r1":self.w_line,"v1":value,"r2":self.w_line,"v2":0})
			self.target,img_view.full_annotation.set_rotated_clip(xform,from_numpy(np.where(fill_vol>0,val,to_numpy(subvol))))
			#self.img_view.full_annotation.set_rotated_clip(xform,subvol*2)
			# full_vol=from_numpy(np.where(fill_vol>0,value,to_numpy(self.img_view.full_annotation)))
			# self.img_view.full_annotation = full_vol.process("xform",{"transform":self.img_view.xform})
		self.target.img_view.force_display_update()
		self.target.img_view.updateGL()



	def do_morp_close(self):
		n_iters = int(self.morp_n_iters_sp.value())
		sels = self.get_target_selected()
		self.target.save_current_state()
		self.target.activateWindow()
		try:
			for sel in sels:
				val=int(sel.text(0))
				mask=self.target.get_annotation().process("threshold.binaryrange",{"high":val+0.1,"low":val-0.1})
				if mask:
					self.target.get_annotation().process_inplace("threshold.rangetozero",{"maxval":(val+0.1),"minval":(val-0.1)})
					self.target.annotate += val*from_numpy(ndi.binary_closing(to_numpy(mask),iterations=n_iters))

					del mask
			self.target.img_view.set_data(self.target.data, self.target.annotate)
		except Exception as e:
			print(e)


	def do_morp_open(self):


		n_iters = int(self.morp_n_iters_sp.value())
		sels = self.get_target_selected()
		self.target.save_current_state()
		self.target.activateWindow()
		try:
			for sel in sels:
				val=int(sel.text(0))
				mask=self.target.get_annotation().process("threshold.binaryrange",{"high":val+0.1,"low":val-0.1})
				if mask:
					self.target.get_annotation().process_inplace("threshold.rangetozero",{"maxval":(val+0.1),"minval":(val-0.1)})
					self.target.annotate += val*from_numpy(ndi.binary_opening(to_numpy(mask),iterations=n_iters))

					del mask
			self.target.img_view.set_data(self.target.data, self.target.annotate)
		except Exception as e:
			print(e)

	def do_morp_dilate(self):
		n_iters = int(self.morp_n_iters_sp.value())
		sels = self.get_target_selected()
		self.target.save_current_state()
		self.target.activateWindow()
		try:
			for sel in sels:
				val=int(sel.text(0))
				mask=self.target.get_annotation().process("threshold.binaryrange",{"high":val+0.1,"low":val-0.1})
				if mask:
					self.target.get_annotation().process_inplace("threshold.rangetozero",{"maxval":(val+0.1),"minval":(val-0.1)})
					self.target.annotate += val*from_numpy(ndi.binary_dilation(to_numpy(mask),iterations=n_iters))

					del mask
			self.target.img_view.set_data(self.target.data, self.target.annotate)
		except Exception as e:
			print(e)



	def do_morp_erode(self):
		n_iters = int(self.morp_n_iters_sp.value())
		self.target.save_current_state()
		self.target.activateWindow()
		try:
			sels = self.get_target_selected()
			for sel in sels:
				val=int(sel.text(0))
				mask=self.target.get_annotation().process("threshold.binaryrange",{"high":val+0.1,"low":val-0.1})

				if mask:
					self.target.get_annotation().process_inplace("threshold.rangetozero",{"maxval":(val+0.1),"minval":(val-0.1)})
					self.target.annotate += val*from_numpy(ndi.binary_erosion(to_numpy(mask),iterations=n_iters))
					del mask
			self.target.img_view.set_data(self.target.data, self.target.annotate)
		except Exception as e:
			print(e)

	def do_morp_label(self):
		sel = self.get_target_selected()
		val = int(sel.text(0))
		mask = self.target.get_segtab().get_whole_annotate(sel)

		if raw_mask:
			self.target.get_annotation().process_inplace("threshold.rangetozero",{"maxval":(val+0.1),"minval":(val-0.1)})
			mask, num = ndi.label(to_numpy(raw_mask))
			t_mask = np.zeros(mask.shape)
			for i in range(1,num+1):
				ind = self.target.get_segtab().get_unused_index()
				name = self.target.get_segtab().tree_set.currentItem().text(1)
				self.target.get_segtab().add_child(child_l=[str(ind),name+"_"+str(i),"-1"])
				self.target.get_segtab().update_sets()
				t_mask += np.where(mask==i,ind,0)

			self.target.annotate += from_numpy(t_mask)
			print("number of object detected:", num)
			self.target.img_view.set_data(self.target.data, self.target.annotate)
			del raw_mask,mask


class Extrapolate_Tab(QtWidgets.QWidget):
	def __init__(self,target) :
		super().__init__()
		self.target = target
		ext_gbl = QtWidgets.QGridLayout(self)
		self.current_z = self.target.zc_spinbox.value()
		# self.rf_path_text = QtWidgets.QLineEdit()
		# self.vol_path_text = QtWidgets.QLineEdit()
		#
		# self.vol_browser_bt = QtWidgets.QPushButton("Browse")
		# self.map_ptcls_bt = QtWidgets.QPushButton("Map ptcls to tomogram")
		# self.show_ptcls_bt = QtWidgets.QPushButton("Show ptcls on tomogram")
		# self.show_ptcls_vs = ValSlider(value=0.5,rng=(0.001,5),rounding=2,label= "Thres")
		#
		self.ext_class_bt = QtWidgets.QPushButton("Extrapolate Class")
		self.ext_all_bt = QtWidgets.QPushButton("Extrapolate All")
		self.ext_from_spinbox = QtWidgets.QSpinBox()
		self.ext_from_spinbox.setMaximum(self.target.get_nz())
		self.ext_from_spinbox.setMinimum(-self.target.get_nz())
		self.ext_from_spinbox.setValue(0)



		# self.n_iters_spinbox.setMinimum(1)
		# self.annotate_ori=None

		ext_gbl.addWidget(QtWidgets.QLabel("Apply to 2D slice. Z-thick set to 0"),0,0,1,2)
		ext_gbl.addWidget(self.ext_class_bt,1,0,1,1)
		ext_gbl.addWidget(self.ext_all_bt,1,1,1,1)
		ext_gbl.addWidget(QtWidgets.QLabel("Extrapolate From"),2,0,1,1)
		ext_gbl.addWidget(self.ext_from_spinbox,2,1,1,1)


		self.ext_class_bt.clicked[bool].connect(self.extrapolate_classes)
		self.ext_all_bt.clicked[bool].connect(self.extrapolate_all)

	def set_target_z(self,value):
		self.target.zt_spinbox.setValue(value)
		print("Setting z-thick to 0, process on a 2D slice at z={}".format(self.target.zc_spinbox.value()))

	def get_target_selected(self):
		sels = []
		for sel in self.target.get_segtab().tree_set.selectedItems():
			for item in self.target.get_segtab().get_whole_branch(sel):
				sels.append(item)

		# if len(sels) == 0:
		# 	print("Select a single class or grouping before extrapolate operation")
		# 	return 0,None

		return sels

	def get_annotate_at_z(self,zcen,sels=[],all=False):

		print("zcen", zcen)
		target_region = Region(self.target.cur_region.x_origin(),self.target.cur_region.x_origin(),zcen, self.target.cur_region.get_width(), self.target.cur_region.get_height(),1)
		print(target_region)
		target_annotate = EMData(self.target.seg_path, 0, False, target_region)
		if all:
			return target_annotate
		else:
			annotation_out=self.target.get_annotation()
			temp = annotation_out.copy_head()

			for ind in sels:
				print("Sel ind", ind)
				num = int(ind.text(0))
				annotation_out.process_inplace("threshold.rangetozero",{"maxval":(num+0.1),"minval":(num-0.1)})
				temp += num*(target_annotate.process("threshold.binaryrange",{"high":num+0.1,"low":num-0.1}))

			return annotation_out + temp


	def extrapolate_classes(self):
		self.set_target_z(0)
		zcen = self.target.zc_spinbox.value()+self.ext_from_spinbox.value()
		self.target.annotate = self.get_annotate_at_z(zcen = zcen,sels=self.get_target_selected())
		self.target.img_view.set_data(self.target.data,self.target.annotate)
		return



	def extrapolate_all(self):
		self.set_target_z(0)
		zcen = self.target.zc_spinbox.value()+self.ext_from_spinbox.value()
		self.target.annotate =  self.get_annotate_at_z(zcen = zcen,all=True)
		self.target.img_view.set_data(self.target.data,self.target.annotate)
		return


class Binary_Tab(QtWidgets.QWidget):
	def __init__(self,target) :
		QtWidgets.QWidget.__init__(self,None)
		self.target = target
		bin_gbl = QtWidgets.QGridLayout(self)
		self.bin_invert_cb = QtWidgets.QCheckBox("Dark Feature")
		self.bin_invert_cb.setChecked(True)
		#valslider set to 0 will paint all the annotate file
		self.bin_detect_bt = QtWidgets.QPushButton("Detect Feature")
		self.bin_fill_bt = QtWidgets.QPushButton("Fill")
		self.bin_trim_bt = QtWidgets.QPushButton("Trim")
		self.bin_low_pass_vs = ValSlider(value=5,rng=(0.001,5),rounding=2,label= "Cut-off Abs")
		self.bin_threshold_vs = ValSlider(value=0.001,rng=(0.001,5),rounding=2,label="Threshold  ")
		self.closing_n_iters =1
		self.opening_n_iters =1
		self.saved_ori = False
		self.quiet = False

		bin_gbl.addWidget(self.bin_invert_cb,0,0,1,1)
		bin_gbl.addWidget(self.bin_detect_bt,0,1,1,1)
		bin_gbl.addWidget(self.bin_low_pass_vs,1,0,1,2)
		bin_gbl.addWidget(self.bin_threshold_vs,2,0,1,2)
		bin_gbl.addWidget(self.bin_fill_bt,3,0,1,1)
		bin_gbl.addWidget(self.bin_trim_bt,3,1,1,1)
		self.bin_detect_bt.clicked[bool].connect(self.bin_detect_bt_clicked)
		self.bin_fill_bt.clicked[bool].connect(self.fill_hole_bin)
		self.bin_trim_bt.clicked[bool].connect(self.do_area_opening)
		self.bin_low_pass_vs.valueChanged.connect(self.update_mask_from_vs)
		self.bin_threshold_vs.valueChanged.connect(self.update_mask_from_vs)

	def get_selected_item(self):
		sels = self.target.get_segtab().tree_set.selectedItems()
		if len(sels) == 0 or len(sels) >1:
			print("Select a single class or group before quantification")
			return
		return sels[0]

	def do_area_opening(self):
		self.target.morp_tab.morp_n_iters_sp.setValue(self.opening_n_iters)
		self.target.morp_tab.do_morp_open()
		self.opening_n_iters +=1

	def fill_hole_bin(self):
		sel = self.get_selected_item()
		val=int(sel.text(0))
		self.target.save_current_state()
		self.target.activateWindow()
		mask=self.target.get_annotation().process("threshold.binaryrange",{"high":val+0.1,"low":val-0.1})
		if mask:
			self.target.get_annotation().process_inplace("threshold.rangetozero",{"maxval":(val+0.1),"minval":(val-0.1)})
			self.target.annotate += val*from_numpy(ndi.binary_fill_holes(to_numpy(mask)))
			del mask
		self.target.img_view.set_data(self.target.data, self.target.annotate)

	def bin_detect_bt_clicked(self):
		self.target.save_current_state()

		sel = self.get_selected_item()
		if not sel:
			print("Choose only one class for detect feature")
			return
		# if len(sel) == 0:
		# 	print("Please select the class feature to detect")
		# 	return


		self.target.get_segtab().get_whole_annotate(sel).write_image("ori_detect.hdf")
		self.saved_ori = True

		self.target.activateWindow()
		lp_v = self.target.get_inspector().maxs.value*0.6
		thres_v = 0.03

		self.bin_low_pass_vs.setValue(lp_v)
		self.bin_threshold_vs.setValue(thres_v)

		self.target.reset_morp_params()
		self.quiet = False
		self.update_mask_from_vs()

	def update_mask_from_vs(self):
		if not self.saved_ori:
			sel = self.get_selected_item()
			self.target.get_segtab().get_whole_annotate(sel).write_image("ori_detect.hdf")
		else:
			pass
		self.closing_n_iters =1
		self.opening_n_iters =1
		sel = self.get_selected_item()

		if self.bin_invert_cb.isChecked():
			mult = -1
		else:
			mult= 1
		if sel:
			val = int(sel.text(0))
		else:
			return

		self.mask = EMData("ori_detect.hdf")
		lp = self.bin_low_pass_vs.value
		thres = self.bin_threshold_vs.value
		lp_masked = mult*self.mask*self.target.get_data().process("filter.lowpass.gauss",{"cutoff_abs":lp})
		mask = lp_masked.process("threshold.binary",{"value":thres})
		self.target.annotate *= (1-self.mask)
		self.target.annotate += val*mask
		self.target.img_view.set_data(self.target.data, self.target.annotate)
		return


class Templ_Match_Tab(QtWidgets.QWidget):
	def __init__(self,target) :
		QtWidgets.QWidget.__init__(self,None)
		self.target = target
		self.template_text = QtWidgets.QLineEdit()
		self.template_browser_bt = QtWidgets.QPushButton("Browse")
		self.create_template_bt = QtWidgets.QPushButton("Make Template")
		self.template_match_bt = QtWidgets.QPushButton("Template Match")
		self.tplt_low_pass_vs = ValSlider(value=1,rng=(0.001,1),rounding=2,label= "Low-pass Filt")
		self.tplt_threshold_vs = ValSlider(value=0.8,rng=(0.001,6),rounding=2,label="Binary Thresh")


		templ_gbl = QtWidgets.QGridLayout()
		templ_gbl.addWidget(self.template_text,0,0,1,3)
		templ_gbl.addWidget(self.template_browser_bt,0,3,1,1)
		templ_gbl.addWidget(self.create_template_bt,1,0,1,2)
		templ_gbl.addWidget(self.template_match_bt,1,2,1,2)
		templ_gbl.addWidget(self.tplt_threshold_vs,2,0,1,4)
		#templ_gbl.addWidget(self.tplt_low_pass_vs,3,0,1,4)
		self.setLayout(templ_gbl)

		self.tplt_match_launcher = QtWidgets.QWidget()
		self.tplt_match_launcher.setWindowTitle("Command for template-matching")
		self.tplt_match_launcher.setMinimumWidth(500)
		self.tplt_match_cmd = QtWidgets.QTextEdit()
		self.tplt_launch_bt = QtWidgets.QPushButton("Launch")

		tml_layout = QtWidgets.QVBoxLayout()
		tml_layout.addWidget(self.tplt_match_cmd)
		tml_layout.addWidget(self.tplt_launch_bt)
		self.tplt_match_launcher.setLayout(tml_layout)

		self.template_browser_bt.clicked[bool].connect(self.load_template)
		self.template_match_bt.clicked[bool].connect(self.tplt_match)
		self.create_template_bt.clicked[bool].connect(self.create_tplt)
		self.tplt_launch_bt.clicked[bool].connect(self.tplt_match_launched)
		self.tplt_threshold_vs.valueChanged.connect(self.binarize_tplt_match)




	def load_template(self):
		self.tplt_browser = EMBrowserWidget(withmodal=True,multiselect=False)
		self.tplt_browser.ok.connect(self.tplt_browser_ok)
		self.tplt_browser.show()
		return

	def tplt_browser_ok(self):
		self.tplt_browser_ret = (self.tplt_browser.getResult())
		self.template_text.setText(self.tplt_browser_ret[0])
		self.template = EMData(self.tplt_browser_ret[0])
		#process template

	def tplt_match(self):
		self.tm_cmd="e2spt_tempmatch.py "+self.target.data_file+" --reference="+self.template_text.text()+" --nptcl=500 --dthr=-1.0 --vthr=0.2 --minvol=-1 --maxvol=-1 --delta=90.0 --sym=c1 --rmedge --rmgold --boxsz=-1 --threads=4 --shrink=1"
		self.tplt_match_cmd.setText(self.tm_cmd)
		self.tplt_match_launcher.show()
		return

	def tplt_match_launched(self):
		self.tplt_match_launcher.close()
		try:
			os.system(self.tplt_match_cmd.toPlainText())
		except:
			pass
		os.system("e2proc3d.py tmp_ccc.hdf tmp_ccc_inv.hdf --mult=-1")
		#tplt_match_map = EMData("tmp_ccc_inv.hdf")
		self.tplt_threshold_vs.setValue(0.81)
		return

	def binarize_tplt_match(self, event):
		tplt_match_map = EMData("tmp_ccc_inv.hdf")
		self.target.save_current_state()
		self.target.activateWindow()
		tplt_match_map_bin=tplt_match_map.process("threshold.binary",{"value":self.tplt_threshold_vs.value})
		self.target.img_view.full_annotation =tplt_match_map_bin
		self.target.img_view.force_display_update(set_clip=0)
		self.target.img_view.updateGL()

	def create_tplt(self):
		return

class Cmd_Launcher(QtWidgets.QWidget):
	def __init__(self,title="") :
		super().__init__(self,None)
		self.setWindowTitle(title)
		self.setMinimumWidth(500)
		self.cmd= QtWidgets.QTextEdit()
		self.launch_bt = QtWidgets.QPushButton("Launch")

		vbl = QtWidgets.QVBoxLayout()
		vb.addWidget(self.cmd)
		vb.addWidget(self.launch_bt)
		self.setLayout(vbl)

	def get_cmd():
		return self.cmd






class Fila_Tab(QtWidgets.QWidget):
	def __init__(self,target) :
		QtWidgets.QWidget.__init__(self,None)
		target = target
		self.tree = QtWidgets.QTreeWidget(self)
		self.tree.setColumnCount(3)
		self.tree.setHeaderLabels(["a","b","c"])
		self.itemflags = Qt.ItemFlags(Qt.ItemIsEditable)|Qt.ItemFlags(Qt.ItemIsSelectable)|Qt.ItemFlags(Qt.ItemIsEnabled)|Qt.ItemFlags(Qt.ItemIsUserCheckable)
		self.tree_root = self.tree.invisibleRootItem()
		data = {'Project A': {'file_a.py': {'cool': {'cool': None}}, 'file_a.txt': None, 'something.xls': None}, 'Project B': {'file_b.csv': None, 'photo.jpg': None}, 'Project C': None}


		def fill_item(item, value):
			def new_item(parent, text, val=None):
				child = QtWidgets.QTreeWidgetItem([text])
				fill_item(child, val)
				parent.addChild(child)
				#child.setExpanded(True)
			if value is None: return
			elif isinstance(value, dict):
				for key, val in sorted(value.items()):
					new_item(item, str(key), val)
			elif isinstance(value, (list, tuple)):
				for val in value:
					text = (str(val) if not isinstance(val, (dict, list, tuple))
							else '[%s]' % type(val).__name__)
					new_item(item, text, val)
			else:
				new_item(item, str(value))
		fill_item(self.tree.invisibleRootItem(),data)
		self.t_button = QtWidgets.QPushButton('test')
		self.p_button = QtWidgets.QPushButton('print')
		self.c_button = QtWidgets.QPushButton('child')
		self.gbl = QtWidgets.QGridLayout(self)
		self.gbl.addWidget(QtWidgets.QLabel("This panel is currently for testing only."),0,0,1,3)
		self.gbl.addWidget(self.tree,1,0,5,3)

		self.gbl.addWidget(self.t_button,0,3,1,1)
		self.gbl.addWidget(self.p_button,1,3,1,1)
		self.gbl.addWidget(self.c_button,2,3,1,1)

		self.t_button.clicked[bool].connect(self.t_button_clicked)
		self.p_button.clicked[bool].connect(self.p_button_clicked)
		self.c_button.clicked[bool].connect(self.c_button_clicked)




	def tree_to_dict(self,parent):
		childCount = parent.childCount()
		if not childCount:
			return
		content = {}
		for row in range(childCount):
			child = parent.child(row)
			content[child.text(0)] = self.tree_to_dict(child)


			values = json.loads(file_in["ann_name"])
			serialize_name = json.dumps(names, default=lambda a: "[%s,%s]" % (str(type(a)), a.pk))

		return content

	def t_button_clicked(self):
		self.tree.insertTopLevelItem(0, QtWidgets.QTreeWidgetItem(['cool name']))
	def p_button_clicked(self):
		json_str = self.tree_to_dict(self.tree.invisibleRootItem())
		print(json_str)

	def change_parent(self,item, new_parent):
		old_parent = item.parent()
		print("Old parent",old_parent)
		try:
			ix = old_parent.indexOfChild(item)
			item_without_parent = old_parent.takeChild(ix)
		except:
			self.tree_root.addChild(item)
			ix = self.tree_root.indexOfChild(item)
			item_without_parent = self.tree_root.takeChild(ix)
		new_parent.addChild(item_without_parent)
	def c_button_clicked(self):
		child = QtWidgets.QTreeWidgetItem(['cool', 'verycool'])
		self.tree.currentItem().addChild(child)


class Subtom_Tab(QtWidgets.QWidget):
	def __init__(self,target) :
		QtWidgets.QWidget.__init__(self,None)
		self.target = target

		self.rf_path_text = QtWidgets.QLineEdit()
		self.vol_path_text = QtWidgets.QLineEdit()

		self.vol_browser_bt = QtWidgets.QPushButton("Browse")
		self.map_ptcls_bt = QtWidgets.QPushButton("Map ptcls to tomogram")
		self.show_ptcls_bt = QtWidgets.QPushButton("Show ptcls on tomogram")
		self.show_ptcls_vs = ValSlider(value=0.5,rng=(0.001,5),rounding=2,label= "Thres")

		self.n_iters_spinbox = QtWidgets.QSpinBox()
		self.n_iters_spinbox.setValue(4)
		self.n_iters_spinbox.setMinimum(1)
		self.annotate_ori=None


		subtom_gbl = QtWidgets.QGridLayout()
		subtom_gbl.addWidget(QtWidgets.QLabel("Refinement folder path"),0,0,1,1)
		subtom_gbl.addWidget(self.rf_path_text,0,1,1,3)
		subtom_gbl.addWidget(QtWidgets.QLabel("3D volume"),1,0,1,1)
		subtom_gbl.addWidget(self.vol_path_text,1,1,1,2)
		subtom_gbl.addWidget(self.vol_browser_bt,1,3,1,1)
		subtom_gbl.addWidget(QtWidgets.QLabel("No iters"),2,0,1,1)
		subtom_gbl.addWidget(self.n_iters_spinbox,2,1,1,1)
		subtom_gbl.addWidget(self.map_ptcls_bt,2,2,1,2)
		subtom_gbl.addWidget(self.show_ptcls_vs,3,0,1,2)
		subtom_gbl.addWidget(self.show_ptcls_bt,3,2,1,2)
		self.setLayout(subtom_gbl)

		self.map_ptcls_launcher = QtWidgets.QWidget()
		self.map_ptcls_launcher.setWindowTitle("Command for mapping ptcles to current tomogram")
		self.map_ptcls_launcher.setMinimumWidth(500)
		self.map_ptcls_cmd = QtWidgets.QTextEdit()
		self.map_ptcls_launch_bt = QtWidgets.QPushButton("Launch")

		mpl_layout = QtWidgets.QVBoxLayout()
		mpl_layout.addWidget(self.map_ptcls_cmd)
		mpl_layout.addWidget(self.map_ptcls_launch_bt)
		self.map_ptcls_launcher.setLayout(mpl_layout)

		self.vol_browser_bt.clicked[bool].connect(self.load_vol)
		self.map_ptcls_bt.clicked[bool].connect(self.show_map_launcher)
		self.show_ptcls_bt.clicked[bool].connect(self.show_ptcls_on_tom)
		self.map_ptcls_launch_bt.clicked[bool].connect(self.map_ptcls_launched)
		self.show_ptcls_vs.valueChanged.connect(self.show_ptcls_vs_changed)

	def load_vol(self):
		self.mp_browser = EMBrowserWidget(withmodal=True,multiselect=False)
		self.mp_browser.ok.connect(self.load_vol_browser_ok)
		self.mp_browser.show()
		return

	def load_vol_browser_ok(self):
		self.mp_browser_ret = (self.mp_browser.getResult())
		self.vol_path_text.setText(self.mp_browser_ret[0])
		#process template

	def show_map_launcher(self):
		if len(self.vol_path_text.text()) == 0:
			if len(self.rf_path_text.text()) == 0:
				print("Must specify 3D volume or spt_xx folder contain 3D volume.")

			else:
				path = self.rf_path_text.text()
				itr = self.n_iters_spinbox.value()
				self.vol_path_text.setText("{}/threed_{:02d}.hdf".format(path, itr))
		else:
			pass
		if len(self.rf_path_text.text()) == 0:
			self.path = "."
		else:
			self.path = self.rf_path_text.text()
		self.mp_cmd="e2spt_mapptclstotomo.py  --tomo={}  --path={}  --avg={} --iter={} --new ".format(self.target.data_file,self.path,self.vol_path_text.text(),self.n_iters_spinbox.value())
		self.map_ptcls_cmd.setText(self.mp_cmd)
		self.map_ptcls_launcher.show()
		return

	def map_ptcls_launched(self,new=True):
		self.vol = EMData(self.vol_path_text.text())
		self.annotate_ori =  self.target.annotate.copy()
		self.map_ptcls_launcher.close()
		try:
			os.system(self.map_ptcls_cmd.toPlainText())
		except:
			print("Error launching e2spt_mapptclstotomo.py program. Abort.")
			pass
		return

	def get_selected_item(self):
		sels = self.target.get_segtab().tree_set.selectedItems()
		if len(sels) == 0 or len(sels) >1:
			print("Select a single class or group for display group")
			return None
		return sels[0]

	def show_ptcls_on_tom(self):
		sel = self.get_selected_item()
		if sel:
			val = int(sel.text(0))
		if len(self.rf_path_text.text()) == 0:
			self.path = "."
		else:
			self.path = self.rf_path_text.text()
		tomo_name = base_name(self.target.data_file)
		self.mask_path = "{}/ptcls_in_tomo_{}_{:02d}.hdf".format(self.path, tomo_name, self.n_iters_spinbox.value())
		if not (os.path.isfile(self.mask_path)):
			print("Map of particles to the tomogram does not exist. Run Map ptcls to tomogram first")
			return
		else:
			self.target.save_current_state()
			self.target.activateWindow()
			if self.annotate_ori is None:
				self.annotate_ori = self.target.annotate.copy()
			try:

				self.mask = EMData(self.mask_path, 0, False, self.target.cur_region).process("threshold.binary",{"value":self.show_ptcls_vs.value})
				self.target.annotate = self.annotate_ori*(1-self.mask) + self.mask*val
				self.target.img_view.set_data(self.target.data,self.target.annotate)
			except Exception as e:
				print(e)
			return

	def show_ptcls_vs_changed(self,val):
		self.show_ptcls_on_tom()
		return

class Statistics_Tab(QtWidgets.QWidget):
	def __init__(self,target) :
		QtWidgets.QWidget.__init__(self,None)
		self.target = target
		stat_gbl = QtWidgets.QVBoxLayout(self)
		self.blob_tab = QtWidgets.QWidget()
		bltlay = QtWidgets.QGridLayout(self.blob_tab)
		self.counted_item = []
		self.n_obj_thres_vs = ValSlider(value=10,rng=(0.001,5000),rounding=2,label= "Area/Vol Thresh")
		self.label_objs_bt = QtWidgets.QPushButton("Label objects")
		self.n_objects_text = QtWidgets.QLabel("---")
		self.stat_bt = QtWidgets.QPushButton("Calc Stat")
		self.convex_bt = QtWidgets.QPushButton("Convex Hull")
		self.apix_vs =  ValSlider(value=1.0,rng=(0.001,5.0),rounding=2,label= "Angstrom/Pixel")
		self.apix_vs.setValue(self.target.apix)



		self.stat_combo = QtWidgets.QComboBox()
		self.stat_combo.addItem('Center of Mass')
		self.stat_combo.addItem('Volume/Area')
		self.stat_combo.addItem('Largest Area')
		self.stat_combo.addItem('Largest Perimeter')
		self.stat_combo.addItem('Custom')

		bltlay.addWidget(self.n_obj_thres_vs,0,0,1,3)
		bltlay.addWidget(QtWidgets.QLabel("Num objects larger than thresh:"),1,0,1,3)
		bltlay.addWidget(self.n_objects_text,1,2,1,1)
		bltlay.addWidget(self.label_objs_bt,2,0,1,1)

		bltlay.addWidget(self.apix_vs,3,0,1,3)
		bltlay.addWidget(self.stat_bt,4,0,1,1)
		bltlay.addWidget(self.stat_combo,4,1,1,2)


		self.fila_tab = Fila_Tab(target=self)
		fillay = QtWidgets.QGridLayout(self.blob_tab)

		self.stat_tab = QtWidgets.QTabWidget()
		self.stat_tab.addTab(self.blob_tab,"Blob-like")
		self.stat_tab.addTab(self.fila_tab,"Filament-like")
		stat_gbl.addWidget(self.stat_tab)

		self.stat_bt.clicked[bool].connect(self.calc_stat)

		#self.convex_bt.clicked[bool].connect(self.show_convex_hull)
		self.n_obj_thres_vs.valueChanged.connect(self.thres_vs_changed)
		self.label_objs_bt.clicked[bool].connect(self.label_objs)


	def get_selected_item(self):
		sels = self.target.get_segtab().tree_set.selectedItems()
		if len(sels) == 0 or len(sels) >1:
			print("Select a single class or group before quantification")
			return None
		return sels[0]
		#return self.target.get_segtab().get_whole_branch(sels[0])



	def do_morp_close(self):
		n_iters = int(self.morp_n_iters_sp.value())
		sels = self.get_target_selected()
		for sel in sels:
			val=int(sel.text(0))
			mask=self.target.get_annotation().process("threshold.binaryrange",{"high":val+0.1,"low":val-0.1})
			if mask:
				self.target.get_annotation().process_inplace("threshold.rangetozero",{"maxval":(val+0.1),"minval":(val-0.1)})
				self.target.annotate += val*from_numpy(ndi.binary_closing(to_numpy(mask),iterations=n_iters))

				del mask
		self.target.img_view.set_data(self.target.data, self.target.annotate)


	#WORK BUT NEED TO BE CONDENSED
	def thres_vs_changed(self):
		thres=self.n_obj_thres_vs.value
		try:
			sels = self.target.get_segtab().get_whole_branch(self.get_selected_item())
		except:
			print ("Must select a class to count objects.")
			return
		if len(sels) == 1:
			sel = sels[0]
			val=int(sel.text(0))
			raw_mask=to_numpy(self.target.get_annotation().process("threshold.binaryrange",{"high":val+0.1,"low":val-0.1}))
			self.target.get_annotation().process_inplace("threshold.rangetozero",{"maxval":(val+0.1),"minval":(val-0.1)})
			mask_op = ndi.binary_opening(raw_mask, iterations=1)
			labeled, num = ndi.label(mask_op)
			areas = np.array(ndi.sum(mask_op,labeled, np.arange(labeled.max()+1)))
			mask = areas > thres
			remove_small = mask[labeled.ravel()].reshape(labeled.shape)
			labeled, num = ndi.label(mask_op)
			total = num
			self.target.annotate += val*from_numpy(remove_small)
			del raw_mask, mask_op, labeled, num, areas, mask
		else:
			total= 0
			for sel in sels[1:]:
				val=int(sel.text(0))
				raw_mask=to_numpy(self.target.get_annotation().process("threshold.binaryrange",{"high":val+0.1,"low":val-0.1}))
				self.target.get_annotation().process_inplace("threshold.rangetozero",{"maxval":(val+0.1),"minval":(val-0.1)})
				mask_op = ndi.binary_opening(raw_mask, iterations=1)
				labeled, num = ndi.label(mask_op)
				areas = np.array(ndi.sum(mask_op,labeled, np.arange(labeled.max()+1)))
				mask = areas > thres
				remove_small = mask[labeled.ravel()].reshape(labeled.shape)
				labeled, num = ndi.label(mask_op)
				total += num
				self.target.annotate += val*from_numpy(remove_small)
				del raw_mask, mask_op, labeled, num, areas, mask
				#self.target.annotate *= from_numpy(raw_mask-remove_small)
		self.n_objects_text.setText(str(total))
		self.target.img_view.set_data(self.target.data, self.target.annotate)


	def label_objs(self,ask_user=True):
		sel = self.get_selected_item()
		try:
			val = int(sel.text(0))
			sels = self.target.get_segtab().get_whole_branch(sel)
		except:
			print ("Must select a class to count objects.")
			return
		raw_mask = self.target.get_segtab().get_whole_annotate(sel)
		relabel = False
		if len(sels)>1:
			if not ask_user:
				relabel = True
			else:
				msg="This class is already labeled into objects. Do you want to redo the labeling?"
				result = QtWidgets.QMessageBox.question(self,"Relabel",msg,QtWidgets.QMessageBox.Yes|QtWidgets.QMessageBox.No)
				if result == QtWidgets.QMessageBox.Yes:
					relabel = True
					pass
				else:
					return
		else:
			pass
		if relabel:
			bg_in_f = 1-raw_mask
			self.target.annotate *= bg_in_f
			self.target.annotate += val*raw_mask
			sel.takeChildren()
			self.target.get_segtab().update_sets()
		self.target.save_current_state()
		self.target.activateWindow()
		self.target.get_annotation().process_inplace("threshold.rangetozero",{"maxval":(val+0.1),"minval":(val-0.1)})
		# 	# 	self.target.annotate += val*from_numpy(ndi.binary_opening(to_numpy(mask),iterations=n_iters))
		self.labeled_ann,self.num = ndi.label(to_numpy(raw_mask))
		self.loc=ndi.find_objects(self.labeled_ann,self.num)
		t_mask = np.zeros(self.labeled_ann.shape)
		count = 0
		self.area_vol = []
		self.objs = []
		open_lab=to_numpy(raw_mask)
		for i in range(1,self.num+1):
			area_temp=len(np.where(open_lab[self.loc[i-1]]>0)[0])
			count = count+1
			self.area_vol.append(area_temp)
			self.objs.append(open_lab[self.loc[i-1]])
			ind = self.target.get_segtab().get_unused_index()
			name = sel.text(1)
			self.target.get_segtab().add_child(child_l=[str(ind),name+"_"+str(i),"-1"])
			self.target.get_segtab().update_sets()
			#self.target.annotate += ind*(raw_mask)
			t_mask += np.where(self.labeled_ann==i,ind,0)
		self.target.annotate += from_numpy(t_mask)
		print("number of object detected:", self.num)
		self.target.img_view.set_data(self.target.data, self.target.annotate)
		del t_mask
		self.n_objects_text.setText(str(count))


	# def count_objs(self):
	# 	thres=self.n_obj_thres_vs.value
	# 	#n_iters = int(self.morp_n_iters_sp.value())
	# 	sel = self.get_selected_item()
	# 	#self.counted_item.append(sel)
	# 	#val,raw_mask = self.get_target_selected()
	# 	if sel:
	# 		val = int(sel.text(0))
	# 		raw_mask = self.target.get_segtab().get_whole_annotate(sel)
	# 	else:
	# 		print("Select a class to count")
	# 		return
	#
	# 	# if raw_mask:
	# 	self.target.get_annotation().process_inplace("threshold.rangetozero",{"maxval":(val+0.1),"minval":(val-0.1)})
	# 	# 	self.target.annotate += val*from_numpy(ndi.binary_opening(to_numpy(mask),iterations=n_iters))
	# 	# mask, num = ndi.label(to_numpy(self.target.get_annotation()))
	# 	self.labeled_ann,self.num = ndi.label(to_numpy(raw_mask))
	# 	self.loc=ndi.find_objects(self.labeled_ann,self.num)
	# 	t_mask = np.zeros(self.labeled_ann.shape)
	# 	count = 0
	# 	self.area_vol = []
	# 	self.objs = []
	# 	open_lab=to_numpy(raw_mask)
	# 	# current_item = self.target.get_segtab().tree_set.currentItem()
	# 	# current_item.takeChildren()
	#
	# 	for i in range(1,self.num+1):
	# 		#area_temp=len(np.where(open_lab[self.loc[i-1]]>0)[0])
	# 		area_temp=len(np.where(open_lab[self.loc[i-1]]>0)[0])
	# 		#print(open_lab[loc[i]].shape[0],open_lab[loc[i]].shape[1])
	# 		if area_temp >= thres:
	# 			count = count+1
	# 			self.area_vol.append(area_temp)
	# 			self.objs.append(open_lab[self.loc[i-1]])
	# 			ind = self.target.get_segtab().get_unused_index()
	# 			name = current_item.text(1)
	# 			self.target.get_segtab().add_child(child_l=[str(ind),name+"_"+str(i),"-1"])
	# 			self.target.get_segtab().update_sets()
	# 	#self.target.annotate += *(raw_mask)
	# 			t_mask += np.where(self.labeled_ann==i,ind,0)
	#
	# 	self.target.annotate += from_numpy(t_mask)
	# 	print("number of object detected:", self.num)
	# 	self.target.img_view.set_data(self.target.data, self.target.annotate)
	# 	del t_mask
	# 	self.n_objects_text.setText(str(count))



	# def calc_convex_hull(self):
	# 	def is_point_in_polygon(point,vertices):
	# 		points = np.array([[x[0],x[1]] for x in vertices])
	# 		poly_path = mplPath.Path(points)
	# 		return poly_path.contains_point(point)
	# 	def polygon(points):
	# 		min_r = int(min([x[0] for x in points]))
	# 		max_r = int(max([x[0] for x in points]))
	# 		min_c = int(min([x[1] for x in points]))
	# 		max_c = int(max([x[1] for x in points]))
	# 		rr = []
	# 		cc = []
	# 		for r in range(min_r,max_r+1):
	# 			for c in range(min_c,max_c+1):
	# 				if is_point_in_polygon([r,c],points):
	# 					rr.append(r)
	# 					cc.append(c)
	# 		return np.array(rr), np.array(cc)
	#
	#
	#
	# 	open_lab=ndi.binary_opening(to_numpy(img),iterations=1)
	# 	labeled,num = ndi.label(open_lab>0.5)
	#
	# 	for i in range(1,num+1):
	# 		o1 = np.ma.masked_where(labeled!=i,labeled)
	# 		p1 = np.nonzero(o1)
	# 		points = np.stack((p1[0], p1[1]), axis=-1)
	# 		hull = ConvexHull(points=points,qhull_options='QG4')
	# 		pts = []
	# 		for simplex in hull.vertices:
	# 			pts.append([points[simplex,0],points[simplex,1]])
	# 		rr,cc = polygon(pts)
	# 		labeled[rr,cc]=i
	# 	return labeled
	#
	# def show_convex_hull(self):
	# 	sel = self.get_selected_item()
	# 	if sel not in self.counted_item:
	# 		self.count_objs()
	# 	convex = self.calc_convex_hull(img=self.target.annotate)
	# 	self.target.img_view.set_data(self.target.data, from_numpy(convex))
	# 	del convex



	def calc_stat(self):
		apix = self.apix_vs.value
		sel = self.get_selected_item()
		#if sel not in self.counted_item:
		#val = int(sel.text(0))
		#raw_mask = self.target.get_segtab().get_whole_annotate(sel)

		#if len(self.target.get_segtab().get_whole_branch(sel))==1:
		self.label_objs(ask_user=False)
		if (self.stat_combo.currentText() == "Center of Mass"):
			self.cent_mass = ndi.center_of_mass(self.target.get_annotation().numpy(),self.labeled_ann,[i+1 for i in range(len(self.objs))])
			self.cent_mass_em = []
			for pair in self.cent_mass:
				if len(pair) == 3:
					self.cent_mass_em.append([pair[2],pair[1],pair[0]])
				elif len(pair) == 2:
					self.cent_mass_em.append([pair[1],pair[0]])
			for  i in range(1,len(self.cent_mass_em)+1):
				print("Center of Mass of Obj",i,":", str(self.cent_mass_em[-i]))
		elif (self.stat_combo.currentText() == "Volume/Area"):
			for  i in range(1,len(self.area_vol)+1):
				print("Volume/Area of Obj",i,":", str(self.area_vol[-i]),"pixels/voxels")
				#print("Largest Perimeter of Obj",i,":", peri_temp*apix,"A")
		elif (self.stat_combo.currentText() == "Largest Area"):
			if  len(self.objs[0].shape)== 2:
				for  i in range(1,len(self.area_vol)+1):
					print("Largest Area of Obj",i,":", str(self.area_vol[-i]),"pixels")
					print("Largest Area of Obj",i,":", str(self.area_vol[-i]*apix*apix)," A^2")
			else:
				#self.projs = [from_numpy(obj).process("xform",{"transform":Transform()}).process("misc.directional_sum",{"axis":"z"}) for obj in self.objs]
				for i in range(1,len(self.objs)+1):
					obj = self.objs[-i]
					a = from_numpy(obj)
					proj = a.process("misc.directional_sum",{"axis":"z"}).numpy()
					area_temp=len(np.where(proj>0)[0])
					print("Largest Area of Obj",i,":", area_temp,"pixel")
					print("Largest Area of Obj",i,":", area_temp*apix*apix,"A^2")

		elif self.stat_combo.currentText() == "Largest Perimeter":

			for i in range(1,len(self.objs)+1):
				obj = self.objs[-i]
				if len(obj.shape) == 2:
					im = obj
				else:
					a = from_numpy(obj)
					proj = a.process("misc.directional_sum",{"axis":"z"}).process("threshold.binary",{"value":0.1})
					im = proj.numpy()
				# sx = ndi.sobel(im, axis=0, mode='constant')
				# sy = ndi.sobel(im, axis=1, mode='constant')
				# sob = np.hypot(sx, sy)
				# from_numpy(sob).write_image('peri_temp_2.hdf')
				peri_temp = self.calc_perimeter(im)
				print("Largest Perimeter of Obj",i,":", peri_temp, "pixel")
				print("Largest Perimeter of Obj",i,":", peri_temp*apix,"A")
				#[data.process("xform",{"transform":xform}).process("misc.directional_sum",{"first":nz//2-layers,"last":nz//2+layers,"axis":"z"}) for data in self.datalist]
			#self.projs = [data.process("xform",{"transform":Transform()}).process("misc.directional_sum",{"first":nz//2-layers,"last":nz//2+layers,"axis":"z"}) for data in self.datalist]
		else:
			return

	def calc_perimeter(self,image):
		(w, h) = image.shape
		data = np.zeros((w + 2, h + 2), dtype=image.dtype)
		data[1:-1, 1:-1] = image
		newdata = np.copy(data)
		for i in range(1, w + 1):
			for j in range(1, h + 1):
				cond = data[i, j] == data[i, j + 1] and \
						data[i, j] == data[i, j - 1] and \
						data[i, j] == data[i + 1, j] and \
						data[i, j] == data[i - 1, j]
				if cond:
					newdata[i, j] = 0
		return np.count_nonzero(newdata)




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
		self.select_anchors_bt = QtWidgets.QPushButton("Select anchors")
		self.detect_plane_bt = QtWidgets.QPushButton("Detect plane")

		self.atlay.addWidget(QtWidgets.QLabel("Ctrl+Click to select 2 points along the filament to detect plane"),0,0,1,2)


		self.atlay.addWidget(self.select_anchors_bt,1,0,1,1)
		self.atlay.addWidget(self.detect_plane_bt,1,1,1,1)


		self.detect_plane_bt.clicked[bool].connect(self.detect_plane)
		self.select_anchors_bt.clicked[bool].connect(self.select_anchors_bt_clicked)


		self.gra_tab = QtWidgets.QWidget()
		self.gtlay = QtWidgets.QVBoxLayout(self.gra_tab)
		self.dark_gra_button = QtWidgets.QPushButton("Dark granule")
		self.alpha_gra_button = QtWidgets.QPushButton("Alpha granule")
		self.ves_button = QtWidgets.QPushButton("Vesicle")
		self.gtlay.addWidget(self.dark_gra_button)
		self.gtlay.addWidget(self.alpha_gra_button)
		self.gtlay.addWidget(self.ves_button)

		self.cell_tab = QtWidgets.QTabWidget()
		self.cell_tab.addTab(self.actin_tab,"Actin/Microtubules")
		self.cell_tab.addTab(self.memb_tab,"Membrane")
		self.cell_tab.addTab(self.gra_tab,"Granule")
		cell_button_l.addWidget(self.cell_tab)
		cell_vbl.addLayout(cell_button_l)

	def detect_plane(self):
		if len(self.target.p1p2) <2:
			print("Must select 2 anchors point to detect plane")
			return
		az = self.target.get_inspector().azs.value/180*np.pi
		p2_p1 = np.subtract(np.array(self.target.p1p2[0]), np.array(self.target.p1p2[1]))
		rot_p2_p1 =  [p2_p1[0]*np.cos(-az) - p2_p1[1]*np.sin(-az),p2_p1[0]*np.sin(-az) + p2_p1[1]*np.cos(-az),p2_p1[2]]
		n = np.cross(rot_p2_p1,np.array([np.cos(az),np.sin(az),0]))
		nn = n / np.linalg.norm(n)
		angles = np.abs(np.arcsin(nn))
		for a in angles:
			print(a/np.pi*180)

		alt = angles[1]/np.pi*180
		#print("ALT", alt)
		self.target.get_inspector().alts.setValue(alt)
		self.target.get_inspector().ns.setValue(0)
		self.target.img_view.del_shape("p1")
		self.target.img_view.del_shape("p2")
		self.target.img_view.updateGL()
		return

	def select_anchors_bt_clicked(self):
		self.target.basic_tab.setCurrentIndex(3)
		print("Select anchors point for detecting plane")
		return

if __name__ == '__main__':
	main()
