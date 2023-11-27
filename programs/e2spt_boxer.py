#!/usr/bin/env python

# LAST update: May/2017 by Muyuan Chen
# Author: Steven Ludtke  2/8/2011 (rewritten)
# Author: Jesus Galaz-Montoya, all command line functionality + updates/enhancements/fixes.
# Author: John Flanagan  9/7/2011 (helixboxer)
# Copyright (c) 2011- Baylor College of Medicine
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




from past.utils import old_div
from builtins import range
from EMAN2 import *
from EMAN2_utils import numpy2pdb
import numpy as np

import weakref
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QSplitter, QHBoxLayout # Erik add for Qsplitter
from PyQt5.QtCore import Qt
from eman2_gui.emapplication import get_application, EMApp
from eman2_gui.emimage2d import EMImage2DWidget
from eman2_gui.emimagemx import EMImageMXWidget
#from emimage3d import EMImage3DWidget
from eman2_gui.emscene3d import EMScene3D
from eman2_gui.emdataitem3d import EMDataItem3D, EMIsosurface
from eman2_gui.emshape import EMShape
from eman2_gui.valslider import ValSlider, ValBox
from sklearn.decomposition import PCA

	
def run(cmd):
	print(cmd)
	launch_childprocess(cmd)
	
def main():
	
	usage="""e2spt_boxer now supports multiple simultaneous features
	
	[prog] <tomogram>
	
	"""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	#parser.add_argument("--path", type=str,help="path", default=None)

	parser.add_pos_argument(name="tomogram",help="Specify a tomogram from which you want to extract particles.", default="", guitype='filebox', browser="EMTomoBoxesTable(withmodal=True,multiselect=False)", row=0, col=0,rowspan=1, colspan=2, mode="box3d,box2d")
	parser.add_argument("--box2d",action="store_true",help="Boxing 2D particls from tomograms.",default=False, guitype='boolbox', row=4, col=0, rowspan=1, colspan=1, mode='box2d[True]')
	parser.add_argument("--box3d",action="store_true",help="Boxing 3D particls from tomograms (default).",default=False, guitype='boolbox', row=4, col=0, rowspan=1, colspan=1, mode='box3d[True]')
	parser.add_header(name="instruction0", help='instruction', title="### Use '~' and '1' to go through slices along Z axis. ###", row=10, col=0, rowspan=1, colspan=2, mode="box3d,box2d")
	parser.add_header(name="instruction1", help='instruction', title="### Hold Shift to delete particles ###", row=11, col=0, rowspan=1, colspan=2, mode="box3d,box2d")
	parser.add_argument("--label", type=str,help="start from viewing particles of the specified label.", default=None)


	#parser.add_argument("--mode", type=str,help="Boxing mode. choose from '2D' and '3D'. Default is '3D'", default="3D",guitype='combobox',choicelist="('2D', '3D')",row=1, col=0, rowspan=1, colspan=1,mode="boxing")

	parser.add_argument("--ppid", type=int,help="ppid", default=-2)

	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)

	if len(args) == 0:
		print("INPUT ERROR: You must specify a tomogram.")
		sys.exit(1)
	
	if options.box2d==True and options.box3d==False:
		options.mode="2D"
	else:
		options.mode="3D"
	
	img = args[0]

	app = EMApp()

	img=args[0]

	imghdr = EMData(img,0,True)
	options.apix = imghdr['apix_x']
	
	boxer=EMTomoBoxer(app,options,datafile=img)

	boxer.show()
	app.execute()
	E2end(logid)
	return

class EMTomoBoxer(QtWidgets.QMainWindow):
	"""This class represents the EMTomoBoxer application instance.  """
	keypress = QtCore.pyqtSignal(QtGui.QKeyEvent)
	module_closed = QtCore.pyqtSignal()

	def __init__(self,application,options,datafile):
		QtWidgets.QWidget.__init__(self)
		self.initialized=False
		self.app=weakref.ref(application)
		self.options=options
		self.apix=options.apix
		self.currentset=0
		self.shrink=1#options.shrink
		self.setWindowTitle("Main Window (e2spt_boxer.py)")
		if options.mode=="3D":
			self.boxshape="circle"
		else:
			self.boxshape="rect"

		self.globalxf=Transform()
		
		# Menu Bar
		self.mfile=self.menuBar().addMenu("File")
		#self.mfile_open=self.mfile.addAction("Open")
		self.mfile_read_boxloc=self.mfile.addAction("Read Box Coord")
		self.mfile_save_boxloc=self.mfile.addAction("Save Box Coord")
		self.mfile_save_boxpdb=self.mfile.addAction("Save Coord as PDB")
		self.mfile_save_boxes_stack=self.mfile.addAction("Save Boxes as Stack")
		self.mfile_save_gif=self.mfile.addAction("Save GIF movie")
		#self.mfile_quit=self.mfile.addAction("Quit")


#old code, Erik commented out
#		self.setCentralWidget(QtWidgets.QWidget())
#		self.gbl = QtWidgets.QGridLayout(self.centralWidget())

#New QSplitter code, Erik addition
		self.setCentralWidget(QtWidgets.QWidget())
		self.gbl = QtWidgets.QGridLayout(self.centralWidget())
		self.gblh = QtWidgets.QHBoxLayout()
		self.gbl.addLayout(self.gblh,0,0,2,2)

		self.splitter_top = QSplitter(Qt.Horizontal) #top panel of images
		self.splitter_bottom = QSplitter(Qt.Horizontal) #bottom panel of images
		self.splitter_wrapper = QSplitter(Qt.Vertical) #will stack the top and bottom image panels

		self.xyview = EMImage2DWidget(sizehint=(1024,1024))
		self.xzview = EMImage2DWidget(sizehint=(1024,256))
		self.zyview	= EMImage2DWidget(sizehint=(256,1024))
		
		self.wdepth = QtWidgets.QSlider()
		self.gbl.addWidget(self.wdepth,0,2)

		self.splitter_top.addWidget(self.zyview)
		self.splitter_top.addWidget(self.xyview)
		
		self.splitter_wrapper.addWidget(self.splitter_top)
		self.splitter_wrapper.addWidget(self.splitter_bottom)
		
		self.gblh.addWidget(self.splitter_wrapper)
		
		#control panel
		self.gbl2 = QtWidgets.QGridLayout()
		self.grid_widget = QtWidgets.QWidget()
		self.grid_widget.setLayout(self.gbl2)
		self.splitter_bottom.addWidget(self.grid_widget)

		self.splitter_bottom.addWidget(self.xzview)
		self.splitter_top.splitterMoved.connect(self.splitter_bottom.moveSplitter)


#########################################################################
		
		# relative stretch factors
		#self.gbl.setColumnMinimumWidth(0,200)
		#self.gbl.setRowMinimumHeight(0,200)
		#self.gbl.setColumnStretch(0,0)
		
		self.wzheight=ValBox(label="Z height:",value=256)
		self.gbl2.addWidget(self.wzheight,1,0)

		# box size
		self.wboxsize=ValBox(label="Box Size:",value=0)
		self.gbl2.addWidget(self.wboxsize,2,0)
		

		# max or mean
		#self.wmaxmean=QtWidgets.QPushButton("MaxProj")
		#self.wmaxmean.setCheckable(True)
		#self.gbl2.addWidget(self.wmaxmean,3,0)

		# number slices
		label0=QtWidgets.QLabel("Thickness")
		self.gbl2.addWidget(label0,3,0)

		self.wnlayers=QtWidgets.QSpinBox()
		self.wnlayers.setMinimum(1)
		self.wnlayers.setMaximum(256)
		self.wnlayers.setValue(1)
		self.gbl2.addWidget(self.wnlayers,3,1)

		# Local boxes in side view
		self.wlocalbox=QtWidgets.QCheckBox("Limit Side Boxes")
		self.gbl2.addWidget(self.wlocalbox,4,0)
		self.wlocalbox.setChecked(True)
		
		self.button_flat = QtWidgets.QPushButton("Flatten")
		self.gbl2.addWidget(self.button_flat,5,0)
		self.button_reset = QtWidgets.QPushButton("Reset")
		self.gbl2.addWidget(self.button_reset,5,1)
		## scale factor
		#self.wscale=ValSlider(rng=(.1,2),label="Sca:",value=1.0)
		#self.gbl2.addWidget(self.wscale,4,0,1,2)

		# 2-D filters
		self.wfilt = ValSlider(rng=(0,150),label="Filt",value=0.0)
		self.gbl2.addWidget(self.wfilt,6,0,1,2)
		
		self.curbox=-1
		
		self.boxes=[]						# array of box info, each is (x,y,z,...)
		self.boxesimgs=[]					# z projection of each box
		self.dragging=-1

		##coordinate display
		self.wcoords=QtWidgets.QLabel("")
		self.gbl2.addWidget(self.wcoords, 0, 0, 1, 2)
		
		self.button_flat.clicked[bool].connect(self.flatten_tomo)
		self.button_reset.clicked[bool].connect(self.reset_flatten_tomo)

		# file menu
		#self.mfile_open.triggered[bool].connect(self.menu_file_open)
		self.mfile_read_boxloc.triggered[bool].connect(self.menu_file_read_boxloc)
		self.mfile_save_boxloc.triggered[bool].connect(self.menu_file_save_boxloc)
		self.mfile_save_boxpdb.triggered[bool].connect(self.menu_file_save_boxpdb)
		
		self.mfile_save_boxes_stack.triggered[bool].connect(self.save_boxes)
		self.mfile_save_gif.triggered[bool].connect(self.save_gif)
		#self.mfile_quit.triggered[bool].connect(self.menu_file_quit)

		# all other widgets
		self.wdepth.valueChanged[int].connect(self.event_depth)
		self.wnlayers.valueChanged[int].connect(self.event_nlayers)
		self.wboxsize.valueChanged.connect(self.event_boxsize)
		
#Erik commented out because QHBoxlayerout has no attribute 'setRowMinimumheight'
#		self.wzheight.valueChanged.connect(self.event_zheight)

		#self.wmaxmean.clicked[bool].connect(self.event_projmode)
		#self.wscale.valueChanged.connect(self.event_scale)
		self.wfilt.valueChanged.connect(self.event_filter)
		self.wlocalbox.stateChanged[int].connect(self.event_localbox)

		self.xyview.mousemove.connect(self.xy_move)
		self.xyview.mousedown.connect(self.xy_down)
		self.xyview.mousedrag.connect(self.xy_drag)
		self.xyview.mouseup.connect(self.mouse_up)
		self.xyview.mousewheel.connect(self.xy_wheel)
		self.xyview.signal_set_scale.connect(self.event_scale)
		self.xyview.origin_update.connect(self.xy_origin)

		self.xzview.mousedown.connect(self.xz_down)
		self.xzview.mousedrag.connect(self.xz_drag)
		self.xzview.mouseup.connect(self.mouse_up)
		self.xzview.mousewheel.connect(self.xz_wheel)
		self.xzview.signal_set_scale.connect(self.event_scale)
		self.xzview.origin_update.connect(self.xz_origin)
		self.xzview.mousemove.connect(self.xz_move)

		self.zyview.mousedown.connect(self.zy_down)
		self.zyview.mousedrag.connect(self.zy_drag)
		self.zyview.mouseup.connect(self.mouse_up)
		self.zyview.mousewheel.connect(self.zy_wheel)
		self.zyview.signal_set_scale.connect(self.event_scale)
		self.zyview.origin_update.connect(self.zy_origin)
		self.zyview.mousemove.connect(self.zy_move)
		
		self.xyview.keypress.connect(self.key_press)
		self.datafilename=datafile
		self.basename=base_name(datafile)
		p0=datafile.find('__')
		if p0>0:
			p1=datafile.rfind('.')
			self.filetag=datafile[p0:p1]
			if self.filetag[-1]!='_':
				self.filetag+='_'
		else:
			self.filetag="__"
			
		data=EMData(datafile)
		self.set_data(data)

		# Boxviewer subwidget (details of a single box)
		#self.boxviewer=EMBoxViewer()
		#self.app().attach_child(self.boxviewer)

		# Boxes Viewer (z projections of all boxes)
		self.boxesviewer=EMImageMXWidget()
		
		#self.app().attach_child(self.boxesviewer)
		self.boxesviewer.show()
		self.boxesviewer.set_mouse_mode("App")
		self.boxesviewer.setWindowTitle("Particle List")
		self.boxesviewer.rzonce=True
		
		self.setspanel=EMTomoSetsPanel(self)

		self.optionviewer=EMTomoBoxerOptions(self)
		self.optionviewer.add_panel(self.setspanel,"Sets")
		
		
		self.optionviewer.show()
		
		self.boxesviewer.mx_image_selected.connect(self.img_selected)
		
		##################
		#### deal with metadata in the _info.json file...
		
		self.jsonfile=info_name(datafile)
		info=js_open_dict(self.jsonfile)
		
		#### read particle classes
		self.sets={}
		self.boxsize={}
		if "class_list" in info:
			clslst=info["class_list"]
			for k in sorted(clslst.keys()):
				if type(clslst[k])==dict:
					self.sets[int(k)]=str(clslst[k]["name"])
					self.boxsize[int(k)]=int(clslst[k]["boxsize"])
				else:
					self.sets[int(k)]=str(clslst[k])
					self.boxsize[int(k)]=64
					
		clr=QtGui.QColor
		self.setcolors=[QtGui.QBrush(clr("blue")),QtGui.QBrush(clr("green")),QtGui.QBrush(clr("red")),QtGui.QBrush(clr("cyan")),QtGui.QBrush(clr("purple")),QtGui.QBrush(clr("orange")), QtGui.QBrush(clr("yellow")),QtGui.QBrush(clr("hotpink")),QtGui.QBrush(clr("gold"))]
		self.sets_visible={}
				
		#### read boxes
		if "boxes_3d" in info:
			box=info["boxes_3d"]
			for i,b in enumerate(box):
				#### X-center,Y-center,Z-center,method,[score,[class #]]
				bdf=[0,0,0,"manual",0.0, 0, 0]
				for j,bi in enumerate(b):  bdf[j]=bi
				
				
				if bdf[5] not in list(self.sets.keys()):
					clsi=int(bdf[5])
					self.sets[clsi]="particles_{:02d}".format(clsi)
					self.boxsize[clsi]=64
				
				self.boxes.append(bdf)
		
		###### this is the new (2018-09) metadata standard..
		### now we use coordinates at full size from center of tomogram so it works for different binning and clipping
		### have to make it compatible with older versions though..
		if "apix_unbin" in info:
			self.apix_unbin=info["apix_unbin"]
			self.apix_cur=apix=data["apix_x"]
			for b in self.boxes:
				b[0]=b[0]/apix*self.apix_unbin+data["nx"]//2
				b[1]=b[1]/apix*self.apix_unbin+data["ny"]//2
				b[2]=b[2]/apix*self.apix_unbin+data["nz"]//2
				
			for k in self.boxsize.keys():
				self.boxsize[k]=int(np.round(self.boxsize[k]*self.apix_unbin/apix))
		else:
			self.apix_unbin=-1
			
		info.close()
		
		E2loadappwin("e2sptboxer","main",self)
		E2loadappwin("e2sptboxer","boxes",self.boxesviewer.qt_parent)
		E2loadappwin("e2sptboxer","option",self.optionviewer)
		
		#### particle classes
		if len(self.sets)==0:
			self.new_set("particles_00")
			
		self.currentset=sorted(self.sets.keys())[0]
		if options.label:
			sid=[i for i in self.sets.keys() if self.sets[i]==options.label]
			if len(sid)>0:
				self.currentset=sid[0]
			else:
				print("cannot find specified label")
				
		self.sets_visible[self.currentset]=0
		self.setspanel.update_sets()
		self.wboxsize.setValue(self.get_boxsize())

		#print(self.sets)
		for i in range(len(self.boxes)):
			self.update_box(i)
		
		self.update_all()
		self.initialized=True
#		self.splitter_bottom.moveSplitter(self.splitter_top.handle(0).pos(),0)
		
	def set_data(self,data):

		self.data=data
		self.apix=data["apix_x"]

		self.datasize=(data["nx"],data["ny"],data["nz"])
		self.x_loc, self.y_loc, self.z_loc=data["nx"]//2,data["ny"]//2,data["nz"]//2

		#self.gbl.setRowMinimumHeight(1,max(250,data["nz"]))
		#self.gbl.setColumnMinimumWidth(0,max(250,data["nz"]))
		#print(data["nx"],data["ny"],data["nz"])
		self.wzheight.setValue(data["nz"])

		self.wdepth.setRange(0,data["nz"]-1)
		self.wdepth.setValue(data["nz"]//2)
		self.boxes=[]
		self.curbox=-1

		if self.initialized:
			self.update_all()

	def eraser_width(self):
		return int(self.optionviewer.eraser_radius.getValue())
		
	def get_cube(self,x,y,z, centerslice=False, boxsz=-1):
		"""Returns a box-sized cube at the given center location"""
		if boxsz<0:
			bs=self.get_boxsize()
		else:
			bs=boxsz
			
		if centerslice:
			bz=1
		else:
			bz=bs
		
		if ((x<-bs//2) or (y<-bs//2) or (z<-bz//2)
			or (x>self.data["nx"]+bs//2) or (y>self.data["ny"]+bs//2) or (z>self.data["nz"]+bz//2) ):
			r=EMData(bs,bs,bz)
		else:
			r=self.data.get_clip(Region(x-bs//2,y-bs//2,z-bz//2,bs,bs,bz))

		if self.apix!=0 :
			r["apix_x"]=r["apix_y"]=r["apix_z"]=self.apix

		return r

	def get_slice(self,idx,thk=1,axis="z"):
		if not self.globalxf.is_identity():
			data=self.dataxf
		elif self.wfilt.getValue()!=0.0:
			data=self.datalp
		else:
			data=self.data
			
		t=int(thk-1)
		idx=int(idx)
		r=data.process("misc.directional_sum",{"axis":axis,"first":idx-t,"last":idx+t})
		r.div(t*2+1)
		
		if self.apix!=0 :
			r["apix_x"]=r["apix_y"]=r["apix_z"]=self.apix
		return r

#Erik commented out because QHBoxLayout has no attribute 'setRowMinimumheight'
#	def event_zheight(self):
#		z=self.wzheight.getValue()
#		self.gbl.setRowMinimumHeight(1,z)
#		self.gbl.setColumnMinimumWidth(0,z)
#		#print(data["nx"],data["ny"],data["nz"])
#		return
	def save_gif(self):
		

		cmd, ok = QtWidgets.QInputDialog.getText(self, 'Command', 'xy:100-200')
		if not ok:
			return
		
		try:
			from PIL import ImageGrab
		except:
			print("require PIL")
			return
		
			
		cmd=str(cmd)
		print(cmd)
		cmd=cmd.split(':')
		win=cmd[0]
		if win=="xy":
			view=self.xyview
		elif win=="xz":
			view=self.xzview
		elif win=="zy":
			view=self.zyview
		else:
			print("choose from xy, xz, zy.")
			return
		
		z=[int(i) for i in cmd[1].split('-')]
		fnames=[]
		for i in range(z[0], z[1]):
			if win=="xy":
				self.wdepth.setValue(i)
			elif win=="xz":
				self.y_loc=i
				self.update_sliceview(['y'])
			else:
				self.x_loc=i
				self.update_sliceview(['x'])
				
			p= view.mapToGlobal(QtCore.QPoint(0, 0))
			ss_region=(p.x(), p.y(), p.x()+view.width(), p.y()+view.height())
			# print(ss_region)
			ss_img = ImageGrab.grab(ss_region)
			fnames.append(f"snapshot_{i:03d}.png")
			ss_img.save(fnames[-1])
# 		
		c="convert {} snapshot.gif".format(' '.join(fnames+fnames[::-1]))
		print(c)
		os.system(c)

	
	def event_boxsize(self):
		if self.get_boxsize()==int(self.wboxsize.getValue()):
			return
		
		self.boxsize[self.currentset]=int(self.wboxsize.getValue())
		
		#cb=self.curbox
		self.initialized=False
		for i in range(len(self.boxes)):
			if self.boxes[i][5]==self.currentset:
				self.update_box(i)
		#self.update_box(cb)
		self.initialized=True
		self.update_all()

	#def event_projmode(self,state):
		#"""Projection mode can be simple average (state=False) or maximum projection (state=True)"""
		#self.update_all()

	def event_scale(self,newscale):
		self.xyview.set_scale(newscale)
		self.xzview.set_scale(newscale)
		self.zyview.set_scale(newscale)

	def event_depth(self):
		if self.z_loc!=self.wdepth.value():
			self.z_loc=self.wdepth.value()
		if self.initialized:
			self.update_sliceview()

	def event_nlayers(self):
		self.update_all()

	def event_filter(self):
		if self.wfilt.getValue()!=0.0:
			print("Filtering tomogram...")
			self.datalp=self.data.process("filter.lowpass.gauss",{"cutoff_freq":1.0/self.wfilt.getValue()})
		else:
			self.datalp=self.data
			
		self.update_all()

	def event_localbox(self,tog):
		self.update_sliceview()

	def get_boxsize(self, clsid=-1):
		if clsid<0:
			return int(self.boxsize[self.currentset])
		else:
			try:
				ret= int(self.boxsize[clsid])
			except:
				print("No box size saved for {}..".format(clsid))
				ret=32
			return ret

	def nlayers(self):
		return int(self.wnlayers.value())

	def menu_file_read_boxloc(self):
		fsp=str(QtWidgets.QFileDialog.getOpenFileName(self, "Select output text file")[0])
		
		if not os.path.isfile(fsp):
			print("file does not exist")
			return

		f=open(fsp,"r")
		for b in f:
			b2=[old_div(int(float(i)),self.shrink) for i in b.split()[:3]]
			bdf=[0,0,0,"manual",0.0, self.currentset]
			for j in range(len(b2)):
				bdf[j]=b2[j]
			self.boxes.append(bdf)
			self.update_box(len(self.boxes)-1)
		f.close()

	def menu_file_save_boxloc(self):
		shrinkf=self.shrink 								#jesus

		fsp=str(QtWidgets.QFileDialog.getSaveFileName(self, "Select output text file")[0])
		if len(fsp)==0:
			return
		
		clsid=list(self.sets_visible.keys())
		if len(clsid)==0:
			print("No visible particles to save")
			return
		
		out=open(fsp,"w")
		for b in self.boxes:
			if int(b[5]) in clsid:
				out.write("%d\t%d\t%d\n"%(b[0]*shrinkf,b[1]*shrinkf,b[2]*shrinkf))
		out.close()
		
	def menu_file_save_boxpdb(self):
		fsp=str(QtWidgets.QFileDialog.getSaveFileName(self, "Select output PDB file", filter="PDB (*.pdb)")[0])
		if len(fsp)==0:
			return
		if fsp[-4:].lower()!=".pdb" :
			fsp+=".pdb"
		clsid=list(self.sets_visible.keys())
		if len(clsid)==0:
			print("No visible particles to save")
			return
		
		bxs=np.array([[b[0], b[1], b[2]] for b in self.boxes if int(b[5]) in clsid])/10
		
		numpy2pdb(bxs, fsp)
		print("PDB saved to {}. Use voxel size 0.1".format(fsp))
		
	def save_boxes(self, clsid=[]):
		if len(clsid)==0:
			defaultname="ptcls.hdf"
		else:
			defaultname="_".join([self.sets[i] for i in clsid])+".hdf"
		
		name,ok=QtWidgets.QInputDialog.getText( self, "Save particles", "Filename suffix:", text=defaultname)
		if not ok:
			return
		name=self.filetag+str(name)
		if name[-4:].lower()!=".hdf" :
			name+=".hdf"
			
			
		if self.options.mode=="3D":
			dr="particles3d"
			is2d=False
		else:
			dr="particles"
			is2d=True
		
		
		if not os.path.isdir(dr):
			os.mkdir(dr)
		
		fsp=os.path.join(dr,self.basename)+name

		print("Saving {} particles to {}".format(self.options.mode, fsp))
		
		if os.path.isfile(fsp):
			print("{} exist. Overwritting...".format(fsp))
			os.remove(fsp)
		
		progress = QtWidgets.QProgressDialog("Saving", "Abort", 0, len(self.boxes),None)
		
		
		boxsz=-1
		for i,b in enumerate(self.boxes):
			if len(clsid)>0:
				if int(b[5]) not in clsid:
					continue
			
			#img=self.get_cube(b[0],b[1],b[2])
			bs=self.get_boxsize(b[5])
			if boxsz<0:
				boxsz=bs
			else:
				if boxsz!=bs:
					print("Inconsistant box size in the particles to save.. Using {:d}..".format(boxsz))
					bs=boxsz
			
			sz=[s//2 for s in self.datasize]
			
			img=self.get_cube(b[0], b[1], b[2], centerslice=is2d, boxsz=bs)
			if is2d==False:
				img.process_inplace('normalize')
			
			img["ptcl_source_image"]=self.datafilename
			img["ptcl_source_coord"]=(b[0]-sz[0], b[1]-sz[1], b[2]-sz[2])
			
			if is2d==False: #### do not invert contrast for 2D images
				img.mult(-1)
			
			img.write_image(fsp,-1)

			progress.setValue(i+1)
			if progress.wasCanceled():
				break

	def update_sliceview(self, axis=['x','y','z']):
		boxes=self.get_rotated_boxes()
		
		allside=(not self.wlocalbox.isChecked())
		
		pms={'z':[2, self.xyview, self.z_loc],
		     'y':[1, self.xzview, self.y_loc],
		     'x':[0, self.zyview, self.x_loc]} 
		
		if self.boxshape=="circle": lwi=7
		else: lwi=8
		
		for ax in axis:
			ia, view, loc=pms[ax]
			shp=view.get_shapes()
			if len(shp)!=len(boxes):
				### something changes the box shapes...
				for i,b in enumerate(boxes):
					self.update_box_shape(i,b)
			
		for ax in axis:
			ia, view, loc=pms[ax]
			
			## update the box shapes
			shp=view.get_shapes()
			for i,b in enumerate(boxes):
				bs=self.get_boxsize(b[5])
				dst=abs(b[ia] - loc)
				
				inplane=dst<bs//2
				rad=bs//2-dst
				
				if allside:
					## display all side view boxes in this mode
					inplane=True
					rad=bs//2
					
				if ax=='z' and self.options.mode=="2D":
					## boxes are 1 slice thick in 2d mode
					inplane=dst<1
				
				
				if inplane and (b[5] in self.sets_visible):
					shp[i][0]=self.boxshape
					## selected box is slightly thicker
					if self.curbox==i:
						shp[i][lwi]=3
					else:
						shp[i][lwi]=2
					if self.options.mode=="3D":
						shp[i][6]=rad
				else:
					shp[i][0]="hidden"
				
			view.shapechange=1
			img=self.get_slice(loc, self.nlayers(), ax)
			#if self.wfilt.getValue()!=0.0:
				#img.process_inplace("filter.lowpass.gauss",{"cutoff_freq":1.0/self.wfilt.getValue(),"apix":self.apix})

			view.set_data(img)
			
		self.update_coords()

	def get_rotated_boxes(self):
		if len(self.boxes)==0:
			return []
		if self.globalxf.is_identity():
			boxes=self.boxes
		else:
			cnt=[self.data["nx"]//2,self.data["ny"]//2,self.data["nz"]//2]
			pts=np.array([b[:3] for b in self.boxes])-cnt
			pts=np.array([self.globalxf.transform(p.tolist()) for p in pts])+cnt
			boxes=[]
			for i,b in enumerate(self.boxes):
				p=pts[i]
				boxes.append([p[0],p[1],p[2],b[3],b[4],b[5]])
		return boxes
		
	def update_all(self):
		"""redisplay of all widgets"""
		if self.data==None:
			return

		self.update_sliceview()
		self.update_boximgs()

	def update_coords(self):
		self.wcoords.setText("X: {:d}\tY: {:d}\tZ: {:d}".format(int(self.x_loc), int(self.y_loc), int(self.z_loc)))

	def inside_box(self,n,x=-1,y=-1,z=-1):
		"""Checks to see if a point in image coordinates is inside box number n. If any value is negative, it will not be checked."""
		box=self.boxes[n]
		if box[5] not in self.sets_visible:
			return False
		bs=self.get_boxsize(box[5])/2
		if self.options.mode=="3D":
			rr=(x>=0)*((box[0]-x)**2) + (y>=0)*((box[1]-y) **2) + (z>=0)*((box[2]-z)**2)
		else:
			rr=(x>=0)*((box[0]-x)**2) + (y>=0)*((box[1]-y) **2) + (z>=0)*(box[2]!=z)*(1e3*bs**2)
		return rr<=bs**2

	def del_box(self, delids):
		
		if type(delids)!=list:
			delids=[delids]
		
		kpids=[i for i,b in enumerate(self.boxes) if i not in delids]
		self.boxes=[self.boxes[i] for i in kpids]
		self.boxesimgs=[self.boxesimgs[i] for i in kpids]
		self.xyview.shapes={i:self.xyview.shapes[k] for i,k in enumerate(kpids)}
		self.xzview.shapes={i:self.xzview.shapes[k] for i,k in enumerate(kpids)}
		self.zyview.shapes={i:self.zyview.shapes[k] for i,k in enumerate(kpids)}
		#print self.boxes, self.xyview.get_shapes()
		self.curbox=-1
		self.update_all()
	
	def update_box_shape(self,n, box):
		bs2=self.get_boxsize(box[5])//2
		if n==self.curbox:
			lw=3
		else:
			lw=2
		color=self.setcolors[box[5]%len(self.setcolors)].color().getRgbF()
		if self.options.mode=="3D":
			self.xyview.add_shape(n,EMShape(["circle",color[0],color[1],color[2],box[0],box[1],bs2,lw]))
			self.xzview.add_shape(n,EMShape(["circle",color[0],color[1],color[2],box[0],box[2],bs2,lw]))
			self.zyview.add_shape(n,EMShape(("circle",color[0],color[1],color[2],box[2],box[1],bs2,lw)))
		else:
			self.xyview.add_shape(n,EMShape(["rect",color[0],color[1],color[2],
				    box[0]-bs2,box[1]-bs2,box[0]+bs2,box[1]+bs2,2]))
			self.xzview.add_shape(n,EMShape(["rect",color[0],color[1],color[2], 
				    box[0]-bs2,box[2]-1,box[0]+bs2,box[2]+1,2]))
			self.zyview.add_shape(n,EMShape(["rect",color[0],color[1],color[2],
				    box[2]-1,box[1]-bs2,box[2]+1,box[1]+bs2,2]))
			
	
	def update_box(self,n,quiet=False):
		"""After adjusting a box, call this"""
#		print "upd ",n,quiet
		if n<0 or n>=len(self.boxes):
			return
		
		box=self.boxes[n]
		
		boxes=self.get_rotated_boxes()
		self.update_box_shape(n,boxes[n])

		if self.initialized: 
			self.update_sliceview()

		# For speed, we turn off updates while dragging a box around. Quiet is set until the mouse-up
		if not quiet:
			# Get the cube from the original data (normalized)
			proj=self.get_cube(box[0], box[1], box[2], centerslice=True, boxsz=self.get_boxsize(box[5]))
			proj.process_inplace("normalize")
			
			for i in range(len(self.boxesimgs),n+1): 
				self.boxesimgs.append(None)
			
			self.boxesimgs[n]=proj

			mm=[m for im,m in enumerate(self.boxesimgs) if self.boxes[im][5] in self.sets_visible]
			
			if self.initialized: self.SaveJson()
			
		if self.initialized:
			self.update_boximgs()
			#if n!=self.curbox:
				#self.boxesviewer.set_selected((n,),True)

		self.curbox=n

	def update_boximgs(self):
		self.boxids=[im for im,m in enumerate(self.boxesimgs) if self.boxes[im][5] in self.sets_visible]
		self.boxesviewer.set_data([self.boxesimgs[i] for i in self.boxids])
		self.boxesviewer.update()
		return

	def img_selected(self,event,lc):
		#print("sel",lc[0])
		lci=self.boxids[lc[0]]
		if event.modifiers()&Qt.ShiftModifier:
			if event.modifiers()&Qt.ControlModifier:
				self.del_box(list(range(lci, len(self.boxes))))
			else:
				self.del_box(lci)
		else:
			#self.update_box(lci)
			self.curbox=lci
			box=self.boxes[lci]
			self.x_loc,self.y_loc,self.z_loc=self.rotate_coord([box[0], box[1], box[2]], inv=False)
			self.scroll_to(self.x_loc,self.y_loc,self.z_loc)
			
			self.update_sliceview()
			
			
	def rotate_coord(self, p, inv=True):
		if not self.globalxf.is_identity():
			cnt=[self.data["nx"]//2,self.data["ny"]//2,self.data["nz"]//2]
			p=[p[i]-cnt[i] for i in range(3)]
			xf=Transform(self.globalxf)
			if inv:
				xf.invert()
			p=xf.transform(p)+cnt
		return p

	def del_region_xy(self, x=-1, y=-1, z=-1, rad=-1):
		if rad<0:
			rad=self.eraser_width()
		
		#print(x,y,z, rad)
		delids=[]
		boxes=self.get_rotated_boxes()
		for i,b in enumerate(boxes):
			if b[5] not in self.sets_visible:
				continue
			
			if (x>=0)*(b[0]-x)**2 + (y>=0)*(b[1]-y)**2 +(z>=0)*(b[2]-z)**2 < rad**2:
				delids.append(i)
		self.del_box(delids)
	
	def scroll_to(self, x,y,z, axis=""):
		if axis!="z": self.xyview.scroll_to(x,y,True)
		if axis!="y": self.xzview.scroll_to(x,self.data["nz"]/2,True)
		if axis!="x": self.zyview.scroll_to(self.data["nz"]/2,y,True)
	
	#### mouse click
	def xy_down(self,event):
		x,y=self.xyview.scr_to_img((event.x(),event.y()))
		self.mouse_down(event, x,y,self.z_loc, "z")
		
	def xz_down(self,event):
		x,z=self.xzview.scr_to_img((event.x(),event.y()))
		self.mouse_down(event,x,self.y_loc,z, "y")
			
	def zy_down(self,event):
		z,y=self.zyview.scr_to_img((event.x(),event.y()))
		self.mouse_down(event,self.x_loc,y,z, "x")
		
	def mouse_down(self,event, x, y, z, axis):
		if min(x,y,z)<0: return
		
		xr,yr,zr=self.rotate_coord((x,y,z))
		#print(x,y,z,xr,yr,zr)
	
		if self.optionviewer.erasercheckbox.isChecked():
			
			side=self.wlocalbox.isChecked()
			xyz={'x':x,'y':y,'z':z}
			if not side:
				xyz[axis]=-1
				
			self.del_region_xy(xyz['x'],xyz['y'],xyz['z'],-1)
			return
			
		for i in range(len(self.boxes)):
			if self.inside_box(i,xr,yr,zr):
				
				if event.modifiers()&Qt.ShiftModifier:  ## delete box
					self.del_box(i)

				else:  ## start dragging
					self.dragging=i
					self.curbox=i
					self.scroll_to(x,y,z,axis)
					
				break
		else:
			if not event.modifiers()&Qt.ShiftModifier: ## add box

				self.x_loc, self.y_loc, self.z_loc=x,y,z
				self.scroll_to(x,y,z,axis)
				self.curbox=len(self.boxes)
				self.boxes.append(([xr,yr,zr, 'manual', 0.0, self.currentset]))
				self.update_box(len(self.boxes)-1)
				self.dragging=len(self.boxes)-1
				
				

	#### eraser mode
	def xy_move(self,event):
		self.mouse_move(event, self.xyview)
			
	def xz_move(self,event):
		self.mouse_move(event, self.xzview)
		
	def zy_move(self,event):
		self.mouse_move(event, self.zyview)
			
	
	def mouse_move(self,event,view):
		
		if self.optionviewer.erasercheckbox.isChecked(): 
			self.xyview.eraser_shape=self.xzview.eraser_shape=self.zyview.eraser_shape=None
			x,y=view.scr_to_img((event.x(),event.y()))
			view.eraser_shape=EMShape(["circle",1,1,1,x,y,self.eraser_width(),2])
			view.shapechange=1
			view.update()
		else:
			view.eraser_shape=None
			
	
	#### dragging...
	def mouse_drag(self,x, y, z):
		if self.dragging<0:
			return
		if min(x,y,z)<0:
			return
		
		self.x_loc, self.y_loc, self.z_loc=x,y,z
		x,y,z=self.rotate_coord((x,y,z))
		self.boxes[self.dragging][:3]= x,y,z
		self.update_box(self.dragging,True)

	def xy_drag(self,event):
		if self.dragging>=0:
			x,y=self.xyview.scr_to_img((event.x(),event.y()))
			self.mouse_drag(x,y,self.z_loc)

	def xz_drag(self,event):
		if self.dragging>=0:
			x,z=self.xzview.scr_to_img((event.x(),event.y()))
			self.mouse_drag(x,self.y_loc,z)
	
	def zy_drag(self,event):
		if self.dragging>=0:
			z,y=self.zyview.scr_to_img((event.x(),event.y()))
			self.mouse_drag(self.x_loc,y,z)
		
	def mouse_up(self,event):
		if self.dragging>=0:
			self.update_box(self.dragging)
		self.dragging=-1

	
	#### keep the same origin for the 3 views
	def xy_origin(self,newor):
		xzo=self.xzview.get_origin()
		self.xzview.set_origin(newor[0],xzo[1],True)

		zyo=self.zyview.get_origin()
		self.zyview.set_origin(zyo[0],newor[1],True)
	
	def xz_origin(self,newor):
		xyo=self.xyview.get_origin()
		self.xyview.set_origin(newor[0],xyo[1],True)

	def zy_origin(self,newor):
		xyo=self.xyview.get_origin()
		self.xyview.set_origin(xyo[0],newor[1],True)


	##### go up/down with shift+wheel
	def xy_wheel(self, event):
		z=int(self.z_loc+ np.sign(event.angleDelta().y()))
		if z>0 and z<self.data["nz"]:
			self.wdepth.setValue(z)
	
	def xz_wheel(self, event):
		y=int(self.y_loc+np.sign(event.angleDelta().y()))
		if y>0 and y<self.data["ny"]:
			self.y_loc=y
			self.update_sliceview(['y'])
		
	def zy_wheel(self, event):
		x=int(self.x_loc+np.sign(event.angleDelta().y()))
		if x>0 and x<self.data["nx"]:
			self.x_loc=x
			self.update_sliceview(['x'])
			

	########
	def set_current_set(self, name):
		
		#print "set current", name
		name=parse_setname(name)
		self.currentset=name
		self.wboxsize.setValue(self.get_boxsize())
		self.update_all()
		return
	
	
	def hide_set(self, name):
		name=parse_setname(name)
		
		if name in self.sets_visible: self.sets_visible.pop(name)
		
		
		if self.initialized: 
			self.update_all()
			self.update_boximgs()
		return
	
	
	def show_set(self, name):
		name=parse_setname(name)
		self.sets_visible[name]=0
		#self.currentset=name
		#self.wboxsize.setValue(self.get_boxsize())
		if self.initialized: 
			self.update_all()
			self.update_boximgs()
		return
	
	
	def delete_set(self, name):
		name=parse_setname(name)
		## idx to keep
		delids=[i for i,b in enumerate(self.boxes) if b[5]==int(name)]
		self.del_box(delids)
		
		if name in self.sets_visible: self.sets_visible.pop(name)
		if name in self.sets: self.sets.pop(name)
		if name in self.boxsize: self.boxsize.pop(name)
		
		self.update_all()
		
		return
	
	def rename_set(self, oldname,  newname):
		name=parse_setname(oldname)
		if name in self.sets: 
			self.sets[name]=newname
		return
	
	
	def new_set(self, name):
		for i in range(len(self.sets)+1):
			if i not in self.sets:
				break
			
		self.sets[i]=name
		self.sets_visible[i]=0
		if self.options.mode=="3D":
			self.boxsize[i]=32
		else:
			self.boxsize[i]=64
		
		return
	
	def save_set(self):
		
		self.save_boxes(list(self.sets_visible.keys()))
		return
	
	
	def key_press(self,event):
		if event.key() == 96: ## "`" to move up a slice since arrow keys are occupied...
			self.wdepth.setValue(self.z_loc+1)

		elif event.key() == 49: ## "1" to move down a slice
			self.wdepth.setValue(self.z_loc-1)
		else:
			self.keypress.emit(event)

	def flatten_tomo(self):
		print("Flatten tomogram by particles coordinates")
		vis=list(self.sets_visible.keys())
		pts=[b[:3] for b in self.boxes if b[5] in vis]
		if len(pts)<3:
			print("Too few visible particles. Cannot flatten tomogram.")
			return
		pts=np.array(pts)
		pca=PCA(3)
		pca.fit(pts);
		c=pca.components_
		t=Transform()
		cc=c[2]
		if cc[2]!=0:
			cc*=np.sign(cc[2])
		
		t.set_rotation(c[2].tolist())
		t.invert()
		xyz=t.get_params("xyz")
		xyz["ztilt"]=0
		print("xtilt {:.02f}, ytilt {:.02f}".format(xyz["xtilt"], xyz["ytilt"]))
		t=Transform(xyz)
		self.globalxf=t
		if self.wfilt.getValue()!=0.0:
			data=self.datalp
		else:
			data=self.data
		self.dataxf=data.process("xform",{"transform":t})
		
		self.xyview.shapes={}
		self.zyview.shapes={}
		self.xzview.shapes={}
		
		boxes=self.get_rotated_boxes()
		for i,b in enumerate(boxes):
			self.update_box_shape(i,b)
		
		self.update_sliceview()
		print("Done")
	
	def reset_flatten_tomo(self, event):
		self.globalxf=Transform()
		self.xyview.shapes={}
		self.zyview.shapes={}
		self.xzview.shapes={}
		
		boxes=self.get_rotated_boxes()
		for i,b in enumerate(boxes):
			self.update_box_shape(i,b)
		
		self.update_sliceview()
		

	def SaveJson(self):
		
		info=js_open_dict(self.jsonfile)
		sx,sy,sz=(self.data["nx"]//2,self.data["ny"]//2,self.data["nz"]//2)
		if "apix_unbin" in info:
			bxs=[]
			for b0 in self.boxes:
				b=[	(b0[0]-sx)*self.apix_cur/self.apix_unbin,
					(b0[1]-sy)*self.apix_cur/self.apix_unbin,
					(b0[2]-sz)*self.apix_cur/self.apix_unbin,
					b0[3], b0[4], b0[5]	]
				bxs.append(b)
				
			bxsz={}
			for k in self.boxsize.keys():
				bxsz[k]=np.round(self.boxsize[k]*self.apix_cur/self.apix_unbin)

				
		else:
			bxs=self.boxes
			bxsz=self.boxsize
				
		info["boxes_3d"]=bxs
		clslst={}
		for key in list(self.sets.keys()):
			clslst[int(key)]={
				"name":self.sets[key],
				"boxsize":int(bxsz[key]),
				}
		info["class_list"]=clslst
		info.close()
	
	def closeEvent(self,event):
		print("Exiting")
		self.SaveJson()
		
		E2saveappwin("e2sptboxer","main",self)
		E2saveappwin("e2sptboxer","boxes",self.boxesviewer.qt_parent)
		E2saveappwin("e2sptboxer","option",self.optionviewer)
		
		#self.boxviewer.close()
		self.boxesviewer.close()
		self.optionviewer.close()
		#self.optionviewer.close()
		try:
			self.xyview.close()
			self.xzview.close()
			self.zyview.close()
		except:
			pass
		
		self.module_closed.emit() # this signal is important when e2ctf is being used by a program running its own event loop

def parse_setname(name):
	p0=name.find('::')
	ret=-1
	if p0>0:
		try:
			ret=int(name[:p0])
		except:
			pass
	
	return ret
			
class EMTomoBoxerOptions(QtWidgets.QWidget):
	def __init__(self,target) :
		QtWidgets.QWidget.__init__(self)
		#print "aaaaaaaa"
		self.setWindowTitle("Options")
		self.target=weakref.ref(target)
		
		self.gbl = QtWidgets.QGridLayout(self)
		#self.gbl.setContentsMargins(2, 2, 2, 2)
		#self.gbl.setSpacing(6)
		self.gbl.setObjectName("gbl")
		
		
		self.erasercheckbox=QtWidgets.QCheckBox("Eraser")
		self.gbl.addWidget(self.erasercheckbox,0,0)
		
		self.eraser_radius=ValBox(label="Radius:",value=64)
		self.gbl.addWidget(self.eraser_radius,0,1)

		self.tabwidget = QtWidgets.QTabWidget()
		self.gbl.addWidget(self.tabwidget,1,0,1,2)
		
	def add_panel(self,widget,name):
		self.tabwidget.addTab(widget,name)

		return 

#### Copied from emimagemx.py since some modification are needed...
class EMTomoSetsPanel(QtWidgets.QWidget):
	'''
	This is the set display panel
	'''
	def __init__(self,target):
		QtWidgets.QWidget.__init__(self)

		self.target = weakref.ref(target) # this should be the EMImageMXWidget
		self.busy = False
		self.initialized=False

		# cached values for speed later
		self.itemflags=	Qt.ItemFlags(Qt.ItemIsEditable)|Qt.ItemFlags(Qt.ItemIsSelectable)|Qt.ItemFlags(Qt.ItemIsEnabled)|Qt.ItemFlags(Qt.ItemIsUserCheckable)

		# now build the interface
		hbl = QtWidgets.QHBoxLayout(self)
		self.setlist=QtWidgets.QListWidget()
		self.setlist.setSizePolicy(QtWidgets.QSizePolicy.Preferred,QtWidgets.QSizePolicy.Expanding)
		hbl.addWidget(self.setlist)

		vbl = QtWidgets.QVBoxLayout()

		self.new_set_button = QtWidgets.QPushButton("New")
		vbl.addWidget(self.new_set_button)
		self.rename_set_button = QtWidgets.QPushButton("Rename")
		vbl.addWidget(self.rename_set_button)
		self.save_set_button = QtWidgets.QPushButton("Save")
		vbl.addWidget(self.save_set_button)
		self.delete_set_button = QtWidgets.QPushButton("Delete")
		vbl.addWidget(self.delete_set_button)

		hbl.addLayout(vbl)

		self.save_set_button.clicked[bool].connect(self.save_set)
		self.new_set_button.clicked[bool].connect(self.new_set)
		self.rename_set_button.clicked[bool].connect(self.rename_set)
		self.delete_set_button.clicked[bool].connect(self.delete_set)
		self.setlist.itemChanged[QtWidgets.QListWidgetItem].connect(self.set_list_item_changed)
		self.setlist.currentRowChanged[int].connect(self.set_list_row_changed)


	def sets_changed(self):
		self.update_sets()
		#keys=sorted(self.target().sets.keys())
		#for i,k in enumerate(keys):
			#try:
				#if k!=str(self.setlist.item(i).text()) : raise Exception
			#except:
				#self.update_sets()
				#break

	def set_list_row_changed(self,i):
		#print(i)
		if not self.initialized: return 
		a = self.setlist.item(i)
		if a==None : return
		name = str(a.text())
		self.target().set_current_set(name)
		#self.update_sets()

	def set_list_item_changed(self,item):
		name=str(item.text())
		#print(name)
		if item.checkState() == Qt.Checked : 
			self.target().show_set(name)
		else: 
			self.target().hide_set(name)
		

	def delete_set(self,unused):
		selections = self.setlist.selectedItems()
		if len(selections)==0 : return
		names=[str(i.text()) for i in selections]
		cancel=QtWidgets.QMessageBox.question(self, "Delete set", "Are you sure to delete {}? This will remove all particles in that class".format(names[0]),QtWidgets.QMessageBox.Yes, QtWidgets.QMessageBox.No)
		#print(cancel, QtWidgets.QMessageBox.Yes, QtWidgets.QMessageBox.No)
		if cancel==QtWidgets.QMessageBox.Yes :
			self.target().delete_set(names[0])
			self.update_sets()


	def new_set(self,unused=None):
		name,ok=QtWidgets.QInputDialog.getText( self, "Set Name", "Enter a name for the new set:")
		if not ok : return
		name=str(name)
		if name in self.target().sets :
			print("Set name exists")
			return

		self.target().new_set(name)
		self.update_sets()

	def rename_set(self,unused=None):
		selections = self.setlist.selectedItems()
		sels=[str(i.text()) for i in selections]
		if len(sels)==0:
			return
		name,ok=QtWidgets.QInputDialog.getText( self, "Set Name", "Enter a name for the new set:")
		if not ok : return
		name=str(name)
		
		if name in self.target().sets :
			print("Set name exists")
			return
		
		self.target().rename_set(sels[0], name)
		self.update_sets()

	def save_set(self):
		self.target().save_set()
		self.update_sets()

	def update_sets(self):
		keys=sorted(self.target().sets.keys())
		viskeys=set(self.target().sets_visible.keys())
		self.setlist.clear()

		for i,k in enumerate(keys):
			
			kname="{:02d} :: {}".format(int(k), self.target().sets[k])
			item=QtWidgets.QListWidgetItem(kname)
			item.setFlags(self.itemflags)
			item.setForeground(self.target().setcolors[i%len(self.target().setcolors)])
			self.setlist.addItem(item)
			if k in viskeys : item.setCheckState(Qt.Checked)
			else : item.setCheckState(Qt.Unchecked)
			
			if not self.initialized:
				if k==self.target().currentset:
					self.setlist.setCurrentItem(item)
					self.initialized=True
		return 


if __name__ == '__main__':
	main()
