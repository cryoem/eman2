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




from __future__ import print_function
from __future__ import division
from past.utils import old_div
from builtins import range
from EMAN2 import *
from EMAN2_utils import numpy2pdb
import numpy as np

import weakref
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtCore import Qt
from eman2_gui.emapplication import get_application, EMApp
from eman2_gui.emimage2d import EMImage2DWidget
from eman2_gui.emimagemx import EMImageMXWidget
#from emimage3d import EMImage3DWidget
from eman2_gui.emscene3d import EMScene3D
from eman2_gui.emdataitem3d import EMDataItem3D, EMIsosurface
from eman2_gui.emshape import EMShape
from eman2_gui.valslider import ValSlider, ValBox

	
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
		self.yshort=False
		self.apix=options.apix
		self.currentset=0
		self.shrink=1#options.shrink
		self.setWindowTitle("Main Window (e2spt_boxer.py)")
		if options.mode=="3D":
			self.boxshape="circle"
		else:
			self.boxshape="rect"


		# Menu Bar
		self.mfile=self.menuBar().addMenu("File")
		self.mfile_open=self.mfile.addAction("Open")
		self.mfile_read_boxloc=self.mfile.addAction("Read Box Coord")
		self.mfile_save_boxloc=self.mfile.addAction("Save Box Coord")
		self.mfile_save_boxpdb=self.mfile.addAction("Save Coord as PDB")
		self.mfile_save_boxes_stack=self.mfile.addAction("Save Boxes as Stack")
		self.mfile_quit=self.mfile.addAction("Quit")


		self.setCentralWidget(QtWidgets.QWidget())
		self.gbl = QtWidgets.QGridLayout(self.centralWidget())

		# relative stretch factors
		self.gbl.setColumnStretch(0,1)
		self.gbl.setColumnStretch(1,4)
		self.gbl.setColumnStretch(2,0)
		self.gbl.setRowStretch(1,1)
		self.gbl.setRowStretch(0,4)

		# 3 orthogonal restricted projection views
		self.xyview = EMImage2DWidget()
		self.gbl.addWidget(self.xyview,0,1)

		self.xzview = EMImage2DWidget()
		self.gbl.addWidget(self.xzview,1,1)

		self.zyview = EMImage2DWidget()
		self.gbl.addWidget(self.zyview,0,0)

		# Select Z for xy view
		self.wdepth = QtWidgets.QSlider()
		self.gbl.addWidget(self.wdepth,1,2)

		### Control panel area in upper left corner
		self.gbl2 = QtWidgets.QGridLayout()
		self.gbl.addLayout(self.gbl2,1,0)

		#self.wxpos = QtWidgets.QSlider(Qt.Horizontal)
		#self.gbl2.addWidget(self.wxpos,0,0)
		
		#self.wypos = QtWidgets.QSlider(Qt.Vertical)
		#self.gbl2.addWidget(self.wypos,0,3,6,1)
		
		# box size
		self.wboxsize=ValBox(label="Box Size:",value=0)
		self.gbl2.addWidget(self.wboxsize,2,0,1,2)

		# max or mean
		#self.wmaxmean=QtWidgets.QPushButton("MaxProj")
		#self.wmaxmean.setCheckable(True)
		#self.gbl2.addWidget(self.wmaxmean,3,0)

		# number slices
		self.wnlayers=QtWidgets.QSpinBox()
		self.wnlayers.setMinimum(1)
		self.wnlayers.setMaximum(256)
		self.wnlayers.setValue(1)
		self.gbl2.addWidget(self.wnlayers,3,1)

		# Local boxes in side view
		self.wlocalbox=QtWidgets.QCheckBox("Limit Side Boxes")
		self.gbl2.addWidget(self.wlocalbox,3,0)
		self.wlocalbox.setChecked(True)

		# scale factor
		self.wscale=ValSlider(rng=(.1,2),label="Sca:",value=1.0)
		self.gbl2.addWidget(self.wscale,4,0,1,2)

		# 2-D filters
		self.wfilt = ValSlider(rng=(0,150),label="Filt:",value=0.0)
		self.gbl2.addWidget(self.wfilt,5,0,1,2)
		
		self.curbox=-1
		
		self.boxes=[]						# array of box info, each is (x,y,z,...)
		self.boxesimgs=[]					# z projection of each box
		self.xydown=self.xzdown=self.zydown=None
		self.firsthbclick = None

		# coordinate display
		self.wcoords=QtWidgets.QLabel("X: " + str(self.get_x()) + "\t\t" + "Y: " + str(self.get_y()) + "\t\t" + "Z: " + str(self.get_z()))
		self.gbl2.addWidget(self.wcoords, 1, 0, 1, 2)

		# file menu
		self.mfile_open.triggered[bool].connect(self.menu_file_open)
		self.mfile_read_boxloc.triggered[bool].connect(self.menu_file_read_boxloc)
		self.mfile_save_boxloc.triggered[bool].connect(self.menu_file_save_boxloc)
		self.mfile_save_boxpdb.triggered[bool].connect(self.menu_file_save_boxpdb)
		
		self.mfile_save_boxes_stack.triggered[bool].connect(self.save_boxes)
		self.mfile_quit.triggered[bool].connect(self.menu_file_quit)

		# all other widgets
		self.wdepth.valueChanged[int].connect(self.event_depth)
		self.wnlayers.valueChanged[int].connect(self.event_nlayers)
		self.wboxsize.valueChanged.connect(self.event_boxsize)
		#self.wmaxmean.clicked[bool].connect(self.event_projmode)
		self.wscale.valueChanged.connect(self.event_scale)
		self.wfilt.valueChanged.connect(self.event_filter)
		self.wlocalbox.stateChanged[int].connect(self.event_localbox)

		self.xyview.mousemove.connect(self.xy_move)
		self.xyview.mousedown.connect(self.xy_down)
		self.xyview.mousedrag.connect(self.xy_drag)
		self.xyview.mouseup.connect(self.xy_up)
		self.xyview.mousewheel.connect(self.xy_wheel)
		self.xyview.signal_set_scale.connect(self.xy_scale)
		self.xyview.origin_update.connect(self.xy_origin)

		self.xzview.mousedown.connect(self.xz_down)
		self.xzview.mousedrag.connect(self.xz_drag)
		self.xzview.mouseup.connect(self.xz_up)
		self.xzview.signal_set_scale.connect(self.xz_scale)
		self.xzview.origin_update.connect(self.xz_origin)

		self.zyview.mousedown.connect(self.zy_down)
		self.zyview.mousedrag.connect(self.zy_drag)
		self.zyview.mouseup.connect(self.zy_up)
		self.zyview.signal_set_scale.connect(self.zy_scale)
		self.zyview.origin_update.connect(self.zy_origin)
		
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
		
		# Average viewer shows results of background tomographic processing
#		self.averageviewer=EMAverageViewer(self)
		#self.averageviewer.show()

		self.boxesviewer.mx_image_selected.connect(self.img_selected)
		
		self.jsonfile=info_name(datafile)
		
		info=js_open_dict(self.jsonfile)
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
				
		if "boxes_3d" in info:
			box=info["boxes_3d"]
			for i,b in enumerate(box):
				#### X-center,Y-center,Z-center,method,[score,[class #]]
				bdf=[0,0,0,"manual",0.0, 0]
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
		if len(self.sets)==0:
			self.new_set("particles_00")
		self.sets_visible[list(self.sets.keys())[0]]=0
		self.currentset=sorted(self.sets.keys())[0]
		self.setspanel.update_sets()
		self.wboxsize.setValue(self.get_boxsize())

		print(self.sets)
		for i in range(len(self.boxes)):
			self.update_box(i)
		
		self.update_all()
		self.initialized=True

	def set_data(self,data):

		self.data=data
		self.apix=data["apix_x"]

		self.datasize=(data["nx"],data["ny"],data["nz"])

		self.wdepth.setRange(0,self.datasize[2]-1)
		self.boxes=[]
		self.curbox=-1

		self.wdepth.setValue(old_div(self.datasize[2],2))
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
			r["apix_x"]=self.apix
			r["apix_y"]=self.apix
			r["apix_z"]=self.apix

		#if options.normproc:
			#r.process_inplace(options.normproc)
		return r

	def get_slice(self,n,xyz):
		"""Reads a slice either from a file or the preloaded memory array.
		xyz is the axis along which 'n' runs, 0=x (yz), 1=y (xz), 2=z (xy)"""

		if xyz==0:
			r=self.data.get_clip(Region(n,0,0,1,self.datasize[1],self.datasize[2]))
			r.set_size(self.datasize[1],self.datasize[2],1)
		elif xyz==1:
			r=self.data.get_clip(Region(0,n,0,self.datasize[0],1,self.datasize[2]))
			r.set_size(self.datasize[0],self.datasize[2],1)
		else:
			r=self.data.get_clip(Region(0,0,n,self.datasize[0],self.datasize[1],1))

		if self.apix!=0 :
			r["apix_x"]=self.apix
			r["apix_y"]=self.apix
			r["apix_z"]=self.apix
		return r

	def event_boxsize(self):
		if self.get_boxsize()==int(self.wboxsize.getValue()):
			return
		
		self.boxsize[self.currentset]=int(self.wboxsize.getValue())
		
		cb=self.curbox
		self.initialized=False
		for i in range(len(self.boxes)):
			if self.boxes[i][5]==self.currentset:
				self.update_box(i)
		self.update_box(cb)
		self.initialized=True
		self.update_all()

	def event_projmode(self,state):
		"""Projection mode can be simple average (state=False) or maximum projection (state=True)"""
		self.update_all()

	def event_scale(self,newscale):
		self.xyview.set_scale(newscale)
		self.xzview.set_scale(newscale)
		self.zyview.set_scale(newscale)

	def event_depth(self):
		if self.initialized:
			self.update_xy()

	def event_nlayers(self):
		self.update_all()

	def event_filter(self):
		self.update_all()

	def event_localbox(self,tog):
		self.update_sides()

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

	def depth(self):
		return int(self.wdepth.value())

	def scale(self):
		return self.wscale.getValue()

	def get_x(self):
		return self.get_coord(0)

	def get_y(self):
		return self.get_coord(1)

	def get_z(self):
		return self.depth()

	def get_coord(self, coord_index):
		if len(self.boxes) > 1:
			if self.curbox:
				return self.boxes[self.curbox][coord_index]
			else:
				return self.boxes[-1][coord_index]
		else:
			return 0


	def menu_file_open(self,tog):
		QtWidgets.QMessageBox.warning(None,"Error","Sorry, in the current version, you must provide a file to open on the command-line.")

	def load_box_yshort(self, boxcoords):
		if options.yshort:
			return [boxcoords[0], boxcoords[2], boxcoords[1]]
		else:
			return boxcoords

	def menu_file_read_boxloc(self):
		fsp=str(QtWidgets.QFileDialog.getOpenFileName(self, "Select output text file")[0])
		
		if not os.path.isfile(fsp):
			print("file does not exist")
			return

		f=file(fsp,"r")
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

		out=file(fsp,"w")
		for b in self.boxes:
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


	def menu_file_quit(self):
		self.close()

	def transform_coords(self, point, xform):
		xvec = xform.get_matrix()
		return [xvec[0]*point[0] + xvec[4]*point[1] + xvec[8]*point[2] + xvec[3], xvec[1]*point[0] + xvec[5]*point[1] + xvec[9]*point[2] + xvec[7], xvec[2]*point[0] + xvec[6]*point[1] + xvec[10]*point[2] + xvec[11]]

	def get_averager(self):
		"""returns an averager of the appropriate type for generating projection views"""
		#if self.wmaxmean.isChecked() : return Averagers.get("minmax",{"max":1})

		return Averagers.get("mean")

	def update_sides(self):
		"""updates xz and yz views due to a new center location"""

		#print "\n\n\n\n\nIn update sides, self.datafile is", self.datafile
		#print "\n\n\n\n"

		if self.data==None:
			return

		if self.curbox==-1 :
			x=self.datasize[0]//2
			y=self.datasize[1]//2
			z=0
		else:
			x,y,z=self.boxes[self.curbox][:3]

		self.cury=y
		self.curx=x

		# update shape display
		if self.wlocalbox.isChecked():
			xzs=self.xzview.get_shapes()
			for i in range(len(self.boxes)):
				bs=self.get_boxsize(self.boxes[i][5])
				if self.boxes[i][1]<self.cury+old_div(bs,2) and self.boxes[i][1]>self.cury-old_div(bs,2) and  self.boxes[i][5] in self.sets_visible:
					xzs[i][0]=self.boxshape
				else:
					xzs[i][0]="hidden"

			zys=self.zyview.get_shapes()
			
			for i in range(len(self.boxes)):
				bs=self.get_boxsize(self.boxes[i][5])
				if self.boxes[i][0]<self.curx+old_div(bs,2) and self.boxes[i][0]>self.curx-old_div(bs,2) and  self.boxes[i][5] in self.sets_visible:
					zys[i][0]=self.boxshape
				else:
					zys[i][0]="hidden"
		else :
			xzs=self.xzview.get_shapes()
			zys=self.zyview.get_shapes()
		
			for i in range(len(self.boxes)):
				bs=self.get_boxsize(self.boxes[i][5])
				if  self.boxes[i][5] in self.sets_visible:
					xzs[i][0]=self.boxshape
					zys[i][0]=self.boxshape
				else:
					xzs[i][0]="hidden"
					zys[i][0]="hidden"

		self.xzview.shapechange=1
		self.zyview.shapechange=1

		# yz
		avgr=self.get_averager()

		for x in range(x-(self.nlayers()//2),x+((self.nlayers()+1)//2)):
			slc=self.get_slice(x,0)
			avgr.add_image(slc)

		av=avgr.finish()
		if not self.yshort:
			av.process_inplace("xform.transpose")

		if self.wfilt.getValue()!=0.0:
			av.process_inplace("filter.lowpass.gauss",{"cutoff_freq":old_div(1.0,self.wfilt.getValue()),"apix":self.apix})

		self.zyview.set_data(av)

		# xz
		avgr=self.get_averager()

		for y in range(y-old_div(self.nlayers(),2),y+old_div((self.nlayers()+1),2)):
			slc=self.get_slice(y,1)
			avgr.add_image(slc)

		av=avgr.finish()
		if self.wfilt.getValue()!=0.0:
			av.process_inplace("filter.lowpass.gauss",{"cutoff_freq":old_div(1.0,self.wfilt.getValue()),"apix":self.apix})

		self.xzview.set_data(av)


	def update_xy(self):
		"""updates xy view due to a new slice range"""

		#print "\n\n\n\n\nIn update_xy, self.datafile is", self.datafile
		#print "\n\n\n\n"

		if self.data==None:
			return

		# Boxes should also be limited by default in the XY view
		if len(self.boxes) > 0:
			zc=self.wdepth.value()
			#print "The current depth is", self.wdepth.value()
			xys=self.xyview.get_shapes()
			for i in range(len(self.boxes)):

				bs=self.get_boxsize(self.boxes[i][5])
				zdist=abs(self.boxes[i][2] - zc)

				if self.options.mode=="3D":
					zthr=bs/2
					xys[i][6]=bs//2-zdist
				else:
					zthr=1
					
				if zdist < zthr and self.boxes[i][5] in self.sets_visible:
					xys[i][0]=self.boxshape
					
				else :
					xys[i][0]="hidden"
			self.xyview.shapechange=1

		#if self.wmaxmean.isChecked():
			#avgr=Averagers.get("minmax",{"max":1})

		#else:
		avgr=Averagers.get("mean")

		slc=EMData()
		for z in range(self.wdepth.value()-self.nlayers()//2,self.wdepth.value()+(self.nlayers()+1)//2):
			slc=self.get_slice(z,2)
			avgr.add_image(slc)

		av=avgr.finish()

		#print "\n\nIn update xy, av and type are", av, type(av)

		if self.wfilt.getValue()!=0.0:

			av.process_inplace("filter.lowpass.gauss",{"cutoff_freq":old_div(1.0,self.wfilt.getValue()),"apix":self.apix})
		if self.initialized:
			self.xyview.set_data(av, keepcontrast=True)
		else:
			self.xyview.set_data(av)

	def update_all(self):
		"""redisplay of all widgets"""

		#print "\n\n\n\n\nIn update all, self.datafile is", self.datafile
		#print "\n\n\n\n"
		if self.data==None:
			return

		self.update_xy()
		self.update_sides()
		self.update_boximgs()

	def update_coords(self):
		self.wcoords.setText("X: " + str(self.get_x()) + "\t\t" + "Y: " + str(self.get_y()) + "\t\t" + "Z: " + str(self.get_z()))

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

	def do_deletion(self, delids):
		
		kpids=[i for i,b in enumerate(self.boxes) if i not in delids]
		self.boxes=[self.boxes[i] for i in kpids]
		self.boxesimgs=[self.boxesimgs[i] for i in kpids]
		self.xyview.shapes={i:self.xyview.shapes[k] for i,k in enumerate(kpids)}
		self.xzview.shapes={i:self.xzview.shapes[k] for i,k in enumerate(kpids)}
		self.zyview.shapes={i:self.zyview.shapes[k] for i,k in enumerate(kpids)}
		#print self.boxes, self.xyview.get_shapes()
		self.curbox=-1
		self.update_all()

	def del_box(self,n):
		"""Delete an existing box by replacing the deleted box with the last box. A bit funny, but otherwise
		update after deletion is REALLY slow."""
#		print "del ",n
		if n<0 or n>=len(self.boxes): return

		#if self.boxviewer.get_data(): self.boxviewer.set_data(None)
		self.curbox=-1
		self.do_deletion([n])




	def update_box(self,n,quiet=False):
		"""After adjusting a box, call this"""
#		print "upd ",n,quiet

		try:
			box=self.boxes[n]
		except IndexError:
			return
		bs2=self.get_boxsize(box[5])//2

		
		color=self.setcolors[box[5]%len(self.setcolors)].color().getRgbF()
		if self.options.mode=="3D":
			self.xyview.add_shape(n,EMShape(["circle",color[0],color[1],color[2],box[0],box[1],bs2,2]))
			self.xzview.add_shape(n,EMShape(["circle",color[0],color[1],color[2],box[0],box[2],bs2,2]))
			self.zyview.add_shape(n,EMShape(("circle",color[0],color[1],color[2],box[2],box[1],bs2,2)))
		else:
			self.xyview.add_shape(n,EMShape(["rect",color[0],color[1],color[2],
				    box[0]-bs2,box[1]-bs2,box[0]+bs2,box[1]+bs2,2]))
			self.xzview.add_shape(n,EMShape(["rect",color[0],color[1],color[2], 
				    box[0]-bs2,box[2]-1,box[0]+bs2,box[2]+1,2]))
			self.zyview.add_shape(n,EMShape(["rect",color[0],color[1],color[2],
				    box[2]-1,box[1]-bs2,box[2]+1,box[1]+bs2,2]))
			
			

		if self.depth()!=box[2]:
			self.wdepth.setValue(box[2])
		else:
			self.xyview.update()
		if self.initialized: self.update_sides()

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
			

			if n!=self.curbox:
				self.boxesviewer.set_selected((n,),True)

		self.curbox=n
		self.update_coords()

	def update_boximgs(self):
		self.boxids=[im for im,m in enumerate(self.boxesimgs) if self.boxes[im][5] in self.sets_visible]
		self.boxesviewer.set_data([self.boxesimgs[i] for i in self.boxids])
		self.boxesviewer.update()
		return

	def img_selected(self,event,lc):
		#print "sel",lc[0]
		lci=self.boxids[lc[0]]
		if event.modifiers()&Qt.ShiftModifier:
			self.del_box(lci)
		else:
			self.update_box(lci)
		if self.curbox>=0 :
			box=self.boxes[self.curbox]
			self.xyview.scroll_to(box[0],box[1])
			self.xzview.scroll_to(None,box[2])
			self.zyview.scroll_to(box[2],None)
			self.currentset=box[5]
			self.setspanel.initialized=False
			self.setspanel.update_sets()

	def del_region_xy(self, x=-1, y=-1, z=-1, rad=-1):
		if rad<0:
			rad=self.eraser_width()
		
		delids=[]
		for i,b in enumerate(self.boxes):
			if b[5] not in self.sets_visible:
				continue
			
			if (x>=0)*(b[0]-x)**2 + (y>=0)*(b[1]-y)**2 +(z>=0)*(b[2]-z)**2 < rad**2:
				delids.append(i)
		self.do_deletion(delids)

	def xy_down(self,event):
		x,y=self.xyview.scr_to_img((event.x(),event.y()))
		x,y=int(x),int(y)
		z=int(self.get_z())
		self.xydown=None
		if x<0 or y<0 : return		# no clicking outside the image (on 2 sides)
		if self.optionviewer.erasercheckbox.isChecked():
			self.del_region_xy(x,y)
			return
			
		for i in range(len(self.boxes)):
			if self.inside_box(i,x,y,z):
				if event.modifiers()&Qt.ShiftModifier:
					self.del_box(i)
					self.firsthbclick = None
				else:
					self.xydown=(i,x,y,self.boxes[i][0],self.boxes[i][1])
					self.update_box(i)
				break
		else:
#			if x>self.get_boxsize()/2 and x<self.datasize[0]-self.get_boxsize()/2 and y>self.get_boxsize()/2 and y<self.datasize[1]-self.get_boxsize()/2 and self.depth()>self.get_boxsize()/2 and self.depth()<self.datasize[2]-self.get_boxsize()/2 :
			if not event.modifiers()&Qt.ShiftModifier:
				self.boxes.append(([x,y,self.depth(), 'manual', 0.0, self.currentset]))
				self.xydown=(len(self.boxes)-1,x,y,x,y)		# box #, x down, y down, x box at down, y box at down
				self.update_box(self.xydown[0])

		if self.curbox>=0:
			box=self.boxes[self.curbox]
			self.xzview.scroll_to(None,box[2])
			self.zyview.scroll_to(box[2],None)

	def xy_drag(self,event):
		
		x,y=self.xyview.scr_to_img((event.x(),event.y()))
		x,y=int(x),int(y)
		if self.optionviewer.erasercheckbox.isChecked():
			self.del_region_xy(x,y)
			self.xyview.eraser_shape=EMShape(["circle",1,1,1,x,y,self.eraser_width(),2])
			self.xyview.shapechange=1
			self.xyview.update()
			return
		
		if self.xydown==None : return


		dx=x-self.xydown[1]
		dy=y-self.xydown[2]

		self.boxes[self.xydown[0]][0]=dx+self.xydown[3]
		self.boxes[self.xydown[0]][1]=dy+self.xydown[4]
		self.update_box(self.curbox,True)

	def xy_up  (self,event):
		if self.xydown!=None: self.update_box(self.curbox)
		self.xydown=None

	def xy_wheel (self,event):
		if event.angleDelta().y() > 0:
			#self.wdepth.setValue(self.wdepth.value()+4)
			self.wdepth.setValue(self.wdepth.value()+1) #jesus

		elif event.angleDelta().y() < 0:
			#self.wdepth.setValue(self.wdepth.value()-4)
			self.wdepth.setValue(self.wdepth.value()-1) #jesus


	def xy_scale(self,news):
		"xy image view has been rescaled"
		self.wscale.setValue(news)
		#self.xzview.set_scale(news,True)
		#self.zyview.set_scale(news,True)

	def xy_origin(self,newor):
		"xy origin change"
		xzo=self.xzview.get_origin()
		self.xzview.set_origin(newor[0],xzo[1],True)

		zyo=self.zyview.get_origin()
		self.zyview.set_origin(zyo[0],newor[1],True)
	
	def xy_move(self,event):
		if self.optionviewer.erasercheckbox.isChecked():
			x,y=self.xyview.scr_to_img((event.x(),event.y()))
			#print x,y
			self.xyview.eraser_shape=EMShape(["circle",1,1,1,x,y,self.eraser_width(),2])
			self.xyview.shapechange=1
			self.xyview.update()
		else:
			self.xyview.eraser_shape=None

	def xz_down(self,event):
		x,z=self.xzview.scr_to_img((event.x(),event.y()))
		x,z=int(x),int(z)
		y=int(self.get_y())
		self.xzdown=None
		if x<0 or z<0 : return		# no clicking outside the image (on 2 sides)
		if self.optionviewer.erasercheckbox.isChecked():
			return
		for i in range(len(self.boxes)):
			if (not self.wlocalbox.isChecked() and self.inside_box(i,x,y,z)) or self.inside_box(i,x,self.cury,z) :
				if event.modifiers()&Qt.ShiftModifier:
					self.del_box(i)
					self.firsthbclick = None
				else :
					self.xzdown=(i,x,z,self.boxes[i][0],self.boxes[i][2])
					self.update_box(i)
				break
		else:
			if not event.modifiers()&Qt.ShiftModifier:
				self.boxes.append(([x,self.cury,z, 'manual', 0.0, self.currentset]))
				self.xzdown=(len(self.boxes)-1,x,z,x,z)		# box #, x down, y down, x box at down, y box at down
				self.update_box(self.xzdown[0])

		if self.curbox>=0 :
			box=self.boxes[self.curbox]
			self.xyview.scroll_to(None,box[1])
			self.zyview.scroll_to(box[2],None)

	def xz_drag(self,event):
		if self.xzdown==None : return

		x,z=self.xzview.scr_to_img((event.x(),event.y()))
		x,z=int(x),int(z)

		dx=x-self.xzdown[1]
		dz=z-self.xzdown[2]

		self.boxes[self.xzdown[0]][0]=dx+self.xzdown[3]
		self.boxes[self.xzdown[0]][2]=dz+self.xzdown[4]
		self.update_box(self.curbox,True)

	def xz_up  (self,event):
		if self.xzdown!=None: self.update_box(self.curbox)
		self.xzdown=None

	def xz_scale(self,news):
		"xy image view has been rescaled"
		self.wscale.setValue(news)
		#self.xyview.set_scale(news,True)
		#self.zyview.set_scale(news,True)

	def xz_origin(self,newor):
		"xy origin change"
		xyo=self.xyview.get_origin()
		self.xyview.set_origin(newor[0],xyo[1],True)

		#zyo=self.zyview.get_origin()
		#self.zyview.set_origin(zyo[0],newor[1],True)


	def zy_down(self,event):
		z,y=self.zyview.scr_to_img((event.x(),event.y()))
		z,y=int(z),int(y)
		x=int(self.get_x())
		self.xydown=None
		if z<0 or y<0 : return		# no clicking outside the image (on 2 sides)

		for i in range(len(self.boxes)):
			if (not self.wlocalbox.isChecked() and self.inside_box(i,x,y,z)) or  self.inside_box(i,self.curx,y,z):
				if event.modifiers()&Qt.ShiftModifier:
					self.del_box(i) 
					self.firsthbclick = None
				else :
					self.zydown=(i,z,y,self.boxes[i][2],self.boxes[i][1])
					self.update_box(i)
				break
		else:
			if not event.modifiers()&Qt.ShiftModifier:
				###########
				self.boxes.append(([self.curx,y,z, 'manual', 0.0, self.currentset]))
				self.zydown=(len(self.boxes)-1,z,y,z,y)		# box #, x down, y down, x box at down, y box at down
				self.update_box(self.zydown[0])

		if self.curbox>=0 :
			box=self.boxes[self.curbox]
			self.xyview.scroll_to(box[0],None)
			self.xzview.scroll_to(None,box[2])

	def zy_drag(self,event):
		if self.zydown==None : return

		z,y=self.zyview.scr_to_img((event.x(),event.y()))
		z,y=int(z),int(y)

		dz=z-self.zydown[1]
		dy=y-self.zydown[2]

		self.boxes[self.zydown[0]][2]=dz+self.zydown[3]
		self.boxes[self.zydown[0]][1]=dy+self.zydown[4]
		self.update_box(self.curbox,True)

	def zy_up  (self,event):
		if self.zydown!=None:
			self.update_box(self.curbox)
		self.zydown=None

	def zy_scale(self,news):
		"xy image view has been rescaled"
		self.wscale.setValue(news)
		#self.xyview.set_scale(news,True)
		#self.xzview.set_scale(news,True)

	def zy_origin(self,newor):
		"xy origin change"
		xyo=self.xyview.get_origin()
		self.xyview.set_origin(xyo[0],newor[1],True)

		#xzo=self.xzview.get_origin()
		#self.xzview.set_origin(xzo[0],newor[1],True)
	
	
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
		self.do_deletion(delids)
		
		if name in self.sets_visible: self.sets_visible.pop(name)
		if name in self.sets: self.sets.pop(name)
		if name in self.boxsize: self.boxsize.pop(name)
		
		self.curbox=-1
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
		if event.key() == 96:
			self.wdepth.setValue(self.wdepth.value()+1)

		elif event.key() == 49:
			self.wdepth.setValue(self.wdepth.value()-1)
		else:
			self.keypress.emit(event)

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
		
		#self.boxviewer.close()
		self.boxesviewer.close()
		self.optionviewer.close()
		self.xyview.close()
		self.xzview.close()
		self.zyview.close()
		
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
			
	

class EMBoxViewer(QtWidgets.QWidget):
	"""This is a multi-paned view showing a single boxed out particle from a larger tomogram"""

	def __init__(self):
		QtWidgets.QWidget.__init__(self)
		self.setWindowTitle("Single Particle View")

		self.resize(300,300)

		self.gbl = QtWidgets.QGridLayout(self)
		self.xyview = EMImage2DWidget()
		self.gbl.addWidget(self.xyview,0,1)

		self.xzview = EMImage2DWidget()
		self.gbl.addWidget(self.xzview,1,1)

		self.zyview = EMImage2DWidget()
		self.gbl.addWidget(self.zyview,0,0)
		self.data = None


		# This puts an isosurface view in the lower left corner, but was causing a lot of segfaults, so switching to 2-D slices for now
		#self.d3view = EMScene3D()
		#self.d3viewdata = EMDataItem3D(test_image_3d(3), transform=Transform())
		#isosurface = EMIsosurface(self.d3viewdata, transform=Transform())
		#self.d3view.insertNewNode('', self.d3viewdata, parentnode=self.d3view)
		#self.d3view.insertNewNode("Iso", isosurface, parentnode=self.d3viewdata )

		self.d3view = EMImage2DWidget()
		self.gbl.addWidget(self.d3view,1,0)

		self.wfilt = ValSlider(rng=(0,50),label="Filter:",value=0.0)
		self.gbl.addWidget(self.wfilt,2,0,1,2)

		self.wfilt.valueChanged.connect(self.event_filter)

		self.gbl.setRowStretch(2,1)
		self.gbl.setRowStretch(0,5)
		self.gbl.setRowStretch(1,5)
		
	def set_data(self,data):
		"""Sets the current volume to display"""

		self.data=data
		self.fdata=data

		self.update()
		self.show()

	def get_data(self):
		return self.data

	def update(self):
		if self.data==None:
			self.xyview.set_data(None)
			self.xzview.set_data(None)
			self.zyview.set_data(None)

			#self.d3viewdata.setData(test_image_3d(3))
			#self.d3view.updateSG()
			self.d3view.set_data(test_image_3d(3))

			return

		if self.wfilt.getValue()>4 :
			self.fdata=self.data.process("filter.lowpass.gauss",{"cutoff_freq":old_div(1.0,self.wfilt.getValue()),"apix":self.data['apix_x']}) #JESUS

		xyd=self.fdata.process("misc.directional_sum",{"axis":"z"})
		xzd=self.fdata.process("misc.directional_sum",{"axis":"y"})
		zyd=self.fdata.process("misc.directional_sum",{"axis":"x"})

		self.xyview.set_data(xyd)
		self.xzview.set_data(xzd)
		self.zyview.set_data(zyd)

		#self.d3viewdata.setData(self.fdata)
		#self.d3view.updateSG()
		self.d3view.set_data(self.fdata)


	def event_filter(self,value):
		self.update()

	def closeEvent(self, event):
		self.d3view.close()
		self.xyview.close()
		self.xzview.close()
		self.zyview.close()


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
