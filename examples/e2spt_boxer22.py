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




from EMAN2 import *
import numpy as np

import weakref
from PyQt4 import QtCore, QtGui
from PyQt4.QtCore import Qt
from emapplication import get_application, EMApp
from emimage2d import EMImage2DWidget
from emimagemx import EMImageMXWidget
from emimage3d import EMImage3DWidget
from emscene3d import EMScene3D
from emdataitem3d import EMDataItem3D, EMIsosurface
from emshape import EMShape
from valslider import ValSlider, ValBox

	
def run(cmd):
	print cmd
	launch_childprocess(cmd)
	
def main():
	
	usage=""" !!! Experimental program....
	Cleaned and modified version of e2spt_boxer.py. Allow multiple type of features in one tomogram and save metadata in the new EMAN2 tomogram framework.
	
	[prog] <tomogram input> --inmemory --invert
	
	"""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	#parser.add_argument("--path", type=str,help="path", default=None)	
	#parser.add_argument("--shrink", type=int,help="shrink", default=1)	
	parser.add_argument("--apix", type=float,help="apix", default=-1)	
	parser.add_argument('--invert', action="store_true", default=False, help='''This means you want the contrast to me inverted while boxing, AND for the extracted sub-volumes.\nRemember that EMAN2 **MUST** work with "white" protein. You can very easily figure out what the original color\nof the protein is in your data by looking at the gold fiducials or the edge of the carbon hole in your tomogram.\nIf they look black you MUST specify this option''', guitype='boolbox', row=4, col=0, rowspan=1, colspan=1, mode="boxing")
	parser.add_argument("--inmemory",action="store_true",default=False,help="This will read the entire tomogram into memory. Much faster, but you must have enough ram !", guitype='boolbox', row=2, col=1, rowspan=1, colspan=1, mode="boxing")

	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)


	img = args[0]

	app = EMApp()

	img=args[0]

	imghdr = EMData(img,0,True)
	if options.apix<=0:
		options.apix = imghdr['apix_x']

	box = 32
	
	boxer=EMTomoBoxer(app,options,datafile=img)


	boxer.show()
	app.execute()
	E2end(logid)
	return

class EMTomoBoxer(QtGui.QMainWindow):
	"""This class represents the EMTomoBoxer application instance.  """

	def __init__(self,application,options,datafile=None):
		QtGui.QWidget.__init__(self)
		self.initialized=False
		self.app=weakref.ref(application)
		self.options=options
		boxsize=32
		self.helixboxer=False
		self.yshort=False
		self.apix=options.apix
		self.currentset=0
		self.shrink=1#options.shrink
		self.invert=options.invert
		self.center=None
		self.normalize=None
		self.setWindowTitle("Main Window (e2spt_boxer.py)")

#		self.setWindowTitle("e2spt_boxer.py")

		# Menu Bar
		self.mfile=self.menuBar().addMenu("File")
		self.mfile_open=self.mfile.addAction("Open")
		self.mfile_read_boxloc=self.mfile.addAction("Read Box Coord")
		self.mfile_save_boxloc=self.mfile.addAction("Save Box Coord")
		self.mfile_save_boxes_stack=self.mfile.addAction("Save Boxes as Stack")
		self.mfile_quit=self.mfile.addAction("Quit")

		self.mwin=self.menuBar().addMenu("Window")
		self.mwin_boxes=self.mwin.addAction("Particles")
		self.mwin_single=self.mwin.addAction("Single Particle")
		self.mwin_average=self.mwin.addAction("Averaging")


		self.setCentralWidget(QtGui.QWidget())
		self.gbl = QtGui.QGridLayout(self.centralWidget())

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
		self.wdepth = QtGui.QSlider()
		self.gbl.addWidget(self.wdepth,1,2)

		### Control panel area in upper left corner
		self.gbl2 = QtGui.QGridLayout()
		self.gbl.addLayout(self.gbl2,1,0)

		# box size
		self.wboxsize=ValBox(label="Box Size:",value=boxsize)
		self.gbl2.addWidget(self.wboxsize,1,0,1,2)

		# max or mean
		self.wmaxmean=QtGui.QPushButton("MaxProj")
		self.wmaxmean.setCheckable(True)
		self.gbl2.addWidget(self.wmaxmean,2,0)

		# number slices
		self.wnlayers=QtGui.QSpinBox()
		self.wnlayers.setMinimum(1)
		self.wnlayers.setMaximum(256)
		self.wnlayers.setValue(1)
		self.gbl2.addWidget(self.wnlayers,2,1)

		# Local boxes in side view
		self.wlocalbox=QtGui.QCheckBox("Limit Side Boxes")
		self.gbl2.addWidget(self.wlocalbox,3,0)

		# scale factor
		self.wscale=ValSlider(rng=(.1,2),label="Sca:",value=1.0)
		self.gbl2.addWidget(self.wscale,4,0,1,2)

		# 2-D filters
		self.wfilt = ValSlider(rng=(0,150),label="Filt:",value=0.0)
		self.gbl2.addWidget(self.wfilt,5,0,1,2)
		
		self.curbox=-1
		
		self.boxes=[]						# array of box info, each is (x,y,z,...)
		self.helixboxes=[]					# array of helix box info. each is (xi, yi, zi, xf, yf, zf)
		self.boxesimgs=[]					# z projection of each box
		self.xydown=self.xzdown=self.zydown=None
		self.firsthbclick = None

		# coordinate display
		self.wcoords=QtGui.QLabel("X: " + str(self.get_x()) + "\t\t" + "Y: " + str(self.get_y()) + "\t\t" + "Z: " + str(self.get_z()))
		self.gbl2.addWidget(self.wcoords, 0, 0, 1, 2)

		# file menu
		QtCore.QObject.connect(self.mfile_open,QtCore.SIGNAL("triggered(bool)")  ,self.menu_file_open  )
		QtCore.QObject.connect(self.mfile_read_boxloc,QtCore.SIGNAL("triggered(bool)")  ,self.menu_file_read_boxloc  )
		QtCore.QObject.connect(self.mfile_save_boxloc,QtCore.SIGNAL("triggered(bool)")  ,self.menu_file_save_boxloc  )
		QtCore.QObject.connect(self.mfile_save_boxes_stack,QtCore.SIGNAL("triggered(bool)")  ,self.save_boxes)
		QtCore.QObject.connect(self.mfile_quit,QtCore.SIGNAL("triggered(bool)")  ,self.menu_file_quit)

		# window menu
		QtCore.QObject.connect(self.mwin_boxes,QtCore.SIGNAL("triggered(bool)")  ,self.menu_win_boxes  )
		QtCore.QObject.connect(self.mwin_single,QtCore.SIGNAL("triggered(bool)")  ,self.menu_win_single  )
#		QtCore.QObject.connect(self.mwin_average,QtCore.SIGNAL("triggered(bool)")  ,self.menu_win_average  )

		# all other widgets
		QtCore.QObject.connect(self.wdepth,QtCore.SIGNAL("valueChanged(int)"),self.event_depth)
		QtCore.QObject.connect(self.wnlayers,QtCore.SIGNAL("valueChanged(int)"),self.event_nlayers)
		QtCore.QObject.connect(self.wboxsize,QtCore.SIGNAL("valueChanged"),self.event_boxsize)
		QtCore.QObject.connect(self.wmaxmean,QtCore.SIGNAL("clicked(bool)"),self.event_projmode)
		QtCore.QObject.connect(self.wscale,QtCore.SIGNAL("valueChanged")  ,self.event_scale  )
		QtCore.QObject.connect(self.wfilt,QtCore.SIGNAL("valueChanged")  ,self.event_filter  )
		QtCore.QObject.connect(self.wlocalbox,QtCore.SIGNAL("stateChanged(int)")  ,self.event_localbox  )

		QtCore.QObject.connect(self.xyview,QtCore.SIGNAL("mousemove"),self.xy_move)
		QtCore.QObject.connect(self.xyview,QtCore.SIGNAL("mousedown"),self.xy_down)
		QtCore.QObject.connect(self.xyview,QtCore.SIGNAL("mousedrag"),self.xy_drag)
		QtCore.QObject.connect(self.xyview,QtCore.SIGNAL("mouseup"),self.xy_up  )
		QtCore.QObject.connect(self.xyview,QtCore.SIGNAL("mousewheel"),self.xy_wheel  )
		QtCore.QObject.connect(self.xyview,QtCore.SIGNAL("set_scale"),self.xy_scale)
		QtCore.QObject.connect(self.xyview,QtCore.SIGNAL("origin_update"),self.xy_origin)

		QtCore.QObject.connect(self.xzview,QtCore.SIGNAL("mousedown"),self.xz_down)
		QtCore.QObject.connect(self.xzview,QtCore.SIGNAL("mousedrag"),self.xz_drag)
		QtCore.QObject.connect(self.xzview,QtCore.SIGNAL("mouseup")  ,self.xz_up  )
		QtCore.QObject.connect(self.xzview,QtCore.SIGNAL("set_scale"),self.xz_scale)
		QtCore.QObject.connect(self.xzview,QtCore.SIGNAL("origin_update"),self.xz_origin)

		QtCore.QObject.connect(self.zyview,QtCore.SIGNAL("mousedown"),self.zy_down)
		QtCore.QObject.connect(self.zyview,QtCore.SIGNAL("mousedrag"),self.zy_drag)
		QtCore.QObject.connect(self.zyview,QtCore.SIGNAL("mouseup")  ,self.zy_up  )
		QtCore.QObject.connect(self.zyview,QtCore.SIGNAL("set_scale"),self.zy_scale)
		QtCore.QObject.connect(self.zyview,QtCore.SIGNAL("origin_update"),self.zy_origin)
		
		QtCore.QObject.connect(self.xyview,QtCore.SIGNAL("keypress"),self.key_press)
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
			
		if options.inmemory:
			data=EMData(datafile)
			self.set_data(data)
		else:
			self.set_datafile(datafile)		# This triggers a lot of things to happen, so we do it last

		# Boxviewer subwidget (details of a single box)
		self.boxviewer=EMBoxViewer()
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

		QtCore.QObject.connect(self.boxesviewer,QtCore.SIGNAL("mx_image_selected"),self.img_selected)
		
		self.jsonfile=info_name(datafile)
		
		info=js_open_dict(self.jsonfile)
		self.sets={}
		self.boxsize={}
		if info.has_key("class_list"):
			clslst=info["class_list"]
			for k in sorted(clslst.keys()):
				if type(clslst[k])==dict:
					self.sets[int(k)]=str(clslst[k]["name"])
					self.boxsize[int(k)]=int(clslst[k]["boxsize"])
				else:
					self.sets[int(k)]=str(clslst[k])
					self.boxsize[int(k)]=boxsize
				
					
		
			
			
		clr=QtGui.QColor
		self.setcolors=[clr("blue"),clr("green"),clr("red"),clr("cyan"),clr("purple"),clr("orange"), clr("yellow"),clr("hotpink"),clr("gold")]
		self.sets_visible={}
		
		if info.has_key("boxes_3d"):
			box=info["boxes_3d"]
			for i,b in enumerate(box):
				#### X-center,Y-center,Z-center,method,[score,[class #]]
				bdf=[0,0,0,"manual",0.0, 0]
				for j in range(len(b)):
					bdf[j]=b[j]
				
				
				if bdf[5] not in self.sets.keys():
					clsi=int(bdf[5])
					self.sets[clsi]="particles_{:02d}".format(clsi)
				
				self.boxes.append(bdf)
				
		
		info.close()
		if len(self.sets)>0:
			self.sets_visible[self.sets.keys()[0]]=0
			self.currentset=sorted(self.sets.keys())[0]
			self.setspanel.update_sets()
		self.e = None
		print self.sets
		for i in range(len(self.boxes)):
			self.update_box(i)
		
		self.update_all()
		self.initialized=True

	def menu_win_boxes(self) : self.boxesviewer.show()
	def menu_win_single(self) : self.boxviewer.show()
#	def menu_win_average(self) : self.averageviewer.show()

	def set_datafile(self,datafile):
		print "\nIn set_datafile, received datafile", datafile
		if datafile==None :
			self.datafile=None
			self.data=None
			self.xyview.set_data(None)
			self.xzview.set_data(None)
			self.zyview.set_data(None)
			return

		self.data=None
		self.datafile=datafile

		#print "\nDatafile set, see!", self.datafile, type(self.datafile)

		imgh=EMData(datafile,0,1)

		if self.yshort:
			self.datasize=(imgh["nx"],imgh["nz"],imgh["ny"])
		else:
			self.datasize=(imgh["nx"],imgh["ny"],imgh["nz"])

		self.wdepth.setRange(0,self.datasize[2]-1)
		self.boxes=[]
		self.curbox=-1

		self.wdepth.setValue(self.datasize[2]/2)
		self.update_all()

	def set_data(self,data):
		if data==None :
			self.datafile=None
			self.data=None
			self.xyview.set_data(None)
			self.xzview.set_data(None)
			self.zyview.set_data(None)
			return

		self.data=data
		self.datafile=None

		if self.yshort:
			self.datasize=(data["nx"],data["nz"],data["ny"])
		else:
			self.datasize=(data["nx"],data["ny"],data["nz"])

		self.wdepth.setRange(0,self.datasize[2]-1)
		self.boxes=[]
		self.curbox=-1

		self.wdepth.setValue(self.datasize[2]/2)
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
		if self.yshort:
			if self.data!=None:
				r=self.data.get_clip(Region(x-bs/2,z-bz/2,y-bs/2,bs,bz,bs))
				if options.normproc:
					r.process_inplace(options.normproc)
				r.process_inplace("xform",{"transform":Transform({"type":"eman","alt":90.0})})
				r.process_inplace("xform.mirror",{"axis":"z"})
			elif self.datafile!=None:
				r=EMData(self.datafile,0,0,Region(x-bs/2,z-bz/2,y-bs/2,bs,bz,bs))
				if options.normproc:
					r.process_inplace(options.normproc)
				r.process_inplace("xform",{"transform":Transform({"type":"eman","alt":90.0})})
				r.process_inplace("xform.mirror",{"axis":"z"})
			else: return None

		else :
			if self.data!=None:
				r=self.data.get_clip(Region(x-bs/2,y-bs/2,z-bz/2,bs,bs,bz))
			elif self.datafile!=None:
				r=EMData(self.datafile,0,0,Region(x-bs/2,y-bs/2,z-bz/2,bs,bs,bz))
			else: return None

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
		if self.yshort:
			if self.data!=None :
				if xyz==0:
					r=self.data.get_clip(Region(n,0,0,1,self.datasize[2],self.datasize[1]))
					r.set_size(self.datasize[2],self.datasize[1],1)
				elif xyz==2:
					r=self.data.get_clip(Region(0,n,0,self.datasize[0],1,self.datasize[1]))
					r.set_size(self.datasize[0],self.datasize[1],1)
				else:
					r=self.data.get_clip(Region(0,0,n,self.datasize[0],self.datasize[2],1))

			elif self.datafile!=None:
				if xyz==0:
					r=EMData()
					r.read_image(self.datafile,0,0,Region(n,0,0,1,self.datasize[2],self.datasize[1]))
					r.set_size(self.datasize[2],self.datasize[1],1)

				elif xyz==2:
					r=EMData()
					r.read_image(self.datafile,0,0,Region(0,n,0,self.datasize[0],1,self.datasize[1]))
					r.set_size(self.datasize[0],self.datasize[1],1)
				else:
					r=EMData()
					r.read_image(self.datafile,0,0,Region(0,0,n,self.datasize[0],self.datasize[2],1))
			else:
				return None

		else :
			if self.data!=None :
				if xyz==0:
					r=self.data.get_clip(Region(n,0,0,1,self.datasize[1],self.datasize[2]))
					r.set_size(self.datasize[1],self.datasize[2],1)
				elif xyz==1:
					r=self.data.get_clip(Region(0,n,0,self.datasize[0],1,self.datasize[2]))
					r.set_size(self.datasize[0],self.datasize[2],1)
				else:
					r=self.data.get_clip(Region(0,0,n,self.datasize[0],self.datasize[1],1))

			elif self.datafile!=None:
				if xyz==0:
					r=EMData()
					r.read_image(self.datafile,0,0,Region(n,0,0,1,self.datasize[1],self.datasize[2]))
					r.set_size(self.datasize[1],self.datasize[2],1)
				elif xyz==1:
					r=EMData()
					r.read_image(self.datafile,0,0,Region(0,n,0,self.datasize[0],1,self.datasize[2]))
					r.set_size(self.datasize[0],self.datasize[2],1)
				else:
					r=EMData()
					r.read_image(self.datafile,0,0,Region(0,0,n,self.datasize[0],self.datasize[1],1))

			else :
				return None

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
		for i in range(len(self.boxes)):
			if self.boxes[i][5]==self.currentset:
				self.update_box(i)
		self.update_box(cb)

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
				print "No box size saved for {}..".format(clsid)
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
		QtGui.QMessageBox.warning(None,"Error","Sorry, in the current version, you must provide a file to open on the command-line.")

	def load_box_yshort(self, boxcoords):
		if options.yshort:
			return [boxcoords[0], boxcoords[2], boxcoords[1]]
		else:
			return boxcoords

	def menu_file_read_boxloc(self):
		fsp=str(QtGui.QFileDialog.getOpenFileName(self, "Select output text file"))

		f=file(fsp,"r")
		if self.helixboxer:
			for b in f:
				b2=[int(float(i))/self.shrink for i in b.split()[:6]]
				self.boxes.append(self.load_box_yshort(b2[3:6]))
				self.update_box(len(self.boxes)-1)
				self.helixboxes.append(b2)
				self.update_helixbox(len(self.helixboxes)-1)
				self.boxes.append(self.load_box_yshort(b2[0:3]))
				self.update_box(len(self.boxes)-1)
		else:
			for b in f:
				b2=[int(float(i))/self.shrink for i in b.split()[:3]]
				bdf=[0,0,0,"manual",0.0, 0]
				for j in range(len(b2)):
					bdf[j]=b2[j]
				self.boxes.append(bdf)
				self.update_box(len(self.boxes)-1)
		f.close()

	def menu_file_save_boxloc(self):
		shrinkf=self.shrink 								#jesus

		fsp=str(QtGui.QFileDialog.getSaveFileName(self, "Select output text file"))

		out=file(fsp,"w")
		if self.helixboxer:
			for b in self.helixboxes:
				out.write("%d\t%d\t%d\t%d\t%d\t%d\n"%(b[0]*shrinkf,b[1]*shrinkf,b[2]*shrinkf,b[3]*shrinkf,b[4]*shrinkf,b[5]*shrinkf))
		else:
			for b in self.boxes:
				out.write("%d\t%d\t%d\n"%(b[0]*shrinkf,b[1]*shrinkf,b[2]*shrinkf))
		out.close()


	def save_boxes(self, clsid=[]):
		if len(clsid)==0:
			defaultname="ptcls.hdf"
		else:
			defaultname="_".join([self.sets[i] for i in clsid])+".hdf"
		
		name,ok=QtGui.QInputDialog.getText( self, "Save particles", "Filename suffix:", text=defaultname)
		if not ok:
			return
		name=self.filetag+str(name)
		if name[-4:].lower()!=".hdf" :
			name+=".hdf"
			
		
		for dr in ["particles3d", "particles"]:
			if not os.path.isdir(dr):
				os.mkdir(dr)
		
		fsp=os.path.join("particles3d",self.basename)+name
		fspprjs=os.path.join("particles",self.basename)+name.replace('.hdf','_prjs.hdf')
		print "Saving 3D particles to {},\n Saving particle projections to {}".format(fsp, fspprjs)
		for f in [fsp, fspprjs]:
			if os.path.isfile(f):
				print "{} exist. Overwritting...".format(f)
				os.remove(f)
		
		progress = QtGui.QProgressDialog("Saving", "Abort", 0, len(self.boxes),None)
		
		
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
					print "Inconsistant box size in the particles to save.. Using {:d}..".format(boxsz)
					bs=boxsz
			
			sz=[s/2 for s in self.datasize]
			img=self.get_cube(b[0], b[1], b[2], boxsz=bs)
			img.process_inplace('normalize')
			
			img["ptcl_source_image"]=self.datafilename
			img["ptcl_source_coord"]=(b[0]-sz[0], b[1]-sz[1], b[2]-sz[2])
			
			if self.invert:
				img.mult(-1)
			prj=img.project("standard", Transform())
			
			prj["ptcl_source_image"]=self.datafilename
			prj["ptcl_source_coord"]=(b[0]-sz[0], b[1]-sz[1], b[2]-sz[2])
			
			img.write_image(fsp,-1)
			prj.write_image(fspprjs,-1)

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
		if self.wmaxmean.isChecked() : return Averagers.get("minmax",{"max":1})

		return Averagers.get("mean")

	def update_sides(self):
		"""updates xz and yz views due to a new center location"""

		#print "\n\n\n\n\nIn update sides, self.datafile is", self.datafile
		#print "\n\n\n\n"

		if self.datafile==None and self.data==None:
			return

		if self.curbox==-1 :
			x=self.datasize[0]/2
			y=self.datasize[1]/2
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
				if self.boxes[i][1]<self.cury+bs/2 and self.boxes[i][1]>self.cury-bs/2 and  self.boxes[i][5] in self.sets_visible:
					xzs[i][0]="rect"
				else:
					xzs[i][0]="hidden"

			zys=self.zyview.get_shapes()
			
			for i in range(len(self.boxes)):
				bs=self.get_boxsize(self.boxes[i][5])
				if self.boxes[i][0]<self.curx+bs/2 and self.boxes[i][0]>self.curx-bs/2 and  self.boxes[i][5] in self.sets_visible:
					zys[i][0]="rect"
				else:
					zys[i][0]="hidden"
		else :
			xzs=self.xzview.get_shapes()
			zys=self.zyview.get_shapes()
		
			for i in range(len(self.boxes)):
				bs=self.get_boxsize(self.boxes[i][5])
				if  self.boxes[i][5] in self.sets_visible:
					xzs[i][0]="rect"
					zys[i][0]="rect"
				else:
					xzs[i][0]="hidden"
					zys[i][0]="hidden"

		self.xzview.shapechange=1
		self.zyview.shapechange=1

		# yz
		avgr=self.get_averager()

		for x in range(x-self.nlayers()/2,x+(self.nlayers()+1)/2):
			slc=self.get_slice(x,0)
			avgr.add_image(slc)

		av=avgr.finish()
		if not self.yshort:
			av.process_inplace("xform.transpose")

		if self.wfilt.getValue()!=0.0:
			av.process_inplace("filter.lowpass.gauss",{"cutoff_freq":1.0/self.wfilt.getValue(),"apix":self.apix})

		self.zyview.set_data(av)

		# xz
		avgr=self.get_averager()

		for y in range(y-self.nlayers()/2,y+(self.nlayers()+1)/2):
			slc=self.get_slice(y,1)
			avgr.add_image(slc)

		av=avgr.finish()
		if self.wfilt.getValue()!=0.0:
			av.process_inplace("filter.lowpass.gauss",{"cutoff_freq":1.0/self.wfilt.getValue(),"apix":self.apix})

		self.xzview.set_data(av)


	def update_xy(self):
		"""updates xy view due to a new slice range"""

		#print "\n\n\n\n\nIn update_xy, self.datafile is", self.datafile
		#print "\n\n\n\n"

		if self.datafile==None and self.data==None:
			return



		# Boxes should also be limited by default in the XY view
		if len(self.boxes) > 0:
			zc=self.wdepth.value()
			#print "The current depth is", self.wdepth.value()
			xys=self.xyview.get_shapes()
			for i in range(len(self.boxes)):
				
				bs=self.get_boxsize(self.boxes[i][5])
				#print "the z coord of box %d is %d" %(i,self.boxes[i][2])
				#print "therefore the criteria to determine whether to display it is", abs(self.boxes[i][2] - zc)
				zdist=abs(self.boxes[i][2] - zc)
				if zdist < bs/2 and self.boxes[i][5] in self.sets_visible:
					#print "Which is less than half the box thus it survives"
					xys[i][0]="circle"
					xys[i][6]=bs/2-zdist
				else :
					xys[i][0]="hidden"
					#print "Which is more than half the box and thus it dies"

			self.xyview.shapechange=1

		if self.wmaxmean.isChecked():
			avgr=Averagers.get("minmax",{"max":1})

		else:
			avgr=Averagers.get("mean")

		slc=EMData()
		for z in range(self.wdepth.value()-self.nlayers()/2,self.wdepth.value()+(self.nlayers()+1)/2):
			slc=self.get_slice(z,2)
			avgr.add_image(slc)

		av=avgr.finish()

		#print "\n\nIn update xy, av and type are", av, type(av)

		if self.wfilt.getValue()!=0.0:

			av.process_inplace("filter.lowpass.gauss",{"cutoff_freq":1.0/self.wfilt.getValue(),"apix":self.apix})
		self.xyview.set_data(av)

	def update_all(self):
		"""redisplay of all widgets"""

		#print "\n\n\n\n\nIn update all, self.datafile is", self.datafile
		#print "\n\n\n\n"
		if self.datafile==None and self.data==None:
			return



		self.update_xy()
		self.update_sides()
		self.update_boximgs()

		#self.xyview.update()
		#self.xzview.update()
		#self.zyview.update()

	def update_coords(self):
		self.wcoords.setText("X: " + str(self.get_x()) + "\t\t" + "Y: " + str(self.get_y()) + "\t\t" + "Z: " + str(self.get_z()))

	def inside_box(self,n,x=-1,y=-1,z=-1):
		"""Checks to see if a point in image coordinates is inside box number n. If any value is negative, it will not be checked."""
		box=self.boxes[n]
		if box[5] not in self.sets_visible:
			return False
		bs=self.get_boxsize(box[5])/2
		rr=(x>=0)*((box[0]-x)**2) + (y>=0)*((box[1]-y) **2) + (z>=0)*((box[2]-z)**2)
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

	def do_helix_deletion(self, n):
		if n==len(self.helixboxes)-1 :
			self.helixboxes.pop()
			self.xyview.del_shape(str(n)+"helix")
			self.xzview.del_shape(str(n)+"helix")
			self.zyview.del_shape(str(n)+"helix")
		else:
			a=self.helixboxes.pop()
			self.helixboxes[n]=a
			self.xyview.del_shape(str(len(self.helixboxes))+"helix")
			self.xzview.del_shape(str(len(self.helixboxes))+"helix")
			self.zyview.del_shape(str(len(self.helixboxes))+"helix")
			self.update_helixbox(n)

	def del_box(self,n):
		"""Delete an existing box by replacing the deleted box with the last box. A bit funny, but otherwise
		update after deletion is REALLY slow."""
#		print "del ",n
		if n<0 or n>=len(self.boxes): return

		if self.boxviewer.get_data(): self.boxviewer.set_data(None)
		self.curbox=-1
		self.do_deletion([n])

	def compute_crossAB(self, a, b):
		c1 = a[1]*b[2] - a[2]*b[1]
		c2 = a[2]*b[0] - a[0]*b[2]
		c3 = a[0]*b[1] - a[1]*b[0]
		return Vec3f(c1,c2,c3)

	def compute_perpZ(self, a):
		# Z axis
		b1 = -a[1]
		b2 = a[0]
		b3 = 0
		return Vec3f(b1,b2,b3)

	def compute_perpY(self, a):
		# Y axis
		b1 = -a[2]
		b2 = 0
		b3 = a[0]

		return Vec3f(b1,b2,b3)

	def get_box_coord_system(self, helixbox):
		"""
		Compute the coordinate system for the box
		"""
		a = Vec3f((helixbox[0]-helixbox[3]), (helixbox[1]-helixbox[4]), (helixbox[2]-helixbox[5]))

		a.normalize()
		b = self.compute_perpZ(a)
		b.normalize()
		c = self.compute_crossAB(a, b)

		return Transform([a[0],a[1],a[2],helixbox[3],b[0],b[1],b[2],helixbox[4],c[0],c[1],c[2],helixbox[5]])


	def get_extended_a_vector(self, helixbox):
		"""
		Extend the A vector to the box ends
		"""
		a = Vec3f((helixbox[3]-helixbox[0]), (helixbox[4]-helixbox[1]), (helixbox[5]-helixbox[2]))
		a.normalize()
		bs = self.get_boxsize()
		return [(helixbox[0] - a[0]*bs/2),(helixbox[1] - a[1]*bs/2),(helixbox[2] - a[2]*bs/2),(helixbox[3] + a[0]*bs/2),(helixbox[4] + a[1]*bs/2),(helixbox[5] + a[2]*bs/2)]



	def update_box(self,n,quiet=False):
		"""After adjusting a box, call this"""
#		print "upd ",n,quiet

		try:
			box=self.boxes[n]
		except IndexError:
			return
		bs2=self.get_boxsize(box[5])/2

		#if self.curbox!=n :
			#self.xzview.scroll_to(None,box[2])
			#self.zyview.scroll_to(box[2],None)


		# Boxes may not extend outside the tomogram
		if box[0]<bs2 : box[0]=bs2
		if box[0]>self.datasize[0]-bs2 : box[0]=self.datasize[0]-bs2
		if box[1]<bs2 : box[1]=bs2
		if box[1]>self.datasize[1]-bs2 : box[1]=self.datasize[1]-bs2
		if box[2]<bs2 : box[2]=bs2
		if box[2]>self.datasize[2]-bs2 : box[2]=self.datasize[2]-bs2
#		print self.boxes

		
		
		color=self.setcolors[box[5]%len(self.setcolors)].getRgbF()
		
		#self.xyview.add_shape(n,EMShape(("rect",.2,.2,.8,box[0]-bs2,box[1]-bs2,box[0]+bs2,box[1]+bs2,2)))
		self.xyview.add_shape(n,EMShape(["circle",color[0],color[1],color[2],box[0],box[1],bs2,2]))
		#self.xyview.add_shape("xl",EMShape(("line",.8,.8,.1,0,box[1],self.datasize[0],box[1],1)))
		#self.xyview.add_shape("yl",EMShape(("line",.8,.8,.1,box[0],0,box[0],self.datasize[1],1)))
		#self.xzview.add_shape(n,EMShape(["circle",.2,.2,.8,box[0],box[2],bs2,2]))
		self.xzview.add_shape(n,EMShape(("rect",color[0],color[1],color[2],box[0]-bs2,box[2]-bs2,box[0]+bs2,box[2]+bs2,2)))
		#self.xzview.add_shape("xl",EMShape(("line",.8,.8,.1,0,box[2],self.datasize[0],box[2],1)))
		#self.xzview.add_shape("zl",EMShape(("line",.8,.8,.1,box[0],0,box[0],self.datasize[2],1)))
		#self.zyview.add_shape(n,EMShape(["circle",.2,.2,.8,box[2],box[1],bs2,2]))
		self.zyview.add_shape(n,EMShape(("rect",color[0],color[1],color[2],box[2]-bs2,box[1]-bs2,box[2]+bs2,box[1]+bs2,2)))
		#self.zyview.add_shape("yl",EMShape(("line",.8,.8,.1,box[2],0,box[2],self.datasize[1],1)))
		#self.zyview.add_shape("zl",EMShape(("line",.8,.8,.1,0,box[1],self.datasize[2],box[1],1)))

		if self.depth()!=box[2]:
			self.wdepth.setValue(box[2])
		else:
			self.xyview.update()
		if self.initialized: self.update_sides()

		# For speed, we turn off updates while dragging a box around. Quiet is set until the mouse-up
		if not quiet and not self.helixboxer:
			# Get the cube from the original data (normalized)
			cube=self.get_cube(box[0], box[1], box[2], centerslice=True, boxsz=self.get_boxsize(box[5]))
			#self.boxviewer.set_data(cube)

			# Make a z projection and store it in the list of all boxes
			#proj=cube.process("misc.directional_sum",{"axis":"z"})
			proj=cube
			proj.process_inplace("normalize")
			
			for i in range(len(self.boxesimgs),n+1): 
				self.boxesimgs.append(None)
			
			self.boxesimgs[n]=proj
			#try: 
				#self.boxesimgs[n]=proj
			#except:
				#for i in range(len(self.boxesimgs),n+1): 
					#self.boxesimgs.append(None)
				#self.boxesimgs[n]=proj
			mm=[m for im,m in enumerate(self.boxesimgs) if self.boxes[im][5] in self.sets_visible]
			
		if self.initialized:
			self.update_boximgs()

			if n!=self.curbox and not self.helixboxer:
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

	def add_helix_box(self, xf, yf, zf, xi, yi, zi):
		print xf, yf, zf, xi, yi, zi
		if options.yshort:
			self.helixboxes.append([xf, zf, yf, xi, zi, yi])
		else:
			self.helixboxes.append([xf, yf, zf, xi, yi, zi])
	
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
					if self.del_box(i) != "DELHELIX": self.firsthbclick = None
				else:
					self.xydown=(i,x,y,self.boxes[i][0],self.boxes[i][1])
					if self.helixboxer: self.update_helixbox(int(i/2))
					self.update_box(i)
				break
		else:
#			if x>self.get_boxsize()/2 and x<self.datasize[0]-self.get_boxsize()/2 and y>self.get_boxsize()/2 and y<self.datasize[1]-self.get_boxsize()/2 and self.depth()>self.get_boxsize()/2 and self.depth()<self.datasize[2]-self.get_boxsize()/2 :
			if not event.modifiers()&Qt.ShiftModifier:
				###########
				if self.helixboxer:	# Only create a helixbox every 2 clicks
					if self.firsthbclick:
						self.add_helix_box(x, y, self.depth(), self.firsthbclick[0], self.firsthbclick[1], self.firsthbclick[2])
						self.firsthbclick = None
						self.update_helixbox(len(self.helixboxes)-1)
					else:
						self.firsthbclick = [x, y, self.depth()]
				###########
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
		if event.delta() > 0:
			#self.wdepth.setValue(self.wdepth.value()+4)
			self.wdepth.setValue(self.wdepth.value()+1) #jesus

		elif event.delta() < 0:
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
					if self.del_box(i) != "DELHELIX": self.firsthbclick = None
				else :
					self.xzdown=(i,x,z,self.boxes[i][0],self.boxes[i][2])
					if self.helixboxer: self.update_helixbox(int(i/2))
					self.update_box(i)
				break
		else:
			if not event.modifiers()&Qt.ShiftModifier:
				###########
				if self.helixboxer:	# Only create a helixbox every 2 clicks
					if self.firsthbclick:
						self.add_helix_box(x, self.cury, z, self.firsthbclick[0], self.firsthbclick[1], self.firsthbclick[2])
						self.firsthbclick = None
						self.update_helixbox(len(self.helixboxes)-1)
					else:
						self.firsthbclick = [x, self.cury, z]
				###########
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
		if self.helixboxer:
			if len(self.boxes) % 2 == 0 or (self.xzdown[0] != len(self.boxes)-1):	# Only update the helix boxer if it is paired, otherwise treat it as a regular box
				hb = self.helixboxes[int(self.xzdown[0]/2)]
				if self.xzdown[0] % 2 == 0:
					hb[3] = dx+self.xzdown[3]
					hb[5] = dz+self.xzdown[4]
				else:
					hb[0] = dx+self.xzdown[3]
					hb[2] = dz+self.xzdown[4]
				self.update_helixbox(int(self.xzdown[0]/2))
			else:
				self.firsthbclick[0] = x
				self.firsthbclick[2] = z

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
					if self.del_box(i) != "DELHELIX": self.firsthbclick = None
				else :
					self.zydown=(i,z,y,self.boxes[i][2],self.boxes[i][1])
					if self.helixboxer: self.update_helixbox(int(i/2))
					self.update_box(i)
				break
		else:
			if not event.modifiers()&Qt.ShiftModifier:
				###########
				if self.helixboxer:	# Only create a helixbox every 2 clicks
					if self.firsthbclick:
						self.add_helix_box(self.curx, y, z, self.firsthbclick[0], self.firsthbclick[1], self.firsthbclick[2])
						self.firsthbclick = None
						self.update_helixbox(len(self.helixboxes)-1)
					else:
						self.firsthbclick = [self.curx, y, z]
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
		self.boxsize[i]=32
		
		return
	
	def save_set(self):
		
		self.save_boxes(self.sets_visible.keys())
		return
	
	
	def key_press(self,event):
		if event.key() == 96:
			self.wdepth.setValue(self.wdepth.value()+1)

		elif event.key() == 49:
			self.wdepth.setValue(self.wdepth.value()-1)
		else:
			self.emit(QtCore.SIGNAL("keypress"), event)

	
	def closeEvent(self,event):
		print "Exiting"
		info=js_open_dict(self.jsonfile)
		info["boxes_3d"]=self.boxes
		clslst={}
		for key in self.sets.keys():
			clslst[int(key)]={
				"name":self.sets[key],
				"boxsize":self.boxsize[key],
				}
		info["class_list"]=clslst
		info.close()
		
		self.boxviewer.close()
		self.boxesviewer.close()
		self.optionviewer.close()
		self.xyview.close()
		self.xzview.close()
		self.zyview.close()
		
#		self.averageviewer.close()
		#event.accept()
		#self.app().close_specific(self)
		self.emit(QtCore.SIGNAL("module_closed")) # this signal is important when e2ctf is being used by a program running its own event loop

	#def closeEvent(self,event):
		#self.target().done()
		

def parse_setname(name):
	p0=name.find('::')
	ret=-1
	if p0>0:
		try:
			ret=int(name[:p0])
		except:
			pass
	
	return ret
			
	

class EMBoxViewer(QtGui.QWidget):
	"""This is a multi-paned view showing a single boxed out particle from a larger tomogram"""

	def __init__(self):
		QtGui.QWidget.__init__(self)
		self.setWindowTitle("Single Particle View")

		self.resize(300,300)

		self.gbl = QtGui.QGridLayout(self)
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

		QtCore.QObject.connect(self.wfilt,QtCore.SIGNAL("valueChanged")  ,self.event_filter  )

		self.gbl.setRowStretch(2,1)
		self.gbl.setRowStretch(0,5)
		self.gbl.setRowStretch(1,5)
		#QtCore.QObject.connect(self.xyview,QtCore.SIGNAL("mousedown"),self.xy_down)
		#QtCore.QObject.connect(self.xyview,QtCore.SIGNAL("mousedrag"),self.xy_drag)
		#QtCore.QObject.connect(self.xyview,QtCore.SIGNAL("mouseup")  ,self.xy_up  )

		#QtCore.QObject.connect(self.xzview,QtCore.SIGNAL("mousedown"),self.xz_down)
		#QtCore.QObject.connect(self.xzview,QtCore.SIGNAL("mousedrag"),self.xz_drag)
		#QtCore.QObject.connect(self.xzview,QtCore.SIGNAL("mouseup")  ,self.xz_up  )

		#QtCore.QObject.connect(self.zyview,QtCore.SIGNAL("mousedown"),self.zy_down)
		#QtCore.QObject.connect(self.zyview,QtCore.SIGNAL("mousedrag"),self.zy_drag)
		#QtCore.QObject.connect(self.zyview,QtCore.SIGNAL("mouseup")  ,self.zy_up  )

#		self.setSizeGripEnabled(True)

#		if get_platform() == "Darwin": # because OpenGL widgets in Qt don't leave room in the bottom right hand corner for the resize tool
#			self.status = QtGui.QStatusBar()
#			self.gbl.addWidget(self.status,3,0,1,2)
#			self.margin = 0

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
			self.fdata=self.data.process("filter.lowpass.gauss",{"cutoff_freq":1.0/self.wfilt.getValue(),"apix":self.data['apix_x']}) #JESUS

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


class EMTomoBoxerOptions(QtGui.QWidget):
	def __init__(self,target) :
		QtGui.QWidget.__init__(self)
		#print "aaaaaaaa"
		self.setWindowTitle("Options")
		self.target=weakref.ref(target)
		
		self.gbl = QtGui.QGridLayout(self)
		#self.gbl.setMargin(2)
		#self.gbl.setSpacing(6)
		self.gbl.setObjectName("gbl")
		
		
		self.erasercheckbox=QtGui.QCheckBox("Eraser")
		self.gbl.addWidget(self.erasercheckbox,0,0)
		
		self.eraser_radius=ValBox(label="Radius:",value=64)
		self.gbl.addWidget(self.eraser_radius,0,1)

		self.tabwidget = QtGui.QTabWidget()
		self.gbl.addWidget(self.tabwidget,1,0,1,2)
		
	def add_panel(self,widget,name):
		self.tabwidget.addTab(widget,name)

		return 

#### Copied from emimagemx.py since some modification are needed...
class EMTomoSetsPanel(QtGui.QWidget):
	'''
	This is the set display panel
	'''
	def __init__(self,target):
		QtGui.QWidget.__init__(self)

		self.target = weakref.ref(target) # this should be the EMImageMXWidget
		self.busy = False
		self.initialized=False

		# cached values for speed later
		self.itemflags=	Qt.ItemFlags(Qt.ItemIsEditable)|Qt.ItemFlags(Qt.ItemIsSelectable)|Qt.ItemFlags(Qt.ItemIsEnabled)|Qt.ItemFlags(Qt.ItemIsUserCheckable)

		# now build the interface
		hbl = QtGui.QHBoxLayout(self)
		self.setlist=QtGui.QListWidget()
		self.setlist.setSizePolicy(QtGui.QSizePolicy.Preferred,QtGui.QSizePolicy.Expanding)
		hbl.addWidget(self.setlist)

		vbl = QtGui.QVBoxLayout()

		self.new_set_button = QtGui.QPushButton("New")
		vbl.addWidget(self.new_set_button)
		self.rename_set_button = QtGui.QPushButton("Rename")
		vbl.addWidget(self.rename_set_button)
		self.save_set_button = QtGui.QPushButton("Save")
		vbl.addWidget(self.save_set_button)
		self.delete_set_button = QtGui.QPushButton("Delete")
		vbl.addWidget(self.delete_set_button)

		hbl.addLayout(vbl)

		QtCore.QObject.connect(self.save_set_button, QtCore.SIGNAL("clicked(bool)"), self.save_set)
		QtCore.QObject.connect(self.new_set_button, QtCore.SIGNAL("clicked(bool)"), self.new_set)
		QtCore.QObject.connect(self.rename_set_button, QtCore.SIGNAL("clicked(bool)"), self.rename_set)
		QtCore.QObject.connect(self.delete_set_button, QtCore.SIGNAL("clicked(bool)"), self.delete_set)
		QtCore.QObject.connect(self.setlist,QtCore.SIGNAL("itemChanged(QListWidgetItem*)"),self.set_list_item_changed)
		QtCore.QObject.connect(self.setlist,QtCore.SIGNAL("currentRowChanged(int)"),self.set_list_row_changed)


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
		if not self.initialized: return 
		a = self.setlist.item(i)
		if a==None : return
		name = str(a.text())
		self.target().set_current_set(name)
		self.update_sets()

	def set_list_item_changed(self,item):
		name=str(item.text())
		if item.checkState() == Qt.Checked : self.target().show_set(name)
		else: self.target().hide_set(name)
		

	def delete_set(self,unused):
		selections = self.setlist.selectedItems()
		names=[str(i.text()) for i in selections]
		cancel=QtGui.QMessageBox.warning(self, "Delete set", "Are you sure to delete {}? This will remove all particles in that class".format(names[0]), "Yes", "No")
		if not cancel:
			self.target().delete_set(names[0])
		self.update_sets()


	def new_set(self,unused=None):
		name,ok=QtGui.QInputDialog.getText( self, "Set Name", "Enter a name for the new set:")
		if not ok : return
		name=str(name)
		if name in self.target().sets :
			print "Set name exists"
			return

		self.target().new_set(name)
		self.update_sets()

	def rename_set(self,unused=None):
		selections = self.setlist.selectedItems()
		sels=[str(i.text()) for i in selections]
		if len(sels)==0:
			return
		name,ok=QtGui.QInputDialog.getText( self, "Set Name", "Enter a name for the new set:")
		if not ok : return
		name=str(name)
		
		if name in self.target().sets :
			print "Set name exists"
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
			item=QtGui.QListWidgetItem(kname)
			item.setFlags(self.itemflags)
			item.setTextColor(self.target().setcolors[i%len(self.target().setcolors)])
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
	