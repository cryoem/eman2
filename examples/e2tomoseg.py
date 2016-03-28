#!/usr/bin/env python

#
# Author: James Michael Bell 5/20/2014 (jmbell@bcm.edu)
# Copyright (c) 2000-2006 Baylor College of Medicine
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
from emapplication import EMApp
from emdataitem3d import EMDataItem3D, EMIsosurface
from emimage2d import EMImage2DWidget
from emimagemx import EMImageMXWidget
from emscene3d import EMScene3D, EMInspector3D, EMQTreeWidget
import os
from valslider import ValSlider, EMANToolButton, EMSpinWidget, EMQTColorWidget
import weakref

from PyQt4 import QtCore
from PyQt4 import QtGui


def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options] <image file> ...
	
	e2tomoseg.py is a simple tomogram segmentation tool with built-in automation.
	"""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	#parser.add_argument("-t","--tomo",type=str,help="<rawptcl>,<classmx> Show particles associated class-averages")
	#parser.add_argument("-s","--seg",type=str,help="A specialized flag that disables auto contrast for the display of particles stacks and 2D images only.")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-2)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	
	(options, args) = parser.parse_args()
	
	logid=E2init(sys.argv)
	
	app = EMApp()
	
	volume = TomoSegVolumeViewer() # window to see 3D segmentation overlayed on top of current slice
	slices = TomoSegSliceViewer() # window to view current slice of tomogram
	tools = TomoSegInspector() # window to hold tools/annotation tree for segmentation
	
	volume.show() 
	slices.show()
	tools.show()
	
	app.exec_()

	E2end(logid)

class TomoSegVolumeViewer(EMScene3D):
	
	def __init__(self):
		super(TomoSegVolumeViewer,self).__init__()
		
		self.setWindowTitle("TomoSeg Subvolume Viewer")

class TomoSegSliceViewer(QtGui.QMainWindow):

	def __init__(self,data=None,datafile=None,yshort=False,apix=0.0,boxsize=32,shrink=1,contrast=None,center=None,mod=False,normalize=False):
		QtGui.QWidget.__init__(self)
		
		self.yshort=yshort
		self.apix=apix

		self.shrink=shrink
		self.contrast=contrast
		self.mod=mod
		self.center=center
		self.normalize=normalize
		self.setWindowTitle("Main Window")

		self.setWindowTitle("TomoSeg Slice Viewer")

		# Menu Bar
		self.mfile=self.menuBar().addMenu("File")
		self.mfile_open=self.mfile.addAction("Open")
		self.mfile_read_boxloc=self.mfile.addAction("Read Box Coord")
		self.mfile_save_boxloc=self.mfile.addAction("Save Box Coord")
		self.mfile_save_boxes=self.mfile.addAction("Save Boxed Data")
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
		#self.wboxsize=ValBox(label="Box Size:",value=boxsize)
		#self.gbl2.addWidget(self.wboxsize,0,0,1,2)
		self.oldboxsize=boxsize

		# max or mean
		self.wmaxmean=QtGui.QPushButton("MaxProj")
		self.wmaxmean.setCheckable(True)
		self.gbl2.addWidget(self.wmaxmean,1,0)

		# number slices
		self.wnlayers=QtGui.QSpinBox()
		self.wnlayers.setMinimum(1)
		self.wnlayers.setMaximum(256)
		self.wnlayers.setValue(1)
		self.gbl2.addWidget(self.wnlayers,1,1)

		# Local boxes in side view
		self.wlocalbox=QtGui.QCheckBox("Limit Side Boxes")
		self.gbl2.addWidget(self.wlocalbox,2,0)

		# scale factor
		self.wscale=ValSlider(rng=(.1,2),label="Sca:",value=1.0)
		self.gbl2.addWidget(self.wscale,3,0,1,2)

		# 2-D filters
		self.wfilt = ValSlider(rng=(0,50),label="Filt:",value=0.0)
		self.gbl2.addWidget(self.wfilt,4,0,1,2)

		self.curbox=-1
		self.boxes=[]						# array of box info, each is (x,y,z,...)
		self.helixboxes=[]					# array of helix box info. each is (xi, yi, zi, xf, yf, zf)
		self.boxesimgs=[]					# z projection of each box
		self.xydown=None
		self.firsthbclick = None

		# file menu
		QtCore.QObject.connect(self.mfile_open,QtCore.SIGNAL("triggered(bool)")  ,self.menu_file_open  )
		QtCore.QObject.connect(self.mfile_read_boxloc,QtCore.SIGNAL("triggered(bool)")  ,self.menu_file_read_boxloc  )
		QtCore.QObject.connect(self.mfile_save_boxloc,QtCore.SIGNAL("triggered(bool)")  ,self.menu_file_save_boxloc  )
		QtCore.QObject.connect(self.mfile_save_boxes,QtCore.SIGNAL("triggered(bool)")  ,self.menu_file_save_boxes  )
		QtCore.QObject.connect(self.mfile_save_boxes_stack,QtCore.SIGNAL("triggered(bool)")  ,self.menu_file_save_boxes_stack)
		QtCore.QObject.connect(self.mfile_quit,QtCore.SIGNAL("triggered(bool)")  ,self.menu_file_quit)

		# window menu
		QtCore.QObject.connect(self.mwin_boxes,QtCore.SIGNAL("triggered(bool)")  ,self.menu_win_boxes  )
		QtCore.QObject.connect(self.mwin_single,QtCore.SIGNAL("triggered(bool)")  ,self.menu_win_single  )
#		QtCore.QObject.connect(self.mwin_average,QtCore.SIGNAL("triggered(bool)")  ,self.menu_win_average  )

		# all other widgets
		QtCore.QObject.connect(self.wdepth,QtCore.SIGNAL("valueChanged(int)"),self.event_depth)
		QtCore.QObject.connect(self.wnlayers,QtCore.SIGNAL("valueChanged(int)"),self.event_nlayers)
		#QtCore.QObject.connect(self.wboxsize,QtCore.SIGNAL("valueChanged"),self.event_boxsize)
		QtCore.QObject.connect(self.wmaxmean,QtCore.SIGNAL("clicked(bool)"),self.event_projmode)
		QtCore.QObject.connect(self.wscale,QtCore.SIGNAL("valueChanged")  ,self.event_scale  )
		QtCore.QObject.connect(self.wfilt,QtCore.SIGNAL("valueChanged")  ,self.event_filter  )
		QtCore.QObject.connect(self.wlocalbox,QtCore.SIGNAL("stateChanged(int)")  ,self.event_localbox  )

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

		if datafile!=None:
			print "\nIn ETomoBoxer, datafile is", datafile
			self.set_datafile(datafile)		# This triggers a lot of things to happen, so we do it last

		if data!=None:
			self.set_data(data)

		# Boxviewer subwidget (details of a single box)
		self.boxviewer=EMBoxViewer()
		#self.app().attach_child(self.boxviewer)

		# Boxes Viewer (z projections of all boxes)
		self.boxesviewer=EMImageMXWidget()
		#self.app().attach_child(self.boxesviewer)
		#self.boxesviewer.show()
		#self.boxesviewer.set_mouse_mode("App")
		#self.boxesviewer.setWindowTitle("Particle List")

		# Average viewer shows results of background tomographic processing
#		self.averageviewer=EMAverageViewer(self)
		#self.averageviewer.show()

		QtCore.QObject.connect(self.boxesviewer,QtCore.SIGNAL("mx_image_selected"),self.img_selected)
		self.e = None

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

		print "\nDatafile set, see!", self.datafile, type(self.datafile)

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
		self.update_all()

	def get_cube(self,x,y,z):
		"""Returns a box-sized cube at the given center location"""
		bs=self.boxsize()

		if self.yshort:
			if self.data!=None:
				r=self.data.get_clip(Region(x-bs/2,z-bs/2,y-bs/2,bs,bs,bs))
				if options.normproc:
					r.process_inplace(options.normproc)
				r.process_inplace("xform",{"transform":Transform({"type":"eman","alt":90.0})})
				r.process_inplace("xform.mirror",{"axis":"z"})
			elif self.datafile!=None:
				r=EMData(self.datafile,0,0,Region(x-bs/2,z-bs/2,y-bs/2,bs,bs,bs))
				if options.normproc:
					r.process_inplace(options.normproc)
				r.process_inplace("xform",{"transform":Transform({"type":"eman","alt":90.0})})
				r.process_inplace("xform.mirror",{"axis":"z"})
			else: return None

		else :
			if self.data!=None:
				r=self.data.get_clip(Region(x-bs/2,y-bs/2,z-bs/2,bs,bs,bs))
			elif self.datafile!=None:
				r=EMData(self.datafile,0,0,Region(x-bs/2,y-bs/2,z-bs/2,bs,bs,bs))
			else: return None

		if self.apix!=0 :
			r["apix_x"]=self.apix
			r["apix_y"]=self.apix
			r["apix_z"]=self.apix

		if options.normproc:
			r.process_inplace(options.normproc)
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
		if self.boxsize()==self.oldboxsize:
			return
		self.oldboxsize=self.boxsize()

		cb=self.curbox
		for i in range(len(self.boxes)):
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
		self.update_xy()

	def event_nlayers(self):
		self.update_all()

	def event_filter(self):
		self.update_all()

	def event_localbox(self,tog):
		self.update_sides()

	def boxsize(self):
		return 32 #int(self.wboxsize.getValue())

	def nlayers(self):
		return int(self.wnlayers.value())

	def depth(self):
		return int(self.wdepth.value())

	def scale(self):
		return self.wscale.getValue()

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
		if options.helixboxer:
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
				self.boxes.append(b2)
				self.update_box(len(self.boxes)-1)
		f.close()

	def menu_file_save_boxloc(self):
		shrinkf=self.shrink 								#jesus

		fsp=str(QtGui.QFileDialog.getSaveFileName(self, "Select output text file"))

		out=file(fsp,"w")
		if options.helixboxer:
			for b in self.helixboxes:
				out.write("%d\t%d\t%d\t%d\t%d\t%d\n"%(b[0]*shrinkf,b[1]*shrinkf,b[2]*shrinkf,b[3]*shrinkf,b[4]*shrinkf,b[5]*shrinkf))
		else:
			for b in self.boxes:
				out.write("%d\t%d\t%d\n"%(b[0]*shrinkf,b[1]*shrinkf,b[2]*shrinkf))
		out.close()

	def menu_file_save_boxes(self):
		fsp=os.path.basename(str(QtGui.QFileDialog.getSaveFileName(self, "Select output file (numbers added)")))

		fspprjs=fsp.replace('.','_prjs.hdf')
		prj=EMData() #Dummy

		progress = QtGui.QProgressDialog("Saving", "Abort", 0, len(self.boxes),None)
		if options.helixboxer:
			for i,b in enumerate(self.helixboxes):
				img = self.extract_subtomo_box(self.get_extended_a_vector(b), cshrink=self.shrink)

				#img['origin_x'] = 0
				#img['origin_y'] = 0
				#img['origin_z'] = 0

				if self.normalize:
					img.process_inplace(normalize)
				#img=img.process('normalize.edgemean')

				#if fsp[:4].lower()=="bdb:":
				#	img.write_image(os.path.join(options.path,"%s_%03d"%(fsp,i)),0)

				if "." in fsp:
					img.write_image(os.path.join(options.path,"%s_%03d.%s"%(fsp.rsplit(".",1)[0],i,fsp.rsplit(".",1)[1])))
				else:
					QtGui.QMessageBox.warning(None,"Error","Please provide a valid image file extension. The numerical sequence will be inserted before the extension.")
					return

				progress.setValue(i+1)
				if progress.wasCanceled() : break
		else:
			for i,b in enumerate(self.boxes):
				#img=self.get_cube(b[0],b[1],b[2])
				bs=self.boxsize()
				shrinkf=self.shrink
				if shrinkf >1:
					bs=bs*shrinkf

				contrast=self.contrast
				center=self.center

				if self.yshort:
					ret = unbinned_extractor(options,bs,b[0],b[2],b[1],shrinkf,contrast,center,args[0])
					img = ret[0]
					prj = ret[1]
				else:
					ret = unbinned_extractor(options,bs,b[0],b[1],b[2],shrinkf,contrast,center,args[0])
					img = ret[0]
					prj = ret[1]

				if "." in fsp:
					img.write_image(os.path.join(options.path,"%s_%03d.%s"%(fsp.rsplit(".",1)[0],i,fsp.rsplit(".",1)[1])))
					prj.write_image(fspprjs,-1)

				else:
					QtGui.QMessageBox.warning(None,"Error","Please provide a valid image file extension. The numerical sequence will be inserted before the extension.")
					return

				progress.setValue(i+1)
				if progress.wasCanceled() : break

	def menu_file_save_boxes_stack(self):

		fsp=os.path.join(options.path,os.path.basename(str(QtGui.QFileDialog.getSaveFileName(self, "Select output file (.hdf supported only)"))))
		#if fsp[:4].lower()!="bdb:" and fsp[-4:].lower()!=".hdf" :


		if fsp[-4:].lower()!=".hdf" :
			QtGui.QMessageBox.warning(None,"Error","3-D stacks supported only for .hdf files")
			return

		fspprjs=fsp.replace('.hdf','_prjs.hdf')
		prj=EMData() #Dummy

		progress = QtGui.QProgressDialog("Saving", "Abort", 0, len(self.boxes),None)
		if options.helixboxer:
			for i,b in enumerate(self.helixboxes):
				img = self.extract_subtomo_box(self.get_extended_a_vector(b), cshrink=self.shrink)

				#img['origin_x'] = 0
				#img['origin_y'] = 0
				#img['origin_z'] = 0
				if self.normalize:
					e.process_inplace(normalize)
				#img=img.process('normalize.edgemean')

				img.write_image(fsp,i)

				progress.setValue(i+1)
				if progress.wasCanceled():
					break
		else:
			for i,b in enumerate(self.boxes):
				#img=self.get_cube(b[0],b[1],b[2])
				bs=self.boxsize()
				shrinkf=self.shrink
				if shrinkf >1:
					bs=bs*shrinkf

				contrast=self.contrast
				center=self.center

				if self.yshort:
					ret = unbinned_extractor(options,bs,b[0],b[2],b[1],shrinkf,contrast,center,args[0])
					img = ret[0]
					prj = ret[1]
				else:
					ret = unbinned_extractor(options,bs,b[0],b[1],b[2],shrinkf,contrast,center,args[0])
					img = ret[0]
					prj = ret[1]

				#img['origin_x'] = 0
				#img['origin_y'] = 0
				#img['origin_z'] = 0
				if self.normalize:
					img.process_inplace(normalize)
				#img=img.process('normalize.edgemean')

				img.write_image(fsp,i)
				prj.write_image(fspprjs,-1)

				progress.setValue(i+1)
				if progress.wasCanceled():
					break


	def menu_file_quit(self):
		self.close()

	def transform_coords(self, point, xform):
		xvec = xform.get_matrix()
		return [xvec[0]*point[0] + xvec[4]*point[1] + xvec[8]*point[2] + xvec[3], xvec[1]*point[0] + xvec[5]*point[1] + xvec[9]*point[2] + xvec[7], xvec[2]*point[0] + xvec[6]*point[1] + xvec[10]*point[2] + xvec[11]]

	def extract_subtomo_box(self, helixbox, cshrink=1, tomogram=None):
		""" Retruns an extracted subtomogram box"""
		if tomogram:
			# Only scale the helix boxer values transiently
			x1 = round(helixbox[0]*cshrink)
			y1 = round(helixbox[1]*cshrink)
			z1 = round(helixbox[2]*cshrink)
			x2 = round(helixbox[3]*cshrink)
			y2 = round(helixbox[4]*cshrink)
			z2 = round(helixbox[5]*cshrink)
	
			bs=self.boxsize()/2
			# Get the extended vector based on boxsize
			a = Vec3f((x2-x1), (y2-y1), (z2-z1))	# Find the a, the long vector
			tcs = self.get_box_coord_system([x1,y1,z1,x2,y2,z2])							# Get the local coord system
			# Get the new coord system
			# First extract a subtomo gram bounding region from the tomogram so we do have to read the whole bloody thing in!
			rv = [self.transform_coords([0, -bs, -bs], tcs), self.transform_coords([0, bs, bs], tcs), self.transform_coords([0, bs, -bs], tcs), self.transform_coords([0, -bs, bs], tcs), self.transform_coords([a.length(), -bs, -bs], tcs), self.transform_coords([a.length(), bs, bs], tcs), self.transform_coords([a.length(), bs, -bs], tcs), self.transform_coords([a.length(), -bs, bs], tcs)]
			rvmin = [int(min([i[0] for i in rv])), int(min([i[1] for i in rv])), int(min([i[2] for i in rv]))]	# Min bounding box extension
			rvmax = [int(max([i[0] for i in rv])), int(max([i[1] for i in rv])), int(max([i[2] for i in rv]))]	# Max bounding box extension
			r = Region(rvmin[0],rvmin[1],rvmin[2],rvmax[0]-rvmin[0],rvmax[1]-rvmin[1],rvmax[2]-rvmin[2])		# Extract the region
			e = EMData()
			e.read_image(tomogram,0,False,r)
			e.set_attr("source_path", tomogram)
			e["ptcl_source_image"]=tomogram
			e["ptcl_source_coord"]=((rvmin[0]+rvmax[0])/2,(rvmin[1]+rvmax[1])/2,(rvmin[2]+rvmax[2])/2)
			# Next adjust the transform matrix to move it to the origin
			origin = self.transform_coords([0,0,0], tcs)
			tcs.set_trans(origin[0] - rvmin[0], origin[1] - rvmin[1], origin[2] - rvmin[2])
	
			return e.extract_box(tcs, Region(0, -bs, -bs, a.length(), bs, bs))
		else:
			return

	def get_averager(self):
		"""returns an averager of the appropriate type for generating projection views"""
		if self.wmaxmean.isChecked() : return Averagers.get("minmax",{"max":1})

		return Averagers.get("mean")

	def update_sides(self):
		"""updates xz and yz views due to a new center location"""

		print "\n\n\n\n\nIn update sides, self.datafile is", self.datafile
		print "\n\n\n\n"

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
		bs=self.boxsize()

		# update shape display
		if self.wlocalbox.isChecked():
			xzs=self.xzview.get_shapes()
			if options.helixboxer:
				for i in range(0, len(self.boxes), 2):
					if i + 1 >= len(self.boxes):
						break
					#if abs(self.boxes[i][1] - zc) < bs/2 or abs(self.boxes[i+1][2] - zc) < bs/2:
					if (self.boxes[i][1]<self.cury+bs/2 and self.boxes[i][1]>self.cury-bs/2) or (self.boxes[i+1][1]<self.cury+bs/2 and self.boxes[i+1][1]>self.cury-bs/2):
						xzs[i][0]="rect"
						xzs[str(i/2)+"helix"][0]="line"
						xzs[i+1][0]="rect"
					else:
						xzs[i][0]="hidden"
						xzs[str(i/2)+"helix"][0]="hidden"
						xzs[i+1][0]="hidden"
			else:
					for i in range(len(self.boxes)):
						if self.boxes[i][1]<self.cury+bs/2 and self.boxes[i][1]>self.cury-bs/2:
							xzs[i][0]="rect"
						else:
							xzs[i][0]="hidden"

			zys=self.zyview.get_shapes()
			if options.helixboxer:
				for i in range(0, len(self.boxes), 2):
					if i + 1 >= len(self.boxes):
						break
					if (self.boxes[i][0]<self.curx+bs/2 and self.boxes[i][0]>self.curx-bs/2) or (self.boxes[i+1][0]<self.curx+bs/2 and self.boxes[i+1][0]>self.curx-bs/2):
						zys[i][0]="rect"
						zys[str(i/2)+"helix"][0]="line"
						zys[i+1][0]="rect"
					else:
						zys[i][0]="hidden"
						zys[str(i/2)+"helix"][0]="hidden"
						zys[i+1][0]="hidden"
			else:
				for i in range(len(self.boxes)):
					if self.boxes[i][0]<self.curx+bs/2 and self.boxes[i][0]>self.curx-bs/2:
						zys[i][0]="rect"
					else:
						zys[i][0]="hidden"
		else :
			xzs=self.xzview.get_shapes()
			zys=self.zyview.get_shapes()
			if options.helixboxer:
				for i in range(0, len(self.boxes), 2):
					if i + 1 >= len(self.boxes):
						break
					xzs[i][0]="rect"
					xzs[str(i/2)+"helix"][0]="line"
					xzs[i+1][0]="rect"
					zys[i][0]="rect"
					zys[str(i/2)+"helix"][0]="line"
					zys[i+1][0]="rect"
			else:
				for i in range(len(self.boxes)):
					xzs[i][0]="rect"
					zys[i][0]="rect"

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

		print "\n\n\n\n\nIn update_xy, self.datafile is", self.datafile
		print "\n\n\n\n"

		if self.datafile==None and self.data==None:
			return



		# Boxes should also be limited by default in the XY view
		if len(self.boxes) > 0:
			zc=self.wdepth.value()
			#print "The current depth is", self.wdepth.value()
			bs=self.boxsize()
			xys=self.xyview.get_shapes()
			if options.helixboxer:
				for i in range(0, len(self.boxes), 2):
					if i + 1 >= len(self.boxes):
						break
					if abs(self.boxes[i][2] - zc) < bs/2 or abs(self.boxes[i+1][2] - zc) < bs/2:
						xys[i][0]="rect"
						xys[str(i/2)+"helix"][0]="line"
						xys[i+1][0]="rect"
					else:
						xys[i][0]="hidden"
						xys[str(i/2)+"helix"][0]="hidden"
						xys[i+1][0]="hidden"
			else:
				for i in range(len(self.boxes)):
					#print "the z coord of box %d is %d" %(i,self.boxes[i][2])
					#print "therefore the criteria to determine whether to display it is", abs(self.boxes[i][2] - zc)
					if abs(self.boxes[i][2] - zc) < bs/2:
						#print "Which is less than half the box thus it survives"
						xys[i][0]="rect"
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

		print "\n\nIn update xy, av and type are", av, type(av)

		if self.wfilt.getValue()!=0.0:

			av.process_inplace("filter.lowpass.gauss",{"cutoff_freq":1.0/self.wfilt.getValue(),"apix":self.apix})
		self.xyview.set_data(av)

	def update_all(self):
		"""redisplay of all widgets"""

		print "\n\n\n\n\nIn update all, self.datafile is", self.datafile
		print "\n\n\n\n"
		if self.datafile==None and self.data==None:
			return



		self.update_xy()
		self.update_sides()

		#self.xyview.update()
		#self.xzview.update()
		#self.zyview.update()


	def inside_box(self,n,x=-1,y=-1,z=-1):
		"""Checks to see if a point in image coordinates is inside box number n. If any value is negative, it will not be checked."""
		box=self.boxes[n]
		if x>=0 and (x<box[0]-self.boxsize()/2 or x>box[0]+self.boxsize()/2) : return False
		if y>=0 and (y<box[1]-self.boxsize()/2 or y>box[1]+self.boxsize()/2) : return False
		if z>=0 and (z<box[2]-self.boxsize()/2 or z>box[2]+self.boxsize()/2) : return False
		return True

	def do_deletion(self, n, delimgs=True):
		""" Helper for del_box"""
		if n==len(self.boxes)-1 :
			self.boxes.pop()
			if delimgs:
				self.boxesimgs.pop()
				self.boxesviewer.set_data(self.boxesimgs)
				self.boxesviewer.update()
			self.xyview.del_shape(n)
			self.xzview.del_shape(n)
			self.zyview.del_shape(n)
			self.curbox=-1
			self.xyview.update()
			self.xzview.update()
			self.zyview.update()
			#self.update()
		else :
			a=self.boxes.pop()
			self.boxes[n]=a
			if delimgs:
				a=self.boxesimgs.pop()
				self.boxesimgs[n]=a
				self.boxesviewer.set_data(self.boxesimgs)
				self.boxesviewer.set_selected([],True)
				self.boxesviewer.update()
			self.xyview.del_shape(len(self.boxes))
			self.xzview.del_shape(len(self.boxes))
			self.zyview.del_shape(len(self.boxes))
			self.update_box(n,True)
#			self.update()

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
		if options.helixboxer:
			if n + 1 == len(self.boxes) and len(self.boxes) % 2 == 1: 	# Delete unpaired box
				self.do_deletion(n, delimgs=False)
			else:								# Delete box pairs
				if n % 2:
					self.do_helix_deletion(int(n/2))
					self.do_deletion(n, delimgs=False)
					self.do_deletion(n-1, delimgs=False)
				else:
					self.do_helix_deletion(int(n/2))
					self.do_deletion(n+1, delimgs=False)
					self.do_deletion(n, delimgs=False)
				return "DELHELIX"	# If we have deleted a pair do not reset the pair toggle/counter
		else:
			self.do_deletion(n)

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
		bs = self.boxsize()
		return [(helixbox[0] - a[0]*bs/2),(helixbox[1] - a[1]*bs/2),(helixbox[2] - a[2]*bs/2),(helixbox[3] + a[0]*bs/2),(helixbox[4] + a[1]*bs/2),(helixbox[5] + a[2]*bs/2)]


	def update_helixbox(self, n, quiet=False):
		"""
		Update a helix box
		"""
		if n > len(self.helixboxes)-1: return	# Some boxes may not be paired
		helixbox = self.get_extended_a_vector(self.helixboxes[n])

		key = str(n)+"helix"
		if options.yshort:
			self.xyview.add_shape(key,EMShape(("line",.2,.2,.8, helixbox[0], helixbox[2], helixbox[3], helixbox[5],2)))
			self.xzview.add_shape(key,EMShape(("line",.2,.2,.8, helixbox[0], helixbox[1], helixbox[3], helixbox[4],2)))
			self.zyview.add_shape(key,EMShape(("line",.2,.2,.8, helixbox[1], helixbox[2], helixbox[4], helixbox[5],2)))
		else:
			self.xyview.add_shape(key,EMShape(("line",.2,.2,.8, helixbox[0], helixbox[1], helixbox[3], helixbox[4],2)))
			self.xzview.add_shape(key,EMShape(("line",.2,.2,.8, helixbox[0], helixbox[2], helixbox[3], helixbox[5],2)))
			self.zyview.add_shape(key,EMShape(("line",.2,.2,.8, helixbox[2], helixbox[1], helixbox[5], helixbox[4],2)))
		self.xyview.update()
		self.xzview.update()
		self.zyview.update()

		if not quiet and options.helixboxer:
			hb = self.extract_subtomo_box(helixbox, cshrink=self.shrink)
			self.boxviewer.set_data(hb)

			proj=hb.process("misc.directional_sum",{"axis":"z"})
			try: self.boxesimgs[n]=proj
			except:
				for i in range(len(self.boxesimgs),n+1): self.boxesimgs.append(None)
				self.boxesimgs[n]=proj
			self.boxesviewer.set_data(self.boxesimgs)
			self.boxesviewer.update()

		if n!=self.curbox and options.helixboxer:
			self.boxesviewer.set_selected((n,),True)

	def update_box(self,n,quiet=False):
		"""After adjusting a box, call this"""
#		print "upd ",n,quiet

		try:
			box=self.boxes[n]
		except IndexError:
			return
		bs2=self.boxsize()/2

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
		self.xyview.add_shape(n,EMShape(("rect",.2,.2,.8,box[0]-bs2,box[1]-bs2,box[0]+bs2,box[1]+bs2,2)))
		self.xyview.add_shape("xl",EMShape(("line",.8,.8,.1,0,box[1],self.datasize[0],box[1],1)))
		self.xyview.add_shape("yl",EMShape(("line",.8,.8,.1,box[0],0,box[0],self.datasize[1],1)))
		self.xzview.add_shape(n,EMShape(("rect",.2,.2,.8,box[0]-bs2,box[2]-bs2,box[0]+bs2,box[2]+bs2,2)))
		self.xzview.add_shape("xl",EMShape(("line",.8,.8,.1,0,box[2],self.datasize[0],box[2],1)))
		self.xzview.add_shape("zl",EMShape(("line",.8,.8,.1,box[0],0,box[0],self.datasize[2],1)))
		self.zyview.add_shape(n,EMShape(("rect",.2,.2,.8,box[2]-bs2,box[1]-bs2,box[2]+bs2,box[1]+bs2,2)))
		self.zyview.add_shape("yl",EMShape(("line",.8,.8,.1,box[2],0,box[2],self.datasize[1],1)))
		self.zyview.add_shape("zl",EMShape(("line",.8,.8,.1,0,box[1],self.datasize[2],box[1],1)))

		if self.depth()!=box[2]:
			self.wdepth.setValue(box[2])
		else:
			self.xyview.update()
		self.update_sides()

		# For speed, we turn off updates while dragging a box around. Quiet is set until the mouse-up
		if not quiet and not options.helixboxer:
			# Get the cube from the original data (normalized)
			cube=self.get_cube(*box)
			self.boxviewer.set_data(cube)

			# Make a z projection and store it in the list of all boxes
			proj=cube.process("misc.directional_sum",{"axis":"z"})
			try: self.boxesimgs[n]=proj
			except:
				for i in range(len(self.boxesimgs),n+1): self.boxesimgs.append(None)
				self.boxesimgs[n]=proj
			self.boxesviewer.set_data(self.boxesimgs)
			self.boxesviewer.update()

		if n!=self.curbox and not options.helixboxer:
			self.boxesviewer.set_selected((n,),True)

		self.curbox=n


	def img_selected(self,event,lc):
#		print "sel",lc[0]
		if event.modifiers()&Qt.ShiftModifier:
			self.del_box(lc[0])
		else:
			self.update_box(lc[0])
		if self.curbox>=0 :
			box=self.boxes[self.curbox]
			self.xyview.scroll_to(box[0],box[1])
			self.xzview.scroll_to(None,box[2])
			self.zyview.scroll_to(box[2],None)

	def add_helix_box(self, xf, yf, zf, xi, yi, zi):
		print xf, yf, zf, xi, yi, zi
		if options.yshort:
			self.helixboxes.append([xf, zf, yf, xi, zi, yi])
		else:
			self.helixboxes.append([xf, yf, zf, xi, yi, zi])

	def xy_down(self,event):
		x,y=self.xyview.scr_to_img((event.x(),event.y()))
		x,y=int(x),int(y)
		self.xydown=None
		if x<0 or y<0 : return		# no clicking outside the image (on 2 sides)

		for i in range(len(self.boxes)):
			if self.inside_box(i,x,y) :
				if event.modifiers()&Qt.ShiftModifier:
					if self.del_box(i) != "DELHELIX": self.firsthbclick = None
				else:
					self.xydown=(i,x,y,self.boxes[i][0],self.boxes[i][1])
					if options.helixboxer: self.update_helixbox(int(i/2))
					self.update_box(i)
				break
		else:
#			if x>self.boxsize()/2 and x<self.datasize[0]-self.boxsize()/2 and y>self.boxsize()/2 and y<self.datasize[1]-self.boxsize()/2 and self.depth()>self.boxsize()/2 and self.depth()<self.datasize[2]-self.boxsize()/2 :
			if not event.modifiers()&Qt.ShiftModifier:
				###########
				if options.helixboxer:	# Only create a helixbox every 2 clicks
					if self.firsthbclick:
						self.add_helix_box(x, y, self.depth(), self.firsthbclick[0], self.firsthbclick[1], self.firsthbclick[2])
						self.firsthbclick = None
						self.update_helixbox(len(self.helixboxes)-1)
					else:
						self.firsthbclick = [x, y, self.depth()]
				###########
				self.boxes.append(([x,y,self.depth()]))
				self.xydown=(len(self.boxes)-1,x,y,x,y)		# box #, x down, y down, x box at down, y box at down
				self.update_box(self.xydown[0])

		if self.curbox>=0:
			box=self.boxes[self.curbox]
			self.xzview.scroll_to(None,box[2])
			self.zyview.scroll_to(box[2],None)

	def xy_drag(self,event):
		if self.xydown==None : return

		x,y=self.xyview.scr_to_img((event.x(),event.y()))
		x,y=int(x),int(y)

		dx=x-self.xydown[1]
		dy=y-self.xydown[2]
		if options.helixboxer:
			if len(self.boxes) % 2 == 0 or (self.xydown[0] != len(self.boxes)-1):	# Only update the helix boxer if it is paired, otherwise treat it as a regular box
				hb = self.helixboxes[int(self.xydown[0]/2)]
				if self.xydown[0] % 2 == 0:
					hb[3] = dx+self.xydown[3]
					hb[4] = dy+self.xydown[4]
				else:
					hb[0] = dx+self.xydown[3]
					hb[1] = dy+self.xydown[4]
				self.update_helixbox(int(self.xydown[0]/2))
			else:
				self.firsthbclick[0] = x
				self.firsthbclick[1] = y

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

	def xz_down(self,event):
		x,z=self.xzview.scr_to_img((event.x(),event.y()))
		x,z=int(x),int(z)
		self.xzdown=None
		if x<0 or z<0 : return		# no clicking outside the image (on 2 sides)

		for i in range(len(self.boxes)):
			if (not self.wlocalbox.isChecked() and self.inside_box(i,x,-1,z)) or self.inside_box(i,x,self.cury,z) :
				if event.modifiers()&Qt.ShiftModifier:
					if self.del_box(i) != "DELHELIX": self.firsthbclick = None
				else :
					self.xzdown=(i,x,z,self.boxes[i][0],self.boxes[i][2])
					if options.helixboxer: self.update_helixbox(int(i/2))
					self.update_box(i)
				break
		else:
			if not event.modifiers()&Qt.ShiftModifier:
				###########
				if options.helixboxer:	# Only create a helixbox every 2 clicks
					if self.firsthbclick:
						self.add_helix_box(x, self.cury, z, self.firsthbclick[0], self.firsthbclick[1], self.firsthbclick[2])
						self.firsthbclick = None
						self.update_helixbox(len(self.helixboxes)-1)
					else:
						self.firsthbclick = [x, self.cury, z]
				###########
				self.boxes.append(([x,self.cury,z]))
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
		if options.helixboxer:
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
		self.xydown=None
		if z<0 or y<0 : return		# no clicking outside the image (on 2 sides)

		for i in range(len(self.boxes)):
			if (not self.wlocalbox.isChecked() and self.inside_box(i,-1,y,z)) or  self.inside_box(i,self.curx,y,z):
				if event.modifiers()&Qt.ShiftModifier:
					if self.del_box(i) != "DELHELIX": self.firsthbclick = None
				else :
					self.zydown=(i,z,y,self.boxes[i][2],self.boxes[i][1])
					if options.helixboxer: self.update_helixbox(int(i/2))
					self.update_box(i)
				break
		else:
			if not event.modifiers()&Qt.ShiftModifier:
				###########
				if options.helixboxer:	# Only create a helixbox every 2 clicks
					if self.firsthbclick:
						self.add_helix_box(self.curx, y, z, self.firsthbclick[0], self.firsthbclick[1], self.firsthbclick[2])
						self.firsthbclick = None
						self.update_helixbox(len(self.helixboxes)-1)
					else:
						self.firsthbclick = [self.curx, y, z]
				###########
				self.boxes.append(([self.curx,y,z]))
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
		if options.helixboxer:
			if len(self.boxes) % 2 == 0 or (self.zydown[0] != len(self.boxes)-1):	# Only update the helix boxer if it is paired, otherwise treat it as a regular box
				hb = self.helixboxes[int(self.zydown[0]/2)]
				if self.zydown[0] % 2 == 0:
					hb[5] = dz+self.zydown[3]
					hb[4] = dy+self.zydown[4]
				else:
					hb[2] =  dz+self.zydown[3]
					hb[1] = dy+self.zydown[4]
				self.update_helixbox(int(self.zydown[0]/2))
			else:
				self.firsthbclick[2] = z
				self.firsthbclick[1] = y

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

	def closeEvent(self,event):
		print "Exiting"
		self.boxviewer.close()
		self.boxesviewer.close()
#		self.averageviewer.close()
		event.accept()
		#self.app().close_specific(self)
		self.emit(QtCore.SIGNAL("module_closed")) # this signal is important when e2ctf is being used by a program running its own event loop

	#def closeEvent(self,event):
		#self.target().done()

class TomoSegInspector(QtGui.QWidget):
	
	def __init__(self):
		super(TomoSegInspector,self).__init__()
		
		self.setWindowTitle("TomoSeg Segmentation Tools")
		
		self.mintreewidth = 250		# minimum width of the tree
		self.mincontrolwidth = 0
		
		vbox = QtGui.QVBoxLayout(self)
		self.inspectortab = QtGui.QTabWidget()
		self.inspectortab.addTab(self.getToolsWidget(), "Tools")
		self.inspectortab.addTab(self.getTreeWidget(), "Annotations")
		self.inspectortab.addTab(self.getUtilsWidget(), "Utils")
		vbox.addWidget(self.inspectortab)
		
		self.setLayout(vbox)
		self.updateGeometry()

	def getToolsWidget(self):
		tooltabs = QtGui.QTabWidget()
		tooltabs.addTab(self.getAutomaticTools(), "Automatic")
		tooltabs.addTab(self.getSemiAutomaticTools(), "Interactive")
		tooltabs.addTab(self.getManualTools(), "Manual")
		return tooltabs
	
	def getAutomaticTools(self):
 		widget = QtGui.QWidget()
		
		hbox = QtGui.QHBoxLayout()
		
		frame = QtGui.QFrame()
		grid = QtGui.QGridLayout(widget)
		
		self.rotatetool = EMANToolButton()
# 		self.rotatetool.setIcon(QtGui.QIcon(QtGui.QPixmap(rotateicon)))
 		self.rotatetool.setToolTip("Description")
 		self.rotatetool_label = QtGui.QLabel("Name")
		
		self.translatetool =EMANToolButton()
# 		self.translatetool.setIcon(QtGui.QIcon(QtGui.QPixmap(crosshairsicon)))
 		self.translatetool.setToolTip("Description")
 		self.translatetool_label = QtGui.QLabel("Name")
		
		self.ztranslate = EMANToolButton()
# 		self.ztranslate.setIcon(QtGui.QIcon(QtGui.QPixmap(ztransicon)))
 		self.ztranslate.setToolTip("Description")
 		self.ztranslatetool_label = QtGui.QLabel("Name")
		
		self.scaletool = EMANToolButton()
# 		self.scaletool.setIcon(QtGui.QIcon(QtGui.QPixmap(scaleicon)))
 		self.scaletool.setToolTip("Description")
 		self.scaletool_label = QtGui.QLabel("Name")
		
		self.rulertool = EMANToolButton()
# 		self.rulertool.setIcon(QtGui.QIcon(QtGui.QPixmap(rulericon)))
 		self.rulertool.setToolTip("Description")
 		self.rulertool_label = QtGui.QLabel("Name")
		
		self.selectiontool = EMANToolButton()
# 		self.selectiontool.setIcon(QtGui.QIcon(QtGui.QPixmap(selectionicon)))
 		self.selectiontool.setToolTip("Description")
 		self.selectiontool_label = QtGui.QLabel("Name")
		
		self.multiselectiontool = EMANToolButton()
# 		self.multiselectiontool.setIcon(QtGui.QIcon(QtGui.QPixmap(multiselectoricon)))
 		self.multiselectiontool.setToolTip("Description")
 		self.multiselectiontool_label = QtGui.QLabel("Name")
		
		self.linetool = EMANToolButton()
# 		self.linetool.setIcon(QtGui.QIcon(QtGui.QPixmap(lineicon)))
		self.linetool.setToolTip("Description")
		self.linetool_label = QtGui.QLabel("Name")
		
		self.cubetool = EMANToolButton()
# 		self.cubetool.setIcon(QtGui.QIcon(QtGui.QPixmap(cubeicon)))
 		self.cubetool.setToolTip("Description")
 		self.cubetool_label = QtGui.QLabel("Name")
		
		self.spheretool = EMANToolButton()
# 		self.spheretool.setIcon(QtGui.QIcon(QtGui.QPixmap(sphereicon)))
 		self.spheretool.setToolTip("Description")
 		self.spheretool_label = QtGui.QLabel("Name")
		
		self.cylindertool = EMANToolButton()
# 		self.cylindertool.setIcon(QtGui.QIcon(QtGui.QPixmap(cylindericon)))
 		self.cylindertool.setToolTip("Description")
 		self.cylindertool_label = QtGui.QLabel("Name")
		
		self.conetool = EMANToolButton()
# 		self.conetool.setIcon(QtGui.QIcon(QtGui.QPixmap(coneicon)))
 		self.conetool.setToolTip("Description")
 		self.conetool_label = QtGui.QLabel("Name")
		
		self.texttool = EMANToolButton()
# 		self.texttool.setIcon(QtGui.QIcon(QtGui.QPixmap(texticon)))
 		self.texttool.setToolTip("Description")
 		self.texttool_label = QtGui.QLabel("Name")
 		
		self.datatool = EMANToolButton()
# 		self.datatool.setIcon(QtGui.QIcon(QtGui.QPixmap(dataicon)))
 		self.datatool.setToolTip("Description")
		self.datatool_label = QtGui.QLabel("Name")
		
		self.apptool = EMANToolButton()
#		self.apptool.setIcon(QtGui.QIcon(QtGui.QPixmap(appicon)))
 		self.apptool.setToolTip("Description")
		self.apptool_label = QtGui.QLabel("Name")

		# buttons
		grid.addWidget(self.selectiontool,1,0)
		grid.addWidget(self.multiselectiontool,2,0)
		grid.addWidget(self.translatetool,3,0)
		grid.addWidget(self.ztranslate,4,0)
		grid.addWidget(self.rotatetool,5,0)
		grid.addWidget(self.scaletool,6,0)
		grid.addWidget(self.rulertool,7,0)
		grid.addWidget(self.linetool,8,0)
		grid.addWidget(self.cubetool,9,0)
		grid.addWidget(self.spheretool,10,0)
		grid.addWidget(self.cylindertool,11,0)
		grid.addWidget(self.conetool,12,0)
		grid.addWidget(self.texttool,13,0)
		grid.addWidget(self.datatool,14,0)
		grid.addWidget(self.apptool,15,0)
		grid.setAlignment(QtCore.Qt.AlignLeft)
		# labels
		#grid.addWidget(toollabel,0,0)
		grid.addWidget(self.selectiontool_label,1,1)
		grid.addWidget(self.multiselectiontool_label,2,1)
		grid.addWidget(self.translatetool_label,3,1)
		grid.addWidget(self.ztranslatetool_label,4,1)
		grid.addWidget(self.rotatetool_label,5,1)
		grid.addWidget(self.scaletool_label,6,1)
		grid.addWidget(self.rulertool_label,7,1)
		grid.addWidget(self.linetool_label,8,1)
		grid.addWidget(self.cubetool_label,9,1)
		grid.addWidget(self.spheretool_label,10,1)
		grid.addWidget(self.cylindertool_label,11,1)
		grid.addWidget(self.conetool_label,12,1)
		grid.addWidget(self.texttool_label,13,1)
		grid.addWidget(self.datatool_label,14,1)
		grid.addWidget(self.apptool_label,15,1)
		grid.setAlignment(QtCore.Qt.AlignLeft)	
		
		frame.setLayout(grid)
		hbox.addWidget(frame)
		
		self.stacked_widget = QtGui.QStackedWidget()
		self.stacked_widget.setFrameShape(QtGui.QFrame.StyledPanel)
		hbox.addWidget(self.stacked_widget)
		
		widget.setLayout(hbox) #grid
		
		return widget

	def getSemiAutomaticTools(self):
 		widget = QtGui.QWidget()
		
		hbox = QtGui.QHBoxLayout()
		
		frame = QtGui.QFrame()
		grid = QtGui.QGridLayout(widget)
		
		self.rotatetool = EMANToolButton()
# 		self.rotatetool.setIcon(QtGui.QIcon(QtGui.QPixmap(rotateicon)))
 		self.rotatetool.setToolTip("Description")
 		self.rotatetool_label = QtGui.QLabel("Name")
		
		self.translatetool =EMANToolButton()
# 		self.translatetool.setIcon(QtGui.QIcon(QtGui.QPixmap(crosshairsicon)))
 		self.translatetool.setToolTip("Description")
 		self.translatetool_label = QtGui.QLabel("Name")
		
		self.ztranslate = EMANToolButton()
# 		self.ztranslate.setIcon(QtGui.QIcon(QtGui.QPixmap(ztransicon)))
 		self.ztranslate.setToolTip("Description")
 		self.ztranslatetool_label = QtGui.QLabel("Name")
		
		self.scaletool = EMANToolButton()
# 		self.scaletool.setIcon(QtGui.QIcon(QtGui.QPixmap(scaleicon)))
 		self.scaletool.setToolTip("Description")
 		self.scaletool_label = QtGui.QLabel("Name")
		
		self.rulertool = EMANToolButton()
# 		self.rulertool.setIcon(QtGui.QIcon(QtGui.QPixmap(rulericon)))
 		self.rulertool.setToolTip("Description")
 		self.rulertool_label = QtGui.QLabel("Name")
		
		self.selectiontool = EMANToolButton()
# 		self.selectiontool.setIcon(QtGui.QIcon(QtGui.QPixmap(selectionicon)))
 		self.selectiontool.setToolTip("Description")
 		self.selectiontool_label = QtGui.QLabel("Name")
		
		self.multiselectiontool = EMANToolButton()
# 		self.multiselectiontool.setIcon(QtGui.QIcon(QtGui.QPixmap(multiselectoricon)))
 		self.multiselectiontool.setToolTip("Description")
 		self.multiselectiontool_label = QtGui.QLabel("Name")
		
		self.linetool = EMANToolButton()
# 		self.linetool.setIcon(QtGui.QIcon(QtGui.QPixmap(lineicon)))
		self.linetool.setToolTip("Description")
		self.linetool_label = QtGui.QLabel("Name")
		
		self.cubetool = EMANToolButton()
# 		self.cubetool.setIcon(QtGui.QIcon(QtGui.QPixmap(cubeicon)))
 		self.cubetool.setToolTip("Description")
 		self.cubetool_label = QtGui.QLabel("Name")
		
		self.spheretool = EMANToolButton()
# 		self.spheretool.setIcon(QtGui.QIcon(QtGui.QPixmap(sphereicon)))
 		self.spheretool.setToolTip("Description")
 		self.spheretool_label = QtGui.QLabel("Name")
		
		self.cylindertool = EMANToolButton()
# 		self.cylindertool.setIcon(QtGui.QIcon(QtGui.QPixmap(cylindericon)))
 		self.cylindertool.setToolTip("Description")
 		self.cylindertool_label = QtGui.QLabel("Name")
		
		self.conetool = EMANToolButton()
# 		self.conetool.setIcon(QtGui.QIcon(QtGui.QPixmap(coneicon)))
 		self.conetool.setToolTip("Description")
 		self.conetool_label = QtGui.QLabel("Name")
		
		self.texttool = EMANToolButton()
# 		self.texttool.setIcon(QtGui.QIcon(QtGui.QPixmap(texticon)))
 		self.texttool.setToolTip("Description")
 		self.texttool_label = QtGui.QLabel("Name")
 		
		self.datatool = EMANToolButton()
# 		self.datatool.setIcon(QtGui.QIcon(QtGui.QPixmap(dataicon)))
 		self.datatool.setToolTip("Description")
		self.datatool_label = QtGui.QLabel("Name")
		
		self.apptool = EMANToolButton()
#		self.apptool.setIcon(QtGui.QIcon(QtGui.QPixmap(appicon)))
 		self.apptool.setToolTip("Description")
		self.apptool_label = QtGui.QLabel("Name")

		# buttons
		grid.addWidget(self.selectiontool,1,0)
		grid.addWidget(self.multiselectiontool,2,0)
		grid.addWidget(self.translatetool,3,0)
		grid.addWidget(self.ztranslate,4,0)
		grid.addWidget(self.rotatetool,5,0)
		grid.addWidget(self.scaletool,6,0)
		grid.addWidget(self.rulertool,7,0)
		grid.addWidget(self.linetool,8,0)
		grid.addWidget(self.cubetool,9,0)
		grid.addWidget(self.spheretool,10,0)
		grid.addWidget(self.cylindertool,11,0)
		grid.addWidget(self.conetool,12,0)
		grid.addWidget(self.texttool,13,0)
		grid.addWidget(self.datatool,14,0)
		grid.addWidget(self.apptool,15,0)
		grid.setAlignment(QtCore.Qt.AlignLeft)
		# labels
		#grid.addWidget(toollabel,0,0)
		grid.addWidget(self.selectiontool_label,1,1)
		grid.addWidget(self.multiselectiontool_label,2,1)
		grid.addWidget(self.translatetool_label,3,1)
		grid.addWidget(self.ztranslatetool_label,4,1)
		grid.addWidget(self.rotatetool_label,5,1)
		grid.addWidget(self.scaletool_label,6,1)
		grid.addWidget(self.rulertool_label,7,1)
		grid.addWidget(self.linetool_label,8,1)
		grid.addWidget(self.cubetool_label,9,1)
		grid.addWidget(self.spheretool_label,10,1)
		grid.addWidget(self.cylindertool_label,11,1)
		grid.addWidget(self.conetool_label,12,1)
		grid.addWidget(self.texttool_label,13,1)
		grid.addWidget(self.datatool_label,14,1)
		grid.addWidget(self.apptool_label,15,1)
		grid.setAlignment(QtCore.Qt.AlignLeft)	
		
		frame.setLayout(grid)
		hbox.addWidget(frame)
		
		self.stacked_widget = QtGui.QStackedWidget()
		self.stacked_widget.setFrameShape(QtGui.QFrame.StyledPanel)
		hbox.addWidget(self.stacked_widget)
		
		widget.setLayout(hbox) #grid
		
		return widget

	def getManualTools(self):
		
 		widget = QtGui.QWidget()
		
		hbox = QtGui.QHBoxLayout()
		
		frame = QtGui.QFrame()
		grid = QtGui.QGridLayout(widget)
		
		self.rotatetool = EMANToolButton()
# 		self.rotatetool.setIcon(QtGui.QIcon(QtGui.QPixmap(rotateicon)))
 		self.rotatetool.setToolTip("Description")
 		self.rotatetool_label = QtGui.QLabel("Name")
		
		self.translatetool =EMANToolButton()
# 		self.translatetool.setIcon(QtGui.QIcon(QtGui.QPixmap(crosshairsicon)))
 		self.translatetool.setToolTip("Description")
 		self.translatetool_label = QtGui.QLabel("Name")
		
		self.ztranslate = EMANToolButton()
# 		self.ztranslate.setIcon(QtGui.QIcon(QtGui.QPixmap(ztransicon)))
 		self.ztranslate.setToolTip("Description")
 		self.ztranslatetool_label = QtGui.QLabel("Name")
		
		self.scaletool = EMANToolButton()
# 		self.scaletool.setIcon(QtGui.QIcon(QtGui.QPixmap(scaleicon)))
 		self.scaletool.setToolTip("Description")
 		self.scaletool_label = QtGui.QLabel("Name")
		
		self.rulertool = EMANToolButton()
# 		self.rulertool.setIcon(QtGui.QIcon(QtGui.QPixmap(rulericon)))
 		self.rulertool.setToolTip("Description")
 		self.rulertool_label = QtGui.QLabel("Name")
		
		self.selectiontool = EMANToolButton()
# 		self.selectiontool.setIcon(QtGui.QIcon(QtGui.QPixmap(selectionicon)))
 		self.selectiontool.setToolTip("Description")
 		self.selectiontool_label = QtGui.QLabel("Name")
		
		self.multiselectiontool = EMANToolButton()
# 		self.multiselectiontool.setIcon(QtGui.QIcon(QtGui.QPixmap(multiselectoricon)))
 		self.multiselectiontool.setToolTip("Description")
 		self.multiselectiontool_label = QtGui.QLabel("Name")
		
		self.linetool = EMANToolButton()
# 		self.linetool.setIcon(QtGui.QIcon(QtGui.QPixmap(lineicon)))
		self.linetool.setToolTip("Description")
		self.linetool_label = QtGui.QLabel("Name")
		
		self.cubetool = EMANToolButton()
# 		self.cubetool.setIcon(QtGui.QIcon(QtGui.QPixmap(cubeicon)))
 		self.cubetool.setToolTip("Description")
 		self.cubetool_label = QtGui.QLabel("Name")
		
		self.spheretool = EMANToolButton()
# 		self.spheretool.setIcon(QtGui.QIcon(QtGui.QPixmap(sphereicon)))
 		self.spheretool.setToolTip("Description")
 		self.spheretool_label = QtGui.QLabel("Name")
		
		self.cylindertool = EMANToolButton()
# 		self.cylindertool.setIcon(QtGui.QIcon(QtGui.QPixmap(cylindericon)))
 		self.cylindertool.setToolTip("Description")
 		self.cylindertool_label = QtGui.QLabel("Name")
		
		self.conetool = EMANToolButton()
# 		self.conetool.setIcon(QtGui.QIcon(QtGui.QPixmap(coneicon)))
 		self.conetool.setToolTip("Description")
 		self.conetool_label = QtGui.QLabel("Name")
		
		self.texttool = EMANToolButton()
# 		self.texttool.setIcon(QtGui.QIcon(QtGui.QPixmap(texticon)))
 		self.texttool.setToolTip("Description")
 		self.texttool_label = QtGui.QLabel("Name")
 		
		self.datatool = EMANToolButton()
# 		self.datatool.setIcon(QtGui.QIcon(QtGui.QPixmap(dataicon)))
 		self.datatool.setToolTip("Description")
		self.datatool_label = QtGui.QLabel("Name")
		
		self.apptool = EMANToolButton()
#		self.apptool.setIcon(QtGui.QIcon(QtGui.QPixmap(appicon)))
 		self.apptool.setToolTip("Description")
		self.apptool_label = QtGui.QLabel("Name")

		# buttons
		grid.addWidget(self.selectiontool,1,0)
		grid.addWidget(self.multiselectiontool,2,0)
		grid.addWidget(self.translatetool,3,0)
		grid.addWidget(self.ztranslate,4,0)
		grid.addWidget(self.rotatetool,5,0)
		grid.addWidget(self.scaletool,6,0)
		grid.addWidget(self.rulertool,7,0)
		grid.addWidget(self.linetool,8,0)
		grid.addWidget(self.cubetool,9,0)
		grid.addWidget(self.spheretool,10,0)
		grid.addWidget(self.cylindertool,11,0)
		grid.addWidget(self.conetool,12,0)
		grid.addWidget(self.texttool,13,0)
		grid.addWidget(self.datatool,14,0)
		grid.addWidget(self.apptool,15,0)
		grid.setAlignment(QtCore.Qt.AlignLeft)
		# labels
		#grid.addWidget(toollabel,0,0)
		grid.addWidget(self.selectiontool_label,1,1)
		grid.addWidget(self.multiselectiontool_label,2,1)
		grid.addWidget(self.translatetool_label,3,1)
		grid.addWidget(self.ztranslatetool_label,4,1)
		grid.addWidget(self.rotatetool_label,5,1)
		grid.addWidget(self.scaletool_label,6,1)
		grid.addWidget(self.rulertool_label,7,1)
		grid.addWidget(self.linetool_label,8,1)
		grid.addWidget(self.cubetool_label,9,1)
		grid.addWidget(self.spheretool_label,10,1)
		grid.addWidget(self.cylindertool_label,11,1)
		grid.addWidget(self.conetool_label,12,1)
		grid.addWidget(self.texttool_label,13,1)
		grid.addWidget(self.datatool_label,14,1)
		grid.addWidget(self.apptool_label,15,1)
		grid.setAlignment(QtCore.Qt.AlignLeft)	
		
		frame.setLayout(grid)
		hbox.addWidget(frame)
		
		self.stacked_widget = QtGui.QStackedWidget()
		self.stacked_widget.setFrameShape(QtGui.QFrame.StyledPanel)
		hbox.addWidget(self.stacked_widget)
		
		widget.setLayout(hbox) #grid
		
		return widget
	
	def getTreeWidget(self):
		"""
		This returns the treeview-control panel widget
		"""
		widget = QtGui.QWidget()
		hbox = QtGui.QHBoxLayout(widget)
		treeframe = QtGui.QFrame()
		treeframe.setFrameShape(QtGui.QFrame.StyledPanel)
		treeframe.setLayout(self._get_tree_layout(widget))
		treeframe.setMinimumWidth(self.mintreewidth)
		hbox.addWidget(treeframe)
		self.stacked_widget = QtGui.QStackedWidget()
		self.stacked_widget.setFrameShape(QtGui.QFrame.StyledPanel)
		hbox.addWidget(self.stacked_widget)
		widget.setLayout(hbox)
		
		return widget
	
	def _get_tree_layout(self, parent):
		"""
		Returns the tree layout
		"""
		tvbox = QtGui.QVBoxLayout()
		self.tree_widget = EMQTreeWidget(parent)
		self.tree_widget.setHeaderLabel("Choose a item")
		tvbox.addWidget(self.tree_widget)
		self.tree_node_button_add = QtGui.QPushButton("Add Object")
		self.tree_node_button_remove = QtGui.QPushButton("Remove Object")
		self.tree_node_slider = ValSlider(label="Seq:")
		self.tree_node_slider.setIntonly(True)
		self.tree_node_slider.setRange(0,1)
		self.tree_node_slider.setValue(0)
		tvbox.addWidget(self.tree_node_button_add)
		tvbox.addWidget(self.tree_node_button_remove)
		tvbox.addWidget(self.tree_node_slider)
		
# 		QtCore.QObject.connect(self.tree_widget, QtCore.SIGNAL("itemClicked(QTreeWidgetItem*,int)"), self._tree_widget_click)
# 		QtCore.QObject.connect(self.tree_widget, QtCore.SIGNAL("visibleItem(QTreeWidgetItem*)"), self._tree_widget_visible)
# 		QtCore.QObject.connect(self.tree_widget, QtCore.SIGNAL("editItem(QTreeWidgetItem*)"), self._tree_widget_edit)
# 		QtCore.QObject.connect(self.tree_node_button_remove, QtCore.SIGNAL("clicked()"), self._tree_widget_remove)
# 		QtCore.QObject.connect(self.tree_node_button_add, QtCore.SIGNAL("clicked()"), self._on_add_button)
# 		QtCore.QObject.connect(self.tree_node_slider, QtCore.SIGNAL("valueChanged"), self._slider_change)
		
		return tvbox

	
	def _recursiveupdatetreeselvis(self, item):
		item.setSelectionStateBox()
		item.getVisibleState()
		for childidx in xrange(item.childCount()):
			self._recursiveupdatetreeselvis(item.child(childidx))
			
	def updateTreeSelVis(self, selecteditem=None):
		"""
		Update the selection and visibility states. Makes the Sel Vis states an observer of the SG
		"""
		# Update the tree
		self._recursiveupdatetreeselvis(self.tree_widget.topLevelItem(0))
		# Set the desired item if desired
		if selecteditem:
			try:
				self.stacked_widget.setCurrentWidget(selecteditem.getItemInspector())
				self.tree_widget.setCurrentItem(selecteditem.EMQTreeWidgetItem)
				#if selecteditem: self.scenegraph().setCurrentSelection(selecteditem)
			except:
				pass
			# Unsure unqiue selection
			self.ensureUniqueTreeLevelSelection(selecteditem)
			
	def ensureUniqueTreeLevelSelection(self, item):
		"""
		Make sure that we don't select both an ancestor and child at the same time
		"""
		for ancestor in item.getSelectedAncestorNodes():
			if ancestor.EMQTreeWidgetItem:			# Not al ancestors are listed on the inspector tree (such as a data node)
				ancestor.EMQTreeWidgetItem.setSelectionState(False)
		for child in item.getAllSelectedNodes()[1:]: 	# Lop the node itself off
			child.EMQTreeWidgetItem.setSelectionState(False)
	
	def _slider_change(self):
		mdl=int(self.tree_node_slider.getValue())
# 		for i,child in enumerate(self.scenegraph().getChildren()):
# 			if i==mdl: child.setVisibleItem(True)
# 			else: child.setVisibleItem(False)
		
		#self.scenegraph().updateSG()

	def _recursiveAdd(self, parentitem, parentnode,depth=0):
		"""
		Helper function to laod the SG
		"""
		for child in parentnode.getChildren():
			if not child.getLabel(): child.setLabel(child.name)
			addeditem = self.addTreeNode(child.getLabel(), child, parentitem)
			self._recursiveAdd(addeditem, child,depth+1)
		# Expand the data items
		if parentitem.childCount() > 0: parentitem.setExpanded(True)
		self.tree_node_slider.setRange(0,len(parentnode.getChildren())-1)

		
	def loadSG(self):
		"""
		Load the SG
		"""
		pass
		#rootitem = self.addTreeNode("All Objects", self.scenegraph())
		#self._recursiveAdd(rootitem, self.scenegraph())
		
	def addTreeNode(self, name, item3d, parentitem=None, insertionindex=-1):
		"""
		Add a node (item3d) to the TreeWidget if not parent node, otherwise add a child to parent node
		We need to get a GUI for the treeitem. The treeitem and the GUI need know each other so they can talk
		The Treeitem also needs to know the node, so it can talk to the node.
		You can think of this as a three way conversation (the alterative it to use a mediator, but that is not worth it w/ only three players)
		"""
		tree_item = EMQTreeWidgetItem(QtCore.QStringList(name), item3d, parentitem)	# Make a QTreeItem widget, and let the TreeItem talk to the scenegraph node and its GUI
		item3d.setEMQTreeWidgetItem(tree_item)				# Reference to the EMQTreeWidgetItem
		item_inspector = item3d.getItemInspector()				# Get the node GUI controls 
		#return tree_item
		item_inspector.setInspector(self)					# Associate the item GUI with the inspector
		self.stacked_widget.addWidget(item_inspector)			# Add a widget to the stack
		item3d.setLabel(name)						# Set the label
		# Set icon status
		tree_item.setSelectionStateBox()
		# Set parent if one exists	
		if not parentitem:
			self.tree_widget.insertTopLevelItem(0, tree_item)
		else:
			if insertionindex >= 0:
				parentitem.insertChild(insertionindex, tree_item)
			else:
				parentitem.addChild(tree_item)
		return tree_item
	
	def removeTreeNode(self, parentitem, childindex):
		# I am using the parent item rather than the item itself b/c the stupid widget has no , remove self function...
		# Remove both the QTreeWidgetItem and the widget from the WidgetStack, otherwise we'll get memory leaks 
		if parentitem.child(childindex).item3d():
			self.stacked_widget.removeWidget(parentitem.child(childindex).item3d().getItemInspector())
		parentitem.takeChild(childindex)
	
	def clearTree(self):
		"""
		Clear the entire tree
		"""
		if self.tree_widget.topLevelItem(0):
			self.tree_widget.topLevelItem(0).removeAllChildren(self)
			self.tree_widget.takeTopLevelItem(0)
		
	def _tree_widget_click(self, item, col, quiet=False):
		"""
		When a user clicks on the selection tree check box
		"""
		self.stacked_widget.setCurrentWidget(item.item3d().getItemInspector())
		item.setSelectionState(item.checkState(0))
		# This code is to prevent both decendents and childer from being selected....
		if item.checkState(0) == QtCore.Qt.Checked: self.ensureUniqueTreeLevelSelection(item.item3d())
		if not item.item3d().isSelectedItem(): item.item3d().getItemInspector().updateItemControls() # This is too update a widget, translation and rotation may change in parent nodes change
		#self.scenegraph().setCurrentSelection(item.item3d())
		#if not quiet: self.updateSceneGraph()
		
	def _tree_widget_visible(self, item):
		"""
		When a user clicks on the visible icon
		"""
		item.toggleVisibleState()
	
	def _tree_widget_edit(self):
		"""
		When a use middle clicks
		"""
		nodedialog = NodeEditDialog(self, self.tree_widget.currentItem())
		nodedialog.exec_()
		self.activateWindow()
	
	def _on_add_button(self):
		nodedialog =  NodeDialog(self, self.tree_widget.currentItem())
		nodedialog.exec_()
		self.activateWindow()
		
	def _tree_widget_remove(self):
		"""
		When a use wants to remove a node_name
		"""
		item = self.tree_widget.currentItem()
		if item.parent:
			self.removeTreeNode(item.parent(), item.parent().indexOfChild(item)) 
			item.parent().item3d().removeChild(item.item3d())
			# In case we delete the currently selected item, we want6 to move selected item to last selection
# 			if self.scenegraph().getCurrentSelection() == item.item3d():
# 				self.scenegraph().setCurrentSelection(self.tree_widget.currentItem().item3d())
# 			self.updateSceneGraph()
		else:
			print "Error cannot remove root node!!"

	def getUtilsWidget(self):
		"""
		Return the utilites widget
		"""
		uwidget = QtGui.QWidget()
		uvbox = QtGui.QVBoxLayout()
		font = QtGui.QFont()
		font.setBold(True)
		
		self.opensession_button = QtGui.QPushButton("Open Session")
		self.savesession_button = QtGui.QPushButton("Save Session")
		self.savebutton = QtGui.QPushButton("Save Image Snapshot")
		
		self.open_tomogram_button = QtGui.QPushButton("Open Tomogram")
		self.open_segmentation_button = QtGui.QPushButton("Open Segmentation")
		self.save_segmentation_button = QtGui.QPushButton("Save Segmentation")
		
		uvbox.addWidget(self.opensession_button)
		uvbox.addWidget(self.savesession_button)
		uvbox.addWidget(self.savebutton)
		
		uvbox.addWidget(self.open_tomogram_button)
		uvbox.addWidget(self.open_segmentation_button)
		uvbox.addWidget(self.save_segmentation_button)
		uwidget.setLayout(uvbox)
		
		QtCore.QObject.connect(self.savebutton, QtCore.SIGNAL("clicked()"),self._on_save)
		QtCore.QObject.connect(self.savesession_button, QtCore.SIGNAL("clicked()"),self._on_save_session)
		QtCore.QObject.connect(self.opensession_button, QtCore.SIGNAL("clicked()"),self._on_open_session)
		
		QtCore.QObject.connect(self.open_tomogram_button, QtCore.SIGNAL("clicked()"),self._on_open_tomogram)
		QtCore.QObject.connect(self.open_segmentation_button, QtCore.SIGNAL("clicked()"),self._on_open_segmentation)
		QtCore.QObject.connect(self.save_segmentation_button, QtCore.SIGNAL("clicked()"),self._on_save_segmentation)
		
		return uwidget

	def _on_open_session(self):
		"""
		Open a session... (might want to add a warning dialog that this will close the current session)
		"""
		filename = QtGui.QFileDialog.getOpenFileName(self, 'Open Session', os.getcwd(), "*.eman")
		
	def _on_save_session(self):
		"""
		Return a list of all the child items (actually a tree of sorts)
		"""
		filename = QtGui.QFileDialog.getSaveFileName(self, 'Save Session', os.getcwd(), "*.eman")

	def _on_save(self):
		"""
		Save a snapshot of the scene
		"""
		filename = QtGui.QFileDialog.getSaveFileName(self, 'Save Image', os.getcwd(), "(*.tiff *.jpeg *.png)")

	def _on_open_tomogram(self):
		"""
		Open a session
		"""
		filename = QtGui.QFileDialog.getOpenFileName(self, 'Open Tomogram', os.getcwd(), "*.hdf,*.mrc")
		#if filename:
		#	self.scenegraph().loadSession(filename)
		
	def _on_open_segmentation(self):
		"""
		Open a session
		"""
		# Open the file
		filename = QtGui.QFileDialog.getOpenFileName(self, 'Open Segmentation', os.getcwd(), "*.hdf,*.xml")
		#if filename:
		#	self.scenegraph().loadSession(filename)
		
	def _on_save_segmentation(self):
		"""
		Save a snapshot of the scene
		"""
		filename = QtGui.QFileDialog.getSaveFileName(self, 'Save Segmentation', os.getcwd(), "(*.hdf,*.xml)")
		#if filename: # if we cancel
		#	self.scenegraph().saveSnapShot(filename)
		
	def updateInspector(self):
		"""
		Update Inspector,is called whenever the scence changes
		"""
		pass
	
	def updateTree(self, currentnode=None):
		"""
		Update the SG tree
		"""
		self.clearTree()
		self.loadSG()
		# either set the current node to the argument or set it to the one PM has selected
		if currentnode:
			self.tree_widget.setCurrentItem(currentnode.EMQTreeWidgetItem)
			idx = self.stacked_widget.indexOf(currentnode.getItemInspector())
			if idx >= 0: self.stacked_widget.setCurrentIndex(idx)
			self.scenegraph().setCurrentSelection(currentnode)
		else:
			node = self.scenegraph().getCurrentSelection()
			self.tree_widget.setCurrentItem(node.EMQTreeWidgetItem)
			idx = self.stacked_widget.indexOf(node.getItemInspector())
			if idx >= 0: self.stacked_widget.setCurrentIndex(idx)
			
	def updateSceneGraph(self):
		""" 
		Updates SG, in the near future this will be improved to allow for slow operations
		"""
		pass

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


		self.d3view = EMScene3D()
		self.d3viewdata = EMDataItem3D(test_image_3d(3), transform=Transform())
		isosurface = EMIsosurface(self.d3viewdata, transform=Transform())
		self.d3view.insertNewNode('', self.d3viewdata, parentnode=self.d3view)
		self.d3view.insertNewNode("Iso", isosurface, parentnode=self.d3viewdata )

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

			self.d3viewdata.setData(test_image_3d(3))
			self.d3view.updateSG()

			return

		if self.wfilt.getValue()!=0.0 :
			self.fdata=self.data.process("filter.lowpass.gauss",{"cutoff_freq":1.0/self.wfilt.getValue(),"apix":self.data['apix_x']}) #JESUS

		xyd=self.fdata.process("misc.directional_sum",{"axis":"z"})
		xzd=self.fdata.process("misc.directional_sum",{"axis":"y"})
		zyd=self.fdata.process("misc.directional_sum",{"axis":"x"})

		self.xyview.set_data(xyd)
		self.xzview.set_data(xzd)
		self.zyview.set_data(zyd)

		self.d3viewdata.setData(self.fdata)
		self.d3view.updateSG()


	def event_filter(self,value):
		self.update()

	def closeEvent(self, event):
		self.d3view.close()

if __name__ == '__main__':
	main()
