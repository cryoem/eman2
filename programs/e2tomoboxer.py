#!/usr/bin/env python

#
# Author: Steven Ludtke  2/8/2011 (rewritten)
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

from PyQt4 import QtCore, QtGui
from PyQt4.QtCore import Qt
from optparse import OptionParser
import sys
import os

from EMAN2 import *
from emapplication import get_application, EMApp
from emimage2d import EMImage2DWidget
from emimage3d import EMImage3DWidget
from emshape import EMShape
from valslider import *

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """%prog [options] <Volume file>

	WARNING: This program still under development.
	
	Tomography 3-D particle picker and annotation tool. Still under development."""

	parser = OptionParser(usage=usage,version=EMANVERSION)

#	parser.add_option("--boxsize","-B",type="int",help="Box size in pixels",default=64)
#	parser.add_option("--shrink",type="int",help="Shrink factor for full-frame view, default=0 (auto)",default=0)
	parser.add_option("--inmemory",action="store_true",default=False,help="This will read the entire tomogram into memory. Much faster, but you must have enough ram !")
	parser.add_option("--verbose", "-v", dest="verbose", action="store", metavar="n", type="int", default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	
	(options, args) = parser.parse_args()
		
#	if len(args) < 0: parser.error("You must sepcify the input argument")
#	if not file_exists(args[0]): parser.error("%s does not exist" %args[0])
#	if options.boxsize < 2: parser.error("The boxsize you specified is too small")
#	# The program will not run very rapidly at such large box sizes anyhow
#	if options.boxsize > 2048: parser.error("The boxsize you specified is too large.\nCurrently there is a hard coded max which is 2048.\nPlease contact developers if this is a problem.")
	
#	logid=E2init(sys.argv)
	
	app = EMApp()
	if options.inmemory : 
		print "Reading tomogram. Please wait."
		img=EMData(args[0],0)
		print "Done !"
		boxer=EMTomoBoxer(data=img)
	else : boxer=EMTomoBoxer(datafile=args[0])
	boxer.show()
	app.execute()
#	E2end(logid)

class EMBoxViewer(QtGui.QWidget):
	"""This is a multi-paned view showing a single boxed out particle from a larger tomogram"""
	
	def __init__(self):
		QtGui.QWidget.__init__(self)
		
		self.gbl = QtGui.QGridLayout(self)
		self.xyview = EMImage2DWidget(noparent=True)
		self.gbl.addWidget(self.xyview,1,0)

		self.xzview = EMImage2DWidget(noparent=True)
		self.gbl.addWidget(self.xzview,0,0)

		self.zyview = EMImage2DWidget(noparent=True)
		self.gbl.addWidget(self.zyview,1,1)

		self.d3view = EMImage3DWidget(noparent=True)
		self.gbl.addWidget(self.d3view,0,1)

		QtCore.QObject.connect(self.xyview,QtCore.SIGNAL("mousedown"),self.xy_down)
		QtCore.QObject.connect(self.xyview,QtCore.SIGNAL("mousedrag"),self.xy_drag)
		QtCore.QObject.connect(self.xyview,QtCore.SIGNAL("mouseup")  ,self.xy_up  )

		QtCore.QObject.connect(self.xzview,QtCore.SIGNAL("mousedown"),self.xz_down)
		QtCore.QObject.connect(self.xzview,QtCore.SIGNAL("mousedrag"),self.xz_drag)
		QtCore.QObject.connect(self.xzview,QtCore.SIGNAL("mouseup")  ,self.xz_up  )

		QtCore.QObject.connect(self.zyview,QtCore.SIGNAL("mousedown"),self.zy_down)
		QtCore.QObject.connect(self.zyview,QtCore.SIGNAL("mousedrag"),self.zy_drag)
		QtCore.QObject.connect(self.zyview,QtCore.SIGNAL("mouseup")  ,self.zy_up  )

#		self.setSizeGripEnabled(True)

#		if get_platform() == "Darwin": # because OpenGL widgets in Qt don't leave room in the bottom right hand corner for the resize tool
#			self.status = QtGui.QStatusBar()
#			self.gbl.addWidget(self.status,3,0,1,2)
#			self.margin = 0
				
	def set_data(self,data,pos,sz):
		"""Sets the current volume to display as well as the center position (x,y,z)  and a size (int) for a 'box' in that region"""
		
		self.boxsize=sz
		self.data=data
		self.update_pos(pos)
		
		self.show()
		
	def update_pos(self,pos):
		self.curpos=pos
		
		xyd=self.data.process("misc.directional_sum",{"axis":"z","first":max(self.curpos[2]-self.boxsize/2,0),"last":self.curpos[2]+self.boxsize/2})
		xzd=self.data.process("misc.directional_sum",{"axis":"y","first":max(self.curpos[1]-self.boxsize/2,0),"last":self.curpos[1]+self.boxsize/2})
		zyd=self.data.process("misc.directional_sum",{"axis":"x","first":max(self.curpos[0]-self.boxsize/2,0),"last":self.curpos[0]+self.boxsize/2})
		d3d=self.data.get_clip(Region(pos[0]-self.boxsize/2,pos[1]-self.boxsize/2,pos[2]-self.boxsize/2,self.boxsize,self.boxsize,self.boxsize))
		
		self.xyview.set_data(xyd)
		self.xyview.add_shape("box",EMShape(("rect",.7,.7,.1,self.curpos[0]-self.boxsize/2,self.curpos[1]-self.boxsize/2,self.curpos[0]+self.boxsize/2,self.curpos[1]+self.boxsize/2,2)))
		
		self.xzview.set_data(xzd)
		self.xzview.add_shape("box",EMShape(("rect",.7,.7,.1,self.curpos[0]-self.boxsize/2,self.curpos[2]-self.boxsize/2,self.curpos[0]+self.boxsize/2,self.curpos[2]+self.boxsize/2,2)))

		self.zyview.set_data(zyd)
		self.zyview.add_shape("box",EMShape(("rect",.7,.7,.1,self.curpos[2]-self.boxsize/2,self.curpos[1]-self.boxsize/2,self.curpos[2]+self.boxsize/2,self.curpos[1]+self.boxsize/2,2)))

		self.d3view.set_data(d3d)
		
	def xy_down(self,me):
		pass
	
	def xy_drag(self,me):
		pass
		
	def xy_up  (self,me):
		pass
		
	def xz_down(self,me):
		pass
		
	def xz_drag(self,me):
		pass
		
	def xz_up  (self,me):
		pass
		
	def zy_down(self,me):
		pass
		
	def zy_drag(self,me):
		pass
		
	def zy_up  (self,me):
		pass
		

class EMTomoBoxer(QtGui.QMainWindow):
	"""This class represents the EMTomoBoxer application instance.  """
	
	def __init__(self,data=None,datafile=None):
		QtGui.QWidget.__init__(self)
		
		# Menu Bar
		self.mfile=self.menuBar().addMenu("File")
		self.mfile_open=self.mfile.addAction("Open")
		self.mfile_quit=self.mfile.addAction("Quit")
		
		self.setCentralWidget(QtGui.QWidget())
		self.gbl = QtGui.QGridLayout(self.centralWidget())

		# relative stretch factors
		self.gbl.setColumnStretch(0,1)
		self.gbl.setColumnStretch(1,5)
		self.gbl.setColumnStretch(2,0)
		self.gbl.setRowStretch(0,1)
		self.gbl.setRowStretch(1,5)
		
		# 3 orthogonal restricted projection views
		self.xyview = EMImage2DWidget()
		self.gbl.addWidget(self.xyview,1,1)

		self.xzview = EMImage2DWidget()
		self.gbl.addWidget(self.xzview,0,1)

		self.zyview = EMImage2DWidget()
		self.gbl.addWidget(self.zyview,1,0)

		# Select layers to integrate for x-y view. For xz and yz, boxsize is used.
		self.wdepth = QtGui.QSlider()
		self.gbl.addWidget(self.wdepth,0,2)
		
		### Control panel area in upper left corner
		self.gbl2 = QtGui.QGridLayout()
		self.gbl.addLayout(self.gbl2,0,0)
		
		# box size
		self.wboxsize=ValBox(label="Box Size:",value=64)
		self.gbl2.addWidget(self.wboxsize,0,0,1,2)
		
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
		
		# scale factor
		self.wzoom=ValSlider(rng=(.1,2),label="Sca:",value=1.0)
		self.gbl2.addWidget(self.wzoom,2,0,1,2)
		
		
		
		self.curbox=-1
		self.boxes=[]						# array of box info, each is (x,y,z,...)
		self.xydown=None
		
		QtCore.QObject.connect(self.wdepth,QtCore.SIGNAL("valueChanged(int)"),self.set_depth)
		QtCore.QObject.connect(self.wnlayers,QtCore.SIGNAL("valueChanged(int)"),self.set_nlayers)
		QtCore.QObject.connect(self.wboxsize,QtCore.SIGNAL("valueChanged"),self.set_boxsize)
		QtCore.QObject.connect(self.wmaxmean,QtCore.SIGNAL("clicked(bool)"),self.set_projmode)
		QtCore.QObject.connect(self.wzoom,QtCore.SIGNAL("valueChanged")  ,self.set_scale  )
		QtCore.QObject.connect(self.mfile_open,QtCore.SIGNAL("triggered(bool)")  ,self.menu_file_open  )


		QtCore.QObject.connect(self.xyview,QtCore.SIGNAL("mousedown"),self.xy_down)
		QtCore.QObject.connect(self.xyview,QtCore.SIGNAL("mousedrag"),self.xy_drag)
		QtCore.QObject.connect(self.xyview,QtCore.SIGNAL("mouseup"),self.xy_up  )
		QtCore.QObject.connect(self.xyview,QtCore.SIGNAL("mousewheel"),self.xy_wheel  )
		QtCore.QObject.connect(self.xyview,QtCore.SIGNAL("set_scale"),self.xy_scale)
		QtCore.QObject.connect(self.xyview,QtCore.SIGNAL("set_origin"),self.xy_origin)

		QtCore.QObject.connect(self.xzview,QtCore.SIGNAL("mousedown"),self.xz_down)
		QtCore.QObject.connect(self.xzview,QtCore.SIGNAL("mousedrag"),self.xz_drag)
		QtCore.QObject.connect(self.xzview,QtCore.SIGNAL("mouseup")  ,self.xz_up  )

		QtCore.QObject.connect(self.zyview,QtCore.SIGNAL("mousedown"),self.zy_down)
		QtCore.QObject.connect(self.zyview,QtCore.SIGNAL("mousedrag"),self.zy_drag)
		QtCore.QObject.connect(self.zyview,QtCore.SIGNAL("mouseup")  ,self.zy_up  )

		if datafile!=None : self.set_datafile(datafile)		# This triggers a lot of things to happen, so we do it last
		if data!=None : self.set_data(data)

	def set_datafile(self,datafile):
		if datafile==None :
			self.datafile=None
			self.data=None
			self.xyview.set_data(None)
			self.xzview.set_data(None)
			self.zyview.set_data(None)
			return
		
		self.data=None
		self.datafile=datafile
		imgh=EMData(datafile,0,1)
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
		self.datasize=(data["nx"],data["ny"],data["nz"])
		self.wdepth.setRange(0,self.datasize[2]-1)
		self.boxes=[]
		self.curbox=-1
		
		self.wdepth.setValue(self.datasize[2]/2)
		self.update_all()

	def get_slice(self,n,xyz):
		"""Reads a slice either from a file or the preloaded memory array. xyz is the axis along which 'n' runs, 0=x (yz), 1=y (xz), 2=z (xy)"""
		if self.data!=None :
			if xyz==0:
				r=self.data.get_clip(Region(n,0,0,1,self.datasize[1],self.datasize[2]))
				r.set_size(self.datasize[1],self.datasize[2],1)
				return r
			elif xyz==1:
				r=self.data.get_clip(Region(0,n,0,self.datasize[0],1,self.datasize[2]))
				r.set_size(self.datasize[0],self.datasize[2],1)
				return r
			else:
				r=self.data.get_clip(Region(0,0,n,self.datasize[0],self.datasize[1],1))
				return r

		elif self.datafile!=None:
			if xyz==0:
				r=EMData()
				r.read_image(self.datafile,0,0,Region(n,0,0,1,self.datasize[1],self.datasize[2]))
				r.set_size(self.datasize[1],self.datasize[2],1)
				return r
			elif xyz==1:
				r=EMData()
				r.read_image(self.datafile,0,0,Region(0,n,0,self.datasize[0],1,self.datasize[2]))
				r.set_size(self.datasize[0],self.datasize[2],1)
				return r
			else:
				r=EMData()
				r.read_image(self.datafile,0,0,Region(0,0,n,self.datasize[0],self.datasize[1],1))
				return r



		else : return None

	def set_boxsize(self):
		pass
	
	def set_projmode(self,state):
		"""Projection mode can be simple average (state=False) or maximum projection (state=True)"""
		self.update_all()
	
	def set_scale(self):
		pass
	
	def set_depth(self):
		self.update_xy()
		
	def set_nlayers(self):
		self.update_all()

	def boxsize(self): return int(self.wboxsize.getValue())
	
	def nlayers(self): return int(self.wnlayers.value())
	
	def depth(self): return int(self.wdepth.value())

	def menu_file_open(self,tog):
		print("file open")

	def get_averager(self):
		"""returns an averager of the appropriate type for generating projection views"""
		if self.wmaxmean.isChecked() : return Averagers.get("minmax",{"max":1})
		
		return Averagers.get("mean")

	def update_sides(self):
		"""updates xz and yz views due to a new center location"""
		
		if self.datafile==None and self.data==None : return
		
		if self.curbox==-1 :
			x=self.datasize[0]/2
			y=self.datasize[1]/2
			z=0
		else:
			x,y,z=self.boxes[self.curbox][:3]
		
		self.cury=y
		self.curx=x
		
		# yz
		avgr=self.get_averager()

		for x in range(x-self.nlayers()/2,x+(self.nlayers()+1)/2):
			slc=self.get_slice(x,0)
			avgr.add_image(slc)
			
		av=avgr.finish()
		av.process_inplace("xform.transpose")
		self.zyview.set_data(av)
		
		# xz
		avgr=self.get_averager()

		for y in range(y-self.nlayers()/2,y+(self.nlayers()+1)/2):
			slc=self.get_slice(y,1)
			avgr.add_image(slc)
			
		self.xzview.set_data(avgr.finish())


	def update_xy(self):
		"""updates xy view due to a new slice range"""
		
		if self.datafile==None and self.data==None: return
		
		if self.wmaxmean.isChecked() : avgr=Averagers.get("minmax",{"max":1})
		else : avgr=Averagers.get("mean")

		slc=EMData()
		for z in range(self.wdepth.value()-self.nlayers()/2,self.wdepth.value()+(self.nlayers()+1)/2):
			slc=self.get_slice(z,2)
			avgr.add_image(slc)
			
		self.xyview.set_data(avgr.finish())

	def update_all(self):
		"""redisplay of all widgets"""
		if self.datafile==None and self.data==None: return
		
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

	def update_box(self,n):
		"""After adjusting a box, call this"""
		self.curbox=n
		box=self.boxes[n]
		bs2=self.boxsize()/2
		
		# Boxes may not extend outside the tomogram
		if box[0]<bs2 : box[0]=bs2
		if box[0]>self.datasize[0]-bs2 : box[0]=self.datasize[0]-bs2
		if box[1]<bs2 : box[1]=bs2
		if box[1]>self.datasize[1]-bs2 : box[1]=self.datasize[1]-bs2
		if box[2]<bs2 : box[2]=bs2
		if box[2]>self.datasize[2]-bs2 : box[2]=self.datasize[2]-bs2
#		print self.boxes
		self.xyview.add_shape(n,EMShape(("rect",.2,.2,.8,box[0]-bs2,box[1]-bs2,box[0]+bs2,box[1]+bs2,1)))
		self.xyview.add_shape("xl",EMShape(("line",.8,.8,.1,0,box[1],self.datasize[0],box[1],1)))
		self.xyview.add_shape("yl",EMShape(("line",.8,.8,.1,box[0],0,box[0],self.datasize[1],1)))
		self.xzview.add_shape(n,EMShape(("rect",.2,.2,.8,box[0]-bs2,box[2]-bs2,box[0]+bs2,box[2]+bs2,1)))
		self.xzview.add_shape("zl",EMShape(("line",.8,.8,.1,box[0],0,box[0],self.datasize[2],1)))
		self.zyview.add_shape(n,EMShape(("rect",.2,.2,.8,box[2]-bs2,box[1]-bs2,box[2]+bs2,box[1]+bs2,1)))
		self.zyview.add_shape("zl",EMShape(("line",.8,.8,.1,0,box[1],self.datasize[2],box[1],1)))
		
		if self.depth()!=box[2] : self.wdepth.setValue(box[2])
		else : self.xyview.update()
		self.update_sides()

	def xy_down(self,event):
		x,y=self.xyview.scr_to_img((event.x(),event.y()))
		x,y=int(x),int(y)
		self.xydown=None
		if x<0 or y<0 : return		# no clicking outside the image (on 2 sides)
		
		for i in range(len(self.boxes)):
			if self.inside_box(i,x,y) :
				self.xydown=(i,x,y,self.boxes[i][0],self.boxes[i][1])
				self.update_box(i)
				break
		else:
#			if x>self.boxsize()/2 and x<self.datasize[0]-self.boxsize()/2 and y>self.boxsize()/2 and y<self.datasize[1]-self.boxsize()/2 and self.depth()>self.boxsize()/2 and self.depth()<self.datasize[2]-self.boxsize()/2 :
			self.boxes.append(([x,y,self.depth()]))
			self.xydown=(len(self.boxes)-1,x,y,x,y)		# box #, x down, y down, x box at down, y box at down
			self.update_box(self.xydown[0])
		
	def xy_drag(self,event):
		if self.xydown==None : return
		
		x,y=self.xyview.scr_to_img((event.x(),event.y()))
		x,y=int(x),int(y)
		
		dx=x-self.xydown[1]
		dy=y-self.xydown[2]
		self.boxes[self.xydown[0]][0]=dx+self.xydown[3]
		self.boxes[self.xydown[0]][1]=dy+self.xydown[4]
		self.update_box(self.curbox)
		
	def xy_up  (self,event):
		self.xydown=None
	
	def xy_wheel (self,event):
		if event.delta() > 0:
			self.wdepth.setValue(self.wdepth.value()+4)
		elif event.delta() < 0:
			self.wdepth.setValue(self.wdepth.value()-4)
	
	def xy_scale(self,news):
		"Image view has been rescaled"
		self.xzview.set_scale(news,True)
		self.zyview.set_scale(news,True)
		
	def xy_origin(self,newor):
		xzo=self.xzview.get_origin()
		xzo.set_origin(newor[0],xzo[1],True)
		
		zyo=self.zyview.get_origin()
		zyo.set_origin(zyo[0],newor[1],True)
	
	def xz_down(self,event):
		x,z=self.xzview.scr_to_img((event.x(),event.y()))
		x,z=int(x),int(z)
		self.xzdown=None
		if x<0 or z<0 : return		# no clicking outside the image (on 2 sides)
		
		for i in range(len(self.boxes)):
			if self.inside_box(i,x,-1,z) :
				self.xzdown=(i,x,z,self.boxes[i][0],self.boxes[i][2])
				self.update_box(i)
				break
		else:
			self.boxes.append(([x,self.cury,z]))
			self.xzdown=(len(self.boxes)-1,x,z,x,z)		# box #, x down, y down, x box at down, y box at down
			self.update_box(self.xzdown[0])
		
	def xz_drag(self,event):
		if self.xzdown==None : return
		
		x,z=self.xzview.scr_to_img((event.x(),event.y()))
		x,z=int(x),int(z)
		
		dx=x-self.xzdown[1]
		dz=z-self.xzdown[2]
		self.boxes[self.xzdown[0]][0]=dx+self.xzdown[3]
		self.boxes[self.xzdown[0]][2]=dz+self.xzdown[4]
		self.update_box(self.curbox)
		
	def xz_up  (self,event):
		self.xzdown=None

	def zy_down(self,event):
		z,y=self.zyview.scr_to_img((event.x(),event.y()))
		z,y=int(z),int(y)
		self.xydown=None
		if z<0 or y<0 : return		# no clicking outside the image (on 2 sides)
		
		for i in range(len(self.boxes)):
			if self.inside_box(i,-1,y,z) :
				self.zydown=(i,z,y,self.boxes[i][2],self.boxes[i][1])
				self.update_box(i)
				break
		else:
			self.boxes.append(([self.curx,y,z]))
			self.zydown=(len(self.boxes)-1,z,y,z,y)		# box #, x down, y down, x box at down, y box at down
			self.update_box(self.zydown[0])
		
	def zy_drag(self,event):
		if self.zydown==None : return
		
		z,y=self.zyview.scr_to_img((event.x(),event.y()))
		z,y=int(z),int(y)
		
		dz=z-self.zydown[1]
		dy=y-self.xydown[2]
		self.boxes[self.zydown[0]][2]=dx+self.xydown[3]
		self.boxes[self.zydown[0]][1]=dy+self.xydown[4]
		self.update_box(self.curbox)
		
	def zy_up  (self,event):
		self.zydown=None

		
	#def closeEvent(self,event):
		#self.target().done()


if __name__ == "__main__":
	main()
		
		
		
