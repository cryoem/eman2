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
	parser.add_option("--shrink",type="int",help="Shrink factor for full-frame view, default=0 (auto)",default=0)
	parser.add_option("--verbose", "-v", dest="verbose", action="store", metavar="n", type="int", default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	
	(options, args) = parser.parse_args()
		
#	if len(args) < 0: parser.error("You must sepcify the input argument")
#	if not file_exists(args[0]): parser.error("%s does not exist" %args[0])
#	if options.boxsize < 2: parser.error("The boxsize you specified is too small")
#	# The program will not run very rapidly at such large box sizes anyhow
#	if options.boxsize > 2048: parser.error("The boxsize you specified is too large.\nCurrently there is a hard coded max which is 2048.\nPlease contact developers if this is a problem.")
	
#	logid=E2init(sys.argv)
	
	app = EMApp()
	boxer=EMTomoBoxer(args[0])
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
	
	def __init__(self,datafile=None):
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
		self.gbl.addWidget(self.xzview,1,0)

		self.zyview = EMImage2DWidget()
		self.gbl.addWidget(self.zyview,0,1)

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
		self.boxes=[]						# array of box info [(x,y,z), ]
		
		QtCore.QObject.connect(self.wdepth,QtCore.SIGNAL("valueChanged(int)"),self.set_depth)
		QtCore.QObject.connect(self.wnlayers,QtCore.SIGNAL("valueChanged(int)"),self.set_nlayers)
		QtCore.QObject.connect(self.wboxsize,QtCore.SIGNAL("valueChanged"),self.set_boxsize)
		QtCore.QObject.connect(self.wmaxmean,QtCore.SIGNAL("clicked(bool)"),self.set_projmode)
		QtCore.QObject.connect(self.wzoom,QtCore.SIGNAL("valueChanged")  ,self.set_scale  )
		QtCore.QObject.connect(self.mfile_open,QtCore.SIGNAL("triggered(bool)")  ,self.menu_file_open  )


		QtCore.QObject.connect(self.xyview,QtCore.SIGNAL("mousedown"),self.xy_down)
		QtCore.QObject.connect(self.xyview,QtCore.SIGNAL("mousedrag"),self.xy_drag)
		QtCore.QObject.connect(self.xyview,QtCore.SIGNAL("mouseup")  ,self.xy_up  )

		#QtCore.QObject.connect(self.xzview,QtCore.SIGNAL("mousedown"),self.xz_down)
		#QtCore.QObject.connect(self.xzview,QtCore.SIGNAL("mousedrag"),self.xz_drag)
		#QtCore.QObject.connect(self.xzview,QtCore.SIGNAL("mouseup")  ,self.xz_up  )

		#QtCore.QObject.connect(self.zyview,QtCore.SIGNAL("mousedown"),self.zy_down)
		#QtCore.QObject.connect(self.zyview,QtCore.SIGNAL("mousedrag"),self.zy_drag)
		#QtCore.QObject.connect(self.zyview,QtCore.SIGNAL("mouseup")  ,self.zy_up  )

		self.set_data(datafile)		# This triggers a lot of things to happen, so we do it last

	def set_data(self,datafile):
		if datafile==None :
			self.datafile=None
			self.xyview.set_data(None)
			self.xzview.set_data(None)
			self.zyview.set_data(None)
			return
		
		self.datafile=datafile
		imgh=EMData(datafile,0,1)
		self.datasize=(imgh["nx"],imgh["ny"],imgh["nz"])
		self.wdepth.setRange(0,self.datasize[2]-1)
		self.boxes=[]
		self.curbox=-1
		
		self.update_all()
		
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
		self.update_xy()
	
	def menu_file_open(self,tog):
		print("file open")

	def get_averager(self):
		"""returns an averager of the appropriate type for generating projection views"""
		if self.wmaxmean.isChecked() : return Averagers.get("minmax",{"max":1})
		
		return Averagers.get("mean")

	def update_sides(self):
		"""updates xz and yz views due to a new center location"""
		
		if self.datafile==None : return
		
		if self.curbox==-1 :
			x=self.datasize[0]/2
			y=self.datasize[1]/2
			z=0
		else:
			x,y,z=self.boxes[self.curbox][0]
		
		
		slc=EMData()
		
		# yz
		avgr=self.get_averager()

		for x in range(x-self.nlayers()/2,x+(self.nlayers()+1)/2):
			slc.read_image(self.datafile,0,0,Region(x,0,0,1,self.datasize[1],self.datasize[2]))
			avgr.add_image(slc)
			
		self.zyview.set_data(avgr.finish())
		
		# xz
		avgr=self.get_averager()

		for y in range(y-self.nlayers()/2,y+(self.nlayers()+1)/2):
			slc.read_image(self.datafile,0,0,Region(0,y,0,self.datasize[0],1,self.datasize[2]))
			avgr.add_image(slc)
			
		self.xzview.set_data(avgr.finish())


	def update_xy(self):
		"""updates xy view due to a new slice range"""
		
		if self.datafile==None: return
		
		if self.wmaxmean.isChecked() : avgr=Averagers.get("minmax",{"max":1})
		else : avgr=Averagers.get("mean")

		slc=EMData()
		for z in range(self.wdepth.value()-self.nlayers()/2,self.wdepth.value()+(self.nlayers()+1)/2):
			slc.read_image(self.datafile,0,0,Region(0,0,z,self.datasize[0],self.datasize[1],1))
			avgr.add_image(slc)
			
		self.xyview.set_data(avgr.finish())

	def update_all(self):
		"""redisplay of all widgets"""
		if self.datafile==None: return
		
		self.update_xy()
		self.update_sides()
		
		#self.xyview.update()
		#self.xzview.update()
		#self.zyview.update()

	def boxsize(self): return self.wboxsize.getValue()
	
	def nlayers(self): return self.wnlayers.value()
	
	def depth(self): return self.wdepth.value()

	def inside_box(self,n,x=-1,y=-1,z=-1):
		"""Checks to see if a point in image coordinates is inside box number n. If any value is negative, it will not be checked."""
		if x>=0 and (x<self.boxes[n][0]-self.boxsize()/2 or x>self.boxes[n][0]+self.boxsize()/2 : return False
		if y>=0 and (y<self.boxes[n][1]-self.boxsize()/2 or y>self.boxes[n][1]+self.boxsize()/2 : return False
		if z>=0 and (z<self.boxes[n][2]-self.boxsize()/2 or z>self.boxes[n][2]+self.boxsize()/2 : return False
		return True

	def xy_down(self,me):
		x,y=self.xyview.scr_to_img((event.x(),event.y()))
		if x<0 or y<0 : return		# no clicking outside the image (on 2 sides)
		
		for i in self.boxes:
			if self.inside_box(i,x,y) :
				self.xydown=(i,x,y)
		elif x>self.boxsize()/2 and x<self.datasize[0]-self.boxsize()/2 and y>self.boxsize()/2 and y<self.datasize[1]-self.boxsize()/2 :
			self.boxes.append((x,y,self.depth()))
			self.xydown(len(self.boxes)-1,x,y))
		else: self.xydown=None
		
	def xy_drag(self,me):
		if self.xydown==None : return
		
		x,y=self.xyview.scr_to_img((event.x(),event.y()))
		dx=x-self.xydown[1]
		dy=y-self.xydown[2]
		self.update_pos((self.curpos[0]+dx,self.curpos[1]+dy,self.curpos[2]))
		
	def xy_up  (self,me):
		self.xydown=None
		
	def xz_down(self,me):
		x,y=self.parent.scr_to_img((event.x(),event.y()))
		if me.x()>self.curpos[0]-self.boxsize/2 and me.x()<self.curpos[0]+self.boxsize/2 and me.y()>self.curpos[2]-self.boxsize/2 and me.y()<self.curpos[2]+self.boxsize/2:
			self.xzdown=(me.x(),me.y())
		else : self.xzdown=None
		
	def xz_drag(self,me):
		x,y=self.parent.scr_to_img((event.x(),event.y()))
		if self.xzdown==None : return
		dx=me.x()-self.xzdown[0]
		dy=me.y()-self.xzdown[1]
		self.update_pos((self.curpos[0]+dx,self.curpos[1],self.curpos[2]+dy))
		
	def xz_up  (self,me):
		self.xzdown=None
		
	def zy_down(self,me):
		x,y=self.parent.scr_to_img((event.x(),event.y()))
		if me.x()>self.curpos[2]-self.boxsize/2 and me.x()<self.curpos[2]+self.boxsize/2 and me.y()>self.curpos[1]-self.boxsize/2 and me.y()<self.curpos[1]+self.boxsize/2:
			self.zydown=(me.x(),me.y())
		else : self.zydown=None
		
	def zy_drag(self,me):
		x,y=self.parent.scr_to_img((event.x(),event.y()))
		if self.zydown==None : return
		dx=me.x()-self.zydown[0]
		dy=me.y()-self.zydown[1]
		self.update_pos((self.curpos[0],self.curpos[1]+dy,self.curpos[2]+dx))
		
	def zy_up  (self,me):
		self.zydown=None
		
		
	#def closeEvent(self,event):
		#self.target().done()


if __name__ == "__main__":
	main()
		
		
		
