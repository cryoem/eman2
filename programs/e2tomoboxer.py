#!/usr/bin/env python

#
# Author: Steven Ludtke
# Copyright (c) 2009-2010 Baylor College of Medicine
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

from optparse import OptionParser
import sys
import os
from EMAN2 import file_exists, gimme_image_dimensions3D,EMANVERSION,EMData,Region,Transform,get_image_directory,db_check_dict,db_open_dict,db_close_dict
from EMAN2 import get_file_tag,check_eman2_type,Processors,parsemodopt,E2init,E2end,E2progress,get_platform
from emapplication import get_application, EMStandAloneApplication,EMQtWidgetModule
from emimage2d import EMImage2DModule
from emimage3d import EMImage3DModule
import weakref
from emshape import EMShape
from PyQt4 import QtCore, QtGui
from PyQt4.QtCore import Qt
import math
import EMAN2db
import time

tomo_db_name = "bdb:e2tomoboxercache#"

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """%prog [options] <Volume file>
	
Tomography 3-D particle picker and annotation tool. Still under development."""

	parser = OptionParser(usage=usage,version=EMANVERSION)

#	parser.add_option("--boxsize","-B",type="int",help="Box size in pixels",default=64)
	parser.add_option("--shrink",type="int",help="Shrink factor for full-frame view, default=0 (auto)",default=0)
	parser.add_option("--verbose", "-v", dest="verbose", action="store", metavar="n", type="int", default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	
	(options, args) = parser.parse_args()
	
	error = check(args,options)
	
	if error:
		sys.exit(1)
#	if len(args) < 0: parser.error("You must sepcify the input argument")
#	if not file_exists(args[0]): parser.error("%s does not exist" %args[0])
#	if options.boxsize < 2: parser.error("The boxsize you specified is too small")
#	# The program will not run very rapidly at such large box sizes anyhow
#	if options.boxsize > 2048: parser.error("The boxsize you specified is too large.\nCurrently there is a hard coded max which is 2048.\nPlease contact developers if this is a problem.")
	
	logid=E2init(sys.argv)
	
	app = EMStandAloneApplication()
	boxer=EMTomoBoxer(args[0])
	boxer.show()
	app.execute()
	E2end(logid)

class EMBoxViewer(QtGui.QWidget):
	"""This is a multi-paned view showing a single boxed out particle from a larger tomogram"""
	
	def __init__(self):
		QtGui.QWidget.__init__(self)
		
		self.gbl = QtGui.QGridLayout(self)
		self.xyview = EMImage2DModule(noparent=True)
		self.gbl.addWidget(self.xyview.get_qt_widget(),1,0)

		self.xzview = EMImage2DModule(noparent=True)
		self.gbl.addWidget(self.xzview.get_qt_widget(),0,0)

		self.zyview = EMImage2DModule(noparent=True)
		self.gbl.addWidget(self.zyview.get_qt_widget(),1,1)

		self.d3view = EMImage3DModule(noparent=True)
		self.gbl.addWidget(self.d3view.get_qt_widget(),0,1)

		QtCore.QObject.connect(self.xyview.emitter(),QtCore.SIGNAL("mousedown"),self.xy_down)
		QtCore.QObject.connect(self.xyview.emitter(),QtCore.SIGNAL("mousedrag"),self.xy_drag)
		QtCore.QObject.connect(self.xyview.emitter(),QtCore.SIGNAL("mouseup")  ,self.xy_up  )

		QtCore.QObject.connect(self.xzview.emitter(),QtCore.SIGNAL("mousedown"),self.xz_down)
		QtCore.QObject.connect(self.xzview.emitter(),QtCore.SIGNAL("mousedrag"),self.xz_drag)
		QtCore.QObject.connect(self.xzview.emitter(),QtCore.SIGNAL("mouseup")  ,self.xz_up  )

		QtCore.QObject.connect(self.zyview.emitter(),QtCore.SIGNAL("mousedown"),self.zy_down)
		QtCore.QObject.connect(self.zyview.emitter(),QtCore.SIGNAL("mousedrag"),self.zy_drag)
		QtCore.QObject.connect(self.zyview.emitter(),QtCore.SIGNAL("mouseup")  ,self.zy_up  )

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
		x,y=self.main_2d_window.scr_to_img((event.x(),event.y()))
		if x>self.curpos[0]-self.boxsize/2 and x<self.curpos[0]+self.boxsize/2 and y>self.curpos[1]-self.boxsize/2 and y<self.curpos[1]+self.boxsize/2:
			self.xydown=(x,y)
		else : self.xydown=None
		
	def xy_drag(self,me):
		if self.xydown==None : return
		x,y=self.main_2d_window.scr_to_img((event.x(),event.y()))
		dx=x-self.xydown[0]
		dy=y-self.xydown[1]
		self.update_pos((self.curpos[0]+dx,self.curpos[1]+dy,self.curpos[2]))
		
	def xy_up  (self,me):
		self.xydown=None
		
	def xz_down(self,me):
		x,y=self.main_2d_window.scr_to_img((event.x(),event.y()))
		if me.x()>self.curpos[0]-self.boxsize/2 and me.x()<self.curpos[0]+self.boxsize/2 and me.y()>self.curpos[2]-self.boxsize/2 and me.y()<self.curpos[2]+self.boxsize/2:
			self.xzdown=(me.x(),me.y())
		else : self.xzdown=None
		
	def xz_drag(self,me):
		x,y=self.main_2d_window.scr_to_img((event.x(),event.y()))
		if self.xzdown==None : return
		dx=me.x()-self.xzdown[0]
		dy=me.y()-self.xzdown[1]
		self.update_pos((self.curpos[0]+dx,self.curpos[1],self.curpos[2]+dy))
		
	def xz_up  (self,me):
		self.xzdown=None
		
	def zy_down(self,me):
		x,y=self.main_2d_window.scr_to_img((event.x(),event.y()))
		if me.x()>self.curpos[2]-self.boxsize/2 and me.x()<self.curpos[2]+self.boxsize/2 and me.y()>self.curpos[1]-self.boxsize/2 and me.y()<self.curpos[1]+self.boxsize/2:
			self.zydown=(me.x(),me.y())
		else : self.zydown=None
		
	def zy_drag(self,me):
		x,y=self.main_2d_window.scr_to_img((event.x(),event.y()))
		if self.zydown==None : return
		dx=me.x()-self.zydown[0]
		dy=me.y()-self.zydown[1]
		self.update_pos((self.curpos[0],self.curpos[1]+dy,self.curpos[2]+dx))
		
	def zy_up  (self,me):
		self.zydown=None
		

class EMTomoBoxer(QtGui.QWidget):
	"""This class represents the EMTomoBoxer application instance. It is a GUI-only object which will make
	use of appropriate module-level functions for specific operations. """
	
	def __init__(self,targfile) :
		QtGui.QWidget.__init__(self)

		self.targfile=targfile
		hdr=EMData(targfile,0,1)
		self.targnx=hdr["nx"]
		self.targny=hdr["ny"]
		self.targnz=hdr["nz"]
		self.apix=hdr["apix_x"]
		
		QtGui.QWidget.__init__(self,None)
		self.setWindowIcon(QtGui.QIcon(get_image_directory() +"eman.png"))
		
		self.vbl = QtGui.QVBoxLayout(self)
		self.vbl.setMargin(2)
		self.vbl.setSpacing(6)
		self.vbl.setObjectName("vbl")
		
		self.thickvs=ValSlider(self,(1,self.targnz),"Thickness:")
		self.thickvs.setValue(16)
		self.vbl.addWidget(self.thickvs)
		
		self.zvs=ValSlider(self,(0,self.targnz),"Z:")
		self.zvs.setValue(self.targnz/2)
		self.vbl.addWidget(self.thickvs)
		
		# This is the tab-bar for mouse mode selection
		self.modetab = QtGui.QTabWidget()
		
		# Box tab
		self.boxtab = QtGui.QWidget()
		self.boxlay = QtGui.QGridLayout(self.boxtab)
		
		self.bsl = QtGui.QHBoxLayout()
		self.boxlay.addLayout(self.bsl,1,0,0,4)
		self.bslbl = QtGui.QLabel("Box Size :")
		self.bssb = QtGui.QSpinBox()
		self.bssb.setRange(8,512)
		self.bssb.setValue(48)
		
		self.bsl.addWidget(self.bslbl)
		self.bsl.addWidget(self.bssb)
		
		
		self.mmtab.addTab(self.boxtab,"Box")
		QtCore.QObject.connect(self.bssb,QtCore.SIGNAL("valueChanged(int)"), self.set_box_size)

		# Filament tab
		self.filtab = QtGui.QWidget()
		self.fillay = QtGui.QGridLayout(self.filtab)

		self.mmtab.addTab(self.boxtab,"Filament")
		
		# Ellipsoid tab
		self.filtab = QtGui.QWidget()
		self.fillay = QtGui.QGridLayout(self.filtab)

		self.mmtab.addTab(self.boxtab,"Ellipsoid")
		
		self.vbl.addWidget(self.modetab)
		# Overall widget connections
		QtCore.QObject.connect(self.modetab, QtCore.SIGNAL("currentChanged(int)"), self.set_mouse_mode)
#		self.connect(self.box_size,QtCore.SIGNAL("editingFinished()"),self.new_box_size)
#		self.connect(self.gen_output_but,QtCore.SIGNAL("clicked(bool)"),self.target().run_output_dialog)
#		self.connect(self.done_but,QtCore.SIGNAL("clicked(bool)"),self.done)

		self.mainview=EMImage2DModule()
		self.curslice=EMData(targfile,0,0,Region(0,0,self.targnz/2,self.targnx,self.targny,1))
		self.mainview.set_data(self.curslice)
		self.mainview.show()
		
		self.boxview=EMBoxViewer()


	def set_mouse_mode(self,mode):
		self.mmode=mode

	def set_box_size(self):
		box_size=int(self.box_size.text())
		self.target().set_box_size(box_size)	
		
	def closeEvent(self,event):
		self.target().done()


if __name__ == "__main__":
	main()
		
		
		