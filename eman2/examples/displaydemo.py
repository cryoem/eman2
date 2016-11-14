#!/usr/bin/env python

#
# Author: Steven Ludtke, 04/10/2003 (sludtke@bcm.edu)
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
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  2111-1307 USA
#
#

# displaydemo.py  09/18/2009  Steven Ludtke


from EMAN2 import *
from math import *
from PyQt4 import QtCore
from emapplication import EMApp
from emimage2d import EMImage2DWidget
from emshape import EMShape


def main():
	# an application
	em_app = EMApp()
	control1=TestControl(em_app)
	control2=TestControl(em_app)

	em_app.execute()

class TestControl():
	def __init__(self,app):
		# the single image display widget
		self.im2d = EMImage2DWidget(application=app)
	
		# get some signals from the window.
		QtCore.QObject.connect(self.im2d,QtCore.SIGNAL("mousedown"),self.down)
		QtCore.QObject.connedisplaydemo.pyct(self.im2d,QtCore.SIGNAL("mousedrag"),self.drag)
		QtCore.QObject.connect(self.im2d,QtCore.SIGNAL("mouseup"),self.up)
	
		#self explanatory
		a=test_image(size=(512,512))
		self.im2d.set_data(a)
		self.im2d.show()


	def down(self,event,lc):
		"""The event contains the x,y coordinates in window space, lc are the coordinates in image space"""

		self.downloc=lc	

	def drag(self,event,lc):
		s=EMShape(["line",0,.7,0,self.downloc[0],self.downloc[1],lc[0],lc[1],1])
		self.im2d.add_shape("mine",s)
		self.im2d.updateGL()
	
	def up(self,event,lc):
		s=EMShape(["line",.7,.2,0,self.downloc[0],self.downloc[1],lc[0],lc[1],1])
		self.im2d.del_shape("mine")
		self.im2d.add_shape("done",s)
		self.im2d.updateGL()
	
if __name__ == "__main__":  main()
