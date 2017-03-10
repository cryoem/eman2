#!/usr/bin/env python

#
# Author: Steven Ludtke, 03/03/2017 (sludtke@bcm.edu)
# Copyright (c) 2000-2017 Baylor College of Medicine
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

from EMAN2 import *
from PyQt4 import QtCore
from emapplication import EMApp
from emimage2d import EMImage2DWidget

def main():
	global cur_data,data,imdisp
	
	# an application
	em_app = EMApp()

	widget=TestDisplay()
	widget.show()

	em_app.execute()

class TestDisplay(EMImage2DWidget):
	def __init__(self):
		EMImage2DWidget.__init__(self)
		# the single image display widget
		self.datatodisp=[test_image(1),test_image(2)]

		self.curdata=0
		self.set_data(self.datatodisp[0])
	
		timer = QtCore.QTimer(self)
		self.connect(timer,QtCore.SIGNAL("timeout()"),self.mytimeout)
		timer.start(1000)

	def mytimeout(self):
		self.curdata=self.curdata^1

		self.set_data(self.datatodisp[self.curdata])
		
main()
