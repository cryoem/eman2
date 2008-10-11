#!/usr/bin/env python

#
# Author: David Woolford (woolford@bcm.edu)
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
#

import PyQt4
from PyQt4 import QtCore, QtGui, QtOpenGL
from PyQt4.QtCore import Qt
from emapplication import EMStandAloneApplication, EMQtWidgetModule

from emselector import EMSelectorDialog

from EMAN2 import EMData

class EMBrowserDialog(EMSelectorDialog):
	def __init__(self,target,application):
		EMSelectorDialog.__init__(self,target,application)

	def set_preview(self,filename):
		try:
			a=EMData.read_images(filename)
		except:
			return
		
		if len(a) == 1: a = a[0]
		
		import emimage
		
		if self.single_preview.isChecked():
			if self.gl_image_preview != None:
				self.application.close_specific(self.gl_image_preview)
				self.gl_image_preview == None
			
			#print self.gl_image_preview
			self.gl_image_preview =emimage.EMImageModule(a,None,False,self.application)
			#print self.gl_image_preview
				
			#self.gl_image_preview.set_data(a,filename)
			#self.gl_image_preview.set_file_name(f)
			self.application.show_specific(self.gl_image_preview)
			self.gl_image_preview.updateGL()
		else:
			preview = emimage.EMImageModule(a,None,False,self.application)
			preview.set_data(a,filename)
			self.application.show_specific(preview)
			preview.updateGL()
			
		self.application.setOverrideCursor(Qt.ArrowCursor)
	
	
	
app = None
def on_done(string_list):
	print "on done"
	if len(string_list) != 0:
		for s in string_list:
			print s,
		print
	app.quit()


if __name__ == '__main__':
	em_app = EMStandAloneApplication()
	app = em_app
	dialog = EMBrowserDialog(None,em_app)
	em_qt_widget = EMQtWidgetModule(dialog,em_app)
	QtCore.QObject.connect(dialog,QtCore.SIGNAL("done"),on_done)
	em_app.show()
	em_app.execute()

