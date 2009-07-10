#!/usr/bin/env python
#
# Author: David Woolford (woolford@bcm.edu)
# Copyright (c) 2000-2006 Baylor College of Medicine


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

from emboxerbase import *

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """%prog [options] <image> <image2>....

a refactoring of e2boxer
"""

	parser = OptionParser(usage=usage,version=EMANVERSION)
	parser.add_option("--boxsize","-B",type="int",help="Box size in pixels",default=128)
	
	(options, args) = parser.parse_args()
	
	error_message = check(options,args)
	if len(error_message) > 0:
		error = "\n"
		for e in error_message:
			error += "Error: "+e +"\n"
		parser.error(error)
		
	args = [abs_path(arg) for arg in args] # always try to use full file names 

	application = EMStandAloneApplication()
#	QtCore.QObject.connect(gui, QtCore.SIGNAL("module_idle"), on_idle)
	
	module = EMBoxerModule(args,options.boxsize)
	module.show_interfaces()
	# this is an example of how to add your own custom tools:
	module.add_2d_window_mouse_tool(SwarmEventHandling,SwarmPanel,particle_diameter=options.boxsize)
	application.execute()


class SwarmPanel:
	def __init__(self,target,particle_diameter=128):
		self.busy = True
		self.particle_diameter = particle_diameter
		self.target = weakref.ref(target)
		self.widget = None
		self.busy = False
		
	def icon(self):
		from PyQt4 import QtGui
		return QtGui.QIcon(get_image_directory() + "swarm_icon.png")
		
	def get_widget(self):
		if self.widget == None:
			from PyQt4 import QtCore, QtGui, Qt
			self.widget = QtGui.QWidget()
			vbl = QtGui.QVBoxLayout(self.widget)
			vbl.setMargin(0)
			vbl.setSpacing(6)
			vbl.setObjectName("vbl")
			
			hbl = QtGui.QHBoxLayout()
			hbl.addWidget(QtGui.QLabel("Particle Diameter:"))
			
			self.ptcl_diam_edit = QtGui.QLineEdit(str(self.particle_diameter))
			hbl.addWidget(self.ptcl_diam_edit)
			
			self.update_template = QtGui.QCheckBox("Update Template")
			self.update_template.setChecked(True)
			
			vbl.addLayout(hbl)
			vbl.addWidget(self.update_template)
			QtCore.QObject.connect(self.ptcl_diam_edit,QtCore.SIGNAL("editingFinished()"),self.new_ptcl_diam)
			QtCore.QObject.connect(self.update_template,QtCore.SIGNAL("clicked(bool)"),self.update_template_checked)
			
		return self.widget
	
	def new_ptcl_diam(self):
		if self.busy: return
		print "yo"
		
	def update_template_checked(self,val):
		if self.busy: return
		print "yo"
		
	def hide(self):
		if self.widget != None:
			self.widget.hide()
			
class SwarmEventHandling:
	'''
	A class that knows how to handle mouse erase events for a GUIBox
	'''
	
	def __init__(self,target,panel_object=None,particle_diameter=128):
		self.target = weakref.ref(target)
		self.particle_diameter = particle_diameter
		
	def unique_name(self): return "Swarm"
	
	def set_panel_object(self,panel): self.panel_object = panel
		
	def mouse_move(self,event):
		pass
		
	def mouse_wheel(self,event):
		pass
	def mouse_down(self,event) :
		pass
	def mouse_drag(self,event) :
		pass
	def mouse_up(self,event) :
		pass

if __name__ == "__main__":
	main()
