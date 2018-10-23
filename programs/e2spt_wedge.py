#!/usr/bin/env python
#
# Author: John Flanagan (jfflanag@bcm.edu); modified extensively by Jesus Galaz-Montoya (jgalaz@gmail.com)
# Copyright (c) 2012- Baylor College of Medicine
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
from __future__ import print_function
from __future__ import division
from past.utils import old_div
from builtins import range
from EMAN2 import *
from PyQt4 import QtCore, QtGui, QtOpenGL
from eman2_gui.emapplication import EMApp
from eman2_gui import emscene3d
from eman2_gui import emdataitem3d

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog 3Dstack [options]
	Visulaizse and compute the mean amplitude and sigma in the missing wedge region. After you are sasified that the missing wedge looks sane, compute missing wedge stats
	on all volumes. These stats are used by the aligner tomo.fsc, for subtomogram alignment and averaging.
	"""
	
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	
	parser.add_pos_argument(name="tdstack",help="The 3D stack to examine.", default="", guitype='filebox',  row=0, col=0,rowspan=1, colspan=2)
	parser.add_header(name="wedgeheader", help='Options below this label are specific to e2wedge', title="### e2wedge options ###", row=1, col=0, rowspan=1, colspan=2)
	parser.add_argument("--wedgeangle",type=float,help="Missing wedge angle",default=60.0, guitype='floatbox', row=2, col=0, rowspan=1, colspan=1)
	parser.add_argument("--wedgei",type=float,help="Missingwedge begining", default=0.05)
	parser.add_argument("--wedgef",type=float,help="Missingwedge ending", default=0.5)
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--nogui", action="store_true", default=False, help="Do not launch the GUI and set the average of the missing wedge statistics on all the volumes.")
	parser.add_argument("--averagestats", action="store_true", default=False, help="Do not launch the GUI and set the average of the missing wedge statistics on all the volumes.")

	(options, args) = parser.parse_args()

	logger = E2init(sys.argv, options.ppid)

	stack=args[0]
	
	if not options.nogui:	
		em_app = EMApp()
		wedgeviewer = MissingWedgeViewer(stack, options.wedgeangle, wedgei=options.wedgei, wedgef=options.wedgef)
		wedgeviewer.show()
		ret=em_app.exec_()
		sys.exit(ret)
	else:
		means=[]
		sigmas=[]
		
		n=EMUtil.get_image_count(stack)
		for i in range(n):
			a = EMData(stack,i)
			retr = wedgestats(a,options.wedgeangle,options.wedgei,options.wedgef)
			mean = retr[0]
			sigma = retr[1]
			
			if options.averagestats:
				means.append(mean)
				sigmas.append(sigma)
			else:
				a['spt_wedge_mean'] = mean
				a['spt_wedge_sigma'] = sigma
				print("The mean and sigma for subvolume %d are: mean=%f, sigma=%f" % (i,mean,sigma))
				a.write_image(stack,i)
		
		if options.averagestats:
			meanavg = old_div(sum(means),len(means))
			sigmaavg = old_div(sum(sigmas),len(sigmas))
			
			print("The average mean and sigma for the wedges in the stack are", meanavg, sigmaavg)
			for i in range(n):
				a = EMData(stack,i)
				a['spt_wedge_mean'] = meanavg
				a['spt_wedge_sigma'] = sigmaavg
				a.write_image(stack,i)
	
	E2end(logger)
	
	return
				
	
def wedgestats(volume,angle, wedgei, wedgef):
	vfft = volume.do_fft()
	wedge = vfft.getwedge(angle, wedgei, wedgef)		
	mean = vfft.get_attr('spt_wedge_mean')
	sigma = vfft.get_attr('spt_wedge_sigma')
	return(mean,sigma)


class MissingWedgeViewer(QtWidgets.QWidget):
	""" Display a missing wedge"""
	def __init__(self, filename, angle, wedgei=0.0, wedgef=1.0):
		QtWidgets.QWidget.__init__(self)
		self.setWindowTitle('The Wedge Viewer')
		self.setMinimumWidth(400)
		self.setMinimumHeight(400)
		self.filename = filename
		self.volumes = EMData.read_images(self.filename)
		self.angle = angle
		self.wedgedata = None
		self.wedgemean = None
		self.wedgesigma = None
		
		self.dataitems= []
		
		grid=QtWidgets.QGridLayout()
		# Make threed widget
		self.widget3d=emscene3d.EMScene3D()
		grid.addWidget(self.widget3d, 0, 0, 1, 4)
		# make contols
		combolabel = QtWidgets.QLabel("Volume Idx")
		self.volcombobox = QtWidgets.QComboBox()
		self.fitbutton = QtWidgets.QPushButton("Fit Wedge")
		self.fitbutton.setEnabled(False)
		labeli = QtWidgets.QLabel("Wedge_i")
		self.wedgei = QtWidgets.QLineEdit(str(wedgei))
		labelf = QtWidgets.QLabel("Wedge_f")
		self.wedgef = QtWidgets.QLineEdit(str(wedgef))
		wmlabel = QtWidgets.QLabel("Wedge Mean")
		wslabel = QtWidgets.QLabel("Wedge Sigma")
		self.wedgemeanwidget = QtWidgets.QLineEdit("0.0")
		self.wedgemeanwidget.setReadOnly(True)
		self.wedgesigmawidget = QtWidgets.QLineEdit("0.0")
		self.wedgesigmawidget.setReadOnly(True)
		self.setwedgestats = QtWidgets.QPushButton("Set Wedege stats, This Vol")
		self.setallwedgestats = QtWidgets.QPushButton("Set Wedege stats, All Vols")
		grid.addWidget(combolabel, 1, 0)
		grid.addWidget(self.volcombobox, 1, 1)
		grid.addWidget(self.fitbutton, 1 ,2, 1, 2)
		grid.addWidget(labeli, 2, 0)
		grid.addWidget(self.wedgei, 2, 1)
		grid.addWidget(labelf, 2, 2)
		grid.addWidget(self.wedgef, 2, 3)
		grid.addWidget(wmlabel, 3, 0)
		grid.addWidget(self.wedgemeanwidget, 3, 1)
		grid.addWidget(wslabel, 3, 2)
		grid.addWidget(self.wedgesigmawidget, 3, 3)
		grid.addWidget(self.setwedgestats, 4, 0, 1, 2) 
		grid.addWidget(self.setallwedgestats, 4, 2, 1, 2)
		
		self.setLayout(grid)
		
		# Fill combox box
		for i in range(len(self.volumes)):
			self.volcombobox.addItem(str(i))
			
		#Add actions
		self.volcombobox.activated[int].connect(self.onVolChange)
		self.wedgei.editingFinished.connect(self.onWedgeChange)
		self.wedgef.editingFinished.connect(self.onWedgeChange)
		self.setwedgestats.clicked.connect(self.onOneVolStats)
		self.setallwedgestats.clicked.connect(self.onManyVolStats)
		
		# Use first volume
		self.onVolChange(0)
		
	
	def onVolChange(self, idx):
		""" Compute wedge on new volume """
		children = tuple(self.widget3d.getChildren())
		for child in children:
			self.widget3d.removeChild(child)
		# Make objects
		wedge = self.makeWedge(self.volumes[idx])
		self.addData(wedge,"Wedge")
		self.findWedge(self.volumes[idx], self.angle)
		self.widget3d.updateSG()
		self.widget3d.updateTree()
	
	def onWedgeChange(self):
		""" Recompute wedge stats """
		if not self.wedgedata: return
		self.widget3d.removeChild(self.wedgedata)
		self.findWedge(self.volumes[int(self.volcombobox.currentIndex())], self.angle)
		self.widget3d.updateSG()
		self.widget3d.updateTree()
	
	def onOneVolStats(self):
		idx = int(self.volcombobox.currentIndex())
		self.volumes[idx].set_attr('spt_wedge_mean',self.wedgemean)
		self.volumes[idx].set_attr('spt_wedge_sigma',self.wedgesigma)
		self.volumes[idx].write_image(self.filename, idx, EMUtil.ImageType.IMAGE_UNKNOWN, True)
		
	def onManyVolStats(self):
		for idx, volume in enumerate(self.volumes):
			volume.set_attr('spt_wedge_mean',self.wedgemean)
			volume.set_attr('spt_wedge_sigma',self.wedgesigma)
			volume.write_image(self.filename, idx, EMUtil.ImageType.IMAGE_UNKNOWN, True)
		
	def addData(self, wedge, name, color=None):
		""" Add Data """
		data=emdataitem3d.EMDataItem3D(wedge)
		self.widget3d.insertNewNode(name, data, parentnode=self.widget3d)
		isosurface = emdataitem3d.EMIsosurface(data)
		self.widget3d.insertNewNode("Iso", isosurface, parentnode=data)
		if color:
			isosurface.setAmbientColor(*color)
			isosurface.setDiffuseColor(*color)
		# add to list
		self.dataitems.append([data, isosurface])
		
		return data
		
		
	def makeWedge(self, volume):
		""" Compute the amps and center """
		vfft = volume.do_fft()
		vfft.ri2ap()
		amps = vfft.amplitude()
		amps.process_inplace('xform.phaseorigin.tocenter')
		symamps = amps.process('xform.mirror', {'axis':'x'})
		return (amps + symamps)
		
	def findWedge(self, volume, angle):
		""" Compute missing wedge stats, and draw missing wedge to check sanity """
		vfft = volume.do_fft()
		wedge = vfft.getwedge(angle, float(self.wedgei.text()), float(self.wedgef.text()))
		
		#vfft.set_attr('is_complex',0)
		#self.addData(vfft,"wedge")
		#self.addData(wedge,"wedegy")
		
		wedge.process_inplace('xform',{'transform':Transform({'type':'eman','tx':old_div(wedge.get_xsize(),2), 'ty':old_div(wedge.get_ysize(),2),'tz':old_div(wedge.get_zsize(),2)})})
		w2 = wedge.process('xform.mirror',{'axis':'x'}) + wedge
		wholewedge = w2.process('xform.mirror',{'axis':'y'}) + w2
		self.wedgedata = self.addData(wholewedge,"Computed Wedge", [0.0625,0.8555,0.9453])
		# Update info
		self.wedgemean = vfft.get_attr('spt_wedge_mean')
		self.wedgesigma = vfft.get_attr('spt_wedge_sigma')
		self.wedgemeanwidget.setText(str(round(self.wedgemean,2)))
		self.wedgesigmawidget.setText(str(round(self.wedgesigma,2)))
	
if __name__ == "__main__":
	main()
