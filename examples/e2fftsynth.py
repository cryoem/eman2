#!/usr/bin/env python

# Author: Steven Ludtke, 10/16/2010 (sludtke@bcm.edu)
# Copyright (c) 2000-2010 Baylor College of Medicine
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
from optparse import OptionParser

try:
	from PyQt4 import QtCore, QtGui, QtOpenGL
	from PyQt4.QtCore import Qt
	from emshape import *
	from valslider import ValSlider,ValBox
	from emimage import EMImageWidget
	from emimage2d import EMImage2DWidget
	from emplot2d import EMPlot2DWidget
	from emapplication import EMApp
except:
	print "Unable to import Python GUI libraries, PYTHONPATH ?"
	sys.exit(1)

def main():
	global debug,logid
	progname = os.path.basename(sys.argv[0])
	usage = """%prog [options]

This program allows the user to play around with Fourier synthesis graphically
	
"""

	parser = OptionParser(usage=usage,version=EMANVERSION)

#	parser.add_option("--gui",action="store_true",help="Start the GUI for interactive fitting",default=False)
	parser.add_option("--verbose", "-v", dest="verbose", action="store", metavar="n", type="int", default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	
	(options, args) = parser.parse_args()

	app=EMApp()
	win=GUIFourierSynth(app)
	win.show()
	try: 
		win.raise_()
		win.synthplot.raise_()
	except: pass
	app.exec_()
	
class GUIFourierSynth(QtGui.QWidget):
	"""This class represents an application for interactive Fourier synthesis"""
	
	def __init__(self,app):
		self.app=app
		QtGui.QWidget.__init__(self,None)

		self.synthplot=EMPlot2DWidget(self.app)
		self.synthplot.show()
		
		# overall layout
		self.vbl1=QtGui.QVBoxLayout()
		self.setLayout(self.vbl1)
		
		# First row contains general purpose controls
		self.hbl1=QtGui.QHBoxLayout()
		self.vbl1.addLayout(self.hbl1)
		
		self.vcell=ValBox(self,(0,128.0),"Cell:",64)
		self.hbl1.addWidget(self.vcell)
		
		self.vncells=ValBox(self,(0,128.0),"n Cells:",1)
		self.hbl1.addWidget(self.vncells)
		
		self.voversamp=ValBox(self,(0,128.0),"Oversample:",1)
		self.hbl1.addWidget(self.voversamp)
		
		self.targfn=None
		
		self.vnsin=ValBox(self,(1,64),"# Sin:",16)
		self.vnsin.intonly=1
		self.hbl1.addWidget(self.vnsin)
		
		self.cbshowall=QtGui.QCheckBox("Show All")
		self.hbl1.addWidget(self.cbshowall)

		self.cbshifted=QtGui.QCheckBox("Shifted")
		self.hbl1.addWidget(self.cbshifted)

		self.cbtargfn=QtGui.QComboBox(self)
		self.cbtargfn.addItem("None")
		self.cbtargfn.addItem("triangle")
		self.cbtargfn.addItem("square")
		self.cbtargfn.addItem("square imp")
		self.cbtargfn.addItem("delta")
		self.cbtargfn.addItem("noise")
		self.cbtargfn.addItem("saw")
		self.cbtargfn.addItem("sin")
		self.cbtargfn.addItem("modsin")
		self.cbtargfn.addItem("modsin2")
		self.cbtargfn.addItem("modsin3")
		self.cbtargfn.addItem("sin low")
		self.cbtargfn.addItem("doubledelta")
		self.hbl1.addWidget(self.cbtargfn)
		
		# Widget containing valsliders
		self.wapsliders=QtGui.QWidget(self)
#		self.wapsliders.setMinimumSize(800,640)
		self.gblap=QtGui.QGridLayout()
		self.gblap.setSizeConstraint(QtGui.QLayout.SetMinAndMaxSize)
		self.gblap.setColumnMinimumWidth(0,250)
		self.gblap.setColumnMinimumWidth(1,250)
		self.wapsliders.setLayout(self.gblap)
		
		# ScrollArea providing view on slider container widget
		self.wapsarea=QtGui.QScrollArea(self)
		self.wapsarea.setWidgetResizable(True)
		self.wapsarea.setWidget(self.wapsliders)
		self.vbl1.addWidget(self.wapsarea)

		QtCore.QObject.connect(self.vcell, QtCore.SIGNAL("valueChanged"), self.recompute)
		QtCore.QObject.connect(self.vncells, QtCore.SIGNAL("valueChanged"), self.recompute)
		QtCore.QObject.connect(self.voversamp, QtCore.SIGNAL("valueChanged"), self.recompute)
		QtCore.QObject.connect(self.vnsin, QtCore.SIGNAL("valueChanged"), self.nsinchange)
		QtCore.QObject.connect(self.cbshowall, QtCore.SIGNAL("stateChanged(int)"), self.recompute)
		QtCore.QObject.connect(self.cbshifted, QtCore.SIGNAL("stateChanged(int)"), self.recompute)
		QtCore.QObject.connect(self.cbtargfn,QtCore.SIGNAL("activated(int)"),self.newtargfn)


		self.wamp=[]
		self.wpha=[]
		self.curves=[]
		self.xvals=[]
		for i in range(65):
			self.wamp.append(ValSlider(self,(0.0,1.0),"%2d:"%i,0.0))
			self.gblap.addWidget(self.wamp[-1],i,0)
			QtCore.QObject.connect(self.wamp[-1], QtCore.SIGNAL("valueChanged"), self.recompute)
			
			self.wpha.append(ValSlider(self,(-180.0,180.0),"%2d:"%i,0.0))
			self.gblap.addWidget(self.wpha[-1],i,1)
			QtCore.QObject.connect(self.wpha[-1], QtCore.SIGNAL("valueChanged"), self.recompute)
		
			self.curves.append(EMData(64,1))
		
			if self.cbshowall.isChecked() :
				self.synthplot
		self.total=EMData(64,1)
	
		self.nsinchange()
		
	def newtargfn(self,index):
		"This function has a number of hardcoded 'target' functions"
	
		if index==0 : 
			self.targfn=None
		else :
			nx=int(self.vcell.getValue())
			self.targfn=EMData(nx,1)
			
			if index==1 : 	# triangle
				for i in xrange(nx/2):
					self.targfn[i]=-1.0+4.0*i/nx
				for i in xrange(nx/2,nx):
					self.targfn[i]=3.0-4.0*i/nx
			
			elif index==2 : # square
				for i in xrange(nx/4): self.targfn[i]=-1.0
				for i in xrange(nx/4,nx*3/4): self.targfn[i]=1.0
				for i in xrange(nx*3/4,nx): self.targfn[i]=-1.0
				
			elif index==3 : # square impulse
				self.targfn.to_zero()
				for i in xrange(nx/4,nx/2): self.targfn[i]=1.0
			
			elif index==4 : # delta
				self.targfn.to_zero()
				self.targfn[nx*2/5]=1.0
				
			elif index==5 : # noise
				self.targfn.process_inplace("testimage.noise.gauss",{"seed":0})
			
			elif index==6 : # saw
				self.targfn.to_zero()
				for i in xrange(nx/4,nx/2): self.targfn[i]=4.0*(i-nx/4.0)/nx
				for i in xrange(nx/2,nx*3/4): self.targfn[i]=-1+4.0*(i-nx/2.0)/nx
				
			elif index==7 : # sin
				for i in xrange(nx): self.targfn[i]=sin(i*pi/4.0)
			
			elif index==8 : # modulated sine
				for i in xrange(nx): self.targfn[i]=sin(i*pi/4.0)*sin(i*pi/32)

			elif index==9 : # modulated sine 2
				for i in xrange(nx): self.targfn[i]=sin(i*pi/4.0)*sin(i*pi/29)

			elif index==10 : # modulated sine 3
				for i in xrange(nx): self.targfn[i]=sin(i*pi/4.0)*sin(i*pi/126)

			elif index==11 : # sin low
				for i in xrange(nx): self.targfn[i]=sin(i*pi/16.0)

			elif index==12 : # double delta
				self.targfn.to_zero()
				self.targfn[nx/16]=4.0
				self.targfn[nx*15/16]=4.0

			self.target2sliders()
		
		self.recompute()

	def target2sliders(self):
		
		nsin=int(self.vnsin.getValue())
#		cp=self.targfn.process("xform.phaseorigin.tocenter")
		cp=self.targfn.copy()
		fft=cp.do_fft()
		fft[0]=fft[0]/2.0
		fft[fft["nx"]/2-1]=fft[fft["nx"]/2-1]/2.0
		fft.ri2ap()
		
		for i in xrange(min(fft["nx"]/2,nsin+1)):
#			print fft[i]
			amp=fft[i].real
			if fabs(amp)<1.0e-5 : amp=0.0
			self.wamp[i].setValue(amp*2/(fft["nx"]-2),quiet=1)
			self.wpha[i].setValue(fft[i].imag*180.0/pi+90.0,quiet=1)


	def nsinchange(self,value=None):
		if value==None : value=int(self.vnsin.getValue())
		
		for i in range(65):
			if i>value: 
				self.wamp[i].hide()
				self.wpha[i].hide()
			else :
				self.wamp[i].show()
				self.wpha[i].show()
		
		if self.targfn!=None :
			self.target2sliders()
			
		self.recompute()
	
	def recompute(self,value=None):
		nsin=int(self.vnsin.getValue())
		cell=int(self.vcell.getValue())
		ncells=self.vncells.getValue()
		oversamp=int(self.voversamp.getValue())
		samples=int(cell*ncells*oversamp)
		
		self.xvals=[xn/float(oversamp) for xn in range(samples)]
		self.total.set_size(samples)
		self.total.to_zero()
		for i in range(nsin+1):
			self.curves[i].set_size(samples)
			if i==0: 
				self.curves[i].to_one()
				if self.wpha[0].getValue()>180.0 : self.curves[i].mult(-1.0)
			else: self.curves[i].process_inplace("testimage.sinewave",{"wavelength":cell*oversamp/float(i),"phase":self.wpha[i].getValue()*pi/180.0})
			self.curves[i].mult(self.wamp[i].getValue())
			
			self.total.add(self.curves[i])
		
		self.synthplot.set_data((self.xvals,self.total.get_data_as_vector()),"Sum",replace=True,quiet=True,linewidth=2)
		
		if self.targfn!=None:
			self.synthplot.set_data(self.targfn,"Target",quiet=True,linewidth=1,linetype=2,symtype=0)
		
		if self.cbshowall.isChecked() :
			csum=self.total["minimum"]*1.1
			for i in range(nsin):
				if self.wamp[i].getValue()==0: continue
				csum-=self.wamp[i].getValue()*1.1
				if self.cbshifted.isChecked() : self.curves[i].add(csum)
				self.synthplot.set_data((self.xvals,self.curves[i].get_data_as_vector()),"%d"%i,quiet=True,linewidth=1,color=2)
				csum-=self.wamp[i].getValue()*1.1

		self.synthplot.updateGL()
		
if __name__ == "__main__":
	main()
	
