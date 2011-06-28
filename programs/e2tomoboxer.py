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

from sys import argv #jesus

from PyQt4 import QtCore, QtGui
from PyQt4.QtCore import Qt
from optparse import OptionParser
import sys
import os
import weakref
import threading

from EMAN2 import *
from emapplication import get_application, EMApp
from emimage2d import EMImage2DWidget
from emimagemx import EMImageMXWidget
from emimage3d import EMImage3DWidget
from emshape import EMShape
from valslider import *


def unbinned_extractor(boxsize,x,y,z,cbin,tomogram=argv[1]):

	boxsize=boxsize*cbin
	x=x*cbin
	y=y*cbin
	z=z*cbin

	r = Region((2*x-boxsize)/2,(2*y-boxsize)/2, (2*z-boxsize)/2, boxsize, boxsize, boxsize)
	e = EMData()
	e.read_image(tomogram,0,False,r)
	return(e)

def commandline_tomoboxer(tomogram,coordinates,subset,boxsize,cbin,output,output_format,swapyz,reverse_contrast):

	clines = open(coordinates,'r').readlines()
	set = len(clines)
	
	if subset:
		if subset > set:
			print "WARNING: The total amount of lines in the coordinates files is LESS that the subset of particles to box you specified; therefore, ALL particles will be extracted"
		else:
			set=subset
	
	print "The size of the set of sub-volumes to extract is", set
	for i in range(set):
	
		#Some people might manually make ABERANT coordinates files with commas, tabs, or more than once space in between coordinates
    	
	       	clines[i] = clines[i].replace(", ",' ')	
		clines[i] = clines[i].replace(",",' ')
		clines[i] = clines[i].replace("x",'')
		clines[i] = clines[i].replace("y",'')
		clines[i] = clines[i].replace("z",'')
		clines[i] = clines[i].replace("=",'')
        	clines[i] = clines[i].replace("_",' ')
		clines[i] = clines[i].replace("\n",' ')
		clines[i] = clines[i].replace("\t",' ')
		clines[i] = clines[i].replace("  ",' ')
		clines[i] = clines[i].split()		
		
		x = int(clines[i][0])
		y = int(clines[i][1])
		z = int(clines[i][2])

		if cbin > 1:
	                x = x * cbin
  			y = y * cbin
			z = z * cbin		
		
		#if cbinx > 1:
		#	x = int(clines[i][0]) * cbinx
		#if cbiny > 1:
		#	y = int(clines[i][0]) * cbiny
		#if cbinz > 1:
		#	z = int(clines[i][0]) * cbinz

		if swapyz:
			print '''You indicated that Y and Z coordinates are flipped in the coordinates file, respect to the orientation of the tomogram; 
			therefore, they will be swapped'''
			aux = y
			y = z
			z = aux

		print "The coordinates for particle#%d are x=%d, y=%d, z=%d " % (i,x,y,z)

		#r = Region((2*x-boxsize)/2,(2*y-boxsize)/2, (2*z-boxsize)/2, boxsize, boxsize, boxsize)
        	#e = EMData()
		#e.read_image(tomogram,0,False,r)
		
		e=unbinned_extractor(boxsize,x,y,z,cbin)

		#IF the boxed out particle is NOT empty, perform BASIC RAW-PARTICLE EDITING: contrast reversak and normalization 
		#Sometimes empty boxes are picked when boxing from the commandline if yshort isn't specified but should have, 
		#or if erroneous binning factors are provided
		
		if e['mean'] != 0:
			
			#It IS CONVENIENT to record any processing done on the particles as header parameters
			#to easily track what has been done to them. All relevant EMAN2 processing steps on sub-volumes are recorded as parameters starting with e2spt
			
			if reverse_contrast:
				e = e*(-1)
				e['e2spt_contrast_reversed'] = 'yes'
			else:
				e['e2spt_contrast_reversed'] = 'no'

			e['e2spt_tomogram'] = tomogram
			e['e2spt_coordx'] = x
			e['e2spt_coordy'] = y
			e['e2spt_coordz'] = z
			
			#The origin WILL most likely be MESSED UP if you don't explicitely set it to ZERO.
			#This can create ANNOYING visualization problems in Chimera
			
			e['origin_x'] = 0
			e['origin_y'] = 0
			e['origin_z'] = 0

			e = e.process('normalize.edgemean')
			e['e2spt_normalize.edgemean'] = 'yes'

			e['xform.align3d'] = Transform({"type":'eman','az':0,'alt':0,'phi':0,'tx':0,'ty':0,'tz':0})
			
			k=i
			if output_format == 'single':
				name = output.split('.')[0] + '_' + str(i).zfill(3) + '.' + output.split('.')[1]
				k=0
			else:
				name = output
			
			e.write_image(name,k)
			
		else:
			print """WARNING! The particle was skipped (and not boxed) because it's mean was ZERO (which often indicates a box is empty).
				Your coordinates file and/or the binning factors specified might be MESSED UP. 
				Also, make sure to specify --yshort if the 'ice thickness' (that is, the shortest legth of your tomogram) lies along the Y direction"""	

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """%prog [options] <Volume file>

	WARNING: This program still under development.
	
	Tomography 3-D particle picker and annotation tool. Still under development."""

	parser = OptionParser(usage=usage,version=EMANVERSION)

	parser.add_option("--boxsize","-B",type="int",help="Box size in pixels",default=32)
#	parser.add_option("--shrink",type="int",help="Shrink factor for full-frame view, default=0 (auto)",default=0)
	parser.add_option("--inmemory",action="store_true",default=False,help="This will read the entire tomogram into memory. Much faster, but you must have enough ram !")
	parser.add_option("--yshort",action="store_true",default=False,help="This means you have a file where y is the short axis")
	parser.add_option("--apix",type="float",help="Override the A/pix value stored in the tomogram header",default=0.0)
	parser.add_option("--verbose", "-v", dest="verbose", action="store", metavar="n", type="int", default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	
	parser.add_option('--bin', type='int', default=1, help="""Specify the binning/shrinking factor you want to use (for X,Y and Z) when opening the tomogram for boxing. 
								Don't worry, the sub-volumes will be extracted from the UNBINNED tomogram. 
								If binx, biny or binz are also specified, they will override the general bin value for the corresponding X, Y or Z directions""")
	
	#parser.add_option('--binx', type='int', default=1, help="Specify the binning/shrinking factor to use in X when opening the tomogram for boxing. Don't worry, the sub-volumes will be extracted from the UNBINNED tomogram")
	#parser.add_option('--biny', type='int', default=1, help="Specify the binning/shrinking factor to use in Y when opening the tomogram for boxing. Don't worry, the sub-volumes will be extracted from the UNBINNED tomogram")
	#parser.add_option('--binz', type='int', default=1, help="Specify the binning/shrinking factor to use in Z when opening the tomogram for boxing. Don't worry, the sub-volumes will be extracted from the UNBINNED tomogram")
	
	parser.add_option('--reverse_contrast', action="store_true", default=False, help='''This means you want the contrast to me inverted while boxing, and for the extracted sub-volumes.
											Remember that EMAN2 **MUST** work with "white" protein. You can very easily figure out what the original color
											of the protein is in your data by looking at the gold fiducials or the edge of the carbon hole in your tomogram.
											If they look black you MUST specify this option''')
	
	#parameters for commandline boxer
	
	parser.add_option('--coords', type='str', default='', help='Provide a coordinates file that contains the center coordinates of the sub-volumes you want to extract, to box from the command line')
	
	parser.add_option('--cbin', type='int', default=1, help='''Specifies the scale of the coordinates respect to the actual size of the tomogram where you want to extract the particles from.
								For example, provide 2 if you recorded the coordinates from a tomogram that was binned by 2, 
								but want to extract the sub-volumes from the UNbinned version of that tomogram''')
	
	#parser.add_option('--cbinx', type='int', default=1, help="""Binning factor of the X coordinates with respect to the actual size of the tomogram from which you want to extract the subvolumes.
	#							Sometimes tomograms are not binned equally in all directions for purposes of recording the coordinates of particles""")

	#parser.add_option('--cbiny', type='int', default=1, help="Binning factor of the Y coordinates with respect to the actual size of the tomogram from which you want to extract the subvolumes")
	#parser.add_option('--cbinz', type='int', default=1, help="Binning factor of the X coordinates with respect to the actual size of the tomogram from which you want to extract the subvolumes")

	parser.add_option('--subset', type='int', default=0, help='''Specify how many sub-volumes from the coordinates file you want to extract; e.g, if you specify 10, the first 10 particles will be boxed.
								0 means "box them all" because it makes no sense to box none''')
	parser.add_option('--output', type='str', default='stack.hdf', help="Specify the name of the stack file where to write the extracted sub-volumes")
	parser.add_option('--output_format', type='str', default='stack', help='''Specify 'single' if you want the sub-volumes to be written to individual files. You MUST still provide an output name in the regular way.
										For example, if you specify --output=myparticles.hdf\n
										but also specify --output_format=single\n
										then the particles will be written as individual files named myparticles_000.hdf myparticles_001.hdf...etc''')
	
	parser.add_option('--swapyz', action="store_true", default=False, help='''This means that the coordinates file and the actual tomogram do not agree regarding which is the "short" direction.
										For example, the coordinates file migh thave a line like this:\n
										1243 3412 45\n
										where clearly the "short" direction is Z; yet, if in the actual tomogram the short direction is Y, as they come out from
										IMOD by default, then the line should have been:\n
										1243 45 3412\n''')
	(options, args) = parser.parse_args()
		
	if len(args) != 1: 
		parser.error("You must specify a single volume data file on the command-line.")		
	
	if not file_exists(args[0]): parser.error("%s does not exist" %args[0])

#	if options.boxsize < 2: parser.error("The boxsize you specified is too small")
#	# The program will not run very rapidly at such large box sizes anyhow
#	if options.boxsize > 2048: parser.error("The boxsize you specified is too large.\nCurrently there is a hard coded max which is 2048.\nPlease contact developers if this is a problem.")	
#	logid=E2init(sys.argv)
	
	img = EMData(args[0],0)
	
	if options.bin > 1:
		img = img.process('math.meanshrink',{'n':options.bin})
		
		if options.coords:
			print "WARNING! You will extract particles from a tomogram that is binned/shrunk by %d, and not its raw version" % (options.bin)
			
	#if options.binx > 1 or options.biny > 1 or options.binz > 1:
	#	img = EMData(args[0],0)

	#	options.binx > 1:

	#	options.biny > 1:

	#	options.binz > 1:

	if options.coords:
		commandline_tomoboxer(args[0],options.coords,options.subset,options.boxsize,options.cbin,options.output,options.output_format,options.swapyz,options.reverse_contrast)
	
	else:	
		if options.reverse_contrast:
			img = img*(-1)
						
		app = EMApp()
		if options.inmemory: 
			print "Reading tomogram. Please wait."
			#img=EMData(args[0],0)
			print "Done !"
			boxer=EMTomoBoxer(app,data=img,yshort=options.yshort,boxsize=options.boxsize,bin=options.bin)
		else : 
	#		boxer=EMTomoBoxer(app,datafile=args[0],yshort=options.yshort,apix=options.apix,boxsize=options.boxsize)
			boxer=EMTomoBoxer(app,datafile=img,yshort=options.yshort,apix=options.apix,boxsize=options.boxsize,bin=options.bin)

		boxer.show()
		app.execute()
	#	E2end(logid)

class EMAverageViewer(QtGui.QWidget):
	"""This is a multi-paned view showing a single boxed out particle from a larger tomogram"""
	
	def __init__(self,parent):
		QtGui.QWidget.__init__(self)
		
		self.setWindowTitle("Particle Average")
		
		self.parent=weakref.ref(parent)
		
		self.resize(300,500)
		
		self.gbl = QtGui.QGridLayout(self)
		#self.xyview = EMImage2DWidget()
		#self.gbl.addWidget(self.xyview,0,1)

		#self.xzview = EMImage2DWidget()
		#self.gbl.addWidget(self.xzview,1,1)

		#self.zyview = EMImage2DWidget()
		#self.gbl.addWidget(self.zyview,0,0)

		self.d3view = EMImage3DWidget()
		self.gbl.addWidget(self.d3view,0,0)
		
		self.gbl2 = QtGui.QGridLayout()
		self.gbl.addLayout(self.gbl2,1,0)
		
		self.wfilt = ValSlider(rng=(0,50),label="Filter:",value=0.0)
		self.gbl2.addWidget(self.wfilt,2,0,1,2)

		self.wmask = ValSlider(rng=(0,100),label="Mask:",value=0.0)
		self.gbl2.addWidget(self.wmask,3,0,1,2)

		self.wsymlbl=QtGui.QLabel("Symmetry:")
		self.gbl2.addWidget(self.wsymlbl,4,0)
		
		self.wsym=QtGui.QLineEdit("C1")
		self.gbl2.addWidget(self.wsym,4,1)
		
		self.wprog=QtGui.QProgressBar()
		self.wprog.setRange(0,100)
		self.gbl2.addWidget(self.wprog,5,0,1,2)
		
		self.wrestart=QtGui.QPushButton("Restart")
		self.gbl2.addWidget(self.wrestart,6,1)
		
		self.needupd=0					# Set by the second thread when a display update is ready, 1 means progress update, 2 means volume update
		self.threadrestart=False		# Set by the GUI thread when the second thread needs to restart from scratch
		self.threadprog=0				# Thread progress (0-100)
		self.threadprogstr=""			# String describing thread action
		self.data=None
		
		# These are values from the widgets, stored so the thread can get at them without making GUI calls
		self.sym="c1"
		self.filt=0.0
		self.mask=0.0
		
		QtCore.QObject.connect(self.wfilt,QtCore.SIGNAL("valueChanged")  ,self.event_filter  )
		QtCore.QObject.connect(self.wmask,QtCore.SIGNAL("valueChanged")  ,self.event_mask  )
		QtCore.QObject.connect(self.wsym,QtCore.SIGNAL("editingFinished()")  ,self.event_symchange  )
		QtCore.QObject.connect(self.wrestart,QtCore.SIGNAL("clicked(bool)")  ,self.event_restart  )
		

		# The timer event handles displaying the results processed by the other thread
		self.timer=QtCore.QTimer(self)
		QtCore.QObject.connect(self.timer,QtCore.SIGNAL("timeout")  ,self.event_timer  )
		self.timer.start(500)

		# The processing is all done in the background by the other thread
		self.bgthread=threading.Thread(target=self.thread_process)
		self.bgthread.daemon=True
		self.bgthread.start()

	def event_timer(self):
		if self.needupd&1 :
			self.wprog.setValue(self.threadprog)
		if self.needupd&2 : self.update()
		self.needupd=0

	def event_symchange(self):
		print "sym"
		self.sym=self.wsym.text()
		self.wrestart.setEnabled(True)

	def event_filter(self,value):
		print "filt"
		self.filt=value
		self.wrestart.setEnabled(True)
		
	def event_mask(self,value):
		print "mask"
		self.mask=value
		self.wrestart.setEnabled(True)
		
	def event_restart(self):
		print "restart"
		self.threadrestart=True
		self.wrestart.setEnabled(False)
		
	#def set_data(self,data):
		#"""Sets the current volume to display"""
		
		#self.data=data
		
		#self.update()
		#self.show()
		
	def update(self):
		#if self.wfilt.getValue()!=0.0 :
			#self.fdata=self.data.process("filter.lowpass.gauss",{"cutoff_freq":1.0/self.wfilt.getValue()})
		
		#xyd=self.fdata.process("misc.directional_sum",{"axis":"z"})
		#xzd=self.fdata.process("misc.directional_sum",{"axis":"y"})
		#zyd=self.fdata.process("misc.directional_sum",{"axis":"x"})
		
		#self.xyview.set_data(xyd)	
		#self.xzview.set_data(xzd)
		#self.zyview.set_data(zyd)
		self.d3view.set_data(self.data)
		
	def thread_process(self):
		
		while 1:
			time.sleep(5)
			

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

		self.d3view = EMImage3DWidget()
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
		
	def update(self):
		if self.data==None:
			self.xyview.set_data(None)
			self.xzview.set_data(None)
			self.zyview.set_data(None)
			self.d3view.set_data(None)
			return
		
		if self.wfilt.getValue()!=0.0 :
			self.fdata=self.data.process("filter.lowpass.gauss",{"cutoff_freq":1.0/self.wfilt.getValue()})
		
		xyd=self.fdata.process("misc.directional_sum",{"axis":"z"})
		xzd=self.fdata.process("misc.directional_sum",{"axis":"y"})
		zyd=self.fdata.process("misc.directional_sum",{"axis":"x"})
		
		self.xyview.set_data(xyd)	
		self.xzview.set_data(xzd)
		self.zyview.set_data(zyd)
		self.d3view.set_data(self.fdata)
		
	def event_filter(self,value):
		self.update()

class EMTomoBoxer(QtGui.QMainWindow):
	"""This class represents the EMTomoBoxer application instance.  """
	
	def __init__(self,application,data=None,datafile=None,yshort=False,apix=0.0,boxsize=32,bin=1): #jesus
		
	
		QtGui.QWidget.__init__(self)
		
		self.app=weakref.ref(application)
		self.yshort=yshort
		self.apix=apix
		
		self.bin=bin			#jesus
		
		self.setWindowTitle("MAIN e2tomoboxer.py")
		
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
		self.wboxsize=ValBox(label="Box Size:",value=boxsize)
		self.gbl2.addWidget(self.wboxsize,0,0,1,2)
		self.oldboxsize=boxsize
		
		# Position boxes showing current position		#Jesus - Coordinates should be displayed
		self.positionx=ValBox(label="X:",value=0)
		self.gbl2.addWidget(self.positionx,1,0)

		self.positiony=ValBox(label="Y:",value=0)
		self.gbl2.addWidget(self.positiony,1,1)

		self.positionz=ValBox(label="Z:",value=0)
		self.gbl2.addWidget(self.positionz,1,2)

		
		# max or mean
		self.wmaxmean=QtGui.QPushButton("MaxProj")
		self.wmaxmean.setCheckable(True)
		#self.gbl2.addWidget(self.wmaxmean,1,0)
		self.gbl2.addWidget(self.wmaxmean,4,0)

		
		# number slices
		self.wnlayers=QtGui.QSpinBox()
		self.wnlayers.setMinimum(1)
		self.wnlayers.setMaximum(256)
		self.wnlayers.setValue(1)
		self.gbl2.addWidget(self.wnlayers,4,1)
		
		# Local boxes in side view
		#self.wlocalbox=QtGui.QCheckBox("Limit Side Boxes")		#jesus - should be true by default; no use for checkbox
		#self.gbl2.addWidget(self.wlocalbox,2,0)
		#self.gbl2.addWidget(self.wlocalbox,5,0)

		
		# scale factor
		self.wscale=ValSlider(rng=(.1,2),label="Sca:",value=1.0)
		#self.gbl2.addWidget(self.wscale,3,0,1,2)
		self.gbl2.addWidget(self.wscale,6,0,1,2)

		
		# 2-D filters
		self.wfilt = ValSlider(rng=(0,50),label="Filt:",value=0.0)
		#self.gbl2.addWidget(self.wfilt,4,0,1,2)
		self.gbl2.addWidget(self.wfilt,7,0,1,2)
		
		
		
		self.curbox=-1
		self.boxes=[]						# array of box info, each is (x,y,z,...)
		self.boxesimgs=[]					# z projection of each box
		self.xydown=None
		
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
		QtCore.QObject.connect(self.wboxsize,QtCore.SIGNAL("valueChanged"),self.event_boxsize)
		QtCore.QObject.connect(self.wmaxmean,QtCore.SIGNAL("clicked(bool)"),self.event_projmode)
		QtCore.QObject.connect(self.wscale,QtCore.SIGNAL("valueChanged")  ,self.event_scale  )
		QtCore.QObject.connect(self.wfilt,QtCore.SIGNAL("valueChanged")  ,self.event_filter  )
		
		#QtCore.QObject.connect(self.wlocalbox,QtCore.SIGNAL("stateChanged(int)")  ,self.event_localbox  )	#jesus
		
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

		if datafile!=None : self.set_datafile(datafile)		# This triggers a lot of things to happen, so we do it last
		if data!=None : self.set_data(data)
		
		# Boxviewer subwidget (details of a single box)
		self.boxviewer=EMBoxViewer()
		#self.app().attach_child(self.boxviewer)
		
		# Boxes Viewer (z projections of all boxes)
		self.boxesviewer=EMImageMXWidget()
		#self.app().attach_child(self.boxesviewer)
		self.boxesviewer.show()
		self.boxesviewer.set_mouse_mode("App")
		self.boxesviewer.setWindowTitle("Particle List")

		# Average viewer shows results of background tomographic processing
#		self.averageviewer=EMAverageViewer(self)
		#self.averageviewer.show()
		
		QtCore.QObject.connect(self.boxesviewer,QtCore.SIGNAL("mx_image_selected"),self.img_selected)
		

	def menu_win_boxes(self) : self.boxesviewer.show()
	def menu_win_single(self) : self.boxviewer.show()
#	def menu_win_average(self) : self.averageviewer.show()


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
		if self.yshort : self.datasize=(imgh["nx"],imgh["nz"],imgh["ny"])
		else: self.datasize=(imgh["nx"],imgh["ny"],imgh["nz"])
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
		if self.yshort: self.datasize=(data["nx"],data["nz"],data["ny"])
		else: self.datasize=(data["nx"],data["ny"],data["nz"])
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
				r.process_inplace("xform",{"transform":Transform({"type":"eman","alt":90.0})})
				r.process_inplace("xform.flip",{"axis":"z"})
			elif self.datafile!=None:
				r=EMData(self.datafile,0,0,Region(x-bs/2,z-bs/2,y-bs/2,bs,bs,bs))
				r.process_inplace("xform",{"transform":Transform({"type":"eman","alt":90.0})})
				r.process_inplace("xform.flip",{"axis":"z"})
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
			
		r.process_inplace("normalize.edgemean")
		return r
	

	def get_slice(self,n,xyz):
		"""Reads a slice either from a file or the preloaded memory array. xyz is the axis along which 'n' runs, 0=x (yz), 1=y (xz), 2=z (xy)"""
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
			else: return None
			
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

			else : return None

		if self.apix!=0 : 
			r["apix_x"]=self.apix
			r["apix_y"]=self.apix
			r["apix_z"]=self.apix
		return r

	def event_boxsize(self):
		if self.boxsize()==self.oldboxsize : return
		self.oldboxsize=self.boxsize()
		
		cb=self.curbox
		for i in range(len(self.boxes)) : self.update_box(i)
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

	def boxsize(self): return int(self.wboxsize.getValue())
	
	def nlayers(self): return int(self.wnlayers.value())
	
	def depth(self): return int(self.wdepth.value())

	def scale(self) : return self.wscale.getValue()

	def menu_file_open(self,tog):
		QtGui.QMessageBox.warning(None,"Error","Sorry, in the current version, you must provide a file to open on the command-line.")

	def menu_file_read_boxloc(self):
		fsp=str(QtGui.QFileDialog.getOpenFileName(self, "Select output text file"))
		
		f=file(fsp,"r")
		for b in f:
			b2=[int(float(i)) for i in b.split()[:3]]
			self.boxes.append(b2)
			self.update_box(len(self.boxes)-1)

	def menu_file_save_boxloc(self):
		binf=self.bin 								#jesus

		fsp=str(QtGui.QFileDialog.getSaveFileName(self, "Select output text file"))
		
		out=file(fsp,"w")
		for b in self.boxes:
			out.write("%d\t%d\t%d\n"%(b[0]*binf,b[1]*binf,b[2]*binf))
		out.close()
		
	def menu_file_save_boxes(self):
		fsp=str(QtGui.QFileDialog.getSaveFileName(self, "Select output file (numbers added)"))
		
		progress = QtGui.QProgressDialog("Saving", "Abort", 0, len(self.boxes),None)
		for i,b in enumerate(self.boxes):
			#img=self.get_cube(b[0],b[1],b[2])
			bs=self.boxsize()
			binf=self.bin
			
			if self.yshort:							#jesus
				img = unbinned_extractor(bs,b[0],b[2],b[1],binf) 	#jesus
			else:
				img = unbinned_extractor(bs,b[0],b[1],b[2],binf) 	#jesus
			
			img['origin_x'] = 0						#jesus
			img['origin_y'] = 0				
			img['origin_z'] = 0

			if fsp[:4].lower()=="bdb:" : 
				img.write_image("%s_%03d"%(fsp,i),0)
			elif "." in fsp: 
				img.write_image("%s_%03d.%s"%(fsp.rsplit(".",1)[0],i,fsp.rsplit(".",1)[1]))
			else :
				QtGui.QMessageBox.warning(None,"Error","Please provide a valid image file extension. The numerical sequence will be inserted before the extension.")
				return
			
			progress.setValue(i+1)
			if progress.wasCanceled() : break
		
	def menu_file_save_boxes_stack(self):
		fsp=str(QtGui.QFileDialog.getSaveFileName(self, "Select output file (bdb and hdf only)"))
		
		if fsp[:4].lower()!="bdb:" and fsp[-4:].lower()!=".hdf" :
			QtGui.QMessageBox.warning(None,"Error","3-D stacks supported only for bdb: and .hdf files")
			return
		
		progress = QtGui.QProgressDialog("Saving", "Abort", 0, len(self.boxes),None)
		for i,b in enumerate(self.boxes):
			#img=self.get_cube(b[0],b[1],b[2])
			bs=self.boxsize()
			binf=self.bin

			if self.yshort:							#jesus
				img = unbinned_extractor(bs,b[0],b[2],b[1],binf) 	#jesus
			else:
				img = unbinned_extractor(bs,b[0],b[1],b[2],binf) 	#jesus
			
			img['origin_x'] = 0						#jesus
			img['origin_y'] = 0				
			img['origin_z'] = 0
			
			img.write_image(fsp,i)
			
			progress.setValue(i+1)
			if progress.wasCanceled() : break
		
	def menu_file_quit(self):
		self.close()
		
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
		bs=self.boxsize()
		
		# update shape display
		#if self.wlocalbox.isChecked():						#jesus- limiting side boxes should be default
		
		xzs=self.xzview.get_shapes()
		for i in range(len(self.boxes)):
			if self.boxes[i][1]<self.cury+bs/2 and self.boxes[i][1]>self.cury-bs/2 : 
				xzs[i][0]="rect"
			else : xzs[i][0]="hidden"

		zys=self.zyview.get_shapes()
		for i in range(len(self.boxes)):
			if self.boxes[i][0]<self.curx+bs/2 and self.boxes[i][0]>self.curx-bs/2 : 
				zys[i][0]="rect"
			else : zys[i][0]="hidden"
		
		#else :
		#	xzs=self.xzview.get_shapes()					#jesus
		#	for i in range(len(self.boxes)): xzs[i][0]="rect"
		#	zys=self.zyview.get_shapes()
		#	for i in range(len(self.boxes)): zys[i][0]="rect"
			
		self.xzview.shapechange=1
		self.zyview.shapechange=1
		
		# yz
		avgr=self.get_averager()

		for x in range(x-self.nlayers()/2,x+(self.nlayers()+1)/2):
			slc=self.get_slice(x,0)
			avgr.add_image(slc)
			
		av=avgr.finish()
		if not self.yshort: av.process_inplace("xform.transpose")
		
		if self.wfilt.getValue()!=0.0 : av.process_inplace("filter.lowpass.gauss",{"cutoff_freq":1.0/self.wfilt.getValue()})
		self.zyview.set_data(av)
		
		# xz
		avgr=self.get_averager()

		for y in range(y-self.nlayers()/2,y+(self.nlayers()+1)/2):
			slc=self.get_slice(y,1)
			avgr.add_image(slc)
			
		av=avgr.finish()
		if self.wfilt.getValue()!=0.0 : av.process_inplace("filter.lowpass.gauss",{"cutoff_freq":1.0/self.wfilt.getValue()})
		self.xzview.set_data(av)


	def update_xy(self):
		"""updates xy view due to a new slice range"""
		
		if self.datafile==None and self.data==None: return
		
		
		# Boxes should also be limited by default in the XY view
		if len(self.boxes) > 0:
			zc=self.wdepth.value()
			#print "The current depth is", self.wdepth.value()
			bs=self.boxsize()
			xys=self.xyview.get_shapes()
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
		
		if self.wmaxmean.isChecked() : avgr=Averagers.get("minmax",{"max":1})
		else : avgr=Averagers.get("mean")

		slc=EMData()
		for z in range(self.wdepth.value()-self.nlayers()/2,self.wdepth.value()+(self.nlayers()+1)/2):
			slc=self.get_slice(z,2)
			avgr.add_image(slc)
		
		av=avgr.finish()
		if self.wfilt.getValue()!=0.0 : av.process_inplace("filter.lowpass.gauss",{"cutoff_freq":1.0/self.wfilt.getValue()})
		self.xyview.set_data(av)

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

	def del_box(self,n):
		"""Delete an existing box by replacing the deleted box with the last box. A bit funny, but otherwise
		update after deletion is REALLY slow."""
#		print "del ",n
		if n<0 or n>=len(self.boxes): return
		self.boxviewer.set_data(None)
		self.curbox=-1
		if n==len(self.boxes)-1 : 
			self.boxes.pop()
			self.boxesimgs.pop()
			self.boxesviewer.set_data(self.boxesimgs)
			self.boxesviewer.update()
			self.xyview.del_shape(n)
			self.xzview.del_shape(n)
			self.zyview.del_shape(n)
			self.curbox=-1
			self.update()
		else :
			a=self.boxes.pop()
			self.boxes[n]=a
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

	def update_box(self,n,quiet=False):
		"""After adjusting a box, call this"""
#		print "upd ",n,quiet
		
		box=self.boxes[n]
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
	

		if self.depth()!=box[2] : self.wdepth.setValue(box[2])
		else : self.xyview.update()
		self.update_sides()

		# For speed, we turn off updates while dragging a box around. Quiet is set until the mouse-up
		if not quiet:
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
			
		if n!=self.curbox:
			self.boxesviewer.set_selected((n,),True)

		self.curbox=n


	def img_selected(self,event,lc):
#		print "sel",lc[0]
		if event.modifiers()&Qt.ShiftModifier: self.del_box(lc[0])
		else : self.update_box(lc[0])
		if self.curbox>=0 :
			box=self.boxes[self.curbox]
			self.xyview.scroll_to(box[0],box[1])
			self.xzview.scroll_to(None,box[2])
			self.zyview.scroll_to(box[2],None)

	def xy_down(self,event):
		x,y=self.xyview.scr_to_img((event.x(),event.y()))
		x,y=int(x),int(y)
		self.xydown=None
		if x<0 or y<0 : return		# no clicking outside the image (on 2 sides)
		
		for i in range(len(self.boxes)):
			if self.inside_box(i,x,y) :
				if event.modifiers()&Qt.ShiftModifier: self.del_box(i)
				else:
				
					self.xydown=(i,x,y,self.boxes[i][0],self.boxes[i][1])
					self.update_box(i)
				break
		else:
#			if x>self.boxsize()/2 and x<self.datasize[0]-self.boxsize()/2 and y>self.boxsize()/2 and y<self.datasize[1]-self.boxsize()/2 and self.depth()>self.boxsize()/2 and self.depth()<self.datasize[2]-self.boxsize()/2 :
			if not event.modifiers()&Qt.ShiftModifier:
				self.boxes.append(([x,y,self.depth()]))
				self.xydown=(len(self.boxes)-1,x,y,x,y)		# box #, x down, y down, x box at down, y box at down
				self.update_box(self.xydown[0])

		if self.curbox>=0 :
			box=self.boxes[self.curbox]
			self.xzview.scroll_to(None,box[2])
			self.zyview.scroll_to(box[2],None)
		
	def xy_drag(self,event):
		if self.xydown==None : return
		
		x,y=self.xyview.scr_to_img((event.x(),event.y()))
		x,y=int(x),int(y)
		
		dx=x-self.xydown[1]
		dy=y-self.xydown[2]
		self.boxes[self.xydown[0]][0]=dx+self.xydown[3]
		self.boxes[self.xydown[0]][1]=dy+self.xydown[4]
		self.update_box(self.curbox,True)
		
	def xy_up  (self,event):
		if self.xydown!=None: self.update_box(self.curbox)
		self.xydown=None
	
	def xy_wheel (self,event):
		if event.delta() > 0:
			self.wdepth.setValue(self.wdepth.value()+2)
		elif event.delta() < 0:
			self.wdepth.setValue(self.wdepth.value()-2)
		
		if len(self.boxes) > 0:			#jesus - boxes should disappear in Z when they're in a Z 1/2 the boxsize away from their center
			self.update_xy()		#xy view needs to be updated to take this into account as you scroll
		
	
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
			#if (not self.wlocalbox.isChecked() and self.inside_box(i,x,-1,z)) or self.inside_box(i,x,self.cury,z) :	#jesus ????
			if self.inside_box(i,x,self.cury,z) :

				if event.modifiers()&Qt.ShiftModifier: 
					self.del_box(i)
				else :
					self.xzdown=(i,x,z,self.boxes[i][0],self.boxes[i][2])
					self.update_box(i)
				break
		else:
			if not event.modifiers()&Qt.ShiftModifier:
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
			#if (not self.wlocalbox.isChecked() and self.inside_box(i,-1,y,z)) or  self.inside_box(i,self.curx,y,z):	#jesus - ????
			if self.inside_box(i,self.curx,y,z):
				if event.modifiers()&Qt.ShiftModifier: self.del_box(i)
				else :
					self.zydown=(i,z,y,self.boxes[i][2],self.boxes[i][1])
					self.update_box(i)
					break
		else:
			if not event.modifiers()&Qt.ShiftModifier:
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
		self.boxes[self.zydown[0]][2]=dz+self.zydown[3]
		self.boxes[self.zydown[0]][1]=dy+self.zydown[4]
		self.update_box(self.curbox,True)
		
	def zy_up  (self,event):
		if self.zydown!=None : self.update_box(self.curbox)
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


if __name__ == "__main__":
	main()
		
		
		
