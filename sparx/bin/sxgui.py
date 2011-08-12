#!/usr/bin/env python

#
# Author: Pawel A.Penczek, 09/09/2006 (Pawel.A.Penczek@uth.tmc.edu)
# Copyright (c) 2000-2006 The University of Texas - Houston Medical School
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
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
#
#

import global_def
from global_def import *
import sys
from PyQt4.Qt import *
from PyQt4 import QtGui
from PyQt4 import QtCore
import os
import subprocess
from subprocess import *
from EMAN2 import *
from sparx import *
from EMAN2_cppwrap import *

#Layout of the Pop Up window infosparx; started by the function info of the main window
class infosparx(QWidget):
    def __init__(self):
        QWidget.__init__(self)
        #Here we just set the window title and  3 different labels, with their positions in the window
	self.setWindowTitle('Sparx GUI Info page')
        title1=QtGui.QLabel('<b>Sparx GUI beta version</b>', self)
	title1.move(10,10)
        title2=QtGui.QLabel('<b>Authors:</b> ', self)
	title2.move(10,40)
        title3=QtGui.QLabel('For more information visit:\nhttp://sparx-em.org/sparxwiki ', self)
	title3.move(10,70)

#Layout of the Pop Up window Popuptwodali (for sxali2d); started by the function twodali of the main window        
class Popuptwodali(QWidget):
    def __init__(self):
        QWidget.__init__(self)
	
	self.setadv=False
	self.cmd = ""
	self.y=90
	self.x = 10
        #Here we just set the window title
	self.setWindowTitle('sxali2d')
        #Here we just set a label and its position in the window
	title1=QtGui.QLabel('<b>sxali2d</b> - performs 2D reference free alignment of an image series', self)
	title1.move(10,10)


	self.repopbtn = QPushButton("Repopulate With Previously Saved Parameters", self)
        self.repopbtn.move(self.x-5,40)
        #sets an infotip for this Pushbutton
        self.repopbtn.setToolTip('Repopulate With Saved Parameters')
        #when this button is clicked, this action starts the subfunction twodali
        self.connect(self.repopbtn, SIGNAL("clicked()"), self.repoparms)

        #Here we create a Button(file_button with title run open .hdf) and its position in the window
	self.file_button = QtGui.QPushButton("Open .hdf", self)
	self.file_button.move(285, self.y-2)
        #Here we define, that when this button is clicked, it starts subfunction choose_file
	QtCore.QObject.connect(self.file_button, QtCore.SIGNAL("clicked()"), self.choose_file)
        #exactly the same as above, but for subfunction choose_file1
	self.file_button1 = QtGui.QPushButton("Open .bdb", self)
	self.file_button1.move(385,self.y-2)
	QtCore.QObject.connect(self.file_button1, QtCore.SIGNAL("clicked()"), self.choose_file1)
	
	
        #not linked to a function yet
        self.activeheader_button = QtGui.QPushButton("activate all images", self)
	self.activeheader_button.move(self.x-5, 340)
	self.connect(self.activeheader_button, SIGNAL("clicked()"), self.setactiveheader)
	self.ali2dheader_button = QtGui.QPushButton("set xform.align2d to zero", self)
	self.ali2dheader_button.move(self.x-5, 370)
	self.connect(self.ali2dheader_button, SIGNAL("clicked()"), self.setali2dheader)
	ort=QtGui.QLabel('<b>or</b>', self)
	ort.move(220,370)
	self.ali2drndheader_button = QtGui.QPushButton("set xform.align2d to random values", self)
	self.ali2drndheader_button.move(self.x-5+240, 370)
	self.connect(self.ali2drndheader_button, SIGNAL("clicked()"), self.setali2drndheader)
	
	self.advbtn = QPushButton("Advanced Parameters", self)
        self.advbtn.move(self.x-5, 420)
        #sets an infotip for this Pushbutton
        self.advbtn.setToolTip('Set Advanced Parameters for ali2d such as center and CTF')
        #when this button is clicked, this action starts the subfunction twodali
        self.connect(self.advbtn, SIGNAL("clicked()"), self.advparams)
	
	self.savepbtn = QPushButton("Save Input Parameters", self)
        self.savepbtn.move(self.x-5, 450)
        #sets an infotip for this Pushbutton
        self.savepbtn.setToolTip('Save Input Parameters')
        #when this button is clicked, this action starts the subfunction twodali
        self.connect(self.savepbtn, SIGNAL("clicked()"), self.saveparms)
	
	self.cmdlinebtn = QPushButton("Generate command line from input parameters", self)
        self.cmdlinebtn.move(self.x-5, 480)
        #sets an infotip for this Pushbutton
        self.cmdlinebtn.setToolTip('Generate command line using input parameters')
        #when this button is clicked, this action starts the subfunction twodali
        self.connect(self.cmdlinebtn, SIGNAL("clicked()"), self.gencmdline)

	 #Here we create a Button(Run_button with title run sxali2d) and its position in the window
	self.RUN_button = QtGui.QPushButton('Run sxali2d', self)
	# make 3D textured push button look
	s = "QPushButton {font: bold; color: #000;border: 1px solid #333;border-radius: 11px;padding: 2px;background: qradialgradient(cx: 0, cy: 0,fx: 0.5, fy:0.5,radius: 1, stop: 0 #fff, stop: 1 #8D0);min-width:90px;margin:5px} QPushButton:pressed {font: bold; color: #000;border: 1px solid #333;border-radius: 11px;padding: 2px;background: qradialgradient(cx: 0, cy: 0,fx: 0.5, fy:0.5,radius: 1, stop: 0 #fff, stop: 1 #084);min-width:90px;margin:5px}"
	
	self.RUN_button.setStyleSheet(s)
	
	
	self.RUN_button.move(230, 520)
        #Here we define, that when this button is clicked, it starts subfunction runsxali2d
        self.connect(self.RUN_button, SIGNAL("clicked()"), self.runsxali2d)
        #Labels and Line Edits for User Input

	# populate with default values
	self.savedparmsdict = {'stackname':'NONE','foldername':'NONE','partradius':'-1','xyrange':'4 2 1 1','trans':'2 1 0.5 0.25','nriter':'3','nproc':'1','maskname':'','center':'-1',"ringstep":"1","innerradius":"1","ctf":"False","snr":"1.0","fourvar":"False", "gpnr":"-1","usrfunc":"ref_ali2d","usrfuncfile":""}
	
        #Example for User input stack name
        #First create the label and define its position
	
	stackname= QtGui.QLabel('Name of your stack', self)
	stackname.move(10,self.y)
        #Now add a line edit and define its position
	self.stacknameedit=QtGui.QLineEdit(self)
        self.stacknameedit.move(140,self.y)
        #Adds a default value for the line edit
	self.stacknameedit.setText(self.savedparmsdict['stackname'])
	#The same as above, but many line edits include Infotips
	foldername= QtGui.QLabel('Output folder', self)
	foldername.move(10,self.y+30)
	self.foldernameedit=QtGui.QLineEdit(self)
	self.foldernameedit.move(140,self.y+30)
	self.foldernameedit.setText(self.savedparmsdict['foldername'])	
	partradius= QtGui.QLabel('Particle radius', self)
	partradius.move(10,self.y+60)
	self.partradiusedit=QtGui.QLineEdit(self)
	self.partradiusedit.move(140,self.y+60)
	self.partradiusedit.setText(self.savedparmsdict['partradius'])
	self.partradiusedit.setToolTip('Parameter ou: Outer radius for rotational correlation \nshould be set to particle radius\nif not sure, set to boxsize/2-2 ')	
	xyrange= QtGui.QLabel('xy range', self)
	xyrange.move(10,self.y+90)
	self.xyrangeedit=QtGui.QLineEdit(self)
	self.xyrangeedit.move(140,self.y+90)
	self.xyrangeedit.setText(self.savedparmsdict['xyrange'])
	self.xyrangeedit.setToolTip('Range for translational search in x, y direction\nif set to 0 only rotational alignment will be performed')
	trans= QtGui.QLabel('translational step', self)
	trans.move(10,self.y+120)
	self.transedit=QtGui.QLineEdit(self)
	self.transedit.move(140,self.y+120)
	self.transedit.setText(self.savedparmsdict['trans'])
	self.transedit.setToolTip('Step of translational search in x, y direction\nlarger values increase the speed but decrease the accuracy')	
	nriter= QtGui.QLabel('Nr. of Iterations', self)
	nriter.move(10,self.y+150)
	self.nriteredit=QtGui.QLineEdit(self)
	self.nriteredit.move(140,self.y+150)
	self.nriteredit.setText(self.savedparmsdict['nriter'])
	self.nriteredit.setToolTip('Maximum number of iterations the program will perform\n Using the default values the program will run 3 rounds with xy-range 4 and translational step 1, 3 rounds with xyrange 2 and translational step 1 and so on..\nif set to 0 maximum iteration number will be 10 and will automatically stop should the criterion falls')
	nproc= QtGui.QLabel('Nr. of Processors', self)
	nproc.move(10,self.y+180)
	self.nprocedit=QtGui.QLineEdit(self)
	self.nprocedit.move(140,self.y+180)
	self.nprocedit.setText(self.savedparmsdict['nproc'])
	self.nprocedit.setToolTip('The number of processors to use. Default is single processor mode')

	header=QtGui.QLabel('Attributes xform.align2d and active parameters must be set in the input stack', self)
	header.move(10,310)
	
     #Function runsxali2d started when  the  RUN_button of the  Poptwodali window is clicked 
   
    def gencmdline(self,writefile=True):
	#Here we just read in all user inputs in the line edits of the Poptwodali window
   	stack = self.stacknameedit.text()
	print "stack defined="+ stack
	output=self.foldernameedit.text()
	print "output folder="+ output
	ou=self.partradiusedit.text()
	print "Particle radius="+ ou
	xr=self.xyrangeedit.text()
	print "xr=" +xr
	yr=self.xyrangeedit.text()
	print "yr=" +yr
	ts=self.transedit.text()
	print "ts=" +ts
	maxit=self.nriteredit.text()
	print "maxit="+maxit
	
	cmd1 = " sxali2d.py "+str(stack) +" "+ str(output)
	
	args = " --ou="+ str(ou)+ " --xr='"+str(xr)+"'"+ " --yr='"+str(yr)+"'"+ " --ts='"+str(ts)+"'"+  " --maxit="+ str(maxit) 
	
	mask=''
	ctr=''
	ringstep=''
	inrad=''
	CTF=''
	snr=''
	fourvar=''
	gpn=''
	userf=''
	userfile=''
	
	if self.setadv:
		mask = self.w.masknameedit.text()
		if len(str(mask))> 1:
			cmd1 = cmd1+" "+str(mask) 
	cmd1 = cmd1 + args
	
	if self.setadv:	
		ctr=self.w.centeredit.text()
		cmd1 = cmd1+" --center=" +str(ctr)
		
		ringstep = self.w.ringstepedit.text()
		cmd1 = cmd1+" --rs="+str(ringstep)
		
		inrad = self.w.innerradiusedit.text()
		cmd1 = cmd1 + " --ir=" + str(inrad)
		
		CTF=str(self.w.ctfedit.text())
		if str(CTF) == 'True':
			cmd1 = cmd1 + " --CTF"
			
		snr = self.w.snredit.text()
		cmd1 = cmd1 + " --snr=" + str(snr)
		
		fourvar = self.w.fourvaredit.text()
		if str(fourvar) == 'True':
			cmd1 = cmd1 + " --Fourvar"
			
		gpn = self.w.gpnredit.text()
		cmd1 = cmd1 + " --Ng=" + str(gpn)
		
		userf = self.w.usrfuncedit.text()
		
		userfile = self.w.usrfuncfileedit.text()
		
		if len(userfile) <= 1:
			cmd1 = cmd1 + " --function="+str(userf)
		else:
			userfile = self.w.usrfuncfileedit.text()
			userfile = str(userfile)
			if len(userfile) > 1:
				# break it down into file name and directory path
				rind = userfile.rfind('/')
				fname = userfile[rind+1:]
				fname, ext = os.path.splitext(fname)
				fdir = userfile[0:rind]
				cmd1 = cmd1 + " --function=\"[" +fdir+","+fname+","+str(userf)+"]\""
	np = self.nprocedit.text()
	
	self.savedparmsdict = {'stackname':str(stack),'foldername':str(output),'partradius':str(ou),'xyrange':str(xr),'trans':str(yr),'nriter':str(maxit),'nproc':str(np),'maskname':str(mask),'center':str(ctr),"ringstep":str(ringstep),"innerradius":str(inrad),"ctf":str(CTF),"snr":str(snr),"fourvar":str(fourvar), "gpnr":str(gpn),"usrfunc":str(userf), "usrfuncfile":str(userfile)}
	
	if self.setadv:
		self.w.savedparmsdict=self.savedparmsdict
	
	if int(str(np)) > 1:
		cmd1="mpirun -np "+ str(np) + " "+ cmd1+" --MPI" 
	
	if writefile:	
		from time import localtime
		a=time.localtime()
		fname = 'ali2dcmd_%02d_%02d_%04d_%02d_%02d_%02d.txt'%(a.tm_mday,a.tm_mon, a.tm_year,a.tm_hour, a.tm_min, a.tm_sec)
		f = open(fname,'a')
		f.write(cmd1)
		f.close()
	
	print cmd1
	self.cmd = cmd1
	
    def runsxali2d(self):
	self.gencmdline(writefile=False)
	outfolder=self.savedparmsdict['foldername']
	if os.path.exists(outfolder):
		print "output folder "+outfolder+" already exists!"
		return
	process = subprocess.Popen(self.cmd,shell=True)
	self.emit(QtCore.SIGNAL("process_started"),process.pid)
	
    def saveparms(self):	
	# save all the parms in a text file so we can repopulate if user requests
	import pickle
	output=open('savedparms.pkl','wb')
	self.gencmdline(writefile=False)
	pickle.dump(self.savedparmsdict,output)
	output.close()
	
    def repoparms(self):	
	# repopulate with saved parms
	import pickle
	pkl = open('savedparms.pkl','rb')
	self.savedparmsdict = pickle.load(pkl)
	print self.savedparmsdict
	self.partradiusedit.setText(self.savedparmsdict['partradius'])
	self.stacknameedit.setText(self.savedparmsdict['stackname'])	
	self.foldernameedit.setText(self.savedparmsdict['foldername'])	
	self.xyrangeedit.setText(self.savedparmsdict['xyrange'])
	self.transedit.setText(self.savedparmsdict['trans'])
	self.nriteredit.setText(self.savedparmsdict['nriter'])
	self.nprocedit.setText(self.savedparmsdict['nproc'])
	if self.setadv:
		self.w.masknameedit.setText(self.savedparmsdict['maskname'])
		self.w.centeredit.setText(self.savedparmsdict['center'])
		self.w.ringstepedit.setText(self.savedparmsdict['ringstep'])
		self.w.innerradiusedit.setText(self.savedparmsdict['innerradius'])
		self.w.ctfedit.setText(self.savedparmsdict['ctf'])
		self.w.snredit.setText(self.savedparmsdict['snr'])
		self.w.fourvaredit.setText(self.savedparmsdict['fourvar'])
		self.w.gpnredit.setText(self.savedparmsdict['gpnr'])
		self.w.usrfuncedit.setText(self.savedparmsdict['usrfunc'])
		self.w.usrfuncfileedit.setText(self.savedparmsdict['usrfuncfile'])
		
    def advparams(self):
        print "Opening a new popup window..."
        self.w = Popupadvparams_ali2d(self.savedparmsdict)
        self.w.resize(500,450)
        self.w.show()
	self.setadv=True
    def setactiveheader(self):
	stack = self.stacknameedit.text()
	print "stack defined="+ stack
	header(str(stack), "active", one=True)
	
    def setali2dheader(self):
	#Here we just read in all user inputs in the line edits of the Poptwodali window
	stack = self.stacknameedit.text()
	print "stack defined="+ stack
	header(str(stack),'xform.align2d',zero=True)

    def setali2drndheader(self):
	#Here we just read in all user inputs in the line edits of the Poptwodali window
	stack = self.stacknameedit.text()
	print "stack defined="+ stack
	header(str(stack),'xform.align2d',rand_alpha=True)
	
	#Function choose_file started when  the  open_file of the  Poptwodali window is clicked
    def choose_file(self):
	#opens a file browser, showing files only in .hdf format
   	file_name = QtGui.QFileDialog.getOpenFileName(self, "Open Data File", "", "HDF files (*.hdf)")
        #after the user selected a file, we obtain this filename as a Qstring
	a=QtCore.QString(file_name)
	print a
        #we convert this Qstring to a string and send it to line edit classed stackname edit of the Poptwodali window
	self.stacknameedit.setText(str(a))
        
	#Function choose_file started when  the  open_file of the  Poptwodali window is clicked (same as above but for bdb files(maybe we can combine these two into one function)
    def choose_file1(self):
	file_name1 = QtGui.QFileDialog.getOpenFileName(self, "Open Data File", "EMAN2DB/", "BDB FILES (*.bdb)" )
	a=QtCore.QString(file_name1)
	b=os.path.basename(str(a))
	c=os.path.splitext(b)[0]
	d="bdb:"+c
	print d
	self.stacknameedit.setText(d)

class Popupadvparams_ali2d(QWidget):
    def __init__(self,savedparms):
        QWidget.__init__(self)
        #Here we just set the window title
	self.setWindowTitle('sxali2d advanced parameter selection')
        #Here we just set a label and its position in the window
	title1=QtGui.QLabel('<b>sxali2d</b> - set advanced params', self)
	title1.move(10,10)
        #Labels and Line Edits for User Input
        #Just a label
	title2= QtGui.QLabel('<b>Advanced</b> parameters', self)
	title2.move(10,40)

	self.savedparmsdict=savedparms
        #Example for User input stack name
        #First create the label and define its position
	maskname= QtGui.QLabel('Mask', self)
	maskname.move(10,60)
        #Now add a line edit and define its position
	self.masknameedit=QtGui.QLineEdit(self)
        self.masknameedit.move(140,60)
        #Adds a default value for the line edit
	self.masknameedit.setText(self.savedparmsdict['maskname'])
	self.masknameedit.setToolTip("Default is a circle mask with radius equal to the particle radius")
	
	center= QtGui.QLabel('Center type', self)
	center.move(10,90)
	self.centeredit=QtGui.QLineEdit(self)
	self.centeredit.move(140,90)
	self.centeredit.setText(self.savedparmsdict['center'])
	self.centeredit.setToolTip('-1 - use average centering method (default),\n0 - if you do not want the average to be centered, \n1 - phase approximation of the center of gravity phase_cog, \n2 - cross-correlate with Gaussian function, \n3 - cross-correlate with donut shape image (e.g. inner radius=2, outer radius=7), \n4 - cross-correlate with reference image provided by user, \n5 - cross-correlate with self-rotated average..\ncentering may fail..use 0 to deactive it')
	
	ringstep= QtGui.QLabel('Ring step', self)
	ringstep.move(10,120)
	self.ringstepedit=QtGui.QLineEdit(self)
	self.ringstepedit.move(140,120)
	self.ringstepedit.setText(self.savedparmsdict['ringstep'])
	self.ringstepedit.setToolTip('step between rings in rotational correlation > 0 (set to 1)')

	innerradius= QtGui.QLabel('Inner radius', self)
	innerradius.move(10,150)
	self.innerradiusedit=QtGui.QLineEdit(self)
	self.innerradiusedit.move(140,150)
	self.innerradiusedit.setText(self.savedparmsdict['innerradius'])
	self.innerradiusedit.setToolTip('inner radius for rotational correlation > 0 (set to 1) ')	
	
	ctf= QtGui.QLabel('CTF', self)
	ctf.move(10,180)
	self.ctfedit=QtGui.QLineEdit(self)
	self.ctfedit.move(140,180)
	self.ctfedit.setText(self.savedparmsdict['ctf'])
	self.ctfedit.setToolTip('if this flag is set, the program will use CTF information provided in file headers')

	snr= QtGui.QLabel('SNR', self)
	snr.move(10,210)
	self.snredit=QtGui.QLineEdit(self)
	self.snredit.move(140,210)
	self.snredit.setText(self.savedparmsdict['snr'])
	self.snredit.setToolTip('signal-to-noise ratio of the data (default SNR=1.0)')		
	fourvar= QtGui.QLabel('Fourvar', self)
	fourvar.move(10,240)
	self.fourvaredit=QtGui.QLineEdit(self)
	self.fourvaredit.move(140,240)
	self.fourvaredit.setText(self.savedparmsdict['fourvar'])
	self.fourvaredit.setToolTip('use Fourier variance to weight the reference (recommended, default False)')

	gpnr= QtGui.QLabel('Number of Groups', self)
	gpnr.move(10,270)
	self.gpnredit=QtGui.QLineEdit(self)
	self.gpnredit.move(140,270)
	self.gpnredit.setText(self.savedparmsdict['gpnr'])
	self.gpnredit.setToolTip('number of groups in the new CTF filteration')	
	
	usrfunc= QtGui.QLabel('User Function Name', self)
	usrfunc.move(10,300)
	self.usrfuncedit=QtGui.QLineEdit(self)
	self.usrfuncedit.move(140,300)
	self.usrfuncedit.setText(self.savedparmsdict['usrfunc'])
	self.usrfuncedit.setToolTip('name of the user-supplied-function that prepares reference image for each iteration')
		
	usrfuncfile= QtGui.QLabel('Enter (full path) name of external file containing user function:', self)
	usrfuncfile.move(10,330)
	usrfuncfile= QtGui.QLabel('(Leave blank if file is not external to Sparx)', self)
	usrfuncfile.move(10,350)
	self.usrfuncfileedit=QtGui.QLineEdit(self)
	self.usrfuncfileedit.move(140,370)
	self.usrfuncfileedit.setText(self.savedparmsdict['usrfuncfile'])
	self.usrfuncfileedit.setToolTip('name of the external file containing user function')	
     #Function runsxali2d started when  the  RUN_button of the  Poptwodali window is clicked 
    	
        #Here we define, that when this button is clicked, it starts subfunction choose_file
	self.usrfile_button = QtGui.QPushButton("Select File", self)
	self.usrfile_button.move(285, 370-2)
	QtCore.QObject.connect(self.usrfile_button, QtCore.SIGNAL("clicked()"), self.choose_usrfile)
	
    def choose_usrfile(self):
	#opens a file browser, showing files only in .hdf format
   	file_name = QtGui.QFileDialog.getOpenFileName(self, "Open File Containing User Fuction", "", "py files (*.py)")
        #after the user selected a file, we obtain this filename as a Qstring
	a=QtCore.QString(file_name)
	print a
        #we convert this Qstring to a string and send it to line edit classed stackname edit of the Poptwodali window
	self.usrfuncfileedit.setText(str(a))
	
#Layout of the Pop Up window Popuptwodali (for sxali2d); started by the function twodali of the main window 
class PopupHelicalRefinement(QWidget):
    def __init__(self):
        
        QtGui.QWidget.__init__(self)
        self.setWindowTitle(self.tr("Helical Refinment"))
        self.__create_ui()
    def __create_ui(self):
	

        input_stack_label = QtGui.QLabel(self.tr("Input Stack:"))
	self.input_stack_line_edit = QtGui.QLineEdit()
	self.input_stack_line_edit.setMinimumWidth(200)
	self.input_stack_browse_button = QtGui.QPushButton(self.tr("Browse"))
	input_stack_layout = QtGui.QHBoxLayout()
	input_stack_layout.addWidget(input_stack_label)
	input_stack_layout.addWidget(self.input_stack_line_edit)
	input_stack_layout.addWidget(self.input_stack_browse_button)
	
	input_parms_label = QtGui.QLabel(self.tr("Input Projection Parameters:"))
	self.input_param_line_edit = QtGui.QLineEdit()
	self.input_param_line_edit.setMinimumWidth(200)
	self.input_param_browse_button = QtGui.QPushButton(self.tr("Browse"))
	input_param_layout = QtGui.QHBoxLayout()
	input_param_layout.addWidget(input_parms_label)
	input_param_layout.addWidget(self.input_param_line_edit)
	input_param_layout.addWidget(self.input_param_browse_button)
	

	self.param_ative_button = QtGui.QPushButton("Open .bdb", self)
	self.param_zeros_button = QtGui.QPushButton("Open .bdb", self)
	param_setting_layout = QtGui.QHBoxLayout()
	param_setting_layout.addWidget(self.param_ative_button)
	param_setting_layout.addWidget(self.param_zeros_button)
	

	self.vbl = QtGui.QVBoxLayout(self)
	self.vbl.setMargin(0)
	self.vbl.setSpacing(6)
	self.vbl.setObjectName("vbl")
	self.vbl.addLayout(input_stack_layout)
	self.vbl.addLayout(input_param_layout)
	self.vbl.addLayout(param_setting_layout)
	   


   
###MAIN WINDOW	(started by class App)
#This class includes the layout of the main window; within each class, i name the main object self, to avoid confusion)    	
class MainWindow(QtGui.QWidget):
    def __init__(self, parent=None):
        QtGui.QWidget.__init__(self, parent)
        #sets the title of the window
	self.setWindowTitle('SPARX GUI')
        #creates a Pushbutton, named sxali2d, defines its position in the window 
        self.btn1 = QPushButton("sxali2d", self)
        self.btn1.move(10, 65)
        #sets an infotip for this Pushbutton
        self.btn1.setToolTip('2D  reference free alignment of an image series ')
        #when this button is clicked, this action starts the subfunction twodali
        self.connect(self.btn1, SIGNAL("clicked()"), self.twodali)
        #another Pushbutton, with a tooltip, not linked to a function yet
	self.btn2 = QPushButton("sxmref_ali2d", self)
        self.btn2.setToolTip('2D MULTI-referencealignment of an image series ')
        self.btn2.move(10, 95)
        #Pushbutton named Info
	self.picbutton = QPushButton(self)
        #when this button is clicked, this action starts the subfunction info
        self.connect(self.picbutton, SIGNAL("clicked()"), self.info)
	#creates a Pushbutton, named sxihrsr defines its position in the window 
        self.btn3 = QPushButton("sxihrsr", self)
        self.btn3.move(10, 125)
        #sets an infotip for this Pushbutton
        self.btn3.setToolTip('Iterative Real Space Helical Refinement ')
        #when this button is clicked, this action starts the subfunction twodali
        self.connect(self.btn3, SIGNAL("clicked()"), self.helicalrefinement)
        #this decorates the button with the sparx image
	icon = QIcon(get_image_directory()+"sparxicon.png")
	self.picbutton.setIcon(icon)
        self.picbutton.move(120, 5)
        self.picbutton.setToolTip('Info Page')
        #Quitbutton
        self.btn3 = QPushButton("Close", self)
        self.btn3.setToolTip('Close SPARX GUI ')
	self.btn3.move(180, 5)
        self.connect(self.btn3, QtCore.SIGNAL('clicked()'),QtGui.qApp, QtCore.SLOT('quit()'))
        #here we set two labels, their position and font style
	title=QtGui.QLabel('<b>SPARX</b> GUI', self)
	title1=QtGui.QLabel('<b>2D alignment</b> applications', self)
	title.move(10,10)
	title1.move(20, 45)
	QtGui.QToolTip.setFont(QtGui.QFont('OldEnglish', 8))
	
    	#here we set the width and height of the main window
        self.resize(300,350)
	
   
    #This is the function two2ali, which is being started when the Pushbutton btn1 of the main window(called sxali2d) is being clicked
    def twodali(self):
        print "Opening a new popup window..."
        #opens the window Poptwodali, and defines its width and height
        #The layout of the Poptwodali window is defined in class Poptwodali(QWidget Window)
        self.w = Popuptwodali()
        self.w.resize(600,600)
        self.w.show()
    def helicalrefinement(self):
        print "Opening a new popup window..."
        #opens the window Poptwodali, and defines its width and height
        #The layout of the Poptwodali window is defined in class Poptwodali(QWidget Window)
        self.w = PopupHelicalRefinement()
        self.w.resize(600,600)
        self.w.show()    
        
    #This is the function info, which is being started when the Pushbutton picbutton of the main window is being clicked
    def info(self):
        print "Opening a new popup window..."
        #opens the window infosparx, and defines its width and height
        #The layout of the infosparx window is defined in class infosparx(QWidget Window)
        self.w = infosparx()
        self.w.resize(250,200)
        self.w.show()

#  this is the main class of the program
#  Here we provide the necessary imports. The basic GUI widgets are located in QtGui module.
class App(QApplication):
    def __init__(self, *args):
        QApplication.__init__(self, *args)
        #here we define the main window (class MainWindow)
        self.main = MainWindow()
        #here we define that when all windows are closed, function byebye of class App will be started
        self.connect(self, SIGNAL("lastWindowClosed()"), self.byebye )
        #hshows main window
        self.main.show()
        
    #function byebye (just quit)  
    def byebye( self ):
        print' bye bye!'
        self.exit(0)

      
#  Necessary for execution of the program
def main(args):
    global app
    app = App(args)
    app.exec_()

if __name__ == "__main__":
    main(sys.argv)
