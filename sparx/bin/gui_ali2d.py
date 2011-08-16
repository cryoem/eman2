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


#Layout of the Pop Up window Popuptwodali (for sxali2d); started by the function twodali of the main window        
class Popuptwodali(QWidget):
    def __init__(self):
        QWidget.__init__(self)
	self.picklename='ali2dsaveparms.pkl'
	self.setadv=False
	self.cmd = ""
	self.y=90
	self.x1 = 10
	self.x2 = self.x1 + 150
	self.x3 = self.x2+145
	self.x4 = self.x3+100
	
        #Here we just set the window title
	self.setWindowTitle('sxali2d')
        #Here we just set a label and its position in the window
	title1=QtGui.QLabel('<b>sxali2d</b> - performs 2D reference free alignment of an image series', self)
	title1.move(self.x1,10)


	self.repopbtn = QPushButton("Repopulate With Saved Parameters", self)
        self.repopbtn.move(self.x1-5,40)
        #sets an infotip for this Pushbutton
        self.repopbtn.setToolTip('Repopulate With Saved Parameters')
        #when this button is clicked, this action starts the subfunction twodali
        self.connect(self.repopbtn, SIGNAL("clicked()"), self.repoparms_ali2d)

        #Here we create a Button(file_button with title run open .hdf) and its position in the window
	self.file_button = QtGui.QPushButton("Open .hdf", self)
	self.file_button.move(self.x3, self.y-2)
        #Here we define, that when this button is clicked, it starts subfunction choose_file
	QtCore.QObject.connect(self.file_button, QtCore.SIGNAL("clicked()"), self.choose_file)
        #exactly the same as above, but for subfunction choose_file1
	self.file_button1 = QtGui.QPushButton("Open .bdb", self)
	self.file_button1.move(self.x4,self.y-2)
	QtCore.QObject.connect(self.file_button1, QtCore.SIGNAL("clicked()"), self.choose_file1)
	
	
        #not linked to a function yet
        self.activeheader_button = QtGui.QPushButton("activate all images", self)
	self.activeheader_button.move(self.x1-5, 340)
	self.connect(self.activeheader_button, SIGNAL("clicked()"), self.setactiveheader)
	
	self.a2dheader_button = QtGui.QPushButton("set xform.align2d", self)
	self.a2dheader_button.move(self.x1-5+180, 340)
	self.connect(self.a2dheader_button, SIGNAL("clicked()"), self.seta2dheader)
	
	self.advbtn = QPushButton("Advanced Parameters", self)
        self.advbtn.move(self.x1-5, 390)
        #sets an infotip for this Pushbutton
        self.advbtn.setToolTip('Set Advanced Parameters for ali2d such as center and CTF')
        #when this button is clicked, this action starts the subfunction twodali
        self.connect(self.advbtn, SIGNAL("clicked()"), self.advparams)
	
	self.savepbtn = QPushButton("Save Input Parameters", self)
        self.savepbtn.move(self.x1-5, 420)
        #sets an infotip for this Pushbutton
        self.savepbtn.setToolTip('Save Input Parameters')
        #when this button is clicked, this action starts the subfunction twodali
        self.connect(self.savepbtn, SIGNAL("clicked()"), self.saveparms)
	
	self.cmdlinebtn = QPushButton("Generate command line from input parameters", self)
        self.cmdlinebtn.move(self.x1-5, 450)
        #sets an infotip for this Pushbutton
        self.cmdlinebtn.setToolTip('Generate command line using input parameters')
        #when this button is clicked, this action starts the subfunction twodali
        self.connect(self.cmdlinebtn, SIGNAL("clicked()"), self.gencmdline)

	 #Here we create a Button(Run_button with title run sxali2d) and its position in the window
	self.RUN_button = QtGui.QPushButton('Run sxali2d', self)
	# make 3D textured push button look
	s = "QPushButton {font: bold; color: #000;border: 1px solid #333;border-radius: 11px;padding: 2px;background: qradialgradient(cx: 0, cy: 0,fx: 0.5, fy:0.5,radius: 1, stop: 0 #fff, stop: 1 #8D0);min-width:90px;margin:5px} QPushButton:pressed {font: bold; color: #000;border: 1px solid #333;border-radius: 11px;padding: 2px;background: qradialgradient(cx: 0, cy: 0,fx: 0.5, fy:0.5,radius: 1, stop: 0 #fff, stop: 1 #084);min-width:90px;margin:5px}"
	
	self.RUN_button.setStyleSheet(s)
	
	
	self.RUN_button.move(230, 500)
        #Here we define, that when this button is clicked, it starts subfunction runsxali2d
        self.connect(self.RUN_button, SIGNAL("clicked()"), self.runsxali2d)
        #Labels and Line Edits for User Input

	# populate with default values
	self.savedparmsdict = {'stackname':'NONE','foldername':'NONE','partradius':'-1','xyrange':'4 2 1 1','trans':'2 1 0.5 0.25','nriter':'3','nproc':'1','maskname':'','center':'-1',"ringstep":"1","innerradius":"1","ctf":"False","snr":"1.0","fourvar":"False", "gpnr":"-1","usrfunc":"ref_ali2d","usrfuncfile":""}
	
        #Example for User input stack name
        #First create the label and define its position
	
	stackname= QtGui.QLabel('Name of input stack', self)
	stackname.move(self.x1,self.y)
        #Now add a line edit and define its position
	self.stacknameedit=QtGui.QLineEdit(self)
        self.stacknameedit.move(self.x2,self.y)
        #Adds a default value for the line edit
	self.stacknameedit.setText(self.savedparmsdict['stackname'])
	#The same as above, but many line edits include Infotips
	foldername= QtGui.QLabel('Output folder', self)
	foldername.move(self.x1,self.y+30)
	self.foldernameedit=QtGui.QLineEdit(self)
	self.foldernameedit.move(self.x2,self.y+30)
	self.foldernameedit.setText(self.savedparmsdict['foldername'])	
	
	self.outinfobtn = QPushButton("Output Info", self)
        self.outinfobtn.move(self.x3,  self.y+30)
        #sets an infotip for this Pushbutton
        self.outinfobtn.setToolTip('Output Info')
        #when this button is clicked, this action starts the subfunction twodali
        self.connect(self.outinfobtn, SIGNAL("clicked()"), self.outputinfo_ali2d)
	
	partradius= QtGui.QLabel('Particle radius', self)
	partradius.move(self.x1,self.y+60)
	self.partradiusedit=QtGui.QLineEdit(self)
	self.partradiusedit.move(self.x2,self.y+60)
	self.partradiusedit.setText(self.savedparmsdict['partradius'])
	self.partradiusedit.setToolTip('Parameter ou: Outer radius for rotational correlation \nshould be set to particle radius\nif not sure, set to boxsize/2-2 ')	
	xyrange= QtGui.QLabel('xy range', self)
	xyrange.move(self.x1,self.y+90)
	self.xyrangeedit=QtGui.QLineEdit(self)
	self.xyrangeedit.move(self.x2,self.y+90)
	self.xyrangeedit.setText(self.savedparmsdict['xyrange'])
	self.xyrangeedit.setToolTip('Range for translational search in x, y direction\nif set to 0 only rotational alignment will be performed')
	trans= QtGui.QLabel('translational step', self)
	trans.move(self.x1,self.y+120)
	self.transedit=QtGui.QLineEdit(self)
	self.transedit.move(self.x2,self.y+120)
	self.transedit.setText(self.savedparmsdict['trans'])
	self.transedit.setToolTip('Step of translational search in x, y direction\nlarger values increase the speed but decrease the accuracy')	
	nriter= QtGui.QLabel('Number of Iterations', self)
	nriter.move(self.x1,self.y+150)
	self.nriteredit=QtGui.QLineEdit(self)
	self.nriteredit.move(self.x2,self.y+150)
	self.nriteredit.setText(self.savedparmsdict['nriter'])
	self.nriteredit.setToolTip('Maximum number of iterations the program will perform\n Using the default values the program will run 3 rounds with xy-range 4 and translational step 1, 3 rounds with xyrange 2 and translational step 1 and so on..\nif set to 0 maximum iteration number will be 10 and will automatically stop should the criterion falls')
	nproc= QtGui.QLabel('Number of Processors', self)
	nproc.move(self.x1,self.y+180)
	self.nprocedit=QtGui.QLineEdit(self)
	self.nprocedit.move(self.x2,self.y+180)
	self.nprocedit.setText(self.savedparmsdict['nproc'])
	self.nprocedit.setToolTip('The number of processors to use. Default is single processor mode')

	header=QtGui.QLabel('Attributes xform.align2d and active parameters must be set in the input stack', self)
	header.move(self.x1,310)
	
     #Function runsxali2d started when  the  RUN_button of the  Poptwodali window is clicked 
    def outputinfo_ali2d(self):
    	QMessageBox.information(self, "ali2d output",'Output files (average of aligned images and Fourier Ring Correlation curve)\nare saved in Output folder. alignment parameters are saved in the attribute \nxform.align2d in each image header. The images themselves are not changed.')
	
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
	
	cmd1 = "sxali2d.py "+str(stack) +" "+ str(output)
	
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
		(fname,stat)= QInputDialog.getText(self,"Generate Command Line","Enter name of file to save command line in",QLineEdit.Normal,"")
		if stat:
			f = open(fname,'a')
			f.write(cmd1)
			f.write('\n')
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
	output=open(self.picklename,'wb')
	self.gencmdline(writefile=False)
	pickle.dump(self.savedparmsdict,output)
	output.close()
	
    def repoparms_ali2d(self):	
	# repopulate with saved parms
	import pickle
	pkl = open(self.picklename,'rb')
	self.savedparmsdict = pickle.load(pkl)
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

    def seta2dheader(self):
	#opens a file browser, showing files only in .hdf format
	ok=False
	zerostr='set xform.align2d to zero'
	randstr='randomize xform.align2d'
	importstr='import parameters from file'
	(item,stat)= QInputDialog.getItem(self,"xform.align2d","choose option",[zerostr,randstr,importstr])
        #we convert this Qstring to a string and send it to line edit classed stackname edit of the Poptwodali window
	#self.stacknameedit.setText(str(a))
   	choice= str(item)
	stack = self.stacknameedit.text()
	if stat:
		if choice == zerostr:
			header(str(stack),'xform.align2d',zero=True)
		if choice == randstr:
			header(str(stack),'xform.align2d',rand_alpha=True)
		if choice == importstr:
			file_name = QtGui.QFileDialog.getOpenFileName(self, "Open Data File", "", "(*)")
			a=str(QtCore.QString(file_name))
			header(str(stack),'xform.align2d',fimport=a)
	
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
	
	self.mskfile_button = QtGui.QPushButton("Open .hdf", self)
	self.mskfile_button.move(285, 60-2)
        #Here we define, that when this button is clicked, it starts subfunction choose_file
	QtCore.QObject.connect(self.mskfile_button, QtCore.SIGNAL("clicked()"), self.choose_mskfile)
	
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
		
	usrfuncfile= QtGui.QLabel('Enter name of external file containing user function:', self)
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
	
    def choose_mskfile(self):
	#opens a file browser, showing files only in .hdf format
   	file_name = QtGui.QFileDialog.getOpenFileName(self, "Open File Containing Mask", "", "(*)")
        #after the user selected a file, we obtain this filename as a Qstring
	a=QtCore.QString(file_name)
	print a
        #we convert this Qstring to a string and send it to line edit classed stackname edit of the Poptwodali window
	self.masknameedit.setText(str(a))
		
