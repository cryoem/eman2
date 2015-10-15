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
import  sys
from    PyQt4.Qt import *
from    PyQt4 import QtGui
from    PyQt4 import QtCore
import  os
import  subprocess
from    subprocess import *
from    EMAN2 import *
from    sparx import *
from    EMAN2_cppwrap import *

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

#helical_start
class PopupHelicalRefinement(QWidget):
	def __init__(self):
		QWidget.__init__(self)
		
		
		# populate with default values
		self.savedparmsdict ={'stackname':'NONE','initialprojectionparameter':'NONE','referencevolume':'NONE','foldername':'NONE','outradius':'-1','xrange':'1.0','xtrans':'1.0','ynumber':"2",'nriter':'3','nproc':'3','dp':'NONE','dphi':'NONE','rmax':'NONE','maskname':'',"delta":"1.0","ringstep":"1","innerradius":"1","ctf":Qt.Unchecked,"snr":"1.0","initial_theta":"90.0", "delta_theta":"1.0","nise":"2","nise":"2","rmin":"0.0","fract":"0.7","dp_step":"0.1","ndp":"12","dphi_step":"0.1","ndphi":"12","psi_max":"15","an":"-1", "npad":"2", "chunk":"-1.0", "sym":"c1","datasym":"symdoc.dat","usrfunc":"helical","usrfuncfile":"",'stackname_prectr':'','outdir_prectr':'','mask_prectr':'','search_rng_prectr':'-1','ou_prectr':'-1','maxit_prectr':'100','snr_prectr':'1','fourvar_prectr':Qt.Unchecked,'ctf_prectr':Qt.Unchecked,'oneDx_prectr':Qt.Unchecked,'nproc_prectr':'2'}
		
		
		self.setadv=False
		self.cmd = ""
		x1 = 10
		x2 = x1 + 170
		x3 = x2 + 145
		x4 = x3 + 100
		x5 = 230 # run button
		#Here we just set the window title
		self.setWindowTitle('sxhelical')
		#Here we just set a label and its position in the window
		title1=QtGui.QLabel('<b>shelical</b> - performs helical refinement', self)
		y = 10
		title1.move(10,y)
		
		
		y = y +30
		self.repopbtn = QPushButton("Retrieve saved parameters", self)
		self.repopbtn.move(x1-5,y)
		#sets an infotip for this Pushbutton
		self.repopbtn.setToolTip('Retrieve saved parameters')
		#when this button is clicked, this action starts the subfunction twodali
		self.connect(self.repopbtn, SIGNAL("clicked()"), self.repoparms_helical)

		
		#Here we create a Button(file_button with title run open .hdf) and its position in the window
		
		y = y +30		
		self.file_button = QtGui.QPushButton("Open .hdf", self)
		self.file_button.move(x3, y)
		#Here we define, that when this button is clicked, it starts subfunction choose_file
		QtCore.QObject.connect(self.file_button, QtCore.SIGNAL("clicked()"), self.choose_file)
		#exactly the same as above, but for subfunction choose_file1
		self.file_button1 = QtGui.QPushButton("Open .bdb", self)
		self.file_button1.move(x4,y)
		QtCore.QObject.connect(self.file_button1, QtCore.SIGNAL("clicked()"), self.choose_file1)
		#Example for User input stack name
		#First create the label and define its position
		stackname= QtGui.QLabel('Name of input stack', self)
		stackname.move(x1,y)
		#Now add a line edit and define its position
		self.stacknameedit=QtGui.QLineEdit(self)
		self.stacknameedit.move(x2,y)
		#Adds a default value for the line edit
		self.stacknameedit.setText(self.savedparmsdict['stackname'])
		
		# file_button2 for input parameters
		
		y = y +30		
		self.file_button2 = QtGui.QPushButton("Open .txt", self)
		self.file_button2.move(x3, y)
		#Here we define, that when this button is clicked, it starts subfunction choose_file
		QtCore.QObject.connect(self.file_button2, QtCore.SIGNAL("clicked()"), self.choose_file2)
		parameterfile= QtGui.QLabel('Projection params file', self)
		parameterfile.move(x1,y)
		#Now add a line edit and define its position
		self.initialprojectionparameteredit=QtGui.QLineEdit(self)
		self.initialprojectionparameteredit.move(x2,y)
		#Adds a default value for the line edit
		self.initialprojectionparameteredit.setText(self.savedparmsdict['initialprojectionparameter'])
		
		# file_button3 for referencevolume
		y = y +30		
		self.file_button3 = QtGui.QPushButton("Open .hdf", self)
		self.file_button3.move(x3, y)
		#Here we define, that when this button is clicked, it starts subfunction choose_file
		QtCore.QObject.connect(self.file_button3, QtCore.SIGNAL("clicked()"), self.choose_file3)
		volumefile= QtGui.QLabel('Initial volume', self)
		volumefile.move(x1,y)
		#Now add a line edit and define its position
		self.referencevolumeedit=QtGui.QLineEdit(self)
		self.referencevolumeedit.move(x2,y)
		#Adds a default value for the line edit
		self.referencevolumeedit.setText(self.savedparmsdict['referencevolume'])
		

		#The same as above, but many line edits include Infotips
		y = y +30		
		foldername= QtGui.QLabel('Output folder', self)
		foldername.move(x1,y)
		self.foldernameedit=QtGui.QLineEdit(self)
		self.foldernameedit.move(x2,y)
		self.foldernameedit.setText(self.savedparmsdict['foldername'])		
		
		self.outinfobtn = QPushButton("Output Info", self)
		self.outinfobtn.move(x3,  y)
		#sets an infotip for this Pushbutton
		self.outinfobtn.setToolTip('Output Info')
		#when this button is clicked, this action starts the subfunction twodali
		self.connect(self.outinfobtn, SIGNAL("clicked()"), self.outputinfo_helical)
		
		y = y +30
		outradius= QtGui.QLabel('Particle outer radius', self)
		outradius.move(x1, y )
		self.outradiusedit=QtGui.QLineEdit(self)
		self.outradiusedit.move(x2,y)
		self.outradiusedit.setText(self.savedparmsdict['outradius'])
		self.outradiusedit.setToolTip('Parameter ou: Outer radius for rotational correlation \nshould be set to particle radius\nif not sure, set to boxsize/2-2 ')
		
		y = y +30		
		xrange= QtGui.QLabel('X search range', self)
		xrange.move(x1,y)
		self.xrangeedit=QtGui.QLineEdit(self)
		self.xrangeedit.move(x2,y)
		self.xrangeedit.setText(self.savedparmsdict['xrange'])
		self.xrangeedit.setToolTip('Range for translational search in x\nif set to 0 no x directional alignment will be performed')
		y = y +30
		xtrans= QtGui.QLabel('Step of x search in pixel', self)
		xtrans.move(x1,y)
		self.xtransedit=QtGui.QLineEdit(self)
		self.xtransedit.move(x2,y)
		self.xtransedit.setText(self.savedparmsdict['xtrans'])
		self.xtransedit.setToolTip('Step of translational search in x direction\nlarger values increase the speed but decrease the accuracy')
		
		y = y +30		
		ynumber= QtGui.QLabel('Number of y search', self)
		ynumber.move(x1,y)
		self.ynumberedit=QtGui.QLineEdit(self)
		self.ynumberedit.move(x2,y)
		self.ynumberedit.setText(self.savedparmsdict['ynumber'])
		self.ynumberedit.setToolTip('number of steps in y direction\n ystep will be dp/ynumber')
				
		y = y +30		
		dp= QtGui.QLabel('Helical rise dz(angstrom)', self)
		dp.move(x1,y)
		self.dpedit=QtGui.QLineEdit(self)
		self.dpedit.move(x2,y)
		self.dpedit.setText(self.savedparmsdict['dp'])
		self.dpedit.setToolTip('helical rise in angstrom')
		
		y = y +30		
		dphi= QtGui.QLabel('Helical angle', self)
		dphi.move(x1,y)
		self.dphiedit=QtGui.QLineEdit(self)
		self.dphiedit.move(x2,y)
		self.dphiedit.setText(self.savedparmsdict['dp'])
		self.dphiedit.setToolTip('helical angle in degree')		
		
		y = y +30		
		rmax= QtGui.QLabel('Outer radius of helix', self)
		rmax.move(x1,y)
		self.rmaxedit=QtGui.QLineEdit(self)
		self.rmaxedit.move(x2,y)
		self.rmaxedit.setText(self.savedparmsdict['rmax'])
		self.rmaxedit.setToolTip('helical out radious')		

		y = y +30
		nriter= QtGui.QLabel('Number of iterations', self)
		nriter.move(x1,y)
		self.nriteredit=QtGui.QLineEdit(self)
		self.nriteredit.move(x2,y)
		self.nriteredit.setText(self.savedparmsdict['nriter'])
		self.nriteredit.setToolTip('Maximum number of iterations the program will perform\n ')
		
		y = y +30
		nproc= QtGui.QLabel('Number of CPUs (>1)', self)
		nproc.move(x1,y)
		self.nprocedit=QtGui.QLineEdit(self)
		self.nprocedit.move(x2,y)
		self.nprocedit.setText(self.savedparmsdict['nproc'])
		self.nprocedit.setToolTip('The number of processors to use, need to be >=2 since we only surrport MPI version')
		
	  
		y = y +30
		header=QtGui.QLabel('Active parameters and xform.projection must be set in the input stack', self)
		header.move(x1,y)
		
		#not linked to a function yet
		
		y = y +30
		self.activeheader_button = QtGui.QPushButton("Activate all images", self)
		self.activeheader_button.move(x1-5, y)
		self.connect(self.activeheader_button, SIGNAL("clicked()"), self.setactiveheader)
		self.projectionheader_button = QtGui.QPushButton("Set xform.projection", self)
		self.projectionheader_button.move(x1+180, y)
		self.projectionheader_button.setToolTip('Depending whether the projection parameters given or not, set xform projection to zero or specific vlues ')
		self.connect(self.projectionheader_button, SIGNAL("clicked()"), self.setprojectionheader)
		
		
		y = y +30
		self.savepbtn = QPushButton("Save Input Parameters", self)
		self.savepbtn.move(x1-5, y)
		#sets an infotip for this Pushbutton
		self.savepbtn.setToolTip('Save Input Parameters')
		#when this button is clicked, this action starts the subfunction twodali
		self.connect(self.savepbtn, SIGNAL("clicked()"), self.saveparms_helical)
		
		y = y +30
		self.cmdlinebtn = QPushButton("Generate command line from input parameters", self)
		self.cmdlinebtn.move(x1-5, y)
		#sets an infotip for this Pushbutton
		self.cmdlinebtn.setToolTip('Generate command line using input parameters')
		#when this button is clicked, this action starts the subfunction twodali
		self.connect(self.cmdlinebtn, SIGNAL("clicked()"), self.gencmdline_helical)

		#Here we create a Button(Run_button with title run sxali2d) and its position in the window
		
		y = y +30
		self.RUN_button = QtGui.QPushButton('Run sxhelical', self)
		# make 3D textured push button look
		s = "QPushButton {font: bold; color: #000;border: 1px solid #333;border-radius: 11px;padding: 2px;background: qradialgradient(cx: 0, cy: 0,fx: 0.5, fy:0.5,radius: 1, stop: 0 #fff, stop: 1 #8D0);min-width:90px;margin:5px} QPushButton:pressed {font: bold; color: #000;border: 1px solid #333;border-radius: 11px;padding: 2px;background: qradialgradient(cx: 0, cy: 0,fx: 0.5, fy:0.5,radius: 1, stop: 0 #fff, stop: 1 #084);min-width:90px;margin:5px}"
		
		self.RUN_button.setStyleSheet(s)


		self.RUN_button.move(x5, y)
		#Here we define, that when this button is clicked, it starts subfunction runsxali2d
		self.connect(self.RUN_button, SIGNAL("clicked()"), self.runsxhelical)
		#Labels and Line Edits for User Input
				
	def outputinfo_helical(self):
		QMessageBox.information(self, "helical output",'Output files (logfile, projection parameter files and pix error files)\nvol--reconstructed volume, volf--reconstructed volume after helicising and user function \n xform.projection was stored in each image header.')
		

		
	#Function runsxali2d started when  the  RUN_button of the  Poptwodali window is clicked 

	def gencmdline_helical(self,writefile=True):
		#Here we just read in all user inputs in the line edits of the Poptwodali window
		stack = self.stacknameedit.text()
		print "stack defined="+ stack
		
		projectionparameters = self.initialprojectionparameteredit.text()
		print "Input parameter files="+ projectionparameters
		referencevolume = self.referencevolumeedit.text()
		print "Initial volume="+ referencevolume
				
		output=self.foldernameedit.text()
		print "output folder="+ output
		ou=self.outradiusedit.text()
		print "Outer radius="+ ou
		xr=self.xrangeedit.text()
		print "xr=" +xr
		ynumber=self.ynumberedit.text()
		print "ynumber=" +ynumber		
		tx=self.xtransedit.text()
		print "tx=" +tx
		dp=self.dpedit.text()
		print "dp=" +dp
		dphi=self.dphiedit.text()
		print "dp=" +dphi
		rmax=self.rmaxedit.text()
		print "rmax==",rmax
		maxit=self.nriteredit.text()
		print "maxit="+maxit
		
		cmd1 = " sxhelical.py "+str(stack) +" "+str(referencevolume)+" " + str(output)
		
		args = " --ou="+ str(ou)+ " --xr="+str(xr)+" " + " --ynumber="+str(ynumber)+ " " + " --txs="+str(tx)+" " +" --rmax="+str(rmax)+" " + " --maxit="+ str(maxit) 
		
		userf=''
		userfile=''
		
		
		mask = self.w1.masknameedit.text()
		if len(str(mask))> 1:
				cmd1 = cmd1+" "+str(mask) 
		cmd1 = cmd1 + args
		
				
		delta=self.w1.deltaedit.text()
		cmd1 = cmd1+" --delta=" +str(delta)
				
		ringstep = self.w1.ringstepedit.text()
		cmd1 = cmd1+" --rs="+str(ringstep)
				
		inrad = self.w1.innerradiusedit.text()
		cmd1 = cmd1 + " --ir=" + str(inrad)
				
								
		snr = self.w1.snredit.text()
		cmd1 = cmd1 + " --snr=" + str(snr)
				
		initial_theta = self.w1.initial_thetaedit.text()
		cmd1 = cmd1 + " --initial_theta=" + str(initial_theta)
						
		delta_theta = self.w1.delta_thetaedit.text()
		cmd1 = cmd1 + " --delta_theta=" + str(delta_theta)
				

				
		an = self.w1.anedit.text()
		cmd1 = cmd1 + " --an=" + str(an)
				
		npad = self.w1.npadedit.text()
		cmd1 = cmd1 + " --npad=" + str(npad)
				
		chunk = self.w1.chunkedit.text()
		cmd1 = cmd1 + " --chunk=" + str(chunk)
		
		nise = self.w2.niseedit.text()
		cmd1 = cmd1 + " --nise=" + str(nise)
				
		rmin = self.w2.rminedit.text()
		cmd1 = cmd1 + " --rmin=" + str(rmin)
				
		fract = self.w2.fractedit.text()
		cmd1 = cmd1 + " --fract=" + str(fract)
		cmd1 = cmd1 +" --dp="+str(dp)		 
		dp_step = self.w2.dp_stepedit.text()
		cmd1 = cmd1 + " --dp_step=" + str(dp_step)
				
		ndp = self.w2.ndpedit.text()
		cmd1 = cmd1 + " --ndp=" + str(ndp)
		
		cmd1 = cmd1 +" --dphi="+str(dphi)		
		dphi_step = self.w2.dphi_stepedit.text()
		cmd1 = cmd1 + " --dphi_step=" + str(dphi_step)
				
		ndphi = self.w2.ndphiedit.text()
		cmd1 = cmd1 + " --ndphi=" + str(ndphi)
				
		psi_max = self.w2.psi_maxedit.text()
		cmd1 = cmd1 + " --psi_max=" + str(psi_max)		
		sym = self.w2.symedit.text()
		cmd1 = cmd1 + " --sym=" + str(sym)
		datasym = self.w2.datasymedit.text()
		cmd1 = cmd1 + " --datasym=" + str(datasym)
		CTF=self.w1.ctfchkbx.checkState()
		if CTF == Qt.Checked:
				cmd1 = cmd1 + " --CTF"		
		userf = self.w1.usrfuncedit.text()
				
		userfile = self.w1.usrfuncfileedit.text()
				
		if len(userfile) <= 1:
				cmd1 = cmd1 + " --function="+str(userf)
		else:
			userfile = self.w1.usrfuncfileedit.text()
			userfile = str(userfile)
			if len(userfile) > 1:
				# break it down into file name and directory path
				rind = userfile.rfind('/')
				fname = userfile[rind+1:]
				fname, ext = os.path.splitext(fname)
				fdir = userfile[0:rind]
				cmd1 = cmd1 + " --function=\"[" +fdir+","+fname+","+str(userf)+"]\""
		
		np = self.nprocedit.text()
		
		(self.savedparmsdict)['stackname']=str(stack)
		(self.savedparmsdict)['initialprojectionparameter']=str(projectionparameters)
		(self.savedparmsdict)['referencevolume']=str(referencevolume)
		(self.savedparmsdict)['foldername']=str(output)
		(self.savedparmsdict)['outradius']=str(ou)
		(self.savedparmsdict)['xrange']=str(xr)
		(self.savedparmsdict)['xtrans']=str(tx)
		(self.savedparmsdict)['ynumber']=str(ynumber)
		(self.savedparmsdict)['dp']=str(dp)
		(self.savedparmsdict)['dphi']=str(dphi)
		(self.savedparmsdict)['rmax']=str(rmax)
		(self.savedparmsdict)['nriter']=str(maxit)
		(self.savedparmsdict)['nproc']=str(np)
		(self.savedparmsdict)['maskname']=str(mask)
		(self.savedparmsdict)['delta']=str(delta)
		(self.savedparmsdict)['ringstep']=str(ringstep)
		(self.savedparmsdict)['innerradius']=str(inrad)
		(self.savedparmsdict)['ctf']=CTF
		(self.savedparmsdict)['snr']=str(snr)
		(self.savedparmsdict)['initial_theta']=str(initial_theta)
		(self.savedparmsdict)['delta_theta']=str(delta_theta)
		(self.savedparmsdict)['nise']=str(nise)
		(self.savedparmsdict)['rmin']=str(rmin)
		(self.savedparmsdict)['fract']=str(fract)
		(self.savedparmsdict)['dp_step']=str(dp_step)
		(self.savedparmsdict)['ndp']=str(ndp)
		(self.savedparmsdict)['dphi_step']=str(dphi_step)
		(self.savedparmsdict)['ndphi']=str(ndphi)
		(self.savedparmsdict)['psi_max']=str(psi_max)
		(self.savedparmsdict)['an']=str(an)
		(self.savedparmsdict)['npad']=str(npad)
		(self.savedparmsdict)['chunk']=str(chunk)
		(self.savedparmsdict)['sym']=str(sym)
		(self.savedparmsdict)['datasym']=str(datasym)
		(self.savedparmsdict)['usrfunc']=str(userf)
		(self.savedparmsdict)['usrfuncfile']=str(userfile)
		
		
		#self.savedparmsdict = {'stackname':str(stack),'initialprojectionparameter':str(projectionparameters),'referencevolume':str(referencevolume),'foldername':str(output),'outradius':str(ou),'xrange':str(xr),'xtrans':str(tx),'ynumber':str(ynumber),'dp':str(dp),'dphi':str(dphi),'rmax':str(rmax),'nriter':str(maxit),'nproc':str(np),'maskname':str(mask),'delta':str(delta),"ringstep":str(ringstep),"innerradius":str(inrad),"ctf":CTF,"snr":str(snr),"initial_theta":str(initial_theta), "delta_theta":str(delta_theta),"nise":str(nise),"rmin":str(rmin),"fract":str(fract),"dp_step":str(dp_step),"ndp":str(ndp),"dphi_step":str(dphi_step),"ndphi":str(ndphi),"psi_max":str(psi_max),"an":str(an), "npad":str(npad), "chunk":str(chunk),"sym":str(sym),"datasym":str(datasym),"usrfunc":str(userf), "usrfuncfile":str(userfile)}
		
		
		self.w1.savedparmsdict=self.savedparmsdict
		self.w2.savedparmsdict=self.savedparmsdict
		self.w3.savedparmsdict=self.savedparmsdict
		
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
		
	def runsxhelical(self):
		self.gencmdline_helical(writefile=False)
		outfolder=self.savedparmsdict['foldername']
		if os.path.exists(outfolder):
			print "output folder "+outfolder+" already exists!"
			return
		process = subprocess.Popen(self.cmd,shell=True)
		self.emit(QtCore.SIGNAL("process_started"),process.pid)
		
	def saveparms_helical(self):		
		# save all the parms in a text file so we can repopulate if user requests

		(fname,stat)= QInputDialog.getText(self,"Save Input Parameters","Enter name of file to save parameters in",QLineEdit.Normal,"")
		if stat:
			import pickle
			output=open(fname,'wb')
			self.gencmdline_helical(writefile=False)
			self.w3.gencmdline_shftali(writefile=False)
			pickle.dump(self.savedparmsdict,output)
			output.close()

	def repoparms_helical(self):
		  
		# repopulate with saved parms
		(fname,stat)= QInputDialog.getText(self,"Retrieve saved parameters","Enter name of file parameters were saved in",QLineEdit.Normal,"")
		if stat:
			import pickle
			pkl = open(fname,'rb')
			#after the user selected a file, we obtain this filename as a Qstring
			self.savedparmsdict = pickle.load(pkl)
			print self.savedparmsdict
			self.outradiusedit.setText(self.savedparmsdict['outradius'])
			self.stacknameedit.setText(self.savedparmsdict['stackname'])
			self.initialprojectionparameteredit.setText(self.savedparmsdict['initialprojectionparameter'])		
			self.referencevolumeedit.setText(self.savedparmsdict['referencevolume'])		
			self.foldernameedit.setText(self.savedparmsdict['foldername'])		
			self.xrangeedit.setText(self.savedparmsdict['xrange'])
			self.xtransedit.setText(self.savedparmsdict['xtrans'])
			self.ynumberedit.setText(self.savedparmsdict['ynumber'])
			self.dpedit.setText(self.savedparmsdict['dp'])
			self.dphiedit.setText(self.savedparmsdict['dphi'])
			self.rmaxedit.setText(self.savedparmsdict['rmax'])
			self.nriteredit.setText(self.savedparmsdict['nriter'])
			self.nprocedit.setText(self.savedparmsdict['nproc'])
			
			self.w1.ctfchkbx.setCheckState(self.savedparmsdict['ctf'])
			self.w1.masknameedit.setText(self.savedparmsdict['maskname'])
			self.w1.deltaedit.setText(self.savedparmsdict['delta'])
			self.w1.ringstepedit.setText(self.savedparmsdict['ringstep'])
			self.w1.innerradiusedit.setText(self.savedparmsdict['innerradius'])
			self.w1.snredit.setText(self.savedparmsdict['snr'])
			self.w1.initial_thetaedit.setText(self.savedparmsdict['initial_theta'])
			self.w1.delta_thetaedit.setText(self.savedparmsdict['delta_theta'])
			self.w1.delta_thetaedit.setText(self.savedparmsdict['delta_theta'])
			self.w1.anedit.setText(self.savedparmsdict['an'])
			self.w1.npadedit.setText(self.savedparmsdict['npad'])
			self.w1.chunkedit.setText(self.savedparmsdict['chunk'])
			self.w1.usrfuncedit.setText(self.savedparmsdict['usrfunc'])
			self.w1.usrfuncfileedit.setText(self.savedparmsdict['usrfuncfile'])
			
			self.w2.niseedit.setText(self.savedparmsdict['nise'])
			self.w2.symedit.setText(self.savedparmsdict['sym'])
			self.w2.datasymedit.setText(self.savedparmsdict['datasym'])
			self.w2.rminedit.setText(self.savedparmsdict['rmin'])
			self.w2.fractedit.setText(self.savedparmsdict['fract'])
			self.w2.dp_stepedit.setText(self.savedparmsdict['dp_step'])
			self.w2.ndpedit.setText(self.savedparmsdict['ndp'])
			self.w2.dphi_stepedit.setText(self.savedparmsdict['dphi_step'])
			self.w2.ndphiedit.setText(self.savedparmsdict['ndphi'])
			self.w2.psi_maxedit.setText(self.savedparmsdict['psi_max'])
			
			
			self.w3.stacknameedit.setText(self.savedparmsdict['stackname_prectr'])
			self.w3.outdiredit.setText(self.savedparmsdict['outdir_prectr'])
			self.w3.maskedit.setText(self.savedparmsdict['mask_prectr'])
			self.w3.search_rngedit.setText(self.savedparmsdict['search_rng_prectr'])
			self.w3.ouedit.setText(self.savedparmsdict['ou_prectr'])
			self.w3.maxitedit.setText(self.savedparmsdict['maxit_prectr'])
			self.w3.snredit.setText(self.savedparmsdict['snr_prectr'])
			self.w3.ctfchkbx.setCheckState(self.savedparmsdict['ctf_prectr'])
			self.w3.fourvarchkbx.setCheckState(self.savedparmsdict['fourvar_prectr'])
			self.w3.oneDxchkbx.setCheckState(self.savedparmsdict['oneDx_prectr'])
			self.w3.nprocedit.setText(self.savedparmsdict['nproc_prectr'])
				

	def setactiveheader(self):
		stack = self.stacknameedit.text()
		print "stack defined="+ stack
		# horatio active_refactoring Jy51i1EwmLD4tWZ9_00000_1
		# header(str(stack), "active", one=True)
		
	def setprojectionheader(self):
		#Here we just read in all user inputs in the line edits of the Poptwodali window
		stack = self.stacknameedit.text()
		input_file = self.initialprojectionparameteredit.text()
		print "stack defined="+ stack
		if( str(input_file) == "NONE" or len( str(input_file) ) <=1  ):
			header(str(stack),'xform.projection',zero=True)
			print "zero projection"
		else:
			header(str(stack),'xform.projection',fimport = str(input_file) )
			print "set projection based on input file"		
		
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
	def choose_file2(self):
		#opens a file browser, showing files only in .hdf format
		file_name = QtGui.QFileDialog.getOpenFileName(self, "Open Data File", "", "HDF files (*.txt)")
		#after the user selected a file, we obtain this filename as a Qstring
		a=QtCore.QString(file_name)
		print a
		#we convert this Qstring to a string and send it to line edit classed stackname edit of the Poptwodali window
		self.initialprojectionparameteredit.setText(str(a))
		
	def choose_file3(self):
		#opens a file browser, showing files only in .hdf format
		file_name = QtGui.QFileDialog.getOpenFileName(self, "Open Reference Volume", "", "hdf files (*.hdf)")
		#after the user selected a file, we obtain this filename as a Qstring
		a=QtCore.QString(file_name)
		print a
		#we convert this Qstring to a string and send it to line edit classed stackname edit of the Poptwodali window
		self.referencevolumeedit.setText(str(a))

class Popupadvparams_helical_1(QWidget):
	def __init__(self,savedparms):
		QWidget.__init__(self)
		#Here we just set the window title
		self.setWindowTitle('sxhelical advanced parameter selection')
		#Here we just set a label and its position in the window
		x1 = 10
		x2 = x1 + 150
		x3 = x2 + 145
		x4 = x3 + 100
		y = 10
		title1=QtGui.QLabel('<b>sxhelical</b> - set advanced params', self)
		title1.move(10,y)
		#Labels and Line Edits for User Input
		#Just a label
		y = y +20
		title2= QtGui.QLabel('<b>Advanced</b> parameters', self)
		title2.move(10,y)

		self.savedparmsdict=savedparms
		#Example for User input stack name
		#First create the label and define its position
		y = y+20
		maskname= QtGui.QLabel('Mask', self)
		maskname.move(x1,y)
		#Now add a line edit and define its position
		self.masknameedit=QtGui.QLineEdit(self)
		self.masknameedit.move(x2,y)
		#Adds a default value for the line edit
		self.masknameedit.setText(self.savedparmsdict['maskname'])
		self.masknameedit.setToolTip("Default is a circle mask with radius equal to the particle radius")
		
		self.mskfile_button = QtGui.QPushButton("Open .hdf", self)
		self.mskfile_button.move(x3, y-2)
		#Here we define, that when this button is clicked, it starts subfunction choose_file
		QtCore.QObject.connect(self.mskfile_button, QtCore.SIGNAL("clicked()"), self.choose_mskfile)
		
		y = y +30
		ctf= QtGui.QLabel('CTF', self)
		ctf.move(x1,y)
		self.ctfchkbx = QtGui.QCheckBox("",self)
		self.ctfchkbx.move(x2, y)
		self.ctfchkbx.setCheckState(self.savedparmsdict['ctf'])
		self.ctfchkbx.setToolTip('Consider CTF correction during the alignment ')
		
		y = y +30
		snr= QtGui.QLabel('SNR', self)
		snr.move(x1,y)
		self.snredit=QtGui.QLineEdit(self)
		self.snredit.move(x2,y)
		self.snredit.setText(self.savedparmsdict['snr'])
		self.snredit.setToolTip('signal-to-noise ratio of the data (default SNR=1.0)')		
		
		y = y+30
		delta= QtGui.QLabel('Angular step', self)
		delta.move(x1,y)
		self.deltaedit=QtGui.QLineEdit(self)
		self.deltaedit.move(x2,y)
		self.deltaedit.setText(self.savedparmsdict["delta"])
		self.deltaedit.setToolTip("angular step of reference projections\n defualt = 1.0")
		
		y = y + 30
		ringstep= QtGui.QLabel('Ring step', self)
		ringstep.move(x1,y)
		self.ringstepedit=QtGui.QLineEdit(self)
		self.ringstepedit.move(x2,y)
		self.ringstepedit.setText(self.savedparmsdict['ringstep'])
		self.ringstepedit.setToolTip('step between rings in rotational correlation > 0 (set to 1)')

		y = y + 30
		innerradius= QtGui.QLabel('Particle inner radius', self)
		innerradius.move(x1,y)
		self.innerradiusedit=QtGui.QLineEdit(self)
		self.innerradiusedit.move(x2,y)
		self.innerradiusedit.setText(self.savedparmsdict['innerradius'])
		self.innerradiusedit.setToolTip('inner radius for rotational correlation > 0 (set to 1) ')		
		
		y = y + 30
		initial_theta= QtGui.QLabel('Initial theta', self)
		initial_theta.move(x1,y)
		self.initial_thetaedit=QtGui.QLineEdit(self)
		self.initial_thetaedit.move(x2,y)
		self.initial_thetaedit.setText(self.savedparmsdict['initial_theta'])
		self.initial_thetaedit.setToolTip('intial theta for reference projection, default 90.0)')
		
		y = y + 30
		delta_theta= QtGui.QLabel('Theta step', self)
		delta_theta.move(x1,y)
		self.delta_thetaedit=QtGui.QLineEdit(self)
		self.delta_thetaedit.move(x2,y)
		self.delta_thetaedit.setText(self.savedparmsdict['delta_theta'])
		self.delta_thetaedit.setToolTip('step of theta for reference projection, default 1.0)')

		y = y + 30
		an = QtGui.QLabel('Angular neighborhood', self)
		an.move(x1,y)
		self.anedit=QtGui.QLineEdit(self)
		self.anedit.move(x2,y)
		self.anedit.setText(self.savedparmsdict['an'])
		self.anedit.setToolTip('angular neighborhood for local searches, defaut -1')		
		
		y = y + 30
		npad = QtGui.QLabel('Npad', self)
		npad.move(x1,y)
		self.npadedit=QtGui.QLineEdit(self)
		self.npadedit.move(x2,y)
		self.npadedit.setText(self.savedparmsdict['npad'])
		self.npadedit.setToolTip('padding size for 3D reconstruction, please use npad = 2')		
		
		y = y + 30
		chunk = QtGui.QLabel('Chunk', self)
		chunk.move(x1,y)
		self.chunkedit=QtGui.QLineEdit(self)
		self.chunkedit.move(x2,y)
		self.chunkedit.setText(self.savedparmsdict['chunk'])
		self.chunkedit.setToolTip('percentage of data used for alignment')		

		y = y + 30
		usrfunc= QtGui.QLabel('User Function Name', self)
		usrfunc.move(x1,y)
		self.usrfuncedit=QtGui.QLineEdit(self)
		self.usrfuncedit.move(x2,y)
		self.usrfuncedit.setText(self.savedparmsdict['usrfunc'])
		self.usrfuncedit.setToolTip('name of the user-supplied-function that prepares reference image for each iteration')
		
		y = y +30		
		usrfuncfile= QtGui.QLabel('Enter (full path) name of external file containing user function:', self)
		usrfuncfile.move(x1,y)
		y = y + 30
		usrfuncfile= QtGui.QLabel('(Leave blank if file is not external to Sparx)', self)		
		usrfuncfile.move(x1,y)
		y = y + 30
		self.usrfuncfileedit=QtGui.QLineEdit(self)
		self.usrfuncfileedit.move(x2,y)
		self.usrfuncfileedit.setText(self.savedparmsdict['usrfuncfile'])
		self.usrfuncfileedit.setToolTip('name of the external file containing user function')		
		#Function runsxali2d started when  the  RUN_button of the  Poptwodali window is clicked 
		#Here we define, that when this button is clicked, it starts subfunction choose_file
		self.usrfile_button = QtGui.QPushButton("Select File", self)
		self.usrfile_button.move(x3, y)
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
#helical_end				

class Popupadvparams_helical_2(QWidget):
	def __init__(self,savedparms):
		QWidget.__init__(self)
		#Here we just set the window title
		self.setWindowTitle('sxhelical advanced parameters related to helix and symmetry')
		#Here we just set a label and its position in the window
		x1 = 10
		x2 = x1 + 300
		x3 = x2 + 145
		x4 = x3 + 100
		y = 10
		title1=QtGui.QLabel('<b>sxhelical</b> - set advanced helical params', self)
		title1.move(10,y)
		#Labels and Line Edits for User Input
		#Just a label
		y = y +20
		title2= QtGui.QLabel('<b>Advanced</b> parameters', self)
		title2.move(10,y)

		self.savedparmsdict=savedparms
		y = y + 30
		nise = QtGui.QLabel('Number of iterations to start dp, dphi search', self)
		nise.move(x1,y)
		self.niseedit=QtGui.QLineEdit(self)
		self.niseedit.move(x2,y)
		self.niseedit.setText(self.savedparmsdict['nise'])
		self.niseedit.setToolTip('start symmetrization searching after nise steps')		
		
		y = y + 30
		rmin = QtGui.QLabel('Inner radius of helix', self)
		rmin.move(x1,y)
		self.rminedit=QtGui.QLineEdit(self)
		self.rminedit.move(x2,y)
		self.rminedit.setText(self.savedparmsdict['rmin'])
		self.rminedit.setToolTip('Inner radius of the helix')		
		
		y = y + 30
		fract = QtGui.QLabel('Fraction of the volume used for helicising', self)
		fract.move(x1,y)
		self.fractedit=QtGui.QLineEdit(self)
		self.fractedit.move(x2,y)
		self.fractedit.setText(self.savedparmsdict['fract'])
		self.fractedit.setToolTip('fraction of the volume used for helicising')		
		
		y = y + 30
		dp_step = QtGui.QLabel('Dz step size in Angstroms during search', self)
		dp_step.move(x1,y)
		self.dp_stepedit=QtGui.QLineEdit(self)
		self.dp_stepedit.move(x2,y)
		self.dp_stepedit.setText(self.savedparmsdict['dp_step'])
		self.dp_stepedit.setToolTip('step size of helicise rise search')
		
		y = y + 30
		ndp = QtGui.QLabel('Number of dz step during search', self)
		ndp.move(x1,y)
		self.ndpedit=QtGui.QLineEdit(self)
		self.ndpedit.move(x2,y)
		self.ndpedit.setText(self.savedparmsdict['ndp'])
		self.ndpedit.setToolTip('In symmetrization search, number of delta z steps equas to 2*ndp+1')		
		
		y = y + 30
		dphi_step = QtGui.QLabel('Dphi step size in degree during search', self)
		dphi_step.move(x1,y)
		self.dphi_stepedit=QtGui.QLineEdit(self)
		self.dphi_stepedit.move(x2,y)
		self.dphi_stepedit.setText(self.savedparmsdict['dphi_step'])
		self.dphi_stepedit.setToolTip('step size of helicise angle search')
		
		y = y + 30
		ndphi = QtGui.QLabel('Number of dphi step during search', self)
		ndphi.move(x1,y)
		self.ndphiedit=QtGui.QLineEdit(self)
		self.ndphiedit.move(x2,y)
		self.ndphiedit.setText(self.savedparmsdict['ndphi'])
		self.ndphiedit.setToolTip('In symmetrization search, number of angular steps equas to 2*ndphi+1')		
		
		y = y + 30
		psi_max = QtGui.QLabel('Maximum psi range to deviate from 90/270', self)
		psi_max.move(x1,y)
		self.psi_maxedit=QtGui.QLineEdit(self)
		self.psi_maxedit.move(x2,y)
		self.psi_maxedit.setText(self.savedparmsdict['psi_max'])
		self.psi_maxedit.setToolTip('maximum psi - how far rotation in plane can can deviate from 90 or 270 degrees')		
		

		y = y + 30
		sym = QtGui.QLabel('Point group symmetry', self)
		sym.move(x1,y)
		self.symedit=QtGui.QLineEdit(self)
		self.symedit.move(x2,y)
		self.symedit.setText(self.savedparmsdict['sym'])
		self.symedit.setToolTip('now support cn and dn group symmetry')
		
		y = y + 30
		datasym = QtGui.QLabel('File to save dp, dphi', self)
		datasym.move(x1,y)
		self.datasymedit=QtGui.QLineEdit(self)
		self.datasymedit.move(x2,y)
		self.datasymedit.setText(self.savedparmsdict['datasym'])
		self.datasymedit.setToolTip('file to save helical parameters of each iteration')		


"""
#Layout of the Pop Up window Popuptwodali (for sxali2d); started by the function twodali of the main window		
class Popuptwodali(QWidget):
	def __init__(self):
		QWidget.__init__(self)
		
		#######################################################################################
		# class variables
		
		self.cmd = ""
		# populate with default values
		self.savedparmsdict = {'stackname':'','foldername':'','partradius':'-1','xyrange':'4 2 1 1','trans':'2 1 0.5 0.25','nriter':'3','nproc':'1','maskname':'','center':'-1',"ringstep":"1","innerradius":"1","ctf":Qt.Unchecked,"snr":"1.0","fourvar":Qt.Unchecked, "gpnr":"-1","usrfunc":"ref_ali2d","usrfuncfile":""}

		#######################################################################################
		# Layout parameters
		
		self.y1 = 10 # title and Repopulate button
		self.y2 = self.y1 + 78 # Text boxes for inputting parameters
		self.y3 = self.y2 + 222 # activate images button and set xform.align2d button
		self.y4 = self.y3 + 80 # Advanced Parameters, Save Input and Generate command line buttons
		self.y5 = self.y4 + 95 # run button 
		self.yspc = 4
		
		self.x1 = 10 # first column (text box labels)
		self.x2 = self.x1 + 150 # second column (text boxes)
		self.x3 = self.x2+145 # third column (Open .hdf button)
		self.x4 = self.x3+100 # fourth column (Open .bdb button)
		self.x5 = 230 # run button
		#######################################################################################
		
		#Here we just set the window title
		self.setWindowTitle('sxali2d')
		#Here we just set a label and its position in the window
		title1=QtGui.QLabel('<b>sxali2d</b> - performs 2D reference free alignment of an image series', self)
		title1.move(self.x1,self.y1)
		self.y1 += 30

		self.repopbtn = QPushButton("Retrieve saved parameters", self)
		self.repopbtn.move(self.x1-5,self.y1)
		#sets an infotip for this Pushbutton
		self.repopbtn.setToolTip('Retrieve saved parameters')
		#when this button is clicked, this action starts the subfunction twodali
		self.connect(self.repopbtn, SIGNAL("clicked()"), self.repoparms_ali2d)

		#######################################################################################
		#Here we create a Button(file_button with title run open .hdf) and its position in the window
		self.file_button = QtGui.QPushButton("Open .hdf", self)
		self.file_button.move(self.x3, self.y2-self.yspc)
		#Here we define, that when this button is clicked, it starts subfunction choose_file
		QtCore.QObject.connect(self.file_button, QtCore.SIGNAL("clicked()"), self.choose_file)
		#exactly the same as above, but for subfunction choose_file1
		self.file_button1 = QtGui.QPushButton("Open .bdb", self)
		self.file_button1.move(self.x4,self.y2-self.yspc)
		QtCore.QObject.connect(self.file_button1, QtCore.SIGNAL("clicked()"), self.choose_file1)
		
		stackname= QtGui.QLabel('Name of input stack', self)
		stackname.move(self.x1,self.y2)
		self.stacknameedit=QtGui.QLineEdit(self)
		self.stacknameedit.move(self.x2,self.y2)
		self.stacknameedit.setText(self.savedparmsdict['stackname'])
		self.y2 += 30
		
		foldername= QtGui.QLabel('Output folder', self)
		foldername.move(self.x1,self.y2)
		self.foldernameedit=QtGui.QLineEdit(self)
		self.foldernameedit.move(self.x2,self.y2)
		self.foldernameedit.setText(self.savedparmsdict['foldername'])		
		
		self.outinfobtn = QPushButton("Output Info", self)
		self.outinfobtn.move(self.x3,  self.y2-self.yspc)
		#sets an infotip for this Pushbutton
		self.outinfobtn.setToolTip('Output Info')
		#when this button is clicked, this action starts the subfunction twodali
		self.connect(self.outinfobtn, SIGNAL("clicked()"), self.outputinfo_ali2d)
		self.y2 += 30
		
		partradius= QtGui.QLabel('Particle radius', self)
		partradius.move(self.x1,self.y2)
		self.partradiusedit=QtGui.QLineEdit(self)
		self.partradiusedit.move(self.x2,self.y2)
		self.partradiusedit.setText(self.savedparmsdict['partradius'])
		self.partradiusedit.setToolTip('Parameter ou: Outer radius for rotational correlation \nshould be set to particle radius\nif not sure, set to boxsize/2-2 ')		
		self.y2 += 30

		xyrange= QtGui.QLabel('xy range', self)
		xyrange.move(self.x1,self.y2)
		self.xyrangeedit=QtGui.QLineEdit(self)
		self.xyrangeedit.move(self.x2,self.y2)
		self.xyrangeedit.setText(self.savedparmsdict['xyrange'])
		self.xyrangeedit.setToolTip('Range for translational search in x, y direction\nif set to 0 only rotational alignment will be performed')
		self.y2 += 30
		
		trans= QtGui.QLabel('translational step', self)
		trans.move(self.x1,self.y2)
		self.transedit=QtGui.QLineEdit(self)
		self.transedit.move(self.x2,self.y2)
		self.transedit.setText(self.savedparmsdict['trans'])
		self.transedit.setToolTip('Step of translational search in x, y direction\nlarger values increase the speed but decrease the accuracy')		
		self.y2 += 30
		
		nriter= QtGui.QLabel('Number of iterations', self)
		nriter.move(self.x1,self.y2)
		self.nriteredit=QtGui.QLineEdit(self)
		self.nriteredit.move(self.x2,self.y2)
		self.nriteredit.setText(self.savedparmsdict['nriter'])
		self.nriteredit.setToolTip('Maximum number of iterations the program will perform\n Using the default values the program will run 3 rounds with xy-range 4 and translational step 1, 3 rounds with xyrange 2 and translational step 1 and so on..\nif set to 0 maximum iteration number will be 10 and will automatically stop should the criterion falls')
		self.y2 += 30
		
		nproc= QtGui.QLabel('MPI processors', self)
		nproc.move(self.x1,self.y2)
		self.nprocedit=QtGui.QLineEdit(self)
		self.nprocedit.move(self.x2,self.y2)
		self.nprocedit.setText(self.savedparmsdict['nproc'])
		self.nprocedit.setToolTip('The number of processors to use. Default is single processor mode')

		##########################################################
		
		header=QtGui.QLabel('Attributes xform.align2d and active parameters must be set in the input stack', self)
		header.move(self.x1,self.y3)
		self.y3 += 30
		
		#not linked to a function yet
		self.activeheader_button = QtGui.QPushButton("activate all images", self)
		self.activeheader_button.move(self.x1-5, self.y3)
		self.connect(self.activeheader_button, SIGNAL("clicked()"), self.setactiveheader)
		
		self.a2dheader_button = QtGui.QPushButton("set xform.align2d", self)
		self.a2dheader_button.move(self.x1-5+180, self.y3)
		self.connect(self.a2dheader_button, SIGNAL("clicked()"), self.seta2dheader)
		
		######################################################################################
		
		self.savepbtn = QPushButton("Save Input Parameters", self)
		self.savepbtn.move(self.x1-5, self.y4)
		#sets an infotip for this Pushbutton
		self.savepbtn.setToolTip('Save Input Parameters')
		#when this button is clicked, this action starts the subfunction twodali
		self.connect(self.savepbtn, SIGNAL("clicked()"), self.saveparms)
		self.y4+=30
		
		self.cmdlinebtn = QPushButton("Generate command line from input parameters", self)
		self.cmdlinebtn.move(self.x1-5, self.y4)
		#sets an infotip for this Pushbutton
		self.cmdlinebtn.setToolTip('Generate command line using input parameters')
		#when this button is clicked, this action starts the subfunction twodali
		self.connect(self.cmdlinebtn, SIGNAL("clicked()"), self.gencmdline_ali2d)

		#######################################################################
		 #Here we create a Button(Run_button with title run sxali2d) and its position in the window
		self.RUN_button = QtGui.QPushButton('Run sxali2d', self)
		# make 3D textured push button look
		s = "QPushButton {font: bold; color: #000;border: 1px solid #333;border-radius: 11px;padding: 2px;background: qradialgradient(cx: 0, cy: 0,fx: 0.5, fy:0.5,radius: 1, stop: 0 #fff, stop: 1 #8D0);min-width:90px;margin:5px} QPushButton:pressed {font: bold; color: #000;border: 1px solid #333;border-radius: 11px;padding: 2px;background: qradialgradient(cx: 0, cy: 0,fx: 0.5, fy:0.5,radius: 1, stop: 0 #fff, stop: 1 #084);min-width:90px;margin:5px}"
		
		self.RUN_button.setStyleSheet(s)
		self.RUN_button.move(self.x5, self.y5)
		#Here we define, that when this button is clicked, it starts subfunction runsxali2d
		self.connect(self.RUN_button, SIGNAL("clicked()"), self.runsxali2d)
		#Labels and Line Edits for User Input		

		
	 #Function runsxali2d started when  the  RUN_button of the  Poptwodali window is clicked 
	def outputinfo_ali2d(self):
			QMessageBox.information(self, "ali2d output",'Output files (average of aligned images and Fourier Ring Correlation curve)\nare saved in Output folder. alignment parameters are saved in the attribute \nxform.align2d in each image header. The images themselves are not changed.')
		
	def gencmdline_ali2d(self,writefile=True):
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
		
		mask = self.w1.masknameedit.text()
				
		if len(str(mask))> 1:
				cmd1 = cmd1+" "+str(mask) 
		cmd1 = cmd1 + args
		
		ctr=self.w1.centeredit.text()
		ringstep = self.w1.ringstepedit.text()
		inrad = self.w1.innerradiusedit.text()
		CTF=self.w1.ctfchkbx.checkState()
		snr = self.w1.snredit.text()
		fourvar = self.w1.fourvarchkbx.checkState()
		gpn = self.w1.gpnredit.text()
		userf = self.w1.usrfuncedit.text()
		userfile = self.w1.usrfuncfileedit.text()
						
		cmd1 = cmd1+" --center=" +str(ctr) +" --rs="+str(ringstep)+ " --ir=" + str(inrad)+ " --snr=" + str(snr)+ " --Ng=" + str(gpn)
		if CTF == Qt.Checked:
				cmd1 = cmd1 + " --CTF"
		if fourvar == Qt.Checked:
				cmd1 = cmd1 + " --Fourvar"		
		if len(userfile) < 1:
				cmd1 = cmd1 + " --function="+str(userf)
		else:
				userfile = str(userfile)
				# break it down into file name and directory path
				rind = userfile.rfind('/')
				
				if rind == -1:
						userfile = os.path.abspath(userfile)
						rind = userfile.rfind('/')
						
				fname = userfile[rind+1:]
				fname, ext = os.path.splitext(fname)
				fdir = userfile[0:rind]
				cmd1 = cmd1 + " --function=\"[" +fdir+","+fname+","+str(userf)+"]\""
				
		np = self.nprocedit.text()
		
		self.savedparmsdict = {'stackname':str(stack),'foldername':str(output),'partradius':str(ou),'xyrange':str(xr),'trans':str(yr),'nriter':str(maxit),'nproc':str(np),'maskname':str(mask),'center':str(ctr),"ringstep":str(ringstep),"innerradius":str(inrad),"ctf":CTF,"snr":str(snr),"fourvar":fourvar, "gpnr":str(gpn),"usrfunc":str(userf), "usrfuncfile":str(userfile)}
		
		self.w1.savedparmsdict=self.savedparmsdict
		
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
		self.gencmdline_ali2d(writefile=False)
		outfolder=self.savedparmsdict['foldername']
		if os.path.exists(outfolder):
			print "output folder "+outfolder+" already exists!"
			return
		process = subprocess.Popen(self.cmd,shell=True)
		self.emit(QtCore.SIGNAL("process_started"),process.pid)
		
	def saveparms(self):		
		(fname,stat)= QInputDialog.getText(self,"Save Input Parameters","Enter name of file to save parameters in",QLineEdit.Normal,"")
		if stat:
			import pickle
			output=open(fname,'wb')
			self.gencmdline_ali2d(writefile=False)
			pickle.dump(self.savedparmsdict,output)
			output.close()
		
	def repoparms_ali2d(self):		
		(fname,stat)= QInputDialog.getText(self,"Retrieve saved parameters","Enter name of file parameters were saved in",QLineEdit.Normal,"")
		if stat:
			import pickle
			pkl = open(fname,'rb')
			self.savedparmsdict = pickle.load(pkl)
			self.partradiusedit.setText(self.savedparmsdict['partradius'])
			self.stacknameedit.setText(self.savedparmsdict['stackname'])		
			self.foldernameedit.setText(self.savedparmsdict['foldername'])		
			self.xyrangeedit.setText(self.savedparmsdict['xyrange'])
			self.transedit.setText(self.savedparmsdict['trans'])
			self.nriteredit.setText(self.savedparmsdict['nriter'])
			self.nprocedit.setText(self.savedparmsdict['nproc'])
			self.w1.masknameedit.setText(self.savedparmsdict['maskname'])
			self.w1.centeredit.setText(self.savedparmsdict['center'])
			self.w1.ringstepedit.setText(self.savedparmsdict['ringstep'])
			self.w1.innerradiusedit.setText(self.savedparmsdict['innerradius'])
			self.w1.ctfchkbx.setCheckState(self.savedparmsdict['ctf'])
			self.w1.snredit.setText(self.savedparmsdict['snr'])
			self.w1.fourvarchkbx.setCheckState(self.savedparmsdict['fourvar'])
			self.w1.gpnredit.setText(self.savedparmsdict['gpnr'])
			self.w1.usrfuncedit.setText(self.savedparmsdict['usrfunc'])
			self.w1.usrfuncfileedit.setText(self.savedparmsdict['usrfuncfile'])
				
	
	def setactiveheader(self):
		stack = self.stacknameedit.text()
		print "stack defined="+ stack
		# horatio active_refactoring Jy51i1EwmLD4tWZ9_00000_1
		# header(str(stack), "active", one=True)


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
						if len(a)>0:
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
		
		self.x1 = 10
		self.x2 = 140
		self.x3 = 285
		
		self.y1 = 10
		self.yspc = 4
		
		#Here we just set the window title
		self.setWindowTitle('sxali2d advanced parameter selection')
		#Here we just set a label and its position in the window
		title1=QtGui.QLabel('<b>sxali2d</b> - set advanced params', self)
		title1.move(self.x1,self.y1)
		#Labels and Line Edits for User Input
		#Just a label
		self.y1 += 30
		
		title2= QtGui.QLabel('<b>Advanced</b> parameters', self)
		title2.move(self.x1,self.y1)
		
		self.y1 += 30
		
		self.savedparmsdict=savedparms
		#Example for User input stack name
		#First create the label and define its position
		maskname= QtGui.QLabel('Mask', self)
		maskname.move(self.x1,self.y1)
		#Now add a line edit and define its position
		self.masknameedit=QtGui.QLineEdit(self)
		self.masknameedit.move(self.x2,self.y1)
		#Adds a default value for the line edit
		self.masknameedit.setText(self.savedparmsdict['maskname'])
		self.masknameedit.setToolTip("Default is a circle mask with radius equal to the particle radius")
		
		self.mskfile_button = QtGui.QPushButton("Open File", self)
		self.mskfile_button.move(self.x3, self.y1-self.yspc)
		#Here we define, that when this button is clicked, it starts subfunction choose_file
		QtCore.QObject.connect(self.mskfile_button, QtCore.SIGNAL("clicked()"), self.choose_mskfile)
		
		self.y1 += 30
		
		center= QtGui.QLabel('Center type', self)
		center.move(self.x1,self.y1)
		self.centeredit=QtGui.QLineEdit(self)
		self.centeredit.move(self.x2,self.y1)
		self.centeredit.setText(self.savedparmsdict['center'])
		self.centeredit.setToolTip('-1 - use average centering method (default),\n0 - if you do not want the average to be centered, \n1 - phase approximation of the center of gravity phase_cog, \n2 - cross-correlate with Gaussian function, \n3 - cross-correlate with donut shape image (e.g. inner radius=2, outer radius=7), \n4 - cross-correlate with reference image provided by user, \n5 - cross-correlate with self-rotated average..\ncentering may fail..use 0 to deactive it')
		
		self.y1 += 30
		
		ringstep= QtGui.QLabel('Ring step', self)
		ringstep.move(self.x1,self.y1)
		self.ringstepedit=QtGui.QLineEdit(self)
		self.ringstepedit.move(self.x2,self.y1)
		self.ringstepedit.setText(self.savedparmsdict['ringstep'])
		self.ringstepedit.setToolTip('step between rings in rotational correlation > 0 (set to 1)')
		
		self.y1 += 30
		
		innerradius= QtGui.QLabel('Inner radius', self)
		innerradius.move(self.x1,self.y1)
		self.innerradiusedit=QtGui.QLineEdit(self)
		self.innerradiusedit.move(self.x2,self.y1)
		self.innerradiusedit.setText(self.savedparmsdict['innerradius'])
		self.innerradiusedit.setToolTip('inner radius for rotational correlation > 0 (set to 1) ')		
		
		self.y1 += 30
		
		ctf= QtGui.QLabel('CTF', self)
		ctf.move(self.x1,self.y1)
		self.ctfchkbx = QtGui.QCheckBox("",self)
		self.ctfchkbx.move(self.x2, self.y1)
		self.ctfchkbx.setCheckState(self.savedparmsdict['ctf'])
		
		self.y1 += 30
		
		snr= QtGui.QLabel('SNR', self)
		snr.move(self.x1,self.y1)
		self.snredit=QtGui.QLineEdit(self)
		self.snredit.move(self.x2,self.y1)
		self.snredit.setText(self.savedparmsdict['snr'])
		self.snredit.setToolTip('signal-to-noise ratio of the data (default SNR=1.0)')		
		
		self.y1 += 30
				
		fourvar= QtGui.QLabel('Fourvar', self)
		fourvar.move(self.x1,self.y1)
		self.fourvarchkbx=QtGui.QCheckBox("",self)
		self.fourvarchkbx.move(self.x2,self.y1)
		self.fourvarchkbx.setCheckState(self.savedparmsdict['fourvar'])
		self.fourvarchkbx.setToolTip('use Fourier variance to weight the reference (recommended, default False)')
		
		self.y1 += 30
		
		gpnr= QtGui.QLabel('Number of Groups', self)
		gpnr.move(self.x1,self.y1)
		self.gpnredit=QtGui.QLineEdit(self)
		self.gpnredit.move(self.x2,self.y1)
		self.gpnredit.setText(self.savedparmsdict['gpnr'])
		self.gpnredit.setToolTip('number of groups in the new CTF filteration')		
		
		self.y1 += 30
		
		usrfunc= QtGui.QLabel('User Function Name', self)
		usrfunc.move(self.x1,self.y1)
		self.usrfuncedit=QtGui.QLineEdit(self)
		self.usrfuncedit.move(self.x2,self.y1)
		self.usrfuncedit.setText(self.savedparmsdict['usrfunc'])
		self.usrfuncedit.setToolTip('name of the user-supplied-function that prepares reference image for each iteration')
		
		self.y1 += 30
				
		usrfuncfile= QtGui.QLabel('Enter name of external file containing user function:', self)
		usrfuncfile.move(self.x1,self.y1)
		
		self.y1 += 20
		
		usrfuncfile= QtGui.QLabel('(Leave blank if file is not external to Sparx)', self)
		usrfuncfile.move(self.x1,self.y1)
		
		self.y1 += 20
		
		self.usrfuncfileedit=QtGui.QLineEdit(self)
		self.usrfuncfileedit.move(self.x2,self.y1)
		self.usrfuncfileedit.setText(self.savedparmsdict['usrfuncfile'])
		self.usrfuncfileedit.setToolTip('name of the external file containing user function')		
	 #Function runsxali2d started when  the  RUN_button of the  Poptwodali window is clicked 
			
		#Here we define, that when this button is clicked, it starts subfunction choose_file
		self.usrfile_button = QtGui.QPushButton("Select File", self)
		self.usrfile_button.move(self.x3, self.y1-self.yspc)
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
"""				

#Layout of the Pop Up window Popuptwodali (for sxmeridien); started by the function twodali of the main window		
class Popupmeridien(QWidget):
	def __init__(self):
		QWidget.__init__(self)
		
		#######################################################################################
		# class variables
		
		self.cmd = ""
		# populate with default values
		self.savedparmsdict = {'stackname':'','foldername':'','inivol':'','radius':'-1','sym':'c1','inires':'25','pwreference':'','nproc':'1','mask':'',"CTF":Qt.Unchecked,"usrfunc":"","usrfuncfile":""}

		#######################################################################################
		# Layout parameters
		
		self.y1 = 10 # title and Repopulate button
		self.y2 = self.y1 + 78 # Text boxes for inputting parameters
		self.y3 = self.y2 + 222 # activate images button and set xform.align2d button
		self.y4 = self.y3 + 80 # Advanced Parameters, Save Input and Generate command line buttons
		self.y5 = self.y4 + 95 # run button 
		self.yspc = 4
		
		self.x1 = 10 # first column (text box labels)
		self.x2 = self.x1 + 150 # second column (text boxes)
		self.x3 = self.x2+145 # third column (Open .hdf button)
		self.x4 = self.x3+100 # fourth column (Open .bdb button)
		self.x5 = 230 # run button
		#######################################################################################
		
		#Here we just set the window title
		self.setWindowTitle('sxmeridien')
		#Here we just set a label and its position in the window
		title1=QtGui.QLabel('<b>sxmeridien</b> - performs 3D structure refinement', self)
		title1.move(self.x1,self.y1)
		self.y1 += 30

		self.repopbtn = QPushButton("Retrieve saved parameters", self)
		self.repopbtn.move(self.x1-5,self.y1)
		#sets an infotip for this Pushbutton
		self.repopbtn.setToolTip('Retrieve saved parameters')
		#when this button is clicked, this action starts the subfunction twodali
		self.connect(self.repopbtn, SIGNAL("clicked()"), self.repoparms_meridien)

		#######################################################################################
		#Here we create a Button(file_button with title run open .hdf) and its position in the window
		self.file_button = QtGui.QPushButton("Open .hdf", self)
		self.file_button.move(self.x3, self.y2-self.yspc)
		#Here we define, that when this button is clicked, it starts subfunction choose_file
		QtCore.QObject.connect(self.file_button, QtCore.SIGNAL("clicked()"), self.choose_file)
		#exactly the same as above, but for subfunction choose_file1
		self.file_button1 = QtGui.QPushButton("Open .bdb", self)
		self.file_button1.move(self.x4,self.y2-self.yspc)
		QtCore.QObject.connect(self.file_button1, QtCore.SIGNAL("clicked()"), self.choose_file1)
		
		stackname= QtGui.QLabel('Name of input stack', self)
		stackname.move(self.x1,self.y2)
		self.stacknameedit=QtGui.QLineEdit(self)
		self.stacknameedit.move(self.x2,self.y2)
		self.stacknameedit.setText(self.savedparmsdict['stackname'])
		self.y2 += 30
		
		foldername= QtGui.QLabel('Output folder', self)
		foldername.move(self.x1,self.y2)
		self.foldernameedit=QtGui.QLineEdit(self)
		self.foldernameedit.move(self.x2,self.y2)
		self.foldernameedit.setText(self.savedparmsdict['foldername'])		
				
		self.outinfobtn = QPushButton("Output Info", self)
		self.outinfobtn.move(self.x3,  self.y2-self.yspc)
		#sets an infotip for this Pushbutton
		self.outinfobtn.setToolTip('Output Info')
		#when this button is clicked, this action starts the subfunction twodali
		self.connect(self.outinfobtn, SIGNAL("clicked()"), self.outputinfo_meridien)
		self.y2 += 30

		inivol= QtGui.QLabel('Initial 3D structure', self)
		inivol.move(self.x1,self.y2)
		self.inivoledit=QtGui.QLineEdit(self)
		self.inivoledit.move(self.x2,self.y2)
		self.inivoledit.setText(self.savedparmsdict['inivol'])
		self.y2 += 30
		
		radius= QtGui.QLabel('Particle radius', self)
		radius.move(self.x1,self.y2)
		self.radiusedit=QtGui.QLineEdit(self)
		self.radiusedit.move(self.x2,self.y2)
		self.radiusedit.setText(self.savedparmsdict['radius'])
		self.radiusedit.setToolTip('Parameter radius: Radius of the structure in pixels\nif not sure, set to boxsize/2-2')		
		self.y2 += 30

		pwreference = QtGui.QLabel('Power spectrum', self)
		pwreference.move(self.x1,self.y2)
		self.pwreferenceedit=QtGui.QLineEdit(self)
		self.pwreferenceedit.move(self.x2,self.y2)
		self.pwreferenceedit.setText(self.savedparmsdict['pwreference'])
		self.pwreferenceedit.setToolTip('1D reference power spectrum of the structure, can be generated using sxprocess.py')		
		self.y2 += 30

		mask = QtGui.QLabel('3D mask', self)
		mask.move(self.x1,self.y2)
		self.maskedit=QtGui.QLineEdit(self)
		self.maskedit.move(self.x2,self.y2)
		self.maskedit.setText(self.savedparmsdict['mask'])
		self.maskedit.setToolTip('3D mask that defines outline of the structure, preferable with soft edges\nif not given, set to spherical mask with radius boxsize/2-1.')		
		self.y2 += 30

		sym = QtGui.QLabel('Point-group symmetry', self)
		sym.move(self.x1,self.y2)
		self.symedit=QtGui.QLineEdit(self)
		self.symedit.move(self.x2,self.y2)
		self.symedit.setText(self.savedparmsdict['sym'])
		self.symedit.setToolTip('Point-group symmetry of the structure: cn, dn, where n is multiplicity (for example c5 or d3)')		
		self.y2 += 30

		inires = QtGui.QLabel('Initial resolution', self)
		inires.move(self.x1,self.y2)
		self.iniresedit=QtGui.QLineEdit(self)
		self.iniresedit.move(self.x2,self.y2)
		self.iniresedit.setText(self.savedparmsdict['inires'])
		self.iniresedit.setToolTip('Initial resolution of the initial 3D structure (default 25A)')		
		self.y2 += 30

		ctf = QtGui.QLabel('CTF', self)
		ctf.move(self.x1,self.y2)
		self.ctfchkbx = QtGui.QCheckBox("",self)
		self.ctfchkbx.move(self.x2, self.y2)
		self.ctfchkbx.setCheckState(self.savedparmsdict['CTF'])		
		self.y2 += 30

		nproc = QtGui.QLabel('MPI processes', self)
		nproc.move(self.x1,self.y2)
		self.nprocedit=QtGui.QLineEdit(self)
		self.nprocedit.move(self.x2,self.y2)
		self.nprocedit.setText(self.savedparmsdict['nproc'])
		self.nprocedit.setToolTip('The number of processors to use. Default is one processor')

		######################################################################################
		
		self.savepbtn = QPushButton("Save Input Parameters", self)
		self.savepbtn.move(self.x1-5, self.y4)
		#sets an infotip for this Pushbutton
		self.savepbtn.setToolTip('Save Input Parameters')
		#when this button is clicked, this action starts the subfunction twodali
		self.connect(self.savepbtn, SIGNAL("clicked()"), self.saveparms)
		self.y4+=30
		
		self.cmdlinebtn = QPushButton("Generate command line from input parameters", self)
		self.cmdlinebtn.move(self.x1-5, self.y4)
		#sets an infotip for this Pushbutton
		self.cmdlinebtn.setToolTip('Generate command line using input parameters')
		#when this button is clicked, this action starts the subfunction twodali
		self.connect(self.cmdlinebtn, SIGNAL("clicked()"), self.gencmdline_meridien)

		#######################################################################
		#Here we create a Button(Run_button with title run sxmeridien) and its position in the window
		self.RUN_button = QtGui.QPushButton('Run sxmeridien', self)
		# make 3D textured push button look
		s = "QPushButton {font: bold; color: #000;border: 1px solid #333;border-radius: 11px;padding: 2px;background: qradialgradient(cx: 0, cy: 0,fx: 0.5, fy:0.5,radius: 1, stop: 0 #fff, stop: 1 #8D0);min-width:90px;margin:5px} QPushButton:pressed {font: bold; color: #000;border: 1px solid #333;border-radius: 11px;padding: 2px;background: qradialgradient(cx: 0, cy: 0,fx: 0.5, fy:0.5,radius: 1, stop: 0 #fff, stop: 1 #084);min-width:90px;margin:5px}"
		
		self.RUN_button.setStyleSheet(s)
		self.RUN_button.move(self.x5, self.y5)
		#Here we define, that when this button is clicked, it starts subfunction runsxmeridien
		self.connect(self.RUN_button, SIGNAL("clicked()"), self.runsxmeridien)
		#Labels and Line Edits for User Input		

		
		#Function runsxmeridien started when  the  RUN_button of the  Poptwodali window is clicked 
	def outputinfo_meridien(self):
		QMessageBox.information(self, "meridien output",'Output files (average of aligned images and Fourier Ring Correlation curve)\nare saved in Output folder. alignment parameters are saved in the attribute \nxform.align2d in each image header. The images themselves are not changed.')
		
	def gencmdline_meridien(self,writefile=True):
		#Here we just read in all user inputs in the line edits of the Poptwodali window
		stack       = self.stacknameedit.text()
		output      = self.foldernameedit.text()
		radius      = self.radiusedit.text()
		inires      = self.iniresedit.text()
		inivol      = self.inivoledit.text()
		sym         = self.symedit.text()
		pwreference = self.pwreferenceedit.text()
		mask        = self.maskedit.text()
		
		cmd1 = "sxmeridien.py "+str(stack) +" "+ str(output)+" "+ str(inivol)
		
		args = " --radius=" + str(radius) +  " --sym=" + str(sym) + " --smear"

		cmd1 += args
		if len(str(mask)) > 0:
				cmd1 += " --mask3D=" +str(mask)
		if len(str(pwreference)) > 0:
			cmd1 += " --pwreference=" + str(pwreference)
		if len(str(inires)) > 0:
			cmd1 += " --inires=" + str(inires)

		CTF        = self.ctfchkbx.checkState()
		userf      = self.w1.usrfuncedit.text()
		userfile   = self.w1.usrfuncfileedit.text()
						
		if CTF == Qt.Checked:
				cmd1 = cmd1 + " --CTF"
		if len(userfile) < 1:
				if(len(userf) > 0):   cmd1 += " --function="+str(userf)
		else:
			userfile = str(userfile)
			# break it down into file name and directory path
			rind = userfile.rfind('/')
			
			if rind == -1:
				userfile = os.path.abspath(userfile)
				rind = userfile.rfind('/')

			fname = userfile[rind+1:]
			fname, ext = os.path.splitext(fname)
			fdir = userfile[0:rind]
			cmd1 += " --function=\"[" +fdir+","+fname+","+str(userf)+"]\""
				
		nproc = self.nprocedit.text()
		self.savedparmsdict = {'stackname':str(stack),'foldername':str(output),'inivol':str(inivol),\
								'radius':str(radius),'sym':str(sym),'inires':str(inires),'nproc':str(nproc),'pwreference':str(pwreference),\
								'maskname':str(mask),"ctf":CTF,"usrfunc":str(userf), "usrfuncfile":str(userfile)}

		self.w1.savedparmsdict = self.savedparmsdict
		
		if int(str(nproc)) > 0:
				cmd1 = "mpirun -np " + str(nproc) + " " + cmd1
		else:  ERROR("sxgui","numper of processors has to be specified",1)

		if writefile:		
			(fname,stat)= QInputDialog.getText(self,"Generate Command Line","Enter name of file to save command line in",QLineEdit.Normal,"")
			if stat:
				f = open(fname,'a')
				f.write(cmd1)
				f.write('\n')
				f.close()

		print cmd1
		self.cmd = cmd1

	def runsxmeridien(self):
		self.gencmdline_meridien(writefile=False)
		outfolder=self.savedparmsdict['foldername']

		process = subprocess.Popen(self.cmd,shell=True)
		self.emit(QtCore.SIGNAL("process_started"),process.pid)
		
	def saveparms(self):		
		(fname,stat)= QInputDialog.getText(self,"Save Input Parameters","Enter name of file to save parameters in",QLineEdit.Normal,"")
		if stat:
			import pickle
			output=open(fname,'wb')
			self.gencmdline_meridien(writefile=False)
			pickle.dump(self.savedparmsdict,output)
			output.close()

	def repoparms_meridien(self):		
		(fname,stat)= QInputDialog.getText(self,"Retrieve saved parameters","Enter name of file parameters were saved in",QLineEdit.Normal,"")
		if stat:
			import pickle
			pkl = open(fname,'rb')
			self.savedparmsdict = pickle.load(pkl)
			self.radiusedit.setText(self.savedparmsdict['radius'])
			self.stacknameedit.setText(self.savedparmsdict['stackname'])		
			self.foldernameedit.setText(self.savedparmsdict['foldername'])		
			self.nprocedit.setText(self.savedparmsdict['nproc'])
			self.maskedit.setText(self.savedparmsdict['maskname'])
			self.iniresedit.setText(self.savedparmsdict['inires'])	
			self.inivoledit.setText(self.savedparmsdict['inivol'])	
			self.symedit.setText(self.savedparmsdict['sym'])	
			self.pwreferenceedit.setText(self.savedparmsdict['pwreference'])	
			self.ctfchkbx.setCheckState(self.savedparmsdict['ctf'])
			self.w1.usrfuncedit.setText(self.savedparmsdict['usrfunc'])
			self.w1.usrfuncfileedit.setText(self.savedparmsdict['usrfuncfile'])


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


class Popupadvparams_meridien(QWidget):
	def __init__(self,savedparms):
		QWidget.__init__(self)
		
		self.x1 = 10
		self.x2 = 140
		self.x3 = 285
		
		self.y1 = 10
		self.yspc = 4
		
		#Here we just set the window title
		self.setWindowTitle('sxmeridien advanced parameter selection')
		#Here we just set a label and its position in the window
		title1=QtGui.QLabel('<b>sxmeridien</b> - set advanced params', self)
		title1.move(self.x1,self.y1)
		#Labels and Line Edits for User Input
		#Just a label
		self.y1 += 30
		
		title2= QtGui.QLabel('<b>Advanced</b> parameters', self)
		title2.move(self.x1,self.y1)
		
		self.y1 += 30
		
		self.savedparmsdict=savedparms
		#Example for User input stack name
		#First create the label and define its position

		usrfunc= QtGui.QLabel('User Function Name', self)
		usrfunc.move(self.x1,self.y1)
		self.usrfuncedit=QtGui.QLineEdit(self)
		self.usrfuncedit.move(self.x2,self.y1)
		self.usrfuncedit.setText(self.savedparmsdict['usrfunc'])
		self.usrfuncedit.setToolTip('name of the user-supplied-function that prepares reference image for each iteration')

		self.y1 += 30
		
		usrfuncfile= QtGui.QLabel('Enter name of external file containing user function:', self)
		usrfuncfile.move(self.x1,self.y1)

		self.y1 += 20

		usrfuncfile= QtGui.QLabel('(Leave blank if file is not external to Sparx)', self)
		usrfuncfile.move(self.x1,self.y1)

		self.y1 += 20

		self.usrfuncfileedit=QtGui.QLineEdit(self)
		self.usrfuncfileedit.move(self.x2,self.y1)
		self.usrfuncfileedit.setText(self.savedparmsdict['usrfuncfile'])
		self.usrfuncfileedit.setToolTip('name of the external file containing user function')		
		#Function runsxmeridien started when  the  RUN_button of the  Poptwodali window is clicked 

		#Here we define, that when this button is clicked, it starts subfunction choose_file
		self.usrfile_button = QtGui.QPushButton("Select File", self)
		self.usrfile_button.move(self.x3, self.y1-self.yspc)
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


class Popupthreedali(QWidget):
	def __init__(self):
		QWidget.__init__(self)
		
		########################################################################################
		# class variables
		self.cmd = ""
		# populate with default values
		self.savedparmsdict = {'stackname':'','refname':'','foldername':'','partradius':'-1','xyrange':'4 2 1 1','trans':'2 1 0.5 0.25', 'delta':'15 5 2','nriter':'3','nproc':'1','maskname':'','center':'-1',"ringstep":"1","innerradius":"1","ctf":Qt.Unchecked,"snr":"1.0","fourvar":Qt.Unchecked,"usrfunc":"ref_ali3d","usrfuncfile":"","an":'-1',"ref_a":'S',"sym":'c1',"npad":'4',"deltapsi":'-1',"startpsi":'-1',"stoprnct":'0.0',"debug":False,"nch":False,"chunk":'0.2',"rantest":False}
		
		if options.demo == DEMO_mpibdbctf:
			self.savedparmsdict['stackname']='bdb:data'
			self.savedparmsdict['refname']='ab_initio_volf.hdf'
			self.savedparmsdict['foldername']='ali3d_d'
			self.savedparmsdict['partradius']='30'
			self.savedparmsdict['xyrange']='2.0 1.0 0.5'
			self.savedparmsdict['trans']='1.0 1 0.5'
			self.savedparmsdict['delta']='31 10 2'
			self.savedparmsdict['an']='-1 12 6'
			self.savedparmsdict['nriter']='2'
			self.savedparmsdict['ref_a']='P'
			self.savedparmsdict['usrfunc']='reference3'
			self.savedparmsdict['ctf']=Qt.Checked
			self.savedparmsdict['nproc']='4'
		
		if options.demo == DEMO_mpibdb:
			self.savedparmsdict['stackname']='bdb:data'
			self.savedparmsdict['refname']='ab_initio_volf.hdf'
			self.savedparmsdict['foldername']='ali3d_d'
			self.savedparmsdict['partradius']='30'
			self.savedparmsdict['xyrange']='2.0 1.0 0.5'
			self.savedparmsdict['trans']='1.0 1 0.5'
			self.savedparmsdict['delta']='31 10 2'
			self.savedparmsdict['an']='-1 12 6'
			self.savedparmsdict['nriter']='2'
			self.savedparmsdict['ref_a']='P'
			self.savedparmsdict['usrfunc']='reference3'
			self.savedparmsdict['nproc']='4'
		########################################################################################
		# layout parameters
		
		self.y1 = 10
		self.y2 = self.y1 + 78
		self.y3 = self.y2 + 282
		self.y4 = self.y3 + 80
		self.y5 = self.y4 + 95
		self.yspc = 4
		
		self.x1 = 10
		self.x2 = self.x1 + 150
		self.x3 = self.x2+145
		self.x4 = self.x3+100
		self.x5 = 230
		########################################################################################
		
		#Here we just set the window title
		self.setWindowTitle('sxali3d')
		#Here we just set a label and its position in the window
		title1=QtGui.QLabel('<b>sxali3d</b> - 3D projection matching given reference structure and an image series', self)
		title1.move(self.x1,self.y1)
		self.y1 += 30

		self.repopbtn = QPushButton("Retrieve saved parameters", self)
		self.repopbtn.move(self.x1-5,self.y1)
		#sets an infotip for this Pushbutton
		self.repopbtn.setToolTip('Retrieve saved parameters')
		#when this button is clicked, this action starts the subfunction twodali
		self.connect(self.repopbtn, SIGNAL("clicked()"), self.repoparms_ali3d)
		
		########################################################################################
		
		stackname= QtGui.QLabel('Name of input stack', self)
		stackname.move(self.x1,self.y2)
		#Now add a line edit and define its position
		self.stacknameedit=QtGui.QLineEdit(self)
		self.stacknameedit.move(self.x2,self.y2)
		#Adds a default value for the line edit
		self.stacknameedit.setText(self.savedparmsdict['stackname'])
		
		#Here we create a Button(file_button with title run open .hdf) and its position in the window
		self.file_button = QtGui.QPushButton("Open .hdf", self)
		self.file_button.move(self.x3, self.y2-self.yspc)
		#Here we define, that when this button is clicked, it starts subfunction choose_file
		QtCore.QObject.connect(self.file_button, QtCore.SIGNAL("clicked()"), self.choose_file)
		#exactly the same as above, but for subfunction choose_file1
		self.file_button1 = QtGui.QPushButton("Open .bdb", self)
		self.file_button1.move(self.x4,self.y2-self.yspc)
		QtCore.QObject.connect(self.file_button1, QtCore.SIGNAL("clicked()"), self.choose_file1)
		
		self.y2 += 30
				
		refname= QtGui.QLabel('Name of reference', self)
		refname.move(self.x1,self.y2)
		#Now add a line edit and define its position
		self.refnameedit=QtGui.QLineEdit(self)
		self.refnameedit.move(self.x2,self.y2)
		#Adds a default value for the line edit
		self.refnameedit.setText(self.savedparmsdict['refname'])
		#Here we create a Button(file_button with title run open .hdf) and its position in the window
		self.ref_button = QtGui.QPushButton("Open files", self)
		self.ref_button.move(self.x3, self.y2-self.yspc)
		#Here we define, that when this button is clicked, it starts subfunction choose_file
		QtCore.QObject.connect(self.ref_button, QtCore.SIGNAL("clicked()"), self.choose_reffile)

		
		self.y2 = self.y2+30
		
		#The same as above, but many line edits include Infotips
		foldername= QtGui.QLabel('Output folder', self)
		foldername.move(self.x1,self.y2)
		self.foldernameedit=QtGui.QLineEdit(self)
		self.foldernameedit.move(self.x2,self.y2)
		self.foldernameedit.setText(self.savedparmsdict['foldername'])		
		
		self.outinfobtn = QPushButton("Output Info", self)
		self.outinfobtn.move(self.x3,  self.y2-self.yspc)
		#sets an infotip for this Pushbutton
		self.outinfobtn.setToolTip('Output Info')
		#when this button is clicked, this action starts the subfunction twodali
		self.connect(self.outinfobtn, SIGNAL("clicked()"), self.outputinfo_ali3d)
		
		self.y2 = self.y2+30
		
		partradius= QtGui.QLabel('Particle radius', self)
		partradius.move(self.x1,self.y2)
		self.partradiusedit=QtGui.QLineEdit(self)
		self.partradiusedit.move(self.x2,self.y2)
		self.partradiusedit.setText(self.savedparmsdict['partradius'])
		self.partradiusedit.setToolTip('Parameter ou: Outer radius for rotational correlation \nshould be set to particle radius\nif not sure, set to boxsize/2-2 ')		

		self.y2 = self.y2+30
		
		xyrange= QtGui.QLabel('xy range', self)
		xyrange.move(self.x1,self.y2)
		self.xyrangeedit=QtGui.QLineEdit(self)
		self.xyrangeedit.move(self.x2,self.y2)
		self.xyrangeedit.setText(self.savedparmsdict['xyrange'])
		self.xyrangeedit.setToolTip('Range for translational search in x, y direction\nif set to 0 only rotational alignment will be performed')

		self.y2 = self.y2+30
		
		trans= QtGui.QLabel('translational step', self)
		trans.move(self.x1,self.y2)
		self.transedit=QtGui.QLineEdit(self)
		self.transedit.move(self.x2,self.y2)
		self.transedit.setText(self.savedparmsdict['trans'])
		self.transedit.setToolTip('Step of translational search in x, y direction\nlarger values increase the speed but decrease the accuracy')		

		self.y2 = self.y2+30
		
		delta= QtGui.QLabel('angular step', self)
		delta.move(self.x1,self.y2)
		self.deltaedit=QtGui.QLineEdit(self)
		self.deltaedit.move(self.x2,self.y2)
		self.deltaedit.setText(self.savedparmsdict['delta'])
		self.deltaedit.setToolTip('angular step for the reference projections in respective iterations')				
		
		self.y2 =self.y2+30

		nriter= QtGui.QLabel('Number of iterations', self)
		nriter.move(self.x1,self.y2)
		self.nriteredit=QtGui.QLineEdit(self)
		self.nriteredit.move(self.x2,self.y2)
		self.nriteredit.setText(self.savedparmsdict['nriter'])
		self.nriteredit.setToolTip('Maximum number of iterations the program will perform\n Using the default values the program will run 3 rounds with xy-range 4 and translational step 1, 3 rounds with xyrange 2 and translational step 1 and so on..\nif set to 0 maximum iteration number will be 10 and will automatically stop should the criterion falls')
		
		self.y2 =self.y2+30
		
		nproc= QtGui.QLabel('MPI processors', self)
		nproc.move(self.x1,self.y2)
		self.nprocedit=QtGui.QLineEdit(self)
		self.nprocedit.move(self.x2,self.y2)
		self.nprocedit.setText(self.savedparmsdict['nproc'])
		self.nprocedit.setToolTip('The number of processors to use. Default is single processor mode')
		
		########################################################################################
				
		header=QtGui.QLabel('Attributes xform.projection and active parameters must be set in the input stack', self)
		header.move(self.x1,self.y3)
		self.y3 = self.y3 + 30
		
		#not linked to a function yet
		self.activeheader_button = QtGui.QPushButton("activate all images", self)
		self.activeheader_button.move(self.x1-5, self.y3)
		self.connect(self.activeheader_button, SIGNAL("clicked()"), self.setactiveheader)
		
		self.projheader_button = QtGui.QPushButton("set xform.projection", self)
		self.projheader_button.move(self.x1-5+180, self.y3)
		self.connect(self.projheader_button, SIGNAL("clicked()"), self.setprojheader)
		
		########################################################################################
		
		self.savepbtn = QPushButton("Save Input Parameters", self)
		self.savepbtn.move(self.x1-5,  self.y4)
		#sets an infotip for this Pushbutton
		self.savepbtn.setToolTip('Save Input Parameters')
		#when this button is clicked, this action starts the subfunction twodali
		self.connect(self.savepbtn, SIGNAL("clicked()"), self.saveparms)
		
		self.y4 = self.y4+30
		
		self.cmdlinebtn = QPushButton("Generate command line from input parameters", self)
		self.cmdlinebtn.move(self.x1-5,  self.y4)
		#sets an infotip for this Pushbutton
		self.cmdlinebtn.setToolTip('Generate command line using input parameters')
		#when this button is clicked, this action starts the subfunction twodali
		self.connect(self.cmdlinebtn, SIGNAL("clicked()"), self.gencmdline_ali3d)
		
		########################################################################################
		
		 #Here we create a Button(Run_button with title run sxali2d) and its position in the window
		self.RUN_button = QtGui.QPushButton('Run sxali3d', self)
		# make 3D textured push button look
		s = "QPushButton {font: bold; color: #000;border: 1px solid #333;border-radius: 11px;padding: 2px;background: qradialgradient(cx: 0, cy: 0,fx: 0.5, fy:0.5,radius: 1, stop: 0 #fff, stop: 1 #8D0);min-width:90px;margin:5px} QPushButton:pressed {font: bold; color: #000;border: 1px solid #333;border-radius: 11px;padding: 2px;background: qradialgradient(cx: 0, cy: 0,fx: 0.5, fy:0.5,radius: 1, stop: 0 #fff, stop: 1 #084);min-width:90px;margin:5px}"
		
		self.RUN_button.setStyleSheet(s)
		
		self.RUN_button.move(self.x5,  self.y5)
		#Here we define, that when this button is clicked, it starts subfunction runsxali2d
		self.connect(self.RUN_button, SIGNAL("clicked()"), self.runsxali3d)
		#Labels and Line Edits for User Input

		
	 #Function runsxali2d started when  the  RUN_button of the  Poptwodali window is clicked 

	def outputinfo_ali3d(self):
			QMessageBox.information(self, "ali3d output",'Output volumes and Fourier Shell Criterion curves are saved in Output folder. Projection parameters are saved in the attribute xform.projection in image headers. The images themselves are not changed.')
		
	def gencmdline_ali3d(self,writefile=True):
		#Here we just read in all user inputs in the line edits of the Poptwodali window
		stack = self.stacknameedit.text()
		print "stack defined="+ stack
		ref = self.refnameedit.text()
		print "reference structure ="+ stack
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
		delta=self.deltaedit.text()
		print "delta=" +delta
		maxit=self.nriteredit.text()
		print "maxit="+maxit
		
		cmd1 = "sxali3d.py "+str(stack) +" "+ str(ref)+" "+ str(output)
		
		args = " --ou="+ str(ou)+ " --xr='"+str(xr)+"'"+ " --yr='"+str(yr)+"'"+ " --ts='"+str(ts)+"'"+ " --delta='"+str(delta)+"'"+" --maxit="+ str(maxit) 
		
		
		mask = self.w1.masknameedit.text()
		if len(str(mask))> 1:
				cmd1 = cmd1+" "+str(mask) 
		cmd1 = cmd1 + args
		
		
		ctr=self.w1.centeredit.text()
		ringstep = self.w1.ringstepedit.text()
		inrad = self.w1.innerradiusedit.text()
		CTF=self.w1.ctfchkbx.checkState()
		snr = self.w1.snredit.text()
		userf = self.w1.usrfuncedit.text()
		userfile = self.w1.usrfuncfileedit.text()
		an = self.w1.anedit.text()
		ref_a = self.w1.refaedit.text()
		sym = self.w1.symedit.text()
		npad = self.w1.npadedit.text()
				
		fourvar = self.w2.fourvarchkbx.checkState()
		rantest = self.w2.rantestchkbx.checkState()
		nch = self.w2.nchchkbx.checkState()
		debug = self.w2.debugchkbx.checkState()
		deltapsi = self.w2.deltapsiedit.text()
		startpsi = self.w2.startpsiedit.text()
		stoprnct = self.w2.stoprnctedit.text()
		chunk = self.w2.chunkedit.text()
				
		cmd1 = cmd1+" --center=" +str(ctr)+" --rs="+str(ringstep)+ " --ir=" + str(inrad)+ " --snr=" + str(snr)+ " --an='" + str(an)+ "'"+" --sym=" + str(sym)+ " --ref_a=" + str(ref_a)+ " --npad=" + str(npad)+ " --deltapsi='" + str(deltapsi)+"'"+ " --startpsi='" + str(startpsi)+ "'"+" --stoprnct=" + str(stoprnct)+ " --chunk=" + str(chunk)
		if CTF == Qt.Checked:
			cmd1 = cmd1 + " --CTF"
		if fourvar == Qt.Checked:
			cmd1 = cmd1 + " --Fourvar"		
		if nch == Qt.Checked:
			cmd1 = cmd1 + " --n"
		if rantest == Qt.Checked:
			cmd1 = cmd1 + " --rantest"
		if debug == Qt.Checked:
			cmd1 = cmd1 + " --debug"
								
		if len(userfile) < 1:
			cmd1 = cmd1 + " --function="+str(userf)
		else:
			userfile = str(userfile)
			# break it down into file name and directory path
			rind = userfile.rfind('/')
			
			if rind == -1:
				userfile = os.path.abspath(userfile)
				rind = userfile.rfind('/')
					
			fname = userfile[rind+1:]
			fname, ext = os.path.splitext(fname)
			fdir = userfile[0:rind]
			cmd1 = cmd1 + " --function=\"[" +fdir+","+fname+","+str(userf)+"]\""
				
		np = self.nprocedit.text()
		
		self.savedparmsdict = {'stackname':str(stack),'refname':str(ref),'foldername':str(output),'partradius':str(ou),'xyrange':str(xr),'trans':str(ts),'delta':str(delta),'nriter':str(maxit),'nproc':str(np),'maskname':str(mask),'center':str(ctr),"ringstep":str(ringstep),"innerradius":str(inrad),"ctf":CTF,"snr":str(snr),"fourvar":fourvar,"usrfunc":str(userf), "usrfuncfile":str(userfile),"an":str(an),"ref_a":str(ref_a),"sym":str(sym),"npad":str(npad),"deltapsi":str(deltapsi),"startpsi":str(startpsi),"stoprnct":str(stoprnct),"debug":debug,"nch":nch,"chunk":str(chunk),"rantest":rantest}

		self.w1.savedparmsdict=self.savedparmsdict
		self.w2.savedparmsdict=self.savedparmsdict
		
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
		
	def runsxali3d(self):
		self.gencmdline_ali3d(writefile=False)
		outfolder=self.savedparmsdict['foldername']
		if os.path.exists(outfolder):
			print "output folder "+outfolder+" already exists!"
			return
		process = subprocess.Popen(self.cmd,shell=True)
		self.emit(QtCore.SIGNAL("process_started"),process.pid)
		
	def saveparms(self):		
		(fname,stat)= QInputDialog.getText(self,"Save Input Parameters","Enter name of file to save parameters in",QLineEdit.Normal,"")
		if stat:
			import pickle
			output=open(fname,'wb')
			self.gencmdline_ali3d(writefile=False)
			pickle.dump(self.savedparmsdict,output)
			output.close()
		
	def repoparms_ali3d(self):		
		# repopulate with saved parms
		(fname,stat)= QInputDialog.getText(self,"Retrieve saved parameters","Enter name of file parameters were saved in",QLineEdit.Normal,"")
		if stat:
			import pickle
			pkl = open(fname,'rb')
			self.savedparmsdict = pickle.load(pkl)
			print self.savedparmsdict
			self.partradiusedit.setText(self.savedparmsdict['partradius'])
			self.stacknameedit.setText(self.savedparmsdict['stackname'])
			self.refnameedit.setText(self.savedparmsdict['refname'])		
			self.foldernameedit.setText(self.savedparmsdict['foldername'])		
			self.xyrangeedit.setText(self.savedparmsdict['xyrange'])
			self.transedit.setText(self.savedparmsdict['trans'])
			self.deltaedit.setText(self.savedparmsdict['delta'])
			self.nriteredit.setText(self.savedparmsdict['nriter'])
			self.nprocedit.setText(self.savedparmsdict['nproc'])
	
			self.w1.masknameedit.setText(self.savedparmsdict['maskname'])
			self.w1.centeredit.setText(self.savedparmsdict['center'])
			self.w1.ringstepedit.setText(self.savedparmsdict['ringstep'])
			self.w1.innerradiusedit.setText(self.savedparmsdict['innerradius'])
			self.w1.ctfchkbx.setCheckState(self.savedparmsdict['ctf'])
			self.w1.snredit.setText(self.savedparmsdict['snr'])
			self.w1.usrfuncedit.setText(self.savedparmsdict['usrfunc'])
			self.w1.usrfuncfileedit.setText(self.savedparmsdict['usrfuncfile'])
			self.w1.anedit.setText(self.savedparmsdict['an'])
			self.w1.refaedit.setText(self.savedparmsdict['ref_a'])
			self.w1.symedit.setText(self.savedparmsdict['sym'])
			self.w1.npadedit.setText(self.savedparmsdict['npad'])
			
			self.w2.fourvarchkbx.setCheckState(self.savedparmsdict['fourvar'])
			self.w2.debugchkbx.setCheckState(self.savedparmsdict['debug'])
			self.w2.nchchkbx.setCheckState(self.savedparmsdict['nch'])
			self.w2.rantestchkbx.setCheckState(self.savedparmsdict['rantest'])
			self.w2.deltapsiedit.setText(self.savedparmsdict['deltapsi'])
			self.w2.startpsiedit.setText(self.savedparmsdict['startpsi'])
			self.w2.stoprnctedit.setText(self.savedparmsdict['stoprnct'])
			self.w2.chunkedit.setText(self.savedparmsdict['chunk'])
				
	def setactiveheader(self):
		stack = self.stacknameedit.text()
		print "stack defined="+ stack
		# horatio active_refactoring Jy51i1EwmLD4tWZ9_00000_1
		# header(str(stack), "active", one=True)


	def setprojheader(self):
		#opens a file browser, showing files only in .hdf format
		ok=False
		zerostr='set xform.projection to zero'
		randstr='randomize xform.projection'
		importstr='import parameters from file'
		(item,stat)= QInputDialog.getItem(self,"xform.projection","choose option",[zerostr,randstr,importstr])
		#we convert this Qstring to a string and send it to line edit classed stackname edit of the Poptwodali window
		#self.stacknameedit.setText(str(a))
		choice= str(item)
		stack = self.stacknameedit.text()
		if stat:
			if choice == zerostr:
				header(str(stack),'xform.projection',zero=True)
			if choice == randstr:
				header(str(stack),'xform.projection',rand_alpha=True)
			if choice == importstr:
				file_name = QtGui.QFileDialog.getOpenFileName(self, "Open Data File", "", "(*)")
				a=str(QtCore.QString(file_name))
				if len(a)>0:
					header(str(stack),'xform.projection',fimport=a)
		
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

	def choose_reffile(self):
		#opens a file browser, showing files only in .hdf format
		file_name = QtGui.QFileDialog.getOpenFileName(self, "Open reference structure file", "", "(*)")
		#after the user selected a file, we obtain this filename as a Qstring
		a=QtCore.QString(file_name)
		print a
		#we convert this Qstring to a string and send it to line edit classed stackname edit of the Poptwodali window
		self.refnameedit.setText(str(a))

class Popupadvparams_ali3d_1(QWidget):
	def __init__(self,savedparms):
		QWidget.__init__(self)
		
		########################################################################################
		# layout parameters
		
		self.y1=10
		self.yspc = 4
		
		self.x1 = 10
		self.x2 = self.x1+220 #140
		self.x3 = self.x2+145 #285
		########################################################################################
		
		#Here we just set the window title
		#self.setWindowTitle('sxali3d advanced parameter selection')
		#Here we just set a label and its position in the window
		title1=QtGui.QLabel('<b>sxali3d</b> - set advanced parameters related to CTF and search', self)
		title1.move(self.x1,self.y1)
		
		self.y1 += 30
		#Labels and Line Edits for User Input
		#Just a label
		#title2= QtGui.QLabel('<b>Advanced</b> parameters', self)
		#title2.move(self.x1,self.y1)
		
		#self.y1 += 30
		
		self.savedparmsdict=savedparms
		#Example for User input stack name
		#First create the label and define its position
		maskname= QtGui.QLabel('Mask', self)
		maskname.move(self.x1,self.y1)
		#Now add a line edit and define its position
		self.masknameedit=QtGui.QLineEdit(self)
		self.masknameedit.move(self.x2,self.y1)
		#Adds a default value for the line edit
		self.masknameedit.setText(self.savedparmsdict['maskname'])
		self.masknameedit.setToolTip("filename of the file containing 3D mask. If not provided, a 3D spherical mask will be created with radius equal to outer_radius. This mask is used to multiply the reference volume for calculation of reference projections.")
		
		self.mskfile_button = QtGui.QPushButton("Open Mask", self)
		self.mskfile_button.move(self.x3, self.y1-self.yspc)
		#Here we define, that when this button is clicked, it starts subfunction choose_file
		QtCore.QObject.connect(self.mskfile_button, QtCore.SIGNAL("clicked()"), self.choose_mskfile)
		
		self.y1 += 30
		
		center= QtGui.QLabel('Center type', self)
		center.move(self.x1,self.y1)
		self.centeredit=QtGui.QLineEdit(self)
		self.centeredit.move(self.x2,self.y1)
		self.centeredit.setText(self.savedparmsdict['center'])
		self.centeredit.setToolTip('-1: average shift method; 0: no centering; 1: center of gravity (default=-1)')
		
		self.y1 += 30
		
		ringstep= QtGui.QLabel('Ring step', self)
		ringstep.move(self.x1,self.y1)
		self.ringstepedit=QtGui.QLineEdit(self)
		self.ringstepedit.move(self.x2,self.y1)
		self.ringstepedit.setText(self.savedparmsdict['ringstep'])
		self.ringstepedit.setToolTip('step between rings in rotational correlation >0  (set to 1)')
		
		self.y1 += 30
		
		innerradius= QtGui.QLabel('Inner radius', self)
		innerradius.move(self.x1,self.y1)
		self.innerradiusedit=QtGui.QLineEdit(self)
		self.innerradiusedit.move(self.x2,self.y1)
		self.innerradiusedit.setText(self.savedparmsdict['innerradius'])
		self.innerradiusedit.setToolTip('inner radius for rotational correlation > 0 (set to 1)')		
		
		self.y1 += 30
		
		ctf= QtGui.QLabel('CTF', self)
		ctf.move(self.x1,self.y1)
		self.ctfchkbx = QtGui.QCheckBox("",self)
		self.ctfchkbx.move(self.x2, self.y1)
		self.ctfchkbx.setCheckState(self.savedparmsdict['ctf'])
		self.ctfchkbx.setToolTip('Consider CTF correction during the alignment ')
		
		self.y1 += 30
		
		snr= QtGui.QLabel('SNR', self)
		snr.move(self.x1,self.y1)
		self.snredit=QtGui.QLineEdit(self)
		self.snredit.move(self.x2,self.y1)
		self.snredit.setText(self.savedparmsdict['snr'])
		self.snredit.setToolTip('Signal-to-Noise Ratio of the data')		
		
		self.y1 += 30
		
		refa= QtGui.QLabel('Method for generating reference \nprojections', self)
		refa.move(self.x1,self.y1)
		#Now add a line edit and define its position
		self.refaedit=QtGui.QLineEdit(self)
		self.refaedit.move(self.x2,self.y1)
		#Adds a default value for the line edit
		self.refaedit.setText(self.savedparmsdict['ref_a'])
		self.refaedit.setToolTip("method for generating the quasi-uniformly distributed projection directions (default S)")
		
		self.y1 += 50
		
		sym= QtGui.QLabel('Symmetry', self)
		sym.move(self.x1,self.y1)
		self.symedit=QtGui.QLineEdit(self)
		self.symedit.move(self.x2,self.y1)
		self.symedit.setText(self.savedparmsdict['sym'])
		self.symedit.setToolTip('symmetry of the refined structure')
		
		self.y1 += 30
		
		pad= QtGui.QLabel('Padding size for 3D \nreconstruction', self)
		pad.move(self.x1,self.y1)
		self.npadedit=QtGui.QLineEdit(self)
		self.npadedit.move(self.x2,self.y1)
		self.npadedit.setText(self.savedparmsdict['npad'])
		self.npadedit.setToolTip('padding size for 3D \nreconstruction')
		
		self.y1 += 50
		
		an= QtGui.QLabel('Angular neighborhood for local \nsearches', self)
		an.move(self.x1,self.y1)
		self.anedit=QtGui.QLineEdit(self)
		self.anedit.move(self.x2,self.y1)
		self.anedit.setText(self.savedparmsdict['an'])
		self.anedit.setToolTip('angular neighborhood for local searches')		

		self.y1 += 50
		
		usrfunc= QtGui.QLabel('User Function Name', self)
		usrfunc.move(self.x1,self.y1)
		self.usrfuncedit=QtGui.QLineEdit(self)
		self.usrfuncedit.move(self.x2,self.y1)
		self.usrfuncedit.setText(self.savedparmsdict['usrfunc'])
		self.usrfuncedit.setToolTip('name of the reference preparation function (ref_ali3d)')
		
		self.y1 += 30
				
		usrfuncfile= QtGui.QLabel('Enter name of external file containing user function:\n(Leave blank if file is not external to Sparx)', self)
		usrfuncfile.move(self.x1,self.y1)
		
		
		self.y1 += 40
		
		self.usrfuncfileedit=QtGui.QLineEdit(self)
		self.usrfuncfileedit.move(self.x2,self.y1)
		self.usrfuncfileedit.setText(self.savedparmsdict['usrfuncfile'])
		self.usrfuncfileedit.setToolTip('name of the external file containing user function')		
		#Function runsxali2d started when  the  RUN_button of the  Poptwodali window is clicked 
			
		#Here we define, that when this button is clicked, it starts subfunction choose_file
		self.usrfile_button = QtGui.QPushButton("Select File", self)
		self.usrfile_button.move(self.x3, self.y1-self.yspc)
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

class Popupadvparams_ali3d_2(QWidget):
	def __init__(self,savedparms):
		QWidget.__init__(self)
		
		########################################################################################
		# layout parameters
		
		self.y1=10
		self.yspc = 4
		
		self.x1 = 10
		self.x2 = self.x1+220 #140
		self.x3 = self.x2+145 # 285
		########################################################################################
		
		#Here we just set the window title
		#self.setWindowTitle('sxali3d advanced parameter selection 2')
		#Here we just set a label and its position in the window
		title1=QtGui.QLabel('<b>sxali3d</b> - set remaining advanced parameters', self)
		title1.move(self.x1,self.y1)
		
		self.y1 += 30
		
		self.savedparmsdict=savedparms
			
		fourvar= QtGui.QLabel('Compute fourier variance', self)
		fourvar.move(self.x1,self.y1)
		self.fourvarchkbx=QtGui.QCheckBox("",self)
		self.fourvarchkbx.move(self.x2,self.y1)
		self.fourvarchkbx.setCheckState(self.savedparmsdict['fourvar'])
		self.fourvarchkbx.setToolTip('compute Fourier variance')
		
		self.y1 += 30
		
		deltapsi= QtGui.QLabel('Delta psi for coarse search', self)
		deltapsi.move(self.x1,self.y1)
		self.deltapsiedit=QtGui.QLineEdit(self)
		self.deltapsiedit.move(self.x2,self.y1)
		self.deltapsiedit.setText(self.savedparmsdict['deltapsi'])
		self.deltapsiedit.setToolTip('Delta psi for coarse search')		

		self.y1 += 30
		
		startpsi= QtGui.QLabel('Start psi for coarse search', self)
		startpsi.move(self.x1,self.y1)
		self.startpsiedit=QtGui.QLineEdit(self)
		self.startpsiedit.move(self.x2,self.y1)
		self.startpsiedit.setText(self.savedparmsdict['startpsi'])
		self.startpsiedit.setToolTip('Start psi for coarse search')		

		self.y1 += 30
		
		stoprnct= QtGui.QLabel('Minimum percentage of particles \nwhich must change orientation, \nbelow which the program stops', self)
		stoprnct.move(self.x1,self.y1)
		self.stoprnctedit=QtGui.QLineEdit(self)
		self.stoprnctedit.move(self.x2,self.y1)
		self.stoprnctedit.setText(self.savedparmsdict['stoprnct'])
		self.stoprnctedit.setToolTip('Minimum percentage of particles that change orientation to stop the program ')		

		self.y1 += 60
		
		nch= QtGui.QLabel('New version of ali3d where \na percentage of data \nis used for alignment', self)
		nch.move(self.x1,self.y1)
		self.nchchkbx=QtGui.QCheckBox("",self)
		self.nchchkbx.move(self.x2,self.y1)
		self.nchchkbx.setCheckState(self.savedparmsdict['nch'])
		self.nchchkbx.setToolTip('new')
		
		self.y1 += 60
		
		chunk= QtGui.QLabel('Percentage of data used \nfor alignment', self)
		chunk.move(self.x1,self.y1)
		self.chunkedit=QtGui.QLineEdit(self)
		self.chunkedit.move(self.x2,self.y1)
		self.chunkedit.setText(self.savedparmsdict['chunk'])
		self.chunkedit.setToolTip('percentage of data used for alignment')		

		self.y1 += 50
		
		dbg= QtGui.QLabel('Debug mode', self)
		dbg.move(self.x1,self.y1)
		self.debugchkbx=QtGui.QCheckBox("",self)
		self.debugchkbx.move(self.x2,self.y1)
		self.debugchkbx.setCheckState(self.savedparmsdict['debug'])
		self.debugchkbx.setToolTip('debug mode')
		
		self.y1 += 30
		
		rant= QtGui.QLabel('rantest', self)
		rant.move(self.x1,self.y1)
		self.rantestchkbx=QtGui.QCheckBox("",self)
		self.rantestchkbx.move(self.x2,self.y1)
		self.rantestchkbx.setCheckState(self.savedparmsdict['rantest'])
		self.rantestchkbx.setToolTip('rantest')
		
		self.y1 += 30
		
	
#################
## kmeans

#Layout of the Pop Up window Popuptwodali (for sxali2d); started by the function twodali of the main window		
class Popupkmeans(QWidget):
	def __init__(self):
		QWidget.__init__(self)
		
		#######################################################################################
		# class variables
		self.cmd = ""
		# populate with default values
		self.savedparmsdict = {'stackname':'','foldername':'','kc':'2','trials':'1','maxiter':'100','nproc':'1','randsd':'-1','maskname':'',"ctf":Qt.Unchecked,"crit":"all","normalize":Qt.Unchecked, "init_method":"rnd"}

		#######################################################################################
		# Layout parameters
		
		self.y1 = 10 # title and Repopulate button
		self.y2 = self.y1 + 78 # Text boxes for inputting parameters
		self.y3 = self.y2 + 222 # activate images button and set xform.align2d button
		self.y4 = self.y3 + 80 # Advanced Parameters, Save Input and Generate command line buttons
		self.y5 = self.y4 + 95 # run button 
		self.yspc = 4
		
		self.x1 = 10 # first column (text box labels)
		self.x2 = self.x1 + 200 # second column (text boxes)
		self.x3 = self.x2+145 # third column (Open .hdf button)
		self.x4 = self.x3+100 # fourth column (Open .bdb button)
		self.x5 = 230 # run button
		#######################################################################################
		
		#Here we just set the window title
		self.setWindowTitle('sxk_means')
		#Here we just set a label and its position in the window
		title1=QtGui.QLabel('<b>sxk_means</b> - K-means classification of a set of images', self)
		title1.move(self.x1,self.y1)
		self.y1 += 30

		self.repopbtn = QPushButton("Retrieve saved parameters", self)
		self.repopbtn.move(self.x1-5,self.y1)
		#sets an infotip for this Pushbutton
		self.repopbtn.setToolTip('Retrieve saved parameters')
		#when this button is clicked, this action starts the subfunction twodali
		self.connect(self.repopbtn, SIGNAL("clicked()"), self.repoparms_kmeans)

		#######################################################################################
		#Here we create a Button(file_button with title run open .hdf) and its position in the window
		self.file_button = QtGui.QPushButton("Open .hdf", self)
		self.file_button.move(self.x3, self.y2-self.yspc)
		#Here we define, that when this button is clicked, it starts subfunction choose_file
		QtCore.QObject.connect(self.file_button, QtCore.SIGNAL("clicked()"), self.choose_file)
		#exactly the same as above, but for subfunction choose_file1
		self.file_button1 = QtGui.QPushButton("Open .bdb", self)
		self.file_button1.move(self.x4,self.y2-self.yspc)
		QtCore.QObject.connect(self.file_button1, QtCore.SIGNAL("clicked()"), self.choose_file1)
		
		stackname= QtGui.QLabel('Name of input stack', self)
		stackname.move(self.x1,self.y2)
		self.stacknameedit=QtGui.QLineEdit(self)
		self.stacknameedit.move(self.x2,self.y2)
		self.stacknameedit.setText(self.savedparmsdict['stackname'])
		self.y2 += 30
		
		foldername= QtGui.QLabel('Output folder', self)
		foldername.move(self.x1,self.y2)
		self.foldernameedit=QtGui.QLineEdit(self)
		self.foldernameedit.move(self.x2,self.y2)
		self.foldernameedit.setText(self.savedparmsdict['foldername'])		
		
		self.outinfobtn = QPushButton("Output Info", self)
		self.outinfobtn.move(self.x3,  self.y2-self.yspc)
		#sets an infotip for this Pushbutton
		self.outinfobtn.setToolTip('Output Info')
		#when this button is clicked, this action starts the subfunction twodali
		self.connect(self.outinfobtn, SIGNAL("clicked()"), self.outputinfo_kmeans)
		self.y2 += 30
		
		kc= QtGui.QLabel('Number of clusters', self)
		kc.move(self.x1,self.y2)
		self.kcedit=QtGui.QLineEdit(self)
		self.kcedit.move(self.x2,self.y2)
		self.kcedit.setText(self.savedparmsdict['kc'])
		self.kcedit.setToolTip('The requested number of clusters')		
		self.y2 += 30

		trials= QtGui.QLabel('Number of trials', self)
		trials.move(self.x1,self.y2)
		self.trialsedit=QtGui.QLineEdit(self)
		self.trialsedit.move(self.x2,self.y2)
		self.trialsedit.setText(self.savedparmsdict['trials'])
		self.trialsedit.setToolTip('number of trials of K-means (see description below) (default one trial). MPI version ignore --trials, the number of trials in MPI version will be the number of cpu used.')
		self.y2 += 30
		
		maxiter= QtGui.QLabel('Number of iterations', self)
		maxiter.move(self.x1,self.y2)
		self.maxiteredit=QtGui.QLineEdit(self)
		self.maxiteredit.move(self.x2,self.y2)
		self.maxiteredit.setText(self.savedparmsdict['maxiter'])
		self.maxiteredit.setToolTip('maximum number of iterations the program will perform (default 100)')
		self.y2 += 30
		
		# make ctf, normalize and init_method radio button...
		
		nproc= QtGui.QLabel('MPI processors', self)
		nproc.move(self.x1,self.y2)
		self.nprocedit=QtGui.QLineEdit(self)
		self.nprocedit.move(self.x2,self.y2)
		self.nprocedit.setText(self.savedparmsdict['nproc'])
		self.nprocedit.setToolTip('The number of processors to use. Default is single processor mode')

		##########################################################
		
		header=QtGui.QLabel('The clustering is performed only using images from the stack that has flag active set to one', self)
		header.move(self.x1,self.y3)
		self.y3 += 30
		
		#not linked to a function yet
		self.activeheader_button = QtGui.QPushButton("activate all images", self)
		self.activeheader_button.move(self.x1-5, self.y3)
		self.connect(self.activeheader_button, SIGNAL("clicked()"), self.setactiveheader)
		
		######################################################################################
		
		self.savepbtn = QPushButton("Save Input Parameters", self)
		self.savepbtn.move(self.x1-5, self.y4)
		#sets an infotip for this Pushbutton
		self.savepbtn.setToolTip('Save Input Parameters')
		#when this button is clicked, this action starts the subfunction twodali
		self.connect(self.savepbtn, SIGNAL("clicked()"), self.saveparms)
		self.y4+=30
		
		self.cmdlinebtn = QPushButton("Generate command line from input parameters", self)
		self.cmdlinebtn.move(self.x1-5, self.y4)
		#sets an infotip for this Pushbutton
		self.cmdlinebtn.setToolTip('Generate command line using input parameters')
		#when this button is clicked, this action starts the subfunction twodali
		self.connect(self.cmdlinebtn, SIGNAL("clicked()"), self.gencmdline_kmeans)

		#######################################################################
		 #Here we create a Button(Run_button with title run sxali2d) and its position in the window
		self.RUN_button = QtGui.QPushButton('Run sxk_means', self)
		# make 3D textured push button look
		s = "QPushButton {font: bold; color: #000;border: 1px solid #333;border-radius: 11px;padding: 2px;background: qradialgradient(cx: 0, cy: 0,fx: 0.5, fy:0.5,radius: 1, stop: 0 #fff, stop: 1 #8D0);min-width:90px;margin:5px} QPushButton:pressed {font: bold; color: #000;border: 1px solid #333;border-radius: 11px;padding: 2px;background: qradialgradient(cx: 0, cy: 0,fx: 0.5, fy:0.5,radius: 1, stop: 0 #fff, stop: 1 #084);min-width:90px;margin:5px}"
		
		self.RUN_button.setStyleSheet(s)
		self.RUN_button.move(self.x5, self.y5)
		#Here we define, that when this button is clicked, it starts subfunction runsxali2d
		self.connect(self.RUN_button, SIGNAL("clicked()"), self.runsxk_means)
		#Labels and Line Edits for User Input		

	 #Function runsxali2d started when  the  RUN_button of the  Poptwodali window is clicked 
	def outputinfo_kmeans(self):
		QMessageBox.information(self, "kmeans output",'The directory to which the averages of K clusters, and the variance. The classification charts are written to the logfile. Warning: If the output directory already exists, the program will crash and an error message will come up. Please change the name of directory and restart the program.')
		
	def gencmdline_kmeans(self,writefile=True):
		#Here we just read in all user inputs in the line edits of the Poptwodali window
		stack = self.stacknameedit.text()
		print "stack defined="+ stack
		output=self.foldernameedit.text()
		print "output folder="+ output
		kc=self.kcedit.text()
		print "Number of requested clusters="+ kc
		trials=self.trialsedit.text()
		print "Number of trials="+ trials
		maxiter=self.maxiteredit.text()
		print "maximum number of iterations=" +maxiter
		
		cmd1 = "sxk_means.py "+str(stack) +" "+ str(output)
		
		args = " --K="+ str(kc)+ " --trials="+str(trials)+ " --maxit="+str(maxiter)
		
		
		mask = self.w1.masknameedit.text()
				
		if len(str(mask))> 1:
				cmd1 = cmd1+" "+str(mask) 
						
		cmd1 = cmd1 + args
		
		randsd=self.w1.randsdedit.text()
		ctf=self.w1.ctfchkbx.checkState()
		normalize=self.w1.normachkbx.checkState()
		init_method=self.w1.init_methodedit.text()
		crit=self.w1.critedit.text()
				
		cmd1 = cmd1+" --rand_seed=" +str(randsd)+" --init_method="+str(init_method)+" --crit="+str(crit)
				
		if ctf == Qt.Checked:
				cmd1 = cmd1 + " --CTF"
		if normalize == Qt.Checked:
				cmd1 = cmd1 + " --normalize"
						
		np = self.nprocedit.text()
		
		self.savedparmsdict = {'stackname':str(stack),'foldername':str(output),'kc':str(kc),'trials':str(trials),'maxiter':str(maxiter),'nproc':str(np),'randsd':str(randsd),'maskname':str(mask),"ctf":ctf,"crit":str(crit),"normalize":normalize, "init_method":str(init_method)}

		self.w1.savedparmsdict=self.savedparmsdict
		
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
		
	def runsxk_means(self):
		self.gencmdline_kmeans(writefile=False)
		outfolder=self.savedparmsdict['foldername']
		if os.path.exists(outfolder):
			print "output folder "+outfolder+" already exists!"
			return
		process = subprocess.Popen(self.cmd,shell=True)
		self.emit(QtCore.SIGNAL("process_started"),process.pid)
		
	def saveparms(self):		
		(fname,stat)= QInputDialog.getText(self,"Save Input Parameters","Enter name of file to save parameters in",QLineEdit.Normal,"")
		if stat:
			import pickle
			output=open(fname,'wb')
			self.gencmdline_kmeans(writefile=False)
			pickle.dump(self.savedparmsdict,output)
			output.close()
		
	def repoparms_kmeans(self):		
		(fname,stat)= QInputDialog.getText(self,"Retrieve saved parameters","Enter name of file parameters were saved in",QLineEdit.Normal,"")
		if stat:
			import pickle
			pkl = open(fname,'rb')
			self.savedparmsdict = pickle.load(pkl)
			self.kcedit.setText(self.savedparmsdict['kc'])
			self.stacknameedit.setText(self.savedparmsdict['stackname'])		
			self.foldernameedit.setText(self.savedparmsdict['foldername'])		
			self.trialsedit.setText(self.savedparmsdict['trials'])
			self.maxiteredit.setText(self.savedparmsdict['maxiter'])
			self.nprocedit.setText(self.savedparmsdict['nproc'])
			self.w1.masknameedit.setText(self.savedparmsdict['maskname'])
			self.w1.randsdedit.setText(self.savedparmsdict['randsd'])
			self.w1.ctfchkbx.setCheckState(self.savedparmsdict['ctf'])
			self.w1.normachkbx.setCheckState(self.savedparmsdict['normalize'])
			self.w1.critedit.setText(self.savedparmsdict['crit'])
			self.w1.init_methodedit.setText(self.savedparmsdict['init_method'])
				
	def setactiveheader(self):
		stack = self.stacknameedit.text()
		print "stack defined="+ stack
		# horatio active_refactoring Jy51i1EwmLD4tWZ9_00000_1
		# header(str(stack), "active", one=True)



	def choose_file(self):
		#opens a file browser, showing files only in .hdf format
		file_name = QtGui.QFileDialog.getOpenFileName(self, "Open Data File", "", "HDF files (*.hdf)")
		#after the user selected a file, we obtain this filename as a Qstring
		a=QtCore.QString(file_name)
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
		
class Popupadvparams_kmeans(QWidget):
	def __init__(self,savedparms):
		QWidget.__init__(self)
		
		self.x1 = 10
		self.x2 = 140
		self.x3 = 285
		
		self.y1 = 10
		self.yspc = 4
		
		#Here we just set the window title
		self.setWindowTitle('sxk_means advanced parameter selection')
		#Here we just set a label and its position in the window
		title1=QtGui.QLabel('<b>sxk_means</b> - set advanced params', self)
		title1.move(self.x1,self.y1)
		#Labels and Line Edits for User Input
		#Just a label
		self.y1 += 30
		
		title2= QtGui.QLabel('<b>Advanced</b> parameters', self)
		title2.move(self.x1,self.y1)
		
		self.y1 += 30
		
		self.savedparmsdict=savedparms
		#Example for User input stack name
		#First create the label and define its position
		maskname= QtGui.QLabel('Mask', self)
		maskname.move(self.x1,self.y1)
		#Now add a line edit and define its position
		self.masknameedit=QtGui.QLineEdit(self)
		self.masknameedit.move(self.x2,self.y1)
		#Adds a default value for the line edit
		self.masknameedit.setText(self.savedparmsdict['maskname'])
		
		self.mskfile_button = QtGui.QPushButton("Open .hdf", self)
		self.mskfile_button.move(self.x3, self.y1-self.yspc)
		#Here we define, that when this button is clicked, it starts subfunction choose_file
		QtCore.QObject.connect(self.mskfile_button, QtCore.SIGNAL("clicked()"), self.choose_mskfile)
		
		self.y1 += 30
		
		randsd= QtGui.QLabel('Random Seed', self)
		randsd.move(self.x1,self.y1)
		self.randsdedit=QtGui.QLineEdit(self)
		self.randsdedit.move(self.x2,self.y1)
		self.randsdedit.setText(self.savedparmsdict['randsd'])
		self.randsdedit.setToolTip('the seed used to generating random numbers (set to -1, means different and pseudo-random each time)')
		
		self.y1 += 30
		
		ctf= QtGui.QLabel('CTF', self)
		ctf.move(self.x1,self.y1)
		self.ctfchkbx = QtGui.QCheckBox("",self)
		self.ctfchkbx.move(self.x2, self.y1)
		self.ctfchkbx.setCheckState(self.savedparmsdict['ctf'])
		
		self.y1 += 30
		
		norma= QtGui.QLabel('normalize', self)
		norma.move(self.x1,self.y1)
		self.normachkbx = QtGui.QCheckBox("",self)
		self.normachkbx.move(self.x2, self.y1)
		self.normachkbx.setCheckState(self.savedparmsdict['normalize'])
		self.y1 += 30
		
		crit= QtGui.QLabel('Criterion', self)
		crit.move(self.x1,self.y1)
		self.critedit=QtGui.QLineEdit(self)
		self.critedit.move(self.x2,self.y1)
		self.critedit.setText(self.savedparmsdict['crit'])
		self.critedit.setToolTip('names of criterion used: \'all\' all criterions, \'C\' Coleman, \'H\' Harabasz or \'D\' Davies-Bouldin, thoses criterions return the values of classification quality, see also sxk_means_groups. Any combination is accepted, i.e., \'CD\', \'HC\', \'CHD\'.')
		
		self.y1 += 30
		
		initmethod= QtGui.QLabel('Criterion', self)
		initmethod.move(self.x1,self.y1)
		self.init_methodedit=QtGui.QLineEdit(self)
		self.init_methodedit.move(self.x2,self.y1)
		self.init_methodedit.setText(self.savedparmsdict['init_method'])
		self.init_methodedit.setToolTip('Method used to initialize partition: "rnd" randomize or "d2w" for d2 weighting initialization (default is rnd)')
		
		
	def choose_mskfile(self):
		#opens a file browser, showing files only in .hdf format
		file_name = QtGui.QFileDialog.getOpenFileName(self, "Open File Containing Mask", "", "(*)")
		#after the user selected a file, we obtain this filename as a Qstring
		a=QtCore.QString(file_name)
		print a
		#we convert this Qstring to a string and send it to line edit classed stackname edit of the Poptwodali window
		self.masknameedit.setText(str(a))

		
class Popupkmeansgroups(QWidget):
	def __init__(self):
		QWidget.__init__(self)
		
		#######################################################################################
		# class variables
		self.cmd = ""
		# populate with default values
		self.savedparmsdict = {'stackname':'','foldername':'','kc1':'2','kc2':'3','trials':'1','nproc':'1','randsd':'-1','maskname':'',"ctf":Qt.Unchecked,'maxiter':'100'}
		
		#######################################################################################
		# Layout parameters
		
		self.y1 = 10 # title and Repopulate button
		self.y2 = self.y1 + 78 # Text boxes for inputting parameters
		self.y3 = self.y2 + 222 # activate images button and set xform.align2d button
		self.y4 = self.y3 + 80 # Advanced Parameters, Save Input and Generate command line buttons
		self.y5 = self.y4 + 95 # run button 
		self.yspc = 4
		
		self.x1 = 10 # first column (text box labels)
		self.x2 = self.x1 + 200 # second column (text boxes)
		self.x3 = self.x2+145 # third column (Open .hdf button)
		self.x4 = self.x3+100 # fourth column (Open .bdb button)
		self.x5 = 230 # run button
		#######################################################################################
		
		#Here we just set the window title
		self.setWindowTitle('sxk_means_groups')
		#Here we just set a label and its position in the window
		title1=QtGui.QLabel('<b>sxk_means_groups</b> - K-means groups classification of a set of images', self)
		title1.move(self.x1,self.y1)
		self.y1 += 30

		self.repopbtn = QPushButton("Retrieve saved parameters", self)
		self.repopbtn.move(self.x1-5,self.y1)
		#sets an infotip for this Pushbutton
		self.repopbtn.setToolTip('Retrieve saved parameters')
		#when this button is clicked, this action starts the subfunction twodali
		self.connect(self.repopbtn, SIGNAL("clicked()"), self.repoparms_kmeans_groups)

		#######################################################################################
		#Here we create a Button(file_button with title run open .hdf) and its position in the window
		self.file_button = QtGui.QPushButton("Open .hdf", self)
		self.file_button.move(self.x3, self.y2-self.yspc)
		#Here we define, that when this button is clicked, it starts subfunction choose_file
		QtCore.QObject.connect(self.file_button, QtCore.SIGNAL("clicked()"), self.choose_file)
		#exactly the same as above, but for subfunction choose_file1
		self.file_button1 = QtGui.QPushButton("Open .bdb", self)
		self.file_button1.move(self.x4,self.y2-self.yspc)
		QtCore.QObject.connect(self.file_button1, QtCore.SIGNAL("clicked()"), self.choose_file1)
		
		stackname= QtGui.QLabel('Name of input stack', self)
		stackname.move(self.x1,self.y2)
		self.stacknameedit=QtGui.QLineEdit(self)
		self.stacknameedit.move(self.x2,self.y2)
		self.stacknameedit.setText(self.savedparmsdict['stackname'])
		self.y2 += 30
		
		foldername= QtGui.QLabel('Output folder', self)
		foldername.move(self.x1,self.y2)
		self.foldernameedit=QtGui.QLineEdit(self)
		self.foldernameedit.move(self.x2,self.y2)
		self.foldernameedit.setText(self.savedparmsdict['foldername'])		
		
		self.outinfobtn = QPushButton("Output Info", self)
		self.outinfobtn.move(self.x3,  self.y2-self.yspc)
		#sets an infotip for this Pushbutton
		self.outinfobtn.setToolTip('Output Info')
		#when this button is clicked, this action starts the subfunction twodali
		self.connect(self.outinfobtn, SIGNAL("clicked()"), self.outputinfo_kmeans_groups)
		self.y2 += 30
		
		kc1= QtGui.QLabel('Minimum number of clusters', self)
		kc1.move(self.x1,self.y2)
		self.kc1edit=QtGui.QLineEdit(self)
		self.kc1edit.move(self.x2,self.y2)
		self.kc1edit.setText(self.savedparmsdict['kc1'])
		self.kc1edit.setToolTip('Minimum requested number of clusters')		
		self.y2 += 30
		
		kc2= QtGui.QLabel('Maximum number of clusters', self)
		kc2.move(self.x1,self.y2)
		self.kc2edit=QtGui.QLineEdit(self)
		self.kc2edit.move(self.x2,self.y2)
		self.kc2edit.setText(self.savedparmsdict['kc2'])
		self.kc2edit.setToolTip('The requested number of clusters')		
		self.y2 += 30

		trials= QtGui.QLabel('Number of trials', self)
		trials.move(self.x1,self.y2)
		self.trialsedit=QtGui.QLineEdit(self)
		self.trialsedit.move(self.x2,self.y2)
		self.trialsedit.setText(self.savedparmsdict['trials'])
		self.trialsedit.setToolTip('number of trials of K-means (see description below) (default one trial). MPI version ignore --trials, the number of trials in MPI version will be the number of cpu used.')
		self.y2 += 30
		
		maxiter= QtGui.QLabel('Number of iterations', self)
		maxiter.move(self.x1,self.y2)
		self.maxiteredit=QtGui.QLineEdit(self)
		self.maxiteredit.move(self.x2,self.y2)
		self.maxiteredit.setText(self.savedparmsdict['maxiter'])
		self.maxiteredit.setToolTip('maximum number of iterations the program will perform (default 100)')
		self.y2 += 30
		# make ctf, normalize and init_method radio button...
		
		nproc= QtGui.QLabel('MPI processors', self)
		nproc.move(self.x1,self.y2)
		self.nprocedit=QtGui.QLineEdit(self)
		self.nprocedit.move(self.x2,self.y2)
		self.nprocedit.setText(self.savedparmsdict['nproc'])
		self.nprocedit.setToolTip('The number of processors to use. Default is single processor mode')

		##########################################################
		
		header=QtGui.QLabel('The clustering is performed only using images from the stack that has flag active set to one', self)
		header.move(self.x1,self.y3)
		self.y3 += 30
		
		#not linked to a function yet
		self.activeheader_button = QtGui.QPushButton("activate all images", self)
		self.activeheader_button.move(self.x1-5, self.y3)
		self.connect(self.activeheader_button, SIGNAL("clicked()"), self.setactiveheader)
		
		######################################################################################
		
		self.savepbtn = QPushButton("Save Input Parameters", self)
		self.savepbtn.move(self.x1-5, self.y4)
		#sets an infotip for this Pushbutton
		self.savepbtn.setToolTip('Save Input Parameters')
		#when this button is clicked, this action starts the subfunction twodali
		self.connect(self.savepbtn, SIGNAL("clicked()"), self.saveparms)
		self.y4+=30
		
		self.cmdlinebtn = QPushButton("Generate command line from input parameters", self)
		self.cmdlinebtn.move(self.x1-5, self.y4)
		#sets an infotip for this Pushbutton
		self.cmdlinebtn.setToolTip('Generate command line using input parameters')
		#when this button is clicked, this action starts the subfunction twodali
		self.connect(self.cmdlinebtn, SIGNAL("clicked()"), self.gencmdline_kmeans_groups)

		#######################################################################
		 #Here we create a Button(Run_button with title run sxali2d) and its position in the window
		self.RUN_button = QtGui.QPushButton('Run sxk_means', self)
		# make 3D textured push button look
		s = "QPushButton {font: bold; color: #000;border: 1px solid #333;border-radius: 11px;padding: 2px;background: qradialgradient(cx: 0, cy: 0,fx: 0.5, fy:0.5,radius: 1, stop: 0 #fff, stop: 1 #8D0);min-width:90px;margin:5px} QPushButton:pressed {font: bold; color: #000;border: 1px solid #333;border-radius: 11px;padding: 2px;background: qradialgradient(cx: 0, cy: 0,fx: 0.5, fy:0.5,radius: 1, stop: 0 #fff, stop: 1 #084);min-width:90px;margin:5px}"
		
		self.RUN_button.setStyleSheet(s)
		self.RUN_button.move(self.x5, self.y5)
		#Here we define, that when this button is clicked, it starts subfunction runsxali2d
		self.connect(self.RUN_button, SIGNAL("clicked()"), self.runsxk_means_groups)
		#Labels and Line Edits for User Input		

				
	 #Function runsxali2d started when  the  RUN_button of the  Poptwodali window is clicked 
	def outputinfo_kmeans_groups(self):
		QMessageBox.information(self, "k_means_groups output",'Output directory will contain two text files: output_file and output_file.p, where output_file is the name of the output folder. \n\n output_file will contain columns according the criteria chosen, for example if crit=\'CHD\', the columns of numbers: (1) number of clusters, (2) values of Coleman criterion, (3) values of Harabasz criterion and (4) values of Davies-Bouldin criterion \n\n output_file.p will contain a gnuplot script, this file allow plot directly the values of all criteria with the same range. Use this command in gnuplot: load \'output_file.p\'\n\n WATCH_GRP_KMEANS or WATCH_MPI_GRP_KMEANS contain the progress of k-means groups. This file can be read in real-time to watch the evolution of criteria')
		
	def gencmdline_kmeans_groups(self,writefile=True):
		#Here we just read in all user inputs in the line edits of the Poptwodali window
		stack = self.stacknameedit.text()
		print "stack defined="+ stack
		output=self.foldernameedit.text()
		print "output folder="+ output
		kc1=self.kc1edit.text()
		print "Minimum number of requested clusters="+ kc1
		kc2=self.kc2edit.text()
		print "Maximum number of requested clusters="+ kc2
		trials=self.trialsedit.text()
		print "Number of trials="+ trials
		maxiter=self.maxiteredit.text()
		print "maxiter="+ maxiter
		
		cmd1 = "sxk_means_groups.py "+str(stack) +" "+ str(output)
		
		args = " --K1="+ str(kc1)+" --K2="+ str(kc2)+ " --trials="+str(trials)+ " --maxit="+str(maxiter)
		
		mask = self.w1.masknameedit.text()
				
		if len(str(mask))> 1:
			cmd1 = cmd1+" "+str(mask) 
						
		cmd1 = cmd1 + args
				
		randsd=self.w1.randsdedit.text()
		ctf=self.w1.ctfchkbx.checkState()
				
		cmd1 = cmd1+" --rand_seed=" +str(randsd)
				
		if ctf == Qt.Checked:
			cmd1 = cmd1 + " --CTF"
						
		np = self.nprocedit.text()
		
		self.savedparmsdict = {'stackname':str(stack),'foldername':str(output),'kc1':str(kc1),'kc2':str(kc2),'trials':str(trials),'nproc':str(np),'randsd':str(randsd),'maskname':str(mask),"ctf":ctf,'maxiter':str(maxiter)}

		self.w1.savedparmsdict=self.savedparmsdict
		
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
		
	def runsxk_means_groups(self):
		self.gencmdline_kmeans_groups(writefile=False)
		outfolder=self.savedparmsdict['foldername']
		if os.path.exists(outfolder):
			print "output folder "+outfolder+" already exists!"
			return
		process = subprocess.Popen(self.cmd,shell=True)
		self.emit(QtCore.SIGNAL("process_started"),process.pid)
		
	def saveparms(self):		
		# save all the parms in a text file so we can repopulate if user requests
		(fname,stat)= QInputDialog.getText(self,"Save Input Parameters","Enter name of file to save parameters in",QLineEdit.Normal,"")
		if stat:
			import pickle
			output=open(fname,'wb')
			self.gencmdline_kmeans_groups(writefile=False)
			pickle.dump(self.savedparmsdict,output)
			output.close()
		
	def repoparms_kmeans_groups(self):		
		(fname,stat)= QInputDialog.getText(self,"Get Input Parameters","Enter name of file parameters were saved in",QLineEdit.Normal,"")
		if stat:
			import pickle
			pkl = open(fname,'rb')
			self.savedparmsdict = pickle.load(pkl)
			self.kc1edit.setText(self.savedparmsdict['kc1'])
			self.kc2edit.setText(self.savedparmsdict['kc2'])
			self.stacknameedit.setText(self.savedparmsdict['stackname'])		
			self.foldernameedit.setText(self.savedparmsdict['foldername'])		
			self.trialsedit.setText(self.savedparmsdict['trials'])
			self.maxiteredit.setText(self.savedparmsdict['maxiter'])
			self.nprocedit.setText(self.savedparmsdict['nproc'])
			self.w1.masknameedit.setText(self.savedparmsdict['maskname'])
			self.w1.randsdedit.setText(self.savedparmsdict['randsd'])
			self.w1.ctfchkbx.setCheckState(self.savedparmsdict['ctf'])
				
	def setactiveheader(self):
		stack = self.stacknameedit.text()
		print "stack defined="+ stack
		# horatio active_refactoring Jy51i1EwmLD4tWZ9_00000_1
		# header(str(stack), "active", one=True)



	def choose_file(self):
		#opens a file browser, showing files only in .hdf format
		file_name = QtGui.QFileDialog.getOpenFileName(self, "Open Data File", "", "HDF files (*.hdf)")
		#after the user selected a file, we obtain this filename as a Qstring
		a=QtCore.QString(file_name)
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

class Popupadvparams_kmeans_groups(QWidget):
	def __init__(self,savedparms):
		QWidget.__init__(self)
		
		self.x1 = 10
		self.x2 = 140
		self.x3 = 285
		
		self.y1 = 10
		self.yspc = 4
		
		#Here we just set the window title
		self.setWindowTitle('sxk_means_groups advanced parameter selection')
		#Here we just set a label and its position in the window
		title1=QtGui.QLabel('<b>sxk_means_groups</b> - set advanced params', self)
		title1.move(self.x1,self.y1)
		#Labels and Line Edits for User Input
		#Just a label
		self.y1 += 30
		
		title2= QtGui.QLabel('<b>Advanced</b> parameters', self)
		title2.move(self.x1,self.y1)
		
		self.y1 += 30
		
		self.savedparmsdict=savedparms
		#Example for User input stack name
		#First create the label and define its position
		maskname= QtGui.QLabel('Mask', self)
		maskname.move(self.x1,self.y1)
		#Now add a line edit and define its position
		self.masknameedit=QtGui.QLineEdit(self)
		self.masknameedit.move(self.x2,self.y1)
		#Adds a default value for the line edit
		self.masknameedit.setText(self.savedparmsdict['maskname'])
		
		self.mskfile_button = QtGui.QPushButton("Open .hdf", self)
		self.mskfile_button.move(self.x3, self.y1-self.yspc)
		#Here we define, that when this button is clicked, it starts subfunction choose_file
		QtCore.QObject.connect(self.mskfile_button, QtCore.SIGNAL("clicked()"), self.choose_mskfile)
		
		self.y1 += 30
		
		randsd= QtGui.QLabel('Random Seed', self)
		randsd.move(self.x1,self.y1)
		self.randsdedit=QtGui.QLineEdit(self)
		self.randsdedit.move(self.x2,self.y1)
		self.randsdedit.setText(self.savedparmsdict['randsd'])
		self.randsdedit.setToolTip('the seed used to generating random numbers (set to -1, means different and pseudo-random each time)')
		
		self.y1 += 30
		
		ctf= QtGui.QLabel('CTF', self)
		ctf.move(self.x1,self.y1)
		self.ctfchkbx = QtGui.QCheckBox("",self)
		self.ctfchkbx.move(self.x2, self.y1)
		self.ctfchkbx.setCheckState(self.savedparmsdict['ctf'])
		
		self.y1 += 30
		
	def choose_mskfile(self):
		#opens a file browser, showing files only in .hdf format
		file_name = QtGui.QFileDialog.getOpenFileName(self, "Open File Containing Mask", "", "(*)")
		#after the user selected a file, we obtain this filename as a Qstring
		a=QtCore.QString(file_name)
		print a
		#we convert this Qstring to a string and send it to line edit classed stackname edit of the Poptwodali window
		self.masknameedit.setText(str(a))


class Popuppdb2em(QWidget):
	def __init__(self):
		QWidget.__init__(self)
		
		#######################################################################################
		# class variables
		self.cmd = ""
		# populate with default values
		self.savedparmsdict = {'pdbfile':'','output':'','apix':'1.0','box':'150','het':Qt.Unchecked,'center':'n','Och':Qt.Unchecked,'quiet':Qt.Unchecked,'tr0':''}
		if (options.demo == DEMO_mpibdbctf) or (options.demo == DEMO_mpibdb):
			self.savedparmsdict['pdbfile'] = '../tteftu_with_tRNA.pdb'
			self.savedparmsdict['output']='tmp.hdf'
			self.savedparmsdict['apix']='5.2'
			self.savedparmsdict['box']='64'
		#######################################################################################
		# Layout parameters
		
		self.y1 = 10 # title and Repopulate button
		self.y2 = self.y1 + 85 # Text boxes for inputting parameters
		#self.y3 = self.y2 + 222 # activate images button and set xform.align2d button
		self.y4 = self.y2 + 290 # Advanced Parameters, Save Input and Generate command line buttons
		self.y5 = self.y4 + 95 # run button 
		self.yspc = 4
		
		self.x1 = 10 # first column (text box labels)
		self.x2 = self.x1 + 260 # second column (text boxes)
		self.x3 = self.x2+145 # third column (Open .hdf button)
		self.x4 = self.x3+100 # fourth column (Open .bdb button)
		self.x5 = 230 # run button
		#######################################################################################

		#Here we just set the window title
		self.setWindowTitle('sxpdb2em')
		#Here we just set a label and its position in the window
		title1=QtGui.QLabel('<b>sxpdb2em</b> - convert atomic model (pdb file) into sampled electron density map', self)
		title1.move(self.x1,self.y1)
		self.y1 += 30

		self.repopbtn = QPushButton("Retrieve saved parameters", self)
		self.repopbtn.move(self.x1-5,self.y1)
		#sets an infotip for this Pushbutton
		self.repopbtn.setToolTip('Retrieve saved parameters')
		#when this button is clicked, this action starts the subfunction twodali
		self.connect(self.repopbtn, SIGNAL("clicked()"), self.repoparms_pdb2em)

		#######################################################################################
		#Here we create a Button(file_button with title run open .hdf) and its position in the window
		self.file_button = QtGui.QPushButton("Open .pdb", self)
		self.file_button.move(self.x3, self.y2-self.yspc)
		#Here we define, that when this button is clicked, it starts subfunction choose_file
		QtCore.QObject.connect(self.file_button, QtCore.SIGNAL("clicked()"), self.choose_file)
		
		pdb= QtGui.QLabel('Name of pdb file', self)
		pdb.move(self.x1,self.y2)
		self.pdbfileedit=QtGui.QLineEdit(self)
		self.pdbfileedit.move(self.x2,self.y2)
		self.pdbfileedit.setText(self.savedparmsdict['pdbfile'])
		self.y2 += 30
		
		output= QtGui.QLabel('EM Output file', self)
		output.move(self.x1,self.y2)
		self.outputedit=QtGui.QLineEdit(self)
		self.outputedit.move(self.x2,self.y2)
		self.outputedit.setText(self.savedparmsdict['output'])		
		self.outinfobtn = QPushButton("Output Info", self)
		self.outinfobtn.move(self.x3,  self.y2-self.yspc)
		#sets an infotip for this Pushbutton
		self.outinfobtn.setToolTip('Output Info')
		#when this button is clicked, this action starts the subfunction twodali
		self.connect(self.outinfobtn, SIGNAL("clicked()"), self.outputinfo_pdb2em)
		self.y2 += 30

		apix= QtGui.QLabel('Pixel size (Angstroms)', self)
		apix.move(self.x1,self.y2)
		self.apixedit=QtGui.QLineEdit(self)
		self.apixedit.move(self.x2,self.y2)
		self.apixedit.setText(self.savedparmsdict['apix'])
		self.apixedit.setToolTip('Angstrom/voxel')		
		self.y2 += 30

		box= QtGui.QLabel('Box size (voxels)', self)
		box.move(self.x1,self.y2)
		self.boxedit=QtGui.QLineEdit(self)
		self.boxedit.move(self.x2,self.y2)
		self.boxedit.setText(self.savedparmsdict['box'])
		self.boxedit.setToolTip('Box size in pixels, <xyz> or <x,y,z>')		
		self.y2 += 30

		het= QtGui.QLabel('Include HET atoms', self)
		het.move(self.x1,self.y2)
		self.hetchkbx = QtGui.QCheckBox("",self)
		self.hetchkbx.move(self.x2, self.y2)
		self.hetchkbx.setCheckState(self.savedparmsdict['het'])
		self.hetchkbx.setToolTip('Include HET atoms in the map.')		
		self.y2 += 30

		center= QtGui.QLabel('Center atomic model', self)
		center.move(self.x1,self.y2)
		self.centeredit=QtGui.QLineEdit(self)
		self.centeredit.move(self.x2,self.y2)
		self.centeredit.setText(self.savedparmsdict['center'])
		self.centeredit.setToolTip('center: c - coordinates; a - center of gravity; <x,y,z> - a vector (in Angstrom) to substract from all coordinates; default: n - no')		
		self.y2 += 30
		
		Och= QtGui.QLabel('Use O system of coordinates', self)
		Och.move(self.x1,self.y2)
		self.Ochchkbx = QtGui.QCheckBox("",self)
		self.Ochchkbx.move(self.x2, self.y2)
		self.Ochchkbx.setCheckState(self.savedparmsdict['Och'])
		self.Ochchkbx.setToolTip('apply additional rotation so the model will appear in O in the same rotation as in chimera.')		
		self.y2 += 30
		
		quiet= QtGui.QLabel('Do not print to monitor', self)
		quiet.move(self.x1,self.y2)
		self.quietchkbx = QtGui.QCheckBox("",self)
		self.quietchkbx.move(self.x2, self.y2)
		self.quietchkbx.setCheckState(self.savedparmsdict['quiet'])
		self.quietchkbx.setToolTip('Verbose is the default')		
		self.y2 += 30
		
		tr0= QtGui.QLabel('Transformation matrix to apply to PDB', self)
		tr0.move(self.x1,self.y2)
		self.tr0edit=QtGui.QLineEdit(self)
		self.tr0edit.move(self.x2,self.y2)
		self.tr0edit.setText(self.savedparmsdict['tr0'])
		self.tr0edit.setToolTip('Filename of initial 3x4 transformation matrix')
		
		self.file_button = QtGui.QPushButton("Open File", self)
		self.file_button.move(self.x3, self.y2-self.yspc)
		#Here we define, that when this button is clicked, it starts subfunction choose_file
		QtCore.QObject.connect(self.file_button, QtCore.SIGNAL("clicked()"), self.choose_file1)
		
		self.y2 += 30
		# make ctf, normalize and init_method radio button...
		
		self.savepbtn = QPushButton("Save Input Parameters", self)
		self.savepbtn.move(self.x1-5, self.y4)
		#sets an infotip for this Pushbutton
		self.savepbtn.setToolTip('Save Input Parameters')
		#when this button is clicked, this action starts the subfunction twodali
		self.connect(self.savepbtn, SIGNAL("clicked()"), self.saveparms)
		self.y4+=30
		
		self.cmdlinebtn = QPushButton("Generate command line from input parameters", self)
		self.cmdlinebtn.move(self.x1-5, self.y4)
		#sets an infotip for this Pushbutton
		self.cmdlinebtn.setToolTip('Generate command line using input parameters')
		#when this button is clicked, this action starts the subfunction twodali
		self.connect(self.cmdlinebtn, SIGNAL("clicked()"), self.gencmdline_pdb2em)

		#######################################################################
		 #Here we create a Button(Run_button with title run sxali2d) and its position in the window
		self.RUN_button = QtGui.QPushButton('Run sxpdb2em', self)
		# make 3D textured push button look
		s = "QPushButton {font: bold; color: #000;border: 1px solid #333;border-radius: 11px;padding: 2px;background: qradialgradient(cx: 0, cy: 0,fx: 0.5, fy:0.5,radius: 1, stop: 0 #fff, stop: 1 #8D0);min-width:90px;margin:5px} QPushButton:pressed {font: bold; color: #000;border: 1px solid #333;border-radius: 11px;padding: 2px;background: qradialgradient(cx: 0, cy: 0,fx: 0.5, fy:0.5,radius: 1, stop: 0 #fff, stop: 1 #084);min-width:90px;margin:5px}"
		
		self.RUN_button.setStyleSheet(s)
		self.RUN_button.move(self.x5, self.y5)
		#Here we define, that when this button is clicked, it starts subfunction runsxali2d
		self.connect(self.RUN_button, SIGNAL("clicked()"), self.runsxpdb2em)
		#Labels and Line Edits for User Input		

				
	 #Function runsxali2d started when  the  RUN_button of the  Poptwodali window is clicked 
	def outputinfo_pdb2em(self):
		QMessageBox.information(self, "sxpdb2em output",'output 3-D electron density map (any EM format).\n Attributes apix_x,apix_y,apix_z will be set to the specified value.')
		
	def gencmdline_pdb2em(self,writefile=True):
		#Here we just read in all user inputs in the line edits of the Poptwodali window
		pdbfile = self.pdbfileedit.text()
		output  = self.outputedit.text()
		apix    = self.apixedit.text()
		box     = self.boxedit.text()
		center  = self.centeredit.text()
		tr0     = self.tr0edit.text()
		het     = self.hetchkbx.checkState()
		quiet   = self.quietchkbx.checkState()
		Och     = self.Ochchkbx.checkState()

		cmd1 = "sxpdb2em.py "+str(pdbfile) +" "+ str(output)

		args = " --apix="+ str(apix)+" --box="+ str(box)+ " --center="+str(center)

		cmd1 += args

		if het == Qt.Checked:
			cmd1 += " --het"
		if quiet == Qt.Checked:
			cmd1 += " --quiet"
		if Och == Qt.Checked:
			cmd1 += " --O"				
		if len(str(tr0)) > 0:
			cmd1 += " --tr0="+str(tr0)
		
		self.savedparmsdict = {'pdbfile':str(pdbfile),'output':str(output),'apix':str(apix),'box':str(box),'het':het,'center':str(center),'Och':Och,'quiet':quiet,'tr0':str(tr0)}

		if writefile:		
			(fname,stat)= QInputDialog.getText(self,"Generate Command Line","Enter name of file to save command line in",QLineEdit.Normal,"")
			if stat:
				f = open(fname,'a')
				f.write(cmd1)
				f.write('\n')
				f.close()

		print cmd1
		self.cmd = cmd1

	def runsxpdb2em(self):
		self.gencmdline_pdb2em(writefile=False)
		process = subprocess.Popen(self.cmd,shell=True)
		self.emit(QtCore.SIGNAL("process_started"),process.pid)
		
	def saveparms(self):		
		# save all the parms in a text file so we can repopulate if user requests
		(fname,stat)= QInputDialog.getText(self,"Save Input Parameters","Enter name of file to save parameters in",QLineEdit.Normal,"")
		if stat:
			import pickle
			output=open(fname,'wb')
			self.gencmdline_pdb2em(writefile=False)
			pickle.dump(self.savedparmsdict,output)
			output.close()
		
	def repoparms_pdb2em(self):		
		(fname,stat)= QInputDialog.getText(self,"Get Input Parameters","Enter name of file parameters were saved in",QLineEdit.Normal,"")
		if stat:
			import pickle
			pkl = open(fname,'rb')
			self.savedparmsdict = pickle.load(pkl)
			self.pdbfileedit.setText(self.savedparmsdict['pdbfile'])
			self.outputedit.setText(self.savedparmsdict['output'])
			self.apixedit.setText(self.savedparmsdict['apix'])		
			self.boxedit.setText(self.savedparmsdict['box'])		
			self.centeredit.setText(self.savedparmsdict['center'])
			self.hetchkbx.setCheckState(self.savedparmsdict['het'])
			self.Ochchkbx.setCheckState(self.savedparmsdict['Och'])
			self.quietchkbx.setCheckState(self.savedparmsdict['quiet'])
			self.tr0edit.setText(self.savedparmsdict['tr0'])
				
	def choose_file(self):
		#opens a file browser, showing files only in .hdf format
		file_name = QtGui.QFileDialog.getOpenFileName(self, "Open PDB File", "", "PDB files (*.pdb)")
		#after the user selected a file, we obtain this filename as a Qstring
		a=QtCore.QString(file_name)
		#we convert this Qstring to a string and send it to line edit classed stackname edit of the Poptwodali window
		self.pdbfileedit.setText(str(a))
	
	
	def choose_file1(self):
		#opens a file browser, showing files only in .hdf format
		file_name = QtGui.QFileDialog.getOpenFileName(self, "Open file containing transformation matrix", "", "(*)")
		#after the user selected a file, we obtain this filename as a Qstring
		a=QtCore.QString(file_name)
		#we convert this Qstring to a string and send it to line edit classed stackname edit of the Poptwodali window
		self.tr0edit.setText(str(a))
		 
#Layout of the Pop Up window Popuptwodali (for sxali2d); started by the function twodali of the main window		
class Popuppca(QWidget):
	def __init__(self):
		QWidget.__init__(self)
		
		#######################################################################################
		# class variables
		
		self.cmd = ""
		# populate with default values
		self.savedparmsdict = {'input_stack_list':'','output_stack':'','subavg':'','rad':'-1','nvec':'1','mask':'','sdir':'.','usebuf':Qt.Unchecked,'nproc':'1','shuffle':Qt.Unchecked}

		#######################################################################################
		# Layout parameters
		
		self.y1 = 10 # title and Repopulate button
		self.y2 = self.y1 + 78 # Text boxes for inputting parameters
		#self.y3 = self.y2 + 222 # activate images button and set xform.align2d button
		self.y4 = self.y2 + 222 # Advanced Parameters, Save Input and Generate command line buttons
		self.y5 = self.y4 + 95 # run button 
		self.yspc = 4
		
		self.x1 = 10 # first column (text box labels)
		self.x2 = self.x1 + 160 # second column (text boxes)
		self.x3 = self.x2+145 # third column (Open .hdf button)
		self.x4 = self.x3+100 # fourth column (Open .bdb button)
		self.x5 = self.x4+100
		self.x6 = 230 # run button
		#######################################################################################
		
		#Here we just set the window title
		self.setWindowTitle('sxpca')
		#Here we just set a label and its position in the window
		title1=QtGui.QLabel('<b>sxpca</b> - Principal Component Analysis of images', self)
		title1.move(self.x1,self.y1)
		self.y1 += 30

		self.repopbtn = QPushButton("Retrieve saved parameters", self)
		self.repopbtn.move(self.x1-5,self.y1)
		#sets an infotip for this Pushbutton
		self.repopbtn.setToolTip('Retrieve saved parameters')
		#when this button is clicked, this action starts the subfunction twodali
		self.connect(self.repopbtn, SIGNAL("clicked()"), self.repoparms_pca)

		#######################################################################################
		#Here we create a Button(file_button with title run open .hdf) and its position in the window
		self.file_button = QtGui.QPushButton("Open .hdf", self)
		self.file_button.move(self.x3, self.y2-self.yspc)
		#Here we define, that when this button is clicked, it starts subfunction choose_file
		QtCore.QObject.connect(self.file_button, QtCore.SIGNAL("clicked()"), self.choose_file)
		#exactly the same as above, but for subfunction choose_file1
		self.file_button1 = QtGui.QPushButton("Open .bdb", self)
		self.file_button1.move(self.x4,self.y2-self.yspc)
		QtCore.QObject.connect(self.file_button1, QtCore.SIGNAL("clicked()"), self.choose_file1)
		
		inpstklist= QtGui.QLabel('Enter input stacks', self)
		inpstklist.move(self.x1,self.y2)
		self.input_stack_listedit=QtGui.QLineEdit(self)
		self.input_stack_listedit.move(self.x2,self.y2)
		self.input_stack_listedit.setText(self.savedparmsdict['input_stack_list'])
		
		self.inputinfobtn = QPushButton("Input Info", self)
		self.inputinfobtn.move(self.x5,  self.y2-self.yspc)
		#sets an infotip for this Pushbutton
		self.inputinfobtn.setToolTip('Input help')
		#when this button is clicked, this action starts the subfunction twodali
		self.connect(self.inputinfobtn, SIGNAL("clicked()"), self.inputinfo_pca)
		
		self.y2 += 30
		
		outstk= QtGui.QLabel('Output Stack', self)
		outstk.move(self.x1,self.y2)
		self.output_stackedit=QtGui.QLineEdit(self)
		self.output_stackedit.move(self.x2,self.y2)
		self.output_stackedit.setText(self.savedparmsdict['output_stack'])		
		
		self.outinfobtn = QPushButton("Output Info", self)
		self.outinfobtn.move(self.x3,  self.y2-self.yspc)
		#sets an infotip for this Pushbutton
		self.outinfobtn.setToolTip('Output Info')
		#when this button is clicked, this action starts the subfunction twodali
		self.connect(self.outinfobtn, SIGNAL("clicked()"), self.outputinfo_pca)
		self.y2 += 30

		sa= QtGui.QLabel('Average to subtract', self)
		sa.move(self.x1,self.y2)
		self.subavgedit=QtGui.QLineEdit(self)
		self.subavgedit.move(self.x2,self.y2)
		self.subavgedit.setText(self.savedparmsdict['subavg'])
		self.subavgedit.setToolTip('subtract average')		
		
		self.file_button2 = QtGui.QPushButton("Open file", self)
		self.file_button2.move(self.x3, self.y2-self.yspc)
		#Here we define, that when this button is clicked, it starts subfunction choose_file
		QtCore.QObject.connect(self.file_button2, QtCore.SIGNAL("clicked()"), self.choose_file2)
		
		self.y2 += 30
		
		radius= QtGui.QLabel('Radius of mask', self)
		radius.move(self.x1,self.y2)
		self.radedit=QtGui.QLineEdit(self)
		self.radedit.move(self.x2,self.y2)
		self.radedit.setText(self.savedparmsdict['rad'])
		self.radedit.setToolTip('Radius of mask')		
		self.y2 += 30

		nv= QtGui.QLabel('Number of eigenvectors', self)
		nv.move(self.x1,self.y2)
		self.nvecedit=QtGui.QLineEdit(self)
		self.nvecedit.move(self.x2,self.y2)
		self.nvecedit.setText(self.savedparmsdict['nvec'])
		self.nvecedit.setToolTip('Number of eigenvectors')
		self.y2 += 30
		
		nproc= QtGui.QLabel('MPI processors', self)
		nproc.move(self.x1,self.y2)
		self.nprocedit=QtGui.QLineEdit(self)
		self.nprocedit.move(self.x2,self.y2)
		self.nprocedit.setText(self.savedparmsdict['nproc'])
		self.nprocedit.setToolTip('The number of processors to use. Default is single processor mode')

		
		######################################################################################
		
		self.savepbtn = QPushButton("Save Input Parameters", self)
		self.savepbtn.move(self.x1-5, self.y4)
		#sets an infotip for this Pushbutton
		self.savepbtn.setToolTip('Save Input Parameters')
		#when this button is clicked, this action starts the subfunction twodali
		self.connect(self.savepbtn, SIGNAL("clicked()"), self.saveparms)
		self.y4+=30
		
		self.cmdlinebtn = QPushButton("Generate command line from input parameters", self)
		self.cmdlinebtn.move(self.x1-5, self.y4)
		#sets an infotip for this Pushbutton
		self.cmdlinebtn.setToolTip('Generate command line using input parameters')
		#when this button is clicked, this action starts the subfunction twodali
		self.connect(self.cmdlinebtn, SIGNAL("clicked()"), self.gencmdline_pca)

		#######################################################################
		 #Here we create a Button(Run_button with title run sxali2d) and its position in the window
		self.RUN_button = QtGui.QPushButton('Run sxpca', self)
		# make 3D textured push button look
		s = "QPushButton {font: bold; color: #000;border: 1px solid #333;border-radius: 11px;padding: 2px;background: qradialgradient(cx: 0, cy: 0,fx: 0.5, fy:0.5,radius: 1, stop: 0 #fff, stop: 1 #8D0);min-width:90px;margin:5px} QPushButton:pressed {font: bold; color: #000;border: 1px solid #333;border-radius: 11px;padding: 2px;background: qradialgradient(cx: 0, cy: 0,fx: 0.5, fy:0.5,radius: 1, stop: 0 #fff, stop: 1 #084);min-width:90px;margin:5px}"
		
		self.RUN_button.setStyleSheet(s)
		self.RUN_button.move(self.x6, self.y5)
		#Here we define, that when this button is clicked, it starts subfunction runsxali2d
		self.connect(self.RUN_button, SIGNAL("clicked()"), self.runsxpca)
		#Labels and Line Edits for User Input		

		
	 #Function runsxali2d started when  the  RUN_button of the  Poptwodali window is clicked 
	def outputinfo_pca(self):
		QMessageBox.information(self, "sxpca output",'output_stack contains the result of the PCA. saved as stack files')
		
		
	def inputinfo_pca(self):
		QMessageBox.information(self, "sxpca input",'input list of filenames delimited by space, e.g., myfile1.hdf /usr/folder/myfile2.hdf myfile3.hdf')
		
	def gencmdline_pca(self,writefile=True):
		#Here we just read in all user inputs in the line edits of the Poptwodali window
		input_stack_list = self.input_stack_listedit.text()
		output_stack=self.output_stackedit.text()
		rad=self.radedit.text()
		nvec=self.nvecedit.text()
		subavg=self.subavgedit.text()
		
		inpstklist=str(input_stack_list).split(' ')
		
		cmd1 = "sxpca.py"
		
		for fn in inpstklist:
			if len(str(fn))>0:
				cmd1 = cmd1+" "+fn
		cmd1 = cmd1+" "+str(output_stack)
		
		args = " --rad="+ str(rad)+ " --subavg='"+str(subavg)+"'"+" --nvec="+str(nvec)
		
		cmd1 = cmd1 + args
		
		mask = self.w1.maskedit.text()
		sdir = self.w1.sdiredit.text()
		usebuf=self.w1.usebufchkbx.checkState()
		shuffle=self.w1.shufflechkbx.checkState()				
						
		cmd1 = cmd1+" --mask='" +str(mask) +"'"+" --sdir="+str(sdir)
		
		if usebuf == Qt.Checked:
			cmd1 = cmd1 + " --usebuf"
		if shuffle == Qt.Checked:
			cmd1 = cmd1 + " --shuffle"		
		
		np = self.nprocedit.text()
		
		self.savedparmsdict = {'input_stack_list':str(input_stack_list),'output_stack':str(output_stack),'subavg':str(subavg),'rad':str(rad),'nvec':str(nvec),'mask':str(mask),'sdir':str(sdir),'usebuf':usebuf,'nproc':str(np),'shuffle':shuffle}

		self.w1.savedparmsdict=self.savedparmsdict
		
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
		
	def runsxpca(self):
		self.gencmdline_pca(writefile=False)
		process = subprocess.Popen(self.cmd,shell=True)
		self.emit(QtCore.SIGNAL("process_started"),process.pid)
		
	def saveparms(self):		
		(fname,stat)= QInputDialog.getText(self,"Save Input Parameters","Enter name of file to save parameters in",QLineEdit.Normal,"")
		if stat:
			import pickle
			output=open(fname,'wb')
			self.gencmdline_pca(writefile=False)
			pickle.dump(self.savedparmsdict,output)
			output.close()
		
	def repoparms_pca(self):		
		(fname,stat)= QInputDialog.getText(self,"Retrieve saved parameters","Enter name of file parameters were saved in",QLineEdit.Normal,"")
		if stat:
			import pickle
			pkl = open(fname,'rb')
			self.savedparmsdict = pickle.load(pkl)
			self.input_stack_listedit.setText(self.savedparmsdict['input_stack_list'])
			self.output_stackedit.setText(self.savedparmsdict['output_stack'])
			self.radedit.setText(self.savedparmsdict['rad'])
			self.subavgedit.setText(self.savedparmsdict['subavg'])
			self.nvecedit.setText(self.savedparmsdict['nvec'])
			self.nprocedit.setText(self.savedparmsdict['nproc'])
			self.w1.maskedit.setText(self.savedparmsdict['mask'])
			self.w1.sdiredit.setText(self.savedparmsdict['sdir'])
			self.w1.usebufchkbx.setCheckState(self.savedparmsdict['usebuf'])
			self.w1.shufflechkbx.setCheckState(self.savedparmsdict['shuffle'])
				
		#Function choose_file started when  the  open_file of the  Poptwodali window is clicked
	def choose_file(self):
		#opens a file browser, showing files only in .hdf format
		file_name = QtGui.QFileDialog.getOpenFileName(self, "Open Data File", "", "HDF files (*.hdf)")
		#after the user selected a file, we obtain this filename as a Qstring
		a=QtCore.QString(file_name)
		print a
		#we convert this Qstring to a string and send it to line edit classed stackname edit of the Poptwodali window
		curlist = self.input_stack_listedit.text()
		curlist = curlist+' '+str(a)
		self.input_stack_listedit.setText(str(curlist))
		
		#Function choose_file started when  the  open_file of the  Poptwodali window is clicked (same as above but for bdb files(maybe we can combine these two into one function)
	def choose_file1(self):
		file_name1 = QtGui.QFileDialog.getOpenFileName(self, "Open Data File", "EMAN2DB/", "BDB FILES (*.bdb)" )
		a=QtCore.QString(file_name1)
		b=os.path.basename(str(a))
		c=os.path.splitext(b)[0]
		d="bdb:"+c
		print d
		curlist = self.input_stack_listedit.text()
		curlist = curlist+' '+str(d)
		self.input_stack_listedit.setText(str(curlist))

	def choose_file2(self):
		#opens a file browser, showing files only in .hdf format
		file_name = QtGui.QFileDialog.getOpenFileName(self, "Open Data File", "", "files (*)")
		#after the user selected a file, we obtain this filename as a Qstring
		a=QtCore.QString(file_name)
		print a
		#we convert this Qstring to a string and send it to line edit classed stackname edit of the Poptwodali window
		self.subavgedit.setText(str(a))
				
class Popupadvparams_pca(QWidget):
	def __init__(self,savedparms):
		QWidget.__init__(self)
		
		self.x1 = 10
		self.x2 = 140
		self.x3 = 285
		
		self.y1 = 10
		self.yspc = 4
		
		#Here we just set the window title
		self.setWindowTitle('sxpca advanced parameter selection')
		#Here we just set a label and its position in the window
		title1=QtGui.QLabel('<b>sxpca</b> - set advanced params', self)
		title1.move(self.x1,self.y1)
		#Labels and Line Edits for User Input
		#Just a label
		self.y1 += 30
		
		title2= QtGui.QLabel('<b>Advanced</b> parameters', self)
		title2.move(self.x1,self.y1)
		
		self.y1 += 30
		
		self.savedparmsdict=savedparms
		#Example for User input stack name
		#First create the label and define its position
		maskname= QtGui.QLabel('Mask', self)
		maskname.move(self.x1,self.y1)
		#Now add a line edit and define its position
		self.maskedit=QtGui.QLineEdit(self)
		self.maskedit.move(self.x2,self.y1)
		#Adds a default value for the line edit
		self.maskedit.setText(self.savedparmsdict['mask'])
		self.maskedit.setToolTip("mask file")
		
		self.mskfile_button = QtGui.QPushButton("Open File", self)
		self.mskfile_button.move(self.x3, self.y1-self.yspc)
		#Here we define, that when this button is clicked, it starts subfunction choose_file
		QtCore.QObject.connect(self.mskfile_button, QtCore.SIGNAL("clicked()"), self.choose_mskfile)
		
		self.y1 += 30
		
		sd= QtGui.QLabel('Scratch directory', self)
		sd.move(self.x1,self.y1)
		self.sdiredit=QtGui.QLineEdit(self)
		self.sdiredit.move(self.x2,self.y1)
		self.sdiredit.setText(self.savedparmsdict['sdir'])
		self.sdiredit.setToolTip('Scratch directory')
		
		self.y1 += 30
		
		useb= QtGui.QLabel('Use existing buffer', self)
		useb.move(self.x1,self.y1)
		self.usebufchkbx = QtGui.QCheckBox("",self)
		self.usebufchkbx.move(self.x2, self.y1)
		self.usebufchkbx.setCheckState(self.savedparmsdict['usebuf'])
		
		self.y1 += 30
				
		shuf= QtGui.QLabel('Use shuffle', self)
		shuf.move(self.x1,self.y1)
		self.shufflechkbx=QtGui.QCheckBox("",self)
		self.shufflechkbx.move(self.x2,self.y1)
		self.shufflechkbx.setCheckState(self.savedparmsdict['shuffle'])
		self.shufflechkbx.setToolTip('use shuffle')
		
		
	def choose_mskfile(self):
		#opens a file browser, showing files only in .hdf format
		file_name = QtGui.QFileDialog.getOpenFileName(self, "Open File Containing Mask", "", "(*)")
		#after the user selected a file, we obtain this filename as a Qstring
		a=QtCore.QString(file_name)
		print a
		#we convert this Qstring to a string and send it to line edit classed stackname edit of the Poptwodali window
		self.maskedit.setText(str(a))
				
class Popupmrefthreedali(QWidget):
	def __init__(self):
		QWidget.__init__(self)
		
		########################################################################################
		# class variables
		self.cmd = ""
		# populate with default values
		self.savedparmsdict = {'stackname':'','refname':'','foldername':'','partradius':'-1','xyrange':'4 2 1 1','trans':'2 1 0.5 0.25', 'delta':'15 5 2','nriter':'3','nproc':'1','maskname':'','focus':'','center':'-1',"ringstep":"1","innerradius":"1","ctf":Qt.Unchecked,"snr":"1.0","fourvar":Qt.Unchecked,"usrfunc":"ref_ali3dm","usrfuncfile":"","an":'-1',"ref_a":'S',"sym":'c1',"npad":'4',"stoprnct":'0.0',"debug":False,'nrefine':'1','nassign':'0'}
		
		########################################################################################
		# layout parameters
		
		self.y1 = 10
		self.y2 = self.y1 + 98
		self.y3 = self.y2 + 282
		self.y4 = self.y3 + 80
		self.y5 = self.y4 + 95
		self.yspc = 4
		
		self.x1 = 10
		self.x2 = self.x1 + 150
		self.x3 = self.x2+145
		self.x4 = self.x3+100
		self.x5 = 230
		########################################################################################
		
		#Here we just set the window title
		self.setWindowTitle('sxmref_ali3d')
		#Here we just set a label and its position in the window
		title1=QtGui.QLabel('<b>sxmref_ali3d</b> -  perform 3-D multireference projection matching given initial<br>reference volumes and image series', self)
		title1.move(self.x1,self.y1)
		self.y1 += 50

		self.repopbtn = QPushButton("Retrieve saved parameters", self)
		self.repopbtn.move(self.x1-5,self.y1)
		#sets an infotip for this Pushbutton
		self.repopbtn.setToolTip('Retrieve saved parameters')
		#when this button is clicked, this action starts the subfunction twodali
		self.connect(self.repopbtn, SIGNAL("clicked()"), self.repoparms_mrefali3d)
		
		########################################################################################
		
		stackname= QtGui.QLabel('Name of input stack', self)
		stackname.move(self.x1,self.y2)
		#Now add a line edit and define its position
		self.stacknameedit=QtGui.QLineEdit(self)
		self.stacknameedit.move(self.x2,self.y2)
		#Adds a default value for the line edit
		self.stacknameedit.setText(self.savedparmsdict['stackname'])
		
		#Here we create a Button(file_button with title run open .hdf) and its position in the window
		self.file_button = QtGui.QPushButton("Open .hdf", self)
		self.file_button.move(self.x3, self.y2-self.yspc)
		#Here we define, that when this button is clicked, it starts subfunction choose_file
		QtCore.QObject.connect(self.file_button, QtCore.SIGNAL("clicked()"), self.choose_file)
		#exactly the same as above, but for subfunction choose_file1
		self.file_button1 = QtGui.QPushButton("Open .bdb", self)
		self.file_button1.move(self.x4,self.y2-self.yspc)
		QtCore.QObject.connect(self.file_button1, QtCore.SIGNAL("clicked()"), self.choose_file1)
		
		self.y2 += 30
				
		refname= QtGui.QLabel('Name of reference', self)
		refname.move(self.x1,self.y2)
		#Now add a line edit and define its position
		self.refnameedit=QtGui.QLineEdit(self)
		self.refnameedit.move(self.x2,self.y2)
		#Adds a default value for the line edit
		self.refnameedit.setText(self.savedparmsdict['refname'])
		#Here we create a Button(file_button with title run open .hdf) and its position in the window
		self.ref_button = QtGui.QPushButton("Open files", self)
		self.ref_button.move(self.x3, self.y2-self.yspc)
		#Here we define, that when this button is clicked, it starts subfunction choose_file
		QtCore.QObject.connect(self.ref_button, QtCore.SIGNAL("clicked()"), self.choose_reffile)

		
		self.y2 = self.y2+30
		
		#The same as above, but many line edits include Infotips
		foldername= QtGui.QLabel('Output folder', self)
		foldername.move(self.x1,self.y2)
		self.foldernameedit=QtGui.QLineEdit(self)
		self.foldernameedit.move(self.x2,self.y2)
		self.foldernameedit.setText(self.savedparmsdict['foldername'])		
		
		self.outinfobtn = QPushButton("Output Info", self)
		self.outinfobtn.move(self.x3,  self.y2-self.yspc)
		#sets an infotip for this Pushbutton
		self.outinfobtn.setToolTip('Output Info')
		#when this button is clicked, this action starts the subfunction twodali
		self.connect(self.outinfobtn, SIGNAL("clicked()"), self.outputinfo_mrefali3d)
		
		self.y2 = self.y2+30
		
		partradius= QtGui.QLabel('Particle radius', self)
		partradius.move(self.x1,self.y2)
		self.partradiusedit=QtGui.QLineEdit(self)
		self.partradiusedit.move(self.x2,self.y2)
		self.partradiusedit.setText(self.savedparmsdict['partradius'])
		self.partradiusedit.setToolTip('Parameter ou: Outer radius for rotational correlation \nshould be set to particle radius\nif not sure, set to boxsize/2-2 ')		

		self.y2 = self.y2+30
		
		xyrange= QtGui.QLabel('xy range', self)
		xyrange.move(self.x1,self.y2)
		self.xyrangeedit=QtGui.QLineEdit(self)
		self.xyrangeedit.move(self.x2,self.y2)
		self.xyrangeedit.setText(self.savedparmsdict['xyrange'])
		self.xyrangeedit.setToolTip('Range for translational search in x, y direction\nif set to 0 only rotational alignment will be performed')

		self.y2 = self.y2+30
		
		trans= QtGui.QLabel('translational step', self)
		trans.move(self.x1,self.y2)
		self.transedit=QtGui.QLineEdit(self)
		self.transedit.move(self.x2,self.y2)
		self.transedit.setText(self.savedparmsdict['trans'])
		self.transedit.setToolTip('Step of translational search in x, y direction\nlarger values increase the speed but decrease the accuracy')		

		self.y2 = self.y2+30
		
		delta= QtGui.QLabel('angular step', self)
		delta.move(self.x1,self.y2)
		self.deltaedit=QtGui.QLineEdit(self)
		self.deltaedit.move(self.x2,self.y2)
		self.deltaedit.setText(self.savedparmsdict['delta'])
		self.deltaedit.setToolTip('angular step for the reference projections in respective iterations')				
		
		self.y2 =self.y2+30

		nriter= QtGui.QLabel('Number of iterations', self)
		nriter.move(self.x1,self.y2)
		self.nriteredit=QtGui.QLineEdit(self)
		self.nriteredit.move(self.x2,self.y2)
		self.nriteredit.setText(self.savedparmsdict['nriter'])
		self.nriteredit.setToolTip('Maximum number of iterations the program will perform\n Using the default values the program will run 3 rounds with xy-range 4 and translational step 1, 3 rounds with xyrange 2 and translational step 1 and so on..\nif set to 0 maximum iteration number will be 10 and will automatically stop should the criterion falls')
		
		self.y2 =self.y2+30
		
		nproc= QtGui.QLabel('MPI processors', self)
		nproc.move(self.x1,self.y2)
		self.nprocedit=QtGui.QLineEdit(self)
		self.nprocedit.move(self.x2,self.y2)
		self.nprocedit.setText(self.savedparmsdict['nproc'])
		self.nprocedit.setToolTip('The number of processors to use. Default is single processor mode')
		
		########################################################################################
				
		header=QtGui.QLabel('Attributes xform.projection and active parameters must be set in the input stack', self)
		header.move(self.x1,self.y3)
		self.y3 = self.y3 + 30
		
		#not linked to a function yet
		self.activeheader_button = QtGui.QPushButton("activate all images", self)
		self.activeheader_button.move(self.x1-5, self.y3)
		self.connect(self.activeheader_button, SIGNAL("clicked()"), self.setactiveheader)
		
		self.projheader_button = QtGui.QPushButton("set xform.projection", self)
		self.projheader_button.move(self.x1-5+180, self.y3)
		self.connect(self.projheader_button, SIGNAL("clicked()"), self.setprojheader)
		
		########################################################################################
		
		self.savepbtn = QPushButton("Save Input Parameters", self)
		self.savepbtn.move(self.x1-5,  self.y4)
		#sets an infotip for this Pushbutton
		self.savepbtn.setToolTip('Save Input Parameters')
		#when this button is clicked, this action starts the subfunction twodali
		self.connect(self.savepbtn, SIGNAL("clicked()"), self.saveparms)
		
		self.y4 = self.y4+30
		
		self.cmdlinebtn = QPushButton("Generate command line from input parameters", self)
		self.cmdlinebtn.move(self.x1-5,  self.y4)
		#sets an infotip for this Pushbutton
		self.cmdlinebtn.setToolTip('Generate command line using input parameters')
		#when this button is clicked, this action starts the subfunction twodali
		self.connect(self.cmdlinebtn, SIGNAL("clicked()"), self.gencmdline_mrefali3d)
		
		########################################################################################
		
		 #Here we create a Button(Run_button with title run sxali2d) and its position in the window
		self.RUN_button = QtGui.QPushButton('Run sxmref_ali3d', self)
		# make 3D textured push button look
		s = "QPushButton {font: bold; color: #000;border: 1px solid #333;border-radius: 11px;padding: 2px;background: qradialgradient(cx: 0, cy: 0,fx: 0.5, fy:0.5,radius: 1, stop: 0 #fff, stop: 1 #8D0);min-width:90px;margin:5px} QPushButton:pressed {font: bold; color: #000;border: 1px solid #333;border-radius: 11px;padding: 2px;background: qradialgradient(cx: 0, cy: 0,fx: 0.5, fy:0.5,radius: 1, stop: 0 #fff, stop: 1 #084);min-width:90px;margin:5px}"
		
		self.RUN_button.setStyleSheet(s)
		
		self.RUN_button.move(self.x5,  self.y5)
		#Here we define, that when this button is clicked, it starts subfunction runsxali2d
		self.connect(self.RUN_button, SIGNAL("clicked()"), self.runsxmrefali3d)
		#Labels and Line Edits for User Input

		
		#Function runsxali2d started when  the  RUN_button of the  Poptwodali window is clicked 

	def outputinfo_mrefali3d(self):
		QMessageBox.information(self, "mref_ali3d output",'Output directory \n\n directory name into which the output files will be written. If it does not exist, the directory will be created. If it does exist, the program will crash and an error message will come up. Please change the name of directory and restart the program . The program will write the low-pass filtered reconstructed volume volf****.spi, and Fourier ring correlation text file resolution****. \n\n header\n\nthe alignment parameters are stored in the headers of input files as xform.projection and assignment to reference volumes as group (it is an index ranging from zero to number of templates minus one).')
		
	def gencmdline_mrefali3d(self,writefile=True):
		#Here we just read in all user inputs in the line edits of the Poptwodali window
		stack = self.stacknameedit.text()
		print "stack defined="+ stack
		ref = self.refnameedit.text()
		print "reference structure ="+ stack
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
		delta=self.deltaedit.text()
		print "delta=" +delta
		maxit=self.nriteredit.text()
		print "maxit="+maxit
		
		cmd1 = "sxmref_ali3d.py "+str(stack) +" "+ str(ref)+" "+ str(output)
		
		args = " --ou="+ str(ou)+ " --xr='"+str(xr)+"'"+ " --yr='"+str(yr)+"'"+ " --ts='"+str(ts)+"'"+ " --delta='"+str(delta)+"'"+" --maxit="+ str(maxit) 
		
		
		mask = self.w1.masknameedit.text()
		if len(str(mask))> 1:
			cmd1 = cmd1+" "+str(mask) 
		cmd1 = cmd1 + args
		
		
		ctr=self.w1.centeredit.text()
		ringstep = self.w1.ringstepedit.text()
		inrad = self.w1.innerradiusedit.text()
		CTF=self.w1.ctfchkbx.checkState()
		snr = self.w1.snredit.text()
		userf = self.w1.usrfuncedit.text()
		userfile = self.w1.usrfuncfileedit.text()
		an = self.w1.anedit.text()
		ref_a = self.w1.refaedit.text()
		sym = self.w1.symedit.text()
		npad = self.w1.npadedit.text()
		nassign = self.w1.nassignedit.text()
		nrefine = self.w1.nrefineedit.text()		
		focus = self.w1.focusedit.text()
				
		fourvar = self.w1.fourvarchkbx.checkState()
		debug = self.w1.debugchkbx.checkState()
		stoprnct = self.w1.stoprnctedit.text()
		
				
		cmd1 = cmd1+" --center=" +str(ctr)+" --rs="+str(ringstep)+ " --ir=" + str(inrad)+ " --snr=" + str(snr)+ " --an='" + str(an)+ "'"+" --sym=" + str(sym)+ " --ref_a=" + str(ref_a)+ " --npad=" + str(npad)+ " --stoprnct=" + str(stoprnct)+ " --nassign=" + str(nassign)+ " --nrefine=" + str(nrefine)
		if CTF == Qt.Checked:
			cmd1 = cmd1 + " --CTF"
		if fourvar == Qt.Checked:
			cmd1 = cmd1 + " --fourvar"		
		
		if debug == Qt.Checked:
			cmd1 = cmd1 + " --debug"
		if len(str(focus)) > 0:
			cmd1 = cmd1 + " --focus='"+str(focus)+"'"
												
		if len(userfile) < 1:
			cmd1 = cmd1 + " --function="+str(userf)
		else:
			userfile = str(userfile)
			# break it down into file name and directory path
			rind = userfile.rfind('/')
				
			if rind == -1:
					userfile = os.path.abspath(userfile)
					rind = userfile.rfind('/')
					
			fname = userfile[rind+1:]
			fname, ext = os.path.splitext(fname)
			fdir = userfile[0:rind]
			cmd1 = cmd1 + " --function=\"[" +fdir+","+fname+","+str(userf)+"]\""
				
		np = self.nprocedit.text()
		
		self.savedparmsdict = {'stackname':str(stack),'refname':str(ref),'foldername':str(output),'partradius':str(ou),'xyrange':str(xr),'trans':str(yr),'delta':str(delta),'nriter':str(maxit),'nproc':str(np),'maskname':str(mask),'focus':str(focus),'center':str(ctr),"ringstep":str(ringstep),"innerradius":str(inrad),"ctf":CTF,"snr":str(snr),"fourvar":fourvar,"usrfunc":str(userf), "usrfuncfile":str(userfile),"an":str(an),"ref_a":str(ref_a),"sym":str(sym),"npad":str(npad),"stoprnct":str(stoprnct),"debug":debug,'nassign':str(nassign),'nrefine':str(nrefine)}

		self.w1.savedparmsdict=self.savedparmsdict
		
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
		
	def runsxmrefali3d(self):
		self.gencmdline_mrefali3d(writefile=False)
		outfolder=self.savedparmsdict['foldername']
		if os.path.exists(outfolder):
			print "output folder "+outfolder+" already exists!"
			return
		process = subprocess.Popen(self.cmd,shell=True)
		self.emit(QtCore.SIGNAL("process_started"),process.pid)
		
	def saveparms(self):		
		(fname,stat)= QInputDialog.getText(self,"Save Input Parameters","Enter name of file to save parameters in",QLineEdit.Normal,"")
		if stat:
			import pickle
			output=open(fname,'wb')
			self.gencmdline_mrefali3d(writefile=False)
			pickle.dump(self.savedparmsdict,output)
			output.close()
		
	def repoparms_mrefali3d(self):		
		# repopulate with saved parms
		(fname,stat)= QInputDialog.getText(self,"Retrieve saved parameters","Enter name of file parameters were saved in",QLineEdit.Normal,"")
		if stat:
			import pickle
			pkl = open(fname,'rb')
			self.savedparmsdict = pickle.load(pkl)
			print self.savedparmsdict
			self.partradiusedit.setText(self.savedparmsdict['partradius'])
			self.stacknameedit.setText(self.savedparmsdict['stackname'])
			self.refnameedit.setText(self.savedparmsdict['refname'])		
			self.foldernameedit.setText(self.savedparmsdict['foldername'])		
			self.xyrangeedit.setText(self.savedparmsdict['xyrange'])
			self.transedit.setText(self.savedparmsdict['trans'])
			self.deltaedit.setText(self.savedparmsdict['delta'])
			self.nriteredit.setText(self.savedparmsdict['nriter'])
			self.nprocedit.setText(self.savedparmsdict['nproc'])
	
			self.w1.masknameedit.setText(self.savedparmsdict['maskname'])
			self.w1.focusedit.setText(self.savedparmsdict['focus'])
			self.w1.centeredit.setText(self.savedparmsdict['center'])
			self.w1.ringstepedit.setText(self.savedparmsdict['ringstep'])
			self.w1.innerradiusedit.setText(self.savedparmsdict['innerradius'])
			self.w1.ctfchkbx.setCheckState(self.savedparmsdict['ctf'])
			self.w1.snredit.setText(self.savedparmsdict['snr'])
			self.w1.usrfuncedit.setText(self.savedparmsdict['usrfunc'])
			self.w1.usrfuncfileedit.setText(self.savedparmsdict['usrfuncfile'])
			self.w1.anedit.setText(self.savedparmsdict['an'])
			self.w1.refaedit.setText(self.savedparmsdict['ref_a'])
			self.w1.symedit.setText(self.savedparmsdict['sym'])
			self.w1.npadedit.setText(self.savedparmsdict['npad'])
			self.w1.nassignedit.setText(self.savedparmsdict['nassign'])
			self.w1.nrefineedit.setText(self.savedparmsdict['nrefine'])
			
			self.w1.fourvarchkbx.setCheckState(self.savedparmsdict['fourvar'])
			self.w1.debugchkbx.setCheckState(self.savedparmsdict['debug'])
			self.w1.stoprnctedit.setText(self.savedparmsdict['stoprnct'])
				
				
	def setactiveheader(self):
		stack = self.stacknameedit.text()
		print "stack defined="+ stack
		# horatio active_refactoring Jy51i1EwmLD4tWZ9_00000_1
		# header(str(stack), "active", one=True)


	def setprojheader(self):
		#opens a file browser, showing files only in .hdf format
		ok=False
		zerostr='set xform.projection to zero'
		randstr='randomize xform.projection'
		importstr='import parameters from file'
		(item,stat)= QInputDialog.getItem(self,"xform.projection","choose option",[zerostr,randstr,importstr])
		#we convert this Qstring to a string and send it to line edit classed stackname edit of the Poptwodali window
		#self.stacknameedit.setText(str(a))
		choice= str(item)
		stack = self.stacknameedit.text()
		if stat:
			if choice == zerostr:
				header(str(stack),'xform.projection',zero=True)
			if choice == randstr:
				header(str(stack),'xform.projection',rand_alpha=True)
			if choice == importstr:
				file_name = QtGui.QFileDialog.getOpenFileName(self, "Open Data File", "", "(*)")
				a=str(QtCore.QString(file_name))
				if len(a)>0:
					header(str(stack),'xform.projection',fimport=a)
		
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

	def choose_reffile(self):
		#opens a file browser, showing files only in .hdf format
		file_name = QtGui.QFileDialog.getOpenFileName(self, "Open reference structure file", "", "(*)")
		#after the user selected a file, we obtain this filename as a Qstring
		a=QtCore.QString(file_name)
		print a
		#we convert this Qstring to a string and send it to line edit classed stackname edit of the Poptwodali window
		self.refnameedit.setText(str(a))


class Popupadvparams_mrefali3d_1(QWidget):
	def __init__(self,savedparms):
		QWidget.__init__(self)
		
		########################################################################################
		# layout parameters
		
		self.y1=10
		self.yspc = 4
		
		self.x1 = 10
		self.x2 = self.x1+280
		self.x3 = self.x2+145
		########################################################################################
		
		#Here we just set the window title
		#self.setWindowTitle('sxali3d advanced parameter selection')
		#Here we just set a label and its position in the window
		title1=QtGui.QLabel('<b>sxmref_ali3d</b> - set advanced parameters related to CTF and search', self)
		title1.move(self.x1,self.y1)
		
		self.y1 += 30
		#Labels and Line Edits for User Input
		#Just a label
		#title2= QtGui.QLabel('<b>Advanced</b> parameters', self)
		#title2.move(self.x1,self.y1)
		
		#self.y1 += 30
		
		self.savedparmsdict=savedparms
		#Example for User input stack name
		#First create the label and define its position
		maskname= QtGui.QLabel('Mask', self)
		maskname.move(self.x1,self.y1)
		#Now add a line edit and define its position
		self.masknameedit=QtGui.QLineEdit(self)
		self.masknameedit.move(self.x2,self.y1)
		#Adds a default value for the line edit
		self.masknameedit.setText(self.savedparmsdict['maskname'])
		self.masknameedit.setToolTip("File name of mask")
		
		self.mskfile_button = QtGui.QPushButton("Open Mask", self)
		self.mskfile_button.move(self.x3, self.y1-self.yspc)
		#Here we define, that when this button is clicked, it starts subfunction choose_file
		QtCore.QObject.connect(self.mskfile_button, QtCore.SIGNAL("clicked()"), self.choose_mskfile)
		
		self.y1 += 30
		
		focus= QtGui.QLabel('3D mask for focused clustering', self)
		focus.move(self.x1,self.y1)
		#Now add a line edit and define its position
		self.focusedit=QtGui.QLineEdit(self)
		self.focusedit.move(self.x2,self.y1)
		#Adds a default value for the line edit
		self.focusedit.setText(self.savedparmsdict['focus'])
		self.focusedit.setToolTip("3D mask for focused clustering")
		
		self.focusfile_button = QtGui.QPushButton("Open Mask", self)
		self.focusfile_button.move(self.x3, self.y1-self.yspc)
		#Here we define, that when this button is clicked, it starts subfunction choose_file
		QtCore.QObject.connect(self.focusfile_button, QtCore.SIGNAL("clicked()"), self.choose_focusfile)
		
		self.y1 += 30
		
		center= QtGui.QLabel('Center: 0 for none, 1 for cog', self)
		center.move(self.x1,self.y1)
		self.centeredit=QtGui.QLineEdit(self)
		self.centeredit.move(self.x2,self.y1)
		self.centeredit.setText(self.savedparmsdict['center'])
		self.centeredit.setToolTip('0 - if you do not want the volume to be centered, 1 - center the volume using cog (default=-1)')
		
		self.y1 += 30
		
		nasg= QtGui.QLabel('Number of reassignment iterations for \neach angular step', self)
		nasg.move(self.x1,self.y1)
		self.nassignedit=QtGui.QLineEdit(self)
		self.nassignedit.move(self.x2,self.y1)
		self.nassignedit.setText(self.savedparmsdict['nassign'])
		self.nassignedit.setToolTip("number of reassignment iterations performed for each angular step (set to 3)")
		
		self.y1 += 50
		
		nrefi= QtGui.QLabel('Number of alignment iterations for \neach angular step', self)
		nrefi.move(self.x1,self.y1)
		self.nrefineedit=QtGui.QLineEdit(self)
		self.nrefineedit.move(self.x2,self.y1)
		self.nrefineedit.setText(self.savedparmsdict['nrefine'])
		self.nrefineedit.setToolTip("number of alignment iterations performed for each angular step (set to 1)")
		
		self.y1 += 50
		
		ringstep= QtGui.QLabel('Step between rings in rotational correlation', self)
		ringstep.move(self.x1,self.y1)
		self.ringstepedit=QtGui.QLineEdit(self)
		self.ringstepedit.move(self.x2,self.y1)
		self.ringstepedit.setText(self.savedparmsdict['ringstep'])
		self.ringstepedit.setToolTip('step between rings in rotational correlation > 0 (set to 1)')
		
		self.y1 += 30
		
		innerradius= QtGui.QLabel('Inner radius for rotational correlation', self)
		innerradius.move(self.x1,self.y1)
		self.innerradiusedit=QtGui.QLineEdit(self)
		self.innerradiusedit.move(self.x2,self.y1)
		self.innerradiusedit.setText(self.savedparmsdict['innerradius'])
		self.innerradiusedit.setToolTip('Inner radius for rotational correlation > 0 (set to 1) ')		
		
		self.y1 += 30
		
		ctf= QtGui.QLabel('Consider CTF correction during alignment', self)
		ctf.move(self.x1,self.y1)
		self.ctfchkbx = QtGui.QCheckBox("",self)
		self.ctfchkbx.move(self.x2, self.y1)
		self.ctfchkbx.setCheckState(self.savedparmsdict['ctf'])
		self.ctfchkbx.setToolTip('Consider CTF correction during the alignment')
		
		self.y1 += 30
		
		snr= QtGui.QLabel('Signal-to-noise ratio of data', self)
		snr.move(self.x1,self.y1)
		self.snredit=QtGui.QLineEdit(self)
		self.snredit.move(self.x2,self.y1)
		self.snredit.setText(self.savedparmsdict['snr'])
		self.snredit.setToolTip('signal-to-noise ratio of the data (default SNR=1.0)')		
		
		self.y1 += 30
		
		refa= QtGui.QLabel('Method for generating the quasi-uniformly \ndistributed projection directions', self)
		refa.move(self.x1,self.y1)
		#Now add a line edit and define its position
		self.refaedit=QtGui.QLineEdit(self)
		self.refaedit.move(self.x2,self.y1)
		#Adds a default value for the line edit
		self.refaedit.setText(self.savedparmsdict['ref_a'])
		self.refaedit.setToolTip("method for generating the quasi-uniformly distributed projection directions (default S)")
		
		self.y1 += 50
		
		sym= QtGui.QLabel('Symmetry of the structure', self)
		sym.move(self.x1,self.y1)
		self.symedit=QtGui.QLineEdit(self)
		self.symedit.move(self.x2,self.y1)
		self.symedit.setText(self.savedparmsdict['sym'])
		self.symedit.setToolTip('symmetry of the structure')
		
		self.y1 += 30
		
		pad= QtGui.QLabel('Padding size for 3D reconstruction', self)
		pad.move(self.x1,self.y1)
		self.npadedit=QtGui.QLineEdit(self)
		self.npadedit.move(self.x2,self.y1)
		self.npadedit.setText(self.savedparmsdict['npad'])
		self.npadedit.setToolTip('padding size for 3D reconstruction')
		
		self.y1 += 30
		
		an= QtGui.QLabel('Angular neighborhood for local searches', self)
		an.move(self.x1,self.y1)
		self.anedit=QtGui.QLineEdit(self)
		self.anedit.move(self.x2,self.y1)
		self.anedit.setText(self.savedparmsdict['an'])
		self.anedit.setToolTip('Angular neighborhood for local searches ')		

		self.y1 += 30

		fourvar= QtGui.QLabel('Compute and use fourier variance', self)
		fourvar.move(self.x1,self.y1)
		self.fourvarchkbx=QtGui.QCheckBox("",self)
		self.fourvarchkbx.move(self.x2,self.y1)
		self.fourvarchkbx.setCheckState(self.savedparmsdict['fourvar'])
		self.fourvarchkbx.setToolTip('Compute and use fourier variance')
		
		self.y1 += 30
		
		stoprnct= QtGui.QLabel('Minimum percentage of assignment \nchange to stop the program', self)
		stoprnct.move(self.x1,self.y1)
		self.stoprnctedit=QtGui.QLineEdit(self)
		self.stoprnctedit.move(self.x2,self.y1)
		self.stoprnctedit.setText(self.savedparmsdict['stoprnct'])
		self.stoprnctedit.setToolTip('Minimum percentage of assignment change to stop the program')		

		self.y1 += 50
		
		
		dbg= QtGui.QLabel('Debug mode', self)
		dbg.move(self.x1,self.y1)
		self.debugchkbx=QtGui.QCheckBox("",self)
		self.debugchkbx.move(self.x2,self.y1)
		self.debugchkbx.setCheckState(self.savedparmsdict['debug'])
		self.debugchkbx.setToolTip('Debug')


		self.y1 += 30
		
		usrfunc= QtGui.QLabel('Name of the reference preparation function', self)
		usrfunc.move(self.x1,self.y1)
		self.usrfuncedit=QtGui.QLineEdit(self)
		self.usrfuncedit.move(self.x2,self.y1)
		self.usrfuncedit.setText(self.savedparmsdict['usrfunc'])
		self.usrfuncedit.setToolTip('name of the reference preparation function')
		
		self.y1 += 30
				
		usrfuncfile= QtGui.QLabel('Enter name of external file containing user function:', self)
		usrfuncfile.move(self.x1,self.y1)
		
		self.y1 += 20
		
		usrfuncfile= QtGui.QLabel('(Leave blank if file is not external to Sparx)', self)
		usrfuncfile.move(self.x1,self.y1)
		
		self.y1 += 20
		
		self.usrfuncfileedit=QtGui.QLineEdit(self)
		self.usrfuncfileedit.move(self.x2,self.y1)
		self.usrfuncfileedit.setText(self.savedparmsdict['usrfuncfile'])
		self.usrfuncfileedit.setToolTip('name of the external file containing user function')		
		#Function runsxali2d started when  the  RUN_button of the  Poptwodali window is clicked 
			
		#Here we define, that when this button is clicked, it starts subfunction choose_file
		self.usrfile_button = QtGui.QPushButton("Select File", self)
		self.usrfile_button.move(self.x3, self.y1-self.yspc)
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
		
	def choose_focusfile(self):
		#opens a file browser, showing files only in .hdf format
		file_name = QtGui.QFileDialog.getOpenFileName(self, "Open File Containing Mask", "", "(*)")
		#after the user selected a file, we obtain this filename as a Qstring
		a=QtCore.QString(file_name)
		print a
		#we convert this Qstring to a string and send it to line edit classed stackname edit of the Poptwodali window
		self.focusedit.setText(str(a))

class Popupcter(QWidget):
	def __init__(self):
		QWidget.__init__(self)
		
		#######################################################################################
		# class variables
		
		self.cmd = ""
		# populate with default values
		self.savedparmsdict = {'stackname':'','partradius':'-1','xyrange':'1','trans':'1','nriter':'30','nproc':'2',"ringstep":"1","innerradius":"1","ctf":Qt.Unchecked,"snr":"1.0","dst":"90.0","FL":"0.1","FH":"0.3","FF":"0.2","init_iter":"3","main_iter":"3","iter_reali":"1","match_first":"1","max_round":"20","match_second":"5","stab_ali":"5","thld_err":"1.0","indep_run":"4","thld_grp":"10","img_per_grp":"100","generation":"1",'stackname_prectr':'','outdir_prectr':'','mask_prectr':'','search_rng_prectr':'-1','ou_prectr':'-1','maxit_prectr':'100','snr_prectr':'1','fourvar_prectr':Qt.Unchecked,'ctf_prectr':Qt.Unchecked,'oneDx_prectr':Qt.Unchecked,'nproc_prectr':'2'}

		if (options.demo == DEMO_mpibdbctf) or (options.demo == DEMO_mpibdb):
			self.savedparmsdict['stackname'] = 'bdb:data'
			self.savedparmsdict['stab_ali']='2'
			self.savedparmsdict['init_iter']='1'
			self.savedparmsdict['main_iter']='1'
			self.savedparmsdict['match_second']='1'
			self.savedparmsdict['partradius']='30'
			self.savedparmsdict['max_round']='5'
			self.savedparmsdict['img_per_grp']='120'
			self.savedparmsdict['thld_err']='0.75'
			self.savedparmsdict['generation']='1'
			self.savedparmsdict['nproc']='4'
		#######################################################################################
		# Layout parameters
		
		self.y1 = 10 # title and Repopulate button
		self.y2 = self.y1 + 78 # Text boxes for inputting parameters
		self.y4 = self.y2 + 500 # Save Input and Generate command line buttons
		self.y5 = self.y4 + 95 # run button 
		self.yspc = 4
		
		self.x1 = 10 # first column (text box labels)
		self.x2 = self.x1 + 330 # second column (text boxes)
		self.x3 = self.x2+145 # third column (Open .hdf button)
		self.x4 = self.x3+100 # fourth column (Open .bdb button)
		self.x5 = 230 # run button
		#######################################################################################
		
		#Here we just set the window title
		self.setWindowTitle('sxisac')
		#Here we just set a label and its position in the window
		title1=QtGui.QLabel('<b>sxisac</b> - Perform Iterative Stable Alignment and Clustering (ISAC) on a 2-D image stack', self)
		title1.move(self.x1,self.y1)
		self.y1 += 30

		self.repopbtn = QPushButton("Retrieve saved parameters", self)
		self.repopbtn.move(self.x1-5,self.y1)
		#sets an infotip for this Pushbutton
		self.repopbtn.setToolTip('Retrieve saved parameters')
		#when this button is clicked, this action starts the subfunction twodali
		self.connect(self.repopbtn, SIGNAL("clicked()"), self.repoparms_isac)

		#######################################################################################
		#Here we create a Button(file_button with title run open .hdf) and its position in the window
		self.file_button = QtGui.QPushButton("Open .hdf", self)
		self.file_button.move(self.x3, self.y2-self.yspc)
		#Here we define, that when this button is clicked, it starts subfunction choose_file
		QtCore.QObject.connect(self.file_button, QtCore.SIGNAL("clicked()"), self.choose_file)
		#exactly the same as above, but for subfunction choose_file1
		self.file_button1 = QtGui.QPushButton("Open .bdb", self)
		self.file_button1.move(self.x4,self.y2-self.yspc)
		QtCore.QObject.connect(self.file_button1, QtCore.SIGNAL("clicked()"), self.choose_file1)
		
		stackname= QtGui.QLabel('Name of input stack', self)
		stackname.move(self.x1,self.y2)
		self.stacknameedit=QtGui.QLineEdit(self)
		self.stacknameedit.move(self.x2,self.y2)
		self.stacknameedit.setText(self.savedparmsdict['stackname'])
		self.y2 += 30
		
		partradius= QtGui.QLabel('Particle radius', self)
		partradius.move(self.x1,self.y2)
		self.partradiusedit=QtGui.QLineEdit(self)
		self.partradiusedit.move(self.x2,self.y2)
		self.partradiusedit.setText(self.savedparmsdict['partradius'])
		self.partradiusedit.setToolTip('Parameter ou: Outer radius for rotational correlation \nshould be set to particle radius\nif not sure, set to boxsize/2-2 ')		
		self.y2 += 30
		
		match_first= QtGui.QLabel('match_first (number of iterations to run 2-way \nmatching in the first phase)', self)
		match_first.move(self.x1,self.y2)
		self.match_firstedit=QtGui.QLineEdit(self)
		self.match_firstedit.move(self.x2,self.y2)
		self.match_firstedit.setText(self.savedparmsdict['match_first'])
		self.match_firstedit.setToolTip('number of iterations to run 2-way matching in the first phase')
		
		self.y2 += 50
		
		
		match_second= QtGui.QLabel('match_second (number of iterations to run \n2-way or 3-way matching in the second phase)', self)
		match_second.move(self.x1,self.y2)
		self.match_secondedit=QtGui.QLineEdit(self)
		self.match_secondedit.move(self.x2,self.y2)
		self.match_secondedit.setText(self.savedparmsdict['match_second'])
		self.match_secondedit.setToolTip('number of iterations to run 2-way (or 3-way) matching in the second phase')
		
		self.y2 += 50
		
		nriter= QtGui.QLabel('maxit (Number of iterations for reference-free \nalignment)', self)
		nriter.move(self.x1,self.y2)
		self.nriteredit=QtGui.QLineEdit(self)
		self.nriteredit.move(self.x2,self.y2)
		self.nriteredit.setText(self.savedparmsdict['nriter'])
		self.nriteredit.setToolTip('number of iterations for reference-free alignment')
		self.y2 += 50
		
		
		img_per_grp= QtGui.QLabel('img_per_grp (number of images per group in \nthe ideal case-essentially maximum size of \nclass)', self)
		img_per_grp.move(self.x1,self.y2)
		self.img_per_grpedit=QtGui.QLineEdit(self)
		self.img_per_grpedit.move(self.x2,self.y2)
		self.img_per_grpedit.setText(self.savedparmsdict['img_per_grp'])
		self.img_per_grpedit.setToolTip('number of images per group in the ideal case (essentially maximum size of class)')
		
		self.y2 += 70
		
		thld_grp= QtGui.QLabel('thld_grp (the threshold of size of reproducible \nclass (essentially minimum size of class))', self)
		thld_grp.move(self.x1,self.y2)
		self.thld_grpedit=QtGui.QLineEdit(self)
		self.thld_grpedit.move(self.x2,self.y2)
		self.thld_grpedit.setText(self.savedparmsdict['thld_grp'])
		self.thld_grpedit.setToolTip('the threshold of size of reproducible class (essentially minimum size of class)')
		self.y2 += 50
		
		thld_err= QtGui.QLabel('thld_err (the threshold of pixel error when \nchecking stability)', self)
		thld_err.move(self.x1,self.y2)
		self.thld_erredit=QtGui.QLineEdit(self)
		self.thld_erredit.move(self.x2,self.y2)
		self.thld_erredit.setText(self.savedparmsdict['thld_err'])
		self.thld_erredit.setToolTip('the threshold of pixel error when checking stability')
		
		self.y2 += 50
		
		generation= QtGui.QLabel('generation (the n-th approach on the dataset)', self)
		generation.move(self.x1,self.y2)
		self.generationedit=QtGui.QLineEdit(self)
		self.generationedit.move(self.x2,self.y2)
		self.generationedit.setText(self.savedparmsdict['generation'])
		self.generationedit.setToolTip('the n-th approach on the dataset')
		
		self.y2 += 30
		
		
		nproc= QtGui.QLabel('MPI processors (default=True, False is not \ncurrently supported)', self)
		nproc.move(self.x1,self.y2)
		self.nprocedit=QtGui.QLineEdit(self)
		self.nprocedit.move(self.x2,self.y2)
		self.nprocedit.setText(self.savedparmsdict['nproc'])
		self.nprocedit.setToolTip('use MPI version (default=True, currently False is not supported)')

		
		######################################################################################
		
		self.savepbtn = QPushButton("Save Input Parameters", self)
		self.savepbtn.move(self.x1-5, self.y4)
		#sets an infotip for this Pushbutton
		self.savepbtn.setToolTip('Save Input Parameters')
		#when this button is clicked, this action starts the subfunction twodali
		self.connect(self.savepbtn, SIGNAL("clicked()"), self.saveparms)
		self.y4+=30
		
		self.cmdlinebtn = QPushButton("Generate command line from input parameters", self)
		self.cmdlinebtn.move(self.x1-5, self.y4)
		#sets an infotip for this Pushbutton
		self.cmdlinebtn.setToolTip('Generate command line using input parameters')
		#when this button is clicked, this action starts the subfunction twodali
		self.connect(self.cmdlinebtn, SIGNAL("clicked()"), self.gencmdline_isac)
		self.y4+=30

		self.outinfobtn = QPushButton("Output Info", self)
		self.outinfobtn.move(self.x1-5, self.y4)
		#sets an infotip for this Pushbutton
		self.outinfobtn.setToolTip('Output Info')
		#when this button is clicked, this action starts the subfunction twodali
		self.connect(self.outinfobtn, SIGNAL("clicked()"), self.outputinfo_isac)
		
		#######################################################################
		 #Here we create a Button(Run_button with title run sxali2d) and its position in the window
		self.RUN_button = QtGui.QPushButton('Run sxisac', self)
		# make 3D textured push button look
		s = "QPushButton {font: bold; color: #000;border: 1px solid #333;border-radius: 11px;padding: 2px;background: qradialgradient(cx: 0, cy: 0,fx: 0.5, fy:0.5,radius: 1, stop: 0 #fff, stop: 1 #8D0);min-width:90px;margin:5px} QPushButton:pressed {font: bold; color: #000;border: 1px solid #333;border-radius: 11px;padding: 2px;background: qradialgradient(cx: 0, cy: 0,fx: 0.5, fy:0.5,radius: 1, stop: 0 #fff, stop: 1 #084);min-width:90px;margin:5px}"
		
		self.RUN_button.setStyleSheet(s)
		self.RUN_button.move(self.x5, self.y5)
		#Here we define, that when this button is clicked, it starts subfunction runsxali2d
		self.connect(self.RUN_button, SIGNAL("clicked()"), self.runsxisac)
		#Labels and Line Edits for User Input		

	def outputinfo_isac(self):
		QMessageBox.information(self, "isac output",'For each generation of running the program, there are two phases. The first phase is an exploratory phase. In this phase, we set the criteria to be very loose and try to find as much candidate class averages as possible. This phase typically should have 10 to 20 rounds (set by max_round, default = 20). The candidate class averages are stored in class_averages_candidate_generation_n.hdf.\n\n The second phase is where the actual class averages are generated, it typically have 3~9 iterations (set by match_second, default = 5) of matching. The first half of iterations are 2-way matching, the second half of iterations are 3-way matching, and the last iteration is 4-way matching. In the second phase, three files will be generated:\n\n class_averages_generation_n.hdf : class averages generated in this generation \n\n generation_n_accounted.txt : IDs of accounted particles in this generation\n\n generation_n_unaccounted.txt : IDs of unaccounted particles in this generation')
		
	def gencmdline_isac(self,writefile=True):
		#Here we just read in all user inputs in the line edits of the Poptwodali window
		stack = self.stacknameedit.text()
		ou=self.partradiusedit.text()
		
		match_first = self.match_firstedit.text()
		match_second = self.match_secondedit.text()
		maxit=self.nriteredit.text()
		thld_err = self.thld_erredit.text()
		thld_grp = self.thld_grpedit.text()
		img_per_grp = self.img_per_grpedit.text()
		generation = self.generationedit.text()
		
		xr=self.w1.xyrangeedit.text()
		yr=self.w1.xyrangeedit.text()
		ts=self.w1.transedit.text()
		dst = self.w1.dstedit.text()
		FL = self.w1.FLedit.text()
		FH = self.w1.FHedit.text()
		FF = self.w1.FFedit.text()
		init_iter = self.w1.init_iteredit.text()
		main_iter = self.w1.main_iteredit.text()
		iter_reali = self.w1.iter_realiedit.text()
		max_round = self.w1.max_roundedit.text()
		stab_ali = self.w1.stab_aliedit.text()
		indep_run = self.w1.indep_runedit.text()
		ringstep = self.w1.ringstepedit.text()
		inrad = self.w1.innerradiusedit.text()
		CTF=self.w1.ctfchkbx.checkState()
		snr = self.w1.snredit.text()
		
		cmd1 = "sxisac.py "+str(stack)+ " --ir=" + str(inrad)+" --ou="+ str(ou) +" --rs="+str(ringstep)+ " --xr="+str(xr)+ " --yr="+str(yr)+ " --ts="+str(ts)+ " --maxit="+ str(maxit)+  " --snr=" + str(snr)
		
		
		if CTF == Qt.Checked:
				cmd1 = cmd1 + " --CTF"
		
		cmd1 = cmd1 + " --dst=" + str(dst) + " --FL="+str(FL)+" --FH="+str(FH)+" --FF="+str(FF)+" --init_iter="+str(init_iter)+" --main_iter="+str(main_iter)+ " --iter_reali="+str(iter_reali)+" --match_first="+str(match_first)+ " --match_second="+str(match_second)+" --max_round="+str(max_round)+" --stab_ali="+str(stab_ali)+" --thld_err="+str(thld_err)+" --indep_run="+str(indep_run)+" --thld_grp="+str(thld_grp)+" --img_per_grp="+str(img_per_grp)+" --generation="+str(generation)
		
		np = self.nprocedit.text()
		
		self.savedparmsdict['stackname']=str(stack)
		self.savedparmsdict['partradius']=str(ou)
		self.savedparmsdict['xyrange']=str(xr)
		self.savedparmsdict['trans']=str(ts)
		self.savedparmsdict['nriter']=str(maxit)
		self.savedparmsdict['nproc']=str(np)
		self.savedparmsdict['ringstep']=str(ringstep)
		self.savedparmsdict['innerradius']=str(inrad)
		self.savedparmsdict['ctf']=CTF
		self.savedparmsdict['snr']=str(snr)
		self.savedparmsdict['dst']=str(dst)
		self.savedparmsdict['FL']=str(FL)
		self.savedparmsdict['FH']=str(FH)
		self.savedparmsdict['FF']=str(FF)
		self.savedparmsdict['init_iter']=str(init_iter)
		self.savedparmsdict['main_iter']=str(main_iter)
		self.savedparmsdict['iter_reali']=str(iter_reali)
		self.savedparmsdict['match_first']=str(match_first)
		self.savedparmsdict['max_round']=str(max_round)
		self.savedparmsdict['match_second']=str(match_second)
		self.savedparmsdict['stab_ali']=str(stab_ali)
		self.savedparmsdict['thld_err']=str(thld_err)
		self.savedparmsdict['indep_run']=str(indep_run)
		self.savedparmsdict['thld_grp']=str(thld_grp)
		self.savedparmsdict['img_per_grp']=str(img_per_grp)
		self.savedparmsdict['generation']=str(generation)
		
		#self.savedparmsdict = {'stackname':str(stack),'partradius':str(ou),'xyrange':str(xr),'trans':str(ts),'nriter':str(maxit),'nproc':str(np),"ringstep":str(ringstep),"innerradius":str(inrad),"ctf":CTF,"snr":str(snr),"dst":str(dst),"FL":str(FL),"FH":str(FH),"FF":str(FF),"init_iter":str(init_iter),"main_iter":str(main_iter),"iter_reali":str(iter_reali),"match_first":str(match_first),"max_round":str(max_round),"match_second":str(match_second),"stab_ali":str(stab_ali),"thld_err":str(thld_err),"indep_run":str(indep_run),"thld_grp":str(thld_grp),"img_per_grp":str(img_per_grp),"generation":str(generation)}
		
		self.w1.savedparmsdict=self.savedparmsdict
		
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
		
	def runsxisac(self):
		self.gencmdline_isac(writefile=False)
		process = subprocess.Popen(self.cmd,shell=True)
		self.emit(QtCore.SIGNAL("process_started"),process.pid)
		
	def saveparms(self):		
		(fname,stat)= QInputDialog.getText(self,"Save Input Parameters","Enter name of file to save parameters in",QLineEdit.Normal,"")
		if stat:
			import pickle
			output=open(fname,'wb')
			self.gencmdline_isac(writefile=False)
			self.w2.gencmdline_shftali(writefile=False)
			pickle.dump(self.savedparmsdict,output)
			output.close()
		
	def repoparms_isac(self):		
		(fname,stat)= QInputDialog.getText(self,"Retrieve saved parameters","Enter name of file parameters were saved in",QLineEdit.Normal,"")
		if stat:
			import pickle
			pkl = open(fname,'rb')
			self.savedparmsdict = pickle.load(pkl)
			self.partradiusedit.setText(self.savedparmsdict['partradius'])
			self.stacknameedit.setText(self.savedparmsdict['stackname'])		
			self.nriteredit.setText(self.savedparmsdict['nriter'])
			self.nprocedit.setText(self.savedparmsdict['nproc'])
			self.match_firstedit.setText(self.savedparmsdict['match_first'])
			self.match_secondedit.setText(self.savedparmsdict['match_second'])
			self.thld_erredit.setText(self.savedparmsdict['thld_err'])
			self.thld_grpedit.setText(self.savedparmsdict['thld_grp'])
			self.img_per_grpedit.setText(self.savedparmsdict['img_per_grp'])
			self.generationedit.setText(self.savedparmsdict['generation'])
			
			self.w1.indep_runedit.setText(self.savedparmsdict['indep_run'])
			self.w1.stab_aliedit.setText(self.savedparmsdict['stab_ali'])
			self.w1.max_roundedit.setText(self.savedparmsdict['max_round'])
			self.w1.xyrangeedit.setText(self.savedparmsdict['xyrange'])
			self.w1.transedit.setText(self.savedparmsdict['trans'])
			self.w1.ringstepedit.setText(self.savedparmsdict['ringstep'])
			self.w1.innerradiusedit.setText(self.savedparmsdict['innerradius'])
			self.w1.ctfchkbx.setCheckState(self.savedparmsdict['ctf'])
			self.w1.snredit.setText(self.savedparmsdict['snr'])
			self.w1.dstedit.setText(self.savedparmsdict['dst'])
			self.w1.FLedit.setText(self.savedparmsdict['FL'])
			self.w1.FHedit.setText(self.savedparmsdict['FH'])
			self.w1.FFedit.setText(self.savedparmsdict['FF'])
			self.w1.init_iteredit.setText(self.savedparmsdict['init_iter'])				
			self.w1.main_iteredit.setText(self.savedparmsdict['main_iter'])
			self.w1.iter_realiedit.setText(self.savedparmsdict['iter_reali'])

			self.w2.stacknameedit.setText(self.savedparmsdict['stackname_prectr'])
			self.w2.outdiredit.setText(self.savedparmsdict['outdir_prectr'])
			self.w2.maskedit.setText(self.savedparmsdict['mask_prectr'])
			self.w2.search_rngedit.setText(self.savedparmsdict['search_rng_prectr'])
			self.w2.ouedit.setText(self.savedparmsdict['ou_prectr'])
			self.w2.maxitedit.setText(self.savedparmsdict['maxit_prectr'])
			self.w2.snredit.setText(self.savedparmsdict['snr_prectr'])
			self.w2.ctfchkbx.setCheckState(self.savedparmsdict['ctf_prectr'])
			self.w2.fourvarchkbx.setCheckState(self.savedparmsdict['fourvar_prectr'])
			self.w2.oneDxchkbx.setCheckState(self.savedparmsdict['oneDx_prectr'])
			self.w2.nprocedit.setText(self.savedparmsdict['nproc_prectr'])
				
				
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


class Popupisac(QWidget):

	def __init__(self):
		QWidget.__init__(self)
		
		########################################################################################
		# class variables
		self.cmd = ""
		# populate with default values
		# self.savedparmsdict = {'nproc': "1"}
		self.savedparmsdict = {'nproc': "1", 'stackname':'','partradius':'-1','xyrange':'1','trans':'1','nriter':'30','nproc':'2',"ringstep":"1","innerradius":"1","ctf":Qt.Unchecked,"snr":"1.0","dst":"90.0","FL":"0.1","FH":"0.3","FF":"0.2","init_iter":"3","main_iter":"3","iter_reali":"1","match_first":"1","max_round":"20","match_second":"5","stab_ali":"5","thld_err":"1.0","indep_run":"4","thld_grp":"10","img_per_grp":"100","generation":"1",'stackname_prectr':'','outdir_prectr':'','mask_prectr':'','search_rng_prectr':'-1','ou_prectr':'-1','maxit_prectr':'100','snr_prectr':'1','fourvar_prectr':Qt.Unchecked,'ctf_prectr':Qt.Unchecked,'oneDx_prectr':Qt.Unchecked,'nproc_prectr':'2'}
		
		self.y1 = 10
		self.y2 = self.y1 + 98
		self.y4 = self.y2 + 450
		self.y5 = self.y4 + 95
		self.yspc = 4
		
		self.x1 = 10
		self.x2 = self.x1 + 200
		self.x3 = self.x2+145
		self.x4 = self.x3+100
		self.x5 = 230
		########################################################################################
		
		#Here we just set the window title
		self.setWindowTitle('sxisac')
		#Here we just set a label and its position in the window
		title1=QtGui.QLabel('<b>sxisac</b> - Perform Iterative Stable Alignment and Clustering (ISAC) on a 2-D image stack', self)
		title1.move(self.x1,self.y1)
		self.y1 += 50

		self.repopbtn = QPushButton("Retrieve saved parameters", self)
		self.repopbtn.move(self.x1-5,self.y1)
		#sets an infotip for this Pushbutton
		self.repopbtn.setToolTip('Retrieve saved parameters')
		#when this button is clicked, this action starts the subfunction twodali
		self.connect(self.repopbtn, SIGNAL("clicked()"), self.repoparms_isac)
		
		args_list_and_parser =  self.get_args_list_and_parser()
		QLabels = [None]*len(args_list_and_parser)
		self.QLineEditsAndChecks = [None]*len(args_list_and_parser)
		# for option_iterator in range(2, len(parser.option_list)):
		for option_iterator in range(len(args_list_and_parser)):
			if args_list_and_parser[option_iterator].help == None or \
				args_list_and_parser[option_iterator].help == "" or \
				args_list_and_parser[option_iterator].help[0] != "<" or \
				args_list_and_parser[option_iterator].help[-10:] == "(advanced)":
				continue
			label = args_list_and_parser[option_iterator].help[1:].split(">")[0]
			# a = QtGui.QCheckBox("",self)
			QLabels[option_iterator] = QtGui.QLabel(label, self)
			QLabels[option_iterator].move(self.x1,self.y2)

			control_is_an_edit_box = args_list_and_parser[option_iterator].action != "store_true"
			control_is_for_a_parameter = args_list_and_parser[option_iterator].action == "LPlFHy3uNTjucRlk"
			if control_is_an_edit_box:
				self.QLineEditsAndChecks[option_iterator] = QtGui.QLineEdit(self)
				# modify_function = self.QLineEditsAndChecks[option_iterator].setText
				self.QLineEditsAndChecks[option_iterator].setText(str(args_list_and_parser[option_iterator].default))
			else:
				self.QLineEditsAndChecks[option_iterator] = QtGui.QCheckBox("",self)
				# modify_function = self.QLineEditsAndChecks[option_iterator].setCheckState
				modify_function = self.QLineEditsAndChecks[option_iterator].setCheckState(args_list_and_parser[option_iterator].default)

			self.QLineEditsAndChecks[option_iterator].move(self.x2,self.y2 - 7)
			# self.QLineEdits[option_iterator].setText(self.savedparmsdict[parser.option_list[option_iterator].dest])
			self.savedparmsdict[args_list_and_parser[option_iterator].dest] = [option_iterator, str(args_list_and_parser[option_iterator].default), control_is_an_edit_box, control_is_for_a_parameter]
			# modify_function(str(parser.option_list[option_iterator].default))
			self.QLineEditsAndChecks[option_iterator].setToolTip(args_list_and_parser[option_iterator].help)		
	
			self.y2 = self.y2+25


		nproc= QtGui.QLabel('MPI processors', self)
		nproc.move(self.x1,self.y2)
		self.nprocedit=QtGui.QLineEdit(self)
		self.nprocedit.move(self.x2,self.y2)
		self.nprocedit.setText(self.savedparmsdict['nproc'])
		self.nprocedit.setToolTip('The number of processors to use. Default is single processor mode')

		self.y2 =self.y2+25
		

		
		########################################################################################
		
		self.savepbtn = QPushButton("Save Input Parameters", self)
		self.savepbtn.move(self.x1-5,  self.y4)
		#sets an infotip for this Pushbutton
		self.savepbtn.setToolTip('Save Input Parameters')
		#when this button is clicked, this action starts the subfunction twodali
		self.connect(self.savepbtn, SIGNAL("clicked()"), self.saveparms)
		
		self.y4 = self.y4+30
		
		self.cmdlinebtn = QPushButton("Generate command line from input parameters", self)
		self.cmdlinebtn.move(self.x1-5,  self.y4)
		#sets an infotip for this Pushbutton
		self.cmdlinebtn.setToolTip('Generate command line using input parameters')
		#when this button is clicked, this action starts the subfunction twodali
		self.connect(self.cmdlinebtn, SIGNAL("clicked()"), self.gencmdline_isac)
		
		########################################################################################
		
		 #Here we create a Button(Run_button with title run sxali2d) and its position in the window
		self.RUN_button = QtGui.QPushButton('Run sxisac', self)
		# make 3D textured push button look
		s = "QPushButton {font: bold; color: #000;border: 1px solid #333;border-radius: 11px;padding: 2px;background: qradialgradient(cx: 0, cy: 0,fx: 0.5, fy:0.5,radius: 1, stop: 0 #fff, stop: 1 #8D0);min-width:90px;margin:5px} QPushButton:pressed {font: bold; color: #000;border: 1px solid #333;border-radius: 11px;padding: 2px;background: qradialgradient(cx: 0, cy: 0,fx: 0.5, fy:0.5,radius: 1, stop: 0 #fff, stop: 1 #084);min-width:90px;margin:5px}"
		
		self.RUN_button.setStyleSheet(s)
		
		self.RUN_button.move(self.x5,  self.y5)
		#Here we define, that when this button is clicked, it starts subfunction runsxali2d
		self.connect(self.RUN_button, SIGNAL("clicked()"), self.runsxisac)

	def get_args_list_and_parser(self):

		from sxisac import main as sxisac_main
		parser = sxisac_main(["aa", "--return_options"])

		# ['__doc__', '__init__', '__module__', '_add_help_option', '_add_version_option', '_check_conflict', '_create_option_list', '_create_option_mappings', '_get_all_options', '_get_args', '_get_encoding', '_init_parsing_state', '_long_opt', '_match_long_opt', '_populate_option_list', '_process_args', '_process_long_opt', '_process_short_opts', '_share_option_mappings', '_short_opt', 'add_option', 'add_option_group', 'add_options', 'allow_interspersed_args', 'check_values', 'conflict_handler', 'defaults', 'description', 'destroy', 'disable_interspersed_args', 'enable_interspersed_args', 'epilog', 'error', 'exit', 'expand_prog_name', 'format_description', 'format_epilog', 'format_help', 'format_option_help', 'formatter', 'get_default_values', 'get_description', 'get_option', 'get_option_group', 'get_prog_name', 'get_usage', 'get_version', 'has_option', 'largs', 'option_class', 'option_groups', 'option_list', 'parse_args', 'print_help', 'print_usage', 'print_version', 'process_default_values', 'prog', 'rargs', 'remove_option', 'set_conflict_handler', 'set_default', 'set_defaults', 'set_description', 'set_process_default_values', 'set_usage', 'standard_option_list', 'usage', 'values', 'version']
		# >>> a.optionlist
		# Traceback (most recent call last):
		#   File "<stdin>", line 1, in <module>
		# AttributeError: OptionParser instance has no attribute 'optionlist'
		# >>> len(a.option_list)
		# 25
		# >>> a.option_list[1].help
		# 'show this help message and exit'
		# >>> a.option_list[11].help
		# 'maximum number of iterations performed for the finishing up part (set to 50) '
		# >>>
		# >>> dir(a.option_list[11])
		# ['ACTIONS', 'ALWAYS_TYPED_ACTIONS', 'ATTRS', 'CHECK_METHODS', 'CONST_ACTIONS', 'STORE_ACTIONS', 'TYPED_ACTIONS', 'TYPES', 'TYPE_CHECKER', '__doc__', '__init__', '__module__', '__repr__', '__str__', '_check_action', '_check_callback', '_check_choice', '_check_const', '_check_dest', '_check_nargs', '_check_opt_strings', '_check_type', '_long_opts', '_set_attrs', '_set_opt_strings', '_short_opts', 'action', 'callback', 'callback_args', 'callback_kwargs', 'check_value', 'choices', 'const', 'container', 'convert_value', 'default', 'dest', 'get_opt_string', 'help', 'metavar', 'nargs', 'process', 'take_action', 'takes_value', 'type']
		# >>> a.option_list[11].dest
		# 'maxit2'
		
		class ImmitateOptionList:
			default=""
			help=""
			dest = ""
			action = "LPlFHy3uNTjucRlk"

		prog_args = parser.usage.split("--")[0].split(" ")
		# print prog_args
		args_list = []
		for a in prog_args[1:]:
			if a == "": continue
			aa = a.strip()
			b = ImmitateOptionList()
			b.help = "<" + aa + ">"
			b.dest = aa[:]
			if a == "stack":
				b.help = "<Projection stack file>"
			args_list.append(b)
		
		return args_list + parser.option_list[2:]


	def outputinfo_isac(self):
		QMessageBox.information(self, "sxisac output",'outdir is the name of the output folder specified by the user. If it does not exist, the directory will be created. If it does exist, the program will crash and an error message will come up. Please change the name of directory and restart the program. \n\noutdir/angles_000: \n\nThis file contains Eulerian angles fount in trial #000 \n\noutdir/plot_agls_000.hdf: \n\nThis image in the hdf format contains visualization of the distribution of projections found during trial #000 (see also sxplot_projs_distrib) \n\noutdir/structure_000.hdf: \n\nCopy of the stack of 2D projections with Eulerian angles found at trial #000 set in the header. In order to examine the structure, one has to do the 3D reconstructions sxrecons3d_n.py outdir/structure_000.hdf myvol.hdf \n\noutdir/structure.hdf: \n\nFor multiple trials, it is a copy of the stack of 2D projections with Eulerian angles found at the best trial set in the header (this feature is no longer supported).')
		
	def gencmdline_isac(self,writefile=True):
		#Here we just read in all user inputs in the line edits of the Poptwodali window
		cmd1 = "sxisac.py "

		args = " "
		for param in self.get_args_list_and_parser():
			if param.action != "LPlFHy3uNTjucRlk": break
			args += "%s "%self.QLineEditsAndChecks[self.savedparmsdict[param.dest][0]].text()
		
		for key in self.savedparmsdict:
			if type(self.savedparmsdict[key]) != list:
				continue
			if self.savedparmsdict[key][2]: ## check to see if this is not a boolean option
				if not self.savedparmsdict[key][3]: ## check to see if this is not a parameter
					args += "--%s=%s "%(key,self.QLineEditsAndChecks[self.savedparmsdict[key][0]].text())
				else:
					args += "%s "%self.QLineEditsAndChecks[self.savedparmsdict[key][0]].text()
				self.savedparmsdict[key][1] = self.QLineEditsAndChecks[self.savedparmsdict[key][0]].text()
			else:
				if self.QLineEditsAndChecks[self.savedparmsdict[key][0]].checkState() == Qt.Checked:
					args += "--%s "%key
					self.savedparmsdict[key][1] = self.QLineEditsAndChecks[self.savedparmsdict[key][0]].checkState()
		cmd1 = cmd1 + args

		args = " "
		for key in self.w1.savedparmsdict:
			if type(self.w1.savedparmsdict[key]) != list:
				continue
			# print self.savedparmsdict
			if self.w1.savedparmsdict[key][2]:
				val_str = self.w1.QLineEditsAndChecks[self.w1.savedparmsdict[key][0]].text()
				self.w1.savedparmsdict[key][1] = val_str
				if val_str != "":
					args += "--%s=%s "%(key,val_str)
			else:
				if self.w1.QLineEditsAndChecks[self.w1.savedparmsdict[key][0]].checkState() == Qt.Checked:
					args += "--%s "%key
					self.w1.savedparmsdict[key][1] = self.w1.QLineEditsAndChecks[self.w1.savedparmsdict[key][0]].checkState()
		cmd1 = cmd1 + args
		self.savedparmsdict['nproc'] = self.nprocedit.text()
		
		# np = self.QLineEditsAndChecks[int(self.savedparmsdict["nproc"])].text()
		np = self.nprocedit.text()
		
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
		
	def runsxisac(self):
		self.gencmdline_isac(writefile=False)
		outfolder=self.savedparmsdict['foldername']
		if os.path.exists(outfolder):
			print "output folder "+outfolder+" already exists!"
			return
		process = subprocess.Popen(self.cmd,shell=True)
		self.emit(QtCore.SIGNAL("process_started"),process.pid)
		
	def saveparms(self):		
		(fname,stat)= QInputDialog.getText(self,"Save Input Parameters","Enter name of file to save parameters in",QLineEdit.Normal,"")
		if stat:
			import pickle
			output=open(fname,'wb')
			self.gencmdline_isac(writefile=False)
			pickle.dump((self.savedparmsdict, self.w1.savedparmsdict),output)
			output.close()
		

	def repoparms_isac(self):		
		# repopulate with saved parms
		(fname,stat)= QInputDialog.getText(self,"Retrieve saved parameters","Enter name of file parameters were saved in",QLineEdit.Normal,"")
		print (fname,stat)
		
		if stat:
			import pickle
			pkl = open(fname,'rb')
			# self.savedparmsdict = pickle.load(pkl)
			(self.savedparmsdict, self.w1.savedparmsdict) = pickle.load(pkl)
			# print self.savedparmsdict
			# self.stacknameedit.setText(self.savedparmsdict['stackname'])
			# self.foldernameedit.setText(self.savedparmsdict['foldername'])
			for key in self.savedparmsdict:
				if type(self.savedparmsdict[key]) != list:
					continue
				if self.savedparmsdict[key][2]:
					self.QLineEditsAndChecks[self.savedparmsdict[key][0]].setText(self.savedparmsdict[key][1])
				else:
					# print self.savedparmsdict[key]
					self.QLineEditsAndChecks[self.savedparmsdict[key][0]].setChecked(self.savedparmsdict[key][1] == Qt.Checked)
			for key in self.w1.savedparmsdict:
				if type(self.w1.savedparmsdict[key]) != list:
					continue
				if self.w1.savedparmsdict[key][2]:
					self.w1.QLineEditsAndChecks[self.w1.savedparmsdict[key][0]].setText(self.w1.savedparmsdict[key][1])
				else:
					self.w1.QLineEditsAndChecks[self.w1.savedparmsdict[key][0]].setChecked(self.w1.savedparmsdict[key][1] == Qt.Checked)
		pass	
		self.nprocedit.setText(self.savedparmsdict['nproc'])

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

class Popupadvparams_isac(QWidget):
	def __init__(self,savedparms):
		QWidget.__init__(self)
		
		self.savedparmsdict=dict()
		
		########################################################################################
		# layout parameters
		
		self.y1=10
		self.yspc = 4
		
		self.x1 = 20
		self.x2 = self.x1+280
		self.x3 = self.x2+145
		########################################################################################
		
		#Here we just set the window title
		#self.setWindowTitle('sxali3d advanced parameter selection')
		#Here we just set a label and its position in the window
		title1=QtGui.QLabel('<b>sxisac</b> - set advanced parameters', self)
		title1.move(self.x1,self.y1)
		self.y1 = self.y1+25
		

		from sxisac import main as sxisac_main
		parser = sxisac_main(["aa", "--return_options"])
		
		
		QLabels = [None]*len(parser.option_list)
		self.QLineEditsAndChecks = [None]*len(parser.option_list)
		for option_iterator in range(2, len(parser.option_list)):
		# for option_iterator in range(2, 6):
			if parser.option_list[option_iterator].help == None or \
				parser.option_list[option_iterator].help == "" or \
				parser.option_list[option_iterator].help[0] != "<" or \
				parser.option_list[option_iterator].help[-10:] != "(advanced)":
				continue
			label = parser.option_list[option_iterator].help[1:].split(">")[0]
		
			QLabels[option_iterator] = QtGui.QLabel(label, self)
			QLabels[option_iterator].move(self.x1,self.y1)
		
			control_is_an_edit_box = parser.option_list[option_iterator].action != "store_true"
			if control_is_an_edit_box:
				self.QLineEditsAndChecks[option_iterator] = QtGui.QLineEdit(self)
				modify_function = self.QLineEditsAndChecks[option_iterator].setText
			else:
				self.QLineEditsAndChecks[option_iterator] = QtGui.QCheckBox("",self)
				modify_function = self.QLineEditsAndChecks[option_iterator].setCheckState
		
			self.QLineEditsAndChecks[option_iterator].move(self.x2,self.y1)
			# self.QLineEdits[option_iterator].setText(self.savedparmsdict[parser.option_list[option_iterator].dest])
			self.savedparmsdict[parser.option_list[option_iterator].dest] = [option_iterator, str(parser.option_list[option_iterator].default), control_is_an_edit_box]
			modify_function(str(parser.option_list[option_iterator].default))
			self.QLineEditsAndChecks[option_iterator].setToolTip(parser.option_list[option_iterator].help)		
		
			self.y1 = self.y1+25
		


class Popupadvparams_cter_1(QWidget):
	def __init__(self,savedparms):
		QWidget.__init__(self)
		
		self.x1 = 10
		self.x2 = self.x1+380 
		self.x3 = self.x2 + 145 
		
		self.y1 = 10
		self.yspc = 4
		
		#Here we just set the window title
		self.setWindowTitle('sxisac advanced parameter selection')
		#Here we just set a label and its position in the window
		title1=QtGui.QLabel('<b>sxisac</b> - set advanced params', self)
		title1.move(self.x1,self.y1)
		#Labels and Line Edits for User Input
		#Just a label
		self.y1 += 30
		
		title2= QtGui.QLabel('<b>Advanced</b> parameters', self)
		title2.move(self.x1,self.y1)
		
		self.y1 += 30
		
		self.savedparmsdict=savedparms
		
		
		xyrange= QtGui.QLabel('xr (x and y range of translational search)', self)
		xyrange.move(self.x1,self.y1)
		self.xyrangeedit=QtGui.QLineEdit(self)
		self.xyrangeedit.move(self.x2,self.y1)
		self.xyrangeedit.setText(self.savedparmsdict['xyrange'])
		self.xyrangeedit.setToolTip('Range for translational search in x, y direction\nif set to 0 only rotational alignment will be performed')
		self.y1 += 30
		
		trans= QtGui.QLabel('ts (search step of translational search)', self)
		trans.move(self.x1,self.y1)
		self.transedit=QtGui.QLineEdit(self)
		self.transedit.move(self.x2,self.y1)
		self.transedit.setText(self.savedparmsdict['trans'])
		self.transedit.setToolTip('search step of translational search')		
		self.y1 += 30
		
		
		ringstep= QtGui.QLabel('rs (ring step of the resampling to polar coordinates)', self)
		ringstep.move(self.x1,self.y1)
		self.ringstepedit=QtGui.QLineEdit(self)
		self.ringstepedit.move(self.x2,self.y1)
		self.ringstepedit.setText(self.savedparmsdict['ringstep'])
		self.ringstepedit.setToolTip('ring step of the resampling to polar coordinates')
		
		self.y1 += 30
		
		innerradius= QtGui.QLabel('ir (Inner ring of the resampling to polar coordinates)', self)
		innerradius.move(self.x1,self.y1)
		self.innerradiusedit=QtGui.QLineEdit(self)
		self.innerradiusedit.move(self.x2,self.y1)
		self.innerradiusedit.setText(self.savedparmsdict['innerradius'])
		self.innerradiusedit.setToolTip('inner ring of the resampling to polar coordinates ')		
		
		self.y1 += 30
		
		ctf= QtGui.QLabel('Use CTF information (default=False, currently \nTrue is not supported)', self)
		ctf.move(self.x1,self.y1)
		self.ctfchkbx = QtGui.QCheckBox("",self)
		self.ctfchkbx.move(self.x2, self.y1)
		self.ctfchkbx.setCheckState(self.savedparmsdict['ctf'])
		
		self.y1 += 50
		
		snr= QtGui.QLabel('Signal-to-noise ratio (only meaningful when CTF \nis enabled, currently not supported)', self)
		snr.move(self.x1,self.y1)
		self.snredit=QtGui.QLineEdit(self)
		self.snredit.move(self.x2,self.y1)
		self.snredit.setText(self.savedparmsdict['snr'])
		self.snredit.setToolTip('signal-to-noise ratio (only meaningful when CTF is enabled, currently not supported)')		
		
		self.y1 += 50
				
		dst= QtGui.QLabel('dst (discrete angle used in within group alignment)', self)
		dst.move(self.x1,self.y1)
		self.dstedit=QtGui.QLineEdit(self)
		self.dstedit.move(self.x2,self.y1)
		self.dstedit.setText(self.savedparmsdict['dst'])
		self.dstedit.setToolTip('discrete angle used in within group alignment')		
		
		self.y1 += 30
		
		FL= QtGui.QLabel('FL (lowest stopband frequency used in the tangent filter)', self)
		FL.move(self.x1,self.y1)
		self.FLedit=QtGui.QLineEdit(self)
		self.FLedit.move(self.x2,self.y1)
		self.FLedit.setText(self.savedparmsdict['FL'])
		self.FLedit.setToolTip('lowest stopband frequency used in the tangent filter')
		
		self.y1 += 30
		
		FH= QtGui.QLabel('FH (highest stopband frequency used in the tangent filter)', self)
		FH.move(self.x1,self.y1)
		self.FHedit=QtGui.QLineEdit(self)
		self.FHedit.move(self.x2,self.y1)
		self.FHedit.setText(self.savedparmsdict['FH'])
		self.FHedit.setToolTip('highest stopband frequency used in the tangent filter')
		
		self.y1 += 30
		
		FF= QtGui.QLabel('FF (fall-off of the tangent filter)', self)
		FF.move(self.x1,self.y1)
		self.FFedit=QtGui.QLineEdit(self)
		self.FFedit.move(self.x2,self.y1)
		self.FFedit.setText(self.savedparmsdict['FF'])
		self.FFedit.setToolTip('fall-off of the tangent filter')
		
		self.y1 += 30
		
		init_iter= QtGui.QLabel('init_iter (number of iterations of ISAC program in \ninitialization)', self)
		init_iter.move(self.x1,self.y1)
		self.init_iteredit=QtGui.QLineEdit(self)
		self.init_iteredit.move(self.x2,self.y1)
		self.init_iteredit.setText(self.savedparmsdict['init_iter'])
		self.init_iteredit.setToolTip('number of iterations of ISAC program in initialization')
		
		self.y1 += 50
		
		main_iter= QtGui.QLabel('main_iter (number of iterations of ISAC program in main \npart)', self)
		main_iter.move(self.x1,self.y1)
		self.main_iteredit=QtGui.QLineEdit(self)
		self.main_iteredit.move(self.x2,self.y1)
		self.main_iteredit.setText(self.savedparmsdict['main_iter'])
		self.main_iteredit.setToolTip('number of iterations of ISAC program in main part')
		
		self.y1 += 50
		
		iter_reali= QtGui.QLabel('iter_reali (number of iterations in ISAC before checking \nstability)', self)
		iter_reali.move(self.x1,self.y1)
		self.iter_realiedit=QtGui.QLineEdit(self)
		self.iter_realiedit.move(self.x2,self.y1)
		self.iter_realiedit.setText(self.savedparmsdict['iter_reali'])
		self.iter_realiedit.setToolTip('number of iterations in ISAC before checking stability')
		
		self.y1 += 50
		
		
		
		max_round= QtGui.QLabel('max_round (maximum rounds of generating candidate \naverages in the first phase)', self)
		max_round.move(self.x1,self.y1)
		self.max_roundedit=QtGui.QLineEdit(self)
		self.max_roundedit.move(self.x2,self.y1)
		self.max_roundedit.setText(self.savedparmsdict['max_round'])
		self.max_roundedit.setToolTip('maximum rounds of generating candidate averages in the first phase')
		
		self.y1 += 50
		
		
		stab_ali= QtGui.QLabel('stab_ali (number of alignments when checking stability)', self)
		stab_ali.move(self.x1,self.y1)
		self.stab_aliedit=QtGui.QLineEdit(self)
		self.stab_aliedit.move(self.x2,self.y1)
		self.stab_aliedit.setText(self.savedparmsdict['stab_ali'])
		self.stab_aliedit.setToolTip('number of alignments when checking stability')
		
		self.y1 += 30
		
		
		
		indep_run= QtGui.QLabel('indep_run (number of indepentdent runs for reproducibility \n(default=4, currently other values not supported)', self)
		indep_run.move(self.x1,self.y1)
		self.indep_runedit=QtGui.QLineEdit(self)
		self.indep_runedit.move(self.x2,self.y1)
		self.indep_runedit.setText(self.savedparmsdict['indep_run'])
		self.indep_runedit.setToolTip('number of indepentdent runs for reproducibility (default=4, currently other values not supported')

class Popupadvparams_isac_1(QWidget):
	def __init__(self,savedparms):
		QWidget.__init__(self)
		
		self.x1 = 10
		self.x2 = self.x1+380 
		self.x3 = self.x2 + 145 
		
		self.y1 = 10
		self.yspc = 4
		
		#Here we just set the window title
		self.setWindowTitle('sxisac advanced parameter selection')
		#Here we just set a label and its position in the window
		title1=QtGui.QLabel('<b>sxisac</b> - set advanced params', self)
		title1.move(self.x1,self.y1)
		#Labels and Line Edits for User Input
		#Just a label
		self.y1 += 30
		
		title2= QtGui.QLabel('<b>Advanced</b> parameters', self)
		title2.move(self.x1,self.y1)
		
		self.y1 += 30
		
		self.savedparmsdict=savedparms
		
		
		xyrange= QtGui.QLabel('xr (x and y range of translational search)', self)
		xyrange.move(self.x1,self.y1)
		self.xyrangeedit=QtGui.QLineEdit(self)
		self.xyrangeedit.move(self.x2,self.y1)
		self.xyrangeedit.setText(self.savedparmsdict['xyrange'])
		self.xyrangeedit.setToolTip('Range for translational search in x, y direction\nif set to 0 only rotational alignment will be performed')
		self.y1 += 30
		
		trans= QtGui.QLabel('ts (search step of translational search)', self)
		trans.move(self.x1,self.y1)
		self.transedit=QtGui.QLineEdit(self)
		self.transedit.move(self.x2,self.y1)
		self.transedit.setText(self.savedparmsdict['trans'])
		self.transedit.setToolTip('search step of translational search')		
		self.y1 += 30
		
		
		ringstep= QtGui.QLabel('rs (ring step of the resampling to polar coordinates)', self)
		ringstep.move(self.x1,self.y1)
		self.ringstepedit=QtGui.QLineEdit(self)
		self.ringstepedit.move(self.x2,self.y1)
		self.ringstepedit.setText(self.savedparmsdict['ringstep'])
		self.ringstepedit.setToolTip('ring step of the resampling to polar coordinates')
		
		self.y1 += 30
		
		innerradius= QtGui.QLabel('ir (Inner ring of the resampling to polar coordinates)', self)
		innerradius.move(self.x1,self.y1)
		self.innerradiusedit=QtGui.QLineEdit(self)
		self.innerradiusedit.move(self.x2,self.y1)
		self.innerradiusedit.setText(self.savedparmsdict['innerradius'])
		self.innerradiusedit.setToolTip('inner ring of the resampling to polar coordinates ')		
		
		self.y1 += 30
		
		ctf= QtGui.QLabel('Use CTF information (default=False, currently \nTrue is not supported)', self)
		ctf.move(self.x1,self.y1)
		self.ctfchkbx = QtGui.QCheckBox("",self)
		self.ctfchkbx.move(self.x2, self.y1)
		self.ctfchkbx.setCheckState(self.savedparmsdict['ctf'])
		
		self.y1 += 50
		
		snr= QtGui.QLabel('Signal-to-noise ratio (only meaningful when CTF \nis enabled, currently not supported)', self)
		snr.move(self.x1,self.y1)
		self.snredit=QtGui.QLineEdit(self)
		self.snredit.move(self.x2,self.y1)
		self.snredit.setText(self.savedparmsdict['snr'])
		self.snredit.setToolTip('signal-to-noise ratio (only meaningful when CTF is enabled, currently not supported)')		
		
		self.y1 += 50
				
		dst= QtGui.QLabel('dst (discrete angle used in within group alignment)', self)
		dst.move(self.x1,self.y1)
		self.dstedit=QtGui.QLineEdit(self)
		self.dstedit.move(self.x2,self.y1)
		self.dstedit.setText(self.savedparmsdict['dst'])
		self.dstedit.setToolTip('discrete angle used in within group alignment')		
		
		self.y1 += 30
		
		FL= QtGui.QLabel('FL (lowest stopband frequency used in the tangent filter)', self)
		FL.move(self.x1,self.y1)
		self.FLedit=QtGui.QLineEdit(self)
		self.FLedit.move(self.x2,self.y1)
		self.FLedit.setText(self.savedparmsdict['FL'])
		self.FLedit.setToolTip('lowest stopband frequency used in the tangent filter')
		
		self.y1 += 30
		
		FH= QtGui.QLabel('FH (highest stopband frequency used in the tangent filter)', self)
		FH.move(self.x1,self.y1)
		self.FHedit=QtGui.QLineEdit(self)
		self.FHedit.move(self.x2,self.y1)
		self.FHedit.setText(self.savedparmsdict['FH'])
		self.FHedit.setToolTip('highest stopband frequency used in the tangent filter')
		
		self.y1 += 30
		
		FF= QtGui.QLabel('FF (fall-off of the tangent filter)', self)
		FF.move(self.x1,self.y1)
		self.FFedit=QtGui.QLineEdit(self)
		self.FFedit.move(self.x2,self.y1)
		self.FFedit.setText(self.savedparmsdict['FF'])
		self.FFedit.setToolTip('fall-off of the tangent filter')
		
		self.y1 += 30
		
		init_iter= QtGui.QLabel('init_iter (number of iterations of ISAC program in \ninitialization)', self)
		init_iter.move(self.x1,self.y1)
		self.init_iteredit=QtGui.QLineEdit(self)
		self.init_iteredit.move(self.x2,self.y1)
		self.init_iteredit.setText(self.savedparmsdict['init_iter'])
		self.init_iteredit.setToolTip('number of iterations of ISAC program in initialization')
		
		self.y1 += 50
		
		main_iter= QtGui.QLabel('main_iter (number of iterations of ISAC program in main \npart)', self)
		main_iter.move(self.x1,self.y1)
		self.main_iteredit=QtGui.QLineEdit(self)
		self.main_iteredit.move(self.x2,self.y1)
		self.main_iteredit.setText(self.savedparmsdict['main_iter'])
		self.main_iteredit.setToolTip('number of iterations of ISAC program in main part')
		
		self.y1 += 50
		
		iter_reali= QtGui.QLabel('iter_reali (number of iterations in ISAC before checking \nstability)', self)
		iter_reali.move(self.x1,self.y1)
		self.iter_realiedit=QtGui.QLineEdit(self)
		self.iter_realiedit.move(self.x2,self.y1)
		self.iter_realiedit.setText(self.savedparmsdict['iter_reali'])
		self.iter_realiedit.setToolTip('number of iterations in ISAC before checking stability')
		
		self.y1 += 50
		
		
		
		max_round= QtGui.QLabel('max_round (maximum rounds of generating candidate \naverages in the first phase)', self)
		max_round.move(self.x1,self.y1)
		self.max_roundedit=QtGui.QLineEdit(self)
		self.max_roundedit.move(self.x2,self.y1)
		self.max_roundedit.setText(self.savedparmsdict['max_round'])
		self.max_roundedit.setToolTip('maximum rounds of generating candidate averages in the first phase')
		
		self.y1 += 50
		
		
		stab_ali= QtGui.QLabel('stab_ali (number of alignments when checking stability)', self)
		stab_ali.move(self.x1,self.y1)
		self.stab_aliedit=QtGui.QLineEdit(self)
		self.stab_aliedit.move(self.x2,self.y1)
		self.stab_aliedit.setText(self.savedparmsdict['stab_ali'])
		self.stab_aliedit.setToolTip('number of alignments when checking stability')
		
		self.y1 += 30
		
		
		
		indep_run= QtGui.QLabel('indep_run (number of indepentdent runs for reproducibility \n(default=4, currently other values not supported)', self)
		indep_run.move(self.x1,self.y1)
		self.indep_runedit=QtGui.QLineEdit(self)
		self.indep_runedit.move(self.x2,self.y1)
		self.indep_runedit.setText(self.savedparmsdict['indep_run'])
		self.indep_runedit.setToolTip('number of indepentdent runs for reproducibility (default=4, currently other values not supported')
		


"""
class Popupcenter(QWidget):
	def __init__(self,winmain,intro_string):
		QWidget.__init__(self)
		
		#######################################################################################
		# class variables
		
		self.cmd = ""
		self.winmain=winmain
		self.intro_string = intro_string
		#######################################################################################
		# Layout parameters
		
		self.y1 = 10 # title and Repopulate button
		self.y2 = self.y1 + 150 # Text boxes for inputting parameters
		self.y3 = self.y2 + 510 # activate images button and set xform.align2d button
		self.y4 = self.y3 + 20 # run button 
		self.yspc = 4
		
		self.x1 = 10 # first column (text box labels)
		self.x2 = self.x1 + 350 # second column (text boxes)
		self.x3 = self.x2+145 # third column (Open .hdf button)
		self.x4 = self.x3+100 # fourth column (Open .bdb button)
		self.x5 = 230 # run button
		#######################################################################################
		
		#Here we just set the window title
		self.setWindowTitle('sxshiftali')
		#Here we just set a label and its position in the window
		title1=QtGui.QLabel(self.intro_string, self)
		title1.move(self.x1,self.y1)
		#self.y1 += 30

		#######################################################################################
		#Here we create a Button(file_button with title run open .hdf) and its position in the window
		self.file_button = QtGui.QPushButton("Open .hdf", self)
		self.file_button.move(self.x3, self.y2-self.yspc)
		#Here we define, that when this button is clicked, it starts subfunction choose_file
		QtCore.QObject.connect(self.file_button, QtCore.SIGNAL("clicked()"), self.choose_file)
		#exactly the same as above, but for subfunction choose_file1
		self.file_button1 = QtGui.QPushButton("Open .bdb", self)
		self.file_button1.move(self.x4,self.y2-self.yspc)
		QtCore.QObject.connect(self.file_button1, QtCore.SIGNAL("clicked()"), self.choose_file1)
		
		stackname= QtGui.QLabel('Name of input stack', self)
		stackname.move(self.x1,self.y2)
		self.stacknameedit=QtGui.QLineEdit(self)
		self.stacknameedit.move(self.x2,self.y2)
		self.stacknameedit.setText(self.winmain.savedparmsdict['stackname_prectr'])
		
		self.y2 += 30
		
		outdir= QtGui.QLabel('Output folder', self)
		outdir.move(self.x1,self.y2)
		self.outdiredit=QtGui.QLineEdit(self)
		self.outdiredit.move(self.x2,self.y2)
		self.outdiredit.setText(self.winmain.savedparmsdict['outdir_prectr'])		
		self.y2 += 30
		
		maskname= QtGui.QLabel('Mask', self)
		maskname.move(self.x1,self.y2)
		#Now add a line edit and define its position
		self.maskedit=QtGui.QLineEdit(self)
		self.maskedit.move(self.x2,self.y2)
		#Adds a default value for the line edit
		self.maskedit.setText(self.winmain.savedparmsdict['mask_prectr'])
		self.maskedit.setToolTip("Default is a circle mask with radius equal to the particle radius")
		
		self.mskfile_button = QtGui.QPushButton("Open File", self)
		self.mskfile_button.move(self.x3, self.y2-self.yspc)
		#Here we define, that when this button is clicked, it starts subfunction choose_file
		QtCore.QObject.connect(self.mskfile_button, QtCore.SIGNAL("clicked()"), self.choose_mskfile)
		
		self.y2 += 30
		
		searchrange= QtGui.QLabel('search_rng (Used to compute the dimension of a \nnwx by nwx section of the 2D ccf which is \nwindowed out for peak search: nwx=2*search_rng+1 (nwx=nx if search_rng is -1))', self)
		searchrange.move(self.x1,self.y2)
		self.search_rngedit=QtGui.QLineEdit(self)
		self.search_rngedit.move(self.x2,self.y2)
		self.search_rngedit.setText(self.winmain.savedparmsdict['search_rng_prectr'])
		self.search_rngedit.setToolTip('Used to compute the dimension of a \nnwx by nwx section of the 2D ccf which is \nwindowed out for peak search: \nnwx=2*search_rng+1 (nwx=nx if search_rng is -1))')		
		self.y2 += 70
		
		partradius= QtGui.QLabel('ou (radius of the particle - used for constructing the \ndefault mask. If ou is -1, then the mask is a circle \nwith radius nx/2 - 2.)', self)
		partradius.move(self.x1,self.y2)
		self.ouedit=QtGui.QLineEdit(self)
		self.ouedit.move(self.x2,self.y2)
		self.ouedit.setText(self.winmain.savedparmsdict['ou_prectr'])
		self.ouedit.setToolTip('Parameter ou: Outer radius for rotational correlation \nshould be set to particle radius\nif not sure, set to boxsize/2-2 ')		
		self.y2 += 70
		
		
		maxit= QtGui.QLabel('maxit (maximum number of iterations program will \nperform)', self)
		maxit.move(self.x1,self.y2)
		self.maxitedit=QtGui.QLineEdit(self)
		self.maxitedit.move(self.x2,self.y2)
		self.maxitedit.setText(self.winmain.savedparmsdict['maxit_prectr'])
		self.maxitedit.setToolTip('')
		self.y2 += 50
		
		
		ctf= QtGui.QLabel('CTF (use CTF correction during centering)', self)
		ctf.move(self.x1,self.y2)
		self.ctfchkbx = QtGui.QCheckBox("",self)
		self.ctfchkbx.move(self.x2, self.y2)
		self.ctfchkbx.setCheckState(self.winmain.savedparmsdict['ctf_prectr'])
		
		self.y2 += 30
		
		snr= QtGui.QLabel('SNR (signal-to-noise ratio of the data)', self)
		snr.move(self.x1,self.y2)
		self.snredit=QtGui.QLineEdit(self)
		self.snredit.move(self.x2,self.y2)
		self.snredit.setText(self.winmain.savedparmsdict['snr_prectr'])
		self.snredit.setToolTip('signal-to-noise ratio of the data (default SNR=1.0)')		
		
		self.y2+= 30
				
		fourvar= QtGui.QLabel('Compute fourier variance', self)
		fourvar.move(self.x1,self.y2)
		self.fourvarchkbx=QtGui.QCheckBox("",self)
		self.fourvarchkbx.move(self.x2,self.y2)
		self.fourvarchkbx.setCheckState(self.winmain.savedparmsdict['fourvar_prectr'])
		self.fourvarchkbx.setToolTip('use Fourier variance to weight the reference (recommended, default False)')
		
		self.y2 += 30
		
		oneDx= QtGui.QLabel('oneDx (Window out central line of 2D cross \ncorrelation for peak search)', self)
		oneDx.move(self.x1,self.y2)
		self.oneDxchkbx = QtGui.QCheckBox("",self)
		self.oneDxchkbx.move(self.x2, self.y2)
		self.oneDxchkbx.setCheckState(self.winmain.savedparmsdict['oneDx_prectr'])
		
		self.y2+= 50
		
		nproc= QtGui.QLabel('MPI Processors (currently only works for n>1 \nprocessors)', self)
		nproc.move(self.x1,self.y2)
		self.nprocedit=QtGui.QLineEdit(self)
		self.nprocedit.move(self.x2,self.y2)
		self.nprocedit.setText(self.winmain.savedparmsdict['nproc_prectr'])
		self.nprocedit.setToolTip('')
		
		self.y2 += 50
		
		self.cmdlinebtn = QPushButton("Generate command line from input parameters", self)
		self.cmdlinebtn.move(self.x1-5, self.y3)
		#sets an infotip for this Pushbutton
		self.cmdlinebtn.setToolTip('Generate command line using input parameters')
		#when this button is clicked, this action starts the subfunction twodali
		self.connect(self.cmdlinebtn, SIGNAL("clicked()"), self.gencmdline_shftali)
		self.y4+=30

		 #Here we create a Button(Run_button with title run sxali2d) and its position in the window
		self.RUN_button = QtGui.QPushButton('Run centering', self)
		# make 3D textured push button look
		s = "QPushButton {font: bold; color: #000;border: 1px solid #333;border-radius: 11px;padding: 2px;background: qradialgradient(cx: 0, cy: 0,fx: 0.5, fy:0.5,radius: 1, stop: 0 #fff, stop: 1 #8D0);min-width:90px;margin:5px} QPushButton:pressed {font: bold; color: #000;border: 1px solid #333;border-radius: 11px;padding: 2px;background: qradialgradient(cx: 0, cy: 0,fx: 0.5, fy:0.5,radius: 1, stop: 0 #fff, stop: 1 #084);min-width:90px;margin:5px}"
		
		self.RUN_button.setStyleSheet(s)
		self.RUN_button.move(self.x5, self.y4)
		#Here we define, that when this button is clicked, it starts subfunction runsxali2d
		self.connect(self.RUN_button, SIGNAL("clicked()"), self.runsxshftali)
		#Labels and Line Edits for User Input		

	def gencmdline_shftali(self,writefile=True):
		#Here we just read in all user inputs in the line edits of the Poptwodali window
		stackname_prectr = self.stacknameedit.text()
		outdir_prectr = self.outdiredit.text()
		mask_prectr = self.maskedit.text()
		search_rng_prectr = self.search_rngedit.text()
		ou_prectr=self.ouedit.text()
		maxit_prectr = self.maxitedit.text()
		ctf_prectr=self.ctfchkbx.checkState()
		snr_prectr = self.snredit.text()
		fourvar_prectr=self.fourvarchkbx.checkState()
		oneDx_prectr=self.oneDxchkbx.checkState()
		nproc_prectr = self.nprocedit.text()
		
		
		cmd1 = "sxshiftali.py "+str(stackname_prectr)+" "+str(outdir_prectr)
		if len(str(mask_prectr))>0:
				print 'has mask'
				cmd1 = cmd1+" "+str(mask_prectr)
				
		cmd1 = cmd1+" --search_rng=" + str(search_rng_prectr)+" --ou="+ str(ou_prectr) +" --maxit="+ str(maxit_prectr)+  " --snr=" + str(snr_prectr)
		
		
		if ctf_prectr == Qt.Checked:
				cmd1 = cmd1 + " --CTF"
		if oneDx_prectr == Qt.Checked:
				cmd1 = cmd1 + " --oneDx"
				
		nproc_prectr = self.nprocedit.text()
		
		(self.winmain.savedparmsdict)['stackname_prectr']=str(stackname_prectr)
		(self.winmain.savedparmsdict)['outdir_prectr']=str(outdir_prectr)
		(self.winmain.savedparmsdict)['mask_prectr']=str(mask_prectr)
		(self.winmain.savedparmsdict)['search_rng_prectr']=str(search_rng_prectr)
		(self.winmain.savedparmsdict)['ou_prectr']=str(ou_prectr)
		(self.winmain.savedparmsdict)['maxit_prectr']=str(maxit_prectr)
		(self.winmain.savedparmsdict)['ctf_prectr']=ctf_prectr
		(self.winmain.savedparmsdict)['snr_prectr']=str(snr_prectr)
		(self.winmain.savedparmsdict)['fourvar_prectr']=fourvar_prectr
		(self.winmain.savedparmsdict)['oneDx_prectr']=oneDx_prectr
		(self.winmain.savedparmsdict)['nproc_prectr']=str(nproc_prectr)
		
		if int(str(nproc_prectr)) > 1:
				cmd1="mpirun -np "+ str(nproc_prectr) + " "+ cmd1+" --MPI" 
		
		if writefile:		
				(fname,stat)= QInputDialog.getText(self,"Generate Command Line","Enter name of file to save command line in",QLineEdit.Normal,"")
				if stat:
						f = open(fname,'a')
						f.write(cmd1)
						f.write('\n')
						f.close()
		
		print cmd1
		self.cmd = cmd1
		
	def runsxshftali(self):
		self.gencmdline_shftali(writefile=False)
		process = subprocess.Popen(self.cmd,shell=True)
		self.emit(QtCore.SIGNAL("process_started"),process.pid)
		
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

	def choose_mskfile(self):
		#opens a file browser, showing files only in .hdf format
		file_name = QtGui.QFileDialog.getOpenFileName(self, "Open File Containing Mask", "", "(*)")
		#after the user selected a file, we obtain this filename as a Qstring
		a=QtCore.QString(file_name)
		#we convert this Qstring to a string and send it to line edit classed stackname edit of the Poptwodali window
		self.maskedit.setText(str(a))		
"""


class Popuplocalthreedali(QWidget):
	def __init__(self):
		QWidget.__init__(self)
		
		########################################################################################
		# class variables
		self.cmd = ""
		# populate with default values
		self.savedparmsdict = {'stackname':'','foldername':'','partradius':'-1', 'delta':'2','trans':'0.25','nriter':'10','nproc':'1','maskname':'','center':'-1',"ctf":Qt.Unchecked,"snr":"1.0","fourvar":Qt.Unchecked,"usrfunc":"ref_ali3d","usrfuncfile":"","sym":'c1',"npad":'4',"debug":False,"chunk":'1.0'}
		
		if options.demo == DEMO_mpibdbctf:
			self.savedparmsdict['stackname']='bdb:data'
			self.savedparmsdict['foldername']='ali3d_e'
			self.savedparmsdict['partradius']='30'
			self.savedparmsdict['delta']='2'
			self.savedparmsdict['trans']='2'
			self.savedparmsdict['nriter']='1'
			self.savedparmsdict['chunk']='0.25'
			self.savedparmsdict['usrfunc']='reference3'
			self.savedparmsdict['center']='0'
			self.savedparmsdict['nproc']='4'
			self.savedparmsdict['ctf']=Qt.Checked
			self.savedparmsdict['debug']=Qt.Checked
				
		if options.demo == DEMO_mpibdb:
			self.savedparmsdict['stackname']='bdb:data'
			self.savedparmsdict['foldername']='ali3d_e'
			self.savedparmsdict['partradius']='30'
			self.savedparmsdict['delta']='2'
			self.savedparmsdict['nriter']='1'
			self.savedparmsdict['chunk']='0.25'
			self.savedparmsdict['usrfunc']='reference4'
			self.savedparmsdict['nproc']='4'
			self.savedparmsdict['debug']=Qt.Checked
		########################################################################################
		# layout parameters
		
		self.y1 = 10
		self.y2 = self.y1 + 98
		self.y3 = self.y2 + 230
		self.y4 = self.y3 + 80
		self.y5 = self.y4 + 95
		self.yspc = 4
		
		self.x1 = 10
		self.x2 = self.x1 + 150
		self.x3 = self.x2+145
		self.x4 = self.x3+100
		self.x5 = 230
		########################################################################################
		
		#Here we just set the window title
		self.setWindowTitle('sxlocal_ali3d')
		#Here we just set a label and its position in the window
		title1=QtGui.QLabel('<b>sxlocal_ali3d</b> - Perform local refinement of 3-D projection alignment of image <br>series using highy accurate gridding method.', self)
		title1.move(self.x1,self.y1)
		self.y1 += 50

		self.repopbtn = QPushButton("Retrieve saved parameters", self)
		self.repopbtn.move(self.x1-5,self.y1)
		#sets an infotip for this Pushbutton
		self.repopbtn.setToolTip('Retrieve saved parameters')
		#when this button is clicked, this action starts the subfunction twodali
		self.connect(self.repopbtn, SIGNAL("clicked()"), self.repoparms_localali3d)
		
		########################################################################################
		
		stackname= QtGui.QLabel('Name of input stack', self)
		stackname.move(self.x1,self.y2)
		#Now add a line edit and define its position
		self.stacknameedit=QtGui.QLineEdit(self)
		self.stacknameedit.move(self.x2,self.y2)
		#Adds a default value for the line edit
		self.stacknameedit.setText(self.savedparmsdict['stackname'])
		
		#Here we create a Button(file_button with title run open .hdf) and its position in the window
		self.file_button = QtGui.QPushButton("Open .hdf", self)
		self.file_button.move(self.x3, self.y2-self.yspc)
		#Here we define, that when this button is clicked, it starts subfunction choose_file
		QtCore.QObject.connect(self.file_button, QtCore.SIGNAL("clicked()"), self.choose_file)
		#exactly the same as above, but for subfunction choose_file1
		self.file_button1 = QtGui.QPushButton("Open .bdb", self)
		self.file_button1.move(self.x4,self.y2-self.yspc)
		QtCore.QObject.connect(self.file_button1, QtCore.SIGNAL("clicked()"), self.choose_file1)
		
		self.y2 += 30
		
		#The same as above, but many line edits include Infotips
		foldername= QtGui.QLabel('Output folder', self)
		foldername.move(self.x1,self.y2)
		self.foldernameedit=QtGui.QLineEdit(self)
		self.foldernameedit.move(self.x2,self.y2)
		self.foldernameedit.setText(self.savedparmsdict['foldername'])		
		
		self.outinfobtn = QPushButton("Output Info", self)
		self.outinfobtn.move(self.x3,  self.y2-self.yspc)
		#sets an infotip for this Pushbutton
		self.outinfobtn.setToolTip('Output Info')
		#when this button is clicked, this action starts the subfunction twodali
		self.connect(self.outinfobtn, SIGNAL("clicked()"), self.outputinfo_localali3d)
		
		self.y2 = self.y2+30
		
		partradius= QtGui.QLabel('Particle radius', self)
		partradius.move(self.x1,self.y2)
		self.partradiusedit=QtGui.QLineEdit(self)
		self.partradiusedit.move(self.x2,self.y2)
		self.partradiusedit.setText(self.savedparmsdict['partradius'])
		self.partradiusedit.setToolTip('outer radius of a circular mask that should encompass the particle< int(nx/2)-1 (set to int(nx/2)-1)')		

		self.y2 = self.y2+30
		
		trans= QtGui.QLabel('Shift bracket (ts)', self)
		trans.move(self.x1,self.y2)
		self.transedit=QtGui.QLineEdit(self)
		self.transedit.move(self.x2,self.y2)
		self.transedit.setText(self.savedparmsdict['trans'])
		self.transedit.setToolTip('shift bracket') 

		self.y2 = self.y2+30
		
		delta= QtGui.QLabel('Angular bracket (delta)', self)
		delta.move(self.x1,self.y2)
		self.deltaedit=QtGui.QLineEdit(self)
		self.deltaedit.move(self.x2,self.y2)
		self.deltaedit.setText(self.savedparmsdict['delta'])
		self.deltaedit.setToolTip('Angular bracket')				
		
		self.y2 =self.y2+30

		nriter= QtGui.QLabel('Number of iterations', self)
		nriter.move(self.x1,self.y2)
		self.nriteredit=QtGui.QLineEdit(self)
		self.nriteredit.move(self.x2,self.y2)
		self.nriteredit.setText(self.savedparmsdict['nriter'])
		self.nriteredit.setToolTip('Maximum number of iterations (set to 10)')
		
		self.y2 =self.y2+30
		
		nproc= QtGui.QLabel('MPI processors', self)
		nproc.move(self.x1,self.y2)
		self.nprocedit=QtGui.QLineEdit(self)
		self.nprocedit.move(self.x2,self.y2)
		self.nprocedit.setText(self.savedparmsdict['nproc'])
		self.nprocedit.setToolTip('The number of processors to use. Default is single processor mode')
		
		########################################################################################
				
		header=QtGui.QLabel('Attributes xform.projection and active parameters must be set in the input stack', self)
		header.move(self.x1,self.y3)
		self.y3 = self.y3 + 30
		
		#not linked to a function yet
		self.activeheader_button = QtGui.QPushButton("activate all images", self)
		self.activeheader_button.move(self.x1-5, self.y3)
		self.connect(self.activeheader_button, SIGNAL("clicked()"), self.setactiveheader)
		
		self.projheader_button = QtGui.QPushButton("set xform.projection", self)
		self.projheader_button.move(self.x1-5+180, self.y3)
		self.connect(self.projheader_button, SIGNAL("clicked()"), self.setprojheader)
		
		########################################################################################
		
		self.savepbtn = QPushButton("Save Input Parameters", self)
		self.savepbtn.move(self.x1-5,  self.y4)
		#sets an infotip for this Pushbutton
		self.savepbtn.setToolTip('Save Input Parameters')
		#when this button is clicked, this action starts the subfunction twodali
		self.connect(self.savepbtn, SIGNAL("clicked()"), self.saveparms)
		
		self.y4 = self.y4+30
		
		self.cmdlinebtn = QPushButton("Generate command line from input parameters", self)
		self.cmdlinebtn.move(self.x1-5,  self.y4)
		#sets an infotip for this Pushbutton
		self.cmdlinebtn.setToolTip('Generate command line using input parameters')
		#when this button is clicked, this action starts the subfunction twodali
		self.connect(self.cmdlinebtn, SIGNAL("clicked()"), self.gencmdline_localali3d)
		
		########################################################################################
		
		#Here we create a Button(Run_button with title run sxali2d) and its position in the window
		self.RUN_button = QtGui.QPushButton('Run sxlocal_ali3d', self)
		# make 3D textured push button look
		s = "QPushButton {font: bold; color: #000;border: 1px solid #333;border-radius: 11px;padding: 2px;background: qradialgradient(cx: 0, cy: 0,fx: 0.5, fy:0.5,radius: 1, stop: 0 #fff, stop: 1 #8D0);min-width:90px;margin:5px} QPushButton:pressed {font: bold; color: #000;border: 1px solid #333;border-radius: 11px;padding: 2px;background: qradialgradient(cx: 0, cy: 0,fx: 0.5, fy:0.5,radius: 1, stop: 0 #fff, stop: 1 #084);min-width:90px;margin:5px}"
		
		self.RUN_button.setStyleSheet(s)
		
		self.RUN_button.move(self.x5,  self.y5)
		#Here we define, that when this button is clicked, it starts subfunction runsxali2d
		self.connect(self.RUN_button, SIGNAL("clicked()"), self.runsxlocal_ali3d)
		#Labels and Line Edits for User Input

		
	#Function runsxali2d started when  the  RUN_button of the  Poptwodali window is clicked 

	def outputinfo_localali3d(self):
		QMessageBox.information(self, "local_ali3d output",'Output directory is the directory name into which the output files will be written. If it does not exist, the directory will be created. If it does exist, the program will crash and an error message will come up; in this case, please change the name of directory and restart the program . The program will write two kind of files: the average of the aligned image series (aqe***.spi) and the Fourier Resolution Criterion curve (dre***). Both files are numbered by the iteration number.\n\n The alignment parameters of each 2D projection (Euler angles and two in-plane translations) are stored in the headers of input files as xform.proj.')
		
	def gencmdline_localali3d(self,writefile=True):
		#Here we just read in all user inputs in the line edits of the Poptwodali window
		stack = self.stacknameedit.text()
		print "stack defined="+ stack
		output=self.foldernameedit.text()
		print "output folder="+ output
		ou=self.partradiusedit.text()
		print "Particle radius="+ ou
		delta=self.deltaedit.text()
		print "delta=" +delta
		ts=self.transedit.text()
		print "ts=" +ts
		maxit=self.nriteredit.text()
		print "maxit="+maxit
		
		cmd1 = "sxlocal_ali3d.py "+str(stack) +" "+ str(output)
		
		args = " --ou="+ str(ou)+ " --ts='"+str(ts)+"'"+ " --delta='"+str(delta)+"'"+" --maxit="+ str(maxit) 
		
		
		mask = self.w1.masknameedit.text()
		if len(str(mask))> 1:
			cmd1 = cmd1+" "+str(mask) 
		cmd1 = cmd1 + args
		
		
		ctr=self.w1.centeredit.text()
		CTF=self.w1.ctfchkbx.checkState()
		snr = self.w1.snredit.text()
		userf = self.w1.usrfuncedit.text()
		userfile = self.w1.usrfuncfileedit.text()
		sym = self.w1.symedit.text()
		npad = self.w1.npadedit.text()
				
		fourvar = self.w1.fourvarchkbx.checkState()
		debug = self.w1.debugchkbx.checkState()
		chunk = self.w1.chunkedit.text()
				
		cmd1 = cmd1+" --center=" +str(ctr)+ " --snr=" + str(snr)+ " --sym=" + str(sym)+ " --npad=" + str(npad)+ " --chunk=" + str(chunk)
		if CTF == Qt.Checked:
			cmd1 = cmd1 + " --CTF"
		if fourvar == Qt.Checked:
			cmd1 = cmd1 + " --fourvar"		
		if debug == Qt.Checked:
			cmd1 = cmd1 + " --debug"
								
		if len(userfile) < 1:
			cmd1 = cmd1 + " --function="+str(userf)
		else:
			userfile = str(userfile)
			# break it down into file name and directory path
			rind = userfile.rfind('/')
			
			if rind == -1:
				userfile = os.path.abspath(userfile)
				rind = userfile.rfind('/')
				
			fname = userfile[rind+1:]
			fname, ext = os.path.splitext(fname)
			fdir = userfile[0:rind]
			cmd1 = cmd1 + " --function=\"[" +fdir+","+fname+","+str(userf)+"]\""
				
		np = self.nprocedit.text()
		
		self.savedparmsdict = {'stackname':str(stack),'foldername':str(output),'partradius':str(ou),'delta':str(delta),'trans':str(ts),'nriter':str(maxit),'nproc':str(np),'maskname':str(mask),'center':str(ctr),"ctf":CTF,"snr":str(snr),"fourvar":fourvar,"usrfunc":str(userf), "usrfuncfile":str(userfile),"sym":str(sym),"npad":str(npad),"debug":debug,"chunk":str(chunk)}

		self.w1.savedparmsdict=self.savedparmsdict
		self.w1.savedparmsdict=self.savedparmsdict
		
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
		
	def runsxlocal_ali3d(self):
		self.gencmdline_localali3d(writefile=False)
		outfolder=self.savedparmsdict['foldername']
		if os.path.exists(outfolder):
			print "output folder "+outfolder+" already exists!"
			return
		process = subprocess.Popen(self.cmd,shell=True)
		self.emit(QtCore.SIGNAL("process_started"),process.pid)
		
	def saveparms(self):		
		(fname,stat)= QInputDialog.getText(self,"Save Input Parameters","Enter name of file to save parameters in",QLineEdit.Normal,"")
		if stat:
			import pickle
			output=open(fname,'wb')
			self.gencmdline_localali3d(writefile=False)
			pickle.dump(self.savedparmsdict,output)
			output.close()
		
	def repoparms_localali3d(self):		
		# repopulate with saved parms
		(fname,stat)= QInputDialog.getText(self,"Retrieve saved parameters","Enter name of file parameters were saved in",QLineEdit.Normal,"")
		if stat:
			import pickle
			pkl = open(fname,'rb')
			self.savedparmsdict = pickle.load(pkl)
			print self.savedparmsdict
			self.partradiusedit.setText(self.savedparmsdict['partradius'])
			self.stacknameedit.setText(self.savedparmsdict['stackname'])
			self.foldernameedit.setText(self.savedparmsdict['foldername'])		
			self.deltaedit.setText(self.savedparmsdict['delta'])
			self.transedit.setText(self.savedparmsdict['trans'])
			self.nriteredit.setText(self.savedparmsdict['nriter'])
			self.nprocedit.setText(self.savedparmsdict['nproc'])
	
			self.w1.masknameedit.setText(self.savedparmsdict['maskname'])
			self.w1.centeredit.setText(self.savedparmsdict['center'])
			self.w1.ctfchkbx.setCheckState(self.savedparmsdict['ctf'])
			self.w1.snredit.setText(self.savedparmsdict['snr'])
			self.w1.usrfuncedit.setText(self.savedparmsdict['usrfunc'])
			self.w1.usrfuncfileedit.setText(self.savedparmsdict['usrfuncfile'])
			self.w1.symedit.setText(self.savedparmsdict['sym'])
			self.w1.npadedit.setText(self.savedparmsdict['npad'])
			
			self.w1.fourvarchkbx.setCheckState(self.savedparmsdict['fourvar'])
			self.w1.debugchkbx.setCheckState(self.savedparmsdict['debug'])
			self.w1.chunkedit.setText(self.savedparmsdict['chunk'])
				
	def setactiveheader(self):
		stack = self.stacknameedit.text()
		print "stack defined="+ stack
		# horatio active_refactoring Jy51i1EwmLD4tWZ9_00000_1
		# header(str(stack), "active", one=True)


	def setprojheader(self):
		#opens a file browser, showing files only in .hdf format
		ok=False
		zerostr='set xform.projection to zero'
		randstr='randomize xform.projection'
		importstr='import parameters from file'
		(item,stat)= QInputDialog.getItem(self,"xform.projection","choose option",[zerostr,randstr,importstr])
		#we convert this Qstring to a string and send it to line edit classed stackname edit of the Poptwodali window
		#self.stacknameedit.setText(str(a))
		choice= str(item)
		stack = self.stacknameedit.text()
		if stat:
			if choice == zerostr:
				header(str(stack),'xform.projection',zero=True)
			if choice == randstr:
				header(str(stack),'xform.projection',rand_alpha=True)
			if choice == importstr:
				file_name = QtGui.QFileDialog.getOpenFileName(self, "Open Data File", "", "(*)")
				a=str(QtCore.QString(file_name))
				if len(a)>0:
					header(str(stack),'xform.projection',fimport=a)
		
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

class Popupadvparams_localali3d_1(QWidget):
	def __init__(self,savedparms):
		QWidget.__init__(self)
		
		########################################################################################
		# layout parameters
		
		self.y1=10
		self.yspc = 4
		
		self.x1 = 10
		self.x2 = self.x1+220 #140
		self.x3 = self.x2+145 #285
		########################################################################################
		
		#Here we just set the window title
		#self.setWindowTitle('sxali3d advanced parameter selection')
		#Here we just set a label and its position in the window
		title1=QtGui.QLabel('<b>sxlocal_ali3d</b> - set advanced parameters related to CTF and search', self)
		title1.move(self.x1,self.y1)
		
		self.y1 += 30
		#Labels and Line Edits for User Input
		#Just a label
		#title2= QtGui.QLabel('<b>Advanced</b> parameters', self)
		#title2.move(self.x1,self.y1)
		
		#self.y1 += 30
		
		self.savedparmsdict=savedparms
		#Example for User input stack name
		#First create the label and define its position
		maskname= QtGui.QLabel('Mask', self)
		maskname.move(self.x1,self.y1)
		#Now add a line edit and define its position
		self.masknameedit=QtGui.QLineEdit(self)
		self.masknameedit.move(self.x2,self.y1)
		#Adds a default value for the line edit
		self.masknameedit.setText(self.savedparmsdict['maskname'])
		self.masknameedit.setToolTip("filename of the file containing 3D mask. If not provided, a 3D spherical mask will be created with radius equal to outer_radius. This mask is used to multiply the reference volume for calculation of reference projections.")
		
		self.mskfile_button = QtGui.QPushButton("Open Mask", self)
		self.mskfile_button.move(self.x3, self.y1-self.yspc)
		#Here we define, that when this button is clicked, it starts subfunction choose_file
		QtCore.QObject.connect(self.mskfile_button, QtCore.SIGNAL("clicked()"), self.choose_mskfile)
		
		self.y1 += 30
		
		center= QtGui.QLabel('Center type', self)
		center.move(self.x1,self.y1)
		self.centeredit=QtGui.QLineEdit(self)
		self.centeredit.move(self.x2,self.y1)
		self.centeredit.setText(self.savedparmsdict['center'])
		self.centeredit.setToolTip('-1: average shift method; 0: no centering; 1: center of gravity (default=-1)')
		
		self.y1 += 30
		
		fourvar= QtGui.QLabel('Compute fourier variance', self)
		fourvar.move(self.x1,self.y1)
		self.fourvarchkbx=QtGui.QCheckBox("",self)
		self.fourvarchkbx.move(self.x2,self.y1)
		self.fourvarchkbx.setCheckState(self.savedparmsdict['fourvar'])
		self.fourvarchkbx.setToolTip('compute Fourier variance')
		
		self.y1 += 30
		
		chunk= QtGui.QLabel('Chunk of data after \nwhich the 3D will be updated \n0<chunk<=1.0 (default 1.0)', self)
		chunk.move(self.x1,self.y1)
		self.chunkedit=QtGui.QLineEdit(self)
		self.chunkedit.move(self.x2,self.y1)
		self.chunkedit.setText(self.savedparmsdict['chunk'])
		self.chunkedit.setToolTip('Chunk of data after which the 3D will be updated 0<chunk<=1.0 (default 1.0)')		

		self.y1 += 70
		
		ctf= QtGui.QLabel('CTF', self)
		ctf.move(self.x1,self.y1)
		self.ctfchkbx = QtGui.QCheckBox("",self)
		self.ctfchkbx.move(self.x2, self.y1)
		self.ctfchkbx.setCheckState(self.savedparmsdict['ctf'])
		self.ctfchkbx.setToolTip('Consider CTF correction during the alignment ')
		
		self.y1 += 30
		
		snr= QtGui.QLabel('SNR', self)
		snr.move(self.x1,self.y1)
		self.snredit=QtGui.QLineEdit(self)
		self.snredit.move(self.x2,self.y1)
		self.snredit.setText(self.savedparmsdict['snr'])
		self.snredit.setToolTip('Signal-to-Noise Ratio of the data. Default is 1.0.')		
		
		self.y1 += 30
		
		sym= QtGui.QLabel('Symmetry', self)
		sym.move(self.x1,self.y1)
		self.symedit=QtGui.QLineEdit(self)
		self.symedit.move(self.x2,self.y1)
		self.symedit.setText(self.savedparmsdict['sym'])
		self.symedit.setToolTip('Symmetry group (default is c1).')
		
		self.y1 += 30
		
		pad= QtGui.QLabel('Padding size for 3D \nreconstruction', self)
		pad.move(self.x1,self.y1)
		self.npadedit=QtGui.QLineEdit(self)
		self.npadedit.move(self.x2,self.y1)
		self.npadedit.setText(self.savedparmsdict['npad'])
		self.npadedit.setToolTip('padding size for 3D \nreconstruction')
		
		self.y1 += 50
		
		usrfunc= QtGui.QLabel('User Function Name', self)
		usrfunc.move(self.x1,self.y1)
		self.usrfuncedit=QtGui.QLineEdit(self)
		self.usrfuncedit.move(self.x2,self.y1)
		self.usrfuncedit.setText(self.savedparmsdict['usrfunc'])
		self.usrfuncedit.setToolTip('name of the reference preparation function (ref_ali3d)')
		
		self.y1 += 30
				
		usrfuncfile= QtGui.QLabel('Enter name of external file containing user function:\n(Leave blank if file is not external to Sparx)', self)
		usrfuncfile.move(self.x1,self.y1)
		
		
		self.y1 += 40
		
		self.usrfuncfileedit=QtGui.QLineEdit(self)
		self.usrfuncfileedit.move(self.x2,self.y1)
		self.usrfuncfileedit.setText(self.savedparmsdict['usrfuncfile'])
		self.usrfuncfileedit.setToolTip('name of the external file containing user function')		
		#Function runsxali2d started when  the  RUN_button of the  Poptwodali window is clicked 
			
		#Here we define, that when this button is clicked, it starts subfunction choose_file
		self.usrfile_button = QtGui.QPushButton("Select File", self)
		self.usrfile_button.move(self.x3, self.y1-self.yspc)
		QtCore.QObject.connect(self.usrfile_button, QtCore.SIGNAL("clicked()"), self.choose_usrfile)
		
		self.y1 += 40
		
		dbg= QtGui.QLabel('Debug mode', self)
		dbg.move(self.x1,self.y1)
		self.debugchkbx=QtGui.QCheckBox("",self)
		self.debugchkbx.move(self.x2,self.y1)
		self.debugchkbx.setCheckState(self.savedparmsdict['debug'])
		self.debugchkbx.setToolTip('debug mode')
		
		self.y1 += 30
		
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

class Popupviper(QWidget):

	def __init__(self):
		QWidget.__init__(self)
		
		########################################################################################
		# class variables
		self.cmd = ""
		# populate with default values
		self.savedparmsdict = {'stackname':'','foldername':'','partradius':'-1','delta':'2.0','maxit':'100','moon_elimination':'','lf':'0.15','hf':'0.25','nproc':'1',"innerradius":"-1","debug":False,"given":False,"rand_seed":"-1","first_zero":False,"noweights":False,"pcross":"0.95","pmut":"0.05","maxgen":"10"}
		
		if options.demo == DEMO_mpibdbctf:
			self.savedparmsdict['stackname']='class_averages_generation_1.hdf'
			self.savedparmsdict['foldername']='structure'
			self.savedparmsdict['partradius']='30'
			self.savedparmsdict['delta']='6'
			self.savedparmsdict['lf']='0.05'
			self.savedparmsdict['hf']='0.25'
			self.savedparmsdict['maxit']='20'
			self.savedparmsdict['rand_seed']='10'
			self.savedparmsdict['nproc']='4'
		
		if options.demo == DEMO_mpibdb:
			self.savedparmsdict['stackname']='class_averages_generation_1.hdf'
			self.savedparmsdict['foldername']='structure'
			self.savedparmsdict['partradius']='30'
			self.savedparmsdict['delta']='6'
			self.savedparmsdict['lf']='0.05'
			self.savedparmsdict['hf']='0.3'
			self.savedparmsdict['maxit']='20'
			self.savedparmsdict['rand_seed']='200'
			self.savedparmsdict['nproc']='4'
		########################################################################################
		# layout parameters
		
		self.y1 = 10
		self.y2 = self.y1 + 98
		self.y4 = self.y2 + 450
		self.y5 = self.y4 + 95
		self.yspc = 4
		
		self.x1 = 10
		self.x2 = self.x1 + 200
		self.x3 = self.x2+145
		self.x4 = self.x3+100
		self.x5 = 230
		########################################################################################
		
		#Here we just set the window title
		self.setWindowTitle('sxviper')
		#Here we just set a label and its position in the window
		title1=QtGui.QLabel('<b>sxviper</b> -  use the "common line" algorithm to assign initial values of phi, theta, <br>psi to 2D average projections', self)
		title1.move(self.x1,self.y1)
		self.y1 += 50

		self.repopbtn = QPushButton("Retrieve saved parameters", self)
		self.repopbtn.move(self.x1-5,self.y1)
		#sets an infotip for this Pushbutton
		self.repopbtn.setToolTip('Retrieve saved parameters')
		#when this button is clicked, this action starts the subfunction twodali
		self.connect(self.repopbtn, SIGNAL("clicked()"), self.repoparms_viper)
		
		########################################################################################
		
		# stackname= QtGui.QLabel('Projection stack file', self)
		# stackname.move(self.x1,self.y2)
		# #Now add a line edit and define its position
		# self.stacknameedit=QtGui.QLineEdit(self)
		# self.stacknameedit.move(self.x2,self.y2)
		# #Adds a default value for the line edit
		# self.stacknameedit.setText(self.savedparmsdict['stackname'])
		# 
		# #Here we create a Button(file_button with title run open .hdf) and its position in the window
		# self.file_button = QtGui.QPushButton("Open .hdf", self)
		# self.file_button.move(self.x3, self.y2-self.yspc)
		# #Here we define, that when this button is clicked, it starts subfunction choose_file
		# QtCore.QObject.connect(self.file_button, QtCore.SIGNAL("clicked()"), self.choose_file)
		# #exactly the same as above, but for subfunction choose_file1
		# self.file_button1 = QtGui.QPushButton("Open .bdb", self)
		# self.file_button1.move(self.x4,self.y2-self.yspc)
		# QtCore.QObject.connect(self.file_button1, QtCore.SIGNAL("clicked()"), self.choose_file1)
		
		# self.y2 += 30
		
		# #The same as above, but many line edits include Infotips
		# foldername= QtGui.QLabel('Output folder', self)
		# foldername.move(self.x1,self.y2)
		# self.foldernameedit=QtGui.QLineEdit(self)
		# self.foldernameedit.move(self.x2,self.y2)
		# self.foldernameedit.setText(self.savedparmsdict['foldername'])		
		
		# self.outinfobtn = QPushButton("Output Info", self)
		# self.outinfobtn.move(self.x3,  self.y2-self.yspc)
		# #sets an infotip for this Pushbutton
		# self.outinfobtn.setToolTip('Output Info')
		# #when this button is clicked, this action starts the subfunction twodali
		# self.connect(self.outinfobtn, SIGNAL("clicked()"), self.outputinfo_viper)
		
		# self.y2 = self.y2+30
		

		args_list_and_parser =  self.get_args_list_and_parser()
		QLabels = [None]*len(args_list_and_parser)
		self.QLineEditsAndChecks = [None]*len(args_list_and_parser)
		# for option_iterator in range(2, len(parser.option_list)):
		for option_iterator in range(len(args_list_and_parser)):
			if args_list_and_parser[option_iterator].help == None or \
				args_list_and_parser[option_iterator].help == "" or \
				args_list_and_parser[option_iterator].help[0] != "<" or \
				args_list_and_parser[option_iterator].help[-10:] == "(advanced)":
				continue
			label = args_list_and_parser[option_iterator].help[1:].split(">")[0]
			# a = QtGui.QCheckBox("",self)
			QLabels[option_iterator] = QtGui.QLabel(label, self)
			QLabels[option_iterator].move(self.x1,self.y2)

			control_is_an_edit_box = args_list_and_parser[option_iterator].action != "store_true"
			control_is_for_a_parameter = args_list_and_parser[option_iterator].action == "LPlFHy3uNTjucRlk"
			if control_is_an_edit_box:
				self.QLineEditsAndChecks[option_iterator] = QtGui.QLineEdit(self)
				# modify_function = self.QLineEditsAndChecks[option_iterator].setText
				self.QLineEditsAndChecks[option_iterator].setText(str(args_list_and_parser[option_iterator].default))
			else:
				self.QLineEditsAndChecks[option_iterator] = QtGui.QCheckBox("",self)
				# modify_function = self.QLineEditsAndChecks[option_iterator].setCheckState
				modify_function = self.QLineEditsAndChecks[option_iterator].setCheckState(args_list_and_parser[option_iterator].default)

			self.QLineEditsAndChecks[option_iterator].move(self.x2,self.y2 - 7)
			# self.QLineEdits[option_iterator].setText(self.savedparmsdict[parser.option_list[option_iterator].dest])
			self.savedparmsdict[args_list_and_parser[option_iterator].dest] = [option_iterator, str(args_list_and_parser[option_iterator].default), control_is_an_edit_box, control_is_for_a_parameter]
			# modify_function(str(parser.option_list[option_iterator].default))
			self.QLineEditsAndChecks[option_iterator].setToolTip(args_list_and_parser[option_iterator].help)		
	
			self.y2 = self.y2+25


		nproc= QtGui.QLabel('MPI processors', self)
		nproc.move(self.x1,self.y2)
		self.nprocedit=QtGui.QLineEdit(self)
		self.nprocedit.move(self.x2,self.y2)
		self.nprocedit.setText(self.savedparmsdict['nproc'])
		self.nprocedit.setToolTip('The number of processors to use. Default is single processor mode')

		self.y2 =self.y2+25
		

		
		########################################################################################
		
		self.savepbtn = QPushButton("Save Input Parameters", self)
		self.savepbtn.move(self.x1-5,  self.y4)
		#sets an infotip for this Pushbutton
		self.savepbtn.setToolTip('Save Input Parameters')
		#when this button is clicked, this action starts the subfunction twodali
		self.connect(self.savepbtn, SIGNAL("clicked()"), self.saveparms)
		
		self.y4 = self.y4+30
		
		self.cmdlinebtn = QPushButton("Generate command line from input parameters", self)
		self.cmdlinebtn.move(self.x1-5,  self.y4)
		#sets an infotip for this Pushbutton
		self.cmdlinebtn.setToolTip('Generate command line using input parameters')
		#when this button is clicked, this action starts the subfunction twodali
		self.connect(self.cmdlinebtn, SIGNAL("clicked()"), self.gencmdline_viper)
		
		########################################################################################
		
		 #Here we create a Button(Run_button with title run sxali2d) and its position in the window
		self.RUN_button = QtGui.QPushButton('Run sxviper', self)
		# make 3D textured push button look
		s = "QPushButton {font: bold; color: #000;border: 1px solid #333;border-radius: 11px;padding: 2px;background: qradialgradient(cx: 0, cy: 0,fx: 0.5, fy:0.5,radius: 1, stop: 0 #fff, stop: 1 #8D0);min-width:90px;margin:5px} QPushButton:pressed {font: bold; color: #000;border: 1px solid #333;border-radius: 11px;padding: 2px;background: qradialgradient(cx: 0, cy: 0,fx: 0.5, fy:0.5,radius: 1, stop: 0 #fff, stop: 1 #084);min-width:90px;margin:5px}"
		
		self.RUN_button.setStyleSheet(s)
		
		self.RUN_button.move(self.x5,  self.y5)
		#Here we define, that when this button is clicked, it starts subfunction runsxali2d
		self.connect(self.RUN_button, SIGNAL("clicked()"), self.runsxviper)



	# def __init__(self):
	# 	QWidget.__init__(self)
	# 	
	# 	########################################################################################
	# 	# class variables
	# 	self.cmd = ""
	# 	# populate with default values
	# 	self.savedparmsdict = {'stackname':'','foldername':'','partradius':'-1','delta':'2.0','maxit':'100','moon_elimination':'','lf':'0.15','hf':'0.25','nproc':'1',"innerradius":"-1","debug":False,"given":False,"rand_seed":"-1","first_zero":False,"noweights":False,"pcross":"0.95","pmut":"0.05","maxgen":"10"}
	# 	
	# 	if options.demo == DEMO_mpibdbctf:
	# 		self.savedparmsdict['stackname']='class_averages_generation_1.hdf'
	# 		self.savedparmsdict['foldername']='structure'
	# 		self.savedparmsdict['partradius']='30'
	# 		self.savedparmsdict['delta']='6'
	# 		self.savedparmsdict['lf']='0.05'
	# 		self.savedparmsdict['hf']='0.25'
	# 		self.savedparmsdict['maxit']='20'
	# 		self.savedparmsdict['rand_seed']='10'
	# 		self.savedparmsdict['nproc']='4'
	# 	
	# 	if options.demo == DEMO_mpibdb:
	# 		self.savedparmsdict['stackname']='class_averages_generation_1.hdf'
	# 		self.savedparmsdict['foldername']='structure'
	# 		self.savedparmsdict['partradius']='30'
	# 		self.savedparmsdict['delta']='6'
	# 		self.savedparmsdict['lf']='0.05'
	# 		self.savedparmsdict['hf']='0.3'
	# 		self.savedparmsdict['maxit']='20'
	# 		self.savedparmsdict['rand_seed']='200'
	# 		self.savedparmsdict['nproc']='4'
	# 	########################################################################################
	# 	# layout parameters
	# 	
	# 	self.y1 = 10
	# 	self.y2 = self.y1 + 98
	# 	self.y4 = self.y2 + 450
	# 	self.y5 = self.y4 + 95
	# 	self.yspc = 4
	# 	
	# 	self.x1 = 10
	# 	self.x2 = self.x1 + 200
	# 	self.x3 = self.x2+145
	# 	self.x4 = self.x3+100
	# 	self.x5 = 230
	# 	########################################################################################
	# 	
	# 	#Here we just set the window title
	# 	self.setWindowTitle('sxviper')
	# 	#Here we just set a label and its position in the window
	# 	title1=QtGui.QLabel('<b>sxviper</b> -  use the "common line" algorithm to assign initial values of phi, theta, <br>psi to 2D average projections', self)
	# 	title1.move(self.x1,self.y1)
	# 	self.y1 += 50
	# 
	# 	self.repopbtn = QPushButton("Retrieve saved parameters", self)
	# 	self.repopbtn.move(self.x1-5,self.y1)
	# 	#sets an infotip for this Pushbutton
	# 	self.repopbtn.setToolTip('Retrieve saved parameters')
	# 	#when this button is clicked, this action starts the subfunction twodali
	# 	self.connect(self.repopbtn, SIGNAL("clicked()"), self.repoparms_viper)
	# 	
	# 	########################################################################################
	# 	
	# 	stackname= QtGui.QLabel('Projection stack file', self)
	# 	stackname.move(self.x1,self.y2)
	# 	#Now add a line edit and define its position
	# 	self.stacknameedit=QtGui.QLineEdit(self)
	# 	self.stacknameedit.move(self.x2,self.y2)
	# 	#Adds a default value for the line edit
	# 	self.stacknameedit.setText(self.savedparmsdict['stackname'])
	# 	
	# 	#Here we create a Button(file_button with title run open .hdf) and its position in the window
	# 	self.file_button = QtGui.QPushButton("Open .hdf", self)
	# 	self.file_button.move(self.x3, self.y2-self.yspc)
	# 	#Here we define, that when this button is clicked, it starts subfunction choose_file
	# 	QtCore.QObject.connect(self.file_button, QtCore.SIGNAL("clicked()"), self.choose_file)
	# 	#exactly the same as above, but for subfunction choose_file1
	# 	self.file_button1 = QtGui.QPushButton("Open .bdb", self)
	# 	self.file_button1.move(self.x4,self.y2-self.yspc)
	# 	QtCore.QObject.connect(self.file_button1, QtCore.SIGNAL("clicked()"), self.choose_file1)
	# 	
	# 	self.y2 += 30
	# 	
	# 	#The same as above, but many line edits include Infotips
	# 	foldername= QtGui.QLabel('Output folder', self)
	# 	foldername.move(self.x1,self.y2)
	# 	self.foldernameedit=QtGui.QLineEdit(self)
	# 	self.foldernameedit.move(self.x2,self.y2)
	# 	self.foldernameedit.setText(self.savedparmsdict['foldername'])		
	# 	
	# 	self.outinfobtn = QPushButton("Output Info", self)
	# 	self.outinfobtn.move(self.x3,  self.y2-self.yspc)
	# 	#sets an infotip for this Pushbutton
	# 	self.outinfobtn.setToolTip('Output Info')
	# 	#when this button is clicked, this action starts the subfunction twodali
	# 	self.connect(self.outinfobtn, SIGNAL("clicked()"), self.outputinfo_viper)
	# 	
	# 	self.y2 = self.y2+30
	# 	
	# 	partradius= QtGui.QLabel('Outer radius of particle \n< int(nx/2)-1', self)
	# 	partradius.move(self.x1,self.y2)
	# 	self.partradiusedit=QtGui.QLineEdit(self)
	# 	self.partradiusedit.move(self.x2,self.y2)
	# 	self.partradiusedit.setText(self.savedparmsdict['partradius'])
	# 	self.partradiusedit.setToolTip('Parameter ou: Outer radius for rotational correlation \nshould be set to particle radius\nif not sure, set to boxsize/2-2 ')		
	# 
	# 	self.y2 = self.y2+50
	# 	
	# 	innerradius= QtGui.QLabel('Inner radius of particle \n(set to -1)', self)
	# 	innerradius.move(self.x1,self.y2)
	# 	self.innerradiusedit=QtGui.QLineEdit(self)
	# 	self.innerradiusedit.move(self.x2,self.y2)
	# 	self.innerradiusedit.setText(self.savedparmsdict['innerradius'])
	# 	self.innerradiusedit.setToolTip('Inner radius of particle (set to -1)')		
	# 	
	# 	self.y2 += 50
	# 	
	# 	delta= QtGui.QLabel('Angular step', self)
	# 	delta.move(self.x1,self.y2)
	# 	self.deltaedit=QtGui.QLineEdit(self)
	# 	self.deltaedit.move(self.x2,self.y2)
	# 	self.deltaedit.setText(self.savedparmsdict['delta'])
	# 	self.deltaedit.setToolTip('angular step for the reference projections in respective iterations')				
	# 	
	# 	self.y2 =self.y2+30
	# 
	# 	moon_eliminationlab= QtGui.QLabel('Moon elimination', self)
	# 	moon_eliminationlab.move(self.x1,self.y2)
	# 	self.moon_elimination_edit=QtGui.QLineEdit(self)
	# 	self.moon_elimination_edit.move(self.x2,self.y2)
	# 	self.moon_elimination_edit.setText(self.savedparmsdict['moon_elimination'])
	# 	self.moon_elimination_edit.setToolTip('mass in KDa and resolution in px/A separated by comma, no space, leave empty to disable')				
	# 	
	# 	self.y2 =self.y2+50
	# 	
	# 	lflab= QtGui.QLabel('Filter, minimum frequency \n(set to 0.0)', self)
	# 	lflab.move(self.x1,self.y2)
	# 	self.lfedit=QtGui.QLineEdit(self)
	# 	self.lfedit.move(self.x2,self.y2)
	# 	self.lfedit.setText(self.savedparmsdict['lf'])
	# 	self.lfedit.setToolTip('Filter, minimum frequency (set to 0.0)')				
	# 	
	# 	self.y2 =self.y2+50
	# 	
	# 	hflab= QtGui.QLabel('Filter, maximum frequency \n(set to 0.5)', self)
	# 	hflab.move(self.x1,self.y2)
	# 	self.hfedit=QtGui.QLineEdit(self)
	# 	self.hfedit.move(self.x2,self.y2)
	# 	self.hfedit.setText(self.savedparmsdict['hf'])
	# 	self.hfedit.setToolTip('Filter, maximum frequency (set to 0.5)')				
	# 	
	# 	self.y2 =self.y2+50
	# 	
	# 	maxit= QtGui.QLabel('Maximum number of iterations', self)
	# 	maxit.move(self.x1,self.y2)
	# 	self.maxitedit=QtGui.QLineEdit(self)
	# 	self.maxitedit.move(self.x2,self.y2)
	# 	self.maxitedit.setText(self.savedparmsdict['maxit'])
	# 	self.maxitedit.setToolTip('Maximum number of iterations')
	# 	
	# 	self.y2 =self.y2+50
	# 	
	# 	nproc= QtGui.QLabel('MPI processors', self)
	# 	nproc.move(self.x1,self.y2)
	# 	self.nprocedit=QtGui.QLineEdit(self)
	# 	self.nprocedit.move(self.x2,self.y2)
	# 	self.nprocedit.setText(self.savedparmsdict['nproc'])
	# 	self.nprocedit.setToolTip('The number of processors to use. Default is single processor mode')
	# 	
	# 	########################################################################################
	# 	
	# 	self.savepbtn = QPushButton("Save Input Parameters", self)
	# 	self.savepbtn.move(self.x1-5,  self.y4)
	# 	#sets an infotip for this Pushbutton
	# 	self.savepbtn.setToolTip('Save Input Parameters')
	# 	#when this button is clicked, this action starts the subfunction twodali
	# 	self.connect(self.savepbtn, SIGNAL("clicked()"), self.saveparms)
	# 	
	# 	self.y4 = self.y4+30
	# 	
	# 	self.cmdlinebtn = QPushButton("Generate command line from input parameters", self)
	# 	self.cmdlinebtn.move(self.x1-5,  self.y4)
	# 	#sets an infotip for this Pushbutton
	# 	self.cmdlinebtn.setToolTip('Generate command line using input parameters')
	# 	#when this button is clicked, this action starts the subfunction twodali
	# 	self.connect(self.cmdlinebtn, SIGNAL("clicked()"), self.gencmdline_viper)
	# 	
	# 	########################################################################################
	# 	
	# 	 #Here we create a Button(Run_button with title run sxali2d) and its position in the window
	# 	self.RUN_button = QtGui.QPushButton('Run sxviper', self)
	# 	# make 3D textured push button look
	# 	s = "QPushButton {font: bold; color: #000;border: 1px solid #333;border-radius: 11px;padding: 2px;background: qradialgradient(cx: 0, cy: 0,fx: 0.5, fy:0.5,radius: 1, stop: 0 #fff, stop: 1 #8D0);min-width:90px;margin:5px} QPushButton:pressed {font: bold; color: #000;border: 1px solid #333;border-radius: 11px;padding: 2px;background: qradialgradient(cx: 0, cy: 0,fx: 0.5, fy:0.5,radius: 1, stop: 0 #fff, stop: 1 #084);min-width:90px;margin:5px}"
	# 	
	# 	self.RUN_button.setStyleSheet(s)
	# 	
	# 	self.RUN_button.move(self.x5,  self.y5)
	# 	#Here we define, that when this button is clicked, it starts subfunction runsxali2d
	# 	self.connect(self.RUN_button, SIGNAL("clicked()"), self.runsxviper)
		#Labels and Line Edits for User Input

		
	#Function runsxali2d started when  the  RUN_button of the  Poptwodali window is clicked 

	def get_args_list_and_parser(self):

		from sxviper import main as sxviper_main
		parser = sxviper_main(["aa", "--return_options"])

		# ['__doc__', '__init__', '__module__', '_add_help_option', '_add_version_option', '_check_conflict', '_create_option_list', '_create_option_mappings', '_get_all_options', '_get_args', '_get_encoding', '_init_parsing_state', '_long_opt', '_match_long_opt', '_populate_option_list', '_process_args', '_process_long_opt', '_process_short_opts', '_share_option_mappings', '_short_opt', 'add_option', 'add_option_group', 'add_options', 'allow_interspersed_args', 'check_values', 'conflict_handler', 'defaults', 'description', 'destroy', 'disable_interspersed_args', 'enable_interspersed_args', 'epilog', 'error', 'exit', 'expand_prog_name', 'format_description', 'format_epilog', 'format_help', 'format_option_help', 'formatter', 'get_default_values', 'get_description', 'get_option', 'get_option_group', 'get_prog_name', 'get_usage', 'get_version', 'has_option', 'largs', 'option_class', 'option_groups', 'option_list', 'parse_args', 'print_help', 'print_usage', 'print_version', 'process_default_values', 'prog', 'rargs', 'remove_option', 'set_conflict_handler', 'set_default', 'set_defaults', 'set_description', 'set_process_default_values', 'set_usage', 'standard_option_list', 'usage', 'values', 'version']
		# >>> a.optionlist
		# Traceback (most recent call last):
		#   File "<stdin>", line 1, in <module>
		# AttributeError: OptionParser instance has no attribute 'optionlist'
		# >>> len(a.option_list)
		# 25
		# >>> a.option_list[1].help
		# 'show this help message and exit'
		# >>> a.option_list[11].help
		# 'maximum number of iterations performed for the finishing up part (set to 50) '
		# >>>
		# >>> dir(a.option_list[11])
		# ['ACTIONS', 'ALWAYS_TYPED_ACTIONS', 'ATTRS', 'CHECK_METHODS', 'CONST_ACTIONS', 'STORE_ACTIONS', 'TYPED_ACTIONS', 'TYPES', 'TYPE_CHECKER', '__doc__', '__init__', '__module__', '__repr__', '__str__', '_check_action', '_check_callback', '_check_choice', '_check_const', '_check_dest', '_check_nargs', '_check_opt_strings', '_check_type', '_long_opts', '_set_attrs', '_set_opt_strings', '_short_opts', 'action', 'callback', 'callback_args', 'callback_kwargs', 'check_value', 'choices', 'const', 'container', 'convert_value', 'default', 'dest', 'get_opt_string', 'help', 'metavar', 'nargs', 'process', 'take_action', 'takes_value', 'type']
		# >>> a.option_list[11].dest
		# 'maxit2'
		
		class ImmitateOptionList:
			default=""
			help=""
			dest = ""
			action = "LPlFHy3uNTjucRlk"

		prog_args = parser.usage.split("--")[0].split(" ")
		# print prog_args
		args_list = []
		for a in prog_args[1:]:
			if a == "": continue
			aa = a.strip()
			b = ImmitateOptionList()
			b.help = "<" + aa + ">"
			b.dest = aa[:]
			if a == "stack":
				b.help = "<Projection stack file>"
			args_list.append(b)
		
		return args_list + parser.option_list[2:]


	def outputinfo_viper(self):
		QMessageBox.information(self, "sxviper output",'outdir is the name of the output folder specified by the user. If it does not exist, the directory will be created. If it does exist, the program will crash and an error message will come up. Please change the name of directory and restart the program. \n\noutdir/angles_000: \n\nThis file contains Eulerian angles fount in trial #000 \n\noutdir/plot_agls_000.hdf: \n\nThis image in the hdf format contains visualization of the distribution of projections found during trial #000 (see also sxplot_projs_distrib) \n\noutdir/structure_000.hdf: \n\nCopy of the stack of 2D projections with Eulerian angles found at trial #000 set in the header. In order to examine the structure, one has to do the 3D reconstructions sxrecons3d_n.py outdir/structure_000.hdf myvol.hdf \n\noutdir/structure.hdf: \n\nFor multiple trials, it is a copy of the stack of 2D projections with Eulerian angles found at the best trial set in the header (this feature is no longer supported).')
		

	def gencmdline_viper(self,writefile=True):
		#Here we just read in all user inputs in the line edits of the Poptwodali window
		# stack = self.stacknameedit.text()
		# output=self.foldernameedit.text()
		# cmd1 = "sxviper.py "+str(stack) +" "+str(output)
		cmd1 = "sxviper.py "

		args = " "
		for param in self.get_args_list_and_parser():
			if param.action != "LPlFHy3uNTjucRlk": break
			args += "%s "%self.QLineEditsAndChecks[self.savedparmsdict[param.dest][0]].text()
		
		for key in self.savedparmsdict:
			if type(self.savedparmsdict[key]) != list:
				continue
			if self.savedparmsdict[key][2]: ## check to see if this is not a boolean option
				if not self.savedparmsdict[key][3]: ## check to see if this is not a parameter
					args += "--%s=%s "%(key,self.QLineEditsAndChecks[self.savedparmsdict[key][0]].text())
				else:
					args += "%s "%self.QLineEditsAndChecks[self.savedparmsdict[key][0]].text()
				self.savedparmsdict[key][1] = self.QLineEditsAndChecks[self.savedparmsdict[key][0]].text()
			else:
				if self.QLineEditsAndChecks[self.savedparmsdict[key][0]].checkState() == Qt.Checked:
					args += "--%s "%key
					self.savedparmsdict[key][1] = self.QLineEditsAndChecks[self.savedparmsdict[key][0]].checkState()
		cmd1 = cmd1 + args

		args = " "
		for key in self.w1.savedparmsdict:
			if type(self.w1.savedparmsdict[key]) != list:
				continue
			# print self.savedparmsdict
			if self.w1.savedparmsdict[key][2]:
				val_str = self.w1.QLineEditsAndChecks[self.w1.savedparmsdict[key][0]].text()
				self.w1.savedparmsdict[key][1] = val_str
				if val_str != "":
					args += "--%s=%s "%(key,val_str)
			else:
				if self.w1.QLineEditsAndChecks[self.w1.savedparmsdict[key][0]].checkState() == Qt.Checked:
					args += "--%s "%key
					self.w1.savedparmsdict[key][1] = self.w1.QLineEditsAndChecks[self.w1.savedparmsdict[key][0]].checkState()
		cmd1 = cmd1 + args


		self.savedparmsdict['nproc'] = self.nprocedit.text()
		
		# np = self.QLineEditsAndChecks[int(self.savedparmsdict["nproc"])].text()
		
		np = self.nprocedit.text()
		
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



	# def gencmdline_viper(self,writefile=True):
	# 	#Here we just read in all user inputs in the line edits of the Poptwodali window
	# 	stack = self.stacknameedit.text()
	# 	output=self.foldernameedit.text()
	# 	ou=self.partradiusedit.text()
	# 	delta=self.deltaedit.text()
	# 	maxit=self.maxitedit.text()
	# 	moon_elimination=self.moon_elimination_edit.text()
	# 	lf=self.lfedit.text()
	# 	hf=self.hfedit.text()
	# 	inrad = self.innerradiusedit.text()
	# 	
	# 	cmd1 = "sxviper.py "+str(stack) +" "+str(output)
	# 	
	# 	args = " --ou="+ str(ou)+ " --ir=" + str(inrad)+" --delta='"+str(delta)+"'"+" --moon_elimination='"+str(moon_elimination)+"'"+" --maxit="+ str(maxit) +" --lf="+ str(lf)+" --hf="+ str(hf)
	# 	
	# 	cmd1 = cmd1 + args
	# 	
	# 	rand_seed=self.w1.rand_seededit.text()
	# 	given = self.w1.givenchkbx.checkState()		
	# 	debug = self.w1.debugchkbx.checkState()
	# 	first_zero = self.w1.first_zerochkbx.checkState()
	# 	noweights = self.w1.noweightschkbx.checkState()
	# 	pcross = self.w1.pcrossedit.text()
	# 	pmut = self.w1.pmutedit.text()
	# 	maxgen = self.w1.maxgenedit.text()
	# 	
	# 	cmd1 = cmd1+" --rand_seed="+str(rand_seed)+" --pcross="+str(pcross)+" --pmut="+str(pmut)+" --maxgen="+str(maxgen)
	# 	
	# 	if debug == Qt.Checked:
	# 		cmd1 = cmd1 + " --debug"
	# 	if given == Qt.Checked:
	# 		cmd1 = cmd1 + " --given"
	# 	if first_zero == Qt.Checked:
	# 		cmd1 = cmd1 + " --first_zero"
	# 	if noweights == Qt.Checked:
	# 		cmd1 = cmd1 + " --noweights"
	# 							
	# 	np = self.nprocedit.text()
	# 	
	# 	self.savedparmsdict = {'stackname':str(stack),'foldername':str(output),'partradius':str(ou),'delta':str(delta),'moon_elimination':str(moon_elimination),'lf':str(lf),'hf':str(hf),'maxit':str(maxit),'nproc':str(np),"innerradius":str(inrad),"debug":debug,"given":given,"rand_seed":str(rand_seed),"first_zero":first_zero,"noweights":noweights,"pcross":str(pcross),"pmut":str(pmut),"maxgen":str(maxgen)}
	# 
	# 	self.w1.savedparmsdict=self.savedparmsdict
	# 	
	# 	if int(str(np)) > 1: 
	# 		cmd1="mpirun -np "+ str(np) + " "+ cmd1+" --MPI" 
	# 	
	# 	if writefile:		
	# 		(fname,stat)= QInputDialog.getText(self,"Generate Command Line","Enter name of file to save command line in",QLineEdit.Normal,"")
	# 		if stat:
	# 			f = open(fname,'a')
	# 			f.write(cmd1)
	# 			f.write('\n')
	# 			f.close()
	# 	
	# 	print cmd1
	# 	self.cmd = cmd1
		
	def runsxviper(self):
		self.gencmdline_viper(writefile=False)
		outfolder=self.savedparmsdict['foldername']
		if os.path.exists(outfolder):
			print "output folder "+outfolder+" already exists!"
			return
		process = subprocess.Popen(self.cmd,shell=True)
		self.emit(QtCore.SIGNAL("process_started"),process.pid)
		
	def saveparms(self):		
		(fname,stat)= QInputDialog.getText(self,"Save Input Parameters","Enter name of file to save parameters in",QLineEdit.Normal,"")
		if stat:
			import pickle
			output=open(fname,'wb')
			self.gencmdline_viper(writefile=False)
			pickle.dump((self.savedparmsdict, self.w1.savedparmsdict),output)
			output.close()
		

	def repoparms_viper(self):		
		# repopulate with saved parms
		(fname,stat)= QInputDialog.getText(self,"Retrieve saved parameters","Enter name of file parameters were saved in",QLineEdit.Normal,"")
		print (fname,stat)
		
		if stat:
			import pickle
			pkl = open(fname,'rb')
			# self.savedparmsdict = pickle.load(pkl)
			(self.savedparmsdict, self.w1.savedparmsdict) = pickle.load(pkl)
			# print self.savedparmsdict
			# self.stacknameedit.setText(self.savedparmsdict['stackname'])
			# self.foldernameedit.setText(self.savedparmsdict['foldername'])
			for key in self.savedparmsdict:
				if type(self.savedparmsdict[key]) != list:
					continue
				if self.savedparmsdict[key][2]:
					self.QLineEditsAndChecks[self.savedparmsdict[key][0]].setText(self.savedparmsdict[key][1])
				else:
					# print self.savedparmsdict[key]
					self.QLineEditsAndChecks[self.savedparmsdict[key][0]].setChecked(self.savedparmsdict[key][1] == Qt.Checked)
			for key in self.w1.savedparmsdict:
				if type(self.w1.savedparmsdict[key]) != list:
					continue
				if self.w1.savedparmsdict[key][2]:
					self.w1.QLineEditsAndChecks[self.w1.savedparmsdict[key][0]].setText(self.w1.savedparmsdict[key][1])
				else:
					self.w1.QLineEditsAndChecks[self.w1.savedparmsdict[key][0]].setChecked(self.w1.savedparmsdict[key][1] == Qt.Checked)
		pass	
		self.nprocedit.setText(self.savedparmsdict['nproc'])





	# def repoparms_viper(self):		
	# 	# repopulate with saved parms
	# 	(fname,stat)= QInputDialog.getText(self,"Retrieve saved parameters","Enter name of file parameters were saved in",QLineEdit.Normal,"")
	# 	if stat:
	# 		import pickle
	# 		pkl = open(fname,'rb')
	# 		self.savedparmsdict = pickle.load(pkl)
	# 		print self.savedparmsdict
	# 		self.partradiusedit.setText(self.savedparmsdict['partradius'])
	# 		self.stacknameedit.setText(self.savedparmsdict['stackname'])
	# 		self.foldernameedit.setText(self.savedparmsdict['foldername'])		
	# 		self.deltaedit.setText(self.savedparmsdict['delta'])
	# 		self.moon_elimination_edit.setText(self.savedparmsdict['moon_elimination'])
	# 		self.lfedit.setText(self.savedparmsdict['lf'])
	# 		self.hfedit.setText(self.savedparmsdict['hf'])
	# 		self.maxitedit.setText(self.savedparmsdict['maxit'])
	# 		self.nprocedit.setText(self.savedparmsdict['nproc'])
	# 		self.innerradiusedit.setText(self.savedparmsdict['innerradius'])
	# 		
	# 		self.w1.rand_seededit.setText(self.savedparmsdict['rand_seed'])
	# 		self.w1.givenchkbx.setCheckState(self.savedparmsdict['given'])
	# 		self.w1.debugchkbx.setCheckState(self.savedparmsdict['debug'])
	# 		self.w1.first_zerochkbx.setCheckState(self.savedparmsdict['first_zero'])
	# 		self.w1.noweightschkbx.setCheckState(self.savedparmsdict['noweights'])
	# 		self.w1.pcrossedit.setText(self.savedparmsdict['pcross'])
	# 		self.w1.pmutedit.setText(self.savedparmsdict['pmut'])
	# 		self.w1.maxgenedit.setText(self.savedparmsdict['maxgen'])
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


class Popupadvparams_viper(QWidget):
	def __init__(self,savedparms):
		QWidget.__init__(self)
		
		# self.savedparmsdict=savedparms
		self.savedparmsdict=dict()
		
		########################################################################################
		# layout parameters
		
		self.y1=10
		self.yspc = 4
		
		self.x1 = 20
		self.x2 = self.x1+280
		self.x3 = self.x2+145
		########################################################################################
		
		#Here we just set the window title
		#self.setWindowTitle('sxali3d advanced parameter selection')
		#Here we just set a label and its position in the window
		title1=QtGui.QLabel('<b>sxviper</b> - set advanced parameters', self)
		title1.move(self.x1,self.y1)
		self.y1 = self.y1+25
		

		from sxviper import main as sxviper_main
		parser = sxviper_main(["aa", "--return_options"])
		
		
		QLabels = [None]*len(parser.option_list)
		self.QLineEditsAndChecks = [None]*len(parser.option_list)
		for option_iterator in range(2, len(parser.option_list)):
		# for option_iterator in range(2, 6):
			if parser.option_list[option_iterator].help == None or \
				parser.option_list[option_iterator].help == "" or \
				parser.option_list[option_iterator].help[0] != "<" or \
				parser.option_list[option_iterator].help[-10:] != "(advanced)":
				continue
			label = parser.option_list[option_iterator].help[1:].split(">")[0]
		
			QLabels[option_iterator] = QtGui.QLabel(label, self)
			QLabels[option_iterator].move(self.x1,self.y1)
		
			control_is_an_edit_box = parser.option_list[option_iterator].action != "store_true"
			if control_is_an_edit_box:
				self.QLineEditsAndChecks[option_iterator] = QtGui.QLineEdit(self)
				modify_function = self.QLineEditsAndChecks[option_iterator].setText
			else:
				self.QLineEditsAndChecks[option_iterator] = QtGui.QCheckBox("",self)
				modify_function = self.QLineEditsAndChecks[option_iterator].setCheckState
		
			self.QLineEditsAndChecks[option_iterator].move(self.x2,self.y1)
			# self.QLineEdits[option_iterator].setText(self.savedparmsdict[parser.option_list[option_iterator].dest])
			self.savedparmsdict[parser.option_list[option_iterator].dest] = [option_iterator, str(parser.option_list[option_iterator].default), control_is_an_edit_box]
			modify_function(str(parser.option_list[option_iterator].default))
			self.QLineEditsAndChecks[option_iterator].setToolTip(parser.option_list[option_iterator].help)		
		
			self.y1 = self.y1+25
		
		
		# self.y1 += 30
		# #Labels and Line Edits for User Input
		# #Just a label
		# #title2= QtGui.QLabel('<b>Advanced</b> parameters', self)
		# #title2.move(self.x1,self.y1)
		# 
		# #self.y1 += 30
		# 
		# self.savedparmsdict=savedparms
		# #Example for User input stack name
		# #First create the label and define its position
		# 
		# 
		# gven= QtGui.QLabel('Start from given projections orientation \n(set to False, means start with randomize \norientations)', self)
		# gven.move(self.x1,self.y1)
		# self.givenchkbx=QtGui.QCheckBox("",self)
		# self.givenchkbx.move(self.x2,self.y1)
		# self.givenchkbx.setCheckState(self.savedparmsdict['given'])
		# self.givenchkbx.setToolTip('Start from given projections orientation (set to False, means start with randomize orientations)')
		# 
		# self.y1 += 70
		# 
		# rs= QtGui.QLabel('Random seed of initial orientations \n(if set to randomly)', self)
		# rs.move(self.x1,self.y1)
		# self.rand_seededit=QtGui.QLineEdit("",self)
		# self.rand_seededit.move(self.x2,self.y1)
		# self.rand_seededit.setText(self.savedparmsdict['rand_seed'])
		# self.rand_seededit.setToolTip('Random seed of initial orientations (if set to randomly)')
		# 
		# self.y1 += 50
		# 
		# fz= QtGui.QLabel('Assign the first projection orientation to 0', self)
		# fz.move(self.x1,self.y1)
		# self.first_zerochkbx=QtGui.QCheckBox("",self)
		# self.first_zerochkbx.move(self.x2,self.y1)
		# self.first_zerochkbx.setCheckState(self.savedparmsdict['first_zero'])
		# self.first_zerochkbx.setToolTip('Assign the first projection orientation to 0')
		# 
		# self.y1 += 30
		# 
		# nw= QtGui.QLabel('Use Voronoi weighting (by default use \nweights)', self)
		# nw.move(self.x1,self.y1)
		# self.noweightschkbx=QtGui.QCheckBox("",self)
		# self.noweightschkbx.move(self.x2,self.y1)
		# self.noweightschkbx.setCheckState(self.savedparmsdict['noweights'])
		# self.noweightschkbx.setToolTip('Use Voronoi weighting (by default use weights)')
		# 
		# self.y1 += 50
		# 
		# 
		# pc= QtGui.QLabel('Cross-over probability (set to 0.95)', self)
		# pc.move(self.x1,self.y1)
		# self.pcrossedit=QtGui.QLineEdit("",self)
		# self.pcrossedit.move(self.x2,self.y1)
		# self.pcrossedit.setText(self.savedparmsdict['pcross'])
		# self.pcrossedit.setToolTip('Cross-over probability (set to 0.95)')
		# 
		# self.y1 += 30
		# 
		# pm= QtGui.QLabel('Mutation probability (set to 0.05)', self)
		# pm.move(self.x1,self.y1)
		# self.pmutedit=QtGui.QLineEdit("",self)
		# self.pmutedit.move(self.x2,self.y1)
		# self.pmutedit.setText(self.savedparmsdict['pmut'])
		# self.pmutedit.setToolTip('Mutation probability (set to 0.05)')
		# 
		# self.y1 += 30
		# 
		# mg= QtGui.QLabel('Maximum number of generations \n(set to 10)', self)
		# mg.move(self.x1,self.y1)
		# self.maxgenedit=QtGui.QLineEdit(self)
		# self.maxgenedit.move(self.x2,self.y1)
		# self.maxgenedit.setText(self.savedparmsdict['maxgen'])
		# self.maxgenedit.setToolTip('Maximum number of generations (set to 10)')
		# 
		# self.y1 += 50
		# 
		# dbg= QtGui.QLabel('Debug mode', self)
		# dbg.move(self.x1,self.y1)
		# self.debugchkbx=QtGui.QCheckBox("",self)
		# self.debugchkbx.move(self.x2,self.y1)
		# self.debugchkbx.setCheckState(self.savedparmsdict['debug'])
		# self.debugchkbx.setToolTip('Debug')
		# 
		# self.y1 += 30



class Popupvariability3d(QWidget):
	def __init__(self):
		QWidget.__init__(self)
		
		#######################################################################################
		# class variables
		self.cmd = ""
		# populate with default values
		self.savedparmsdict = {'pdbfile':'','output':'','apix':'1.0','box':'150','het':Qt.Unchecked,'center':'n','Och':Qt.Unchecked,'quiet':Qt.Unchecked,'tr0':''}
		#######################################################################################
		# Layout parameters
		
		self.y1 = 10 # title and Repopulate button
		self.y2 = self.y1 + 85 # Text boxes for inputting parameters
		#self.y3 = self.y2 + 222 # activate images button and set xform.align2d button
		self.y4 = self.y2 + 290 # Advanced Parameters, Save Input and Generate command line buttons
		self.y5 = self.y4 + 95 # run button 
		self.yspc = 4
		
		self.x1 = 10 # first column (text box labels)
		self.x2 = self.x1 + 260 # second column (text boxes)
		self.x3 = self.x2+145 # third column (Open .hdf button)
		self.x4 = self.x3+100 # fourth column (Open .bdb button)
		self.x5 = 230 # run button
		#######################################################################################
		
		#Here we just set the window title
		self.setWindowTitle('sxvariability3d')
		#Here we just set a label and its position in the window
		title1=QtGui.QLabel('<b>sxvariability3d</b> - convert atomic model (pdb file) into sampled electron density map', self)
		title1.move(self.x1,self.y1)
		self.y1 += 30

		self.repopbtn = QPushButton("Retrieve saved parameters", self)
		self.repopbtn.move(self.x1-5,self.y1)
		#sets an infotip for this Pushbutton
		self.repopbtn.setToolTip('Retrieve saved parameters')
		#when this button is clicked, this action starts the subfunction twodali
		self.connect(self.repopbtn, SIGNAL("clicked()"), self.repoparms_variability3d)

		#######################################################################################
		#Here we create a Button(file_button with title run open .hdf) and its position in the window
		self.file_button = QtGui.QPushButton("Open .pdb", self)
		self.file_button.move(self.x3, self.y2-self.yspc)
		#Here we define, that when this button is clicked, it starts subfunction choose_file
		QtCore.QObject.connect(self.file_button, QtCore.SIGNAL("clicked()"), self.choose_file)
		
		pdb= QtGui.QLabel('Name of pdb file', self)
		pdb.move(self.x1,self.y2)
		self.pdbfileedit=QtGui.QLineEdit(self)
		self.pdbfileedit.move(self.x2,self.y2)
		self.pdbfileedit.setText(self.savedparmsdict['pdbfile'])
		self.y2 += 30
		
		output= QtGui.QLabel('EM Output file', self)
		output.move(self.x1,self.y2)
		self.outputedit=QtGui.QLineEdit(self)
		self.outputedit.move(self.x2,self.y2)
		self.outputedit.setText(self.savedparmsdict['output'])		
		self.outinfobtn = QPushButton("Output Info", self)
		self.outinfobtn.move(self.x3,  self.y2-self.yspc)
		#sets an infotip for this Pushbutton
		self.outinfobtn.setToolTip('Output Info')
		#when this button is clicked, this action starts the subfunction twodali
		self.connect(self.outinfobtn, SIGNAL("clicked()"), self.outputinfo_variability3d)
		self.y2 += 30
		
		apix= QtGui.QLabel('Pixel size (Angstroms)', self)
		apix.move(self.x1,self.y2)
		self.apixedit=QtGui.QLineEdit(self)
		self.apixedit.move(self.x2,self.y2)
		self.apixedit.setText(self.savedparmsdict['apix'])
		self.apixedit.setToolTip('Angstrom/voxel')		
		self.y2 += 30
		
		box= QtGui.QLabel('Box size (voxels)', self)
		box.move(self.x1,self.y2)
		self.boxedit=QtGui.QLineEdit(self)
		self.boxedit.move(self.x2,self.y2)
		self.boxedit.setText(self.savedparmsdict['box'])
		self.boxedit.setToolTip('Box size in pixels, <xyz> or <x,y,z>')		
		self.y2 += 30

		het= QtGui.QLabel('Include HET atoms', self)
		het.move(self.x1,self.y2)
		self.hetchkbx = QtGui.QCheckBox("",self)
		self.hetchkbx.move(self.x2, self.y2)
		self.hetchkbx.setCheckState(self.savedparmsdict['het'])
		self.hetchkbx.setToolTip('Include HET atoms in the map.')		
		self.y2 += 30
		
		center= QtGui.QLabel('Center atomic model', self)
		center.move(self.x1,self.y2)
		self.centeredit=QtGui.QLineEdit(self)
		self.centeredit.move(self.x2,self.y2)
		self.centeredit.setText(self.savedparmsdict['center'])
		self.centeredit.setToolTip('center: c - coordinates; a - center of gravity; <x,y,z> - a vector (in Angstrom) to substract from all coordinates; default: n - no')		
		self.y2 += 30
		
		Och= QtGui.QLabel('Use O system of coordinates', self)
		Och.move(self.x1,self.y2)
		self.Ochchkbx = QtGui.QCheckBox("",self)
		self.Ochchkbx.move(self.x2, self.y2)
		self.Ochchkbx.setCheckState(self.savedparmsdict['Och'])
		self.Ochchkbx.setToolTip('apply additional rotation so the model will appear in O in the same rotation as in chimera.')		
		self.y2 += 30
		
		quiet= QtGui.QLabel('Do not print to monitor', self)
		quiet.move(self.x1,self.y2)
		self.quietchkbx = QtGui.QCheckBox("",self)
		self.quietchkbx.move(self.x2, self.y2)
		self.quietchkbx.setCheckState(self.savedparmsdict['quiet'])
		self.quietchkbx.setToolTip('Verbose is the default')		
		self.y2 += 30
		
		tr0= QtGui.QLabel('Transformation matrix to apply to PDB', self)
		tr0.move(self.x1,self.y2)
		self.tr0edit=QtGui.QLineEdit(self)
		self.tr0edit.move(self.x2,self.y2)
		self.tr0edit.setText(self.savedparmsdict['tr0'])
		self.tr0edit.setToolTip('Filename of initial 3x4 transformation matrix')
		
		self.file_button = QtGui.QPushButton("Open File", self)
		self.file_button.move(self.x3, self.y2-self.yspc)
		#Here we define, that when this button is clicked, it starts subfunction choose_file
		QtCore.QObject.connect(self.file_button, QtCore.SIGNAL("clicked()"), self.choose_file1)
		
		self.y2 += 30
		# make ctf, normalize and init_method radio button...
		
		self.savepbtn = QPushButton("Save Input Parameters", self)
		self.savepbtn.move(self.x1-5, self.y4)
		#sets an infotip for this Pushbutton
		self.savepbtn.setToolTip('Save Input Parameters')
		#when this button is clicked, this action starts the subfunction twodali
		self.connect(self.savepbtn, SIGNAL("clicked()"), self.saveparms)
		self.y4+=30
		
		self.cmdlinebtn = QPushButton("Generate command line from input parameters", self)
		self.cmdlinebtn.move(self.x1-5, self.y4)
		#sets an infotip for this Pushbutton
		self.cmdlinebtn.setToolTip('Generate command line using input parameters')
		#when this button is clicked, this action starts the subfunction twodali
		self.connect(self.cmdlinebtn, SIGNAL("clicked()"), self.gencmdline_variability3d)

		#######################################################################
		 #Here we create a Button(Run_button with title run sxali2d) and its position in the window
		self.RUN_button = QtGui.QPushButton('Run sxvariability3d', self)
		# make 3D textured push button look
		s = "QPushButton {font: bold; color: #000;border: 1px solid #333;border-radius: 11px;padding: 2px;background: qradialgradient(cx: 0, cy: 0,fx: 0.5, fy:0.5,radius: 1, stop: 0 #fff, stop: 1 #8D0);min-width:90px;margin:5px} QPushButton:pressed {font: bold; color: #000;border: 1px solid #333;border-radius: 11px;padding: 2px;background: qradialgradient(cx: 0, cy: 0,fx: 0.5, fy:0.5,radius: 1, stop: 0 #fff, stop: 1 #084);min-width:90px;margin:5px}"
		
		self.RUN_button.setStyleSheet(s)
		self.RUN_button.move(self.x5, self.y5)
		#Here we define, that when this button is clicked, it starts subfunction runsxali2d
		self.connect(self.RUN_button, SIGNAL("clicked()"), self.runsxvariability3d)
		#Labels and Line Edits for User Input		

				
	 #Function runsxali2d started when  the  RUN_button of the  Poptwodali window is clicked 
	def outputinfo_variability3d(self):
		QMessageBox.information(self, "sxvariability3d output",'output 3-D electron density map (any EM format). Attribute pixel_size will be set to the specified value.')
		
	def gencmdline_variability3d(self,writefile=True):
		#Here we just read in all user inputs in the line edits of the Poptwodali window
		pdbfile = self.pdbfileedit.text()
		output=self.outputedit.text()
		apix=self.apixedit.text()
		box=self.boxedit.text()
		center=self.centeredit.text()
		tr0=self.tr0edit.text()
		het=self.hetchkbx.checkState()
		quiet=self.quietchkbx.checkState()
		Och=self.Ochchkbx.checkState()
		
		cmd1 = "sxvariability3d.py "+str(pdbfile) +" "+ str(output)
		
		args = " --apix="+ str(apix)+" --box="+ str(box)+ " --center="+str(center)
		
		cmd1 = cmd1 + args
				
		if het == Qt.Checked:
			cmd1 = cmd1 + " --het"
		if quiet == Qt.Checked:
			cmd1 = cmd1 + " --quiet"
		if Och == Qt.Checked:
			cmd1 = cmd1 + " --O"				
		
		if len(str(tr0)) > 0:
				cmd1 = cmd1 + " --tr0="+str(tr0)
		
		self.savedparmsdict = {'pdbfile':str(pdbfile),'output':str(output),'apix':str(apix),'box':str(box),'het':het,'center':str(center),'Och':Och,'quiet':quiet,'tr0':str(tr0)}
		
		if writefile:		
			(fname,stat)= QInputDialog.getText(self,"Generate Command Line","Enter name of file to save command line in",QLineEdit.Normal,"")
			if stat:
				f = open(fname,'a')
				f.write(cmd1)
				f.write('\n')
				f.close()
		
		print cmd1
		self.cmd = cmd1
		
	def runsxvariability3d(self):
		self.gencmdline_variability3d(writefile=False)
		process = subprocess.Popen(self.cmd,shell=True)
		self.emit(QtCore.SIGNAL("process_started"),process.pid)
		
	def saveparms(self):		
		# save all the parms in a text file so we can repopulate if user requests
		(fname,stat)= QInputDialog.getText(self,"Save Input Parameters","Enter name of file to save parameters in",QLineEdit.Normal,"")
		if stat:
			import pickle
			output=open(fname,'wb')
			self.gencmdline_variability3d(writefile=False)
			pickle.dump(self.savedparmsdict,output)
			output.close()
		
	def repoparms_variability3d(self):		
		(fname,stat)= QInputDialog.getText(self,"Get Input Parameters","Enter name of file parameters were saved in",QLineEdit.Normal,"")
		if stat:
			import pickle
			pkl = open(fname,'rb')
			self.savedparmsdict = pickle.load(pkl)
			self.pdbfileedit.setText(self.savedparmsdict['pdbfile'])
			self.outputedit.setText(self.savedparmsdict['output'])
			self.apixedit.setText(self.savedparmsdict['apix'])		
			self.boxedit.setText(self.savedparmsdict['box'])		
			self.centeredit.setText(self.savedparmsdict['center'])
			self.hetchkbx.setCheckState(self.savedparmsdict['het'])
			self.Ochchkbx.setCheckState(self.savedparmsdict['Och'])
			self.quietchkbx.setCheckState(self.savedparmsdict['quiet'])
			self.tr0edit.setText(self.savedparmsdict['tr0'])
				
	def choose_file(self):
		#opens a file browser, showing files only in .hdf format
		file_name = QtGui.QFileDialog.getOpenFileName(self, "Open PDB File", "", "PDB files (*.pdb)")
		#after the user selected a file, we obtain this filename as a Qstring
		a=QtCore.QString(file_name)
		#we convert this Qstring to a string and send it to line edit classed stackname edit of the Poptwodali window
		self.pdbfileedit.setText(str(a))
	
	
	def choose_file1(self):
		#opens a file browser, showing files only in .hdf format
		file_name = QtGui.QFileDialog.getOpenFileName(self, "Open file containing transformation matrix", "", "(*)")
		#after the user selected a file, we obtain this filename as a Qstring
		a=QtCore.QString(file_name)
		#we convert this Qstring to a string and send it to line edit classed stackname edit of the Poptwodali window
		self.tr0edit.setText(str(a))



class Popuplocres(QWidget):
	def __init__(self):
		QWidget.__init__(self)

		#######################################################################################
		# class variables
		self.cmd = ""
		# populate with default values
		self.savedparmsdict = {'vol1':'','vol2':'','mask':'','output':'','wn':'11','step':'1','cutoff':'0.5','radius':'-1','nproc':'4','fsc':''}
		#######################################################################################
		# Layout parameters
		self.y1 = 10 # title and Repopulate button
		self.y2 = self.y1 + 85 # Text boxes for inputting parameters
		#self.y3 = self.y2 + 222 # activate images button and set xform.align2d button
		self.y4 = self.y2 + 290 # Advanced Parameters, Save Input and Generate command line buttons
		self.y5 = self.y4 + 95 # run button 
		self.yspc = 4
		
		self.x1 = 10 # first column (text box labels)
		self.x2 = self.x1 + 260 # second column (text boxes)
		self.x3 = self.x2+145 # third column (Open .hdf button)
		self.x4 = self.x3+100 # fourth column (Open .bdb button)
		self.x5 = 230 # run button
		#######################################################################################
		
		#Here we just set the window title
		self.setWindowTitle('sxlocres')
		#Here we just set a label and its position in the window
		title1=QtGui.QLabel('<b>sxlocres</b> - compute local FSC resolution in real space unsing half-maps, within area outlined by the maskfile and within regions wn x wn x wn', self)
		title1.move(self.x1,self.y1)
		self.y1 += 30

		self.repopbtn = QPushButton("Retrieve saved parameters", self)
		self.repopbtn.move(self.x1-5,self.y1)
		#sets an infotip for this Pushbutton
		self.repopbtn.setToolTip('Retrieve saved parameters')
		#when this button is clicked, this action starts the subfunction twodali
		self.connect(self.repopbtn, SIGNAL("clicked()"), self.repoparms_locres)

		#######################################################################################
		#Here we create a Button(file_button with title run open .hdf) and its position in the window
		self.file_button = QtGui.QPushButton("Open first half volume", self)
		self.file_button.move(self.x3, self.y2-self.yspc)
		#Here we define, that when this button is clicked, it starts subfunction choose_file
		QtCore.QObject.connect(self.file_button, QtCore.SIGNAL("clicked()"), self.choose_file)
		
		vol1 = QtGui.QLabel('Name of volume 1', self)
		vol1.move(self.x1,self.y2)
		self.vol1edit=QtGui.QLineEdit(self)
		self.vol1edit.move(self.x2,self.y2)
		self.vol1edit.setText(self.savedparmsdict['vol1'])
		self.y2 += 30
		
		vol2 = QtGui.QLabel('Name of volume 2', self)
		vol2.move(self.x1,self.y2)
		self.vol2edit=QtGui.QLineEdit(self)
		self.vol2edit.move(self.x2,self.y2)
		self.vol2edit.setText(self.savedparmsdict['vol2'])
		self.y2 += 30
		
		mask = QtGui.QLabel('Name of 3D mask', self)
		mask.move(self.x1,self.y2)
		self.maskedit=QtGui.QLineEdit(self)
		self.maskedit.move(self.x2,self.y2)
		self.maskedit.setText(self.savedparmsdict['mask'])
		self.y2 += 30

		output = QtGui.QLabel('Output local res volume', self)
		output.move(self.x1,self.y2)
		self.outputedit=QtGui.QLineEdit(self)
		self.outputedit.move(self.x2,self.y2)
		self.outputedit.setText(self.savedparmsdict['output'])		
		self.outinfobtn = QPushButton("Output", self)
		self.outinfobtn.move(self.x3,  self.y2-self.yspc)
		#sets an infotip for this Pushbutton
		self.outinfobtn.setToolTip('Output Info')
		#when this button is clicked, this action starts the subfunction twodali
		self.connect(self.outinfobtn, SIGNAL("clicked()"), self.outputinfo_locres)
		self.y2 += 30
		
		wn = QtGui.QLabel('Window size (voxels)', self)
		wn.move(self.x1,self.y2)
		self.wnedit=QtGui.QLineEdit(self)
		self.wnedit.move(self.x2,self.y2)
		self.wnedit.setText(self.savedparmsdict['wn'])
		self.wnedit.setToolTip('Window size defines local FSC region')		
		self.y2 += 30
		
		step = QtGui.QLabel('Resolution step (voxels)', self)
		step.move(self.x1,self.y2)
		self.stepedit = QtGui.QLineEdit(self)
		self.stepedit.move(self.x2,self.y2)
		self.stepedit.setText(self.savedparmsdict['step'])
		self.stepedit.setToolTip('Resolution step size defines shell width')		
		self.y2 += 30

		cutoff = QtGui.QLabel('Resolution cutoff', self)
		cutoff.move(self.x1,self.y2)
		self.cutoffedit = QtGui.QLineEdit(self)
		self.cutoffedit.move(self.x2,self.y2)
		self.cutoffedit.setText(self.savedparmsdict['cutoff'])
		self.cutoffedit.setToolTip('Resolution cut-off for FSC.  The local resolution volume will contain absolute frequencies corresponding to the selected cut-off values.  Default is 0.5.  Note excesively small values (as 0.14) will result in a choppy local resolution volume.')		
		self.y2 += 30

		nproc = QtGui.QLabel('Number of processors', self)
		nproc.move(self.x1,self.y2)
		self.nprocedit = QtGui.QLineEdit(self)
		self.nprocedit.move(self.x2, self.y2)
		self.nprocedit.setText(self.savedparmsdict['nproc'])
		self.nprocedit.setToolTip('Number of processors for MPI processing.  If set to 1 or 0, a single processor non-MPI version will be used.')		
		self.y2 += 30

		radius = QtGui.QLabel('Radius of the structure', self)
		radius.move(self.x1,self.y2)
		self.radiusedit = QtGui.QLineEdit(self)
		self.radiusedit.move(self.x2, self.y2)
		self.radiusedit.setText(self.savedparmsdict['radius'])
		self.radiusedit.setToolTip('Radius of the structure in pixels.  Only required if mask is not given.')		
		self.y2 += 30

		fsc = QtGui.QLabel('FSC file name', self)
		fsc.move(self.x1,self.y2)
		self.fscedit = QtGui.QLineEdit(self)
		self.fscedit.move(self.x2, self.y2)
		self.fscedit.setText(self.savedparmsdict['fsc'])
		self.fscedit.setToolTip('Optional. Name of the output text file to contain overall 1D resolution curve. It is a rotational average of local FSC resolution values')		
		self.y2 += 30
		"""

		center= QtGui.QLabel('Center atomic model', self)
		center.move(self.x1,self.y2)
		self.centeredit=QtGui.QLineEdit(self)
		self.centeredit.move(self.x2,self.y2)
		self.centeredit.setText(self.savedparmsdict['center'])
		self.centeredit.setToolTip('center: c - coordinates; a - center of gravity; <x,y,z> - a vector (in Angstrom) to substract from all coordinates; default: n - no')		
		self.y2 += 30
		
		Och= QtGui.QLabel('Use O system of coordinates', self)
		Och.move(self.x1,self.y2)
		self.Ochchkbx = QtGui.QCheckBox("",self)
		self.Ochchkbx.move(self.x2, self.y2)
		self.Ochchkbx.setCheckState(self.savedparmsdict['Och'])
		self.Ochchkbx.setToolTip('apply additional rotation so the model will appear in O in the same rotation as in chimera.')		
		self.y2 += 30
		
		quiet= QtGui.QLabel('Do not print to monitor', self)
		quiet.move(self.x1,self.y2)
		self.quietchkbx = QtGui.QCheckBox("",self)
		self.quietchkbx.move(self.x2, self.y2)
		self.quietchkbx.setCheckState(self.savedparmsdict['quiet'])
		self.quietchkbx.setToolTip('Verbose is the default')		
		self.y2 += 30
		tr0= QtGui.QLabel('Transformation matrix to apply to PDB', self)
		tr0.move(self.x1,self.y2)
		self.tr0edit=QtGui.QLineEdit(self)
		self.tr0edit.move(self.x2,self.y2)
		self.tr0edit.setText(self.savedparmsdict['tr0'])
		self.tr0edit.setToolTip('Filename of initial 3x4 transformation matrix')
		
		self.file_button = QtGui.QPushButton("Open File", self)
		self.file_button.move(self.x3, self.y2-self.yspc)
		#Here we define, that when this button is clicked, it starts subfunction choose_file
		QtCore.QObject.connect(self.file_button, QtCore.SIGNAL("clicked()"), self.choose_file1)
		
		self.y2 += 30
		"""
		
		self.savepbtn = QPushButton("Save Input Parameters", self)
		self.savepbtn.move(self.x1-5, self.y4)
		#sets an infotip for this Pushbutton
		self.savepbtn.setToolTip('Save Input Parameters')
		#when this button is clicked, this action starts the subfunction twodali
		self.connect(self.savepbtn, SIGNAL("clicked()"), self.saveparms)
		self.y4+=30
		
		self.cmdlinebtn = QPushButton("Generate command line from input parameters", self)
		self.cmdlinebtn.move(self.x1-5, self.y4)
		#sets an infotip for this Pushbutton
		self.cmdlinebtn.setToolTip('Generate command line using input parameters')
		#when this button is clicked, this action starts the subfunction twodali
		self.connect(self.cmdlinebtn, SIGNAL("clicked()"), self.gencmdline_locres)

		#######################################################################
		#Here we create a Button(Run_button with title run sxali2d) and its position in the window
		self.RUN_button = QtGui.QPushButton('Run sxlocres', self)
		# make 3D textured push button look
		s = "QPushButton {font: bold; color: #000;border: 1px solid #333;border-radius: 11px;padding: 2px;background: qradialgradient(cx: 0, cy: 0,fx: 0.5, fy:0.5,radius: 1, stop: 0 #fff, stop: 1 #8D0);min-width:90px;margin:5px} QPushButton:pressed {font: bold; color: #000;border: 1px solid #333;border-radius: 11px;padding: 2px;background: qradialgradient(cx: 0, cy: 0,fx: 0.5, fy:0.5,radius: 1, stop: 0 #fff, stop: 1 #084);min-width:90px;margin:5px}"
		
		self.RUN_button.setStyleSheet(s)
		self.RUN_button.move(self.x5, self.y5)
		#Here we define, that when this button is clicked, it starts subfunction runsxali2d
		self.connect(self.RUN_button, SIGNAL("clicked()"), self.runsxlocres)
		#Labels and Line Edits for User Input		

				
	 #Function runsxali2d started when  the  RUN_button of the  Poptwodali window is clicked 
	def outputinfo_locres(self):
		QMessageBox.information(self, "sxlocres output",'output 3D local resolution map (any EM format).')
		
	def gencmdline_locres(self,writefile=True):
		#Here we just read in all user inputs in the line edits of the Poptwodali window
		vol1     = self.vol1edit.text()
		vol2     = self.vol2edit.text()
		mask     = self.maskedit.text()
		output   = self.outputedit.text()
		wn       = self.wnedit.text()
		step     = self.stepedit.text()
		cutoff   = self.cutoffedit.text()
		radius   = self.radiusedit.text()
		fsc      = self.fscedit.text()
		nproc    = self.nprocedit.text()

		cmd1 = "sxlocres.py "+str(vol1) +" "+ str(vol2)
		if(len(str(mask)) > 0): cmd1 += " "+ str(mask)
		cmd1 += " "+ str(output)
		if int(str(nproc)) >1:
				cmd1 = "mpirun -np " + str(nproc) + " " + cmd1
		cmd1 += " --wn="+ str(wn)+" --step="+ str(step) + " --cutoff="+str(cutoff)+ " --radius="+str(radius)+ " --fsc="+str(fsc)
		if int(str(nproc)) >1:
			cmd1 += "  --MPI"

		
		self.savedparmsdict = {'vol1':str(vol1),'vol2':str(vol2),'mask':str(mask),'output':str(output),'wn':str(wn),\
						'step':str(step),'cutoff':cutoff,'radius':str(radius),'fsc':str(fsc),'nproc':nproc}
		
		
		if writefile:		
			(fname,stat)= QInputDialog.getText(self,"Generate Command Line","Enter name of file to save command line in",QLineEdit.Normal,"")
			if stat:
				f = open(fname,'a')
				f.write(cmd1)
				f.write('\n')
				f.close()
		
		print cmd1
		self.cmd = cmd1
		
	def runsxlocres(self):
		self.gencmdline_locres(writefile=False)
		process = subprocess.Popen(self.cmd,shell=True)
		self.emit(QtCore.SIGNAL("process_started"),process.pid)
		
	def saveparms(self):		
		# save all the parms in a text file so we can repopulate if user requests
		(fname,stat)= QInputDialog.getText(self,"Save Input Parameters","Enter name of file to save parameters in",QLineEdit.Normal,"")
		if stat:
			import pickle
			output=open(fname,'wb')
			self.gencmdline_locres(writefile=False)
			pickle.dump(self.savedparmsdict,output)
			output.close()
		
	def repoparms_locres(self):		
		(fname,stat)= QInputDialog.getText(self,"Get Input Parameters","Enter name of file parameters were saved in",QLineEdit.Normal,"")
		if stat:
			import pickle
			pkl = open(fname,'rb')
			self.savedparmsdict = pickle.load(pkl)
			self.vol1edit.setText(self.savedparmsdict['vol1'])
			self.vol2edit.setText(self.savedparmsdict['vol2'])
			self.maskedit.setText(self.savedparmsdict['mask'])
			self.outputedit.setText(self.savedparmsdict['output'])
			self.wnedit.setText(self.savedparmsdict['wn'])		
			self.stepedit.setText(self.savedparmsdict['step'])		
			self.cutoffedit.setText(self.savedparmsdict['cutoff'])
			self.radiusedit.setText(self.savedparmsdict['radius'])
			self.fscedit.setText(self.savedparmsdict['fsc'])
			self.nprocedit.setText(self.savedparmsdict['nproc'])

				
	def choose_file(self):
		#opens a file browser, showing files only in .hdf format
		file_name = QtGui.QFileDialog.getOpenFileName(self, "Open volume", "", "PDB files (*.pdb)")
		#after the user selected a file, we obtain this filename as a Qstring
		a=QtCore.QString(file_name)
		#we convert this Qstring to a string and send it to line edit classed stackname edit of the Poptwodali window
		self.vol1edit.setText(str(a))

	def choose_file1(self):
		#opens a file browser, showing files only in .hdf format
		file_name = QtGui.QFileDialog.getOpenFileName(self, "Open file containing transformation matrix", "", "(*)")
		#after the user selected a file, we obtain this filename as a Qstring
		a=QtCore.QString(file_name)
		#we convert this Qstring to a string and send it to line edit classed stackname edit of the Poptwodali window
		self.tr0edit.setText(str(a))




class Popupfilterlocal(QWidget):
	def __init__(self):
		QWidget.__init__(self)
		
		#######################################################################################
		# class variables
		self.cmd = ""
		# populate with default values
		self.savedparmsdict = {'pdbfile':'','output':'','apix':'1.0','box':'150','het':Qt.Unchecked,'center':'n','Och':Qt.Unchecked,'quiet':Qt.Unchecked,'tr0':''}
		if (options.demo == DEMO_mpibdbctf) or (options.demo == DEMO_mpibdb):
			self.savedparmsdict['pdbfile'] = '../tteftu_with_tRNA.pdb'
			self.savedparmsdict['output']='tmp.hdf'
			self.savedparmsdict['apix']='5.2'
			self.savedparmsdict['box']='64'
		#######################################################################################
		# Layout parameters
		
		self.y1 = 10 # title and Repopulate button
		self.y2 = self.y1 + 85 # Text boxes for inputting parameters
		#self.y3 = self.y2 + 222 # activate images button and set xform.align2d button
		self.y4 = self.y2 + 290 # Advanced Parameters, Save Input and Generate command line buttons
		self.y5 = self.y4 + 95 # run button 
		self.yspc = 4
		
		self.x1 = 10 # first column (text box labels)
		self.x2 = self.x1 + 260 # second column (text boxes)
		self.x3 = self.x2+145 # third column (Open .hdf button)
		self.x4 = self.x3+100 # fourth column (Open .bdb button)
		self.x5 = 230 # run button
		#######################################################################################
		
		#Here we just set the window title
		self.setWindowTitle('sxfilterlocal')
		#Here we just set a label and its position in the window
		title1=QtGui.QLabel('<b>sxfilterlocal</b> - convert atomic model (pdb file) into sampled electron density map', self)
		title1.move(self.x1,self.y1)
		self.y1 += 30

		self.repopbtn = QPushButton("Retrieve saved parameters", self)
		self.repopbtn.move(self.x1-5,self.y1)
		#sets an infotip for this Pushbutton
		self.repopbtn.setToolTip('Retrieve saved parameters')
		#when this button is clicked, this action starts the subfunction twodali
		self.connect(self.repopbtn, SIGNAL("clicked()"), self.repoparms_filterlocal)

		#######################################################################################
		#Here we create a Button(file_button with title run open .hdf) and its position in the window
		self.file_button = QtGui.QPushButton("Open .pdb", self)
		self.file_button.move(self.x3, self.y2-self.yspc)
		#Here we define, that when this button is clicked, it starts subfunction choose_file
		QtCore.QObject.connect(self.file_button, QtCore.SIGNAL("clicked()"), self.choose_file)
		
		pdb= QtGui.QLabel('Name of pdb file', self)
		pdb.move(self.x1,self.y2)
		self.pdbfileedit=QtGui.QLineEdit(self)
		self.pdbfileedit.move(self.x2,self.y2)
		self.pdbfileedit.setText(self.savedparmsdict['pdbfile'])
		self.y2 += 30
		
		output= QtGui.QLabel('EM Output file', self)
		output.move(self.x1,self.y2)
		self.outputedit=QtGui.QLineEdit(self)
		self.outputedit.move(self.x2,self.y2)
		self.outputedit.setText(self.savedparmsdict['output'])		
		self.outinfobtn = QPushButton("Output Info", self)
		self.outinfobtn.move(self.x3,  self.y2-self.yspc)
		#sets an infotip for this Pushbutton
		self.outinfobtn.setToolTip('Output Info')
		#when this button is clicked, this action starts the subfunction twodali
		self.connect(self.outinfobtn, SIGNAL("clicked()"), self.outputinfo_filterlocal)
		self.y2 += 30
		
		apix= QtGui.QLabel('Pixel size (Angstroms)', self)
		apix.move(self.x1,self.y2)
		self.apixedit=QtGui.QLineEdit(self)
		self.apixedit.move(self.x2,self.y2)
		self.apixedit.setText(self.savedparmsdict['apix'])
		self.apixedit.setToolTip('Angstrom/voxel')		
		self.y2 += 30
		
		box= QtGui.QLabel('Box size (voxels)', self)
		box.move(self.x1,self.y2)
		self.boxedit=QtGui.QLineEdit(self)
		self.boxedit.move(self.x2,self.y2)
		self.boxedit.setText(self.savedparmsdict['box'])
		self.boxedit.setToolTip('Box size in pixels, <xyz> or <x,y,z>')		
		self.y2 += 30

		het= QtGui.QLabel('Include HET atoms', self)
		het.move(self.x1,self.y2)
		self.hetchkbx = QtGui.QCheckBox("",self)
		self.hetchkbx.move(self.x2, self.y2)
		self.hetchkbx.setCheckState(self.savedparmsdict['het'])
		self.hetchkbx.setToolTip('Include HET atoms in the map.')		
		self.y2 += 30
		
		center= QtGui.QLabel('Center atomic model', self)
		center.move(self.x1,self.y2)
		self.centeredit=QtGui.QLineEdit(self)
		self.centeredit.move(self.x2,self.y2)
		self.centeredit.setText(self.savedparmsdict['center'])
		self.centeredit.setToolTip('center: c - coordinates; a - center of gravity; <x,y,z> - a vector (in Angstrom) to substract from all coordinates; default: n - no')		
		self.y2 += 30
		
		Och= QtGui.QLabel('Use O system of coordinates', self)
		Och.move(self.x1,self.y2)
		self.Ochchkbx = QtGui.QCheckBox("",self)
		self.Ochchkbx.move(self.x2, self.y2)
		self.Ochchkbx.setCheckState(self.savedparmsdict['Och'])
		self.Ochchkbx.setToolTip('apply additional rotation so the model will appear in O in the same rotation as in chimera.')		
		self.y2 += 30
		
		quiet= QtGui.QLabel('Do not print to monitor', self)
		quiet.move(self.x1,self.y2)
		self.quietchkbx = QtGui.QCheckBox("",self)
		self.quietchkbx.move(self.x2, self.y2)
		self.quietchkbx.setCheckState(self.savedparmsdict['quiet'])
		self.quietchkbx.setToolTip('Verbose is the default')		
		self.y2 += 30
		
		tr0= QtGui.QLabel('Transformation matrix to apply to PDB', self)
		tr0.move(self.x1,self.y2)
		self.tr0edit=QtGui.QLineEdit(self)
		self.tr0edit.move(self.x2,self.y2)
		self.tr0edit.setText(self.savedparmsdict['tr0'])
		self.tr0edit.setToolTip('Filename of initial 3x4 transformation matrix')
		
		self.file_button = QtGui.QPushButton("Open File", self)
		self.file_button.move(self.x3, self.y2-self.yspc)
		#Here we define, that when this button is clicked, it starts subfunction choose_file
		QtCore.QObject.connect(self.file_button, QtCore.SIGNAL("clicked()"), self.choose_file1)
		
		self.y2 += 30
		# make ctf, normalize and init_method radio button...
		
		self.savepbtn = QPushButton("Save Input Parameters", self)
		self.savepbtn.move(self.x1-5, self.y4)
		#sets an infotip for this Pushbutton
		self.savepbtn.setToolTip('Save Input Parameters')
		#when this button is clicked, this action starts the subfunction twodali
		self.connect(self.savepbtn, SIGNAL("clicked()"), self.saveparms)
		self.y4+=30
		
		self.cmdlinebtn = QPushButton("Generate command line from input parameters", self)
		self.cmdlinebtn.move(self.x1-5, self.y4)
		#sets an infotip for this Pushbutton
		self.cmdlinebtn.setToolTip('Generate command line using input parameters')
		#when this button is clicked, this action starts the subfunction twodali
		self.connect(self.cmdlinebtn, SIGNAL("clicked()"), self.gencmdline_filterlocal)

		#######################################################################
		 #Here we create a Button(Run_button with title run sxali2d) and its position in the window
		self.RUN_button = QtGui.QPushButton('Run sxfilterlocal', self)
		# make 3D textured push button look
		s = "QPushButton {font: bold; color: #000;border: 1px solid #333;border-radius: 11px;padding: 2px;background: qradialgradient(cx: 0, cy: 0,fx: 0.5, fy:0.5,radius: 1, stop: 0 #fff, stop: 1 #8D0);min-width:90px;margin:5px} QPushButton:pressed {font: bold; color: #000;border: 1px solid #333;border-radius: 11px;padding: 2px;background: qradialgradient(cx: 0, cy: 0,fx: 0.5, fy:0.5,radius: 1, stop: 0 #fff, stop: 1 #084);min-width:90px;margin:5px}"
		
		self.RUN_button.setStyleSheet(s)
		self.RUN_button.move(self.x5, self.y5)
		#Here we define, that when this button is clicked, it starts subfunction runsxali2d
		self.connect(self.RUN_button, SIGNAL("clicked()"), self.runsxfilterlocal)
		#Labels and Line Edits for User Input		

				
	 #Function runsxali2d started when  the  RUN_button of the  Poptwodali window is clicked 
	def outputinfo_filterlocal(self):
		QMessageBox.information(self, "sxfilterlocal output",'output 3-D electron density map (any EM format). Attribute pixel_size will be set to the specified value.')
		
	def gencmdline_filterlocal(self,writefile=True):
		#Here we just read in all user inputs in the line edits of the Poptwodali window
		pdbfile = self.pdbfileedit.text()
		output=self.outputedit.text()
		apix=self.apixedit.text()
		box=self.boxedit.text()
		center=self.centeredit.text()
		tr0=self.tr0edit.text()
		het=self.hetchkbx.checkState()
		quiet=self.quietchkbx.checkState()
		Och=self.Ochchkbx.checkState()
		
		cmd1 = "sxfilterlocal.py "+str(pdbfile) +" "+ str(output)
		
		args = " --apix="+ str(apix)+" --box="+ str(box)+ " --center="+str(center)
		
		cmd1 = cmd1 + args
				
		if het == Qt.Checked:
			cmd1 = cmd1 + " --het"
		if quiet == Qt.Checked:
			cmd1 = cmd1 + " --quiet"
		if Och == Qt.Checked:
			cmd1 = cmd1 + " --O"				
		
		if len(str(tr0)) > 0:
			cmd1 = cmd1 + " --tr0="+str(tr0)
		
		self.savedparmsdict = {'pdbfile':str(pdbfile),'output':str(output),'apix':str(apix),'box':str(box),'het':het,'center':str(center),'Och':Och,'quiet':quiet,'tr0':str(tr0)}
		
		
		if writefile:		
			(fname,stat)= QInputDialog.getText(self,"Generate Command Line","Enter name of file to save command line in",QLineEdit.Normal,"")
			if stat:
				f = open(fname,'a')
				f.write(cmd1)
				f.write('\n')
				f.close()
		
		print cmd1
		self.cmd = cmd1
		
	def runsxfilterlocal(self):
		self.gencmdline_filterlocal(writefile=False)
		process = subprocess.Popen(self.cmd,shell=True)
		self.emit(QtCore.SIGNAL("process_started"),process.pid)
		
	def saveparms(self):		
		# save all the parms in a text file so we can repopulate if user requests
		(fname,stat)= QInputDialog.getText(self,"Save Input Parameters","Enter name of file to save parameters in",QLineEdit.Normal,"")
		if stat:
			import pickle
			output=open(fname,'wb')
			self.gencmdline_filterlocal(writefile=False)
			pickle.dump(self.savedparmsdict,output)
			output.close()
		
	def repoparms_filterlocal(self):		
		(fname,stat)= QInputDialog.getText(self,"Get Input Parameters","Enter name of file parameters were saved in",QLineEdit.Normal,"")
		if stat:
			import pickle
			pkl = open(fname,'rb')
			self.savedparmsdict = pickle.load(pkl)
			self.pdbfileedit.setText(self.savedparmsdict['pdbfile'])
			self.outputedit.setText(self.savedparmsdict['output'])
			self.apixedit.setText(self.savedparmsdict['apix'])		
			self.boxedit.setText(self.savedparmsdict['box'])		
			self.centeredit.setText(self.savedparmsdict['center'])
			self.hetchkbx.setCheckState(self.savedparmsdict['het'])
			self.Ochchkbx.setCheckState(self.savedparmsdict['Och'])
			self.quietchkbx.setCheckState(self.savedparmsdict['quiet'])
			self.tr0edit.setText(self.savedparmsdict['tr0'])

	def choose_file(self):
		#opens a file browser, showing files only in .hdf format
		file_name = QtGui.QFileDialog.getOpenFileName(self, "Open PDB File", "", "PDB files (*.pdb)")
		#after the user selected a file, we obtain this filename as a Qstring
		a=QtCore.QString(file_name)
		#we convert this Qstring to a string and send it to line edit classed stackname edit of the Poptwodali window
		self.pdbfileedit.setText(str(a))
	
	
		def choose_file1(self):
			#opens a file browser, showing files only in .hdf format
			file_name = QtGui.QFileDialog.getOpenFileName(self, "Open file containing transformation matrix", "", "(*)")
			#after the user selected a file, we obtain this filename as a Qstring
			a=QtCore.QString(file_name)
			#we convert this Qstring to a string and send it to line edit classed stackname edit of the Poptwodali window
			self.tr0edit.setText(str(a))



###MAIN WINDOW		(started by class App)
#This class includes the layout of the main window; within each class, i name the main object self, to avoid confusion)			
class MainWindow(QtGui.QWidget):
	def __init__(self, parent=None):
		QtGui.QWidget.__init__(self, parent)
		# self.setStyleSheet('background-image: url("1.png")')
		#sets the title of the window
		self.setWindowTitle('SPARX GUI')
		self.setAutoFillBackground(True)				
		palette = QPalette(self)				
		palette.setBrush(QPalette.Background, QBrush(QPixmap(get_image_directory()+"sxgui.py_main_window_background_image.png")))				
		# palette.setBrush(QPalette.Background, QBrush(QPixmap("Fig6.png")))
		# palette.setBrush(QPalette.Background, QBrush(QPixmap("spaxgui02.png")))
		self.setPalette(palette)				


		self.y2 = 65

		self.btn7 = QPushButton("sxpdb2em", self)
		self.btn7.move(10, self.y2)
		#sets an infotip for this Pushbutton
		self.btn7.setToolTip('Convert atomic model (pdb file) into sampled electron density map')
		self.connect(self.btn7, SIGNAL("clicked()"), self.pdb2em)
		"""
		self.y2 += 30

		#creates a Pushbutton, named sxali2d, defines its position in the window 
		self.btn1 = QPushButton("sxali2d", self)
		self.btn1.move(10, self.y2)
		#sets an infotip for this Pushbutton
		self.btn1.setToolTip('2D  reference free alignment of an image series ')
		#when this button is clicked, this action starts the subfunction twodali
		self.connect(self.btn1, SIGNAL("clicked()"), self.twodali)

		self.y2 += 30

		#another Pushbutton, with a tooltip, not linked to a function yet
		self.btn2 = QPushButton("sxmref_ali2d", self)
		self.btn2.setToolTip('2D MULTI-referencealignment of an image series ')
		self.btn2.move(10, self.y2)


		self.y2 += 30

		self.btn6 = QPushButton("sxk_means_groups", self)
		self.btn6.move(10, self.y2)
		#sets an infotip for this Pushbutton
		self.btn6.setToolTip('determine \'best\' number of clusters in the data using K-means classification of a set of images')
		self.connect(self.btn6, SIGNAL("clicked()"), self.kmeansgroups)

		self.y2 += 30

		self.btn5 = QPushButton("sxk_means", self)
		self.btn5.move(10, self.y2)
		#sets an infotip for this Pushbutton
		self.btn5.setToolTip('K-means classification of a set of images')
		self.connect(self.btn5, SIGNAL("clicked()"), self.kmeans)

		self.y2 += 30


		self.btn8 = QPushButton("sxpca", self)
		self.btn8.move(10, self.y2)
		#sets an infotip for this Pushbutton
		self.btn8.setToolTip('Principal Component Analysis of images')
		self.connect(self.btn8, SIGNAL("clicked()"), self.pca)
		"""		
		self.y2 += 30

		self.btn10 = QPushButton("sxcter", self)
		self.btn10.move(10, self.y2)
		#sets an infotip for this Pushbutton
		self.btn10.setToolTip('Automated estimation of CTF parameters with error assessment')
		self.connect(self.btn10, SIGNAL("clicked()"), self.cter)

		self.y2 += 30

		self.btn10 = QPushButton("sxisac", self)
		self.btn10.move(10, self.y2)
		#sets an infotip for this Pushbutton
		self.btn10.setToolTip('Perform Iterative Stable Alignment and Clustering (ISAC) on a 2-D image stack')
		self.connect(self.btn10, SIGNAL("clicked()"), self.isac)

		self.y2 += 30

		self.btn10 = QPushButton("sxviper", self)
		self.btn10.move(10, self.y2)
		#sets an infotip for this Pushbutton
		self.btn10.setToolTip('Use the "common line" algorithm to assign initial values of phi, theta, psi to 2D average projections.')
		# self.btn10.mouseDoubleClickEvent(lambda : os.system("python -m webbrowser http://sparx-em.org/sparxwiki/sxviper"))
		self.connect(self.btn10, SIGNAL("clicked()"), self.viper)
		# def mmm(x):
		#		 os.system("python -m webbrowser http://sparx-em.org/sparxwiki/sxviper")
		# self.btn10.mouseDoubleClickEvent = mmm



		self.y2 += 30

		"""		
		self.y2 += 30

		self.btn4 = QPushButton("sxali3d", self)
		self.btn4.move(10, self.y2)
		#sets an infotip for this Pushbutton
		self.btn4.setToolTip('Perform 3-D projection matching given initial reference volume and image series')
		self.connect(self.btn4, SIGNAL("clicked()"), self.ali3d)
		"""

		self.btn4 = QPushButton("sxmeridien", self)
		self.btn4.move(10, self.y2)
		#sets an infotip for this Pushbutton
		self.btn4.setToolTip('3D structure refinement')
		self.connect(self.btn4, SIGNAL("clicked()"), self.meridien)

		self.y2 += 30

		self.btn4 = QPushButton("sx3dvariability", self)
		self.btn4.move(10, self.y2)
		#sets an infotip for this Pushbutton
		self.btn4.setToolTip('3D local variability (real space)')
		self.connect(self.btn4, SIGNAL("clicked()"), self.variability3d)

		self.y2 += 30

		self.btn4 = QPushButton("sxlocres", self)
		self.btn4.move(10, self.y2)
		#sets an infotip for this Pushbutton
		self.btn4.setToolTip('Local resolution (FSC)')
		self.connect(self.btn4, SIGNAL("clicked()"), self.locres)

		self.y2 += 30

		self.btn4 = QPushButton("sxfilterlocal", self)
		self.btn4.move(10, self.y2)
		#sets an infotip for this Pushbutton
		self.btn4.setToolTip('3D local filter')
		self.connect(self.btn4, SIGNAL("clicked()"), self.filterlocal)

		self.y2 += 30
		"""
		self.btn9 = QPushButton("sxmref_ali3d", self)
		self.btn9.move(10, self.y2)
		#sets an infotip for this Pushbutton
		self.btn9.setToolTip('perform 3-D multireference projection matching given initial reference volumes and image series')
		self.connect(self.btn9, SIGNAL("clicked()"), self.mref_ali3d)

		self.y2 += 30

		self.btn11 = QPushButton("sxlocal_ali3d", self)
		self.btn11.move(10, self.y2)
		#sets an infotip for this Pushbutton
		self.btn11.setToolTip('Perform local refinement of 3-D projection alignment of image series using highy accurate gridding method')
		self.connect(self.btn11, SIGNAL("clicked()"), self.localali3d)
		"""

		self.btn11 = QPushButton("sxsort3d", self)
		self.btn11.move(10, self.y2)
		#sets an infotip for this Pushbutton
		self.btn11.setToolTip('Sorts out possible conformations from one heterogenous data set whose xform.projection parameters are already determined using K-means, and Equal K-means method.')
		# self.connect(self.btn11, SIGNAL("clicked()"), self.localali3d)

		self.y2 += 30

		self.btn3 = QPushButton("sxhelical", self)
		self.btn3.move(10, self.y2)
		#sets an infotip for this Pushbutton
		self.btn3.setToolTip('Iterative Real Space Helical Refinement ')
		#when this button is clicked, this action starts the subfunction twodali
		self.connect(self.btn3, SIGNAL("clicked()"), self.helicalrefinement)


		#Pushbutton named Info
		self.picbutton = QPushButton(self)
		#when this button is clicked, this action starts the subfunction info
		self.connect(self.picbutton, SIGNAL("clicked()"), self.info)
		#creates a Pushbutton, named sxhelical defines its position in the window 
		#this decorates the button with the sparx image
		icon = QIcon(get_image_directory()+"sparxicon.png")
		self.picbutton.setIcon(icon)
		self.picbutton.move(5, 5)
		self.picbutton.setToolTip('Info Page')
		#Quitbutton
		self.btn3 = QPushButton("Close", self)
		self.btn3.setToolTip('Close SPARX GUI ')
		self.btn3.move(65, 5)
		self.connect(self.btn3, QtCore.SIGNAL('clicked()'),QtGui.qApp, QtCore.SLOT('quit()'))
		#here we set two labels, their position and font style
		title=QtGui.QLabel("<span style='font-size:18pt; font-weight:600; color:#aa0000;'><b>PROGRAMS </b></span><span style='font-size:12pt; font-weight:60; color:#aa0000;'>(shift-click for wiki)</span>", self)
		title.move(17,47)
		QtGui.QToolTip.setFont(QtGui.QFont('OldEnglish', 8))

		#here we set the width and height of the main window
		self.resize(300,400)


	#This is the function two2ali, which is being started when the Pushbutton btn1 of the main window(called sxali2d) is being clicked
	def twodali(self):
		#print "Opening a new popup window..."
		#opens the window Poptwodali, and defines its width and height
		#The layout of the Poptwodali window is defined in class Poptwodali(QWidget Window)
		self.w = Popuptwodali()
		self.w1 = Popupadvparams_ali2d(self.w.savedparmsdict)
		self.w.w1 = self.w1
		self.TabWidget = QtGui.QTabWidget()
		self.TabWidget.insertTab(0,self.w,'Main')
		self.TabWidget.insertTab(1,self.w1,'Advanced')
		self.TabWidget.resize(550,570)
		self.TabWidget.show()
	
	def helicalrefinement(self):
		#print "Opening a new popup window..."
		#opens the window Poptwodali, and defines its width and height
		#self.w.show()
		self.w = PopupHelicalRefinement()
		self.w1 = Popupadvparams_helical_1(self.w.savedparmsdict)
		self.w2 = Popupadvparams_helical_2(self.w.savedparmsdict)

		intro_string = "Place holder text....fill this in"
		self.w3 = Popupcenter(self.w,intro_string)

		self.w.w1 = self.w1
		self.w.w2 = self.w2
		self.w.w3 = self.w3
		self.TabWidget = QtGui.QTabWidget()
		self.TabWidget.insertTab(0,self.w,'Main')
		self.TabWidget.insertTab(1,self.w1,'Advanced CTF and Search')
		self.TabWidget.insertTab(2,self.w2,'Advanced Symmetry')
		self.TabWidget.insertTab(3,self.w3,'Pre-center input stack (Recommended)')
		self.TabWidget.resize(730,800)
		self.TabWidget.show()
	

	def ali3d(self):
			#print "Opening a new popup window..."
			#opens the window Poptwodali, and defines its width and height
			#The layout of the Poptwodali window is defined in class Poptwodali(QWidget Window)
			self.w = Popupthreedali()
			self.w1 = Popupadvparams_ali3d_1(self.w.savedparmsdict)
			self.w2 = Popupadvparams_ali3d_2(self.w.savedparmsdict)
			self.w.w1 = self.w1
			self.w.w2 = self.w2
			self.TabWidget = QtGui.QTabWidget()
			self.TabWidget.insertTab(0,self.w,'Main')
			self.TabWidget.insertTab(1,self.w1,'Advanced CTF and Search')
			self.TabWidget.insertTab(2,self.w2,'Advanced MISC')
			self.TabWidget.resize(550,650)
			self.TabWidget.show()


	def localali3d(self):
			#print "Opening a new popup window..."
			#opens the window Poptwodali, and defines its width and height
			#The layout of the Poptwodali window is defined in class Poptwodali(QWidget Window)
			self.w = Popuplocalthreedali()
			self.w1 = Popupadvparams_localali3d_1(self.w.savedparmsdict)
			self.w.w1 = self.w1
			self.TabWidget = QtGui.QTabWidget()
			self.TabWidget.insertTab(0,self.w,'Main')
			self.TabWidget.insertTab(1,self.w1,'Advanced CTF')
			self.TabWidget.resize(550,600)
			self.TabWidget.show()		 

	def kmeans(self):
			#print "Opening a new popup window..."
			#opens the window Poptwodali, and defines its width and height
			#The layout of the Poptwodali window is defined in class Poptwodali(QWidget Window)
			self.w = Popupkmeans()
			self.w1 = Popupadvparams_kmeans(self.w.savedparmsdict)
			self.w.w1 = self.w1
			self.TabWidget = QtGui.QTabWidget()
			self.TabWidget.insertTab(0,self.w,'Main')
			self.TabWidget.insertTab(1,self.w1,'Advanced')
			self.TabWidget.resize(580,570)
			self.TabWidget.show()
	"""
def kmeansgroups(self):
			#print "Opening a new popup window..."		
			self.w = Popupkmeansgroups()
			self.w1 = Popupadvparams_kmeans_groups(self.w.savedparmsdict)
			self.w.w1 = self.w1
			self.TabWidget = QtGui.QTabWidget()
			self.TabWidget.insertTab(0,self.w,'Main')
			self.TabWidget.insertTab(1,self.w1,'Advanced')
			self.TabWidget.resize(580,570)
			self.TabWidget.show()
	"""
	def pdb2em(self):
			#print "Opening a new popup window..."
			#opens the window Poptwodali, and defines its width and height
			#The layout of the Poptwodali window is defined in class Poptwodali(QWidget Window)
			self.w = Popuppdb2em()
			self.w.resize(580,570)
			self.w.show()
	
	def pca(self):
			#print "Opening a new popup window..."
			#opens the window Poptwodali, and defines its width and height
			#The layout of the Poptwodali window is defined in class Poptwodali(QWidget Window)
			self.w = Popuppca()
			self.w1 = Popupadvparams_pca(self.w.savedparmsdict)
			self.w.w1 = self.w1
			self.TabWidget = QtGui.QTabWidget()
			self.TabWidget.insertTab(0,self.w,'Main')
			self.TabWidget.insertTab(1,self.w1,'Advanced')
			self.TabWidget.resize(630,500)
			self.TabWidget.show()
	
	def mref_ali3d(self):
			#print "Opening a new popup window..."
			#opens the window Poptwodali, and defines its width and height
			#The layout of the Poptwodali window is defined in class Poptwodali(QWidget Window)
			self.w = Popupmrefthreedali()
			self.w1 = Popupadvparams_mrefali3d_1(self.w.savedparmsdict)
			self.w.w1 = self.w1
			self.TabWidget = QtGui.QTabWidget()
			self.TabWidget.insertTab(0,self.w,'Main')
			self.TabWidget.insertTab(1,self.w1,'Advanced Parameters')
			self.TabWidget.resize(570,740)
			self.TabWidget.show()

	def cter(self):
			##print "Opening a new popup window..."
			#opens the window Poptwodali, and defines its width and height
			#The layout of the Poptwodali window is defined in class Poptwodali(QWidget Window)
			self.w = Popupisac()
			self.w1 = Popupadvparams_cter_1(self.w.savedparmsdict)
			#intro_string = "Before running sxisac.py, it is recommended that the stack be centered. The centering is performed \nusing sxshftali.py. The alignment parameters calculated by the centering procedure is stored in the \nheaders of the input stack as xform.align2d. \n\nTo apply orientation parameters stored in the file headers, check the 'Apply calculated centering parameters \nto input stack' box below and enter the name of the output stack. The resulting output stack will be the input \nstack after applying the shifts calculated by the centering procedure. The orientation parameters 'sx' and \n'sy' in the header of the transformed stack will be set to 0."
			#self.w2 = Popupcenter(self.w,intro_string)
			self.w.w1 = self.w1
			#self.w.w2 = self.w2
			self.TabWidget = QtGui.QTabWidget()
			self.TabWidget.insertTab(0,self.w,'Main')
			self.TabWidget.insertTab(1,self.w1,'Advanced')
			#self.TabWidget.insertTab(2,self.w2,'Pre-center input stack (Recommended)')
			self.TabWidget.resize(730,800)
			self.TabWidget.show()

	def isac(self):
		modifiers = QtGui.QApplication.keyboardModifiers()
		if modifiers == QtCore.Qt.ShiftModifier:
			os.system("python -m webbrowser %ssxisac"%SPARX_DOCUMENTATION_WEBSITE)
			return
		##print "Opening a new popup window..."
		#opens the window Poptwodali, and defines its width and height
		#The layout of the Poptwodali window is defined in class Poptwodali(QWidget Window)
		self.w = Popupisac()
		# self.w1 = Popupadvparams_isac_1(self.w.savedparmsdict)
		self.w1 = Popupadvparams_isac(self.w.savedparmsdict)
		intro_string = "Before running sxisac.py, it is recommended that the stack be centered. The centering is performed \nusing sxshftali.py. The alignment parameters calculated by the centering procedure is stored in the \nheaders of the input stack as xform.align2d. \n\nTo apply orientation parameters stored in the file headers, check the 'Apply calculated centering parameters \nto input stack' box below and enter the name of the output stack. The resulting output stack will be the input \nstack after applying the shifts calculated by the centering procedure. The orientation parameters 'sx' and \n'sy' in the header of the transformed stack will be set to 0."
		# self.w2 = Popupcenter(self.w,intro_string)
		self.w.w1 = self.w1
		# self.w.w2 = self.w2
		self.TabWidget = QtGui.QTabWidget()
		self.TabWidget.insertTab(0,self.w,'Main')
		self.TabWidget.insertTab(1,self.w1,'Advanced')
		# self.TabWidget.insertTab(2,self.w2,'Pre-center input stack (Recommended)')
		self.TabWidget.resize(730,800)
		self.TabWidget.show()


	def viper(self):
			modifiers = QtGui.QApplication.keyboardModifiers()
			if modifiers == QtCore.Qt.ShiftModifier:
					os.system("python -m webbrowser %ssxviper"%SPARX_DOCUMENTATION_WEBSITE)
					return
			#print "Opening a new popup window..."
			#print "Opening a new popup window..."
			#opens the window Poptwodali, and defines its width and height
			#The layout of the Poptwodali window is defined in class Poptwodali(QWidget Window)
			self.w = Popupviper()
			self.w1 = Popupadvparams_viper(self.w.savedparmsdict)
			self.w.w1 = self.w1
			self.TabWidget = QtGui.QTabWidget()
			self.TabWidget.insertTab(0,self.w,'Main')
			self.TabWidget.insertTab(1,self.w1,'Advanced Parameters')
			self.TabWidget.resize(570,740)
			self.TabWidget.show()


	def meridien(self):
			modifiers = QtGui.QApplication.keyboardModifiers()
			if modifiers == QtCore.Qt.ShiftModifier:
					os.system("python -m webbrowser %ssxmeridien"%SPARX_DOCUMENTATION_WEBSITE)
					return
			#print "Opening a new popup window..."
			#opens the window Poptwodali, and defines its width and height
			#The layout of the Poptwodali window is defined in class Poptwodali(QWidget Window)
			self.w = Popupmeridien()
			self.w1 = Popupadvparams_meridien(self.w.savedparmsdict)
			self.w.w1 = self.w1
			self.TabWidget = QtGui.QTabWidget()
			self.TabWidget.insertTab(0,self.w,'Main')
			self.TabWidget.insertTab(1,self.w1,'Advanced Parameters')
			self.TabWidget.resize(570,740)
			self.TabWidget.show()


	def variability3d(self):
			modifiers = QtGui.QApplication.keyboardModifiers()
			if modifiers == QtCore.Qt.ShiftModifier:
					os.system("python -m webbrowser %ssx3dvariability"%SPARX_DOCUMENTATION_WEBSITE)
					return
			#print "Opening a new popup window..."
			#opens the window Poptwodali, and defines its width and height
			#The layout of the Poptwodali window is defined in class Poptwodali(QWidget Window)
			self.w = Popupvariability3d()
			#self.w1 = Popupadvparams_variability3d(self.w.savedparmsdict)
			#self.w.w1 = self.w1
			self.TabWidget = QtGui.QTabWidget()
			self.TabWidget.insertTab(0,self.w,'Main')
			#self.TabWidget.insertTab(1,self.w1,'Advanced Parameters')
			self.TabWidget.resize(570,740)
			self.TabWidget.show()



	def locres(self):
			#print "Opening a new popup window..."
			#opens the window Poptwodali, and defines its width and height
			#The layout of the Poptwodali window is defined in class Poptwodali(QWidget Window)
			self.w = Popuplocres()
			#self.w1 = Popupadvparams_locres(self.w.savedparmsdict)
			#self.w.w1 = self.w1
			self.TabWidget = QtGui.QTabWidget()
			self.TabWidget.insertTab(0,self.w,'Main')
			#self.TabWidget.insertTab(1,self.w1,'Advanced Parameters')
			self.TabWidget.resize(570,740)
			self.TabWidget.show()



	def filterlocal(self):
			#print "Opening a new popup window..."
			#opens the window Poptwodali, and defines its width and height
			#The layout of the Poptwodali window is defined in class Poptwodali(QWidget Window)
			self.w = Popupfilterlocal()
			#self.w1 = Popupadvparams_filterlocal(self.w.savedparmsdict)
			#self.w.w1 = self.w1
			self.TabWidget = QtGui.QTabWidget()
			self.TabWidget.insertTab(0,self.w,'Main')
			#self.TabWidget.insertTab(1,self.w1,'Advanced Parameters')
			self.TabWidget.resize(570,740)
			self.TabWidget.show()

								
	#This is the function info, which is being started when the Pushbutton picbutton of the main window is being clicked
	def info(self):
			#print "Opening a new popup window..."
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
				# self.main.resize(400,450)
				self.main.resize(1000, 755)
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
	from optparse import OptionParser
	progname = os.path.basename(sys.argv[0])
	usage = """%prog [options] <image>
		
Automatic and manual particle selection. This version is specifically aimed at square boxes
for single particle analysis."""
	parser = OptionParser(usage=usage,version=EMANVERSION)

	parser.add_option("--demo", type="string", default="",   help="Name of the demo whose input parameters should be used to initialize the GUI fields: --demo=mpibdb means the input parameters in demo/mpi_bdb will be used, --demo=mpibdbctf means the input parameters in demo/mpi_bdb_ctf will be used")
	global options
	(options, args) = parser.parse_args()
	global DEMO_mpibdbctf
	DEMO_mpibdbctf = 'mpibdbctf'
	global DEMO_mpibdb
	DEMO_mpibdb = 'mpibdb'
	
	global app
	app = App(args)
	app.setWindowIcon(QIcon(get_image_directory()+"sparxicon.png"))
	
	app.main.show()
	app.main.raise_()
	app.exec_()

if __name__ == "__main__":
	main(sys.argv)
