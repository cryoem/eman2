#!/usr/bin/env python

# Author: Muthu Alagappan, 07/22/09
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

from em3dmodule import *
from EMAN2 import PDBReader

class EMPDBViewer(EM3DModule):
	def __init__(self, application=None):
		EM3DModule.__init__(self,application)
		#self.fName = raw_input ("Enter the file name of a pdb file: ")
		self.fName = ""
		self.text = self.fName
		self.dl = None
	
	def current_text(self): return self.text
	
	def set_current_text(self,text):
		self.text = text
	
	def get_inspector(self):
		if self.inspector == None:
			self.inspector = EMPDBInspector(self)
		return self.inspector
		
	def draw_objects(self):

		if (self.text == ""): return
		
		if (self.text != self.fName): 
			self.dl=None
			self.fName = self.text

		if (self.dl == None):
			self.dl=glGenLists(1)
			glNewList(self.dl,GL_COMPILE)
			self.buildResList()
			
			for res in self.allResidues:
				for i in range (0, len(res[0])):
					glPushMatrix()
					glTranslate(res[0][i], res[1][i], res[2][i])
					glScale(1,1,1)
					if (str(res[3][i])[0] == 'C'): self.load_gl_color("white")
					elif (str(res[3][i])[0] == 'N'): self.load_gl_color("green")
					elif (str(res[3][i])[0] == 'O'): self.load_gl_color("blue")
					elif (str(res[3][i])[0] == 'S'): self.load_gl_color("red")
					else: self.load_gl_color("silver")
					glCallList(self.spheredl)
					glPopMatrix()
		
			for k in range (0, len(self.allResidues)):
				res = self.allResidues[k]
				if res[4][0] == "ALA":
					#t1 = res[3].index('N')
					#t2 = res[3].index('CA')
					#t3 = res[3].index('C')
					#t4 = res[3].index('O')
					t1 = res[3].index('CB')

					self.makeStick(res, 0, 1)
					self.makeStick(res, 1, 2)
					self.makeStick(res, 2, 3)

					self.makeStick(res, 1, t1)
				
			
				elif res[4][0] == "ARG":
					t1 = res[3].index('CB')
					t2 = res[3].index('CG')
					t3 = res[3].index('CD')
					t4 = res[3].index('NE')
					t5 = res[3].index('CZ')
					t6 = res[3].index('NH1')
					t7 = res[3].index('NH2')


					self.makeStick(res, 0, 1)
					self.makeStick(res, 1, 2)
					self.makeStick(res, 2, 3)

					self.makeStick(res, 1, t1)
					self.makeStick(res, t1, t2)
					self.makeStick(res, t2, t3)
					self.makeStick(res, t3, t4)
					self.makeStick(res, t4, t5)
					self.makeStick(res, t5, t6)
					self.makeStick(res, t5, t7)

			
				elif res[4][0] == "ASP":
					t1 = res[3].index('CB')
					t2 = res[3].index('CG')
					t3 = res[3].index('OD1')
					t4 = res[3].index('OD2')


					self.makeStick(res, 0, 1)
					self.makeStick(res, 1, 2)
					self.makeStick(res, 2, 3)

					self.makeStick(res, 1, t1)
					self.makeStick(res, t1, t2)
					self.makeStick(res, t2, t3)
					self.makeStick(res, t2, t4)
				
			
				elif res[4][0] == "ASN":
					t1 = res[3].index('CB')
					t2 = res[3].index('CG')
					t3 = res[3].index('OD1')
					t4 = res[3].index('ND2')


					self.makeStick(res, 0, 1)
					self.makeStick(res, 1, 2)
					self.makeStick(res, 2, 3)

					self.makeStick(res, 1, t1)
					self.makeStick(res, t1, t2)
					self.makeStick(res, t2, t3)
					self.makeStick(res, t2, t4)

				elif res[4][0] == "CYS":
					t1 = res[3].index('CB')
					t2 = res[3].index('SG')

					self.makeStick(res, 0, 1)
					self.makeStick(res, 1, 2)
					self.makeStick(res, 2, 3)

					self.makeStick(res, 1, t1)
					self.makeStick(res, t1, t2)


				elif res[4][0] == "GLY":
					self.makeStick(res, 0, 1)
					self.makeStick(res, 1, 2)
					self.makeStick(res, 2, 3)

				elif res[4][0] == "GLN":
					t1 = res[3].index('CB')
					t2 = res[3].index('CG')
					t3 = res[3].index('CD')
					t4 = res[3].index('OE1')
					t5 = res[3].index('NE2')


					self.makeStick(res, 0, 1)
					self.makeStick(res, 1, 2)
					self.makeStick(res, 2, 3)

					self.makeStick(res, 1, t1)
					self.makeStick(res, t1, t2)
					self.makeStick(res, t2, t3)
					self.makeStick(res, t3, t4)
					self.makeStick(res, t3, t5)

				elif res[4][0] == "GLU":
					t1 = res[3].index('CB')
					t2 = res[3].index('CG')
					t3 = res[3].index('CD')
					t4 = res[3].index('OE1')
					t5 = res[3].index('OE2')


					self.makeStick(res, 0, 1)
					self.makeStick(res, 1, 2)
					self.makeStick(res, 2, 3)

					self.makeStick(res, 1, t1)
					self.makeStick(res, t1, t2)
					self.makeStick(res, t2, t3)
					self.makeStick(res, t3, t4)
					self.makeStick(res, t3, t5)

			
				elif res[4][0] == "HIS":
					t1 = res[3].index('CB')
					t2 = res[3].index('CG')
					t3 = res[3].index('CD2')
					t4 = res[3].index('ND1')
					t5 = res[3].index('NE2')
					t6 = res[3].index('CE1')


					self.makeStick(res, 0, 1)
					self.makeStick(res, 1, 2)
					self.makeStick(res, 2, 3)

					self.makeStick(res, 1, t1)
					self.makeStick(res, t1, t2)
					self.makeStick(res, t2, t3)
					self.makeStick(res, t2, t4)
					self.makeStick(res, t3, t5)
					self.makeStick(res, t5, t6)
					self.makeStick(res, t4, t6)


				elif res[4][0] == "ILE":	
					t1 = res[3].index('CB')
					t2 = res[3].index('CG1')
					t3 = res[3].index('CG2')
					t4 = res[3].index('CD1')


					self.makeStick(res, 0, 1)
					self.makeStick(res, 1, 2)
					self.makeStick(res, 2, 3)

					self.makeStick(res, 1, t1)
					self.makeStick(res, t1, t2)
					self.makeStick(res, t1, t3)
					self.makeStick(res, t2, t4)

			
				elif res[4][0] == "LEU":
					t1 = res[3].index('CB')
					t2 = res[3].index('CG')
					t3 = res[3].index('CD1')
					t4 = res[3].index('CD2')


					self.makeStick(res, 0, 1)
					self.makeStick(res, 1, 2)
					self.makeStick(res, 2, 3)

					self.makeStick(res, 1, t1)
					self.makeStick(res, t1, t2)
					self.makeStick(res, t2, t3)
					self.makeStick(res, t2, t4)

				elif res[4][0] == "LYS":
					t1 = res[3].index('CB')
					t2 = res[3].index('CG')
					t3 = res[3].index('CD')
					t4 = res[3].index('CE')
					t5 = res[3].index('NZ')


					self.makeStick(res, 0, 1)
					self.makeStick(res, 1, 2)
					self.makeStick(res, 2, 3)

					self.makeStick(res, 1, t1)
					self.makeStick(res, t1, t2)
					self.makeStick(res, t2, t3)
					self.makeStick(res, t3, t4)
					self.makeStick(res, t4, t5)

			
				elif res[4][0] == "MET":
					t1 = res[3].index('CB')
					t2 = res[3].index('CG')
					t3 = res[3].index('SD')
					t4 = res[3].index('CE')


					self.makeStick(res, 0, 1)
					self.makeStick(res, 1, 2)
					self.makeStick(res, 2, 3)

					self.makeStick(res, 1, t1)
					self.makeStick(res, t1, t2)
					self.makeStick(res, t2, t3)
					self.makeStick(res, t3, t4)

				elif res[4][0] == "PHE":
					t1 = res[3].index('CB')
					t2 = res[3].index('CG')
					t3 = res[3].index('CD1')
					t4 = res[3].index('CD2')
					t5 = res[3].index('CE1')
					t6 = res[3].index('CE2')
					t7 = res[3].index('CZ')


					self.makeStick(res, 0, 1)
					self.makeStick(res, 1, 2)
					self.makeStick(res, 2, 3)

					self.makeStick(res, 1, t1)
					self.makeStick(res, t1, t2)
					self.makeStick(res, t2, t3)
					self.makeStick(res, t2, t4)
					self.makeStick(res, t3, t5)
					self.makeStick(res, t4, t6)
					self.makeStick(res, t5, t7)
					self.makeStick(res, t6, t7)

				elif res[4][0] == "PRO":
					t1 = res[3].index('CB')
					t2 = res[3].index('CG')
					t3 = res[3].index('CD')
					t4 = res[3].index('N')


					self.makeStick(res, 0, 1)
					self.makeStick(res, 1, 2)
					self.makeStick(res, 2, 3)

					self.makeStick(res, 1, t1)
					self.makeStick(res, t1, t2)
					self.makeStick(res, t2, t3)
					self.makeStick(res, t3, t4)

				elif res[4][0] == "SER":
					t1 = res[3].index('CB')
					t2 = res[3].index('OG')

					self.makeStick(res, 0, 1)
					self.makeStick(res, 1, 2)
					self.makeStick(res, 2, 3)

					self.makeStick(res, 1, t1)
					self.makeStick(res, t1, t2)

			
				elif res[4][0] == "THR":
					t1 = res[3].index('CB')
					t2 = res[3].index('CG2')
					t3 = res[3].index('OG1')

					self.makeStick(res, 0, 1)
					self.makeStick(res, 1, 2)
					self.makeStick(res, 2, 3)

					self.makeStick(res, 1, t1)
					self.makeStick(res, t1, t2)
					self.makeStick(res, t1, t3)

				elif res[4][0] == "TRP":
					t1 = res[3].index('CB')
					t2 = res[3].index('CG')
					t3 = res[3].index('CD1')
					t4 = res[3].index('CD2')
					t5 = res[3].index('NE1')
					t6 = res[3].index('CE2')
					t7 = res[3].index('CE3')
					t8 = res[3].index('CZ3')
					t9 = res[3].index('CH2')
					t10 = res[3].index('CZ2')

					self.makeStick(res, 0, 1)
					self.makeStick(res, 1, 2)
					self.makeStick(res, 2, 3)

					self.makeStick(res, 1, t1)
					self.makeStick(res, t1, t2)
					self.makeStick(res, t2, t3)
					self.makeStick(res, t2, t4)
					self.makeStick(res, t3, t5)
					self.makeStick(res, t5, t6)
					self.makeStick(res, t4, t6)
					self.makeStick(res, t4, t7)
					self.makeStick(res, t7, t8)
					self.makeStick(res, t8, t9)
					self.makeStick(res, t10, t9)

				elif res[4][0] == "TYR":
					t1 = res[3].index('CB')
					t2 = res[3].index('CG')
					t3 = res[3].index('CD1')
					t4 = res[3].index('CD2')
					t5 = res[3].index('CE1')
					t6 = res[3].index('CE2')
					t7 = res[3].index('CZ')
					t8 = res[3].index('OH')


					self.makeStick(res, 0, 1)
					self.makeStick(res, 1, 2)
					self.makeStick(res, 2, 3)

					self.makeStick(res, 1, t1)
					self.makeStick(res, t1, t2)
					self.makeStick(res, t2, t3)
					self.makeStick(res, t2, t4)
					self.makeStick(res, t3, t5)
					self.makeStick(res, t4, t6)
					self.makeStick(res, t5, t7)
					self.makeStick(res, t6, t7)
					self.makeStick(res, t7, t8)

				elif res[4][0] == "VAL":
					t1 = res[3].index('CB')
					t2 = res[3].index('CG2')
					t3 = res[3].index('CG1')

					self.makeStick(res, 0, 1)
					self.makeStick(res, 1, 2)
					self.makeStick(res, 2, 3)

					self.makeStick(res, 1, t1)
					self.makeStick(res, t1, t2)
					self.makeStick(res, t1, t3)

				if (k!=0):
				
					nt = [0,0,0]
					pt = [0,0,0]
					nt[0] = res[0][0]
					nt[1] = res[1][0]
					nt[2] = res[2][0]

					pt[0] = self.allResidues[(k-1)][0][2]
					pt[1] = self.allResidues[(k-1)][1][2]
					pt[2] = self.allResidues[(k-1)][2][2]
					self.cylinder_to_from(nt, pt, 0.2)
			glEndList()

		glCallList(self.dl)


	def makeStick (self, res, index1, index2):
		n = [0,0,0]
		p = [0,0,0]
		p[0] = res[0][index1]
		p[1] = res[1][index1]
		p[2] = res[2][index1]

		n[0] = res[0][index2]
		n[1] = res[1][index2]
		n[2] = res[2][index2]
		self.cylinder_to_from(n, p, 0.2)	

	def buildResList (self):

		self.allResidues = []
		
		try:
			f = open(self.fName)
			f.close()
		except IOError:	
			print "Sorry, the file name \"" + str(self.fName) + "\" does not exist"
			sys.exit()
		
   		self.a = PDBReader()
    		self.a.read_from_pdb(self.fName)
    		point_x = self.a.get_x()
   		point_y = self.a.get_y()
	        point_z = self.a.get_z()
		point_atomName = self.a.get_atomName()
		point_resName = self.a.get_resName()
		point_resNum = self.a.get_resNum()
		x =[]
		y =[]
		z =[]
		atomName =[]
		resName = []
		amino = []
		currentRes = 1

    		for i in range(0, len(point_x)):
        		if (point_resNum[i]==currentRes):
           			x.append(point_x[i])
            			y.append(point_y[i])
            			z.append(point_z[i])
				temp = point_atomName[i]
				temp2 = temp.strip()
				atomName.append(temp2)
            			resName.append(point_resName[i])
       			else:
            			currentRes = point_resNum[i]
				amino.append(x[:])
				amino.append(y[:])
				amino.append(z[:])
				amino.append(atomName[:])
				amino.append(resName[:])
				self.allResidues.append(amino[:])
				del amino[:]
            			del x[:]
            			del y[:]
            			del z[:]
            			del atomName[:]
            			del resName[:]
           			x.append(point_x[i])
            			y.append(point_y[i])
            			z.append(point_z[i])
				temp = point_atomName[i]
				temp2 = temp.strip()
				atomName.append(temp2)
            			resName.append(point_resName[i])
			if (i == (len(point_x)-1)): 
				amino.append(x[:])
				amino.append(y[:])
				amino.append(z[:])
				amino.append(atomName[:])
				amino.append(resName[:])
				self.allResidues.append(amino[:])
				break


	def cylinder_to_from(self,next,prev,scale=0.5):

		dx = next[0] - prev[0]
		dy = next[1] - prev[1]
		dz = next[2] - prev[2]
		try:
			length = sqrt(dx**2 + dy**2 + dz**2)
		except: return
		if length == 0: return
		

		alt = acos(dz/length)*180.0/pi
		phi = atan2(dy,dx)*180.0/pi
		
		glPushMatrix()
		glTranslatef(prev[0],prev[1],prev[2] )
		glRotatef(90+phi,0,0,1)
		glRotatef(alt,1,0,0)
		glScalef(scale,scale,length)
		self.load_gl_color("silver")
		glCallList(self.cylinderdl)

		glPopMatrix()		
		
class EMPDBInspector(EM3DInspector):
	def __init__(self,target,enable_advanced=False):
		EM3DInspector.__init__(self,target,enable_advanced)
		self.tabwidget.insertTab(0,self.get_example_tab(),"Main")
		self.tabwidget.setCurrentIndex(0)
			
	def get_example_tab(self):
		'''
		@return an QtGui.QWidget, i.e. for insertion into a tab widget, or layour, etc
		'''
		widget = QtGui.QWidget()
		vbl = QtGui.QVBoxLayout(widget)
		vbl.setMargin(0)
		vbl.setSpacing(6)

		hbl1 = QtGui.QHBoxLayout()
		self.text = QtGui.QLineEdit()
		self.text.setText(self.target().current_text())
		#text_label = QtGui.QLabel("Enter Text:",self)
		#hbl1.addWidget(text_label)
		hbl1.addWidget(self.text)
		self.browse = QtGui.QPushButton("Browse")
		hbl1.addWidget(self.browse)
		vbl.addLayout(hbl1)
		
		QtCore.QObject.connect(self.text, QtCore.SIGNAL("textChanged(const QString&)"), self.on_text_change)
		QtCore.QObject.connect(self.browse, QtCore.SIGNAL("clicked(bool)"), self.on_browse)
		
		return widget
	
	def on_text_change(self,text):
		self.target().set_current_text(str(text))
		self.target().updateGL()

	def on_browse(self):
		self.fileName = QtGui.QFileDialog.getOpenFileName(self, "open file", "/home", "Text files (*.pdb)")
		if (self.fileName == ""): return
		self.target().set_current_text(str(self.fileName)) #self.target().text and self.text are what the user sees. 
		self.text.setText(self.fileName) #if self.text changes, then self.fName becomes self.text and the image regenerates
		self.target().updateGL()
	
if __name__ == '__main__':
	from emapplication import EMStandAloneApplication
	em_app = EMStandAloneApplication()
	window = EMPDBViewer()
	em_app.show()
	em_app.execute()

