#!/usr/bin/env python

# Author: Muthu Alagappan, m.alagappan901@gmail.com, 07/22/09
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

from EMAN2 import PDBReader, get_image_directory
from libpyGLUtils2 import *
from OpenGL.GL import *
from OpenGL.GLU import *
from emglobjects import EM3DModel, get_default_gl_colors, EMViewportDepthTools, Camera2
from PyQt4 import QtGui, QtCore
from PyQt4.QtCore import Qt
import sys
import weakref
from emimageutil import EMTransformPanel



class AlaRenderer:
	def __init__(self): pass
		
	def __call__(self,res,target):
		try: t1 = res[3].index('CB')
		except: pass

		try: target.makeStick(res, 0, 1)
		except: pass		
		try: target.makeStick(res, 1, 2)
		except: pass
		try: target.makeStick(res, 2, 3)
		except: pass

		try: target.makeStick(res, 1, t1)
		except: pass

class ArgRenderer:
	def __init__(self): pass
		
	def __call__(self,res,target):

		try: t1 = res[3].index('CB')
		except: pass
		try: t2 = res[3].index('CG')
		except: pass
		try: t3 = res[3].index('CD')
		except: pass
		try: t4 = res[3].index('NE')
		except: pass
		try: t5 = res[3].index('CZ')
		except: pass
		try: t6 = res[3].index('NH1')
		except: pass
		try: t7 = res[3].index('NH2')
		except: pass
		


		try: target.makeStick(res, 0, 1)
		except: pass
		try: target.makeStick(res, 1, 2)
		except: pass
		try: target.makeStick(res, 2, 3)
		except: pass

		try: target.makeStick(res, 1, t1)
		except: pass
		try: target.makeStick(res, t1, t2)
		except: pass
		try: target.makeStick(res, t2, t3)
		except: pass
		try: target.makeStick(res, t3, t4)
		except: pass
		try: target.makeStick(res, t4, t5)
		except: pass
		try: target.makeStick(res, t5, t6)
		except: pass
		try: target.makeStick(res, t5, t7)
		except: pass

class AspRenderer:
	def __init__(self): pass
		
	def __call__(self,res,target):
		try: t1 = res[3].index('CB')
		except: pass
		try: t2 = res[3].index('CG')
		except: pass
		try: t3 = res[3].index('OD1')
		except: pass
		try: t4 = res[3].index('OD2')
		except: pass


		try: target.makeStick(res, 0, 1)
		except: pass
		try: target.makeStick(res, 1, 2)
		except: pass
		try: target.makeStick(res, 2, 3)
		except: pass

		try: target.makeStick(res, 1, t1)
		except: pass
		try: target.makeStick(res, t1, t2)
		except: pass
		try: target.makeStick(res, t2, t3)
		except: pass
		try: target.makeStick(res, t2, t4)
		except: pass

class AsnRenderer:
	def __init__(self): pass
		
	def __call__(self,res,target):
		try: t1 = res[3].index('CB')
		except: pass
		try: t2 = res[3].index('CG')
		except: pass
		try: t3 = res[3].index('OD1')
		except: pass
		try: t4 = res[3].index('ND2')
		except: pass


		try: target.makeStick(res, 0, 1)
		except: pass
		try: target.makeStick(res, 1, 2)
		except: pass
		try: target.makeStick(res, 2, 3)
		except: pass

		try: target.makeStick(res, 1, t1)
		except: pass
		try: target.makeStick(res, t1, t2)
		except: pass
		try: target.makeStick(res, t2, t3)
		except: pass
		try: target.makeStick(res, t2, t4)
		except: pass

class CysRenderer:
	def __init__(self): pass
		
	def __call__(self,res,target):
		try: t1 = res[3].index('CB')
		except: pass
		try: t2 = res[3].index('SG')
		except: pass

		try: target.makeStick(res, 0, 1)
		except: pass
		try: target.makeStick(res, 1, 2)
		except: pass
		try: target.makeStick(res, 2, 3)
		except: pass

		try: target.makeStick(res, 1, t1)
		except: pass
		try: target.makeStick(res, t1, t2)
		except: pass

class GlyRenderer:
	def __init__(self): pass
		
	def __call__(self,res,target):
		try: target.makeStick(res, 0, 1)
		except: pass
		try: target.makeStick(res, 1, 2)
		except: pass
		try: target.makeStick(res, 2, 3)
		except: pass

class GlnRenderer:
	def __init__(self): pass
		
	def __call__(self,res,target):
		try: t1 = res[3].index('CB')
		except: pass
		try: t2 = res[3].index('CG')
		except: pass
		try: t3 = res[3].index('CD')
		except: pass
		try: t4 = res[3].index('OE1')
		except: pass
		try: t5 = res[3].index('NE2')
		except: pass


		try: target.makeStick(res, 0, 1)
		except: pass
		try: target.makeStick(res, 1, 2)
		except: pass
		try: target.makeStick(res, 2, 3)
		except: pass

		try: target.makeStick(res, 1, t1)
		except: pass
		try: target.makeStick(res, t1, t2)
		except: pass
		try: target.makeStick(res, t2, t3)
		except: pass
		try: target.makeStick(res, t3, t4)
		except: pass
		try: target.makeStick(res, t3, t5)
		except: pass

class GluRenderer:
	def __init__(self): pass
		
	def __call__(self,res,target):
		try: t1 = res[3].index('CB')
		except: pass
		try: t2 = res[3].index('CG')
		except: pass
		try: t3 = res[3].index('CD')
		except: pass
		try: t4 = res[3].index('OE1')
		except: pass
		try: t5 = res[3].index('OE2')
		except: pass


		try: target.makeStick(res, 0, 1)
		except: pass
		try: target.makeStick(res, 1, 2)
		except: pass
		try: target.makeStick(res, 2, 3)
		except: pass

		try: target.makeStick(res, 1, t1)
		except: pass
		try: target.makeStick(res, t1, t2)
		except: pass
		try: target.makeStick(res, t2, t3)
		except: pass
		try: target.makeStick(res, t3, t4)
		except: pass
		try: target.makeStick(res, t3, t5)
		except: pass

class HisRenderer:
	def __init__(self): pass
		
	def __call__(self,res,target):
		try: t1 = res[3].index('CB')
		except: pass
		try: t2 = res[3].index('CG')
		except: pass
		try: t3 = res[3].index('CD2')
		except: pass
		try: t4 = res[3].index('ND1')
		except: pass
		try: t5 = res[3].index('NE2')
		except: pass
		try: t6 = res[3].index('CE1')
		except: pass


		try: target.makeStick(res, 0, 1)
		except: pass
		try: target.makeStick(res, 1, 2)
		except: pass
		try: target.makeStick(res, 2, 3)
		except: pass

		try: target.makeStick(res, 1, t1)
		except: pass
		try: target.makeStick(res, t1, t2)
		except: pass
		try: target.makeStick(res, t2, t3)
		except: pass
		try: target.makeStick(res, t2, t4)
		except: pass
		try: target.makeStick(res, t3, t5)
		except: pass
		try: target.makeStick(res, t5, t6)
		except: pass
		try: target.makeStick(res, t4, t6)
		except: pass

class IleRenderer:
	def __init__(self): pass
		
	def __call__(self,res,target):
		try: t1 = res[3].index('CB')
		except: pass
		try: t2 = res[3].index('CG1')
		except: pass
		try: t3 = res[3].index('CG2')
		except: pass
		try: t4 = res[3].index('CD1')
		except: pass


		try: target.makeStick(res, 0, 1)
		except: pass
		try: target.makeStick(res, 1, 2)
		except: pass
		try: target.makeStick(res, 2, 3)
		except: pass

		try: target.makeStick(res, 1, t1)
		except: pass
		try: target.makeStick(res, t1, t2)
		except: pass
		try: target.makeStick(res, t1, t3)
		except: pass
		try: target.makeStick(res, t2, t4)
		except: pass

class LeuRenderer:
	def __init__(self): pass
		
	def __call__(self,res,target):
		try: t1 = res[3].index('CB')
		except: pass
		try: t2 = res[3].index('CG')
		except: pass
		try: t3 = res[3].index('CD1')
		except: pass
		try: t4 = res[3].index('CD2')
		except: pass


		try: target.makeStick(res, 0, 1)
		except: pass
		try: target.makeStick(res, 1, 2)
		except: pass
		try: target.makeStick(res, 2, 3)
		except: pass

		try: target.makeStick(res, 1, t1)
		except: pass
		try: target.makeStick(res, t1, t2)
		except: pass
		try: target.makeStick(res, t2, t3)
		except: pass
		try: target.makeStick(res, t2, t4)
		except: pass

class LysRenderer:
	def __init__(self): pass
		
	def __call__(self,res,target):
		try: t1 = res[3].index('CB')
		except: pass
		try: t2 = res[3].index('CG')
		except: pass
		try: t3 = res[3].index('CD')
		except: pass
		try: t4 = res[3].index('CE')
		except: pass
		try: t5 = res[3].index('NZ')
		except: pass


		try: target.makeStick(res, 0, 1)
		except: pass
		try: target.makeStick(res, 1, 2)
		except: pass
		try: target.makeStick(res, 2, 3)
		except: pass

		try: target.makeStick(res, 1, t1)
		except: pass
		try: target.makeStick(res, t1, t2)
		except: pass
		try: target.makeStick(res, t2, t3)
		except: pass
		try: target.makeStick(res, t3, t4)
		except: pass
		try: target.makeStick(res, t4, t5)
		except: pass

class MetRenderer:
	def __init__(self): pass
		
	def __call__(self,res,target):
		try: t1 = res[3].index('CB')
		except: pass
		try: t2 = res[3].index('CG')
		except: pass
		try: t3 = res[3].index('SD')
		except: pass
		try: t4 = res[3].index('CE')
		except: pass


		try: target.makeStick(res, 0, 1)
		except: pass
		try: target.makeStick(res, 1, 2)
		except: pass
		try: target.makeStick(res, 2, 3)
		except: pass

		try: target.makeStick(res, 1, t1)
		except: pass
		try: target.makeStick(res, t1, t2)
		except: pass
		try: target.makeStick(res, t2, t3)
		except: pass
		try: target.makeStick(res, t3, t4)
		except: pass

class PheRenderer:
	def __init__(self): pass
		
	def __call__(self,res,target):
		try: t1 = res[3].index('CB')
		except: pass
		try: t2 = res[3].index('CG')
		except: pass
		try: t3 = res[3].index('CD1')
		except: pass
		try: t4 = res[3].index('CD2')
		except: pass
		try: t5 = res[3].index('CE1')
		except: pass
		try: t6 = res[3].index('CE2')
		except: pass
		try: t7 = res[3].index('CZ')
		except: pass


		try: target.makeStick(res, 0, 1)
		except: pass
		try: target.makeStick(res, 1, 2)
		except: pass
		try: target.makeStick(res, 2, 3)
		except: pass

		try: target.makeStick(res, 1, t1)
		except: pass
		try: target.makeStick(res, t1, t2)
		except: pass
		try: target.makeStick(res, t2, t3)
		except: pass
		try: target.makeStick(res, t2, t4)
		except: pass
		try: target.makeStick(res, t3, t5)
		except: pass
		try: target.makeStick(res, t4, t6)
		except: pass
		try: target.makeStick(res, t5, t7)
		except: pass
		try: target.makeStick(res, t6, t7)
		except: pass

class ProRenderer:
	def __init__(self): pass
		
	def __call__(self,res,target):
		try: t1 = res[3].index('CB')
		except: pass
		try: t2 = res[3].index('CG')
		except: pass
		try: t3 = res[3].index('CD')
		except: pass
		try: t4 = res[3].index('N')
		except: pass


		try: target.makeStick(res, 0, 1)
		except: pass
		try: target.makeStick(res, 1, 2)
		except: pass
		try: target.makeStick(res, 2, 3)
		except: pass

		try: target.makeStick(res, 1, t1)
		except: pass
		try: target.makeStick(res, t1, t2)
		except: pass
		try: target.makeStick(res, t2, t3)
		except: pass
		try: target.makeStick(res, t3, t4)
		except: pass

class SerRenderer:
	def __init__(self): pass
		
	def __call__(self,res,target):
		try: t1 = res[3].index('CB')
		except: pass
		try: t2 = res[3].index('OG')
		except: pass

		try: target.makeStick(res, 0, 1)
		except: pass
		try: target.makeStick(res, 1, 2)
		except: pass
		try: target.makeStick(res, 2, 3)
		except: pass

		try: target.makeStick(res, 1, t1)
		except: pass
		try: target.makeStick(res, t1, t2)
		except: pass

class ThrRenderer:
	def __init__(self): pass
		
	def __call__(self,res,target):
		try: t1 = res[3].index('CB')
		except: pass
		try: t2 = res[3].index('CG2')
		except: pass
		try: t3 = res[3].index('OG1')
		except: pass

		try: target.makeStick(res, 0, 1)
		except: pass
		try: target.makeStick(res, 1, 2)
		except: pass
		try: target.makeStick(res, 2, 3)
		except: pass

		try: target.makeStick(res, 1, t1)
		except: pass
		try: target.makeStick(res, t1, t2)
		except: pass
		try: target.makeStick(res, t1, t3)
		except: pass

class TrpRenderer:
	def __init__(self): pass
		
	def __call__(self,res,target):
		try: t1 = res[3].index('CB')
		except: pass
		try: t2 = res[3].index('CG')
		except: pass
		try: t3 = res[3].index('CD1')
		except: pass
		try: t4 = res[3].index('CD2')
		except: pass
		try: t5 = res[3].index('NE1')
		except: pass
		try: t6 = res[3].index('CE2')
		except: pass
		try: t7 = res[3].index('CE3')
		except: pass
		try: t8 = res[3].index('CZ3')
		except: pass
		try: t9 = res[3].index('CH2')
		except: pass
		try: t10 = res[3].index('CZ2')
		except: pass

		try: target.makeStick(res, 0, 1)
		except: pass
		try: target.makeStick(res, 1, 2)
		except: pass
		try: target.makeStick(res, 2, 3)
		except: pass

		try: target.makeStick(res, 1, t1)
		except: pass
		try: target.makeStick(res, t1, t2)
		except: pass
		try: target.makeStick(res, t2, t3)
		except: pass
		try: target.makeStick(res, t2, t4)
		except: pass
		try: target.makeStick(res, t3, t5)
		except: pass
		try: target.makeStick(res, t5, t6)
		except: pass
		try: target.makeStick(res, t4, t6)
		except: pass
		try: target.makeStick(res, t4, t7)
		except: pass
		try: target.makeStick(res, t7, t8)
		except: pass
		try: target.makeStick(res, t8, t9)
		except: pass
		try: target.makeStick(res, t10, t9)
		except: pass

class TyrRenderer:
	def __init__(self): pass
		
	def __call__(self,res,target):
		try: t1 = res[3].index('CB')
		except: pass
		try: t2 = res[3].index('CG')
		except: pass
		try: t3 = res[3].index('CD1')
		except: pass
		try: t4 = res[3].index('CD2')
		except: pass
		try: t5 = res[3].index('CE1')
		except: pass
		try: t6 = res[3].index('CE2')
		except: pass
		try: t7 = res[3].index('CZ')
		except: pass
		try: t8 = res[3].index('OH')
		except: pass


		try: target.makeStick(res, 0, 1)
		except: pass
		try: target.makeStick(res, 1, 2)
		except: pass
		try: target.makeStick(res, 2, 3)
		except: pass

		try: target.makeStick(res, 1, t1)
		except: pass
		try: target.makeStick(res, t1, t2)
		except: pass
		try: target.makeStick(res, t2, t3)
		except: pass
		try: target.makeStick(res, t2, t4)
		except: pass
		try: target.makeStick(res, t3, t5)
		except: pass
		try: target.makeStick(res, t4, t6)
		except: pass
		try: target.makeStick(res, t5, t7)
		except: pass
		try: target.makeStick(res, t6, t7)
		except: pass
		try: target.makeStick(res, t7, t8)
		except: pass

class ValRenderer:
	def __init__(self): pass
		
	def __call__(self,res,target):
		try: t1 = res[3].index('CB')
		except: pass
		try: t2 = res[3].index('CG2')
		except: pass
		try: t3 = res[3].index('CG1')
		except: pass

		try: target.makeStick(res, 0, 1)
		except: pass
		try: target.makeStick(res, 1, 2)
		except: pass
		try: target.makeStick(res, 2, 3)
		except: pass

		try: target.makeStick(res, 1, t1)
		except: pass
		try: target.makeStick(res, t1, t2)
		except: pass
		try: target.makeStick(res, t1, t3)
		except: pass



class EMPDBModel(EM3DModel):
	def __init__(self, gl_widget):
		self.fName = ""
		self.text = self.fName
		self.dl = None
		EM3DModel.__init__(self, gl_widget)
		# basic shapes will be stored in these lists
		self.gq = None # will be a glu quadric
		self.cylinderdl = 0 # will be a cylinder with no caps
		self.diskdl = 0 # this will be a flat disk
		self.spheredl = 0 # this will be a low resolution sphere
		self.highresspheredl = 0 # high resolution sphere
		self.cappedcylinderdl = 0 # a capped cylinder
		self.first_render_flag = True # this is used to catch the first call to the render function - so you can do an GL context sensitive initialization when you know there is a valid context
	
		self.inspector = None # will be the inspector
		self.radius = 100
		self.perspective = False
		self.colors = get_default_gl_colors()
		self.vdtools = EMViewportDepthTools(self)
		self.cam = Camera2(self)
		self.cam.basicmapping = True #new by Ross
		
		self.side_chain_renderer = {}
		self.side_chain_renderer["ALA"] = AlaRenderer()
		self.side_chain_renderer["ARG"] = ArgRenderer()
		self.side_chain_renderer["ASP"] = AspRenderer()
		self.side_chain_renderer["ASN"] = AsnRenderer()
		self.side_chain_renderer["CYS"] = CysRenderer()
		self.side_chain_renderer["GLY"] = GlyRenderer()
		self.side_chain_renderer["GLN"] = GlnRenderer()
		self.side_chain_renderer["GLU"] = GluRenderer()
		self.side_chain_renderer["HIS"] = HisRenderer()
		self.side_chain_renderer["ILE"] = IleRenderer()
		self.side_chain_renderer["LEU"] = LeuRenderer()
		self.side_chain_renderer["LYS"] = LysRenderer()
		self.side_chain_renderer["MET"] = MetRenderer()
		self.side_chain_renderer["PHE"] = PheRenderer()
		self.side_chain_renderer["PRO"] = ProRenderer()
		self.side_chain_renderer["SER"] = SerRenderer()
		self.side_chain_renderer["THR"] = ThrRenderer()
		self.side_chain_renderer["TRP"] = TrpRenderer()
		self.side_chain_renderer["TYR"] = TyrRenderer()
		self.side_chain_renderer["VAL"] = ValRenderer()
	def buildResList (self): # calls PDBReader to read the given pdb file and create a list (self.allResidues) of lists (x,y,z,atom name, residue name) of lists (all the values for that residue)

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
		currentRes = point_resNum[0]


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
	def createDefault(self):
		return          #display a default pdb here, currently not done
	def current_text(self): return self.text
	def cylinder_to_from(self,next,prev,scale=0.5):
		dx = next[0] - prev[0]
		dy = next[1] - prev[1]
		dz = next[2] - prev[2]
		from math import sqrt,acos,atan2,pi
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
	def draw_objects(self):
		self.init_basic_shapes() # only does something the first time you call it
		if (self.text == ""): 
			#default drawing
			self.createDefault()
			return
			
		
		#self.get_gl_widget().makeCurrent()
		if (self.text != self.fName): 
			if (self.dl != None): glDeleteLists(self.dl, 1)
			self.dl=None
			self.fName = self.text

		if (self.dl == None): #self.dl is the display list, every time a new file is added, this is changed back to None
			self.dl=glGenLists(1)
			glNewList(self.dl,GL_COMPILE)
			self.buildResList()

			for res in self.allResidues: #goes through self.allResidues and displays a sphere for every atom in the pdb
				for i in range (0, len(res[0])):
					glPushMatrix()
					glTranslate(res[0][i], res[1][i], res[2][i])
					glScale(1,1,1)
					if (str(res[3][i])[0] == 'C'): self.load_gl_color("white")
					elif (str(res[3][i])[0] == 'N'): self.load_gl_color("green")
					elif (str(res[3][i])[0] == 'O'): self.load_gl_color("blue")
					elif (str(res[3][i])[0] == 'S'): self.load_gl_color("red")
					else: self.load_gl_color("silver")
					glCallList(self.highresspheredl)
					glPopMatrix()
			
#			self.load_gl_color("silver")
			for k in range (0, len(self.allResidues)):
				
				res = self.allResidues[k]
				key =  res[4][0]
				if self.side_chain_renderer.has_key(key): #goes through each residue and draws the newtwork of sticks connecting atoms
					self.side_chain_renderer[key](res,self)
					continue


				if (k!=0): #connects residues together from the nitrogen of one residue to the O of the next residue
				
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

		try: glCallList(self.dl)
		except: 
			print "call list failed",self.dl
			glDeleteLists(self.dl,1)
			self.dl = None
	def init_basic_shapes(self):
		#self.get_gl_widget().makeCurrent()
		if self.gq == None:
			
			self.gq=gluNewQuadric() # a quadric for general use
			gluQuadricDrawStyle(self.gq,GLU_FILL)
			gluQuadricNormals(self.gq,GLU_SMOOTH)
			gluQuadricOrientation(self.gq,GLU_OUTSIDE)
			gluQuadricTexture(self.gq,GL_FALSE)
		
		if ( self.cylinderdl == 0 ):
			self.cylinderdl=glGenLists(1)
				
			glNewList(self.cylinderdl,GL_COMPILE)
			glPushMatrix()
			gluCylinder(self.gq,1.0,1.0,1.0,12,2)
			glPopMatrix()
				
			glEndList()
		
		if self.diskdl == 0:
			self.diskdl=glGenLists(1)
				
			glNewList(self.diskdl,GL_COMPILE)
			gluDisk(self.gq,0,1,12,2)
			glEndList()
		
		if self.spheredl == 0:
			self.spheredl=glGenLists(1)
				
			glNewList(self.spheredl,GL_COMPILE)
			gluSphere(self.gq,.5,4,2)
			glEndList()

		
		if self.highresspheredl == 0:
			self.highresspheredl=glGenLists(1)
				
			glNewList(self.highresspheredl,GL_COMPILE)
			gluSphere(self.gq,.5,16,16)
			glEndList()
			
		if ( self.cappedcylinderdl == 0 ):
			self.cappedcylinderdl=glGenLists(1)
			glNewList(self.cappedcylinderdl,GL_COMPILE)
			glCallList(self.cylinderdl)
			glPushMatrix()
			glTranslate(0,0,1)
			glCallList(self.diskdl)
			glPopMatrix()
			glPushMatrix()
			glRotate(180,0,1,0)
			glCallList(self.diskdl)
			glPopMatrix()
			glEndList()
	def get_inspector(self):
		if self.inspector == None:
			self.inspector = EMPDBInspector(self)
		return self.inspector
	def get_pdb_file(self):
		return self.fName
	def get_type(self):
		return "EMPDBModel"
	def load_gl_color(self,name):
		color = self.colors[name]
		glColor(color["ambient"])
		glMaterial(GL_FRONT,GL_AMBIENT,color["ambient"])
		glMaterial(GL_FRONT,GL_DIFFUSE,color["diffuse"])
		glMaterial(GL_FRONT,GL_SPECULAR,color["specular"])
		glMaterial(GL_FRONT,GL_EMISSION,color["emission"])
		glMaterial(GL_FRONT,GL_SHININESS,color["shininess"])
	def makeStick (self, res, index1, index2): #draws a cylinder between two atoms once the index for start and stop is given
		n = [0,0,0]
		p = [0,0,0]
		p[0] = res[0][index1]
		p[1] = res[1][index1]
		p[2] = res[2][index1]

		n[0] = res[0][index2]
		n[1] = res[1][index2]
		n[2] = res[2][index2]
		self.cylinder_to_from(n, p, 0.2)
	def render(self):
		if self.first_render_flag:
			if not self.perspective:self.get_gl_widget().load_orthographic()
			else: self.get_gl_widget().load_perspective()
			self.first_render_flag = False
					
		#self.vdtools.set_update_P_inv()
		glPushMatrix()
		self.cam.position(True)
		# the ones are dummy variables atm... they don't do anything
		self.vdtools.update(1,1)
		glPopMatrix()
		
		glPushMatrix()
		self.cam.position() #FIXME: figure out why translation doesn't work
		self.draw_objects()
		glPopMatrix()
	def set_current_text(self,text):  #changes self.text and updatesGL, when self.text changes, it redisplays using the new file
		self.text = text
		self.get_inspector().text.setText(self.text)
		self.updateGL()
		
class EMPDBInspector(QtGui.QWidget):
	def __init__(self,target,enable_advanced=False):
		QtGui.QWidget.__init__(self)
		self.target = weakref.ref(target)

		self.rotation_sliders = EMTransformPanel(target,self)
		
		self.text = QtGui.QLineEdit()
		text_value = self.target().current_text()
		if text_value:
			self.text.setText(text_value)
		self.browse = QtGui.QPushButton("Browse")

		hbl1 = QtGui.QHBoxLayout()
		hbl1.addWidget(self.text)
		hbl1.addWidget(self.browse)

		vbl = QtGui.QVBoxLayout()
		vbl.setMargin(0)
		vbl.setSpacing(6)
		vbl.addLayout(hbl1)
		
		self.rotation_sliders.addWidgets(vbl)

		self.setLayout(vbl)
		
		QtCore.QObject.connect(self.text, QtCore.SIGNAL("textEdited(const QString&)"), self.on_text_change)
		QtCore.QObject.connect(self.browse, QtCore.SIGNAL("clicked(bool)"), self.on_browse)
	
	def on_text_change(self,text):
		print "Use the Browse button to update the pdb file"

	def on_browse(self):
		import os
		self.fileName = QtGui.QFileDialog.getOpenFileName(self, "open file", os.getcwd(), "Text files (*.pdb)")
		if (self.fileName == ""): return
		self.target().set_current_text(str(self.fileName)) #self.target().text and self.text are what the user sees. 
		self.text.setText(self.fileName) #if self.text changes, then self.fName becomes self.text and the image regenerates	
		self.target().updateGL()
	
	def update_rotations(self,t3d):
		self.rotation_sliders.update_rotations(t3d)
	
	def set_scale(self,val):
		self.rotation_sliders.set_scale(val)
	
	def set_xy_trans(self, x, y):
		self.rotation_sliders.set_xy_trans(x,y)
	
	def set_xyz_trans(self,x,y,z):
		self.rotation_sliders.set_xyz_trans(x,y,z)

if __name__ == '__main__':
	from emapplication import EMApp
	from emimage3d import EMImage3DWidget
	em_app = EMApp()

	window = EMImage3DWidget()
	pdb_model = EMPDBModel(window)
	window.add_model(pdb_model)

	em_app.show()
	em_app.execute()

