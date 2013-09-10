#!/usr/bin/env python
#
# Author: Steven Ludtke, 11/01/2007 (sludtke@bcm.edu)
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

from OpenGL import GL,GLUT
from math import *
from EMAN2 import get_3d_font_renderer, Util
import warnings

from libpyGLUtils2 import *


def initCircle():
	"""Call this static function once to initialize necessary display lists"""
	if EMShape.dlists>=0 and GL.glIsList(EMShape.dlists): return
	EMShape.dlists=GL.glGenLists(1)
	GL.glNewList(EMShape.dlists,GL.GL_COMPILE)
	GL.glBegin(GL.GL_LINE_LOOP)
	d2r=pi/180.0
	for i in range(90):
		GL.glVertex(sin(i*d2r*4.0),cos(i*d2r*4.0))
	GL.glEnd()
	GL.glEndList()

def shidentity(x,y) : return x,y

class EMShape:
	
	"""This class represents a geometric shape which can be used to annotate
	the various data display widgets in EMAN2. The 'scr' shapes are in screen
	coordinates, and the others are in data coordinates. Shapes are initialized
	and read out in the form of short lists or tuples. Note that no validation is
	performed. A shape that is not understood will simply not be rendered, meaning
	the programmer may create 'invisible' shapes for out-of-band use. Colors
	are on the range 0-1.0 
	
		0               1  2  3  4   5     6     7     8        9  10  11  12  13  14  15  16  17  19  20
		"rect"          R  G  B  x0  y0    x1    y1    linew
		"rectpoint"     R  G  B  x0  y0    x1    y1    linew
		"rectline"      R  G  B  x0  y0    x1    y1    boxw     linew
		"rcircle"       R  G  B  x0  y0    x1    y1    linew
		"rcirclepoint"  R  G  B  x0  y0    x1    y1    linew
		"line"          R  G  B  x0  y0    x1    y1    linew
		"label"         R  G  B  x0  y0    text  size  linew
		"circle"        R  G  B  x0  y0    r     linew
		"ellipse"       R  G  B  x0  y0    r1    r2    ang1     linew
		"scrrect"       R  G  B  x0  y0    x1    y1    linew
		"scrline"       R  G  B  x0  y0    x1    y1    linew
		"scrlabel"      R  G  B  x0  y0    text  size  linew
		"scrcircle"     R  G  B  x0  y0    r     linew
		"point"         R  G  B  x0  y0    r
		"mask"		R  G  B  xi0 yi0   xo0   yo0   xi1      yi1 xo1 yo1 xi2 yi2 xo2 yo2 xi3 yi3 xo3 yo3
		"linemask"	R  G  B  x0  y0    x1    y1    x2       y2  x3  y3  x4  y4  linew
		"hidden"		anything, not rendered
	"""
	font_renderer = None
#	glutinit = True
	dlists=-1

	def __init__(self,init=None) :
		"""init is a list/tuple containing the above parameters describing the shape"""
#		from emglobjects import init_glut
		#if EMShape.glutinit:
			#init_glut()
			#EMShape.glutint = False
		if init : self.shape=list(init)
		else : self.shape=["None",0,0,0,0,0,0,0,0]
		
		initCircle()
		self.blend = 1.0
		self.blendinc = 0.2
		self.isanimated = False
		
		#These are now allocated at the window level, not the shape level !
		#if EMShape.font_renderer == None:
			
			#try:
				#EMShape.font_renderer = get_3d_font_renderer()
				#EMShape.font_renderer.set_face_size(16)
				#EMShape.font_renderer.set_font_mode(FTGLFontMode.TEXTURE)
			#except:
				#EMShape.font_renderer = 1
	
	def __getitem__(self,key): return self.shape[key]
		
	def __setitem__(self,key,value):
		self.shape[key]=value
	
	def draw(self,d2s=None,col=None):
		"""This function causes the shape to render itself into the current GL context.
		d2s is a function of x,y which will convert data coordinates to screen
		coordinates. For data coordinate shapes, only the positional information
		is in data coordinates, font size, line width, etcare in screen units.
		col can be used to override the current shape's color."""
#		global glut_inited
		s=self.shape
		
		if s[0]=="hidden": return
		
		if col==None: col=self.shape[1:4]
		if d2s==None : d2s=shidentity
		
		v=d2s(s[4],s[5])
		v2=d2s(s[4]+1,s[5]+1)
		sc=v2[0]-v[0]
		
		if  self.isanimated:
#			print self.blend, "was the blend value"
			GL.glEnable(GL.GL_BLEND);
			depth_testing_was_on = GL.glIsEnabled(GL.GL_DEPTH_TEST);
			GL.glDisable(GL.GL_DEPTH_TEST);
			try:GL.glBlendEquation(GL.GL_FUNC_ADD)
			except:pass
			GL.glBlendFunc(GL.GL_SRC_ALPHA,GL.GL_ONE_MINUS_SRC_ALPHA);
			
			col=[self.shape[1],self.shape[2],self.shape[3],self.blend]
		
		if s[0]=="rect":
			assert self.shape[8] >= 0
			GL.glLineWidth(s[8])
			#if  self.isanimated:
				#v1 = d2s(s[4],s[5])
				#v2 = d2s(s[6],s[7])
				#dx = (v1[0]+v2[0])/2.0
				#dy = (v1[1]+v2[1])/2.0
				#GL.glPushMatrix()
				#GL.glTranslate(dx,dy,0)
				#GL.glScale(self.blend,self.blend,1)
				#GL.glTranslate(-dx,-dy,0)
			
			GL.glBegin(GL.GL_LINE_LOOP)
			GL.glColor(*col)
			GL.glVertex(*d2s(s[4],s[5]))
			GL.glVertex(*d2s(s[6],s[5]))
			GL.glVertex(*d2s(s[6],s[7]))
			GL.glVertex(*d2s(s[4],s[7]))
			GL.glEnd()
		
		if s[0]=="rect4point":
			assert self.shape[12] >= 0
			GL.glLineWidth(s[12])
			
			GL.glBegin(GL.GL_LINE_LOOP)
			GL.glColor(*col)
			GL.glVertex(*d2s(s[4],s[5]))
			GL.glVertex(*d2s(s[6],s[7]))
			GL.glVertex(*d2s(s[8],s[9]))
			GL.glVertex(*d2s(s[10],s[11]))
			GL.glEnd()
			
		elif s[0]=="rectline":
			#This makes a rectangle based on a width and two points
			#The two points are the endpoints for the longer axis of symmetry
			#Those two points implicitly define the length
			#The width is a parameter and together with the two coordinates, defines the rectangle
			
			assert s[4] != s[6] or s[5] != s[7], "The endpoints cannot coincide"
			assert self.shape[8] >= 0, "The width must be positive"
			assert self.shape[9] >= 0, "The line-width must be positive"
			
			pt0 = (s[4], s[5]) #first coordinate
			pt1 = (s[6], s[7]) #second coordinate
			width = s[8] #width
			GL.glLineWidth(s[9])
			#l_vect = (pt1[0]-pt0[0], pt1[1]-pt2[0]) #vector parallel to a longer side -- length vector
			w_vect = ( -(pt1[1]-pt0[1]), pt1[0]-pt0[0] ) #vector parallel to a shorter side -- width vector
			mag = sqrt(w_vect[0]**2 + w_vect[1]**2)
			w_uvect = (w_vect[0]/mag, w_vect[1]/mag)  #unit vector parallel to a short side
			
			#vertices - add/subtract a vector of half the box width with w_uvect direction
			v1 = ( pt0[0]-w_uvect[0]*width/2.0, pt0[1]-w_uvect[1]*width/2.0 )
			v2 = ( pt0[0]+w_uvect[0]*width/2.0, pt0[1]+w_uvect[1]*width/2.0 )
			v3 = ( pt1[0]+w_uvect[0]*width/2.0, pt1[1]+w_uvect[1]*width/2.0 )
			v4 = ( pt1[0]-w_uvect[0]*width/2.0, pt1[1]-w_uvect[1]*width/2.0 )

			#rectangle
			GL.glLineWidth(1)
			GL.glBegin(GL.GL_LINE_LOOP)
			GL.glColor(*col)
			GL.glVertex(*d2s(*v1))
			GL.glVertex(*d2s(*v2))
			GL.glVertex(*d2s(*v3))
			GL.glVertex(*d2s(*v4))
			GL.glEnd()
			#line
			GL.glBegin(GL.GL_LINES)
			GL.glColor(*col)
			GL.glVertex(*d2s(*pt0))
			GL.glVertex(*d2s(*pt1))
			GL.glEnd()
					
		elif s[0]=="rectpoint":
			assert self.shape[8] >= 0
			GL.glLineWidth(s[8])
			GL.glPointSize(s[8])
			
			#rectangle
			GL.glBegin(GL.GL_LINE_LOOP)
			GL.glColor(*col)
			GL.glVertex(*d2s(s[4],s[5]))
			GL.glVertex(*d2s(s[6],s[5]))
			GL.glVertex(*d2s(s[6],s[7]))
			GL.glVertex(*d2s(s[4],s[7]))
			GL.glEnd()
			
			#midpoint
			GL.glBegin(GL.GL_POINTS)
			p1 = d2s(s[4],s[5])
			p2 = d2s(s[6],s[7])
			GL.glVertex((p1[0]+p2[0])/2,(p1[1]+p2[1])/2)
			GL.glEnd()
			
			
		elif s[0]=="rcirclepoint":
			assert self.shape[8] >= 0
			GL.glLineWidth(s[8])
			GL.glPointSize(s[8])
			
			#midpoint
			GL.glBegin(GL.GL_POINTS)
			p1 = d2s(s[4],s[5])
			p2 = d2s(s[6],s[7])
			GL.glVertex((p1[0]+p2[0])/2,(p1[1]+p2[1])/2)
			GL.glEnd()
			v2=d2s(s[6],s[7])
			GL.glPushMatrix()
			GL.glColor(*col)
			
			#circle inscribed in the rectangle
			GL.glTranslate((v[0]+v2[0])/2.0,(v[1]+v2[1])/2.0,0)
			GL.glScalef((v2[0]-v[0])/2.0,(v2[1]-v[1])/2.0,1.0)
			GL.glCallList(EMShape.dlists)
			GL.glPopMatrix()
			
		elif s[0]=="rcircle":
			v2=d2s(s[6],s[7])
			assert self.shape[8] >= 0
			GL.glLineWidth(s[8])
			
			GL.glPushMatrix()
			GL.glColor(*col)
			
			GL.glTranslate((v[0]+v2[0])/2.0,(v[1]+v2[1])/2.0,0)
			GL.glScalef((v2[0]-v[0])/2.0,(v2[1]-v[1])/2.0,1.0)
			GL.glCallList(EMShape.dlists)
			GL.glPopMatrix()
		
		elif s[0]=="point":
			GL.glColor(*col)
			GL.glPointSize(s[6])
			p1 = d2s(s[4],s[5])
			GL.glBegin(GL.GL_POINTS)
			GL.glVertex(p1[0],p1[1])
			GL.glEnd()
			
		elif s[0]=="line":
#			print "A line ",s[4],s[5],s[6],s[7]
#			print "A line ",d2s(s[4],s[5]),d2s(s[6],s[7])
			GL.glColor(*col)
			assert self.shape[8] >= 0
			GL.glLineWidth(s[8])
			GL.glBegin(GL.GL_LINES)
			GL.glVertex(*d2s(s[4],s[5]))
			GL.glVertex(*d2s(s[6],s[7]))
			GL.glEnd()
		
		elif s[0]=="mask":
#			print "A line ",s[4],s[5],s[6],s[7]
#			print "A line ",d2s(s[4],s[5]),d2s(s[6],s[7])
			GL.glColor(*col)
			GL.glBegin(GL.GL_TRIANGLE_STRIP)
			GL.glVertex(*d2s(s[4],s[5]))
			GL.glVertex(*d2s(s[6],s[7]))
			GL.glVertex(*d2s(s[8],s[9]))
			GL.glVertex(*d2s(s[10],s[11]))
			GL.glVertex(*d2s(s[12],s[13]))
			GL.glVertex(*d2s(s[14],s[15]))
			GL.glVertex(*d2s(s[16],s[17]))
			GL.glVertex(*d2s(s[18],s[19]))
			#First Vertex
			GL.glVertex(*d2s(s[4],s[5]))
			GL.glVertex(*d2s(s[6],s[7]))
			GL.glEnd()
		
		elif s[0]=="linemask":
			assert self.shape[12] >= 0
			GL.glLineWidth(s[12])
			GL.glColor(*col)
			GL.glBegin(GL.GL_LINE_LOOP)
			GL.glVertex(*d2s(s[4],s[5]))
			GL.glVertex(*d2s(s[6],s[7]))
			GL.glVertex(*d2s(s[8],s[9]))
			GL.glVertex(*d2s(s[10],s[11]))
			GL.glEnd()
			
		elif s[0]=="label":
			if len(col)==4 : col=col[:3]
			GL.glPushMatrix()
			if EMShape.font_renderer != None :
				GL.glPushAttrib(GL.GL_ALL_ATTRIB_BITS)
				GL.glTranslate(v[0],v[1],.2)
				EMShape.font_renderer.set_face_size(int(s[7]))
				if s[8]<0 :
					# try for a sensible background color
					if col[0]+col[1]+col[2]>.4 and col[0]+col[1]+col[2]<.6 : bgcol=(0.0,0.0,0.0)
					else : bgcol = (1.0-col[0],1.0-col[1],1.0-col[2])
					bbox = EMShape.font_renderer.bounding_box(s[6])
					GLUtil.mx_bbox(bbox,col,bgcol)
				GL.glEnable(GL.GL_TEXTURE_2D)
				GL.glTexEnvi (GL.GL_TEXTURE_ENV, GL.GL_TEXTURE_ENV_MODE, GL.GL_REPLACE)
				GL.glColor(*col)
				EMShape.font_renderer.render_string(s[6])
				GL.glPopAttrib()
			else :
				pass
			#if s[8]<0 :
				#GL.glColor(1.,1.,1.)
				#GL.glTranslate(v[0],v[1],0)
				#GL.glScalef(s[7]/100.0/sc,s[7]/100.0/sc,s[7]/100.0/sc)
				#GL.glLineWidth(-s[8])
				#w=104.76*len(s[6])
				#GL.glBegin(GL.GL_QUADS)
				#GL.glVertex(-10.,-33.0)
				#GL.glVertex(w+10.,-33.0)
				#GL.glVertex(w+10.,119.05)
				#GL.glVertex(-10.,119.05)
				#GL.glEnd()
				#GL.glColor(*col)
				#for i in s[6]:
					#GLUT.glutStrokeCharacter(GLUT.GLUT_STROKE_MONO_ROMAN,ord(i))
			#else:
				#GL.glColor(*col)
				#GL.glTranslate(v[0],v[1],0)
##				GL.glScalef(s[7]/100.0,s[7]/100.0,s[7]/100.0)
				#GL.glScalef(s[7]/100.0/sc,s[7]/100.0/sc,s[7]/100.0/sc)
				#GL.glLineWidth(fabs(s[8]))
				#for i in s[6]:
					#GLUT.glutStrokeCharacter(GLUT.GLUT_STROKE_ROMAN,ord(i))
			GL.glPopMatrix()
			
		elif s[0]=="circle":
#			print s[6],v,v2
			GL.glPushMatrix()
			GL.glColor(*col)
			GL.glLineWidth(s[7])
			GL.glTranslate(v[0],v[1],0)
			GL.glScalef(s[6]*(v2[0]-v[0]),s[6]*(v2[1]-v[1]),1.0)
			GL.glCallList(EMShape.dlists)
			GL.glPopMatrix()
		elif s[0]=="ellipse":
#			print s[6],v,v2
			GL.glPushMatrix()
			GL.glColor(*col)
			GL.glLineWidth(s[9])
			GL.glTranslate(v[0],v[1],0)
			GL.glScalef((v2[0]-v[0]),(v2[1]-v[1]),1.0)
			GL.glRotatef(s[8],0,0,1.0)
			GL.glScalef(s[6],s[7],1.0)
			GL.glCallList(EMShape.dlists)
			GL.glPopMatrix()
		else:
			#mx=GL.glGetFloatv(GL.GL_MODELVIEW_MATRIX)
			GL.glPopMatrix()
			GL.glPushMatrix()
			if s[0]=="scrrect":
				x1 = int( round(s[4]) )
				y1 = int( round(s[5]) )
				x2 = int( round(s[6]) )
				y2 = int( round(s[7]) )
				GL.glLineWidth(s[8])
				GL.glBegin(GL.GL_LINE_LOOP)
				GL.glColor(*col)
				GL.glVertex(x1,y1)
				GL.glVertex(x2,y1)
				GL.glVertex(x2,y2)
				GL.glVertex(x1,y2)
				GL.glEnd()
			elif s[0]=="scrline":
				x1 = int( round(s[4]) )
				y1 = int( round(s[5]) )
				x2 = int( round(s[6]) )
				y2 = int( round(s[7]) )
				GL.glColor(*col)
				assert self.shape[8] >= 0
				GL.glLineWidth(s[8])
				GL.glBegin(GL.GL_LINES)
				GL.glVertex(x1,y1)
				GL.glVertex(x2,y2)
				GL.glEnd()
			elif s[0]=="scrlabel":
				if len(col)==4 : col=col[:3]
				x1 = int( round(s[4]) )
				y1 = int( round(s[5]) )
#				if s[8]<0 :
				if EMShape.font_renderer != None :
					GL.glPushAttrib(GL.GL_ALL_ATTRIB_BITS)
					GL.glTranslate(x1,y1,.2)
					if s[8]<0 :
						# try for a sensible background color
						if col[0]+col[1]+col[2]>.4 and col[0]+col[1]+col[2]<.6 : bgcol=(0.0,0.0,0.0)
						else : bgcol = (1.0-col[0],1.0-col[1],1.0-col[2])
						bbox = EMShape.font_renderer.bounding_box(s[6])
						GLUtil.mx_bbox(bbox,col,(1.0-col[0],1.0-col[1],1.0-col[2]))
					GL.glEnable(GL.GL_TEXTURE_2D)
					GL.glTexEnvi (GL.GL_TEXTURE_ENV, GL.GL_TEXTURE_ENV_MODE, GL.GL_REPLACE)
					GL.glColor(*col)
					EMShape.font_renderer.render_string(s[6])
					GL.glPopAttrib()
				else :
					pass
				
					#else:
						#GL.glColor(1.,1.,1.)
						#GL.glTranslate(x1,y1,0)
						##GL.glScalef(s[7]/1500.0/sc,s[7]/1500.0/sc,s[7]/1500.0/sc)
						#GL.glScalef(s[7]/1500.0,s[7]/1500.0,1)
						#GL.glLineWidth(-s[8])
						#w=104.76*len(text)
						#GL.glBegin(GL.GL_QUADS)
						#GL.glVertex(-10.,-33.0)
						#GL.glVertex(w+10.,-33.0)
						#GL.glVertex(w+10.,119.05)
						#GL.glVertex(-10.,119.05)
						#GL.glEnd()
						#GL.glTranslate(0,0,.1)
						#GL.glColor(*col)
						#for i in text:
							#GLUT.glutStrokeCharacter(GLUT.GLUT_STROKE_MONO_ROMAN,ord(i))
				#else:
					#GL.glColor(*col)
					#GL.glTranslate(x1,y1,-1)
	##				GL.glScalef(y2/100.0,y2/100.0,y2/100.0)
					#GL.glScalef(s[7]/1500.0/sc,s[7]/1500.0/sc,1)
					#GL.glLineWidth(fabs(s[8]))
					#for i in text:
						#GLUT.glutStrokeCharacter(GLUT.GLUT_STROKE_ROMAN,ord(i))
			elif s[0]=="scrcircle":
				x1 = int( round(s[4]) )
				y1 = int( round(s[5]) )

				GL.glColor(*col)
				GL.glLineWidth(s[7])
				GL.glTranslate(x1,y1,0)
				GL.glScalef(s[6],s[6],s[6])
				GL.glCallList(EMShape.dlists)
			#GL.glLoadMatrixf(mx)
		
		if self.isanimated:
			GL.glDisable( GL.GL_BLEND);
			if (depth_testing_was_on):
				GL.glEnable(GL.GL_DEPTH_TEST)			
	
	def getShape(self):
		return self.shape
	
	def setShape(self,shape):
		"""sets the shape to a new tuple/list"""
		self.shape=list(shape)
	
	def setColor(self,r,g=-1,b=-1):
		"""Sets the color of the shape. You may call this either as s.setColor(r,g,b) or s.setcolor([r,g,b])"""
		
		try:
			if g<0 :
				self.shape[1]=r[0]
				self.shape[2]=r[1]
				self.shape[3]=r[2]
			else:
				self.shape[1]=r
				self.shape[2]=g
				self.shape[3]=b
		except:
			print "Invalid color set in shape ",self.shape
			return
			
	def translate(self,dx,dy):
		"""This translates the shape without resizing it"""
		self.shape[4]+=dx
		self.shape[5]+=dy
		if self.shape[0] in ("rect","rectpoint","rectline","line","rcircle","rcirclepoint","scrrect","scrline") :
			self.shape[6]+=dx
			self.shape[7]+=dy
		
	def incblend(self):
		self.blend += self.blendinc
		if self.blend >= 1:
			self.isanimated = False
			return 1
		elif self.blend <= 0:
			# remove me please
			return 2
		else: return 0
		
	def set_blend(self,blend):
		if not self.isanimated: return
		self.blend = blend
		if self.blend >= 1: self.isanimated = False
		
	def is_animated(self):
		return self.isanimated
	def collision(self, x, y, fuzzy=False):
		s = self.shape
		
		if s[0] == "rect":
			warnings.warn("Not yet implemented.")
		elif s[0] == "rectpoint":
			warnings.warn("Not yet implemented.")
		elif s[0] == "rectline":
			#      0     1  2  3  4   5   6   7   8
			# s = [type, r, g, b, x1, y1, x2, y2, width]
			centroid = ( (s[4]+s[6])/2.0, (s[5]+s[7])/2.0 )
			l_vect = (s[6]-s[4], s[7]-s[5])
			length = sqrt(l_vect[0]**2+l_vect[1]**2)
			l_uvect = (l_vect[0]/length, l_vect[1]/length)
			w_uvect = (l_uvect[1], -l_uvect[0])
			width = s[8]
			#New coordinate system (w, l) with origin at centroid
			translated = (x-centroid[0], y-centroid[1]) #translate the origin to the centroid
			w = translated[0]*w_uvect[0] + translated[1]*w_uvect[1] #projection onto w_vect
			l = translated[0]*l_uvect[0] + translated[1]*l_uvect[1] #projection onto l_vect
			if abs(w) <= width/2.0:
				if fuzzy and abs(l) <= 5*length/8.0: #fuzzy = True includes area L/8 units beyond the ends
					return True
				elif abs(l) <= length/2.0:
					return True
				else:
					return False
			else:
				return False
		elif s[0] == "rcircle":
			warnings.warn("Not yet implemented.")
		elif s[0] == "rcirclepoint":
			warnings.warn("Not yet implemented.")
		elif s[0] == "line":
			warnings.warn("Not yet implemented.")
		elif s[0] == "label":
			warnings.warn("Not yet implemented.")
		elif s[0] == "circle":
			warnings.warn("Not yet implemented.")
		elif s[0] == "scrrect":
			warnings.warn("Not yet implemented.")
		elif s[0] == "scrline":
			warnings.warn("Not yet implemented.")
		elif s[0] == "scrlabel":
			warnings.warn("Not yet implemented.")
		elif s[0] == "scrcircle":
			warnings.warn("Not yet implemented.")
		elif s[0] == "point":
			warnings.warn("Not yet implemented.")
		else:
			raise LookupError
	def control_pts(self):
		s = self.shape
		
		if s[0] == "rect":
			warnings.warn("Not yet implemented.")
		elif s[0] == "rectpoint":
			warnings.warn("Not yet implemented.")
		elif s[0] == "rectline":
			#      0     1  2  3  4   5   6   7   8
			# s = [type, r, g, b, x1, y1, x2, y2, width]
			return ( (s[4],s[5]), (s[6],s[7]), ((s[4]+s[6])/2.0, (s[5]+s[7])/2.0) )			
		elif s[0] == "rcircle":
			warnings.warn("Not yet implemented.")
		elif s[0] == "rcirclepoint":
			warnings.warn("Not yet implemented.")
		elif s[0] == "line":
			warnings.warn("Not yet implemented.")
		elif s[0] == "label":
			warnings.warn("Not yet implemented.")
		elif s[0] == "circle":
			warnings.warn("Not yet implemented.")
		elif s[0] == "scrrect":
			warnings.warn("Not yet implemented.")
		elif s[0] == "scrline":
			warnings.warn("Not yet implemented.")
		elif s[0] == "scrlabel":
			warnings.warn("Not yet implemented.")
		elif s[0] == "scrcircle":
			warnings.warn("Not yet implemented.")
		elif s[0] == "point":
			warnings.warn("Not yet implemented.")
		else:
			raise LookupError
	def control_pt_min_distance(self, x, y):
		control_points = self.control_pts()
		squared_distances = [point[0]**2+point[1]**2 for point in control_points]
		return sqrt(min(squared_distances))

class EMShapeList(list):
	def collisions(self,x,y,fuzzy=False):
		"""
		This returns a list of shapes that enclose the point (x,y).	If fuzzy
		is set to True, the point need not be exactly inside the shape.
		"""
		ret = []
		for shape in self:
			if shape.collision(x,y,fuzzy):
				ret.append[shape]
		return ret
	
	def closest_collision(self, x, y, fuzzy=False):
		"""
		This searches the list returned by collisions() to find the shape
		with a control point closest to x,y. The parameter fuzzy is used by
		collisions(). If there are no collisions, this returns none.
		"""
		shapes = self.collisions(x, y, fuzzy)
		if shapes:
			min_dist = None
			for shape in shapes:
				dist = shape.control_pt_min_distance(x,y)
				if min_dist == None:
					min_dist = dist
					closest_shape = shape
				elif dist < min_dist:
					min_dist = dist
					closest_shape = shape
			return closest_shape		
		else:
			return None

class EMShapeDict(dict):
	def collisions(self,x,y,fuzzy=False):
		"""
		This returns a list of keys to the shapes that enclose the point (x,y).
		If fuzzy is set to True, the point need not be exactly inside the shape.
		"""
		ret = []
		for k in self:
			if self.get(k).collision(x,y, fuzzy):
				ret.append(k)
		return ret
	
	def closest_collision(self,x,y,fuzzy=False):
		"""
		This searches the list returned by collisions() to find the shape
		with a control point closest to x,y. The key to that shape is returned.
		The parameter fuzzy is used by collisions(). If there are no collisions, 
		this returns None.
		"""
		shape_keys = self.collisions(x, y, fuzzy)
		if shape_keys:
			min_dist = None
			for key in shape_keys:
				shape = self.get(key)
				
				dist = shape.control_pt_min_distance(x,y)
				if min_dist == None:
					min_dist = dist
					closest_shape_key = key
				elif dist < min_dist:
					min_dist = dist
					closest_shape_key = key
			return closest_shape_key
		else:
			return None