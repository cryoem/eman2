#!/usr/bin/env python

#
# Author: Steven Ludtke, 04/10/2003 (sludtke@bcm.edu)
# and David Woolford 10/26/2007 (woolford@bcm.edu)
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

from EMAN2 import *
from OpenGL import GL, GLU, GLUT
from OpenGL.GL import *
from OpenGL.GLU import *
from OpenGL.GLUT import *
from PyQt4 import QtCore, QtGui, QtOpenGL
from PyQt4.QtCore import Qt
from emglobjects import Camera2, get_default_gl_colors, EMViewportDepthTools2, get_RGB_tab, get_gl_lights_vector, init_glut, EM3DModel
from emimageutil import EMTransformPanel # for EMLightsInspector
from math import *
from time import *
from valslider import ValSlider
import weakref # for EMLightsInspector



MAG_INCREMENT_FACTOR = 1.1

class EMLightsDrawer:
	'''
	Base clase, works with EMLightsInspectorBase
	'''
	def __init__(self):
		
		self.gl_lights = get_gl_lights_vector() # [GL_LIGHTS0, GL_LIGHTS1, ... GL_LIGHTS<MAX>]
		self.colors = get_default_gl_colors()  # a dictionary of dictionaries
		
		self.gq=gluNewQuadric()
		gluQuadricDrawStyle(self.gq,GLU_FILL)
		gluQuadricNormals(self.gq,GLU_SMOOTH)
		gluQuadricOrientation(self.gq,GLU_OUTSIDE)
		gluQuadricTexture(self.gq,GL_FALSE)
		
		self.mouse_modes = [None, "directional", "point source"]
		self.current_mouse_mode = None
		self.current_light = GL_LIGHT0 # change this to change which light is drawn
		self.display_lights = False
		self.ball_dl = 0 # makes a funny kind of ball
		self.light_dl = 0
		self.cylinderdl = 0
		self.arc_t = 16
		self.radius = 10 #  This is an important variable. Set this value to adapt the way the lights appear in your scene.
		
		self.inspector = None
		self.mpressx = None
		self.mpressy = None
	
	def set_gl_widget(self, gl_widget):
		self.gl_widget = weakref.ref(gl_widget)
	def set_current_light(self,light):
		self.current_light = light
		if self.current_mouse_mode != None:
			pos = glGetLightfv(self.current_light,GL_POSITION)
			if pos[3] == 0:
				self.current_mouse_mode = "directional"
			else:
				self.current_mouse_mode = "point source"
			
	def generate_arcs(self,points,n,halt=0):
		for i in range(0,n-halt):
			p1 = points[i]
			if ( i == n-1 ): p2 = points[0]
			else: p2 = points[i+1]
			angle = acos(p2.dot(p1))
			sinangle = sin(angle)
			prev = self.radius*Vec3f(p1[0],p1[1],p1[2])
			for t in range(1, self.arc_t+1):
				timeangle = float(t)/float(self.arc_t)*angle
				p1Copy = self.radius*Vec3f(p1[0],p1[1],p1[2])
				p2Copy = self.radius*Vec3f(p2[0],p2[1],p2[2])
				next = (sin(angle-timeangle)*p1Copy + sin(timeangle)*p2Copy)/sinangle
				
				self.cylinder_to_from(next,prev)
				prev = Vec3f(next[0],next[1],next[2])
				
	
	def cylinder_to_from(self,next,prev):
		dx = next[0] - prev[0]
		dy = next[1] - prev[1]
		dz = next[2] - prev[2]
		
		length = sqrt(dx**2 + dy**2 + dz**2)
		
		if length == 0: return
		
		alt = acos(dz/length)*180.0/pi
		phi = atan2(dy,dx)*180.0/pi
		
		glPushMatrix()
		glTranslatef(prev[0],prev[1],prev[2] )
		#print "positioned at", prev[0],prev[1],prev[2]
		glRotatef(90+phi,0,0,1)
		glRotatef(alt,1,0,0)
		
		glMaterial(GL_FRONT, GL_AMBIENT, self.colors["emerald"]["ambient"])
		glMaterial(GL_FRONT, GL_DIFFUSE, self.colors["emerald"]["diffuse"])
		glMaterial(GL_FRONT, GL_SPECULAR, self.colors["emerald"]["specular"])
		glMaterial(GL_FRONT, GL_SHININESS, self.colors["emerald"]["shininess"])
		glutSolidTorus(.05,0.25,16,16)
		glScalef(0.2,0.2,length)
		glMaterial(GL_FRONT, GL_AMBIENT, self.colors["ruby"]["ambient"])
		glMaterial(GL_FRONT, GL_DIFFUSE, self.colors["ruby"]["diffuse"])
		glMaterial(GL_FRONT, GL_SPECULAR, self.colors["ruby"]["specular"])
		glMaterial(GL_FRONT, GL_SHININESS, self.colors["ruby"]["shininess"])
		glCallList(self.cylinderdl) #FIXME: this cylinder isn't visible on Mac OS
		glPopMatrix()

	
	
	
	def draw_lights(self):
		init_glut()

		
		if ( self.cylinderdl == 0 ):
			self.cylinderdl=glGenLists(1)
				
			glNewList(self.cylinderdl,GL_COMPILE)
			glPushMatrix()
			gluCylinder(self.gq,1.0,1.0,1.0,12,2)
			glPopMatrix()
				
			glEndList()
			
		if self.ball_dl == 0:
			self.ball_dl = glGenLists(1)
			glNewList(self.ball_dl,GL_COMPILE)
			glPushMatrix()
			glScale(self.radius/10.0,self.radius/10.0,self.radius/10.0)
			self.draw_light_cocoon()
			self.draw_inside_light()
			glPopMatrix()
			
			glMaterial(GL_FRONT, GL_AMBIENT, self.colors["obsidian"]["ambient"])
			glMaterial(GL_FRONT, GL_DIFFUSE, self.colors["obsidian"]["diffuse"])
			glMaterial(GL_FRONT, GL_SPECULAR, self.colors["obsidian"]["specular"])
			glMaterial(GL_FRONT, GL_EMISSION, self.colors["obsidian"]["emission"])
			glMaterial(GL_FRONT, GL_SHININESS, self.colors["obsidian"]["shininess"])
			
			glPushMatrix()
			n = 12
			
			for i in range(n):
				if i % 2 == 0: color = "gold"
				else: color = "silver"
				
				glMaterial(GL_FRONT, GL_AMBIENT, self.colors[color]["ambient"])
				glMaterial(GL_FRONT, GL_DIFFUSE, self.colors[color]["diffuse"])
				glMaterial(GL_FRONT, GL_SPECULAR, self.colors[color]["specular"])
				glMaterial(GL_FRONT, GL_EMISSION, self.colors[color]["emission"])
				glMaterial(GL_FRONT, GL_SHININESS, self.colors[color]["shininess"])
				rot = 180.0*i/(n-1)
				glPushMatrix()
				glTranslate(0,0,-.5)
				glRotate(rot,0,0,1)
				glRotate(90,0,1,0)
				glScale(self.radius/10.0,self.radius/10.0,self.radius/10.0)
				glutSolidTorus(.1,0.71,32,32)
				glPopMatrix()
			
			glPopMatrix()
			
			glEndList()
			
			
		if self.light_dl == 0:
			self.light_dl = glGenLists(1)
			glNewList(self.light_dl,GL_COMPILE)
			glPushMatrix()
			glTranslate(0,0,self.radius-1.42)
			glRotate(180,1,0,0)
			
			glCallList(self.ball_dl)
			
			glPopMatrix()
			
			glPushMatrix()
			points = [Vec3f(1,0,0),Vec3f(0,1,0),Vec3f(-1,0,0),Vec3f(0,-1,0)]
			self.generate_arcs(points,len(points))
			glPopMatrix()
			
			glPushMatrix()
			points = [Vec3f(1,0,0),Vec3f(0,1,0),Vec3f(-1,0,0)]
			glRotate(90,1,0,0)
			self.generate_arcs(points,len(points),1)
			glPopMatrix()
			glEndList()
		
		
		glPushMatrix()
		self.position_light()
		
		glPopMatrix()
	
	def position_light(self):
		if not self.gl_lights:
			self.gl_lights = get_gl_lights_vector()
		for light in self.gl_lights:
			
			if not glIsEnabled(light):continue
			pos = glGetLightfv(light,GL_POSITION)
			glPushMatrix()
			if pos[3] == 1:
				glLoadIdentity() # The position of a point source is stored in eye coordinates!
				sp = glGetLightfv(light,GL_SPOT_DIRECTION)
				sc = glGetLightfv(light,GL_SPOT_CUTOFF)
				#print sp
				#print sc
			
				v = Vec3f(float(sp[0]),float(sp[1]),float(sp[2]))
				if v.length() == 0:
					glPopMatrix()
					return
				t = Transform()
				t.set_rotation(v)
				d = t.get_rotation("eman")
				
				glTranslate(pos[0],pos[1],pos[2])
				#print pos[0],pos[1],pos[2]
				#print d
				#print d
				glRotate(-d["phi"],0,0,1)
				glRotate(-d["alt"],1,0,0)
				glRotate(-d["az"],0,0,1)


				glCallList(self.ball_dl)
				
				
			else:
				v = Vec3f(float(pos[0]),float(pos[1]),float(pos[2]))
				t = Transform()
				t.set_rotation(v)
				d = t.get_rotation("eman")
				
				#print d
				glRotate(-d["phi"],0,0,1)
				glRotate(-d["alt"],1,0,0)
				glRotate(-d["az"],0,0,1)
			
		
				glCallList(self.light_dl)
		
			glPopMatrix()
	def draw_light_cocoon(self):
		'''
		draw a nice looking torch
		'''

		n = 20 # height discretizations
		bot = 0.5 # square width at bottom
		top = 1.0 # square width at top
		
		dz = (top-bot)/(n-1)
		yellow = [1,1,0,0.5]
		glMaterial(GL_FRONT, GL_AMBIENT, yellow)
		glMaterial(GL_FRONT, GL_DIFFUSE, yellow)
		glMaterial(GL_FRONT, GL_SPECULAR, yellow)
		glMaterial(GL_FRONT, GL_EMISSION, [0,0,0,0])
		glMaterial(GL_FRONT, GL_SHININESS, 32)
		glColor(*yellow)
		
		glEnable(GL_BLEND)
		glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA)
		
		glPushMatrix()
		self.square_from_points([bot,bot,0],[top,top,1],[top,-top,1],[bot,-bot,0])
		glPopMatrix()
		glPushMatrix()
		glRotate(90,0,0,1)
		self.square_from_points([bot,bot,0],[top,top,1],[top,-top,1],[bot,-bot,0])
		glPopMatrix()
		glPushMatrix()
		glRotate(180,0,0,1)
		self.square_from_points([bot,bot,0],[top,top,1],[top,-top,1],[bot,-bot,0])
		glPopMatrix()
		glPushMatrix()
		glRotate(270,0,0,1)
		self.square_from_points([bot,bot,0],[top,top,1],[top,-top,1],[bot,-bot,0])
		glPopMatrix()
		glDisable(GL_BLEND)
		
	def square_from_points(self,p1,p2,p3,p4):
		'''
		draws a square from the given points, calculates the normal automatically
		assumes counter clockwise direction
		'''
		
		v1 = Vec3f(p1)
		v1.normalize()
		v2 = Vec3f(p2)
		v2.normalize()
		v3 = Vec3f(p4)
		v3.normalize()
		
		normal = (v2-v1).cross(v3-v1)
		
		glBegin(GL_QUADS)
		glNormal(normal[0],normal[1],normal[2])
		glVertex(*p1)
		glVertex(*p2)
		glVertex(*p3)
		glVertex(*p4)
		glEnd()
	
	
	def draw_inside_light(self):
		'''
		draw a nice looking torch
		'''

		n = 20 # height discretizations
		bot = 0.5 # square width at bottom
		top = 1.0 # square width at top
		
		dz = (top-bot)/(n-1)
		yellow = [1,1,0,1.0/n]
		glMaterial(GL_FRONT, GL_AMBIENT, yellow)
		glMaterial(GL_FRONT, GL_DIFFUSE, yellow)
		glMaterial(GL_FRONT, GL_SPECULAR, yellow)
		glMaterial(GL_FRONT, GL_EMISSION, [0,0,0,0])
		glMaterial(GL_FRONT, GL_SHININESS, 32)
		glColor(*yellow)
		
		glEnable(GL_BLEND)
		glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA)
		
		
		for i in range(n):
			width = i*dz + bot
			glPushMatrix()
			glTranslate(0,0,2*(width-bot))
			glScale(2*width,2*width,1)
			self.square()
			glPopMatrix()
		
		glDisable(GL_BLEND)
	
	def square(self):
		'''
		a 1x1 square centered on the origin looking along positive z
		'''
		glBegin(GL_QUADS)
		glNormal(0,0,1)
		glVertex(-0.5,-0.5,0)
		glVertex(0.5,-0.5,0)
		glVertex(0.5,0.5,0)
		glVertex(-0.5,0.5,0)
		glEnd()
		
	def light_manipulation_toggled(self,state):
		self.light_manip = state
		if state:
			pos = glGetLightfv(self.current_light,GL_POSITION)
			if pos[3] == 0:
				self.current_mouse_mode = "directional"
			else:
				self.current_mouse_mode = "point source"
		else:
			self.current_mouse_mode = None
			
	def show_lights(self,state):
		self.display_lights = state
		self.updateGL()
		
	def draw(self):
		if self.display_lights:
			self.draw_lights()
	def motion_translate_z_only(self,prev_x,prev_y,event):
		[dx,dy] = [event.x()-prev_x,prev_y-event.y()]
		dx /= 10.0
		dy /= 10.0
		d = abs(dx) + abs(dy)
		if dy > 0: d = -d 
		
		pos = glGetLightfv(self.current_light,GL_POSITION)
		
		test_pos = [0,0,0,1]
		glLightfv(self.current_light,GL_POSITION,test_pos)
		test_pos_out = glGetLightfv(self.current_light,GL_POSITION)
		# reset to the correction position
		pos = [ (pos[i] - test_pos_out[i]) for i in range(3)]
		pos[2] += d

		pos.append(1)
		glLightfv(self.current_light,GL_POSITION,pos)
		if self.inspector:
			self.inspector.set_positional_light_pos(pos)
		
	def motion_translate(self,prev_x,prev_y,event):
		
		[dx,dy] = [event.x()-prev_x,prev_y-event.y()]
		dx /= 10.0
		dy /= 10.0
		pos = glGetLightfv(self.current_light,GL_POSITION)
		
		test_pos = [0,0,0,1]
		glLightfv(self.current_light,GL_POSITION,test_pos)
		test_pos_out = glGetLightfv(self.current_light,GL_POSITION)
		# reset to the correction position
		pos = [ (pos[i] - test_pos_out[i]) for i in range(3)]
		pos[0] += dx
		pos[1] += dy
		
		pos.append(1)
		glLightfv(self.current_light,GL_POSITION,pos)
		if self.inspector:
			self.inspector.set_positional_light_pos(pos)
	def motion_rotate(self,x,y,fac=1.0):
		# this function implements mouse interactive rotation
		# [x,y] is the vector generating by the mouse movement (in the plane of the screen)
		# Rotation occurs about the vector 90 degrees to [x,y,0]
		# The amount of rotation is linealy proportional to the length of [x,y]
		
		if self.current_mouse_mode not in ["directional", "point source"]:
			raise NotImplementedError
		
		if ( x == 0 and y == 0): return
		
		theta = atan2(-y,x)

		rotaxis_x = -sin(theta)
		rotaxis_y = cos(theta)
		rotaxis_z = 0
		length = sqrt(x*x + y*y)
		# motiondull is a magic number - things rotate more if they are closer and slower if they are far away in this appproach
		# This magic number could be overcome using a strategy based on the results of get_render_dims_at_depth
		angle = 4*fac*length/pi # the four is just because I liked the feel of it
		
		t = Transform()
		quaternion = {}
		quaternion["omega"] = angle
		quaternion["n1"] = rotaxis_x
		quaternion["n2"] = rotaxis_y
		quaternion["n3"] = rotaxis_z
		quaternion["type"] = "spin"
		
		t.set_rotation( quaternion )
		
		if self.current_mouse_mode == "point source":
			dr = glGetLightfv(self.current_light,GL_SPOT_DIRECTION)
			v = Vec3f(float(dr[0]),float(dr[1]),float(dr[2]))
		elif self.current_mouse_mode == "directional":
			pos = glGetLightfv(self.current_light,GL_POSITION)
			v = Vec3f(float(pos[0]),float(pos[1]),float(pos[2]))

		torig = Transform()
		torig.set_rotation(v)
		
		
		t2 = t*torig
		if self.current_mouse_mode == "point source":
			new_dr = t2*Vec3f(0,0,1)
			new_dr_list = [new_dr[i] for i in range(3)]
	
			glLightfv(self.current_light,GL_SPOT_DIRECTION,new_dr_list)
			if self.inspector != None: self.inspector.set_positional_light_dir(new_dr_list)
		elif self.current_mouse_mode == "directional":
			new_pos = t2*Vec3f(0,0,1)
			new_pos_list = [new_pos[i] for i in range(3)]
			new_pos_list.append(0)
			glLightfv(self.current_light,GL_POSITION,new_pos_list)
			if self.inspector != None: self.inspector.set_directional_light_dir(new_pos_list)

	def mousePressEvent(self,event):
		if self.current_mouse_mode != None:
			self.mpressx = event.x()
			self.mpressy = event.y()
			self.updateGL()
	
	def mouseMoveEvent(self,event):
		if self.current_mouse_mode == "point source":
			if event.buttons()&Qt.RightButton and event.modifiers()&Qt.ShiftModifier:
				
				self.motion_translate_z_only(self.mpressx, self.mpressy,event)
				self.mpressx = event.x()
				self.mpressy = event.y()
	
			elif event.buttons()&Qt.RightButton:
				self.motion_translate(self.mpressx, self.mpressy,event)
				
				self.mpressx = event.x()
				self.mpressy = event.y()
			else:
				self.motion_rotate(self.mpressx - event.x(), self.mpressy - event.y())
				self.mpressx = event.x()
				self.mpressy = event.y()
			self.updateGL()
		elif self.current_mouse_mode == "directional":
			self.motion_rotate(self.mpressx - event.x(), self.mpressy - event.y())
			self.mpressx = event.x()
			self.mpressy = event.y()
			self.updateGL()
	
	def mouseReleaseEvent(self,event):
		if self.current_mouse_mode != None:
			pass
	
class EMLights(EMLightsDrawer,EM3DModel):
	def eye_coords_dif(self,x1,y1,x2,y2,mdepth=True):
		return self.vdtools.eye_coords_dif(x1,y1,x2,y2,mdepth)
	
	def __init__(self, gl_widget):
		EM3DModel.__init__(self, gl_widget)
		EMLightsDrawer.__init__(self)
		self.display_lights = True
		self.set_gl_widget(gl_widget)

		self.init()
		self.initialized = True
		
		self.cam=Camera2(self)
		
		self.brightness = 0
		self.contrast = 10
		self.glcontrast = 1.0
		self.glbrightness = 0.0
		self.rank = 1
		self.inspector=None
		self.currentcolor = "emerald"
		self.vdtools = None	
		self.highresspheredl = 0 # display list id
		self.draw_dl = 0

	def set_gl_widget(self, gl_widget):
		EMLightsDrawer.set_gl_widget(self, gl_widget)
		self.vdtools = EMViewportDepthTools2(self.gl_widget())
		self.gl_widget().cam.default_z = -25	 # this is me hacking
		self.gl_widget().cam.cam_z = -25 # this is me hacking		
	def get_type(self):
		return "lights"

	def render(self):
		init_glut()
		#if (not isinstance(self.data,EMData)): return
		lighting = glIsEnabled(GL_LIGHTING)
		cull = glIsEnabled(GL_CULL_FACE)
		depth = glIsEnabled(GL_DEPTH_TEST)
		polygonmode = glGetIntegerv(GL_POLYGON_MODE)

		glDisable(GL_CULL_FACE)
		glEnable(GL_DEPTH_TEST)
		
		if ( self.wire ):
			glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
		else:
			glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
		
		if self.light:
			glEnable(GL_LIGHTING)
		else:
			glDisable(GL_LIGHTING)

		glPushMatrix()
		self.cam.position()
			
		glShadeModel(GL_SMOOTH)

		glStencilFunc(GL_EQUAL,self.rank,0)
		glStencilOp(GL_KEEP,GL_KEEP,GL_REPLACE)
		glMaterial(GL_FRONT, GL_AMBIENT, self.colors[self.currentcolor]["ambient"])
		glMaterial(GL_FRONT, GL_DIFFUSE, self.colors[self.currentcolor]["diffuse"])
		glMaterial(GL_FRONT, GL_SPECULAR, self.colors[self.currentcolor]["specular"])
		glMaterial(GL_FRONT, GL_EMISSION, self.colors[self.currentcolor]["emission"])
		glMaterial(GL_FRONT, GL_SHININESS, self.colors[self.currentcolor]["shininess"])
		glColor(self.colors[self.currentcolor]["ambient"])
		
		glEnable(GL_NORMALIZE)
		#HERE
		
		if self.highresspheredl == 0:
			self.highresspheredl=glGenLists(1)
				
			glNewList(self.highresspheredl,GL_COMPILE)
			gluSphere(self.gq,.5,16,16)
			glEndList()
		
		glScale(2,2,2)
		
		if ( self.draw_dl == 0 ):
			self.draw_dl=glGenLists(1)
				
			glNewList(self.draw_dl,GL_COMPILE)
			z = [-10,-4,-2,2,4,10]
			t = [[2,0,0],[-2,0,0],[0,2,0],[0,-2,0],[0,0,0]]
			for r in z:
				
				for s in t:
					glPushMatrix()
					glTranslate(0,0,r)
					glTranslate(*s)
					v = Util.get_irand(0,7)
					if v == 0:
						glCallList(self.highresspheredl)
					elif v == 1:
						try:
							glutSolidSierpinskiSponge(3,0,1.0)
						except: pass # maybe some versions don't have it?
						#glutSolidRhombicDodecahedron()
					elif v == 2:
						glScale(0.5,0.5,0.5)
						glutSolidDodecahedron()
					elif v == 3:
						#glScale(0.5,0.5,0.5)
						glutSolidIcosahedron()
					elif v == 4:
						#glScale(0.5,0.5,0.5)
						glutSolidTetrahedron()
						
					elif v == 5:
						#glScale(0.5,0.5,0.5)
						glutSolidCube(1)
					
					elif v == 6:
						#glScale(0.5,0.5,0.5)
						glutSolidTorus(.25,.5,8,32)
						
					elif v == 7:
						#glScale(0.5,0.5,0.5)
						glutSolidTeapot(0.5)
					
					glPopMatrix()
			glEndList(self.draw_dl)
			
		glCallList(self.draw_dl)
		
		glPopMatrix()
		EMLightsDrawer.draw(self)
		
		self.draw_bc_screen()


		if ( lighting ): glEnable(GL_LIGHTING)
		else: glDisable(GL_LIGHTING)
		if ( cull ): glEnable(GL_CULL_FACE)
		else: glDisable(GL_CULL_FACE)
		if ( depth ): glEnable(GL_DEPTH_TEST)
		else : glDisable(GL_DEPTH_TEST)
		
		if ( polygonmode[0] == GL_LINE ): glPolygonMode(GL_FRONT, GL_LINE)
		else: glPolygonMode(GL_FRONT, GL_FILL)
		if ( polygonmode[1] == GL_LINE ): glPolygonMode(GL_BACK, GL_LINE)
		else: glPolygonMode(GL_BACK, GL_FILL)
		
	def init(self):
		self.mmode = 0
		self.wire = False
		self.light = True
	
	def setInit(self):

		#self.cam.default_z = -1.25*32
		#self.cam.cam_z = -1.25*32
		
		if not self.inspector or self.inspector ==None:
			self.inspector=EMLightsInspector(self)
		
		self.load_colors()
		
	
	def setColor(self,val):
		#print val
		self.currentcolor = str(val)
		self.updateGL()
	
	def toggle_wire(self,val):
		self.wire = not self.wire
		self.updateGL()
		
	def toggle_light(self,val):
		self.light = not self.light
		self.updateGL()
	
	def update_inspector(self,t3d):
		if not self.inspector or self.inspector ==None:
			self.inspector=EMLightsInspector(self)
		self.inspector.update_rotations(t3d)
	
	def get_inspector(self):
		if not self.inspector : 
			self.inspector=EMLightsInspector(self)
			self.inspector.set_colors(self.colors,self.currentcolor)
		return self.inspector
	def mouseMoveEvent(self, event):
		if self.current_mouse_mode:
			EMLightsDrawer.mouseMoveEvent(self, event)
		else:
			EM3DModel.mouseMoveEvent(self, event)
	def mousePressEvent(self, event):
		if self.current_mouse_mode:
			EMLightsDrawer.mousePressEvent(self, event)
		else:
			EM3DModel.mousePressEvent(self, event)
	def mouseReleaseEvent(self, event):
		if self.current_mouse_mode:
			EMLightsDrawer.mouseReleaseEvent(self, event)
		else:
			EM3DModel.mouseReleaseEvent(self, event)


class EMLightsInspectorBase:
	'''
	Inherit from this if you want its functionality
	'''

	def __init__(self):
		self.gl_lights = get_gl_lights_vector()
		self.quiet = False


	def set_positional_light_pos(self,pos):
		self.quiet = True
		p = [self.light_x_pos,self.light_y_pos,self.light_z_pos]
		for i,w in enumerate(p):
			w.setValue(pos[i])
		
		self.quiet = False

	def set_positional_light_dir(self,direction):
		self.quiet = True
		p = [self.light_ps_xdir,self.light_ps_ydir,self.light_ps_zdir]
		for i,w in enumerate(p):
			w.setValue(direction[i])
		
		self.quiet = False

	def set_directional_light_dir(self,direction):
		
		self.quiet = True
		p = [self.light_x_dir,self.light_y_dir,self.light_z_dir]
		for i,w in enumerate(p):
			w.setValue(direction[i])
		
		self.quiet = False
	
	def update_light(self):
		if self.quiet: return
		
		l = self.get_current_light()
		if l == None: return
		
		attribs = ["r","g","b"]
		amb = [ getattr(self.light_ambient, a).getValue() for a in attribs]
		amb.append(1.0) # alpha
		dif = [ getattr(self.light_diffuse, a).getValue() for a in attribs]
		dif.append(1.0) # alpha
		spec = [ getattr(self.light_specular, a).getValue() for a in attribs]
		spec.append(1.0) # alpha
		
		glLightfv(l, GL_AMBIENT, amb)
		glLightfv(l, GL_DIFFUSE, dif)
		glLightfv(l, GL_SPECULAR, spec)
		
		pos = glGetLightfv(l,GL_POSITION)
		if pos[3] == 0:
			# directional light
			p = [self.light_x_dir,self.light_y_dir,self.light_z_dir]
			pos = [float(a.value()) for a in p]
			pos.append(0) # i.e. directional light
			glLightfv(l, GL_POSITION, pos)
		else:
			p = [self.light_x_pos,self.light_y_pos,self.light_z_pos]
			pos = [float(a.value()) for a in p]
			pos.append(1) # i.e.point light
			glLightfv(l, GL_POSITION, pos)
			
			d = [self.light_ps_xdir,self.light_ps_ydir,self.light_ps_zdir]
			dr = [float(a.value()) for a in d]
			glLightfv(l, GL_SPOT_DIRECTION, dr)
			
			glLightfv(l,GL_CONSTANT_ATTENUATION,self.const_atten.getValue())
			glLightfv(l,GL_LINEAR_ATTENUATION,self.linear_atten.getValue())
			glLightfv(l,GL_QUADRATIC_ATTENUATION,self.quad_atten.getValue())
			glLightfv(l,GL_SPOT_CUTOFF,self.spot_cutoff.getValue())
			glLightfv(l,GL_SPOT_EXPONENT,self.spot_exponent.getValue())

		
		self.target().updateGL()
		
	def get_light_tab(self):
		self.light_tab = QtGui.QWidget()
		light_tab = self.light_tab
		
		vbl = QtGui.QVBoxLayout(self.light_tab )
		vbl.setMargin(0)
		vbl.setSpacing(6)
		vbl.setObjectName("Lights")
		
		self.light_manip_check = QtGui.QCheckBox("Mouse moves lights")
		self.local_viewer_check = QtGui.QCheckBox("Local light model")
		self.local_viewer_check.setChecked(glGetInteger(GL_LIGHT_MODEL_LOCAL_VIEWER))
		show_lights = QtGui.QCheckBox("Show lights")
		show_lights.setChecked(self.target().display_lights)
		max_lights_label = QtGui.QLabel()
		max_lights_label.setText("Max lights : " + str(glGetInteger(GL_MAX_LIGHTS)))
		
		hdl_l = QtGui.QHBoxLayout()
		hdl_t = QtGui.QHBoxLayout()
		
		hdl_l.addWidget(self.light_manip_check)
		hdl_l.addWidget(self.local_viewer_check)
		hdl_t.addWidget(show_lights)
		hdl_t.addWidget(max_lights_label)
		
		vbl.addLayout(hdl_l)
		vbl.addLayout(hdl_t)
		
		self.light_tab_widget = QtGui.QTabWidget()
		
		
		self.light_tab_widget.addTab(self.get_directional_light_tab(), "Directional")
		self.light_tab_widget.addTab(self.get_pointsource_light_tab(), "Point source")
		
		vbl.addWidget(self.light_tab_widget)
		
		light_material_tab_widget = QtGui.QTabWidget()
		self.light_ambient = get_RGB_tab(self,"ambient")
		light_material_tab_widget.addTab(self.light_ambient, "Ambient")
		self.light_diffuse = get_RGB_tab(self,"diffuse")
		light_material_tab_widget.addTab(self.light_diffuse, "Diffuse")
		self.light_specular = get_RGB_tab(self,"specular")
		light_material_tab_widget.addTab(self.light_specular, "Specular")
		
		
		self.refresh_light_states()
		
		vbl.addWidget(light_material_tab_widget)

		QtCore.QObject.connect(self.light_ambient.r, QtCore.SIGNAL("valueChanged"), self.update_light)
		QtCore.QObject.connect(self.light_ambient.g, QtCore.SIGNAL("valueChanged"), self.update_light)
		QtCore.QObject.connect(self.light_ambient.b, QtCore.SIGNAL("valueChanged"), self.update_light)
		QtCore.QObject.connect(self.light_diffuse.r, QtCore.SIGNAL("valueChanged"), self.update_light)
		QtCore.QObject.connect(self.light_diffuse.g, QtCore.SIGNAL("valueChanged"), self.update_light)
		QtCore.QObject.connect(self.light_diffuse.b, QtCore.SIGNAL("valueChanged"), self.update_light)
		QtCore.QObject.connect(self.light_specular.r, QtCore.SIGNAL("valueChanged"), self.update_light)
		QtCore.QObject.connect(self.light_specular.g, QtCore.SIGNAL("valueChanged"), self.update_light)
		QtCore.QObject.connect(self.light_specular.b, QtCore.SIGNAL("valueChanged"), self.update_light)
		QtCore.QObject.connect(self.light_x_dir, QtCore.SIGNAL("valueChanged(double)"), self.update_light)
		QtCore.QObject.connect(self.light_y_dir, QtCore.SIGNAL("valueChanged(double)"), self.update_light)
		QtCore.QObject.connect(self.light_z_dir, QtCore.SIGNAL("valueChanged(double)"), self.update_light)
		QtCore.QObject.connect(self.light_manip_check, QtCore.SIGNAL("valueChanged(double)"), self.update_light)
		QtCore.QObject.connect(self.light_manip_check, QtCore.SIGNAL("stateChanged(int)"), self.target().light_manipulation_toggled)
		QtCore.QObject.connect(show_lights, QtCore.SIGNAL("stateChanged(int)"), self.target().show_lights)
		QtCore.QObject.connect(self.local_viewer_check, QtCore.SIGNAL("stateChanged(int)"), self.local_viewer_checked)
		#QtCore.QObject.connect(self.light_w_pos, QtCore.SIGNAL("valueChanged(int)"), self.update_light)
	 
		return light_tab
	
	def get_current_light(self):
		selected_items = self.light_list.selectedItems()
		if len(selected_items) == 0: selected_items = self.point_light_list.selectedItems()
			
		if len(selected_items) == 0:
			return None
		if len(selected_items) > 1:
			print "write a message box to say there is a function that assumes only one light is selected"
			return None
		
		item = selected_items[0]

		#print "selected lights is",item.text()
		idx = int(str(item.text())[-1]) # this shouldn't really ever fail
		return self.gl_lights[idx]
	
	
	def refresh_light_states(self):
		self.quiet = True
		l = self.get_current_light()
		if l == None: return
		
		self.target().set_current_light(l)
		amb = glGetLightfv(l,GL_AMBIENT)
		dif =  glGetLightfv(l,GL_DIFFUSE)
		spec = glGetLightfv(l,GL_SPECULAR)
		self.light_ambient.r.setValue(float(amb[0]))
		self.light_ambient.g.setValue(float(amb[1]))
		self.light_ambient.b.setValue(float(amb[2]))
		
		self.light_diffuse.r.setValue(float(dif[0]))
		self.light_diffuse.g.setValue(float(dif[1]))
		self.light_diffuse.b.setValue(float(dif[2]))
		
		self.light_specular.r.setValue(float(spec[0]))
		self.light_specular.g.setValue(float(spec[1]))
		self.light_specular.b.setValue(float(spec[2]))
		
		pos = glGetLightfv(l,GL_POSITION)
		if pos[3] == 0:
			self.light_x_dir.setValue(float(pos[0]))
			self.light_y_dir.setValue(float(pos[1]))
			self.light_z_dir.setValue(float(pos[2]))
		else:
			# it's a point source
			
			# have to figure out the model coordinates, do this
			# by setting the light position to zero and then 
			# getting the light position. Calculate the difference...
			test_pos = [0,0,0,1]
			glLightfv(l,GL_POSITION,test_pos)
			test_pos_out = glGetLightfv(l,GL_POSITION)
			# reset to the correction position
			pos = [ (pos[i] - test_pos_out[i]) for i in range(3)]
			pos.append(1)
			glLightfv(l,GL_POSITION,pos)
			
			self.light_x_pos.setValue(float(pos[0]))
			self.light_y_pos.setValue(float(pos[1]))
			self.light_z_pos.setValue(float(pos[2]))
			
			dr = glGetLightfv(l,GL_SPOT_DIRECTION)
			self.light_ps_xdir.setValue(float(dr[0]))
			self.light_ps_ydir.setValue(float(dr[1]))
			self.light_ps_zdir.setValue(float(dr[2]))
			
			self.const_atten.setValue(float(glGetLightfv(l,GL_CONSTANT_ATTENUATION)))
			self.linear_atten.setValue(float(glGetLightfv(l,GL_LINEAR_ATTENUATION)))
			self.quad_atten.setValue(float(glGetLightfv(l,GL_QUADRATIC_ATTENUATION)))
			
			self.spot_cutoff.setValue(float(glGetLightfv(l,GL_SPOT_CUTOFF)))
			self.spot_exponent.setValue(float(glGetLightfv(l,GL_SPOT_EXPONENT)))
			
			
		self.target().updateGL()
		self.quiet = False
	def local_viewer_checked(self,i):
		glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER,i)
		self.target().updateGL()
	
	
	def del_directional_light(self):
		glDisable(self.get_current_light())
		self.redo_directional_light_list()
		
		self.refresh_light_states()
		self.target().updateGL()
		
	def del_pointsource_light(self):
		glDisable(self.get_current_light())
		self.redo_pointsource_light_list()
		self.refresh_light_states()
		self.target().updateGL()
		
	def redo_directional_light_list(self):
		self.light_list.clear()
		for i,l in enumerate(self.gl_lights):
			if glIsEnabled(l):
				pos = glGetLightfv(l,GL_POSITION)
				if pos[3] == 0:
					a = QtGui.QListWidgetItem("Light "+str(i),self.light_list)
					if len(self.light_list.selectedItems()) == 0:
						a.setSelected(True)
						
	def redo_pointsource_light_list(self):
		self.point_light_list.clear()
		for i,l in enumerate(self.gl_lights):
			if glIsEnabled(l):
				pos = glGetLightfv(l,GL_POSITION)
				if pos[3] == 1:
					a = QtGui.QListWidgetItem("Light "+str(i),self.light_list)
					if len(self.light_list.selectedItems()) == 0 and len(self.point_light_list.selectedItems()) == 0:
						a.setSelected(True)
		
	
	def new_directional_light(self):
		self.new_light()
	
	def new_pointsource_light(self):
		self.new_light(point_source=True)
		
	def new_light(self,point_source=False):
		for i,l in enumerate(self.gl_lights):
			if not glIsEnabled(l):
				glEnable(l)
				pos = glGetLightfv(l,GL_POSITION)
				
				if point_source:
					# make sure that it's a point source
					if pos[3] != 1:
						pos[3] = 1
						glLightfv(l,GL_POSITION,pos)
					
					spot_cutoff = glGetLightfv(l,GL_SPOT_CUTOFF)
					if spot_cutoff > 90: glLightfv(l,GL_SPOT_CUTOFF,90)
				else:
					# make sure that it's directionaly
					if pos[3] != 0:
						pos[3] = 0
						glLightfv(l,GL_POSITION,pos)
						
					
				new_label = "Light "+str(i)
				
				if not point_source: 
					a = QtGui.QListWidgetItem(new_label,self.light_list)
					for item in self.point_light_list.selectedItems(): item.setSelected(False)
						
				else:
					a = QtGui.QListWidgetItem(new_label,self.point_light_list)
					for item in self.light_list.selectedItems(): item.setSelected(False)
					
				a.setSelected(True)
				self.refresh_light_states()
				break
		else:
			print "write a message box to say that there are no available lights!"
	
	def new_point_source_light(self):
		for i,l in enumerate(self.gl_lights):
			if not glIsEnabled(l):
				glEnable(l)
				new_label = "Light "+str(i)
				
				a = QtGui.QListWidgetItem(new_label,self.light_list)
				a.setSelected(True)
				
				self.refresh_light_states()
				break
		else:
			print "write a message box to say that there are no available lights!"
			
		#print "new directional",glGetInteger(GL_MAX_LIGHTS)
		#for i in range(glGetInteger(GL_MAX_LIGHTS)):
			
			#if glIsEnabled(
	
	def light_list_clicked(self,item):
		for item in self.point_light_list.selectedItems(): item.setSelected(False)
		self.refresh_light_states()
		
	def point_light_list_clicked(self,item):
		for item in self.light_list.selectedItems(): item.setSelected(False)
		self.refresh_light_states()
		
	
	def get_directional_light_tab(self):
		
		self.directional_light_widget = QtGui.QWidget()
		
		vbl = QtGui.QVBoxLayout(self.directional_light_widget)
		vbl.setMargin(0)
		vbl.setSpacing(6)
		vbl.setObjectName("Lights")
		
		hbl = QtGui.QHBoxLayout()
		hbl.setMargin(0)
		hbl.setSpacing(6)
		hbl.setObjectName("hbl")
		vbl.addLayout(hbl)

		self.light_list = QtGui.QListWidget(None)
		self.light_list.setMouseTracking(True)
		self.redo_directional_light_list()

		hbl.addWidget(self.light_list)

		vbl2 = QtGui.QVBoxLayout()
		vbl2.setMargin(0)
		vbl2.setSpacing(6)
		vbl2.setObjectName("vbl2")
		hbl.addLayout(vbl2)
		
		new_light = QtGui.QPushButton("New")
		vbl2.addWidget(new_light)
		#copy_light = QtGui.QPushButton("Copy")
		#vbl2.addWidget(copy_light)
		del_light = QtGui.QPushButton("Delete")
		vbl2.addWidget(del_light)
		
		
		x_label = QtGui.QLabel()
		x_label.setText('x')
		
		self.light_x_dir = QtGui.QDoubleSpinBox(self)
		self.light_x_dir.setMinimum(-100000)
		self.light_x_dir.setMaximum(100000)
		self.light_x_dir.setValue(0.0)
	
		y_label = QtGui.QLabel()
		y_label.setText('y')
		
		self.light_y_dir = QtGui.QDoubleSpinBox(self)
		self.light_y_dir.setMinimum(-100000)
		self.light_y_dir.setMaximum(100000)
		self.light_y_dir.setValue(0.0)
		
		z_label = QtGui.QLabel()
		z_label.setText('z')
		
		self.light_z_dir = QtGui.QDoubleSpinBox(self)
		self.light_z_dir.setMinimum(-100000)
		self.light_z_dir.setMaximum(100000)
		self.light_z_dir.setValue(0.0)
		
		
		hbl_trans = QtGui.QHBoxLayout()
		hbl_trans.setMargin(0)
		hbl_trans.setSpacing(6)
		hbl_trans.setObjectName("Trans")
		hbl_trans.addWidget(x_label)
		hbl_trans.addWidget(self.light_x_dir)
		hbl_trans.addWidget(y_label)
		hbl_trans.addWidget(self.light_y_dir)
		hbl_trans.addWidget(z_label)
		hbl_trans.addWidget(self.light_z_dir)
		
		vbl.addLayout(hbl_trans)
		
		QtCore.QObject.connect(new_light, QtCore.SIGNAL("clicked()"), self.new_directional_light)
		QtCore.QObject.connect(del_light, QtCore.SIGNAL("clicked()"), self.del_directional_light)
		QtCore.QObject.connect(self.light_list, QtCore.SIGNAL("itemPressed(QListWidgetItem*)"), self.light_list_clicked)
		
		return self.directional_light_widget
	
	
	def get_pointsource_light_tab(self):
		
		self.pointsource_light_widget = QtGui.QWidget()
		
		vbl = QtGui.QVBoxLayout(self.pointsource_light_widget)
		vbl.setMargin(0)
		vbl.setSpacing(6)
		vbl.setObjectName("Lights")
		
		hbl = QtGui.QHBoxLayout()
		hbl.setMargin(0)
		hbl.setSpacing(6)
		hbl.setObjectName("hbl")
		vbl.addLayout(hbl)

		self.point_light_list = QtGui.QListWidget(None)
		self.point_light_list.setMouseTracking(True)
		
		self.redo_pointsource_light_list()
		#a = QtGui.QListWidgetItem(str("Light 0"),self.light_list)
		#a.setSelected(True)
		hbl.addWidget(self.point_light_list)

		vbl2 = QtGui.QVBoxLayout()
		vbl2.setMargin(0)
		vbl2.setSpacing(6)
		vbl2.setObjectName("vbl2")
		hbl.addLayout(vbl2)
		
		new_light = QtGui.QPushButton("New")
		vbl2.addWidget(new_light)
		#copy_light = QtGui.QPushButton("Copy")
		#vbl2.addWidget(copy_light)
		del_light = QtGui.QPushButton("Delete")
		vbl2.addWidget(del_light)
		
		
		pos_label = QtGui.QLabel()
		pos_label.setText('Pos: ')
		
		x_label = QtGui.QLabel()
		x_label.setText('x')
		
		self.light_x_pos = QtGui.QDoubleSpinBox(self)
		self.light_x_pos.setMinimum(-100000)
		self.light_x_pos.setMaximum(100000)
		self.light_x_pos.setValue(0.0)
	
		y_label = QtGui.QLabel()
		y_label.setText('y')
		
		self.light_y_pos = QtGui.QDoubleSpinBox(self)
		self.light_y_pos.setMinimum(-100000)
		self.light_y_pos.setMaximum(100000)
		self.light_y_pos.setValue(0.0)
		
		z_label = QtGui.QLabel()
		z_label.setText('z')
		
		self.light_z_pos = QtGui.QDoubleSpinBox(self)
		self.light_z_pos.setMinimum(-100000)
		self.light_z_pos.setMaximum(100000)
		self.light_z_pos.setValue(0.0)
		
		
		hbl_trans = QtGui.QHBoxLayout()
		hbl_trans.setMargin(0)
		hbl_trans.setSpacing(6)
		hbl_trans.setObjectName("Trans")
		hbl_trans.addWidget(pos_label)
		hbl_trans.addWidget(x_label)
		hbl_trans.addWidget(self.light_x_pos)
		hbl_trans.addWidget(y_label)
		hbl_trans.addWidget(self.light_y_pos)
		hbl_trans.addWidget(z_label)
		hbl_trans.addWidget(self.light_z_pos)
		
		vbl.addLayout(hbl_trans)
	
		
		self.const_atten = ValSlider(self.pointsource_light_widget,(0.0,5.0),"Const atten.:")
		self.const_atten.setValue(1.0)
		vbl.addWidget(self.const_atten)
		
		self.linear_atten = ValSlider(self.pointsource_light_widget,(0.0,0.5),"Linear atten.:")
		self.linear_atten.setValue(0.5)# why o why?
		self.linear_atten.setValue(0.0)
		vbl.addWidget(self.linear_atten)
		
		self.quad_atten = ValSlider(self.pointsource_light_widget,(0.0,0.05),"Quad. atten.:")
		self.quad_atten.setValue(0.05) # why o why?
		self.quad_atten.setValue(0.0)
		vbl.addWidget(self.quad_atten)
		
		
		dir_label = QtGui.QLabel()
		dir_label.setText('Dir: ')
		
		self.light_ps_xdir = QtGui.QDoubleSpinBox(self)
		self.light_ps_xdir.setMinimum(-100000)
		self.light_ps_xdir.setMaximum(100000)
		self.light_ps_xdir.setValue(0.0)
	
		y_label = QtGui.QLabel()
		y_label.setText('y')
		
		self.light_ps_ydir = QtGui.QDoubleSpinBox(self)
		self.light_ps_ydir.setMinimum(-100000)
		self.light_ps_ydir.setMaximum(100000)
		self.light_ps_ydir.setValue(0.0)
		
		z_label = QtGui.QLabel()
		z_label.setText('z')
		
		self.light_ps_zdir = QtGui.QDoubleSpinBox(self)
		self.light_ps_zdir.setMinimum(-100000)
		self.light_ps_zdir.setMaximum(100000)
		self.light_ps_zdir.setValue(0.0)
		
		hbl_trans2 = QtGui.QHBoxLayout()
		hbl_trans2.setMargin(0)
		hbl_trans2.setSpacing(6)
		hbl_trans2.setObjectName("Trans")
		hbl_trans2.addWidget(dir_label)
		hbl_trans2.addWidget(x_label)
		hbl_trans2.addWidget(self.light_ps_xdir)
		hbl_trans2.addWidget(y_label)
		hbl_trans2.addWidget(self.light_ps_ydir)
		hbl_trans2.addWidget(z_label)
		hbl_trans2.addWidget(self.light_ps_zdir)
		
		vbl.addLayout(hbl_trans2)
		
			
		self.spot_cutoff = ValSlider(self.pointsource_light_widget,(0.0,90.0),"Spot cutoff:")
		self.spot_cutoff.setValue(90)
		vbl.addWidget(self.spot_cutoff)
		
		self.spot_exponent = ValSlider(self.pointsource_light_widget,(0,10.0),"Spot exponent:")
		self.spot_exponent.setValue(1.0) # why o why?
		self.spot_exponent.setValue(0.0)
		vbl.addWidget(self.spot_exponent)
		
		
		QtCore.QObject.connect(new_light, QtCore.SIGNAL("clicked()"), self.new_pointsource_light)
		QtCore.QObject.connect(self.point_light_list, QtCore.SIGNAL("itemPressed(QListWidgetItem*)"), self.point_light_list_clicked)
		QtCore.QObject.connect(self.light_x_pos, QtCore.SIGNAL("valueChanged(double)"), self.update_light)
		QtCore.QObject.connect(self.light_y_pos, QtCore.SIGNAL("valueChanged(double)"), self.update_light)
		QtCore.QObject.connect(self.light_z_pos, QtCore.SIGNAL("valueChanged(double)"), self.update_light)
		QtCore.QObject.connect(self.light_ps_xdir, QtCore.SIGNAL("valueChanged(double)"), self.update_light)
		QtCore.QObject.connect(self.light_ps_ydir, QtCore.SIGNAL("valueChanged(double)"), self.update_light)
		QtCore.QObject.connect(self.light_ps_zdir, QtCore.SIGNAL("valueChanged(double)"), self.update_light)
		QtCore.QObject.connect(self.spot_cutoff, QtCore.SIGNAL("valueChanged"), self.update_light)
		QtCore.QObject.connect(self.spot_exponent, QtCore.SIGNAL("valueChanged"), self.update_light)
		QtCore.QObject.connect(self.const_atten, QtCore.SIGNAL("valueChanged"), self.update_light)
		QtCore.QObject.connect(self.linear_atten, QtCore.SIGNAL("valueChanged"), self.update_light)
		QtCore.QObject.connect(self.quad_atten, QtCore.SIGNAL("valueChanged"), self.update_light)

		QtCore.QObject.connect(del_light, QtCore.SIGNAL("clicked()"), self.del_pointsource_light)
		
		return self.pointsource_light_widget
	

class EMLightsInspector(QtGui.QWidget,EMLightsInspectorBase):
	def __init__(self,target) :
		QtGui.QWidget.__init__(self,None)
		EMLightsInspectorBase.__init__(self)
		self.target=weakref.ref(target)
		
		self.vbl = QtGui.QVBoxLayout(self)
		self.vbl.setMargin(0)
		self.vbl.setSpacing(6)
		self.vbl.setObjectName("vbl")
		
		self.hbl = QtGui.QHBoxLayout()
		self.hbl.setMargin(0)
		self.hbl.setSpacing(6)
		self.hbl.setObjectName("hbl")
		self.vbl.addLayout(self.hbl)

		self.vbl2 = QtGui.QVBoxLayout()
		self.vbl2.setMargin(0)
		self.vbl2.setSpacing(6)
		self.vbl2.setObjectName("vbl2")
		self.hbl.addLayout(self.vbl2)
		
		self.wiretog = QtGui.QPushButton("Wire")
		self.wiretog.setCheckable(1)
		self.vbl2.addWidget(self.wiretog)
		
		self.lighttog = QtGui.QPushButton("Light")
		self.lighttog.setCheckable(1)
		self.lighttog.setChecked(True)
		self.vbl2.addWidget(self.lighttog)
		
		self.tabwidget = QtGui.QTabWidget()
		self.maintab = None
		self.tabwidget.addTab(self.get_light_tab(), "Lights")
		self.tabwidget.addTab(self.get_main_tab(), "Transform")
		self.tabwidget.addTab(self.get_GL_tab(),"GL")
		self.vbl.addWidget(self.tabwidget)
		self.n3_showing = False
		self.quiet = False

		QtCore.QObject.connect(self.cbb, QtCore.SIGNAL("currentIndexChanged(QString)"), target.setColor)
		QtCore.QObject.connect(self.wiretog, QtCore.SIGNAL("toggled(bool)"), target.toggle_wire)
		QtCore.QObject.connect(self.lighttog, QtCore.SIGNAL("toggled(bool)"), target.toggle_light)
		QtCore.QObject.connect(self.glcontrast, QtCore.SIGNAL("valueChanged"), target.set_GL_contrast)
		QtCore.QObject.connect(self.glbrightness, QtCore.SIGNAL("valueChanged"), target.set_GL_brightness)

	def update_rotations(self,t3d):
		self.rotation_sliders.update_rotations(t3d)
	
	def set_scale(self,val):
		self.rotation_sliders.set_scale(val)
	
	def set_xy_trans(self, x, y):
		self.rotation_sliders.set_xy_trans(x,y)
	
	def set_xyz_trans(self,x,y,z):
		self.rotation_sliders.set_xyz_trans(x,y,z)
		#return self.advanced_tab
	
	def get_GL_tab(self):
		self.gltab = QtGui.QWidget()
		gltab = self.gltab
		
		gltab.vbl = QtGui.QVBoxLayout(self.gltab )
		gltab.vbl.setMargin(0)
		gltab.vbl.setSpacing(6)
		gltab.vbl.setObjectName("Main")
		
		self.hbl_color = QtGui.QHBoxLayout()
		self.hbl_color.setMargin(0)
		self.hbl_color.setSpacing(6)
		gltab.vbl.addLayout(self.hbl_color)

		self.color_label = QtGui.QLabel()
		self.color_label.setText('Material')
		self.hbl_color.addWidget(self.color_label)
		
		self.cbb = QtGui.QComboBox(gltab)
		self.hbl_color.addWidget(self.cbb)
		
		self.glcontrast = ValSlider(gltab,(1.0,5.0),"GLShd:")
		self.glcontrast.setObjectName("GLShade")
		self.glcontrast.setValue(1.0)
		gltab.vbl.addWidget(self.glcontrast)
		
		self.glbrightness = ValSlider(gltab,(-1.0,0.0),"GLBst:")
		self.glbrightness.setObjectName("GLBoost")
		self.glbrightness.setValue(0.1)
		self.glbrightness.setValue(0.0)
		gltab.vbl.addWidget(self.glbrightness)
	
		return gltab
	
	def get_main_tab(self):
		if ( self.maintab == None ):
			self.maintab = QtGui.QWidget()
			maintab = self.maintab
			maintab.vbl = QtGui.QVBoxLayout(self.maintab)
			maintab.vbl.setMargin(0)
			maintab.vbl.setSpacing(6)
			maintab.vbl.setObjectName("Main")
			
			self.rotation_sliders = EMTransformPanel(self.target(),self)
			self.rotation_sliders.addWidgets(maintab.vbl)
		
		return self.maintab
	
	def set_colors(self,colors,current_color):
		a = 0
		for i in colors:
			self.cbb.addItem(i)
			if ( i == current_color):
				self.cbb.setCurrentIndex(a)
			a += 1

		
# This is just for testing, of course
if __name__ == '__main__':
	from emapplication import EMApp
	from emglobjects import EM3DGLWidget
	em_app = EMApp()
	window = EM3DGLWidget()
	em_lights = EMLights(window)
	window.set_model(em_lights) 
	em_app.show()
	em_app.execute()
