#!/usr/bin/env python

#Author: Ross Coleman

from EMAN2 import test_image#,EMData
from emapplication import EMStandAloneApplication#, get_application
from emimage2d import EMImage2DModule
from emshape import EMShape
from math import *
import weakref

import sys
from PyQt4 import QtGui, QtCore

def counterGen():
	"""
	Calling this function will create a counter object.
	Ex: "counter = counterGen()"
	Then calling "counter.next()" will return the next number in {1, 2, 3, ...}
	"""
	i = 0
	while True:
		i += 1
		yield i

class Im2DModEventsHandler:
	"""Adds helix-boxing specific mouse interactivity to
	an EMImage2DModule object"""
	def __init__(self, target, boxwidth=50):
		"""The target should be an EMImage2DModule instance and the
		boxwidth is the width to draw the boxes."""
		self.counter = counterGen() #Used to make keys for the helix objects
		self.startPoint = None #This is the location of the first click when making a box
		self.endPoint = None #This is the location of the other endpoint for the box
		self.currentKey = None #This is the key for the "rectline" EMShape we're working on
		#The current key is created during a "mousedrag" if one doesn't exist and is set to None upon "mouseup"
		self.boxwidth = boxwidth
		self.editmode = None #Can be None, 'n' for new, 'm' for move, 'f' for move first point, or 's' for move second point
		self.old_box_key = None #Gives the key in EMImage2DModule's dictionary for the "rectline" EMShape to edit
		self.old_box_list = None #Gives the list used to create the old EMShape box
		self.color = (1, 1, 1)
		self.target = weakref.ref(target)
		QtCore.QObject.connect( self.target().emitter(), QtCore.SIGNAL("mousedown"), self.mouse_down)
		QtCore.QObject.connect( self.target().emitter(), QtCore.SIGNAL("mousedrag"), self.mouse_drag)
		QtCore.QObject.connect( self.target().emitter(), QtCore.SIGNAL("mouseup"), self.mouse_up)

	def get_helix_keys(self):
		"""This returns the keys for the "rectline" EMShape objects
		in a dictionary beloning to the EMImage2DModule object"""
		shapesdict = self.target().get_shapes()
		helix_keys = [shapekey for shapekey in shapesdict if shapesdict[shapekey].getShape()[0] == "rectline"]
		return helix_keys

	def get_nearest_helix_key(self, loc):
		"""
		loc is the coordinates for the mouse click.
		The function returns EMImage2DModule's key for the "rectline" EMShape
		object closest to the mouse click. Specifically, the nearest shape
		is determined by finding the closest endpoint or midpoint of a
		"rectline" EMShape object, and the key to that object is returned.
		"""
		helix_keys = self.get_helix_keys()
		shapesdict = self.target().get_shapes()
		pointsdict = {}
		for key in helix_keys:
			helix = shapesdict[key].getShape() #This is the list used to create the EMShape object
			pt0 = (helix[4], helix[5])
			pt1 = (helix[6], helix[7])
			midpoint = ( (pt0[0]+pt1[0])/2.0, (pt0[1]+pt1[1])/2.0 )
			if not pt0 in pointsdict:
			#If that same point belongs to another object, we won't replace it in pointsdict
				pointsdict[pt0] = key
			if not pt1 in pointsdict:
				pointsdict[pt1] = key
			if not midpoint in pointsdict:
				pointsdict[midpoint] = key

		closestpoint = None
		mindist_squared = None
		for point in pointsdict.keys():
			distance_squared = (point[0]-loc[0])**2 + (point[1]-loc[1])**2
			if (closestpoint == None) or (distance_squared < mindist_squared):
				mindist_squared = distance_squared
				closestpoint = point

		nearest_helix_key = pointsdict[closestpoint]
		return nearest_helix_key

	def mouse_down(self, event, loc):
		"""
		If the shift key is pressed and the click is inside a box, delete it.
		Otherwise, either create a new box or edit an existing one depending on click location.
		Imagine drawing two (infinite) lines through the long sides of each box.
		If the click is not between two of the lines for a box, we will create a new box.
		Then the behavior depends on distance from the shorter axis of symmetry--in other
		words, how far up or down the length of the box. Clicking in the middle 3/4 of the box
		(3/8 L from the shorter axis of symmetry) will result in moving the entire box.
		Clicking on a point betwen 3/8 L and 5/8 L from the shorter axis of symmetry
		results in moving that end of the box while keeping the midpoint of the other end fixed.
		"""
		self.startPoint = loc
		if self.get_helix_keys(): #helix boxes already exist
			closest_helix_key = self.get_nearest_helix_key(loc)

			box = self.target().get_shapes()[closest_helix_key].getShape() #list used to create the "rectline" EMShape
			mid = ( (box[4]+box[6])/2.0, (box[5]+box[7])/2.0 ) #midpoint
			length = sqrt( (box[6]-box[4])**2 + (box[7]-box[5])**2 )
			#Converting to a coordinate system with origin at the centroid and axes along the width and length
			#Unit vector along the length pointing from x0,y0 to x1,y1
			l_uvect = ( (box[6]-box[4])/length, (box[7]-box[5])/length )
			#Unit vector along the width such that w_uvect (cross product) l_uvect = k_hat
			w_uvect = ( l_uvect[1], -1*l_uvect[0] )
			click_vect = ( loc[0]-mid[0], loc[1] - mid[1] ) #vector from midpoint to click location
			#dot product of click_vect on w_uvect
			w_coord = click_vect[0]*w_uvect[0]+click_vect[1]*w_uvect[1]
			half_width = box[8]/2.0

			if abs(w_coord) > half_width: #outside box
				if event.modifiers()&QtCore.Qt.ShiftModifier: #nothing to delete
					self.editmode = None
					self.old_box_key = None
				else: #new box
					self.editmode = 'n' #new box
					self.old_box_key = None
			else:
				#dot product of click_vect on l_uvect
				l_coord = click_vect[0]*l_uvect[0]+click_vect[1]*l_uvect[1]
				if event.modifiers()&QtCore.Qt.ShiftModifier: #Shift key means delete
					self.editmode = None
					self.old_box_key = None
					if abs(l_coord) <= length /2.0: #inside box, so delete it
						self.target().del_shape(closest_helix_key)
						self.target().updateGL()
				else: #either new box or edit current box
					if l_coord > 5*length/8.0:
						self.editmode = 'n' #new box
						self.old_box_key = None
					elif l_coord >= 3*length/8.0:
						self.editmode = 's' #second point
						self.old_box_key = closest_helix_key
					elif l_coord > -3*length/8.0:
						self.editmode = 'm' #move entire box
						self.old_box_key = closest_helix_key
					elif l_coord >= -5*length/8.0:
						self.editmode = 'f' #first point
						self.old_box_key = closest_helix_key
					elif l_coord < -5*length/8.0:
						self.editmode = 'n' #new box
						self.old_box_key = None
					else:
						raise ValueError

		else: #no boxes exist
			if event.modifiers()&QtCore.Qt.ShiftModifier: #nothing to delete
				self.editmode = None
			else:
				self.editmode = 'n' #new box
			self.old_box_key = None #there is no old box

		if self.old_box_key:
			self.old_box_list = self.target().get_shapes()[self.old_box_key].getShape() #The list used to make the EMShape
		else:
			self.old_box_list = None

	def mouse_drag(self, event, location):
		"""
		Boxes are deleted in mouse_down, and the decision of how to edit is made there.
		However, new boxes are made and existing boxes are edited here.
		"""
		if self.startPoint and self.editmode: #self.startPoint and self.editmode are set in mouse_down
			if self.editmode == 'n': #new
				box = EMShape(["rectline", self.color[0], self.color[1], self.color[2],
					self.startPoint[0], self.startPoint[1],
					location[0], location[1], self.boxwidth])
				if not self.currentKey: 
					self.currentKey = self.counter.next() #We only make boxes if there is a mouse drag, so we create it here

				self.target().add_shape(self.currentKey, box)
			else:
				color = self.old_box_list[1:4]
				first = self.old_box_list[4:6]
				second = self.old_box_list[6:8]
				width = self.old_box_list[8]
				move = (location[0] - self.startPoint[0], location[1]-self.startPoint[1])

				if self.old_box_key in self.target().get_shapes().keys():
					self.target().del_shape(self.old_box_key)
				if self.editmode == 'm': #move entire box
					newbox = EMShape(["rectline", color[0], color[1], color[2],
						move[0]+first[0], move[1]+first[1],
						move[0]+second[0], move[1]+second[1], width])
				elif self.editmode == 'f': #move first point
					newbox = EMShape(["rectline", color[0], color[1], color[2],
						move[0]+first[0], move[1]+first[1], second[0], second[1], width])
				elif self.editmode == 's': #move second point
					newbox = EMShape(["rectline", color[0], color[1], color[2],
						first[0], first[1], move[0]+second[0], move[1]+second[1], width])

				self.target().add_shape(self.old_box_key, newbox) #replaces the previous box
				
			self.target().updateGL()

	def mouse_up(self, event, location):
		"""
		Once the mouse button comes back up, creating a new box, or editing
		an existing box is complete, so we need only clear variables relevant
		to creating or editing boxes.
		"""
		self.startPoint = None
		self.endPoint = None
		self.currentKey = None
		self.editmode = None
		self.old_box_key = None
		self.old_box_list = None

	def set_boxwidth(self, boxwidth):
		"""
		Set the width to draw the boxes in image coordinates, not screen coordinates.
		Thus, box width doesn't depend on your zoom level when drawing a box.
		"""
		self.boxwidth = boxwidth

def main():
	app = EMStandAloneApplication()
	module = EMImage2DModule(test_image(), app)
	eventsHandler = Im2DModEventsHandler(module)

	from OpenGL import GL
	print str(GL.glGetIntegerv(GL.GL_LINE_WIDTH_RANGE))
	#The line width setting seems to have no effect for "rectpoint" shapes
	rect1 = EMShape(["rectpoint", 1, 0, 0, 10, 10, 30, 30, 1])
	rect2 = EMShape(["rectpoint", 0, 0, 1, 10, 35, 30, 55, 10])
	module.add_shapes({"first": rect1, "second": rect2})

	module.show()
	sys.exit(app.exec_())

if __name__ == '__main__':
	main()
