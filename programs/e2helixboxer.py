#!/usr/bin/env python
from EMAN2 import test_image#,EMData
from emapplication import EMStandAloneApplication#, get_application
from emimage2d import EMImage2DModule
from emshape import EMShape
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
		print helix_keys
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
#			print "start %s, end %s, mid %s" % (pt0, pt1, midpoint)
			if not pt0 in pointsdict:
			#If that same point belongs to another object, we won't replace it in pointsdict
				pointsdict[pt0] = key
			if not pt1 in pointsdict:
				pointsdict[pt1] = key
			if not midpoint in pointsdict:
				pointsdict[midpoint] = key
		print pointsdict

		closestpoint = None
		mindist_squared = None
		for point in pointsdict.keys():
			distance_squared = (point[0]-loc[0])**2 + (point[1]-loc[1])**2
			if (closestpoint == None) or (distance_squared < mindist_squared):
				mindist_squared = distance_squared
				closestpoint = point
		
		print "Click location: %s, nearest point: %s" % (loc, closestpoint)
		nearest_helix_key = pointsdict[closestpoint]
		return nearest_helix_key

	def mouse_down(self, event, location):
		#print "You clicked on image coordinates: " + str(location)
		if self.get_helix_keys(): #helix boxes already exist
			closest_helix_key = self.get_nearest_helix_key(location)
			print closest_helix_key
		if not self.startPoint:
			self.startPoint = location

	def mouse_drag(self, event, location):
		if self.startPoint:
			box = EMShape(["rectline", self.color[0], self.color[1], self.color[2],
					self.startPoint[0], self.startPoint[1],
					location[0], location[1], self.boxwidth])
			if not self.currentKey: 
				self.currentKey = self.counter.next() #We only make boxes if there is a mouse drag, so we create it here
			self.target().add_shape(self.currentKey, box)
			self.target().updateGL()

	def mouse_up(self, event, location):
		if self.startPoint:
			self.endPoint = location
			if self.currentKey: #This is true only if a "mousedrag" event occurred
				box = EMShape(["rectline", self.color[0], self.color[1], self.color[2],
					self.startPoint[0], self.startPoint[1],
					self.endPoint[0], self.endPoint[1], self.boxwidth])
				self.target().add_shape(self.currentKey, box)
				self.target().updateGL()
		#Once the mouse comes up, we finish creating a box
		self.startPoint = None
		self.endPoint = None
		self.currentKey = None

	def set_boxwidth(self, boxwidth):
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
