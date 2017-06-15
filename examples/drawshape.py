#!/usr/bin/env python
# Muyuan Chen 2017-03
from EMAN2 import *
import numpy as np
import weakref
from PyQt4 import QtCore, QtGui
from PyQt4.QtCore import Qt
from emapplication import get_application, EMApp
from emimage2d import EMImage2DWidget
from emshape import EMShape

def main():
	
	usage=" "
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	#parser.add_argument("--path", type=str,help="path", default="het_3d")
	parser.add_argument("--noupdate", action="store_true",help="do not erase shapes", default=False)
	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)
	img = EMData(args[0])

	app = EMApp()

	drawer=EMDrawWindow(app,options,datafile=img)


	drawer.show()
	app.execute()
	E2end(logid)
	
def run(cmd):
	print cmd
	launch_childprocess(cmd)
	
	
class EMDrawWindow(QtGui.QMainWindow):

	def __init__(self,application,options,datafile=None):
		QtGui.QWidget.__init__(self)
		self.imgview = EMImage2DWidget()
		self.setCentralWidget(QtGui.QWidget())
		self.gbl = QtGui.QGridLayout(self.centralWidget())
		
		self.gbl.addWidget(self.imgview,0,0)
		self.options=options
		self.app=weakref.ref(application)
		
		self.datafile=datafile
		self.imgview.set_data(datafile)
		
		self.shape=["hidden",1,0,0,0,0,2,2]
		self.state=0
		self.shape_index = 0
		self.imgview.shapes = {}

		print("x,y,major,minor,angle")
		
		QtCore.QObject.connect(self.imgview,QtCore.SIGNAL("mouseup"),self.mouseup  )
		QtCore.QObject.connect(self.imgview,QtCore.SIGNAL("mousemove"),self.mousemv)
		
	def update_view(self):
		if self.options.noupdate:
			self.imgview.shapes[self.shape_index] = EMShape(self.shape)
		else:
			self.imgview.shapes={0:EMShape(self.shape)}
		self.imgview.shapechange=1
		self.imgview.updateGL()
		
	def mouseup(self, event):
		x,y=self.imgview.scr_to_img((event.x(),event.y()))
		#x,y=int(x),int(y)
		if self.state==0:
			self.shape=["ellipse",1,0,0,x,y,3,3, 0,2]
			self.state=1
		elif self.state==1:
			self.state=2
		else:
			a = max(self.shape[6],self.shape[7])
			b = min(self.shape[6],self.shape[7])
			print("{},{},{},{},{}".format(self.shape[4],self.shape[5],a,b,self.shape[8]))
			self.shape_index += 1
			self.state=0

		
		self.update_view()
		
		
	def mousemv(self, event):
		x,y=self.imgview.scr_to_img((event.x(),event.y()))
		#x,y=int(x),int(y)
		if self.state==1:
			x0=self.shape[4]
			y0=self.shape[5]
			r=np.sqrt((x-x0)**2+(y-y0)**2)
			ang=np.arctan2((y-y0), (x-x0)) *180./np.pi
			
			self.shape=["ellipse",1,0,0,x0,y0,r,3,ang,2]
			
			self.update_view()
		elif self.state==2:
			x0=self.shape[4]
			y0=self.shape[5]
			r0=self.shape[6]
			ang=self.shape[8]
			r=np.sqrt((x-x0)**2+(y-y0)**2)
			
			self.shape=["ellipse",1,0,0,x0,y0,r0,r,ang,2]
			
			self.update_view()
	
if __name__ == '__main__':
	main()
	