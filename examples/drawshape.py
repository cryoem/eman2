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
		
		self.shape=[0]*6
		self.all_shapes=[]
		self.state=0
		self.shape_index = 0
		self.imgview.shapes = {}

		print("x,y,major,minor,angle")
		
		QtCore.QObject.connect(self.imgview,QtCore.SIGNAL("mouseup"),self.mouseup  )
		QtCore.QObject.connect(self.imgview,QtCore.SIGNAL("mousemove"),self.mousemv)
		
		
	def update_view(self):
		shps={}
		for i,s in enumerate(self.all_shapes):
			shp=["ellipse",1,0,0]+s[1:]+[2]+[s[0]]
			shps[i]=EMShape(shp)
			
		shps[len(self.all_shapes)]=EMShape(["ellipse",1,0,0]+self.shape[1:]+[2]+[self.shape[0]])
		self.imgview.shapes=shps
		self.imgview.shapechange=1
		self.imgview.updateGL()
		
	def mouseup(self, event):
		x,y=self.imgview.scr_to_img((event.x(),event.y()))
		#x,y=int(x),int(y)
		if self.state==0:
			self.shape=[self.imgview.list_idx, x,y,3,3,0]
			self.state=1
		elif self.state==1:
			self.state=2
		else:
			a = max(self.shape[3],self.shape[4])
			b = min(self.shape[3],self.shape[4])
			print("{},{},{},{},{}".format(self.shape[1],self.shape[2],a,b,self.shape[5]))
			#self.shape_index += 1
			if self.options.noupdate:
				self.all_shapes.append(self.shape)
			self.state=0
			#print self.shape

		
		self.update_view()
		
		
	def mousemv(self, event):
		x,y=self.imgview.scr_to_img((event.x(),event.y()))
		#x,y=int(x),int(y)
		if self.state==1:
			x0=self.shape[1]
			y0=self.shape[2]
			r=np.sqrt((x-x0)**2+(y-y0)**2)
			ang=np.arctan2((y-y0), (x-x0)) *180./np.pi
			
			self.shape=[self.imgview.list_idx, x0,y0,r,3,ang]
			
			self.update_view()
		elif self.state==2:
			x0=self.shape[1]
			y0=self.shape[2]
			r0=self.shape[3]
			ang=self.shape[5]
			r=np.sqrt((x-x0)**2+(y-y0)**2)
			
			self.shape=[self.imgview.list_idx, x0,y0,r0,r,ang]
			
			self.update_view()
	
if __name__ == '__main__':
	main()
	