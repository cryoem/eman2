#!/usr/bin/env python
from __future__ import print_function
from __future__ import division
# Muyuan Chen 2017-03
from past.utils import old_div
from EMAN2 import *
import numpy as np
import weakref
from PyQt5 import QtCore, QtGui, QtWidgets
from eman2_gui.emapplication import get_application, EMApp
from eman2_gui.emimage2d import EMImage2DWidget
from eman2_gui.emshape import EMShape

def main():
	
	usage=" "
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	#parser.add_argument("--path", type=str,help="path", default="het_3d")
	parser.add_argument("--noupdate", action="store_true",help="do not erase shapes", default=False)
	parser.add_argument("--fromedge", action="store_true",help="draw from edge", default=False)
	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)
	img = EMData(args[0])

	app = EMApp()

	drawer=EMDrawWindow(app,options,datafile=img)


	drawer.show()
	drawer.raise_()
	app.execute()
	E2end(logid)
	
def run(cmd):
	print(cmd)
	launch_childprocess(cmd)
	
	
class EMDrawWindow(QtWidgets.QMainWindow):

	def __init__(self,application,options,datafile=None):
		QtWidgets.QWidget.__init__(self)
		self.imgview = EMImage2DWidget()
		self.setCentralWidget(QtWidgets.QWidget())
		self.gbl = QtWidgets.QGridLayout(self.centralWidget())
		
		self.gbl.addWidget(self.imgview,0,0)
		self.options=options
		self.app=weakref.ref(application)
		
		self.datafile=datafile
		self.imgview.set_data(datafile)
		
		self.shape=[0]*6
		self.all_shapes=[]
		self.state=0
		self.shape_index = 0
		self.origin=[0,0]
		self.imgview.shapes = {}

		print("imgnum,x,y,major,minor,angle")
		
		self.imgview.mouseup.connect(self.on_mouseup)
		self.imgview.mousemove.connect(self.mousemv)
		
		
	def update_view(self):
		shps={}
		for i,s in enumerate(self.all_shapes):
			shp=["ellipse",1,0,0]+s[1:]+[2]+[s[0]]
			shps[i]=EMShape(shp)
		
		shps[len(self.all_shapes)]=EMShape(["ellipse",1,0,0]+self.shape[1:]+[2]+[self.shape[0]])
		self.imgview.shapes=shps
		self.imgview.shapechange=1
		self.imgview.updateGL()
		
	def on_mouseup(self, event):
		x,y=self.imgview.scr_to_img((event.x(),event.y()))
		#x,y=int(x),int(y)
		if self.state==0:
			self.shape=[self.imgview.list_idx, x,y,3,3,0]
			self.origin=[x,y]
			self.state=1
		elif self.state==1:
			self.state=2
		else:
			a = max(self.shape[3],self.shape[4])
			b = min(self.shape[3],self.shape[4])
			print(("{},{},{},{},{},{}".format(self.imgview.list_idx,self.shape[1],self.shape[2],a,b,self.shape[5])))
			#self.shape_index += 1
			if self.options.noupdate:
				self.all_shapes.append(self.shape)
			self.state=0

		self.update_view()
		
		
	def mousemv(self, event):
		x,y=self.imgview.scr_to_img((event.x(),event.y()))
		#x,y=int(x),int(y)
		if self.state==1:
			if self.options.fromedge:
				dx=x-self.origin[0]
				dy=y-self.origin[1]
				r=old_div(np.sqrt((dx)**2+(dy)**2),2.)
				x0=x-old_div(dx,2.)
				y0=y-old_div(dy,2.)
			else:
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
	
