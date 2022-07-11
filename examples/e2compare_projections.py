#!/usr/bin/env python
import sys
from PyQt5 import QtGui, QtWidgets, QtCore
from EMAN2 import *
from eman2_gui.emapplication import get_application, EMApp
from eman2_gui.emimage import EMImageWidget
from eman2_gui.emimage2d import EMImage2DWidget
from eman2_gui.emimagemx import EMImageMXWidget
from eman2_gui.emscene3d import EMScene3D
from eman2_gui.emdataitem3d import EMDataItem3D, EMIsosurface

from eman2_gui.valslider import ValSlider
import weakref


def main():
	usage="""display projections of input 3D files"

    [prog] <input>

    """
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_pos_argument(name="input", help="Specify an image name to display")
	(options, args) = parser.parse_args()
	#logid=E2init(sys.argv)

	if len(args) == 0:
		print("INPUT ERROR: You must specify an image to display.")
		sys.exit(1)


	app = EMApp()
	cp = ControlPanel(app,args)
	cp.show()


	x=app.exec_()
	#E2end(logid)
	sys.exit(0)

class ControlPanel(QtWidgets.QMainWindow):
	def __init__(self, application,fname_l=[]):
		super().__init__()
		self.app = weakref.ref(application)
		self.setWindowTitle("Control Panel")
		self.fname_l=fname_l
		self.centralWidget=QtWidgets.QWidget()
		self.setCentralWidget(self.centralWidget) #set central widget
		self.datalist = [EMData(f) for f in self.fname_l]
		self.d_sz=[self.datalist[0]['nx'],self.datalist[0]['ny'],self.datalist[0]['nz']] #list of data size x,y,z

		#Create image widgets
		self.im2d=[EMImage2DWidget() for f in self.fname_l]
		self.im3d=EMScene3D()
		self.im3d.setWindowTitle(self.fname_l[0])
		self.im3ddata=EMDataItem3D(self.datalist[0])
		self.isosurface = EMIsosurface(self.im3ddata)
		self.im3d.insertNewNode("Item", self.im3ddata, parentnode=self.im3d)
		self.im3d.insertNewNode("Iso", self.isosurface, parentnode=self.im3ddata)
		self.im3d.initialViewportDims(self.im3ddata.getData().get_xsize())

		#Create ValSliders for Euler Angles and Value for Filters
		self.vsalt=ValSlider(self,(0,360),label='alt',labelwidth=85,rounding=1)
		self.vsaz=ValSlider(self,(0,360),label='az',labelwidth=85,rounding=1)
		self.vsphi=ValSlider(self,(0,360),label='phi',labelwidth=85,rounding=1)
		self.vslow=ValSlider(self,(0,50),label='Low Pass (A)',labelwidth=85,rounding=1)
		self.vshigh=ValSlider(self,(0,50),label='High Pass (A)',labelwidth=85,rounding=1)

		#Create SumLayers Spinbox Widget
		self.sumLayers = QtWidgets.QWidget()
		sl_layout = QtWidgets.QHBoxLayout()
		sl_layout.addWidget(QtWidgets.QLabel("Sum Layers (0->1, 1->3, 2->5,...)"))
		self.sumLayersBox=QtWidgets.QSpinBox(self.centralWidget)
		self.sumLayersBox.setMinimum(-1)
		self.sumLayersBox.setValue(-1)
		sl_layout.addWidget(self.sumLayersBox)
		self.sumLayers.setLayout(sl_layout)

		#Add sliders and spinbox to the central widget
		self.vbl=QtWidgets.QVBoxLayout(self.centralWidget)
		self.vbl.addWidget(self.vsalt)
		self.vbl.addWidget(self.vsaz)
		self.vbl.addWidget(self.vsphi)
		self.vbl.addWidget(self.vslow)
		self.vbl.addWidget(self.vshigh)
		self.vbl.addWidget(self.sumLayers)

		#Connect VS and Spinbox with appropriate event functions
		self.vsalt.valueChanged.connect(self.neworientation_fromslider)
		self.vsaz.valueChanged.connect(self.neworientation_fromslider)
		self.vsphi.valueChanged.connect(self.neworientation_fromslider)
		self.vslow.valueChanged.connect(self.neworientation_fromslider)
		self.vshigh.valueChanged.connect(self.neworientation_fromslider)
		self.sumLayersBox.valueChanged.connect(self.newdata_fromspinbox)
		self.im3d.sgtransform.connect(self.neworientation_from3d)

		#Show all image widgets
		self.im3d.show()
		for i in range(len(self.im2d)):
			self.im2d[i].setWindowTitle(self.fname_l[i])
			self.im2d[i].show()
		#Initiate the Image2DWidgets with projections of the non-transfrom volumes
		self.newOrientation_2d(Transform())

	def neworientation_fromslider(self,event):
		"""Listen to change from ValSlider and update image widgets"""
		self.xform = Transform({"type":"eman","alt":self.vsalt.value,"az":self.vsaz.value,"phi":self.vsphi.value})
		self.newOrientation_3d(xform=self.xform)
		self.newOrientation_2d(xform=self.xform,lp=self.vslow.value,hp=self.vshigh.value)

	def newdata_fromspinbox(self,event):
		"""Listen to change from SumLayers Spinbox, update data_clip list and 2D image widgets accordingly"""
		self.newOrientation_2d(xform=self.xform,lp=self.vslow.value,hp=self.vshigh.value)
		
			#first = self.d_sz[2]//2-n_layers
			#last = self.d_sz[2]//2+n_layers+1
			#self.projs = [data.process("misc.directional_sum",{"axis":"z","first":first,"last":last}) for data in self.datalist]
		#print(self.data_clip[0]['nz'])

	def neworientation_from3d(self,x,y):
		self.newOrientation_2d(y)
		ort=y.get_rotation()
		self.vsalt.setValue(ort["alt"],1)
		self.vsaz.setValue(ort["az"],1)
		self.vsphi.setValue(ort["phi"],1)

	def newOrientation_3d(self,xform):
		"""Update EMScene3D"""
		self.im3d.transform=xform
		self.im3d.updateSG()

	def newOrientation_2d(self,xform,lp=0,hp=0):
		"""Update EMImage2DWidgets"""
		layers=int(self.sumLayersBox.value())
		#nx=self.datalist[0]["nx"]
		#ny=self.datalist[0]["ny"]
		nz=self.datalist[0]["nz"]
		if layers<0: layers=nz//2
		self.projs = [data.process("xform",{"transform":xform}).process("misc.directional_sum",{"first":nz//2-layers,"last":nz//2+layers,"axis":"z"}) for data in self.datalist]
		
		if lp > 0:
			for im in self.projs: im.process_inplace("filter.lowpass.gauss",{"cutoff_freq":1.0/lp})
		if hp > 0:
			for im in self.projs: im.process_inplace("filter.highpass.gauss",{"cutoff_freq":1.0/hp})

		for i in range(len(self.im2d)):
			self.im2d[i].set_data(self.projs[i])

	def closeEvent(self,event):
		"""Close everything when close the Control Panel"""
		#print("Bye!!!")
		self.im3d.close()
		for im in self.im2d:
			im.close()

if __name__ == '__main__':
	main()

