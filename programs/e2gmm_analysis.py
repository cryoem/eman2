#!/usr/bin/env python
# Muyuan Chen 2020-05
from EMAN2 import *
import numpy as np
from PyQt5 import QtGui, QtWidgets, QtCore
from PyQt5.QtCore import Qt
from eman2_gui.emapplication import get_application, EMApp
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure
from sklearn.decomposition import PCA

os.environ["CUDA_VISIBLE_DEVICES"]='0' 
os.environ["TF_FORCE_GPU_ALLOW_GROWTH"]='true' 
import tensorflow as tf

def main():
	
	usage=" "
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--midinput", type=str,help="middle layer input", default=None)
	parser.add_argument("--decoder", type=str,help="decoder input", default=None)

	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)

	app = EMApp()
	win=EMSptEval(app, options)
	win.show()
	app.execute()
	
	E2end(logid)
	


class MplCanvas(FigureCanvasQTAgg):

    def __init__(self, parent=None):
        fig = Figure(figsize=(4, 4))
        self.axes = fig.add_subplot(111)
        super(MplCanvas, self).__init__(fig)
        
        a=np.random.randn(1000,2)
        self.axes.plot(a[:,0], a[:,1],'.')


class EMSptEval(QtWidgets.QMainWindow):

	def __init__(self,application,options):
		QtWidgets.QWidget.__init__(self)
		
		self.options=options
		
		self.setMinimumSize(1000,200)
		self.setCentralWidget(QtWidgets.QWidget())
		self.gbl = QtWidgets.QGridLayout(self.centralWidget())
		
		
		#self.bt_save=QtWidgets.QPushButton("Save")
		#self.gbl.addWidget(self.bt_save, 1,1,1,1)
		##self.bt_save.clicked[bool].connect(self.save_json)

		
		self.plotwinx=MplCanvas(self)
		self.plotwinx.setMinimumSize(500, 500)
		self.plotx=self.plotwinx.axes
		self.gbl.addWidget(self.plotwinx, 1,2,10,1)
		self.plotwinx.mpl_connect('button_press_event', self.onclick_plotx)
	
		self.plotwiny=MplCanvas(self)
		self.plotwiny.setMinimumSize(500, 500)
		self.ploty=self.plotwiny.axes
		self.gbl.addWidget(self.plotwiny, 1,3,10,1)
		
		self.mid=np.loadtxt(options.midinput)
		print(self.mid.shape)
		self.decoder=tf.keras.models.load_model(options.decoder)
		self.decoder.summary()

		
		self.pca=PCA(2)
		p2=self.pca.fit_transform(self.mid[:,1:])
		p2-=np.mean(p2, axis=0)
		p2/=np.std(p2, axis=0)
		self.pcaout=p2.copy()
		self.plotx.cla()
		self.plotx.scatter(p2[:,0], p2[:,1],s=5, alpha=.2)
		self.plotwinx.draw()
		
		self.traj=[]


	def onclick_plotx(self,event):
		print(event.xdata, event.ydata)
		pos=[event.xdata, event.ydata]
		if len(self.traj)==1:
			self.traj.append(pos)
			self.draw_traj()
			
		else:
			self.traj=[pos]
			
		print(self.traj)
		
	def draw_traj(self):
		print(self.traj)
		cf=self.pca.inverse_transform(self.traj).astype(np.float32)
		print(cf)
		p=self.decoder(cf*0)[0]-.5
		print(p.shape)
		pc=self.decoder(cf)
		pc=pc[0]-pc[1]
		
		plt=self.ploty
		plt.cla()
		plt.quiver(p[::2,0], -p[::2,1], pc[::2,0], -pc[::2,1],scale=.6, width=3e-3)
		
		plt.axis("square")
		#plt.xlim(-.4, .4)
		#plt.ylim(-.4, .4)
		self.plotwiny.draw()
	
if __name__ == '__main__':
	main()
	