#!/usr/bin/env python
# Muyuan Chen 2020-03
from EMAN2 import *
from EMAN2_utils import *
import numpy as np
import weakref
import OpenGL
OpenGL.ERROR_CHECKING = False
from OpenGL.GL import *
from OpenGL.GLU import *
from PyQt5 import QtGui, QtWidgets, QtCore
from PyQt5.QtCore import Qt
from eman2_gui.emapplication import get_application, EMApp
from eman2_gui.emimage2d import EMImage2DWidget
from eman2_gui.emimagemx import EMImageMXWidget
from eman2_gui.emshape import EMShape
import scipy.spatial.distance as scipydist
from scipy import ndimage

def main():

	usage=" "
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--gpuid", type=str,help="Specify the gpu to use", default="")
	#parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-2)
	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)
	options.setname=args[0]
	app = EMApp()

	drawer=EMPtclClassify(app,options)


	drawer.show()
	app.execute()
	E2end(logid)
	

class NNet:
	def __init__(self, boxsize=96):
		
		pad='valid'

		ki=tf.keras.initializers.TruncatedNormal(stddev=0.01)
		layers=[
			tf.keras.layers.Conv2D(8, 5, activation="relu", padding=pad,input_shape=(boxsize, boxsize,1)),
			tf.keras.layers.MaxPooling2D(),
			tf.keras.layers.Dropout(0.2),
			tf.keras.layers.Conv2D(16, 5, activation="relu", padding=pad),
			tf.keras.layers.MaxPooling2D(),
			tf.keras.layers.Dropout(0.2),
			tf.keras.layers.Conv2D(32, 5, activation="relu", padding=pad),
			tf.keras.layers.MaxPooling2D(),
			tf.keras.layers.Dropout(0.2),
			tf.keras.layers.Flatten(),
			tf.keras.layers.Dense(128, activation="relu"),
			tf.keras.layers.BatchNormalization(),
			tf.keras.layers.Dense(1, activation="sigmoid"),
		]
		self.model=model=tf.keras.Sequential(layers)
		self.outsz=self.model.layers[-1].output_shape#[1]
		
		self.boxsize=boxsize
		self.layers=layers
		
	def predict_class(self, inp, fromdata=True, usemax=False):
		if fromdata:
			inp=self.model(inp)
		inp=tf.math.minimum(inp,1.0)
		if usemax:
			out0=tf.math.reduce_max(inp*self.mask2, axis=(1,2))
			out1=1-tf.math.reduce_max(inp*self.mask, axis=(1,2))
		else:
			out0=tf.math.reduce_sum(inp*self.mask2, axis=(1,2))
			out1=tf.math.reduce_sum(abs(inp-self.mask)*self.mask2, axis=(1,2))
		return out0, out1
	
	def do_training(self, dataset, learnrate=1e-5, niter=10, posmult=0.5):
		#niter=1; learnrate=0
	
		lossbce=tf.keras.losses.BinaryCrossentropy()
		opt=tf.keras.optimizers.Adam(learning_rate=learnrate)
		self.model.compile(optimizer=opt, loss=lossbce)
		print("Training...")
		for it in range(niter):
			#self.model.reset_metrics()
			cost=[]
			for image, label in dataset:
				result = self.model.train_on_batch(image, label)
				cost.append(result)
			print("iteration {}, cost {:.3f}".format(it, np.mean(cost)))
		
	def apply_network(self, imgs):
		out=self.model.predict(imgs)
		return out
		
	def save_network(self, fname, options=None):
		if os.path.isfile(fname):
			os.remove(fname)
		print("Saving the nnet to {}...".format(fname))
		
		weights=[]
		bias=[]
		for ly in self.model.layers:
			if len(ly.weights)==2:
				w,b=ly.weights
				weights.append(w.numpy())
				bias.append(b.numpy())
		
		hdr=EMData(7,7)

		nlayer=len(self.model.layers)

		hdr.write_image(fname,0)

		k=1
		for i in range(len(weights)):
			w=weights[i]
			b=bias[i]
			w=w.transpose(3,2,0,1).copy()
			s=w.shape

			e=from_numpy(b)
			e["w_shape"]=s
			e.write_image(fname,k)
			k+=1
			w=w.reshape(s[0]*s[1], s[2], s[3])
			for wi in w:
				ws=wi.T.copy()
				e=from_numpy(ws)
				e.write_image(fname,k)
				k+=1	

	@staticmethod
	def load_network(fname, **kwargs):
		print("Loading nnet from {}...".format(fname))
		hdr=EMData(fname,0)
		
		nnet=NNet(**kwargs)
		k=1
		for ly in nnet.model.layers:
			
			wts=ly.weights
			if len(wts)==0:
				continue
			
			e=EMData(fname,k)
			s=e["w_shape"]
			b=e.numpy().copy()
			k+=1
			wts[1].assign(b, read_value=False)
			ks=wts[0].shape[1]
			allw=np.zeros((s[0]*s[1], ks, ks))
			for wi in range(s[0]*s[1]):
				e=EMData(fname,k)
				sw=e["nx"]
				#e=e.get_clip(Region((sw-ks)//2,(sw-ks)//2,ks,ks))
				k+=1
				w=e.numpy().copy()
				allw[wi]=w
			allw=allw.reshape([s[0], s[1], ks, ks]).transpose(3,2,1,0)
			wts[0].assign(allw, read_value=False)
			
		return nnet	

def get_image(ptcl, nimgs):
	ncopy=ceil(nimgs/len(ptcl))
	imgs=[]
	for p in ptcl:
		for i in range(ncopy):
			if i>0:
				pc=p.copy()
				pc.rotate(0,0,np.random.rand()*360)
			else:
				pc=p
			imgs.append(pc.numpy().copy())
	
	imgs=np.array(imgs, dtype=np.float32)
	if ncopy>1: 
		np.random.shuffle(imgs)
	imgs=imgs[:nimgs,:,:,None]
	
	return imgs


class TextBox(QtWidgets.QWidget):
	
	def __init__(self, label="", text=""):
		
		QtWidgets.QWidget.__init__(self)
		self.layout=QtWidgets.QHBoxLayout(self)
		self.layout.setContentsMargins(0,0,0,0)
		self.layout.addWidget(QtWidgets.QLabel(label))
		
		self.text=QtWidgets.QLineEdit()
		#self.text.setValidator()
		self.text.setText(str(text))
		self.layout.addWidget(self.text)
	
	def getval(self):
		return float(self.text.text())
	def setval(self, val):
		self.text.setText(str(val))
		
		
class EMPtclClassify(QtWidgets.QMainWindow):

	def __init__(self,application,options,datafile=None):
		QtWidgets.QWidget.__init__(self)
		self.setMinimumSize(150,100)
		self.setCentralWidget(QtWidgets.QWidget())
		self.gbl = QtWidgets.QGridLayout(self.centralWidget())

		#self.bt_new=QtWidgets.QPushButton("New")
		#self.bt_new.setToolTip("Build new neural network")
		#self.gbl.addWidget(self.bt_new, 0,0,1,2)
		
		#self.bt_load=QtWidgets.QPushButton("Load")
		#self.bt_load.setToolTip("Load neural network")
		#self.gbl.addWidget(self.bt_load, 1,0,1,2)
		
		self.bt_train=QtWidgets.QPushButton("Train")
		self.bt_train.setToolTip("Train neural network")
		self.gbl.addWidget(self.bt_train, 0,0,1,2)
		
		self.bt_save=QtWidgets.QPushButton("Save")
		self.bt_save.setToolTip("Save neural network")
		self.gbl.addWidget(self.bt_save, 1,0,1,2)
		
		#self.bt_apply=QtWidgets.QPushButton("Apply")
		#self.bt_apply.setToolTip("Apply neural network")
		#self.gbl.addWidget(self.bt_apply, 4,0,1,2)
				
		#self.bt_new.clicked[bool].connect(self.new_nnet)
		#self.bt_load.clicked[bool].connect(self.load_nnet)
		self.bt_train.clicked[bool].connect(self.train_nnet)
		self.bt_save.clicked[bool].connect(self.save_set)
		#self.bt_apply.clicked[bool].connect(self.apply_nnet)
		#self.bt_chgbx.clicked[bool].connect(self.change_boxsize)
		#self.box_display.currentIndexChanged.connect(self.do_update)

		self.val_learnrate=TextBox("LearnRate", 1e-4)
		self.gbl.addWidget(self.val_learnrate, 0,2,1,1)
		
		self.val_ptclthr=TextBox("PtclThresh", 200)
		self.gbl.addWidget(self.val_ptclthr, 1,2,1,1)
		
		self.val_niter=TextBox("Niter", 10)
		self.gbl.addWidget(self.val_niter, 2,2,1,1)
		
		self.options=options
		self.app=weakref.ref(application)
		self.nnet=None
		self.trainset=[]
		#self.nnetsize=96
		
		self.particles=EMData.read_images(options.setname)
		self.boxsz=self.particles[0]["nx"]
		#for p in self.particles:
			#p.process_inplace("math.fft.resample",{"n":float(boxsz/self.nnetsize)})
			#p.process_inplace("filter.lowpass.gauss", {"cutoff_abs":.4})
		
		self.ptclviewer=EMImageMXWidget()
		self.ptclviewer.setWindowTitle("Particles")
		self.sortidx=list(range(len(self.particles)))
		
		for boxview in [self.ptclviewer]:
			boxview.usetexture=False
			boxview.show()
			boxview.set_mouse_mode("Sets")
			boxview.rzonce=True
		
		
		self.ptclviewer.set_data(self.particles)
		self.ptclviewer.update()
		self.ptclimg=get_image(self.particles, len(self.particles))
		global tf
		tf=import_tensorflow(options.gpuid)
		
	def train_nnet(self):
		
		if self.nnet==None:
			self.nnet=NNet(self.boxsz)
			
		if 1:#len(self.trainset)==0:
			print("Preparing training set...")
			sets=self.ptclviewer.sets
			#print(sets)
			if (not "bad_particles" in sets) or  len(sets["bad_particles"])==0:
				print("No references.")
				return
				
			bids=sets["bad_particles"]
			ptcldata=self.ptclviewer.data
			#print(bids)
			badrefs=[ptcldata[i] for i in bids]
			
			gids=[i for i in range(len(ptcldata)) if i not in bids]
			np.random.shuffle(gids)
			nsample=512
			gids=gids[:nsample]
			#print(gids)
			goodrefs=[ptcldata[i] for i in gids]
			gimgs=get_image(goodrefs, nsample)
			bimgs=get_image(badrefs, nsample)
			#print(gimgs.shape, bimgs.shape)
			
			imgs=np.concatenate([gimgs, bimgs], axis=0)
			
			print(imgs.shape)
			labs=np.zeros(nsample*2, dtype=np.float32)
			labs[:nsample]=1.0
			
			self.trainset=(imgs, labs)
			
		dataset = tf.data.Dataset.from_tensor_slices(self.trainset)
		dataset=dataset.shuffle(500).batch(64)
		self.nnet.do_training(
			dataset, 
			learnrate=self.val_learnrate.getval(), 
			niter=int(self.val_niter.getval()),
			)
		
		
		bids=self.ptclviewer.sets["bad_particles"]
		#bids=[self.sortidx.index(i) for i in bids]
		bids=[self.sortidx[i] for i in bids]
		#print(bids)
		score=self.nnet.apply_network(self.ptclimg).flatten()
		#print(score)
		sid=np.argsort(score).tolist()
		self.sortidx=sid
		ptcls=[self.particles[i] for i in sid]
		self.ptclviewer.set_data(ptcls)
		idx=[sid.index(i) for i in bids]
		self.ptclviewer.enable_set("bad_particles",idx,update=True)
		self.ptclviewer.commit_sets()
		#print(idx)
		#print('------------------')
		#self.ptclviewer.sets={"bad_particles":s}
		self.ptclviewer.set_mouse_mode("Sets")
		self.ptclviewer.update()
		
	def save_set(self):
		fname=self.options.setname
		oname=fname[:fname.rfind('.')]+"_good.lst"
		thr=int(self.val_ptclthr.getval())
		#print(oname, thr)
		badi=self.sortidx[:thr]
		#print(badi)
		if os.path.isfile(oname):
			os.remove(oname)
		lst=LSXFile(fname, True)
		lout=LSXFile(oname, False)
		nn=lst.n
		for i in range(nn):
			if i in badi:
				continue
			l=lst.read(i)
			lout.write(-1, l[0], l[1], l[2])
			
		lst=lout=None
		print("{} particles written to {}".format(nn-thr, oname))
		
	def closeEvent(self, event):
		for b in  [self.ptclviewer]:
			b.close()
			
if __name__ == '__main__':
	main()
	
