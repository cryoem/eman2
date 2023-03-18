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
from eman2_gui.embrowser import EMSliceParamDialog,makeOrthoProj

def main():

	usage=" "
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--gpuid", type=str,help="Specify the gpu to use", default="")
	parser.add_argument("--keep", type=float,help="fraction to keep", default=.9)
	parser.add_argument("--nogui", action="store_true", default=False ,help="apply to a set of particles using trained network without the gui.")
	parser.add_argument("--batchsz", type=int,help="batch size for nogui mode", default=50000)

	#parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-2)
	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)
	options.setname=args[0]
	
	if options.nogui:
		print("loading network...")
		global tf
		tf=import_tensorflow(options.gpuid)
		ptcl=EMData(options.setname,0,True)
		boxsz=ptcl["nx"]		
		nnet=NNet(boxsz)
		nnet.model=tf.keras.models.load_model("nnet_classifycnn.h5",compile=False)
		bsz=options.batchsz
		nptcl=EMUtil.get_image_count(options.setname)
		nbatch=nptcl//bsz+1
		allscore=[]
		for ib in range(nbatch):
			last=min((ib+1)*bsz, nptcl-1)
			idx=list(range(ib*bsz,last))
			ptcls=EMData.read_images(options.setname, idx)
			ptclimg=get_image(ptcls, len(ptcls))
			
			score=nnet.apply_network(ptclimg).flatten()
			print("particle {} - {} : min {:.3f}, max {:.3f}".format(ib*bsz, last, np.min(score), np.max(score)))
			allscore.extend(score.tolist())
			
		print("sorting particles...")
		sortid=np.argsort(allscore).tolist()
		sortid=sortid[int(nptcl*(1-options.keep)):]
		lst=load_lst_params(options.setname)
		lstgood=[l for i,l in enumerate(lst) if i in sortid]
		oname=options.setname[:options.setname.rfind('.')]+"_good.lst"
		save_lst_params(lstgood, oname)
		print("output written to {}".format(oname))
		
		lstbad=[l for i,l in enumerate(lst) if i not in sortid]
		oname=options.setname[:options.setname.rfind('.')]+"_bad.lst"
		save_lst_params(lstbad, oname)
		
		
	
	else:
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
			tf.keras.layers.Conv2D(32, 5, activation="relu", padding=pad,input_shape=(boxsize, boxsize,1)),
			tf.keras.layers.MaxPooling2D(),
			tf.keras.layers.Dropout(0.2),
			tf.keras.layers.Conv2D(64, 5, activation="relu", padding=pad),
			tf.keras.layers.MaxPooling2D(),
			tf.keras.layers.Dropout(0.2),
			tf.keras.layers.Conv2D(128, 5, activation="relu", padding=pad),
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
		try: os.remove("nnet_classifycnn.h5")
		except:pass
		self.model.save("nnet_classifycnn.h5")
		
	def apply_network(self, imgs):
		out=self.model.predict(imgs)
		return out
		
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
		
		self.bt_train=QtWidgets.QPushButton("Train")
		self.bt_train.setToolTip("Train neural network")
		self.gbl.addWidget(self.bt_train, 0,0,1,2)
		
		self.bt_save=QtWidgets.QPushButton("Save")
		self.bt_save.setToolTip("Save particle set")
		self.gbl.addWidget(self.bt_save, 1,0,1,2)
		
		self.bt_load=QtWidgets.QPushButton("Load")
		self.bt_load.setToolTip("Load neural network")
		self.gbl.addWidget(self.bt_load, 2,0,1,2)
		
		
		#self.bt_apply=QtWidgets.QPushButton("Apply")
		#self.bt_apply.setToolTip("Apply neural network")
		#self.gbl.addWidget(self.bt_apply, 4,0,1,2)
				
		#self.bt_new.clicked[bool].connect(self.new_nnet)
		self.bt_load.clicked[bool].connect(self.load_nnet)
		self.bt_train.clicked[bool].connect(self.train_nnet)
		self.bt_save.clicked[bool].connect(self.save_set)
		#self.bt_apply.clicked[bool].connect(self.apply_nnet)
		#self.bt_chgbx.clicked[bool].connect(self.change_boxsize)
		#self.box_display.currentIndexChanged.connect(self.do_update)

		self.val_learnrate=TextBox("LearnRate", 1e-4)
		self.gbl.addWidget(self.val_learnrate, 0,2,1,1)
		
		nptcl=EMUtil.get_image_count(options.setname)
		self.val_ptclthr=TextBox("PtclThresh", int(nptcl*(1-options.keep)))
		self.gbl.addWidget(self.val_ptclthr, 1,2,1,1)
		
		self.val_niter=TextBox("Niter", 10)
		self.gbl.addWidget(self.val_niter, 2,2,1,1)
		
		self.options=options
		self.app=weakref.ref(application)
		self.nnet=None
		self.trainset=[]
		#self.nnetsize=96
		e=EMData(options.setname, 0, True)
		is3d=e["nz"]>1
		
		if is3d:
			nimg=EMUtil.get_image_count(options.setname)
			self.secparm=EMSliceParamDialog(self,nimg)
			ret=self.secparm.exec_()
			print(ret)
			
			layers=self.secparm.wspinlayers.value()
			center=self.secparm.wspincenter.value()
			lowpass=float(self.secparm.wlelp.text())
			highpass=float(self.secparm.wlehp.text())
			self.particles=[]
			for i in range(nimg):
				first=e["nx"]/2+center-layers
				ptcl=EMData(options.setname, i,False,Region(0,0,first, e["nx"], e["ny"], layers*2))
				x=ptcl.process("misc.directional_sum",{"axis":"z"})
				if lowpass>0 : x.process_inplace("filter.lowpass.gauss",{"cutoff_freq":1.0/lowpass})
				if highpass>0 : x.process_inplace("filter.highpass.gauss",{"cutoff_freq":1.0/highpass})
				x.process_inplace("normalize.edgemean")
				self.particles.append(x)
				sys.stdout.write("\r{}/{} finished.".format(i, nimg))
				sys.stdout.flush()
			print()
		else:
			self.particles=EMData.read_images(options.setname)
		self.boxsz=self.particles[0]["nx"]
		
		self.ptclviewer=EMImageMXWidget()
		self.ptclviewer.setWindowTitle("Particles")
		self.sortidx=list(range(len(self.particles)))
		
		for boxview in [self.ptclviewer]:
			boxview.usetexture=False
			boxview.show()
			boxview.set_mouse_mode("Sets")
			boxview.rzonce=True
		
		
		self.ptclviewer.set_data(self.particles)
		self.ptclviewer.sets["bad_particles"]=set()
		self.ptclviewer.update()
		self.ptclimg=get_image(self.particles, len(self.particles))
		global tf
		tf=import_tensorflow(options.gpuid)
		
	def train_nnet(self):
		
		if self.nnet==None:
			self.nnet=NNet(self.boxsz)
			
		if int(self.val_niter.getval())>0:
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
			thr=int(self.val_ptclthr.getval())
			gids=self.sortidx[thr:]
			gids=[i for i in gids if i not in bids]
			#gids=[i for i in range(len(ptcldata)) if i not in bids]
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
		oname2=fname[:fname.rfind('.')]+"_bad.lst"
		thr=int(self.val_ptclthr.getval())
		#print(oname, thr)
		badi=self.sortidx[:thr]
		#print(badi)
		if os.path.isfile(oname):
			os.remove(oname)
		lst=LSXFile(fname, True)
		lout=LSXFile(oname, False)
		lout2=LSXFile(oname2, False)
		nn=lst.n
		for i in range(nn):
			l=lst.read(i)
			if i in badi:
				lout2.write(-1, l[0], l[1], l[2])
			else:
				lout.write(-1, l[0], l[1], l[2])
			
		lst=lout=lout2=None
		print("{} particles written to {}".format(nn-thr, oname))
		
	def load_nnet(self):
		if self.nnet==None:
			self.nnet=NNet(self.boxsz)
			
		self.nnet.model=tf.keras.models.load_model("nnet_classifycnn.h5",compile=False)
		
	def closeEvent(self, event):
		for b in  [self.ptclviewer]:
			b.close()
			
if __name__ == '__main__':
	main()
	
