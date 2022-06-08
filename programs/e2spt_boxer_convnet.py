#!/usr/bin/env python
# Muyuan Chen 2020-01
from builtins import range
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
from scipy.spatial import KDTree

def main():

	usage=" "
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--label", type=str,help="Load previous contour segmentation.", default="tomobox",guitype='strbox',row=0, col=0, rowspan=1, colspan=1)
	parser.add_argument("--gpuid", type=str,help="Specify the gpu to use", default="0",guitype='strbox',row=1, col=0, rowspan=1, colspan=1)
	parser.add_argument("--mult", type=float,help="multiply data by factor. useful for vpp data...", default=1)
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-2)
	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)
	
	app = EMApp()

	test = EMImage2DWidget()
	drawer=EMTomobox(app,options)


	drawer.show()
	app.execute()
	E2end(logid)
	

def get_circle(p,r):
	t=np.arange(0, np.pi*2, np.pi/10.)
	pts=r*np.vstack([np.cos(t), np.sin(t)]).T
	pts+=np.array(p[:2])
	return pts

class BoxShapes(EMShape):
	def __init__(self, img, points):
		self.points=points
		self.isanimated = False
		self.shape=["tomo_boxes",0]
		self.image=img
		self.triangles=[]
		self.select=0
		self.classid=0
		self.circlesize=32
		
	def add_point(self, newpt=[], cls=0):
		zpos=self.image.list_idx
		#ci=self.select
		self.points.append([newpt[0], newpt[1], zpos, int(cls>0)])
		#print(self.points)
	

	def draw(self,d2s=None,col=None):
		zpos=self.image.list_idx
		for p in self.points:
			dz=abs(p[2]-zpos)
			if dz>=self.circlesize:
				continue
			
			lns=get_circle(p, self.circlesize-dz)
			if p[3]==1:
				glColor3f( .2, .2, 1 );
			else:
				glColor3f( 1, .2, .2 );
				
			glLineWidth(3.)
			glEnableClientState(GL_VERTEX_ARRAY)
			glVertexPointerf(lns)
			glDrawArrays(GL_LINES, 0, len(lns))

		return
	
class NNet:
	def __init__(self, boxsize=96, thick=9):
		
		pad='same'

		ki=tf.keras.initializers.TruncatedNormal(stddev=0.01)
		layers=[
			tf.keras.layers.Conv2D(48, 9, activation="relu", padding=pad, kernel_initializer=ki,input_shape=(boxsize, boxsize, thick//3)),
			tf.keras.layers.MaxPooling2D(),
			tf.keras.layers.Conv2D(48, 9, activation="relu", padding=pad, kernel_initializer=ki),
			tf.keras.layers.MaxPooling2D(),
			tf.keras.layers.Conv2D(1, 9, padding=pad, activation="relu", kernel_initializer=ki),
		]
		self.model=model=tf.keras.Sequential(layers)
		self.outsz=self.model.layers[-1].output_shape[1]
		
		self.thick=thick
		self.boxsize=boxsize
		self.layers=layers
		
	#def predict_class(self, inp, usemax=False):
		#inp=tf.math.minimum(inp,1.0)
		#if usemax:
			#out0=tf.math.reduce_max(inp*self.mask2, axis=(1,2))
			#out1=1-tf.math.reduce_sum(inp*self.mask, axis=(1,2))
		#else:
			#out0=tf.reduce_sum(inp*self.mask2, axis=(1,2))
			#out1=tf.reduce_sum(abs(inp-self.mask)*self.mask2, axis=(1,2))
		#return out0, out1
	
	def do_training(self, dataset, learnrate=1e-5, niter=10, tarsz=2, usemax=False, posmult=0.5):
		#niter=1; learnrate=0
		#### prepare the output mask
		print(usemax)
		ksize=(1./tarsz)**2
		sz=self.outsz
		m=np.exp(-ksize*np.sum((np.indices((sz,sz))-(sz-1)/2.)**2, axis=0))
		m/=np.max(m)
		self.mask=m[None,:,:,None]
		m2=np.exp(-.25*ksize*np.sum((np.indices((sz,sz))-(sz-1)/2.)**2, axis=0))
		m2/=np.max(m2)
		self.mask2=m2[None,:,:,None]
		posmult=posmult/(1+posmult)
		
		opt=tf.keras.optimizers.Adam(learning_rate=1e-5) 
		wts=self.model.trainable_variables
		print("Training...")
		for it in range(niter):
			#self.model.reset_metrics()
			cost=[]
			for image, label in dataset:
				with tf.GradientTape() as gt:
					py=self.model(image, training=True)
					if usemax:
						py1=tf.math.minimum(py,1.0)
						out0=1-tf.math.reduce_max(py1*self.mask2, axis=(1,2))[:,0]
						out1=tf.math.reduce_sum((py*self.mask)**2, axis=(1,2))[:,0]
						#out0=tf.math.maximum(out0,.0)
					else:
						out0=tf.reduce_sum((py-self.mask2)**2, axis=(1,2))[:,0]
						out1=tf.reduce_sum(py**2, axis=(1,2))[:,0]
					
					
					loss=tf.reduce_mean(out0*label*posmult+out1*(1-label)*(1-posmult))
					
				grad=gt.gradient(loss, wts)
				#print(loss, label, out0[:,0], out1[:,0])
				opt.apply_gradients(zip(grad, wts))
					
				#result = self.model.train_on_batch(image, label)
				cost.append(loss)
			print("iteration {}, cost {:.3f}".format(it, np.mean(cost)))
			#print(out0, out1, label, )
		
	def apply_network(self, img):
		bigbox=img.shape[1]
		ly=tf.keras.layers.Conv2D(48, 9, activation="relu", input_shape=[bigbox,bigbox,self.thick//3])
		modelbig=tf.keras.Sequential([ly]+self.layers[1:])
		ly.weights[0].assign(self.layers[0].weights[0], read_value=False)
		ly.weights[1].assign(self.layers[0].weights[1], read_value=False)
		bsz=int(2e7/(img.shape[1]*img.shape[2]))
		bsz=min(bsz,64)
		nb=len(img)//bsz+1
		out=[]
		print("image shape: {}, batch size: {}".format(img.shape, bsz))
		for i in range(nb):
			m=img[i*bsz:(i+1)*bsz]
			if m.shape[0]==0: continue
			#print(m.shape)
			o=modelbig.predict(m)
			out.append(o)
			
		out=np.concatenate(out,axis=0)
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
		
		#shp=hdr["imageshape"]
		#if imgsz>0: sz=imgsz
		#else:  sz=shp[-1]
		
		#if bsz==0:  batchsize=shp[0]
		#else: batchsize=bsz

		#kernels=[(k,ksize[i],poolsz[i]) for i,k in enumerate(nkernel)]
		
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

def get_image(img, pos=[], boxsz=-1, thick=9, ncopy=1):
	
	if boxsz<0:
		bx=img["nx"]
		by=img["ny"]
		b=max(bx,by)
		
	else:
		b=bx=by=boxsz
		
	t=thick
	
	if len(pos)==0:
		pos=[img["nx"]//2, img["ny"]//2, img["nz"]//2]
	
	m=img.get_clip(Region(pos[0]-b//2, pos[1]-b//2, pos[2]-t//2, b, b, t))
	
	if ncopy>1:
		m.process_inplace("mask.soft",{"outer_radius":b//2-6,"width":4})
		rots=np.random.rand(ncopy)*360
	else:
		rots=[0]
		
	ptcls=[]
	for i in range(t//3):
		ptcls.append(m.process("misc.directional_sum", {"axis":'z', "first":i*3, "last":i*3+3}))
		
	imgs=[]
	for r in rots:
		ms=[]
		for ptcl in ptcls:
			if r!=0:
				pt=ptcl.copy()
				pt.rotate(0,0,r)
			else:
				pt=ptcl
			p=pt.numpy().copy()
			ms.append(p)
	
		m=np.array(ms)/3.
		m=m.transpose((2,1,0))
		#print(np.max(m), np.min(m), np.std(m))
		m=np.clip(m, -4,4)
		imgs.append(m)
		
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
		
class EMImageList(QtWidgets.QWidget):
	def __init__(self, options):
		QtWidgets.QWidget.__init__(self)
		self.win_size=[720,480]
		self.setMinimumSize(self.win_size[0], self.win_size[1])
		self.gbl = QtWidgets.QGridLayout(self)
		# Micrograph list
		self.imglst=QtWidgets.QTableWidget(1, 5, self)
		#self.imglst.verticalHeader().hide()
		self.imglst.setColumnWidth(0,50)
		self.imglst.setColumnWidth(1,200)
		self.gbl.addWidget(self.imglst, 0,0)
		self.options=options
		self.update_list()
	

class EMTomobox(QtWidgets.QMainWindow):

	def __init__(self,application,options,datafile=None):
		QtWidgets.QWidget.__init__(self)
		self.setMinimumSize(700,200)
		
		#### load references first
		self.reffile="info/boxrefs3d.hdf"
		if os.path.isfile(self.reffile):
			self.references=EMData.read_images(self.reffile)
		else:
			self.references=[]
			
		self.path="tomograms/"
		self.setCentralWidget(QtWidgets.QWidget())
		self.gbl = QtWidgets.QGridLayout(self.centralWidget())
		
		self.imglst=QtWidgets.QTableWidget(1, 5, self)
		#self.imglst.verticalHeader().hide()
		for i,w in enumerate([50,200,70,70,70]):
			self.imglst.setColumnWidth(i,w)
		self.imglst.setMinimumSize(450, 100)
		self.gbl.addWidget(self.imglst, 0,4,10,10)
		
		self.bt_new=QtWidgets.QPushButton("New")
		self.bt_new.setToolTip("Build new neural network")
		self.gbl.addWidget(self.bt_new, 0,0,1,2)
		
		self.bt_load=QtWidgets.QPushButton("Load")
		self.bt_load.setToolTip("Load neural network")
		self.gbl.addWidget(self.bt_load, 1,0,1,2)
		
		self.bt_train=QtWidgets.QPushButton("Train")
		self.bt_train.setToolTip("Train neural network")
		self.gbl.addWidget(self.bt_train, 2,0,1,2)
		
		self.bt_save=QtWidgets.QPushButton("Save")
		self.bt_save.setToolTip("Save neural network")
		self.gbl.addWidget(self.bt_save, 3,0,1,2)
		
		self.bt_apply=QtWidgets.QPushButton("Apply")
		self.bt_apply.setToolTip("Apply neural network")
		self.gbl.addWidget(self.bt_apply, 4,0,1,2)
		
		self.bt_chgbx=QtWidgets.QPushButton("ChangeBx")
		self.bt_chgbx.setToolTip("Change box size")
		self.gbl.addWidget(self.bt_chgbx, 5,0,1,2)
		
		self.bt_applyall=QtWidgets.QPushButton("ApplyAll")
		self.bt_applyall.setToolTip("Apply to all tomograms")
		self.gbl.addWidget(self.bt_applyall, 6,0,1,2)
		
		self.box_display = QtWidgets.QComboBox()
		self.box_display.addItem("References")
		self.box_display.addItem("Particles")
		self.gbl.addWidget(self.box_display, 0,2,1,1)
		
		
		self.bt_new.clicked[bool].connect(self.new_nnet)
		self.bt_load.clicked[bool].connect(self.load_nnet)
		self.bt_train.clicked[bool].connect(self.train_nnet)
		self.bt_save.clicked[bool].connect(self.save_nnet)
		self.bt_apply.clicked[bool].connect(self.apply_nnet)
		self.bt_chgbx.clicked[bool].connect(self.change_boxsize)
		self.bt_applyall.clicked[bool].connect(self.apply_nnet_all)
		self.box_display.currentIndexChanged.connect(self.do_update)

		self.val_targetsize=TextBox("TargetSize", 1)
		self.gbl.addWidget(self.val_targetsize, 1,2,1,1)
		
		self.val_learnrate=TextBox("LearnRate", 1e-4)
		self.gbl.addWidget(self.val_learnrate, 2,2,1,1)
		
		self.val_ptclthr=TextBox("PtclThresh", 0.8)
		self.gbl.addWidget(self.val_ptclthr, 3,2,1,1)
		
		self.val_circlesize=TextBox("CircleSize", 24)
		self.gbl.addWidget(self.val_circlesize, 4,2,1,1)
		
		self.val_niter=TextBox("Niter", 20)
		self.gbl.addWidget(self.val_niter, 5,2,1,1)
		
		self.val_posmult=TextBox("PosMult", 1)
		self.gbl.addWidget(self.val_posmult, 6,2,1,1)
		
		self.val_lossfun = QtWidgets.QComboBox()
		self.val_lossfun.addItem("Sum")
		self.val_lossfun.addItem("Max")
		self.gbl.addWidget(self.val_lossfun, 7,2,1,1)
		
		self.options=options
		self.app=weakref.ref(application)
		
		self.nnet=None
		global tf
		tf=import_tensorflow(options.gpuid)
		
		
		self.nnetsize=96
		if len(self.references)==0:
			self.boxsize=self.nnetsize
		else:
			self.boxsize=self.references[0]["boxsz"]
		self.thick=9
		self.datafile=""
		self.data=None
		
		if not os.path.isdir("neuralnets"):
			os.mkdir("neuralnets")
		
		self.imgview = EMImage2DWidget()
		self.boxesviewer=[EMImageMXWidget(), EMImageMXWidget()]
		self.boxesviewer[0].setWindowTitle("Negative")
		self.boxesviewer[1].setWindowTitle("Positive")
		self.ptclviewer=EMImageMXWidget()
		self.ptclviewer.setWindowTitle("Particles")
		
		for boxview in (self.boxesviewer+[self.ptclviewer]):
			boxview.usetexture=False
			boxview.show()
			boxview.set_mouse_mode("App")
			boxview.rzonce=True
		
		self.boximages=[]
		self.ptclimages=[]
		for img in self.references:
			self.add_boximage(img)
		
		self.imgview.mouseup.connect(self.on_tomo_mouseup)
		self.imgview.keypress.connect(self.key_press)
		self.boxesviewer[0].mx_image_selected.connect(self.on_boxpos_selected)
		self.boxesviewer[1].mx_image_selected.connect(self.on_boxneg_selected)
		self.ptclviewer.mx_image_selected.connect(self.on_ptcl_selected)
		self.imglst.cellClicked[int, int].connect(self.on_list_selected)
		
		self.boxshapes=BoxShapes(img=self.imgview, points=[] )
		self.imgview.shapes = {0:self.boxshapes}

		glEnable(GL_POINT_SMOOTH)
		glEnable(GL_LINE_SMOOTH );
		glEnable(GL_POLYGON_SMOOTH );
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		self.trainset=[]
		self.tomoinp=[]
		self.segout=None
		self.curinfo=None
		
		self.update_list()
		self.do_update()
			
	
	def update_list(self):
		#### update file list
		files=natural_sort([os.path.join(self.path,f) for f in os.listdir(self.path)])
		self.imginfo=[]
		for i,name in enumerate(files):
			basename=base_name(name)
			jsfile=info_name(basename)
			if not os.path.isfile(jsfile):
				continue
			info={"id":len(self.imginfo), "name":name, "basename":basename, "nptcl":0, "gref":0, "bref":0, "clsid":-1}
			js=js_open_dict(jsfile)
			clsid=-1
			if ("class_list" in js) and ("boxes_3d" in js):
				cls=js["class_list"]
				
				for k in list(cls.keys()):
					vname=str(cls[k]["name"])
					if vname==self.options.label:
						clsid=int(k)
						break
					
				ptcls=[p for p in js["boxes_3d"] if p[5]==clsid]
				info["nptcl"]=len(ptcls)
			
			refs=[r for r in self.references if base_name(r["fromtomo"])==basename]
			info["gref"]=len([r for r in refs if r["label"]==1])
			info["bref"]=len([r for r in refs if r["label"]==0])
			info["clsid"]=clsid
				
			self.imginfo.append(info)
			js.close()
			
		#print(self.imginfo)
		self.imglst.clear()
		self.imglst.setRowCount(len(self.imginfo))
		self.imglst.setColumnCount(5)
		self.imglst.setHorizontalHeaderLabels(["ID", "FileName", "Ptcls", "GoodRef", "BadRef"])
		self.imglst.setColumnHidden(0, True)
		for i,info in enumerate(self.imginfo):
			#### use Qt.EditRole so we can sort them as numbers instead of strings
			it=QtWidgets.QTableWidgetItem()
			it.setData(Qt.EditRole, int(info["id"]))
			self.imglst.setItem(i,0,it)
			
			self.imglst.setItem(i,1,QtWidgets.QTableWidgetItem(str(info["basename"])))
			
			it=QtWidgets.QTableWidgetItem()
			it.setData(Qt.EditRole, info["nptcl"])
			self.imglst.setItem(i,2, it)
			
			
			it=QtWidgets.QTableWidgetItem()
			it.setData(Qt.EditRole, info["gref"])
			self.imglst.setItem(i,3, it)
			
			it=QtWidgets.QTableWidgetItem()
			it.setData(Qt.EditRole, info["bref"])
			self.imglst.setItem(i,4, it)
			
			
		self.imglst.setVerticalHeaderLabels([str(i) for i in range(len(self.imginfo))])
	
	
	def on_list_selected(self, row, col):
		if self.curinfo:
			self.save_points()
			self.update_list()
		
		idx=self.imglst.item(row, 0).text()
		self.curinfo=self.imginfo[int(idx)]
		self.set_data(self.curinfo["name"])
		
		#hdr=EMData(info["filename"], 0,True)
		#iz=hdr["nz"]//2
		#e=EMData(info["filename"], 0, False, Region(0,0,iz, hdr["nx"], hdr["ny"],1))
		#fac=float(hdr["nx"])/self.bt_show2d.width()*1.01
		#e.process_inplace('math.fft.resample',{"n":fac})
		#self.wg_thumbnail.set_data(e)
		
		
	def set_data(self, datafile):
		if self.datafile==datafile:
			return
		print("Reading {}...".format(datafile))
		self.datafile=datafile
		self.data=EMData(datafile)
		self.data.mult(self.options.mult)
		self.imgview.setWindowTitle(base_name(datafile))
		self.imgview.list_idx=self.data["nz"]//2
		self.imgview.set_data(self.data)
		self.imgview.show()
		self.infofile=info_name(datafile)
		
		js=js_open_dict(self.infofile)
		apix_cur=apix=self.data["apix_x"]
		apix_unbin=js["apix_unbin"]
		self.apix_scale=apix_cur/apix_unbin
		self.tomocenter= np.array([self.data["nx"],self.data["ny"],self.data["nz"]])/2
		self.ptclimages=[]
		if "boxes_3d" in js:
			ptcls=np.array([p[:3] for p in js["boxes_3d"] if p[5]==self.curinfo["clsid"]])
			if len(ptcls)>0:
				ptcls= ptcls  / self.apix_scale + self.tomocenter
				for p in ptcls.tolist():
					self.add_ptcls(p[0], p[1], p[2])
		
		js.close()
		self.tomoinp=[]
		self.segout=None
		self.do_update()
	
	def change_boxsize(self):
		size,ok=QtWidgets.QInputDialog.getText( self, "Box size", "Enter a new box size:")
		if not ok : return
		size=int(size)
		if size<self.nnetsize:
			print("Cannot set box size ({}) smaller than Network input size ({}). Stop.".format(size, self.nnetsize))
			return
		self.boxsize=size
		print("Updating references...")
		oldref=[r for r in self.references]
		self.references=[]
		self.boximages=[]
		curdatafile=self.datafile
		
		tomolist={r["fromtomo"] for r in oldref}
		print("  {} references from {} tomograms...".format(len(oldref), len(tomolist)))
		for tomo in tomolist:
			imgs=[r for r in oldref if r["fromtomo"]==tomo]
			self.datafile=tomo
			self.data=EMData(tomo)
			self.data.mult(self.options.mult)
			for m in imgs:
				p=m["pos"]
				self.add_reference(m["label"], p[0],p[1],p[2])
				
				
		oldref=[]
		self.datafile=curdatafile
		self.do_update()
		self.trainset= []
		self.tomoinp=[]
		if len(self.datafile)>0:
			self.curdata=EMData(self.datafile)
		else:
			self.curdata=None
		
	def new_nnet(self):
		print("New network..")
		self.nnet=NNet(boxsize=self.nnetsize, thick=self.thick)
		
	def train_nnet(self):
		if self.nnet==None:
			self.new_nnet()
			
		if len(self.trainset)==0:
			if len(self.references)==0:
				print("No references.")
				return
			
			print("Preparing training set...")
			labels=[b["label"] for b in self.references]
			ngood=np.mean(labels)
			nc=500/len(labels)
			ncopy=[int(round(nc/(1-ngood)+1)), int(round(nc/ngood+1))]
			#print(ncopy)
			imgs=[]
			labs=[]
			for i,p in enumerate(self.references):
				lb=labels[i]
				m=get_image(p, thick=self.thick, ncopy=ncopy[lb])
				imgs.extend(m)
				labs.extend([lb]*ncopy[lb])
			imgs=np.array(imgs, dtype=np.float32)
			labs=np.array(labs, dtype=np.float32)
			#print(len(labs),np.sum(labs==0), np.sum(labs==1), np.mean(labels))
			self.trainset=(imgs, labs)
			
		dataset = tf.data.Dataset.from_tensor_slices(self.trainset)
		dataset=dataset.shuffle(500).batch(32)
		usemax=(self.val_lossfun.currentText()=="Max")
		self.nnet.do_training(
			dataset, 
			learnrate=self.val_learnrate.getval(), 
			niter=int(self.val_niter.getval()),
			tarsz=self.val_targetsize.getval(),
			usemax=usemax, 
			posmult=self.val_posmult.getval()
			)
		
		self.segout=None
		print("Generating output...")
		
		imgs,labs=self.trainset
		idx=np.arange(len(imgs))
		np.random.shuffle(idx)
		idx=idx[:100]
		out=self.nnet.apply_network(imgs[idx])[:,:,:,0]
		
		#outval=self.nnet.predict_class(np.array(imgs[idx], dtype=np.float32), usemax=usemax)
		#outval=np.array(outval)[:,:,0]
		
		fname="neuralnets/trainouts.hdf"
		if os.path.isfile(fname):
			os.remove(fname)
			
		bx=self.nnetsize
		sz=out.shape[1]
		for i,o in enumerate(out):
			m=np.mean(imgs[idx[i]], axis=-1)/3.0
			m=from_numpy(m)
			m["score"]=[int(labs[idx[i]]),int(labs[idx[i]])]
			m.write_image(fname, -1)
			m=from_numpy(o)
			m=m.get_clip(Region(sz//2-bx//2, sz//2-bx//2, bx, bx))
			m.scale(4)
			#m["score"]=[float(outval[0,i]),float(outval[1,i])]
			m.write_image(fname, -1)
		print("Output written to {}...".format(fname))
		
	
	def apply_nnet_all(self):
		if self.nnet==None:
			print("Neural network not initialized...")
			return
		
		modifiers = QtWidgets.QApplication.keyboardModifiers()
		skipexist=(modifiers == QtCore.Qt.ShiftModifier)
		if skipexist:
			print("Skipping tomograms with particles")
			
		for i,info in enumerate(self.imginfo):
			self.curinfo=info
			if skipexist and info["nptcl"]>0:
				print("Skipping {}..".format(info["name"]))
				continue
			self.set_data(info["name"])
			self.apply_nnet()
		
	def apply_nnet(self):
		bx=self.nnetsize//8
		thk=self.thick//8
		if self.data==None:
			print("No data loaded")
			return
		if self.nnet==None:
			print("Neural network not initialized...")
			return
		
		ts=[self.data["nx"], self.data["ny"], self.data["nz"]]
		if self.nnetsize==self.boxsize:
			data=self.data
		else:
			data=self.data.copy()
			scale=self.nnetsize/self.boxsize
			data.scale(scale)
			ts2=[int(t*scale) for t in ts]
			#print(ts, ts2)
			data=data.get_clip(Region(ts[0]//2-ts2[0]//2, ts[1]//2-ts2[1]//2, ts[2]//2-ts2[2]//2,  ts2[0], ts2[1], ts2[2]))
			ts=ts2
			
		if len(self.tomoinp)==0:
			print("Preparing input...")
			
			m=[]
			for i in range(0, ts[2],4):
				m.extend(get_image(data, pos=[ts[0]//2, ts[1]//2, i], thick=self.thick, ncopy=1))
				
			self.tomoinp=np.array(m)

		if self.segout==None:
			print("Applying...")
			out=self.nnet.apply_network(self.tomoinp)[:,:,:,0]
			#o=out.reshape((-1, 4, out.shape[1], out.shape[2]))
			#o=np.mean(o, axis=1)
			o=out.transpose(0,2,1).copy()
			o=from_numpy(o)
			o=o.get_clip(Region(o["nx"]//2-ts[0]//8, o["ny"]//2-ts[1]//8, o["nz"]//2-ts[2]//8, ts[0]//4, ts[1]//4, ts[2]//4))
			#o.scale(self.boxsize/self.nnetsize)
			o.process_inplace("mask.zeroedge3d",{"x0":bx,"x1":bx,"y0":bx,"y1":bx,"z0":thk,"z1":thk})
			o.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.5})
			self.segout=o.copy()
			o.write_image("neuralnets/segout.hdf")
		
		print("Finding boxes...")
		o=self.segout.copy()
		self.ptclimages=[]
		o.process_inplace("mask.onlypeaks",{"npeaks":5})
		img=o.numpy().copy()
		img[img<self.val_ptclthr.getval()]=0
		
		pnew=np.array(np.where(img>self.val_ptclthr.getval())).T
		val=img[pnew[:,0], pnew[:,1], pnew[:,2]]
		pnew=pnew[np.argsort(-val)]

		pnew=pnew[:, ::-1]*4#+4
		pnew=pnew*self.boxsize/self.nnetsize

		dthr=self.val_circlesize.getval()
		
		
		tree=KDTree(pnew)

		tokeep=np.ones(len(pnew), dtype=bool)
		for i in range(len(pnew)):
			if tokeep[i]:
				k=tree.query_ball_point(pnew[i], dthr)
				tokeep[k]=False
				tokeep[i]=True
			
		#print(np.sum(tokeep))
		pts=pnew[tokeep].tolist()
		#scr=pkscore[tokeep]
		
		
		#dst=scipydist.cdist(pnew, pnew)+(np.eye(len(pnew))*dthr*100)
		#tokeep=np.ones(len(dst), dtype=bool)
		#for i in range(len(dst)):
			#if tokeep[i]:
				#tokeep[dst[i]<dthr]=False

		#pts=pnew[tokeep].tolist()

		for p in pts:
			self.add_ptcls(p[0], p[1], p[2])
			
		print("Found {} particles...".format(len(pts)))
		self.do_update()
		self.save_points()
		self.update_list()
		#print(self.boxshapes.points)
		
	def load_nnet(self):
		self.nnet=NNet.load_network("neuralnets/nnet_save.hdf", boxsize=self.nnetsize, thick=self.thick)
		
	def save_nnet(self):
		self.nnet.save_network("neuralnets/nnet_save.hdf")
		
	def key_press(self, event):
		#print(event.key())
		if event.key()==96:
			self.imgview.increment_list_data(1)
		elif event.key()==49:	
			self.imgview.increment_list_data(-1)
		self.imgview.shapechange=1
		self.imgview.updateGL()
		return
	
	def on_tomo_mouseup(self, event):
		x,y=self.imgview.scr_to_img((event.x(),event.y()))		
		x,y =np.round(x), np.round(y)
		
		if not event.button()&Qt.LeftButton:
			return
		
		if  event.modifiers()&Qt.ShiftModifier:
			### delete point
			z=self.imgview.list_idx
			pos=np.array([x,y,z])
			mode=self.box_display.currentText()
			if mode=="Particles":
				pts=np.array([b["pos"] for b in self.ptclimages])
			else:
				pts=np.array([img["pos"] for img in self.references])
			if len(pts)==0:
				return
			dst=np.linalg.norm(pts[:,:3]-pos, axis=1)
			
			if np.sum(dst<self.val_circlesize.getval())>0:
				idx=np.argsort(dst)[0]
				if mode=="Particles":
					self.ptclimages=[p for i,p in enumerate(self.ptclimages) if i!=idx]
					self.save_points()
				else:
					self.boxshapes.points.pop(idx)
					self.references=[b for ib,b in enumerate(self.references) if ib!=idx]
					self.boximages=[b for ib,b in enumerate(self.boximages) if ib!=idx]
					
			self.do_update()
			
		else:
			#### hold control to add negative references
			label=int(event.modifiers()&Qt.ControlModifier)
			label=int(label==0)
			self.add_reference(label, x, y)
			self.do_update()
		
	def on_boxpos_selected(self,event,lc):
		self.on_box_selected(event, lc, ib=0)
		
	def on_boxneg_selected(self,event,lc):
		self.on_box_selected(event, lc, ib=1)
		
	def on_box_selected(self,event, lc, ib):
		ic=int(1-ib)
		if event.modifiers()&Qt.ShiftModifier:
			if event.modifiers()&Qt.ControlModifier:
				torm=len(self.boxesviewer[ib].data)-lc[0]
			else:
				torm=1
			
			for i in range(torm):
				img=self.boxesviewer[ib].data[lc[0]+i]
				pos=np.array(img["pos"])
				self.rm_reference(pos)
				self.boximages=[b for b in self.boximages if b!=img]
			self.do_update()
			
	def on_ptcl_selected(self,event,lc):
		
		
		if  event.modifiers()&Qt.ShiftModifier:
			## delete ptcl
			if event.modifiers()&Qt.ControlModifier:
				self.ptclimages=[p for i,p in enumerate(self.ptclimages) if i<lc[0]]
			else:
				self.ptclimages=[p for i,p in enumerate(self.ptclimages) if i!=lc[0]]
			self.do_update()
			
			
		else:
			## add to good/bad refs
			if event.modifiers()&Qt.ControlModifier:
				ic=0
			else:# event.modifiers()&Qt.ShiftModifier:
				ic=1
			
			img=self.ptclimages[lc[0]]
			pos=img["pos"]
			self.add_reference(ic, pos[0], pos[1], pos[2])
			self.do_update()
			
			
	def rm_reference(self, pos):
		pts=np.array([b["pos"] for b in self.references])
		dst=np.linalg.norm(pts-pos, axis=1)
		if np.sum(dst<1)>0:
			idx=int(np.where(dst<1)[0][0])
			self.references=[b for i,b in enumerate(self.references) if i!=idx]
			self.trainset=[]
	
	def add_reference(self, label, x, y, z=None):
		if z==None:
			z=self.imgview.list_idx
		
		
		
		if self.nnetsize==self.boxsize:
			b=self.boxsize
			t=self.thick
			img=self.data.get_clip(Region(x-b//2, y-b//2, z-t//2, b, b, t))
			
		else:
			b=self.boxsize
			img=self.data.get_clip(Region(x-b//2, y-b//2, z-b//2, b, b, b))
			img.scale(self.nnetsize/self.boxsize)
			b2=self.nnetsize
			t=self.thick
			img=img.get_clip(Region(b//2-b2//2, b//2-b2//2, b//2-t//2,  b2, b2, t))
			
		img["pos"]=[x,y,z]
		img["fromtomo"]=self.datafile
		img["label"]=label
		img["boxsz"]=self.boxsize
		self.references.append(img)
		self.add_boximage(img)
		self.trainset=[]
		
	def add_boximage(self, img):
		img2d=img.get_clip(Region(0, 0, img["nz"]//2, img["nx"], img["nx"], 1))
		self.boximages.append(img2d)
			
	def add_ptcls(self, x, y, z):
		b=self.boxsize
		img=self.data.get_clip(Region(x-b//2, y-b//2, z, b, b, 1))
		img["pos"]=[x,y,z]
		self.ptclimages.append(img)

	def save_points(self):
		if self.data==None:
			return
		#print("Saving particles...")
		js=js_open_dict(self.infofile)
		label=self.options.label
		pts=np.array([b["pos"] for b in self.ptclimages])
		if len(pts)>0:
			pts=(pts - self.tomocenter) * self.apix_scale
		pts=np.round(pts).astype(int)
		clsid=self.curinfo["clsid"]
		print('save particles to {} : {}'.format(clsid, label))
		if not "class_list" in js:
			js["class_list"]={}
		if not "boxes_3d" in js:
			js["boxes_3d"]=[]
			
		clslst=js["class_list"]
		if clsid==-1:
			if len(clslst.keys())==0:
				clsid=0
			else:
				clsid=max([int(k) for k in clslst.keys()])+1
			clslst[str(clsid)]={"name":label, "boxsize": int(self.boxshapes.circlesize*8)}
			js["class_list"]=clslst
			self.curinfo["clsid"]=clsid
			
		boxes=[b for b in js["boxes_3d"] if b[5]!=clsid]
		boxes=boxes+[[p[0], p[1], p[2], "tomobox", 0.0, clsid] for p in pts.tolist()]
		js["boxes_3d"]=boxes
		js.close()
		#self.update_list()
		
	def clear_points(self):
		return
	
	#def box_display_changed(self):
		#print(self.box_display.currentText())
		
	
	def do_update(self):
		if self.box_display.currentText()=="Particles":
			pts=[b["pos"]+[1] for b in self.ptclimages]
		else:
			pts=[b["pos"]+[b["label"]] for b in self.references if b["fromtomo"]==self.datafile]
			
		
		#print(pts)
		self.boxshapes.points=pts
		self.boxshapes.circlesize=self.val_circlesize.getval()
		self.imgview.shapechange=1
		self.imgview.updateGL()
		
		self.ptclviewer.set_data(self.ptclimages)
		self.ptclviewer.update()
		for i in [0,1]:
			self.boxesviewer[i].set_data([b for b in self.boximages if b["label"]==i])
			self.boxesviewer[i].update()
		
	
	def save_references(self):
		if os.path.isfile(self.reffile):
			os.remove(self.reffile)
		
		for ref in self.references:
			ref.write_image(self.reffile, -1)
		
		
	def closeEvent(self, event):
		self.save_references()
		self.imgview.close()
		self.save_points()
		#self.imagelist.close()
		for b in  (self.boxesviewer+[self.ptclviewer]):
			b.close()
if __name__ == '__main__':
	main()
	
