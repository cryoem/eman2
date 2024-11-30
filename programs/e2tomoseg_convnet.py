#!/usr/bin/env python
# Muyuan July 2015
from past.utils import old_div
from future import standard_library
standard_library.install_aliases()
from builtins import range
from builtins import object
import sys
import random
import numpy as np
from EMAN2 import *
import pickle
import time
import threading
from queue import Queue
from multiprocessing import Array

def import_theano():
	global theano,T,conv2d,pool
	import theano
	import theano.tensor as T
	try:
		from theano.tensor.nnet import conv2d
	except:
		from theano.tensor.nnet.conv import conv2d
	from theano.tensor.signal import pool

def import_tensorflow(gpuid=None):
	import os
	global tf
	if gpuid!=None: #### decide which gpu to use
		os.environ["CUDA_VISIBLE_DEVICES"]=str(gpuid)
	import tensorflow as tf
	os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3' #### reduce log output


def main():

	usage="""Segment a tomograph using convolutional neural network. Please run this program from the GUI in e2projectmanager.py."""
	#print usage
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	#parser.add_header(name="tmpheader", help='temp label', title="### This program is NOT avaliable yet... ###", row=0, col=0, rowspan=1, colspan=2, mode="train,test")
	parser.add_argument("--trainset",help="Training set.", default=None, guitype='filebox', browser="EMBrowserWidget(withmodal=True, startpath='particles', dirregex='*_trainset.hdf')",  row=1, col=0,rowspan=1, colspan=3, mode="train")
	parser.add_argument("--from_trained", type=str,help="Train from an existing network", default=None,guitype='filebox',browser="EMBrowserWidget(withmodal=True, startpath='neuralnets', dirregex='nnet_save*')", row=2, col=0, rowspan=1, colspan=3, mode="train")
	parser.add_argument("--nnet", type=str,help="Trained network input (nnet_save_xx.hdf)", default=None,guitype='filebox',browser="EMBrowserWidget(withmodal=True, startpath='neuralnets', dirregex='nnet_save*')", row=2, col=0, rowspan=1, colspan=3, mode="test")
	#parser.add_argument("--netout", type=str,help="Output neural net file name", default="nnet_save.hdf",guitype='strbox', row=3, col=0, rowspan=1, colspan=3, mode="train")
	parser.add_argument("--nettag", type=str,help="Tag of the output neural net file. Will use the tag of good particles in training set by default.", default="",guitype='strbox', row=3, col=0, rowspan=1, colspan=3, mode="train")
	
	parser.add_argument("--learnrate", type=float,help="Learning rate ", default=1e-4, guitype='floatbox', row=4, col=0, rowspan=1, colspan=1, mode="train")
	parser.add_argument("--niter", type=int,help="Training iterations", default=20, guitype='intbox', row=4, col=1, rowspan=1, colspan=1, mode="train")
	parser.add_argument("--ncopy", type=int,help="Number of copies for each particle", default=1, guitype='intbox', row=5, col=0, rowspan=1, colspan=1, mode="train")
	parser.add_argument("--batch", type=int,help="Batch size for the stochastic gradient descent. Default is 20.", default=20, guitype='intbox', row=5, col=1, rowspan=1, colspan=1, mode="train")
	parser.add_argument("--nkernel", type=str,help="Number of kernels for each layer, from input to output. The number of kernels in the last layer must be 1. ", default="40,40,1", guitype='strbox', row=6, col=0, rowspan=1, colspan=1, mode="train")
	parser.add_argument("--ksize", type=str,help="Width of kernels of each layer, the numbers must be odd. Note the number of layers should be the same as the nkernel option. ", default="15,15,15", guitype='strbox', row=6, col=1, rowspan=1, colspan=1, mode="train")
	parser.add_argument("--poolsz", type=str,help="Pooling size for each layer. Note the number of layers should be the same as the nkernel option. ", default="2,1,1", guitype='strbox', row=7, col=0, rowspan=1, colspan=1, mode="train")
	#parser.add_argument("--weightdecay", type=float,help="Weight decay. Used for regularization.", default=1e-6, guitype='floatbox', row=7, col=1, rowspan=1, colspan=1, mode="train")
	parser.add_argument("--trainout", action="store_true", default=False ,help="Output the result of the training set", guitype='boolbox', row=9, col=0, rowspan=1, colspan=1, mode='train[True]')
	parser.add_argument("--training", action="store_true", default=False ,help="Doing training", guitype='boolbox', row=9, col=1, rowspan=1, colspan=1, mode='train[True]')
	parser.add_argument("--tomograms", type=str,help="Tomograms input.", default=None,guitype='filebox',browser="EMBrowserWidget(withmodal=True, startpath='tomograms', multiselect=True)", row=1, col=0, rowspan=1, colspan=3, mode="test")
	parser.add_argument("--applying", action="store_true", default=False ,help="Applying the neural network on tomograms", guitype='boolbox', row=4, col=0, rowspan=1, colspan=1, mode='test[True]')
	#parser.add_argument("--dream", action="store_true", default=False ,help="Iterativly applying the neural network on noise")
	#parser.add_argument("--to3d", action="store_true", default=True ,help="convert to result to 3D.", guitype='boolbox', row=5, col=1, rowspan=1, colspan=1, mode='test')
	parser.add_argument("--outtag", type=str,help="Tag of the segmentation output. When left empty, the segmentation will be saved to 'segmentations/<tomogram name>__<neural network tag>_seg.hdf'. When set, the output will be written to 'segmentations/<tomogram name>__<outtag>.hdf'", default="", guitype='strbox', row=3, col=0, rowspan=1, colspan=1, mode="test")
	parser.add_argument("--threads", type=int,help="Number of thread to use when applying neural net on test images. Not used during trainning", default=12, guitype='intbox', row=10, col=0, rowspan=1, colspan=1, mode="test")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--device", type=str, help="For Convnet training only. Pick a device to use. chose from cpu, gpu, or gpuX (X=0,1,...) when multiple gpus are available. default is cpu",default="cpu",guitype='strbox', row=7, col=1, rowspan=1, colspan=1,mode="train,test")


	(options, args) = parser.parse_args()
	E2n=E2init(sys.argv,options.ppid)
	time00=time.time()
	
	if options.from_trained==None and options.nnet!=None:
		options.from_trained=options.nnet
	
	#### parse the options.
	nkernel=[int(i) for i in options.nkernel.split(',')]
	ksize=[int(i) for i in options.ksize.split(',')]
	poolsz=[int(i) for i in options.poolsz.split(',')]
	options.labelshrink=np.prod(poolsz)
	
		
	if "CUDA_VISIBLE_DEVICES" in os.environ:
		print("CUDA_VISIBLE_DEVICES is already set as environment variable. This will override the --device option...")
	else:
		try: options.device=options.device.lower()
		except: pass
		if options.device=="gpu":
			print("Using GPU...")
		elif options.device.startswith("gpu"):
			try:
				gid=int(options.device[3:])
				print("Using GPU #{}..".format(gid))
				os.environ["CUDA_VISIBLE_DEVICES"]="{}".format(gid)
			except:
				print("Cannot parse {}, will use CPU instead...".format(options.device))
				os.environ["CUDA_VISIBLE_DEVICES"]=""
			
		elif options.device=="cpu":
			print("Using CPU...")
			os.environ["CUDA_VISIBLE_DEVICES"]=""
		else:
			print("Cannot parse {}, will use CPU instead...".format(options.device))
			os.environ["CUDA_VISIBLE_DEVICES"]=""

	import_tensorflow()
	
	
	if options.applying:
		dirname="segmentations"
		try: os.mkdir(dirname)
		except: pass
	
		tomos=options.tomograms.split(',')
		### there might be a minor issue if the list of tomograms have different x-y sizes...
		print("Loading the Neural Net...")
		e=EMData(tomos[0], 0, True)
		tsz=max(e["nx"], e["ny"])
		convnet=StackedConvNet_tf.load_network(options.from_trained, imgsz=tsz, bsz=1)
		
		for tm in tomos:
			print("Starting on {}...".format(tm))
			bn=base_name(tm)
			if options.outtag=="":	
				nn=options.nnet
				nn=nn[nn.rfind("__")+2:-4]
				outname="segmentations/{}__{}_seg.hdf".format(bn, nn)
			else:
				outname="segmentations/{}__{}.hdf".format(bn, options.outtag)
				
			segout=apply_neuralnet(convnet, options, tm)
			segout.write_image(outname)
			print("Output segmentation written to {}...".format(outname))
		print("Done.")
		print("Total time: ", time.time()-time00)
	
		E2end(E2n)
		return
	
	
	if options.nettag=="":
		tag=options.trainset
		options.nettag=tag[tag.rfind('__')+2:-4].replace("_trainset","")
	
	
	if options.trainset==None:
		print("No training set input...exit.")
		return
	
	
	#session=tf.Session()
	
	print("loading particles...")
	data, labels, shape, ntrain=load_particles(options.trainset,options.labelshrink,options.ncopy)
	batch_size=options.batch
	
	if options.from_trained!=None:
		convnet=StackedConvNet_tf.load_network(options.from_trained, imgsz=shape[0], bsz=batch_size)

	else:
		print("Setting up model...")
		kernels=[(k,ksize[i],poolsz[i]) for i,k in enumerate(nkernel)]
		convnet = StackedConvNet_tf(kernels, shape[0], batch_size)
		#session.run(tf.global_variables_initializer())

	
	
	#convnet.save_network(options.netout, session, options)
	
	if (options.niter>0):	
		
		convnet.do_training(data, labels, shuffle=False, learnrate=options.learnrate, niter=options.niter)
		
		
		
	dirname="neuralnets"
	try: os.mkdir(dirname)
	except: pass
	options.netout="{}/nnet_save__{}.hdf".format(dirname, options.nettag)

	if options.trainout:
		fname="{}/trainout_nnet_save__{}.hdf".format(dirname, options.nettag)
		convnet.write_output_train(fname, writelabel=True)
		
	convnet.save_network(options.netout, options)
	print("Done")
	print("Total time: {:.1f} s".format(time.time()-time00))
	E2end(E2n)

def run(cmd):
	print(cmd)
	launch_childprocess(cmd)

	
def apply_neuralnet(convnet, options, tomogram=None):
	
	if tomogram==None:
		tomogram=options.tomograms
	
	nframe=EMUtil.get_image_count(tomogram)
	is3d=False
	### deal with 3D volume or image stack
	e=EMData(tomogram, 0, True)
	apix=e["apix_x"]
	if nframe==1:
		nframe=e["nz"]
		if nframe>1:
			#### input data is 3D volume
			is3d=True
			
	enx,eny=e["nx"], e["ny"]
	tsz=max(enx,eny)
	
	output=EMData(e["nx"], e["ny"], nframe)
	output["tomogram_src"]=tomogram
	output["nnet_src"]=options.from_trained
	output["apix_x"]=apix
	output["apix_y"]=apix
	output["apix_z"]=apix
	
	#####################
	
	
	#print("Loading tomogram...")
	tomo_in=[]
	for nf in range(nframe):
		if is3d:
			e0=EMData(tomogram, 0, False, Region((enx-tsz)//2,(eny-tsz)//2,nf,tsz,tsz,1))
		else:
			e0=EMData(tomogram, nf, False, Region((enx-tsz)//2,(eny-tsz)//2,tsz,tsz))
		tomo_in.append(e0)
	
	for idx, img in enumerate(tomo_in):
		m=img.numpy()
		p=convnet.model.predict(m[None, :, :, None]/3.)
		p[p<0]=0
		cout=from_numpy(p[0,:,:,0])
		cout=cout.get_clip(Region((cout["nx"]-enx)//2,(cout["ny"]-eny)//2 ,enx, eny))
		cout.scale(int(options.labelshrink))
		output.insert_clip(cout, [0,0,idx])
		
		sys.stdout.write("\r  {}/{} finished.".format(idx+1, len(tomo_in)))
		sys.stdout.flush()
		
	return output
	
	
	
	
def load_particles(ptcls,labelshrink,ncopy=5, rng=None):
	if rng==None:
		rng=random
	num=old_div(EMUtil.get_image_count(ptcls),2)
	
	data=[]
	label=[]
	ntrain=-1
	for i in range(num):
		for nc in range(ncopy):
			ptl=EMData(ptcls,i*2)
			if ntrain<0 and ptl.get_attr_default("valid_set", 0)==1:
				ntrain=len(data)
			#ptl.process_inplace("threshold.belowtozero")
			if ncopy>1:
				tr=Transform()
				tr.set_rotation({"type":"2d","alpha":rng.random()*360.0})
				ptl.process_inplace("xform",{"transform":tr})
			
			
			ar=ptl.numpy().copy()
			#shp=np.shape(ar)
			data.append(ar)
			
			ptl=EMData(ptcls,i*2+1)
			#ptl.process_inplace("threshold.belowtozero")
			if ncopy>1:
				ptl.process_inplace("xform",{"transform":tr})
			if labelshrink>1:
				ptl.process_inplace("math.meanshrink",{'n':int(labelshrink)})
			ar=ptl.numpy().copy()
			#shp=np.shape(ar)
			label.append(ar)
	
	if ntrain<0: ntrain=len(data)
	print("{:d} particles loaded, {:d} in training set, {:d} in validation set".format(len(data), ntrain, len(data)-ntrain))
	data=np.asarray(data,dtype=np.float32)
	print(data.shape)
	print("Std of particles: ",np.std(data))
	#data/=np.std(data.flatten())*3  #np.max(np.abs(data))
	data/=3. ### so the range is roughly (-1,1)
	label=np.asarray(label,dtype=np.float32)
	label/=np.max(np.abs(label))
	
	
	header=EMData(ptcls,0,True)
	shape=[header["nx"],header["ny"],header["nz"]]
	return data, label, shape, ntrain


class StackedConvNet_tf(object):
		
	def __init__(self, kernels, imgsz=64, batchsz=10, meanout=False):
		
		#convnet = type('convnet', (), {})() ### an empty object
		nlayer=len(kernels)
		sz=imgsz
			
		outsz=sz//np.prod([k[2] for k in kernels])
		shrink=1
		layers=[]
		
		def min_act(x):
			#return tf.math.maximum(tf.math.minimum(x,1.0),0.0)
			return tf.math.minimum(x,1.0)
		
		for i in range(nlayer):
			nk, ks, pl=kernels[i]
				
			if i<nlayer-1:
				act="relu"
			else:
				act="linear"
			
			if i==0:
				kwg={"input_shape":(imgsz, imgsz, 1)}
			else:
				kwg={}
			
			layers.append(tf.keras.layers.Conv2D(nk, ks, padding="same", activation=act, kernel_initializer=tf.keras.initializers.TruncatedNormal(stddev=0.01), **kwg))
			if pl==2:
				layers.append(tf.keras.layers.MaxPooling2D())
				shrink*=2
		if meanout:
			layers.append(tf.keras.layers.MaxPooling2D(pool_size=(8,8)))
			layers.append(tf.keras.layers.Flatten())
			layers.append(tf.keras.layers.Dense(1))
				
		self.model=tf.keras.Sequential(layers)
		self.imgsz=imgsz
		self.outsz=self.model.get_layer(index=-1).output.shape[1]
		self.batchsize=batchsz
		self.kernels=kernels
		self.meanout=meanout
		self.labelshrink=shrink
		#print(layers[1].weights)
		#exit()

		

	def do_training(self, data, label, shuffle=False, learnrate=1e-4, niter=10):
		print("Training...")
		data=data.reshape((-1, self.imgsz, self.imgsz, 1))
		if not self.meanout:
			label=label.reshape((-1, self.outsz, self.outsz, 1))
		dataset = tf.data.Dataset.from_tensor_slices((data, label))
		dataset=dataset.shuffle(500).batch(self.batchsize)

		def calc_loss(yt, yp):
			if self.meanout:
				#ayp=tf.reduce_min(yp, axis=(1,2))
				ayp=yp
				ayp=tf.math.minimum(1.0, tf.math.maximum(-1.0, ayp))
			else:
				ayp=tf.math.minimum(yp,1.0)
			loss=tf.math.log(tf.math.reduce_mean((yt-ayp)**2))
			#loss+=sum([tf.nn.l2_loss(t) for t in self.model.weights])*1e-2
			return loss
		
		self.model.compile(optimizer=tf.keras.optimizers.Adam(learning_rate=learnrate), loss=calc_loss)
		for it in range(niter):
			self.model.reset_metrics()
			cost=[]
			for image, label in dataset:
				
				y=self.model.predict(image)
				result = self.model.train_on_batch(image, label)
				cost.append(result)
			print("iteration {}, cost {:.3f}".format(it, np.mean(cost)))

		
		self.dataset=dataset
		
		
	def save_network(self, fname, options=None):
		try: os.remove(fname)
		except: pass
		print("Saving the trained net to {}...".format(fname))
		
		weights=[]
		bias=[]
		for ly in self.model.layers:
			if len(ly.weights)==2:
				w,b=ly.weights
				weights.append(w.numpy())
				bias.append(b.numpy())
		
		
		wsz=int(weights[0].shape[0])

		hdr=EMData(wsz,wsz)
		hdr["nkernel"]=[k[0] for k in self.kernels]
		hdr["ksize"]=[k[1] for k in self.kernels]
		hdr["poolsz"]=[k[2] for k in self.kernels]
		hdr["imageshape"]=(self.batchsize,1,self.imgsz,self.imgsz)
		#hdr["meanout"]=self.meanout
		if options!=None:
			hdr["trained_from"]=options.trainset
		nlayer=len(self.kernels)

		hdr.write_image(fname,0)

		k=1
		for i in range(nlayer):
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
	def load_network(fname, imgsz=-1, bsz=0):
		
		hdr=EMData(fname,0)
		
		shp=hdr["imageshape"]
		if imgsz>0: sz=imgsz
		else:  sz=shp[-1]
		
		if bsz==0:  batchsize=shp[0]
		else: batchsize=bsz
		
		nkernel=hdr["nkernel"]
		ksize=hdr["ksize"]
		poolsz=hdr["poolsz"]

		bf=1

		kernels=[(k,ksize[i],poolsz[i]) for i,k in enumerate(nkernel)]
		
		nnet=StackedConvNet_tf(kernels, sz, batchsize)
		k=1
		for ly in nnet.model.layers:
			
			wts=ly.weights
			if len(wts)==0:
				continue
			
			e=EMData(fname,k)
			s=e["w_shape"]
			b=e.numpy().copy()
			k+=1
			# wts[1].assign(b, read_value=False)
			wts[1].assign(b) #read_value not used in Keras 3+

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
			# wts[0].assign(allw, read_value=False)
			wts[0].assign(allw) #read_value not used in Keras 3+
			
		return nnet	
		
	def write_output_train(self, outfile, ncopy=10, writelabel=False):
		sz=self.imgsz
		outsz=self.outsz
		try: os.remove(outfile)
		except: pass
		
		print("Writting network output of training set to {}...".format(outfile))
	
		k=0
		if self.meanout:
			print(self.model.layers[:-3])
			modelpred=tf.keras.Sequential(self.model.layers[:-3])
		else:
			modelpred=self.model
			
		for image, label in self.dataset:
			k+=1
			if k>ncopy: break
			
			result = modelpred.predict_on_batch(image)
			#print(result.shape)

			for i in range(len(image)):
				im=np.array(image[i][:,:,0])
				
				rs=np.array(result[i][:,:,0])
				#rs[rs<0]=0
				#print(np.std(rs))

				e0=from_numpy(im.copy())
				e0.process_inplace("normalize")
				e0.write_image(outfile, -1)
				
				if not self.meanout:
					lb=np.array(label[i][:,:,0])
					e1=from_numpy(lb.copy())
					e1=e1.get_clip(Region(-(sz-outsz)//2,-(sz-outsz)//2,sz,sz))
					e1.scale(sz/outsz)
					e1.write_image(outfile, -1)
				
				e1=from_numpy(rs.copy())
				e1=e1.get_clip(Region(-(sz-outsz)//2,-(sz-outsz)//2,sz,sz))
				e1.scale(sz/outsz)
				e1.write_image(outfile, -1)
				
			
		#print(np.mean(amp))

if __name__ == '__main__':
    main()
