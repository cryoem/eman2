#!/usr/bin/env python
from __future__ import print_function
from __future__ import division
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
	import tensorflow
	if tensorflow.__version__[0]=="2" : 
		import tensorflow.compat.v1 as tf
		tf.disable_eager_execution()
	else: import tensorflow as tf
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
	parser.add_argument("--device", type=str, help="For Convnet training only. Pick a device to use. chose from cpu, gpu, or gpuX (X=0,1,...) when multiple gpus are available. default is cpu",default="cpu",guitype='strbox', row=7, col=1, rowspan=1, colspan=1,mode="train")


	(options, args) = parser.parse_args()
	E2n=E2init(sys.argv,options.ppid)
	time00=time.time()
	
	if options.from_trained==None and options.nnet!=None:
		options.from_trained=options.nnet
	
	#### parse the options.
	nkernel=[int(i) for i in options.nkernel.split(',')]
	ksize=[int(i) for i in options.ksize.split(',')]
	poolsz=[int(i) for i in options.poolsz.split(',')]
	
	
	if options.applying:
		dirname="segmentations"
		try: os.mkdir(dirname)
		except: pass
	
		tomos=options.tomograms.split(',')
		
		for tm in tomos:
			print("Starting on {}...".format(tm))
			bn=base_name(tm)
			if options.outtag=="":	
				nn=options.nnet
				nn=nn[nn.rfind("__")+2:-4]
				outname="segmentations/{}__{}_seg.hdf".format(bn, nn)
			else:
				outname="segmentations/{}__{}.hdf".format(bn, options.outtag)
				
			print("Output segmentation will be written to {}...".format(outname))
			segout=apply_neuralnet(options, tm)
			segout.write_image(outname)
		print("Done.")
		print("Total time: ", time.time()-time00)
	
		E2end(E2n)
		return
	
	
	if options.trainset==None:
		print("No training set input...exit.")
		return
	
	
	if options.nettag=="":
		tag=options.trainset
		options.nettag=tag[tag.rfind('__')+2:-4].replace("_trainset","")
		
	if "CUDA_VISIBLE_DEVICES" in os.environ:
		print("CUDA_VISIBLE_DEVICES is already set as environment variable. This will overwrite the device option...")
	else:
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
	
	
	session=tf.Session()
	labelshrink=np.prod(poolsz)
	print("loading particles...")
	data, labels, shape, ntrain=load_particles(options.trainset,labelshrink,options.ncopy)
	batch_size=options.batch
	
	if options.from_trained!=None:
		convnet=StackedConvNet_tf.load_network(options.from_trained, session, imgsz=shape[0], bsz=batch_size)

	else:
		print("Setting up model...")
		kernels=[(k,ksize[i],poolsz[i]) for i,k in enumerate(nkernel)]
		convnet = StackedConvNet_tf(kernels, shape[0], batch_size, meanout=False)
		session.run(tf.global_variables_initializer())

	
	
	#convnet.save_network(options.netout, session, options)
	
	if (options.niter>0):	
		
		convnet.do_training(data, labels, session, shuffle=False, learnrate=options.learnrate, niter=options.niter)
		
		
		
	dirname="neuralnets"
	try: os.mkdir(dirname)
	except: pass
	options.netout="{}/nnet_save__{}.hdf".format(dirname, options.nettag)

	if options.trainout:
		fname="{}/trainout_nnet_save__{}.hdf".format(dirname, options.nettag)
		convnet.write_output_train(fname, session, writelabel=True)
		
	convnet.save_network(options.netout, session, options)
	print("Done")
	print("Total time: ", time.time()-time00)
	E2end(E2n)

def run(cmd):
	print(cmd)
	launch_childprocess(cmd)

	
def apply_neuralnet(options, tomogram=None):
	
	if tomogram==None:
		tomogram=options.tomograms
	tt0=time.time()
	
	nframe=EMUtil.get_image_count(tomogram)
	is3d=False
	### deal with 3D volume or image stack
	e=EMData(tomogram, 0, True)
	apix=e["apix_x"]
	if nframe==1:
		nframe=e["nz"]
		if nframe>1:
			#### input data is 3D volume
			#esz=max(e["nx"],e["ny"])
			is3d=True
			
	enx=e["nx"]
	eny=e["ny"]
	shape=[enx,eny]
	
	output=EMData(enx,eny, nframe)
	output["tomogram_src"]=tomogram
	output["nnet_src"]=options.from_trained
	output["apix_x"]=apix
	output["apix_y"]=apix
	output["apix_z"]=apix
	
	#####################
	print("Loading the Neural Net...")
	
	fname=options.from_trained
	hdr=EMData(fname,0)
	amplitude=hdr.get_attr_default("amplitude", 1.0)
	if amplitude<=0: amplitude=1.
	ksize=hdr["ksize"]
	poolsz=hdr["poolsz"]
	if hdr.has_attr("trained_from"):
		output["trainset"]=hdr["trained_from"]
	labelshrink=np.prod(poolsz)
	k=1
	#global layers
	layers=[]
	for i in range(len(ksize)):
		layer={}
		b=EMData(fname,k)
		s=b["w_shape"]
		k+=1
		allw=[]
		layer["b"]=b
		layer["shp"]=s
		layer["pool"]=poolsz[i]
		for wi in range(s[0]*s[1]):
			w=EMData(fname,k)
			allw.append(w)
			k+=1
			
		allw=np.asarray(allw).reshape((s[0],s[1]))
		for wi in range(s[0]):
			
			for mi in range(s[1]):
				sw=allw[wi][mi]["nx"]
				allw[wi][mi]=allw[wi][mi].get_clip(Region(((sw-enx)//2),((sw-eny)//2),enx,eny))
				
				allw[wi][mi].process_inplace("xform.phaseorigin.tocenter")
				#allw[wi][mi].do_fft_inplace()
				
		enx/=poolsz[i]
		eny/=poolsz[i]
		layer["allw"]=allw
		layers.append(layer)
	
	
	################
	enx=e["nx"]
	eny=e["ny"]
	
	print("Loading tomogram...")
	#global tomo_in
	tomo_in=[]
	for nf in range(nframe):
		if is3d:
			e0=EMData(tomogram, 0, False, Region(0,0,nf,enx,eny,1))
		else:
			e0=EMData(tomogram, nf, False, Region(0,0,enx,eny))
		tomo_in.append(e0)
	
	#####################
	print("Doing covolution...")
	
	try: os.remove(options.output)
	except: pass
	
	jobs=[]
	for nf in range(nframe):
		jobs.append((tomo_in, nf, layers))
		
		
	######### threading copied from e2spt_align.py
	jsd=Queue(0)
	NTHREADS=max(options.threads+1,2)
	thrds=[threading.Thread(target=do_convolve,args=(jsd,job)) for job in jobs]
	thrtolaunch=0
	while thrtolaunch<len(thrds) or threading.active_count()>1:
		#print thrtolaunch, len(thrds)
		if thrtolaunch<len(thrds) :
			 
			while (threading.active_count()==NTHREADS ) : time.sleep(.1)
			thrds[thrtolaunch].start()
			print("starting: ", thrtolaunch)#, e0["nx"]
			thrtolaunch+=1
		else: time.sleep(1)
	
		while not jsd.empty():
			idx,cout=jsd.get()
			#cout.write_image("tmp.hdf",-1)
			cout=cout.get_clip(Region(((cout["nx"]-enx)//2),((cout["ny"]-eny)//2) ,enx, eny))
			cout.scale(labelshrink)
			cout.div(amplitude)
			output.insert_clip(cout, [0,0,idx])
	return output
	
def do_convolve(jsd, job):
	tomo_in, idx, layers= job
	#idx=job
	
	e0=tomo_in[idx]
	e0.div(3.)
	imgs=[e0]
	
		
	for layer in layers:
		
		s0=imgs[0]["nx"]
		s1=imgs[0]["ny"]
		
		imgout=[]
		allw=layer["allw"]
		s=layer["shp"]
		poolsz=layer["pool"]
		b=layer["b"]
		#print s,poolsz,s0,s1
		for wi in range(s[0]):
			
			cv=EMData(imgs[0])
			cv.to_zero()
			for mi in range(s[1]):
				ww=allw[wi][mi]
				#print ww.numpy().shape
				cv.add(imgs[mi].process("math.convolution",{"with":ww}))
			
			if poolsz>1:
				cv=cv.process("math.maxshrink",{"n":poolsz})
			cv.add(b[wi])
			cv.process_inplace("threshold.belowtozero")
			imgout.append(cv)
		
		imgs=imgout
	#imgs[0].process_inplace("xform.phaseorigin.tocenter")

	jsd.put((idx,imgs[0]))
	#imgs[0].write_image(outname, idx)
	return# (idx,imgs[0])
	
	
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
			
			
			ar=EMNumPy.em2numpy(ptl)
			#shp=np.shape(ar)
			data.append(ar.flatten())
			
			ptl=EMData(ptcls,i*2+1)
			#ptl.process_inplace("threshold.belowtozero")
			if ncopy>1:
				ptl.process_inplace("xform",{"transform":tr})
			if labelshrink>1:
				ptl.process_inplace("math.meanshrink",{'n':labelshrink})
			ar=EMNumPy.em2numpy(ptl)
			#shp=np.shape(ar)
			label.append(ar.flatten())
	
	if ntrain<0: ntrain=len(data)
	## randomize
	rndid=list(range(ntrain))
	rng.shuffle(rndid)	
	rndid=rndid+list(range(ntrain, len(data)))
	data=[data[i] for i in rndid]
	label=[label[i] for i in rndid]
	
	print("{:d} particles loaded, {:d} in training set, {:d} in validation set".format(len(data), ntrain, len(data)-ntrain))
	data=np.asarray(data,dtype=np.float32)
	print("Std of particles: ",np.std(data.flatten()))
	#data/=np.std(data.flatten())*3  #np.max(np.abs(data))
	data/=3.
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
			
		outsz=(old_div(sz,np.prod([k[2] for k in kernels])))

		tf_data = tf.placeholder(tf.float32, shape=[None, sz*sz], name="tfdata")
		if meanout:
			tf_label = tf.placeholder(tf.float32, shape=[None], name="tflabel")
		else:
			tf_label = tf.placeholder(tf.float32, shape=[None, outsz*outsz], name="tflabel")

		#### since they are changing this after tf 1.5
		try:
			dataset = tf.data.Dataset.from_tensor_slices((tf_data, tf_label))
		except:
			dataset = tf.contrib.data.Dataset.from_tensor_slices((tf_data, tf_label))

		dataset_batch=dataset.batch(batchsz)
		iterator = dataset_batch.make_initializable_iterator()
		data_inp, data_tar=iterator.get_next()
		
		k0=1
		wts=[]
		bss=[]
		layerout=[tf.reshape(data_inp, (-1, sz, sz, 1), name="reshape_input")]
		for i in range(nlayer):
			x=layerout[-1]
			nk, ks, pl=kernels[i]
			w=tf_make_weight((ks,ks,k0,nk), pl)
			b=tf_make_bias([nk])
			convout=tf_conv2d(x, w)
			if pl==2:
				poolout=max_pool_2x2(convout)
			else:
				poolout=convout
				
			if i<nlayer-1:
				actout=tf.nn.relu(poolout+b)
			else:
		#		 actout=tf.minimum(1.,poolout+b)
				actout=poolout+b
			
			layerout.append(actout)
			k0=nk
			wts.append(w)
			bss.append(b)
			
		if meanout:
			yrl=y=tf.reshape(layerout[-1], (-1, outsz*outsz), name="reshape_output")
			yout=tf.maximum(-1.,(tf.reduce_mean(y, axis=1)))
			yout=tf.minimum(1., yout)
			loss=tf.reduce_mean((data_tar-yout)**2)
		else:
			y=tf.reshape(tf.minimum(1.,layerout[-1]), (-1, outsz*outsz), name="reshape_output")
			yrl=tf.reshape(tf.nn.relu(layerout[-1]), (-1, outsz*outsz), name="reshape_output_relu")
			loss=tf.log(tf.reduce_mean((data_tar-y)**2))
			
			
		self.meanout=meanout
		self.dataset=dataset
		self.kernels=kernels
		self.batchsize=batchsz
		self.wts=wts
		self.bss=bss
		self.outsz=outsz
		self.imgsz=sz
		self.layerout=layerout
		self.data=tf_data
		self.label=tf_label
		
		self.batch_in=data_inp
		self.batch_out=y
		self.batch_out_rl=yrl
		self.batch_tar=data_tar
		
		self.loss=loss
		self.iterator=iterator
		self.amplitude=1.
		

	def do_training(self, data, label, session, shuffle=False, learnrate=1e-4, niter=10):
		print("Training...")
		optimizer=tf.train.AdamOptimizer(learnrate)
		train_step=optimizer.minimize(self.loss)
		session.run(tf.variables_initializer(optimizer.variables()))
		
		for it in range(niter):
			if shuffle:
				self.dataset.shuffle(10000)
			session.run(self.iterator.initializer,
				feed_dict={self.data: data, self.label: label})
			cost=[]
			while True:
				try:
					cost.append(session.run((train_step, self.loss))[1])
				except tf.errors.OutOfRangeError:
					break
					
			print("iteration {}, cost {:.3f}".format(it, np.mean(cost)))
		
		session.run(self.iterator.initializer, feed_dict={self.data: data, self.label: label})
		
		
	def save_network(self, fname, session, options=None):
		try: os.remove(fname)
		except: pass
		print("Saving the trained net to {}...".format(fname))
		weights=session.run(self.wts)
		bias=session.run(self.bss)
		wsz=int(weights[0].shape[0])

		hdr=EMData(wsz,wsz)
		hdr["nkernel"]=[k[0] for k in self.kernels]
		hdr["ksize"]=[k[1] for k in self.kernels]
		hdr["poolsz"]=[k[2] for k in self.kernels]
		hdr["imageshape"]=(self.batchsize,1,self.imgsz,self.imgsz)
		hdr["amplitude"]=self.amplitude
		hdr["meanout"]=self.meanout
		if options!=None:
			hdr["trained_from"]=options.trainset
		nlayer=len(self.kernels)



		hdr.write_image(fname,0)

		k=1
		for i in range(nlayer):
			w=weights[i].copy()
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
	def load_network(fname, session,  imgsz=-1, bsz=0):
		
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
		
		meanout=hdr.get_attr_default("meanout", True)
		nnet=StackedConvNet_tf(kernels, sz, batchsize, meanout)
		nnet.amplitude=hdr.get_attr_default("amplitude", 1.0)
		k=1
		for i in range(len(nnet.kernels)):
			e=EMData(fname,k)
			s=e["w_shape"]
			b=e.numpy().copy().astype(np.float32)
			k+=1
			session.run(tf.assign(nnet.bss[i], b))
			ks=ksize[i]
			allw=np.zeros((s[0]*s[1], ks, ks))
			for wi in range(s[0]*s[1]):
				e=EMData(fname,k)
				sw=e["nx"]
				#e=e.get_clip(Region((sw-ks)//2,(sw-ks)//2,ks,ks))
				k+=1
				w=e.numpy().copy().astype(np.float32)
				allw[wi]=w
			allw=allw.reshape([s[0], s[1], ks, ks]).astype(np.float32).transpose(3,2,1,0)
			session.run(tf.assign(nnet.wts[i], allw))
			
		return nnet	
		
	def write_output_train(self, outfile, session, ncopy=10, writelabel=False):
		sz=self.imgsz
		outsz=self.outsz
		try: os.remove(outfile)
		except: pass
		
		print("Writting network output of training set to {}...".format(outfile))
	
		amp=[]
		for nc in range(ncopy):
			if writelabel:
				outx, outtar, outy =session.run((self.batch_in,self.batch_tar,self.batch_out_rl))
			else:
				outx, outy =session.run((self.batch_in,self.batch_out_rl))
			for i in range(len(outx)):
				ox=outx[i].reshape((sz,sz))
				
				oy=outy[i].reshape((outsz, outsz))
				
				

				e0=from_numpy(ox.copy())
				e0.process_inplace("normalize")
				e0.write_image(outfile, -1)
				
				if writelabel:
					ot=outtar[i].reshape((outsz, outsz))
					lb=ot>.7
					if np.sum(lb)>0:
						amp.append(np.mean(oy[lb]))
					e1=from_numpy(ot.copy())
					e1=e1.get_clip(Region(old_div(-(sz-outsz),2),old_div(-(sz-outsz),2),sz,sz))
					e1.scale(old_div(sz,outsz))

					e1.write_image(outfile, -1)
					

				e1=from_numpy(oy.copy())
				e1=e1.get_clip(Region(old_div(-(sz-outsz),2),old_div(-(sz-outsz),2),sz,sz))
				e1.scale(old_div(sz,outsz))

				e1.write_image(outfile, -1)

		if len(amp)>0:
			#print(np.mean(amp), amp)
			self.amplitude=float(np.mean(amp))
	
def tf_make_weight(shape, pz=1, small=.01):
	initial = tf.truncated_normal(shape, stddev=small)
	return tf.Variable(initial)

def tf_make_bias(shape, small=0.0):
	initial = tf.constant(small, shape=shape)
	return tf.Variable(initial)

def tf_conv2d(x, W):
	return tf.nn.conv2d(x, W, strides=[1, 1, 1, 1], padding='SAME')

def max_pool_2x2(x):
	return tf.nn.max_pool(x, ksize=[1, 2, 2, 1], strides=[1, 2, 2, 1], padding='SAME')

	
	

####################### theano version






#class StackedConvNet(object):
	
	#def __init__(self,rng,nkernel,ksize,poolsz,imageshape):
		#self.amplitude =1.
		#self.n_kernel=nkernel
		#self.ksize=ksize
		#self.n_convlayers=len(self.ksize)
		#self.x = T.matrix(name='input')
		#self.image_shape=imageshape
		#input_shape=self.image_shape
		#print("Shape of neural networl input: ",input_shape)
		#self.convlayers=[]
		#self.params=[]
		#if poolsz:
			#self.poolsz=poolsz
		#else:
			#self.poolsz=[2]*self.n_convlayers
			
		##self.poolsz[-1]=1
		#poolsz=self.poolsz
		#convin=self.x
		#self.labelshrink=1
		#for i in range(self.n_convlayers):
			#pz=poolsz[i]
			#if pz<0:
				#pz=1.0/abs(pz)
			#convlayer = LeNetConvPoolLayer(
				#rng,
				#image_shape=input_shape,
				#filter_shape=(self.n_kernel[i], input_shape[1], self.ksize[i], self.ksize[i]),
				#poolsize=poolsz[i],
				#xin=convin
			#)
			#self.convlayers.append(convlayer)
			##self.weights.append(convlayer.W)
			
			#self.labelshrink=int(self.labelshrink*pz)#poolsz[i]
			#input_shape=(input_shape[0],self.n_kernel[i],input_shape[2]/pz,input_shape[3]/pz)
			#convin=convlayer.hidden
			
			#self.params.extend(convlayer.params)
		
		
		
		#self.clslayer=self.convlayers[-1]
		#self.outsize=int(input_shape[2])
		##print 
		##self.labelshrink=2**(self.n_convlayers-1)
		
	#def get_pretrain_func(self,data,batch_size):
		#learning_rate = T.scalar('lr')  # learning rate to use
		#weight_decay = T.scalar('wd')  # learning rate to use
		#index = T.lscalar()	# index to a [mini]batch
		
		#cost, updates = self.convlayer.get_cost_updates(learning_rate,weight_decay)
		#train_model = theano.function(
					#inputs=[
						#index,
						#theano.In(learning_rate, value=0.1),
						#theano.In(weight_decay, value=1e-5)
					#],
					#outputs=cost,
					#updates=updates,
					#givens={
						#self.convlayer.x: data[index * batch_size: (index+1) * batch_size]
					#}
				#)
		#return train_model
					
	#def get_classify_func(self,data,lab,batch_size):
		#learning_rate = T.scalar('lr')  # learning rate to use
		#weight_decay = T.scalar('wd')  # learning rate to use
		#index = T.lscalar()	# index to a [mini]batch
		
		#label=T.matrix(name='label')
		
		#cost = self.clslayer.get_cost_hidden(label)
		#for c in self.convlayers:
			#cost+=weight_decay*T.sum(abs(c.W))
		#gparams = T.grad(cost, self.params)
		#updates = [
			#(param, param - learning_rate * gparam)
			#for param, gparam in zip(self.params, gparams)
		#]
		##cost, updates = self.clslayer.get_cost_updates_hidden(cls,learning_rate,weight_decay)
		#train_classify = theano.function(
					#inputs=[
						#index,
						#theano.In(learning_rate, value=0.1),
						#theano.In(weight_decay, value=1e-5)
					#],
					#outputs=cost,
					#updates=updates,
					#givens={
						#self.x: data[index * batch_size: (index+1) * batch_size],
						#label: lab[index * batch_size: (index+1) * batch_size]
					#}
					#)
		#return train_classify
	
	#def update_shape(self, imgshape):
		#input_shape=imgshape
		#if isinstance(self.poolsz,list):
			#poolsz=self.poolsz
		#else:
			#poolsz=[2]*(self.n_convlayers-1)+[1]
		#for i in range(self.n_convlayers):
			#pz=poolsz[i]
			#if pz<0:
				#pz=1.0/abs(pz)
			
			#self.convlayers[i].image_shape.set_value(input_shape, borrow=True)		
			#input_shape=(input_shape[0],self.n_kernel[i],input_shape[2]/pz,input_shape[3]/pz)
			
		#self.outsize=int(input_shape[2])
		
		
#class LeNetConvPoolLayer(object):

	#def __init__(self, rng, filter_shape, image_shape, poolsize=2, xin=None):
		
		#assert image_shape[1] == filter_shape[1]
		#self.image_shape=theano.shared(
			#value=np.asarray(image_shape,dtype='int16'),borrow=True)
		#self.poolsize=(poolsize,poolsize)
		##self.input = input
		#if xin:
			#self.x=xin
		#else:
			#self.x = T.matrix(name='input')
		#self.x1=self.x.reshape(self.image_shape,ndim=4)
		#self.filter_shape=filter_shape
		
		## there are "num input feature maps * filter height * filter width"
		## inputs to each hidden unit
		#fan_in = np.prod(filter_shape[1:])
		## each unit in the lower layer receives a gradient from:
		## "num output feature maps * filter height * filter width" /
		##   pooling size
		#fan_out = (filter_shape[0] * np.prod(filter_shape[2:]) /
			#np.prod(self.poolsize))
		## initialize weights with random weights
		#W_bound = np.sqrt(6. / (fan_in + fan_out))
		#self.W = theano.shared(
			#np.asarray(
				#rng.uniform(low=-W_bound, high=W_bound, size=filter_shape),
				#dtype=theano.config.floatX
			#),
			#borrow=True
		#)
		#self.W_prime=self.W[:,:,::-1,::-1]
		#self.W_prime=self.W_prime.dimshuffle(1,0,2,3)
		##self.W_prime=self.W_prime[:,::-1]
		##print self.W.get_value()
		##print self.W_prime.eval()
		## the bias is a 1D tensor -- one bias per output feature map
		#b_values = np.zeros((filter_shape[0],), dtype=theano.config.floatX)
		#bp_values = np.zeros((filter_shape[1],), dtype=theano.config.floatX)
		#self.b = theano.shared(value=b_values, borrow=True)
		#self.b_prime = theano.shared(value=bp_values, borrow=True)
		
		#if poolsize<-1:
			#self.x1=self.x1.repeat(int(-poolsize), axis=2).repeat(int(-poolsize), axis=3)

		## convolve input feature maps with filters
		#conv_out = conv2d(
			#input=self.x1,
			#filters=self.W,
			#filter_shape=filter_shape,
			##image_shape=self.image_shape.eval(),
			#border_mode='full'
		#)
		#bp=(filter_shape[2]-1)/2
		
		#conv_out=conv_out[:,:,bp:-bp,bp:-bp]
		
		## downsample each feature map individually, using maxpooling
		#if poolsize>1:
			#try:
				#self.pooled_out = pool.pool_2d(
					#input=conv_out,
					#ws=self.poolsize,
					#ignore_border=True
				#)
			#except:
				
				#self.pooled_out = pool.pool_2d(
					#input=conv_out,
					#ds=self.poolsize,
					#ignore_border=True
				#)
		#else:
			#self.pooled_out=conv_out
		
		#self.hidden = T.maximum(0,(self.pooled_out + self.b.dimshuffle('x', 0, 'x', 'x')))

		## store parameters of this layer
		#self.params = [self.W,self.b]
		##self.params = [self.W]
	
	
	#def get_reconstructed_input(self):
		#""" Computes the reconstructed input given the values of the hidden layer """

		#repeated_conv = conv2d(input = self.hidden,
			      #filters = self.W_prime,
			      #border_mode='full')
		##repeated_conv=repeated_conv[:,:,1:-1,1:-1]
		#bp=(self.filter_shape[2]-1)/2
		#repeated_conv=repeated_conv[:,:,bp:-bp,bp:-bp]
		
		#multiple_conv_out = [repeated_conv.flatten()] * np.prod(self.poolsize)
		#stacked_conv_neibs = T.stack(*multiple_conv_out).T

		#stretch_unpooling_out =  T.nnet.neighbours.neibs2images(stacked_conv_neibs, 
							 #self.poolsize, self.x1.shape)
		##return self.hidden
		#z=T.tanh(stretch_unpooling_out + self.b_prime.dimshuffle('x', 0, 'x', 'x'))
		##return T.sum(T.sum((self.x1-z)**2,axis=3),axis=2)
		##rectified_linear_activation = lambda x: T.maximum(0.0, x)
		#return z

	#def get_image(self, relu=True):
		##return T.tanh(self.hidden)*2-1
		##return T.maximum(0,(self.pooled_out + self.b.dimshuffle('x', 0, 'x', 'x')))
		#if relu:
			#return T.minimum(1,(self.pooled_out + self.b.dimshuffle('x', 0, 'x', 'x')))
		#else:
			#return self.pooled_out + self.b.dimshuffle('x', 0, 'x', 'x')
	
	#def get_cost_hidden(self,  label):
		
		
		##z = self.hidden
		##shp=z.shape
		##z = z.reshape((shp[0],shp[2]*shp[3]))
		##z = T.mean(z,axis=1)
		##x=label-.5
		####cost=T.mean((z-x)**2)
		###cost=1-T.mean(z*x)
		##cost=1-T.mean(z*x)
		##cost+= T.sum(abs(self.W))*weight_decay
		
		
		
		#xin=label
		##xin=xin.reshape(self.hidden.shape,ndim=4)
		#z = self.get_image()
		
		#shp=z.shape
		
		#xin=xin.reshape((shp[0],-1))
		#z=z.reshape((shp[0],-1))
		
		##L = T.sum(xin*z,axis=1)/T.sqrt(T.sum(xin**2,axis=1))/T.sqrt(T.sum(z**2,axis=1))
		##cost = 1-T.mean(L)
		
		##L = T.mean(xin*z,axis=1)
		##cost = 1-T.mean(L)
		
		##cost=T.log(T.mean(((xin-z)**2)*T.exp(0.7*(xin-z)) ))
		#cost=T.log(T.mean(((xin-z)**2) ))
		##cost+= T.sum(abs(self.W))*weight_decay
		

		#return cost
	
	#def get_cost_updates(self, learning_rate,weight_decay,xin=None):
		
			#""" This function computes the cost and the updates for one trainng
			#step of the dA """
			#if xin==None:
				#xin=self.x1
				
			#z = self.get_reconstructed_input()
			
			#shp=self.image_shape
			#xin=xin.reshape((shp[0],shp[2]*shp[3]))
			#z=z.reshape((shp[0],shp[2]*shp[3]))
			##L = T.sum(T.sum((xin-z)**2,axis=3),axis=2)
			#L = T.sum(xin*z,axis=1)/T.sqrt(T.sum(xin**2,axis=1))/T.sqrt(T.sum(z**2,axis=1))
			#cost = 1-T.mean(L)
			##L = T.sum(xin*z,axis=1)/T.sqrt(T.sum(xin**2,axis=1))/T.sqrt(T.sum(z**2,axis=1)
			## note : L is now a vector, where each element is the
			##	cross-entropy cost of the reconstruction of the
			##	corresponding example of the minibatch. We need to
			##	compute the average of all these to get the cost of
			##	the minibatch
			##cost = T.mean(L)
			#cost+= T.sum(abs(self.W))*weight_decay
			## compute the gradients of the cost of the `dA` with respect
			## to its parameters
			#gparams = T.grad(cost, self.params)
			## generate the list of updates
			#updates = [
				#(param, param - learning_rate * gparam)
				#for param, gparam in zip(self.params, gparams)
			#]

			#return (cost, updates)

if __name__ == '__main__':
    main()
