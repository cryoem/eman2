#!/usr/bin/env python
# Muyuan July 2015
import sys
import random
import numpy as np
from EMAN2 import *
import cPickle
import time
import threading
from Queue import Queue
from multiprocessing import Array

def import_theano():
	global theano,T,conv,pool
	import theano
	import theano.tensor as T
	from theano.tensor.nnet import conv
	from theano.tensor.signal import pool


def main():

	usage="""Segment a tomograph using convolutional neural network. This program is still experimental. Please consult the developers before using."""
	print usage
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_header(name="tmpheader", help='temp label', title="### This program is NOT avaliable yet... ###", row=0, col=0, rowspan=1, colspan=2, mode="train,test")
	parser.add_argument("--trainset",help="Training set.", default=None, guitype='filebox', browser="EMParticlesTable(withmodal=True)",  row=1, col=0,rowspan=1, colspan=3, mode="train")
	parser.add_argument("--from_trained", type=str,help="Start from pre-trained neural network", default=None,guitype='filebox',browser="EMBrowserWidget(withmodal=True)", row=2, col=0, rowspan=1, colspan=3, mode="train,test")
	parser.add_argument("--netout", type=str,help="Output neural net file name", default="nnet_save.hdf",guitype='strbox', row=3, col=0, rowspan=1, colspan=3, mode="train")
	
	parser.add_argument("--learnrate", type=float,help="Learning rate ", default=.01, guitype='floatbox', row=4, col=0, rowspan=1, colspan=1, mode="train")
	parser.add_argument("--niter", type=int,help="Training iterations", default=20, guitype='intbox', row=4, col=1, rowspan=1, colspan=1, mode="train")
	parser.add_argument("--ncopy", type=int,help="Number of copies for each particle", default=1, guitype='intbox', row=5, col=0, rowspan=1, colspan=1, mode="train")
	parser.add_argument("--batch", type=int,help="Batch size for the stochastic gradient descent. Default is 20.", default=20, guitype='intbox', row=5, col=1, rowspan=1, colspan=1, mode="train")
	parser.add_argument("--nkernel", type=str,help="Number of kernels for each layer, from input to output. The number of kernels in the last layer must be 1. ", default="40,40,1", guitype='strbox', row=6, col=0, rowspan=1, colspan=1, mode="train")
	parser.add_argument("--ksize", type=str,help="Width of kernels of each layer, the numbers must be odd. Note the number of layers should be the same as the nkernel option. ", default="15,15,15", guitype='strbox', row=6, col=1, rowspan=1, colspan=1, mode="train")
	parser.add_argument("--poolsz", type=str,help="Pooling size for each layer. Note the number of layers should be the same as the nkernel option. ", default="2,1,1", guitype='strbox', row=7, col=0, rowspan=1, colspan=1, mode="train")
	parser.add_argument("--weightdecay", type=float,help="Weight decay. Used for regularization.", default=1e-6, guitype='floatbox', row=7, col=1, rowspan=1, colspan=1, mode="train")
	parser.add_argument("--trainout", action="store_true", default=False ,help="Output the result of the training set", guitype='boolbox', row=8, col=0, rowspan=1, colspan=1, mode='train[True]')
	parser.add_argument("--training", action="store_true", default=False ,help="Doing training", guitype='boolbox', row=8, col=1, rowspan=1, colspan=1, mode='train[True]')
	parser.add_argument("--tomograms", type=str,help="Tomograms input.", default=None,guitype='filebox',browser="EMBrowserWidget(withmodal=True)", row=1, col=0, rowspan=1, colspan=3, mode="test")
	parser.add_argument("--applying", action="store_true", default=False ,help="Applying the neural network on tomograms", guitype='boolbox', row=4, col=0, rowspan=1, colspan=1, mode='test[True]')
	parser.add_argument("--dream", action="store_true", default=False ,help="Iterativly applying the neural network on noise")
	parser.add_argument("--to3d", action="store_true", default=True ,help="convert to result to 3D.", guitype='boolbox', row=5, col=1, rowspan=1, colspan=1, mode='test')
	parser.add_argument("--output", type=str,help="Segmentation out file name", default="tomosegresult.hdf", guitype='strbox', row=3, col=0, rowspan=1, colspan=1, mode="test")
	parser.add_argument("--threads", type=int,help="Number of thread to use when applying neural net on test images. Not used during trainning", default=12, guitype='intbox', row=10, col=0, rowspan=1, colspan=1, mode="test")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)

	(options, args) = parser.parse_args()
	E2n=E2init(sys.argv,options.ppid)
	
	#### parse the options.
	options.nkernel=[int(i) for i in options.nkernel.split(',')]
	options.ksize=[int(i) for i in options.ksize.split(',')]
	if options.poolsz:
		options.poolsz=[int(i) for i in options.poolsz.split(',')]
	
	#### This is supposed to test the overfitting of the network by applying it on pure noise repeatly
	#### The function is no longer maintained so it may or may not work..
	if options.dream:
		print "This function is no longer supported.. exit."
		return
		#os.environ["THEANO_FLAGS"]="optimizer=None"
		#print "Testing on big images, Theano optimizer disabled"
		#import_theano()
		#convnet=load_model(options.from_trained)
		#dream(convnet,options)
		#E2end(E2n)
		#exit()
		
	
	if options.applying:
		apply_neuralnet(options)
		E2end(E2n)
		exit()
	
	
	os.environ["THEANO_FLAGS"]="optimizer=fast_run"
	import_theano()

	
	batch_size=options.batch
	#### Train da with particles first.
	
	if options.trainset==None:
		print "No training set input...exit."
		exit()
	
	
	rng = np.random.RandomState(123)

	labelshrink=np.prod(options.poolsz)
	print "loading particles..."
	particles=load_particles(options.trainset,labelshrink,options.ncopy, rng)

	train_set_x= particles[0]
	labels=particles[1]
	shape=particles[2]
	ntrain=particles[3]
	#print "Number of particles: {}".format(train_set_x.shape.eval()[0])
	
	# allocate symbolic variables for the data
	index = T.lscalar()	# index to a [mini]batch
	x = T.matrix('x')  # the data is presented as rasterized images
	image_shape=(batch_size, shape[2], shape[0],shape[1])
	
	if options.from_trained!=None:
		convnet=load_model(options.from_trained)
		convnet.update_shape(image_shape)
	else:
		print "setting up model"
		convnet = StackedConvNet(
			rng,
			nkernel=options.nkernel,
			ksize=options.ksize,
			poolsz=options.poolsz,
			imageshape=image_shape
		)
	
	
	
	#print shape
	
	
	if (options.niter>0):	
		print "training the convolutional network..."
		
		classify=convnet.get_classify_func(train_set_x,labels,batch_size)
			
		learning_rate=options.learnrate
		n_train_batches = train_set_x.get_value(borrow=True).shape[0] / batch_size
		v0=np.inf
		nbad=0
		for epoch in xrange(options.niter):
		# go through the training set
			
			c = []   #### train set loss
			v = []   #### valid set loss
			if epoch==0:
				print classify(0,lr=learning_rate,wd=options.weightdecay)
			for batch_index in xrange(n_train_batches):
				if batch_index*batch_size < ntrain:
					err=classify(batch_index,
						lr=learning_rate,
						wd=options.weightdecay)
					c.append(err)
					if epoch==0 and batch_index<5:
						print err
				else:
					err=classify(batch_index,
						lr=0,
						wd=options.weightdecay)
					v.append(err)
					
			#print len(v), len(c)
			learning_rate*=.9
			print "Training epoch {:d}, train loss {:.3f}, learning rate {:.3f}".format(epoch, np.mean(c),  learning_rate),
			if len(v)>0:
				print "valid loss {:.3f}".format(np.mean(v)), 
				if np.mean(v)>v0 and np.mean(v)>np.mean(c):
					nbad+=1
					print '*'
				else:
					nbad=0
					print 
				v0=np.mean(v)
				if nbad>2:
					print "loss increase in validation set. Overfitting. Stop."
					break
			else:
				print 

		
		
		
	#######################################
	#print convnet.clslayer.W.get_value()
	#print convnet.clslayer.b.get_value()
			
	if options.trainout:
		print "Generating results ..."
		nsample=100
		convnet.update_shape((nsample, shape[2], shape[0],shape[1]))
		test_cls = theano.function(
			inputs=[],
			outputs=convnet.clslayer.get_image(False),
			givens={
				convnet.x: train_set_x[:nsample]
			}
		)
		if options.netout.endswith(".hdf"):
			fname="trainout_{}".format(options.netout)
		else:
			fname="trainout_{}.hdf".format(options.netout)
		try:os.remove(fname)
		except: pass
		#print convnet.outsize,shape
		mid=test_cls()
		
		ipt= train_set_x[:nsample]
		ipt= ipt.eval()
		
		lb= labels[:nsample].eval()
		amp=[]
		for t in range(nsample):
			
			#### raw image
			if shape[2]==1:
				img=ipt[t].reshape(shape[0],shape[1])
			else:
				img=ipt[t].reshape(shape[2],shape[0],shape[1])
				img=img[shape[2]/2]
			e0 = from_numpy(img.astype("float32"))
			e0.write_image(fname,-1)
			
			#### manual annotation
			img=lb[t].reshape(convnet.outsize,convnet.outsize)
			e1 = from_numpy(img.astype("float32"))
			e1=e1.get_clip(Region((convnet.outsize-shape[0])/2,(convnet.outsize-shape[0])/2,shape[0],shape[0]))
			e1.scale(float(shape[0])/float(convnet.outsize))
			e1.process_inplace("threshold.binary", {"value":.67})
			e1.write_image(fname,-1)
			
			#### neural net output
			img=mid[t].reshape(convnet.outsize,convnet.outsize)
			e2 = from_numpy(img.astype("float32"))
			e2=e2.get_clip(Region((convnet.outsize-shape[0])/2,(convnet.outsize-shape[0])/2,shape[0],shape[0]))
			#print float(shape[0])/float(convnet.outsize)
			e2.scale(float(shape[0])/float(convnet.outsize))
			e2.write_image(fname,-1)
			
			#### measure the amplitude of the neural network output by comparing it to the label
			e2.mult(e1)
			amp.append(e2["mean_nonzero"])
		print "amplitude: ", np.mean(amp)
		convnet.amplitude=np.mean(amp)
		print "Writing output on training set in {}".format(fname)
		
	save_model(convnet, options.netout, options)
	
	print "Done"
	E2end(E2n)

def run(cmd):
	print cmd
	launch_childprocess(cmd)

def save_model(convnet, fname, options=None):
	print "Saving the trained net to {}...".format(fname)
	#fname="nnet_save.hdf"
	sz=int(convnet.convlayers[0].W.shape[-1].eval())

	hdr=EMData(sz,sz)
	hdr["nkernel"]=convnet.n_kernel
	hdr["ksize"]=convnet.ksize
	hdr["poolsz"]=convnet.poolsz
	hdr["imageshape"]=convnet.image_shape
	hdr["amplitude"]=convnet.amplitude
	if options:
		if options.trainset:
			hdr["trainset_src"]=options.trainset
	
	
	hdr.write_image(fname,0)

	k=1
	for i in range(convnet.n_convlayers):
		w=convnet.convlayers[i].W.get_value()        
		b=convnet.convlayers[i].b.get_value()
		s=w.shape
		
		e=from_numpy(b)
		e["w_shape"]=s
		e.write_image(fname,k)
		k+=1
		w=w.reshape(s[0]*s[1], s[2], s[3])
		for wi in w:
			e=from_numpy(wi)
			e.write_image(fname,k)
			k+=1


	
def load_model(fname):
	print "loading model from {}...".format(fname)
	try:
		f = file(fname, 'rb')
		convnet = cPickle.load(f)
		f.close()
		return convnet
	except:
		print "Reading weight matrix from hdf file..."
		
	hdr=EMData(fname,0)
	
	nkernel=hdr["nkernel"]
	ksize=hdr["ksize"]
	poolsz=hdr["poolsz"]
	imageshape=hdr["imageshape"]
	rng = np.random.RandomState(123)
	
	savenet= StackedConvNet(
		rng,
		nkernel=nkernel,
		ksize=ksize,
		poolsz=poolsz,
		imageshape=imageshape
	)
	k=1
	for i in range(savenet.n_convlayers):
		e=EMData(fname,k)
		s=e["w_shape"]
		b=e.numpy().copy().astype(theano.config.floatX)
		k+=1
		savenet.convlayers[i].b.set_value(b)
		allw=np.zeros((s[0]*s[1], s[2],s[3]))
		for wi in range(s[0]*s[1]):
			e=EMData(fname,k)
			k+=1
			w=e.numpy()
			allw[wi]=w.copy()
		allw=allw.reshape(s).astype(theano.config.floatX)
		savenet.convlayers[i].W.set_value(allw)
		#print w.shape, b.shape
	
	savenet.amplitude=hdr.get_attr_default("amplitude", 1.0)
	return savenet

#def dream(convnet,options):
	#convz=convnet.image_shape[1]
	
	#print "Dreaming....."
	#try: 
		#os.remove(options.output)
		#print "Overwriting the output.."
	#except: pass
	
	#sz=500
	#shape=[sz,sz,convz]
	#img=np.random.randn(shape[2],shape[0],shape[1])
	#convnet.update_shape((1, shape[2], shape[0],shape[1]))
	#e=from_numpy(img)
	#e.write_image(options.output,-1)
	
	#for it in range(50):
		
		
		#print "Iteration {}...".format(it)
		
		
		##### prepare inputs
		#img[np.abs(img)>1.0]=1.0
		##print img.shape
		##print convz,sz
		#ar=img.reshape((convz,sz*sz))
		
		#data=theano.shared(np.asarray(ar,dtype=theano.config.floatX),borrow=True)

		##### write input when testing...
		##### 
		##print np.shape(data.get_value())
		##print "Applying the convolution net..."
		#test_imgs = theano.function(
			#inputs=[],
			#outputs=convnet.clslayer.get_image(),
			#givens={
				#convnet.x: data
			#}
		#)
			
		#img=test_imgs()
		##print np.shape(img)
		#img=img.reshape(convnet.outsize,convnet.outsize).T
		
		#e = from_numpy(img)
		#e.process_inplace("math.fft.resample",{"n":float(1./convnet.labelshrink)})
		#e.process_inplace("normalize")
		#e.write_image(options.output,-1)
		#if convz>1:
			#img=np.random.randn(shape[2],shape[0],shape[1])
			#img[convz/2]=e.numpy().copy()
		#else:
			#img=e.numpy().copy()
	
	#print "Output written to {}.".format(options.output)
	
	
def apply_neuralnet(options):
	
	tt0=time.time()
	
	nframe=EMUtil.get_image_count(options.tomograms)
	is3d=False
	### deal with 3D volume or image stack
	e=EMData(options.tomograms, 0, True)
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
	output["tomogram_src"]=options.tomograms
	output["nnet_src"]=options.from_trained
	output["apix_x"]=apix
	output["apix_y"]=apix
	output["apix_z"]=apix
	
	#####################
	print "Loading the Neural Net..."
	
	fname=options.from_trained
	hdr=EMData(fname,0)
	amplitude=hdr.get_attr_default("amplitude", 1.0)
	if amplitude<=0: amplitude=1.
	ksize=hdr["ksize"]
	poolsz=hdr["poolsz"]
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
				allw[wi][mi]=allw[wi][mi].get_clip(Region((sw-enx)/2,(sw-eny)/2,enx,eny))
				
				allw[wi][mi].process_inplace("xform.phaseorigin.tocenter")
				#allw[wi][mi].do_fft_inplace()
				
		enx/=poolsz[i]
		eny/=poolsz[i]
		layer["allw"]=allw
		layers.append(layer)
	
	
	################
	enx=e["nx"]
	eny=e["ny"]
	
	print "Loading tomogram..."
	#global tomo_in
	tomo_in=[]
	for nf in range(nframe):
		if is3d:
			e0=EMData(options.tomograms, 0, False, Region(0,0,nf,enx,eny,1))
		else:
			e0=EMData(options.tomograms, nf, False, Region(0,0,enx,eny))
		tomo_in.append(e0)
	
	#####################
	print "Doing covolution..."
	
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
			print "starting: ", thrtolaunch#, e0["nx"]
			thrtolaunch+=1
		else: time.sleep(1)
	
		while not jsd.empty():
			idx,cout=jsd.get()
			cout=cout.get_clip(Region((cout["nx"]-enx)/2,(cout["ny"]-eny)/2 ,enx, eny))
			cout.scale(labelshrink)
			cout.div(amplitude)
			output.insert_clip(cout, [0,0,idx])
			
	#if is3d and options.to3d:
		#to3d=" --twod2threed"
	#else:
		#to3d=""
		
		
	#fout="__tomoseg_tmp.hdf"
	#run("e2proc2d.py {} {} --process math.fft.resample:n={} {} --apix {} ".format(options.output,fout,float(1./labelshrink), to3d, apix))
	##os.rename(fout, options.output)
	output.write_image(options.output)
	print "Done."
	print "Total time: ", time.time()-tt0
	
def do_convolve(jsd, job):
	tomo_in, idx, layers= job
	#idx=job
	
	e0=tomo_in[idx]
	e0.div(3.)
	imgs=[e0]
	
	#if len(rg)>4:
		#e0=EMData(tomo, 0, False, Region(rg[0],rg[1],rg[2],rg[3],rg[4],rg[5]))
	#else:
		#e0=EMData(tomo, idx, False, Region(rg[0],rg[1],rg[2],rg[3]))
	
		
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
	
#def do_saveimg(job):
	#job[1].write_image(

#def apply_neuralnet_theano(options):
	
	#tt0=time.time()
	#import_theano()
	
	#convnet=load_model(options.from_trained)
	
	#convz=convnet.image_shape[1]
	
	#try: 
		#os.remove(options.output)
		#print "Overwriting the output.."
	#except: pass
	
	
	#nframe=EMUtil.get_image_count(options.tomograms)
	#is3d=False
	#### deal with 3D volume or image stack
	#e=EMData(options.tomograms, 0, True)
	#if nframe==1:
		#nframe=e["nz"]
		#if nframe>1:
			##### input data is 3D volume
			##esz=max(e["nx"],e["ny"])
			#is3d=True
			
	#esz=max(e["nx"],e["ny"])
	#enx=e["nx"]
	#eny=e["ny"]
	#shape=[esz,esz,convz]
	#convnet.update_shape((1, shape[2], shape[0],shape[1]))
	#for nf in range(nframe):
		#if is3d:
			#e=EMData(options.tomograms,0,False,Region(0,0,nf,esz,esz,convz))
		#else:
			#e=EMData(options.tomograms,nf, False, Region(0,0,esz,esz))
		
		#print "Applying the convolution net on image No {}...".format(nf)

		##### prepare inputs
		##e.process_inplace("normalize")
		#enp=e.numpy()
		##print enp.shape
		##exit()
		#enp[enp>3.0]=3.0
		#enp/=3
		#ar=enp.reshape((convz,shape[0]*shape[1]))
		
		
		#data=theano.shared(np.asarray(ar,dtype=theano.config.floatX),borrow=True)
		
	
		##### input is one 2D slice: just test the performance.
		##### write input when testing...
		#if nframe==1:
			#img=data[0].eval().reshape(shape[0],shape[1]).T
			#e = EMNumPy.numpy2em(img.astype("float32"))
			##e.process_inplace("normalize")
			#newshp=convnet.outsize
			#scl=float(convnet.outsize)/float(esz)
			#e.scale(scl)
			#e.clip_inplace(Region((esz-newshp)/2,(esz-newshp)/2,enx*scl,eny*scl))
			##e=e.get_clip(Region((shape[0]-newshp)/2,(shape[0]-newshp)/2,newshp,newshp))
			#e.process_inplace("normalize")
			#e.write_image(options.output,-1)
			##print (esz-newshp)/2,enx,eny, float(convnet.outsize)/float(esz), img.shape
		##exit()
		##### Applying the convolution net...
		#test_imgs = theano.function(
			#inputs=[],
			#outputs=convnet.clslayer.get_image(),
			#givens={
				#convnet.x: data
			#}
		#)
		
		#img=test_imgs()
		##print np.shape(img)
		#img=img.reshape(convnet.outsize,convnet.outsize).T
		#e = EMNumPy.numpy2em(img.astype("float32"))
		#scl=float(convnet.outsize)/float(esz)
		#e.clip_inplace(Region(0,0,enx*scl,eny*scl))
		##print e["nx"], e["ny"], img.shape
		#if nframe==1:
			#e.process_inplace("normalize")
		#e.write_image(options.output,-1)
		
		##exit()
	#if nframe>1:
		#ss=options.output
		#fout=ss[:ss.rfind('.')]+"_pp.hdf"
		##### unbin the output and copy the header
		
		#try: os.remove(fout)
		#except: pass
		#if is3d:
			#run("e2proc2d.py {} {} --process math.fft.resample:n={} --twod2threed".format(ss,fout,float(1./convnet.labelshrink)))
			
			#e=EMData(options.tomograms,0,True)
			#e.write_image(fout)
			
		#else:
			#run("e2proc2d.py {} {} --process math.fft.resample:n={}".format(ss,fout,float(1./convnet.labelshrink)))
			#for ni in range(nframe):
				#e=EMData(options.tomograms,ni,True)
				##e.set_size(shape[0],shape[1],1)
				#a=EMData(fout,ni,True)
				#a.set_attr_dict(e.get_attr_dict())
				#a.write_image(fout,ni)
			
		#print "Output written to {}.".format(fout)
	#else:
		#print "Output written to {}.".format(options.output)
	
	#print "Total time:", time.time()-tt0
	
def load_particles(ptcls,labelshrink,ncopy=5, rng=None):
	if rng==None:
		rng=random
	num=EMUtil.get_image_count(ptcls)/2
	
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
	rndid=range(ntrain)
	rng.shuffle(rndid)	
	rndid=rndid+range(ntrain, len(data))
	data=[data[i] for i in rndid]
	label=[label[i] for i in rndid]
	
	print "{:d} particles loaded, {:d} in training set, {:d} in validation set".format(len(data), ntrain, len(data)-ntrain)
	data=np.asarray(data,dtype=theano.config.floatX)
	print "Std of particles: ",np.std(data.flatten())
	#data/=np.std(data.flatten())*3  #np.max(np.abs(data))
	data/=3.
	label=np.asarray(label,dtype=theano.config.floatX)
	label/=np.max(np.abs(label))
	
	
	shared_x = theano.shared(data,borrow=True)
	shared_y = theano.shared(label,borrow=True)
	header=EMData(ptcls,0,True)
	shape=[header["nx"],header["ny"],header["nz"]]
	return shared_x,shared_y,shape, ntrain

class StackedConvNet(object):
	
	def __init__(self,rng,nkernel,ksize,poolsz,imageshape):
		self.amplitude =1.
		self.n_kernel=nkernel
		self.ksize=ksize
		self.n_convlayers=len(self.ksize)
		self.x = T.matrix(name='input')
		self.image_shape=imageshape
		input_shape=self.image_shape
		print "Shape of neural networl input: ",input_shape
		self.convlayers=[]
		self.params=[]
		if poolsz:
			self.poolsz=poolsz
		else:
			self.poolsz=[2]*self.n_convlayers
			
		#self.poolsz[-1]=1
		poolsz=self.poolsz
		convin=self.x
		self.labelshrink=1
		for i in range(self.n_convlayers):
			pz=poolsz[i]
			if pz<0:
				pz=1.0/abs(pz)
			convlayer = LeNetConvPoolLayer(
				rng,
				image_shape=input_shape,
				filter_shape=(self.n_kernel[i], input_shape[1], self.ksize[i], self.ksize[i]),
				poolsize=poolsz[i],
				xin=convin
			)
			self.convlayers.append(convlayer)
			#self.weights.append(convlayer.W)
			
			self.labelshrink=int(self.labelshrink*pz)#poolsz[i]
			input_shape=(input_shape[0],self.n_kernel[i],input_shape[2]/pz,input_shape[3]/pz)
			convin=convlayer.hidden
			
			self.params.extend(convlayer.params)
		
		
		
		self.clslayer=self.convlayers[-1]
		self.outsize=int(input_shape[2])
		#print 
		#self.labelshrink=2**(self.n_convlayers-1)
		
	def get_pretrain_func(self,data,batch_size):
		learning_rate = T.scalar('lr')  # learning rate to use
		weight_decay = T.scalar('wd')  # learning rate to use
		index = T.lscalar()	# index to a [mini]batch
		
		cost, updates = self.convlayer.get_cost_updates(learning_rate,weight_decay)
		train_model = theano.function(
					inputs=[
						index,
						theano.In(learning_rate, value=0.1),
						theano.In(weight_decay, value=1e-5)
					],
					outputs=cost,
					updates=updates,
					givens={
						self.convlayer.x: data[index * batch_size: (index+1) * batch_size]
					}
				)
		return train_model
					
	def get_classify_func(self,data,lab,batch_size):
		learning_rate = T.scalar('lr')  # learning rate to use
		weight_decay = T.scalar('wd')  # learning rate to use
		index = T.lscalar()	# index to a [mini]batch
		
		label=T.matrix(name='label')
		
		cost = self.clslayer.get_cost_hidden(label)
		for c in self.convlayers:
			cost+=weight_decay*T.sum(abs(c.W))
		gparams = T.grad(cost, self.params)
		updates = [
			(param, param - learning_rate * gparam)
			for param, gparam in zip(self.params, gparams)
		]
		#cost, updates = self.clslayer.get_cost_updates_hidden(cls,learning_rate,weight_decay)
		train_classify = theano.function(
					inputs=[
						index,
						theano.In(learning_rate, value=0.1),
						theano.In(weight_decay, value=1e-5)
					],
					outputs=cost,
					updates=updates,
					givens={
						self.x: data[index * batch_size: (index+1) * batch_size],
						label: lab[index * batch_size: (index+1) * batch_size]
					}
					)
		return train_classify
	
	def update_shape(self, imgshape):
		input_shape=imgshape
		if isinstance(self.poolsz,list):
			poolsz=self.poolsz
		else:
			poolsz=[2]*(self.n_convlayers-1)+[1]
		for i in range(self.n_convlayers):
			pz=poolsz[i]
			if pz<0:
				pz=1.0/abs(pz)
			
			self.convlayers[i].image_shape.set_value(input_shape, borrow=True)		
			input_shape=(input_shape[0],self.n_kernel[i],input_shape[2]/pz,input_shape[3]/pz)
			
		self.outsize=int(input_shape[2])
		
		
class LeNetConvPoolLayer(object):

	def __init__(self, rng, filter_shape, image_shape, poolsize=2, xin=None):
		
		assert image_shape[1] == filter_shape[1]
		self.image_shape=theano.shared(
			value=np.asarray(image_shape,dtype='int16'),borrow=True)
		self.poolsize=(poolsize,poolsize)
		#self.input = input
		if xin:
			self.x=xin
		else:
			self.x = T.matrix(name='input')
		self.x1=self.x.reshape(self.image_shape,ndim=4)
		self.filter_shape=filter_shape
		
		# there are "num input feature maps * filter height * filter width"
		# inputs to each hidden unit
		fan_in = np.prod(filter_shape[1:])
		# each unit in the lower layer receives a gradient from:
		# "num output feature maps * filter height * filter width" /
		#   pooling size
		fan_out = (filter_shape[0] * np.prod(filter_shape[2:]) /
			np.prod(self.poolsize))
		# initialize weights with random weights
		W_bound = np.sqrt(6. / (fan_in + fan_out))
		self.W = theano.shared(
			np.asarray(
				rng.uniform(low=-W_bound, high=W_bound, size=filter_shape),
				dtype=theano.config.floatX
			),
			borrow=True
		)
		self.W_prime=self.W[:,:,::-1,::-1]
		self.W_prime=self.W_prime.dimshuffle(1,0,2,3)
		#self.W_prime=self.W_prime[:,::-1]
		#print self.W.get_value()
		#print self.W_prime.eval()
		# the bias is a 1D tensor -- one bias per output feature map
		b_values = np.zeros((filter_shape[0],), dtype=theano.config.floatX)
		bp_values = np.zeros((filter_shape[1],), dtype=theano.config.floatX)
		self.b = theano.shared(value=b_values, borrow=True)
		self.b_prime = theano.shared(value=bp_values, borrow=True)
		
		if poolsize<-1:
			self.x1=self.x1.repeat(int(-poolsize), axis=2).repeat(int(-poolsize), axis=3)

		# convolve input feature maps with filters
		conv_out = conv.conv2d(
			input=self.x1,
			filters=self.W,
			filter_shape=filter_shape,
			#image_shape=self.image_shape.eval(),
			border_mode='full'
		)
		bp=(filter_shape[2]-1)/2
		
		conv_out=conv_out[:,:,bp:-bp,bp:-bp]
		
		# downsample each feature map individually, using maxpooling
		if poolsize>1:
			self.pooled_out = pool.pool_2d(
				input=conv_out,
				ws=self.poolsize,
				ignore_border=True
			)
		else:
			self.pooled_out=conv_out
		
		self.hidden = T.maximum(0,(self.pooled_out + self.b.dimshuffle('x', 0, 'x', 'x')))

		# store parameters of this layer
		self.params = [self.W,self.b]
		#self.params = [self.W]
	
	
	def get_reconstructed_input(self):
		""" Computes the reconstructed input given the values of the hidden layer """

		repeated_conv = conv.conv2d(input = self.hidden,
			      filters = self.W_prime,
			      border_mode='full')
		#repeated_conv=repeated_conv[:,:,1:-1,1:-1]
		bp=(self.filter_shape[2]-1)/2
		repeated_conv=repeated_conv[:,:,bp:-bp,bp:-bp]
		
		multiple_conv_out = [repeated_conv.flatten()] * np.prod(self.poolsize)
		stacked_conv_neibs = T.stack(*multiple_conv_out).T

		stretch_unpooling_out =  T.nnet.neighbours.neibs2images(stacked_conv_neibs, 
							 self.poolsize, self.x1.shape)
		#return self.hidden
		z=T.tanh(stretch_unpooling_out + self.b_prime.dimshuffle('x', 0, 'x', 'x'))
		#return T.sum(T.sum((self.x1-z)**2,axis=3),axis=2)
		#rectified_linear_activation = lambda x: T.maximum(0.0, x)
		return z

	def get_image(self, relu=True):
		#return T.tanh(self.hidden)*2-1
		#return T.maximum(0,(self.pooled_out + self.b.dimshuffle('x', 0, 'x', 'x')))
		if relu:
			return T.minimum(1,(self.pooled_out + self.b.dimshuffle('x', 0, 'x', 'x')))
		else:
			return self.pooled_out + self.b.dimshuffle('x', 0, 'x', 'x')
	
	def get_cost_hidden(self,  label):
		
		
		#z = self.hidden
		#shp=z.shape
		#z = z.reshape((shp[0],shp[2]*shp[3]))
		#z = T.mean(z,axis=1)
		#x=label-.5
		###cost=T.mean((z-x)**2)
		##cost=1-T.mean(z*x)
		#cost=1-T.mean(z*x)
		#cost+= T.sum(abs(self.W))*weight_decay
		
		
		
		xin=label
		#xin=xin.reshape(self.hidden.shape,ndim=4)
		z = self.get_image()
		
		shp=z.shape
		
		xin=xin.reshape((shp[0],-1))
		z=z.reshape((shp[0],-1))
		
		#L = T.sum(xin*z,axis=1)/T.sqrt(T.sum(xin**2,axis=1))/T.sqrt(T.sum(z**2,axis=1))
		#cost = 1-T.mean(L)
		
		#L = T.mean(xin*z,axis=1)
		#cost = 1-T.mean(L)
		
		#cost=T.log(T.mean(((xin-z)**2)*T.exp(0.7*(xin-z)) ))
		cost=T.log(T.mean(((xin-z)**2) ))
		#cost+= T.sum(abs(self.W))*weight_decay
		

		return cost
	
	def get_cost_updates(self, learning_rate,weight_decay,xin=None):
		
			""" This function computes the cost and the updates for one trainng
			step of the dA """
			if xin==None:
				xin=self.x1
				
			z = self.get_reconstructed_input()
			
			shp=self.image_shape
			xin=xin.reshape((shp[0],shp[2]*shp[3]))
			z=z.reshape((shp[0],shp[2]*shp[3]))
			#L = T.sum(T.sum((xin-z)**2,axis=3),axis=2)
			L = T.sum(xin*z,axis=1)/T.sqrt(T.sum(xin**2,axis=1))/T.sqrt(T.sum(z**2,axis=1))
			cost = 1-T.mean(L)
			#L = T.sum(xin*z,axis=1)/T.sqrt(T.sum(xin**2,axis=1))/T.sqrt(T.sum(z**2,axis=1)
			# note : L is now a vector, where each element is the
			#	cross-entropy cost of the reconstruction of the
			#	corresponding example of the minibatch. We need to
			#	compute the average of all these to get the cost of
			#	the minibatch
			#cost = T.mean(L)
			cost+= T.sum(abs(self.W))*weight_decay
			# compute the gradients of the cost of the `dA` with respect
			# to its parameters
			gparams = T.grad(cost, self.params)
			# generate the list of updates
			updates = [
				(param, param - learning_rate * gparam)
				for param, gparam in zip(self.params, gparams)
			]

			return (cost, updates)

if __name__ == '__main__':
    main()
