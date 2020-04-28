#!/usr/bin/env python
# Muyuan Chen 2017-10
from past.utils import old_div
from future import standard_library
standard_library.install_aliases()
from builtins import range
from EMAN2 import *
import numpy as np
import threading
import queue
import theano
import theano.tensor as T
from theano.tensor.nnet import conv
from theano.tensor.signal import pool
from e2tomoseg_convnet import import_theano, StackedConvNet, save_model, load_model
import_theano()

def main():
	
	usage=""" 
	train network:
		e2convnet_classifyptcls.py [good_refs.hdf] [bad_refs.hdf]
	
	sort particle set (output to xx/xx_sort.lst)
		e2covnet_classifyptcls.py [all_ptcls.lst]
	
	"""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	#parser.add_argument("--autobox", action="store_true", default=False ,help="box particles")
	#parser.add_argument("--boxsz", type=int,help="box size of input particles. Necessery for applying network", default=256)
	#parser.add_argument("--threads", type=int,help="number of threads", default=10)
	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)
	
	if len(args)==2:
		do_training(args)
	else:
		classify_on(args[0])
	
	E2end(logid)
	
def run(cmd):
	print(cmd)
	launch_childprocess(cmd)
	
def classify_on(fname):
	n=EMUtil.get_image_count(fname)
	tstsz=1024
	nnet_savename="nnet_classify.hdf"
	e0=EMData(fname,0)
	bxsz=e0["nx"]
	sz=64
	shrinkfac=old_div(float(bxsz),float(sz))
	convnet=load_model(nnet_savename)
	
	nbatch=n//tstsz+1
	#nbatch=1
	allscore=[]
	for ib in range(nbatch):
		sn=min(tstsz, n-ib*tstsz)
		data=[]
		
		for i in range(sn):
			e=EMData(fname, i+ib*tstsz)
			ptl=e.process("math.fft.resample",{"n":shrinkfac})
			ptl.clip_inplace(Region(0,0, sz, sz))
			#ptl.process_inplace("filter.highpass.gauss",{"cutoff_freq":0.005})
			#ptl.process_inplace("filter.lowpass.gauss",{"cutoff_freq":0.05})
			ptl.process_inplace("filter.highpass.gauss",{"cutoff_pixels":2})
			ptl.process_inplace("filter.lowpass.gauss",{"cutoff_abs":0.25})
			ar=ptl.numpy().copy()
			data.append(ar.flatten())
		data=np.array(data, dtype=theano.config.floatX)
		print(data.shape, ib)
	
		div=np.mean(np.std(data,axis=1))
		data/=div#np.std(data)#*2.
		mx=4.
		data[data>mx]=mx
		data[data<-mx]=-mx
		testset = theano.shared(data,borrow=True)
		
		
		
		convnet.update_shape((sn, 1, sz, sz))
		convnet.sumout=T.mean(convnet.clslayer.get_image().reshape((sn, -1)), axis=1)
		test_cls = theano.function(
			inputs=[],
			outputs=convnet.sumout,
			givens={
				convnet.x: testset
			}
		)
		tstout=test_cls()
		#print tstout
		allscore.append(tstout)
		
	allscore=np.hstack(allscore)
	srt=np.argsort(-allscore)
	np.savetxt("score_srt.txt",np.vstack([np.arange(n), allscore[srt]]).T)
	lstin=LSXFile(fname, True)
	
	lname=fname[:-4]+"_sort.lst"
	try: os.remove(lname)
	except: pass
	lst=LSXFile(lname,False)
	for i in srt:
		ln=lstin.read(i)
		#print i, ln
		lst.write(-1, ln[0], ln[1])
	lst=None
	print(lname)
	lstin=None

def do_training(args=None):
	
	goodrefs, badrefs=[EMData.read_images(a) for a in args]
	
		
	
	nnet_savename="nnet_classify.hdf"
	bxsz=goodrefs[0]["nx"]
	sz=64
	shrinkfac=old_div(float(bxsz),float(sz))
	
	print("Setting up model ...")
	rng = np.random.RandomState(123)
	nkernel=[20,20,1]
	ksize=[15,15,15]
	poolsz=[2,1,1]
	batch_size=10
	
	image_shape=(batch_size, 1, sz, sz)
	convnet = StackedConvNet(
		rng,
		nkernel=nkernel,
		ksize=ksize,
		poolsz=poolsz,
		imageshape=image_shape
	)
	
	convnet.sumout=T.mean(convnet.clslayer.get_image().reshape((batch_size, -1)), axis=1)
	convnet.sumout=T.maximum(0,convnet.sumout)
	#convnet.sumout=T.minimum(1,convnet.sumout)
	
	
	print("Pre-processing particles...")
	#### here we shrink the particles so they are 64x64
	#### and duplicate so there are more than 500 good and 500 bad particles
	
	nref_target=500
	
	data=[] ### particles in flattened numpy array
	lbs=[]  ### labels in flattened numpy array
	
	for label, refs in enumerate([badrefs,goodrefs]):
		nref=len(refs)
		if nref<5:
			print("Not enough references. Please box at least 5 good and 5 bad reference...")
			return []
		ncopy=old_div(nref_target,nref) + 1
		
		for pp in refs:
			ptl=pp.process("math.fft.resample",{"n":shrinkfac})
			ptl.clip_inplace(Region(0,0, sz, sz))
			#ptl.process_inplace("filter.highpass.gauss",{"cutoff_freq":0.005})
			#ptl.process_inplace("filter.lowpass.gauss",{"cutoff_freq":0.05})
			ptl.process_inplace("filter.highpass.gauss",{"cutoff_pixels":2})
			ptl.process_inplace("filter.lowpass.gauss",{"cutoff_abs":0.25})
			for c in range(ncopy):
				tr=Transform()
				tr.set_rotation({"type":"2d","alpha":np.random.random()*360.0})
				img=ptl.process("xform",{"transform":tr})
				ar=img.numpy().copy()
				data.append(ar.flatten())
				lbs.append(label)
	#print shrinkfac
	rndid=list(range(len(data)))
	np.random.shuffle(rndid)
	data=[data[i] for i in rndid]
	lbs=[lbs[i] for i in rndid]
	data=np.asarray(data,dtype=theano.config.floatX)
	#print np.std(data), np.mean(np.std(data,axis=1))
	div=np.mean(np.std(data,axis=1))
	data/=div#np.std(data)#*2.
	mx=4.
	data[data>mx]=mx
	data[data<-mx]=-mx
	lbs=np.asarray(lbs,dtype=theano.config.floatX)
	train_set_x = theano.shared(data,borrow=True)
	
	#### make target output
	#img=EMData(sz/2,sz/2)
	#img.process_inplace("testimage.gaussian",{'sigma':5.})
	#img.div(img["maximum"])
	#gaus=img.numpy().copy().flatten()
	
	#gaus=gaus.astype(theano.config.floatX)
	#lbarrs=np.array([np.zeros_like(gaus, dtype=theano.config.floatX), gaus])
	#label_np=lbarrs[lbs]
	labels=theano.shared(lbs.astype(theano.config.floatX), borrow=True)
	#print lbs.astype(theano.config.floatX)
	#print gaus.shape, data.shape, label_np.shape
	
	print("Now Training...")
	classify=get_classify_func(convnet, train_set_x,labels,batch_size)
	learning_rate=0.005
	weightdecay=1e-5
	n_train_batches = old_div(len(data), batch_size)
	for epoch in range(20):
	# go through the training set
		c = []
		for batch_index in range(n_train_batches):
			err=classify(batch_index,
				lr=learning_rate,
				wd=weightdecay)
			c.append(err)

		learning_rate*=.96
		print('Training epoch %d, cost ' % ( epoch), end=' ')
		print(np.mean(c),", learning rate",learning_rate)

	
	save_model(convnet, nnet_savename)
	
	tstsz=100
	convnet.update_shape((tstsz, 1, sz, sz))
	test_cls = theano.function(
		inputs=[],
		outputs=convnet.clslayer.get_image(),
		givens={
			convnet.x: train_set_x[:tstsz]
		}
	)
	tstout=test_cls()
	trainoutfile="trainout_nnet_classify.hdf"
	if os.path.isfile(trainoutfile):
		os.remove(trainoutfile)
	
	for i,mm in enumerate(tstout):
		t=train_set_x[i].eval().reshape(sz,sz)
		img=from_numpy(t)
		#img.process_inplace("normalize")
		img.write_image(trainoutfile, -1)
		
		for m in mm:
			img=from_numpy(m)
			nx=img["nx"]
			img=img.get_clip(Region(old_div(-nx,2),old_div(-nx,2),nx*2,nx*2))
			img.scale(2.)
			#img.process_inplace("math.fft.resample",{"n":.5})
			#img.mult(5)
			img.process_inplace("threshold.clampminmax.nsigma", {"nsigma":4})
			img.write_image(trainoutfile, -1)
		
def get_classify_func(convnet, data,lab,batch_size):
		learning_rate = T.scalar('lr')  # learning rate to use
		weight_decay = T.scalar('wd')  # learning rate to use
		index = T.lscalar()	# index to a [mini]batch
		
		label=T.vector(name='label')
		cost=T.mean((convnet.sumout-label)**2)
		#cost = convnet.clslayer.get_cost_hidden(label)
		
		
		for c in convnet.convlayers:
			cost+=weight_decay*T.sum(abs(c.W))
		gparams = T.grad(cost, convnet.params)
		updates = [
			(param, param - learning_rate * gparam)
			for param, gparam in zip(convnet.params, gparams)
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
						convnet.x: data[index * batch_size: (index+1) * batch_size],
						label: lab[index * batch_size: (index+1) * batch_size]
					}
					)
		return train_classify

if __name__ == '__main__':
	main()
	
