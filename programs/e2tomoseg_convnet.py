#!/usr/bin/env python
# Muyuan July 2015
import sys
import random
import numpy as np
from EMAN2 import *
import cPickle

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
	parser.add_argument("--netout", type=str,help="Output neural net file name", default="conv.save",guitype='strbox', row=3, col=0, rowspan=1, colspan=3, mode="train")
	
	parser.add_argument("--learnrate", type=float,help="Learning rate ", default=.01, guitype='floatbox', row=4, col=0, rowspan=1, colspan=1, mode="train")
	parser.add_argument("--niter", type=int,help="Training iterations", default=10, guitype='intbox', row=4, col=1, rowspan=1, colspan=1, mode="train")
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
	parser.add_argument("--dream", action="store_true", default=False ,help="Iterativly applying the neural network on noise", guitype='boolbox', row=5, col=0, rowspan=1, colspan=1, mode='test')
	parser.add_argument("--output", type=str,help="Segmentation out file name", default="tomosegresult.mrcs", guitype='strbox', row=3, col=0, rowspan=1, colspan=1, mode="test")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)

	(options, args) = parser.parse_args()
	E2n=E2init(sys.argv,options.ppid)
	
	options.nkernel=[int(i) for i in options.nkernel.split(',')]
	options.ksize=[int(i) for i in options.ksize.split(',')]
	if options.poolsz:
		options.poolsz=[int(i) for i in options.poolsz.split(',')]
	
	if options.applying:
		os.environ["THEANO_FLAGS"]="optimizer=None"
		print "Testing on big images, Theano optimizer disabled"
		import_theano()
		convnet=load_model(options.from_trained)
		if options.dream:
			dream(convnet,options)
		else:
			apply_neuralnet(convnet,options)
		print "Done"
		E2end(E2n)
		exit()
	else:
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
	if options.trainset.endswith(".pkl"):
		f = file(options.trainset, 'rb')
		particles=cPickle.load(f)
		f.close()
	else:
		particles=load_particles(options.trainset,labelshrink,options.ncopy)
		f = file("data_training.pkl", 'wb')
		cPickle.dump(particles, f, protocol=cPickle.HIGHEST_PROTOCOL)
		f.close()

	train_set_x= particles[0]
	labels=particles[1]
	shape=particles[2]
	print "Number of particles: {}".format(train_set_x.shape.eval()[0])
	
	# allocate symbolic variables for the data
	index = T.lscalar()	# index to a [mini]batch
	x = T.matrix('x')  # the data is presented as rasterized images
	image_shape=(batch_size, shape[2], shape[0],shape[1])
	print image_shape,shape
	
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
		for epoch in xrange(options.niter):
		# go through the training set
			c = []
			if epoch==0:
				print classify(0,lr=learning_rate,wd=options.weightdecay)
			for batch_index in xrange(n_train_batches):
				err=classify(batch_index,
					lr=learning_rate,
					wd=options.weightdecay)
				c.append(err)
				if epoch==0 and batch_index<5:
					print err
			learning_rate*=.9
			print 'Training epoch %d, cost ' % ( epoch),
			print np.mean(c),", learning rate",learning_rate

		print "Saving the trained net to {}...".format(options.netout)
		f = file(options.netout, 'wb')
		cPickle.dump(convnet, f, protocol=cPickle.HIGHEST_PROTOCOL)
		f.close()
		
	#######################################
	#print convnet.clslayer.W.get_value()
	#print convnet.clslayer.b.get_value()
			
	if options.trainout:
		print "Generating results ..."
		test_imgs = theano.function(
			inputs=[index],
			outputs=convnet.convlayers[0].get_reconstructed_input(),
			givens={
				convnet.x: train_set_x[index * batch_size: (index+1) * batch_size]
			}
		)
		test_cls = theano.function(
			inputs=[index],
			outputs=convnet.clslayer.get_image(),
			givens={
				convnet.x: train_set_x[index * batch_size: (index+1) * batch_size]
			}
		)
		fname="result_conv_{}.hdf".format(options.netout)
		try:os.remove(fname)
		except: pass
		for idi in range(2):
			rt=test_imgs(idi)
			mid=test_cls(idi)
			mid_cent=mid
			mid_mean=np.mean(mid_cent)
			mid_std=np.std(mid_cent)
			#print "mean:",mid_mean,"std:",mid_std
			#print np.shape(test_imgs(0))

			ipt= train_set_x[idi * batch_size: (idi + 1) * batch_size]
			ipt= ipt.eval()
			
			lb= labels[idi * batch_size: (idi + 1) * batch_size].eval()
			
			for t in range(len(rt)):
				#img=ipt[t].reshape(lth,lth)
				#print ipt[t].shape,shape
				if shape[2]==1:
					img=ipt[t].reshape(shape[0],shape[1])
				else:
					img=ipt[t].reshape(shape[2],shape[0],shape[1])
					img=img[shape[2]/2]
				e = from_numpy(img.astype("float32"))
				#e.mult(-1)
				e.process_inplace("normalize")
				e.write_image(fname,-1)
				
				img=lb[t].reshape(convnet.outsize,convnet.outsize)
				e = from_numpy(img.astype("float32"))
				#e.process_inplace("normalize")
				e.mult(2)
				e=e.get_clip(Region((convnet.outsize-shape[0])/2,(convnet.outsize-shape[0])/2,shape[0],shape[0]))
				e.scale(float(shape[0])/float(convnet.outsize))
				e.write_image(fname,-1)
				
				img=mid[t].reshape(convnet.outsize,convnet.outsize)
				df=(np.mean(img)-mid_mean)/mid_std
				e = from_numpy(img.astype("float32"))
				#e.process_inplace("normalize")
				e.div(float(mid_std))
				
				e=e.get_clip(Region((convnet.outsize-shape[0])/2,(convnet.outsize-shape[0])/2,shape[0],shape[0]))
				#print float(shape[0])/float(convnet.outsize)
				e.scale(float(shape[0])/float(convnet.outsize))
				e.write_image(fname,-1)
		print "Writing output on training set in {}".format(fname)
	print "Done"
	E2end(E2n)

def run(cmd):
	print cmd
	launch_childprocess(cmd)
	
def load_model(fname):
	print "loading model from {}...".format(fname)
	f = file(fname, 'rb')
	convnet = cPickle.load(f)
	f.close()
	return convnet

def dream(convnet,options):
	convz=convnet.image_shape[1]
	
	print "Dreaming....."
	try: 
		os.remove(options.output)
		print "Overwriting the output.."
	except: pass
	
	sz=500
	shape=[sz,sz,convz]
	img=np.random.randn(shape[2],shape[0],shape[1])
	convnet.update_shape((1, shape[2], shape[0],shape[1]))
	e=from_numpy(img)
	e.write_image(options.output,-1)
	
	for it in range(50):
		
		
		print "Iteration {}...".format(it)
		
		
		#### prepare inputs
		img[np.abs(img)>1.0]=1.0
		#print img.shape
		#print convz,sz
		ar=img.reshape((convz,sz*sz))
		
		data=theano.shared(np.asarray(ar,dtype=theano.config.floatX),borrow=True)

		#### write input when testing...
		#### 
		#print np.shape(data.get_value())
		#print "Applying the convolution net..."
		test_imgs = theano.function(
			inputs=[],
			outputs=convnet.clslayer.get_image(),
			givens={
				convnet.x: data
			}
		)
			
		img=test_imgs()
		#print np.shape(img)
		img=img.reshape(convnet.outsize,convnet.outsize).T
		
		e = from_numpy(img)
		e.process_inplace("math.fft.resample",{"n":float(1./convnet.labelshrink)})
		e.process_inplace("normalize")
		e.write_image(options.output,-1)
		if convz>1:
			img=np.random.randn(shape[2],shape[0],shape[1])
			img[convz/2]=e.numpy().copy()
		else:
			img=e.numpy().copy()
	
	print "Output written to {}.".format(options.output)
	
def apply_neuralnet(convnet,options):
	
	
	convz=convnet.image_shape[1]
	
	try: 
		os.remove(options.output)
		print "Overwriting the output.."
	except: pass
	
	
	nframe=EMUtil.get_image_count(options.tomograms)
	is3d=False
	### deal with 3D volume or image stack
	e=EMData(options.tomograms, 0, True)
	if nframe==1:
		nframe=e["nz"]
		if nframe>1:
			#### input data is 3D volume
			esz=e["nx"]
			is3d=True
			
			
	shape=[e["nx"],e["ny"],convz]		
	convnet.update_shape((1, shape[2], shape[0],shape[1]))
	for nf in range(nframe):
		if is3d:
			e=EMData(options.tomograms,0,False,Region(0,0,nf,esz,esz,convz))
		else:
			e=EMData(options.tomograms,nf)
		
		print "Applying the convolution net on image No {}...".format(nf)

		#### prepare inputs
		#e.process_inplace("normalize")
		enp=e.numpy()
		#print enp.shape
		#exit()
		enp[enp>3.0]=3.0
		enp/=3
		ar=enp.reshape((convz,shape[0]*shape[1]))
		
		
		data=theano.shared(np.asarray(ar,dtype=theano.config.floatX),borrow=True)
		
	
		#### input is one 2D slice: just test the performance.
		#### write input when testing...
		if nframe==1:
			img=data[0].eval().reshape(shape[0],shape[1]).T
			e = EMNumPy.numpy2em(img.astype("float32"))
			#e.process_inplace("normalize")
			newshp=convnet.outsize
			e.scale(float(newshp)/float(shape[0]))
			e=e.get_clip(Region((shape[0]-newshp)/2,(shape[0]-newshp)/2,newshp,newshp))
			e.process_inplace("normalize")
			e.write_image(options.output,-1)
		
		#### Applying the convolution net...
		test_imgs = theano.function(
			inputs=[],
			outputs=convnet.clslayer.get_image(),
			givens={
				convnet.x: data
			}
		)
			
		img=test_imgs()
		#print np.shape(img)
		img=img.reshape(convnet.outsize,convnet.outsize).T
		e = EMNumPy.numpy2em(img.astype("float32"))
		
		if nframe==1:
			e.process_inplace("normalize")
		e.write_image(options.output,-1)
		
		#exit()
	if nframe>1:
		ss=options.output
		fout=ss[:ss.rfind('.')]+"_pp.hdf"
		#### unbin the output and copy the header
		
		try: os.remove(fout)
		except: pass
		if is3d:
			run("e2proc2d.py {} {} --process math.fft.resample:n={} --twod2threed".format(ss,fout,float(1./convnet.labelshrink)))
			
			e=EMData(options.tomograms,0,True)
			e.write_image(fout)
			
		else:
			run("e2proc2d.py {} {} --process math.fft.resample:n={}".format(ss,fout,float(1./convnet.labelshrink)))
			for ni in range(nframe):
				e=EMData(options.tomograms,ni,True)
				#e.set_size(shape[0],shape[1],1)
				a=EMData(fout,ni,True)
				a.set_attr_dict(e.get_attr_dict())
				a.write_image(fout,ni)
			
		print "Output written to {}.".format(fout)
	else:
		print "Output written to {}.".format(options.output)

def load_particles(ptcls,labelshrink,ncopy=5):
	num=EMUtil.get_image_count(ptcls)/2
	
	data=[]
	label=[]
	for nc in range(ncopy):
		for i in range(num):
			ptl=EMData(ptcls,i*2)
			#ptl.process_inplace("threshold.belowtozero")
			if ncopy>1:
				tr=Transform()
				tr.set_rotation({"type":"2d","alpha":random.random()*360.0})
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
	
	## randomize
	rndid=range(len(data))
	random.shuffle(rndid)	
	data=[data[i] for i in rndid]
	label=[label[i] for i in rndid]
	
	data=np.asarray(data,dtype=theano.config.floatX)
	data/=np.std(data.flatten())*3  #np.max(np.abs(data))
	label=np.asarray(label,dtype=theano.config.floatX)
	label/=np.max(np.abs(label))
	
	
	shared_x = theano.shared(data,borrow=True)
	shared_y = theano.shared(label,borrow=True)
	header=EMData(ptcls,0,True)
	shape=[header["nx"],header["ny"],header["nz"]]
	return shared_x,shared_y,shape

class StackedConvNet(object):
	
	def __init__(self,rng,nkernel,ksize,poolsz,imageshape):
		
		self.n_kernel=nkernel
		self.ksize=ksize
		self.n_convlayers=len(self.ksize)
		self.x = T.matrix(name='input')
		self.image_shape=imageshape
		input_shape=self.image_shape
		print input_shape
		self.convlayers=[]
		self.params=[]
		if poolsz:
			self.poolsz=poolsz
		else:
			self.poolsz=[2]*self.n_convlayers
			
		self.poolsz[-1]=1
		poolsz=self.poolsz
		convin=self.x
		self.labelshrink=1
		for i in range(self.n_convlayers):
			convlayer = LeNetConvPoolLayer(
				rng,
				image_shape=input_shape,
				filter_shape=(self.n_kernel[i], input_shape[1], self.ksize[i], self.ksize[i]),
				poolsize=(poolsz[i], poolsz[i]),
				xin=convin
			)
			self.convlayers.append(convlayer)
			#self.weights.append(convlayer.W)
			self.labelshrink*=poolsz[i]
			input_shape=(input_shape[0],self.n_kernel[i],input_shape[2]/poolsz[i],input_shape[3]/poolsz[i])
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
		
		#negnum=int(batch_size*.3) # number of negative samples
		#theano_rng=RandomStreams(1)
		#tdata=data[index * batch_size: ((index+1) * batch_size)-negnum]
		#noise=theano_rng.normal(size=tdata[:negnum].shape,dtype=theano.config.floatX)
		#tdata=T.concatenate([tdata,noise])
		#cls=T.concatenate([cls[:-negnum], T.zeros_like(cls[:negnum])-1])
		#print cls.shape.eval()
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
			self.convlayers[i].image_shape.set_value(input_shape, borrow=True)		
			input_shape=(input_shape[0],self.n_kernel[i],input_shape[2]/poolsz[i],input_shape[3]/poolsz[i])
			
		self.outsize=int(input_shape[2])
		
		
class LeNetConvPoolLayer(object):

	def __init__(self, rng, filter_shape, image_shape, poolsize=(2, 2), xin=None):
		
		assert image_shape[1] == filter_shape[1]
		self.image_shape=theano.shared(
			value=np.asarray(image_shape,dtype='int16'),borrow=True)
		self.poolsize=poolsize
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
			np.prod(poolsize))
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
		self.pooled_out = pool.pool_2d(
			input=conv_out,
			ds=poolsize,
			ignore_border=True
		)
		#shp=conv_out.shape
		#y = T.nnet.neighbours.images2neibs(conv_out, poolsize,mode='ignore_borders')
		#pooled_out=y.mean(axis=-1).reshape((shp[0],shp[1],shp[2]/poolsize[0],shp[3]/poolsize[1]))

		# add the bias term. Since the bias is a vector (1D array), we first
		# reshape it to a tensor of shape (1, n_filters, 1, 1). Each bias will
		# thus be broadcasted across mini-batches and feature map
		# width & height
		#self.hidden = T.tanh(self.pooled_out + self.b.dimshuffle('x', 0, 'x', 'x'))
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

	def get_image(self):
		#return T.tanh(self.hidden)*2-1
		#return T.maximum(0,(self.pooled_out + self.b.dimshuffle('x', 0, 'x', 'x')))
		return T.minimum(1,(self.pooled_out + self.b.dimshuffle('x', 0, 'x', 'x')))
	
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
		
		xin=xin.reshape((shp[0],shp[2]*shp[3]))
		z=z.reshape((shp[0],shp[2]*shp[3]))
		
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
