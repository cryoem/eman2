#!/usr/bin/env python
# Muyuan July 2015
import sys
import random
import numpy as np
from EMAN2 import *
import cPickle

def import_theano():
	global theano,T,conv,downsample
	import theano
	import theano.tensor as T
	from theano.tensor.nnet import conv
	from theano.tensor.signal import downsample


def main():

	usage="""convnet_pickptcls.py  [options]
	*** This program requires Theano and cPickle in addition to normal EMAN2 requirments.
	Start with "convnet_pickptcls.py [particles] --ngtvs [negative samples] --trainout" to train the neural net.
	Look at result_conv0.hdf to see the training result. The images are arranged as (preprocessed particle / first layer output / last layer output ). Idealy last layer output should be gaussian balls for particles, dark for negative samples.
	If the training results looks good, run "convnet_pickptcls.py --teston [micrograph]" and view testresult.hdf to see the result. Also the boxes can be viewed in e2boxer.py
	Use "convnet_pickptcls.py --teston [micrograph folder]" to box all micrographs in a folder
	"""
	
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--ngtvs", type=str,help="Negative samples for training", default=None)
	parser.add_argument("--ncopy", type=int,help="Number of copies made for input particles", default=10)
	parser.add_argument("--learnrate", type=float,help="Learning rate for the auto-encoder", default=.05)
	parser.add_argument("--niter", type=int,help="Training iterations", default=10)
	parser.add_argument("--nkernel", type=str,help="Number of kernels for each layer, from input to output. The number of kernels in the last layer must be 1. Default is 10,5,1", default="10,5,1")
	parser.add_argument("--ksize", type=str,help="Width of kernels of each layer, the numbers must be odd. Note the number of layers should be the same as the nkernel option. Default is 11,11,5", default="11,11,5")
	parser.add_argument("--weightdecay", type=float,help="Weight decay. Used for regularization.", default=1e-6)
	parser.add_argument("--batch", type=int,help="Batch size for the stochastic gradient descent. Default is 30.", default=30)
	parser.add_argument("--tarsize", type=float,help="Sigma of the target gaussian output,default is 15", default=13.)
	parser.add_argument("--pretrainnet", type=str,help="Output pre-trained neural net temp file name", default="conv.save")
	parser.add_argument("--trainout", action="store_true", default=False ,help="Output the result of the training set")
	parser.add_argument("--fromlast", action="store_true", default=False ,help="Start from previous pretrained network")
	parser.add_argument("--teston", type=str,help="test on one image. The image must be square", default=None)
	parser.add_argument("--boxptcl_noedge", type=int,help="Ignore edge of the micrograph while boxing particles", default=50)
	parser.add_argument("--boxptcl_boxsep", type=int,help="Minimum seperation for particle picking", default=25)
	parser.add_argument("--boxptcl_minscore", type=int,help="Minimum score for particle picking", default=.1)
	parser.add_argument("--shrink", type=int,help="Shrink particles", default=1)
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)

	(options, args) = parser.parse_args()
	E2n=E2init(sys.argv,options.ppid)
	
	options.nkernel=[int(i) for i in options.nkernel.split(',')]
	options.ksize=[int(i) for i in options.ksize.split(',')]
	
	if options.teston!=None:
		os.environ["THEANO_FLAGS"]="optimizer=None"
		print "Testing on big images, Theano optimizer disabled"
		import_theano()
		convnet=load_model(options.pretrainnet)
		box_particles(convnet,options)
		exit()
	else:
		os.environ["THEANO_FLAGS"]="optimizer=fast_run"
		import_theano()
	
	
	batch_size=options.batch
	#### Train da with particles first.
	
	if len(args)==0:
		print "no particle input...exit."
		exit()
	
	print "loading particles..."
	if args[0].endswith(".pkl"):
		f = file(args[0], 'rb')
		particles=cPickle.load(f)
		f.close()
	else:
		particles=load_particles(args[0],options)
		f = file("data_training.pkl", 'wb')
		cPickle.dump(particles, f, protocol=cPickle.HIGHEST_PROTOCOL)
		f.close()

	train_set_x= particles[0]
	labels=particles[1]
	shape=particles[2]
	

	# allocate symbolic variables for the data
	index = T.lscalar()	# index to a [mini]batch
	x = T.matrix('x')  # the data is presented as rasterized images
	image_shape=(batch_size, 1, shape[0],shape[1])
	
	
	rng = np.random.RandomState(123)
	#print shape
	if options.fromlast:
		convnet=load_model(options.pretrainnet)
	
	else:
		print "setting up model"
		convnet = StackedConvNet(
			rng,
			image_shape=image_shape,
			nkernel=options.nkernel,
			ksize=options.ksize
		)
	
	convnet.update_shape(image_shape)
	
	if (options.niter>0):	
		print "training the convolutional network..."
		
		classify=convnet.get_classify_func(train_set_x,labels,batch_size,options.tarsize)
			
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

		print "Saving the trained net to file..."
		f = file(options.pretrainnet, 'wb')
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
		fname="result_conv0.hdf"
		try:os.remove(fname)
		except: pass
		for idi in range(3):
			rt=test_imgs(idi)
			mid=test_cls(idi)
			mid_cent=mid
			mid_mean=np.mean(mid_cent)
			mid_std=np.std(mid_cent)
			print "mean:",mid_mean,"std:",mid_std
			#print np.shape(test_imgs(0))

			ipt= train_set_x[idi * batch_size: (idi + 1) * batch_size]
			ipt= ipt.eval()
			
			lb= labels[idi * batch_size: (idi + 1) * batch_size].eval()
			
			for t in range(len(rt)):
				#img=ipt[t].reshape(lth,lth)
				
				img=ipt[t].reshape(shape[0],shape[1])
				e = EMNumPy.numpy2em(img.astype("float32"))
				e.process_inplace("normalize")
				e["label"]=float(lb[t])
				e.write_image(fname,-1)
				
				img=rt[t].reshape(shape[0],shape[1])
				e = EMNumPy.numpy2em(img.astype("float32"))
				e.process_inplace("normalize")
				e["label"]=float(lb[t])
				e.write_image(fname,-1)
				
				img=mid[t].reshape(convnet.outsize,convnet.outsize)
				df=(np.mean(img)-mid_mean)/mid_std
				#print lb[t],df,df*lb[t]
				#print img
				e = EMNumPy.numpy2em(img.astype("float32"))
				e.sub(float(mid_mean))
				e.div(float(mid_std))
				e=e.get_clip(Region((convnet.outsize-shape[0])/2,(convnet.outsize-shape[0])/2,shape[0],shape[0]))
				#print float(shape[0])/float(convnet.outsize)
				e.scale(float(shape[0]+8)/float(convnet.outsize))
				e["label"]=float(lb[t])
				e.write_image(fname,-1)

def load_model(fname):
	print "loading model from {}...".format(fname)
	f = file(fname, 'rb')
	convnet = cPickle.load(f)
	f.close()
	return convnet
	
def box_particles(convnet,options):
	
	filelist=[]
	isfolder=False
	if os.path.isdir(options.teston):
		print "Check files..."
		isfolder=True
		lst=sorted(os.listdir(options.teston))
		lst=[i for i in lst if i[0]!="."]
		for mg in lst:
			print "  ",mg
			filelist.append(os.path.join(options.teston,mg))
	else:
		filelist.append(options.teston)
	
	if not isfolder:
		try: os.remove("testresult.hdf")
		except: pass
	
	lastshape=-1
	
	for tarfile in filelist:
		print "Boxing particles on ",tarfile
		try:
			e=EMData(tarfile)
		except:
			print "Cannot read file, continue"
			continue
		
		print "\tPreprocessing image..."

		e.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.2})
		e.process_inplace("filter.highpass.gauss",{"cutoff_abs":.005})
		e.process_inplace("normalize.edgemean")
		e.process_inplace("threshold.belowtominval",{"minval":-3, "newval":-3})
		e.mult(-1)
		e.process_inplace("threshold.belowtominval",{"minval":-3, "newval":-3})
		e.mult(-1)
		
		if (options.shrink>1):
			e.process_inplace("math.meanshrink",{"n":options.shrink})
		
		eg=options.boxptcl_noedge
		e.process_inplace("mask.zeroedge2d",{"x0":eg,"x1":eg,"y0":eg,"y1":eg})
		
		## make the image square
		ori_shape=[e["nx"],e["ny"]]
		if (e["nx"]!=e["ny"]):
			shape=max(e["nx"],e["ny"])
			e=e.get_clip(Region((e["nx"]-shape)/2,(e["ny"]-shape)/2,shape,shape))
		else:
			shape=e["nx"]
		
		e.div(e["maximum"])
		ar=[EMNumPy.em2numpy(e)]
		
		if shape!=lastshape:	## only update when shape changes
			convnet.update_shape((1, 1, shape,shape))
			
			
		data=theano.shared(np.asarray(ar,dtype=theano.config.floatX),borrow=True)
		img=data[0].eval().reshape(shape,shape).T
		e = EMNumPy.numpy2em(img.astype("float32"))
		e.process_inplace("normalize")
		newshp=convnet.outsize
		e.scale(float(newshp)/float(shape))
		ori_shape=[i*(float(newshp)/float(shape)) for i in ori_shape]
		e=e.get_clip(Region((shape-newshp)/2,(shape-newshp)/2,newshp,newshp))
		e=e.get_clip(Region((e["nx"]-ori_shape[0])/2,(e["ny"]-ori_shape[1])/2,ori_shape[0],ori_shape[1]))
		if not isfolder: e.write_image("testresult.hdf",-1)

		
		print "\tApplying the convolution net..."
		test_imgs = theano.function(
			inputs=[],
			outputs=convnet.clslayer.hidden,
			givens={
				convnet.x: data[0]
			}
		)
			
		img=test_imgs()
		#print np.shape(img)
		img=img.reshape(convnet.outsize,convnet.outsize).T
		e = EMNumPy.numpy2em(img.astype("float32"))
		#print "Post-processing..."
		#eg=20
		#e.process_inplace("mask.zeroedge2d",{"x0":eg,"x1":eg,"y0":eg,"y1":eg})
		e.process_inplace("normalize")
		#e.process_inplace("threshold.belowtozero",{"minval":0})
		#e.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.05})
		#e.div(e["maximum"])
		e=e.get_clip(Region((e["nx"]-ori_shape[0])/2,(e["ny"]-ori_shape[1])/2,ori_shape[0],ori_shape[1]))
		if not isfolder: e.write_image("testresult.hdf",-1)
		
		###########
		print "\tBoxing particles.."
		#e=EMData("testresult.hdf",1)
		e.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.2})
		#shape=[e["nx"],e["ny"],e["nz"]]
		#convnet.update_shape((1, 1, shape[0],shape[1]))
		#newshp=convnet.outsize
		pks=e.peak_ccf(options.boxptcl_boxsep)
		#print pks
		
		box=[]
		scale=float(shape)/float(newshp)
		for i in range(0,len(pks),3):
			#print i//3,pks[i]
			if pks[i]<options.boxptcl_minscore:
				break
			box.append([pks[i+1]*scale,pks[i+2]*scale,"manual"])
		
		print "\t{} particles.".format(i//3)
		import json
		jn={}
		jn["boxes"]=box
		#print jn
		#options.teston="testimg.hdf"
		mn=tarfile.replace('/','-')
		#n1=mn.rfind('/')+1
		#mn=mn[n1:]
		#print mn
		f = open('info/{n}_info.json'.format(n=mn[:-4]), 'w')
		json.dump(jn, f, indent=0)
		f.close()
		lastshape=shape

def load_particles(ptcls, options):
	
	nptcls=EMUtil.get_image_count(ptcls)
	if options.ngtvs==None:
		label=[1]
	else:
		label=[1,0]
		
	source=[options.ngtvs,ptcls]
	lbs=[] # ptcls are labeled 1, 
	data=[]
	for l in label:
		num=EMUtil.get_image_count(source[l])
		if l==0:
			nc=options.ncopy*nptcls/num
		else:
			nc=options.ncopy
		for i in range(num):
			ptl=EMData(source[l],i)
			ptl.process_inplace("normalize.edgemean")
			ptl.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.2})
			ptl.process_inplace("filter.highpass.gauss",{"cutoff_abs":.005})
			if (options.shrink>1):
				ptl.process_inplace("math.meanshrink",{"n":options.shrink})
			for c in range(nc):
				tr=Transform()
				tr.set_rotation({"type":"2d","alpha":random.random()*360.0})
				img=ptl.process("xform",{"transform":tr})
				
				img.write_image("tst.hdf",-1)
				ar=EMNumPy.em2numpy(img)
				data.append(ar.flatten())
				lbs.append(l)
	
	rndid=range(len(data))
	random.shuffle(rndid)	
	data=[data[i] for i in rndid]
	lbs=[lbs[i] for i in rndid]
	data=np.asarray(data,dtype=theano.config.floatX)
	data/=np.max(np.abs(data))
	lbs=np.asarray(lbs,dtype=theano.config.floatX)
	
	#header=EMData(ptcls,0,True)
	shape=[ptl["nx"],ptl["ny"],ptl["nz"]]
	shared_data = [theano.shared(data,borrow=True),theano.shared(lbs,borrow=True),shape]

	return shared_data


def gen_classify_target(img_size,sigma=10):
	img=EMData(img_size[0],img_size[1])
	img.process_inplace("testimage.gaussian",{'sigma':sigma})
	img.div(img["maximum"])
	#img.mult(.5)
	#img.add(.5)
	img.write_image("gaussian_example.hdf")
	#img.mult(5)
		
	ar=EMNumPy.em2numpy(img)
	shp=np.shape(ar)
	data=ar.flatten()
	
	shared_x = theano.shared(np.asarray(data,dtype=theano.config.floatX),borrow=True)
	
	return shared_x

class StackedConvNet(object):
	
	def __init__(self,rng,image_shape,nkernel,ksize):
		
		self.n_kernel=nkernel
		self.ksize=ksize
		self.n_convlayers=len(self.ksize)
		self.x = T.matrix(name='input')
		self.image_shape=image_shape
		input_shape=image_shape
		self.convlayers=[]
		self.params=[]
		self.poolsz=2
		poolsz=self.poolsz
		convin=self.x
		for i in range(self.n_convlayers):
			if i==self.n_convlayers-1:
				poolsz=1
			convlayer = LeNetConvPoolLayer(
				rng,
				image_shape=input_shape,
				filter_shape=(self.n_kernel[i], input_shape[1], self.ksize[i], self.ksize[i]),
				poolsize=(poolsz, poolsz),
				xin=convin
			)
			self.convlayers.append(convlayer)
			#self.weights.append(convlayer.W)
			input_shape=(input_shape[0],self.n_kernel[i],input_shape[2]/poolsz,input_shape[3]/poolsz)
			convin=convlayer.hidden
			
			self.params.extend(convlayer.params)
		
		
		
		self.clslayer=self.convlayers[-1]
		self.outsize=int(input_shape[2])
		
	def get_pretrain_func(self,data,batch_size):
		learning_rate = T.scalar('lr')  # learning rate to use
		weight_decay = T.scalar('wd')  # learning rate to use
		index = T.lscalar()	# index to a [mini]batch
		
		cost, updates = self.convlayer.get_cost_updates(learning_rate,weight_decay)
		train_model = theano.function(
					inputs=[
						index,
						theano.Param(learning_rate, default=0.1),
						theano.Param(weight_decay, default=1e-5)
					],
					outputs=cost,
					updates=updates,
					givens={
						self.convlayer.x: data[index * batch_size: (index+1) * batch_size]
					}
				)
		return train_model
					
	def get_classify_func(self,data,lab,batch_size,sigma=15):
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
		label=T.vector(name='label')
		
		target=gen_classify_target((self.outsize,self.outsize),sigma)
		cost = self.clslayer.get_cost_hidden(target,label)
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
						theano.Param(learning_rate, default=0.1),
						theano.Param(weight_decay, default=1e-5)
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
		poolsz=self.poolsz
		for i in range(self.n_convlayers):
			if i==self.n_convlayers-1:
				poolsz=1
			self.convlayers[i].image_shape.set_value(input_shape, borrow=True)		
			input_shape=(input_shape[0],self.n_kernel[i],input_shape[2]/poolsz,input_shape[3]/poolsz)
			
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
		self.pooled_out = downsample.max_pool_2d(
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
		return T.minimum(1,(self.hidden))
	
	def get_cost_hidden(self, tar, label):
		
		
		#z = self.hidden
		#shp=z.shape
		#z = z.reshape((shp[0],shp[2]*shp[3]))
		#z = T.mean(z,axis=1)
		#x=label-.5
		###cost=T.mean((z-x)**2)
		##cost=1-T.mean(z*x)
		#cost=1-T.mean(z*x)
		#cost+= T.sum(abs(self.W))*weight_decay
		
		
		
		
		
		xin=T.outer(label,tar)
		#xin=xin.reshape(self.hidden.shape,ndim=4)
		z = self.get_image()
		
		shp=z.shape
		
		xin=xin.reshape((shp[0],shp[2]*shp[3]))
		z=z.reshape((shp[0],shp[2]*shp[3]))
		
		#L = T.sum(xin*z,axis=1)/T.sqrt(T.sum(xin**2,axis=1))/T.sqrt(T.sum(z**2,axis=1))
		#cost = 1-T.mean(L)
		
		#L = T.mean(xin*z,axis=1)
		#cost = 1-T.mean(L)
		cost=T.mean((xin-z)**2)
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
