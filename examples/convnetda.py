#!/usr/bin/env python
# Muyuan July 2015
import sys
import random
import numpy as np
from EMAN2 import *
import theano
import theano.tensor as T
from theano.tensor.shared_randomstreams import RandomStreams
import cPickle
from theano.tensor.nnet import conv
from theano.tensor.signal import downsample
def main():


	usage="""datest.py <particles> [options]
	Run a stacked denoise auto-encoder on the data set."""
	
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--learnrate", type=float,help="Learning rate for the auto-encoder", default=.0001)
	parser.add_argument("--niter", type=int,help="Training iterations", default=10)
	parser.add_argument("--append", type=int,help="Append a layer", default=-1)
	parser.add_argument("--layer", type=str,help="train only these layers", default=None)
	parser.add_argument("--weightdecay", type=float,help="Weight decay. Used for regularization.", default=0)
	parser.add_argument("--batch", type=int,help="Batch size for the stochastic gradient descent.", default=50)
	parser.add_argument("--shrink", type=int,help="Shrink the image.", default=1)
	parser.add_argument("--pretrainnet", type=str,help="Output pretrained autoencoder file name", default="conv.save")
	parser.add_argument("--hidden", type=str,help="Size of hidden layers: n1,n2,n3...", default="1000,500")
	parser.add_argument("--sdaout", action="store_true", default=False ,help="Output the result of the autoencoder")
	parser.add_argument("--checkmid", action="store_true", default=False ,help="check the middle layer output")
	parser.add_argument("--fromlast", action="store_true", default=False ,help="Start from previous pretrained network")
	parser.add_argument("--writeall", type=str, default=None ,help="write all output")
	parser.add_argument("--weights", type=str,help="output file name for the weights", default=None)
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)

	(options, args) = parser.parse_args()
	E2n=E2init(sys.argv,options.ppid)
	training_epochs=options.niter
	
	batch_size=options.batch
	hidden=[int(i) for i in options.hidden.split(',')]
	#### Train da with particles first.
	

	print "loading particles..."
	#particles = load_mat(args[1])
	##print train_set_x.get_value()
	particles = load_particles(args[0],shrink=options.shrink)
	
	#f = file("ptcls_test.pkl", 'wb')
	#cPickle.dump(particles, f, protocol=cPickle.HIGHEST_PROTOCOL)
	##particles=cPickle.load(f)
	#f.close()
	
	
	train_set_x= particles
	
	ltrain=train_set_x.get_value().shape[1]
	lth=int(sqrt(ltrain))
	
	ptl=EMData(args[0],0)
	shape=[ptl["nx"],ptl["ny"],ptl["nz"]]
	shape=[i/options.shrink for i in shape]
	print shape
	# compute number of minibatches 
	

	# allocate symbolic variables for the data
	index = T.lscalar()	# index to a [mini]batch
	x = T.matrix('x')  # the data is presented as rasterized images
	image_shape=(batch_size, 1, shape[0],shape[1])
	
	print "setting up model"
	rng = np.random.RandomState(123)
	#print shape
	if options.fromlast:
		print "loading {}...".format(options.pretrainnet)
		f = file(options.pretrainnet, 'rb')
		convnet = cPickle.load(f)
		f.close()		
		x=convnet.x
	
	else:
		convnet = LeNetConvPoolLayer(
			rng,
			image_shape=image_shape,
			filter_shape=(5, 1, 15, 15),
			poolsize=(2, 2)
		)
	convnet.image_shape.set_value(image_shape, borrow=True)
	print convnet.image_shape.eval()
	#test_imgs = theano.function(
			#inputs=[index],
			#outputs=convnet.get_reconstructed_input(),
			#givens={
				#convnet.x: train_set_x[index * batch_size: (index+1) * batch_size]
			#}
		#)
	##print test_imgs(0)
	#print np.shape(test_imgs(0))
	#f = file(options.pretrainnet, 'wb')
	#cPickle.dump(convnet, f, protocol=cPickle.HIGHEST_PROTOCOL)
	#f.close()
	
	#exit()
	#if options.append>0:
		#convnet.appendlayer(options.append)
	learning_rate = T.scalar('lr')  # learning rate to use
	weight_decay = T.scalar('wd')  # learning rate to use
	#layer0_input = x.reshape((batch_size, 1,  shape[0],shape[1]))
	cost, updates = convnet.get_cost_updates(learning_rate,weight_decay)
	train_model = theano.function(
				inputs=[
					index,
					theano.Param(learning_rate, default=0.1),
					theano.Param(weight_decay, default=1e-5)
				],
				outputs=cost,
				updates=updates,
				givens={
					convnet.x: train_set_x[index * batch_size: (index+1) * batch_size]
				}
			)
	
	
	print '... pre-training the model'
	### Pre-train layer-wise
	#if options.layer==None:
		#totrain=range(convnet.n_layers)
	#else:
		#totrain=[int(i) for i in options.layer.split(',')]
	totrain=[0]
	for i in totrain:
		learning_rate=options.learnrate

		n_train_batches = train_set_x.get_value(borrow=True).shape[0] / batch_size
		# go through pretraining epochs
		#training_epochs=eplist[i]
		for epoch in xrange(training_epochs):
		# go through the training set
			c = []
			for batch_index in xrange(n_train_batches):
				err=train_model(batch_index,
					lr=learning_rate,
					wd=options.weightdecay)
				c.append(err)
				#print err
			learning_rate*=.999
			#print test_imgs[convnet.n_layers-1](index=0)
			print 'Pre-training layer %i, epoch %d, cost ' % (i, epoch),
			print np.mean(c)#,", learning rate",learning_rate

	f = file(options.pretrainnet, 'wb')
	cPickle.dump(convnet, f, protocol=cPickle.HIGHEST_PROTOCOL)
	f.close()
	
	
	if options.checkmid:
		print "checking mid layer ..."
		
		test_imgs = theano.function(
			inputs=[index],
			outputs=convnet.hidden,
			givens={
				convnet.x: train_set_x[index * batch_size: (index+1) * batch_size]
			}
		)
		for ni in [0]:
			fname="midout_conv{:d}.hdf".format(ni)
			try:os.remove(fname)
			except: pass
			for idi in range(1):
				rt=test_imgs(0)[0]
				print np.shape(rt)

				ipt= train_set_x[idi * batch_size: (idi + 1) * batch_size]
				ipt= ipt.eval()
				for t in range(len(rt)):
					img=rt[t]
					e = EMNumPy.numpy2em(img.astype("float32"))
					e.write_image(fname,-1)
	
	if options.sdaout:
		print "Generating results ..."
		
		test_imgs = theano.function(
			inputs=[index],
			outputs=convnet.get_reconstructed_input(),
			givens={
				convnet.x: train_set_x[index * batch_size: (index+1) * batch_size]
			}
		)
		for ni in [0]:
			fname="result_conv{:d}.hdf".format(ni)
			try:os.remove(fname)
			except: pass
			for idi in range(1):
				rt=test_imgs(0)
				print np.shape(test_imgs(0))

				ipt= train_set_x[idi * batch_size: (idi + 1) * batch_size]
				ipt= ipt.eval()
				for t in range(len(rt)):
					#img=ipt[t].reshape(lth,lth)
					
					if shape[2]>1:
						img=ipt[t].reshape(shape[0],shape[1],shape[2])
					else:
						img=ipt[t].reshape(shape[0],shape[1])
					e = EMNumPy.numpy2em(img.astype("float32"))
					e.write_image(fname,-1)
					#print len(rt[t]),len(ipt[t])
					
					if shape[2]>1:
						img=rt[t].reshape(shape[0],shape[1],shape[2])
					else:
						img=rt[t].reshape(shape[0],shape[1])
					e = EMNumPy.numpy2em(img.astype("float32"))
					e.write_image(fname,-1)

def load_particles(ptcls,shrink=3,mask=12):
	num=EMUtil.get_image_count(ptcls)
	mx=0
	mn=1e10
	data=[]
	for i in range(num):
		ptl=EMData(ptcls,i)
		if ptl["maximum"]>mx:
			mx=ptl["maximum"]
		if ptl["minimum"]<mn:
			mn=ptl["minimum"]
			
	for i in range(num):
		ptl=EMData(ptcls,i)
		if shrink>1:
			ptl.process_inplace("math.meanshrink",{"n":shrink})

		ptl.div(mx)
		ar=EMNumPy.em2numpy(ptl)
		shp=np.shape(ar)
		##print shp
		#ar=ar.reshape(1,shp[0],shp[1])
		#print ar
		data.append(ar.flatten())
	
	shared_x = theano.shared(np.asarray(data,dtype=theano.config.floatX),borrow=True)
	
	return shared_x

class LeNetConvPoolLayer(object):

	def __init__(self, rng, filter_shape, image_shape, poolsize=(2, 2)):
		
		assert image_shape[1] == filter_shape[1]
		self.image_shape=theano.shared(
			value=np.asarray(image_shape,dtype='int16'),borrow=True)
		self.poolsize=poolsize
		#self.input = input
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
		pooled_out = downsample.max_pool_2d(
			input=conv_out,
			ds=poolsize,
			ignore_border=True
		)

		# add the bias term. Since the bias is a vector (1D array), we first
		# reshape it to a tensor of shape (1, n_filters, 1, 1). Each bias will
		# thus be broadcasted across mini-batches and feature map
		# width & height
		self.hidden = T.tanh(pooled_out + self.b.dimshuffle('x', 0, 'x', 'x'))

		# store parameters of this layer
		self.params = [self.W, self.b]
	
	
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

		stretch_unpooling_out = T.nnet.neighbours.neibs2images(stacked_conv_neibs, self.poolsize, self.x1.shape)
		#return self.hidden
		#z=T.tanh(stretch_unpooling_out + self.b_prime.dimshuffle('x', 0, 'x', 'x'))
		#return T.sum((self.x-z)**2,axis=3)
		#rectified_linear_activation = lambda x: T.maximum(0.0, x)
		return T.tanh(stretch_unpooling_out + self.b_prime.dimshuffle('x', 0, 'x', 'x'))

	def get_cost_updates(self, learning_rate,weight_decay):
		
			""" This function computes the cost and the updates for one trainng
			step of the dA """

			z = self.get_reconstructed_input()
			# note : we sum over the size of a datapoint; if we are using
			#	minibatches, L will be a vector, with one entry per
			#		example in minibatch
			#L = - T.sum(self.x * T.log(z) + (1 - self.x) * T.log(1 - z), axis=1)
			L = T.sum((self.x1-z)*(self.x1-z),axis=3)
			# note : L is now a vector, where each element is the
			#	cross-entropy cost of the reconstruction of the
			#	corresponding example of the minibatch. We need to
			#	compute the average of all these to get the cost of
			#	the minibatch
			cost = T.mean(L)
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
