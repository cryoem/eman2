#!/usr/bin/env python
# Muyuan June 2015
import os
import sys
import time

import random
import numpy as np
from EMAN2 import *
import theano
import theano.tensor as T
from theano.tensor.shared_randomstreams import RandomStreams
import cPickle

def main():


	usage="""datest.py <particles> [options]
	Run a stacked denoise auto-encoder on the data set."""
	
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--learnrate", type=float,help="Learning rate for the auto-encoder", default=.01)
	parser.add_argument("--niter", type=int,help="Training iterations", default=10)
	parser.add_argument("--weightdecay", type=float,help="Weight decay. Used for regularization.", default=1e-5)
	parser.add_argument("--batch", type=int,help="Batch size for the stochastic gradient descent.", default=50)
	parser.add_argument("--shrink", type=int,help="Shrink the image.", default=1)
	parser.add_argument("--pretrainnet", type=str,help="Output pretrained autoencoder file name", default="da.save")
	parser.add_argument("--hidden", type=str,help="Size of hidden layers: n1,n2,n3...", default="1000,500")
	parser.add_argument("--sdaout", action="store_true", default=False ,help="Output the result of the autoencoder")
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
	train_set_x= particles
	
	ltrain=train_set_x.get_value().shape[1]
	#lth=int(sqrt(ltrain))
	
	
	ptl=EMData(args[0],0)
	shape=[ptl["nx"],ptl["ny"],ptl["nz"]]
	#print shape
	# compute number of minibatches 
	

	# allocate symbolic variables for the data
	index = T.lscalar()	# index to a [mini]batch
	x = T.matrix('x')  # the data is presented as rasterized images
	
	
	print "setting up model"
	rng = np.random.RandomState(123)
	
	if options.fromlast:
		print "loading {}...".format(options.pretrainnet)
		f = file(options.pretrainnet, 'rb')
		sda = cPickle.load(f)
		f.close()		
		x=sda.x
	
	else:
		sda = SdA(
			numpy_rng=rng,
			n_ins=ltrain,
			hidden_layers_sizes=hidden
		)
	pretraining_fns = sda.pretraining_functions(train_set_x=train_set_x,batch_size=batch_size)
	
	#test_imgs = sda.pretraining_get_result(train_set_x=train_set_x,batch_size=1)
	
	
	print '... pre-training the model'
	### Pre-train layer-wise
	for i in xrange(sda.n_layers):
		learning_rate=options.learnrate

		n_train_batches = train_set_x.get_value(borrow=True).shape[0] / batch_size
		# go through pretraining epochs
		for epoch in xrange(training_epochs):
		# go through the training set
			c = []
			for batch_index in xrange(n_train_batches):
				err=pretraining_fns[i](index=batch_index,
					lr=learning_rate,
					wd=options.weightdecay)
				c.append(err)
				#print err
			learning_rate*=.999
			#print test_imgs[sda.n_layers-1](index=0)
			print 'Pre-training layer %i, epoch %d, cost ' % (i, epoch),
			print np.mean(c),", learning rate",learning_rate

			
	
	### testing..
	if options.sdaout:
		print "Generating results ..."
		
		if options.writeall:
			n_imgs= train_set_x.get_value(borrow=True).shape[0] 
			test_imgs = sda.pretraining_get_result(train_set_x=train_set_x,batch_size=1)
			
			fname=options.writeall
			try:os.remove(fname)
			except: pass
			for i in range(n_imgs):
				rt=test_imgs[sda.n_layers-1](index=i)[0]
			
				if shape[2]>1:
					img=rt.reshape(shape[0],shape[1],shape[2])
				else:
					img=rt.reshape(shape[0],shape[1])
				e = EMNumPy.numpy2em(img.astype("float32"))
				e.process_inplace("normalize")
				e.process_inplace("filter.highpass.gauss",{"cutoff_abs":.02})
				e.write_image(fname,i)
				hdr=EMData(args[0],i,True)
				hdr.write_image(fname,i)
		else:
			batch_size=100
			test_imgs = sda.pretraining_get_result(train_set_x=train_set_x,batch_size=batch_size)
			
			for ni in xrange(sda.n_layers):
				fname="result_sda{:d}.hdf".format(ni)
				try:os.remove(fname)
				except: pass
				for idi in range(1):
					rt=test_imgs[ni](index=idi)

					ipt= train_set_x[idi * batch_size: (idi + 1) * batch_size]
					ipt= ipt.eval()
					for t in range(len(rt)):
						#img=ipt[t].reshape(lth,lth)
						#print np.shape(rt[t])
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
		
		
		f = file(options.pretrainnet, 'wb')
		cPickle.dump(sda, f, protocol=cPickle.HIGHEST_PROTOCOL)
		f.close()

	### save the weights (first layer only)
	if options.weights!=None:
		try:os.remove(options.weights)
		except: pass
		wt=sda.dA_layers[0].W.get_value(borrow=True).T
		for w in wt:
			img=w.reshape(lth,lth)
			img=img.astype("float32")
			em = EMNumPy.numpy2em(img)
			em.write_image(options.weights,-1)
	
	E2progress(E2n,1.0)
	E2end(E2n)
	

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
		ar=ar.flatten()
		data.append(ar)
	
	shared_x = theano.shared(np.asarray(data,dtype=theano.config.floatX),borrow=True)
	return shared_x
	

###############################################################
###  deep learning classes borrowed from the tutorial
###########################################################


class dA(object):
	"""Denoising Auto-Encoder class (dA)

	A denoising autoencoders tries to reconstruct the input from a corrupted
	version of it by projecting it first in a latent space and reprojecting
	it afterwards back in the input space. Please refer to Vincent et al.,2008
	for more details. If x is the input then equation (1) computes a partially
	destroyed version of x by means of a stochastic mapping q_D. Equation (2)
	computes the projection of the input into the latent space. Equation (3)
	computes the reconstruction of the input, while equation (4) computes the
	reconstruction error.

	.. math::

		\tilde{x} ~ q_D(\tilde{x}|x)									 (1)

		y = s(W \tilde{x} + b)										   (2)

		x = s(W' y  + b')												(3)

		L(x,z) = -sum_{k=1}^d [x_k \log z_k + (1-x_k) \log( 1-z_k)]	  (4)

	"""

	def __init__(
		self,
		numpy_rng,
		theano_rng=None,
		input=None,
		n_visible=784,
		n_hidden=500,
		W=None,
		bhid=None,
		bvis=None
	):
		"""
		Initialize the dA class by specifying the number of visible units (the
		dimension d of the input ), the number of hidden units ( the dimension
		d' of the latent or hidden space ) and the corruption level. The
		constructor also receives symbolic variables for the input, weights and
		bias. Such a symbolic variables are useful when, for example the input
		is the result of some computations, or when weights are shared between
		the dA and an MLP layer. When dealing with SdAs this always happens,
		the dA on layer 2 gets as input the output of the dA on layer 1,
		and the weights of the dA are used in the second stage of training
		to construct an MLP.

		:type numpy_rng: np.random.RandomState
		:param numpy_rng: number random generator used to generate weights

		:type theano_rng: theano.tensor.shared_randomstreams.RandomStreams
		:param theano_rng: Theano random generator; if None is given one is
					 generated based on a seed drawn from `rng`

		:type input: theano.tensor.TensorType
		:param input: a symbolic description of the input or None for
					  standalone dA

		:type n_visible: int
		:param n_visible: number of visible units

		:type n_hidden: int
		:param n_hidden:  number of hidden units

		:type W: theano.tensor.TensorType
		:param W: Theano variable pointing to a set of weights that should be
				  shared belong the dA and another architecture; if dA should
				  be standalone set this to None

		:type bhid: theano.tensor.TensorType
		:param bhid: Theano variable pointing to a set of biases values (for
					 hidden units) that should be shared belong dA and another
					 architecture; if dA should be standalone set this to None

		:type bvis: theano.tensor.TensorType
		:param bvis: Theano variable pointing to a set of biases values (for
					 visible units) that should be shared belong dA and another
					 architecture; if dA should be standalone set this to None


		"""
		self.n_visible = n_visible
		self.n_hidden = n_hidden

		# create a Theano random generator that gives symbolic random values
		# note : W' was written as `W_prime` and b' as `b_prime`
		if not W:
			# W is initialized with `initial_W` which is uniformely sampled
			# from -4*sqrt(6./(n_visible+n_hidden)) and
			# 4*sqrt(6./(n_hidden+n_visible))the output of uniform if
			# converted using asarray to dtype
			# theano.config.floatX so that the code is runable on GPU
			#print np.sqrt(6. / (n_hidden + n_visible))
			initial_W = np.asarray(
				numpy_rng.uniform(
					low=-1*np.sqrt(6. / (n_hidden + n_visible)) ,
					high=1*np.sqrt(6. / (n_hidden + n_visible)),
					size=(n_visible, n_hidden)
				),
				dtype=theano.config.floatX
			)
			W = theano.shared(value=initial_W, name='W', borrow=True)
			#initial_W2 = np.asarray(
				#numpy_rng.uniform(
					#low=-.4 * np.sqrt(6. / (n_hidden + n_visible)),
					#high=.4 * np.sqrt(6. / (n_hidden + n_visible)),
					#size=(n_hidden,n_visible)
				#),
				#dtype=theano.config.floatX
			#)
			#Wt = theano.shared(value=initial_W2, name='Wt', borrow=True)

		if not bvis:
			bvis = theano.shared(
				value=np.zeros(
					n_visible,
					dtype=theano.config.floatX
				),
				borrow=True
			)

		if not bhid:
			bhid = theano.shared(
				value=np.zeros(
					n_hidden,
					dtype=theano.config.floatX
				),
				name='b',
				borrow=True
			)

		self.W = W
		# b corresponds to the bias of the hidden
		self.b = bhid
		# b_prime corresponds to the bias of the visible
		self.b_prime = bvis
		# tied weights, therefore W_prime is W transpose
		self.W_prime = W.T
		#self.W_prime = Wt
		# if no input is given, generate a variable representing the input
		if input is None:
			# we use a matrix because we expect a minibatch of several
			# examples, each example being a row
			self.x = T.dmatrix(name='input')
		else:
			self.x = input

		self.params = [self.W, self.b, self.b_prime]


	def get_hidden_values(self, input):
		""" Computes the values of the hidden layer """
		#return T.nnet.sigmoid(T.dot(input, self.W) + self.b)
		return T.tanh(T.dot(input, self.W) + self.b)

	def get_reconstructed_input(self, hidden):
		"""Computes the reconstructed input given the values of the
		hidden layer

		"""
		#return T.nnet.sigmoid(T.dot(hidden, self.W_prime) + self.b_prime)
		return T.tanh(T.dot(hidden, self.W_prime) + self.b_prime)
	
	def get_cost_updates(self, learning_rate,weight_decay):
		""" This function computes the cost and the updates for one trainng
		step of the dA """

		y = self.get_hidden_values(self.x)
		z = self.get_reconstructed_input(y)
		# note : we sum over the size of a datapoint; if we are using
		#	minibatches, L will be a vector, with one entry per
		#		example in minibatch
		#L = - T.sum(self.x * T.log(z) + (1 - self.x) * T.log(1 - z), axis=1)
		L = T.sum((self.x-z)*(self.x-z),axis=1)
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
	
	def get_result(self):
		""" This function computes the cost and the updates for one trainng
		step of the dA """

		y = self.get_hidden_values(self.x)
		z = self.get_reconstructed_input(y)
		
		return z

# start-snippet-1
class HiddenLayer(object):
	def __init__(self, rng, input, n_in, n_out, W=None, b=None,
				 activation=T.tanh):
		"""
		Typical hidden layer of a MLP: units are fully-connected and have
		sigmoidal activation function. Weight matrix W is of shape (n_in,n_out)
		and the bias vector b is of shape (n_out,).

		NOTE : The nonlinearity used here is tanh

		Hidden unit activation is given by: tanh(dot(input,W) + b)

		:type rng: np.random.RandomState
		:param rng: a random number generator used to initialize weights

		:type input: theano.tensor.dmatrix
		:param input: a symbolic tensor of shape (n_examples, n_in)

		:type n_in: int
		:param n_in: dimensionality of input

		:type n_out: int
		:param n_out: number of hidden units

		:type activation: theano.Op or function
		:param activation: Non linearity to be applied in the hidden
						   layer
		"""
		self.input = input
		# end-snippet-1

		# `W` is initialized with `W_values` which is uniformely sampled
		# from sqrt(-6./(n_in+n_hidden)) and sqrt(6./(n_in+n_hidden))
		# for tanh activation function
		# the output of uniform if converted using asarray to dtype
		# theano.config.floatX so that the code is runable on GPU
		# Note : optimal initialization of weights is dependent on the
		#		activation function used (among other things).
		#		For example, results presented in [Xavier10] suggest that you
		#		should use 4 times larger initial weights for sigmoid
		#		compared to tanh
		#		We have no info for other function, so we use the same as
		#		tanh.
		if W is None:
			W_values = np.asarray(
				rng.uniform(
					low=-np.sqrt(6. / (n_in + n_out)),
					high=np.sqrt(6. / (n_in + n_out)),
					size=(n_in, n_out)
				),
				dtype=theano.config.floatX
			)
			if activation == theano.tensor.nnet.sigmoid:
				W_values *= 4

			W = theano.shared(value=W_values, name='W', borrow=True)

		if b is None:
			b_values = np.zeros((n_out,), dtype=theano.config.floatX)
			b = theano.shared(value=b_values, name='b', borrow=True)

		self.W = W
		self.b = b

		lin_output = T.dot(input, self.W) + self.b
		self.output = (
			lin_output if activation is None
			else activation(lin_output)
		)
		# parameters of the model
		self.params = [self.W, self.b]

# start-snippet-1
class SdA(object):
	"""Stacked denoising auto-encoder class (SdA)

	A stacked denoising autoencoder model is obtained by stacking several
	dAs. The hidden layer of the dA at layer `i` becomes the input of
	the dA at layer `i+1`. The first layer dA gets as input the input of
	the SdA, and the hidden layer of the last dA represents the output.
	Note that after pretraining, the SdA is dealt with as a normal MLP,
	the dAs are only used to initialize the weights.
	"""

	def __init__(
		self,
		numpy_rng,
		theano_rng=None,
		n_ins=784,
		hidden_layers_sizes=[500, 500],
		n_outs=10,
		corruption_levels=[0.1, 0.1]
	):
		""" This class is made to support a variable number of layers.

		:type numpy_rng: np.random.RandomState
		:param numpy_rng: numpy random number generator used to draw initial
					weights

		:type theano_rng: theano.tensor.shared_randomstreams.RandomStreams
		:param theano_rng: Theano random generator; if None is given one is
						   generated based on a seed drawn from `rng`

		:type n_ins: int
		:param n_ins: dimension of the input to the sdA

		:type n_layers_sizes: list of ints
		:param n_layers_sizes: intermediate layers size, must contain
							   at least one value

		:type n_outs: int
		:param n_outs: dimension of the output of the network

		:type corruption_levels: list of float
		:param corruption_levels: amount of corruption to use for each
								  layer
		"""

		self.sigmoid_layers = []
		self.dA_layers = []
		self.params = []
		self.n_layers = len(hidden_layers_sizes)

		assert self.n_layers > 0

		if not theano_rng:
			theano_rng = RandomStreams(numpy_rng.randint(2 ** 30))
		# allocate symbolic variables for the data
		self.x = T.matrix('x')  # the data is presented as rasterized images
		self.y = T.ivector('y')  # the labels are presented as 1D vector of
								 # [int] labels
		# end-snippet-1

		# The SdA is an MLP, for which all weights of intermediate layers
		# are shared with a different denoising autoencoders
		# We will first construct the SdA as a deep multilayer perceptron,
		# and when constructing each sigmoidal layer we also construct a
		# denoising autoencoder that shares weights with that layer
		# During pretraining we will train these autoencoders (which will
		# lead to chainging the weights of the MLP as well)
		# During finetunining we will finish training the SdA by doing
		# stochastich gradient descent on the MLP

		# start-snippet-2
		for i in xrange(self.n_layers):
			# construct the sigmoidal layer

			# the size of the input is either the number of hidden units of
			# the layer below or the input size if we are on the first layer
			if i == 0:
				input_size = n_ins
			else:
				input_size = hidden_layers_sizes[i - 1]

			# the input to this layer is either the activation of the hidden
			# layer below or the input of the SdA if you are on the first
			# layer
			if i == 0:
				layer_input = self.x
			else:
				layer_input = self.sigmoid_layers[-1].output

			

			# Construct a denoising autoencoder that shared weights with this
			# layer
			dA_layer = dA(numpy_rng=numpy_rng,
						  theano_rng=theano_rng,
						  input=layer_input,
						  n_visible=input_size,
						  n_hidden=hidden_layers_sizes[i])
			
			sigmoid_layer = HiddenLayer(rng=numpy_rng,
						input=layer_input,
						W=dA_layer.W,
						b=dA_layer.b,
						n_in=input_size,
						n_out=hidden_layers_sizes[i],
						#activation=T.nnet.sigmoid)
						activation=T.tanh)
			# add the layer to our list of layers
			self.sigmoid_layers.append(sigmoid_layer)
			# its arguably a philosophical question...
			# but we are going to only declare that the parameters of the
			# sigmoid_layers are parameters of the StackedDAA
			# the visible biases in the dA are parameters of those
			# dA, but not the SdA
			self.params.extend(sigmoid_layer.params)
			self.dA_layers.append(dA_layer)
	
	def pretraining_get_result(self,train_set_x, batch_size):
		
		index = T.lscalar('index')  # index to a minibatch
		corruption_level = T.scalar('corruption')  # % of corruption to use
		
		# begining of a batch, given `index`
		batch_begin = index * batch_size
		# ending of a batch given `index`
		batch_end = batch_begin + batch_size
		
		ret_imgs=[]
		for nl in xrange(self.n_layers):
			da=self.dA_layers[nl]
			result_da = da.get_result()
			### feed forward for result
			for li in xrange(nl):
				nda=self.dA_layers[nl-1-li]
				result_da=nda.get_reconstructed_input(result_da)
				
			test_da = theano.function(
				[index],
				result_da,
				givens={
					self.x: train_set_x[batch_begin: batch_end]
				}
			)
			ret_imgs.append(test_da)
		
		return ret_imgs
		
	def pretraining_functions(self, train_set_x, batch_size):
		''' Generates a list of functions, each of them implementing one
		step in trainnig the dA corresponding to the layer with same index.
		The function will require as input the minibatch index, and to train
		a dA you just need to iterate, calling the corresponding function on
		all minibatch indexes.

		:type train_set_x: theano.tensor.TensorType
		:param train_set_x: Shared variable that contains all datapoints used
							for training the dA

		:type batch_size: int
		:param batch_size: size of a [mini]batch

		:type learning_rate: float
		:param learning_rate: learning rate used during training for any of
							  the dA layers
		'''

		# index to a [mini]batch
		index = T.lscalar('index')  # index to a minibatch
		learning_rate = T.scalar('lr')  # learning rate to use
		weight_decay = T.scalar('wd')  # learning rate to use
		pretrain_fns = []
		ct=0
		for dA in self.dA_layers:
				
			batch_begin = index * batch_size
			batch_end = batch_begin + batch_size
			ct+=1
			# get the cost and the updates list
			cost, updates = dA.get_cost_updates(learning_rate,weight_decay)
			# compile the theano function
			fn = theano.function(
				inputs=[
					index,
					theano.Param(learning_rate, default=0.1),
					theano.Param(weight_decay, default=1e-5)
				],
				outputs=cost,
				updates=updates,
				givens={
					self.x: train_set_x[batch_begin: batch_end]
				}
			)
			# append `fn` to the list of functions
			pretrain_fns.append(fn)

		return pretrain_fns
	
	
		

class LogisticRegression(object):
	"""Multi-class Logistic Regression Class

	The logistic regression is fully described by a weight matrix :math:`W`
	and bias vector :math:`b`. Classification is done by projecting data
	points onto a set of hyperplanes, the distance to which is used to
	determine a class membership probability.
	"""

	def __init__(self, input, n_in, n_out):
		""" Initialize the parameters of the logistic regression

		:type input: theano.tensor.TensorType
		:param input: symbolic variable that describes the input of the
					  architecture (one minibatch)

		:type n_in: int
		:param n_in: number of input units, the dimension of the space in
					 which the datapoints lie

		:type n_out: int
		:param n_out: number of output units, the dimension of the space in
					  which the labels lie

		"""
		# start-snippet-1
		# initialize with 0 the weights W as a matrix of shape (n_in, n_out)
		#rng = np.random.RandomState(123)
		#W_values = np.asarray(
				#rng.uniform(
					#low=-np.sqrt(1. / (n_in + n_out)),
					#high=np.sqrt(1. / (n_in + n_out)),
					#size=(n_in, n_out)
				#),
				#dtype=theano.config.floatX
			#)
		self.W = theano.shared(
			#value=W_values,
			value=np.zeros(
				(n_in, n_out),
				dtype=theano.config.floatX
			),
			name='W',
			borrow=True
		)
		# initialize the baises b as a vector of n_out 0s
		self.b = theano.shared(
			value=np.zeros(
				(n_out,),
				dtype=theano.config.floatX
			),
			name='b',
			borrow=True
		)

		# symbolic expression for computing the matrix of class-membership
		# probabilities
		# Where:
		# W is a matrix where column-k represent the separation hyper plain for
		# class-k
		# x is a matrix where row-j  represents input training sample-j
		# b is a vector where element-k represent the free parameter of hyper
		# plain-k
		self.p_y_given_x = T.nnet.softmax(T.dot(input, self.W) + self.b)

		# symbolic description of how to compute prediction as class whose
		# probability is maximal
		self.y_pred = T.argmax(self.p_y_given_x, axis=1)
		# end-snippet-1

		# parameters of the model
		self.params = [self.W, self.b]

	def negative_log_likelihood(self, y):
		"""Return the mean of the negative log-likelihood of the prediction
		of this model under a given target distribution.

		.. math::

			\frac{1}{|\mathcal{D}|} \mathcal{L} (\theta=\{W,b\}, \mathcal{D}) =
			\frac{1}{|\mathcal{D}|} \sum_{i=0}^{|\mathcal{D}|}
				\log(P(Y=y^{(i)}|x^{(i)}, W,b)) \\
			\ell (\theta=\{W,b\}, \mathcal{D})

		:type y: theano.tensor.TensorType
		:param y: corresponds to a vector that gives for each example the
				  correct label

		Note: we use the mean instead of the sum so that
			  the learning rate is less dependent on the batch size
		"""
		# start-snippet-2
		# y.shape[0] is (symbolically) the number of rows in y, i.e.,
		# number of examples (call it n) in the minibatch
		# T.arange(y.shape[0]) is a symbolic vector which will contain
		# [0,1,2,... n-1] T.log(self.p_y_given_x) is a matrix of
		# Log-Probabilities (call it LP) with one row per example and
		# one column per class LP[T.arange(y.shape[0]),y] is a vector
		# v containing [LP[0,y[0]], LP[1,y[1]], LP[2,y[2]], ...,
		# LP[n-1,y[n-1]]] and T.mean(LP[T.arange(y.shape[0]),y]) is
		# the mean (across minibatch examples) of the elements in v,
		# i.e., the mean log-likelihood across the minibatch.
		return -T.mean(T.log(self.p_y_given_x)[T.arange(y.shape[0]), y])
		# end-snippet-2

	def errors(self, y):
		"""Return a float representing the number of errors in the minibatch
		over the total number of examples of the minibatch ; zero one
		loss over the size of the minibatch

		:type y: theano.tensor.TensorType
		:param y: corresponds to a vector that gives for each example the
				  correct label
		"""

		# check if y has same dimension of y_pred
		if y.ndim != self.y_pred.ndim:
			raise TypeError(
				'y should have the same shape as self.y_pred',
				('y', y.type, 'y_pred', self.y_pred.type)
			)
		# check if y is of the correct datatype
		if y.dtype.startswith('int'):
			# the T.neq operator returns a vector of 0s and 1s, where 1
			# represents a mistake in prediction
			return T.mean(T.neq(self.y_pred, y))
		else:
			raise NotImplementedError()
		
	def result(self):
		return self.y_pred




if __name__ == '__main__':
	main()
