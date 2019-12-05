#!/usr/bin/env python
# Muyuan Chen 2017-10
from __future__ import division
from past.utils import old_div
from builtins import range
import os
os.environ["CUDA_VISIBLE_DEVICES"]="0"
import tensorflow as tf
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3' 
import numpy as np
from EMAN2 import *

def main():
	
	usage=" "
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--goodrefs", type=str,help="good reference", default="info/boxrefs.hdf")
	parser.add_argument("--badrefs", type=str,help="bad reference", default="info/boxrefsbad.hdf")
	parser.add_argument("--bgrefs", type=str,help="background reference", default="info/boxrefsbg.hdf")
	parser.add_argument("--shuffle",action="store_true",help="shuffle dataset at each iteration.",default=False)
	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)

	session=tf.Session()
	
	goodrefs=EMData.read_images(options.goodrefs)
	bgrefs=EMData.read_images(options.bgrefs)
	badrefs=EMData.read_images(options.badrefs)
	
	sz=64
	data, label=load_ptcls(bgrefs, goodrefs, sz, True)
	
	### number of kernel, kernel size, pooling size(2/1)
	kernels=[(20,15,2), (20,15,1), (1,15,1)]
	batchsize=10
	
	nnet0=tf_build_cnn(kernels, sz, batchsize, meanout=False)
	tf_do_training(nnet0, data, label, session, options.shuffle, learnrate=1e-4, niter=10)
	write_output_train(nnet0, 'trainout_pickptcl.hdf', session)
	tf_save_cnn(nnet0, "nnet_pickptcls.hdf", session)
	
	
	data, label=load_ptcls(badrefs, goodrefs, sz, False)
	
	### number of kernel, kernel size, pooling size(2/1)
	kernels=[(20,15,2), (20,15,1), (1,15,1)]
	batchsize=10
	
	nnet1=tf_build_cnn(kernels, sz, batchsize, meanout=True)
	tf_do_training(nnet1, data, label, session, options.shuffle, learnrate=1e-4, niter=15)
	write_output_train(nnet1, 'trainout_classify.hdf', session)
	tf_save_cnn(nnet1, "nnet_classify.hdf", session)
	
	E2end(logid)
	
def load_ptcls(ref0, ref1, sz=64, makegaussian=True):
	print("Pre-processing particles...")
	#### here we shrink the particles so they are 64x64
	#### and duplicate so there are more than 500 good and 500 bad particles

	nref_target=500
	bxsz=ref0[0]["nx"]
	shrinkfac=old_div(float(bxsz),float(sz))

	data=[] ### particles in flattened numpy array
	lbs=[]  ### labels in flattened numpy array

	for label, refs in enumerate([ref0, ref1]):
		nref=len(refs)
		if nref<5:
			print("Not enough references. Please box at least 5 good and 5 background reference...")
		ncopy=old_div(nref_target,nref) + 1

		for pp in refs:
			ptl=pp.process("math.fft.resample",{"n":shrinkfac})
			ptl.clip_inplace(Region(0,0, sz, sz))
	#		 ptl.process_inplace("filter.highpass.gauss",{"cutoff_pixels":3})
	#		 ptl.process_inplace("filter.lowpass.gauss",{"cutoff_abs":0.25})
			ptl.process_inplace("filter.highpass.gauss",{"cutoff_freq":0.005})
			ptl.process_inplace("filter.lowpass.gauss",{"cutoff_freq":0.05})
			for c in range(ncopy):

				tr=Transform()
				tr.set_rotation({"type":"2d","alpha":np.random.random()*360.0})
				img=ptl.process("xform",{"transform":tr})
				ar=img.numpy().copy()
				data.append(ar.flatten())
				lbs.append(label)
				
	rndid=list(range(len(data)))
	np.random.shuffle(rndid)
	data=[data[i] for i in rndid]
	lbs=[lbs[i] for i in rndid]
	data=np.asarray(data,dtype=np.float32)
	div=np.mean(np.std(data,axis=1))
	data/=div#np.std(data)#*2.
	mx=4.
	data[data>mx]=mx
	data[data<-mx]=-mx
	lbs=np.asarray(lbs,dtype=int)
	
	if makegaussian:
		#### make target output
		img=EMData(old_div(sz,2),old_div(sz,2))
		img.process_inplace("testimage.gaussian",{'sigma':5.})
		img.div(img["maximum"])
		gaus=img.numpy().copy().flatten()
		gaus=gaus.astype(np.float32)
		lbarrs=np.array([np.zeros_like(gaus, dtype=np.float32), gaus])
		label_np=lbarrs[lbs]
		return data, label_np
	else:
		return data, lbs

#### network training.
def tf_do_training(convnet, data, label, session, shuffle=False, learnrate=1e-4, niter=10):
	train_step = tf.train.AdamOptimizer(learnrate).minimize(convnet.loss)
	session.run(tf.global_variables_initializer())
	
	for it in range(niter):
		if shuffle:
			convnet.dataset.shuffle(10000)
		session.run(convnet.iterator.initializer,
			feed_dict={convnet.data: data, convnet.label: label})
		cost=[]
		while True:
			try:
				cost.append(session.run((train_step, convnet.loss))[1])
			except tf.errors.OutOfRangeError:
				break
				
		print("iteration {}, cost {:.3f}".format(it, np.mean(cost)))
	
	session.run(convnet.iterator.initializer, feed_dict={convnet.data: data, convnet.label: label})
	
#### save network to hdf file
def tf_save_cnn(convnet, fname, session):
	try: os.remove(fname)
	except: pass
	print("Saving the trained net to {}...".format(fname))
	weights=session.run(convnet.wts)
	bias=session.run(convnet.bss)
	wsz=int(weights[0].shape[0])

	hdr=EMData(wsz,wsz)
	hdr["nkernel"]=[k[0] for k in convnet.kernels]
	hdr["ksize"]=[k[1] for k in convnet.kernels]
	hdr["poolsz"]=[k[2] for k in convnet.kernels]
	hdr["imageshape"]=(10,1,64,64)
	hdr["amplitude"]=1
	nlayer=len(convnet.kernels)



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
			e=from_numpy(wi)
			e.write_image(fname,k)
			k+=1

def write_output_train(convnet,outfile, session, ncopy=10):
	sz=convnet.imgsz
	outsz=convnet.outsz
	try: os.remove(outfile)
	except: pass
	for nc in range(ncopy):
		outx, outy =session.run((convnet.batch_in, convnet.batch_out))
		for i in range(len(outx)):
			ox=outx[i].reshape((sz,sz))
			oy=outy[i].reshape((outsz, outsz))

			e0=from_numpy(ox.copy())
			e0.process_inplace("normalize")
			e0.write_image(outfile, -1)

			e1=from_numpy(oy.copy())
			e1=e1.get_clip(Region(old_div(-(sz-outsz),2),old_div(-(sz-outsz),2),sz,sz))
			e1.scale(old_div(sz,outsz))

			e1.write_image(outfile, -1)

	
def tf_build_cnn(kernels, imgsz=64, batchsz=10, meanout=False):
	
	convnet = type('convnet', (), {})() ### an empty object
	nlayer=len(kernels)
	outsz=(old_div(imgsz,np.prod([k[2] for k in kernels])))
	tf_data = tf.placeholder(tf.float32, shape=[None, imgsz*imgsz], name="tfdata")
	if meanout:
		tf_label = tf.placeholder(tf.float32, shape=[None], name="tflabel")
	else:
		tf_label = tf.placeholder(tf.float32, shape=[None, outsz**2], name="tflabel")
	dataset = tf.contrib.data.Dataset.from_tensor_slices((tf_data, tf_label))
	dataset_batch=dataset.batch(batchsz)
	iterator = dataset_batch.make_initializable_iterator()
	data_inp, data_tar=iterator.get_next()
	
	k0=1
	wts=[]
	bss=[]
	layerout=[tf.reshape(data_inp, (-1, imgsz,imgsz, 1), name="reshape_input")]
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
		y=tf.reshape(layerout[-1], (-1, outsz**2), name="reshape_output")
		yout=tf.maximum(-1.,(tf.reduce_mean(y, axis=1)))
		yout=tf.minimum(1., yout)
		loss=tf.reduce_mean((data_tar-yout)**2)
	else:
		y=tf.reshape(tf.minimum(1.,layerout[-1]), (-1, outsz**2), name="reshape_output")
		loss=tf.log(tf.reduce_mean((data_tar-y)**2))
		
		
	convnet.dataset=dataset
	convnet.kernels=kernels
	convnet.batchsize=batchsz
	convnet.wts=wts
	convnet.bss=bss
	convnet.outsz=outsz
	convnet.imgsz=imgsz
	convnet.layerout=layerout
	convnet.data=tf_data
	convnet.label=tf_label
	
	convnet.batch_in=data_inp
	convnet.batch_out=y
	convnet.batch_tar=data_tar
	
	convnet.loss=loss
	convnet.iterator=iterator
	
	return convnet
	
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

	
def run(cmd):
	print(cmd)
	launch_childprocess(cmd)
	
	
if __name__ == '__main__':
	main()

