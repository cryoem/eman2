#!/usr/bin/env python
# Muyuan Chen 2019-04

import numpy as np
from EMAN2 import *

def main():
	
	usage=" "
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--gpuid", type=str,help="gpu to use", default=None)
	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)
	
	if options.gpuid:
		os.environ["CUDA_VISIBLE_DEVICES"]=options.gpuid
		
	import tensorflow as tf
	
	print("Testing basic operations...")
	#session=tf.Session(); 
	a=tf.Variable(1); 
	#session.run(tf.global_variables_initializer()); 
	print("\t1 + 1 = {}".format(a+a))


	print("Testing matrix multiplication...")
	n=3
	a=np.random.randint(0,5,(n,n)).astype(np.float32)
	b=np.random.randint(0,5,(n,n)).astype(np.float32)
	ta=tf.Variable(a)
	tb=tf.Variable(b)
	#session.run(tf.global_variables_initializer())
	tc=tf.tensordot(ta, tb, 1)
	print(tc)
		
		
	print("Testing convolution...")
	n1=4
	n2=2
	a=np.random.randint(0,5,(1,n1,n1,1)).astype(np.float32)
	b=np.random.randint(0,5,(n2,n2,1,1)).astype(np.float32)
	ta=tf.Variable(a)
	tb=tf.Variable(b)
	#session.run(tf.global_variables_initializer())
	tc=tf.nn.conv2d(ta, tb, strides=[1, 1, 1, 1], padding='SAME')
	print(tc[0,:,:,0])



	print("Testing training set...")
	a=np.random.randn(128, 3).astype(np.float32)
	trainset=tf.data.Dataset.from_tensor_slices((a))
	trainset=trainset.batch(32)
	n=0
	for t in trainset:
		print(tf.reduce_mean(t,axis=0))
		n=n+1
	
	print(n)


	print("Testing training...")
	bp=b+np.random.randn(n2,n2,1,1).astype(np.float32)
	tbp=tf.Variable(bp)
	tcp=tf.nn.conv2d(ta, tbp, strides=[1, 1, 1, 1], padding='SAME')
	loss=tf.reduce_mean((tcp-tc)**2)

	optimizer=tf.keras.optimizers.Adam(0.05)
	#train_step=optimizer.minimize(loss, var_list=tbp)
	#session.run(tf.global_variables_initializer())
	for i in range(51):
		with tf.GradientTape() as gt:
			tcp=tf.nn.conv2d(ta, tbp, strides=[1, 1, 1, 1], padding='SAME')
			loss=tf.reduce_mean((tcp-tc)**2)
			
		grad=gt.gradient(loss, [tbp])
		optimizer.apply_gradients(zip(grad, [tbp]))
		if i%5==0:
			gd=np.mean(abs(grad[0][0]))
			print("  iter {}, loss {:.2f}, mean grad {:.2f}".format(i, loss, gd))
		
		
	print("Truth:\n{}\nEstimate:\n{}".format(b[:,:,0,0],tbp[:,:,0,0]))
	
	
	E2end(logid)
	
	
if __name__ == '__main__':
	main()
	

