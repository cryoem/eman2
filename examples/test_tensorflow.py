#!/usr/bin/env python
# Muyuan Chen 2019-04

import numpy as np
from EMAN3 import *
from EMAN3tensor import *
import tensorflow as tf
from time import time


def main():
	
	usage=" "
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--gpuid", type=int,help="gpu to use, -1 for CPU", default=0)
	(options, args) = parser.parse_args()
#	logid=E2init(sys.argv)
	
	dev=tf_set_device(options.gpuid,8192)

	with dev:
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

		start=time()
		optimizer=tf.keras.optimizers.Adam(0.05)
		#train_step=optimizer.minimize(loss, var_list=tbp)
		#session.run(tf.global_variables_initializer())
		for i in range(100):
			with tf.GradientTape() as gt:
				tcp=tf.nn.conv2d(ta, tbp, strides=[1, 1, 1, 1], padding='SAME')
				loss=tf.reduce_mean((tcp-tc)**2)

			grad=gt.gradient(loss, [tbp])
			optimizer.apply_gradients(zip(grad, [tbp]))
			if i%5==0:
				gd=np.mean(abs(grad[0][0]))
				print("  iter {}, loss {:.2f}, mean grad {:.2f}".format(i, loss, gd))


		print("Truth: {}\tEstimate:{}\nTime:{:0.2f} s\n".format(b[:,:,0,0],tbp[:,:,0,0],time()-start))

		print("Testing FFTs")
		imgs=[test_image(size=(1024,1024)) for i in range(256)]
		imgstf=to_tf(imgs)
		print(f"{imgstf.shape} Allocated")
		start=time()
		if options.gpuid==-1:
			ffts=tf_fft2d(imgstf)
			print(f"Done: {(time()-start)/.002048:0.3f} us/fft (one core)")
		else:
			for i in range(1000): ffts=tf_fft2d(imgstf)
			print(f"Done: {(time()-start)/2.048:0.3f} us/fft")

		imgs=[test_image(size=(256,256)) for i in range(4096)]
		imgstf=to_tf(imgs)
		print(f"{imgstf.shape} Allocated")
		start=time()
		if options.gpuid==-1:
			ffts=tf_fft2d(imgstf)
			print(f"Done: {(time()-start)/.004096:0.3f} us/fft (one core)")
		else:
			for i in range(1000): ffts=tf_fft2d(imgstf)
			print(f"Done: {(time()-start)/4.096:0.3f} us/fft")

		imgs=[test_image(size=(256,256)) for i in range(1024)]
		imgstf=to_tf(imgs)
		print(f"{imgstf.shape} Allocated")
		start=time()
		if options.gpuid==-1:
			ffts=tf_fft2d(imgstf)
			print(f"Done: {(time()-start)/.001024:0.3f} us/fft (one core)")
		else:
			for i in range(1000): ffts=tf_fft2d(imgstf)
			print(f"Done: {(time()-start)/1.024:0.3f} us/fft")

		imgs=[test_image(size=(256,256)) for i in range(512)]
		imgstf=to_tf(imgs)
		print(f"{imgstf.shape} Allocated")
		start=time()
		if options.gpuid==-1:
			ffts=tf_fft2d(imgstf)
			print(f"Done: {(time()-start)/.000512:0.3f} us/fft (one core)")
		else:
			for i in range(1000): ffts=tf_fft2d(imgstf)
			print(f"Done: {(time()-start)/0.512:0.3f} us/fft")


#	E2end(logid)
	
	
if __name__ == '__main__':
	main()
	

