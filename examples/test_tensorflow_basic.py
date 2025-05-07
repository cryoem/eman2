#!/usr/bin/env python

# This program doesn't have any EMAN dependencies, does some basic tensorflow testing
print("Ignore any warnings you see here:")

import tensorflow as tf
import numpy as np
import time

print("\n\n\nThis performs a very basic test to insure TensorFlow is installed and running. If you see 'test successful', then TensorFlow should be installed correctly. Otherwise, a careful reading of the error may help.")

# make a base image (noise) and a then 1024 unique noise images. Add the base to each unique one
imgnoise=np.random.normal(0,1,(1024,256,256))
imgbase=np.random.normal(0,1,(256,256))
imgnoise+=imgbase		# 1024 images all with a common base noise with additional random noise

imgbase=tf.constant(imgbase)
imgnoise=tf.constant(imgnoise)

discard=tf.signal.rfft2d(imgbase)	# make sure everything is loaded, so timing is consistent

# Time to compute translational alignment by cross-correlation of 1024 256x256 images (should all have 0 translation, so a peak at the origin)
t0=time.time()
imgbasefc=tf.math.conj(tf.signal.rfft2d(imgbase))	# complex conjugate of the base tensor
imgnoisef=tf.signal.rfft2d(imgnoise)				# ffts of the individual noisy images
ccfsf=imgnoisef*imgbasefc							# cross correlation in Fourier space
ccfs=tf.signal.irfft2d(ccfsf)						# in real space
t1=time.time()
print(t1-t0)
print("\ntest successful")
