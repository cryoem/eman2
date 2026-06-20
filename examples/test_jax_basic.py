#!/usr/bin/env python

# This program doesn't have any EMAN dependencies, does some basic tensorflow testing

import numpy as np
import jax
import jax.numpy as jnp
import time

def test_fn(imgbase,imgnoise):
	imgbasefc=jnp.conj(jnp.fft.rfft2(imgbase))	# complex conjugate of the base tensor
	imgnoisef=jnp.fft.rfft2(imgnoise)				# ffts of the individual noisy images
	ccfsf=imgnoisef*imgbasefc							# cross correlation in Fourier space
	ccfs=jnp.fft.irfft2(ccfsf)						# in real space
	return ccfs


print("This performs a very basic test to insure JAX is installed and running. If you see 'test successful', then JAX should be installed correctly. Otherwise, a careful reading of the error may help.")

# make a base image (noise) and a then 8192 unique noise images. Add the base to each unique one
imgnoise=jnp.array(np.random.normal(0,1,(8192,256,256)))
imgbase=jnp.array(np.random.normal(0,1,(256,256)))
imgnoise+=imgbase		# 8192 images all with a common base noise with additional random noise

discard=jnp.fft.rfft2(imgbase)	# make sure everything is loaded on the GPU (hopefully), so timing is consistent

# Time to compute translational alignment by cross-correlation of 8192 256x256 images (should all have 0 translation, so a peak at the origin)
t0=time.process_time()
imgbasefc=jnp.conj(jnp.fft.rfft2(imgbase))	# complex conjugate of the base tensor
imgnoisef=jnp.fft.rfft2(imgnoise)				# ffts of the individual noisy images
ccfsf=imgnoisef*imgbasefc							# cross correlation in Fourier space
ccfs=jnp.fft.irfft2(ccfsf)						# in real space
t1=time.process_time()
print(t1-t0)

print("test JIT compile")
t0=time.process_time()
test_jit=jax.jit(test_fn)
t1=time.process_time()
print(t1-t0)


t0=time.process_time()
ccfs=test_jit(imgbase,imgnoise)
t1=time.process_time()
print(t1-t0)

print("test successful")

