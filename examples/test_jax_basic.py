#!/usr/bin/env python

# This program doesn't have any EMAN dependencies, does some basic tensorflow testing

import numpy as np
import jax
import jax.numpy as jnp
import time

print("This performs a very basic test to insure JAX is installed and running. If you see 'test successful', then JAX should be installed correctly. Otherwise, a careful reading of the error may help.")

# make a base image (noise) and a then 1024 unique noise images. Add the base to each unique one
imgnoise=jnp.array(np.random.normal(0,1,(1024,256,256)))
imgbase=jnp.array(np.random.normal(0,1,(256,256)))
imgnoise+=imgbase		# 1024 images all with a common base noise with additional random noise

discard=jnp.fft.rfft2(imgbase)	# make sure everything is loaded, so timing is consistent

# Time to compute translational alignment by cross-correlation of 1024 256x256 images (should all have 0 translation, so a peak at the origin)
t0=time.time()
imgbasefc=jnp.conj(jnp.fft.rfft2(imgbase))	# complex conjugate of the base tensor
imgnoisef=jnp.fft.rfft2(imgnoise)				# ffts of the individual noisy images
ccfsf=imgnoisef*imgbasefc							# cross correlation in Fourier space
ccfs=jnp.fft.irfft2(ccfsf)						# in real space
t1=time.time()
print(t1-t0)
print("test successful")
