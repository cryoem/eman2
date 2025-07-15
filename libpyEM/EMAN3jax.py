#!/usr/bin/env python
#
# Author: Steven Ludtke, 09/10/2024 (sludtke@bcm.edu)
# Copyright (c) 2000-2024 Baylor College of Medicine
#
# This software is issued under a joint BSD/GPL license. You may use the
# source code in this file under either license. However, note that the
# complete EMAN2 and SPARX software packages have some GPL dependencies,
# so you are responsible for compliance with the licenses of these packages
# if you opt to use BSD licensing. The warranty disclaimer below holds
# in either instance.
#
# This complete copyright notice must be included in any revised version of the
# source code. Additional authorship citations may be added, but existing
# author citations must be preserved.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston MA 02111-1307 USA
#

"""
This module is exclusive with EMAN3tensor. Do not try and import both in a single program.

ONLY import this file if you will be working with JAX in your program, otherwise the JAX initialization may add unreasonable startup delays

There are several key classes for data representation in this module:
EMStack3D, 2D, 1D - A set of 3 classes to represent stacks of images of different dimensionality with seamless interconversion among EMData, NumPy and Tensorflow.
	Implemented as 3 separate classes to avoid validation overhead and provide dimensionality-specific routines. All share a common generic interface.

Orientations - an array {N,X,Y,Z} of orientations, interconvertable to EMAN Transform objects. The main representation is an XYZ vector where the length of the vector represents
	the amount of rotation, and the direction is the plane of rotation. This representation is particularly conventient for work with deep learning.

Gaussians - an array of Gaussian objects {N,X,Y,Z,A} with amplitude, but no width

VERY important to note that when indexing EMData objects it is emd[x,y,z], whereas indexing numpy/tensorflow objects, the last index is the fastest varying ary[z,y,x] !

====
tf.signal.rfft2d    - returns a new constant tensor which is the FFT of the last 2 indices of the tensor.
                       ie - if used on a 3-D tensor (N,X,Y) will return (N,FX,FY)
                       FFT tensors are complex, and padded in X, ie (NX,NY) -> (NX/2+1,NY)

tf.signal.irfft2d   - the inverse operation

np.fromfunction(lambda x,y: np.hypot(x,y),(nx,ny)) - for example

"""

from EMAN3 import *
import numpy as np
import jax
import jax.numpy as jnp
from jax import grad, jit
from jax import lax
from jax import random
import jaxlib

# TODO
class StackCache():
	"""This object serves as a cache of EMStack objects which can be conveniently read back in. This provides
	methods for easy/efficient sampling of subsets of data sets which may be too large for RAM. Caches are not persistent across sessions
	Writing images to the cache will coerce them to tensorflow representation."""

	def __init__(self,filename,n):
		"""Specify filename and number of images to be cached."""
		print("cache: ",filename)
		self.filename=filename

		self.fp=open(filename,"wb+")		# erase file!
		self.locs=np.zeros(n+1,dtype=np.int64)			# list of seek locations in the binary file for each image
		self.orts=np.zeros((n,3),dtype=np.float32)		# orientations in spinvec format
		self.tytx=np.zeros((n,3),dtype=np.float32)		# image shifts in absolute [-0.5,0.5] format, third column used for per particle defocus
		self.frcs=np.zeros((n),dtype=np.float32)+2.0	# FRCs from previous round, initialize to 2.0 (> 1.0 normal max)
		self.cloc=0
		self.locked=False

	def __del__(self):
		"""free all resources if possible"""
		self.fp=None
		os.unlink(self.filename)
		self.locs=None

	def write(self,stack,n0,ortss=None,tytxs=None):
		"""writes stack of images starting at n0 to cache"""
		while self.locked: time.sleep(0.1)
		self.locked=True
		stack.coerce_numpy()
		if ortss is not None:
			try: self.orts[n0:n0+len(stack)]=ortss
			except: self.orts[n0:n0+len(stack)]=np.array(ortss)
		if tytxs is not None:
			try: self.tytx[n0:n0+len(stack)]=tytxs
			except:
#				print(tytxs,tytxs.shape)
				self.tytx[n0:n0+len(stack)]=np.array(tytxs)

		# we go through the images one at a time, serialze, and write to a file with a directory
		self.fp.seek(self.cloc)
		for i in range(len(stack)):
			im=stack[i]
			self.locs[n0+i]=self.cloc
			np.save(self.fp,im)
			# self.fp.write(tf.io.serialize_tensor(im).numpy())	#TODO
			self.cloc=self.fp.tell()
			self.locs[n0+i+1]=self.cloc

		self.locked=False

	def add_orts(self,nlist,dorts=None,dtytxs=None):
		"""adds dorts and dtytxs to existing arrays at locations described by nlist, used to take a gradient step
		on a subset of the data."""
		if dorts is not None: self.orts[nlist]+=dorts
		if dtytxs is not None:
			try: self.tytx[nlist]+=dtytxs
			except: self.tytx[nlist,:2]+=dtytxs

	def set_frcs(self,nlist,frcs):
		self.frcs[nlist]=frcs


	def read(self,nlist):
		while self.locked: time.sleep(0.1)

		self.locked=True

		stack=[]
		for i in nlist:
			try:
				self.fp.seek(self.locs[i])
				stack.append(np.load(self.fp))
#				stack.append(tf.io.parse_tensor(self.fp.read(self.locs[i+1]-self.locs[i]),out_type=tf.complex64))  # TODO
			except:
				raise Exception(f"Error reading cache {self.filename}: {i} -> {self.locs[i]}")

		self.locked=False
		ret=EMStack2D(jnp.stack(stack))
		orts=Orientations(self.orts[nlist])
		tytx=np.array(self.tytx[nlist])
		return ret,orts,tytx

class EMStack():
	"""This class represents a stack of images in either an EMData, NumPy or Tensorflow representation, with easy interconversion
	- Shape of the array is {N,Z,Y,X}, as EMData, it is a list of N x EMData(X,Y,Z)
	- All images in the stack must have the same dimensions.
	- Only one representation exists at a time. Coercing to a new type is relatively expensive in terms of time.
	- Coerce routines will insure that the required representation exists, but others will be lost to insure self-consistency.
	- The convenience method numpy_list() will return a python list of N {Z,X,Y} NumPy arrays sharing memory with the EMData objects in the EMDataStack,
	  but beware, as coercing the EMDataStack to a differnt type will invalidate these NumPy arrays, and could create a crash!

	Individual images in the stack may be accessed using [n]
	"""

	def __init__(self,imgs=None):
		"""	imgs - one of:
		None
		filename, with optional ":" range specifier (see https://eman2.org/ImageFormats)
		single EMData object
		list or tuple of EMData objects
		numpy array, with first axis being image number {N,Z,Y,X} | {N,Y,X} | {N,X}
		Tensor, with first axis being image number {N,Z,Y,X} ...
		"""
		self._data=None	# representation in whatever the current format is
		self._npy_list=None # list of NumPy arrays sharing memory with EMData objects. RISKY - never copy this list, only use in place!
		self.set_data(imgs)
		self.apix=None

	def set_data(self,imgs):
		raise Exception("EMStack should not be used directly, please use EMStack3D, EMStack2D or EMStack1d")

	def __len__(self): return len(self._data)

	@property
	def shape(self):
		raise Exception("EMStack should not be used directly, please use EMStack3D, EMStack2D or EMStack1d")

	def __getitem__(self,key): return self._data[key]

	def __setitem__(self,key,value):
		raise Exception("Cannot set individual elements")

	@property
	def jax(self):
		self.coerce_jax()
		return self._data

	@jax.setter
	def jax(self,value):
		self.set_data(value)

	@property
	def numpy(self):
		self.coerce_numpy()
		return self._data

	@numpy.setter
	def numpy(self,value):
		self.set_data(value)

	@property
	def emdata(self):
		self.coerce_emdata()
		return self._data

	@emdata.setter
	def emdata(self,value):
		self.set_data(value)

	@property
	def numpy_list(self):
		self.coerce_emdata()
		if self._npy_list is not None: return
		self._npy_list=[i.numpy() for i in self._data]

	def coerce_emdata(self):
		"""Forces the current representation to EMData/NumPy"""
		if isinstance(self._data,list): return
		elif isinstance(self._data,np.ndarray): self._data=[from_numpy(i) for i in self._data]
		elif isinstance(self._data,jax.Array): self._data=from_jax(self._data,True)
		else: raise Exception(f"Invalid data in EMStack3D: {type(self._data)}")
		self._npy_list=None		# not necessary if already EMData list

	def coerce_numpy(self):
		if isinstance(self._data,np.ndarray): return
		elif self._npy_list is not None: self._data=np.stack(self._npy_list)
		elif isinstance(self._data,list): self._data=np.stack([i.numpy() for i in self._data])
		elif isinstance(self._data,jax.Array): self._data=np.array(self._data)
		else: raise Exception(f"Invalid data in EMStack3D: {type(self._data)}")
		self._npy_list=None		# not necessary if already EMData list

	def coerce_jax(self):
		if isinstance(self._data,jax.Array): return
		elif isinstance(self._data,list): self._data=to_jax(self._data)
		elif isinstance(self._data,np.ndarray): self._data=jnp.array(self._data)
		else: raise Exception(f"Invalid data in EMStack3D: {type(self._data)}")

	def center_clip(self,size):
		raise Exception("EMStack should not be used directly, please use EMStack3D, EMStack2D or EMStack1d")

	def do_fft(self,keep_type=False):
		raise Exception("EMStack should not be used directly, please use EMStack3D, EMStack2D or EMStack1d")

	def do_ift(self,keep_type=False):
		raise Exception("EMStack should not be used directly, please use EMStack3D, EMStack2D or EMStack1d")

	def calc_ccf(self,target,offset=0):
		"""Compute the cross correlation between each image in the stack and target, which may be a single image or another EMStack of the same size"""

		if isinstance(target,EMStack3D):
			return self.jax*jnp.conj(target)

	def mult(self,img):
		"""multiply each image in the stack by img"""
		if isinstance(img,jax.Array):
			self._data=self._data*img
		else: raise Exception("Only JAX data currently supported")

	def align_translate(ref,maxshift=-1):
		"""compute translational alignment of a stack of images to a same sized stack or single reference image.
		returns array of shifts the same size as the input stack. maxshift limits the maximum search area to +-maxshift
		on each axis"""
		pass

	def write_images(self,fsp=None,bits=12,n0=0):
		self.coerce_emdata()
		im_write_compressed(self._data,fsp,n0,bits)

	def downsample(self,newsize):
		"""Downsamples each image/volume in Fourier space such that its real-space dimensions after downsampling
		are "newsize" in all 2/3 dimensions. Downsampled images/volumes will be in Fourier space regardless of whether
		current stack is in real or Fourier space. This cannot be used to upsample (make images larger) and should
		not be used on rectangular images/volumes."""
		pass

class EMStack3D(EMStack):
	"""This class represents a stack of 3-D Volumes in either an EMData, NumPy or Tensorflow representation, with easy interconversion
	- Shape of the array is {N,Z,Y,X}, as EMData, it is a list of N x EMData(X,Y,Z)
	- All images in the stack must have the same dimensions.
	- Only one representation exists at a time. Coercing to a new type is relatively expensive in terms of time.
	- Coerce routines will insure that the required representation exists, but others will be lost to insure self-consistency.
	- The convenience method numpy_list() will return a python list of N {Z,X,Y} NumPy arrays sharing memory with the EMData objects in the EMDataStack,
	  but beware, as coercing the EMDataStack to a differnt type will invalidate these NumPy arrays, and could create a crash!

	Individual images in the stack may be accessed using emdata[n], tensor[n], numpy[n]
	"""

	def set_data(self,imgs):
		""" """
		if imgs is None:
			self._data=None
			self._npy_list=None
		elif isinstance(imgs,EMData):
			if imgs.get_ndim()!=3: raise Exception("EMStack3D only supports 3-D data")
			self._data=[imgs]
			self._npy_list=None
			self.apix=imgs[0]["apix_x"]
		elif isinstance(imgs,jax.Array) or isinstance(imgs,np.ndarray):
			if len(imgs.shape)==3:
				imgs=jnp.expand_dims(imgs,0)
			elif len(imgs.shape)!=4: raise Exception(f"EMStack3D only supports stacks of 3-D data, the provided images were {len(imgs.shape)}-D")
			self._data=imgs
			self._npy_list=None
		elif isinstance(imgs,str):
			self._data=EMData.read_images(imgs)
			self.apix=self._data[0]["apix_x"]
			if imgs[0].get_ndim()!=3: raise Exception(f"EMStack3D only supports stacks of 3-D data. {imgs} is {imgs[0].get_ndim()}-D")
			self._npy_list=None
		else:
			try:
				if not isinstance(imgs[0],EMData): raise Exception(f"EMDataStack cannot be initialized with a list of {type(imgs[0])}")
				self._data=list(imgs)		# copy the list, not the elements of the list
				self._npy_list=None
			except: raise Exception("EMDataStack may be initialized with None, a filename, an EMData object, a list/tuple of EMData objects, a NumPy array or a Tensor {N,Z,Y,X}")

	def __len__(self): return len(self._data)

	@property
	def shape(self):
		# note that the returned shape is N,Z,Y,X regardless of representation
		if isinstance(self._data,list): return(np.array((len(self._data),self._data[0]["nz"],self._data[0]["ny"],self._data[0]["nx"])))
		return(np.array(self._data.shape))

	def center_clip(self,size):
		size=int(size)
		if size<1: raise Exception("center_clip(size) must be called with a positive integer")
		shp=(self.shape-size)//2
		if isinstance(self._data,list):
			newlst=[im.get_clip(Region(int(shp[1]),int(shp[2]),int(shp[2]),size,size,size)) for im in self._data]
			return EMStack3D(newlst)
		elif isinstance(self._data,np.ndarray) or isinstance(self._data,jax.Array):
			newary=self._data[:,shp[1]:shp[1]+size,shp[2]:shp[2]+size,shp[3]:shp[3]+size]
			return EMStack3D(newary)

	def do_fft(self,keep_type=False):
		"""Computes the FFT of each image and returns a new EMStack3D. If keep_type is not set, will convert to Tensor before computing FFT."""
		if keep_type: raise Exception("do_fft: keep_type not functional yet")
		self.coerce_jax()

		return EMStack3D(jax_fft3d(self._data))

	def do_ift(self,keep_type=False):
		"""Computes the IFT of each image and returns a new EMStack3D. If keep_type is not set, will convert to Tensor before computing."""
		if keep_type: raise Exception("do_ift: keep_type not functional yet")
		self.coerce_jax()

		return jax_ift3d(self._data)

	def calc_ccf(self,target):
		"""Compute the cross correlation between each image in the stack and target, which may be a single image or another EMStack of the same size"""

		if isinstance(target,EMStack3D):
			return EMStack3D(self.jax*jnp.conj(target.jax))
		elif isinstance(target,jax.Array):
			return EMStack3D(self.jax*jnp.conj(target))
		else: raise Exception("calc_ccf: target must be either EMStack2D or single Tensor")

	def calc_flcf(self,target):
		""" Compute the fast local cross correlation between each image in the stack and a target, which may be a single image or another EMStack of the same size. Both image and target(s) need to be real, target need to have a boxsize smaller than the one of image"""
		if isinstance(target,EMStack3D):
			target = target.jax
		elif isinstance(target,jax.Array):
			pass
		else: raise Exception("calc_flcf: target must be either EMStack3D or single Tensor")
		tomof = jax_fft3d(self.jax)
		tomosqr_f = jax_fft3d(self.jax*self.jax)
		r2 = (target.shape[1]//2)**2
		cz,cy,cx = tomof.shape[1]//2, tomof.shape[2]//2, tomof.shape[3]
		mask = jax_fft3d(jnp.fromfunction(lambda z,y,x:jnp.where(((z-cz)**2+(y-cy)**2+(x-cx)**2)>r2,0,1), (cz*2,cy*2,(cx-1)*2),dtype=jnp.int32))
		norm = jax_ift3d(jax_phaseorigin3d((tomosqr_f)*jnp.conj(mask)))
		norm = norm/jnp.max(norm)
		pz,py,px = cz-target.shape[1]//2,cy-target.shape[2]//2,(cx-1)-target.shape[2]//2
		targetf = jax_fft3d(jnp.pad(target,((0,0),(pz,pz),(py,py),(px,px))))# pad to the same size with self
		ccf = jax_ift3d(jax_phaseorigin3d(tomof*jnp.conj(targetf)))
		return EMStack3D(ccf/(norm))


	def downsample(self,newsize):
		"""Downsamples each image/volume in Fourier space such that its real-space dimensions after downsampling
		are "newsize" in all 2/3 dimensions. Downsampled images/volumes will be in Fourier space regardless of whether
		current stack is in real or Fourier space. This cannot be used to upsample (make images larger) and should
		not be used on rectangular images/volumes."""

		return EMStack3D(jax_downsample_3d(self.jax,newsize))	# TODO: for now we're forcing this to be a tensor, probably better to leave it in the current format


class EMStack2D(EMStack):
	"""This class represents a stack of 2-D Images in either an EMData, NumPy or Tensorflow representation, with easy interconversion
	- Shape of the array is {N,Y,X}, as EMData, it is a list of N x EMData(X,Y)
	- All images in the stack must have the same dimensions.
	- Only one representation exists at a time. Coercing to a new type is relatively expensive in terms of time.
	- Coerce routines will insure that the required representation exists, but others will be lost to insure self-consistency.
	- The convenience method numpy_list() will return a python list of N {Z,X,Y} NumPy arrays sharing memory with the EMData objects in the EMDataStack,
	  but beware, as coercing the EMDataStack to a differnt type will invalidate these NumPy arrays, and could create a crash!

	Note that EMData headers are lost if converted to jax or numpy!

	Individual images in the stack may be accessed using emdata[n], jax[n], numpy[n]
	"""

	def set_data(self,imgs):
		""" """
		self._xforms=None
		self._df=None
		if imgs is None:
			self._data=None
			self._npy_list=None
		elif isinstance(imgs,EMData):
			if imgs.get_ndim()!=2: raise Exception("EMStack2D only supports 2-D data")
			self._data=[imgs]
			try: self._xforms=[imgs["xform.projection"]]
			except: pass
			try: self._df=np.array([imgs["ctf"].to_dict()["defocus"]])
			except:
				try:
					js=js_open_dict(info_name(imgs["ptcl_source_image"]))
					self._df=np.array([js["ctf"][0].to_dict()["defocus"]])
					js.close()
				except: pass
			self._npy_list=None
			self.apix=imgs[0]["apix_x"]
		elif isinstance(imgs,jax.Array) or isinstance(imgs,np.ndarray):
			if len(imgs.shape)!=3: raise Exception(f"EMStack2D only supports stacks of 2-D data, the provided images were {len(imgs.shape)}-D")
			self._data=imgs
			self._npy_list=None
		elif isinstance(imgs,str):
			self._data=EMData.read_images(imgs)
			self.apix=self._data[0]["apix_x"]
			try: self._xforms=[im["xform.projection"] for im in self._data]
			except: pass
			try: self._df=np.array([im["ctf"].to_dict()["defocus"] for im in self._data])
			except:
				try:
					df=[]
					for im in self._data:
						js=js_open_dict(info_name(im["ptcl_source_image"]))
						df.append(js["ctf"][0].to_dict()["defocus"])
						js.close()
					self._df=np.array(df)
				except: pass
			if self._data[0].get_ndim()!=2:
				if len(self._data)!=1 : raise Exception(f"EMStack2D only supports stacks of 2-D data or a single volume. {imgs} is a stack of {self._data[0].get_ndim()}-D")
				self._data=to_jax(self._data[0])
			self._npy_list=None
		else:
			try:
				if not isinstance(imgs[0],EMData): raise Exception(f"EMDataStack cannot be initialized with a list of {type(imgs[0])}")
				self._data=list(imgs)		# copy the list, not the elements of the list
				try: self._xforms=[im["xform.projection"] for im in self._data]
				except: pass
				try: self._df=np.array([im["ctf"].to_dict()["defocus"] for im in self._data])
				except:
					try:
						df=[]
						for im in self._data:
							js=js_open_dict(info_name(im["ptcl_source_image"]))
							df.append(js["ctf"][0].to_dict()["defocus"])
							js.close()
						self._df=np.array(df)
					except: pass
				self._npy_list=None
			except: raise Exception("EMStack2D may be initialized with None, a filename, an EMData object, a list/tuple of EMData objects, a NumPy array or a Tensor {N,Y,X}")

	def __len__(self): return len(self._data)

	@property
	def shape(self):
		# note that the returned shape is N,Y,X regardless of representation
		if isinstance(self._data,list): return(np.array((len(self._data),self._data[0]["ny"],self._data[0]["nx"])))
		return(np.array(self._data.shape))

	@property
	def orientations(self):
		"""returns an Orientations object and tytx array for the current images if available or None if not. If postxf is a Transform
		object, it will be applied to each Transform before generating the Orientations object"""
		if self._xforms is None: return None,None
		orts=Orientations()
		tytx=orts.init_from_transforms(self._xforms)
		if self._df is not None: tytx=jnp.stack([tytx[:,0],tytx[:,1], self._df], axis=-1)
		else: tytx = jnp.stack([tytx[:,0],tytx[:,1], jnp.zeros(tytx.shape[0])], axis=-1)
		return orts,tytx

	def orientations_withxf(self,prexf=None):
		"""Will apply prexf rotatation to each orientation before returning. Typically used to modify orientations for
		symmetry. eg = s=Symmetries.get("c4"), orientations_withxf(s.get_sym(1)) will return projections oriented to
		the first (non identity) symmetric orientation"""
		if self._xforms is None: return None,None
		if prexf is None: return self.orientations

		xfc=[xf*prexf for xf in self._xforms]
		orts=Orientations()
		tytx=orts.init_from_transforms(xfc)
		if self._df is not None: tytx=jnp.stack([tytx[:,0],tytx[:,1], self._df], axis=-1)
		else: tytx = jnp.stack([tytx[:,0],tytx[:,1], jnp.zeros(tytx.shape[0])], axis=-1)
		return orts,tytx


	def center_clip(self,size):
		try: size=np.array((int(size),int(size)))
		except: size=(int(size[0]),int(size[1]))
		if size[0]<1 or size[1]<1: raise Exception("center_clip(size) must be called with a positive integer")
		shp=(self.shape[1:]-size)//2
		if isinstance(self._data,list):
			newlst=[im.get_clip(Region(int(shp[0]),int(shp[1]),int(size[0]),int(size[1]))) for im in self._data]
			return EMStack2D(newlst)
		elif isinstance(self._data,np.ndarray) or isinstance(self._data,jax.Array):
			newary=self._data[:,shp[0]:shp[0]+size[0],shp[1]:shp[1]+size[1]]
			return EMStack2D(newary)

	def do_fft(self,keep_type=False):
		"""Computes the FFT of each image and returns a new EMStack3D. If keep_type is not set, will convert to Tensor before computing FFT."""
		if keep_type: raise Exception("do_fft: keep_type not functional yet")
		self.coerce_jax()

		return EMStack2D(jax_fft2d(self._data))

	def do_ift(self,keep_type=False):
		"""Computes the IFT of each image and returns a new EMStack3D. If keep_type is not set, will convert to Tensor before computing."""
		if keep_type: raise Exception("do_ift: keep_type not functional yet")
		self.coerce_jax()

		return EMStack2D(jax_ift2d(self._data))

	def calc_ccf(self,target,center=True,offset=0):
		"""Compute the cross correlation between each image in the stack and target, which may be a single image or another EMStack of the same size.
	If center is True, will shift the phase origin so zero shift corresponds to the middle of the image"""

		if center:
			if isinstance(target,EMStack2D) and offset!=0:
				return EMStack2D(jax_phaseorigin2d(self.jax[:-offset]*jnp.conj(target.jax[offset:])))
			elif isinstance(target,EMStack2D):
				return EMStack2D(jax_phaseorigin2d(self.jax*jnp.conj(target.jax)))
			elif isinstance(target,jax.Array) and offset==0:
				return EMStack2D(jax_phaseorigin2d(self.jax*jnp.conj(target)))
			else: raise Exception("calc_ccf: target must be either EMStack2D or single Tensor")
		else:
			if isinstance(target,EMStack2D) and offset!=0:
				return EMStack2D(self.jax[:-offset]*jnp.conj(target.jax[offset:]))
			elif isinstance(target,EMStack2D):
				return EMStack2D(self.jax*jnp.conj(target.jax))
			elif isinstance(target,jax.Array) and offset==0:
				return EMStack2D(self.jax*jnp.conj(target))
			else: raise Exception("calc_ccf: target must be either EMStack2D or single Tensor")

	def convolve(self,target):
		"""Compute the convolution between each image in the stack and target, which may be a single image or another EMStack of the same size"""

		if isinstance(target,EMStack2D):
			return EMStack2D(self.jax*target.jax)
		elif isinstance(target,jax.Array):
			return EMStack2D(self.jax*target)
		else: raise Exception("calc_ccf: target must be either EMStack2D or single Tensor")

	def downsample(self,newsize):
		"""Downsamples each image/volume in Fourier space such that its real-space dimensions after downsampling
		are "newsize" in all 2/3 dimensions. Downsampled images/volumes will be in Fourier space regardless of whether
		current stack is in real or Fourier space. This cannot be used to upsample (make images larger) and should
		not be used on rectangular images/volumes."""

		if newsize==self.shape[1]: return EMStack2D(self.jax) # this won't copy, but since the tensor is constant should be ok?
		return EMStack2D(jax_downsample_2d(self.jax,newsize))	# TODO: for now we're forcing this to be a tensor, probably better to leave it in the current format

	def normalize_standard(self):
		"""linear transform of values such that the mean of each image is 0 and the standard deviation is 1"""
		return EMStack2D((self.jax-jnp.mean(self.jax,axis=(1,2)))[:,jnp.newaxis,jnp.newaxis]/jnp.std(self.jax,axis=(1,2))[:,jnp.newaxis,jnp.newaxis])

	def normalize_edgemean(self):
		"""linear transform of values such that the mean of each image over a 3-pixel wide border around the edge
		of each image is set to zero and the standard deviation is 1"""
		N,ny,nx=self.shape
		edgemean=(jnp.mean(self.jax[:,0:ny-3,0:3],axis=(1,2))+
			jnp.mean(self.jax[:,ny-3:ny,0:nx-3],axis=(1,2))+
			jnp.mean(self.jax[:,3:ny,nx-3:nx],axis=(1,2))+
			jnp.mean(self.jax[:,0:3,3:nx],axis=(1,2)))/4.0

		return EMStack2D((self.jax-edgemean[:,jnp.newaxis,jnp.newaxis])/jnp.std(self.jax,axis=(1,2))[:,jnp.newaxis,jnp.newaxis])

	def center_align_seq(self,region_size=-1):
		"""Aligns a stack of (real space) images using the middle 1/2 (in x/y) of each image, and the middle image as a starting point.
		designed to do a rough alignment of tilt series. region_size is the size in pixels of the region around the center of each
		image to use for alignment. Default is ~1/2 the box size. All sizes adjusted to a good size"""

		#TODO - While functional, this whole method seems inefficient, particularly in terms of JAX
		if region_size<=32 : region_size=self.shape[1]//2
		region_size=good_size(region_size)

		cens=self.center_clip(region_size).do_fft()
		nc=cens.shape[0]//2
		atop0=EMStack2D(cens.jax[nc:-1])
		atop1=EMStack2D(cens.jax[nc+1:])
		alitop=atop0.align_translate(atop1)

		abot0=EMStack2D(cens.jax[jnp.arange(nc,0,-1)])
		abot1=EMStack2D(cens.jax[jnp.arange(nc-1,-1,-1)])
		alibot=abot0.align_translate(abot1)

#		print(alibot,alitop)

		self.coerce_emdata()
		dx,dy=0,0
		for i,sh in enumerate(alibot):
			dx+=sh[1]
			dy+=sh[0]
			self.emdata[nc-i-1].translate(-int(dx),-int(dy),0)

		dx,dy=0,0
		for i,sh in enumerate(alitop):
			dx+=sh[1]
			dy+=sh[0]
			self.emdata[nc+i+1].translate(-int(dx),-int(dy),0)


	def align_translate(self,ref,maxshift=-1):
		"""compute translational alignment of a stack of images to a same sized stack (or single) of reference images.
		returns array of (dy,dx) the same size as the input stack required to bring each "this" image into alignment with "ref". maxshift limits the maximum search area to +-maxshift
		on each axis. If maxshift is unspecified -> box size //4"""

		ny,nx=self.shape[1:]
		if self.jax.dtype==jnp.complex64 :
			nx=(nx-1)*2
			data=self
		else:
			data=self.do_fft()
			ref=ref.do_fft()

		if maxshift<=0: maxshift=ny//4

		ccfs=data.calc_ccf(ref)
		ccfsr=ccfs.do_ift()
		ccfsrc=ccfsr[:,ny//2-maxshift:ny//2+maxshift,nx//2-maxshift:nx//2+maxshift]		# only search a region around the center defined by maxshift

		# reshaped CCF so we can use reduce_max on it
		ccfsrs=jnp.reshape(ccfsrc,(ccfsrc.shape[0],maxshift*2*maxshift*2))

		# The y,x coordinates of the peak location
		peaks=maxshift-jnp.unravel_index(jnp.argmax(ccfsrs,1),(maxshift*2,maxshift*2))

		return jnp.transpose(peaks)

class EMAN3Ctf():
	# TODO: Change this to be accurate
	"""This class represents CTF conditions for an image."""

	def __init__(self, ctf=None, dfmajor = None, dfminor=None, dfang=None, voltage=None, cs=None, ampcont=None,apix=None, bfactor=None):
		"""	imgs - one of:
		None
		filename, with optional ":" range specifier (see https://eman2.org/ImageFormats)
		single EMData object
		list or tuple of EMData objects
		numpy array, with first axis being image number {N,Z,Y,X} | {N,Y,X} | {N,X}
		Tensor, with first axis being image number {N,Z,Y,X} ...
		"""

		if ctf != None: # TODO: make sure you can pass either EMAN2 or EMAN3 ctf
			self.voltage=ctf.voltage
			self.cs=ctf.cs
			self.ampcont=ctf.ampcont
			self.apix=ctf.apix
			self.bfactor=ctf.bfactor
			self.dfang=ctf.dfang
			self.dfmajor=ctf.defocus+ctf.dfdiff/2
			self.dfminor=ctf.defocus-ctf.dfdiff/2
		if voltage !=None: self.voltage=voltage
		if cs != None: self.cs=cs
		if ampcont!= None: self.ampcont=ampcont # TODO: Store internally as phase shift
		if apix!=None: self.apix=apix
		if dfmajor!=None: self.dfmajor=dfmajor
		if dfminor!=None: self.dfminor=dfminor
		if dfang!=None: self.dfang=dfang

	@property
	def wavelength(self):
		"""Calculates relativistic wavelength of electron (in A) for a given voltage"""
		return 12.2639/np.sqrt(self.voltage*1000.0+0.97845*self.voltage*self.voltage)

	@property
	def defocus_step(self):
		"""Calculates the step in defocus in um that results in a phase shift of no more than 90deg through Nyquist.
		voltage- the voltage of the microscope, in keV
		apix-The apix of the original image, to know Nyquist"""
		return 2*self.apix*self.apix/(self.wavelength*10000)

	def get_defocus(self):
		"""Returns the EMAN style defocus"""
		return (self.dfmajor+self.dfminor)/2

	def df(self, ang):
		"""Returns defocus as a function of angle, ang should be in radians"""
		if self.dfmajor == self.dfminor: return self.dfmajor
		else: return self.get_defocus + ((self.dfmajor - self.dfminor)/2)*cos(2*ang-(2*pi/180)*self.dfang)

	def df_range(self, ang, defocus, dfang):
		"""Returns defocus as a function of angle using a provided defocus instead of the one in self
		Inputs:
		ang--The angle to calculate the defocus at. It should be in radians
		defocus--The defocus to use instead of calculating average defocus from self
		dfang--The astigmatism angle to use instead of dfang from self (should be in degrees)"""
		if self.dfmajor == self.dfminor: return defocus
		else: return defocus + ((self.dfmajor - self.dfminor)/2)*cos(2*ang-(2*pi/180)*dfang)

	def df_img(self, ny, defocus, dfang):
		return jnp.array(jnp.vstack((jnp.fromfunction(lambda y,x: self.df_range(jnp.arctan2(x,y), defocus, dfang),(ny//2,ny//2+1)),jnp.fromfunction(lambda y,x: self.df_range(jnp.arctan2(x,(ny//2-y)),defocus,dfang),(ny//2,ny//2+1)))))

	df_defocus_img = jax.vmap(df_img, in_axes=(None, None, 0, None))
	df_dfang_img = jax.vmap(df_img, in_axes=(None, None, None, 0))

	def get_phase(self):
		"""Returns the phase shift based on the percent amplitude contrast"""
		if self.ampcont>-100.0 and self.ampcont<=100.0: return asin(self.ampcont/100)
		elif ampcont>100: return pi-asin(2-self.ampcont/100)
		else: return -pi-asin(-2-self.ampcont/100)

	def set_phase(self, phase):
		"""Sets the amplitude contrast depending on the given phase shift. Phase should be given in radians."""
		if phase>=-pi/2.0 and phase < pi/2.0: self.ampcont=sin(phase)*100
		elif phase >=pi/2.0 and phase<=pi: self.ampcont=(2-sin(phase))*100
		elif phase>-pi and phase<-pi/2: self.ampcont=(-2-sin(phase))*100
		else: print(f"Phase={phase*180/pi} degrees out of range")

	def to_emanCTF(self):
		"""Returns an EMAN CTF object with the information in the object"""
		ctf = EMAN2Ctf()
		ctf.from_dict({"voltage": self.voltage, "cs": self.cs, "ampcont": self.ampcont, "apix": self.apix, "bfactor": self.bfactor, "dfang": self.dfang, "defocus": self.get_defocus, "dfdiff": self.dfmajor-self.dfminor})
		return ctf

	def compute_2d_stack_complex(self, ny, ctf_type, var_range, var_name, apix=None):
	def to_dict(self):
		"""Returns a dictionary representation of the object"""
		return {"voltage": self.voltage, "cs": self.cs, "ampcont": self.ampcont, "apix":self.apix, "defocus":self.defocus, "dfang":self.dfang, "dfdiff":self.dfdiff}

	def to_jsondict(self): # Included here instead of EMAN3jsondb.py because of recursive imports
		"""Creates a dict representation of the object for JSON serialization purposes."""
		ret = self.to_dict()
		ret["__class__"] = "EMAN3Ctf"
		return ret

		"""Returns a stack of CTF images, with size ny, of type ctf_type. The stack will vary var_name in the range var_range.
			Inputs:
			ny--Size of images to make. Will assume it corresponds to the apix unless overriden
			ctf_type--One of 'amplitude', 'sign'
			var_range--The range to let the variable vary
			var_name--Which variable to vary, right now 'defocus' and 'dfang' are supported
			apix--Override in case ny is downsampled from the original
			Returns:
			ctfary--A jax array of the CTF images
			step--The step that was used to fill the range"""
		if apix==None: apix=self.apix
		wl = self.wavelength
		g1=(jnp.pi/2.0)*self.cs*1.0e7*pow(wl,3)
		g2=jnp.pi*wl*10000.0
		phase=jnp.pi/2.0-self.get_phase()
		rad2 = rad2_img(ny)/(apix*apix*ny*ny)
		if var_name == "defocus":
			step = self.defocus_step
			dfimg = self.df_defocus_img(ny, jnp.arange(var_range[0], var_range[1], step, jnp.complex64), self.dfang)
		elif var_name == "dfang":
			step = 1 # TODO: If I can, calculate a better angular step eg based on defocus_step
			dfimg = self.df_dfang_img(ny, self.get_defocus(), jnp.arange(var_range[0], var_range[1], step, jnp.complex64))
		else:
			print("var_name was not recognized. It should be one of 'defocus' or 'dfang'")
			return
		if ctf_type == "amplitude":
			ctf = jnp.cos(-g1*rad2*rad2+g2*rad2*dfimg-phase)
		elif ctf_type == "sign":
			ctf = jnp.sign(jnp.cos(-g1*rad2*rad2+g2*rad2*dfimg-phase))
		else:
			print("The ctf_type was not recognized, it should be one of 'amplitude', 'sign'.")
			return
		return ctf, step


class Orientations():
	"""This represents a set of orientations, with a standard representation of an XYZ vector where the vector length indicates the amount
		of rotation with a length of 0.5 corresponding to 180 degrees. This form is a good representation for deep learning minimization strategies
		which conventionally span a range of ~1.0. This form can be readily interconverted to EMAN2 Transform objects or transformation matrices
		for use with Gaussians. Note that this object handles only orientation, not translation.
		NOTE: unlike some other objects, the Nx3 array is XYZ order, NOT ZYX order! This matches the vector ordering of Gaussians
	"""

	def __init__(self,xyzs=None):
		"""Initialize with the number of orientations, a N x 3 matrix or a symmetry"""
		if isinstance(xyzs,int):
			if xyzs<=0: self._data=None
			else: self._data=np.zeros((xyzs,3),np.float32)
		elif isinstance(xyzs,str): self.init_symmetry(xyzs)
		else:
			try: self._data=np.array(xyzs,np.float32)
			except: raise Exception("Orientations must be initialized with an integer (number of orientations) or a N x 3 numpy array")


	def __len__(self): return len(self._data)

	def __getitem__(self,key):
		"""Return the keyed Gaussian parameter, may return a tensor or numpy array. G[i] returns the 3-vector for the i'th Orientation"""
		return self._data[key]

	def __setitem__(self,key,value):
		# if the Gaussians are a tensor, we turn it back into numpy for modification
		self.coerce_numpy()
		self._data[key]=value

	def __len__(self): return self._data.shape[0]

	def coerce_jax(self):
		if not isinstance(self._data,jax.Array): self._data=jnp.array(self._data,jnp.float32)

	def coerce_numpy(self):
		if isinstance(self._data,jax.Array): self._data=np.array(self._data)

	@property
	def jax(self):
		self.coerce_jax()
		return self._data

	@property
	def numpy(self):
		self.coerce_numpy()
		return self._data

	def init_from_transforms(self,xformlist):
		"""Replaces current contents of Orientations object with orientations from a list of Transform objects,
		returns tytx array with any translations (not stored within Orientations)"""
		self._data=np.zeros((len(xformlist),3))
		tytx=[]
		for i,x in enumerate(xformlist):
			r=x.get_rotation("spinvec")
			self._data[i]=(r["v1"],r["v2"],r["v3"])
			tytx.append((x.get_trans_2d()[1],x.get_trans_2d()[0]))

		return(np.array(tytx))

	def init_symmetry(self,sym="c1"):
		"""Replaces current orientations with those necessary to symmetrize a volume or set of points. No translation,
	so will not work for helical symmetries"""
		s=Symmetries.get(sym)
		self._data=np.zeros((s.get_nsym(),3))
		for i in range(s.get_nsym()):
			r=s.get_sym(i).get_rotation("spinvec")
			self._data[i]=(r["v1"],r["v2"],r["v3"])

		return

	def transforms(self,tytx=None):
		"""converts the current orientations to a list of Transform objects"""
		if tytx is not None:
			return [Transform({"type":"spinvec","v1":float(self._data[i][0]),"v2":float(self._data[i][1]),"v3":float(self._data[i][2]),"tx":float(tytx[i][1]),"ty":float(tytx[i][0])}) for i in range(len(self._data))]

		return [Transform({"type":"spinvec","v1":self._data[i][0],"v2":self._data[i][1],"v3":self._data[i][2]}) for i in range(len(self._data))]

	def to_mx2d(self,swapxy=False):
		"""Returns the current set of orientations as a 2 x 3 x N matrix which will transform a set of 3-vectors to a set of
		2-vectors, ignoring the resulting Z component. Typically used with Gaussians to generate projections.

		To apply to a set of vectors:
		mx=self.to_mx2d()
		vecs=tf.constant(((1,0,0),(0,1,0),(0,0,1),(2,2,2),(1,1,0),(0,1,1)),dtype=tf.float32)

		tf.transpose(tf.matmul(mx[:,:,0],tf.transpose(vecs)))
		or
		tf.einsum("ij,kj->ki",mx[:,:,0],vecs)

		if swapxy is set, then the input vector is XYZ, but output is YX. This corrects for the fact that Gaussians are XYZA,
		but images are YX.
		"""

		self.coerce_jax()

		# Adding a tiny value avoids the issue with zero rotations. While it would be more correct to use a conditional
		# it is much slower, and the tiny pertutbation should not significantly impact the math.
		l=jnp.linalg.norm(self._data,axis=1)+1.0e-37

		w=jnp.cos(pi*l)  # cos "real" component of quaternion
		s=jnp.sin(-pi*l)/l
		q=jnp.transpose(self._data)*s		# transpose makes the vectorized math below work properly
		if swapxy :
			mx=jnp.array(((2*q[0]*q[1]+2*q[2]*w,1-(2*q[0]*q[0]+2*q[2]*q[2]),2*q[1]*q[2]-2*q[0]*w),
			(1-2*(q[1]*q[1]+q[2]*q[2]),2*q[0]*q[1]-2*q[2]*w,2*q[0]*q[2]+2*q[1]*w)))
		else:
			mx=jnp.array(((1-2*(q[1]*q[1]+q[2]*q[2]),2*q[0]*q[1]-2*q[2]*w,2*q[0]*q[2]+2*q[1]*w),
			(2*q[0]*q[1]+2*q[2]*w,1-(2*q[0]*q[0]+2*q[2]*q[2]),2*q[1]*q[2]-2*q[0]*w)))
		return jnp.array(mx)

	def to_mx3d(self):
		"""Returns the current set of orientations as a 3 x 3 x N matrix which will transform a set of 3-vectors to a set of
		rotated 3-vectors. Can be used to symmetrize a set of 3-vectors

		To apply to a set of vectors:
		ort=Orientations("c4")
		mx=self.to_mx3d()
		vecs=tf.constant(((1,0,0),(0,1,0),(0,0,1),(2,2,2)),dtype=tf.float32)
		vecs_c4_0=jnp.matmul(vecs,mx[:,:,0])    # this is all vectors rotated by the first C4 transformation
		vecs_c4_1=jnp.matmul(vecs,mx[:,:,1])	# all vectors rotated by the second C4 transformation
		"""

		self.coerce_jax()

		# Adding a tiny value avoids the issue with zero rotations. While it would be more correct to use a conditional
		# it is much slower, and the tiny pertutbation should not significantly impact the math.
		l=jnp.linalg.norm(self._data,axis=1)+1.0e-37

		w=jnp.cos(pi*l)  # cos "real" component of quaternion
		s=jnp.sin(-pi*l)/l
		q=jnp.transpose(self._data)*s		# transpose makes the vectorized math below work properly

		mx=jnp.array(((1-2*(q[1]*q[1]+q[2]*q[2]),2*q[0]*q[1]-2*q[2]*w,2*q[0]*q[2]+2*q[1]*w),
		(2*q[0]*q[1]+2*q[2]*w,1-(2*q[0]*q[0]+2*q[2]*q[2]),2*q[1]*q[2]-2*q[0]*w),
		(2*q[0]*q[2]-2*q[1]*w,2*q[1]*q[2]+2*q[0]*w,1-(2*q[0]*q[0]+2*q[1]*q[1]))))
		return jnp.array(mx)

	def to_mx3da(self):
		"""Returns the current set of orientations as a 4 x 4 x N matrix which will transform a set of 3-vectors + amplitude to a set of
		rotated 3-vectors with unchanged amplitudes. Typically used with Gaussians to symmetrize or rotate corresponding to a volume.

		To apply to a set of vectors:
		ort=Orientations("c4")
		mx=self.to_mx3d()
		vecs=tf.constant(((1,0,0),(0,1,0),(0,0,1),(2,2,2)),dtype=tf.float32)
		vecs_c4_0=jnp.matmul(vecs,mx[:,:,0])    # this is all vectors rotated by the first C4 transformation
		vecs_c4_1=jnp.matmul(vecs,mx[:,:,1])	# all vectors rotated by the second C4 transformation
		"""

		self.coerce_jax()

		# Adding a tiny value avoids the issue with zero rotations. While it would be more correct to use a conditional
		# it is much slower, and the tiny pertutbation should not significantly impact the math.
		l=jnp.linalg.norm(self._data,axis=1)+1.0e-37

		w=jnp.cos(pi*l)  # cos "real" component of quaternion
		s=jnp.sin(-pi*l)/l
		q=jnp.transpose(self._data)*s		# transpose makes the vectorized math below work properly
		c0=jnp.zeros(w.shape)
		c1=c0+1.0

		mx=jnp.array(((1-2*(q[1]*q[1]+q[2]*q[2]), 2*q[0]*q[1]-2*q[2]*w, 2*q[0]*q[2]+2*q[1]*w, c0),
		(2*q[0]*q[1]+2*q[2]*w, 1-(2*q[0]*q[0]+2*q[2]*q[2]), 2*q[1]*q[2]-2*q[0]*w, c0),
		(2*q[0]*q[2]-2*q[1]*w, 2*q[1]*q[2]+2*q[0]*w, 1-(2*q[0]*q[0]+2*q[1]*q[1]), c0),
		(c0,c0,c0,c1)))
		return jnp.array(mx)


class Gaussians():
	"""This represents a set of Gaussians with x,y,z,amp parameters (but no width). Representation is a N x 4 numpy array or tensor (x,y,z,amp) ],
x,y,z are ~-0.5 to ~0.5 (typ) and amp is 0 to ~1. A scaling factor (value -> pixels) is applied when generating projections. """

	def __init__(self,gaus=0):
		if isinstance(gaus,int):
			if gaus<=0: self._data=None
			else: self._data=np.zeros((gaus,4),np.float32)
		elif isinstance(gaus,str):
			self._data=np.loadtxt(gaus,delimiter="\t")
			if self._data.shape[1] != 4:
				raise Exception(f"File {gaus} has shape {self._data.shape()}. Should be Nx4")
		else:
			try: self._data=np.array(gaus,np.float32)
			except: raise Exception("Gaussians must be initialized with an integer (number of Gaussians) or N x 4 matrix")

	def __getitem__(self,key):
		"""Return the keyed Gaussian parameter, may return a tensor or numpy array. G[i] returns the 4-vector for the i'th Gaussian"""
		return self._data[key]

	def __setitem__(self,key,value):
		# if the Gaussians are a tensor, we turn it back into numpy for modification
		self.coerce_numpy()
		self._data[key]=value

	def __len__(self): return len(self._data)

	def coerce_jax(self):
		if not isinstance(self._data,jax.Array): self._data=jnp.array(self._data,jnp.float32)

	def coerce_numpy(self):
		if isinstance(self._data,jax.Array): self._data=np.array(self._data)

	def add_array(self,array):
		self._data+=array

	@property
	def jax(self):
		self.coerce_jax()
		return self._data

	@property
	def numpy(self):
		self.coerce_numpy()
		return self._data

	def init_from_map(self,vol,res,minratio=0.1,apix=None):
		"""Replace the current set of Gaussians with a set of Gaussians generated from a 3-D map by progressive Gaussian decomposition.
		The map is filtered to res, then the highest amplitude peak is assigned to the first Gaussian. After subtracting that Gaussian from the
		map the process is repeated until the ratio of the next amplitude to the first amplitude falls below the minratio limit.
		vol - a single EMData, numpy or tensorflow volume
		res - lowpass filter resolution and the FWHM of the Gaussians to be subtracted, specified in A, use apix for non-EMData volumes
		minratio - minimum peak ratio, >0
		apix - A/pix override"""

		if isinstance(vol,jax.Array): emd=from_jax(vol)
		elif isinstance(vol,np.ndarray): emd=from_numpy(vol)
		elif isinstance(vol,EMData): emd=vol
		else: raise Exception("init_from_map: vol must be EMData, Tensor or NumPy Array")

		if apix is not None: emd["apix_x"],emd["apix_y"],emd["apix_z"]=apix,apix,apix

		# The actual Gaussian segmentation
		seg=emd.process("segment.gauss",{"minratio":minratio,"width":res,"skipseg":1})
		amps=np.array(seg["segment_amps"])
		centers=np.array(seg["segment_centers"]).reshape(len(amps),3)
		centers/=(emd["nx"],emd["ny"],emd["nz"])
		centers-=(0.5,0.5,0.5)
		amps/=max(amps)
		self._data=np.concatenate((centers.transpose(),amps.reshape((1,len(amps))))).transpose()

	def replicate(self,n=2,dev=0.01):
		"""Makes n copies of the current Gaussians shifted by a small random amount to improve the level of detail without
significantly altering the spatial distribution. Note that amplitudes are also perturbed by the same distribution. Default
stddev=0.01"""
		if n<=1 : return
		rng=np.random.default_rng()
		self.coerce_jax()
		dups=[self._data+rng.normal(0,dev,self._data.shape) for i in range(n)]
		self._data=jnp.concat(dups,axis=0)

	def replicate_abs(self,ngaus,dev=0.01):
		"""Makes copies of the current Gaussians shifted by a small random amount to improve the level of detail without
significantly altering the spatial distribution. ngaus specifies the total number of desired gaussians after replication. Note that amplitudes are also perturbed by the same distribution. Default stddev=0.01"""
		if len(self)==ngaus : return
		if len(self)>ngaus :
			print("Warning, replicate targeted fewer gaussians than currently exist! Unchanged")
			return
		rng=np.random.default_rng()
		self.coerce_jax()
		dups=[self._data]
		for i in range(ngaus//len(self)-1): dups.append(self._data+rng.normal(0,dev,self._data.shape))
		tail=ngaus%len(self)
		if tail>0: dups.append((self._data+rng.normal(0,dev,self._data.shape))[:tail])
		self._data=jnp.concat(dups,axis=0)

	def replicate_sym(self,sym="c1"):
		"""Makes symmertry related copies of each Gaussian so the full Gaussians object matches the specified symmetry."""
		symorts=Orientations(sym)
		mx=np.moveaxis(symorts.to_mx3da(),2,0)	# gets matrix as Nsym x 4 x 4

		self._data=jnp.matmul(self._data,mx).reshape(mx.shape[0]*self._data.shape[0],4)	# actual matrix multiplication becomes (Nsym,Ngau,4)

	def norm_filter(self,sig=0.5,rad_downweight=-1,cyl_mask=-1):
		"""Rescale the amplitudes so the maximum is 1, with amplitude below mean+sig*sigma removed.
		rad_downweight, if >0 will apply a radial linear amplitude decay beyond the specified radius to the corner of the cube. eg - 0.5 will downweight the corners. Downweighting only works if Gaussian coordinate range follows the -0.5 - 0.5 standard range for the box.
		cyl_mask, if >0 will apply a cylindrical mask along xz to account for what space in 3d the tilts cover. The value is assumed to be the Gaussian corrdinate equivalent of 1px in the final boxsize """
		self.coerce_jax()
		self._data=self._data*jnp.array((1.0,1.0,1.0,1.0/jnp.max(jnp.absolute(self._data[:,3]))))		# "normalize" amplitudes so max amplitude is scaled to 1.0, not sure how necessary this really is
		if cyl_mask>0:
			self._data = self._data[jnp.linalg.norm(self._data[:,(0,2)], axis=1)<0.5-cyl_mask*10]
			self._data = self._data[(self._data[:,1]>=-0.5+cyl_mask*3) & (self._data[:,1]<=0.5-cyl_mask*3)]
		if rad_downweight>0:
			famp=jnp.absolute(self._data[:,3])*(1.0-jax.nn.relu(jnp.linalg.vector_norm(self._data[:,:3],axis=1)-rad_downweight))
		else: famp=jnp.absolute(self._data[:,3])
		thr=jnp.mean(famp)+sig*jnp.std(famp)
		self._data=self._data[famp>thr]			# remove any gaussians with amplitude below threshold

	def mask(self,mask):
		"""
		"""
		if isinstance(mask,EMStack2D) : mask=mask.jax
		elif isinstance(mask,EMData) : mask=to_jax(mask)
		elif isinstance(mask,np.ndarray) : mask=jnp.array(mask)
		elif not isinstance(mask,jax.Array) : raise Exception("invalid data type for Gaussians.mask")

		coords=((self.jax+0.5)*mask.shape[1]).astype(jnp.int32)
		cmask=(mask[coords[:,2],coords[:,1],coords[:,0]]).astype(jnp.bool_)
		return Gaussians(self.jax[cmask])


	def project_simple(self,orts,boxsize,tytx=None):
		"""Generates a tensor containing a simple 2-D projection (interpolated delta functions) of the set of Gaussians for each of N Orientations in orts.
		orts - must be an Orientations object
		tytx =  is an (optional) N x 2+ vector containing an in-plane translation in unit (-0.5 - 0.5) coordinates to be applied to the set of Gaussians for each Orientation.
		boxsize in pixels. Scaling factor is equal to boxsize, such that -0.5 to 0.5 range covers the box.

		With these definitions, Gaussian coordinates are sampling-independent as long as no box size alterations are performed. That is, raw projection data
		used for comparisons should be resampled without any "clip" operations.
		"""
		self.coerce_jax()

#		proj=tf.zeros((len(orts),boxsize,boxsize))		# projections
		proj2=[]
		mx=orts.to_mx2d(swapxy=True)
		if tytx is None: tytx=jnp.zeros
		return EMStack2D(gauss_project_simple_fn(self._data,mx,boxsize,tytx))

	def project_ctf(self,orts,ctf_stack,boxsize,dfrange,dfstep,tytx=None):
		"""Generates a tensor containing a phase-flipped 2-D projection of the set of Gaussians for each of N orientations in orts. Phase flipping is off a single
		defocus value.

		orts-must be an Orientations object
		ctf_stack-A stack of ctf correction images at different defocuses
		tytx-is a Nx3+ vector containing an in-plane translation in units (-0.5,0.5) coordinates to be applied to the set of Gaussians in the first two columns and
			per particle defocus in the third
		boxsize-Value in pixels. Scaling factor is equal to boxsize, such that -0.5 to 0.5 range covers the box
		dfrange-a tuple of (min defocus,max defocus) in the project
		dfstep-The step between two correction images in ctf_stack
		With these definitions, Gaussian coordinates are sampling-independent as long as no box size alterations are performed. That is, raw projection data is
		used for comparisons should be resampled without any "clip" operations.
		"""
		self.coerce_jax()
		mx=orts.to_mx2d(swapxy=True) # 2d is fine becasue don't need z position when we apply one ctf regardless
		if tytx is None: tytx=jnp.zeros(3)
		return EMStack2D(gauss_project_ctf_fn(self._data,mx,ctf_stack._data,boxsize,dfrange[0],dfrange[1],dfstep,tytx))

	def project_layered_ctf(self,orts,ctf_stack,boxsize,apix,dfrange,dfstep,tytx=None):
		"""Generates a tensor containing a phase-flipped 2-D projection accounting for defocus levels of the set of Gaussians for each of N Orientations in orts
		orts-must be an Orientations object
		ctf_stack-A stack of ctf correction images at different defocuses
		tytx-is a Nx3+ vector containing an in-plane translation in units (-0.5,0.5) coordinates to be applied to the set of Gaussians in the first two columns and
			per particle defocus in the third
		boxsize-Value in pixels. Scaling factor is equal to boxsize, such that -0.5 to 0.5 range covers the box
		dfrange-a tuple of (min defocus,max defocus) in the project
		dfstep-The step between two correction images in ctf_stack

		With these definitions, Gaussian coordinates are sampling-independent as long as no box size alterations are performed. That is, raw projection data
		is used for comparisons should be resampled without any "clip" operations.
		"""
		self.coerce_jax()
		mx=orts.to_mx3d()
		if tytx is None: tytx=jnp.zeros(3)
		return EMStack2D(gauss_project_layered_ctf_fn(self._data,mx,ctf_stack._data,boxsize,dfrange[0],dfrange[1],dfstep,tytx))

	def volume(self,boxsize,zaspect=0.5):
		self.coerce_jax()

		if zaspect<=0 : zsize=boxsize
		else: zsize=good_size(boxsize*zaspect*2.0)

		return EMStack3D(gauss_volume_fn(self._data,boxsize,zsize))

		vol=jnp.zeros((zsize,boxsize,boxsize),dtype=jnp.float32)		# output

	def volume_np(self, boxsize, zaspect=0.5):
		"""A numpy implementation of volume since the jax at method makes a copy instead of updating in place when not in jit compiled and was causing OOM errors"""
		zsize=good_size(boxsize*zaspect*2.0)
		vol=np.zeros((zsize,boxsize,boxsize))
		xfgauss=np.flip((self[:,:3]+jnp.array((0.5,0.5,zaspect)))*boxsize,-1)		# shift and scale both x and y the same, reverse handles the XYZ -> ZYX EMData->Tensorflow issue

		xfgaussf=np.floor(xfgauss)
		xfgaussi=np.ndarray.astype(xfgaussf,np.int32)   # integer index
		xfgaussf=xfgauss-xfgaussf			       # remainder used for bilinear interpolation

		# messy trilinear interpolation
		bamp000=self._data[:,3]*(1.0-xfgaussf[:,0])*(1.0-xfgaussf[:,1])*(1.0-xfgaussf[:,2])
		bamp001=self._data[:,3]*(1.0-xfgaussf[:,0])*(1.0-xfgaussf[:,1])*(    xfgaussf[:,2])
		bamp010=self._data[:,3]*(1.0-xfgaussf[:,0])*(    xfgaussf[:,1])*(1.0-xfgaussf[:,2])
		bamp011=self._data[:,3]*(1.0-xfgaussf[:,0])*(    xfgaussf[:,1])*(    xfgaussf[:,2])
		bamp100=self._data[:,3]*(    xfgaussf[:,0])*(1.0-xfgaussf[:,1])*(1.0-xfgaussf[:,2])
		bamp101=self._data[:,3]*(    xfgaussf[:,0])*(1.0-xfgaussf[:,1])*(    xfgaussf[:,2])
		bamp110=self._data[:,3]*(    xfgaussf[:,0])*(    xfgaussf[:,1])*(1.0-xfgaussf[:,2])
		bamp111=self._data[:,3]*(    xfgaussf[:,0])*(    xfgaussf[:,1])*(    xfgaussf[:,2])
		bampall=np.concatenate([bamp000,bamp001,bamp010,bamp011,bamp100,bamp101,bamp110,bamp111],0)
		bposall=np.concatenate([xfgaussi,xfgaussi+(0,0,1),xfgaussi+(0,1,0),xfgaussi+(0,1,1),xfgaussi+(1,0,0),xfgaussi+(1,0,1),xfgaussi+(1,1,0),xfgaussi+(1,1,1)],0)

		# Jax doesn't give OOB errors it just puts them at the closest available index, but numpy will so need an extra check
		mask = np.logical_and(np.logical_and(np.logical_and(bposall[:,0]>=0, bposall[:,0]<zsize), np.logical_and(bposall[:,1]>=0, bposall[:,1]<boxsize)),np.logical_and(bposall[:,2]>=0, bposall[:,2]<boxsize))
		bampall=bampall[mask]
		bposall=bposall[mask]
		vol[bposall[:,0], bposall[:,1],bposall[:,2]]+=bampall

		return EMStack3D(vol)


	def volume_tiled(self,xsize,ysize,boxz,xtiles,ytiles,zaspect=0.5):
		"""A numpy version of creating a volume because the GPU can't allocate the full tomogram volume in memory. Returns an EMData object not a EMStack3D object."""
		self.coerce_numpy()

		zsize=good_size(boxz*zaspect*2.0)
		vol=np.zeros((zsize,ysize,xsize)) # TODO: Add a try/except state that if this fails to allocate it makes the binned tomogram? That way the refinement isn't wasted. Or add check beforehand

		xfgauss=np.flip((self._data[:,:3]+(0.25*xtiles,0.25*ytiles,zaspect))*(xsize/xtiles*2,ysize/ytiles*2,zsize),[-1])	# shift and scale both x and y the same, flip handles the XYZ -> ZYX EMData->Jax/Numpy issue

		xfgaussf=np.floor(xfgauss)
		xfgaussi=np.ndarray.astype(xfgaussf,np.int32)   # integer index
		xfgaussf=xfgauss-xfgaussf			       # remainder used for bilinear interpolation

		# messy trilinear interpolation
		bamp000=self._data[:,3]*(1.0-xfgaussf[:,0])*(1.0-xfgaussf[:,1])*(1.0-xfgaussf[:,2])
		bamp001=self._data[:,3]*(1.0-xfgaussf[:,0])*(1.0-xfgaussf[:,1])*(    xfgaussf[:,2])
		bamp010=self._data[:,3]*(1.0-xfgaussf[:,0])*(    xfgaussf[:,1])*(1.0-xfgaussf[:,2])
		bamp011=self._data[:,3]*(1.0-xfgaussf[:,0])*(    xfgaussf[:,1])*(    xfgaussf[:,2])
		bamp100=self._data[:,3]*(    xfgaussf[:,0])*(1.0-xfgaussf[:,1])*(1.0-xfgaussf[:,2])
		bamp101=self._data[:,3]*(    xfgaussf[:,0])*(1.0-xfgaussf[:,1])*(    xfgaussf[:,2])
		bamp110=self._data[:,3]*(    xfgaussf[:,0])*(    xfgaussf[:,1])*(1.0-xfgaussf[:,2])
		bamp111=self._data[:,3]*(    xfgaussf[:,0])*(    xfgaussf[:,1])*(    xfgaussf[:,2])
		bampall=np.concatenate([bamp000,bamp001,bamp010,bamp011,bamp100,bamp101,bamp110,bamp111],0)
		bposall=np.concatenate([xfgaussi,xfgaussi+(0,0,1),xfgaussi+(0,1,0),xfgaussi+(0,1,1),xfgaussi+(1,0,0),xfgaussi+(1,0,1),xfgaussi+(1,1,0),xfgaussi+(1,1,1)],0)

		# Jax doesn't give OOB errors it just puts them at the closest available index, but numpy will so need an extra check
		mask = np.logical_and(np.logical_and(np.logical_and(bposall[:,0]>=0, bposall[:,0]<zsize), np.logical_and(bposall[:,1]>=0, bposall[:,1]<ysize)),np.logical_and(bposall[:,2]>=0, bposall[:,2]<xsize))
		bampall=bampall[mask]
		bposall=bposall[mask]
		vol[bposall[:,0], bposall[:,1],bposall[:,2]]+=bampall

		return from_numpy(vol) # I think this makes a copy which isn't ideal but I don't know how to get it to an EMData object another way

# def gauss_project_simple_fn(gausary,mx,boxsize,tytx):
#	"""This exists as a function separate from the Gaussian class to better support JAX optimization. It is called by the corresponding Gaussian method.
#
#	Generates an array containing a simple 2-D projection (interpolated delta functions) of the set of Gaussians for each of N Orientations in orts.
#	gausary - a Gaussians.jax array
#	mx - an Orientations object converted to a stack of 2d matrices
#	tytx =  a N x 2+ vector containing an in-plane translation in unit (-0.5 - 0.5) coordinates to be applied to the set of Gaussians for each Orientation.
#	boxsize in pixels. Scaling factor is equal to boxsize, such that -0.5 to 0.5 range covers the box.
#
#	With these definitions, Gaussian coordinates are sampling-independent as long as no box size alterations are performed. That is, raw projection data
#	used for comparisons should be resampled without any "clip" operations.
#	"""
#
#	proj2=[]
#
#	# iterate over projections
#	# TODO - at some point this outer loop should be converted to a tensor axis for (potentially) better performance
#	# note that the mx dimensions have N as the 3rd not 1st component!
#	gpsf=jax.jit(gauss_project_single_fn,static_argnames=["boxsize"])
#
#	for j in range(mx.shape[2]):
#		proj2.append(gpsf(gausary,mx[:,:,j],boxsize,tytx[j]))
#
#	return jnp.stack(proj2)
#	#proj=tf.stack([tf.tensor_scatter_nd_add(proj[i],bposall[i],bampall[i]) for i in range(proj.shape[0])])


def gauss_project_single_fn(gausary,mx,boxsize,tytx):
	"""This exists as a function separate from the Gaussian class to better support JAX optimization. It is called by the corresponding Gaussian method.

	Generates an array containing a simple 2-D projection (interpolated delta functions) of the set of Gaussians for each of N Orientations in orts.
	gausary - a Gaussians.jax array
	mx - an Orientations object converted to a stack of 2d matrices
	tytx =  a N x 2+ vector containing an in-plane translation in unit (-0.5 - 0.5) coordinates to be applied to the set of Gaussians for each Orientation.
	boxsize in pixels. Scaling factor is equal to boxsize, such that -0.5 to 0.5 range covers the box.

	With these definitions, Gaussian coordinates are sampling-independent as long as no box size alterations are performed. That is, raw projection data
	used for comparisons should be resampled without any "clip" operations.
	"""

	proj2=[]
	shift10=jnp.array((1,0))
	shift01=jnp.array((0,1))
	shift11=jnp.array((1,1))
#	print("t1")

	xfgauss=jnp.einsum("ij,kj->ki",mx,gausary[:,:3])	# changed to ik instead of ki due to y,x ordering in tensorflow
	xfgauss+=tytx[:2]	# translation, ignore z or any other variables which might be used for per particle defocus, etc
	xfgauss=(xfgauss+0.5)*boxsize			# shift and scale both x and y the same

	xfgaussf=jnp.floor(xfgauss)
	xfgaussi=xfgaussf.astype(jnp.int32)		# integer index
	xfgaussf=xfgauss-xfgaussf				# remainder used for bilinear interpolation

		# messy tensor math here to implement bilinear interpolation
	bamp0=gausary[:,3]*(1.0-xfgaussf[:,0])*(1.0-xfgaussf[:,1])	#0,0
	bamp1=gausary[:,3]*(xfgaussf[:,0])*(1.0-xfgaussf[:,1])		#1,0
	bamp2=gausary[:,3]*(xfgaussf[:,0])*(xfgaussf[:,1])		#1,1
	bamp3=gausary[:,3]*(1.0-xfgaussf[:,0])*(xfgaussf[:,1])		#0,1
	bampall=jnp.concat([bamp0,bamp1,bamp2,bamp3],axis=0)  			# TODO: this would be ,1 with the loop subsumed
	bposall=jnp.concat([xfgaussi,xfgaussi+shift10,xfgaussi+shift11,xfgaussi+shift01],axis=0).transpose() # TODO: this too

		# note: tried this using advanced indexing, but JAX wouldn't accept the syntax for 2-D arrays
	proj=jnp.zeros((boxsize,boxsize),dtype=jnp.float32)
	proj=proj.at[bposall[0],bposall[1]].add(bampall, mode="drop")		# projection
	return proj

gauss_project_simple_fn=jax.jit(jax.vmap(gauss_project_single_fn, in_axes=[None, 2, None, 0]), static_argnames=["boxsize"])


#def gauss_project_ctf_fn(gausary,mx,ctfary,dfmin,dfmax,dfstep,boxsize,tytx):
#	"""This exists as a function separate from the Gaussian class to better support JAX optimization. It is called by the corresponding Gaussian method.
#
#	Generates an array containing a simple 2-D projection (interpolated delta functions) of the set of Gaussians for each of N Orientations in orts.
#	gausary - a Gaussians.jax array
#	mx - an Orientations object converted to a stack of 2d matrices
#	ctfary - A jax array of ctf corrections for defocuses given by dfmin, dfmax, and dfstep
#	dfmin - the min defocus value used in ctfary
#	dfmax - the max defocus value used in ctfary
#	dfstep - the defocus step between images in ctfary
#	tytx =  a N x 2+ vector containing an in-plane translation in unit (-0.5 - 0.5) coordinates to be applied to the set of Gaussians for each Orientation.
#	boxsize in pixels. Scaling factor is equal to boxsize, such that -0.5 to 0.5 range covers the box.
#
#	With these definitions, Gaussian coordinates are sampling-independent as long as no box size alterations are performed. That is, raw projection data
#	used for comparisons should be resampled without any "clip" operations.
#	"""
#	proj2=[]
#
#	# iterate over projections
#	# TODO - at some point this outer loop should be converted to a tensor axis for better performance
#	# note that the mx dimensions have N as the 3rd not 1st component!
#	# TODO: I think we can just use jax.vmap() instead of having this whole function. It is made to vectorize a function and works well with jit
#	gpcsf=jax.jit(gauss_project_ctf_single_fn,static_argnames=["boxsize"])
#
#	for j in range(mx.shape[2]):
#		proj2.append(gpcsf(gausary,mx[:,:,j],ctfary,dfmin,dfstep,boxsize,tytx[j]))
#
#	return jnp.stack(proj2)

def gauss_project_ctf_single_fn(gausary,mx,ctfary,dfmin,dfstep,boxsize,tytx):
	"""This exists as a function separate from the Gaussian class to better support JAX optimization. It is called by the corresponding Gaussian method.

	Generates an array containing a simple 2-D projection (interpolated delta functions) of the set of Gaussians for each of N Orientations in orts.
	gausary - a Gaussians.jax array
	mx - an Orientations object converted to a stack of 2d matrices
	tytx =  a N x 2+ vector containing an in-plane translation in unit (-0.5 - 0.5) coordinates to be applied to the set of Gaussians for each Orientation.
	boxsize in pixels. Scaling factor is equal to boxsize, such that -0.5 to 0.5 range covers the box.

	With these definitions, Gaussian coordinates are sampling-independent as long as no box size alterations are performed. That is, raw projection data
	used for comparisons should be resampled without any "clip" operations.
	"""
	shift10=jnp.array((1,0))
	shift01=jnp.array((0,1))
	shift11=jnp.array((1,1))

	xfgauss=jnp.einsum("ij,kj->ki",mx,gausary[:,:3])	# changed to ik instead of ki due to y,x ordering in tensorflow
	xfgauss+=tytx[:2]	# translation, ignore z or any other variables which might be used for per particle defocus, etc
	xfgauss=(xfgauss+0.5)*boxsize			# shift and scale both x and y the same

	xfgaussf=jnp.floor(xfgauss)
	xfgaussi=xfgaussf.astype(jnp.int32)		# integer index
	xfgaussf=xfgauss-xfgaussf				# remainder used for bilinear interpolation

	# messy tensor math here to implement bilinear interpolation
	bamp0=gausary[:,3]*(1.0-xfgaussf[:,0])*(1.0-xfgaussf[:,1])	#0,0
	bamp1=gausary[:,3]*(xfgaussf[:,0])*(1.0-xfgaussf[:,1])		#1,0
	bamp2=gausary[:,3]*(xfgaussf[:,0])*(xfgaussf[:,1])			#1,1
	bamp3=gausary[:,3]*(1.0-xfgaussf[:,0])*(xfgaussf[:,1])		#0,1
	bampall=jnp.concat([bamp0,bamp1,bamp2,bamp3],axis=0)  			# TODO: this would be ,1 with the loop subsumed
	bposall=jnp.concat([xfgaussi,xfgaussi+shift10,xfgaussi+shift11,xfgaussi+shift01],axis=0).transpose() # TODO: this too

	# note: tried this using advanced indexing, but JAX wouldn't accept the syntax for 2-D arrays
	proj=jnp.zeros((boxsize,boxsize),dtype=jnp.float32)
	proj=proj.at[bposall[0],bposall[1]].add(bampall, mode="drop")
	return jit_apply_ctf(ctfary, proj, tytx[2], dfmin, dfstep)

gauss_project_ctf_fn=jax.jit(jax.vmap(gauss_project_ctf_single_fn, in_axes=[None, 2, None, None, None, None, 0]) ,static_argnames=["boxsize"])

#def gauss_project_layered_ctf_fn(gausary,mx,ctfary,boxsize,dfmin,dfmax,dfstep,apix,tytx):
#	proj2=[]
#
#	# iterate over projections
#	# TODO - at some point this outer loop should be converted to a tensor axis for better performance
#	# note that the mx dimensions have N as the 3rd not 1st component!
#	gplcsf=jax.jit(gauss_project_layered_ctf_single_fn,static_argnames=["boxsize","apix","dfstep"])
#
#	for j in range(mx.shape[2]):
#		proj2.append(gplcsf(gausary,mx[:,:,j],ctfary,dfmin,dfmax,dfstep,apix,boxsize,tytx[j]))
#
#	return jnp.stack(proj2)

def gauss_project_layered_ctf_single_fn(gausary, mx, ctfary,dfmin,dfstep,apix,boxsize,tytx):
	shift10=jnp.array((1,0))
	shift01=jnp.array((0,1))
	shift11=jnp.array((1,1))

	boxstep=dfstep*10000.0/apix
#	offset=ceil((boxsize*jnp.sqrt(3)-boxstep)/(2*boxstep))
	offset = ceil((boxsize*1.733-boxstep)/(2*boxstep))

	#TODO: make sure einsum is doing what I want given ordering changes from 2d -> 3d mx. In tf I specified axes=[-1]
	xfgauss=jnp.fliplr(jnp.einsum("ij,kj->ki",mx,gausary[:,:3]))	# changed to ik instead of ki due to y,x ordering in tensorflow
	xfgaussz = xfgauss[:,0]
	xfgauss = xfgauss[:,1:]+tytx[:2]	# translation, ignore z or any other variables which might be used for per particle defocus, etc
	xfgauss=(xfgauss+0.5)*boxsize			# shift and scale both x and y the same
	xfgaussz *= boxsize
	xfgaussz = jnp.round(xfgaussz/boxstep).astype(jnp.int32)
#	xfgaussz = jnp.clip(xfgaussz, 0, 2*offset) # Clip z such that if any Gaussians are outside the expected range but in the box they will be in the outermost sections

	xfgaussf=jnp.floor(xfgauss)
	xfgaussi=xfgaussf.astype(jnp.int32)		# integer index
	xfgaussf=xfgauss-xfgaussf				# remainder used for bilinear interpolation

	# messy tensor math here to implement bilinear interpolation
	bamp0=gausary[:,3]*(1.0-xfgaussf[:,0])*(1.0-xfgaussf[:,1])	#0,0
	bamp1=gausary[:,3]*(xfgaussf[:,0])*(1.0-xfgaussf[:,1])		#1,0
	bamp2=gausary[:,3]*(xfgaussf[:,0])*(xfgaussf[:,1])			#1,1
	bamp3=gausary[:,3]*(1.0-xfgaussf[:,0])*(xfgaussf[:,1])		#0,1
	bampall=jnp.concat([bamp0,bamp1,bamp2,bamp3],axis=0)  			# TODO: this would be ,1 with the loop subsumed
	bposall=jnp.concat([xfgaussi,xfgaussi+shift10,xfgaussi+shift11,xfgaussi+shift01],axis=0).transpose() # TODO: this too

	# note: tried this using advanced indexing, but JAX wouldn't accept the syntax for 2-D arrays
	proj=jnp.zeros((2*offset+1,boxsize,boxsize),dtype=jnp.float32)
	proj=proj.at[jnp.tile(xfgaussz,4),bposall[0],bposall[1]].add(bampall, mode="drop")
	proj = jit_apply_ctf(ctfary, proj, jnp.linspace(tytx[2]-offset*dfstep, tytx[2]+offset*dfstep, 2*offset+1), dfmin, dfstep)
	return jnp.sum(proj, axis=0)

gauss_project_layered_ctf_fn=jax.jit(jax.vmap(gauss_project_layered_ctf_single_fn, in_axes=[None, 2, None, None, None, None, None, 0]) ,static_argnames=["boxsize","apix","dfstep"])

def gauss_volume_fn(gausary,boxsize,zsize):
	"""This exists as a function separate from the Gaussian class to better support JAX optimization. It is called by the corresponding Gaussian method."""

#		xfgauss=tf.reverse((gausary[:,:3]+(0.5,0.5,zaspect))*boxsize,[-1])		# shift and scale both x and y the same, reverse handles the XYZ -> ZYX EMData->Tensorflow issue
	zaspect=zsize/(2.0*boxsize)
	xfgauss=jnp.flip((gausary[:,:3]+jnp.array((0.5,0.5,zaspect)))*boxsize,-1)		# shift and scale both x and y the same, reverse handles the XYZ -> ZYX EMData->Tensorflow issue
	pos_mask = jnp.all(jnp.logical_and(xfgauss>0.0, xfgauss<boxsize-1.0001),axis=1).astype(float)
	xfgauss=jnp.clip(xfgauss,0.0,boxsize-1.0001)

	xfgaussf=jnp.floor(xfgauss)
	xfgaussi=xfgaussf.astype(jnp.int32)	# integer index
	xfgaussf=xfgauss-xfgaussf				# remainder used for bilinear interpolation
	xfamps=gausary[:,3]*pos_mask

	shift001=jnp.array((0,0,1))
	shift010=jnp.array((0,1,0))
	shift011=jnp.array((0,1,1))
	shift100=jnp.array((1,0,0))
	shift101=jnp.array((1,0,1))
	shift110=jnp.array((1,1,0))
	shift111=jnp.array((1,1,1))


	# messy trilinear interpolation
	bamp000=xfamps*(1.0-xfgaussf[:,0])*(1.0-xfgaussf[:,1])*(1.0-xfgaussf[:,2])
	bamp001=xfamps*(1.0-xfgaussf[:,0])*(1.0-xfgaussf[:,1])*(    xfgaussf[:,2])
	bamp010=xfamps*(1.0-xfgaussf[:,0])*(    xfgaussf[:,1])*(1.0-xfgaussf[:,2])
	bamp011=xfamps*(1.0-xfgaussf[:,0])*(    xfgaussf[:,1])*(    xfgaussf[:,2])
	bamp100=xfamps*(    xfgaussf[:,0])*(1.0-xfgaussf[:,1])*(1.0-xfgaussf[:,2])
	bamp101=xfamps*(    xfgaussf[:,0])*(1.0-xfgaussf[:,1])*(    xfgaussf[:,2])
	bamp110=xfamps*(    xfgaussf[:,0])*(    xfgaussf[:,1])*(1.0-xfgaussf[:,2])
	bamp111=xfamps*(    xfgaussf[:,0])*(    xfgaussf[:,1])*(    xfgaussf[:,2])
	bampall=jnp.concat([bamp000,bamp001,bamp010,bamp011,bamp100,bamp101,bamp110,bamp111],axis=0)
	bposall=jnp.concat([xfgaussi,xfgaussi+shift001,xfgaussi+shift010,xfgaussi+shift011,xfgaussi+shift100,xfgaussi+shift101,xfgaussi+shift110,xfgaussi+shift111],axis=0).transpose()

	vol=jnp.zeros((zsize,boxsize,boxsize),dtype=jnp.float32).at[bposall[0],bposall[1],bposall[2]].add(bampall)

	return vol

def jit_apply_ctf(ctfary, projs, dfary, dfmin, dfstep):
	"""jitable version of apply_ctf function in CTFStack class. Called in projection code"""
	jax.debug.print("defocus value(s): {df}\nindex in array: {i}", df=dfary, i=jnp.round((dfary-dfmin)/dfstep).astype(jnp.int32))
	return jnp.fft.irfft2(jnp.fft.rfft2(projs) * ctfary[jnp.round((dfary-dfmin)/dfstep).astype(jnp.int32)])

def jax_to_mx2d(ortary,swapxy=False):
	"""Returns the current set of orientations as a 2 x 3 x N matrix which will transform a set of 3-vectors to a set of
	2-vectors, ignoring the resulting Z component. Typically used with Gaussians to generate projections.

	To apply to a set of vectors:
	mx=self.to_mx2d()
	vecs=tf.constant(((1,0,0),(0,1,0),(0,0,1),(2,2,2),(1,1,0),(0,1,1)),dtype=tf.float32)

	tf.transpose(tf.matmul(mx[:,:,0],tf.transpose(vecs)))
	or
	tf.einsum("ij,kj->ki",mx[:,:,0],vecs)

	if swapxy is set, then the input vector is XYZ, but output is YX. This corrects for the fact that Gaussians are XYZA,
	but images are YX.
	"""
	# Find the zero rotation values
	l=jnp.linalg.norm(ortary, axis=1)
	is_zero = l < 1e-12
	# Filter them out so there won't be nans when I take the gradient
	l = jnp.linalg.norm(jnp.where(is_zero[:,None], jnp.ones_like(ortary), ortary), axis=1)

	# jnp.where calls from here on are to fix second derivative at the zero (fixed) values
	w = jnp.where(is_zero, -pi*pi*(l**2)/2, jnp.cos(pi*l)) # cos "real" component of quaternion
	s = jnp.where(is_zero, pi*pi*pi*(l**2)/6, jnp.sin(-pi*l)/l)
	q=jnp.transpose(ortary)*s		# transpose makes the vectorized math below work properly
	if swapxy :
		mx=jnp.array(((2*q[0]*q[1]+2*q[2]*w,1-(2*q[0]*q[0]+2*q[2]*q[2]),2*q[1]*q[2]-2*q[0]*w),
		(1-2*(q[1]*q[1]+q[2]*q[2]),2*q[0]*q[1]-2*q[2]*w,2*q[0]*q[2]+2*q[1]*w)))
	else:
		mx=jnp.array(((1-2*(q[1]*q[1]+q[2]*q[2]),2*q[0]*q[1]-2*q[2]*w,2*q[0]*q[2]+2*q[1]*w),
		(2*q[0]*q[1]+2*q[2]*w,1-(2*q[0]*q[0]+2*q[2]*q[2]),2*q[1]*q[2]-2*q[0]*w)))
	return jnp.array(mx)

JAXDEV=jax.devices()[0]
def jax_set_device(dev=0,maxmem=4096):
	"""Sets maximum memory for a specific Tensorflow device and returns a device to use with "with:"
	dev - GPU number or -1 for CPU (CPU doesn't actually permit memory size allocation)
	maxmem - maximum memory to allocate in megabytes

	dev=tf_set_device(gpuid,6144)
	with dev:
		# tensorflow operations, "with" block optional
	"""
	global JAXDEV
	if dev<0 :
		raise Exception("CPU not currently supported on this JAX installation")
		# pdevice=tf.config.list_physical_devices('CPU')[0]
		# tf.config.set_logical_device_configuration(pdevice,[tf.config.LogicalDeviceConfiguration()])
		# return tf.device('/CPU:0')
	else:
		JAXDEV=jax.devices()[dev]
		# pdevice=tf.config.list_physical_devices('GPU')[dev]
		# tf.config.set_logical_device_configuration(pdevice,[tf.config.LogicalDeviceConfiguration(memory_limit=maxmem)])
		# return tf.device(f'/GPU:{dev}')

def from_jax(jaxtensor,stack=False):
	"""Convert a specified tensor to an EMData object
	If stack is set, then the first axis of the tensor will be unpacked to form a list. ie a 3D tensor would become a list of 2D EMData objects"""

	if stack:
		return [EMNumPy.numpy2em(np.array(jaxtensor[i])) for i in range(jaxtensor.shape[0])]
	return EMNumPy.numpy2em(np.array(jaxtensor))

def to_jax(emdata):
	"""Convert a specified EMData object or list of EMData objects into a JAX tensor. The tensor is immutable."""

	if isinstance(emdata,EMData):
		return jnp.array(EMNumPy.em2numpy(emdata))

	if isinstance(emdata,list) or isinstance(emdata,tuple):
		npstack=np.stack([to_numpy(im) for im in emdata],axis=0)
		return jnp.array(npstack)

def jax_fft2d(imgs):
	if isinstance(imgs,EMData) or ((isinstance(imgs,list) or isinstance(imgs,tuple)) and isinstance(imgs[0],EMData)): imgs=to_jax(imgs)

	if imgs.dtype==jnp.complex64: raise Exception("Data type must be real")

	return jnp.fft.rfft2(imgs)

def jax_fft3d(imgs):
	if isinstance(imgs,EMData) or ((isinstance(imgs,list) or isinstance(imgs,tuple)) and isinstance(imgs[0],EMData)): imgs=to_jax(imgs)

	if imgs.dtype==jnp.complex64: raise Exception("Data type must be real")

	return jnp.fft.rfftn(imgs,axes=(-3,-2,-1))

def jax_ift2d(imgs):
	if isinstance(imgs,EMData) or ((isinstance(imgs,list) or isinstance(imgs,tuple)) and isinstance(imgs[0],EMData)): imgs=to_jax(imgs)

	if imgs.dtype!=jnp.complex64: raise Exception("Data type must be complex")

	return jnp.fft.irfft2(imgs)

def jax_ift3d(imgs):
	if isinstance(imgs,EMData) or ((isinstance(imgs,list) or isinstance(imgs,tuple)) and isinstance(imgs[0],EMData)): imgs=to_jax(imgs)

	if imgs.dtype!=jnp.complex64: raise Exception("Data type must be complex")

	return jnp.fft.irfftn(imgs,axes=[-3,-2,-1])

POF2D=None
def jax_phaseorigin2d(imgs):
	global POF2D
	if len(imgs.shape)==3: shp=imgs.shape[1:]
	else: shp=imgs.shape

	if POF2D is None or shp!=POF2D.shape:
#		POF2D=jnp.fromfunction(lambda y,x: ((x+y)%2)*-2+1,shp,dtype=jnp.complex64)
		POF2D=jnp.fromfunction(lambda y,x: ((x+y)%2)*-2+1,shp,dtype=jnp.int32)

	return imgs*POF2D

POF3D=None
def jax_phaseorigin3d(imgs):
	global POF3D
	if len(imgs.shape)==4: shp=imgs.shape[1:]
	else: shp=imgs.shape

	if POF3D is None or shp!=POF3D.shape:
		#POF3D=jnp.fromfunction(lambda z,y,x: ((z+x+y)%2)*-2+1,shp,dtype=jnp.complex64)
		POF3D=jnp.fromfunction(lambda z,y,x: ((z+x+y)%2)*-2+1,shp,dtype=jnp.int32)

	return imgs*POF3D

def jax_gaussfilt_2d(boxsize,halfwidth):
	"""create a (multiplicative) Gaussian lowpass filter for boxsize with halfwidth (0.5=Nyquist)"""
	coef=-1.0/(halfwidth*boxsize)**2
	r2img=rad2_img(boxsize)
	filt=jnp.exp(r2img*coef)

	return filt


def jax_downsample_2d(imgs,newx,stack=False):
	"""Fourier downsamples a tensorflow 2D image or stack of 2D images (similar to math.fft.resample processor conceptually)
	return will always be a stack (3d tensor) even if the first dimension is 1
	passed image/stack may be real or complex (FFT), return is always complex!
	final image will be a square/cube with the (real space) size nx on all axes. Should not be used to downsample rectangular images.
	newx specifies the real-space image size after downsampling, MUST be even, and the input image must have even dimensions in real space
	note that complex conjugate relationships aren't enforced in the cropped Fourier volume in redundant locations
	"""

	if newx%2!=0 : raise Exception("newx must be an even number")

	if isinstance(imgs,EMData) or ((isinstance(imgs,list) or isinstance(imgs,tuple)) and isinstance(imgs[0],EMData)): imgs=to_jax(imgs)

	if imgs.dtype!=jnp.complex64: imgs=jnp.fft.rfft2(imgs)

	if imgs.ndim==2: imgs=jnp.expand_dims(imgs,0)	# we need a 3 rank tensor

	# Note: tried replacing the 2-step process with a single 2-D mask, but that led to a flattened intermediate,
	# which led to reshaping, so wound up back here for transparency, and (maybe) speed
	imgs=jnp.concatenate((imgs[:,:newx//2,:],imgs[:,imgs.shape[1]-newx//2:,:]),axis=1)
	return imgs[:,:,:newx//2+1]

def jax_downsample_3d(imgs,newx,stack=False):
	"""Fourier downsamples a tensorflow 3D image or stack of 3D images (similar to math.fft.resample processor conceptually)
	return will always be a stack (3d tensor) even if the first dimension is 1
	passed image/stack may be real or complex (FFT), return is always complex!
	final image will be a square/cube with the size nx on all axes. Should not be used to downsample rectangular images.
	newx specifies the real-space image size after downsampling,MUST be even, and the input image must have even dimensions in real space
	note that complex conjugate relationships aren't enforced in the cropped Fourier volume
	"""

	if newx%2!=0 : raise Exception("newx must be an even number")

	if isinstance(imgs,EMData) or ((isinstance(imgs,list) or isinstance(imgs,tuple)) and isinstance(imgs[0],EMData)): imgs=to_jax(imgs)

	if imgs.dtype!=jnp.complex64: imgs=jnp.fft.rfftn(imgs,axes=(-3,-2,-1))

	if imgs.ndim==3: imgs=jnp.expand_dims(imgs,0)	# we need a 3 rank tensor

	imgs=jnp.concatenate((imgs[:,:newx//2,:,:],imgs[:,imgs.shape[1]-newx//2:,:,:]),axis=1)		#Z
	imgs=jnp.concatenate((imgs[:,:,:newx//2,:],imgs[:,:,imgs.shape[1]-newx//2:,:]),axis=2)		#Y
	return imgs[:,:,:,:newx//2+1]		# X

	# cropz=lax.gather(imgs,jnp.concatenate((jnp.arange(newx//2),jnp.arange(imgs.shape[1]-newx//2,imgs.shape[1]))),axis=1)
	# cropy=lax.gather(cropz,jnp.concatenate((jnp.arange(newx//2),jnp.arange(imgs.shape[2]-newx//2,imgs.shape[1]))),axis=2)
	# return cropy[:,:,:,:newx//2+1]

# def tf_ccf_2d(ima,imb):
# 	"""Compute the cross correlation between a stack of 2-D images (ima) and either a single 2-D image or a 1-1 stack of 2-D images (imb)"""
#
# 	if ima.dtype!=jnp.complex64 or imb.dtype!=jnp.complex64 : raise Exception("tf_frc requires FFTs")

FRC_RADS={}		# dictionary (cache) of constant tensors of size ny/2+1,ny containing the integer Fourier radius to each point in the image
def rad_img_int(ny):
	global FRC_RADS
	try: return FRC_RADS[ny]
	except:
		rad_img=jnp.array(jnp.vstack((jnp.fromfunction(lambda y,x: jnp.int32(jnp.hypot(x,y)),(ny//2,ny//2+1)),jnp.fromfunction(lambda y,x: jnp.int32(jnp.hypot(x,ny//2-y)),(ny//2,ny//2+1)))))
#		rad_img=jnp.expand_dims(jnp.array(jnp.vstack((jnp.fromfunction(lambda y,x: jnp.int32(jnp.hypot(x,y)),(ny//2,ny//2+1)),jnp.fromfunction(lambda y,x: jnp.int32(jnp.hypot(x,ny//2-y)),(ny//2,ny//2+1))))),2)
		FRC_RADS[ny]=rad_img
		return rad_img

GEN_RAD2={}		# dictionary (cache) of constant tensors of size ny/2+1,ny containing the floating point radius at each Fourier pixel location
def rad2_img(ny):
	"""Returns a complex tensor ny/2+2,ny containing the (real value) Fourier radius**2 (squared) in each pixel location. Tensors for a
given size are cached for reuse. """
	global GEN_RAD2
	try: return GEN_RAD2[ny]
	except:
		rad2_img=jnp.array(jnp.vstack((jnp.fromfunction(lambda y,x: jnp.complex64(x**2+y**2),(ny//2,ny//2+1)),jnp.fromfunction(lambda y,x: jnp.complex64((x**2+(ny//2-y)**2)),(ny//2,ny//2+1)))))
		GEN_RAD2[ny]=rad2_img
		return rad2_img

FSC_RADS={}		# dictionary (cache) of constant tensors of size ny/2+1,ny containing the integer Fourier radius to each point in the image
def rad_vol_int(ny):
	global FSC_RADS
	try: return FSC_RADS[ny]
	except:
		rad_img=jnp.array(np.fromfunction(lambda z,y,x: np.int32(np.sqrt((x**2+np.min((y,ny-y),axis=0)**2+np.min((z,ny-z),axis=0)**2))),(ny,ny,ny//2+1)))
		FSC_RADS[ny]=rad_img
		return rad_img

#FRC_NORM={}		# dictionary (cache) of constant tensors of size ny/2*1.414 (we don't actually need this for anything)
#TODO iterating over the images is handled with a python for loop. This may not be taking great advantage of the GPU (just don't know)
# two possible approaches would be to add an extra dimension to rad_img to cover image number, and handle the scatter_nd as a single operation
# or to try making use of DataSet. I started a DataSet implementation, but decided it added too much design complexity
def jax_frc(ima,imb,avg=0,weight=1.0,minfreq=0):
	"""Computes the pairwise FRCs between two stacks of complex images. Returns a list of 1D FSC tensors or if avg!=0
	then the average of the first 'avg' values. If -1, averages through Nyquist. Weight permits a frequency based weight
	(only for avg>0): 1-2 will upweight low frequencies, 0-1 will upweight high frequencies"""
	if ima.dtype!=jnp.complex64 or imb.dtype!=jnp.complex64 : raise Exception("jax_frc requires FFTs")
#	if tf.rank(ima)<3 or tf.rank(imb)<3 or ima.shape != imb.shape: raise Exception("tf_frc works on stacks of FFTs not individual images, and the shape of both inputs must match")

	global FRC_RADS
#	global FRC_NORM		# we don't actually need this unless we want to compute uncertainties (number of points at each radius)
	ny=ima.shape[1]
	nimg=ima.shape[0]
	nr=int(ny*0.70711)+1	# max radius we consider
	rad_img=rad_img_int(ny)
#	try:
	imar=jnp.real(ima) # if you do the dot product with complex math the processor computes the cancelling cross-terms. Want to avoid the waste
	imai=jnp.imag(ima)
	imbr=jnp.real(imb)
	imbi=jnp.imag(imb)

	imabr=imar*imbr		# compute these before squaring for normalization
	imabi=imai*imbi

	imar=imar*imar		# just need the squared versions, not the originals now
	imai=imai*imai
	imbr=imbr*imbr
	imbi=imbi*imbi
	# except:
	# 	raise Exception(f"failed in FRC with sizes {ima.shape} {imb.shape} {imar.shape} {imbr.shape}")

	frc=[]
	zero=jnp.zeros([nr])
	for i in range(nimg):
		cross=zero.at[rad_img].add(imabr[i]+imabi[i])
		aprd=zero.at[rad_img].add(imar[i]+imai[i])
		bprd=zero.at[rad_img].add(imbr[i]+imbi[i])
		frc.append(cross/jnp.sqrt(aprd*bprd))

	frc=jnp.stack(frc)
	if avg>len(frc[0]): avg=-1
	if avg>0:
		frc=jnp.stack(frc)
		if weight!=1.0:
			w=jnp.linspace(weight,2.0-weight,nr)
			frc=frc*w
		return frc[:,minfreq:avg].mean()
#		return tf.math.reduce_mean(frc[:,minfreq:avg],1)
	elif avg==-1: return frc.mean(1)
#	elif avg==-1: return tf.math.reduce_mean(frc,1)
	else: return frc

def jax_frc_allvs1(ima,imb,avg=0,weight=1.0,minfreq=0):
	"""Computes the pairwise FRCs between a stack of complex images and a single image. Returns a list of 1D FSC tensors or if avg!=0
	then the average of the first 'avg' values. If -1, averages through Nyquist. Weight permits a frequency based weight
	(only for avg>0): 1-2 will upweight low frequencies, 0-1 will upweight high frequencies"""
	if ima.dtype!=jnp.complex64 or imb.dtype!=jnp.complex64 : raise Exception("jax_frc requires FFTs")
#	if tf.rank(ima)<3 or tf.rank(imb)<3 or ima.shape != imb.shape: raise Exception("tf_frc works on stacks of FFTs not individual images, and the shape of both inputs must match")

	global FRC_RADS
#	global FRC_NORM		# we don't actually need this unless we want to compute uncertainties (number of points at each radius)
	ny=ima.shape[1]
	nimg=ima.shape[0]
	nr=int(ny*0.70711)+1	# max radius we consider
	rad_img=rad_img_int(ny)
#	try:
	imar=jnp.real(ima) # if you do the dot product with complex math the processor computes the cancelling cross-terms. Want to avoid the waste
	imai=jnp.imag(ima)
	imbr=jnp.real(imb)
	imbi=jnp.imag(imb)

	imabr=imar*imbr		# compute these before squaring for normalization
	imabi=imai*imbi

	imar=imar*imar		# just need the squared versions, not the originals now
	imai=imai*imai
	imbr=imbr*imbr
	imbi=imbi*imbi
	# except:
	# 	raise Exception(f"failed in FRC with sizes {ima.shape} {imb.shape} {imar.shape} {imbr.shape}")

	frc=[]
	zero=jnp.zeros([nr])
	for i in range(nimg):
		cross=zero.at[rad_img].add(imabr[i]+imabi[i])
		aprd=zero.at[rad_img].add(imar[i]+imai[i])
		bprd=zero.at[rad_img].add(imbr+imbi)
		frc.append(cross/jnp.sqrt(aprd*bprd))

	frc=jnp.stack(frc)
	if avg>len(frc[0]): avg=-1
	if avg>0:
		frc=jnp.stack(frc)
		if weight!=1.0:
			w=jnp.linspace(weight,2.0-weight,nr)
			frc=frc*w
		return frc[:,minfreq:avg].mean()
#		return tf.math.reduce_mean(frc[:,minfreq:avg],1)
	elif avg==-1: return frc.mean(1)
#	elif avg==-1: return tf.math.reduce_mean(frc,1)
	else: return frc

def __jax_frc_jit(ima,imb,weight=1.0,minfreq=0,frc_Z=-3):
	"""Simplified jax_frc with fewer options to permit JIT compilation. Computes averaged FRCs to ny//2. Note that rad_img_int(ny) MUST
	be called with the appropriate size prior to using this function!"""
	global FRC_RADS

	ny=ima.shape[1]
	nimg=ima.shape[0]
	nr=int(ny*0.70711)+1	# max radius we consider
	rad_img=FRC_RADS[ny]	# NOTE: This is unsafe to permit JIT, appropriate FRC_RADS MUST be precomputed to use this

	imar=jnp.real(ima) # if you do the dot product with complex math the processor computes the cancelling cross-terms. Want to avoid the waste
	imai=jnp.imag(ima)
	imbr=jnp.real(imb)
	imbi=jnp.imag(imb)

	imabr=imar*imbr		# compute these before squaring for normalization
	imabi=imai*imbi

	imar=imar*imar		# just need the squared versions, not the originals now
	imai=imai*imai
	imbr=imbr*imbr
	imbi=imbi*imbi

	frc=[]
	zero=jnp.zeros([nr])
	for i in range(nimg):
		cross=zero.at[rad_img].add(imabr[i]+imabi[i])
		aprd=zero.at[rad_img].add(imar[i]+imai[i])
		bprd=zero.at[rad_img].add(imbr[i]+imbi[i])
		frc.append(cross/jnp.sqrt(aprd*bprd))

	frc=jnp.stack(frc)
	w=jnp.linspace(weight,2.0-weight,nr)
#	frc=frc*w
	ret=jax.lax.dynamic_slice(frc, (0,minfreq), (nimg,ny//2-minfreq-1)).mean(axis=1) # average over frequencies # TODO: This has shape ny//2 not stop at ny//2. Same as nimg but that seems fine
	return jnp.clip(ret,ret.mean()-ret.std()*frc_Z,1.0).mean()
#	return jnp.square(jnp.clip(ret,0.0,1.0)).mean()   # Experimental to bias gradients towards better FRCs
#	return jnp.pow(jnp.clip(ret,0.0,1.0),1.5).mean()   # Experimental to bias gradients towards better FRCs

jax_frc_jit=jax.jit(__jax_frc_jit, static_argnames="minfreq")

FSC_REFS={}
def __jax_fsc_jit(ima,imb):
	"""Computes the FSC between a stack of complex volumes and a single reference volume. Returns a stack of 1D FSC curves. Note that rad_vol_int(ny) MUST be called prior to using this function"""
	if ima.dtype!=jnp.complex64 or imb.dtype!=jnp.complex64 : raise Exception("tf_fsc requires FFTs")

	global FSC_RADS

	ny=ima.shape[1]
	nimg=ima.shape[0]
	nr=int(ny*0.86602)+1	# max radius we consider
	rad_img=FSC_RADS[ny]	# NOTE: This is unsafe to permit JIT, appropriate FRC_RADS MUST be precomputed to use this

	imar=jnp.real(ima) # if you do the dot product with complex math the processor computes the cancelling cross-terms. Want to avoid the waste
	imai=jnp.imag(ima)
	imbr=jnp.real(imb)
	imbi=jnp.imag(imb)

	imabr=imar*imbr		# compute these before squaring for normalization
	imabi=imai*imbi

	imar=imar*imar		# just need the squared versions, not the originals now
	imai=imai*imai
	imbr=imbr*imbr
	imbi=imbi*imbi

	fsc=[]
	zero=jnp.zeros([nr])
	for i in range(nimg):
		cross=zero.at[rad_img].add(imabr[i]+imabi[i])
		aprd=zero.at[rad_img].add(imar[i]+imai[i])
		bprd=zero.at[rad_img].add(imbr[i]+imbi[i])
		fsc.append(cross/jnp.sqrt(aprd*bprd))

	return jnp.stack(fsc)

jax_fsc_jit=jax.jit(__jax_fsc_jit)

#TODO: consider implementing this. Issue is we need a radial filter for each image being processed, which may be expensive
# def jax_subtract_filt(ptcl,proj,masked_proj):
# 	"""ima and imb should both be jax arrays. Computes ptcl_n-filt(masked_proj_n) where filt() is computed from the FRC between ptcl_n and proj_n,
# 	to "optimally" subtract. The typical usage would be for single particle analysis to subtract a masked version of a projection from a
# 	(phase flipped) single particle image. """
#
# 	if ptcl.dtype!=jnp.complex64 or proj.dtype!=jnp.complex64 or masked_proj.dtype!=jnp.complex64: raise Exception("jax_subtract_filt requires FFTs")
# #	if tf.rank(ima)<3 or tf.rank(imb)<3 or ima.shape != imb.shape: raise Exception("tf_frc works on stacks of FFTs not individual images, and the shape of both inputs must match")
#
# 	global FRC_RADS
# #	global FRC_NORM		# we don't actually need this unless we want to compute uncertainties (number of points at each radius)
# 	ny=ima.shape[1]
# 	nimg=ima.shape[0]
# 	nr=int(ny*0.70711)+1	# max radius we consider
# 	rad_img=FRC_RADS[ny]	# NOTE: This is unsafe to permit JIT, appropriate FRC_RADS MUST be precomputed to use this
# #	try:
# 	imar=jnp.real(ptcl) # if you do the dot product with complex math the processor computes the cancelling cross-terms. Want to avoid the waste
# 	imai=jnp.imag(ptcl)
# 	imbr=jnp.real(proj)
# 	imbi=jnp.imag(proj)
#
# 	imabr=imar*imbr		# compute these before squaring for normalization
# 	imabi=imai*imbi
#
# 	imar=imar*imar		# just need the squared versions, not the originals now
# 	imai=imai*imai
# 	#imbr=imbr*imbr		# don't want to normalize
# 	#imbi=imbi*imbi
#
# 	frc=[]
# 	zero=jnp.zeros([nr])
# 	for i in range(nimg):
# 		cross=zero.at[rad_img].add(imabr[i]+imabi[i])
# 		aprd=zero.at[rad_img].add(imar[i]+imai[i])
# 		bprd=zero.at[rad_img].add(imbr[i]+imbi[i])
# 		frc.append(cross/aprd)
#
