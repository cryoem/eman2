#!/usr/bin/env python
#
# Author: Steven Ludtke, 05/19/2023 (sludtke@bcm.edu)
# Copyright (c) 2000-2023 Baylor College of Medicine
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
This module contains EMAN3-specific operations. We don't generate EMAN-like aliases for monolithic tensorflow operations.

ONLY import this file if you will be working with tensorflow in your program, otherwise the tensorflow initialization may add unreasonable startup delays

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
import tensorflow as tf
import numpy as np

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
		self.tytx=np.zeros((n,2),dtype=np.float32)		# image shifts in absolute [-0.5,0.5] format
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
		stack.coerce_tensor()
		if ortss is not None:
			try: self.orts[n0:n0+len(stack)]=ortss
			except: self.orts[n0:n0+len(stack)]=ortss.numpy()
		if tytxs is not None:
			try: self.tytx[n0:n0+len(stack)]=tytxs
			except: self.tytx[n0:n0+len(stack)]=tytxs.numpy()

		# we go through the images one at a time, serialze, and write to a file with a directory
		self.fp.seek(self.cloc)
		for i in range(len(stack)):
			im=stack[i]
			self.locs[n0+i]=self.cloc
			self.fp.write(tf.io.serialize_tensor(im).numpy())
			self.cloc=self.fp.tell()
			self.locs[n0+i+1]=self.cloc

		self.locked=False

	def add_orts(self,nlist,dorts=None,dtytxs=None):
		"""adds dorts and dtytxs to existing arrays at locations described by nlist, used to take a gradient step
		on a subset of the data."""
		if dorts is not None: self.orts[nlist]+=dorts
		if dtytxs is not None: self.tytx[nlist]+=dtytxs

	def set_frcs(self,nlist,frcs):
		self.frcs[nlist]=frcs


	def read(self,nlist):
		while self.locked: time.sleep(0.1)

		self.locked=True

		stack=[]
		for i in nlist:
			try:
				self.fp.seek(self.locs[i])
				stack.append(tf.io.parse_tensor(self.fp.read(self.locs[i+1]-self.locs[i]),out_type=tf.complex64))
			except:
				raise Exception(f"Error reading cache {self.filename}: {i} -> {self.locs[i]}")

		self.locked=False
		ret=EMStack2D(tf.stack(stack))
		orts=Orientations(self.orts[nlist])
		tytx=tf.constant(self.tytx[nlist])
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
	def tensor(self):
		self.coerce_tensor()
		return self._data

	@tensor.setter
	def tensor(self,value):
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
		elif isinstance(self._data,tf.Tensor): self._data=from_tf(self._data,True)
		else: raise Exception(f"Invalid data in EMStack3D: {type(self._data)}")
		self._npy_list=None		# not necessary if already EMData list

	def coerce_numpy(self):
		if isinstance(self._data,np.ndarray): return
		elif self._npy_list is not None: self._data=np.stack(self._npy_list)
		elif isinstance(self._data,list): self._data=np.stack([i.numpy() for i in self._data])
		elif isinstance(self._data,tf.Tensor): self._data=self._data.numpy()
		else: raise Exception(f"Invalid data in EMStack3D: {type(self._data)}")
		self._npy_list=None		# not necessary if already EMData list

	def coerce_tensor(self):
		if isinstance(self._data,tf.Tensor): return
		elif isinstance(self._data,list): self._data=to_tf(self._data)
		elif isinstance(self._data,np.ndarray): self._data=tf.constant(self._data)
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
			return self.tensor*tf.math.conj(target)

	def align_translate(ref,maxshift=-1):
		"""compute translational alignment of a stack of images to a same sized stack or single reference image.
		returns array of shifts the same size as the input stack. maxshift limits the maximum search area to +-maxshift
		on each axis"""
		pass

	def write_images(self,fsp=None,bits=12):
		self.coerce_emdata()
		im_write_compressed(self._data,fsp,0,bits)

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
		elif isinstance(imgs,tf.Tensor) or isinstance(imgs,np.ndarray):
			if len(imgs.shape)==3:
				imgs=tf.expand_dims(imgs,0)
			elif len(imgs.shape)!=4: raise Exception(f"EMStack3D only supports stacks of 3-D data, the provided images were {len(imgs.shape)}-D")
			self._data=imgs
			self._npy_list=None
		elif isinstance(imgs,str):
			self._data=EMData.read_images(imgs)
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
		return(self._data.shape)

	def center_clip(self,size):
		size=int(size)
		if size<1: raise Exception("center_clip(size) must be called with a positive integer")
		shp=(self.shape-size)//2
		if isinstance(self._data,list):
			newlst=[im.get_clip(Region(int(shp[1]),int(shp[2]),int(shp[2]),size,size,size)) for im in self._data]
			return EMStack2D(newlst)
		elif isinstance(self._data,np.ndarray) or isinstance(self._data,tf.Tensor):
			newary=self._data[:,shp[1]:shp[1]+size,shp[2]:shp[2]+size,shp[3]:shp[3]+size]
			return EMStack2D(newary)

	def do_fft(self,keep_type=False):
		"""Computes the FFT of each image and returns a new EMStack3D. If keep_type is not set, will convert to Tensor before computing FFT."""
		if keep_type: raise Exception("do_fft: keep_type not functional yet")
		self.coerce_tensor()

		return tf_fft3d(self._data)

	def do_ift(self,keep_type=False):
		"""Computes the IFT of each image and returns a new EMStack3D. If keep_type is not set, will convert to Tensor before computing."""
		if keep_type: raise Exception("do_ift: keep_type not functional yet")
		self.coerce_tensor()

		return tf_ift3d(self._data)

	def calc_ccf(self,target):
		"""Compute the cross correlation between each image in the stack and target, which may be a single image or another EMStack of the same size"""

		if isinstance(target,EMStack3D):
			return EMStack3D(self.tensor*tf.math.conj(target.tensor))
		elif isinstance(target,tf.Tensor):
			return EMStack3D(self.tensor*tf.math.conj(target))
		else: raise Exception("calc_ccf: target must be either EMStack2D or single Tensor")

	def downsample(self,newsize):
		"""Downsamples each image/volume in Fourier space such that its real-space dimensions after downsampling
		are "newsize" in all 2/3 dimensions. Downsampled images/volumes will be in Fourier space regardless of whether
		current stack is in real or Fourier space. This cannot be used to upsample (make images larger) and should
		not be used on rectangular images/volumes."""

		return EMStack3D(tf_downsample_3d(self.tensor,newsize))	# TODO: for now we're forcing this to be a tensor, probably better to leave it in the current format


class EMStack2D(EMStack):
	"""This class represents a stack of 2-D Images in either an EMData, NumPy or Tensorflow representation, with easy interconversion
	- Shape of the array is {N,Y,X}, as EMData, it is a list of N x EMData(X,Y)
	- All images in the stack must have the same dimensions.
	- Only one representation exists at a time. Coercing to a new type is relatively expensive in terms of time.
	- Coerce routines will insure that the required representation exists, but others will be lost to insure self-consistency.
	- The convenience method numpy_list() will return a python list of N {Z,X,Y} NumPy arrays sharing memory with the EMData objects in the EMDataStack,
	  but beware, as coercing the EMDataStack to a differnt type will invalidate these NumPy arrays, and could create a crash!

	Individual images in the stack may be accessed using emdata[n], tensor[n], numpy[n]
	"""

	def set_data(self,imgs):
		""" """
		self._xforms=None
		if imgs is None:
			self._data=None
			self._npy_list=None
		elif isinstance(imgs,EMData):
			if imgs.get_ndim()!=2: raise Exception("EMStack2D only supports 2-D data")
			self._data=[imgs]
			try: self._xforms=[imgs["xform.projection"]]
			except: pass
			self._npy_list=None
		elif isinstance(imgs,tf.Tensor) or isinstance(imgs,np.ndarray):
			if len(imgs.shape)!=3: raise Exception(f"EMStack2D only supports stacks of 2-D data, the provided images were {len(imgs.shape)}-D")
			self._data=imgs
			self._npy_list=None
		elif isinstance(imgs,str):
			self._data=EMData.read_images(imgs)
			try: self._xforms=[im["xform.projection"] for im in self._data]
			except: pass
			if self._data[0].get_ndim()!=2: raise Exception(f"EMStack2D only supports stacks of 2-D data. {imgs} is {self._data[0].get_ndim()}-D")
			self._npy_list=None
		else:
			try:
				if not isinstance(imgs[0],EMData): raise Exception(f"EMDataStack cannot be initialized with a list of {type(imgs[0])}")
				self._data=list(imgs)		# copy the list, not the elements of the list
				try: self._xforms=[im["xform.projection"] for im in self._data]
				except: pass
				self._npy_list=None
			except: raise Exception("EMDataStack may be initialized with None, a filename, an EMData object, a list/tuple of EMData objects, a NumPy array or a Tensor {N,Z,Y,X}")

	def __len__(self): return len(self._data)

	@property
	def shape(self):
		# note that the returned shape is N,Y,X regardless of representation
		if isinstance(self._data,list): return(np.array((len(self._data),self._data[0]["ny"],self._data[0]["nx"])))
		return(self._data.shape)

	@property
	def orientations(self):
		"""returns an Orientations object and tytx array for the current images if available or None if not"""
		if self._xforms is None: return None,None
		orts=Orientations()
		tytx=orts.init_from_transforms(self._xforms)
		return orts,tytx

	def center_clip(self,size):
		try: size=np.array((int(size),int(size)))
		except: size=(int(size[0]),int(size[1]))
		if size[0]<1 or size[1]<1: raise Exception("center_clip(size) must be called with a positive integer")
		shp=(self.shape[1:]-size)//2
		if isinstance(self._data,list):
			newlst=[im.get_clip(Region(int(shp[0]),int(shp[1]),int(size[0]),int(size[1]))) for im in self._data]
			return EMStack2D(newlst)
		elif isinstance(self._data,np.ndarray) or isinstance(self._data,tf.Tensor):
			newary=self._data[:,shp[0]:shp[0]+size[0],shp[1]:shp[1]+size[1]]
			return EMStack2D(newary)

	def do_fft(self,keep_type=False):
		"""Computes the FFT of each image and returns a new EMStack3D. If keep_type is not set, will convert to Tensor before computing FFT."""
		if keep_type: raise Exception("do_fft: keep_type not functional yet")
		self.coerce_tensor()

		return EMStack2D(tf_fft2d(self._data))

	def do_ift(self,keep_type=False):
		"""Computes the IFT of each image and returns a new EMStack3D. If keep_type is not set, will convert to Tensor before computing."""
		if keep_type: raise Exception("do_ift: keep_type not functional yet")
		self.coerce_tensor()

		return EMStack2D(tf_ift2d(self._data))

	def calc_ccf(self,target,center=True,offset=0):
		"""Compute the cross correlation between each image in the stack and target, which may be a single image or another EMStack of the same size.
	If center is True, will shift the phase origin so zero shift corresponds to the middle of the image"""

		if center:
			if isinstance(target,EMStack2D) and offset!=0:
				return EMStack2D(tf_phaseorigin2d(self.tensor[:-offset]*tf.math.conj(target.tensor[offset:])))
			elif isinstance(target,EMStack2D):
				return EMStack2D(tf_phaseorigin2d(self.tensor*tf.math.conj(target.tensor)))
			elif isinstance(target,tf.Tensor) and offset==0:
				return EMStack2D(tf_phaseorigin2d(self.tensor*tf.math.conj(target)))
			else: raise Exception("calc_ccf: target must be either EMStack2D or single Tensor")
		else:
			if isinstance(target,EMStack2D) and offset!=0:
				return EMStack2D(self.tensor[:-offset]*tf.math.conj(target.tensor[offset:]))
			elif isinstance(target,EMStack2D):
				return EMStack2D(self.tensor*tf.math.conj(target.tensor))
			elif isinstance(target,tf.Tensor) and offset==0:
				return EMStack2D(self.tensor*tf.math.conj(target))
			else: raise Exception("calc_ccf: target must be either EMStack2D or single Tensor")

	def convolve(self,target):
		"""Compute the convolution between each image in the stack and target, which may be a single image or another EMStack of the same size"""

		if isinstance(target,EMStack2D):
			return EMStack2D(self.tensor*target.tensor)
		elif isinstance(target,tf.Tensor):
			return EMStack2D(self.tensor*target)
		else: raise Exception("calc_ccf: target must be either EMStack2D or single Tensor")

	def downsample(self,newsize):
		"""Downsamples each image/volume in Fourier space such that its real-space dimensions after downsampling
		are "newsize" in all 2/3 dimensions. Downsampled images/volumes will be in Fourier space regardless of whether
		current stack is in real or Fourier space. This cannot be used to upsample (make images larger) and should
		not be used on rectangular images/volumes."""

		if newsize==self.shape[1]: return EMStack2D(self.tensor) # this won't copy, but since the tensor is constant should be ok?
		return EMStack2D(tf_downsample_2d(self.tensor,newsize))	# TODO: for now we're forcing this to be a tensor, probably better to leave it in the current format

	def align_translate(self,ref,maxshift=-1):
		"""compute translational alignment of a stack of images to a same sized stack (or single) of reference images.
		returns array of (dy,dx) the same size as the input stack required to bring each "this" image into alignment with "ref". maxshift limits the maximum search area to +-maxshift
		on each axis. If maxshift is unspecified -> box size //4"""

		ny,nx=self.shape[1:]
		if self.tensor.dtype==tf.complex64 :
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
		ccfsrs=tf.reshape(ccfsrc,(ccfsrc.shape[0],maxshift*2*maxshift*2))

		# The y,x coordinates of the peak location
		peaks=maxshift-tf.unravel_index(tf.argmax(ccfsrs,1),(maxshift*2,maxshift*2))

		return tf.transpose(peaks)

class Orientations():
	"""This represents a set of orientations, with a standard representation of an XYZ vector where the vector length indicates the amount
		of rotation with a length of 0.5 corresponding to 180 degrees. This form is a good representation for deep learning minimization strategies
		which conventionally span a range of ~1.0. This form can be readily interconverted to EMAN2 Transform objects or transformation matrices
		for use with Gaussians. Note that this object handles only orientation, not translation.
	"""

	def __init__(self,xyzs=None):
		"""Initialize with either the number of orientations or a N x 3 matrix"""
		if isinstance(xyzs,int):
			if xyzs<=0: self._data=None
			else: self._data=np.zeros((xyzs,3))
		else:
			try: self._data=np.array(xyzs)
			except: raise Exception("Orientations must be initialized with an integer (number of orientations) or a N x 3 numpy array")


	def __len__(self): return len(self._data)

	def __getitem__(self,key):
		"""Return the keyed Gaussian parameter, may return a tensor or numpy array. G[i] returns the 4-vector for the i'th Gaussian"""
		return self._data[key]

	def __setitem__(self,key,value):
		# if the Gaussians are a tensor, we turn it back into numpy for modification
		self.coerce_numpy()
		self._data[key]=value

	def __len__(self): return self._data.shape[0]

	def coerce_tensor(self):
		if not isinstance(self._data,tf.Tensor): self._data=tf.constant(self._data,tf.float32)

	def coerce_numpy(self):
		if isinstance(self._data,tf.Tensor): self._data=self._data.numpy()

	@property
	def tensor(self):
		self.coerce_tensor()
		return self._data

	@property
	def numpy(self):
		self.coerce_numpy()
		return self._data

	def init_from_transforms(self,xformlist):
		"""Replaces current contents of Orientations object with a list of Transform objects,
		returns tytx array with any translations (not stored within Orientations)"""
		self._data=np.zeros((len(xformlist),3))
		tytx=[]
		for i,x in enumerate(xformlist):
			r=x.get_rotation("spinvec")
			self._data[i]=(r["v1"],r["v2"],r["v3"])
			tytx.append((x.get_trans_2d()[1],x.get_trans_2d()[0]))

		return(tf.constant(tytx))

	def transforms(self,tytx=None):
		"""converts the current orientations to a list of Transform objects"""

		if tytx is not None:
			return [Transform({"type":"spinvec","v1":self._data[i][0],"v2":self._data[i][1],"v3":self._data[i][2],"tx":tytx[i][1],"ty":tytx[i][0]}) for i in range(len(self._data))]

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

		self.coerce_tensor()

		# Adding a tiny value avoids the issue with zero rotations. While it would be more correct to use a conditional
		# it is much slower, and the tiny pertutbation should not significantly impact the math.
		l=tf.norm(self._data,axis=1)+1.0e-37

		w=tf.cos(pi*l)  # cos "real" component of quaternion
		s=tf.sin(-pi*l)/l
		q=tf.transpose(self._data)*s		# transpose makes the vectorized math below work properly

		if swapxy :
			mx=tf.stack(((2*q[0]*q[1]+2*q[2]*w,1-(2*q[0]*q[0]+2*q[2]*q[2]),2*q[1]*q[2]-2*q[0]*w),
			(1-2*(q[1]*q[1]+q[2]*q[2]),2*q[0]*q[1]-2*q[2]*w,2*q[0]*q[2]+2*q[1]*w)))
		else:
			mx=tf.stack(((1-2*(q[1]*q[1]+q[2]*q[2]),2*q[0]*q[1]-2*q[2]*w,2*q[0]*q[2]+2*q[1]*w),
			(2*q[0]*q[1]+2*q[2]*w,1-(2*q[0]*q[0]+2*q[2]*q[2]),2*q[1]*q[2]-2*q[0]*w)))
		return mx

	def to_mx3d(self):
		"""Returns the current set of orientations as a 3 x 3 x N matrix which will transform a set of 3-vectors to a set of
		rotated 3-vectors, ignoring the resulting Z component. Typically used with Gaussians to generate projections.

		To apply to a set of vectors:
		mx=self.to_mx2d()
		vecs=tf.constant(((1,0,0),(0,1,0),(0,0,1),(2,2,2),(1,1,0),(0,1,1)),dtype=tf.float32)

		tf.transpose(tf.matmul(mx[:,:,0],tf.transpose(vecs)))
		or
		tf.einsum("ij,kj->ki",mx[:,:,0],vecs)"""

		self.coerce_tensor()

		# Adding a tiny value avoids the issue with zero rotations. While it would be more correct to use a conditional
		# it is much slower, and the tiny pertutbation should not significantly impact the math.
		l=tf.norm(self._data,axis=1)+1.0e-37

		w=tf.cos(pi*l)  # cos "real" component of quaternion
		s=tf.sin(pi*l)/l
		q=tf.transpose(self._data)*s		# transpose makes the vectorized math below work properly

		mx=tf.stack(((1-2*(q[1]*q[1]+q[2]*q[2]),2*q[0]*q[1]-2*q[2]*w,2*q[0]*q[2]+2*q[1]*w),
		(2*q[0]*q[1]+2*q[2]*w,1-(2*q[0]*q[0]+2*q[2]*q[2]),2*q[1]*q[2]-2*q[0]*w),
		(2*q[0]*q[2]+2*q[1]*w,2*q[1]*q[2]+2*q[0]*w,1-(2*q[0]*q[0]+2*q[1]*q[1]))))
		return mx


class Gaussians():
	"""This represents a set of Gaussians with x,y,z,amp parameters (but no width). Representation is a N x 4 numpy array or tensor (x,y,z,amp) ],
x,y,z are ~-0.5 to ~0.5 (typ) and amp is 0 to ~1. A scaling factor (value -> pixels) is applied when generating projections. """

	def __init__(self,gaus=0):
		if isinstance(gaus,int):
			if gaus<=0: self._data=None
			else: self._data=np.zeros((gaus,4))
		else:
			try: self._data=np.array(gaus)
			except: raise Exception("Gaussians must be initialized with an integer (number of Gaussians) or N x 4 matrix")

	def __getitem__(self,key):
		"""Return the keyed Gaussian parameter, may return a tensor or numpy array. G[i] returns the 4-vector for the i'th Gaussian"""
		return self._data[key]

	def __setitem__(self,key,value):
		# if the Gaussians are a tensor, we turn it back into numpy for modification
		self.coerce_numpy()
		self._data[key]=value

	def __len__(self): return len(self._data)

	def coerce_tensor(self):
		if not isinstance(self._data,tf.Tensor): self._data=tf.constant(self._data,tf.float32)

	def coerce_numpy(self):
		if isinstance(self._data,tf.Tensor): self._data=self._data.numpy()

	def add_tensor(self,tensor):
		self._data+=tensor

	@property
	def tensor(self):
		self.coerce_tensor()
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

		if isinstance(vol,tf.Tensor): emd=from_tf(vol)
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
		self.coerce_tensor()
		dups=[self._data+tf.random.normal(self._data.shape,stddev=dev) for i in range(n)]
		self._data=tf.concat(dups,0)

	def norm_filter(self,sig=0.5,rad_downweight=-1):
		"""Rescale the amplitudes so the maximum is 1, with amplitude below mean+sig*sigma removed. rad_downweight, if >0 will apply a radial linear amplitude decay beyond the specified radius to the corner of the cube. eg - 0.5 will downweight the corners. Downweighting only works if Gaussian coordinate range follows the -0.5 - 0.5 standard range for the box. """
		self.coerce_tensor()
		self._data=self._data*(1.0,1.0,1.0,1.0/tf.reduce_max(self._data[:,3]))		# "normalize" amplitudes so max amplitude is scaled to 1.0, not sure how necessary this really is
		if rad_downweight>0:
			famp=self._data[:,3]*(1.0-tf.nn.relu(tf.math.reduce_euclidean_norm(self._data[:,:3],1)-rad_downweight))
		else: famp=self._data[:,3]
		thr=tf.math.reduce_mean(famp)+sig*tf.math.reduce_std(famp)
		self._data=tf.boolean_mask(self._data,famp>thr)					# remove any gaussians with amplitude below threshold

	def project_simple(self,orts,boxsize,tytx=None):
		"""Generates a tensor containing a simple 2-D projection (interpolated delta functions) of the set of Gaussians for each of N Orientations in orts.
		orts - must be an Orientations object
		tytx =  is an (optional) N x 2+ vector containing an in-plane translation in unit (-0.5 - 0.5) coordinates to be applied to the set of Gaussians for each Orientation.
		boxsize in pixels. Scaling factor is equal to boxsize, such that -0.5 to 0.5 range covers the box.

		With these definitions, Gaussian coordinates are sampling-independent as long as no box size alterations are performed. That is, raw projection data
		used for comparisons should be resampled without any "clip" operations.
		"""
		self.coerce_tensor()

#		proj=tf.zeros((len(orts),boxsize,boxsize))		# projections
		proj=tf.zeros((boxsize,boxsize))		# projections
		proj2=[]
		mx=orts.to_mx2d(swapxy=True)

		# iterate over projections
		# TODO - at some point this outer loop should be converted to a tensor axis for better performance
		for j in range(len(orts)):
			xfgauss=tf.einsum("ij,kj->ki",mx[:,:,j],self._data[:,:3])	# changed to ik instead of ki due to y,x ordering in tensorflow
			if tytx is not None:
				#print(xfgauss.shape,tytx.shape)
				xfgauss+=tytx[j,:2]	# translation, ignore z or any other variables which might be used for per particle defocus, etc
			xfgauss=(xfgauss+0.5)*boxsize		# shift and scale both x and y the same

			xfgaussf=tf.floor(xfgauss)
			xfgaussi=tf.cast(xfgaussf,tf.int32)	# integer index
			xfgaussf=xfgauss-xfgaussf				# remainder used for bilinear interpolation

			# print(xfgaussf,xfgaussi)
			# print(xfgaussf.shape,xfgaussi.shape,self._data.shape,self._data[:,3].shape)
			# messy tensor math here to implement bilinear interpolation
			bamp0=self._data[:,3]*(1.0-xfgaussf[:,0])*(1.0-xfgaussf[:,1])	#0,0
			bamp1=self._data[:,3]*(xfgaussf[:,0])*(1.0-xfgaussf[:,1])		#1,0
			bamp2=self._data[:,3]*(xfgaussf[:,0])*(xfgaussf[:,1])			#1,1
			bamp3=self._data[:,3]*(1.0-xfgaussf[:,0])*(xfgaussf[:,1])		#0,1
			bampall=tf.concat([bamp0,bamp1,bamp2,bamp3],0)  # TODO: this would be ,1 with the loop subsumed
			bposall=tf.concat([xfgaussi,xfgaussi+(1,0),xfgaussi+(1,1),xfgaussi+(0,1)],0) # TODO: this too
			proj2.append(tf.tensor_scatter_nd_add(proj,bposall,bampall))

		return EMStack2D(tf.stack(proj2))
		#proj=tf.stack([tf.tensor_scatter_nd_add(proj[i],bposall[i],bampall[i]) for i in range(proj.shape[0])])

	def volume(self,boxsize,zaspect=0.5):
		self.coerce_tensor()

		zsize=good_size(boxsize*zaspect*2.0)
		vol=tf.zeros((zsize,boxsize,boxsize))		# output

		xfgauss=tf.reverse((self._data[:,:3]+(0.5,0.5,zaspect))*boxsize,[-1])		# shift and scale both x and y the same, reverse handles the XYZ -> ZYX EMData->Tensorflow issue

		xfgaussf=tf.floor(xfgauss)
		xfgaussi=tf.cast(xfgaussf,tf.int32)	# integer index
		xfgaussf=xfgauss-xfgaussf				# remainder used for bilinear interpolation

		# messy trilinear interpolation
		bamp000=self._data[:,3]*(1.0-xfgaussf[:,0])*(1.0-xfgaussf[:,1])*(1.0-xfgaussf[:,2])
		bamp001=self._data[:,3]*(1.0-xfgaussf[:,0])*(1.0-xfgaussf[:,1])*(    xfgaussf[:,2])
		bamp010=self._data[:,3]*(1.0-xfgaussf[:,0])*(    xfgaussf[:,1])*(1.0-xfgaussf[:,2])
		bamp011=self._data[:,3]*(1.0-xfgaussf[:,0])*(    xfgaussf[:,1])*(    xfgaussf[:,2])
		bamp100=self._data[:,3]*(    xfgaussf[:,0])*(1.0-xfgaussf[:,1])*(1.0-xfgaussf[:,2])
		bamp101=self._data[:,3]*(    xfgaussf[:,0])*(1.0-xfgaussf[:,1])*(    xfgaussf[:,2])
		bamp110=self._data[:,3]*(    xfgaussf[:,0])*(    xfgaussf[:,1])*(1.0-xfgaussf[:,2])
		bamp111=self._data[:,3]*(    xfgaussf[:,0])*(    xfgaussf[:,1])*(    xfgaussf[:,2])
		bampall=tf.concat([bamp000,bamp001,bamp010,bamp011,bamp100,bamp101,bamp110,bamp111],0)
		bposall=tf.concat([xfgaussi,xfgaussi+(0,0,1),xfgaussi+(0,1,0),xfgaussi+(0,1,1),xfgaussi+(1,0,0),xfgaussi+(1,0,1),xfgaussi+(1,1,0),xfgaussi+(1,1,1)],0)
		vol=tf.tensor_scatter_nd_add(vol,bposall,bampall)

		return EMStack3D(vol)
		#proj=tf.stack([tf.tensor_scatter_nd_add(proj[i],bposall[i],bampall[i]) for i in range(proj.shape[0])])


def tf_set_device(dev=0,maxmem=4096):
	"""Sets maximum memory for a specific Tensorflow device and returns a device to use with "with:"
	dev - GPU number or -1 for CPU (CPU doesn't actually permit memory size allocation)
	maxmem - maximum memory to allocate in megabytes

	dev=tf_set_device(gpuid,6144)
	with dev:
		# tensorflow operations, "with" block optional
	"""
	if dev<0 :
		pdevice=tf.config.list_physical_devices('CPU')[0]
		tf.config.set_logical_device_configuration(pdevice,[tf.config.LogicalDeviceConfiguration()])
		return tf.device('/CPU:0')
	else:
		pdevice=tf.config.list_physical_devices('GPU')[dev]
		tf.config.set_logical_device_configuration(pdevice,[tf.config.LogicalDeviceConfiguration(memory_limit=maxmem)])
		return tf.device(f'/GPU:{dev}')

def from_tf(tftensor,stack=False):
	"""Convert a specified tensor to an EMData object
	If stack is set, then the first axis of the tensor will be unpacked to form a list. ie a 3D tensor would become a list of 2D EMData objects"""

	if stack:
		return [EMNumPy.numpy2em(tftensor[i].numpy()) for i in range(tftensor.shape[0])]
	return EMNumPy.numpy2em(tftensor.numpy())

def to_tfvar(emdata):
	"""Convert a specified EMData object or list of EMData objects into a TensorFlow Variable. WARNING many tensorflow operations are very inefficient with Variable tensors!"""
	if isinstance(emdata,EMData):
		return tf.Variable(EMNumPy.em2numpy(emdata))

	if isinstance(emdata,list) or isinstance(emdata,tuple):
		npstack=np.stack([to_numpy(im) for im in emdata],axis=0)
		return tf.Variable(npstack)

def to_tf(emdata):
	"""Convert a specified EMData object or list of EMData objects into a Tensorflow constant tensor. The tensor is immutable, but will have much better performance for most operations."""

	if isinstance(emdata,EMData):
		return tf.constant(EMNumPy.em2numpy(emdata))

	if isinstance(emdata,list) or isinstance(emdata,tuple):
		npstack=np.stack([to_numpy(im) for im in emdata],axis=0)
		return tf.constant(npstack)

def tf_fft2d(imgs):
	if isinstance(imgs,EMData) or ((isinstance(imgs,list) or isinstance(imgs,tuple)) and isinstance(imgs[0],EMData)): imgs=to_tf(imgs)

	if imgs.dtype==tf.complex64: raise Exception("Data type must be real")

	return tf.signal.rfft2d(imgs)

def tf_fft3d(imgs):
	if isinstance(imgs,EMData) or ((isinstance(imgs,list) or isinstance(imgs,tuple)) and isinstance(imgs[0],EMData)): imgs=to_tf(imgs)

	if imgs.dtype==tf.complex64: raise Exception("Data type must be real")

	return tf.signal.rfft3d(imgs)

def tf_ift2d(imgs):
	if isinstance(imgs,EMData) or ((isinstance(imgs,list) or isinstance(imgs,tuple)) and isinstance(imgs[0],EMData)): imgs=to_tf(imgs)

	if imgs.dtype!=tf.complex64: raise Exception("Data type must be complex")

	return tf.signal.irfft2d(imgs)

def tf_ift3d(imgs):
	if isinstance(imgs,EMData) or ((isinstance(imgs,list) or isinstance(imgs,tuple)) and isinstance(imgs[0],EMData)): imgs=to_tf(imgs)

	if imgs.dtype!=tf.complex64: raise Exception("Data type must be complex")

	return tf.signal.irfft3d(imgs)

POF2D=None
def tf_phaseorigin2d(imgs):
	global POF2D
	if len(imgs.shape)==3: shp=imgs.shape[1:]
	else: shp=imgs.shape

	if POF2D is None or shp!=POF2D.shape:
		POF2D=tf.constant(np.fromfunction(lambda y,x: ((x+y)%2)*-2+1,shp),dtype=tf.complex64)

	return imgs*POF2D

POF3D=None
def tf_phaseorigin3d(imgs):
	global POF3D
	if len(imgs.shape)==3: shp=imgs.shape[1:]
	else: shp=imgs.shape

	if POF3D is None or shp!=POF3D.shape:
		POF3D=tf.constant(np.fromfunction(lambda z,y,x: ((z+x+y)%2)*-2+1,shp),dtype=tf.complex64)

	return imgs*POF3D

def tf_gaussfilt_2d(boxsize,halfwidth):
	"""create a (multiplicative) Gaussian lowpass filter for boxsize with halfwidth (0.5=Nyquist)"""
	coef=1.0/(halfwidth*boxsize)**2
	rad2_img=tf.expand_dims(tf.constant(np.vstack((np.fromfunction(lambda y,x: np.float32(x**2+y**2),(ny//2,ny//2+1)),np.fromfunction(lambda y,x: np.float32((x**2+(ny//2-y)**2)),(ny//2,ny//2+1))))),2)
	filt=tf.math.exp(rad2_img*coef)

	return filt


def tf_downsample_2d(imgs,newx,stack=False):
	"""Fourier downsamples a tensorflow 2D image or stack of 2D images (similar to math.fft.resample processor conceptually)
	return will always be a stack (3d tensor) even if the first dimension is 1
	passed image/stack may be real or complex (FFT), return is always complex!
	final image will be a square/cube with the (real space) size nx on all axes. Should not be used to downsample rectangular images.
	newx specifies the real-space image size after downsampling, MUST be even, and the input image must have even dimensions in real space
	note that complex conjugate relationships aren't enforced in the cropped Fourier volume in redundant locations
	"""

	if newx%2!=0 : raise Exception("newx must be an even number")

	if isinstance(imgs,EMData) or ((isinstance(imgs,list) or isinstance(imgs,tuple)) and isinstance(imgs[0],EMData)): imgs=to_tf(imgs)

	if imgs.dtype!=tf.complex64: imgs=tf.signal.rfft2d(imgs)

	if imgs.ndim==2: imgs=tf.expand_dims(imgs,0)	# we need a 3 rank tensor

	cropy=tf.gather(imgs,np.concatenate((np.arange(newx//2),np.arange(imgs.shape[1]-newx//2,imgs.shape[1]))),axis=1)
	return cropy[:,:,:newx//2+1]

def tf_downsample_3d(imgs,newx,stack=False):
	"""Fourier downsamples a tensorflow 3D image or stack of 3D images (similar to math.fft.resample processor conceptually)
	return will always be a stack (3d tensor) even if the first dimension is 1
	passed image/stack may be real or complex (FFT), return is always complex!
	final image will be a square/cube with the size nx on all axes. Should not be used to downsample rectangular images.
	newx specifies the real-space image size after downsampling,MUST be even, and the input image must have even dimensions in real space
	note that complex conjugate relationships aren't enforced in the cropped Fourier volume
	"""

	if newx%2!=0 : raise Exception("newx must be an even number")

	if isinstance(imgs,EMData) or ((isinstance(imgs,list) or isinstance(imgs,tuple)) and isinstance(imgs[0],EMData)): imgs=to_tf(imgs)

	if imgs.dtype!=tf.complex64: imgs=tf.signal.rfft3d(imgs)

	if imgs.ndim==3: imgs=tf.expand_dims(imgs,0)	# we need a 3 rank tensor

	cropz=tf.gather(imgs,np.concatenate((np.arange(newx//2),np.arange(imgs.shape[1]-newx//2,imgs.shape[1]))),axis=1)
	cropy=tf.gather(cropz,np.concatenate((np.arange(newx//2),np.arange(imgs.shape[2]-newx//2,imgs.shape[1]))),axis=2)
	return cropy[:,:,:,:newx//2+1]

def tf_ccf_2d(ima,imb):
	"""Compute the cross correlation between a stack of 2-D images (ima) and either a single 2-D image or a 1-1 stack of 2-D images (imb)"""

	if ima.dtype!=tf.complex64 or imb.dtype!=tf.complex64 : raise Exception("tf_frc requires FFTs")



FRC_RADS={}		# dictionary (cache) of constant tensors of size ny/2+1,ny containing the Fourier radius to each point in the image
def rad_img(ny):
	try: return FRC_RADS[ny]
	except:
		rad_img=tf.expand_dims(tf.constant(np.vstack((np.fromfunction(lambda y,x: np.int32(np.hypot(x,y)),(ny//2,ny//2+1)),np.fromfunction(lambda y,x: np.int32(np.hypot(x,ny//2-y)),(ny//2,ny//2+1))))),2)
		FRC_RADS[ny]=rad_img
		return rad_img


#FRC_NORM={}		# dictionary (cache) of constant tensors of size ny/2*1.414 (we don't actually need this for anything)
#TODO iterating over the images is handled with a python for loop. This may not be taking great advantage of the GPU (just don't know)
# two possible approaches would be to add an extra dimension to rad_img to cover image number, and handle the scatter_nd as a single operation
# or to try making use of DataSet. I started a DataSet implementation, but decided it added too much design complexity
def tf_frc(ima,imb,avg=0,weight=1.0,minfreq=0):
	"""Computes the pairwise FRCs between two stacks of complex images. imb may alternatively be a single image. Returns a list of 1D FSC tensors or if avg!=0
	then the average of the first 'avg' values. If -1, averages through Nyquist. Weight permits a frequency based weight
	(only for avg>0): 1-2 will upweight low frequencies, 0-1 will upweight high frequencies"""
	if ima.dtype!=tf.complex64 or imb.dtype!=tf.complex64 : raise Exception("tf_frc requires FFTs")
#	if tf.rank(ima)<3 or tf.rank(imb)<3 or ima.shape != imb.shape: raise Exception("tf_frc works on stacks of FFTs not individual images, and the shape of both inputs must match")

	global FRC_RADS
#	global FRC_NORM		# we don't actually need this unless we want to compute uncertainties (number of points at each radius)
	ny=ima.shape[1]
	nimg=ima.shape[0]
	nr=int(ny*0.70711)+1	# max radius we consider
	rad_img=rad_img(ny)
	try:
		imar=tf.math.real(ima) # if you do the dot product with complex math the processor computes the cancelling cross-terms. Want to avoid the waste
		imai=tf.math.imag(ima)
		imbr=tf.math.real(imb)
		imbi=tf.math.imag(imb)

		imabr=imar*imbr		# compute these before squaring for normalization
		imabi=imai*imbi

		imar=imar*imar		# just need the squared versions, not the originals now
		imai=imai*imai
		imbr=imbr*imbr
		imbi=imbi*imbi
	except:
		raise Exception(f"failed in FRC with sizes {ima.shape} {imb.shape} {imar.shape} {imbr.shape}")

	if len(imbr.shape)==3: single=False
	else: single=True
	frc=[]
	for i in range(nimg):
		zero=tf.zeros([nr])
#		print(zero.shape,rad_img.shape,imabr.shape)
		cross=tf.tensor_scatter_nd_add(zero,rad_img,imabr[i])	#start with zero when we add the real component
		cross=tf.tensor_scatter_nd_add(cross,rad_img,imabi[i])	#add the imaginary component to the real

		aprd=tf.tensor_scatter_nd_add(zero,rad_img,imar[i])
		aprd=tf.tensor_scatter_nd_add(aprd,rad_img,imai[i])

		if single:
			bprd=tf.tensor_scatter_nd_add(zero,rad_img,imbr)
			bprd=tf.tensor_scatter_nd_add(bprd,rad_img,imbi)
		else:
			bprd=tf.tensor_scatter_nd_add(zero,rad_img,imbr[i])
			bprd=tf.tensor_scatter_nd_add(bprd,rad_img,imbi[i])

		frc.append(cross/tf.sqrt(aprd*bprd))

	if avg>len(frc[0]): avg=-1
	if avg>0:
		frc=tf.stack(frc)
		if weight!=1.0:
			w=np.linspace(weight,2.0-weight,nr)
			frc=frc*w
		return tf.math.reduce_mean(frc[:,minfreq:avg],1)
	elif avg==-1: return tf.math.reduce_mean(frc,1)
	else: return frc

FSC_REFS={}
def tf_fsc(ima,imb):
	"""Computes the FSC between a stack of complex volumes and a single reference volume. Returns a stack of 1D FSC curves."""
	if ima.dtype!=tf.complex64 or imb.dtype!=tf.complex64 : raise Exception("tf_fsc requires FFTs")


#### Project 3d Gaussian coordinates based on transforms to make projection
##   input:  pts - ( batch size, number of Gaussian, 3 (x,y,z) )
##                 ( number of Gaussian, 3) should also work
##           ang - ( batch size, 5 (az, alt, phi, tx, ty) )
#@tf.function
def xf2pts(pts, ang):

	#### input EMAN style euler angle (az, alt, phi) and make projection matrix
	##   note we need to be able to deal with a batch of particles at once
	##   so everything is in matrix form
	azp=-ang[:,0]
	altp=ang[:,1]
	phip=-ang[:,2]

	matrix=tf.stack([(tf.cos(phip)*tf.cos(azp) - tf.cos(altp)*tf.sin(azp)*tf.sin(phip)),
	(tf.cos(phip)*tf.sin(azp) + tf.cos(altp)*tf.cos(azp)*tf.sin(phip)),
	(tf.sin(altp)*tf.sin(phip)),

	(-tf.sin(phip)*tf.cos(azp) - tf.cos(altp)*tf.sin(azp)*tf.cos(phip)),
	(-tf.sin(phip)*tf.sin(azp) + tf.cos(altp)*tf.cos(azp)*tf.cos(phip)),
	(tf.sin(altp)*tf.cos(phip)),

	(tf.sin(altp)*tf.sin(azp)),
	(-tf.sin(altp)*tf.cos(azp)),
	tf.cos(altp)], 0)

	matrix=tf.transpose(matrix)
	matrix=tf.reshape(matrix, shape=[-1, 3,3]) #### Here we get a batch_size x 3 x 3 matrix

	#### rotate Gaussian positions
	##   here we try to make it also work when pts contains only the neutral model
	if len(pts.shape)>2:
		pts_rot=tf.tensordot(pts, matrix, [[2],[2]])
		pts_rot=tf.transpose(pts_rot, (0,2,1,3))

		#### the eye matrix here is mathematically unnecessary
		##   but somehow tensorflow 2.0 does not track gradient properly without it...
		##   shouldn't do much damage on the performance anyway
		e=tf.eye(pts.shape[0], dtype=bool)#.flatten()
		pts_rot=pts_rot[e]

	else:
		pts_rot=tf.tensordot(pts, matrix, [[1],[2]])
		pts_rot=tf.transpose(pts_rot, [1,0,2])

	#### finally do the translation
	tx=ang[:,3][:,None]
	ty=ang[:,4][:,None]
#	pts_rot_trans=tf.stack([(pts_rot[:,:,0]+tx), (-pts_rot[:,:,1])+ty], 2)
	pts_rot_trans=tf.stack([(-pts_rot[:,:,1])+ty,(pts_rot[:,:,0]+tx)], 2)

	#pts_rot_trans=pts_rot_trans*sz+sz/2
	return pts_rot_trans


#### make 2D projections from Gaussian coordinates in Fourier space
##   input:  pts - ( batch size, number of Gaussian, 5 (x,y,z,amp,sigma) )
##                 ( number of Gaussian, 3) should also work
##           ang - ( batch size, 5 (az, alt, phi, tx, ty) )
##        params - a dictionary of some Fourier indices for slicing
##                 sz - Fourier box size
##                 idxft - Fourier indices
##                 rrft - radial Fourier indices
##            lp - lowpass filter applied to the images
##                 this should not be necessary since we use FRC for loss
##                 but the dynamic range of values in Fourier space can sometimes be too high...
##           sym - symmetry string
#@tf.function
def pts2img(pts, ang, params, lp=.1, sym="c1"):
	bsz=ang.shape[0]
	sz, idxft, rrft=params["sz"], params["idxft"], params["rrft"]
	xfo=params["xforigin"]

	### initialize output and parse input
	imgs=tf.zeros((bsz, sz,sz), dtype=floattype)
	if len(pts.shape)>2 and pts.shape[0]>1:
		ni=pts.shape[1]
		pts=tf.reshape(pts, (-1, pts.shape[-1]))
		bamp=tf.reshape(pts[:, 3], (bsz,-1))
		multmodel=True

	else:
		bamp=pts[:, 3][None, :]
		multmodel=False

	### when a non c1 symmetry is provided, this will return a list of points
	##  one for each asymmetrical unit so we loop thro`ugh them and sum the images
	p0=get_sym_pts(sym, pts)
	for p in p0:
		p=tf.transpose(p)
		if multmodel:
			p=tf.reshape(p, (bsz, ni, -1))

		## need to change from (-0.5, 0.5) to actual image coordinates
		bpos=xf2pts(p,ang)
		bpos=bpos*sz+sz/2

		bposf=tf.floor(bpos)
		bposi=tf.cast(bposf,tf.int32)	# integer index
		bposf=bpos-bposf				# remainder used for bilinear interpolation

		# messy tensor math here to implement bilinear interpolation
		bamp0=bamp*(1.0-bposf[:,:,0])*(1.0-bposf[:,:,1])	#0,0
		bamp1=bamp*(bposf[:,:,0])*(1.0-bposf[:,:,1])	#1,0
		bamp2=bamp*(bposf[:,:,0])*(bposf[:,:,1])		#1,1
		bamp3=bamp*(1.0-bposf[:,:,0])*(bposf[:,:,1])	#0,1
		bampall=tf.concat([bamp0,bamp1,bamp2,bamp3],1)
		bposall=tf.concat([bposi,bposi+(1,0),bposi+(1,1),bposi+(0,1)],1)
		imgs=tf.stack([tf.tensor_scatter_nd_add(imgs[i],bposall[i],bampall[i]) for i in range(imgs.shape[0])])

		#try: imgs=tf.tensor_scatter_nd_add(imgs,bposi,bamp)
		#except:
			#print(imgs.shape,bposi.shape,bamp.shape)
			#raise Exception

	fimgs=tf.signal.rfft2d(imgs)

	return (tf.math.real(fimgs)*xfo,tf.math.imag(fimgs)*xfo)
