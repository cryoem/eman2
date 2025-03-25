#!/usr/bin/env python
#
# Author: Steven Ludtke  06/07/2023
# Copyright (c) 2023- Baylor College of Medicine
#
# This software is issued under a joint BSD/GNU license. You may use the
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

from EMAN3 import *
from EMAN3tensor import *
import random
from sys import argv

def main():

	usage="""e3movie.py <movie stack> ...

This will estimate gain for counting mode cameras, and has a "first draft" of frame alignment.
	"""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--est_gain", type=str,help="specify output file for gain image. Estimates a gain image when given a set of many movies via hierarchical median estimation", default=None)
	parser.add_argument("--seqavg",type=int,default=-1,help="Average N frames in sequence before alignment")
	parser.add_argument("--shrinkout",type=float,default=-1,help="Fourier downsample the output movie stack by a factor of N (may be fractional)")
	parser.add_argument("--alignbyccf",action="store_true",default=False,help="Performs movie alignment via sequential ccfs ")
	parser.add_argument("--alignbyacfccf",action="store_true",default=False,help="Performs movie alignment via progressive ")
	parser.add_argument("--align_gain",type=str,help="Gain image for correcting movie images before alignment. Applied correction is to divide by the gain image.",default=None)
	parser.add_argument("--clip",type=str,default=None,help="nx,ny output image size. Trims both edges equally as necessary")
	parser.add_argument("--frames",type=str,default=None,help="<first>,<last+1> movie frames to use, first frame is 0, '0,3' will use frames 0,1,2")
	parser.add_argument("--acftest",action="store_true",default=False,help="compute ACF images for input stack")
	parser.add_argument("--ccftest",action="store_true",default=False,help="compute CCF between each image and the middle image in the movie")
	parser.add_argument("--ccfdtest",type=str,default=None,help="compute the CCF between each image and the next image in the movie, length n-1, provide the filename of the gain correction image")
	parser.add_argument("--ccftiletest",action="store_true",default=False,help="test on tiled average of CCF")
	parser.add_argument("--gpudev",type=int,help="GPU Device, default 0", default=0)
	parser.add_argument("--gpuram",type=int,help="Maximum GPU ram to allocate in MB, default=4096", default=4096)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higher number means higher level of verboseness")

	(options, args) = parser.parse_args()
	tf_set_device(dev=0,maxmem=options.gpuram)

	pid=E3init(argv)
	nmov=len(args)
	if options.clip is not None: clip=(int(options.clip.split(",")[0]),int(options.clip.split(",")[1]))
	else: clip=None
	if options.frames is not None: options.frames=(int(options.frames.split(",")[0]),int(options.frames.split(",")[1]))
	else: options.frames=(0,100000)

	if options.est_gain is not None:

		sig=EMData()
		avgr=Averagers.get("mean",{"sigma":sig})
		for i in range(nmov):
			if (i%10==0): E3progress(pid,i/nmov)
			nimg=EMUtil.get_image_count(args[i])
			print(f"{args[i]}: {nimg}")
			#if (nimg>500) : raise Exception("Can't deal with movies with >500 frames at present")
			for j in range(0,nimg,50):
				imgs=EMData.read_images(f"{args[i]}:{j}:{j+50}")
				for img in imgs: avgr.add_image(img)

		avg=avgr.finish()
		avg.mult(1.0/avg["mean"])	# this insures that on average e- counts are still counts
		sig.mult(1.0/avg["mean"])	# scale the same as the average
		sigm=sig.process("threshold.binaryrange",{"low":max(0.1,sig["mean"]-sig["sigma"]*4.0),"high":sig["mean"]+sig["sigma"]*4.0}) 	#pixels with small or large standard deviations are likely bad
		avg.mult(sigm)
		avg.write_image(options.est_gain,0)
		sig.write_image(options.est_gain,1)

	if options.alignbyccf :
		rsz=128				# 1/2 the size of the CCF regions to correlate
		gain=EMData(options.align_gain,0)
		for mi in range(nmov):
			base=base_name(args[mi])[0]
			nimg=EMUtil.get_image_count(args[mi])
			frames=(options.frames[0],min(options.frames[1],nimg))
			if options.seqavg>1:
				framesteps=(frames[1]-frames[0])//options.seqavg
				imgs=EMStack2D([sum(EMData.read_images(f"{args[mi]}:{frames[0]+i*options.seqavg}:{frames[0]+(i+1)*options.seqavg}")) for i in range(framesteps)])
			else: imgs=EMStack2D(EMData.read_images(f"{args[mi]}:{frames[0]}:{frames[1]}"))
			if options.align_gain is not None:
				for i in imgs.emdata: i.process_inplace("math.fixgain.counting",{"gain":gain,"gainmin":2,"gainmax":2})
			# normalization only for alignment
			for i in imgs.emdata: i.process_inplace("normalize.edgemean")
			imgs_clip=imgs.center_clip(3072)	# this should fit the dimensions of pretty much any currently used movie mode sensor, and still use enough of the image
			ffts=imgs_clip.do_fft()
			imgs_clip=None

			# CCF between each image and the subsequent image
			ccfs=ffts.calc_ccf(ffts,offset=1)
			ccfsr=ccfs.do_ift()
			ccfsr=ccfsr.center_clip(rsz*2)
			#ccfsr=EMStack2D(ccfsr.tensor[:,nx//2-rsz:nx//2+rsz,ny//2-rsz:ny//2+rsz])
			ccfsr.numpy[:,rsz,rsz]=(ccfsr[:,rsz-1,rsz]+ccfsr[:,rsz+1,rsz]+ccfsr[:,rsz,rsz-1]+ccfsr[:,rsz,rsz+1])/4.0
			for i in ccfsr.numpy:
				i-=np.mean(i)
				i/=np.std(i)
#			ccfsr.write_images("ccfsr.hdf")
			ccfsr.coerce_tensor()
			ccfsf=ccfsr.do_fft()

			# make a gaussian low pass filter to smooth the ccfs after removing the 0 translation value
			filt=EMData(rsz*2,rsz*2,1)
			filt.to_one()
			filt.process_inplace("mask.gaussian",{"outer_radius":2})
			filt.process_inplace("xform.phaseorigin.tocorner")
			filttf=tf_fft2d(to_tf(filt))

			ccfsfilt=ccfsf.convolve(filttf)
			ccfsreal=ccfsfilt.do_ift()

			# Extract the peak locations
			peaks=tf.math.argmax(tf.reshape(ccfsreal.tensor,[len(ccfsreal),rsz*2*rsz*2]),axis=1)
			seqoff=np.insert((tf.unravel_index(peaks,dims=[rsz*2,rsz*2])-rsz).numpy(),0,(0,0),1)		# insert adds a 0 shift element for the first image
			# Convert the n -> n+1 shifts into absolute shifts
			print(seqoff)
			for i in range(1,len(seqoff[0])):
				seqoff[0][i]+=seqoff[0][i-1]
				seqoff[1][i]+=seqoff[1][i-1]

			# recenter the result to minimize edge effects
			seqoff[0]-=seqoff[0].mean().astype("int64")
			seqoff[1]-=seqoff[1].mean().astype("int64")

			print(seqoff)
#			ccfsreal.write_images("ccfs.hdf")

			avgorig=Averagers.get("mean")
			avgali=Averagers.get("mean")
			for i in range(frames[0],frames[1],max(options.seqavg,1)):
				if options.seqavg>1: im=sum(EMData.read_images(f"{args[mi]}:{i}:{i+options.seqavg}"))
				else: im=EMData(f"{args[mi]}",i)
				if options.align_gain is not None:
					im.process_inplace("math.fixgain.counting",{"gain":gain,"gainmin":2,"gainmax":2})
				if clip is not None: im=im.get_clip(Region((im["nx"]-clip[0])//2,(im["ny"]-clip[1])//2,clip[0],clip[1]))
				avgorig.add_image(im)
				im.translate(int(seqoff[1][i-frames[0]]),int(seqoff[0][i-frames[0]]),0)		# note the apparent x/y swap here due to the numpy/tf N/Y/X indices
				avgali.add_image(im)

			avgo=avgorig.finish()
			avga=avgali.finish()
			if options.shrinkout>1.0: 
				avgo.process_inplace("math.fft.resample",{"n":options.shrinkout})
				avga.process_inplace("math.fft.resample",{"n":options.shrinkout})
			avgo.write_image(f"micrographs/{base}_{frames[0]}-{frames[1]}_unali.hdf:6")
			avga.write_image(f"micrographs/{base}_{frames[0]}-{frames[1]}.hdf:6")
			js=js_open_dict(info_name(args[mi]))
			js["movie_frames"]=frames
			js["movie_align_x"]=seqoff[1]
			js["movie_align_y"]=seqoff[0]
			js=None
	#		for i in ccfacfr.tensor: print(tf.math.argmax(tf.reshape(i,[rsz*2*rsz*2])))
	#		for im in ccfacfr.numpy: im/=np.std(im)
	#		ccfacfr.write_images("ccfacfs.hdf")

	if options.alignbyacfccf :
		rsz=128				# 1/2 the size of the CCF regions to correlate
		gain=EMData(options.align_gain,0)
		for mi in range(nmov):
			base=base_name(args[mi])[0]
			nimg=EMUtil.get_image_count(args[mi])
			frames=(options.frames[0],min(options.frames[1],nimg))
			imgs=EMStack2D(EMData.read_images(f"{args[mi]}:{frames[0]}:{frames[1]}"))
			if options.align_gain is not None:
				for i in imgs.emdata: i.process_inplace("math.fixgain.counting",{"gain":gain,"gainmin":2,"gainmax":2})
			for i in imgs.emdata: i.process_inplace("normalize.edgemean")
			imgs_clip=imgs.center_clip(3072)	# this should fit the dimensions of pretty much any currently used movie mode sensor, and still use enough of the image
			ffts=imgs_clip.do_fft()
			imgs_clip=None

			# compute the ACF of each image in the movie, then average them
			acfs=ffts.calc_ccf(ffts,offset=0)	# ACF stack
			acfsr=acfs.do_ift()					# real space version
			_,nx,ny=acfsr.shape
			acfsrc=acfsr.tensor[:,nx//2-rsz:nx//2+rsz,ny//2-rsz:ny//2+rsz]	# pull out the central region limiting shifts to +-rsz pixels
			acfav=EMStack2D(tf.math.reduce_mean(acfsrc,0,keepdims=True))	# average acf over all images
			a=acfav.numpy[0]		# this is to permit manipulation, since TF object is a constant

			# we zero out the region near the origin
			a[rsz-1,rsz-1]=(a[rsz-2,rsz-1]+a[rsz-1,rsz-2])/2.0
			a[rsz,rsz-1]=a[rsz,rsz-2]
			a[rsz+1,rsz-1]=(a[rsz+1,rsz-2]+a[rsz+2,rsz-1])/2.0
			a[rsz-1,rsz]=a[rsz-2,rsz]
			a[rsz+1,rsz]=a[rsz+2,rsz]
			a[rsz-1,rsz+1]=(a[rsz-2,rsz+1]+a[rsz-1,rsz+2])/2.0
			a[rsz,rsz+1]=a[rsz,rsz+2]
			a[rsz+1,rsz+1]=(a[rsz+2,rsz+1]+a[rsz+1,rsz+2])/2.0
			a[rsz,rsz]=(a[rsz,rsz-1]+a[rsz-1,rsz]+a[rsz+1,rsz]+a[rsz,rsz+1])/4.0
			print(np.std(acfav),np.mean(acfav))
			a-=np.mean(a)
			a/=np.std(a)
#			print(np.std(acfav),np.mean(acfav))

			acfavf=acfav.do_fft()	# this is the ACF FFT we will align the CCFS to
#			acfav.write_images("acfref.hdf")

			# CCF between each image and the subsequent image
			ccfs=ffts.calc_ccf(ffts,offset=1)
			ccfsr=ccfs.do_ift()
			ccfsr=EMStack2D(ccfsr.tensor[:,nx//2-rsz:nx//2+rsz,ny//2-rsz:ny//2+rsz])
			ccfsr.numpy[:,rsz,rsz]=(ccfsr[:,rsz-1,rsz]+ccfsr[:,rsz+1,rsz]+ccfsr[:,rsz,rsz-1]+ccfsr[:,rsz,rsz+1])/4.0
			for i in ccfsr.numpy:
				i-=np.mean(i)
				i/=np.std(i)
#			ccfsr.write_images("ccfsr.hdf")
			ccfsr.coerce_tensor()
			ccfsf=ccfsr.do_fft()

			# CCF between CCFS and ACF ref
			ccfacf=ccfsf.calc_ccf(acfavf)
			ccfacfr=ccfacf.do_ift()
#			ccfacfr.write_images("ccfacfr.hdf")

			# Extract the peak locations
			peaks=tf.math.argmax(tf.reshape(ccfacfr.tensor,[len(ccfacfr),rsz*2*rsz*2]),axis=1)
			seqoff=np.insert((tf.unravel_index(peaks,dims=[rsz*2,rsz*2])-rsz).numpy(),0,(0,0),1)		# insert adds a 0 shift element for the first image
			# Convert the n -> n+1 shifts into absolute shifts
			print(seqoff)
			for i in range(1,len(seqoff[0])):
				seqoff[0][i]+=seqoff[0][i-1]
				seqoff[1][i]+=seqoff[1][i-1]

			# recenter the result to minimize edge effects
			seqoff[0]-=seqoff[0].mean().astype("int64")
			seqoff[1]-=seqoff[1].mean().astype("int64")

			print(seqoff)

			avgorig=Averagers.get("mean")
			avgali=Averagers.get("mean")
			for i in range(frames[0],frames[1]):
				im=EMData(f"{args[mi]}",i)
				if options.align_gain is not None:
					im.process_inplace("math.fixgain.counting",{"gain":gain,"gainmin":2,"gainmax":2})
				if clip is not None: im=im.get_clip(Region((im["nx"]-clip[0])//2,(im["ny"]-clip[1])//2,clip[0],clip[1]))
				avgorig.add_image(im)
				im.translate(int(seqoff[1][i-frames[0]]),int(seqoff[0][i-frames[0]]),0)		# note the apparent x/y swap here due to the numpy/tf N/Y/X indices
				avgali.add_image(im)

			avgorig.finish().write_image(f"micrographs/{base}_{frames[0]}-{frames[1]}_unali.hdf:6")
			avgali.finish().write_image(f"micrographs/{base}_{frames[0]}-{frames[1]}.hdf:6")
			js=js_open_dict(info_name(args[mi]))
			js["movie_frames"]=frames
			js["movie_align_x"]=seqoff[1]
			js["movie_align_y"]=seqoff[0]
			js=None
	#		for i in ccfacfr.tensor: print(tf.math.argmax(tf.reshape(i,[rsz*2*rsz*2])))
	#		for im in ccfacfr.numpy: im/=np.std(im)
	#		ccfacfr.write_images("ccfacfs.hdf")


	if options.ccftest:
		avg=EMData("average.hdf",0)
		avg.div(avg["mean"])
		#avg.add(-avg["mean"])
		nimg=EMUtil.get_image_count(args[0])
		imgs=EMStack2D(EMData.read_images(f"{args[0]}:0:{min(50,nimg)}"))
		for im in imgs:
			im.div(avg)
		ffts=imgs.do_fft()
		ccfs=ffts.calc_ccf(ffts[len(imgs)//2])
		ccfsr=ccfs.do_ift()
		_,nx,ny=ccfsr.shape

		cens=EMStack2D(ccfsr.tensor[:,nx//2-64:nx//2+64,ny//2-64:ny//2+64])
		cens.write_images("ccfs.hdf")

	if options.acftest:
		avg=EMData("average.hdf",0)
		avg.div(avg["mean"])
		#avg.add(-avg["mean"])
		nimg=EMUtil.get_image_count(args[0])
		imgs=EMStack2D(EMData.read_images(f"{args[0]}:0:{min(50,nimg)}"))
		for im in imgs:
			im.div(avg)
		ffts=imgs.do_fft()
		ccfs=ffts.calc_ccf(ffts,offset=0)
		ccfsr=ccfs.do_ift()
		_,nx,ny=ccfsr.shape

		cens=EMStack2D(ccfsr.tensor[:,nx//2-64:nx//2+64,ny//2-64:ny//2+64])
		for im in cens.emdata: im.process_inplace("normalize.edgemean")
		cens.write_images("ccfs.hdf")

	if options.ccfdtest is not None:
		try: os.unlink("ccfs.hdf")
		except: pass
		try: os.unlink("ccfs2k.hdf")
		except: pass
		try: os.unlink("ccfs1k.hdf")
		except: pass
		try: os.unlink("ccfeo5.hdf")
		except: pass
		try: os.unlink("ccfeo5_imgs.hdf")
		except: pass

		avg=EMData(options.ccfdtest,0)
		avg.div(avg["mean"])
		#avg.add(-avg["mean"])
		nimg=file_image_count(args[0])

		print(f"0:{min(50,nimg)}")
		aimgs=EMStack2D(EMData.read_images(f"{args[0]}:0:{min(50,nimg)}"))
		for im in aimgs:
			im.process_inplace("math.fixgain.counting",{"gain":avg,"gainmin":3,"gainmax":3})
		imgs=EMStack2D(np.stack([sum(aimgs.numpy[i:i+10:2]) for i in range(40)]))
		imgs.coerce_tensor()
		print(imgs.shape)
		ffts=imgs.do_fft()
		ccfs=ffts.calc_ccf(ffts,offset=1)
		ccfsr=ccfs.do_ift()
		_,nx,ny=ccfsr.shape
		cens=ccfsr.center_clip(64)
		for im in cens.emdata: im.process_inplace("normalize.edgemean")
		cens.write_images("ccfeo5.hdf",bits=0)
		imgs.write_images("ccfeo5_imgs.hdf",bits=5)

		for i in range(1,nimg,25):
			print(f"{i-1}:{min(i+25,nimg)}")
			imgs=EMStack2D(EMData.read_images(f"{args[0]}:{i-1}:{min(i+25,nimg)}"))
			for im in imgs:
				im.process_inplace("math.fixgain.counting",{"gain":avg,"gainmin":3,"gainmax":3})
			imgs.coerce_tensor()
#		imgs=imgs.downsample(4096)
			ffts=imgs.do_fft()
			ccfs=ffts.calc_ccf(ffts,offset=1)
			ccfsr=ccfs.do_ift()
			_,nx,ny=ccfsr.shape
			cens=ccfsr.center_clip(64)
			#cens=EMStack2D(ccfsr.tensor[:,nx//2-64:nx//2+64,ny//2-64:ny//2+64])
			for im in cens.emdata: im.process_inplace("normalize.edgemean")
			cens.write_images("ccfs.hdf",bits=0,n_start=i-1)

			imgs=imgs.center_clip(2048)
			ffts=imgs.do_fft()
			ccfs=ffts.calc_ccf(ffts,offset=1)
			ccfsr=ccfs.do_ift()
			_,nx,ny=ccfsr.shape
			cens=ccfsr.center_clip(64)
			for im in cens.emdata: im.process_inplace("normalize.edgemean")
			cens.write_images("ccfs2k.hdf",bits=0,n_start=i-1)

			imgs=imgs.center_clip(1024)
			ffts=imgs.do_fft()
			ccfs=ffts.calc_ccf(ffts,offset=1)
			ccfsr=ccfs.do_ift()
			_,nx,ny=ccfsr.shape
			cens=ccfsr.center_clip(64)
			for im in cens.emdata: im.process_inplace("normalize.edgemean")
			cens.write_images("ccfs1k.hdf",bits=0,n_start=i-1)


	if options.ccftiletest:
		avg=EMData("average.hdf",0)
		avg.div(avg["mean"])
		#avg.add(-avg["mean"])
		nimg=EMUtil.get_image_count(args[0])
		imgs=EMStack2D(EMData.read_images(f"{args[0]}:0:{min(50,nimg)}"))
		_,nx,ny=imgs.shape
		for im in imgs: im.div(avg)

		for x in range(0,nx-1024,1024):
			for y in range(0,ny-1024,1024):
				imgt=EMStack2D(imgs.tensor[:,x:x+1024,y:y+1024])

				ffts=imgt.do_fft()
				ccfs=ffts.calc_ccf(ffts,offset=1)
				ccfsr=ccfs.do_ift()

				try: ccfsum=ccfsum+ccfsr
				except: ccfsum=ccfsr

#		cens=EMStack2D(ccfsr.tensor[:,nx//2-64:nx//2+64,ny//2-64:ny//2+64])
		for im in ccfsum.emdata: im.process_inplace("normalize.edgemean")
		ccfsum.write_images("ccfs.hdf")


	E3end(argv)

if __name__ == '__main__':
	main()
