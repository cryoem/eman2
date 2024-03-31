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

At the moment this program provides only an option for estimating the gain image from a large set of counting-mode images. Eventually this will include movie alignment.
	"""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--est_gain", type=str,help="specify output file for gain image. Estimates a gain image when given a set of many movies via hierarchical median estimation", default=None)
	parser.add_argument("--align",action="store_true",default=False,help="Performs movie alignment via progressive ")
	parser.add_argument("--align_gain",type=str,help="Gain image for correcting movie images before alignment. Applied correction is to divide by the gain image.",default=None)
	parser.add_argument("--acftest",action="store_true",default=False,help="compute ACF images for input stack")
	parser.add_argument("--ccftest",action="store_true",default=False,help="compute CCF between each image and the middle image in the movie")
	parser.add_argument("--ccfdtest",action="store_true",default=False,help="compute the CCT between each image and the next image in the movie, length n-1")
	parser.add_argument("--ccftiletest",action="store_true",default=False,help="test on tiled average of CCF")
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higher number means higher level of verboseness")

	(options, args) = parser.parse_args()

	pid=E3init(argv)
	nmov=len(args)

	if options.est_gain is not None:

		sig=EMData()
		avgr=Averagers.get("mean",{"sigma":sig})
		for i in range(nmov):
			if (i%10==0): E3progress(pid,i/nmov)
			nimg=EMUtil.get_image_count(args[i])
			print(f"{args[i]}: {nimg}")
			if (nimg>500) : raise Exception("Can't deal with movies with >500 frames at present")
			for j in range(0,nimg,50):
				imgs=EMData.read_images(f"{args[i]}:{j}:{j+50}")
				for img in imgs: avgr.add_image(img)

		avg=avgr.finish()
		avg.write_image("average.hdf",0)
		sig.write_image("average.hdf",1)
		avg2=avg.copy()
		avg2.mult(sig.process("math.reciprocal"))
		avg2.write_image("average.hdf",2)

	if options.align :
		for mi in range(nmov):
			base=base_name(args[mi])
			nimg=EMUtil.get_image_count(args[mi])
			imgs=EMStack2D(EMData.read_images(f"{args[mi]}:0:{min(32,nimg)}"))
			if options.align_gain is not None:
				gain=EMData(options.align_gain,0)
				gain.process_inplace("math.reciprocal")
				for i in imgs.emdata: i.mult(gain)
			for i in imgs.emdata: i.process_inplace("normalize.edgemean")
			imgs_clip=imgs.center_clip(3072)	# this should fit the dimensions of pretty much any currently used movie mode sensor
			ffts=imgs_clip.do_fft()
			imgs_clip=None

			# compute the ACF of each image in the movie, then average them
			acfs=ffts.calc_ccf(ffts,offset=0)	# ACF stack
			acfsr=acfs.do_ift()					# real space version
			_,nx,ny=acfsr.shape
			acfsrc=acfsr.tensor[:,nx//2-64:nx//2+64,ny//2-64:ny//2+64]	# pull out the central region limiting shifts to +-64 pixels
			acfav=EMStack2D(tf.math.reduce_mean(acfsrc,0,keepdims=True))	# average acf over all images
			a=acfav.numpy[0]		# this is to permit manipulation, since TF object is a constant

			# we zero out the region near the origin
			a[63,63]=(a[62,63]+a[63,62])/2.0
			a[64,63]=a[64,62]
			a[65,63]=(a[65,62]+a[66,63])/2.0
			a[63,64]=a[62,64]
			a[65,64]=a[66,64]
			a[63,65]=(a[62,65]+a[63,66])/2.0
			a[64,65]=a[64,66]
			a[65,65]=(a[66,65]+a[65,66])/2.0
			a[64,64]=(a[64,63]+a[63,64]+a[65,64]+a[64,65])/4.0
			print(np.std(acfav),np.mean(acfav))
			a/=np.std(a)
			print(np.std(acfav),np.mean(acfav))

			acfavf=acfav.do_fft()	# this is the ACF FFT we will align the CCFS to
			#acfav.write_images("acfref.hdf")

			# CCF between each image and the subsequent image
			ccfs=ffts.calc_ccf(ffts,offset=1)
			ccfsr=ccfs.do_ift()
			ccfsr=EMStack2D(ccfsr.tensor[:,nx//2-64:nx//2+64,ny//2-64:ny//2+64])
			ccfsr.numpy[:,64,64]=(ccfsr[:,63,64]+ccfsr[:,65,64]+ccfsr[:,64,63]+ccfsr[:,64,65])/4.0
			for i in ccfsr.numpy: i/=np.std(i)
			#ccfsr.write_images("ccfsr.hdf")
			ccfsr.coerce_tensor()
			ccfsf=ccfsr.do_fft()

			# CCF between CCFS and ACF ref
			ccfacf=ccfsf.calc_ccf(acfavf)
			ccfacfr=ccfacf.do_ift()

			# Extract the peak locations
			peaks=tf.math.argmax(tf.reshape(ccfacfr.tensor,[len(ccfacfr),128*128]),axis=1)
			seqoff=np.insert((tf.unravel_index(peaks,dims=[128,128])-64).numpy(),0,(0,0),1)		# insert adds a 0 shift element for the first image
			# Convert the n -> n+1 shifts into absolute shifts
			for i in range(1,len(seqoff[0])):
				seqoff[0][i]+=seqoff[0][i-1]
				seqoff[1][i]+=seqoff[1][i-1]

			# recenter the result to minimize edge effects
			seqoff[0]-=seqoff[0].mean().astype("int64")
			seqoff[1]-=seqoff[1].mean().astype("int64")

			print(seqoff)

			avgorig=Averagers.get("mean")
			avgali=Averagers.get("mean")
			for i in range(min(32,nimg)):
				im=EMData(f"{args[mi]}",i)
				im.mult(gain)
				avgorig.add_image(im)
				im.translate(int(seqoff[0][i]),int(seqoff[1][i]),0)
				avgali.add_image(im)

			avgorig.finish().write_image(f"micrographs/{base}_unali.hdf:6")
			avgali.finish().write_image(f"micrographs/{base}.hdf:6")
	#		for i in ccfacfr.tensor: print(tf.math.argmax(tf.reshape(i,[128*128])))
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

	if options.ccfdtest:
		avg=EMData("average.hdf",0)
		avg.div(avg["mean"])
		#avg.add(-avg["mean"])
		nimg=EMUtil.get_image_count(args[0])
		imgs=EMStack2D(EMData.read_images(f"{args[0]}:0:{min(50,nimg)}"))
		for im in imgs:
			im.div(avg)
		ffts=imgs.do_fft()
		ccfs=ffts.calc_ccf(ffts,offset=1)
		ccfsr=ccfs.do_ift()
		_,nx,ny=ccfsr.shape

		cens=EMStack2D(ccfsr.tensor[:,nx//2-64:nx//2+64,ny//2-64:ny//2+64])
		for im in cens.emdata: im.process_inplace("normalize.edgemean")
		cens.write_images("ccfs.hdf")

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
