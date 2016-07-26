#!/usr/bin/env python

from EMAN2 import *
import numpy as np
from scipy import ndimage

def main():
progname = os.path.basename(sys.argv[0])
usage = """prog [options] stack1.hdf stack2.mrcs ...

Program to erase gold from DDD movies. Requires scipy.
"""

parser = EMArgumentParser(usage=usage,version=EMANVERSION)

parser.add_argument("--writenoise", default=False, action="store_true", help="Save noise image(s).")
parser.add_argument("--average", default=False, action="store_true", help="Erase gold from average of input stack(s).")
parser.add_argument("--downsamp", default=1, type=int, help="Downsample the input stack(s). Default is 1, i.e. no downsampling.")
parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-2)
(options, args) = parser.parse_args()

for arg in args:
	if options.verbose: print(arg)
	frames = load(arg,ds=options.downsample,inv=True)
	if options.average:
		avgr = Averagers.get("mean")
		avgr.add_image_list(frames)
		img = avgr.finish()
		img.process_inplace("normalize")
		sharp_msk, soft_msk = generate_masks(img)
		mskd_sharp = sharp_msk*img
		sub_sharp = img-mskd_sharp
		noise = local_noise(sub_sharp)
		if options.writenoise:
			noise.write_image("{}_noise.hdf")
		mskd_soft = soft_msk*img
		sub_soft = img-mskd_soft
		result = sub_soft + noise * soft_msk
		result *= -1
		result.write_image("{}_proc.hdf")
	else:
		for i,f in enumerate(frames):
			f.process_inplace("normalize")
			sharp_msk, soft_msk = generate_masks(img)
			mskd_sharp = sharp_msk*img
			sub_sharp = img-mskd_sharp
			noise = local_noise(sub_sharp)
			if options.writenoise:
				noise.write_image("{}_noise.hdf",i)
			mskd_soft = soft_msk*img
			sub_soft = img-mskd_soft
			result = sub_soft + noise * soft_msk
			result *= -1
			result.write_image("{}_proc.hdf",i)

def load(fn,inv=False,ds=1):
	frames = []
	for i in range(EMUtil.get_image_count(fn)):
		f = EMData(fn,i)
		if ds != 1: f.process_inplace("math.fft.resample",{"n":ds})
		if inv: f*=-1
		frames.append(f)
	return frames

def generate_masks(img,goldsize=30):
	img.process_inplace("normalize")
	# create sharp mask
	msk = img.process("filter.highpass.gauss",{"cutoff_pixels":25})
	msk.process_inplace("filter.lowpass.gauss",{"cutoff_abs":0.1})
	msk.process_inplace("threshold.clampminmax",{"maxval":msk["maximum"],"minval":msk["mean"]+3*msk["sigma"],"tozero":True})
	msk.process_inplace("threshold.binary",{"value":msk["mean"]})
	# remove dust
	nproc = msk.numpy().copy()
	s = np.sqrt(goldsize*2).astype(int)
	se2=np.ones((s,s))
	nproc = ndimage.binary_closing(nproc,structure=se2).astype(int)
	nproc = ndimage.binary_opening(nproc,structure=se2).astype(int)
	sharp_msk = from_numpy(nproc)
	# grow slightly and create soft mask
	sharp_msk = sharp_msk.process("mask.addshells.gauss",{"val1":8,"val2":0})
	soft_msk = sharp_msk.process("mask.addshells.gauss",{"val1":0,"val2":8})
	return sharp_msk,soft_msk

def local_noise(img,bs=128,oversamp=6):
	localnoise = EMData(img["nx"],img["ny"])
	localnoise.to_zero()
	nbxs = len(np.arange(-bs,img['nx']+bs,bs))*oversamp
	nbys = len(np.arange(-bs,img['ny']+bs,bs))*oversamp
	mx = np.linspace(0,img["nx"],nbxs).astype(int)
	my = np.linspace(0,img["ny"],nbys).astype(int)
	for x in mx:
		for y in my:
			r = img.get_clip(Region(x-bs/2,y-bs/2,bs,bs))
			n = EMData(bs,bs)
			n.to_zero()
			n.process_inplace("math.addnoise",{"noise":r["sigma_nonzero"]})
			n.process_inplace("filter.highpass.gauss",{"cutoff_abs":0.01})
			n *= r["sigma_nonzero"]/n["sigma_nonzero"]
			n += r["mean_nonzero"]
			localnoise.insert_clip(n,(x-bs/2,y-bs/2))
	return localnoise


if __name__ == "__main__":
	main()
