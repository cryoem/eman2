#!/usr/bin/env python

#====================
#Author: Michael Bell July, 2016 (minor edits, Jesus Galaz-Montoya). Last update: July, 2016
#====================
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
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  2111-1307 USA



from EMAN2 import *
import numpy as np
from scipy import ndimage

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options] stack1.hdf stack2.mrcs ...

	Program to erase gold fiducials and other high-density features from images, such as frames in DDD movies or images in tiltseries. Requires scipy.
	"""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_argument("--writenoise", default=False, action="store_true", help="Save noise image(s).")
	parser.add_argument("--average", default=False, action="store_true", help="Erase gold from average of input stack(s).")
	parser.add_argument("--downsample", default=1, type=int, help="Downsample the input stack(s). Default is 1, i.e. no downsampling.")
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-2)
	
	(options, args) = parser.parse_args()

	nfiles = len(args)
	
	kk=0
	for arg in args:
		newarg=''
		originalarg = arg
		if '.ali' == arg[-4:] or '.mrc' == arg[-4:]:
			newarg = arg.replace(arg[-4:],'.hdf') #intermediate stack in converting .ali or .mrc images into a stack of 2-D images that can be processed by this program. It will be deleted in the end.
			cmd = 'e2proc2d.py ' + arg + ' dummy_stack.hdf --threed2twod'
			
			#Unfortunately, e2proc2d.py appends to existing files instead of overwriting them. If you run this two consecutive times, you'll find your stack growing.
			#To prevent this, we write the output to a 'dummy' file, which gets cleared at the beginning. If the program runs successfully to the end, the dummy file gets renamed.
			try:
				os.remove('dummy_stack.hdf')
			except:
				pass

			runcmd(options,cmd)
			os.rename('dummy_stack.hdf',newarg)
			arg = newarg

		if options.verbose: 
			#print(arg)
			print "\nprocessing file %s with %d images in it" %(arg, EMUtil.get_image_count(arg))
		
		frames = load(arg,ds=options.downsample,inv=True)
		
		if options.average:
			avgr = Averagers.get("mean")
			avgr.add_image_list(frames)
			img = avgr.finish()
			img.process_inplace("normalize")
			sharp_msk, soft_msk = generate_masks(img)
			mskd_sharp = sharp_msk*img
			sub_sharp = img-mskd_sharp
			noise = local_noise(options,sub_sharp)
			if options.writenoise:
				noise.write_image("{}_noise.hdf")
			mskd_soft = soft_msk*img
			sub_soft = img-mskd_soft
			result = sub_soft + noise * soft_msk
			result *= -1
			result.write_image("{}_proc.hdf")
		else:
			for i,f in enumerate(frames):
				#if i > 0:
				#	break
				print "processing image", i
				f.process_inplace("normalize")
				sharp_msk, soft_msk = generate_masks(f)
				mskd_sharp = sharp_msk*f
				sub_sharp = f-mskd_sharp

				
				noise = local_noise(options,sub_sharp)
				
				apix = f['apix_x']
				noise['apix_x'] = apix
				noise['apix_y'] = apix
				#noise['apix_z'] = apix
				nyquistres=apix*2.0
				filtres=nyquistres*10.0/9.0
				filtfreq=1.0/filtres
				print "apix is",apix
				print "resolution and corresponding frequency to filter at",filtres,filtfreq
				noise.process_inplace('filter.lowpass.tanh', {'cutoff_freq':filtfreq,'apix':apix})
				outputnoise = arg.replace('.hdf','_noise.hdf')

				if options.writenoise:
					noise.write_image(outputnoise,i)
				mskd_soft = soft_msk*f
				sub_soft = f-mskd_soft
				sub_soft.write_image( arg.replace('.hdf','_masked.hdf'),i )
				
				result = sub_soft + noise * soft_msk
				result *= -1
				output = arg.replace('.hdf','_proc.hdf')
				result.write_image(output,i)
		
		if newarg:
			os.remove(newarg)
		
		if '.ali' == originalarg[-4:] or '.mrc' == originalarg[-4:]:
			#intermediate = arg.replace('.hdf','.mrcs')
			finaloutput = arg.replace('.hdf',originalarg[-4:])
			cmd = 'e2proc2d.py ' + arg + ' ' + finaloutput + ' --twod2threed --outmode int16'
			runcmd(options,cmd)
			os.remove(arg)

		kk+=1


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


def local_noise(options,img,bs=128,oversamp=6):
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
			try:
				n *= r["sigma_nonzero"]/n["sigma_nonzero"]
			except:
				if options.verbose > 9:
					print "WARNING: division by zero, from n['sigma_nonzero']", n["sigma_nonzero"]
			n += r["mean_nonzero"]
			localnoise.insert_clip(n,(x-bs/2,y-bs/2))
	return localnoise


def runcmd(options,cmd):
	if options.verbose > 8:
		print "(erase_gold)(runcmd) running command", cmd
	
	p=subprocess.Popen( cmd, shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	text=p.communicate()	
	p.stdout.close()
	
	if options.verbose > 8:
		print "(erase_gold)(runcmd) done"

if __name__ == "__main__":
	main()
