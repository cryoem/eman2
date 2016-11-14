#!/usr/bin/env python

#====================
#Author: Michael Bell July, 2016 (edits, Jesus Galaz-Montoya). Last update: September, 2016
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
try: from scipy import ndimage
except: ndimage=None
import subprocess
import os
import time
from EMAN2jsondb import JSTask,jsonclasses

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options] stack1.hdf stack2.mrcs ...

	Program to erase gold fiducials and other high-density features from images, such as frames in DDD movies or images in tiltseries. Requires scipy.
	"""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	#parser.add_argument("--average", default=False, action="store_true", help="Erase gold from average of input stack(s).")
	parser.add_argument("--apix", default=None, type=float, help="Override Apix in image header.")
	parser.add_argument("--lowpass", default=False, action="store_true", help="Also lowpass filter noise based on local properties. Useful for processing tomographic tilt series.")
	parser.add_argument("--keepdust", default=False, action="store_true", help="Do not remove 'dust' from mask (include objects smaller than gold fiducials).")
	parser.add_argument("--goldsize", default=30, type=float, help="Diameter (in pixels) of gold fiducials to erase.")
	#parser.add_argument("--downsample", default=1.0, type=float, help="Downsample the input stack(s). Default is 1, i.e. no downsampling.")
	parser.add_argument("--oversample", default=4, type=int, help="Oversample noise image to smooth transitions from regions with different noise.")
	parser.add_argument("--boxsize", default=128, type=int, help="Box size to use when computing local noise.")
	parser.add_argument("--debug", default=False, action="store_true", help="Save noise and mask/masked image(s).")
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-2)
	parser.add_argument("--parallel",type=str, default=None, help="""Default=None (not used). Parallelism. See http://blake.bcm.edu/emanwiki/EMAN2/Parallel""")

	parser.add_argument("--subset", default=0, type=int, help="Default=0 (not used). Apply algorithm to only a subset of images in each stack file.")
	parser.add_argument("--nsigmas", default=3.0,type=float, help="Default=3.0. Number of standard deviations above the mean to determine pixels to mask out (erase).")


	(options, args) = parser.parse_args()

	nfiles = len(args)

	if options.parallel:
		from EMAN2PAR import EMTaskCustomer
		etc=EMTaskCustomer(options.parallel)

	for argnum,arg in enumerate(args):

		t0 = time.time()

		newarg=''
		originalarg = arg

		hdr = EMData(arg,0,True) #load header only to get parameters used below
		if options.apix: apix = options.apix
		else: apix = hdr['apix_x']
		nx=hdr['nx']
		ny=hdr['ny']

		if '.ali' == arg[-4:] or '.mrc' == arg[-4:]:

			#Unfortunately, e2proc2d.py appends to existing files instead of overwriting them. If you run this program two consecutive times and the first one failed for whatever reason,
			#you'll find your stack growing.
			#To prevent this, we create a 'dummy' file, but first remove any dummy files from previous failed runs. (If the program runs successfully to the end, the dummy file gets renamed).
			try: os.remove('dummy_stack.hdf')
			except: pass

			#turn .ali or .mrc 3D images into a stack of 2D images that can be processed by this program.
			cmd = 'e2proc2d.py ' + arg + ' dummy_stack.hdf --threed2twod'
			if options.subset:
				cmd += ' --first 0 --last ' + str(options.subset-1)

			runcmd(options,cmd)

			#make the new stack of 2D images (dummy_stack.hdf) the new input (the name of the input file but with .hdf format); this intermediate file will be deleted in the end.
			newarg = arg.replace(arg[-4:],'.hdf')
			os.rename('dummy_stack.hdf',newarg)
			arg = newarg

		outf = "{}_proc.hdf".format( os.path.splitext(arg)[0] )
		if os.path.isfile(outf):
			print("Results are already stored in {}. Please erase or move and try again.".format(outf))
			sys.exit(1)

		nfs = EMUtil.get_image_count(arg)

		tasks=[]
		results=[]
		results=None

		#parallelized tasks don't run "in order"; therefore, a dummy stack needs to be pre-created with as many images as the final stack will have
		#(otherwise, writing output images to stack indexes randomly makes the program crash or produces garbage output)
		dummy=EMData(8,8)
		dummy.to_one()
		dummy['apix_x']=apix
		dummy['apix_y']=apix
		for j in range(nfs):
			dummy.write_image(outf,j)

		#EMAN2 does not allow stacks of images with different size; this, and possibly some bug, prevent images written from the parallelization task from
		#having the corret size if the pre-created dummy doesn't have the correct size to begin with. No point in writing big images for the dummy from the start.
		#re-writing the index=0 image will change the size of all images in the stack to the correct size
		dummy_correct_size = EMData(nx,ny)
		dummy_correct_size.to_one()
		dummy_correct_size['apix_x']=apix
		dummy_correct_size['apix_y']=apix
		dummy.write_image(outf,0)

		for i in range(nfs):
			if options.verbose:
				sys.stdout.write("\rstaging images ({}/{})".format(i+1,nfs))
				sys.stdout.flush()

			if options.parallel:
				#print "parallelism started"
				task = EraseGold2DTask( options, arg, i, outf)
				tasks.append(task)
			else:
				results=fiximage( options, arg, i, outf)

		if options.parallel:
			if tasks:
				tids = etc.send_tasks(tasks)
				if options.verbose:
					print "\n(erase_gold) %d tasks queued" % (len(tids))

				results = get_results( etc, tids, options )

		#if results:
		#	#pass
		#
		#	if '.ali' == originalarg[-4:] or '.mrc' == originalarg[-4:]:
		#		#intermediate = arg.replace('.hdf','.mrcs')
		#		finaloutput = arg.replace('.hdf',originalarg[-4:])
		#		cmd = 'e2proc2d.py ' + arg + ' ' + finaloutput + ' --twod2threed --outmode int16'
		#		runcmd(options,cmd)
		#		os.remove(arg)
		#
		#	if newarg: os.remove(newarg)

		if results:
			#pass

			if options.parallel:
				#outfstem = outf.replace('.hdf','')
				cmdbuildstack = 'e2buildstacks.py erasegold_tmp-*_proc.hdf --stackname ' + outf
				runcmd(options,cmdbuildstack)

				if options.debug:
					outfmasked = outf.replace('.hdf','_masked.hdf')
					cmdbuildstack = 'e2buildstacks.py erasegold_tmp-*_masked.hdf --stackname ' + outfmasked
					runcmd(options,cmdbuildstack)

					outfnoise= outf.replace('.hdf','_noise.hdf')
					cmdbuildstack = 'e2buildstacks.py erasegold_tmp-*_noise.hdf --stackname ' + outfnoise
					runcmd(options,cmdbuildstack)

			if '.ali' == originalarg[-4:] or '.mrc' == originalarg[-4:]:
				#intermediate = arg.replace('.hdf','.mrcs')
				finaloutput = outf.replace('.hdf',originalarg[-4:])
				cmd = 'e2proc2d.py ' + outf + ' ' + finaloutput + ' --twod2threed --outmode int16'
				
				#print "\ncomand to generate finaloutput",cmd
				runcmd(options,cmd)
				os.remove(arg)

			if newarg: 
				try:
					os.remove(newarg)
				except:
					try:
						#print "would have removed",newarg.replace('.hdf','_proc.hdf')
						os.remove(newarg.replace('.hdf','_proc.hdf'))
					except:
						pass
		try:
			filelist = [ tmpf for tmpf in os.listdir(".") if 'erasegold_tmp' in tmpf ]
			for tf in filelist:
			    os.remove(tf)
		except:
			print "WARNING: cleanup failed."


		dt = time.time() - t0
		if options.verbose:
			print("\n")
			sys.stdout.write("Erased fiducials from {} ({} minutes)\n".format(arg,round(dt/60.,2)))
	return


def fiximage(options,imgfile,imgindx,outf):
	#sys.stdout.write("\r{}/{}".format(i+1,nfs))
	#sys.stdout.flush()

	f = EMData(imgfile,imgindx) * -1
	f.process_inplace("normalize")

	sharp_msk, soft_msk = generate_masks(options,f)
	mskd_sharp = sharp_msk*f
	sub_sharp = f-mskd_sharp

	noise = local_noise(options,sub_sharp)
	outn = "{}_noise.hdf".format(imgfile)

	if options.debug: noise.write_image(outn,imgindx)

	mskd_soft = soft_msk*f
	sub_soft = f-mskd_soft

	if options.debug:
		sub_soft.write_image("{}_masked.hdf".format( os.path.splitext(imgfile), imgindx))

	result=sub_soft + noise * soft_msk
	result *= -1
	result.write_image(outf,imgindx)

	f *= -1
	#f.write_image("{}_compare.hdf".format( os.path.splitext(imgfile) ,imgindx)
	#result.write_image("{}_compare.hdf".format( os.path.splitext(imgfile) ,imgindx+1)
	#ctr+=2

	return


def generate_masks(options,img):
	img.process_inplace("normalize")
	# create sharp mask
	
	'''
	random cutoff values or int multiplication factors for thresholding ('3' below) work only on few datasets
	it's best to estimate these things based on the size of the feature to mask, and decide on a per-dataset basis
	how harshly to threshold
	'''
	#msk = img.process("filter.highpass.gauss",{"cutoff_pixels":25})
	#msk.process_inplace("filter.lowpass.gauss",{"cutoff_abs":0.1})
	#msk.process_inplace("threshold.clampminmax",{"maxval":msk["maximum"],"minval":msk["mean"]+3*msk["sigma"],"tozero":True})
	#msk.process_inplace("threshold.binary",{"value":msk["mean"]})]

	fourierpixels = img['nx']/2
	cutoffpixels = fourierpixels - options.goldsize/2
	
	msk = img.process("filter.highpass.gauss",{"cutoff_pixels":cutoffpixels})
	
	apix = img['apix_x']
	goldsizeinangstroms = apix*options.goldsize
	freq = 1.0/goldsizeinangstroms

	msk.process_inplace("filter.lowpass.tanh",{"cutoff_freq":freq})	#c:lowpass shouldn't be arbitrary; rather, use gold size to derive it.
	msk.process_inplace("threshold.clampminmax",{"maxval":msk["maximum"],"minval":msk["mean"]+options.nsigmas*msk["sigma"],"tomean":True})

	# remove dust
	if options.keepdust: sharp_msk = msk.copy()
	else:
		nproc = msk.numpy().copy()
		s = np.sqrt(options.goldsize*2).astype(int)
		se2=np.ones((s,s))
		nproc = ndimage.binary_closing(nproc,structure=se2).astype(int)
		nproc = ndimage.binary_opening(nproc,structure=se2).astype(int)
		sharp_msk = from_numpy(nproc)
	# grow slightly and create soft mask
	sharp_msk = sharp_msk.process("mask.addshells.gauss",{"val1":8,"val2":0})
	soft_msk = sharp_msk.process("mask.addshells.gauss",{"val1":0,"val2":8})
	return sharp_msk,soft_msk


def local_noise(options,img):
	localnoise = EMData(img["nx"],img["ny"])
	localnoise.to_zero()
	bs = options.boxsize
	nbxs = len(np.arange(-bs,img['nx']+bs,bs))*options.oversample
	nbys = len(np.arange(-bs,img['ny']+bs,bs))*options.oversample
	mx = np.linspace(0,img["nx"],nbxs).astype(int)
	my = np.linspace(0,img["ny"],nbys).astype(int)
	for x in mx:
		for y in my:
			r = img.get_clip(Region(x-bs/2,y-bs/2,bs,bs))
			n = EMData(bs,bs)
			n.to_zero()
			n.process_inplace("math.addnoise",{"noise":r["sigma_nonzero"]})
			
			#n.process_inplace("filter.highpass.gauss",{"cutoff_abs":0.01})
			
			fourierpixels = n['nx']/2
			cutoffpixels = fourierpixels - options.goldsize/2
			n.process_inplace("filter.highpass.gauss",{"cutoff_pixels":cutoffpixels})

			try: 
				n *= r["sigma_nonzero"]/n["sigma_nonzero"]
			except:
				if options.verbose > 8:
					print("WARNING: division by zero, from n['sigma_nonzero']={}".format(n["sigma_nonzero"]))
			n += r["mean_nonzero"]
			localnoise.insert_clip(n,(x-bs/2,y-bs/2))
	
	if options.lowpass:
		apix = img['apix_x']
		localnoise['apix_x'] = apix
		localnoise['apix_y'] = apix
		nyquistres=apix*2.0
		filtres=nyquistres*10.0/9.0
		filtfreq=1.0/filtres
		#if options.verbose: print("Apix {}\tResolution {}\tFrequency {}".format(apix,filtres,filtfreq))
		localnoise.process_inplace('filter.lowpass.tanh', {'cutoff_freq':filtfreq,'apix':apix})
	return localnoise


def runcmd(options,cmd):
	if options.verbose > 8: print("(erase_gold)(runcmd) running command: {}".format(cmd))
	p=subprocess.Popen( cmd, shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	text=p.communicate()
	p.stdout.close()
	if options.verbose > 8: print("(erase_gold)(runcmd) done")


def get_results(etc,tids,options):
	"""This will get results for a list of submitted tasks. Won't return until it has all requested results.
	aside from the use of options["ptcl"] this is fairly generalizable code. """

	# wait for them to finish and get the results
	# results for each will just be a list of (qual,Transform) pairs
	results=[0]*len(tids)		# storage for results
	ncomplete=0
	tidsleft=tids[:]
	while 1:
		time.sleep(5)
		proglist=etc.check_task(tidsleft)
		nwait=0
		for i,prog in enumerate(proglist):
			if prog==-1 : nwait+=1
			if prog==100 :
				r=etc.get_results(tidsleft[i])				#Results for a completed task

				if r:
					#print "r is", r
					#ptcl=r[0].classoptions["ptclnum"]		#Get the particle number from the task rather than trying to work back to it
					#results[ptcl] = r[1]
					results[i] = 1
				ncomplete+=1

		tidsleft=[j for i,j in enumerate(tidsleft) if proglist[i]!=100]		# remove any completed tasks from the list we ask about
		if options.verbose:
			sys.stdout.write("\r{} tasks\t{} complete\t{} in queue".format(len(tids),ncomplete,nwait))
			sys.stdout.flush()

		if len(tidsleft)==0: break

	return results


class EraseGold2DTask(JSTask):
	"""This is a task object for the parallelism system."""

	def __init__(self, options, imgfile, imgindx, outf):

		JSTask.__init__(self,"TomoBoxer3d",{},"")
		self.classoptions={"options":options,"imgfile":imgfile,"imgindx":imgindx,"outf":outf}

	def execute(self,callback=None):

		#print "insdie class EraseGold2DTask"
		options = self.classoptions['options']
		#print "options",options
		imgfile = self.classoptions['imgfile']
		imgindx = self.classoptions['imgindx']
		outf = self.classoptions['outf']

		#print "(EraseGold2DTask) imgindx",imgindx
		fiximage( options, imgfile, imgindx, outf)

		return



if __name__ == "__main__":
	main()



 # Redundant code. Staged for removal.

 		#if options.verbose: print("processing {} ({} images)".format(arg, EMUtil.get_image_count(arg)))

		#Averaging can be outsorced to e2proc2d via the command line, and the average can be read in as the new input
		#if options.average:

		#	newarg = arg.replace('.hdf','_avg.hdf')
		#	cmdavg = 'e2proc2d.py ' + arg + ' ' + newarg + ' --average'
		#	if ds > 1.0:
		#		cmdavg += ' --process math.fft.resample:n=' + str(ds)
		#	cmdavg += ' --process normalize'
		#	runcmd(options,cmdavg)
		#	arg = newarg

		#The code to operate on frame averages seems to be the same as that to operate on single images; no need for redundancy.
		# '''
		# 	avgr = Averagers.get("mean")
		# 	for i in range(EMUtil.get_image_count(fn)):
		# 		f = EMData(fn,i) * -1
		# 		if ds > 1.0: f.process_inplace("math.fft.resample",{"n":ds})
		# 		avgr.add_image(f)
		# 	img = avgr.finish()
		# 	img.process_inplace("normalize")

		# 	sharp_msk, soft_msk = generate_masks(options,img)
		# 	mskd_sharp = sharp_msk*img
		# 	sub_sharp = img-mskd_sharp
		# 	noise = local_noise(options,sub_sharp)

		# 	if options.debug: noise.write_image("{}_noise.hdf".format(arg))

		# 	mskd_soft = soft_msk*img
		# 	sub_soft = img-mskd_soft
		# 	result = sub_soft + noise * soft_msk
		# 	result *= -1

		# 	print("Writing result to {}".format(outf))

		# 	result.write_image(outf,0)
		# 	avg.write_image("{}_compare.hdf".format(arg),0)
		# 	result.write_image("{}_compare.hdf".format(arg),1)
		# '''
		#else:
		#ctr = 0
