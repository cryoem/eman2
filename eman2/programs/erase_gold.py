#!/usr/bin/env python

#====================
#Author: Michael Bell July, 2016 (edits, Jesus Galaz-Montoya). Last update: August, 2016
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
import subprocess
import os

from EMAN2jsondb import JSTask,jsonclasses


def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options] stack1.hdf stack2.mrcs ...

	Program to erase gold fiducials and other high-density features from images, such as frames in DDD movies or images in tiltseries. Requires scipy.
	"""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_argument("--average", default=False, action="store_true", help="Erase gold from average of input stack(s).")
	parser.add_argument("--lowpass", default=False, action="store_true", help="Also lowpass filter noise based on local properties. Useful for processing tomographic tilt series.")
	parser.add_argument("--keepdust", default=False, action="store_true", help="Do not remove 'dust' from mask (include objects smaller than gold fiducials).")
	parser.add_argument("--goldsize", default=30, type=float, help="Diameter (in pixels) of gold fiducials to erase.")
	parser.add_argument("--downsample", default=1.0, type=float, help="Downsample the input stack(s). Default is 1, i.e. no downsampling.")
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

	logger = E2init(sys.argv, options.ppid)
	print "\n(e2tomopreproc)(main) started log"	

	if options.parallel == 'None' or options.parallel == 'none':
		options.parallel == None

	if options.parallel:
		from EMAN2PAR import EMTaskCustomer
		etc=EMTaskCustomer(options.parallel)

	for arg in args:
		newarg=''
		originalarg = arg

		hdr = EMData(arg,0,True) #load header only to get parameters used below
		apix = hdr['apix_x']
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

			#make the new stack of 2D images (dumy_stack.hdf) the new input (the name of the input file but with .hdf format); this intermediate file will be deleted in the end.
			newarg = arg.replace(arg[-4:],'.hdf')
			os.rename('dummy_stack.hdf',newarg)
			arg = newarg

		if options.verbose: print("processing {} ({} images)".format(arg, EMUtil.get_image_count(arg)))
		
		#Averaging can be outsorced to e2proc2d via the command line, and the average can be read in as the new input
		if options.average:
			
			newarg = arg.replace('.hdf','_avg.hdf')
			
			cmdavg = 'e2proc2d.py ' + arg + ' ' + newarg + ' --average'

			if ds > 1.0:
				cmdavg += ' --process math.fft.resample:n=' + str(ds)

			cmdavg += ' --process normalize'

			runcmd(options,cmdavg)

			arg = newarg

		#The code to operate on frame averages seems to be the same as that to operate on single images; no need for redundancy.
		'''
			avgr = Averagers.get("mean")
			for i in range(EMUtil.get_image_count(fn)):
				f = EMData(fn,i) * -1
				if ds > 1.0: f.process_inplace("math.fft.resample",{"n":ds})
				avgr.add_image(f)
			img = avgr.finish()
			img.process_inplace("normalize")

			sharp_msk, soft_msk = generate_masks(options,img)
			mskd_sharp = sharp_msk*img
			sub_sharp = img-mskd_sharp
			noise = local_noise(options,sub_sharp)

			if options.debug: noise.write_image("{}_noise.hdf".format(arg))

			mskd_soft = soft_msk*img
			sub_soft = img-mskd_soft
			result = sub_soft + noise * soft_msk
			result *= -1

			print("Writing result to {}".format(outf))

			result.write_image(outf,0)
			avg.write_image("{}_compare.hdf".format(arg),0)
			result.write_image("{}_compare.hdf".format(arg),1)
		'''
		#else:
		#ctr = 0

		outf = "{}_proc.hdf".format( os.path.splitext(arg)[0] )

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

		print "outf",outf

		if options.parallel:
			cmdunstacking = 'e2proc2d.py ' + arg + ' erasegold_tmp.hdf --unstacking'
			runcmd(options,cmdunstacking)

		if options.subset:
			nfs=options.subset

		for i in range(nfs):
			
				#if i > options.subset -1:
				#	break

			if options.verbose: print "processing image {}/{}".format(i,nfs)
			
			if options.parallel:
				print "parallelism started"
				thisimg = 'erasegold_tmp-' + str(i+1).zfill(len(str(nfs))) + '.hdf'			#c: when e2proc2d.py unstacks images, it starts from 1, not from 0
				thisoutf = 'erasegold_tmp-' + str(i+1).zfill(len(str(nfs))) + '_proc.hdf'
				task = EraseGold2DTask( options, thisimg, 0, thisoutf,nfs)
				tasks.append(task)
			else:
				results=fiximage( options, arg, i, outf,nfs)

		if options.parallel:	
			if tasks:
				tids = etc.send_tasks(tasks)
				if options.verbose: 
					print "\n(erase_gold)(main) preprocessing %d tasks queued" % (len(tids)) 

				results = get_results( etc, tids, options )

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

	
	E2end(logger)

	return


def fiximage(options,imgfile,imgindx,outf,nfs):
	#sys.stdout.write("\r{}/{}".format(i+1,nfs))
	#sys.stdout.flush()
	#print "fiximage function received imgfile, outf, imgindx",imgfile, outf, imgindx
	
	f = EMData(imgfile,imgindx) * -1
	#print "loaded image"
	
	if options.downsample > 1.0:
		f.process_inplace("math.fft.resample",{"n":options.downsample})
	f.process_inplace("normalize")

	#print "downsampled and normalized"

	sharp_msk, soft_msk = generate_masks(options,f)
	mskd_sharp = sharp_msk*f
	sub_sharp = f-mskd_sharp

	#print "generated masks"

	noise = local_noise(options,sub_sharp)

	outnoise = "{}_noise.hdf".format( os.path.splitext(imgfile)[0] )

	#print "generated noise"

	if options.debug: noise.write_image(outnoise,imgindx)

	mskd_soft = soft_msk*f
	sub_soft = f-mskd_soft

	if options.debug:
		outmasked = "{}_masked.hdf".format( os.path.splitext(imgfile)[0] )
		sub_soft.write_image(outmasked, imgindx)

	#print "wrote out noise and mask"
	result=sub_soft + noise * soft_msk
	result *= -1
	#print "generated result"

	#if options.parallel:
		#outfindividual=outf.replace('.hdf', str(imgindx).zfill(len(str(nfs))) + '_tmp.hdf')
		#print "writing out result to", outfindividual
		#result.write_image(outfindividual,0)
		#print "wrote out image", outfindividual
	#else:
	result.write_image(outf,imgindx)

	#print "wrote out result"
	f *= -1
	#f.write_image("{}_compare.hdf".format( os.path.splitext(imgfile) ,imgindx)
	#result.write_image("{}_compare.hdf".format( os.path.splitext(imgfile) ,imgindx+1)
	#ctr+=2

	return 1


def generate_masks(options,img):
	img.process_inplace("normalize")
	# create sharp mask
	
	#msk = img.process("filter.highpass.gauss",{"cutoff_pixels":25})
	fourierpixels = img['nx']/2
	pixels = fourierpixels - options.goldsize/2
	
	msk = img.process("filter.highpass.gauss",{"cutoff_pixels":pixels})
	apix = img['apix_x']
	
	goldinAngs=apix*options.goldsize
	freq=1.0/goldinAngs

	#msk.process_inplace("filter.lowpass.gauss",{"cutoff_abs":0.1})
	msk.process_inplace("filter.lowpass.tanh",{"cutoff_freq":freq})	#c:lowpass shouldn't be arbitrary; rather, use gold size to derive it.

	#msk.process_inplace("threshold.clampminmax",{"maxval":msk["maximum"],"minval":msk["mean"]+3*msk["sigma"],"tozero":True})
	#msk.process_inplace("threshold.clampminmax",{"maxval":msk["mean"]+3*msk["sigma"],"minval":msk["mean"]-3*msk["sigma"],"tomean":True})
	
	msk.process_inplace("threshold.clampminmax",{"maxval":msk["maximum"],"minval":msk["mean"]+options.nsigmas*msk["sigma"],"tomean":True})

	msk.process_inplace("threshold.binary",{"value":msk["mean"]})
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
	sharp_msk = sharp_msk.process("mask.addshells.gauss",{"val1":1,"val2":0})
	soft_msk = sharp_msk.process("mask.addshells.gauss",{"val1":0,"val2":4})
	return sharp_msk,soft_msk


def local_noise(options,img):
	localnoise = EMData(img["nx"],img["ny"])
	apix = img['apix_x']
	
	localnoise['apix_x']=apix
	localnoise['apix_y']=apix
	
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

			rth=r.process("threshold.clampminmax",{"maxval":r["mean"]+3*r['sigma'],"minval":r["mean"]-3*r["sigma"],"tozero":True})

			n.process_inplace("math.addnoise",{"noise":rth["sigma_nonzero"]})
			#n.process_inplace("normalize.toimage",{"to":rth})

			#n.process_inplace("filter.highpass.gauss",{"cutoff_abs":0.01})
			fourierpixels = n['nx']/2
			pixels = fourierpixels - options.goldsize/2
			n.process_inplace("filter.highpass.gauss",{"cutoff_pixels":pixels})

			try: n *= r["sigma_nonzero"]/n["sigma_nonzero"]
			except:
				if options.verbose > 9:
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
	
	localnoise.process_inplace("normalize.toimage",{"to":img})

	return localnoise


def runcmd(options,cmd):
	if options.verbose > 8:
		print "\n(erase_gold)(runcmd) running command", cmd

	p=subprocess.Popen( cmd, shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	text=p.communicate()
	p.stdout.close()

	if options.verbose > 8:
		print "\n(erase_gold)(runcmd) done"


'''
CLASS TO PARALLELIZE GOLD ERASING
'''
class EraseGold2DTask(JSTask):
	'''This is a task object for the parallelism system.'''

	def __init__(self, options, imgfile, imgindx, outf,nfs):
	
		JSTask.__init__(self,"TomoBoxer3d",{},"")
		self.classoptions={"options":options,"imgfile":imgfile,"imgindx":imgindx,"outf":outf,"nfs":nfs}
	
	def execute(self,callback=None):
		
		#print "insdie class EraseGold2DTask"
		options = self.classoptions['options']
		#print "options",options
		imgfile = self.classoptions['imgfile']
		imgindx = self.classoptions['imgindx']
		outf = self.classoptions['outf']
		nfs = self.classoptions['nfs']
		#print "BEFORE fiximage in class EraseGold2DTask, input, output, imgindx", imgfile,outf,imgindx

		#print "(EraseGold2DTask) imgindx",imgindx
		ret=fiximage( options, imgfile, imgindx, outf,nfs)
		#print "AFTER fiximage in class EraseGold2DTask, ret=",ret
		
		return


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
			print ("\n%d tasks, %d complete, %d waiting to start \r" % (len(tids),ncomplete,nwait))
			sys.stdout.flush()
	
		if len(tidsleft)==0: break
		
	return results





if __name__ == "__main__":
	main()
