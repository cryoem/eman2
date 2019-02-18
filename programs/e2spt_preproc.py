#!/usr/bin/env python
#
# Author: Jesus Galaz-Montoya 03/2011, 
# (based on Steven Ludtke's initial implementation [02/15/2011] of Jesus's older scripts, from M.F.Schmid's methods).
# Last modification: July/08/2015
#
# Copyright (c) 2011 Baylor College of Medicine
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
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  2111-1307 USA
#
#
from __future__ import print_function
from __future__ import division
from builtins import range
from EMAN2 import *
from EMAN2_utils import *
from EMAN2jsondb import JSTask,jsonclasses
import os
from sys import argv

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog <output> [options]

	This program is used to preprocess subtomograms before aligning them. The same can be accomplished with 
	e2proc3d, except that this program is parallelized and thus should be substantially faster for large subtomograms.
	"""
			
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	
	parser.add_argument("--apix",type=float,default=0.0,help="""Default=0.0 (not used). Use this apix value where relevant instead of whatever is in the header of the reference and the particles. Will overwrite particle header as well.""")

	parser.add_argument("--clip",type=int,default=0,help="""Default=0 (which means it's not used). Boxsize to clip particles. For example, the boxsize of the particles might be 100 pixels, but the particles are only 50 pixels in diameter. Aliasing effects are not always as deleterious for all specimens, and sometimes 2x padding isn't necessary.""")

	parser.add_argument("--highpass",type=str,default='',help="""Default=None. A highpass filtering processor (see 'e2help.py processors -v 10' at the command line) to be applied to each volume prior to COARSE alignment. Not applied to aligned particles before averaging.""")

	parser.add_argument("--input", type=str, default='',help="""Default=None. The name of the input volume stack. MUST be HDF since volume stack support is required.""")
	
	parser.add_argument("--lowpass",type=str,default='',help="""Default=None. A lowpass filtering processor (see 'e2help.py processors -v 10' at the command line) to be applied to each volume prior to COARSE alignment. Not applied to aligned particles before averaging.""")

	parser.add_argument("--mask",type=str,default='', help="""Default=None. Masking processor applied to particles before alignment. IF using --clip, make sure to express outer mask radii as negative pixels from the edge.""")
	parser.add_argument("--maskfile",type=str,default='',help="""Default=None. Mask file (3D IMAGE) applied to particles before alignment. Must be in HDF format. Default is None.""")
	
	parser.add_argument("--nopath",action='store_true',default=False,help="""If supplied, this option will save results in the directory where the command is run. A directory to store the results will not be made.""")
	parser.add_argument("--normproc",type=str, default='',help="""Default=None (see 'e2help.py processors -v 10' at the command line). Normalization processor applied to particles before alignment. If normalize.mask is used, results of the mask option will be passed in automatically. If you want to turn this option off specify \'None\'""")
	
	parser.add_argument("--output", type=str, default='',help="""Default=None. Specific name of HDF file to write processed particles to.""")

	parser.add_argument("--path",type=str,default='spt_preproc',help="""Default=spt. Directory to store results in. The default is a numbered series of directories containing the prefix 'sptpreproc'; for example, sptpreproc_02 will be the directory by default if 'sptpreproc_01' already exists.""")
	parser.add_argument("--parallel",type=str, default='', help="""default=None. Parallelism. See http://blake.bcm.edu/emanwiki/EMAN2/Parallel""")
	parser.add_argument("--ppid", type=int, help="""Default=-1. Set the PID of the parent process, used for cross platform PPID""",default=-1)
	parser.add_argument("--preprocess",type=str,default='',help="""Any processor (see 'e2help.py processors -v 10' at the command line) to be applied to each volume prior to COARSE alignment. Not applied to aligned particles before averaging.""")

	parser.add_argument("--shrink", type=int,default=0,help="""Default=0 (no shrinking). Optionally shrink the input volumes by an integer amount for coarse alignment.""")
	parser.add_argument("--subset",type=int,default=0,help="""Default=0 (not used). Refine only this substet of particles from the stack provided through --input""")
	
	parser.add_argument("--threshold",type=str,default='',help="""Default=None. A threshold applied to the subvolumes after normalization. For example, --threshold=threshold.belowtozero:minval=0 makes all negative pixels equal 0, so that they do not contribute to the correlation score.""")
	
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="""Default=0. Verbose level [0-9], higner number means higher level of verboseness; 10-11 will trigger many messages that might make little sense since this level of verboseness corresponds to 'debugging mode'""")

	(options, args) = parser.parse_args()
	
	logger = E2init(sys.argv, options.ppid)
	
	print("\n(e2spt_preproc)(main) started log")
		
	if options.path and not options.nopath:
	
		options = makepath(options,'spt_preproc')

	if options.parallel=='None' or options.parallel=='none':
		options.parallel=None
	
	options = checkinput( options )
			
	if options.mask or options.maskfile or options.threshold or options.clip or options.threshold or options.normproc or options.preprocess or options.lowpass or options.highpass or int(options.shrink) > 1:
		
		preprocstack =  str(os.path.basename(options.input).replace('.hdf','_preproc.hdf'))
		
		if options.path and not options.nopath:
			preprocstack = options.path + '/' + preprocstack
		
		if options.output:
			if '.hdf' in options.output[-4:]:
				preprocstack = options.output
			else:
				print("\n(e2spt_preproc)(main) ERROR: '.hdf' must be the last four characters of the output filename.")
			
		print("\n(e2spt_preproc)(main) output stack will be %s" %( preprocstack))

		n = 0
		try:
			n = EMUtil.get_image_count( options.input )
		except:
			print("\n(e2spt_preproc)(main) ERROR: --input stack seems to be invalid")
			sys.exit()
		
		print("\n(e2spt_preproc)(main) number of particles is %d" %( n)) 
		
		
		c = os.getcwd()
		
		findir = os.listdir( c )
		
		if preprocstack not in findir:
		
			dimg = EMData(8,8,8)
			dimg.to_one()

			for i in range(n):
				dimg.write_image( preprocstack, i )
		
		else:
			print("\n(e2spt_preproc)(main) WARNING: a file with the name of the output stack %s is already in the current directory and will be DELETED" %( preprocstack ))
			os.remove( preprocstack )
			
			dimg = EMData(8,8,8)
			dimg.to_one()

			for i in range(n):
				dimg.write_image( preprocstack, i )
		

		finalbox = EMData(options.input,0,True)['nx']
		if options.clip:
			finalbox=options.clip

		
		#dimglarge=EMData(finalbox,finalbox,finalbox)
		#dimglarge.to_one()
		#dimglarge.write_image(preprocstack,0)
		#n=EMUtil.get_image_count(options.input)
		#if options.subset:
		#	n=options.subset
		#dimglarge.write_image(preprocstack,n-1)

		if options.verbose:
			print("\n(e2spt_preproc)(main) wrote dummy ptcls to %s" %( preprocstack))
	
		
		print("\n(e2spt_preproc)(main) - INITIALIZING PARALLELISM!\n")
		
		if options.parallel:
			from EMAN2PAR import EMTaskCustomer
			etc=EMTaskCustomer(options.parallel)
			pclist=[options.input]

			etc.precache(pclist)
			print("\n(e2spt_preproc)(main) - precaching --input")

			tasks=[]
			results=[]
				
		options = sptOptionsParser( options )
		
		
		for j in range(n):
			#print "processing  particle", j
			
			img = EMData( options.input, j )
			
			if options.parallel:
				#task = Preproc3DTask( ["cache",options.input,j], options, j, preprocstack )
				task = Preproc3DTask( img, options, j, preprocstack )
				tasks.append(task)
		
			else:
				img = EMData( options.input, j )
				pimg = preprocfunc( img, options, j, preprocstack)
								
		
		
		if options.parallel and tasks:
			tids = etc.send_tasks(tasks)
			if options.verbose: 
				print("\n(e2spt_preproc)(main) preprocessing %d tasks queued" % (len(tids))) 

	
			results = get_results( etc, tids, options )
		#print "\n(e2spt_preproc)(main) preprocessing results are", results	
		
		
		#print "\n(e2spt_preproc)(main) input changing to preprocstack"
		#options.input = preprocstack

		#cache needs to be reloaded with the new options.input		
		
	else:
		print("\n(e2spt_preproc)(main) Nothing to do. No preprocessing parameters specified.")
		
	E2end(logger)
	
	return
	


'''
CLASS TO PARALLELIZE PREPROCESSING
'''


class Preproc3DTask(JSTask):
	'''This is a task object for the parallelism system. It is responsible for preprocessing one 3-D volume, with a variety of options'''


	def __init__(self,image,options,ptclnum,outname):
	#def __init__(self,options,ptclnum,outname):
		data={"image":image}
		
		JSTask.__init__(self,"Preproc3d",data,{},"")

		self.classoptions={"options":options,"ptclnum":ptclnum,"outname":outname}
	
	def execute(self,callback=None):
		classoptions=self.classoptions
		options=self.classoptions['options']
		outname = self.classoptions['outname']
		
		i = self.classoptions['ptclnum']
		#image=self.data["image"]
		#print "simage is", simagefile
		print("\n(e2spt_preproc)(Preproc3DTask) ptclnum is i",i)
		
		if isinstance(self.data["image"],EMData):
			simage=self.data["image"]
		else: 
			simage=EMData(self.data["image"][1],self.data["image"][2])
			print("\n(e2spt_preproc)(Preproc3DTask) image was actually file",image)
			
		#simage = self.data['image']
		
		if options.verbose:
			print("\n(e2spt_preproc)(Preproc3DTask) outname %s and i=%d " %(outname, i))
		
		
		preprocfunc( simage, options, i, outname )
		#preprocfunc( options, i, outname )
		
		return


def preprocfunc( simage, options, i, outname, simulation=False, resizeonly=False, ref=False ):
	#def preprocfunc( options, i, outname, simulation=False, resizeonly=False ):
		
	#simage = EMData(options.input,i)
	#if fyle:
	#	simage = fyle.copy()

	if options.verbose:
		print("\n(e2spt_preproc) preprocessing particle", i)


	if not resizeonly:
	
		apix = simage['apix_x']
	
		if options.verbose > 9:
			print("\n(e2spt_preproc)(preprocfunc)apix is %f" % (apix))
	
		'''
		Make the mask first 
		'''
		if options.verbose > 9:
			print("\n(e2spt_preproc)(preprocfunc) masking particle", i)
	
		maskimg = EMData( int(simage["nx"]), int(simage["ny"]), int(simage["nz"]) )
		maskimg.to_one()
		print("\n(e2spt_preproc)(preprocfunc) done creating mask")

	
		if options.mask and options.mask != 'None' and options.mask != 'none':
			#if options.verbose:
			#print "\n(e2spt_preproc)(preprocfunc) this is the mask I will apply: %s,%s" %(options.mask[0],options.mask[1]) 
			print("masking")
			maskimg.process_inplace(options.mask[0],options.mask[1])
	
			#print ("\n(e2spt_preproc)(preprocfunc) --mask provided: %s" %( options.mask))
			#mask.write_image(options.path + '/mask.hdf',-1)

		if ref:
			try:
				if options.maskfile:
					maskfileimg = EMData(options.maskfile,0)
			
					if maskfileimg['nx'] !=  maskimg['nx'] or maskfileimg['ny'] !=  maskimg['ny'] or maskfileimg['nz'] !=  maskimg['nz']:
						maskfileimg = clip3d( maskfileimg, maskimg['nx'] )
		
					maskimg.mult( maskfileimg )
				
					if options.verbose > 9:
						print("\n(e2spt_preproc)(preprocfunc) including --maskfile in mask")
					#print "\n(e2spt_preproc)(preprocfunc)a maskfile was multiplied by the mask %s" %( options.maskfile) 
				else:
					if options.verbose > 9:
						print("\n(e2spt_preproc)(preprocfunc) apparently, there was no --maskfile")
					#pass
			except:
				pass
	
	
		'''
		Set the 'mask' parameter for --normproc if normalize.mask is being used
		'''
		if options.normproc and options.normproc != 'None' and options.normproc != 'none':
			#print "normproc is %s" %(options.normproc[0]) 
			if options.normproc[0]=="normalize.mask": 
				options.normproc[1]["mask"]=maskimg

		'''
		Normalize-Mask
		'''	
		
		if options.normproc and options.normproc != 'None' and options.normproc != 'none':
			simage.process_inplace(options.normproc[0],options.normproc[1])
			#simage.write_image(options.path + '/imgMsk1norm.hdf',-1)
			print("normalizing")
			#print "\n(e2spt_preproc)(preprocfunc) --normproc provided: %s" %( str(options.normproc[0]) ) 

		try:
			#if mask and mask != 'None' and mask != 'none' or options.maskfile:
			if options.mask or options.maskfile:
				if options.verbose > 9:
					print("\n(e2spt_preproc)(preprocfunc) masking again after normalizing")
				simage.mult(maskimg)
				#simage.write_image(options.path + '/imgMsk1normMsk2.hdf',-1)
		except:
			pass
	
		'''
		Any threshold picked visually was based on the original pixels give the original
		box size; therefore, it needs to be applied before clipping the box
		'''
		if options.threshold and options.threshold != 'None' and options.threshold != 'none':
			if options.verbose > 9:
				print("\n(e2spt_preproc)(preprocfunc) --thresholding provided")
			simage.process_inplace( options.threshold[0], options.threshold[1] )

		if options.shrink  == 'None' or options.shrink == 'none':
			options.shrink = None
	
		if options.lowpass  == 'None' or options.lowpass == 'none':
			options.lowpass = None

		if options.highpass == 'None' or options.highpass == 'none':
			options.highpass = None

		if options.preprocess == 'None' or options.preprocess == 'none':
			options.preprocess = None

	
	
		if options.lowpass:
			if options.verbose > 9:
				print("lowpassing particle",i)
			#print "\n(e2spt_preproc)(preprocfunc) --lowpass provided: %s, %s" %( options.lowpass[0],options.lowpass[1]) 
			simage.process_inplace(options.lowpass[0],options.lowpass[1])
			#fimage.write_image(options.path + '/imgPrepLp.hdf',-1)

		if options.highpass:
			if options.verbose > 9:
				print("highpassing particle",i) 
			#print "\n(e2spt_preproc)(preprocfunc) --highpass provided: %s, %s" %( options.highpass[0],options.highpass[1]) 
			simage.process_inplace(options.highpass[0],options.highpass[1])
			#fimage.write_image(options.path + '/imgPrepLpHp.hdf',-1)

		if options.preprocess:
			if options.verbose > 9:
				print("preprocessing particle",i)
			#print "\n(e2spt_preproc)(preprocfunc) --preprocess provided: %s, %s" %( options.preprocess[0],options.preprocess[1])
			simage.process_inplace(options.preprocess[0],options.preprocess[1])
			#fimage.write_image(options.path + '/imgPrep.hdf',-1)
	
	
	
	'''
	If the box is clipped, you need to make sure rotations won't induce aberrant blank
	corners; therefore, mask again softly at radius -4, independently of whatever --mask might have been provided or not
	'''
	if options.clip and options.clip != 'None' and options.clip != 'none':
		#if simage['nx'] != options.clip or simage['ny'] != options.clip or simage['nz'] != options.clip:
		
		clipf = options.clip
		
		#if not simulation:
		#	if options.shrink > 1:
		#		clipf /= options.shrink
		
		#	simage = clip3D( simage, clipf)
		#else:
		
		simage = clip3d( simage, clipf)
			
		simage.process_inplace('mask.soft',{'outer_radius':-4})
	
	
	
	if options.shrink and int( options.shrink  ) > 1:
		if options.verbose > 9:
			print("shrinking particle", i)
		#print "\n(e2spt_preproc)(preprocfunc)--shrink provided: %d" %(options.shrink)
		simage.process_inplace("math.fft.resample",{"n":options.shrink })
	
	
	try:
		if options.output:
			if options.verbose > 9:
				print("\n(e2spt_preproc)(preprocfunc) outname  %s and i=%d"  % (options.output, i))
			simage.write_image( outname, i )
		else:
			if options.verbose > 9:
				print("\n(e2spt_preproc)(preprocfunc) no --output provided. Default outname %s and i=%d  " % ( outname, i)) 
			simage.write_image( outname, i )
	except:
		if options.verbose > 9:	
			print("\n(e2spt_preproc)(preprocfunc) parameter --output probably doesn't exist in program calling this function. default outname %s and i=%d " %( outname , i)) 
		
		if outname:
			simage.write_image( outname, i )
		else:
			print("no output name")
		
		
	#del simage		
	return simage


'''
def clip3D( vol, sizex, sizey=0, sizez=0 ):
	
	if not sizey:
		sizey=sizex
	
	if not sizez:
		sizez=sizex
	
	volxc = vol['nx']/2
	volyc = vol['ny']/2
	volzc = vol['nz']/2
	
	print("clipping volume to size", sizex, sizey, sizez)
	
	Rvol =  Region( (2*volxc - sizex)/2, (2*volyc - sizey)/2, (2*volzc - sizez)/2, sizex , sizey , sizez)
	vol.clip_inplace( Rvol )
	#vol.process_inplace('mask.sharp',{'outer_radius':-1})
	
	return vol
'''


"""
def get_results(etc,tids,verbose):
	'''This will get results for a list of submitted tasks. Won't return until it has all requested results.
	aside from the use of options["ptcl"] this is fairly generalizable code.'''
	
	print "(e2spt_classaverage)(get_results_preproc)"
	
	# wait for them to finish and get the results
	# results for each will just be a list of (qual,Transform) pairs
	results=[0]*len(tids)		# storage for results
	ncomplete=0
	tidsleft=tids[:]
	
	print "results len is", len(results)
	
	while 1:
		time.sleep(5)
		proglist=etc.check_task(tidsleft)
		nwait=0
		for i,prog in enumerate(proglist):
			if prog==-1 : nwait+=1
			if prog==100 :
				r = etc.get_results(tidsleft[i])				#Results for a completed task
				
				print "r is", r
				
				if r:
					#print "r is", r
					ptcl=r[0].classoptions["ptclnum"]		#Get the particle number from the task rather than trying to work back to it
					results[ptcl] = r[1]
					
					print "ptcl is", ptcl
					print "results inside get_results are", results
										
				ncomplete+=1
		
		tidsleft=[j for i,j in enumerate(tidsleft) if proglist[i]!=100]		# remove any completed tasks from the list we ask about
		if verbose:
			print "  %d tasks, %d complete, %d waiting to start        \r"%(len(tids),ncomplete,nwait)
			sys.stdout.flush()
	
		if len(tidsleft)==0: break
		
	return results
"""

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
					ptcl=r[0].classoptions["ptclnum"]		#Get the particle number from the task rather than trying to work back to it
					results[ptcl] = r[1]					
				ncomplete+=1
		
		tidsleft=[j for i,j in enumerate(tidsleft) if proglist[i]!=100]		# remove any completed tasks from the list we ask about
		if options.verbose:
			print(("%d tasks, %d complete, %d waiting to start \r" % (len(tids),ncomplete,nwait)))
			sys.stdout.flush()
	
		if len(tidsleft)==0: break
		
	return results


	
if __name__ == "__main__":
    main()
    
