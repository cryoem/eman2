#!/usr/bin/env python

#
# Author: Jesus G. Galaz 9/1/2010 - Last Update 11/14/2012 
# Copyright (c) 2011- Baylor College of Medicine
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

import sys
import os, commands
from EMAN2 import *
import subprocess

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options] <Volume file>
	Files to be processed must NOT contain ANY dots '.' other than the one preceding the file's format.
	Prepares one volume (the "reference") to be at the same apix and boxsize than another (a "particle" you might want to align to the reference).
	Takes PDB, MRC and HDF file formats, but outputs HDF or MRC only. The reference can also be subject to symmetrization, filters, and other processing options. 
	"""
	
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument('--stack2process', type=str, default='', help="""Path to the stack 
		that needs to be processed to match stack2match. If you want to process multiple stacks 
		or files, just separate them by commas --stacks2process=vo1.mrc,vol2.hdf,file.pdb""")
		
	parser.add_argument('--stack2match', type=str, default='', help="Path to the image or stack of images which stack2process will match after processing.")
	parser.add_argument('--output', type=str, default='', help="Specify the name of the file to which the edited stack2process will be written.")
	parser.add_argument("--boxsize",type=int,default=0,help="If NOT specified, the reference's box size will match that of the data. If specified, both the reference and the data will be resized.")
	parser.add_argument("--sym", type=str, default='', help='Will apply the specified symmetry to the edited stack2process.')
	
	parser.add_argument("--apix",type=float,default=0.0,help="""If specified, the program 
		will assume this is the correct apix for the stack2match, which will be written to 
		the data's header if the data is in .hdf format, and will also be used to scale 
		the files in stack2process.""")
	
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	
	parser.add_argument("--preprocess",type=str,help="Any processor (as in e2proc3d.py) to be applied to the edited stack2process.", default='')
	parser.add_argument("--lowpass",type=str,help="A lowpass filtering processor (as in e2proc3d.py) be applied to the edited stack2process.", default='')
	parser.add_argument("--highpass",type=str,help="A highpass filtering processor (as in e2proc3d.py) to be applied to the edited stack2process.", default='')
	parser.add_argument("--mask",type=str,help="Mask processor applied to particles before alignment. Default is None", default='')
	parser.add_argument("--normproc",type=str,default='',help="""Normalization processor applied to 
		particles before alignment. Default is None. If normalize.mask is used, 
		results of the mask option will be passed in automatically.""")

	parser.add_argument("--threshold",type=str,default='',help="""EMAN2 processor to be used
		to threshold the stack2process. See available thresholds by typing 
		e2help.py processors --verbose=10
		at the command line.""")
		
	parser.add_argument("--shrink",type=float,default=0.0,help="""This will scale the stack2process
		by the factor specified. This does NOT need to be an integer. You can enter 1.5, 2.7, etc.
		any decimal number > 1.0 is valid.""")	
	
	parser.add_argument("--mirror", action="store_true", help="Will generate a mirrored copy of the edited stack2process.", default=False)
	parser.add_argument("--sharpfiltres",type=float,default=0.0,help="If specified, the edited stack2process will be sharply filtered to this resolution.")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)

	(options, args) = parser.parse_args()

	logger = E2init(sys.argv, options.ppid)
	
	'''
	Check for sane formats for stack2process and output
	'''
	
	stacks2process = options.stack2process.split(',')
	
	for stack2p in stacks2process:
		
		options.stack2process = stack2p
		
		if len(stacks2process) > 1 or not options.output:
			options.output = stack2p.split('.')[0] + '_PREP.hdf'
			
		print "Options.output is", options.output

		if '.hdf' not in options.stack2process and '.pdb' not in options.stack2process and '.mrc' not in options.stack2process:
			print "(e2spt_refprep.py) ERROR: The format of the reference to edit must be .pdb, .mrc or .hdf"
			sys.exit()
	
		if '.hdf' not in options.output and '.mrc' not in options.output:
			print "(e2spt_refprep.py) ERROR: The output reference must be written to either an mrc or an hdf file."
			sys.exit() 
	
		check=0
	
		current = os.getcwd()
		findir = os.listdir(current)
	
		if options.stack2process.split('/')[-1] not in findir:
			os.system('cp ' + options.stack2process + ' ' + current)
			options.stack2process = current + '/' + options.stack2process.split('/')[-1]
		
		if '.pdb' in options.stack2process:
			pdbmodel = options.stack2process
			pdbmodel = pdbmodel.split('/')[-1]
			mrcmodel = pdbmodel.replace('.pdb','.mrc')
		
			cmd = 'e2pdb2mrc.py ' + pdbmodel + ' ' + mrcmodel
			p=subprocess.Popen( cmd, shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			text=p.communicate()	
			p.stdout.close()

			options.stack2process = mrcmodel
			check=1
			if options.verbose:
				print "(e2spt_refprep.py) I've converted the .pdf reference to .mrc"
	
		if '.mrc' in options.stack2process:
			mrcmodel = options.stack2process
			if check==0:
				mrcmodel = mrcmodel.split('/')[-1]
			hdfmodel = mrcmodel.replace('.mrc','.hdf')
		
			cmd = 'e2proc3d.py ' + mrcmodel + ' ' + hdfmodel + ' && rm ' + mrcmodel
			p=subprocess.Popen( cmd, shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			text=p.communicate()	
			p.stdout.close()
		
			options.stack2process = hdfmodel
			check=1
			if options.verbose:
				print "(e2spt_refprep.py) I've converted the .mrc reference to .hdf"
	
	
		nStack2process = EMUtil.get_image_count( options.stack2process )
	
		if options.verbose:
			print "\n\n(e2spt_refprep.py) There are these many particles in stack2process to be edited",nStack2process
	
		'''
		Make all sides of the reference box equally sized and even, if they're not
		'''
	
		stack2processEd = options.stack2process.replace('.hdf','_ed.hdf')
		if options.output:
			stack2processEd = options.output
	
	
		#print "(e2spt_refprep.py) stack2processEd is", stack2processEd
	
		stack2processSample = EMData( options.stack2process, 0, True)
	
		cmd = 'e2proc3d.py ' + options.stack2process + ' ' + stack2processEd + ' && e2fixheaderparam.py --input=' + stack2processEd + ' --stem=origin --stemval=0' 
		
		
		#print "(e2spt_refprep.py) This command will copy and create it", cmd
	
		p=subprocess.Popen( cmd, shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		text=p.communicate()	
		p.stdout.close()
	
		#print "(e2spt_refprep.py) Command executed. Feedback was", text
	
		if stack2processSample['nx'] != stack2processSample['ny'] or stack2processSample['nx'] != stack2processSample['nz'] or stack2processSample['ny'] != stack2processSample['nz']:	
			x0 = stack2processSample['nx']
			y0 = stack2processSample['ny']
			z0 = stack2processSample['nz']
	
			side = max(x0,y0,z0)
			if side % 2:
				side += 1
		
			cmd = 'e2proc3d.py ' + stack2processEd + ' ' + stack2processEd + ' --clip=' + str(side)
			p=subprocess.Popen( cmd, shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			text=p.communicate()	
			p.stdout.close()
		
			#R = Region((x0 - side)/2, (y0 - side)/2, (z0 - side)/2, side, side, side)
			#ref.clip_inplace(R)
			
			if options.verbose:
				print "\n(e2spt_refprep.py) I have clipped the particles in stack2process to have a cubical and even-sided box size."

		'''
		If options.apix is supplied, first change the apix of stack2match to that; then make stack2process match it
		'''
		stack2matchhdr=EMData(options.stack2match,0)
		targetApix = round(stack2matchhdr['apix_x'],4)
		targetBox = stack2matchhdr['nx']
		
		if options.stack2match:
			print "\n\n\n(e2spt_refprep.py) Options stack2match is", options.stack2match
			
			if options.apix:		
				cmd = 'e2fixheaderparam.py --input=' + str (options.stack2match) + ' --stem=apix --stemval=' + str( options.apix )
			
				p = subprocess.Popen( cmd, shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
				text=p.communicate()
				p.stdout.close()
			
			print "\n\nWill call ref match\n\n"
			print "before refmatch, size is", EMData(stack2processEd,0)['nx']
			refmatch(options, stack2processEd, stack2processSample)
			print "after refmatch, size is", EMData(stack2processEd,0)['nx']
			print "because targetbox is", targetBox
			
			refpostprocessing(options.stack2process, stack2processEd, options )
		
		else:
			print "(e2spt_refprep.py) You have NOT supplied stack2match !!!"
			
			if options.shrink > 1.0:
				targetApix = round(stack2processSample['apix_x'],4) * options.shrink
				targetBox = stack3processSample['nx'] / options.shrink
				print "There's shrink >1.0"
				
			elif options.shrink:
				print """ERROR: If supplying --shrink, it needs to be a decimal number
				larger than 1.0"""
				sys.exit()
		
			print "test apix"
			if options.apix:
				print "(e2spt_refprep.py) There's options.apix", options.apix
				targetApix = float( options.apix )
			
			print "Test bpx"
			if options.boxsize:
				print "There's boxsize"
				targetBox = options.boxsize	
			
			print "test preproc opts"
			if not options.lowpass and not options.highpass and not options.shrink and not options.preprocess and not options.mask and not options.normproc and not options.threshold and not options.mirror and not options.sym and not options.boxsize and not options.apix:
				print """(e2spt_refprep.py) ERROR: You must supply at least one of the following preprocessing parameters,
				if you don't supply a stack to match via --stack2match:
				--lowpass, --highpass, --shrink, --preprocess, --mask, --threshold, --normproc, --mirror, --sym"""
				sys.exit()
			
			else:
				print "You have supplied at least one processing parameter and therefore you'll get an output file"
				refpostprocessing(stack2processSample, stack2processEd, options )
			
				if options.shrink:
					preciseShrink(options, stack2processSample, stack2processEd, targetApix, targetBox)
			
			print "Will proceed"
			
		print "The file whose header needs fixing is", stack2processEd
			
		cmd = 'e2fixheaderparam.py --input=' + stack2processEd + ' --stem=origin --stemval=0'
		print "The first command to fix the header is", cmd
		
		p = subprocess.Popen( cmd, shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		text=p.communicate()
		p.stdout.close()
		
		hdr = EMData( stack2processEd, 0, True )
		nx = hdr['nx']
		ny = hdr['ny']
		nz = hdr['nz']
		
		cmd2 = 'e2fixheaderparam.py --input=' + stack2processEd + ' --params=MRC.mx:' + str(targetBox) + ',MRC.my:' + str(targetBox) + ',MRC.mz:' + str(targetBox) + ',MRC.nx:' + str(targetBox) + ',MRC.ny:' + str(targetBox) + ',MRC.nz:' + str(targetBox)
		cmd2 += ',MRC.xlen:' + str(targetBox) + ',MRC.ylen:' + str(targetBox) + ',MRC.zlen:' + str(targetBox) + ',MRC.nxstart:0,MRC.nystart:0,MRC.nzstart:0'	
		print "The second command to fix the header is", cmd2
		
		p = subprocess.Popen( cmd2, shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		text=p.communicate()
		p.stdout.close()
		
	E2end(logger)
		
	return


def refmatch( options , stack2processEd, stack2processSample):
	'''
	Get the target boxsize for the reference. If options.boxsize is supplied, first change the boxsize of the data the reference is supposed to match.
	'''	
	#if options.refStack:

	print "\nInside refmatch"
	targetApix = round(EMData(options.stack2match,0,True)['apix_x'],4)
	
	if options.apix:
		targetApix = float( options.apix )

	targetBox = EMData(options.stack2match,0,True)['nx']
	
	print "\n\nTargetBox is", targetBox
	cmd = ''
	if options.boxsize:
		if int( options.boxsize ) != int( targetBox ):
			refStackEd = options.stack2match.replace('.hdf','_ed.hdf')
		
			print "(e2spt_refprep.py) refStackEd is", refStackEd
		
			cmd ='e2proc3d.py ' + str(options.stack2match) + ' ' + str(refStackEd) + ' --clip=' + str( options.boxsize )
			p=subprocess.Popen( cmd, shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			text=p.communicate()	
			p.stdout.close()
	
			options.stack2match = refStackEd

		targetBox = options.boxsize

	#if not options.stack2match and not options.boxsize:
	#	print """(e2spt_refprep.py) ERROR: You must supply an image through --refStack to get the target boxsize for the edited stack2process, or specify this number through --boxsize. If all you want to do
	#	is scale the data in refStack or stack2process by an integer amount, you ought to use e2proc3d.py <input> <output> --process=math.meanshrink:n="""
	#	sys.exit()

	cmd = 'e2fixheaderparam.py --input=' + stack2processEd + ' --stem=origin --stemval=0 && e2fixheaderparam.py --input=' + stack2processEd + " --params=MRC.nxstart:0,MRC.nystart:0,MRC.nzstart:0"
	p=subprocess.Popen( cmd, shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	text=p.communicate()	
	p.stdout.close()
	#print "\n(e2spt_refprep.py) Feedback from fix original is", text

	#stack2processEdHdr = EMData(stack2processEd,0,True)
	#stack2processEdApix = float( stack2processEdHdr['apix_x'])

	#refStackHdr = EMData(options.stack2match,0,True)
	#refStackApix = float( refStackHdr['apix_x'])
	#print "\n(e2spt_refprep.py) stack2processEdApix and refStackApix are", stack2processEdApix,refStackApix
	#if stack2processEdApix != refStackApix:
	#	cmd = 'e2fixheaderparam.py --input=' + stack2processEd + ' --stem=apix --stemval=' + str( targetApix ) 
	#	p=subprocess.Popen( cmd, shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	#	text=p.communicate()	
	#	p.stdout.close()
	#	print "Feedback from fix apix is", text


	#ref.process_inplace("xform.scale",{"scale":scalefactor,"clip":target_boxsize})

	preciseShrink( options, stack2processSample, stack2processEd, targetApix, targetBox)
	if options.verbose:
		print "\n(e2spt_refprep.py) The stack2process %shas now been precisely shrunk" %(stack2processEd)
	
	return
	

def preciseShrink( options, stack2processSample, stack2processEd, targetApix, targetBox ):
	'''
	Calculate the scale between the reference model and the data, round to the nearest integer to be able to use math.meanshrink (which preserves features),
	then calculate the scale factor again and use xform.scale, to get the exact scaling, even if features might become a bit distorted by xform.scale
	'''	
	
	print "inside preciseshrink targetbox is", targetBox
	
	stack2processApix = round(float( stack2processSample['apix_x'] ),4)
	if options.verbose:
		print "\n(e2spt_refprep.py) I've read the apix of the particles in stack2process, which is", stack2processApix

	meanshrinkfactor = float( targetApix )/float(stack2processApix)
	#meanshrinkfactor_int = int(round(meanshrinkfactor))

	meanshrinkfactor_int = int(meanshrinkfactor)

	if options.verbose:
		print "\n(e2spt_refprep.py) The refStack apix is", round(EMData( options.stack2match, 0, True)['apix_x'],4)
		print "(e2spt_refprep.py) And the target apix is", targetApix
		print "(e2spt_refprep.py) Therefore, the meanshrink factor is", meanshrinkfactor
		print "(e2spt_refprep.py) Which, for the first step of shrinking (using math.meanshrink), will be rounded to", meanshrinkfactor_int

	cmd = ''
	if meanshrinkfactor_int > 1:
		print "\n\n\n\n(e2spt_refprep.py) About to MEAN shrink becuase meanshrink factor is", meanshrinkfactor
		print "(e2spt_refprep.py) The type of stack2process is", type( stack2processSample )
	
		cmd = 'e2proc3d.py ' + str(stack2processEd) + ' ' + str(stack2processEd) + ' --process=math.meanshrink:n=' + str(meanshrinkfactor_int)
		p=subprocess.Popen( cmd, shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		text=p.communicate()	
		p.stdout.close()
		#print "\n(e2spt_refprep.py) Feedback from SHRINK is", text
	
	
		#ref.process_inplace("math.meanshrink",{"n":meanshrinkfactor_int})
		if options.verbose:
			print "(e2spt_refprep.py) The stack2process was shrunk, to a first approximation, see", EMData(options.stack2process,0,True)['nx']
	
		stack2processApix = round(EMData( stack2processEd, 0 , True)['apix_x'],4)
	
	else:
		#targetbox = EMData( options.stack2match,0,True )['nx']
		
		cmd = 'e2proc3d.py ' + str(stack2processEd) + ' ' + str(stack2processEd) + ' --clip=' + str(targetBox)

		#cmd += ' && e2proc3d.py ' + str(stack2processEd) + ' ' + str(stack2processEd) + ' --process=math.meanshrink:n=' + str(meanshrinkfactor_int)
		p=subprocess.Popen( cmd, shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		text=p.communicate()	
		p.stdout.close()
		
	
	scalefactor = round(float( stack2processApix ),4)/round(float( targetApix),4)

	print "\n\n\n\n(e2spt_refprep.py) The finer scale factor to apply is", scalefactor
	
	print "Right before, apix is", round(EMData(stack2processEd,0,True)['apix_x'],4)
	print "(e2spt_refprep.py) The final clip box is", targetBox

	cmd = 'e2proc3d.py ' + str(stack2processEd) + ' ' + str(stack2processEd) + ' --process=xform.scale:scale=' + str(scalefactor) 
	
	'''
	Only need to add clipping to scaling if the box wasn't clipped before. Boxes need to be clipped first when you're "scaling up the data" (making it bigger) rather than shrinking it
	'''
	if int(meanshrinkfactor_int) > 1:
		cmd += ':clip=' + str(targetBox) + ' --apix=' + str(targetApix)
		cmd += ' && e2fixheaderparam.py --input=' + str(stack2processEd) + ' --stem apix --valtype float --stemval ' + str(targetApix)
		print "meanshrinkfactor_int > 1, it is", meanshrinkfactor_int
		print "target apix is", targetApix
		print "cmd is", cmd
	elif int(meanshrinkfactor_int) < 1:
		cmd += ' && e2proc3d.py ' + str(stack2processEd) + ' ' + str(stack2processEd) + ' --clip=' + str(targetBox)
		cmd += ' && e2fixheaderparam.py --input=' + str(stack2processEd) + ' --stem apix --valtype float --stemval ' + str(targetApix)
		print "meanshrinkfactor_int < 1, it is", meanshrinkfactor_int
		print "target apix is", targetApix
		print "cmd is", cmd
	
	p=subprocess.Popen( cmd, shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	text=p.communicate()	
	p.stdout.close()
	
		
	print "Right after, apix is", round(EMData(stack2processEd,0,True)['apix_x'],4)

	#print "(e2spt_refprep.py) Feedback from scale and clip is", text


	print "\n\n!!!!!!!!!!!!\n(e2spt_refprep.py) stack2porcessEd should have been clipped by now, let's see", EMData(stack2processEd,0,True)['nx']
	
	return


def refpostprocessing(stack2process, stack2processEd, options ):

	''' 
	Apply processing options. First, parse those that need to be parsed, then apply
	'''	
	if options.normproc: 
		options.normproc=parsemodopt(options.normproc)

	if options.mask: 
		options.mask=parsemodopt(options.mask)

	if options.preprocess: 
		options.preprocess=parsemodopt(options.preprocess)
	
	if options.lowpass: 
		options.lowpass=parsemodopt(options.lowpass)

	if options.highpass: 
		options.highpass=parsemodopt(options.highpass)
		
	if options.threshold: 
		options.highpass=parsemodopt(options.highpass)

	if options.verbose:
		print "I have parsed the processing options."
	
	
	nptcls = EMUtil.get_image_count(stack2processEd)
	
	
	for i in range(nptcls):
		
		ref = EMData(stack2processEd,i)
	
		'''
		Make the mask first (if specified), use it to normalize (optionally), then apply it.
		'''
		mask=EMData(ref["nx"],ref["ny"],ref["nz"])
		mask.to_one()
	
		if options.mask:
			mask.process_inplace(options.mask[0],options.mask[1])

		if options.normproc:
			if options.normproc[0]=="normalize.mask": 
				options.normproc[1]["mask"]=mask
		
			ref.process_inplace(options.normproc[0],options.normproc[1])
	
		ref.mult(mask)
	
		'''
		If normalizing, it's best to do normalize->mask->normalize>mask
		'''
		if options.normproc:		
			ref.process_inplace(options.normproc[0],options.normproc[1])
			ref.mult(mask)
	
		'''
		Apply any specified filters (preprocess, lowpass and/or highpass)
		'''
		if options.preprocess:
			ref.process_inplace(options.preprocess[0],options.preprocess[1])
		
		if options.lowpass:
			ref.process_inplace(options.lowpass[0],options.lowpass[1])
		
		if options.highpass:
			ref.process_inplace(options.highpass[0],options.highpass[1])

		'''
		If a sharp filter for resolution wants to be applied, build it
		'''

		if options.sharpfiltres:
			stack2processHdr = EMData(stack2processEd,0,True)
			stack2processBox = stack2processHdr['nx']
			stack2processApix = round(stack2processHdr['apix_x'],4)
	
			resfac = 1.0/float(options.sharpfiltres)
			npixels = int(round(float( ref_box * ref_apix * res_fac )))

			actual_res = float(ref_box * ref_apix) / npixels
	
			if options.verbose:
				print "The sharp lowpass filter will be actually applied at this resolution", actual_res
				print "Becuase these many pixels in Fourier space will be zeroed out", stack2processBox/2 - npixels
	
			ref_table = [1.0] * npixels + [0.0] * (( stack2processBox/2) - npixels )
	
			ref.process_inplace("filter.radialtable",{"table":ref_table})

		ref['origin_x']=0
		ref['origin_y']=0
		ref['origin_z']=0
	
		ref['MRC.mx']= ref['nx']
		ref['MRC.my']= ref['ny']
		ref['MRC.mz']= ref['nz']
	
		ref['MRC.nxstart']= 0
		ref['MRC.nystart']= 0
		ref['MRC.nzstart']= 0
	
		ref['MRC.nx']= ref['nx']
		ref['MRC.ny']= ref['nx']
		ref['MRC.nz']= ref['nx']
		
		ref['MRC.xlen']= ref['nx']
		ref['MRC.ylen']= ref['ny']
		ref['MRC.zlen']= ref['nz']
	
		ref.write_image(stack2processEd,i)

		if options.sym and options.sym != 'c1':
			cmd = ''
			cmd = 'e2proc3d.py ' + str(stack2processEd) + ' ' + str(stack2processEd) + ' --sym=' + options.sym
			p=subprocess.Popen( cmd, shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			text=p.communicate()	
			p.stdout.close()

		if options.mirror:
	
			stack2processEdMirror = stack2processEd.replace('.hdf','_mirror.hdf')
	
			cmd = 'e2proc3d.py ' + stack2processEd + ' ' + stack2processEdMirror + ' --process=xform.mirror:axis=x'
			p=subprocess.Popen( cmd, shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			text=p.communicate()	
			p.stdout.close()

	return


if __name__ == '__main__':
	main()













