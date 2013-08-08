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
	Prepares one volume (the "reference") to be at the same apix and boxsize than another (a "particle" you might want to align to the reference).
	Takes PDB, MRC and HDF file formats, but outputs HDF or MRC only. The reference can also be subject to symmetrization, filters, and other processing options. 
	"""
	
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument('--stack2process', type=str, default='', help="Path to the stack that needs to be processed to match stack2.")
	parser.add_argument('--refStack', type=str, default='', help="Path to the image or stack of images which stack2process will match after processing.")
	parser.add_argument('--output', type=str, default='', help="Specify the name of the file to which the edited stack2process will be written.")
	parser.add_argument("--boxsize",type=int,help="If NOT specified, the reference's box size will match that of the data. If specified, both the reference and the data will be resized.",default=0)
	parser.add_argument("--sym", type=str, default='', help='Will apply the specified symmetry to the edited stack2process.')
	parser.add_argument("--apix",type=float,default=None,help="If specified, the program will assume this is the correct apix for the DATA, which will be written to the data's header if the data is in .hdf format, and will also be used to scale the edited stack2.")
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	parser.add_argument("--preprocess",type=str,help="Any processor (as in e2proc3d.py) to be applied to the edited stack2process.", default='')
	parser.add_argument("--lowpass",type=str,help="A lowpass filtering processor (as in e2proc3d.py) be applied to the edited stack2process.", default='')
	parser.add_argument("--highpass",type=str,help="A highpass filtering processor (as in e2proc3d.py) to be applied to the edited stack2process.", default='')
	parser.add_argument("--mask",type=str,help="Mask processor applied to particles before alignment. Default is mask.sharp:outer_radius=-2", default="mask.sharp:outer_radius=-2")
	parser.add_argument("--normproc",type=str,help="Normalization processor applied to particles before alignment. Default is to use normalize.mask. If normalize.mask is used, results of the mask option will be passed in automatically. If you want to turn this option off specify \'None\'", default="normalize.mask")
	parser.add_argument("--mirror", action="store_true", help="Will generate a mirrored copy of the edited stack2process.", default=False)
	parser.add_argument("--sharplowpassresolution",type=float,default=None,help="If specified, the edited stack2process will be sharply filtered to this resolution.")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)

	(options, args) = parser.parse_args()

	logger = E2init(sys.argv, options.ppid)
	
	'''
	Check for sane formats for stack2process and output
	'''
	
	if '.hdf' not in options.stack2process and '.pdb' not in options.stack2process and '.mrc' not in options.stack2process:
		print "ERROR: The format of the reference to edit must be .pdb, .mrc or .hdf"
		sys.exit()
	
	if '.hdf' not in options.output and '.mrc' not in options.output:
		print "ERROR: The output reference must be written to either an mrc or an hdf file."
		sys.exit() 
	
	check=0
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
			print "I've converted the .pdf reference to .mrc"
	
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
			print "I've converted the .mrc reference to .hdf"
	
	
	nStack2process = EMUtil.get_image_count( options.stack2process )
	
	if options.verbose:
		print "There are these many particles in stack2process to be edited",nStack2process
	
	'''
	Make all sides of the reference box equally sized and even, if they're not
	'''
	
	stack2processEd = options.stack2process.replace('.hdf','_ed.hdf')
	if options.output:
		stack2processEd = options.output
	
	
	print "stack2processEd is", stack2processEd
	
	stack2processSample = EMData( options.stack2process, 0, True)
	
	cmd = 'e2proc3d.py ' + options.stack2process + ' ' + stack2processEd + ' && e2fixheaderparam.py --input=' + stack2processEd + ' --stem=origin --stemval=0' 
	
	print "This command will copy and create it", cmd
	
	p=subprocess.Popen( cmd, shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	text=p.communicate()	
	p.stdout.close()
	
	print "Command executed. Feedback was", text
	
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
			print "I have clipped the particles in stack2process to have a cubical and even-sided box size."

	'''
	If options.apix is supplied, first change the apix of refStack to that; then make stack2process match it
	'''
	targetApix=0.0
	cmd = ''
	
	print "Options refstack is", options.refStack


	if options.refStack:
		targetApix = EMData(options.refStack,0,True)['apix_x']
		
		if options.apix:		
			cmd = 'e2fixheaderparam.py --input=' + str (options.refStack) + ' --stem=apix --stemval=' + str( options.apix )
			
			p = subprocess.Popen( cmd, shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			text=p.communicate()
			p.stdout.close()
			
			targetApix = float( options.apix )
			
	else:
		print "You have NOT supplied refStack see", options.refStack
		if not options.apix:
			print """ERROR: You must supply an image through --refStack to get the target apix for the edited stack2process, or specify this number through --apix. If all you want to do
			is clip the box of stack2process, you ought to use e2proc3d.py <input> <output> --clip=newboxsize. If all you want to do is change the apix value on the header, you ought to
			use e2fixheaderparam.py --input=stack.hdf --stem=apix --stemval=whatever_apix_value."""
			sys.exit()
		else:
			print "There'e options.apix", options.apix
			targetApix = float( options.apix )
	
	'''
	Get the target boxsize for the reference. If options.boxsize is supplied, first change the boxsize of the data the reference is supposed to match.
	'''	
	#if options.refStack:
	
	targetBox = EMData(options.refStack,0,True)['nx']
	print "\n\nTargetBox is", targetBox
	cmd = ''
	if options.boxsize:
		if int( options.boxsize ) != int( targetBox ):
			refStackEd = options.refStack.replace('.hdf','_ed.hdf')
			
			print "refStackEd is", refStackEd
			
			cmd ='e2proc3d.py ' + str(options.refStack) + ' ' + str(refStackEd) + ' --clip=' + str( options.boxsize )
			p=subprocess.Popen( cmd, shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			text=p.communicate()	
			p.stdout.close()
		
			otpions.refStack = refStackEd

		targetBox = options.boxsize
	
	if not options.refStack and not options.boxsize:
		print """ERROR: You must supply an image through --refStack to get the target boxsize for the edited stack2process, or specify this number through --boxsize. If all you want to do
		is scale the data in refStack or stack2process by an integer amount, you ought to use e2proc3d.py <input> <output> --process=math.meanshrink:n="""
		sys.exit()
	
	'''
	Calculate the scale between the reference model and the data, round to the nearest integer to be able to use math.meanshrink (which preserves features),
	then calculate the scale factor again and use xform.scale, to get the exact scaling, even if features might become a bit distorted by xform.scale
	'''	
	
	stack2processApix = float( stack2processSample['apix_x'] )
	if options.verbose:
		print "I've read the apix of the particles in stack2process, which is", stack2processApix

	meanshrinkfactor = float( targetApix )/float(stack2processApix)
	meanshrinkfactor_int = int(round(meanshrinkfactor))
	
	if options.verbose:
		print "The refStack apix is", EMData( options.refStack, 0, True)['apix_x']
		print "And the target apix is", targetApix
		print "Therefore, the meanshrink factor is", meanshrinkfactor
		print "Which, for the first step of shrinking (using math.meanshrink), will be rounded to", meanshrinkfactor_int
	
	cmd = ''
	if meanshrinkfactor_int > 1:
		print "bout to shrink"
		print "The type of stack2process is", type( stack2processSample )
		
		cmd = 'e2proc3d.py ' + str(stack2processEd) + ' ' + str(stack2processEd) + ' --process=math.meanshrink:n=' + str(meanshrinkfactor_int)
		p=subprocess.Popen( cmd, shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		text=p.communicate()	
		p.stdout.close()
		print "\nFeedback from SHRINK is", text
		
		
		#ref.process_inplace("math.meanshrink",{"n":meanshrinkfactor_int})
		if options.verbose:
			print "The stack2process was shrunk, to a first approximation, see", EMData(options.stack2process,0,True)['nx']
		
		stack2processApix = EMData( stack2processEd, 0 , True)['apix_x']
	
	scalefactor = float( stack2processApix )/float( targetApix)
	
	
	
	#if options.verbose:
	
	print "The finer scale factor to apply is", scalefactor
	
	print "The final clip box is", targetBox
	cmd = 'e2proc3d.py ' + str(stack2processEd) + ' ' + str(stack2processEd) + ' --process=xform.scale:scale=' + str(scalefactor) + ':clip=' + str(targetBox)
	p=subprocess.Popen( cmd, shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	text=p.communicate()	
	p.stdout.close()
	print "Feedback from scale and clip is", text

	
	print "\n\n!!!!!!!!!!!!\nstack2porcessEd should have been clipped by now, let's see", EMData(stack2processEd,0,True)['nx']
	
	cmd = 'e2fixheaderparam.py --input=' + stack2processEd + ' --stem=origin --stemval=0 && e2fixheaderparam.py --input=' + stack2processEd + " --params=MRC.nxstart:0,MRC.nystart:0,MRC.nzstart:0"
	p=subprocess.Popen( cmd, shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	text=p.communicate()	
	p.stdout.close()
	print "Feedback from fix origina is", text

	stack2processEdHdr = EMData(stack2processEd,0,True)
	stack2processEdApix = float( stack2processEdHdr['apix_x'])
	
	refStackHdr = EMData(options.refStack,0,True)
	refStackApix = float( refStackHdr['apix_x'])
	print "stack2processEdApix and refStackApix are", stack2processEdApix,refStackApix
	if stack2processEdApix != refStackApix:
		cmd = 'e2fixheaderparam.py --input=' + stack2processEd + ' --stem=apix --stemval=' + str( targetApix ) 
		p=subprocess.Popen( cmd, shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		text=p.communicate()	
		p.stdout.close()
		print "Feedback from fix apix is", text

	
	#ref.process_inplace("xform.scale",{"scale":scalefactor,"clip":target_boxsize})

	if options.verbose:
		print "The stack2process has now been precisely shrunk."






	"""
	
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
	
	if options.verbose:
		print "I have parsed the processing options."

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
	if options.preprocess != None:
		ref.process_inplace(options.preprocess[0],options.preprocess[1])
			
	if options.lowpass != None:
		ref.process_inplace(options.lowpass[0],options.lowpass[1])
			
	if options.highpass != None:
		ref.process_inplace(options.highpass[0],options.highpass[1])
	"""
	
	'''
	If a sharp filter for resolution wants to be applied, build it
	'''
	
	if options.sharplowpassresolution:
		stack2processHdr = EMData(stack2process,0,True)
		stack2processBox = stack2processHdr['nx']
		stack2processApix = stack2processHdr['apix_x']
		
		resfac = 1.0/float(options.sharplowpassresolution)
		npixels = int(round(float( ref_box * ref_apix * res_fac )))

		actual_res = float(ref_box * ref_apix) / npixels
		
		if options.verbose:
			print "The sharp lowpass filter will be actually applied at this resolution", actual_res
			print "Becuase these many pixels in Fourier space will be zeroed out", stack2processBox/2 - npixels
		
		ref_table = [1.0] * npixels + [0.0] * (( stack2processBox/2) - npixels )
		
		
		for i in range(nStack2process):
			a=EMData( stack2processEd,i)
			a = a.process("filter.radialtable",{"table":ref_table})
			a['origin_x']=0
			a['origin_y']=0
			a['origin_z']=0
			a.write_image(stack2processEd,i)
	
	if options.sym and options.sym != 'c1':
		#refs=ref.process('xform.applysym',{'sym':options.refsym})
		#refsymname=options.output.replace('.','_' + options.refsym + '.')
		#refs.write_image(refsymname)
		cmd = ''
		cmd = 'e2proc3d.py ' + str(stack2processEd) + ' ' + str(stack2processEd) + ' --sym=' + options.sym
		p=subprocess.Popen( cmd, shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		text=p.communicate()	
		p.stdout.close()

	if options.mirror:
		#refm=ref.process('xform.mirror',{'axis':'x'})
		#refmname = options.output.replace('.','_mirror.')
		#refm.write_image(refmname)
		
		stack2processEdMirror = stack2processEd.replace('.hdf','_mirror.hdf')
		
		cmd = 'e2proc3d.py ' + stack2processEd + ' ' + stack2processEdMirror + ' --process=xform.mirror:axis=x'
		p=subprocess.Popen( cmd, shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		text=p.communicate()	
		p.stdout.close()


	E2end(logger)
		
	return

	
	

if __name__ == '__main__':
	main()













