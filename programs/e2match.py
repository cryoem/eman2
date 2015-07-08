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
	parser.add_argument('--img2process', type=str, default='', help="""Path to the stack 
		that needs to be processed to match img2match. If you want to process multiple stacks 
		or files, just separate them by commas --imgs2process=vo1.mrc,vol2.hdf,file.pdb""")
		
	parser.add_argument('--img2match', type=str, default='', help="Path to the image or stack of images which --img2process will match after processing. Not compulsory if --apix is provided.")
	parser.add_argument('--output', type=str, default='', help="Specify the name of the file to which the edited img2process will be written.")
	parser.add_argument("--boxsize",type=int,default=0,help="If NOT specified, the reference's box size will match that of the data. If specified, both the reference and the data will be resized.")
	parser.add_argument("--sym", type=str, default='', help='Will apply the specified symmetry to the edited img2process.')
	
	parser.add_argument("--apix",type=float,default=0.0,help="""If specified, the program 
		will assume this is the correct apix for --img2match if --img2match is provided,
		so the current value will be overwritten if --img2match is in .hdf format. 
		This value will also be used to scale the images in --img2process.""")
	
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	
	parser.add_argument("--preprocess",type=str,help="Any processor (as in e2proc3d.py) to be applied to the edited img2process.", default='')
	parser.add_argument("--lowpass",type=str,help="A lowpass filtering processor (as in e2proc3d.py) be applied to the edited img2process.", default='')
	parser.add_argument("--highpass",type=str,help="A highpass filtering processor (as in e2proc3d.py) to be applied to the edited img2process.", default='')
	parser.add_argument("--mask",type=str,help="Mask processor applied to particles before alignment. Default is None", default='')
	parser.add_argument("--normproc",type=str,default='',help="""Normalization processor applied to 
		particles before alignment. Default is None. If normalize.mask is used, 
		results of the mask option will be passed in automatically.""")

	parser.add_argument("--threshold",type=str,default='',help="""EMAN2 processor to be used
		to threshold the img2process. See available thresholds by typing 
		e2help.py processors --verbose=10
		at the command line.""")
		
	parser.add_argument("--shrink",type=float,default=0.0,help="""This will scale the img2process
		by the factor specified. This does NOT need to be an integer. You can enter 1.5, 2.7, etc.
		any decimal number > 1.0 is valid.""")	
	
	parser.add_argument("--mirror", action="store_true", help="Will generate a mirrored copy of the edited img2process.", default=False)
	parser.add_argument("--sharpfiltres",type=float,default=0.0,help="If specified, the edited img2process will be sharply filtered to this resolution.")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)

	(options, args) = parser.parse_args()

	logger = E2init(sys.argv, options.ppid)
	
	'''
	make sure apix values are consistently float and rounded to 4 significant digits
	'''
	if options.apix:
		options.apix = round(float(options.apix),4)
			
	'''
	if --apix and/or --boxsize are supplied, first change these values for --img2match if supplied; further down --img2process will be made to match them
	'''
	if options.img2match:
		'''
		check that --img2match has a valid format
		'''
		if '.mrc' not in options.img2match and '.hdf' not in options.img2match:
			print "\nERROR: --img2match must be in .hdf or .mrc format. Use e2proc3d.py or e2pdb2mrc.py to convert the file provided, %s, to a valid image format." %( options.img2match )
			sys.exit()
		else:
			
			img2matchhdr = EMData(options.img2match,0, True)
			img2matchhdrx = img2matchhdr['nx']
			img2matchhdry = img2matchhdr['ny']
			img2matchhdrz = img2matchhdr['nz']
			
			if options.boxsize or options.apix or img2matchhdrx!= img2matchhdry or img2matchhdrx != img2matchhdrz or img2matchhdry != img2matchhdrz:
				img2matchED = os.path.basename( options.img2match ).replace( '.','_ED.')
				cmd = 'cp ' + options.img2match + ' ' + img2matchED
				runcmd( options, cmd )
				options.img2match = img2matchED
			
			'''
			make all sides of --img2match equally sized and even, if they're not
			'''
			if img2matchhdrx!= img2matchhdry or img2matchhdrx != img2matchhdrz or img2matchhdry != img2matchhdrz:	
	
				newsize = max(img2matchhdrx,img2matchhdry,img2matchhdrz)
				if newsize % 2:
					newsize += 1
		
				cmd = 'e2proc3d.py ' + options.img2match + ' ' + options.img2match + ' --clip=' + str( newsize )
				runcmd( options, cmd )
			
			img2matchapix = round( float(img2matchhdr['apix_x']),4 )
		
			if options.boxsize:
				if int( options.boxsize ) != int( img2matchhdrx ) or int( options.boxsize ) != int( img2matchhdry ) or int( options.boxsize ) != int( img2matchhdrz ):
					cmd ='e2proc3d.py ' + options.img2match + ' ' + options.img2match + ' --clip=' + str( options.boxsize )
					runcmd( options, cmd )
				else:
					print "\nWARNING: the box side lengths of --img2match are already equal to --boxsize."
	
			if options.apix:
				if options.apix != img2matchapix:
					cmd = 'e2fixheaderparam.py --input=' + options.img2match + ' --stem=apix --stemval=' + str( options.apix ) + ' --valtype=float'
					runcmd( options, cmd )
				else:
					print "\nWARNING: the apix of --img2match is already equal to --apix."
		
			'''
			fix the origin values for --img2match
			'''		
			cmd = 'e2fixheaderparam.py --input=' + options.img2match + ' --stem=origin --stemval=0 && e2fixheaderparam.py --input=' + options.img2match + " --params=MRC.nxstart:0,MRC.nystart:0,MRC.nzstart:0"
			
				
	'''
	iterate through img2process files
	'''
	imgs2process = options.img2process.split(',')
	
	for img2p in imgs2process:
		
		img2p
		
		if len(imgs2process) > 1 or not options.output:
			options.output = img2p.split('.')[0] + '_PREP.hdf'
		
		if options.verbose:
			print "\n--output is", options.output

		if '.hdf' not in img2p and '.pdb' not in img2p and '.mrc' not in img2p:
			print "\n(e2match) ERROR: The format of --img2process %s must be .pdb, .mrc or .hdf" %( img2p )
			sys.exit()
	
		if '.hdf' not in options.output:
			print "\n(e2match) ERROR: --output must be in .hdf format."
			sys.exit() 
	
		current = os.getcwd()
		findir = os.listdir(current)
	
		if '.pdb' in img2p:
			pdbmodel = img2p
			#pdbmodel = pdbmodel.split('/')[-1]
			hdfmodel = os.path.basename( pdbmodel ).replace('.pdb','.hdf')
		
			cmd = 'e2pdb2mrc.py ' + pdbmodel + ' ' + hdfmodel
			runcmd( options, cmd )
		
			img2p = hdfmodel
			if options.verbose:
				print "\n(e2match) I've converted the .pdb reference to .hdf"
	
		elif '.mrc' in img2p:
			mrcmodel = img2p
			hdfmodel = os.path.basename( mrcmodel ).replace('.mrc','.hdf')
		
			cmd = 'e2proc3d.py ' + mrcmodel + ' ' + hdfmodel
			runcmd( options, cmd )
			
			img2p = hdfmodel
			
			if options.verbose:
				print "\n(e2match) I've converted the .mrc reference to .hdf"
		
		if '.hdf' in img2p:
			
			img2processED = os.path.basename( img2p ).replace( '.','_ED.')
			cmd = 'cp ' + img2p+ ' ' + img2processED
			runcmd( options, cmd )
			img2p = img2processED
			
			nStack2process = EMUtil.get_image_count( img2p )
	
			if options.verbose:
				print "\n(e2match) There are these many images in --img2process file %s" %( img2p )
	
			'''
			fix origin on header of --img2process current file img2p
			'''	
			cmd = ' e2fixheaderparam.py --input=' + img2p+ ' --stem=origin --stemval=0 && e2fixheaderparam.py --input=' + img2p + " --params=MRC.nxstart:0,MRC.nystart:0,MRC.nzstart:0" 
			runcmd( options, cmd )
		
			img2phdr = EMData( img2p, 0, True )
			
			
			'''
			make all sides of --img2process img2p equally sized and even, if they're not.
			In principle, you shouldn't ned to, but --process=scale:clip in e2proc3d crashes 
			if all the dimensions of an img are not equal
			'''
			if img2phdr['nx']!= img2phdr['ny'] or img2phdr['nx'] != img2phdr['nz'] or img2phdr['ny'] != img2phdr['nz']:	
	
				newsize = max(img2phdr['nx'],img2phdr['ny'],img2phdr['nz'])
				if newsize % 2:
					newsize += 1
		
				cmd = 'e2proc3d.py ' + img2p + ' ' + img2p+ ' --clip=' + str( newsize )
				runcmd( options, cmd )
					
				img2phdr = EMData( img2p, 0, True )
			
			targetApix = None
			targetBox = max(img2phdr['nx'],img2phdr['ny'],img2phdr['nz'])
			if targetBox % 2:
				targetBox += 1
			if options.boxsize:
				targetBox = int( options.boxsize )
			
			if options.img2match:
				img2matchhdr = EMData(options.img2match,0, True)
				targetBox = img2matchhdr['nx']
				targetApix = img2matchhdr['apix_x']	
			else:
				if options.apix:
					targetApix = round(float( options.apix ),4 )
				else:
					print """\nERROR: --apix required if img2match is not provided.
					If the apix of --img2process is already correct and you just want to
					match the boxsize, use e2proc3d.py input output --clip=targetBoxSize"""
					sys.exit()
				
			print "\n\ncalling preciseshrink function"
			print "before, size is", img2phdr['nx']		
						
			preciseshrink( options, img2p, targetApix, targetBox)
			
			img2processhdr = EMData( img2p, 0, True )
			print "after preciseshrink, size is", img2processhdr['nx']
			
			if options.normproc or options.mask or options.preprocess or options.lowpass or options.highpass or options.threshold or options.mirror or options.sym:			
				refpostprocessing( options, img2p )		
				
			cmd = 'e2fixheaderparam.py --input=' + img2p + ' --stem=origin --stemval=0 && e2fixheaderparam.py --input=' + img2p + " --params=MRC.nxstart:0,MRC.nystart:0,MRC.nzstart:0"
			if options.verbose:
				print "The first command to fix the header is", cmd
		
			hdr = EMData( img2p, 0, True )
			nx = hdr['nx']
			ny = hdr['ny']
			nz = hdr['nz']
		
			cmd2 = ' && e2fixheaderparam.py --input=' + img2p + ' --params=MRC.mx:' + str(targetBox) + ',MRC.my:' + str(targetBox) + ',MRC.mz:' + str(targetBox) + ',MRC.nx:' + str(targetBox) + ',MRC.ny:' + str(targetBox) + ',MRC.nz:' + str(targetBox)
			cmd2 += ',MRC.xlen:' + str(targetBox) + ',MRC.ylen:' + str(targetBox) + ',MRC.zlen:' + str(targetBox) + ',MRC.nxstart:0,MRC.nystart:0,MRC.nzstart:0'	
			if options.verbose:
				print "The second command to fix the header is", cmd2
		
			cmdf = cmd + cmd2
			runcmd( options, cmdf )
			
	E2end(logger)
		
	return
	

def runcmd(options,cmd):
	if options.verbose > 9:
		print "(e2spt_autoboxer.py)(runcmd) Running command", cmd
	
	p=subprocess.Popen( cmd, shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	text=p.communicate()	
	p.stdout.close()
	
	if options.verbose > 9:
		print "(e2spt_autoboxer.py)(runcmd) Done"
	return
	

def preciseshrink( options, img2processEd, targetApix, targetBox ):
	'''
	Calculate the scale between the reference model and the data, round to the nearest 
	integer to be able to use math.meanshrink (which preserves features), then calculate 
	the scale factor again and use xform.scale, to get the exact scaling, even if features 
	might become a bit distorted by xform.scale.
	
	Process img2proccessEd via the command line calling e2proc3d.py because it might be a stack,
	so we don't bother to iterate over the images in here.
	'''
	
	print "\n(e2match)(preciseshrink)"
	img2processEdhdr = EMData(img2processEd,0, True )
	img2processApix = round(float( img2processEdhdr['apix_x'] ),4)
	if options.verbose:
		print "\n(e2match) I've read the apix of the particles in img2process, which is", img2processApix

	meanshrinkfactor = float( targetApix )/float(img2processApix)
	#meanshrinkfactor_int = int(round(meanshrinkfactor))

	meanshrinkfactor_int = int(meanshrinkfactor)

	if options.verbose:
		print "\n(e2match) the refStack or --img2match apix is", round(EMData( options.img2match, 0, True)['apix_x'],4)
		print "(e2match) and the target apix is", targetApix
		print "(e2match) therefore, the meanshrink factor is", meanshrinkfactor
		print "(e2match) which, for the first step of shrinking (using math.meanshrink), will be rounded to", meanshrinkfactor_int

	cmd = ''
	
	if meanshrinkfactor_int > 1:
		print "\n\n(e2match) about to MEAN shrink becuase meanshrink factor is > 1", meanshrinkfactor
		print "(e2match) the type of img2process is", type( img2processEd )
	
		cmd = 'e2proc3d.py ' + img2processEd + ' ' + img2processEd + ' --process=math.meanshrink:n=' + str(meanshrinkfactor_int)
		runcmd( options, cmd )
	
		#ref.process_inplace("math.meanshrink",{"n":meanshrinkfactor_int})
		if options.verbose:
			print "(e2match) the img2process was shrunk, to a first approximation, see", EMData(options.img2process,0,True)['nx']
	
		img2processApix = round(EMData( img2processEd, 0 , True)['apix_x'],4)
	
	else:
		#targetbox = EMData( options.img2match,0,True )['nx']
		
		cmd = 'e2proc3d.py ' + img2processEd + ' ' + img2processEd + ' --clip=' + str(targetBox)
		runcmd( options, cmd )
	
	scalefactor = round(float( img2processApix ),4)/round(float( targetApix),4)

	print "\n\n\n(e2match) the finer scale factor to apply is", scalefactor
	
	print "right before, apix is", round(EMData(img2processEd,0,True)['apix_x'],4)
	print "(e2match) the final clip box is", targetBox

	cmd = 'e2proc3d.py ' + img2processEd + ' ' + img2processEd + ' --process=xform.scale:scale=' + str(scalefactor) 
	
	'''
	Only need to add clipping to scaling if the box wasn't clipped before. Boxes need to be clipped first when you're "scaling up the data" (making it bigger) rather than shrinking it
	'''
	if int(meanshrinkfactor_int) > 1:
		cmd += ':clip=' + str(targetBox) + ' --apix=' + str(targetApix)
		cmd += ' && e2fixheaderparam.py --input=' + img2processEd + ' --stem apix --valtype float --stemval ' + str(targetApix)
		print "meanshrinkfactor_int > 1, it is", meanshrinkfactor_int
		print "target apix is", targetApix
		print "cmd is", cmd
	elif int(meanshrinkfactor_int) < 1:
		cmd += ' && e2proc3d.py ' + img2processEd + ' ' + img2processEd + ' --clip=' + str(targetBox)
		cmd += ' && e2fixheaderparam.py --input=' + str(img2processEd) + ' --stem apix --valtype float --stemval ' + str(targetApix)
		print "meanshrinkfactor_int < 1, it is", meanshrinkfactor_int
		print "target apix is", targetApix
		print "cmd is", cmd
	
	runcmd( options, cmd )
	
		
	print "Right after, apix is", round(EMData(img2processEd,0,True)['apix_x'],4)

	#print "(e2match) Feedback from scale and clip is", text


	print "\n\n!!!!!!!!!!!!\n(e2match) img2porcessEd should have been clipped by now, let's see", EMData(img2processEd,0,True)['nx']
	
	return


def refpostprocessing( options, img2processEd ):

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
	
	
	nptcls = EMUtil.get_image_count(img2processEd)
	
	
	for i in range(nptcls):
		
		ref = EMData(img2processEd,i)
	
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
		
			#ref.process_inplace(options.normproc[0],options.normproc[1])
	
		'''
		If normalizing, it's best to do mask->normalize>mask
		'''
		
		ref.mult(mask)
		
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
			img2processHdr = EMData(img2processEd,0,True)
			img2processBox = img2processHdr['nx']
			img2processApix = round(img2processHdr['apix_x'],4)
	
			resfac = 1.0/float(options.sharpfiltres)
			npixels = int(round(float( ref_box * ref_apix * res_fac )))

			actual_res = float(ref_box * ref_apix) / npixels
	
			if options.verbose:
				print "The sharp lowpass filter will be actually applied at this resolution", actual_res
				print "Becuase these many pixels in Fourier space will be zeroed out", img2processBox/2 - npixels
	
			ref_table = [1.0] * npixels + [0.0] * (( img2processBox/2) - npixels )
	
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
	
		ref.write_image(img2processEd,i)

		if options.sym and options.sym != 'c1':
			cmd = 'e2proc3d.py ' + str(img2processEd) + ' ' + str(img2processEd) + ' --sym=' + options.sym
			p=subprocess.Popen( cmd, shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			text=p.communicate()	
			p.stdout.close()

		if options.mirror:
			img2processEdMirror = img2processEd.replace('.hdf','_mirror.hdf')
	
			cmd = 'e2proc3d.py ' + img2processEd + ' ' + img2processEdMirror + ' --process=xform.mirror:axis=z'
			p=subprocess.Popen( cmd, shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			text=p.communicate()	
			p.stdout.close()

	return


if __name__ == '__main__':
	main()













