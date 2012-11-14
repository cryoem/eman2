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

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options] <Volume file>
	Prepares one volume (the "reference") to be at the same apix and boxsize than another (a "particle" you might want to align to the reference).
	Takes PDB, MRC and HDF file formats, but outputs HDF or MRC only. The reference can also be subject to symmetrization, filters, and other processing options. 
	"""
	
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument('--ref', type=str, default='', help="Path to the reference that needs to be processed to match the data.")
	parser.add_argument('--input', type=str, default='', help="Path to the image or stack of images which the reference will match after processing.")
	parser.add_argument('--output', type=str, default='', help="Specify the name of the file to which the edited reference will be written.")
	parser.add_argument("--boxsize",type=int,help="If NOT specified, the reference's box size will match that of the data. If specified, both the reference and the data will be resized.",default=None)
	parser.add_argument("--refsym", type=str, default='', help='Will apply the specified symmetry to the edited reference.')
	parser.add_argument("--apix",type=float,default=None,help="If specified, the program will assume this is the correct apix for the DATA, which will be written to the data's header if the data is in .hdf format, and will also be used to scale the reference.")
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	parser.add_argument("--preprocess",type=str,help="Any processor (as in e2proc3d.py) to be applied to the edited reference.", default=None)
	parser.add_argument("--lowpass",type=str,help="A lowpass filtering processor (as in e2proc3d.py) be applied to the edited reference.", default=None)
	parser.add_argument("--highpass",type=str,help="A highpass filtering processor (as in e2proc3d.py) to be applied to the edited reference.", default=None)
	parser.add_argument("--mask",type=str,help="Mask processor applied to particles before alignment. Default is mask.sharp:outer_radius=-2", returnNone=True, default="mask.sharp:outer_radius=-2", guitype='comboparambox', choicelist='re_filter_list(dump_processors_list(),\'mask\')', row=11, col=0, rowspan=1, colspan=3, mode='alignment,breaksym')
	parser.add_argument("--normproc",type=str,help="Normalization processor applied to particles before alignment. Default is to use normalize.mask. If normalize.mask is used, results of the mask option will be passed in automatically. If you want to turn this option off specify \'None\'", default="normalize.mask")
	parser.add_argument("--mirrorref", action="store_true", help="Will generate a mirrored copy of the reference.", default=False)
	parser.add_argument("--sharplowpassresolution",type=float,default=None,help="If specified, the reference will be sharply filtered to this resolution.")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)

	(options, args) = parser.parse_args()

	logger = E2init(sys.argv, options.ppid)
	
	'''
	Parse processing options, if supplied
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
	Check for sane formats for ref and output
	'''
	
	if '.hdf' not in options.ref and '.pdb' not in options.ref and '.mrc' not in options.ref:
		print "ERROR: The format of the reference to edit must be .pdb, .mrc or .hdf"
		sys.exit()
	
	if '.hdf' not in options.output and '.mrc' not in options.output:
		print "ERROR: The output reference must be written to either an mrc or an hdf file."
		sys.exit() 
	
	check=0
	if '.pdb' in options.ref:
		pdbmodel = options.ref
		#os.system('cp ' + pdbmodel + ' ' + options.path)
		pdbmodel = pdbmodel.split('/')[-1]
		mrcmodel = pdbmodel.replace('.pdb','.mrc')
		os.system('e2pdb2mrc.py ' + pdbmodel + ' ' + mrcmodel)
		options.ref = mrcmodel
		check=1
		if options.verbose:
			print "I've converted the .pdf reference to .mrc"
	
	if '.mrc' in options.ref and '.hdf' in options.output:
		mrcmodel = options.ref
		if check==0:
			#os.system('cp ' + mrcmodel + ' ' + options.path)
			mrcmodel = mrcmodel.split('/')[-1]
		hdfmodel = mrcmodel.replace('.mrc','.hdf')
		os.system('e2proc3d.py ' + mrcmodel + ' ' + hdfmodel + ' && rm ' + mrcmodel)
		options.ref = hdfmodel
		check=1
		if options.verbose:
			print "I've converted the .mrc reference to .hdf"
		
	ref = EMData(options.ref,0)
	
	if options.verbose:
		print "I have loaded the reference."
	
	'''
	Make all sides of the reference box equally sized and even, if they're not
	'''
	if ref['nx'] != ref['ny'] or ref['nx'] != ref['nz'] or ref['ny'] != ref['nz']:	
		x0 = ref['nx']
		y0 = ref['ny']
		z0 = ref['nz']
	
		side = max(x0,y0,z0)
		if side % 2:
			side += 1
		
		R = Region((x0 - side)/2, (y0 - side)/2, (z0 - side)/2, side, side, side)
		ref.clip_inplace(R)
		
		if options.verbose:
			print "I have clipped the reference to have a cubical and even-sided box size."

	ref_apix = ref['apix_x']
	if options.verbose:
		print "I've read the reference's apix, which is", ref_apix

	'''
	Get the target apix for the reference. If options.apix is supplied, first change the apix of the data the reference is supposed to match.
	'''	
	if options.input:
		target_apix = EMData(options.input,0,True)['apix_x']
	
	if options.apix:
		if options.input and '.hdf' in options.input:
			nptcls = EMUtil.get_image_count(options.input)
			for i in range(nptcls):
				a =EMData(options.input,i)
				a['apix_x'] = options.apix
				a['apix_y'] = options.apix
				a['apix_z'] = options.apix
				a.write_image(options.input.replace('.hdf','_ed.hdf'),i)
		
		target_apix = options.apix
	
	if not options.input and not options.apix:
		print """ERROR: You must supply an image through --input= to get the target apix for the edited reference, or specify this number through --apix= . If all you want to do
		is clip the box of the reference, you ought to use e2proc3d.py <input> <output> --clip=newboxsize"""
		sys.exit()
	
	'''
	Get the target boxsize for the reference. If options.boxsize is supplied, first change the boxsize of the data the reference is supposed to match.
	'''	
	if options.input:
		target_boxsize = EMData(options.input,0,True)['nx']
	
	if options.boxsize:
		if options.input:
			nptcls = 1
			if '.hdf' in options.input:
				nptcls = EMUtil.get_image_count(options.input)
				aux=0
				if options.apix:
					options.input=options.input.replace('.hdf','_ed.hdf')
					aux=1
					
			for i in range(nptcls):
				print "\n\nOPTIONSINPUT is!!!!!!!\n", options.input
				
				a = EMData(options.input,i)
				a.process_inplace('xform.scale',{'scale':1,'clip':options.boxsize})
				outname=options.input
				if aux==0:
					outname=options.input.replace('.hdf','_ed.hdf')
				a.write_image(outname,i)

		target_boxsize = options.boxsize
	
	if not options.input and not options.boxsize:
		print """ERROR: You must supply an image through --input= to get the target boxsize for the edited reference, or specify this number through --boxsize= . If all you want to do
		is scale the box of the reference, you ought to use e2proc3d.py <input> <output> --process=math.meanshrink:n="""
		sys.exit()
	
	'''
	Calculate the scale between the reference model and the data, round to the nearest integer to be able to use math.meanshrink (which preserves features),
	then calculate the scale factor again and use xform.scale, to get the exact scaling, even if features might become a bit distorted by xform.scale
	'''	
	meanshrinkfactor = target_apix/ref_apix
	meanshrinkfactor_int = int(round(meanshrinkfactor))
	
	if options.verbose:
		print "The ref's apix is", ref_apix
		print "And the target apix is", target_apix
		print "Therefore, the meanshrink factor is", meanshrinkfactor
		print "Which, for the first step of shrinking (using math.meanshrink), will be rounded to", meanshrinkfactor_int
	
	if meanshrinkfactor_int > 1:
		print "bout to shrink"
		print "The type of ref is", type(ref)
		ref.process_inplace("math.meanshrink",{"n":meanshrinkfactor_int})
		if options.verbose:
			print "The reference was shrunk, to a first approximation."
		
		ref_apix = ref['apix_x']
		
	scalefactor = ref_apix/target_apix
	
	if options.verbose:
		print "The finer scale factor to apply is", scalefactor
	
	ref.process_inplace("xform.scale",{"scale":scalefactor,"clip":target_boxsize})

	if options.verbose:
		print "The reference has now been precisely shrunk."

	''' 
	Apply processing options. Make the mask first (if specified), use it to normalize (optionally), then apply it. Then apply any specified filters.
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
	#If normalizing, it's best to do normalize->mask->normalize>mask
	'''
	if options.normproc:		
		ref.process_inplace(options.normproc[0],options.normproc[1])
		ref.mult(mask)
		
	'''
	#Preprocess, lowpass and/or highpass
	'''
	if options.preprocess != None:
		ref.process_inplace(options.preprocess[0],options.preprocess[1])
			
	if options.lowpass != None:
		ref.process_inplace(options.lowpass[0],options.lowpass[1])
			
	if options.highpass != None:
		ref.process_inplace(options.highpass[0],options.highpass[1])
	
	'''
	If a sharp filter for resolution wants to be applied, build it
	'''
	if options.sharplowpassresolution:
		ref_box = ref['nx']
		ref_apix = ref['apix_x']
		resfac = 1.0/float(options.sharplowpassresolution)
		npixels = int(round(float( ref_box * ref_apix * res_fac )))

		actual_res = float(ref_box * ref_apix) / npixels
		
		if options.verbose:
			print "The sharp lowpass filter will be actually applied at this resolution", actual_res
			print "Becuase these many pixels in Fourier space will be zeroed out", ref_box/2 - npixels
		
		ref_table = [1.0] * npixels + [0.0] * (( ref_box/2) - npixels )
		
		ref = ref.process("filter.radialtable",{"table":ref_table})

	ref.write_image(options.output,0)
	
	if options.refsym and options.refsym != 'c1':
		os.system( 'e2proc3d.py --sym=' + options.refsym + ' ' + options.output + ' ' + options.output )
	
	if options.mirrorref:
		ref = EMData(options.output)
	
		t = Transform({'type':'eman','mirror':True})
		refm = ref.transform(t)
		
		if '.hdf' in options.output:
			refmname = options.output.replace('.hdf','_mirror.hdf')
			refm.write_image(refmname)
	
		elif '.mrc' in options.output:
			refmname = options.output.replace('.mrc','_mirror.mrc')
			refm.write_image(refmname)

	E2end(logger)
		
	return()

if __name__ == '__main__':
	main()













