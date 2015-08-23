#!/usr/bin/env python

#
# Author: Jesus Galaz, 28/March/2013. Updated: 25/September/2013
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

import os
from EMAN2 import *
		 
import sys
import numpy

import math	 
	 
def main():
	
	progname = os.path.basename(sys.argv[0])
	usage = """Must be run in the directory containing the stack(s)/particle(s) whose header 
	is to be modified. e2fixheaderparam.py --input=stack_to_fix --output=fixed_stack_name 
	--params=param1:value1,param2:value2... --type=type_of_parameters. 
	This programs fix values for any parameter on the header of an HDF file or a stack of 
	HDF files or an MRC file. ADDING new parameters, such as through --addfilename, or any
	parameters supplied through --params=parameter:value,parameter:value,.... that are not
	already in the header, will force HDF format on the output."""
			
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	
	parser.add_argument("--refheader",type=str,help="""If supplied, the header of this image
		will be copied to the header of all images in --input.""")
	
	parser.add_argument("--input", type=str, help="""File or stack for which to fix header 
		parameters. To indicate multiple files with a common string, use *.
		For example, if you want to process all the .mrc files in a given directory, 
		supply options.input=*.mrc. To process all files containing 'this' and 'that' in any part of the filename, supply --input=this*that.
		To process all files containing 'esto', 'eso' and 'aquello' in any part of the filename, supply --input=esto*eso*aquello.
		So on and so forth.""", default='')
	
	parser.add_argument("--output", type=str, help="""File to write the fixed stack to. 
		If not provided, the stack in --input will be overwritten.""", default='')
	
	parser.add_argument("--params", type=str, help="""Comma separated pairs of parameter:value. 
		The parameter will be changed to the value specified.""", default='')
	
	parser.add_argument("--stem", type=str, help="""Some parameters have common stems. 
		For example, 'origin_x', 'origin_y', 'origin"x'. 
		Supply the stem and all parameters containing it will be modified.""",default='')
	
	parser.add_argument("--stemval", type=str, help="""New value for all parameters 
		containing --stem.""",default='')
	
	parser.add_argument("--valtype", type=str, help="""Type of the value to enforce. 
		It can be: str, float, int, list, or transform.""",default='str')
	
	#parser.add_argument("--addparam", action='store_true', help="""If you want to add a new parameter to the header opposed to overwriting an existing one, turn this option on.""",default=False) 
	
	parser.add_argument("--addfilename", action='store_true', help="""Adds 
		the original filename of a file or stack to the header of each particle.
		This only works for .hdf files.""",default=False)
		
	parser.add_argument

	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n",type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness.")

	(options, args) = parser.parse_args()
	
	
	if not options.input:
		print "ERROR: You must supply an input image."
		sys.exit()
	
	t=[]
	tags=[]
	outputbase = os.path.basename( options.input ).split('.')[0]
	
	formats=['.hdf','.mrc','.mrcs','.st','.ali','.rec']
	
	originaloutput = options.output
	
	if options.output:
		outputbase = options.output
		
		print "\n(e2fixheader)(main)--output is", options.output
	
	files2process = []
	
	if '*' in options.input:
		ts = options.input.split('*')
		for t in ts:
			if t:
				tags.append( t )		
		
		current = os.getcwd()
		findir = os.listdir( current )
		
		ntags = len(tags)
	
		
		for f in findir:
			proceed=0
			for tag in tags:
				if tag in f:
					proceed+=1
		
			if '*' in options.input and int(proceed) == int(ntags):
				#options.input = f
				inputextension = f.split('.')[-1]
				if inputextension in f:
					files2process.append( f )
					print "\nFile appended!", f
				else:
					print "\nERROR: invalid image %s. You must supply .hdf,.mrc,.mrcs,.st,.ali or .rec files" %(f)
					sys.exit(1)
	else:
		inputextension = options.input.split('.')[-1]
		if inputextension in options.input:
			files2process.append( options.input )
		else:
			print "\nERROR: invalid image %s. You must supply .hdf,.mrc,.mrcs,.st,.ali or .rec files" %(options.input)
			sys.exit(1)
			
	k=0
	for fyle in files2process:
		if originaloutput:
			if len(files2process) > 1:
				options.output = outputbase.split('.')[0] + str(k).zfill( len(files2process) ) + '.' +  outputbase.split('.')[-1]
		else:
			#options.output = fyle.replace('.','_hdrEd.')
			options.output = fyle
			print "(e2fixheaderparam.py line 134) Defining options.output as", fyle
			
		if options.addfilename:
			
			if options.params:
				options.params += ',tag_originalfile:' + os.path.basename(fyle)
			else:
				options.params += 'tag_originalfile:' + os.path.basename(fyle)			
	
		
			newinput = fyle.split('.')[0] + '.hdf'
			os.system('e2proc3d.py ' + fyle + ' ' + newinput )
		
			options.output = options.output.split('.')[0] + '.hdf'
			
			fyle = newinput
			
			print "\n\n\n\n\n\n\n\n\nNEW input is", fyle
			print "\n\n\n\n\n\n\n\n\n"
		
		
		print "Sending this file for fixing", fyle
		fixer( fyle, options )
		k+=1
				
	return
		
	
def fixer(fyle, options):	
	formats=['.hdf','.mrc','.mrcs','.st','.ali','.rec']
	nonhdfformats = ['.mrc','.mrcs','.st','.ali','.rec']
	
	if options.output:
		#print "options output isa", options.output
		#print "lets see if .hdf is in it", '.hdf' in options.output[-4:0]
		#print "Therefore not in it should be false, and it is...", '.hdf' not in options.output[-4:0]
		#if '.hdf' not in options.output[-4:] and '.mrc' not in options.output[-4:] and 'bdb:' not in options.output[:5] and '.st' not in options.output[-4:] and '.ali' not in options.output[-4:] and '.rec' not in options.output[-4:]:
		
		print "\n(e2fixheader)(fixer) output is", options.output
		if options.output[-4:] in formats:
			pass
		elif options.output[-5:] in formats:
			pass
		elif options.output[-3:] in formats:
			pass
		else:	
			print "\n(e2fixheader)(fixer) ERROR: The output filename must be in .hdf or .mrc format (.mrcs, .st, .ali and .rec extensions are also allowed for MRC files)."
			sys.exit(1)
	

		
	aux1=aux2=aux3=0
	
	refheader = None
	if options.refheader:
		refheader = EMData( options.refheader, 0, True ).get_attr_dict()
		print "\n(e2fixheader)(fixer) reading --refheader from", options.refheader
	
	n = 1
	if fyle[-4:] == '.hdf' or fyle[-4:] == '.mrcs':
		n = EMUtil.get_image_count( fyle )
	
	for i in range(n):
		aux1=aux2=aux3=aux4=0
		
		indx=i
		if fyle[-4:] == '.hdf' or fyle[-4:] == '.mrcs':
			print "\n(e2fixheader)(fixer) fixing the header of img %d in stack %s" %( indx, fyle )
		else:
			indx=0
				
		imgHdr = EMData(fyle,indx,True)
		print "\n(e2fixheader)(fixer) type of imgHdr is", type(imgHdr)
		print "\n(e2fixheader)(fixer) and imgHdr is", imgHdr
		
		existingps = imgHdr.get_attr_dict()
		print "\n(e2fixheader)(fixer) existingps are", existingps
		
		aux1 = 0
		aux2 = 0
		
		if options.refheader:
			if refheader:
				img = EMData( fyle, indx )
				img.set_attr_dict( refheader )
			
				outfile = fyle
				if options.output:
					outfile = options.output 
			
				img.write_image( outfile, indx )
				print "overwriting header of img", indx
				aux4+=1
			else:
				print "ERROR: could not read --refheader", refheader
		
		
		elif not options.refheader:
		
			if options.params:
				paramValPairs=options.params.split(',')
			
				p2add = []
				print "!!!!!!!!!!!!!!!!!!!!\n!!!!!!!!!!!!!!!!!!!!Param pair vals are", paramValPairs
				for pair in paramValPairs:
					p = pair.split(':')[0]
					print "\nParsed parameter", p
					if p not in existingps:
						print "latter was not present already!"
						p2add.append( p )
			
				if len(p2add) > 0 and fyle[-4:] in nonhdfformats:
					tmp = options.input.split('.')[0] + '.hdf'
					cmd = 'e2proc3d.py ' + fyle + ' ' + tmp
					p=subprocess.Popen( cmd, shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
					text=p.communicate()
					p.stdout.close()
					fyle = tmp
				
					print "\nThere are parameters to add and the format is non-hdf. Therefore, this image will be created", fyle
				
					print """\nWARNING: You are trying to add parameters to a file format that does not allow this.
							You can only add new parameters to .hdf files. 
							A copy of the image will be created and saved as .hdf and the parameters will be added to it.
							The parameters that are not originally on the header of the image are:"""
				
					for p2a in p2add:
						print p2a
						
				for param in paramValPairs:
					p=param.split(':')[0]
					v=param.split(':')[-1]
					if options.valtype:
						v=valtyper(options,v)
					
					if p and p not in p2add:
						print "Will consider changing this existing param", p
						previousParam = imgHdr[p]
					
						if v != previousParam:
							aux1=1
							imgHdr[p] = v
							print "\nNew value %s for previous parameter %s" %(v,p)
				
					print "\n\nThe params to add are", p2add
					if p and p in p2add:
						print "Will add this new parameter", p
						aux2 = 1
						imgHdr.set_attr(p,v)
						print "\nNew PARAMETER %s added with this value %s" %(p,v)
						print "Let's see if we can read it", p, imgHdr[p]
						
				if not options.output:
				
					outputformat = fyle.split('.')[-1]
					if '.hdf' in fyle[-4:]:
						imgHdr.write_image(fyle,indx,EMUtil.ImageType.IMAGE_HDF,True)
					
					elif fyle[-4:] in nonhdfformats:
						imgHdr.write_image(fyle,-1,EMUtil.ImageType.IMAGE_MRC, True, None, EMUtil.EMDataType.EM_SHORT)
				
					else:
						print "ERROR: Only MRC (.mrc, .rec, .ali, .st) and HDF (.hdf) formats supported."
						sys.exit()
				
					print """\n\nOutput format will the be same as the input format (changed 
						to HDF by default if new parameters here added), which is""", outputformat
			
				else:
					outindx = 0
					if '.hdf' in options.output[-4:]:
						outindx=indx
					img = EMData(fyle,outindx)
					img.set_attr_dict(imgHdr.get_attr_dict())
					img.write_image(options.output,outindx)	
		
		
		
			if options.stem:
				try:
					v=options.stemval
				except:
					if not options.stemval:
						print "ERROR: If supplying --stem, you must also supply --stemval."
						sys.exit()	
			
			
				for param in existingps:
					#print "param being analyzed and its type are", param, type(param)
					#print "for stem", options.stem
					#print "param and stem are", param, options.stem
					if str(options.stem) in str(param):
						#print "Found stem in param!"
						if options.valtype:
							v=valtyper(options,v)
				
						imgHdr[param]=v
						aux3=1
				
				if aux3 !=0:
					if not options.output:
						outputformat = fyle.split('.')[-1]
					
						if '.hdf' in fyle[-4:]:
							imgHdr.write_image(fyle,indx,EMUtil.ImageType.IMAGE_HDF,True)				
						elif fyle[-4:] in nonhdfformats:
							imgHdr.write_image(fyle,-1,EMUtil.ImageType.IMAGE_MRC, True, None, EMUtil.EMDataType.EM_SHORT)
						else:
							print "ERROR: Only MRC (.mrc, .rec, .ali, .st) and HDF (.hdf) formats supported."
							sys.exit()
				
					else:
						img = EMData(fyle,indx)
						print "\nType of imgHdr is", type(imgHdr)
						print "\n\n\nand imgHdr is", imgHdr
						img.set_attr_dict(imgHdr.get_attr_dict())
						img.write_image(options.output,indx)	
				else:
					print "Couldn't find any parameters with the stem", options.stem
	
	if aux1 == n:
		print "\nformer parameter value(s) changed successfully!"
	if aux2 == n:
		print "\nnew parameter(s) added successfully!"	
	if aux3 == n:
		print "\nstem used successfully to change former parameter(s)!"
	if aux4 == n:
		print "\n--refheader copied onto --input"
	elif aux1 == 0 and aux2 == 0 and aux3 == 0 and aux4 == 0:
		print "\nno parameter(s) changed."	
	
	return	


def valtyper(options,v):
	if options.valtype == 'str':
		v=str(v)
	if options.valtype == 'int':
		v=int(v)
	if options.valtype == 'float':
		v=float(v)
	if options.valtype == 'list':
		v=list(v)
	return(v)

if __name__ == '__main__':
	main()
