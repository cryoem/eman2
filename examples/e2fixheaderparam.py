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
	usage = """Must be run in the directory containing the stack(s)/particle(s) whose header is to be modified.
	e2fixheaderparam.py --input=stack_to_fix --output=fixed_stack_name --params=param1:value1,param2:value2... --type=type_of_parameters. 
	This programs fix values for any parameter on the header of an HDF file or a stack of HDF files."""
			
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	
	parser.add_argument("--input", type=str, help="""File or stack for which to fix header parameters. To indicate multiple files with a common string, use *.
													For example, if you want to process all the .mrc files in a given directory, supply options.input=*.mrc""", default='')
	parser.add_argument("--output", type=str, help="""File to write the fixed stack to. If not provided, the stack in --input will be overwritten.""", default='')
	parser.add_argument("--params", type=str, help="""Comma separated pairs of parameter:value. The parameter will be changed to the value specified.""", default='')
	parser.add_argument("--stem", type=str, help="""Some parameters have common stems. For example, 'origin_x', 'origin_y', 'origin"x'. 
												Supply the stem and all parameters containing it will be modified.""",default='')
	parser.add_argument("--stemval", type=str, help="""New value for all parameters containing --stem.""",default='')

	parser.add_argument("--valtype", type=str, help="""Type of the value to enforce. It can be: str, float, int, list, or transform.""",default='str')
	#parser.add_argument("--addparam", action='store_true', help="""If you want to add a new parameter to the header opposed to overwriting an existing one, turn this option on.""",default=False) 
	parser.add_argument("--addfilename", action='store_true', help="""Automatically adds the original filename of a file or stack to the header of each particle.
																	--params will be overwritten if this option is on.""",default=False) 

	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n",type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness.")

	(options, args) = parser.parse_args()
	
	if not options.input:
		print "ERROR: You must supply an input stack and it must be in HDF format."
		sys.exit()
	
	t=[]
	tags=[]
	outputbase = options.input.split('.')[0]
	if options.output:
		outputbase = options.output
		
		print "Output is", options.output
	
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
				files2process.append( f )
	else:
		files2process.append( options.input )
	
	k=0
	for fyle in files2process:
		if options.output:
			if len(files2process) > 1:
				options.output = outputbase.split('.')[0] + str(k).zfill( len(files2process) ) + '.' +  outputbase.split('.')[-1]
		else:
			options.output = fyle.replace('.','_hdrEd.')
			
		if options.addfilename:
			if options.params:
				options.params += ',originalfile:' + fyle.split('/')[0]
			else:
				options.params += 'originalfile:' + fyle.split('/')[0]			
	
		print "Sending this file for fixing", fyle
		fixer( options )
		k+=1
				
	return
		
	
def fixer(options):	
	formats=['.hdf','.mrc','.st','.ali','.rec']
	nonhdfformats = ['.mrc','.st','.ali','.rec']
	
	if options.output:
		#print "options output isa", options.output
		#print "lets see if .hdf is in it", '.hdf' in options.output[-4:0]
		#print "Therefore not in it should be false, and it is...", '.hdf' not in options.output[-4:0]
		#if '.hdf' not in options.output[-4:] and '.mrc' not in options.output[-4:] and 'bdb:' not in options.output[:5] and '.st' not in options.output[-4:] and '.ali' not in options.output[-4:] and '.rec' not in options.output[-4:]:
		if options.output[-4:] not in formats:	
			print "ERROR: The output filename must be in .hdf or .mrc format (.st, .ali and .rec format endings are also allowed for MRC files)."
			sys.exit()
	
	#if options.addparam:
	

			
		#if '.mrc' in options.output[-4:] or '.rec' in options.output[-4:] or '.st' in options.output[-4:] or '.ali' in options.output[-4:] or '.hdf' not in options.output[-4:]:
		#	print "ERROR: To add parameters to the header the OUTPUT format must be .hdf"
		#	sys.exit()
	
	
	
	if options.input[-4:] == '.hdf':
		n=EMUtil.get_image_count(options.input)
	else:
		n=1
		
	aux1=aux2=aux3=0
	for i in range(n):
		aux1=aux2=aux3=0
		
		indx=i
		if options.input[-4:] == '.hdf':
			print "Fixing the header of particle %d in the stack %s" %( indx, options.input )
		else:
			indx=0
				
		imgHdr = EMData(options.input,indx,True)
		print "\nType of imgHdr is", type(imgHdr)
		print "\n\n\nand imgHdr is", imgHdr
		
		existingps=imgHdr.get_attr_dict()
		print "\nExistingps are", existingps
		
		aux1 = 0
		aux2 = 0
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
			
			if len(p2add) > 0 and options.input[-4:] in nonhdfformats:
				tmp = options.input.split('.')[0] + '.hdf'
				cmd = 'e2proc3d.py ' + options.input + ' ' + tmp
				p=subprocess.Popen( cmd, shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
				text=p.communicate()
				p.stdout.close()
				options.input = tmp
				
				print "\nThere are parameters to add and the format is non-hdf. Therefore, this image will be create", options.input
				
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
						
			if not options.output:
				if '.hdf' in options.input[-4:]:
					imgHdr.write_image(options.input,indx,EMUtil.ImageType.IMAGE_HDF,True)
	
				elif options.input[-4:] in nonhdfformats:
					imgHdr.write_image(options.input,-1,EMUtil.ImageType.IMAGE_MRC, True, None, EMUtil.EMDataType.EM_SHORT)
				
				else:
					print "ERROR: Only MRC (.mrc, .rec, .ali, .st) and HDF (.hdf) formats supported."
					sys.exit()
			else:
				outindx = 0
				if '.hdf' in options.output[-4:]:
					outindx=indx
				img = EMData(options.input,outindx)
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
					if '.hdf' in options.input[-4:]:
						imgHdr.write_image(options.input,indx,EMUtil.ImageType.IMAGE_HDF,True)				
					elif options.input[-4:] in nonhdfformats:
						imgHdr.write_image(options.input,-1,EMUtil.ImageType.IMAGE_MRC, True, None, EMUtil.EMDataType.EM_SHORT)
				
				else:
					img = EMData(options.input,indx)
					print "\nType of imgHdr is", type(imgHdr)
					print "\n\n\nand imgHdr is", imgHdr
					img.set_attr_dict(imgHdr.get_attr_dict())
					img.write_image(options.output,indx)	
			else:
				print "Couldn't find any parameters with the stem", options.stem
	
	if aux1 !=0:
		print "Former parameter value changed successfully!"
	if aux2 !=0:
		print "A new parameter added successfully!"	
	if aux3 !=0:
		print "Stem used successfully to change former parameters!"
	elif aux1 == 0 and aux2 == 0 and aux3 == 0:
		print "No parameters seem to have changed."	
				
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
