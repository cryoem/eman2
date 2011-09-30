#!/usr/bin/env python

#
# Author: David Woolford 04/16/2009 (woolford@bcm.edu)
# Copyright (c) 2000-2006 Baylor College of Medicine
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

from optparse import OptionParser
from EMAN2 import file_exists,EMData,E2init,E2progress,E2end,EMANVERSION,check_eman2_type_string,EMUtil,remove_file,get_file_tag
import EMAN2
from EMAN2db import EMTask
from e2tomoaverage import *

def check_options(options,args):
	'''
	check function for e2tomoallvall
	'''
	error = []
	error.extend(check_tomo_options(options))
	
	if len(args) < 2:
		error.append("Error - to average you must supply atleast two images")
	else:
		for arg in args:
			if not file_exists(arg):
				error.append("%s does not exist" %arg)
	
	if file_exists(options.output):
		if not options.force:
			error.append("Error - output file exists. Remove it or supply the force argument")
		else:
			remove_file(options.output)
	
	return error
	
class EMTomoAllVAll:
	'''
	This class breaks the jobs of the boot-strapped average generation procedure
	so that they can be run in parallel. 
	'''
	def __init__(self,files,options,logger):
		'''
		@param options - the options returned by the call to (options, args) = parser.parse_args() 
		@param args - a list of image names - args is that which is returned by the call to (options, args) = parser.parse_args()
		'''
		self.options = options
		self.files = files
		self.logger = logger
		self.using_cuda = EMUtil.cuda_available() and options.cuda
		
		self.output_images = None # these will be a list of images to write the output to
		self.init_output()
		
	def init_output(self):
		'''
		Called internally to initialize output images - they are stored in memory
		'''
		
		self.output_images = []
		e = EMData(len(self.files),len(self.files))
		e.to_zero()
		
		# need to store correlation score (or cmp score)
		self.output_images.append(e)
		
		# need to store translation x
		self.output_images.append(e.copy())
		# need to store translation y
		self.output_images.append(e.copy())
		# need to store translation z
		self.output_images.append(e.copy())
		
		# need to store translation az
		self.output_images.append(e.copy())
		# need to store translation alt
		self.output_images.append(e.copy())
		# need to store translation phi
		self.output_images.append(e.copy())
		
		
	def execute(self):
		'''
		Function to perform the main work
		'''
		alignment_jobs = []
		for i in range(len(self.files)):
			for j in range(i+1,len(self.files)):
				alignment_jobs.append([i,j])
		from e2tomoaverage import EMTomoAlignments
		alignments_manager = EMTomoAlignments(self.options)
		alignments_manager.execute(alignment_jobs, self.files,self)
				
		self.finalize_writing()

	def process_output(self,results):
		'''
		store the result of an alignment
		'''
		from e2tomoaverage import EMTomoOutputWriter
		output_writer = EMTomoOutputWriter() # note - has a relationship
		output_writer.process_output(results, self.output_images, self.files)

	def finalize_writing(self):
		'''
		Writes the output images to disk
		'''
		
		for i,image in enumerate(self.output_images):
			image.set_attr("file_names",self.files)
			image.write_image(self.options.output,i)
		


import os,sys
def main():
	progname = os.path.basename(sys.argv[0])
	usage = """%prog [options] <image1> <image2> <image3> <image4> ....

	WARNING: Experimental program. Contact jgmontoy@bcm.edu for more info.
	
	Boot straps an initial probe doing all versus all alignment of the input images

"""

	parser = OptionParser(usage=usage,version=EMANVERSION)

	parser.add_option("--align",type="string",help="The aligner and its parameters. e.g. --align=rt.3d.grid:ralt=180:dalt=10:dphi=10:rphi=180:search=5", default="rt.3d.grid")
	parser.add_option("--aligncmp",type="string",help="The comparator used for determing the best initial alignment", default="dot.tomo:threshold=0")
	parser.add_option("--cmp",type="string",help="The comparator used to obtain the final similarity", default="dot.tomo:threshold=0")
	parser.add_option("--ralign",type="string",help="This is the second stage aligner used to refine the first alignment. This is usually the \'refine\' aligner.", default=None)
	parser.add_option("--raligncmp",type="string",help="The comparator used for determing the refined alignment", default="dot.tomo:threshold=0")
	parser.add_option("--output",type="string",default="e2tomoallvall.hdf",help="The output image which will store the results matrix")
	parser.add_option("--force",action="store_true",default=False,help="Force overwrite the output file if it already exists")
	parser.add_option("--parallel",type="string",default=None,help="Use parallelism")
	parser.add_option("--nsoln", default=1, type="int",help="If supplied and greater than 1, the nsoln-best alignments will be written to a text file. This is useful for debug but may be left unspecified")
	parser.add_option("--shrink",type="int",help="Shrink the data as part of the alignment - for speed purposes but at the potential loss of accuracy",default=None)
	parser.add_option("--filter",type="string",help="The name and parameters of an EMAN2 processor. Will be applied prior to shrinking.",default=None)
	parser.add_option("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_option("--verbose", "-v", dest="verbose", action="store", metavar="n", type="int", default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	if EMUtil.cuda_available():
		parser.add_option("--cuda",action="store_true",help="GPU acceleration using CUDA. Experimental", default=False)
 
	(options, args) = parser.parse_args()
	
	error_messages = check_options(options,args)
	if len(error_messages) != 0:
		msg = "\n"
		for error in error_messages:
			msg += error +"\n"
		parser.error(msg)
		exit(1)
	
	logger=E2init(sys.argv,options.ppid)

	module = EMTomoAllVAll(args,options,logger)
	module.execute()
	
	E2end(logger)
	
	
# If executed as a program
if __name__ == '__main__':
	main() 
	
