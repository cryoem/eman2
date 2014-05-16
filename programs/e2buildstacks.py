#!/usr/bin/env python
#
# Author: John Flanagan (jfflanag@bcm.edu)
# Copyright (c) 2000-2011 Baylor College of Medicine


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
import os, re
from EMAN2 import *

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options] --stackname myfile.hdf <image1> <image2> <image3> ...
	This program will combine many image files into a single output file. 
	
	If the output name has a ".lst" extension:
	the output is a formatted text file, one line per image, describing the file containing the actual
	image data in a searchable form. .lst files can be used as if they conatined actual images in any
	EMAN2 programs.
	
	If the output is a normal image file (.hdf, .spi, etc.) then the images will be copied into the
	output file sequentially in the order provided on the command-line. Some file formats will not
	support multiple images, or multiple volumes. Appropriate errors will be raised in these cases.
	HDF is the only format supporting full metadata retention for stacks of images or volumes.
	
	The output file will be emptied and overwritten!
	"""
	
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	
	parser.add_pos_argument(name="stack_files",help="List building material (sets) here.", default="", guitype='filebox', browser="EMParticlesEditTable(withmodal=True,multiselect=True)",  row=0, col=0,rowspan=1, colspan=2, nosharedb=True)
	parser.add_argument("--stackname",type=str,help="Name of the stack to build", default=None, guitype='strbox',row=2, col=0, rowspan=1, colspan=1)
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, help="verbose level [0-9], higner number means higher level of verboseness",default=1)
	
	(options, args) = parser.parse_args()
	
	if options.stackname==None :
		print "--stackname is required (output file)"
		sys.exit(1)
	
	# remove existing output file
	if os.path.exists(options.stackname) :
		try: os.unlink(options.stackname)
		except:
			print "ERROR: Unable to remove ",options.stackname,". Cannot proceed"
			sys.exit(1)
			
	# if output is LSX format, we handle it differently, with a specific object for these files
	if options.stackname[-4:].lower()==".lst" :
		outfile=LSXFile(options.stackname)
	else: outfile=None
	
	n=0		# number of images in output file
	for infile in args:
		nimg = EMUtil.get_image_count(infile)		# number of images in each input file as it is processed
		
		if options.verbose : 
			if nimg==1 : print infile
			else : print infile,nimg

		for i in xrange(nimg):
			if outfile!=None:
				outfile.write(n,i,infile)
			else:
				img=EMData(infile,i)
				img.write_image(options.stackname,n)
			n+=1
			
	if options.verbose : print n," total images written to ",options.stackname
			
			
if __name__ == "__main__":
	main()
