#!/usr/bin/env python

#
# Author: David Woolford, 11/19/2007 (woolford@bcm.edu)
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
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  2111-1307 USA
#
#

from EMAN2 import *
from optparse import OptionParser
import os
import sys

#constants
HEADER_ONLY=True
HEADER_AND_DATA=False


def main():
	progname = os.path.basename(sys.argv[0])
	usage = """%prog [options] <input file> <tlt file> <output file>"""

	parser = OptionParser(usage=usage,version=EMANVERSION)
	(options, args) = parser.parse_args()
	if len(args)<3 : 
		print usage
		parser.error("Specify all input and output files")

	angles_filename=args[1]
	ali_filename=args[0]
	output_stack=args[2]

	logid=E2init(sys.argv)
	f=file(angles_filename,'r')
	lines=f.readlines()
	angles=[]
	for line in lines:
		angles.append(float(line))
		
	

	input_image=EMData()
	input_image.read_image(str(ali_filename),0,HEADER_ONLY)
	nx = input_image.get_attr('nx')
	ny = input_image.get_attr('ny')
	nz = input_image.get_attr('nz')
	
	print len(angles),nz
	if len(angles) != nz:
		print len(angles),nz,ny,nx
		raise RuntimeError("The number of angles in the tlt file does not match the number of images in the input image")

	print nx,ny,nz
	for z_index in range(0,nz):
		roi=Region(0,0,z_index,nx,ny,1)
		input_image=EMData()
		input_image.read_image(ali_filename,0, HEADER_AND_DATA, roi)
		print angles[z_index]
		input_image.set_attr("xform.projection",Transform({"type":"eman","az":90,"alt":angles[z_index],"phi":90}))
		#input_image.set_rotation(90,angles[z_index],90)
		input_image.set_attr('ptcl_repr',1)
		input_image.write_image(output_stack,-1)
		
	E2end(logid)
# If executed as a program
if __name__ == '__main__':
	main()	

