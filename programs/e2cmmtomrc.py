#!/usr/bin/env python

#
# Author: Steven Ludtke, 04/10/2003 (sludtke@bcm.edu)
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
from math import *
import xml.sax
import os
import sys

from xml.sax.handler import ContentHandler

class myhandler(ContentHandler):
	def startDocument(self):
		self.parsed=[]
	
	def startElement(self,name,attrs):
		if name=="marker" :
			self.parsed.append((float(attrs["x"]),float(attrs["y"]),float(attrs["z"])))

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options] input.cmm output.mrc

	This program will read a 'marker file' produced by UCSF Chimera and turn it into an density map.
	"""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_argument("--apix", "-A", type=float, help="A/voxel", default=1.0)
	parser.add_argument("--res", "-R", type=float, help="Resolution in A, equivalent to Gaussian lowpass with 1/e width at 1/res",default=2.8)
	parser.add_argument("--box", "-B", type=int, help="Box size in pixels",default=0)
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	
	(options, args) = parser.parse_args()
	if len(args)<2 : parser.error("Input and output files required")
	try: chains=options.chains
	except: chains=None

	if options.box==0: 
		print("Box size required")
		sys.exit(1)
	
	if options.res<=options.apix : print "Warning: res<=apix. Generally res should be 2x apix or more"
	handler=myhandler()
	xml.sax.parse(args[0],handler)

	print "%d markers in CMM file"%len(handler.parsed)

	pa=PointArray()
	pa.set_number_points(len(handler.parsed))

	for i,j in enumerate(handler.parsed):
		pa.set_vector_at(i,Vec3f(j[0],j[1],j[2]),1.0)

	out=pa.pdb2mrc_by_summation(options.box,options.apix,options.res)
	out.write_image(args[1])
	
if __name__ == "__main__":
    main()
