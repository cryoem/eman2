#!/usr/bin/env python
# e2pdb2mrc.py  07/23/3004  Steven Ludtke
# This program will generate an electron density map from an MRC file.

from EMAN2 import *
from optparse import OptionParser
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
	usage = """%prog [options] input.cmm output.mrc
	"""

	parser = OptionParser(usage=usage,version=EMANVERSION)

	parser.add_option("--apix", "-A", type="float", help="A/voxel", default=1.0)
	parser.add_option("--res", "-R", type="float", help="Resolution in A, equivalent to Gaussian lowpass with 1/e width at 1/res",default=2.8)
	parser.add_option("--box", "-B", type="int", help="Box size in pixels",default=0)
	
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
