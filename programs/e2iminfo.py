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

from optparse import OptionParser
import pprint
from EMAN2 import *

def get_data_type_string(datatype):
	dtstring = {
		0 : "UNKNOWN",
		1 : "CHAR",
		2 : "UNSIGNED CHAR",
		3 : "SHORT",
		4 : "UNSIGNED SHORT",
		5 : "INT",
		6 : "UNSIGNED INT",
		7 : "FLOAT",
		8 : "DOUBLE",
		9 : "SHORT_COMPLEX",
		10 : "USHORT_COMPLEX",
		11 : "FLOAT_COMPLEX"
	}
	return dtstring.get(datatype, 'UNKNOWN')

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """%prog [options] imagefile ...
This program will print out some information about the image.
	"""
	
	parser = OptionParser(usage=usage,version=EMANVERSION)
	
	parser.add_option("-H", "--header", action="store_true",help="show all header information",default=False)
	parser.add_option("-N", "--number", type="int", help="Image number for single image info",default=0)
	parser.add_option("-s", "--stat", action="store_true",help="Show statistical information about the image(s).",default=False)
	parser.add_option("-E", "--euler", action="store_true",help="Show Euler angles from header",default=False)
	parser.add_option("-a", "--all", action="store_true",help="Show info for all images in file",default=False)
	parser.add_option("--verbose", "-v", dest="verbose", action="store", metavar="n", type="int", default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	
	(options, args) = parser.parse_args()
	
	if len(args)<1:
		print usage
		parser.error("Specify image file")

	if options.header and options.stat : 
		print "Only one of --header and --stat may be specified"
		sys.exit(1)
	
	for imagefile in args:
	
		if imagefile.lower()[:4]!="bdb:" :
			nimg = EMUtil.get_image_count(imagefile)
			imgtype = EMUtil.get_image_type(imagefile)
			imgtypename = EMUtil.get_imagetype_name(imgtype)
			if imgtype == EMUtil.ImageType.IMAGE_SPIDER and not options.stat: 
				image_index = -1
				
			d=EMData()
			d.read_image(imagefile, options.number, True)
			if d["nz"]==1 : print "%s\t %d images in %s format\t%d x %d"%(imagefile,nimg,imgtypename,d["nx"],d["ny"])
			else : print "%s\t %d images in %s format\t%d x %d x %d"%(imagefile,nimg,imgtypename,d["nx"],d["ny"],d["nz"])
			
			
		else :
			dct=db_open_dict(imagefile,ro=True)
			d=dct.get_header(0)
			nimg=len(dct)
			
			d=EMData(imagefile, options.number, True)
			if d["nz"]==1 : print "%s\t %d images in BDB format\t%d x %d"%(imagefile,len(dct),d["nx"],d["ny"])
			else : print "%s\t %d images in BDB format\t%d x %d x %d"%(imagefile,len(dct),d["nx"],d["ny"],d["nz"])

		if options.all : imgn=xrange(nimg)
		else : imgn=[options.number]

		if options.stat :
			for i in imgn:
				d=EMData(imagefile,i)
				print "%d.min=%1.4g\tmax=%1.4g\tmean=%1.4g\tsigma=%1.4g\tskewness=%1.4g \tkurtosis=%1.4g"%(i,d["minimum"],d["maximum"],d["mean"],d["sigma"],d["skewness"],d["kurtosis"]),
				try:
					c=d["ctf"]
					print "\tdefocus=%1.2f\tB=%1.0f"%(c.defocus,c.bfactor)
				except:
					print " "
		if options.euler:
			for i in imgn:
				d=EMData(imagefile,i)
				print "%d. %s"%(i,str(d["xform.projection"]))

		if options.header :
			for i in imgn:
				d=EMData(imagefile,i)
				print "%d."%i,
				keys=d.get_attr_dict().keys()
				keys.sort()
				for k in keys:
					print "\t%s: %s"%(k,str(d[k]))
				print
			
if __name__ == "__main__":
	main()
