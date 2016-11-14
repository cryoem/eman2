#!/usr/bin/env python
# This program performs simple processing of XML files

# Author: Stephen Murray, 10/17/14 (scmurray@bcm.edu)
# Copyright (c) 2014- Baylor College of Medicine
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

from EMAN2 import *
from math import *
import os
import sys
import xml.etree.ElementTree as ET

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """Usage:\nprocxml.py [options] <input file> <output file> \nSimple manipulations of XML files.\n."""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	####################
	parser.add_argument("--fscxml2txt",action="store_true",help="Convert the input XML FSC file into a text file")
	parser.add_argument("--fsctxt2xml",type=str,help="Convert the input text FSC file into an xml file. The value provided is the title of the graph. For example: --fsctxt2xml=\"My Graph Title\"", default="FSC Plot")
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, help="verbose level [0-9], higner number means higher level of verboseness",default=1)

	(options, args) = parser.parse_args()
	
	if len(args)!=2 : 
		parser.error("Input and Output files required")
		sys.exit(1)
	input_file, output_file = args[0],args[1]

	if options.fscxml2txt and options.fsctxt2xml != "FSC Plot":
		print "Invalid Options: Either --fscxml2txt or --fsctxt2xml can be used but not both"
		sys.exit(-1)

	if options.fsctxt2xml != "FSC Plot":
		#not ideal to do it by hand - should use the xml parser....
		input_file = open(str(input_file),'r')
		output_file = open(str(output_file),'w')
		output_file.write("<fsc title=\"" + options.fsctxt2xml+"\" xaxis=\"Resolution (A-1)\" yaxis=\"Correlation Coefficient\">\n")
		data = input_file.readlines()
		for i in range(len(data)):
			output_file.write("  <coordinate>\n    <x>"+str(data[i].split()[0])+"</x>\n    <y>"+str(data[i].split()[1])+"</y>\n  </coordinate>\n")
		output_file.write("</fsc>\n")
		input_file.close()
		output_file.close()

	if options.fscxml2txt:
		output_file = open(str(output_file),'w')
		element = ET.parse(input_file)
		fsc = element.getroot()
		for datapoint in fsc:
			for coordinate in datapoint:
				if coordinate.tag =="y":
					s = coordinate.text + "\n"
				else:
					s = coordinate.text + "\t"
				output_file.write(s)
		output_file.close()
		











































if __name__ == "__main__":
	main()
