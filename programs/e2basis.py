#!/usr/bin/env python

#
# Author: Steven Ludtke, 1/21/2008 (sludtke@bcm.edu)
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
# e2basis.py  01/21/2008  Steven Ludtke

from EMAN2 import *
from optparse import OptionParser
from math import *
import time
import os
import sys

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """%prog [options] <command> <file A> ...
	
Performs various options with basis sets, such as the orthogonal
basis produced by e2msa.py. 

varimax <basis input> <basis output>
	Performs a varimax rotation on an input basis set
	
project <basis input> <image input> <projection output>
	Projects a set of images into the input basis subspace"""

	parser = OptionParser(usage=usage,version=EMANVERSION)

	parser.add_option("--normproj",action="store_true",help="Normalize the projections resulting from 'project', such that the length of each vector is 1",default=False)
	
	#parser.add_option("--gui",action="store_true",help="Start the GUI for interactive boxing",default=False)
	#parser.add_option("--boxsize","-B",type="int",help="Box size in pixels",default=-1)
	#parser.add_option("--dbin","-D",type="string",help="Filename to read an existing box database from",default=None)
	
	(options, args) = parser.parse_args()
#	if len(args)>0 : parser.error("e2basis.py takes no arguments, only options")

	logid=E2init(sys.argv)
	
	basis=EMData.read_images(args[0])
	
	# Project an image stack into a basis subspace
	if args[0]=="project" :
		
	elif args[0]=="varimax" :
		
		
	
