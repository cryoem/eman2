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

# e2help.py  07/23/3004  Steven Ludtke
# This program will provide a variety of EMAN2 help

from EMAN2 import *
from math import *
import os
import sys
from sys import exit


def main():
	progname = os.path.basename(sys.argv[0])
	helpstring =  """Help is available on the following topics:
processors, cmps, aligners, averagers, projectors, reconstructors, analyzers, symmetries, orientgens"""
	usage = """prog <topic>
	
Interactive help on a variety of topics."""
	usage += " "+helpstring

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	#parser.add_argument("--res", "-R", type=float, help="Resolution in A, equivalent to Gaussian lowpass with 1/e width at 1/res",default=2.8)
	#parser.add_argument("--box", "-B", type=str, help="Box size in pixels, <xyz> or <x>,<y>,<z>")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-2)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	
	(options, args) = parser.parse_args()
	
	if len(args)<1 : 
		print helpstring
		exit(0)
	
	l=None
	if args[0] in ("cmp","cmps") :
		print "Available comparitors:"
		l=dump_cmps_list()
	elif args[0] in ("analyzer","analyzers") :
		print "Available analysers:"
		l=dump_analyzers_list()
	elif args[0] in ("averager","averagers") :
		print "Available averagers:"
		l=dump_averagers_list()
	elif args[0] in ("processor","processors") :
		print "Available processors:"
		l=dump_processors_list()
	elif args[0] in ("projector","projectors") :
		print "Available projectors:"
		l=dump_projectors_list()
	elif args[0] in ("reconstructor","reconstructors") :
		print "Available reconstructors:"
		l=dump_reconstructors_list()
	elif args[0] in ("aligner","aligners") :
		print "Available aligners:"
		l=dump_aligners_list()
	elif args[0] in ("sym","symmetry","symmetries") :
		print "Available symmetries:"
		l=dump_symmetries_list()
	elif args[0] in ("orientgen","orientationgen","orientgens","orientationgens","orientationgenerators") :
		print "Available orientation generators:"
		l=dump_orientgens_list()
	elif args[0] in ("version"):
	   print EMANVERSION + ' (CVS' + CVSDATESTAMP[6:-2] +')' 
	else:
		print helpstring
		print "unknown option:",args[0]
		
	if l:
		if options.verbose>0:
			k=l.keys()
			k.sort()
			for i in k:
				print "%s : %s"%(i, l[i][0])
				for j in range(1,len(l[i]),3): 
					print "\t%s(%s) - %s"%(l[i][j],l[i][j+1],l[i][j+2])
		else :
			k=l.keys()
			maxk=max([len(ii) for ii in k])
			fmt="%%-%0ds : "%maxk
			k.sort()
			for i in k:
				print fmt%i,
				for j in range(1,len(l[i]),3): 
					print "%s(%s)  "%(l[i][j],l[i][j+1]),
				if len(k)>1: print ""

if __name__ == "__main__":
    main()
