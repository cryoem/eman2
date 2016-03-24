#!/usr/bin/env python

#
# Author: Steve Ludtke, 03/24/2016 (sludtke@bcm.edu)
# Copyright (c) 2000-2007 Baylor College of Medicine
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
from math import *
from os import remove
import sys

fprefs=[]	# these are persistent references for generating invariants. Specific to one box size, but we assume that we only deal with one size per run

def footprint(im):
	"""This routine is similar in function to the make_footprint method in concept, but with a specific implementation for this purpose.
	implemented in Python rather than C++ as a convenience for easy "tweaking"."""
	global fprefs
	
	if len(fprefs)==0:
		apix=im["apix_x"]
		

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options]

	This program is used to produce reference-free class averages from a population of mixed,
	unaligned particle images. Unlike e2refine2d, this program performs only a single pass of classification,
	based entirely on specific low-resolution invariants. Its primary aim is to identify groups of bad particles,
	or rather to separate them from the general population. This program may also tend to split particles by
	defocus if the range within a project is too large."""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	# we grab all relevant options from e2refine.py for consistency
	# and snag a bunch of related code from David

	parser.add_argument("--path",type=str,default=None,help="Path for the refinement, default=auto")
	parser.add_argument("--input", default="start.hdf",type=str, help="The name of the file containing the particle data", browser='EMSetsTable(withmodal=True,multiselect=False)', guitype='filebox', row=0, col=0, rowspan=1, colspan=3, mode="spr")
	parser.add_argument("--ncls", default=32, type=int, help="Number of classes to generate", guitype='intbox', row=1, col=0, rowspan=1, colspan=1, mode="spr")
	parser.add_argument("--threads", default=4,type=int,help="Number of threads to run in parallel on a single computer when multi-computer parallelism isn't useful", guitype='intbox', row=24, col=2, rowspan=1, colspan=1, mode="refinement")
#	parser.add_argument("--parallel","-P",type=str,help="Run in parallel, specify type:<option>=<value>:<option>:<value>",default=None, guitype='strbox', row=4, col=0, rowspan=1, colspan=3, mode="spr")
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)

	global options
	(options, args) = parser.parse_args()

	if options.path!=None and ("/" in options.path) :
		print "Path specifier should be the name of a subdirectory to use in the current directory. Neither '/' or '#' can be included. "
		sys.exit(1)
	if options.path == None:
		fls=[int(i[-2:]) for i in os.listdir(".") if i[:4]=="dnc_" and len(i)==6]
		if len(fls)==0 : fls=[0]
		options.path = "dnc_{:02d}".format(max(fls)+1)
		try: os.mkdir(options.path)
		except: pass



if __name__ == "__main__":
    main()
