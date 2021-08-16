#!/usr/bin/env python
#
# Author: Steven Ludtke, 03/06/2009 (sludtke@bcm.edu)
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

# e2parallel.py Steven Ludtke
# This program implements, via various options, the parallelism system for EMAN2

from past.utils import old_div
from future import standard_library
standard_library.install_aliases()
from builtins import range
from builtins import object
import traceback

import sys
from EMAN2 import *
from EMAN2PAR import EMMpiClient

debug=False
logid=None


def load_module(module):
	fname, cls=module.split('.')
	#print("import {} from {}".format(cls, fname))
	sys.path.append(os.path.join(e2getinstalldir(),"bin"))
	mod=__import__(fname, fromlist=[cls])
	setattr(sys.modules["__main__"], cls, getattr(mod,cls))

def main():
	
	usage=" "
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--scratchdir", type=str,help="Internal use only. Used by the MPI client",default="/tmp")
	parser.add_argument("--mode", type=str,help="choose from thread and mpi",default="thread")
	parser.add_argument("--taskin", type=str,help="Internal use only. Used when executing local threaded tasks.")
	parser.add_argument("--taskout", type=str,help="Internal use only. Used when executing local threaded tasks.")
	parser.add_argument("--loadmodule", type=str,help="load module",default="")
	parser.add_argument("--usethreads", type=int,help="max thread to use. only used for producing occupancy in mpi mode. default is the same as threads/mpi option given",default=-1)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	
	(options, args) = parser.parse_args()
	
	if options.mode=="mpi":
		client=EMMpiClient(options.scratchdir)
		if options.loadmodule!="":
			load_module(options.loadmodule)
		client.test(options.verbose)		# don't skip this. It's necessary to identify node names
		client.run(options.verbose, usethreads=options.usethreads)
		
	elif options.mode=="thread":
		from pickle import load,dump
		if options.loadmodule!="":
			load_module(options.loadmodule)
			

		task=load(open(options.taskin,"rb"))
		try: 
			ret=task.execute(empty_func)
			dump(ret,open(options.taskout,"wb"),-1)
		except:
			#### print to both stdout and text file
			print("Error executing: ",task)
			traceback.print_exc()
			err=traceback.format_exc()
			print(err)
			sys.exit(1)
		
	else:
		print("Cannot recognize mode. select from thread and mpi")

	return

def empty_func(val):
	return True



if __name__ == '__main__':
	main()
