#!/usr/bin/env python
# This program performs simple processing of .LST files

# Author: Steven Ludtke, 5/19/2016 (sludtke@bcm.edu)
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

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """Usage:\nproctxt.py [options] <txt 1> <txt 2> ... 
Simple manipulations of text files conatining multi-column data, as would be used with plotting programs."""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	####################
	parser.add_argument("--merge",type=str,help="Merge several files into a single output. All inputs must have the same number of rows. Row comments stripped.",default=None)
	parser.add_argument("--sortcomment",action="store_true",default=False,help="Sorts rows based on per-row comment (after #) before merging")
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, help="verbose level [0-9], higner number means higher level of verboseness",default=1)
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)


	(options, args) = parser.parse_args()
	
	if len(args)<1 : 
		parser.error("At least one lst file required")
		sys.exit(1)

	logid=E2init(sys.argv,options.ppid)


	if options.merge!=None:
		# read all files. data_sets is a list of files. each file is a list of rows. each row is a list of str values
		data_sets=[]
		for filename in args:
			fin=open(filename,"r")
			data=[]
			for line in fin:
				if "#" in line : 
					comment=line.split("#")[1].strip()
					line=line.split("#")[0].strip()
				if len(line)==0 : continue
				if "," in line : line=line.split(",")
				elif ";" in line : line=line.split(";")
				else : line=line.split()
				if options.sortcomment : line.insert(0,comment)
				data.append(line)
			
			if options.sortcomment:
				data.sort()
				data=[i[1:] for i in data]
				
			data_sets.append(data)
			
		# merge all of the columns into data_sets[0]
		for i in range(len(args)-1):
			if len(data_sets[i])!=len(data_sets[i+1]) :
				print "Error: {} has {} rows and {} has {}".format(args[i],len(data_sets[i]),args[i+1],len(data_sets[i]))
				sys.exit(1)
			
			for row in xrange(len(data_sets[i+1])): 
				data_sets[0][row].extend(data_sets[i+1][row])

		out=open(options.merge,"w")
		for row in data_sets[0]:
			out.write("\t".join(row))
			out.write("\n")
			

	E2end(logid)

if __name__ == "__main__":
	main()
