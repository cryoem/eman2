#!/usr/bin/env python
from __future__ import print_function
from __future__ import division
#
# Author: Jesus Galaz, 06/05/2012 - Modifid 08/Jun/2017
# Copyright (c) 2011 Baylor College of Medicine
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

from future import standard_library
standard_library.install_aliases()
from builtins import range
import os
import sys 
import subprocess
from EMAN2 import *
from sys import argv

def main():
	progname = os.path.basename(sys.argv[0])
	
	usage = """Processes coordinate files (a text file with x y z coordinates) in various ways. May supply first argument instead. e2spt_coordsrefactor.py coords.txt. Otherwise, e2spt_coordsrefactor.py --input=coords.txt"""
		
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	
	parser.add_argument("--input", type=str, help="""Default=None. Text file containing the coordinates for SPT subvolumes of a tomogram, in a sane format X, Y, Z, or X, Z, Y, with ONE set of X, Y, Z coordinates per line, and NO EXTRA CHARACTERS, except blank spaces (a tab is fine too) in between coordinates.""")
	
	parser.add_argument("--mult", type=float,default=1.0,help="""Default=1.0. Factor to multiply the coordinates by. This can allow to expand or shrink the coordinates, but note that the resulting number will be rounded to the nearest integer.""")	

	parser.add_argument("--output", type=str, default='', help="Default=None. Output name for the refactored coordinates file.")
	
	parser.add_argument("--randomize", action="store_true",default=False, help="Default=false. Randomizes the coordinates so that they are in no preferential order.")
	parser.add_argument("--rotx", action="store_true", default=False, help="""Default=false. Will produce the effects of rotating the coordinates about the x axis. --tomothickness needs to be supplied if --rotx is on.""")

	parser.add_argument("--sort", type=int, default=0,help="""Will sort the coordinates in the file, by the order provided; it can sort by 1 coordinate, or 2, or 3); for example, 'z' will sort by z only, so that all the coordinates in the same 'slice' z-height, will be together in the file [assuming 'z' is the shortest dimension]; 'zx' would leave y unsorted; zxy will sort by z, and then at each z height, it will sort by x, then by y. You can provide ANY combination of sorting: xyz, xy, xz, yx, yzx... etc.""")
	parser.add_argument("--subset", type=int,default=0,help="--subset=n will select a subset of coordinate lines, from 1 to n, to write into the refactored file.")
	parser.add_argument("--subtract", type=str, default='', help="Default=None. Name of a coordinates .txt file to subtract from the --input coordinates file. Any coordinates found in --subtract that are also present in --input will be subtracted from --input.")
	parser.add_argument("--swapyz", action="store_true", help="This will swap the Y and Z coordinates.")

	parser.add_argument("--tomothickness", type=int, help="""Lentgth of the 'short' direction in the tomogram whence the cooridates came from. Must be in the same scale as the coordinates file supplied through --coords. This parameter is required if --rotx is used.""", default=0)

	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n",type=int, default=0, help="verbose level [0-9], higher number means higher level of verboseness")
	
	(options, args) = parser.parse_args()
	
	logger = E2init(sys.argv, options.ppid)

	cfile=''
	try:
		if len(argv) > 0:
			if '.txt' in argv[1]:
				cfile = argv[1]
			else:

				print("ERROR: coordinates file must be in .txt format")
				sys.exit(1)
	except:
		pas

	if options.input:
		cfile = options.input
	

	if options.rotx and not options.tomothickness:
		print("ERROR: You must supply --tomothickness if you intend to use --rotx")
		sys.exit()
		
	if not cfile:
		print("ERROR: Must provide a coordinates file")
		sys.exit(1)
	
	if options.randomize and options.sort:
		print("ERROR: Cannot randomize and sort at the same time; the functions are contradictory. Chooe one, please.")
	

	sanelines = loadlines(cfile)
	sanelines = fixlines(sanelines)
	
	if options.subtract:
		subtractlines = loadlines(options.subtract)
		subtractlines = fixlines(subtractlines)
		finalsetlines = [ line for line in sanelines if line not in subtractlines ]
		sanelines = finalsetlines
	
	n=len(sanelines)
	if options.verbose:
		print("\nyou have these many potentially sane lines in your coordinates file {} ".format(n))
	
	if options.subset and not options.randomize and not options.sort:
		n = options.subset

	#newlines = ['']*n
	newlines = []

	'''Some people might manually make ABERRANT coordinates files with commas, tabs, or more than once space in between coordinates'''
	for i in range(n):                             
				
		sanelines[i] = sanelines[i].split()

		x = str(round( float(sanelines[i][0]) ))
		y = str(round( float(sanelines[i][1]) ))
		z = str(round( float(sanelines[i][2]) ))
		
	
		if options.swapyz:

			print("You indicated Y and Z are flipped in the coords file, respect to the tomogram's orientation; therefore, they will be swapped")
			aux = y
			y = z
			z = aux
		
		if options.rotx and options.tomothickness:
			z = float( options.tomothickness ) - float(z)
		
		if options.mult:
			x = str( int(round( float(x) * options.mult )))
			y = str( int(round( float(y) * options.mult )))
			z = str( int(round( float(z) * options.mult )))

		newline = x + '\t' + y + '\t' + z + '\n'
		newlines.append(newline)
	
	output = cfile.replace('.txt','_refactored.txt')
	if options.output:
		output=options.output

	f=open(output,'w')
	f.writelines(newlines)
	f.close()
	#else:
	#	print "ERROR: Terminating. You must specify the output file in .txt format"
	#	sys.exit()
	
	if options.randomize and not options.sort:
		pass
		
	if options.sort and not options.randomze:
		pass
	
	E2end(logger)

	return


def fixlines(inlines):
	n=len(inlines)
	newlines=[]
	for i in range(0,n):
		inlines[i] = inlines[i].replace(", ",' ')	
		inlines[i] = inlines[i].replace(",",' ')
		inlines[i] = inlines[i].replace("x",'')
		inlines[i] = inlines[i].replace("y",'')
		inlines[i] = inlines[i].replace("z",'')
		inlines[i] = inlines[i].replace("=",'')
		inlines[i] = inlines[i].replace("_",' ')
		inlines[i] = inlines[i].replace("\n",'')
		inlines[i] = inlines[i].replace("\t",' ')
		inlines[i] = inlines[i].replace("  ",' ')

		if inlines[i]:
			newlines.append(inlines[i])

	return newlines


def loadlines(infile):
	f = open(infile,'r')
	lines =f.readlines()
	f.close()
	
	outlines = []
	for line in lines:
		#print "The len of this line is", len(line)
		if len(line)<5 or len(line) > 30:

			print("This line is insane and therefore will be removed", line)
		else:
			outlines.append(line)

	return outlines


if '__main__' == __name__:
	main()
