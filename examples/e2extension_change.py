#!/usr/bin/env python
# Author: Jesus Galaz-Montoya 2018 (jgalaz@gmail.com); last update Feb/9
# Copyright (c) 2000-2011 Baylor College of Medicine
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
from __future__ import print_function
from __future__ import division

from EMAN2_utils import *
from EMAN2 import *
import os
import sys
from sys import argv

from EMAN2jsondb import JSTask,jsonclasses

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options] 
	This program changes the extension of all files in a directory by default (or a selected subset) to the desired image format/extension. This could be accomplished
	with e2proc3d.py, except that the user would have to run a command for each file, or know how to make a bash script to do it, and e2proc3d.py isn't parallelized.
	This wrapper performs the task without manual intervention by the user or bash scripting knowledge.
	"""
	
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_argument('--apix', type=float, default=0.0, help="""default=0.0 (not used). Reset the apix on the header of each image to this value.""")
	
	parser.add_argument('--deleteoriginal', action="store_true", default=False, help="""default=False. Deletes the original input images to avoid having redundant copies and save disk space.""")

	parser.add_argument('--inputstring', type=str, default='', help="""default=empty string (not used to filter/select input files). String common to all images to be processed. E.g., if --inputstring=.mrc, all files in the directory ending in '.mrc' (or containing this as part of the filename) will be analyzed.""")

	parser.add_argument('--outputextension', type=str, default=None, help="""default=None. Required, and must be different from the extension of the input images. Possible options are: '.hdf','.mrc','.mrcs','.st','.ali','.rec','.tiff'""")
	#parser.add_argument("--parallel","-P",type=str,help="Run in parallel, specify type:<option>=<value>:<option>:<value>",default=None, guitype='strbox', row=8, col=0, rowspan=1, colspan=2, mode="align")

	parser.add_argument("--ppid", type=int, default=-1, help="Set the PID of the parent process, used for cross platform PPID")

	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n",type=int, default=0, help="verbose level [0-9], higher number means higher level of verboseness.")	

	(options, args) = parser.parse_args()
	
	if not options.outputextension:
		print("\nERROR: parameter --outputextension required.")
		sys.exit(1)
	
	logid=E2init(sys.argv,options.ppid)
	
	extensions=['.hdf','.mrc','.mrcs','.st','.ali','.rec','.tiff','.dm3']
	#dotlessextensions=['hdf','mrc','mrcs','st','ali','rec','tiff','dm3']

	directory=os.getcwd()
	filesindir=os.listdir(directory)
	
	nf=len(filesindir)
	
	filesindir.sort()

	if options.verbose:
		print("\nexamining nf={} files in directory={}".format(nf,directory))
		if options.verbose >9:
			print("\n which are")
			for fyle in filesindir:
				print("\n{}".format(fyle))
	
	files=set([])

	for fyle in filesindir:
		#dotlessextension = fyle[-5:].split('.')[-1]
		#extension = '.' + fyle[-5:].split('.')[-1]
		extension = os.path.splitext(fyle)[-1]
		if extension in extensions and options.inputstring in fyle: #or dotlessextension in dotlessextensions:
			if extension != options.outputextension:
				try:
					EMData(fyle,0,True)
					files.add(fyle)
				except:
					print("\nERROR: Skipping file. Could not read file={}".format(fyle))
			else:
				print("\nWARNING: Skipping file. Extension={} in input file={} is the same as outputextension={}".format(extension,fyle,options.outputextension))
		else:
			print("\nWARNING: Skipping file. invalid extension={} for input file={}".format(extension,fyle))

	files=list(files)
	files.sort()
	n=len(files)
	if n > 0:
		print("\nfound n={} valid files".format(n))
		print("\nwhich are={}".format(files))
	else:
		print("\nERROR: no files found with any of the following valid input formats/extensions={}".format(extensions))
		sys.exit(1)

	i=0
	for fyle in files:
		print("\nconverting file {}/{}, which is file={}".format(i+1,n,fyle))
		extension = '.' + fyle[-5:].split('.')[-1]
		stem = fyle.replace(extension,'')
		outfile = stem + options.outputextension

		cmd = 'e2proc3d.py ' + fyle + ' ' + outfile

		if options.apix:
			cmd += ' --apix ' + str(options.apix)
		#	print("\napix will be changed to {}".format(options.apix))
		#	cmd += ' && e2procheader.py --input ' + outfile + ' --stem apix --stemval ' + str(options.apix) + ' --valtype float'

		print('\ncmd to run is {}'.format(cmd))
		feedback = runcmd(options,cmd)
		print("\nfeedback is".format(feedback))
		print("\ndone")
		
		if options.deleteoriginal:
			
			c = os.getcwd()
			findir = os.listdir(c)
			if outfile in findir:
				delete=False
				try:
					EMData(outfile,0,True)
					delete=True
				except:
					print("\nWARNING: could not read outfile={}".format(outfile))

				if delete:
					print("\ndeleting original image, as requested")
					os.remove(fyle)
				else:
					print("\nWARNING: original file={} not removed because outfile={} seems to be faulty".format(fyle,outfile))
			else:
				print("\nWARNING: intended outfile={} not present in directory".format(outfile))
		i+=1

	E2end(logid)

	return


if __name__ == '__main__':
	main()

