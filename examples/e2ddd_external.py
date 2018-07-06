#!/usr/bin/env python
from __future__ import print_function

#
# Author: Steven Ludtke, 02/12/2013 (sludtke@bcm.edu). Updated on 08/28/16.
# Modified by James Michael Bell, 03/27/2017 (jmbell@bcm.edu)
# Copyright (c) 2000-2013 Baylor College of Medicine
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
from numpy import *
import pprint
import sys
import os
from sys import argv
from time import sleep,time,ctime
import threading
import Queue
import numpy as np
from sklearn import linear_model
from scipy import optimize

def which(program):
	# found at https://stackoverflow.com/questions/377017/test-if-executable-exists-in-python
	def is_exe(fpath) : return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

	fpath, fname = os.path.split(program)
	if fpath:
		if is_exe(program): return program
	else:
		for path in os.environ["PATH"].split(os.pathsep):
			for f in os.listdir(path):
				if program in f:
					exe_file = os.path.join(path, f)
					if is_exe(exe_file): return exe_file
	return None

def main():

	progname = os.path.basename(sys.argv[0])
	usage = """prog [options] <ddd_movie_stacks>

	This program is a wrapper for various external DDD alignment routines including:
		MotionCor2
		alignframes (IMOD)

	Note, this is a simple script that uses the default alignment parameters for these
	external programs. Programs must be installed and accessible via the PATH 
	environment variable. To customize alignment, you will need to run these programs 
	from the command line independently and import the aligned averages into your project.
	"""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_pos_argument(name="movies",help="List the movies to align.", default="", guitype='filebox', browser="EMMovieDataTable(withmodal=True,multiselect=True)",  row=0, col=0,rowspan=1, colspan=2, mode="movie")

	parser.add_header(name="orblock1", help='Just a visual separation', title="Options", row=1, col=0, rowspan=1, colspan=1, mode="tomo")
	
	parser.add_argument("--device",  default = "GPU", type=str, choices=["CPU","GPU"], help="When possible, use this device to process movie frames. Default is GPU.",guitype='combobox', choicelist='["GPU","CPU"]', row=1, col=0, rowspan=1, colspan=1, mode="movie")
	parser.add_argument("--program",  default = "IMOD_alignframes", type=str, choices=["IMOD_alignframes","MotionCor2"], help="Use this external program to align frames. Choose between IMOD_alignframes and MotionCor2. Note, programs must be accessible from your PATH environment variable.",guitype='combobox', choicelist='["IMOD_alignframes","MotionCor2"]', row=1, col=1, rowspan=1, colspan=1, mode="movie")
	
	parser.add_argument("--dark",  default = None, type=str, help="Use this dark reference.",guitype='filebox',  browser="EMMovieDataTable(withmodal=True,multiselect=True)",  row=2, col=0,rowspan=1, colspan=2, mode="align",mode="movie")
	parser.add_argument("--gain",  default = None, type=str, help="Use this gain reference.",guitype='filebox',  browser="EMMovieDataTable(withmodal=True,multiselect=True)",  row=3, col=0,rowspan=1, colspan=2, mode="align",mode="movie")

	# rotation? IMOD
	parser.add_argument("--flipgain",  default = "1 1", type=str, help="Use this many patches with MotionCor2. Format is 'X Y'. Default is: '1 1'",guitype='strbox', row=4, col=0, rowspan=1, colspan=1, mode="movie")
	parser.add_argument("--rotgain",  default = "1 1", type=str, help="Use this many patches with MotionCor2. Format is 'X Y'. Default is: '1 1'",guitype='strbox', row=4, col=0, rowspan=1, colspan=1, mode="movie")

	parser.add_argument("--group",  default = "1 1", type=str, help="Use this many patches with MotionCor2. Format is 'X Y'. Default is: '1 1'",guitype='strbox', row=4, col=0, rowspan=1, colspan=1, mode="movie")
	parser.add_argument("--iters",  default = "1 1", type=str, help="Use this many patches with MotionCor2. Format is 'X Y'. Default is: '1 1'",guitype='strbox', row=4, col=0, rowspan=1, colspan=1, mode="movie")
	parser.add_argument("--binby",  default = "1 1", type=str, help="Use this many patches with MotionCor2. Format is 'X Y'. Default is: '1 1'",guitype='strbox', row=4, col=0, rowspan=1, colspan=1, mode="movie")

	# throw
	parser.add_argument("--first",  default = "1 1", type=str, help="Use this many patches with MotionCor2. Format is 'X Y'. Default is: '1 1'",guitype='strbox', row=4, col=0, rowspan=1, colspan=1, mode="movie")
	# trunc
	parser.add_argument("--last",  default = "1 1", type=str, help="Use this many patches with MotionCor2. Format is 'X Y'. Default is: '1 1'",guitype='strbox', row=4, col=0, rowspan=1, colspan=1, mode="movie")

	parser.add_header(name="orblock1", help='Just a visual separation', title="Options", row=1, col=0, rowspan=1, colspan=1, mode="tomo")

	parser.add_argument("--mdoc",  default = "1 1", type=str, help="Use this many patches with MotionCor2. Format is 'X Y'. Default is: '1 1'",guitype='strbox', row=4, col=0, rowspan=1, colspan=1, mode="movie")

	parser.add_argument("--smooth",  default = "1 1", type=str, help="Use this many patches with MotionCor2. Format is 'X Y'. Default is: '1 1'",guitype='strbox', row=4, col=0, rowspan=1, colspan=1, mode="movie")
	parser.add_argument("--maxshift",  default = "1 1", type=str, help="Use this many patches with MotionCor2. Format is 'X Y'. Default is: '1 1'",guitype='strbox', row=4, col=0, rowspan=1, colspan=1, mode="movie")
	parser.add_argument("--kfactor",  default = "1 1", type=str, help="Use this many patches with MotionCor2. Format is 'X Y'. Default is: '1 1'",guitype='strbox', row=4, col=0, rowspan=1, colspan=1, mode="movie")
	parser.add_argument("--pairwiseframes",  default = "1 1", type=str, help="Use this many patches with MotionCor2. Format is 'X Y'. Default is: '1 1'",guitype='strbox', row=4, col=0, rowspan=1, colspan=1, mode="movie")

	parser.add_argument("--hybrid",  default = "1 1", type=str, help="Use this many patches with MotionCor2. Format is 'X Y'. Default is: '1 1'",guitype='strbox', row=4, col=0, rowspan=1, colspan=1, mode="movie")
	parser.add_argument("--antialias",  default = "1 1", type=str, help="Use this many patches with MotionCor2. Format is 'X Y'. Default is: '1 1'",guitype='strbox', row=4, col=0, rowspan=1, colspan=1, mode="movie")

	parser.add_header(name="orblock1", help='Just a visual separation', title="Options", row=1, col=0, rowspan=1, colspan=1, mode="tomo")

	parser.add_argument("--patch",  default = "1 1", type=str, help="Use this many patches with MotionCor2. Format is 'X Y'. Default is: '1 1'",guitype='strbox', row=4, col=0, rowspan=1, colspan=1, mode="movie")
	parser.add_argument("--mag",  action="store_true", help="Use this many patches with MotionCor2. Format is 'X Y'. Default is: '1 1'",guitype='boolbox', row=4, col=0, rowspan=1, colspan=1, mode="movie")
	parser.add_argument("--voltage",  default = "1 1", type=str, help="Use this many patches with MotionCor2. Format is 'X Y'. Default is: '1 1'",guitype='strbox', row=4, col=0, rowspan=1, colspan=1, mode="movie")
	parser.add_argument("--apix",  default = "1 1", type=str, help="Use this many patches with MotionCor2. Format is 'X Y'. Default is: '1 1'",guitype='strbox', row=4, col=0, rowspan=1, colspan=1, mode="movie")
	parser.add_argument("--dose",  default = "1 1", type=str, help="Use this many patches with MotionCor2. Format is 'X Y'. Default is: '1 1'",guitype='strbox', row=4, col=0, rowspan=1, colspan=1, mode="movie")
	parser.add_argument("--initdose",  default = "1 1", type=str, help="Use this many patches with MotionCor2. Format is 'X Y'. Default is: '1 1'",guitype='strbox', row=4, col=0, rowspan=1, colspan=1, mode="movie")
	parser.add_argument("--bfactor",  default = "1 1", type=str, help="Use this many patches with MotionCor2. Format is 'X Y'. Default is: '1 1'",guitype='strbox', row=4, col=0, rowspan=1, colspan=1, mode="movie")
	parser.add_argument("--tolerance",  default = "1 1", type=str, help="Use this many patches with MotionCor2. Format is 'X Y'. Default is: '1 1'",guitype='strbox', row=4, col=0, rowspan=1, colspan=1, mode="movie")

	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	(options, args) = parser.parse_args()

	device = options.device

	if options.program == "IMOD_alignframes":
		program = which("alignframes")
	
	if options.program == "MotionCor2":
		try: patchx,patchy = options.patch.split(" ")
		except:
			print("Could not interpret --patch for use with MotionCor2. Please input integer values in the form: 'X,Y'")
			sys.exit(1)
		if device == "CPU":
			print("INPUT ERROR: Cannot use --device CPU with MotionCor2. This program requires 1+ GPU.")
			sys.exit(1)
		program = which("MotionCor2")

	if program == None:
		print("Could not locate {}. Please check your PATH environment variable.".format(options.program.replace("_"," ")))
		sys.exit(1)

	if args[0].split(".")[-1] not in [".mrc",".mrcs"]:
		print("Input files must be in .mrc format. Please reformat these files and try again.")
		sys.exit(1)

	for input_file in sorted(args):

		output_file = input_file.split(".")[0]+"_ali.mrc"

		if options.program == "IMOD_alignframes":
			cmd = "alignframes -input {} -output {}".format(input_file,output_file)
			if device == "GPU": cmd += " -gpu 0 "
			if options.dark != None: cmd += " -dark {} ".format(options.dark)
			if options.gain != None: cmd += " -gain {} ".format(options.gain)

		elif options.program == "MotionCor2":
			cmd = "{} -InMrc {} -OutMrc {} ".format(program,input_file,output_file)
			if options.patch != 1 1: cmd += " -Patch {} ".format(options.patch)
			if options.dark != None: cmd += " -Dark {} ".format(options.dark)
			if options.gain != None: cmd += " -Gain {} ".format(options.gain)
		
		run(cmd,shell=True)

def run(cmd,shell=False,cwd=None,exe="/bin/sh"):
	if cwd == None: cwd = os.getcwd()
	if shell == False: cmd = cmd.split()
	p = psutil.Popen("/usr/bin/time -p "+cmd, shell=shell, cwd=cwd, executable=exe,stderr=subprocess.PIPE) 
	_,err=p.communicate()
	return