#!/usr/bin/env python
# This program performs simple processing of .LST files

# Author: Steven Ludtke, 11/16/25 (sludtke@bcm.edu)
# Copyright (c) 2025- Baylor College of Medicine
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

from EMAN3 import *
from math import *
import numpy as np
import os
import sys

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """Usage:\ne3proclst.py [options] <lst 1> <lst 2> ... \nManipulate .lst files in various ways\nMany of the options in this program act in-place rather than creating a new file, please read the help carefully!\nIf your goal is to produce an actual image file rather than the sort of virtual stack represented by .lst files,\nuse e2proc2d.py or e2proc3d.py instead. Those other programs will treat LST files as normal image files for input.\n."""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	####################
#	parser.add_argument("--average", action="store_true", help="Averages all input images (without alignment) and writes a single output image")
	
	parser.add_argument("--output", type=str, default=None, help="The input file(s) may be other .lst files or image files. All inputs will be merged sequentially into this single output. If this .lst file already exists, new images will be appended to it.")
	parser.add_argument("--eman2to3", action="store_true", help="This will convert an EMAN2 .lst file referencing phase flipped particle images into an EMAN3 style .lst file containing the original particles (if available) with embedded CTF information")

# 	parser.add_argument("--eosplit", action="store_true", help="Will generate _even and _odd .lst files for each specified input .lst file")
# 	parser.add_argument("--split", type=int,default=0, help="Will put every nth particle in a separate numbered .lst file based on --create name. Ignores other subset selection options! Single input only!")
#
# 	parser.add_argument("--dereforig", type=str, default=None, help="Extract the data_source and data_n parameters from each image in the file and create a new .lst file referencing the original image(s)")
#
# 	parser.add_argument("--exclude", type=str, default=None, help="only works if --create is supplied. comma-separated list of indexes from the input file(s) to EXCLUDE from the created .lst file.")
#
# 	#parser.add_argument("--first", type=int, default=0, help="Default=0 (first image index in input(s)). This will be the first particle index in the input images to put in the output lsx/lst file.")
#
# 	parser.add_argument("--include", type=str, default=None, help="only works if --create is supplied. comma-separated list of indexes to take from the input file(s) to INCLUDE in the created .lst file. if you have the list of indexes to include in a .txt file, you can provide it through --list.")
# 	parser.add_argument("--inplace", action="store_true", default=False, help="only works with --create. if the stack specified in --create already exists, this will prevent appending to it. rather, the file will be modified in place.")
# #	parser.add_argument("--force", action="store_true", default=False, help="only works with --create. if the stack specified in --create already exists, it will be removed and rewritten.")
#
# 	#parser.add_argument("--last", type=str, default=-1, help="Default=-1 (last image index in input (s)). This will be the first particle index in the input images to put in the output lsx/lst file.")
# 	parser.add_argument("--list", type=str, default=None, help="only works if --create is supplied. .txt file with a list of indexes (one per line/row) to take from the input file(s) to INCLUDE in the created .lst file.")
#
# 	parser.add_argument("--merge", type=str, default=None, help="Specify the output name here. This will concatenate all of the input .lst files into a single output")
# 	parser.add_argument("--mergesort", type=str, default=None, help="Specify the output name here. This will merge all of the input .lst files into a single (resorted) output")
# 	parser.add_argument("--mergeinterleave", type=str, default=None, help="Specify the output name here. Interleaves images from input .lst files, eg - A0,B0,C0,A1,B1,C1,... truncates based on size of smallest input, eg- 1000,500,300 -> 900")
# 	parser.add_argument("--mergeeo", action="store_true", default=False, help="Merge even odd lst.")
# 	parser.add_argument("--mergeref", type=str, default=None, help="reference lst file to determine the order for --create")
# 	parser.add_argument("--replacexf", type=str, default=None, help="replace xform.projection or xform.align3d in the particle set with the same attribute from another set. Require --create")
# 	parser.add_argument("--replaceclass", action="store_true", default=False, help="following --replacexf only. also take the class entry from the given list.")
#
# 	parser.add_argument("--minhisnr", type=float, help="Integrated SNR from 1/10-1/4 1/A must be larger than this",default=-1,guitype='floatbox', row=8, col=1)
# 	parser.add_argument("--minlosnr", type=float, help="Integrated SNR from 1/200-1/20 1/A must be larger than this",default=-1,guitype='floatbox', row=8, col=0)
# 	parser.add_argument("--mindf", type=float, help="Minimum defocus",default=-1,guitype='floatbox', row=8, col=1)
# 	parser.add_argument("--maxdf", type=float, help="Maximum defocus",default=-1,guitype='floatbox', row=8, col=0)
# 	parser.add_argument("--maxscore", type=float, help="Per-particle score used to exclude particles above a specified value. Works only with --create, with .lst inputs. Some other options may also prevent it. ",default=1)
#
# 	parser.add_argument("--numaslist", type=str, default=None, help="extract the particle indexes (numbers) only from an lst file into a text file (one number per line).")
#
#
# 	parser.add_argument("--keeptilt", type=int, default=-1, help="Keep X subtilt close to the center tilt. Require ptcl3d_id and tilt_id.")
#
# 	parser.add_argument("--range", type=str, default=None, help="Range of particles to use. Works only with --create option. Input of 0,10,2 means range(0,10, step=2).")
# 	parser.add_argument("--retype", type=str, default=None, help="If a lst file is referencing a set of particles from particles/imgname__oldtype.hdf, this will change oldtype to the specified string in-place (modifies input files), use 'none' to remove the type.")
# 	parser.add_argument("--refile", type=str, default=None, help="similar to retype, but replaces the full filename of the source image file with the provided string")
# 	parser.add_argument("--shuffle", action="store_true", default=False, help="shuffle list inplace.")
# 	parser.add_argument("--flip", action="store_true", default=False, help="flip xform.")
# 	parser.add_argument("--unique", action="store_true", default=False, help="keep only unique particles.")
# 	parser.add_argument("--sym", type=str, default=None, help="WARNING: operates in-place, modifying input files!  Apply symmetry to .lst files of particles with xform.projection by duplicating each particle N times.")
# 	parser.add_argument("--extractattr", type=str, default=None, help="Extract an attribute from the particle image file and add it as an attribute in the .lst file")
# 	parser.add_argument("--getclass", type=int, help="select a class when --create",default=-1)
#
# 	parser.add_argument("--getptcls", action="store_true", default=False, help="Get particles from input 2D class averages and put them in a lst specified in --create")
# 	parser.add_argument("--nocomments", action="store_true", default=False, help="Removes the comments from each line of the lst file.")
# 	parser.add_argument("--scalexf", type=float, help="scale the translation in xform in header.",default=-1)
# 	parser.add_argument("--applyxf", type=str, help="apply the transform to a list of particles with xform.projection or align3d.",default=None)

	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, help="verbose level [0-9], higher number means higher level of verboseness",default=1)
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)


	(options, args) = parser.parse_args()

	if len(args)<1 :
		parser.error("At least one lst file required")
		sys.exit(1)

	logid=E3init(sys.argv,options.ppid)

	lut={}
	if options.output is not None:
		out=LSXFile(options.output)

		for infile in args:
			# if the input is already a .lst file
			if infile.endswith(".lst"):
				try: lsxin=LSXFile(infile,ifexists=True)
				except: error_exit(f"{infile} does not appear to be an existing .lst file in LST or LSX format")
				N=len(lsxin)
			else:
				lsxin=None
				try: N=file_image_count(infile)
				except: error_exit(f"{infile} does not appear to be a valid image file")

			# N should now be the number of images in this input files
			for i in range(N):
				if lsxin is not None:
					nextf,extf,dct=lsxin[i]

				else:
					nextf=i
					extf=infile
					dct={}

				if options.eman2to3:
					# needs to happen before we update extf
					if "ctf" not in dct:
						try: ctf=EMData(extf,nextf,True)["ctf"]
						except: error_exit(f"No CTF information available in original image {extf} or .lst file {infile}")
						ctf.background=[]
						ctf.snr=[]
						dct["ctf"]=ctf

					if "__" in extf:
						if extf not in lut:
							if os.path.exists(extf.split("__")[0]+".hdf"): lut[extf]=extf.split("__")[0]+".hdf"
							elif os.path.exists(extf.split("__")[0]+"_ptcls.hdf"): lut[extf]=extf.split("__")[0]+"_ptcls.hdf"
							else: error_exit(f"Cannot find non-ctf flipped file for {extf}")
						extf=lut[extf]

				out[-1]=(nextf,extf,dct)	# Append to the end of the output LST file


	E3end(logid)


if __name__ == "__main__":
	main()
