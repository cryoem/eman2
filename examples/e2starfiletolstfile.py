#!/usr/bin/env python
# Author: Zhao Wang
# Edited by Stephen Murray (scmurray@bcm.edu), 12/12/11
# Copyright (c) 2000-2011 Baylor College of Medicine

# Official copyright notice. EMAN2 is distributed under a joint GPL/BSD license. Please copy
# this statement from one of the other programs. You must agree to use this license if your
# code is distributed with EMAN2. While you may use your own institution for the copyright notice
# the terms of the GPL/BSD license permit us to redistribute it.

# This script will convert a relion star file to an eman2 lst file, looking specifically for the particle number in the image and the image name

from EMAN2 import *
import sys
import os
from EMAN2star import *
import pyemtbx.options

progname = os.path.basename(sys.argv[0])
usage = """ prog [options] <input star file> <output lst file>

 This script will convert a relion star file to an eman2 lst file, looking specifically for the particle number in the image and the image name
Examples: e2starfiletolstfile.py <input star file> <output lst file>


"""

print "Running e2refinetorelion3d.py"
# Required Program Options and Parameters (GUI and Command Line)
parser = EMArgumentParser(usage, version=EMANVERSION)
parser.add_pos_argument(name="input_star", type=str, help="The name of the input star file", default="")
parser.add_pos_argument(name="output_lst", type=str, help="The name of the output .lst file", default="")
parser.add_argument("--path", type=str, default=None, help="Path for the files. Use this if you want to change the output path for the files in the .lst file. For example, instead of the file being located at /raid1/myaccount/data/image.hdf you could do --path=data/particles and the new path written into the .lst file will be data/particles/image.hdf")
optionList = pyemtbx.options.get_optionlist(sys.argv[1:])
(options, args) = parser.parse_args()

if len(args) != 2:
	exit(-1)

#does the input star file exist?
if os.path.exists(args[0]):
	starf = StarFile(args[0])
else:
	exit(-1)

lstfile = LSXFile(args[1])

# in the below line, -1 tells it to append the line, the next value is the image number in that image, and the final value is the image name

for i in range(len(starf[starf.keys()[0]])):
	if options.path:
		filepath1=starf['rlnImageName'][i].split("@")[1].replace(".mrcs",".hdf").split("/")
		filepath = options.path + "/" + filepath1[len(filepath1)-1]
	else:
		filepath = starf['rlnImageName'][i].split("@")[1].replace(".mrcs",".hdf")
	lstfile.write(-1,str(int(starf['rlnImageName'][i].split("@")[0])-1) ,filepath )


# normalize the length of the lines in the LST file
lstfile.normalize() 

print "e2starfiletolstfile.py complete"