#!/usr/bin/env python
#********************************************************************************
# Author: Stephen Murray (scmurray@bcm.edu), 7/1/14
# Copyright (c) 2000-2013 Baylor College of Medicine

# Official copyright notice. EMAN2 is distributed under a joint GPL/BSD license. Please copy
# this statement from one of the other programs. You must agree to use this license if your
# code is distributed with EMAN2. While you may use your own institution for the copyright notice
# the terms of the GPL/BSD license permit us to redistribute it.
#********************************************************************************

from EMAN2 import *
from EMAN2db import db_open_dict
import pyemtbx.options
import os
import sys
import shutil
from subprocess import *
from distutils.spawn import find_executable

EMAN_MASK="mask.hdf"
EMAN_ODD="threed_odd_unmasked.hdf"
EMAN_EVEN="threed_even_unmasked.hdf"
refine_dir = sys.argv[1]
progname = os.path.basename(sys.argv[0])
usage = """add stuff


"""

print "Running e2runresmap.py"
# Required Program Options and Parameters (GUI and Command Line)
parser = EMArgumentParser(usage, version=EMANVERSION)
parser.add_header(name="resmap", help="Options in this GUI pertain to using ResMap", title="---ResMap Options---", row=0, col=0, rowspan=1, colspan=3)
#parser.add_pos_argument(name="split_map",help="The name of one of the unfiltered map halves", default="", guitype='filebox', browser="EMBrowserWidget(withmodal=True,multiselect=False)", row=1, col=0, rowspan=2, colspan=2)
parser.add_pos_argument(name="dir",help="The refinement directory to ResMap on.", default="", guitype='dirbox', dirbasename='refine',  row=1, col=0,rowspan=1, colspan=1)
parser.add_argument("--use_mask", action="store_true",help="Use the final EMAN2 mask?",default=False, guitype='boolbox',row=1, col=1,rowspan=1, colspan=1)
parser.add_argument("--res_step", help="Step size (in Angstroms)", default=1.0,  guitype='floatbox', row=2, col=0, rowspan=1, colspan=1)
parser.add_argument("--res_min", help="Minimum resolution (in Angstroms)", default=0.0,  guitype='floatbox', row=2, col=1, rowspan=1, colspan=1)
parser.add_argument("--res_max", help="Maximum resolution (in Angstroms)", default=0.0,  guitype='floatbox', row=2, col=2, rowspan=1, colspan=1)
parser.add_argument("--p_value", help="Confidence level (Usually between .01 and .05) ", default=0.05,  guitype='floatbox', row=3, col=0, rowspan=1, colspan=1)
#parser.add_argument("--optional_mask",type=str, help="Location of the optional mask to be used", guitype='filebox',default="", browser="EMBrowserWidget(withmodal=True,multiselect=False)", filecheck=False,row=15, col=0, rowspan=2, colspan=2)
parser.add_argument("--verbose", type=int, help="Set the level of verbosity for the code", default=1, guitype='combobox', choicelist='0,1,2,3,4,5,6,7,8,9', row=26, col=1, rowspan=1, colspan=1, expert=True)
parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
optionList = pyemtbx.options.get_optionlist(sys.argv[1:])
(options, args) = parser.parse_args()

# Create the E2ResMap directory structure if it does not exist
i = 1
found = 1
mask=False
while found == 1:
	if i < 10:
		res_run = '0' + str(i)
	else:
		res_run = str(i)
	found = os.path.exists("resmap_" + res_run)
	i = i+1
E2RES = "resmap_" + res_run
os.mkdir(E2RES)

#get the two map halves and the mask (if needed)
if EMAN_EVEN in os.listdir(refine_dir) and EMAN_ODD in os.listdir(refine_dir):
	apix = EMUtil.get_all_attributes(refine_dir+"/"+EMAN_ODD,'apix_x')[0]
	s = "e2proc3d.py " + refine_dir + "/" + EMAN_ODD + " " + E2RES + "/" + EMAN_ODD.replace("hdf","mrc")
	call(s,shell=True)
	s = "e2proc3d.py " + refine_dir + "/" + EMAN_EVEN + " " + E2RES + "/" + EMAN_EVEN.replace("hdf","mrc")
	call(s,shell=True)
else:
	print "Missing one of the two halves! Make sure you are trying to run this program on a completed EMAN2 refinement"
	exit(-1)
if "use_mask" in optionList:
	if not EMAN_MASK in os.listdir(refine_dir):
		print "Mask not present! Make sure you are trying to run this program on a completed EMAN2 refinement" 
		exit(-1)
	else:
		s = "e2proc3d.py " + refine_dir + "/" + EMAN_MASK + " " + E2RES + "/mask.mrc"
		call(s,shell=True)
		mask = True

#Can we find the Resmap Executable?
for path in os.environ["PATH"].split(os.pathsep):
	if os.path.isfile(str(path)+"/ResMap-1.1.4-linux64"):
		resmap=str(path)+"/ResMap-1.1.4-linux64"
		exe = True
if not exe:
	print "Could not find the ResMap executable. Make sure it is in the $PATH environment variable."
	exit(-1)
	
s = resmap+ " --noguiSplit " + E2RES + "/" + EMAN_ODD.replace("hdf","mrc") + " " + E2RES + "/" + EMAN_EVEN.replace("hdf","mrc") + " --vxSize=" + str(apix) + " --pVal=" + str(options.p_value)+ " --stepRes=" + str(options.res_step)

if float(options.res_min) != 0.0:
	s = s + " --minRes=" + str(options.res_min)
if float(options.res_max) != 0.0:
	s = s + " --maxRes=" + str(options.res_max)
if mask:
	s = s + " --maskVol=" + E2RES + "/mask.mrc"
s = s + ">" + E2RES + "/resmap_output_text.txt"
print s
call(s,shell=True)

print "e2runresmap complete"