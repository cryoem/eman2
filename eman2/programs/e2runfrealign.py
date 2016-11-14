#!/usr/bin/env python
# Author: Stephen Murray (scmurray@bcm.edu), 2/14/11
# Copyright (c) 2000-2011 Baylor Colelge of Medicine

# Official copyright notice. EMAN2 is distributed under a joint GPL/BSD license. Please copy
# this statement from one of the other programs. You must agree to use this license if your
# code is distributed with EMAN2. While you may use your own institution for the copyright notice
# the terms of the GPL/BSD license permit us to redistribute it.

# import block
import os
import sys
from subprocess import *
import shutil
from EMAN2 import *

usage = """prog
...
"""

parser = EMArgumentParser(usage,version=EMANVERSION)

parser.add_header(name="runfrealign", help='Click Launch to run Frealign', title="### Click Launch to run e2runfrealign ###", row=0, col=0, rowspan=1, colspan=1)
parser.add_pos_argument(name="dir",help="The refinement directory to use for FreAlign.", default="", guitype='dirbox', dirbasename='frealign',  row=1, col=0,rowspan=1, colspan=2)
parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)

(options, args) = parser.parse_args()

if len(args) != 1:
   print "Please use e2runfrealign.py with only one argument - which FreAlign directory to run in"
   exit(-1)
os.chdir(str(args[0]))
E2n=E2init(args,options.ppid)

shutil.copy('3DMapInOut.mrc', '3DMapInOut.mrc.old')
shutil.copy('ptcl_meta_data', 'ptcl_meta_data.old')

try:
	retcode = call('frealign_v9.exe < card.txt', shell=True)
	if retcode < 0:
		print >>sys.stderr, "Child was terminated by signal", -retcode
	else:
		print >>sys.stderr, "child returned", retcode
except OSError, e:
	print >>sys.stderr, "Execution Failed:", e

shutil.copy('3DMapInOut.mrc', 'OutputMap.mrc')
shutil.move('3DMapInOut.mrc.old', '3DMapInOut.mrc')
shutil.move('ptcl_meta_data.old', 'ptcl_meta_data')
print "e2runfrealign.py finished"

E2end(E2n)
