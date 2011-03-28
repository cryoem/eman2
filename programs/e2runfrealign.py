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

if len(sys.argv) != 1:
   print "Please use e2runfrealign.py"
   exit(-1)

E2n=E2init(sys.argv)

os.system('clear')
shutil.copy('3DMapInOut.mrc', '3DMapInOut.mrc.old')
shutil.copy('ptcl_meta_data', ' ptcl_meta_data.old')

try:
	retcode = call('frealign_v8.exe < card.txt', shell=True)
	if retcode < 0:
		print >>sys.stderr, "Child was terminated by signal", -retcode
	else:
		print >>sys.stderr, "child returned", retcode
except OSError, e:
	print >>sys.stderr, "Execution Failed:", e

shutil.copy('3DMapInOut.mrc', ' OutputMap.mrc')
shutil.move('3DMapInOut.mrc.old', '3DMapInOut.mrc')
shutil.move('ptcl_meta_data.old', ' ptcl_meta_data')
print "e2runfrealign.py finished"

E2end(E2n)
