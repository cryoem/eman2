#!/usr/bin/env python

from EMAN2 import *
from mpi_eman import *
import os

proc,nproc=mpi_init()

print "%d\t"%(proc),
os.system('hostname')

mpi_finalize()

if proc==0:
	print "done"


