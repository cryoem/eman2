#!/usr/bin/env python

from EMAN2 import *
from mpi_eman import *
import socket
import os

proc,nproc=mpi_init()

print "%d\t%s"%(proc,socket.gethostname())

mpi_finalize()

if proc==0:
	print "done"


