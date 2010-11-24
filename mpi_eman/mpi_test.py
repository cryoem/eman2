from EMAN2 import *
from mpi_eman import *

proc=mpi_init()

print "running on CPU ",proc

mpi_finalize()

