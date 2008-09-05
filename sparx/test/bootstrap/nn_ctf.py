#!/usr/bin/env python

from EMAN2 import *
from sparx import *

from sys import argv
from sys import exit
from string import atoi

from time import time

if len(argv) != 3 and len(argv) != 6 :
    print "usage: nn_ctf proj_stack vol_stack [start end step]"
    exit(-1)


proj_stack = argv[1]
vol_stack  = argv[2]

if len(argv) == 3 :
    start = 0
    end = EMUtil.get_image_count( proj_stack )
    step = 1
else :
    start = atoi( argv[3] )
    end   = atoi( argv[4] )
    step  = atoi( argv[5] )

time_start = time()

v = recons3d_4nn_ctf( proj_stack, range(start,end,step), snr=100.0 )

dropImage( v, vol_stack )

print "Time:", time() - time_start


