#!/bin/env python

from sys import argv
import os
from EMAN2 import *

img=EMData()
for fsp in argv[1:]:
	img.read_image(fsp,0,1)
	print "%f\t%f"%(-img.get_ctf().defocus,img.get_ctf().bfactor)

