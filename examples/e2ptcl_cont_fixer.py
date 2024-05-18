#!/usr/bin/env python
# Steve Ludtke 2020
from EMAN2 import *
import numpy as np
from EMAN2_utils import *
from sys import argv
import numpy as np

for f in argv[1:]:
	imgs=EMData.read_images(f)
	r=imgs[0]["nx"]//3
	means=np.array([i.process("mask.sharp",{"outer_radius":r})["mean"] for i in imgs])
	print(f"{f}\t{np.mean(means)} +- {np.std(means)}")
	
