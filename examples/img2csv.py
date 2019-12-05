#!/usr/bin/env python
from __future__ import division

### This program will read the first image from an image file
### and write the image data in a csv file

from EMAN2 import *
import numpy as np
from sys import argv

im=EMData(argv[1],0)
data=to_numpy(im).copy()
np.savetxt(argv[1].rsplit(".",1)[0]+".csv",data,delimiter=",")

