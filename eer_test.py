import os

from EMAN2 import *

print(os.getpid())

import pdb;pdb.set_trace()
fname='FoilHole_24015405_Data_24016401_24016403_20200225_0014_Fractions.mrc.eer'

a = EMData()
imgs = a.read_images(fname)

print(len(imgs))