from EMAN2 import *

# print(test_image())
im = EMData()
# ims = im.read_images('FoilHole_29355520_Data_29330528_29330530_20200329_234551_Fractions.mrc.eer')
im.sum_images('FoilHole_29355520_Data_29330528_29330530_20200329_234551_Fractions.mrc.eer')
im.write_image('zsum5.hdf')
# print(len(ims))
