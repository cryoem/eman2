from EMAN2 import *


eer = 'FoilHole_29355520_Data_29330528_29330530_20200329_234551_Fractions.mrc.eer'
# print(test_image())
im = EMData()
# ims = im.read_images('FoilHole_29355520_Data_29330528_29330530_20200329_234551_Fractions.mrc.eer')
im.sum_images(eer)
im.write_image('zsum.hdf')
# print(len(ims))

frames = EMData.read_images(eer)


# :TODO
# mrcs write_images
# hdf write image header only
# eer
# e2proc write compressed
