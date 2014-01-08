from EMAN2 import *
from math import *
import random

a=PointArray()
a.set_number_points(336)
a.set_pot_parms(3.3,10.0,10.0,35.9*pi/180.0,100.0,0,None)

# random distribution
for i in xrange(336):
	a.set_vector_at(i,Vec3f(random.uniform(-250,250),random.uniform(-250,250),random.uniform(-250,250)),1.0)


for i in xrange(400):
	a.minstep(10.0)
#	a.center_to_zero()
	if i%8==0:
		img=a.pdb2mrc_by_summation(160,4.52,40.,-1)
		img.write_image("x.hdf",i/8)

for i in xrange(50):
	a.minstep(1.0)
	img=a.pdb2mrc_by_summation(160,4.52,40.,-1)
	img.write_image("x.hdf",i+50)
	
