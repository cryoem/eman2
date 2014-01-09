from EMAN2 import *
from math import *
import random

a=PointArray()
a.set_number_points(336)
a.set_pot_parms(3.3,100.0,10.0,35.9*pi/180.0,100.0,0,None)

# random distribution
for i in xrange(336):
	#a.set_vector_at(i,Vec3f(random.uniform(-250,250),random.uniform(-250,250),random.uniform(-250,250)),1.0)
	ang=2.0*pi*i/336.0
	a.set_vector_at(i,Vec3f(sin(ang)*200,cos(ang)*150,random.uniform(-2,2)),1.0)

for i in xrange(100000):
	a.minstep(0.1)
	a.center_to_zero()
	print a.potential()
	if i%2000==0:
		img=a.pdb2mrc_by_summation(160,4.52,30.,-1)
		img.write_image("x.hdf",i/2000)

for i in xrange(100000):
	a.minstep(0.01)
	a.center_to_zero()
	if i%2000==0:
		img=a.pdb2mrc_by_summation(160,4.52,30.,-1)
		img.write_image("x.hdf",i/2000+50)

dsts=[]
for r in xrange(336-1):
	x=a.get_vector_at(r)
	y=a.get_vector_at(r+1)
	dsts.append((x-y).length())
	print dsts[-1],

print "\n",sum(dsts)/335
#for i in xrange(50):
	#a.minstep(1.0)
	#a.center_to_zero()
	#img=a.pdb2mrc_by_summation(160,4.52,40.,-1)
	#img.write_image("x.hdf",i+50)
	
