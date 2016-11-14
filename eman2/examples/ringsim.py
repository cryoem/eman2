from EMAN2 import *
from math import *
import random

a=PointArray()
a.set_number_points(336)
a.sim_set_pot_parms(3.3,10.0,1000.0,35.9*pi/180.0,20.0,0,None)

# random distribution
for i in xrange(336):
	#a.set_vector_at(i,Vec3f(random.uniform(-250,250),random.uniform(-250,250),random.uniform(-250,250)),1.0)
	ang=2.0*pi*i/336.0
	a.set_vector_at(i,Vec3f(sin(ang)*300,cos(ang)*200,random.uniform(5,5)),1.0)

for i in xrange(40000):
	a.sim_minstep_seq(0.2)
	a.center_to_zero()
	if i%1000==0:
#		a.sim_rescale()
		print a.sim_potential()
		a.sim_printstat()
		img=a.pdb2mrc_by_summation(160,4.52,10.,-1)
		img.write_image("x.hdf",i/1000)

#a.sim_rescale()
a.sim_set_pot_parms(3.3,500.0,1000.0,35.9*pi/180.0,500.0,0,None)
print "------------"
for i in xrange(40000):
	a.sim_minstep_seq(0.02)
	a.center_to_zero()
	if i%1000==0:
#		a.sim_rescale()
		print a.sim_potential()
		a.sim_printstat()
		img=a.pdb2mrc_by_summation(160,4.52,10.,-1)
		img.write_image("x.hdf",i/1000+40)

# Compute all distances
#dsts=[]
#for r in xrange(336-1):
	#x=a.get_vector_at(r)
	#y=a.get_vector_at(r+1)
	#dsts.append((x-y).length())
	#print dsts[-1],

#print "\n",sum(dsts)/335
	
