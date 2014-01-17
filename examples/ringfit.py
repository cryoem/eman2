from EMAN2 import *
import sys
import numpy.linalg as LA
import random

p=EMData(sys.argv[1],0)

print "Preprocessing"
# This will (hopefully) isolate the mini-circle
p.process_inplace("normalize.edgemean")
p.process_inplace("mask.dust3d",{"voxels":8000,"threshold":1.6})
p.process_inplace("mask.auto3d",{"nshells":2,"nshellsgauss":3,"radius":35,"threshold":1.6})
p.process_inplace("threshold.belowtozero",{"minval":0})
p.process_inplace("xform.centerofmass")
comxf=p["xform.align3d"]
p.process_inplace("normalize.unitsum")
p.mult(10000.0)

print "Inertia matrix and alignment"
# compute the resulting inertia matrix
an=Analyzers.get("inertiamatrix",{"verbose":0})
an.insert_image(p)
mxi=an.analyze()
mx=EMNumPy.em2numpy(mxi[0])

# Compute the eigenvalues/vectors
eigvv=LA.eig(mx)		# a 3-vector with eigenvalues and a 3x3 with the vectors
if min(eigvv[0])==0 :
	print "error on ",pf
	drn-=1
	sys.exit(1)

eig=[(1.0/eigvv[0][i],eigvv[1][:,i]) for i in xrange(3)]  # extract for sorting
eig=sorted(eig)		# now eig is sorted in order from major to minor axes
#T=array([eig[0][1],eig[1][1],eig[2][1]])            # reassemble sorted matrix

T=Transform((float(i) for i in (eig[0][1][0],eig[0][1][1],eig[0][1][2],0,eig[1][1][0],eig[1][1][1],eig[1][1][2],0,eig[2][1][0],eig[2][1][1],eig[2][1][2],0)))
#T=Transform((float(i) for i in (eig[0][1][0],eig[1][1][0],eig[2][1][0],0,eig[0][1][1],eig[1][1][1],eig[2][1][1],0,eig[0][1][2],eig[1][1][2],eig[2][1][2],0)))

#T.printme()
#print "------------------"
p.transform(T)
#p.write_image("xf/"+dr+".hdf",drn)

### ok, p now contains the masked/aligned mini-circle we want to fit

# start by making its density uniform and smearing it out a bit
m=p.process("threshold.binary",{"value":p["mean"]+p["sigma"]*5.0})
m.process_inplace("filter.lowpass.gauss",{"cutoff_abs":0.04})

pa=PointArray()
pa.set_number_points(336)

for i in xrange(336):
	ang=2.0*pi*i/336.0
	pa.set_vector_at(i,Vec3f(random.uniform(5,5),sin(ang)*200,cos(ang)*300),1.0)

pa.sim_set_pot_parms(3.3,50.0,500.0,35.9*pi/180.0,20.0,5000.0,m)
#pa.sim_set_pot_parms(3.3,10.0,1000.0,35.9*pi/180.0,20.0,0.0,m)

for i in xrange(60000):
	pa.sim_minstep_seq(0.2)
#	pa.center_to_zero()
	if i==10000: pa.sim_set_pot_parms(3.3,50.0,500.0,35.9*pi/180.0,20.0,1500.0,m)
	elif i==15000: pa.sim_set_pot_parms(3.3,50.0,500.0,35.9*pi/180.0,20.0,2500.0,m)
	elif i==20000: pa.sim_set_pot_parms(3.3,50.0,500.0,35.9*pi/180.0,20.0,1500.0,m)
	elif i==25000: pa.sim_set_pot_parms(3.3,50.0,500.0,35.9*pi/180.0,20.0,2500.0,m)
	elif i==30000: pa.sim_set_pot_parms(3.3,50.0,500.0,35.9*pi/180.0,20.0,1500.0,m)
	elif i==35000: pa.sim_set_pot_parms(3.3,50.0,500.0,35.9*pi/180.0,20.0,2500.0,m)
	elif i==40000: pa.sim_set_pot_parms(3.3,50.0,500.0,35.9*pi/180.0,20.0,1500.0,m)
	elif i==45000: pa.sim_set_pot_parms(3.3,50.0,500.0,35.9*pi/180.0,20.0,2500.0,m)
	elif i==50000: pa.sim_set_pot_parms(3.3,50.0,500.0,35.9*pi/180.0,20.0,1500.0,m)
	elif i==55000: pa.sim_set_pot_parms(3.3,50.0,500.0,35.9*pi/180.0,20.0,2500.0,m)

	if i%1000==0:
#		pa.sim_rescale()
		print pa.sim_potential()
		
		print i,") ",
		pa.sim_printstat()
		img=pa.pdb2mrc_by_summation(160,4.52,10.,-1)
		img.write_image("x.hdf",i/1000)
