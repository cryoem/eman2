# 12/26/2013	Steven Ludtke
# This script is designed to quantitatively analyze DNA minicircles. Could be used for any other small objects as well
# It first filters and normalizes the data to try to isolate the particles under consideration
# Next it aligns the particles so the longest axis is along Z, second longest, on Y, and shortest on X. It does this by computing
# the eigenvectors of the Intertia matrix.

from EMAN2 import *
import os
import sys
from numpy import *
import numpy.linalg as LA

#ptcls=[i for i in os.listdir(".") if ".hdf" in i]

try: os.mkdir("xf")
except: pass

ldr=""
for pf in sorted(sys.argv[1:]):
	try: p=EMData(pf)
	except:
		print "Couldn't read:",pf
		continue

	# identify when we've started on a new directory
	dr=pf.split("/")[0]
	if ldr!=dr :
		out=file("stat_%s.txt"%dr,"w")
		ldr=dr
		drn=0
	else: drn+=1

	# This will (hopefully) isolate the mini-circle
	p.process_inplace("normalize.edgemean")
	p.process_inplace("mask.dust3d",{"voxels":8000,"threshold":1.6})
	p.process_inplace("mask.auto3d",{"nshells":2,"nshellsgauss":3,"radius":35,"threshold":1.6})
	p.process_inplace("threshold.belowtozero",{"minval":0})
	p.process_inplace("xform.centerofmass")
	comxf=p["xform.align3d"]
	p.process_inplace("normalize.unitsum")
	p.mult(10000.0)

	# compute the resulting inertia matrix
	an=Analyzers.get("inertiamatrix",{"verbose":0})
	an.insert_image(p)
	mxi=an.analyze()
	mx=EMNumPy.em2numpy(mxi[0])

	#print mx

	# Compute the eigenvalues/vectors
	eigvv=LA.eig(mx)		# a 3-vector with eigenvalues and a 3x3 with the vectors
	if min(eigvv[0])==0 :
		print "error on ",pf
		drn-=1
		continue
	#print eigvv[0]
	#print eigvv[1]
	eig=[(1.0/eigvv[0][i],eigvv[1][:,i]) for i in xrange(3)]  # extract for sorting
	#eig=sorted(eig,reverse=True)		# now eig is sorted in order from major to minor axes
	eig=sorted(eig)		# now eig is sorted in order from major to minor axes
	T=array([eig[0][1],eig[1][1],eig[2][1]])            # reassemble sorted matrix
	#print eig[0][0],eig[1][0],eig[2][0]
	#print T
	#print LA.inv(T)
	#print "============================"
	
	#out.write("%1.3g\t%1.3g\t%1.3g\t# %s\n"%(1.0/eig[0][0],1.0/eig[1][0],1.0/eig[2][0],pf.split("/")[-1]))
	#print "%1.3g\t%1.3g\t%1.3g\t# %s"%(1.0/eig[0][0],1.0/eig[1][0],1.0/eig[2][0],pf)
	
	T=Transform((float(i) for i in (eig[0][1][0],eig[0][1][1],eig[0][1][2],0,eig[1][1][0],eig[1][1][1],eig[1][1][2],0,eig[2][1][0],eig[2][1][1],eig[2][1][2],0)))
	#T=Transform((float(i) for i in (eig[0][1][0],eig[1][1][0],eig[2][1][0],0,eig[0][1][1],eig[1][1][1],eig[2][1][1],0,eig[0][1][2],eig[1][1][2],eig[2][1][2],0)))
	
	#T.printme()
	#print "------------------"
	p.transform(T)
	p.write_image("xf/"+dr+".hdf",drn)

	# write the original unfiltered, spherically masked, centered and oriented particles
	p2=EMData(pf)
	p2.process_inplace("normalize.edgemean")
	p2.transform(comxf)		# center
	p2.transform(T)			# reorient
	p2.process_inplace("mask.sharp",{"outer_radius":p2["nx"]/2-1})
	p2.write_image("xf/orig_"+dr+".hdf",drn)

	# now the shape is aligned to Z/Y/X so the greatest axial extent should be along Z
	an=Analyzers.get("shape",{"verbose":0})
	an.insert_image(p)
	shp=an.analyze()[0]
	
	# Z/Y - should always be >1, Y/X, Z/X
	out.write("%1.3g\t%1.3g\t%1.3g\t# %s\n"%(shp[2]/shp[1],shp[1]/shp[0],shp[2]/shp[0],pf.split("/")[-1]))
	print "%1.3g\t%1.3g\t%1.3g\t# %s"%(shp[2]/shp[1],shp[1]/shp[0],shp[2]/shp[0],pf)
	
