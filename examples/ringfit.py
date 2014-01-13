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
