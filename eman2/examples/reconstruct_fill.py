# proj_complete_test.py		Steven Ludtke	6/2013
# this program will generate various numbers of projections and then test completeness in Fourier space during reconstruction with a specified reconstructor

from EMAN2 import *

from sys import argv
from os import system
import random
import os

size=512
recon=Reconstructors.get("fourier",{"size":(size,size,size),"sym":"c1","savenorm":"test.hdf"})

ptcl=test_image(size=(size,size,1))
ptclf=recon.preprocess_slice(ptcl,Transform())

sym_object=parsesym("c1")

#nslices=[16,32,64,128,256,512,1024,2048,4096]
angs=[4,5,6,9,10,12,15,18,20,24,30,36,45]
for da in angs:
	eulers = sym_object.gen_orientations("eman",{"delta":90.0/da,"perturb":True})
	sl=len(eulers)
	
	recon.setup()
	
	for i,e in enumerate(eulers):
		print i,e
		recon.insert_slice(ptclf,e,1.0)
	
	final=recon.finish(True)

	img=EMData("test.hdf",0)
	os.rename ("test.hdf","test.{:1.1f}.hdf".format(90.0/da))
	
	r1=[0]*(size/2)
	r2=[0]*(size/2)
	n=[0]*(size/2)
	for z in xrange(-size/2,size/2):
		for y in xrange(-size/2,size/2):
			for x in xrange(size/2+1):
				r=int(sqrt(x*x+y*y+z*z))
				if r>=size/2 : continue
				n[r]+=1
				if img[x,y+size/2,z+size/2]>0.5 : r1[r]+=1
				if img[x,y+size/2,z+size/2]>1.5 : r2[r]+=1
				
	out=file("fill_ptrb_%04d.txt"%sl,"w")
	for i in xrange(size/2):
		out.write("%d\t%1.2f\t%1.2f\n"%(i,100.0*r1[i]/n[i],100.0*r2[i]/n[i]))
	out=None
				
