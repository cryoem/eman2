#!/usr/bin/env python
# 01/13/2014		Steven Ludtke
# This program tries to fit a closed loop of blobs to DNA-minicircle density pattern using
# a simple distance, angle, dihedral potential with a closed linear chain of balls

from EMAN2 import *
import sys
import numpy.linalg as LA
import random
import math
from numpy import *
import os
from os import listdir
from os.path import isfile, join

def main():
	
	
	usage=""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--output", type=str,help="output file name", default="circlestat")
	(options, args) = parser.parse_args()
	
	mypath=args[0]
	if mypath.endswith(".pdb") or mypath.endswith(".mrc"):
		process_image(mypath,"")
	else:
		filetype=args[1]
		
		if filetype=="all":
			mrcf = sorted([ f for f in listdir(mypath) if f.endswith("finalimg.mrc")])
			pdbf = sorted([ f for f in listdir(mypath) if f.endswith("_result.pdb")])
			#print mrcf
			
			for i in range(len(mrcf)):
				print mrcf[i],pdbf[i]
				pdbshp=process_image(join(mypath,pdbf[i]),pdbf[i])
				mrcshp=process_image(join(mypath,mrcf[i]),mrcf[i])
				outfile=open(options.output,"a")
				outfile.write("%s\t%1.3g\t%1.3g\t%1.3g\t%1.5g\t"%(pdbf[i],pdbshp[0],pdbshp[1],pdbshp[2],pdbshp[3]))
				outfile.write("||\t%1.3g\t%1.3g\t%1.3g\t%1.5g\n"%(mrcshp[0],mrcshp[1],mrcshp[2],mrcshp[3]))
				#outfile.write("%s\t%1.5g\n"%(fname,totalen))
				outfile.close()
			
			
		else:			
			if filetype=="mrc":
				files = sorted([ f for f in listdir(mypath) if f.endswith("finalimg.mrc")])
			elif filetype=="pdb":
				files = sorted([ f for f in listdir(mypath) if f.endswith("_result.pdb")])
			else:
				print "pdb or mrc only.."
				exit()
			
			print files
			
			for fname in files:
				print join(mypath,fname)
				shape=process_image(join(mypath,fname),fname)
				#totalen=process_image(join(mypath,fname),fname)
				
				#outfile=open("circlestat_m","a")
				outfile=open(options.output,"a")
				outfile.write("%s"%(fname))
				for s in shape:
					outfile.write("\t%f"%(s))
				outfile.write("\n")
				#outfile.write("%s\t%1.5g\n"%(fname,totalen))
				outfile.close()

	
def process_image(imgname,imgprefix):
	#compute the resulting inertia matrix
	
	if imgname.endswith(".pdb"):
		ispdb=1
	elif imgname.endswith(".mrc"):
		ispdb=0
		
	if ispdb:
		pa=PointArray()
		pa.read_from_pdb(imgname)
		totalen=pa.calc_total_length()
		pa.center_to_zero()
		#pa.sim_add_point_double()
		#pa.sim_add_point_double()
		#pa.sim_add_point_double()
		nn=pa.get_number_points()
		p=empty((nn,3),float)
		for i in range(nn):
			p[i,:]=array(pa.get_vector_at(i))
		p=delete(p,0,axis=0)
		tp=transpose(p)
		mx=dot(tp,p)
		### align the polygon 
		eigvv=LA.eig(mx)		# a 3-vector with eigenvalues and a 3x3 with the vectors
		eig=[(1.0/eigvv[0][i],eigvv[1][:,i]) for i in xrange(3)]  # extract for sorting
		eig=sorted(eig)		# now eig is sorted in order from major to minor axes
		et=matrix([[eig[0][1][0],eig[0][1][1],eig[0][1][2]],[eig[1][1][0],eig[1][1][1],eig[1][1][2]],[eig[2][1][0],eig[2][1][1],eig[2][1][2]]])
		
		pl=dot(p,et.T)
		shp=(pl.max(0)-pl.min(0))
		
		
		
		mn=mean(pl,axis=0)
		for j in range(20):
			pl[j]-=mn
		
		#print pl
		area=0
		for ip in range(len(pl)):
			### x1*y2-x2*y1
			p1= pl[ip].A1
			p2=pl[(ip+1)%20].A1
			area+=(p1[1]*p2[2]-p2[1]*p1[2])
			#print area

		area=abs(area/2)
		#print area
		
		pmax=argmax(abs(pl),axis=0)
		pm=pmax.A1[0]
		#print pm[0]
		d=0
		for j in range(10):
			n1=(pm+j)%20
			n2=(pm-j)%20
			p1=pl[n1]
			p2=pl[n2]
			t=linalg.norm(p1[0]-p2[0])
			d+=t
		d=d/10
		#shp=append(shp.A1,d)
		shp=append(shp.A1,totalen)
		shp=append(shp,area)
		print shp
		#for i in range(nn):
			#pa.set_vector_at(i,Vec3f(pl[i,0],pl[i,1],pl[i,2]),1.0)
		#pa.save_to_pdb(imgprefix+"aaa.pdb")
		return shp
		
		
	else:
		
		finalimg=EMData(imgname,0)
		finalimg.process_inplace("normalize.edgemean")
		finalimg.process_inplace("threshold.belowtozero",{"minval":.5})
		mnz=finalimg["mean_nonzero"]
		finalimg.process_inplace("threshold.belowtozero",{"minval":mnz})
		#finalimg.process_inplace("normalize.edgemean")
		#wtname="m_"+imgprefix
		#finalimg.write_image(wtname,0)
		finalimg.process_inplace("xform.centerofmass",{"threshold":0})
		
		an=Analyzers.get("inertiamatrix",{"verbose":0})
		an.insert_image(finalimg)
		mxi=an.analyze()
		mx=EMNumPy.em2numpy(mxi[0])
		
		# Compute the eigenvalues/vectors
		eigvv=LA.eig(mx)		# a 3-vector with eigenvalues and a 3x3 with the vectors
		eig=[(1.0/eigvv[0][i],eigvv[1][:,i]) for i in xrange(3)]  # extract for sorting
		eig=sorted(eig)		# now eig is sorted in order from major to minor axes
		#T=array([eig[0][1],eig[1][1],eig[2][1]])            # reassemble sorted matrix

		T=Transform((float(i) for i in (eig[0][1][0],eig[0][1][1],eig[0][1][2],0,eig[1][1][0],eig[1][1][1],eig[1][1][2],0,eig[2][1][0],eig[2][1][1],eig[2][1][2],0)))
		
		
		
		finalimg.transform(T)
		#nm=imgprefix+"ttt.mrc"
		#finalimg.write_image(nm,0)
		an=Analyzers.get("shape",{"verbose":0})
		an.insert_image(finalimg)
		shp=an.analyze()[0]
		#shp=EMNumPy.em2numpy(shp)
		
		print shp[0],shp[1],shp[2],shp[3],shp[2]/shp[1],shp[1]/shp[0]
		for i in range(4):
			shp[i]=sqrt(shp[i]/finalimg["mean"])*finalimg["apix_x"]
		
		return shp
		
		##finalimg.transform(T)
		##finalimg.write_image("trans_finalimg.mrc",0)
		## now the shape is aligned to Z/Y/X so the greatest axial extent should be along Z	
		#SX=finalimg.get_xsize()
		#SY=finalimg.get_ysize()
		#SZ=finalimg.get_zsize()
		#xsta=SX;xend=-10000;ysta=SY;yend=-10000;zsta=SZ;zend=-10000
		#print finalimg["mean"],finalimg["sigma"]
		#for i in range(SX):
			#for j in range(SY):
				#for k in range(SZ):
					#if(finalimg.get_value_at(i,j,k)>5):
						#v=T.transform(i,j,k)
						#x=v[0]
						#y=v[1]
						#z=v[2]
						#if(x<xsta): xsta=x
						#if(x>xend): xend=x
						#if(y<ysta): ysta=y
						#if(y>yend): yend=y
						#if(z<zsta): zsta=z
						#if(z>zend): zend=z
							
		##print xsta,xend,ysta,yend,zsta,zend
		#shp=[float(xend-xsta),float(yend-ysta),float(zend-zsta)]
		## Z/Y - should always be >1, Y/X, Z/X
		##out.write("%1.3g\t%1.3g\t%1.3g\t# %s\n"%(shp[2]/shp[1],shp[1]/shp[0],shp[2]/shp[0],pf.split("/")[-1]))
		#shape=sorted([abs(shp[0]),abs(shp[1]),abs(shp[2])])
		
		#print shape
		#print "%1.3g\t%1.3g\t%1.3g\t#"%(shape[2]/shape[1],shape[1]/shape[0],shape[2]/shape[0])
		#finalimg.transform(T)
		#nm=imgprefix+"ttt.mrc"
		#finalimg.write_image(nm,0)
		#return shape,totalen


if __name__ == '__main__':
    main() 
