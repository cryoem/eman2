#!/usr/bin/env python

#
# Author: Steven Ludtke, 10/21/2011 (sludtke@bcm.edu)
# Copyright (c) 2000- Baylor College of Medicine
#
# This software is issued under a joint BSD/GNU license. You may use the
# source code in this file under either license. However, note that the
# complete EMAN2 and SPARX software packages have some GPL dependencies,
# so you are responsible for compliance with the licenses of these packages
# if you opt to use BSD licensing. The warranty disclaimer below holds
# in either instance.
#
# This complete copyright notice must be included in any revised version of the
# source code. Additional authorship citations may be added, but existing
# author citations must be preserved.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  2111-1307 USA
#
#



from EMAN2 import *
from optparse import OptionParser
from math import *
import os
import sys
import time


def main():
	progname = os.path.basename(sys.argv[0])
	usage = """e2ptclzvssim.py [options] input1 input2 ...
	
Reads a full similarity matrix. Computes a Z score for each particle and plots vs actual
similarity score.

If there are N particles and M input files, the output file is:
an M*4 or M*6 column text file with N rows (lines). M*6 columns are present only if --refs is specified.
The columns are:

input1-Q input2-Q ... input1-QA input2-QA ... input1-Z input2-Z ... input1-B ... input1-alt ... input1-az

where:
Z - Z score
Q - raw quality value of best particle orientation
QA - Average quality for this particle across all measured references
B - reference number of best particle orientation
alt - altitude of best particle orientation
az - azimuth of best particle orientation

Note also that Z scores will not be accurate if the similarity matrix was computed using --twostage classification,
as not all elements are computed.
"""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

#	parser.add_option("--input",type=str,help="Similarity matrix to analyze",default=None)
	parser.add_argument("--refine",type=str,default=None,help="Automatically get parameters for a refine directory")
	parser.add_argument("--output",type=str,help="Output text file",default="zvssim.txt")
	parser.add_argument("--refs",type=str,help="Reference images from the similarity matrix (projections)",default=None)
	parser.add_argument("--inimgs",type=str,help="Input image file",default=None)
	parser.add_argument("--outimgs",type=str,help="Output image file",default="imgs.hdf")
	parser.add_argument("--filtimgs",type=str,help="A python expression using Z[n], Q[n] and N[n] for selecting specific particles to output. n is the 0 indexed number of the input file",default=None)
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	
	(options, args) = parser.parse_args()

		
	if options.refine!=None:
		if len(args)<1 :
			args.extend(sorted([ "bdb:%s#"%options.refine+i for i in db_list_dicts("bdb:"+options.refine) if "simmx" in i and len(i)==8]))
			print ", ".join(args)
		options.refs="bdb:%s#projections_%s"%(options.refine,args[-1][-2:])
		print "refs: ",options.refs

		db=db_open_dict("bdb:%s#register"%options.refine,True)
		options.inimgs=db["cmd_dict"]["input"]
		print "inimgs: ",options.inimgs

		
	if len(args)<1 : 
		print "Please specify input files"
		sys.exit(1)


	tmp=EMData(args[0],0,True)		# header info from first input, assume others are same
	nx=tmp["nx"]
	ny=tmp["ny"]

	logid=E2init(sys.argv,options.ppid)

# read in projection Euler angles
	if options.refs:
		ALTs=[]
		AZs=[]
		for i in xrange(nx):	
			# this reads the header, gets the orientation, and reads it out EMAN style
			ort=EMData(options.refs,i,True)["xform.projection"].get_rotation("eman")
			ALTs.append(ort["alt"])
			if ort["az"]>180.0 : az=ort["az"]-360.0		# we do this to reduce problems with symmetric structures with angle perturbation
			else: az=ort["az"]
			AZs.append(az)
		print nx," projections read"

	out=file(options.output,"w")
	outkey=file(options.output.rsplit(".",1)[0]+"_key.txt","w")

	colh=1
	t0=time.time()
	outkey.write( "0 - N\n")
	angcols=[]	# list of columns containing angles (used in image conversion at end)
	ncols=[]	# list of ptcl number cols

	# We read one line of the simmx at a time. The line represents all values
	# for a single particle
	for y in xrange(ny):
		if time.time()-t0>0.3 :
			print " %d/%d\r"%(y+1,ny),
			sys.stdout.flush()
			t0=time.time()
		Qs=[]		# Quality
		QAs=[]		# Average Quality for particle
		Q0s=[]		# Quality - Worst Quality
		Zs=[]		# Z score for classification
		Ns=[]		# classified best class
		DTAs=[]		# angle of translation (best orient)
		DTRs=[]		# translation distance
		DAs=[]		# rotation angle
		for cm in args:
			im=EMData(cm,0,False,Region(0,y,nx,1))
			im2=im.copy()										# here we make a copy with the max value-> zero. Used with 2-stage classification to get a better sigma
			im2.add(-im2["maximum"])
			try : Z=(im["mean"]-im["minimum"])/im2["sigma_nonzero"]
			except: Z=0
			Zs.append(Z)
			Q=im["minimum"]
			Qs.append(Q)
			QAs.append(im["mean"])
			Q0s.append(im2["minimum"])	# im2 has worst quality -> 0
			N=im.calc_min_index()
			Ns.append(N)
			
			# a little silly to read the whole line, but also a little silly to have a 1 pixel image :^/
			imdx=EMData(cm,1,False,Region(0,y,nx,1))
			imdy=EMData(cm,2,False,Region(0,y,nx,1))
			imda=EMData(cm,3,False,Region(0,y,nx,1))
			
			DTAs.append(atan2(imdy[N],imdx[N])*180.0/pi)
			DTRs.append(hypot(imdx[N],imdy[N]))
			DAs.append(imda[N])
			
		out.write("%d\t"%y)

		# note that we are multiplying qualities by -1 so BIGGER is better
		for i,q in enumerate(Qs) : 
			out.write("%1.4g\t"%-q)			# Quality of best alignment
			if y==0:
				outkey.write( "%d - (Qs) Best Ali Quality *-1 (%s)\n"%(colh,args[i]))
				colh+=1
		for i,qa in enumerate(QAs) : 
			out.write("%1.4g\t"%-qa)			#  Average quality for this particle
			if y==0:
				outkey.write( "%d - (QAs) Avg Quality *-1 (%s)\n"%(colh,args[i]))
				colh+=1
		for i,q0 in enumerate(Q0s) : 
			out.write("%1.4g\t"%-q0)			#  Average quality for this particle
			if y==0:
				outkey.write( "%d - (Q0s) Quality - Worst Quality *-1 (%s)\n"%(colh,args[i]))
				colh+=1
		for i,z in enumerate(Zs) : 
			out.write("%1.4g\t"%z)			# Z for this alignment
			if y==0:
				outkey.write( "%d - (Zs) Z score (%s)\n"%(colh,args[i]))
				colh+=1
		for i,n in enumerate(Ns) : 
			out.write("%d\t"%n)				# class number
			if y==0:
				outkey.write( "%d - (Ns) Best Ali Class (%s)\n"%(colh,args[i]))
				ncols.append(colh)
				colh+=1
		for i,dta in enumerate(DTAs) : 
			out.write("%1.4g\t"%dta)			# Trans align angle for best alignment
			if y==0: 
				outkey.write( "%d - (DTAs) Best Ali trans angle (%s)\n"%(colh,args[i]))
				angcols.append(colh)
				colh+=1
		for i,dtr in enumerate(DTRs) : 
			out.write("%1.4g\t"%dtr)			# Trans align dist for best alignment
			if y==0: 
				outkey.write( "%d - (DTRs) Best Ali trans dist (%s)\n"%(colh,args[i]))
				colh+=1
		for i,da in enumerate(DAs) : 
			out.write("%1.4g\t"%da)			# Rot align angle for best alignment
			if y==0: 
				outkey.write( "%d - (DAs) Best Ali rot angle (%s)\n"%(colh,args[i]))
				angcols.append(colh)
				colh+=1


		# if refs were provided we also write out Euler angle columns
		if options.refs :
			for i,n in enumerate(Ns) : 
				out.write("%1.5g\t"%ALTs[n])	# Altitude Euler for best alignment
				if y==0:
					outkey.write( "%d - (ALTs) Alt (%s)\n"%(colh,args[i]))
					angcols.append(colh)
					colh+=1
			for i,n in enumerate(Ns) : 
				out.write("%1.5g\t"%AZs[n])	# Azimuth Euler for best alignment
				if y==0:
					outkey.write( "%d - (AZs) Az (%s)\n"%(colh,args[i]))
					angcols.append(colh)
					colh+=1

		# if input values were provided, we write per-particle info
		if options.inimgs:
			ptcl=EMData(options.inimgs,y,True)		# read image header
			try :
				ctf=ptcl["ctf"]
				out.write("%1.4g\t%1.4g\t"%(ctf.defocus,ctf.bfactor))
				snr=ctf.snr
				snrw=[n*i for n,i in enumerate(snr)]
				nsnr=len(snr)
				# This gives integrated radial weighted SSNR over 3 resolution ranges
#				out.write("%1.3g\t%1.3g\t%1.3g\t"%(sum(snrw[1:nsnr/8])/(nsnr/8),sum(snrw[nsnr/16:nsnr/3])/(nsnr/3-nsnr/16),sum(snrw[nsnr/3:nsnr*2/3])/(nsnr*2/3-nsnr/3)));
				out.write("%1.3g\t%1.3g\t%1.3g\t"%(sum(snr[3:nsnr/6])/(nsnr/6-3),sum(snr[nsnr/8:nsnr/3])/(nsnr/3-nsnr/8),sum(snr[nsnr/3:nsnr*2/3])/(nsnr*2/3-nsnr/3)));
				if y==0:
					outkey.write( "%d - defocus\n"%(colh))
					colh+=1
					outkey.write( "%d - bfactor\n"%(colh))
					colbfac=colh
					colh+=1
					outkey.write( "%d - snr low\n"%(colh))
					colh+=1
					outkey.write( "%d - snr middle\n"%(colh))
					colh+=1
					outkey.write( "%d - snr high\n"%(colh))
					colh+=1
				
			except:
				pass								# no particle CTF info

		out.write("\n")

		if options.filtimgs!=None :
			e=eval(options.filtimgs)
			#print y,e,DTRs[0]-DTRs[1],fabs(DAs[0]-DAs[1]),fabs(ALTs[0]-ALTs[1]),Qs[1]
			if options.inimgs!=None and e:
				EMData(options.inimgs,y).write_image(options.outimgs,-1)

	print " %d/%d\n"%(ny,ny),
	print "Output in ",options.output
	print "Key in ",options.output.rsplit(".",1)[0]+"_key.txt"
	out.close()
	
	# we convert the output we just generated to an image file for other potential analysis techniques
	import numpy
	ary=numpy.loadtxt(options.output).transpose()[1:]		# load the entire text file, rotate so we can manipulate columns, and throw away the row number column
	for i in angcols: ary[i-1]*=pi/180.0					# convert angles to radians to better match scale of other parameters
	ary[colbfac-1]=numpy.sqrt(ary[colbfac-1])/100.0				# B-factor -> sqrt(B)/100.0
	print ncols
	print max(ary[ncols[0]-1])
	for i in ncols: ary[i-1]/=max(ary[i-1])					# class numbers -> 0-1 range
	print max(ary[ncols[0]-1])
	ary=ary.transpose()										# transpose again
	
	imgpath=options.output.rsplit(".",1)[0]+".hdf"
	print "Image stack output in ",imgpath
	for i in xrange(ary.shape[0]):							# we write the output to a stack of images
		im=EMNumPy.numpy2em(ary[i])
		im.write_image(imgpath,i)
	
	E2end(logid)


if __name__ == "__main__":  main()
