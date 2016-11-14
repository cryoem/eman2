#!/usr/bin/env python

#
# Author: Steven Ludtke, 2/18/15 (sludtke@bcm.edu)
# Copyright (c) 2015 Baylor College of Medicine
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
import random
from math import *
from numpy import *
import os
import sys
from e2simmx import cmponetomany
from EMAN2jsondb import JSTask,jsonclasses

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options]
	This program will take a set of reference-free class-averages (or other projections) and generate a set of possible
	3-D initial models. Unlike e2initialmodel.py this program relies on particles having high symmetry, and uses this
	to generate a more direct solution. It is typically used with icosahedral symmetry, but could also work with other
	high symmetries."""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_header(name="initialmodelheader", help='Options below this label are specific to e2initialmodel', title="### e2initialmodel options ###", row=1, col=0, rowspan=1, colspan=3)
	parser.add_argument("--input", dest="input", default=None,type=str, help="This file should contain good class-averages to use in constructing the initial model", browser='EMBrowserWidget(withmodal=True,multiselect=False)', guitype='filebox', row=0, col=0, rowspan=1, colspan=3)
	parser.add_argument("--shrink", dest="shrink", type = int, default=0, help="Optionally shrink the input particles by an integer factor prior to reconstruction. Default=0, no shrinking", guitype='shrinkbox', row=2, col=2, rowspan=1, colspan=1)
	parser.add_argument("--sym", dest = "sym", help = "Specify symmetry - choices are: c<n>, d<n>, h<n>, tet, oct, icos",default="c1", guitype='symbox', row=4, col=0, rowspan=1, colspan=2)
	parser.add_argument("--orientgen",type=str, default="eman:delta=3.0:inc_mirror=0:perturb=0",help="The type of orientation generator. Default is eman:delta=3.0:inc_mirror=0. See e2help.py orientgens", guitype='strbox', expert=True, row=4, col=2, rowspan=1, colspan=1)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	parser.add_argument("--threads", default=1,type=int,help="Number of threads to run in parallel on a single computer when multi-computer parallelism isn't useful", guitype='intbox', row=24, col=2, rowspan=1, colspan=1, mode="refinement[4]")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)

	# Database Metadata storage
	#parser.add_argument("--dbls",type=str,default=None,help="data base list storage, used by the workflow. You can ignore this argument.")

	(options, args) = parser.parse_args()
	verbose=options.verbose

	try: ptcls=EMData.read_images(options.input)
	except:
		print "Error: bad input file"
		exit(1)

	for i in range(len(ptcls)):
		if options.shrink>1 :
			ptcls[i].process_inplace("math.meanshrink",{"n":options.shrink})
		ptcls[i].process_inplace("normalize.edgemean",{})
	if ptcls[0]["nx"]>160 : print "WARNING: using a large box size may be slow. Suggest trying --shrink="
	if not ptcls or len(ptcls)==0 : parser.error("Bad input file")
	boxsize=ptcls[0].get_xsize()
	apix=ptcls[0]["apix_x"]
	if verbose>0 : print "%d particles %dx%d"%(len(ptcls),boxsize,boxsize)
	print "Models will be %1.3f A/pix"%apix

	# pad the particles with zeroes
	padsize=good_size(boxsize*1.25)
	recon=Reconstructors.get("fourier",{"size":padsize,"sym":options.sym,"mode":"gauss_2","verbose":max(options.verbose-2,0)})
	ptclsf=[]			# particles that have been preprocessed for reconstruction
	for i in range(len(ptcls)):
		padded=ptcls[i].get_clip(Region((boxsize-padsize)/2,(boxsize-padsize)/2,padsize,padsize))
		ptclsf.append(recon.preprocess_slice(padded,Transform()))

	# angles to use for testing
	[og_name,og_args] = parsemodopt(options.orientgen)
	sym_object = parsesym(options.sym)
	orts = sym_object.gen_orientations(og_name,og_args)

	logid=E2init(sys.argv,options.ppid)
	results=[]

	try: os.mkdir("initial_models")
	except: pass

	curmap=EMData(boxsize,boxsize,boxsize)
	curmap.to_zero()
	
	#daz=2.0*360.0/(boxsize*pi)		# 1 pixel step at r/2
	recon=Reconstructors.get("fourier",{"size":(padsize,padsize,padsize),"sym":options.sym,"mode":"gauss_2","verbose":max(options.verbose-3,0)})
	daz=3.0							# matches default orientation generator, good enough for an initial model?
	eulers=[]
#	out=file("dbg.txt","w")
	allbest=[]
	for n in xrange(len(ptcls)):
		best=(1e100,None,None,None)
		best2=(1e100,None,None,None)
		if options.verbose : print "Particle/average: ",n
		for ort in orts:
			for phi in arange(0,359.9,daz):
				if options.verbose>1: print ort.get_rotation()["az"],ort.get_rotation()["alt"],phi, " : ",
				
				recon.setup()
				
				# we reconstruct the map using a single projection in some test orientation, with symmetry
				ortins=Transform({"type":"eman","az":ort.get_rotation()["az"],"alt":ort.get_rotation()["alt"],"phi":phi})
				recon.insert_slice(ptclsf[n],ortins,1.0)
				trymapf=recon.finish(False)
				trymap=trymapf.do_ift()
				trymap.process_inplace("xform.phaseorigin.tocenter")
				trymap=trymap.get_clip(Region((padsize-boxsize)/2,(padsize-boxsize)/2,(padsize-boxsize)/2,boxsize,boxsize,boxsize))
	
				# then we reproject that map to see how well it matches the original image
				# effectively this is just a very expensive way of doing self-common-lines until we add multiple projections
				proj=trymap.project("standard",{"transform":ortins})
#				display((proj,ptcls[n]))
#				sim=proj.cmp("optsub",ptcls[n])
				sim=proj.cmp("ccc",ptcls[n])
			
				if sim<best[0] :
					if mapcmp(trymap,best[2])<-0.7:						# this means we have a better version of 'best'
						best=(sim,ortins,trymap,trymapf,proj,ptcls[n])
					else:												# the better map doesn't look like the existing map
						best2=best
						best=(sim,ortins,trymap,trymapf,proj,ptcls[n])
				elif sim<best2[0]:
					if mapcmp(trymap,best[2])>-0.7:						# better than the current second, but not similar to the first
						best2=(sim,ortins,trymap,trymapf,proj,ptcls[n])
				
				if options.verbose>1 : print sim
#				out.write("{}\t{}\t# {},{}\n".format(phi,sim,ort.get_rotation()["az"],ort.get_rotation()["alt"]))
		if options.verbose: print best[:2]
		
		if (options.verbose>2) :
			#per particle results
			best[2].write_image("best.hdf",-1)
			best[4].write_image("bestcmp.hdf",-1)
			best[5].write_image("bestcmp.hdf",-1)
			if best2[2]!=None:
				best2[2].write_image("best.hdf",-1)
				best2[4].write_image("bestcmp.hdf",-1)
				best2[5].write_image("bestcmp.hdf",-1)

		allbest.append(best)
		if best2[2]!=None:
			allbest.append(best2)

	## similarity matrix among per particle results (2 maps/particle)
	#simmap=EMData(len(allbest),len(allbest),1)
	#simmap.to_one()
	#for i in xrange(1,len(allbest)):
		#for j in xrange(i):
			#c=mapcmp(allbest[i][2],allbest[j][2])
			#simmap[i,j]=c
			#simmap[j,i]=c
	
	#simmap.write_image("simmap.hdf")
	
	#bi,bj,bk=simmap.calc_min_location()
	#cursum=mapsum(allbest[bi][2],allbest[bj][2])
	#used=[bi,bj]

	cursum=min(allbest)[2]			# start with the best matching map we got
	used=[allbest.index(min(allbest))]
	if options.verbose : print "Best map: ",used[0]
	
	# we add in 1/3 more of the best matching volumes
	for n in xrange(len(allbest)/3):
		# Find the best match to the current sum and add it in
		best=(1.0,None)
		for i in xrange(len(allbest)):
			if i in used : continue
			best=min(best,(mapcmp(cursum,allbest[i][2]),allbest[i][2],i))
		
		if options.verbose : print "{}.  {}\t{}".format(i,best[2],best[0])
		used.append(best[2])
		cursum=mapsum(cursum,best[1])
	
	if options.verbose: print used
	cursum.process_inplace("normalize.edgemean")
	cursum.write_image("final.hdf")
	
	# write projection comparisons
	for i in used:
		allbest[i][-1].write_image("final_cmp.hdf",-1)
		cursum.project("standard",{"transform":allbest[i][1]}).process("normalize.edgemean").write_image("final_cmp.hdf",-1)
	
	E2end(logid)

def mapsum(m1,m2):
	"""Adds 2 maps taking into account handedness flips and 5-fold orientation uncertainty, picking the best solution for m2"""
	if m1==None or m2==None : return 2.0
	c1=m1.cmp("ccc",m2)
	mb=m2
	m3=m2.process("xform.flip",{"axis":"z"})
	c2=m1.cmp("ccc",m3)
	if c2<c1:
		mb=m3.copy()
		c1=c2
	m3.rotate(180.0,0,0)
	c2=m1.cmp("ccc",m3)
	if c2<c1:
		mb=m3.copy()
		c1=c2
	m3.process_inplace("xform.flip",{"axis":"z"})
	c2=m1.cmp("ccc",m3)
	if c2<c1:
		mb=m3.copy()
		c1=c2

	return m1+mb

def mapcmp(m1,m2):
	"""Compares 2 maps taking into account handedness flips and 5-fold orientation uncertainty"""
	if m1==None or m2==None : return 2.0
	c1=m1.cmp("ccc",m2)
	m3=m2.process("xform.flip",{"axis":"z"})
	c1=min(c1,m1.cmp("ccc",m3))
	m3.rotate(180.0,0,0)
	c1=min(c1,m1.cmp("ccc",m3))
	m3.process_inplace("xform.flip",{"axis":"z"})
	c1=min(c1,m1.cmp("ccc",m3))

	return c1

if __name__ == "__main__":
    main()

