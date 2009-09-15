#!/usr/bin/env python

#
# Author: Steven Ludtke, 11/30/2008 (ludtke@bcm.edu)
# Copyright (c) 2000-2007 Baylor College of Medicine
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
import random
from math import *
import os
import sys
from e2simmx import cmponetomany
from e2make3d import fourier_reconstruction

class silly:
	"""This class exists so we can pass an object in to fourier_reconstruction"""
	pass

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """%prog [options] 
	Initial model generator"""
	parser = OptionParser(usage=usage,version=EMANVERSION)

	parser.add_option("--input", dest="input", default=None,type="string", help="The name of the image containing the particle data")
	parser.add_option("--iter", type = "int", default=8, help = "The total number of refinement iterations to perform")
	parser.add_option("--tries", type="int", default=1, help="The number of different initial models to generate in search of a good one")
	parser.add_option("--sym", dest = "sym", help = "Specify symmetry - choices are: c<n>, d<n>, h<n>, tet, oct, icos",default="c1")
	parser.add_option("--savemore",action="store_true",help="Will cause intermediate results to be written to flat files",default=False)
	parser.add_option("--verbose","-v", type="int", default=0,help="Verbosity of output (1-9)")
	parser.add_option("--orientgen",type="string", default="eman",help="The type of orientation generator. Default is safe. See e2help.py orientgens")

	# Database Metadata storage
	parser.add_option("--dbls",type="string",default=None,help="data base list storage, used by the workflow. You can ignore this argument.")
	
	(options, args) = parser.parse_args()
	verbose=options.verbose

	ptcls=EMData.read_images(options.input)
	for i in ptcls: i.process_inplace("normalize.edgemean",{})
	if not ptcls or len(ptcls)==0 : parser.error("Bad input file")
	boxsize=ptcls[0].get_xsize()
	if verbose : print "%d particles %dx%d"%(len(ptcls),boxsize,boxsize)

	# angles to use for refinement
	sym_object = parsesym(options.sym)
#	orts = sym_object.gen_orientations("eman",{"delta":7.5})
	#orts = sym_object.gen_orientations("rand",{"n":15})
	if options.orientgen == "rand":
		orts = sym_object.gen_orientations(options.orientgen,{"n":15})
	else:
   		orts = sym_object.gen_orientations(options.orientgen,{"delta":9.0})

	logid=E2init(sys.argv)
	results=[]
	
	try: os.mkdir("initial_models")
	except: pass
	dcts=db_list_dicts("bdb:initial_models")
	for ii in range(1,100):
		for jj in dcts: 
			if "model_%02d"%ii in jj : break
		else: break
	results_name="bdb:initial_models#model_%02d"%ii
	print results_name

	# We make one new reconstruction for each loop of t 
	for t in range(options.tries):
		if verbose: print "Try %d"%t
		threed=[make_random_map(boxsize)]		# initial model
		apply_sym(threed[0],options.sym)		# with the correct symmetry
		
		# This is the refinement loop
		for it in range(options.iter):
			E2progress(logid,(it+t*options.iter)/float(options.tries*options.iter))
			if verbose : print "Iteration %d"%it
			if options.savemore : threed[it].write_image("imdl.%02d.%02d.mrc"%(t,it))
			projs=[threed[it].project("standard",ort) for ort in orts]		# projections
			for i in projs : i.process_inplace("normalize.edgemean")
			if verbose>2: print "%d projections"%len(projs)
			
			# determine particle orientation
			bss=0.0
			bslst=[]
			for i in range(len(ptcls)):
				sim=cmponetomany(projs,ptcls[i],align=("rotate_translate_flip",{}),alicmp=("dot",{}),cmp=("frc",{}))
				bs=min(sim)
				#print bs[0]
				bss+=bs[0]
				bslst.append((bs[0],i))
				if verbose>2 : print "align %d \t(%1.3f)\t%1.3g"%(i,bs[0],bss)
				n=sim.index(bs)
				ptcls[i]["match_n"]=n
				ptcls[i]["match_qual"]=bs[0]
				ptcls[i]["xform.projection"]=orts[n]	# best orientation set in the original particle
			
			bslst.sort()					# sorted list of all particle qualities
			bslst.reverse()
			aptcls=[]
#			for i in range(len(ptcls)*3/4):		# We used to include 3/4 of the particles
			for i in range(len(ptcls)*7/8):
				n=ptcls[bslst[i][1]]["match_n"]
				aptcls.append(ptcls[bslst[i][1]].align("rotate_translate_flip",projs[n],{},"dot",{}))
				if it<2 : aptcls[-1].process_inplace("xform.centerofmass",{})
			
			bss/=len(ptcls)

			# 3-D reconstruction
			pad=(boxsize*3/2)
			pad-=pad%8
			recon=Reconstructors.get("fourier", {"quiet":True,"sym":options.sym,"x_in":pad,"y_in":pad})
			recon.setup()
			qual=0
			for ri in range(3):
				if ri>0 :
					for i,p in enumerate(aptcls):
						p2=p.get_clip(Region(-(pad-boxsize)/2,-(pad-boxsize)/2,pad,pad))
						recon.determine_slice_agreement(p2,p["xform.projection"],1)
						if ri==2 : qual+=recon.get_score(i)
				for p in aptcls:
					p2=p.get_clip(Region(-(pad-boxsize)/2,-(pad-boxsize)/2,pad,pad))
					recon.insert_slice(p2,p["xform.projection"])
			threed.append(recon.finish())
			
			if verbose>1 : print "Iter %d \t %1.4g (%1.4g)"%(it,bss,qual)
			
			threed[-1].process_inplace("xform.centerofmass")
			threed[-1]=threed[-1].get_clip(Region((pad-boxsize)/2,(pad-boxsize)/2,(pad-boxsize)/2,boxsize,boxsize,boxsize))
			threed[-1].process_inplace("mask.gaussian",{"inner_radius":boxsize/3.0,"outer_radius":boxsize/12.0})
			threed[-1].process_inplace("filter.lowpass.gauss",{"cutoff_abs":.2})
			threed[-1].process_inplace("normalize.edgemean")
			if it>1 : threed[-1].process_inplace("mask.auto3d",{"radius":boxsize/6,"threshold":2.5,"nshells":boxsize/20,"nshellsgauss":boxsize/20})
#			threed[-1]["quality"]=bss
			threed[-1]["quality"]=qual

#			threed[-1]["quality_projmatch"]=
#			display(threed[0])
#			display(threed[-1])
			#debugging output
			#for i in range(len(aptcls)):
				#projs[aptcls[i]["match_n"]].write_image("x.%d.hed"%t,i*2) 
				#aptcls[i].write_image("x.%d.hed"%t,i*2+1)
		#display(threed[-1])
		#threed[-1].write_image("x.mrc")
		if verbose : print "Model %d complete. Quality = %1.4f (%1.4f)"%(t,bss,qual)

		results.append((bss,threed[-1],aptcls,projs))
		results.sort()
		
		#dct=db_open_dict(results_name)
		#for i,j in enumerate(results): dct[i]=j[1]
		#dct.close()
		
		for i,j in enumerate(results):
			out_name = results_name+"_%02d"%(i+1)
			j[1].write_image(out_name,0)
			print out_name,j[1]["quality"],j[0]
			if options.dbls: # database list storage
				pdb = db_open_dict("bdb:project")
				old_data = pdb.get(options.dbls,dfl={})
				if isinstance(old_data,list): # in June 2009 we decided we'd transition the lists to dicts - so this little loop is for back compatibility
					d = {}
					for name in old_data:
						s = {}
						s["Original Data"] = name
						d[name] = s
					old_data = d
				
				s = {}
				s["Original Data"] = out_name
				old_data[out_name] = s
				pdb[options.dbls] = old_data
			for k,l in enumerate(j[3]): l.write_image(results_name+"_%02d_proj"%(i+1),k)	# set of projection images
			for k,l in enumerate(j[2]): 
				l.write_image(results_name+"_%02d_aptcl"%(i+1),k*2)						# set of aligned particles
				j[3][l["match_n"]].write_image(results_name+"_%02d_aptcl"%(i+1),k*2+1)	# set of projections matching aligned particles
			
#		threed[-1].write_image("x.%d.mrc"%t)
		
#		display(aptcls)
			
			
	E2end(logid)


def make_random_map(boxsize):
	"""This will make a map consisting of random noise, low-pass filtered and center-weighted for use
	as a random starting model in initial model generation. Note that the mask is eliptical and has random aspect."""
	
	ret=EMData(boxsize,boxsize,boxsize)
	ret.process_inplace("testimage.noise.gauss",{"mean":0.02,"sigma":1.0})
	ret.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.05})
	ret.process_inplace("xform.centerofmass",{})
#	ret.process_inplace("mask.gaussian",{"inner_radius":boxsize/3.0,"outer_radius":boxsize/12.0})
	ret.process_inplace("mask.gaussian.nonuniform",{"radius_x":boxsize/random.uniform(2.0,5.0),"radius_y":boxsize/random.uniform(2.0,5.0),"radius_z":boxsize/random.uniform(2.0,5.0)})
	
	return ret
	
def apply_sym(data,sym):
	"""applies a symmetry to a 3-D volume in-place"""
	xf = Transform()
	xf.to_identity()
	nsym=xf.get_nsym(sym)
	ref=data.copy()
	for i in range(1,nsym):
		dc=ref.copy()
		dc.transform(xf.get_sym(sym,i))
		data.add(dc)
	data.mult(1.0/nsym)	
	

if __name__ == "__main__":
    main()

