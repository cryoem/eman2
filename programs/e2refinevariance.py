#!/usr/bin/env python

#
# Author: Steven Ludtke, 04/10/2003 (sludtke@bcm.edu)
# Copyright (c) 2000-2006 Baylor College of Medicine
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
from math import *
import os
import sys

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options] 
	This program will compute a variance map, using the results from an iteration of e2refine_easy.py as a basis
	for the computation. It uses a bootstrap resampling technique, where random particles are excluded from
	the model with replacement. This is repeated N times, producing a new 3-D model each time. The variance of
	the 3-D models is then computed. Note that this program requires a fair bit of memory. If running the entire program
	on a single machine, you will need enough memory to hold ~5 copies of a 3-D map + whatever is required
	by e2classaverage.py. It may be best-used on downsampled data (also for speed)."""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	
	# Updated to EMAN2.1 by Muyuan Chen, on Dec 8th, 2014
	
	#options associated with e2refinevariance.py
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	parser.add_argument("--input", dest="input", default=None,type=str, help="The name of the image containing the particle data")
	parser.add_argument("--usefilt", dest="usefilt", type=str,default=None, help="Specify a particle data file that has been low pass or Wiener filtered. Has a one to one correspondence with your particle data. If specified will be used in projection matching routines, and elsewhere.")
	parser.add_argument("--path", default=None, type=str,help="The name of a directory where results of e2refine_easy.py are placed. If unspecified will generate one automatically of type refine_??.")
	parser.add_argument("--output", default=None, type=str,help="The name of a directory where the variance calculated should be placed. If unspecified will generate one automatically of type refinevar_??.")
	parser.add_argument("--mass", default=None, type=float,help="The mass of the particle in kilodaltons, used to run normalize.bymass. If unspecified nothing happens. Requires the --apix argument.")
	parser.add_argument("--apix", default=None, type=float,help="The angstrom per pixel of the input particles. This argument is required if you specify the --mass argument. If unspecified, the convergence plot is generated using either the project apix, or if not an apix of 1.")
	parser.add_argument("--automask3d", default=None, type=str,help="The 5 parameters of the mask.auto3d processor, applied after 3D reconstruction. These paramaters are, in order, isosurface threshold,radius,nshells, ngaussshells and nmaxseed. From e2proc3d.py you could achieve the same thing using --process=mask.auto3d:threshold=1.1:radius=30:nmaxseed=10:nshells=5:ngaussshells=5.")
	parser.add_argument("--nmodels", type=int, help="The number of different bootstrap models to generate for the variance computation. Default=10", default=10)
	parser.add_argument("--iteration", type=int, help="The refinement iteration to use as a basis for the variance map", default=-1)
	parser.add_argument("--volfiles",action="store_true",help="This will bypass the construction of the individual resampled models, and use files previously generated with the --keep3d options",default=False)
	parser.add_argument("--threads", default=1, type=int, help="Number of threads to run in parallel on a single computer when multi-computer parallelism isn't useful")


	# options associated with e2project3d.py
	parser.add_argument("--sym", dest = "sym", help = "Specify symmetry - choices are: c<n>, d<n>, h<n>, tet, oct, icos",default='c1')
		
	# options associated with e2classaverage.py
	parser.add_argument("--classkeep",type=float,help="The fraction of particles to keep in each class, based on the similarity score generated by the --cmp argument.")
	parser.add_argument("--classkeepsig", default=False, action="store_true", help="Change the keep (\'--keep\') criterion from fraction-based to sigma-based.")
	parser.add_argument("--classiter", type=int, help="The number of iterations to perform. Default is 1.", default=1)
	parser.add_argument("--classalign",type=str,help="If doing more than one iteration, this is the name and parameters of the 'aligner' used to align particles to the previous class average.", default="rotate_translate_flip")
	parser.add_argument("--classaligncmp",type=str,help="This is the name and parameters of the comparitor used by the fist stage aligner  Default is dot.",default="phase")
	parser.add_argument("--classralign",type=str,help="The second stage aligner which refines the results of the first alignment in class averaging. Default is None.", default=None)
	parser.add_argument("--classraligncmp",type=str,help="The comparitor used by the second stage aligner in class averageing. Default is dot:normalize=1.",default="dot:normalize=1")
	parser.add_argument("--classaverager",type=str,help="The averager used to generate the class averages. Default is \'mean\'.",default="mean")
	parser.add_argument("--classcmp",type=str,help="The name and parameters of the comparitor used to generate similarity scores, when class averaging. Default is \'dot:normalize=1\'", default="dot:normalize=1")
	parser.add_argument("--classnormproc",type=str,default="normalize.edgemean",help="Normalization applied during class averaging")
	parser.add_argument("--classrefsf",default=False, action="store_true", help="Use the setsfref option in class averaging to produce better filtered averages.")
	parser.add_argument("--prefilt",action="store_true",help="Filter each reference (c) to match the power spectrum of each particle (r) before alignment and comparison",default=False)
	
	#options associated with e2make3d.py
	parser.add_argument("--pad", type=int, dest="pad", help="To reduce Fourier artifacts, the model is typically padded by ~25 percent - only applies to Fourier reconstruction", default=0)
	parser.add_argument("--recon", dest="recon", default="fourier", help="Reconstructor to use see e2help.py reconstructors -v")
	parser.add_argument("--m3dkeep", type=float, help="The percentage of slices to keep in e2make3d.py")
	parser.add_argument("--m3dkeepsig", default=False, action="store_true", help="The standard deviation alternative to the --m3dkeep argument")
	parser.add_argument("--m3dsetsf", default=False, action="store_true", help="The standard deviation alternative to the --m3dkeep argument")
	parser.add_argument("--m3diter", type=int, default=4, help="The number of times the 3D reconstruction should be iterated")
	parser.add_argument("--m3dpreprocess", type=str, default="normalize.edgemean", help="Normalization processor applied before 3D reconstruction")
	parser.add_argument("--m3dpostprocess", type=str, default=None, help="Post processor to be applied to the 3D volume once the reconstruction is completed")
	parser.add_argument("--shrink3d",type=int,help="Shrink the class-averages and make a downsampled variance map",default=0)
	parser.add_argument("--keep3d",action="store_true",help="Keep all of the individual 3-D models used to make the variance map. This make take substantial disk space.")

	#lowmem!
	parser.add_argument("--lowmem", default=False, action="store_true",help="Make limited use of memory when possible - useful on lower end machines")
	parser.add_argument("--parallel","-P",type=str,help="Run in parallel, specify type:<option>=<value>:<option>:<value>",default=None)
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	
	(options, args) = parser.parse_args()
	
	# read json file of the refinement for some parameters
	db=js_open_dict(options.path+"/0_refine_parms.json")
	options.apix=db["apix"]
	options.astep=db["astep"]
	
	if options.sym.lower()=="none" : options.sym="c1"
	
	if options.iteration<0 :
		print "You must specify a refinement iteration to use as a basis for the variance map."
		sys.exit(1)
	
	if options.automask3d: 
		vals = options.automask3d.split(",")
		mapping = ["threshold","radius","nshells","nshellsgauss","nmaxseed"]
		s = "mask.auto3d"
		try:
			for i,p in enumerate(mapping):
				s += ":"+p+"="+vals[i]
		except:
			print "Error: automask3d requires 5 parameters now. See the Wiki."
			sys.exit(1)
		s+= ":return_mask=1"
		options.automask3d = s

	logid=E2init(sys.argv,options.ppid)
	
	if options.output == None:
		fls=[int(i[-2:]) for i in os.listdir(".") if i[:10]=="refinevar_" and len(i)==12 and str.isdigit(i[-2:])]
		if len(fls)==0 : fls=[0]
		options.output = "refinevar_{:02d}".format(max(fls)+1)
	try: os.makedirs(options.output)
	except: pass
	
	# Combine even/odd class matrixes
	clsmx=EMData.read_images("%s/classmx_%02d_odd.hdf"%(options.path,options.iteration))
	clsmx_odd=EMData.read_images("%s/classmx_%02d_odd.hdf"%(options.path,options.iteration))
	clsmx_even=EMData.read_images("%s/classmx_%02d_even.hdf"%(options.path,options.iteration))
	
	for i in range(len(clsmx)):
		oddsize=clsmx_odd[i].get_ysize()
		evensize=clsmx_even[i].get_ysize()
		ysize=oddsize+evensize
		clsmx[i].set_size(1,ysize)
		for j in range(clsmx[i].get_ysize()):
			if(j%2==0): #even
				clsmx[i].set_value_at(0,j,clsmx_even[i].get_value_at(0,j/2))
			else: #odd
				clsmx[i].set_value_at(0,j,clsmx_odd[i].get_value_at(0,(j-1)/2))
			
		
		clsmx[i].write_image("%s/classify.hdf"%(options.output),i)

	
	
	# Generating projection from the combination of the even/odd set
	projectcmd = "e2project3d.py {path}/threed_{itrm1:02d}.hdf  --outfile {outpath}/projections.hdf -f --orientgen {orient} --sym {sym} --parallel thread:{threads} --verbose {verbose}".format(
		path=options.path,itrm1=options.iteration-1,itr=options.iteration,orient=db["orientgen"],sym=options.sym,threads=options.threads,verbose=options.verbose,outpath=options.output)
	launch_childprocess(projectcmd)
	
	nprogress=options.nmodels*2.0+1.0
	
	# this loops over each of the n models we create to compute the variance
	for mod in xrange(options.nmodels) :
		if not options.volfiles :
			if options.verbose : print "Class-averaging"
			# Compute class-averages with the --resample option
			
			options.classifyfile="%s/classify.hdf"%(options.output)
			options.projfile="%s/projections.hdf"%(options.output)
			options.cafile="%s/variance_classes_tmp.hdf"%(options.output)
			print get_classaverage_cmd(options)
			if ( launch_childprocess(get_classaverage_cmd(options)) != 0 ):
				print "Failed to execute %s" %get_classaverage_cmd(options)
				sys.exit(1)
			E2progress(logid,(mod*2.0+1)/nprogress)
			
			# deal with shrink3d
			if options.shrink3d : 
				if options.verbose : print "Shrinking"
				print "e2proc2d.py %s %s --meanshrink=%d --inplace --writejunk"%(options.cafile,"%s/variance_classes_tmp_s.hdf"%(options.output),options.shrink3d)
				if ( launch_childprocess("e2proc2d.py %s %s --meanshrink=%d --inplace --writejunk"%(options.cafile,"%s/variance_classes_tmp_s.hdf"%(options.output),options.shrink3d)) != 0 ):
					print "Failed to execute CA shrink"
					sys.exit(1)
				options.cafile="%s/variance_classes_tmp_s.hdf"%(options.output)
			
			if options.verbose : print "3-D Reconstruction"
			# build a new 3-D map
			options.model="%s/variance_threed_tmp.hdf"%(options.output)
			print get_make3d_cmd(options)
			if ( launch_childprocess(get_make3d_cmd(options)) != 0 ):
				print "Failed to execute %s" %get_make3d_cmd(options)
				sys.exit(1)
			E2progress(logid,(mod*2.0+2.0)/nprogress)

			# enforce symmetry
			if options.sym.lower()!="c1" :
				launch_childprocess("e2proc3d.py %s %s --sym=%s"%(options.model,options.model,options.sym))

			if options.verbose : print "Post-processing"

			cur_map=EMData(options.model,0)
			nx=cur_map["nx"]
			ny=cur_map["ny"]
			nz=cur_map["nz"]
			apix=cur_map["apix_x"]
			
			if options.mass:
				# if options.mass is not none, the check function has already ascertained that it's postivie non zero, and that the 
				# apix argument has been specified.
				cur_map.process_inplace("normalize.bymass",{"apix":options.apix, "mass":options.mass})
				if options.automask3d: 

					automask_parms = parsemodopt(options.automask3d) # this is just so we only ever have to do it

					mask = cur_map.process(automask_parms[0],automask_parms[1])
					cur_map.mult(mask)
			
			# Write the individual models to MRC files
			if options.keep3d and not options.volfiles : cur_map.write_image("%s/var3d_maps.hdf"%(options.output),mod)
	
		if options.volfiles : cur_map=EMData("%s/var3d_maps.hdf"%(options.output),mod)
	
		# now keep a sum of all of the maps and all of the maps^2
		if mod==0 :
			mean_map=cur_map
			sqr_map=cur_map.process("math.squared")
		else :
			mean_map.add(cur_map)
			cur_map.process_inplace("math.squared")
			sqr_map.add(cur_map)
			cur_map=None
		
	# Ok, all done, now compute the mean and standard deviation. The mean map should look nearly identical to
	# the original results from the same iteration
	mean_map.mult(1.0/options.nmodels)
	mean_map.write_image("%s/threed_mean.hdf"%(options.output),0)
	
	# Now compute the variance from the two maps
	sqr_map.mult(1.0/options.nmodels)			# pixels are mean of x^2
	mean_map.process_inplace("math.squared")	
	sqr_map.sub(mean_map)						
	
	### Symmetry downweighting
	weight=EMData(nx,ny,nz)
	weight.to_zero()
	weight["apix_x"]=apix
	weight["apix_y"]=apix
	weight["apix_z"]=apix
	
	w=weight.copy()
	# A line along Z
	for i in xrange(0,nz): w[nx/2,ny/2,i]=1.0

	# replicate the line under symmetry
	t=Transform()
	ns=t.get_nsym(options.sym)
	for i in xrange(ns):
		t2=t.get_sym(options.sym,i)
		wc=w.process("xform",{"transform":t2})		# transformed version of weight
		weight.add(wc)
	
	weight.process_inplace("threshold.clampminmax",{"minval":1.0,"maxval":500.0})	# 60 would really normally be the max here, but with helical symmetry supported...

	# This filters the duplication map similarly to the actual map, approximating the exaggeration in variance
	if options.m3dpostprocess :
		tv=weight[nx/2,ny/2,nz/4]
		(processorname, param_dict) = parsemodopt(options.m3dpostprocess)
		weight.process_inplace(processorname, param_dict)
	
		# Now this is a bit tricky to argue. Along the axis we have n-fold redundancy, so that will define our normalization
		# so, we  a point 1/2 way up the Z axis to its value before filtration. 1/2 way up Z is to deal with things like icosahedral symmetry reasonably
		if weight[nx/2,ny/2,nz/4]>1.0 : rescale=(tv-1.0)/(weight[nx/2,ny/2,nz/4]-1.0)
		else : rescale=1.0
#		print rescale
		
		# now we reprocess with the weighting factor, but don't reweight the 1.0 region...
		weight.add(-1.0)
		weight.mult(rescale)
		weight.add(1.0)
	
#	display(weight)
	weight.process_inplace("math.reciprocal",{"zero_to":1.0})		# this gives us 1/duplication on the symmetry axes and 1 off of the axes


#	display(weight)
	sqr_map.mult(weight)
	
	# Finally write the result
	sqr_map.write_image("%s/threed_variance.hdf"%(options.output),0)
	
	E2progress(logid,nprogress/nprogress)
	
	E2end(logid)

def get_make3d_cmd(options,check=False,nofilecheck=False):
	
	e2make3dcmd="e2make3dpar.py --input {inputfile} --sym {sym} --output {outfile} --apix {apix} --threads {threads} --fillangle {astep}".format(
	inputfile=options.cafile,  sym=options.sym, outfile=options.model,threads=options.threads, apix=options.apix, astep=options.astep)

	
	#e2make3dcmd += " --recon=%s --output=%s" %(options.recon,options.model)

	if str(options.m3dpreprocess) != "None":
		e2make3dcmd += " --preprocess=%s" %options.m3dpreprocess
		
	if str(options.m3dpostprocess) != "None":
		e2make3dcmd += " --postprocess=%s" %options.m3dpostprocess

	
	if (options.m3dkeep):
		e2make3dcmd += " --keep=%f" %options.m3dkeep
		if (options.m3dkeepsig): e2make3dcmd += " --keepsig"
	
	if options.m3dsetsf :
		e2make3dcmd += " --setsf=auto"
	
	if (options.lowmem): e2make3dcmd += " --lowmem"

	if (options.pad != 0):
		e2make3dcmd += " --pad=%d" %options.pad
		
	#if (options.verbose):
		#e2make3dcmd += " --verbose=" + str(options.verbose - 1)
	
	if ( check ):
		e2make3dcmd += " --check"	
			
	if ( nofilecheck ):
		e2make3dcmd += " --nofilecheck"
	
	return e2make3dcmd

def get_classaverage_cmd(options,check=False,nofilecheck=False):
	
	e2cacmd = "e2classaverage.py --resample --input=%s --classmx=%s --output=%s --storebad" %(options.input,options.classifyfile,options.cafile)
	
	e2cacmd += " --ref=%s --iter=%d -f --normproc=%s --averager=%s %s" %(options.projfile,options.classiter,options.classnormproc,options.classaverager,options.classrefsf)
	
	e2cacmd += " --dbpath=%s" %options.path
	
	if (options.classkeep):
		e2cacmd += " --keep=%f" %options.classkeep
		
	if (options.classkeepsig):
		e2cacmd += " --keepsig"
	
	if (options.classiter >= 1 ):
		e2cacmd += " --cmp=%s --align=%s --aligncmp=%s" %(options.classcmp,options.classalign,options.classaligncmp)

		if (options.classralign != None):
			e2cacmd += " --ralign=%s --raligncmp=%s" %(options.classralign,options.classraligncmp)
	
	if options.usefilt != None:
		e2cacmd += " --usefilt=%s" %options.usefilt
	
	#if (options.verbose):
		#e2cacmd += " --verbose=" + str(options.verbose - 1)
	
	if options.prefilt:
		e2cacmd += " --prefilt"
	
	if options.parallel: e2cacmd += " --parallel=%s" %options.parallel

		
	#lowmem becamoe the only supportable behaviour as of May 5th 2009
#	if (options.lowmem): e2cacmd += " --lowmem"
	
	if ( check ):
		e2cacmd += " --check"	
			
	if ( nofilecheck ):
		e2cacmd += " --nofilecheck"
	
	return e2cacmd

if __name__ == "__main__":
    main()
