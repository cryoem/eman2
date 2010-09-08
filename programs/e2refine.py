#!/usr/bin/env python

#
# Author: David Woolford, 10/19/2007 (woolford@bcm.edu)
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
from math import *
import os
import sys

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """%prog [options] 
	
	Single particle reconstruction refinement program. This is the main program used to perform
	iterative single-model single particle reconstruction in EMAN2. It has MANY options, many
	of which are passed on to other programs called as part of this process. For more information
	on the parameters and using this program, suggest reading the tutorial or using the
	e2workflow.py interface.
"""
	parser = OptionParser(usage=usage,version=EMANVERSION)
		
	#options associated with e2refine.py
	parser.add_option("--iter", dest = "iter", type = "int", default=0, help = "The total number of refinement iterations to perform")
	parser.add_option("--check", "-c", dest="check", default=False, action="store_true",help="Checks the contents of the current directory to verify that e2refine.py command will work - checks for the existence of the necessary starting files and checks their dimensions. Performs no work ")
	parser.add_option("--verbose", "-v", dest="verbose", action="store", metavar="n", type="int", default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	parser.add_option("--nomirror", dest="nomirror", default=False, action="store_true",help="Turn projection over the mirror portion of the asymmetric unit off")
	parser.add_option("--input", dest="input", default=None,type="string", help="The name of the image containing the particle data")
	parser.add_option("--model", dest="model", type="string",default="threed.0a.mrc", help="The name 3D image that will seed the refinement")
	parser.add_option("--usefilt", dest="usefilt", type="string",default=None, help="Specify a particle data file that has been low pass or Wiener filtered. Has a one to one correspondence with your particle data. If specified will be used in projection matching routines, and elsewhere.")
	parser.add_option("--path", default=None, type="string",help="The name of a directory where results are placed. If unspecified will generate one automatically of type refine_??.")
	parser.add_option("--mass", default=None, type="float",help="The mass of the particle in kilodaltons, used to run normalize.bymass. If unspecified nothing happens. Requires the --apix argument.")
	parser.add_option("--apix", default=None, type="float",help="The angstrom per pixel of the input particles. This argument is required if you specify the --mass argument. If unspecified, the convergence plot is generated using either the project apix, or if not an apix of 1.")
	parser.add_option("--automask3d", default=None, type="string",help="The 5 parameters of the mask.auto3d processor, applied after 3D reconstruction. These paramaters are, in order, isosurface threshold,radius,nshells and ngaussshells. From e2proc3d.py you could achieve the same thing using --process=mask.auto3d:threshold=1.1:radius=30:nshells=5:ngaussshells=5.")

	# options associated with e2project3d.py
	parser.add_option("--sym", dest = "sym", help = "Specify symmetry - choices are: c<n>, d<n>, h<n>, tet, oct, icos")
	parser.add_option("--projector", dest = "projector", default = "standard",help = "Projector to use")
	parser.add_option("--orientgen", type = "string",help = "The orientation generation argument for e2project3d.py")
		
	# options associated with e2simmx.py
	parser.add_option("--simalign",type="string",help="The name of an 'aligner' to use prior to comparing the images", default="rotate_translate_flip")
	parser.add_option("--simaligncmp",type="string",help="Name of the aligner along with its construction arguments",default="dot")
	parser.add_option("--simralign",type="string",help="The name and parameters of the second stage aligner which refines the results of the first alignment", default=None)
	parser.add_option("--simraligncmp",type="string",help="The name and parameters of the comparitor used by the second stage aligner. Default is dot.",default="dot")
	parser.add_option("--simcmp",type="string",help="The name of a 'cmp' to be used in comparing the aligned images", default="dot:normalize=1")
	parser.add_option("--simmask",type="string",help="A file containing a single 0/1 image to apply as a mask before comparison but after alignment", default=None)
	parser.add_option("--shrink", dest="shrink", type = "int", default=None, help="Optionally shrink the input particles by an integer amount prior to computing similarity scores. For speed purposes.")
	parser.add_option("--twostage", dest="twostage", type = "int", help="Optionally run a faster 2-stage similarity matrix, ~5-10x faster, generally same accuracy. Value specifies shrink factor for first stage, typ 1-3",default=0)
	parser.add_option("--prefilt",action="store_true",help="Filter each reference (c) to match the power spectrum of each particle (r) before alignment and comparison",default=False)
	
	# options associated with e2classify.py
	parser.add_option("--sep", type="int", help="The number of classes a particle can contribute towards (default is 1)", default=1)
	
	# options associated with e2classaverage.py
	parser.add_option("--classkeep",type="float",help="The fraction of particles to keep in each class, based on the similarity score generated by the --cmp argument.")
	parser.add_option("--classkeepsig", default=False, action="store_true", help="Change the keep (\'--keep\') criterion from fraction-based to sigma-based.")
	parser.add_option("--classiter", type="int", help="The number of iterations to perform. Default is 1.", default=3)
	parser.add_option("--classalign",type="string",help="If doing more than one iteration, this is the name and parameters of the 'aligner' used to align particles to the previous class average.", default="rotate_translate_flip")
	parser.add_option("--classaligncmp",type="string",help="This is the name and parameters of the comparitor used by the fist stage aligner  Default is dot.",default="phase")
	parser.add_option("--classralign",type="string",help="The second stage aligner which refines the results of the first alignment in class averaging. Default is None.", default=None)
	parser.add_option("--classraligncmp",type="string",help="The comparitor used by the second stage aligner in class averageing. Default is dot:normalize=1.",default="dot:normalize=1")
	parser.add_option("--classaverager",type="string",help="The averager used to generate the class averages. Default is \'mean\'.",default="mean")
	parser.add_option("--classcmp",type="string",help="The name and parameters of the comparitor used to generate similarity scores, when class averaging. Default is \'dot:normalize=1\'", default="dot:normalize=1")
	parser.add_option("--classnormproc",type="string",default="normalize.edgemean",help="Normalization applied during class averaging")
	parser.add_option("--classrefsf",default=False, action="store_true", help="Use the setsfref option in class averaging to produce better filtered averages.")
	parser.add_option("--classautomask",default=False, action="store_true", help="This will apply an automask to the class-average during iterative alignment for better accuracy. The final class averages are unmasked.")
	
	
	#options associated with e2make3d.py
	parser.add_option("--pad", type=int, dest="pad", help="To reduce Fourier artifacts, the model is typically padded by ~25% - only applies to Fourier reconstruction", default=0)
	parser.add_option("--recon", dest="recon", default="fourier", help="Reconstructor to use see e2help.py reconstructors -v")
	parser.add_option("--m3dkeep", type=float, help="The percentage of slices to keep in e2make3d.py")
	parser.add_option("--m3dkeepsig", default=False, action="store_true", help="The standard deviation alternative to the --m3dkeep argument")
	parser.add_option("--m3dsetsf", default=False, action="store_true", help="The standard deviation alternative to the --m3dkeep argument")
	parser.add_option("--m3diter", type=int, default=4, help="The number of times the 3D reconstruction should be iterated")
	parser.add_option("--m3dpreprocess", type="string", default="normalize.edgemean", help="Normalization processor applied before 3D reconstruction")
	parser.add_option("--m3dpostprocess", type="string", default=None, help="Post processor to be applied to the 3D volume once the reconstruction is completed")

	#lowmem!
	parser.add_option("--lowmem", default=False, action="store_true",help="Make limited use of memory when possible - useful on lower end machines")
	parser.add_option("--parallel","-P",type="string",help="Run in parallel, specify type:<option>=<value>:<option>:<value>",default=None)

	
	(options, args) = parser.parse_args()
	error = False
	if check(options,True) == True : 
		error = True
	if check_projection_args(options) == True : 
		error = True
	if check_simmx_args(options,True) == True :
		error = True
	if check_classify_args(options,True) == True :
		error = True
	if check_classaverage_args(options,True) == True :
		error = True
#	if check_make3d_args(options,True) == True:
#		error = True
	
	logid=E2init(sys.argv)
	
	
	if options.classrefsf : options.classrefsf=" --setsfref"
	else : options.classrefsf=" "
	
	if options.classautomask : options.classautomask=" --automask"
	else: options.classautomask=" "
	
	if error:
		print "Error encountered while checking command line, bailing"
		exit_refine(1,logid)
	
	if (options.check):
		exit_refine(0,logid)
	
	if options.path == None:
		options.path = numbered_path("refine",True)

	# store the inputs arguments forever in the refinement directory
	db = db_open_dict("bdb:"+options.path+"#register")
	db["cmd_dict"] = options.__dict__
		
	# this is the main refinement loop
	
	progress = 0.0
	total_procs = 5*options.iter
	
	if options.automask3d: automask_parms = parsemodopt(options.automask3d) # this is just so we only ever have to do it
	
	apix = get_apix_used(options)
	
	for i in range(0,options.iter) :
		
		number_options_file(i,"projections",options,"projfile")
		if ( os.system(get_projection_cmd(options)) != 0 ):
			print "Failed to execute %s" %get_projection_cmd(options)
			exit_refine(1,logid)
		progress += 1.0
		E2progress(logid,progress/total_procs)
		
		number_options_file(i,"simmx",options,"simmxfile")
		if options.twostage>0 :
			number_options_file(i,"proj_simmx",options,"simmxprojfile")
			number_options_file(i,"proj_stg1",options,"projstg1file")
			number_options_file(i,"simmx_stg1",options,"simmxstg1file")
			
		if ( os.system(get_simmx_cmd(options)) != 0 ):
			print "Failed to execute %s" %get_simmx_cmd(options)
			exit_refine(1,logid)
		progress += 1.0
		E2progress(logid,progress/total_procs)
			
		number_options_file(i,"classify",options,"classifyfile")
		if ( os.system(get_classify_cmd(options)) != 0 ):
			print "Failed to execute %s" %get_classify_cmd(options)
			exit_refine(1,logid)
		progress += 1.0
		E2progress(logid,progress/total_procs)
			
		number_options_file(i,"classes",options,"cafile")
		number_options_file(i,"cls_result",options,"resultfile")
		if ( os.system(get_classaverage_cmd(options)) != 0 ):
			print "Failed to execute %s" %get_classaverage_cmd(options)
			exit_refine(1,logid)
		progress += 1.0
		E2progress(logid,progress/total_procs)
			
		try : previous_model = options.filtered_model
		except : previous_model = options.model
		number_options_file(i,"threed",options,"model")
		new_model = options.model
		if ( os.system(get_make3d_cmd(options)) != 0 ):
			print "Failed to execute %s" %get_make3d_cmd(options)
			exit_refine(1,logid)
		progress += 1.0
		
		
		a = EMData(previous_model,0)
		b = EMData(new_model,0)
		
		if options.mass:
			# if options.mass is not none, the check function has already ascertained that it's postivie non zero, and that the 
			# apix argument has been specified.
			b.process_inplace("normalize.bymass",{"apix":options.apix, "mass":options.mass})
			if options.automask3d:
				number_options_file(i,"threed_mask",options,"mask_image")
				mask = b.process(automask_parms[0],automask_parms[1])
				mask.write_image(options.mask_image,0)
				db_close_dict(options.mask_image)
				b.mult(mask)
			
			number_options_file(i,"threed_filt",options,"filtered_model")
			b.write_image(options.filtered_model,0)
			db_close_dict(options.filtered_model)
		
		fsc = a.calc_fourier_shell_correlation(b)
		third = len(fsc)/3
		xaxis = fsc[0:third]
		plot = fsc[third:2*third]
		error = fsc[2*third:]
		
		if i == 0:
			s = "init_00"
		else:
			s1 = str(i-1)
			s2 = str(i)
			if len(s1) == 1: s1 = "0" + s1
			if len(s2) == 1: s2 = "0" + s2
			s = s1 + "_" + s2
		
		convergence_db_name = "bdb:"+options.path+"#convergence.results"
		db = db_open_dict(convergence_db_name)
		
		tmpaxis = [x/apix for x in xaxis]
		db[s+"_fsc"] = [tmpaxis,plot]
		#db["error_"+s+"_fsc"] = [xaxis,error] #we're not plotting the errors
		db_close_dict(convergence_db_name)
		
		E2progress(logid,progress/total_procs)
	E2end(logid)

def get_apix_used(options):
	'''
	Just an encapsulation apix retrieval
	Basically, if the apix argument is in the options, that's what you get
	Else the project db is checked for the global.apix parameter
	Else you just get 1
	'''
	apix = 1.0
	if options.apix: apix = options.apix # check function checks whether or not this value is positive, non zero
	elif db_check_dict("bdb:project"):
		# this is a temporary workaround to get things working in the workflow
		db2 = db_open_dict("bdb:project")
		apix = db2.get("global.apix",dfl=1)
		db_close_dict("bdb:project")
		if apix == 0: apix = 1
	#else apix is 1 !
	return apix
		
def number_options_file(i,file,options,attr):
	name = "bdb:"+options.path+"#" + file+"_"
	if i < 10:
		name += "0"
	
	name += str(i)
	setattr(options,attr,name)
	
def exit_refine(n,logid):
	E2end(logid)
	exit(n)

def get_make3d_cmd(options,check=False,nofilecheck=False):
	e2make3dcmd = "e2make3d.py --input=%s --sym=%s --iter=%d -f" %(options.cafile,options.sym,options.m3diter)
	
	e2make3dcmd += " --recon=%s --output=%s" %(options.recon,options.model)

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
		
	if (options.verbose):
		e2make3dcmd += " --verbose=" + str(options.verbose - 1)
	
	if ( check ):
		e2make3dcmd += " --check"	
			
	if ( nofilecheck ):
		e2make3dcmd += " --nofilecheck"
	
	return e2make3dcmd

def check_make3d_args(options, nofilecheck=False):
	
	cmd = get_make3d_cmd(options,True,nofilecheck)
	print ""
	print "#### Test executing make3d command: %s" %cmd
	return ( os.system(cmd) != 0)

def get_classaverage_cmd(options,check=False,nofilecheck=False):
	
	e2cacmd = "e2classaverage.py --input=%s --classmx=%s --storebad --output=%s" %(options.input,options.classifyfile,options.cafile)
	
	e2cacmd += " --ref=%s --iter=%d -f --resultmx=%s --normproc=%s --averager=%s %s %s" %(options.projfile,options.classiter,options.resultfile,options.classnormproc,options.classaverager,options.classrefsf,options.classautomask)
	
	e2cacmd += " --idxcache --dbpath=%s" %options.path
	
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
	
	if (options.verbose):
		e2cacmd += " --verbose=" + str(options.verbose - 1)
	
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

def check_classaverage_args(options, nofilecheck=False):
	if not hasattr(options,"cafile"): setattr(options,"cafile","dummy")
	if not hasattr(options,"resultfile"): setattr(options,"resultfile","dummy")
	cmd = get_classaverage_cmd(options,True,nofilecheck)
	print ""
	print "#### Test executing classaverage command: %s" %cmd
	return ( os.system(cmd) != 0)

def get_classify_cmd(options,check=False,nofilecheck=False):
	e2classifycmd = "e2classify.py %s %s --sep=%d -f" %(options.simmxfile,options.classifyfile,options.sep)
	
	if (options.verbose):
		e2classifycmd += " --verbose=" + str(options.verbose - 1)
	
	if ( check ):
		e2classifycmd += " --check"	
			
	if ( nofilecheck ):
		e2classifycmd += " --nofilecheck"
	
	return e2classifycmd

def check_classify_args(options, nofilecheck=False):
	if not hasattr(options,"classifyfile"): setattr(options,"classifyfile","dummy")
	cmd = get_classify_cmd(options,True,nofilecheck)
	print ""
	print "#### Test executing classify command: %s" %cmd
	return ( os.system(cmd) != 0)

def get_simmx_cmd(options,check=False,nofilecheck=False):
	
	if options.usefilt != None:
		image = options.usefilt
	else:
		image = options.input
	
	if options.twostage>0 : 
		try : e2simmxcmd = "e2simmx2stage.py %s %s %s %s %s %s -f --saveali --cmp=%s --align=%s --aligncmp=%s --shrinks1=%d"  %(options.projfile, image,options.simmxfile,options.simmxprojfile,options.projstg1file,options.simmxstg1file,options.simcmp,options.simalign,options.simaligncmp,options.twostage)
		except: print options
	else :
		e2simmxcmd = "e2simmx.py %s %s %s -f --saveali --cmp=%s --align=%s --aligncmp=%s"  %(options.projfile, image,options.simmxfile,options.simcmp,options.simalign,options.simaligncmp)

	if options.prefilt : e2simmxcmd+=" --prefilt"


	if options.simmask!=None : e2simmxcmd += " --mask=%s"%options.simmask
	
	if ( options.simralign != None ):
		e2simmxcmd += " --ralign=%s --raligncmp=%s" %(options.simralign,options.simraligncmp)
	
	if (options.verbose):
		e2simmxcmd += " --verbose=" + str(options.verbose - 1)
	
	if options.parallel: e2simmxcmd += " --parallel=%s" %options.parallel
	
	#if (options.lowmem): e2simmxcmd += " --lowmem"	
	
	if (options.shrink):
		e2simmxcmd += " --shrink="+str(options.shrink)
		
	if ( check ):
		e2simmxcmd += " --check"	
			
	if ( nofilecheck ):
		e2simmxcmd += " --nofilecheck"
		
	
	return e2simmxcmd

def check_simmx_args(options, nofilecheck=False):
	if not hasattr(options,"simmxfile"): setattr(options,"simmxfile","dummy")
	number_options_file(0,"proj_simmx",options,"simmxprojfile")
	number_options_file(0,"proj_stg1",options,"projstg1file")
	number_options_file(0,"simmx_stg1",options,"simmxstg1file")
	cmd = get_simmx_cmd(options,True,nofilecheck)
	print ""
	print "#### Test executing simmx command: %s" %cmd
	return (False)
	return ( os.system(cmd) != 0)

def get_projection_cmd(options,check=False):
	
	model = options.model
	if hasattr(options,"filtered_model"): # sometimes there is a filtered model, i.e. if mass or automask3d is specified
		model = options.filtered_model
		
	e2projcmd = "e2project3d.py %s -f --sym=%s --projector=%s --outfile=%s --orientgen=%s --postprocess=normalize.circlemean" %(model,options.sym,options.projector,options.projfile,options.orientgen)
	
	if options.parallel: e2projcmd += " --parallel=%s" %options.parallel
	
	if ( check ):
		e2projcmd += " --check"	
		
	if (options.verbose):
		e2projcmd += " --verbose=" + str(options.verbose - 1)
	
	return e2projcmd
	
def check_projection_args(options):
	if not hasattr(options,"projfile"): setattr(options,"projfile","dummy")
	cmd = get_projection_cmd(options,True)
	print ""
	print "#### Test executing projection command: %s" %cmd
	return ( os.system(cmd) != 0 )

def check(options,verbose=0):
	if (options.verbose>0):
		print ""
		print "#### Testing directory contents and command line arguments for e2refine.py"
	
	error = False
	if options.input == None or not file_exists(options.input):
		print "Error: failed to find input file %s" %options.input
		error = True
	
	if options.usefilt != None:
		if not file_exists(options.usefilt):
			print "Error: failed to find usefilt file %s" %options.usefilt
			error = True
		n1 = EMUtil.get_image_count(options.usefilt)
		n2 = EMUtil.get_image_count(options.input)
		if n1 != n2:
			print "Error, the number of images in the starting particle set:",n2,"does not match the number in the usefilt set:",n1
			error = True
		read_header_only=True
		img1 = EMData()
		img1.read_image(options.input,0,read_header_only)
		img2 = EMData()
		img2.read_image(options.usefilt,0,read_header_only)
		
		nx1 = img1.get_attr("nx") 
		nx2 = img2.get_attr("nx") 
		
		ny1 = img1.get_attr("ny") 
		ny2 = img2.get_attr("ny") 
		
		if nx1 != nx2 or ny1 != ny2:
			error = True
			if options.verbose>0: print "Error, the dimensions of particle data (%i x %i) and the usefilt data (%i x %i) do not match" %(nx1,ny1,nx2,ny2)
	
	if not file_exists(options.model):
		print "Error: 3D image %s does not exist" %options.model
		error = True
	
	if options.mass:
		if options.mass <=0:
			print "If specified, the mass argument must be greater than 0"
			error = True
		if not options.apix:
			print "If you specify the mass argument, you must also specify the apix argument"
			error = True
			
	if options.apix:
		if options.apix <=0:
			print "If specified, the apix argument must be greater than 0"
			error = True
		
	if options.automask3d:
		vals = options.automask3d.split(",")
		mapping = ["threshold","radius","nshells","nshellsgauss","nmaxseed"]
		if len(vals) != 5:
			print "If specified, the automask3d options must provide 5 parameters (threshold,radius,nshells,nshellsgauss,nmaxseed), for example --automask3d=1.7,0,5,5,3"
			error = True
		else:
			# Here I turn options.automask3d into what we would have expected if the user was supplying the whole processor argument,
			# e.g. --automask3d=mask.auto3d:threshold=1.7:radi.... etc. I also add the return_mask=1 parameters - this could be misleading for future
			# programmers, who will potentially wander where it came from
			s = "mask.auto3d"
			for i,p in enumerate(mapping):
				s += ":"+p+"="+vals[i]
			s+= ":return_mask=1"
			options.automask3d = s
			
			if not check_eman2_type(options.automask3d,Processors,"Processors"): error = True
	
	if not options.iter:
		print "Error: you must specify the --it argument"
		error = True
		
	if ( file_exists(options.model) and options.input != None and file_exists(options.input)):
		(xsize, ysize ) = gimme_image_dimensions2D(options.input);
		(xsize3d,ysize3d,zsize3d) = gimme_image_dimensions3D(options.model)
		
		if (options.verbose>0):
			print "%s contains %d images of dimensions %dx%d" %(options.input,EMUtil.get_image_count(options.input),xsize,ysize)
			print "%s has dimensions %dx%dx%d" %(options.model,xsize3d,ysize3d,zsize3d)
		
		if ( xsize != ysize ):
			if ( ysize == zsize3d and xsize == ysize3d and ysize3D == xsize3d ):
				print "Error: it appears as though you are trying to do helical reconstruction. This is not supported"
				error = True
			else:
				print "Error: images dimensions (%d x %d) of %s are not identical. This mode of operation is not supported" %(xsize,ysize,options.input)
				error = True
		
		if ( xsize3d != ysize3d or ysize3d != zsize3d ):
			print "Error: image dimensions (%dx%dx%d) of %s are not equal" %(xsize3d,ysize3d,zsize3d,options.model)
			error = True
			
		if ( xsize3d != xsize ) :
			print "Error: the dimensions of %s (%d) do not match the dimension of %s (%d)" %(options.input,xsize,options.model,xsize3d)
			error = True
			
	if options.path != None:
		if not os.path.exists(options.path):
			print "Error: the path %s does not exist" %options.path
			error = True
	else:
		options.path = numbered_path("refine",True)
			
	if hasattr(options,"parallel") and options.parallel != None:
  		if len(options.parallel) < 2:
  			print "The parallel option %s does not make sense" %options.parallel
  			error = True
 
	if (options.verbose>0):
		if (error):
			s = "FAILED"
		else:
			s = "PASSED"
			
		print "e2refine.py test.... %s" %s

	return error == True
	
if __name__ == "__main__":
    main()
