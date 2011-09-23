#!/usr/bin/env python

#
# Author: Steve Ludtke, 1/18/2008 (sludtke@bcm.edu)
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
from EMAN2db import db_open_dict, db_list_dicts
from optparse import OptionParser
from math import *
from os import system,remove
import sys

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options] 
	
	This program is used to produce reference-free class averages from a population of mixed,
	unaligned particle images. These averages can be used to generate initial models or assess
	the structural variability of the data. They are not normally themselves used as part of
	the single particle reconstruction refinement process, which uses the raw particles in a
	reference-based classification approach. However, with a good structure, projections of
	the final 3-D model should be consistent with the results of this reference-free analysis.
	
	This program uses a fully automated iterative alignment/MSA approach. You should normally
	target a minimum of 10-20 particles per class-average, though more is fine.
	
	Default parameters should give a good start, but are likely not optimal for any given system.
	
	Note that it does have the --parallel option, but a few steps of the iterative process
	are not parallelized, so don't be surprised if multiple cores are not always active."""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	# we grab all relevant options from e2refine.py for consistency
	# and snag a bunch of related code from David
	
	#options associated with e2refine2d.py
	parser.add_header(name="refine2dheader", help='Options below this label are specific to e2refine2d', title="### e2refine2d options ###", row=1, col=0, rowspan=1, colspan=3)
	parser.add_argument("--path",type=str,default=None,help="Path for the refinement, default=auto")
	parser.add_argument("--iter", type=int, default=8, help = "The total number of refinement iterations to perform", guitype='intbox', row=3, col=0, rowspan=1, colspan=1)
	parser.add_argument("--automask", default=False, action="store_true",help="This will perform a 2-D automask on class-averages to help with centering. May be useful for negative stain data particularly.")
	parser.add_argument("--check", "-c",default=False, action="store_true",help="Checks the contents of the current directory to verify that e2refine2d.py command will work - checks for the existence of the necessary starting files and checks their dimensions. Performs no work ")
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	parser.add_argument("--input", default="start.hdf",type=str, help="The name of the file containing the particle data", guitype='filebox', row=0, col=0, rowspan=1, colspan=3)
	parser.add_argument("--ncls", default=32, type=int, help="Number of classes to generate", guitype='intbox', row=2, col=0, rowspan=1, colspan=1)
	parser.add_argument("--maxshift", default=-1, type=int, help="Maximum particle translation in x and y")
	parser.add_argument("--naliref", default=5, type=int, help="Number of alignment references to when determining particle orientations", guitype='intbox', row=4, col=0, rowspan=1, colspan=1)
	parser.add_argument("--exclude", type=str,default=None,help="The named file should contain a set of integers, each representing an image from the input file to exclude.")
	parser.add_argument("--resume", default=False, action="store_true",help="This will cause a check of the files in the current directory, and the refinement will resume after the last completed iteration. It's ok to alter other parameters.")
	
	#options associated with generating initial class-averages
	parser.add_argument("--initial",type=str,default=None,help="File containing starting class-averages. If not specified, will generate starting averages automatically", guitype='strbox', row=3, col=1, rowspan=1, colspan=2)
	parser.add_argument("--nbasisfp",type=int,default=5,help="Number of MSA basis vectors to use when classifying particles", guitype='intbox', row=4, col=1, rowspan=1, colspan=1)
	parser.add_argument("--minchange", type=int,default=-1,help="Minimum number of particles that change group before deicding to terminate. Default = -1 (auto)")
	parser.add_argument("--fastseed", action="store_true", default=False,help="Will seed the k-means loop quickly, but may produce less consistent results.")

	parser.add_header(name="simmxheader", help='Options below this label are specific to simmx', title="### simmx options ###", row=6, col=0, rowspan=1, colspan=3)
	# options associated with e2simmx.py
	parser.add_argument("--simalign",type=str,help="The name of an 'aligner' to use prior to comparing the images (default=rotate_translate_flip)", default="rotate_translate_flip", guitype='comboparambox', choicelist='re_filter_list(dump_aligners_list(),\'refine|3d\', 1)', row=9, col=0, rowspan=1, colspan=3)
	parser.add_argument("--simaligncmp",type=str,help="Name of the aligner along with its construction arguments (default=ccc)",default="ccc", guitype='comboparambox', choicelist='re_filter_list(dump_cmps_list(),\'tomo\', True)', row=10, col=0, rowspan=1, colspan=3)
	parser.add_argument("--simralign",type=str,help="The name and parameters of the second stage aligner which refines the results of the first alignment", default=None, guitype='comboparambox', choicelist='re_filter_list(dump_aligners_list(),\'refine$\')', row=11, col=0, rowspan=1, colspan=3)
	parser.add_argument("--simraligncmp",type=str,help="The name and parameters of the comparitor used by the second stage aligner. (default=dot).",default="dot", guitype='comboparambox', choicelist='re_filter_list(dump_cmps_list(),\'tomo\', 1)', row=12, col=0, rowspan=1, colspan=3)
	parser.add_argument("--simcmp",type=str,help="The name of a 'cmp' to be used in comparing the aligned images (default=ccc)", default="ccc", guitype='comboparambox', choicelist='re_filter_list(dump_cmps_list(),\'tomo\', 1)', row=8, col=0, rowspan=1, colspan=3)
	parser.add_argument("--shrink", dest="shrink", type=int, default=0, help="Optionally shrink the input particles by an integer amount prior to computing similarity scores. For speed purposes.", guitype='intbox', row=7, col=0, rowspan=1, colspan=1)
	
	parser.add_header(name="caheader", help='Options below this label are specific to class averageing', title="### Class averaging options ###", row=13, col=0, rowspan=1, colspan=3)
	# options associated with e2classaverage
	parser.add_argument("--classkeep",type=float,default=0.85,help="The fraction of particles to keep in each class, based on the similarity score generated by the --cmp argument (default=0.85).", guitype='floatbox', row=14, col=0, rowspan=1, colspan=1)
	parser.add_argument("--classkeepsig", default=False, action="store_true", help="Change the keep (\'--keep\') criterion from fraction-based to sigma-based.", guitype='boolbox', row=14, col=2, rowspan=1, colspan=1)
	parser.add_argument("--classiter", default=5, type=int, help="Number of iterations to use when making class-averages (default=5)", guitype='intbox', row=14, col=1, rowspan=1, colspan=1)
	parser.add_argument("--classalign",type=str,help="If doing more than one iteration, this is the name and parameters of the 'aligner' used to align particles to the previous class average.", default="rotate_translate_flip", guitype='comboparambox', choicelist='re_filter_list(dump_aligners_list(),\'refine|3d\', 1)', row=18, col=0, rowspan=1, colspan=3)
	parser.add_argument("--classaligncmp",type=str,help="This is the name and parameters of the comparitor used by the fist stage aligner  Default is dot.",default="ccc", guitype='comboparambox', choicelist='re_filter_list(dump_cmps_list(),\'tomo\', True)', row=19, col=0, rowspan=1, colspan=3)
	parser.add_argument("--classralign",type=str,help="The second stage aligner which refines the results of the first alignment in class averaging. Default is None.", default=None, guitype='comboparambox', choicelist='re_filter_list(dump_aligners_list(),\'refine$\')', row=20, col=0, rowspan=1, colspan=3)
	parser.add_argument("--classraligncmp",type=str,help="The comparitor used by the second stage aligner in class averageing. Default is dot:normalize=1.",default="ccc", guitype='comboparambox', choicelist='re_filter_list(dump_cmps_list(),\'tomo\', True)', row=21, col=0, rowspan=1, colspan=3)
	parser.add_argument("--classaverager",type=str,help="The averager used to generate the class averages. Default is \'mean\'.",default="mean", guitype='combobox', choicelist='dump_averagers_list()', row=16, col=0, rowspan=1, colspan=2)
	parser.add_argument("--classcmp",type=str,help="The name and parameters of the comparitor used to generate similarity scores, when class averaging. Default is ccc'", default="ccc", guitype='comboparambox', choicelist='re_filter_list(dump_cmps_list(),\'tomo\', 1)', row=17, col=0, rowspan=1, colspan=3)
	parser.add_argument("--classnormproc",type=str,default="normalize.edgemean",help="Normalization applied during class averaging", guitype='comboparambox', choicelist='re_filter_list(dump_processors_list(),\'normalize\')', row=15, col=0, rowspan=1, colspan=3)
	parser.add_argument("--classrefsf",default=False, action="store_true", help="Use the setsfref option in class averaging to produce better filtered averages.", guitype='boolbox', row=16, col=2, rowspan=1, colspan=1)
	
	 
	#options associated with e2basis.py
	parser.add_argument("--normproj", default=False, action="store_true",help="Normalizes each projected vector into the MSA subspace. Note that this is different from normalizing the input images since the subspace is not expected to fully span the image", guitype='boolbox', row=2, col=1, rowspan=1, colspan=1)

	# Parallelism
	parser.add_argument("--parallel","-P",type=str,help="Run in parallel, specify type:<option>=<value>:<option>:<value>",default=None, guitype='strbox', row=5, col=0, rowspan=1, colspan=3)
	
	# Database Metadata storage
	parser.add_argument("--dbls",type=str,help="data base list storage, used by the workflow. You can ignore this argument.",default=None)
	
		
	global options
	(options, args) = parser.parse_args()
	subverbose=options.verbose-1
	if subverbose<0: subverbose=0

	if options.exclude : excludestr="exclude="+options.exclude
	else: excludestr=""

	if options.maxshift<0 : 
		tmp=EMData()
		tmp.read_image(options.input,0)
		options.maxshift=tmp.get_xsize()/3	
	
	if options.parallel :
		parstr="--parallel="+options.parallel
	else : parstr=""

	if options.normproj :
		options.normproj="--normproj"
	else : options.normproj=""
	
	if options.classrefsf : 
		print "Warning: classrefsf option has no effect on e2refine2d.py"
	
	if options.path and ("/" in options.path or "#" in options.path) :
		print "Path specifier should be the name of a subdirectory to use in the current directory. Neither '/' or '#' can be included. "
		sys.exit(1)
	if options.path and options.path[:4].lower()!="bdb:" : options.path="bdb:"+options.path
	if not options.path : options.path="bdb:"+numbered_path("r2d",not options.resume)
	
	logid=E2init(sys.argv)

	fit=1
	dcts=db_list_dicts(options.path)
	if options.resume :
		while "classes_%02d"%fit in dcts: fit+=1

		options.initial=options.path+"#classes_%02d"%fit
		fit+=1
		print "starting at iteration ",fit
	
	total_procs = options.iter*7 + 5 # one for every run command
	proc_tally = 0.0
	
	if options.dbls:
		most_recent_classes = None

	if options.fastseed : fastseed="--fastseed"
	else : fastseed=""

	# if we aren't given starting class-averages, make some
#	if not options.initial and not "classes_init" in dcts:
	if not options.initial:
		print "Building initial averages"
		
		# make footprint images (rotational/translational invariants)
		fpfile=options.path+"#input_fp"
		#run("e2proc2d.py %s %s --fp --verbose=%d --inplace %s"%(options.input,fpfile,subverbose,parstr)) # parallel doesn't work in e2proc2d.py as of November 28 2008 - d.woolford
		run("e2proc2d.py %s %s --fp=0 --verbose=%d --inplace"%(options.input,fpfile,subverbose))
		proc_tally += 1.0
		if logid : E2progress(logid,proc_tally/total_procs)
		
		# MSA on the footprints
		fpbasis=options.path+"#input_fp_basis"
		fpbasis=fpbasis.split("/")[-1]
		run("e2msa.py %s %s --normalize --nbasis=%0d"%(fpfile,fpbasis,options.nbasisfp))
		proc_tally += 1.0
		if logid : E2progress(logid,proc_tally/total_procs)
#			run("e2msa.py %s %s --nbasis=%0d --varimax"%(fpfile,fpbasis,options.nbasisfp))
	
		# reproject the particle footprints into the basis subspace
		inputproj=options.path+"#input_fp_basis_proj"
		run("e2basis.py project %s %s %s --oneout --mean1 --verbose=%d %s"%(fpbasis,fpfile,inputproj,subverbose,options.normproj))
		proc_tally += 1.0
		if logid : E2progress(logid,proc_tally/total_procs)
		# classify the subspace vectors
#		try: db_remove_dict(path+"#classmx_00")
#		except: pass
		run("e2classifykmeans.py %s --original=%s --mininclass=2 --ncls=%d --clsmx=%s#classmx_00 --minchange=%d --onein %s %s"%(inputproj,options.input,options.ncls,options.path,options.minchange,excludestr,fastseed))
		
		proc_tally += 1.0
		if logid : E2progress(logid,proc_tally/total_procs)
		# make class-averages
		#run("e2classaverage.py %s %s#classmx_%02d %s#classes_%02d --iter=%d --align=%s:maxshift=%d --averager=%s -vf  --keep=%f --cmp=%s --aligncmp=%s"%(options.input,options.path,it,options.path,it,options.classiter,options.classalign,options.maxshift,options.classaverager,options.classkeep,options.classcmp,options.classaligncmp))
		
		cls_cmd = "e2classaverage.py --input=%s --classmx=%s#classmx_00 --output=%s#classes_init --iter=8 --force --bootstrap" %(options.input,options.path,options.path)
		cls_cmd += get_classaverage_extras(options)
		
		#run("e2classaverage.py %s %s#classmx_00 %s#classes_init --iter=6 --align=%s:maxshift=%d --averager=%s -vf --bootstrap --keep=%f --cmp=%s --aligncmp=%s --normproc=%s"%(options.input,options.path,options.path,options.classalign,options.maxshift,options.classaverager,options.classkeep,options.classcmp,options.classaligncmp,options.classnormproc))
		run (cls_cmd)
		proc_tally += 1.0
		if logid : E2progress(logid,proc_tally/total_procs)
		
		if options.dbls:
			pdb = db_open_dict("bdb:project")
			tmp_data = pdb.get(options.dbls, dfl={})
			if isinstance(tmp_data,list):
				d = {}
				for name in tmp_data: # this is for back compatibility it could be removed in 2010
					s = {}
					s["Original Data"] = name
					d[name] = s
				tmp_data = d	
			most_recent_classes = "%s#classes_init" %options.path
			s = {}
			s["Original Data"] = most_recent_classes
			tmp_data[most_recent_classes]= s
			# global.spr_ref_free_class_aves
			pdb[options.dbls] = tmp_data
			
	if not options.initial : options.initial=options.path+"#classes_init"
		
	print "Using references from ",options.initial
	# this is the main refinement loop
	for it in range(fit,options.iter+1) :		
		# first we sort and align the class-averages from the last step
		run("e2proc2d.py %s %s#allrefs_%02d --inplace --process=filter.highpass.gauss:cutoff_abs=.02 --process=normalize.edgemean --process=xform.centerofmass:threshold=1"%(options.initial,options.path,it))
		# now we try for mutual alignment of particle orientations
		run("e2stacksort.py %s#allrefs_%02d %s#allrefs_%02d --simcmp=sqeuclidean:normto=1 --simalign=rotate_translate_flip --useali --iterative"%(options.path,it,options.path,it))
		# however we don't want things off-center, so we do a final recentering
		run("e2proc2d.py %s#allrefs_%02d %s#allrefs_%02d --inplace --process=xform.centerofmass:threshold=1"%(options.path,it,options.path,it))
		proc_tally += 1.0
		if logid : E2progress(logid,proc_tally/total_procs)
		# Compute a classification basis set
		run("e2msa.py %s#allrefs_%02d %s#basis_%02d --normalize --nbasis=%d "%(options.path,it,options.path,it,options.nbasisfp))
		proc_tally += 1.0
		if logid : E2progress(logid,proc_tally/total_procs)
#		run("e2msa.py allrefs.%02d.hdf basis.%02d.hdf --nbasis=%d --varimax"%(it,it,options.nbasisfp))
		
		# extract the most different references for alignment
#		run("e2stacksort.py %s aliref.%02d.hdf --simcmp=sqeuclidean --reverse --nsort=%d"%(options.initial,it,options.naliref))

		# sort by particles
		run("e2stacksort.py %s#allrefs_%02d %s#allrefs_%02d --bykurtosis"%(options.path,it,options.path,it))

       	# now extract most different refs. ninput eliminates 25% particles with lowest kurtosis from consideration
		run("e2stacksort.py %s#allrefs_%02d %s#aliref_%02d --reverse --ninput=%d --nsort=%d --simcmp=ccc --simalign=rotate_translate_flip"%(options.path,it,options.path,it,options.ncls*3/4,options.naliref))
		proc_tally += 1.0
		if logid : E2progress(logid,proc_tally/total_procs)
		# We use e2simmx to compute the optimal particle orientations
		# ERROR e2simmxy doesn't support parallel
		#e2simmxcmd = "e2simmx.py %s#aliref_%02d %s %s#simmx_%02d -f --saveali --cmp=%s --align=%s --aligncmp=%s --verbose=%d %s %s"%(options.path,it, options.input,options.path,it,options.simcmp,options.simalign,options.simaligncmp,subverbose,excludestr,parstr) # e2simmx doesn't do parallel
		e2simmxcmd = "e2simmx.py %s#aliref_%02d %s %s#simmx_%02d -f --saveali --cmp=%s --align=%s --aligncmp=%s --verbose=%d %s"%(options.path,it, options.input,options.path,it,options.simcmp,options.simalign,options.simaligncmp,subverbose,excludestr)
		if options.simralign : e2simmxcmd += " --ralign=%s --raligncmp=%s" %(options.simralign,options.simraligncmp)
		if options.parallel: e2simmxcmd += " --parallel=%s" %options.parallel
		if options.shrink: e2simmxcmd += " --shrink=%d" %options.shrink
		run(e2simmxcmd)
		proc_tally += 1.0
		if logid : E2progress(logid,proc_tally/total_procs)
		
		# e2basis projectrot here
		inputproj=options.path+"#input_%02d_proj"%it
		run("e2basis.py projectrot %s#basis_%02d %s %s#simmx_%02d %s --oneout --mean1 --normproj --verbose=%d %s"%(options.path,it,options.input,options.path,it,inputproj,subverbose,options.normproj))
		
		# classify the subspace vectors
		run("e2classifykmeans.py %s --original=%s --mininclass=2 --ncls=%d --clsmx=%s#classmx_%02d --minchange=%d --oneinali %s %s"%(inputproj,options.input,options.ncls,options.path,it,options.minchange,excludestr,fastseed))
		proc_tally += 1.0
		if logid : E2progress(logid,proc_tally/total_procs)
		
		# make class-averages
		cls_cmd = "e2classaverage.py --input=%s --classmx=%s#classmx_%02d --output=%s#classes_%02d --force --iter=%d " %(options.input,options.path,it,options.path,it,options.classiter)
		cls_cmd += get_classaverage_extras(options)
		#run("e2classaverage.py %s %s#classmx_%02d %s#classes_%02d --iter=%d --align=%s:maxshift=%d --averager=%s -vf  --keep=%f --cmp=%s --aligncmp=%s"%(options.input,options.path,it,options.path,it,options.classiter,options.classalign,options.maxshift,options.classaverager,options.classkeep,options.classcmp,options.classaligncmp))
		run(cls_cmd)
		
		if options.dbls:
			pdb = db_open_dict("bdb:project")
			tmp_data = pdb.get(options.dbls, dfl={})
			if most_recent_classes != None:
				try:
					tmp_data.pop(most_recent_classes)
				except: pass # the user removed it in the workflow!
			most_recent_classes = "%s#classes_%02d" %(options.path,it)
			s = {}
			s["Original Data"] = most_recent_classes
			tmp_data[most_recent_classes] = s
			# global.spr_ref_free_class_aves
			pdb[options.dbls] = tmp_data
		
		proc_tally += 1.0
		if logid : E2progress(logid,proc_tally/total_procs)
		
		options.initial=options.path+"#classes_%02d"%it
			
	E2end(logid)
	

def get_classaverage_extras(options):
	s = " --align=%s:maxshift=%d --averager=%s  --keep=%f --cmp=%s --aligncmp=%s --normproc=%s" %(options.classalign,options.maxshift,options.classaverager,options.classkeep,options.classcmp,options.classaligncmp,options.classnormproc)

	if options.automask : s+= " --automask"
	if options.classkeepsig:
		s += " --keepsig"
	if options.classralign != None:
		s += " --ralign=%s --raligncmp=%s" %(options.classralign,options.classraligncmp)
	if options.parallel != None:
		s += " --parallel=%s" %options.parallel
	
	return s

def run(command):
	"Execute a command with optional verbose output"
	global options
	if options.verbose>0 : print "***************",command
	error = system(command)
	if error==11 :
		pass
#		print "Segfault running %s\nNormal on some platforms, ignoring"%command
	elif error : 
		print "Error running:\n%s"%command
		exit(1)

def get_simmx_cmd(options,refs,simmx,check=False,nofilecheck=False):
	
		
	
	return e2simmxcmd

def get_classaverage_cmd(options,check=False,nofilecheck=False):
	
	e2cacmd = "e2classaverage.py --input=%s --classmx=%s --force --output=%s --keep=.75" %(options.startimg,options.classifyfile,options.cafile)
	
	e2cacmd += " --ref=%s --iter=%d -f" %(options.projfile,options.classiter)
	
	#if (options.classkeepsig):
		#e2cacmd += " --keepsig=%f" %options.classkeepsig
	#elif (options.classkeep):
		#e2cacmd += " --keep=%f" %options.classkeep
	
	if (options.classiter > 1 ):
		e2cacmd += " --cmp=%s --align=%s --aligncmp=%s" %(options.classcmp,options.classalign,options.classaligncmp)

		if (options.classralign != None):
			e2cacmd += " --ralign=%s --raligncmp=%s" %(options.classralign,options.classraligncmp)
	
	if (options.verbose>0):
		e2cacmd += " -verbose=%d"%(options.verbose-1)
	
	if ( check ):
		e2cacmd += " --check"	
			
	if ( nofilecheck ):
		e2cacmd += " --nofilecheck"
	
	# We need to tell e2classaverage.py to bootstrap the original class average, because there are is no alignment
	print 'using bootstrap'
	e2cacmd += " --bootstrap"	
	
	return e2cacmd

def check_e2refin2d_args(options): # this function is required by the workflow, it is a little specialized, but it makes sense to be here
	'''
	Returns a list of error messages based on the input options
	List is empty if there are no messages
	Currently doesn't check the input file, because it's used from the workflow, not from e2refine2d iteself
	'''
	error_message = []
#	if len(options.filenames) == 0: # this is the specialized part - the workflow creates starting data sets from a list of filenames
#		error_message.append("Please choose the file(s) that you want to use as as input data for e2refine2d")
# 		 
 	if options.shrink < 1:
 		error_message.append("Shrink must be atleast 1")
 		 
  	if options.initial != None and len(options.initial) > 0:
 		if not file_exists(options.initial):
 			error_message.append("The initial class averages file you specified (%s) does not exist." %(options.initial))
 		
 	if options.iter < 0:
 		error_message.append("The number of e2refine2d iterations must be atleast 0.")
 		
 	if options.classiter < 0:
 		error_message.append("The number of class average iterations iteration must be atleast 0.")

  	if options.naliref < 1:
  		error_message.append("The number alignment references must be atleast 1.")
  		
  	if options.nbasisfp < 1:
  		error_message.append("The number of MSA basis vectors must be atleast 1.")
  	
  	if options.ncls < 2:
  		error_message.append("The number of classes must be atleast 2.")
  	
  	if options.parallel < 1:
  		error_message.append("The number CPUs availables must be atleast 1.")
  		
  	if  not check_eman2_type(options.simalign,Aligners,"Aligner",False):
  		error_message.append("There is problem with the aligner arguments.")
  		
  	if not check_eman2_type(options.simaligncmp,Cmps,"Cmp",False):
  		error_message.append("There is problem with main aligner comparitor arguments.")
  	
  	if not check_eman2_type(options.simcmp,Cmps,"Cmp",False):
  		error_message.append("There is problem with main comparitor arguments.")
  	
  	if options.simralign != None and not check_eman2_type(options.simralign,Aligners,"Aligner"):
  		error_message.append("There is problem with the refine aligner arguments.")
  		
  	if options.simraligncmp != None and  not check_eman2_type(options.simraligncmp,Cmps,"Cmps"):
  		error_message.append("There is problem with the refine aligner comparitor arguments.")

  	if hasattr(options,"parallel") and options.parallel != None:
  		if len(options.parallel) < 2:
  			error_message.append("The parallel option %s does not make sense" %options.parallel)
  		
  	return error_message

if __name__ == "__main__":
    main()
