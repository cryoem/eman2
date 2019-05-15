#!/usr/bin/env python
from __future__ import print_function
from __future__ import division

#
# Author: Steve Ludtke, 07/26/17 (sludtke@bcm.edu)
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


from past.utils import old_div
from builtins import range
from EMAN2 import *
from EMAN2db import db_open_dict, db_list_dicts
from math import *
from os import remove
import time
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

	This variant of the program uses rotational/translational invariants derived from the bispectrum
	of each particle."""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	# we grab all relevant options from e2refine.py for consistency
	# and snag a bunch of related code from David

	#options associated with e2refine2d.py
	parser.add_argument("--path",type=str,default=None,help="Path for the refinement, default=auto")
	parser.add_argument("--input", default=None,type=str, help="The name of the file containing the particle data", browser='EMSetsTable(withmodal=True,multiselect=False)', guitype='filebox', row=0, col=0, rowspan=1, colspan=3, mode="spr")
	parser.add_argument("--ncls", default=32, type=int, help="Number of classes to generate", guitype='intbox', row=1, col=0, rowspan=1, colspan=1, mode="spr")
	parser.add_argument("--alignsort", default=False, action="store_true",help="This will align and sort the final class-averages based on mutual similarity.", guitype='boolbox', row=1, col=1, rowspan=1, colspan=1, mode="spr[True]")
	parser.add_argument("--msamode",default="pca",type=str,help="e2msa can use a variety of different dimensionality reduction algorithms, the default is Principal Component Analysis (PCA), but others are available, see e2msa.py")
#	parser.add_argument("--normproj", default=False, action="store_true",help="Normalizes each projected vector into the MSA subspace. Note that this is different from normalizing the input images since the subspace is not expected to fully span the image", guitype='boolbox', row=1, col=1, rowspan=1, colspan=1, mode="spr[True]")
#	parser.add_argument("--fastseed", action="store_true", default=False,help="Will seed the k-means loop quickly, but may produce less consistent results. Always use this when generating >~100 classes.",guitype='boolbox', row=1, col=2, rowspan=1, colspan=1, mode="spr[True]")
	parser.add_argument("--iter", type=int, default=0, help = "The total number of refinement iterations to perform")  #, guitype='intbox', row=2, col=0, rowspan=1, colspan=1, mode="spr")
	parser.add_argument("--nbasisfp",type=int,default=8,help="Number of MSA basis vectors to use when classifying particles", guitype='intbox', row=2, col=1, rowspan=1, colspan=1, mode="spr")
#	parser.add_argument("--automask",default=False, action="store_true",help="Automasking during class-averaging to help with centering when particle density is high",guitype="boolbox", row=2,col=2,rowspan=1,colspan=1,mode="spr")
#	parser.add_argument("--naliref", default=5, type=int, help="Number of alignment references to when determining particle orientations", guitype='intbox', row=3, col=0, rowspan=1, colspan=1, mode="spr")
	parser.add_argument("--parallel","-P",type=str,help="Run in parallel, specify type:<option>=<value>:<option>:<value>",default="thread:4", guitype='strbox', row=4, col=0, rowspan=1, colspan=3, mode="spr")
	parser.add_argument("--threads", default=4,type=int,help="Number of threads to run in parallel on a single computer when multi-computer parallelism isn't useful", guitype='intbox', row=30, col=2, rowspan=1, colspan=1, mode="refinement[4]")
#	parser.add_argument("--centeracf", default=False, action="store_true",help="This option has been removed in favor of a new centering algorithm")
	parser.add_argument("--center",type=str,default="xform.center",help="If the default centering algorithm (xform.center) doesn't work well, you can specify one of the others here (e2help.py processor center)",guitype='comboparambox', choicelist='dict(re_filter_list(dump_processors_list(),"xform.center").items()+[("nocenter",["Do not center class averages. (similar to what relion does)"])])', row=3, col=1, rowspan=1, colspan=2, mode="spr")
#	parser.add_argument("--check", "-c",default=False, action="store_true",help="Checks the contents of the current directory to verify that e2refine2d.py command will work - checks for the existence of the necessary starting files and checks their dimensions. Performs no work ")
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higher number means higher level of verboseness")
#	parser.add_argument("--maxshift", default=-1, type=int, help="Maximum particle translation in x and y")
#	parser.add_argument("--exclude", type=str,default=None,help="The named file should contain a set of integers, each representing an image from the input file to exclude.")
#	parser.add_argument("--resume", default=False, action="store_true",help="This will cause a check of the files in the current directory, and the refinement will resume after the last completed iteration. It's ok to alter other parameters.")

	parser.add_header(name="advancedheader", help='*** Advanced Users Only (most projects will not need to change) ***', title="Advanced Options:", row=5, col=0, rowspan=1, colspan=3, mode="spr")
	parser.add_header(name="refine2dheader", help='General Options', title="Similarity Matrix:", row=8, col=0, rowspan=1, colspan=3, mode="spr")

	#options associated with generating initial class-averages
#	parser.add_argument("--initial",type=str,default=None,help="File containing starting class-averages. If not specified, will generate starting averages automatically", guitype='strbox', row=6, col=0, rowspan=1, colspan=2, mode="spr")
#	parser.add_argument("--minchange", type=int,default=-1,help="Minimum number of particles that change group before deicding to terminate. Default = -1 (auto)")

	# options associated with e2simmx.py
#	parser.add_argument("--simalign",type=str,help="The name of an 'aligner' to use prior to comparing the images (default=rotate_translate_tree)", default="rotate_translate_tree", guitype='comboparambox', choicelist='re_filter_list(dump_aligners_list(),\'refine|3d\', 1)', row=11, col=0, rowspan=1, colspan=3, mode="spr")
#	parser.add_argument("--simaligncmp",type=str,help="Name of the aligner along with its construction arguments (default=ccc)",default="ccc", guitype='comboparambox', choicelist='re_filter_list(dump_cmps_list(),\'tomo\', True)', row=12, col=0, rowspan=1, colspan=3, mode="spr")
#	parser.add_argument("--simralign",type=str,help="The name and parameters of the second stage aligner which refines the results of the first alignment", default=None, guitype='comboparambox', choicelist='re_filter_list(dump_aligners_list(),\'refine$\')', row=13, col=0, rowspan=1, colspan=3, mode="spr")
#	parser.add_argument("--simraligncmp",type=str,help="The name and parameters of the comparitor used by the second stage aligner. (default=dot).",default="dot", guitype='comboparambox', choicelist='re_filter_list(dump_cmps_list(),\'tomo\', 1)', row=14, col=0, rowspan=1, colspan=3, mode="spr")
#	parser.add_argument("--simcmp",type=str,help="The name of a 'cmp' to be used in comparing the aligned images (default=ccc)", default="ccc", guitype='comboparambox', choicelist='re_filter_list(dump_cmps_list(),\'tomo\', 1)', row=10, col=0, rowspan=1, colspan=3, mode="spr")
#	parser.add_argument("--shrink", dest="shrink", type=int, default=0, help="Optionally shrink the input particles by an integer amount prior to computing similarity scores. For speed purposes. default=0, no shrinking", guitype='shrinkbox', row=9, col=0, rowspan=1, colspan=1, mode="spr")

	parser.add_header(name="caheader", help='Options below this label are specific to class averageing', title="Class averaging:", row=15, col=0, rowspan=1, colspan=3, mode="spr")
	# options associated with e2classaverage
	parser.add_argument("--classkeep",type=float,default=0.8,help="The fraction of particles to keep in each class, based on the similarity score generated by the --cmp argument (default=0.8).", guitype='floatbox', row=16, col=0, rowspan=1, colspan=1, mode="spr")
	parser.add_argument("--classkeepsig", default=False, action="store_true", help="Change the keep (\'--keep\') criterion from fraction-based to sigma-based.", guitype='boolbox', row=16, col=2, rowspan=1, colspan=1, mode="spr")
	parser.add_argument("--classiter", default=4, type=int, help="Number of iterations to use when making class-averages (default=4)", guitype='intbox', row=16, col=1, rowspan=1, colspan=1, mode="spr")
	parser.add_argument("--classalign",type=str,help="If doing more than one iteration, this is the name and parameters of the 'aligner' used to align particles to the previous class average.", default="rotate_translate_tree:flip=1", guitype='comboparambox', choicelist='re_filter_list(dump_aligners_list(),\'refine|3d\', 1)', row=20, col=0, rowspan=1, colspan=3, mode="spr")
	parser.add_argument("--classaligncmp",type=str,help="This is the name and parameters of the comparitor used by the fist stage aligner  Default is dot.",default="ccc", guitype='comboparambox', choicelist='re_filter_list(dump_cmps_list(),\'tomo\', True)', row=21, col=0, rowspan=1, colspan=3, mode="spr")
	parser.add_argument("--classralign",type=str,help="The second stage aligner which refines the results of the first alignment in class averaging. Default is None.", default=None, guitype='comboparambox', choicelist='re_filter_list(dump_aligners_list(),\'refine$\')', row=22, col=0, rowspan=1, colspan=3, mode="spr")
	parser.add_argument("--classraligncmp",type=str,help="The comparitor used by the second stage aligner in class averageing. Default is dot:normalize=1.",default="ccc", guitype='comboparambox', choicelist='re_filter_list(dump_cmps_list(),\'tomo\', True)', row=23, col=0, rowspan=1, colspan=3, mode="spr")
	parser.add_argument("--classaverager",type=str,help="The averager used to generate the class averages. Default is \'mean\'.",default="ctf.weight.autofilt", guitype='combobox', choicelist='dump_averagers_list()', row=18, col=0, rowspan=1, colspan=2, mode="spr")
	parser.add_argument("--classcmp",type=str,help="The name and parameters of the comparitor used to generate similarity scores, when class averaging. Default is ccc'", default="ccc", guitype='comboparambox', choicelist='re_filter_list(dump_cmps_list(),\'tomo\', 1)', row=19, col=0, rowspan=1, colspan=3, mode="spr")
	parser.add_argument("--classnormproc",type=str,default="normalize.edgemean",help="Normalization applied during class averaging", guitype='comboparambox', choicelist='re_filter_list(dump_processors_list(),\'normalize\')', row=17, col=0, rowspan=1, colspan=3, mode="spr")


	#options associated with e2basis.py


	# Database Metadata storage
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)

	global options
	(options, args) = parser.parse_args()
	subverbose=options.verbose-1
	if subverbose<0: subverbose=0

	msamode=options.msamode
	
	if options.parallel :
		parstr="--parallel="+options.parallel
		if options.parallel[:6]=="thread" :
			options.threads=int(options.parallel.split(":")[-1])
			print("--threads set to match --parallel")
	else : parstr=""

	if options.path and ("/" in options.path or "#" in options.path) :
		print("Path specifier should be the name of a subdirectory to use in the current directory. Neither '/' or '#' can be included. ")
		sys.exit(1)

	if options.path == None:
		fls=[int(i[-2:]) for i in os.listdir(".") if i[:5]=="r2db_" and len(i)==7]
		if len(fls)==0 : fls=[0]
		options.path = "r2db_{:02d}".format(max(fls)+1)
		try: os.mkdir(options.path)
		except: pass


	fit=1
	dcts=os.listdir(options.path)

	total_procs = options.iter*4 + 4 # one for every run command
	proc_tally = 0.0

	# store the input arguments forever in the refinement directory
	db = js_open_dict(options.path+"/0_refine2d_parms.json")
	db.update(vars(options))
	db["commandline"]=" ".join(sys.argv)
	db["timestamp"]=str(time.ctime())
	db.close()

	print("Building initial averages")

	n=EMUtil.get_image_count(options.input)

	# make footprint images (rotational/translational invariants)
	fpfile=options.input.split("__")[0]+"__ctf_flip_invar.lst"
	if not os.path.exists(fpfile):
		print("ERROR: no _invar file found. Please run e2ctf_auto.py in your standard EMAN2 project folder. This program will not work with standalone stack files.")
		sys.exit(1)
		#print("WARNING: ",fpfile," not found. Computing bispectra. This will slow processing. ")
		#fpfile=options.path+"/input_bispec.hdf"
		#run("e2proc2dpar.py {} {} --process filter.highpass.gauss:cutoff_pixels=2 --process math.bispectrum.slice:fp={}:size={} --threads {}".format(options.input,fpfile,bispec_invar_parm[1],bispec_invar_parm[0],options.threads))
	else:
		tmp1=EMData(fpfile,0)
		tmp2=EMData(options.input,0)
		if tmp1.has_attr("is_harmonic_fp") : invmode="harmonic"
		else:
			tmp2=tmp2.process("math.bispectrum.slice",{"fp":bispec_invar_parm[1],"size":bispec_invar_parm[0]})
			if tmp1["nx"]!=tmp2["nx"] or tmp1["ny"]!=tmp2["ny"] :
				print("ERROR: images in ",fpfile," have the wrong dimensions. It is likely that you ran e2ctf_auto.py on an older version of EMAN2. Please rerun CTF autoprocessing, with the 'outputonly' option set to generate the correct bispectra.")
				sys.exit(1)
			invmode="bispec"

	logid=E2init(sys.argv,options.ppid)

	# MSA on the footprints
	fpbasis=options.path+"/basis_00.hdf"
	inputproj=options.path+"/basis_proj_00.hdf"
	if n>10000 : step="--step 0,{}".format((n+10000)//20000)
	else: step=""
	#run("e2msa.py %s %s --normalize --nbasis=%0d --scratchfile=%s/msa_scratch.bin %s"%(fpfile,fpbasis,options.nbasisfp,options.path,step))
	run("e2msa.py %s %s %s --nbasis %0d %s --mode %s --nomean"%(fpfile,fpbasis,inputproj,options.nbasisfp,step,msamode))
	proc_tally += 1.0
	if logid : E2progress(logid,old_div(proc_tally,total_procs))

	# reproject the particle footprints into the basis subspace
#	run("e2basis.py project %s %s %s --oneout --mean1 --verbose=%d"%(fpbasis,fpfile,inputproj,subverbose))
	proc_tally += 1.0
	if logid : E2progress(logid,old_div(proc_tally,total_procs))

	# Classification
	run("e2classifykmeans.py %s --original %s --mininclass 2 --ncls %d --clsmx %s/classmx_00.hdf --onein --fastseed --axes 0-%d"%(inputproj,options.input,options.ncls,options.path,options.nbasisfp-1))

	proc_tally += 1.0
	if logid : E2progress(logid,old_div(proc_tally,total_procs))

	# Make class averages
	cls_cmd = "e2classaverage.py --input=%s --classmx=%s/classmx_00.hdf --output=%s/classes_00.hdf --iter=%d --force --storebad --center=%s" %(options.input,options.path,options.path,options.classiter,options.center)
	cls_cmd += get_classaverage_extras(options)
	run (cls_cmd)

	class_postproc(options,0,invmode)

	proc_tally += 1.0
	if logid : E2progress(logid,old_div(proc_tally,total_procs))

	# this is the iterative refinement loop, but the default with bispectra is to not run any iterative refinement
	for it in range(1,options.iter+1) :
		# first we sort and align the class-averages from the last step

		# MSA on class-average bispectra
		run("e2msa.py {path}/classes_fp_{it1:02d}.hdf {path}/basis_{it:02d}.hdf {path}/basis_proj_{it:02d}.hdf --projin {fpfile} --nbasis {nbasis} --mode {mode} --nomean".format(path=options.path,it1=it-1,it=it,nbasis=options.nbasisfp,mode=msamode,fpfile=fpfile))
		#run("e2msa.py %s/classes_fp_%02d.hdf %s/basis_%02d.hdf  --normalize --nbasis=%0d --scratchfile=%s/msa_scratch.bin"%(options.path,it-1,options.path,it,options.nbasisfp,options.path))
		proc_tally += 1.0
		if logid : E2progress(logid,old_div(proc_tally,total_procs))

		# now project original image bispectra into class-average basis space
		#run("e2basis.py project %s/basis_%02d.hdf %s %s/basis_proj_%02d.hdf --oneout --mean1 --verbose=%d"%(options.path,it,fpfile,options.path,it,subverbose))
		#proc_tally += 1.0
		#if logid : E2progress(logid,proc_tally/total_procs)

		# Classification
		run("e2classifykmeans.py %s/basis_proj_%02d.hdf --original=%s --mininclass=2 --ncls=%d --clsmx=%s/classmx_%02d.hdf --onein --fastseed"%(options.path,it,options.input,options.ncls,options.path,it))
		proc_tally += 1.0
		if logid : E2progress(logid,old_div(proc_tally,total_procs))


		# Make class averages
		cls_cmd = "e2classaverage.py --input=%s --classmx=%s/classmx_%02d.hdf --output=%s/classes_%02d.hdf --iter=%d --force --storebad --center=%s" %(options.input,options.path,it,options.path,it,options.classiter,options.center)
		cls_cmd += get_classaverage_extras(options)
		run (cls_cmd)
		proc_tally += 1.0
		if logid : E2progress(logid,old_div(proc_tally,total_procs))

		class_postproc(options,it,invmode)

	print("e2refine2d.py complete")
	E2end(logid)

def class_postproc(options,it,invmode):

	# bispectra of class-averages

	run("e2proc2d.py %s/classes_%02d.hdf %s/classes_%02d.hdf --inplace --calccont --process=filter.highpass.gauss:cutoff_pixels=2 --process=normalize.circlemean:radius=-5"%(options.path,it,options.path,it))

	if invmode=="bispec" :
		run("e2proc2dpar.py {}/classes_{:02d}.hdf {}/classes_fp_{:02d}.hdf --process math.bispectrum.slice:fp={}:size={} --threads {}".format(options.path,it,options.path,it,bispec_invar_parm[1],bispec_invar_parm[0],options.threads))
	else:
		run("e2proc2dpar.py {}/classes_{:02d}.hdf {}/classes_fp_{:02d}.hdf --process math.harmonicpow:fp=1 --threads {}".format(options.path,it,options.path,it,options.threads))
		
	run("e2stacksort.py %s/classes_fp_%02d.hdf %s/classes_fp_%02d.hdf %s/classes_%02d.hdf %s/classes_%02d.hdf --simcmp=ccc --seqalicen"%(options.path,it,options.path,it,options.path,it,options.path,it))
	
	if options.alignsort:
		run("e2stacksort.py %s/classes_%02d.hdf %s/classes_sort_%02d.hdf --bispec"%(options.path,it,options.path,it))


def get_classaverage_extras(options):
	s = " --align=%s --averager=%s  --keep=%f --cmp=%s --aligncmp=%s --normproc=%s" %(options.classalign,options.classaverager,options.classkeep,options.classcmp,options.classaligncmp,options.classnormproc)

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
	if options.verbose>0 : print("*** ",command)
	if options.verbose>1 : tm=time.time()
	error = launch_childprocess(command)

	if error==11 :
		pass
#	#	print "Segfault running %s\nNormal on some platforms, ignoring"%command
	elif error :
		print("Error running:\n%s"%command)
		exit(1)

	if options.verbose>1 : print(time.time()-tm," seconds to complete")

def get_simmx_cmd(options,refs,simmx,check=False,nofilecheck=False):



	return e2simmxcmd

def get_classaverage_cmd(options,check=False,nofilecheck=False):

	e2cacmd = "e2classaverage.py --input=%s --classmx=%s --force --output=%s --keep=.75 --center %s" %(options.startimg,options.classifyfile,options.cafile,options.center)

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
	print('using bootstrap')
	e2cacmd += " --bootstrap"

	return e2cacmd

if __name__ == "__main__":
    main()
