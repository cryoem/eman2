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
from optparse import OptionParser
from math import *
from os import system,remove
import sys

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """%prog [options] 
	THIS PROGRAM IS NOT YET COMPLETE

	EMAN2 iterative single particle reconstruction"""
	parser = OptionParser(usage=usage,version=EMANVERSION)

	# we grab all relevant options from e2refine.py for consistency
	# and snag a bunch of related code from David
	
	#options associated with e2refine2d.py
	parser.add_option("--path",type="string",default=None,help="Path for the refinement, default=auto")
	parser.add_option("--iter", type = "int", default=0, help = "The total number of refinement iterations to perform")
	parser.add_option("--check", "-c",default=False, action="store_true",help="Checks the contents of the current directory to verify that e2refine2d.py command will work - checks for the existence of the necessary starting files and checks their dimensions. Performs no work ")
	parser.add_option("--verbose","-v", type="int", default=0,help="Verbosity of output (1-9)")
	parser.add_option("--input", default="start.hdf",type="string", help="The name of the file containing the particle data")
	parser.add_option("--ncls", default=32, type="int", help="Number of classes to generate")
	parser.add_option("--maxshift", default=-1, type="int", help="Maximum particle translation in x and y")
	parser.add_option("--iterclassav", default=2, type="int", help="Number of iterations when making class-averages")
	parser.add_option("--naliref", default=8, type="int", help="Number of alignment references to when determining particle orientations")
	parser.add_option("--exclude", type="string",default=None,help="The named file should contain a set of integers, each representing an image from the input file to exclude.")
	parser.add_option("--resume", default=False, action="store_true",help="This will cause a check of the files in the current directory, and the refinement will resume after the last completed iteration. It's ok to alter other parameters.")
	
	#options associated with generating initial class-averages
	parser.add_option("--initial",type="string",default=None,help="File containing starting class-averages. If not specified, will generate starting averages automatically")
	parser.add_option("--nbasisfp",type="int",default=5,help="Number of MSA basis vectors to use when classifiying based on invariants for making starting class-averages")

	# options associated with e2simmx.py
	parser.add_option("--simalign",type="string",help="The name of an 'aligner' to use prior to comparing the images (default=rotate_translate_flip)", default="rotate_translate_flip")
	parser.add_option("--simaligncmp",type="string",help="Name of the aligner along with its construction arguments (default=dot)",default="frc")
	parser.add_option("--simralign",type="string",help="The name and parameters of the second stage aligner which refines the results of the first alignment", default=None)
	parser.add_option("--simraligncmp",type="string",help="The name and parameters of the comparitor used by the second stage aligner. (default=dot).",default="dot")
	parser.add_option("--simcmp",type="string",help="The name of a 'cmp' to be used in comparing the aligned images (default=dot:normalize=1)", default="dot:normalize=1")

	# options associated with e2basis.py
	parser.add_option("--normproj", default=False, action="store_true",help="Normalizes each projected vector into the MSA subspace. Note that this is different from normalizing the input images since the subspace is not expected to fully span the image")

	# Parallelism
	parser.add_option("--parallel","-P",type="string",help="Run in parallel, specify type:n=<proc>:option:option",default=None)
	
		
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

	# if we aren't given starting class-averages, make some
	if not options.initial and not "classes_init" in dcts:
		print "Building initial averages"
		
		# make footprint images (rotational/translational invariants)
		fpfile=options.path+"#input_fp"
		run("e2proc2d.py %s %s --fp --verbose=%d --inplace %s"%(options.input,fpfile,subverbose,parstr))
		
		# MSA on the footprints
		fpbasis=options.path+"#input_fp_basis"
		fpbasis=fpbasis.split("/")[-1]
		run("e2msa.py %s %s --nbasis=%0d"%(fpfile,fpbasis,options.nbasisfp))
#			run("e2msa.py %s %s --nbasis=%0d --varimax"%(fpfile,fpbasis,options.nbasisfp))
	
		# reproject the particle footprints into the basis subspace
		inputproj=options.path+"#input_fp_basis_proj"
		run("e2basis.py project %s %s %s --oneout --verbose=%d %s"%(fpbasis,fpfile,inputproj,subverbose,options.normproj))
		
		# classify the subspace vectors
#		try: db_remove_dict(path+"#classmx_00")
#		except: pass
		run("e2classifykmeans.py %s --original=%s --ncls=%d --clsmx=%s#classmx_00 --onein %s"%(inputproj,options.input,options.ncls,options.path,excludestr))
		
		# make class-averages
		run("e2classaverage.py %s %s#classmx_00 %s#classes_init --iter=6 --align=rotate_translate_flip:maxshift=%d --averager=image -vf --bootstrap --keep=.9 --cmp=frc --aligncmp=frc"%(options.input,options.path,options.path,options.maxshift))
	if not options.initial : options.initial=options.path+"#classes_init"
		
	print "Using references from ",options.initial
	# this is the main refinement loop
	for it in range(fit,options.iter+1) :		
		# first we sort and align the class-averages from the last step
		run("e2stacksort.py %s %s#allrefs_%02d --simcmp=optvariance:matchfilt=1 --simalign=rotate_translate:maxshift==%d --center --useali --iterative"%
		    (options.initial,options.path,it,options.maxshift))
		
		# Compute a classification basis set
		run("e2msa.py %s#allrefs_%02d %s#basis_%02d --nbasis=%d"%(options.path,it,options.path,it,options.nbasisfp))
#		run("e2msa.py allrefs.%02d.hdf basis.%02d.hdf --nbasis=%d --varimax"%(it,it,options.nbasisfp))
		
		# extract the most different references for alignment
#		run("e2stacksort.py %s aliref.%02d.hdf --simcmp=sqeuclidean --reverse --nsort=%d"%(options.initial,it,options.naliref))

		# extract the averages with the most particles
#		run("e2stacksort.py allrefs.%02d.hdf aliref.%02d.hdf --byptcl --nsort=%d"%(it,it,options.naliref))
		run("e2stacksort.py %s#allrefs_%02d %s#aliref_%02d --reverse --nsort=%d --simcmp=sqeuclidean"%(options.path,it,options.path,it,options.naliref))
		
		# We use e2simmx to compute the optimal particle orientations
		e2simmxcmd = "e2simmx.py %s#aliref_%02d %s %s#simmx_%02d -f --saveali --cmp=%s --align=%s --aligncmp=%s --verbose=%d %s %s"%(options.path,it, options.input,options.path,it,options.simcmp,options.simalign,options.simaligncmp,subverbose,excludestr,parstr)
		if options.simralign : e2simmxcmd += " --ralign=%s --raligncmp=%s" %(options.simralign,options.simraligncmp)
		run(e2simmxcmd)
		
		# e2basis projectrot here
		inputproj=options.path+"#input_%02d_proj"%it
		run("e2basis.py projectrot %s#basis_%02d %s %s#simmx_%02d %s --oneout --verbose=%d %s"%(options.path,it,options.input,options.path,it,inputproj,subverbose,options.normproj))
		
		# classify the subspace vectors
		run("e2classifykmeans.py %s --original=%s --ncls=%d --clsmx=%s#classmx_%02d --oneinali %s"%(inputproj,options.input,options.ncls,options.path,it,excludestr))
		
		# make class-averages
		run("e2classaverage.py %s %s#classmx_%02d %s#classes_%02d --iter=%d --align=rotate_translate_flip:maxshift=%d --averager=image -vf  --keep=.9 --cmp=frc --aligncmp=frc"%(options.input,options.path,it,options.path,it,options.iterclassav,options.maxshift))
		
		options.initial=options.path+"#classes_%02d"%it
			
	E2end(logid)
	
	
def run(command):
	"Execute a command with optional verbose output"
	global options
	if options.verbose : print "***************",command
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
	
	e2cacmd = "e2classaverage.py %s %s %s" %(options.startimg,options.classifyfile,options.cafile)
	
	e2cacmd += " --ref=%s --iter=%d -f" %(options.projfile,options.classiter)
	
	if (options.classkeepsig):
		e2cacmd += " --keepsig=%f" %options.classkeepsig
	elif (options.classkeep):
		e2cacmd += " --keep=%f" %options.classkeep
	
	if (options.classiter > 1 ):
		e2cacmd += " --cmp=%s --align=%s --aligncmp=%s" %(options.classcmp,options.classalign,options.classaligncmp)

		if (options.classralign != None):
			e2cacmd += " --ralign=%s --raligncmp=%s" %(options.classralign,options.classraligncmp)
	
	if (options.verbose):
		e2cacmd += " -v"
	
	if ( check ):
		e2cacmd += " --check"	
			
	if ( nofilecheck ):
		e2cacmd += " --nofilecheck"
	
	# We need to tell e2classaverage.py to bootstrap the original class average, because there are is no alignment
	print 'using bootstrap'
	e2cacmd += " --bootstrap"	
	
	return e2cacmd

	
if __name__ == "__main__":
    main()
