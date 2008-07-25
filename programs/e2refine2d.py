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
	parser.add_option("--simalign",type="string",help="The name of an 'aligner' to use prior to comparing the images", default="rotate_translate")
	parser.add_option("--simaligncmp",type="string",help="Name of the aligner along with its construction arguments",default="dot")
	parser.add_option("--simralign",type="string",help="The name and parameters of the second stage aligner which refines the results of the first alignment", default=None)
	parser.add_option("--simraligncmp",type="string",help="The name and parameters of the comparitor used by the second stage aligner. Default is dot.",default="dot")
	parser.add_option("--simcmp",type="string",help="The name of a 'cmp' to be used in comparing the aligned images", default="dot:normalize=1")

	# options associated with e2basis.py
	parser.add_option("--normproj", default=False, action="store_true",help="Normalizes each projection vector. Note that this is different from normalizing the input images since the subspace is not expected to fully span the image")

	# Parallelism
	parser.add_option("--parallel","-P",type="string",help="Run in parallel, specify type:n=<proc>:option:option",default=None)
	
	## options associated with e2classify.py
	#parser.add_option("--sep", type="int", help="The number of classes a particle can contribute towards (default is 2)", default=2)
	#parser.add_option("--classifyfile", dest="classifyfile", type = "string", default="e2classify.img", help="The file that will store the classification matrix")
	
	## options associated with e2classaverage.py
	#parser.add_option("--classkeep",type="float",help="The fraction of particles to keep in each class, based on the similarity score generated by the --cmp argument.")
	#parser.add_option("--classkeepsig", type="float",default=1.0, help="Change the keep (\'--keep\') criterion from fraction-based to sigma-based.")
	#parser.add_option("--classiter", type="int", help="The number of iterations to perform. Default is 1.", default=3)
	#parser.add_option("--classalign",type="string",help="If doing more than one iteration, this is the name and parameters of the 'aligner' used to align particles to the previous class average.", default="rotate_translate")
	#parser.add_option("--classaligncmp",type="string",help="This is the name and parameters of the comparitor used by the fist stage aligner  Default is dot.",default="phase")
	#parser.add_option("--classralign",type="string",help="The second stage aligner which refines the results of the first alignment in class averaging. Default is None.", default=None)
	#parser.add_option("--classraligncmp",type="string",help="The comparitor used by the second stage aligner in class averageing. Default is dot:normalize=1.",default="dot:normalize=1")
	#parser.add_option("--classaverager",type="string",help="The averager used to generate the class averages. Default is \'image\'.",default="image")
	#parser.add_option("--classcmp",type="string",help="The name and parameters of the comparitor used to generate similarity scores, when class averaging. Default is \'dot:normalize=1\'", default="dot:normalize=1")
		
	global options
	(options, args) = parser.parse_args()
	subverbose=options.verbose-1
	if subverbose<0: subverbose=0
	#0
	#error = False
	#if check(options,True) == True : 
		#error = True
	#if check_projection_args(options) == True : 
		#error = True
	#if check_simmx_args(options,True) == True :
		#error = True
	#if check_classify_args(options,True) == True :
		#error = True
	#options.cafile = "e2classes.1.img"
	#if check_classaverage_args(options,True) == True :
		#error = True
	#if check_make3d_args(options,True) == True:
		#error = True
	
#	if error:
#		print "Error encountered while, bailing"
#		exit(1)
	
	logid=E2init(sys.argv)

	if options.exclude : excludestr="exclude="+options.exclude
	else: excludestr=""

	if options.maxshift<0 : 
		tmp=EMData(options.input,0)
		options.maxshift=tmp.get_xsize()/3	
	
	if options.parallel :
		parstr="--parallel="+options.parallel
	else : parstr=""

	if options.normproj :
		options.normproj="--normproj"
	else : options.normproj=""
	
	fit=1
	if options.resume :
		while os.access("classes.%02d.hdf"%fit,os.F_OK): fit+=1
		options.initial="classes.%02d.hdf"%fit
		fit+=1
		print "starting at iteration ",fit

	# if we aren't given starting class-averages, make some
	if not options.initial and not os.access("classes.init.hdf",os.R_OK):
		print "Building initial averages"
		
		# make footprint images (rotational/translational invariants)
		fpfile=options.input[:options.input.rfind(".")]+".fp.hdf"
		fpfile=fpfile.split("/")[-1]
		if not os.access(fpfile,os.R_OK) :
			run("e2proc2d.py %s %s --fp --verbose=%d %s"%(options.input,fpfile,subverbose,parstr))
		
		# MSA on the footprints
		fpbasis=options.input[:options.input.rfind(".")]+".fp.basis.hdf"
		fpbasis=fpbasis.split("/")[-1]
		if not os.access(fpbasis,os.R_OK) :
			run("e2msa.py %s %s --nbasis=%0d"%(fpfile,fpbasis,options.nbasisfp))
#			run("e2msa.py %s %s --nbasis=%0d --varimax"%(fpfile,fpbasis,options.nbasisfp))
	
		# reproject the particle footprints into the basis subspace
		inputproj=options.input[:options.input.rfind(".")]+".fp.proj.hdf"
		inputproj=inputproj.split("/")[-1]
		if not os.access(inputproj,os.R_OK) :
			run("e2basis.py project %s %s %s --oneout --verbose=%d %s"%(fpbasis,fpfile,inputproj,subverbose,options.normproj))
		
		# classify the subspace vectors
		try: remove("classmx.00.hdf")
		except: pass
		run("e2classifykmeans.py %s --original=%s --ncls=%d --clsmx=classmx.00.hdf --onein %s"%(inputproj,options.input,options.ncls,excludestr))
		
		# make class-averages
		try: re0move("classes.init.hdf")
		except: pass
		run("e2classaverage.py %s classmx.00.hdf classes.init.hdf --iter=3 --align=rotate_translate_flip:maxshift=%d --averager=image -vf --bootstrap --keep=.9 --cmp=optvariance --aligncmp=optvariance"%(options.input,options.maxshift))
	if not options.initial : options.initial="classes.init.hdf"
		
	print "Using references from ",options.initial
	# this is the main refinement loop
	for it in range(fit,options.iter+1) :		
		# first we sort and align the class-averages from the last step
		run("e2stacksort.py %s allrefs.%02d.hdf --simcmp=optvariance:matchfilt=1 --simalign=rotate_translate:maxshift==%d --center --useali --iterative"%(options.initial,it,options.maxshift))
		
		# Compute a classification basis set
		try: remove("basis.%02d.hdf"%it)
		except: pass
		run("e2msa.py allrefs.%02d.hdf basis.%02d.hdf --nbasis=%d"%(it,it,options.nbasisfp))
#		run("e2msa.py allrefs.%02d.hdf basis.%02d.hdf --nbasis=%d --varimax"%(it,it,options.nbasisfp))
		
		# extract the most different references for alignment
#		run("e2stacksort.py %s aliref.%02d.hdf --simcmp=sqeuclidean --reverse --nsort=%d"%(options.initial,it,options.naliref))

		# extract the averages with the most particles
#		run("e2stacksort.py allrefs.%02d.hdf aliref.%02d.hdf --byptcl --nsort=%d"%(it,it,options.naliref))
		run("e2stacksort.py allrefs.%02d.hdf aliref.%02d.hdf --reverse --nsort=%d --simcmp=sqeuclidean"%(it,it,options.naliref))
		
		# We use e2simmx to compute the optimal particle orientations
		try: remove("simmx.%02d.hdf"%it)
		except: pass
		e2simmxcmd = "e2simmx.py aliref.%02d.hdf %s simmx.%02d.hdf -f --saveali --cmp=%s --align=%s --aligncmp=%s --verbose=%d %s %s"  %(it, options.input,it,options.simcmp,options.simalign,options.simaligncmp,subverbose,excludestr,parstr)
		if options.simralign : e2simmxcmd += " --ralign=%s --raligncmp=%s" %(options.simralign,options.simraligncmp)
		run(e2simmxcmd)
		
		# e2basis projectrot here
		inputproj=options.input[:options.input.rfind(".")]+".%02d.proj.hdf"%it
		inputproj=inputproj.split("/")[-1]
		try: remove(inputproj)
		except: pass
		run("e2basis.py projectrot basis.%02d.hdf %s simmx.%02d.hdf %s --oneout --verbose=%d %s"%(it,options.input,it,inputproj,subverbose,options.normproj))
		
		# classify the subspace vectors
		try: remove("classmx.%02d.hdf"%it)
		except: pass
		run("e2classifykmeans.py %s --original=%s --ncls=%d --clsmx=classmx.%02d.hdf --oneinali %s"%(inputproj,options.input,options.ncls,it,excludestr))
		
		# make class-averages
		try: remove("classes.%02d.hdf"%it)
		except: pass
		run("e2classaverage.py %s classmx.%02d.hdf classes.%02d.hdf --iter=%d --align=rotate_translate:maxshift=%d --averager=image -vf  --keep=.9 --cmp=optvariance --aligncmp=optvariance"%(options.input,it,it,options.iterclassav,options.maxshift))
		
		options.initial="classes.%02d.hdf"%it
			
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
