#!/usr/bin/env python

#
# Author: Steven Ludtke (sludtke@bcm.edu) 02/06/2012
# Copyright (c) 2000-2012 Baylor College of Medicine
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
import traceback
import e2refine

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options] 
	
	WARNING: THIS PROGRAM IS OBSOLETE. E2REFINE_EASY SHOULD ALWAYS BE USED INSTEAD.

	This program is used for robust resolution evaluation. It avoids noise and initial 
	model bias by using an initial model with phases randomized beyond some specified 
	resolution, splitting the data into even/odd halves, then refining to near convergence
	before comparing the maps. This is a self-consistency check in addition to a resolution test,
	so in the presence of structural variability in the data, you will determine the resolution at
	which the two resultant structures are self-consistent. In the presence of, say, discrete heterogentity
	results may vary from one run to the next.
	
	Command-line parameters are identical to e2refine.py, with a couple of optional additions. Normally you would
	run this command with identical options to the corresponding e2refine.py. The 'path' and other parameters will
	be internally modified for purposes of the even/odd split.
"""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
		
	#options associated with e2refine.py
	parser.add_header(name="refineheader", help='Options below this label are specific to e2refine', title="### e2refine options ###", row=1, col=0, rowspan=1, colspan=3, mode="refinement")
	parser.add_header(name="modelheader", help='Options below this label are specific to e2refine Model', title="### e2refine model options ###", row=4, col=0, rowspan=1, colspan=3, mode="refinement")
	parser.add_argument("--iter", dest = "iter", type = int, default=0, help = "The total number of refinement iterations to perform", guitype='intbox', row=2, col=2, rowspan=1, colspan=1, mode="refinement")
	parser.add_argument("--startiter", dest = "startiter", type = int, default=0, help = "If a refinement crashes, this can be used to pick it up where it left off. This should NOT be used to change parameters, but only to resume an incomplete run.")
	parser.add_argument("--check", "-c", dest="check", default=False, action="store_true",help="Checks the contents of the current directory to verify that e2refine.py command will work - checks for the existence of the necessary starting files and checks their dimensions. Performs no work ")
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	parser.add_argument("--input", dest="input", default=None,type=str, help="The name of the set containing the particle data", browser='EMSetsTable(withmodal=True,multiselect=False)', guitype='filebox', row=0, col=0, rowspan=1, colspan=3, mode="refinement")
	parser.add_argument("--model", dest="model", type=str,default="threed.0a.mrc", help="The name of the 3D image that will seed the refinement", guitype='filebox', browser='EMModelsTable(withmodal=True,multiselect=False)', row=5, col=0, rowspan=1, colspan=3, mode="refinement")
	parser.add_argument("--randomres",type=float,default=25.0,help="Resolution for the lowpass phase-randomization to apply to the initial model. Specify in A. (default=25)")
	parser.add_argument("--usefilt", dest="usefilt", type=str,default=None, help="Specify a particle data file that has been low pass or Wiener filtered. Has a one to one correspondence with your particle data. If specified will be used in projection matching routines, and elsewhere.")
	parser.add_argument("--path", default=None, type=str,help="The name of a directory where results are placed. If unspecified will generate one automatically of type refine_??.")
	parser.add_argument("--mass", default=0, type=float,help="The mass of the particle in kilodaltons, used to run normalize.bymass. If unspecified (set to 0) nothing happens. Requires the --apix argument.", guitype='floatbox', row=2, col=1, rowspan=1, colspan=1, mode="refinement['self.pm().getMass()']")
	parser.add_argument("--apix", default=0, type=float,help="The angstrom per pixel of the input particles. This argument is required if you specify the --mass argument. \n If unspecified (set to 0), the convergence plot is generated using either the project apix, or if not an apix of 1.", guitype='floatbox', row=2, col=0, rowspan=1, colspan=1, mode="refinement['self.pm().getAPIX()']")
	parser.add_argument("--automask3d", default=None, type=str,help="The 5 parameters of the mask.auto3d processor, applied after 3D reconstruction. \n These parameters are, in order, isosurface threshold, radius, nshells, and ngaussshells. \n From e2proc3d.py you could achieve the same thing using: \n --process=mask.auto3d:threshold=1.1:radius=30:nshells=5:ngaussshells=5.", guitype='automask3d', row=6, col=0, rowspan=1, colspan=3, mode="refinement")
	parser.add_argument("--automaskalign",action="store_true",help="This will use the automask to improve 2-D alignments and classification.",default=False, guitype='boolbox', row=7, col=0, rowspan=1, colspan=1,mode="refinement")

	# options associated with e2project3d.py
	parser.add_header(name="projectheader", help='Options below this label are specific to e2project', title="### e2project options ###", row=7, col=0, rowspan=1, colspan=3, mode="refinement")
	parser.add_argument("--sym", dest = "sym", default="c1", help = "Specify symmetry - choices are: c<n>, d<n>, h<n>, tet, oct, icos. For asymmetric reconstruction omit this option or specify c1.", guitype='symbox', row=10, col=0, rowspan=1, colspan=3, mode="refinement")
	parser.add_argument("--projector", dest = "projector", default = "standard",help = "Projector to use", guitype='comboparambox', choicelist='dump_projectors_list()', row=8, col=0, rowspan=1, colspan=3, mode="refinement")
	parser.add_argument("--orientgen", type = str, default='eman:delta=5.0:inc_mirror=0:perturb=1',help = "The orientation generation argument for e2project3d.py", guitype='comboparambox', choicelist='dump_orientgens_list()', row=9, col=0, rowspan=1, colspan=3, mode="refinement")
		
	# options associated with e2simmx.py
	parser.add_header(name="simmxheader", help='Options below this label are specific to e2simmx', title="### e2simmx options ###", row=11, col=0, rowspan=1, colspan=3, mode="refinement")
	parser.add_argument("--simalign",type=str,help="The name of an 'aligner' to use prior to comparing the images", default="rotate_translate_flip", guitype='comboparambox', choicelist='re_filter_list(dump_aligners_list(),\'refine|3d\', 1)', row=14, col=0, rowspan=1, colspan=3, mode="refinement")
	parser.add_argument("--simaligncmp",type=str,help="Name and parameters of the comparator used by the first stage aligner",default="dot", guitype='comboparambox', choicelist='re_filter_list(dump_cmps_list(),\'tomo\', True)', row=15, col=0, rowspan=1, colspan=3, mode="refinement")
	parser.add_argument("--simralign",type=str,help="The name and parameters of the second stage aligner which refines the results of the first alignment, currently 'refine' or not specified.", default=None, guitype='comboparambox', choicelist='re_filter_list(dump_aligners_list(),\'refine|3d\', 1)', row=16, col=0, rowspan=1, colspan=3, mode="refinement")
	parser.add_argument("--simraligncmp",type=str,help="The name and parameters of the comparator used by the second stage aligner. Default is dot.",default="dot", guitype='comboparambox', choicelist='re_filter_list(dump_cmps_list(),\'tomo\', True)', row=17, col=0, rowspan=1, colspan=3, mode="refinement")
	parser.add_argument("--simcmp",type=str,help="The name of a comparator to be used in comparing the aligned images", default="frc:zeromask=1:snrweight=1", guitype='comboparambox', choicelist='re_filter_list(dump_cmps_list(),\'tomo\', True)', row=13, col=0, rowspan=1, colspan=3, mode="refinement")
	parser.add_argument("--simmask",type=str,help="A file containing a single 0/1 image to apply as a mask before comparison but after alignment", default=None)
	parser.add_argument("--shrink", dest="shrink", type = int, default=0, help="Optionally shrink the input particles by an integer amount prior to computing similarity scores, for speed purposes. Default=0, no shrinking", guitype='shrinkbox', row=12, col=0, rowspan=1, colspan=1, mode="refinement")
	parser.add_argument("--twostage", dest="twostage", type = int, help="Optionally run a faster 2-stage similarity matrix, ~5-10x faster, generally same accuracy. Value specifies shrink factor for first stage, typ 1-3",default=0, guitype='intbox', row=12, col=1, rowspan=1, colspan=1, mode="refinement")
	parser.add_argument("--prefilt",action="store_true",help="Filter each reference (c) to match the power spectrum of each particle (r) before alignment and comparison",default=False, guitype='boolbox', row=12, col=2, rowspan=1, colspan=1)
	
	# options associated with e2classify.py
	parser.add_argument("--simvec", action="store_true",help="Causes the classification algorithm to use patterns rather than peak values",default=False)
	parser.add_argument("--sep", type=int, help="The number of classes a particle can contribute towards (default is 1)", default=1, guitype='intbox', row=19, col=2, rowspan=1, colspan=1, mode="refinement")
	
	# options associated with e2classaverage.py
	parser.add_header(name="caheader", help='Options below this label are specific to e2classaverage', title="### e2classaverage options ###", row=18, col=0, rowspan=1, colspan=3, mode="refinement")
	parser.add_argument("--classkeep",type=float,help="The fraction of particles to keep in each class, based on the similarity score generated by the --cmp argument.",default=0.8, guitype='floatbox', row=20, col=0, rowspan=1, colspan=2, mode="refinement")
	parser.add_argument("--classkeepsig", default=False, action="store_true", help="Change the keep (\'--keep\') criterion from fraction-based to sigma-based.", guitype='boolbox', row=20, col=2, rowspan=1, colspan=1, mode="refinement")
	parser.add_argument("--classiter", type=int, help="The number of iterations to perform. Default is 1.", default=3, guitype='intbox', row=19, col=0, rowspan=1, colspan=2, mode="refinement")
	parser.add_argument("--classalign",type=str,help="If doing more than one iteration, this is the name and parameters of the 'aligner' used to align particles to the previous class average.", default="rotate_translate_flip", guitype='comboparambox', choicelist='re_filter_list(dump_aligners_list(),\'refine|3d\', 1)', row=24, col=0, rowspan=1, colspan=3, mode="refinement")
	parser.add_argument("--classaligncmp",type=str,help="This is the name and parameters of the comparator used by the first stage aligner  Default is dot.",default="phase", guitype='comboparambox', choicelist='re_filter_list(dump_cmps_list(),\'tomo\', True)', row=25, col=0, rowspan=1, colspan=3, mode="refinement")
	parser.add_argument("--classralign",type=str,help="The second stage aligner which refines the results of the first alignment in class averaging. Either specify 'refine' or omit the option.", default=None, guitype='comboparambox', choicelist='re_filter_list(dump_aligners_list(),\'refine|3d\', 1)', row=26, col=0, rowspan=1, colspan=3, mode="refinement")
	parser.add_argument("--classraligncmp",type=str,help="The comparator used by the second stage aligner in class averageing. Default is ccc.",default="ccc", guitype='comboparambox', choicelist='re_filter_list(dump_cmps_list(),\'tomo\', True)', row=27, col=0, rowspan=1, colspan=3, mode="refinement")
	parser.add_argument("--classaverager",type=str,help="The averager used to generate the class averages. Default is \'mean\'.",default="mean", guitype='combobox', choicelist='dump_averagers_list()', row=22, col=0, rowspan=1, colspan=2, mode="refinement")
	parser.add_argument("--classcmp",type=str,help="The name and parameters of the comparator used to generate similarity scores, when class averaging. Default is ccc", default="frc:snrweight=1", guitype='comboparambox', choicelist='re_filter_list(dump_cmps_list(),\'tomo\', True)', row=23, col=0, rowspan=1, colspan=3, mode="refinement")
	parser.add_argument("--classnormproc",type=str,default="normalize.edgemean",help="Normalization applied during class averaging", guitype='combobox', choicelist='re_filter_list(dump_processors_list(),\'normalize\')', row=21, col=0, rowspan=1, colspan=3, mode="refinement")
	parser.add_argument("--classrefsf",default=False, action="store_true", help="Use the setsfref option in class averaging to produce better filtered averages.", guitype='boolbox', row=22, col=2, rowspan=1, colspan=1, mode="refinement")
	parser.add_argument("--classautomask",default=False, action="store_true", help="This will apply an automask to the class-average during iterative alignment for better accuracy. The final class averages are unmasked.")
	
	
	#options associated with e2make3d.py
	parser.add_header(name="make3dheader", help='Options below this label are specific to e2make3d', title="### e2make3d options ###", row=28, col=0, rowspan=1, colspan=3, mode="refinement")
	parser.add_argument("--pad", type=int, dest="pad", help="To reduce Fourier artifacts, the model is typically padded by ~25 percent - only applies to Fourier reconstruction", default=0, guitype='intbox', row=30, col=2, rowspan=1, colspan=1, mode="refinement")
	parser.add_argument("--recon", dest="recon", default="fourier", help="Reconstructor to use. See 'e2help.py reconstructors -v' for more information", guitype='combobox', choicelist='dump_reconstructors_list()', row=29, col=0, rowspan=1, colspan=2, mode="refinement")
	parser.add_argument("--m3dkeep", type=float, help="The percentage of slices to keep in e2make3d.py", default=0.8, guitype='floatbox', row=31, col=0, rowspan=1, colspan=1, mode="refinement")
	parser.add_argument("--m3dkeepsig", default=False, action="store_true", help="The standard deviation alternative to the --m3dkeep argument", guitype='boolbox', row=31, col=1, rowspan=1, colspan=1, mode="refinement")
	parser.add_argument("--m3dsetsf", default=False, action="store_true", help="The standard deviation alternative to the --m3dkeep argument", guitype='boolbox', row=31, col=2, rowspan=1, colspan=1, mode="refinement")
	parser.add_argument("--m3dsffile", default=None, type=str, help="If specified, will use the structure factor from specified file rather than project default")
	parser.add_argument("--m3diter", type=int, default=4, help="The number of times the 3D reconstruction should be iterated", guitype='intbox', row=29, col=2, rowspan=1, colspan=1, mode="refinement")
	parser.add_argument("--m3dpreprocess", type=str, default="normalize.edgemean", help="Normalization processor applied before 3D reconstruction", guitype='combobox', choicelist='re_filter_list(dump_processors_list(),\'normalize\')', row=30, col=0, rowspan=1, colspan=2, mode="refinement")
	parser.add_argument("--m3dpostprocess", type=str, default=None, help="Post processor to be applied to the 3D volume once the reconstruction is completed", guitype='comboparambox', choicelist='re_filter_list(dump_processors_list(),\'filter.lowpass|filter.highpass\')', row=32, col=0, rowspan=1, colspan=3, mode="refinement")
	parser.add_argument("--m3dpostprocess2", type=str, default=None, help="A second post processor to be applied to the 3D volume once the reconstruction is completed")
	
	#lowmem!
	parser.add_argument("--lowmem", default=False, action="store_true",help="Make limited use of memory when possible - useful on lower end machines", guitype='boolbox', row=3, col=2, rowspan=1, colspan=1, mode="refinement")
	parser.add_argument("--parallel","-P",type=str,help="Run in parallel, specify type:<option>=<value>:<option>:<value> EX thread:4",default=None, guitype='strbox', row=3, col=0, rowspan=1, colspan=2, mode="refinement")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	
	(options, args) = parser.parse_args()
	
	logid=E2init(sys.argv,options.ppid)
	
	# make directories for the two jobs
	try: os.mkdir(options.path+"_even")
	except: pass
	try: os.mkdir(options.path+"_odd")
	except: pass
	try: s.mkdir(options.path) # we at least need a place to put refinement results
	except: pass

	if options.path==None:
		print "ERROR: --path=refine_xx is a required option"
		sys.exit(1)

	# create the even and odd data sets
	eset,oset=image_eosplit(options.input)

	if options.usefilt!=None :
		efset,ofset==image_eosplit(options.usefilt)

	# Prepare the starting models for each run
	# each model will have different random phases beyond the specified resolution
	print "### Preparing initial models for refinement, phase-randomized at %1.1f A resolution"%options.randomres
	launch_childprocess("e2proc3d.py %s %s_even/initial_model.hdf --process=filter.lowpass.randomphase:cutoff_freq=%1.4f"%(options.model,options.path,1.0/options.randomres))
	launch_childprocess("e2proc3d.py %s %s_odd/initial_model.hdf --process=filter.lowpass.randomphase:cutoff_freq=%1.4f"%(options.model,options.path,1.0/options.randomres))
	
	# Ok, now we're ready to run the actual refinements !
	argv=sys.argv[1:]
	for a in argv:
		if a[:11]=="--randomres" : argv.remove(a)

	iuf=0
	for i,a in enumerate(argv):
		if a[:6]=="--path" : ipath=i
		if a[:7]=="--model" : imodel=i
		if a[:7]=="--input" : iinp=i
		if a[:9]=="--usefilt" : iuf=i

	# run even refinement
	print "### Starting even data refinement"
	argv[ipath]="--path=%s"%(options.path+"_even")
	argv[imodel]="--model=%s_even/initial_model.hdf"%options.path
	argv[iinp]="--input=%s"%eset
	if iuf>0: argv[iuf]="--usefilt=%s"%efset
	launch_childprocess("e2refine.py "+" ".join(argv))

	# run odd refinement
	print "### Starting odd data refinement"
	argv[ipath]="--path=%s"%(options.path+"_odd")
	argv[imodel]="--model=%s_odd/initial_model.hdf"%options.path
	argv[iinp]="--input=%s"%oset
	if iuf>0: argv[iuf]="--usefilt=%s"%ofset
	launch_childprocess("e2refine.py "+" ".join(argv))
	
	# compute convergence results for even odd test
	for i in xrange(options.startiter,options.iter):
		# do a refine alignment of each odd map to the corresponding even map before resolution calc
		try:
			print "aligning iteration %d"%i
			launch_childprocess("e2proc3d.py %s_odd/threed_filt_%02d.hdf tmp1.hdf --alignref=%s_even/threed_filt_%02d.hdf --align=refine_3d"%(options.path,i,options.path,i))
			launch_childprocess("e2proc3d.py %s_odd/threed_filt_%02d.hdf tmp2.hdf --process=xform.flip:axis=z --alignref=%s_even/threed_filt_%02d.hdf --align=refine_3d"%(options.path,i,options.path,i))
		except:
			print "Alignment failed"
			
		# Pick the best handedness
		a=EMData("tmp1.hdf",0)
		b=EMData("tmp2.hdf",0)
		c=EMData("%s_even/threed_filt_%02d.hdf"%(options.path,i))
		ca=c.cmp(a,"ccc")
		cb=c.cmp(b,"ccc")
		if ca<cb :
			print "correct hand detected"
			os.unlink("%s_odd/threed_filt_%02d.hdf"%(options.path,i))
			os.rename("tmp1.hdf","%s_odd/threed_filt_%02d.hdf"%(i,options.path))
			os.unlink("tmp2.hdf")
		else :
			print "handedness flip required"
			os.unlink("%s_odd/threed_filt_%02d.hdf"%(options.path,i))
			os.rename("tmp2.hdf","%s_odd/threed_filt_%02d.hdf"%(i,options.path))
			os.unlink("tmp1.hdf")

		# Compute FSC convergence plot
		com="e2proc3d.py {path}_even/threed_filt_{iter:02d}.hdf {path}/fsc_eo_{iter:02d}.txt --apix={apix} --calcfsc={path}_odd/threed_filt_{iter:02d}.hdf".format(path=options.path,iter=i,apix=apix)
		if ( launch_childprocess(com) != 0 ):
			print "Failed to execute %s" %com
			exit_refine(1,logid)

	# measure resolution curve
	com="e2proc3d.py {path}_even/threed_filt_{iter:02d}.hdf {path}/fsc_gold.txt --apix={apix} --calcfsc={path}_odd/threed_filt_{iter:02d}.hdf".format(path=options.path,iter=options.iter-1,apix=apix)
	if ( launch_childprocess(com) != 0 ):
		print "Failed to execute %s" %com
		exit_refine(1,logid)


	E2end(logid)

if __name__ == "__main__":
    main()
