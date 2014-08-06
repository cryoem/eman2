#!/usr/bin/env python

#
# Author: John Flanagan Sept 2011 (jfflanag@bcm.edu); last update by Jesus Galaz-Montoya on March/20/2014
# Modified by Jesus Galaz-Montoya to enable iteration over particle stacks and particle preprocessing
# Copyright (c) 2000-2011 Baylor College of Medicine
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


from EMAN2 import *
import math
import os
from EMAN2jsondb import JSTask,jsonclasses
import sys

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options] 
	This program aligns a paricle to its symmetry axis. There are two algorithmic modes. A coarse search followed by simplex
	minimization(not yet implimented) OR monte carlo course search followed by simplex minimization. The Goal is to align the paricle to its 
	symmetry axis so symmetry can be applied for avergaing and for alignment speed up(it is only necessary to search over the
	assymetric unit!
	"""
	
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	
	parser.add_header(name="refineheader", help='Options below this label are specific to e2refine', title="### e2refine options ###", row=3, col=0, rowspan=1, colspan=2, mode="align")
	parser.add_argument("--input", dest="input", default=None,type=str, help="The name of input volume or hdf stack of volumes", guitype='filebox', browser="EMBrowserWidget(withmodal=True,multiselect=False)", row=0, col=0, rowspan=1, colspan=2, mode="align")
	parser.add_argument("--output", dest="output", default="e2symsearch3d_OUTPUT.hdf",type=str, help="The name of the output volume", guitype='strbox', filecheck=False, row=1, col=0, rowspan=1, colspan=2, mode="align")
	parser.add_argument("--path",type=str,help="Name of path for output file",default='', guitype='strbox', row=2, col=0, rowspan=1, colspan=2, mode="align['initial_models']")
	
	
	
	parser.add_argument("--sym", dest = "sym", default="c1", help = "Specify symmetry - choices are: c<n>, d<n>, h<n>, tet, oct, icos. For asymmetric reconstruction omit this option or specify c1.", guitype='symbox', row=4, col=0, rowspan=1, colspan=2, mode="align")
	
	
	parser.add_argument("--shrink", dest="shrink", type = int, default=0, help="Optionally shrink the input particles by an integer amount prior to computing similarity scores. For speed purposes. Default=0, no shrinking", guitype='shrinkbox', row=5, col=0, rowspan=1, colspan=1, mode="align")

	parser.add_argument("--mask",type=str,help="""Mask processor applied to particles before alignment. 
		Default is mask.sharp:outer_radius=-2. IF using --clipali, make sure to express outer mask radii as negative 
		pixels from the edge.""", returnNone=True, default="mask.sharp:outer_radius=-2", guitype='comboparambox', choicelist='re_filter_list(dump_processors_list(),\'mask\')', row=11, col=0, rowspan=1, colspan=3, mode='alignment,breaksym')
	
	
	parser.add_argument("--maskfile",type=str,default=None,help="""Mask file (3D IMAGE) applied to particles 
		before alignment. Must be in HDF format. Default is None.""")
	
	parser.add_argument("--normproc",type=str,default='',help="Normalization processor applied to particles before alignment. Default is to use normalize. If normalize.mask is used, results of the mask option will be passed in automatically. If you want to turn this option off specify \'None\'")
	
	parser.add_argument("--threshold",default='',type=str,help="""A threshold applied to the subvolumes after normalization. 
													For example, --threshold=threshold.belowtozero:minval=0 makes all negative pixels equal 0, so that they do not contribute to the correlation score.""", guitype='comboparambox', choicelist='re_filter_list(dump_processors_list(),\'filter\')', row=10, col=0, rowspan=1, colspan=3, mode='alignment,breaksym')
	
	parser.add_argument("--preprocess",default='',type=str,help="Any processor (as in e2proc3d.py) to be applied to each volume prior to COARSE alignment. Not applied to aligned particles before averaging.", guitype='comboparambox', choicelist='re_filter_list(dump_processors_list(),\'filter\')', row=10, col=0, rowspan=1, colspan=3, mode='alignment,breaksym')
	#parser.add_argument("--preprocessfine",type=str,help="Any processor (as in e2proc3d.py) to be applied to each volume prior to FINE alignment. Not applied to aligned particles before averaging.", default=None)
	
	parser.add_argument("--lowpass",type=str,default='',help="A lowpass filtering processor (as in e2proc3d.py) to be applied to each volume prior to COARSE alignment. Not applied to aligned particles before averaging.", guitype='comboparambox', choicelist='re_filter_list(dump_processors_list(),\'filter\')', row=17, col=0, rowspan=1, colspan=3, mode='alignment,breaksym')
	#parser.add_argument("--lowpassfine",type=str,help="A lowpass filtering processor (as in e2proc3d.py) to be applied to each volume prior to FINE alignment. Not applied to aligned particles before averaging.", default=None)

	parser.add_argument("--highpass",type=str,default='',help="A highpass filtering processor (as in e2proc3d.py) to be applied to each volume prior to COARSE alignment. Not applied to aligned particles before averaging.", guitype='comboparambox', choicelist='re_filter_list(dump_processors_list(),\'filter\')', row=18, col=0, rowspan=1, colspan=3, mode='alignment,breaksym')
	#parser.add_argument("--highpassfine",type=str,help="A highpass filtering processor (as in e2proc3d.py) to be applied to each volume prior to FINE alignment. Not applied to aligned particles before averaging.", default=None)

	parser.add_argument("--clipali",type=int,default=0,help="""Boxsize to clip particles as part of preprocessing
		to speed up alignment. For example, the boxsize of the particles might be 100 pixels, but the particles are only 50 pixels 
		in diameter. Aliasing effects are not always as deleterious for all specimens, and sometimes 2x padding isn't necessary;
		still, there are some benefits from 'oversampling' the data during averaging; so you might still want an average of size
		2x, but perhaps particles in a box of 1.5x are sufficiently good for alignment. In this case, you would supply --clipali=75""")
	
	#parser.add_argument("--postprocess",type=str,help="A processor to be applied to the FINAL volume after averaging the raw volumes in their FINAL orientations, after all iterations are done.",default=None, guitype='comboparambox', choicelist='re_filter_list(dump_processors_list(),\'filter\')', row=16, col=0, rowspan=1, colspan=3, mode='alignment,breaksym')
	
	
	parser.add_argument("--savepreprocessed",action="store_true", help="Will save stacks of preprocessed particles (one for coarse alignment and one for fine alignment if preprocessing options are different).", default=False)

	
	
	parser.add_argument("--steps", dest="steps", type = int, default=10, help="Number of steps (for the MC)", guitype='intbox', row=5, col=1, rowspan=1, colspan=1, mode="align")
	parser.add_argument("--symmetrize", default=False, action="store_true", help="Symmetrize volume after alignment.", guitype='boolbox', row=6, col=0, rowspan=1, colspan=1, mode="align")
	parser.add_argument("--cmp",type=str,help="The name of a 'cmp' to be used in comparing the symmtrized object to unsymmetrized", default="ccc", guitype='comboparambox', choicelist='re_filter_list(dump_cmps_list(),\'tomo\', True)', row=7, col=0, rowspan=1, colspan=2, mode="align")
	parser.add_argument("--parallel","-P",type=str,help="Run in parallel, specify type:<option>=<value>:<option>:<value>",default=None, guitype='strbox', row=8, col=0, rowspan=1, colspan=2, mode="align")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)


	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n",type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness.")


	#parser.add_argument("--mask",type=str,help="Mask processor applied to particles before alignment. Default is mask.sharp:outer_radius=-2", default="mask.sharp:outer_radius=-2")
	

	(options, args) = parser.parse_args()


	
	if ".hdf" not in options.output and ".mrc" not in options.output:
		print "ERROR. The output must contain a valid format ending, for example '.hdf.' TERMINATING!"
		sys.exit()
	
	if not options.input:
		parser.print_help()
		sys.exit(0)
	
	#If no failures up until now, initialize logger
	logid=E2init(sys.argv,options.ppid)
	
	#inimodeldir = os.path.join(".",options.path)
	#if not os.access(inimodeldir, os.R_OK):
	#	os.mkdir(options.path)
	
	#Make directory to save results
	from e2spt_classaverage import sptmakepath
	options = sptmakepath(options,'symsearch')
	
	#Import preprocessing function
	from e2spt_classaverage import preprocessing	
	
	#Import parallelization class
	from EMAN2PAR import EMTaskCustomer
	
	#Determine number of particles in the stack
	n = EMUtil.get_image_count( options.input )

	#Iterate over particles in stack
	for i in range(n):
	
		print "\nI'll look for symmetry in particle number",i
		#Load particle
		volume = EMData(options.input,i)
		
		preproc = 0
		preprocvol = volume.copy()
		
		#Preprocess volume if any preprocessing options are specified
		if options.shrink or options.mask or options.lowpass or options.highpass or options.normproc or options.preprocess or options.threshold or options.clipali:
			print "\nHowever, I will first preprocess particle number",i
			#Parse preprocessing parameters in a way that the preprocessing function will be able to apply them
			
			if options.normproc and i == 0: 
				print "parsing normproc"
				options.normproc=parsemodopt(options.normproc)
	
			if options.mask and i == 0: 
				print "parsing mask"

				options.mask=parsemodopt(options.mask)
	
			if options.preprocess and i == 0:
				print "parsing preproc"
 
				options.preprocess=parsemodopt(options.preprocess)
		
			if options.threshold and i == 0: 
				print "parsing thresh"

				options.threshold=parsemodopt(options.threshold)

			if options.lowpass and i == 0:
				print "parsing lowpass"
 
				options.lowpass=parsemodopt(options.lowpass)
	
			if options.highpass and i == 0: 
				print "parsing highpass"

				options.highpass=parsemodopt(options.highpass)
		
			
			print "\nWill call preprocessing on ptcl",i
			#preprocvol = preprocessing(options, preprocvol)		
			preprocvol = preprocessing(preprocvol,options,options.mask,options.clipali,options.normproc,options.shrink,options.lowpass,options.highpass,options.preprocess,options.threshold,i)
			
			print "\nDone preprocessing on ptcl",i
			preproc = 1
			
			#Save the preprocessed volume for inspection if desired
			#if options.savepreprocessed:
			#	print "Saving preproc ptcl",i
			#	outpreproc = options.path + '/' + options.input.replace('.hdf','_preproc.hdf')
			#	preprocvol.write_image(outpreproc,-1)
				
		#Initialize parallelism if being used
		
		if options.parallel :
			etc=EMTaskCustomer(options.parallel)
		else:
			etc=EMTaskCustomer("thread:1")
		
		symalgorithm = SymALignStrategy( preprocvol, options.sym, options.steps, options.cmp, etc)
		symxform = symalgorithm.execute()
	
		print "\nWriting out put for best alignment found for particle number",i
		
		output = None
		if preproc:
			trans = symxform.get_trans()
			symxform.set_trans(trans[0]*options.shrink, trans[1]*options.shrink, trans[2]*options.shrink)
			output = EMData(options.input,i)
			output.process_inplace('xform',{'transform':symxform})
			print "\nApplying this transform to particle",symxform
			output.set_attr('symxform', symxform)	# Obviously only HDF or BDB files will contain this metadata
			if options.symmetrize:
				output = output.process('xform.applysym',{'sym':options.sym})
		else:
			output = volume.process('xform',{'transform':symxform})
			output.set_attr('symxform', symxform)	# Obviously only HDF or BDB files will contain this metadata
			print "\nApplying this transform to particle",symxform
			if options.symmetrize:
				output = output.process('xform.applysym',{'sym':options.sym})
			
		#if inimodeldir =='./' or inimodeldir.upper == './NONE': # That is to say no directory is wanted (--path='')
		#	output.write_image(options.output)
		#else:
		#	print "WHT%sZZZ"%inimodeldir
		#	output.write_image(os.path.join(inimodeldir, options.output))
		
			
		if output:
			print "\nWrittng to output ptcl",i
			output.write_image(options.path + '/' + options.output,-1)
	
	
	
	
	
	
	"""	
	if options.mask:
		volume.mult(mask)
	if options.shrink:
		volume.process_inplace("math.meanshrink",{"n":options.shrink})
		
	symalgorithm = SymALignStrategy(volume, options.sym, options.steps, options.cmp, etc)
	symxform = symalgorithm.execute()
	
	print "Found Best alignment, now writting output..."
	if options.shrink:
		trans = symxform.get_trans()
		symxform.set_trans(trans[0]*options.shrink, trans[1]*options.shrink, trans[2]*options.shrink)
		output = EMData(options.input)
		output.process_inplace('xform',{'transform':symxform})
		output.set_attr('symxform', symxform)	# Obviously only HDF or BDB files will contain this metadata
		if options.symmetrize:
			print "symmetrize output"
			output = output.process('xform.applysym',{'sym':options.sym})
	else:
		output = volume.process('xform',{'transform':symxform})
		output.set_attr('symxform', symxform)	# Obviously only HDF or BDB files will contain this metadata
		if options.symmetrize:
			print "symmetrize output"
			output = output.process('xform.applysym',{'sym':options.sym})
			
	if inimodeldir =='./' or inimodeldir.upper == './NONE': # That is to say no directory is wanted (--path='')
		output.write_image(options.output)
	else:
		print "WHT%sZZZ"%inimodeldir
		output.write_image(os.path.join(inimodeldir, options.output))
	"""
	
	E2end(logid)
	
	return


# Use strategy pattern here. Any new stategy needs to inherit this
class Strategy:
	def __init__(self, volume, sym, steps, comp, etc):
		self.volume = volume
		self.sym = sym
		self.steps = steps
		self.cmp=parsemodopt(comp)
		self.etc = etc
		
	def execute(self):
		raise NotImplementedError("Subclass must implement abstract method")


class SymALignStrategy(Strategy):
	""" MC followed by minimization """
	def __init__(self, volume, sym, steps, comp, etc):
		Strategy.__init__(self, volume, sym, steps, comp, etc)
		
	def execute(self):
		Util.set_randnum_seed(Util.get_randnum_seed())
		tasks=[]
		for i in xrange(self.steps):
			az = Util.get_frand(0,360) 					# theta
			alt  = math.degrees(math.acos(2*Util.get_frand(0,1) - 1))	# phi
			phi = Util.get_frand(0,360)					# kappa
			t = Transform({"type":"eman","az":az,"alt":alt,"phi":phi})
			# Now do the simplex alignment
			tasks.append(SymAlignTask(self.volume, self.sym, self.cmp, t))
		tids=self.etc.send_tasks(tasks)
		
		solns = []
		# Hang till tasks are done
		while 1:
			time.sleep(5)
			proglist=self.etc.check_task(tids)
			
			for i,prog in enumerate(proglist):
				if prog==100:
					print "Finished MC trial number", i
					r=self.etc.get_results(tids[i])
					solns.append(r[1]["symalign"])
					
			tids=[j for i,j in enumerate(tids) if proglist[i]!=100]		# remove any completed tasks from the list we ask about
			if len(tids)==0: break
		
		print "\nFinished alignments...\n"
		# Search for the best scoring
		bestscore = 0
		bestxform = Transform()
		for i,soln in enumerate(solns):
			print "score for MC trial %d is %f"%(i,soln.get_attr('score'))
			if soln.get_attr('score') < bestscore:
				bestscore = soln.get_attr('score')
				bestxform = soln.get_attr('xform.align3d')
				
		return bestxform
		
		
class SymAlignTask(JSTask):
	def __init__(self, volume, sym, comp, xform):
		data = {"volume":volume}
		JSTask.__init__(self,"CmpTilt",data,{},"")
		
		self.sym = sym
		self.cmp=comp
		self.xform=xform
		
	def execute(self,callback=None):
		symalign = self.data['volume'].align('symalignquat',self.data['volume'],{"sym":self.sym,"xform.align3d":self.xform},self.cmp[0],self.cmp[1])
		return {"symalign":symalign}
	
jsonclasses["SymAlignTask"]=SymAlignTask.from_jsondict

if __name__ == "__main__":
    main()
