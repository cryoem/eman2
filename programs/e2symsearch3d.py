#!/usr/bin/env python

#
# Author: John Flanagan Sept 2011 (jfflanag@bcm.edu)
# Modified by Jesus Galaz-Montoya (jgalaz@gmail.com)
# Last modification: 19/Feb/2015
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
from e2spt_classaverage import preprocessing


def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options] 
	This program aligns a paricle to its symmetry axis. There are two algorithmic modes. 
	A coarse search followed by simplex minimization (not yet implimented) OR monte carlo course 
	search followed by simplex minimization. The Goal is to align the paricle to its 
	symmetry axis so symmetry can be applied for avergaing and for alignment speed up 
	(it is only necessary to search over one asymmetric unit!
	"""
	
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	
	parser.add_header(name="symsearch3dheader", help="""Options below this label are specific to e2symsearch3d""", title="### e2symsearch3d options ###", row=3, col=0, rowspan=1, colspan=2)
	
	parser.add_argument("--input", dest="input", default='',type=str, help="""The name of input volume or hdf stack of volumes""", guitype='filebox', browser="EMBrowserWidget(withmodal=True,multiselect=False)", row=0, col=0, rowspan=1, colspan=2)
	
	#parser.add_argument("--output", dest="output", default="""e2symsearch3d_OUTPUT.hdf""", type=str, help="The name of the output volume", guitype='strbox', filecheck=False, row=1, col=0, rowspan=1, colspan=2)
	
	parser.add_argument("--ref",type=str,default='',help="""Default=None. If provided and --average is also provided and --keep < 1.0 or --keepsig is specified, 'good particles' will be determined by correlation to --ref.""")
	
	parser.add_argument("--mirror",type=str,default='',help="""Axis across of which to generate a mirrored copy of --ref. All particles will be compared to it in addition to the unmirrored image in --ref if --keepsig is provided or if --keep < 1.0.""")
	
	parser.add_argument("--path",type=str, default='', help="""Name of path for output file""", guitype='strbox', row=2, col=0, rowspan=1, colspan=2)

	parser.add_argument("--plots", action='store_true', default=False,help="""Default=False. Turn this option on to generate a plot of the ccc scores if --average is supplied. Running on a cluster or via ssh remotely might not support plotting.""")

	parser.add_argument("--sym", dest = "sym", default="c1", help = """Specify symmetry -choices are: c<n>, d<n>, h<n>, tet, oct, icos. For asymmetric reconstruction ommit this option or specify c1.""", guitype='symbox', row=4, col=0, rowspan=1, colspan=2)
	
	parser.add_argument("--shrink", dest="shrink", type = int, default=0, help="""Optionally shrink the input particles by an integer amount prior to computing similarity scores. For speed purposes. Default=0, no shrinking""", guitype='shrinkbox', row=5, col=0, rowspan=1, colspan=1)

	parser.add_argument("--mask",type=str,help="""Mask processor applied to particles before alignment. Default is mask.sharp:outer_radius=-2. IF using --clipali, make sure to express outer mask radii as negative pixels from the edge.""", returnNone=True, default="mask.sharp:outer_radius=-2", guitype='comboparambox', choicelist='re_filter_list(dump_processors_list(),\'mask\')', row=11, col=0, rowspan=1, colspan=3)
	
	parser.add_argument("--maskfile",type=str,default='',help="""Mask file (3D IMAGE) applied to particles before alignment. Must be in HDF format. Default is None.""")
	
	parser.add_argument("--normproc",type=str,default='',help="""Normalization processor applied to particles before alignment. Default is to use normalize. If normalize.mask is used, results of the mask option will be passed in automatically. If you want to turn this option off specify \'None\'""")
	
	parser.add_argument("--threshold",default='',type=str,help="""A threshold applied to the subvolumes after normalization. For example, --threshold=threshold.belowtozero:minval=0 makes all negative pixels equal 0, so that they do not contribute to the correlation score.""", guitype='comboparambox', choicelist='re_filter_list(dump_processors_list(),\'filter\')', row=10, col=0, rowspan=1, colspan=3)
	
	parser.add_argument("--preprocess",default='',type=str,help="""Any processor (as in e2proc3d.py) to be applied to each volume prior to COARSE alignment. Not applied to aligned particles before averaging.""", guitype='comboparambox', choicelist='re_filter_list(dump_processors_list(),\'filter\')', row=10, col=0, rowspan=1, colspan=3)
		
	parser.add_argument("--lowpass",type=str,default='',help="""A lowpass filtering processor (from e2proc3d.py; see e2help.py processors) to be applied to each volume prior to COARSE alignment. Not applied to aligned particles before averaging.""", guitype='comboparambox', choicelist='re_filter_list(dump_processors_list(),\'filter\')', row=17, col=0, rowspan=1, colspan=3)
	
	parser.add_argument("--highpass",type=str,default='',help="""A highpass filtering processor (from e2proc3d.py, see e2help.py processors) to be applied to each volume prior to COARSE alignment. Not applied to aligned particles before averaging.""", guitype='comboparambox', choicelist='re_filter_list(dump_processors_list(),\'filter\')', row=18, col=0, rowspan=1, colspan=3)
	
	parser.add_argument("--clipali",type=int,default=0,help="""Boxsize to clip particles as part of preprocessing to speed up alignment. For example, the boxsize of the particles might be 100 pixels, but the particles are only 50 pixels in diameter. Aliasing effects are not always as deleterious for all specimens, and sometimes 2x padding isn't necessary; still, there are some benefits from 'oversampling' the data during averaging; so you might still want an average of size 2x, but perhaps particles in a box of 1.5x are sufficiently good for alignment. In this case, you would supply --clipali=75""")
		
	parser.add_argument("--savepreprocessed",action="store_true", default=False, help="""Default=False. Will save stacks of preprocessed particles (one for coarse alignment and one for fine alignment if preprocessing options are different).""")

	parser.add_argument("--average",action='store_true',default=False,help="""Default=False. If supplied and a stack is provided through --input, the average of the aligned and/or symmetrized stack will also be saved.""")
	
	parser.add_argument("--averager",type=str,default="mean.tomo",help="""Default=mean.tomo. The type of averager used to produce the class average. Default=mean.tomo.""")
	
	parser.add_argument("--keep",type=float,default=1.0,help="""Fraction of particles to include if --average is on, after correlating the particles with the average.""")
	
	parser.add_argument("--keepsig", action="store_true", default=False,help="""Default=False. Causes theoptions.keep argument to be interpreted in standard deviations.""", guitype='boolbox', row=6, col=1, rowspan=1, colspan=1, mode='alignment,breaksym')
	
	parser.add_argument("--iter",type=int,default=1,help="""Default=1. If --keep is different from 1.0 and --average is on, an  average will be computed with all the particles, but then the percent specified by --keep (or determined from --keepsigma) will be kept (the rest thrown away) and a new average will be computed. If --iter > 1, this new average will be compared again against all the particles. The procedure will be repeated for however many iterations --iter is given, or the process will stop automatically if in two consecutive iterations exactly the same particles are kept""") 
	
	parser.add_argument('--subset',type=int,default=0,help="""Number of particles in a subset of particles from the --input stack of particles to run the alignments on.""")
	
	parser.add_argument("--steps", dest="steps", type = int, default=10, help="""Number of steps (for the MC). Default=10.""", guitype='intbox', row=5, col=1, rowspan=1, colspan=1)
	
	parser.add_argument("--symmetrize", default=False, action="store_true", help="""Symmetrize volume after alignment.""", guitype='boolbox', row=6, col=0, rowspan=1, colspan=1)
	
	parser.add_argument("--cmp",type=str,help="""The name of a 'cmp' to be used in comparing the symmtrized object to unsymmetrized""", default="ccc", guitype='comboparambox', choicelist='re_filter_list(dump_cmps_list(),\'tomo\', True)', row=7, col=0, rowspan=1, colspan=2)
	
	parser.add_argument("--parallel","-P",type=str,help="""Run in parallel, specify type:<option>=<value>:<option>:<value>""",default=None, guitype='strbox', row=8, col=0, rowspan=1, colspan=2)
	
	parser.add_argument("--ppid", type=int, help="""Set the PID of the parent process, used for cross platform PPID.""",default=-1)

	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="""verbose level [0-9], higner number means higher level ofoptions.verboseness.""")
	
	parser.add_argument("--nopath",action='store_true',default=False,help="""If supplied, this option will save results in the directory where the command is run. A directory to store the results will not be made.""")
		
	parser.add_argument("--nolog",action='store_true',default=False,help="""If supplied, this option will prevent logging the command run in .eman2log.txt.""")
	
	parser.add_argument("--saveali",action='store_true',default=False,help="""Save the stack of aligned/symmetrized particles.""")
	
	parser.add_argument("--savesteps",action='store_true',default=False,help="""If --iter > 1, save all intermediate averages and intermediate aligned kept stacks.""")
	
	#parser.add_argument("--notmatchimgs",action='store_true',default=False,help="""Default=True. This option prevents applying filter.match.to to one image so that it matches the other's spectral profile during preprocessing for alignment purposes.""")
	
	parser.add_argument("--preavgproc1",type=str,default='',help="""Default=None. A processor (see 'e2help.py processors -v 10' at the command line) to be applied to the raw particle after alignment but before averaging (for example, a threshold to exclude extreme values, or a highphass filter if you have phaseplate data.)""")
	
	parser.add_argument("--preavgproc2",type=str,default='',help="""Default=None. A processor (see 'e2help.py processors -v 10' at the command line) to be applied to the raw particle after alignment but before averaging (for example, a threshold to exclude extreme values, or a highphass filter if you have phaseplate data.)""")

	parser.add_argument("--weighbytiltaxis",type=str,default='',help="""Default=None. A,B, where A is an integer number and B a decimal. A represents the location of the tilt axis in the tomogram in pixels (eg.g, for a 4096x4096xZ tomogram, this value should be 2048), and B is the weight of the particles furthest from the tomogram. For example, --weighbytiltaxis=2048,0.5 means that praticles at the tilt axis (with an x coordinate of 2048) will have a weight of 1.0 during averaging, while the distance in the x coordinates of particles not-on the tilt axis will be used to weigh their contribution to the average, with particles at the edge(0+radius or 4096-radius) weighing 0.5, as specified by the value provided for B.""")
	
	parser.add_argument("--weighbyscore",action='store_true',default=False,help="""Default=False. This option will weigh the contribution of each subtomogram to the average by score/bestscore.""")
	
	parser.add_argument("--reconsiderbad",action='store_true',default=False,help="""This option will trigger the re-alignment of particles excluded from the average if --iter > 1 AND ( --keep < 1.0 OR if --keepsig is provided).""")
	
	
	(options, args) = parser.parse_args()
	
	if not options.input:
		parser.print_help()
		sys.exit(0)
	
	#If no failures up until now, initialize logger
	log = 0
	if not options.nolog:
		logid=E2init(sys.argv,options.ppid)
		log = 1

	#inimodeldir = os.path.join(".",options.path)
	#if not os.access(inimodeldir, os.R_OK):
	#	os.mkdir(options.path)
	
	#Make directory to save results
	from e2spt_classaverage import sptmakepath
	options = sptmakepath(options,'symsearch')
	
	if options.nopath:
		options.path = '.'
	
	
	rootpath =os.getcwd()
	
	if rootpath not in options.path:
		options.path = rootpath + '/' + options.path
	
	from e2spt_classaverage import sptOptionsParser

	options = sptOptionsParser( options )
	
	outputstack = options.path + '/aliptcls.hdf'
	
	rets = symsearcher( options, options.input, options.average )
	final_avg = rets[0]
	results = rets[1]
	
	if options.keep == 1.0 and not options.keepsig and options.iter < 2:	
			final_avg.write_image( options.path + '/final_avg.hdf' , 0)
			if options.iter > 1:
				print """ERROR: --iter > 1 must be accompanied by --keepsing, or by --keep < 1.0"""
				sys.exit(1)
	
	elif options.keep < 1.0 or options.keepsig:
		if options.savesteps:
			final_avg.write_image( options.path + '/final_avg_all.hdf' , 0)
		
		if options.ref:
			ref = EMData( options.ref, 0 )
			refComp( options, outputstack, ref, results, '' )
			
			if options.mirror:
				ref.process_inplace('xform.mirror',{'axis': options.mirror })
				refComp( options, outputstack, ref, results, '_vsMirror')
		else:
			ref2compare =  final_avg
			refComp( options, outputstack, final_avg, results, '')	
		

	if log:
		E2end(logid)
	
	return


def symsearcher( options, infile, averageflag=False, results={} ):
	from EMAN2PAR import EMTaskCustomer	

	avgr = Averagers.get( options.averager[0], options.averager[1 ])
	#results = {}
	scores=[]
	
	outputstack = options.path + '/aliptcls.hdf'
	
	#Determine number of particles in the stack
	n = EMUtil.get_image_count( options.input )
	if options.subset and options.subset < n:
		n = options.subset
	
	sptmultinit=0
	for i in range(n):
	
		print "\nI'll look for symmetry in particle number",i
		#Load particle and make a copy to modify if preprocessing options are specified
		volume = EMData(options.input,i)
		preprocvol = volume.copy()
		
		#Preprocess volume if any preprocessing options are specified
		if (options.shrink and options.shrink > 1) or options.mask or options.maskfile or options.lowpass or options.highpass or options.normproc or options.preprocess or options.threshold or options.clipali:
			print "\nHowever, I will first preprocess particle number",i
			
			print "\nWill call preprocessing on ptcl",i
			preprocvol = preprocessing(preprocvol,options,i)
			#preprocessing(s2image,options, ptclindx, savetagp ,'no',round)
			
			print "\nDone preprocessing on ptcl",i
		
		if options.parallel :
			etc=EMTaskCustomer(options.parallel)
		else:
			etc=EMTaskCustomer("thread:1")
		
		symalgorithm = SymALignStrategy( preprocvol, options.sym, options.steps, options.cmp, etc)
		ret = symalgorithm.execute()
		symxform = ret[0]
		score = ret[1]
		
		print "score returned from symalgorithm", score
		
		scores.append( score )
		
		try:																		#This is used when --reconsiderbad is supplied, in which case the 'spt_symsearch_indx' parameter should exist in the header of the "bad" particles.
			results.update( { volume['spt_symsearch_indx']:[symxform,score] } )
			print "writing updated results for low-scoring particle", volume['spt_symsearch_indx']
		except:
			results.update( { i:[symxform,score] } )
	
		print "\nwriting output for best alignment found for particle number",i
		
		if options.shrink and options.shrink > 1:
			trans = symxform.get_trans()
			symxform.set_trans(trans[0]*options.shrink, trans[1]*options.shrink, trans[2]*options.shrink)	
				
		print "\nWrittng to output ptcl",i
	
		#Rotate volume to the best orientation found, set the orientation in the header, apply symmetry if specified and write out the aligned (and symmetrized) particle to the output stack
		output = volume.process('xform',{'transform':symxform})
		output.set_attr('symxform', symxform)
		print "\nApplying this transform to particle",symxform
		if options.symmetrize:
			output = output.process('xform.applysym',{'sym':options.sym})
			
		try:																	#This is used when --reconsiderbad is supplied, in which case the 'spt_symsearch_indx' parameter should exist in the header of the "bad" particles.
			output.write_image( outputstack, volume['spt_symsearch_indx'] )
			print "writing reconsidered low-scoring particle %d to aliptcls.hdf stack" %( volume['spt_symsearch_indx'] )
		except:
			output.write_image( outputstack, -1)
		
		#Averaging here only makes sense if all particles are going to be kept. Otherwise, different code is needed (below)
		if averageflag:
			
			weight = 1.0
			if options.weighbyscore:
				scoreweight = float( score ) / max( scores )
				print "the score weight is %f because score was %f and the best score was %f" % (scoreweight, score, max(scores) )
				weight = weight * scoreweight
			
			if weight < 1.0:
				ptcl.mult( weight )
		
			avgr.add_image( output )
	
			sptmultinit+=1
			
	#Finalize average of all particles. Then, if --keep or --keepsig provided, determine the discrimination threshold and then average the particles that pass it.
	if averageflag: 
			
		final_avg = avgr.finish()
		final_avg['spt_multiplicity']=sptmultinit
		final_avg['origin_x']=0
		final_avg['origin_y']=0		#The origin needs to be reset to ZERO to avoid display issues in Chimera
		final_avg['origin_z']=0
		final_avg['xform.align3d'] = Transform()
	
		return [ final_avg, results ]
	else:
		return results
	
	
def refComp( options, outputstack, ref2compare, results, mirrortag ):
	
	lastexcluded=[]
	meanscores = []
	
	stop = 0
	for it in range( options.iter ):
		#print "Averaging iteration", it
				
		ret = calcScores( outputstack, ref2compare, results )
		scores = ret[0]
		results = ret[1]
		print "in averaging iteration %d scores are" % ( it )
		print scores
		
		meanscore = sum(scores)/len(scores)
		print "for it %d, avg mean score is %.3f" %( it, meanscore )
		meanscores.append( meanscore )
		
		retm = makeSsaAverage( options, scores, results, it )
		ref2compare = retm[0]
		excluded = retm[-1]
	
		if lastexcluded == excluded:
			stop = 1
	
		if options.savesteps and it < options.iter -1 and not stop:
		
			ref2compare.write_image( options.path + '/avgs' + mirrortag + '.hdf', -1)
		
			it == options.iter -1
		
			cmd = 'e2proc3d.py ' + options.path + '/avgs' + mirrortag + '.hdf ' + options.path + '/final_avg' + mirrortag + '.hdf --first -1 --last -1'
		
		elif it == options.iter -1 or stop:
			print "Final mean score, in iteration %d, is %.4f" % ( it, meanscore )
			ref2compare.write_image( options.path + '/final_avg' + mirrortag + '.hdf', 0)
	
		print "\n\nin it %d, excluded is" %( it )
		print excluded
		print "lastexcluded is" 
		print lastexcluded
		if stop:
			print "TERMINATING. The algorithm has converged"
			break;
			
		lastexcluded = list(excluded)	
	
		if options.reconsiderbad:
			print "reconsidering bad particles!"
			badptclsfile = options.path + '/tmpbadptcls.hdf' 
			results = symsearcher( options, badptclsfile, False )
	
	
	lines=[]
	for s in meanscores:
		line = str(s) + '\n'
		lines.append( line )
	
	print "meanscores are", lines
	
	mscoresfile = options.path + '/meanscores.txt' 
	f = open( mscoresfile,'w' )
	f.writelines(lines)
	f.close()
	
	print "written to", mscoresfile
	return


def calcScores( stack, avg, results):
	
	newresults = {}
	
	scores = []
	print "Stack is", stack
	n = EMUtil.get_image_count( stack )

	for i in results:
	
		indx = i 
		#results[r][-1]
	
		img = EMData( stack, indx )
	
		ccmap = avg.calc_ccf( img )
		ccmap.process_inplace('normalize')
		maxloc = ccmap.calc_max_location()
		
		scoredefault = ccmap.get_value_at(0,0,0)
		
		score = ccmap.get_value_at( maxloc[0], maxloc[1], maxloc[2] )
		
		print "scoredefault %.3f and score %.3f" %( scoredefault, score )
		scores.append( score )
		
		t = results[i][0]
		
		newresults.update( { i:[t,score] } )
		
	return [ scores, newresults ]


def makeSsaAverage( options, scores, results, it ):
	
	includedptcls = []
	excludedptcls = []
	
	thresh = None
	
	if options.keep < 1.0 and options.average:
		print "len of scores is", len(scores)
	
		vals=[]
		for p in scores:
			vals.append( p )
	
		vals.sort()
		thresh = vals[ int(options.keep * len(vals) ) - 1]
		if options.verbose: 
			print "Keep threshold : %f (min=%f  max=%f)"%(thresh,vals[0],vals[-1])

	if options.keepsig and options.average:		
		vals=[]
		vals2=[]
		for p in scores:
			vals.append( p )
			vals2.append(p ** 2)
		
		val = sum(vals)
		val2 = sum(vals2)

		mean = val/len( scores )
		sig = sqrt(val2/len( scores )-mean*mean )
		thresh = mean+sig* options.keep
		if options.verbose: 
			print "\nKeep threshold : %f (mean=%f  sigma=%f)"%(thresh,mean,sig)	
	
	ssaavg = None
	if thresh :	
															
		avgr = Averagers.get( options.averager[0], options.averager[1] )
	
		print "Threshold is", thresh
		print "and scores are", scores
		
		sptmult = 0
		for r in results:
			indx =  r
			score = results[r][-1]
			if score > thresh:
				includedptcls.append( indx )
				print "in iteration %d ptcl %d KEPT because its score %.4f was above the threshold %.4f" %(it, indx, score, thresh)
				
				a = EMData( options.input, indx )
			
				if it == options.iter -1:
					a.write_image( options.path + '/keptPtclsRaw_' + str(it).zfill( len( str( options.iter))) + '.hdf', -1 )
			
				t = results[ r ][0]
				a.transform( t )
				#a = a.process('xform',{'transform':symxform})
				
				#print "adding image", a, type(a)
				avgr.add_image( a )
				
				sptmult+=1
				
				if options.symmetrize:
					b = a.copy()
					b.process_inplace('xform.applysym',{'sym':options.sym})
			
					if options.saveali and it == options.iter -1:
						b.write_image( options.path + '/keptPtclsAli_' + str(it).zfill( len( str( options.iter))) + '.hdf', -1 )
			else:
				excludedptcls.append( indx )
				print "in iteration %d ptcl %d DISCARDED because its score %.4f was under the threshold %.4f" % ( it, indx, score, thresh )

				if options.reconsiderbad:
					c = EMData( options.input, indx )
					c['spt_symsearch_indx'] = indx
					badptcls = options.path + '/tmpbad.hdf'
					c.write_image( badptcls, -1 )
					print "however, it will be reconsidered"	
			
		ssaavg = avgr.finish()
		#print "finalizing average", ssaavg, type(ssaavg)
		if ssaavg:
			ssaavg['origin_x']=0
			ssaavg['origin_y']=0		#The origin needs to be reset to ZERO to avoid display issues in Chimera
			ssaavg['origin_z']=0
			ssaavg['xform.align3d'] = Transform()
			ssaavg['spt_multiplicity'] = sptmult

			if options.symmetrize:
				ssaavg.process_inplace('xform.applysym',{'sym':options.sym})
			
			print "returning ssaavg"
			return [ssaavg,includedptcls,excludedptcls]
		else:
			print "ERROR! ssaavg = None!!!"
			sys.exit()
		#else:
		#	print "ERROR: in iteration it %d none of the particles had a score higher than the threshold %.4f" %(it, thresh)
		#	print "scores", scores
		#	sys.exit()

	else:
		print "ERROR: threshold=None."
		sys.exit(1)

		return None
	


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
				
		return [bestxform, bestscore]
		
		
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
    
    
    
"""
PROGRAM SCRAPS

	'''
		output = None
		if preproc:
			trans = symxform.get_trans()
			symxform.set_trans(trans[0]*options.shrink, trans[1]*options.shrink, trans[2]*options.shrink)
			output = EMData(options.input,i)
			output.process_inplace('xform',{'transform':symxform})
			print "\nApplying this transform to particle",symxform
			output.set_attr('symxform', symxform)	# Only HDF files will contain this metadata
			if options.symmetrize:
				output = output.process('xform.applysym',{'sym':options.sym})
		else:
			output = volume.process('xform',{'transform':symxform})
			output.set_attr('symxform', symxform)	# Only HDF files will contain this metadata
			print "\nApplying this transform to particle",symxform
			if options.symmetrize:
				output = output.process('xform.applysym',{'sym':options.sym})
		'''

"""
