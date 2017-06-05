#!/usr/bin/env python
#
# Author: Jesus Galaz, 22/Jan/2017; last update Jan/22/2017
# Copyright (c) 2011 Baylor College of Medicine
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

import os
from EMAN2 import *

import sys
import numpy
import math
import collections

#from e2spt_intrafsc import genOddAndEvenVols, fscOddVsEven	
	
def main():
	
	progname = os.path.basename(sys.argv[0])
	usage = """run e2spt_bfactorplot.py --inputstack --ref <other_options>, to see whether resolution against a known structure plateaus. to see whether 
	gold-standard resolution has plateaued, run e2spt_bfactorplot.py --inputeven --inputodd <other_options>. this program helps determine whether resolution for
	subtomgoram alignment is particle-limited, or whetehr you should add more particles"""
			
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	
	parser.add_argument("--apix",type=float,default=0.0, help="""Default=0.0 (not used). Use this apix value where relevant instead of whatever is in the header of the reference and the particles.""")

	#parser.add_argument("--automask",action='store_true',default=False,help="""Applies loose automask at threshold=2.0""")
	
	parser.add_argument("--averager",type=str,default="mean.tomo", help="""Default=mean.tomo. The type of averager used to produce the class average.""")

	parser.add_argument("--clip",type=float,default=0.0,help="""Default=0.0 (not used). Size to clip the box to before calculating FSCs. If particle of 64 pixels in diameter are in a 128x128x128 box, it's safe to provide --clip 70 or so""")

	#parser.add_argument("--fsccutoff",type=float,default=0.143,help="""Default=0.0 (not used).""")

	parser.add_argument("--inputeven", type=str, default='', help="""Default=None. 'EVEN' stack of aligned subvolumes after gold-standard SPT refinement. MUST be HDF format, since volume stack support is required.""")
	parser.add_argument("--inputodd", type=str, default='', help="""Default=None. 'ODD' stack of aligned subvolumes after gold-standard SPT refinement. MUST be HDF format, since volume stack support is required.""")
	
	parser.add_argument("--inputstack", type=str, default='', help="""Default=None. stack aligned subvolumes. MUST be HDF format, since volume stack support is required. needs to be in the same orientation/symmetry-axis as --ref, and have the same apix and box size, otherwise the results will be invalid.""")

	parser.add_argument("--mask1", type=str,default='', help="""Default=None. Mask processor to apply to averages before FSC computation.""")
	parser.add_argument("--mask2", type=str,default='', help="""Default=None. Mask processor to apply to averages before FSC computation.""")
	parser.add_argument("--maskfile",type=str,default='', help="""Default=None. Mask file to multiply the averages by before FSC computation.""")

	parser.add_argument("--path",type=str,default='sptbfactor', help="""Default=spt. Directory to store results in. The default is a numbered series of directories containing the prefix 'spt'; for example, spt_02 will be the directory by default if 'spt_01' already exists.""")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)

	parser.add_argument("--ref", type=str, default='', help="""Default=None. known structure to compare --inputstack against. needs to be in the same orientation/symmetry-axis as --inputstack, and have the same box size and apix, otherwise the results will be invalid.""")

	parser.add_argument("--runningavg",action='store_true', default=False, help="""Computes running average of particles weighing properly, instead of computing each average (as N grows) from scratch.""")


	parser.add_argument("--step", type=int, default=1, help="""Number of particles to increase from one data point to the next. For example, --step=10 will compute the B-factor averaging 10 particles from the even set and 10 from the odd set; then 20; then 40; etc.""")
	parser.add_argument("--sym", type=str, default='', help="""Default=None (equivalent to c1). Symmetry to impose -choices are: c<n>, d<n>, h<n>, tet, oct, icos""")
	
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n",type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	
	(options, args) = parser.parse_args()	#c:this parses the options or "arguments" listed 
											#c:above so that they're accesible in the form of option.argument; 
											#c:for example, the input for --template would be accesible as options.template
		
	logger = E2init(sys.argv, options.ppid)	#c:this initiates the EMAN2 logger such that the execution
											#of the program with the specified parameters will be logged
											#(written) in the invisible file .eman2log.txt
	
	if options.inputstack:
		if options.inputeven or options.inputodd:
			print "\nERROR: cannot supply --inputodd or --inputeven with --inputstack"
			sys.exit(1)
		if not options.ref:
			print "\nERROR: --ref required with --inputstack"

	if options.inputodd and not options.inputeven:
		print "\nERROR: --inputodd requires --inputeven"

	if options.inputeven and not options.inputodd:
		print "\nERROR: --inputeven requires --inputodd"

	if options.ref:
		if not options.inputstack:
			print "\nERROR: --inputstack required with --ref"
		if options.inputeven or options.inputodd:
			print "\nERROR: do not supply --inputeven nor --inputodd with --ref and --inputstack; these are mutually exclusive options."

	if options.averager:
		options.averager=parsemodopt(options.averager)
	
	if options.mask1:
		options.mask1=parsemodopt(options.mask1)
	
	if options.mask2:
		options.mask2=parsemodopt(options.mask2)
	
	from e2spt_classaverage import sptmakepath
	options = sptmakepath(options,'sptbfactor')


	if options.inputeven and options.inputodd:
		bfactorfuncgold(options)

	elif options.inputstack and options.inputref:
		bfactorfunc(options)
	
	'''	
	ccc
	ccc.tomo.thresh
	ccc.tomo
	fsc.tomo.auto
	'''
	
	return


def averagerfunc(stack,options,n):
	print "\ncomputing average %d, from file %s" %(n,stack)
	
	avgr = Averagers.get( options.averager[0], options.averager[1])
	print '\navgr is',avgr
	print "\nn is",n

	for x in range(n):
		print "\nadded ptcl %d/%d" %(x,n)
		ptcl = EMData( stack, x )
		
		print "\napix is",ptcl['apix_x']
		avgr.add_image( ptcl )
	
	avg = avgr.finish()
	if avg:
		#avg.process_inplace('normalize.edgemean')

		return avg
	else:
		print "average failed"
		sys.exit(1)

		return


def runningavg(stack,previousavg,options,start,stop,count):
	print "\ncomputing running avg %d, from file %s" %(count,stack)
	
	avgr = Averagers.get( options.averager[0], options.averager[1])
		
	for i in xrange(start,stop,1):
		ptcl = EMData( stack, i )
		#ptcl.process_inplace('normalize.edgemean')	
	
		avgr.add_image( ptcl )
	
	#previousavg.process_inplace('normalize.edgemean')
	
	previousavg.mult(count+1)
	avgr.add_image( previousavg )
	
	avg = avgr.finish()
	#if avg:
	#avg.process_inplace('normalize.edgemean')
		
	return avg
	#else:
	#	print "average failed"
	#	return

	#return avg

	
def calcfsc( options, img1, img2, gold=True ):
	
	img1fsc = img1.copy()
	img2fsc = img2.copy()
	
	apix = img1['apix_x']
	if options.apix:
		apix=options.apix
	
	#if options.clip:
		#img1fsc = clip3D( img1fsc, options.clip )
		#img1fsc.process_inpl
		
	#	img1fsc.write_image(options.path +'/vol4fsc1.hdf',0)
		
	#	img2fsc = clip3D( img2fsc, options.clip )
	#	img2fsc.write_image(options.path +'/vol4fsc2.hdf',0)
		
	fsc = img1fsc.calc_fourier_shell_correlation( img2fsc )
	third = len( fsc )/3
	xaxis = fsc[0:third]
	fsc = fsc[third:2*third]
	saxis = [x/apix for x in xaxis]

	fscfile = options.path + '/tmpfsc.txt'
	Util.save_data( saxis[1], saxis[1]-saxis[0], fsc[1:-1], fscfile )

	f=open(fscfile,'r')
	lines=f.readlines()
	fscarea = sum( [ float(line.split()[-1].replace('\n','')) for line in lines ])

	fscvals=[]
	resolutions=[]
	for line in lines:
		linesplit = line.split()
		fscval = float(linesplit[-1].replace('\n',''))
		resolution = float(linesplit[-1].replace('\n',''))
		
		fscvals.append( fscval )
		resolutions.append( resolution )

		fsccutoff = 0.143

		if not gold:
			fsccutoff = 0.5



	f.close()
	
	return fscarea
	

def bfactorfuncgold( options ):
	ne = EMUtil.get_image_count( options.inputeven )
	no = EMUtil.get_image_count( options.inputodd )

	hdr = EMData( options.inputeven, 0, True )
	
	#box = hdr['nx']
	#radius = box/4
	
	#print "\nradius is", radius
	#radius_expanded = radius + math.ceil(0.2*radius)
	
	#print "\nradius expanded is", radius_expanded
	
	nfinal = (ne+no)/2
	if ne < nfinal:
		nfinal = ne
	elif no < nfinal:
		nfinal = no
	
	fscareas = []
	fscareasdict={}
	
	preve = EMData( options.inputeven, 0 )
	prevo = EMData( options.inputodd, 0 )
	
	avge = preve.copy()
	avgo = prevo.copy()

	if options.runningavg:
			
		for i in range(nfinal):
			avgre = Averagers.get( options.averager[0], options.averager[1])
			avgro = Averagers.get( options.averager[0], options.averager[1])

			ptcle = EMData( options.inputeven, i )
			ptclo = EMData( options.inputodd, i )
			
			avgre.add_image( ptcle )
			avgro.add_image( ptclo )

			if i > 0:
				#pavge.process_inplace('normalize.edgemean')
				#pavgo.process_inplace('normalize.edgemean')
				
				#pavge.mult(i)
				#pavgo.mult(i)
				
				avgre.add_image( pavge )
				avgro.add_image( pavgo )

			avge = avgre.finish()
			avgo = avgro.finish()

			pavge = avge.copy()
			pavgo = avgo.copy()
			
			avge_w = avge.copy()
			avgo_w = avgo.copy()

			if options.sym:
				avge_w.process_inplace('xform.applysym',{'sym':options.sym})
				avgo_w.process_inplace('xform.applysym',{'sym':options.sym})
		
			if options.mask1:
				print "parsed mask1 is",options.mask1
				avge_w.process_inplace(options.mask1[0], options.mask1[1])
				avgo_w.process_inplace(options.mask1[0], options.mask1[1])
			
			if options.mask2:
				print "parsed mask2 is",options.mask2
				avge_w.process_inplace(options.mask2[0], options.mask2[1])
				avgo_w.process_inplace(options.mask2[0], options.mask2[1])
		
			if options.maskfile:
				mask = EMData(options.maskfile)
				avge_w.mult(mask)
				avgo_w.mult(mask)
			
			fscarea = calcfsc( options, avge_w, avgo_w, gold=True )

			fscareas.append( fscarea )
			fscareasdict.update({i:fscarea})
			
			f = open( options.path +'/n_vs_fsc.txt','w' )
			fsclines = []
			
			x=0
			sortedfscareasdict = collections.OrderedDict(sorted(fscareasdict.items()))
			
			#for fscarea in fscareas:
			for k in sortedfscareasdict:
				fscareak = sortedfscareasdict[k]
				fscline = str(x) + '\t'+ str(fscareak) + '\n'
				fsclines.append( fscline )
				x+=1
		
			f.writelines( fsclines )
			f.close()
			print "\ndone with fsc %d/%d" %(i,nfinal)

		g = open( options.path +'/n_vs_fsc_unsorted.txt','w' )
		fsclines_unsorted = []
		xx=0
		for area in fscareas:
			fscline = str(xx) + '\t'+ str(area) + '\n'
			fsclines.append( fscline )
			xx+=1
		
		g.writelines( fsclines )
		g.close()

	else:
		ngroups = nfinal/options.step
		excedent = nfinal%options.step
		count=0
		
		h = open( options.path +'/n_vs_fsc_unsorted.txt','a' )

		for thisn in xrange( 0, nfinal, options.step ):
			
			startindx = thisn
			stopindx = options.step*(count+1)
			if count == ngroups -1 :
				stopindx = nfinal
				
			avge = averagerfunc( options.inputeven, options, stopindx )
			avgo = averagerfunc( options.inputodd, options, stopindx )
			
			avge_w = avge.copy()
			avgo_w = avgo.copy()
			if options.sym:
				avge_w.process_inplace('xform.applysym',{'sym':options.sym})
				avgo_w.process_inplace('xform.applysym',{'sym':options.sym})	
			
			if options.mask1:
				print "parsed mask1 is",options.mask1
				avge_w.process_inplace(options.mask1[0], options.mask1[1])
				avgo_w.process_inplace(options.mask1[0], options.mask1[1])
			
			if options.mask2:
				print "parsed mask2 is",options.mask2
				avge_w.process_inplace(options.mask2[0], options.mask2[1])
				avgo_w.process_inplace(options.mask2[0], options.mask2[1])
		
			if options.maskfile:
				mask = EMData( options.maskfile, 0 )
				avge_w.mult(mask)
				avgo_w.mult(mask)
			
			fscarea = calcfsc( options, avge_w, avgo_w, gold=True )

			fscareas.append( fscarea )
			fscareasdict.update({count:fscarea})
			
			f = open( options.path +'/n_vs_fsc.txt','w' )
			fsclines = []
			
			
			x=0
			sortedfscareasdict = collections.OrderedDict(sorted(fscareasdict.items()))
			
			for k in sortedfscareasdict:
				fscareak = sortedfscareasdict[k]
				fscline = str(x) + '\t'+ str(fscareak) + '\n'
				fsclines.append( fscline )

				x+=1
		
			f.writelines( fsclines )
			f.close()
			
			hline = str(count) + '\t'+ str(fscarea) + '\n'

			h.write(hline)


			print "\ndone with fsc %d/%d" %(count,ngroups)

			count+=1

		h.close()

		g = open( options.path +'/n_vs_fsc_unsorted_final.txt','w' )
		fsclines_unsorted = []
		xx=0
		for area in fscareas:
			fscline = str(xx) + '\t'+ str(area) + '\n'
			fsclines.append( fscline )
			xx+=1
		
		g.writelines( fsclines )
		g.close()
	
		print "\ndone with fsc {}/{}".format(count,ngroups)
	
	return


def bfactorfunc( options ):
	n = EMUtil.get_image_count( options.inputstack )

	hdr = EMData( options.inputstack, 0, True )
	
	fscareas = []
	fscareasdict={}
	
	preve = EMData( options.ref, 0 )
	prevo = EMData( options.inputstack, 0 )
	
	avgo = prevo.copy()

	if options.runningavg:
			
		for i in range(n):
			avgro = Averagers.get( options.averager[0], options.averager[1])

			ptclo = EMData( options.inputstack, i )
			
			avgro.add_image( ptclo )

			if i > 0:
				avgro.add_image( pavgo )

			avgo = avgro.finish()

			pavgo = avgo.copy()
			
			avgo_w = avgo.copy()

			if options.sym:
				avgo_w.process_inplace('xform.applysym',{'sym':options.sym})
		
			if options.mask1:
				print "parsed mask1 is",options.mask1
				avgo_w.process_inplace(options.mask1[0], options.mask1[1])
			
			if options.mask2:
				print "parsed mask2 is",options.mask2
				avgo_w.process_inplace(options.mask2[0], options.mask2[1])
		
			if options.maskfile:
				mask = EMData( options.maskfile, 0 )
				avgo_w.mult( mask )
			
			fscarea = calcfsc( options, preve, avgo_w, gold=False )

			fscareas.append( fscarea )
			fscareasdict.update({i:fscarea})
			
			f = open( options.path +'/n_vs_fsc.txt','w' )
			fsclines = []
			
			x=0
			sortedfscareasdict = collections.OrderedDict(sorted(fscareasdict.items()))
			
			#for fscarea in fscareas:
			for k in sortedfscareasdict:
				fscareak = sortedfscareasdict[k]
				fscline = str(x) + '\t'+ str(fscareak) + '\n'
				fsclines.append( fscline )
				x+=1
		
			f.writelines( fsclines )
			f.close()
			print "\ndone with fsc %d/%d" %(i,nfinal)

	else:
		ngroups = nfinal/options.step
		excedent = nfinal%options.step
		count=0
		for thisn in xrange( 0, nfinal, options.step ):
			
			startindx = thisn
			stopindx = options.step*(count+1)
			if count == ngroups -1 :
				stopindx = nfinal
				
			avgo = averagerfunc( options.inputodd, options, stopindx )
			
			avgo_w = avgo.copy()
			if options.sym:
				avgo_w.process_inplace('xform.applysym',{'sym':options.sym})	
			
			if options.mask1:
				print "parsed mask1 is",options.mask1
				avgo_w.process_inplace(options.mask1[0], options.mask1[1])
			
			if options.mask2:
				print "parsed mask2 is",options.mask2
				avgo_w.process_inplace(options.mask2[0], options.mask2[1])
		
			if options.maskfile:
				mask = EMData(options.maskfile,0)
				avgo_w.mult(mask)
			
			fscarea = calcfsc( options, preve, avgo_w, gold=False )

			fscareas.append( fscarea )
			fscareasdict.update({count:fscarea})
			
			f = open( options.path +'/n_vs_fsc.txt','w' )
			fsclines = []
			
			x=0
			sortedfscareasdict = collections.OrderedDict(sorted(fscareasdict.items()))
			
			for k in sortedfscareasdict:
				fscareak = sortedfscareasdict[k]
				fscline = str(x) + '\t'+ str(fscareak) + '\n'
				fsclines.append( fscline )
				x+=1
		
			f.writelines( fsclines )
			f.close()
			print "\ndone with fsc %d/%d" %(count,ngroups)

			count+=1

	return

	
if '__main__' == __name__:
	main()
	



	
