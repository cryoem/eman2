#!/usr/bin/env python

#
# Author: Jesus Galaz, 04/28/2012 - using code and concepts drawn from Jesus Galaz's scripts
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

from EMAN2 import *
import os
import sys
import time
import numpy
import pylab
#from operator import itemgetter
from matplotlib.ticker import MaxNLocator
from pylab import figure, show	
import matplotlib.pyplot as plt				 

def main():
	
	progname = os.path.basename(sys.argv[0])
	usage = """Aligns a 3d volume to another by executing e2classaverage3d.py and then calculates the FSC between them by calling e2proc3d.py . It returns both a number for the resolution based on the FSC0.5 
	criterion(on the screen) and a plot as an image in .png format."""
			
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	
	parser.add_argument("--vol1", type=str, help="File for volume 1. Note that volume 1 will be aligned TO volume 2; that is, volume 1 will be 'moved' to find its best fit to volume 2. The format MUST be '.hdf' or '.mrc' ", default=None)
	parser.add_argument("--vol2", type=str, help="File for volume 2. Volume 2 to will be 'static' (the 'reference' to which volume 1 will be aligned to). The format MUST be '.hdf' or '.mrc' ", default=None)
	parser.add_argument("--output", type=str, help="Name for the .txt file that will contain the FSC data. If not specified, a default name will be used, vol1_VS_vol2, where vol1 and vol2 are taken from --vol1 and --vol2 without the format.", default=None)

	parser.add_argument("--sym", type=str, default='c1', help = """Will symmetrize the reference (--vol1) and limit alignment of other structure (--vol2) against the model to searching the asymmetric unit only. 
								Then, after alignment, vol2 will be symmetrized as well, and the FSC will be calculated. 
								Note that this will only work IF the reference (--vol1) is ALREADY aligned to the symmetry axis as defined by EMAN2.""")
	parser.add_argument("--mask",type=str,help="Mask processor applied to particles before alignment. Default is mask.sharp:outer_radius=-2", default="mask.sharp:outer_radius=-2")
	parser.add_argument("--normproc",type=str,help="Normalization processor applied to particles before alignment. Default is to use normalize.mask. If normalize.mask is used, results of the mask option will be passed in automatically. If you want to turn this option off specify \'None\'", default="normalize.mask")
	parser.add_argument("--preprocess",type=str,help="Any processor (as in e2proc3d.py) to be applied to each volume prior to alignment. Not applied to aligned particles before averaging.", default=None)
	parser.add_argument("--saveallalign",action="store_true", help="If set, will save the alignment parameters after each iteration",default=True)

	parser.add_argument("--mirror",action="store_true", help="If set, it will generate a mirrored version of --vol2, align that against --vol1 (the reference), and then compute FSCs. This will be done IN ADDITION to aligment and FSC of the unmirrored, original version of --vol2",default=False)

	parser.add_argument("--lowpass",type=str,help="A lowpass filtering processor (as in e2proc3d.py) to be applied to each volume prior to alignment. Not applied to aligned particles before averaging.", default=None)
	parser.add_argument("--highpass",type=str,help="A highpass filtering processor (as in e2proc3d.py) to be applied to each volume prior to alignment. Not applied to aligned particles before averaging.", default=None)
	
	parser.add_argument("--npeakstorefine", type=int, help="The number of best 'coarse peaks' from 'coarse alignment' to refine in search for the best final alignment. Default=4.", default=4)
	
	parser.add_argument("--align",type=str,help="This is the aligner used to align particles to the previous class average.", default=None)
	parser.add_argument("--aligncmp",type=str,help="The comparator used for the --align aligner. Default is the internal tomographic ccc. Do not specify unless you need to use another specific aligner.",default="ccc.tomo")
	
	parser.add_argument("--ralign",type=str,help="This is the second stage aligner used to refine the first alignment. Default is refine_3d:search=2:delta=3:range=12", default="refine_3d:search=2:delta=3:range=12")
	parser.add_argument("--raligncmp",type=str,help="The comparator used by the second stage aligner. Default is the internal tomographic ccc",default="ccc.tomo")
	
	parser.add_argument("--postprocess",type=str,help="A processor to be applied to the volume after averaging the raw volumes, before subsequent iterations begin.",default=None)
		
	parser.add_argument("--shrink", type=int,default=1,help="Optionally shrink the input volumes by an integer amount for coarse alignment.")
	parser.add_argument("--shrinkrefine", type=int,default=1,help="Optionally shrink the input volumes by an integer amount for refine alignment.")
	
	parser.add_argument("--parallel",  help="Parallelism. See http://blake.bcm.edu/emanwiki/EMAN2/Parallel", default="thread:1")
	
	parser.add_argument("--apix", type=float, help="Provide --apix for automatic FSC calculation (one of two methods) if you supply --curves and no volumes for --vol1 and --vol2", default=0)
	parser.add_argument("--boxsize", type=float, help="Provide --boxsize for automatic FSC calculation (one of two methods) if you supply --curves and no volumes for --vol1 and --vol2", default=0)

	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n",type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	#parser.add_argument("--plotonly",action="store_true", help="Assumes vol1 and vol2 are already aligned and fsc files will be provided through --curves; thus skips alignment and fsc curve generation", default=False)
	parser.add_argument("--fsconly",action="store_true", help="Assumes vol1 and vol2 are already aligned with respect to each other and thus skips alignment", default=False)
	parser.add_argument("--curves",help="FSC curves to plot in separate plots. Skips alignment and fsc curve generation. Provide .txt. files separated by commas --curves=file1.txt,file2.txt,file3.txt etc...", default=None)
	
	#print "Before parsing, options are", parser.args()
	
	(options, args) = parser.parse_args()
	
	print "options are", options
	
	logger = E2init(sys.argv, options.ppid)

	vol1 = options.vol1
	vol2 = options.vol2
	vol2final = ''
	vol2mirrorfinal = ''
	
	if vol1:
		vol1hdr = EMData(vol1,0,True)
		apix = vol1hdr['apix_x']
		boxsize = vol1hdr['nx']
		options.apix = apix
		options.boxsize = boxsize
	
	fscoutputname = options.output
	fscoutputnamemirror = ''

	if not options.output:
		fscoutputname = 'FSC_' + vol1.split('.')[0] + '_VS_' + vol2.split('.')[0] + '.txt'
	
	if '.hdf' in fscoutputname:
		#output = options.output
		fscoutputname = fscoutputname.replace('.hdf','.txt')
		#options.output = output
		
	if '.txt' not in fscoutputname:
		fscoutputname = fscoutputname + '.txt'

	if options.mirror:
		vol2mirror=vol2.replace('.hdf','_mirror.hdf')
		os.system('e2proc3d.py ' + vol2 + ' ' + vol2mirror  + ' --process=xform.mirror:axis=x')
		fscoutputnamemirror = fscoutputname.replace('.txt','_mirror.txt')
		
	if options.sym is not 'c1' and options.sym is not "C1" and not options.align and not options.curves:
		vol1symname = vol1.replace('.hdf','_' + options.sym + '.hdf')
		os.system('e2proc3d.py ' + vol1 + ' ' + vol1symname + ' --sym=' + options.sym)
		
		vol1 = vol1symname
		
		vol2sym = vol2.replace('.','_' + options.sym + '.')
		fscoutputname = fscoutputname.replace('.','_' + options.sym + '.')
		os.system('e2proc3d.py ' + vol1 + ' ' + vol2sym + ' --sym=' + options.sym)
		vol2final = vol2sym
		
		vol2mirrorsym = vol2mirror.replace('.','_' + options.sym + '.')
		fscoutputnamemirror = fscoutputnamemirror.replace('.','_' + options.sym + '.')
		os.system('e2proc3d.py ' + vol1 + ' ' + vol2mirrorsym + ' --sym=' + options.sym)
		vol2mirrorfinal = vol2mirrorsy
	else:
		vol2final = vol2
		if options.mirror:			
			vol2mirrorfinal = vol2mirror

	alicmd = 'echo'
	fsccmd = 'echo'

	if not options.curves and not options.fsconly and options.align:
		print "Will compute all alignment, fsc and plotting"
		aligner=options.align
		aligner=aligner.split(':')
		print "The split aligner is", aligner
		newaligner=''
		
		if 'sym='+options.sym not in aligner:
			for element in aligner:
				if 'sym' in element:
					element='sym=' +options.sym
				newaligner=newaligner + element + ":"
			if newaligner[-1] == ':':
				newaligner = newaligner[:-1]
				options.align=newaligner
		
		print "The reformed aligner is", options.align
			
		vol2ALIname = vol2.replace('.hdf','_ALI.hdf')
		alicmd = 'e2classaverage3d.py --path=. --input=' + vol2 + ' --output=' + vol2ALIname + ' --ref=' + vol1 + ' --npeakstorefine=' + str(options.npeakstorefine) + ' --verbose=' + str(options.verbose) + ' --mask=' + options.mask + ' --preprocess=' + options.preprocess + ' --align=' + options.align + ' --parallel=' + options.parallel + ' --ralign=' + options.ralign + ' --aligncmp=' + options.aligncmp + ' --raligncmp=' + options.raligncmp + ' --shrink=' + str(options.shrink) + ' --shrinkrefine=' + str(options.shrinkrefine) + ' --saveali' + ' --normproc=' + options.normproc + ' --sym=' + options.sym + ' --breaksym'
		vol2final = vol2ALIname
		
		if options.mirror:
			vol2mirrorALIname = vol2mirror.replace('.hdf','_ALI.hdf')
			alicmd += ' && e2classaverage3d.py --path=. --input=' + vol2mirror + ' --output=' + vol2mirrorALIname + ' --ref=' + vol1 + ' --npeakstorefine=' + str(options.npeakstorefine) + ' --verbose=' + str(options.verbose) + ' --mask=' + options.mask + ' --preprocess=' + options.preprocess + ' --align=' + options.align + ' --parallel=' + options.parallel + ' --ralign=' + options.ralign + ' --aligncmp=' + options.aligncmp + ' --raligncmp=' + options.raligncmp + ' --shrink=' + str(options.shrink) + ' --shrinkrefine=' + str(options.shrinkrefine) + ' --saveali' + ' --normproc=' + options.normproc + ' --sym=' + options.sym + ' --breaksym'
			vol2mirrorfinal = vol2mirrorALIname
		
		if options.sym is not 'c1' and options.sym is not "C1":
			print "\n@@@@@@@@@@@@@@@@@@@@@@@@@@@\nSym is NOT C1, therefore, the vol2 must be symmetrized after alignment\n"
			
			vol2ALIsym = vol2ALIname.replace('.','_' + options.sym + '.')
			alicmd = alicmd + ' && e2proc3d.py ' + vol2ALIname + ' ' + vol2ALIsym + ' --sym=' + options.sym
			vol2final = vol2ALIsym
			
			if options.mirror:
				vol2mirrorALIsym = vol2mirrorALIname.replace('.','_' + options.sym + '.')
				alicmd +=  ' && e2proc3d.py ' + vol2mirrorALIname + ' ' + vol2mirrorALIsym + ' --sym=' + options.sym
				vol2mirrorfinal = vol2mirrorALIsym
			
		print "\nFINAL alicmd INSIDE alignment is", alicmd		
	
	if not options.curves:
		
		print "\n\nWIll compute the FSC now!!"
		print "\nBecuase options.plotonly is NEGATIVE, see", options.curves
		print "Where vol2final is", vol2final
		
		fsccmd = 'e2proc3d.py ' + vol2final + ' ' + fscoutputname + ' --calcfsc=' + vol1
		if options.mirror:
			print "Where vol2mirrorfinal is", vol2mirrorfinal
			fsccmd += ' && e2proc3d.py ' + vol2mirrorfinal + ' ' + fscoutputnamemirror + ' --calcfsc=' + vol1
			
		cmd = alicmd + ' && ' + fsccmd 
		print "\nFinal command ali and fsc in FSC calculation is\n", cmd

		if options.verbose:
			print "\n\nAli command is\n", alicmd
			print "\n\nFSC command is\n", fsccmd	
			print "\n\nThus together they are\n", cmd
		
		os.system(cmd)
		
		findir=os.listdir(os.getcwd())
		print "\n\n\n\nFiles in dir are", findir
		
		if fscoutputname not in findir:
			print "\n\nERROR: Could not find the FSC file", fscoutputname
			sys.exit()
		elif fscoutputname:
			fscplotter(fscoutputname,options.apix,options.boxsize)
			
		if options.mirror and fscoutputnamemirror not in findir:
			print "\n\nERROR: Could not fine the FSC file for the mirrored version of things", fscoutputnamemirror
			sys.exit()
		elif options.mirror and fscoutputnamemirror:
			fscplotter(fscoutputnamemirror,options.apix,options.boxsize)
			
	elif options.curves:
		curves = options.curves
		print "I will compute the plots only!"
		curves = curves.replace('--curves=','')
		curves = curves.split(',')

		#flag=options.sameplot

		for curve in curves:
			fscplotter(curve,options.apix,options.boxsize)
			
	E2end(logger)
	return()
			
def fscplotter(fscoutputname,apix,boxsize):
	print "I found THIS FSC file and will thus plot it", fscoutputname
	
	f= open(fscoutputname,'r')
	lines = f.readlines()
	
	if not boxsize:
		boxsize = len(lines)
	
	k=0
	x = []

	values = []
	inversefreqs = []
	inversefreqslabels = []

	xticks = []
	nlines=len(lines)
	factorOfTicks = round(nlines/(round(nlines/10)))
	print "\n\n\nTHEE factor of ticks is!!!!\n", factorOfTicks
	for line in lines:
		x.append(float(k))
		values.append( float( line.split()[-1] ) )
		inverse =  float( line.split()[0] )
		#print "inverse is"
		inversefreqs.append(inverse)

		element = ''
		if inverse:
			element = '1/' + str(int(round( 1.0/inverse  )))
		else:
			element='0'
		#print "Therefore, turned into a number it is", element
		inversefreqslabels.append(element)
		if not k % factorOfTicks:

			print "I have appended this tick!", element
			print "Because k is", k
			print "And k mod factorOfTicks is", k % factorOfTicks
			xticks.append(element)
		k += 1

	plot_name = fscoutputname.replace('.txt','_PLOT.png')

	fullinfo = {}

	difs0p5 = []
	difs0p143 = []


	polycoeffs = numpy.polyfit(x, values, 30)
	yfit = numpy.polyval(polycoeffs, x)

	for i in range(len(yfit)):
		dif0p5 = abs(yfit[i] - 0.5)
		difs0p5.append(dif0p5)
		fullinfo.update({dif0p5:i})

		dif0p143 = abs(yfit[i] - 0.143)
		difs0p143.append(dif0p143)
		fullinfo.update({dif0p143:i})


	difs0p5min1=min(difs0p5)
	difs0p5.remove(difs0p5min1)
	difs0p5min2=min(difs0p5)

	difs0p143min1=min(difs0p143)
	difs0p143.remove(difs0p143min1)
	difs0p143min2=min(difs0p143)

	fsc0p5minpixel1 = fullinfo[difs0p5min1]
	fsc0p5minpixel2 = fullinfo[difs0p5min2]
	fsc0p5pixel = (fsc0p5minpixel1 + fsc0p5minpixel2)/2
	fsc0p5freq1 = inversefreqs[fsc0p5minpixel1]
	fsc0p5freq2 = inversefreqs[fsc0p5minpixel2]

	fsc0p143minpixel1 = fullinfo[difs0p143min1]
	fsc0p143minpixel2 = fullinfo[difs0p143min2]
	fsc0p143pixel = (fsc0p5minpixel1 + fsc0p143minpixel2)/2
	fsc0p143freq1 = inversefreqs[fsc0p143minpixel1]
	fsc0p143freq2 = inversefreqs[fsc0p143minpixel2]

	fsc0p5freqavg = (fsc0p5freq1 + fsc0p5freq2)/2.0

	fsc0p143freqavg = (fsc0p143freq1 + fsc0p143freq2)/2.0

	fsc0p5resolution1 = ''
	fsc0p5resolution1label=''

	fsc0p143resolution1 = ''
	fsc0p143resolution1label=''

	if fsc0p5pixel and apix and boxsize:		
		fsc0p5resolution1 = (float(apix) * float(boxsize)) / float(fsc0p5pixel)
		fsc0p5resolution1label = "%.1f" % ( fsc0p5resolution1 )
	else:
		print "Method 1 for resolution calculation failed (there was a division by zero somewhere, or you forgot to provide --boxsize or --apix)"


	if fsc0p143pixel and apix and boxsize:		
		fsc0p143resolution1 = (float(apix) * float(boxsize)) / float(fsc0p143pixel)
		fsc0p143resolution1label = "%.1f" % ( fsc0p143resolution1 )

	elif not fsc0p5pixel:
		print "Method 1 for resolution calculation failed (there was a division by zero somewhere)"

	fsc0p5resolution2=''
	fsc0p5resolution2label=''

	if fsc0p5freqavg:	
		fsc0p5resolution2 = 1/fsc0p5freqavg
		fsc0p5resolution2label = "%.1f" % ( fsc0p5resolution2 )
	else:
		print "Method 2 for resolution calculation failed (there was a division by zero somewhere)"


	fsc0p143resolution2=''
	fsc0p143resolution2label=''

	if fsc0p143freqavg:	
		fsc0p143resolution2 = 1/fsc0p143freqavg
		fsc0p143resolution2label = "%.1f" % ( fsc0p143resolution2 )

	elif not fsc0p5resolution2:
		print "Method 2 for resolution calculation failed (there was a division by zero somewhere)"


	if sum(values)/len(values) == 1.0:
		print "The particles you are aligning are exactly the same; cannot compute a reasonable FSC curve for a particle with itself! (but I'll plot the straight line anyway)."
	else:
		print "FSC0.5 resolution calculations 1 and 2 are", fsc0p5resolution1, fsc0p5resolution2
		print "FSC0.143 resolution calculations 1 and 2 are", fsc0p143resolution1, fsc0p143resolution2


	fig = figure()
	curve = pylab.plot(x, values, 'k',marker='o')

	#fit = pylab.plot(x, yfit, 'r-')

	ax = fig.add_subplot(111)
	plt.xticks(x,inversefreqslabels)
	print "Len of x and inversefreqslabels are", len(x), len(inversefreqslabels)

	print "len of ticks is", len(xticks)
	print "ticks are", xticks
	ax.xaxis.set_major_locator(MaxNLocator(nbins=len(xticks)))
	pylab.setp(ax, xticklabels=xticks)

	final0p5='NA'
	if fsc0p5resolution2label:
		final0p5 = fsc0p5resolution2label
	elif fsc0p5resolution1label:
		final0p5 = fsc0p5resolution1label

	final0p143='NA'
	if fsc0p143resolution2label:
		final0p143 = fsc0p143resolution2label
	elif fsc0p143resolution1label:
		final0p143 = fsc0p143resolution1label	

	pylab.title(plot_name)
	pylab.ylabel('FSC')
	pylab.xlabel('Frequency 1/Angstroms')
	pylab.grid(True)
	#pylab.legend( (curve, fit), ('Values', 'Fit'))
	pylab.ylim([0,1.2])

	ww=open(plot_name.replace('.png','_RESvalues.txt'),'w')
	resline1 = 'fsc0.5=' +fsc0p5resolution1label + ' , fsc0.143=' +fsc0p143resolution1label + '\n'
	resline2 = 'fsc0.5=' +fsc0p5resolution2label + ' , fsc0.143=' +fsc0p143resolution2label+ '\n'
	resline3 = 'fsc0.5=' +str(final0p5) + ' , fsc0.143=' + str(final0p143) + '\n'
	reslines = [resline1, resline2, resline3]
	ww.writelines(reslines)
	ww.close()

	#pylab.annotate("FSC0.5 = "+str(final0p5)+" A", xy=(0, 1), xytext=(300, -30), xycoords='data', textcoords='offset points',bbox=dict(boxstyle="round", fc="0.8"))
	#pylab.annotate("FSC0.143 = "+str(final0p143)+" A", xy=(0, 1), xytext=(300, -55), xycoords='data', textcoords='offset points',bbox=dict(boxstyle="round", fc="0.8"))
	#pylab.annotate("Sampling ="+"%.2f"%(apix) + " A/pixel", xy=(0, 1), xytext=(300, -80), xycoords='data', textcoords='offset points',bbox=dict(boxstyle="round", fc="0.8"))

	pylab.savefig(plot_name)
	
	return()
	
if __name__ == "__main__":
	main()
