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

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """Aligns a 3d volume to another by executing e2classaverage3d.py and then calculates the FSC between them by calling e2proc3d.py . It returns both a number for the resolution based on the FSC0.5 
	criterion(on the screen) and a plot as an image in .png format."""
			
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	
	parser.add_argument("--vol1", type=str, help="File for volume 1. Note that volume 1 will be aligned TO volume 2; that is, volume 1 will be 'moved' to find its best fit to volume 2. The format MUST be '.hdf' or '.mrc' ", default=None)
	parser.add_argument("--vol2", type=str, help="File for volume 2. Volume 2 to will be 'static' (the 'reference' to which volume 1 will be aligned to). The format MUST be '.hdf' or '.mrc' ", default=None)
	parser.add_argument("--output", type=str, help="Name for the .txt file that will contain the FSC data. If not specified, a default name will be used, vol1_VS_vol2, where vol1 and vol2 are taken from --vol1 and --vol2 without the format.", default=None)

	parser.add_argument("--sym", type=str, default='c1', help = "Asymmetric unit to limit the alignment search to. Note that this should only be on IF the reference (--vol1) is ALREADY aligned to the symmetry axis.")
	parser.add_argument("--mask",type=str,help="Mask processor applied to particles before alignment. Default is mask.sharp:outer_radius=-2", default="mask.sharp:outer_radius=-2")
	parser.add_argument("--normproc",type=str,help="Normalization processor applied to particles before alignment. Default is to use normalize.mask. If normalize.mask is used, results of the mask option will be passed in automatically. If you want to turn this option off specify \'None\'", default="normalize.mask")
	parser.add_argument("--preprocess",type=str,help="Any processor (as in e2proc3d.py) to be applied to each volume prior to alignment. Not applied to aligned particles before averaging.", default=None)
	parser.add_argument("--saveallalign",action="store_true", help="If set, will save the alignment parameters after each iteration",default=True)

	parser.add_argument("--lowpass",type=str,help="A lowpass filtering processor (as in e2proc3d.py) to be applied to each volume prior to alignment. Not applied to aligned particles before averaging.", default=None)
	parser.add_argument("--highpass",type=str,help="A highpass filtering processor (as in e2proc3d.py) to be applied to each volume prior to alignment. Not applied to aligned particles before averaging.", default=None)
	
	parser.add_argument("--npeakstorefine", type=int, help="The number of best 'coarse peaks' from 'coarse alignment' to refine in search for the best final alignment. Default=4.", default=4)
	parser.add_argument("--align",type=str,help="This is the aligner used to align particles to the previous class average. Default is rotate_translate_3d:search=10:delta=10:dphi=10, specify 'None' to disable", default="rotate_translate_3d:search=10:delta=10:dphi=10")
	parser.add_argument("--aligncmp",type=str,help="The comparator used for the --align aligner. Default is the internal tomographic ccc. Do not specify unless you need to use another specific aligner.",default="ccc.tomo")
	parser.add_argument("--ralign",type=str,help="This is the second stage aligner used to refine the first alignment. Default is refine_3d:search=2:delta=3:range=12", default="refine_3d:search=2:delta=3:range=12")
	parser.add_argument("--raligncmp",type=str,help="The comparator used by the second stage aligner. Default is the internal tomographic ccc",default="ccc.tomo")
	
	parser.add_argument("--postprocess",type=str,help="A processor to be applied to the volume after averaging the raw volumes, before subsequent iterations begin.",default=None)
		
	parser.add_argument("--shrink", type=int,default=1,help="Optionally shrink the input volumes by an integer amount for coarse alignment.")
	parser.add_argument("--shrinkrefine", type=int,default=1,help="Optionally shrink the input volumes by an integer amount for refine alignment.")
	
	parser.add_argument("--parallel",  help="Parallelism. See http://blake.bcm.edu/emanwiki/EMAN2/Parallel", default="thread:1")
	
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n",type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	
	#print "Before parsing, options are", parser.args()
	
	(options, args) = parser.parse_args()
	
	print "options are", options
	
	logger = E2init(sys.argv, options.ppid)

	vol1 = options.vol1
	vol2 = options.vol2
	
	vol1hdr = EMData(vol1,0,True)
	apix = vol1hdr['apix_x']
	boxsize = vol1hdr['nx']
	
	vol1ALIname =  vol1.replace('.','_ALI.')
	
	fscoutputname = options.output
	if not options.output:
		fscoutputname = 'FSC_' + vol1.split('.')[0] + '_VS_' + vol2.split('.')[0] + '.txt'
	
	alicmd = 'e2classaverage3d.py --path=. --input=' + vol1 + ' --output=' + vol1ALIname + ' --ref=' + vol1 + ' --npeakstorefine=' + str(options.npeakstorefine) + ' --verbose=' + str(options.verbose) + ' --mask=' + options.mask + ' --preprocess=' + options.preprocess + ' --align=' + options.align + ' --parallel=' + options.parallel + ' --ralign=' + options.ralign + ' --aligncmp=' + options.aligncmp + ' --raligncmp=' + options.raligncmp + ' --shrink=' + str(options.shrink) + ' --shrinkrefine=' + str(options.shrinkrefine) + ' --saveali' + ' --normproc=' + options.normproc + ' --sym=' + options.sym + ' --breaksym'
	
	print "FSCOUTPUTNAME is", fscoutputname
	print "Because options.output was", options.output
	
	fsccmd = 'e2proc3d.py ' + vol2 + ' ' + fscoutputname + ' --calcfsc=' + vol1ALIname
	
	print "Ali command is\n", alicmd
	print "\n\nFSC command is\n", fsccmd
	
	#cmd = alicmd + ' && ' + 'mv res/' + vol1ALIname + ' . && ' + fsccmd + ' && rm -r res' 
	
	cmd = alicmd + ' && ' + fsccmd + ' && rm -r res' 

		
	print "\n\nThus together they are\n", cmd
		
	os.system(cmd)
	
	findir=os.listdir(os.getcwd())
	
	if fscoutputname not in findir:
		print "I am waiting for the fsc file to be generated"
		while fscoutputname not in findir:
			print "Can't find it yet!"
			findir=os.listdir(os.getcwd())
			time.sleep=(10)
	else:
		print "I found the fsc file", fscoutputname
	
		f= open(fscoutputname,'r')
		lines = f.readlines()
		
		k=0
		x = []
		values = []
		inversefreqs = []
		for line in lines:
			x.append(float(k))
			values.append( float( line.split()[-1] ) )
			inversefreqs.append( float( line.split()[0] ))
			k += 1
		
		polycoeffs = numpy.polyfit(x, values, 10)
		yfit = numpy.polyval(polycoeffs, x)
		
		#print "\n\n\n\n\n Y FIT IS!!!\n", yfit
		
		plot_name = fscoutputname.replace('.txt','_PLOT.png')
		
		fullinfo = {}
		#fullinfo2 = {}

		difs = []
		
		for i in range(len(yfit)):
			dif = abs(yfit[i] - 0.5)
			difs.append(dif)
			fullinfo.update({dif:i})
			#fullinfo2.update({inversefreqs[i]:dif})
		
		difsmin1=min(difs)
		difs.remove(difsmin1)
		difsmin2=min(difs)
		
		
		minpixel1 = fullinfo[difsmin1]
		minpixel2 = fullinfo[difsmin2]
		
		fsc0p5pixel = (minpixel1+minpixel2)/2
		fscfreq1 = inversefreqs[minpixel1]
		fscfreq2 = inversefreqs[minpixel2]
		
		fscfreqavg = (fscfreq1 + fscfreq2)/2.0
		
		resolution1 = (float(apix) * float(boxsize)) / float(fsc0p5pixel)
		
		resolution2 = 1/fscfreq2
		
		print "The resolution calculations 1 and 2 are", resolution1, resolution2
		#difs = sorted(bestfinal, key=itemgetter('score'))
		
		pylab.plot(x, values, 'k.')
		pylab.plot(x, yfit, 'r-')
				
		pylab.title(plot_name)
		pylab.ylabel('FSC value')
		pylab.xlabel('pixel number')
		
		pylab.savefig(plot_name)
		
		#pylab.plot(x, values, linewidth=1, marker='o', linestyle='--', color='r')	
		
		#a = plt.gca()
		#a.set_xlim(1,mults[-1])
		#a.set_ylim(0,max(x))

	E2end(logger)
	
if __name__ == "__main__":
	main()
