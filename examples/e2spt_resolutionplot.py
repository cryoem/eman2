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
import scipy
import pylab

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """Aligns a 3d volume to another by executing e2classaverage3d.py and then calculates the FSC between them by calling e2proc3d.py . It returns both a number for the resolution based on the FSC0.5 
	criterion(on the screen) and a plot as an image in .png format."""
			
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	
	parser.add_argument("--vol1", type=str, help="File for volume 1. Note that volume 2 will be aligned TO volume 1; that is, volume 1 is 'the reference'. The format MUST be '.hdf' or '.mrc' ", default=None)
	parser.add_argument("--vol2", type=str, help="File for volume 2. Volume 2 to will be 'moved' to find its best alignment to volume 1. The format MUST be '.hdf' or '.mrc' ", default=None)
	parser.add_argument("--output", type=str, help="Name for the .txt file that will contain the FSC data. If not specified, a default name will be used, vol1_VS_vol2, where vol1 and vol2 are taken from --vol1 and --vol2 without the format"., default=None)

	parser.add_argument("--sym", type=str, default=None, help = "Asymmetric unit to limit the alignment search to. Note that this should only be on IF the reference (--vol1) is ALREADY aligned to the symmetry axis.")
	parser.add_argument("--mask",type=str,help="Mask processor applied to particles before alignment. Default is mask.sharp:outer_radius=-2", default="mask.sharp:outer_radius=-2")
	parser.add_argument("--normproc",type=str,help="Normalization processor applied to particles before alignment. Default is to use normalize.mask. If normalize.mask is used, results of the mask option will be passed in automatically. If you want to turn this option off specify \'None\'", default="normalize.mask")
	parser.add_argument("--preprocess",type=str,help="Any processor (as in e2proc3d.py) to be applied to each volume prior to alignment. Not applied to aligned particles before averaging.", default=None)
	parser.add_argument("--saveallalign",action="store_true", help="If set, will save the alignment parameters after each iteration",default=True)

	parser.add_argument("--lowpass",type=str,help="A lowpass filtering processor (as in e2proc3d.py) to be applied to each volume prior to alignment. Not applied to aligned particles before averaging.", default=None)
	parser.add_argument("--highpass",type=str,help="A highpass filtering processor (as in e2proc3d.py) to be applied to each volume prior to alignment. Not applied to aligned particles before averaging.", default=None)
	
	parser.add_argument("--npeakstorefine", type=int, help="The number of best 'coarse peaks' from 'coarse alignment' to refine in search for the best final alignment. Default=4.", default=4)
	parser.add_argument("--align",type=str,help="This is the aligner used to align particles to the previous class average. Default is rotate_translate_3d:search=10:delta=10:dphi=10, specify 'None' to disable", returnNone=True, default="rotate_translate_3d:search=10:delta=10:dphi=10")
	parser.add_argument("--aligncmp",type=str,help="The comparator used for the --align aligner. Default is the internal tomographic ccc. Do not specify unless you need to use another specific aligner.",default="ccc.tomo")
	parser.add_argument("--ralign",type=str,help="This is the second stage aligner used to refine the first alignment. Default is refine_3d:search=2:delta=3:range=12", default="refine_3d:search=2:delta=3:range=12", returnNone=True)
	parser.add_argument("--raligncmp",type=str,help="The comparator used by the second stage aligner. Default is the internal tomographic ccc",default="ccc.tomo")
	
	parser.add_argument("--postprocess",type=str,help="A processor to be applied to the volume after averaging the raw volumes, before subsequent iterations begin.",default=None)
		
	parser.add_argument("--shrink", type=int,default=1,help="Optionally shrink the input volumes by an integer amount for coarse alignment.")
	parser.add_argument("--shrinkrefine", type=int,default=1,help="Optionally shrink the input volumes by an integer amount for refine alignment.")
	
	parser.add_argument("--parallel",  help="Parallelism. See http://blake.bcm.edu/emanwiki/EMAN2/Parallel", default="thread:1")
	
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n",type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")

	(options, args) = parser.parse_args()

	vol1 = options.vol1
	vol2 = options.vol2
	
	vol1ALIname =  vol1.replace('.',+'_ALI.')
	
	alicmd = """e2classaverage3d.py --input=""" + vol1ALIname + """ --output=""" + + 
	""" --ref=""" + vol2 + """ --npeakstorefine=""" + options.npeakstorefine + """ --verbose=""" + options.verbose + """ --mask=""" + options.mask + 
	""" --preprocess=""" + options.preprocess + """ --align=""" + options.align + """--parallel=""" + options.parallel + """ --ralign=""" + options.ralign + 
	""" --aligncmp=""" + options.aligncmp + """ --raligncmp=""" + """ --shrink=""" + options.shrink + """ --shrinkrefine=""" + options.shrinkrefine + 
	""" --saveali""" + """ --normproc=""" + options.normproc + """ --sym=""" + options.sym + """ --breaksym"""
	
	fscoutputname = options.output
	if not options.output:
		fscoutputname = 'FSC_' + vol1.split('.')[0] + '_VS_' + vol2.split('.')[0] + '.txt'
	
	fsccmd = """e2proc3d.py """ + vol2 + """ """ + fscoutputname + """ --calcfsc=""" + vol1ALIname
	
	pltcmd = 
	
	options.system(alicmd + ' && ' + fsccmd)
	
	findir=os.listdir(os.getcwd())
	
	if fscoutputname not in findir:
		print "I am waiting for the fsc file to be generated"
		while fscoutputname not in findir:
			findir=os.listdir(os.getcwd())
			time.sleep=(10)
	else:
		print "I found the fsc file", fscoutputname
	
		f= open(fscoutputname,'r')
		lines = f.readlines()
		
		k=0
		x = []
		values = []
		for line in lines:
			x.append(float(k))
			values.append(float( line.split()[-1] ) )
			k += 1
		
		polycoeffs = scipy.polyfit(x, values, 3)
		yfit = scipy.polyval(polycoeffs, x)
		
		pylab.plot(x, values, 'k.')
		pylab.plot(x, yfit, 'r-')
				
		pylab.title(plot_name)
		pylab.ylabel('FSC value')
		pylab.xlabel('pixel number')
		
		pylab.savefig(plot_name)

		plot_name = fscoutputname.replace('.txt','_PLOT.png')
		
		#pylab.plot(x, values, linewidth=1, marker='o', linestyle='--', color='r')	
		
		#a = plt.gca()
		#a.set_xlim(1,mults[-1])
		#a.set_ylim(0,max(x))


if __name__ == "__main__":
	main()
