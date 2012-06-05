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
	usage = """Plots the variation of correlation of a volume with itself as it is rotated in azimuth or altitude"""
			
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	
	parser.add_argument("--input", type=str, help="Filename of the .hdf volume whose rotational correlation plot you want to compute.", default=None)
	parser.add_argument("--output", type=str, help="Name for the .txt file with the results and the corresponding .png plot")
	
	parser.add_argument("--mask",type=str,help="Mask processor applied to particles before alignment. Default is mask.sharp:outer_radius=-2", default="mask.sharp:outer_radius=-2")
	parser.add_argument("--preprocess",type=str,help="Any processor (as in e2proc3d.py) to be applied to each volume prior to alignment. Not applied to aligned particles before averaging.", default=None)

	parser.add_argument("--lowpass",type=str,help="A lowpass filtering processor (as in e2proc3d.py) to be applied to each volume prior to alignment. Not applied to aligned particles before averaging.", default=None)
	parser.add_argument("--highpass",type=str,help="A highpass filtering processor (as in e2proc3d.py) to be applied to each volume prior to alignment. Not applied to aligned particles before averaging.", default=None)
		
	parser.add_argument("--shrink", type=int,default=1,help="Optionally shrink the input volumes by an integer amount for coarse alignment.")
	parser.add_argument("--shrinkrefine", type=int,default=1,help="Optionally shrink the input volumes by an integer amount for refine alignment.")
	
	parser.add_argument("--parallel",  help="Parallelism. See http://blake.bcm.edu/emanwiki/EMAN2/Parallel", default="thread:1")
	
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n",type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	parser.add_argument("--plotonly",action="store_true", help="Assumes vol1 and vol2 are already aligned with respect to each other and thus skips alignment", default=False)
	
	(options, args) = parser.parse_args()
		
	logger = E2init(sys.argv, options.ppid)

	if options.mask: 
		options.mask=parsemodopt(options.mask)
	
	if options.preprocess: 
		options.preprocess=parsemodopt(options.preprocess)
		
	if options.lowpass: 
		options.lowpass=parsemodopt(options.lowpass)
	
	if options.highpass: 
		options.highpass=parsemodopt(options.highpass)
		
	if options.postprocess: 
		options.postprocess=parsemodopt(options.postprocess)

	vol1 = EMData(ptions.input,0)
	
	# Make the mask first, use it to normalize (optionally), then apply it 
	mask = EMData(vol1["nx"],vol1["ny"],vol1["nz"])
	mask.to_one()
	if options["mask"]:
		#print "This is the mask I will apply: mask.process_inplace(%s,%s)" %(options["mask"][0],options["mask"][1]) 
		mask.process_inplace(options["mask"][0],options["mask"][1])

	# normalize
	if options["normproc"]:
		#if options["normproc"][0]=="normalize.mask": 
		#	options["normproc"][1]["mask"]=mask
		vol1.process_inplace('normalize')

	vol1.mult(mask)
	vol1.process_inplace('normalize')
	vol1.mult(mask)

	# preprocess
	if options["preprocess"] != None:
		vol1.process_inplace(options["preprocess"][0],options["preprocess"][1])

	# lowpass
	if options["lowpass"] != None:
		vol1.process_inplace(options["lowpass"][0],options["lowpass"][1])

	# highpass
	if options["highpass"] != None:
		vol1.process_inplace(options["highpass"][0],options["highpass"][1])

	# Shrinking both for initial alignment and reference
	if options["shrink"]!=None and options["shrink"]>1 :
		vol1=vol1.process("math.meanshrink",{"n":options["shrink"]})
	
	vol1hdr = EMData(vol1,0,True)
	apix = vol1hdr['apix_x']
	boxsize = vol1hdr['nx']
	
	az=0
	values = []
	azs=[]
	while az <= 360:
		vol2 = EMAN2.EMData()
		vol2 = vol1.copy()

		ccf = vol1.calc_ccf(vol2)
		ccf.process_inplace("xform.phaseorigin.tocorner") 
		ccf.process_inplace('normalize')

		#box = ccf.get_zsize()
		#r =  Region((box/2) - int(parameters['searchx']), (box/2) - int(parameters['searchy']), (box/2) - int(parameters['searchz']), 2 * int(parameters['searchx']) + 1, 2 * int(parameters['searchy']) + 1, 2 * int(parameters['searchz']) + 1) 
		#sub_ccf = ccf.get_clip(r)

		loc_sub = ccf.calc_max_location()
		
		#xbest = loc_sub[0]
		#ybest = loc_sub[1]
		#zbest = loc_sub[2]

		best_value = ccf.get_value_at(loc_sub[0],loc_sub[1],loc_sub[2])
		values.append(best_value)
		
		azs.append(az)
		
		az += options.delta
		vol2.rotate(az,0,0)
	
	fileoutputname = options.output
	
	f = open(fscoutputname,'w')
	f.writelines(values)
	f.close()
	
	#polycoeffs = numpy.polyfit(x, values, 30)
	#yfit = numpy.polyval(polycoeffs, x)

	plot_name = filecoutputname.replace('.txt','_PLOT.png')

	pylab.plot(azs, values, 'k-', marker='o')
	#fit = pylab.plot(x, yfit, 'r-')

	pylab.title(plot_name)
	pylab.ylabel('Correlation')
	pylab.xlabel('Azimuth')
	pylab.grid(True)
	#pylab.legend( (curve, fit), ('Values', 'Fit'))
	#pylab.yaxis([0,boxsize/2 + 10,0,1.2])

	#pylab.annotate("FSC0.5 = "+str(final0p5)+" A", xy=(0, 1), xytext=(300, -30), xycoords='data', textcoords='offset points',bbox=dict(boxstyle="round", fc="0.8"))
	#pylab.annotate("FSC0.143 = "+str(final0p143)+" A", xy=(0, 1), xytext=(300, -55), xycoords='data', textcoords='offset points',bbox=dict(boxstyle="round", fc="0.8"))
	#pylab.annotate("Sampling ="+"%.2f"%(apix) + " A/pixel", xy=(0, 1), xytext=(300, -80), xycoords='data', textcoords='offset points',bbox=dict(boxstyle="round", fc="0.8"))

	pylab.savefig(plot_name)

	#pylab.plot(x, values, linewidth=1, marker='o', linestyle='--', color='r')	
	#a = plt.gca()
	#a.set_xlim(1,mults[-1])
	#a.set_ylim(0,max(x))

	E2end(logger)
	
if __name__ == "__main__":
	main()
