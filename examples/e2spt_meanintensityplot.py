#!/usr/bin/env python

import os
import matplotlib.pyplot as plt
from EMAN2 import *
from sys import argv
from optparse import OptionParser
import sys
import time
import numpy
import pylab
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

def main():

	progname = os.path.basename(sys.argv[0])
	usage = """Produces mean intensity histograms of stack of sub-volumes"""
			
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	
	parser.add_argument("--input",type=str,default=None,help="Stack you want to analyze.")
	parser.add_argument("--output",type=str,default=None,help="Stack to save the edited volumes (masked, lowpassed, etc).")
	
	parser.add_argument("--shrink", type=int,default=1,help="Optionally shrink the input volumes by an integer amount.")
	#parser.add_argument("--bins", type=int,default=5,help="Number of bins for histogram.")
	
	parser.add_argument("--sym", type=str, default='c1', help = "Symmetry to enforce before computing mean intensity in the box. Note that this should only be used if the particles are properly aligned to the symmetry axis.")
	parser.add_argument("--mask",type=str,help="Mask processor applied to particles before alignment. Default is mask.sharp:outer_radius=-2", default="mask.sharp:outer_radius=-2")
	parser.add_argument("--maskfile",type=str,default=None,help="A maskfile that will be multiplied by the image.")
	
	#parser.add_argument("--normproc",type=str,help="Normalization processor applied to particles before alignment. Default is to use normalize.mask. If normalize.mask is used, results of the mask option will be passed in automatically. If you want to turn this option off specify \'None\'", default="normalize.mask")
	
	parser.add_argument("--preprocess",type=str,help="Any processor (as in e2proc3d.py) to be applied to each volume prior to alignment. Not applied to aligned particles before averaging.", default=None)
	parser.add_argument("--lowpass",type=str,help="A lowpass filtering processor (as in e2proc3d.py) to be applied to each volume prior to alignment. Not applied to aligned particles before averaging.", default=None)
	parser.add_argument("--highpass",type=str,help="A highpass filtering processor (as in e2proc3d.py) to be applied to each volume prior to alignment. Not applied to aligned particles before averaging.", default=None)

	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)

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
	
	cavity_scan(options)
	
	E2end(logger)

def cavity_scan(options):
	print "I am in the intensity scanner"
	intensities = []
	data = options.input
	
	n = EMUtil.get_image_count(data)
	
	for i in range(n):
		a=EMData(data,i)
		print "I am analyzing particle number", i

		# Make the mask first, use it to normalize (optionally), then apply it 
		mask=EMData(a["nx"],a["ny"],a["nz"])
		mask.to_one()
		
		if options.mask != None:
			#print "This is the mask I will apply: mask.process_inplace(%s,%s)" %(options.mask[0],options.mask[1]) 
			mask.process_inplace(options.mask[0],options.mask[1])
			print "Options mask 0, and 1 are", options.mask[0], options.mask[1]
		
		# normalize
		#if options["normproc"] != None:
		#	if options["normproc"][0]=="normalize.mask" : options["normproc"][1]["mask"]=mask
		#	fixedimage.process_inplace(options["normproc"][0],options["normproc"][1])
		#	image.process_inplace(options["normproc"][0],options["normproc"][1])
		
		a.mult(mask)
		
		if options.maskfile:
			m=EMData(options.maskfile,0)
			a.mult(m)		

		# preprocess
		if options.preprocess != None:
			a.process_inplace(options.preprocess[0],options.preprocess[1])
			
		# lowpass
		if options.lowpass != None:
			a.process_inplace(options.lowpass[0],options.lowpass[1])
			
		# highpass
		if options.highpass != None:
			a.process_inplace(options.highpass[0],options.highpass[1])
		
		# Shrink
		if options.shrink !=None and options.shrink>1 :
			a=a.process("math.meanshrink",{"n":options.shrink})

		intensities.append(a['mean_nonzero']*1000)
		
		if options.output:
			if options.output == 'same':
				a.write_image(data.replace('.','_EDITED.'),i)
			else:
				a.write_image(options.output,i)
	
	msktag=''
	if options.mask:	
		msktag = str(options.mask[1]).replace(':','_').replace(' ','').replace('}','').replace('{','').replace("'",'').replace(',','')
	
	if options.maskfile:
		msktag = options.maskfile
	
	msktag=msktag.replace('.hdf','')
	msktag=msktag.split('.')[0]

	print "\n\n\n\n\n\n\nMSK TAG is\n", msktag
	plot_name = data.replace('.hdf', '_MEANintensityPLOT_MSK' + msktag + '.png')

	print "The total number of particles is", len(intensities)

	#plt.hist(y, bins=5)

	std = numpy.std(intensities)
	mean = numpy.mean(intensities)
	statistics = 'mean='+str(mean) + ' std='+str(std) + '\n'
	
	statsfile = plot_name.replace('.png','_STATS.txt')
	f=open(statsfile,'w')
	f.writelines(statistics)
	f.close()

	intensitiesfile = plot_name.replace('.png','_INTENSITIES.txt')
	
	intensitieslines = []
	for i in intensities:	
		intensitieslines.append(str(i)+'\n')

	g=open(intensitiesfile,'w')
	g.writelines(statistics)
	g.close()

	print "The standard deviation is", std
	cuberoot = numpy.power(len(intensities),1.0/3.0)
	print "The cuberoot of n is", cuberoot
	width = (3.5*std)/cuberoot
	print "Therefore the bin width according to Scott's normal reference rule is", width
	calcbins = (max(intensities) - min(intensities)) / width
	print "And the number of bins should be", calcbins

	plt.hist(intensities, calcbins, facecolor='blue', alpha=0.75)	

	'''
	FITTING A CURVE
	binmids = []
	
	for i in range( len(bins) ):
		print bins[i]
		if i < len(bins)-1:
			#print i
			#print "i was smaller, because len bins is", len(bins)
			binmids.append( (bins[i] + bins[i+1])/2 )
	
	
	polycoeffs = numpy.polyfit(binmids, n, 3)
	yfit = numpy.polyval(polycoeffs, binmids)
	
	print "\n\nBinmids are\n", binmids
	print "\n\nThe length of binmids is\n", len(binmids)
	
	#hist, bins = numpy.histogram(y, options.bins)
	'''

	plt.title("Mean density for " + str(len(intensities)) + " particles")
	plt.ylabel("Number of particles")
	plt.xlabel("Density")
	#plt.axis([min(intensities), max(intensities)])	
	#plt.tick_params(axis='both',which='both',direction='out',length=1,width=2)
	
	#plt.plot(binmids, yfit, 'r--', linewidth=1)

	plt.grid(True)
	#plt.show()
	
	plt.savefig(plot_name)
	#plt.clf()
	
	return()

if __name__ == '__main__':
	main()
