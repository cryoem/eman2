#!/usr/bin/env/ python

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
	
	parser.add_argument("--input",type=str,default=None,help="Image you want to apply the cylindrical mask to.")
	parser.add_argument("--output",type=str,default=None,help="File where you want to save the masked image.")
	
	parser.add_argument("--shrink", type=int,default=1,help="Optionally shrink the input volumes by an integer amount.")
	parser.add_argument("--bins", type=int,default=5,help="Number of bins for histogram.")
	
	parser.add_argument("--sym", type=str, default='c1', help = "Asymmetric unit to limit the alignment search to. Note that this should only be on IF the reference (--vol1) is ALREADY aligned to the symmetry axis.")
	parser.add_argument("--mask",type=str,help="Mask processor applied to particles before alignment. Default is mask.sharp:outer_radius=-2", default="mask.sharp:outer_radius=-2")
	
	#parser.add_argument("--normproc",type=str,help="Normalization processor applied to particles before alignment. Default is to use normalize.mask. If normalize.mask is used, results of the mask option will be passed in automatically. If you want to turn this option off specify \'None\'", default="normalize.mask")
	
	parser.add_argument("--preprocess",type=str,help="Any processor (as in e2proc3d.py) to be applied to each volume prior to alignment. Not applied to aligned particles before averaging.", default=None)
	parser.add_argument("--lowpass",type=str,help="A lowpass filtering processor (as in e2proc3d.py) to be applied to each volume prior to alignment. Not applied to aligned particles before averaging.", default=None)
	parser.add_argument("--highpass",type=str,help="A highpass filtering processor (as in e2proc3d.py) to be applied to each volume prior to alignment. Not applied to aligned particles before averaging.", default=None)
	parser.add_argument("--saveeditedvol",action="store_true",help="Saves the edited stack, mean that it is masked, normalized, lowpass filtered and symmetrized, and whatever other processing option was selected", default=False)
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)

	(options, args) = parser.parse_args()
	
	logger = E2init(sys.argv, options.ppid)
	
	#params_dict = parameters.__dict__
	#parameters = params_dict
	#cavity_scan(parameters)
	
	#if options.normproc: 
	#	options.normproc=parsemodopt(options.normproc)
	
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
	x = []
	y = []
	#k = 0
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
	
		#b=b.process('xform.scale',{'scale':1,'clip':clipf})
		
		#b.process_inplace('testimage.cylinder',{'radius':r,'height':parameters['height']})

		y.append(a['mean_nonzero'])
		
		if options.saveeditedvol:
			a.write_image(data.replace('.','_EDITED.'),i)
	
	msktag = str(options.mask[1]).split(':')[-1].replace(' ','').replace('}','').zfill(3)
	print "\n\n\n\n\n\n\nMSK TAG is\n", msktag
	plot_name = data.replace('.hdf', '_MEANintensityPLOT_MSK' + msktag + '.png')

	print "The total number of particles is", len(y)

	#plt.hist(y, bins=5)
	n, bins, patches = plt.hist(y, options.bins, normed=1, facecolor='blue', alpha=0.75)

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

	plt.title("Cavity densities for " + str(len(y)) + " particles")
	plt.ylabel("Number of particles")
	plt.xlabel("Density")
	
	#print "\n\nBins are\n", bins
	#print "Bin length is", len(bins)
	#print "N length is", len(n)
	
	print "\n\n n is\n", n
	y = [4]*11
	print "\n\ny is\n", y
	print "Y length is", len(y)
	
	#width=0.7*(bins[1]-bins[0])
	#center=(bins[:-1]+bins[1:])/2
	#plt.bar(center,hist,align='center',width=width)

	l = plt.plot(binmids, yfit, 'r--', linewidth=1)

	plt.grid(True)
	#plt.show()
	
	plt.savefig(plot_name)
	#plt.clf()
	
	return()

if __name__ == '__main__':
	main()
