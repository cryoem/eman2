#!/usr/bin/env python
#
# Author: Jesus Galaz, 2011?2012? - Last change 16/Sep/2014
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

import os
from EMAN2 import *
from sys import argv
from optparse import OptionParser
import sys

import numpy as np


from scipy.stats import norm


def main():
	#import pylab
	#import matplotlib.mlab as mlab
	import matplotlib.pyplot as plt

	progname = os.path.basename(sys.argv[0])
	usage = """Produces mean intensity histograms of stack of sub-volumes"""
			
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	
	parser.add_argument("--input",type=str,default='',help="""Comma separated stacks of images
		whose mean intensity distribution you want to plot.""")

	parser.add_argument("--path",type=str,default='',help="Directory to store results in. The default is a numbered series of directories containing the prefix 'sptsim'; for example, sptsim_02 will be the directory by default if 'sptsim_01' already exists.")

	parser.add_argument("--output",type=str,default='',help="""Name of output plot inf comparing
		two populations or more.""")
	
	parser.add_argument("--shrink", type=int,default=1,help="Optionally shrink the input volumes by an integer amount.")
	parser.add_argument("--bins", type=int,default=0,help="""Number of bins for histogram.
		If not provided, the optimal bin number will be calculated.""")
	
	#parser.add_argument("--sym", type=str, default='c1', help = "Symmetry to enforce before computing mean intensity in the box. Note that this should only be used if the particles are properly aligned to the symmetry axis.")
	
	parser.add_argument("--mask",type=str,help="Mask processor applied to particles before alignment. Default is mask.sharp:outer_radius=-2", default="mask.sharp:outer_radius=-2")
	parser.add_argument("--maskfile",type=str,default='',help="A maskfile that will be multiplied by the image.")
	
	#parser.add_argument("--normproc",type=str,help="Normalization processor applied to particles before alignment. Default is to use normalize.mask. If normalize.mask is used, results of the mask option will be passed in automatically. If you want to turn this option off specify \'None\'", default="normalize.mask")
	
	parser.add_argument("--clip",type=int,default=0,help="""Boxsize to clip particles before 
		computing mean and std values for each image. Default=0, which means no resizing.""")	
	
	parser.add_argument("--preprocess",type=str,help="""Any processor (as in e2proc3d.py) 
		to be applied to each image before computing mean and std values.""", default=None)
	
	parser.add_argument("--lowpass",type=str,help="""A lowpass filtering processor (as in e2proc3d.py) 
		to be applied before computing mean and std values for each image""", default=None)
	
	parser.add_argument("--highpass",type=str,help="""A highpass filtering processor (as in e2proc3d.py) 
		to be applied before computing mean and std values for each image.""", default=None)

	parser.add_argument("--threshold",type=str,help="""A thresholding processor (as in from e2proc3d.py) 
		to be applied before computing mean and std values for each image.""", default=None)
		
	parser.add_argument("--normproc",type=str,help="""Normalization processor applied to particles 
		before computing mean and std values for each iamge. Default is 'normalize.edgemean'.
		If normalize.mask is used, results of the mask option will be passed in automatically. 
		If you want to turn this option off specify \'None\'""", default="normalize.edgemean")

	parser.add_argument("--savepreprocessed",action="store_true", help="""Will save image stacks 
		after preprocessing options (lowpass, highpass, preprocess, masking, etc) have been 
		applied.""", default=False)
		
	parser.add_argument("--normalizeplot",action="store_true",help="""This will normalize the
		intensity values of the distribution to be between 0 and 1""")
	parser.add_argument("--removesigma",type=int,default=0,help="""Provide a value for the
		number of standard deviations away from the mean to consider values to exclude.
		For example, if --removesigma=3, values further than 3 standard deviations away from the 
		mean will be excluded.""")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)

	parser.add_argument("--verbose", "-v", type=int, dest="verbose", action="store", metavar="n", default=0, help="verbose level [0-9], higner number means higher level of verboseness")

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
	
	if options.threshold: 
		options.threshold=parsemodopt(options.threshold)
		
	if options.normproc: 
		options.normproc=parsemodopt(options.normproc)
	
	datafiles = options.input.split(',')
	
	
	from e2spt_classaverage import sptmakepath
	options = sptmakepath( options, 'meanintensityplots')
	
	
	intensitiesSeveral = []
	means = []
	stds = []
	
	for datafile in datafiles:
	
		intensitiesSingle = calcintensities( options, datafile )
	
		intensitiesSeveral.append( [ datafile, list(intensitiesSingle) ] )
		
		if options.normalizeplot:
			intensitiesSingleNorm = normintensities( intensitiesSingle, 0, 0 )
		
		ret = plotintensities( intensitiesSingleNorm, options, datafile )
		mean = ret[0]
		std = ret[1]
		means.append(mean)
		stds.append(std)
		
	print "\nIntensities several len is", len( intensitiesSeveral )
	if len( intensitiesSeveral ) > 1:
		
		datafile1 = intensitiesSeveral[0][0]
		datafile2 = intensitiesSeveral[1][0]
		
		intensities1 = intensitiesSeveral[0][1]
		intensities2 = intensitiesSeveral[1][1]
		n1 = len( intensities1 )
		n2 = len( intensities2 )
		
		zscore = ( means[0]-means[1] )/ np.sqrt( (stds[0]*stds[0])/n1 + (stds[1]*stds[1])/n2 )
		
		g = open(options.path + '/MIboth_INFO.txt','w')
		zscoreline = 'zscore=' + str(zscore)+' for ' + datafile1 + ' vs ' + datafile2 + ' \n'
		lines=[ zscoreline ]
		g.writelines(lines)
		g.close()
		
		print "\nzzzzzzz\n%s" %( zscoreline )
		
		absmax = absmin = 0
		if options.normalizeplot:
		
			minses = []
			maxes = []
			for intenS in intensitiesSeveral:
				minS = float(min( intenS[1] ))
				maxS = float(max( intenS[1] ))
				
				minses.append( minS )
				maxes.append( maxS )
	
			absmin = min( minses )
			absmax = max( maxes ) - absmin
		
		
		for intensities in intensitiesSeveral:	
			print "Type and len of intensities is", type(intensities[1]), len(intensities[1])
						
			if options.normalizeplot:
				print "normalize on"	
				intensitiesNorm = normintensities( intensities[1], absmin, absmax )
				
			plotintensities( intensitiesNorm, options, datafile, 'no' )
	
		plt.savefig(options.path + '/MIbothPlot.png')
		plt.clf()
		
	E2end(logger)


def normintensities( intensitiesR, minval=0, maxval=0 ):
	print "normalize function"
	intensitiesNormalizedMin = []
	
	print "Intensities R type and length are", type(intensitiesR), len(intensitiesR)
	imin = min( intensitiesR )
	if minval:
		imin = minval
	
	print "\nMin is", imin
	
	for x in range(len( intensitiesR )):
		intenNormMin = ( float(intensitiesR[x]) - imin )
		intensitiesNormalizedMin.append( intenNormMin )
	
	print "After minnorm, intensitiesNormalizedMin[10] is", intensitiesNormalizedMin[10]
	print "Len of intensitiesNormalizedMin is", len(intensitiesNormalizedMin)
	imax = float(max( intensitiesNormalizedMin ))
	print "max is", imax
	
	if maxval:
		imax = maxval
	
	#print "intensitiesNoramlizedMin are", intensitiesNormalizedMin
	
	intensitiesNormalizedMax = []
	i=0
	for x in range(len(intensitiesNormalizedMin)):
		intenNormMax = float(intensitiesNormalizedMin[x]) / imax
		#print "for intenNorm %d value is %f" %(i,intenNormMax)
		intensitiesNormalizedMax.append( intenNormMax )
		i+=1
	
	#print "intensitiesNormalizedMax are", intensitiesNormalizedMax

	finalIntensities = list( intensitiesNormalizedMax )
	#print "After maxnorm, Intensities[10] is", finalIntensities[10]
	
	#print "normalized intensities are", finalIntensities
	return finalIntensities


def calcintensities( options, datafile ):

	print "\n(e2spt_meanintensityplot) (calcintensities)"
	intensities = []
	
	n = EMUtil.get_image_count( datafile )
	
	hdr = EMData( datafile, 0, True)
	dimensionality = 3
	mask=EMData(hdr["nx"],hdr["ny"],hdr["nz"])
	
	if int(hdr['nz']) <= 1:
		dimensionality = 2
		mask=EMData(hdr["nx"],hdr["ny"])
	
	mask.to_one()
		
	if options.mask:
		#print "This is the mask I will apply: mask.process_inplace(%s,%s)" %(options.mask[0],options.mask[1]) 
		mask.process_inplace(options.mask[0],options.mask[1])
		#print "Options mask 0, and 1 are", options.mask[0], options.mask[1]
	
	
	for i in range(n):
		a=EMData(datafile,i)
		print "\nAnalyzing particle number %d of %d for stack %s" % (i,n,datafile)

		# Make the mask first, use it to normalize (optionally), then apply it 
		# normalize
		#if options["normproc"] != None:
		#	if options["normproc"][0]=="normalize.mask" : options["normproc"][1]["mask"]=mask
		#	fixedimage.process_inplace(options["normproc"][0],options["normproc"][1])
		#	image.process_inplace(options["normproc"][0],options["normproc"][1])
		
		#print "\nImg size", a['nx']
		#print "Mask size",mask['nx']
		
		
		if options.clip:						
			sx = a['nx']
			sy = a['ny']
			sz = a['nz']
	
			xc = sx/2
			yc = sy/2
			zc = sz/2
	
			newsize = options.clip
	
			a['origin_x'] = 0
			a['origin_y'] = 0
			a['origin_z'] = 0
			
			print "\nBEFORE CLIP, intensity is", a['mean_nonzero']
	
			r=Region( (2*xc - newsize)/2, (2*yc - newsize)/2, (2*zc - newsize)/2, newsize , newsize , newsize)		
			
			a.clip_inplace( r )
			if i == 0:
				mask.clip_inplace( r )
		
		print "BEFORE normalizing, intensity is", a['mean_nonzero']
		#a.process_inplace('normalize.edgemean')
		
		if options.normproc != None:
			a.process_inplace(options.normproc[0],options.normproc[1])
			print "AFTER normalizing, intensity is", a['mean_nonzero']	

		# preprocess
		if options.preprocess != None:
			a.process_inplace(options.preprocess[0],options.preprocess[1])
		
		# lowpass
		if options.lowpass != None:
			a.process_inplace(options.lowpass[0],options.lowpass[1])
			print "After LOWPASS, intensity is", a['mean_nonzero']

		# highpass
		if options.highpass != None:
			a.process_inplace(options.highpass[0],options.highpass[1])
			print "After HIGHPASS, intensity is", a['mean_nonzero']

		if options.threshold != None:
			a.process_inplace(options.threshold[0],options.threshold[1])
			print "After THRESHOLD, intensity is", a['mean_nonzero']
		
		a.mult(mask)
		print "After MASKING, intensity is", a['mean_nonzero']
		
		if options.maskfile:
			m=EMData(options.maskfile,0)
			a.mult(m)	
		
		# Shrink
		if options.shrink !=None and options.shrink>1 :
			a.process_inplace("math.meanshrink",{"n":options.shrink})
			print "After SHRINKING, intensity is", a['mean_nonzero']
		
		intensities.append(a['mean_nonzero']*1000)
		
		if options.savepreprocessed:
			a.write_image(options.path + '/' + datafile.replace('.','_EDITED.'),i)

	finalvalues = []
	stddin = np.std( intensities )
	meanin = np.mean( intensities )
	
	for val in intensities:
		if not options.removesigma:
			finalvalues.append( val )
		elif options.removesigma:
			topthreshold = float( options.removesigma ) * stddin + meanin
			bottomthreshold = meanin - float( options.removesigma ) * stddin
			if float( val ) < topthreshold and float(val) > bottomthreshold:
				finalvalues.append( val )
			else:
				print """Value %f excluded because bottom and top thresholds to include are bottomthresh=%f, topthresh=%f
				""" %( float(val), bottomthreshold, topthreshold )	
	
	intensitiestxt = options.path + '/' + datafile.replace('.hdf','_INTENSITIES.txt')
	
	
	#if options.normalizeplot:
	#	intensitiesNormalized = []
	#	imax = max( intensities )
	#	imin = min( intensities )
	#	for inten in intensities:
	#		intenNorm = ( inten - imin ) / imax
	#		intensitiesNormalized.append( intenNorm )
	#	
	#	intensities = list( intensitiesNormalized )
	
	f=open( intensitiestxt, 'w')
	lines = []
	k=0
	for inten in finalvalues:
		lines.append( str(k) + ' ' + str(inten) + '\n')
		k+=1
		
	f.writelines(lines)
	f.close()
	
	plotintensities( finalvalues, options, datafile, 'yes' )
	
	return finalvalues
	
	
def plotintensities( intensities, options, datafile, onefile='yes' ):	
	import matplotlib.pyplot as plt
	
	#msktag=''
	#if options.mask:	
	#	msktag = str(options.mask[1]).replace(':','_').replace(' ','').replace('}','').replace('{','').replace("'",'').replace(',','')
	
	#if options.maskfile:
	#	msktag = options.maskfile
	
	#msktag=msktag.replace('.hdf','')
	#msktag=msktag.split('.')[0]

	#print "\n\n\n\n\n\n\nMSK TAG is\n", msktag
	plotname = datafile.replace('.hdf', '_MIplotMSK.png')

	print "The total number of particles is", len(intensities)

	#plt.hist(y, bins=5)

	std = np.std(intensities)
	mean = np.mean(intensities)
	statistics = ['mean='+str(mean) + ' std='+str(std) + '\n']

	#intensitiesfile = plotname.replace('.png','_INTENSITIES.txt')
	
	#intensitieslines = []
	#for i in intensities:	
	#	intensitieslines.append(str(i)+'\n')

	#g=open(options.path + '/' + intensitiesfile,'w')
	#g.writelines(statistics)
	#g.close()


	print "\nThe standard deviation is", std
	cuberoot = np.power(len(intensities),1.0/3.0)
	print "The cuberoot of n is", cuberoot
	width = (3.5*std)/cuberoot
	print "Therefore the bin width according to Scott's normal reference rule is", width
	
	calcbins = (max(intensities) - min(intensities)) / width
	
	if options.bins:
		calcbins = options.bins
	
	print "\nAnd the number of bins should be", calcbins
	
	statistics.append( 'bins=' + str( calcbins ) + ' , binwidth=' + str( width ) + '\n')
	
	statsfile = plotname.replace('.png','_INFO.txt')
	f=open(options.path + '/' + statsfile,'w')
	f.writelines(statistics)
	f.close()
	
	
	#calcbins=50
	plt.hist(intensities, calcbins, alpha=0.30, label=datafile)	
	
	intensities.sort()
	#hmean = np.mean(h)
	#hstd = np.std(h)
	
	pdf = norm.pdf(intensities, mean, std)
	plt.plot(intensities, pdf)

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

	plt.title("Mean density distribution")
	plt.ylabel("Number of particles")
	plt.xlabel("Density")
	
	#plt.axis([min(intensities), max(intensities)])	
	#plt.tick_params(axis='both',which='both',direction='out',length=1,width=2)
	
	#plt.plot(binmids, yfit, 'r--', linewidth=1)

	#plt.grid(True)
	#plt.show()
	
	if onefile == 'yes':
		plt.savefig(options.path + '/' +plotname)
		plt.clf()
	
	return [mean,std]

if __name__ == '__main__':
	main()
