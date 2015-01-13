#!/usr/bin/env python
#
# Author: Jesus Galaz, 2011?2012? - Last change 12/Jan/2015
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
	
	parser.add_argument("--input",type=str,default='',help="""Comma separated stacks of images whose mean intensity distribution you want to plot.""")

	parser.add_argument("--path",type=str,default='',help="Directory to store results in. The default is a numbered series of directories containing the prefix 'sptsim'; for example, sptsim_02 will be the directory by default if 'sptsim_01' already exists.")

	parser.add_argument("--output",type=str,default='',help="""Name of output plot if comparing two populations or more.""")
	
	parser.add_argument("--shrink", type=int,default=1,help="Optionally shrink the input volumes by an integer amount.")
	
	parser.add_argument("--bins", type=int,default=0,help="""Number of bins for histogram. If not provided, the optimal bin number will be calculated.""")
	
	#parser.add_argument("--sym", type=str, default='c1', help = "Symmetry to enforce before computing mean intensity in the box. Note that this should only be used if the particles are properly aligned to the symmetry axis.")
	
	parser.add_argument("--mask",type=str,default="mask.sharp:outer_radius=-2",help="Mask processor applied to particles before alignment. Default is mask.sharp:outer_radius=-2")
	
	parser.add_argument("--maskfile",type=str,default='',help="A maskfile that will be multiplied by the image.")
		
	parser.add_argument("--clip",type=int,default=0,help="""Boxsize to clip particles before computing mean and std values for each image. Default=0, which means no resizing.""")	
	
	parser.add_argument("--preprocess",type=str,default='',help="""Any processor (as in e2proc3d.py) to be applied to each image before computing mean and std values.""")
	
	parser.add_argument("--lowpass",type=str,default='',help="""A lowpass filtering processor (as in e2proc3d.py) to be applied before computing mean and std values for each image""")
	
	parser.add_argument("--highpass",type=str,default='',help="""A highpass filtering processor (as in e2proc3d.py) to be applied before computing mean and std values for each image.""")

	parser.add_argument("--threshold",type=str,default='',help="""A thresholding processor (as in from e2proc3d.py) to be applied before computing mean and std values for each image.""")
		
	parser.add_argument("--normproc",type=str,default="normalize.edgemean",help="""Normalization processor applied to particles before computing mean and std values for each iamge. Default is 'normalize.edgemean'. If normalize.mask is used, results of the mask option will be passed in automatically. If you want to turn this option off specify \'None\'""")

	parser.add_argument("--savepreprocessed",action="store_true",default=False,help="""Default=False. Will save image stacks after preprocessing options (lowpass, highpass, preprocess, masking, etc) have been applied.""")
		
	parser.add_argument("--normalizeplot",action="store_true",default=False,help="""Default=False. This will normalize the intensity values of the distribution to be between 0 and 1""")

	parser.add_argument("--removesigma",type=int,default=0,help="""Default=0. Provide a value for the number of standard deviations away from the mean to consider values to exclude. For example, if --removesigma=3, values further than 3 standard deviations away from the mean will be excluded.""")
	
	parser.add_argument("--ppid", type=int, help="Default=1. Set the PID of the parent process, used for cross platform PPID",default=-1)

	parser.add_argument("--verbose", "-v", type=int, default=0, help="Default 0. Verbose level [0-9], higner number means higher level of verboseness",dest="verbose", action="store", metavar="n")

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
		n = EMUtil.get_image_count(datafile)
		if n < 3:
			print "ERROR: All stacks must have at least 3 partilces. This one doesn't:", datafile
			sys.exit(1)
	
	for datafile in datafiles:
		intensitiesSingle = calcintensities( options, datafile )
	
		intensitiesSeveral.append( [ datafile, list(intensitiesSingle) ] )
		
		intensitiesSingleNorm = intensitiesSingle
		
		if options.normalizeplot:
			intensitiesSingleNorm = normintensities( intensitiesSingle, 0, 0 )
		
		#print "\]n\\n\nIntensities before plotting are", intensitiesSingleNorm
		#print "\n\n\n\n\n"
		
		ret = plotintensities( intensitiesSingleNorm, options, datafile )
		mean = ret[0]
		std = ret[1]
		means.append(mean)
		stds.append(std)
		
	#print "\nIntensities several len is", len( intensitiesSeveral )
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
				print "normalizeplot on"	
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


def clip3D( vol, size ):
	
	volxc = vol['nx']/2
	volyc = vol['ny']/2
	volzc = vol['nz']/2
	
	Rvol =  Region( (2*volxc - size)/2, (2*volyc - size)/2, (2*volzc - size)/2, size , size , size)
	vol.clip_inplace( Rvol )
	#vol.process_inplace('mask.sharp',{'outer_radius':-1})
	
	return vol


def clip2D( img, size ):
	
	imgxc = img['nx']/2
	imgyc = img['ny']/2
	#imgzc = img['nz']/2
	
	Rimg =  Region( (2*imgxc - size)/2, (2*imgyc - size)/2, 0, size , size , 1)
	img.clip_inplace( Rimg )
	#img.process_inplace('mask.sharp',{'outer_radius':-1})
	
	return img


def calcintensities( options, datafile ):

	print "\n(e2spt_meanintensityplot) (calcintensities)"
	intensities = []
	
	n = EMUtil.get_image_count( datafile )
	
	hdr = EMData( datafile, 0, True)
	dimensionality = 3
	
	mask = EMData(hdr["nx"],hdr["ny"],hdr["nz"])
	mask.to_one()
	
	if int(hdr['nz']) <= 1:
		dimensionality = 2
		mask=EMData(hdr["nx"],hdr["ny"])
	
	if options.clip:
		if dimensionality == 3:
			mask = clip3D( mask, options.clip)
		if dimensionality == 2:
			mask = clip2D( mask, options.clip)
	
	if options.mask:
		mask.process_inplace(options.mask[0],options.mask[1])
	
	for i in range(n):
		a = EMData(datafile,i)
		
		print "\nAnalyzing particle number %d/%d for stack %s" % (i,n,datafile)
		
		if options.clip:
			if options.verbose:
				print "\nBEFORE CLIPPING, mean intensity is", a['mean_nonzero']
			if dimensionality == 3:
				a = clip3D( a, options.clip)
			if dimensionality == 2:
				a = clip2D( a, options.clip)
			if options.verbose:
				print "\nAFTER CLIPPING with clipsize=%d, mean intensity is %f" %( options.clip, a['mean_nonzero'] )
		
		if options.normproc:
			if options.verbose:
				print "\nBEFORE NORMALIZING, mean intensity is", a['mean_nonzero']
			
			if options.normproc[0]=="normalize.mask":
				if options.mask:			
					options.normproc[1]["mask"]=mask
				else:
					print "ERROR: To use normalize.mask you also need to specify a mask"
					sys.exit(1)
			
			a.process_inplace(options.normproc[0],options.normproc[1])
			
			if options.verbose:
				print "\nAFTER NORMALIZING with %s, intensity is %f" %( options.normproc, a['mean_nonzero'])	

		if options.preprocess:
			if options.verbose:
				print "\nBEFORE PREPROCESS, intensity is", a['mean_nonzero']
			a.process_inplace(options.preprocess[0],options.preprocess[1])
			if options.verbose:
					print "\nATER PREPROCESS, intensity is", a['mean_nonzero']
		
		if options.lowpass:
			if options.verbose:
				print "\nBEFORE LOWPASSING, intensity is", a['mean_nonzero']
			a.process_inplace(options.lowpass[0],options.lowpass[1])
			if options.verbose:
				print "\nAFTER LOWPASSING, intensity is", a['mean_nonzero']

		if options.highpass:
			if options.verbose:
				print "\nBEFORE HIGHPASSING, intensity is", a['mean_nonzero']
			a.process_inplace(options.highpass[0],options.highpass[1])
			if options.verbose:
				print "\nAFTER HIGHPASSING, intensity is", a['mean_nonzero']

		if options.threshold:
			if options.verbose:
				print "\nBEFORE THRESHOLDING, intensity is", a['mean_nonzero']
			a.process_inplace(options.threshold[0],options.threshold[1])
			if options.verbose:
				print "\nAFTER THRESHOLDING, intensity is", a['mean_nonzero']
		
		if options.mask:
			if options.verbose:
				print "\nBEFORE MASKING, intensity is", a['mean_nonzero']
			a.mult(mask)
			if options.verbose:
				print "\nAFTER MASKING, intensity is", a['mean_nonzero']
		
		if options.maskfile:
			m = EMData(options.maskfile,0)
			if options.clip:
				if dimensionality == 2:
					m = clip2D( m, options.clip )
				if dimensionality == 3:
					m = clip3D( m, options.clip )
	
			if options.verbose:	
				print "\nBEFORE MASKING with MASKFILE, intensity is", a['mean_nonzero']
			a.mult(m)
			if options.verbose:	
				print "\nAFTER MASKING with MASKFILE, intensity is", a['mean_nonzero']
			
		if options.shrink and options.shrink > 1 :
			if options.verbose:
				print "\nBEFORE SHRINKING, intensity is", a['mean_nonzero']
			a.process_inplace("math.meanshrink",{"n":options.shrink})
			if options.verbose:
				print "\nAFTER SHRINKING, intensity is", a['mean_nonzero']
		
		#intensities.append(a['mean_nonzero']*1000)
		finalval=1+a['mean_nonzero']
		intensities.append(finalval)
		
		print "Value added to 1!!!!!",finalval
		
		if options.savepreprocessed:
			a.write_image(options.path + '/' + datafile.replace('.','_EDITED.'),i)

	finalvalues = []
	stddin = np.std( intensities )
	meanin = np.mean( intensities )
	
	for val in intensities:
		
		if not options.removesigma:
			finalvalues.append( val )
			print "finalval added", val
		elif options.removesigma:
			topthreshold = float( options.removesigma ) * stddin + meanin
			bottomthreshold = meanin - float( options.removesigma ) * stddin
			if float( val ) < topthreshold and float(val) > bottomthreshold:
				finalvalues.append( val )
			else:
				print """Value %f excluded because bottom and top thresholds to include are bottomthresh=%f, topthresh=%f
				""" %( float(val), bottomthreshold, topthreshold )	
	
	intensitiestxt = options.path + '/' + datafile.replace('.hdf','_INTENSITIES.txt')
	
	f=open( intensitiestxt, 'w')
	lines = []
	k=0
	for inten in finalvalues:
		lines.append( str(k) + ' ' + str(inten) + '\n')
		k+=1
		
	f.writelines(lines)
	f.close()
	
	#plotintensities( finalvalues, options, datafile, 'yes' )
	
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
	
	if options.verbose==10:
		print "\n\n\nThe intensities to plot are", intensities
		print "\n\n\n\n\n"
	
	std = np.std(intensities)
	mean = np.mean(intensities)
	statistics = ['mean='+str(mean) + ' std='+str(std) + '\n']

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
	
	#pdf = norm.pdf(intensities, mean, std)
	#plt.plot(intensities, pdf)

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
