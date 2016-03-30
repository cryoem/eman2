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
#from scipy.stats import norm

def main():
	#import pylab
	#import matplotlib.mlab as mlab
	import matplotlib.pyplot as plt

	progname = os.path.basename(sys.argv[0])
	usage = """Produces mean intensity histograms of stack of sub-volumes"""
			
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	
	parser.add_argument("--input",type=str,default='',help="""Default=None. Comma-separated stacks of images whose mean intensity distribution you want to plot.""")

	parser.add_argument("--subset",type=int,default=0,help="""Default=0 (not used). N > 2 number of particles to from each stack provided through --input to consider.""")

	parser.add_argument("--path",type=str,default='',help="Directory to store results in. The default is a numbered series of directories containing the prefix 'sptsim'; for example, sptsim_02 will be the directory by default if 'sptsim_01' already exists.")

	#parser.add_argument("--output",type=str,default='',help="""Name of output plot if comparing two populations or more.""")
	
	parser.add_argument("--shrink", type=int,default=1,help="Default=1 (no shrinking). Optionally shrink the input volumes by an integer amount n > 1.")
	
	parser.add_argument("--bins", type=int,default=0,help="""Default=0 (not used). Number of bins for histogram. If not provided, the optimal bin number will be automatically calculated based on bin-width, computed using Scott's normal reference rule, width = (3.5*std)/cuberoot(n), where 'std' is the standard deviation of the mean intensity distribution of population and n is the number of mean intensity values considered (this is affected by --removesigma). Then, bins will be nbins = (max(intensities) - min(intensities)) / width.""")
	
	#parser.add_argument("--sym", type=str, default='c1', help = "Symmetry to enforce before computing mean intensity in the box. Note that this should only be used if the particles are properly aligned to the symmetry axis.")
	
	parser.add_argument("--mask",type=str,default="mask.sharp:outer_radius=-2",help="Default=mask.sharp:outer_radius=-2. Mask processor applied to the particles before alignment. (See 'e2help.py processors' at the command line for a list of processors that can be applied through e2proc3d.py).")
	
	parser.add_argument("--maskfile",type=str,default='',help="""Default=None. An image file containing an additional mask to apply besides --mask.""")
		
	parser.add_argument("--clip",type=int,default=0,help="""Default=0 (not used). Boxsize to clip particles to before computing mean and standard deviation values for each image. (This can act as a mask, as you'd want to clip the boxes to a smaller size than their current, original size, excluding neighboring particles and background pixels/voxels).""")	
	
	parser.add_argument("--preprocess",type=str,default='',help="""Any processor to be applied to each image before computing mean and standard deviation values. (See 'e2help.py processors' at the command line for a list of processors that can be applied through e2proc3d.py).""")
	
	parser.add_argument("--lowpass",type=str,default='',help="""Default=None. A lowpass filtering processor to be applied before computing mean and standard deviation values for each image. (See 'e2help.py processors' at the command line for a list of processors that can be applied through e2proc3d.py).""")
	
	parser.add_argument("--highpass",type=str,default='',help="""Default=None. A highpass filtering processor to be applied before computing mean and standard deviation values for each image. (See 'e2help.py processors' at the command line for a list of processors that can be applied through e2proc3d.py).""")

	parser.add_argument("--threshold",type=str,default='',help="""A thresholding processor to be applied before computing mean and standard deviation values for each image. (See 'e2help.py processors' at the command line for a list of processors that can be applied through e2proc3d.py).""")
		
	parser.add_argument("--normproc",type=str,default="normalize.edgemean",help="""Default=normalize.edgemean. Normalization processor applied to particles before computing mean and standard deviation values for each iamge. If normalize.mask is used, --mask will be passed in automatically. If you want to turn normalization off specify \'None\'. (See 'e2help.py processors' at the command line for a list of processors that can be applied through e2proc3d.py).""")

	parser.add_argument("--savepreprocessed",action="store_true",default=False,help="""Default=False. If provided, this option will save the image stacks in --input after all preprocessing options (lowpass, highpass, preprocess, masking, etc.) have been applied.""")
		
	parser.add_argument("--normalizeplot",action="store_true",default=False,help="""Default=False. This will normalize the intensity values of the distribution to be between 0 and 1""")

	parser.add_argument("--removesigma",type=int,default=0,help="""Default=0. Provide a value for the number of standard deviations away from the mean to consider values to exclude. For example, if --removesigma=3, values further than 3 standard deviations away from the mean will be excluded.""")
	
	parser.add_argument("--ppid", type=int, help="Default=1. Set the PID of the parent process, used for cross platform PPID",default=-1)

	parser.add_argument("--verbose", "-v", type=int, default=0, help="Default 0. Verbose level [0-9], higner number means higher level of verboseness",dest="verbose", action="store", metavar="n")

	(options, args) = parser.parse_args()
	
	logger = E2init(sys.argv, options.ppid)
	
	'''
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
	'''
	
	from e2spt_classaverage import sptOptionsParser
	options = sptOptionsParser( options )
	
	datafiles = options.input.split(',')
	
	
	from e2spt_classaverage import sptmakepath
	options = sptmakepath( options, 'meanintensityplots')
	
	
	intensitiesSeveral = []
	iwzSeveral = []
	iminsSeveral = []
	imaxsSeveral = []
	istdsSeveral = []
	
	means = []
	stds = []
	
	from e2spt_classaverage import writeParameters
	cmdwp = writeParameters(options,'e2spt_meanintensityplot.py', 'sptmeanintensity')
	
	for datafile in datafiles:
		n = EMUtil.get_image_count(datafile)
		
		if options.subset:
			if options.subset < 3:
				print "ERROR:Subset must be > 2."
				sys.exit(1)
			
			n = options.subset
			
		if n < 3:
			print "ERROR: All stacks must have at least 3 particles in them. This one doesn't:", datafile
			sys.exit(1)
	
	for datafile in datafiles:
		ret = calcintensities( options, datafile )
		
		intensitiesSingle = ret[0]
		iwz = ret[1]
		imins = ret[2]
		imaxs = ret[3]
		istds = ret[4]
		
		intensitiesSeveral.append( [ datafile, list( intensitiesSingle ) ] )
		
		iwzSeveral.append( [ datafile, list( iwz ) ] )
		iminsSeveral.append( [ datafile, list( imins ) ] ) 
		imaxsSeveral.append( [ datafile, list( imaxs ) ] ) 
		istdsSeveral.append( [ datafile, list( istds ) ] ) 
				
		
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
		
		
		ret = plotintensities( iwz, options, datafile,'wz' )
		ret = plotintensities( imins, options, datafile,'mins' )
		ret = plotintensities( imaxs, options, datafile,'maxs' )
		ret = plotintensities( istds, options, datafile,'stds' )
		
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
			
			intensitiesNorm = intensities[1]			
			if options.normalizeplot:
				print "Normalizeplot on"	
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
	intensitiesWzeros = []
	intensitiesMins = []
	intensitiesMaxs = []
	intensitiesStds = []
	
	n = EMUtil.get_image_count( datafile )
	
	if options.subset and options.subset > 2:
		n = options.subset
		print "Taking subset", n
	
	hdr = EMData( datafile, 0, True)
	dimensionality = 3
	
	mask = EMData(hdr["nx"],hdr["ny"],hdr["nz"])
	
	if int(hdr['nz']) <= 1:
		dimensionality = 2
		
		mask=EMData(hdr["nx"],hdr["ny"])
	
	mask.to_one()
	
	print "dimensionality", dimensionality
	
	if options.clip:
		if dimensionality == 3:
			mask = clip3D( mask, options.clip)
		if dimensionality == 2:
			mask = clip2D( mask, options.clip)
	
	if options.mask:
		mask.process_inplace(options.mask[0],options.mask[1])
	
	print "Created mask of size", mask['nx'],mask['ny'],mask['nz']
	mask.write_image(options.path + '/mask.hdf',0)
	
	for i in range(n):
		a = EMData(datafile,i)
		
		print "\nAnalyzing particle number %d/%d for stack %s with mean %f and non-zero mean %f before preprocessing" % (i,n,datafile, a['mean'], a['mean_nonzero'])
		
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
				print "\nAFTER NORMALIZING with %s, non-zero mean intensity is %f" %( options.normproc, a['mean_nonzero'])	

		if options.preprocess:
			if options.verbose:
				print "\nBEFORE PREPROCESS, non-zero mean intensity is", a['mean_nonzero']
			a.process_inplace(options.preprocess[0],options.preprocess[1])
			if options.verbose:
					print "\nATER PREPROCESS, non-zero mean intensity is", a['mean_nonzero']
		
		if options.lowpass:
			if options.verbose:
				print "\nBEFORE LOWPASSING, non-zero mean intensity is", a['mean_nonzero']
			a.process_inplace(options.lowpass[0],options.lowpass[1])
			if options.verbose:
				print "\nAFTER LOWPASSING, non-zero mean intensity is", a['mean_nonzero']

		if options.highpass:
			if options.verbose:
				print "\nBEFORE HIGHPASSING, non-zero mean intensity is", a['mean_nonzero']
			a.process_inplace(options.highpass[0],options.highpass[1])
			if options.verbose:
				print "\nAFTER HIGHPASSING, non-zero mean intensity is", a['mean_nonzero']

		if options.threshold:
			if options.verbose:
				print "\nBEFORE THRESHOLDING, non-zero mean intensity is", a['mean_nonzero']
			a.process_inplace(options.threshold[0],options.threshold[1])
			if options.verbose:
				print "\nAFTER THRESHOLDING, non-zero mean intensity is", a['mean_nonzero']
		
		if options.mask:
			if options.verbose:
				print "\nBEFORE MASKING, non-zero mean intensity is", a['mean_nonzero']
			a.mult(mask)
			if options.verbose:
				print "\nAFTER MASKING, non-zero mean intensity is", a['mean_nonzero']
		
		if options.maskfile:
			m = EMData(options.maskfile,0)
			if options.clip:
				if dimensionality == 2:
					m = clip2D( m, options.clip )
				if dimensionality == 3:
					m = clip3D( m, options.clip )
	
			if options.verbose:	
				print "\nBEFORE MASKING with MASKFILE, non-zero mean intensity is", a['mean_nonzero']
			a.mult(m)
			if options.verbose:	
				print "\nAFTER MASKING with MASKFILE, non-zero mean intensity is", a['mean_nonzero']
			
		if options.shrink and options.shrink > 1 :
			if options.verbose:
				print "\nBEFORE SHRINKING, non-zero mean intensity is", a['mean_nonzero']
			a.process_inplace("math.meanshrink",{"n":options.shrink})
			if options.verbose:
				print "\nAFTER SHRINKING, non-zero mean intensity is", a['mean_nonzero']
		
		#intensities.append(a['mean_nonzero']*1000)
		#finalval=1+a['mean_nonzero']
		
		#intensities.append(finalval)
		
		#print "Value added to 1!!!!!",finalval
			
		intensities.append(a['mean_nonzero'])
		
		intensitiesWzeros.append(a['mean'])
		intensitiesMins.append(a['minimum'])
		intensitiesMaxs.append(a['maximum'])
		intensitiesStds.append(a['sigma'])
	
			
		if not a['mean_nonzero']:
			print "WARNING: mean intensity appended is zero"
		else:
			print "appended non-zero mean intensity of", a['mean_nonzero']
			
		if options.savepreprocessed:
			a.write_image(options.path + '/' + datafile.replace('.','_EDITED.'),i)

	

	
	
	
	finalvalues = []
	stdin = np.std( intensities )
	meanin = np.mean( intensities )
	finalvalues = prunevals( options, intensities, meanin, stdin )	
	
	intensitiestxt = options.path + '/' + datafile.replace('.hdf','_INTENSITIES.txt')		
	
	
	finalvaluesWz = []
	stdinWz = np.std( intensitiesWzeros )
	meaninWz = np.mean( intensitiesWzeros )
	finalvaluesWz = prunevals( options, intensitiesWzeros, meaninWz, stdinWz )
	
	intensitiestxtWz = options.path + '/' + datafile.replace('.hdf','_INTENSITIESwZeros.txt')	
	
	finalvaluesMins = []
	stdinMins = np.std( intensitiesMins  )
	meaninMins  = np.mean( intensitiesMins  )
	finalvaluesMins  = prunevals( options, intensitiesMins , meaninMins , stdinMins  )
	intensitiestxtMin = options.path + '/' + datafile.replace('.hdf','_MIN.txt')
	
	finalvaluesMaxs = []
	stdinMaxs = np.std( intensitiesMaxs )
	meaninMaxs = np.mean( intensitiesMaxs )
	finalvaluesMaxs = prunevals( options, intensitiesMaxs, meaninMaxs, stdinMaxs )
	intensitiestxtMax = options.path + '/' + datafile.replace('.hdf','_MAX.txt')
	
	finalvaluesStds = []
	stdinStds = np.std( intensitiesStds )
	meaninStds = np.mean( intensitiesStds )
	finalvaluesStds = prunevals( options, intensitiesStds, meaninStds, stdinStds )
	intensitiestxtStd = options.path + '/' + datafile.replace('.hdf','_STD.txt')
	
	
	#print "Type of intensities is", type(finalvalues)
	#stddinpruned = np.std( finalvalues )
	#meaninpruned = np.mean( finalvalues )
	
	#plotintensities( finalvalues, options, datafile, 'yes' )
	
	return [finalvalues, finalvaluesWz, finalvaluesMins, finalvaluesMaxs, finalvaluesStds]


def writetxt():

	f=open( intensitiestxt, 'w')
	lines = []
	k=0
	for inten in finalvalues:
		lines.append( str(k) + ' ' + str(inten) + '\n')
		k+=1
		
	f.writelines(lines)
	f.close()

	return


def prunevals( options, intensities, mean, std ):
	
	finalvalues = []
	for val in intensities:
		if not options.removesigma:
			finalvalues.append( val )
			print "finalval added", val
		elif options.removesigma:
			print "\nBecause --removesigma is non-zero", options.removesigma
			topthreshold = float( options.removesigma ) * stddin + meanin
			print "the top threshold for keeping particles is", topthreshold
			bottomthreshold = meanin - float( options.removesigma ) * stddin
			print "whereas the bottom threshold is", bottomthreshold			
			if float( val ) < topthreshold and float(val) > bottomthreshold:
				finalvalues.append( val )
				print """Value %f included""" %( float(val) )
				
			else:
				print """Value %f EXCLUDED""" %( float(val) )
	
	return finalvalues
	
	
def plotintensities( intensities, options, datafile, tag='', onefile='yes' ):	
	import matplotlib.pyplot as plt
	
	#msktag=''
	#if options.mask:	
	#	msktag = str(options.mask[1]).replace(':','_').replace(' ','').replace('}','').replace('{','').replace("'",'').replace(',','')
	
	#if options.maskfile:
	#	msktag = options.maskfile
	
	#msktag=msktag.replace('.hdf','')
	#msktag=msktag.split('.')[0]

	#print "\n\n\n\n\n\n\nMSK TAG is\n", msktag
	
	if tag:
		tag ='_' + tag
	plotname = datafile.replace('.hdf', '_MIplotMSK' + tag + '.png')

	print "The total number of particles is", len(intensities)
	#print "because intensities are", intensities

	#plt.hist(y, bins=5)
	
	if options.verbose==10:
		print "\n\n\nThe intensities to plot are", intensities
		print "\n\n\n\n\n"
	
	std = np.std(intensities)
	mean = np.mean(intensities)
	statistics = ['mean='+str(mean) + ' std='+str(std) + '\n']


	print "The number of particles kept is", len(intensities)
	print "The standard deviation of the mean intensity distribution for this population is", std
	
	if not std:
		print "ERROR: std=0, which means all intensity values are the same."
		sys.exit()
		
	
	
	cuberoot = np.power(len(intensities),1.0/3.0)
	#print "The cuberoot of n is", cuberoot
	width = (3.5*std)/cuberoot
	print "Therefore, according to Scott's normal reference rule, width = (3.5*std)/cuberoot(n), the width of the histogram bins will be", width
	
	calcbins = (max(intensities) - min(intensities)) / width
	
	if options.bins:
		calcbins = options.bins
	
	print "\nAnd the number of bins n = ( max(intensities) - min(intensities) ) / width will thus be", calcbins
	calcbins = round(calcbins)
	print "rounding to", calcbins
	
	statistics.append( 'bins=' + str( calcbins ) + ' , binwidth=' + str( width ) + '\n')
	
	print "statistics are", statistics
	
	if not calcbins:
		print "WARNING: nins=0, which means max and min intensity are the same, which probably means all intensities are zero. Defaulting nbins to number of partilces."
		calcbins = len(intensities)
			
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
