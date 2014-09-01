#!/usr/bin/env python

import os
from EMAN2 import *
from sys import argv
from optparse import OptionParser
import sys

import numpy as np
import pylab
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

from scipy.stats import norm


def main():

	progname = os.path.basename(sys.argv[0])
	usage = """Produces mean intensity histograms of stack of sub-volumes"""
			
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	
	parser.add_argument("--input",type=str,default='',help="""Comma separated stacks of images
		whose mean intensity distribution you want to plot.""")

	parser.add_argument("--path",type=str,default='',help="Directory to store results in. The default is a numbered series of directories containing the prefix 'sptsim'; for example, sptsim_02 will be the directory by default if 'sptsim_01' already exists.")

	parser.add_argument("--output",type=str,default='',help="""Name of output plot inf comparing
		two populations or more.""")
	
	parser.add_argument("--shrink", type=int,default=1,help="Optionally shrink the input volumes by an integer amount.")
	#parser.add_argument("--bins", type=int,default=5,help="Number of bins for histogram.")
	
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
	
		intensities = calcintensities( options, datafile )
	
		intensitiesSeveral.append( intensities )
	
		ret = plotintensities( intensities, options, datafile )
		mean = ret[0]
		std = ret[1]
		means.append(mean)
		stds.append(std)
		
	print "\nIntensities several len is", len( intensitiesSeveral )
	if len( intensitiesSeveral ) > 1:
		
		n1 = len( intensitiesSeveral[0] )
		n2 = len( intensitiesSeveral[1] )
		
		zscore = ( means[0]-means[1] )/ np.sqrt( (stds[0]*stds[0])/n1 + (stds[1]*stds[1])/n2 )
		
		g = open(options.path + '/zscore.txt','w')
		lines=['zscore=' + str(zscore)+'\n']
		g.writelines(lines)
		g.close()
		
		print "\nzzzzzzzzz\nZ-score is", zscore
		for inten in intensitiesSeveral:	
			plotintensities( inten, options, datafile, 'no' )
	
		plt.savefig(options.path + '/MIbothPlot.png')
		plt.clf()
		
	E2end(logger)


def calcintensities( options, datafile ):

	print "\nI am in the intensity scanner"
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
		print "\nI am analyzing particle number %d of %d for stack %s" % (i,n,datafile)

		# Make the mask first, use it to normalize (optionally), then apply it 
		# normalize
		#if options["normproc"] != None:
		#	if options["normproc"][0]=="normalize.mask" : options["normproc"][1]["mask"]=mask
		#	fixedimage.process_inplace(options["normproc"][0],options["normproc"][1])
		#	image.process_inplace(options["normproc"][0],options["normproc"][1])
		
		#print "\nImg size", a['nx']
		#print "Mask size",mask['nx']
		
		
		if options.clip:				
			#sys.exit()
		
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
	
	intensitiestxt = options.path + '/' + datafile.replace('.hdf','_Intensities.txt')
	
	f=open( intensitiestxt, 'w')
	lines = []
	k=0
	for inten in intensities:
		lines.append( str(k) + ' ' + str(inten) + '\n')
		k+=1
		
	f.writelines(lines)
	f.close()
	
	plotintensities( intensities, options, datafile, 'yes' )
	
	return intensities
	
	
def plotintensities( intensities, options, datafile, onefile='yes' ):	
	
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
	statistics = 'mean='+str(mean) + ' std='+str(std) + '\n'
	
	statsfile = plotname.replace('.png','_STATS.txt')
	f=open(options.path + '/' + statsfile,'w')
	f.writelines(statistics)
	f.close()

	intensitiesfile = plotname.replace('.png','_INTENSITIES.txt')
	
	intensitieslines = []
	for i in intensities:	
		intensitieslines.append(str(i)+'\n')

	g=open(options.path + '/' + intensitiesfile,'w')
	g.writelines(statistics)
	g.close()


	
	print "The standard deviation is", std
	cuberoot = np.power(len(intensities),1.0/3.0)
	print "The cuberoot of n is", cuberoot
	width = (3.5*std)/cuberoot
	print "Therefore the bin width according to Scott's normal reference rule is", width
	calcbins = (max(intensities) - min(intensities)) / width
	print "And the number of bins should be", calcbins
	
	
	#calcbins=50
	plt.hist(intensities, calcbins, alpha=0.30)	
	
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
