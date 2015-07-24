#!/usr/bin/env python

#
# Author: Jesus Galaz, 04/01/2012 - Last update July/23/2015
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


from EMAN2 import *
from operator import itemgetter	

import sys
import numpy as np


def main():
	print "I have entered main"
	#progname = os.path.basename(sys.argv[0])
	#usage = """Aligns a 3d volume to another by executing e2spt_classaverage.py and then calculates the FSC between them by calling e2proc3d.py . It returns both a number for the resolution based on the FSC0.5 
	#criterion(on the screen) and a plot as an image in .png format."""
	
	progname = os.path.basename(sys.argv[0])
	usage = """Aligns a 3d volume to another by executing e2spt_classaverage.py and then calculates the FSC between them by calling e2proc3d.py . It returns both a number for the resolution based on the FSC0.5 
	criterion(on the screen) and a plot as an image in .png format."""
			
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
		
	parser.add_argument("--alistack", type=str, default='', help="""Default=None. Stack of particles aligned with e2spt_classaverage.p.y""")
	
	parser.add_argument("--scores",type=str,default='',help="""Default=None. Text file containing a list of correlations scores.""")
	
	parser.add_argument("--nbins", type=int,default=0,help="""Used for histogram plot. Default=0 (not used). Number of bins for histogram. If not provided, the optimal bin number will be automatically calculated based on bin-width, computed using Scott's normal reference rule, width = (3.5*std)/cuberoot(n), where 'std' is the standard deviation of the distribution of scores and n is the number of values considered. Then, bins will be nbins = (max(scores) - min(scores)) / width.""")
	
	parser.add_argument("--cutoff", type=float, help="""Fraction of particles (as a decimal, where 1.0 is the entire set, 0.8 is 80 percent. 0.5 is 50 percent, etc); 
														where to make the cutoff to divide the set into two groups. For example, if you specify 
														--cutoff=0.2, the 20 percent of particles with the highest correlation scores will be bundled into the first group, 
														and the remaining 80 percent into the second group.""", default=None)
	
	parser.add_argument("--lowpass",type=str,default='default',help="""Filter applied to averages when --groups is used. Default=filter.lowpass.gauss:cutoff_freq=1/F, where F is Nyquist frequency, automatically calculated. To disable type --lowpass=None""")
	
	parser.add_argument("--groups", type=int, help="Number of groups you want the data to be divided into based on correlation.", default=None)
	
	parser.add_argument("--topn", type=int, help="Number of particles from best to worst that you want to be written as a substack, averaged, and generate a coordinates .txt file with their coordinates.", default=None)

	parser.add_argument("--sigmaprune", type=float, default=0.0, help = """Number of standard deviations below the mean to cut off particles; 
																		that is, the mean cross correlation coefficient of all particles will be computed, and those that are 
																		--sigmaprune=N standard deviations below the mean will not be considered.""")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n",type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	
	parser.add_argument("--path", type=str, help="Results directory. If not specified, defaults to e2sptcoeff/", default='sptcoeff')
	
	parser.add_argument("--normalizeplot",action="store_true", help="Make maximum correlation value on plot equal to 1 and the minimum equal to 0, and scale all other values accordingly.", default=False)

	
	(options, args) = parser.parse_args()
	
	if options.groups and (options.cutoff or options.topn):
		print "ERROR: you cannot specify --cutoff, --groups and --topn all at the same time. Choose one."
		sys.exit()

	if options.cutoff and (options.groups or options.topn):
		print "ERROR: you cannot specify --cutoff, --groups and --topn all at the same time. Choose one."
		sys.exit()
		
	if options.topn and (options.cutoff or options.groups):
		print "ERROR: you cannot specify --cutoff, --groups and --topn all at the same time. Choose one."
		sys.exit()
		
	logger = E2init(sys.argv, options.ppid)
	
	from e2spt_classaverage import sptmakepath
	options = sptmakepath( options, 'sptcoeff')
	
	print "\nI have read the parameters"
	
	scores=[]
	dataset={}
	#x=[]
	
	k = 0
	usenewstack = False
	newstack = options.alistack.replace('.hdf','_prepruned.hdf')
	
	if options.alistack:
		n = EMUtil.get_image_count(options.alistack)
		
		for i in range(n):
			ahdr = EMData(options.alistack,i,True)
			#print "ahdr is", ahdr	
			score = None

			if options.verbose:
				print  "\nanalyzing header for particle", i
			
			if 'spt_score' in ahdr.get_attr_dict():
				score=-1*ahdr['spt_score']
				print "spt_score is", score
			elif 'spt_coefficient' in ahdr.get_attr_dict():
				score=ahdr['spt_coefficient']*-1
				print "spt_coefficient is", score
			else:
				print "\nERROR: no score found in header for particle %d. Skipping it! A prepruned stack will be made with all the particles that did have score info in their header" %(i)
				a = EMData( options.alistack, i )
				a.write_image( newstack, k )
				usenewstack = True
			
			if score:
				scores.append( score )
		
	elif options.scores and not options.alistack:
		f = open( options.scores,'r')
		lines=f.readlines()
		f.close()
		scores = [ round(float(line.replace('\n','')),6) for line in lines]
	
	elif not options.scores and not options.alistack:
		print "\n(e2spt_coeffplot)(main) ERROR: you need to supply either --alistack or --scores, with the former taking precedence over the latter."
		sys.exit()
	
	
	if usenewstack:
		options.alistack = newstack
	
	n = len(scores)
	
	nstack = EMUtil.get_image_count( options.alistack )
	
	if n != nstack:
		print "\n!!!! WARNING: the number of scores %d does not match the number of images %d in the stack %s" %( n, nstack, options.alistack )

	if scores and n > 1:
		print "\nThe set has these many particles", len( scores )
	else:
		print "\nERROR: There seems to be no information on particle scores. The number of scores must be larger than one to do any statistics with them."
		sys.exit()
		
	for i in range( n ):
		#dataset.append({'score':float(score), 'particle':a})
		dataset.update({i:scores[i]})
		#x.append(i)

	dataset = sorted( dataset.items(), key=itemgetter(1), reverse=True )
	
	if options.normalizeplot:
		#for s in scores:
			#val = values[ele]
		minv1 = min(scores)
		#maxv1 = max(scores)
		#print "Min max before normalization was", minv,maxv1
		for k in range(len(scores)):
			scores[k] = scores[k] - minv1
		
		#minv2 = min(scores)
		maxv2 = max(scores)
		#print "After subtracting min, the are", minv2,maxv
		#print "Max before normalization was", maxv
		for k in range(len(scores)):
			scores[k] = scores[k] / maxv2
	
	
	scores.sort()
	scores.reverse()
	#scores=scores[:-2]

	'''
	c:plot the distribution of scores
	'''
	import matplotlib.pyplot as plt
	import pylab
	import matplotlib
	
	std = np.std( scores )
	mean = np.mean( scores )
	statistics = ['mean='+str(mean) + ' std='+str(std) + '\n']

	print "\nthe standard deviation %.6f, mean %.6f" %( std, mean )
	
	if not std:
		print "\nERROR: std=0, which means all intensity values are the same."
		sys.exit()
	
	cuberoot = np.power(len( scores ),1.0/3.0)
	width = (3.5*std)/cuberoot
	print "\naccording to Scott's normal reference rule, width = (3.5*std)/cuberoot(n), the width of the histogram bins will be", width
	
	calcbins = ( max(scores) - min( scores )) / width
	
	if options.nbins:
		calcbins = options.nbins
		print "\overwriting number of bins to be", options.nbins
	
	print "\nand the number of bins n = ( max(scores) - min(scores) ) / width will thus be", calcbins
	calcbins = round(calcbins)
	print "rounding to", calcbins
	
	statistics.append( 'bins=' + str( calcbins ) + ' , binwidth=' + str( width ) + '\n')
		
	if not calcbins:
		print "WARNING: nbins=0, which means max and min intensity are the same, which probably means all scores are zero. Defaulting nbins to number of partilces."
		calcbins = len( scores )
	
	datafile = ''
	if options.alistack:
		datafile = os.path.basename( options.alistack )
	elif options.scores:
		datafile = os.path.basename( options.scores )
		
	f=open(options.path + '/stats.txt','w')
	f.writelines(statistics)
	f.close()
	
	plt.hist( scores, calcbins, alpha=0.30, label=datafile)	

	plottitle = os.path.basename( options.alistack ).replace('.hdf','') + ' CC scores distribution histogram'
	plt.title( plottitle )
	plt.ylabel("Number of particles")
	plt.xlabel("Score (au)")
  
  	matplotlib.rc('xtick', labelsize=16) 
	matplotlib.rc('ytick', labelsize=16) 
  		 	
  	font = {'weight':'bold','size':16}
	matplotlib.rc('font', **font)
  		 	
	pylab.rc("axes", linewidth=2.0)
		
	pylab.xlabel('CC score', fontsize=16, fontweight='bold')
  	pylab.ylabel('Number of particles', fontsize=16, fontweight='bold')
  	
  	plt.savefig( options.path + '/scores_histogram.png' )
  	plt.clf()
  	
  	
  	'''
  	c:plot decay in ranked correlation scores
  	'''
  	x = [i for i in range(len(scores))]
  	
	plt.plot(x, scores, color='k', linewidth=3)
	plottitle = os.path.basename( options.alistack ).replace('.hdf','') + ' CC scores decay' 
	plt.title( plottitle )
	
	pylab.xlim([min(x),max(x)+0.1*max(x)])
	font = {'weight':'bold','size':16}
	matplotlib.rc('font', **font)

	pylab.xlabel('Particle index', fontsize=16, fontweight='bold')
  	pylab.ylabel('CC score (au)', fontsize=16, fontweight='bold')

	plt.savefig( options.path + '/scores_decay.png')
	plt.clf()
	
	
	'''
  	c:plot the distance to mean value in # of standard deviations
  	'''
	distancestomean = [ (scores[i]-mean)/std for i in range(len(scores))]
	plt.plot(x, distancestomean, color='b', linewidth=2,)
	
	plottitle = os.path.basename( options.alistack ).replace('.hdf','') + ' CC scores distance to mean'
	plt.title( plottitle )
	
	pylab.xlim([min(x),max(x)+0.1*max(x)])
	
	font = {'weight':'bold','size':16}
	matplotlib.rc('font', **font)
	
	pylab.xlabel('Particle index', fontsize=16, fontweight='bold')
  	pylab.ylabel('Distance from mean (sigmas)', fontsize=16, fontweight='bold')
	
	plt.savefig( options.path + '/scores_distance2mean.png')
	plt.clf()
	
	
	'''
  	c:plot the distance to max value in # of standard deviations
  	'''
	maxi = max(scores)
	distancestomax = [ (maxi-scores[i])/std for i in range(len(scores))]
	plt.plot(x, distancestomax, color='b', linewidth=2,)
	
	plottitle = os.path.basename( options.alistack ).replace('.hdf','') + ' CC scores distance to max'
	plt.title( plottitle )
	
	pylab.xlim([min(x),max(x)+0.1*max(x)])
	
	font = {'weight':'bold','size':16}
	matplotlib.rc('font', **font)
	
	pylab.xlabel('Particle index', fontsize=16, fontweight='bold')
  	pylab.ylabel('Distance from max (sigmas)', fontsize=16, fontweight='bold')
	
	plt.savefig( options.path + '/scores_distance2max.png')
	plt.clf()
	
	
	print 'sorted dataset is', dataset
	
	'''
	c:prune and/or divide into groups
	'''		
	newscores=[]
	if options.sigmaprune:

		#print "Will analyze scores to remove aberrantly low ones"
	
		#mu=numpy.mean(scores)
		#sigma=numpy.std(scores)
			
		mu=mean
		sigma=std
		
		print "\nMean is", mu
		print "Std is", sigma
		filter = mu - sigma * (options.sigmaprune)
		print "Therefore, filter is", filter

		for d in dataset:
			if float(d['score']) < float(filter):
				dataset.remove(d)
				print "I have removed this aberrant particle from the dataset due to its low score", d['score']
			else:
				newscores.append(float(d['score']))
				print "This score is high enough to survive", d['score']
		
		newscores.sort()
		newscores.reverse()
		x=range(len(newscores))
	
		plottitle = os.path.basename(options.alistack).replace('.hdf', '_prunedSCORES')
		plt.plot(x, newscores,color='k', linewidth=2)
		pylab.xlim([min(x),max(x)+0.1*max(x)])
		#plt.title(plottitle)
		plt.ylabel('CC coefficient')
		plt.xlabel('Particle number')
		#a = plt.gca()
		plotfile = options.path + '/' + plottitle + '.png'
		plt.savefig(plotfile)
		plt.clf()
		
	newN=len(dataset)
	
	if options.lowpass:
		if options.lowpass == "default":
			hdr = EMData( options.alistack, 0, True )
			apix=hdr['apix_x']
			nyquist = 1.0/2.0*apix
			options.lowpass = 'filter.lowpass.tanh:cutoff_freq='+str(nyquist)+':apix='+str(apix)
			if apix =='1.0':
				print "\nWARNING: apix is 1.0, most likely wrong (default empty value). You can fix/change it with e2fixheaderparam.py"
		
		if options.lowpass and options.lowpass != 'None' and options.lowpass != 'none': 
			options.lowpass=parsemodopt(options.lowpass)
		elif 'None' in options.lowpass or 'none' in options.lowpass:
			options.lowpass=None
		
	if options.groups:
		
		if not options.alistack:
			print "\nERROR: --groups requires --alistack"
			sys.exit()
	
		halfptclnum= int(round(n/2.0))
		halfptcl = dataset[halfptclnum][0]
		print "half ptclnum is", halfptclnum
		print "which happens to be ptcl indx", halfptcl
		halfscore = dataset[halfptcl][1]
		print "with halfscore being", halfscore
		
		subdatasetsSET=[]
		N=len(dataset)
		#print "THe type of dataset is", type(dataset)
		print "The len of the dataset is", N
		subN = N/options.groups
		print "The len of each subset, except the last, should be", subN
		for g in range(options.groups):
			#subscores = scores[g*subN : (g+1)*subN]
			subdataset = dataset[g*subN : (g+1)*subN]			
			if g == options.groups -1:
				subdataset = dataset[g*subN :]
			subdatasetsSET.append(subdataset)
	
		kk=0
		for subdataset in subdatasetsSET:
			print "The len of subset %d is %d" % (kk, len(subdataset))
			jj = 0		
			groupname = options.path + '/' + os.path.basename(options.alistack).replace('.', '_'+ str(kk).zfill(len(str(options.groups))) + '.')
			particleLIST = []			
			lines = []
			
			avgr = Averagers.get('mean.tomo')
			for element in subdataset:
				#particle = element['particle']
				
				ptclnum = element[0]
				img = EMData( options.alistack, ptclnum )
				img.write_image(groupname,jj)
				jj+=1
				particleLIST.append( ptclnum )
				
				avgr.add_image( img )
				
			averageNAME = groupname.replace('.hdf','_AVG.hdf')		
			
			average = avgr.finish()
			
			average['origin_x'] = 0
			average['origin_y'] = 0
			average['origin_z'] = 0
			
			average.process_inplace('normalize.edgemean')
			if options.lowpass:
				print "(e2spt_coeffplot)(main) --lowpass provided:", options.lowpass
				average.process_inplace(options.lowpass[0],options.lowpass[1])
			
			average.write_image(averageNAME,0)
			#print "\nThe group average has been written to", averageNAME
			kk+=1	
			
	if options.cutoff:
	
		if not options.alistack:
			print "\nERROR: --cutoff requires --alistack"
			sys.exit()
	
		threshptclnum=int(round(newN/(1/options.cutoff)))
		threshptcl=dataset[threshptclnum]	
		print "The new threshptcl is", threshptcl
		threshscore = dataset[threshptclnum][1]

		threshlabel="%0.2f" %(options.cutoff)
		threshlabel=str(threshlabel).replace('.','p')	
		group1=[]
		group2=[]
		k1=0
		k2=0
		g1name = options.path + '/' + os.path.basename(options.alistack).replace('.', '_G1.')
		g2name = options.path + '/' + os.path.basename(options.alistack).replace('.', '_G2.')
		avgr1 = avgr = Averagers.get('mean.tomo')
		avgr2 = avgr = Averagers.get('mean.tomo')
		for i in dataset:
			img = EMData( options.alistack, i[0] )
			
			if i[1] >= threshscore:
				group1.append( i[0] )
				img.write_image(g1name,k1)
				avgr1.add_image( img )
				k1+=1
			else:
				group2.append( i[0] )
				img.write_image(g2name,k2)
				avgr2.add_image( img )
				k2+=1
				#else:
				#	print "\n\n@@@@ Found a garbage particle, and thus did not consider it!\n\n"	

		if group1:
			g1avg = avgr1.finish()
			g1avg.write_image(g1name.replace('.hdf','_avg.hdf'),0)

		if group2:
			g2avg = avgr1.finish()
			g2avg.write_image(g2name.replace('.hdf','_avg.hdf'),0)
	
	if options.topn:
	
		if not options.alistack:
			print "\nERROR: --topn requires --alistack"
			sys.exit()
			
		topndataset = dataset[0:options.topn]
		bottomdataset = dataset[options.topn:]
		
		outnamestack = options.path + '/' + os.path.basename(options.alistack).replace('.','top' + str(options.topn) + '.')
		
		coordsname = options.path + '/' + os.path.basename(outnamestack).split('.')[0] + '_coords.txt'
		
		indxsname = options.path + '/' + os.path.basename(coordsname).replace('coords','indxs')
		
		k=0
		linescoords=[]
		linesindxs=[]
		topptcls=[]
		
		avgr = Averagers.get('mean.tomo')
		for ptcl in topndataset:
			p = ptcl[0]
			img = EMData( options.alistack, p )
			img.write_image(outnamestack,k)
			avgr.add_imag(img)
			topptcls.append(p)
			linescoords.append(str(p['ptcl_source_coord'])+ '\n')
			linesindxs.append(str(p['source_n'])+ '\n')
			k+=1
		
		avg = avgr.finish()
		avgname = options.path + '/' + os.path.basename(outnamestack).replace('.','_avg.') 
		
		avg.write_image(avgname,0)
		f=open(coordsname,'w')
		f.writelines(linescoords)
		f.close()
	
		f=open(indxsname,'w')
		f.writelines(linesindxs)
		f.close()
	
	#x=range(len(scores))
	
	#plottitle = os.path.basename(options.alistack).replace('.hdf', '_prunedSCORES')
	#plt.plot(x, scores, marker='o', color='r', linewidth=2)
	#plt.title(plottitle)
	#plt.ylabel('score')
	#plt.xlabel('n ptcl')
	#a = plt.gca()
	
	#plotfile = options.path + '/' + plottitle + '.png'
	#plt.savefig(plotfile)
	#plt.clf()
	
	E2end(logger)

	return()

if __name__ == '__main__':
	main()
