#!/usr/bin/env python

#
# Author: Jesus Galaz, 04/01/2012 - Last update 05/01/2012
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
import matplotlib.pyplot as plt
import sys
import numpy

def main():
	print "I have entered main"
	progname = os.path.basename(sys.argv[0])
	usage = """Aligns a 3d volume to another by executing e2spt_classaverage.py and then calculates the FSC between them by calling e2proc3d.py . It returns both a number for the resolution based on the FSC0.5 
	criterion(on the screen) and a plot as an image in .png format."""
	
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--input", type=str,help="Stack of 3D volumes that have been aligned with EMAN2 SPT programs. The format MUST be .hdf, and the files must contain the header parameter spt_coefficient.", default=None)
	parser.add_argument("--cutoff", type=float,help="Fraction of particles (as a decimal, where 1.0 is the entire set, 0.8 is 80%. 0.5 is 50%, etc); where to make the cutoff to divide the set into two groups. For example, if you specify --cutoff=0.2, the 20% of particles with the highest correlation scores will be bundled into the first group, and the remaining 80% into the second group.", default=None)
	parser.add_argument("--groups", type=int,help="Number of groups you want the data to be divided into.", default=None)
	parser.add_argument("--topn", type=int,help="Number of particles from best to worst that you want to be written as a substack,averaged, and generate a coordinates .txt file with their coordinates.", default=None)
	parser.add_argument("--sigmaprune", type=float,help="Number of standard deviations to below the mean cut off particles; that is, the mean cross correlation coefficient of all particles will be computed, and those that are --sigmaprune=N standard deviations below the mean will not be considered.",default=None)
	(options,args)=parser.parse_args()
	
	
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

	data=options.input
	cutoff=options.cutoff
		
	print "I have read the parameters"
	
	scores=[]
	dataset=[]
	#x=[]

	n=EMUtil.get_image_count(data)
	
	print "The set has these many particles", n
	
	for i in range(n):
		a=EMData(data,i)
		hdr=a.get_attr_dict()
		print  "analyzing header for particle", i
		if 'spt_score' in hdr:
			score=-1*a['spt_score']
		elif 'spt_coefficient' in hdr:
			score=a['spt_coefficient']*-1
		else:
			print "No score found in header. Terminating!"
			sys.exit()
		scores.append(score)
		dataset.append({'score':float(score), 'particle':a})
		#x.append(i)

	dataset = sorted(dataset, key=itemgetter('score'), reverse=True)

	newscores=[]
	if options.sigmaprune:

		print "Will analyze scores to remove aberrantly low ones"
	
		mu=numpy.mean(scores)
		sigma=numpi.std(scores)
	
		print "Mean is", mu
		print "Std is", sigma
		filter = mean - sigma * (options.prune)
		
		#for s in scores:
		#	if s < halfscore/2.0:
		#		scores.remove(s)
				#x=x[0:-1]
				#print "X is", x
		#		print "I have removed this aberrantly low score", s

		for d in dataset:
			if float(d['score']) < float(filter):
				dataset.remove(d)
				print "I have removed this aberrant particle from the dataset due to its low score", d['score']
			else:
				newscores.append(float(d['score']))
				print "This score is high enough to survive", d['score']
		#for d in dataset:
		#	if d['score'] < halfscore/2.0:
		#		dataset.remove(d)
		#		print "I have removed this aberrant particle from the dataset due to its low score", d['score']

	newN=len(dataset)
	
	scores = newscores
	scores.sort()
	scores.reverse()
	
	if options.groups:
		halfptclnum=int(round(n/2.0))
		halfptcl=dataset[halfptclnum]
		print "halfptcl is", halfptcl
		halfscore=halfptcl['score']
		
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
			groupname = data.replace('.hdf', '_'+ str(kk).zfill(len(str(options.groups))) + '.hdf')
			particleLIST = []			
			for element in subdataset:
				particle = element['particle']
				particle.write_image(groupname,jj)
				jj+=1
				particleLIST.append(particle)
			averageNAME = groupname.replace('.hdf','_AVG.hdf')		
			average = sum(particleLIST)/len(particleLIST)
			average.process_inplace('normalize')
			print "\nThe group %d has been written to %s, and has these many particles in it %d" % ( kk, groupname, len(particleLIST) )
			average.write_image(averageNAME,0)
			print "\nThe group average has been written to", averageNAME
			kk+=1	
			
	if options.cutoff:
		threshptclnum=int(round(newN/(1/options.cutoff)))
		threshptcl=dataset[threshptclnum]	
		print "The new threshptcl is", threshptcl
		threshscore=threshptcl['score']

		threshlabel="%0.2f" %(cutoff)
		threshlabel=str(threshlabel).replace('.','p')	
		group1=[]
		group2=[]
		kk=0
		for i in dataset:
				
			if i['score'] >= threshscore:
				group1.append(i['particle'])
				i['particle'].write_image(g1name,k1)
				k1+=1
			else:
				if ['score'] > halfscore/2.0:
					group2.append(i['particle'])
					i['particle'].write_image(g2name,k2)
					k2+=1
				else:
					print "\n\n@@@@ Found a garbage particle, and thus did not consider it!\n\n"	

		if group1:	
			g1avg = sum(group1)/len(group1)	
			g1avg.write_image(g1name.replace('.hdf','_avg.hdf'),0)

		if group2:
			g2avg = sum(group2)/len(group2)
			g2avg.write_image(g2name.replace('.hdf','_avg.hdf'),0)
	
	if options.topn:
		topndataset = dataset[0:options.topn]
		outnamestack = os.path.basename(options.input).replace('.','top' + str(options.topn) + '.')
		coordsname = outnamestack.split('.')[0] + '_coords.txt'
		indxsname = coordsname.replace('coords','indxs')
		k=0
		linescoords=[]
		linesindxs=[]
		for ptcl in topndataset:
			p = topndataset['particle']
			p.write_image(outnamestack,k)
			linescoords.append(str(p['ptcl_source_coord'])+ '\n')
			linescoords.append(str(p['ptcl_indx'])+ '\n')
			k+=1
		f=open(coordsname,'w')
		f.writelines(linescoords)
		f.close()
	
		f=open(linesindxs,'w')
		f.writelines()
		f.close()
	
	x=range(len(scores))
	
	plot_name = data.replace('.hdf', '_SCORES.png')
	plt.plot(x, scores, marker='+', color='r', linewidth=2)
	plt.title(plot_name)
	plt.ylabel('score')
	plt.xlabel('n ptcl')
	a = plt.gca()

	#a.set_xlim(1,mults[-1])
	#a.set_ylim(0,max(x))

	plt.savefig(plot_name)
	plt.clf()
	
	E2end(logger)

	return()

if __name__ == '__main__':
	main()
