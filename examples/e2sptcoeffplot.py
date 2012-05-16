#!/usr/bin/env python

#
# Author: Jesus Galaz, 04/xx/2012 - using code and concepts drawn from Jesus Galaz's scripts
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
from sys import argv
from operator import itemgetter	
import matplotlib.pyplot as plt
import sys

def main():
	print "I have entered main"

	progname = os.path.basename(sys.argv[0])
	usage = """Aligns a 3d volume to another by executing e2classaverage3d.py and then calculates the FSC between them by calling e2proc3d.py . It returns both a number for the resolution based on the FSC0.5 
	criterion(on the screen) and a plot as an image in .png format."""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_argument("--input", type=str, help="Stack of 3D volumes that have been aligned with EMAN2 SPT programs. The format MUST be '.hdf', and the files must contain the header parameter 'spt_coefficient'.", default=None)
	parser.add_argument("--cutoff", type=float, help="Fraction of particles (as a decimal, where 1.0 is the entire set, 0.8 is 80%. 0.5 is 50%, etc); where to make the cutoff to divide the set into two groups. For example, if you specify --cutoff=0.2, the 20% of particles with the highest correlation scores will be bundled into the first group, and the remaining 80% into the second group.", default=None)
	parser.add_argument("--groups", type=int, help="Number of groups you want the data to be divided into.", default=None)

	(options, args) = parser.parse_args()
	
	data=options.input
	cutoff=options.cutoff
	
	print "I have read the parameters"
	
	scores=[]
	dataset=[]
	x=[]

	n=EMUtil.get_image_count(data)
	
	print "The set has these many particles", n
	
	for i in range(n):
		a=EMData(data,i)
		hdr=a.get_attr_dict()
		print  "analyzing header for particle", i
		if 'spt_score' in hdr:
			score=-1*a['spt_score']
		elif 'spt_coefficient' in hdr:
			score=a['spt_coefficient']
		else:
			print "No score found in header. Terminating!"
			sys.exit()
		scores.append(score)
		dataset.append({'score':float(score), 'particle':a})
		x.append(i)

	dataset = sorted(dataset, key=itemgetter('score'), reverse=True)

	halfptclnum=int(round(n/2.0))
	halfptcl=dataset[halfptclnum]
	print "halfptcl is", halfptcl
	halfscore=halfptcl['score']

	print "Will analyze scores to remove aberrantly low ones"
	
	for s in scores:
		if s < halfscore/2.0:
			scores.remove(s)
			x=x[0:-1]
			print "I have removed this aberrantly low score", s

	for d in dataset:
		if d['score'] < halfscore/2.0:
			dataset.remove(d)
			print "I have removed this aberrant particle from the dataset due to its low score", d['score']

	for d in dataset:
		if d['score'] < halfscore/2.0:
			dataset.remove(d)
			print "I have removed this aberrant particle from the dataset due to its low score", d['score']

	newN=len(dataset)

	scores.sort()
	scores.reverse()
	
	if options.groups and options.cutoff:
		print "ERROR; you cannot specify both --cutoff and --groups at the same time"
		sys.exit()

	if options.groups and not options.cutoff:
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
			
	if options.cutoff and not options.group:
	
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

	plot_name = data.replace('.hdf', '_SCORES.png')
	plt.plot(x, scores, linewidth=1, marker='o', linestyle='--', color='r')
	plt.title(plot_name)
	plt.ylabel('score')
	plt.xlabel('n ptcl')
	a = plt.gca()

	#a.set_xlim(1,mults[-1])
	#a.set_ylim(0,max(x))

	plt.savefig(plot_name)
	plt.clf()

if __name__ == '__main__':
	main()
