#!/usr/bin/env python

'''
====================
Author: Jesus Galaz - 2012, Last update: 08/2012
====================

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
'''

from optparse import OptionParser
from EMAN2 import *
from sys import argv
import sys
import os
import EMAN2
#import heapq
#import operator
#import random
import numpy

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """%prog [options]

	This program calculates the angular distance between the orientation of simulated subtomograms that have been generated with e2spt_simulation.py,
	and the 'solution' alignment found when aligning them to the model where they came from, using e2spt_classaverage.py 
	"""
			
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)	
	
	parser.add_argument("--sym",type=str,default=None,help="""If the particles are symmetric, provide the symmetry. 
								If the particle is symmetric and you didn't select the --sym option when simulating subtomograms from it,
								then, multiple solutions would be valid, but they would have a different absolute distance to the orientation
								in which the particle was randomized. 
								e2spt_classaverage.py would have aligned the particles to the closest asymmetric unit 
								and this needs to be taken into account to determine whether it found 'a' correct solution or not.""")

	parser.add_argument("--input", type=str, help="""The name of the input stack in .hdf format containing the ALIGNED PARTICLES after aligning the simulated 
							subtomograms with e2spt_classaverage.py to the reference from which they were generated. The particle must contain
							the header parameters 'spt_randT' and 'spt_ali_params' """, default=None)
	parser.add_argument("--output", type=str, help="""The name of the output text file in .txt format containing the average angular distance and translational distance of the solutions proposed
							by e2spt_classaverage3d.py respect to the simulated transforms.""", default=None) 
	
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n",type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	
	
	(options, args) = parser.parse_args()
	
	logger = E2init(sys.argv, options.ppid)
	
	if  not options.input:
		print "ERRO. You must supply and input stack."
		sys.exit()
	
	nptcls = EMUtil.get_image_count(options.input)
	
	print "I will calculate the distance of these many particles", nptcls
	
	angles=[]
	translations=[]
	
	for i in range(nptcls):
		print "I am working on ptcl number", i
		ptcl = EMData(options.input,i, True)
		simT = ptcl['sptsim_randT']
		solutionT = ptcl['spt_ali_param']	#If the solution IS an aproximate solution, its inverse should be almost equal to simT, and thus
							#the distance between simT and solutionT.inverse() should be small.
		solutionTinverse = solutionT.inverse()
		
		angles.append(angdist(simT,solutionTinverse))
	
		transX = pow(simT.get_trans()[0] - solutionTinverse.get_trans()[0],2)
		transY = pow(simT.get_trans()[1] - solutionTinverse.get_trans()[1],2)
		transZ = pow(simT.get_trans()[2] - solutionTinverse.get_trans()[2],2)
		trans = sqrt(transX + transY + transZ)
		print "The translational distance is", trans
		translations.append(trans)	
	
	avgA = sum(angles)/len(angles)
	avgT = sum(translations)/len(translations)
	print "The average angular and translational ditances are", avgA, avgT
	resultsfile = options.input.replace('.hdf', '_RESULTS.txt')
	if options.output:
		resultsfile = options.output
	
	f=open(resultsfile,'w')
	line='average angular distance=' + str(avgA) + '\n'+ 'average translational error=' + str(avgT) +'\n'
	f.write(line)
	f.close()
	
	
	return()
	
'''
ANGDIST - Calculates the distance between two transforms. Used to verify the answers provided by the aligner for simmulated data,
and to determine the proximity of the 10 best peaks proposed, for experimental and simmulated data.
'''
def angdist(t1,t2):
	t2inv = t2.inverse()
	product = t2inv * t1
	product_SPIN = product.get_rotation('spin')
	angular_distance = round(float(product_SPIN["Omega"]),2)
	print "The angular distance is", angular_distance
	
	return(angular_distance)


if __name__ == '__main__':
	main()
