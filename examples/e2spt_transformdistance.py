#!/usr/bin/env python

'''
====================
Author: Jesus Galaz - 2012, Last update: 08/2013
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
	usage = """e2spt_transformdistance.py [options]

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
							the header parameters 'spt_randT' and 'xform.align3d' """, default=None)
							
							
	#parser.add_argument("--finalref", type=str, help="""The last reference used in e2spt_classaverage.py if --iter was more than 1.
	#						If the particles were aligned through multiple iterations of refinement in e2spt_classaverage.py
	#						the alignment transform on their headers will be respect to the average that served as the final reference.
	#						You need to align this final reference to the model that was used to simulated the particles, and take THAT relative alignment into account
	#						when determining whether the final alignment of the particles undid the 'random' transform they were rotated into
	#						when simulating the subvolumes with e2spt_simulation.py.""", default=None)
							
							
	#parser.add_argument("--refused", type=str, help="""Original reference/model from which the particles were simulated. """, default=None)
	
	parser.add_argument("--angles",type=str,default='',help="""Two sets of 3 euler angles 
		representing two rotational orientations between which you want to calculate the
		absolute angular distance (the one rotation that, using quaternions, would bring
		you from the first orientation, to the second). The 6 euler angles are provides as 
		--angles=az0,alt0,phi0,az1,alt1,phi1""")
	
	parser.add_argument("--transform",type=str,default='',help="""Transform to apply to particles 
		before measuring error respect to randT on their header.
		You can supply the parameters in the form --transform=type:eman,az:azval,alt:altval,phi:phival,tx:txval,ty:tyval,tz:tzval, 
		with appropriate values for all 6 variables (the right side after az:,alt:,phi:,tx:,ty: and tz:). 
		Note that you only have to provide non-zero variables and type:eman.
		
		Alternatively the transform can be read from the first line of a tomo_xforms.json file, --transform=tomo_xforms.json.
		If the particles were aligned through multiple iterations of refinement in e2spt_classaverage.py
		the alignment transform on their headers will be respect to the average that served as the final reference.
		In that case, you need to supply here the Transform that aligns this 'FINAL REFERENCE' (the penultimate image in the class_00.hdf file
		generated by e2pst_classaverage.py if the --savesteps parameter is on) to the ORIGINAL REFERENCE that was used to simulate the particles, 
		so that such relative alignment is taken into account when determining whether the final alignment of the particles 'undid' the random transforms 
		they were rotated into when simulating the subvolumes with e2spt_simulation.py (these random transforms are stored as 'spt_randT' in the headers of the particles.""")
		
	parser.add_argument("--output", type=str, help="""The name of the output text file in .txt format containing the average angular distance and translational distance of the solutions proposed
							by e2spt_classaverage3d.py respect to the simulated transforms.""", default='results.txt') 
	
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n",type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	
	#parser.add_argument("--shrink", type=int,default=0,help="""If the 'solution transform'
	#	storied in the header parameter 'xform.align3d' was derived from particles that 
	#	were shrunk during alignment, providing --shrinkalign in e2spt_tomosimjobs.py,
	#	then this program needs to know what the shrink factor was, to properly compute
	#	the error between 'xform.align3d' and 'sptsim_randT', the latter containing the
	#	'randomized Transform' that was used to rotate the simulated subtomogram. 
	#	""")
	
	parser.add_argument("--nolog",action="store_true",help="Turn off recording of the command ran for this program onto the .eman2log.txt file",default=False) 
	
	(options, args) = parser.parse_args()
	
	log=0
	if not options.nolog:
		logger = E2init(sys.argv, options.ppid)
		log=1
		
	if  not options.input and not options.angles:
		print """ERROR. You must supply either a 3-D image stack (previously ran through 
		e2spt_classaverage.py using --input, or two rotational orientations using --angles."""
		sys.exit()
	
	if options.input:
		nptcls = EMUtil.get_image_count(options.input)
	
		print "I will calculate the distance of these many particles", nptcls
	
		angularDistances=[]
		translations=[]
	
		compensatoryT = None
	
		if options.transform:
			if '.json' in options.transform:
				js = js_open_dict(options.transform)
				tnums = len(js)
				tomoID = "tomo_" + str(0).zfill( len(str( tnums )) )
				compensatoryT = js[tomoID]
					
				print "compensatoryT is", compensatoryT
				
			elif '.json' not in options.transform:
				params = options.transform.split(',')
				print "These are the supplied transform params", params
				pdict={}
				for p in params:	
					key = p.split(':')[0]
					value = p.split(':')[-1]
					if 'type' not in key:
						value = float( value )
					pdict.update({ key : value })
					
				compensatoryT = Transform(pdict)				#Compensatory T should be the alignment that takes the average/model
																#back to the symmetry axis.
	
		for i in range(nptcls):
			print "\nIn e2spt_transformdistance, I am working on ptcl number", i
			ptcl = EMData(options.input,i, True)
			simT = ptcl['sptsim_randT']
			#solutionT = ptcl['spt_ali_param']	#If the solution IS an aproximate solution, its inverse should be almost equal to simT (or 'spt_randT' in particle headers), and thus
								#the distance between simT and solutionT.inverse() should be small.
			solutionT = ptcl['xform.align3d']
			
			
			
			
			#'''
			#If solutionT was found for particles shrunk during alignment, but not during
			#suntomogram simulation, then the solution is also 'shrunk' respect to the real
			#answer. Therefore, translations must be multiplied by the shrinking factor.
			#'''
			#if options.shrink:
			#	solutionT.set_trans( solutionT.get_trans() * options.shrink )
			
			originalSolutionT = solutionT
			originalSolutionTinverse = originalSolutionT.inverse()
			
			if compensatoryT:
				solutionT = compensatoryT * solutionT
			
			if options.sym:
				solutionT = accountForSym( options, solutionT, simT )
			
			solutionTinverse = solutionT.inverse()
		
			angularDistances.append(angdist(simT,solutionTinverse))
	
			transX = pow(simT.get_trans()[0] - originalSolutionTinverse.get_trans()[0],2)
			transY = pow(simT.get_trans()[1] - originalSolutionTinverse.get_trans()[1],2)
			transZ = pow(simT.get_trans()[2] - originalSolutionTinverse.get_trans()[2],2)
			trans = sqrt(transX + transY + transZ)
			
			print "The translational distance is", trans
			translations.append(trans)	
	
		avgA = sum(angularDistances)/len(angularDistances)
		avgT = sum(translations)/len(translations)
		print "The average angular and translational distances are", avgA, avgT
		resultsfile = options.input.replace('.hdf', '_RESULTS.txt')
		if options.output:
			resultsfile = options.output
	
		f=open(resultsfile,'w')
		line='average angular distance=' + str(avgA) + '\n'+ 'average translational error=' + str(avgT) +'\n'
		f.write(line)
		f.close()
	
	elif options.angles:
		angles=options.angles.split(',')
		
		t1 = Transform({'type':'eman','az':float(angles[0]),'alt':float(angles[1]),'phi':float(angles[2]) })
		t2 = Transform({'type':'eman','az':float(angles[3]),'alt':float(angles[4]),'phi':float(angles[5]) })
		#t2i = t2.inverse()
		aDistance = angdist(t1,t2)
		print "The angular distance between t1=" + str(t1) + " and t2=" + str(t2) + "is: " + str(aDistance)
	
	if not options.nolog and log:
		E2end(logger)
	
	return()
	
	
'''
ANGDIST - Calculates the distance between two transforms. Used to verify the answers provided by the aligner for simmulated data,
and to determine the proximity of the 10 best peaks proposed, for experimental and simmulated data.
'''
def angdist(t1,t2, num=0):
	t2inv = t2.inverse()
	product = t2inv * t1
	product_SPIN = product.get_rotation('spin')
	angular_distance = round(float(product_SPIN["omega"]),2)
	
	print "For orientation %d the angular distance %f" %(num, angular_distance)
	
	return(angular_distance)


def accountForSym( options, solutionT, simT ):
	
	symDistances = {}
	
	symnames = ['oct','OCT','icos','ICOS','tet','TET']
	
	orientations = []
	
	symletter = 'c'
	symnum = 1
	
	if options.sym:
		print "\nsym found", options.sym
		
		if options.sym not in symnames:
			
			symletter = options.sym
			symnum = options.sym
			
			for x in options.sym:
				if x.isalpha():
					symnum = symnum.replace(x,'')
					
			for x in options.sym:
				if x.isdigit():
					symletter = symletter.replace(x,'')
			
			print "\nThe letter for sym is", symletter
			print "\nThe num for sym is", symnum
		
		if options.sym == 'oct' or options.sym == 'OCT':
			symnum = 8
			
		if options.sym == 'tet' or options.sym == 'TET':
			symnum = 12	
			
		if options.sym == 'icos' or options.sym == 'ICOS':
			symnum = 60	
		
		symnum = int(symnum)
			
		t = Transform()
		
		if symnum:
			print "\nSymnum determined",symnum
			
			if symletter == 'd' or symletter == 'D':
				symnum*=2
			
			#if options.save
			
			for i in range(symnum):
				orientation = t.get_sym( options.sym , i )
				orientations.append( orientation )
	
	
	
	groundT = Transform()
	kk=0
	for o in orientations:
		
		oi = o.inverse()
			
		diff = solutionT * simT				#The solution applied to the known random transformation
											#will give you a large "difference" if the solution found
											#is far from the symmetry axis. Otherwise, they'll cancel each other out 
		
		distance = angdist( o,  diff,kk+1 )	 	#Find the distance between this difference and all symmetry-related positions
		
		solutionTsym = oi * solutionT		#For each possible symmetry related position, you calculate a
											#symmetry-related "variant" of the solution by multiplying the 
											#inverse of the symmetry-related orientation times the solution 
											#(in an ideal case, the solution by itself would be the inverse of random transform "simT")
		
		symDistances.update( {distance: solutionTsym } )	#Keep track of all symmetry-related solutions and their transforms
		
		kk+=1
		
	mindist = min( symDistances.keys() )	#Find the symmetry-related solution that yields the minimum distance to a symmetry-related axis
	
	print "\nThe minimum distance is", mindist
	
	bestT = symDistances[mindist]
	print "\nTherefore the best transform is", bestT
	
	return bestT
	

if __name__ == '__main__':
	main()
