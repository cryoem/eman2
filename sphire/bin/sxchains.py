#!/usr/bin/env python
from __future__ import print_function

#
# Author: Steven Ludtke, 04/10/2003 (sludtke@bcm.edu)
# Copyright (c) 2000-2006 Baylor College of Medicine
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
#

import EMAN2_cppwrap
import alignment
import copy
import fundamentals
import global_def
import numpy
import numpy.random
import optparse
import os
import statistics
import sys
import time
import utilities
pass#IMPORTIMPORTIMPORT import EMAN2
pass#IMPORTIMPORTIMPORT import EMAN2_cppwrap
pass#IMPORTIMPORTIMPORT import EMAN2jsondb
pass#IMPORTIMPORTIMPORT import alignment
pass#IMPORTIMPORTIMPORT import copy
pass#IMPORTIMPORTIMPORT import fundamentals
pass#IMPORTIMPORTIMPORT import global_def
pass#IMPORTIMPORTIMPORT import math
pass#IMPORTIMPORTIMPORT import numpy
pass#IMPORTIMPORTIMPORT import numpy.random
pass#IMPORTIMPORTIMPORT import optparse
pass#IMPORTIMPORTIMPORT import os
pass#IMPORTIMPORTIMPORT import pyemtbx.options
pass#IMPORTIMPORTIMPORT import random
pass#IMPORTIMPORTIMPORT import statistics
pass#IMPORTIMPORTIMPORT import sys
pass#IMPORTIMPORTIMPORT import time
pass#IMPORTIMPORTIMPORT import utilities
from builtins import range
pass#IMPORTIMPORTIMPORT import	global_def
pass#IMPORTIMPORTIMPORT from	global_def 	import *
pass#IMPORTIMPORTIMPORT from	EMAN2 		import EMUtil, parsemodopt, EMAN2Ctf
pass#IMPORTIMPORTIMPORT from    EMAN2jsondb import js_open_dict

pass#IMPORTIMPORTIMPORT from	utilities 	import *
pass#IMPORTIMPORTIMPORT from    statistics import mono
pass#IMPORTIMPORTIMPORT import numpy.random
pass#IMPORTIMPORTIMPORT import  os

"""
	Traveling salesman problem solved using Simulated Annealing.
"""
#from scipy import *
#from pylab import *
pass#IMPORTIMPORTIMPORT from numpy import *

def Distance(i1, i2, lccc):
	return max(1.0 - lccc[statistics.mono(i1,i2)][0], 0.0)
# 	return sqrt((R1[0]-R2[0])**2+(R1[1]-R2[1])**2)

def TotalDistance(city, lccc):
	dist = 0.0
	for i in range(len(city)-1):
		dist += Distance(city[i], city[i+1], lccc)
	dist += Distance(city[-1], city[0], lccc)
	return dist
    
def reverse(city, n):
    nct = len(city)
    nn = (1+ ((n[1]-n[0]) % nct))/2 # half the lenght of the segment to be reversed
    # the segment is reversed in the following way n[0]<->n[1], n[0]+1<->n[1]-1, n[0]+2<->n[1]-2,...
    # Start at the ends of the segment and swap pairs of cities, moving towards the center.
    for j in range(nn):
        k = (n[0]+j) % nct
        l = (n[1]-j) % nct
        (city[k],city[l]) = (city[l],city[k])  # swap

def transpt(city, n):
    nct = len(city)
    
    newcity=[]
    # Segment in the range n[0]...n[1]
    for j in range( (n[1]-n[0])%nct + 1):
        newcity.append(city[ (j+n[0])%nct ])
    # is followed by segment n[5]...n[2]
    for j in range( (n[2]-n[5])%nct + 1):
        newcity.append(city[ (j+n[5])%nct ])
    # is followed by segment n[3]...n[4]
    for j in range( (n[4]-n[3])%nct + 1):
        newcity.append(city[ (j+n[3])%nct ])
    return newcity

"""Multiline Comment0"""

def tsp(lccc):

	#     ncity = 100        # Number of cities to visit
	pass#IMPORTIMPORTIMPORT from math import sqrt
	ncity = int( (1+numpy.sqrt(1+8*len(lccc)))/2 )        # Number of cities to visit
    #  sanity check
	if( ncity*(ncity-1)/2 != len(lccc) ): return [-1]

	maxTsteps = 100    # Temperature is lowered not more than maxTsteps
	Tstart = 0.2       # Starting temperature - has to be high enough
	fCool = 0.9        # Factor to multiply temperature at each cooling step
	maxSteps = 100*ncity     # Number of steps at constant temperature
	maxAccepted = 10*ncity   # Number of accepted steps at constant temperature

	Preverse = 0.5      # How often to choose reverse/transpose trial move


	# The index table -- the order the cities are visited.
	city = list(range(ncity))
	# Distance of the travel at the beginning
	dist = TotalDistance(city, lccc)

	#  Not clear what is nct
	nct = ncity
	# Stores points of a move
	n = numpy.zeros(6, dtype=int)

	T = Tstart # temperature

	#     Plot(city, R, dist)

	for t in range(maxTsteps):  # Over temperature

		accepted = 0
		for i in range(maxSteps): # At each temperature, many Monte Carlo steps
            
			while True: # Will find two random cities sufficiently close by
				# Two cities n[0] and n[1] are choosen at random
				n[0] = int((nct)*numpy.random.rand())     # select one city
				n[1] = int((nct-1)*numpy.random.rand())   # select another city, but not the same
				if (n[1] >= n[0]): n[1] += 1   #
				if (n[1] < n[0]): (n[0],n[1]) = (n[1],n[0]) # swap, because it must be: n[0]<n[1]
				nn = (n[0]+nct -n[1]-1) % nct  # number of cities not on the segment n[0]..n[1]
				if nn>=3: break
        
			# We want to have one index before and one after the two cities
			# The order hence is [n2,n0,n1,n3]
			n[2] = (n[0]-1) % nct  # index before n0  -- see figure in the lecture notes
			n[3] = (n[1]+1) % nct  # index after n2   -- see figure in the lecture notes
            
			if Preverse > numpy.random.rand(): 
				# Here we reverse a segment
				# What would be the cost to reverse the path between city[n[0]]-city[n[1]]?
				de = Distance(city[n[2]], city[n[1]], lccc) + Distance(city[n[3]], city[n[0]], lccc)\
					 - Distance(city[n[2]], city[n[0]], lccc) - Distance(city[n[3]] ,city[n[1]], lccc)
                
				if de<0 or numpy.exp(-de/T)>numpy.random.rand(): # Metropolis
					accepted += 1
					dist += de
					reverse(city, n)
			else:
				# Here we transpose a segment
				nc = (n[1]+1+ int(numpy.random.rand()*(nn-1)))%nct  # Another point outside n[0],n[1] segment. See picture in lecture nodes!
				n[4] = nc
				n[5] = (nc+1) % nct

				# Cost to transpose a segment
				de = -Distance( city[n[1]], city[n[3]], lccc) - Distance( city[n[0]], city[n[2]], lccc) \
						- Distance( city[n[4]], city[n[5]], lccc)
				de += Distance( city[n[0]], city[n[4]], lccc) + Distance( city[n[1]], city[n[5]], lccc) \
						+ Distance( city[n[2]], city[n[3]], lccc)

				if de<0 or numpy.exp(-de/T)>numpy.random.rand(): # Metropolis
					accepted += 1
					dist += de
					city = transpt(city, n)
                    
			if accepted > maxAccepted: break

		# Plot
		#         Plot(city, R, dist)
            
		###print "T=%10.5f , distance= %10.5f , accepted steps= %d" %(T, dist, accepted)
		T *= fCool             # The system is cooled down
		if accepted == 0: break  # If the path does not want to change any more, we can stop

        
#     Plot(city, R, dist)
	return city




def main():
	pass#IMPORTIMPORTIMPORT import sys
	pass#IMPORTIMPORTIMPORT import os
	pass#IMPORTIMPORTIMPORT import math
	pass#IMPORTIMPORTIMPORT import random
	pass#IMPORTIMPORTIMPORT import pyemtbx.options
	pass#IMPORTIMPORTIMPORT from   random   import random, seed, randint
	pass#IMPORTIMPORTIMPORT from   optparse import OptionParser
	pass#IMPORTIMPORTIMPORT from global_def import ERROR

	progname = os.path.basename(sys.argv[0])
	usage = progname + """ [options] <inputfile> <outputfile>

	Forms chains of 2D images based on their similarities.

	Functionality:


	Order a 2-D stack of image based on pair-wise similarity (computed as a cross-correlation coefficent).
		Options 1-3 require image stack to be aligned.  The program will apply orientation parameters if present in headers.
	    The ways to use the program:
	   1.  Use option initial to specify which image will be used as an initial seed to form the chain.
	        sxchains.py input_stack.hdf output_stack.hdf --initial=23 --radius=25

	   2.  If options initial is omitted, the program will determine which image best serves as initial seed to form the chain
	        sxchains.py input_stack.hdf output_stack.hdf --radius=25

	   3.  Use option circular to form a circular chain.
	        sxchains.py input_stack.hdf output_stack.hdf --circular--radius=25

	   4.  New circular code based on pairwise alignments
			sxchains.py aclf.hdf chain.hdf circle.hdf --align  --radius=25 --xr=2 --pairwiseccc=lcc.txt

	   5.  Circular ordering based on pairwise alignments
			sxchains.py vols.hdf chain.hdf mask.hdf --dd  --radius=25


"""

	parser = optparse.OptionParser(usage,version=global_def.SPARXVERSION)
	parser.add_option("--dd", action="store_true", help="Circular ordering without adjustment of orientations", default=False)
	parser.add_option("--circular", action="store_true", help="Select circular ordering (first image has to be similar to the last)", default=False)
	parser.add_option("--align", action="store_true", help="Compute all pairwise alignments and from the table of image similarities find the best chain", default=False)
	parser.add_option("--initial", type="int", default=-1, help="Specifies which image will be used as an initial seed to form the chain. (default = 0, means the first image)")
	parser.add_option("--radius", type="int", default=-1, help="Radius of a circular mask for similarity based ordering")
	#  import params for 2D alignment
	parser.add_option("--ou",           type="int",    default=-1,          help="outer radius for 2D alignment < nx/2-1 (set to the radius of the particle)")
	parser.add_option("--xr",           type="int",    default=0,     		help="range for translation search in x direction, search is +/xr (0)")
	parser.add_option("--yr",           type="int",    default=0,          	help="range for translation search in y direction, search is +/yr (0)")
	#parser.add_option("--nomirror",     action="store_true", default=False,   help="Disable checking mirror orientations of images (default False)")
	parser.add_option("--pairwiseccc",  type="string",	default= " ",      help="Input/output pairwise ccc file")


	(options, args) = parser.parse_args()

	global_def.BATCH = True

					
	if options.dd:
		nargs = len(args)
		if nargs != 3:
			print("must provide name of input and two output files!")
			return
		stack = args[0]
		new_stack = args[1]


		pass#IMPORTIMPORTIMPORT from utilities import model_circle
		pass#IMPORTIMPORTIMPORT from statistics import ccc
		pass#IMPORTIMPORTIMPORT from statistics import mono
		lend = EMAN2_cppwrap.EMUtil.get_image_count(stack)
		lccc = [None]*(lend*(lend-1)/2)

		for i in range(lend-1):
			v1 = utilities.get_im( stack, i )
			if( i == 0 and nargs == 2):
				nx = v1.get_xsize()
				ny = v1.get_ysize()
				nz = v1.get_ysize()
				if options.ou < 1 : radius = nx//2-2
				else:  radius = options.ou
				mask = utilities.model_circle(radius, nx, ny, nz)
			else:
				mask = utilities.get_im(args[2])
				
			for j in range(i+1, lend):
				lccc[statistics.mono(i,j)] = [statistics.ccc(v1, utilities.get_im( stack, j ), mask), 0]


		order = tsp(lccc)
		if(len(order) != lend):
			print(" problem with data length")
			pass#IMPORTIMPORTIMPORT from sys import exit
			exit()
		print("Total sum of cccs :",TotalDistance(order, lccc))
		print("ordering :",order)
		for i in range(lend):  utilities.get_im(stack, order[i]).write_image( new_stack, i )

	elif options.align:
		nargs = len(args)
		if nargs != 3:
			print("must provide name of input and two output files!")
			return

		pass#IMPORTIMPORTIMPORT from utilities import get_params2D, model_circle
		pass#IMPORTIMPORTIMPORT from fundamentals import rot_shift2D
		pass#IMPORTIMPORTIMPORT from statistics import ccc
		pass#IMPORTIMPORTIMPORT from time import time
		pass#IMPORTIMPORTIMPORT from alignment import align2d, align2d_scf
		
		stack = args[0]
		new_stack = args[1]
		
		d = EMAN2_cppwrap.EMData.read_images(stack)
		if(len(d)<6):
			global_def.ERROR("Chains requires at least six images in the input stack to be executed", "sxchains.py", 1)

		"""Multiline Comment1"""

		nx = d[0].get_xsize()
		ny = d[0].get_ysize()
		if options.ou < 1 : radius = nx//2-2
		else:  radius = options.ou
		mask = utilities.model_circle(radius, nx, ny)

		if(options.xr < 0):	xrng = 0
		else:				xrng = options.xr
		if(options.yr < 0):	yrng = xrng
		else:				yrng = options.yr
			
		initial = max(options.initial, 0)

		pass#IMPORTIMPORTIMPORT from statistics import mono
		lend = len(d)
		lccc = [None]*(lend*(lend-1)/2)
		pass#IMPORTIMPORTIMPORT from utilities import read_text_row

		if  options.pairwiseccc == " " or not os.path.exists(options.pairwiseccc) :
			st = time.time()
			for i in range(lend-1):
				for j in range(i+1, lend):
					#  j>i meaning mono entry (i,j) or (j,i) indicates T i->j (from smaller index to larger)
					#alpha, sx, sy, mir, peak = align2d(d[i],d[j], xrng, yrng, step=options.ts, first_ring=options.ir, last_ring=radius, mode = "F")
					alpha, sx, sy, mir, peak = alignment.align2d_scf(d[i],d[j], xrng, yrng, ou=radius)
					lccc[statistics.mono(i,j)] = [statistics.ccc(d[j], fundamentals.rot_shift2D(d[i], alpha, sx, sy, mir, 1.0), mask), alpha, sx, sy, mir]
				#print "  %4d   %10.1f"%(i,time()-st)

			if((not os.path.exists(options.pairwiseccc)) and (options.pairwiseccc != " ")):
				pass#IMPORTIMPORTIMPORT from utilities import write_text_row
				utilities.write_text_row([[initial,0,0,0,0]]+lccc, options.pairwiseccc)
		elif(os.path.exists(options.pairwiseccc)):
			lccc = utilities.read_text_row(options.pairwiseccc)
			initial = int(lccc[0][0] + 0.1)
			del lccc[0]


		for i in range(len(lccc)):
			T = EMAN2_cppwrap.Transform({"type":"2D","alpha":lccc[i][1],"tx":lccc[i][2],"ty":lccc[i][3],"mirror":int(lccc[i][4]+0.1)})
			lccc[i] = [lccc[i][0],T]

		tdummy = EMAN2_cppwrap.Transform({"type":"2D"})
		maxsum = -1.023
		for m in range(0,lend):#initial, initial+1):
			indc = list(range( lend))
			lsnake = [[m, tdummy, 0.0]]
			del indc[m]

			lsum = 0.0
			while len(indc) > 1:
				maxcit = -111.
				for i in range(len(indc)):
						cuc = lccc[statistics.mono(indc[i], lsnake[-1][0])][0]
						if cuc > maxcit:
								maxcit = cuc
								qi = indc[i]
								#  Here we need transformation from the current to the previous,
								#     meaning indc[i] -> lsnake[-1][0]
								T = lccc[statistics.mono(indc[i], lsnake[-1][0])][1]
								#  If direction is from larger to smaller index, the transformation has to be inverted
								if( indc[i] > lsnake[-1][0] ):  T = T.inverse()
								
								
				lsnake.append([qi,T, maxcit])
				lsum += maxcit

				del indc[indc.index(qi)]

			T = lccc[statistics.mono(indc[-1], lsnake[-1][0])][1]
			if( indc[-1] > lsnake[-1][0]):  T = T.inverse()
			lsnake.append([indc[-1], T, lccc[statistics.mono(indc[-1], lsnake[-1][0])][0]])
			print(" initial image and lsum  ",m,lsum)
			#print lsnake
			if(lsum > maxsum):
				maxsum = lsum
				init = m
				snake = [lsnake[i] for i in range(lend)]
		print("  Initial image selected : ",init,maxsum,"    ",TotalDistance([snake[m][0] for m in range(lend)], lccc))
		#for q in snake: print q

		pass#IMPORTIMPORTIMPORT from copy import deepcopy
		trans=copy.deepcopy([snake[i][1] for i in range(len(snake))])
		print([snake[i][0] for i in range(len(snake))])
		"""Multiline Comment2"""
		for k in range(lend-2,0,-1):
			T = snake[k][1]
			for i in range(k+1, lend):
 					trans[i] = T*trans[i]
		#  To add - apply all transformations and do the overall centering.
		for m in range(lend):
			prms = trans[m].get_params("2D")
			#print " %3d   %7.1f   %7.1f   %7.1f   %2d  %6.2f"%(snake[m][0], prms["alpha"], prms["tx"], prms["ty"], prms["mirror"], snake[m][2])
			#rot_shift2D(d[snake[m][0]], prms["alpha"], prms["tx"], prms["ty"], prms["mirror"]).write_image(new_stack, m)
			fundamentals.rot_shift2D(d[snake[m][0]], prms["alpha"], 0.0,0.0, prms["mirror"]).write_image(new_stack, m)

		order = tsp(lccc)
		if(len(order) != lend):
			print(" problem with data length")
			pass#IMPORTIMPORTIMPORT from sys import exit
			exit()
		print(TotalDistance(order, lccc))
		print(order)
		ibeg = order.index(init)
		order = [order[(i+ibeg)%lend] for i in range(lend)]
		print(TotalDistance(order, lccc))
		print(order)


		snake = [tdummy]
		for i in range(1,lend):
			#  Here we need transformation from the current to the previous,
			#     meaning order[i] -> order[i-1]]
			T = lccc[statistics.mono(order[i], order[i-1])][1]
			#  If direction is from larger to smaller index, the transformation has to be inverted
			if( order[i] > order[i-1] ):  T = T.inverse()
			snake.append(T)
		assert(len(snake) == lend)
		pass#IMPORTIMPORTIMPORT from copy import deepcopy
		trans = copy.deepcopy(snake)
		for k in range(lend-2,0,-1):
			T = snake[k]
			for i in range(k+1, lend):
 					trans[i] = T*trans[i]

		#  Try to smooth the angles - complicated, I am afraid one would have to use angles forward and backwards
		#     and find their average??
		#  In addition, one would have to recenter them
		"""Multiline Comment3"""

		for m in range(lend):
			prms = trans[m].get_params("2D")
			#rot_shift2D(d[order[m]], prms["alpha"], prms["tx"], prms["ty"], prms["mirror"]).write_image("metro.hdf", m)
			fundamentals.rot_shift2D(d[order[m]], prms["alpha"], 0.0,0.0, prms["mirror"]).write_image(args[2], m)

		"""Multiline Comment4"""
		"""Multiline Comment5"""

		"""Multiline Comment6"""
	else:
		nargs = len(args)
		if nargs != 2:
			print("must provide name of input and output file!")
			return
		
		pass#IMPORTIMPORTIMPORT from utilities import get_params2D, model_circle
		pass#IMPORTIMPORTIMPORT from fundamentals import rot_shift2D
		pass#IMPORTIMPORTIMPORT from statistics import ccc
		pass#IMPORTIMPORTIMPORT from time import time
		pass#IMPORTIMPORTIMPORT from alignment import align2d
		
		stack = args[0]
		new_stack = args[1]
		
		d = EMAN2_cppwrap.EMData.read_images(stack)
		try:
			print("Using 2D alignment parameters from header.")
			ttt = d[0].get_attr('xform.params2d')
			for i in range(len(d)):
				alpha, sx, sy, mirror, scale = utilities.get_params2D(d[i])
				d[i] = fundamentals.rot_shift2D(d[i], alpha, sx, sy, mirror)
		except:
			pass

		nx = d[0].get_xsize()
		ny = d[0].get_ysize()
		if options.radius < 1 : radius = nx//2-2
		else:  radius = options.radius
		mask = utilities.model_circle(radius, nx, ny)

		init = options.initial
		
		if init > -1 :
			print("      initial image: %d" % init)
			temp = d[init].copy()
			temp.write_image(new_stack, 0)
			del d[init]
			k = 1
			lsum = 0.0
			while len(d) > 1:
				maxcit = -111.
				for i in range(len(d)):
						cuc = statistics.ccc(d[i], temp, mask)
						if cuc > maxcit:
								maxcit = cuc
								qi = i
				# 	print k, maxcit
				lsum += maxcit
				temp = d[qi].copy()
				del d[qi]
				temp.write_image(new_stack, k)
				k += 1
			print(lsum)
			d[0].write_image(new_stack, k)
		else:			
			if options.circular :
				print("Using options.circular, no alignment")
				#  figure the "best circular" starting image
				maxsum = -1.023
				for m in range(len(d)):
					indc = list(range(len(d)))
					lsnake = [-1]*(len(d)+1)
					lsnake[0]  = m
					lsnake[-1] = m
					del indc[m]
					temp = d[m].copy()
					lsum = 0.0
					direction = +1
					k = 1
					while len(indc) > 1:
						maxcit = -111.
						for i in range(len(indc)):
								cuc = statistics.ccc(d[indc[i]], temp, mask)
								if cuc > maxcit:
										maxcit = cuc
										qi = indc[i]
						lsnake[k] = qi
						lsum += maxcit
						del indc[indc.index(qi)]
						direction = -direction
						for i in range( 1,len(d) ):
							if( direction > 0 ):
								if(lsnake[i] == -1):
									temp = d[lsnake[i-1]].copy()
									#print  "  forw  ",lsnake[i-1]
									k = i
									break
							else:
								if(lsnake[len(d) - i] == -1):
									temp = d[lsnake[len(d) - i +1]].copy()
									#print  "  back  ",lsnake[len(d) - i +1]
									k = len(d) - i
									break

					lsnake[lsnake.index(-1)] = indc[-1]
					#print  " initial image and lsum  ",m,lsum
					#print lsnake
					if(lsum > maxsum):
						maxsum = lsum
						init = m
						snake = [lsnake[i] for i in range(len(d))]
				print("  Initial image selected : ",init,maxsum)
				print(lsnake)
				for m in range(len(d)):  d[snake[m]].write_image(new_stack, m)
			else:
				#  figure the "best" starting image
				print("Straight chain, no alignment")
				maxsum = -1.023
				for m in range(len(d)):
					indc = list(range(len(d)))
					lsnake = [m]
					del indc[m]
					temp = d[m].copy()
					lsum = 0.0
					while len(indc) > 1:
						maxcit = -111.
						for i in range(len(indc)):
								cuc = statistics.ccc(d[indc[i]], temp, mask)
								if cuc > maxcit:
										maxcit = cuc
										qi = indc[i]
						lsnake.append(qi)
						lsum += maxcit
						temp = d[qi].copy()
						del indc[indc.index(qi)]

					lsnake.append(indc[-1])
					#print  " initial image and lsum  ",m,lsum
					#print lsnake
					if(lsum > maxsum):
						maxsum = lsum
						init = m
						snake = [lsnake[i] for i in range(len(d))]
				print("  Initial image selected : ",init,maxsum)
				print(lsnake)
				for m in range(len(d)):  d[snake[m]].write_image(new_stack, m)
		

if __name__ == "__main__":
	main()
