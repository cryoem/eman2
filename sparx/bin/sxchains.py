#!/usr/bin/env python

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

import	global_def
from	global_def 	import *
from	EMAN2 		import EMUtil, parsemodopt, EMAN2Ctf
from    EMAN2jsondb import js_open_dict

from	utilities 	import *
from    statistics import mono
import  os

"""
	Traveling salesman problem solved using Simulated Annealing.
"""
from scipy import *
#from pylab import *

def Distance(i1, i2, lccc):
	return max(1.0 - lccc[mono(i1,i2)][0], 0.0)
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

"""
def Plot(city, R, dist):
    # Plot
    Pt = [R[city[i]] for i in range(len(city))]
    Pt += [R[city[0]]]
    Pt = array(Pt)
    title('Total distance='+str(dist))
    plot(Pt[:,0], Pt[:,1], '-o')
    show()
"""

def tsp(lccc):

	#     ncity = 100        # Number of cities to visit
	from math import sqrt
	ncity = int( (1+sqrt(1+8*len(lccc)))/2 )        # Number of cities to visit
    #  sanity check
	if( ncity*(ncity-1)/2 != len(lccc) ): return [-1]

	maxTsteps = 100    # Temperature is lowered not more than maxTsteps
	Tstart = 0.2       # Starting temperature - has to be high enough
	fCool = 0.9        # Factor to multiply temperature at each cooling step
	maxSteps = 100*ncity     # Number of steps at constant temperature
	maxAccepted = 10*ncity   # Number of accepted steps at constant temperature

	Preverse = 0.5      # How often to choose reverse/transpose trial move


	# The index table -- the order the cities are visited.
	city = range(ncity)
	# Distance of the travel at the beginning
	dist = TotalDistance(city, lccc)

	#  Not clear what is nct
	nct = ncity
	# Stores points of a move
	n = zeros(6, dtype=int)

	T = Tstart # temperature

	#     Plot(city, R, dist)

	for t in range(maxTsteps):  # Over temperature

		accepted = 0
		for i in range(maxSteps): # At each temperature, many Monte Carlo steps
            
			while True: # Will find two random cities sufficiently close by
				# Two cities n[0] and n[1] are choosen at random
				n[0] = int((nct)*rand())     # select one city
				n[1] = int((nct-1)*rand())   # select another city, but not the same
				if (n[1] >= n[0]): n[1] += 1   #
				if (n[1] < n[0]): (n[0],n[1]) = (n[1],n[0]) # swap, because it must be: n[0]<n[1]
				nn = (n[0]+nct -n[1]-1) % nct  # number of cities not on the segment n[0]..n[1]
				if nn>=3: break
        
			# We want to have one index before and one after the two cities
			# The order hence is [n2,n0,n1,n3]
			n[2] = (n[0]-1) % nct  # index before n0  -- see figure in the lecture notes
			n[3] = (n[1]+1) % nct  # index after n2   -- see figure in the lecture notes
            
			if Preverse > rand(): 
				# Here we reverse a segment
				# What would be the cost to reverse the path between city[n[0]]-city[n[1]]?
				de = Distance(city[n[2]], city[n[1]], lccc) + Distance(city[n[3]], city[n[0]], lccc)\
					 - Distance(city[n[2]], city[n[0]], lccc) - Distance(city[n[3]] ,city[n[1]], lccc)
                
				if de<0 or exp(-de/T)>rand(): # Metropolis
					accepted += 1
					dist += de
					reverse(city, n)
			else:
				# Here we transpose a segment
				nc = (n[1]+1+ int(rand()*(nn-1)))%nct  # Another point outside n[0],n[1] segment. See picture in lecture nodes!
				n[4] = nc
				n[5] = (nc+1) % nct

				# Cost to transpose a segment
				de = -Distance( city[n[1]], city[n[3]], lccc) - Distance( city[n[0]], city[n[2]], lccc) \
						- Distance( city[n[4]], city[n[5]], lccc)
				de += Distance( city[n[0]], city[n[4]], lccc) + Distance( city[n[1]], city[n[5]], lccc) \
						+ Distance( city[n[2]], city[n[3]], lccc)

				if de<0 or exp(-de/T)>rand(): # Metropolis
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




def pca(cov):
	from numpy import  linalg, argsort
	""" assume one sample per column """
	values, vecs = linalg.eigh(cov)
	perm = argsort(-values)  # sort in descending order
	return values[perm], vecs[:, perm]


def main():
	import sys
	import os
	import math
	import random
	import pyemtbx.options
	import time
	from   random   import random, seed, randint
	from   optparse import OptionParser

	progname = os.path.basename(sys.argv[0])
	usage = progname + """ [options] <inputfile> <outputfile>

	Forms chains of 2D images based on their similarities.

	Functionality:


	Order a 2-D stack of image based on pair-wise similarity (computed as a cross-correlation coefficent).
		Options 1-3 require image stack to be aligned.  The program will apply orientation parameters if present in headers.
	    The ways to use the program:
	   1  Use option initial to specify which image will be used as an initial seed to form the chain.
	        sxchains.py input_stack.hdf output_stack.hdf --initial=23 --radius=25
	   2  If options initial is omitted, the program will determine which image best serves as initial seed to form the chain
	        sxchains.py input_stack.hdf output_stack.hdf --radius=25
	   3  Use option circular to form a circular chain.
	        sxchains.py input_stack.hdf output_stack.hdf --circular--radius=25
	   4  New circular code based on pairwise alignments
			sxchains.py aclf.hdf chain.hdf circle.hdf --align  --radius=25 --xr=2 --pairwiseccc=lcc.txt

	   5  Circular ordering based on pairwise alignments
			sxchains.py vols.hdf chain.hdf mask.hdf --dd  --radius=25


"""

	parser = OptionParser(usage,version=SPARXVERSION)
	parser.add_option("--dd", action="store_true", help="Circular ordering without adjustment of orientations", default=False)
	parser.add_option("--circular", action="store_true", help="Select circular ordering (first image has to be similar to the last)", default=False)
	parser.add_option("--align", action="store_true", help="Compute all pairwise alignments and for the table of their similarities find the best chain", default=False)
	parser.add_option("--initial", type="int", default=-1, help="Specifies which image will be used as an initial seed to form the chain. (default = 0, means the first image)")
	parser.add_option("--radius", type="int", default=-1, help="Radius of a circular mask for similarity based ordering")
	#  import params for 2D alignment
	parser.add_option("--ou",           type="int",    default=-1,          help="outer radius for 2D alignment < nx/2-1 (set to the radius of the particle)")
	parser.add_option("--xr",           type="int",    default=0,     		help="range for translation search in x direction, search is +/xr (0)")
	parser.add_option("--yr",           type="int",    default=0,          	help="range for translation search in y direction, search is +/yr (0)")
	#parser.add_option("--nomirror",     action="store_true", default=False,   help="Disable checking mirror orientations of images (default False)")
	parser.add_option("--pairwiseccc",  type="string",	default= None,      help="Input/output pairwise ccc file")


 	(options, args) = parser.parse_args()

	global_def.BATCH = True

					
	if options.dd:
		nargs = len(args)
		if nargs != 3:
			print "must provide name of input and two output files!"
			return
		stack = args[0]
		new_stack = args[1]


		from utilities import model_circle
		from statistics import ccc
		from statistics import mono
		lend = EMUtil.get_image_count(stack)
		lccc = [None]*(lend*(lend-1)/2)

		for i in xrange(lend-1):
			v1 = get_im( stack, i )
			if( i == 0 and nargs == 2):
				nx = v1.get_xsize()
				ny = v1.get_ysize()
				nz = v1.get_ysize()
				if options.ou < 1 : radius = nx//2-2
				else:  radius = options.ou
				mask = model_circle(radius, nx, ny, nz)
			else:
				mask = get_im(args[2])
				
			for j in xrange(i+1, lend):
				lccc[mono(i,j)] = [ccc(v1, get_im( stack, j ), mask), 0]


		order = tsp(lccc)
		if(len(order) != lend):
			print  " problem with data length"
			from sys import exit
			exit()
		print  "Total sum of cccs :",TotalDistance(order, lccc)
		print "ordering :",order
		for i in xrange(lend):  get_im(stack, order[i]).write_image( new_stack, i )

	elif options.align:
		nargs = len(args)
		if nargs != 3:
			print "must provide name of input and two output files!"
			return

		from utilities import get_params2D, model_circle
		from fundamentals import rot_shift2D
		from statistics import ccc
		from time import time
		from alignment import align2d, align2d_scf
		from multi_shc import mult_transform 
		
		stack = args[0]
		new_stack = args[1]
		
		d = EMData.read_images(stack)

		"""
		# will align anyway
		try:
			ttt = d[0].get_attr('xform.params2d')
			for i in xrange(len(d)):
				alpha, sx, sy, mirror, scale = get_params2D(d[i])
				d[i] = rot_shift2D(d[i], alpha, sx, sy, mirror)
		except:
			pass
		"""

		nx = d[0].get_xsize()
		ny = d[0].get_ysize()
		if options.ou < 1 : radius = nx//2-2
		else:  radius = options.ou
		mask = model_circle(radius, nx, ny)

		if(options.xr < 0):	xrng = 0
		else:				xrng = options.xr
		if(options.yr < 0):	yrng = xrng
		else:				yrng = options.yr
			
		initial = max(options.initial, 0)

		from statistics import mono
		lend = len(d)
		lccc = [None]*(lend*(lend-1)/2)
		from utilities import read_text_row

		if  options.pairwiseccc == None or not os.path.exists(options.pairwiseccc) :
			st = time()
			for i in xrange(lend-1):
				for j in xrange(i+1, lend):
					#  j>i meaning mono entry (i,j) or (j,i) indicates T i->j (from smaller index to larger)
					#alpha, sx, sy, mir, peak = align2d(d[i],d[j], xrng, yrng, step=options.ts, first_ring=options.ir, last_ring=radius, mode = "F")
					alpha, sx, sy, mir, peak = align2d_scf(d[i],d[j], xrng, yrng, ou=radius)
					lccc[mono(i,j)] = [ccc(d[j], rot_shift2D(d[i], alpha, sx, sy, mir, 1.0), mask), alpha, sx, sy, mir]
				#print "  %4d   %10.1f"%(i,time()-st)

			if(not os.path.exists(options.pairwiseccc)):
				from utilities import write_text_row
				write_text_row([[initial,0,0,0,0]]+lccc,options.pairwiseccc)
		elif(os.path.exists(options.pairwiseccc)):
			lccc = read_text_row(options.pairwiseccc)
			initial = int(lccc[0][0] + 0.1)
			del lccc[0]


		for i in xrange(len(lccc)):
			T = Transform({"type":"2D","alpha":lccc[i][1],"tx":lccc[i][2],"ty":lccc[i][3],"mirror":int(lccc[i][4]+0.1)})
			lccc[i] = [lccc[i][0],T]

		tdummy = Transform({"type":"2D"})
		maxsum = -1.023
		for m in xrange(0,lend):#initial, initial+1):
			indc = range( lend )
			lsnake = [[m, tdummy, 0.0]]
			del indc[m]

			lsum = 0.0
			while len(indc) > 1:
				maxcit = -111.
				for i in xrange(len(indc)):
						cuc = lccc[mono(indc[i], lsnake[-1][0])][0]
						if cuc > maxcit:
								maxcit = cuc
								qi = indc[i]
								#  Here we need transformation from the current to the previous,
								#     meaning indc[i] -> lsnake[-1][0]
								T = lccc[mono(indc[i], lsnake[-1][0])][1]
								#  If direction is from larger to smaller index, the transformation has to be inverted
								if( indc[i] > lsnake[-1][0] ):  T = T.inverse()
								
								
				lsnake.append([qi,T, maxcit])
				lsum += maxcit

				del indc[indc.index(qi)]

			T = lccc[mono(indc[-1], lsnake[-1][0])][1]
			if( indc[-1] > lsnake[-1][0]):  T = T.inverse()
			lsnake.append([indc[-1], T, lccc[mono(indc[-1], lsnake[-1][0])][0]])
			print  " initial image and lsum  ",m,lsum
			#print lsnake
			if(lsum > maxsum):
				maxsum = lsum
				init = m
				snake = [lsnake[i] for i in xrange(lend)]
		print  "  Initial image selected : ",init,maxsum,"    ",TotalDistance([snake[m][0] for m in xrange(lend)], lccc)
		#for q in snake: print q

		from copy import deepcopy
		trans=deepcopy([snake[i][1] for i in xrange(len(snake))])
		print  [snake[i][0] for i in xrange(len(snake))]
		"""
		for m in xrange(lend):
			prms = trans[m].get_params("2D")
			print " %3d   %7.1f   %7.1f   %7.1f   %2d  %6.2f"%(snake[m][0], prms["alpha"], prms["tx"], prms["ty"], prms["mirror"], snake[m][2])
		"""
		for k in xrange(lend-2,0,-1):
			T = snake[k][1]
			for i in xrange(k+1, lend):
 					trans[i] = T*trans[i]
		#  To add - apply all transformations and do the overall centering.
		for m in xrange(lend):
			prms = trans[m].get_params("2D")
			#print " %3d   %7.1f   %7.1f   %7.1f   %2d  %6.2f"%(snake[m][0], prms["alpha"], prms["tx"], prms["ty"], prms["mirror"], snake[m][2])
			#rot_shift2D(d[snake[m][0]], prms["alpha"], prms["tx"], prms["ty"], prms["mirror"]).write_image(new_stack, m)
			rot_shift2D(d[snake[m][0]], prms["alpha"], 0.0,0.0, prms["mirror"]).write_image(new_stack, m)

		order = tsp(lccc)
		if(len(order) != lend):
			print  " problem with data length"
			from sys import exit
			exit()
		print  TotalDistance(order, lccc)
		print order
		ibeg = order.index(init)
		order = [order[(i+ibeg)%lend] for i in xrange(lend)]
		print  TotalDistance(order, lccc)
		print order


		snake = [tdummy]
		for i in xrange(1,lend):
			#  Here we need transformation from the current to the previous,
			#     meaning order[i] -> order[i-1]]
			T = lccc[mono(order[i], order[i-1])][1]
			#  If direction is from larger to smaller index, the transformation has to be inverted
			if( order[i] > order[i-1] ):  T = T.inverse()
			snake.append(T)
		assert(len(snake) == lend)
		from copy import deepcopy
		trans = deepcopy(snake)
		for k in xrange(lend-2,0,-1):
			T = snake[k]
			for i in xrange(k+1, lend):
 					trans[i] = T*trans[i]

		#  Try to smooth the angles - complicated, I am afraid one would have to use angles forward and backwards
		#     and find their average??
		#  In addition, one would have to recenter them
		"""
		trms = []
		for m in xrange(lend):
			prms = trans[m].get_params("2D")
			trms.append([prms["alpha"], prms["mirror"]])
		for i in xrange(3):
			for m in xrange(lend):
				mb = (m-1)%lend
				me = (m+1)%lend
				#  angles order mb,m,me
				# calculate predicted angles mb->m 
		"""

		for m in xrange(lend):
			prms = trans[m].get_params("2D")
			#rot_shift2D(d[order[m]], prms["alpha"], prms["tx"], prms["ty"], prms["mirror"]).write_image("metro.hdf", m)
			rot_shift2D(d[order[m]], prms["alpha"], 0.0,0.0, prms["mirror"]).write_image(args[2], m)

		"""
		#  This was an effort to get number of loops, inconclusive, to say the least
		from numpy import outer, zeros, float32, sqrt
		lend = len(d)
 		cor = zeros(lend,float32)
 		cor = outer(cor, cor)
		for i in xrange(lend):  cor[i][i] = 1.0
		for i in xrange(lend-1):
			for j in xrange(i+1, lend):
				cor[i,j] = lccc[mono(i,j)][0]
				cor[j,i] = cor[i,j]

		lmbd, eigvec = pca(cor)

		from utilities import write_text_file

		nvec=20
		print  [lmbd[j] for j in xrange(nvec)]
		print  " G"
		mm = [-1]*lend
		for i in xrange(lend):  # row
			mi = -1.0e23
			for j in xrange(nvec):
				qt = eigvec[j][i]
				if(abs(qt)>mi):
					mi = abs(qt)
					mm[i] = j
			for j in xrange(nvec):
				qt = eigvec[j][i]
				print  round(qt,3),   #  eigenvector
			print  mm[i]
		print
		for j in xrange(nvec):
			qt = []
			for i in xrange(lend):
				if(mm[i] == j):  qt.append(i)
			if(len(qt)>0):  write_text_file(qt,"loop%02d.txt"%j)
		"""
		"""
		print  [lmbd[j] for j in xrange(nvec)]
		print  " B"
		mm = [-1]*lend
		for i in xrange(lend):  # row
			mi = -1.0e23
			for j in xrange(nvec):
				qt = eigvec[j][i]/sqrt(lmbd[j])
				if(abs(qt)>mi):
					mi = abs(qt)
					mm[i] = j
			for j in xrange(nvec):
				qt = eigvec[j][i]/sqrt(lmbd[j])
				print  round(qt,3),   #  eigenvector
			print  mm[i]
		print
		"""

		"""
		lend=3
 		cor = zeros(lend,float32)
 		
 		cor = outer(cor, cor)
 		
 		
 		cor[0][0] =136.77
 		cor[0][1] = 79.15
 		cor[0][2] = 37.13
 		
 		cor[1][0] = 79.15
 		cor[2][0] = 37.13
 		
 		
 		cor[1][1] = 50.04
 		cor[1][2] = 21.65
 		
 		cor[2][1] = 21.65
 		
 		
 		cor[2][2] = 13.26

		lmbd, eigvec = pca(cor)
		print  lmbd
		print  eigvec
		for i in xrange(lend):  # row
			for j in xrange(lend):  print  eigvec[j][i],   #  eigenvector
			print
		print  " B"
		for i in xrange(lend):  # row
			for j in xrange(lend):  print  eigvec[j][i]/sqrt(lmbd[j]),   #  eigenvector
			print
		print  " G"
		for i in xrange(lend):  # row
			for j in xrange(lend):  print  eigvec[j][i]*sqrt(lmbd[j]),   #  eigenvector
			print
		"""
	else:
		nargs = len(args)
		if nargs != 2:
			print "must provide name of input and output file!"
			return
		
		from utilities import get_params2D, model_circle
		from fundamentals import rot_shift2D
		from statistics import ccc
		from time import time
		from alignment import align2d
		from multi_shc import mult_transform 
		
		stack = args[0]
		new_stack = args[1]
		
		d = EMData.read_images(stack)
		try:
			print "Using 2D alignment parameters from header.
			ttt = d[0].get_attr('xform.params2d')
			for i in xrange(len(d)):
				alpha, sx, sy, mirror, scale = get_params2D(d[i])
				d[i] = rot_shift2D(d[i], alpha, sx, sy, mirror)
		except:
			pass

		nx = d[0].get_xsize()
		ny = d[0].get_ysize()
		if options.radius < 1 : radius = nx//2-2
		else:  radius = options.radius
		mask = model_circle(radius, nx, ny)

		init = options.initial
		
		if init > -1 :
			print "      initial image: %d" % init
			temp = d[init].copy()
			temp.write_image(new_stack, 0)
			del d[init]
			k = 1
			lsum = 0.0
			while len(d) > 1:
				maxcit = -111.
				for i in xrange(len(d)):
						cuc = ccc(d[i], temp, mask)
						if cuc > maxcit:
								maxcit = cuc
								qi = i
				# 	print k, maxcit
				lsum += maxcit
				temp = d[qi].copy()
				del d[qi]
				temp.write_image(new_stack, k)
				k += 1
			print  lsum
			d[0].write_image(new_stack, k)
		else:			
			if options.circular :
				print "Using options.circular, no alignment"
				#  figure the "best circular" starting image
				maxsum = -1.023
				for m in xrange(len(d)):
					indc = range(len(d) )
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
						for i in xrange(len(indc)):
								cuc = ccc(d[indc[i]], temp, mask)
								if cuc > maxcit:
										maxcit = cuc
										qi = indc[i]
						lsnake[k] = qi
						lsum += maxcit
						del indc[indc.index(qi)]
						direction = -direction
						for i in xrange( 1,len(d) ):
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
						snake = [lsnake[i] for i in xrange(len(d))]
				print  "  Initial image selected : ",init,maxsum
				print lsnake
				for m in xrange(len(d)):  d[snake[m]].write_image(new_stack, m)
			else:
				#  figure the "best" starting image
				print "Straight chain, no alignment"
				maxsum = -1.023
				for m in xrange(len(d)):
					indc = range(len(d) )
					lsnake = [m]
					del indc[m]
					temp = d[m].copy()
					lsum = 0.0
					while len(indc) > 1:
						maxcit = -111.
						for i in xrange(len(indc)):
								cuc = ccc(d[indc[i]], temp, mask)
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
						snake = [lsnake[i] for i in xrange(len(d))]
				print  "  Initial image selected : ",init,maxsum
				print lsnake
				for m in xrange(len(d)):  d[snake[m]].write_image(new_stack, m)
		

if __name__ == "__main__":
	main()
