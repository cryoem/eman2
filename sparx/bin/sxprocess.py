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
 rotate_shift_params(paramsin, transf) has been moved to utilities
"""

from utilities import rotate_shift_params

"""
	Traveling salesman problem solved using Simulated Annealing.
"""
#from scipy import *
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

		print "T=%10.5f , distance= %10.5f , accepted steps= %d" %(T, dist, accepted)
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

	Generic 2-D image processing programs.

	Functionality:

	1.  Phase flip a stack of images and write output to new file:
		sxprocess.py input_stack.hdf output_stack.hdf --phase_flip

	2.  Resample (decimate or interpolate up) images (2D or 3D) in a stack to change the pixel size.
	    The window size will change accordingly.
		sxprocess input.hdf output.hdf  --changesize --ratio=0.5

	3.  Compute average power spectrum of a stack of 2D images with optional mask and/or padding (option wn) with zeroes or a 3-D volume.
		sxprocess.py input_stack.hdf powerspectrum.hdf  <mask.hdf> --pw [--wn=1024]

	4.  Generate a stack of projections bdb:data and micrographs with prefix mic (i.e., mic0.hdf, mic1.hdf etc) from structure input_structure.hdf, with CTF applied to both projections and micrographs:
		sxprocess.py input_structure.hdf data mic --generate_projections format="bdb":apix=5.2:CTF=True:boxsize=64

    5.  Retrieve original image numbers in the selected ISAC group (here group 12 from generation 3):
    	sxprocess.py  bdb:test3 class_averages_generation_3.hdf  list3_12.txt --isacgroup=12 --params=originalid

    6.  Retrieve original image numbers of images listed in ISAC output stack of averages:
    	sxprocess.py  select1.hdf  ohk.txt

    7.  Adjust rotationally averaged power spectrum of an image to that of a reference image or a reference 1D power spectrum stored in an ASCII file.
    	Optionally use a tangent low-pass filter.  Also works for a stack of images, in which case the output is also a stack.
    	sxprocess.py  vol.hdf ref.hdf  avol.hdf < 0.25 0.2> --adjpw
   	 	sxprocess.py  vol.hdf pw.txt   avol.hdf < 0.25 0.2> --adjpw

    8.  Generate a 1D rotationally averaged power spectrum of an image.
		sxprocess.py  vol.hdf --rotwp=rotpw.txt
    	# Output will contain three columns:
       (1) rotationally averaged power spectrum
       (2) logarithm of the rotationally averaged power spectrum
       (3) integer line number (from zero to approximately to half the image size)

    9.  Apply 3D transformation (rotation and/or shift) to a set of orientation parameters associated with projection data.
    	sxprocess.py  --transfromparams=phi,theta,psi,tx,ty,tz      input.txt  output.txt
    	The output file is then imported and 3D transformed volume computed:
    	sxheader.py  bdb:p  --params=xform.projection  --import=output.txt
    	mpirun -np 2 sxrecons3d_n.py  bdb:p tvol.hdf --MPI
    	The reconstructed volume is in the position of the volume computed using the input.txt parameters and then
    	transformed with rot_shift3D(vol, phi,theta,psi,tx,ty,tz)

   10.  Import ctf parameters from the output of sxcter into windowed particle headers.
	    There are three possible input files formats:  (1) all particles are in one stack, (2 aor 3) particles are in stacks, each stack corresponds to a single micrograph.
	    In each case the particles should contain a name of the micrograph of origin stores using attribute name 'ptcl_source_image'.
        Normally this is done by e2boxer.py during windowing.
	    Particles whose defocus or astigmatism error exceed set thresholds will be skipped, otherwise, virtual stacks with the original way preceded by G will be created.
		sxprocess.py  --input=bdb:data  --importctf=outdir/partres  --defocuserror=10.0  --astigmatismerror=5.0
		#  Output will be a vritual stack bdb:Gdata
		sxprocess.py  --input="bdb:directory/stacks*"  --importctf=outdir/partres  --defocuserror=10.0  --astigmatismerror=5.0
		To concatenate output files:
		cd directory
		e2bdb.py . --makevstack=bdb:allparticles  --filt=G
		IMPORTANT:  Please do not move (or remove!) any input/intermediate EMAN2DB files as the information is linked between them.

   11. Scale 3D shifts.  The shifts in the input five columns text file with 3D orientation parameters will be DIVIDED by the scale factor
		sxprocess.py  orientationparams.txt  scaledparams.txt  scale=0.5

   12. Generate soft-edged 3D mask from input 3D volume automatically or using the user-provided threshold.
        Automatically compute the threshold to intially obtain the largest density cluster.
        sxprocess.py  vol3d.hdf  mask3d.hdf  --adaptive_mask  --nsigma=3.0  --ndilation=1  --kernel_size=9  --gauss_standard_dev=5
        
        Use the user-provided threshold to intially obtain the largest density cluster.
        sxprocess.py  vol3d.hdf  mask3d.hdf  --adaptive_mask --threshold=0.05  -ndilation=0  --kernel_size=9  --gauss_standard_dev=5

   13. Generate binary 3D mask from input 3D volume using the user-provided threshold.
        sxprocess.py  vol3d.hdf  mask3d.hdf  --binary_mask  --threshold=0.05  --ne=3  --nd==3

   14. Postprocess 3-D or 2-D images:
   
	for 3-D volumes: 
		a. calculate FSC with provided mask; 
		b. sum two volume; 
		c. apply mask
		d. apply MTF correction;
		e. adjust power spectrum by 2*FSC/(1+FSC);  
		f. estimate B-factor from 10 Angstrom to resolution; 
		g. apply negative B-factor to enhance the volume;
		h. low_pass filter the volume
		options are independent of each others. However, if there is only one input map, do not use --fsc_adj option. 	
		--low_pass_filter: =0.0, low_pass filter to resolution; =-1., no low_pass filter; =5.8 low_pass filter to 5.8 Angstrom; =.2 low_pass filter to 0.2  
		--B_enhance:       =-1, B-factor is not applied; =0, program estimates B-factor from options.B_start(usually set as 10 Angstrom)to the resolution determined by FSC 0.143; =128., program use the given value 128. to enhance map.
		--mtf:             =aa.txt, for those high resolution maps, mtf corrections would significantly enhance structural features.
		--fsc_adj:         fsc adjustment of power spectrum is inclined to increase the slope of power spectrum of the summed volume.
		--do_adaptive_mask =True when it is restored, the program adaptively creates surface mask file using summed two volumes. This takes a couple of minutes. For map with dimension of 384*384*384, it takes 6 minutes.
		--output           output volume 
										
		sxprocess.py vol_0_unfil_026.hdf vol_1_unfil_026.hdf  --mask=mask15.mrc --postprocess   --pixel_size=1.2     --low_pass_filter =-1  --mtf=aa.txt  --fsc_adj --output=vol_post.hdf 
		sxprocess.py vol_0_unfil_026.hdf  --mask=mask15.mrc --postprocess   --pixel_size=1.2     --mtf=aa.txt        --output=vol_0_post.hdf
		sxprocess.py vol_0_unfil_026.hdf vol_1_unfil_026.hdf  --mask=mask15.mrc --postprocess   --pixel_size=1.2     --low_pass_filter=4.7  --mtf=aa.txt --fsc_adj
		sxprocess.py vol_0_unfil_026.hdf vol_1_unfil_026.hdf  --do_adaptive_mask   --postprocess   --pixel_size=1.2   --mtf=aa.txt --fsc_adj
		
	 for 2-D images:       calculate B-factor and apply negative B-factor to 2-D images.
		
   15. Window stack file -reduce the size of images without changing the pixel size.

   16. Create angular distribution .build file
        sxprocess.py --angular_distribution  inputfile=example/path/params.txt --pixel_size=1.0  --round_digit=5  --box_size=500  --particle_radius=175  --cylinder_width=1  --cylinder_length=10000
        

"""

	parser = OptionParser(usage,version=SPARXVERSION)
	parser.add_option("--order", 				action="store_true", help="Two arguments are required: name of input stack and desired name of output stack. The output stack is the input stack sorted by similarity in terms of cross-correlation coefficent.", default=False)
	parser.add_option("--order_lookup", 		action="store_true", help="Test/Debug.", default=False)
	parser.add_option("--order_metropolis", 	action="store_true", help="Test/Debug.", default=False)
	parser.add_option("--order_pca", 			action="store_true", help="Test/Debug.", default=False)
	parser.add_option("--initial",				type="int", 		default=-1, help="Specifies which image will be used as an initial seed to form the chain. (default = 0, means the first image)")
	parser.add_option("--circular", 			action="store_true", help="Select circular ordering (fisr image has to be similar to the last", default=False)
	parser.add_option("--radius", 				type="int", 		default=-1, help="Radius of a circular mask for similarity based ordering")
	parser.add_option("--changesize", 			action="store_true", help="resample (decimate or interpolate up) images (2D or 3D) in a stack to change the pixel size.", default=False)
	parser.add_option("--ratio", 				type="float", 		default=1.0, help="The ratio of new to old image size (if <1 the pixel size will increase and image size decrease, if>1, the other way round")
	parser.add_option("--pw", 					action="store_true", help="compute average power spectrum of a stack of 2-D images with optional padding (option wn) with zeroes", default=False)
	parser.add_option("--wn", 					type="int", 		default=-1, help="Size of window to use (should be larger/equal than particle box size, default padding to max(nx,ny))")
	parser.add_option("--phase_flip", 			action="store_true", help="Phase flip the input stack", default=False)
	parser.add_option("--makedb", 				metavar="param1=value1:param2=value2", type="string",
					action="append",  help="One argument is required: name of key with which the database will be created. Fill in database with parameters specified as follows: --makedb param1=value1:param2=value2, e.g. 'gauss_width'=1.0:'pixel_input'=5.2:'pixel_output'=5.2:'thr_low'=1.0")
	parser.add_option("--generate_projections", metavar="param1=value1:param2=value2", type="string",
					action="append", help="Three arguments are required: name of input structure from which to generate projections, desired name of output projection stack, and desired prefix for micrographs (e.g. if prefix is 'mic', then micrographs mic0.hdf, mic1.hdf etc will be generated). Optional arguments specifying format, apix, box size and whether to add CTF effects can be entered as follows after --generate_projections: format='bdb':apix=5.2:CTF=True:boxsize=100, or format='hdf', etc., where format is bdb or hdf, apix (pixel size) is a float, CTF is True or False, and boxsize denotes the dimension of the box (assumed to be a square). If an optional parameter is not specified, it will default as follows: format='bdb', apix=2.5, CTF=False, boxsize=64.")
	parser.add_option("--isacgroup", 			type="int", 		        help="Retrieve original image numbers in the selected ISAC group. See ISAC documentation for details.", default=-1)
	parser.add_option("--isacselect", 			action="store_true", 		help="Retrieve original image numbers of images listed in ISAC output stack of averages. See ISAC documentation for details.", default=False)
	parser.add_option("--params",	   			type="string",              default=None,    help="Name of header of parameter, which one depends on specific option")
	parser.add_option("--adjpw", 				action="store_true",	    help="Adjust rotationally averaged power spectrum of an image", default=False)
	parser.add_option("--rotpw", 				type="string",   	        default=None,    help="Name of the text file to contain rotationally averaged power spectrum of the input image.")
	parser.add_option("--transformparams",		type="string",   	        default=None,    help="Transform 3D projection orientation parameters using six 3D parameters (phi, theta,psi,sx,sy,sz).  Input: --transformparams=45.,66.,12.,-2,3,-5.5 desired six transformation of the reconstructed structure. Output: file with modified orientation parameters.")


	# import ctf estimates done using cter
	parser.add_option("--input",              	type="string",		default= None,     		  help="Input particles.")
	parser.add_option("--importctf",          	type="string",		default= None,     		  help="Name of the file containing CTF parameters produced by sxcter.")
	parser.add_option("--defocuserror",       	type="float",  		default=1000000.0,        help="Exclude micrographs whose relative defocus error as estimated by sxcter is larger than defocuserror percent.  The error is computed as (std dev defocus)/defocus*100%")
	parser.add_option("--astigmatismerror",   	type="float",  		default=360.0,            help="Set to zero astigmatism for micrographs whose astigmatism angular error as estimated by sxcter is larger than astigmatismerror degrees.")

	# import ctf estimates done using cter
	parser.add_option("--scale",              	type="float", 		default=-1.0,      		  help="Divide shifts in the input 3D orientation parameters text file by the scale factor.")

	# Generate soft-edged 3D mask from input 3D volume and Generate binarized version of input 3D volume
	parser.add_option("--adaptive_mask",        action="store_true",                      help="generate soft-edged 3D mask from input 3D volume", default= False)
	parser.add_option("--nsigma",               type="float",        default= 1.0,        help="number of times of sigma of the input volume to intially obtain the largest density cluster")
	parser.add_option("--threshold",            type="float",        default= -9999.0,    help="threshold provided by user to intially obtain the largest density cluster")
	parser.add_option("--ndilation",            type="int",          default= 3,          help="number of times of dilation applied to the largest cluster of density")
	parser.add_option("--kernel_size",          type="int",          default= 11,         help="convolution kernel for smoothing the edge of the mask")
	parser.add_option("--gauss_standard_dev",   type="int",          default= 9,          help="stanadard deviation value to generate Gaussian edge")
	
	# Generate soft-edged 3D mask from input 3D volume and Generate binarized version of input 3D volume
	parser.add_option("--binary_mask",          action="store_true",                      help="generate binary 3D mask from input 3D volume", default=False)
	parser.add_option("--bin_threshold",        type="float",        default= 0.0,        help="threshold provided by user to binarize input volume")
	parser.add_option("--ne",                   type="int",          default= 0,          help="number of times to erode binarized volume")
	parser.add_option("--nd",                   type="int",          default= 0,          help="number of times to dilate binarized volume")

	# Postprocess 3-D or 2-D images
	parser.add_option("--postprocess",          action="store_true",                      help="postprocess unfiltered odd, even 3-D volumes",default=False)
	parser.add_option("--mtf",                  type="string",        default= None,      help="mtf text file of camera")
	parser.add_option("--fsc_adj",              action="store_true",                      help="adjust the power spectrum of summed volume by their FSC", default=False)
	parser.add_option("--B_enhance",            type="float",         default=0.0,        help="apply Bfactor to enhance map or not")
	parser.add_option("--low_pass_filter",      type="float",         default=0.0,        help="=0.0, low_pass filter to resolution limit; =some value, low_pass filter to some valume; =-1, not low_pass filter ")
	parser.add_option("--aa",                   type="float",         default=.1,         help="low pass filter falloff" )
	parser.add_option("--mask",                 type="string",        help="input mask file",  default= None)
	parser.add_option("--output",               type="string",        help="output file name", default = "vol_postrefine_masked.hdf")
	parser.add_option("--pixel_size",           type="float",         help="pixel size of the data", default=0.0)
	parser.add_option("--B_start",              type="float",         help="starting frequency in Angstrom for B-factor estimation", default=10.)
	parser.add_option("--B_stop",               type="float",         help="cutoff frequency in Angstrom for B-factor estimation, cutoff is set to the frequency where fsc < 0.0", default=0.0)
	parser.add_option("--do_adaptive_mask",     action="store_true",  help="generate adaptive mask with the given threshold ", default= False)
	parser.add_option("--mask_threshold",       type="float",         help=" the threshold for adaptive_mask", default= 0.02)
	parser.add_option("--consine_edge", 	    type="float",		  help="the width for cosine transition area ", default= 6.0)
	parser.add_option("--dilation", 			type="float",		  help="the pixels for dilate or erosion of binary mask ", default= 3.0)
	#parser.add_option("--randomphasesafter", 	type="float",		  help=" set Fourier pixels random phases after FSC value ", default= 0.8)	
	# 
	parser.add_option("--window_stack",         action="store_true",                      help="window stack images using a smaller window size", default=False)
	parser.add_option("--box",                  type="int",		      default= 0,         help="the new window size ")
	
	# Options for angular distribution
	parser.add_option('--angular_distribution',    	action="store_true",  	default=False,        	help='create an angular distribution file based on a project3d.txt')
	parser.add_option('--round_digit',             	type='int',          	default=5,           	help='accuracy of the loaded angle (default 5)')
	parser.add_option('--box_size',                	type='int',          	default=500,         	help='box size in pixel used for calculating the center of the particle [px] (default 500)')
	parser.add_option('--particle_radius',     		type='int',          	default=175,         	help='particle radius [Pixels] (default 175)')
	parser.add_option('--cylinder_width',      		type='int',          	default=1,           	help='width of the cylinder (default 1)')
	parser.add_option('--cylinder_length',     		type='int',          	default=10000,       	help='length of the cylinder (default 10000)')
	(options, args) = parser.parse_args()

	global_def.BATCH = True

	if options.phase_flip:
		nargs = len(args)
		if nargs != 2:
			print "must provide name of input and output file!"
			return
		from EMAN2 import Processor
		instack = args[0]
		outstack = args[1]
		nima = EMUtil.get_image_count(instack)
		from filter import filt_ctf
		for i in xrange(nima):
			img = EMData()
			img.read_image(instack, i)
			try:
				ctf = img.get_attr('ctf')
			except:
				print "no ctf information in input stack! Exiting..."
				return

			dopad = True
			sign = 1
			binary = 1  # phase flip

			assert img.get_ysize() > 1
			dict = ctf.to_dict()
			dz = dict["defocus"]
			cs = dict["cs"]
			voltage = dict["voltage"]
			pixel_size = dict["apix"]
			b_factor = dict["bfactor"]
			ampcont = dict["ampcont"]
			dza = dict["dfdiff"]
			azz = dict["dfang"]

			if dopad and not img.is_complex(): ip = 1
			else:                             ip = 0


			params = {"filter_type": Processor.fourier_filter_types.CTF_,
	 			"defocus" : dz,
				"Cs": cs,
				"voltage": voltage,
				"Pixel_size": pixel_size,
				"B_factor": b_factor,
				"amp_contrast": ampcont,
				"dopad": ip,
				"binary": binary,
				"sign": sign,
				"dza": dza,
				"azz":azz}

			tmp = Processor.EMFourierFilter(img, params)
			tmp.set_attr_dict({"ctf": ctf})

			tmp.write_image(outstack, i)

	elif options.changesize:
		nargs = len(args)
		if nargs != 2:
			ERROR("must provide name of input and output file!", "change size", 1)
			return
		from utilities import get_im
		instack = args[0]
		outstack = args[1]
		sub_rate = float(options.ratio)

		nima = EMUtil.get_image_count(instack)
		from fundamentals import resample
		for i in xrange(nima):
			resample(get_im(instack, i), sub_rate).write_image(outstack, i)

	elif options.isacgroup>-1:
		nargs = len(args)
		if nargs != 3:
			ERROR("Three files needed on input!", "isacgroup", 1)
			return
		from utilities import get_im
		instack = args[0]
		m=get_im(args[1],int(options.isacgroup)).get_attr("members")
		l = []
		for k in m:
			l.append(int(get_im(args[0],k).get_attr(options.params)))
		from utilities import write_text_file
		write_text_file(l, args[2])

	elif options.isacselect:
		nargs = len(args)
		if nargs != 2:
			ERROR("Two files needed on input!", "isacgroup", 1)
			return
		from utilities import get_im
		nima = EMUtil.get_image_count(args[0])
		m = []
		for k in xrange(nima):
			m += get_im(args[0],k).get_attr("members")
		m.sort()
		from utilities import write_text_file
		write_text_file(m, args[1])

	elif options.pw:
		nargs = len(args)
		if nargs < 2:
			ERROR("must provide name of input and output file!", "pw", 1)
			return
		from utilities import get_im, write_text_file
		from fundamentals import rops_table
		d = get_im(args[0])
		ndim = d.get_ndim()
		if ndim ==3:
			pw = rops_table(d)
			write_text_file(pw, args[1])
		else:
			nx = d.get_xsize()
			ny = d.get_ysize()
			if nargs ==3: mask = get_im(args[2])
			wn = int(options.wn)
			if wn == -1:
				wn = max(nx, ny)
			else:
				if( (wn<nx) or (wn<ny) ):  ERROR("window size cannot be smaller than the image size","pw",1)
			n = EMUtil.get_image_count(args[0])
			from utilities import model_blank, model_circle, pad
			from EMAN2 import periodogram
			p = model_blank(wn,wn)

			for i in xrange(n):
				d = get_im(args[0], i)
				if nargs==3:
					d *=mask
				st = Util.infomask(d, None, True)
				d -= st[0]
				p += periodogram(pad(d, wn, wn, 1, 0.))
			p /= n
			p.write_image(args[1])

	elif options.adjpw:

		if len(args) < 3:
			ERROR("filt_by_rops input target output fl aa (the last two are optional parameters of a low-pass filter)","adjpw",1)
			return
		img_stack = args[0]
		from math         import sqrt
		from fundamentals import rops_table, fft
		from utilities    import read_text_file, get_im
		from filter       import  filt_tanl, filt_table
		if(  args[1][-3:] == 'txt'):
			rops_dst = read_text_file( args[1] )
		else:
			rops_dst = rops_table(get_im( args[1] ))

		out_stack = args[2]
		if(len(args) >4):
			fl = float(args[3])
			aa = float(args[4])
		else:
			fl = -1.0
			aa = 0.0

		nimage = EMUtil.get_image_count( img_stack )

		for i in xrange(nimage):
			img = fft(get_im(img_stack, i) )
			rops_src = rops_table(img)

			assert len(rops_dst) == len(rops_src)

			table = [0.0]*len(rops_dst)
			for j in xrange( len(rops_dst) ):
				table[j] = sqrt( rops_dst[j]/rops_src[j] )

			if( fl > 0.0):
				img = filt_tanl(img, fl, aa)
			img = fft(filt_table(img, table))
			img.write_image(out_stack, i)

	elif options.rotpw != None:

		if len(args) != 1:
			ERROR("Only one input permitted","rotpw",1)
			return
		from utilities import write_text_file, get_im
		from fundamentals import rops_table
		from math import log10
		t = rops_table(get_im(args[0]))
		x = range(len(t))
		r = [0.0]*len(x)
		for i in x:  r[i] = log10(t[i])
		write_text_file([t,r,x],options.rotpw)

	elif options.transformparams != None:
		if len(args) != 2:
			ERROR("Please provide names of input and output files with orientation parameters","transformparams",1)
			return
		from utilities import read_text_row, write_text_row
		transf = [0.0]*6
		spl=options.transformparams.split(',')
		for i in xrange(len(spl)):  transf[i] = float(spl[i])

		write_text_row( rotate_shift_params(read_text_row(args[0]), transf)	, args[1])

	elif options.makedb != None:
		nargs = len(args)
		if nargs != 1:
			print "must provide exactly one argument denoting database key under which the input params will be stored"
			return
		dbkey = args[0]
		print "database key under which params will be stored: ", dbkey
		gbdb = js_open_dict("e2boxercache/gauss_box_DB.json")

		parmstr = 'dummy:'+options.makedb[0]
		(processorname, param_dict) = parsemodopt(parmstr)
		dbdict = {}
		for pkey in param_dict:
			if (pkey == 'invert_contrast') or (pkey == 'use_variance'):
				if param_dict[pkey] == 'True':
					dbdict[pkey] = True
				else:
					dbdict[pkey] = False
			else:
				dbdict[pkey] = param_dict[pkey]
		gbdb[dbkey] = dbdict

	elif options.generate_projections:
		nargs = len(args)
		if nargs != 3:
			ERROR("Must provide name of input structure(s) from which to generate projections, name of output projection stack, and prefix for output micrographs."\
			"sxprocess - generate projections",1)
			return
		inpstr  = args[0]
		outstk  = args[1]
		micpref = args[2]

		parmstr = 'dummy:'+options.generate_projections[0]
		(processorname, param_dict) = parsemodopt(parmstr)

		parm_CTF    = False
		parm_format = 'bdb'
		parm_apix   = 2.5

		if 'CTF' in param_dict:
			if param_dict['CTF'] == 'True':
				parm_CTF = True

		if 'format' in param_dict:
			parm_format = param_dict['format']

		if 'apix' in param_dict:
			parm_apix = float(param_dict['apix'])

		boxsize = 64
		if 'boxsize' in param_dict:
			boxsize = int(param_dict['boxsize'])

		print "pixel size: ", parm_apix, " format: ", parm_format, " add CTF: ", parm_CTF, " box size: ", boxsize

		scale_mult      = 2500
		sigma_add       = 1.5
		sigma_proj      = 30.0
		sigma2_proj     = 17.5
		sigma_gauss     = 0.3
		sigma_mic       = 30.0
		sigma2_mic      = 17.5
		sigma_gauss_mic = 0.3

		if 'scale_mult' in param_dict:
			scale_mult = float(param_dict['scale_mult'])
		if 'sigma_add' in param_dict:
			sigma_add = float(param_dict['sigma_add'])
		if 'sigma_proj' in param_dict:
			sigma_proj = float(param_dict['sigma_proj'])
		if 'sigma2_proj' in param_dict:
			sigma2_proj = float(param_dict['sigma2_proj'])
		if 'sigma_gauss' in param_dict:
			sigma_gauss = float(param_dict['sigma_gauss'])
		if 'sigma_mic' in param_dict:
			sigma_mic = float(param_dict['sigma_mic'])
		if 'sigma2_mic' in param_dict:
			sigma2_mic = float(param_dict['sigma2_mic'])
		if 'sigma_gauss_mic' in param_dict:
			sigma_gauss_mic = float(param_dict['sigma_gauss_mic'])

		from filter import filt_gaussl, filt_ctf
		from utilities import drop_spider_doc, even_angles, model_gauss, delete_bdb, model_blank,pad,model_gauss_noise,set_params2D, set_params_proj
		from projection import prep_vol,prgs
		seed(14567)
		delta = 29
		angles = even_angles(delta, 0.0, 89.9, 0.0, 359.9, "S")
		nangle = len(angles)

		modelvol = []
		nvlms = EMUtil.get_image_count(inpstr)
		from utilities import get_im
		for k in xrange(nvlms):  modelvol.append(get_im(inpstr,k))

		nx = modelvol[0].get_xsize()

		if nx != boxsize:
			ERROR("Requested box dimension does not match dimension of the input model.", \
			"sxprocess - generate projections",1)
		nvol = 10
		volfts = [[] for k in xrange(nvlms)]
		for k in xrange(nvlms):
			for i in xrange(nvol):
				sigma = sigma_add + random()  # 1.5-2.5
				addon = model_gauss(sigma, boxsize, boxsize, boxsize, sigma, sigma, 38, 38, 40 )
				scale = scale_mult * (0.5+random())
				vf, kb = prep_vol(modelvol[k] + scale*addon)
				volfts[k].append(vf)
		del vf, modelvol

		if parm_format == "bdb":
			stack_data = "bdb:"+outstk
			delete_bdb(stack_data)
		else:
			stack_data = outstk + ".hdf"
		Cs      = 2.0
		pixel   = parm_apix
		voltage = 120.0
		ampcont = 10.0
		ibd     = 4096/2-boxsize
		iprj    = 0

		width = 240
		xstart = 8 + boxsize/2
		ystart = 8 + boxsize/2
		rowlen = 17
		from random import randint
		params = []
		for idef in xrange(3, 8):

			irow = 0
			icol = 0

			mic = model_blank(4096, 4096)
			defocus = idef * 0.5#0.2
			if parm_CTF:
				astampl=defocus*0.15
				astangl=50.0
				ctf = generate_ctf([defocus, Cs, voltage,  pixel, 0.0, ampcont, astampl, astangl])

			for i in xrange(nangle):
				for k in xrange(12):
					dphi = 8.0*(random()-0.5)
					dtht = 8.0*(random()-0.5)
					psi  = 360.0*random()

					phi = angles[i][0]+dphi
					tht = angles[i][1]+dtht

					s2x = 4.0*(random()-0.5)
					s2y = 4.0*(random()-0.5)

					params.append([phi, tht, psi, s2x, s2y])

					ivol = iprj % nvol
					#imgsrc = randint(0,nvlms-1)
					imgsrc = iprj % nvlms
					proj = prgs(volfts[imgsrc][ivol], kb, [phi, tht, psi, -s2x, -s2y])

					x = xstart + irow * width
					y = ystart + icol * width

					mic += pad(proj, 4096, 4096, 1, 0.0, x-2048, y-2048, 0)

					proj = proj + model_gauss_noise( sigma_proj, nx, nx )
					if parm_CTF:
						proj = filt_ctf(proj, ctf)
						proj.set_attr_dict({"ctf":ctf, "ctf_applied":0})

					proj = proj + filt_gaussl(model_gauss_noise(sigma2_proj, nx, nx), sigma_gauss)
					proj.set_attr("origimgsrc",imgsrc)
					proj.set_attr("test_id", iprj)
					proj.set_attr("ptcl_source_image",micpref + "%1d.hdf" % (idef-3))
					# flags describing the status of the image (1 = true, 0 = false)
					set_params2D(proj, [0.0, 0.0, 0.0, 0, 1.0])
					set_params_proj(proj, [phi, tht, psi, s2x, s2y])

					proj.write_image(stack_data, iprj)

					icol += 1
					if icol == rowlen:
						icol = 0
						irow += 1

					iprj += 1

			mic += model_gauss_noise(sigma_mic,4096,4096)
			if parm_CTF:
				#apply CTF
				mic = filt_ctf(mic, ctf)
			mic += filt_gaussl(model_gauss_noise(sigma2_mic, 4096, 4096), sigma_gauss_mic)

			mic.write_image(micpref + "%1d.hdf" % (idef-3), 0)

		drop_spider_doc("params.txt", params)

	elif options.importctf != None:
		print ' IMPORTCTF  '
		from utilities import read_text_row,write_text_row
		from random import randint
		import subprocess
		grpfile = 'groupid%04d'%randint(1000,9999)
		ctfpfile = 'ctfpfile%04d'%randint(1000,9999)
		cterr = [options.defocuserror/100.0, options.astigmatismerror]
		ctfs = read_text_row(options.importctf)
		for kk in xrange(len(ctfs)):
			root,name = os.path.split(ctfs[kk][-1])
			ctfs[kk][-1] = name[:-4]
		if(options.input[:4] != 'bdb:'):
			ERROR('Sorry, only bdb files implemented','importctf',1)
		d = options.input[4:]
		#try:     str = d.index('*')
		#except:  str = -1
		from string import split
		import glob
		uu = os.path.split(d)
		uu = os.path.join(uu[0],'EMAN2DB',uu[1]+'.bdb')
		flist = glob.glob(uu)
		for i in xrange(len(flist)):
			root,name = os.path.split(flist[i])
			root = root[:-7]
			name = name[:-4]
			fil = 'bdb:'+os.path.join(root,name)
			sourcemic = EMUtil.get_all_attributes(fil,'ptcl_source_image')
			nn = len(sourcemic)
			gctfp = []
			groupid = []
			for kk in xrange(nn):
				junk,name2 = os.path.split(sourcemic[kk])
				name2 = name2[:-4]
				ctfp = [-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
				for ll in xrange(len(ctfs)):
					if(name2 == ctfs[ll][-1]):
						#  found correct
						if(ctfs[ll][8]/ctfs[ll][0] <= cterr[0]):
							#  acceptable defocus error
							ctfp = ctfs[ll][:8]
							if(ctfs[ll][10] > cterr[1] ):
								# error of astigmatism exceed the threshold, set astigmatism to zero.
								ctfp[6] = 0.0
								ctfp[7] = 0.0
							gctfp.append(ctfp)
							groupid.append(kk)
						break
			if(len(groupid) > 0):
				write_text_row(groupid, grpfile)
				write_text_row(gctfp, ctfpfile)
				cmd = "{} {} {} {}".format('e2bdb.py',fil,'--makevstack=bdb:'+root+'G'+name,'--list='+grpfile)
				#print cmd
				subprocess.call(cmd, shell=True)
				cmd = "{} {} {} {}".format('sxheader.py','bdb:'+root+'G'+name,'--params=ctf','--import='+ctfpfile)
				#print cmd
				subprocess.call(cmd, shell=True)
			else:
				print  ' >>>  Group ',name,'  skipped.'

		cmd = "{} {} {}".format("rm -f",grpfile,ctfpfile)
		subprocess.call(cmd, shell=True)

	elif options.scale > 0.0:
		from utilities import read_text_row,write_text_row
		scale = options.scale
		nargs = len(args)
		if nargs != 2:
			print "Please provide names of input and output file!"
			return
		p = read_text_row(args[0])
		for i in xrange(len(p)):
			p[i][3] /= scale
			p[i][4] /= scale
		write_text_row(p, args[1])

	elif options.adaptive_mask:
		from utilities import get_im
		from morphology import adaptive_mask1
		nargs = len(args)
		if nargs ==0:
			print " Generate soft-edged 3D mask from input 3D volume automatically or using the user provided threshold."
			return
		elif nargs > 2:
			ERROR( "Too many arguments are given, try again!", "options.adaptive_mask")
			return
		
		print "Started sxprocess.py  --adaptive_mask"
		inputvol = get_im(args[0]) # args[0]: input 3D volume file path
		input_path, input_file_name = os.path.split(args[0])
		input_file_name_root,ext=os.path.splitext(input_file_name)
		if nargs == 2:  mask_file_name = args[1] # args[1]: output 3D mask file path
		else:           mask_file_name = "adaptive_mask_for_" + input_file_name_root + ".hdf" # Only hdf file is output.
		mask3d, density_stats = adaptive_mask1(inputvol, options.nsigma, options.threshold, options.ndilation, options.kernel_size, options.gauss_standard_dev)
		mask3d.write_image(mask_file_name)
		print "  Applied threshold for binarize: %f" % density_stats[0]
		print "  Background density average    : %f" % density_stats[1]
		print "  Background density sigma      : %f" % density_stats[2]
		print "  Sigma factor (nsigma)         : %f" % density_stats[3]
		print "Finished sxprocess.py  --adaptive_mask"
	
	elif options.binary_mask:
		from utilities import get_im
		from morphology import binarize, erosion, dilation
		nargs = len(args)
		if nargs == 0:
			print " Generate binary 3D mask from input 3D volume using the user-provided threshold."
			return
		elif nargs > 2:
			print "Too many arguments are given, try again!"
			return
		
		print "Started sxprocess.py  --binary_mask"
		inputvol = get_im(args[0])
		input_path, input_file_name = os.path.split(args[0])
		input_file_name_root,ext=os.path.splitext(input_file_name)
		if nargs == 2:  mask_file_name = args[1]
		else:           mask_file_name = "binary_mask_for_" + input_file_name_root + ".hdf" # Only hdf file is output.
		mask3d = binarize(inputvol, options.bin_threshold)
		for i in xrange(options.ne): mask3d = erosion(mask3d)
		for i in xrange(options.nd): mask3d = dilation(mask3d)
		mask3d.write_image(mask_file_name)
		print "Applied threshold value for binarization is %f" % options.bin_threshold
		print "Finished sxprocess.py  --binary_mask"

	elif options.postprocess:
		from logger import Logger,BaseLogger_Files
		if os.path.exists("log.txt"): os.system(" rm log.txt")
		log_main=Logger(BaseLogger_Files())
		log_main.prefix="./"
		print_msg ="--------------------------------------------"
		#line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
		log_main.add(print_msg)
		print_msg=" Sphire postprocess"
		log_main.add(print_msg)
		from utilities    import get_im
		from fundamentals import rot_avg_table
		from morphology   import compute_bfactor,power
		from statistics   import fsc
		from filter       import filt_table, filt_gaussinv
		from EMAN2 import periodogram
		if len(args)<1 or len(args)>2:
			ERROR(" number of inputs is incorrection", " --postprocess option")
			exit()
		if options.pixel_size ==0:
			ERROR(" Set pixel_size value first! There is no default value for pixel_size", " --postprocess option")
			exit()
		try:
			e1   = get_im(args[0],0)
		except:
			ERROR(args[0]+" does not exist", " --postprocess option")
			exit()
		nx = e1.get_xsize()
		ny = e1.get_ysize()
		nz = e1.get_zsize()
		log_main.add(" The sphire postprocess command:  ")
		line=" "
		for a in sys.argv:
			line +=" "+a
		log_main.add(line)
		log_main.add(" ---------- Settings given by all options")
		log_main.add("pixle_size  "+str(options.pixel_size))
		log_main.add("mask        "+str(options.mask))
		log_main.add("fsc_adj     "+str(options.fsc_adj))
		log_main.add("B_enhance   "+str(options.B_enhance))
		log_main.add("low_pass_filter  "+str(options.low_pass_filter))
		log_main.add("B_start  "+str(options.B_start))
		log_main.add("B_stop   "+str(options.B_stop))
		log_main.add("mtf     "+str(options.mtf))
		log_main.add("output  "+str(options.output))
		log_main.add("do_adaptive_mask  "+str(options.do_adaptive_mask))
		log_main.add("cosine_edge    "+str(options.consine_edge))
		log_main.add("dilation    "+str(options.dilation))
		#log_main.add("randomphasesafter    "+str(options.randomphasesafter))
		log_main.add("-----------")
		if e1.get_zsize() == 1:  # 2D case
			log_main.add("2-D postprocess for ISAC averaged images")
			nimage = EMUtil.get_image_count(args[0])
			if options.mask !=None:
				try:
					m = get_im(options.mask)
					log_main.add("user provided mask is %s"%options.mask)
				except:
					ERROR(" mask image %s does not exists"%options.mask, " --postprocess for 2-D")
					exit()
			else:
				m = None
				log_main.add("mask is not used")
			log_main.add("total number of average images is %d"%nimage)
			for i in xrange(nimage):
				e1 = get_im(args[0],i)
				if m: e1 *=m
				if options.B_enhance ==0.0 or options.B_enhance == -1.:
					guinerline = rot_avg_table(power(periodogram(e1),.5))
					if options.B_stop:
						freq_max   =  1./(2.*options.pixel_size)
					else:
						freq_max =1./options.B_stop
					freq_min   =  1./options.B_start
					log_main.add("B-factor exp(-B*s^2) is estimated from %f Angstrom to %f Angstrom"%(options.B_start, 2*options.pixel_size))
					b,junk,ifreqmin, ifreqmax =compute_bfactor(guinerline, freq_min, freq_max, options.pixel_size)
					global_b = b*4
					log_main.add( "the estimated slope of rotationally averaged Fourier factors  of the summed volumes is %f"%round(-b,2))
				else:
					global_b = option.B_enhance
					log_main.add( "User provided B_factor is %f"%global_b)
				sigma_of_inverse=sqrt(2./global_b)
				e1 = filt_gaussinv(e1,sigma_of_inverse)
				if options.low_pass_filter>0.0 and options.low_pass_filte<0.5:
					log_main.add("low-pass filter ff %   aa  %f"%(options.low_pass_filter, options.aa))
					from filter import filt_tanl
					e1 =filt_tanl(e1,options.low_pass_filter, options.aa)
				elif options.low_pass_filte>0.5:
					e1 =filt_tanl(e1,options.low_pass_filter/option.pixel_size, options.aa)
				e1.write_image(options.output)
		else:   # 3D case
			log_main.add( "3-D refinement postprocess ")
			nargs     = len(args)
			if nargs >=3:
				ERROR(" Too many inputs!", "--postprocess option for 3-D")
			log_main.add("The first input volume: %s"%args[0])
			try: 
				map1    = get_im(args[0])
			except:
				ERROR(" Sphire postprocess fails to read the first map "+args[0], "--postprocess option for 3-D")
				exit()
			if nargs >1:
				log_main.add("The second input volume: %s"%args[1])
				try:
					map2  = get_im(args[1])
				except:
					ERROR(" fail to read the second map "+args[1], "--postprocess option for 3-D")
					exit()
			if (map2.get_xsize() != map1.get_xsize()) or (map2.get_ysize() != map1.get_ysize()) or (map2.get_zsize() != map1.get_zsize()):
				ERROR(" Two input maps have different image size", "--postprocess option for 3-D")
			## prepare mask 
			if options.mask != None and options.do_adaptive_mask:
				ERROR("Wrong options, use either adaptive_mask or user provided mask", " options.mask and options.do_adaptive_mask ")
			if options.mask != None:
				log_main.add("User provided mask: %s"%options.mask)
				try:
					m = get_im(options.mask)
				except:
					ERROR("Sphire postprocess fails to read mask file "+options.mask, "--postprocess option for 3-D")
					exit()
				if (m.get_xsize() != map1.get_xsize()) or (m.get_ysize() != map1.get_ysize()) or (m.get_zsize() != map1.get_zsize()):
					ERROR(" mask file  "+options.mask+" has different size with input image  ", "--postprocess for mask "+options.mask)
			elif options.do_adaptive_mask:
				if nargs >1 :
					map1 +=map2
					map1 /=2.
				m = Util.surface_mask(map1, options.mask_threshold, options.dilation, options.consine_edge)
				m.write_image("vol_adaptive_mask.hdf")
				map1 = get_im(args[0]) # re-read map1
			else:
				m = None
				log_main.add(" No mask is applied ")
			## prepare FSC
			from math import sqrt
			resolution_FSC143   	= 0.5 # for single volume, this is the default resolution
			resolution_FSChalf  	= 0.5
			frc_without_mask    	= None
			frc_with_mask       	= None
			frc_with_random_phases 	= None 
			if nargs >1: 
				frc_without_mask 		= fsc(map1, map2,1)
				if m: 
					frc_with_mask     = fsc(map1*m, map2*m,1)
					"""
					# determine random_phase_cutoff_pixel
					random_phase_cutoff_pixel =-1
					for ifreq in xrange(1,len(frc_without_mask[1])):
						if frc_without_mask[1][ifreq] < options.randomphasesafter: 
							random_phase_cutoff_pixel = ifreq
							break
					print("random phase %f  %d"%(options.randomphasesafter,random_phase_cutoff_pixel))
					from fundamentals import fft
					if random_phase_cutoff_pixel >0:
						map1 = Util.randomizedphasesafter(fft(map1),60)
						map2 = Util.randomizedphasesafter(fft(map2),60)
						map1 = fft(map1)*m
						map2 = fft(map2)*m
						frc_random_phase_masked = fsc(map1, map2,1)
						from utilities import write_text_file
						write_text_file(frc_random_phase_masked, "fsc_random_phase.txt")
						write_text_file(frc_with_mask, "fsc_with_mask.txt")
						write_text_file(frc_without_mask, "fsc_without_mask.txt")
					else:
						ERROR(" negative random phase cutoff pixel ", "random phase ", 1)
						#outfrc = [ frc_with_mask[1]]
						for ifreq in xrange(len(frc_with_mask[0])):
							if ifreq==0:
								outfrc[1].append(frc_with_mask[0][ifreq])
							else:	
								outfrc[1].append(options.pixel_size/frc_with_mask[0][ifreq])
					# FSC_RH see Richard Handson's paper
					frc_RH =[frc_with_mask[0],len(frc_with_mask[0])*[None]]
					frc_RH =[frc_with_mask[0],frc_with_mask[1]]
					#print(" frc_RH  %f"%frc_RH[1][5])
					for ifreq in xrange(len(frc_with_mask[1])):
						if ifreq <= random_phase_cutoff_pixel+2: frc_RH[1][ifreq] = frc_with_mask[1][ifreq]
						else:   
							if frc_random_phase_masked[1][ifreq] > frc_with_mask[1][ifreq]: frc_RH[1][ifreq] =0.0
							else:     frc_RH[1][ifreq] = (frc_with_mask[1][ifreq]-frc_random_phase_masked[1][ifreq])/(1.-frc_random_phase_masked[1][ifreq])
					"""
					frc_RH =[frc_with_mask[0],frc_with_mask[1]]
				else:
					frc_RH = fsc_without_mask
				from utilities import write_text_file
				write_text_file(frc_RH[1], "fsc.txt")
				for ifreq in xrange(len(frc_RH[1])):
					if frc_RH[1][ifreq] < 0.143:
						resolution_FSC143   = frc_RH[0][ifreq-1]
						break
				for ifreq in xrange(len(frc_RH[1])):
					if frc_RH[1][ifreq] < 0.5:
						resolution_FSChalf  = frc_RH[0][ifreq-1]
						break															
				from utilities import write_text_file
				#write_text_file(outfrc, "fsc_with_mask.txt")
			map1 = get_im(args[0])
			if nargs >1: 
				map2 = get_im(args[1])
				map1 +=map2
				map1 /=2.0
			outtext = [["Squaredfreqencies"],[ "LogOrignal"]]
			guinerline = rot_avg_table(power(periodogram(map1),.5))
			from math import log
			for ig in xrange(len(guinerline)):
				x = ig*.5/float(len(guinerline))/options.pixel_size
				outtext[0].append(x*x)
				outtext[1].append(log(guinerline[ig]))
			# starts adjustment of powerspectrum
			if options.mtf: # divided by the mtf   #1
				from fundamentals import fft
				log_main.add("MTF correction is applied")
				from utilities import read_text_file
				log_main.add("MTF file is %s"%options.mtf)
				try:
					mtf_core  = read_text_file(options.mtf, -1)
				except:
					ERROR(" Sphire postprocess fails to read MTF file "+options.mtf, "--postprocess option for 3-D")
					exit()
				map1 = fft(Util.divide_mtf(fft(map1), mtf_core[1], mtf_core[0]))
				outtext.append(["LogMTFdiv"])
				guinerline   = rot_avg_table(power(periodogram(map1),.5))
				for ig in xrange(len(guinerline)): outtext[-1].append(log(guinerline[ig]))
			if options.fsc_adj:    #2
				log_main.add("(2*FSC)/(1+FSC) is applied to adjust power spectrum ")
				if nargs==1:
					print("WARNING! there is only one input map,  and FSC adjustment cannot be done! Skip and continue...", "--postprocess  for 3-D")					
				else:
					#### FSC adjustment ((2.*fsc)/(1+fsc)) to the powerspectrum;
					fil = len(frc_RH[1])*[None]
					for i in xrange(len(fil)):
						if frc_RH[1][i]>=0.0: 	tmp = frc_RH[1][i]
						else: 					tmp = 0.0
						fil[i] = sqrt(2.*tmp/(1.+tmp))
					map1=filt_table(map1,fil)
					guinerline   = rot_avg_table(power(periodogram(map1),.5))
					outtext.append(["LogFSCadjusted"])
					for ig in xrange(len(guinerline)):
						outtext[-1].append(log(guinerline[ig]))
			
			if options.B_enhance !=-1:  #3
				if options.B_enhance == 0.0: # auto mode 
					cutoff_by_fsc = 0
					for ifreq in xrange(len(frc_RH[1])):
						if frc_RH[1][ifreq]<0.143:
							break
					cutoff_by_fsc = float(ifreq)
					#print(" cutoff_by_fsc ", cutoff_by_fsc)
					freq_max     =cutoff_by_fsc/(2.*len(frc_RH[0]))/options.pixel_size
					guinerline   = rot_avg_table(power(periodogram(map1),.5))
					logguinerline = []
					for ig in xrange(len(guinerline)):
						logguinerline.append(log(guinerline[ig]))
					freq_min     = 1./options.B_start # given frequencies with unit of Angstrom, say 10 Angstrom, 15  Angstrom
					if options.B_stop!=0.0: freq_max     = 1./options.B_stop 
					if freq_min>=freq_max:
						log_main.add("Your B_start is too high! Decrease it and rerun the program!")
						ERROR("your B_start is too high! Decrease it and re-run the program!", "--postprocess option")
						exit()
					from utilities import write_text_file
					#log_main.add("B-factor exp(-B*s^2) is estimated from %f Angstrom to %f Angstrom"%(round(1./freq_min,2), round(1./freq_max,2)))
					b,junk , ifreqmin, ifreqmax  =  compute_bfactor(guinerline, freq_min, freq_max, options.pixel_size)
					#log_main.add("The used pixels are from %d to %d"%(ifreqmin, ifreqmax))
					global_b     =  4.*b
					from statistics import pearson
					cc =pearson(junk[1],logguinerline)
					log_main.add("Similiarity between the fitted line and 1-D rotationally average power spectrum within [%d, %d] is %f"% \
					                                                  (ifreqmin, ifreqmax, pearson(junk[1][ifreqmin:ifreqmax],logguinerline[ifreqmin:ifreqmax])))
					log_main.add("The slope is %f Angstrom^2 "%(round(-b,2)))
					sigma_of_inverse = sqrt(2./(global_b/options.pixel_size**2))

				else: # User provided value
					#log_main.add( " apply user provided B-factor to enhance map!")
					log_main.add("User provided B-factor is %f Angstrom^2   "%options.B_enhance)
					sigma_of_inverse = sqrt(2./((abs(options.B_enhance))/options.pixel_size**2))
					global_b = options.B_enhance
				map1  = filt_gaussinv(map1,sigma_of_inverse)
				guinerline   = rot_avg_table(power(periodogram(map1),.5))
				outtext.append([" LogBfactorsharpened"])
				last_non_zero = -999.0
				for ig in xrange(len(guinerline)):
					if guinerline[ig]>0: 
						outtext[-1].append(log(guinerline[ig]))
						last_non_zero = log(guinerline[ig])
					else:
						outtext[-1].append(last_non_zero)			
			if options.low_pass_filter !=-1.: # User provided low-pass filter #4.
				from filter       import filt_tanl
				if options.low_pass_filter>0.5: # Input is in Angstrom 
					map1 =filt_tanl(map1,options.pixel_size/options.low_pass_filter, min(options.aa,.1))
					cutoff = options.low_pass_filter
				elif options.low_pass_filter>0.0 and options.low_pass_filter<0.5:  # input is in absolution frequency
					map1 =filt_tanl(map1,options.low_pass_filter, min(options.aa,.1))
					cutoff = options.pixel_size/options.low_pass_filter
				else: # low-pass filter to resolution determined by FSC0.143
					map1 = filt_tanl(map1,resolution_FSC143, options.aa)
					cutoff = options.pixel_size/resolution_FSC143
			map1.write_image("vol_postprocess_nomask.hdf")
			if m: map1 *=m
			map1.write_image(options.output)
			log_main.add("------ Summary -------")
			log_main.add("Resolution at criteria 0.143 is %f Angstrom"%round((options.pixel_size/resolution_FSC143),3))
			log_main.add("Resolution at criteria 0.5   is %f Angstrom"%round((options.pixel_size/resolution_FSChalf),3))
			if options.B_enhance !=-1:  log_main.add( " B-factor is  %f Angstrom^2  "%(round((-global_b),2)))
			else:                       log_main.add( " B-factor is not applied  ")
			log_main.add( "FSC curve is saved in fsc.txt  ")
			log_main.add( "Final processed volume is "+options.output)
			log_main.add("guinierlines in logscale are saved in guinierlines.txt")
			if options.low_pass_filter !=-1:  	log_main.add(" Top hat low-pass filter is applied to cut off high frequencies from resolution 1./%f Angstrom" %round(cutoff,2))
			else: 								log_main.add(" The final volume is not low_pass filtered. ")
			write_text_file(outtext, "guinierlines.txt")
				
	elif options.window_stack:
		nargs = len(args)
		if nargs ==0:
			print "  reduce image size of a stack"
			return
		else:
			output_stack_name = None
			inputstack = args[0]
			if nargs ==2:output_stack_name = args[1]
			input_path,input_file_name     = os.path.split(inputstack)
			input_file_name_root,ext       = os.path.splitext(input_file_name)
			if input_file_name_root[0:3]=="bdb":stack_is_bdb = True
			else:                               stack_is_bdb = False
			if output_stack_name is None:
				if stack_is_bdb: output_stack_name  = "bdb:reduced_"+input_file_name_root[4:]
				else: output_stack_name = "reduced_"+input_file_name_root+".hdf" # Only hdf file is output.
			nimage = EMUtil.get_image_count(inputstack)
			from fundamentals import window2d
			from utilities import get_im
			for i in xrange(nimage): window2d(get_im(inputstack,i),options.box,options.box).write_image(output_stack_name,i)

	elif options.angular_distribution:
		from utilities import angular_distribution
		nargs = len(args)
		if nargs > 1:
			print 'Too many inputs are given, see usage and restart the program!'
		else:
			if not os.path.exists(args[0]):
				ERROR(
					'Params file does not exists! Please rename and restart the program.', 1
					)
			strInput = args[0]
			strOutput = strInput[:-len(strInput.split('/')[-1])] + 'distribution.bild'
			angular_distribution(inputfile=strInput, options=options, output=strOutput)
	else:  ERROR("Please provide option name","sxprocess.py",1)

if __name__ == "__main__":
	main()
