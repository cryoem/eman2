#
from __future__ import print_function
# Author: Markus Stabrin 2019 (markus.stabrin@mpi-dortmund.mpg.de)
# Author: Fabian Schoenfeld 2019 (fabian.schoenfeld@mpi-dortmund.mpg.de)
# Author: Thorsten Wagner 2019 (thorsten.wagner@mpi-dortmund.mpg.de)
# Author: Tapu Shaikh 2019 (tapu.shaikh@mpi-dortmund.mpg.de)
# Author: Adnan Ali 2019 (adnan.ali@mpi-dortmund.mpg.de)
# Author: Luca Lusnig 2019 (luca.lusnig@mpi-dortmund.mpg.de)
# Author: Toshio Moriya 2019 (toshio.moriya@kek.jp)
# Author: Pawel A.Penczek, 09/09/2006 (Pawel.A.Penczek@uth.tmc.edu)
#
# Copyright (c) 2019 Max Planck Institute of Molecular Physiology
# Copyright (c) 2000-2006 The University of Texas - Houston Medical School
#
# This software is issued under a joint BSD/GNU license. You may use the
# source code in this file under either license. However, note that the
# complete EMAN2 and SPARX software packages have some GPL dependencies,
# so you are responsible for compliance with the licenses of these packages
# if you opt to use BSD licensing. The warranty disclaimer below holfds
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
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
#

# Comment by Zhengfan Yang on 06/11/10
# I decided to move everything related to pixel error to this file. Otherwise, they are
# scattered all over the places and there are a lot of duplications and confusions.
#
# This file contains the following functions:
#  1. pixel_error_2D (originally max_pixel_error, also I changed its input format)
#  2. max_3D_pixel_error
#  3. angle_diff
#  4. align_diff_params
#  5. align_diff
#  6. align_diff_textfile (new)
#  7. ave_ali_err_params
#  8. ave_ali_err
#  9. ave_ali_err_textfile (new)
# 10. multi_align_diff_params (new)
# 11. calc_connect_list (new)
# 12. ali_stable_list (new)
#
# Update on 09/01/10, we have decided for multiple alignment function 10 and 11 are not best way
# to determine the stability. We have instead written a new function to to this.
# 13. multi_align_stability (latest)
#
# Here, align_diff_params() and ave_ali_err_params() takes two lists of alignmnet parameters
# in the following format:
# [alpha1, sx1, sy1, mirror1, alpha2, sx2, sy2, mirror2, ...]
# align_diff() and ave_ali_err() takes two lists of images
# align_diff_textfile() and ave_ali_err_textfile() takes two textfiles, which are usually the output
# file of header() or sxheader.py.
# Per previous discussion, I decided not to support two stacks of images.

def pixel_error_2D(ali_params1, ali_params2, r = 1.0):
	"""
	Compute average squared 2D pixel error
	"""
	from math import radians, sin, pi, sqrt
	return (sin(radians(ali_params1[0]-ali_params2[0])/2)*(2*r+1))**2 / 2 + (ali_params1[1]-ali_params2[1])**2 + (ali_params1[2]-ali_params2[2])**2


def max_3D_pixel_error(t1, t2, r=1.0):
	"""
	Compute maximum pixel error between two sets of orientation parameters
	assuming object has radius r, t1 is the projection transformation
	of the first projection and t2 of the second one, respectively:
		t = Transform({"type":"spider","phi":phi,"theta":theta,"psi":psi})
		t.set_trans(Vec2f(-tx, -ty))
	Note the function is symmetric in t1, t2.
	"""
	from EMAN2 import Vec2f
	import types
	dummy = Transform()
	if( type(dummy) != type(t1)):
		t = Transform({"type":"spider","phi":t1[0],"theta":t1[1],"psi":t1[2]})
		t.set_trans(Vec2f(-t1[3], -t1[4]))
	else: t = t1
	if( type(dummy) != type(t2)):
		u = Transform({"type":"spider","phi":t2[0],"theta":t2[1],"psi":t2[2]})
		u.set_trans(Vec2f(-t2[3], -t2[4]))
	else: u = t2

	return EMData().max_3D_pixel_error(t, u, r)


def angle_ave(angle1):
	'''
	This function computes average angle of a set of angles.
	It also computes a measure of dispersion (incorrect).
	'''
	from math import cos, sin, pi, atan2, degrees, radians, sqrt

	nima = len(angle1)

	cosi = 0.0
	sini = 0.0
	for i in range(nima):
		qt = radians( angle1[i] )
		cosi += cos(qt)
		sini += sin(qt)
	alphai = degrees(atan2(sini, cosi))%360.0
	# what follows is not correct, it is just to give a measure of dispersion
	stdv = 0.0
	for i in range(nima):
		qt = angle1[i] - alphai
		if   qt >  180.0:   qt -= 360.
		elif qt < -180.0:   qt += 360.
		stdv += qt*qt
	stdv = sqrt(stdv/nima)

	return alphai, stdv


def angle_diff(angle1, angle2):
	'''
	This function determines the relative angle between two sets of angles.
	The resulting angle has to be added (modulo 360) to the first set.
	'''
	from math import cos, sin, pi, atan2, degrees, radians
	
	nima  = len(angle1)
	nima2 = len(angle2)
	if nima2 != nima:
		ERROR("Error: List lengths do not agree!","angle_diff",1)
	else:
		del nima2

	cosi = 0.0
	sini = 0.0
	for i in range(nima):
		qt = radians(angle2[i]-angle1[i])
		cosi += cos( qt )
		sini += sin( qt )
	alphai = degrees(atan2(sini, cosi))%360.0

	return alphai

def angle_diff_sym(angle1, angle2, simi=1):
	'''
	This function determines the relative angle around Z axis (phi) between two sets of angles
	   taking into account point group symmetry with multiplicity simi.
	The input has to be in the form [[phi0,theta0], [phi1,theta1], ...]
	  Only sets that have theta in the same range (0,90), or (90,180) are included in calculation.
	The resulting angle has to be added (modulo 360/simi) to the first set.
	'''
	from math import cos, sin, pi, atan2, degrees, radians
	
	nima  = len(angle1)
	if len(angle2) != nima:
		ERROR( "List lengths do not agree!", "angle_diff_sym",1)

	cosi = 0.0
	sini = 0.0
	agree = 0
	for i in range(nima):
		if( ( (angle2[i][1] <90.0) and (angle1[i][1] <90.0) ) or ( (angle2[i][1] >90.0) and (angle1[i][1] >90.0) ) ):
			qt = radians((angle2[i][0]-angle1[i][0])*simi)
			cosi += cos( qt )
			sini += sin( qt )
			agree += 1
	if(agree == 0):  return 0.0
	else:            return degrees(atan2(sini, cosi)/simi)%(360.0/simi)

def angle_error(ang1, ang2, delta_ang=0.0):
	'''
	This function calculates the error (variance) between two sets of angles after delta_ang (angle difference) is added to the
	first sets. When the angle difference (delta_ang) is the true difference, this function will return maximum error.
	'''
	from math import cos, sin, pi, radians
	
	erra = 0.0
	errb = 0.0
	delta_ang = radians( delta_ang )
	for i in range(len(ang1)):
		p2   = radians( ang2[i] )
		p2_x = cos(p2)
		p2_y = sin(p2)
		p1   = radians( ang1[i] )
		p1_x = cos(p1)
		p1_y = sin(p1)

		erra += p2_x*p1_x+p2_y*p1_y
		errb += p2_y*p1_x-p2_x*p1_y

	return erra*cos(delta_ang) + errb*sin(delta_ang)


def align_diff_params(ali_params1, ali_params2):
	'''
	This function determines the relative angle, shifts and mirrorness between
	two sets of alignment parameters.	
	'''
	from math import cos, sin, pi
	from sp_utilities import combine_params2
	
	nima = len(ali_params1)
	nima2 = len(ali_params2)
	if nima2 != nima:
		sxprint("Error: Number of images do not agree!")
		return 0.0, 0.0, 0.0, 0
	else:
		nima/=4
		del nima2

	# Read the alignment parameters and determine the relative mirrorness
	mirror_same = 0
	for i in range(nima):
		if ali_params1[i*4+3] == ali_params2[i*4+3]: mirror_same += 1
	if mirror_same > nima/2:
		mirror = 0
	else:
		mirror_same = nima-mirror_same
		mirror = 1

	# Determine the relative angle
	cosi = 0.0
	sini = 0.0
	angle1 = []
	angle2 = []
	for i in range(nima):
		mirror1 = ali_params1[i*4+3]
		mirror2 = ali_params2[i*4+3]
		if abs(mirror1-mirror2) == mirror:
			alpha1 = ali_params1[i*4]
			alpha2 = ali_params2[i*4]
			if mirror1 == 1:
				alpha1 = -alpha1
				alpha2 = -alpha2
			angle1.append(alpha1)
			angle2.append(alpha2)
	alphai = angle_diff(angle1, angle2)

	# Determine the relative shift
	sxi = 0.0
	syi = 0.0
	for i in range(nima):
		mirror1 = ali_params1[i*4+3]
		mirror2 = ali_params2[i*4+3]
		if abs(mirror1-mirror2) == mirror:
			alpha1 = ali_params1[i*4]
			#alpha2 = ali_params2[i*4]
			sx1 = ali_params1[i*4+1]
			sx2 = ali_params2[i*4+1]
			sy1 = ali_params1[i*4+2]
			sy2 = ali_params2[i*4+2]
			alpha12, sx12, sy12, mirror12 = combine_params2(alpha1, sx1, sy1, int(mirror1), alphai, 0.0, 0.0, 0)
			if mirror1 == 0: sxi += sx2-sx12
			else: sxi -= sx2-sx12
			syi += sy2-sy12

	sxi /= mirror_same
	syi /= mirror_same

	return alphai, sxi, syi, mirror


def align_diff(data1, data2=None, suffix="_ideal"):
	
	'''
	This function determines the relative angle, shifts and mirrorness between
	two list of data
	'''
	from sp_utilities import get_params2D
	
	nima = len(data1)

	if data2 != None: 
		nima2 = len(data2)
		if nima2 != nima:
			sxprint("Error: Number of images don't agree!")
			return 0.0, 0.0, 0.0, 0
		else:
			del nima2

	# Read the alignment parameters and determine the relative mirrorness
	ali_params1 = []
	ali_params2 = []
	for i in range(nima):
		alpha1, sx1, sy1, mirror1, scale1 = get_params2D(data1[i])
		if data2 != None:
			alpha2, sx2, sy2, mirror2, scale2 = get_params2D(data2[i])
		else:
			alpha2, sx2, sy2, mirror2, scale2 = get_params2D(data1[i], "xform.align2d"+suffix)
		ali_params1.extend([alpha1, sx1, sy1, mirror1])
		ali_params2.extend([alpha2, sx2, sy2, mirror2])

	return align_diff_params(ali_params1, ali_params2)


def align_diff_textfile(textfile1, textfile2):
	
	'''
	This function (2D) determines the relative angle, shifts and mirrorness between
	the two textfile of alignment parameters
	'''
	from sp_utilities import read_text_row
	
	ali1 = read_text_row(textfile1, "", "")
	ali2 = read_text_row(textfile2, "", "")

	nima = len(ali1)
	nima2 = len(ali2)
	if nima2 != nima:
		sxprint("Error: Number of images don't agree!")
		return 0.0, 0.0, 0.0, 0
	else:
		del nima2

	# Read the alignment parameters and determine the relative mirrorness
	ali_params1 = []
	ali_params2 = []
	for i in range(nima):
		ali_params1.extend(ali1[i][0:4])
		ali_params2.extend(ali2[i][0:4])

	return align_diff_params(ali_params1, ali_params2)


def ave_ali_err(data1, data2=None, r=25, suffix="_ideal"):
	'''
	This function determines the relative angle, shifts and mirrorness between
	the two lists of data. It also calculates the mirror consistent
	rate and average pixel error between two sets of parameters.
	'''
	from sp_utilities import get_params2D, combine_params2
	from math import sqrt, sin, pi
	
	# Determine relative angle, shifts and mirror
	alphai, sxi, syi, mirror = align_diff(data1, data2, suffix)

	# Determine the average pixel error
	err = 0.0
	nima = len(data1)
	mirror_same = 0
	for i in range(nima):
		alpha1, sx1, sy1, mirror1, scale1 = get_params2D(data1[i])
		if data2 != None:
			alpha2, sx2, sy2, mirror2, scale2 = get_params2D(data2[i])
		else:
			alpha2, sx2, sy2, mirror2, scale2 = get_params2D(data1[i], "xform.align2d"+suffix)
		
		if abs(mirror1-mirror2) == mirror: 
			mirror_same += 1
			alpha12, sx12, sy12, mirror12 = combine_params2(alpha1, sx1, sy1, int(mirror1), alphai, sxi, syi, 0)
			err += pixel_error_2D([alpha12, sx12, sy12], [alpha2, sx2, sy2], r)
	
	return alphai, sxi, syi, mirror, float(mirror_same)/nima, err/mirror_same


def ave_ali_err_params(ali_params1, ali_params2, r=25):
	'''
	This function determines the relative angle, shifts and mirrorness between
	the two sets of alignment parameters. It also calculates the mirror consistent
	rate and average pixel error between two sets of parameters.
	'''
	from sp_utilities import combine_params2
	from math import sqrt, sin, pi

	# Determine relative angle, shift and mirror
	alphai, sxi, syi, mirror = align_diff_params(ali_params1, ali_params2)

	# Determine the average pixel error
	nima = len(ali_params1)/4
	mirror_same = 0
	err = 0.0
	for i in range(nima):
		alpha1, sx1, sy1, mirror1 = ali_params1[i*4:i*4+4]
		alpha2, sx2, sy2, mirror2 = ali_params2[i*4:i*4+4]

		if abs(mirror1-mirror2) == mirror: 
			mirror_same += 1
			alpha12, sx12, sy12, mirror12 = combine_params2(alpha1, sx1, sy1, int(mirror1), alphai, sxi, syi, 0)
			err += pixel_error_2D([alpha12, sx12, sy12], [alpha2, sx2, sy2], r)

	return alphai, sxi, syi, mirror, float(mirror_same)/nima, err/mirror_same


def ave_ali_err_textfile(textfile1, textfile2, r=25):
	'''
	This function determines the relative angle, shifts and mirrorness between
	the two sets of alignment parameters. It also calculates the mirror consistent
	rate and average pixel error between two sets of parameters.
	'''
	from sp_utilities import combine_params2
	from math import sqrt, sin, pi
	from sp_utilities import read_text_row
	
	ali1 = read_text_row(textfile1, "", "")
	ali2 = read_text_row(textfile2, "", "")

	nima = len(ali1)
	nima2 = len(ali2)
	if nima2 != nima:
		sxprint("Error: Number of images don't agree!")
		return 0.0, 0.0, 0.0, 0, 0.0, 0.0
	else:
		del nima2

	# Read the alignment parameters
	ali_params1 = []
	ali_params2 = []
	for i in range(nima):
		ali_params1.extend(ali1[i][0:4])
		ali_params2.extend(ali2[i][0:4])

	# Determine relative angle, shift and mirror
	alphai, sxi, syi, mirror = align_diff_params(ali_params1, ali_params2)

	# Determine the average pixel error
	nima = len(ali_params1)/4
	mirror_same = 0
	err = 0.0
	for i in range(nima):
		alpha1, sx1, sy1, mirror1 = ali_params1[i*4:i*4+4]
		alpha2, sx2, sy2, mirror2 = ali_params2[i*4:i*4+4]
		
		if abs(mirror1-mirror2) == mirror: 
			mirror_same += 1
			alpha12, sx12, sy12, mirror12 = combine_params2(alpha1, sx1, sy1, int(mirror1), alphai, sxi, syi, 0)
			err += pixel_error_2D([alpha12, sx12, sy12], [alpha2, sx2, sy2], r)
	
	return alphai, sxi, syi, mirror, float(mirror_same)/nima, err/mirror_same


def multi_align_diff_params(ali_params, verbose=0):
	"""
	Calculate the mirror consistency and pixel error between different runs of alignment
	The input is in the following format (here n is the number of alignments done):
		[[alpha1_r1, sx1_s1, sy1_r1, mirror1_r1, alpha2_r1, sx2_r1, sy2_r1, mirror2_r1, ...],
		 [alpha1_r2, sx1_s2, sy1_r2, mirror1_r2, alpha2_r2, sx2_r2, sy2_r2, mirror2_r2, ...],
		 [alpha1_r3, sx1_s3, sy1_r3, mirror1_r3, alpha2_r3, sx2_r3, sy2_r3, mirror2_r3, ...],
		 ...
		 [alpha1_rn, sx1_sn, sy1_rn, mirror1_rn, alpha2_rn, sx2_rn, sy2_rn, mirror2_rn, ...]]
	The output is in the following format (here k=n*(n+1)/2):
		[[pixel_error_1, mirror_consistency_1, i_1, j_1, alpha_1, sx_1, sy_1, mirror_1],
		 [pixel_error_2, mirror_consistency_2, i_2, j_2, alpha_2, sx_2, sy_2, mirror_2],
		 [pixel_error_3, mirror_consistency_3, i_3, j_3, alpha_3, sx_3, sy_3, mirror_3],
		...
		 [pixel_error_k, mirror_consistency_k, i_k, j_k, alpha_k, sx_k, sy_k, mirror_k]]
	"""
	num_ali = len(ali_params)
	multi_align_results = []
	for i in range(num_ali-1):
		for j in range(i+1, num_ali):
			alpha, sx, sy, mirror, stab_mirror, pixel_error = ave_ali_err_params(ali_params[i], ali_params[j])
			if verbose == 1:
				sxprint("Between trial %d and %d: mirror stability = %6.3f   pixel error = %6.3f"%(i, j, stab_mirror, pixel_error))
			multi_align_results.append([pixel_error, stab_mirror, i, j, alpha, sx, sy, mirror])
	return multi_align_results
	

def calc_connect_list(multi_align_results, pixel_error_threshold = 5.0, mirror_consistency_threshold = 0.8):
	"""
	Generate the connection list from the multi_align_results, which generally comes from multi_align_diff_params()
	The connection list will have the following format:
		[[1, 2, 5], [4, 6], [0, 7]]
	You will also get the size of the largest connection in the list.
	"""
	import sets
	
	k = len(multi_align_results)
	multi_align_results.sort()
	connect_list = []
	for i in range(k):
		if multi_align_results[i][0] <= pixel_error_threshold:
			if multi_align_results[i][1] >= mirror_consistency_threshold: 
				connect_list.append([multi_align_results[i][2], multi_align_results[i][3]])
		else:	break
	to_break = True
	while to_break:
		l = len(connect_list)
		to_break = False
		for i in range(l-1):
			for j in range(i+1, l):
				set1 = set(connect_list[i])
				set2 = set(connect_list[j])
				if list(set1.intersection(set2)) != []:
					connect_list[i] = list(set1.union(set2))
					del connect_list[j]
					to_break = True
					break
			if to_break: break
	max_connect = 0
	for l in connect_list: max_connect = max(max_connect, len(l))
	return connect_list, max_connect


def ali_stable_list(ali_params1, ali_params2, pixel_error_threshold, r=25):
	'''
	This function first determines the relative angle, shifts and mirrorness between
	the two sets of alignment parameters. It then determines whether each image is
	stable or not and return this information as an int list. (1 is stable and 0 is unstable)
	'''
	from sp_utilities import combine_params2
	from math import sqrt, sin, pi
	
	# Determine relative angle, shift and mirror
	alphai, sxi, syi, mirror = align_diff_params(ali_params1, ali_params2)

	# Determine the average pixel error
	nima = len(ali_params1)/4
	ali_list = []
	for i in range(nima):
		alpha1, sx1, sy1, mirror1 = ali_params1[i*4:i*4+4]
		alpha2, sx2, sy2, mirror2 = ali_params2[i*4:i*4+4]
		if abs(mirror1-mirror2) == mirror:
			alpha12, sx12, sy12, mirror12 = combine_params2(alpha1, sx1, sy1, int(mirror1), alphai, sxi, syi, 0)
			if pixel_error_2D([alpha12, sx12, sy12], [alpha2, sx2, sy2], r) < pixel_error_threshold: ali_list.append(1)
			else: ali_list.append(0)
		else: ali_list.append(0)

	return ali_list



def multi_align_stability(ali_params, mir_stab_thld = 0.0, grp_err_thld = 10000.0, err_thld = 1.732, print_individual = False, d = 64):

	def sqerr(a):
		n = len(a)
		avg = sum(a)
		sq = 0.0
		for i in range(n): sq += a[i]**2
		return (sq-avg*avg/n)/n

	# args - G, data - [T, d]
	def func(args, data, return_avg_pixel_error=True):

		from math import pi, sin, cos, radians, degrees

		ali_params = data[0]
		d = data[1]
		L = len(ali_params)
		N = len(ali_params[0])/4

		args_list= [0.0]*(L*3)
		for i in range(L*3-3):  args_list[i] = args[i]

		cosa = [0.0]*L
		sina = [0.0]*L
		for i in range(L):
			cosa[i] = cos(radians(args_list[i*3]))
			sina[i] = sin(radians(args_list[i*3]))
		sqr_pixel_error = [0.0]*N
		ave_params = []
		for i in range(N):
			sum_cosa = 0.0
			sum_sina = 0.0
			sx = [0.0]*L
			sy = [0.0]*L
			for j in range(L):
				if int(ali_params[j][i*4+3]) == 0:
					sum_cosa += cos(radians(args_list[j*3]+ali_params[j][i*4]))
					sum_sina += sin(radians(args_list[j*3]+ali_params[j][i*4]))
					sx[j] =  args_list[j*3+1] + ali_params[j][i*4+1]*cosa[j] + ali_params[j][i*4+2]*sina[j]
					sy[j] =  args_list[j*3+2] - ali_params[j][i*4+1]*sina[j] + ali_params[j][i*4+2]*cosa[j]
				else:
					sum_cosa += cos(radians(-args_list[j*3]+ali_params[j][i*4]))
					sum_sina += sin(radians(-args_list[j*3]+ali_params[j][i*4]))
					sx[j] = -args_list[j*3+1] + ali_params[j][i*4+1]*cosa[j] - ali_params[j][i*4+2]*sina[j]
					sy[j] =  args_list[j*3+2] + ali_params[j][i*4+1]*sina[j] + ali_params[j][i*4+2]*cosa[j]
			sqrtP = sqrt(sum_cosa**2+sum_sina**2)
			sqr_pixel_error[i] = max( 0.0, d*d/4.*(1-sqrtP/L) + sqerr(sx) + sqerr(sy) )
			# Get ave transform params
			H = Transform({"type":"2D"})
			H.set_matrix([sum_cosa/sqrtP, sum_sina/sqrtP, 0.0, sum(sx)/L, -sum_sina/sqrtP, sum_cosa/sqrtP, 0.0, sum(sy)/L, 0.0, 0.0, 1.0, 0.0])
			dd = H.get_params("2D")
			#  We are using here mirror of the LAST SET.
			H = Transform({"type":"2D","alpha":dd[ "alpha" ],"tx":dd[ "tx" ],"ty": dd[ "ty" ],"mirror":int(ali_params[L-1][i*4+3]),"scale":1.0})
			dd = H.get_params("2D")
			ave_params.append([dd[ "alpha" ], dd[ "tx" ], dd[ "ty" ], dd[ "mirror" ]])
		# Warning: Whatever I return here is squared pixel error, this is for the easy expression of derivative
		# Don't forget to square root it after getting the value
		if return_avg_pixel_error:
			return sum(sqr_pixel_error)/N
		else:
			return sqr_pixel_error, ave_params

	"""Multiline Comment0"""
	#MULTILINEMULTILINEMULTILINE 0

	        #MULTILINEMULTILINEMULTILINE 0
	        #MULTILINEMULTILINEMULTILINE 0

	        #MULTILINEMULTILINEMULTILINE 0

	        #MULTILINEMULTILINEMULTILINE 0
	        #MULTILINEMULTILINEMULTILINE 0

	        #MULTILINEMULTILINEMULTILINE 0
	        #MULTILINEMULTILINEMULTILINE 0

	        #MULTILINEMULTILINEMULTILINE 0
	        #MULTILINEMULTILINEMULTILINE 0
	        #MULTILINEMULTILINEMULTILINE 0
	        #MULTILINEMULTILINEMULTILINE 0
	        #MULTILINEMULTILINEMULTILINE 0
	        	#MULTILINEMULTILINEMULTILINE 0
	        	#MULTILINEMULTILINEMULTILINE 0
	        #MULTILINEMULTILINEMULTILINE 0
	        	#MULTILINEMULTILINEMULTILINE 0
	        	#MULTILINEMULTILINEMULTILINE 0
	        	#MULTILINEMULTILINEMULTILINE 0
	        	#MULTILINEMULTILINEMULTILINE 0
	        	#MULTILINEMULTILINEMULTILINE 0
	        		#MULTILINEMULTILINEMULTILINE 0
	        			#MULTILINEMULTILINEMULTILINE 0
	        			#MULTILINEMULTILINEMULTILINE 0
	        			#MULTILINEMULTILINEMULTILINE 0
	        			#MULTILINEMULTILINEMULTILINE 0
	        		#MULTILINEMULTILINEMULTILINE 0
	        			#MULTILINEMULTILINEMULTILINE 0
	        			#MULTILINEMULTILINEMULTILINE 0
	        			#MULTILINEMULTILINEMULTILINE 0
	        			#MULTILINEMULTILINEMULTILINE 0
	        	#MULTILINEMULTILINEMULTILINE 0
	        	#MULTILINEMULTILINEMULTILINE 0
	        	#MULTILINEMULTILINEMULTILINE 0

	        	#MULTILINEMULTILINEMULTILINE 0
	        		#MULTILINEMULTILINEMULTILINE 0
	        		#MULTILINEMULTILINEMULTILINE 0
	        		#MULTILINEMULTILINEMULTILINE 0
	        		#MULTILINEMULTILINEMULTILINE 0
	        		#MULTILINEMULTILINEMULTILINE 0
	        		#MULTILINEMULTILINEMULTILINE 0
	        		#MULTILINEMULTILINEMULTILINE 0
	        		#MULTILINEMULTILINEMULTILINE 0
	        			#MULTILINEMULTILINEMULTILINE 0
	        					#MULTILINEMULTILINEMULTILINE 0
	        					#MULTILINEMULTILINEMULTILINE 0
	        					#MULTILINEMULTILINEMULTILINE 0
	        			#MULTILINEMULTILINEMULTILINE 0
	        			#MULTILINEMULTILINEMULTILINE 0
	        		#MULTILINEMULTILINEMULTILINE 0
	        			#MULTILINEMULTILINEMULTILINE 0
	        					    #MULTILINEMULTILINEMULTILINE 0
	        					 #MULTILINEMULTILINEMULTILINE 0
	        					 #MULTILINEMULTILINEMULTILINE 0
	        			#MULTILINEMULTILINEMULTILINE 0
	        			#MULTILINEMULTILINEMULTILINE 0
	        #MULTILINEMULTILINEMULTILINE 0

	        #MULTILINEMULTILINEMULTILINE 0
	#MULTILINEMULTILINEMULTILINE 0

	from sp_statistics import k_means_stab_bbenum
	from sp_utilities import combine_params2
	from numpy import array
	from math import sqrt

	# I decided not to use scipy in order to reduce the dependency, I wrote the C++ code instead
	# from scipy import array, int32
	# from scipy.optimize.lbfgsb import fmin_l_bfgs_b

	# Find out the subset which is mirror stable over all runs
	all_part = []
	num_ali = len(ali_params)
	nima = len(ali_params[0])/4
	for i in range(num_ali):
		mirror0 = []
		mirror1 = []
		for j in range(nima):
			if ali_params[i][j*4+3] == 0: mirror0.append(j)
			else: mirror1.append(j)
		mirror0 = array(mirror0, 'int32')
		mirror1 = array(mirror1, 'int32')
		all_part.append([mirror0, mirror1])
	match, stab_part, CT_s, CT_t, ST, st = k_means_stab_bbenum(all_part, T=0, nguesses=1)
	mir_stab_part = stab_part[0] + stab_part[1]

	mir_stab_rate = len(mir_stab_part)/float(nima)
	if mir_stab_rate <= mir_stab_thld: return [], mir_stab_rate, -1.0
	mir_stab_part.sort()


	del all_part, match, stab_part, CT_s, CT_t, ST, st	

	#for i in xrange(len(mir_stab_part)):  print i, mir_stab_part[i]

	nima2 = len(mir_stab_part)

	#print "  mirror stable  ",nima2


	# Keep the alignment parameters of mirror stable particles
	ali_params_mir_stab = [[] for i in range(num_ali)]
	for i in range(num_ali):
		for j in mir_stab_part:
			ali_params_mir_stab[i].extend(ali_params[i][j*4:j*4+4])
	# Find out the alignment parameters for each iteration against the last one
	args = []
	for i in range(num_ali-1):
		alpha, sx, sy, mirror = align_diff_params(ali_params_mir_stab[i], ali_params_mir_stab[num_ali-1])
		args.extend([alpha, sx, sy])

	# Do an initial analysis, purge all outlier particles, whose pixel error are larger than three times the threshold
	data = [ali_params_mir_stab, d]
	pixel_error_before, ave_params = func(array(args), data, return_avg_pixel_error=False)
	#  We have to return mir_stab_rate (see above), even if the group does not survive it and the pixel error before cleaning,
	#   see below, and in ISAC print the user the overall statistics (histograms?) 
	#   so one can get on overall idea how good/bad data is.  PAP  01/25/2015
	#print  " >>> ",sqrt(sum(pixel_error_before)/nima2)
	ali_params_cleaned = [[] for i in range(num_ali)]
	cleaned_part = []
	for j in range(nima2):
		pixel_error_before[j] = max(0.0, pixel_error_before[j])  # prevent sqrt of 0
		if sqrt(pixel_error_before[j]) > 3*err_thld:
			pass  #print "  removed ",3*err_thld,j,sqrt(pixel_error_before[j])
		else:
			cleaned_part.append(mir_stab_part[j])
			for i in range(num_ali):
				ali_params_cleaned[i].extend(ali_params_mir_stab[i][j*4:j*4+4])
	nima3 = len(cleaned_part)
	prever = sqrt(sum(pixel_error_before)/nima2)
	if nima3 <= 1:  return [], mir_stab_rate, prever

	#print "  cleaned part  ",nima3


	# Use LBFGSB to minimize the sum of pixel errors
	data = [ali_params_cleaned, d]
	# Use Python code
	#ps_lp, val, d = fmin_l_bfgs_b(func, array(args), args=[data], fprime=dfunc, bounds=None, m=10, factr=1e3, pgtol=1e-4, iprint=-1, maxfun=100)
	# Use C++ code
	ali_params_cleaned_list = []
	for params in ali_params_cleaned: ali_params_cleaned_list.extend(params)
	results = Util.multi_align_error(args, ali_params_cleaned_list, d)
	ps_lp = results[:-1]

	# Negative val can happen in some rare cases, it should be due to rounding errors, 
	# because all results show the val is about 1e-13.
	#print "Strange results"
	#print "args =", args
	#print "ali_params_cleaned_list =", ali_params_cleaned_list
	#print "results = ", results
	val = max(0.0, results[-1])

	del ali_params_cleaned_list

	if sqrt(val) > grp_err_thld: return [], mir_stab_rate, sqrt(val)

	pixel_error_after, ave_params = func(ps_lp, data, return_avg_pixel_error=False)

	stable_set = []
	val = 0.0
	for i in range(nima):
		if i in cleaned_part:
			j = cleaned_part.index(i)
			err = sqrt(pixel_error_after[j])
			if err < err_thld:
				stable_set.append([err, i, ave_params[j]])
				val += err
				if print_individual:  sxprint("Particle %4d :  pixel error = %18.4f"%(i, err))
			else:
				if print_individual:  sxprint("Particle %4d :  pixel error = %18.4f  unstable"%(i, err))
		else:
			if print_individual:  sxprint("Particle %4d :  Mirror unstable"%i)
	#  return average pixel error before pruning as it is more informative
	return stable_set, mir_stab_rate, prever# sqrt(val/len(stable_set))


# args - G, data - [T, d]
def ave2dtransform(args, data, return_avg_pixel_error=False):

	from math import pi, sin, cos, radians, degrees

	ali_params = data[0]
	d = data[1]
	L = len(ali_params)
	N = len(ali_params[0])/4

	args_list= [0.0]*(L*3)
	for i in range(L*3-3):  args_list[i] = args[i]

	cosa = [0.0]*L
	sina = [0.0]*L
	for i in range(L):
		cosa[i] = cos(radians(args_list[i*3]))
		sina[i] = sin(radians(args_list[i*3]))
	sqr_pixel_error = [0.0]*N
	ave_params = []
	for i in range(N):
		sum_cosa = 0.0
		sum_sina = 0.0
		sx = [0.0]*L
		sy = [0.0]*L
		for j in range(L):
			if int(ali_params[j][i*4+3]) == 0:
				sum_cosa += cos(radians(args_list[j*3]+ali_params[j][i*4]))
				sum_sina += sin(radians(args_list[j*3]+ali_params[j][i*4]))
				sx[j] =  args_list[j*3+1] + ali_params[j][i*4+1]*cosa[j] + ali_params[j][i*4+2]*sina[j]
				sy[j] =  args_list[j*3+2] - ali_params[j][i*4+1]*sina[j] + ali_params[j][i*4+2]*cosa[j]
			else:
				sum_cosa += cos(radians(-args_list[j*3]+ali_params[j][i*4]))
				sum_sina += sin(radians(-args_list[j*3]+ali_params[j][i*4]))
				sx[j] = -args_list[j*3+1] + ali_params[j][i*4+1]*cosa[j] - ali_params[j][i*4+2]*sina[j]
				sy[j] =  args_list[j*3+2] + ali_params[j][i*4+1]*sina[j] + ali_params[j][i*4+2]*cosa[j]
		sqrtP = sqrt(sum_cosa**2+sum_sina**2)
		sqr_pixel_error[i] = max( 0.0, d*d/4.*(1-sqrtP/L) + sqerr(sx) + sqerr(sy) )
		# Get ave transform params
		H = Transform({"type":"2D"})
		H.set_matrix([sum_cosa/sqrtP, sum_sina/sqrtP, 0.0, sum(sx)/L, -sum_sina/sqrtP, sum_cosa/sqrtP, 0.0, sum(sy)/L, 0.0, 0.0, 1.0, 0.0])
		dd = H.get_params("2D")
		#  We are using here mirror of the LAST SET.
		H = Transform({"type":"2D","alpha":dd[ "alpha" ],"tx":dd[ "tx" ],"ty": dd[ "ty" ],"mirror":int(ali_params[L-1][i*4+3]),"scale":1.0})
		dd = H.get_params("2D")
		ave_params.append([dd[ "alpha" ], dd[ "tx" ], dd[ "ty" ], dd[ "mirror" ]])
	# Warning: Whatever I return here is squared pixel error, this is for the easy expression of derivative
	# Don't forget to square root it after getting the value
	if return_avg_pixel_error:
		return sum(sqr_pixel_error)/N
	else:
		return sqr_pixel_error, ave_params

def rotate_angleset_to_match(agls1, agls2):
	"""
	  Finds rotation between two sets of angles, agls2 is the template
	  It will also establish whether mirror is required
	  Rotation is applied to agsl1 and the set of rotated angles is returned
	  Rotation itself is not returned.
	  Makes sense only for no symmetry
	"""
	from sp_utilities    import rotation_between_anglesets
	from sp_fundamentals import rotate_params

	t1 = rotation_between_anglesets(agls1, agls2)

	return rotate_params(agls1,[-t1[2],-t1[1],-t1[0]])

def ordersegments(infilaments, ptclcoords):
	'''
	Input:
	
	stack: Input stack of images whose headers contain filament membership information and particle coordinates in the original micrograph (stored under attribute ptcl_source_coord).
	filament_attr: Attribute under which filament membership ID is stored.
	It is assumed the prtl coords are nonnegative
	
	Output: 
	
	Returns a list of lists, where each inner list consists of IDs of segments windowed from
	a single filament ordered according to their relative positions on the filament.
	
	'''

	def orderbymodule(xxp,yyp):
		from math import atan,sin,cos,pi, atan2
		from sp_statistics import linreg
		nq = len(xxp)
		xs = sum(xxp)/nq
		ys = sum(yyp)/nq
		xp = [0.0]*nq
		yp = [0.0]*nq
		for i in range(nq):
			xp[i] = xxp[i] - xs
			yp[i] = yyp[i] - ys
		try:
			a,b = linreg(xp,yp)
			alpha = pi/4-atan(a)
		except:
			a,b = linreg([(xp[i]-yp[i]) for i in range(nq)], [(xp[i]+yp[i]) for i in range(nq)])
			alpha = atan(a)
			#print "except"

		cs = cos(alpha)
		ss = sin(alpha)
		qm = 1.e23
		dd = [[0.0, 0] for i in range(nq)]
		for i in range(nq):
			xt =  cs*xp[i] - ss*yp[i]
			yt =  ss*xp[i] + cs*yp[i]
			xp[i] = xt; yp[i] = yt
		xs = min(xp)
		ys = min(yp)
		for i in range(nq):
			dd[i] = [(xp[i]-xs)**2+(yp[i]-ys)**2, i]
		dd.sort()
		return [dd[i][1] for i in range(nq)]

	allfilaments = [None]*len(infilaments)
	for i in range(len(infilaments)):
		allfilaments[i] = [infilaments[i],i]
	allfilaments.sort()
	filaments = []
	current = allfilaments[0][0]
	temp = [allfilaments[0][1]]
	for i in range(1,len(allfilaments)):
		if( allfilaments[i][0] == current ):
			temp.extend([allfilaments[i][1]])
		else:
			filaments.append(temp)
			current = allfilaments[i][0]
			temp = [allfilaments[i][1]]
	filaments.append(temp)

	del allfilaments,temp

	nfil = len(filaments)

	for i in range(nfil):
		nsegs = len(filaments[i])
		if(nsegs > 1):
			ord = orderbymodule([ptclcoords[filaments[i][ii]][0] for ii in range(nsegs)],[ptclcoords[filaments[i][ii]][1] for ii in range(nsegs)])
			filaments[i] = [filaments[i][ord[ii]] for ii in range(nsegs)]
			# To preserve the original order check indexes and invert if it appears to be inverted
			if(filaments[i][0] > filaments[i][-1]):
				for k in range(nsegs//2):
					temp = filaments[i][k]
					filaments[i][k] = filaments[i][nsegs-1-k]
					filaments[i][nsegs-1-k]=temp
	return filaments	


def mapcoords(x, y, r, nx, ny):
	from math 			import ceil, floor
	from sp_utilities 	import get_dist
	import sys
	'''
	Input:
	
	(x,y): Coordinate in old image. 
	r: ratio by which old image is resampled by. If r < 1, then pixel size of resampled image is original pixel size divided by r.
	nx, ny: dimensions of old image
	
	Assumes coordinates are positive and run from 0 to nx-1 and 0 to ny-1, where nx and ny
	are the x and y dimensions of the micrograph respectively.
	
	Output:
	
	x'',y'':	The pixel coordinate in the resampled image where
	
					(x',y') = (Util.round(x''/r), Util.round(y''/r))
					
				and (x',y') is the closest point to (x,y) over all other points (a,b) in the
				old image where (a,b)=  (Util.round(a''/r), Util.round(b''/r)) for some
				pixel coordinate (a'', b'') in resampled image.
	'''	
	
	# Neighborhood of (x,y) in old image which contains a point (x',y') such that
	# (x',y') = (Util.round(x''/r), Util.round(y''/r)) for some (x'', y'') in resampled image
	if r > 1:
		nbrhd = 1
	else:
		nbrhd = int(ceil(1.0/r))+1
			
	allxnew = []
	allynew = []
	
	for i in range(-nbrhd, nbrhd+1):
		xold = Util.round(x + i)
		if xold < 0 or xold >= nx:
			continue
		# See if there is xnew in new image where xold == int(Util.round(xnew/r))
		# If there is such a xnew, then xold == int(Util.round(xnew/r)) implies r*(xold-0.5) <= xnew < r*(xold+0.5)
		lxnew = int(floor(r*(xold - 0.5)))
		uxnew = int(ceil(r*(xold + 0.5))) 
		for xn in range(lxnew, uxnew + 1):
			if xold == Util.round(xn/r):
				allxnew.append(xn)
				
	for j in range(-nbrhd, nbrhd+1):
		yold = Util.round(y + j)
		if yold < 0 or yold >= ny:
			continue
		lynew = int(floor(r*(yold - 0.5)))
		uynew = int(ceil(r*(yold + 0.5)))
		for yn in range(lynew, uynew + 1):
			if yold == Util.round(yn/r):
				allynew.append(yn)
				
	if len(allxnew) == 0 or len(allynew) == 0:
		ERROR("Could not find mapping")
	
	mindist = -1
	minxnew = -1
	minynew = -1
	
	for xnew in allxnew:
		for ynew in allynew:
			xold = Util.round(xnew/r)
			yold = Util.round(ynew/r)
			dst = get_dist([x,y],[xold,yold])
			if dst > mindist:
				mindist = dst
				minxnew = int(xnew)
				minynew = int(ynew)
					
	return minxnew, minynew

def consistency_params(stack, dphi, dp, pixel_size, phithr=2.5, ythr=1.5, THR=3):
	'''
		stack        - contains coding of filaments and coordinates of segments ptcl_source_coord
		fname_params - parameters whose consistency is tested
	'''
	from sp_utilities import read_text_row, write_text_row, get_dist
	from sp_applications import ordersegments
	from sp_pixel_error import angle_diff

	filaments = ordersegments(stack)
	ptclcoords = EMUtil.get_all_attributes(stack, 'ptcl_source_coord')
	params     = EMUtil.get_all_attributes(stack, 'xform.projection')
	for i in range(len(params)):
		d = params[i].get_params("spider")
		params[i] = [d["phi"], d["theta"], d["psi"], -d["tx"], -d["ty"] ]

	N = len(filaments)
	sxprint("N: ", N)
	totsegs    = 0
	totpsicons = 0
	totphicons = 0
	tot_nopred = 0
	allphier = []
	imi = 0
	for mic in filaments:
		#for kkk in xrange(1):
		#mic = filaments[0]
		imi += 1
		#mic = mic[:4]
		if len(mic) < THR:
			#print "less than threshold"
			allphier.append([[mic[0],mic[-1],0],[]])
		else:

			totsegs += len(mic)

			a90  = [] # segments in this filament with psi ~ 90
			a270 = [] # segments in this filament with psi ~ 270

			for iseg in mic:
				if abs(params[iseg][2] - 90.0) <45.0: a90.append(iseg)
				else:                                 a270.append(iseg)

			if (len(a90) == len(a270)):
				# these cannot be predicted
				tot_nopred += len(mic)
				continue

			if( len(a90) > len(a270) ):
				thetapsi1 =  a90
				flip =  1
			else:
				thetapsi1 = a270
				flip = -1

			ns = len(thetapsi1)
			# given phi = phig
			phig = [0.0]*ns
			for j in range(ns): phig[j] = params[thetapsi1[j]][0]
			totpsicons += ns
			distances = [0.0]*ns
			for i in range(1,ns):  distances[i] = get_dist( ptclcoords[thetapsi1[i]], ptclcoords[thetapsi1[0]] )
			ganger = [0.0]*ns
			terr = 1.e23
			for idir in range(-1,2,2):
				phierr = []
				#  get phi's
				ddphi = pixel_size/dp*idir*dphi
				phis = [0.0]*ns
				#print "  MIC  ",mic
				for i in range(ns):
					yy = distances[i] + params[thetapsi1[i]][4]
					phis[i] = (yy*ddphi)%360.0
					#print " %7.3f   %7.3f   %7.3f  %7.3f   "%(ycoords[i],yy,phis[i], phig[i]),params[thetapsi1[i]]

				# find the overall angle
				angdif = angle_diff(phis,phig)
				#print " angdif ",angdif
				lerr = 0.0
				for i in range(ns):
					anger = (phis[i]+angdif - phig[i] + 360.0)%360.0
					if( anger > 180.0 ): anger -= 360.0
					lerr += abs(anger)
					phierr.append(anger)
					#print  " %7.3f   %7.3f   %7.3f"%((phis[i]+angdif+360.0)%360.0 , phig[i],anger)
				if(lerr < terr):
					terr = lerr
					for j in range(ns):  ganger[j] = phierr[j]
			allphier.append([[mic[0], mic[-1], flip, terr/ns], ganger])

	sxprint("number of segments belonging to filaments from which at least %i segments were windowed: "%THR, totsegs)
	sxprint("number of segments oriented 50/50 wrt psi (and therefore could not be predicted):       ", tot_nopred)
	sxprint("segments whose psi agreed with the majority of segments in its filament:                ", totpsicons)
	return  allphier

def getnewhelixcoords(hcoordsname, outdir, ratio,nx,ny, newpref="resampled_", boxsize=-1):
	"""
	Input
	
		helixcoordsfile: Full path name of file with coordinates of boxed helices
		
		outdir: Full path name of directory in which to put new helix coordinates file.
		
		ratio: factor by which new image (micrograph) is resampled from old
		
		nx, ny: dimensions of old image (micrograph)
		
		newpref: prefix for attaching to fname to get name of new helix coordinates file
	
	Output:
		Returns full path name of file containing new box coordinates
	"""
	import os
	from sp_utilities 		import read_text_row
	from sp_pixel_error	import mapcoords
	
	fname = (hcoordsname.split('/'))[-1] # name of old coordinates files minus the path
	newhcoordsname = os.path.join(outdir , newpref+fname) # full path name of new coordinates file to be created
	f = open(newhcoordsname, 'w')
	coords = read_text_row(hcoordsname) # old coordinates
	ncoords = len(coords)
	newcoords=[]
	w = coords[0][2]
	new_w = boxsize
	if new_w < 0:
		new_w = w*ratio
	for i in range(ncoords):
		xold = coords[i][0] + w/2
		yold = coords[i][1] + w/2
		xnew, ynew = mapcoords(xold,yold,ratio,nx,ny)
		s = '%d\t%d\t%d\t%d\t%d\n'%(xnew-new_w/2,ynew-new_w/2, new_w, new_w, coords[i][4])
		f.write(s)
	return newhcoordsname	

def helical_params_err(params1, params2, fil_list):
	'''
	Input:
		params1: First set of projection orientation parameters
		params2: Second set of projection orientation parameters
		fil_list: List of filament IDs, where i-th element on list is filament of i-th image.
				  Assume fil_list is ordered, i.e., i-th element of fil_list is the filament ID
				  of the i-th image in the stack.
	
	Output:
		The phi angle difference is computed between params1 and params2 between those parameters
		whose psi agree, and the angle is applied to params2.
		
		For each filament, the program computes the phi error averaged over the number 
		of segments between params1 and the aligned params2.
		
		Returns a list of lists, where each inner list is [fil_ID, avg_phi_err], where fil_ID
		is filament name, and avg_phi_err is the average phi error for the filament.
	'''
	from sp_pixel_error import angle_diff
	from EMAN2 import Vec2f
	nima = len(params1)
	# Identify those where psi agrees
	phi1 = []
	phi2 = []
	fil_psisame = []
	pos = []
	pref = []
	for i in range(nima):
		if abs(params1[i][2] - params2[i][2]) < 90:
			phi1.append(params1[i][0])
			phi2.append(params2[i][0])
			fil_psisame.append(fil_list[i])
			pos.append(i)
			pref.append(params2[i])

	fag = len(fil_psisame)

	tflip = Transform({"type":"spider","theta":180.0})
	# Identify those where psi agrees
	phi1 = []
	phi2 = []
	fil_psisame = []
	pos  = []
	prot = []
	pref = []
	for i in range(nima):
		t2 = Transform({"type":"spider","phi":params1[i][0],"theta":params1[i][1],"psi":params1[i][2]})
		t2.set_trans( Vec2f( -params1[i][3], -params1[i][4] ) )
		t2 = t2*tflip
		d = t2.get_params("spider")
		p1r = [d["phi"], d["theta"], d["psi"], -d["tx"], -d["ty"]]
		if abs(p1r[2] - params2[i][2]) < 90:
			phi1.append(p1r[0])
			phi2.append(params2[i][0])
			fil_psisame.append(fil_list[i])
			pos.append(i)
			prot.append(p1r)
			pref.append(params2[i])
			
	nima2 = len(fil_psisame)
	if(fag > nima2):
		phi1 = []
		phi2 = []
		fil_psisame = []
		pos = []
		prot = []
		pref = []
		for i in range(nima):
			if abs(params1[i][2] - params2[i][2]) < 90:
				phi1.append(params1[i][0])
				phi2.append(params2[i][0])
				fil_psisame.append(fil_list[i])
				pos.append(i)
				prot.append(params1[i])
				pref.append(params2[i])
		nima2 = len(fil_psisame)
	else:
		sxprint("better agreement afer upside-down rotation")
		
	# agls1psi and agls2psi agree in psi. 
	sxprint("Number of images which agree on psi between params1 and params2: ",nima2)
	sxprint("Percentage of total number images being compared: ", nima2/float(len(params1))*100)
	
	angdif = angle_diff(phi1, phi2)
	sxprint("angdif: ", angdif)
	
	phierr_byfil = []

	start_fil = None
	i = 0
	while (i < nima2):
		ibeg = i
		iend = i
		cur_fil = fil_psisame[ibeg]
		for k in range(ibeg+1,nima2):
			if(cur_fil != fil_psisame[k]):
				iend = k
				break
		if( ibeg == iend ):  iend = nima2-1
		sum_phierr = 0.
		for k in range(ibeg,iend):
			phidf = abs(phi2[k] - (phi1[k] + angdif)%360.)
			phidf = min(phidf, abs(360.0-phidf))
			sum_phierr += phidf
			prot[k][0] =  (phi1[k] + angdif)%360.
		avg_phierr = sum_phierr/(iend-ibeg+1)
		phierr_byfil.append([avg_phierr, cur_fil,ibeg,iend,pos[ibeg],pos[iend]])
		i = iend+1

	return phierr_byfil,prot,pref

"""Multiline Comment1"""
#MULTILINEMULTILINEMULTILINE 1
	#MULTILINEMULTILINEMULTILINE 1
	#MULTILINEMULTILINEMULTILINE 1
	#MULTILINEMULTILINEMULTILINE 1
	#MULTILINEMULTILINEMULTILINE 1


	#MULTILINEMULTILINEMULTILINE 1
	#MULTILINEMULTILINEMULTILINE 1
	#MULTILINEMULTILINEMULTILINE 1

	#MULTILINEMULTILINEMULTILINE 1
	#MULTILINEMULTILINEMULTILINE 1
	#MULTILINEMULTILINEMULTILINE 1
	#MULTILINEMULTILINEMULTILINE 1
		#MULTILINEMULTILINEMULTILINE 1
		#MULTILINEMULTILINEMULTILINE 1

		#MULTILINEMULTILINEMULTILINE 1
		#MULTILINEMULTILINEMULTILINE 1
		#MULTILINEMULTILINEMULTILINE 1
		#MULTILINEMULTILINEMULTILINE 1
			#MULTILINEMULTILINEMULTILINE 1
				#MULTILINEMULTILINEMULTILINE 1
				#MULTILINEMULTILINEMULTILINE 1
				#MULTILINEMULTILINEMULTILINE 1
				#MULTILINEMULTILINEMULTILINE 1
				#MULTILINEMULTILINEMULTILINE 1
				#MULTILINEMULTILINEMULTILINE 1
				#MULTILINEMULTILINEMULTILINE 1
				#MULTILINEMULTILINEMULTILINE 1
				#MULTILINEMULTILINEMULTILINE 1
					#MULTILINEMULTILINEMULTILINE 1
					#MULTILINEMULTILINEMULTILINE 1
					#MULTILINEMULTILINEMULTILINE 1
					#MULTILINEMULTILINEMULTILINE 1
						#MULTILINEMULTILINEMULTILINE 1
						#MULTILINEMULTILINEMULTILINE 1
						#MULTILINEMULTILINEMULTILINE 1
						#MULTILINEMULTILINEMULTILINE 1
						#MULTILINEMULTILINEMULTILINE 1
						#MULTILINEMULTILINEMULTILINE 1
						#MULTILINEMULTILINEMULTILINE 1
						#MULTILINEMULTILINEMULTILINE 1
						#MULTILINEMULTILINEMULTILINE 1
						#MULTILINEMULTILINEMULTILINE 1
						#MULTILINEMULTILINEMULTILINE 1
						#MULTILINEMULTILINEMULTILINE 1

			#MULTILINEMULTILINEMULTILINE 1
				#MULTILINEMULTILINEMULTILINE 1


	#MULTILINEMULTILINEMULTILINE 1
	#MULTILINEMULTILINEMULTILINE 1

	#MULTILINEMULTILINEMULTILINE 1
	#MULTILINEMULTILINEMULTILINE 1
	#MULTILINEMULTILINEMULTILINE 1
		#MULTILINEMULTILINEMULTILINE 1
			#MULTILINEMULTILINEMULTILINE 1

	#MULTILINEMULTILINEMULTILINE 1
#MULTILINEMULTILINEMULTILINE 1


# These are some obsolete codes, we retain them just in case.
"""Multiline Comment2"""

#MULTILINEMULTILINEMULTILINE 2
	#MULTILINEMULTILINEMULTILINE 2
	  #MULTILINEMULTILINEMULTILINE 2
	  #MULTILINEMULTILINEMULTILINE 2

		 #MULTILINEMULTILINEMULTILINE 2

	  #MULTILINEMULTILINEMULTILINE 2
	         #MULTILINEMULTILINEMULTILINE 2
	         #MULTILINEMULTILINEMULTILINE 2
		 #MULTILINEMULTILINEMULTILINE 2
		 #MULTILINEMULTILINEMULTILINE 2
		 #MULTILINEMULTILINEMULTILINE 2

	  #MULTILINEMULTILINEMULTILINE 2
	  	      #MULTILINEMULTILINEMULTILINE 2
		      #MULTILINEMULTILINEMULTILINE 2
		       #MULTILINEMULTILINEMULTILINE 2
	#MULTILINEMULTILINEMULTILINE 2
	#MULTILINEMULTILINEMULTILINE 2
	#MULTILINEMULTILINEMULTILINE 2
	#MULTILINEMULTILINEMULTILINE 2

	#MULTILINEMULTILINEMULTILINE 2
	#MULTILINEMULTILINEMULTILINE 2
		#MULTILINEMULTILINEMULTILINE 2
		#MULTILINEMULTILINEMULTILINE 2
	#MULTILINEMULTILINEMULTILINE 2
		#MULTILINEMULTILINEMULTILINE 2
		#MULTILINEMULTILINEMULTILINE 2


	#MULTILINEMULTILINEMULTILINE 2
	#MULTILINEMULTILINEMULTILINE 2

	#MULTILINEMULTILINEMULTILINE 2
	#MULTILINEMULTILINEMULTILINE 2
	#MULTILINEMULTILINEMULTILINE 2
	#MULTILINEMULTILINEMULTILINE 2
	#MULTILINEMULTILINEMULTILINE 2

	#MULTILINEMULTILINEMULTILINE 2
	#MULTILINEMULTILINEMULTILINE 2
		#MULTILINEMULTILINEMULTILINE 2
			#MULTILINEMULTILINEMULTILINE 2
			#MULTILINEMULTILINEMULTILINE 2
	#MULTILINEMULTILINEMULTILINE 2
		#MULTILINEMULTILINEMULTILINE 2
			#MULTILINEMULTILINEMULTILINE 2
	#MULTILINEMULTILINEMULTILINE 2
	#MULTILINEMULTILINEMULTILINE 2
	#MULTILINEMULTILINEMULTILINE 2

	#MULTILINEMULTILINEMULTILINE 2
	#MULTILINEMULTILINEMULTILINE 2
	#MULTILINEMULTILINEMULTILINE 2
		#MULTILINEMULTILINEMULTILINE 2
			#MULTILINEMULTILINEMULTILINE 2
			#MULTILINEMULTILINEMULTILINE 2
				#MULTILINEMULTILINEMULTILINE 2
	#MULTILINEMULTILINEMULTILINE 2
		#MULTILINEMULTILINEMULTILINE 2
			#MULTILINEMULTILINEMULTILINE 2
			#MULTILINEMULTILINEMULTILINE 2
				#MULTILINEMULTILINEMULTILINE 2

	#MULTILINEMULTILINEMULTILINE 2


#MULTILINEMULTILINEMULTILINE 2

#MULTILINEMULTILINEMULTILINE 2
	#MULTILINEMULTILINEMULTILINE 2
	#MULTILINEMULTILINEMULTILINE 2
	#MULTILINEMULTILINEMULTILINE 2
	#MULTILINEMULTILINEMULTILINE 2
	#MULTILINEMULTILINEMULTILINE 2
		#MULTILINEMULTILINEMULTILINE 2
		#MULTILINEMULTILINEMULTILINE 2
		#MULTILINEMULTILINEMULTILINE 2
			#MULTILINEMULTILINEMULTILINE 2
			#MULTILINEMULTILINEMULTILINE 2
	#MULTILINEMULTILINEMULTILINE 2
		#MULTILINEMULTILINEMULTILINE 2
		#MULTILINEMULTILINEMULTILINE 2
		#MULTILINEMULTILINEMULTILINE 2
			#MULTILINEMULTILINEMULTILINE 2
			#MULTILINEMULTILINEMULTILINE 2
			#MULTILINEMULTILINEMULTILINE 2
			#MULTILINEMULTILINEMULTILINE 2
				#MULTILINEMULTILINEMULTILINE 2
				#MULTILINEMULTILINEMULTILINE 2
				#MULTILINEMULTILINEMULTILINE 2
				#MULTILINEMULTILINEMULTILINE 2
					#MULTILINEMULTILINEMULTILINE 2
					#MULTILINEMULTILINEMULTILINE 2
					#MULTILINEMULTILINEMULTILINE 2
					#MULTILINEMULTILINEMULTILINE 2
					#MULTILINEMULTILINEMULTILINE 2
					#MULTILINEMULTILINEMULTILINE 2
						#MULTILINEMULTILINEMULTILINE 2
						#MULTILINEMULTILINEMULTILINE 2
						#MULTILINEMULTILINEMULTILINE 2
						#MULTILINEMULTILINEMULTILINE 2
						#MULTILINEMULTILINEMULTILINE 2
						#MULTILINEMULTILINEMULTILINE 2
						#MULTILINEMULTILINEMULTILINE 2
						#MULTILINEMULTILINEMULTILINE 2
				#MULTILINEMULTILINEMULTILINE 2
		#MULTILINEMULTILINEMULTILINE 2
			#MULTILINEMULTILINEMULTILINE 2
				#MULTILINEMULTILINEMULTILINE 2
				#MULTILINEMULTILINEMULTILINE 2
				#MULTILINEMULTILINEMULTILINE 2
				#MULTILINEMULTILINEMULTILINE 2
					#MULTILINEMULTILINEMULTILINE 2
					#MULTILINEMULTILINEMULTILINE 2
					#MULTILINEMULTILINEMULTILINE 2
					#MULTILINEMULTILINEMULTILINE 2
					#MULTILINEMULTILINEMULTILINE 2
					#MULTILINEMULTILINEMULTILINE 2

					#MULTILINEMULTILINEMULTILINE 2
						#MULTILINEMULTILINEMULTILINE 2
						#MULTILINEMULTILINEMULTILINE 2
						#MULTILINEMULTILINEMULTILINE 2
						#MULTILINEMULTILINEMULTILINE 2
						#MULTILINEMULTILINEMULTILINE 2
						#MULTILINEMULTILINEMULTILINE 2
						#MULTILINEMULTILINEMULTILINE 2
	#MULTILINEMULTILINEMULTILINE 2
		#MULTILINEMULTILINEMULTILINE 2

#MULTILINEMULTILINEMULTILINE 2
	#MULTILINEMULTILINEMULTILINE 2
	#MULTILINEMULTILINEMULTILINE 2
	#MULTILINEMULTILINEMULTILINE 2
	#MULTILINEMULTILINEMULTILINE 2
	#MULTILINEMULTILINEMULTILINE 2
	#MULTILINEMULTILINEMULTILINE 2
	#MULTILINEMULTILINEMULTILINE 2
		#MULTILINEMULTILINEMULTILINE 2
			#MULTILINEMULTILINEMULTILINE 2
	#MULTILINEMULTILINEMULTILINE 2
		#MULTILINEMULTILINEMULTILINE 2
		#MULTILINEMULTILINEMULTILINE 2
			#MULTILINEMULTILINEMULTILINE 2
				#MULTILINEMULTILINEMULTILINE 2
				#MULTILINEMULTILINEMULTILINE 2
				#MULTILINEMULTILINEMULTILINE 2
				#MULTILINEMULTILINEMULTILINE 2
				#MULTILINEMULTILINEMULTILINE 2
				#MULTILINEMULTILINEMULTILINE 2
	#MULTILINEMULTILINEMULTILINE 2


#MULTILINEMULTILINEMULTILINE 2

	#MULTILINEMULTILINEMULTILINE 2
		#MULTILINEMULTILINEMULTILINE 2
		#MULTILINEMULTILINEMULTILINE 2
		#MULTILINEMULTILINEMULTILINE 2
		#MULTILINEMULTILINEMULTILINE 2

	#MULTILINEMULTILINEMULTILINE 2
		#MULTILINEMULTILINEMULTILINE 2
		#MULTILINEMULTILINEMULTILINE 2
		#MULTILINEMULTILINEMULTILINE 2
		#MULTILINEMULTILINEMULTILINE 2
		#MULTILINEMULTILINEMULTILINE 2

	#MULTILINEMULTILINEMULTILINE 2
		#MULTILINEMULTILINEMULTILINE 2
		#MULTILINEMULTILINEMULTILINE 2
		#MULTILINEMULTILINEMULTILINE 2
		#MULTILINEMULTILINEMULTILINE 2
		#MULTILINEMULTILINEMULTILINE 2

	#MULTILINEMULTILINEMULTILINE 2
		#MULTILINEMULTILINEMULTILINE 2
	        #MULTILINEMULTILINEMULTILINE 2
	        #MULTILINEMULTILINEMULTILINE 2
        	#MULTILINEMULTILINEMULTILINE 2
	        #MULTILINEMULTILINEMULTILINE 2

        	#MULTILINEMULTILINEMULTILINE 2
	        #MULTILINEMULTILINEMULTILINE 2
	        #MULTILINEMULTILINEMULTILINE 2
		#MULTILINEMULTILINEMULTILINE 2

	        #MULTILINEMULTILINEMULTILINE 2
	        	#MULTILINEMULTILINEMULTILINE 2
	        	#MULTILINEMULTILINEMULTILINE 2
	        	#MULTILINEMULTILINEMULTILINE 2
	        	#MULTILINEMULTILINEMULTILINE 2
	        	#MULTILINEMULTILINEMULTILINE 2
	        	#MULTILINEMULTILINEMULTILINE 2
	        	#MULTILINEMULTILINEMULTILINE 2
	        	#MULTILINEMULTILINEMULTILINE 2
	        	#MULTILINEMULTILINEMULTILINE 2
	        		#MULTILINEMULTILINEMULTILINE 2
	        		#MULTILINEMULTILINEMULTILINE 2
	        		#MULTILINEMULTILINEMULTILINE 2
	        		#MULTILINEMULTILINEMULTILINE 2
	        		#MULTILINEMULTILINEMULTILINE 2

	        	#MULTILINEMULTILINEMULTILINE 2

	        #MULTILINEMULTILINEMULTILINE 2

	#MULTILINEMULTILINEMULTILINE 2
		#MULTILINEMULTILINEMULTILINE 2
	        #MULTILINEMULTILINEMULTILINE 2
	        #MULTILINEMULTILINEMULTILINE 2
        	#MULTILINEMULTILINEMULTILINE 2
	        #MULTILINEMULTILINEMULTILINE 2

        	#MULTILINEMULTILINEMULTILINE 2
	        #MULTILINEMULTILINEMULTILINE 2
	        #MULTILINEMULTILINEMULTILINE 2
		#MULTILINEMULTILINEMULTILINE 2

		#MULTILINEMULTILINEMULTILINE 2
	        #MULTILINEMULTILINEMULTILINE 2
	        	#MULTILINEMULTILINEMULTILINE 2
	        	#MULTILINEMULTILINEMULTILINE 2
	        	#MULTILINEMULTILINEMULTILINE 2
	        	#MULTILINEMULTILINEMULTILINE 2
	        	#MULTILINEMULTILINEMULTILINE 2
	        	#MULTILINEMULTILINEMULTILINE 2
	        	#MULTILINEMULTILINEMULTILINE 2
	        	#MULTILINEMULTILINEMULTILINE 2
	        	#MULTILINEMULTILINEMULTILINE 2
	        		#MULTILINEMULTILINEMULTILINE 2
	        		#MULTILINEMULTILINEMULTILINE 2
	        		#MULTILINEMULTILINEMULTILINE 2
	        		#MULTILINEMULTILINEMULTILINE 2
	        		#MULTILINEMULTILINEMULTILINE 2

	        	#MULTILINEMULTILINEMULTILINE 2

		#MULTILINEMULTILINEMULTILINE 2
        	#MULTILINEMULTILINEMULTILINE 2
	        #MULTILINEMULTILINEMULTILINE 2
	        #MULTILINEMULTILINEMULTILINE 2
		#MULTILINEMULTILINEMULTILINE 2

	        #MULTILINEMULTILINEMULTILINE 2
	        	#MULTILINEMULTILINEMULTILINE 2
	        	#MULTILINEMULTILINEMULTILINE 2
	        	#MULTILINEMULTILINEMULTILINE 2
	        	#MULTILINEMULTILINEMULTILINE 2
	        	#MULTILINEMULTILINEMULTILINE 2
	        	#MULTILINEMULTILINEMULTILINE 2
	        	#MULTILINEMULTILINEMULTILINE 2
	        	#MULTILINEMULTILINEMULTILINE 2
	        	#MULTILINEMULTILINEMULTILINE 2
	        		#MULTILINEMULTILINEMULTILINE 2
	        		#MULTILINEMULTILINEMULTILINE 2
	        		#MULTILINEMULTILINEMULTILINE 2
	        		#MULTILINEMULTILINEMULTILINE 2
	        		#MULTILINEMULTILINEMULTILINE 2

	        	#MULTILINEMULTILINEMULTILINE 2

        	#MULTILINEMULTILINEMULTILINE 2
	        #MULTILINEMULTILINEMULTILINE 2
	        #MULTILINEMULTILINEMULTILINE 2
		#MULTILINEMULTILINEMULTILINE 2

	        #MULTILINEMULTILINEMULTILINE 2
	        	#MULTILINEMULTILINEMULTILINE 2
	        	#MULTILINEMULTILINEMULTILINE 2
	        	#MULTILINEMULTILINEMULTILINE 2
	        	#MULTILINEMULTILINEMULTILINE 2
	        	#MULTILINEMULTILINEMULTILINE 2
	        	#MULTILINEMULTILINEMULTILINE 2
	        	#MULTILINEMULTILINEMULTILINE 2
	        	#MULTILINEMULTILINEMULTILINE 2
	        	#MULTILINEMULTILINEMULTILINE 2
	        		#MULTILINEMULTILINEMULTILINE 2
	        		#MULTILINEMULTILINEMULTILINE 2
	        		#MULTILINEMULTILINEMULTILINE 2
	        		#MULTILINEMULTILINEMULTILINE 2
	        		#MULTILINEMULTILINEMULTILINE 2

	        	#MULTILINEMULTILINEMULTILINE 2

	        #MULTILINEMULTILINEMULTILINE 2

	#MULTILINEMULTILINEMULTILINE 2
	#MULTILINEMULTILINEMULTILINE 2
	#MULTILINEMULTILINEMULTILINE 2

	#MULTILINEMULTILINEMULTILINE 2
	#MULTILINEMULTILINEMULTILINE 2
	#MULTILINEMULTILINEMULTILINE 2
	#MULTILINEMULTILINEMULTILINE 2
	#MULTILINEMULTILINEMULTILINE 2
		#MULTILINEMULTILINEMULTILINE 2
		#MULTILINEMULTILINEMULTILINE 2
		#MULTILINEMULTILINEMULTILINE 2
			#MULTILINEMULTILINEMULTILINE 2
			#MULTILINEMULTILINEMULTILINE 2
		#MULTILINEMULTILINEMULTILINE 2
		#MULTILINEMULTILINEMULTILINE 2
		#MULTILINEMULTILINEMULTILINE 2
	#MULTILINEMULTILINEMULTILINE 2
	#MULTILINEMULTILINEMULTILINE 2
	#MULTILINEMULTILINEMULTILINE 2
	#MULTILINEMULTILINEMULTILINE 2
	#MULTILINEMULTILINEMULTILINE 2
	#MULTILINEMULTILINEMULTILINE 2

	#MULTILINEMULTILINEMULTILINE 2
	#MULTILINEMULTILINEMULTILINE 2
	#MULTILINEMULTILINEMULTILINE 2
	#MULTILINEMULTILINEMULTILINE 2
		#MULTILINEMULTILINEMULTILINE 2
		#MULTILINEMULTILINEMULTILINE 2
		#MULTILINEMULTILINEMULTILINE 2
		#MULTILINEMULTILINEMULTILINE 2

	#MULTILINEMULTILINEMULTILINE 2
	#MULTILINEMULTILINEMULTILINE 2
	#MULTILINEMULTILINEMULTILINE 2
		#MULTILINEMULTILINEMULTILINE 2
		#MULTILINEMULTILINEMULTILINE 2

	#MULTILINEMULTILINEMULTILINE 2
	#MULTILINEMULTILINEMULTILINE 2
	#MULTILINEMULTILINEMULTILINE 2

	#MULTILINEMULTILINEMULTILINE 2
	#MULTILINEMULTILINEMULTILINE 2

	#MULTILINEMULTILINEMULTILINE 2
		#MULTILINEMULTILINEMULTILINE 2
			#MULTILINEMULTILINEMULTILINE 2
			#MULTILINEMULTILINEMULTILINE 2
				#MULTILINEMULTILINEMULTILINE 2
				#MULTILINEMULTILINEMULTILINE 2
			#MULTILINEMULTILINEMULTILINE 2
	#MULTILINEMULTILINEMULTILINE 2
	#MULTILINEMULTILINEMULTILINE 2
		#MULTILINEMULTILINEMULTILINE 2
	#MULTILINEMULTILINEMULTILINE 2

	#MULTILINEMULTILINEMULTILINE 2


#MULTILINEMULTILINEMULTILINE 2
	#MULTILINEMULTILINEMULTILINE 2
	#MULTILINEMULTILINEMULTILINE 2
	#MULTILINEMULTILINEMULTILINE 2
	#MULTILINEMULTILINEMULTILINE 2
	#MULTILINEMULTILINEMULTILINE 2
	#MULTILINEMULTILINEMULTILINE 2

	#MULTILINEMULTILINEMULTILINE 2
	#MULTILINEMULTILINEMULTILINE 2
	#MULTILINEMULTILINEMULTILINE 2
	#MULTILINEMULTILINEMULTILINE 2
	#MULTILINEMULTILINEMULTILINE 2
	#MULTILINEMULTILINEMULTILINE 2

	#MULTILINEMULTILINEMULTILINE 2
	#MULTILINEMULTILINEMULTILINE 2
	#MULTILINEMULTILINEMULTILINE 2
	#MULTILINEMULTILINEMULTILINE 2
	#MULTILINEMULTILINEMULTILINE 2
		#MULTILINEMULTILINEMULTILINE 2
		#MULTILINEMULTILINEMULTILINE 2
			#MULTILINEMULTILINEMULTILINE 2
			#MULTILINEMULTILINEMULTILINE 2
		#MULTILINEMULTILINEMULTILINE 2

	#MULTILINEMULTILINEMULTILINE 2
	#MULTILINEMULTILINEMULTILINE 2
		#MULTILINEMULTILINEMULTILINE 2
	#MULTILINEMULTILINEMULTILINE 2
		#MULTILINEMULTILINEMULTILINE 2

	#MULTILINEMULTILINEMULTILINE 2
	#MULTILINEMULTILINEMULTILINE 2
		#MULTILINEMULTILINEMULTILINE 2
	#MULTILINEMULTILINEMULTILINE 2
		#MULTILINEMULTILINEMULTILINE 2

	#MULTILINEMULTILINEMULTILINE 2
	#MULTILINEMULTILINEMULTILINE 2

	#MULTILINEMULTILINEMULTILINE 2
	#MULTILINEMULTILINEMULTILINE 2
	#MULTILINEMULTILINEMULTILINE 2
		#MULTILINEMULTILINEMULTILINE 2
		#MULTILINEMULTILINEMULTILINE 2

		#MULTILINEMULTILINEMULTILINE 2

		#MULTILINEMULTILINEMULTILINE 2
			#MULTILINEMULTILINEMULTILINE 2
			#MULTILINEMULTILINEMULTILINE 2
			#MULTILINEMULTILINEMULTILINE 2

	#MULTILINEMULTILINEMULTILINE 2



#MULTILINEMULTILINEMULTILINE 2
	#MULTILINEMULTILINEMULTILINE 2
	  #MULTILINEMULTILINEMULTILINE 2
	  #MULTILINEMULTILINEMULTILINE 2
	  #MULTILINEMULTILINEMULTILINE 2
		#MULTILINEMULTILINEMULTILINE 2
		#MULTILINEMULTILINEMULTILINE 2
	  #MULTILINEMULTILINEMULTILINE 2
	#MULTILINEMULTILINEMULTILINE 2
	#MULTILINEMULTILINEMULTILINE 2
	#MULTILINEMULTILINEMULTILINE 2
	#MULTILINEMULTILINEMULTILINE 2
	#MULTILINEMULTILINEMULTILINE 2
		#MULTILINEMULTILINEMULTILINE 2
			#MULTILINEMULTILINEMULTILINE 2
			#MULTILINEMULTILINEMULTILINE 2
			#MULTILINEMULTILINEMULTILINE 2
			#MULTILINEMULTILINEMULTILINE 2
	#MULTILINEMULTILINEMULTILINE 2

#MULTILINEMULTILINEMULTILINE 2
	#MULTILINEMULTILINEMULTILINE 2
	  #MULTILINEMULTILINEMULTILINE 2
	  #MULTILINEMULTILINEMULTILINE 2
	  #MULTILINEMULTILINEMULTILINE 2
		#MULTILINEMULTILINEMULTILINE 2
		#MULTILINEMULTILINEMULTILINE 2
	  #MULTILINEMULTILINEMULTILINE 2
	#MULTILINEMULTILINEMULTILINE 2
	#MULTILINEMULTILINEMULTILINE 2
	#MULTILINEMULTILINEMULTILINE 2
	#MULTILINEMULTILINEMULTILINE 2
	#MULTILINEMULTILINEMULTILINE 2
		#MULTILINEMULTILINEMULTILINE 2
			#MULTILINEMULTILINEMULTILINE 2
			#MULTILINEMULTILINEMULTILINE 2
			#MULTILINEMULTILINEMULTILINE 2
			#MULTILINEMULTILINEMULTILINE 2
	#MULTILINEMULTILINEMULTILINE 2
#MULTILINEMULTILINEMULTILINE 2


from builtins import range
from sp_global_def import *
