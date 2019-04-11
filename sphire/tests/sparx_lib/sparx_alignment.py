#
from __future__ import print_function
# Author: Pawel A.Penczek, 09/09/2006 (Pawel.A.Penczek@uth.tmc.edu)
# Copyright (c) 2000-2006 The University of Texas - Houston Medical School
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
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
#

import EMAN2_cppwrap
import sparx_applications
import sparx_filter
import sparx_fundamentals
import sparx_global_def
import sparx_logger
import sparx_morphology
import mpi
import sparx_multi_shc
import numpy
import numpy as np
import numpy.random
import operator
import sparx_pixel_error
import sparx_projection
import scipy
import sparx_statistics
import sys
import time
import types
import sparx_utilities
pass#IMPORTIMPORTIMPORT import EMAN2
pass#IMPORTIMPORTIMPORT import EMAN2_cppwrap
pass#IMPORTIMPORTIMPORT import alignment
pass#IMPORTIMPORTIMPORT import applications
pass#IMPORTIMPORTIMPORT import cPickle as pickle
pass#IMPORTIMPORTIMPORT import collections
pass#IMPORTIMPORTIMPORT import development
pass#IMPORTIMPORTIMPORT import filter
pass#IMPORTIMPORTIMPORT import fundamentals
pass#IMPORTIMPORTIMPORT import gc
pass#IMPORTIMPORTIMPORT import glob
pass#IMPORTIMPORTIMPORT import global_def
pass#IMPORTIMPORTIMPORT import inspect
pass#IMPORTIMPORTIMPORT import logger
pass#IMPORTIMPORTIMPORT import math
pass#IMPORTIMPORTIMPORT import morphology
pass#IMPORTIMPORTIMPORT import mpi
pass#IMPORTIMPORTIMPORT import multi_shc
pass#IMPORTIMPORTIMPORT import numpy
pass#IMPORTIMPORTIMPORT import numpy as np
pass#IMPORTIMPORTIMPORT import numpy.random
pass#IMPORTIMPORTIMPORT import operator
pass#IMPORTIMPORTIMPORT import os
pass#IMPORTIMPORTIMPORT import pixel_error
pass#IMPORTIMPORTIMPORT import projection
pass#IMPORTIMPORTIMPORT import random
pass#IMPORTIMPORTIMPORT import scipy
pass#IMPORTIMPORTIMPORT import socket
pass#IMPORTIMPORTIMPORT import sparx
pass#IMPORTIMPORTIMPORT import statistics
pass#IMPORTIMPORTIMPORT import subprocess
pass#IMPORTIMPORTIMPORT import sys
pass#IMPORTIMPORTIMPORT import time
pass#IMPORTIMPORTIMPORT import types
pass#IMPORTIMPORTIMPORT import utilities
from builtins import range
pass#IMPORTIMPORTIMPORT from global_def import *
pass#IMPORTIMPORTIMPORT import numpy.random

#  06-12-14 code lifted
def ali2d_single_iter(data, numr, wr, cs, tavg, cnx, cny, \
						xrng, yrng, step, nomirror = False, mode="F", CTF=False, \
						random_method="", T=1.0, ali_params="xform.align2d", delta = 0.0):
	"""
		single iteration of 2D alignment using ormq
		if CTF = True, apply CTF to data (not to reference!)
	"""
	pass#IMPORTIMPORTIMPORT from utilities import combine_params2, inverse_transform2, get_params2D, set_params2D
	pass#IMPORTIMPORTIMPORT from alignment import ormq, ornq

	if CTF:
		pass#IMPORTIMPORTIMPORT from filter  import filt_ctf

	maxrin = numr[-1]  #  length
	ou = numr[-3]  #  maximum radius
	if random_method == "SCF":
		pass#IMPORTIMPORTIMPORT from fundamentals import fft, scf
		pass#IMPORTIMPORTIMPORT from alignment import multalign2d_scf
		frotim = [sparx_fundamentals.fft(tavg)]
		xrng = int(xrng+0.5)
		yrng = int(yrng+0.5)
		cimage = EMAN2_cppwrap.Util.Polar2Dm(sparx_fundamentals.scf(tavg), cnx, cny, numr, mode)
		EMAN2_cppwrap.Util.Frngs(cimage, numr)
		EMAN2_cppwrap.Util.Applyws(cimage, numr, wr)
	else:
		# 2D alignment using rotational ccf in polar coords and quadratic interpolation
		cimage = EMAN2_cppwrap.Util.Polar2Dm(tavg, cnx, cny, numr, mode)
		EMAN2_cppwrap.Util.Frngs(cimage, numr)
		EMAN2_cppwrap.Util.Applyws(cimage, numr, wr)

	sx_sum = 0.0
	sy_sum = 0.0
	sxn = 0.
	syn = 0.
	mn = 0
	nope = 0
	mashi = cnx-ou-2
	for im in range(len(data)):
		if CTF:
			#Apply CTF to image
			ctf_params = data[im].get_attr("ctf")
			ima = sparx_filter.filt_ctf(data[im], ctf_params, True)
		else:
			ima = data[im]

		if( random_method == "PCP"):
			sxi = data[im][0][0].get_attr('sxi')
			syi = data[im][0][0].get_attr('syi')
			nx = ny = data[im][0][0].get_attr('inx')
		else:
			nx = ima.get_xsize()
			ny = ima.get_ysize()
			alpha, sx, sy, mirror, dummy = sparx_utilities.get_params2D(data[im], ali_params)
			alpha, sx, sy, dummy         = sparx_utilities.combine_params2(alpha, sx, sy, mirror, 0.0, -cs[0], -cs[1], 0)
			alphai, sxi, syi, scalei     = sparx_utilities.inverse_transform2(alpha, sx, sy)
			#  introduce constraints on parameters to accomodate use of cs centering
			sxi = min(max(sxi,-mashi),mashi)
			syi = min(max(syi,-mashi),mashi)

		#  The search range procedure was adjusted for 3D searches, so since in 2D the order of operations is inverted, we have to invert ranges
		txrng = search_range(nx, ou, sxi, xrng, "ali2d_single_iter")
		txrng = [txrng[1],txrng[0]]
		tyrng = search_range(ny, ou, syi, yrng, "ali2d_single_iter")
		tyrng = [tyrng[1],tyrng[0]]
		#print im, "B",cnx,sxi,syi,txrng, tyrng
		# align current image to the reference
		if random_method == "SHC":
			"""Multiline Comment0"""
			#  For shc combining of shifts is problematic as the image may randomly slide away and never come back.
			#  A possibility would be to reject moves that results in too large departure from the center.
			#  On the other hand, one cannot simply do searches around the proper center all the time,
			#    as if xr is decreased, the image cannot be brought back if the established shifts are further than new range
			olo = EMAN2_cppwrap.Util.shc(ima, [cimage], txrng, tyrng, step, -1.0, mode, numr, cnx+sxi, cny+syi, "c1")
			##olo = Util.shc(ima, [cimage], xrng, yrng, step, -1.0, mode, numr, cnx, cny, "c1")
			if(data[im].get_attr("previousmax")<olo[5]):
				#[angt, sxst, syst, mirrort, peakt] = ormq(ima, cimage, xrng, yrng, step, mode, numr, cnx+sxi, cny+syi, delta)
				#print  angt, sxst, syst, mirrort, peakt,olo
				angt = olo[0]
				sxst = olo[1]
				syst = olo[2]
				mirrort = int(olo[3])
				# combine parameters and set them to the header, ignore previous angle and mirror
				[alphan, sxn, syn, mn] = sparx_utilities.combine_params2(0.0, -sxi, -syi, 0, angt, sxst, syst, mirrort)
				sparx_utilities.set_params2D(data[im], [alphan, sxn, syn, mn, 1.0], ali_params)
				##set_params2D(data[im], [angt, sxst, syst, mirrort, 1.0], ali_params)
				data[im].set_attr("previousmax",olo[5])
			else:
				# Did not find a better peak, but we have to set shifted parameters, as the average shifted
				sparx_utilities.set_params2D(data[im], [alpha, sx, sy, mirror, 1.0], ali_params)
				nope += 1
				mn = 0
				sxn = 0.0
				syn = 0.0
		elif random_method == "SCF":
			sxst,syst,iref,angt,mirrort,totpeak = multalign2d_scf(data[im], [cimage], frotim, numr, xrng, yrng, ou = ou)
			[alphan, sxn, syn, mn] = sparx_utilities.combine_params2(0.0, -sxi, -syi, 0, angt, sxst, syst, mirrort)
			sparx_utilities.set_params2D(data[im], [alphan, sxn, syn, mn, 1.0], ali_params)
		elif random_method == "PCP":
			[angt, sxst, syst, mirrort, peakt] = ormq_fast(data[im], cimage, txrng, tyrng, step, numr, mode, delta)
			sxst = rings[0][0][0].get_attr("sxi")
			syst = rings[0][0][0].get_attr("syi")
			print(sxst, syst,sx,sy)
			dummy,sxs,sys, dummy = sparx_utilities.inverse_transform2(-angt,sx+sxst,sy+syst)
			sparx_utilities.set_params2D(data[im][0][0], [angt, sxs, sys, mirrort, 1.0], ali_params)
		else:
			if nomirror:  [angt, sxst, syst, mirrort, peakt] = ornq(ima, cimage, txrng, tyrng, step, mode, numr, cnx+sxi, cny+syi)
			else:	      [angt, sxst, syst, mirrort, peakt] = ormq(ima, cimage, txrng, tyrng, step, mode, numr, cnx+sxi, cny+syi, delta)
			# combine parameters and set them to the header, ignore previous angle and mirror
			[alphan, sxn, syn, mn] = sparx_utilities.combine_params2(0.0, -sxi, -syi, 0, angt, sxst, syst, mirrort)
			sparx_utilities.set_params2D(data[im], [alphan, sxn, syn, mn, 1.0], ali_params)

		if mn == 0: sx_sum += sxn
		else:       sx_sum -= sxn
		sy_sum += syn

	return sx_sum, sy_sum, nope


"""Multiline Comment1"""

def ang_n(tot, mode, maxrin):
	"""
	  Calculate angle based on the position of the peak
	"""
	pass#IMPORTIMPORTIMPORT from math import fmod
	if (mode == 'f' or mode == 'F'): return numpy.fmod(((tot-1.0) / maxrin+1.0)*360.0, 360.0)
	else:                            return numpy.fmod(((tot-1.0) / maxrin+1.0)*180.0, 180.0)

# Copy of this function is implemented in C++ in Util (Util.Applyws). It works much faster than this one.
"""Multiline Comment2"""

def log2(n):
	""" Returns the smallest power by which 2 has to be raised to obtain
	    an integer less equal n
	"""
	m = 1
	k =-1
	while (m <= n):
		i = m
		k +=1
		m = 2*i
	return k
    
def Numrinit(first_ring, last_ring, skip=1, mode="F"):
	"""This function calculates the necessary information for the 2D 
	   polar interpolation. For each ring, three elements are recorded:
	   numr[i*3]:  Radius of this ring
	   numr[i*3+1]: Total number of samples of all inner rings+1
	   		(Or, the beginning point of this ring)
	   numr[i*3+2]: Number of samples of this ring. This number is an 
	   		FFT-friendly power of the 2.
			
	   "F" means a full circle interpolation
	   "H" means a half circle interpolation
	"""
	MAXFFT = 32768
	pass#IMPORTIMPORTIMPORT from math import pi

	if (mode == 'f' or mode == 'F'): dpi = 2*numpy.pi
	else:                            dpi = numpy.pi
	numr = []
	lcirc = 1
	for k in range(first_ring, last_ring+1, skip):
		numr.append(k)
		jp = int(dpi * k+0.5)
		ip = 2**(log2(jp)+1)  # two times oversample each ring
		if (k+skip <= last_ring and jp > ip+ip//2): ip=min(MAXFFT,2*ip)
		if (k+skip  > last_ring and jp > ip+ip//5): ip=min(MAXFFT,2*ip)

		numr.append(lcirc)
		numr.append(ip)
		lcirc += ip

	return  numr
"""Multiline Comment3"""

def ringwe(numr, mode="F"):
	"""
	   Calculate ring weights for rotational alignment
	   The weights are r*delta(r)*delta(phi).
	"""
	pass#IMPORTIMPORTIMPORT from math import pi
	if (mode == 'f' or mode == 'F'):
		dpi = 2*numpy.pi
	else:
		dpi = numpy.pi
	nring = len(numr)/3
	wr=[0.0]*nring
	maxrin = float(numr[len(numr)-1])
	for i in range(0,nring): wr[i] = numr[i*3]*dpi/float(numr[2+i*3])*maxrin/float(numr[2+i*3])
	return wr

def ornq(image, crefim, xrng, yrng, step, mode, numr, cnx, cny, deltapsi = 0.0):
	"""Determine shift and rotation between image and reference image (refim)
	   no mirror
		quadratic interpolation
		cnx, cny in FORTRAN convention
	"""
	pass#IMPORTIMPORTIMPORT from math import pi, cos, sin, radians
	pass#IMPORTIMPORTIMPORT from alignment import ang_n
	#from utilities import info
	#print "ORNQ"
	peak = -1.0E23

	lkx = int(xrng[0]/step)
	rkx = int(xrng[-1]/step)

	lky = int(yrng[0]/step)
	rky = int(yrng[-1]/step)

	for i in range(-lky, rky+1):
		iy = i*step
		for j in range(-lkx, rkx+1):
			ix = j*step
			cimage = EMAN2_cppwrap.Util.Polar2Dm(image, cnx+ix, cny+iy, numr, mode)
			EMAN2_cppwrap.Util.Frngs(cimage, numr)
			retvals = EMAN2_cppwrap.Util.Crosrng_e(crefim, cimage, numr, 0, deltapsi)
			qn = retvals["qn"]
			if qn >= peak:
				sx = -ix
				sy = -iy
				ang = ang_n(retvals["tot"], mode, numr[-1])
				peak = qn
	# mirror is returned as zero for consistency
	mirror = 0
	co =  numpy.cos(numpy.radians(ang))
	so = -numpy.sin(numpy.radians(ang))
	sxs = sx*co - sy*so
	sys = sx*so + sy*co
	return  ang, sxs, sys, mirror, peak


def ormq(image, crefim, xrng, yrng, step, mode, numr, cnx, cny, delta = 0.0):
	"""Determine shift and rotation between image and reference image (crefim)
		crefim should be as FT of polar coords with applied weights
	        consider mirror
		quadratic interpolation
		cnx, cny in FORTRAN convention
	"""
	pass#IMPORTIMPORTIMPORT from math import pi, cos, sin, radians
	#print "ORMQ"
	peak = -1.0E23

	lkx = int(xrng[0]/step)
	rkx = int(xrng[-1]/step)

	lky = int(yrng[0]/step)
	rky = int(yrng[-1]/step)

	for i in range(-lky, rky+1):
		iy = i*step
		for j in range(-lkx, rkx+1):
			ix = j*step
			cimage = EMAN2_cppwrap.Util.Polar2Dm(image, cnx+ix, cny+iy, numr, mode)
			EMAN2_cppwrap.Util.Frngs(cimage, numr)
			# The following code it used when mirror is considered
			if delta == 0.0:
				retvals = EMAN2_cppwrap.Util.Crosrng_ms(crefim, cimage, numr, 0.0)
			else:
				retvals = EMAN2_cppwrap.Util.Crosrng_ms_delta(crefim, cimage, numr, 0.0, delta)
			qn = retvals["qn"]
			qm = retvals["qm"]
			if (qn >= peak or qm >= peak):
				sx = -ix
				sy = -iy
				if (qn >= qm):
					ang = ang_n(retvals["tot"], mode, numr[-1])
					peak = qn
					mirror = 0
				else:
					ang = ang_n(retvals["tmt"], mode, numr[-1])
					peak = qm
					mirror = 1
			"""Multiline Comment4"""
	co  =  numpy.cos(numpy.radians(ang))
	so  = -numpy.sin(numpy.radians(ang))
	sxs = sx*co - sy*so
	sys = sx*so + sy*co
	return  ang, sxs, sys, mirror, peak

def ormq_fast(dimage, crefim, xrng, yrng, step, numr, mode, delta = 0.0):
	"""Determine shift and rotation between image and reference image (crefim)
		crefim should be as FT of polar coords with applied weights
	        consider mirror
		cnx, cny in FORTRAN convention
	"""
	#from math import pi, cos, sin, radians
	#print "ORMQ_FAST"
	maxrange = len(dimage)//2
	#istep = int(2*step)
	istep = int(step)
	"""Multiline Comment5"""
	lkx = rkx = int(xrng*istep)

	lky = rky = int(yrng*istep)

	peak = -1.0E23
	for j in range(-lky, rky+1, istep):
		for i in range(-lkx, rkx+1, istep):
			if delta == 0.0: retvals = EMAN2_cppwrap.Util.Crosrng_ms(crefim, dimage[i+maxrange][j+maxrange], numr, 0.0)
			else:            retvals = EMAN2_cppwrap.Util.Crosrng_ms_delta(crefim, dimage[i+maxrange][j+maxrange], numr, delta)
			qn = retvals["qn"]
			qm = retvals["qm"]
			if (qn >= peak or qm >= peak):
				sx = i
				sy = j
				if (qn >= qm):
					ang = ang_n(retvals["tot"], mode, numr[-1])
					peak = qn
					mirror = 0
				else:
					ang = ang_n(retvals["tmt"], mode, numr[-1])
					peak = qm
					mirror = 1
	"""Multiline Comment6"""
	if( peak < -1.0e20): sparx_global_def.ERROR("ormq_fast","failed, most likely due to search ranges",1)
	#return  ang, sx/2.0, sy/2.0, mirror, peak
	return  ang, sx, sy, mirror, peak
			

def prepref(data, maskfile, cnx, cny, numr, mode, maxrangex, maxrangey, step):
	pass#IMPORTIMPORTIMPORT from utilities import get_params2D, combine_params2
	pass#IMPORTIMPORTIMPORT from EMAN2 import Util
	#step = 1
	mashi = cnx -numr[-3] -2
	nima = len(data)
	istep = int(1.0/step)
	dimage = [[[None for j in range(2*maxrangey*istep+1)] for i in range(2*maxrangex*istep+1)] for im in range(nima) ]
	for im in range(nima):
		sts = EMAN2_cppwrap.Util.infomask(data[im], maskfile, False)
		data[im] -= sts[0]
		data[im] /= sts[1]
		alpha, sx, sy, mirror, dummy = sparx_utilities.get_params2D(data[im])
		#alpha, sx, sy, dummy         = combine_params2(alpha, sx, sy, mirror, 0.0, -cs[0], -cs[1], 0)
		alphai, sxi, syi, dummy      = sparx_utilities.combine_params2(0.0, sx, sy, 0, -alpha, 0,0, 0)
		#  introduce constraints on parameters to accomodate use of cs centering
		sxi = min(max(sxi,-mashi),mashi)
		syi = min(max(syi,-mashi),mashi)	
		for j in range(-maxrangey*istep, maxrangey*istep+1):
			iy = j*step
			for i in range(-maxrangex*istep, maxrangex*istep+1):
				ix = i*step
				dimage[im][i+maxrangex][j+maxrangey] = EMAN2_cppwrap.Util.Polar2Dm(data[im], cnx+sxi+ix, cny+syi+iy, numr, mode)
				#print ' prepref  ',j,i,j+maxrangey,i+maxrangex
				EMAN2_cppwrap.Util.Frngs(dimage[im][i+maxrangex][j+maxrangey], numr)
		dimage[im][0][0].set_attr("sxi",sxi)
		dimage[im][0][0].set_attr("syi",syi)

	return dimage

def prepare_refrings( volft, kb, nz = -1, delta = 2.0, ref_a = "P", sym = "c1", numr = None, MPI=False, \
						phiEqpsi = "Zero", kbx = None, kby = None, initial_theta = None, \
						delta_theta = None, initial_phi = None):
	"""
		Generate quasi-evenly distributed reference projections converted to rings
		ref_a can be a list of angles, in which case it is used instead of being generated
	"""
	pass#IMPORTIMPORTIMPORT from projection   import prep_vol, prgs
	pass#IMPORTIMPORTIMPORT from applications import MPI_start_end
	pass#IMPORTIMPORTIMPORT from utilities    import even_angles, getfvec
	pass#IMPORTIMPORTIMPORT from types        import BooleanType



	# mpi communicator can be sent by the MPI parameter
	if type(MPI) is types.BooleanType:
		if MPI:
			pass#IMPORTIMPORTIMPORT from mpi import MPI_COMM_WORLD
			mpi_comm = mpi.MPI_COMM_WORLD
	else:
		mpi_comm = MPI
		MPI = True

	mode = "F"

	pass#IMPORTIMPORTIMPORT from types import ListType
	if(type(ref_a) is types.ListType):
		# if ref_a is  list, it has to be a list of projection directions, use it
		ref_angles = ref_a
	else:
		# generate list of Eulerian angles for reference projections
		#  phi, theta, psi
		if initial_theta and initial_phi :
			ref_angles = sparx_utilities.even_angles(delta, theta1 = initial_theta, phi1 = initial_phi, symmetry=sym, method = ref_a, phiEqpsi = phiEqpsi)
		else:
			if initial_theta is None:
				if(sym[:1] == "c" or sym[:1] == "d"):
					ref_angles = sparx_utilities.even_angles(delta, symmetry=sym, method = ref_a, phiEqpsi = phiEqpsi)
				else:
					pass#IMPORTIMPORTIMPORT from fundamentals import symclass
					psp = sparx_fundamentals.symclass(sym)
					ref_angles = psp.even_angles(delta)
					del psp
			else:
				if delta_theta is None: delta_theta = 1.0
				ref_angles = sparx_utilities.even_angles(delta, theta1 = initial_theta, theta2 = delta_theta, symmetry=sym, method = ref_a, phiEqpsi = phiEqpsi)


	wr_four  = ringwe(numr, mode)
	cnx = nz//2 + 1
	cny = nz//2 + 1
	num_ref = len(ref_angles)

	if MPI:
		pass#IMPORTIMPORTIMPORT from mpi import mpi_comm_rank, mpi_comm_size
		myid = mpi.mpi_comm_rank( mpi_comm )
		ncpu = mpi.mpi_comm_size( mpi_comm )
	else:
		ncpu = 1
		myid = 0

	if(nz <1):  sparx_global_def.ERROR("Data size has to be given (nz)", "prepare_refrings", 1, myid)
	
	ref_start, ref_end = sparx_applications.MPI_start_end(num_ref, ncpu, myid)

	refrings = []     # list of (image objects) reference projections in Fourier representation

	sizex = numr[len(numr)-2] + numr[len(numr)-1]-1

	for i in range(num_ref):
		prjref = EMAN2_cppwrap.EMData()
		prjref.set_size(sizex, 1, 1)
		refrings.append(prjref)

	if kbx is None:
		for i in range(ref_start, ref_end):
			prjref = sparx_projection.prgs(volft, kb, [ref_angles[i][0], ref_angles[i][1], ref_angles[i][2], 0.0, 0.0])
			cimage = EMAN2_cppwrap.Util.Polar2Dm(prjref, cnx, cny, numr, mode)  # currently set to quadratic....
			EMAN2_cppwrap.Util.Normalize_ring(cimage, numr, 0 )
			EMAN2_cppwrap.Util.Frngs(cimage, numr)
			EMAN2_cppwrap.Util.Applyws(cimage, numr, wr_four)
			refrings[i] = cimage
	else:
		for i in range(ref_start, ref_end):
			prjref = sparx_projection.prgs(volft, kb, [ref_angles[i][0], ref_angles[i][1], ref_angles[i][2], 0.0, 0.0], kbx, kby)
			cimage = EMAN2_cppwrap.Util.Polar2Dm(prjref, cnx, cny, numr, mode)  # currently set to quadratic....
			EMAN2_cppwrap.Util.Normalize_ring(cimage, numr, 0 )
			EMAN2_cppwrap.Util.Frngs(cimage, numr)
			EMAN2_cppwrap.Util.Applyws(cimage, numr, wr_four)
			refrings[i] = cimage

	if MPI:
		pass#IMPORTIMPORTIMPORT from utilities import bcast_compacted_EMData_all_to_all
		sparx_utilities.bcast_compacted_EMData_all_to_all(refrings, myid, comm=mpi_comm)

	for i in range(len(ref_angles)):
		n1,n2,n3 = sparx_utilities.getfvec(ref_angles[i][0], ref_angles[i][1])
		refrings[i].set_attr_dict( {"phi":ref_angles[i][0], "theta":ref_angles[i][1], "psi":ref_angles[i][2], "n1":n1, "n2":n2, "n3":n3} )

	return refrings

def proj_ali_incore(data, refrings, numr, xrng, yrng, step, finfo=None, sym = "c1", delta_psi = 0.0, rshift = 0.0):
	pass#IMPORTIMPORTIMPORT from alignment import search_range
	pass#IMPORTIMPORTIMPORT from EMAN2 import Vec2f

	if finfo:
		pass#IMPORTIMPORTIMPORT from utilities    import get_params_proj
		phi, theta, psi, s2x, s2y = sparx_utilities.get_params_proj(data)
		finfo.write("Old parameters: %9.4f %9.4f %9.4f %9.4f %9.4f\n"%(phi, theta, psi, s2x, s2y))
		finfo.flush()

	mode = "F"
	#  center is in SPIDER convention
	nx   = data.get_xsize()
	ny   = data.get_ysize()
	cnx  = nx//2 + 1
	cny  = ny//2 + 1

	#phi, theta, psi, sxo, syo = get_params_proj(data)
	t1 = data.get_attr("xform.projection")
	dp = t1.get_params("spider")
	ou = numr[-3]
	sxi = round(-dp["tx"]+rshift,2)
	syi = round(-dp["ty"]+rshift,2)
	txrng = search_range(nx, ou, sxi, xrng)
	tyrng = search_range(ny, ou, syi, yrng)

	[ang, sxs, sys, mirror, iref, peak] = EMAN2_cppwrap.Util.multiref_polar_ali_3d(data, refrings, txrng, tyrng, step, mode, numr, cnx-sxi, cny-syi, delta_psi)
	#print ang, sxs, sys, mirror, iref, peak
	iref = int(iref)
	#  What that means is that one has to change the the Eulerian angles so they point into mirrored direction: phi+180, 180-theta, 180-psi
	#  rotation has to be reversed
	if mirror:
		phi   = (refrings[iref].get_attr("phi")+540.0)%360.0
		theta = 180.0-refrings[iref].get_attr("theta")
		psi   = (540.0-refrings[iref].get_attr("psi")-ang)%360.0
	else:
		phi   = refrings[iref].get_attr("phi")
		theta = refrings[iref].get_attr("theta")
		psi   = (360.0+refrings[iref].get_attr("psi")-ang)%360.0
	s2x   = sxs + sxi
	s2y   = sys + syi
	#set_params_proj(data, [phi, theta, psi, s2x, s2y])
	t2 = EMAN2_cppwrap.Transform({"type":"spider","phi":phi,"theta":theta,"psi":psi})
	t2.set_trans(EMAN2_cppwrap.Vec2f(-s2x, -s2y))
	data.set_attr("xform.projection", t2)
	data.set_attr("referencenumber", iref)
	pass#IMPORTIMPORTIMPORT from pixel_error import max_3D_pixel_error
	ts = t2.get_sym_proj(sym)
	if(len(ts) > 1):
		# only do it if it is not c1
		pixel_error = +1.0e23
		for ut in ts:
			# we do not care which position minimizes the error
			pixel_error = min(sparx_pixel_error.max_3D_pixel_error(t1, ut, numr[-3]), pixel_error)
	else:
		pixel_error = sparx_pixel_error.max_3D_pixel_error(t1, t2, numr[-3])
	

	if finfo:
		finfo.write( "New parameters: %9.4f %9.4f %9.4f %9.4f %9.4f %10.5f  %11.3e\n\n" %(phi, theta, psi, s2x, s2y, peak, pixel_error))
		finfo.flush()

	return peak, pixel_error

def proj_ali_incore_local(data, refrings, list_of_reference_angles, numr, xrng, yrng, step, an, finfo=None, sym='c1', delta_psi = 0.0, rshift = 0.0):
	pass#IMPORTIMPORTIMPORT from alignment    import search_range
	#from utilities    import set_params_proj, get_params_proj
	pass#IMPORTIMPORTIMPORT from math         import cos, sin, pi, radians
	pass#IMPORTIMPORTIMPORT from EMAN2        import Vec2f

	mode = "F"
	nx   = data.get_xsize()
	ny   = data.get_ysize()
	#  center is in SPIDER convention
	cnx  = nx//2 + 1
	cny  = ny//2 + 1

	ant = numpy.cos(numpy.radians(an))

	return True

	#phi, theta, psi, sxo, syo = get_params_proj(data)
	t1 = data.get_attr("xform.projection")
	dp = t1.get_params("spider")
	ou = numr[-3]
	sxi = round(-dp["tx"] + rshift, 2)
	syi = round(-dp["ty"] + rshift, 2)
	txrng = search_range(nx, ou, sxi, xrng)
	tyrng = search_range(ny, ou, syi, yrng)
	if finfo:
		finfo.write("Old parameters: %6.2f %6.2f %6.2f %6.2f %6.2f\n"%(dp["phi"], dp["theta"], dp["psi"], -dp["tx"], -dp["ty"]))
		finfo.write("ou, nx, ny, xrng, yrng, cnx, cny, sxi, syi, txrng[0],txrng[1],tyrng[0],tyrng[1] : %3d  %3d  %3d    %4.1f  %4.1f %3d %3d   %4.1f  %4.1f     %4.1f  %4.1f %4.1f %4.1f\n"%(ou, nx, ny, xrng, yrng, cnx, cny, sxi, syi, txrng[0],txrng[1],tyrng[0],tyrng[1]))
		finfo.flush()

	[ang, sxs, sys, mirror, iref, peak] = EMAN2_cppwrap.Util.multiref_polar_ali_3d_local(data, refrings, list_of_reference_angles, txrng, tyrng, step, ant, mode, numr, cnx-sxi, cny-syi, sym, delta_psi)

	iref=int(iref)
	if iref > -1:
		# What that means is that one has to change the the Eulerian angles so they point into mirrored direction: phi+180, 180-theta, 180-psi
		if  mirror:
			phi   = (refrings[iref].get_attr("phi")+540.0)%360.0
			theta = 180.0-refrings[iref].get_attr("theta")
			psi   = (540.0-refrings[iref].get_attr("psi")-ang)%360.0
		else:
			phi   = refrings[iref].get_attr("phi")
			theta = refrings[iref].get_attr("theta")
			psi   = (360.0+refrings[iref].get_attr("psi")-ang)%360.0
		s2x   = sxs + sxi
		s2y   = sys + syi

		#set_params_proj(data, [phi, theta, psi, s2x, s2y])
		t2 = EMAN2_cppwrap.Transform({"type":"spider","phi":phi,"theta":theta,"psi":psi})
		t2.set_trans(EMAN2_cppwrap.Vec2f(-s2x, -s2y))
		data.set_attr("xform.projection", t2)
		pass#IMPORTIMPORTIMPORT from pixel_error import max_3D_pixel_error
		ts = t2.get_sym_proj(sym)
		if(len(ts) > 1):
			# only do it if it is not c1
			pixel_error = +1.0e23
			for ut in ts:
				# we do not care which position minimizes the error
				pixel_error = min(sparx_pixel_error.max_3D_pixel_error(t1, ut, numr[-3]), pixel_error)
		else:
			pixel_error = sparx_pixel_error.max_3D_pixel_error(t1, t2, numr[-3])
		#print phi, theta, psi, s2x, s2y, peak, pixel_error
		if finfo:
			pass#IMPORTIMPORTIMPORT from utilities import get_params_proj
			phi, theta, psi, s2x, s2y = sparx_utilities.get_params_proj(data)
			finfo.write( "New parameters: %6.2f %6.2f %6.2f %6.2f %6.2f   %10.5f  %11.3e\n\n" %(phi, theta, psi, s2x, s2y, peak, pixel_error))
			finfo.flush()
		return peak, pixel_error
	else:
		return -1.0e23, 0.0


def ali_vol_func(params, data):
	pass#IMPORTIMPORTIMPORT from utilities    import model_gauss
	pass#IMPORTIMPORTIMPORT from fundamentals import rot_shift3D, cyclic_shift
	pass#IMPORTIMPORTIMPORT from morphology   import binarize
	#print  params
	#print  data[3]
	#cphi, ctheta, cpsi, cs2x, cs2y, cs2z, cscale= compose_transform3(data[3][0], data[3][1], data[3][2], data[3][3], data[3][4], data[3][5], data[3][6], params[0], params[1], params[2],params[3], params[4], params[5],1.0)
	#print  cphi, ctheta, cpsi, cs2x, cs2y, cs2z, cscale
	x = sparx_fundamentals.rot_shift3D(data[0], params[0], params[1], params[2], params[3], params[4], params[5], 1.0)

	res = -x.cmp("ccc", data[1], {"mask":data[2]})
	#print  " %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f  %10.5f" %(params[0], params[1], params[2],params[3], params[4], params[5], -res)
	return res

def align2d(image, refim, xrng=[0, 0], yrng=[0, 0], step=1, first_ring=1, last_ring=0, rstep=1, mode = "F"):
	"""  Determine shift and rotation between image and reference image
	     quadratic interpolation
	     Output: ang, sxs, sys, mirror, peak
	"""
	#from utilities import print_col
	pass#IMPORTIMPORTIMPORT from alignment import Numrinit, ringwe
	step = float(step)
	nx = refim.get_xsize()
	ny = refim.get_ysize()
	if(last_ring == 0):  last_ring = nx/2-2-int(max(max(xrng),max(yrng)))
	# center in SPIDER convention
	cnx = nx//2+1
	cny = ny//2+1
	#precalculate rings
	numr = Numrinit(first_ring, last_ring, rstep, mode)
	wr   = ringwe(numr, mode)
	#cimage=Util.Polar2Dmi(refim, cnx, cny, numr, mode, kb)
	crefim = EMAN2_cppwrap.Util.Polar2Dm(refim, cnx, cny, numr, mode)
	#crefim = Util.Polar2D(refim, numr, mode)
	#print_col(crefim)
	EMAN2_cppwrap.Util.Frngs(crefim, numr)
	EMAN2_cppwrap.Util.Applyws(crefim, numr, wr)
	return ormq(image, crefim, xrng, yrng, step, mode, numr, cnx, cny)
	
"""Multiline Comment14"""

def align2d_scf(image, refim, xrng=-1, yrng=-1, ou = -1):
	pass#IMPORTIMPORTIMPORT from fundamentals import scf, rot_shift2D, ccf, mirror
	pass#IMPORTIMPORTIMPORT from fundamentals import fft
	pass#IMPORTIMPORTIMPORT from utilities import peak_search
	pass#IMPORTIMPORTIMPORT from math import radians, sin, cos
	nx = image.get_xsize()
	ny = image.get_xsize()
	if(ou<0):  ou = min(nx//2-1,ny//2-1)
	if(yrng < 0):  yrng = xrng
	if(ou<2):
		sparx_global_def.ERROR('Radius of the object (ou) has to be given','align2d_scf',1)
	sci = sparx_fundamentals.scf(image)
	scr = sparx_fundamentals.scf(refim)
	first_ring = 1

	#alpha1, sxs, sys, mirr, peak1 = align2d_no_mirror(scf(image), scr, last_ring=ou, mode="H")
	#alpha2, sxs, sys, mirr, peak2 = align2d_no_mirror(scf(mirror(image)), scr, last_ring=ou, mode="H")
	#alpha1, sxs, sys, mirr, peak1 = align2d_no_mirror(sci, scr, first_ring = 1, last_ring=ou, mode="H")
	#alpha2, sxs, sys, mirr, peak2 = align2d_no_mirror(mirror(sci), scr,  first_ring = 1, last_ring=ou, mode="H")


	pass#IMPORTIMPORTIMPORT from alignment import Numrinit, ringwe, ornq
	# center in SPIDER convention
	cnx = nx//2+1
	cny = ny//2+1
	#precalculate rings
	numr = Numrinit(first_ring, ou, 1, "H")
	wr   = ringwe(numr, "H")
	crefim = EMAN2_cppwrap.Util.Polar2Dm(scr, cnx, cny, numr, "H")
	EMAN2_cppwrap.Util.Frngs(crefim, numr)
	EMAN2_cppwrap.Util.Applyws(crefim, numr, wr)
	alpha1, sxs, sys, mirr, peak1 = ornq(sci, crefim, [0.0], [0.0], 1.0, "H", numr, cnx, cny)
	alpha2, sxs, sys, mirr, peak2 = ornq(sparx_fundamentals.mirror(sci), crefim, [0.0], [0.0], 1.0, "H", numr, cnx, cny)


	if(peak1>peak2):
		mirr = 0
		alpha = alpha1
	else:
		mirr = 1
		alpha = -alpha2
	nrx = min( 2*(xrng+1)+1, (((nx-2)//2)*2+1) )
	nry = min( 2*(yrng+1)+1, (((ny-2)//2)*2+1) )
	frotim = sparx_fundamentals.fft( refim )
	ccf1 = EMAN2_cppwrap.Util.window(sparx_fundamentals.ccf(sparx_fundamentals.rot_shift2D(image, alpha, 0.0, 0.0, mirr), frotim),nrx,nry)
	p1 = sparx_utilities.peak_search(ccf1)
	
	ccf2 = EMAN2_cppwrap.Util.window(sparx_fundamentals.ccf(sparx_fundamentals.rot_shift2D(image, alpha+180.0, 0.0, 0.0, mirr), frotim),nrx,nry)
	p2 = sparx_utilities.peak_search(ccf2)
	#print p1
	#print p2

	peak_val1 = p1[0][0]
	peak_val2 = p2[0][0]
	
	if peak_val1 > peak_val2:
		sxs = -p1[0][4]
		sys = -p1[0][5]
		cx = int(p1[0][1])
		cy = int(p1[0][2])
		peak = peak_val1
	else:
		alpha += 180.0
		sxs = -p2[0][4]
		sys = -p2[0][5]
		peak = peak_val2
		cx = int(p2[0][1])
		cy = int(p2[0][2])
		ccf1 = ccf2
	pass#IMPORTIMPORTIMPORT from utilities import model_blank
	#print cx,cy
	z = sparx_utilities.model_blank(3,3)
	for i in range(3):
		for j in range(3):
			z[i,j] = ccf1[i+cx-1,j+cy-1]
	#print  ccf1[cx,cy],z[1,1]
	XSH, YSH, PEAKV = parabl(z)
	#print sxs, sys, XSH, YSH, PEAKV, peak
	if(mirr == 1):  	sx = -sxs+XSH
	else:               sx =  sxs-XSH
	return alpha, sx, sys-YSH, mirr, PEAKV



def multalign2d_scf(image, refrings, frotim, numr, xrng=-1, yrng=-1, ou = -1):
	pass#IMPORTIMPORTIMPORT from fundamentals import scf, rot_shift2D, ccf, mirror
	pass#IMPORTIMPORTIMPORT from utilities import peak_search, model_blank
	pass#IMPORTIMPORTIMPORT from math import radians, sin, cos
	pass#IMPORTIMPORTIMPORT from alignment import ang_n

	nx = image.get_xsize()
	ny = image.get_xsize()
	if(ou<0):  ou = min(nx//2-1,ny//2-1)
	if(yrng < 0):  yrng = xrng
	if(ou<2):
		sparx_global_def.ERROR('Radius of the object (ou) has to be given','align2d_scf',1)
	sci = sparx_fundamentals.scf(image)
	first_ring = 1
	# center in SPIDER convention
	cnx = nx//2+1
	cny = ny//2+1

	cimage = EMAN2_cppwrap.Util.Polar2Dm(sci, cnx, cny, numr, "H")
	EMAN2_cppwrap.Util.Frngs(cimage, numr)
	mimage = EMAN2_cppwrap.Util.Polar2Dm(sparx_fundamentals.mirror(sci), cnx, cny, numr, "H")
	EMAN2_cppwrap.Util.Frngs(mimage, numr)

	nrx = min( 2*(xrng+1)+1, (((nx-2)//2)*2+1) )
	nry = min( 2*(yrng+1)+1, (((ny-2)//2)*2+1) )

	totpeak = -1.0e23

	for iki in range(len(refrings)):
		#print  "TEMPLATE  ",iki
		#  Find angle
		retvals = EMAN2_cppwrap.Util.Crosrng_e(refrings[iki], cimage, numr, 0, 0.0)
		alpha1  = ang_n(retvals["tot"], "H", numr[-1])
		peak1 	= retvals["qn"]
		retvals = EMAN2_cppwrap.Util.Crosrng_e(refrings[iki], mimage, numr, 0, 0.0)
		alpha2  = ang_n(retvals["tot"], "H", numr[-1])
		peak2 	= retvals["qn"]
		#print  alpha1, peak1
		#print  alpha2, peak2

		if(peak1>peak2):
			mirr = 0
			alpha = alpha1
		else:
			mirr = 1
			alpha = -alpha2

		ccf1 = EMAN2_cppwrap.Util.window(sparx_fundamentals.ccf(sparx_fundamentals.rot_shift2D(image, alpha, 0.0, 0.0, mirr), frotim[iki]), nrx, nry)
		p1 = sparx_utilities.peak_search(ccf1)
	
		ccf2 = EMAN2_cppwrap.Util.window(sparx_fundamentals.ccf(sparx_fundamentals.rot_shift2D(image, alpha+180.0, 0.0, 0.0, mirr), frotim[iki]), nrx, nry)
		p2 = sparx_utilities.peak_search(ccf2)
		#print p1
		#print p2

		peak_val1 = p1[0][0]
		peak_val2 = p2[0][0]
	
		if peak_val1 > peak_val2:
			sxs = -p1[0][4]
			sys = -p1[0][5]
			cx = int(p1[0][1])
			cy = int(p1[0][2])
			peak = peak_val1
		else:
			alpha += 180.0
			sxs = -p2[0][4]
			sys = -p2[0][5]
			peak = peak_val2
			cx = int(p2[0][1])
			cy = int(p2[0][2])
			ccf1 = ccf2
		#print cx,cy
		z = sparx_utilities.model_blank(3,3)
		for i in range(3):
			for j in range(3):
				z[i,j] = ccf1[i+cx-1,j+cy-1]
		#print  ccf1[cx,cy],z[1,1]
		XSH, YSH, PEAKV = parabl(z)
		#print  PEAKV
		if(PEAKV > totpeak):
			totpeak = PEAKV
			iref = iki
			if(mirr == 1):  	sx = -sxs+XSH
			else:               sx =  sxs-XSH
			sy = sys-YSH
			talpha = alpha
			tmirr = mirr
			#print "BETTER",sx,sy,iref,talpha,tmirr,totpeak
			#return alpha, sx, sys-YSH, mirr, PEAKV
	return sx,sy,iref,talpha,tmirr,totpeak

def parabl(Z):
	#  parabolic fit to a peak, C indexing
	C1 = (26.*Z[0,0] - Z[0,1] + 2*Z[0,2] - Z[1,0] - 19.*Z[1,1] - 7.*Z[1,2] + 2.*Z[2,0] - 7.*Z[2,1] + 14.*Z[2,2])/9.

	C2 = (8.* Z[0,0] - 8.*Z[0,1] + 5.*Z[1,0] - 8.*Z[1,1] + 3.*Z[1,2] +2.*Z[2,0] - 8.*Z[2,1] + 6.*Z[2,2])/(-6.)

	C3 = (Z[0,0] - 2.*Z[0,1] + Z[0,2] + Z[1,0] -2.*Z[1,1] + Z[1,2] + Z[2,0] - 2.*Z[2,1] + Z[2,2])/6.

	C4 = (8.*Z[0,0] + 5.*Z[0,1] + 2.*Z[0,2] -8.*Z[1,0] -8.*Z[1,1] - 8.*Z[1,2] + 3.*Z[2,1] + 6.*Z[2,2])/(-6.)

	C5 = (Z[0,0] - Z[0,2] - Z[2,0] + Z[2,2])/4.

	C6 = (Z[0,0] + Z[0,1] + Z[0,2] - 2.*Z[1,0] - 2.*Z[1,1] -2.*Z[1,2] + Z[2,0] + Z[2,1] + Z[2,2])/6.

	DENOM = 4. * C3 * C6 - C5 * C5
	if(DENOM == 0.):
		return 0.0, 0.0, 0.0

	YSH   = (C4*C5 - 2.*C2*C6) / DENOM - 2.
	XSH   = (C2*C5 - 2.*C4*C3) / DENOM - 2.

	PEAKV = 4.*C1*C3*C6 - C1*C5*C5 - C2*C2*C6 + C2*C4*C5 - C4*C4*C3
	PEAKV = PEAKV / DENOM
	#print  "  in PARABL  ",XSH,YSH,Z[1,1],PEAKV

	XSH = min(max( XSH, -1.0), 1.0)
	YSH = min(max( YSH, -1.0), 1.0)

	return XSH, YSH, PEAKV
"""Multiline Comment15"""

def shc(data, refrings, list_of_reference_angles, numr, xrng, yrng, step, an = -1.0, sym = "c1", finfo=None):
	pass#IMPORTIMPORTIMPORT from alignment import search_range
	pass#IMPORTIMPORTIMPORT from math         import cos, sin, degrees, radians
	pass#IMPORTIMPORTIMPORT from EMAN2 import Vec2f

	number_of_checked_refs = 0

	mode = "F"
	nx   = data.get_xsize()
	ny   = data.get_ysize()
	#  center is in SPIDER convention
	cnx  = nx//2 + 1
	cny  = ny//2 + 1

	if( an>= 0.0):  ant = numpy.cos(numpy.radians(an))
	else:           ant = -1.0
	#phi, theta, psi, sxo, syo = get_params_proj(data)
	t1 = data.get_attr("xform.projection")
	dp = t1.get_params("spider")
	ou = numr[-3]
	sxi = round(-dp["tx"],2)
	syi = round(-dp["ty"],2)
	txrng = search_range(nx, ou, sxi, xrng)
	tyrng = search_range(ny, ou, syi, yrng)

	if finfo:
		finfo.write("Old parameters: %9.4f %9.4f %9.4f %9.4f %9.4f\n"%(dp["phi"], dp["theta"], dp["psi"], -dp["tx"], -dp["ty"]))
		finfo.flush()
		pass#IMPORTIMPORTIMPORT from utilities import get_params_proj
		z1,z2,z3,z4,z5 = sparx_utilities.get_params_proj(data, "xform.anchor")
		finfo.write("Anc parameters: %9.4f %9.4f %9.4f %9.4f %9.4f\n"%(z1,z2,z3,-z4,-z5))
		finfo.flush()

	previousmax = data.get_attr("previousmax")
	[ang, sxs, sys, mirror, iref, peak, checked_refs] = EMAN2_cppwrap.Util.shc(data, refrings, list_of_reference_angles, txrng, tyrng, step, ant, mode, numr, cnx-sxi, cny-syi, sym)
	iref=int(iref)
	number_of_checked_refs += int(checked_refs)
	if peak <= previousmax:
		return -1.0e23, 0.0, number_of_checked_refs, -1
		"""Multiline Comment50"""
	else:
		# The ormqip returns parameters such that the transformation is applied first, the mirror operation second.
		# What that means is that one has to change the the Eulerian angles so they point into mirrored direction: phi+180, 180-theta, 180-psi
		if  mirror:
			phi   = (refrings[iref].get_attr("phi")+540.0)%360.0
			theta = 180.0-refrings[iref].get_attr("theta")
			psi   = (540.0-refrings[iref].get_attr("psi")-ang)%360.0
		else:
			phi   = refrings[iref].get_attr("phi")
			theta = refrings[iref].get_attr("theta")
			psi   = (360.0+refrings[iref].get_attr("psi")-ang)%360.0
		s2x   = sxs + sxi
		s2y   = sys + syi

		#set_params_proj(data, [phi, theta, psi, s2x, s2y])
		t2 = EMAN2_cppwrap.Transform({"type":"spider","phi":phi,"theta":theta,"psi":psi})
		t2.set_trans(EMAN2_cppwrap.Vec2f(-s2x, -s2y))
		data.set_attr("xform.projection", t2)
		data.set_attr("previousmax", peak)
		#  Find the pixel error that is minimum over symmetry transformations
		pass#IMPORTIMPORTIMPORT from pixel_error import max_3D_pixel_error
		if(sym == "nomirror" or sym == "c1"):
			pixel_error = sparx_pixel_error.max_3D_pixel_error(t1, t2, numr[-3])
		else:		
			ts = t2.get_sym_proj(sym)
			# only do it if it is not c1
			pixel_error = +1.0e23
			for ut in ts:
				# we do not care which position minimizes the error
				pixel_error = min(sparx_pixel_error.max_3D_pixel_error(t1, ut, numr[-3]), pixel_error)
		if finfo:
			finfo.write( "New parameters: %9.4f %9.4f %9.4f %9.4f %9.4f %10.5f  %11.3e\n\n" %(phi, theta, psi, s2x, s2y, peak, pixel_error))
			finfo.flush()
		return peak, pixel_error, number_of_checked_refs, iref


# parameters: list of (all) projections | reference volume is optional, if provided might be shrank| ...
#  This functions centers projections using an self-correlation-based exhaustive search
#  It only returns shifts
#  Data is assumed to be shrunk and CTF-applied
#  The input volume is assumed to be shrunk but not filtered, if not provided, it will be reconstructed and shrunk
#  We apply ali3d_options.fl
def search_range(n, radius, shift, range, location = ""):
	"""
		Find permissible ranges for translational searches by resampling into polar coordinates
		n - image size; radius - particle radius, the circle has to fit into the square image;
		shift - current particle shift; range - desired maximum range search
		Output: a list of two elements:
		  left range (positive)
		  right range
		NOTE - ranges are with respect to the point n//2+1-shift within image (in 3D)
	"""
	cn = n//2 +1
	ql = cn+shift-radius -2   # lower end is positive
	qe = n - cn-shift-radius    # upper end
	if( ql < 0 or qe < 0 ):
		sparx_global_def.ERROR("Shift of particle too large, results may be incorrect:  %4d   %3d   %f  %f  %f  %f  %f"%(n, cn, radius, shift, range, ql, qe),"search_range  "+location,0)
		ql = max(ql,0)
		qe = max(qe,0)
	# ???for mysterious reasons it has to be this way as C code changes the order of searches.
	return  [ min( qe, range), min(ql, range) ]


def generate_list_of_reference_angles_for_search(input_angles, sym):
	"""
	  Generate full set of reference angles, including mirror and symmetry related
	  from a unique subrange generated by even_angles and stored in refrings.
	  Input - input_angles [[angles],[angles]]
	  Output - [[angles], [angles]] (no shifts)
			Blocks - [[basic][mirrored basic]] [[basic sym1][mirrored basic sym1]] ...
	"""
	pass#IMPORTIMPORTIMPORT from EMAN2 import Transform
	t2   = EMAN2_cppwrap.Transform()
	nsym = t2.get_nsym(sym)

	original_number_of_angles = len(input_angles)
	# original_number_of_angles is the same as the number of refrings
	
	list_of_reference_angles = [None]*original_number_of_angles
	for i in range(original_number_of_angles): 
		list_of_reference_angles[i] = [input_angles[i][0],input_angles[i][1], 0]

	#  add mirror related
	list_of_reference_angles += [[0.0,0.0,0.0] for i in range(original_number_of_angles)]
	for i in range(original_number_of_angles):
		list_of_reference_angles[i+original_number_of_angles][0] = (list_of_reference_angles[i][0]+180.0)%360.0
		list_of_reference_angles[i+original_number_of_angles][1] = 180.0-list_of_reference_angles[i][1]
		list_of_reference_angles[i+original_number_of_angles][2] =  list_of_reference_angles[i][2]

	#  add symmetry related
	if(nsym>1):	
		number_of_angles_original_and_mirror = len(list_of_reference_angles)
		for l in range(1,nsym):
			list_of_reference_angles += [[0.0,0.0,0.0] for i in range(number_of_angles_original_and_mirror)]

		for i in range(number_of_angles_original_and_mirror):
			t2 = EMAN2_cppwrap.Transform({"type":"spider","phi":list_of_reference_angles[i][0],"theta":list_of_reference_angles[i][1]})
			ts = t2.get_sym_proj(sym)
			for ll in range(1,nsym,1):
				d = ts[ll].get_params("spider")
				list_of_reference_angles[i+ll*number_of_angles_original_and_mirror][0] = round(d["phi"],5)
				list_of_reference_angles[i+ll*number_of_angles_original_and_mirror][1] = round(d["theta"],5)
				list_of_reference_angles[i+ll*number_of_angles_original_and_mirror][2] = round(d["psi"],5)  #  Not needed?

	return list_of_reference_angles


