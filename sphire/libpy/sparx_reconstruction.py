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
import sparx_filter
import sparx_fundamentals
import sparx_global_def
import sparx_morphology
import mpi
import numpy
import numpy.random
import os
import random
import sparx_statistics
import string
import sys
import time
import sparx_utilities
pass#IMPORTIMPORTIMPORT import EMAN2
pass#IMPORTIMPORTIMPORT import EMAN2_cppwrap
pass#IMPORTIMPORTIMPORT import datetime
pass#IMPORTIMPORTIMPORT import filter
pass#IMPORTIMPORTIMPORT import fundamentals
pass#IMPORTIMPORTIMPORT import global_def
pass#IMPORTIMPORTIMPORT import math
pass#IMPORTIMPORTIMPORT import morphology
pass#IMPORTIMPORTIMPORT import mpi
pass#IMPORTIMPORTIMPORT import numpy
pass#IMPORTIMPORTIMPORT import numpy.random
pass#IMPORTIMPORTIMPORT import os
pass#IMPORTIMPORTIMPORT import projection
pass#IMPORTIMPORTIMPORT import random
pass#IMPORTIMPORTIMPORT import reconstruction
pass#IMPORTIMPORTIMPORT import statistics
pass#IMPORTIMPORTIMPORT import string
pass#IMPORTIMPORTIMPORT import subprocess
pass#IMPORTIMPORTIMPORT import sys
pass#IMPORTIMPORTIMPORT import time
pass#IMPORTIMPORTIMPORT import types
pass#IMPORTIMPORTIMPORT import utilities
from builtins import range
from builtins import object
pass#IMPORTIMPORTIMPORT from global_def import *

def insert_slices(reconstructor, proj):
	xforms = [ proj.get_attr("xform.projection") ]
	weights = [ proj.get_attr_default("weight", 1.0) ]
	ixform = 0
	while True:
		ixform += 1
		xform_proj = proj.get_attr_default("xform.projection" + str(ixform), None)
		if xform_proj == None:
			break
		# putting params in a list does not seem to be necessary, one could call reconstructor as one goes.
		xforms.append(xform_proj)
		#weights.append(proj.get_attr_default("weight" + str(ixform), 1.0))
		weights.append(1.0)
	for i in range(len(xforms)):
		reconstructor.insert_slice( proj, xforms[i], weights[i] )

def insert_slices_pdf(reconstructor, proj):
	xforms =   proj.get_attr("xform.projection") 
	weights =  proj.get_attr_default("weight", 1.0) 
	reconstructor.insert_slice( proj, xforms, weights )
	ixform = 0
	while True:
		ixform += 1
		xform_proj = proj.get_attr_default("xform.projection" + str(ixform), None)
		if xform_proj == None:
			return 
		weights = proj.get_attr_default("weight" + str(ixform), 1.0)
		reconstructor.insert_slice( proj, xforms, weights)

def recons3d_4nn_MPI(myid, prjlist, symmetry="c1", finfo=None, snr = 1.0, npad=2, xysize=-1, zsize=-1, mpi_comm=None):
	pass#IMPORTIMPORTIMPORT from utilities  import reduce_EMData_to_root, pad
	pass#IMPORTIMPORTIMPORT from EMAN2      import Reconstructors
	pass#IMPORTIMPORTIMPORT from utilities  import iterImagesList
	pass#IMPORTIMPORTIMPORT from mpi        import MPI_COMM_WORLD
	pass#IMPORTIMPORTIMPORT import types

	if mpi_comm == None:
		mpi_comm = mpi.MPI_COMM_WORLD

	if type(prjlist) == list:
		prjlist = sparx_utilities.iterImagesList(prjlist)

	if not prjlist.goToNext():
		sparx_global_def.ERROR("empty input list","recons3d_4nn_MPI",1)

	imgsize = prjlist.image().get_xsize()
	if prjlist.image().get_ysize() != imgsize:
		imgsize = max(imgsize, prjlist.image().get_ysize())
		dopad = True
	else:
		dopad = False
	prjlist.goToPrev()

	fftvol = EMAN2_cppwrap.EMData()		
	weight = EMAN2_cppwrap.EMData()
	if (xysize == -1 and zsize == -1 ):
		params = {"size":imgsize, "npad":npad, "symmetry":symmetry, "fftvol":fftvol, "weight":weight, "snr":snr}
		r = EMAN2_cppwrap.Reconstructors.get( "nn4", params )
	else:
		if ( xysize != -1 and zsize != -1):
			rx = float(xysize)/imgsize
			ry = float(xysize)/imgsize
			rz = float(zsize)/imgsize
		elif( xysize != -1):
			rx = float(xysize)/imgsize
			ry = float(xysize)/imgsize
			rz = 1.0
		else:
			rx = 1.0
			ry = 1.0
			rz = float(zsize)/imgsize
		params = {"sizeprojection":imgsize, "npad":npad, "symmetry":symmetry, "fftvol":fftvol,"weight":weight,"xratio":rx,"yratio":ry,"zratio":rz}
		r = EMAN2_cppwrap.Reconstructors.get( "nn4_rect", params )
	r.setup()

	if not (finfo is None): nimg = 0
	while prjlist.goToNext():
		prj = prjlist.image()
		if dopad:
			prj = sparx_utilities.pad(prj, imgsize,imgsize, 1, "circumference")
		insert_slices(r, prj)
		if( not (finfo is None) ):
			nimg += 1
			finfo.write("Image %4d inserted.\n" %(nimg) )
			finfo.flush()

	if not (finfo is None): 
		finfo.write( "Begin reducing ...\n" )
		finfo.flush()

	sparx_utilities.reduce_EMData_to_root(fftvol, myid, comm=mpi_comm)
	sparx_utilities.reduce_EMData_to_root(weight, myid, comm=mpi_comm)

	if myid == 0:  dummy = r.finish(True)
	else:
		pass#IMPORTIMPORTIMPORT from utilities import model_blank
		if ( xysize == -1 and zsize == -1 ):
			fftvol = sparx_utilities.model_blank(imgsize, imgsize, imgsize)
		else:
			if zsize == -1:
				fftvol = sparx_utilities.model_blank(xysize, xysize, imgsize)
			elif xysize == -1:
				fftvol = sparx_utilities.model_blank(imgsize, imgsize, zsize)
			else:
				fftvol = sparx_utilities.model_blank(xysize, xysize, zsize)
	return fftvol

"""Multiline Comment0"""
"""Multiline Comment1"""

def recons3d_trl_struct_MPI(myid, main_node, prjlist, paramstructure, refang, rshifts_shrank, delta, upweighted = True, mpi_comm=None, CTF = True, target_size=-1, avgnorm = 1.0, norm_per_particle = None):
	"""
		recons3d_4nn_ctf - calculate CTF-corrected 3-D reconstruction from a set of projections using three Eulerian angles, two shifts, and CTF settings for each projeciton image
		Input
			list_of_prjlist: list of lists of projections to be included in the reconstruction
	"""
	pass#IMPORTIMPORTIMPORT from utilities  import reduce_EMData_to_root, random_string, get_im, findall
	pass#IMPORTIMPORTIMPORT from EMAN2      import Reconstructors
	pass#IMPORTIMPORTIMPORT from utilities  import model_blank
	pass#IMPORTIMPORTIMPORT from filter	import filt_table
	pass#IMPORTIMPORTIMPORT from fundamentals import fshift
	pass#IMPORTIMPORTIMPORT from mpi        import MPI_COMM_WORLD, mpi_barrier
	pass#IMPORTIMPORTIMPORT import types
	pass#IMPORTIMPORTIMPORT import datetime
	
	if mpi_comm == None: mpi_comm = mpi.MPI_COMM_WORLD

	refvol = sparx_utilities.model_blank(target_size)
	refvol.set_attr("fudge", 1.0)

	if CTF: do_ctf = 1
	else:   do_ctf = 0

	fftvol = EMAN2_cppwrap.EMData()
	weight = EMAN2_cppwrap.EMData()

	pass#IMPORTIMPORTIMPORT from utilities import info
	params = {"size":target_size, "npad":2, "snr":1.0, "sign":1, "symmetry":"c1", "refvol":refvol, "fftvol":fftvol, "weight":weight, "do_ctf": do_ctf}
	r = EMAN2_cppwrap.Reconstructors.get( "nn4_ctfw", params )
	r.setup()
	
	if norm_per_particle == None: norm_per_particle = len(prjlist)*[1.0]

	nnx = prjlist[0].get_xsize()
	nny = prjlist[0].get_ysize()
	nshifts = len(rshifts_shrank)
	for im in range(len(prjlist)):
		#  parse projection structure, generate three lists:
		#  [ipsi+iang], [ishift], [probability]
		#  Number of orientations for a given image
		numbor = len(paramstructure[im][2])
		ipsiandiang = [ paramstructure[im][2][i][0]/1000  for i in range(numbor) ]
		allshifts   = [ paramstructure[im][2][i][0]%1000  for i in range(numbor) ]
		probs       = [ paramstructure[im][2][i][1] for i in range(numbor) ]
		#  Find unique projection directions
		tdir = list(set(ipsiandiang))
		bckgn = prjlist[im].get_attr("bckgnoise")
		ct = prjlist[im].get_attr("ctf")
		#  For each unique projection direction:
		data = [None]*nshifts
		for ii in range(len(tdir)):
			#  Find the number of times given projection direction appears on the list, it is the number of different shifts associated with it.
			lshifts = sparx_utilities.findall(tdir[ii], ipsiandiang)
			toprab  = 0.0
			for ki in range(len(lshifts)):  toprab += probs[lshifts[ki]]
			recdata = EMAN2_cppwrap.EMData(nny,nny,1,False)
			recdata.set_attr("is_complex",0)
			for ki in range(len(lshifts)):
				lpt = allshifts[lshifts[ki]]
				if( data[lpt] == None ):
					data[lpt] = sparx_fundamentals.fshift(prjlist[im], rshifts_shrank[lpt][0], rshifts_shrank[lpt][1])
					data[lpt].set_attr("is_complex",0)
				EMAN2_cppwrap.Util.add_img(recdata, EMAN2_cppwrap.Util.mult_scalar(data[lpt], probs[lshifts[ki]]/toprab))
			recdata.set_attr_dict({"padffted":1, "is_fftpad":1,"is_fftodd":0, "is_complex_ri":1, "is_complex":1})
			if not upweighted:  recdata = sparx_filter.filt_table(recdata, bckgn )
			recdata.set_attr_dict( {"bckgnoise":bckgn, "ctf":ct} )
			ipsi = tdir[ii]%100000
			iang = tdir[ii]/100000
			r.insert_slice( recdata, EMAN2_cppwrap.Transform({"type":"spider","phi":refang[iang][0],"theta":refang[iang][1],"psi":refang[iang][2]+ipsi*delta}), toprab*avgnorm/norm_per_particle[im])
	#  clean stuff
	del bckgn, recdata, tdir, ipsiandiang, allshifts, probs


	sparx_utilities.reduce_EMData_to_root(fftvol, myid, main_node, comm=mpi_comm)
	sparx_utilities.reduce_EMData_to_root(weight, myid, main_node, comm=mpi_comm)

	if myid == main_node:
		dummy = r.finish(True)
	mpi.mpi_barrier(mpi_comm)

	if myid == main_node: return fftvol, weight, refvol
	else: return None, None, None


def recons3d_4nn_ctf(stack_name, list_proj = [], snr = 1.0, sign=1, symmetry="c1", verbose=0, npad=2, xysize = -1, zsize = -1 ):
	"""Perform a 3-D reconstruction using Pawel's FFT Back Projection algoritm.

	   Input:
	    stack_name - name of the stack file on a disk,
	                 each image has to have the following attributes set:
			 psi, theta, phi, sx, sy, defocus, 
	    list_proj - list of images from stack_name to be included in the reconstruction
	    symmetry	 -- Point group of the target molecule (defaults to "C1")

	   Return:  3d reconstructed volume image

	   Usage:
	     
	     anglelist = getAngles("myangles.txt") # not yet written
	     vol = do_reconstruction(filepattern, start, end, anglelist, symmetry)
	"""
	pass#IMPORTIMPORTIMPORT import types
	pass#IMPORTIMPORTIMPORT from EMAN2     import Reconstructors
	pass#IMPORTIMPORTIMPORT from utilities import pad

	# read first image to determine the size to use
	if list_proj == []:	
		if type(stack_name) == bytes: nima = EMAN2_cppwrap.EMUtil.get_image_count(stack_name)
		else : nima = len(stack_name)
		list_proj = list(range(nima)) 
	# read first image to determine the size to use
	if type(stack_name) == bytes:
		proj = EMAN2_cppwrap.EMData()
		proj.read_image(stack_name, list_proj[0])
	else:    proj = stack_name[list_proj[0]].copy()

	# convert angles to transform (rotation) objects
	# horatio active_refactoring Jy51i1EwmLD4tWZ9_00000_1
	# active = proj.get_attr_default('active', 1)
	size   = proj.get_xsize()
	if proj.get_ysize() != size:
		size = max(size, proj.get_ysize())
		dopad = True
	else:
		dopad = False

	# reconstructor
	fftvol = EMAN2_cppwrap.EMData()
	weight = EMAN2_cppwrap.EMData()
	params = {"npad":npad, "symmetry":symmetry, "snr":snr, "sign":sign, "fftvol":fftvol, "weight":weight}
	if ( xysize == -1 and zsize == -1 ):
		params["size"] = size
		r = EMAN2_cppwrap.Reconstructors.get("nn4_ctf", params)
	else:
		if ( xysize != -1 and zsize != -1):
			rx = float(xysize)/size
			ry = float(xysize)/size
			rz = float(zsize)/size
		elif( xysize != -1):
			rx = float(xysize)/size
			ry = float(xysize)/size
			rz = 1.0
		else:
			rx = 1.0
			ry = 1.0
			rz = float(zsize)/size

		params["sizeprojection"] = size
		params["xratio"] = rx
		params["yratio"] = ry
		params["zratio"] = rz
		r = EMAN2_cppwrap.Reconstructors.get("nn4_ctf_rect", params)
	r.setup()

	if type(stack_name) == bytes:
		for i in range(len(list_proj)):
			proj.read_image(stack_name, list_proj[i])
			if dopad: 
				proj = sparx_utilities.pad(proj, size, size, 1, "circumference")
			insert_slices(r, proj)
	else:
		for i in range(len(list_proj)):
			insert_slices(r, stack_name[list_proj[i]])
	dummy = r.finish(True)
	return fftvol


def recons3d_4nn_ctf_MPI(myid, prjlist, snr = 1.0, sign=1, symmetry="c1", finfo=None, npad=2, xysize=-1, zsize=-1, mpi_comm=None, smearstep = 0.0):
	"""
		recons3d_4nn_ctf - calculate CTF-corrected 3-D reconstruction from a set of projections using three Eulerian angles, two shifts, and CTF settings for each projeciton image
		Input
			stack: name of the stack file containing projection data, projections have to be squares
			list_proj: list of projections to be included in the reconstruction or image iterator
			snr: Signal-to-Noise Ratio of the data 
			sign: sign of the CTF 
			symmetry: point-group symmetry to be enforced, each projection will enter the reconstruction in all symmetry-related directions.
	"""
	pass#IMPORTIMPORTIMPORT from utilities  import reduce_EMData_to_root, pad
	pass#IMPORTIMPORTIMPORT from EMAN2      import Reconstructors
	pass#IMPORTIMPORTIMPORT from utilities  import iterImagesList, set_params_proj
	pass#IMPORTIMPORTIMPORT from mpi        import MPI_COMM_WORLD
	pass#IMPORTIMPORTIMPORT import types

	if mpi_comm == None:
		mpi_comm = mpi.MPI_COMM_WORLD

	if type(prjlist) == list:
		prjlist = sparx_utilities.iterImagesList(prjlist)
	if not prjlist.goToNext():
		sparx_global_def.ERROR("empty input list","recons3d_4nn_ctf_MPI",1)
	imgsize = prjlist.image().get_xsize()
	if prjlist.image().get_ysize() != imgsize:
		imgsize = max(imgsize, prjlist.image().get_ysize())
		dopad = True
	else:
		dopad = False
	prjlist.goToPrev()

	fftvol = EMAN2_cppwrap.EMData()

	if( smearstep > 0.0 ):
		#if myid == 0:  print "  Setting smear in prepare_recons_ctf"
		ns = 1
		smear = []
		for j in range(-ns,ns+1):
			if( j != 0):
				for i in range(-ns,ns+1):
					for k in range(-ns,ns+1):
						smear += [i*smearstep,j*smearstep,k*smearstep,1.0]
		# Deal with theta = 0.0 cases
		prj = []
		for i in range(-ns,ns+1):
			for k in range(-ns,ns+1):
				prj.append(i+k)
		for i in range(-2*ns,2*ns+1,1):
			smear += [i*smearstep,0.0,0.0,float(prj.count(i))]
		#if myid == 0:  print "  Smear  ",smear
		fftvol.set_attr("smear", smear)

	weight = EMAN2_cppwrap.EMData()
	if (xysize == -1 and zsize == -1 ):
		params = {"size":imgsize, "npad":npad, "snr":snr, "sign":sign, "symmetry":symmetry, "fftvol":fftvol, "weight":weight}
		r = EMAN2_cppwrap.Reconstructors.get( "nn4_ctf", params )
	else:
		if ( xysize != -1 and zsize != -1):
			rx = float(xysize)/imgsize
			ry = float(xysize)/imgsize
			rz = float(zsize)/imgsize
		elif( xysize != -1):
			rx = float(xysize)/imgsize
			ry = float(xysize)/imgsize
			rz = 1.0
		else:
			rx = 1.0
			ry = 1.0
			rz = float(zsize)/imgsize
		#  There is an error here with sizeprojection  PAP 10/22/2014
		params = {"size":sizeprojection, "npad":npad, "snr":snr, "sign":sign, "symmetry":symmetry, "fftvol":fftvol, "weight":weight,"xratio":rx,"yratio":ry,"zratio":rz}
		r = EMAN2_cppwrap.Reconstructors.get( "nn4_ctf_rect", params )
	r.setup()

	#if not (finfo is None):
	nimg = 0
	while prjlist.goToNext():
		prj = prjlist.image()
		if dopad:
			prj = sparx_utilities.pad(prj, imgsize, imgsize, 1, "circumference")
		#if params:
		insert_slices(r, prj)
		if not (finfo is None):
			nimg += 1
			finfo.write(" %4d inserted\n" %(nimg) )
			finfo.flush()
	del sparx_utilities.pad
	if not (finfo is None): 
		finfo.write( "begin reduce\n" )
		finfo.flush()

	sparx_utilities.reduce_EMData_to_root(fftvol, myid, comm=mpi_comm)
	sparx_utilities.reduce_EMData_to_root(weight, myid, comm=mpi_comm)

	if not (finfo is None): 
		finfo.write( "after reduce\n" )
		finfo.flush()

	if myid == 0 :
		dummy = r.finish(True)
	else:
		pass#IMPORTIMPORTIMPORT from utilities import model_blank
		if ( xysize == -1 and zsize == -1 ):
			fftvol = sparx_utilities.model_blank(imgsize, imgsize, imgsize)
		else:
			if zsize == -1:
				fftvol = sparx_utilities.model_blank(xysize, xysize, imgsize)
			elif xysize == -1:
				fftvol = sparx_utilities.model_blank(imgsize, imgsize, zsize)
			else:
				fftvol = sparx_utilities.model_blank(xysize, xysize, zsize)
	return fftvol


def recons3d_nn_SSNR_MPI(myid, prjlist, mask2D, ring_width=1, npad =1, sign=1, symmetry="c1", CTF = False, random_angles = 0, mpi_comm = None):
	pass#IMPORTIMPORTIMPORT from utilities import reduce_EMData_to_root
	pass#IMPORTIMPORTIMPORT from EMAN2 import Reconstructors
	pass#IMPORTIMPORTIMPORT from mpi import MPI_COMM_WORLD

	if mpi_comm == None:
		mpi_comm = mpi.MPI_COMM_WORLD

	if( len(prjlist) == 0 ):    sparx_global_def.ERROR("empty input list","recons3d_nn_SSNR_MPI",1)
	imgsize = prjlist[0].get_xsize()
	if prjlist[0].get_ysize() != imgsize:  sparx_global_def.ERROR("input data has to be square","recons3d_nn_SSNR_MPI",1)
	fftvol   = EMAN2_cppwrap.EMData()
	weight   = EMAN2_cppwrap.EMData()
	weight2  = EMAN2_cppwrap.EMData()
	SSNR     = EMAN2_cppwrap.EMData()
	vol_ssnr = EMAN2_cppwrap.EMData()
	params = {"size":imgsize, "npad":npad, "symmetry":symmetry, "SSNR":SSNR, "fftvol":fftvol, "weight":weight, "weight2":weight2, "vol_ssnr":vol_ssnr, "w":ring_width }
	if CTF:
		weight3  = EMAN2_cppwrap.EMData()
		params["sign"] = sign
		params["weight3"] = weight3
		r = EMAN2_cppwrap.Reconstructors.get("nnSSNR_ctf", params)
	else:
		r = EMAN2_cppwrap.Reconstructors.get("nnSSNR", params)
	r.setup()

	if prjlist[0].get_xsize() != imgsize or prjlist[0].get_ysize() != imgsize: sparx_global_def.ERROR("inconsistent image size","recons3d_nn_SSNR_MPI",1)
	for prj in prjlist:
		# horatio active_refactoring Jy51i1EwmLD4tWZ9_00000_1
		# active = prj.get_attr_default('active', 1)
		# if active == 1:
		pass#IMPORTIMPORTIMPORT import numpy.random
		if random_angles  == 2:
			pass#IMPORTIMPORTIMPORT from  random import  random
			phi	 = 360.0*numpy.random.random()
			theta    = 180.0*numpy.random.random()
			psi	 = 360.0*numpy.random.random()
			xform_proj = EMAN2_cppwrap.Transform( {"type":"spider", "phi":phi, "theta":theta, "psi":psi} )
		elif random_angles  == 3:
			pass#IMPORTIMPORTIMPORT from  random import  random
			phi    = 360.0*numpy.random.random()
			theta  = 180.0*numpy.random.random()
			psi    = 360.0*numpy.random.random()
			tx     = 6.0*(numpy.random.random() - 0.5)
			ty     = 6.0*(numpy.random.random() - 0.5)
			xform_proj = EMAN2_cppwrap.Transform( {"type":"spider", "phi":phi, "theta":theta, "psi":psi, "tx":tx, "ty":ty} )
		elif random_angles  == 1:
			pass#IMPORTIMPORTIMPORT from  random import  random
			old_xform_proj = prj.get_attr( "xform.projection" )
			dict = old_xform_proj.get_rotation( "spider" )
			dict["psi"] = 360.0*numpy.random.random()
			xform_proj = EMAN2_cppwrap.Transform( dict )
		else:
			xform_proj = prj.get_attr( "xform.projection" )
		if mask2D:
			stats = EMAN2_cppwrap.Util.infomask(prj, mask2D, True)
			prj -= stats[0]
			prj *= mask2D
		r.insert_slice(prj, xform_proj )
		# horatio active_refactoring Jy51i1EwmLD4tWZ9_00000_1 END

	#from utilities import info
	sparx_utilities.reduce_EMData_to_root(weight,  myid, 0, comm=mpi_comm)
	sparx_utilities.reduce_EMData_to_root(fftvol,  myid, 0, comm=mpi_comm)
	sparx_utilities.reduce_EMData_to_root(weight2, myid, 0, comm=mpi_comm)
	if CTF:
		sparx_utilities.reduce_EMData_to_root(weight3, myid, 0, comm=mpi_comm)
	if myid == 0 :
		dummy = r.finish(True)		
		outlist = [[] for i in range(6)]
		nn = SSNR.get_xsize()
		for i in range(1,nn): outlist[0].append((float(i)-0.5)/(float(nn-1)*2))
		for i in range(1,nn):
			if SSNR(i,1,0) > 0.0:
				outlist[1].append(max(0.0,(SSNR(i,0,0)/SSNR(i,1,0)-1.)))     # SSNR
			else:
				outlist[1].append(0.0)
		for i in range(1,nn): 
			if SSNR(i,2,0) > 0.0:
				outlist[2].append(SSNR(i,1,0)/SSNR(i,2,0))	          # variance
			else:
				outlist[2].append(0.0)
		for i in range(1,nn): outlist[3].append(SSNR(i,2,0))				  # number of points in the shell
		for i in range(1,nn): outlist[4].append(SSNR(i,3,0))				  # number of added Fourier points
		for i in range(1,nn): outlist[5].append(SSNR(i,0,0))				  # square of signal
		return [outlist, vol_ssnr]


def prepare_recons(data, symmetry, myid, main_node_half, half_start, step, index, finfo=None, npad = 2, mpi_comm=None):
	pass#IMPORTIMPORTIMPORT from random     import randint
	pass#IMPORTIMPORTIMPORT from utilities  import reduce_EMData_to_root
	pass#IMPORTIMPORTIMPORT from mpi        import mpi_barrier, MPI_COMM_WORLD
	pass#IMPORTIMPORTIMPORT from EMAN2 import Reconstructors

	if mpi_comm == None:
		mpi_comm = mpi.MPI_COMM_WORLD

	nx = data[0].get_xsize()

	fftvol_half = EMAN2_cppwrap.EMData()
	weight_half = EMAN2_cppwrap.EMData()
	half_params = {"size":nx, "npad":npad, "symmetry":symmetry, "fftvol":fftvol_half, "weight":weight_half}
	half = EMAN2_cppwrap.Reconstructors.get( "nn4", half_params )
	half.setup()

	group = -1
	for i in range(half_start, len(data), step):
		if(index >-1 ):  group = data[i].get_attr('group')
		if(group == index):
			# horatio active_refactoring Jy51i1EwmLD4tWZ9_00000_1
			# if( data[i].get_attr_default('active',1) == 1):
			# 	xform_proj = data[i].get_attr( "xform.projection" )
			# 	half.insert_slice(data[i], xform_proj )
			xform_proj = data[i].get_attr( "xform.projection" )
			half.insert_slice(data[i], xform_proj )

	if not(finfo is None):
		finfo.write( "begin reduce half\n" )
		finfo.flush()

	sparx_utilities.reduce_EMData_to_root(fftvol_half, myid, main_node_half, mpi_comm)
	sparx_utilities.reduce_EMData_to_root(weight_half, myid, main_node_half, mpi_comm)

	if not(finfo is None):
		finfo.write( "after reduce half\n" )
		finfo.flush()

	if myid == main_node_half:
		tmpid = random.randint(0, 1000000)
		fftvol_half_file = ("fftvol_half%d.hdf" % tmpid)
		weight_half_file = ("weight_half%d.hdf" % tmpid)
		fftvol_half.write_image(fftvol_half_file)
		weight_half.write_image(weight_half_file)
	mpi.mpi_barrier(mpi_comm)

	fftvol_half = None
	weight_half = None

	if myid == main_node_half:  return fftvol_half_file, weight_half_file

	return None, None

"""Multiline Comment11"""

def prepare_recons_ctf(nx, data, snr, symmetry, myid, main_node_half, half_start, step, finfo=None, npad = 2, mpi_comm=None, smearstep = 0.0):
	pass#IMPORTIMPORTIMPORT from random     import randint
	pass#IMPORTIMPORTIMPORT from utilities  import reduce_EMData_to_root
	pass#IMPORTIMPORTIMPORT from mpi        import mpi_barrier, MPI_COMM_WORLD
	pass#IMPORTIMPORTIMPORT from EMAN2 import Reconstructors

	if mpi_comm == None:
		mpi_comm = mpi.MPI_COMM_WORLD

	fftvol_half = EMAN2_cppwrap.EMData()

	if( smearstep > 0.0 ):
		#if myid == 0:  print "  Setting smear in prepare_recons_ctf"
		ns = 1
		smear = []
		for j in range(-ns,ns+1):
			if( j != 0):
				for i in range(-ns,ns+1):
					for k in range(-ns,ns+1):
						smear += [i*smearstep,j*smearstep,k*smearstep,1.0]
		# Deal with theta = 0.0 cases
		prj = []
		for i in range(-ns,ns+1):
			for k in range(-ns,ns+1):
				prj.append(i+k)
		for i in range(-2*ns,2*ns+1,1):
			smear += [i*smearstep,0.0,0.0,float(prj.count(i))]
		#if myid == 0:  print "  Smear  ",smear
		fftvol_half.set_attr("smear", smear)

	weight_half = EMAN2_cppwrap.EMData()
	half_params = {"size":nx, "npad":npad, "snr":snr, "sign":1, "symmetry":symmetry, "fftvol":fftvol_half, "weight":weight_half}
	half = EMAN2_cppwrap.Reconstructors.get( "nn4_ctf", half_params )
	half.setup()

	for i in range(half_start, len(data), step):
		xform_proj = data[i].get_attr( "xform.projection" )
		half.insert_slice(data[i], xform_proj )

	if not(finfo is None):
		finfo.write( "begin reduce half\n" )
		finfo.flush()

	sparx_utilities.reduce_EMData_to_root(fftvol_half, myid, main_node_half, mpi_comm)
	sparx_utilities.reduce_EMData_to_root(weight_half, myid, main_node_half, mpi_comm)

	if not(finfo is None):
		finfo.write( "after reduce half\n" )
		finfo.flush()

	if myid == main_node_half:
		tmpid = random.randint(0, 1000000) 
		fftvol_half_file = ("fftvol_half%d.hdf" % tmpid)
		weight_half_file = ("weight_half%d.hdf" % tmpid)
		fftvol_half.write_image(fftvol_half_file)
		weight_half.write_image(weight_half_file)
	mpi.mpi_barrier(mpi_comm)

	fftvol_half = None
	weight_half = None

	if myid == main_node_half:
		return fftvol_half_file, weight_half_file

	return None,None


def recons_from_fftvol(size, fftvol, weight, symmetry, npad = 2):
	pass#IMPORTIMPORTIMPORT from EMAN2 import Reconstructors

	params = {"size":size, "npad":npad, "symmetry":symmetry, "fftvol":fftvol, "weight":weight}
	r = EMAN2_cppwrap.Reconstructors.get("nn4", params)
	r.setup()
	dummy = r.finish(True)
	return fftvol


def recons_ctf_from_fftvol(size, fftvol, weight, snr, symmetry, weighting=1, npad = 2):
	pass#IMPORTIMPORTIMPORT from EMAN2 import Reconstructors

	params = {"size":size, "npad":npad, "snr":snr, "sign":1, "symmetry":symmetry, "fftvol":fftvol, "weight":weight, "weighting":weighting}
	r = EMAN2_cppwrap.Reconstructors.get("nn4_ctf", params)
	r.setup()
	dummy = r.finish(True)
	return fftvol

def get_image_size( imgdata, myid ):
	pass#IMPORTIMPORTIMPORT from mpi import mpi_gather, mpi_bcast, MPI_COMM_WORLD, MPI_INT
	nimg = len(imgdata)

	nimgs = mpi.mpi_gather( nimg, 1, mpi.MPI_INT, 1, mpi.MPI_INT, 0, mpi.MPI_COMM_WORLD )

	if myid==0:
		src = -1
		for i in range( len(nimgs) ):
			if int(nimgs[i]) > 0 :
				src = i
				break
		if src==-1:
			return 0
	else:
		src = -1

	size_src = mpi.mpi_bcast( src, 1, mpi.MPI_INT, 0, mpi.MPI_COMM_WORLD )

	if myid==int(size_src[0]):
		assert nimg > 0
		size = imgdata[0].get_xsize()
	else:
		size = -1

	nx = mpi.mpi_bcast( size, 1, mpi.MPI_INT, size_src[0], mpi.MPI_COMM_WORLD )
	return int(nx[0])


def rec3D_MPI(data, snr = 1.0, symmetry = "c1", mask3D = None, fsc_curve = None, \
		myid = 0, main_node = 0, rstep = 1.0, odd_start=0, eve_start=1, finfo=None, \
		index=-1, npad = 2, mpi_comm=None, smearstep = 0.0):
	'''
	  This function is to be called within an MPI program to do a reconstruction on a dataset kept 
	  in the memory, computes reconstruction and through odd-even, in order to get the resolution
	'''
	pass#IMPORTIMPORTIMPORT import os
	pass#IMPORTIMPORTIMPORT from statistics import fsc_mask
	pass#IMPORTIMPORTIMPORT from utilities  import model_blank, model_circle, get_image, send_EMData, recv_EMData
	pass#IMPORTIMPORTIMPORT from mpi        import mpi_comm_size, MPI_COMM_WORLD
	
	if mpi_comm == None:
		mpi_comm = mpi.MPI_COMM_WORLD
	
	nproc = mpi.mpi_comm_size(mpi_comm)

	if nproc==1:
		assert main_node==0
		main_node_odd = main_node
		main_node_eve = main_node
		main_node_all = main_node
	elif nproc==2:
		main_node_odd = main_node
		main_node_eve = (main_node+1)%2
		main_node_all = main_node

		tag_voleve     = 1000
		tag_fftvol_eve = 1001
		tag_weight_eve = 1002
	else:
		#spread CPUs between different nodes to save memory
		main_node_odd = main_node
		main_node_eve = (int(main_node)+nproc-1)%int(nproc)
		main_node_all = (int(main_node)+nproc//2)%int(nproc)

		tag_voleve     = 1000
		tag_fftvol_eve = 1001
		tag_weight_eve = 1002

		tag_fftvol_odd = 1003
		tag_weight_odd = 1004
		tag_volall     = 1005


	if index != -1 :
		grpdata = []
		for i in range(len(data)):
			if data[i].get_attr('group') == index:
				grpdata.append(data[i])
		imgdata = grpdata
	else:
		imgdata = data

	nx = get_image_size(imgdata, myid)
	if nx == 0:
		sparx_global_def.ERROR("Warning: no images were given for reconstruction, this usually means there is an empty group, returning empty volume", "rec3D", 0)
		return sparx_utilities.model_blank( 2, 2, 2 ), None

	fftvol_odd_file, weight_odd_file = prepare_recons_ctf(nx, imgdata, snr, symmetry, myid, main_node_odd, odd_start, 2, finfo, npad, mpi_comm=mpi_comm, smearstep = smearstep)
	fftvol_eve_file, weight_eve_file = prepare_recons_ctf(nx, imgdata, snr, symmetry, myid, main_node_eve, eve_start, 2, finfo, npad, mpi_comm=mpi_comm, smearstep = smearstep)
	del imgdata



	if nproc == 1:
		fftvol = sparx_utilities.get_image(fftvol_odd_file)
		weight = sparx_utilities.get_image(weight_odd_file)
		volodd = recons_ctf_from_fftvol(nx, fftvol, weight, snr, symmetry, npad = npad)

		fftvol = sparx_utilities.get_image(fftvol_eve_file)
		weight = sparx_utilities.get_image(weight_eve_file)
		voleve = recons_ctf_from_fftvol(nx, fftvol, weight, snr, symmetry, npad = npad)

		if( not mask3D ):
			nx = volodd.get_xsize()
			ny = volodd.get_ysize()
			nz = volodd.get_zsize()
			mask3D = sparx_utilities.model_circle(min(nx,ny,nz)//2 - 2, nx,ny,nz)
		fscdat = sparx_statistics.fsc_mask( volodd, voleve, mask3D, rstep, fsc_curve)
		del  volodd, voleve, mask3D

		fftvol = sparx_utilities.get_image( fftvol_odd_file )
		fftvol_tmp = sparx_utilities.get_image(fftvol_eve_file)
		fftvol += fftvol_tmp
		fftvol_tmp = None

		weight = sparx_utilities.get_image( weight_odd_file )
		weight_tmp = sparx_utilities.get_image(weight_eve_file)
		weight += weight_tmp
		weight_tmp = None

		volall = recons_ctf_from_fftvol(nx, fftvol, weight, snr, symmetry, npad = npad)
		os.system( "rm -f " + fftvol_odd_file + " " + weight_odd_file )
		os.system( "rm -f " + fftvol_eve_file + " " + weight_eve_file )

		return volall,fscdat

	if nproc == 2:
		if myid == main_node_odd:
			fftvol = sparx_utilities.get_image( fftvol_odd_file )
			weight = sparx_utilities.get_image( weight_odd_file )
			volodd = recons_ctf_from_fftvol(nx, fftvol, weight, snr, symmetry, npad = npad)
			voleve = sparx_utilities.recv_EMData(main_node_eve, tag_voleve, mpi_comm)
			
			if( not mask3D ):
				nx = volodd.get_xsize()
				ny = volodd.get_ysize()
				nz = volodd.get_zsize()
				mask3D = sparx_utilities.model_circle(min(nx,ny,nz)//2 - 2, nx,ny,nz)
			fscdat = sparx_statistics.fsc_mask( volodd, voleve, mask3D, rstep, fsc_curve)
			del  volodd, voleve, mask3D
		else:
			assert myid == main_node_eve
			fftvol = sparx_utilities.get_image( fftvol_eve_file )
			weight = sparx_utilities.get_image( weight_eve_file )
			voleve = recons_ctf_from_fftvol(nx, fftvol, weight, snr, symmetry, npad = npad)
			sparx_utilities.send_EMData(voleve, main_node_odd, tag_voleve, mpi_comm)

		if myid == main_node_odd:
			fftvol = sparx_utilities.get_image( fftvol_odd_file )
			fftvol_tmp = sparx_utilities.recv_EMData( main_node_eve, tag_fftvol_eve, mpi_comm)
			fftvol += fftvol_tmp
			fftvol_tmp = None

			weight = sparx_utilities.get_image( weight_odd_file )
			weight_tmp = sparx_utilities.recv_EMData( main_node_eve, tag_weight_eve, mpi_comm)
			weight += weight_tmp
			weight_tmp = None

			volall = recons_ctf_from_fftvol(nx, fftvol, weight, snr, symmetry, npad = npad)
			os.system( "rm -f " + fftvol_odd_file + " " + weight_odd_file )

			return volall,fscdat
		else:
			assert myid == main_node_eve
			fftvol = sparx_utilities.get_image( fftvol_eve_file )
			weight = sparx_utilities.get_image( weight_eve_file )
			sparx_utilities.send_EMData(fftvol, main_node_odd, tag_fftvol_eve, mpi_comm)
			sparx_utilities.send_EMData(weight, main_node_odd, tag_weight_eve, mpi_comm)
			os.system( "rm -f " + fftvol_eve_file + " " + weight_eve_file )
			return sparx_utilities.model_blank(nx,nx,nx), None

	# cases from all other number of processors situations
	if myid == main_node_odd:
		fftvol = sparx_utilities.get_image( fftvol_odd_file )
		sparx_utilities.send_EMData(fftvol, main_node_eve, tag_fftvol_odd, mpi_comm)

		if not(finfo is None):
			finfo.write("fftvol odd sent\n")
			finfo.flush()

		weight = sparx_utilities.get_image( weight_odd_file )
		sparx_utilities.send_EMData(weight, main_node_all, tag_weight_odd, mpi_comm)

		if not(finfo is None):
			finfo.write("weight odd sent\n")
			finfo.flush()

		volodd = recons_ctf_from_fftvol(nx, fftvol, weight, snr, symmetry, npad = npad)
		del fftvol, weight
		voleve = sparx_utilities.recv_EMData(main_node_eve, tag_voleve, mpi_comm)

		if( not mask3D ):
			nx = volodd.get_xsize()
			ny = volodd.get_ysize()
			nz = volodd.get_zsize()
			mask3D = sparx_utilities.model_circle(min(nx,ny,nz)//2 - 2, nx,ny,nz)

		fscdat = sparx_statistics.fsc_mask(volodd, voleve, mask3D, rstep, fsc_curve)
		del  volodd, voleve, mask3D
		volall = sparx_utilities.recv_EMData(main_node_all, tag_volall, mpi_comm)
		os.system( "rm -f " + fftvol_odd_file + " " + weight_odd_file )
		return volall, fscdat

	if myid == main_node_eve:
		ftmp = sparx_utilities.recv_EMData(main_node_odd, tag_fftvol_odd, mpi_comm)
		fftvol = sparx_utilities.get_image( fftvol_eve_file )
		EMAN2_cppwrap.Util.add_img( ftmp, fftvol )
		sparx_utilities.send_EMData(ftmp, main_node_all, tag_fftvol_eve, mpi_comm)
		del ftmp

		weight = sparx_utilities.get_image( weight_eve_file )
		sparx_utilities.send_EMData(weight, main_node_all, tag_weight_eve, mpi_comm)

		voleve = recons_ctf_from_fftvol(nx, fftvol, weight, snr, symmetry, npad = npad)
		sparx_utilities.send_EMData(voleve, main_node_odd, tag_voleve, mpi_comm)
		os.system( "rm -f " + fftvol_eve_file + " " + weight_eve_file );

		return sparx_utilities.model_blank(nx,nx,nx), None


	if myid == main_node_all:
		fftvol = sparx_utilities.recv_EMData(main_node_eve, tag_fftvol_eve, mpi_comm)
		if not(finfo is None):
			finfo.write( "fftvol odd received\n" )
			finfo.flush()

		weight = sparx_utilities.recv_EMData(main_node_odd, tag_weight_odd, mpi_comm)
		weight_tmp = sparx_utilities.recv_EMData(main_node_eve, tag_weight_eve, mpi_comm)
		EMAN2_cppwrap.Util.add_img( weight, weight_tmp )
		weight_tmp = None

		volall = recons_ctf_from_fftvol(nx, fftvol, weight, snr, symmetry, npad = npad)
		sparx_utilities.send_EMData(volall, main_node_odd, tag_volall, mpi_comm)

		return sparx_utilities.model_blank(nx,nx,nx),None

	return sparx_utilities.model_blank(nx,nx,nx),None


def rec3D_MPI_noCTF(data, symmetry = "c1", mask3D = None, fsc_curve = None, myid = 2, main_node = 0, \
		rstep = 1.0, odd_start=0, eve_start=1, finfo=None, index = -1, npad = 2, mpi_comm=None):
	'''
	  This function is to be called within an MPI program to do a reconstruction on a dataset kept in the memory 
	  Computes reconstruction and through odd-even, in order to get the resolution
	  if index > -1, projections should have attribute group set and only those whose group matches index will be used in the reconstruction
	    this is for multireference alignment
	'''
	pass#IMPORTIMPORTIMPORT import os
	pass#IMPORTIMPORTIMPORT from statistics import fsc_mask
	pass#IMPORTIMPORTIMPORT from utilities  import model_blank, get_image,send_EMData, recv_EMData
	pass#IMPORTIMPORTIMPORT from mpi        import mpi_comm_size, MPI_COMM_WORLD
	
	if mpi_comm == None:
		mpi_comm = mpi.MPI_COMM_WORLD
	
	nproc = mpi.mpi_comm_size(mpi_comm)

	if nproc==1:
		assert main_node==0
		main_node_odd = main_node
		main_node_eve = main_node
		main_node_all = main_node
	elif nproc==2:
		main_node_odd = main_node
		main_node_eve = (main_node+1)%2
		main_node_all = main_node

		tag_voleve     = 1000
		tag_fftvol_eve = 1001
		tag_weight_eve = 1002
	else:
		#spread CPUs between different nodes to save memory
		main_node_odd = main_node
		main_node_eve = (int(main_node)+nproc-1)%int(nproc)
		main_node_all = (int(main_node)+nproc//2)%int(nproc)

		tag_voleve     = 1000
		tag_fftvol_eve = 1001
		tag_weight_eve = 1002

		tag_fftvol_odd = 1003
		tag_weight_odd = 1004
		tag_volall     = 1005

	nx = data[0].get_xsize()

	fftvol_odd_file,weight_odd_file = prepare_recons(data, symmetry, myid, main_node_odd, odd_start, 2, index, finfo, npad, mpi_comm=mpi_comm)
	fftvol_eve_file,weight_eve_file = prepare_recons(data, symmetry, myid, main_node_eve, eve_start, 2, index, finfo, npad, mpi_comm=mpi_comm) 

	if nproc == 1:
		fftvol = sparx_utilities.get_image( fftvol_odd_file )
		weight = sparx_utilities.get_image( weight_odd_file )
		volodd = recons_from_fftvol(nx, fftvol, weight, symmetry, npad)

		fftvol = sparx_utilities.get_image( fftvol_eve_file )
		weight = sparx_utilities.get_image( weight_eve_file )
		voleve = recons_from_fftvol(nx, fftvol, weight, symmetry, npad)

		fscdat = sparx_statistics.fsc_mask( volodd, voleve, mask3D, rstep, fsc_curve)
		del  volodd, voleve

		fftvol = sparx_utilities.get_image( fftvol_odd_file )
		EMAN2_cppwrap.Util.add_img( fftvol, sparx_utilities.get_image(fftvol_eve_file) )

		weight = sparx_utilities.get_image( weight_odd_file )
		EMAN2_cppwrap.Util.add_img( weight, sparx_utilities.get_image(weight_eve_file) )

		volall = recons_from_fftvol(nx, fftvol, weight, symmetry, npad)
		os.system( "rm -f " + fftvol_odd_file + " " + weight_odd_file );
		os.system( "rm -f " + fftvol_eve_file + " " + weight_eve_file );
		return volall,fscdat

	if nproc == 2:
		if myid == main_node_odd:
			fftvol = sparx_utilities.get_image( fftvol_odd_file )
			weight = sparx_utilities.get_image( weight_odd_file )
			volodd = recons_from_fftvol(nx, fftvol, weight, symmetry, npad)
			voleve = sparx_utilities.recv_EMData(main_node_eve, tag_voleve, mpi_comm)
			fscdat = sparx_statistics.fsc_mask( volodd, voleve, mask3D, rstep, fsc_curve)
			del  volodd, voleve
		else:
			assert myid == main_node_eve
			fftvol = sparx_utilities.get_image( fftvol_eve_file )
			weight = sparx_utilities.get_image( weight_eve_file )
			voleve = recons_from_fftvol(nx, fftvol, weight, symmetry, npad)
			sparx_utilities.send_EMData(voleve, main_node_odd, tag_voleve, mpi_comm)

		if myid == main_node_odd:
			fftvol = sparx_utilities.get_image( fftvol_odd_file )
			fftvol_tmp = sparx_utilities.recv_EMData( main_node_eve, tag_fftvol_eve, mpi_comm)
			EMAN2_cppwrap.Util.add_img( fftvol, fftvol_tmp )
			fftvol_tmp = None

			weight = sparx_utilities.get_image( weight_odd_file )
			weight_tmp = sparx_utilities.recv_EMData( main_node_eve, tag_weight_eve, mpi_comm)
			EMAN2_cppwrap.Util.add_img( weight, weight_tmp )
			weight_tmp = None
			volall = recons_from_fftvol(nx, fftvol, weight, symmetry, npad)
			os.system( "rm -f " + fftvol_odd_file + " " + weight_odd_file );
			return volall,fscdat
		else:
			assert myid == main_node_eve
			fftvol = sparx_utilities.get_image( fftvol_eve_file )
			sparx_utilities.send_EMData(fftvol, main_node_odd, tag_fftvol_eve, mpi_comm)

			weight = sparx_utilities.get_image( weight_eve_file )
			sparx_utilities.send_EMData(weight, main_node_odd, tag_weight_eve, mpi_comm)
			os.system( "rm -f " + fftvol_eve_file + " " + weight_eve_file );
			return sparx_utilities.model_blank(nx,nx,nx), None
	# cases from all other number of processors situations
	if myid == main_node_odd:
		fftvol = sparx_utilities.get_image( fftvol_odd_file )
		sparx_utilities.send_EMData(fftvol, main_node_eve, tag_fftvol_odd, mpi_comm)

		if not(finfo is None):
			finfo.write("fftvol odd sent\n")
			finfo.flush()

		weight = sparx_utilities.get_image( weight_odd_file )
		sparx_utilities.send_EMData(weight, main_node_all, tag_weight_odd, mpi_comm)

		if not(finfo is None):
			finfo.write("weight odd sent\n")
			finfo.flush()

		volodd = recons_from_fftvol(nx, fftvol, weight, symmetry, npad)
		del fftvol, weight
		voleve = sparx_utilities.recv_EMData(main_node_eve, tag_voleve, mpi_comm)
		fscdat = sparx_statistics.fsc_mask(volodd, voleve, mask3D, rstep, fsc_curve)
		del  volodd, voleve
		volall = sparx_utilities.recv_EMData(main_node_all, tag_volall, mpi_comm)
		os.system( "rm -f " + fftvol_odd_file + " " + weight_odd_file );
		return volall,fscdat

	if myid == main_node_eve:
		ftmp = sparx_utilities.recv_EMData(main_node_odd, tag_fftvol_odd, mpi_comm)
		fftvol = sparx_utilities.get_image( fftvol_eve_file )
		EMAN2_cppwrap.Util.add_img( ftmp, fftvol )
		sparx_utilities.send_EMData(ftmp, main_node_all, tag_fftvol_eve, mpi_comm)
		del ftmp

		weight = sparx_utilities.get_image( weight_eve_file )
		sparx_utilities.send_EMData(weight, main_node_all, tag_weight_eve, mpi_comm)

		voleve = recons_from_fftvol(nx, fftvol, weight, symmetry, npad)
		sparx_utilities.send_EMData(voleve, main_node_odd, tag_voleve, mpi_comm)
		os.system( "rm -f " + fftvol_eve_file + " " + weight_eve_file );

		return sparx_utilities.model_blank(nx,nx,nx), None


	if myid == main_node_all:
		fftvol = sparx_utilities.recv_EMData(main_node_eve, tag_fftvol_eve, mpi_comm)
		if not(finfo is None):
			finfo.write( "fftvol odd received\n" )
			finfo.flush()

		weight = sparx_utilities.recv_EMData(main_node_odd, tag_weight_odd, mpi_comm)
		weight_tmp = sparx_utilities.recv_EMData(main_node_eve, tag_weight_eve, mpi_comm)
		EMAN2_cppwrap.Util.add_img( weight, weight_tmp )
		weight_tmp = None

		volall = recons_from_fftvol(nx, fftvol, weight, symmetry, npad)
		sparx_utilities.send_EMData(volall, main_node_odd, tag_volall, mpi_comm)

		return sparx_utilities.model_blank(nx,nx,nx),None


	return sparx_utilities.model_blank(nx,nx,nx),None
	
def prepare_recons_ctf_two_chunks(nx,data,snr,symmetry,myid,main_node_half,chunk_ID,finfo=None,npad=2,mpi_comm=None,smearstep = 0.0):
	pass#IMPORTIMPORTIMPORT from random     import randint
	pass#IMPORTIMPORTIMPORT from utilities  import reduce_EMData_to_root
	pass#IMPORTIMPORTIMPORT from mpi        import mpi_barrier, MPI_COMM_WORLD
	pass#IMPORTIMPORTIMPORT from EMAN2 import Reconstructors

	if mpi_comm == None:
		mpi_comm = mpi.MPI_COMM_WORLD

	fftvol_half = EMAN2_cppwrap.EMData()

	if( smearstep > 0.0 ):
		#if myid == 0:  print "  Setting smear in prepare_recons_ctf"
		ns = 1
		smear = []
		for j in range(-ns,ns+1):
			if( j != 0):
				for i in range(-ns,ns+1):
					for k in range(-ns,ns+1):
						smear += [i*smearstep,j*smearstep,k*smearstep,1.0]
		# Deal with theta = 0.0 cases
		prj = []
		for i in range(-ns,ns+1):
			for k in range(-ns,ns+1):
				prj.append(i+k)
		for i in range(-2*ns,2*ns+1,1):
			 smear += [i*smearstep,0.0,0.0,float(prj.count(i))]
		#if myid == 0:  print "  Smear  ",smear
		fftvol_half.set_attr("smear", smear)

	weight_half = EMAN2_cppwrap.EMData()
	half_params = {"size":nx, "npad":npad, "snr":snr, "sign":1, "symmetry":symmetry, "fftvol":fftvol_half, "weight":weight_half}
	half = EMAN2_cppwrap.Reconstructors.get( "nn4_ctf", half_params )
	half.setup()
	for i in range(len(data)):
		if data[i].get_attr("chunk_id") == chunk_ID:
			xform_proj = data[i].get_attr( "xform.projection" )
			half.insert_slice(data[i], xform_proj )
	if not(finfo is None):
		finfo.write( "begin reduce half\n" )
		finfo.flush()

	sparx_utilities.reduce_EMData_to_root(fftvol_half, myid, main_node_half, mpi_comm)
	sparx_utilities.reduce_EMData_to_root(weight_half, myid, main_node_half, mpi_comm)

	if not(finfo is None):
		finfo.write( "after reduce half\n" )
		finfo.flush()

	if myid == main_node_half:
		tmpid = random.randint(0, 1000000) 
		fftvol_half_file = ("fftvol_half%d.hdf" % tmpid)
		weight_half_file = ("weight_half%d.hdf" % tmpid)
		fftvol_half.write_image(fftvol_half_file)
		weight_half.write_image(weight_half_file)
	mpi.mpi_barrier(mpi_comm)

	fftvol_half = None
	weight_half = None

	if myid == main_node_half:
		return fftvol_half_file, weight_half_file

	return None,None
	
def rec3D_two_chunks_MPI(data, snr = 1.0, symmetry = "c1", mask3D = None, fsc_curve = None, \
		myid = 0, main_node = 0, rstep = 1.0, finfo=None, \
		index=-1, npad = 2, mpi_comm=None, smearstep = 0.0):
	'''
	  This function is to be called within an MPI program to do a reconstruction on a dataset kept 
	  in the memory, computes reconstruction and through odd-even, in order to get the resolution
	'''
	pass#IMPORTIMPORTIMPORT import os
	pass#IMPORTIMPORTIMPORT from statistics import fsc_mask
	pass#IMPORTIMPORTIMPORT from utilities  import model_blank, model_circle, get_image, send_EMData, recv_EMData
	pass#IMPORTIMPORTIMPORT from mpi        import mpi_comm_size, MPI_COMM_WORLD
	
	if mpi_comm == None:
		mpi_comm = mpi.MPI_COMM_WORLD
	
	nproc = mpi.mpi_comm_size(mpi_comm)

	if nproc==1:
		assert main_node==0
		main_node_odd = main_node
		main_node_eve = main_node
		main_node_all = main_node
	elif nproc==2:
		main_node_odd = main_node
		main_node_eve = (main_node+1)%2
		main_node_all = main_node

		tag_voleve     = 1000
		tag_fftvol_eve = 1001
		tag_weight_eve = 1002
	else:
		#spread CPUs between different nodes to save memory
		main_node_odd = main_node
		main_node_eve = (int(main_node)+nproc-1)%int(nproc)
		main_node_all = (int(main_node)+nproc//2)%int(nproc)

		tag_voleve     = 1000
		tag_fftvol_eve = 1001
		tag_weight_eve = 1002

		tag_fftvol_odd = 1003
		tag_weight_odd = 1004
		tag_volall     = 1005


	if index != -1 :
		grpdata = []
		for i in range(len(data)):
			if data[i].get_attr('group') == index:
				grpdata.append(data[i])
		imgdata = grpdata
	else:
		imgdata = data

	nx = get_image_size(imgdata, myid)
	if nx == 0:
		sparx_global_def.ERROR("Warning: no images were given for reconstruction, this usually means there is an empty group, returning empty volume", "rec3D", 0)
		return sparx_utilities.model_blank( 2, 2, 2 ), None

	fftvol_odd_file,weight_odd_file = prepare_recons_ctf_two_chunks(nx, imgdata, snr, symmetry, myid, main_node_odd, 0, finfo, npad, mpi_comm=mpi_comm, smearstep = smearstep)
	fftvol_eve_file,weight_eve_file = prepare_recons_ctf_two_chunks(nx, imgdata, snr, symmetry, myid, main_node_eve, 1, finfo, npad, mpi_comm=mpi_comm, smearstep = smearstep)
	del imgdata

	if nproc == 1:
		fftvol = sparx_utilities.get_image(fftvol_odd_file)
		weight = sparx_utilities.get_image(weight_odd_file)
		volodd = recons_ctf_from_fftvol(nx, fftvol, weight, snr, symmetry, npad = npad)

		fftvol = sparx_utilities.get_image(fftvol_eve_file)
		weight = sparx_utilities.get_image(weight_eve_file)
		voleve = recons_ctf_from_fftvol(nx, fftvol, weight, snr, symmetry, npad = npad)

		if( not mask3D ):
			nx = volodd.get_xsize()
			ny = volodd.get_ysize()
			nz = volodd.get_zsize()
			mask3D = sparx_utilities.model_circle(min(nx,ny,nz)//2 - 2, nx,ny,nz)
		fscdat = sparx_statistics.fsc_mask( volodd, voleve, mask3D, rstep, fsc_curve)
		del  volodd, voleve, mask3D

		fftvol = sparx_utilities.get_image( fftvol_odd_file )
		fftvol_tmp = sparx_utilities.get_image(fftvol_eve_file)
		fftvol += fftvol_tmp
		fftvol_tmp = None

		weight = sparx_utilities.get_image( weight_odd_file )
		weight_tmp = sparx_utilities.get_image(weight_eve_file)
		weight += weight_tmp
		weight_tmp = None

		volall = recons_ctf_from_fftvol(nx, fftvol, weight, snr, symmetry, npad = npad)
		os.system( "rm -f " + fftvol_odd_file + " " + weight_odd_file )
		os.system( "rm -f " + fftvol_eve_file + " " + weight_eve_file )

		return volall,fscdat

	if nproc == 2:
		if myid == main_node_odd:
			fftvol = sparx_utilities.get_image( fftvol_odd_file )
			weight = sparx_utilities.get_image( weight_odd_file )
			volodd = recons_ctf_from_fftvol(nx, fftvol, weight, snr, symmetry, npad = npad)
			voleve = sparx_utilities.recv_EMData(main_node_eve, tag_voleve, mpi_comm)
			
			if( not mask3D ):
				nx = volodd.get_xsize()
				ny = volodd.get_ysize()
				nz = volodd.get_zsize()
				mask3D = sparx_utilities.model_circle(min(nx,ny,nz)//2 - 2, nx,ny,nz)
			fscdat = sparx_statistics.fsc_mask( volodd, voleve, mask3D, rstep, fsc_curve)
			del  volodd, voleve, mask3D
		else:
			assert myid == main_node_eve
			fftvol = sparx_utilities.get_image( fftvol_eve_file )
			weight = sparx_utilities.get_image( weight_eve_file )
			voleve = recons_ctf_from_fftvol(nx, fftvol, weight, snr, symmetry, npad = npad)
			sparx_utilities.send_EMData(voleve, main_node_odd, tag_voleve, mpi_comm)

		if myid == main_node_odd:
			fftvol = sparx_utilities.get_image( fftvol_odd_file )
			fftvol_tmp = sparx_utilities.recv_EMData( main_node_eve, tag_fftvol_eve, mpi_comm)
			fftvol += fftvol_tmp
			fftvol_tmp = None

			weight = sparx_utilities.get_image( weight_odd_file )
			weight_tmp = sparx_utilities.recv_EMData( main_node_eve, tag_weight_eve, mpi_comm)
			weight += weight_tmp
			weight_tmp = None

			volall = recons_ctf_from_fftvol(nx, fftvol, weight, snr, symmetry, npad = npad)
			os.system( "rm -f " + fftvol_odd_file + " " + weight_odd_file )

			return volall,fscdat
		else:
			assert myid == main_node_eve
			fftvol = sparx_utilities.get_image( fftvol_eve_file )
			weight = sparx_utilities.get_image( weight_eve_file )
			sparx_utilities.send_EMData(fftvol, main_node_odd, tag_fftvol_eve, mpi_comm)
			sparx_utilities.send_EMData(weight, main_node_odd, tag_weight_eve, mpi_comm)
			os.system( "rm -f " + fftvol_eve_file + " " + weight_eve_file )
			return sparx_utilities.model_blank(nx,nx,nx), None

	# cases from all other number of processors situations
	if myid == main_node_odd:
		fftvol = sparx_utilities.get_image( fftvol_odd_file )
		sparx_utilities.send_EMData(fftvol, main_node_eve, tag_fftvol_odd, mpi_comm)

		if not(finfo is None):
			finfo.write("fftvol odd sent\n")
			finfo.flush()

		weight = sparx_utilities.get_image( weight_odd_file )
		sparx_utilities.send_EMData(weight, main_node_all, tag_weight_odd, mpi_comm)

		if not(finfo is None):
			finfo.write("weight odd sent\n")
			finfo.flush()

		volodd = recons_ctf_from_fftvol(nx, fftvol, weight, snr, symmetry, npad = npad)
		del fftvol, weight
		voleve = sparx_utilities.recv_EMData(main_node_eve, tag_voleve, mpi_comm)

		if( not mask3D ):
			nx = volodd.get_xsize()
			ny = volodd.get_ysize()
			nz = volodd.get_zsize()
			mask3D = sparx_utilities.model_circle(min(nx,ny,nz)//2 - 2, nx,ny,nz)

		fscdat = sparx_statistics.fsc_mask(volodd, voleve, mask3D, rstep, fsc_curve)
		del  volodd, voleve, mask3D
		volall = sparx_utilities.recv_EMData(main_node_all, tag_volall, mpi_comm)
		os.system( "rm -f " + fftvol_odd_file + " " + weight_odd_file )
		return volall, fscdat

	if myid == main_node_eve:
		ftmp = sparx_utilities.recv_EMData(main_node_odd, tag_fftvol_odd, mpi_comm)
		fftvol = sparx_utilities.get_image( fftvol_eve_file )
		EMAN2_cppwrap.Util.add_img( ftmp, fftvol )
		sparx_utilities.send_EMData(ftmp, main_node_all, tag_fftvol_eve, mpi_comm)
		del ftmp

		weight = sparx_utilities.get_image( weight_eve_file )
		sparx_utilities.send_EMData(weight, main_node_all, tag_weight_eve, mpi_comm)

		voleve = recons_ctf_from_fftvol(nx, fftvol, weight, snr, symmetry, npad = npad)
		sparx_utilities.send_EMData(voleve, main_node_odd, tag_voleve, mpi_comm)
		os.system( "rm -f " + fftvol_eve_file + " " + weight_eve_file );

		return sparx_utilities.model_blank(nx,nx,nx), None


	if myid == main_node_all:
		fftvol = sparx_utilities.recv_EMData(main_node_eve, tag_fftvol_eve, mpi_comm)
		if not(finfo is None):
			finfo.write( "fftvol odd received\n" )
			finfo.flush()

		weight = sparx_utilities.recv_EMData(main_node_odd, tag_weight_odd, mpi_comm)
		weight_tmp = sparx_utilities.recv_EMData(main_node_eve, tag_weight_eve, mpi_comm)
		EMAN2_cppwrap.Util.add_img( weight, weight_tmp )
		weight_tmp = None

		volall = recons_ctf_from_fftvol(nx, fftvol, weight, snr, symmetry, npad = npad)
		sparx_utilities.send_EMData(volall, main_node_odd, tag_volall, mpi_comm)

		return sparx_utilities.model_blank(nx,nx,nx),None

	return sparx_utilities.model_blank(nx,nx,nx),None

