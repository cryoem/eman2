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
import sparx_fundamentals
import sparx_morphology
import mpi
import numpy
import sparx_utilities
pass#IMPORTIMPORTIMPORT import EMAN2
pass#IMPORTIMPORTIMPORT import EMAN2_cppwrap
pass#IMPORTIMPORTIMPORT import filter
pass#IMPORTIMPORTIMPORT import fundamentals
pass#IMPORTIMPORTIMPORT import global_def
pass#IMPORTIMPORTIMPORT import math
pass#IMPORTIMPORTIMPORT import morphology
pass#IMPORTIMPORTIMPORT import mpi
pass#IMPORTIMPORTIMPORT import numpy
pass#IMPORTIMPORTIMPORT import utilities
from builtins import range
pass#IMPORTIMPORTIMPORT from global_def import *

def filt_tophatb(e, freql, freqh, pad = False):
	"""
		Name
			filt_tophatb - top-hat band-pass Fourier filter (truncation of a Fourier series)
		Input
			e: input image (can be either real or Fourier)
			freql: low-end frequency of the filter pass-band
			freqh: high-end frequency of the filter pass-band
			pad: logical flag specifying whether before filtering the image should be padded with zeroes in real space to twice the size (this helps avoiding aliasing artifacts). (Default pad = False).
			All frequencies are in absolute frequency units fa and their valid range is [0:0.5].
		Output
			filtered image. Output image is real when input image is real or Fourier when input image is Fourier
	"""
	pass#IMPORTIMPORTIMPORT from EMAN2 import Processor
	params = {"filter_type" : EMAN2_cppwrap.Processor.fourier_filter_types.TOP_HAT_BAND_PASS,
		  "low_cutoff_frequency" : freql, "high_cutoff_frequency" : freqh, "dopad" : pad}
	return EMAN2_cppwrap.Processor.EMFourierFilter(e, params)
    
def filt_gaussl(e, sigma, pad = False):
	"""
		Name
			filt_gaussl - Gaussian low-pass Fourier filter
		Input
			e: input image (can be either real or Fourier)
			sigma: standard deviation of the Gaussian function in absolute frequency units fa.
			pad: logical flag specifying whether before filtering the image should be padded with zeroes in real space to twice the size (this helps avoiding aliasing artifacts). (Default pad = False).
			Note: for sigma = 0.5 (the Nyquist frequency) the value of the filter at the maximum frequency is G(fN)=1e=0.61.
		Output
			filtered image. Output image is real when input image is real or Fourier when input image is Fourier
	"""
	pass#IMPORTIMPORTIMPORT from EMAN2 import Processor
	params = {"filter_type" : EMAN2_cppwrap.Processor.fourier_filter_types.GAUSS_LOW_PASS,
		  "cutoff_abs" : sigma, "dopad" : pad}
	return EMAN2_cppwrap.Processor.EMFourierFilter(e, params)
    
def filt_gaussinv(e, sigma, pad = False):
	"""
		Name
			filt_gaussinv - inverse Gaussian (high-pass) Fourier filter (division by a Gaussian function in Fourier space)
		Input
			e: input image (can be either real or Fourier)
			sigma: standard deviation of the Gaussian function in absolute frequency units fa.
			pad: logical flag specifying whether before filtering the image should be padded with zeroes in real space to twice the size (this helps avoiding aliasing artifacts). (Default pad = False).
			Note: for sigma = 0.5 (the Nyquist frequency) the value of the filter at the maximum frequency is G(fN)=e=1.65.
		Output
			filtered image. Output image is real when input image is real or Fourier when input image is Fourier
	"""
	pass#IMPORTIMPORTIMPORT from EMAN2 import Processor
	params = {"filter_type" : EMAN2_cppwrap.Processor.fourier_filter_types.GAUSS_INVERSE,
		  "cutoff_abs" : sigma, "dopad" : pad}
	return EMAN2_cppwrap.Processor.EMFourierFilter(e, params)

def filt_gaussh(e, sigma, pad = False):
	"""
		Name
			filt_gaussh - Gaussian high-pass Fourier filter
		Input
			e: input image (can be either real or Fourier)
			sigma: standard deviation of the Gaussian function in absolute frequency units fa.
			pad: logical flag specifying whether before filtering the image should be padded with zeroes in real space to twice the size (this helps avoiding aliasing artifacts). (Default pad = False).
			Note: for sigma = 0.5 (the Nyquist frequency) the value of the filter at the maximum frequency is 1-G(0)=1-1e=0.39
		Output
			filtered image. Output image is real when input image is real or Fourier when input image is Fourier
	"""
	pass#IMPORTIMPORTIMPORT from EMAN2 import Processor
	params = {"filter_type" : EMAN2_cppwrap.Processor.fourier_filter_types.GAUSS_HIGH_PASS,
		  "cutoff_abs" : sigma, "dopad" : pad}
	return EMAN2_cppwrap.Processor.EMFourierFilter(e, params)

def filt_tanl(e, freq, fall_off, pad = False):
	"""
		Name
			filt_tanl - hyperbolic tangent low-pass Fourier filter
		Input
			e: input image (can be either real or Fourier)
			freq: stop-band frequency
			fall_off: fall off of the filter
			pad: logical flag specifying whether before filtering the image should be padded with zeroes in real space to twice the size (this helps avoiding aliasing artifacts). (Default pad = False).
			All frequencies are in absolute frequency units fa and their valid range is [0:0.5].
		Output
			filtered image. Output image is real when input image is real or Fourier when input image is Fourier
	"""
	pass#IMPORTIMPORTIMPORT from EMAN2 import Processor
	params = {"filter_type" : EMAN2_cppwrap.Processor.fourier_filter_types.TANH_LOW_PASS,
		  "cutoff_abs" : freq, "fall_off": fall_off, "dopad" : pad}
	return EMAN2_cppwrap.Processor.EMFourierFilter(e, params)
    
def filt_table(e, table):
	"""
		Name
			filt_table - filter image using a user-provided (as a list) of Fourier filter values
		Input
			e: input image (real or Fourier)
			table: user-provided 1-D table of filter values.
		Output
			image Fourier-filtered using coefficients in table (real or Fourier, according the input image)
		Options
			fast: use the fast method; may combust certain computers.
			huge: gobble memory; there is plenty of it, anyway.
	"""
	pass#IMPORTIMPORTIMPORT from EMAN2 import Processor
	params = {"filter_type" : EMAN2_cppwrap.Processor.fourier_filter_types.RADIAL_TABLE,
			"table" : table}
	return EMAN2_cppwrap.Processor.EMFourierFilter(e, params)

def filt_ctf(img, ctf, dopad=True, sign=1, binary = 0):
	"""
		Name
			filt_ctf - apply Contrast Transfer Function (CTF) to an image in Fourier space
		Input
			image: input image, it can be real or complex
			ctf: an CTF object, please see CTF_info for description.
			pad: apply padding with zeroes in real space to twice the size before CTF application (Default is True, if set to False, no padding, if input image is Fourier, the flag has no effect).
			sign: sign of the CTF. If cryo data had inverted contrast, i.e., images were multiplied by -1 and particle projections appear bright on dark background, it has to be set to -1). (Default is 1).
			binary: phase flipping if set to 1 (default is 0).
		Output
			image multiplied in Fourier space by the CTF, the output image has the same format as the input image.
	"""
	pass#IMPORTIMPORTIMPORT from EMAN2 import Processor
	assert img.get_ysize() > 1
	dict        = ctf.to_dict()
	dz          = dict["defocus"]
	cs          = dict["cs"]
	voltage     = dict["voltage"]
	pixel_size  = dict["apix"]
	b_factor    = dict["bfactor"]
	ampcont     = dict["ampcont"]
	dza         = dict["dfdiff"]
	azz         = dict["dfang"]

	if dopad and not img.is_complex():  ip = 1
	else:                               ip = 0

	params = {"filter_type": EMAN2_cppwrap.Processor.fourier_filter_types.CTF_,
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
	tmp = EMAN2_cppwrap.Processor.EMFourierFilter(img, params)
	tmp.set_attr_dict({"ctf":ctf})
	return tmp

def filterlocal(ui, vi, m, falloff, myid, main_node, number_of_proc):
	pass#IMPORTIMPORTIMPORT from mpi 	  	  import mpi_init, mpi_comm_size, mpi_comm_rank, MPI_COMM_WORLD
	pass#IMPORTIMPORTIMPORT from mpi 	  	  import mpi_reduce, mpi_bcast, mpi_barrier, mpi_gatherv, mpi_send, mpi_recv
	pass#IMPORTIMPORTIMPORT from mpi 	  	  import MPI_SUM, MPI_FLOAT, MPI_INT
	pass#IMPORTIMPORTIMPORT from utilities import bcast_number_to_all, bcast_list_to_all, model_blank, bcast_EMData_to_all, reduce_EMData_to_root
	pass#IMPORTIMPORTIMPORT from morphology import threshold_outside
	pass#IMPORTIMPORTIMPORT from filter import filt_tanl
	pass#IMPORTIMPORTIMPORT from fundamentals import fft, fftip

	if(myid == main_node):

		nx = vi.get_xsize()
		ny = vi.get_ysize()
		nz = vi.get_zsize()
		#  Round all resolution numbers to two digits
		for x in range(nx):
			for y in range(ny):
				for z in range(nz):
					ui.set_value_at_fast( x,y,z, round(ui.get_value_at(x,y,z), 2) )
		dis = [nx,ny,nz]
	else:
		falloff = 0.0
		radius  = 0
		dis = [0,0,0]
	falloff = sparx_utilities.bcast_number_to_all(falloff, main_node)
	dis = sparx_utilities.bcast_list_to_all(dis, myid, source_node = main_node)

	if(myid != main_node):
		nx = int(dis[0])
		ny = int(dis[1])
		nz = int(dis[2])

		vi = sparx_utilities.model_blank(nx,ny,nz)
		ui = sparx_utilities.model_blank(nx,ny,nz)

	sparx_utilities.bcast_EMData_to_all(vi, myid, main_node)
	sparx_utilities.bcast_EMData_to_all(ui, myid, main_node)

	sparx_fundamentals.fftip(vi)  #  volume to be filtered

	st = EMAN2_cppwrap.Util.infomask(ui, m, True)


	filteredvol = sparx_utilities.model_blank(nx,ny,nz)
	cutoff = max(st[2] - 0.01,0.0)
	while(cutoff < st[3] ):
		cutoff = round(cutoff + 0.01, 2)
		#if(myid == main_node):  print  cutoff,st
		pt = EMAN2_cppwrap.Util.infomask( sparx_morphology.threshold_outside(ui, cutoff - 0.00501, cutoff + 0.005), m, True)  # Ideally, one would want to check only slices in question...
		if(pt[0] != 0.0):
			#print cutoff,pt[0]
			vovo = sparx_fundamentals.fft( filt_tanl(vi, cutoff, falloff) )
			for z in range(myid, nz, number_of_proc):
				for x in range(nx):
					for y in range(ny):
						if(m.get_value_at(x,y,z) > 0.5):
							if(round(ui.get_value_at(x,y,z),2) == cutoff):
								filteredvol.set_value_at_fast(x,y,z,vovo.get_value_at(x,y,z))

	mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
	sparx_utilities.reduce_EMData_to_root(filteredvol, myid, main_node, mpi.MPI_COMM_WORLD)
	return filteredvol
