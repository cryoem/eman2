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
import EMAN2db
import sparx_alignment
import sparx_applications
import configparser
import copy
import sparx_filter
import sparx_fundamentals
import sparx_global_def
import logging
import sparx_morphology
import mpi
import numpy
import numpy.random
import os
import sparx_pixel_error
import sparx_projection
import random
import random as rdq
import sparx_reconstruction
import scipy.stats
import string
import sys
import time
import sparx_utilities
pass#IMPORTIMPORTIMPORT import EMAN2
pass#IMPORTIMPORTIMPORT import EMAN2_cppwrap
pass#IMPORTIMPORTIMPORT import EMAN2db
pass#IMPORTIMPORTIMPORT import alignment
pass#IMPORTIMPORTIMPORT import applications
pass#IMPORTIMPORTIMPORT import array
pass#IMPORTIMPORTIMPORT import configparser
pass#IMPORTIMPORTIMPORT import copy
pass#IMPORTIMPORTIMPORT import development
pass#IMPORTIMPORTIMPORT import filter
pass#IMPORTIMPORTIMPORT import fundamentals
pass#IMPORTIMPORTIMPORT import global_def
pass#IMPORTIMPORTIMPORT import logging
pass#IMPORTIMPORTIMPORT import math
pass#IMPORTIMPORTIMPORT import morphology
pass#IMPORTIMPORTIMPORT import mpi
pass#IMPORTIMPORTIMPORT import numpy
pass#IMPORTIMPORTIMPORT import numpy.random
pass#IMPORTIMPORTIMPORT import os
pass#IMPORTIMPORTIMPORT import pickle
pass#IMPORTIMPORTIMPORT import pixel_error
pass#IMPORTIMPORTIMPORT import projection
pass#IMPORTIMPORTIMPORT import random
pass#IMPORTIMPORTIMPORT import random as rdq
pass#IMPORTIMPORTIMPORT import reconstruction
pass#IMPORTIMPORTIMPORT import scipy.stats
pass#IMPORTIMPORTIMPORT import statistics
pass#IMPORTIMPORTIMPORT import string
pass#IMPORTIMPORTIMPORT import sys
pass#IMPORTIMPORTIMPORT import time
pass#IMPORTIMPORTIMPORT import types
pass#IMPORTIMPORTIMPORT import utilities
from future import standard_library
standard_library.install_aliases()
from builtins import range
from builtins import object
pass#IMPORTIMPORTIMPORT from global_def import *

def sum_oe(data, mode = "a", CTF = False, ctf_2_sum = None, ctf_eo_sum = False, return_params = False):
	"""
		Calculate average of an image series
		mode - "a": use current alignment parameters
		CTF  - if True, use CTF for calculations of the average.
		In addition, calculate odd and even sums, these are not divided by the ctf^2
		If ctf_2_sum not provided, sum of ctf^2 will be calculated and returned
		If ctf_eo_sum is True, then compute ctf^2 in odd and even form
		If return_params is True, then return ali2d.xform parameters
	"""
	pass#IMPORTIMPORTIMPORT from utilities    import    model_blank, get_params2D, same_ctf
	pass#IMPORTIMPORTIMPORT from fundamentals import    rot_shift2D, fft
	pass#IMPORTIMPORTIMPORT from copy import deepcopy
	if CTF: sparx_global_def.ERROR("This function was disabled as it does not treat astigmatism properly","sum_oe",1)
	n      = len(data)
	if return_params: params_list = [None]*n
	if CTF:
		origin_size = data[0].get_xsize()
		ave1   = EMAN2_cppwrap.EMData(origin_size, origin_size, 1, False) 
		ave2   = EMAN2_cppwrap.EMData(origin_size, origin_size, 1, False)
		pass#IMPORTIMPORTIMPORT from morphology import ctf_img
		ctf_2_sumo = EMAN2_cppwrap.EMData(origin_size, origin_size, 1, False)
		ctf_2_sume = EMAN2_cppwrap.EMData(origin_size, origin_size, 1, False)

		if data[0].get_attr_default('ctf_applied', 1) == 1:  sparx_global_def.ERROR("data cannot be ctf-applied", "sum_oe", 1)
		if ctf_2_sum:  get_ctf2 = False
		else:          get_ctf2 = True
		for im in range(n):
			current_ctf = data[im].get_attr("ctf")
			if im == 0: 
				myctf = copy.deepcopy(current_ctf)
				ctt = sparx_morphology.ctf_img(origin_size, myctf)
			else:
				if not sparx_utilities.same_ctf(current_ctf, myctf):
					myctf = copy.deepcopy(current_ctf)
					ctt   = sparx_morphology.ctf_img(origin_size, myctf)
			if mode == "a":
				alpha, sx, sy, mirror, scale = sparx_utilities.get_params2D(data[im])
				ima = sparx_fundamentals.rot_shift2D(data[im], alpha, sx, sy, mirror, scale, "quadratic")
				if return_params:
					params_list[im]= [alpha, sx, sy, mirror, scale]
			else:
				ima = data[im]
				if return_params:
					alpha, sx, sy, mirror, scale = sparx_utilities.get_params2D(data[im])
					params_list[im]= [alpha, sx, sy, mirror, scale]
			ima = sparx_fundamentals.fft(ima)
			EMAN2_cppwrap.Util.mul_img(ima, ctt)
			if im%2 == 0:	
				EMAN2_cppwrap.Util.add_img(ave1, ima)
				EMAN2_cppwrap.Util.add_img2(ctf_2_sume, ctt)
			else:	        
				EMAN2_cppwrap.Util.add_img(ave2, ima)
				EMAN2_cppwrap.Util.add_img2(ctf_2_sumo, ctt)
	else:
		nx     = data[0].get_xsize()
		ny     = data[0].get_ysize()
		ave1   = sparx_utilities.model_blank(nx, ny, 1) 
		ave2   = sparx_utilities.model_blank(nx, ny, 1)
		for im in range(n):
			if mode == "a":
				alpha, sx, sy, mirror, scale = sparx_utilities.get_params2D(data[im])
				ima = sparx_fundamentals.rot_shift2D(data[im], alpha, sx, sy, mirror, scale, "quadratic")
				if return_params:
					params_list[im]= [alpha, sx, sy, mirror, scale]
			else:
				ima = data[im]
				if return_params:
					alpha, sx, sy, mirror, scale = sparx_utilities.get_params2D(data[im])
					params_list[im]= [alpha, sx, sy, mirror, scale]
			if im%2 == 0:	EMAN2_cppwrap.Util.add_img(ave1, ima)
			else:	        EMAN2_cppwrap.Util.add_img(ave2, ima)

	if  CTF:
		if get_ctf2:
			if not ctf_eo_sum:# Old usage
				ctf_2_sum  = EMAN2_cppwrap.Util.addn_img(ctf_2_sume, ctf_2_sumo)
				return sparx_fundamentals.fft(ave1), sparx_fundamentals.fft(ave2), ctf_2_sum
			else: # return Fourier images
				if return_params: return ave1, ave2, ctf_2_sume, ctf_2_sumo, params_list 
				else: return ave1, ave2, ctf_2_sume, ctf_2_sumo
		else:
			if return_params: return  ave1, ave2, params_list
			else: return ave1, ave2
	else:
		if not return_params: return  ave1, ave2
		else: return  ave1, ave2, params_list
		
def ave_var(data, mode = "a", listID=None):
	"""
		Calculate average and variance of a 2D or 3D image series
		with optional application of orientation parameters
		data can be either in-core stack or a disk file
	"""
	pass#IMPORTIMPORTIMPORT from utilities import model_blank, get_im
	if  type(data) == type(""): n = EMAN2_cppwrap.EMUtil.get_image_count(data)
	else:                       n = len(data)
	if listID == None:
		listID = list(range(n))
	img = sparx_utilities.get_im(data, 0)
	nx = img.get_xsize()
	ny = img.get_ysize()
	nz = img.get_zsize()
	if(mode == "a"):
		if(nz > 1):
			ali_params = "xform.align3d"
			pass#IMPORTIMPORTIMPORT from fundamentals import rot_shift3D
			pass#IMPORTIMPORTIMPORT from utilities import get_params3D
		else:
			ali_params = "xform.align2d"
			pass#IMPORTIMPORTIMPORT from fundamentals import rot_shift2D
			pass#IMPORTIMPORTIMPORT from utilities import get_params2D

	ave = sparx_utilities.model_blank(nx,ny,nz)
	var = sparx_utilities.model_blank(nx,ny,nz)
	nlistID = len(listID)
	for i in range(nlistID):
		img = sparx_utilities.get_im(data,listID[i])
		if(mode == "a"):
			if(nz > 1):
				phi, theta, psi, s3x, s3y, s3z, mirror, scale = sparx_utilities.get_params3D(img)
				img = sparx_fundamentals.rot_shift3D(img, phi, theta, psi, s3x, s3y, s3z, scale)
			else:
				angle, sx, sy, mirror, scale = sparx_utilities.get_params2D(img)
				img = sparx_fundamentals.rot_shift2D(img, angle, sx, sy, mirror, scale)
		EMAN2_cppwrap.Util.add_img(ave, img)
		EMAN2_cppwrap.Util.add_img2(var, img)
	EMAN2_cppwrap.Util.mul_scalar(ave, 1.0 /float(nlistID) )

	return ave, (var - ave*ave*nlistID)/(nlistID-1)

def ave_series(data, pave = True, mask = None):
	"""
		Calculate average of a image series using current alignment parameters
		data - real space image series
	"""
	pass#IMPORTIMPORTIMPORT from utilities    import model_blank, get_params2D
	pass#IMPORTIMPORTIMPORT from fundamentals import rot_shift2D
	n = len(data)
	nx = data[0].get_xsize()
	ny = data[0].get_ysize()
	ave = sparx_utilities.model_blank(nx, ny)
	for i in range(n):
		alpha, sx, sy, mirror, scale = sparx_utilities.get_params2D(data[i])
		temp = sparx_fundamentals.rot_shift2D(data[i], alpha, sx, sy, mirror)
		EMAN2_cppwrap.Util.add_img(ave, temp)
	if mask: EMAN2_cppwrap.Util.mul_img(ave, mask)
	if pave:  EMAN2_cppwrap.Util.mul_scalar(ave, 1.0/float(n))
	return ave

def ccc(img1, img2, mask=None):
	"""Cross-correlation coefficient.	   
	   Usage: result = ccc(image1, image2 [, mask])
	"""
	return img1.cmp("ccc", img2, {"mask":mask,"negative":0})

def fsc(img1, img2, w = 1.0, filename=None):
	"""Fourier Shell (or Ring) Correlation.

	   Usage: [frsc =] fsc(image1, image2 [, w, filename])

	   Computes the Fourier Shell (3d) or Fourier Ring (2d) correlation
	   function of two images.  If a filename is provided, then the 
	   result is saved using that filename.
	"""
	result = img1.calc_fourier_shell_correlation(img2, w)
	# repack results as a list of (freq, fsc, n) triplets
	size = len(result)/3
	frsc = []
	for i in range(3):
		frsc.append(result[i*size:(i+1)*size])
	if filename:
		outf = open(filename, "w")
		for i in range(size):
			datstrings = []
			datstrings.append("  %12f" % (frsc[0][i]))
			datstrings.append("  %12f" % (frsc[1][i]))
			datstrings.append("  %12f" % (frsc[2][i]))
			datstrings.append("\n")
			outf.write("".join(datstrings))
		outf.close()
	return frsc

def locres(vi, ui, m, nk, cutoff, step, myid, main_node, number_of_proc):
	pass#IMPORTIMPORTIMPORT from mpi 	  	  import mpi_init, mpi_comm_size, mpi_comm_rank, MPI_COMM_WORLD
	pass#IMPORTIMPORTIMPORT from mpi 	  	  import mpi_reduce, mpi_bcast, mpi_barrier, mpi_gatherv, mpi_send, mpi_recv
	pass#IMPORTIMPORTIMPORT from mpi 	  	  import MPI_SUM, MPI_FLOAT, MPI_INT
	pass#IMPORTIMPORTIMPORT from fundamentals import fft
	pass#IMPORTIMPORTIMPORT from utilities import model_blank, bcast_EMData_to_all, recv_EMData, send_EMData, bcast_number_to_all, info
	pass#IMPORTIMPORTIMPORT from filter import filt_tophatb
	pass#IMPORTIMPORTIMPORT from EMAN2 import rsconvolution
	pass#IMPORTIMPORTIMPORT from morphology import square_root, threshold
	

	nx = m.get_xsize()
	ny = m.get_ysize()
	nz = m.get_zsize()

	mc = sparx_utilities.model_blank(nx,ny,nz,1.0)-m

	if(myid == main_node):
		st = EMAN2_cppwrap.Util.infomask(vi,m,True)
		vi -= st[0]

		st = EMAN2_cppwrap.Util.infomask(ui,m,True)
		ui -= st[1]

	sparx_utilities.bcast_EMData_to_all(vi, myid, main_node)
	sparx_utilities.bcast_EMData_to_all(ui, myid, main_node)

	vf = sparx_fundamentals.fft(vi)
	uf = sparx_fundamentals.fft(ui)

	if(myid == 0):
		freqvol = sparx_utilities.model_blank(nx,ny,nz)
		resolut = []
	lp = int(max(nx,ny,nz)/2/step+0.5)
	step = 0.5/lp
	lt = lp//number_of_proc
	lp = (lt+1)*number_of_proc
	bailout = 0
	for i in range(myid,lp,number_of_proc):
		fl = step*i
		fh = fl+step
		freq=(fl+fh)/2.0
		#print number_of_proc,myid,lp,i,step,fl,fh,freq

		if i>0 :
			v = sparx_fundamentals.fft(sparx_filter.filt_tophatb( vf, fl, fh))
			u = sparx_fundamentals.fft(sparx_filter.filt_tophatb( uf, fl, fh))
			tmp1 = EMAN2_cppwrap.Util.muln_img(v,v)
			tmp2 = EMAN2_cppwrap.Util.muln_img(u,u)
			tmp3 = EMAN2_cppwrap.Util.muln_img(u,v)
			do = EMAN2_cppwrap.Util.infomask(sparx_morphology.square_root(sparx_morphology.threshold(EMAN2_cppwrap.Util.muln_img(tmp1,tmp2))),m,True)[0]
			dp = EMAN2_cppwrap.Util.infomask(tmp3,m,True)[0]
			#print "dpdo   ",myid,dp,do
			if do == 0.0: dis = [freq, 0.0]
			else:  dis = [freq, dp/do]
		else:
			tmp1 = sparx_utilities.model_blank(nx,ny,nz,1.0)
			tmp2 = sparx_utilities.model_blank(nx,ny,nz,1.0)
			tmp3 = sparx_utilities.model_blank(nx,ny,nz,1.0)
			dis = [freq, 1.0]


		tmp1 = EMAN2_cppwrap.Util.box_convolution(tmp1, nk)
		tmp2 = EMAN2_cppwrap.Util.box_convolution(tmp2, nk)
		tmp3 = EMAN2_cppwrap.Util.box_convolution(tmp3, nk)

		EMAN2_cppwrap.Util.mul_img(tmp1,tmp2)

		tmp1 = sparx_morphology.square_root(sparx_morphology.threshold(tmp1))

		EMAN2_cppwrap.Util.mul_img(tmp1,m)
		EMAN2_cppwrap.Util.add_img(tmp1,mc)


		EMAN2_cppwrap.Util.mul_img(tmp3,m)
		EMAN2_cppwrap.Util.add_img(tmp3,mc)

		EMAN2_cppwrap.Util.div_img(tmp3,tmp1)

		EMAN2_cppwrap.Util.mul_img(tmp3,m)

		mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

		if(myid == main_node):
			for k in range(number_of_proc):
				if(k != main_node):
					#print " start receiving",myid,i
					tag_node = k+1001
					dis = mpi.mpi_recv(2, mpi.MPI_FLOAT, k, sparx_global_def.SPARX_MPI_TAG_UNIVERSAL, mpi.MPI_COMM_WORLD)
					#print  "received ",myid,dis
					tmp3 = sparx_utilities.recv_EMData(k, tag_node)
					#print  "received ",myid
				if(dis[0] <=0.5):  resolut.append(dis)
				fl = step*(i+k)
				fh = fl+step
				freq=(fl+fh)/2.0
				#print k,dis,Util.infomask(tmp3,m,True)
				#if(k == number_of_proc-1):  bailout = 1
				bailout = 0
				#print  "setting freqvol  ",k
				EMAN2_cppwrap.Util.set_freq(freqvol,tmp3,m,cutoff, freq)
				"""Multiline Comment5"""
			"""Multiline Comment6"""

		else:
			tag_node = myid+1001
			#print   "sent from", myid,dis
			mpi.mpi_send(dis, 2, mpi.MPI_FLOAT, main_node, sparx_global_def.SPARX_MPI_TAG_UNIVERSAL, mpi.MPI_COMM_WORLD)
			#print   "sending EMD from", myid
			sparx_utilities.send_EMData(tmp3, main_node, tag_node)
			#print   "sent EMD from",myid

		bailout = sparx_utilities.bcast_number_to_all(bailout, main_node)
		if(bailout == 1):  break

	mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
	if( myid == main_node ):  return freqvol, resolut
	else:  return None, None


def histogram(image, mask = None, nbins = 0, hmin = 0.0, hmax = 0.0):
	"""
		Name
			histogram - calculate a histogram of the image pixel values
		Input
			input: input image
			mask: optional binary mask
			nbins:number of bins. Optional. 
			hmin, hmax: Optional. 
		Output
			h:a list containining 2*nbins elements.
	"""
	return  EMAN2_cppwrap.Util.histogram(image, mask, nbins, hmin, hmax)

class Munkres(object):
    """
Calculate the Munkres solution to the classical assignment problem.
See the module documentation for usage.
    """   
    def __init__(self):
        """Create a new instance"""
        self.C = None
        self.row_covered = []
        self.col_covered = []
        self.n = 0
        self.Z0_r = 0
        self.Z0_c = 0
        self.marked = None
        self.path = None

    def make_cost_matrix(profit_matrix, inversion_function):
        """
        Create a cost matrix from a profit matrix by calling
        'inversion_function' to invert each value. The inversion
        function must take one numeric argument (of any type) and return
        another numeric argument which is presumed to be the cost inverse
        of the original profit.

        This is a static method. Call it like this:

        cost_matrix = Munkres.make_cost_matrix(matrix, inversion_func)

        For example:

        cost_matrix = Munkres.make_cost_matrix(matrix, lambda x : sys.maxint - x)
        """
        cost_matrix = []
        for row in profit_matrix:
            cost_row = []
            for value in row:
                cost_row += [inversion_function(value)]
            cost_matrix += [cost_row]
        return cost_matrix

    make_cost_matrix = staticmethod(make_cost_matrix)

    def compute(self, cost_matrix):
        """
        Compute the indexes for the lowest-cost pairings between rows and
        columns in the database. Returns a list of (row, column) tuples
        that can be used to traverse the matrix.

        The matrix must be square.
        """
        self.C = self.__copy_matrix(cost_matrix)
        self.n = len(cost_matrix)
        self.row_covered = [False for i in range(self.n)]
        self.col_covered = [False for i in range(self.n)]
        self.Z0_r = 0
        self.Z0_c = 0
        self.path = self.__make_matrix(self.n * 2, 0)
        self.marked = self.__make_matrix(self.n, 0)

        done = False
        step = 1

        steps = { 1 : self.__step1,
                  2 : self.__step2,
                  3 : self.__step3,
                  4 : self.__step4,
                  5 : self.__step5,
                  6 : self.__step6 }

        while not done:
            try:
                func = steps[step]
                #print 'calling ' + str(func)
                step = func()
            except KeyError:
                done = True

        # Look for the starred columns
        results = []
        for i in range(self.n):
            for j in range(self.n):
                if self.marked[i][j] == 1:
                    results += [(i, j)]
        assert(len(results) == self.n)

        return results

    def __copy_matrix(self, matrix):
        """Return an exact copy of the supplied matrix"""
        copy = []
        for row in matrix:
            new_row = []
            for item in row:
                new_row += [item]
            copy += [new_row]
        return copy

    def __make_matrix(self, n, val):
        """Create an NxN matrix, populating it with the specific value."""
        matrix = []
        for i in range(n):
            matrix += [[val for j in range(n)]]
        return matrix

    def __step1(self):
        """
        For each row of the matrix, find the smallest element and
        subtract it from every element in its row. Go to Step 2.
        """
        C = self.C
        n = self.n
        for i in range(n):
            minval = self.C[i][0]
            # Find the minimum value for this row
            for j in range(n):
                if minval > self.C[i][j]:
                    minval = self.C[i][j]

            # Subtract that minimum from every element in the row.
            for j in range(n):
                self.C[i][j] -= minval

        return 2

    def __step2(self):
        """
        Find a zero (Z) in the resulting matrix. If there is no starred
        zero in its row or column, star Z. Repeat for each element in the
        matrix. Go to Step 3.
        """
        n = self.n
        for i in range(n):
            for j in range(n):
                if (self.C[i][j] == 0) and \
                   (not self.col_covered[j]) and \
                   (not self.row_covered[i]):
                    self.marked[i][j] = 1
                    self.col_covered[j] = True
                    self.row_covered[i] = True

        self.__clear_covers()
        return 3

    def __step3(self):
        """
        Cover each column containing a starred zero. If K columns are
        covered, the starred zeros describe a complete set of unique
        assignments. In this case, Go to DONE, otherwise, Go to Step 4.
        """
        n = self.n
        count = 0
        for i in range(n):
            for j in range(n):
                if self.marked[i][j] == 1:
                    self.col_covered[j] = True
                    count += 1

        if count >= n:
            step = 7 # done
        else:
            step = 4

        return step

    def __step4(self):
        """
        Find a noncovered zero and prime it. If there is no starred zero
        in the row containing this primed zero, Go to Step 5. Otherwise,
        cover this row and uncover the column containing the starred
        zero. Continue in this manner until there are no uncovered zeros
        left. Save the smallest uncovered value and Go to Step 6.
        """
        step = 0
        done = False
        row = -1
        col = -1
        star_col = -1
        while not done:
            (row, col) = self.__find_a_zero()
            if row < 0:
                done = True
                step = 6
            else:
                self.marked[row][col] = 2
                star_col = self.__find_star_in_row(row)
                if star_col >= 0:
                    col = star_col
                    self.row_covered[row] = True
                    self.col_covered[col] = False
                else:
                    done = True
                    self.Z0_r = row
                    self.Z0_c = col
                    step = 5

        return step

    def __step5(self):
        """
        Construct a series of alternating primed and starred zeros as
        follows. Let Z0 represent the uncovered primed zero found in Step 4.
        Let Z1 denote the starred zero in the column of Z0 (if any).
        Let Z2 denote the primed zero in the row of Z1 (there will always
        be one). Continue until the series terminates at a primed zero
        that has no starred zero in its column. Unstar each starred zero
        of the series, star each primed zero of the series, erase all
        primes and uncover every line in the matrix. Return to Step 3
        """
        count = 0
        path = self.path
        path[count][0] = self.Z0_r
        path[count][1] = self.Z0_c
        done = False
        while not done:
            row = self.__find_star_in_col(path[count][1])
            if row >= 0:
                count += 1
                path[count][0] = row
                path[count][1] = path[count-1][1]
            else:
                done = True

            if not done:
                col = self.__find_prime_in_row(path[count][0])
                count += 1
                path[count][0] = path[count-1][0]
                path[count][1] = col

        self.__convert_path(path, count)
        self.__clear_covers()
        self.__erase_primes()
        return 3

    def __step6(self):
        """
        Add the value found in Step 4 to every element of each covered
        row, and subtract it from every element of each uncovered column.
        Return to Step 4 without altering any stars, primes, or covered
        lines.
        """
        minval = self.__find_smallest()
        for i in range(self.n):
            for j in range(self.n):
                if self.row_covered[i]:
                    self.C[i][j] += minval
                if not self.col_covered[j]:
                    self.C[i][j] -= minval
        return 4

    def __find_smallest(self):
        """Find the smallest uncovered value in the matrix."""
        minval = sys.maxsize
        for i in range(self.n):
            for j in range(self.n):
                if (not self.row_covered[i]) and (not self.col_covered[j]):
                    if minval > self.C[i][j]:
                        minval = self.C[i][j]
        return minval

    def __find_a_zero(self):
        """Find the first uncovered element with value 0"""
        row = -1
        col = -1
        i = 0
        n = self.n
        done = False

        while not done:
            j = 0
            while True:
                if (self.C[i][j] == 0) and \
                   (not self.row_covered[i]) and \
                   (not self.col_covered[j]):
                    row = i
                    col = j
                    done = True
                j += 1
                if j >= n:
                    break
            i += 1
            if i >= n:
                done = True

        return (row, col)

    def __find_star_in_row(self, row):
        """
        Find the first starred element in the specified row. Returns
        the column index, or -1 if no starred element was found.
        """
        col = -1
        for j in range(self.n):
            if self.marked[row][j] == 1:
                col = j
                break

        return col

    def __find_star_in_col(self, col):
        """
        Find the first starred element in the specified row. Returns
        the row index, or -1 if no starred element was found.
        """
        row = -1
        for i in range(self.n):
            if self.marked[i][col] == 1:
                row = i
                break

        return row

    def __find_prime_in_row(self, row):
        """
        Find the first prime element in the specified row. Returns
        the column index, or -1 if no starred element was found.
        """
        col = -1
        for j in range(self.n):
            if self.marked[row][j] == 2:
                col = j
                break

        return col

    def __convert_path(self, path, count):
        for i in range(count+1):
            if self.marked[path[i][0]][path[i][1]] == 1:
                self.marked[path[i][0]][path[i][1]] = 0
            else:
                self.marked[path[i][0]][path[i][1]] = 1

    def __clear_covers(self):
        """Clear all covered matrix cells"""
        for i in range(self.n):
            self.row_covered[i] = False
            self.col_covered[i] = False

    def __erase_primes(self):
        """Erase all prime markings"""
        for i in range(self.n):
            for j in range(self.n):
                if self.marked[i][j] == 2:
                    self.marked[i][j] = 0

# NEEDS TO BE FIX
"""Multiline Comment11"""

# K-means main stability stream command line
def k_means_match_clusters_asg_new(asg1, asg2, T=0):
	# asg1 and asg2 are numpy array
	pass#IMPORTIMPORTIMPORT from numpy      import zeros, array
	pass#IMPORTIMPORTIMPORT from statistics import Munkres
	pass#IMPORTIMPORTIMPORT import sys

	K        = len(asg1)
	MAT      = [[0] * K for i in range(K)] 
	cost_MAT = [[0] * K for i in range(K)]
	dummy    = numpy.array([0], 'int32')
	for k1 in range(K):
		for k2 in range(K):
			MAT[k1][k2] = EMAN2_cppwrap.Util.k_means_cont_table(asg1[k1], asg2[k2], dummy, asg1[k1].size, asg2[k2].size, 0)
			if MAT[k1][k2] <= T:
				MAT[k1][k2] = 0
	for i in range(K):
		for j in range(K):
			cost_MAT[i][j] = sys.maxsize - MAT[i][j]
	m = Munkres()
	indexes = m.compute(cost_MAT)
	newindexes = []
	list_stable = []
	nb_tot_objs = 0
	for r, c in indexes:
		cont = MAT[r][c]
		if cont <= T:
			list_stable.append(numpy.array([], 'int32'))
			continue
		nb_tot_objs += cont
		objs = numpy.zeros(cont, 'int32')
		dummy = EMAN2_cppwrap.Util.k_means_cont_table(asg1[r], asg2[c], objs, asg1[r].size, asg2[c].size, 1)
		list_stable.append(objs)
		newindexes.append([r,c])

	return newindexes, list_stable, nb_tot_objs

# Hierarchical stability between partitions given by k-means
def hist_list(data, nbin = -1, fminiu = None, fmaxiu = None):
	"""
	  Calculate histogram of the list elements
	  nbin - number of bins, if not provided it will be set such that in average there is 10 elements per bin
	  fminiu - user provided minimum value for the histogram, it has to be smaller than the smallest element in data
	  fmaxiu - user provided maximum value for the histogram, it has to be larger than the largest element in data
	"""
	if nbin < 0:  nbin = len(data)/10
	fmaxi = max(data)
	fmini = min(data)

	if fmaxi == fmini:
		hist = [0]*nbin
		hist[0] = len(data)
		return [fmaxi]*nbin, hist
	if fminiu != None:
		if fminiu < fmini : fmini = fminiu
	if fmaxiu != None:
		if fmaxiu > fmaxi : fmaxi = fmaxiu

	binsize_i = (fmaxi-fmini)/float(nbin)
	start_i = fmini

	region = [None]*nbin
	hist = [None]*nbin
	for i in range(nbin):
		region[i] = start_i + i*binsize_i
		hist[i] = 0

	for d in data:
		i = min(int((d-start_i)/binsize_i), nbin-1)
		hist[i] += 1

	return region, hist

def pearson(X, Y):
	"""
	  Pearson correlation coefficient between two lists
	"""
	pass#IMPORTIMPORTIMPORT from math import sqrt
	Sx = Sy = Sxx = Syy = Sxy = 0.0
	N = len(X)
	for x, y in map(None, X, Y):
		Sx  += x
		Sy  += y
		Sxx += x*x
		Syy += y*y
		Sxy += x*y
	return (Sxy - Sx * Sy / N) / numpy.sqrt((Sxx - Sx*Sx/N)*(Syy - Sy*Sy/N))

def table_stat(X):
	"""
	  Basic statistics of numbers stored in a list: average, variance, minimum, maximum
	"""
	av = X[0]
	va = X[0]*X[0]
	mi = X[0]
	ma = X[0]
	N = len(X)
	for i in range(1,N):
		av += X[i]
		va += X[i]*X[i]
		mi = min(mi, X[i])
		ma = max(ma, X[i])
	return  av/N,(va - av*av/N)/float(N-1) , mi, ma

def mono(k1,k2):
	"""
	get index of a square nxn matrix stored in a triangular form
	for i in xrange(1,n):
	    for j in xrange(i):
		print  i,j,mono(i,j)

	"""
	mk = max(k1,k2)
	return  min(k1,k2) + mk*(mk-1)/2
	
def scale_fsc_datasetsize(fsc_to_be_adjusted, nfsc, nnew):
	s = float(nfsc)/float(nnew)
	fsc_sub = [0.0]*len(fsc_to_be_adjusted)
	for i,q in enumerate(fsc_to_be_adjusted):  fsc_sub[i] = q/(q*(1.0-s)+s)
	return fsc_sub
