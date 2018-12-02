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
import sparx_applications
import collections
import copy
import datetime
import sparx_filter
import fractions
import sparx_fundamentals
import sparx_global_def
import heapq
import inspect
import json
import math
import sparx_morphology
import mpi
import numpy
import operator
import os
import pickle
import random
import re
import sparx_reconstruction
import sets
import six
import socket
import sparx_statistics
import string
import struct
import sys
import time
import traceback
import types
import zlib
pass#IMPORTIMPORTIMPORT import EMAN2
pass#IMPORTIMPORTIMPORT import EMAN2_cppwrap
pass#IMPORTIMPORTIMPORT import EMAN2db
pass#IMPORTIMPORTIMPORT import applications
pass#IMPORTIMPORTIMPORT import collections
pass#IMPORTIMPORTIMPORT import copy
pass#IMPORTIMPORTIMPORT import datetime
pass#IMPORTIMPORTIMPORT import filter
pass#IMPORTIMPORTIMPORT import fractions
pass#IMPORTIMPORTIMPORT import functools
pass#IMPORTIMPORTIMPORT import fundamentals
pass#IMPORTIMPORTIMPORT import global_def
pass#IMPORTIMPORTIMPORT import heapq
pass#IMPORTIMPORTIMPORT import inspect
pass#IMPORTIMPORTIMPORT import json
pass#IMPORTIMPORTIMPORT import math
pass#IMPORTIMPORTIMPORT import morphology
pass#IMPORTIMPORTIMPORT import mpi
pass#IMPORTIMPORTIMPORT import numpy
pass#IMPORTIMPORTIMPORT import operator
pass#IMPORTIMPORTIMPORT import os
pass#IMPORTIMPORTIMPORT import pickle
pass#IMPORTIMPORTIMPORT import random
pass#IMPORTIMPORTIMPORT import re
pass#IMPORTIMPORTIMPORT import reconstruction
pass#IMPORTIMPORTIMPORT import sets
pass#IMPORTIMPORTIMPORT import six
pass#IMPORTIMPORTIMPORT import socket
pass#IMPORTIMPORTIMPORT import statistics
pass#IMPORTIMPORTIMPORT import string
pass#IMPORTIMPORTIMPORT import struct
pass#IMPORTIMPORTIMPORT import subprocess
pass#IMPORTIMPORTIMPORT import sys
pass#IMPORTIMPORTIMPORT import time
pass#IMPORTIMPORTIMPORT import traceback
pass#IMPORTIMPORTIMPORT import types
pass#IMPORTIMPORTIMPORT import utilities
pass#IMPORTIMPORTIMPORT import zlib
from future import standard_library
standard_library.install_aliases()
from builtins import range
from builtins import object
pass#IMPORTIMPORTIMPORT from global_def import *
pass#IMPORTIMPORTIMPORT from functools import reduce

def amoeba(var, scale, func, ftolerance=1.e-4, xtolerance=1.e-4, itmax=500, data=None):
	'''Use the simplex method to maximize a function of 1 or more variables.

	   Input:
		  var = the initial guess, a list with one element for each variable
		  scale = the search scale for each variable, a list with one
			  element for each variable.
		  func = the function to maximize.

	   Optional Input:
		  ftolerance = convergence criterion on the function values (default = 1.e-4)
		  xtolerance = convergence criterion on the variable values (default = 1.e-4)
		  itmax = maximum number of iterations allowed (default = 500).
		  data = data to be passed to func (default = None).

	   Output:
		  (varbest,funcvalue,iterations)
		  varbest = a list of the variables at the maximum.
		  funcvalue = the function value at the maximum.
		  iterations = the number of iterations used.

	   - Setting itmax to zero disables the itmax check and the routine will run
	     until convergence, even if it takes forever.
	   - Setting ftolerance or xtolerance to 0.0 turns that convergence criterion
	     off.  But do not set both ftolerance and xtolerance to zero or the routine
	     will exit immediately without finding the maximum.
	   - To check for convergence, check if (iterations < itmax).

	   The function should be defined like func(var,data) where
	   data is optional data to pass to the function.

	   Example:

	       pass#IMPORTIMPORTIMPORT import amoeba
	       def afunc(var,data=None): return 1.0-var[0]*var[0]-var[1]*var[1]
	       print amoeba.amoeba([0.25,0.25],[0.5,0.5],afunc)

	   Version 1.0 2005-March-28 T. Metcalf
		   1.1 2005-March-29 T. Metcalf - Use scale in simsize calculation.
						- Use func convergence *and* x convergence
						  rather than func convergence *or* x
						  convergence.
		   1.2 2005-April-03 T. Metcalf - When contracting, contract the whole
						  simplex.
	   '''

	nvar = len(var)       # number of variables in the minimization
	nsimplex = nvar + 1   # number of vertices in the simplex

	# first set up the simplex

	simplex = [0]*(nvar+1)  # set the initial simplex
	simplex[0] = var[:]
	for i in range(nvar):
		simplex[i+1] = var[:]
		simplex[i+1][i] += scale[i]

	fvalue = []
	for i in range(nsimplex):  # set the function values for the simplex
		fvalue.append(func(simplex[i],data=data))

	# Ooze the simplex to the maximum

	iteration = 0

	while 1:
		# find the index of the best and worst vertices in the simplex
		ssworst = 0
		ssbest  = 0
		for i in range(nsimplex):
			if fvalue[i] > fvalue[ssbest]:
				ssbest = i
			if fvalue[i] < fvalue[ssworst]:
				ssworst = i

		# get the average of the nsimplex-1 best vertices in the simplex
		pavg = [0.0]*nvar
		for i in range(nsimplex):
			if i != ssworst:
				for j in range(nvar): pavg[j] += simplex[i][j]
		for j in range(nvar): pavg[j] = pavg[j]/nvar # nvar is nsimplex-1
		simscale = 0.0
		for i in range(nvar):
			simscale += abs(pavg[i]-simplex[ssworst][i])/scale[i]
		simscale = simscale/nvar

		# find the range of the function values
		fscale = (abs(fvalue[ssbest])+abs(fvalue[ssworst]))/2.0
		if fscale != 0.0:
			frange = abs(fvalue[ssbest]-fvalue[ssworst])/fscale
		else:
			frange = 0.0  # all the fvalues are zero in this case

		# have we converged?
		if (((ftolerance <= 0.0 or frange < ftolerance) and    # converged to maximum
		 (xtolerance <= 0.0 or simscale < xtolerance)) or  # simplex contracted enough
		 (itmax and iteration >= itmax)):	     # ran out of iterations
			return simplex[ssbest],fvalue[ssbest],iteration

		# reflect the worst vertex
		pnew = [0.0]*nvar
		for i in range(nvar):
			pnew[i] = 2.0*pavg[i] - simplex[ssworst][i]
		fnew = func(pnew,data=data)
		if fnew <= fvalue[ssworst]:
			# the new vertex is worse than the worst so shrink
			# the simplex.
			for i in range(nsimplex):
				if i != ssbest and i != ssworst:
					for j in range(nvar):
						simplex[i][j] = 0.5*simplex[ssbest][j] + 0.5*simplex[i][j]
					fvalue[i] = func(simplex[i],data=data)
			for j in range(nvar):
				pnew[j] = 0.5*simplex[ssbest][j] + 0.5*simplex[ssworst][j]
			fnew = func(pnew,data=data)
		elif fnew >= fvalue[ssbest]:
			# the new vertex is better than the best so expand
			# the simplex.
			pnew2 = [0.0]*nvar
			for i in range(nvar):
				pnew2[i] = 3.0*pavg[i] - 2.0*simplex[ssworst][i]
			fnew2 = func(pnew2,data=data)
			if fnew2 > fnew:
				# accept the new vertex in the simplexe
				pnew = pnew2
				fnew = fnew2
		# replace the worst vertex with the new vertex
		for i in range(nvar):
			simplex[ssworst][i] = pnew[i]
		fvalue[ssworst] = fnew
		iteration += 1
		#print "Iteration:",iteration,"  ",ssbest,"  ",fvalue[ssbest]

def compose_transform2(alpha1, sx1, sy1, scale1, alpha2, sx2, sy2, scale2):
	"""Print the composition of two transformations  T2*T1
		Here  if v's are vectors:   vnew = T2*T1 vold
		     with T1 described by alpha1, sx1, scale1 etc.

	  Combined parameters correspond to image first transformed by set 1 followed by set 2.

	    Usage: compose_transform2(alpha1,sx1,sy1,scale1,alpha2,sx2,sy2,scale2)
	       angles in degrees
	"""

	t1 = EMAN2_cppwrap.Transform({"type":"2D","alpha":alpha1,"tx":sx1,"ty":sy1,"mirror":0,"scale":scale1})
	t2 = EMAN2_cppwrap.Transform({"type":"2D","alpha":alpha2,"tx":sx2,"ty":sy2,"mirror":0,"scale":scale2})
	tt = t2*t1
	d = tt.get_params("2D")
	return d[ "alpha" ], d[ "tx" ], d[ "ty" ], d[ "scale" ]


def compose_transform3(phi1,theta1,psi1,sx1,sy1,sz1,scale1,phi2,theta2,psi2,sx2,sy2,sz2,scale2):
	"""
	  Compute the composition of two transformations  T2*T1
		Here  if v's are vectors:	vnew = T2*T1 vold
		with T1 described by phi1, sx1,  scale1 etc.

		Usage: compose_transform3(phi1,theta1,psi1,sx1,sy1,sz1,scale1,phi2,theta2,psi2,sx2,sy2,sz2,scale2)
		   angles in degrees
	"""

	R1 = EMAN2_cppwrap.Transform({"type":"spider","phi":float(phi1),"theta":float(theta1),"psi":float(psi1),"tx":float(sx1),"ty":float(sy1),"tz":float(sz1),"mirror":0,"scale":float(scale1)})
	R2 = EMAN2_cppwrap.Transform({"type":"spider","phi":float(phi2),"theta":float(theta2),"psi":float(psi2),"tx":float(sx2),"ty":float(sy2),"tz":float(sz2),"mirror":0,"scale":float(scale2)})
	Rcomp=R2*R1
	d = Rcomp.get_params("spider")
	return d["phi"],d["theta"],d["psi"],d["tx"],d["ty"],d["tz"],d["scale"]

def combine_params2(alpha1, sx1, sy1, mirror1, alpha2, sx2, sy2, mirror2):
	"""
	  Combine 2D alignment parameters including mirror: Tnew = T2*T1
	  Combined parameters correspond to image first transformed by set 1 followed by set 2.
	"""

	t1 = EMAN2_cppwrap.Transform({"type":"2D","alpha":alpha1,"tx":sx1,"ty":sy1,"mirror":mirror1,"scale":1.0})
	t2 = EMAN2_cppwrap.Transform({"type":"2D","alpha":alpha2,"tx":sx2,"ty":sy2,"mirror":mirror2,"scale":1.0})
	tt = t2*t1
	"""Multiline Comment1"""
	d = tt.get_params("2D")
	return d[ "alpha" ], d[ "tx" ], d[ "ty" ], int(d[ "mirror" ]+0.1)

def inverse_transform2(alpha, tx = 0.0, ty = 0.0, mirror = 0):
	"""Returns the inverse of the 2d rot and trans matrix

	    Usage: nalpha, ntx, nty, mirror = inverse_transform2(alpha,tx,ty,mirror)
	"""

	t = EMAN2_cppwrap.Transform({"type":"2D","alpha":alpha,"tx":tx,"ty":ty,"mirror":mirror,"scale":1.0})
	t = t.inverse()
	t = t.get_params("2D")
	return t[ "alpha" ], t[ "tx" ], t[ "ty" ], int(t[ "mirror" ]+0.1)

def drop_image(imagename, destination, itype="h"):
	"""Write an image to the disk.

	Usage:  drop_image(name_of_existing_image, "path/to/image",
			  type = <type>)
	<type> is "h" (hdf) or "s" (spider)
	"""

	if type(destination) == type(""):
		if(itype == "h"):    imgtype = EMAN2_cppwrap.EMUtil.ImageType.IMAGE_HDF
		elif(itype == "s"):  imgtype = EMAN2_cppwrap.EMUtil.ImageType.IMAGE_SINGLE_SPIDER
		else:  sparx_global_def.ERROR("unknown image type","drop_image",1)
		imagename.write_image(destination, 0, imgtype)
	else:
		sparx_global_def.ERROR("destination is not a file name","drop_image",1)

def drop_spider_doc(filename, data, comment = None):
	"""Create a spider-compatible "Doc" file.

	   filename: name of the Doc file
	   data: List of lists, with the inner list being a list of floats
	         and is written as a line into the doc file.
	"""
	outf = open(filename, "w")
	pass#IMPORTIMPORTIMPORT from datetime import datetime
	outf.write(" ;   %s   %s   %s\n" % (datetime.datetime.now().ctime(), filename, comment))
	count = 1  # start key from 1; otherwise, it is confusing...
	for dat in data:
		try:
			nvals = len(dat)
			if nvals <= 5: datstrings = ["%5d %d" % (count, nvals)]
			else         : datstrings = ["%6d %d" % (count, nvals)]
			for num in dat:
				datstrings.append("%12.5g" % (num))
		except TypeError:
			# dat is a single number
			datstrings = ["%5d 1%12.5g" % (count, dat)]
		datstrings.append("\n")
		outf.write("".join(datstrings))
		count += 1
	outf.close()

def even_angles(delta = 15.0, theta1=0.0, theta2=90.0, phi1=0.0, phi2=359.99, \
				method = 'S', phiEqpsi = "Minus", symmetry='c1', ant = 0.0):
	"""Create a list of Euler angles suitable for projections.
	   method is either 'S' - for Saff algorithm
	               or   'P' - for Penczek '94 algorithm
			         'S' assumes phi1<phi2 and phi2-phi1>> delta ;
	   symmetry  - if this is set to point-group symmetry (cn or dn) or helical symmetry with point-group symmetry (scn or sdn), \
	                 it will yield angles from the asymmetric unit, not the specified range;
	   ant - neighborhood for local searches.  I believe the fastest way to deal with it in case of point-group symmetry
		        it to generate cushion projections within an/2 of the border of unique zone
	"""

	pass#IMPORTIMPORTIMPORT from math      import pi, sqrt, cos, acos, tan, sin
	pass#IMPORTIMPORTIMPORT from utilities import even_angles_cd
	pass#IMPORTIMPORTIMPORT from string    import lower,split
	angles = []
	symmetryLower = symmetry.lower()
	symmetry_string = symmetry.split()[0]
	if(symmetry_string[0]  == "c"):
		if(phi2 == 359.99):
			angles = even_angles_cd(delta, theta1, theta2, phi1-ant, phi2/int(symmetry_string[1:])+ant, method, phiEqpsi)
		else:
			angles = even_angles_cd(delta, theta1, theta2, phi1-ant, phi2+ant, method, phiEqpsi)
		if(int(symmetry_string[1:]) > 1):
			if( int(symmetry_string[1:])%2 ==0):
				qt = 360.0/int(symmetry_string[1:])
			else:
				qt = 180.0/int(symmetry_string[1:])
			n = len(angles)
			for i in range(n):
				t = n-i-1
				if(angles[t][1] == 90.0):
					if(angles[t][0] >= qt+ant):  del angles[t]
	elif(symmetry_string[0]  == "d"):
		if(phi2 == 359.99):
			angles = even_angles_cd(delta, theta1, theta2, phi1, 360.0/int(symmetry_string[1:]), method, phiEqpsi)
		else:
			angles = even_angles_cd(delta, theta1, theta2, phi1, phi2, method, phiEqpsi)
		n = len(angles)
		badb = 360.0/int(symmetry_string[1:])/4 + ant
		bade = 2*badb -ant
		bbdb = badb + 360.0/int(symmetry_string[1:])/2 + ant
		bbde = bbdb + 360.0/int(symmetry_string[1:])/4 - ant
		for i in range(n):
			t = n-i-1
			qt = angles[t][0]
			if((qt>=badb and qt<bade) or (qt>=bbdb and qt<bbde)):  del angles[t]

		if (int(symmetry_string[1:])%2 == 0):
			qt = 360.0/2/int(symmetry_string[1:])
		else:
			qt = 180.0/2/int(symmetry_string[1:])
		n = len(angles)
		for i in range(n):
			t = n-i-1
			if(angles[t][1] == 90.0):
				if(angles[t][0] >= qt + ant ):  del angles[t]
	elif(symmetry_string[0]  == "s"):

		#if symetry is "s", deltphi=delta, theata intial=theta1, theta end=90, delttheta=theta2
		# for helical, theta1 cannot be 0.0
		if theta1 > 90.0:
			sparx_global_def.ERROR('theta1 must be less than 90.0 for helical symmetry', 'even_angles', 1)
		if theta1 == 0.0: theta1 =90.0
		theta_number = int((90.0 - theta1)/theta2)
		#for helical, symmetry = s or scn
		cn = int(symmetry_string[2:])
		for j in range(theta_number,-1, -1):

			if( j == 0):
				if (symmetry_string[1] =="c"):
					if cn%2 == 0:
						k=int(359.99/cn/delta)
					else:
						k=int(359.99/2/cn/delta)
				elif (symmetry_string[1] =="d"):
					if cn%2 == 0:
						k=int(359.99/2/cn/delta)
					else:
						k=int(359.99/4/cn/delta)
				else:
					sparx_global_def.ERROR("For helical strucutre, we only support scn and sdn symmetry","even_angles",1)

			else:
				if (symmetry_string[1] =="c"):
					k=int(359.99/cn/delta)
				elif (symmetry_string[1] =="d"):
					k=int(359.99/2/cn/delta)

			for i in range(k+1):
					angles.append([i*delta,90.0-j*theta2,90.0])

	else:
		sparx_global_def.ERROR("even_angles","Symmetry not supported: "+symmetry_string,0)
	"""Multiline Comment2"""


	return angles

def even_angles_cd(delta, theta1=0.0, theta2=90.0, phi1=0.0, phi2=359.99, method = 'P', phiEQpsi='Minus'):
	"""Create a list of Euler angles suitable for projections.
	   method is either 'S' - for Saff algorithm
	                  or   'P' - for Penczek '94 algorithm
			  'S' assumes phi1<phi2 and phi2-phi1>> delta ;
	   phiEQpsi  - set this to 'Minus', if you want psi=-phi;
	"""
	pass#IMPORTIMPORTIMPORT from math import pi, sqrt, cos, acos
	angles = []
	if (method == 'P'):
		temp = EMAN2_cppwrap.Util.even_angles(delta, theta1, theta2, phi1, phi2)
		#		                                              phi, theta, psi
		for i in range(len(temp)/3): angles.append([temp[3*i],temp[3*i+1],temp[3*i+2]]);
	else:              #elif (method == 'S'):
		Deltaz  = numpy.cos(theta2*numpy.pi/180.0)-numpy.cos(theta1*numpy.pi/180.0)
		s       = delta*numpy.pi/180.0
		NFactor = 3.6/s
		wedgeFactor = abs(Deltaz*(phi2-phi1)/720.0)
		NumPoints   = int(NFactor*NFactor*wedgeFactor)
		angles.append([phi1, theta1, 0.0])
		z1 = numpy.cos(theta1*numpy.pi/180.0); 	phi=phi1            # initialize loop
		for k in range(1,(NumPoints-1)):
			z=z1 + Deltaz*k/(NumPoints-1)
			r= numpy.sqrt(1-z*z)
			phi = phi1+(phi + delta/r -phi1)%(abs(phi2-phi1))
			#[k, phi,180*acos(z)/pi, 0]
			angles.append([phi, 180*math.acos(z)/numpy.pi, 0.0])
		#angles.append([p2,t2,0])  # This is incorrect, as the last angle is really the border, not the element we need. PAP 01/15/07
	if (phiEQpsi == 'Minus'):
		for k in range(len(angles)): angles[k][2] = (720.0 - angles[k][0])%360.0
	if( theta2 == 180.0 or (theta2 >180. and delta == 180.0)):  angles.append( [0.0, 180.0, 0.0] )

	return angles

def find(vv, cmp_str, n):
	jFoundVec= [];
	for jFound in range(len(vv)):
		if (cmp_str=='lt'):
			if (vv[jFound]<n):
				jFoundVec.append(jFound);
		if (cmp_str=='le'):
			if (vv[jFound]<=n):
				jFoundVec.append(jFound);
		if (cmp_str=='eq'):
			if (vv[jFound]==n):
				jFoundVec.append(jFound);
		if (cmp_str=='ge'):
			if (vv[jFound]>=n):
				jFoundVec.append(jFound);
		if (cmp_str=='gt'):
			if (vv[jFound]>n):
				jFoundVec.append(jFound);
	return jFoundVec;

def gauss_edge(sharp_edge_image, kernel_size = 7, gauss_standard_dev =3):
	"""
		smooth sharp_edge_image with Gaussian function
		1. The sharp-edge image is convoluted with a gassian kernel
		2. The convolution normalized
	"""
	pass#IMPORTIMPORTIMPORT from utilities import model_gauss
	pass#IMPORTIMPORTIMPORT from EMAN2 import rsconvolution
	nz = sharp_edge_image.get_ndim()
	if(nz == 3):   kern = model_gauss(gauss_standard_dev, kernel_size , kernel_size, kernel_size)
	elif(nz == 2):  kern = model_gauss(gauss_standard_dev, kernel_size , kernel_size)
	else:          kern = model_gauss(gauss_standard_dev, kernel_size)
	aves = EMAN2_cppwrap.Util.infomask(kern, None, False)
	nx = kern.get_xsize()
	ny = kern.get_ysize()
	nz = kern.get_zsize()

	kern /= (aves[0]*nx*ny*nz)
	return  EMAN2_cppwrap.rsconvolution(sharp_edge_image, kern)

def get_image(imagename, nx = 0, ny = 1, nz = 1, im = 0):
	"""Read an image from the disk or assign existing object to the output.

	Usage: myimage = readImage("path/to/image")
	or     myimage = readImage(name_of_existing_image)
	"""
	if type(imagename) == type(""):
	    e = EMAN2_cppwrap.EMData()
	    e.read_image(imagename, im)
	elif not imagename:
	    e = EMAN2_cppwrap.EMData()
	    if (nx > 0): e.set_size(nx, ny, nz)
	else:
	    e = imagename
	return e

def get_im(stackname, im = 0):
	"""Read an image from the disk stack, or return im's image from the list of images

	Usage: myimage = get_im("path/to/stack", im)
	   or: myimage = get_im( data, im )
	"""
	if type(stackname) == type(""):
		e = EMAN2_cppwrap.EMData()
		e.read_image(stackname, im)
		return  e
	else:
		return  stackname[im].copy()

def get_image_data(img):
	"""
		Return a NumPy array containing the image data.
		Note: The NumPy array and the image data share the same memory,
		so if the NumPy array is altered then the image is altered
		as well (and vice versa).
	"""
	pass#IMPORTIMPORTIMPORT from EMAN2 import EMNumPy
	return EMAN2_cppwrap.EMNumPy.em2numpy(img)

def get_symt(symmetry):
	"""
	get a list of point-group symmetry transformations, symmetry="c3"
	"""
	pass#IMPORTIMPORTIMPORT from fundamentals import symclass
	scl = sparx_fundamentals.symclass(symmetry)
	trans = []
	for q in scl.symangles:
		trans.append(EMAN2_cppwrap.Transform({"type":"spider","phi":q[0],"theta":q[1],"psi":q[2]}))
	return trans

def get_input_from_string(str_input):
	"""
		Extract input numbers from a given string
	"""
	pass#IMPORTIMPORTIMPORT import re
	qq = re.split(" |,",str_input)
	for i in range(len(qq)-1, -1, -1):
		if(qq[i] == ""):  del qq[i]
	o = []
	for i in range(len(qq)):
		if(qq[i].find(".") >= 0):  o.append(float(qq[i]))
		else:  o.append(int(qq[i]))
	return o

def model_circle(r, nx, ny, nz=1):
	"""
	Create a centered circle (or sphere) having radius r.
	"""
	e = EMAN2_cppwrap.EMData()
	e.set_size(nx, ny, nz)
	e.process_inplace("testimage.circlesphere", {"radius":r, "fill":1})
	return e

def model_gauss(xsigma, nx, ny=1, nz=1, ysigma=None, zsigma=None, xcenter=None, ycenter=None, zcenter=None):
	"""
	Create an image of a Gaussian function with standard deviation "xsigma,ysigma,zsigma"
	 and centered at (xcenter,ycenter,zcenter), by default the center is image center.
	"""
	e = EMAN2_cppwrap.EMData()
	e.set_size(nx, ny, nz)
	if( ysigma  == None ) : ysigma = xsigma
	if( zsigma  == None ) : zsigma = xsigma
	if( xcenter == None ) : xcenter = nx//2
	if( ycenter == None ) : ycenter = ny//2
	if( zcenter == None ) : zcenter = nz//2
	e.process_inplace("testimage.puregaussian", {"x_sigma":xsigma,"y_sigma":ysigma,"z_sigma":zsigma,"x_center":xcenter,"y_center":ycenter,"z_center":zcenter} )
	return e

def model_gauss_noise(sigma, nx, ny=1, nz=1):
	"""
	Create an image of noise having standard deviation "sigma",
	and average 0.
	"""
	e = EMAN2_cppwrap.EMData()
	e.set_size(nx, ny, nz)
	if(sigma == 0.0): e.to_zero()
	else: e.process_inplace("testimage.noise.gauss", {"sigma":sigma})
	return e

def model_blank(nx, ny=1, nz=1, bckg = 0.0):
	"""
	Create a blank image.
	"""
	e = EMAN2_cppwrap.EMData()
	e.set_size(nx, ny, nz)
	e.to_zero()
	if( bckg != 0.0):  e+=bckg
	return e

def peak_search(e, npeak = 1, invert = 1, print_screen = 0):
	peaks    = e.peak_search(npeak, invert)
	ndim     = peaks[0]
	nlist    = int((len(peaks)-1)/((ndim+1)*2))
	if(nlist > 0):
		outpeaks = []
		if(print_screen):
			if  ndim == 1 :
				print('%10s%10s%10s%10s%10s'%("Index  "," Peak_value","X   ",		     "Peak/P_max", "X-NX/2"))
				print_list_format(peaks[1:], 4)
			elif ndim == 2 :
				print('%10s%10s%10s%10s%10s%10s%10s'%("Index  ", "Peak_value","X   ","Y   ",	     "Peak/P_max", "X-NX/2", "Y-NY/2"))
				print_list_format(peaks[1:], 6)
			elif ndim == 3 :
				print('%10s%10s%10s%10s%10s%10s%10s%10s%10s'%("Index  ", "Peak_value","X   ","Y   ","Z   ", "Peak/P_max", "X-NX/2", "Y-NY/2", "Z-NZ/2"))
				print_list_format(peaks[1:], 8)
			else:	sparx_global_def.ERROR("Image dimension extracted in peak_search is wrong", "Util.peak_search", 1)
		for i in range(nlist):
			k=int((ndim+1)*i*2)
			if   ndim == 1 :  p=[peaks[k+1], peaks[k+2], peaks[k+3], peaks[k+4]]
			elif ndim == 2 :  p=[peaks[k+1], peaks[k+2], peaks[k+3], peaks[k+4], peaks[k+5], peaks[k+6]]
			elif ndim == 3 :  p=[peaks[k+1], peaks[k+2], peaks[k+3], peaks[k+4], peaks[k+5], peaks[k+6], peaks[k+7], peaks[k+8]]
			outpeaks.append(p)
	else:
		ndim = e.get_ndim()
		#ERROR("peak search fails to find any peaks, returns image center as a default peak position","peak_search",0)
		if  ndim == 1 :
			nx = e.get_xsize()
			outpeaks = [[1.0, float(nx/2), 1.0, 0.0]]
		elif ndim == 2 :
			nx = e.get_xsize()
			ny = e.get_ysize()
			outpeaks = [[1.0, float(nx/2), float(ny/2), 1.0, 0.0, 0.0]]
		elif ndim == 3 :
			nx = e.get_xsize()
			ny = e.get_ysize()
			nz = e.get_ysize()
			outpeaks = [[1.0, float(nx/2), float(ny/2), float(nz/2), 1.0, 0.0, 0.0, 0.0]]
	return outpeaks

####--------------------------------------------------------------------------------------------------#########
def print_list_format(m, narray = 0):
	pass#IMPORTIMPORTIMPORT from string 	import split
	pass#IMPORTIMPORTIMPORT from math 	import sqrt
	pass#IMPORTIMPORTIMPORT import string
	pass#IMPORTIMPORTIMPORT import types
	"""
		Print formated elements in a list to screen
		The screen output is in the form of narray*int(len(m)/narray)
		Or when narray is zero, int(sqrt(len(m)))*int(sqrt(len(m)))
	"""
	flist = []
	for i in range(len(m)):
		if   type(m[i])  is float: flist.append('%10.3g'%(m[i]))
		elif type(m[i])  is int :  flist.append(  '%10d'%(m[i]))
		else				   : flist.append(  '%10s'%(m[i]))
	if(narray > len(m)):
		narray = 0
		sparx_global_def.ERROR("improper input narray number, use default value", "print_list_foramt",0)
	if(narray == 0 ):
		num = int(numpy.sqrt(len(m)))
		if( len(m) % num != 0): lnum = int(len(m)/num) + 1
		else: 			lnum = int(len(m)/num)
	else:
		num = narray
		if( len(m) % num == 0): lnum = int(len(m)/num)
		else: 			lnum = int(len(m)/num) + 1
	ncount = -1
	plist  = []
	for i in range(lnum):
		qlist = ""
		for j in range(num):
			ncount += 1
			if ncount <= len(m) - 1: qlist=qlist+flist[ncount]
			else:			 break
		plist.append(qlist)
	for i in range(lnum):
		print('%6d '%(i+1),plist[i])

def pad(image_to_be_padded, new_nx, new_ny = 1,	new_nz = 1, background = "average", off_center_nx = 0, off_center_ny = 0, off_center_nz = 0):
	pass#IMPORTIMPORTIMPORT import types
	if type(background) != bytes: background = str(background)
	if   background == "average"       :     image_padded = EMAN2_cppwrap.Util.pad(image_to_be_padded, new_nx, new_ny, new_nz, off_center_nx, off_center_ny, off_center_nz, "average")
	elif background == "circumference" :     image_padded = EMAN2_cppwrap.Util.pad(image_to_be_padded, new_nx, new_ny, new_nz, off_center_nx, off_center_ny, off_center_nz, "circumference")
	else:                                    image_padded = EMAN2_cppwrap.Util.pad(image_to_be_padded, new_nx, new_ny, new_nz, off_center_nx, off_center_ny, off_center_nz,  background  )
	return  image_padded

def chooseformat(t):
    ee = "%12.5f"%t
    if(len(ee)>12):  return "e"
    df1 = float(ee)
    df2 = float("%12.4e"%t)
    if(abs(t-df1) <= abs(t-df2)):  ouo = "f"
    else: ouo = "e"
    return ouo

def read_text_row(fnam, format="", skip=";"):
	"""
	 	Read a column-listed txt file.
		INPUT: filename: name of the Doc file
	 	OUTPUT:
	    	nc : number of entries in each lines (number of columns)
	    	len(data)/nc : number of lines (rows)
	    	data: List of numbers from the doc file
 	"""
	pass#IMPORTIMPORTIMPORT from string import split

	inf  = open(fnam, "r")
	strg = inf.readline()
	x    = []
	data = []
	while (len(strg) > 0):
		com_line = False
		for j in range(len(strg)):
			if(strg[j] == skip):	com_line = True
		if com_line == False:
			word=strg.split()
			if format == "s" :
				key = int(word[1])
				if key != len(word) - 2:
					del word
					word = []
					word.append(strg[0 : 5])
					word.append(strg[6 : 7])
					for k in range(key):
						k_start = 7       + k*13
						k_stop  = k_start + 13
						word.append(strg[k_start : k_stop])
			line=[]
			for i in range(len(word)):
				try:  line.append(int(word[i]))
				except:
					try:  	line.append(float(word[i]))
					except:	line.append(word[i])
			data.append(line)
		strg=inf.readline()
	inf.close
	return data


def write_text_row(data, file_name):
	"""
	   Write to an ASCII file a list of lists containing floats.

	   filename: name of the text file
	   data: List of lists, with the inner list being a list of floats, i.e., [ [first list], [second list], ...]
	         First list will be written as a first line, second as a second, and so on...
		 If only one list is given, the file will contain one line
	"""
	pass#IMPORTIMPORTIMPORT import types
	outf = open(file_name, "w")
	if (type(data[0]) == list):
		# It is a list of lists
		for i in range(len(data)):
			for j in range(len(data[i])):
				tpt = data[i][j]
				qtp = type(tpt)
				if qtp == int:		outf.write("  %12d"%tpt)
				elif qtp == float:
					frmt = chooseformat(tpt)
					if( frmt == "f" ):			outf.write("  %12.5f"%tpt)
					else:						outf.write("  %12.5e"%tpt)
				else:                   		outf.write("  %s"%tpt)
			outf.write("\n")
	else:
		# Single list
		for j in range(len(data)):
			tpt = data[j]
			qtp = type(tpt)
			if qtp == int :			outf.write("  %12d\n"%tpt)
			elif qtp == float:
				frmt = chooseformat(tpt)
				if( frmt == "f" ):				outf.write("  %12.5f\n"%tpt)
				else:							outf.write("  %12.5e\n"%tpt)
			else:								outf.write("  %s\n"%tpt)
	outf.flush()
	outf.close()


def read_text_file(file_name, ncol = 0):
	"""
		Read data from text file, if ncol = -1, read all columns
		if ncol >= 0, just read the (ncol)-th column.
	"""

	pass#IMPORTIMPORTIMPORT from string import split
	inf = open(file_name, "r")
	line = inf.readline()
	data = []
	while len(line) > 0:
		if ncol == -1:
			vdata = line.split()
			if data == []:
				for i in range(len(vdata)):
					try:     data.append([int(vdata[i])])
					except:
						try:  	 data.append([float(vdata[i])])
						except:  data.append([vdata[i]])
			else:
				for i in range(len(vdata)):
					try:  data[i].append(int(vdata[i]))
					except:
						try:  data[i].append(float(vdata[i]))
						except:  data[i].append(vdata[i])
		else:
			vdata = line.split()[ncol]
			try:     data.append(int(vdata))
			except:
				try:  	data.append(float(vdata))
				except:  data.append(vdata)
		line = inf.readline()
	return data

def write_text_file(data, file_name):
	"""
	   Write to an ASCII file a list of lists containing floats.

	   filename: name of the text file
	   data: List of lists, with the inner list being a list of floats, i.e., [ [first list], [second list], ...]
	         First list will be written as a first column, second as a second, and so on...
		 If only one list is given, the file will contain one column
	"""

	if data == []:
		outf = open(file_name, "w")
		outf.close()
		return

	pass#IMPORTIMPORTIMPORT import types
	outf = open(file_name, "w")
	if (type(data[0]) == list):
		# It is a list of lists
		for i in range(len(data[0])):
			for j in range(len(data)):
				tpt = data[j][i]
				qtp = type(tpt)
				if qtp == int:			outf.write("  %12d"%tpt)
				elif qtp == float:
					frmt = chooseformat(tpt)
					if( frmt == "f" ):				outf.write("  %12.5f"%tpt)
					else:							outf.write("  %12.5e"%tpt)
				else:                   			outf.write("  %s"%tpt)
			outf.write("\n")
	else:
		# Single list
		for j in range(len(data)):
			tpt = data[j]
			qtp = type(tpt)
			if qtp == int :			outf.write("  %12d\n"%tpt)
			elif qtp == float:
				frmt = chooseformat(tpt)
				if( frmt == "f" ):				outf.write("  %12.5f\n"%tpt)
				else:							outf.write("  %12.5e\n"%tpt)
			else:                   			outf.write("  %s\n"%tpt)
	outf.close()

def rotate_shift_params(paramsin, transf):
	# moved from sxprocess.py
	if len(paramsin[0])>3 :
		pass#IMPORTIMPORTIMPORT from EMAN2 import Vec2f
		t = EMAN2_cppwrap.Transform({"type":"spider","phi":transf[0],"theta":transf[1],"psi":transf[2],"tx":transf[3],"ty":transf[4],"tz":transf[5],"mirror":0,"scale":1.0})
		t = t.inverse()
		cpar = []
		for params in paramsin:
			d = EMAN2_cppwrap.Transform({"type":"spider","phi":params[0],"theta":params[1],"psi":params[2]})
			d.set_trans(EMAN2_cppwrap.Vec2f(-params[3], -params[4]))
			c = d*t
			u = c.get_params("spider")
			cpar.append([u["phi"],u["theta"],u["psi"],-u["tx"],-u["ty"]])
	else:
		t = EMAN2_cppwrap.Transform({"type":"spider","phi":transf[0],"theta":transf[1],"psi":transf[2]})
		t = t.inverse()
		cpar = []
		for params in paramsin:
			d = EMAN2_cppwrap.Transform({"type":"spider","phi":params[0],"theta":params[1],"psi":params[2]})
			c = d*t
			u = c.get_params("spider")
			cpar.append([u["phi"],u["theta"],u["psi"]])
			#cpar.append([u["phi"],u["theta"],u["psi"],-u["tx"],-u["ty"]])
	return cpar

def reshape_1d(input_object, length_current=0, length_interpolated=0, Pixel_size_current = 0., Pixel_size_interpolated = 0.):
	"""
		linearly interpolate a 1D power spectrum to required length with required Pixel size
		input_object - a 1D list with a 1D curve to be interpolated
		length_current - half size of the image size (in case of power spectrum, it can be different from the length of the input_object)
		length_interpolated - length of the interpolated 1D curve
		Pixel_size_current - pixel size of the input 1D list
		Pixel_size_interpolated - pixel size of the target 1D list
		One can either input the two lengths or two respective pixel sizes
	"""

	interpolated = []
	if length_current == 0: length_current = len(input_object)
	lt = len(input_object) - 2
	if length_interpolated == 0:
		if( Pixel_size_interpolated != Pixel_size_current):
			length_interpolated = int(length_current*Pixel_size_current/Pixel_size_interpolated + 0.5)
		else:
			sparx_global_def.ERROR("Incorrect input parameters","reshape_1d",1)
			return []

	if  Pixel_size_current == 0.:
		Pixel_size_current = 1.
		Pixel_size_interpolated = Pixel_size_current*float(length_current)/float(length_interpolated)
	qt =Pixel_size_interpolated/Pixel_size_current

	for i in range(length_interpolated):
		xi = float(i)*qt
		ix = min(int(xi),lt)
		df = xi -ix
		xval = input_object[ix] + df*(input_object[ix+1]-input_object[ix])
		interpolated.append(xval)
	return interpolated

def estimate_3D_center_MPI(data, nima, myid, number_of_proc, main_node, mpi_comm=None):
	pass#IMPORTIMPORTIMPORT from math import cos, sin, radians
	pass#IMPORTIMPORTIMPORT from numpy import matrix
	pass#IMPORTIMPORTIMPORT from numpy import linalg
	pass#IMPORTIMPORTIMPORT from mpi import MPI_COMM_WORLD
	pass#IMPORTIMPORTIMPORT from mpi import mpi_recv, mpi_send, MPI_FLOAT
	pass#IMPORTIMPORTIMPORT from applications import MPI_start_end

	if mpi_comm == None:
		mpi_comm = mpi.MPI_COMM_WORLD

	ali_params_series = []
	for im in data:
		phi, theta, psi, s2x, s2y = get_params_proj(im)
		ali_params_series.append(phi)
		ali_params_series.append(theta)
		ali_params_series.append(psi)
		ali_params_series.append(s2x)
		ali_params_series.append(s2y)

	if myid == main_node:
		for proc in range(number_of_proc):
			if proc != main_node:
				image_start_proc, image_end_proc = sparx_applications.MPI_start_end(nima, number_of_proc, proc)
				n_params = (image_end_proc - image_start_proc)*5
				temp = mpi.mpi_recv(n_params, mpi.MPI_FLOAT, proc, proc, mpi_comm)
				for nn in range(n_params):
					ali_params_series.append(float(temp[nn]))

		ali_params = []
		N = len(ali_params_series)/5
		for im in range(N):
			ali_params.append([ali_params_series[im*5], ali_params_series[im*5+1], ali_params_series[im*5+2], ali_params_series[im*5+3], ali_params_series[im*5+4]])

		A = []
		b = []

		for i in range(N):
			phi_rad   = numpy.radians(ali_params[i][0])
			theta_rad = numpy.radians(ali_params[i][1])
			psi_rad   = numpy.radians(ali_params[i][2])
			A.append([numpy.cos(psi_rad)*numpy.cos(theta_rad)*numpy.cos(phi_rad)-numpy.sin(psi_rad)*numpy.sin(phi_rad),
				numpy.cos(psi_rad)*numpy.cos(theta_rad)*numpy.sin(phi_rad)+numpy.sin(psi_rad)*numpy.cos(phi_rad), -numpy.cos(psi_rad)*numpy.sin(theta_rad), 1, 0])
			A.append([-numpy.sin(psi_rad)*numpy.cos(theta_rad)*numpy.cos(phi_rad)-numpy.cos(psi_rad)*numpy.sin(phi_rad),
				-numpy.sin(psi_rad)*numpy.cos(theta_rad)*numpy.sin(phi_rad)+numpy.cos(psi_rad)*numpy.cos(phi_rad), numpy.sin(psi_rad)*numpy.sin(theta_rad), 0, 1])
			b.append([ali_params[i][3]])
			b.append([ali_params[i][4]])

		A_matrix = numpy.matrix(A)
		b_matrix = numpy.matrix(b)

		K = numpy.linalg.solve(A_matrix.T*A_matrix, A_matrix.T*b_matrix)
		return float(K[0][0]), float(K[1][0]), float(K[2][0]), float(K[3][0]), float(K[4][0])

	else:
		image_start_proc, image_end_proc = sparx_applications.MPI_start_end(nima, number_of_proc, myid)
		n_params = (image_end_proc - image_start_proc)*5
		mpi.mpi_send(ali_params_series, n_params, mpi.MPI_FLOAT, main_node, myid, mpi_comm)

		return 0.0, 0.0, 0.0, 0.0, 0.0


def rotate_3D_shift(data, shift3d):

	t = EMAN2_cppwrap.Transform({"type":"spider","phi":0.0,"theta":0.0,"psi":0.0,"tx":-shift3d[0],"ty":-shift3d[1],"tz":-shift3d[2],"mirror":0,"scale":1.0})

	for i in range(len(data)):
		d = data[i].get_attr('xform.projection')
		c = d*t
		data[i].set_attr('xform.projection', c)


def set_arb_params(img, params, par_str):

	"""
		filling arbitary headers
	"""
	for i in range(len(par_str)): img.set_attr_dict({par_str[i]:params[i]})

def get_arb_params(img, par_str):

	"""
		reading arbitary headers
	"""
	params=[]
	for i in range(len(par_str)): params.append(img.get_attr(par_str[i]))
	return params

###------------------------------------------------------------------------------------------

def reduce_EMData_to_root(data, myid, main_node = 0, comm = -1):
	pass#IMPORTIMPORTIMPORT from numpy import shape, reshape
	pass#IMPORTIMPORTIMPORT from mpi   import mpi_reduce, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD, mpi_barrier

	if comm == -1 or comm == None:  comm = mpi.MPI_COMM_WORLD

	array = get_image_data(data)
	n     = numpy.shape(array)
	ntot  = 1
	for i in n: ntot *= i
	count = (75*4+2)*(75*4)**2
	array1d = numpy.reshape( array, (ntot,))
	ntime = (ntot-1) /count + 1
	for i in range(ntime):
		block_begin = i*count
		block_end   = min(block_begin + count, ntot)
		block_size  = block_end - block_begin
		tmpsum = mpi.mpi_reduce(array1d[block_begin:block_begin+block_size], block_size, mpi.MPI_FLOAT, mpi.MPI_SUM, main_node, comm)
		mpi.mpi_barrier(comm)
		if myid == main_node:
			array1d[block_begin:block_end] = tmpsum[0:block_size]

def bcast_compacted_EMData_all_to_all(list_of_em_objects, myid, comm=-1):

	"""
	The assumption in <<bcast_compacted_EMData_all_to_all>> is that each processor
	calculates part of the list of elements and then each processor sends
	its results to the other ones.

	Therefore, each processor has access to the header. If we assume that the
	attributes of interest from the header are the same for all elements then
	we can copy the header and no mpi message is necessary for the
	header.

	"""
	pass#IMPORTIMPORTIMPORT from applications import MPI_start_end
	pass#IMPORTIMPORTIMPORT from EMAN2 import EMNumPy
	pass#IMPORTIMPORTIMPORT from numpy import concatenate, shape, array, split
	pass#IMPORTIMPORTIMPORT from mpi import mpi_comm_size, mpi_bcast, MPI_FLOAT, MPI_COMM_WORLD
	pass#IMPORTIMPORTIMPORT from numpy import reshape

	if comm == -1 or comm == None: comm = mpi.MPI_COMM_WORLD

	num_ref = len(list_of_em_objects)
	ncpu = mpi.mpi_comm_size(comm)	# Total number of processes, passed by --np option.

	ref_start, ref_end = sparx_applications.MPI_start_end(num_ref, ncpu, myid)

	for first_myid_process_that_has_em_elements in range(ncpu):
		sim_start, sim_ref_end = sparx_applications.MPI_start_end(num_ref, ncpu, first_myid_process_that_has_em_elements)
		if sim_start != sim_ref_end:
			break
	else:
		raise ValueError("No processor contains em objects!")

	if myid == first_myid_process_that_has_em_elements:
		# used for copying the header and other info

		reference_em_object = list_of_em_objects[ref_start]
		data = EMAN2_cppwrap.EMNumPy.em2numpy(reference_em_object)
		size_of_one_refring_assumed_common_to_all = data.size

		nx = reference_em_object.get_xsize()
		ny = reference_em_object.get_ysize()
		nz = reference_em_object.get_zsize()

		em_dict = reference_em_object.get_attr_dict()
		dict_to_send = {"size_of_one_refring_assumed_common_to_all":size_of_one_refring_assumed_common_to_all, \
						"em_dict":em_dict, "nx":nx, "ny":ny, "nz":nz}
	else:
		dict_to_send = None

	dict_received = wrap_mpi_bcast(dict_to_send, first_myid_process_that_has_em_elements, comm)

	em_dict = dict_received["em_dict"]
	nx = dict_received["nx"]
	ny = dict_received["ny"]
	nz = dict_received["nz"]
	size_of_one_refring_assumed_common_to_all = dict_received["size_of_one_refring_assumed_common_to_all"]

	if size_of_one_refring_assumed_common_to_all*(ref_end-ref_start) > (2**31-1):
		print("Sending refrings: size of data to broadcast is greater than 2GB")

	for sender_id in range(ncpu):
		sender_ref_start, sender_ref_end = sparx_applications.MPI_start_end(num_ref, ncpu, sender_id)

		if sender_id == myid:
			if ref_start == ref_end:
				continue
			data = EMAN2_cppwrap.EMNumPy.em2numpy(list_of_em_objects[ref_start])  #array([], dtype = 'float32')
			for i in range(ref_start+1,ref_end):
				data = numpy.concatenate([data, EMAN2_cppwrap.EMNumPy.em2numpy(list_of_em_objects[i])])
		else:
			if sender_ref_start == sender_ref_end:
				continue
			data = numpy.array([], dtype = 'float32')

		sender_size_of_refrings = (sender_ref_end - sender_ref_start)*size_of_one_refring_assumed_common_to_all

		# size_of_refrings = mpi_bcast(size_of_refrings, 1, MPI_INT, sender_id, comm)
		data = mpi.mpi_bcast(data, sender_size_of_refrings, mpi.MPI_FLOAT, sender_id, comm)
		# print "Just sent %d float32 elements"%data.size

		if myid != sender_id:
			for i in range(sender_ref_start, sender_ref_end):
				offset_ring = sender_ref_start
				start_p = (i-offset_ring)*size_of_one_refring_assumed_common_to_all
				end_p   = (i+1-offset_ring)*size_of_one_refring_assumed_common_to_all
				image_data = data[start_p:end_p]

				if int(nz) != 1:
					image_data = numpy.reshape(image_data, (nz, ny, nx))
				elif ny != 1:
					image_data = numpy.reshape(image_data, (ny, nx))

				em_object = EMAN2_cppwrap.EMNumPy.numpy2em(image_data)
				em_object.set_attr_dict(em_dict)
				list_of_em_objects[i] = em_object

def gather_compacted_EMData_to_root(number_of_all_em_objects_distributed_across_processes, list_of_em_objects_for_myid_process, myid, comm=-1):

	"""

	The assumption in <<gather_compacted_EMData_to_root>> is that each processor
	calculates part of the list of elements and then each processor sends
	its results to the root

	Therefore, each processor has access to the header. If we assume that the
	attributes of interest from the header are the same for all elements then
	we can copy the header and no mpi message is necessary for the
	header.

	"""
	pass#IMPORTIMPORTIMPORT from applications import MPI_start_end
	pass#IMPORTIMPORTIMPORT from EMAN2 import EMNumPy
	pass#IMPORTIMPORTIMPORT from numpy import concatenate, shape, array, split
	pass#IMPORTIMPORTIMPORT from mpi import mpi_comm_size, mpi_bcast, MPI_FLOAT, MPI_COMM_WORLD
	pass#IMPORTIMPORTIMPORT from numpy import reshape
	pass#IMPORTIMPORTIMPORT from mpi import mpi_recv, mpi_send, mpi_barrier

	if comm == -1 or comm == None: comm = mpi.MPI_COMM_WORLD

	ncpu = mpi.mpi_comm_size(comm)	# Total number of processes, passed by --np option.

	ref_start, ref_end = sparx_applications.MPI_start_end(number_of_all_em_objects_distributed_across_processes, ncpu, myid)
	ref_end -= ref_start
	ref_start = 0
	tag_for_send_receive = 123456

	# used for copying the header
	reference_em_object = list_of_em_objects_for_myid_process[ref_start]
	nx = reference_em_object.get_xsize()
	ny = reference_em_object.get_ysize()
	nz = reference_em_object.get_zsize()
	# is_complex = reference_em_object.is_complex()
	is_ri = reference_em_object.is_ri()
	changecount = reference_em_object.get_attr("changecount")
	is_complex_x = reference_em_object.is_complex_x()
	is_complex_ri = reference_em_object.get_attr("is_complex_ri")
	apix_x = reference_em_object.get_attr("apix_x")
	apix_y = reference_em_object.get_attr("apix_y")
	apix_z = reference_em_object.get_attr("apix_z")
	is_complex = reference_em_object.get_attr_default("is_complex",1)
	is_fftpad = reference_em_object.get_attr_default("is_fftpad",1)
	is_fftodd = reference_em_object.get_attr_default("is_fftodd", nz%2)

	data = EMAN2_cppwrap.EMNumPy.em2numpy(list_of_em_objects_for_myid_process[ref_start])
	size_of_one_refring_assumed_common_to_all = data.size

	# n = shape(data)
	# size_of_one_refring_assumed_common_to_all = 1
	# for i in n: size_of_one_refring_assumed_common_to_all *= i

	if size_of_one_refring_assumed_common_to_all*(ref_end-ref_start) > (2**31-1):
		print("Sending refrings: size of data to broadcast is greater than 2GB")

	for sender_id in range(1,ncpu):
		if sender_id == myid:
			data = EMAN2_cppwrap.EMNumPy.em2numpy(list_of_em_objects_for_myid_process[ref_start])  #array([], dtype = 'float32')
			for i in range(ref_start+1,ref_end):
				data = numpy.concatenate([data, EMAN2_cppwrap.EMNumPy.em2numpy(list_of_em_objects_for_myid_process[i])])
		else:
			data = numpy.array([], dtype = 'float32')

		sender_ref_start, sender_ref_end = sparx_applications.MPI_start_end(number_of_all_em_objects_distributed_across_processes, ncpu, sender_id)

		sender_size_of_refrings = (sender_ref_end - sender_ref_start)*size_of_one_refring_assumed_common_to_all

		if myid == 0:
			# print "root, receiving from ", sender_id, "  sender_size_of_refrings = ", sender_size_of_refrings
			data = mpi.mpi_recv(sender_size_of_refrings,mpi.MPI_FLOAT, sender_id, tag_for_send_receive, mpi.MPI_COMM_WORLD)
		elif sender_id == myid:
			# print "sender_id = ", sender_id, "sender_size_of_refrings = ", sender_size_of_refrings
			mpi.mpi_send(data, sender_size_of_refrings, mpi.MPI_FLOAT, 0, tag_for_send_receive, mpi.MPI_COMM_WORLD)

		mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

		# if myid != sender_id:
		if myid == 0:
			for i in range(sender_ref_start, sender_ref_end):
				offset_ring = sender_ref_start
				start_p = (i-offset_ring)*size_of_one_refring_assumed_common_to_all
				end_p   = (i+1-offset_ring)*size_of_one_refring_assumed_common_to_all
				image_data = data[start_p:end_p]

				if int(nz) != 1:
					image_data = numpy.reshape(image_data, (nz, ny, nx))
				elif ny != 1:
					image_data = numpy.reshape(image_data, (ny, nx))

				em_object = EMAN2_cppwrap.EMNumPy.numpy2em(image_data)

				# em_object.set_complex(is_complex)
				em_object.set_ri(is_ri)
				em_object.set_attr_dict({
				"changecount":changecount,
				"is_complex_x":is_complex_x,
				"is_complex_ri":is_complex_ri,
				"apix_x":apix_x,
				"apix_y":apix_y,
				"apix_z":apix_z,
				'is_complex':is_complex,
				'is_fftodd':is_fftodd,
				'is_fftpad':is_fftpad})

				# list_of_em_objects[i] = em_object
				list_of_em_objects_for_myid_process.append(em_object)

		mpi.mpi_barrier(mpi.MPI_COMM_WORLD)


def bcast_EMData_to_all(tavg, myid, source_node = 0, comm = -1):
	pass#IMPORTIMPORTIMPORT from EMAN2 import EMNumPy
	pass#IMPORTIMPORTIMPORT from numpy import array, shape, reshape
	pass#IMPORTIMPORTIMPORT from mpi   import mpi_bcast, MPI_FLOAT, MPI_COMM_WORLD

	if comm == -1 or comm == None: comm = mpi.MPI_COMM_WORLD
	tavg_data = EMAN2_cppwrap.EMNumPy.em2numpy(tavg)
	n = numpy.shape(tavg_data)
	ntot = 1
	for i in n: ntot *= i
	tavg_tmp = mpi.mpi_bcast(tavg_data, ntot, mpi.MPI_FLOAT, source_node, comm)
	if(myid != source_node):
		tavg_data1d = numpy.reshape(tavg_data,(ntot,))
		tavg_data1d[0:ntot] = tavg_tmp[0:ntot]

"""Multiline Comment4"""

def send_EMData(img, dst, tag, comm=-1):
	pass#IMPORTIMPORTIMPORT from mpi import mpi_send, MPI_INT, MPI_FLOAT, MPI_COMM_WORLD

	if comm == -1: comm = mpi.MPI_COMM_WORLD
	img_head = []
	img_head.append(img.get_xsize())
	img_head.append(img.get_ysize())
	img_head.append(img.get_zsize())
	img_head.append(img.is_complex())
	img_head.append(img.is_ri())
	img_head.append(img.get_attr("changecount"))
	img_head.append(img.is_complex_x())
	img_head.append(img.get_attr("is_complex_ri"))
	img_head.append(int(img.get_attr("apix_x")*10000))
	img_head.append(int(img.get_attr("apix_y")*10000))
	img_head.append(int(img.get_attr("apix_z")*10000))

	head_tag = 2*tag
	mpi.mpi_send(img_head, 11, mpi.MPI_INT, dst, head_tag, comm)

	img_data = get_image_data(img)
	data_tag = 2*tag+1
	ntot = img_head[0]*img_head[1]*img_head[2]
	mpi.mpi_send(img_data, ntot, mpi.MPI_FLOAT, dst, data_tag, comm)

	"""Multiline Comment5"""

def recv_EMData(src, tag, comm=-1):
	pass#IMPORTIMPORTIMPORT from mpi import mpi_recv, MPI_INT, MPI_FLOAT, MPI_COMM_WORLD
	pass#IMPORTIMPORTIMPORT from numpy import reshape
	pass#IMPORTIMPORTIMPORT from EMAN2 import EMNumPy


	if comm==-1: comm = mpi.MPI_COMM_WORLD
	head_tag = 2*tag
	img_head = mpi.mpi_recv(11, mpi.MPI_INT, src, head_tag, comm)

	nx = int(img_head[0])
	ny = int(img_head[1])
	nz = int(img_head[2])
	is_complex = int(img_head[3])
	is_ri = int(img_head[4])

	data_tag = 2*tag+1
	ntot = nx*ny*nz

	img_data = mpi.mpi_recv(ntot, mpi.MPI_FLOAT, src, data_tag, comm)
	if nz != 1:
		img_data = numpy.reshape(img_data, (nz, ny, nx))
	elif ny != 1:
		img_data = numpy.reshape(img_data, (ny, nx))
	else:
		pass

	img = EMAN2_cppwrap.EMNumPy.numpy2em(img_data)
	img.set_complex(is_complex)
	img.set_ri(is_ri)
	img.set_attr_dict({"changecount":int(img_head[5]),  "is_complex_x":int(img_head[6]),  "is_complex_ri":int(img_head[7]),  "apix_x":int(img_head[8])/10000.0,  "apix_y":int(img_head[9])/10000.0,  "apix_z":int(img_head[10])/10000.0})
	return img

	"""Multiline Comment6"""

def bcast_number_to_all(number_to_send, source_node = 0, mpi_comm = -1):
	"""
		number_to_send has to be pre-defined in each node
	"""
	pass#IMPORTIMPORTIMPORT from mpi import mpi_bcast, MPI_INT, MPI_COMM_WORLD, MPI_FLOAT
	pass#IMPORTIMPORTIMPORT import types
	if mpi_comm == -1:  mpi_comm = mpi.MPI_COMM_WORLD
	if    type(number_to_send) is int:
		TMP = mpi.mpi_bcast(number_to_send, 1, mpi.MPI_INT,   source_node, mpi_comm)
		return int(TMP[0])
	elif  type(number_to_send) is float:
		TMP = mpi.mpi_bcast(number_to_send, 1, mpi.MPI_FLOAT, source_node, mpi_comm)
		return float(TMP[0])
	elif  type(number_to_send) is bool:
		if number_to_send: number_to_send = 1
		else: number_to_send = 0
		TMP = mpi.mpi_bcast(number_to_send, 1, mpi.MPI_INT, source_node, mpi_comm)
		if TMP == 1:  return True
		else:         return False
	else:
		print(" ERROR in bcast_number_to_all")

def bcast_list_to_all(list_to_send, myid, source_node = 0, mpi_comm = -1):
	pass#IMPORTIMPORTIMPORT from mpi import mpi_bcast, MPI_COMM_WORLD, MPI_FLOAT, MPI_INT
	pass#IMPORTIMPORTIMPORT import   types
	if mpi_comm == -1:  mpi_comm = mpi.MPI_COMM_WORLD
	if(myid == source_node):
		n = len(list_to_send)
		# we will also assume all elements on the list are of the same type
		if( type(list_to_send[0]) == int ): tp = 0
		elif( type(list_to_send[0]) == float ): tp = 1
		else: tp = 2
	else:
		n = 0
		tp = 0
	n = bcast_number_to_all(n, source_node = source_node, mpi_comm = mpi_comm)
	tp = bcast_number_to_all(tp, source_node = source_node, mpi_comm = mpi_comm)
	if( tp == 2 ): 	sparx_global_def.ERROR("Only list of the same type numbers can be brodcasted","bcast_list_to_all",1, myid)
	if(myid != source_node): list_to_send = [0]*n

	if( tp == 0 ):
		list_to_send = mpi.mpi_bcast(list_to_send, n, mpi.MPI_INT, source_node, mpi_comm)
		return [int(n) for n in list_to_send]
	else:
		list_to_send = mpi.mpi_bcast(list_to_send, n, mpi.MPI_FLOAT, source_node, mpi_comm)
		return [float(n) for n in list_to_send]


def recv_attr_dict(main_node, stack, data, list_params, image_start, image_end, number_of_proc, comm = -1):
	pass#IMPORTIMPORTIMPORT import types
	pass#IMPORTIMPORTIMPORT from  utilities import  get_arb_params, set_arb_params
	pass#IMPORTIMPORTIMPORT from  mpi 	import mpi_recv
	pass#IMPORTIMPORTIMPORT from  mpi 	import MPI_FLOAT, MPI_INT, MPI_COMM_WORLD

	#   hdf version!
	# This is done on the main node, so for images from the main node, simply write headers

	if comm == -1:  comm = mpi.MPI_COMM_WORLD

	TransType = type(EMAN2_cppwrap.Transform())
	# prepare keys for float/int
	value = get_arb_params(data[0], list_params)
	ink = []
	len_list = 0
	for il in range(len(list_params)):
		if type(value[il]) is int:
			ink.append(1)
			len_list += 1
		elif type(value[il]) is float:
			ink.append(0)
			len_list += 1
		elif type(value[il]) is TransType:
			ink.append(2)
			len_list += 12
	ldis = []
	headers = []
	for n in range(number_of_proc):
		if n != main_node:
			dis = mpi.mpi_recv(2, mpi.MPI_INT, n, sparx_global_def.SPARX_MPI_TAG_UNIVERSAL, comm)
			value = mpi.mpi_recv(len_list*(dis[1]-dis[0]), mpi.MPI_FLOAT, n, sparx_global_def.SPARX_MPI_TAG_UNIVERSAL, comm)
			ldis.append([dis[0], dis[1]])
			headers.append(value)
			del  dis
	del  value
	for im in range(image_start, image_end):
		data[im-image_start].write_image(stack, data[im-image_start].get_attr_default('ID', im), EMAN2_cppwrap.EMUtil.ImageType.IMAGE_HDF, True)

	for n in range(len(ldis)):
		img_begin = ldis[n][0]
		img_end = ldis[n][1]
		for im in range(img_begin, img_end):
			par_begin = (im-img_begin)*len_list
			nvalue = []
			header = headers[n]
			ilis = 0
			for il in range(len(list_params)):
				if(ink[il] == 1):
					nvalue.append(int(header[par_begin+ilis]))
					ilis += 1
				elif ink[il]==0:
					nvalue.append(float(header[par_begin+ilis]))
					ilis += 1
				else:
					assert ink[il]==2
					t = EMAN2_cppwrap.Transform()
					tmp = []
					for iii in range(par_begin+ilis, par_begin+ilis+12):
						tmp.append(float(header[iii]))
					t.set_matrix(tmp)
					ilis += 12
					nvalue.append(t)
			ISID = list_params.count('ID')
			if(ISID == 0):
				imm = im
			else:
				imm = nvalue[list_params.index('ID')]
			# read head, set params, and write it
			dummy = EMAN2_cppwrap.EMData()
			dummy.read_image(stack, imm, True)
			set_arb_params(dummy, nvalue, list_params)
			dummy.write_image(stack, dummy.get_attr_default('ID', im), EMAN2_cppwrap.EMUtil.ImageType.IMAGE_HDF, True)

def send_attr_dict(main_node, data, list_params, image_start, image_end, comm = -1):
	pass#IMPORTIMPORTIMPORT import types
	pass#IMPORTIMPORTIMPORT from utilities import get_arb_params
	pass#IMPORTIMPORTIMPORT from mpi 	   import mpi_send
	pass#IMPORTIMPORTIMPORT from mpi 	   import MPI_FLOAT, MPI_INT, MPI_COMM_WORLD

	#  This function is called from a node other than the main node

	if comm == -1: comm = mpi.MPI_COMM_WORLD
	TransType = type(EMAN2_cppwrap.Transform())
	mpi.mpi_send([image_start, image_end], 2, mpi.MPI_INT, main_node, sparx_global_def.SPARX_MPI_TAG_UNIVERSAL, comm)
	nvalue = []
	for im in range(image_start, image_end):
		value = get_arb_params(data[im-image_start], list_params)
		for il in range(len(value)):
			if    type(value[il]) is int:  nvalue.append(float(value[il]))
			elif  type(value[il]) is float: nvalue.append(value[il])
			elif  type(value[il]) is TransType:
				m = value[il].get_matrix()
				assert (len(m)==12)
				for f in m: nvalue.append(f)
	mpi.mpi_send(nvalue, len(nvalue), mpi.MPI_FLOAT, main_node, sparx_global_def.SPARX_MPI_TAG_UNIVERSAL, comm)

def recv_attr_dict_bdb(main_node, stack, data, list_params, image_start, image_end, number_of_proc, comm = -1):
	pass#IMPORTIMPORTIMPORT import types
	pass#IMPORTIMPORTIMPORT from  utilities import  get_arb_params, set_arb_params
	pass#IMPORTIMPORTIMPORT from  mpi 	import mpi_recv
	pass#IMPORTIMPORTIMPORT from  mpi 	import MPI_FLOAT, MPI_INT, MPI_COMM_WORLD
	pass#IMPORTIMPORTIMPORT from EMAN2db import db_open_dict
	#  bdb version!
	# This is done on the main node, so for images from the main node, simply write headers

	if comm == -1: comm = mpi.MPI_COMM_WORLD

	DB = EMAN2db.db_open_dict(stack)
	TransType = type(EMAN2_cppwrap.Transform())
	# prepare keys for float/int
	value = get_arb_params(data[0], list_params)
	ink = []
	len_list = 0
	ISID = -1
	for il in range(len(list_params)):
		if(list_params[il] == 'ID'):  ISID = il
		if type(value[il]) is int:
			ink.append(1)
			len_list += 1
		elif type(value[il]) is float:
			ink.append(0)
			len_list += 1
		elif type(value[il]) is TransType:
			ink.append(2)
			len_list += 12
	ldis = []
	headers = []
	for n in range(number_of_proc):
		if n != main_node:
			dis = mpi.mpi_recv(2, mpi.MPI_INT, n, sparx_global_def.SPARX_MPI_TAG_UNIVERSAL, comm)
			img_begin = int(dis[0])
			img_end = int(dis[1])
			header = mpi.mpi_recv(len_list*(img_end-img_begin), mpi.MPI_FLOAT, n, sparx_global_def.SPARX_MPI_TAG_UNIVERSAL, comm)
			for im in range(img_begin, img_end):
				par_begin = (im-img_begin)*len_list
				nvalue = []
				ilis = 0
				for il in range(len(list_params)):
					if(ink[il] == 1):
						nvalue.append(int(header[par_begin+ilis]))
						ilis += 1
					elif ink[il]==0:
						nvalue.append(float(header[par_begin+ilis]))
						ilis += 1
					else:
						assert ink[il]==2
						t = EMAN2_cppwrap.Transform()
						tmp = []
						for iii in range(par_begin+ilis, par_begin+ilis+12):
							tmp.append(float(header[iii]))
						t.set_matrix(tmp)
						ilis += 12
						nvalue.append(t)
				if(ISID == -1):
					imm = im
				else:
					imm = nvalue[ISID]
				for i in range(len(list_params)):
					if(list_params[i] != "ID"):  DB.set_attr(imm, list_params[i], nvalue[i])
		else:
			for n in range(image_start, image_end):
				ID = data[n-image_start].get_attr_default('ID', n)
				for param in list_params:
					if(param != "ID"):  DB.set_attr(ID, param, data[n-image_start].get_attr(param))
	DB.close()

def print_begin_msg(program_name, onscreen=False):
	pass#IMPORTIMPORTIMPORT from time import localtime, strftime
	t = 100
	stars = '*'*t
	string = "Beginning of the program " + program_name + ": " + time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
	s = (t-len(string))/2
	spacing = ' '*s
	if onscreen:
		print(stars)
		print(spacing+string)
		print(stars)
	else:
		print_msg(stars+"\n")
		print_msg(spacing+string+"\n")
		print_msg(stars+"\n")

def print_end_msg(program_name, onscreen=False):
	pass#IMPORTIMPORTIMPORT from time import localtime, strftime
	t = 100
	stars = '*'*t
	string = "End of the program " + program_name + ": " + time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
	s = (t-len(string))/2
	spacing = ' '*s
	if onscreen:
		print(stars)
		print(spacing+string)
		print(stars)
	else:
		print_msg(stars+"\n")
		print_msg(spacing+string+"\n")
		print_msg(stars+"\n")

def print_msg(msg):
	pass#IMPORTIMPORTIMPORT import sys
	pass#IMPORTIMPORTIMPORT import global_def
	if (sparx_global_def.IS_LOGFILE_OPEN == False):
		sparx_global_def.LOGFILE_HANDLE = open(sparx_global_def.LOGFILE,"w")
		sparx_global_def.IS_LOGFILE_OPEN = True
	if (sparx_global_def.BATCH):
		sparx_global_def.LOGFILE_HANDLE.write(msg)
	else:
		sys.stdout.write(msg)
		sparx_global_def.LOGFILE_HANDLE.write(msg)
	sparx_global_def.LOGFILE_HANDLE.flush()

def read_fsc( filename ):
	pass#IMPORTIMPORTIMPORT from string import split, atof
	f = open( filename, 'r' )
	fscc = None
	line = f.readline()
	while len(line) > 0:
		items =  line.split()
		if fscc is None:
			fscc = [None]*len(items)
			for i in range( len(items) ):
				fscc[i] = []

		for i in range( len(items) ) :
			fscc[i].append( string.atof(items[i]) )

		line = f.readline()

	return fscc
"""Multiline Comment7"""

def circumference( img, inner = -1, outer = -1):
	nx = img.get_xsize()
	ny = img.get_ysize()
	nz = img.get_zsize()
	if( inner == -1):
		inner = nx//2 -2
		if( outer <= inner ):  outer = inner + 1
	else:
		if( outer <= inner ):  outer = inner + 1
	inner_sphere = model_circle(inner, nx, ny, nz)

	[mean_a,sigma,imin,imax] = EMAN2_cppwrap.Util.infomask(img, model_circle(outer, nx, ny, nz) - inner_sphere, True)
	inner_rest   = model_blank(nx, ny, nz, 1.0) - inner_sphere
	EMAN2_cppwrap.Util.mul_img(inner_sphere, img)
	return EMAN2_cppwrap.Util.addn_img(inner_sphere, EMAN2_cppwrap.Util.mult_scalar(inner_rest, mean_a ) )

def write_headers(filename, data, lima):
	"""
	  write headers from files in data into a disk file called filename.
	  The filename has to be either hdf or bdb.
	  lima - list with positions in the disk files into which headers will be written,
	    i.e., header from data[k] will be written into file number lima[k]
	  WARNING: this function will open and close DB library!
	"""
	pass#IMPORTIMPORTIMPORT from utilities import file_type
	pass#IMPORTIMPORTIMPORT from EMAN2db import db_open_dict

	ftp = file_type(filename)
	if ftp == "bdb":
		#  For unknown reasons this does not work on Linux, but works on Mac ??? Really?
		DB = EMAN2db.db_open_dict(filename)
		for i in range(len(lima)):
			DB.set_header(lima[i], data[i])
		DB.close()
		#for i in range(len(lima)):
		#	data[i].write_image(filename, lima[i])
	elif ftp == "hdf":
		for i in range(len(lima)):
			data[i].write_image(filename, lima[i], EMAN2_cppwrap.EMUtil.ImageType.IMAGE_HDF, True)
	else:
		sparx_global_def.ERROR("Unacceptable file format","write_headers",1)

def write_header(filename, data, lima):
	"""
	  write header from a single file data into a disk file called filename.
	  The filename has to be either hdf or bdb.
	  lima - position in the disk files into which header will be written,
	    i.e., header from data will be written into file number lima
	  WARNING: this function assums DB library is opened and will NOT close it!
	"""
	pass#IMPORTIMPORTIMPORT from utilities import file_type
	pass#IMPORTIMPORTIMPORT from EMAN2db import db_open_dict

	ftp = file_type(filename)
	if ftp == "bdb":
		DB = EMAN2db.db_open_dict(filename)
		DB.set_header(lima, data)
	elif ftp == "hdf":
		data.write_image(filename, lima, EMAN2_cppwrap.EMUtil.ImageType.IMAGE_HDF, True)
	else:
		sparx_global_def.ERROR("Unacceptable file format","write_headers",1)

def file_type(name):
	if(len(name)>4):
		if(name[:4] == "bdb:"): return "bdb"
		elif(name[-4:-3] == "."):  return name[-3:]
	sparx_global_def.ERROR("Unacceptable file format","file_type",1)

def get_params2D(ima, xform = "xform.align2d"):
	"""
	  retrieve 2D alignment parameters from the header
	  alpha, tx, ty, mirror, scale
	"""
	d = EMAN2_cppwrap.Util.get_transform_params(ima, xform, "2D")
	return d["alpha"],d["tx"],d["ty"],d["mirror"],d["scale"]

def set_params2D(ima, p, xform = "xform.align2d"):
	"""
	  set 2D alignment parameters in the header
	  p = [alpha, tx, ty, mirror, scale]
	"""
	t = EMAN2_cppwrap.Transform({"type":"2D","alpha":p[0],"tx":p[1],"ty":p[2],"mirror":p[3],"scale":p[4]})
	ima.set_attr(xform, t)

def get_params3D(ima, xform = "xform.align3d"):
	"""
	  retrieve 3D alignment parameters from the header
	  phi,theta, psi, tx, ty, tz, mirror,scale
	"""
	d = EMAN2_cppwrap.Util.get_transform_params(ima, xform, "spider")
	return  d["phi"],d["theta"],d["psi"],d["tx"],d["ty"],d["tz"],d["mirror"],d["scale"]

def set_params3D(ima, p, xform = "xform.align3d"):
	"""
	  set 3D alignment parameters in the header
	  p = [phi,theta, psi, tx, ty, tz, mirror,scale]
	"""
	t = EMAN2_cppwrap.Transform({"type":"spider","phi":p[0],"theta":p[1],"psi":p[2],"tx":p[3],"ty":p[4],"tz":p[5],"mirror":p[6],"scale":p[7]})
	ima.set_attr(xform, t)

def get_params_proj(ima, xform = "xform.projection"):
	"""
	  retrieve projection alignment parameters from the header
	  phi, theta, psi, s2x, s2y
	"""
	d = EMAN2_cppwrap.Util.get_transform_params(ima, xform, "spider")
	return  d["phi"],d["theta"],d["psi"],-d["tx"],-d["ty"]

def set_params_proj(ima, p, xform = "xform.projection"):
	"""
	  set projection alignment parameters in the header
	  p = [phi, theta, psi, s2x, s2y]
	"""
	pass#IMPORTIMPORTIMPORT from EMAN2 import Vec2f
	t = EMAN2_cppwrap.Transform({"type":"spider","phi":p[0],"theta":p[1],"psi":p[2]})
	t.set_trans(EMAN2_cppwrap.Vec2f(-p[3], -p[4]))
	ima.set_attr(xform, t)


def get_ctf(ima):
	"""
	  recover numerical values of CTF parameters from EMAN2 CTF object stored in a header of the input image
	  order of returned parameters:
        defocus, cs, voltage, apix, bfactor, ampcont, astigmatism amplitude, astigmatism angle
	"""
	pass#IMPORTIMPORTIMPORT from EMAN2 import EMAN2Ctf
	ctf_params = ima.get_attr("ctf")
	return ctf_params.defocus, ctf_params.cs, ctf_params.voltage, ctf_params.apix, ctf_params.bfactor, ctf_params.ampcont, ctf_params.dfdiff, ctf_params.dfang

def same_ctf(c1,c2):
	"""
	  Compare two CTF objects and return True if they are the same
	"""
	return c1.to_string()==c2.to_string()

def generate_ctf(p):
	"""
	  generate EMAN2 CTF object using values of CTF parameters given in the list p
	  order of parameters:
        p = [defocus, cs, voltage, apix, bfactor, ampcont, astigmatism_amplitude, astigmatism_angle]
	    [ microns, mm, kV, Angstroms, A^2, microns, microns, radians]
	"""
	pass#IMPORTIMPORTIMPORT from EMAN2 import EMAN2Ctf

	defocus      = p[0]
	cs           = p[1]
	voltage      = p[2]
	pixel_size   = p[3]
	bfactor      = p[4]
	amp_contrast = p[5]

	if defocus > 100:  # which means it is very likely in Angstrom, therefore we are using the old convention
		defocus *= 1.0e-4

	ctf = EMAN2_cppwrap.EMAN2Ctf()
	if(len(p) == 6):
		ctf.from_dict({"defocus":defocus, "cs":cs, "voltage":voltage, "apix":pixel_size, "bfactor":bfactor, "ampcont":amp_contrast})
	elif(len(p) == 8):
		ctf.from_dict({"defocus":defocus, "cs":cs, "voltage":voltage, "apix":pixel_size, "bfactor":bfactor, "ampcont":amp_contrast,'dfdiff':p[6],'dfang':p[7]})
	else:
		sparx_global_def.ERROR("Incorrect number of entries on a list, cannot generate CTF","generate_ctf",0)
		return None
	return ctf

def delete_bdb(name):
	"""
	  Delete bdb stack
	"""
	pass#IMPORTIMPORTIMPORT from EMAN2db import db_open_dict, db_remove_dict
	a = EMAN2db.db_open_dict(name)
	EMAN2db.db_remove_dict(name)

def disable_bdb_cache():
	pass#IMPORTIMPORTIMPORT import EMAN2db
	EMAN2db.BDB_CACHE_DISABLE = True

def getvec( phi, tht ):
	pass#IMPORTIMPORTIMPORT from math import radians,cos,sin

	if tht > 180.0:
		tht -= 180.0
		phi += 180.0
	if tht > 90.0:
		tht = 180.0 - tht
		phi += 180.0

	assert tht <=90.0

	qt = numpy.radians(tht)
	qp = numpy.radians(phi)
	qs = numpy.sin(qt)

	x = qs*numpy.cos(qp)
	y = qs*numpy.sin(qp)
	z = numpy.cos(qt)

	return (x,y,z)

def getfvec( phi, tht ):
	pass#IMPORTIMPORTIMPORT from math import radians,cos,sin
	qt = numpy.radians(tht)
	qp = numpy.radians(phi)
	qs = numpy.sin(qt)
	x = qs*numpy.cos(qp)
	y = qs*numpy.sin(qp)
	z = numpy.cos(qt)

	return (x,y,z)

def nearest_fang( vecs, phi, tht ):
	"""
		vecs = [ [x0,y0,z0], [x1,y1,z1], ...]
	"""
	pass#IMPORTIMPORTIMPORT from utilities import getfvec
	vec = getfvec( phi, tht )
	return  EMAN2_cppwrap.Util.nearest_fang(vecs, vec[0],vec[1],vec[2])[0]

def nearest_many_full_k_projangles(reference_normals, angles, howmany = 1, sym_class=None):
	# 
	pass#IMPORTIMPORTIMPORT from utilities import getfvec, angles_to_normals
	#refnormal = normals[:]
	assignments = [-1]*len(angles)
	if( sym_class.sym[:2] == "c1"):
		for i,q in enumerate(angles):
			ref = getfvec(q[0],q[1])
			assignments[i] = EMAN2_cppwrap.Util.nearest_fang_select(reference_normals, ref[0],ref[1],ref[2], howmany)
	else:
		for i,q in enumerate(angles):
			ancordir = angles_to_normals(sym_class.symmetry_neighbors([q[:3]]))
			assignments[i] = EMAN2_cppwrap.Util.nearest_fang_sym(ancordir, reference_normals, len(ancordir), howmany)

	return assignments


def assign_projdirs_f(projdirs, refdirs, neighbors):
	#  projdirs - data
	#  refdirs  - templates, each template has neighbors related copies 
	#  output - list of lists, ofr each of refdirs/neighbors there is a list of projdirs indexes that are closest to it
	"""
	qsti = [-1]*len(projdirs)
	for i,q in enumerate(projdirs):
		dn = -2.0
		for l in xrange(len(refdirs)/neighbors):
			sq = -2.0
			for k in xrange(l*neighbors,(l+1)*neighbors):
				sq = max(q[0]*refdirs[k][0] + q[1]*refdirs[k][1] + q[2]*refdirs[k][2], sq)
			if(sq > dn):
				dn = sq
				this = l
		qsti[i] = this
	"""
	#  Create a list that for each projdirs contains an index of the closest refdirs/neighbors
	qsti = EMAN2_cppwrap.Util.assign_projdirs_f(projdirs, refdirs, neighbors)
	assignments = [[] for i in range(len(refdirs)/neighbors)]
	for i in range(len(projdirs)):
		assignments[qsti[i]].append(i)

	return assignments



def angles_to_normals(angles):
	temp = EMAN2_cppwrap.Util.angles_to_normals(angles)
	return [[temp[l*3+i] for i in range(3)] for l in range(len(angles)) ]
"""Multiline Comment15"""

def angular_occupancy(params, angstep = 15., sym= "c1", method='S'):
	pass#IMPORTIMPORTIMPORT from fundamentals import symclass
	pass#IMPORTIMPORTIMPORT from utilities import nearest_fang, angles_to_normals

	smc  = sparx_fundamentals.symclass(sym)
	eah  = smc.even_angles(angstep, inc_mirror=0, method=method)

	leah = len(eah)
	u = []
	for q in eah:
		#print("q",q)
		m = smc.symmetry_related([(180.0+q[0])%360.0,180.0-q[1],0.0])
		#print("m",m)
		itst = len(u)
		for c in m:
			#print("c",c)
			if smc.is_in_subunit(c[0],c[1],1) :
				#print(" is in 1")
				if not smc.is_in_subunit(c[0],c[1],0) :
					#print("  outside")
					u.append(c)
					break
		if(len(u) != itst+1):
			u.append(q)  #  This is for exceptions that cannot be easily handled
			"""Multiline Comment16"""
	seaf = []
	for q in eah+u:  seaf += smc.symmetry_related(q)

	lseaf = len(seaf)/(2*leah)
	#print(lseaf)
	#for i,q in enumerate(seaf):  print(" seaf  ",i,q)
	#print(seaf)
	seaf = angles_to_normals(seaf)

	occupancy = [[] for i in range(leah)]

	for i,q in enumerate(params):
		l = nearest_fang(seaf,q[0],q[1])
		l = l/lseaf
		if(l>=leah):  l = l-leah
		occupancy[l].append(i)
	#for i,q in enumerate(occupancy):
	#	if q:
	#		print("  ",i,q,eah[i])
	return occupancy, eah


def angular_histogram(params, angstep = 15., sym= "c1", method='S'):
	occupancy, eah = angular_occupancy(params, angstep, sym, method)
	return  [len(q) for q in occupancy], eah

def balance_angular_distribution(params, max_occupy = -1, angstep = 15., sym= "c1"):
	pass#IMPORTIMPORTIMPORT from fundamentals import symclass
	occupancy,eah = angular_occupancy(params, angstep, sym, method='S')

	if(max_occupy > 0):
		outo = []
		pass#IMPORTIMPORTIMPORT from random import shuffle
		for l,q in enumerate(occupancy):
			random.shuffle(q)
			q = q[:max_occupy]
			outo += q
			print("  %10d   %10d        %6.1f   %6.1f"%(l,len(q),eah[l][0],eah[l][1]))
			#print(l,len(q),q)
		outo.sort()
		#write_text_file(outo,"select.txt")
		return outo
	else:
		#for l,q in enumerate(occupancy):
		#	print("  %10d   %10d        %6.1f   %6.1f"%(l,len(q),eah[l][0],eah[l][1]))
		return occupancy


def symmetry_neighbors(angles, symmetry):
	#  input is a list of lists  [[phi0,theta0,psi0],[phi1,theta1,psi1],...]
	#  output is [[phi0,theta0,psi0],[phi0,theta0,psi0]_SYM1,...,[phi1,theta1,psi1],[phi1,theta1,psi1]_SYM1,...]
	temp = EMAN2_cppwrap.Util.symmetry_neighbors(angles, symmetry)
	nt = len(temp)/3
	return [[temp[l*3+i] for i in range(3)] for l in range(nt) ]
	#  We could make it a list of lists
	#mt = len(angles)
	#nt = len(temp)/mt/3
	#return [[ [temp[m*3*nt+l*3+i] for i in xrange(3)] for l in xrange(nt)] for m in xrange(mt) ]

#def nearest_angular_direction(normals, vect, symmetry):

def rotation_between_anglesets(agls1, agls2):
	"""
	  Find an overall 3D rotation (phi theta psi) between two sets of Eulerian angles. (psi irrelevant)
	  The two sets have to have the same number of elements and it is assumed that k'th element on the first
	  list corresponds to the k'th element on the second list.
	  Input: two lists [[phi1, theta1, psi1], [phi2, theta2, psi2], ...].  Second list is considered reference.
	  Output: overall rotation phi, theta, psi that has to be applied to the first list (agls1) so resulting
	    angles will agree with the second list.
	  Note: all angles have to be in spider convention.
	  For details see: Appendix in Penczek, P., Marko, M., Buttle, K. and Frank, J.:  Double-tilt electron tomography.  Ultramicroscopy 60:393-410, 1995.
	"""
	pass#IMPORTIMPORTIMPORT from math  import sin, cos, pi, sqrt, atan2, acos, atan, radians
	pass#IMPORTIMPORTIMPORT from numpy import array, linalg, matrix
	pass#IMPORTIMPORTIMPORT import types

	def ori2xyz(ori):
		if(type(ori) == list):
			phi, theta, psi = ori[:3]
		else:
			# it has to be Transformation object
			d = ori.get_params("spider")
			phi   = d["phi"]
			theta = d["theta"]
			#psi   = d["psi"]

		phi   = numpy.radians(phi)
		theta = numpy.radians(theta)
		sint = numpy.sin(theta)
		x = sint * numpy.sin(phi)
		y = sint * numpy.cos(phi)
		z = numpy.cos(theta)

		return [x, y, z]

	N = len(agls1)
	if N != len(agls2):
		sparx_global_def.ERROR('rotation_between_anglesets', 'Both lists must have the same length',1)
		return -1
	if N < 2:
		sparx_global_def.ERROR('rotation_between_anglesets',  'At least two orientations are required in each list',1)
		return -1

	U1 = [ori2xyz(q) for q in agls1]
	U2 = [ori2xyz(q) for q in agls2]

	# compute all Suv with uv = {xx, xy, xz, yx, ..., zz}
	Suv   = [0] * 9

	nbori = len(U1)
	for i in range(3):
		for j in range(3):
			for s in range(nbori):
				Suv[j+3*i] += (U2[s][i] * U1[s][j])

	# create matrix N
	N = numpy.array([[Suv[0]+Suv[4]+Suv[8], Suv[5]-Suv[7],    Suv[6]-Suv[2],                 Suv[1]-Suv[3]],
		   [Suv[5]-Suv[7],        Suv[0]-Suv[4]-Suv[8], Suv[1]+Suv[3],                 Suv[6]+Suv[2]],
		   [Suv[6]-Suv[2],        Suv[1]+Suv[3],        -Suv[0]+Suv[4]-Suv[8],         Suv[5]+Suv[7]],
		   [Suv[1]-Suv[3],        Suv[6]+Suv[2],        Suv[5]+Suv[7],         -Suv[0]-Suv[4]+Suv[8]]])

	# eigenvector corresponding to the most positive eigenvalue
	val, vec = numpy.linalg.eig(N)
	q0, qx, qy, qz = vec[:, val.argmax()]
	# create quaternion Rot matrix
	r = [
		[q0*q0-qx*qx+qy*qy-qz*qz,         2*(qy*qx+q0*qz),          2*(qy*qz-q0*qx)],
		[2*(qx*qy-q0*qz),                 q0*q0+qx*qx-qy*qy-qz*qz,  2*(qx*qz+q0*qy)],
		[2*(qz*qy+q0*qx),                 2*(qz*qx-q0*qy),          q0*q0-qx*qx-qy*qy+qz*qz],
		]

	pass#IMPORTIMPORTIMPORT from fundamentals import recmat
	return  sparx_fundamentals.recmat(r)

def angle_between_projections_directions(proj1, proj2):
	"""
	  It returns angle between two projections directions.
	  INPUT: two lists: [phi1, theta1] , [phi2, theta2]
	  OUTPUT: angle (in degrees)
	"""
	pass#IMPORTIMPORTIMPORT from math import sin, cos, acos, radians, degrees
	pass#IMPORTIMPORTIMPORT from utilities import lacos
	theta1 = numpy.radians(proj1[1])
	theta2 = numpy.radians(proj2[1])
	cp1cp2_sp1sp2 = numpy.cos(numpy.radians(proj1[0]) - numpy.radians(proj2[0]))
	temp = numpy.sin(theta1) * numpy.sin(theta2) * cp1cp2_sp1sp2 + numpy.cos(theta1) * numpy.cos(theta2)
	return lacos( temp )

getang3 = angle_between_projections_directions

def get_pixel_size(img):
	"""
	  Retrieve pixel size from the header.
	  We check attribute Pixel_size and also pixel size from ctf object, if exisits.
	  If the two are different or if the pixel size is not set, return -1.0 and print a warning.
	"""
	p1 = img.get_attr_default("apix_x", -1.0)
	cc = img.get_attr_default("ctf", None)
	if cc == None:
		p2 = -1.0
	else:
		p2 = round(cc.apix, 3)
	if p1 == -1.0 and p2 == -1.0:
		#ERROR("Pixel size not set", "get_pixel_size", 0)
		return -1.0
	elif p1 > -1.0 and p2 > -1.0:
		#if abs(p1-p2) >= 0.001:
		#	ERROR("Conflict between pixel size in attribute and in ctf object", "get_pixel_size", 0)
		# pixel size is positive, so what follows omits -1 problem
		return max(p1, p2)
	else:
		return max(p1, p2)

def set_pixel_size(img, pixel_size):
	"""
	  Set pixel size in the header.
	  Set attribute Pixel_size and also pixel size in ctf object, if exists.
	"""
	nz = img.get_zsize()
	img.set_attr("apix_x", round(pixel_size, 3))
	img.set_attr("apix_y", round(pixel_size, 3))
	img.set_attr("apix_z", round(pixel_size, 3))
	cc = img.get_attr_default("ctf", None)
	if(cc):
		cc.apix = pixel_size
		img.set_attr("ctf", cc)

def lacos(x):
	"""
		compute acos(x) in degrees after enforcing -1<=x<=1
	"""
	pass#IMPORTIMPORTIMPORT from math import degrees, acos
	return  numpy.degrees(math.acos(max(-1.0,min(1.0,x))))

def nearest_proj(proj_ang, img_per_grp=100, List=[]):
	pass#IMPORTIMPORTIMPORT from utilities import getfvec
	pass#IMPORTIMPORTIMPORT from math import exp, pi
	pass#IMPORTIMPORTIMPORT from sets import Set
	pass#IMPORTIMPORTIMPORT from time import time
	pass#IMPORTIMPORTIMPORT from random import randint

	def ang_diff(v1, v2):
		# The first return value is the angle between two vectors
		# The second return value is whether we need to mirror one of them (0 - no need, 1 - need)
		pass#IMPORTIMPORTIMPORT from math import acos, degrees
		pass#IMPORTIMPORTIMPORT from utilities import lacos

		v = v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2]
		if v >= 0: return lacos(v), 0
		else:  return lacos(-v), 1

	def get_ref_ang_list(delta, sym):
		ref_ang = even_angles(delta, symmetry=sym)
		ref_ang_list = [0.0]*(len(ref_ang)*2)
		for i in range(len(ref_ang)):
			ref_ang_list[2*i] = ref_ang[i][0]
			ref_ang_list[2*i+1] = ref_ang[i][1]
		return ref_ang_list, len(ref_ang)

	def binary_search(a, x):
		N = len(a)
		begin = 0
		end = N-1
		while begin <= end:
			mid = (begin+end)/2
			if a[mid] == x: return mid
			if a[mid] < x: begin = mid+1
			else: end = mid-1
		return -1

	def binary_search_l(a, x):
		# This function returns an index i such that i is the smallest number
		# such that when t >= i, a[t] >= x
		N = len(a)
		t = binary_search(a, x)
		if t != -1:
			while t-1 >= 0 and a[t-1] == a[t]: t -= 1
			return t
		else:
			if x > a[N-1]: return -1;
			if x < a[0]: return 0;
			begin = 0
			end = N-2
			while end >= begin:
				mid = (begin+end)/2
				if x > a[mid] and x < a[mid+1]: break;
				if x < a[mid]: end = mid-1
				else: begin = mid+1
			return mid+1

	def binary_search_r(a, x):
		# This function returns an index i such that i is the largest number
		# such that when t <= i, a[t] <= x
		N = len(a)
		t = binary_search(a, x)
		if t != -1:
			while t+1 <= N-1 and a[t+1] == a[t]: t += 1
			return t
		else:
			if x > a[N-1]: return N-1;
			if x < a[0]: return -1;
			begin = 0
			end = N-2
			while end >= begin:
				mid = (begin+end)/2
				if x > a[mid] and x < a[mid+1]: break;
				if x < a[mid]: end = mid-1
				else: begin = mid+1
			return mid

	N = len(proj_ang)
	if len(List) == 0: List = list(range(N))
	if N < img_per_grp:
		print("Error: image per group larger than the number of particles!")
		exit()
	phi_list   = [[0.0, 0] for i in range(N)]
	theta_list = [[0.0, 0] for i in range(N)]
	vec = [None]*N
	for i in range(N):
		phi = proj_ang[i][0]
		theta = proj_ang[i][1]
		vec[i] = getfvec(phi, theta)
		if theta > 90.0:
			theta = 180.0-theta
			phi  += 180.0
		phi = phi%360.0
		phi_list[i][0]   = phi
		phi_list[i][1]   = i
		theta_list[i][0] = theta
		theta_list[i][1] = i
	theta_list.sort()
	phi_list.sort()
	theta_list_l = [0.0]*N
	phi_list_l = [0.0]*N
	for i in range(N):
		theta_list_l[i] = theta_list[i][0]
		phi_list_l[i] = phi_list[i][0]

	g = [[360.0, 0, 0] for i in range(N)]
	proj_list = []
	mirror_list = []
	neighbor   = [0]*img_per_grp
	#neighbor2 = [0]*img_per_grp
	dis        = [0.0]*img_per_grp
	#dis2      = [0.0]*img_per_grp
	mirror     = [0]*img_per_grp
	S = sets.Set()
	T = sets.Set()
	#tt1 = time()
	for i in range(len(List)):
		k = List[i]
		#print "\nCase #%3d: Testing projection %6d"%(i, k)
		#t1 = time()
		phi = proj_ang[k][0]
		theta = proj_ang[k][1]
		if theta > 90.0:
			theta = 180.0-theta
			phi += 180.0
		phi = phi%360.0
		delta = 0.01
		while True:
			min_theta = max( 0.0, theta-delta)
			max_theta = min(90.0, theta+delta)
			if min_theta == 0.0:
				min_phi = 0.0
				max_phi = 360.0
			else:
				dphi = min(delta/(2*min_theta)*180.0, 180.0)
				min_phi = phi - dphi
				max_phi = phi + dphi
				if min_phi < 0.0: min_phi += 360.0
				if max_phi > 360.0: max_phi -= 360.0
				if theta+delta > 90.0:
					phi_mir = (phi+180.0)%360.0
					min_phi_mir = phi_mir - dphi
					max_phi_mir = phi_mir + dphi
					if min_phi_mir < 0.0: min_phi_mir += 360.0
					if max_phi_mir > 360.0: max_phi_mir -= 360.0

			phi_left_bound    = binary_search_l(phi_list_l, min_phi)
			phi_right_bound   = binary_search_r(phi_list_l, max_phi)
			theta_left_bound  = binary_search_l(theta_list_l, min_theta)
			theta_right_bound = binary_search_r(theta_list_l, max_theta)
			if theta+delta > 90.0:
				phi_mir_left_bound = binary_search_l(phi_list_l, min_phi_mir)
				phi_mir_right_bound = binary_search_r(phi_list_l, max_phi_mir)
			#print delta
			#print min_phi, max_phi, min_theta, max_theta
			#print phi_left_bound, phi_right_bound, theta_left_bound, theta_right_bound
			if phi_left_bound < phi_right_bound:
				for j in range(phi_left_bound, phi_right_bound+1):
					S.add(phi_list[j][1])
			else:
				for j in range(phi_right_bound+1):
					S.add(phi_list[j][1])
				for j in range(phi_left_bound, N):
					S.add(phi_list[j][1])
			if theta+delta > 90.0:
				if phi_mir_left_bound < phi_mir_right_bound:
					for j in range(phi_mir_left_bound, phi_mir_right_bound+1):
						S.add(phi_list[j][1])
				else:
					for j in range(phi_mir_right_bound+1):
						S.add(phi_list[j][1])
					for j in range(phi_mir_left_bound, N):
						S.add(phi_list[j][1])
			for j in range(theta_left_bound, theta_right_bound+1):
				T.add(theta_list[j][1])
			v = list(T.intersection(S))
			S.clear()
			T.clear()
			if len(v) >= min(1.5*img_per_grp, N): break
			delta *= 2
			del v

		for j in range(len(v)):
			d = ang_diff(vec[v[j]], vec[k])
			g[j][0] = d[0]
			if v[j] == k: g[j][0] = -1.  # To ensure the image itself is always included in the group
			g[j][1] = d[1]
			g[j][2] = v[j]
		g[:len(v)] = sorted(g[:len(v)])
		for j in range(img_per_grp):
			neighbor[j] = g[j][2]
			dis[j] = g[j][0]
			mirror[j] = (g[j][1] == 1)
		proj_list.append(neighbor[:])
		mirror_list.append(mirror[:])
		#t2 = time()

		"""Multiline Comment18"""
	#tt2 = time()
	#print tt2-tt1
	return proj_list, mirror_list


def findall(value, L, start=0):
	"""
	 return a list of all indices of a value on the list L beginning from position start
	"""
	positions = []
	lL = len(L) - 1
	i = start - 1
	while( i < lL  ):
		i += 1
		try:
			i = L.index(value, i)
			positions.append(i)
		except:
			pass
	return positions

"""Multiline Comment19"""

class iterImagesList(object):
	images = []
	imagesIndexes = []
	position = -1
	def __init__(self, list_of_images, list_of_indexes = None):
		if list_of_indexes == None:
			self.images = list_of_images[:]
			self.imagesIndexes = list(range(len(self.images)))
		else:
			for i in list_of_indexes:
				self.images.append(list_of_images[i])
			self.imagesIndexes = list_of_indexes[:]
	def iterNo(self):
		return self.position
	def imageIndex(self):
		return self.imagesIndexes[self.position]
	def image(self):
		return self.images[self.position]
	def goToNext(self):
		if len(self.imagesIndexes) <= self.position:
			return False
		self.position += 1
		return (self.position < len(self.imagesIndexes))
	def goToPrev(self):
		if 0 > self.position:
			return False
		self.position -= 1
		return (self.position >= 0)

# ================ Iterator for stack of images
def pack_message(data):
	"""Convert data for transmission efficiently"""

	if isinstance(data,str):
		if len(data)>256 : return "C"+zlib.compress(data,1)
		else : return "S"+data
	else :
		d2x=pickle.dumps(data,-1)
		if len(d2x)>256 : return "Z"+zlib.compress(d2x,1)
		else : return "O"+d2x


def unpack_message(msg):
	"""Unpack a data payload prepared by pack_message"""

	if msg[0]=="C" : return zlib.decompress((msg[1:]).tostring())
	elif msg[0]=="S" : return (msg[1:]).tostring()
	elif msg[0]=="Z" : return pickle.loads(zlib.decompress((msg[1:]).tostring()))
	elif msg[0]=="O" : return pickle.loads((msg[1:]).tostring())
	else :
		print("ERROR: Invalid MPI message. Please contact developers. (%s)"%str(msg[:20]))
		raise Exception("unpack_message")


statistics_send_recv = dict()

def update_tag(communicator, target_rank):   # TODO - it doesn't work when communicators are destroyed and recreated
	return 123456
	global statistics_send_recv
	if communicator not in statistics_send_recv:
		pass#IMPORTIMPORTIMPORT from mpi import mpi_comm_size
		statistics_send_recv[communicator] = [0] * mpi.mpi_comm_size(communicator)
	statistics_send_recv[communicator][target_rank] += 1
	return statistics_send_recv[communicator][target_rank]

# ===================================== WRAPPER FOR MPI

def wrap_mpi_send(data, destination, communicator = None):
	pass#IMPORTIMPORTIMPORT from mpi import mpi_send, MPI_COMM_WORLD, MPI_CHAR

	if communicator == None:
		communicator = mpi.MPI_COMM_WORLD

	msg = pack_message(data)
	tag = update_tag(communicator, destination)
	#from mpi import mpi_comm_rank
	#print communicator, mpi_comm_rank(communicator), "send to", destination, tag
	mpi.mpi_send(msg, len(msg), mpi.MPI_CHAR, destination, tag, communicator) # int MPI_Send( void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm )


def wrap_mpi_recv(source, communicator = None):
	pass#IMPORTIMPORTIMPORT from mpi import mpi_recv, MPI_COMM_WORLD, MPI_CHAR, mpi_probe, mpi_get_count

	if communicator == None:
		communicator = mpi.MPI_COMM_WORLD

	tag = update_tag(communicator, source)
	#from mpi import mpi_comm_rank
	#print communicator, mpi_comm_rank(communicator), "recv from", source, tag
	mpi.mpi_probe(source, tag, communicator)
	n = mpi.mpi_get_count(mpi.MPI_CHAR)
	msg = mpi.mpi_recv(n, mpi.MPI_CHAR, source, tag, communicator)
	return unpack_message(msg)


def wrap_mpi_bcast(data, root, communicator = None):
	pass#IMPORTIMPORTIMPORT from mpi import mpi_bcast, MPI_COMM_WORLD, mpi_comm_rank, MPI_CHAR

	if communicator == None:
		communicator = mpi.MPI_COMM_WORLD

	rank = mpi.mpi_comm_rank(communicator)

	if rank == root:
		msg = pack_message(data)
		n = struct.pack("I",len(msg))
	else:
		msg = None
		n = None

	n = mpi.mpi_bcast(n, 4, mpi.MPI_CHAR, root, communicator)  # int MPI_Bcast ( void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm )
	n=struct.unpack("I",n)[0]
	msg = mpi.mpi_bcast(msg, n, mpi.MPI_CHAR, root, communicator)  # int MPI_Bcast ( void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm )
	return unpack_message(msg)


# data must be a python list (numpy array also should be implemented)
def wrap_mpi_gatherv(data, root, communicator = None):
	pass#IMPORTIMPORTIMPORT from mpi import mpi_comm_rank, mpi_comm_size, MPI_COMM_WORLD

	if communicator == None:
		communicator = mpi.MPI_COMM_WORLD

	rank = mpi.mpi_comm_rank(communicator)
	procs = mpi.mpi_comm_size(communicator)

	out_array = None
	if rank == root:
		if type(data) is list:
			out_array = []
			for p in range(procs):
				if p == rank:
					out_array.extend(data)
				else:
					recv_data = wrap_mpi_recv(p, communicator)
					out_array.extend(recv_data)
		else:
			raise Exception("wrap_mpi_gatherv: type of data not supported")
	else:
		wrap_mpi_send(data, root, communicator)

	return out_array

# def rearrange_ranks_of_processors_to_fit_round_robin_assignment():
def get_colors_and_subsets(main_node, mpi_comm, my_rank, shared_comm, sh_my_rank, masters_from_groups_vs_everything_else_comm):
	"""
	It is assumed that this code or equivalent is ran before calling this function

	mpi_init(0, [])
	
	mpi_comm = MPI_COMM_WORLD
	main_node = 0
	
	my_rank = mpi_comm_rank(mpi_comm)
	mpi_size = mpi_comm_size(mpi_comm)
	
	shared_comm = mpi_comm_split_type(mpi_comm, MPI_COMM_TYPE_SHARED,  0, MPI_INFO_NULL)
	
	sh_my_rank = mpi_comm_rank(shared_comm)
	key = sh_my_rank
	sh_mpi_size = mpi_comm_size(shared_comm)
	
	masters_from_groups_vs_everything_else_comm = mpi_comm_split(mpi_comm, sh_my_rank == main_node, my_rank)
	"""

	pass#IMPORTIMPORTIMPORT from mpi import mpi_barrier

	sh_my_ranks = [my_rank]
	sh_my_ranks = wrap_mpi_gatherv(sh_my_ranks, main_node, shared_comm)

	group_infos = None
	if sh_my_rank == main_node:
		# executed only by masters from groups
		group_infos = [my_rank, sh_my_ranks]
		group_infos = wrap_mpi_gatherv(group_infos, main_node, masters_from_groups_vs_everything_else_comm)

	mpi.mpi_barrier(mpi_comm)

	group_infos = wrap_mpi_bcast(group_infos, main_node, mpi_comm)

	number_of_groups = len(group_infos)/2

	for i in range(number_of_groups):
		if my_rank in group_infos[2*i+1]:
			color = i#group_infos[2*i]
			break

	number_of_processes_in_each_group = []
	for i in range(number_of_groups):
		number_of_processes_in_each_group.append(len(group_infos[2*i+1]))

	balanced_processor_load_on_nodes = len(set(number_of_processes_in_each_group)) == 1 

	return color, number_of_groups, balanced_processor_load_on_nodes


def wrap_mpi_split(comm, no_of_groups):
	"""

	Takes the processes of a communicator (comm) and splits them in groups (no_of_groups).
	Each subgroup of processes has ids generated from 0 to number of processes per group - 1.
	Consecutive global process ids have consecutive subgroup process ids.

	"""
	pass#IMPORTIMPORTIMPORT from mpi import mpi_comm_size, mpi_comm_rank, mpi_comm_split
	nproc = mpi.mpi_comm_size(comm)
	myid = mpi.mpi_comm_rank(comm)

	no_of_proc_per_group = nproc / no_of_groups
	color = myid / no_of_proc_per_group
	key = myid % no_of_proc_per_group

	return mpi.mpi_comm_split(comm, color, key)
	

def get_dist(c1, c2):
	pass#IMPORTIMPORTIMPORT from math import sqrt
	d = numpy.sqrt((c1[0] - c2[0])**2 + (c1[1] - c2[1])**2)
	return d


def eliminate_moons(my_volume, moon_elimination_params):
	"""
	moon_elimination_params[0] - mass in KDa
	moon_elimination_params[1] - pixel size in A
	"""

	pass#IMPORTIMPORTIMPORT from morphology import binarize
	histogram_threshold  =  my_volume.find_3d_threshold(moon_elimination_params[0], moon_elimination_params[1])*1.1
	# clean11 88x88,  4.84 px/A 750 kDa

	my_volume_binarized = sparx_morphology.binarize(my_volume, histogram_threshold)
	# my_volume_binarized.write_image ("my_volume_binarized.hdf")
	my_volume_binarized_with_no_moons = EMAN2_cppwrap.Util.get_biggest_cluster(my_volume_binarized)
	# my_volume_binarized_with_no_moons.write_image("my_volume_binarized_with_no_moons.hdf")
	volume_difference = my_volume_binarized - my_volume_binarized_with_no_moons
	# volume_difference.write_image("volume_difference.hdf")

	if volume_difference.get_value_at(volume_difference.calc_max_index()) == 0 and \
		volume_difference.get_value_at(volume_difference.calc_min_index()) == 0:
		return my_volume
	else:
		pass#IMPORTIMPORTIMPORT from utilities import gauss_edge
		return gauss_edge(my_volume_binarized_with_no_moons) * my_volume

		# from utilities   import model_blank
		# # mask = model_blank(my_volume_binarized_with_no_moons.get_xsize(), my_volume_binarized_with_no_moons.get_ysize(), my_volume_binarized_with_no_moons.get_zsize())
		# # mask.to_one()
	# this is only in master

def combinations_of_n_taken_by_k(n, k):
	pass#IMPORTIMPORTIMPORT from fractions import Fraction
	return int(reduce(lambda x, y: x * y, (fractions.Fraction(n-i, i+1) for i in range(k)), 1))

def cmdexecute(cmd, printing_on_success = True):
	pass#IMPORTIMPORTIMPORT from   time import localtime, strftime
	pass#IMPORTIMPORTIMPORT import subprocess
	pass#IMPORTIMPORTIMPORT import os
	#outcome = subprocess.call(cmd, shell=True)
	outcome = os.system(cmd)
	line = time.strftime("%Y-%m-%d_%H:%M:%S", time.localtime()) + " =>"
	if(outcome != 0):
		print(line,"ERROR!!   Command failed:  ", cmd, " return code of failed command: ", outcome)
		return 0
	elif printing_on_success:
		print(line,"Executed successfully: ",cmd)
		return 1

def string_found_in_file(myregex, filename):
	pass#IMPORTIMPORTIMPORT import re
	pattern = re.compile(myregex)
	for line in open(filename):
		if re.findall(pattern, line) != []:
			return True
	return False

def get_latest_directory_increment_value(directory_location, directory_name, start_value = 1, myformat = "%03d"):
	pass#IMPORTIMPORTIMPORT import os
	dir_count = start_value
	while os.path.isdir(directory_location + directory_name + myformat%(dir_count)):
		dir_count += 1
	if dir_count == start_value:
		return start_value
	return dir_count - 1

def if_error_then_all_processes_exit_program(error_status):
	pass#IMPORTIMPORTIMPORT import sys, os
	pass#IMPORTIMPORTIMPORT from utilities import print_msg

	# if "OMPI_COMM_WORLD_SIZE" not in os.environ:
	if len({"OMPI_COMM_WORLD_SIZE", "PMI_RANK", "PMI_ID", "SLURM_PROCID", "LAMRANK", "MPI_RANKID", "MP_CHILD", "MP_CHILD", "MP_CHILD"}.intersection(set(os.environ))) == 0:
		def my_mpi_comm_rank(n): return 0
		def my_mpi_bcast(*largs):
			return [largs[0]]
		def my_mpi_finalize():
			return None
		MY_MPI_INT, MY_MPI_COMM_WORLD = 0, 0
	else:
		pass#IMPORTIMPORTIMPORT from mpi import mpi_comm_rank, mpi_bcast, mpi_finalize, MPI_INT, MPI_COMM_WORLD
		my_mpi_comm_rank = mpi.mpi_comm_rank
		my_mpi_bcast = mpi.mpi_bcast
		my_mpi_finalize = mpi.mpi_finalize
		MY_MPI_INT = mpi.MPI_INT
		MY_MPI_COMM_WORLD = mpi.MPI_COMM_WORLD

	myid = my_mpi_comm_rank(MY_MPI_COMM_WORLD)
	if error_status != None and error_status != 0:
		error_status_info = error_status
		error_status = 1
	else:
		error_status = 0

	error_status = my_mpi_bcast(error_status, 1, MY_MPI_INT, 0, MY_MPI_COMM_WORLD)
	error_status = int(error_status[0])

	if error_status > 0:
		if myid == 0:
			if type(error_status_info) == type((1,1)):
				if len(error_status_info) == 2:
					frameinfo = error_status_info[1]
					print_msg("***********************************\n")
					print_msg("** Error: %s\n"%error_status_info[0])
					print_msg("***********************************\n")
					print_msg("** Location: %s\n"%(frameinfo.filename + ":" + str(frameinfo.lineno)))
					print_msg("***********************************\n")
		sys.stdout.flush()
		my_mpi_finalize()
		exit() # sys.exit(1)

def get_shrink_data_huang(Tracker, nxinit, partids, partstack, myid, main_node, nproc, preshift = False):
	# The function will read from stack a subset of images specified in partids
	#   and assign to them parameters from partstack with optional CTF application and shifting of the data.
	# So, the lengths of partids and partstack are the same.
	#  The read data is properly distributed among MPI threads.
	# 10142015 --- preshift is set to True when doing 3-D sorting.
	# chunk_id are set when data is read in

	pass#IMPORTIMPORTIMPORT from fundamentals import resample, fshift
	pass#IMPORTIMPORTIMPORT from filter import filt_ctf
	pass#IMPORTIMPORTIMPORT from applications import MPI_start_end

	"""Multiline Comment21"""
	if( myid == main_node ): lpartids = read_text_file(partids)
	else:  lpartids = 0
	lpartids = wrap_mpi_bcast(lpartids, main_node)
	ndata = len(lpartids)
	if( myid == main_node ):  partstack = read_text_row(partstack)
	else:  partstack = 0
	partstack = wrap_mpi_bcast(partstack, main_node)
	if( ndata < nproc):
		if(myid<ndata):
			image_start = myid
			image_end   = myid+1
		else:
			image_start = 0
			image_end   = 1
	else:
		image_start, image_end = sparx_applications.MPI_start_end(ndata, nproc, myid)
	lpartids  = lpartids[image_start:image_end]
	#partstack = partstack[image_start:image_end]
	#  Preprocess the data
	mask2D  = model_circle(Tracker["constants"]["radius"],Tracker["constants"]["nnxo"],Tracker["constants"]["nnxo"])
	nima = image_end - image_start
	oldshifts = [[0.0,0.0]]#*nima
	data = [None]*nima
	shrinkage = nxinit/float(Tracker["constants"]["nnxo"])
	radius = int(Tracker["constants"]["radius"] * shrinkage +0.5)
	#  Note these are in Fortran notation for polar searches
	#txm = float(nxinit-(nxinit//2+1) - radius -1)
	#txl = float(2 + radius - nxinit//2+1)
	txm = float(nxinit-(nxinit//2+1) - radius)
	txl = float(radius - nxinit//2+1)
	for im in range(nima):
		data[im] = get_im(Tracker["constants"]["stack"], lpartids[im])
		if im ==0:
			if data[im].get_xsize() > Tracker["constants"]["nnxo"]:
				window_particle =True
				pass#IMPORTIMPORTIMPORT from EMAN2 import Region
			else:
				window_particle =False
		phi,theta,psi,sx,sy = partstack[lpartids[im]][0], partstack[lpartids[im]][1], partstack[lpartids[im]][2], partstack[lpartids[im]][3], partstack[lpartids[im]][4]
		if( Tracker["constants"]["CTF"] and Tracker["applyctf"] ):
			ctf_params = data[im].get_attr("ctf")
			st = EMAN2_cppwrap.Util.infomask(data[im], mask2D, False)
			data[im] -= st[0]
			data[im] = sparx_filter.filt_ctf(data[im], ctf_params)
			data[im].set_attr('ctf_applied', 1)
		if preshift:# always true
			data[im] = sparx_fundamentals.fshift(data[im], sx, sy)
			set_params_proj(data[im],[phi,theta,psi,0.0,0.0])
			sx = 0.0
			sy = 0.0
		if window_particle:
			mx = data[im].get_xsize()//2-Tracker["constants"]["nnxo"]//2
			my = data[im].get_ysize()//2-Tracker["constants"]["nnxo"]//2
			data[im] = data[im].get_clip(EMAN2_cppwrap.Region(mx,my,Tracker["constants"]["nnxo"],Tracker["constants"]["nnxo"]))
			data[im].set_attr('ctf_applied', 1)
			set_params_proj(data[im],[phi,theta,psi,0.0,0.0])
		#oldshifts[im] = [sx,sy]
		#  resample will properly adjusts shifts and pixel size in ctf
		data[im] = sparx_fundamentals.resample(data[im], shrinkage)
		#  We have to make sure the shifts are within correct range, shrinkage or not
		set_params_proj(data[im],[phi,theta,psi,max(min(sx*shrinkage,txm),txl),max(min(sy*shrinkage,txm),txl)])
		#  For local SHC set anchor
		#if(nsoft == 1 and an[0] > -1):
		#  We will always set it to simplify the code
		set_params_proj(data[im],[phi,theta,psi,0.0,0.0], "xform.anchor")
		chunk_id_state = data[im].get_attr_default("chunk_id",None)
		if chunk_id_state == None: data[im].set_attr("chunk_id",Tracker["chunk_dict"][ lpartids[im]])
	assert( nxinit == data[0].get_xsize() )  #  Just to make sure.
	#oldshifts = wrap_mpi_gatherv(oldshifts, main_node, MPI_COMM_WORLD)
	return data, oldshifts

"""Multiline Comment22"""
def getindexdata(stack, partids, partstack, myid, nproc):
	# The function will read from stack a subset of images specified in partids
	#   and assign to them parameters from partstack
	# So, the lengths of partids and partstack are the same.
	#  The read data is properly distributed among MPI threads.

	pass#IMPORTIMPORTIMPORT from applications import MPI_start_end

	lpartids = read_text_file(partids)
	ndata = len(lpartids)
	partstack = read_text_row(partstack)

	if( ndata < nproc):
		if(myid<ndata):
			image_start = myid
			image_end   = myid+1
		else:
			image_start = 0
			image_end   = 1
	else:
		image_start, image_end = sparx_applications.MPI_start_end(ndata, nproc, myid)
	lpartids  = lpartids[image_start:image_end]
	partstack = partstack[image_start:image_end]
	data = EMAN2_cppwrap.EMData.read_images(stack, lpartids)

	for i in range(len(partstack)):
		set_params_proj(data[i], partstack[i])
	return data


def store_value_of_simple_vars_in_json_file(filename, local_vars, exclude_list_of_vars = [], write_or_append = "w",
	vars_that_will_show_only_size = []):

	pass#IMPORTIMPORTIMPORT import json, types, collections

	allowed_types = [type(None), bool, int, int, float, complex,
					 str, bytes]

	local_vars_keys = list(local_vars.keys())

	my_vars = dict()
	for key in set(local_vars_keys) - set(exclude_list_of_vars):
		if type(local_vars[key]) in allowed_types:
			my_vars[key] = local_vars[key]
		elif type(local_vars[key]) in [list, tuple, type(set())]:
			if len({type(i) for i in local_vars[key]} - set(allowed_types)) == 0:
				if key in vars_that_will_show_only_size:
					my_vars[key] = "%s with length: %d"%(str(type(local_vars[key])),len(local_vars[key]))
				else:
					if	type(local_vars[key]) == type(set()):
						my_vars[key] = list(local_vars[key])
					else:
						my_vars[key] = local_vars[key]
		elif type(local_vars[key]) == dict:
			if len({type(local_vars[key][i]) for i in local_vars[key]} - set(allowed_types)) == 0:
					my_vars[key] = local_vars[key]

	ordered_my_vars = collections.OrderedDict(sorted(my_vars.items()))

	with open(filename, write_or_append) as fp:
		json.dump(ordered_my_vars, fp, indent = 2)
	fp.close()


def convert_json_fromunicode(data):
	pass#IMPORTIMPORTIMPORT import  collections
	pass#IMPORTIMPORTIMPORT import six
	if isinstance(data, six.string_types):
		return str(data)
	elif isinstance(data, collections.Mapping):
		return dict(list(map(convert_json_fromunicode, iter(list(data.items())))))
	elif isinstance(data, collections.Iterable):
		return type(data)(list(map(convert_json_fromunicode, data)))
	else:
		return data

"""Multiline Comment23"""

def get_sorting_attr_stack(data_stack):
	pass#IMPORTIMPORTIMPORT from utilities import get_params_proj
	attr_value_list = []
	for idat in range(len(data_stack)):
		group                 = data_stack[idat].get_attr("group")
		phi,theta,psi,s2x,s2y = get_params_proj(data_stack[idat],xform = "xform.projection")
		attr_value_list.append([group, phi, theta, psi, s2x, s2y])
	return attr_value_list

def get_sorting_params_refine(Tracker,data,ndata):
	pass#IMPORTIMPORTIMPORT from mpi import mpi_barrier, MPI_COMM_WORLD
	pass#IMPORTIMPORTIMPORT from utilities import read_text_row,wrap_mpi_bcast,even_angles
	pass#IMPORTIMPORTIMPORT from applications import MPI_start_end
	myid       = Tracker["constants"]["myid"]
	main_node  = Tracker["constants"]["main_node"]
	nproc      = Tracker["constants"]["nproc"]
	#ndata     = Tracker["total_stack"]
	mpi_comm   = mpi.MPI_COMM_WORLD
	if myid == main_node:
		total_attr_value_list = []
		for n in range(ndata):
			total_attr_value_list.append([])
	else:
		total_attr_value_list = 0
	for inode in range(nproc):
		attr_value_list = get_sorting_attr_stack(data)
		attr_value_list = wrap_mpi_bcast(attr_value_list,inode)
		if myid == main_node:
			image_start,image_end = sparx_applications.MPI_start_end(ndata,nproc,inode)
			total_attr_value_list = fill_in_mpi_list(total_attr_value_list, attr_value_list, image_start,image_end)
		mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
	total_attr_value_list = wrap_mpi_bcast(total_attr_value_list, main_node)
	return total_attr_value_list

def parsing_sorting_params(sorting_params_list):
	group_list        = []
	ali3d_params_list = []
	for element in sorting_params_list:
		group_list.append(element[0])
		ali3d_params_list.append(element[1:])
	return group_list, ali3d_params_list

def fill_in_mpi_list(mpi_list,data_list,index_start,index_end):
	for index in range(index_start, index_end):
		mpi_list[index] = data_list[index-index_start]
	return mpi_list

def sample_down_1D_curve(nxinit, nnxo, pspcurv_nnxo_file):
	shrinkage=float(nnxo)/float(nxinit)
	curv_orgn = read_text_file(pspcurv_nnxo_file)
	new_curv=int(1.5*len(curv_orgn))*[0.0]
	for index in range(len(curv_orgn)):
		new_index = int(index/shrinkage)
		fraction  = index/shrinkage-new_index
		if fraction <=0:
			new_curv[new_index] +=curv_orgn[index]
		else:
			new_curv[new_index]  +=(1.-fraction)*curv_orgn[index]
			new_curv[new_index+1] += fraction*curv_orgn[index]
	return new_curv

def get_initial_ID(part_list, full_ID_dict):
	part_initial_id_list = []
	#new_dict = {}
	for iptl in range(len(part_list)):
		part_initial_id_list.append(full_ID_dict[part_list[iptl]])
		#new_dict[iptl] = id
	return part_initial_id_list#, new_dict

def print_upper_triangular_matrix(data_table_dict,N_indep,log_main):
		msg =""
		for i in range(N_indep):
			msg +="%7d"%i
		log_main.add(msg)
		for i in range(N_indep):
			msg ="%5d "%i
			for j in range(N_indep):
				if i<j:
					msg +="%5.2f "%data_table_dict[(i,j)]
				else:
					msg +="      "
			log_main.add(msg)

def convertasi(asig,K):
	pass#IMPORTIMPORTIMPORT from numpy import array
	p = []
	for k in range(K):
		l = []
		for i in range(len(asig)):
			if( asig[i ]== k ): l.append(i)
		l = numpy.array(l,"int32")
		l.sort()
		p.append(l)
	return p

def prepare_ptp(data_list, K):
	num_of_pt = len(data_list)
	ptp=[]
	for ipt in range(num_of_pt):
		ptp.append([])
	for ipt in range(num_of_pt):
		nc = len(data_list[ipt])
		asig  =[-1]*nc
		for i in range(nc):
			asig[i] = data_list[ipt][i]
		ptp[ipt] = convertasi(asig, K)
	return ptp

def print_dict(dict,theme):
		line = time.strftime("%Y-%m-%d_%H:%M:%S", time.localtime()) + " =>"
		print((line+theme))
		spaces = "                           "
		for key, value in sorted( dict.items() ):
			if(key != "constants"):
				print(("                    => "+key+spaces[len(key):]+":  "+str(value)))

def get_resolution_mrk01(vol, radi, nnxo, fscoutputdir, mask_option):
	# this function is single processor
	#  Get updated FSC curves, user can also provide a mask using radi variable
	pass#IMPORTIMPORTIMPORT import types
	pass#IMPORTIMPORTIMPORT from statistics import fsc
	pass#IMPORTIMPORTIMPORT from utilities import model_circle, get_im
	pass#IMPORTIMPORTIMPORT from filter import fit_tanh1
	pass#IMPORTIMPORTIMPORT import os
	if(type(radi) == int):
		if(mask_option is None):  mask = model_circle(radi,nnxo,nnxo,nnxo)
		else:                     mask = get_im(mask_option)
	else:  mask = radi
	nfsc = sparx_statistics.fsc(vol[0]*mask,vol[1]*mask, 1.0,os.path.join(fscoutputdir,"fsc.txt") )
	currentres = -1.0
	ns = len(nfsc[1])
	#  This is actual resolution, as computed by 2*f/(1+f)
	for i in range(1,ns-1):
		if ( nfsc[1][i] < 0.333333333333333333333333):
			currentres = nfsc[0][i-1]
			break
		#if(currentres < 0.0):
			#print("  Something wrong with the resolution, cannot continue")
		currentres = nfsc[0][i-1]

	"""Multiline Comment24"""
	lowpass, falloff = sparx_filter.fit_tanh1(nfsc, 0.01)
	return  round(lowpass,4), round(falloff,4), round(currentres,2)

def partition_to_groups(alist, K):
	res =[]
	for igroup in range(K):
		this_group =[]
		for imeb in range(len(alist)):
			if( alist[imeb] == igroup ):   this_group.append(imeb)
		this_group.sort()
		res.append(this_group)
	return res

def partition_independent_runs(run_list, K):
	indep_runs_groups = {}
	for indep in range(len(run_list)):
		indep_runs_groups[indep] = partition_to_groups(run_list[indep], K)
	return indep_runs_groups

def merge_groups(stable_members_list):
	alist=[]
	for i in range(len(stable_members_list)):
		for j in range(len(stable_members_list[i])):alist.append(stable_members_list[i][j])
	return alist

def save_alist(Tracker,name_of_the_text_file,alist):
	pass#IMPORTIMPORTIMPORT from utilities import write_text_file
	pass#IMPORTIMPORTIMPORT import os
	log       =Tracker["constants"]["log_main"]
	myid      =Tracker["constants"]["myid"]
	main_node =Tracker["constants"]["main_node"]
	dir_to_save_list =Tracker["this_dir"]
	if myid==main_node:
		file_name=os.path.join(dir_to_save_list,name_of_the_text_file)
		write_text_file(alist, file_name)

def margin_of_error(P, size_of_this_sampling):
	# margin of an error, or radius of an error for a percentage
	pass#IMPORTIMPORTIMPORT from math import sqrt
	return numpy.sqrt(P*(1.-P)/size_of_this_sampling)

def do_two_way_comparison(Tracker):
	pass#IMPORTIMPORTIMPORT from mpi import mpi_barrier, MPI_COMM_WORLD
	pass#IMPORTIMPORTIMPORT from utilities import read_text_file,write_text_file
	pass#IMPORTIMPORTIMPORT from statistics import k_means_match_clusters_asg_new
	pass#IMPORTIMPORTIMPORT import os
	######
	myid              = Tracker["constants"]["myid"]
	main_node         = Tracker["constants"]["main_node"]
	log_main          = Tracker["constants"]["log_main"]
	total_stack       = Tracker["this_total_stack"]
	workdir           = Tracker["this_dir"]
	number_of_groups  = Tracker["number_of_groups"]
	######
	if myid ==main_node:
		msg="-------Two_way comparisons analysis of %3d independent runs of equal Kmeans-------"%Tracker["constants"]["indep_runs"]
		log_main.add(msg)
	total_partition = []
	if( Tracker["constants"]["indep_runs"]<2 ):
		if myid ==main_node:
			log_main.add(" Error! One single run cannot make two-way comparison")
		pass#IMPORTIMPORTIMPORT from mpi import mpi_finalize
		pass#IMPORTIMPORTIMPORT from sys import exit
		mpi.mpi_finalize()
		exit()
	else:
		for iter_indep in range(Tracker["constants"]["indep_runs"]):  total_partition.append(Tracker["partition_dict"][iter_indep])
		### Two-way comparision is carried out on all nodes
		ptp = prepare_ptp(total_partition, number_of_groups)
		indep_runs_to_groups = partition_independent_runs(total_partition, number_of_groups)
		###### Check margin of error
		if myid ==main_node:
			log_main.add("--------------------------margin of error--------------------------------------------")
		for indep in range(len(indep_runs_to_groups)):
			for index_of_class in range(len(indep_runs_to_groups[indep])):
				one_group_in_old_ID = get_initial_ID(indep_runs_to_groups[indep][index_of_class], Tracker["full_ID_dict"])
				rate1, rate2, size_of_this_group = count_chunk_members(Tracker["chunk_dict"], one_group_in_old_ID)
				error = margin_of_error(Tracker["P_chunk0"], size_of_this_group)
				if myid ==main_node:
					log_main.add(" margin of error for chunk0 is %f   %f    %d"%((Tracker["P_chunk0"]-error),(Tracker["P_chunk0"]+error),size_of_this_group))
					log_main.add(" actual percentage is %f"%rate1)
		if myid ==main_node:
			log_main.add("------------------------------------------------------------------------------")
		total_pop=0
		two_ways_stable_member_list = {}
		avg_two_ways                = 0.0
		avg_two_ways_square         = 0.0
		scores                      = {}
		for iptp in range(len(ptp)):
			for jptp in range(len(ptp)):
				newindeces, list_stable, nb_tot_objs = sparx_statistics.k_means_match_clusters_asg_new(ptp[iptp], ptp[jptp])
				tt = 0.0
				if myid ==main_node and iptp<jptp:
					aline="Two-way comparison between independent run %3d and %3d"%(iptp,jptp)
					log_main.add(aline)
				for m in range(len(list_stable)):
					tt +=len(list_stable[m])
				if( (myid == main_node) and (iptp<jptp) ):
					unaccounted = total_stack-tt
					ratio_unaccounted  = 100.-tt/total_stack*100.
					ratio_accounted    = tt/total_stack*100
				rate = tt/total_stack*100.0
				scores[(iptp,jptp)]    = rate
				if iptp<jptp :
					avg_two_ways 	    += rate
					avg_two_ways_square += rate**2
					total_pop += 1
					new_list=[]
					for any in list_stable:
						any.tolist()
						new_list.append(any)
					two_ways_stable_member_list[(iptp,jptp)] = new_list[:]
					del new_list
		if myid ==main_node:
			log_main.add("two_way comparison is done!")
		#### Score each independent run by pairwise summation
		summed_scores = []
		two_way_dict  = {}
		for ipp in range(len(ptp)):
			avg_scores =0.0
			for jpp in range(len(ptp)):
				if ipp!=jpp:
					avg_scores += scores[(ipp,jpp)]
			avg_rate =avg_scores/(len(ptp)-1)
			summed_scores.append(avg_rate)
			two_way_dict[avg_rate] =ipp
		#### Select two independent runs that have the first two highest scores
		run1, run2,rate1,rate2 = select_two_runs(summed_scores,two_way_dict)
		Tracker["two_way_stable_member"]      = two_ways_stable_member_list[(run1,run2)]
		Tracker["pop_size_of_stable_members"] = 1
		if myid == main_node:
			log_main.add("Get outliers of the selected comparison")
		####  Save both accounted ones and unaccounted ones
		if myid == main_node:
			log_main.add("Save outliers")
		stable_class_list = []
		small_group_list  = []
		if myid ==main_node:
			log_main.add("------------------margin of error--------------------------------------------")
		for istable in range(len(Tracker["two_way_stable_member"])):
			new_one_class                    = get_initial_ID(Tracker["two_way_stable_member"][istable], Tracker["full_ID_dict"])
			rate1, rate2, size_of_this_group = count_chunk_members(Tracker["chunk_dict"], new_one_class)
			error=margin_of_error(Tracker["P_chunk0"],size_of_this_group)
			if myid ==main_node:
				log_main.add(" margin of error for chunk0 is %f    %f    %d"%((Tracker["P_chunk0"]-error),(Tracker["P_chunk0"]+error),size_of_this_group))
				log_main.add(" actual percentage is %f"%rate1)
			if( len(new_one_class)>= Tracker["constants"]["smallest_group"] ):  stable_class_list.append(new_one_class)
			else:                                                               small_group_list.append(new_one_class)
		if myid ==main_node:
			log_main.add("----------------------------------------------------------------------------")
		accounted_list = merge_groups(stable_class_list)
		Tracker["this_accounted_list"]   =  accounted_list
		Tracker["two_way_stable_member"] =  stable_class_list
		mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
		save_alist(Tracker,"Accounted.txt", accounted_list)
		update_full_dict(accounted_list,Tracker)# Update full_ID_dict for Kmeans
		mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
		Tracker["this_unaccounted_dir"]     = workdir
		Tracker["this_unaccounted_text"]    = os.path.join(workdir,"Unaccounted.txt")
		Tracker["this_accounted_text"]      = os.path.join(workdir,"Accounted.txt")
		Tracker["ali3d_of_outliers"]        = os.path.join(workdir,"ali3d_params_of_outliers.txt")
		Tracker["ali3d_of_accounted"]       = os.path.join(workdir,"ali3d_params_of_accounted.txt")
		if myid==main_node:
			log_main.add(" Selected indepedent runs      %5d and  %5d"%(run1,run2))
			log_main.add(" Their pair-wise averaged rates are %5.2f  and %5.2f "%(rate1,rate2))
		pass#IMPORTIMPORTIMPORT from math import sqrt
		avg_two_ways        = avg_two_ways/total_pop
		two_ways_std        = numpy.sqrt(avg_two_ways_square/total_pop-avg_two_ways**2)
		net_rate            = avg_two_ways-1./number_of_groups*100.
		Tracker["net_rate"] = net_rate
		if myid == main_node:
			msg="average of two-way comparison  %5.3f"%avg_two_ways
			log_main.add(msg)
			msg="net rate of two-way comparison  %5.3f"%net_rate
			log_main.add(msg)
			msg="std of two-way comparison %5.3f"%two_ways_std
			log_main.add(msg)
			msg ="Score table of two_way comparison when Kgroup =  %5d"%number_of_groups
			log_main.add(msg)
			print_upper_triangular_matrix(scores,Tracker["constants"]["indep_runs"],log_main)
		del two_ways_stable_member_list
		Tracker["score_of_this_comparison"]=(avg_two_ways,two_ways_std,net_rate)
		mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

def select_two_runs(summed_scores,two_way_dict):
	summed_scores.sort()
	rate1 = summed_scores[-1]
	rate2 = None
	for index in range(2,len(summed_scores)+1):
		rate2 =summed_scores[-index]
		if rate2 !=rate1:
			break
	if rate2 != None:
		if rate1 != rate2:
			tmp_run1= two_way_dict[rate1]
			tmp_run2= two_way_dict[rate2]
			run1 = min(tmp_run1,tmp_run2)
			run2 = max(tmp_run1,tmp_run2)
		else:
			run1 = 0
			run2 = 1
			rate2=rate1
	else:
		run1 =0
		run2 =1
		rate2 = rate1
	return run1, run2, rate1, rate2

def counting_projections(delta, ali3d_params, image_start):
	pass#IMPORTIMPORTIMPORT from utilities import even_angles,angle_between_projections_directions
	sampled_directions = {}
	angles=even_angles(delta,0,180)
	for a in angles:
		[phi0, theta0, psi0]=a
		sampled_directions[(phi0,theta0)]=[]
	pass#IMPORTIMPORTIMPORT from math import sqrt
	for i in range(len(ali3d_params)):
		[phi, theta, psi, s2x, s2y] = ali3d_params[i]
		dis_min    = 9999.
		this_phi   = 9999.
		this_theta = 9999.
		this_psi   = 9999.
		prj1       =[phi,theta]
		for j in range(len(angles)):
			[phi0, theta0, psi0] = angles[j]
			prj2 =[phi0,theta0]
			dis = angle_between_projections_directions(prj1, prj2)
			if dis<dis_min:
				dis_min    =dis
				this_phi   =phi0
				this_theta =theta0
				this_psi   =psi0
		alist = sampled_directions[(this_phi,this_theta)]
		alist.append(i+image_start)
		sampled_directions[(this_phi,this_theta)]=alist
	return sampled_directions

def unload_dict(dict_angles):
	dlist =[]
	for a in dict_angles:
		tmp=[a[0],a[1]]
		tmp_list=dict_angles[a]
		for b in tmp_list:
			tmp.append(b)
		dlist.append(tmp)
	return dlist

def load_dict(dict_angle_main_node, unloaded_dict_angles):
	for ang_proj in unloaded_dict_angles:
		if len(ang_proj)>2:
			for item in range(2,len(ang_proj)):
				dict_angle_main_node[(ang_proj[0],ang_proj[1])].append(item)
	return dict_angle_main_node

def get_stat_proj(Tracker,delta,this_ali3d):
	pass#IMPORTIMPORTIMPORT from mpi import mpi_barrier, MPI_COMM_WORLD
	pass#IMPORTIMPORTIMPORT from utilities import read_text_row,wrap_mpi_bcast,even_angles
	pass#IMPORTIMPORTIMPORT from applications import MPI_start_end
	myid      = Tracker["constants"]["myid"]
	main_node = Tracker["constants"]["main_node"]
	nproc     = Tracker["constants"]["nproc"]
	mpi_comm  = mpi.MPI_COMM_WORLD
	if myid ==main_node:
		ali3d_params=read_text_row(this_ali3d)
		lpartids    = list(range(len(ali3d_params)))
	else:
		lpartids      = 0
		ali3d_params  = 0
	lpartids = wrap_mpi_bcast(lpartids, main_node)
	ali3d_params = wrap_mpi_bcast(ali3d_params, main_node)
	ndata=len(ali3d_params)
	image_start, image_end = sparx_applications.MPI_start_end(ndata, nproc, myid)
	ali3d_params=ali3d_params[image_start:image_end]
	sampled=counting_projections(delta,ali3d_params,image_start)
	for inode in range(nproc):
		if myid ==inode:
			dlist=unload_dict(sampled)
		else:
			dlist =0
		dlist=wrap_mpi_bcast(dlist,inode)
		if myid ==main_node and inode != main_node:
			sampled=load_dict(sampled,dlist)
		mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
	return sampled

def create_random_list(Tracker):
	pass#IMPORTIMPORTIMPORT import copy
	pass#IMPORTIMPORTIMPORT import random
	pass#IMPORTIMPORTIMPORT from utilities import wrap_mpi_bcast

	myid        = Tracker["constants"]["myid"]
	main_node   = Tracker["constants"]["main_node"]
	total_stack = Tracker["total_stack"]

	if Tracker["constants"]["seed"] ==- 1: random.seed()
	else:                                  random.seed(Tracker["constants"]["seed"])

	indep_list  = []
	for irandom in range(Tracker["constants"]["indep_runs"]):
		ll = copy.copy(Tracker["this_data_list"])
		random.shuffle(ll)
		ll = wrap_mpi_bcast(ll, main_node)
		indep_list.append(ll)
	Tracker["this_indep_list"] = indep_list

def recons_mref(Tracker):
	pass#IMPORTIMPORTIMPORT from mpi import mpi_barrier, MPI_COMM_WORLD
	pass#IMPORTIMPORTIMPORT import os
	pass#IMPORTIMPORTIMPORT from time import sleep
	pass#IMPORTIMPORTIMPORT from reconstruction import recons3d_4nn_ctf_MPI
	pass#IMPORTIMPORTIMPORT from utilities import get_shrink_data_huang
	myid             = Tracker["constants"]["myid"]
	main_node        = Tracker["constants"]["main_node"]
	nproc            = Tracker["constants"]["nproc"]
	number_of_groups = Tracker["number_of_groups"]
	particle_list    = Tracker["this_particle_list"]
	nxinit           = Tracker["nxinit"]
	partstack        = Tracker["constants"]["partstack"]
	total_data       = len(particle_list)
	ref_list = []
	number_of_ref_class = []
	for igrp in range(number_of_groups):
		a_group_list = particle_list[(total_data*igrp)//number_of_groups:(total_data*(igrp+1))//number_of_groups]
		a_group_list.sort()
		Tracker["this_data_list"] = a_group_list
		pass#IMPORTIMPORTIMPORT from utilities import write_text_file
		particle_list_file = os.path.join(Tracker["this_dir"], "iclass%d.txt"%igrp)
		if myid ==main_node:
			write_text_file(Tracker["this_data_list"],particle_list_file)
		mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
		data, old_shifts =  get_shrink_data_huang(Tracker,nxinit,particle_list_file,partstack,myid,main_node,nproc,preshift=True)
		#vol=reconstruct_3D(Tracker,data)
		mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
		vol = sparx_reconstruction.recons3d_4nn_ctf_MPI(myid=myid,prjlist=data,symmetry=Tracker["constants"]["sym"],finfo=None)
		if myid ==main_node:
			print("reconstructed %3d"%igrp)
		ref_list.append(vol)
		number_of_ref_class.append(len(Tracker["this_data_list"]))
	Tracker["number_of_ref_class"] = number_of_ref_class
	return ref_list

def apply_low_pass_filter(refvol,Tracker):
	pass#IMPORTIMPORTIMPORT from filter import filt_tanl
	for iref in range(len(refvol)):
		refvol[iref]=sparx_filter.filt_tanl(refvol[iref],Tracker["low_pass_filter"],.1)
	return refvol

def get_groups_from_partition(partition, initial_ID_list, number_of_groups):
	# sort out Kmref results to individual groups that has initial IDs
	# make a dictionary
	dict = {}
	for iptl in range(len(initial_ID_list)):
		dict[iptl] = initial_ID_list[iptl]
	res = []
	for igrp in range(number_of_groups):
		class_one = []
		for ipt in range(len(partition)):
			if partition[ipt] == igrp:
				orginal_id = dict[ipt]
				class_one.append(orginal_id)
		res.append(class_one)
	return res

def get_complementary_elements(total_list,sub_data_list):
	if len(total_list)<len(sub_data_list):
		print("Wrong input list!")
		return []
	else:
		sub_data_dict     = {}
		complementary     = []
		for index in range(len(sub_data_list)):sub_data_dict[sub_data_list[index]]=index
		for any in total_list:
			if (any in sub_data_dict) is False:complementary.append(any)
		return complementary

def update_full_dict(leftover_list, Tracker):
	full_dict = {}
	for iptl in range(len(leftover_list)):
		full_dict[iptl]     = leftover_list[iptl]
	Tracker["full_ID_dict"] = full_dict

def count_chunk_members(chunk_dict, one_class):
	N_chunk0 = 0
	N_chunk1 = 0
	for a in one_class:
		if chunk_dict[a] == 0: N_chunk0 += 1
		else:                  N_chunk1 += 1
	n = len(one_class)
	if n <=1:  return 0.0, 0.0, n
	else: return  float(N_chunk0)/n, float(N_chunk1)/n, n

def remove_small_groups(class_list,minimum_number_of_objects_in_a_group):
	new_class  = []
	final_list = []
	for one_class in class_list:
		if len(one_class)>=minimum_number_of_objects_in_a_group:
			new_class.append(one_class)
			for element in one_class:
				final_list.append(element)
	final_list.sort()
	return final_list, new_class

def get_number_of_groups(total_particles,number_of_images_per_group):
	#minimum_number_of_members = 1000
	number_of_groups   =  float(total_particles)/number_of_images_per_group
	if number_of_groups - int(number_of_groups)<.4:number_of_groups = int(number_of_groups)
	else:number_of_groups = int(number_of_groups)+1
	return number_of_groups

def angular_distribution(inputfile, options, output):
	sparx_global_def.ERROR("Code disabled, please use sxplot_projs_distrib.py instead","angular_distribution",1)
	"""
	pass#IMPORTIMPORTIMPORT import numpy

	#print('Loading data')
	# Import data
	listDType = [
		('Phi', '<f8'),
		('Theta', '<f8'),
	]
	#arrData = numpy.genfromtxt(inputfile, dtype=listDType, usecols=(0, 1))
	# The following two lines are in case there is no symmetry option in sxprocess, to be simiplified
	try: sym = options.symmetry
	except: sym = "c1"
	pass#IMPORTIMPORTIMPORT from fundamentals import symclass
	pass#IMPORTIMPORTIMPORT from utilities import read_text_row
	if( sym == "c0" ):
		angs = read_text_row(inputfile)
	else:
		scs = symclass(sym)
		angs = scs.reduce_anglesets(read_text_row(inputfile),0)
		del scs
	nang = len(angs)

	# Load angle Data
	# below makes no sense
	#arrPhi = numpy.round(arrData['Phi'], options.round_digit)
	#arrTheta = numpy.round(arrData['Theta'], options.round_digit)

	arrPhi   = numpy.array([numpy.round(angs[i][0], options.round_digit) for i in xrange(nang)])
	arrTheta = numpy.array([numpy.round(angs[i][1], options.round_digit) for i in xrange(nang)])
	del angs

	# Set the vectors for transformation and plotting
	vectorInital = numpy.array([0, 0, 1])
	vectorCenter = 0.5 * numpy.array([
		options.box_size,
		options.box_size,
		options.box_size
	])

	#print('Calculate vector length')
	# Create array for the angles
	dtype = [
		('alpha', '<f8'),
		('beta', '<f8')
	]
	arrayAngles = numpy.empty(nang, dtype=dtype)
	arrayAngles['alpha'] = numpy.radians(arrTheta)
	arrayAngles['beta'] = numpy.radians(arrPhi)

	# Create length of the vectors. One angstrom is one particle.
	uniqueArray, allArray = numpy.unique(
		arrayAngles, return_inverse=True
	)
	arrayRadius = numpy.histogram(allArray, bins=len(uniqueArray))[0]

	# Calculate the overall number of particles for the normalisation.
	# Normalise the radius and calculate
	# how many times there is the same radius.
	particleNumber = len(arrayAngles)
	arrayRadius = arrayRadius / float(particleNumber)
	uniqueRadius, indicesRadius = numpy.unique(
		arrayRadius, return_index=True
	)

	# Set the right colour to the right radius
	uniqueRadiusNumber = len(uniqueRadius)
	rangeGreen = numpy.linspace(0, 1, uniqueRadiusNumber)
	rangeBlue = numpy.linspace(1, 0, uniqueRadiusNumber)

	sortRadius = numpy.sort(uniqueRadius)
	dictColor = {}
	for number, radius in enumerate(sortRadius):
		dictColor.update(
			{
				radius:
				str(rangeGreen[number]) +
				' ' +
				str(rangeBlue[number])
			}
		)

	# Merge all unique data and the related radius into one array
	dtype = [
		('alpha', '<f8'),
		('beta', '<f8'),
		('radius', '<f8')
	]
	arrayAnglesRadius = numpy.empty(len(uniqueArray['alpha']), dtype=dtype)
	arrayAnglesRadius['alpha'] = uniqueArray['alpha']
	arrayAnglesRadius['beta'] = uniqueArray['beta']
	arrayAnglesRadius['radius'] = arrayRadius

	#print('Write output')
	# Create vectors for chimera
	with open(output, 'w') as f:
		for vector in arrayAnglesRadius:
			arrayVector1 = numpy.empty(3)
			arrayVector2 = numpy.empty(3)
			arrayVectorSphere = numpy.empty(3)

			arrayVectorSphere[0] = numpy.sin(vector[0]) * numpy.cos(vector[1])
			arrayVectorSphere[1] = numpy.sin(vector[0]) * numpy.sin(vector[1])
			arrayVectorSphere[2] = numpy.cos(vector[0])

			arrayVector1 = vectorCenter
			arrayVector2 = vectorCenter

			arrayVector1 = arrayVector1 + \
				options.particle_radius * arrayVectorSphere / options.pixel_size
			arrayVector2 = arrayVector2 + \
				(
					options.particle_radius / options.pixel_size +
					0.01 + vector[2] * options.cylinder_length
				) * \
				arrayVectorSphere
			f.write('.color 0 {:s} \n'.format(dictColor[vector[2]]))
			f.write(
				'.cylinder {:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {:.3f} \n'.format(
					arrayVector1[0],
					arrayVector1[1],
					arrayVector1[2],
					arrayVector2[0],
					arrayVector2[1],
					arrayVector2[2],
					options.cylinder_width
				)
			)
	print(('All done! Saved output to: {0}'.format(output)))
	"""

#####---------------------------------------------------
# used in new meridien
def tabessel(nx, nnxo, nbel = 5000):
	beltab = [0.0]*nbel
	radius = 1.9
	alpha = 15
	#order = 0
	normk = EMAN2_cppwrap.Util.bessel0(0., radius, alpha)
	for i in range(nbel):
		rr = i/float(nbel-1)/2.0
		beltab[i] = EMAN2_cppwrap.Util.bessel0(rr, radius, alpha)/normk
	return beltab

####

