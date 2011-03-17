#
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

from global_def import *

def params_2D_3D(alpha, sx, sy, mirror):
	"""
		Convert 2D alignment parameters (alpha, sx, sy, mirror) into
		3D alignment parameters (phi, theta, psi, s2x, s2y, mirror)
	"""
	phi = 0
	psi = 0
	theta = 0
	alphan, s2x, s2y, scalen = compose_transform2(0, sx, sy, 1, -alpha, 0, 0, 1)
	if mirror > 0:
		phi   = (540.0 + phi)%360.0
		theta = 180.0  - theta
		psi   = (540.0 - psi + alphan)%360.0
	else:
		psi   = (psi   + alphan)%360.0
	return phi, theta, psi, s2x, s2y

def params_3D_2D(phi, theta, psi, s2x, s2y):
	"""
		Convert 3D alignment parameters (phi, theta, psi, s2x, s2y)  # there is no mirror in 3D! 
		into 2D alignment parameters (alpha, sx, sy, mirror)
	"""
	if theta > 90.0:
		mirror = 1
		alpha, sx, sy, scalen = compose_transform2(0, s2x, s2y, 1.0, 540.0-psi, 0, 0, 1.0)
	else:
		mirror = 0
		alpha, sx, sy, scalen = compose_transform2(0, s2x, s2y, 1.0, 360.0-psi, 0, 0, 1.0)
	return  alpha, sx, sy, mirror

# Amoeba uses the simplex method of Nelder and Mead to maximize a
# function of 1 or more variables.
#
#   Copyright (C) 2005  Thomas R. Metcalf
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
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
	   
	       import amoeba
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
	for i in xrange(nvar):
	    simplex[i+1] = var[:]
	    simplex[i+1][i] += scale[i]

	fvalue = []
	for i in xrange(nsimplex):  # set the function values for the simplex
	    fvalue.append(func(simplex[i],data=data))

	# Ooze the simplex to the maximum

	iteration = 0

	while 1:
	    # find the index of the best and worst vertices in the simplex
	    ssworst = 0
	    ssbest  = 0
	    for i in xrange(nsimplex):
		if fvalue[i] > fvalue[ssbest]:
		    ssbest = i
		if fvalue[i] < fvalue[ssworst]:
		    ssworst = i
		    
	    # get the average of the nsimplex-1 best vertices in the simplex
	    pavg = [0.0]*nvar
	    for i in xrange(nsimplex):
		if i != ssworst:
		    for j in range(nvar): pavg[j] += simplex[i][j]
	    for j in xrange(nvar): pavg[j] = pavg[j]/nvar # nvar is nsimplex-1
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
	    for i in xrange(nvar):
		pnew[i] = 2.0*pavg[i] - simplex[ssworst][i]
	    fnew = func(pnew,data=data)
	    if fnew <= fvalue[ssworst]:
		# the new vertex is worse than the worst so shrink
		# the simplex.
		for i in xrange(nsimplex):
		    if i != ssbest and i != ssworst:
			for j in xrange(nvar):
			    simplex[i][j] = 0.5*simplex[ssbest][j] + 0.5*simplex[i][j]
			fvalue[i] = func(simplex[i],data=data)
		for j in xrange(nvar):
		    pnew[j] = 0.5*simplex[ssbest][j] + 0.5*simplex[ssworst][j]
		fnew = func(pnew,data=data)
	    elif fnew >= fvalue[ssbest]:
		# the new vertex is better than the best so expand
		# the simplex.
		pnew2 = [0.0]*nvar
		for i in xrange(nvar):
		    pnew2[i] = 3.0*pavg[i] - 2.0*simplex[ssworst][i]
		fnew2 = func(pnew2,data=data)
		if fnew2 > fnew:
		    # accept the new vertex in the simplexe
		    pnew = pnew2
		    fnew = fnew2
	    # replace the worst vertex with the new vertex
	    for i in xrange(nvar):
		simplex[ssworst][i] = pnew[i]
	    fvalue[ssworst] = fnew
	    iteration += 1
	    #print "Iteration:",iteration,"  ",ssbest,"  ",fvalue[ssbest]

def amoeba_multi_level(var, scale, func, ftolerance=1.e-4, xtolerance=1.e-4, itmax=500, data=None):
	"""
	Commented by Zhengfan Yang on 05/01/07

	I made some change to the original amoeba so that it can now pass out some values
	calculated by func other than the criteria. This is important in multi-level
	amoeba refinement because otherwise, upper level refinement will lose the information
	of lower level refinement.
	"""
	#print " ENTER AMOEBA MULTI LEVEL"
	nvar = len(var)       # number of variables in the minimization
	nsimplex = nvar + 1   # number of vertices in the simplex

	# first set up the simplex

	simplex = [0]*(nvar+1)  # set the initial simplex
	simplex[0] = var[:]
	for i in xrange(nvar):
	    simplex[i+1] = var[:]
	    simplex[i+1][i] += scale[i]

	fvalue = []
	for i in xrange(nsimplex):  # set the function values for the simplex
	    result, passout = func(simplex[i], data=data)
	    #print  " amoeba setting ",i,simplex[i],result, passout
	    fvalue.append([result, passout])

	# Ooze the simplex to the maximum

	iteration = 0

	while 1:
	    # find the index of the best and worst vertices in the simplex
	    ssworst = 0
	    ssbest  = 0
	    for i in xrange(nsimplex):
		if fvalue[i][0] > fvalue[ssbest][0]:
		    ssbest = i
		if fvalue[i][0] < fvalue[ssworst][0]:
		    ssworst = i
		    
	    # get the average of the nsimplex-1 best vertices in the simplex
	    pavg = [0.0]*nvar
	    for i in xrange(nsimplex):
		if i != ssworst:
		    for j in range(nvar): pavg[j] += simplex[i][j]
	    for j in xrange(nvar): pavg[j] = pavg[j]/nvar # nvar is nsimplex-1
	    simscale = 0.0
	    for i in range(nvar):
		simscale += abs(pavg[i]-simplex[ssworst][i])/scale[i]
	    simscale = simscale/nvar

	    # find the range of the function values
	    fscale = (abs(fvalue[ssbest][0])+abs(fvalue[ssworst][0]))/2.0
	    if fscale != 0.0:
		frange = abs(fvalue[ssbest][0]-fvalue[ssworst][0])/fscale
	    else:
		frange = 0.0  # all the fvalues are zero in this case
		
	    # have we converged?
	    if (((ftolerance <= 0.0 or frange < ftolerance) and    # converged to maximum
		 (xtolerance <= 0.0 or simscale < xtolerance)) or  # simplex contracted enough
		(itmax and iteration >= itmax)):	     # ran out of iterations
		return simplex[ssbest],fvalue[ssbest][0],iteration,fvalue[ssbest][1]

	    # reflect the worst vertex
	    pnew = [0.0]*nvar
	    for i in xrange(nvar):
		pnew[i] = 2.0*pavg[i] - simplex[ssworst][i]
	    fnew = func(pnew,data=data)
	    if fnew[0] <= fvalue[ssworst][0]:
		# the new vertex is worse than the worst so shrink
		# the simplex.
		for i in xrange(nsimplex):
		    if i != ssbest and i != ssworst:
			for j in xrange(nvar):
			    simplex[i][j] = 0.5*simplex[ssbest][j] + 0.5*simplex[i][j]
	    		fvalue[i]  = func(simplex[i],data=data)
		for j in xrange(nvar):
		    pnew[j] = 0.5*simplex[ssbest][j] + 0.5*simplex[ssworst][j]  	    
		fnew = func(pnew, data=data)
	    elif fnew[0] >= fvalue[ssbest][0]:
		# the new vertex is better than the best so expand
		# the simplex.
		pnew2 = [0.0]*nvar
		for i in xrange(nvar):
		    pnew2[i] = 3.0*pavg[i] - 2.0*simplex[ssworst][i]
		fnew2 = func(pnew2,data=data)
		if fnew2[0] > fnew[0]:
		    # accept the new vertex in the simplexe
		    pnew = pnew2
		    fnew = fnew2
	    # replace the worst vertex with the new vertex
	    for i in xrange(nvar):
		simplex[ssworst][i] = pnew[i]
	    fvalue[ssworst] = fnew
	    iteration += 1
	    #print "Iteration:",iteration,"  ",ssbest,"  ",fvalue[ssbest]

def golden(func, args=(), brack=None, tol=1.e-4, full_output=0):
	""" Given a function of one-variable and a possible bracketing interval,
	return the minimum of the function isolated to a fractional precision of
	tol. A bracketing interval is a triple (a,b,c) where (a<b<c) and
	func(b) < func(a),func(c).  If bracket is two numbers then they are
	assumed to be a starting interval for a downhill bracket search
	(see bracket)

	Uses analog of bisection method to decrease the bracketed interval.
	"""
	if brack is None:
		xa,xb,xc,fa,fb,fc,funcalls = bracket(func, args=args)
	elif len(brack) == 2:
		xa,xb,xc,fa,fb,fc,funcalls = bracket(func, xa=brack[0], xb=brack[1], args=args)
	elif len(brack) == 3:
		xa,xb,xc = brack
		if (xa > xc):  # swap so xa < xc can be assumed
			dum = xa; xa=xc; xc=dum
		assert ((xa < xb) and (xb < xc)), "Not a bracketing interval."
		fa = apply(func, (xa,)+args)
		fb = apply(func, (xb,)+args)
		fc = apply(func, (xc,)+args)
		assert ((fb<fa) and (fb < fc)), "Not a bracketing interval."
		funcalls = 3
	else:
		raise ValueError, "Bracketing interval must be length 2 or 3 sequence."

	_gR = 0.61803399
	_gC = 1.0-_gR
	x3 = xc
	x0 = xa
	if (abs(xc-xb) > abs(xb-xa)):
		x1 = xb
		x2 = xb + _gC*(xc-xb)
	else:
		x2 = xb
		x1 = xb - _gC*(xb-xa)
	f1 = apply(func, (x1,)+args)
	f2 = apply(func, (x2,)+args)
	funcalls += 2
	while (abs(x3-x0) > tol*(abs(x1)+abs(x2))):
		if (f2 < f1):
			x0 = x1; x1 = x2; x2 = _gR*x1 + _gC*x3
			f1 = f2; f2 = apply(func, (x2,)+args)
		else:
			x3 = x2; x2 = x1; x1 = _gR*x2 + _gC*x0
			f2 = f1; f1 = apply(func, (x1,)+args)
		funcalls += 1
	if (f1 < f2):
		xmin = x1
		fval = f1
	else:
		xmin = x2
		fval = f2
	if full_output:
		return xmin, fval, funcalls
	else:
		return xmin


def bracket(func, xa=0.0, xb=1.0, args=(), grow_limit=110.0, maxiter=1000):
	"""Given a function and distinct initial points, search in the downhill
	direction (as defined by the initital points) and return new points
	xa, xb, xc that bracket the minimum of the function:
	f(xa) > f(xb) < f(xc)
	"""
	_gold = 1.618034
	_verysmall_num = 1e-21
	fa = apply(func, (xa,)+args)
	fb = apply(func, (xb,)+args)
	if (fa < fb):			   # Switch so fa > fb
	    dum = xa; xa = xb; xb = dum
	    dum = fa; fa = fb; fb = dum
	xc = xb + _gold*(xb-xa)
	fc = apply(func, (xc,)+args)
	funcalls = 3
	iter = 0
	while (fc < fb):
	    tmp1 = (xb - xa)*(fb-fc)
	    tmp2 = (xb - xc)*(fb-fa)
	    val = tmp2-tmp1
	    if abs(val) < _verysmall_num:
		denom = 2.0*_verysmall_num
	    else:
		denom = 2.0*val
	    w = xb - ((xb-xc)*tmp2-(xb-xa)*tmp1)/denom
	    wlim = xb + grow_limit*(xc-xb)
	    if iter > maxiter:
		raise RuntimeError, "Too many iterations."
	    iter += 1
	    if (w-xc)*(xb-w) > 0.0:
		fw = apply(func, (w,)+args)
		funcalls += 1
		if (fw < fc):
		    xa = xb; xb=w; fa=fb; fb=fw
		    return xa, xb, xc, fa, fb, fc, funcalls
		elif (fw > fb):
		    xc = w; fc=fw
		    return xa, xb, xc, fa, fb, fc, funcalls
		w = xc + _gold*(xc-xb)
		fw = apply(func, (w,)+args)
		funcalls += 1
	    elif (w-wlim)*(wlim-xc) >= 0.0:
		w = wlim
		fw = apply(func, (w,)+args)
		funcalls += 1
	    elif (w-wlim)*(xc-w) > 0.0:
		fw = apply(func, (w,)+args)
		funcalls += 1
		if (fw < fc):
		    xb=xc; xc=w; w=xc+_gold*(xc-xb)
		    fb=fc; fc=fw; fw=apply(func, (w,)+args)
		    funcalls += 1
	    else:
		w = xc + _gold*(xc-xb)
		fw = apply(func, (w,)+args)
		funcalls += 1
	    xa=xb; xb=xc; xc=w
	    fa=fb; fb=fc; fc=fw
	return xa, xb, xc, fa, fb, fc, funcalls

def ce_fit(inp_image, ref_image, mask_image):
	""" Fit the histogram of the input image under mask with the reference image.
		    
	     Usage : ce_fit(inp_image,ref_image,mask_image):
	    	 A and B, number of iterations and the chi-square
	"""
	hist_res = Util.histc(ref_image,inp_image,mask_image)
	args=hist_res["args"]
	scale = hist_res["scale"]
	data = [hist_res['data'],inp_image,hist_res["ref_freq_bin"],mask_image,hist_res['size_img'],hist_res['hist_len']]
	res=amoeba(args,scale,hist_func,1.e-4,1.e-4,500,data)
	resu = ["Final Parameter [A,B]:", res[0],"Final Chi-square :",-1*res[1],"Number of Iteration :",res[2]]
	corrected_image = inp_image*res[0][0] + res[0][1]
	result=[resu,"Corrected Image :",corrected_image]
	del data[:],args[:],scale[:]
	return result
	
def center_2D(image_to_be_centered, center_method = 1, searching_range = -1, Gauss_radius_inner = 2, Gauss_radius_outter = 7, self_defined_reference = None):
	"""
		Put an input image into image center (nx/2, ny/2) using method :
		1. phase_cog
		2. cross-correlate with Gaussian function
		3. cross-correlate with donut shape image
		4. cross-correlate with reference image provided by user
		5. cross-correlate with self-rotated average
	        The function will return centered_image, and shifts
	"""
	from   utilities import peak_search
	from   fundamentals import fshift
	import types
	if type(image_to_be_centered) == types.StringType: image_to_be_centered = get_im(image_to_be_centered)
	if    center_method == 0 :  return  image_to_be_centered,0.,0.
	elif  center_method == 1 :
		cs = image_to_be_centered.phase_cog()
		if searching_range > 0 :
			if(abs(cs[0]) > searching_range):  cs[0]=0.0
			if(abs(cs[1]) > searching_range):  cs[1]=0.0
		return fshift(image_to_be_centered, -cs[0], -cs[1]), cs[0], cs[1]

	elif center_method == 5:
		from fundamentals import rot_avg_image,ccf
		from math import sqrt
		not_centered = True
		tmp_image = image_to_be_centered.copy()
		shiftx = 0
		shifty = 0
		while (not_centered):
			reference = rot_avg_image(tmp_image)			
			ccmap = ccf(tmp_image, reference)
			if searching_range > 0:  ccmap = Util.window(ccmap, searching_range, searching_range, 1, 0, 0, 0)
			peak  = peak_search(ccmap)
			centered_image = fshift(tmp_image, -peak[0][4], -peak[0][5])
			if sqrt(peak[0][4]**2 + peak[0][5]**2) < 1. : not_centered = False
			else : tmp_image = centered_image.copy()
			shiftx += peak[0][4]
			shifty += peak[0][5]
		return centered_image, shiftx, shifty
		
	elif center_method == 6:
		from morphology import threshold_to_minval
		nx = image_to_be_centered.get_xsize()
		ny = image_to_be_centered.get_ysize()
		r = nx//2-2
		mask = model_circle(r, nx, ny)
		[mean, sigma, xmin, xmax] = Util.infomask(image_to_be_centered, mask, True)
		new_image = threshold_to_minval(image_to_be_centered, mean+sigma)
		cs = new_image.phase_cog()
		if searching_range > 0 :
			if(abs(cs[0]) > searching_range):  cs[0]=0.0
			if(abs(cs[1]) > searching_range):  cs[1]=0.0
		return fshift(image_to_be_centered, -cs[0], -cs[1]), cs[0], cs[1]

	else :
		nx = image_to_be_centered.get_xsize()
		ny = image_to_be_centered.get_ysize()
		from fundamentals import ccf	
		if center_method == 2 :
			reference = model_gauss(Gauss_radius_inner, nx, ny)
		if center_method == 3 :
			do1 = model_gauss(Gauss_radius_outter, nx, ny)
			do2 = model_gauss(Gauss_radius_inner,  nx, ny)
			s = Util.infomask(do1, None, True)
			do1/= s[3]
			s = Util.infomask(do2, None, True)
			do2/=s[3]
			reference = do1 - do2
		if center_method == 4:	reference = self_defined_reference
		ccmap = ccf(image_to_be_centered, reference)
		if searching_range > 1: ccmap = Util.window(ccmap, searching_range, searching_range, 1, 0, 0, 0)
		peak  = peak_search(ccmap)
		return fshift(image_to_be_centered, -peak[0][4], -peak[0][5]), peak[0][4], peak[0][5]

def common_line_in3D(phiA,thetaA,phiB,thetaB):
	"""Find the position of the commone line in 3D
           Formula is   (RB^T zhat)   cross  (RA^T zhat)
	   Returns phi, theta of the common line in degrees. theta always < 90
	   Notice you don't need to enter psi's; they are irrelevant
	"""

	from math import pi, sqrt, cos, sin, asin, atan2

	piOver=pi/180.0;
	ph1 = phiA*piOver; 
	th1 = thetaA*piOver; 
	ph2 = phiB*piOver; 
	th2 = thetaB*piOver;
	
 	#nx = cos(thetaBR)*sin(thetaAR)*sin(phiAR) - cos(thetaAR)*sin(thetaBR)*sin(phiBR) ;
	#ny = cos(thetaAR)*sin(thetaBR)*cos(phiBR) - cos(thetaBR)*sin(thetaAR)*cos(phiAR) ;
	#nz = sin(thetaAR)*sin(thetaBR)*sin(phiAR-phiBR);


	nx = sin(th1)*cos(ph1)*sin(ph2)-sin(th2)*sin(ph1)*cos(ph2)
	ny = sin(th1)*cos(th2)*cos(ph1)*cos(ph2)-cos(th1)*sin(th2)*cos(ph1)*cos(ph2)
	nz = cos(th2)*sin(ph1)*cos(ph2)-cos(th1)*cos(ph1)*sin(ph2)

	norm = nx*nx + ny*ny + nz*nz
 
	if norm < 1e-5:
		#print 'phiA,thetaA,phiB,thetaB:', phiA, thetaA, phiB, thetaB
		return 0.0, 0.0

	if nz<0: nx=-nx; ny=-ny; nz=-nz;

	#thetaCom = asin(nz/sqrt(norm))
	phiCom    = asin(nz/sqrt(norm))
	#phiCom   = atan2(ny,nx)
	thetaCom  = atan2(ny, nx)
	
	return phiCom*180.0/pi , thetaCom*180.0/pi

def compose_transform2(alpha1, sx1, sy1, scale1, alpha2, sx2, sy2, scale2):
	"""Print the composition of two transformations  T2*T1
		Here  if v's are vectors:   vnew = T2*T1 vold
		     with T1 described by alpha1, sx1, scale1 etc.

	    Usage: compose_transform2(alpha1,sx1,sy1,scale1,alpha2,sx2,sy2,scale2)
	       angles in degrees
	"""
	t1 = Transform({"type":"2D","alpha":alpha1,"tx":sx1,"ty":sy1,"mirror":0,"scale":scale1})
	t2 = Transform({"type":"2D","alpha":alpha2,"tx":sx2,"ty":sy2,"mirror":0,"scale":scale2})
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
	R1 = Transform({"type":"spider","phi":float(phi1),"theta":float(theta1),"psi":float(psi1),"tx":float(sx1),"ty":float(sy1),"tz":float(sz1),"mirror":0,"scale":float(scale1)})
	R2 = Transform({"type":"spider","phi":float(phi2),"theta":float(theta2),"psi":float(psi2),"tx":float(sx2),"ty":float(sy2),"tz":float(sz2),"mirror":0,"scale":float(scale2)})
	Rcomp=R2*R1
	d = Rcomp.get_params("spider")
	return d["phi"],d["theta"],d["psi"],d["tx"],d["ty"],d["tz"],d["scale"]
   
def combine_params2(alpha1, sx1, sy1, mirror1, alpha2, sx2, sy2, mirror2):
	"""
	  Combine 2D alignent parameters including mirror
	"""
	t1 = Transform({"type":"2D","alpha":alpha1,"tx":sx1,"ty":sy1,"mirror":mirror1,"scale":1.0})
	t2 = Transform({"type":"2D","alpha":alpha2,"tx":sx2,"ty":sy2,"mirror":mirror2,"scale":1.0})
	tt = t2*t1
	d = tt.get_params("2D")
	return d[ "alpha" ], d[ "tx" ], d[ "ty" ], d[ "mirror" ]

def create_spider_doc(fname,spiderdoc):
	"""Convert a text file that is composed of columns of numbers into spider doc file
	"""
	from string import atoi,atof
	infile = open(fname,"r")
	lines  = infile.readlines()
	infile.close()
	nmc  = len(lines[0].split())
	table=[]
	for line in lines:
		data = line.split()
	for i in xrange(0,nmc):
		data[i] = atof(data[i])
		table.append(data)
	drop_spider_doc(spiderdoc ,table)

def drop_image(imagename, destination, itype="h"):
	"""Write an image to the disk.

	Usage:  drop_image(name_of_existing_image, "path/to/image",
			  type = <type>)
	<type> is "h" (hdf) or "s" (spider)
	"""

	if type(destination) == type(""):
		if(itype == "h"):    imgtype = EMUtil.ImageType.IMAGE_HDF
		elif(itype == "s"):  imgtype = EMUtil.ImageType.IMAGE_SINGLE_SPIDER
		else:  ERROR("unknown image type","drop_image",1)
		imagename.write_image(destination, 0, imgtype)
	else:
		ERROR("destination is not a file name","drop_image",1)

def drop_png_image(im, trg):
	"""Write an image with the proper png save
	Usage: drop_png_image(name_of_existing_image, 'path/to/image.png')
	"""

	if trg[-4:] != '.png':
		ERROR('destination name must be png extension', 'drop_png_image', 1)

	if isinstance(trg, basestring):
		im['render_min'] = im['minimum']
		im['render_max'] = im['maximum']
		im.write_image(trg, 0)
	else:
		ERROR('destination is not a file name', 'drop_png_image', 1)

def drop_spider_doc(filename, data, comment = None):
	"""Create a spider-compatible "Doc" file.

	   filename: name of the Doc file
	   data: List of lists, with the inner list being a list of floats
	         and is written as a line into the doc file.
	"""
	outf = open(filename, "w")
	from datetime import datetime
	outf.write(" ;   %s   %s   %s\n" % (datetime.now().ctime(), filename, comment))
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

def dump_row(input, fname, ix=0, iz=0):
	"""Output the data in slice iz, row ix of an image to standard out.

	Usage: dump_row(image, ix, iz)
	   or
	       dump_row("path/to/image", ix, iz)
	"""
	fout = open(fname, "w")
	image=get_image(input)
	nx = image.get_xsize()
	ny = image.get_ysize()
	nz = image.get_zsize()
	fout.write("# z = %d slice, x = %d row)\n" % (iz, ix))
	line = []
	for iy in xrange(ny):
		fout.write("%d\t%12.5g\n" % (iy, image.get_value_at(ix,iy,iz)))
	fout.close()

def even_angles(delta = 15.0, theta1=0.0, theta2=90.0, phi1=0.0, phi2=359.99, method = 'S', phiEqpsi = "Minus", symmetry='c1'):
	"""Create a list of Euler angles suitable for projections.
	   method is either 'S' - for Saff algorithm
	                  or   'P' - for Penczek '94 algorithm
			  'S' assumes phi1<phi2 and phi2-phi1>> delta ;
	   symmetry  - if this is set to c3, will yield angles from the asymmetric unit, not the specified range;
	"""

	from math      import pi, sqrt, cos, acos, tan, sin
	from utilities import even_angles_cd
	from string    import lower,split
	angles = []
	symmetryLower = symmetry.lower()
	symmetry_string = split(symmetry)[0]
	if  (symmetry_string[0]  == "c"):
		if(phi2 == 359.99):
			angles = even_angles_cd(delta, theta1, theta2, phi1, phi2/int(symmetry_string[1:]), method, phiEqpsi)
			if(int(symmetry_string[1:]) > 1):
				qt = 180.0/int(symmetry_string[1:])
				n = len(angles)
				for i in xrange(n):
					t = n-i-1
					if(angles[t][1] == 90.0):
						if(angles[t][0] >= qt):  del angles[t]
		else:
			angles = even_angles_cd(delta, theta1, theta2, phi1, phi2, method, phiEqpsi)
	elif(symmetry_string[0]  == "d"):
		if(phi2 == 359.99):
			angles = even_angles_cd(delta, theta1, theta2, phi1, 360.0/2/int(symmetry_string[1:]), method, phiEqpsi)
			qt = 180.0/2/int(symmetry_string[1:])
			n = len(angles)
			for i in xrange(n):
				t = n-i-1
				if(angles[t][1] == 90.0):
					if(angles[t][0] >= qt):  del angles[t]
		else:
			angles = even_angles_cd(delta, theta1, theta2, phi1, phi2, method, phiEqpsi)
	elif(symmetry_string[0]  == "s"):
	#if symetry is "s", deltphi=delta, theata intial=theta1, theta end=90, delttheta=theta2
				
		theta_number = int((90.0 - theta1)/theta2)
		for j in xrange(theta_number,-1, -1):
			if( j == 0):
				k=int( ((phi2/2) - phi1)/delta ) + 1
			else:
				k=int( ((phi2) - phi1)/delta ) + 1 
			for i in xrange(k):
					t = phi1 +i*delta
					angles.append([t,90.0-j*theta2,90.0])
	
	else : # This is very close to the Saff even_angles routine on the asymmetric unit;
		# the only parameters used are symmetry and delta
		# The formulae are given in the Transform Class Paper
		# The symmetric unit 		nVec=[]; # x,y,z triples
		# is defined by three points b,c, v of Fig 2 of the paper
		# b is (0,0,1)
		# c is (sin(thetac),0,cos(thetac))
		# a is (sin(thetac)cos(Omega),sin(thetac)cos(Omega),cos(thetac))
		# f is the normalized sum of all 3
		
		# The possible symmetries are in list_syms
		# The symmetry determines thetac and Omega
		# The spherical area is Omega - pi/3; 
		#  should be equal to 4 *pi/(3*# Faces)
		#		
		# symmetry ='tet';   delta  = 6;

		scrunch = 0.9  # closeness factor to eliminate oversampling corners
		#nVec=[]       # x,y,z triples

		piOver = pi/180.0
		Count=0   # used to count the number of angles
		
		if (symmetryLower[0:3] =="tet"):  m=3.0; fudge=0.9 # fudge is a factor used to adjust phi steps
		elif (symmetryLower[0:3] =="oct"):  m=4.0; fudge=0.8
		elif (symmetryLower[0:3] =="ico"):  m=5.0; fudge=0.95
		else: ERROR("allowable symmetries are cn, dn, tet, oct, icos","even_angles",1)

		n=3.0
		OmegaR = 2.0*pi/m; cosOmega= cos(OmegaR)
		Edges  = 2.0*m*n/(2.0*(m+n)-m*n)
		Faces  = 2*Edges/n
		Area   = 4*pi/Faces/3.0; # also equals  2*pi/3 + Omega
		costhetac = cosOmega/(1-cosOmega)
		deltaRad= delta*pi/180
		NumPoints = int(Area/(deltaRad*deltaRad))
		fheight = 1/sqrt(3)/ (tan(OmegaR/2.0))

		z0      = costhetac  # initialize loop	
		z       = z0
		phi     = 0
		Deltaz  = (1-costhetac)/(NumPoints-1)

		#[1, phi,180.0*acos(z)/pi,0.]
		anglesLast = [phi,180.0*acos(z)/pi,0.]
		angles.append(anglesLast)
		nLast=  [ sin(acos(z))*cos(phi*piOver) ,  sin(acos(z))*sin(phi*piOver) , z]
		nVec = []
		nVec.append(nLast)

		Count +=1

		for k in xrange(1,(NumPoints-1)):
			z=z0 + Deltaz*k  # Is it higher than fhat or lower
			r= sqrt(1-z*z)
			if (z > fheight): phiRmax= OmegaR/2.0
			if (z<= fheight):
				thetaR   = acos(z); 
				cosStuff = (cos(thetaR)/sin(thetaR))*sqrt(1. - 2 *cosOmega);
				phiMax   =  180.0*( OmegaR - acos(cosStuff))/pi
			angleJump = fudge* delta/r
			phi = (phi + angleJump)%(phiMax)
			anglesNew = [phi,180.0*acos(z)/pi,0.];
			nNew = [ sin(acos(z))*cos(phi*piOver) ,  sin(acos(z))*sin(phi*piOver) , z]
			diffangleVec = [acos(nNew[0]*nVec[k][0] +   nNew[1]*nVec[k][1] +    nNew[2]*nVec[k][2] ) for k in xrange(Count)] 
			diffMin = min(diffangleVec)
			if (diffMin>angleJump*piOver *scrunch):
				Count +=1
				angles.append(anglesNew)
				nVec.append(nNew)
				#[Count, phi,180*acos(z)/pi,0.]
			anglesLast = anglesNew
			nLast=nNew

		angles.append( [0.0, 0.0, 0.0] )
		nLast=  [ 0., 0. , 1.]
		nVec.append(nLast)
		if(theta2 == 180.0):   angles.append( [0.0, 180.0, 0.0] )
		
		angles.reverse()
		if(phiEqpsi == "Minus"):
			for i in xrange(len(angles)):  angles[i][2] = (720.0-angles[i][0])%360.0
		#print(Count,NumPoints)
		
#		look at the distribution
#		Count =len(angles); piOver= pi/180.0;
#		phiVec    =  [ angles[k][0] for k in range(Count)] ;
#		thetaVec  =  [ angles[k][1] for k in range(Count)] ;
#		xVec = [sin(piOver * angles[k][1]) * cos(piOver * angles[k][0]) for k in range(Count) ]
#		yVec = [sin(piOver * angles[k][1])* sin(piOver * angles[k][0]) for k in range(Count) ]
#		zVec = [cos(piOver *  angles[k][1]) for k in range(Count) ]
#		pylab.plot(yVec,zVec,'.'); pylab.show()


	return angles

def even_angles_cd(delta, theta1=0.0, theta2=90.0, phi1=0.0, phi2=359.99, method = 'P', phiEQpsi='Minus'):
	"""Create a list of Euler angles suitable for projections.
	   method is either 'S' - for Saff algorithm
	                  or   'P' - for Penczek '94 algorithm
			  'S' assumes phi1<phi2 and phi2-phi1>> delta ;
	   phiEQpsi  - set this to 'Minus', if you want psi=-phi;
	"""
	from math import pi, sqrt, cos, acos
	angles = []
	if (method == 'P'):
		temp = Util.even_angles(delta, theta1, theta2, phi1, phi2)
		#		                                              phi, theta, psi
		for i in xrange(len(temp)/3): angles.append([temp[3*i],temp[3*i+1],temp[3*i+2]]);
	else:              #elif (method == 'S'):
		Deltaz  = cos(theta2*pi/180.0)-cos(theta1*pi/180.0)
		s       = delta*pi/180.0
		NFactor = 3.6/s
		wedgeFactor = abs(Deltaz*(phi2-phi1)/720.0)
		NumPoints   = int(NFactor*NFactor*wedgeFactor)
		angles.append([phi1, theta1, 0.0])
		z1 = cos(theta1*pi/180.0); 	phi=phi1            # initialize loop
		for k in xrange(1,(NumPoints-1)):
			z=z1 + Deltaz*k/(NumPoints-1)
			r= sqrt(1-z*z)
			phi = phi1+(phi + delta/r -phi1)%(abs(phi2-phi1))
			#[k, phi,180*acos(z)/pi, 0]
			angles.append([phi, 180*acos(z)/pi, 0.0])
		#angles.append([p2,t2,0])  # This is incorrect, as the last angle is really the border, not the element we need. PAP 01/15/07
	if (phiEQpsi == 'Minus'):
		for k in xrange(len(angles)): angles[k][2] = (720.0 - angles[k][0])%360.0
	if( theta2 == 180.0 ):  angles.append( [0.0, 180.0, 0.0] )

	return angles
	
def eigen_images_get(stack, eigenstack, mask, num, avg):
	"""
		Perform PCA on stack file 
		and Get eigen images
	"""
	
	from utilities import get_image
	
	a = Analyzers.get('pca_large')
	e = EMData()
	if(avg == 1): s = EMData()
	nima = EMUtil.get_image_count(stack)
	for im in xrange(nima):
		e.read_image(stack,im)
		e *= mask
		a.insert_image(e)
		if( avg==1):
			if(im==0): s  = a
			else:      s += a
	if(avg == 1): a -= s/nima
	eigenimg = a.analyze()
	if(num>= EMUtil.get_image_count(eigenimg)):
		num=EMUtil.get_image_count(eigenimg)	
	for  i in xrange(num): eigenimg.write_image(eigenstack,i)
		
def find_inplane_to_match(phiA,thetaA,phiB,thetaB,psiA=0,psiB=0):
	"""Find the z rotation such that
	    ZA  RA is as close as possible to RB
	        this maximizes trace of ( RB^T ZA RA) = trace(ZA RA RB^T)
	"""
	#from math import pi, sqrt, cos, acos, sin

	RA   = Transform({'type': 'spider', 'phi': phiA, 'theta': thetaA, 'psi': psiA})
	RB   = Transform({'type': 'spider', 'phi': phiB, 'theta': thetaB, 'psi': psiB})
	RBT  = RB.transpose()
	RABT = RA * RBT

	RABTeuler = RABT.get_rotation('spider')
	RABTphi   = RABTeuler['phi']
	RABTtheta = RABTeuler['theta']
	RABTpsi   = RABTeuler['psi']

	#deg_to_rad = pi/180.0
	#thetaAR = thetaA*deg_to_rad
	#thetaBR = thetaB*deg_to_rad
	#phiAR   = phiA*deg_to_rad
	#phiBR   = phiB *deg_to_rad

	#d12=cos(thetaAR)*cos(thetaBR) + sin(thetaAR)*sin(thetaBR)*cos(phiAR-phiBR)
	return (-RABTpsi-RABTphi),RABTtheta #  180.0*acos(d12)/pi;

def find(vv, cmp_str, n):
	jFoundVec= [];
	for jFound in xrange(len(vv)):
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
	from utilities import model_gauss
	nz = sharp_edge_image.get_ndim()
	if(nz == 3):   kern = model_gauss(gauss_standard_dev, kernel_size , kernel_size, kernel_size)
	elif(nz == 2):  kern = model_gauss(gauss_standard_dev, kernel_size , kernel_size)
	else:          kern = model_gauss(gauss_standard_dev, kernel_size)
	aves = Util.infomask(kern, None, False)
	nx = kern.get_xsize()
	ny = kern.get_ysize()
	nz = kern.get_zsize()

	kern /= (aves[0]*nx*ny*nz)
	return  rsconvolution(sharp_edge_image, kern)

def get_image(imagename, nx = 0, ny = 1, nz = 1, im = 0):
	"""Read an image from the disk or assign existing object to the output.

	Usage: myimage = readImage("path/to/image")
	or     myimage = readImage(name_of_existing_image)
	"""
	if type(imagename) == type(""):
	    e = EMData()
	    e.read_image(imagename, im)
	elif not imagename:
	    e = EMData()
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
		e = EMData()
		e.read_image(stackname, im)
		return  e
	else:
		return  stackname[im].copy()

def get_image_data(img):
    """
    Return a NumPY array containing the image data. 
    Note: The NumPY array and the image data share the same memory,
          so if the NumPY array is altered then the image is altered
          as well (and vice versa).
    """
    return EMNumPy.em2numpy(img)

def get_inplane_angle(ima,ref, iring=1, fring=-1, ringstep=1, xtransSearch=0, ytransSearch=0, stp=1, center=1):
	""" 
		Get the in_plane angle from two images
		and output the crosss correlation value
		The function won't destroy input two images
		This is the angle that rotates the first image, ima, into the second image, ref.
		The sense of the rotation is clockwise.
		center=1 means image is first centered, then rotation angle is found
	"""

	from alignment import Numrinit, ringwe, Applyws, ormq
	from filter import fshift

	first_ring=int(iring); last_ring=int(fring); rstep=int(ringstep); xrng=int(xtransSearch); yrng=int(ytransSearch); step=int(stp)	
	nx=ima.get_xsize()
	if(last_ring == -1): last_ring=int(nx/2)-2
	cnx = int(nx/2)+1
 	cny = cnx
 	mode = "F"
 	#precalculate rings
	numr = Numrinit(first_ring, last_ring, rstep, mode)
 	wr = ringwe(numr, mode)
	if(center==1):
		cs = [0.0]*2  # additio
		cs = ref.phase_cog()
		ref1 = fshift(ref, -cs[0], -cs[1])
		cimage=Util.Polar2Dm(ref1, cnx, cny, numr, mode)
		cs = ima.phase_cog()
		ima1 = fshift(ima, -cs[0], -cs[1])
	else:
		ima1=ima.copy()
		cimage=Util.Polar2Dm(ref, cnx, cny, numr, mode)
	Util.Frngs(cimage, numr)
	Applyws(cimage, numr, wr)
	[angt, sxst, syst, mirrort, peakt]=ormq(ima1, cimage, xrng, yrng, step, mode, numr, cnx, cny)
	return angt,sxst, syst, mirrort, peakt

def get_sym(symmetry):
	RA   = Transform()
	NTot = RA.get_nsym(symmetry)
	angs = []
	for j in xrange(NTot):
		RNow  = RA.get_sym(symmetry, j)
		RNowE = RNow.get_rotation('spider')
		angs.append([RNowE['phi'], RNowE['theta'], RNowE['psi']])

	return angs

def get_textimage(fname):
	"""	
		Return an image created from a text file.  The first line of
		the image should contain "nx ny nz" (separated by whitespace)
		All subsequent lines contain "ix iy iz val", where ix, iy, 
		and iz are the integer x, y, and z coordinates of the point
		and val is the floating point value of that point.  All points
		not explicitly listed are set to zero.
	"""
    	from string import atoi,atof
    	infile = open(fname)
	lines = infile.readlines()
	infile.close()
	data = lines[0].split()
	nx = atoi(data[0])
	ny = atoi(data[1])
	nz = atoi(data[2])
	e = EMData()
	e.set_size(nx, ny, nz)
	e.to_zero()
	for line in lines[1:]:
		data = line.split()
		ix = atoi(data[0])
		iy = atoi(data[1])
		iz = atoi(data[2])
		val = atof(data[3])
		e[ix,iy,iz] = val
	return e

def get_input_from_string(str_input):
	"""
		Extract input numbers from given string, 
	"""
	from string import split
	res        = []
	list_input = split(str_input)
	for i in xrange(len(list_input)):
		res.append(float(list_input[i]))
	return res

def hist_func(args, data):
	#Util.hist_comp_freq(float PA,float PB,int size_img, int hist_len, float *img_ptr, float *ref_freq_bin, float *mask_ptr, float ref_h_diff, float ref_h_min)
	return Util.hist_comp_freq(args[0],args[1],data[4],data[5],data[1],data[2],data[3],data[0][0],data[0][1])
    	
def info(image, mask=None, Comment=""):
	"""Calculate and print the descriptive statistics of an image.

	Usage: [mean, sigma, xmin, xmax, nx, ny, nz =] info(image object)
	       or
	       [mean, sigma, xmin, xmax, nx, ny, nz =] info("path/image")

	Purpose: calculate basic statistical characteristics of an image.
	"""
	if(Comment):  print  " ***  ", Comment
	e = get_image(image)
	[mean, sigma, imin, imax] = Util.infomask(e, mask, True)
	nx = e.get_xsize()
	ny = e.get_ysize()
	nz = e.get_zsize()
	if (e.is_complex()):
		s = ""
		if e.is_shuffled():
			s = " (shuffled)"
		if (e.is_fftodd()):
			print "Complex odd image%s: nx = %i, ny = %i, nz = %i" % (s, nx, ny, nz)
		else:
			print "Complex even image%s: nx = %i, ny = %i, nz = %i" % (s, nx, ny, nz)

	else:
		print "Real image: nx = %i, ny = %i, nz = %i" % (nx, ny, nz)

	print "avg = %g, std dev = %g, min = %g, max = %g" % (mean, sigma, imin, imax)
	return mean, sigma, imin, imax, nx, ny, nz

def image_decimate(img,decimation=2, fit_to_fft=1,frequency_low=0, frequency_high=0):
	from filter import filt_btwl
	from fundamentals import smallprime, window2d
	from utilities import get_image
	"""
		Window image to FFT-friendly size, apply Butterworth low pass filter,
		and decimate 2D image 
	"""
	if type(img)     == str :	img=get_image(img)
	if decimation    <= 1   :  	ERROR("Improper decimation ratio","image_decimation",1)
	if frequency_low <= 0   :	
		frequency_low  = .5/decimation- .05
		if(frequency_low <= 0): ERRROR("Butterworth pass_band frequency is too low","image_decimation",1)			
		frequency_high = .5/decimation+ .05
	if fit_to_fft :
		nx_d = (img.get_xsize())/int(decimation)
		ny_d = (img.get_ysize())/int(decimation)
		nx_fft_d = smallprime(int(nx_d))
		ny_fft_d = smallprime(int(ny_d))
		nx_fft_m = nx_fft_d*int(decimation)
		ny_fft_m = ny_fft_d*int(decimation)
		e   = window2d(img, nx_fft_m, ny_fft_m, "l")
		e1  = filt_btwl(e, frequency_low, frequency_high)
		img = Util.decimate(e1, int(decimation), int(decimation), 1)
	else:
		
		e1  = filt_btwl(img, frequency_low, frequency_high)
		img = Util.decimate(e1, int(decimation), int(decimation), 1)
	return  img

def inverse_transform2(alpha, tx = 0.0, ty = 0.0, mirror = 0):
	"""Returns the inverse of the 2d rot and trans matrix

	    Usage: nalpha, ntx, nty, mirror = inverse_transform2(alpha,tx,ty,mirror)
	"""
	t = Transform({"type":"2D","alpha":alpha,"tx":tx,"ty":ty,"mirror":mirror,"scale":1.0})
	t = t.inverse()
	t = t.get_params("2D")
	return t[ "alpha" ], t[ "tx" ], t[ "ty" ], t[ "mirror" ]

def inverse_transform3(phi, theta=0.0, psi=0.0, tx=0.0, ty=0.0, tz=0.0, mirror = 0, scale=1.0):
	"""Returns the inverse of the 3d rot and trans matrix

	    Usage: nphi,ntheta,npsi,ntx,nty,ntz,nmirror,nscale = inverse_transform3(phi,theta,psi,tx,ty,tz,mirror,scale)
	       angles in degrees
	"""
	d = Transform({'type': 'spider', 'phi': phi, 'theta': theta, 'psi': psi, 'tx': tx, 'ty': ty, 'tz': tz, "mirror":mirror,"scale":scale})
	d = d.inverse()
	d = d.get_params("spider")
	return  d["phi"],d["theta"],d["psi"],d["tx"],d["ty"],d["tz"],d["mirror"],d["scale"]

def list_syms():
	"""Create a list of available symmetries
	"""
	SymStringVec=[];
	SymStringVec.append("CSYM");
	SymStringVec.append("DSYM");
	SymStringVec.append("TET_SYM");
	SymStringVec.append("OCT_SYM");
	SymStringVec.append("ICOS_SYM");
	SymStringVec.append("ISYM");
	return SymStringVec

#### -----M--------
def model_circle(r, nx, ny, nz=1):
	"""
	Create a centered circle (or sphere) having radius r.
	"""
	e = EMData()
	e.set_size(nx, ny, nz)
	e.process_inplace("testimage.circlesphere", {"radius":r, "fill":1})
	return e

def model_square(d, nx, ny, nz=1):
	"""
	Create a centered square (or cube) with edge length of d.
	"""
	e = EMData()
	e.set_size(nx, ny, nz)
	e.process_inplace("testimage.squarecube", {"edge_length":d, "fill":1})
	return e

def model_gauss(xsigma, nx, ny=1, nz=1, ysigma=None, zsigma=None, xcenter=None, ycenter=None, zcenter=None):
	"""
	Create a centered Gaussian image having standard deviation "sigma".
	"""
	e = EMData()
	e.set_size(nx, ny, nz)
	if( ysigma  == None ) : ysigma = xsigma
	if( zsigma  == None ) : zsigma = xsigma
	if( xcenter == None ) : xcenter = nx//2
	if( ycenter == None ) : ycenter = ny//2
	if( zcenter == None ) : zcenter = nz//2
	e.process_inplace("testimage.puregaussian", {"x_sigma":xsigma,"y_sigma":ysigma,"z_sigma":zsigma,"x_center":xcenter,"y_center":ycenter,"z_center":zcenter} )
	return e

def model_cylinder(radius, nx, ny, nz):
	"""
	 create a cylinder along z axis
	"""
	e = EMData()
	e.set_size(nx, ny, nz)
	e.process_inplace("testimage.cylinder", {"radius":radius})
	return  e

def model_gauss_noise(sigma, nx, ny=1, nz=1):
	"""
	Create an image of noise having standard deviation "sigma", 
	and average 0.
	"""
	e = EMData()
	e.set_size(nx, ny, nz)
	e.process_inplace("testimage.noise.gauss", {"sigma":sigma})
	return e

def model_blank(nx, ny=1, nz=1, bckg = 0.0):
	"""
	Create a blank image.
	"""
	e = EMData()
	e.set_size(nx, ny, nz)
	e.to_zero()
	if( bckg != 0.0):  e+=bckg
	return e

def set_seed(sde):
	from random import seed
	seed(int(sde))
	e = EMData()
	e.set_size(1,1,1)
	e.process_inplace("testimage.noise.gauss", {"sigma":1.0, "seed":int(sde)})

###----P-------
def parse_spider_fname(mystr, *fieldvals):
	"""
	Parse a Spider filename string and insert parameters.

	Example input: "foo{***}/img{****}.mrc"
	This string has two fields that should be replaced by integers,
	and the number of '*'s determines how "wide" that field should be.
	So, if the parameters to be inserted are 10 and 3, then the resulting
	filename should be "foo010/img0003.mrc".

	Note: If the image is a stack file, the last character in the string
	must be a '@' (except for possible extraneous whitespace, which is
	ignored).  This stack symbol will be stripped in the output filename.

	Example:

	   In [1]: mystr = "foo{***}/img{****}.mrc"
	   In [2]: parse_spider_fname(mystr, 10, 3)
	   Out[2]: 'foo010/img0003.mrc'

	@param mystr Spider filename string to be parsed
	@param fieldvals Integer values to be placed into the fields

	@return Parsed filename
	"""
	# helper functions and classes
	def rm_stack_char(mystr):
		"Helper function to remove a stack character if it exists"
		stackloc = mystr.find("@")
		if stackloc != -1: 
			# there's an '@' somewhere
			if len(mystr) - 1 == stackloc:
				# It's at the end of the string
				return mystr[:-1]
			else:
				# '@' not at the end, so it's an error
				raise ValueError, "Invalid format: misplaced '@'."
		else:
			# no '@' at all
			return mystr
	class Fieldloc:
		"Helper class to store description of a field"
		def __init__(self, begin, end):
			self.begin = begin
			self.end = end
		def count(self):
			"Size of the field (including braces)"
			return self.end - self.begin + 1
	def find_fields(mystr):
		"Helper function to identify and validate fields in a string"
		fields = []
		loc = 0
		while True:
			begin = mystr.find('{', loc)
			if begin == -1: break
			end = mystr.find('}', begin)
			field = Fieldloc(begin, end)
			# check validity
			asterisks = mystr[begin+1:end]
			if asterisks.strip("*") != "":
			    raise ValueError, "Malformed {*...*} field: %s" % \
				mystr[begin:end+1]
			fields.append(Fieldloc(begin, end))
			loc = end
		return fields
	# remove leading whitespace
	mystr.strip()
	# remove stack character (if it exists)
	mystr = rm_stack_char(mystr)
	# locate fields to replace
	fields = find_fields(mystr)
	if len(fields) != len(fieldvals):
		# wrong number of fields?
		raise ValueError, "Number of field values provided differs from" \
			"the number of {*...*} fields."
	newstrfrags = []
	loc = 0
	for i, field in enumerate(fields):
		# text before the field
		newstrfrags.append(mystr[loc:field.begin])
		# replace the field with the field value
		fieldsize = field.count() - 2
		fielddesc = "%0" + str(fieldsize) + "d"
		newstrfrags.append(fielddesc % fieldvals[i])
		loc = field.end + 1
	newstrfrags.append(mystr[loc:])
	return "".join(newstrfrags)

def peak_search(e, npeak = 1, invert = 1, print_screen = 0):
	peaks    = e.peak_search(npeak, invert)
	ndim     = peaks[0]
	nlist    = int((len(peaks)-1)/((ndim+1)*2))
	if(nlist > 0):
		outpeaks = []
		if(print_screen):
			if  ndim == 1 : 
				print		      '%10s%10s%10s%10s%10s'%("Index  "," Peak_value","X   ",		     "Peak/P_max", "X-NX/2")
				print_list_format(peaks[1:], 4)
			elif ndim == 2 :
				print	      '%10s%10s%10s%10s%10s%10s%10s'%("Index  ", "Peak_value","X   ","Y   ",	     "Peak/P_max", "X-NX/2", "Y-NY/2")
				print_list_format(peaks[1:], 6)
			elif ndim == 3 : 
				print '%10s%10s%10s%10s%10s%10s%10s%10s%10s'%("Index  ", "Peak_value","X   ","Y   ","Z   ", "Peak/P_max", "X-NX/2", "Y-NY/2", "Z-NZ/2")
				print_list_format(peaks[1:], 8)
			else:	ERROR("Image dimension extracted in peak_search is wrong", "Util.peak_search", 1)
		for i in xrange(nlist):
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
def print_row(input, ix=0, iz=0):
	"""Print the data in slice iz, row ix of an image to standard out.

	Usage: print_row(image, ix, iz)
	   or
	       print_row("path/to/image", ix, iz)
	"""
	image=get_image(input)
	nx = image.get_xsize()
	ny = image.get_ysize()
	nz = image.get_zsize()
	print "(z = %d slice, x = %d row)" % (iz, ix)
	line = []
	for iy in xrange(ny):
		line.append("%12.5g  " % (image.get_value_at(ix,iy,iz)))
		if ((iy + 1) % 5 == 0): line.append("\n   ")
	line.append("\n")
	print "".join(line)

def print_col(input, iy=0, iz=0):
	"""Print the data in slice iz, column iy of an image to standard out.

	   Usage: print_col(image, iy, iz)
	      or
		  print_col("path/to/image", iy, iz)
	"""
	image=get_image(input)
	nx = image.get_xsize()
	ny = image.get_ysize()
	nz = image.get_zsize()
	print "(z = %d slice, y = %d col)" % (iz, iy)
	line = []
	for ix in xrange(nx):
		line.append("%12.5g  " % (image.get_value_at(ix,iy,iz)))
		if ((ix + 1) % 5 == 0): line.append("\n   ")
	line.append("\n")
	print "".join(line)

def print_slice(input, iz=0):
	"""Print the data in slice iz of an image to standard out.

	Usage: print_image(image, int)
	   or
	       print_image("path/to/image", int)
	"""
	image=get_image(input)
	nx = image.get_xsize()
	ny = image.get_ysize()
	nz = image.get_zsize()
	print "(z = %d slice)" % (iz)
	line = []
	for iy in xrange(ny):
		line.append("Row ")
		line.append("%4i " % iy)
		for ix in xrange(nx):
			line.append("%12.5g  " % (image.get_value_at(ix,iy,iz)))
			if ((ix + 1) % 5 == 0): 
				line.append("\n   ")
				line.append("      ")
	    	line.append("\n")
	    	if(nx%5 != 0): line.append("\n")
	print "".join(line)

def print_image(input):
	"""Print the data in an image to standard out.

	Usage: print_image(image)
	   or
	       print_image("path/to/image")
	"""
	image=get_image(input)
	nz = image.get_zsize()
	for iz in xrange(nz): print_slice(input, iz)


def print_image_col(input, ix=0, iz=0):
	"""Print the data in slice iz, row ix of an image to standard out.

	Usage: print_image_col(image, ix, iz)
	   or
	       print_image_col("path/to/image", ix, iz)
	"""
	image=get_image(input)
	nx = image.get_xsize()
	ny = image.get_ysize()
	nz = image.get_zsize()
	print "(z = %d slice, x = %d row)" % (iz, ix)
	line = []
	for iy in xrange(ny):
		line.append("%12.5g  " % (image.get_value_at(ix,iy,iz)))
		if ((iy + 1) % 5 == 0): line.append("\n   ")
	line.append("\n")
	print "".join(line)

def print_image_row(input, iy=0, iz=0):
	"""Print the data in slice iz, column iy of an image to standard out.

	   Usage: print_image_row(image, iy, iz)
	      or
		  print_image_row("path/to/image", iy, iz)
	"""
	image=get_image(input)
	nx = image.get_xsize()
	ny = image.get_ysize()
	nz = image.get_zsize()
	print "(z = %d slice, y = %d col)" % (iz, iy)
	line = []
	for ix in xrange(nx):
		line.append("%12.5g  " % (image.get_value_at(ix,iy,iz)))
		if ((ix + 1) % 5 == 0): line.append("\n   ")
	line.append("\n")
	print "".join(line)

def print_image_slice(input, iz=0):
	"""Print the data in slice iz of an image to standard out in a format that agrees with v2

	Usage: print_image_slice(image, int)
	   or
	       print_image_slice("path/to/image", int)
	"""
	image=get_image(input)
	nx = image.get_xsize()
	ny = image.get_ysize()
	nz = image.get_zsize()
	print "(z = %d slice)" % (iz)
	line = []
	for iy in xrange(ny-1,-1,-1):
		line.append("Row ")
		line.append("%4i " % iy)
		for ix in xrange(nx):
			line.append("%12.5g  " % (image.get_value_at(ix,iy,iz)))
			if ((ix + 1) % 5 == 0): 
				line.append("\n   ")
				line.append("      ")
	    	line.append("\n")
	    	if(nx%5 != 0): line.append("\n")
	print "".join(line)
	
def print_image_slice_3d(input, num=0,direction="z"):
	"""Print the data in slice iz of an image to standard out in a format that agrees with v2

	Usage: print_image_slice(image, int)
	   or
	       print_image_slice("path/to/image", int)
	"""
	#print "print slice at 3 directions"
	image=get_image(input)
	nx = image.get_xsize()
	ny = image.get_ysize()
	nz = image.get_zsize()
	if(direction=="x"):
		#print "xxxxx"
		ix=num
		print "(x = %d slice)" % (ix)
		line = []
		for iz in xrange(nz-1,-1,-1):
			line.append("Z ")
			line.append("%4i " % iz)
			for iy in xrange(ny):
				line.append("%12.5g  " % (image.get_value_at(ix,iy,iz)))
				if ((iy + 1) % 5 == 0): 
					line.append("\n   ")
					line.append("      ")
	    		line.append("\n")
	    		if(ny%5 != 0): line.append("\n")
		print "".join(line)
	elif(direction=="y"):
		#print "yyy"
		iy=num
		print "(y = %d slice)" % (iy)
		line = []
		for iz in xrange(nz-1,-1,-1):
			line.append("Z ")
			line.append("%4i " % iz)
			for ix in xrange(nx):
				line.append("%12.5g  " % (image.get_value_at(ix,iy,iz)))
				if ((ix + 1) % 5 == 0): 
					line.append("\n   ")
					line.append("      ")
	    		line.append("\n")
	    		if(nx%5 != 0): line.append("\n")
		print "".join(line)
	else:
		#print "zzzz"
		iz=num
		print "(z = %d slice)" % (iz)
		line = []
		for iy in xrange(ny-1,-1,-1):
			line.append("Row ")
			line.append("%4i " % iy)
			for ix in xrange(nx):
				line.append("%12.5g  " % (image.get_value_at(ix,iy,iz)))
				if ((ix + 1) % 5 == 0): 
					line.append("\n   ")
					line.append("      ")
	    		line.append("\n")
	    		if(nx%5 != 0): line.append("\n")
		print "".join(line)
		

def print_list_format(m, narray = 0):
	from string 	import split
	from math 	import sqrt
	import string
	import types
	"""
		Print formated elements in a list to screen
		The screen output is in the form of narray*int(len(m)/narray)
		Or when narray is zero, int(sqrt(len(m)))*int(sqrt(len(m)))
	"""
	flist = []
	for i in xrange(len(m)):
		if   type(m[i])  is types.FloatType: flist.append('%10.3g'%(m[i]))
		elif type(m[i])  is types.IntType :  flist.append(  '%10d'%(m[i]))
		else				   : flist.append(  '%10s'%(m[i]))
	if(narray > len(m)):
		narray = 0
		ERROR("improper input narray number, use default value", "print_list_foramt",0)
	if(narray == 0 ):
		num = int(sqrt(len(m)))
		if( len(m) % num != 0): lnum = int(len(m)/num) + 1
		else: 			lnum = int(len(m)/num)
	else:
		num = narray
		if( len(m) % num == 0): lnum = int(len(m)/num)
		else: 			lnum = int(len(m)/num) + 1
	ncount = -1
	plist  = []
	for i in xrange(lnum):
		qlist = ""
		for j in xrange(num):
			ncount += 1
			if ncount <= len(m) - 1: qlist=qlist+flist[ncount]
			else:			 break
		plist.append(qlist)
	for i in xrange(lnum):
		print '%6d '%(i+1),plist[i]

def pad(image_to_be_padded, new_nx, new_ny = 1,	new_nz = 1, background = "average", off_center_nx = 0, off_center_ny = 0, off_center_nz = 0):
	import types
	if type(background) != types.StringType: background = str(background)
	if   background == "average"       :     image_padded = Util.pad(image_to_be_padded, new_nx, new_ny, new_nz, off_center_nx, off_center_ny, off_center_nz, "average")
	elif background == "circumference" :     image_padded = Util.pad(image_to_be_padded, new_nx, new_ny, new_nz, off_center_nx, off_center_ny, off_center_nz, "circumference")
	else:                                    image_padded = Util.pad(image_to_be_padded, new_nx, new_ny, new_nz, off_center_nx, off_center_ny, off_center_nz,  background  )
	return  image_padded

def read_spider_doc(fnam):
	from string import atof, atoi
	"""
		spider doc file format:
		key nrec record ...
		5   2    12     ...(when key <=99999)
		6   2    12     ...(when key >99999)
	"""
	inf = file(fnam, "r")
	comment_line=inf.readline() # assume there is only one comment line
	docf_in = []
	data    = []
	line    = inf.readline()
	while len(line) > 0:
		line_data=[]
		if(line[5:6] == " "): ibeg = 6
		else:  ibeg = 7
		for irec in xrange(atoi(line[ibeg:ibeg+2])):
		 	start= ibeg+2+irec*12
		 	end  = ibeg+2+(irec+1)*12
		 	line_data.append(atof(line[start:end]))
		data.append(line_data)
		line = inf.readline()
	return data
		
def read_text_row(fnam, format="", skip=";"):
	"""
	 	Read a column-listed txt file.
		INPUT: filename: name of the Doc file
	 	OUTPUT:
	    	nc : number of entries in each lines (number of columns)
	    	len(data)/nc : number of lines (rows)
	    	data: List of numbers from the doc file
 	"""
	from string import split

	inf  = file(fnam, "r")
	strg = inf.readline()
	x    = []
	data = []
	while (len(strg) > 0):
		com_line = False
		for j in xrange(len(strg)):
			if(strg[j] == skip):	com_line = True
		if com_line == False:
			word=split(strg)
			if format == "s" :
				key = int(word[1])
				if key != len(word) - 2:
					del word
					word = []
					word.append(strg[0 : 5])
					word.append(strg[6 : 7])
					for k in xrange(key):
						k_start = 7       + k*13
						k_stop  = k_start + 13
						word.append(strg[k_start : k_stop])				
			line=[]
			for i in xrange(len(word)):
				line.append(float(word[i]))
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
	import types
	outf = open(file_name, "w")
	if (type(data[0]) == types.ListType):
		# It is a list of lists
		for i in xrange(len(data)):
			for j in xrange(len(data[i])):
				outf.write("  %12.5g"%data[i][j])
			outf.write("\n")
	else:
		# Single list
		for j in xrange(len(data)):
			outf.write("  %12.5g"%data[j])
		outf.write("  \n")
	outf.close()


def read_text_file(file_name, ncol = 0):
	"""
		Read data from text file, if ncol = -1, read all columns
		if ncol >= 0, just read the (ncol+1)-th column.
	"""
	
	from string import split
	inf = file(file_name, "r")
	line = inf.readline()
	data = []
	while len(line) > 0:
		if ncol == -1:
			vdata = split(line)
			if data == []:
				for i in xrange(len(vdata)):
					data.append([float(vdata[i])])
			else:
				for i in xrange(len(vdata)):
					data[i].append(float(vdata[i]))			
		else:
			vdata = float(split(line)[ncol])
			data.append(vdata)
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
	import types
	outf = open(file_name, "w")
	if (type(data[0]) == types.ListType):
		# It is a list of lists
		for i in xrange(len(data[0])):
			for j in xrange(len(data)):
				outf.write("  %12.5g"%data[j][i])
			outf.write("\n")
	else:
		# Single list
		for j in xrange(len(data)):
			outf.write("  %12.5g\n"%data[j])
	outf.close()

def reconstitute_mask(image_mask_applied_file, new_mask_file, save_file_on_disk = True, saved_file_name = "image_in_reconstituted_mask.hdf"):
	import types
	"""
		Substitute masked area value with image average 
	"""
	if type(image_mask_applied_file) == types.StringType:
		nima = EMUtil.get_image_count(image_mask_applied_file)
		if (nima > 1):
			image_mask_applied = []
		 	for ima in xrange(nima):
				e = EMData()
				e.read_image(image_mask_applied_file, ima)
				image_mask_applied.append(e)
		else:
			image_mask_applied = get_im(image_mask_applied_file)
	elif  type(image_mask_applied_file) == types.ListType:   
		nima =  len( image_mask_applied )
		image_mask_applied = image_mask_applied_file
	if type(new_mask_file) == types.StringType:     
		new_mask = get_im( new_mask_file )
        elif type(new_mask_file) == types.IntType or type( new_mask_file ) == types.floatType:
		if nima > 1:
			e = image_mask_applied[0]
			nx = e.get_xsize()
			ny = e.get_ysize()
			nz = e.get_zsize()
		else :
			nx = image_mask_applied.get_xsize()
			ny = image_mask_applied.get_ysize()
			nz = image_mask_applied.get_zsize()
		new_mask = model_circle(new_mask_file, nx, ny, nz)
	if  nima > 1:
		image_in_reconstituted_mask = []
		for i in xrange(nima):
			tmp_image = Util.reconstitute_image_mask(image_mask_applied[i], new_mask)
			image_in_reconstituted_mask.append (tmp_image)
			if (save_file_on_disk ):  image_in_reconstituted_mask[i].write_image(saved_file_name, i)
		if(not save_file_on_disk):  return  image_in_reconstituted_mask
	else :  
		if(save_file_on_disk ):
			image_in_reconstituted_mask = Util.reconstitute_image_mask(image_mask_applied, new_mask)
			image_in_reconstituted_mask.write_image(saved_file_name)
		else:	return  Util.reconstitute_image_mask(image_mask_applied, new_mask)

def rotate_about_center(alpha, cx, cy):
	"""Rotate about a different center

	    Usage: rotate_about_center(alpha,cx,cy):
	       angles in degrees
	"""

	cmp1 = compose_transform2(0, -cx, -cy, 1, alpha, 0, 0, 1)
	cmp2 = compose_transform2(cmp1[0], cmp1[1], cmp1[2], cmp1[3], 0, cx, cy, 1)

	#   return compalpha, comptrans.at(0),comptrans.at(1), compscale
	return cmp2[0], cmp2[1], cmp2[2], cmp2[3]
	
def reshape_1dpw2(input_object, rimage_size_to_be, rimage_size_interpolated, Pixel_size_to_be = 0, Pixel_size_interpolated = 0):
	"""
		linearly interpolate a 1D power spectrum to required length with required Pixel size
	""" 
	import types
	from utilities import read_spider_doc
	from fundamentals import rops_table
	if  type(input_object) == types.ClassType:
		to_be_interpolated = []
		if input_object.get_ndim() > 1: 
			to_be_interpolated = rops_table(input_object)
			rimage_size_to_be  = input_object.get_xsize()
		else       : 
			for ir in xrange(input_object.get_xsize()): to_be_interpolated.append(input_object.get_value_at(ir))
	elif type(input_object) == types.StringType:
		to_be_interpolated = [] 
		tmp_table = read_spider_doc(input_object) # assume 1dpw2 is saved in SPIDER format
		for i in xrange(len(tmp_table)):
			to_be_interpolated.append(tmp_table[i][0])
	elif type(input_object) == types.ListType or type(input_object) == types.TupleType: to_be_interpolated = input_object
	else:   ERROR("Input type is wrong ","linear_interpolation_1dpw2",1)
	interpolated = []
	if  Pixel_size_to_be == 0 :
		Pixel_size_to_be = 1
		Pixel_size_interpolated = Pixel_size_to_be 
	if  Pixel_size_interpolated == 0: Pixel_size_interpolated = Pixel_size_to_be				
	for i in xrange(int(rimage_size_interpolated/2)):
		xi =  float(i)*(rimage_size_to_be*Pixel_size_to_be)/(Pixel_size_interpolated*rimage_size_interpolated)
		ix = int(xi)
		df = xi -ix
		xval = (1.0-df)*to_be_interpolated[ix] + df*to_be_interpolated[ix+1]
		interpolated.append(xval)
	return interpolated
	
def rops_dir(indir, output_dir = "1dpw2_dir"):
	"""
		Calculate 1D rotationally averaged power spectra from 
		image stack listed in a directory
	"""
	import os
	flist = os.listdir(indir)
	print flist
	if os.path.exists(output_dir) is False: os.mkdir(output_dir)
	for i, v in enumerate(flist):
		(filename, filextension) = os.path.splitext(v)
		nima = EMUtil.get_image_count(os.path.join(indir,v))
		print nima
		for im in xrange(nima):
			e = EMData()
			file_name = os.path.join(indir,v)
			e.read_image(file_name, im)
			tmp1 = periodogram(e)
			tmp  = tmp1.rotavg()
			if im == 0:
				sum_ima  = model_blank(tmp.get_xsize())
				sum_ima += tmp
			else :  sum_ima += tmp
		table = []
		nr = sum_ima.get_xsize()
		for ir in xrange(nr):  table.append([sum_ima.get_value_at(ir)])
		drop_spider_doc(os.path.join(output_dir, "1dpw2_"+filename+".txt"), table)


def estimate_3D_center(data):
	from math import cos, sin, pi
	from numpy import matrix
	from numpy import linalg
	
	ali_params = []
	for im in data:
		phi, theta, psi, s2x, s2y = get_params_proj(im)
		ali_params.append([phi, theta, psi, s2x, s2y])

	N = len(ali_params)
	A = []
	b = []
	
	for i in xrange(N):
		phi_rad = ali_params[i][0]/180*pi
		theta_rad = ali_params[i][1]/180*pi
		psi_rad = ali_params[i][2]/180*pi
		A.append([cos(psi_rad)*cos(theta_rad)*cos(phi_rad)-sin(psi_rad)*sin(phi_rad), 
			cos(psi_rad)*cos(theta_rad)*sin(phi_rad)+sin(psi_rad)*cos(phi_rad), -cos(psi_rad)*sin(theta_rad), 1, 0])
		A.append([-sin(psi_rad)*cos(theta_rad)*cos(phi_rad)-cos(psi_rad)*sin(phi_rad), 
			-sin(psi_rad)*cos(theta_rad)*sin(phi_rad)+cos(psi_rad)*cos(phi_rad), sin(psi_rad)*sin(theta_rad), 0, 1])	
		b.append([ali_params[i][3]])
		b.append([ali_params[i][4]])
	
	A_matrix = matrix(A)
	b_matrix = matrix(b)

	K = linalg.solve(A_matrix.T*A_matrix, A_matrix.T*b_matrix)

	return float(K[0][0]), float(K[1][0]), float(K[2][0]), float(K[3][0]), float(K[4][0])


def estimate_3D_center_MPI(data, nima, myid, number_of_proc, main_node):
	from math import cos, sin, pi
	from numpy import matrix
	from numpy import linalg
	from mpi import MPI_COMM_WORLD
	from mpi import mpi_recv, mpi_send, MPI_FLOAT
	from applications import MPI_start_end
	
	ali_params_series = []
	for im in data:
		phi, theta, psi, s2x, s2y = get_params_proj(im)
		ali_params_series.append(phi)
		ali_params_series.append(theta)
		ali_params_series.append(psi)
		ali_params_series.append(s2x)
		ali_params_series.append(s2y)

	if myid == main_node:
		for proc in xrange(number_of_proc):
			if proc != main_node:
				image_start_proc, image_end_proc = MPI_start_end(nima, number_of_proc, proc)
				n_params = (image_end_proc - image_start_proc)*5
				temp = mpi_recv(n_params, MPI_FLOAT, proc, proc, MPI_COMM_WORLD)
				for nn in xrange(n_params): 	ali_params_series.append(float(temp[nn]))
					
		ali_params = []
		N = len(ali_params_series)/5
		for im in xrange(N):
			ali_params.append([ali_params_series[im*5], ali_params_series[im*5+1], ali_params_series[im*5+2], ali_params_series[im*5+3], ali_params_series[im*5+4]])

		A = []
		b = []
		DEG_to_RAD = pi/180.0
	
		for i in xrange(N):
			phi_rad = ali_params[i][0]*DEG_to_RAD
			theta_rad = ali_params[i][1]*DEG_to_RAD
			psi_rad = ali_params[i][2]*DEG_to_RAD
			A.append([cos(psi_rad)*cos(theta_rad)*cos(phi_rad)-sin(psi_rad)*sin(phi_rad), 
				cos(psi_rad)*cos(theta_rad)*sin(phi_rad)+sin(psi_rad)*cos(phi_rad), -cos(psi_rad)*sin(theta_rad), 1, 0])
			A.append([-sin(psi_rad)*cos(theta_rad)*cos(phi_rad)-cos(psi_rad)*sin(phi_rad), 
				-sin(psi_rad)*cos(theta_rad)*sin(phi_rad)+cos(psi_rad)*cos(phi_rad), sin(psi_rad)*sin(theta_rad), 0, 1])	
			b.append([ali_params[i][3]])
			b.append([ali_params[i][4]])
		
		A_matrix = matrix(A)
		b_matrix = matrix(b)

		K = linalg.solve(A_matrix.T*A_matrix, A_matrix.T*b_matrix)
		return float(K[0][0]), float(K[1][0]), float(K[2][0]), float(K[3][0]), float(K[4][0])

	else:
		image_start_proc, image_end_proc = MPI_start_end(nima, number_of_proc, myid)
		n_params = (image_end_proc - image_start_proc)*5
		mpi_send(ali_params_series, n_params, MPI_FLOAT, main_node, myid, MPI_COMM_WORLD)
		
		return 0.0, 0.0, 0.0, 0.0, 0.0	

def rotate_3D_shift(data, shift3d):

	t = Transform({"type":"spider","phi":0.0,"theta":0.0,"psi":0.0,"tx":-shift3d[0],"ty":-shift3d[1],"tz":-shift3d[2],"mirror":0,"scale":1.0})

	for i in xrange(len(data)):
		d = data[i].get_attr('xform.projection')
		c = d*t
		data[i].set_attr('xform.projection', c)


def sym_vol(image, symmetry="c1"):
	" Symmetrize a volume"
	if(symmetry == "c1"):  return  image.copy()
	else:                  return  image.symvol(symmetry)

##----------------------------------HDF headers related code --------------------------
def set_arb_params(img, params, par_str):

	"""
		filling arbitary headers 
	"""
	for i in xrange(len(par_str)): img.set_attr_dict({par_str[i]:params[i]})

def get_arb_params(img, par_str):

	"""
		reading arbitary headers 
	"""
	params=[]
	for i in xrange(len(par_str)): params.append(img.get_attr(par_str[i]))
	return params

###------------------------------------------------------------------------------------------	

def start_time():
	import time
	start_time = time.time()
	return start_time

def finish_time(start_time):
	import time
	finish_time = time.time()
	print ("Running time is"), finish_time-start_time
	return finish_time

def ttime():
	import time
	now = time.localtime(time.time())
	return time.asctime(now)

def running_time(start_time):
	from utilities import print_msg
	from time import time
	time_run = int(time() - start_time)
	time_h   = time_run / 3600
	time_m   = (time_run % 3600) / 60
	time_s   = (time_run % 3600) % 60
	print_msg('\nRunning time is: %s h %s min %s s\n\n' % (str(time_h).rjust(2, '0'), str(time_m).rjust(2, '0'), str(time_s).rjust(2, '0')))

def running_time_txt(start_time):
	from time import time
	time_run = int(time() - start_time)
	time_h   = time_run / 3600
	time_m   = (time_run % 3600) / 60
	time_s   = (time_run % 3600) % 60
	return 'Running time is: %s h %s min %s s' % (str(time_h).rjust(2, '0'), str(time_m).rjust(2, '0'), str(time_s).rjust(2, '0'))

'''
def reduce_array_to_root(data, myid, main_node = 0, comm = -1):
	from numpy import array, shape, reshape
	from mpi import MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD, mpi_reduce, mpi_barrier
	
	if comm == -1:  comm = MPI_COMM_WORLD
	n = shape(data)
	ntot = 1
	for i in xrange(len(n)):  ntot *= n[i]
        count = 500000
        array1d = reshape(data, (ntot,))
        ntime = (ntot-1) /count + 1
        for i in xrange(ntime):
        	block_begin = i*count
        	block_end   = i*count + count
        	if block_end > ntot: block_end = ntot
        	block_size  = block_end - block_begin
        	tmpsum = mpi_reduce(array1d[block_begin], block_size, MPI_FLOAT, MPI_SUM, main_node, comm)
		mpi_barrier(comm)
        	if myid == main_node:
        		array1d[block_begin:block_end] = tmpsum[0:block_size]
'''

def reduce_EMData_to_root(data, myid, main_node = 0, comm = -1):
	from numpy import array, shape, reshape
	from mpi   import mpi_reduce, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD, mpi_barrier
	from utilities import get_image_data
	
	if comm == -1:  comm = MPI_COMM_WORLD

	array = get_image_data(data)
	n     = shape(array)
	ntot  = 1
	for i in n: ntot *= i
	count = (75*4+2)*(75*4)**2
	array1d = reshape( array, (ntot,))
	ntime = (ntot-1) /count + 1
	for i in xrange(ntime):
		block_begin = i*count
		block_end   = min(block_begin + count, ntot)
		block_size  = block_end - block_begin
		tmpsum = mpi_reduce(array1d[block_begin:block_begin+block_size], block_size, MPI_FLOAT, MPI_SUM, main_node, comm)
		mpi_barrier(comm)
		if myid == main_node:
			array1d[block_begin:block_end] = tmpsum[0:block_size]

def bcast_EMData_to_all(tavg, myid, source_node = 0, comm = -1):
	from numpy import array, shape, reshape
	from mpi   import mpi_bcast, MPI_FLOAT, MPI_COMM_WORLD

	if comm == -1: comm = MPI_COMM_WORLD
	tavg_data = EMNumPy.em2numpy(tavg)
	n = shape(tavg_data)
	ntot = 1
	for i in n: ntot *= i
	tavg_tmp = mpi_bcast(tavg_data, ntot, MPI_FLOAT, source_node, comm)
	if(myid != source_node):
		tavg_data1d = reshape(tavg_data,(ntot,))
		tavg_data1d[0:ntot] = tavg_tmp[0:ntot]
	
'''
def bcast_EMData_to_all(img, myid, main_node = 0, comm = -1):

	# Comment by Zhengfan Yang on 01/05/10
	#
	# Notice: 
	# (1) one should use this new version of broadcasting EMData in the following way:
	# 	img = bcast_EMData_to_all(img, myid, main_node, comm)
	# instead of
	# 	bcast_EMData_to_all(img, myid, main_node, comm)
	# The latter is inconsistent with mpi_bcast() and difficult to implement efficiently
	#
	# (2) To be consistent with send_EMData() and recv_EMData(), we assume that the node 
	# other than the broadcasting node know nothing about the EMData(). Therefore, one 
	# need to broadcast the size of EMData() and two attributes: is_complex and is_ri. 
	# For all other attributes, you are on your own.

	from numpy import reshape
	from mpi import mpi_bcast, MPI_INT, MPI_FLOAT, MPI_COMM_WORLD

	if comm == -1: comm = MPI_COMM_WORLD
	
	img_head = []
	if myid == main_node:
		img_head.append(img.get_xsize())
		img_head.append(img.get_ysize())
		img_head.append(img.get_zsize())
		img_head.append(img.is_complex())
		img_head.append(img.is_ri())
	img_head = mpi_bcast(img_head, 5, MPI_INT, main_node, comm)
	nx = int(img_head[0])
	ny = int(img_head[1])
	nz = int(img_head[2])
	is_complex = int(img_head[3])
	is_ri = int(img_head[4])

	ntot = nx*ny*nz
	img_data = EMNumPy.em2numpy(img)
	img_data = mpi_bcast(img_data, ntot, MPI_FLOAT, main_node, comm)
	if nz != 1:
		img_data = reshape(img_data, (nz, ny, nx))     # For some reason, the order should be like this -- Zhengfan Yang
	elif ny != 1:
		img_data = reshape(img_data, (ny, nx))
	else:
		pass
	img = EMNumPy.numpy2em(img_data)
	img.set_complex(is_complex)
	img.set_ri(is_ri)

	return img


def reduce_EMData_to_root(img, myid, main_node = 0, comm = -1):

	# Comment by Zhengfan Yang on 01/05/10
	#
	# Notice: 
	# (1) one should use this new version of reducing EMData in the following way:
	# 	img = reduce_EMData_to_root(img, myid, main_node, comm)
	# instead of
	# 	reduce_EMData_to_root(img, myid, main_node, comm)
	# The latter is inconsistent with mpi_bcast() and difficult to implement efficiently

	from numpy import reshape
	from mpi   import mpi_reduce, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD
	
	if comm == -1: comm = MPI_COMM_WORLD
	
	nx = img.get_xsize()
	ny = img.get_ysize()
	nz = img.get_zsize()
	is_complex = img.is_complex()
	is_ri = img.is_ri()
	ntot = nx*ny*nz
	
	img_data = EMNumPy.em2numpy(img)
	img_data = mpi_reduce(img_data, ntot, MPI_FLOAT, MPI_SUM, main_node, comm)

	if myid == main_node:
		if nz!=1:
			img_data = reshape(img_data, (nz, ny, nx))
		elif ny!=1:
			img_data = reshape(img_data, (ny, nx))
		else:
			pass

		img = EMNumPy.numpy2em(img_data)
		img.set_complex(is_complex)
		img.set_ri(is_ri)
		return img
	else:
		return img
'''

def send_EMData(img, dst, tag, comm=-1):
	from mpi import mpi_send, MPI_INT, MPI_FLOAT, MPI_COMM_WORLD
	
	if comm == -1: comm = MPI_COMM_WORLD
	img_head = []
	img_head.append(img.get_xsize())
	img_head.append(img.get_ysize())
	img_head.append(img.get_zsize())
	img_head.append(img.is_complex())
	img_head.append(img.is_ri())

	head_tag = 2*tag
	mpi_send(img_head, 5, MPI_INT, dst, head_tag, comm)

	img_data = get_image_data(img)
	data_tag = 2*tag+1
	ntot = img_head[0]*img_head[1]*img_head[2]
	mpi_send(img_data, ntot, MPI_FLOAT, dst, data_tag, comm)

        '''
	count = 100000
	data1d = reshape(img_data, (ntot,))
	ntime = (ntot-1) /count + 1

	for i in xrange(ntime):
		block_begin = i*count
		block_end   = i*count + count
		if block_end > ntot:
		    block_end = ntot
		block_size  = block_end - block_begin
		mpi_send(data1d[block_begin], block_size, MPI_FLOAT, dst, data_tag*ntime+i, comm)
        '''

def recv_EMData(src, tag, comm=-1):
	from mpi import mpi_recv, MPI_INT, MPI_FLOAT, MPI_COMM_WORLD
	from numpy import reshape
	
	if comm==-1: comm = MPI_COMM_WORLD
	head_tag = 2*tag
	img_head = mpi_recv(5, MPI_INT, src, head_tag, comm)

	nx = int(img_head[0])
	ny = int(img_head[1])
	nz = int(img_head[2])
	is_complex = int(img_head[3])
	is_ri = int(img_head[4])

	data_tag = 2*tag+1
	ntot = nx*ny*nz
	
	img_data = mpi_recv(ntot, MPI_FLOAT, src, data_tag, comm)
	if nz != 1:
		img_data = reshape(img_data, (nz, ny, nx))
	elif ny != 1:
		img_data = reshape(img_data, (ny, nx))
	else:
		pass

	img = EMNumPy.numpy2em(img_data)
	img.set_complex(is_complex)
	img.set_ri(is_ri)
	return img
		
	'''
	#construct a EMData by taking the ownership of numpy array, no memory copying  --Grant Tang
	#recv_data_numeric = mpi_recv(ntot, MPI_FLOAT, src, data_tag, comm)
	#recv_data_numpy = numpy.array(recv_data_numeric)
	#numpy_data = recv_data.reshape(recv_data, (nz,ny,nx))
	#img = EMNumPy.numpy2em(numpy_data)

	#comment out Wei's original code, which makes memory copy to construct EMData from numpy array  --Grant Tang
	img = EMData()
	img.set_size(nx, ny, nz)
	if( complex > 0 ):
		img.set_complex(True)
	else:
		img.set_complex(False)

	data1d = reshape( get_image_data(img), (ntot,) )
	tmp_data = mpi_recv(ntot, MPI_FLOAT, src, data_tag, comm)
        data1d[0:ntot] = tmp_data[0:ntot] 
        
        
	count = 100000
	ntime = (ntot-1)/count + 1

	for i in xrange(ntime):
		block_begin = i*count
		block_end   = i*count + count
		if block_end > ntot:
		    block_end = ntot
		block_size  = block_end - block_begin
		tmp_data = mpi_recv(block_size, MPI_FLOAT, src, data_tag*ntime+i, comm)
		data1d[block_begin:block_end] = tmp_data[0:block_size]
        
	return img
	'''

def gather_EMData(data, number_of_proc, myid, main_node):
	"""
	Gather the a list of EMData on all nodes to the main node, we assume the list has the same length on each node.
	"""
	
	l = len(data)
	gathered_data = []
	if myid == main_node:
		for i in xrange(number_of_proc):
			if i == main_node:
				for k in xrange(l):
					gathered_data.append(data[k])
			else:
				for k in xrange(l):
					gathered_data.append(recv_EMData(i, i*l+k))
	else:
		for k in xrange(l):
			send_EMData(data[k], main_node, myid*l+k)
	return gathered_data
	

def bcast_string_to_all(str_to_send, source_node = 0):
	from mpi import mpi_bcast, MPI_INT, MPI_COMM_WORLD
	str_tmp = ""
	str_TMP = mpi_bcast(str_to_send, len(str_to_send), MPI_INT, source_node, MPI_COMM_WORLD)
	for i in xrange(len(str_to_send)):  str_tmp += chr(str_TMP[i])
	return str_tmp
	
def bcast_number_to_all(number_to_send, source_node = 0):
	"""
		number_to_send has to be pre-defined in each node
	"""
	from mpi import mpi_bcast, MPI_INT, MPI_COMM_WORLD, MPI_FLOAT
	import types
	if    type(number_to_send) is types.IntType: 
		TMP = mpi_bcast(number_to_send, 1, MPI_INT,   source_node, MPI_COMM_WORLD)
		return int(TMP[0])
	elif  type(number_to_send) is types.FloatType:
		TMP = mpi_bcast(number_to_send, 1, MPI_FLOAT, source_node, MPI_COMM_WORLD)
		return float(TMP[0])
	else:
		print  " ERROR "
	
def bcast_list_to_all(list_to_send, source_node = 0):
	from mpi import mpi_bcast, MPI_COMM_WORLD, MPI_FLOAT
	import   types
	list_tmp = mpi_bcast(list_to_send, len(list_to_send), MPI_FLOAT, source_node, MPI_COMM_WORLD)
	list_to_bcast = []
	for i in xrange(len(list_to_send)):
		if (type(list_to_send[i]) == types.IntType ): list_to_bcast.append( int(  list_tmp[i] ) )
		else:                                        list_to_bcast.append( float( list_tmp[i] ) )
	return list_to_bcast

def recv_attr_dict(main_node, stack, data, list_params, image_start, image_end, number_of_proc, comm = -1):
	import types
	from  utilities import  get_arb_params, set_arb_params
	from  mpi 	import mpi_recv
	from  mpi 	import MPI_FLOAT, MPI_INT, MPI_TAG_UB, MPI_COMM_WORLD
	#   hdf version!
	# This is done on the main node, so for images from the main node, simply write headers
	
	if comm == -1:  comm = MPI_COMM_WORLD
	
	TransType = type(Transform())
	# prepare keys for float/int
	value = get_arb_params(data[0], list_params)
	ink = []
	len_list = 0
	for il in xrange(len(list_params)):
		if type(value[il]) is types.IntType:     
			ink.append(1)
			len_list += 1
		elif type(value[il]) is types.FloatType:  
			ink.append(0)
			len_list += 1
		elif type(value[il]) is TransType:
			ink.append(2)
			len_list += 12
	ldis = []
	headers = []
	for n in xrange(number_of_proc):
		if n != main_node:
			dis = mpi_recv(2, MPI_INT, n, MPI_TAG_UB, comm)
			value = mpi_recv(len_list*(dis[1]-dis[0]), MPI_FLOAT, n, MPI_TAG_UB, comm)
			ldis.append([dis[0], dis[1]])
			headers.append(value)
			del  dis
	del  value
	for im in xrange(image_start, image_end):
		data[im-image_start].write_image(stack, data[im-image_start].get_attr_default('ID', im), EMUtil.ImageType.IMAGE_HDF, True)

	for n in xrange(len(ldis)):
		img_begin = ldis[n][0]
		img_end = ldis[n][1]
		for im in xrange(img_begin, img_end):
			par_begin = (im-img_begin)*len_list
			nvalue = []
			header = headers[n]
			ilis = 0
			for il in xrange(len(list_params)):
				if(ink[il] == 1):
					nvalue.append(int(header[par_begin+ilis]))
					ilis += 1
				elif ink[il]==0:
					nvalue.append(float(header[par_begin+ilis]))
					ilis += 1
				else:
					assert ink[il]==2
					t = Transform()
					tmp = []
					for iii in xrange(par_begin+ilis, par_begin+ilis+12):
						tmp.append(float(header[iii]))
					t.set_matrix(tmp)
					ilis += 12
					nvalue.append(t)
			ISID = list_params.count('ID')
			if(ISID == 0):
				imm = im
			else:
				imm = nvalue[ISID]
			# read head, set params, and write it
			dummy = EMData()
			dummy.read_image(stack, imm, True)
			set_arb_params(dummy, nvalue, list_params)
			dummy.write_image(stack, dummy.get_attr_default('ID', im), EMUtil.ImageType.IMAGE_HDF, True)

def send_attr_dict(main_node, data, list_params, image_start, image_end, comm = -1):
	import types
	from utilities import get_arb_params
	from mpi 	   import mpi_send
	from mpi 	   import MPI_FLOAT, MPI_INT, MPI_TAG_UB, MPI_COMM_WORLD
	#  This function is called from a node other than the main node

	if comm == -1: comm = MPI_COMM_WORLD
	TransType = type(Transform())
	mpi_send([image_start, image_end], 2, MPI_INT, main_node, MPI_TAG_UB, comm)
	nvalue = []
	for im in xrange(image_start, image_end):
		value = get_arb_params(data[im-image_start], list_params)
		for il in xrange(len(value)):
			if    type(value[il]) is types.IntType:  nvalue.append(float(value[il]))
			elif  type(value[il]) is types.FloatType: nvalue.append(value[il])
			elif  type(value[il]) is TransType: 
				m = value[il].get_matrix()
				assert (len(m)==12)
				for f in m: nvalue.append(f)
	mpi_send(nvalue, len(nvalue), MPI_FLOAT, main_node, MPI_TAG_UB, comm)

def recv_attr_dict_bdb(main_node, stack, data, list_params, image_start, image_end, number_of_proc, comm = -1):
	import types
	from  utilities import  get_arb_params, set_arb_params
	from  mpi 	import mpi_recv
	from  mpi 	import MPI_FLOAT, MPI_INT, MPI_TAG_UB, MPI_COMM_WORLD
	#  bdb version!
	# This is done on the main node, so for images from the main node, simply write headers
	
	if comm == -1: comm = MPI_COMM_WORLD
	
	DB = db_open_dict(stack)
	TransType = type(Transform())
	# prepare keys for float/int
	value = get_arb_params(data[0], list_params)
	ink = []
	len_list = 0
	ISID = -1
	for il in xrange(len(list_params)):
		if(list_params[il] == 'ID'):  ISID = il
		if type(value[il]) is types.IntType:
			ink.append(1)
			len_list += 1
		elif type(value[il]) is types.FloatType:
			ink.append(0)
			len_list += 1
		elif type(value[il]) is TransType:
			ink.append(2)
			len_list += 12
	ldis = []
	headers = []
	for n in xrange(number_of_proc):
		if n != main_node:
			dis = mpi_recv(2, MPI_INT, n, MPI_TAG_UB, comm)
			img_begin = int(dis[0])
			img_end = int(dis[1])
			header = mpi_recv(len_list*(img_end-img_begin), MPI_FLOAT, n, MPI_TAG_UB, comm)
			for im in xrange(img_begin, img_end):
				par_begin = (im-img_begin)*len_list
				nvalue = []
				ilis = 0
				for il in xrange(len(list_params)):
					if(ink[il] == 1):
						nvalue.append(int(header[par_begin+ilis]))
						ilis += 1
					elif ink[il]==0:
						nvalue.append(float(header[par_begin+ilis]))
						ilis += 1
					else:
						assert ink[il]==2
						t = Transform()
						tmp = []
						for iii in xrange(par_begin+ilis, par_begin+ilis+12):
							tmp.append(float(header[iii]))
						t.set_matrix(tmp)
						ilis += 12
						nvalue.append(t)
				if(ISID == -1):
					imm = im
				else:
					imm = nvalue[ISID]
				for i in xrange(len(list_params)):
					if(list_params[i] != "ID"):  DB.set_attr(imm, list_params[i], nvalue[i])
		else:
			for n in xrange(image_start, image_end):
				ID = data[n-image_start].get_attr_default('ID', n)
				for param in list_params:
					if(param != "ID"):  DB.set_attr(ID, param, data[n-image_start].get_attr(param))
	DB.close()

def check_attr(ima, num, params, default_value, action="Warning"):
	from sys import exit
	attr_list = ima.get_attr_dict()	
	if attr_list.has_key(params) == False:
		if action=="Warning":
			print "WARNING: In image %i, cannot find attribute \'%s\' in the header, set it to the default value" %(num, params), default_value
			ima.set_attr_dict({params:default_value})
		elif action=="Error":
			print "ERROR:   In image %i, cannot find attribute \'%s\' in the header, the program has to terminate" %(num, params)
			exit()
		return False
	else: return True

def print_begin_msg(program_name, onscreen=False):
	from time import localtime, strftime
	t = 100
	stars = '*'*t
	string = "Beginning of the program " + program_name + ": " + strftime("%a, %d %b %Y %H:%M:%S", localtime())
	s = (t-len(string))/2
	spacing = ' '*s
	if onscreen:
		print stars
		print spacing+string
		print stars
	else:
		print_msg(stars+"\n")
		print_msg(spacing+string+"\n")
		print_msg(stars+"\n")

def print_end_msg(program_name, onscreen=False):
	from time import localtime, strftime
	t = 100
	stars = '*'*t
	string = "End of the program " + program_name + ": " + strftime("%a, %d %b %Y %H:%M:%S", localtime())
	s = (t-len(string))/2
	spacing = ' '*s
	if onscreen:
		print stars
		print spacing+string
		print stars
	else:
		print_msg(stars+"\n")
		print_msg(spacing+string+"\n")
		print_msg(stars+"\n")

def print_msg(msg):
	import sys
	import global_def
	if (global_def.IS_LOGFILE_OPEN == False):
		global_def.LOGFILE_HANDLE = open(LOGFILE,"w")
		global_def.IS_LOGFILE_OPEN = True
	if (global_def.BATCH):
		global_def.LOGFILE_HANDLE.write(msg)		
	else:
		sys.stdout.write(msg)
		global_def.LOGFILE_HANDLE.write(msg)
	global_def.LOGFILE_HANDLE.flush()

def read_fsc( filename ):
	from string import split, atof
	f = open( filename, 'r' )
	fscc = None
	line = f.readline()
	while len(line) > 0:
		items = split( line )
		if fscc is None:
			fscc = [None]*len(items)
			for i in xrange( len(items) ):
				fscc[i] = []

		for i in xrange( len(items) ) :
			fscc[i].append( atof(items[i]) )

		line = f.readline()

	return fscc
"""
#  This would not work on windows
def memory_usage():
	import os
	from string import split
	return 0
	file = "/proc/%d/status" % os.getpid()
	f = open(file, 'r')
	line = f.readline()
	while len(line) > 0 :
		items = split( line )
		if items[0]=='VmSize:':
			return items[1]+items[2]
		line = f.readline()
"""

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

	[mean_a,sigma,imin,imax] = Util.infomask(img, model_circle(outer, nx, ny, nz) - inner_sphere, True)
	inner_rest   = model_blank(nx, ny, nz, 1.0) - inner_sphere
	Util.mul_img(inner_sphere, img)
	return Util.addn_img(inner_sphere, Util.mult_scalar(inner_rest, mean_a ) )

def copy_attr( pin, name, pot ):
	pot.set_attr( name, pin.get_attr(name) )
	pass

def write_headers(filename, data, lima):
	"""
	  write headers from files in data into a disk file called filename.
	  The filename has to be either hdf or bdb.
	  lima - list with positions in the disk files into which headers will be written,
	    i.e., header from data[k] will be written into file number lima[k]
	  WARNING: this function will open and close DB library!
	"""
	from utilities import file_type
	ftp = file_type(filename)
	if ftp == "bdb":
		#  For unknown reasons this does not work on Linux, but works on Mac ??? Really?
		DB = db_open_dict(filename)
		for i in range(len(lima)):
			DB.set_header(lima[i], data[i])
		DB.close()
		#for i in range(len(lima)):
		#	data[i].write_image(filename, lima[i])
	elif ftp == "hdf":
		for i in range(len(lima)):
			data[i].write_image(filename, lima[i], EMUtil.ImageType.IMAGE_HDF, True)
	else:
		ERROR("Unacceptable file format","write_headers",1)

def write_header(filename, data, lima):
	"""
	  write header from a single file data into a disk file called filename.
	  The filename has to be either hdf or bdb.
	  lima - position in the disk files into which header will be written,
	    i.e., header from data will be written into file number lima
	  WARNING: this function assums DB library is opened and will NOT close it!
	"""
	from utilities import file_type
	ftp = file_type(filename)
	if ftp == "bdb":
		DB = db_open_dict(filename)
		DB.set_header(lima, data)
	elif ftp == "hdf":
		data.write_image(filename, lima, EMUtil.ImageType.IMAGE_HDF, True)
	else:
		ERROR("Unacceptable file format","write_headers",1)

def file_type(name):
	if(len(name)>4):
		if(name[:4] == "bdb:"): return "bdb"
		elif(name[-4:-3] == "."):  return name[-3:]
	ERROR("Unacceptable file format","file_type",1)

def get_params2D(ima, xform = "xform.align2d"):
	"""
	  retrieve 2D alignment parameters from the header
	  alpha tx ty mirror scale
	"""
	t = ima.get_attr(xform)
	d = t.get_params("2D")
	return d["alpha"],d["tx"],d["ty"],d["mirror"],d["scale"]

def set_params2D(ima, p, xform = "xform.align2d"):
	"""
	  set 2D alignment parameters in the header
	  alpha tx ty mirror scale
	"""
	t = Transform({"type":"2D","alpha":p[0],"tx":p[1],"ty":p[2],"mirror":p[3],"scale":p[4]})
	ima.set_attr(xform, t)

def get_params3D(ima, xform = "xform.align3d"):
	"""
	  retrieve 3D alignment parameters from the header
	  phi  theta  psi  tx  ty  tz mirror scale
	"""
	t = ima.get_attr(xform)
	d = t.get_params("spider")
	return  d["phi"],d["theta"],d["psi"],d["tx"],d["ty"],d["tz"],d["mirror"],d["scale"]

def set_params3D(ima, p, xform = "xform.align3d"):
	"""
	  set 3D alignment parameters in the header
	  phi  theta  psi  tx  ty  tz mirror scale
	"""
	t = Transform({"type":"spider","phi":p[0],"theta":p[1],"psi":p[2],"tx":p[3],"ty":p[4],"tz":p[5],"mirror":p[6],"scale":p[7]})
	ima.set_attr(xform, t)

def get_params_proj(ima, xform = "xform.projection"):
	"""
	  retrieve projection alignment parameters from the header
	  phi  theta  psi  s2x  s2y
	"""
	t = ima.get_attr(xform)
	d = t.get_params("spider")
	return  d["phi"],d["theta"],d["psi"],-d["tx"],-d["ty"]

def set_params_proj(ima, p, xform = "xform.projection"):
	"""
	  set projection alignment parameters in the header
	  phi  theta  psi  s2x  s2y
	"""
	t = Transform({"type":"spider","phi":p[0],"theta":p[1],"psi":p[2]})
	t.set_trans(Vec2f(-p[3], -p[4]))
	ima.set_attr(xform, t)

def get_ctf(ima):
	"""
	  recover numerical values of CTF parameters from EMAN2 CTF object stored in a header of the input image
	  order of returned parameters:
        [defocus, cs, voltage, apix, bfactor, ampcont]
	"""

	from EMAN2 import EMAN2Ctf
	
	ctf_params = ima.get_attr("ctf")
	
	return ctf_params.defocus, ctf_params.cs, ctf_params.voltage, ctf_params.apix, ctf_params.bfactor, ctf_params.ampcont

def generate_ctf(p):
	"""
	  generate EMAN2 CTF object using values of CTF parameters given in the list p
	  order of parameters:
        [defocus, cs, voltage, apix, bfactor, ampcont]
	[ microns, mm, kV, Angstroms, A^2, %]
	"""
	from EMAN2 import EMAN2Ctf

	defocus = p[0]
	cs = p[1]
	voltage = p[2]
	pixel_size = p[3]
	bfactor = p[4]
	amp_contrast = p[5]
	
	if defocus > 100:  # which means it is very likely in Angstrom, therefore we are using the old convention
		defocus *= 1e-4
	
	if amp_contrast < 1.0:
		from math import sqrt
		amp_contrast = amp_contrast*100/sqrt(2*amp_contrast**2-2*amp_contrast+1)

	ctf = EMAN2Ctf()
	ctf.from_dict({"defocus":defocus, "cs":cs, "voltage":voltage, "apix":pixel_size, "bfactor":bfactor, "ampcont":amp_contrast})
	return ctf

def set_ctf(ima, p):
	"""
	  set EMAN2 CTF object in the header of input image using values of CTF parameters given in the list p
	  order of parameters:
        [defocus, cs, voltage, apix, bfactor, ampcont]
	"""
	from utilities import generate_ctf
	ctf = generate_ctf( p )
	ima.set_attr( "ctf", ctf )

def delete_bdb(name):
	"""
	  Delete bdb stack
	"""
	a = db_open_dict(name)
	db_remove_dict(name)


# parse user function parses the --function option. this option
#    can be either a single function name (i.e. --function=ali3d_e)
#    or a list of names, specifying the location of a user-defined
#    function; this will have the form --function=/path/module/function

def parse_user_function(opt_string):

    # check if option string is a string and return None if not. this
    #    will cause the user function to be set to default value
    #    "ref_ali3d" in the ali functions....
    if not(type(opt_string) is str):
        return None

    # check opt_string for format:
    if (opt_string.startswith("[") and opt_string.endswith("]")):
        # options string is [path,file,function]
        opt_list = opt_string[1:-1].split(",")
        if (2 == len(opt_list)):
            # options are [file,function]
            return [opt_list[0],opt_list[1]]
        elif (3 == len(opt_list)):
            # options are [path,file,function]. note the order!
            return [opt_list[1],opt_list[2],opt_list[0]]
        else:
            # neither. assume this is an error and return default
            return None
    else:
        # no list format used, so we assume this is a function name
        #    defined (and referenced) in user_functions.
        return opt_string


def getvec( phi, tht ):
	from math import pi,cos,sin
	angle_to_rad = pi/180.0

	if tht > 180.0:
		tht -= 180.0
		phi += 180.0
	elif tht > 90.0:
		tht = 180.0 - tht
		phi += 180.0

	assert tht <=90.0

	tht *= angle_to_rad
	phi *= angle_to_rad

	x = sin(tht)*cos(phi) 
	y = sin(tht)*sin(phi)
	z = cos(tht)

	return (x,y,z)

def nearest_ang( vecs, phi, tht ) :
	from utilities import getvec
	vec = getvec( phi, tht )

	best_s = -1.0
	best_i = -1

	for i in xrange( len(vecs) ):
		s = abs(vecs[i][0]*vec[0] + vecs[i][1]*vec[1] + vecs[i][2]*vec[2])
		if s > best_s:
			best_s = s
			best_i = i

	return best_i

def assign_projangles(projangles, refangles):
	refnormal = [None]*len(refangles)
	for i in xrange(len(refangles)):
		refnormal[i] = getvec( refangles[i][0], refangles[i][1] )
	assignments = [[] for i in xrange(len(refangles))]
	for i in xrange(len(projangles)):
		mm = nearest_ang( refnormal, projangles[i][0], projangles[i][1] )
		assignments[mm].append(i)
	return assignments

def cone_ang( projangles, phi, tht, ant ) :
	from utilities import getvec
	from math import cos, pi
	vec = getvec( phi, tht )

	cone = cos(ant*pi/180.0)
	la = []
	for i in xrange( len(projangles) ):
		vec  = getvec( phi, tht )
		vecs = getvec( projangles[i][0], projangles[i][1] )
		s = abs(vecs[0]*vec[0] + vecs[1]*vec[1] + vecs[2]*vec[2])
		if s >= cone:
			la.append(i)

	return la

def findall(lo,val):
	"""
	  Find all occurences of val on list lo
	  Returns a list of indices of val on lo.
	"""
	u = []
	i = -1
	while( i < len(lo)-1):
		try:
			i = lo.index(val,i+1)
			u.append(i)
		except:
			i += 1
	return  u

def disable_bdb_cache():

	import EMAN2db

	EMAN2db.BDB_CACHE_DISABLE = True

def enable_bdb_cache():

	import EMAN2db
	EMAN2db.BDB_CACHE_DISABLE = False

def helical_consistency(p2i, p1):
	"""
	  Find overall phi angle and z shift difference between two sets of projection parameters for helical structure.
	  The two sets have to be of the same length and it is assume that k'th element on the first
	  list corresponds to the k'th element on the second list.
	  Input: two lists [ [phi2, theta2, psi2, sx2, sy2], [phi1, theta1, psi1, sx1, sy1], ...].  Second list is considered reference.
	  	 parametes for helical symmetry-- dp, pixel_size, dphi
	  Output: , 3D_error between the two sets after adjustment 
	  Note: all angles have to be in spider convention.
	  """
	from pixel_error import angle_diff
	from math import cos,pi
	from utilities import getvec
	from pixel_error import angle_error
	n =len(p1[0])
	print n
	qtm = -1.0e10
	for lf in xrange(0,181,180):
		p2 = []
		p2.extend(p2i)
		if( lf == 180):
			tflip = Transform({"type":"spider","theta":180.0})
			for j in xrange(n):
				t2 = Transform({"type":"spider","phi":p2[0][j],"theta":p2[1][j],"psi":p2[2][j]})
				t2.set_trans( Vec2f( -p2[3][j], -p2[4][j] ) )
				t2 = t2*tflip
				d = t2.get_params("spider")
				p2[0][j] = d["phi"]
				p2[1][j] = d["theta"]
				p2[2][j] = d["psi"]
				p2[3][j] = -d["tx"]
				p2[4][j] = -d["ty"]
		tt1 = [0.0]*n
		tt2 = [0.0]*n
		mirror = [False]*n
		ln = 0
		for j in xrange( n ):
			t1 = getvec(p1[0][j],p1[1][j])
			t2 = getvec(p2[0][j],p2[1][j])
			tm = getvec(180.0+p2[0][j],180.0-p2[1][j])
			tt1[j] = t1[0]*t2[0]+t1[1]*t2[1]+t1[2]*t2[2]
			tt2[j] = t1[0]*tm[0]+t1[1]*tm[1]+t1[2]*tm[2]
			if(abs(tt1[j])<1.0e-7): tt1[j] = 0.0
			if(abs(tt2[j])<1.0e-7): tt2[j] = 0.0
			if(tt1[j]>tt2[j]):
				mirror[j] = True
				ln+=1
		print  " FLIP  ",lf
		if(ln < n//2):
			print "mirror  ",ln
			for j in xrange( n ):
				p2[0][j] += 180.0
				p2[1][j]  = 180.0-p2[1][j]
				p2[2][j]  = -p2[2][j]
				p2[4][j]  = -p2[4][j]
				mirror[j] = not(mirror[j])
		else:
			print " straight", ln
		phi1 = []
		phi2 = []
		agree = []
		for j in xrange(n):
			if(mirror[j]):
				phi1.append(p1[0][j])
				phi2.append(p2[0][j])
				agree.append(j)
		print  len(phi1)
		delta_phi = angle_diff( phi2, phi1 )
		print "close form diff===", delta_phi

		phi1 = []
		phi2 = []
		errorm = []
		for j in xrange( len( p1[0]) ):
			p2[0][j] = (p2[0][j] + delta_phi + 360)%360.0
			if(mirror[j]):
				phi1.append(p1[0][j])
				phi2.append(p2[0][j])
				errorm.append(angle_error( [ p2[0][j] ], [ p1[0][j] ]))
		qt = sum(errorm)/len(errorm)
		print  len(errorm),qt
		if(qt > qtm):
			qtm = qt
			p2o = []
			p2o.extend(p2)
			errormo = []
			phi1o = []
			phi2o = []
			errormo.extend(errorm)
			phi1o.extend(phi1)
			phi2o.extend(phi2)
		
	return p2o, errormo, agree, delta_phi, phi1o, phi2o

# according two lists of orientation or marker (phi, theta, psi for each one)
# return the global rotation (dphi, dtheta, dpsi) between the two systems
def rotation_between_anglesets(agls1, agls2):
	"""
	  Find overall 3D rotation (phi theta psi) between two sets of Eulerian angles.
	  The two sets have to be of the same length and it is assume that k'th element on the first
	  list corresponds to the k'th element on the second list.
	  Input: two lists [[phi1, theta1, psi1], [phi2, theta2, psi2], ...].  Second list is considered reference.
	  Output: overall rotation phi, theta, psi that has to be applied to the first list (agls1) so resulting
	    angles will agree with the second list.
	  Note: all angles have to be in spider convention.
	  For details see: Appendix in Penczek, P., Marko, M., Buttle, K. and Frank, J.:  Double-tilt electron tomography.  Ultramicroscopy 60:393-410, 1995.
	"""
	from math  import sin, cos, pi, sqrt, atan2, acos, atan
	from numpy import array, linalg, matrix
	import types

	rad2deg = 180.0 / pi
	deg2rad = 1.0 / rad2deg

	def ori2xyz(ori):
		if(type(ori) == types.ListType):
			phi, theta, psi = ori[:3]
		else:
			# it has to be Transformation object
			d = ori.get_params("spider")
			phi = d["phi"]
			theta = d["theta"]
			psi = d["psi"]

		mapt = False
		if theta > 90:
			theta = theta - 2 * (theta - 90.0)
			mapt   = True
		phi   *= deg2rad
		theta *= deg2rad
		x = sin(theta) * sin(phi)
		y = sin(theta) * cos(phi)
		val = 1.0 - x*x - y*y
		if val < 0.0: val = 0.0
		z = sqrt(val)
		if mapt: z = -z

		return [x, y, z]

	N = len(agls1)
	if N != len(agls2):
		print 'Both lists must have the same length'
		return -1
	if N < 2:
		print 'At least two orientations are required in each list'
		return -1
	U1, U2 = [], []
	for n in xrange(N):
		p1 = ori2xyz(agls1[n])
		p2 = ori2xyz(agls2[n])
		U1.append(p1)
		U2.append(p2)

	# compute all Suv with uv = {xx, xy, xz, yx, ..., zz}
	Suv   = [0] * 9
	c     = 0
	nbori = len(U1)
	for i in xrange(3):
		for j in xrange(3):
			for s in xrange(nbori):
				Suv[c] += (U2[s][i] * U1[s][j])
			c += 1

        # create matrix N
	N = array([[Suv[0]+Suv[4]+Suv[8], Suv[5]-Suv[7],        Suv[6]-Suv[2],                 Suv[1]-Suv[3]], 
		   [Suv[5]-Suv[7],        Suv[0]-Suv[4]-Suv[8], Suv[1]+Suv[3],                 Suv[6]+Suv[2]], 
		   [Suv[6]-Suv[2],        Suv[1]+Suv[3],        -Suv[0]+Suv[4]-Suv[8],         Suv[5]+Suv[7]],
		   [Suv[1]-Suv[3],        Suv[6]+Suv[2],        Suv[5]+Suv[7],         -Suv[0]-Suv[4]+Suv[8]]])

        # eigenvector corresponding to the most positive eigenvalue
	val, vec = linalg.eig(N)
	q0, qx, qy, qz = vec[:, val.argmax()]

        # create quaternion Rot matrix 
	r = [q0*q0-qx*qx+qy*qy-qz*qz,         2*(qy*qx+q0*qz),          2*(qy*qz-q0*qx),          0.0,
	     2*(qx*qy-q0*qz),                 q0*q0+qx*qx-qy*qy-qz*qz,  2*(qx*qz+q0*qy),          0.0,
	     2*(qz*qy+q0*qx),                 2*(qz*qx-q0*qy),          q0*q0-qx*qx-qy*qy+qz*qz,  0.0]
	
	R = Transform(r)
	dictR = R.get_rotation('SPIDER')

	return dictR['phi'], dictR['theta'], dictR['psi']


def get_pixel_size(img):
	"""
	  Retrieve pixel size from the header.
	  We check attribute Pixel_size and also pixel size from ctf object, if exisits.
	  If the two are different or if the pixel size is not set, return -1.0 and print a warning.
	"""
	p1 = img.get_attr_default("Pixel_size", -1.0)
	cc = img.get_attr_default("ctf", None)
	if cc == None:
		p2 = -1.0
	else:
		p2 = round(cc.apix, 3)
	if p1 == -1.0 and p2 == -1.0:
		ERROR("Pixel size not set", "get_pixel_size", 0)
		return -1.0
	elif p1 > -1.0 and p2 > -1.0:
		if abs(p1-p2) >= 0.001:
			ERROR("Conflict between pixel size in attribute and in ctf object", "get_pixel_size", 0)
		# pixel size is positive, so what follows omits -1 problem
		return max(p1, p2)
	else:
		return max(p1, p2)

def set_pixel_size(img, pixel_size):
	"""
	  Set pixel size in the header.
	  Set attribute Pixel_size and also pixel size in ctf object, if exists.
	"""
	img.set_attr("Pixel_size", round(pixel_size, 3))
	cc = img.get_attr_default("ctf", None)
	if(cc):
		cc.apix = pixel_size
		img.set_attr("ctf", cc)
