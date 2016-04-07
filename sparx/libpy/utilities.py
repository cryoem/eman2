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
		3D alignment parameters (phi, theta, psi, s2x, s2y)
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
				#### <--------->
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

'''
def golden(func, args=(), brack=None, tol=1.e-4, full_output=0):
	""" Given a function of one-variable and a possible bracketing interval,
	return the minimum of the function isolated to a fractional precision of
	tol. A bracketing interval is a triple (a,b,c) where (a<b<c) and
	func(b) < func(a),func(c).  If bracket is two numbers then they are
	assumed to be a starting interval for a downhill bracket search
	(see bracketing)

	Uses analog of bisection method to decrease the bracketed interval.
	"""
	from utilities import bracketing
	if brack is None:
		xa,xb,xc,fa,fb,fc,funcalls = bracketing(func, args=args)
	elif len(brack) == 2:
		xa,xb,xc,fa,fb,fc,funcalls = bracketing(func, xa=brack[0], xb=brack[1], args=args)
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


def bracketing(func, xa=0.0, xb=1.0, args=(), grow_limit=110.0, maxiter=1000):
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
'''


def ce_fit(inp_image, ref_image, mask_image):
	""" Fit the histogram of the input image under mask with the reference image.
		    
	     Usage : ce_fit(inp_image,ref_image,mask_image):
	    	 A and B, number of iterations and the chi-square
	"""
	hist_res = Util.histc(ref_image, inp_image, mask_image)
	args = hist_res["args"]
	scale = hist_res["scale"]
	data = [hist_res['data'], inp_image, hist_res["ref_freq_bin"], mask_image, int(hist_res['size_img']), hist_res['hist_len']]
	res = amoeba(args, scale, hist_func, 1.e-4, 1.e-4, 500, data)
	resu = ["Final Parameter [A,B]:", res[0], "Final Chi-square :", -1*res[1], "Number of Iteration :", res[2]]
	corrected_image = inp_image*res[0][0] + res[0][1]
	result = [resu,"Corrected Image :",corrected_image]
	del data[:], args[:], scale[:]
	return result

def center_2D(image_to_be_centered, center_method = 1, searching_range = -1, Gauss_radius_inner = 2, Gauss_radius_outter = 7, self_defined_reference = None):
	"""
		Put an input image into image center (nx/2, ny/2) using method :
		1. phase_cog
		2. cross-correlate with Gaussian function
		3. cross-correlate with donut shape image
		4. cross-correlate with reference image provided by user
		5. cross-correlate with self-rotated average
		7. binarize at ave+sigma and cross-correlate with a circle
	        The function will return centered_image, and shifts
	"""
	from   utilities    import peak_search
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

	elif center_method == 7:
		from fundamentals import ccf, cyclic_shift
		from morphology   import binarize
		from utilities    import model_blank
		from EMAN2        import rsconvolution
		p = Util.infomask(image_to_be_centered,None,True)
		cc = binarize(rsconvolution(binarize(image_to_be_centered,p[0]+p[1]),model_blank(5,5,1,1.0/(5.0*5.0))),0.5)	
		c = ccf(cc, self_defined_reference)
		p = Util.infomask(c,None,True)[3]
		nx = c.get_xsize()
		ny = c .get_ysize()
		cx = nx//2
		cy = ny//2
		n = 0
		x=0
		y=0
		for i in xrange(nx):
			for j in xrange(ny):
				if c.get_value_at(i,j) == p :
					x+=(i-cx)
					y+=(j-cy)
					n+=1
		shiftx = x/n
		shifty = y/n
		if searching_range > 0 :
			if(abs(shiftx) > searching_range):  shiftx=0
			if(abs(shifty) > searching_range):  shifty=0
		return cyclic_shift(image_to_be_centered, -shiftx, -shifty), shiftx, shifty

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


def compose_transform2m(alpha1=0.0, sx1=0., sy1=0.0, mirror1=0, scale1=1.0, alpha2=0.0, sx2=0.0, sy2=0.0, mirror2=0, scale2=1.0):
	"""Print the composition of two transformations  T2*T1
		Here  if v's are vectors:   vnew = T2*T1 vold
		     with T1 described by alpha1, sx1, scale1 etc.

	    Usage: compose_transform2(alpha1,sx1,sy1,mirror1,scale1,alpha2,sx2,sy2,mirror2,scale2)
	       angles in degrees
	"""

	t1 = Transform({"type":"2D","alpha":alpha1,"tx":sx1,"ty":sy1,"mirror":mirror1,"scale":scale1})
	t2 = Transform({"type":"2D","alpha":alpha2,"tx":sx2,"ty":sy2,"mirror":mirror2,"scale":scale2})
	tt = t2*t1
	d = tt.get_params("2D")
	return d[ "alpha" ], d[ "tx" ], d[ "ty" ], int(d[ "mirror" ]+0.1), d[ "scale" ]

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
	  Combine 2D alignent parameters including mirror: tt = t2*t1
	"""

	t1 = Transform({"type":"2D","alpha":alpha1,"tx":sx1,"ty":sy1,"mirror":mirror1,"scale":1.0})
	t2 = Transform({"type":"2D","alpha":alpha2,"tx":sx2,"ty":sy2,"mirror":mirror2,"scale":1.0})
	tt = t2*t1
	d = tt.get_params("2D")
	return d[ "alpha" ], d[ "tx" ], d[ "ty" ], int(d[ "mirror" ]+0.1)

def inverse_transform2(alpha, tx = 0.0, ty = 0.0, mirror = 0):
	"""Returns the inverse of the 2d rot and trans matrix

	    Usage: nalpha, ntx, nty, mirror = inverse_transform2(alpha,tx,ty,mirror)
	"""

	t = Transform({"type":"2D","alpha":alpha,"tx":tx,"ty":ty,"mirror":mirror,"scale":1.0})
	t = t.inverse()
	t = t.get_params("2D")
	return t[ "alpha" ], t[ "tx" ], t[ "ty" ], int(t[ "mirror" ]+0.1)

def inverse_transform3(phi, theta=0.0, psi=0.0, tx=0.0, ty=0.0, tz=0.0, mirror = 0, scale=1.0):
	"""Returns the inverse of the 3d rot and trans matrix

	    Usage: nphi,ntheta,npsi,ntx,nty,ntz,nmirror,nscale = inverse_transform3(phi,theta,psi,tx,ty,tz,mirror,scale)
	       angles in degrees
	"""

	d = Transform({'type': 'spider', 'phi': phi, 'theta': theta, 'psi': psi, 'tx': tx, 'ty': ty, 'tz': tz, "mirror":mirror,"scale":scale})
	d = d.inverse()
	d = d.get_params("spider")
	return  d["phi"],d["theta"],d["psi"],d["tx"],d["ty"],d["tz"],int(d["mirror"]+0.1),d["scale"]

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

	from math      import pi, sqrt, cos, acos, tan, sin
	from utilities import even_angles_cd
	from string    import lower,split
	angles = []
	symmetryLower = symmetry.lower()
	symmetry_string = split(symmetry)[0]
	if  (symmetry_string[0]  == "c"):
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
			for i in xrange(n):
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
		for i in xrange(n):
			t = n-i-1
			qt = angles[t][0]
			if((qt>=badb and qt<bade) or (qt>=bbdb and qt<bbde)):  del angles[t]
			
		if (int(symmetry_string[1:])%2 == 0):
			qt = 360.0/2/int(symmetry_string[1:])
		else:
			qt = 180.0/2/int(symmetry_string[1:])
		n = len(angles)
		for i in xrange(n):
			t = n-i-1
			if(angles[t][1] == 90.0):
				if(angles[t][0] >= qt + ant ):  del angles[t]
	elif(symmetry_string[0]  == "s"):
	
	#if symetry is "s", deltphi=delta, theata intial=theta1, theta end=90, delttheta=theta2
		# for helical, theta1 cannot be 0.0
		if theta1 > 90.0:
			ERROR('theta1 must be less than 90.0 for helical symmetry', 'even_angles', 1)
		if theta1 == 0.0: theta1 =90.0
		theta_number = int((90.0 - theta1)/theta2)
		#for helical, symmetry = s or scn
		cn = int(symmetry_string[2:])
		for j in xrange(theta_number,-1, -1):

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
					ERROR("For helical strucutre, we only support scn and sdn symmetry","even_angles",1)

			else:
				if (symmetry_string[1] =="c"):
					k=int(359.99/cn/delta)
				elif (symmetry_string[1] =="d"):
					k=int(359.99/2/cn/delta)
						
			for i in xrange(k+1):
					angles.append([i*delta,90.0-j*theta2,90.0])


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
	from EMAN2 import rsconvolution
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
		Return a NumPy array containing the image data. 
		Note: The NumPy array and the image data share the same memory,
		so if the NumPy array is altered then the image is altered
		as well (and vice versa).
	"""
	from EMAN2 import EMNumPy
	return EMNumPy.em2numpy(img)

def get_sym(symmetry):
	"""
	get a list of point-group symmetry angles, symmetry="c3"
	"""
	RA   = Transform()
	NTot = RA.get_nsym(symmetry)
	angs = []
	for j in xrange(NTot):
		RNow  = RA.get_sym(symmetry, j)
		RNowE = RNow.get_rotation('spider')
		angs.append([RNowE['phi'], RNowE['theta'], RNowE['psi']])

	return angs

def get_symt(symmetry):
	"""
	get a list of point-group symmetry transformations, symmetry="c3"
	"""

	RA   = Transform()
	NTot = RA.get_nsym(symmetry)
	angs = []
	for j in xrange(NTot):
		angs.append(RA.get_sym(symmetry, j))

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
		Extract input numbers from a given string
	"""
	from re import split
	qq = split(" |,",str_input)
	for i in xrange(len(qq)-1, -1, -1):
		if(qq[i] == ""):  del qq[i]
	o = []
	for i in xrange(len(qq)):
		if(qq[i].find(".") >= 0):  o.append(float(qq[i]))
		else:  o.append(int(qq[i]))
	return o

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

def image_decimate(img, decimation=2, fit_to_fft=1,frequency_low=0, frequency_high=0):
	from filter import filt_btwl
	from fundamentals import smallprime, window2d
	from utilities import get_image
	"""
		Window image to FFT-friendly size, apply Butterworth low pass filter,
		and decimate 2D image 
	"""
	if type(img)     == str :	img=get_image(img)
	if decimation    <= 1   :  	ERROR("Improper decimation ratio", "image_decimation", 1)
	if frequency_low <= 0   :	
		frequency_low  = .5/decimation- .05
		if frequency_low <= 0: ERROR("Butterworth passband frequency is too low", "image_decimation", 1)
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
	Create an image of a Gaussian function with standard deviation "xsigma,ysigma,zsigma"
	 and centered at (xcenter,ycenter,zcenter), by default the center is image center.
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
		if(line[11:12]==" " and line[8:10] != " "):			# new data format
			start= 13
			end  = 15
			#line_data.append(atoi(line[start:end]))		# 03/21/12 Anna: according to the previous version of this function this value was omitted
			start= end+3
			end  = start+6
			line_data.append(atof(line[start:end]))
			colNo = (len(line)-end)/12 - 1
			for i in xrange(colNo):
				start= end+6
				end  = start+7
				line_data.append(atof(line[start:end]))
			data.append(line_data)
			line = inf.readline()
		else:												# old data format
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
	import types
	outf = open(file_name, "w")
	if (type(data[0]) == types.ListType):
		# It is a list of lists
		for i in xrange(len(data)):
			for j in xrange(len(data[i])):
				qtp = type(data[i][j])
				if qtp == type(0):
					outf.write("  %12d"%data[i][j])
				elif qtp == type(0.0):
					# NOTE: 2015/05/27 Toshio Moriya
					# Since %12.5g does not work as the Python spec,
					# we manually mimic the spec.  
					# %12.5g results in int when a float value has only one non-zero digit after decimal point 
					a = data[i][j]
					z = (a>0.00001 and a<1.e6) or (a<-0.00001 and a>-1.e4)
					if z: outf.write("  %12.5f"%a)
					else: outf.write("  %12.5e"%a)
				else:
					outf.write("  %s"%data[i][j])
			outf.write("\n")
	else:
		# Single list
		for j in xrange(len(data)):
			qtp = type(data[j])
			if qtp == type(0):
				outf.write("  %12d"%data[j])
			elif qtp == type(0.0):
				# NOTE: 2015/05/27 Toshio Moriya
				# Since %12.5g does not work as the Python spec,
				# we manually mimic the spec.  
				# %12.5g results in int when a float value has only one non-zero digit after decimal point 
				a = data[j]
				z = (a>0.00001 and a<1.e6) or (a<-0.00001 and a>-1.e4)
				if z: outf.write("  %12.5f\n"%a)
				else: outf.write("  %12.5e\n"%a)
			else:
				outf.write("  %s"%data[j])
		outf.write("  \n")
	outf.flush()
	outf.close()


def read_text_file(file_name, ncol = 0):
	"""
		Read data from text file, if ncol = -1, read all columns
		if ncol >= 0, just read the (ncol)-th column.
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
					try:     data.append([int(vdata[i])])
					except:
						try:  	 data.append([float(vdata[i])])
						except:  data.append([vdata[i]])
			else:
				for i in xrange(len(vdata)):
					try:  data[i].append(int(vdata[i]))
					except:
						try:  data[i].append(float(vdata[i]))
						except:  data[i].append(vdata[i])
		else:
			vdata = split(line)[ncol]
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
	
	import types
	outf = open(file_name, "w")
	if (type(data[0]) == types.ListType):
		# It is a list of lists
		for i in xrange(len(data[0])):
			for j in xrange(len(data)):
				qtp = type(data[j][i])
				if qtp == type(0):
					outf.write("  %12d"%data[j][i])
				elif qtp == type(0.0):
					# NOTE: 2015/05/27 Toshio Moriya
					# Since %12.5g does not work as the Python spec,
					# we manually mimic the spec.  
					# %12.5g results in int when a float value has only one non-zero digit after decimal point 
					a = data[j][i]
					z = (a>0.00001 and a<1.e6) or (a<-0.00001 and a>-1.e4)
					if z: outf.write("  %12.5f"%a)
					else: outf.write("  %12.5e"%a)
				else:
					outf.write("  %s"%data[j][i])
			outf.write("\n")
	else:
		# Single list 
		for j in xrange(len(data)):
			qtp = type(data[j])
			if qtp == type(0):
				outf.write("  %12d\n"%data[j])
			elif qtp == type(0.0):
				# NOTE: 2015/05/27 Toshio Moriya
				# Since %12.5g does not work as the Python spec,
				# we manually mimic the spec.  
				# %12.5g results in int when a float value has only one non-zero digit after decimal point 
				a = data[j]
				z = (a>0.00001 and a<1.e6) or (a<-0.00001 and a>-1.e4)
				if z: outf.write("  %12.5f\n"%a)
				else: outf.write("  %12.5e\n"%a)
			else:
				outf.write("  %s\n"%data[j])
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

def rotate_shift_params(paramsin, transf):
	# moved from sxprocess.py
	if len(paramsin[0])>3 :
		from EMAN2 import Vec2f
		t = Transform({"type":"spider","phi":transf[0],"theta":transf[1],"psi":transf[2],"tx":transf[3],"ty":transf[4],"tz":transf[5],"mirror":0,"scale":1.0})
		t = t.inverse()
		cpar = []
		for params in paramsin:
			d = Transform({"type":"spider","phi":params[0],"theta":params[1],"psi":params[2]})
			d.set_trans(Vec2f(-params[3], -params[4]))
			c = d*t
			u = c.get_params("spider")
			cpar.append([u["phi"],u["theta"],u["psi"],-u["tx"],-u["ty"]])
	else:
		t = Transform({"type":"spider","phi":transf[0],"theta":transf[1],"psi":transf[2]})
		t = t.inverse()
		cpar = []
		for params in paramsin:
			d = Transform({"type":"spider","phi":params[0],"theta":params[1],"psi":params[2]})
			c = d*t
			u = c.get_params("spider")
			cpar.append([u["phi"],u["theta"],u["psi"]])
			#cpar.append([u["phi"],u["theta"],u["psi"],-u["tx"],-u["ty"]])
	return cpar


#  01/06/2016 - This is my recoding of old FORTRAN code with the hope that python's double precission
#                 will fix the problem of rotation of a 0,0,0 direction.  It does not as one neeeds psi
#                 in this case as well.  So, the only choice is to use small theta instead of exact 0,0,0 direction
def rotate_params(params, transf):
	matinv = rotmatrix( -transf[2], -transf[1], -transf[0] )
	n = len(params)
	cpar = [None]*n
	for i in xrange(n):
		d = rotmatrix( params[i][0], params[i][1], params[i][2] )
		c = mulmat(d,matinv)
		phi, theta, psi = recmat(c)
		cpar[i] = [phi, theta, psi]
	return cpar


def rotmatrix(phi,theta,psi):
	from math import sin,cos,radians
	rphi   = radians(phi)
	rtheta = radians(theta)
	rpsi   = radians(psi)
	mat = [[0.0]*3,[0.0]*3,[0.0]*3]

	mat[0][0] =  cos(rpsi)*cos(rtheta)*cos(rphi) - sin(rpsi)*sin(rphi)
	mat[1][0] = -sin(rpsi)*cos(rtheta)*cos(rphi) - cos(rpsi)*sin(rphi)
	mat[2][0] =            sin(rtheta)*cos(rphi)


	mat[0][1] =  cos(rpsi)*cos(rtheta)*sin(rphi) + sin(rpsi)*cos(rphi)
	mat[1][1] = -sin(rpsi)*cos(rtheta)*sin(rphi) + cos(rpsi)*cos(rphi)
	mat[2][1] =            sin(rtheta)*sin(rphi)


	mat[0][2] = -cos(rpsi)*sin(rtheta)
	mat[1][2] =  sin(rpsi)*sin(rtheta)
	mat[2][2] =            cos(rtheta)
	return mat

def mulmat(m1,m2):
	mat = [[0.0]*3,[0.0]*3,[0.0]*3]
	for i in xrange(3):
		for j in xrange(3):
			for k in xrange(3):
				mat[i][j] += m1[i][k]*m2[k][j]
			mat[i][j] = round(mat[i][j],8)
	return mat

def recmat(mat):
	from math import acos,asin,atan2,degrees,pi
	def sign(x):
		if( x >= 0.0 ): return 1
		else:  return -1
	"""
	mat = [[0.0]*3,[0.0]*3,[0.0]*3]
	# limit precision
	for i in xrange(3):
		for j in xrange(3):
			mat[i][j] = inmat[i][j]	
			#if(abs(inmat[i][j])<1.0e-8):  mat[i][j] = 0.0
			#else: mat[i][j] = inmat[i][j]
	for i in xrange(3):
		for j in xrange(3):  print  "     %14.8f"%mat[i][j],
		print ""
	"""
	if(mat[2][2] == 1.0):
		theta = 0.0
		psi = 0.0
		if( mat[0][0] == 0.0 ):
			phi = asin(mat[0][1])
		else:
			phi = atan2(mat[0][1],mat[0][0])
	elif(mat[2][2] == -1.0):
		theta = pi
		psi = 0.0
		if(mat[0][0] == 0.0):
			phi = asin(-mat[0][1])
		else:
			phi = atan2(-mat[0][1],-mat[0][0])
	else:
		theta = acos(mat[2][2])
		st = sign(theta)
		#print theta,st,mat[2][0]
		if(mat[2][0] == 0.0):
			if( st != sign(mat[2][1]) ):
				phi = 1.5*pi
			else:
				phi = 0.5*pi
		else:
			phi = atan2(st*mat[2][1], st*mat[2][0])

		#print theta,st,mat[0][2],mat[1][2]
		if(mat[0][2] == 0.0):
			if( st != sign(mat[1][2]) ):
				psi = 1.5*pi
			else:
				psi = 0.5*pi
		else:
			psi = atan2(st*mat[1][2], -st*mat[0][2])
	pi2 = 2*pi
	return  degrees(round(phi,10)%pi2),degrees(round(theta,10)%pi2),degrees(round(psi,10)%pi2)
	#return  degrees(round(phi,10)%pi2)%360.0,degrees(round(theta,10)%pi2)%360.0,degrees(round(psi,10)%pi2)%360.0
	#return  degrees(phi)%360.0,degrees(theta)%360.0,degrees(psi)%360.0

def reduce2asymmetric_D(angles_list, symmetry="d2"):
	sym_number =int(symmetry[1:])
	sym_angle  =360./sym_number
	badb      = 360.0/int(symmetry[1:])/4
	bade      = 2*badb
	bbdb      = badb + 360.0/int(symmetry[1:])/2
	bbde      = bbdb + 360.0/int(symmetry[1:])/4
	print badb, bade
	print bbdb, bbde
	alist =[]
	for index in xrange(len(angles_list)):
		if len(angles_list[index]) == 3 :
			[phi,theta,psi]= angles_list[index]
		else:
			[phi,theta,psi,sx,sy] = angles_list[index]
		if (phi > badb and phi < bade) or phi> bbdb: # only mirror those fall in non-permitted zones
			phi   = (phi+540.0)%360.0
			theta = 180.0-theta
			psi   = (540.0-psi)%360.0
		t2 = Transform({"type":"spider","phi":phi,"theta":theta,"psi":psi})
		ts = t2.get_sym_proj(symmetry)
		dlist ={}
		for item in xrange(len(ts)):
			a         = ts[item]
			u         = a.get_params("spider")
			phi       = u["phi"]
			theta     = u["theta"]
			psi       = u["psi"]
			#print qt, badb, bade, bbdb, bbde
			if theta <=90:
				if phi<=badb or (phi>=bade and phi<=bbdb):
					#print phi
					phi = phi%sym_angle
					phi = round(phi,3)
					dlist[phi] = [theta,psi]
		tlist = dlist.keys()
		if len(tlist) == 1 :
			alist.append([index,tlist[0], dlist[tlist[0]][0], dlist[tlist[0]][1]])
	return alist

def reduce2asymmetric_C(angles_list,symmetry="c1"):
	sym_angle  = 360.0/int(symmetry[1:])
	alist      = []
	for index in xrange(len(angles_list)):
		if len(angles_list[index]) == 3:
			[phi,theta,psi] = angles_list[index]
		else:
			[phi,theta,psi,sx,sy] = angles_list[index]
		phi = phi%sym_angle
		alist.append([phi,theta,psi])			
	return alist

def reduce_to_asymmetric_unit(angles_list,symmetry):
	from string import atoi
	sym_type   = symmetry[0:1].lower()
	if sym_type=="c":
		flist = reduce2asymmetric_C(angles_list, symmetry=symmetry)
		return flist
	elif sym_type=="d":
		new_angleslist1 = reduce2asymmetric_D(angles_list, symmetry=symmetry)
		tlist = {}
		for a in new_angleslist1:
			tlist[a[0]] = [a[1],a[2],a[3]]
		flist =[]
		for i in xrange(len(tlist)):
			flist.append(tlist[i])
		return flist


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
			ERROR("Incorrect input parameters","reshape_1d",1)
			return []

	if  Pixel_size_current == 0.:
		Pixel_size_current = 1.
		Pixel_size_interpolated = Pixel_size_current*float(length_current)/float(length_interpolated)
	qt =Pixel_size_interpolated/Pixel_size_current

	for i in xrange(length_interpolated):
		xi = float(i)*qt
		ix = min(int(xi),lt)
		df = xi -ix
		xval = (1.0-df)*input_object[ix] + df*input_object[ix+1]
		interpolated.append(xval)
	return interpolated

def estimate_3D_center(data):
	from math import cos, sin, radians
	from numpy import matrix
	from numpy import linalg
	import types
	if(type(data[0]) is types.ListType):
		ali_params = data
	else:
		ali_params = []
		for im in data:
			phi, theta, psi, s2x, s2y = get_params_proj(im)
			ali_params.append([phi, theta, psi, s2x, s2y])

	N = len(ali_params)
	A = []
	b = []
	
	for i in xrange(N):
		phi_rad   = radians(ali_params[i][0])
		theta_rad = radians(ali_params[i][1])
		psi_rad   = radians(ali_params[i][2])
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


def estimate_3D_center_MPI(data, nima, myid, number_of_proc, main_node, mpi_comm=None):
	from math import cos, sin, radians
	from numpy import matrix
	from numpy import linalg
	from mpi import MPI_COMM_WORLD
	from mpi import mpi_recv, mpi_send, MPI_FLOAT
	from applications import MPI_start_end

	if mpi_comm == None:
		mpi_comm = MPI_COMM_WORLD

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
				temp = mpi_recv(n_params, MPI_FLOAT, proc, proc, mpi_comm)
				for nn in xrange(n_params):
					ali_params_series.append(float(temp[nn]))

		ali_params = []
		N = len(ali_params_series)/5
		for im in xrange(N):
			ali_params.append([ali_params_series[im*5], ali_params_series[im*5+1], ali_params_series[im*5+2], ali_params_series[im*5+3], ali_params_series[im*5+4]])

		A = []
		b = []
	
		for i in xrange(N):
			phi_rad   = radians(ali_params[i][0])
			theta_rad = radians(ali_params[i][1])
			psi_rad   = radians(ali_params[i][2])
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
		mpi_send(ali_params_series, n_params, MPI_FLOAT, main_node, myid, mpi_comm)
		
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
	from numpy import shape, reshape
	from mpi   import mpi_reduce, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD, mpi_barrier
	
	if comm == -1 or comm == None:  comm = MPI_COMM_WORLD

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

def bcast_compacted_EMData_all_to_all(list_of_em_objects, myid, comm=-1):

	"""
	The assumption in <<bcast_compacted_EMData_all_to_all>> is that each processor
	calculates part of the list of elements and then each processor sends
	its results to the other ones. I

	Therefore, each processor has access to the header. If we assume that the
	attributes of interest from the header are the same for all elements then
	we can copy the header and no mpi message is necessary for the
	header.

	"""
	from applications import MPI_start_end
	from EMAN2 import EMNumPy
	from numpy import concatenate, shape, array, split
	from mpi import mpi_comm_size, mpi_bcast, MPI_FLOAT, MPI_COMM_WORLD
	from numpy import reshape

	if comm == -1 or comm == None: comm = MPI_COMM_WORLD

	num_ref = len(list_of_em_objects)
	ncpu = mpi_comm_size(comm)	# Total number of processes, passed by --np option.

	ref_start, ref_end = MPI_start_end(num_ref, ncpu, myid)

	for first_myid_process_that_has_em_elements in range(ncpu):
		sim_start, sim_ref_end = MPI_start_end(num_ref, ncpu, first_myid_process_that_has_em_elements)
		if sim_start != sim_ref_end:
			break
	else:
		raise ValueError("No processor contains em objects!")

	if myid == first_myid_process_that_has_em_elements:
		# used for copying the header and other info
		
		reference_em_object = list_of_em_objects[ref_start]
		data = EMNumPy.em2numpy(reference_em_object)
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
		print "Sending refrings: size of data to broadcast is greater than 2GB"

	for sender_id in range(ncpu):
		sender_ref_start, sender_ref_end = MPI_start_end(num_ref, ncpu, sender_id)		

		if sender_id == myid:
			if ref_start == ref_end:
				continue
			data = EMNumPy.em2numpy(list_of_em_objects[ref_start])  #array([], dtype = 'float32')
			for i in xrange(ref_start+1,ref_end):
				data = concatenate([data, EMNumPy.em2numpy(list_of_em_objects[i])])
		else:
			if sender_ref_start == sender_ref_end:
				continue
			data = array([], dtype = 'float32')

		sender_size_of_refrings = (sender_ref_end - sender_ref_start)*size_of_one_refring_assumed_common_to_all

		# size_of_refrings = mpi_bcast(size_of_refrings, 1, MPI_INT, sender_id, comm)
		data = mpi_bcast(data, sender_size_of_refrings, MPI_FLOAT, sender_id, comm)
		# print "Just sent %d float32 elements"%data.size

		if myid != sender_id:
			for i in xrange(sender_ref_start, sender_ref_end):
				offset_ring = sender_ref_start
				start_p = (i-offset_ring)*size_of_one_refring_assumed_common_to_all
				end_p   = (i+1-offset_ring)*size_of_one_refring_assumed_common_to_all
				image_data = data[start_p:end_p]

				if int(nz) != 1:
					image_data = reshape(image_data, (nz, ny, nx))
				elif ny != 1:
					image_data = reshape(image_data, (ny, nx))

				em_object = EMNumPy.numpy2em(image_data)
				em_object.set_attr_dict(em_dict)
				list_of_em_objects[i] = em_object

def bcast_compacted_EMData_all_to_all___original(list_of_em_objects, myid, comm=-1):

	"""
	The assumption in <<bcast_compacted_EMData_all_to_all>> is that each processor
	calculates part of the list of elements and then each processor sends
	its results to the other ones. I

	Therefore, each processor has access to the header. If we assume that the
	attributes of interest from the header are the same for all elements then
	we can copy the header and no mpi message is necessary for the
	header.

	"""
	from applications import MPI_start_end
	from EMAN2 import EMNumPy
	from numpy import concatenate, shape, array, split
	from mpi import mpi_comm_size, mpi_bcast, MPI_FLOAT, MPI_COMM_WORLD
	from numpy import reshape

	if comm == -1 or comm == None: comm = MPI_COMM_WORLD

	num_ref = len(list_of_em_objects)
	ncpu = mpi_comm_size(comm)	# Total number of processes, passed by --np option.

	ref_start, ref_end = MPI_start_end(num_ref, ncpu, myid)
	
	# used for copying the header
	reference_em_object = list_of_em_objects[ref_start]
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
	
	data = EMNumPy.em2numpy(list_of_em_objects[ref_start])
	size_of_one_refring_assumed_common_to_all = data.size
	
	# n = shape(data)
	# size_of_one_refring_assumed_common_to_all = 1
	# for i in n: size_of_one_refring_assumed_common_to_all *= i

	if size_of_one_refring_assumed_common_to_all*(ref_end-ref_start) > (2**31-1):
		print "Sending refrings: size of data to broadcast is greater than 2GB"

	for sender_id in range(ncpu):
		if sender_id == myid:
			data = EMNumPy.em2numpy(list_of_em_objects[ref_start])  #array([], dtype = 'float32')
			for i in xrange(ref_start+1,ref_end):
				data = concatenate([data, EMNumPy.em2numpy(list_of_em_objects[i])])
		else:
			data = array([], dtype = 'float32')

		sender_ref_start, sender_ref_end = MPI_start_end(num_ref, ncpu, sender_id)
		sender_size_of_refrings = (sender_ref_end - sender_ref_start)*size_of_one_refring_assumed_common_to_all

		# size_of_refrings = mpi_bcast(size_of_refrings, 1, MPI_INT, sender_id, comm)
		data = mpi_bcast(data, sender_size_of_refrings, MPI_FLOAT, sender_id, comm)
		# print "Just sent %d float32 elements"%data.size

		if myid != sender_id:
			for i in xrange(sender_ref_start, sender_ref_end):
				offset_ring = sender_ref_start
				start_p = (i-offset_ring)*size_of_one_refring_assumed_common_to_all
				end_p   = (i+1-offset_ring)*size_of_one_refring_assumed_common_to_all
				image_data = data[start_p:end_p]

				if int(nz) != 1:
					image_data = reshape(image_data, (nz, ny, nx))
				elif ny != 1:
					image_data = reshape(image_data, (ny, nx))
					
				em_object = EMNumPy.numpy2em(image_data)

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

				list_of_em_objects[i] = em_object



def gather_compacted_EMData_to_root_with_header_info_for_each_image(number_of_all_em_objects_distributed_across_processes, list_of_em_objects_for_myid_process, myid, comm=-1):

	"""
	
	The assumption in <<gather_compacted_EMData_to_root>> is that each processor
	calculates part of the list of elements and then each processor sends
	its results to the root

	Therefore, each processor has access to the header. If we assume that the
	attributes of interest from the header are the same for all elements then
	we can copy the header and no mpi message is necessary for the
	header.

	"""
	from applications import MPI_start_end
	from EMAN2 import EMNumPy
	from numpy import concatenate, shape, array, split
	from mpi import mpi_comm_size, mpi_bcast, MPI_FLOAT, MPI_COMM_WORLD
	from numpy import reshape

	if comm == -1 or comm == None: comm = MPI_COMM_WORLD

	ncpu = mpi_comm_size(comm)	# Total number of processes, passed by --np option.

	ref_start, ref_end = MPI_start_end(number_of_all_em_objects_distributed_across_processes, ncpu, myid)
	ref_end -= ref_start
	ref_start = 0
	
	# used for copying the header
	reference_em_object = list_of_em_objects_for_myid_process[ref_start]
	try:
		
		str_to_send = str(reference_em_object.get_attr_dict())

		# print "\nFFFFFFFFFFFFFF\n\n", reference_em_object.get_attr_dict(), "\n\n\n"
		# from mpi import mpi_finalize
		# mpi_finalize()
		# import sys
		# sys.stdout.flush()
		# sys.exit()
	except:
		raise ValueError("Could not convert em attribute dictionary to string s that gather_compacted_EMData_to_root_with_header_info_for_each_image can be used.")


	nx = reference_em_object.get_xsize()
	ny = reference_em_object.get_ysize()
	nz = reference_em_object.get_zsize()

	# # is_complex = reference_em_object.is_complex()
	# is_ri = reference_em_object.is_ri()
	# changecount = reference_em_object.get_attr("changecount")
	# is_complex_x = reference_em_object.is_complex_x()
	# is_complex_ri = reference_em_object.get_attr("is_complex_ri")
	# apix_x = reference_em_object.get_attr("apix_x")
	# apix_y = reference_em_object.get_attr("apix_y")
	# apix_z = reference_em_object.get_attr("apix_z")
	# is_complex = reference_em_object.get_attr_default("is_complex",1)
	# is_fftpad = reference_em_object.get_attr_default("is_fftpad",1)
	# is_fftodd = reference_em_object.get_attr_default("is_fftodd", nz%2)
	
	data = EMNumPy.em2numpy(list_of_em_objects_for_myid_process[ref_start])
	size_of_one_refring_assumed_common_to_all = data.size
	
	# n = shape(data)
	# size_of_one_refring_assumed_common_to_all = 1
	# for i in n: size_of_one_refring_assumed_common_to_all *= i

	if size_of_one_refring_assumed_common_to_all*(ref_end-ref_start) > (2**31-1):
		print "Sending refrings: size of data to broadcast is greater than 2GB"

	for sender_id in range(1,ncpu):
		if sender_id == myid:
			data = EMNumPy.em2numpy(list_of_em_objects_for_myid_process[ref_start])  #array([], dtype = 'float32')
			str_to_send = str(list_of_em_objects_for_myid_process[ref_start].get_attr_dict())
			em_dict_to_send_list = [list_of_em_objects_for_myid_process[ref_start].get_attr_dict()]
			for i in xrange(ref_start+1,ref_end):
				data = concatenate([data, EMNumPy.em2numpy(list_of_em_objects_for_myid_process[i])])
				# str_to_send += str(list_of_em_objects_for_myid_process[i].get_attr_dict())
				em_dict_to_send_list.append(list_of_em_objects_for_myid_process[i].get_attr_dict())
				
		else:
			data = array([], dtype = 'float32')

		sender_ref_start, sender_ref_end = MPI_start_end(number_of_all_em_objects_distributed_across_processes, ncpu, sender_id)

		sender_size_of_refrings = (sender_ref_end - sender_ref_start)*size_of_one_refring_assumed_common_to_all
		
		from mpi import mpi_recv, mpi_send, mpi_barrier
		if myid == 0:
			# print "root, receiving from ", sender_id, "  sender_size_of_refrings = ", sender_size_of_refrings
			str_to_receive = wrap_mpi_recv(sender_id)
			em_dict_list = eval(str_to_receive)
			# print "em_dict_list", em_dict_list
			data = mpi_recv(sender_size_of_refrings, MPI_FLOAT, sender_id, SPARX_MPI_TAG_UNIVERSAL, MPI_COMM_WORLD)
		elif sender_id == myid:
			wrap_mpi_send(str(em_dict_to_send_list), 0)
			# print "sender_id = ", sender_id, "sender_size_of_refrings = ", sender_size_of_refrings
			mpi_send(data, sender_size_of_refrings, MPI_FLOAT, 0, SPARX_MPI_TAG_UNIVERSAL, MPI_COMM_WORLD)
		
		mpi_barrier(MPI_COMM_WORLD)

		# if myid != sender_id:
		if myid == 0:
			for i in xrange(sender_ref_start, sender_ref_end):
				offset_ring = sender_ref_start
				start_p = (i-offset_ring)*size_of_one_refring_assumed_common_to_all
				end_p   = (i+1-offset_ring)*size_of_one_refring_assumed_common_to_all
				image_data = data[start_p:end_p]

				if int(nz) != 1:
					image_data = reshape(image_data, (nz, ny, nx))
				elif ny != 1:
					image_data = reshape(image_data, (ny, nx))
					
				em_object = EMNumPy.numpy2em(image_data)
				em_object.set_attr_dict(em_dict_list[i - sender_ref_start])

				# # em_object.set_complex(is_complex)
				# em_object.set_ri(is_ri)
				# em_object.set_attr_dict({
				# "changecount":changecount,
				# "is_complex_x":is_complex_x,
				# "is_complex_ri":is_complex_ri,
				# "apix_x":apix_x,
				# "apix_y":apix_y,
				# "apix_z":apix_z,
				# 'is_complex':is_complex,
				# 'is_fftodd':is_fftodd,
				# 'is_fftpad':is_fftpad})

				# list_of_em_objects[i] = em_object
				list_of_em_objects_for_myid_process.append(em_object)
				
		mpi_barrier(MPI_COMM_WORLD)

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
	from applications import MPI_start_end
	from EMAN2 import EMNumPy
	from numpy import concatenate, shape, array, split
	from mpi import mpi_comm_size, mpi_bcast, MPI_FLOAT, MPI_COMM_WORLD
	from numpy import reshape
	from mpi import mpi_recv, mpi_send, mpi_barrier

	if comm == -1 or comm == None: comm = MPI_COMM_WORLD

	ncpu = mpi_comm_size(comm)	# Total number of processes, passed by --np option.

	ref_start, ref_end = MPI_start_end(number_of_all_em_objects_distributed_across_processes, ncpu, myid)
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
	
	data = EMNumPy.em2numpy(list_of_em_objects_for_myid_process[ref_start])
	size_of_one_refring_assumed_common_to_all = data.size
	
	# n = shape(data)
	# size_of_one_refring_assumed_common_to_all = 1
	# for i in n: size_of_one_refring_assumed_common_to_all *= i

	if size_of_one_refring_assumed_common_to_all*(ref_end-ref_start) > (2**31-1):
		print "Sending refrings: size of data to broadcast is greater than 2GB"

	for sender_id in range(1,ncpu):
		if sender_id == myid:
			data = EMNumPy.em2numpy(list_of_em_objects_for_myid_process[ref_start])  #array([], dtype = 'float32')
			for i in xrange(ref_start+1,ref_end):
				data = concatenate([data, EMNumPy.em2numpy(list_of_em_objects_for_myid_process[i])])
		else:
			data = array([], dtype = 'float32')

		sender_ref_start, sender_ref_end = MPI_start_end(number_of_all_em_objects_distributed_across_processes, ncpu, sender_id)

		sender_size_of_refrings = (sender_ref_end - sender_ref_start)*size_of_one_refring_assumed_common_to_all
		
		if myid == 0:
			# print "root, receiving from ", sender_id, "  sender_size_of_refrings = ", sender_size_of_refrings
			data = mpi_recv(sender_size_of_refrings,MPI_FLOAT, sender_id, tag_for_send_receive, MPI_COMM_WORLD)
		elif sender_id == myid:
			# print "sender_id = ", sender_id, "sender_size_of_refrings = ", sender_size_of_refrings
			mpi_send(data, sender_size_of_refrings, MPI_FLOAT, 0, tag_for_send_receive, MPI_COMM_WORLD)
		
		mpi_barrier(MPI_COMM_WORLD)

		# if myid != sender_id:
		if myid == 0:
			for i in xrange(sender_ref_start, sender_ref_end):
				offset_ring = sender_ref_start
				start_p = (i-offset_ring)*size_of_one_refring_assumed_common_to_all
				end_p   = (i+1-offset_ring)*size_of_one_refring_assumed_common_to_all
				image_data = data[start_p:end_p]

				if int(nz) != 1:
					image_data = reshape(image_data, (nz, ny, nx))
				elif ny != 1:
					image_data = reshape(image_data, (ny, nx))
					
				em_object = EMNumPy.numpy2em(image_data)

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
				
		mpi_barrier(MPI_COMM_WORLD)


def bcast_EMData_to_all(tavg, myid, source_node = 0, comm = -1):
	from EMAN2 import EMNumPy
	from numpy import array, shape, reshape
	from mpi   import mpi_bcast, MPI_FLOAT, MPI_COMM_WORLD

	if comm == -1 or comm == None: comm = MPI_COMM_WORLD
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
	img_head.append(img.get_attr("changecount"))
	img_head.append(img.is_complex_x())
	img_head.append(img.get_attr("is_complex_ri"))
	img_head.append(int(img.get_attr("apix_x")*10000))
	img_head.append(int(img.get_attr("apix_y")*10000))
	img_head.append(int(img.get_attr("apix_z")*10000))

	head_tag = 2*tag
	mpi_send(img_head, 11, MPI_INT, dst, head_tag, comm)

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
	from EMAN2 import EMNumPy

	
	if comm==-1: comm = MPI_COMM_WORLD
	head_tag = 2*tag
	img_head = mpi_recv(11, MPI_INT, src, head_tag, comm)

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
	img.set_attr_dict({"changecount":int(img_head[5]),  "is_complex_x":int(img_head[6]),  "is_complex_ri":int(img_head[7]),  "apix_x":int(img_head[8])/10000.0,  "apix_y":int(img_head[9])/10000.0,  "apix_z":int(img_head[10])/10000.0})
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
											It is a dangerous assumption, it will have to be changed  07/10/2015
	"""
	from mpi import MPI_COMM_WORLD, MPI_INT
	from mpi import mpi_send, mpi_recv	

	l = len(data)
	gathered_data = []
	inc = 1   # A temp measure
	if myid == main_node:
		for i in xrange(0, number_of_proc*inc, inc):
			if i == main_node:
				for k in xrange(l):
					gathered_data.append(data[k])
			else:
				for k in xrange(l):
					im = recv_EMData(i, i*l+k)
					mem_len = mpi_recv(1, MPI_INT, i, SPARX_MPI_TAG_UNIVERSAL, MPI_COMM_WORLD)
					members = mpi_recv(int(mem_len[0]), MPI_INT, i, SPARX_MPI_TAG_UNIVERSAL, MPI_COMM_WORLD)
					members = map(int, members)
					im.set_attr('members', members)
					gathered_data.append(im)
	else:
		for k in xrange(l):
			send_EMData(data[k], main_node, myid*l+k)
			mem = data[k].get_attr('members')
			mpi_send(len(mem), 1, MPI_INT, main_node, SPARX_MPI_TAG_UNIVERSAL, MPI_COMM_WORLD)
			mpi_send(mem, len(mem), MPI_INT, main_node, SPARX_MPI_TAG_UNIVERSAL, MPI_COMM_WORLD)
	return gathered_data

def send_string_to_all(str_to_send, source_node = 0):
	from mpi import MPI_COMM_WORLD, MPI_INT, MPI_CHAR, mpi_bcast, mpi_comm_rank 

	myid = mpi_comm_rank(MPI_COMM_WORLD)
	str_to_send_len  = len(str_to_send)*int(myid == source_node)
	str_to_send_len = mpi_bcast(str_to_send_len,1,MPI_INT,source_node,MPI_COMM_WORLD)[0]
	str_to_send = mpi_bcast(str_to_send,str_to_send_len,MPI_CHAR,source_node,MPI_COMM_WORLD)
	return "".join(str_to_send)


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
		print  " ERROR in bcast_number_to_all"
	
def bcast_list_to_all(list_to_send, myid, source_node = 0):
	from mpi import mpi_bcast, MPI_COMM_WORLD, MPI_FLOAT, MPI_INT
	import   types
	if(myid == source_node):
		n = len(list_to_send)
		# we will also assume all elements on the list are of the same type
		if( type(list_to_send[0]) == types.IntType ): tp = 0
		elif( type(list_to_send[0]) == types.FloatType ): tp = 1
		else: tp = 2
	else:
		n = 0
		tp = 0
	n = bcast_number_to_all(n, source_node = source_node)
	tp = bcast_number_to_all(tp, source_node = source_node)
	if( tp == 2 ): 	ERROR("Only list of the same type numbers can be brodcasted","bcast_list_to_all",1, myid)
	if(myid != source_node): list_to_send = [0]*n

	if( tp == 0 ):
		list_to_send = mpi_bcast(list_to_send, n, MPI_INT, source_node, MPI_COMM_WORLD)
		return [int(n) for n in list_to_send]
	else:
		list_to_send = mpi_bcast(list_to_send, n, MPI_FLOAT, source_node, MPI_COMM_WORLD)
		return [float(n) for n in list_to_send]


def recv_attr_dict(main_node, stack, data, list_params, image_start, image_end, number_of_proc, comm = -1):
	import types
	from  utilities import  get_arb_params, set_arb_params
	from  mpi 	import mpi_recv
	from  mpi 	import MPI_FLOAT, MPI_INT, MPI_COMM_WORLD

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
			dis = mpi_recv(2, MPI_INT, n, SPARX_MPI_TAG_UNIVERSAL, comm)
			value = mpi_recv(len_list*(dis[1]-dis[0]), MPI_FLOAT, n, SPARX_MPI_TAG_UNIVERSAL, comm)
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
				imm = nvalue[list_params.index('ID')]
			# read head, set params, and write it
			dummy = EMData()
			dummy.read_image(stack, imm, True)
			set_arb_params(dummy, nvalue, list_params)
			dummy.write_image(stack, dummy.get_attr_default('ID', im), EMUtil.ImageType.IMAGE_HDF, True)

def send_attr_dict(main_node, data, list_params, image_start, image_end, comm = -1):
	import types
	from utilities import get_arb_params
	from mpi 	   import mpi_send
	from mpi 	   import MPI_FLOAT, MPI_INT, MPI_COMM_WORLD

	#  This function is called from a node other than the main node

	if comm == -1: comm = MPI_COMM_WORLD
	TransType = type(Transform())
	mpi_send([image_start, image_end], 2, MPI_INT, main_node, SPARX_MPI_TAG_UNIVERSAL, comm)
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
	mpi_send(nvalue, len(nvalue), MPI_FLOAT, main_node, SPARX_MPI_TAG_UNIVERSAL, comm)

def recv_attr_dict_bdb(main_node, stack, data, list_params, image_start, image_end, number_of_proc, comm = -1):
	import types
	from  utilities import  get_arb_params, set_arb_params
	from  mpi 	import mpi_recv
	from  mpi 	import MPI_FLOAT, MPI_INT, MPI_COMM_WORLD
	from EMAN2db import db_open_dict
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
			dis = mpi_recv(2, MPI_INT, n, SPARX_MPI_TAG_UNIVERSAL, comm)
			img_begin = int(dis[0])
			img_end = int(dis[1])
			header = mpi_recv(len_list*(img_end-img_begin), MPI_FLOAT, n, SPARX_MPI_TAG_UNIVERSAL, comm)
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
		global_def.LOGFILE_HANDLE = open(global_def.LOGFILE,"w")
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
	from EMAN2db import db_open_dict

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
	from EMAN2db import db_open_dict

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
	d = Util.get_transform_params(ima, xform, "2D")
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
	d = Util.get_transform_params(ima, xform, "spider")
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
	d = Util.get_transform_params(ima, xform, "spider")
	return  d["phi"],d["theta"],d["psi"],-d["tx"],-d["ty"]

def set_params_proj(ima, p, xform = "xform.projection"):
	"""
	  set projection alignment parameters in the header
	  phi  theta  psi  s2x  s2y
	"""
	from EMAN2 import Vec2f
	t = Transform({"type":"spider","phi":p[0],"theta":p[1],"psi":p[2]})
	t.set_trans(Vec2f(-p[3], -p[4]))
	ima.set_attr(xform, t)


def get_ctf(ima):
	"""
	  recover numerical values of CTF parameters from EMAN2 CTF object stored in a header of the input image
	  order of returned parameters:
        [defocus, cs, voltage, apix, bfactor, ampcont, astigmatism amplitude, astigmatism angle]
	"""
	from EMAN2 import EMAN2Ctf
	ctf_params = ima.get_attr("ctf")	
	return ctf_params.defocus, ctf_params.cs, ctf_params.voltage, ctf_params.apix, ctf_params.bfactor, ctf_params.ampcont, ctf_params.dfdiff, ctf_params.dfang

def generate_ctf(p):
	"""
	  generate EMAN2 CTF object using values of CTF parameters given in the list p
	  order of parameters:
        [defocus, cs, voltage, apix, bfactor, ampcont, astigmatism_amplitude, astigmatism_angle]
	    [ microns, mm, kV, Angstroms, A^2, microns, radians]
	"""
	from EMAN2 import EMAN2Ctf

	defocus      = p[0]
	cs           = p[1]
	voltage      = p[2]
	pixel_size   = p[3]
	bfactor      = p[4]
	amp_contrast = p[5]
	
	if defocus > 100:  # which means it is very likely in Angstrom, therefore we are using the old convention
		defocus *= 1e-4
	
	if amp_contrast < 1.0:
		from math import sqrt
		amp_contrast = amp_contrast*100/sqrt(2*amp_contrast**2-2*amp_contrast+1)

	ctf = EMAN2Ctf()
	if(len(p) == 6):
		ctf.from_dict({"defocus":defocus, "cs":cs, "voltage":voltage, "apix":pixel_size, "bfactor":bfactor, "ampcont":amp_contrast})
	else:
		ctf.from_dict({"defocus":defocus, "cs":cs, "voltage":voltage, "apix":pixel_size, "bfactor":bfactor, "ampcont":amp_contrast,'dfdiff':p[6],'dfang':p[7]})

	return ctf

def set_ctf(ima, p):
	"""
	  set EMAN2 CTF object in the header of input image using values of CTF parameters given in the list p
	  order of parameters:
        [defocus, cs, voltage, apix, bfactor, ampcont, astigmatism amplitude, astigmatism angle]
	"""
	from utilities import generate_ctf
	ctf = generate_ctf( p )
	ima.set_attr( "ctf", ctf )

def delete_bdb(name):
	"""
	  Delete bdb stack
	"""
	from EMAN2db import db_open_dict, db_remove_dict
	a = db_open_dict(name)
	db_remove_dict(name)


# parse user function parses the --function option. this option
#    can be either a single function name (i.e. --function=ali3d_e)
#    or a list of names, specifying the location of a user-defined
#    function; this will have the form --function=/path/module/function

def parse_user_function(opt_string):
	# check if option string is a string and return None if not. this
	# will cause the user function to be set to default value
	# "ref_ali3d" in the ali functions....
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
		# defined (and referenced) in user_functions.
		return opt_string

def getang(n):
	from math import atan2, acos, degrees
	return degrees(atan2(n[1],n[0]))%360.0, degrees(acos(n[2]))%360.0

def getang3(p1,p2):
	from utilities import getfvec
	from math import acos, degrees
	n1 = getfvec(p1[0],p1[1])
	n2 = getfvec(p2[0],p2[1])
	return degrees(acos(max(min((n1[0]*n2[0]+n1[1]*n2[1]+n1[2]*n2[2]),1.0),-1.0)))

def getvec( phi, tht ):
	from math import radians,cos,sin

	if tht > 180.0:
		tht -= 180.0
		phi += 180.0
	if tht > 90.0:
		tht = 180.0 - tht
		phi += 180.0

	assert tht <=90.0

	qt = radians(tht)
	qp = radians(phi)
	qs = sin(qt)

	x = qs*cos(qp) 
	y = qs*sin(qp)
	z = cos(qt)

	return (x,y,z)

def getfvec( phi, tht ):
	from math import radians,cos,sin
	qt = radians(tht)
	qp = radians(phi)
	qs = sin(qt)
	x = qs*cos(qp) 
	y = qs*sin(qp)
	z = cos(qt)

	return (x,y,z)

def nearest_ang( vecs, phi, tht ) :
	from utilities import getvec
	vec = getvec( phi, tht )
	return  Util.nearest_ang(vecs, vec[0],vec[1],vec[2])
	"""
	best_s = -1.0
	best_i = -1

	for i in xrange( len(vecs) ):
		s = abs(vecs[i][0]*vec[0] + vecs[i][1]*vec[1] + vecs[i][2]*vec[2])
		if s > best_s:
			best_s = s
			best_i = i

	return best_i
	"""
"""
Util.nearest_ang(vecs, vec[0],vec[1],vec[2])
def closest_ang( vecs, vec) :
	best_s = -1.0
	best_i = -1

	for i in xrange( len(vecs) ):
		s = abs(vecs[i][0]*vec[0] + vecs[i][1]*vec[1] + vecs[i][2]*vec[2])
		if s > best_s:
			best_s = s
			best_i = i

	return best_i
"""
# This is in python, it is very slow, we keep it just for comparison, use Util.assign_projangles instead
def assign_projangles_slow(projangles, refangles):
	refnormal = [None]*len(refangles)
	for i in xrange(len(refangles)):
		refnormal[i] = getvec(refangles[i][0], refangles[i][1])
	assignments = [[] for i in xrange(len(refangles))]
	for i in xrange(len(projangles)):
		best_i = nearest_ang(refnormal, projangles[i][0], projangles[i][1])
		assignments[best_i].append(i)
	return assignments

def nearestk_projangles(projangles, whichone = 0, howmany = 1, sym="c1"):
	# In both cases mirrored should be treated the same way as straight as they carry the same structural information
	lookup = range(len(projangles))
	if( sym == "c1"):
		from utilities import getvec
		refnormal = [None]*(len(projangles)*3)
		for i in xrange(len(projangles)):
			ref = getvec(projangles[i][0], projangles[i][1])
			for k in xrange(3):
				refnormal[3*i+k] = ref[k]
		# remove the reference projection from the list
		ref = [0.0,0.0,0.0]
		for k in xrange(3):
			ref[k] = refnormal[3*whichone+k]
		for k in xrange(3): del refnormal[3*whichone+2-k]
		del lookup[whichone]
		assignments = [-1]*howmany
		for i in xrange(howmany):
			k = Util.nearest_ang(refnormal, ref[0],ref[1],ref[2])
			assignments[i] = lookup[k]
			for l in xrange(3): del refnormal[3*k+2-l]
			del lookup[k]

	elif( sym[:1] == "d" ):
		from utilities import get_symt, getvec
		from EMAN2 import Vec2f, Transform
		t = get_symt(sym)
		phir = 360.0/int(sym[1:])
		for i in xrange(len(t)):  t[i] = t[i].inverse()
		a = Transform({"type":"spider","phi":projangles[whichone][0], "theta":projangles[whichone][1]})
		for l in xrange(len(t)):
			q = a*t[l]
			q = q.get_params("spider")
			if(q["phi"]<phir and q["theta"] <= 90.0): break
		refvec = getfvec(q["phi"], q["theta"])
		#print  "refvec   ",q["phi"], q["theta"]

		tempan =  [None]*len(projangles)
		for i in xrange(len(projangles)): tempan[i] = projangles[i]
		del tempan[whichone], lookup[whichone]
		assignments = [-1]*howmany

		for i in xrange(howmany):
			best = -1
			for j in xrange(len(tempan)):
				nearest = -1.
				a = Transform({"type":"spider","phi":tempan[j][0], "theta":tempan[j][1]})
				for l in xrange(len(t)):
					q = a*t[l]
					q = q.get_params("spider")
					vecs = getfvec(q["phi"], q["theta"])
					s = abs(vecs[0]*refvec[0] + vecs[1]*refvec[1] + vecs[2]*refvec[2])
					if( s > nearest ):
						nearest = s
						#ttt = (q["phi"], q["theta"])
				if( nearest > best ):
					best = nearest
					best_j = j
					#print  j,tempan[j][0], tempan[j][1],best,lookup[j],ttt
			assignments[i] = lookup[best_j]
			del tempan[best_j], lookup[best_j]

	elif( sym[:1] == "c" ):
		from utilities import get_symt, getvec
		from EMAN2 import Vec2f, Transform
		t = get_symt(sym)
		phir = 360.0/int(sym[1:])

		tempan =  [None]*len(projangles)
		for i in xrange(len(projangles)): tempan[i] = projangles[i]
		del tempan[whichone], lookup[whichone]
		assignments = [-1]*howmany

		for i in xrange(howmany):
			best = -1
			for j in xrange(len(tempan)):
				nearest = -1.
				a = Transform({"type":"spider","phi":tempan[j][0], "theta":tempan[j][1]})
				for l in xrange(len(t)):
					q = a*t[l]
					q = q.get_params("spider")
					vecs = getfvec(q["phi"], q["theta"])
					s = abs(vecs[0]*refvec[0] + vecs[1]*refvec[1] + vecs[2]*refvec[2])
					if( s > nearest ):
						nearest = s
						#ttt = (q["phi"], q["theta"])
				if( nearest > best ):
					best = nearest
					best_j = j
					#print  j,tempan[j][0], tempan[j][1],best,lookup[j],ttt
			assignments[i] = lookup[best_j]
			del tempan[best_j], lookup[best_j]

	else:
		print  "  ERROR:  symmetry not supported  ",sym
		assignments = []


	return assignments


def nearest_full_k_projangles(anormals, refang, howmany = 1, sym="c1"):
	# We assume refang can be on the list of normals
	from utilities import getfvec
	lookup = range(len(anormals))
	#refnormal = normals[:]
	assignments = [-1]*howmany

	if( sym == "c1"):
		refnormal = []
		for i,q in enumerate(anormals):
			refnormal += getfvec(q[0],q[1])
		ref = getfvec(refang[0],refang[1])
		for i in xrange(howmany):
			tmp = Util.nearest_fang(refnormal, ref[0],ref[1],ref[2])
			k = tmp[0]
			assignments[i] = lookup[k]
			for l in xrange(3): del refnormal[3*k+2-l]
			del lookup[k]

	elif( sym[:1] == "c" ):
		from utilities import get_symt, getfvec
		from EMAN2 import Vec2f, Transform
		phin = int(sym[1:])

		refnormal = []
		for i,q in enumerate(anormals):
			refnormal += getfvec(q[0]*phin,q[1])
		ref = getfvec(refang[0]*phin,refang[1])
		for i in xrange(howmany):
			tmp = Util.nearest_fang(refnormal, ref[0],ref[1],ref[2])
			k = tmp[0]
			assignments[i] = lookup[k]
			for l in xrange(3): del refnormal[3*k+2-l]
			del lookup[k]

	elif( sym[:1] == "d" ):
		from utilities import get_symt, getfvec
		from EMAN2 import Vec2f, Transform
		t = get_symt(sym)
		nt = len(t)
		a = Transform({"type":"spider","phi":refang[0], "theta":refang[1]})
		refvec = [None]*nt
		for i in xrange(nt):
			qt = a*(t[i].inverse())
			qt = qt.get_params("spider")
			refvec[i] = getfvec(qt["phi"], qt["theta"])
			print i,qt["phi"], qt["theta"],["psi"],refvec[i]

		refnormal = []
		for i,q in enumerate(anormals):
			refnormal += getfvec(q[0],q[1])

		for i in xrange(howmany):
			best_i = -1
			best_v = -10000000
			for l in xrange(nt):
				tmp = Util.nearest_fang(refnormal, refvec[l][0],refvec[l][1],refvec[l][2])
				if(tmp[1] > best_v):
					best_i = tmp[0]
					best_v = tmp[1]
					print i,l,best_i,best_v

			assignments[i] = lookup[best_i]
			for l in xrange(3): del refnormal[3*best_i+2-l]
			del lookup[best_i]


	else:
		ERROR("  ERROR:  symmetry not supported  "+sym,"nearest_full_k_projangles",1)
		assignments = []

	return assignments


def nearestk_to_refdir(refnormal, refdir, howmany = 1):
	lookup = range(len(refnormal))
	assignments = [-1]*howmany
	for i in xrange(howmany):
		k = Util.nearest_ang(refnormal, refdir[0],refdir[1],refdir[2])
		assignments[i] = lookup[k]
		del refnormal[3*k+2], refnormal[3*k+1], refnormal[3*k+0], lookup[k]
	return assignments


def nearestk_to_refdirs(refnormal, refdir, howmany = 1):
	lookup = range(len(refnormal))
	assignments = []
	for j in xrange(len(refdir)):
		assignment = [-1]*howmany
		for i in xrange(howmany):
			k = Util.nearest_ang(refnormal, refdir[j][0],refdir[j][1],refdir[j][2])
			assignment[i] = lookup[k]
			del refnormal[3*k+2], refnormal[3*k+1], refnormal[3*k+0], lookup[k]
		assignments.append(assignment)
	return assignments



'''
def assign_projangles(projangles, refangles, return_asg = False):

	if len(refangles) > 10000:
		if len(refangles) > 100000:
			coarse_refangles = even_angles(1.5)   # 9453 angles 
		else:
			coarse_refangles = even_angles(5.0)   # 849 angles
		coarse_asg = assign_projangles(projangles, coarse_refangles, True)
		ref_asg = assign_projangles(refangles, coarse_refangles, True)
	else:
		coarse_refangles = []
		coarse_asg = []
		ref_asg = []

	nproj = len(projangles)
	nref = len(refangles)
	proj_ang = [0.0]*(nproj*2)
	ref_ang = [0.0]*(nref*2)
	for i in xrange(nproj):
		proj_ang[i*2] = projangles[i][0]
		proj_ang[i*2+1] = projangles[i][1]
	for i in xrange(nref):
		ref_ang[i*2] = refangles[i][0]
		ref_ang[i*2+1] = refangles[i][1]

	asg = Util.assign_projangles(proj_ang, ref_ang, coarse_asg, ref_asg, len(coarse_refangles))
	if return_asg: return asg
	assignments = [[] for i in xrange(nref)]
	for i in xrange(nproj):
		assignments[asg[i]].append(i)

	return assignments
'''

def assign_projangles(projangles, refangles, return_asg = False):

	nproj = len(projangles)
	nref = len(refangles)
	proj_ang = [0.0]*(nproj*2)
	ref_ang = [0.0]*(nref*2)
	for i in xrange(nproj):
		proj_ang[i*2] = projangles[i][0]
		proj_ang[i*2+1] = projangles[i][1]
	for i in xrange(nref):
		ref_ang[i*2] = refangles[i][0]
		ref_ang[i*2+1] = refangles[i][1]

	asg = Util.assign_projangles(proj_ang, ref_ang)
	if return_asg: return asg
	assignments = [[] for i in xrange(nref)]
	for i in xrange(nproj):
		assignments[asg[i]].append(i)

	return assignments

def assign_projangles_f(projangles, refangles, return_asg = False):

	asg = Util.assign_projangles_f(projangles, refangles)
	if return_asg: return asg
	assignments = [[] for i in xrange(len(refangles))]
	for i in xrange(len(projangles)):
		assignments[asg[i]].append(i)

	return assignments


def cone_ang( projangles, phi, tht, ant ):
	from utilities import getvec
	from math import cos, pi, degrees, radians
	vec = getvec( phi, tht )

	cone = cos(radians(ant))
	la = []
	for i in xrange( len(projangles) ):
		vecs = getvec( projangles[i][0], projangles[i][1] )
		s = abs(vecs[0]*vec[0] + vecs[1]*vec[1] + vecs[2]*vec[2])
		if s >= cone:
			la.append(projangles[i])

	return la

def cone_ang_f( projangles, phi, tht, ant ):
	from utilities import getvec
	from math import cos, pi, degrees, radians
	# vec = getfvec( phi, tht )
	vec = getfvec( phi, tht )

	cone = cos(radians(ant))
	la = []
	for i in xrange( len(projangles) ):
		# vecs = getfvec( projangles[i][0], projangles[i][1] )
		vecs = getfvec( projangles[i][0], projangles[i][1] )
		s = vecs[0]*vec[0] + vecs[1]*vec[1] + vecs[2]*vec[2]
		if s >= cone:
			la.append(projangles[i])
	return la

def cone_ang_f_with_index( projangles, phi, tht, ant ):
	from utilities import getvec
	from math import cos, pi, degrees, radians
	# vec = getvec( phi, tht )
	vec = getfvec( phi, tht )

	cone = cos(radians(ant))
	la = []
	index = []
	for i in xrange( len(projangles) ):
		# vecs = getvec( projangles[i][0], projangles[i][1] )
		vecs = getfvec( projangles[i][0], projangles[i][1] )
		s = vecs[0]*vec[0] + vecs[1]*vec[1] + vecs[2]*vec[2]
		if s >= cone:
		# if abs(s) >= cone:
			la.append(projangles[i])
			index.append(i)
	return la, index

def cone_ang_with_index( projangles, phi, tht, ant ):
	from utilities import getvec
	from math import cos, pi, degrees, radians
	# vec = getvec( phi, tht )
	vec = getfvec( phi, tht )

	cone = cos(radians(ant))
	la = []
	index = []
	for i in xrange( len(projangles) ):
		# vecs = getvec( projangles[i][0], projangles[i][1] )
		vecs = getfvec( projangles[i][0], projangles[i][1] )
		s = vecs[0]*vec[0] + vecs[1]*vec[1] + vecs[2]*vec[2]
		# if s >= cone:
		if abs(s) >= cone:
			la.append(projangles[i] + [i])
			index.append(i)
	
	return la, index

def cone_vectors( normvectors, phi, tht, ant ):
	from utilities import getvec
	from math import cos, pi, degrees, radians
	vec = getvec( phi, tht )

	cone = cos(radians(ant))
	la = []
	for i in xrange( len(normvectors) ):
		s = abs(normvectors[i][0]*vec[0] + normvectors[i][1]*vec[1] + normvectors[i][2]*vec[2])
		if s >= cone:
			la.append(normvectors[i])

	return la

def disable_bdb_cache():
	import EMAN2db
	EMAN2db.BDB_CACHE_DISABLE = True

def enable_bdb_cache():
	import EMAN2db
	EMAN2db.BDB_CACHE_DISABLE = False

def rotation_between_anglesets(agls1, agls2):
	"""
	  Find an overall 3D rotation (phi theta psi) between two sets of Eulerian angles.
	  The two sets have to have the same number of elements and it is assumed that k'th element on the first
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

	deg2rad = pi/180.0

	def ori2xyz(ori):
		if(type(ori) == types.ListType):
			phi, theta, psi = ori[:3]
		else:
			# it has to be Transformation object
			d = ori.get_params("spider")
			phi   = d["phi"]
			theta = d["theta"]
			psi   = d["psi"]
		"""
		#  This makes no sense here!  PAP 09/2011
		if theta > 90.0:
			phi += 180.0
			theta = 180.0-theta
		"""
		phi   *= deg2rad
		theta *= deg2rad
		x = sin(theta) * sin(phi)
		y = sin(theta) * cos(phi)
		z = cos(theta)

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


def angle_between_projections_directions(proj1, proj2):
	"""
	  It returns angle between two projections directions.
	  INPUT: two lists: [phi1, theta1] , [phi2, theta2]
	  OUTPUT: angle (in degrees)
	"""
	from math import sin, cos, acos, radians, degrees
	phi1   = radians(proj1[0])
	phi2   = radians(proj2[0])
	theta1 = radians(proj1[1])
	theta2 = radians(proj2[1])
	st1 = sin(theta1)
	st2 = sin(theta2)
	ct1 = cos(theta1)
	ct2 = cos(theta2)
	cp1cp2_sp1sp2 = cos(phi1 - phi2)
	temp = st1 * st2 * cp1cp2_sp1sp2 + ct1 * ct2
	if temp < -1.0:
		temp = -1.0
	if temp > 1.0:
		temp = 1.0
	return degrees( acos( temp ) )


def angles_between_anglesets(angleset1, angleset2, indexes=None):
	"""
	  It returns list of angles describing differences between the given anglesets (rotations of anglesets don't matter).
	  The function works as follow:
	  1. use the rotation_between_anglesets function to find the transformation between the anglesets
	  2. apply the transformation to the first angleset (to fit it into the second angleset)
	  3. calculate angles between corresponding projections directions from the anglesets 
	  INPUT: each of the parameters should be a list of pairs [phi, theta] (additional parameters (pdi, shifts) don't influence on the results)
	  OUTPUT: list of floats - angles in degrees (the n-th element of the list equals the angle between n-th projections directions from the anglesets)
	  The third parameter (indexes) is optional and may be set to list of indexes. In that case only elements from given list are taken into account.
	"""
	from EMAN2 import Transform

	if indexes != None:
		new_ang1 = []
		new_ang2 = []
		for i in indexes:
			new_ang1.append(angleset1[i])
			new_ang2.append(angleset2[i])
		angleset1 = new_ang1
		angleset2 = new_ang2

	rot = rotation_between_anglesets(angleset1, angleset2)
	T2 = Transform({"type":"spider","phi":rot[0],"theta":rot[1],"psi":rot[2],"tx":0.0,"ty":0.0,"tz":0.0,"mirror":0,"scale":1.0})
	angles = []
	angle_errors = []
	for i in xrange(len(angleset1)):
		T1 = Transform({"type":"spider","phi":angleset1[i][0],"theta":angleset1[i][1],"psi":angleset1[i][2],"tx":0.0,"ty":0.0,"tz":0.0,"mirror":0,"scale":1.0})
		T  = T1*T2
		d  = T.get_params("spider")
		angles.append([d["phi"], d["theta"], d["psi"]])
		angle_errors.append( angle_between_projections_directions(angles[i], angleset2[i]) )
	return angle_errors


def phi_theta_to_xyz(ang):
	from math import sin, cos, pi, radians
	phi   = radians( ang[0] )
	theta = radians( ang[1] )
	z = cos(theta)
	x = sin(theta) * cos(phi)
	y = sin(theta) * sin(phi)
	return [x, y, z]


def xyz_to_phi_theta(xyz):
	from math import pi, acos, sqrt, degrees, atan2
	theta = acos(xyz[2])
	phi   = atan2(xyz[1], xyz[0])
	return [ degrees(phi), degrees(theta), 0.0]
	
# input: list of triplets (phi, theta, psi)
# output: average triplet: (phi, theta, psi)
def average_angles(angles):
	from math import sqrt
	# convert to x, y, z
	ex = 0.0
	ey = 0.0
	ez = 0.0
	sum_psi = 0.0
	for ang in angles:
		xyz = phi_theta_to_xyz(ang)
		ex += xyz[0]
		ey += xyz[1]
		ez += xyz[2]
		sum_psi += ang[2]
	ex /= len(ang)
	ey /= len(ang)
	ez /= len(ang)
	divider = sqrt(ex*ex + ey*ey + ez*ez)
	xyz = [ ex/divider, ey/divider, ez/divider ]
	r = xyz_to_phi_theta(xyz)
	r[2] = sum_psi / len(angles)  # TODO - correct psi calculations !!!!
	return r


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

def group_proj_by_phitheta_slow(proj_ang, symmetry = "c1", img_per_grp = 100, verbose = False):
	from time import time
	from math import exp, pi

	def get_ref_ang_list(delta, sym):
		ref_ang = even_angles(delta, symmetry=sym)
		ref_ang_list = [0.0]*(len(ref_ang)*2)
		for i in xrange(len(ref_ang)):
			ref_ang_list[2*i] = ref_ang[i][0]
			ref_ang_list[2*i+1] = ref_ang[i][1]
		return ref_ang_list, len(ref_ang)

	def gv(phi, theta):
		from math import pi, cos, sin
		angle_to_rad = pi/180.0

		theta *= angle_to_rad
		phi *= angle_to_rad

		x = sin(theta)*cos(phi) 
		y = sin(theta)*sin(phi)
		z = cos(theta)

		return (x, y, z)

	def ang_diff(v1, v2):
		# The first return value is the angle between two vectors
		# The second return value is whether we need to mirror one of them (0 - no need, 1 - need)
		from math import acos, pi, degrees

		v = v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2]
		v = max(min(v,1.0),-1.0)
		if v >= 0: return degrees(acos(v)), 0
		else:      return degrees(acos(-v)), 1

	t0 = time()
	proj_list = []   
	angles_list = []
	N = len(proj_ang)
	if len(proj_ang[0]) == 3:       # determine whether it has shifts provided, make the program more robust
		for i in xrange(N):
			proj_ang[i].append(i)
			proj_ang[i].append(True)
			vec = gv(proj_ang[i][0], proj_ang[i][1])     # pre-calculate the vector for each projection angles
			proj_ang[i].append(vec)
	else:
		for i in xrange(N):
			proj_ang[i][3] = i
			proj_ang[i][4] = True
			vec = gv(proj_ang[i][0], proj_ang[i][1])     # pre-calculate the vector for each projection angles
			proj_ang[i].append(vec)

	ref_ang_list1, nref1 = get_ref_ang_list(20.0, sym = symmetry)   
	ref_ang_list2, nref2 = get_ref_ang_list(10.0, sym = symmetry)   
	ref_ang_list3, nref3 = get_ref_ang_list(5.0, sym = symmetry)    
	ref_ang_list4, nref4 = get_ref_ang_list(2.5, sym = symmetry)

	c = 100
	L = max(100, img_per_grp)
	# This is to record whether we are considering the same group as before
	# If we are, we are only going to read the table and avoid calculating the distance again.
	previous_group = -1
	previous_zone = 5
	for grp in xrange(N/img_per_grp):
		print grp,
		N_remain = N-grp*img_per_grp
		# The idea here is that if each group has more than 100 images in average, 
		# we consider it crowded enough to just consider the most crowded group.
		if N_remain >= nref4*L:
			ref_ang_list = ref_ang_list4
			nref = nref4
			if previous_zone > 4:
				previous_group = -1
				previous_zone = 4
		elif N_remain >= nref3*L:
			ref_ang_list = ref_ang_list3
			nref = nref3
			if previous_zone > 3:
				previous_group = -1
				previous_zone = 3
		elif N_remain >= nref2*L:
			ref_ang_list = ref_ang_list2
			nref = nref2
			if previous_zone > 2:
				previous_group = -1
				previous_zone = 2
		elif N_remain >= nref1*L:
			ref_ang_list = ref_ang_list1
			nref = nref1
			if previous_zone > 1:
				previous_group = -1
				previous_zone = 1
		else:
			if previous_zone > 0:
				previous_group = -1
				previous_zone = 0

		t1 = time()
		v = []
		index = []
		if N_remain >= nref1*L:
			# In this case, assign all projection to groups and only consider the most crowded group.
			proj_ang_list = [0.0]*(N_remain*2)
			nn = 0
			remain_index = [0]*N_remain
			for i in xrange(N):
				if proj_ang[i][4]:
					proj_ang_list[nn*2] = proj_ang[i][0]
					proj_ang_list[nn*2+1] = proj_ang[i][1]
					remain_index[nn] = i
					nn += 1
			asg = Util.assign_projangles(proj_ang_list, ref_ang_list)
			assignments = [[] for i in xrange(nref)]
			for i in xrange(N_remain):
				assignments[asg[i]].append(i)
			# find the largest group and record the group size and group number 
			max_group_size = 0
			max_group = -1
			for i in xrange(nref):
				if len(assignments[i]) > max_group_size:
					max_group_size = len(assignments[i])
					max_group = i
			print max_group_size, max_group, previous_group,
			for i in xrange(len(assignments[max_group])):
				ind = remain_index[assignments[max_group][i]]
				v.append(proj_ang[ind][5])
				index.append(ind)
		else:
			# In this case, use all the projections available
			for i in xrange(N):
				if proj_ang[i][4]:
					v.append(proj_ang[i][5])
					index.append(i)
			max_group = 0

		t2 = time()
		Nn = len(index)
		density = [[0.0, 0] for i in xrange(Nn)]
		if max_group != previous_group:
			diff_table = [[0.0 for i in xrange(Nn)] for j in xrange(Nn)]
			for i in xrange(Nn-1):
				for j in xrange(i+1, Nn):
					diff = ang_diff(v[i], v[j])
					q = exp(-c*(diff[0]/180.0*pi)**2)
					diff_table[i][j] = q
					diff_table[j][i] = q
			diff_table_index = dict()
			for i in xrange(Nn): diff_table_index[index[i]] = i
			print Nn, True, 
		else:
			print Nn, False,

		t21 = time()
		for i in xrange(Nn):
			density[i][0] = sum(diff_table[diff_table_index[index[i]]])
			density[i][1] = i
		t22 = time()
		density.sort(reverse=True)

		t3 = time()
		dang = [[0.0, 0] for i in xrange(Nn)]
		most_dense_point = density[0][1]
		for i in xrange(Nn):
			diff = ang_diff(v[i], v[most_dense_point])
			dang[i][0] = diff[0]
			dang[i][1] = i
		dang[most_dense_point][0] = -1.
		dang.sort()

		t4 = time()
		members = [0]*img_per_grp
		for i in xrange(img_per_grp):
			idd = index[dang[i][1]]
			for j in xrange(len(diff_table)):
				diff_table[diff_table_index[idd]][j] = 0.0
				diff_table[j][diff_table_index[idd]] = 0.0
			members[i] = idd
			proj_ang[members[i]][4] = False
		proj_list.append(members)
		center_i = index[dang[0][1]]
		angles_list.append([proj_ang[center_i][0], proj_ang[center_i][1], dang[img_per_grp-1][0]])
		previous_group = max_group
		print t2-t1, t3-t2, t22-t21, t3-t22, t4-t3

	if N%img_per_grp*3 >= 2*img_per_grp:
		members = []
		for i in xrange(N):
			if proj_ang[i][4]:
				members.append(i)
		proj_list.append(members)
		angles_list.append([proj_ang[members[0]][0], proj_ang[members[0]][1], 90.0])
	elif N%img_per_grp != 0:
		for i in xrange(N):
			if proj_ang[i][4]:
				proj_list[-1].append(i)
	print "Total time used = ", time()-t0

	return proj_list, angles_list

def group_proj_by_phitheta(proj_ang, symmetry = "c1", img_per_grp = 100, verbose = False):
	from math import exp, pi

	def gv(phi, theta):
		from math import pi, cos, sin
		angle_to_rad = pi/180.0

		theta *= angle_to_rad
		phi *= angle_to_rad

		x = sin(theta)*cos(phi) 
		y = sin(theta)*sin(phi)
		z = cos(theta)

		return (x, y, z)

	def ang_diff(v1, v2):
		# The first return value is the angle between two vectors
		# The second return value is whether we need to mirror one of them (0 - no need, 1 - need)
		from math import acos, pi

		v = v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2]
		if v > 1: v = 1
		if v < -1: v = -1
		if v >= 0: return acos(v)*180/pi, 0
		else:  return acos(-v)*180/pi, 1

	
	def get_ref_ang_list(delta, sym):
		ref_ang = even_angles(delta, symmetry=sym)
		ref_ang_list = [0.0]*(len(ref_ang)*2)
		for i in xrange(len(ref_ang)):
			ref_ang_list[2*i] = ref_ang[i][0]
			ref_ang_list[2*i+1] = ref_ang[i][1]
		return ref_ang_list, len(ref_ang)
	
	N = len(proj_ang)
	proj_ang_list = [0]*(N*2)
	for i in xrange(N):
		proj_ang_list[i*2] = proj_ang[i][0]
		proj_ang_list[i*2+1] = proj_ang[i][1]
	
	ref_ang_list1, nref1 = get_ref_ang_list(20.0, sym = symmetry)   
	ref_ang_list2, nref2 = get_ref_ang_list(10.0, sym = symmetry)   
	ref_ang_list3, nref3 = get_ref_ang_list(5.0, sym = symmetry)    
	ref_ang_list4, nref4 = get_ref_ang_list(2.5, sym = symmetry)

	ref_ang_list = []
	ref_ang_list.extend(ref_ang_list1)
	ref_ang_list.extend(ref_ang_list2)
	ref_ang_list.extend(ref_ang_list3)
	ref_ang_list.extend(ref_ang_list4)
	ref_ang_list.append(nref1)
	ref_ang_list.append(nref2)
	ref_ang_list.append(nref3)
	ref_ang_list.append(nref4)
	proj_list = Util.group_proj_by_phitheta(proj_ang_list, ref_ang_list, img_per_grp)
	
	proj_list2 = proj_list[:]
	for i in xrange(len(proj_list2)): proj_list2[i] = abs(proj_list2[i])
	proj_list2.sort()
	assert N == len(proj_list2)
	for i in xrange(N): assert i == proj_list2[i]
	
	Ng = N/img_per_grp
	proj_list_new = [[] for i in xrange(Ng)]
	mirror_list = [[] for i in xrange(Ng)]
	angles_list = []
	for i in xrange(Ng):
		for j in xrange(img_per_grp):
			proj_list_new[i].append(abs(proj_list[i*img_per_grp+j]));
			mirror_list[i].append(proj_list[i*img_per_grp+j] >= 0)
		phi1 = proj_ang[proj_list_new[i][0]][0];
		theta1 = proj_ang[proj_list_new[i][0]][1];
		phi2 = proj_ang[proj_list_new[i][-1]][0];
		theta2 = proj_ang[proj_list_new[i][-1]][1];
		angles_list.append([phi1, theta1, ang_diff(gv(phi1, theta1), gv(phi2, theta2))[0]]);

	if N%img_per_grp*3 >= 2*img_per_grp:
		proj_list_new.append([])
		mirror_list.append([])
		for i in xrange(Ng*img_per_grp, N):
			proj_list_new[-1].append(abs(proj_list[i]));
			mirror_list[-1].append(proj_list[i] >= 0)
		phi1 = proj_ang[proj_list_new[Ng][0]][0];
		theta1 = proj_ang[proj_list_new[Ng][0]][1];
		phi2 = proj_ang[proj_list_new[Ng][-1]][0];
		theta2 = proj_ang[proj_list_new[Ng][-1]][1];
		angles_list.append([phi1, theta1, ang_diff(gv(phi1, theta1), gv(phi2, theta2))[0]]);
	elif N%img_per_grp != 0:
		for i in xrange(Ng*img_per_grp, N):
			proj_list_new[-1].append(abs(proj_list[i]))
			mirror_list[-1].append(proj_list[i] >= 0)

	return proj_list_new, angles_list, mirror_list

def nearest_proj(proj_ang, img_per_grp=100, List=[]):
	from math import exp, pi
	from sets import Set
	from time import time
	from random import randint

	def gv(phi, theta):
		from math import radians, cos, sin

		theta = radians(theta)
		phi   = radians(phi)

		x = sin(theta)*cos(phi) 
		y = sin(theta)*sin(phi)
		z = cos(theta)

		return (x, y, z)

	def ang_diff(v1, v2):
		# The first return value is the angle between two vectors
		# The second return value is whether we need to mirror one of them (0 - no need, 1 - need)
		from math import acos, degrees

		v = v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2]
		if v > 1: v = 1
		if v < -1: v = -1
		if v >= 0: return degrees(acos(v)), 0
		else:  return degrees(acos(-v)), 1
	
	def get_ref_ang_list(delta, sym):
		ref_ang = even_angles(delta, symmetry=sym)
		ref_ang_list = [0.0]*(len(ref_ang)*2)
		for i in xrange(len(ref_ang)):
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
	if len(List) == 0: List = range(N)
	if N < img_per_grp:
		print "Error: image per group larger than the number of particles!"
		exit()
	phi_list   = [[0.0, 0] for i in xrange(N)]
	theta_list = [[0.0, 0] for i in xrange(N)]
	vec = [None]*N
	for i in xrange(N):
		phi = proj_ang[i][0]
		theta = proj_ang[i][1]
		vec[i] = gv(phi, theta)
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
	for i in xrange(N):
		theta_list_l[i] = theta_list[i][0]
		phi_list_l[i] = phi_list[i][0]
	
	g = [[360.0, 0, 0] for i in xrange(N)]
	proj_list = []
	mirror_list = []
	neighbor   = [0]*img_per_grp
	#neighbor2 = [0]*img_per_grp
	dis        = [0.0]*img_per_grp
	#dis2      = [0.0]*img_per_grp
	mirror     = [0]*img_per_grp
	S = Set()
	T = Set()
	#tt1 = time()
	for i in xrange(len(List)):
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
				for j in xrange(phi_left_bound, phi_right_bound+1):
					S.add(phi_list[j][1])
			else:
				for j in xrange(phi_right_bound+1):
					S.add(phi_list[j][1])
				for j in xrange(phi_left_bound, N):
					S.add(phi_list[j][1])
			if theta+delta > 90.0:
				if phi_mir_left_bound < phi_mir_right_bound:
					for j in xrange(phi_mir_left_bound, phi_mir_right_bound+1):
						S.add(phi_list[j][1])
				else:
					for j in xrange(phi_mir_right_bound+1):
						S.add(phi_list[j][1])
					for j in xrange(phi_mir_left_bound, N):
						S.add(phi_list[j][1])
			for j in xrange(theta_left_bound, theta_right_bound+1):
				T.add(theta_list[j][1])
			v = list(T.intersection(S))
			S.clear()
			T.clear()
			if len(v) >= min(1.5*img_per_grp, N): break
			delta *= 2
			del v

		for j in xrange(len(v)):
			d = ang_diff(vec[v[j]], vec[k])
			g[j][0] = d[0]
			if v[j] == k: g[j][0] = -1.  # To ensure the image itself is always included in the group
			g[j][1] = d[1]
			g[j][2] = v[j]
		g[:len(v)] = sorted(g[:len(v)])
		for j in xrange(img_per_grp):
			neighbor[j] = g[j][2]
			dis[j] = g[j][0]
			mirror[j] = (g[j][1] == 1)
		proj_list.append(neighbor[:])
		mirror_list.append(mirror[:])
		#t2 = time()

		'''
		for j in xrange(N):
			d = ang_diff(vec[j], vec[k])
			g[j][0] = d[0]
			g[j][1] = d[1]
			g[j][2] = j
		g.sort()
		for j in xrange(img_per_grp):
			neighbor2[j] = g[j][2]
			dis2[j] = g[j][0]
		t3 = time()
		print "Members in common = %3d   extra delta = %6.3f   time1 = %5.2f   time2 = %5.2f"%(len(Set(neighbor).intersection(Set(neighbor2))),
		dis[-1]-dis2[-1], t2-t1, t3-t2)
		'''
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
"""
def findall(val, lo):
	'''
	  Find all occurrences of val in list lo
	  Returns a list of indices of val in lo.
	'''
	u = []
	i = -1
	while( i < len(lo)-1):
		try:
			i = lo.index(val,i+1)
			u.append(i)
		except:
			i += 1
	return  u
"""

# parameters: list of integers, number of processors
def chunks_distribution(chunks, procs):
	from heapq import heappush, heappop

	# sort chunks in descending order
	chunks.sort(reverse=True)

	# create heap and list with solution    
	results = []
	heap = []
	for p in xrange(procs):
		results.append([])
		heappush(heap, (0, p))

	# main calculations
	# following chunks are added to the least loaded processors
	for c in chunks:
		s, p = heappop(heap)
		results[p].append(c)
		s += c[0]
		heappush(heap, (s, p))

	return results

"""
Iterators over sequence of images. They work for lists and stacks of images.
Additional time cost: about 1 second per 10^6 iterations on 3 GHz processor.

Usage:

it = iterImagesList(list_of_images)  <-or->  it = iterImagesStack(stack_with_images)
while it.goToNext():
	do_something(it.image())
	
"""
# ================ Iterator for list of images
class iterImagesList:
	images = []
	imagesIndexes = []
	position = -1
	def __init__(self, list_of_images, list_of_indexes = None):
		if list_of_indexes == None:
			self.images = list_of_images[:]
			self.imagesIndexes = range(len(self.images))
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
class iterImagesStack:
	stackName = ""
	currentImage = None
	imagesIndexes = []
	position = -1
	def __init__(self, stack_name, list_of_indexes = None):
		if list_of_indexes == None:
			self.imagesIndexes = range(EMUtil.get_image_count(stack_name))
		else:
			self.imagesIndexes = list_of_indexes[:]
		self.stackName = stack_name
	def iterNo(self):
		return self.position
	def imageIndex(self):
		return self.imagesIndexes[self.position]
	def image(self):
		if self.currentImage == None:
			self.currentImage = EMData()
			self.currentImage.read_image(self.stackName, self.imagesIndexes[self.position])
		return self.currentImage
	def goToNext(self):
		self.currentImage = None
		if len(self.imagesIndexes) <= self.position:
			return False
		self.position += 1
		return (self.position < len(self.imagesIndexes))
	def goToPrev(self):
		self.currentImage = None
		if 0 > self.position:
			return False
		self.position -= 1
		return (self.position >= 0)


from cPickle import dumps,loads
from zlib import compress,decompress
from struct import pack,unpack

def pack_message(data):
	"""Convert data for transmission efficiently"""

	if isinstance(data,str):
		if len(data)>256 : return "C"+compress(data,1)
		else : return "S"+data
	else :
		d2x=dumps(data,-1)
		if len(d2x)>256 : return "Z"+compress(d2x,1)
		else : return "O"+d2x

	
def unpack_message(msg):
	"""Unpack a data payload prepared by pack_message"""
	
	if msg[0]=="C" : return decompress((msg[1:]).tostring())
	elif msg[0]=="S" : return (msg[1:]).tostring()
	elif msg[0]=="Z" : return loads(decompress((msg[1:]).tostring()))
	elif msg[0]=="O" : return loads((msg[1:]).tostring())
	else :
		print "ERROR: Invalid MPI message. Please contact developers. (%s)"%str(msg[:20])
		raise Exception("unpack_message")


statistics_send_recv = dict()

def update_tag(communicator, target_rank):   # TODO - it doesn't work when communicators are destroyed and recreated
	return 123456
	global statistics_send_recv
	if not statistics_send_recv.has_key(communicator):
		from mpi import mpi_comm_size
		statistics_send_recv[communicator] = [0] * mpi_comm_size(communicator)
	statistics_send_recv[communicator][target_rank] += 1
	return statistics_send_recv[communicator][target_rank]

# ===================================== WRAPPER FOR MPI

def wrap_mpi_send(data, destination, communicator = None):
	from mpi import mpi_send, MPI_COMM_WORLD, MPI_CHAR

	if communicator == None:
		communicator = MPI_COMM_WORLD
		
	msg = pack_message(data)
	tag = update_tag(communicator, destination)
	#from mpi import mpi_comm_rank
	#print communicator, mpi_comm_rank(communicator), "send to", destination, tag
	mpi_send(msg, len(msg), MPI_CHAR, destination, tag, communicator) # int MPI_Send( void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm )


def wrap_mpi_recv(source, communicator = None):
	from mpi import mpi_recv, MPI_COMM_WORLD, MPI_CHAR, mpi_probe, mpi_get_count

	if communicator == None:
		communicator = MPI_COMM_WORLD
	
	tag = update_tag(communicator, source)
	#from mpi import mpi_comm_rank
	#print communicator, mpi_comm_rank(communicator), "recv from", source, tag
	mpi_probe(source, tag, communicator)
	n = mpi_get_count(MPI_CHAR)
	msg = mpi_recv(n, MPI_CHAR, source, tag, communicator)
	return unpack_message(msg)


def wrap_mpi_bcast(data, root, communicator = None):
	from mpi import mpi_bcast, MPI_COMM_WORLD, mpi_comm_rank, MPI_CHAR
	
	if communicator == None:
		communicator = MPI_COMM_WORLD

	rank = mpi_comm_rank(communicator)

	if rank == root:
		msg = pack_message(data)
		n = pack("I",len(msg))
	else:
		msg = None
		n = None

	n = mpi_bcast(n, 4, MPI_CHAR, root, communicator)  # int MPI_Bcast ( void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm ) 
	n=unpack("I",n)[0]
	msg = mpi_bcast(msg, n, MPI_CHAR, root, communicator)  # int MPI_Bcast ( void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm ) 
	return unpack_message(msg)


# data must be a python list (numpy array also should be implemented)
def wrap_mpi_gatherv(data, root, communicator = None):
	from mpi import mpi_comm_rank, mpi_comm_size, MPI_COMM_WORLD

	if communicator == None:
		communicator = MPI_COMM_WORLD

	rank = mpi_comm_rank(communicator)
	procs = mpi_comm_size(communicator)

	out_array = None
	if rank == root:
		if type(data) is list:
			out_array = []
			for p in xrange(procs):
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

def wrap_mpi_split(comm, no_of_groups):
	"""

	Takes the processes of a communicator (comm) and splits them in groups (no_of_groups).
	Each subgroup of processes has ids generated from 0 to number of processes per group - 1.
	Consecutive global process ids have consecutive subgroup process ids.

	"""
	from mpi import mpi_comm_size, mpi_comm_rank, mpi_comm_split
	nproc = mpi_comm_size(comm)
	myid = mpi_comm_rank(comm)

	no_of_proc_per_group = nproc / no_of_groups
	color = myid / no_of_proc_per_group
	key = myid % no_of_proc_per_group

	return mpi_comm_split(comm, color, key)

def get_dist(c1, c2):
	from math import sqrt
	d = sqrt((c1[0] - c2[0])**2 + (c1[1] - c2[1])**2)
	return d


def eliminate_moons(my_volume, moon_elimination_params):
	"""
	moon_elimination_params[0] - mass in KDa
	moon_elimination_params[1] - resolution in px/A
	"""

	from morphology import binarize
	histogram_threshold  =  my_volume.find_3d_threshold(moon_elimination_params[0], moon_elimination_params[1])*1.1
	# clean11 88x88,  4.84 px/A 750 kDa

	my_volume_binarized = binarize(my_volume, histogram_threshold)
	# my_volume_binarized.write_image ("my_volume_binarized.hdf")
	my_volume_binarized_with_no_moons = Util.get_biggest_cluster(my_volume_binarized)
	# my_volume_binarized_with_no_moons.write_image("my_volume_binarized_with_no_moons.hdf")
	volume_difference = my_volume_binarized - my_volume_binarized_with_no_moons
	# volume_difference.write_image("volume_difference.hdf")

	if volume_difference.get_value_at(volume_difference.calc_max_index()) == 0 and \
		volume_difference.get_value_at(volume_difference.calc_min_index()) == 0:
		return my_volume
	else:
		from utilities import gauss_edge
		return gauss_edge(my_volume_binarized_with_no_moons) * my_volume

		# from utilities   import model_blank
		# # mask = model_blank(my_volume_binarized_with_no_moons.get_xsize(), my_volume_binarized_with_no_moons.get_ysize(), my_volume_binarized_with_no_moons.get_zsize())
		# # mask.to_one()
	# this is only in master
	
def combinations_of_n_taken_by_k(n, k):
	from fractions import Fraction
	return int(reduce(lambda x, y: x * y, (Fraction(n-i, i+1) for i in range(k)), 1))
	
def cmdexecute(cmd):
	from   time import localtime, strftime
	import subprocess
	outcome = subprocess.call(cmd, shell=True)
	line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
	if(outcome == 1):
		print  line,"ERROR!!   Command failed:  ", cmd
		from sys import exit
		exit()
	else:  print line,"Executed successfully: ",cmd

def string_found_in_file(myregex, filename):
	import re
	pattern = re.compile(myregex)
	for line in open(filename):
		if re.findall(pattern, line) <> []:
			return True
	return False

def random_string(length_of_randomstring = 16):
	import random
	chars=map(chr, range(97, 123)) # a..z
	chars.extend(map(chr, range(65, 91))) # A..Z
	chars.extend(map(chr, range(48, 58))) # 0..9
	random_string = ""
	for i in xrange(length_of_randomstring):
		random_string += chars[random.randint(0,len(chars)-1)]
	return random_string

def get_latest_directory_increment_value(directory_location, directory_name, start_value = 1, myformat = "%03d"):
	import os
	dir_count = start_value
	while os.path.isdir(directory_location + directory_name + myformat%(dir_count)):
		dir_count += 1
	if dir_count == start_value:
		return start_value
	return dir_count - 1

def get_nonexistent_directory_increment_value(directory_location, directory_name, start_value = 1, myformat = "%03d"):
	import os
	dir_count = start_value
	while os.path.isdir(directory_location + directory_name + myformat%(dir_count)):
		dir_count += 1
	return dir_count

def print_with_time_info(msg):
	from time import localtime, strftime
	line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>" + msg
	print line

def if_error_then_all_processes_exit_program(error_status):
	import sys, os
	from utilities import print_msg
	
	if "OMPI_COMM_WORLD_SIZE" not in os.environ:
		def mpi_comm_rank(n): return 0
		def mpi_bcast(*largs):
			return [largs[0]]
		def mpi_finalize():
			return None
		MPI_INT, MPI_COMM_WORLD = 0, 0
	else:
		from mpi import mpi_comm_rank, mpi_bcast, mpi_finalize, MPI_INT, MPI_COMM_WORLD
	
	myid = mpi_comm_rank(MPI_COMM_WORLD)
	if error_status != None and error_status != 0:
		error_status_info = error_status
		error_status = 1
	else:
		error_status = 0

	error_status = mpi_bcast(error_status, 1, MPI_INT, 0, MPI_COMM_WORLD)
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
		mpi_finalize()
		sys.exit(1)

def get_shrink_data_huang(Tracker, nxinit, partids, partstack, myid, main_node, nproc, preshift = False):
	# The function will read from stack a subset of images specified in partids
	#   and assign to them parameters from partstack with optional CTF application and shifting of the data.
	# So, the lengths of partids and partstack are the same.
	#  The read data is properly distributed among MPI threads.
	# 10142015 --- preshift is set to True when doing 3-D sorting. 
	# chunk_id are set when data is read in	

	from fundamentals import resample, fshift
	from filter import filt_ctf
	from applications import MPI_start_end

	'''
	if( myid == main_node ):
		print "  "
		line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
		print  line, "Reading data  onx: %3d, nx: %3d, CTF: %s, applyctf: %s, preshift: %s."%(Tracker["constants"]["nnxo"], nxinit, Tracker["constants"]["CTF"], Tracker["applyctf"], preshift)
		print  "                       stack:      %s\n                       partids:     %s\n                       partstack: %s\n"%(Tracker["constants"]["stack"], partids, partstack)
	'''
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
		image_start, image_end = MPI_start_end(ndata, nproc, myid)
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
	for im in xrange(nima):
		data[im] = get_im(Tracker["constants"]["stack"], lpartids[im])
		if im ==0:
			if data[im].get_xsize() > Tracker["constants"]["nnxo"]:
				window_particle =True
				from EMAN2 import Region
			else:
				window_particle =False
		phi,theta,psi,sx,sy = partstack[lpartids[im]][0], partstack[lpartids[im]][1], partstack[lpartids[im]][2], partstack[lpartids[im]][3], partstack[lpartids[im]][4]
		if( Tracker["constants"]["CTF"] and Tracker["applyctf"] ):
			ctf_params = data[im].get_attr("ctf")
			st = Util.infomask(data[im], mask2D, False)
			data[im] -= st[0]
			data[im] = filt_ctf(data[im], ctf_params)
			data[im].set_attr('ctf_applied', 1)
		if preshift:# always true
			data[im] = fshift(data[im], sx, sy)
			set_params_proj(data[im],[phi,theta,psi,0.0,0.0])
			sx = 0.0
			sy = 0.0
		if window_particle:
			mx = data[im].get_xsize()//2-Tracker["constants"]["nnxo"]//2
			my = data[im].get_ysize()//2-Tracker["constants"]["nnxo"]//2
			data[im] = data[im].get_clip(Region(mx,my,Tracker["constants"]["nnxo"],Tracker["constants"]["nnxo"]))
			data[im].set_attr('ctf_applied', 1)
			set_params_proj(data[im],[phi,theta,psi,0.0,0.0])
		#oldshifts[im] = [sx,sy]
		#  resample will properly adjusts shifts and pixel size in ctf
		data[im] = resample(data[im], shrinkage)
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

'''
def get_shrink_data(Tracker, nxinit, partids, partstack, bckgdata = None, myid = 0, main_node = 0, nproc = 1, \
					original_data = None, return_real = False, preshift = False, apply_mask = True, large_memory = True):
	"""
	This function will read from stack a subset of images specified in partids
	   and assign to them parameters from partstack with optional CTF application and shifting of the data.
	So, the lengths of partids and partstack are the same.
	  The read data is properly distributed among MPI threads.
	
	Flow of data:
	1. Read images, if there is enough memory, keep them as original_data.
	2. Read current params
	3.  Apply shift
	4.  Normalize outside of the radius
	5.  Do noise substitution and cosine mask.  (Optional?)
	6.  Shrink data.
	7.  Apply CTF.
	
	"""
	#from fundamentals import resample
	from utilities    import get_im, model_gauss_noise, set_params_proj, get_params_proj
	from fundamentals import fdecimate, fshift, fft
	from filter       import filt_ctf, filt_table
	from applications import MPI_start_end
	from math         import sqrt
	

	if( myid == main_node ):
		print "  "
		line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
		if(original_data == None or not large_memory):
			print  line, "Reading data  onx: %3d, nx: %3d, CTF: %s, applymask: %s, applyctf: %s, preshift: %s."%(Tracker["constants"]["nnxo"], nxinit, Tracker["constants"]["CTF"], apply_mask, Tracker["applyctf"], preshift)
		else:
			print  line, "Processing data  onx: %3d, nx: %3d, CTF: %s, applymask: %s, applyctf: %s, preshift: %s."%(Tracker["constants"]["nnxo"], nxinit, Tracker["constants"]["CTF"], apply_mask, Tracker["applyctf"], preshift)
		print  "                       stack:      %s\n                       partids:     %s\n                       partstack: %s\n"%(Tracker["constants"]["stack"], partids, partstack)
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
		image_start, image_end = MPI_start_end(ndata, nproc, myid)
	lpartids  = lpartids[image_start:image_end]
	partstack = partstack[image_start:image_end]
	#  Preprocess the data
	mask2D  = model_circle(Tracker["constants"]["radius"],Tracker["constants"]["nnxo"],Tracker["constants"]["nnxo"])
	nima = image_end - image_start
	oldshifts = [[0.0,0.0]]*nima
	data = [None]*nima
	if(original_data == None or not large_memory): original_data = [None]*nima
	shrinkage = nxinit/float(Tracker["constants"]["nnxo"])


	#  Note these are in Fortran notation for polar searches
	#txm = float(nxinit-(nxinit//2+1) - radius -1)
	#txl = float(2 + radius - nxinit//2+1)
	radius = int(Tracker["constants"]["radius"]*shrinkage + 0.5)
	txm = float(nxinit-(nxinit//2+1) - radius)
	txl = float(radius - nxinit//2+1)

	if bckgdata :
		nnx = bckgdata[0].get_xsize()
		nny = bckgdata[0].get_ysize()
		bckgnoise = []
		oneover = []
		for i in xrange(nny):
			prj = [0.0]*nnx
			prj[0] = 1.0
			for k in xrange(1,nnx): prj[k] = bckgdata[0].get_value_at(k,i)
			oneover.append(prj)
			for k in xrange(1,nnx):
				if( prj[k] > 0.0 ):  prj[k] = 1.0/sqrt(prj[k])
			bckgnoise.append(prj)

		datastamp = bckgdata[1]

	for im in xrange(nima):
		if(original_data[im] == None or not large_memory):
			original_data[im] = get_im(Tracker["constants"]["stack"], lpartids[im])

		phi,theta,psi,sx,sy = partstack[im][0], partstack[im][1], partstack[im][2], partstack[im][3], partstack[im][4]
		if preshift:
			data[im] = fshift(original_data[im], sx, sy)
			oldshifts[im] = [sx,sy]
			sx = 0.0
			sy = 0.0
		else:  data[im] = original_data[im].copy()
		st = Util.infomask(data[im], mask2D, False)
		data[im] -= st[0]
		data[im] /= st[1]
		if data[im].get_attr_default("bckgnoise", None) :  data[im].delete_attr("bckgnoise")
		#  Do bckgnoise if exists
		if bckgdata:
			try:
				stmp = data[im].get_attr("ptcl_source_image")
			except:
				try:
					stmp = data[im].get_attr("ctf")
					stmp = round(stmp.defocus,4)
				except:
					ERROR("Either ptcl_source_image or ctf has to be present in the header.","get_shrink_data",1, myid)
			try:
				indx = datastamp.index(stmp)
			except:
				ERROR("Problem with indexing ptcl_source_image.","get_shrink_data",1, myid)

			data[im].set_attr("bckgnoise",oneover[indx])
			if apply_mask:
				bckg = model_gauss_noise(1.0,Tracker["constants"]["nnxo"]+2,Tracker["constants"]["nnxo"])
				bckg.set_attr("is_complex",1)
				bckg.set_attr("is_fftpad",1)
				bckg = fft(filt_table(bckg,bckgnoise[indx]))
				#  Normalize bckg noise in real space, only region actually used.
				st = Util.infomask(bckg, mask2D, False)
				bckg -= st[0]
				bckg /= st[1]
				from morphology import cosinemask
				data[im] = cosinemask(data[im],radius = Tracker["constants"]["radius"], bckg = bckg)
		else:
			#  if no bckgnoise, do simple masking instead
			if apply_mask:  data[im] = cosinemask(data[im],radius = Tracker["constants"]["radius"] )

		if( Tracker["constants"]["CTF"] and Tracker["applyctf"] ):
			data[im] = filt_ctf(data[im], data[im].get_attr("ctf"))
			data[im].set_attr('ctf_applied', 1)
		else:  apix = Tracker["constants"]["pixel_size"]

		#  resample will properly adjusts shifts and pixel size in ctf
		#data[im] = resample(data[im], shrinkage)
		#  return Fourier image
		data[im] = fdecimate(data[im], nxinit, nxinit, 1, return_real)
		try:
			ctf_params = original_data[im].get_attr("ctf")
			ctf_params.apix = ctf_params.apix/shrinkage
			data[im].set_attr('ctf', ctf_params)
		except:  pass
		if( Tracker["constants"]["CTF"] and Tracker["applyctf"] ):
			pass
		else:  data[im].set_attr('apix', apix/shrinkage)
		#  We have to make sure the shifts are within correct range, shrinkage or not
		set_params_proj(data[im],[phi,theta,psi,max(min(sx*shrinkage,txm),txl),max(min(sy*shrinkage,txm),txl)])
		#  For local SHC set anchor
		#if(nsoft == 1 and an[0] > -1):
		#  We will always set it to simplify the code
		###set_params_proj(data[im],[phi,theta,psi,0.0,0.0], "xform.anchor")
	assert( nxinit == data[0].get_ysize() )  #  Just to make sure.
	#oldshifts = wrap_mpi_gatherv(oldshifts, main_node, MPI_COMM_WORLD)
	return data, oldshifts, original_data
'''
def getindexdata(stack, partids, partstack, myid, nproc):
	# The function will read from stack a subset of images specified in partids
	#   and assign to them parameters from partstack
	# So, the lengths of partids and partstack are the same.
	#  The read data is properly distributed among MPI threads.
	
	from applications import MPI_start_end

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
		image_start, image_end = MPI_start_end(ndata, nproc, myid)
	lpartids  = lpartids[image_start:image_end]
	partstack = partstack[image_start:image_end]
	data = EMData.read_images(stack, lpartids)

	for i in xrange(len(partstack)):
		set_params_proj(data[i], partstack[i])
	return data


def store_value_of_simple_vars_in_json_file(filename, local_vars, exclude_list_of_vars = [], write_or_append = "w", 
	vars_that_will_show_only_size = []):
	
	import json, types, collections
	 
	allowed_types = [types.NoneType, types.BooleanType, types.IntType, types.LongType, types.FloatType, types.ComplexType,
					 types.UnicodeType, types.StringType]
	
	local_vars_keys = local_vars.keys()

	my_vars = dict()
	for key in set(local_vars_keys) - set(exclude_list_of_vars):
		if type(local_vars[key]) in allowed_types:
			my_vars[key] = local_vars[key]
		elif type(local_vars[key]) in [types.ListType, types.TupleType, type(set())]:
			if len({type(i) for i in local_vars[key]} - set(allowed_types)) == 0:
				if key in vars_that_will_show_only_size:
					my_vars[key] = "%s with length: %d"%(str(type(local_vars[key])),len(local_vars[key]))
				else:
					if	type(local_vars[key]) == type(set()):
						my_vars[key] = list(local_vars[key])
					else:
						my_vars[key] = local_vars[key]
		elif type(local_vars[key]) == types.DictType:
			if len({type(local_vars[key][i]) for i in local_vars[key]} - set(allowed_types)) == 0:	
					my_vars[key] = local_vars[key]

	ordered_my_vars = collections.OrderedDict(sorted(my_vars.items()))
	
	with open(filename, write_or_append) as fp:
		json.dump(ordered_my_vars, fp, indent = 2)
	fp.close()


def print_program_start_information():
	
	from mpi import MPI_COMM_WORLD, mpi_comm_rank, mpi_comm_size, mpi_barrier
	import os
	from socket import gethostname

	myid = mpi_comm_rank(MPI_COMM_WORLD)
	mpi_size = mpi_comm_size(MPI_COMM_WORLD)	# Total number of processes, passed by --np option.

	if(myid == 0):
		print "Location: " + os.getcwd()
		
	print "MPI Rank: %03d/%03d "%(myid, mpi_size) + "Hostname: " + gethostname() +  " proc_id: " + str(os.getpid())



def store_program_state(filename, state, stack):
	import json
	with open(filename, "w") as fp:
		json.dump(zip(stack, state), fp, indent = 2)
	fp.close()

def restore_program_stack_and_state(file_name_of_saved_state):
	import json; f = open(file_name_of_saved_state, 'r')
	saved_state_and_stack = json.load(f); f.close()
	return list(zip(*saved_state_and_stack)[0]), list(zip(*saved_state_and_stack)[1])


def program_state_stack(full_current_state, frameinfo, file_name_of_saved_state=None, last_call="", force_starting_execution = False):

	"""

	When used it needs: from inspect import currentframe, getframeinfo
	Also: from utilities import program_state_stack

	This function is used for restarting time consuming data processing programs/steps from the last saved point. 

	This static variable must be defined before the first call:
	program_state_stack.PROGRAM_STATE_VARIABLES = {"isac_generation", "i", "j"}
	It contains local variables at any level of the stack that define uniquely the state(flow/logic) of the program.
	
	It is assumed that the processed data is saved at each step and it is independent from the variables that uniquely define 
	the state(flow/logic) of the program. All the variables that are used in more than one step must be calculated before
	the "if program_state_stack(locals(), getframeinfo(currentframe())):" call. It is assumed that they are not time consuming.
	Passing processed data from one step to the next is done only through files. 
	
	First call needs to contain "file_name_of_saved_state".
	Then, the next calls are "if program_state_stack(locals(), getframeinfo(currentframe())):" to demarcate the blocks of 
	processing steps that take a long time (hours/days).
	
	Example of initialization:
	program_state_stack.PROGRAM_STATE_VARIABLES = {"isac_generation", "i", "j"}
	program_state_stack(locals(), getframeinfo(currentframe()), "my_state.json")
	
	Then regular usage in the program:
	
	if program_state_stack(locals(), getframeinfo(currentframe())):
	# if 1:
		pass
	
	"""

	from traceback import extract_stack
	from mpi import mpi_comm_rank, mpi_bcast, MPI_COMM_WORLD, MPI_INT
	from utilities import if_error_then_all_processes_exit_program
	import os

	def get_current_stack_info():
		return [[x[0], x[2]] for x in extract_stack()[:-2]]

	START_EXECUTING_FALSE = 0
	START_EXECUTING_TRUE = 1
	START_EXECUTING_ONLY_ONE_TIME_THEN_REVERT = 2
	
	# error_status = 1
	# if_error_then_all_processes_exit_program(error_status)
	
	current_state = dict()
	for var in program_state_stack.PROGRAM_STATE_VARIABLES & set(full_current_state) :
		current_state[var] =  full_current_state[var]

	if "restart_location_title" in program_state_stack.__dict__:
		# location_in_program = frameinfo.filename + "___" + program_state_stack.restart_location_title + "___" + last_call
		location_in_program = frameinfo.filename + "___" + program_state_stack.restart_location_title
		del program_state_stack.restart_location_title
	else:
		location_in_program = frameinfo.filename + "___" + str(frameinfo.lineno) + "_" + last_call
		# location_in_program = frameinfo.filename + "___" + last_call
		
	current_state["location_in_program"] = location_in_program
	
	current_stack = get_current_stack_info()

	error_status = 0

	# not a real while, an if with the possibility of jumping with break
	while mpi_comm_rank(MPI_COMM_WORLD) == 0:
		if "file_name_of_saved_state" not in program_state_stack.__dict__:
			if type(file_name_of_saved_state) != type(""):
				print "Must provide the file name of saved state as a string in the first call of the function!"
				error_status = 1
				break

			program_state_stack.file_name_of_saved_state = os.getcwd() + os.sep + file_name_of_saved_state
			program_state_stack.counter = 0
			program_state_stack.track_stack = get_current_stack_info()
			program_state_stack.track_state = [dict() for i in xrange(len(program_state_stack.track_stack))]
			program_state_stack.track_state[-1] = current_state

			file_name_of_saved_state_contains_information = False
			if (os.path.exists(file_name_of_saved_state)):
				statinfo = os.stat(file_name_of_saved_state)
				file_name_of_saved_state_contains_information = statinfo.st_size > 0 
			if file_name_of_saved_state_contains_information:
				program_state_stack.saved_stack, \
				program_state_stack.saved_state = restore_program_stack_and_state(file_name_of_saved_state)
				program_state_stack.start_executing = START_EXECUTING_FALSE
			else:
				# check to see if file can be created
				f = open(file_name_of_saved_state, "w"); f.close()
				program_state_stack.start_executing = START_EXECUTING_TRUE
		else:
			program_state_stack.counter += 1
			# print "counter: ", program_state_stack.counter
			# if program_state_stack.counter == program_state_stack.CCC:
			# 	error_status = 1
			# 	break

			if program_state_stack.start_executing == START_EXECUTING_ONLY_ONE_TIME_THEN_REVERT:
				program_state_stack.start_executing = START_EXECUTING_FALSE
			
			# correct track_state to reflect track_stack 
			for i in xrange(len(current_stack)):
				if i < len(program_state_stack.track_state):
					if program_state_stack.track_stack[i] != current_stack[i]:
						program_state_stack.track_state[i] = dict()
				else:
					# print "i:", i, len(program_state_stack.track_state), len(current_stack), current_stack
					program_state_stack.track_state.append(dict())
			program_state_stack.track_state[i] = current_state
			
			# correct track_stack to reflect current_stack
			program_state_stack.track_stack = current_stack
			
			# if program_state_stack.counter == 68:
			# 	print range(len(current_stack), len(program_state_stack.track_state))
				
			# delete additional elements in track_state so that size of track_state is the same as current_stack  				
			program_state_stack.track_state[len(current_stack):len(program_state_stack.track_state)] = []
			
			if program_state_stack.start_executing == START_EXECUTING_TRUE or last_call != "" or force_starting_execution:
				store_program_state(program_state_stack.file_name_of_saved_state, program_state_stack.track_state, current_stack)
				program_state_stack.start_executing = START_EXECUTING_TRUE
			else:
				if len(program_state_stack.saved_state) >= len(current_stack):
					for i in range(len(program_state_stack.saved_state)):
						if i < len(current_stack):
							if program_state_stack.track_stack[i] == current_stack[i]:
								if program_state_stack.track_state[i] == program_state_stack.saved_state[i]:
									continue
							break
						else:
							program_state_stack.start_executing = START_EXECUTING_ONLY_ONE_TIME_THEN_REVERT
							# print "////////////////////////////" 
							# print "Entering function: ", location_in_program
							# print "////////////////////////////"
							break
					else:
						program_state_stack.start_executing = START_EXECUTING_TRUE
						# print "////////////////////////////" 
						# print "Start executing: ", location_in_program
						# print "////////////////////////////"
		break
	else:
		program_state_stack.start_executing = START_EXECUTING_FALSE
		
	if_error_then_all_processes_exit_program(error_status)	
		
	program_state_stack.start_executing = mpi_bcast(program_state_stack.start_executing, 1, MPI_INT, 0, MPI_COMM_WORLD)
	program_state_stack.start_executing = int(program_state_stack.start_executing[0])

	# print "program_state_stack.start_executing ", program_state_stack.start_executing

	return program_state_stack.start_executing

def qw(s):
	s = s.replace("\n"," ")
	s = s.replace("\t"," ")
	return tuple(s.split())



def debug_mpi_barrier(comm):
	from mpi import mpi_barrier, mpi_comm_rank, mpi_bcast
	from traceback import extract_stack
	import sys


	# if mpi_comm_rank(comm) in range(4):
	print "Stack info::0::", extract_stack()[-3:]
	
	sys.stdout.flush()
	sys.stderr.flush()
	return mpi_barrier(comm)


def debug_mpi_bcast(newv, s, t, m, comm):
	from mpi import mpi_comm_rank, mpi_bcast
	from traceback import extract_stack	
	import sys

	
	rrr = mpi_bcast(newv, s, t, m, comm)
	# if mpi_comm_rank(comm) in range(4):
	print "Stack info::0::", extract_stack()[-3:], "****************", newv, "####", rrr 

	sys.stdout.flush()
	sys.stderr.flush()

	# return mpi_bcast(newv, s, t, m, comm)
	return rrr

def print_from_process(process_rank, message):
	from mpi import MPI_COMM_WORLD, mpi_comm_rank, mpi_comm_size, mpi_barrier
	import os, sys
	from socket import gethostname

	myid = mpi_comm_rank(MPI_COMM_WORLD)
	mpi_size = mpi_comm_size(MPI_COMM_WORLD)	# Total number of processes, passed by --np option.

	if(myid == process_rank):
		print "MPI Rank: %03d/%03d "%(myid, mpi_size) + "Hostname: " + gethostname() +  " proc_id: " + str(os.getpid()) +\
			"message:::", message
	sys.stdout.flush()
	
def mpi_exit():
	from mpi import mpi_finalize
	import sys

	mpi_finalize()
	sys.stdout.flush()
	sys.exit()
	
### from sort3d

def get_attr_stack(data_stack,attr_string):
	attr_value_list = []
	for idat in xrange(len(data_stack)):
		attr_value = data_stack[idat].get_attr(attr_string)
		attr_value_list.append(attr_value)
	return attr_value_list
	
def get_sorting_attr_stack(data_stack):
	from utilities import get_params_proj
	attr_value_list = []
	for idat in xrange(len(data_stack)):
		group = data_stack[idat].get_attr("group")
		phi,theta,psi,s2x,s2y=get_params_proj(data_stack[idat],xform = "xform.projection")
		attr_value_list.append([group, phi, theta, psi, s2x, s2y])
	return attr_value_list
	
def get_sorting_params(Tracker,data):
	from mpi import mpi_barrier, MPI_COMM_WORLD
	from utilities import read_text_row,wrap_mpi_bcast,even_angles
	from applications import MPI_start_end
	myid      = Tracker["constants"]["myid"]
	main_node = Tracker["constants"]["main_node"]
	nproc     = Tracker["constants"]["nproc"]
	ndata     = Tracker["total_stack"]
	mpi_comm  = MPI_COMM_WORLD
	if myid == main_node:
		total_attr_value_list = []
		for n in xrange(ndata):
			total_attr_value_list.append([])
	else:
		total_attr_value_list = 0
	for inode in xrange(nproc):
		attr_value_list =get_attr_stack(data,"group")
		attr_value_list =wrap_mpi_bcast(attr_value_list,inode)
		if myid ==main_node:
			image_start,image_end=MPI_start_end(ndata,nproc,inode)
			total_attr_value_list=fill_in_mpi_list(total_attr_value_list,attr_value_list,image_start,image_end)		
		mpi_barrier(MPI_COMM_WORLD)
	total_attr_value_list = wrap_mpi_bcast(total_attr_value_list,main_node)
	return total_attr_value_list
	
def get_sorting_params_refine(Tracker,data,ndata):
	from mpi import mpi_barrier, MPI_COMM_WORLD
	from utilities import read_text_row,wrap_mpi_bcast,even_angles
	from applications import MPI_start_end
	myid      = Tracker["constants"]["myid"]
	main_node = Tracker["constants"]["main_node"]
	nproc     = Tracker["constants"]["nproc"]
	#ndata     = Tracker["total_stack"]
	mpi_comm  = MPI_COMM_WORLD
	if myid == main_node:
		total_attr_value_list = []
		for n in xrange(ndata):
			total_attr_value_list.append([])
	else:
		total_attr_value_list = 0
	for inode in xrange(nproc):
		attr_value_list = get_sorting_attr_stack(data)
		attr_value_list = wrap_mpi_bcast(attr_value_list,inode)
		if myid == main_node:
			image_start,image_end = MPI_start_end(ndata,nproc,inode)
			total_attr_value_list = fill_in_mpi_list(total_attr_value_list, attr_value_list, image_start,image_end)		
		mpi_barrier(MPI_COMM_WORLD)
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
	for index in xrange(index_start, index_end):
		mpi_list[index] = data_list[index-index_start]
	return mpi_list
	
def get_groups_from_partition(partition, initial_ID_list, number_of_groups):
	# sort out Kmref results to individual groups that has initial IDs
	# make a dictionary
	dict = {}
	for iptl in xrange(len(initial_ID_list)):
		dict[iptl] = initial_ID_list[iptl]
	res = []
	for igrp in xrange(number_of_groups):
		class_one = []
		for ipt in xrange(len(partition)):
			if partition[ipt] == igrp:
				orginal_id = dict[ipt]
				class_one.append(orginal_id)
		res.append(class_one)
	return res

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
	
#### Used in the main programm

def sample_down_1D_curve(nxinit, nnxo, pspcurv_nnxo_file):
	shrinkage=float(nnxo)/float(nxinit)
	curv_orgn = read_text_file(pspcurv_nnxo_file)
	new_curv=int(1.5*len(curv_orgn))*[0.0]
	for index in xrange(len(curv_orgn)):
		new_index = int(index/shrinkage)
		fraction  =  index/shrinkage-new_index
		if fraction <=0:
			new_curv[new_index] +=curv_orgn[index]
		else:
			new_curv[new_index]  +=(1.-fraction)*curv_orgn[index]
			new_curv[new_index+1] += fraction*curv_orgn[index]
	return new_curv
	
def get_initial_ID(part_list, full_ID_dict):
	part_initial_id_list = []
	#new_dict = {}
	for iptl in xrange(len(part_list)):
		part_initial_id_list.append(full_ID_dict[part_list[iptl]])
		#new_dict[iptl] = id
	return part_initial_id_list#, new_dict

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
	
def print_upper_triangular_matrix(data_table_dict,N_indep,log_main):
		msg =""
		for i in xrange(N_indep):
			msg +="%7d"%i
		log_main.add(msg)
		for i in xrange(N_indep):
			msg ="%5d "%i
			for j in xrange(N_indep):
				if i<j:
					msg +="%5.2f "%data_table_dict[(i,j)]
				else:
					msg +="      "
			log_main.add(msg)
			
def print_a_line_with_timestamp(string_to_be_printed ):                 
	line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
 	print(line,string_to_be_printed)
	return string_to_be_printed

def convertasi(asig,K):
	from numpy import array
	p = []
	for k in xrange(K):
		l = []
		for i in xrange(len(asig)):
			if( asig[i ]== k ): l.append(i)
		l = array(l,"int32")
		l.sort()
		p.append(l)
	return p

def prepare_ptp(data_list, K):
	num_of_pt = len(data_list)
	ptp=[]
	for ipt in xrange(num_of_pt):
		ptp.append([])
	for ipt in xrange(num_of_pt):
		nc = len(data_list[ipt])
		asig  =[-1]*nc
		for i in xrange(nc):
			asig[i] = data_list[ipt][i]
		ptp[ipt] = convertasi(asig, K)
	return ptp

def print_dict(dict,theme):
		line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
		print(line+theme)
		spaces = "                           "
		for key, value in sorted( dict.items() ):
			if(key != "constants"):  
				print("                    => "+key+spaces[len(key):]+":  "+str(value))
"""
def checkstep(item, keepchecking, myid, main_node):
	from utilities import bcast_number_to_all
        if(myid == main_node):
                if keepchecking:
                        if(os.path.exists(item)):
                                doit = 0
                        else:
                                doit = 1
                                keepchecking = False
                else:
                        doit = 1
        else:
                doit = 1
        doit = bcast_number_to_all(doit, source_node = main_node)
        return doit, keepchecking
"""

def get_resolution_mrk01(vol, radi, nnxo, fscoutputdir, mask_option):
        # this function is single processor
        #  Get updated FSC curves, user can also provide a mask using radi variable
	import types
	from statistics import fsc
	from utilities import model_circle, get_im
	from filter import fit_tanh1
	import os
	if(type(radi) == int):
		if(mask_option is None):  mask = model_circle(radi,nnxo,nnxo,nnxo)
		else:                     mask = get_im(mask_option)
	else:  mask = radi
	nfsc = fsc(vol[0]*mask,vol[1]*mask, 1.0,os.path.join(fscoutputdir,"fsc.txt") )
	currentres = -1.0
	ns = len(nfsc[1])
	#  This is actual resolution, as computed by 2*f/(1+f)
	for i in xrange(1,ns-1):
		if ( nfsc[1][i] < 0.333333333333333333333333):
			currentres = nfsc[0][i-1]
			break
		#if(currentres < 0.0):
			#print("  Something wrong with the resolution, cannot continue")
		currentres = nfsc[0][i-1]
        
        """ this commented previously
		lowpass = 0.5
		ns = len(nfsc[1])
        #  This is resolution used to filter half-volumes
        for i in xrange(1,ns-1):
                if ( nfsc[1][i] < 0.5 ):
                        lowpass = nfsc[0][i-1]
                        break
        """  
	lowpass, falloff = fit_tanh1(nfsc, 0.01)
	return  round(lowpass,4), round(falloff,4), round(currentres,2)
        
def partition_to_groups(alist, K):
	res =[]
	for igroup in xrange(K):
		this_group =[]
		for imeb in xrange(len(alist)):
			if( alist[imeb] == igroup ):   this_group.append(imeb)
		this_group.sort()
		res.append(this_group)
	return res

def partition_independent_runs(run_list, K):
	indep_runs_groups = {}
	for indep in xrange(len(run_list)):
		indep_runs_groups[indep] = partition_to_groups(run_list[indep], K)
	return indep_runs_groups

def get_outliers(total_number,plist):
	tlist={}
	for i in xrange(total_number):tlist[i]=i
	for a in plist:   del tlist[a]
	out =[]
	for a in tlist:   out.append(a)
	return out

def merge_groups(stable_members_list):
	alist=[]
	for i in xrange(len(stable_members_list)):
		for j in xrange(len(stable_members_list[i])):alist.append(stable_members_list[i][j])
	return alist

def save_alist(Tracker,name_of_the_text_file,alist):
	from utilities import write_text_file
	import os
	log       =Tracker["constants"]["log_main"]
	myid      =Tracker["constants"]["myid"]
	main_node =Tracker["constants"]["main_node"]
	dir_to_save_list =Tracker["this_dir"]
	if myid==main_node:
		file_name=os.path.join(dir_to_save_list,name_of_the_text_file)
		write_text_file(alist, file_name)

def margin_of_error(P, size_of_this_sampling):
	# margin of an error, or radius of an error for a percentage
	from math import sqrt
	return sqrt(P*(1.-P)/size_of_this_sampling)
	
def get_margin_of_error(this_group_of_data,Tracker):
	ratio = margin_of_error(Tracker["P_chunk0"],len(this_group_of_data))
	rate1, rate2, size_of_this_sampling = count_chunk_members(Tracker["chunk_dict"],this_group_of_data)
	return abs(rate1-Tracker["P_chunk0"]),ratio,abs(rate2-Tracker["P_chunk1"]),ratio
	
def do_two_way_comparison(Tracker):
	from mpi import mpi_barrier, MPI_COMM_WORLD
	from utilities import read_text_file,write_text_file
	from statistics import k_means_match_clusters_asg_new
	import os
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
		from mpi import mpi_finalize
		from sys import exit
		mpi_finalize()
		exit()
	else:
		for iter_indep in xrange(Tracker["constants"]["indep_runs"]):  total_partition.append(Tracker["partition_dict"][iter_indep])
		### Two-way comparision is carried out on all nodes 
		ptp = prepare_ptp(total_partition, number_of_groups)
		indep_runs_to_groups = partition_independent_runs(total_partition, number_of_groups)
		###### Check margin of error
		if myid ==main_node:
			log_main.add("--------------------------margin of error--------------------------------------------")
		for indep in xrange(len(indep_runs_to_groups)):
			for index_of_class in xrange(len(indep_runs_to_groups[indep])):
				one_group_in_old_ID = get_initial_ID(indep_runs_to_groups[indep][index_of_class], Tracker["full_ID_dict"])
				#if myid ==main_node:
				#	print  " chunk_dict ",Tracker["chunk_dict"]
				#	print "one_group_in_old_ID",one_group_in_old_ID
				rate1, rate2, size_of_this_group = count_chunk_members(Tracker["chunk_dict"], one_group_in_old_ID)
				error = margin_of_error(Tracker["P_chunk0"], size_of_this_group)
				if myid ==main_node:
					#log_main.add("    %f     %f     %d"%(rate1,rate2,size_of_this_group))
					log_main.add(" margin of error for chunk0 is %f   %f    %d"%((Tracker["P_chunk0"]-error),(Tracker["P_chunk0"]+error),size_of_this_group))
					log_main.add(" actual percentage is %f"%rate1)
					#log_main.add(" margin of error for chunk1 is %f"%margin_of_error(Tracker["P_chunk1"],size_of_this_group))
					#log_main.add(" actual error is %f"%abs(rate2-Tracker["P_chunk1"]))
		if myid ==main_node:
			log_main.add("------------------------------------------------------------------------------")
		total_pop=0
		two_ways_stable_member_list = {}
		avg_two_ways                = 0.0
		avg_two_ways_square         = 0.0
		scores                      = {}
		for iptp in xrange(len(ptp)):
			for jptp in xrange(len(ptp)):
				newindeces, list_stable, nb_tot_objs = k_means_match_clusters_asg_new(ptp[iptp], ptp[jptp])
				tt = 0.0
				if myid ==main_node and iptp<jptp:
					aline="Two-way comparison between independent run %3d and %3d"%(iptp,jptp)
					log_main.add(aline)
				for m in xrange(len(list_stable)):
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
		for ipp in xrange(len(ptp)):
			avg_scores =0.0
			for jpp in xrange(len(ptp)):
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
		for istable in xrange(len(Tracker["two_way_stable_member"])):
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
		mpi_barrier(MPI_COMM_WORLD)
		save_alist(Tracker,"Accounted.txt", accounted_list)
		update_full_dict(accounted_list,Tracker)# Update full_ID_dict for Kmeans
		mpi_barrier(MPI_COMM_WORLD)
		Tracker["this_unaccounted_dir"]     = workdir
		Tracker["this_unaccounted_text"]    = os.path.join(workdir,"Unaccounted.txt")
		Tracker["this_accounted_text"]      = os.path.join(workdir,"Accounted.txt")
		Tracker["ali3d_of_outliers"]        = os.path.join(workdir,"ali3d_params_of_outliers.txt")
		Tracker["ali3d_of_accounted"]       = os.path.join(workdir,"ali3d_params_of_accounted.txt")
		if myid==main_node:
			log_main.add(" Selected indepedent runs      %5d and  %5d"%(run1,run2))
			log_main.add(" Their pair-wise averaged rates are %5.2f  and %5.2f "%(rate1,rate2))		
		from math import sqrt
		avg_two_ways        = avg_two_ways/total_pop
		two_ways_std        = sqrt(avg_two_ways_square/total_pop-avg_two_ways**2)
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
		mpi_barrier(MPI_COMM_WORLD)

def select_two_runs(summed_scores,two_way_dict):
	summed_scores.sort()
	rate1 = summed_scores[-1]
	rate2 = None
	for index in xrange(2,len(summed_scores)+1):
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

def get_ali3d_params(ali3d_old_text_file,shuffled_list):
	from utilities import read_text_row
	ali3d_old = read_text_row(ali3d_old_text_file)
	ali3d_new = []
	for iptl in xrange(len(shuffled_list)):
		ali3d_new.append(ali3d_old[shuffled_list[iptl]])
	return ali3d_new

def counting_projections(delta, ali3d_params, image_start):
	from utilities import even_angles,angle_between_projections_directions
	sampled_directions = {}
	angles=even_angles(delta,0,180)
	for a in angles:
		[phi0, theta0, psi0]=a
		sampled_directions[(phi0,theta0)]=[]
	from math import sqrt
	for i in xrange(len(ali3d_params)):
		[phi, theta, psi, s2x, s2y] = ali3d_params[i]
		dis_min    = 9999.
		this_phi   = 9999.
		this_theta = 9999.
		this_psi   = 9999.
		prj1       =[phi,theta]
		for j in xrange(len(angles)):
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
			for item in xrange(2,len(ang_proj)):
				dict_angle_main_node[(ang_proj[0],ang_proj[1])].append(item)
	return dict_angle_main_node

def get_stat_proj(Tracker,delta,this_ali3d):
	from mpi import mpi_barrier, MPI_COMM_WORLD
	from utilities import read_text_row,wrap_mpi_bcast,even_angles
	from applications import MPI_start_end
	myid      = Tracker["constants"]["myid"]
	main_node = Tracker["constants"]["main_node"]
	nproc     = Tracker["constants"]["nproc"]
	mpi_comm  = MPI_COMM_WORLD
	if myid ==main_node:
		ali3d_params=read_text_row(this_ali3d)
		lpartids    = range(len(ali3d_params))
	else:
		lpartids      = 0
		ali3d_params  = 0
	lpartids = wrap_mpi_bcast(lpartids, main_node)
	ali3d_params = wrap_mpi_bcast(ali3d_params, main_node)
	ndata=len(ali3d_params)
	image_start, image_end = MPI_start_end(ndata, nproc, myid)
	ali3d_params=ali3d_params[image_start:image_end]
	sampled=counting_projections(delta,ali3d_params,image_start)
	for inode in xrange(nproc):
		if myid ==inode:
			dlist=unload_dict(sampled)
		else:
			dlist =0
		dlist=wrap_mpi_bcast(dlist,inode)
		if myid ==main_node and inode != main_node:
			sampled=load_dict(sampled,dlist)
		mpi_barrier(MPI_COMM_WORLD)
	return sampled
	
def create_random_list(Tracker):
	import copy
	import random
	from utilities import wrap_mpi_bcast

	myid        = Tracker["constants"]["myid"]
	main_node   = Tracker["constants"]["main_node"]
	total_stack = Tracker["total_stack"]

	if Tracker["constants"]["seed"] ==- 1: random.seed()
	else:                                  random.seed(Tracker["constants"]["seed"])

	indep_list  = []
	for irandom in xrange(Tracker["constants"]["indep_runs"]):
		ll = copy.copy(Tracker["this_data_list"])
		random.shuffle(ll)
		ll = wrap_mpi_bcast(ll, main_node)
		indep_list.append(ll)
	Tracker["this_indep_list"] = indep_list
	
def get_number_of_groups(total_particles,number_of_images_per_group, round_off=.2):
	number_of_groups=float(total_particles)/number_of_images_per_group
	if number_of_groups - int(number_of_groups)<round_off:
		number_of_groups = int(number_of_groups)
	else:
		number_of_groups = int(number_of_groups)+1
	return number_of_groups

def recons_mref(Tracker):
	from mpi import mpi_barrier, MPI_COMM_WORLD
	import os
	from time import sleep
	from reconstruction import recons3d_4nn_ctf_MPI
	from utilities import get_shrink_data_huang
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
	for igrp in xrange(number_of_groups):
		a_group_list = particle_list[(total_data*igrp)//number_of_groups:(total_data*(igrp+1))//number_of_groups]
		a_group_list.sort()
		Tracker["this_data_list"] = a_group_list
		from utilities import write_text_file
		particle_list_file = os.path.join(Tracker["this_dir"], "iclass%d.txt"%igrp)
		if myid ==main_node:
			write_text_file(Tracker["this_data_list"],particle_list_file)
		mpi_barrier(MPI_COMM_WORLD)
		data, old_shifts =  get_shrink_data_huang(Tracker,nxinit,particle_list_file,partstack,myid,main_node,nproc,preshift=True)
		#vol=reconstruct_3D(Tracker,data)
		mpi_barrier(MPI_COMM_WORLD)
		vol = recons3d_4nn_ctf_MPI(myid=myid,prjlist=data,symmetry=Tracker["constants"]["sym"],finfo=None)
		if myid ==main_node:
			print "reconstructed %3d"%igrp
		ref_list.append(vol)
		number_of_ref_class.append(len(Tracker["this_data_list"]))
	Tracker["number_of_ref_class"] = number_of_ref_class
	return ref_list

def apply_low_pass_filter(refvol,Tracker):
	from filter import filt_tanl
	for iref in xrange(len(refvol)):
		refvol[iref]=filt_tanl(refvol[iref],Tracker["low_pass_filter"],.1)
	return refvol
	
def get_groups_from_partition(partition, initial_ID_list, number_of_groups):
	# sort out Kmref results to individual groups that has initial IDs
	# make a dictionary
	dict = {}
	for iptl in xrange(len(initial_ID_list)):
		dict[iptl] = initial_ID_list[iptl]
	res = []
	for igrp in xrange(number_of_groups):
		class_one = []
		for ipt in xrange(len(partition)):
			if partition[ipt] == igrp:
				orginal_id = dict[ipt]
				class_one.append(orginal_id)
		res.append(class_one)
	return res

def get_number_of_groups(total_particles,number_of_images_per_group):
	number_of_groups=float(total_particles)/number_of_images_per_group
	if number_of_groups - int(number_of_groups)<.4:
		number_of_groups = int(number_of_groups)
	else:
		number_of_groups = int(number_of_groups)+1
	return number_of_groups
	
def get_complementary_elements(total_list,sub_data_list):
	if len(total_list)<len(sub_data_list):
		print "Wrong input list!"
		return []
	else:
		sub_data_dict     = {}
		complementary     = []
		for index in xrange(len(sub_data_list)):sub_data_dict[sub_data_list[index]]=index
		for any in total_list:
			if sub_data_dict.has_key(any) is False:complementary.append(any)
		return complementary

def get_complementary_elements_total(total_stack, data_list):
	data_dict    ={}
	complementary     = []
	for index in xrange(len(data_list)):data_dict[data_list[index]]=index
	for index in xrange(total_stack):
		if data_dict.has_key(index) is False:complementary.append(index)
	return complementary

def update_full_dict(leftover_list, Tracker):
	full_dict = {}
	for iptl in xrange(len(leftover_list)):
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
	
def get_two_chunks_from_stack(Tracker):
	total_chunk = EMUtil.get_all_attributes(Tracker["orgstack"],"chunk_id")
	chunk_one = []
	chunk_two = []
	for index_of_total_chunk in xrange(len(total_chunk)):
		if total_chunk[index_of_total_chunk]==0:chunk_one.append(index_of_total_chunk)
		else:chunk_two.append(index_of_total_chunk)
	return chunk_one, chunk_two

def adjust_fsc_down(fsc,n1,n2):
	# fsc curve:  frequencies   cc values  number of the sampling points
	# n1 total data n2 subset
	from utilities import read_text_file
	import types
	if type(fsc) == types.StringType:fsc=read_text_file(fsc,-1)
	N_bins =  len(fsc[0])
	adjusted_fsc = N_bins*[None]
	for index_of_cc in xrange(N_bins):
		adjusted_fsc[index_of_cc] = (float(n2)/float(n1))*fsc[1][index_of_cc]/(1.-(1.-float(n2)/float(n1))*fsc[1][index_of_cc])
	calibrated_fsc=[fsc[0], adjusted_fsc, fsc[2]]
	return calibrated_fsc
	
def set_filter_parameters_from_adjusted_fsc(n1,n2,Tracker):
	fsc_cutoff   = 1.0/3.0
	adjusted_fsc = adjust_fsc_down(Tracker["global_fsc"],n1,n2)
	currentres   = -1.0
	ns           = len(adjusted_fsc)
	for i in xrange(1,ns-1):
		if adjusted_fsc[1][i] < fsc_cutoff:
			currentres = adjusted_fsc[0][i-1]
			break
	lowpass, falloff    = fit_tanh1(adjusted_fsc, 0.01)
	lowpass             = round(lowpass,4)
	falloff    =min(.1,falloff)
	falloff             = round(falloff,4)
	currentres          = round(currentres,2)	
	Tracker["lowpass"]  = lowpass
	Tracker["falloff"]  = falloff
##### from RSORT

def get_class_members(sort3d_dir):
	import os
	from utilities import read_text_file
	maximum_generations = 100
	maximum_groups      = 100
	class_list = []
	igen = 0
	while( igen < maximum_generations ):
		gendir = os.path.join(sort3d_dir, "generation%03d"%igen)
		if os.path.exists(gendir):
			igrp = 0
			while( igrp < maximum_groups):
				Class_file = os.path.join(gendir, "Kmref/Class%d.txt"%igrp)
				if os.path.exists(Class_file):
					class_one = read_text_file(Class_file)
					class_list.append(class_one)
					igrp += 1
				else:
					igrp = maximum_groups
			igen += 1
		else:
			igen = maximum_generations
	return class_list

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
	number_of_groups=float(total_particles)/number_of_images_per_group
	if number_of_groups - int(number_of_groups)<.4:number_of_groups = int(number_of_groups)
	else:number_of_groups = int(number_of_groups)+1
	return number_of_groups
	
def get_stable_members_from_two_runs(SORT3D_rootdirs, ad_hoc_number, log_main):
	#SORT3D_rootdirs                       =sys.argv[1]
	# ad_hoc_number would be a number larger than the id simply for handling two_way comparison of non-equal number of groups from two partitions.
	########
	from string import split
	from statistics import k_means_match_clusters_asg_new
	from numpy import array
	
	sort3d_rootdir_list = split(SORT3D_rootdirs)
	dict1              = []
	maximum_elements   = 0
	for index_sort3d in xrange(len(sort3d_rootdir_list)):
		sort3d_dir       = sort3d_rootdir_list[index_sort3d]
		all_groups       = get_class_members(sort3d_dir)
		dict1.append(all_groups)
		if maximum_elements <len(all_groups):
			maximum_elements = len(all_groups)
	TC = ad_hoc_number + 1
	for indep in xrange(len(dict1)):
		alist = dict1[indep] 
		while len(alist)<maximum_elements:
			alist.append([TC])
			TC += 1
		dict1[indep] = alist
		TC += 1
	for a in dict1:   log_main.add(len(a))
	dict = {}
	for index_sort3d in xrange(len(sort3d_rootdir_list)):
		sort3d_dir       = sort3d_rootdir_list[index_sort3d]
		dict[sort3d_dir] = dict1[index_sort3d]
	###### Conduct two-way comparison
	for isort3d in xrange(0,1): #len(sort3d_rootdir_list)):
		li = dict[sort3d_rootdir_list[isort3d]]
		new_li = []
		for ili in xrange(len(li)):
			li[ili].sort()
			t= array(li[ili],'int32')
			new_li.append(t)
		avg_list = {}
		total    = {}
		for ii in xrange(len(li)):
			avg_list[ii]=0.0
			total[ii]=0.0
		for jsort3d in xrange(len(sort3d_rootdir_list)):
			if isort3d != jsort3d:
				new_lj = []
				lj = dict[sort3d_rootdir_list[jsort3d]]
				for a in lj:
					log_main.add("the size is  %d"%len(a))
				for jlj in xrange(len(lj)):
					lj[jlj].sort()
					t= array(lj[jlj],'int32')
					new_lj.append(t)
				ptp=[new_li,new_lj]
				newindeces, list_stable, nb_tot_objs = k_means_match_clusters_asg_new(ptp[0],ptp[1])
				log_main.add("*************************************************************")
				log_main.add("the results of two P1 runs are: ")
				for index in xrange(len(newindeces)):
					log_main.add("  %d of %s matches  %d of %s"%(newindeces[index][0],sort3d_rootdir_list[isort3d],newindeces[index][1],sort3d_rootdir_list[jsort3d]))
				for index in xrange(len(list_stable)):
					log_main.add("%d   stable memebers"%len(list_stable[index]))
				new_stable = []
				for ilist in xrange(len(list_stable)):
					if len(list_stable[ilist])!= 0:
						new_stable.append(list_stable[ilist])
				for istable in xrange(len(new_stable)):
					stable = new_stable[istable]
					if len(stable)>0: 
						group_A =  li[newindeces[istable][0]]
						group_B =  lj[newindeces[istable][1]]
						log_main.add(" %d %d %d   "%(len(group_A),len(group_B),len(stable)))
		return new_stable
		
def two_way_comparison_single(partition_A, partition_B,Tracker):
	###############
	from statistics import k_means_match_clusters_asg_new
	from utilities import count_chunk_members, margin_of_error
	from numpy import array
	#two_way_comparison_single
	total_stack = Tracker["constants"]["total_stack"]
	log_main    = Tracker["constants"]["log_main"]
	myid        = Tracker["constants"]["myid"]
	main_node   = Tracker["constants"]["main_node"]
	numpy32_A   = []
	numpy32_B   = []
	total_A     = 0
	total_B     = 0
	###--------------------------------------
	
	if myid == main_node:
		log_main.add(" the first run has number of particles %d"%len(partition_A))
		log_main.add(" the second run has number of particles %d"%len(partition_B))
	for A in partition_A:
		total_A +=len(A)
	for B in partition_B:
		total_B +=len(B)
	nc_zero = 1
	if len(partition_A) < len(partition_B):
		while len(partition_A) <len(partition_B):
			partition_A.append([nc_zero+total_stack])
			nc_zero +=1
	elif len(partition_A) > len(partition_B):
		while len(partition_B) <len(partition_A):
			partition_B.append([nc_zero+total_stack])
			nc_zero +=1
	number_of_class = len(partition_A)
	for index_of_class in xrange(number_of_class):
		A = partition_A[index_of_class]
		A.sort()
		A= array(A,'int32')
		numpy32_A.append(A)
		B= partition_B[index_of_class]
		B.sort()
		B = array(B,'int32')
		numpy32_B.append(B)
		if myid ==main_node:
			log_main.add("group %d  %d   %d"%(index_of_class,len(A), len(B))) 
	ptp    = [[],[]]
	ptp[0] = numpy32_A
	ptp[1] = numpy32_B
	newindexes, list_stable, nb_tot_objs = k_means_match_clusters_asg_new(ptp[0],ptp[1])
	if myid == main_node:
		log_main.add(" reproducible percentage of the first partition %f"%(nb_tot_objs/float(total_A)*100.))
		log_main.add(" reproducible percentage of the second partition %f"%(nb_tot_objs/float(total_B)*100.))
		for index in xrange(len(newindexes)):
			log_main.add("%d of A match %d of B "%(newindexes[index][0],newindexes[index][1]))
		for index in xrange(len(list_stable)):
			log_main.add("%d number of reproduced objects are found in group %d"%(len(list_stable[index]),index))
		log_main.add(" %d number of objects are reproduced "%nb_tot_objs)
		log_main.add(" margin of error")
	large_stable = []
	for index_of_stable in xrange(len(list_stable)):
		rate1,rate2,size_of_this_group = count_chunk_members(Tracker["chunk_dict"], list_stable[index_of_stable])
		if size_of_this_group>=Tracker["constants"]["smallest_group"]:
			error                          = margin_of_error(Tracker["P_chunk0"],size_of_this_group)
			if myid == main_node:
				log_main.add(" chunk0  lower bound %f  upper bound  %f  for sample size  %d     group id %d"%((Tracker["P_chunk0"]- error),(Tracker["P_chunk0"]+error),size_of_this_group, index_of_stable))
				log_main.add(" actual percentage is %f"%rate1)
			large_stable.append(list_stable[index_of_stable])
		else:
			if myid==main_node:
				log_main.add("%d  group is too small"%index_of_stable)
	return large_stable
	
def get_leftover_from_stable(stable_list, N_total,smallest_group):
	tmp_dict={}
	for i in xrange(N_total):
		tmp_dict[i]=i
	new_stable =[]
	for alist in stable_list:
		if len(alist) > smallest_group:
			for index_of_list in xrange(len(alist)):
				del tmp_dict[alist[index_of_list]]
			new_stable.append(alist)
	leftover_list = []
	for one_element in tmp_dict:
		leftover_list.append(one_element)
	return leftover_list, new_stable
		
def Kmeans_exhaustive_run(ref_vol_list,Tracker):
	from applications import ali3d_mref_Kmeans_MPI
	from utilities import write_text_file
	from reconstruction import rec3D_two_chunks_MPI
	from morphology import get_shrink_3dmask
	from utilities import wrap_mpi_bcast
	import os
	from mpi import MPI_COMM_WORLD, mpi_barrier
	# npad 2 ---------------------------------------
	npad                  = 2
	myid                  = Tracker["constants"]["myid"]
	main_node             = Tracker["constants"]["main_node"]
	log_main              = Tracker["constants"]["log_main"]
	nproc                 = Tracker["constants"]["nproc"]
	final_list_text_file  = Tracker["this_data_list_file"] ## id text file for get_shrink_data_huang
	snr  =1.
	Tracker["total_stack"]= len(Tracker["this_data_list"])
	if myid ==main_node:
		log_main.add("start exhaustive Kmeans")
		log_main.add("total data is %d"%len(Tracker["this_data_list"]))
		log_main.add("final list file is "+final_list_text_file)
	workdir = Tracker["this_dir"]
	####----------------------------------------------
	empty_group = 1
	kmref =0
	while empty_group ==1 and kmref<=5:## In case pctn of Kmeans jumps between 100% to 0%, stop the program
		if myid ==main_node:
			log_main.add(" %d     Kmref run"%kmref) 
		outdir =os.path.join(workdir, "Kmref%d"%kmref)
		empty_group, res_classes, data_list = ali3d_mref_Kmeans_MPI(ref_vol_list, outdir, final_list_text_file, Tracker)
		kmref +=1
		if empty_group ==1:
			if myid ==main_node:
				log_main.add("empty gorup appears, next round of Kmeans requires rebuilding reference volumes!")
				log_main.add(" the number of classes for next round before cleaning is %d"%len(res_classes))
			final_list   = []
			new_class    = []
			for a in res_classes:
				if len(a)>=Tracker["constants"]["smallest_group"]:
					for b in a:
						final_list.append(b)
					new_class.append(a)
			final_list.sort()
			### reset variables of Kmeans run 
			Tracker["total_stack"]    = len(final_list)
			Tracker["this_data_list"] = final_list
			final_list_text_file      = os.path.join(workdir, "final_list%d.txt"%kmref)
			if myid == main_node:
				log_main.add("number of classes for next round is %d"%len(new_class))
				write_text_file(final_list, final_list_text_file)
			mpi_barrier(MPI_COMM_WORLD)
			
			if myid == main_node:
				number_of_ref_class = []
				for igrp in xrange(len(new_class)):
					write_text_file(new_class[igrp],os.path.join(workdir,"final_class%d.txt"%igrp))
					number_of_ref_class.append(len(new_class[igrp]))
			else:   number_of_ref_class = 0
			number_of_ref_class = wrap_mpi_bcast(number_of_ref_class,main_node)
			mpi_barrier(MPI_COMM_WORLD)
			
			ref_vol_list = []
			if  Tracker["constants"]["mask3D"]: mask3D = get_shrink_3dmask(Tracker["constants"]["nxinit"],Tracker["constants"]["mask3D"])
			else: mask3D = None
			Tracker["number_of_ref_class"] = number_of_ref_class
			for igrp in xrange(len(new_class)):
				data,old_shifts = get_shrink_data_huang(Tracker,Tracker["nxinit"],os.path.join(workdir,"final_class%d.txt"%igrp),Tracker["constants"]["partstack"],myid,main_node,nproc,preshift = True)
				#volref = recons3d_4nn_ctf_MPI(myid=myid, prjlist = data, symmetry=Tracker["constants"]["sym"], finfo=None)
				#volref = filt_tanl(volref, Tracker["low_pass_filter"],.1)
				volref, fsc_kmref = rec3D_two_chunks_MPI(data,snr,Tracker["constants"]["sym"],mask3D,\
			 os.path.join(outdir, "resolution_%02d_Kmref%04d"%(igrp,kmref)), myid, main_node, index=-1, npad=npad, finfo = None)
			 	if myid !=main_node:
			 		volref = model_blank(Tracker["nxinit"], Tracker["nxinit"], Tracker["nxinit"])
			 	bcast_EMData_to_all(volref, myid, main_node, MPI_COMM_WORLD)
				ref_vol_list.append(volref)
				mpi_barrier(MPI_COMM_WORLD)
		else:
			new_class    = []
			for a in res_classes:
				if len(a)>=Tracker["constants"]["smallest_group"]:new_class.append(a)
	if myid==main_node:
		log_main.add("Exhaustive Kmeans ends")
		log_main.add(" %d groups are selected out"%len(new_class))
	return new_class
	
def print_a_line_with_timestamp(string_to_be_printed ):                 
	line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
 	print(line,string_to_be_printed)
	return string_to_be_printed
		
def split_a_group(workdir,list_of_a_group,Tracker):
	### Using EQ-Kmeans and Kmeans to split a group
	from utilities import wrap_mpi_bcast
	from random import shuffle
	from mpi import MPI_COMM_WORLD, mpi_barrier
	from utilities import get_shrink_data_huang 
	from reconstructions import recons3d_4nn_ctf_MPI
	from filter import filt_tanl
	from applications import mref_ali3d_EQ_Kmeans
	################
	myid        = Tracker["constants"]["myid"]
	main_node   = Tracker["constants"]["main_node"]
	nproc       = Tracker["constants"]["nproc"]
	total_stack = len(list_of_a_group)
	################
	import copy
	data_list = copy.deepcopy(list_of_a_group)
	update_full_dict(data_list,Tracker)
	this_particle_text_file = os.path.join(workdir,"full_class.txt")
	if myid ==main_node:
		write_text_file(data_list,"full_class.txt")
	# Compute the resolution of leftover 
	if myid ==main_node:
		shuffle(data_list)
		l1=data_list[0:total_stack//2]
		l2=data_list[total_stack//2:]
		l1.sort()
		l2.sort()
	else:
		l1 = 0
		l2 = 0
	l1 = wrap_mpi_bcast(l1, main_node)
	l2 = wrap_mpi_bcast(l2, main_node)
	llist = []
	llist.append(l1)
	llist.append(l2)
	if myid ==main_node:
		for index in xrange(2): 
			partids = os.path.join(workdir,"Class_%d.txt"%index)
			write_text_file(llist[index],partids)
	mpi_barrier(MPI_COMM_WORLD)
	################ create references for EQ-Kmeans
	ref_list = []
	for index in xrange(2):
		partids = os.path.join(workdir,"Class_%d.txt"%index)
		while not os.path.exists(partids):
			#print  " my_id",myid
			sleep(2)
		mpi_barrier(MPI_COMM_WORLD)
		data,old_shifts = get_shrink_data_huang(Tracker,Tracker["constants"]["nxinit"],partids,Tracker["constants"]["partstack"],myid,main_node,nproc,preshift = True)
		vol = recons3d_4nn_ctf_MPI(myid=myid,prjlist = data,symmetry=Tracker["constants"]["sym"],finfo=None)
		vol = filt_tanl(vol,Tracker["constants"]["low_pass_filter"],.1)
		ref_list.append(vol)
	mpi_barrier(MPI_COMM_WORLD)
	### EQ-Kmeans
	outdir = os.path.join(workdir,"EQ-Kmeans")
	mref_ali3d_EQ_Kmeans(ref_list,outdir,this_particle_text_file,Tracker)
	res_EQ = partition_to_groups(Tracker["this_partition"],K=2)
	new_class = []
	for index in xrange(len(res_EQ)):
		new_ID = get_initial_ID(res_EQ(index), Tracker["full_ID_dict"])
		new_class.append(new_ID)
		if myid ==main_node:
			new_class_file = os.path.join(workdir,"new_class%d.txt"%index)
			write_text_file(new_ID,new_class_file)
	mpi_barrier(MPI_COMM_WORLD)
	############# create references for Kmeans
	ref_list = []
	for index in xrange(2):
		partids = os.path.join(workdir,"new_class%d.txt"%index)
		while not os.path.exists(partids):
			#print  " my_id",myid
			sleep(2)
		mpi_barrier(MPI_COMM_WORLD)
		data,old_shifts = get_shrink_data_huang(Tracker,Tracker["constants"]["nxinit"],partids,Tracker["constants"]["partstack"],myid,main_node,nproc,preshift = True)
		vol = recons3d_4nn_ctf_MPI(myid=myid,prjlist = data,symmetry=Tracker["constants"]["sym"],finfo=None)
		vol = filt_tanl(vol,Tracker["constants"]["low_pass_filter"],.1)
		ref_list.append(vol)
	mpi_barrier(MPI_COMM_WORLD)
	#### Kmeans
	
def search_lowpass(fsc):
	fcutoff =.5
	for i in xrange(len(fsc[1])):
		if fsc[0][i]<.5:
			break
	if i<len(fsc[1])-1:
		fcutoff=fsc[0][i-1]
	else:
		fcutoff=.5
	fcutoff=min(.45,fcutoff)
	return fcutoff
#####---------------------------------------------------