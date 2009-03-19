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
#
from EMAN2_cppwrap import *
from global_def import *

def project(volume, params, radius):
        # angles phi, theta, psi
	from fundamentals import rot_shift2D
	from utilities import set_params_proj
        myparams = {"transform":Transform({"type":"spider","phi":params[0],"theta":params[1],"psi":params[2]}), "radius":radius}
        proj = volume.project("pawel", myparams)
	if(params[3]!=0. or params[4]!=0.): 
		params2 = {"filter_type" : Processor.fourier_filter_types.SHIFT, "x_shift" : params[3], "y_shift" : params[4], "z_shift" : 0.0}
		proj=Processor.EMFourierFilter(proj, params2)
		#proj = rot_shift2D(proj, sx = params[3], sy = params[4], interpolation_method = "linear")
	set_params_proj(proj, [params[0], params[1], params[2], -params[3], -params[4]])
	proj.set_attr_dict({'active':1, 'ctf_applied':0})
	return  proj

"""
Temporarily disabled as list cannot be passed to projector.
def prl(vol, params, radius, stack = None):
	from fundamentals import rot_shift2D
	for i in xrange(len(params)):
        	myparams = {"angletype":"SPIDER", "anglelist":params[i][0:3], "radius":radius}
        	proj = vol.project("pawel", myparams)
		if(params[i][3]!=0. or params[i][4]!=0.): proj = rot_shift2D(proj, sx = params[i][3], sy = params[i][4], interpolation_method = "linear")
		proj.set_attr_dict({'phi':params[i][0], 'theta':params[i][1], 'psi':params[i][2], 's2x':-params[i][3], 's2y':-params[i][4]})
		proj.set_attr_dict({'active':1, 'ctf_applied':0})
		if(stack):
			proj.write_image(stack, i)
		else:
			if(i == 0): out= []
			out.append(proj)
	if(stack): return
	else:      return out
"""
def prj(vol, params, stack = None):
	from utilities import set_params_proj
	volft,kb = prep_vol(vol)
	for i in xrange(len(params)):
		proj = prgs(volft, kb, params[i])
		set_params_proj(proj, [params[i][0], params[i][1], params[i][2], -params[i][3], -params[i][4]])
		proj.set_attr_dict({'active':1, 'ctf_applied':0})
		if(stack):
			proj.write_image(stack, i)
		else:
			if(i == 0): out= []
			out.append(proj)
	if(stack):  return
	else:       return out

def prgs(volft, kb, params):
	#  params:  phi, theta, psi, sx, sy
	from fundamentals import fft
	from utilities import set_params_proj
	R = Transform({"type":"spider", "phi":params[0], "theta":params[1], "psi":params[2]})
	temp = volft.extract_plane(R,kb)
	temp.fft_shuffle()
	temp.center_origin_fft()

	if(params[3]!=0. or params[4]!=0.):
		filt_params = {"filter_type" : Processor.fourier_filter_types.SHIFT,
				  "x_shift" : params[3], "y_shift" : params[4], "z_shift" : 0.0}
		temp=Processor.EMFourierFilter(temp, filt_params)
	temp.do_ift_inplace()
	set_params_proj(temp, [params[0], params[1], params[2], -params[3], -params[4]])
	temp.set_attr_dict({'active':1, 'ctf_applied':0, 'npad':2})
	temp.depad()
	return temp

def prgs1d( prjft, kb, params ):
	from fundamentals import fft
	from math import cos, sin, pi

	alpha = params[0]
	shift = params[1]

	tmp = alpha/180.0*pi

	nuxnew =  cos(tmp)
	nuynew = -sin(tmp)
	
	line = prjft.extractline(kb, nuxnew, nuynew)
	line = fft(line)

	M = line.get_xsize()/2
	Util.cyclicshift( line, {"dx":M, "dy":0, "dz":0} )
	line = Util.window( line, M, 1, 1, 0, 0, 0 )

	if shift!=0:
		filt_params = {"filter_type" : Processor.fourier_filter_types.SHIFT,
		   	       "x_shift" : shift, "y_shift" : 0.0, "z_shift" : 0.0}
		line = Processor.EMFourierFilter(temp, filt_params)

	line.set_attr_dict( {'alpha':alpha, 's1x':shift} )
	return line

def prg(volume, params):
	"""Given a volume, a set of projection angles, and Kaiser-Bessel
	   window parameters, use gridding to generate projection
	"""
	volft,kb = prep_vol(volume)
	return  prgs(volft,kb,params)

def prep_vol(vol):
	# prepare the volume
	M     = vol.get_xsize()
	# padd two times
	npad  = 2
	N     = M*npad
	# support of the window
	K     = 6
	alpha = 1.75
	r     = M/2
	v     = K/2.0/N
	kb    = Util.KaiserBessel(alpha, K, r, K/(2.*N), N)
	volft = vol.copy()
	volft.divkbsinh(kb)
	volft = volft.norm_pad(False, npad)
	volft.do_fft_inplace()
	volft.center_origin_fft()
	volft.fft_shuffle()
	return  volft,kb

def next_type(typlst, max):
	nlen = len(typlst)
	typlst[nlen-1] += 1
	for i in xrange(nlen-1, 0, -1):
		if typlst[i] == max:
			typlst[i] = 0
			typlst[i-1] += 1
	return typlst[0] < max

def change_ang(ang, type):
	if type == 0:
		return ang

	if type == 1:
		r = 180.0 - ang
		if r < 0.0: r = r + 360.0
		return r

	if type == 2:
		r = ang - 180.0
		if r < 0.0: r = r + 360.0
		return r

	assert type==3
	return 360.0 - ang

def exhaust_ang(angs, anglst, func, data=None, info=None):
	from utilities import amoeba
	origin = angs[:]
	optpsi = angs[:]
	funval = func(angs, data)

	typlst = [0]*len(anglst)

	lowers = [0.0]*len(angs)
	uppers = [360.0]*len(angs)

	while( next_type(typlst, 4) ):
		work_angs = optpsi[:]
		for i in xrange( len(typlst) ):
			iang = anglst[i]
			ityp = typlst[i]
			work_angs[iang] = change_ang(origin[iang], ityp)
		#curt_angs = work_angs
		#curt_funval = func(work_angs, data)
	
		curt_angs,curt_funval = rollover(work_angs, lowers, uppers, func, data=data)

		#print 'typlst:', typlst
		info.write( '	     typlst,funval,best:' + str(typlst) + " " + str(curt_funval) + " " + str(funval) + "\n" )
		info.flush()

		if( curt_funval > funval ):
			optpsi = curt_angs[:]
			funval = curt_funval
	return optpsi,funval

def drop_angle_doc(filename, phis, thetas, psis, comment=None):
	from utilities import drop_spider_doc
	table = []
	for i in xrange( len(phis) ): table.append( [ phis[i], thetas[i], psis[i] ] )
	drop_spider_doc(filename, table, comment)


###############################################################################################
## COMMON LINES NEW VERSION ###################################################################
# 11/11/08 
# TODO: refinement, trials and MPI

# plot angles, map on half-sphere
# agls: [[phi0, theta0, psi0], [phi1, theta1, psi1], ..., [phin, thetan, psin]]
def plot_angles(agls):
	from math      import cos, sin, fmod, pi
	from utilities import model_blank

	# var
	nx = 256
	im = model_blank(nx, nx)
	c  = 2
	kc = 10

	# draw reperes
	for i in xrange(nx):
		im.set_value_at(i, int(nx / 2.0), 0.006)
		im.set_value_at(int(nx / 2.0), i, 0.006)
	
	# draw the circles
	lth = range(0, 90, kc)
	lth.append(90)
	
	for th in lth:
		
		if th == 90: color = 0.03
		else:        color = 0.006

		rc  = sin((float(th) / 180.0) * pi)
		rc *= (nx - 1)
		
		for n in xrange(3600):
			a  = (n / 1800.0) * pi
			px = nx / 2.0 + (rc - 1) / 2.0 * cos(a)
			py = nx / 2.0 + (rc - 1) / 2.0 * sin(a)
			im.set_value_at(int(px), int(py), color)
	
	# for each angles plot on circle area
	# agsl: [phi, theta, phi]
	for i in xrange(len(agls)):
		if agls[i][1] > 90:
			agls[i][0] = (float(agls[i][0]) + 180) % 360
			agls[i][1] = 180 - float(agls[i][1])
		
		rc  = sin((agls[i][1] / 180.0) * pi)
		rc *= ((nx - 1) / 2)

		px  = nx / 2.0 + rc * cos((fmod(agls[i][0], 360.0) / 180.0) * pi)
		py  = nx / 2.0 + rc * sin((fmod(agls[i][0], 360.0) / 180.0) * pi)

		if px > nx - 1: px = nx - 1
		if px < 0:  px = 0
		px = int(px)

		if py > nx - 1: py = nx - 1
		if py < 0:  py = 0
		py = int(py)

		#if agls[i][1] > 90: style = 2
		#else:               style = 1
		style = 1
	
		for cx in xrange(px - c, px + c + 1, style):
			if cx > nx - 1: cx = nx - 1
			if cx < 0:  cx = 0
			im.set_value_at(cx, py, 1.0)

		for cy in xrange(py - c, py + c + 1, style):
			if cy > nx - 1: cy = nx - 1
			if cy < 0:  cy = 0
			im.set_value_at(px, cy, 1.0)

	return im

# to utilities: get_line
def get_line(im, li):
	from utilities import model_blank
	nx = im.get_xsize()
	e  = model_blank(nx)
	for n in xrange(nx): e.set_value_at(n, 0, im.get_value_at(n, li))
	return e






# transform an image to sinogram (mirror include)
def cml_sinogram_(image2D, diameter = 53, d_psi = 1):
	from math         import cos, sin
	from fundamentals import fft
	from utilities    import model_blank
	M_PI  = 3.141592653589793238462643383279502884197
	
	diameter = int(diameter)
	ri = diameter // 2
	diameter = 2*ri + 1
	ri2 = ri*ri
	nx = image2D.get_xsize()
	ny = image2D.get_ysize()
	nxc = nx//2
	nyc = ny//2

	# get line projection
	nangle = int(360 / d_psi)     
	dangle = 2 * M_PI / float(nangle)
	e = model_blank(diameter,nangle)
	for j in xrange(nangle):
		cs =  cos(dangle * j)
		si = -sin(dangle * j)
		for iy in xrange(ny):
			oiy = iy - nyc
			for ix in xrange(nx):
				oix = ix - nxc
				if( (oiy*oiy + oix*oix) < ri2):
					xb   = oix*cs+oiy*si + ri
					ixb  = int(xb)
					dipx = xb - ixb
					val  = image2D.get_value_at(ix,iy)
					e.set_value_at(ixb, j, e.get_value_at(ixb, j) + (1.0-dipx)*val)
					e.set_value_at(ixb+1, j, e.get_value_at(ixb+1, j) + dipx*val)
				
	return Util.window(e, diameter-1, nangle, 1, 0, 0, 0)

# transform an image to sinogram (mirror include)
def cml_sinogram(image2D, diameter, d_psi = 1):
	from math         import cos, sin
	from fundamentals import fft
	
	M_PI  = 3.141592653589793238462643383279502884197
	
	# prepare 
	M = image2D.get_xsize()
	# padd two times
	npad  = 2
	N     = M * npad
	# support of the window
	K     = 6
	alpha = 1.75
	r     = M / 2
	v     = K / 2.0 / N

	kb     = Util.KaiserBessel(alpha, K, r, K / (2. * N), N)
	volft  = image2D.average_circ_sub()  	# ASTA - in spider
	volft.divkbsinh(kb)		  	# DIVKB2 - in spider
	volft  = volft.norm_pad(False, npad)
	volft.do_fft_inplace()
	volft.center_origin_fft()
	volft.fft_shuffle()

	# get line projection
	nangle = int(360 / d_psi)     
	dangle = 2 * M_PI / float(nangle)
	data   = []
	for j in xrange(nangle):
		nuxnew =  cos(dangle * j)
		nuynew = -sin(dangle * j)
		line   = volft.extractline(kb, nuxnew, nuynew)
		rlines = fft(line)
		data.append(rlines.copy())

	# copy each line in the same im
	e = EMData()
	e.set_size(data[0].get_xsize() ,len(data), 1)
	for n in xrange(len(data)):
		nx = data[n].get_xsize()
		for i in xrange(nx): e.set_value_at(i, n, data[n].get_value_at(i))

	Util.cyclicshift(e, {"dx":M, "dy":0, "dz":0} )

	return Util.window(e, diameter, len(data), 1, 0, 0, 0)

def common_line_in3D(ph1, th1, ph2, th2):
	from math import pi, sqrt, cos, sin, asin, acos

	deg_rad = pi / 180.0
	ph1 *= deg_rad 
	th1 *= deg_rad 
	ph2 *= deg_rad 
	th2 *= deg_rad

	# cross-product between normal vector of projections
	nx = sin(th1)*sin(ph1)*cos(th2) - cos(th1)*sin(th2)*sin(ph2)
	ny = cos(th1)*sin(th2)*cos(ph2) - cos(th2)*sin(th1)*cos(ph1)
	nz = sin(th1)*cos(ph1)*sin(th2)*sin(ph2) - sin(th1)*sin(ph1)*sin(th2)*cos(ph2)

	# normalize
	norm    = nx**2 + ny**2 + nz**2
	rt_norm = sqrt(norm)
	nx /= rt_norm
	ny /= rt_norm
	nz /= rt_norm

	# if theta > 90, apply mirror 
	if nz < 0: nx = -nx; ny = -ny; nz = -nz
	
	# calculate phi and theta (deg)
	thetaCom  = acos(nz)

	if    thetaCom == 0: phiCom = 0
	else:
		val = ny / sin(thetaCom)
		if val > 1.0:  val = 1.0
		if val < -1.0: val = -1.0
		phiCom = asin(val)
	
		phiCom    = (phiCom * 180 / pi + 360)%360
		thetaCom *= (180 / pi)

	return phiCom , thetaCom

def cml_weights_full(Ori):
	from projection   import common_line_in3D

	# gbl vars
	global g_n_prj, g_n_lines, g_anglst

	
	# gbl vars
	l_phs  = [0.0] * g_n_lines  # angle phi of the common lines
	l_ths  = [0.0] * g_n_lines  # angle theta of the common lines
	n      = 0
	for i in xrange(g_n_prj - 1):
		for j in xrange(i + 1, g_n_prj):
			l_phs[n], l_ths[n] = common_line_in3D(Ori[4*i], Ori[4*i+1], Ori[4*j], Ori[4*j+1])
			n+= 1

	tol = 3

	# search the closer cml lines
	ocp_same   = [-1] * g_n_lines
	num_agl    = 0
	for i in xrange(g_n_lines):
	    if ocp_same[i] == -1:
		ocp_same[i] = num_agl
		for j in xrange(i + 1, g_n_lines):
		    if ocp_same[j] == -1:
			dist = (l_phs[i] - l_phs[j])**2 + (l_ths[i] - l_ths[j])**2
			if dist < tol: ocp_same[j] = num_agl

		num_agl += 1

	if num_agl > 2:

		# create the new vector n_phi n_theta without closer
		n_phi   = [0.0] * num_agl
		n_theta = [0.0] * num_agl
		nb_same = [0]   * num_agl
		num_agl = 0
		for n in xrange(g_n_lines):
		    nb_same[ocp_same[n]] += 1
		    if ocp_same[n] == num_agl:
			n_phi[num_agl]   = l_phs[n]
			n_theta[num_agl] = l_ths[n]
			num_agl += 1

		# Voronoi
		n_weights = Util.vrdg(n_phi, n_theta)

		weights = [0.0] * g_n_lines
		for i in xrange(g_n_lines):
			if nb_same[ocp_same[i]] > 1:
				weights[i] = n_weights[ocp_same[i]] / float(nb_same[ocp_same[i]])
			else:
				weights[i] = n_weights[ocp_same[i]]

	else:
		weights = [6.28/float(g_n_lines)] * g_n_lines

	return weights
	

# compute the weight of the common lines
def cml_weights_iagl(Ori, iagl, iprj):
	from projection   import common_line_in3D

	# gbl vars
	global g_n_prj, g_n_lines, g_anglst

	
	# gbl vars
	l_phs  = [0.0] * g_n_lines  # angle phi of the common lines
	l_ths  = [0.0] * g_n_lines  # angle theta of the common lines
	n      = 0
	for i in xrange(g_n_prj - 1):
		for j in xrange(i + 1, g_n_prj):
			if i == iprj:   l_phs[n], l_ths[n] = common_line_in3D(g_anglst[iagl][0], g_anglst[iagl][1], Ori[4*j], Ori[4*j+1])
			elif j == iprj:	l_phs[n], l_ths[n] = common_line_in3D(Ori[4*i], Ori[4*i+1], g_anglst[iagl][0], g_anglst[iagl][1])
			else:		l_phs[n], l_ths[n] = common_line_in3D(Ori[4*i], Ori[4*i+1], Ori[4*j], Ori[4*j+1])
			n+= 1

	tol = 3

	# search the closer cml lines
	ocp_same   = [-1] * g_n_lines
	num_agl    = 0
	for i in xrange(g_n_lines):
	    if ocp_same[i] == -1:
		ocp_same[i] = num_agl
		for j in xrange(i + 1, g_n_lines):
		    if ocp_same[j] == -1:
			dist = (l_phs[i] - l_phs[j])**2 + (l_ths[i] - l_ths[j])**2
			#print i, j, dist
			if dist < tol: ocp_same[j] = num_agl

		num_agl += 1

	if num_agl > 2:

		# create the new vector n_phi n_theta without closer
		n_phi   = [0.0] * num_agl
		n_theta = [0.0] * num_agl
		nb_same = [0]   * num_agl
		num_agl = 0
		for n in xrange(g_n_lines):
		    nb_same[ocp_same[n]] += 1
		    if ocp_same[n] == num_agl:
			n_phi[num_agl]   = l_phs[n]
			n_theta[num_agl] = l_ths[n]
			num_agl += 1

		# Voronoi
		n_weights = Util.vrdg(n_phi, n_theta)

		weights = [0.0] * g_n_lines
		for i in xrange(g_n_lines):
			if nb_same[ocp_same[i]] > 1:
				weights[i] = n_weights[ocp_same[i]] / float(nb_same[ocp_same[i]])
			else:
				weights[i] = n_weights[ocp_same[i]]

	else:
		weights = [6.28/float(g_n_lines)] * g_n_lines

	return weights

# open and transform projections
def cml_open_proj(stack, ir, ou, d_psi, lf, hf):
	from projection   import cml_sinogram
	from utilities    import model_circle, get_params_proj, model_blank, get_im
	from fundamentals import fft

	nprj = EMUtil.get_image_count(stack)                # number of projections
	Prj = []                                            # list of projections
	Ori = [-1] * 4 * nprj                              # orientation intial (phi, theta, psi, index) for each projection

	for i in xrange(nprj):
		image = get_im(stack, i)

		# read initial angles if given
		try:	Ori[4*i], Ori[4*i+1], Ori[4*i+2], s2x, s2y = get_params_proj(image)
		except:	pass
		
		if(i == 0):
			nx = image.get_xsize()
			if(ou < 1):
				ou = nx // 2 - 1
			#else:
			#	ou = int(ou) // 2
			#	ou = 2 * ou +1
			diameter = 2 * ou + 1
			diameter = int(diameter)
			mask2D   = model_circle(ou, nx, nx)
			circ     = mask2D.copy()
			if ou > 1:  circ   -= model_circle(ou - 1, nx, nx)
			if ir > 0:  mask2D -= model_circle(ir, nx, nx)

		# normalize under the mask
		[mean_a, sigma, imin, imax] = Util.infomask(image, circ, True)
		image = (image - mean_a) / sigma
		Util.mul_img(image, mask2D)

		# sinogram
		sino = cml_sinogram(image, diameter, d_psi)

		# prepare the cut positions in order to filter (lf: low freq; hf: high freq)
		ihf = min(int(2 * hf * diameter), diameter + (diameter + 1) % 2)
		ihf = ihf + (ihf + 1) % 2    # index ihf must be odd to take the img part
		ilf = max(int(2 * lf * diameter), 0)
		ilf = ilf + ilf % 2          # index ilf must be even to fall in the real part
		bdf = ihf - ilf + 1

		# process lines
		nxe = sino.get_xsize()
		nye = sino.get_ysize()
		prj = model_blank(bdf, nye)
		prj.set_complex(True)
		for li in xrange(nye):

			# get the line li
			line = model_blank(nxe)
			for ci in xrange(nxe): line.set_value_at(ci, 0, sino.get_value_at(ci, li))

			# normalize this line
			[mean_l, sigma_l, imin, imax] = Util.infomask(line, None, True)
			line = (line - mean_l) / sigma_l

			# fft
			line = fft(line)
	
			# filter (cut part of coef)
			ct = 0
			for ci in xrange(ilf, ihf + 1):
				prj.set_value_at(ct, li, line.get_value_at(ci, 0))
				ct += 1
	
		# store the projection
		Prj.append(prj)

	return Prj, Ori

# export result obtain by the function find_struct
def cml_export_struc(stack, outseed, Ori, BDB):
	from projection import plot_angles
	from utilities  import set_params_proj, get_im

	global g_n_prj
	
	pagls = []
	for i in xrange(g_n_prj):
		data = get_im(stack, i)
		p = [Ori[4*i], Ori[4*i+1], Ori[4*i+2], 0.0, 0.0]
		set_params_proj(data, p)
		data.set_attr('active', 1)
		if BDB:	data.write_image(stack, i)
		else:	data.write_image(stack, i)

		# prepare angles to plot
		pagls.append([Ori[4*i], Ori[4*i+1], Ori[4*i+2]])

	# plot angles
	im = plot_angles(pagls)
	if BDB: im.write_image('bdb:%s_plot_agls' % outseed, 0)
	else:   im.write_image(outseed + 'plot_agls.hdf')

# init the global average used for lot of function to cml
def cml_init_global_var(dpsi, delta, nprj, debug):
	from utilities import even_angles
	global g_anglst, g_d_psi, g_n_psi, g_i_prj, g_n_lines, g_n_prj, g_n_anglst, g_debug
	
	g_anglst   = even_angles(delta, 0.0, 179.9, 0.0, 359.9, 'P')
	g_n_anglst = len(g_anglst)
	g_d_psi    = dpsi
	g_n_psi    = int(360 / dpsi)
	g_i_prj    = -1
	g_n_lines  = ((nprj - 1) * nprj) / 2
	g_n_prj    = nprj
	g_debug    = debug
	

# write the head of the logfile
def cml_head_log(stack, outdir, delta, ir, ou, lf, hf, rand_seed, maxit, given):
	from utilities import print_msg

	# call global var
	global g_anglst, g_n_prj, g_d_psi, g_n_anglst
	
	print_msg('Input stack                  : %s\n'     % stack)
	print_msg('Number of projections        : %d\n'     % g_n_prj)
	print_msg('Output directory             : %s\n'     % outdir)
	print_msg('Angular step                 : %5.2f\n'  % delta)
	print_msg('Sinogram angle accuracy      : %5.2f\n'  % g_d_psi)
	print_msg('Inner particle radius        : %5.2f\n'  % ir)	
	print_msg('Outer particle radius        : %5.2f\n'  % ou)
	print_msg('Filter, minimum frequency    : %5.3f\n'  % lf)
	print_msg('Filter, maximum frequency    : %5.3f\n'  % hf)
	print_msg('Random seed                  : %i\n'     % rand_seed)
	print_msg('Number of maximum iterations : %d\n'     % maxit)
	print_msg('Start from given orientations: %s\n'     % given)
	
	
	#print_msg('Number of trials            : %d\n'     % trials)
	#print_msg('Number of cpus              : %i\n'     % ncpu)
	#if refine:
	#	print_msg('Refinement                  : True\n')
	#else:
	#	print_msg('Refinement                  : False\n')
	print_msg('Number of angles             : %i\n\n'   % g_n_anglst)

# write the end of the logfile
def cml_end_log(Ori, disc, disc_nw, ite):
	from utilities import print_msg
	global g_n_prj
	print_msg('\n\n')
	for i in xrange(g_n_prj): print_msg('Projection #%s: phi %10.5f    theta %10.5f    psi %10.5f\n' % (str(i).rjust(3, '0'), Ori[4*i], Ori[4*i+1], Ori[4*i+2]))
	print_msg('\nNumber of iterations: %d\n' % ite)
	print_msg('Discrepancy: %10.3f\n' % abs(disc))
	print_msg('Discrepancy without weigths: %10.3f\n' % abs(disc_nw))

# display the list of angles for each iterations
def cml_export_txtagls(outdir, Ori, disc, title):
	import time
	global g_n_prj, g_i_prj

	angfile = open(outdir + 'angles', 'a')

	angfile.write('|%s|-----------------------------------------------%s---------\n' % (title, time.ctime()))
	for i in xrange(g_n_prj): angfile.write('%10.3f\t%10.3f\t%10.3f\n' % (Ori[4*i], Ori[4*i+1], Ori[4*i+2]))
			
	angfile.write('\nDiscrepancy: %10.3f\n\n' % disc)
	angfile.close()

# export the progress of the find_struc function
def cml_export_progress(outdir, ite, iprj, iagl, psi, disc, cmd):
	infofile = open(outdir + 'progress', 'a')
	global g_anglst

	if cmd == 'progress':
		txt_ite = str(ite).rjust(3, '0')
		txt_i   = str(iprj).rjust(3, '0')
		txt_a   = str(iagl).rjust(3, '0')
		txt     = 'Ite: %s Prj: %s Agls: %s >> Agls (phi, theta, psi): %10.3f %10.3f %10.3f   Disc: %10.7f' % (txt_ite, txt_i, txt_a, g_anglst[iagl][0], g_anglst[iagl][1], psi, disc)
		

	elif cmd == 'choose':
		txt   = 'Ite: %s  Select Agls: %s >> Agls (phi, theta, psi): %10.3f %10.3f %10.3f   Disc: %10.7f\n' % (str(ite).rjust(3, '0'), str(iagl).rjust(3, '0'), g_anglst[iagl][0], g_anglst[iagl][1], psi, disc)

	infofile.write(txt + '\n')
	infofile.close()

# compute the common lines in sino
def get_common_line_angles(phi1, theta1, psi1, phi2, theta2, psi2, nangle, STOP=False):
	from math import fmod
	R1    = Transform({"type":"spider", "phi":phi1, "theta":theta1, "psi":psi1})
	R2    = Transform({"type":"spider", "phi":phi2, "theta":theta2, "psi":psi2})
	R2T   = R2.inverse()
	R2to1 = R1*R2T

	eulerR2to1 = R2to1.get_rotation("spider")
	phiR2to1   = eulerR2to1["phi"]
	thetaR2to1 = eulerR2to1["theta"]
	psiR2to1   = eulerR2to1["psi"]

	alphain1 = fmod(psiR2to1  + 270.0, 360.0)
	alphain2 = fmod(-phiR2to1 + 270.0, 360.0)

	n1 = int(nangle * fmod(alphain1 + 360, 360) / 360.0)
	n2 = int(nangle * fmod(alphain2 + 360, 360) / 360.0)
	
	return n1, n2

# compute discrepancy according the projections and orientations
def cml_disc(Prj, Ori, flag_weights):
	from projection  import get_common_line_angles, cml_weights_full
	from math        import pi, fmod

	# gbl vars
	global g_n_prj, g_n_psi, g_n_lines

	if flag_weights:
		cml = Util.cml_line_in3d_full(Ori)    # c-code
		weights = Util.cml_weights(cml)       # c-code
		#weights = cml_weights_full_dev(Ori)  # py-code
		
	else:   weights = [1.0] * g_n_lines

	#com = [[] for i in xrange(g_n_lines)]
	com = [0] * 2 * g_n_lines

	# compute the common lines
	count = 0
	for i in xrange(g_n_prj - 1):
		for j in xrange(i + 1, g_n_prj):
			#com[count], com[count + 1] = get_common_line_angles(Ori[4*i], Ori[4*i+1], Ori[4*i+2], Ori[4*j], Ori[4*j+1], Ori[4*j+2], g_n_psi)    # py code
			[com[count], com[count + 1]] = Util.cml_line_pos(Ori[4*i], Ori[4*i+1], Ori[4*i+2], Ori[4*j], Ori[4*j+1], Ori[4*j+2], g_n_psi)        # c  code
			count += 2

	n = 0
	L_tot = 0.0

	# compute the discrepancy for all sinograms
	for i in xrange(g_n_prj - 1):
		for j in xrange(i + 1, g_n_prj):
			L      = Prj[i].cm_euc(Prj[j], com[n], com[n + 1])
			L_tot += (L * weights[int(n/2)])
			n     += 2

	return L_tot

# interface between the simplex function to refine the angles and the function to compute the discrepancy
def cml_refine_agls_wrap(vec_in, data):
	# vec_in: [phi_i, theta_i, psi_i]
	# data:   [Prj, Ori, iprj]

	# unpack
	phi, theta, psi = vec_in
	Prj, Ori, iprj  = data

	# prepare the variables
	Ori[4*iprj]   = phi
	Ori[4*iprj+1] = theta
	Ori[4*iprj+2] = psi

	# compute the discrepancy
	disc = cml_disc(Prj, Ori, True)

	return -disc

# cml refines angles
def cml_refine_agls(Prj, Ori, delta):
	from copy      import deepcopy
	from utilities import amoeba
	global g_n_prj
	
	scales = [delta] * (g_n_prj + 2)

	for iprj in xrange(g_n_prj):
		# init vec_in
		vec_in   = [Ori[4*iprj], Ori[4*iprj+1], Ori[4*iprj+2]]

		# prepare vec_data
		vec_data = [Prj, deepcopy(Ori), iprj]

		# simplex
		optvec, disc, niter = amoeba(vec_in, scales, cml_refine_agls_wrap_dev, data = vec_data)

		# assign new angles refine
		Ori[4*iprj]   = (optvec[0]+360)%360
		Ori[4*iprj+1] = optvec[1]
		Ori[4*iprj+2] = optvec[2]

		print 'refine:', iprj, 'angles:', Ori[4*iprj:4*iprj+4], 'disc:', -disc

	return Ori

# cml spin function for one orientation
def cml_spin(Prj, iprj, Ori, iagl, weights, flag):
	# gbl vars
	global g_n_prj, g_n_psi, g_n_lines, g_anglst
	com = [0] * 2 * g_n_lines

	# compute the common line (only for iprj)
	com = Util.cml_list_line_pos(Ori, g_anglst[iagl][0], g_anglst[iagl][1], iprj, g_n_prj, g_n_psi, g_n_lines)
        #print com
	#import sys
	#sys.exit()
	'''
	if g_anglst[iagl][0] > 179 and g_anglst[iagl][0] < 200 and g_anglst[iagl][1] == 80:
		print iprj, iagl, g_anglst[iagl][0], g_anglst[iagl][1]
		#print Ori
		for i in xrange(0, 2*g_n_lines, 2):
			print com[i], com[i+1]
	'''

	#import sys
	#sys.exit()
	# 337 100
	# 251 40
	#if g_anglst[iagl][0] > 251 and g_anglst[iagl][0] < 252 and g_anglst[iagl][1] == 40 and iprj == 0:
	#	print '====', iprj, iagl, g_anglst[iagl][0], g_anglst[iagl][1]
	
	# do spin over all psi
	#if iprj == 01 and iagl == 256:
	#	print 'iprj 01 iagl 256'
	res = Util.cml_spin(g_n_psi, iprj, g_n_prj, weights, com, Prj, flag)
	#else:
	#	res = Util.cml_spin(g_n_psi, iprj, g_n_prj, weights, com, Prj, 0)
	#else:
	#	res = [1e9, 0]

	#return best_disc, best_psi
	return res[0], int(res[1])

def cml_norm_weights(w):
	wm = 6.28 #max(w)
	nw = []
	for i in xrange(len(w)): nw.append(wm - w[i])
	sw = sum(nw)
	for i in xrange(len(w)): nw[i] /= sw
	return nw

# find structure
def cml_find_structure(Prj, Ori, outdir, maxit, first_zero, flag_weights):
	from projection import cml_spin, cml_export_progress, cml_weights_iagl
	import time
	import sys
	
	# global vars
	global g_i_prj, g_n_prj, g_n_anglst, g_anglst, g_d_psi, g_debug, g_n_lines

	# list of free orientation
	ocp = [-1] * g_n_anglst

	if first_zero:
		listprj = range(1, g_n_prj)
		ocp[0]  = 0 
	else:   listprj = range(g_n_prj)

	# iteration loop
	for ite in xrange(maxit):
		t_start = time.time()

		# loop over i prj
		change = False
		for iprj in listprj:
			

			# Store current index of angles assign
			cur_agl      = Ori[4*iprj+3]
			ocp[cur_agl] = -1

			# loop over all angles
			best_disc = 1e20
			best_psi  = -1
			best_iagl = -1
			for iagl in xrange(g_n_anglst):
				#if g_debug: print 'iprj:', iprj, 'iagl', iagl

				# if agls free
				if ocp[iagl] == -1:
					# weights
					if flag_weights:
						cml = Util.cml_line_in3d_iagl(Ori, g_anglst[iagl][0], g_anglst[iagl][1], iprj)   # c-code
						
						#if g_debug: print ite, iprj, iagl
						'''
						if iprj == 11 and iagl == 177 and ite == 1:
							for i in xrange(g_n_prj):
								print '%10.3f\t%10.3f\t%10.3f' % (Ori[4*i], Ori[4*i+1], Ori[4*i+2])
							print 'projection:', iprj, g_anglst[iagl]
							print '-------------------------'
							for i in xrange(0, len(cml), 2):
								print cml[i], cml[i+1]
							sys.exit()
						'''
						weights = Util.cml_weights(cml)                                                  # c-code
						weights = cml_norm_weights(weights)
						
                                                #if sum(weights) < 6.28:
						#	print 'WARNING ite', ite, 'iprj', iprj, 'iagl', iagl 
                                                #weights = cml_weights_iagl(Ori, iagl, iprj)                                     # py-code
						
					else:   weights = [1.0] * g_n_lines
					
					# spin
					#if g_anglst[iagl][1] == 100:
					#print g_anglst[iagl][0], g_anglst[iagl][1]
					#if iagl == 0 and iprj == 1 and ite == 0:
					#	disc, ind_psi = cml_spin(Prj, iprj, Ori, iagl, weights, 1)
					#	sys.exit()
					#else:
					disc, ind_psi = cml_spin(Prj, iprj, Ori, iagl, weights, 0)
					#else:
					#	disc, ind_psi = 1e9, 0

					# select the best
					if disc < best_disc:
						best_disc = disc
						best_psi  = ind_psi
						best_iagl = iagl

					if g_debug: cml_export_progress(outdir, ite, iprj, iagl, ind_psi * g_d_psi, disc, 'progress')
				else:
					if g_debug: cml_export_progress(outdir, ite, iprj, iagl, -1, -1, 'progress')


			# if change, assign
			if best_iagl != cur_agl:
				ocp[best_iagl] = iprj
				Ori[4*iprj]    = g_anglst[best_iagl][0] # phi
				Ori[4*iprj+1]  = g_anglst[best_iagl][1] # theta
				Ori[4*iprj+2]  = best_psi * g_d_psi     # psi
				Ori[4*iprj+3]  = best_iagl              # index

				change = True
			else:
				ocp[cur_agl]   = iprj

			if g_debug: cml_export_progress(outdir, ite, iprj, best_iagl, best_psi * g_d_psi, best_disc, 'choose')

		# if one change, compute new full disc
		disc = cml_disc(Prj, Ori, flag_weights)

		if g_debug:
			disc2 = cml_disc(Prj, Ori, False)
			print 'Ite: ', disc, '           %6.2f s     ' % (time.time() - t_start), disc2

		# display in the progress file
		cml_export_txtagls(outdir, Ori, disc, 'Ite: %s' % str(ite + 1).rjust(3, '0'))

		if not change: break

	#sys.exit()

	return Ori, disc, ite
	
# cml init for MPI version
def cml_init_MPI(trials):
	from mpi 	  import mpi_init, mpi_comm_size, mpi_comm_rank, mpi_barrier, MPI_COMM_WORLD
	from utilities    import bcast_number_to_all
	from random       import randint
	import sys
	
	# init
	sys.argv       = mpi_init(len(sys.argv),sys.argv)
	number_of_proc = mpi_comm_size(MPI_COMM_WORLD)
	myid           = mpi_comm_rank(MPI_COMM_WORLD)

	# chose a random node as a main one
	main_node = 0
	if myid  == 0:	main_node = randint(0, number_of_proc - 1)
	main_node = bcast_number_to_all(main_node, 0)
	mpi_barrier(MPI_COMM_WORLD)

	# define the number of loop per node
	loop = trials / number_of_proc
	if trials % number_of_proc != 0: loop += 1

	return main_node, myid, number_of_proc, loop

# cml init list of rand_seed for trials version
def cml_init_rnd(trials, rand_seed):
	from random import seed, randrange

	if trials == 1: return [rand_seed]
	
	if rand_seed > 0: seed(rand_seed)
	else:             seed()

	r_min = 100
	r_max = 1000000
	f_min = 1
	f_max = 100

	rnd     = []
	itrials = 0
	while itrials < trials:
		val_rnd = randrange(r_min, r_max)
		val_f   = randrange(f_min, f_max)
		val_o   = randrange(0, 2)
		if val_o: val_rnd = int(val_rnd * val_f)
		else:     val_rnd = int(val_rnd / float(val_f))
		if val_rnd not in rnd:
			rnd.append(val_rnd)
			itrials += 1

	return rnd

# ATTENTION: works well only without Voronoi weights
def find_structure_MPI(stack, trials, ir, ou, delta, dpsi, lf, hf, rand_seed, maxit):
	from projection  import find_struct, cml_init_rnd, cml_init_MPI
	from mpi         import mpi_barrier, mpi_reduce, mpi_bcast
	from mpi         import MPI_COMM_WORLD, MPI_FLOAT, MPI_INT, MPI_SUM
	import os

	main_node, myid, ncpu, loop = cml_init_MPI(trials)
	Rnd = cml_init_rnd(loop * ncpu, rand_seed)
	print Rnd

	for i in xrange(loop):
		outdir = 'trials_%s' % str(myid * loop + i).rjust(4, '0')
		ret = 0
		try:
			ret = find_struct(stack, outdir, ir, ou, delta, dpsi, lf, hf, Rnd[myid * loop + i], maxit, False, False, True, False)
		except SystemExit:
			ret = 0

		f = open('report_node%s' % str(myid).rjust(2, '0'), 'a')
		if ret: txt = 'passed'
		else:   txt = 'failed'
		f.write('node %s   trials %s: %s\n' % (str(myid).rjust(2, '0'), str(myid * loop + i).rjust(4, '0'), txt))
		f.close()


## END COMMON LINES NEW VERSION ###############################################################
###############################################################################################

