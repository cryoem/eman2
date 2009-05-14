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

# cml init for MPI version
def cml_init_MPI(trials):
	from mpi 	  import mpi_init, mpi_comm_size, mpi_comm_rank, mpi_barrier, MPI_COMM_WORLD
	from utilities    import bcast_number_to_all
	from random       import randint
	import sys
	
	# init
	sys.argv  = mpi_init(len(sys.argv),sys.argv)
	ncpu      = mpi_comm_size(MPI_COMM_WORLD)
	myid      = mpi_comm_rank(MPI_COMM_WORLD)
	main_node = 0
	mpi_barrier(MPI_COMM_WORLD)

	N_start = int(round(float(trials) / ncpu * myid))
	N_stop  = int(round(float(trials) / ncpu * (myid + 1)))

	return main_node, myid, ncpu, N_start, N_stop

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

# calculate the discrepancy allong all common-lines 
def cml_disc(Prj, Ori, Rot, flag_weights):
	from math        import pi, fmod

	# gbl vars
	global g_n_prj, g_n_psi, g_n_lines, g_seq

	if flag_weights:
		cml = Util.cml_line_in3d(Ori, g_seq, g_n_prj, g_n_lines)
		weights = Util.cml_weights(cml)       # c-code
		wm = max(weights)
		for i in xrange(g_n_lines): weights[i]  = wm - weights[i]
		sw = sum(weights)
		for i in xrange(g_n_lines): weights[i] /= sw
				
	else:   weights = [1.0] * g_n_lines

	com  = Util.cml_line_insino_all(Rot, g_seq, g_n_prj, g_n_lines)
	disc = Util.cml_disc(Prj, com, g_seq, weights, g_n_lines)
	
	return disc

# export the progress of the find_struc function
def cml_export_progress(outdir, ite, iprj, iagl, psi, mir, disc, cmd):
	infofile = open(outdir + '/progress', 'a')
	global g_anglst

	if cmd == 'progress':
		txt_ite = str(ite).rjust(3, '0')
		txt_i   = str(iprj).rjust(3, '0')
		txt_a   = str(iagl).rjust(3, '0')
		if mir:
			txt     = 'Ite: %s Prj: %s Agls: %s >> Agls (phi, theta, psi): %10.3f %10.3f %10.3f * Disc: %10.7f' % (txt_ite, txt_i, txt_a, g_anglst[iagl][0], g_anglst[iagl][1], psi, disc)
		else:
			txt     = 'Ite: %s Prj: %s Agls: %s >> Agls (phi, theta, psi): %10.3f %10.3f %10.3f   Disc: %10.7f' % (txt_ite, txt_i, txt_a, g_anglst[iagl][0], g_anglst[iagl][1], psi, disc)

	elif cmd == 'choose':
		if mir:
			txt   = 'Ite: %s  Select Agls: %s >> Agls (phi, theta, psi): %10.3f %10.3f %10.3f * Disc: %10.7f\n' % (str(ite).rjust(3, '0'), str(iagl).rjust(3, '0'), g_anglst[iagl][0], g_anglst[iagl][1], psi, disc)
		else:
			txt   = 'Ite: %s  Select Agls: %s >> Agls (phi, theta, psi): %10.3f %10.3f %10.3f   Disc: %10.7f\n' % (str(ite).rjust(3, '0'), str(iagl).rjust(3, '0'), g_anglst[iagl][0], g_anglst[iagl][1], psi, disc)

	infofile.write(txt + '\n')
	infofile.close()


# display the list of angles for each iterations
def cml_export_txtagls(outdir, outname, Ori, disc, title):
	import time
	global g_n_prj, g_i_prj

	angfile = open(outdir + '/' + outname, 'a')

	angfile.write('|%s|-----------------------------------------------%s---------\n' % (title, time.ctime()))
	for i in xrange(g_n_prj): angfile.write('%10.3f\t%10.3f\t%10.3f\n' % (Ori[4*i], Ori[4*i+1], Ori[4*i+2]))
			
	angfile.write('\nDiscrepancy: %10.3f\n\n' % disc)
	angfile.close()

# init global variables used toa quick acces with many function of common-lines
def cml_init_global_var(dpsi, delta, nprj, debug):
	from utilities import even_angles

	global g_anglst, g_d_psi, g_n_psi, g_i_prj, g_n_lines, g_n_prj, g_n_anglst, g_debug, g_seq
	
	g_anglst   = even_angles(delta, 0.0, 179.9, 0.0, 359.9, 'P')
	g_n_anglst = len(g_anglst)
	g_d_psi    = dpsi
	g_n_psi    = int(360 / dpsi)
	g_i_prj    = -1
	g_n_lines  = (nprj - 1) * nprj / 2
	g_n_prj    = nprj
	g_debug    = debug
	g_seq      = [0] * 2 * g_n_lines
	c          = 0
	# prepare pairwise indexes ij
	for i in xrange(g_n_prj):
		for j in xrange(i+1, g_n_prj):
			g_seq[c]   = i
			g_seq[c+1] = j
			c += 2

# export result obtain by the function find_struct
def cml_export_struc(stack, outdir, irun, Ori):
	from projection import plot_angles
	from utilities  import set_params_proj, get_im

	global g_n_prj
	
	pagls = []
	for i in xrange(g_n_prj):
		data = get_im(stack, i)
		p = [Ori[4*i], Ori[4*i+1], Ori[4*i+2], 0.0, 0.0]
		set_params_proj(data, p)
		data.set_attr('active', 1)
		data.write_image(outdir + '/structure_%03i.hdf' % irun, i)

		# prepare angles to plot
		pagls.append([Ori[4*i], Ori[4*i+1], Ori[4*i+2]])

	# plot angles
	im = plot_angles(pagls)
	im.write_image(outdir + '/plot_agls_%03i.hdf' % irun)

# open and transform projections to sinogram
def cml_open_proj(stack, ir, ou, lf, hf, dpsi = 1):
	from projection   import cml_sinogram
	from utilities    import model_circle, get_params_proj, model_blank, get_im
	from fundamentals import fft
	from filter       import filt_tanh

	nprj = EMUtil.get_image_count(stack)               # number of projections
	Prj  = []                                          # list of projections
	Ori  = [-1] * 4 * nprj                             # orientation intial (phi, theta, psi, index) for each projection

	for i in xrange(nprj):
		image = get_im(stack, i)

		# read initial angles if given
		try:	Ori[4*i], Ori[4*i+1], Ori[4*i+2], s2x, s2y = get_params_proj(image)
		except:	pass
		
		if(i == 0):
			nx = image.get_xsize()
			if(ou < 1): ou = nx // 2 - 1
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
		sino = cml_sinogram(image, diameter, dpsi) 

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

			# u2 (not improve the results)
			#line = filt_tanh(line, ou / float(nx), ou / float(nx))

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

# write the head of the logfile
def cml_head_log(stack, outdir, delta, ir, ou, lf, hf, rand_seed, maxit, given, flag_weights, trials):
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
	print_msg('Number of angles             : %i\n'     % g_n_anglst)
	print_msg('Number of trials             : %i\n'     % trials)
	print_msg('Use Voronoi weights          : %s\n\n'   % flag_weights)

# write the end of the logfile
def cml_end_log(Ori, disc, ite):
	from utilities import print_msg
	global g_n_prj
	print_msg('\n\n')
	for i in xrange(g_n_prj): print_msg('Projection #%03i: phi %10.5f    theta %10.5f    psi %10.5f\n' % (i, Ori[4*i], Ori[4*i+1], Ori[4*i+2]))
	print_msg('\nNumber of iterations: %d\n' % ite)
	print_msg('Discrepancy: %10.3f\n' % abs(disc))

# find structure
def cml_find_structure(Prj, Ori, Rot, outdir, outname, maxit, first_zero, flag_weights):
	from projection import cml_export_progress, cml_disc, cml_export_txtagls
	import time, sys
	
	# global vars
	global g_i_prj, g_n_prj, g_n_anglst, g_anglst, g_d_psi, g_debug, g_n_lines, g_seq

	# list of free orientation
	ocp = [-1] * g_n_anglst

	if first_zero:
		listprj = range(1, g_n_prj)
		ocp[0]  = 0 
	else:   listprj = range(g_n_prj)
	
	# to stop when the solution oscillate
	period_disc = [0, 0, 0]
	period_ct   = 0
	period_th   = 2

	# iteration loop
	for ite in xrange(maxit):
		t_start = time.time()

		# loop over i prj
		change = False
		for iprj in listprj:

			# Store current the current orientation
			ind          = 4*iprj
			store_phi    = Ori[ind]
			store_theta  = Ori[ind+1]
			store_psi    = Ori[ind+2]
			cur_agl      = Ori[ind+3]
			if cur_agl  != -1: ocp[cur_agl] = -1

			# TODO optimize that
			iw = [0] * (g_n_prj - 1)
			c  = 0
			ct = 0
			for i in xrange(g_n_prj):
				for j in xrange(i+1, g_n_prj):
					if i == iprj or j == iprj:
						iw[ct] = c
						ct += 1
					c += 1
			
			# loop over all angles
			best_disc = 1e20
			best_psi  = -1
			best_iagl = -1
			for iagl in xrange(g_n_anglst):
				# if orientation is free
				if ocp[iagl] == -1:
					# assign new orientation
					Ori[ind]   = g_anglst[iagl][0]
					Ori[ind+1] = g_anglst[iagl][1]
					Rot        = Util.cml_update_rot(Rot, iprj, Ori[ind], Ori[ind+1], 0.0)
					# weights
					if flag_weights:
						cml = Util.cml_line_in3d(Ori, g_seq, g_n_prj, g_n_lines)
						weights = Util.cml_weights(cml)
						minw = min(weights)
						sc   = max(weights) - minw
						if sc == 0.0: sc = 1.0
						for i in xrange(g_n_lines):
							weights[i]  = (weights[i] - minw) / sc
							weights[i]  = 1 - weights[i]
							weights[i] *= weights[i]
						# TODO optimize that
						#wm = max(weights)
						#for i in xrange(g_n_lines): weights[i]  = wm - weights[i]
						#sw = sum(weights)
						#if sw == 0:
						#	weights = [1.0] * g_n_lines
						#else:
						#	for i in xrange(g_n_lines):
						#		weights[i] /= sw
						#		weights[i] *= weights[i]
					else:   weights = [1.0] * g_n_lines

					# spin all psi
					com = Util.cml_line_insino(Rot, iprj, g_n_prj)
					res = Util.cml_spin_psi(Prj, com, weights, iprj, iw, g_n_psi, g_d_psi, g_n_prj)

					# select the best
					if res[0] < best_disc:
						best_disc = res[0]
						best_psi  = res[1]
						best_iagl = iagl
				
					if g_debug: cml_export_progress(outdir, ite, iprj, iagl, res[1], res[0], 'progress')
				else:
					if g_debug: cml_export_progress(outdir, ite, iprj, iagl, -1, -1, 'progress')

			# if change, assign
			if best_iagl != cur_agl:
				ocp[best_iagl] = iprj
				Ori[ind]       = g_anglst[best_iagl][0] # phi
				Ori[ind+1]     = g_anglst[best_iagl][1] # theta
				Ori[ind+2]     = best_psi * g_d_psi     # psi
				Ori[ind+3]     = best_iagl              # index
				change = True
			else:
				if cur_agl != -1: ocp[cur_agl] = iprj
				Ori[ind]    = store_phi
				Ori[ind+1]  = store_theta
				Ori[ind+2]  = store_psi
				Ori[ind+3]  = cur_agl

			Rot = Util.cml_update_rot(Rot, iprj, Ori[ind], Ori[ind+1], Ori[ind+2])
	
			if g_debug: cml_export_progress(outdir, ite, iprj, best_iagl, best_psi * g_d_psi, best_disc, 'choose')

		# if one change, compute new full disc
		disc = cml_disc(Prj, Ori, Rot, flag_weights)
		
		# display in the progress file
		cml_export_txtagls(outdir, outname, Ori, disc, 'Ite: %03i' % (ite + 1))

		if not change: break

		# to stop when the solution oscillate
		period_disc.pop(0)
		period_disc.append(disc)
		if period_disc[0] == period_disc[2]:
			period_ct += 1
			if period_ct >= period_th and min(period_disc) == disc:
				angfile = open(outdir + '/' + outname, 'a')
				angfile.write('\nSTOP SOLUTION UNSTABLE\n')
				angfile.write('Discrepancy periode: %s\n' % period_disc)
				angfile.close()
				break
		else:
			period_ct = 0

	'''
	cml = Util.cml_line_in3d(Ori, g_seq, g_n_prj, g_n_lines)
	w1 = Util.cml_weights(cml)
	wm = max(w1)
	w2 = [0.0] * g_n_lines
	for i in xrange(g_n_lines): w2[i]  = wm - w1[i]
	sw = sum(w2)
	if sw == 0:
		w2 = [1.0] * g_n_lines
	else:
		for i in xrange(g_n_lines):
			w2[i] /= sw
			w2[i] *= w2[i]

	f = open('weights.txt', 'w')
	for i in xrange(g_n_lines):
		f.write('%f\t%f\n' % (w1[i], w2[i]))
	f.close()

	la = []
	for i in xrange(g_n_lines // 2):
		ind = i*2
		la.append([cml[ind], cml[ind+1], 0.0])
	im = plot_angles(la)
	im.write_image('plot_cml.hdf', 0)
	'''

	return Ori, disc, ite

## END COMMON LINES NEW VERSION ###############################################################
###############################################################################################

