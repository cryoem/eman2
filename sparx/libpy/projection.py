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

def project(volume, params, radius=-1):
	"""
		Name
			project - calculate 2-D projection of a 3-D volume using trilinear interpolation
		Input
			vol: input volume, all dimensions have to be the same
			params: input parameters given as a list [phi, theta, psi, s2x, s2y], projection in calculated using the three Eulerian angles and then shifted by s2x,s2y
		radius: radius of a sphere within which the projection of the volume will be calculated
		Output
		proj: generated 2-D projection
	"""
        # angles phi, theta, psi
	from fundamentals import rot_shift2D
	from utilities import set_params_proj
	from EMAN2 import Processor

	if(radius>0):	myparams = {"transform":Transform({"type":"spider","phi":params[0],"theta":params[1],"psi":params[2]}), "radius":radius}
	else:			myparams = {"transform":Transform({"type":"spider","phi":params[0],"theta":params[1],"psi":params[2]})}
	proj = volume.project("pawel", myparams)
	if(params[3]!=0. or params[4]!=0.):
		params2 = {"filter_type" : Processor.fourier_filter_types.SHIFT, "x_shift" : params[3], "y_shift" : params[4], "z_shift" : 0.0}
		proj=Processor.EMFourierFilter(proj, params2)
		#proj = rot_shift2D(proj, sx = params[3], sy = params[4], interpolation_method = "linear")
	set_params_proj(proj, [params[0], params[1], params[2], -params[3], -params[4]])
	# horatio active_refactoring Jy51i1EwmLD4tWZ9_00000_1
	# proj.set_attr_dict({'active':1, 'ctf_applied':0})
	proj.set_attr_dict({ 'ctf_applied':0})
	return  proj

'''
Temporarily disabled as list cannot be passed to projector.
def prl(vol, params, radius, stack = None):
	"""
		Name
			prl - calculate a set of 2-D projection of a 3-D volume
		Input
			vol: input volume, all dimensions have to be the same (nx=ny=nz)
			params: a list of input parameters given as a list [i][phi, theta, psi, s2x, s2y], projection in calculated using the three Eulerian angles and then shifted by s2x,s2y
			radius: integer radius of a sphere - projection is calculated using voxels within this sphere.
		Output
			proj
				either: an in-core stack of generated 2-D projection
			stack
	"""
	from fundamentals import rot_shift2D
	for i in xrange(len(params)):
        	myparams = {"angletype":"SPIDER", "anglelist":params[i][0:3], "radius":radius}
        	proj = vol.project("pawel", myparams)
		if(params[i][3]!=0. or params[i][4]!=0.): proj = rot_shift2D(proj, sx = params[i][3], sy = params[i][4], interpolation_method = "linear")
		proj.set_attr_dict({'phi':params[i][0], 'theta':params[i][1], 'psi':params[i][2], 's2x':-params[i][3], 's2y':-params[i][4]})
		# horatio active_refactoring Jy51i1EwmLD4tWZ9_00000_1
		# proj.set_attr_dict({'active':1, 'ctf_applied':0})
		proj.set_attr_dict({ 'ctf_applied':0})
		
		if(stack):
			proj.write_image(stack, i)
		else:
			if(i == 0): out= []
			out.append(proj)
	if(stack): return
	else:      return out
'''
def prj(vol, params, stack = None):
	"""
		Name
			prj - calculate a set of 2-D projection of a 3-D volume
		Input
			vol: input volume, all dimensions have to be the same (nx=ny=nz)
			params: a list of input parameters given as a list [i][phi, theta, psi, sx, sy], projection in calculated using the three Eulerian angles and then shifted by sx,sy
		Output
			proj
				either: an in-core stack of generated 2-D projections
			stack
	"""
	from utilities  import set_params_proj
	from projection import prep_vol
	volft,kb = prep_vol(vol)
	for i in xrange(len(params)):
		proj = prgs(volft, kb, params[i])
		set_params_proj(proj, [params[i][0], params[i][1], params[i][2], -params[i][3], -params[i][4]])
		# horatio active_refactoring Jy51i1EwmLD4tWZ9_00000_1
		# proj.set_attr_dict({'active':1, 'ctf_applied':0})
		proj.set_attr_dict({ 'ctf_applied':0})
		
		if(stack):
			proj.write_image(stack, i)
		else:
			if(i == 0): out= []
			out.append(proj)
	if(stack):  return
	else:       return out

def prgs(volft, kb, params, kbx=None, kby=None):
	"""
		Name
			prg - calculate 2-D projection of a 3-D volume
		Input
			vol: input volume, the volume can be either cubic or rectangular
			kb: interpolants generated using prep_vol (tabulated Kaiser-Bessel function). If the volume is cubic, kb is the only interpolant.
			    Otherwise, kb is the for caculating weigthing along z direction.
			kbx,kby: interpolants generated using prep_vol used to calculae weighting aling x and y directin. Default is none when the volume is cubic. 
				 If the volume is rectangular, kbx and kby must be given.
			params: input parameters given as a list [phi, theta, psi, s2x, s2y], projection in calculated using the three Eulerian angles and then shifted by sx,sy
		Output
			proj: generated 2-D projection
	"""
	#  params:  phi, theta, psi, sx, sy
	from fundamentals import fft
	from utilities import set_params_proj
	from EMAN2 import Processor

	R = Transform({"type":"spider", "phi":params[0], "theta":params[1], "psi":params[2]})
	if kbx is None:
		temp = volft.extract_plane(R,kb)
	else:
		temp = volft.extract_plane_rect(R,kbx,kby,kb)
		
	temp.fft_shuffle()
	temp.center_origin_fft()

	if(params[3]!=0. or params[4]!=0.):
		filt_params = {"filter_type" : Processor.fourier_filter_types.SHIFT,
				  "x_shift" : params[3], "y_shift" : params[4], "z_shift" : 0.0}
		temp=Processor.EMFourierFilter(temp, filt_params)
	temp.do_ift_inplace()
	set_params_proj(temp, [params[0], params[1], params[2], -params[3], -params[4]])
	# horatio active_refactoring Jy51i1EwmLD4tWZ9_00000_1	
	# temp.set_attr_dict({'active':1, 'ctf_applied':0, 'npad':2})
	temp.set_attr_dict({'ctf_applied':0, 'npad':2})
	temp.depad()
	return temp
	

def gen_rings_ctf( prjref, nx, ctf, numr):
	"""
	  Convert set of ffts of projections to Fourier rings with additional multiplication by a ctf
	  The command returns list of rings
	"""
        from math         import sin, cos, pi
	from fundamentals import fft
	from alignment    import ringwe
	from filter       import filt_ctf
	mode = "F"
	wr_four  = ringwe(numr, "F")
	cnx = nx//2 + 1
	cny = nx//2 + 1
	qv = pi/180.0

	refrings = []     # list of (image objects) reference projections in Fourier representation

        for i in xrange( len(prjref) ):
		cimage = Util.Polar2Dm(filt_ctf(prjref[i], ctf, True) , cnx, cny, numr, mode)  # currently set to quadratic....
		Util.Normalize_ring(cimage, numr)

		Util.Frngs(cimage, numr)
		Util.Applyws(cimage, numr, wr_four)
		refrings.append(cimage)
		phi   = prjref[i].get_attr('phi')
		theta = prjref[i].get_attr('theta')
		psi   = prjref[i].get_attr('psi')
		n1 = sin(theta*qv)*cos(phi*qv)
		n2 = sin(theta*qv)*sin(phi*qv)
		n3 = cos(theta*qv)
		refrings[i].set_attr_dict( {"n1":n1, "n2":n2, "n3":n3, "phi": phi, "theta": theta,"psi": psi} )

	return refrings

def prgq( volft, kb, nx, delta, ref_a, sym, MPI=False):
	"""
	  Generate set of projections based on even angles
	  The command returns list of ffts of projections
	"""
	from projection   import prep_vol, prgs
	from applications import MPI_start_end
	from utilities    import even_angles, model_blank
	from fundamentals import fft
	# generate list of Eulerian angles for reference projections
	#  phi, theta, psi
	mode = "F"
	ref_angles = even_angles(delta, symmetry=sym, method = ref_a, phiEqpsi = "Minus")
	cnx = nx//2 + 1
	cny = nx//2 + 1
        num_ref = len(ref_angles)

	if MPI:
		from mpi import mpi_comm_rank, mpi_comm_size, MPI_COMM_WORLD
		myid = mpi_comm_rank( MPI_COMM_WORLD )
		ncpu = mpi_comm_size( MPI_COMM_WORLD )
	else:
		ncpu = 1
		myid = 0
	from applications import MPI_start_end
	ref_start,ref_end = MPI_start_end( num_ref, ncpu, myid )

	prjref = []     # list of (image objects) reference projections in Fourier representation

        for i in xrange(num_ref):
		prjref.append(model_blank(nx, nx))  # I am not sure why is that necessary, why not put None's??

        for i in xrange(ref_start, ref_end):
		prjref[i] = prgs(volft, kb, [ref_angles[i][0], ref_angles[i][1], ref_angles[i][2], 0.0, 0.0])

	if MPI:
		from utilities import bcast_EMData_to_all
		for i in xrange(num_ref):
			for j in xrange(ncpu):
				ref_start,ref_end = MPI_start_end(num_ref,ncpu,j)
				if i >= ref_start and i < ref_end: rootid = j
			bcast_EMData_to_all( prjref[i], myid, rootid )

	for i in xrange(len(ref_angles)):
		prjref[i].set_attr_dict({"phi": ref_angles[i][0], "theta": ref_angles[i][1],"psi": ref_angles[i][2]})

	return prjref


def prgs1d( prjft, kb, params ):
	from fundamentals import fft
	from math import cos, sin, pi
	from EMAN2 import Processor

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
	Mx=volume.get_xsize()
	My=volume.get_ysize()
	Mz=volume.get_zsize()
	if(Mx==Mz&My==Mz):
		volft,kb = prep_vol(volume)
		return  prgs(volft,kb,params)
	else:
		volft,kbx,kby,kbz = prep_vol(volume)
		return  prgs(volft,kbz,params,kbx,kby) 
	

def prep_vol(vol):
	"""
		Name
			prep_vol - prepare the volume for calculation of gridding projections and generate the interpolants.
		Input
			vol: input volume for which projections will be calculated using prgs
		Output
			volft: volume prepared for gridding projections using prgs
			kb: interpolants (tabulated Kaiser-Bessel function) when the volume is cubic.
			kbx,kby: interpolants along x, y and z direction (tabulated Kaiser-Bessel function) when the volume is rectangular 
	"""
	# prepare the volume
	Mx=vol.get_xsize()
	My=vol.get_ysize()
	Mz=vol.get_zsize()
	K     = 6
	alpha = 1.75
	npad  = 2
	if(Mx==Mz&My==Mz):
		M     = vol.get_xsize()
		# padd two times
		N     = M*npad
		# support of the window
		kb    = Util.KaiserBessel(alpha, K, M/2, K/(2.*N), N)
		volft = vol.copy()
		volft.divkbsinh(kb)
		volft = volft.norm_pad(False, npad)
		volft.do_fft_inplace()
		volft.center_origin_fft()
		volft.fft_shuffle()
		return  volft,kb
	else:
	
		# padd two times
		
		Nx     = Mx*npad
		Ny     = My*npad
		Nz     = Mz*npad
		# support of the window
		kbx    = Util.KaiserBessel(alpha, K, Mx/2, K/(2.*Nx), Nx)
		kby    = Util.KaiserBessel(alpha, K, My/2, K/(2.*Ny), Ny)
		kbz    = Util.KaiserBessel(alpha, K, Mz/2, K/(2.*Nz), Nz)
		volft = vol.copy()
		volft.divkbsinh_rect(kbx,kby,kbz)
		volft = volft.norm_pad(False, npad)
		volft.do_fft_inplace()
		volft.center_origin_fft()
		volft.fft_shuffle()
		return  volft,kbx,kby,kbz
		

###############################################################################################
## COMMON LINES NEW VERSION ###################################################################

# plot angles, map on half-sphere
# agls: [[phi0, theta0, psi0], [phi1, theta1, psi1], ..., [phin, thetan, psin]]
def plot_angles(agls, nx = 256):
	from math      import cos, sin, fmod, pi, radians
	from utilities import model_blank

	# var
	im = model_blank(nx, nx)
	"""
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
	"""
	# for each angles plot on circle area
	# agsl: [phi, theta, psi]
	ri = nx//2
	rr = ri-1
	conv = pi/180.0
	for i in xrange(len(agls)):
		if agls[i][1] > 90.0:
			agls[i][0] = agls[i][0] + 180.0
			agls[i][1] = 180.0 - float(agls[i][1])

		rc  = rr*sin( radians(agls[i][1]))
		rd  = radians(agls[i][0])
		px  = ri + rc * cos( rd )
		py  = ri + rc * sin( rd )

		px = min(max(int(px+0.5),0), nx-1)

		py = min(max(int(py+0.5),0), nx-1)

		im.set_value_at(px, py, 1.0 + im.get_value_at(px, py))

	return im

# interface between the simplex function to refine the angles and the function to compute the discrepancy
# not used yet
def cml_refine_agls_wrap(vec_in, data, flag_weights = False):
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
	disc = cml_disc(Prj, Ori, True, flag_weights)

	return -disc

# cml refines angles with simplex
# not used yet
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
def cml_disc(Prj, Ori, Rot, flag_weights=True):
	global g_n_prj, g_n_psi, g_n_lines, g_seq
	if flag_weights:
		cml = Util.cml_line_in3d(Ori, g_seq, g_n_prj, g_n_lines)
		weights = Util.cml_weights(cml)
		mw  = max(weights)
		for i in xrange(g_n_lines): weights[i]  = mw - weights[i]
		sw = sum(weights)
		if sw == 0:
			weights = [6.28 / float(g_n_lines)] * g_n_lines
		else:
			for i in xrange(g_n_lines):
				weights[i] /= sw
				weights[i] *= weights[i]
	else:   weights = [1.0] * g_n_lines

	com  = Util.cml_line_insino_all(Rot, g_seq, g_n_prj, g_n_lines)
	disc = Util.cml_disc(Prj, com, g_seq, weights, g_n_lines)

	return disc

# export the progress of the find_struc function
def cml_export_progress(outdir, ite, iprj, iagl, psi, disc, cmd):
	infofile = open(outdir + '/progress', 'a')
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


# display the list of angles for each iterations
def cml_export_txtagls(outdir, outname, Ori, disc, title):
	import time
	global g_n_prj, g_i_prj

	angfile = open(outdir + '/' + outname, 'a')

	angfile.write('|%s|-----------------------------------------------%s---------\n' % (title, time.ctime()))
	for i in xrange(g_n_prj): angfile.write('%10.3f\t%10.3f\t%10.3f\n' % (Ori[4*i], Ori[4*i+1], Ori[4*i+2]))
			
	angfile.write('\nDiscrepancy: %12.5e\n\n' % disc)
	angfile.close()

# init global variables used to a quick acces with many function of common-lines
def cml_init_global_var(dpsi, delta, nprj, debug):
	from utilities import even_angles

	global g_anglst, g_d_psi, g_n_psi, g_i_prj, g_n_lines, g_n_prj, g_n_anglst, g_debug, g_seq
	# TO FIX
	v = 180.0 / float(dpsi)
	if v != int(v):
		v    = int(v + 0.5)
		dpsi = 180 // v

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
		# horatio active_refactoring Jy51i1EwmLD4tWZ9_00000_1
		# data.set_attr('active', 1)
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
	from fundamentals import fftip
	from filter       import filt_tanh

	# number of projections
	if  type(stack) == type(""): nprj = EMUtil.get_image_count(stack)
	else:                       nprj = len(stack)
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
			diameter = int(2 * ou)
			mask2D   = model_circle(ou, nx, nx)
			if ir > 0:  mask2D -= model_circle(ir, nx, nx)

		# normalize under the mask
		[mean_a, sigma, imin, imax] = Util.infomask(image, mask2D, True)
		image -= mean_a
		Util.mul_scalar(image, 1.0/sigma)
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
		prj = model_blank(bdf, 2*nye)
		pp = model_blank(nxe, 2*nye)
		for li in xrange(nye):
			# get the line li
			line = Util.window(sino, nxe, 1, 1, 0, li-nye//2, 0)
			# u2 (not improve the results)
			#line = filt_tanh(line, ou / float(nx), ou / float(nx))
			# normalize this line
			[mean_l, sigma_l, imin, imax] = Util.infomask(line, None, True)
			line = (line - mean_l) / sigma_l
			# fft
			fftip(line)
			# filter (cut part of coef) and create mirror line
			Util.cml_prepare_line(prj, line, ilf, ihf, li, nye)

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
	nangle = int(180.0 / d_psi)
	dangle = M_PI / float(nangle)
	e = EMData()
	e.set_size(diameter, nangle, 1)
	offset = M - diameter // 2
	for j in xrange(nangle):
		nuxnew =  cos(dangle * j)
		nuynew = -sin(dangle * j)
		line   = fft(volft.extractline(kb, nuxnew, nuynew))
		Util.cyclicshift(line, {"dx":M, "dy":0, "dz":0})
		Util.set_line(e, j, line, offset, diameter)

	return e 

# transform an image to sinogram (mirror include)
def cml_sinogram_shift(image2D, diameter, shifts = [0.0, 0.0], d_psi = 1):
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
	#  Apply shift
	from EMAN2 import Processor
	params2 = {"filter_type" : Processor.fourier_filter_types.SHIFT, "x_shift" : 2*shifts[0], "y_shift" : 2*shifts[1], "z_shift" : 0.0}
	volft = Processor.EMFourierFilter(volft, params2)

	volft.center_origin_fft()
	volft.fft_shuffle()

	# get line projection
	nangle = int(180.0 / d_psi)
	dangle = M_PI / float(nangle)
	e = EMData()
	e.set_size(diameter, nangle, 1)
	offset = M - diameter // 2
	for j in xrange(nangle):
		nuxnew =  cos(dangle * j)
		nuynew = -sin(dangle * j)
		line   = fft(volft.extractline(kb, nuxnew, nuynew))
		Util.cyclicshift(line, {"dx":M, "dy":0, "dz":0})
		Util.set_line(e, j, line, offset, diameter)

	return e 

# write the head of the logfile
def cml_head_log(stack, outdir, delta, ir, ou, lf, hf, rand_seed, maxit, given, flag_weights, trials, ncpu):
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
	print_msg('Number of cpus               : %i\n'     % ncpu)
	print_msg('Use Voronoi weights          : %s\n\n'   % flag_weights)

# write the end of the logfile
def cml_end_log(Ori):
	from utilities import print_msg
	global g_n_prj
	print_msg('\n\n')
	for i in xrange(g_n_prj): print_msg('Projection #%03i: phi %10.5f    theta %10.5f    psi %10.5f\n' % (i, Ori[4*i], Ori[4*i+1], Ori[4*i+2]))

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

	# to stop when the solution oscillates
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

			# prepare active index of cml for weighting in order to earn time later
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
			best_disc = 1.0e20
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
						mw  = max(weights)
						for i in xrange(g_n_lines): weights[i]  = mw - weights[i]
						sw = sum(weights)
						if sw == 0:
							weights = [6.28 / float(g_n_lines)] * g_n_lines
						else:
							for i in xrange(g_n_lines):
								weights[i] /= sw
								weights[i] *= weights[i]
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

		# to stop when the solution oscillates
		period_disc.pop(0)
		period_disc.append(disc)
		if period_disc[0] == period_disc[2]:
			period_ct += 1
			if period_ct >= period_th and min(period_disc) == disc:
				angfile = open(outdir + '/' + outname, 'a')
				angfile.write('\nSTOP SOLUTION UNSTABLE\n')
				angfile.write('Discrepancy period: %s\n' % period_disc)
				angfile.close()
				break
		else:
			period_ct = 0

	return Ori, disc, ite

# find structure
def cml_find_structure2(Prj, Ori, Rot, outdir, outname, maxit, first_zero, flag_weights, myid, main_node, number_of_proc):
	from projection import cml_export_progress, cml_disc, cml_export_txtagls
	import time, sys
	from random import shuffle,random

	from mpi import MPI_FLOAT, MPI_INT, MPI_SUM, MPI_COMM_WORLD
	from mpi import mpi_reduce, mpi_bcast, mpi_barrier

	# global vars
	global g_i_prj, g_n_prj, g_n_anglst, g_anglst, g_d_psi, g_debug, g_n_lines, g_seq

	# list of free orientation
	ocp = [-1] * g_n_anglst

	if first_zero:
		listprj = range(1, g_n_prj)
		ocp[0]  = 0 
	else:   listprj = range(g_n_prj)

	# to stop when the solution oscillates
	period_disc = [0, 0, 0]
	period_ct   = 0
	period_th   = 2
	#if not flag_weights:   weights = [1.0] * g_n_lines

	# iteration loop
	for ite in xrange(maxit):
		#print ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>    ite = ", ite, "   myid = ", myid
		t_start = time.time()

		# loop over i prj
		change = False
		tlistprj = listprj[:]
		shuffle(tlistprj)
		nnn = len(tlistprj)
		tlistprj = mpi_bcast(tlistprj, nnn, MPI_INT, main_node, MPI_COMM_WORLD)
		tlistprj = map(int, tlistprj)
		"""
		if(ite>1 and ite%5 == 0  and ite<140):
			if(myid == main_node):
				for i in xrange(0,len(tlistprj),5):
					ind          = 4*i
					Ori[ind]      =  360.*random()
					Ori[ind+1]    =  180.*random()
					Ori[ind+2]    =  360.*random()
					Ori[ind+3]    =  -1
				for i in xrange(len(tlistprj)):
					ind          = 4*i
					Ori[ind+3]    = float(Ori[ind+3])
			nnn = len(Ori)
			Ori = mpi_bcast(Ori, nnn, MPI_FLOAT, main_node, MPI_COMM_WORLD)
			Ori = map(float, Ori)
			for i in xrange(len(tlistprj)):
				ind          = 4*i
				Ori[ind+3]    = int(Ori[ind+3])
		"""

		for iprj in tlistprj:
			#print "**********************************  iprj = ", iprj, g_n_anglst

			# Store current the current orientation
			ind          = 4*iprj
			store_phi    = Ori[ind]
			store_theta  = Ori[ind+1]
			store_psi    = Ori[ind+2]
			cur_agl      = Ori[ind+3]
			if cur_agl  != -1: ocp[cur_agl] = -1

			# prepare active index of cml for weighting in order to earn time later
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
			best_disc_list = [0]*g_n_anglst
			best_psi_list  = [0]*g_n_anglst
			for iagl in xrange(myid, g_n_anglst, number_of_proc):
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
						mw  = max(weights)
						for i in xrange(g_n_lines): weights[i]  = mw - weights[i]
						sw = sum(weights)
						if sw == 0:
							weights = [6.28 / float(g_n_lines)] * g_n_lines
						else:
							for i in xrange(g_n_lines):
								weights[i] /= sw
								weights[i] *= weights[i]

					# spin all psi
					com = Util.cml_line_insino(Rot, iprj, g_n_prj)
					if flag_weights:
						res = Util.cml_spin_psi(Prj, com, weights, iprj, iw, g_n_psi, g_d_psi, g_n_prj)
					else:
						res = Util.cml_spin_psi_now(Prj, com, iprj, iw, g_n_psi, g_d_psi, g_n_prj)

					# select the best
					best_disc_list[iagl] = res[0]
					best_psi_list[iagl]  = res[1]

					if g_debug: cml_export_progress(outdir, ite, iprj, iagl, res[1], res[0], 'progress')
				else:
					if g_debug: cml_export_progress(outdir, ite, iprj, iagl, -1, -1, 'progress')
			best_disc_list = mpi_reduce(best_disc_list, g_n_anglst, MPI_FLOAT, MPI_SUM, main_node, MPI_COMM_WORLD)
			best_psi_list = mpi_reduce(best_psi_list, g_n_anglst, MPI_FLOAT, MPI_SUM, main_node, MPI_COMM_WORLD)

			best_psi = -1
			best_iagl = -1

			if myid == main_node:
				best_disc = 1.0e20
				for iagl in xrange(g_n_anglst):
					if best_disc_list[iagl] > 0.0 and best_disc_list[iagl] < best_disc:
						best_disc = best_disc_list[iagl]
						best_psi = best_psi_list[iagl]
						best_iagl = iagl
			best_psi = mpi_bcast(best_psi, 1, MPI_FLOAT, main_node, MPI_COMM_WORLD)
			best_iagl = mpi_bcast(best_iagl, 1, MPI_INT, main_node, MPI_COMM_WORLD)
			best_psi = float(best_psi[0])
			best_iagl =  int(best_iagl[0])
			
			#print "xxxxx myid = ", myid, "    best_psi = ", best_psi, "   best_ialg = ", best_iagl

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
		if myid == main_node:
			cml_export_txtagls(outdir, outname, Ori, disc, 'Ite: %03i' % (ite + 1))

		if not change: break

		# to stop when the solution oscillates
		period_disc.pop(0)
		period_disc.append(disc)
		if period_disc[0] == period_disc[2]:
			period_ct += 1
			if period_ct >= period_th and min(period_disc) == disc and myid == main_node:
				angfile = open(outdir + '/' + outname, 'a')
				angfile.write('\nSTOP SOLUTION UNSTABLE\n')
				angfile.write('Discrepancy period: %s\n' % period_disc)
				angfile.close()
				break
		else:
			period_ct = 0
		mpi_barrier(MPI_COMM_WORLD)

	return Ori, disc, ite


# this function return the degree of colinearity of the orientations found (if colinear the value is close to zero)
def cml2_ori_collinearity(Ori):
	from math  import sin, cos, pi
	from numpy import array, linalg, matrix, zeros, power

	# ori 3d sphere map to 2d plan
	rad2deg = 180.0 / pi
	deg2rad = 1.0 / rad2deg
	nori    = len(Ori) // 4
	lx, ly  = [], []
	for n in xrange(nori):
		ind     = n * 4
		ph, th  = Ori[ind:ind+2]
		ph     *= deg2rad
		th     *= deg2rad
		lx.append(sin(th) * sin(ph))
		ly.append(sin(th) * cos(ph))

	# direct least square fitting of ellipse (IEEE ICIP Pilu96)
	D = zeros((nori, 6))
	for c in xrange(nori):
		D[c, :] = [lx[c]*lx[c], lx[c]*ly[c], ly[c]*ly[c], lx[c], ly[c], 1.0]
	D = matrix(D)
	S = D.transpose() * D
	C = zeros((6, 6))
	C[5, 5] =  0
	C[0, 2] = -2
	C[1, 2] =  1
	C[2, 0] = -2
	C = matrix(C)
	val, vec = linalg.eig(S.getI() * C)
	ell = vec[:, val.argmin()]
	verr = D * ell
	verr = power(verr, 2)
	serr = sum(verr)

	# sum squared error
	return serr.getA()[0][0]

## END COMMON LINES NEW VERSION ###############################################################
###############################################################################################

