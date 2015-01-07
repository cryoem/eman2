#!/usr/bin/env python

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

# clean up the code, make documentation

import os
import global_def
from   global_def     import *
from   user_functions import *
from   optparse       import OptionParser
import sys

def main():
	progname = os.path.basename(sys.argv[0])
	usage = progname + " stack <maskfile> --search_rng=10 --maxit=max_iteration --CTF --snr=SNR --Fourvar=Fourier_variance --oneDx --MPI"
	parser = OptionParser(usage,version=SPARXVERSION)
	parser.add_option("--search_rng",       type="int",           default=-1,      help="Search range for x-shift")
	parser.add_option("--search_ang",       type="int",           default=-1,      help="Search range for inplane rotation angle")
	parser.add_option("--search_rng_y",     type="int",           default=-1,      help="Search range for x-shift. Not used for 1D search (oneDx flag set).")
	parser.add_option("--maxit",            type="int",           default=100,     help="Maximum number of iterations program will perform")
	parser.add_option("--CTF",              action="store_true",  default=False,   help="Use CTF correction")
	parser.add_option("--snr",              type="float",         default=1.0,     help="signal-to-noise ratio of the data (default is 1.0)")
	parser.add_option("--Fourvar",          action="store_true",  default=False,   help="compute Fourier variance")
	parser.add_option("--oneDx",            action="store_true",  default=False,   help="1D search along x-axis")
	parser.add_option("--MPI",              action="store_true",  default=False,   help="use MPI")
	parser.add_option("--curvature",        action="store_true",  default=False,   help="for curved filament alignment")
	(options, args) = parser.parse_args()
	
	if not(options.MPI):
		print "Only MPI version is currently implemented."
		print "Please run '" + progname + " -h' for detailed options"
		return
			
	if len(args) < 1 or len(args) > 2:
		print "usage: " + usage
		print "Please run '" + progname + " -h' for detailed options"
	else:
	
		if len(args) == 1: mask = None
		else:              mask = args[1]
			
		if global_def.CACHE_DISABLE:
			from utilities import disable_bdb_cache
			disable_bdb_cache()
		
		from mpi import mpi_init
		sys.argv = mpi_init(len(sys.argv),sys.argv)

		global_def.BATCH = True
		if options.oneDx:
			if options.curvature: 
				snakehelicalshiftali_MPI(args[0], mask, options.maxit, options.CTF, options.snr, options.Fourvar, options.search_rng, options.search_ang)
				##curhelicalshiftali_MPI(args[0], mask, options.maxit, options.CTF, options.snr, options.Fourvar, options.search_rng)
			else:
				helicalshiftali_MPI(args[0], mask, options.maxit, options.CTF, options.snr, options.Fourvar, options.search_rng)		
		else:
			shiftali_MPI(args[0], mask, options.maxit, options.CTF, options.snr, options.Fourvar,options.search_rng,options.oneDx,options.search_rng_y)
		global_def.BATCH = False
		
		from mpi import mpi_finalize
		mpi_finalize()

def shiftali_MPI(stack, maskfile=None, maxit=100, CTF=False, snr=1.0, Fourvar=False, search_rng=-1, oneDx=False, search_rng_y=-1):  
	from applications import MPI_start_end
	from utilities    import model_circle, model_blank, get_image, peak_search, get_im
	from utilities    import reduce_EMData_to_root, bcast_EMData_to_all, send_attr_dict, file_type, bcast_number_to_all, bcast_list_to_all
	from statistics   import varf2d_MPI
	from fundamentals import fft, ccf, rot_shift3D, rot_shift2D
	from utilities    import get_params2D, set_params2D
	from utilities    import print_msg, print_begin_msg, print_end_msg
	import os
	import sys
	from mpi 	  	  import mpi_init, mpi_comm_size, mpi_comm_rank, MPI_COMM_WORLD
	from mpi 	  	  import mpi_reduce, mpi_bcast, mpi_barrier, mpi_gatherv
	from mpi 	  	  import MPI_SUM, MPI_FLOAT, MPI_INT
	from EMAN2	  	  import Processor
	from time         import time	
	
	number_of_proc = mpi_comm_size(MPI_COMM_WORLD)
	myid = mpi_comm_rank(MPI_COMM_WORLD)
	main_node = 0
		
	ftp = file_type(stack)

	if myid == main_node:
		print_begin_msg("shiftali_MPI")

	max_iter=int(maxit)

	if myid == main_node:
		if ftp == "bdb":
			from EMAN2db import db_open_dict
			dummy = db_open_dict(stack, True)
		active = EMUtil.get_all_attributes(stack, 'active')
		list_of_particles = []
		for im in xrange(len(active)):
			if active[im]:  list_of_particles.append(im)
		del active
		nima = len(list_of_particles)
	else:
		nima = 0
	nima = bcast_number_to_all(nima, source_node = main_node)

	if myid != main_node:
		list_of_particles = [-1]*nima
	list_of_particles = bcast_list_to_all(list_of_particles, source_node = main_node)

	image_start, image_end = MPI_start_end(nima, number_of_proc, myid)
	list_of_particles = list_of_particles[image_start: image_end]

	# read nx and ctf_app (if CTF) and broadcast to all nodes
	if myid == main_node:
		ima = EMData()
		ima.read_image(stack, list_of_particles[0], True)
		nx = ima.get_xsize()
		ny = ima.get_ysize()
		if CTF:	ctf_app = ima.get_attr_default('ctf_applied', 2)
		del ima
	else:
		nx = 0
		ny = 0
		if CTF:	ctf_app = 0
	nx = bcast_number_to_all(nx, source_node = main_node)
	ny = bcast_number_to_all(ny, source_node = main_node)
	if CTF:
		ctf_app = bcast_number_to_all(ctf_app, source_node = main_node)
		if ctf_app > 0:	ERROR("data cannot be ctf-applied", "shiftali_MPI", 1, myid)

	if maskfile == None:
		mrad = min(nx, ny)
		mask = model_circle(mrad//2-2, nx, ny)
	else:
		mask = get_im(maskfile)

	if CTF:
		from filter import filt_ctf
		from morphology   import ctf_img
		ctf_abs_sum = EMData(nx, ny, 1, False)
		ctf_2_sum = EMData(nx, ny, 1, False)
	else:
		ctf_2_sum = None

	from global_def import CACHE_DISABLE
	if CACHE_DISABLE:
		data = EMData.read_images(stack, list_of_particles)
	else:
		for i in xrange(number_of_proc):
			if myid == i:
				data = EMData.read_images(stack, list_of_particles)
			if ftp == "bdb": mpi_barrier(MPI_COMM_WORLD)


	for im in xrange(len(data)):
		data[im].set_attr('ID', list_of_particles[im])
		st = Util.infomask(data[im], mask, False)
		data[im] -= st[0]
		if CTF:
			ctf_params = data[im].get_attr("ctf")
			ctfimg = ctf_img(nx, ctf_params, ny=ny)
			Util.add_img2(ctf_2_sum, ctfimg)
			Util.add_img_abs(ctf_abs_sum, ctfimg)

	if CTF:
		reduce_EMData_to_root(ctf_2_sum, myid, main_node)
		reduce_EMData_to_root(ctf_abs_sum, myid, main_node)
	else:  ctf_2_sum = None
	if CTF:
		if myid != main_node:
			del ctf_2_sum
			del ctf_abs_sum
		else:
			temp = EMData(nx, ny, 1, False)
			for i in xrange(0,nx,2):
				for j in xrange(ny):
					temp.set_value_at(i,j,snr)
			Util.add_img(ctf_2_sum, temp)
			del temp

	total_iter = 0

	# apply initial xform.align2d parameters stored in header
	init_params = []
	for im in xrange(len(data)):
		t = data[im].get_attr('xform.align2d')
		init_params.append(t)
		p = t.get_params("2d")
		data[im] = rot_shift2D(data[im], p['alpha'], sx=p['tx'], sy=p['ty'], mirror=p['mirror'], scale=p['scale'])

	# fourier transform all images, and apply ctf if CTF
	for im in xrange(len(data)):
		if CTF:
			ctf_params = data[im].get_attr("ctf")
			data[im] = filt_ctf(fft(data[im]), ctf_params)
		else:
			data[im] = fft(data[im])

	sx_sum=0
	sy_sum=0
	sx_sum_total=0
	sy_sum_total=0
	shift_x = [0.0]*len(data)
	shift_y = [0.0]*len(data)
	ishift_x = [0.0]*len(data)
	ishift_y = [0.0]*len(data)

	for Iter in xrange(max_iter):
		if myid == main_node:
			start_time = time()
			print_msg("Iteration #%4d\n"%(total_iter))
		total_iter += 1
		avg = EMData(nx, ny, 1, False)
		for im in data:  Util.add_img(avg, im)

		reduce_EMData_to_root(avg, myid, main_node)

		if myid == main_node:
			if CTF:
				tavg = Util.divn_filter(avg, ctf_2_sum)
			else:	 tavg = Util.mult_scalar(avg, 1.0/float(nima))
		else:
			tavg = EMData(nx, ny, 1, False)                               

		if Fourvar:
			bcast_EMData_to_all(tavg, myid, main_node)
			vav, rvar = varf2d_MPI(myid, data, tavg, mask, "a", CTF)

		if myid == main_node:
			if Fourvar:
				tavg    = fft(Util.divn_img(fft(tavg), vav))
				vav_r	= Util.pack_complex_to_real(vav)

			# normalize and mask tavg in real space
			tavg = fft(tavg)
			stat = Util.infomask( tavg, mask, False ) 
			tavg -= stat[0]
			Util.mul_img(tavg, mask)
			# For testing purposes: shift tavg to some random place and see if the centering is still correct
			#tavg = rot_shift3D(tavg,sx=3,sy=-4)
			tavg = fft(tavg)

		if Fourvar:  del vav
		bcast_EMData_to_all(tavg, myid, main_node)

		sx_sum=0 
		sy_sum=0 
		if search_rng > 0: nwx = 2*search_rng+1
		else:              nwx = nx
		
		if search_rng_y > 0: nwy = 2*search_rng_y+1
		else:                nwy = ny

		not_zero = 0
		for im in xrange(len(data)):
			if oneDx:
				ctx = Util.window(ccf(data[im],tavg),nwx,1)
				p1  = peak_search(ctx)
				p1_x = -int(p1[0][3])
				ishift_x[im] = p1_x
				sx_sum += p1_x
			else:
				p1 = peak_search(Util.window(ccf(data[im],tavg), nwx,nwy))
				p1_x = -int(p1[0][4])
				p1_y = -int(p1[0][5])
				ishift_x[im] = p1_x
				ishift_y[im] = p1_y
				sx_sum += p1_x
				sy_sum += p1_y

			if not_zero == 0:
				if (not(ishift_x[im] == 0.0)) or (not(ishift_y[im] == 0.0)):
					not_zero = 1

		sx_sum = mpi_reduce(sx_sum, 1, MPI_INT, MPI_SUM, main_node, MPI_COMM_WORLD)  

		if not oneDx:
			sy_sum = mpi_reduce(sy_sum, 1, MPI_INT, MPI_SUM, main_node, MPI_COMM_WORLD)

		if myid == main_node:
			sx_sum_total = int(sx_sum[0])
			if not oneDx:
				sy_sum_total = int(sy_sum[0])
		else:
			sx_sum_total = 0	
			sy_sum_total = 0

		sx_sum_total = bcast_number_to_all(sx_sum_total, source_node = main_node)

		if not oneDx:
			sy_sum_total = bcast_number_to_all(sy_sum_total, source_node = main_node)

		sx_ave = round(float(sx_sum_total)/nima)
		sy_ave = round(float(sy_sum_total)/nima)
		for im in xrange(len(data)): 
			p1_x = ishift_x[im] - sx_ave
			p1_y = ishift_y[im] - sy_ave
			params2 = {"filter_type" : Processor.fourier_filter_types.SHIFT, "x_shift" : p1_x, "y_shift" : p1_y, "z_shift" : 0.0}
			data[im] = Processor.EMFourierFilter(data[im], params2)
			shift_x[im] += p1_x
			shift_y[im] += p1_y
		# stop if all shifts are zero
		not_zero = mpi_reduce(not_zero, 1, MPI_INT, MPI_SUM, main_node, MPI_COMM_WORLD)  
		if myid == main_node:
			not_zero_all = int(not_zero[0])
		else:
			not_zero_all = 0
		not_zero_all = bcast_number_to_all(not_zero_all, source_node = main_node)

		if myid == main_node:
			print_msg("Time of iteration = %12.2f\n"%(time()-start_time))
			start_time = time()

		if not_zero_all == 0:  break

	#for im in xrange(len(data)): data[im] = fft(data[im])  This should not be required as only header information is used
	# combine shifts found with the original parameters
	for im in xrange(len(data)):		
		t0 = init_params[im]
		t1 = Transform()
		t1.set_params({"type":"2D","alpha":0,"scale":t0.get_scale(),"mirror":0,"tx":shift_x[im],"ty":shift_y[im]})
		# combine t0 and t1
		tt = t1*t0
		data[im].set_attr("xform.align2d", tt)  

	# write out headers and STOP, under MPI writing has to be done sequentially
	mpi_barrier(MPI_COMM_WORLD)
	par_str = ["xform.align2d", "ID"]
	if myid == main_node:
		from utilities import file_type
		if(file_type(stack) == "bdb"):
			from utilities import recv_attr_dict_bdb
			recv_attr_dict_bdb(main_node, stack, data, par_str, image_start, image_end, number_of_proc)
		else:
			from utilities import recv_attr_dict
			recv_attr_dict(main_node, stack, data, par_str, image_start, image_end, number_of_proc)
		
	else:           send_attr_dict(main_node, data, par_str, image_start, image_end)
	if myid == main_node: print_end_msg("shiftali_MPI")				


		
def helicalshiftali_MPI(stack, maskfile=None, maxit=100, CTF=False, snr=1.0, Fourvar=False, search_rng=-1):
	from applications import MPI_start_end
	from utilities    import model_circle, model_blank, get_image, peak_search, get_im, pad
	from utilities    import reduce_EMData_to_root, bcast_EMData_to_all, send_attr_dict, file_type, bcast_number_to_all, bcast_list_to_all
	from statistics   import varf2d_MPI
	from fundamentals import fft, ccf, rot_shift3D, rot_shift2D, fshift
	from utilities    import get_params2D, set_params2D, chunks_distribution
	from utilities    import print_msg, print_begin_msg, print_end_msg
	import os
	import sys
	from mpi 	  	  import mpi_init, mpi_comm_size, mpi_comm_rank, MPI_COMM_WORLD
	from mpi 	  	  import mpi_reduce, mpi_bcast, mpi_barrier, mpi_gatherv
	from mpi 	  	  import MPI_SUM, MPI_FLOAT, MPI_INT
	from time         import time	
	from pixel_error  import ordersegments
	from math         import sqrt, atan2, tan, pi
	
	nproc = mpi_comm_size(MPI_COMM_WORLD)
	myid = mpi_comm_rank(MPI_COMM_WORLD)
	main_node = 0
		
	ftp = file_type(stack)

	if myid == main_node:
		print_begin_msg("helical-shiftali_MPI")

	max_iter=int(maxit)
	if( myid == main_node):
		infils = EMUtil.get_all_attributes(stack, "filament")
		ptlcoords = EMUtil.get_all_attributes(stack, 'ptcl_source_coord')
		filaments = ordersegments(infils, ptlcoords)
		total_nfils = len(filaments)
		inidl = [0]*total_nfils
		for i in xrange(total_nfils):  inidl[i] = len(filaments[i])
		linidl = sum(inidl)
		nima = linidl
		tfilaments = []
		for i in xrange(total_nfils):  tfilaments += filaments[i]
		del filaments
	else:
		total_nfils = 0
		linidl = 0
	total_nfils = bcast_number_to_all(total_nfils, source_node = main_node)
	if myid != main_node:
		inidl = [-1]*total_nfils
	inidl = bcast_list_to_all(inidl, source_node = main_node)
	linidl = bcast_number_to_all(linidl, source_node = main_node)
	if myid != main_node:
		tfilaments = [-1]*linidl
	tfilaments = bcast_list_to_all(tfilaments, source_node = main_node)
	filaments = []
	iendi = 0
	for i in xrange(total_nfils):
		isti = iendi
		iendi = isti+inidl[i]
		filaments.append(tfilaments[isti:iendi])
	del tfilaments,inidl

	if myid == main_node:
		print "total number of filaments: ", total_nfils
	if total_nfils< nproc:
		ERROR('number of CPUs (%i) is larger than the number of filaments (%i), please reduce the number of CPUs used'%(nproc, total_nfils), "ehelix_MPI", 1,myid)

	#  balanced load
	temp = chunks_distribution([[len(filaments[i]), i] for i in xrange(len(filaments))], nproc)[myid:myid+1][0]
	filaments = [filaments[temp[i][1]] for i in xrange(len(temp))]
	nfils     = len(filaments)

	#filaments = [[0,1]]
	#print "filaments",filaments
	list_of_particles = []
	indcs = []
	k = 0
	for i in xrange(nfils):
		list_of_particles += filaments[i]
		k1 = k+len(filaments[i])
		indcs.append([k,k1])
		k = k1
	data = EMData.read_images(stack, list_of_particles)
	ldata = len(data)
	print "ldata=", ldata
	nx = data[0].get_xsize()
	ny = data[0].get_ysize()
	if maskfile == None:
		mrad = min(nx, ny)//2-2
		mask = pad( model_blank(2*mrad+1, ny, 1, 1.0), nx, ny, 1, 0.0)
	else:
		mask = get_im(maskfile)

	# apply initial xform.align2d parameters stored in header
	init_params = []
	for im in xrange(ldata):
		t = data[im].get_attr('xform.align2d')
		init_params.append(t)
		p = t.get_params("2d")
		data[im] = rot_shift2D(data[im], p['alpha'], p['tx'], p['ty'], p['mirror'], p['scale'])

	if CTF:
		from filter import filt_ctf
		from morphology   import ctf_img
		ctf_abs_sum = EMData(nx, ny, 1, False)
		ctf_2_sum = EMData(nx, ny, 1, False)
	else:
		ctf_2_sum = None
		ctf_abs_sum = None



	from utilities import info

	for im in xrange(ldata):
		data[im].set_attr('ID', list_of_particles[im])
		st = Util.infomask(data[im], mask, False)
		data[im] -= st[0]
		if CTF:
			ctf_params = data[im].get_attr("ctf")
			qctf = data[im].get_attr("ctf_applied")
			if qctf == 0:
				data[im] = filt_ctf(fft(data[im]), ctf_params)
				data[im].set_attr('ctf_applied', 1)
			elif qctf != 1:
				ERROR('Incorrectly set qctf flag', "helicalshiftali_MPI", 1,myid)
			ctfimg = ctf_img(nx, ctf_params, ny=ny)
			Util.add_img2(ctf_2_sum, ctfimg)
			Util.add_img_abs(ctf_abs_sum, ctfimg)
		else:  data[im] = fft(data[im])

	del list_of_particles		

	if CTF:
		reduce_EMData_to_root(ctf_2_sum, myid, main_node)
		reduce_EMData_to_root(ctf_abs_sum, myid, main_node)
	if CTF:
		if myid != main_node:
			del ctf_2_sum
			del ctf_abs_sum
		else:
			temp = EMData(nx, ny, 1, False)
			tsnr = 1./snr
			for i in xrange(0,nx+2,2):
				for j in xrange(ny):
					temp.set_value_at(i,j,tsnr)
					temp.set_value_at(i+1,j,0.0)
			#info(ctf_2_sum)
			Util.add_img(ctf_2_sum, temp)
			#info(ctf_2_sum)
			del temp

	total_iter = 0
	shift_x = [0.0]*ldata

	for Iter in xrange(max_iter):
		if myid == main_node:
			start_time = time()
			print_msg("Iteration #%4d\n"%(total_iter))
		total_iter += 1
		avg = EMData(nx, ny, 1, False)
		for im in xrange(ldata):
			Util.add_img(avg, fshift(data[im], shift_x[im]))

		reduce_EMData_to_root(avg, myid, main_node)

		if myid == main_node:
			if CTF:  tavg = Util.divn_filter(avg, ctf_2_sum)
			else:    tavg = Util.mult_scalar(avg, 1.0/float(nima))
		else:
			tavg = model_blank(nx,ny)

		if Fourvar:
			bcast_EMData_to_all(tavg, myid, main_node)
			vav, rvar = varf2d_MPI(myid, data, tavg, mask, "a", CTF)

		if myid == main_node:
			if Fourvar:
				tavg    = fft(Util.divn_img(fft(tavg), vav))
				vav_r	= Util.pack_complex_to_real(vav)
			# normalize and mask tavg in real space
			tavg = fft(tavg)
			stat = Util.infomask( tavg, mask, False )
			tavg -= stat[0]
			Util.mul_img(tavg, mask)
			tavg.write_image("tavg.hdf",Iter)
			# For testing purposes: shift tavg to some random place and see if the centering is still correct
			#tavg = rot_shift3D(tavg,sx=3,sy=-4)

		if Fourvar:  del vav
		bcast_EMData_to_all(tavg, myid, main_node)
		tavg = fft(tavg)

		sx_sum = 0.0
		nxc = nx//2
		
		for ifil in xrange(nfils):
			"""
			# Calculate filament average
			avg = EMData(nx, ny, 1, False)
			filnima = 0
			for im in xrange(indcs[ifil][0], indcs[ifil][1]):
				Util.add_img(avg, data[im])
				filnima += 1
			tavg = Util.mult_scalar(avg, 1.0/float(filnima))
			"""
			# Calculate 1D ccf between each segment and filament average
			nsegms = indcs[ifil][1]-indcs[ifil][0]
			ctx = [None]*nsegms
			pcoords = [None]*nsegms
			for im in xrange(indcs[ifil][0], indcs[ifil][1]):
				ctx[im-indcs[ifil][0]] = Util.window(ccf(tavg, data[im]), nx, 1)
				pcoords[im-indcs[ifil][0]] = data[im].get_attr('ptcl_source_coord')
				#ctx[im-indcs[ifil][0]].write_image("ctx.hdf",im-indcs[ifil][0])
				#print "  CTX  ",myid,im,Util.infomask(ctx[im-indcs[ifil][0]], None, True)
			# search for best x-shift
			cents = nsegms//2
			
			dst = sqrt(max((pcoords[cents][0] - pcoords[0][0])**2 + (pcoords[cents][1] - pcoords[0][1])**2, (pcoords[cents][0] - pcoords[-1][0])**2 + (pcoords[cents][1] - pcoords[-1][1])**2))
			maxincline = atan2(ny//2-2-float(search_rng),dst)
			kang = int(dst*tan(maxincline)+0.5)
			#print  "  settings ",nsegms,cents,dst,search_rng,maxincline,kang
			
			# ## C code for alignment. @ming
 			results = [0.0]*3;
 			results = Util.helixshiftali(ctx, pcoords, nsegms, maxincline, kang, search_rng,nxc)
			sib = int(results[0])
 			bang = results[1]
 			qm = results[2]
			#print qm, sib, bang
			
			# qm = -1.e23	
# 				
# 			for six in xrange(-search_rng, search_rng+1,1):
# 				q0 = ctx[cents].get_value_at(six+nxc)
# 				for incline in xrange(kang+1):
# 					qt = q0
# 					qu = q0
# 					if(kang>0):  tang = tan(maxincline/kang*incline)
# 					else:        tang = 0.0
# 					for kim in xrange(cents+1,nsegms):
# 						dst = sqrt((pcoords[cents][0] - pcoords[kim][0])**2 + (pcoords[cents][1] - pcoords[kim][1])**2)
# 						xl = dst*tang+six+nxc
# 						ixl = int(xl)
# 						dxl = xl - ixl
# 						#print "  A  ", ifil,six,incline,kim,xl,ixl,dxl
# 						qt += (1.0-dxl)*ctx[kim].get_value_at(ixl) + dxl*ctx[kim].get_value_at(ixl+1)
# 						xl = -dst*tang+six+nxc
# 						ixl = int(xl)
# 						dxl = xl - ixl
# 						qu += (1.0-dxl)*ctx[kim].get_value_at(ixl) + dxl*ctx[kim].get_value_at(ixl+1)
# 					for kim in xrange(cents):
# 						dst = sqrt((pcoords[cents][0] - pcoords[kim][0])**2 + (pcoords[cents][1] - pcoords[kim][1])**2)
# 						xl = -dst*tang+six+nxc
# 						ixl = int(xl)
# 						dxl = xl - ixl
# 						qt += (1.0-dxl)*ctx[kim].get_value_at(ixl) + dxl*ctx[kim].get_value_at(ixl+1)
# 						xl =  dst*tang+six+nxc
# 						ixl = int(xl)
# 						dxl = xl - ixl
# 						qu += (1.0-dxl)*ctx[kim].get_value_at(ixl) + dxl*ctx[kim].get_value_at(ixl+1)
# 					if( qt > qm ):
# 						qm = qt
# 						sib = six
# 						bang = tang
# 					if( qu > qm ):
# 						qm = qu
# 						sib = six
# 						bang = -tang
					#if incline == 0:  print  "incline = 0  ",six,tang,qt,qu
			#print qm,six,sib,bang
			#print " got results   ",indcs[ifil][0], indcs[ifil][1], ifil,myid,qm,sib,tang,bang,len(ctx),Util.infomask(ctx[0], None, True)
			for im in xrange(indcs[ifil][0], indcs[ifil][1]):
				kim = im-indcs[ifil][0]
				dst = sqrt((pcoords[cents][0] - pcoords[kim][0])**2 + (pcoords[cents][1] - pcoords[kim][1])**2)
				if(kim < cents):  xl = -dst*bang+sib
				else:             xl =  dst*bang+sib
				shift_x[im] = xl
							
			# Average shift
			sx_sum += shift_x[indcs[ifil][0]+cents]
			
			
		# #print myid,sx_sum,total_nfils
		sx_sum = mpi_reduce(sx_sum, 1, MPI_FLOAT, MPI_SUM, main_node, MPI_COMM_WORLD)
		if myid == main_node:
			sx_sum = float(sx_sum[0])/total_nfils
			print_msg("Average shift  %6.2f\n"%(sx_sum))
		else:
			sx_sum = 0.0
		sx_sum = 0.0
		sx_sum = bcast_number_to_all(sx_sum, source_node = main_node)
		for im in xrange(ldata):
			shift_x[im] -= sx_sum
			#print  "   %3d  %6.3f"%(im,shift_x[im])
		#exit()


			
	# combine shifts found with the original parameters
	for im in xrange(ldata):		
		t1 = Transform()
		##import random
		##shix=random.randint(-10, 10)
		##t1.set_params({"type":"2D","tx":shix})
		t1.set_params({"type":"2D","tx":shift_x[im]})
		# combine t0 and t1
		tt = t1*init_params[im]
		data[im].set_attr("xform.align2d", tt)
	# write out headers and STOP, under MPI writing has to be done sequentially
	mpi_barrier(MPI_COMM_WORLD)
	par_str = ["xform.align2d", "ID"]
	if myid == main_node:
		from utilities import file_type
		if(file_type(stack) == "bdb"):
			from utilities import recv_attr_dict_bdb
			recv_attr_dict_bdb(main_node, stack, data, par_str, 0, ldata, nproc)
		else:
			from utilities import recv_attr_dict
			recv_attr_dict(main_node, stack, data, par_str, 0, ldata, nproc)
	else:           send_attr_dict(main_node, data, par_str, 0, ldata)
	if myid == main_node: print_end_msg("helical-shiftali_MPI")				


def snakehelicalshiftali_MPI(stack, maskfile=None, maxit=100, CTF=False, snr=1.0, Fourvar=False, search_rng=-1, search_ang=-1):
	from applications import MPI_start_end
	from utilities    import model_circle, model_blank, get_image, peak_search, get_im, pad
	from utilities    import reduce_EMData_to_root, bcast_EMData_to_all, send_attr_dict, file_type, bcast_number_to_all, bcast_list_to_all
	from statistics   import varf2d_MPI
	from fundamentals import fft, ccf, rot_shift3D, rot_shift2D, fshift
	from utilities    import get_params2D, set_params2D, chunks_distribution
	from utilities    import print_msg, print_begin_msg, print_end_msg
	import os
	import sys
	from mpi 	  	  import mpi_init, mpi_comm_size, mpi_comm_rank, MPI_COMM_WORLD
	from mpi 	  	  import mpi_reduce, mpi_bcast, mpi_barrier, mpi_gatherv
	from mpi 	  	  import MPI_SUM, MPI_FLOAT, MPI_INT
	from time         import time	
	from pixel_error  import ordersegments
	from math         import sqrt, atan2, tan, pi, asin
	
	nproc = mpi_comm_size(MPI_COMM_WORLD)
	myid = mpi_comm_rank(MPI_COMM_WORLD)
	main_node = 0
		
	ftp = file_type(stack)

	if myid == main_node:
		print_begin_msg("snake-shiftali_MPI")

	max_iter=int(maxit)
	if( myid == main_node):
		infils = EMUtil.get_all_attributes(stack, "filament")
		ptlcoords = EMUtil.get_all_attributes(stack, 'ptcl_source_coord')
		filaments = ordersegments(infils, ptlcoords)
		total_nfils = len(filaments)
		inidl = [0]*total_nfils
		for i in xrange(total_nfils):  inidl[i] = len(filaments[i])
		linidl = sum(inidl)
		nima = linidl
		tfilaments = []
		for i in xrange(total_nfils):  tfilaments += filaments[i]
		del filaments
	else:
		total_nfils = 0
		linidl = 0
	total_nfils = bcast_number_to_all(total_nfils, source_node = main_node)
	if myid != main_node:
		inidl = [-1]*total_nfils
	inidl = bcast_list_to_all(inidl, source_node = main_node)
	linidl = bcast_number_to_all(linidl, source_node = main_node)
	if myid != main_node:
		tfilaments = [-1]*linidl
	tfilaments = bcast_list_to_all(tfilaments, source_node = main_node)
	filaments = []
	iendi = 0
	for i in xrange(total_nfils):
		isti = iendi
		iendi = isti+inidl[i]
		filaments.append(tfilaments[isti:iendi])
	del tfilaments,inidl

	if myid == main_node:
		print "total number of filaments: ", total_nfils
	if total_nfils< nproc:
		ERROR('number of CPUs (%i) is larger than the number of filaments (%i), please reduce the number of CPUs used'%(nproc, total_nfils), "ehelix_MPI", 1,myid)

	#  balanced load
	temp = chunks_distribution([[len(filaments[i]), i] for i in xrange(len(filaments))], nproc)[myid:myid+1][0]
	filaments = [filaments[temp[i][1]] for i in xrange(len(temp))]
	nfils     = len(filaments)

	#filaments = [[0,1]]
	#print "filaments",filaments
	list_of_particles = []
	indcs = []
	k = 0
	for i in xrange(nfils):
		list_of_particles += filaments[i]
		k1 = k+len(filaments[i])
		indcs.append([k,k1])
		k = k1
	data = EMData.read_images(stack, list_of_particles)
	ldata = len(data)
	print "ldata=", ldata
	nx = data[0].get_xsize()
	ny = data[0].get_ysize()
	if maskfile == None:
		mrad = min(nx, ny)//2-2
		mask = pad( model_blank(2*mrad+1, ny, 1, 1.0), nx, ny, 1, 0.0)
	else:
		mask = get_im(maskfile)

	# apply initial xform.align2d parameters stored in header
	init_params = []
	for im in xrange(ldata):
		t = data[im].get_attr('xform.align2d')
		init_params.append(t)
		p = t.get_params("2d")
		#data[im] = rot_shift2D(data[im], p['alpha'], p['tx'], p['ty'], p['mirror'], p['scale'])  @ming

	if CTF:
		from filter import filt_ctf
		from morphology   import ctf_img
		ctf_abs_sum = EMData(nx, ny, 1, False)
		ctf_2_sum = EMData(nx, ny, 1, False)
	else:
		ctf_2_sum = None
		ctf_abs_sum = None



	from utilities import info

	for im in xrange(ldata):
		data[im].set_attr('ID', list_of_particles[im])
		st = Util.infomask(data[im], mask, False)
		data[im] -= st[0]
		if CTF:
			ctf_params = data[im].get_attr("ctf")
			qctf = data[im].get_attr("ctf_applied")
			if qctf == 0:
				data[im] = filt_ctf(fft(data[im]), ctf_params)
				data[im].set_attr('ctf_applied', 1)
			elif qctf != 1:
				ERROR('Incorrectly set qctf flag', "helicalshiftali_MPI", 1,myid)
			ctfimg = ctf_img(nx, ctf_params, ny=ny)
			Util.add_img2(ctf_2_sum, ctfimg)
			Util.add_img_abs(ctf_abs_sum, ctfimg)
		else:  data[im] = fft(data[im])

	del list_of_particles		

	if CTF:
		reduce_EMData_to_root(ctf_2_sum, myid, main_node)
		reduce_EMData_to_root(ctf_abs_sum, myid, main_node)
	if CTF:
		if myid != main_node:
			del ctf_2_sum
			del ctf_abs_sum
		else:
			temp = EMData(nx, ny, 1, False)
			tsnr = 1./snr
			for i in xrange(0,nx+2,2):
				for j in xrange(ny):
					temp.set_value_at(i,j,tsnr)
					temp.set_value_at(i+1,j,0.0)
			#info(ctf_2_sum)
			Util.add_img(ctf_2_sum, temp)
			#info(ctf_2_sum)
			del temp

	total_iter = 0
	shift_x = [0.0]*ldata
	shift_ang = [0.0]*ldata
	
	# center in SPIDER convention
	cnx = nx//2+1
	cny = ny//2+1
	rstep=1
	mode = "F"
	T = 1.0
	xrng = search_rng
	yrng = search_rng
	first_ring=1
	last_ring = nx/2-2-int(max(xrng,yrng))
	from alignment import Numrinit, ringwe
	numr = Numrinit(first_ring, last_ring, rstep, mode)
	wr   = ringwe(numr, mode)
	maxrin = numr[(last_ring-first_ring)*3+2]
	dang = 2*pi/maxrin
	kx = int(2*xrng/rstep+0.5)/2
	vol = maxrin*(2*kx+1)
	angnxc = maxrin//2
	
	resamp_dang = 2*asin(1.0/nx)
	resamp_maxrin = int(2*pi / resamp_dang+0.5)
	resamp_vol = resamp_maxrin * (2*kx+1)
	resamp_angnxc = resamp_maxrin//2
	print "dang=%f resamp_dfang=%f resamp_maxrin=%d"%(dang,resamp_dang, resamp_maxrin)
	
	print "max ring points=%d"%numr[len(numr)-1]
	
	total_iter = 0
	##tavg = fft(tavg)                                       #transform tavg    into real space.
	tttt = EMData(nx, ny, 1, False)
	
	## compute 2D ccfs.
	#crefim = Util.Polar2Dm(tavg, cnx, cny, numr, mode)
	#Util.Frngs(crefim, numr)
	##Util.Applyws(crefim, numr, wr)
	CCF2d = []
	for ifil in xrange(nfils):
		# test Calculate 2D ccf between each segment and filament average
		nsegms = indcs[ifil][1]-indcs[ifil][0]
		ccf2d = [None]*vol
		ctx2d = [None]*nsegms
		resamp_ccf2d = [None]*resamp_vol
		for im in xrange(indcs[ifil][0], indcs[ifil][1]):
			ttavg = get_im("bdb:vdata",im)
			stat = Util.infomask( ttavg, mask, False )
			ttavg -= stat[0]
			Util.mul_img(ttavg, mask)
			crefim = Util.Polar2Dm(ttavg, cnx, cny, numr, mode)
			Util.Frngs(crefim, numr)
			#Util.Applyws(crefim, numr, wr)
			tttt = fft(data[im])             # transform data[im] into real space.
			ccf2d = Util.ali2d_ccf_list_snake(tttt, crefim,  wr, xrng, yrng, rstep, mode, numr, cnx, cny, T)
			for i in xrange(2*kx+1):
				for j in xrange(resamp_maxrin):
					j_old = int(j * resamp_dang/dang + 0.5)
					resamp_ccf2d[i*resamp_maxrin+j] = ccf2d[i*maxrin+j_old]
			# aaaa=resamp_ccf2d.index(max(resamp_ccf2d))
# 			iaaa = aaaa/resamp_maxrin
# 			jaaa = aaaa%resamp_maxrin
# 			shift_x[im] = iaaa-(2*kx+1)//2
# 			shift_ang[im] = jaaa
# 			print "im=%d rotang=%d  maxrin=%d shift=%f"%(im, jaaa, resamp_maxrin, -shift_x[im] )	
			# if im == 7 :
# 						print "rotang=%f  shift=%f"%(jaaa*resamp_dang*180/pi, -shift_x[im] )	
# 						tttt = fft(data[im])
# 						tttt.write_image("image7.hdf")
# 						tttt=rot_shift2D(tttt, jaaa*resamp_dang*180/pi, -shift_x[im], 0, 0, 1)
# 						tttt.write_image("imagerot7.hdf")
			ccf2dimg = model_blank(2*kx+1, resamp_maxrin) ##EMData(2*kx+1, resamp_maxrin, 1, False)			
			for i in xrange(2*kx+1):
				for j in xrange(resamp_maxrin):
					ccf2dimg.set_value_at(i,j,resamp_ccf2d[i*resamp_maxrin+j])
			ctx2d[im-indcs[ifil][0]] = 	ccf2dimg	 
		CCF2d.append(ctx2d)
			
	
	## ---------------------------------------------- ##
	## step 1: find a straight line for initial guess ##
	## ---------------------------------------------- ##
	# for Iter in xrange(max_iter):
# 		nxc = nx//2
# 		sx_sum  = 0
# 		CCF = [] ##added@ming
# 		for ifil in xrange(nfils):
# 			# Calculate 1D ccf between each segment and filament average
# 			nsegms = indcs[ifil][1]-indcs[ifil][0]
# 			#sccf = []      ##added@ming
# 			ctx = [None]*nsegms
# 			pcoords = [None]*nsegms
# 			for im in xrange(indcs[ifil][0], indcs[ifil][1]):
# 				ttavg = get_im("bdb:vdata",im)
# 				stat = Util.infomask( ttavg, mask, False )
# 				ttavg -= stat[0]
# 				Util.mul_img(ttavg, mask)
# 				ttavg = fft(ttavg)
# 				ctx[im-indcs[ifil][0]] = Util.window(ccf(ttavg, data[im]), nx, 1)
# 				#sccf.append(ctx[im-indcs[ifil][0]])
# 				pcoords[im-indcs[ifil][0]] = data[im].get_attr('ptcl_source_coord')
# 				#ctx[im-indcs[ifil][0]].write_image("ctx.hdf",im-indcs[ifil][0])
# 				#print "  CTX  ",myid,im,Util.infomask(ctx[im-indcs[ifil][0]], None, True)
# 			CCF.append(ctx)	
# 			# search for best x-shift
# 			cents = nsegms//2
# 		
# 			dst = sqrt(max((pcoords[cents][0] - pcoords[0][0])**2 + (pcoords[cents][1] - pcoords[0][1])**2, (pcoords[cents][0] - pcoords[-1][0])**2 + (pcoords[cents][1] - pcoords[-1][1])**2))
# 			maxincline = atan2(ny//2-2-float(search_rng),dst)
# 			kang = int(dst*tan(maxincline)+0.5)
# 			#print  "  settings ",nsegms,cents,dst,search_rng,maxincline,kang
# 		
# 			# ## C code for alignment. @ming
# 			results = [0.0]*3;
# 			results = Util.helixshiftali(ctx, pcoords, nsegms, maxincline, kang, search_rng,nxc)
# 			sib = int(results[0])
# 			bang = results[1]
# 			qm = results[2]
# 			#print qm, sib, bang
# 		
# 			# qm = -1.e23	
# 
# 			#print qm,six,sib,bang
# 			#print " got results   ",indcs[ifil][0], indcs[ifil][1], ifil,myid,qm,sib,tang,bang,len(ctx),Util.infomask(ctx[0], None, True)
# 			for im in xrange(indcs[ifil][0], indcs[ifil][1]):
# 				kim = im-indcs[ifil][0]
# 				dst = sqrt((pcoords[cents][0] - pcoords[kim][0])**2 + (pcoords[cents][1] - pcoords[kim][1])**2)
# 				if(kim < cents):  xl = -dst*bang+sib
# 				else:             xl =  dst*bang+sib
# 				shift_x[im] = xl
# 			# Average shift
# 			sx_sum += shift_x[indcs[ifil][0]+cents]
# 		
# 		
# 		# #print myid,sx_sum,total_nfils
# 		sx_sum = mpi_reduce(sx_sum, 1, MPI_FLOAT, MPI_SUM, main_node, MPI_COMM_WORLD)
# 		if myid == main_node:
# 			sx_sum = float(sx_sum[0])/total_nfils
# 			print_msg("Average shift  %6.2f\n"%(sx_sum))
# 		else:
# 			sx_sum = 0.0
# 		sx_sum = 0.0
# 		sx_sum = bcast_number_to_all(sx_sum, source_node = main_node)
# 		for im in xrange(ldata):
# 			shift_x[im] -= sx_sum
			#print  "   %3d  shift_x=%6.3f rotang=%f"%(im,shift_x[im],shift_ang[im]*resamp_dang*180/pi)

		

	##  --------------------------  ##
	##  step 2. - fitting a snake.  ##
	##  --------------------------- ##
	#use straight line as the initial guess to refine the snake.
	paramsline = shift_x

	## for polynomial
	pordera = 2
	porderb = pordera-1
	#a0=[[0.0]*(pordera+1)]*nfils
	#a=[[0.0]*(pordera+1)]*nfils
	#b0=[[0.0]*(porderb+1)]*nfils
	#b=[[0.0]*(porderb+1)]*nfils
	
	
	#for b-spline
	nknots = 9
	a0=[]
	a=[]
	b0=[]
	b=[]
	
	for ifil in xrange(nfils):
		nsegs =  indcs[ifil][1]-indcs[ifil][0]
		tttt=[0.0]*nsegs
		nperiod = nsegs//2
		for j in xrange(-nperiod,0):
			a = 300.0/(nperiod*nperiod*nperiod) 
			import random
			shix=random.randint(-15, 15)
			tttt[j+nperiod]= a * j*j*j  + shix
		for j in xrange(1, nperiod):
			import random
			shix=random.randint(-15, 15)
			tttt[j+nperiod] = 300.0/nperiod*j+shix
			
		out_file = open("initcubic.txt", "w")
		
		for i in xrange(nsegs):
			out_file.write( "%f\n" % (tttt[i]) )
		out_file.close()
		
		#a0[ifil],b0[ifil] = interpoLinecoeffs([paramsline[indcs[ifil][0]], 0], [paramsline[indcs[ifil][0]+1],0], pordera, porderb, nsegs)
		
		#for b-spline
		at=[0.0]*nknots
		at = tttt
		#at=paramsline[indcs[ifil][0]:indcs[ifil][1]]
		#Util.convertTocubicbsplineCoeffs(at, nknots, 1.e-10)
		#for i in xrange(46):
		#	at[i]=bspline(at,i-nperiod,46)
		from scipy import interpolate
		iii=0
		U=[0.0]*nknots
		W=[0.0]*nknots
		AT=[0.0]*nknots
		for i in xrange(nknots):
			U[i] = -nperiod+i*(nsegs-1)*1.0/(nknots-1)
			id   = int(U[i]+nperiod)
			idx  = U[i]+nperiod-id
			if id < nsegs-1:
				att = at[id+1]
			else:
				att = 0	
			AT[i]= (1-idx)*at[id]+idx*att
			W[i] = 1.0
		
		u=[i-nperiod for i in xrange(nsegs)] 
		
		tck=interpolate.splrep(U,AT,W, k=3,s=100)
		at=interpolate.splev(u, tck, der=0, ext=0)
		out_file = open("approxbs.txt", "w")
	
		for i in xrange(nsegs):
			out_file.write( "%f\n" % (at[i]) )
		out_file.close()
	
		# a0.append(at)
# 		a.append(at)
# 		b0.append(at)
# 		b.append(at)
		#a = a0
		#b = b0
			
	## refine using amoeba
	from utilities import amoeba
	for Iter in xrange(max_iter):
		for ifil in xrange(nfils):
			sccf=CCF2d[ifil]
			params0 = a0[ifil]+b0[ifil]   #[paramsline[indcs[ifil][0]:indcs[ifil][1]]
			nsegs =  indcs[ifil][1]-indcs[ifil][0]
			pordera = nknots-1
			nsegsc = nsegs//2
			params =  a[ifil]+b[ifil] 
			fval0 = snakehelicalali(params,[sccf,params0, pordera, porderb, 0.0,2*kx+1,resamp_maxrin])
			#print params0
			im = 1
			print "Iter=%d nsegs=%d"%(Iter, nsegs)
			print "before amoeba im=%d shift_x=%f  rotang=%f "%(im,-shift_x[im], shift_ang[im]*resamp_dang*180/pi)
			
			im = 45
			print "before amoeba im=%d shift_x=%f  rotang=%f fval=%f "%(im,-shift_x[im], shift_ang[im]*resamp_dang*180/pi, fval0)

			## work for parabola filament. pordera=2
			# ftol = 1.e-6
# 			xtol = 1.e-6
# 			maxi = 500
# 			scale = [8.0]*(pordera+1)+[4.0]*(porderb+1)
			
			# ## work for cubic filament. pordera=3			
# 			ftol = 1.e-16
# 			xtol = 1.e-16
# 			maxi = 500
# 			scale = [0.001, 0.001, 0.001, 500.0, 0.001,0.001,50.0]
			#newparams,fval, numit=amoeba(params, scale, snakehelicalali, ftol, xtol, maxi, [sccf,params0, pordera, porderb, 0.0,2*kx+1,resamp_maxrin])		
			## for b-spline
			ftol = 1.e-8
			xtol = 1.e-8

			maxi = 200
			scale=[5.0]*(nknots//2)+[5.0]+[-5.0]*(nknots//2)+[1.0]*nknots		#+[0.001,0.001,50.0] #[4.0]*(porderb+1)

			
			newparams,fval, numit=amoeba(params, scale, snakehelicalali, ftol, xtol, maxi, [sccf,params0, pordera, porderb, 0.0,2*kx+1,resamp_maxrin])
			
			a[ifil] = newparams[0:pordera+1]
			b[ifil] = newparams[pordera+1:pordera+1+porderb+1]
			for iseg in xrange(indcs[ifil][0], indcs[ifil][1]):
				point=bspline(a[ifil],iseg-indcs[ifil][0]-nsegsc, nsegs)		#parabolaf(a[ifil],b[ifil],(iseg-indcs[ifil][0])*1.0/nsegs-nsegsc*1.0/nsegs, pordera, porderb)
				shift_x[iseg] = point
				shift_ang[iseg] = bsplinedu(b[ifil],iseg-indcs[ifil][0]-nsegsc, nsegs) #parabolafb(b[ifil],iseg-indcs[ifil][0], porderb)			#   point[1]
				
				if shift_ang[iseg] < 0.0:
					shift_ang[iseg] = resamp_maxrin - 1.0 + shift_ang[iseg]
					
			im = 1
			print "after amoeba i m=%d shift_x=%f  rotang=%f max_it=%d angid=%f"%(im,-shift_x[im], shift_ang[im]*resamp_dang*180/pi, numit, shift_ang[im])
			im = 45
			print "after amoeba i m=%d shift_x=%f  rotang=%f  max_it=%d angid=%f fval=%f"%(im,-shift_x[im], shift_ang[im]*resamp_dang*180/pi, numit, shift_ang[im], fval)
		
			print newparams
				
	# ## compute new average.	
	if myid == main_node:
		start_time = time()
		print_msg("Iteration #%4d\n"%(total_iter))
	total_iter += 1
	########
	avg = model_blank(nx,ny)     #EMData(nx, ny, 1, False)
	for im in xrange(ldata):
		tttt = fft(data[im])
		#tttt.write_image("image%03d.hdf"%im)
		tttt = rot_shift2D(tttt, shift_ang[im]*resamp_dang*180/pi, -shift_x[im], 0, 0, 1)
		Util.add_img(avg, tttt)
		#tttt.write_image("rot_image%03d.hdf"%im)
	reduce_EMData_to_root(avg, myid, main_node)
	if myid == main_node:
		if CTF:  tavg = Util.divn_filter(avg, ctf_2_sum)
		else:
			tavg = model_blank(nx,ny)    
			tavg = Util.mult_scalar(avg, 1.0/float(nima))
	else:
		tavg = model_blank(nx,ny)
	if Fourvar:
		bcast_EMData_to_all(tavg, myid, main_node)
		vav, rvar = varf2d_MPI(myid, data, tavg, mask, "a", CTF)
	if myid == main_node:
		if Fourvar:
			tavg    = fft(Util.divn_img(fft(tavg), vav))
			vav_r	= Util.pack_complex_to_real(vav)
		# normalize and mask tavg in real space
		#tavg = fft(tavg)
		stat = Util.infomask( tavg, mask, False )
		tavg -= stat[0]
		Util.mul_img(tavg, mask)
		#print "Iter=%d"%(Iter+20)
		tavg.write_image("tavgttt.hdf",-1)
		
	if Fourvar:  del vav
	 							
	# combine shifts found with the original parameters
	for im in xrange(ldata):		
		t1 = Transform()
		##import random
		##shix=random.randint(-10, 10)
		##t1.set_params({"type":"2D","tx":shix})
		t1.set_params({"type":"2D","tx":-shift_x[im]})
		t1.set_rotation({"type":"2d", "alpha":shift_ang[im]*resamp_dang*180.0/pi})
		# combine t0 and t1
		tt = t1 #*init_params[im] ##@ming
		data[im].set_attr("xform.align2d", tt)
	# write out headers and STOP, under MPI writing has to be done sequentially
	mpi_barrier(MPI_COMM_WORLD)
	par_str = ["xform.align2d", "ID"]
	if myid == main_node:
		from utilities import file_type
		if(file_type(stack) == "bdb"):
			from utilities import recv_attr_dict_bdb
			recv_attr_dict_bdb(main_node, stack, data, par_str, 0, ldata, nproc)
		else:
			from utilities import recv_attr_dict
			recv_attr_dict(main_node, stack, data, par_str, 0, ldata, nproc)
	else:           send_attr_dict(main_node, data, par_str, 0, ldata)
	if myid == main_node: print_end_msg("snakehelical-shiftali_MPI")				



def snakehelicalali(params,data):
	sccf    = data[0]
	sccfn = len(sccf)
	sccfnc = sccfn//2
	pordera = data[2]
	porderb = data[3]
	a0 = data[1] 	#[0:pordera+1]
	#b0 = data[1][pordera+1:pordera+1+porderb+1]
	
	lambw   = data[4]
	nx     = data[5]
	angnx  = data[6]
	nxc = nx//2
	angnxc = angnx//2
	
	
	a = params[0:pordera+1]
	b = params[pordera+1:pordera+1+porderb+1]
	#print "lambw", lambw
	sx_sum=0.0
	
	for id in xrange(sccfn):
		point=bspline(a,id-sccfnc, sccfn)	#parabolaf(a,b,id*1.0/sccfn-sccfnc*1.0/sccfn, pordera, porderb)
		xl = point+nxc
		ixl = int(xl)
		dxl = xl - ixl
		al = bsplinedu(b,id-sccfnc, sccfn) #parabolafb(b,id*1.0/sccfn-sccfnc*1.0/sccfn, porderb)		# point[1]
		bl = al
		if al < 0.0:
			al = (angnx-1+al)
		ial = int(al)
		dal = al - ial
		if al < 0.0: print "bl=%f al=%f ial=%d "%(bl, al, ial)
		#print "sx_sum, xl, ixl, dxl", sx_sum, xl,ixl,dxl
		if ixl < 0:
			ixl = 0          
			#print "ixl=%d xl=%f params[id]=%f, ial=%d "%(ixl,xl,xparam[id], ial)
			dxl = 0
		if ixl >= nx-1:
			ixl = nx - 2
			dxl = 0
		# if ial >= angnx-1:
# 			ial = angnx-2 	
# 			dal = 0
		if ixl < 0 or ixl+1 > nx -1 or ial < 0 or (ial +1)%angnx > angnx-1:
			print "AAAAA ixl=%d ial=%d"%(ixl, ial)
		sx_sum += (1.0-dxl)*(1.0-dal)*sccf[id].get_value_at(ixl,ial%angnx) + dxl*(1.0-dal)*sccf[id].get_value_at(ixl+1,ial%angnx) + (1.0-dxl)*dal*sccf[id].get_value_at(ixl,(ial+1)%angnx) + dxl*dal*sccf[id].get_value_at(ixl+1,(ial+1)%angnx)  
		#print "ixl=%d maxia=%d ial=%d get_value_at(ixl,ial)=%f  get_value_at(ixl+1,ial)=%f  get_value_at(ixl,ial+1)=%f  get_value_at(ixl+1,ial+1)=%f"%(ixl, angnx, ial, sccf[id].get_value_at(ixl,ial),sccf[id].get_value_at(ixl+1,ial),sccf[id].get_value_at(ixl,ial+1),sccf[id].get_value_at(ixl+1,ial+1))	
	#print "part 1", sx_sum
	
	part2_sum=0	
	# for id in xrange(sccfn):
# 		point0=parabolaf(a0,b0,id*1.0/sccfn-sccfnc*1.0/sccfn, pordera, porderb)
# 		point1=parabolaf(a,b,id*1.0/sccfn-sccfnc*1.0/sccfn, pordera, porderb)
# 		part2_sum += lambw*((point0[0]-point1[0])**2+(point0[1]-point1[1])**2)
				
	sx_sum -= part2_sum 
	
	# part3_sum = 0
# 	for id in xrange(sccfn-1):
# 		part3_sum += abs(angpar[id]-angpar[id+1])  # + abs(xparam[id] - xparam[id+1])
# 		
# 	sx_sum -= 10.0 * part3_sum
	return sx_sum
					
def interpoLinecoeffs(p0,p1, pordera, porderb, m):
	a=[0.0]*(pordera+1)
	b=[0.0]*(porderb+1)
	c = m//2
	#z=0/m-c, z=1/m-c
	a[1] = (p1[0] - p0[0])*m
	a[0] = p0[0] + a[1]*c
	
	b[1] = (p1[1] - p0[1])*m
	b[0] = p0[1] + b[1]*c
	
	return a, b 

def parabolaf(a,b,z, pordera, porderb):
	point=[0.0]*2
	for i in xrange(pordera+1):
		point[0] += a[i]*pow(z,i)
	for i in xrange(porderb+1):	
		point[1] += b[i]*pow(z,i)
	
	return point

def parabolafb(b,z, porderb):
	point=0.0
	for i in xrange(porderb+1):	
		point += b[i]*pow(z,i)
	
	return point
		
					
def bspline(coefs,z,m):
	ncpoints = len(coefs)
	cents = m//2
	h = m*1.0/ncpoints
		
	val = 0.0
	for i in xrange(ncpoints):
		val += coefs[i]*Util.bsplineBase((z+cents)/h-i)   #  (z,len(coefs),i,3,U)

	return val

def bsplinedu(coefs,z,m):
	ncpoints = len(coefs)
	cents = m//2
	h = m*1.0/ncpoints
		
	val = 0.0
	for i in xrange(ncpoints):
		val += coefs[i]*Util.bsplineBasedu((z+cents)/h-i)   #  (z,len(coefs),i,3,U)

	return val


# get the interval index of u belonging to	
def get_k(U, u):
    for i in range(len(U)-1):
        if U[i] <= u and U[i+1] > u:
            return i
            
        if U[len(U)-2] == u:
            return len(U)-2
            
    raise Exception("Bounding interval for knot not found.")

def get_multiplicity(U, u, i):

    multiplicity = 0
    
    for i in range(i+1,len(U)):
        if U[i] == u:
            multiplicity += 1
        else:
            break
            
    for i in range(i-1,0,-1):
        if U[i] == u:
            multiplicity += 1
        else:
            break
        
    return multiplicity

# deBoor algorithm for parametric b-spline/Bezier curve
# _P: control points sequence.
# U: knot sequence.
# u: compute the b-spline at point u
# p: degree of the b-spline. For cubic b-spline, p=3.
def deboor(_P, U, u, p):
    k = get_k(U, u)
    s = get_multiplicity(U, u, k)
    
    if u != U[k]:
        h = p
    else:
        h = p - s
    
    P = []
    for _p in _P:
        P.append([_p])
    
    for r in range(1, h+1):
        for i in range(k-p+r, k-s+1):
            a = (u - U[i]) / (U[i+p-r+1] - U[i])
            
            # Each point is a tuple
            x = (1 - a) * P[i-1][r-1][0] + a * P[i][r-1][0]
            y = (1 - a) * P[i-1][r-1][1] + a * P[i][r-1][1]
            z = (1 - a) * P[i-1][r-1][2] + a * P[i][r-1][2]     ## added.
            
            P[i].append((x,y,z))
            
    return P[k-s][p-s]	
		
if __name__ == "__main__":
	main()
