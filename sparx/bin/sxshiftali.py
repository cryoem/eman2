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
	print "dang=%f resamp_dfang=%f"%(dang,resamp_dang)
	for wIter in xrange(1):
		## ---------------------------------------------- ##
		## step 1: find a straight line for initial guess ##
		## ---------------------------------------------- ##
		# for Iter in xrange(15):
# 			if myid == main_node:
# 				start_time = time()
# 				print_msg("Iteration #%4d\n"%(total_iter))
# 			total_iter += 1
# 			avg = EMData(nx, ny, 1, False)
# 			for im in xrange(ldata):
# 				Util.add_img(avg, fshift(data[im], shift_x[im]))
# 
# 			reduce_EMData_to_root(avg, myid, main_node)
# 
# 			if myid == main_node:
# 				if CTF:  tavg = Util.divn_filter(avg, ctf_2_sum)
# 				else:    tavg = Util.mult_scalar(avg, 1.0/float(nima))
# 			else:
# 				tavg = model_blank(nx,ny)
# 
# 			if Fourvar:
# 				bcast_EMData_to_all(tavg, myid, main_node)
# 				vav, rvar = varf2d_MPI(myid, data, tavg, mask, "a", CTF)
# 
# 			if myid == main_node:
# 				if Fourvar:
# 					tavg    = fft(Util.divn_img(fft(tavg), vav))
# 					vav_r	= Util.pack_complex_to_real(vav)
# 				# normalize and mask tavg in real space
# 				tavg = fft(tavg)
# 				stat = Util.infomask( tavg, mask, False )
# 				tavg -= stat[0]
# 				Util.mul_img(tavg, mask)
# 				tavg.write_image("tavg.hdf",-1)
# 				# For testing purposes: shift tavg to some random place and see if the centering is still correct
# 				#tavg = rot_shift3D(tavg,sx=3,sy=-4)
# 
# 			# if Fourvar:  del vav
# # 			bcast_EMData_to_all(tavg, myid, main_node)
# # 			tavg = fft(tavg)
# 		
# 			nxc = nx//2
# 			sx_sum  = 0
# 			CCF = [] ##added@ming
# 			for ifil in xrange(nfils):
# 				# Calculate 1D ccf between each segment and filament average
#  				nsegms = indcs[ifil][1]-indcs[ifil][0]
#  				#sccf = []      ##added@ming
#  				ctx = [None]*nsegms
#  				pcoords = [None]*nsegms
#  				for im in xrange(indcs[ifil][0], indcs[ifil][1]):
#  					ttavg = get_im("bdb:vdata",im)
#  					stat = Util.infomask( ttavg, mask, False )
#  					ttavg -= stat[0]
#  					ttavg = fft(ttavg)
#  					ctx[im-indcs[ifil][0]] = Util.window(ccf(ttavg, data[im]), nx, 1)
#  					#sccf.append(ctx[im-indcs[ifil][0]])
#  					pcoords[im-indcs[ifil][0]] = data[im].get_attr('ptcl_source_coord')
#  					#ctx[im-indcs[ifil][0]].write_image("ctx.hdf",im-indcs[ifil][0])
#  					#print "  CTX  ",myid,im,Util.infomask(ctx[im-indcs[ifil][0]], None, True)
#  				CCF.append(ctx)	
#  				# search for best x-shift
#  				cents = nsegms//2
#  			
#  				dst = sqrt(max((pcoords[cents][0] - pcoords[0][0])**2 + (pcoords[cents][1] - pcoords[0][1])**2, (pcoords[cents][0] - pcoords[-1][0])**2 + (pcoords[cents][1] - pcoords[-1][1])**2))
#  				maxincline = atan2(ny//2-2-float(search_rng),dst)
#  				kang = int(dst*tan(maxincline)+0.5)
#  				#print  "  settings ",nsegms,cents,dst,search_rng,maxincline,kang
#  			
#  				# ## C code for alignment. @ming
#  				results = [0.0]*3;
#  				results = Util.helixshiftali(ctx, pcoords, nsegms, maxincline, kang, search_rng,nxc)
#  				sib = int(results[0])
#  				bang = results[1]
#  				qm = results[2]
#  				#print qm, sib, bang
#  			
#  				# qm = -1.e23	
#  
#  				#print qm,six,sib,bang
#  				#print " got results   ",indcs[ifil][0], indcs[ifil][1], ifil,myid,qm,sib,tang,bang,len(ctx),Util.infomask(ctx[0], None, True)
#  				for im in xrange(indcs[ifil][0], indcs[ifil][1]):
#  					kim = im-indcs[ifil][0]
#  					dst = sqrt((pcoords[cents][0] - pcoords[kim][0])**2 + (pcoords[cents][1] - pcoords[kim][1])**2)
#  					if(kim < cents):  xl = -dst*bang+sib
#  					else:             xl =  dst*bang+sib
#  					shift_x[im] = xl
# 				# Average shift
# 				sx_sum += shift_x[indcs[ifil][0]+cents]
# 			
# 			
# 			# #print myid,sx_sum,total_nfils
# 			sx_sum = mpi_reduce(sx_sum, 1, MPI_FLOAT, MPI_SUM, main_node, MPI_COMM_WORLD)
# 			if myid == main_node:
# 				sx_sum = float(sx_sum[0])/total_nfils
# 				print_msg("Average shift  %6.2f\n"%(sx_sum))
# 			else:
# 				sx_sum = 0.0
# 			sx_sum = 0.0
# 			sx_sum = bcast_number_to_all(sx_sum, source_node = main_node)
# 			for im in xrange(ldata):
# 				shift_x[im] -= sx_sum
				#print  "   %3d  shift_x=%6.3f rotang=%f"%(im,shift_x[im],shift_ang[im]*resamp_dang*180/pi)

########## exhausting search.  #####################
		# for im in xrange(ldata):
# 			shift_x[im] = -1* shift_x[im]						
# 		## for 2D case alignment
# 		for Iter2 in xrange(max_iter):
# 			if myid == main_node:
# 				start_time = time()
# 				print_msg("Iteration #%4d\n"%(total_iter))
# 			total_iter += 1
# 			
# 			####  compute avg.
# 			avg = EMData(nx, ny, 1, False)
# 			for im in xrange(ldata):
# 				tttt = fft(data[im])
# 				tttt = rot_shift2D(tttt, shift_ang[im]*resamp_dang*180/pi, -shift_x[im], 0, 0, 1)
# 				tttt = fft(tttt)
# 				Util.add_img(avg, tttt)
# 
# 			reduce_EMData_to_root(avg, myid, main_node)
# 
# 			if myid == main_node:
# 				if CTF:  tavg = Util.divn_filter(avg, ctf_2_sum)
# 				else:    tavg = Util.mult_scalar(avg, 1.0/float(nima))
# 			else:
# 				tavg = model_blank(nx,ny)
# 
# 			if Fourvar:
# 				bcast_EMData_to_all(tavg, myid, main_node)
# 				vav, rvar = varf2d_MPI(myid, data, tavg, mask, "a", CTF)
# 
# 			if myid == main_node:
# 				if Fourvar:
# 					tavg    = fft(Util.divn_img(fft(tavg), vav))
# 					vav_r	= Util.pack_complex_to_real(vav)
# 				# normalize and mask tavg in real space
# 				tavg = fft(tavg)
# 				stat = Util.infomask( tavg, mask, False )
# 				tavg -= stat[0]
# 				Util.mul_img(tavg, mask)
# 				tavg.write_image("tavg.hdf",-1)
# 				#tavg.write_image("tavg_%d.hdf"%Iter)
# 				# For testing purposes: shift tavg to some random place and see if the centering is still correct
# 				#tavg = rot_shift3D(tavg,sx=3,sy=-4)
# 				
# 			if Fourvar:  del vav
# 			bcast_EMData_to_all(tavg, myid, main_node)
# 			##### end computing avg.
# 			
# 			##tavg = fft(tavg)
# 
# 			
# 			print "max ring points=%d nxc=%d angnxc=%d resamp_maxrin=%d"%(numr[len(numr)-1],2*kx+1,maxrin,resamp_maxrin)
# 			total_iter = 0
# 			##tavg = fft(tavg)                                       #transform tavg    into real space.
# 			
# 			tttt = EMData(nx, ny, 1, False)		
# 			## compute 2D ccfs.
# 			crefim = Util.Polar2Dm(tavg, cnx, cny, numr, mode)
# 			Util.Frngs(crefim, numr)
# 			#Util.Applyws(crefim, numr, wr)
# 			sx_sum = 0
# 			CCF2d = []
# 			for ifil in xrange(nfils):
# 				# test Calculate 2D ccf between each segment and filament average
# 				nsegms = indcs[ifil][1]-indcs[ifil][0]
# 				ccf2d = [None]*vol
# 				ctx2d = [None]*nsegms
# 				resamp_ccf2d = [None]*resamp_vol
# 				pcoords = [None]*nsegms
# 				for im in xrange(indcs[ifil][0], indcs[ifil][1]):
# 					tttt = fft(data[im])             # transform data[im] into real space.
# 					#tttt.write_image("ttt%d%d%d.hdf"%(myid,ifil,im))
# 					ccf2d = Util.ali2d_ccf_list_snake(tttt, crefim,  xrng, yrng, rstep, mode, numr, cnx, cny, T)
# 					
# 					for i in xrange(2*kx+1):
# 						for j in xrange(resamp_maxrin):
# 							j_old = int(j * resamp_dang/dang + 0.5)
# 							resamp_ccf2d[i*resamp_maxrin+j] = ccf2d[i*maxrin+j_old]
# 					ddd=resamp_ccf2d.index(max(resamp_ccf2d))
# 					iii = ddd//resamp_maxrin
# 					jjj = ddd%resamp_maxrin
# 					shift_x[im] = iii-(2*kx+1)//2
# 					shift_ang[im] = jjj
# 					#print "ifil=%d iii=%d jjj=%f"%(ifil,iii-(2*kx+1)//2,jjj*resamp_dang)
# 					ccf2dimg = model_blank(2*kx+1, resamp_maxrin)  #EMData(2*kx+1, resamp_maxrin, 1, False)			
# 					for i in xrange(2*kx+1):
# 						for j in xrange(resamp_maxrin):
# 							ccf2dimg.set_value_at(i,j,resamp_ccf2d[i*resamp_maxrin+j])
# 					ctx2d[im-indcs[ifil][0]] = 	ccf2dimg
# 					pcoords[im-indcs[ifil][0]] = data[im].get_attr('ptcl_source_coord')	 
# 					
# 				CCF2d.append(ctx2d)
# 				## search for best x-shift and in-plane rotation angle
# 				cents = nsegms//2
# 			
# 				#dst = sqrt(max((pcoords[cents][0] - pcoords[0][0])**2 + (pcoords[cents][1] - pcoords[0][1])**2, (pcoords[cents][0] - pcoords[-1][0])**2 + (pcoords[cents][1] - pcoords[-1][1])**2))
# 				maxinprot = search_ang*pi/180
# 				kang = int(maxinprot/resamp_dang+0.5)
# 				#print  "  settings ",nsegms,cents,dst,search_rng,maxincline,kang
# 			
# 				# ## C code for alignment. @ming
# 				results = [0.0]*3;
# 				results = Util.snakeshiftali(ctx2d, pcoords, nsegms, resamp_dang, kang, kx,(2*kx+1)//2, resamp_maxrin)
# 				sib = int(results[0])
# 				bang = int(results[1])
# 				qm = results[2]
# 				print "maxinprot=%f kang=%d sib=%d bang=%d qm=%f"%(maxinprot, kang, sib,bang,qm)
# 				for im in xrange(indcs[ifil][0], indcs[ifil][1]):
# 					kim = im-indcs[ifil][0]
# 					dst = sqrt((pcoords[cents][0] - pcoords[kim][0])**2 + (pcoords[cents][1] - pcoords[kim][1])**2)
# 					if(kim < cents):  xl = -dst*tan(resamp_dang*bang)+sib
# 					else:             xl =  dst*tan(resamp_dang*bang)+sib
# 					shift_x[im] = xl
# 					shift_ang[im] = bang
				# Average shift
# 				sx_sum += shift_x[indcs[ifil][0]+cents]
# 			sx_sum = mpi_reduce(sx_sum, 1, MPI_FLOAT, MPI_SUM, main_node, MPI_COMM_WORLD)
# 			if myid == main_node:
# 				sx_sum = float(sx_sum[0])/total_nfils
# 				print_msg("Average shift  %6.2f\n"%(sx_sum))
# 			else:
# 				sx_sum = 0.0
# 			sx_sum = 0.0
# 			sx_sum = bcast_number_to_all(sx_sum, source_node = main_node)
# 			for im in xrange(ldata):
# 				shift_x[im] -= sx_sum	
# 								
# 		for im in xrange(ldata):
# 			shift_x[im] = -1* shift_x[im]		
			#exit()
###### end exhausting search. ##############			

		##  --------------------------  ##
		##  step 2. - fitting a snake.  ##
		##  --------------------------- ##
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
		print "dang=%f resamp_dfang=%f"%(dang,resamp_dang)
		ccf2dimg = EMData(2*kx+1, resamp_maxrin, 1, False)	
		print "max ring points=%d"%numr[len(numr)-1]
		#use straight line as the initial guess to refine the snake.
		paramsline = shift_x
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
				tttt = fft(data[im])             # transform data[im] into real space.
				ccf2d = Util.ali2d_ccf_list_snake(tttt, crefim,  xrng, yrng, rstep, mode, numr, cnx, cny, T)
				for i in xrange(2*kx+1):
					for j in xrange(resamp_maxrin):
						j_old = int(j * resamp_dang/dang + 0.5)
						resamp_ccf2d[i*resamp_maxrin+j] = ccf2d[i*maxrin+j_old]
				# aaaa=resamp_ccf2d.index(max(resamp_ccf2d))
# 					iaaa = aaaa/resamp_maxrin
# 					jaaa = aaaa%resamp_maxrin
# 					shift_x[im] = iaaa-(2*kx+1)//2
# 					shift_ang[im] = jaaa
# 					print "im=%d rotang=%d  maxrin=%d shift=%f"%(im, jaaa, resamp_maxrin, -shift_x[im] )	
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
		## refine using amoeba
		from utilities import amoeba
		
		for Iter in xrange(max_iter):
			for ifil in xrange(nfils):
				sccf=CCF2d[ifil]
				xparams=shift_x[indcs[ifil][0]:indcs[ifil][1]]
				angpams=shift_ang[indcs[ifil][0]:indcs[ifil][1]] 
				# for im in xrange (indcs[ifil][0], indcs[ifil][1]):
# 					angpams[im] = (2*pi-2.81/180*pi)/resamp_dang
							##method 2. maximize the sum
				if ifil == 0:
				#	for im in xrange (indcs[ifil][0], indcs[ifil][1]):
					im = 32
					print "before im=%d shift_x=%f  rotang=%f"%(im,shift_x[im], angpams[im]*resamp_dang*180/pi)
				params0 = paramsline[indcs[ifil][0]:indcs[ifil][1]]
				nsegs = len(xparams)
				params = xparams+angpams 
				fval0 = snakehelicalali(params,[sccf,params0, 0.0,2*kx+1,resamp_maxrin])
				#print "len x =%d len ang = %d len param0=%d"%(nsegs, len(angpams), len(params0))
				scale = [1.5]*nsegs+[1.5]*nsegs
				newparams,fval, numit=amoeba(params, scale, snakehelicalali, 1.e-4, 1.e-8, 500, [sccf,params0, 0.0,2*kx+1,resamp_maxrin])
				print "ifil=%d Iter: %d before  amoeba: %f, after amoeba: %f  max_it=%d"%(ifil,Iter,fval0, fval,numit)
				im = 32
				print "after im=%d shift_x=%f  rotang=%f"%(im,shift_x[im], shift_ang[im]*resamp_dang*180/pi)
				newxpar = newparams[0:nsegs]
				for iseg in xrange(nsegs):
					if newxpar[iseg] > search_rng or newxpar[iseg] < -search_rng:
						newxpar[iseg] = xparams[iseg]
				shift_x[indcs[ifil][0]:indcs[ifil][1]] = newxpar
				shift_ang[indcs[ifil][0]:indcs[ifil][1]] = newparams[nsegs:2*nsegs]

			
					
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
			#tavg.write_image("tavg_%d.hdf"%(Iter+20))
			# For testing purposes: shift tavg to some random place and see if the centering is still correct
			#tavg = rot_shift3D(tavg,sx=3,sy=-4)
		if Fourvar:  del vav
		bcast_EMData_to_all(tavg, myid, main_node)
			###########
			#tavg=get_im("tavg97.hdf")
			##tavg = fft(tavg)

		# for im in xrange(ldata):
# 			shift_x[im] = -1* shift_x[im]	
		
## ----------------------------- ##		
## snake method for x-shift case ##
## ----------------------------- ##		
	# paramsline = shift_x
# 	for ifil in xrange(nfils):
# 		params00=paramsline[indcs[ifil][0]:indcs[ifil][1]]
# 		params00=[int(params00[ip0]+nxc) for ip0 in xrange(len(params00))]
# 		from utilities import write_text_file
# 		#print params00
# 		write_text_file([params00, range(len(params00))],"snakep%df%d.txt"%(myid,ifil))
# 	total_iter = 0
# 	for Iter in xrange(max_iter):
# 		##using b-spline to refine.
# 		from utilities import amoeba
# 		for ifil in xrange(nfils):
# 			sccf=CCF[ifil]
# 			params=shift_x[indcs[ifil][0]:indcs[ifil][1]]
# 			##method 1 maximize every segmentation
# 	# 			for im in xrange(indcs[ifil][0], indcs[ifil][1]):
# 	# 				iseg = im - indcs[ifil][0]
# 	# 				params0=params[iseg]
# 	# 				scale=[0.1]
# 	# 				sccfi = sccf[iseg]
# 	# 				fval0=flexhelicalali1([params0],[sccfi,[params0], 0.5,nxc])
# 	# 				#print "myid ifil fval before amoeba", myid,ifil, fval
# 	# 				[newparams],fval, numit=amoeba([params[iseg]], scale, flexhelicalali1, 1.e-8, 1.e-8, 5000, [sccfi,[params0], 1.0,nxc])
# 	# 				if myid == 0 and im == 9: print "before  amoeba: %f, after amoeba: %f"%(fval0, fval)
# 	# 				shift_x[im] = newparams
# 
# 			##method 2. maximize the sum
# 			iiid=51
# 			params0 = paramsline[indcs[ifil][0]:indcs[ifil][1]]
# 			if ifil == 1:
# 				sccfi = sccf[iiid]
# 				param0177 = params0[iiid]
# 				fval0_177=flexhelicalali1([params[iiid]],[sccfi,[param0177], 0.0,nxc])
# 			nsegs = len(params)
# 			scale = [0.5]*nsegs
# 			fval0 = flexhelicalali(params,[sccf,params0, 0.0,nxc]) 
# 			newparams,fval, numit=amoeba(params, scale, flexhelicalali, 1.e-6, 1.e-6, 700, [sccf,params0, 0.1,nxc])
# 			#print "ifil=%d Iter: %d before  amoeba: %f, after amoeba: %f"%(ifil,Iter,fval0, fval)
# 			for iseg in xrange(nsegs):
# 				if newparams[iseg] > search_rng or newparams[iseg] < -search_rng:
# 					#print "newparams", newparams[iseg]
# 					newparams[iseg] = params[iseg]
# 			shift_x[indcs[ifil][0]:indcs[ifil][1]] = newparams
# 			if ifil == 1:
# 				param177 = newparams[iiid]
# 				fval_177=flexhelicalali1([param177],[sccfi,[param0177], 0.0,nxc])
# 				print "Iter=%d before amoeba:%f, shift=%f; after amoeba:%f shift=%f"%(Iter,fval0_177,param0177, fval_177,param177)
# 			
# 		##compute new average.	
# 		if myid == main_node:
# 			start_time = time()
# 			print_msg("Iteration #%4d\n"%(total_iter))
# 		total_iter += 1
# 		avg = EMData(nx, ny, 1, False)
# 		for im in xrange(ldata):
# 			Util.add_img(avg, fshift(data[im], shift_x[im]))
# 		reduce_EMData_to_root(avg, myid, main_node)
# 		if myid == main_node:
# 			if CTF:  tavg = Util.divn_filter(avg, ctf_2_sum)
# 			else:    tavg = Util.mult_scalar(avg, 1.0/float(nima))
# 		else:
# 			tavg = model_blank(nx,ny)
# 		if Fourvar:
# 			bcast_EMData_to_all(tavg, myid, main_node)
# 			vav, rvar = varf2d_MPI(myid, data, tavg, mask, "a", CTF)
# 		if myid == main_node:
# 			if Fourvar:
# 				tavg    = fft(Util.divn_img(fft(tavg), vav))
# 				vav_r	= Util.pack_complex_to_real(vav)
# 			# normalize and mask tavg in real space
# 			tavg = fft(tavg)
# 			stat = Util.infomask( tavg, mask, False )
# 			tavg -= stat[0]
# 			Util.mul_img(tavg, mask)
# 			tavg.write_image("tavg.hdf",Iter+20)
# 			# For testing purposes: shift tavg to some random place and see if the centering is still correct
# 			#tavg = rot_shift3D(tavg,sx=3,sy=-4)
# 		if Fourvar:  del vav
# 		bcast_EMData_to_all(tavg, myid, main_node)
# 		tavg = fft(tavg)
# 
# 		##compute new CCF.
# 		CCF = [] ##added@ming
# 		for ifil in xrange(nfils):
# 			# Calculate 1D ccf between each segment and filament average
# 			nsegms = indcs[ifil][1]-indcs[ifil][0]
# 			#sccf = []      ##added@ming
# 			ctx = [None]*nsegms
# 			for im in xrange(indcs[ifil][0], indcs[ifil][1]):
# 				ctx[im-indcs[ifil][0]] = Util.window(ccf(tavg, data[im]), nx, 1)
# 				# if im == 177:
# # 					valll=[]
# # 					freq=range(nx)
# # 					[valll.append(ctx[im-indcs[ifil][0]].get_value_at(iii)) for iii in xrange(nx)]
# # 					from utilities import write_text_file
# # 					write_text_file([freq,valll],"ccf%2dit%02d.txt"%(im,Iter))
# # 					ttttt=get_im(stack,im)
# # 					data[im].write_image("datafft.hdf")
# # 					datareal=fft(data[im])
# # 					datareal.write_image("datareal.hdf")
# # 					tavg.write_image("tavgfft.hdf")
# # 					ttttt.write_image("%2d.hdf"%im)
# 			CCF.append(ctx)	
# 			
# 		##save the snake
# 	for ifil in xrange(nfils):
# 		params=shift_x[indcs[ifil][0]:indcs[ifil][1]]
# 		params=[int(params[ip0]+nxc) for ip0 in xrange(len(params))]
# 		from utilities import write_text_file
# 		write_text_file([params,range(len(params))],"snakep%df%dit%d.txt"%(myid,ifil,Iter+1 ))
#################  end 1D snake #####################################################################

 							
	# combine shifts found with the original parameters
	for im in xrange(ldata):		
		t1 = Transform()
		##import random
		##shix=random.randint(-10, 10)
		##t1.set_params({"type":"2D","tx":shix})
		t1.set_params({"type":"2D","tx":shift_x[im]})
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


def flexhelicalali1(params,data):
	sccf    = data[0]
	params0 = data[1]
	lambw   = data[2]
	nxc     = data[3]
	#print "lambw", lambw
	sx_sum=0.0
	sccfn=len(params)
	for id in xrange(sccfn):
		xl = params[id]+nxc
		ixl = int(xl)
		dxl = xl - ixl
		#print "sx_sum, xl, ixl, dxl", sx_sum, xl,ixl,dxl
		sx_sum += (1.0-dxl)*sccf.get_value_at(ixl) + dxl*sccf.get_value_at(ixl+1)
	#print "part 1", sx_sum
	part2_sum=0	
	for id in xrange(sccfn):
		part2_sum += lambw*(params0[id]-params[id])**2
	#print "part 2", part2_sum	
	sx_sum -= part2_sum 	
	return sx_sum
		
def flexhelicalali(params,data):
	sccf    = data[0]
	params0 = data[1]
	lambw   = data[2]
	nxc     = data[3]
	#print "lambw", lambw
	sx_sum=0.0
	sccfn=len(params)
	for id in xrange(sccfn):
		xl = params[id]+nxc
		ixl = int(xl)
		dxl = xl - ixl
		#print "sx_sum, xl, ixl, dxl", sx_sum, xl,ixl,dxl
		if ixl < 0:
			print "ixl=%d xl=%f params[id]=%f"%(ixl,xl,params[id])
		sx_sum += (1.0-dxl)*sccf[id].get_value_at(ixl) + dxl*sccf[id].get_value_at(ixl+1)
	#print "part 1", sx_sum
	part2_sum=0	
	for id in xrange(sccfn):
		part2_sum += lambw*(params0[id]-params[id])**2
	#print "part 2", part2_sum	
	sx_sum -= part2_sum 	
	return sx_sum
	
def snakehelicalali(params,data):
	sccf    = data[0]
	xparam0 = data[1]
	lambw   = data[2]
	nx     = data[3]
	angnx  = data[4]
	nxc = nx//2
	angnxc = angnx//2
	sccfn=len(xparam0)
	
	xparam = params[0:sccfn]
	angpar = params[sccfn:2*sccfn]
	#print "lambw", lambw
	sx_sum=0.0
	
	for id in xrange(sccfn):
		xl = xparam[id]+nxc
		ixl = int(xl)
		dxl = xl - ixl
		al = angpar[id]
		if al < 0:
			al = (angnx-1+al)
		ial = int(al)
		dal = al - ial
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
		sx_sum += (1.0-dxl)*(1.0-dal)*sccf[id].get_value_at(ixl,ial) + dxl*(1.0-dal)*sccf[id].get_value_at(ixl+1,ial) + (1.0-dxl)*dal*sccf[id].get_value_at(ixl,ial+1) + dxl*dal*sccf[id].get_value_at(ixl+1,ial+1)  
		#print "ixl=%d maxia=%d ial=%d get_value_at(ixl,ial)=%f  get_value_at(ixl+1,ial)=%f  get_value_at(ixl,ial+1)=%f  get_value_at(ixl+1,ial+1)=%f"%(ixl, angnx, ial, sccf[id].get_value_at(ixl,ial),sccf[id].get_value_at(ixl+1,ial),sccf[id].get_value_at(ixl,ial+1),sccf[id].get_value_at(ixl+1,ial+1))	
	#print "part 1", sx_sum
	part2_sum=0	
	for id in xrange(sccfn):
		part2_sum += lambw*(xparam0[id]-xparam[id])**2 
	#print "part 2", part2_sum	
	sx_sum -= part2_sum 
	return sx_sum
					
def curhelicalshiftali_MPI(stack, maskfile=None, maxit=100, CTF=False, snr=1.0, Fourvar=False, search_rng=-1):
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
	from math         import sqrt, atan2, tan, sin, cos
	
	nproc = mpi_comm_size(MPI_COMM_WORLD)
	myid = mpi_comm_rank(MPI_COMM_WORLD)
	main_node = 0
		
	ftp = file_type(stack)

	if myid == main_node:
		print_begin_msg("curhelical-shiftali_MPI")

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
				ERROR('Incorrectly set qctf flag', "curhelicalshiftali_MPI", 1,myid)
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
		circnum = nx//2         #number of circles.  
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
			
			## C code for alignment. @ming
 			results = [0.0]*6;
 			results = Util.curhelixshiftali(ctx, pcoords, nsegms, search_rng, nx, ny)
 			sib = int(results[0])
 			bang = results[1]
 			CircR = results[2]
 			RSinang = results[3]
 			RCosang = results[4]
 			qm = results[5]
# 			print "in C", sib, bang, CircR, RSinang, RCosang, qm
#  			dst1 = sqrt((pcoords[cents][0] - pcoords[0][0])**2 + (pcoords[cents][1] - pcoords[0][1])**2)
#  			dst2 = sqrt((pcoords[cents][0] - pcoords[-1][0])**2 + (pcoords[cents][1] - pcoords[-1][1])**2)
#  			qm = -1.e23
#  			for circ in xrange(circnum):
#  				circd = 4*(circ+1)           ## circd = n*(circ+1), by adjusting n to obtain different number of circles. 
#  				circR = circd/2.0 + (ny*nsegms)**2/(8.0*circd)
#  				x1 = circR - sqrt(circR**2-dst1**2)
#  				x2 = circR - sqrt(circR**2-dst2**2)
#  				xmax = max(x1,x2)
#  				if ( xmax == x1 ):
#  					dst0 = dst1
#  				else:
#  					dst0 = dst2
#  				for six in xrange(-search_rng, search_rng+1,1):
#  					if (nx-(six+nxc+xmax)-2 <= 0): continue
#  					xmin = min(nx-(six+nxc+xmax)-2, six+nxc-2)
#  					maxincline = atan2(xmin,dst0)
#  					kang = int(dst0*tan(maxincline)+0.5)
#  					q0 = ctx[cents].get_value_at(six+nxc)
#  					#for incline in xrange(kang+1):
#  					for incline in xrange(-kang,kang+1):
#  						qt = q0
#  						#qu = q0
#  						if(kang>0):  
#  							tang = tan(maxincline/kang*incline)
#  							Rsinang = circR*sin(maxincline/kang*incline)
#  							Rcosang = circR*cos(maxincline/kang*incline) 
#  						else:        
#  							tang = 0.0
#  							Rsinang = 0.0
#  							Rcosang = 1.0
#  						for kim in xrange(cents+1,nsegms):
#  							dst = sqrt((pcoords[cents][0] - pcoords[kim][0])**2 + (pcoords[cents][1] - pcoords[kim][1])**2)
#  							if ( tang >= 0 ):                                   ## case 1 and 4.								
#  								#xl = dst*tang+six+nxc
#  								xl = Rcosang - sqrt(Rcosang**2-2*dst*Rsinang-dst**2)+six+nxc
#  								ixl = int(xl)
#  								dxl = xl - ixl
#  								#print "  A  ", ifil,six,incline,kim,xl,ixl,dxl
#  								qt += (1.0-dxl)*ctx[kim].get_value_at(ixl) + dxl*ctx[kim].get_value_at(ixl+1)
#  							if ( tang < 0 ):
#  								if ( -dst*tang > circR - sqrt(circR**2-dst**2)):	             ## case 3
#  									#xl = -dst*tang+six+nxc
#  									if ( Rcosang**2+2*dst*Rsinang-dst**2 <0 ):
#  										ERROR('Rcosang**2+2*dst*Rsinang-dst**2 less than 0', "curhelicalshiftali_MPI", 1,myid) 
#  									xl = -Rcosang + sqrt(Rcosang**2+2*dst*Rsinang-dst**2)+six+nxc
#  									#print "incline=%f Rcos=%f Rsin=%f dst=%f ny=%f nxc=%f circd=%f circR=%f kim=%f cents+1=%d xl=%f six=%d"%(incline, Rcosang, Rsinang, dst, ny, nxc, circd, circR, kim, cents+1, xl,six)
#  									ixl = int(xl)
#  									dxl = xl - ixl
#  									if ( ixl < 0 ):
#  										print "kim=%d, cents=%d nsegms=%d incline=%d ixl=%f"%(kim, cents,nsegms, incline,ixl) 
#  									qt += (1.0-dxl)*ctx[kim].get_value_at(ixl) + dxl*ctx[kim].get_value_at(ixl+1)
#  								else:                                           ## case 2
#  									xl = Rcosang - sqrt(Rcosang**2+2*dst*Rsinang-dst**2)+six+nxc
#  									ixl = int(xl)
#  									dxl = xl - ixl
#  									qt += (1.0-dxl)*ctx[kim].get_value_at(ixl) + dxl*ctx[kim].get_value_at(ixl+1)										
#  						for kim in xrange(cents):
#  							dst = sqrt((pcoords[cents][0] - pcoords[kim][0])**2 + (pcoords[cents][1] - pcoords[kim][1])**2)
#  							if ( tang >= 0 ):  
#  								if ( dst*tang > circR - sqrt(circR**2-dst**2)):	             
#  									#xl = -dst*tang+six+nxc
# 									xl = -Rcosang + sqrt(Rcosang**2+2*dst*Rsinang-dst**2)+six+nxc
# 									#print "incline=%f Rcos=%f Rsin=%f dst=%f ny=%f nxc=%f circd=%f circR=%f kim=%f cents+1=%d xl=%f six=%d"%(incline, Rcosang, Rsinang, dst, ny, nxc, circd, circR, kim, cents+1, xl,six)
#  									ixl = int(xl)
#  									dxl = xl - ixl
#  									qt += (1.0-dxl)*ctx[kim].get_value_at(ixl) + dxl*ctx[kim].get_value_at(ixl+1)
#  								else:                                           
#  									xl = Rcosang - sqrt(Rcosang**2+2*dst*Rsinang-dst**2)+six+nxc
#  									ixl = int(xl)
#  									dxl = xl - ixl
#  									qt += (1.0-dxl)*ctx[kim].get_value_at(ixl) + dxl*ctx[kim].get_value_at(ixl+1)	
#  							if ( tang < 0 ):
#  								xl = Rcosang - sqrt(Rcosang**2-2*dst*Rsinang-dst**2)+six+nxc
#  								ixl = int(xl)
#  								dxl = xl - ixl
#  								#print "  A  ", ifil,six,incline,kim,xl,ixl,dxl
#  								qt += (1.0-dxl)*ctx[kim].get_value_at(ixl) + dxl*ctx[kim].get_value_at(ixl+1)															
#  						if( qt > qm ):
#  							qm = qt
#  							sib = six
#  							bang = tang
#  							CircR = circR
#  							RSinang = Rsinang
#  							RCosang = Rcosang
#  						#if incline == 0:  print  "incline = 0  ",six,tang,qt,qu
#  					#print qm,six,sib,bang
#  				#print " got results   ",indcs[ifil][0], indcs[ifil][1], ifil,myid,qm,sib,tang,bang,len(ctx),Util.infomask(ctx[0], None, True)
#  			print "in python", sib, bang, CircR, RSinang, RCosang, qm
			for im in xrange(indcs[ifil][0], indcs[ifil][1]):
				kim = im-indcs[ifil][0]
				dst = sqrt((pcoords[cents][0] - pcoords[kim][0])**2 + (pcoords[cents][1] - pcoords[kim][1])**2)
				if(kim < cents):  
					if ( bang >= 0 ):
						if ( dst*bang > CircR - sqrt(CircR**2-dst**2)):
							xl = -RCosang + sqrt(RCosang**2+2*dst*RSinang-dst**2)+sib
						else:
							xl = RCosang - sqrt(RCosang**2+2*dst*RSinang-dst**2)+sib
					if ( bang < 0 ):
							xl = RCosang - sqrt(RCosang**2-2*dst*RSinang-dst**2)+sib
				else:
					if ( bang >= 0 ):					
						xl = RCosang - sqrt(RCosang**2-2*dst*RSinang-dst**2)+sib
					if ( bang < 0 ):
						if ( -dst*bang > CircR - sqrt(CircR**2-dst**2)):	             ## case 3
							xl = -RCosang + sqrt(RCosang**2+2*dst*RSinang-dst**2)+sib
						else:
							xl = RCosang - sqrt(RCosang**2+2*dst*RSinang-dst**2)+sib				 
				shift_x[im] = xl
			# Average shift
			sx_sum += shift_x[indcs[ifil][0]+cents]
		#print myid,sx_sum,total_nfils
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
	if myid == main_node: print_end_msg("curhelical-shiftali_MPI")				

					


if __name__ == "__main__":
	main()
