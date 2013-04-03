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
	parser.add_option("--search_rng",       type="int",  default=-1,             help="Used to compute the dimension of a \nnwx by nwy section of the 2D ccf which is \nwindowed out for peak search: \nnwx=2*search_rng+1 (nwx=nx if search_rng is -1))")
	parser.add_option("--search_rng_y",       type="int",  default=-1,             help="Used to compute the dimension of a \nnwx by nwy section of the 2D ccf which is \nwindowed out for peak search: \nnwy=2*search_rng_y+1 (nwy=ny if search_rng_y is -1)). This is ignored when oneDx flag is activated for x-shift search.")
	parser.add_option("--maxit",    type="float",  default=100,             help="maximum number of iterations program will perform")
	parser.add_option("--CTF",      action="store_true", default=False,   help="use CTF correction during centering ")
	parser.add_option("--snr",      type="float",  default=1.0,           help="signal-to-noise ratio of the data (default is 1.0)")
	parser.add_option("--Fourvar",  action="store_true", default=False,   help="compute Fourier variance")
	parser.add_option("--oneDx",  action="store_true", default=False,   help="Window out central line of 2D cross correlation for peak search")
	parser.add_option("--MPI",      action="store_true", default=False,   help="use multiple processors ")
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

	total_iter = 0
	cs = [0.0]*2


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
			cs = [0.0]*2

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
		else:              nwx=nx
		
		if search_rng_y > 0: nwy = 2*search_rng_y+1
		else:              nwy=ny
		
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
			
		if not_zero_all == 0:
			break
		
	for im in xrange(len(data)):
		data[im] = fft(data[im])

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

def helishiftali_MPI(stack, maskfile=None, maxit=100, CTF=False, snr=1.0, Fourvar=False, search_rng=-1):  
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
	from pixel_error  import ordersegments
	from development  import chunks_distribution
	
	nproc = mpi_comm_size(MPI_COMM_WORLD)
	myid = mpi_comm_rank(MPI_COMM_WORLD)
	main_node = 0
		
	ftp = file_type(stack)

	if myid == main_node:
		print_begin_msg("shiftali_MPI")

	max_iter=int(maxit)

		if( myid == 0):
		infils = EMUtil.get_all_attributes(stack, "filament")
		ptlcoords = EMUtil.get_all_attributes(stack, 'ptcl_source_coord')
		filaments = ordersegments(infils, ptlcoords)
		total_nfils = len(filaments)
		inidl = [0]*total_nfils
		for i in xrange(total_nfils):  inidl[i] = len(filaments[i])
		linidl = sum(inidl)
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
	nima = len(data)
	nx = data[0].get_xsize()
	ny = data[0].get_ysize()
	if maskfile == None:
		mrad = min(nx, ny)
		mask = model_circle(mrad//2-2, nx, ny)
	else:
		mask = get_im(maskfile)

	for im in xrange(nima):
		data[im].set_attr('ID', list_of_particles[im])
		if CTF:
			ctf_params = data[im].get_attr("ctf")
			st = Util.infomask(data[im], mask, False)
			data[im] -= st[0]
			data[im] = filt_ctf(data[im], ctf_params)
			data[im].set_attr('ctf_applied', 1)

	del list_of_particles		

	total_iter = 0
	cs = [0.0]*2

	# fourier transform all images, and apply ctf if CTF
	for im in xrange(len(data)):
		if CTF:
			ctf_params = data[im].get_attr("ctf")
			data[im] = filt_ctf(fft(data[im]), ctf_params)
		else:
			data[im] = fft(data[im])
	
	for Iter in xrange(max_iter):
	
		# Calculate filament average
		for ifil in xrange(nfils):
			avg = EMData(nx, ny, 1, False)
			filnima = 0
			for im in xrange(indcs[ifil][0], indcs[ifil][1]):
				Util.add_img(avg, data[im])
				filnima += 1
			tavg = Util.mult_scalar(avg, 1.0/float(filnima))
			
			# Calculate 1D ccf between each segment and filament average
			for im in xrange(indcs[ifil][0], indcs[ifil][1]):
				ctx = Util.window(ccf(data[im],tavg),nwx,1)

if __name__ == "__main__":
	main()
