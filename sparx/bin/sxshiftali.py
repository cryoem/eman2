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
	usage = progname + " stack outdir <maskfile> --search_rng --ou=outer_radius --maxit=max_iteration --CTF --snr=SNR --Fourvar=Fourier_variance --oneDx --MPI"
	parser = OptionParser(usage,version=SPARXVERSION)
	parser.add_option("--search_rng",       type="int",  default=-1,             help="Used to compute the dimension of a \nnwx by nwx section of the 2D ccf which is \nwindowed out for peak search: \nnwx=2*search_rng+1 (nwx=nx if search_rng is -1))")
	parser.add_option("--ou",       type="float",  default=-1,            help="radius of the particle - used for constructing the default mask. If ou is -1, then the mask is a circle with radius nx/2 - 2")
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
			
	if len(args) < 2 or len(args) > 3:
    		print "usage: " + usage
    		print "Please run '" + progname + " -h' for detailed options"
	else:
		if args[1] == 'None': outdir = None
		else:		      outdir = args[1]

		if len(args) == 2: mask = None
		else:              mask = args[2]
			
		if global_def.CACHE_DISABLE:
			from utilities import disable_bdb_cache
			disable_bdb_cache()
		
		from mpi import mpi_init
		sys.argv = mpi_init(len(sys.argv),sys.argv)

		global_def.BATCH = True
		shiftali_MPI(args[0], outdir, mask,options.ou, options.maxit, options.CTF, options.snr, options.Fourvar,options.search_rng,options.oneDx)
		global_def.BATCH = False
		
		from mpi import mpi_finalize
		mpi_finalize()

def shiftali_MPI(stack, outdir, maskfile=None, ou=-1, maxit=100, CTF=False, snr=1.0, Fourvar=False, search_rng=-1, oneDx=False):  
	from applications import MPI_start_end
	from utilities    import model_circle, model_blank, get_image, peak_search
	from utilities    import reduce_EMData_to_root, bcast_EMData_to_all, send_attr_dict, file_type, bcast_number_to_all, bcast_list_to_all
	from statistics   import varf2d_MPI
	from fundamentals import fft, ccf, rot_shift3D, rot_shift2D
	from utilities    import get_params2D, set_params2D
	from utilities    import print_msg, print_begin_msg, print_end_msg
	import os
	import sys
	from mpi 	  import mpi_init, mpi_comm_size, mpi_comm_rank, MPI_COMM_WORLD
	from mpi 	  import mpi_reduce, mpi_bcast, mpi_barrier, mpi_gatherv
	from mpi 	  import MPI_SUM, MPI_FLOAT, MPI_INT
	from EMAN2	  import Processor
	
	number_of_proc = mpi_comm_size(MPI_COMM_WORLD)
	myid = mpi_comm_rank(MPI_COMM_WORLD)
	main_node = 0
		
	ftp = file_type(stack)
	if outdir:
		if os.path.exists(outdir):  ERROR('Output directory exists, please change the name and restart the program', "shftali_MPI", 1, myid)
		mpi_barrier(MPI_COMM_WORLD)

	if myid == main_node:
		print_begin_msg("shiftali_MPI")
		if outdir:	os.mkdir(outdir)

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
		if CTF:	ctf_app = ima.get_attr_default('ctf_applied', 2)
		del ima
	else:
		nx = 0
		if CTF:	ctf_app = 0
	nx = bcast_number_to_all(nx, source_node = main_node)
	if CTF:
		ctf_app = bcast_number_to_all(ctf_app, source_node = main_node)
		if ctf_app > 0:	ERROR("data cannot be ctf-applied", "shiftali_MPI", 1, myid)

	mask = model_circle(nx//2-2, nx, nx)

	cnx  = nx/2+1
 	cny  = cnx

	if CTF:
		from filter import filt_ctf
		from morphology   import ctf_img
		ctf_abs_sum = EMData(nx, nx, 1, False)
		ctf_2_sum = EMData(nx, nx, 1, False)
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
			ctfimg = ctf_img(nx, ctf_params)
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
		t = data[im].get_attr(data[im], 'xform.align2d')
		init_params.append(t)
		!!alpha, sx, sy, mirror, scale =  t.get_params({"type":"2d","alpha":alpha,"scale":scale,"mirror":mirror,"tx":sx,"ty":sy})  #  PLEASE CHECK HOW TO DO IT
		data[im] = rot_shift2D(data[im], alpha, sx, sy, mirror)		
	
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
		total_iter += 1
		avg = EMData(nx, nx, 1, False)
		for im in data:  Util.add_img(avg, im)

		reduce_EMData_to_root(avg, myid, main_node)

		if myid == main_node:
			print_msg("Iteration #%4d\n"%(total_iter))
			if CTF:
				tavg = Util.divn_filter(avg, ctf_2_sum)
			else:	 tavg = Util.mult_scalar(avg, 1.0/float(nima))
		else:
			tavg =  fft(model_blank(nx, nx))                                  #  shouldn't it read: tavg = EMData(nx, nx, 1, False) ???
			cs = [0.0]*2

		if Fourvar:
			bcast_EMData_to_all(tavg, myid, main_node)
			vav, rvar = varf2d_MPI(myid, data, tavg, mask, "a", CTF)

		if myid == main_node:
			if Fourvar:
				tavg    = fft(Util.divn_img(fft(tavg), vav))
				vav_r	= Util.pack_complex_to_real(vav)
				if outdir:
					vav_r.write_image(os.path.join(outdir, "varf.hdf"), total_iter-1)

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

		for im in xrange(len(data)):
			if oneDx:
				ctx = Util.window(ccf(data[im],tavg),nwx,1)
				p1  = peak_search(ctx)
				p1_x = -int(p1[0][3])
				ishift_x[im] = p1_x
				sx_sum += p1_x

			else:
				p1 = peak_search(Util.window(ccf(data[im],tavg), nwx,nwx))
				p1_x = -int(p1[0][4])
				p1_y = -int(p1[0][5])
				ishift_x[im] = p1_x
				ishift_y[im] = p1_y
				sx_sum += p1_x
				sy_sum += p1_y

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

	for im in xrange(len(data)):
		data[im] = fft(data[im])

	scale = ((init_params[0]).get_params("2D"))["scale"]        #  I DO NOT UNDERSTAND THIS LINE, LOOKS FUNKY.  INIT_PARAMS CONTAINS TRANFORM OBJECTS, SO THERE IS A PROPER WAY TO RETRIEVE PARAMS
	for im in xrange(len(data)):		
		t0 = init_params[im]
		t1 = Transform()
		t1.set_params({"type":"2D","alpha":0,"scale":scale,"mirror":0,"tx":shift_x[im],"ty":shift_y[im]})
		# combine t0 and t1
		tt = t1*t0
		data[im].set_attr("xform.align2d", tt)   #  I CHANGED HERE,  TRANS OBJ SHOULD BE SET< READ DIRECTLY

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

if __name__ == "__main__":
	main()
