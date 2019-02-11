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

#  This file contains code under development or not currently used.


'''  From fundamentals  '''

def select_stack(stack, select_stack):
	nima = EMUtil.get_image_count(stack)
	ima = EMData()
	ima.read_image(stack,0)
	nclass = ima.get_attr('nclass')
	class2 = []
	for im in xrange(nima):
		ima.read_image(stack,im)
		n = ima.get_attr('assign')
		if n not in class2:
			class2.append(n)		
			ima.write_image(select_stack, n)

'''
def rtsh(image, angle, sx=0., sy=0.):
	m=image.get_xsize()
	# padd two times
	npad=2
	n=m*npad
	# support of the window
	k=6
	alpha=1.75
	kb = Util.KaiserBessel(alpha, k, m/2, k/(2.*n), n)
	return grotshift2d(image, angle, kb, npad, sx, sy)
'''
def list_send_im(list, dst):
	from utilities import send_EMData, model_blank
	N    = len(list)
	data = model_blank(N)
	for n in xrange(N): data.set_value_at(n, float(list[n]))
	send_EMData(data, dst, 0)

def list_recv_im(org):
	from utilities import recv_EMData
	data = recv_EMData(org, 0)
	N    = data.get_xsize()
	list = [0] * N
	for n in xrange(N): list[n] = int(data.get_value_at(n))
	return list

# list must be allocated
def list_bcast_im(list, myid, main_node):
	from utilities import bcast_EMData_to_all, model_blank
	N   = len(list)
	data = model_blank(N)
	for n in xrange(N): data.set_value_at(n, float(list[n]))
	bcast_EMData_to_all(data, myid, main_node)
	for n in xrange(N): list[n] = int(data.get_value_at(n))
	return list


###############################################################################################
# draft ali2d_rac_MPI
###############################################################################################

def ali2d_rac_MPI(stack, maskfile = None, kmeans = 'None', ir = 1, ou = -1, rs = 1, nclass = 2, maxit = 10, maxin = 10, check_mirror = False, rand_seed = 10):
	from global_def import MPI
	from utilities  import bcast_EMData_to_all, reduce_EMData_to_root, bcast_number_to_all
	from utilities  import send_EMData, recv_EMData
	from utilities  import model_circle, combine_params2
	from utilities  import get_arb_params, set_arb_params
	from statistics import MPIlogfile_init, MPIlogfile_print, MPIlogfile_end
	from statistics import kmnr, kmn
	from random     import randint, seed, shuffle
	from alignment  import Numrinit, ringwe, ang_n
	from copy       import deepcopy
	import time
	import sys
	from mpi 	  import mpi_comm_size, mpi_comm_rank, MPI_COMM_WORLD
	from mpi 	  import mpi_reduce, mpi_bcast, mpi_barrier, mpi_recv, mpi_send
	from mpi 	  import MPI_SUM, MPI_FLOAT, MPI_INT, MPI_LOR

	# [id]   part of code different for each node
	# [sync] synchronise each node
	# [main] part of code just for the main node
	# [all]  code write for all node

	# debug mode
	DEBUG = False
	K_th  = 1

	# To work
	if kmeans == 'None':
		start_k = False
	else:
		start_k = True
	
	# ------------------- MPI init ----------------------------- #
	# init
	number_of_proc = mpi_comm_size(MPI_COMM_WORLD)
	myid           = mpi_comm_rank(MPI_COMM_WORLD)

	# chose a random node as a main one
	main_node = 0
	if myid  == 0:	main_node = randint(0,number_of_proc-1)
	main_node = bcast_number_to_all(main_node,0)
	mpi_barrier(MPI_COMM_WORLD)

	if number_of_proc > nclass:
		if myid == main_node:
			print 'need number of cpus > K'
		return	
		
	t_start = time.time()
	
	# ---------------------------------------------------------- #
	
	seed(rand_seed)
	first_ring=int(ir); last_ring=int(ou); rstep=int(rs); K=int(nclass); max_iter=int(maxit); max_internal=int(maxin);
	nima = EMUtil.get_image_count(stack)
	temp = EMData()
	temp.read_image(stack, 0)
	nx = temp.get_xsize()
	ny = nx

	# default value for the last ring
	if(last_ring==-1): last_ring=nx//2-2

	# somebody dedicated could write a version with option "H" for half rings that would work for ACF functions.
	mode = "F"

	# [all] data structure
	grp_asg    = [] # assignment of groups to node
	GRP        = [] # list of local groups, each group contain index images
	IM         = [] # list of local groups, each group contain images
	AVE        = [] # list of averages (for all groups)

	# [main] 
	if myid == main_node:
		# define assignment of groups to node
		grp_asg = []
		for k in xrange(K): grp_asg.append(k % 	number_of_proc)
		grp_asg.sort()

		if start_k == False:
			# randomly list of index
			tmp = range(0, nima)
			shuffle(tmp)

			list_i = []
			for n in xrange(nima): list_i.append(n % K)
			list_i.sort()

			list_index = [[] for n in xrange(K)]
			n = 0
			for i in list_i:
			    list_index[i].append(tmp[n])
			    n += 1
		else:
			# pick-up index from ave of k-means
			list_index = [[] for n in xrange(K)]
			im = EMData()
			for k in xrange(K):
				im.read_image(kmeans, k, True)
				listim = im.get_attr('members')
				for index in listim: list_index[k].append(int(index))
				if k == 0: print list_index

	## to test
	#return

	# [main] broadcast list of assignment of groups to node
	grp_asg = mpi_bcast(grp_asg, K, MPI_INT, main_node, MPI_COMM_WORLD)
	grp_asg = grp_asg.tolist()

	# [all] display infos
	MPIlogfile_init(myid)
	MPIlogfile_print(myid, '************* ali2d_rac MPI *************\n')
	MPIlogfile_print(myid, 'Input stack                          : %s\n' % stack)
	MPIlogfile_print(myid, 'Mask file                            : %s\n' % maskfile)
	MPIlogfile_print(myid, 'Inner radius                         : %i\n' % first_ring)
	MPIlogfile_print(myid, 'Outer radius                         : %i\n' % last_ring)
	MPIlogfile_print(myid, 'Ring step                            : %i\n' % rstep)
	MPIlogfile_print(myid, 'Maximum iteration                    : %i\n' % max_iter)
	MPIlogfile_print(myid, 'Maximum internal iteration           : %i\n' % max_internal)
	MPIlogfile_print(myid, 'Consider mirror                      : %s\n' % check_mirror)
	MPIlogfile_print(myid, 'Random seed                          : %i\n' % rand_seed)
	MPIlogfile_print(myid, 'Number of classes                    : %i\n' % K)
	MPIlogfile_print(myid, 'Number of images                     : %i\n' % nima)
	MPIlogfile_print(myid, 'Number of cpus                       : %i\n' % number_of_proc)
	MPIlogfile_print(myid, 'ID cpu                               : %i\n' % myid)

	nb_grp = 0
	for k in xrange(K):
		if grp_asg[k] == myid: nb_grp += 1

	MPIlogfile_print(myid, 'Number of classes in this cpu        : %i\n' % nb_grp)
	MPIlogfile_print(myid, 'Output file                          : %s\n\n' % stack)

	# [all] generate LUT group, ex global indice grp 4 5 6 -> local indice in the node 0 1 2
	# [4, 5, 6] -> [-1, -1, -1, -1, 0, 1, 2]
	lut_grp = [-1] * K
	i       = 0
	for k in xrange(K):
		if myid == grp_asg[k]:
			lut_grp[k] = i
			i += 1
	
	# [main] send list_index to each individual group
	for k in xrange(K):
		size = []
		dst  = -1
		tmp  = []

		if myid == main_node:
			if grp_asg[k] != main_node:
				size = len(list_index[k])
				mpi_send(size, 1, MPI_INT, grp_asg[k], 0, MPI_COMM_WORLD)
				mpi_send(list_index[k], size, MPI_INT, grp_asg[k], 0, MPI_COMM_WORLD)
			else:
				GRP.append(list_index[k])
		else:
			if myid == grp_asg[k]:
				size = mpi_recv(1, MPI_INT, main_node, 0, MPI_COMM_WORLD)
				tmp  = mpi_recv(size, MPI_INT, main_node, 0, MPI_COMM_WORLD)
				tmp  = tmp.tolist()
				GRP.append(tmp)

		mpi_barrier(MPI_COMM_WORLD)

	# [all] precalculate rings
	numr     = Numrinit(first_ring, last_ring, rstep, mode)
	wr       = ringwe(numr ,mode)
	lnumr    = numr[len(numr)-1]
	norm_rsd = 0
	for n in xrange(1, len(numr), 3): norm_rsd += numr[n]
	
	# [all] prepare 2-D mask for normalization
	if maskfile:
		import  types
		if(type(maskfile) is types.StringType):
			mask2D = get_image(maskfile)
		else:
			mask2D = maskfile
	else:
		mask2D = model_circle(last_ring,nx,nx)
	if(first_ring > 0):
		tave = model_circle(first_ring-1, nx, ny)
		mask2D -= tave

	#  center is in SPIDER convention
	cnx = int(nx/2) + 1
	cny = int(ny/2) + 1
	
	# [id] read images in each local groups and resample them into polar coordinates
	if DEBUG: MPIlogfile_print(myid, 'prepare images\n')
	for grp in xrange(nb_grp):
		data = []
		im   = 0
		for index in GRP[grp]:
			temp = EMData()
			temp.read_image(stack, index)
			sx = temp.get_attr('sx')
			sy = temp.get_attr('sy')
			alpha_original = temp.get_attr('alpha')
			miri = temp.get_attr('mirror')
			[mean, sigma, qn, qm] = Util.infomask(temp, mask2D, True)
			temp = (temp - mean)/sigma
			alpha_original_n,sxn,syn,mir = combine_params2(0, -sx, -sy, 0, -alpha_original,0,0,0)
			cimage = Util.Polar2Dm(temp, cnx+sxn, cny+syn, numr, mode)
			Util.Frngs(cimage, numr)
			data.append(cimage)
			data[im].set_attr_dict({'alpha':1.0, 'alpha_original':alpha_original, 'sx':sx, 'sy':sy, 'mirror': 0})

			im += 1

		IM.append(deepcopy(data))

	del temp
	del data
	del mask2D

	# use image to 0 for define a blank image (use to reset the average)
	blank = IM[0][0]
	blank.to_zero()

	# prepare the syncronize sequence (for the com between the nodes)
	sync_org = []
	sync_dst = []
	for i in xrange(number_of_proc):
		tmp = range(number_of_proc)
		sync_dst.extend(tmp)
		for j in xrange(number_of_proc):
			sync_org.append(i)

	# debug
	if DEBUG:
		flow_ctrl = [0] * nima
		flow_cnt  = 0
	
	# [id] classification
	again =  1
	it    = -1
	while again and it < (max_iter - 1):
		it += 1
	
		mpi_barrier(MPI_COMM_WORLD)
		
		if DEBUG:
			flow_ite = 0

			if myid == main_node:
				tg = time.time()
				print '== ITE %d ==' % it

		# [id] compute average send to main_node
		AVE = []
		for k in xrange(K): AVE.append(blank)
		
		if DEBUG:
			MPIlogfile_print(myid, 'after init AVE\n')
			if myid == main_node: t1 = time.time()
		
		for k in xrange(K):
			if myid == grp_asg[k]:
				temp    = kmnr(IM[lut_grp[k]], -1, len(GRP[lut_grp[k]]), -1, numr, wr, check_mirror, max_internal, rand_seed, myid)
				temp   /= len(GRP[lut_grp[k]])
				AVE[k]  = temp.copy()

		if DEBUG and myid == main_node: print 'time average: %d s' % (time.time() - t1)

		if DEBUG: MPIlogfile_print(myid, 'compute ave\n')

		mpi_barrier(MPI_COMM_WORLD)

		for k in xrange(K):
			if myid != main_node:
				if myid == grp_asg[k]:
					send_EMData(AVE[k], main_node, 0)
			else:
				if grp_asg[k] != main_node and grp_asg[k] != -1:
					AVE[k] = recv_EMData(grp_asg[k], 0)

			mpi_barrier(MPI_COMM_WORLD)

		if DEBUG: MPIlogfile_print(myid, 'recv average\n')

		# [main] broadcast each average to other node
		for k in xrange(K):
			tmp = AVE[k].copy()
			bcast_EMData_to_all(tmp, myid, main_node)
			AVE[k] = tmp.copy() # need to use .copy() to store the new im
			mpi_barrier(MPI_COMM_WORLD)

		if DEBUG: MPIlogfile_print(myid, 'after send AVE\n')

		# [all] compute norm for each average
		norm = []
		for k in xrange(K):
			q   = Util.ener(AVE[k], numr)
			res = Util.Crosrng_ew(AVE[k], AVE[k], numr, wr, 0)
			norm.append((2 * q) / float(res['qn']))

		# [id] info display
		MPIlogfile_print(myid, '\n___________ Iteration %d _____________%s\n' % (it + 1, time.strftime('%a_%d_%b_%Y_%H_%M_%S', time.localtime())))
		for k in xrange(K):
			if lut_grp[k] != -1:
				MPIlogfile_print(myid, 'objects in class %3d : %d \n' % (k, len(GRP[lut_grp[k]])))

		MPIlogfile_print(myid, '\n')

		if DEBUG and myid == main_node: t1 = time.time()
		
		# [id] classification
		list_exg = []
		list_dst = []
		again    = 0
		Je_rsd   = 0
		for ngrp in xrange(K):
			if myid == grp_asg[ngrp]:
				exg = []
				dst = []
				for index in xrange(len(GRP[lut_grp[ngrp]])):
					dmin    = 1.0e20
					old_grp = ngrp
					new_grp = -1
					for k in xrange(K):
						if grp_asg[k] != -1: # if the group k is not empty
							if (check_mirror):
								retvals = Util.Crosrng_ew(AVE[k], IM[lut_grp[ngrp]][index], numr, wr, 0)
								qn      = retvals["qn"]
								retvals = Util.Crosrng_ew(AVE[k], IM[lut_grp[ngrp]][index], numr, wr, 1)
								qm      = retvals["qn"]
								q1      = Util.ener(AVE[k], numr)
								q2      = Util.ener(IM[lut_grp[ngrp]][index], numr)
								qn      = max(qn,qm)
								qn      = q1 + q2 - (qn * norm[k])
							else:
								retvals = Util.Crosrng_ew(AVE[k], IM[lut_grp[ngrp]][index], numr, wr, 0)
								qn      = retvals["qn"]
								q1      = Util.ener(AVE[k], numr)
								q2      = Util.ener(IM[lut_grp[ngrp]][index], numr)
								qn      = q1 + q2 - (qn * norm[k])
							if(qn < dmin):
								dmin    = qn
								new_grp = k

					Je_rsd += (dmin / float(norm_rsd))
					
					if new_grp < 0:
						print  'Error in assignment of objects to averages'
						break
					
					## TO TEST
					#if ngrp    == 0: new_grp = 1
					#if new_grp == 0: new_grp = 1
				
					if old_grp != new_grp:
						again = 1
						exg.append(GRP[lut_grp[ngrp]][index])
						dst.append(new_grp)
					
				# store list exg and dst for each loc groups
				list_exg.append(exg)
				list_dst.append(dst)

		if DEBUG and myid == main_node: print 'time classification: %d s' % (time.time() - t1)

		mpi_barrier(MPI_COMM_WORLD)

		if DEBUG: MPIlogfile_print(myid, 'after classification\n')

		# [sync] gather value of criterion Je_rsd
		Je_rsd = mpi_reduce(Je_rsd, 1, MPI_FLOAT, MPI_SUM, main_node, MPI_COMM_WORLD)
		Je_rsd = mpi_bcast(Je_rsd, 1, MPI_FLOAT, main_node, MPI_COMM_WORLD)
		Je_rsd = Je_rsd.tolist()
		Je_rsd = Je_rsd[0]

		MPIlogfile_print(myid, 'Criterion %11.4e\n' % Je_rsd)

		# [sync] with the other node
		again = mpi_reduce(again, 1, MPI_INT, MPI_LOR, main_node, MPI_COMM_WORLD)
		again = mpi_bcast(again, 1, MPI_INT, main_node, MPI_COMM_WORLD)
		again = again.tolist()
		again = again[0]
		
		if DEBUG:
			for grp in xrange(len(list_exg)):
				for index in list_exg[grp]:
					flow_ctrl[index] += 1
			if myid == main_node: t1 = time.time()

		# COMMUNICATION
		mpi_barrier(MPI_COMM_WORLD)
		for node in xrange(number_of_proc * number_of_proc):
			# define the origin and destination node
			org_node = sync_org[node]
			dst_node = sync_dst[node]

			# define list of group inside org node and dst node
			org_grp = []
			dst_grp = []
			for k in xrange(K):
				if grp_asg[k] == org_node:
					org_grp.append(k)
				if grp_asg[k] == dst_node:
					dst_grp.append(k)

			bag           = []
			pack          = []
			pack_size     = []
			bag_im        = []
			pack_im       = []
			lst_remove    = [[] for n in xrange(len(org_grp))]
			if myid == org_node:
				# prepare data and pack
				for outgrp in dst_grp:
					bag        = []
					bag_im     = []
					for ingrp in org_grp:
						remove   = []
						size_exg = len(list_exg[lut_grp[ingrp]])
						for n in xrange(size_exg):
							if list_dst[lut_grp[ingrp]][n] == outgrp:
								index = list_exg[lut_grp[ingrp]][n]
								bag.append(index)
								pos = GRP[lut_grp[ingrp]].index(index)
								bag_im.append(IM[lut_grp[ingrp]][pos])
								remove.append(index)

						lst_remove[lut_grp[ingrp]].extend(remove)
			
					pack.extend(bag)
					pack_im.extend(bag_im)
					pack_size.append(len(bag))
						
			# first send the pack_size to all
			pack_size = mpi_bcast(pack_size, len(dst_grp), MPI_INT, org_node, MPI_COMM_WORLD)
			pack_size = pack_size.tolist()

			mpi_barrier(MPI_COMM_WORLD)

			if DEBUG:
				data_send = 0
				data_recv = 0
				data_in   = 0
			
			# send-recv
			if myid != org_node:
				pack    = []
				pack_im = []
			for n in xrange(sum(pack_size)):
				if org_node != dst_node:
					if myid == org_node:
						# send index
						mpi_send(pack[n], 1, MPI_INT, dst_node, 0, MPI_COMM_WORLD)
						# send header [alpha, sx, sy, mirror]
						head = get_arb_params(pack_im[n], ['alpha', 'alpha_original', 'sx', 'sy', 'mirror'])
						mpi_send(head, 5, MPI_FLOAT, dst_node, 0, MPI_COMM_WORLD)
						# send image
						send_EMData(pack_im[n], dst_node, 1)
						
						if DEBUG: data_send += 1
					
					if myid == dst_node:
						# recv index
						index = mpi_recv(1, MPI_INT, org_node, 0, MPI_COMM_WORLD)
						index = index.tolist()
						pack.append(index[0])
						# recv header
						head    = mpi_recv(5, MPI_FLOAT, org_node, 0, MPI_COMM_WORLD)
						head    = head.tolist()
						head[4] = int(head[4]) # mirror is integer
						# recv image
						img   = recv_EMData(org_node, 1)
						set_arb_params(img, head, ['alpha', 'alpha_original', 'sx', 'sy', 'mirror'])
						pack_im.append(img)
						
						if DEBUG: data_recv += 1
				else:
					if DEBUG: data_in += 1

				mpi_barrier(MPI_COMM_WORLD)
	
			if DEBUG:
				flow_cnt += data_send
				flow_ite += data_send
				if dst_node != org_node:
					if myid == org_node:
						MPIlogfile_print(myid, 'send data %s %d ---> %d %2s: exp [%d] int [%d]\n' % (org_grp, org_node, dst_node, str(dst_grp).ljust(50, ' '), data_send, data_in))
					if myid == dst_node:
						MPIlogfile_print(myid, 'recv data %s %d <--- %d %2s: exp [%d] int [%d]\n' % (dst_grp, dst_node, org_node, str(org_grp).ljust(50, ' '), data_recv, data_in))
				else:
					if myid == org_node:
						MPIlogfile_print(myid, 'loca data %s %d <--> %d %2s: exp [%d] int [%d]\n' % (dst_grp, dst_node, org_node, str(org_grp).ljust(50, ' '), data_recv, data_in))

			if myid == org_node:
				# remove index and images in each group
				for ingrp in org_grp:
					for index in lst_remove[lut_grp[ingrp]]:
						pos = GRP[lut_grp[ingrp]].index(index)
						trash = GRP[lut_grp[ingrp]].pop(pos)
						trash = IM[lut_grp[ingrp]].pop(pos)
						del trash
				
			if myid == dst_node:
				i = 0
				j = 0
				for ingrp in dst_grp:
					bag    = []
					bag_im = []
					for n in xrange(j, j + pack_size[i]):
						bag.append(pack[n])
						bag_im.append(pack_im[n])
						j += 1
						
					i   += 1

					# add index and images in each group
					GRP[lut_grp[ingrp]].extend(bag)
					IM[lut_grp[ingrp]].extend(bag_im)

			# [sync]
		        mpi_barrier(MPI_COMM_WORLD)

		if DEBUG and myid == main_node:
			print 'time comm: %d s' % (time.time() - t1)
			print '\n'
			print 'time iteration: %d s\n\n' % (time.time() - tg)
		
		# [sync]
		mpi_barrier(MPI_COMM_WORLD)

		if DEBUG:
			flow_ite = mpi_reduce(flow_ite, 1, MPI_INT, MPI_SUM, main_node, MPI_COMM_WORLD)
			flow_ite = flow_ite.tolist()

		flag = True
		if flag:
			# search the groups to remove
			grp_rmv = [0] * K
			nbr_obj = [0] * K
			for k in xrange(K):
				if myid == grp_asg[k]:
					if len(GRP[lut_grp[k]]) < K_th:
						grp_rmv[k] = 1
						nbr_obj[k] = len(GRP[lut_grp[k]])

			for k in xrange(K):
				# if the group k need to removed
				if grp_rmv[k] > 0:
					if nbr_obj[k] != 0:
						MPIlogfile_print(myid, 'Warning: class %d not enough images, remove this group\n' % k)
						org_mv = k
						dst_mv = -1

						## TO TEST
						#MPIlogfile_print(myid, 'list lut_grp: %s\n' % lut_grp)
						#MPIlogfile_print(myid, 'list grp rmv: %s\n' % grp_rmv)
						
						# search if in the node there is another group to move objects
						for q in xrange(K):
							if lut_grp[q] != -1 and grp_rmv[q] == 0:
								dst_mv = q

						## TO TEST
						#MPIlogfile_print(myid, 'grp destination %d\n' % dst_mv)
								
						# move object if there is a group destination
						if dst_mv != -1:
							#MPIlogfile_print(myid, 'from %d -> %d\n' % (org_mv, dst_mv))
							GRP[lut_grp[dst_mv]].extend(GRP[lut_grp[org_mv]])
							IM[lut_grp[dst_mv]].extend(IM[lut_grp[dst_mv]])
							MPIlogfile_print(myid, 'move objects grp %d -> grp %d\n' % (org_mv, dst_mv))
					else:
						MPIlogfile_print(myid, 'Warning: class %d is empty, remove this group.\n' % k)
	
			# clean groups
			for k in xrange(K - 1, -1, -1):
				if grp_rmv[k] > 0:
					#MPIlogfile_print(myid, 'lut_grp[k] and k: %d %d\n' % (lut_grp[k], k))
					del GRP[lut_grp[k]]
					del IM[lut_grp[k]]
							
			# broadcast the grp remove
			grp_rmv = mpi_reduce(grp_rmv, K, MPI_INT, MPI_SUM, main_node, MPI_COMM_WORLD)
			grp_rmv = mpi_bcast(grp_rmv, K, MPI_INT, main_node, MPI_COMM_WORLD)
			grp_rmv = grp_rmv.tolist()

			# new group assign node
			for k in xrange(K):
				if grp_rmv[k] > 0: grp_asg[k] = -1			

			# new list loc group assign
			lut_grp = [-1] * K
			i       = 0
			for k in xrange(K):
				if myid == grp_asg[k]:
					lut_grp[k] = i
					i += 1

			# check empty node
			kill = 0
			if len(GRP) == 0:
				print 'Error: cpu %d is empty (no group), kill all process.\n' % myid
				MPIlogfile_print(myid, 'Error: cpu %d is empty (no group), kill all process.\n' % myid)
				kill = 1
			
			mpi_barrier(MPI_COMM_WORLD)
			kill = mpi_reduce(kill, 1, MPI_INT, MPI_SUM, main_node, MPI_COMM_WORLD)
			kill = mpi_bcast(kill, 1, MPI_INT, main_node, MPI_COMM_WORLD)
			kill = kill.tolist()
			kill = kill[0]

			if kill > 0: return

		else:
			# check empty class
			kill = 0
			for k in xrange(K):
				if myid == grp_asg[k]:
					if len(GRP[lut_grp[k]]) < 1:
						print 'Error: class %d is empty' % k
						MPIlogfile_print(myid, 'Error: class %d is empty, kill all process.\n' % k)
						kill = 1

				mpi_barrier(MPI_COMM_WORLD)
			kill = mpi_reduce(kill, 1, MPI_INT, MPI_SUM, main_node, MPI_COMM_WORLD)
			kill = mpi_bcast(kill, 1, MPI_INT, main_node, MPI_COMM_WORLD)
			kill = kill.tolist()
			kill = kill[0]

			if kill > 0: return

	if DEBUG:
		flow_ctrl = mpi_reduce(flow_ctrl, nima, MPI_INT, MPI_SUM, main_node, MPI_COMM_WORLD)
		flow_ctrl = flow_ctrl.tolist()

		print myid, ': flow count [%d]' % flow_cnt

		if myid == main_node:
			print '\n ############ flow control ####################\n'
			#print TXT_flow
			#print flow_ctrl
			print '\n'
			print 'sum: %d' % sum(flow_ctrl)
	
	mpi_barrier(MPI_COMM_WORLD)

	# display the number of objects in each group
	ngrp = [0] * K
	for k in xrange(K):
		if myid == grp_asg[k]: ngrp[k] = len(GRP[lut_grp[k]])
	ngrp = mpi_reduce(ngrp, K, MPI_INT, MPI_SUM, main_node, MPI_COMM_WORLD)
	ngrp = mpi_bcast(ngrp, K, MPI_INT, main_node, MPI_COMM_WORLD)
	ngrp = ngrp.tolist()
	MPIlogfile_print(myid, '\n== Object in each Class ==\n')
	count = 0
	for k in xrange(K):
		if ngrp[k] == 0:
			MPIlogfile_print(myid, 'class %3d : %d [removed]\n' % (k, ngrp[k]))
		else:
			MPIlogfile_print(myid, 'class %3d : %d\n' % (k, ngrp[k]))
			count += 1

	MPIlogfile_print(myid, 'number of classes: %d\n' % count)
	
	MPIlogfile_print(myid, '\nPREPARE RESULTS\n')

	'''
	# update ave according the empty group
	newAVE = []
	L      = 0
	for k in xrange(K):
		if ngrp[k] != 0:
			tmp = AVE[k].copy()
			newAVE.append(tmp)
			L += 1
	K = L
	'''

	# [main] align class averages and transfer parameters to individual images
	talpha = [0] * K
	tmir   = [0] * K
	if myid == main_node:
		for k in xrange(K): AVE[k].set_attr_dict({'alpha':1.0, 'mirror':0})
		kmn(AVE, numr, wr, check_mirror, max_iter, rand_seed)
		for k in xrange(K):
			talpha[k] = ang_n(AVE[k].get_attr('alpha'), mode, lnumr)
			tmir[k]   = AVE[k].get_attr('mirror')

		del AVE
	
	# [all] broadcast talpha and tmir
	talpha = mpi_bcast(talpha, K, MPI_FLOAT, main_node, MPI_COMM_WORLD)
	talpha = talpha.tolist()
	tmir   = mpi_bcast(tmir, K, MPI_INT, main_node, MPI_COMM_WORLD)
	tmir   = tmir.tolist()
	
	mpi_barrier(MPI_COMM_WORLD)

	#  write out the alignment parameters to headers
	del temp
	temp = EMData()

	if DEBUG and myid == main_node: t1 = time.time()

	for k in xrange(K):
		if myid == grp_asg[k]:
			for n in xrange(len(GRP[lut_grp[k]])):
				#  First combine with angle of the average
				alpha  = ang_n(IM[lut_grp[k]][n].get_attr('alpha'), mode, lnumr)
				mirror =  IM[lut_grp[k]][n].get_attr('mirror')
				alpha, tmp, it, mirror = combine_params2(alpha, 0, 0, mirror, talpha[k], 0, 0, tmir[k])
			
				#  Second combine with given alignment
				alpha_original = IM[lut_grp[k]][n].get_attr('alpha_original')
				sx    =  IM[lut_grp[k]][n].get_attr('sx')
				sy    =  IM[lut_grp[k]][n].get_attr('sy')
				alpha_original_n, sxn, syn, mir = combine_params2(0, -sx, -sy, 0, -alpha_original, 0, 0, 0)
				alphan, sxn, syn, mir           = combine_params2(0, -sxn, -syn, 0, alpha, 0, 0, mirror)
				temp.read_image(stack, GRP[lut_grp[k]][n], True)
				temp.set_attr_dict({'alpha':alphan, 'sx':sxn, 'sy':syn, 'mirror': mir, 'nclass':K, 'ref_num':k})
				temp.write_image(stack, GRP[lut_grp[k]][n], EMUtil.ImageType.IMAGE_HDF, True)

		mpi_barrier(MPI_COMM_WORLD)

	if DEBUG and myid == main_node: print 'time write header: %d s' % (time.time() - t1)

	del temp
	del IM

	MPIlogfile_print(myid, '\n Time: %d s\n' % (time.time() - t_start))
	MPIlogfile_end(myid)



	
def SSNR_func(args, data):
	from utilities import print_msg
	
	img_data = data[0]
	nima = data[1]
	maskI = data[2]
	kb = data[3]
	CTF = data[6]
	SSNR_fit = data[7]
	if CTF:
		index_list = data[8]
		ctfimg_list = data[9]
		ctfimg2 = data[10]
		ctf_data = 3
	else:
		ctf_data = 0
	if SSNR_fit:
		SSNR_img = data[8+ctf_data]
		SSNR_sig = data[9+ctf_data]
		
	N = img_data[0].get_ysize()
	avgimg  = EMData(N, N, 1, False)
	var     = EMData(N, N, 1, False)
	
	for im in xrange(nima):
		imgft = img_data[im].copy()
		alpha = args[im*3]
		sx = args[im*3+1]
		sy = args[im*3+2]
		#imgft = gridrot_shift2D(imgft, kb, alpha, sx, sy)
		imgft = imgft.fouriergridrot_shift2d(alpha, sx, sy, kb)
		imgft.center_origin_fft()
		Util.add_img2(var, imgft)
		if CTF:
			ctfimg = ctfimg_list[index_list[im]]
			Util.mul_img(imgft, ctfimg)
		Util.add_img(avgimg, imgft)

		#print "%2d %20.12e %20.12e"%(im, imgft.get_value_at(2, 0), imgft.get_value_at(3, 0))

	'''
	yyy1 = avgimg.get_value_at(2, 0)
	yyy2 = avgimg.get_value_at(3, 0)
	yyy3 = var.get_value_at(2, 0)
	print yyy1, yyy2, yyy3
	yyy1 /= nima
	yyy2 /= nima
	yyy = yyy1**2+yyy2**2
	print yyy*nima, (yyy3-nima*yyy)/(nima-1), yyy*nima/((yyy3-nima*yyy)/(nima-1))
	'''
	
	if CTF:
		Util.div_filter(avgimg, ctfimg2)
	else:
		Util.mul_scalar(avgimg, 1.0/float(nima))		
	avgimg_2 = avgimg.copy()
	Util.mul_img(avgimg_2, avgimg_2.conjg())
	if CTF:
		Util.mul_img(avgimg_2, ctfimg2)
	else:
		Util.mul_scalar(avgimg_2, float(nima))
	Util.sub_img(var, avgimg_2)
	Util.mul_scalar(var, 1.0/float(nima-1))
	
	SSNR = avgimg_2.copy()
	Util.div_filter(SSNR, var)
	
	#print avgimg_2.get_value_at(2, 0), var.get_value_at(2, 0), SSNR.get_value_at(2, 0)
	#print " "
	
	a0 = Util.infomask(SSNR, maskI, True)
	sum_SSNR = a0[0]
	print_msg("SSNR = %20.7f\n"%(sum_SSNR))
	
	if SSNR_fit:
		beta = 10.0
		Util.sub_img(SSNR, maskI)
		SSNR = SSNR.process("threshold.belowtozero", {"minval": 0.0})
		Util.sub_img(SSNR, SSNR_img)
		Util.div_filter(SSNR, SSNR_sig)
		Util.mul_img(SSNR, SSNR)
		a0 = Util.infomask(SSNR, maskI, True)
		sum_SSNR -= beta*a0[0]
	
	return -sum_SSNR


def SSNR_grad(args, data):
	from numpy import zeros, array, float64

	img_data = data[0]
	nima = data[1]
	maskI = data[2]
	kb = data[3]
	_jX = data[4]
	_jY = data[5]
	CTF = data[6]
	SSNR_fit = data[7]
	if CTF:
		index_list = data[8]
		ctfimg_list = data[9]
		ctfimg2 = data[10]
		ctf_data = 3
	else:
		ctf_data = 0
	if SSNR_fit:
		SSNR_img = data[8+ctf_data]
		SSNR_sig = data[9+ctf_data]
	
	N = img_data[0].get_ysize()
	avgimg  = EMData(N, N, 1, False)
	var     = EMData(N, N, 1, False)
	img_data_new = []
	d_img = []

	for im in xrange(nima):
		imgft = img_data[im].copy()
		alpha = args[im*3]
		sx = args[im*3+1]
		sy = args[im*3+2]
		
		dalpha = 0.05
		#imgft0 = gridrot_shift2D(imgft, kb, alpha, sx, sy)
		#imgft1 = gridrot_shift2D(imgft, kb, alpha-dalpha, sx, sy)
		#imgft2 = gridrot_shift2D(imgft, kb, alpha+dalpha, sx, sy)
		imgft0 = imgft.fouriergridrot_shift2d(alpha, sx, sy, kb)
		imgft0.center_origin_fft()
		imgft1 = imgft.fouriergridrot_shift2d(alpha-dalpha, sx, sy, kb)
		imgft1.center_origin_fft()
		imgft2 = imgft.fouriergridrot_shift2d(alpha+dalpha, sx, sy, kb)
		imgft2.center_origin_fft()
		Util.sub_img(imgft2, imgft1)
		Util.mul_scalar(imgft2, 1/(2*dalpha))
				
		img_data_new.append(imgft0)
		d_img.append(imgft2)
		Util.add_img2(var, imgft0)
		if CTF:
			ctfimg = ctfimg_list[index_list[im]]
			Util.mul_img(imgft0, ctfimg)
		Util.add_img(avgimg, imgft0)
		
	if CTF:
		Util.div_filter(avgimg, ctfimg2)
	else:
		Util.mul_scalar(avgimg, 1.0/float(nima))
	avgimg_2 = avgimg.copy()
	Util.mul_img(avgimg_2, avgimg_2.conjg())
	if CTF:
		Util.mul_img(avgimg_2, ctfimg2)
	else:
		Util.mul_scalar(avgimg_2, float(nima))
	sumimg2 = var.copy()
	Util.sub_img(var, avgimg_2)
	Util.mul_scalar(var, 1.0/float(nima-1))
	
	avgimg_conj = avgimg.conjg()
	dSSNR = avgimg_conj.copy()
	Util.div_filter(dSSNR, var)

	if SSNR_fit:
		beta = 10.0
		SSNR = avgimg_2.copy()
		Util.div_filter(SSNR, var)
		Util.sub_img(SSNR, maskI)
		SSNR = SSNR.process("threshold.belowtozero", {"minval": 0.0})
		Util.sub_img(SSNR, SSNR_img)
		Util.div_filter(SSNR, SSNR_sig)
		Util.div_filter(SSNR, SSNR_sig)
		Util.mul_scalar(SSNR, 2*beta)
		C = maskI.copy()
		Util.sub_img(C, SSNR)
		Util.mul_img(dSSNR, C)
	
	g = zeros(args.shape, float64)
	accurate = True

	for im in xrange(nima):
		img_new = img_data_new[im].copy()
		dSSNR_copy = dSSNR.copy()

		if accurate: 
			
			img_new_copy = img_new.copy()
			Util.sub_img(img_new_copy, avgimg)
			Util.mul_img(img_new_copy, dSSNR)
			img_new_copy = img_new_copy.conjg()
			Util.mul_scalar(img_new_copy, nima/float(nima-1))
			C = maskI.copy()
			Util.sub_img(C, img_new_copy)
			Util.mul_img(dSSNR_copy, C)
			'''
			C = sumimg2.copy()
			img_new_congj = img_new.conjg()
			Util.mul_img(img_new_congj, avgimg)
			Util.mul_scalar(img_new_congj, float(nima))
			Util.sub_img(C, img_new_congj)
			Util.div_filter(C, var)
			Util.mul_scalar(C, 1.0/float(nima-1))
			Util.mul_img(dSSNR_copy, C)
			'''
			
		Util.mul_img(dSSNR_copy, d_img[im])
		if CTF:
			Util.mul_img(dSSNR_copy, ctfimg_list[index_list[im]])
		Util.add_img(dSSNR_copy, dSSNR_copy.conjg())
		a0 = Util.infomask(dSSNR_copy, maskI, True)
		g[im*3] = -a0[0]
		
		dSSNR_copy = dSSNR.copy()
		if accurate:
			Util.mul_img(dSSNR_copy, C)
		Util.mul_img(dSSNR_copy, img_new)
		if CTF:
			Util.mul_img(dSSNR_copy, ctfimg_list[index_list[im]])
		dSSNR_fft = dSSNR_copy.copy()		
		Util.mul_img(dSSNR_fft, _jX)
		Util.add_img(dSSNR_fft, dSSNR_fft.conjg())
		
		#if im == 0: 	print dSSNR_fft.get_value_at(2, 0)
		
		a0 = Util.infomask(dSSNR_fft, maskI, True)
		g[im*3+1] = -a0[0]
		
		dSSNR_fft = dSSNR_copy.copy()		
		Util.mul_img(dSSNR_fft, _jY)
		Util.add_img(dSSNR_fft, dSSNR_fft.conjg())
		a0 = Util.infomask(dSSNR_fft, maskI, True)
		g[im*3+2] = -a0[0]
	return g


def SSNR_func_MPI(args, data):
	from applications import MPI_start_end
	from utilities import reduce_EMData_to_root, print_msg
	from mpi import mpi_comm_size, mpi_comm_rank, mpi_bcast, MPI_COMM_WORLD, MPI_FLOAT

	number_of_proc = mpi_comm_size(MPI_COMM_WORLD)
	myid = mpi_comm_rank(MPI_COMM_WORLD)
	main_node = 0

	img_data = data[0]
	nima = data[1]
	maskI = data[2]
	kb = data[3]
	CTF = data[6]
	SSNR_fit = data[7]
	if CTF:
		index_list = data[8]
		ctfimg_list = data[9]
		ctfimg2 = data[10]
		ctf_data = 3
	else:
		ctf_data = 0
	if SSNR_fit:
		"""
		SSNR_img = data[8+ctf_data]
		SSNR_sig = data[9+ctf_data]
		"""
		SSNR_r = data[8+ctf_data]
	
	image_start, image_end = MPI_start_end(nima, number_of_proc, myid)
	
	N = img_data[0].get_ysize()
	avgimg  = EMData(N, N, 1, False)
	var     = EMData(N, N, 1, False)

	for im in xrange(image_start, image_end):
		imgft = img_data[im-image_start].copy()
		alpha = args[im*3]
		sx = args[im*3+1]
		sy = args[im*3+2]
		#imgft = gridrot_shift2D(imgft, kb, alpha, sx, sy)			
		imgft = imgft.fouriergridrot_shift2d(alpha, sx, sy, kb)
		imgft.center_origin_fft()
		Util.add_img2(var, imgft)
		if CTF:
			ctfimg = ctfimg_list[index_list[im]]
			Util.mul_img(imgft, ctfimg)
		Util.add_img(avgimg, imgft)
		
	reduce_EMData_to_root(avgimg, myid, main_node)
	reduce_EMData_to_root(var, myid, main_node)
	if myid == main_node:
		if CTF:
			Util.div_filter(avgimg, ctfimg2)
		else:
			Util.mul_scalar(avgimg, 1.0/float(nima))
		avgimg_2 = avgimg.copy()
		Util.mul_img(avgimg_2, avgimg_2.conjg())
		if CTF:
			Util.mul_img(avgimg_2, ctfimg2)
		else:
			Util.mul_scalar(avgimg_2, float(nima))
		Util.sub_img(var, avgimg_2)
		Util.mul_scalar(var, 1.0/float(nima-1))
	
		SSNR = avgimg_2.copy()
		Util.div_filter(SSNR, var)
		
		a0 = Util.infomask(SSNR, maskI, True)
		sum_SSNR = a0[0]
		print_msg("SSNR = %20.7f\n"%(sum_SSNR))

		if SSNR_fit:
			from fundamentals import rot_avg_table
			beta = 10.0
			avgimg_2 = Util.pack_complex_to_real(avgimg_2)
			var   = Util.pack_complex_to_real(var)
			ravgimg_2 = rot_avg_table(avgimg_2)
			rvar = rot_avg_table(var)
			SSNR_diff = 0.0
			for i in xrange(N/2+1):
				qt = max(0.0, ravgimg_2[i]/rvar[i] - 1.0)
				diff = (qt-SSNR_r[i])/max(1.0, SSNR_r[i])
				print "In bin %3d  Target SSNR = %10.4f  Actual SSNR= %10.4f difference in pct is %7.3f"%(i,SSNR_r[i], qt, diff)
				print "Numerator: ", ravgimg_2[i], "Denominator: ", rvar[i]
				SSNR_diff += diff**2
			sum_SSNR -= beta*SSNR_diff				
			"""
			Util.sub_img(SSNR, maskI)
			SSNR = SSNR.process("threshold.belowtozero", {"minval": 0.0})
			Util.sub_img(SSNR, SSNR_img)
			Util.div_filter(SSNR, SSNR_sig)
			Util.mul_img(SSNR, SSNR)
			a0 = Util.infomask(SSNR, maskI, True)
			sum_SSNR -= beta*a0[0]
			"""
	else: 	
		sum_SSNR = 0.0
	
	sum_SSNR = mpi_bcast(sum_SSNR, 1, MPI_FLOAT, main_node, MPI_COMM_WORLD)
	return -sum_SSNR[0]


def SSNR_grad_MPI(args, data):
	from numpy import zeros, array, float32, float64
	from applications import MPI_start_end
	from utilities import reduce_EMData_to_root, bcast_EMData_to_all 
	from utilities import reduce_array_to_root, bcast_array_to_all
	from mpi import mpi_comm_size, mpi_comm_rank, mpi_bcast, MPI_COMM_WORLD, MPI_FLOAT

	number_of_proc = mpi_comm_size(MPI_COMM_WORLD)
	myid = mpi_comm_rank(MPI_COMM_WORLD)
	main_node = 0

	img_data = data[0]
	nima = data[1]
	maskI = data[2]
	kb = data[3]
	_jX = data[4]
	_jY = data[5]
	CTF = data[6]
	SSNR_fit = data[7]
	if CTF:
		index_list = data[8]
		ctfimg_list = data[9]
		ctfimg2 = data[10]
		ctf_data = 3
	else:
		ctf_data = 0
	if SSNR_fit:
		"""
		SSNR_img = data[8+ctf_data]
		SSNR_sig = data[9+ctf_data]
		"""
		SSNR_r = data[8+ctf_data]
	
	image_start, image_end = MPI_start_end(nima, number_of_proc, myid)

	N = img_data[0].get_ysize()
	avgimg  = EMData(N, N, 1, False)
	var     = EMData(N, N, 1, False)
	img_data_new = []
	d_img = []

	for im in xrange(image_start, image_end):
		imgft = img_data[im-image_start].copy()
		alpha = args[im*3]
		sx = args[im*3+1]
		sy = args[im*3+2]
		
		dalpha = 0.05
		#imgft0 = gridrot_shift2D(imgft, kb, alpha, sx, sy)
		#imgft1 = gridrot_shift2D(imgft, kb, alpha-dalpha, sx, sy)
		#imgft2 = gridrot_shift2D(imgft, kb, alpha+dalpha, sx, sy)
		imgft0 = imgft.fouriergridrot_shift2d(alpha, sx, sy, kb)
		imgft0.center_origin_fft()		
		imgft1 = imgft.fouriergridrot_shift2d(alpha-dalpha, sx, sy, kb)
		imgft1.center_origin_fft()
		imgft2 = imgft.fouriergridrot_shift2d(alpha+dalpha, sx, sy, kb)
		imgft2.center_origin_fft()
		Util.sub_img(imgft2, imgft1)
		Util.mul_scalar(imgft2, 1/(2*dalpha))
				
		img_data_new.append(imgft0)
		d_img.append(imgft2)
		Util.add_img2(var, imgft0)
		if CTF:
			ctfimg = ctfimg_list[index_list[im]]
			Util.mul_img(imgft0, ctfimg)
		Util.add_img(avgimg, imgft0)
		
	reduce_EMData_to_root(avgimg, myid, main_node)
	reduce_EMData_to_root(var, myid, main_node)
	bcast_EMData_to_all(avgimg, myid, main_node)
	bcast_EMData_to_all(var, myid, main_node)

	if CTF:
		Util.div_filter(avgimg, ctfimg2)
	else:
		Util.mul_scalar(avgimg, 1.0/float(nima))
	avgimg_2 = avgimg.copy()
	Util.mul_img(avgimg_2, avgimg_2.conjg())
	if CTF:
		Util.mul_img(avgimg_2, ctfimg2)
	else:
		Util.mul_scalar(avgimg_2, float(nima))
	Util.sub_img(var, avgimg_2)
	Util.mul_scalar(var, 1.0/float(nima-1))
	
	avgimg_conj = avgimg.conjg()
	dSSNR = avgimg_conj.copy()
	Util.div_filter(dSSNR, var)

	if SSNR_fit:
		from fundamentals import rot_avg_table
		from math import sqrt, pi
		beta = 10.0
		avgimg_2 = Util.pack_complex_to_real(avgimg_2)
		var   = Util.pack_complex_to_real(var)
		ravgimg_2 = rot_avg_table(avgimg_2)
		rvar = rot_avg_table(var)
		C = EMData(N, N, 1, False)
		S = pi*(0.49*N)**2
		for x in xrange((N+2)/2):
			for y in xrange(N):
 				if y > N/2-1: yy = y-N
				else: yy = y
				r = sqrt(x**2+yy**2)
				if r < 0.49*N:
					i = int(r+0.5)
					qt = max(0.0, ravgimg_2[i]/rvar[i] - 1.0)
					temp = 1-2*beta*(qt-SSNR_r[i])/max(1.0, SSNR_r[i])**2*S/max(1, 2*pi*i)
					C.set_value_at(x*2, y, temp)
		Util.mul_img(dSSNR, C)
		
		"""
		SSNR = avgimg_2.copy()
		Util.div_filter(SSNR, var)
		Util.sub_img(SSNR, maskI)
		SSNR = SSNR.process("threshold.belowtozero", {"minval": 0.0})
		Util.sub_img(SSNR, SSNR_img)
		Util.div_filter(SSNR, SSNR_sig)
		Util.div_filter(SSNR, SSNR_sig)
		Util.mul_scalar(SSNR, 2*beta)
		C = maskI.copy()
		Util.sub_img(C, SSNR)
		Util.mul_img(dSSNR, C)
		"""
	
	h = zeros(args.shape, float32)
	accurate = True
	
	for im in xrange(image_start, image_end):
		img_new = img_data_new[im-image_start].copy()
		dSSNR_copy = dSSNR.copy()

		if accurate: 
			img_new_copy = img_new.copy()
			Util.sub_img(img_new_copy, avgimg)
			Util.mul_img(img_new_copy, dSSNR)
			img_new_copy = img_new_copy.conjg()
			Util.mul_scalar(img_new_copy, nima/float(nima-1))
			C = maskI.copy()
			Util.sub_img(C, img_new_copy)
			Util.mul_img(dSSNR_copy, C)

		Util.mul_img(dSSNR_copy, d_img[im-image_start])
		if CTF:
			Util.mul_img(dSSNR_copy, ctfimg_list[index_list[im]])
		Util.add_img(dSSNR_copy, dSSNR_copy.conjg())
		a0 = Util.infomask(dSSNR_copy, maskI, True)
		h[im*3] = -a0[0]
		
		dSSNR_copy = dSSNR.copy()
		if accurate:
			Util.mul_img(dSSNR_copy, C)
		Util.mul_img(dSSNR_copy, img_new)
		if CTF:
			Util.mul_img(dSSNR_copy, ctfimg_list[index_list[im]])
		
		dSSNR_fft = dSSNR_copy.copy()		
		Util.mul_img(dSSNR_fft, _jX)
		Util.add_img(dSSNR_fft, dSSNR_fft.conjg())
		a0 = Util.infomask(dSSNR_fft, maskI, True)
		h[im*3+1] = -a0[0]
		
		dSSNR_fft = dSSNR_copy.copy()		
		Util.mul_img(dSSNR_fft, _jY)
		Util.add_img(dSSNR_fft, dSSNR_fft.conjg())
		a0 = Util.infomask(dSSNR_fft, maskI, True)
		h[im*3+2] = -a0[0]
		
	reduce_array_to_root(h, myid, main_node)
	bcast_array_to_all(h, myid, main_node)
	
	g = zeros(args.shape, float64)
	g = float64(h)
	return g


def ali_SSNR(stack, maskfile=None, ou=-1, maxit=10, CTF=False, opti_method="CG", SSNR_fit=False, SSNR=[], MPI=False):

	if MPI:
		ali_SSNR_MPI(stack, maskfile, ou, maxit, CTF, opti_method, SSNR_fit, SSNR)
		return
		
	from math import pi, sqrt
	from fundamentals import fftip, mirror
	from numpy import Inf
	from scipy.optimize.lbfgsb import fmin_l_bfgs_b
	from scipy.optimize.optimize import fmin_cg
	from utilities import get_image, get_params2D, set_params2D, print_begin_msg, print_end_msg, print_msg
	
	if CTF:
		from utilities import get_arb_params
		from morphology import ctf_img

	print_begin_msg("ali_SSNR")

	if opti_method!="CG" and opti_method!="LBFGSB":
		 print "Unknown optimization method!"
		 return
		 
	nima = EMUtil.get_image_count(stack)
	ima = EMData()
	ima.read_image(stack, 0)
	nx = ima.get_xsize()
	
	last_ring = int(ou);	max_iter = int(maxit)
	if last_ring == -1:	last_ring = nx//2-2
	
	print_msg("Input stack                 : %s\n"%(stack))
	print_msg("Outer radius                : %i\n"%(last_ring))
	print_msg("Maximum iteration           : %i\n"%(max_iter))
	print_msg("Data with CTF               : %s\n"%(CTF))
	print_msg("Optimization method         : %s\n"%(opti_method))

	if maskfile:
		import  types
		if type(maskfile) is types.StringType:  
			print_msg("Maskfile                    : %s\n\n"%(maskfile))
			mask = get_image(maskfile)
		else:
			print_msg("Maskfile                    : user provided in-core mask\n\n")
			mask = maskfile
	else :
		from utilities import model_circle 
		print_msg("Maskfile                    : default, a circle with radius %i\n\n"%(last_ring))
		mask = model_circle(last_ring, nx, nx)
	
	
	# prepare kb for gridding interpolation
	npad = 2
	N = nx*npad
	K = 6
	alpha = 1.75
	r = nx/2
	v = K/2.0/N
	kb = Util.KaiserBessel(alpha, K, r, v, N)
	del alpha, K, r, v

	# generate the mask in Fourier space
	maskI = EMData(N, N, 1, False)
	for x in xrange((N+2)/2):
		for y in xrange(N):
 			if y > N/2-1: yy = y-N
			else: yy = y
			if x**2+yy**2 < (N*0.49)**2:
				maskI.set_value_at(x*2, y, 1) 
	maskI.set_value_at(0, 0, 0)
	maskI.set_value_at(1, 0, 0)
	
	if CTF:
		defocus_list = []
		ctfimg_list = []
		index_list = []
		ctfimg2 = EMData(N, N, 1, False)

	# There two matrixes in Fourier space is necessary for calculating derivatives
	_jX = EMData(N, N, 1, False)
	_jY = EMData(N, N, 1, False)
	for x in xrange((N+2)/2):
		for y in xrange(N):
 			if y > N/2-1: yy = y-N
			else: yy = y
		 	_jX.set_value_at(x*2+1, y, -x*2*pi/N) 
 			_jY.set_value_at(x*2+1, y, -yy*2*pi/N) 

	img_data= []
	x0 = [0.0]*(nima*3)
	mir = [0]*nima	
	if opti_method == "LBFGSB":	bounds = []

	# pre-processing, get initial parameters and boundaries (for LBFGSB)
	for im in xrange(nima):
		if im>0:
			ima = EMData()
			ima.read_image(stack, im)
		if CTF:
			ctf_params = ima.get_attr('ctf')

		x0[im*3], x0[im*3+1], x0[im*3+2], mir[im], dummy = get_params2D(ima)
		
		if mir[im] == 1:
			x0[im*3] = 360.0-x0[im*3]
			x0[im*3+1] = -x0[im*3+1]
			ima = mirror(ima)

		st = Util.infomask(ima, mask, False)
		ima -= st[0]	
		ima.divkbsinh(kb)
		ima = ima.norm_pad(False, npad)
		fftip(ima)
		ima.center_origin_fft()
		img_data.append(ima)
		if CTF:
			if not (ctf_params.defocus in defocus_list):
				defocus_list.append(ctf_params.defocus)
				ctfimg = ctf_img(N, ctf_params, ny = N, nz = 1)
				ctfimg_list.append(ctfimg)
				index_list.append(len(defocus_list)-1)
			else:
				index = defocus_list.index(ctf_params.defocus)
				ctfimg = ctfimg_list[index]
				index_list.append(index)
			Util.add_img2(ctfimg2, ctfimg)			

		if opti_method == "LBFGSB":
			bounds.append((x0[im*3]-2.0, x0[im*3]+2.0))
			bounds.append((x0[im*3+1]-1.0, x0[im*3+1]+1.0))
			bounds.append((x0[im*3+2]-1.0, x0[im*3+2]+1.0))

	# Use a gradient method here
	data = []
	data.append(img_data)
	data.append(nima)
	data.append(maskI)
	data.append(kb)
	data.append(_jX)
	data.append(_jY)
	data.append(CTF)
	data.append(SSNR_fit)

	if CTF:
		data.append(index_list)
		data.append(ctfimg_list)
		data.append(ctfimg2)
	if SSNR_fit:
		"""
		SSNR_img = EMData(N, N, 1, False)
		SSNR_sig = EMData(N, N, 1, False)
		for x in xrange((N+2)/2):
			for y in xrange(N):
 				if y > N/2-1: yy = y-N
				else: yy = y
				r = sqrt(x**2+yy**2)
				if r < N*0.49:
					i = int(r/2)
					j = r/2-i
					temp = (SSNR[i]*(1-j)+SSNR[i+1]*j)*nima
					SSNR_img.set_value_at(x*2, y, temp) 
					if temp < 1.0: temp = 1.0
					SSNR_sig.set_value_at(x*2, y, temp)					
		data.append(SSNR_img)
		data.append(SSNR_sig)
		"""
		SSNR_r = []
		for i in xrange(N/2+1):
			if i%2==0: SSNR_r.append(SSNR[i/2]*nima)
			else:	SSNR_r.append((SSNR[i/2]+SSNR[i/2+1])*0.5*nima)
		data.append(SSNR_r)
	
	from numpy import array
	x1 = array(x0)
	aaa = SSNR_func(x1, data)
	bbb = SSNR_grad(x1, data)

	ccc = [0.0]*(nima*3)
	for i in xrange(nima*3):
		x1[i] += 0.1
		ccc[i] = SSNR_func(x1, data)
		x1[i] -= 0.1
	t1 = []
	t2 = []
	t3 = []
	for i in xrange(nima*3):
		f1 = ccc[i]-aaa
		f2 = bbb[i]*0.1
		if f1!=0.0: pct = abs(f1-f2)/abs(f1)
		else: pct = 0.0
		if i%3==0: t1.append(pct)
		elif i%3==1: t2.append(pct)
		else: t3.append(pct)
		print i, f1, f2, pct
	t1.sort(), t2.sort(), t3.sort()
	print "Median error = ", t1[nima/2], t2[nima/2], t3[nima/2]
	exit()
	
	if opti_method == "CG":
		ps = fmin_cg(SSNR_func, x0, fprime=SSNR_grad, args=([data]), gtol=1e-3, norm=Inf, epsilon=1e-5,
        	      maxiter=max_iter, full_output=0, disp=0, retall=0, callback=None)		
	else:
		ps, val, d = fmin_l_bfgs_b(SSNR_func, x0, args=[data], fprime=SSNR_grad, bounds=bounds, m=10, 
			factr=1e3, pgtol=1e-4, epsilon=1e-2, iprint=-1, maxfun=max_iter)
	
	# write the result into the header
	for im in xrange(nima):
		ima = EMData()
		ima.read_image(stack, im)
		if mir[im] == 0:
			set_params2D(ima, [ps[im*3], ps[im*3+1], ps[im*3+2], 0, 1.0])
		else:
			set_params2D(ima, [360.0-ps[im*3], -ps[im*3+1], ps[im*3+2], 1, 1.0])
		ima.write_image(stack, im)	

	print_end_msg("ali_SSNR")


def ali_SSNR_MPI(stack, maskfile=None, ou=-1, maxit=10, CTF=False, opti_method="CG", SSNR_fit=False, SSNR=[]):

	from applications import MPI_start_end
	from math import pi, sqrt
	from fundamentals import fftip, mirror
	from numpy import Inf
	from scipy.optimize.lbfgsb import fmin_l_bfgs_b
	from scipy.optimize.optimize import fmin_cg
	from mpi import mpi_comm_size, mpi_comm_rank, mpi_bcast, MPI_COMM_WORLD, mpi_barrier
	from utilities import get_image, set_params2D, get_params2D, print_begin_msg, print_end_msg, print_msg
	from utilities import model_circle
	
	if CTF:
		from utilities import get_arb_params
		from morphology import ctf_img

	number_of_proc = mpi_comm_size(MPI_COMM_WORLD)
	myid = mpi_comm_rank(MPI_COMM_WORLD)
	main_node = 0
	
	if myid == main_node:	print_begin_msg("ali_SSNR_MPI")

	if opti_method!="CG" and opti_method!="LBFGSB":
		 print "Unknown optimization method!"
		 return
		 
	nima = EMUtil.get_image_count(stack)
	image_start, image_end = MPI_start_end(nima, number_of_proc, myid)
	ima = EMData()
	ima.read_image(stack, 0)
	nx = ima.get_xsize()
	
	last_ring = int(ou);	max_iter = int(maxit)
	if last_ring == -1:	last_ring = nx//2-2
	
	if myid == main_node:
		print_msg("Input stack                 : %s\n"%(stack))
		print_msg("Outer radius                : %i\n"%(last_ring))
		print_msg("Maximum iteration           : %i\n"%(max_iter))
		print_msg("Data with CTF               : %s\n"%(CTF))
		print_msg("Optimization method         : %s\n"%(opti_method))
		print_msg("Number of processors used   : %d\n"%(number_of_proc))

	if maskfile:
		import  types
		if type(maskfile) is types.StringType:  
			if myid == main_node:		print_msg("Maskfile                    : %s\n\n"%(maskfile))
			mask = get_image(maskfile)
		else:
			if myid == main_node: 		print_msg("Maskfile                    : user provided in-core mask\n\n")
			mask = maskfile
	else : 
		if myid==main_node: 	print_msg("Maskfile                    : default, a circle with radius %i\n\n"%(last_ring))
		mask = model_circle(last_ring, nx, nx)
	
	
	# prepare kb for gridding interpolation
	npad = 2
	N = nx*npad
	K = 6
	alpha = 1.75
	r = nx/2
	v = K/2.0/N
	kb = Util.KaiserBessel(alpha, K, r, v, N)
	del alpha, K, r, v

	# generate the mask in Fourier space
	maskI = EMData(N, N, 1, False)
	for x in xrange((N+2)/2):
		for y in xrange(N):
 			if y > N/2-1: yy = y-N
			else: yy = y
			if x**2+yy**2 < (N*0.49)**2:
				maskI.set_value_at(x*2, y, 1) 
	maskI.set_value_at(0, 0, 0)
	maskI.set_value_at(1, 0, 0)
	
	if CTF:
		defocus_list = []
		ctfimg_list = []
		index_list = []
		ctfimg2 = EMData(N, N, 1, False)

	# There two matrixes in Fourier space is necessary for calculating derivatives
	_jX = EMData(N, N, 1, False)
	_jY = EMData(N, N, 1, False)
	for x in xrange((N+2)/2):
		for y in xrange(N):
 			if y > N/2-1: yy = y-N
			else: yy = y
		 	_jX.set_value_at(x*2+1, y, -x*2*pi/N) 
 			_jY.set_value_at(x*2+1, y, -yy*2*pi/N) 

	img_data= []
	x0 = [0.0]*(nima*3)
	mir = [0]*nima
	if opti_method == "LBFGSB":	bounds = []

	# pre-processing, get initial parameters and boundaries (for LBFGSB)
	for im in xrange(nima):
		if im>0:
			ima = EMData()
			ima.read_image(stack, im)
		if CTF:
			ctf_params = ima.get_attr('ctf')

		x0[im*3], x0[im*3+1], x0[im*3+2], mir[im], dummy = get_params2D(ima)
		
		if mir[im] == 1:
			x0[im*3] = 360-x0[im*3]
			x0[im*3+1] = -x0[im*3+1]
			ima = mirror(ima)		

		if (im >= image_start) and (im < image_end):
			st = Util.infomask(ima, mask, False)
			ima -= st[0]	
			ima.divkbsinh(kb)
			ima = ima.norm_pad(False, npad)
			fftip(ima)
			ima.center_origin_fft()
			img_data.append(ima)
		if CTF:
			if not (ctf_params.defocus in defocus_list):
				defocus_list.append(ctf_params.defocus)
				ctfimg = ctf_img(N, ctf_params, ny = N, nz = 1)
				ctfimg_list.append(ctfimg)
				index_list.append(len(defocus_list)-1)
			else:
				index = defocus_list.index(ctf_params.defocus)
				ctfimg = ctfimg_list[index]
				index_list.append(index)
			Util.add_img2(ctfimg2, ctfimg)			
		
		if opti_method == "LBFGSB":
			bounds.append((x0[im*3]-2.0, x0[im*3]+2.0))
			bounds.append((x0[im*3+1]-1.0, x0[im*3+1]+1.0))
			bounds.append((x0[im*3+2]-1.0, x0[im*3+2]+1.0))

	# Use a gradient method here
	data = []
	data.append(img_data)
	data.append(nima)
	data.append(maskI)
	data.append(kb)
	data.append(_jX)
	data.append(_jY)
	data.append(CTF)
	data.append(SSNR_fit)

	if CTF:
		data.append(index_list)
		data.append(ctfimg_list)
		data.append(ctfimg2)
	if SSNR_fit:
		"""
		SSNR_img = EMData(N, N, 1, False)
		SSNR_sig = EMData(N, N, 1, False)
		for x in xrange((N+2)/2):
			for y in xrange(N):
 				if y > N/2-1: yy = y-N
				else: yy = y
				r = sqrt(x**2+yy**2)
				if r < N*0.49:
					i = int(r/2)
					j = r/2-i
					temp = (SSNR[i]*(1-j)+SSNR[i+1]*j)*nima
					SSNR_img.set_value_at(x*2, y, temp) 
					if temp < 1.0: temp = 1.0
					SSNR_sig.set_value_at(x*2, y, temp)					
		data.append(SSNR_img)
		data.append(SSNR_sig)
		"""
		SSNR_r = []
		for i in xrange(N/2+1):
			if i%2==0: SSNR_r.append(SSNR[i/2]*nima)
			else:	SSNR_r.append((SSNR[i/2]+SSNR[i/2+1])*0.5*nima)
		data.append(SSNR_r)
	
	'''
	from numpy import array
	x0 = array(x0)
	aaa = SSNR_func_MPI(x0, data)
	bbb = SSNR_grad_MPI(x0, data)
	ccc = [0.0]*(len(x0))
	for i in xrange(len(x0)):
		x0[i] += 0.1
		ccc[i] = SSNR_func_MPI(x0, data)
		x0[i] -= 0.1
		mpi_barrier(MPI_COMM_WORLD)
	
	if myid == 0:
		t = 0
		for i in xrange(len(x0)):
			x1 = ccc[i]-aaa
			x2 = bbb[i]*0.1
			if x1!=0.0: pct = abs(x1-x2)/abs(x1)
			else: pct = 0.0
			t += pct	
			print i, x1, x2, pct
		print "Average error = ", t/len(x0)
	exit()
	'''

	if opti_method == "CG":
		ps = fmin_cg(SSNR_func_MPI, x0, fprime=SSNR_grad_MPI, args=([data]), gtol=1e-3, norm=Inf, epsilon=1e-5,
        	      maxiter=max_iter, full_output=0, disp=0, retall=0, callback=None)
	else:		
		ps, val, d = fmin_l_bfgs_b(SSNR_func_MPI, x0, args=[data], fprime=SSNR_grad_MPI, bounds=bounds, m=10, 
			factr=1e3, pgtol=1e-4, epsilon=1e-2, iprint=-1, maxfun=max_iter)
	
	# write the result into the header
	if myid == main_node: 
		for im in xrange(nima):
			ima = EMData()
			ima.read_image(stack, im)
			if mir[im] == 0:
				set_params2D(ima, [ps[im*3], ps[im*3+1], ps[im*3+2], 0, 1.0])
			else:
				set_params2D(ima, [360.0-ps[im*3], -ps[im*3+1], ps[im*3+2], 1, 1.0])
			ima.write_image(stack, im)

	mpi_barrier(MPI_COMM_WORLD)	
	if myid == main_node:	print_end_msg("ali_SSNR_MPI")



def ali3d_eB_MPI_LAST_USED(stack, ref_vol, outdir, maskfile, ou=-1,  delta=2, maxit=10, CTF = None, snr=1.0, sym="c1", chunk = -1.0, user_func_name="ref_aliB_cone"):
	'''
		Cone
	'''
	from alignment	    import proj_ali_incore_cone
	from filter         import filt_ctf, filt_params, filt_table, filt_from_fsc, filt_btwl, filt_tanl
	from fundamentals   import fshift, rot_avg_image
	from projection     import prep_vol, prgs
	from utilities      import model_circle, get_arb_params, set_arb_params, drop_spider_doc
	from utilities      import get_image, drop_image, bcast_EMData_to_all, send_attr_dict, recv_attr_dict
	from utilities      import read_spider_doc, get_im
	from reconstruction import rec3D_MPI
	from statistics     import ccc
	from math           import pi, sqrt
	from string         import replace
	import os
	import sys
	from mpi 	    import mpi_comm_size, mpi_comm_rank, MPI_COMM_WORLD
	from mpi 	    import mpi_barrier

	number_of_proc = mpi_comm_size(MPI_COMM_WORLD)
	myid = mpi_comm_rank(MPI_COMM_WORLD)
	main_node = 0
	if(myid == main_node):
		if os.path.exists(outdir):  os.system('rm -rf '+outdir)
		os.mkdir(outdir)
		import user_functions
		user_func = user_functions.factory[user_func_name]
	mpi_barrier(MPI_COMM_WORLD)

	nima = EMUtil.get_image_count(stack)

	image_start, image_end = MPI_start_end(nima, number_of_proc, myid)

	# figure the size of the chunk (3D is updated after each chunk).  Chunk should be given as 0.0< chunk <= 1.0.  1.0 means all projections
	if(chunk <= 0.0):  chunk = 1.0
	n_in_chunk  = max(int(chunk * (image_end-image_start+1)), 1)
	n_of_chunks = (image_end-image_start+1)//n_in_chunk + min((image_end-image_start+1)%n_in_chunk,1)
	
	outf = file(replace("progress%4d"%myid,' ','0'), "w")
	outf.write("  chunk = "+str(chunk)+"   ")
	outf.write("\n")
	outf.flush()
	outf.write("  n_in_chunk = "+str(n_in_chunk)+"   ")
	outf.write("  n_of_chunks = "+str(n_of_chunks)+"   ")
	outf.write("\n")
	outf.flush()
	#  Here we assume that reference volume exists
	vol = EMData()
	vol.read_image(ref_vol)
	nx  = vol.get_xsize()
	#ima = EMData()
	#ima.read_image(stack)
	#nx  = ima.get_xsize()
	if(ou <= 0):  ou = nx//2-2
	#if(myid == main_node):  vol.write_image(os.path.join(outdir,"ref_volf00.hdf"))
	if maskfile:
		import  types
		if(type(maskfile) is types.StringType):  mask3D = get_image(maskfile)
		else:                                   mask3D = maskfile
	else:
		mask3D = model_circle(ou, nx, nx, nx)
	mask2D = model_circle(ou, nx, nx)



	#from utilities import read_spider_doc, set_arb_params
	#prm = read_spider_doc("params_new.doc")
	#prm_dict = ["phi", "theta", "psi", "s2x", "s2y"]
	#prm_dict = ["phi", "theta", "psi"]
	#from random import seed,gauss
	#seed()
	dataim = EMData.read_images(image_start, image_end)
	for im in xrange(image_start, image_end):
		dataim[im].set_attr('ID', im)
		#angn = get_arb_params(ima, prm_dict)
		#for ian in xrange(3):  angn[ian] += gauss(0.0, 1.0)
		#set_arb_params(ima, angn, prm_dict)
		#set_arb_params(ima, prm[im], prm_dict)
		if(CTF):
			ctf_params = data[im].get_attr('ctf')
			if(im == image_start): data_had_ctf = data[im].get_attr('ctf_applied')
			if data[im].get_attr('ctf_applied') == 0:
				st = Util.infomask(data[im], mask2D, False)
				data[im] -= st[0]
				from filter import filt_ctf
				data[im] = filt_ctf(data[im], ctf_params)
				data[im].set_attr('ctf_applied', 1)

	outf.write("  data read = "+str(image_start)+"   ")
	outf.write("\n")
	outf.flush()


	if (myid == main_node):
		# initialize data for the reference preparation function
		from utilities import read_text_file
		ref_data = []
		ref_data.append( mask3D )
		ref_data.append( read_text_file("pwpdb.txt", 1) )

	from utilities      import bcast_number_to_all
	from morphology     import threshold, threshold_to_minval

	par_str=["phi", "theta", "psi", "s2x", "s2y"]

	#  this is needed for gathering and scattering of cccfs
	disps = []
	recvcount = []
	for im in xrange(number_of_proc):
		if( im == main_node ):  disps.append(0)
		else:                  disps.append(disps[im-1] + recvcount[im-1])
		ib, ie = MPI_start_end(nima, number_of_proc, im)
		recvcount.append( ie - ib )

	from utilities import even_angles
	xrng = 0.0
	yrng = 0.0
	step = 1.0
	template_angles = even_angles(0.2, 0.0, 0.4, phiEqpsi = 'Zero', method = 'P')
	from filter import  filt_tophatb, filt_gaussl
	fifi = True
	for iteration in xrange(maxit):
		outf.write("  iteration = "+str(iteration)+"   ")
		outf.write("\n")
		outf.flush()
		for  ic  in xrange(n_of_chunks):
			# compute updated 3D after each chunk
			if(fifi):
				# resolution
				outf.write("  begin reconstruction = "+str(image_start)+"   ")
				outf.write("\n")
				outf.flush()
				vol, fscc = rec3D_MPI(dataim, snr, sym, mask3D, os.path.join(outdir, replace("resolution%3d_%3d"%(iteration, ic),' ','0') ), myid, main_node)
				outf.write("  done reconstruction = "+str(image_start)+"   ")
				outf.write("\n")
				outf.flush()

				mpi_barrier(MPI_COMM_WORLD)
				if(myid == main_node):
					drop_image(vol, os.path.join(outdir, replace("vol%3d_%3d.hdf"%(iteration, ic),' ','0') ))
					
					ref_data.append( vol )
					ref_data.append( fscc )
					#  call user-supplied function to prepare reference image, i.e., filter it
					vol = user_func( ref_data )
					#  HERE CS SHOULD BE USED TO MODIFY PROJECTIONS' PARAMETERS  !!!
					del ref_data[2]
					del ref_data[2]
					drop_image(vol, os.path.join(outdir, replace("volf%3d_%3d.hdf"%(iteration, ic),' ','0') ))
				#from sys import exit
				#exit()

				bcast_EMData_to_all(vol, myid, main_node)
			fifi = True
			volft,kb  = prep_vol(vol)

			image_start_in_chunk = image_start + ic*n_in_chunk
			image_end_in_chunk   = min(image_start_in_chunk + n_in_chunk, image_end)
			outf.write("ic "+str(ic)+"   image_start "+str(image_start)+"   n_in_chunk "+str(n_in_chunk)+"   image_end "+str(image_end)+"\n")
			outf.write("image_start_in_chunk "+str(image_start_in_chunk)+"  image_end_in_chunk "+str(image_end_in_chunk)+"\n")
			outf.flush()
			for imn in xrange(image_start_in_chunk, image_end_in_chunk):
				#from utilities import start_time, finish_time
				#t3=start_time()
				proj_ali_incore_cone(volft, kb, template_angles, dataim[imn-image_start], 1, ou, 1, xrng, yrng, step, outf)

			soto = []
			for imn in xrange(image_start, image_end):
				from utilities import set_params_proj, get_params_proj
				phi,theta,psi,s2x,s2y = get_params_proj( dataim[imn-image_start] )
				soto.append([phi,theta,psi,s2x,s2y,imn])
			drop_spider_doc(os.path.join(outdir, replace("new_params%3d_%3d_%3d"%(iteration, ic, myid),' ','0')), soto," phi, theta, psi, s2x, s2y, image number")
			del soto

			#from sys import exit
			#exit()
		#  here we should write header info, just in case the program crashes...
	del vol
	del volft
	# write out headers  and STOP, under MPI writing has to be done sequentially
	mpi_barrier(MPI_COMM_WORLD)

def ali3d_eB_CCC(stack, ref_vol, outdir, maskfile, ou=-1,  delta=2, maxit=10, CTF = None, snr=1.0, sym="c1", chunk = -1.0, user_func_name="ref_aliB_cone"):
	'''
		Cone, modified version to test CCC
		single processor version
	'''
	from utilities      import print_begin_msg, print_end_msg, print_msg

	from alignment	    import proj_ali_incore_cone
	from filter         import filt_ctf, filt_params, filt_table, filt_from_fsc, filt_btwl, filt_tanl
	from fundamentals   import fshift, rot_avg_image
	from projection     import prep_vol, prgs
	from utilities      import model_circle, get_arb_params, set_arb_params, drop_spider_doc
	from utilities      import get_image, drop_image, send_attr_dict, recv_attr_dict
	from utilities      import read_spider_doc, get_im
	from statistics     import ccc
	from statistics     import fsc_mask
	from math           import pi, sqrt
	from string         import replace
	import os
	import sys

	number_of_proc = 1
	myid = 0
	main_node = 0
	
	if os.path.exists(outdir):  os.system('rm -rf '+outdir)
	os.mkdir(outdir)
	import user_functions
	user_func = user_functions.factory[user_func_name]

	if CTF :from reconstruction import recons3d_4nn_ctf
	else   : from reconstruction import recons3d_4nn


	nima = EMUtil.get_image_count(stack)

	image_start = 0
	image_end   = nima

	# figure the size of the chunk (3D is updated after each chunk).  Chunk should be given as 0.0< chunk <= 1.0.  1.0 means all projections
	if(chunk <= 0.0):  chunk = 1.0
	n_in_chunk  = max(int(chunk * (image_end-image_start+1)), 1)
	n_of_chunks = (image_end-image_start+1)//n_in_chunk + min((image_end-image_start+1)%n_in_chunk,1)
	
	
	outf = file(os.path.join(outdir, "progress"), "w")
	outf.write("  chunk = "+str(chunk)+"   ")
	outf.write("\n")
	outf.flush()
	outf.write("  n_in_chunk = "+str(n_in_chunk)+"   ")
	outf.write("  n_of_chunks = "+str(n_of_chunks)+"   ")
	outf.write("\n")
	outf.flush()
	#  Here we assume that reference volume exists
	vol = EMData()
	vol.read_image(ref_vol)
	nx  = vol.get_xsize()
	#ima = EMData()
	#ima.read_image(stack)
	#nx  = ima.get_xsize()
	ou = int(ou)
	if(ou <= 0):  ou = nx//2-2
	#if(myid == main_node):  vol.write_image(os.path.join(outdir,"ref_volf00.hdf"))
	if maskfile:
		import  types
		if(type(maskfile) is types.StringType):  mask3D = get_image(maskfile)
		else:                                   mask3D = maskfile
	else:
		mask3D = model_circle(ou, nx, nx, nx)
	mask2D = model_circle(ou, nx, nx)

	dataim = []
	#from utilities import read_spider_doc, set_arb_params
	#prm = read_spider_doc("params_new.doc")
	#prm_dict = ["phi", "theta", "psi", "s2x", "s2y"]
	#prm_dict = ["phi", "theta", "psi"]
	#from random import seed,gauss
	#seed()
	for im in xrange(image_start, image_end):
		ima = EMData()
		ima.read_image(stack, im)
		#angn = get_arb_params(ima, prm_dict)
		#for ian in xrange(3):  angn[ian] += gauss(0.0, 1.0)
		#set_arb_params(ima, angn, prm_dict)
		#set_arb_params(ima, prm[im], prm_dict)
		if(CTF):
			ctf_params = ima.get_attr('ctf')
			if(im == image_start): data_had_ctf = ima.get_attr('ctf_applied')
			if ima.get_attr('ctf_applied') == 0:
				st = Util.infomask(ima, mask2D, False)
				ima -= st[0]
				from filter import filt_ctf
				ima = filt_ctf(ima, ctf_params)
				ima.set_attr('ctf_applied', 1)
		dataim.append(ima)
	outf.write("  data read = "+str(image_start)+"   ")
	outf.write("\n")
	outf.flush()


	# initialize data for the reference preparation function
	from utilities import read_text_file
	ref_data = []
	ref_data.append( mask3D )
	ref_data.append( read_text_file("pwpdb.txt", 1) )
	from utilities import read_text_file
	fscc = [read_text_file("resolution000_000",0), read_text_file("resolution000_000",1)]
 	jtep = 0

	par_str=["phi", "theta", "psi", "s2x", "s2y"]

	from utilities import even_angles
	xrng = 0.0
	yrng = 0.0
	step = 1.0
	template_angles = even_angles(0.2, 0.0, 3.5, phiEqpsi = 'Zero', method = 'P')
	print  len(template_angles)
	for iteration in xrange(maxit):
		msg = "ITERATION #%3d\n"%(iteration+1)
		print_msg(msg)
		for  ic  in xrange(n_of_chunks):
			jtep += 1
			drop_image(vol, os.path.join(outdir, replace("vol%3d_%3d.hdf"%(iteration, ic),' ','0') ))
			ref_data.append( vol )
			ref_data.append( fscc )
			#  call user-supplied function to prepare reference image, i.e., filter it
			vol = user_func( ref_data )
			#  HERE CS SHOULD BE USED TO MODIFY PROJECTIONS' PARAMETERS  !!!
			del ref_data[2]
			del ref_data[2]
			drop_image(vol, os.path.join(outdir, replace("volf%3d_%3d.hdf"%(iteration, ic),' ','0') ))

			volft,kb  = prep_vol(vol)

			image_start_in_chunk = image_start + ic*n_in_chunk
			image_end_in_chunk   = min(image_start_in_chunk + n_in_chunk, image_end)
			outf.write("ic "+str(ic)+"   image_start "+str(image_start)+"   n_in_chunk "+str(n_in_chunk)+"   image_end "+str(image_end)+"\n")
			outf.write("image_start_in_chunk "+str(image_start_in_chunk)+"  image_end_in_chunk "+str(image_end_in_chunk)+"\n")
			outf.flush()
			for imn in xrange(image_start_in_chunk, image_end_in_chunk):
			#for imn in xrange(2,3):
				proj_ali_incore_cone(volft, kb, template_angles, dataim[imn-image_start], 1, ou, 1, xrng, yrng, step, outf)

			soto = []
			for imn in xrange(image_start, image_end):
				from utilities import set_params_proj, get_params_proj
				phi,theta,psi,s2x,s2y = get_params_proj( dataim[imn-image_start] )
				soto.append([phi,theta,psi,s2x,s2y,imn])
			drop_spider_doc(os.path.join(outdir, replace("new_params%3d_%3d"%(iteration, ic),' ','0')), soto," phi, theta, psi, s2x, s2y, image number")
			del soto
			from sys import exit
			exit()

			# compute updated 3D after each chunk
 	    		# resolution
			#print  " start reconstruction",image_start,image_end
			#  3D stuff
			list_p = range(0,nima,2)
 			if(CTF): vol1 = recons3d_4nn_ctf(stack, list_p, snr, 1, sym)
			else:	 vol1 = recons3d_4nn(stack, list_p, sym)

			list_p = range(1,nima,2)
			if(CTF): vol2 = recons3d_4nn_ctf(stack, list_p, snr, 1, sym)
			else:	 vol2 = recons3d_4nn(stack, list_p, sym)

			fscc = fsc_mask(vol1, vol2, mask3D, 1.0, os.path.join(outdir, replace("resolution%4d"%(iteration*n_of_chunks+ic+1),' ','0')))
			del vol1
			del vol2

			# calculate new and improved 3D
			list_p = range(nima)
			if(CTF): vol = recons3d_4nn_ctf(stack, list_p, snr, 1, sym)
			else:	 vol = recons3d_4nn(stack, list_p, sym)
			# store the reference volume
			#drop_image(vol,os.path.join(outdir, replace("vol%4d.spi"%(N_step*max_iter+Iter+1),' ','0')), "s")
			drop_image(vol,os.path.join(outdir, replace("vol%4d.hdf"%(iteration*n_of_chunks+ic+1),' ','0')), "s")

def ali3d_eB_MPI_conewithselect(stack, ref_vol, outdir, maskfile, ou=-1,  delta=2, maxit=10, CTF = None, snr=1.0, sym="c1", chunk = -1.0):
	'''
		Cone
	'''
	from alignment	    import proj_ali_incore_cone
	from filter         import filt_ctf, filt_params, filt_table, filt_from_fsc, filt_btwl, filt_tanl
	from fundamentals   import fshift, rot_avg_image
	from projection     import prep_vol, prgs
	from utilities      import model_circle, get_arb_params, set_arb_params, drop_spider_doc
	from utilities      import get_image, drop_image, bcast_EMData_to_all, send_attr_dict, recv_attr_dict
	from utilities      import read_spider_doc, get_im
	from reconstruction import rec3D_MPI
	from statistics     import ccc
	from math           import pi, sqrt
	from string         import replace
	import os
	import sys
	#from development    import ali_G3
	from mpi 	    import mpi_comm_size, mpi_comm_rank, MPI_COMM_WORLD
	from mpi 	    import mpi_barrier, mpi_gatherv, mpi_scatterv
	from mpi 	    import MPI_FLOAT, MPI_INT, MPI_SUM

	number_of_proc = mpi_comm_size(MPI_COMM_WORLD)
	myid = mpi_comm_rank(MPI_COMM_WORLD)
	main_node = 0
	if(myid == main_node):
		if os.path.exists(outdir):  os.system('rm -rf '+outdir)
		os.mkdir(outdir)
		from utilities import read_text_file
		pwpdb = read_text_file("pwpdb.txt", 1)
	mpi_barrier(MPI_COMM_WORLD)

	nima = EMUtil.get_image_count(stack)

	image_start, image_end = MPI_start_end(nima, number_of_proc, myid)

	# figure the size of the chunk (3D is updated after each chunk).  Chunk should be given as 0.0< chunk <= 1.0.  1.0 means all projections
	if(chunk <= 0.0):  chunk = 1.0
	n_in_chunk  = max(int(chunk * (image_end-image_start+1)), 1)
	n_of_chunks = (image_end-image_start+1)//n_in_chunk + min((image_end-image_start+1)%n_in_chunk,1)
	
	outf = file(replace("progress%4d"%myid,' ','0'), "w")
	outf.write("  chunk = "+str(chunk)+"   ")
	outf.write("\n")
	outf.flush()
	outf.write("  n_in_chunk = "+str(n_in_chunk)+"   ")
	outf.write("  n_of_chunks = "+str(n_of_chunks)+"   ")
	outf.write("\n")
	outf.flush()
	#  Here we assume that reference volume exists
	vol = EMData()
	vol.read_image(ref_vol)
	nx  = vol.get_xsize()
	#ima = EMData()
	#ima.read_image(stack)
	#nx  = ima.get_xsize()
	if(ou <= 0):  ou = nx//2-2
	#if(myid == main_node):  vol.write_image(os.path.join(outdir,"ref_volf00.hdf"))
	if maskfile:
		import  types
		if(type(maskfile) is types.StringType):  mask3D = get_image(maskfile)
		else:                                   mask3D = maskfile
	else:
		mask3D = model_circle(ou, nx, nx, nx)
	mask2D = model_circle(ou, nx, nx)

	dataim = []
	#from utilities import read_spider_doc, set_arb_params
	#prm = read_spider_doc("params_new.doc")
	#prm_dict = ["phi", "theta", "psi", "s2x", "s2y"]
	#prm_dict = ["phi", "theta", "psi"]
	#from random import seed,gauss
	#seed()
	for im in xrange(image_start, image_end):
		ima = EMData()
		ima.read_image(stack, im)
		#angn = get_arb_params(ima, prm_dict)
		#for ian in xrange(3):  angn[ian] += gauss(0.0, 1.0)
		#set_arb_params(ima, angn, prm_dict)
		#set_arb_params(ima, prm[im], prm_dict)
		if(CTF):
			ctf_params = ima.get_attr('ctf')
			if(im == image_start): data_had_ctf = ima.get_attr('ctf_applied')
			if ima.get_attr('ctf_applied') == 0:
				st = Util.infomask(ima, mask2D, False)
				ima -= st[0]
				from filter import filt_ctf
				ima = filt_ctf(ima, ctf_params)
				ima.set_attr('ctf_applied', 1)
		dataim.append(ima)
	outf.write("  data read = "+str(image_start)+"   ")
	outf.write("\n")
	outf.flush()

	from utilities      import bcast_number_to_all
	from morphology     import threshold, threshold_to_minval

	par_str=["phi", "theta", "psi", "s2x", "s2y"]

	#  this is needed for gathering and scattering of cccfs
	disps = []
	recvcount = []
	for im in xrange(number_of_proc):
		if( im == main_node ):  disps.append(0)
		else:                  disps.append(disps[im-1] + recvcount[im-1])
		ib, ie = MPI_start_end(nima, number_of_proc, im)
		recvcount.append( ie - ib )

	from utilities import even_angles
	xrng = 0.0
	yrng = 0.0
	step = 1.0
	template_angles = even_angles(0.2, 0.0, 0.4, phiEqpsi = 'Zero', method = 'P')
	from filter import  filt_tophatb, filt_gaussl
	volftb = vol.copy()
	for iteration in xrange(maxit):
		outf.write("  iteration = "+str(iteration)+"   ")
		outf.write("\n")
		outf.flush()
		for  ic  in xrange(n_of_chunks):
			# compute updated 3D after each chunk
			
			outf.write("  generate projections = "+str(image_start)+"   ")
			outf.write("\n")
			outf.flush()
			# SORT PROJECTIONS
			volftb,kb  = prep_vol( filt_tophatb(volftb, 0.28, 0.44, False) )

			qcc = []
			for imn in xrange(image_start, image_end):
				atparams = get_arb_params(dataim[imn-image_start], par_str)
				projt = prgs(volftb, kb, [atparams[0], atparams[1], atparams[2], -atparams[3], -atparams[4]])
				qcc.append(ccc(projt, dataim[imn-image_start], mask2D))
				#qqcc=ccc(projt, dataim[imn-image_start], mask2D)
				#dataim[imn-image_start] /= (1.0-qqcc*qqcc)
			del projt
			del volftb

			recvbuf = mpi_gatherv(qcc, len(dataim), MPI_FLOAT, recvcount, disps, MPI_FLOAT, main_node, MPI_COMM_WORLD)
			del qcc
			mpi_barrier(MPI_COMM_WORLD)
			if(myid == main_node):
				templ = []
				for im in xrange(len(recvbuf)):    templ.append([float(recvbuf[im]), im])
				del recvbuf
				'''
				for im in xrange(len(templ)):
					outf.write('ccc, image %12.5f  %07d'%( templ[im][0], templ[im][1]  ))
					outf.write("\n")
				outf.flush()
				'''
				
				templ.sort()
				ilow = int(0.25*len(templ))  # reject 25% worst images.
				for im in xrange(ilow):  templ[im] = [templ[im][1], 0]
				for im in xrange(ilow, len(templ)):  templ[im] = [templ[im][1], 1]
				templ.sort()
				sendbuf = []
				for im in xrange(len(templ)):	sendbuf.append(templ[im][1])
				del templ
				'''
				qb = -1.0
				qs = 10.0
				for im in xrange(len(recvbuf)):
					qt = float(recvbuf[im])
					qb = max(qb,qt)
					qs = min(qs,qt)
				qs -= 1.0e-3
				qb -= qs
				sendbuf = []
				for im in xrange(len(recvbuf)):
					sendbuf.append((float(recvbuf[im])-qs)/qb)
				del recvbuf
				'''
			else:
				sendbuf = []
			mpi_barrier(MPI_COMM_WORLD)
			#recvbuf = mpi_scatterv(sendbuf, recvcount, disps, MPI_FLOAT, recvcount[myid], MPI_FLOAT, main_node, MPI_COMM_WORLD)
			recvbuf = mpi_scatterv(sendbuf, recvcount, disps, MPI_INT, recvcount[myid], MPI_INT, main_node, MPI_COMM_WORLD)
			del sendbuf

			# horatio active_refactoring Jy51i1EwmLD4tWZ9_00000_1
			# for imn in xrange(image_start, image_end):
			# 	dataim[imn-image_start].set_attr_dict({'active': int(recvbuf[imn-image_start])})
			# 	#dataim[imn-image_start] /= float(recvbuf[imn-image_start])


			'''
			nact = 0
			for imn in xrange(image_start, image_end):
				nact += dataim[imn-image_start].get_attr('active')
			nact = float(nact)
			tn = mpi_reduce(nact, 1, MPI_FLOAT, MPI_SUM, main_node, MPI_COMM_WORLD)
			if(myid == main_node):
				outf.write('total number of used images %12.2f  '%(float(tn)))
				outf.write("\n")
				outf.flush()
			'''
			# resolution
			outf.write("  begin reconstruction = "+str(image_start)+"   ")
			outf.write("\n")
			outf.flush()
			vol, fscc = rec3D_MPI(dataim, snr, sym, mask3D, os.path.join(outdir, replace("resolution%3d_%3d"%(iteration, ic),' ','0') ), myid, main_node)
			outf.write("  done reconstruction = "+str(image_start)+"   ")
			outf.write("\n")
			outf.flush()
			volftb = vol.copy()

			#  restore original normalization
			#for imn in xrange(image_start, image_end):
			#	#dataim[imn-image_start].set_attr_dict({'active': int(recvbuf[imn-image_start])})
			#	dataim[imn-image_start] *= float(recvbuf[imn-image_start])
			del recvbuf
			mpi_barrier(MPI_COMM_WORLD)
			if(myid == main_node):
				drop_image(vol, os.path.join(outdir, replace("vol%3d_%3d.hdf"%(iteration, ic),' ','0') ))
				stat = Util.infomask(vol, mask3D, False)
				vol -= stat[0]
				vol /= stat[1]
				vol = threshold(vol)
				from  fundamentals  import  rops_table
				pwem = rops_table(vol)
				ftb = []
				for idum in xrange(len(pwem)):
					ftb.append(sqrt(pwpdb[idum]/pwem[idum]))
				from filter import filt_table, fit_tanh
				vol = filt_table(vol, ftb)
				del ftb, stat
				Util.mul_img(vol, mask3D)
				fl, aa = fit_tanh(fscc)
				vol = filt_tanl(vol, fl, aa)
				#vol = filt_gaussl(filt_tanl(vol, fl, aa),  0.2)
				outf.write('tanh params %8.4f  %8.4f '%(fl, aa))
				outf.write("\n")
				outf.flush()
				drop_image(vol, os.path.join(outdir, replace("volf%3d_%3d.hdf"%(iteration, ic),' ','0') ))
			#from sys import exit
			#exit()

			bcast_EMData_to_all(volftb, myid, main_node)
			bcast_EMData_to_all(vol, myid, main_node)
			volft,kb  = prep_vol(vol)

			image_start_in_chunk = image_start + ic*n_in_chunk
			image_end_in_chunk   = min(image_start_in_chunk + n_in_chunk, image_end)
			outf.write("ic "+str(ic)+"   image_start "+str(image_start)+"   n_in_chunk "+str(n_in_chunk)+"   image_end "+str(image_end)+"\n")
			outf.write("image_start_in_chunk "+str(image_start_in_chunk)+"  image_end_in_chunk "+str(image_end_in_chunk)+"\n")
			outf.flush()
			for imn in xrange(image_start_in_chunk, image_end_in_chunk):
				#from utilities import start_time, finish_time
				#t3=start_time()
				proj_ali_incore_cone(volft, kb, template_angles, dataim[imn-image_start], 1, ou, 1, xrng, yrng, step, outf)

			soto = []
			for imn in xrange(image_start, image_end):
				from utilities import set_params_proj, get_params_proj
				phi,theta,psi,s2x,s2y = get_params_proj( dataim[imn-image_start] )
				soto.append([phi,theta,psi,s2x,s2y,imn])
			drop_spider_doc(os.path.join(outdir, replace("new_params%3d_%3d_%3d"%(iteration, ic, myid),' ','0')), soto," phi, theta, psi, s2x, s2y, image number")
			del soto

			#from sys import exit
			#exit()
		#  here we should write header info, just in case the program crashes...
	del vol
	del volft
	# write out headers  and STOP, under MPI writing has to be done sequentially
	mpi_barrier(MPI_COMM_WORLD)


"""
"""
def proj_ali_incore_index(volref, iref, mask3D, projdata, first_ring, last_ring, rstep, xrng, yrng, step, delta, ref_a, symmetry, MPI):
	from utilities    import even_angles, model_circle, compose_transform2, get_params_proj, set_params_proj
	from alignment    import prepare_refprojs
	#  DO NOT USE THIS ONE< WILL BE OBSOLETED SOON  PAP 01/25/08
	mode    = "F"
	# generate list of Eulerian angles for reference projections
	#  phi, theta, psi
	ref_angles = even_angles(delta, symmetry = symmetry, method = ref_a, phiEqpsi = "Minus")
	#from utilities import drop_spider_doc
	#drop_spider_doc("angles.txt",ref_angles)
	#  begin from applying the mask, i.e., subtract the average outside the mask and multiply by the mask
	if(mask3D):
		[mean, sigma, xmin, xmax ] =  Util.infomask(volref, mask3D, False)
		volref -= mean
		Util.mul_img(volref, mask3D)
	#drop_image(volref, "volref.spi", "s")
	#exit()
	nx   = volref.get_xsize()
	ny   = volref.get_ysize()
	#  center is in SPIDER convention
	cnx  = nx//2 + 1
	cny  = ny//2 + 1
	#precalculate rings
	numr = Numrinit(first_ring, last_ring, rstep, mode)
	wr   = ringwe(numr ,mode)

	# prepare 2-D mask for normalization
	mask2D = model_circle(last_ring, nx, ny)
	if(first_ring > 0): mask2D -= model_circle(first_ring, nx, ny)

	# generate reference projections in polar coords
	ref_proj_rings = prepare_refprojs( volref, ref_angles, last_ring, mask2D, cnx, cny, numr, mode, wr, MPI )
	#soto = []
	for imn in xrange(len(projdata)):
		peako = projdata[imn].get_attr('peak')
		from utilities import set_params_proj, get_params_proj
		phi,theta,psi,sxo,syo = get_params_proj( projdata[imn] )
		[ang, sxs, sys, mirror, nref, peak] = Util.multiref_polar_ali_2d(projdata[imn].process("normalize.mask", {"mask":mask2D, "no_sigma":1}), ref_proj_rings, xrng, yrng, step, mode, numr, cnx-sxo, cny-syo)
		if(peak > peako):
			numref=int(nref)
			projdata[imn].set_attr_dict({'peak':peak})
			projdata[imn].set_attr_dict({'group':iref})
			#[ang,sxs,sys,mirror,peak,numref] = apmq(projdata[imn], ref_proj_rings, xrng, yrng, step, mode, numr, cnx-sxo, cny-syo)
			#ang = (ang+360.0)%360.0
			# The ormqip returns parameters such that the transformation is applied first, the mirror operation second.
			#  What that means is that one has to change the the Eulerian angles so they point into mirrored direction: phi+180, 180-theta, 180-psi
			angb, sxb, syb, ct = compose_transform2(0.0, sxs, sys, 1,  -ang, 0.,0.,1)
			if  mirror:
                		phi = (ref_angles[numref][0]+540.0)%360.0
                		theta = 180.0-ref_angles[numref][1]
                		psi = (540.0-ref_angles[numref][2]+angb)%360.0
                		s2x = sxb + sxo
                		s2y = syb + syo
                	else:
                		phi = ref_angles[numref][0]
                		theta = ref_angles[numref][1]
                		psi = (ref_angles[numref][2]+angb+360.0)%360.0
                		s2x = sxb + sxo
                		s2y = syb + syo
			from utilities import set_params_proj, get_params_proj
                	set_params_proj( projdata[imn], [phi, theta, psi, s2x, s2y] )
			#soto.append( [phi, theta,psi,s2x, s2y, mirror, numref, peak] )
	#from utilities import drop_spider_doc
	#drop_spider_doc("ali_s_params.txt",soto)

def proj_ali_incore_localB(volref, mask3D, projdata, first_ring, last_ring, rstep, xrng, yrng, step, delta, an, ref_a, symmetry, info=None, MPI=False):
	#This is for Berlin only
	from utilities    import even_angles, model_circle, compose_transform2, bcast_EMData_to_all
	from alignment    import prepare_refprojs
	from math         import cos, sin, pi
	qv = pi/180.
        
	mode    = "F"
	# generate list of Eulerian angles for reference projections
	#  phi, theta, psi
	ref_angles = even_angles(delta, symmetry = symmetry, method = ref_a, phiEqpsi = "Minus")
	#from utilities import drop_spider_doc
	#drop_spider_doc("angles.txt",ref_angles)
	#  begin from applying the mask, i.e., subtract the average outside the mask and multiply by the mask
	if(mask3D):
		[mean, sigma, xmin, xmax ] =  Util.infomask(volref, mask3D, False)
		volref -= mean
		Util.mul_img(volref, mask3D)
	#drop_image(volref, "volref.spi", "s")
	#exit()
	nx   = volref.get_xsize()
	ny   = volref.get_ysize()
	#  center is in SPIDER convention
	cnx  = nx//2 + 1
	cny  = ny//2 + 1
	#precalculate rings
	numr = Numrinit(first_ring, last_ring, rstep, mode)
	wr   = ringwe(numr ,mode)

	# prepare 2-D mask for normalization
	mask2D = model_circle(last_ring, nx, ny)
	if(first_ring > 0): mask2D -= model_circle(first_ring, nx, ny)

	# generate reference projections in polar coords
	ref_proj_rings = prepare_refprojs( volref, ref_angles, last_ring, mask2D, cnx, cny, numr, mode, wr, MPI )

	for i in xrange(len(ref_angles)):
		n1 = sin(ref_angles[i][1]*qv)*cos(ref_angles[i][0]*qv)
		n2 = sin(ref_angles[i][1]*qv)*sin(ref_angles[i][0]*qv)
		n3 = cos(ref_angles[i][1]*qv)
		ref_proj_rings[i].set_attr_dict( {"n1":n1, "n2":n2, "n3":n3} )

	ant = abs(cos(an*qv))
	for imn in xrange(len(projdata)):
		from utilities import set_params_proj, get_params_proj
		phi,theta,psi,sxo,syo = get_params_proj( projdata[imn] )
		if not(info is None):
			info.write( "prj %4d old params: %8.3f %8.3f %8.3f %8.3f %8.3f\n" %(imn, phi, theta, psi, sxo, syo) )
			info.flush()
		# This is for Berlin only
		from utilities import get_arb_params
		ctf_params = projdata[imn].get_attr('ctf')
		from morphology import ctf_2
		ctf2 = ctf_2(nx, ctf_params)
		nct = len(ctf2)
		from math import exp
		envt = []
		for i in xrange(nct):
			# e(x)*h(x)/(bckg(x)+e(x)**2*ctf(x)**2/20)
			xs = float(i)/2.22/nx
			et = exp(-70.0*xs**2)
			bckgt = exp(-0.8-120.*xs**2)+0.01
			ht = 1.0-0.6*exp(-xs**2/2.0/0.012**2)
			fmt = et/(bckgt + ctf2[i]*et**2/10.0)*ht
			envt.append(fmt)
		from filter import filt_table
		ima = filt_table(projdata[imn], envt)
		#from utilities import drop_spider_doc
		#if(myid == main_node):
		#	if(im == image_start):  drop_spider_doc("matc.doc", envt)
		[ang, sxs, sys, mirror, nref, peak] = Util.multiref_polar_ali_2d_local(ima.process("normalize.mask", {"mask":mask2D, "no_sigma":1}), ref_proj_rings, xrng, yrng, step, ant, mode, numr, cnx-sxo, cny-syo)
		#[ang, sxs, sys, mirror, nref, peak] = Util.multiref_polar_ali_2d_local(projdata[imn].process("normalize.mask", {"mask":mask2D, "no_sigma":1}), ref_proj_rings, xrng, yrng, step, ant, mode, numr, cnx-sxo, cny-syo)
		numref=int(nref)
		#[ang,sxs,sys,mirror,peak,numref] = apmq_local(projdata[imn], ref_proj_rings, xrng, yrng, step, ant, mode, numr, cnx-sxo, cny-syo)
		#ang = (ang+360.0)%360.0
		if(numref > -1):
			# The ormqip returns parameters such that the transformation is applied first, the mirror operation second.
			#  What that means is that one has to change the the Eulerian angles so they point into mirrored direction: phi+180, 180-theta, 180-psi
			angb, sxb, syb, ct = compose_transform2(0.0, sxs, sys, 1,  -ang, 0.,0.,1)
			if  mirror:
				phi = (ref_angles[numref][0]+540.0)%360.0
				theta = 180.0-ref_angles[numref][1]
				psi = (540.0-ref_angles[numref][2]+angb)%360.0
				s2x = sxb + sxo
				s2y = syb + syo
			else:
				phi = ref_angles[numref][0]
				theta = ref_angles[numref][1]
				psi = (ref_angles[numref][2]+angb+360.0)%360.0
				s2x = sxb+sxo
				s2y = syb+syo

			from utilities import set_params_proj, get_params_proj
			set_params_proj( projdata[imn], [phi, theta, psi, s2x, s2y])

			# if -1 local search did not have any neighbors, simply skip it
			if not(info is None):
				info.write( "prj %4d new params: %8.3f %8.3f %8.3f %8.3f %8.3f\n" %(imn, phi, theta, psi, s2x, s2y) )
				info.flush()
	return

def proj_ali_incore_cone(volref, kb, template_angles, projdata, first_ring, last_ring, rstep, xrng, yrng, step, finfo=None):
	#alignment within a cone, no mirror considered

	from utilities    import even_angles, model_circle, compose_transform2, print_msg
	from alignment    import refprojs
	mode    = "F"
	#  Volume is prepared earlier
	nx   = projdata.get_xsize()
	ny   = projdata.get_ysize()
	#  center is in SPIDER convention
	cnx  = nx//2 + 1
	cny  = ny//2 + 1
	#precalculate rings
	numr = Numrinit(first_ring, last_ring, rstep, mode)
	wr   = ringwe(numr ,mode)

	# prepare 2-D mask for normalization
	mask2D = model_circle(last_ring, nx, ny)
	if(first_ring > 0): mask2D -= model_circle(first_ring, nx, ny)
	#  Generate ref_angles according to angles of the current projection

	from utilities import set_params_proj, get_params_proj
	phi, theta, psi, sxo, syo = get_params_proj( projdata )
	R2  = Transform({"type":"spider", "phi":phi,"theta":theta,"psi":0.0})
	ref_angles = []
	for i in xrange(len(template_angles)):
		R1  = Transform({"type":"spider", "phi":template_angles[i][0], "theta":template_angles[i][1], "psi":0.0})
		RR = R1*R2
		Euler = RR.get_rotation("spider")
		ref_angles.append( [ Euler['phi'], Euler['theta'], 0.0])

	# generate reference projections in polar coords
	ref_proj_rings = refprojs( volref, kb, ref_angles, last_ring, mask2D, cnx, cny, numr, mode, wr )

	#if(imn%10 == 0):  print_msg("%d  "%(imn))
	from utilities import set_params_proj, get_params_proj
	phi,theta,psi,sxo,syo = get_params_proj( projdata )
	if not(finfo is None):
		finfo.write( "old params: %8.3f %8.3f %8.3f %8.3f %8.3f\n" %(phi, theta, psi, sxo, syo) )
		finfo.flush()
	[ang, sxs, sys, nref, peak] = Util.multiref_polar_ali_2d_nom(projdata.process("normalize.mask", {"mask":mask2D, "no_sigma":1}), ref_proj_rings, xrng, yrng, step, mode, numr, cnx-sxo, cny-syo)
	numref=int(nref)
	#[ang,sxs,sys,mirror,peak,numref] = apmq(projdata, ref_proj_rings, xrng, yrng, step, mode, numr, cnx-sxo, cny-syo)
	#ang = (ang+360.0)%360.0
	# The ormqip returns parameters such that the transformation is applied first, the mirror operation second.
	#  What that means is that one has to change the the Eulerian angles so they point into mirrored direction: phi+180, 180-theta, 180-psi
	angb, sxb, syb, ct = compose_transform2(0.0, sxs, sys, 1,  -ang, 0.,0.,1)
        phi   = ref_angles[numref][0]
        theta = ref_angles[numref][1]
        psi   = (ref_angles[numref][2]+angb+360.0)%360.0
        s2x   = sxb + sxo
        s2y   = syb + syo
	projdata.set_attr( "peak", peak )

	from utilities import set_params_proj, get_params_proj
	set_params_proj( projdata, [phi, theta, psi, s2x, s2y] )
        if not(finfo is None):
		finfo.write( "new params: %8.3f %8.3f %8.3f %8.3f %8.3f\n" %(phi, theta, psi, s2x, s2y) )
		finfo.flush()

def realid(stack, averages, out_averages, output_dir, ou, xr, ts, maxit, fun, snr, CTF, Fourvar, Ng, num_ali, th_mir, th_err, dst, center_type, CUDA, GPUID, MPI):
	
	if MPI:
		realid_MPI(stack, averages, out_averages, output_dir, ou, xr, ts, maxit, fun, snr, CTF, Fourvar, Ng, num_ali, th_mir, th_err, dst, center_type, CUDA, GPUID)
		return
	from applications import header, cpy
	from utilities    import get_im, file_type, write_header, get_params2D, set_params2D, get_input_from_string
	from statistics   import aves, aves_adw
	from string	  import replace
	from pixel_error  import multi_align_stability
	import os, logging, sys


	number_of_proc = 1
	myid = 0
	main_node = 0

	averages_copy = "averages_org.hdf"
	if myid == main_node:
		if not os.path.exists(output_dir): os.mkdir(output_dir)
		extract_class(stack, averages, output_dir, 'hdf')
		cpy(averages, os.path.join(output_dir, averages_copy))

	os.chdir(output_dir)
	logging.basicConfig(filename = 'main_log_proc_%02d.txt'%(myid), format = '%(asctime)s %(message)s', level = logging.INFO)

	ct = 0
	lrej = []

	K = EMUtil.get_image_count(averages_copy)
	out_averages_myid = replace(out_averages, ".hdf", "_%03d.hdf"%myid)

	if CUDA:
		GPUID = get_input_from_string(GPUID)
		GPUID = map(int, GPUID)
		nGPU = len(GPUID)
		GPUID = GPUID[myid%nGPU]

	for iclass in xrange(K):
		if iclass%number_of_proc == myid:
			logging.info("... Testing the stability of Group %d"%iclass)

			class_name = 'class_%03i.hdf'%iclass

			# Perform reference-free alignment for num_ali times
			all_ali_params = []
			for i in xrange(num_ali):
				header(class_name, 'xform.align2d', rand_alpha=True)

				cmd = "mpirun -np 8 --host node5 --wdir `pwd`  sxali2d.py %s None --MPI --ou=%d --xr='%s' --ts='%s' --maxit=%d --Ng=%d --function=%s --snr=%f --dst=%f --center=%d"%(class_name, ou, xr, ts, maxit, Ng, fun, snr, dst, center_type)
				if CTF: cmd += " --CTF"
				if Fourvar: cmd += " --Fourvar"
				if CUDA:
					cmd += " --CUDA"
					cmd += " --GPUID='%s'"%(GPUID)
				os.system(cmd)

				data = EMData.read_images(class_name)
				ali_params = []
				for im in xrange(len(data)):
					alpha, sx, sy, mirror, scale = get_params2D(data[im])
					ali_params.extend([alpha, sx, sy, mirror])
				all_ali_params.append(ali_params)

			stable_set, mirror_consistent_rate, err = multi_align_stability(all_ali_params, th_mir, th_err, th_err)
			#print "th_mir = ", th_mir
			#print "th_err = ", th_err
			#print "stable_set =", stable_set
			#print mirror_consistent_rate
			#print err
			
			for im in xrange(len(data)):
				set_params2D(data[im], [all_ali_params[0][im*4], all_ali_params[0][im*4+1], all_ali_params[0][im*4+2], all_ali_params[0][im*4+3], 1.0])

			ref = get_im(averages_copy, iclass)
			dic = ref.get_attr_dict()

			if stable_set == []:
				if err == -1:
					logging.info('...... Group rejected due to low mirror cosistent rate %6.3f (threshold = %6.3f)'%(mirror_consistent_rate, th_mir))
				else:
					logging.info('...... mirror consistent rate = %6.3f    error = %6.3f '%(mirror_consistent_rate, err))
					logging.info('...... Group rejected due to high error %6.3f (threshold = %6.3f)'%(err, th_err))
				lrej.extend(dic['members'])
				logging.info('...... All %d particles are rejected and returned to the pool.'%(len(dic['members'])))
			else:
				members = dic['members']
				class_name_stable = "class_%03i_stable.hdf"%iclass

				stable_set_id = []
				for im in xrange(len(stable_set)):  stable_set_id.append(stable_set[im][1])

				stable_members = []
				n_stable = 0
				for im in xrange(len(data)):
					if im in stable_set_id:
						data[im].write_image(class_name_stable, n_stable)
						n_stable += 1
						stable_members.append(members[im])
					else: 
						lrej.append(members[im])
					
				dic['members'] = stable_members
				dic['n_objects'] = len(stable_members)
				dic['err_mir'] = mirror_consistent_rate
				dic['err_pix'] = err
				dic['class_num'] = iclass			

				if CTF:	ave, var = aves_adw(class_name_stable, Ng=-1)
				else:	ave, var = aves(class_name_stable, 'a')

				ave.set_attr_dict(dic)
				ave.write_image(out_averages_myid, ct)
				ct += 1
				logging.info('...... mirror consistent rate = %6.3f    error = %6.3f '%(mirror_consistent_rate, err))
				logging.info('...... Keep this group, but reject %d among %d particles.'%(len(data)-len(stable_set), len(data)))
				os.system('rm -rf '+class_name_stable)

			logging.info('  ')
			os.system('rm -rf '+class_name)

	ln = len(lrej)

	if myid == main_node:
		lrej = map(int, lrej)	
		logging.basicConfig(filename = 'main_log.txt', format = '%(asctime)s %(message)s', level = logging.INFO)
		n_out_avg = 0
		for i in xrange(number_of_proc):
			out_averages_myid = replace(out_averages, ".hdf", "_%03d.hdf"%i)
			if os.path.exists(out_averages_myid):
				data = EMData.read_images(out_averages_myid)
				for im in xrange(len(data)):
					data[im].write_image(out_averages, n_out_avg)
					n_out_avg += 1

		# Set the flag 'active' of rejected images back to 1 
		logging.info('... Keep %d/%d averages, reject %4d images (to be processed again in next iteration). ' % (n_out_avg, K, len(lrej)))
		os.chdir('..')
		if n_out_avg != 0:
			cpy(os.path.join(output_dir, out_averages), out_averages)

		# horatio active_refactoring Jy51i1EwmLD4tWZ9_00000_1
		# if file_type(stack) == 'bdb':
		# 	DB = db_open_dict(stack)
		# 	for idd in lrej: DB.set_attr(idd, 'active', 1)
		# 	DB.close()
		# else:
		# 	im = EMData()
		# 	for idd in lrej:
		# 		im.read_image(stack, idd)
		# 		im.set_attr('active', 1)
		# 		write_header(stack, im, idd)
	
def realid_MPI(stack, averages, out_averages, output_dir, ou, xr, ts, maxit, fun, snr, CTF, Fourvar, Ng, num_ali, th_mir, th_err, dst, center_type, CUDA, GPUID):

	from applications import header, cpy
	from utilities    import get_im, file_type, write_header, get_params2D, set_params2D, get_input_from_string
	from statistics   import aves, aves_adw
	from string	  import replace
	from pixel_error  import multi_align_stability
	import os, logging, sys
	from subprocess   import call

	from mpi import mpi_comm_size, mpi_comm_rank, MPI_COMM_WORLD
	from mpi import mpi_send, mpi_recv, mpi_gatherv, mpi_bcast, mpi_barrier, MPI_INT

	number_of_proc = mpi_comm_size(MPI_COMM_WORLD)
	myid = mpi_comm_rank(MPI_COMM_WORLD)
	main_node = 0
	mpi_barrier(MPI_COMM_WORLD)

	averages_copy = "averages_org.hdf"
	if myid == main_node:
		if not os.path.exists(output_dir): os.mkdir(output_dir)
		extract_class(stack, averages, output_dir, 'hdf')
		cpy(averages, os.path.join(output_dir, averages_copy))
	
	mpi_barrier(MPI_COMM_WORLD)
	os.chdir(output_dir)
	logging.basicConfig(filename = 'main_log_proc_%02d.txt'%(myid), format = '%(asctime)s %(message)s', level = logging.INFO)

	ct = 0
	lrej = []

	K = EMUtil.get_image_count(averages_copy)
	out_averages_myid = replace(out_averages, ".hdf", "_%03d.hdf"%myid)

	if CUDA:
		GPUID = get_input_from_string(GPUID)
		GPUID = map(int, GPUID)
		nGPU = len(GPUID)
		GPUID = GPUID[myid%nGPU]

	for iclass in xrange(K):
		if iclass%number_of_proc == myid:
			logging.info("... Testing the stability of Group %d"%iclass)

			class_name = 'class_%03i.hdf'%iclass

			# Perform reference-free alignment for num_ali times
			'''
			all_ali_params = []
			for i in xrange(num_ali):
				header(class_name, 'xform.align2d', rand_alpha=True)

				cmd = "sxali2d.py %s None --ou=%d --xr='%s' --ts='%s' --maxit=%d --Ng=%d --function=%s --snr=%f --dst=%f --center=%d"%(class_name, ou, xr, ts, maxit, Ng, fun, snr, dst, center_type)
				if CTF: cmd += " --CTF"
				if Fourvar: cmd += " --Fourvar"
				if CUDA:
					cmd += " --CUDA"
					cmd += " --GPUID='%s'"%(GPUID)
				os.system(cmd)

				data = EMData.read_images(class_name)
				ali_params = []
				for im in xrange(len(data)):
					alpha, sx, sy, mirror, scale = get_params2D(data[im])
					ali_params.extend([alpha, sx, sy, mirror])
				all_ali_params.append(ali_params)
			'''
			cmd = ['sxmulti_ali2d.py', class_name, 'None', '--ou=%d'%ou, '--xr=%s'%xr, '--ts=%s'%ts, '--maxit=%d'%maxit, '--Ng=%d'%Ng, 
				'--num_ali=%d'%num_ali, '--function=%s'%fun, '--snr=%f'%snr, '--dst=%f'%dst, '--center=%d'%center_type]
			if CTF: cmd.append('--CTF')
			if Fourvar: cmd.append('--Fourvar')
			if CUDA:
				cmd.append('--CUDA')
				cmd.append('--GPUID=%s'%GPUID)
			status = call(cmd)
			if status!=0: exit()

			data = EMData.read_images(class_name)
			call(['rm', '-f', class_name])

			all_ali_params = []
			for i in xrange(num_ali):
				ali_params = []
				for im in xrange(len(data)):
					alpha, sx, sy, mirror, scale = get_params2D(data[im], "xform.align2d_%02d"%i)
					ali_params.extend([alpha, sx, sy, mirror])
				all_ali_params.append(ali_params)

			stable_set, mirror_consistent_rate, err = multi_align_stability(all_ali_params, th_mir, th_err, th_err)
			#print "th_mir = ", th_mir
			#print "th_err = ", th_err
			#print "stable_set =", stable_set
			#print mirror_consistent_rate
			#print err
			
			for im in xrange(len(data)):
				set_params2D(data[im], [all_ali_params[0][im*4], all_ali_params[0][im*4+1], all_ali_params[0][im*4+2], all_ali_params[0][im*4+3], 1.0])

			ref = get_im(averages_copy, iclass)
			dic = ref.get_attr_dict()

			if stable_set == []:
				if err == -1:
					logging.info('...... Group rejected due to low mirror cosistent rate %6.3f (threshold = %6.3f)'%(mirror_consistent_rate, th_mir))
				else:
					logging.info('...... mirror consistent rate = %6.3f    error = %6.3f '%(mirror_consistent_rate, err))
					logging.info('...... Group rejected due to high error %6.3f (threshold = %6.3f)'%(err, th_err))
				lrej.extend(dic['members'])
				logging.info('...... All %d particles are rejected and returned to the pool.'%(len(dic['members'])))
			else:
				members = dic['members']
				class_name_stable = "class_%03i_stable.hdf"%iclass

				stable_set_id = []
				for im in xrange(len(stable_set)):  stable_set_id.append(stable_set[im][1])

				stable_members = []
				n_stable = 0
				for im in xrange(len(data)):
					if im in stable_set_id:
						data[im].write_image(class_name_stable, n_stable)
						n_stable += 1
						stable_members.append(members[im])
					else: 
						lrej.append(members[im])
					
				dic['members'] = stable_members
				dic['n_objects'] = len(stable_members)
				dic['err_mir'] = mirror_consistent_rate
				dic['err_pix'] = err
				dic['class_num'] = iclass			

				if CTF:	ave, var = aves_adw(class_name_stable, Ng=-1)
				else:	ave, var = aves(class_name_stable, 'a')

				ave.set_attr_dict(dic)
				ave.write_image(out_averages_myid, ct)
				ct += 1
				logging.info('...... mirror consistent rate = %6.3f    error = %6.3f '%(mirror_consistent_rate, err))
				logging.info('...... Keep this group, but reject %d among %d particles.'%(len(data)-len(stable_set), len(data)))
				os.system('rm -rf '+class_name_stable)

			logging.info('  ')
			os.system('rm -rf '+class_name)

	# Bring all rejected particles together
	mpi_barrier(MPI_COMM_WORLD)
	ln = len(lrej)
	if myid == main_node:
		recvcount = [ln]
		disp = []
		for n in xrange(number_of_proc):
			if n != main_node:
				lnn = mpi_recv(1, MPI_INT, n, n*100, MPI_COMM_WORLD)
				lnn = int(lnn[0])
				recvcount.append(lnn)
			if n == 0: disp.append(0)
			else: disp.append(disp[n-1]+recvcount[n-1])
	else:
		mpi_send(ln, 1, MPI_INT, main_node, myid*100, MPI_COMM_WORLD)
		recvcount = [0.0]*number_of_proc
		disp = [0.0]*number_of_proc
	recvcount = mpi_bcast(recvcount, number_of_proc, MPI_INT, main_node, MPI_COMM_WORLD)
	disp = mpi_bcast(disp, number_of_proc, MPI_INT, main_node, MPI_COMM_WORLD)
	recvcount = map(int, recvcount)
	disp = map(int, disp)
	lrej = mpi_gatherv(lrej, ln, MPI_INT, recvcount, disp, MPI_INT, main_node, MPI_COMM_WORLD)

	if myid == main_node:
		lrej = map(int, lrej)		
		logging.basicConfig(filename = 'main_log.txt', format = '%(asctime)s %(message)s', level = logging.INFO)
		n_out_avg = 0
		for i in xrange(number_of_proc):
			out_averages_myid = replace(out_averages, ".hdf", "_%03d.hdf"%i)
			if os.path.exists(out_averages_myid):
				data = EMData.read_images(out_averages_myid)
				for im in xrange(len(data)):
					data[im].write_image(out_averages, n_out_avg)
					n_out_avg += 1

		# Set the flag 'active' of rejected images back to 1 
		logging.info('... Keep %d/%d averages, reject %4d images (to be processed again in next iteration). ' % (n_out_avg, K, len(lrej)))
		os.chdir('..')
		if n_out_avg != 0:
			cpy(os.path.join(output_dir, out_averages), out_averages)

		# horatio active_refactoring Jy51i1EwmLD4tWZ9_00000_1
		# if file_type(stack) == 'bdb':
		# 	DB = db_open_dict(stack)
		# 	for idd in lrej: DB.set_attr(idd, 'active', 1)
		# 	DB.close()
		# else:
		# 	im = EMData()
		# 	for idd in lrej:
		# 		im.read_image(stack, idd)
		# 		im.set_attr('active', 1)
		# 		write_header(stack, im, idd)
	mpi_barrier(MPI_COMM_WORLD)
# ----------------------------------------------------------------------------------------------



def ali3d_n(stack, ref_vol, outdir, maskfile = None, ir = 1, ou = -1, rs = 1, 
            xr = "4 2 2 1", yr = "-1", ts = "1 1 0.5 0.25", delta = "10 6 4 4", an = "-1", deltapsi = "-1", startpsi = "-1",
	    center = -1, maxit = 5, CTF = False, snr = 1.0,  ref_a = "S", sym = "c1",
	    user_func_name = "ref_ali3d", fourvar = True, npad = 4, debug = False, MPI = False, termprec = 0.0, chunk = 0.2):
	"""
		Name
			ali3d - Perform 3-D projection matching given initial reference volume and image series
		Input
			stack: set of 2-D images in a stack file, images have to be squares
			ref_vol: initial reference volume
			outdir: directory name into which the results will be written
			maskfile: filename of the file containing 3D mask.
			ir: inner radius for rotational correlation > 0 
			ou: outer radius for rotational correlation <int(nx/2)-1 
			rs: steps between rings in rotational correlation >0
			xr: range for translation search in x direction in each iteration, search is +/xr
			yr: range for translation search in y direction in each iteration, search is +/yr
			ts: step size of the translation search in both directions, search is -xr, -xr+ts, 0, xr-ts, xr, can be fractional.
			delta: angular step for the reference projections in respective iterations
			an: angular neighborhood for local searches
			center: average center method
			max_iter: maximum iterations at each angle step
			CTF: if the flag is present, program will use the CTF information stored in file headers
			snr: signal noise ratio used in the 3D reconstruction
			ref_a: method for creating quasi-uniform distribution of the projection directions of reference projections: "S" - spiral
			sym: symmetry of the refined structure
			function: name of the user-supplied-function
			MPI: if presetm use MPI version
		Output
			output_directory: directory name into which the output files will be written.
			header: the alignment parameters are stored in the headers of input files as Transform Object xform.proj
	"""
	if MPI:
		from development import ali3d_new_MPI
		ali3d_new_MPI(stack, ref_vol, outdir, maskfile, ir, ou, rs, xr, yr, ts,
	        	delta, an, deltapsi, startpsi, center, maxit, CTF, snr, ref_a, sym, user_func_name,
			fourvar, npad, debug, termprec, chunk)
		return

	from alignment      import proj_ali_incore, proj_ali_incore_local
	from utilities      import model_circle, drop_image, get_image, get_input_from_string
	from utilities      import get_params_proj
	from utilities      import estimate_3D_center, rotate_3D_shift
	from filter         import filt_params, fit_tanh, filt_tanl, filt_ctf
	from statistics     import fsc_mask
	import os
	import types
	from utilities      import print_begin_msg, print_end_msg, print_msg
	print_begin_msg("ali3d")

	from alignment      import Numrinit, prepare_refrings
	from projection     import prep_vol

	import user_functions
	user_func = user_functions.factory[user_func_name]

	if os.path.exists(outdir):  ERROR('Output directory exists, please change the name and restart the program', "ali3d", 1)
	os.mkdir(outdir)
	import global_def
	global_def.LOGFILE =  os.path.join(outdir, global_def.LOGFILE)
	print_begin_msg("ali3d")

	xrng        = get_input_from_string(xr)
	if  yr == "-1":  yrng = xrng
	else          :  yrng = get_input_from_string(yr)
	step        = get_input_from_string(ts)
	delta       = get_input_from_string(delta)
	lstp = min(len(xrng), len(yrng), len(step), len(delta))
	if an == "-1":
		an = [-1] * lstp
	else:
		an = get_input_from_string(an)
	first_ring  = int(ir)
	rstep       = int(rs)
	last_ring   = int(ou)
	max_iter    = int(maxit)
	center      = int(center)

	print_msg("Input stack                 : %s\n"%(stack))
	print_msg("Reference volume            : %s\n"%(ref_vol))	
	print_msg("Output directory            : %s\n"%(outdir))
	print_msg("Maskfile                    : %s\n"%(maskfile))
	print_msg("Inner radius                : %i\n"%(first_ring))

	vol     = EMData()
	vol.read_image(ref_vol)
	nx      = vol.get_xsize()
	if last_ring == -1:	last_ring = nx/2 - 2

	print_msg("Outer radius                : %i\n"%(last_ring))
	print_msg("Ring step                   : %i\n"%(rstep))
	print_msg("X search range              : %s\n"%(xrng))
	print_msg("Y search range              : %s\n"%(yrng))
	print_msg("Translational step          : %s\n"%(step))
	print_msg("Angular step                : %s\n"%(delta))
	print_msg("Angular search range        : %s\n"%(an))
	print_msg("Maximum iteration           : %i\n"%(max_iter))
	print_msg("Center type                 : %i\n"%(center))
	print_msg("CTF correction              : %s\n"%(CTF))
	print_msg("Signal-to-Noise Ratio       : %f\n"%(snr))
	print_msg("Reference projection method : %s\n"%(ref_a))
	print_msg("Symmetry group              : %s\n\n"%(sym))

	if maskfile :
		if type(maskfile) is types.StringType: mask3D = get_image(maskfile)
		else                                  : mask3D = maskfile
	else          :   mask3D = model_circle(last_ring, nx, nx, nx)
	mask2D = model_circle(last_ring, nx, nx) - model_circle(first_ring, nx, nx)
	numr   = Numrinit(first_ring, last_ring, rstep, "F")

	if CTF:
		from reconstruction import recons3d_4nn_ctf
		from filter         import filt_ctf
	else: from reconstruction import recons3d_4nn

	if debug:  outf = file(os.path.join(outdir, "progress"), "w")
	else:      outf = None

	# horatio active_refactoring Jy51i1EwmLD4tWZ9_00000_1
	# active = EMUtil.get_all_attributes(stack, 'active')
	# list_of_particles = []
	# for im in xrange(len(active)):
	# 	if(active[im]):  list_of_particles.append(im)
	# del active

	nima = EMUtil.get_image_count(stack)	
	list_of_particles = range(nima)
	
	data = EMData.read_images(stack, list_of_particles)
	for im in xrange(len(data)):
		data[im].set_attr('ID', list_of_particles[im])
		if CTF:
			ctf_params = data[im].get_attr("ctf")
			st = Util.infomask(data[im], mask2D, False)
			data[im] -= st[0]
			data[im] = filt_ctf(data[im], ctf_params)
			data[im].set_attr('ctf_applied', 1)

	nima = len(data)
	# initialize data for the reference preparation function
	ref_data = [ mask3D, max(center,0), None, None ]#  for center -1 switch of centering by user function

	cs = [0.0]*3
	# do the projection matching
	for N_step in xrange(lstp):
 		for Iter in xrange(max_iter):
			print_msg("\nITERATION #%3d\n"%(N_step*max_iter+Iter+1))

			volft,kb = prep_vol( vol )
			refrings = prepare_refrings( volft, kb, nx, delta[N_step], ref_a, sym, numr, MPI=False)
			del volft,kb

			for im in xrange( nima ):

				if an[N_step] == -1:	
					peak, pixel_error = proj_ali_incore(data[im],refrings,numr,xrng[N_step],yrng[N_step],step[N_step])
				else:
					peak, pixel_error = proj_ali_incore_local(data[im],refrings,numr,xrng[N_step],yrng[N_step],step[N_step],an[N_step])
			if center == -1:
				cs[0], cs[1], cs[2], dummy, dummy = estimate_3D_center(data)
				msg = "Average center x = %10.3f        Center y = %10.3f        Center z = %10.3f\n"%(cs[0], cs[1], cs[2])
				print_msg(msg)
				rotate_3D_shift(data, [-cs[0], -cs[1], -cs[2]])

			if CTF:   vol1 = recons3d_4nn_ctf(data, range(0, nima, 2), snr, 1, sym)
			else:	   vol1 = recons3d_4nn(data, range(0, nima, 2), sym)
			if CTF:   vol2 = recons3d_4nn_ctf(data, range(1, nima, 2), snr, 1, sym)
			else:	   vol2 = recons3d_4nn(data, range(1, nima, 2), sym)

			fscc = fsc_mask(vol1, vol2, mask3D, 1.0, os.path.join(outdir, "resolution%04d"%(N_step*max_iter+Iter+1)))
			del vol1
			del vol2

			# calculate new and improved 3D
			if CTF:  vol = recons3d_4nn_ctf(data, range(nima), snr, 1, sym)
			else:	 vol = recons3d_4nn(data, range(nima), sym)
			# store the reference volume
			drop_image(vol, os.path.join(outdir, "vol%04d.hdf"%(N_step*max_iter+Iter+1)))
			ref_data[2] = vol
			ref_data[3] = fscc

			#  call user-supplied function to prepare reference image, i.e., center and filter it
			vol, dummy = user_func(ref_data)

			drop_image(vol, os.path.join(outdir, "volf%04d.hdf"%(N_step*max_iter+Iter+1)))
			#  here we write header info
			from utilities import write_headers
			#from utilities import write_select_headers
			if CTF:
				for dat in data:  dat.set_attr('ctf_applied',0)
			write_headers(stack, data, list_of_particles)
			#list_params= ['ID','xform.projection']
			#write_select_headers(stack, data, list_params)
			if CTF:
				for dat in data:  dat.set_attr('ctf_applied', 1)
	print_end_msg("ali3d")

def ali3d_new_MPI(stack, ref_vol, outdir, maskfile = None, ir = 1, ou = -1, rs = 1, 
            xr = "4 2 2 1", yr = "-1", ts = "1 1 0.5 0.25", delta = "10 6 4 4", an = "-1", deltapsi = "-1", startpsi = "-1",
	    center = -1, maxit = 5, CTF = False, snr = 1.0,  ref_a = "S", sym = "c1",  user_func_name = "ref_ali3d",
	    fourvar = True, npad = 4, debug = False, termprec = 0.0, chunk = 0.2):

	from alignment       import Numrinit, prepare_refrings, proj_ali_incore, proj_ali_incore_local
	from utilities       import model_circle, get_image, drop_image, get_input_from_string
	from utilities       import bcast_list_to_all, bcast_number_to_all, reduce_EMData_to_root, bcast_EMData_to_all
	from utilities       import send_attr_dict
	from utilities       import get_params_proj, file_type
	from fundamentals    import rot_avg_image
	import os
	import types
	from utilities       import print_begin_msg, print_end_msg, print_msg
	from mpi             import mpi_bcast, mpi_comm_size, mpi_comm_rank, MPI_FLOAT, MPI_COMM_WORLD, mpi_barrier
	from mpi             import mpi_reduce, MPI_INT, MPI_SUM
	from filter          import filt_ctf
	from projection      import prep_vol, prgs
	from statistics      import hist_list, varf3d_MPI
	from random import shuffle
	from applications import MPI_start_end


	number_of_proc = mpi_comm_size(MPI_COMM_WORLD)
	myid           = mpi_comm_rank(MPI_COMM_WORLD)
	main_node = 0
	
	if os.path.exists(outdir):  ERROR('Output directory exists, please change the name and restart the program', "ali3d_MPI", 1, myid)
	mpi_barrier(MPI_COMM_WORLD)
	
	if myid == main_node:
		os.mkdir(outdir)
		import global_def
		global_def.LOGFILE =  os.path.join(outdir, global_def.LOGFILE)
		print_begin_msg("ali3d_MPI")
	mpi_barrier(MPI_COMM_WORLD)

	if debug:
		from time import sleep
		while not os.path.exists(outdir):
			print  "Node ",myid,"  waiting..."
			sleep(5)

		info_file = os.path.join(outdir, "progress%04d"%myid)
		finfo = open(info_file, 'w')
	else:
		finfo = None

	xrng        = get_input_from_string(xr)
	if  yr == "-1":  yrng = xrng
	else          :  yrng = get_input_from_string(yr)
	step        = get_input_from_string(ts)
	delta       = get_input_from_string(delta)
	lstp = min(len(xrng), len(yrng), len(step), len(delta))
	if an == "-1":
		an = [-1] * lstp
	else:
		an = get_input_from_string(an)

	if deltapsi == "-1":
		deltapsi = [-1] * lstp
	else:
		deltapsi = get_input_from_string(deltapsi)

	if startpsi == "-1":
		startpsi = [-1] * lstp
	else:
		startpsi = get_input_from_string(startpsi)

	first_ring  = int(ir)
	rstep       = int(rs)
	last_ring   = int(ou)
	max_iter    = int(maxit)
	center      = int(center)

	vol     = EMData()
	vol.read_image(ref_vol)
	nx      = vol.get_xsize()
	
	if last_ring < 0:	last_ring = int(nx/2) - 2

	if myid == main_node:
		import user_functions
		user_func = user_functions.factory[user_func_name]

		print_msg("Input stack                 : %s\n"%(stack))
		print_msg("Reference volume            : %s\n"%(ref_vol))	
		print_msg("Output directory            : %s\n"%(outdir))
		print_msg("Maskfile                    : %s\n"%(maskfile))
		print_msg("Inner radius                : %i\n"%(first_ring))
		print_msg("Outer radius                : %i\n"%(last_ring))
		print_msg("Ring step                   : %i\n"%(rstep))
		print_msg("X search range              : %s\n"%(xrng))
		print_msg("Y search range              : %s\n"%(yrng))
		print_msg("Translational step          : %s\n"%(step))
		print_msg("Angular step                : %s\n"%(delta))
		print_msg("Angular search range        : %s\n"%(an))
		print_msg("Delta psi                   : %s\n"%(deltapsi))
		print_msg("Start psi                   : %s\n"%(startpsi))
		print_msg("Maximum iteration           : %i\n"%(max_iter))
		print_msg("Percentage of change for termination: %f\n"%(termprec))
		print_msg("Center type                 : %i\n"%(center))
		print_msg("CTF correction              : %s\n"%(CTF))
		print_msg("Signal-to-Noise Ratio       : %f\n"%(snr))
		print_msg("Reference projection method : %s\n"%(ref_a))
		print_msg("Symmetry group              : %s\n\n"%(sym))
		print_msg("chunk                       : %f\n"%(chunk))

	if maskfile:
		if type(maskfile) is types.StringType: mask3D = get_image(maskfile)
		else:                                  mask3D = maskfile
	else: mask3D = model_circle(last_ring, nx, nx, nx)

	numr	= Numrinit(first_ring, last_ring, rstep, "F")
	mask2D  = model_circle(last_ring,nx,nx) - model_circle(first_ring,nx,nx)

	fscmask = mask3D  #model_circle(last_ring,nx,nx,nx)  For a fancy mask circle would work better  PAP 7/21/11
	if CTF:
		from reconstruction import rec3D_MPI
		from filter         import filt_ctf
	else:	 from reconstruction import rec3D_MPI_noCTF

	if myid == main_node:
		# horatio active_refactoring Jy51i1EwmLD4tWZ9_00000_1
		# if(file_type(stack) == "bdb"):
		# 	from EMAN2db import db_open_dict
		# 	dummy = db_open_dict(stack, True)
		# active = EMUtil.get_all_attributes(stack, 'active')
		# list_of_particles = []
		# for im in xrange(len(active)):
		# 	if active[im]:  list_of_particles.append(im)
		# del active
		
		nima = EMUtil.get_image_count(stack)
		list_of_particles = range(nima)
	else:
		nima = 0
	total_nima = bcast_number_to_all(nima, source_node = main_node)

	if myid != main_node:
		list_of_particles = [-1]*total_nima
	list_of_particles = bcast_list_to_all(list_of_particles, source_node = main_node)

	image_start, image_end = MPI_start_end(total_nima, number_of_proc, myid)
	# create a list of images for each node
	list_of_particles = list_of_particles[image_start: image_end]
	nima = len(list_of_particles)
	if debug:
		finfo.write("image_start, image_end: %d %d\n" %(image_start, image_end))
		finfo.flush()

	data = EMData.read_images(stack, list_of_particles)
	if fourvar:  original_data = []
	for im in xrange(nima):
		data[im].set_attr('ID', list_of_particles[im])
		if fourvar: original_data.append(data[im].copy())
		if CTF:
			ctf_params = data[im].get_attr("ctf")
			st = Util.infomask(data[im], mask2D, False)
			data[im] -= st[0]
			data[im] = filt_ctf(data[im], ctf_params)
			data[im].set_attr('ctf_applied', 1)

	if debug:
		finfo.write( '%d loaded  \n' % nima )
		finfo.flush()
	if myid == main_node:
		# initialize data for the reference preparation function
		ref_data = [ mask3D, max(center,0), None, None, None, None ]
		# for method -1, switch off centering in user function

	from time import time	

	#  this is needed for determining chunk, each chunk data used for aligment--anima, data used for reconstruction--rnima
	nchunk = int(1.0/chunk+0.5) 
	chunk_list=[0]*2*nchunk
	anima = [0]*nchunk
	rnima =[0]*nchunk
	for i in xrange(nchunk):
		chunk_list[2*i+0] = int(round(float(nima)/nchunk*i))
		chunk_list[2*i+1]   = int(round(float(nima)/nchunk*(i+1)))
		anima[i] = chunk_list[2*i+1]-chunk_list[2*i+0]
		rnima[i] = nima - anima[i]

	#  this is needed for gathering of pixel errors
	disps = []
	recvcount = []
	for im in xrange(number_of_proc):
		if( im == main_node ):  disps.append(0)
		else:                  disps.append(disps[im-1] + recvcount[im-1])
		ib, ie = MPI_start_end(total_nima, number_of_proc, im)
		recvcount.append( ie - ib )

	pixer = [0.0]*nima
	cs = [0.0]*3
	total_iter = 0
	# do the projection matching
	for N_step in xrange(lstp):
		terminate = 0
		Iter = -1
 		while(Iter < max_iter-1 and terminate == 0):
			Iter += 1
			total_iter += 1
			if myid == main_node:
				start_time = time()
				print_msg("\nITERATION #%3d,  inner iteration #%3d\nDelta = %4.1f, an = %5.2f, xrange = %5.2f, yrange = %5.2f, step = %5.2f, delta psi = %5.2f, start psi = %5.2f\n"%(total_iter, Iter, delta[N_step], an[N_step], xrng[N_step],yrng[N_step],step[N_step],deltapsi[N_step],startpsi[N_step]))
			
			c1 = range(nima)
			shuffle(c1)
			c1_set = set(c1)
			for nch in xrange(nchunk):
				# data is all the input data of each processor, len(data) = nima
				# rdata is data used for reconstruction, len(rdata) = rnima
				# adata is the data used for aligment, len(data) = anima
				list1 = []
				adata = [None]* anima[nch]
				rdata = [None]*rnima[nch]
	
				for k in xrange( len(adata) ):
					adata[k] = data[c1[chunk_list[2*nch]+k]]
					list1.append( c1[chunk_list[2*nch]+k] )
	
				list1_set = set( list1 )
				list2 = list(  c1_set.difference(list1_set) )
				for k in xrange( len(rdata) ):   
					rdata[k]= data[ list2[k] ]
				if myid == main_node:
					start_time = time()
					
									
				if CTF: vol, fscc = rec3D_MPI(rdata, snr, sym, fscmask, os.path.join(outdir, "resolution%04d_%04d"%(total_iter,nch)), myid, main_node, npad = npad)
				else:    vol, fscc = rec3D_MPI_noCTF(rdata, sym, fscmask, os.path.join(outdir, "resolution%04d_%4d"%(total_iter, nch)), myid, main_node, npad = npad)
				del rdata,list2
								
				if myid == main_node:
					ref_data[2] = vol
					ref_data[3] = fscc
					ref_data[4] = None
					#  call user-supplied function to prepare reference image, i.e., center and filter it
					vol, cs = user_func(ref_data)
				
				bcast_EMData_to_all(vol, myid, main_node)
				# write out headers, under MPI writing has to be done sequentially
				mpi_barrier(MPI_COMM_WORLD)
					
				volft,kb = prep_vol( vol )
				refrings = prepare_refrings( volft, kb, nx, delta[N_step], ref_a, sym, numr, True)
				del volft,kb
				
				for im in xrange( anima[nch] ):

					if deltapsi[N_step] > 0.0:
						from alignment import proj_ali_incore_delta
						peak, pixer[list1[im]] = proj_ali_incore_delta(adata[im],refrings,numr,xrng[N_step],yrng[N_step],step[N_step],startpsi[N_step],deltapsi[N_step],finfo)						
					elif an[N_step] == -1:
						peak, pixer[list1[im]] = proj_ali_incore(adata[im],refrings,numr,xrng[N_step],yrng[N_step],step[N_step],finfo)
					else:
						peak, pixer[list1[im]] = proj_ali_incore_local(adata[im],refrings,numr,xrng[N_step],yrng[N_step],step[N_step],an[N_step],finfo)

				del adata,list1,list1_set
			del c1, c1_set
			#output pixel errors
			from mpi import mpi_gatherv
			recvbuf = mpi_gatherv(pixer, nima, MPI_FLOAT, recvcount, disps, MPI_FLOAT, main_node, MPI_COMM_WORLD)
			mpi_barrier(MPI_COMM_WORLD)
			terminate = 0
			if(myid == main_node):
				recvbuf = map(float, recvbuf)
				from statistics import hist_list
				lhist = 20
				region, histo = hist_list(recvbuf, lhist)
				if(region[0] < 0.0):  region[0] = 0.0
				msg = "      Histogram of pixel errors\n      ERROR       number of particles\n"
				print_msg(msg)
				for lhx in xrange(lhist):
					msg = " %10.3f     %7d\n"%(region[lhx], histo[lhx])
					print_msg(msg)
				# Terminate if 95% within 1 pixel error
				im = 0
				for lhx in xrange(lhist):
					if(region[lhx] > 1.0): break
					im += histo[lhx]
				precn = 100*float(total_nima-im)/float(total_nima)
				msg = " Number of particles that changed orientations %7d, percentage of total: %5.1f\n"%(total_nima-im, precn)
				print_msg(msg)
				if(precn <= termprec):  terminate = 1
				del region, histo
			del recvbuf
			terminate = mpi_bcast(terminate, 1, MPI_INT, 0, MPI_COMM_WORLD)
			terminate = int(terminate[0])

			if(center == -1):
				from utilities      import estimate_3D_center_MPI, rotate_3D_shift
				cs[0], cs[1], cs[2], dummy, dummy = estimate_3D_center_MPI(data, total_nima, myid, number_of_proc, main_node)
				if myid == main_node:
					msg = " Average center x = %10.3f        Center y = %10.3f        Center z = %10.3f\n"%(cs[0], cs[1], cs[2])
					print_msg(msg)
				cs = mpi_bcast(cs, 3, MPI_FLOAT, main_node, MPI_COMM_WORLD)
				cs = [-float(cs[0]), -float(cs[1]), -float(cs[2])]
				rotate_3D_shift(data, cs)

			if CTF: vol, fscc = rec3D_MPI(data, snr, sym, fscmask, os.path.join(outdir, "resolution%04d"%(total_iter)), myid, main_node, npad = npad)
			else:    vol, fscc = rec3D_MPI_noCTF(data, sym, fscmask, os.path.join(outdir, "resolution%04d"%(total_iter)), myid, main_node, npad = npad)

			if myid == main_node:
				print_msg("3D reconstruction time = %d\n"%(time()-start_time))
				start_time = time()

			if fourvar:
			#  Compute Fourier variance
				for im in xrange(nima):
					original_data[im].set_attr( 'xform.projection', data[im].get_attr('xform.projection') )
				varf = varf3d_MPI(original_data, ssnr_text_file = os.path.join(outdir, "ssnr%04d"%(total_iter)), mask2D = None, reference_structure = vol, ou = last_ring, rw = 1.0, npad = 1, CTF = CTF, sign = 1, sym =sym, myid = myid)
				if myid == main_node:
					print_msg("Time to calculate 3D Fourier variance= %d\n"%(time()-start_time))
					start_time = time()
					varf = 1.0/varf
			else:  varf = None

			if myid == main_node:
				drop_image(vol, os.path.join(outdir, "vol%04d.hdf"%(total_iter)))
				ref_data[2] = vol
				ref_data[3] = fscc
				ref_data[4] = varf
				#  call user-supplied function to prepare reference image, i.e., center and filter it
				vol, cs = user_func(ref_data)
				drop_image(vol, os.path.join(outdir, "volf%04d.hdf"%(total_iter)))

			del varf
			#bcast_EMData_to_all(vol, myid, main_node)
			# write out headers, under MPI writing has to be done sequentially
			mpi_barrier(MPI_COMM_WORLD)
			par_str = ['xform.projection', 'ID']
			if myid == main_node:
	   			if(file_type(stack) == "bdb"):
	        			from utilities import recv_attr_dict_bdb
	        			recv_attr_dict_bdb(main_node, stack, data, par_str, image_start, image_end, number_of_proc)
	        		else:
	        			from utilities import recv_attr_dict
	        			recv_attr_dict(main_node, stack, data, par_str, image_start, image_end, number_of_proc)
				print_msg("Time to write header information= %d\n"%(time()-start_time))
				start_time = time()
	        	else:	       send_attr_dict(main_node, data, par_str, image_start, image_end)
	if myid == main_node: print_end_msg("ali3d_new_MPI")




'''def ali3d_n(stack, ref_vol, outdir, maskfile = None, ir = 1, ou = -1, rs = 1, 
            xr = "4 2 2 1", yr = "-1", ts = "1 1 0.5 0.25", delta = "10 6 4 4", an = "-1", deltapsi = "-1", startspi = "-1",
	    center = -1, maxit = 5, CTF = False, snr = 1.0,  ref_a = "S", sym = "c1",
	    user_func_name = "ref_ali3d", fourvar = True, debug = False, MPI = False):
	if MPI:
		ali3d_d_new_MPI(stack, ref_vol, outdir, maskfile, ir, ou, rs, xr, yr, ts,
	        	delta, an, deltapsi, startspi, center, maxit, CTF, snr, ref_a, sym, user_func_name,
			fourvar, debug)
		return

	from alignment      import proj_ali_incore, proj_ali_incore_local
	from utilities      import model_circle, drop_image, get_image, get_input_from_string
	from utilities      import get_params_proj, even_angles, assign_projangles
	from utilities      import estimate_3D_center, rotate_3D_shift
	from filter         import filt_params, fit_tanh, filt_tanl, filt_ctf
	from statistics     import fsc_mask
	import os
	import types
	from utilities      import print_begin_msg, print_end_msg, print_msg
	print_begin_msg("ali3d_n")

	from alignment      import Numrinit, ringwe, prepare_refrings
	from projection     import prep_vol

	import user_functions
	user_func = user_functions.factory[user_func_name]

	if os.path.exists(outdir):  ERROR('Output directory exists, please change the name and restart the program', " ", 1)
	os.mkdir(outdir)

	xrng        = get_input_from_string(xr)
	if  yr == "-1":  yrng = xrng
	else          :  yrng = get_input_from_string(yr)
	step        = get_input_from_string(ts)
	delta       = get_input_from_string(delta)
	lstp = min(len(xrng), len(yrng), len(step), len(delta))
	if an == "-1":
		an = [-1] * lstp
	else:
		an = get_input_from_string(an)
	first_ring  = int(ir)
	rstep       = int(rs)
	last_ring   = int(ou)
	max_iter    = int(maxit)
	center      = int(center)

	print_msg("Input stack                 : %s\n"%(stack))
	print_msg("Reference volume            : %s\n"%(ref_vol))	
	print_msg("Output directory            : %s\n"%(outdir))
	print_msg("Maskfile                    : %s\n"%(maskfile))
	print_msg("Inner radius                : %i\n"%(first_ring))

	vol     = EMData()
	vol.read_image(ref_vol)
	nx      = vol.get_xsize()
	cnx = nx//2+1
	if last_ring == -1:	last_ring = nx/2 - 2

	print_msg("Outer radius                : %i\n"%(last_ring))
	print_msg("Ring step                   : %i\n"%(rstep))
	print_msg("X search range              : %s\n"%(xrng))
	print_msg("Y search range              : %s\n"%(yrng))
	print_msg("Translational step          : %s\n"%(step))
	print_msg("Angular step                : %s\n"%(delta))
	print_msg("Angular search range        : %s\n"%(an))
	print_msg("Maximum iteration           : %i\n"%(max_iter))
	print_msg("Center type                 : %i\n"%(center))
	print_msg("CTF correction              : %s\n"%(CTF))
	print_msg("Signal-to-Noise Ratio       : %f\n"%(snr))
	print_msg("Reference projection method : %s\n"%(ref_a))
	print_msg("Symmetry group              : %s\n\n"%(sym))

	if maskfile :
		if type(maskfile) is types.StringType: mask3D = get_image(maskfile)
		else                                  : mask3D = maskfile
	else          :   mask3D = model_circle(last_ring, nx, nx, nx)
	mask2D = model_circle(last_ring, nx, nx) - model_circle(first_ring, nx, nx)
	numr   = Numrinit(first_ring, last_ring, rstep, "F")
	wr     = ringwe(numr,"F")

	if CTF:
		from reconstruction import recons3d_4nn_ctf
		from filter         import filt_ctf
	else: from reconstruction import recons3d_4nn

	if debug:  outf = file(os.path.join(outdir, "progress"), "w")
	else:      outf = None


	# horatio active_refactoring Jy51i1EwmLD4tWZ9_00000_1
	# active = EMUtil.get_all_attributes(stack, 'active')
	# list_of_particles = []
	# for im in xrange(len(active)):
	# 	if(active[im]):  list_of_particles.append(im)
	# del active
	
	nima = EMUtil.get_image_count(stack)
	list_of_particles = range(nima)

	data = EMData.read_images(stack, list_of_particles)
        for im in xrange(len(data)):
                data[im].set_attr('ID', list_of_particles[im])
		if CTF:
			ctf_params = data[im].get_attr("ctf")
			st = Util.infomask(data[im], mask2D, False)
			data[im] -= st[0]
			data[im] = filt_ctf(data[im], ctf_params)
			data[im].set_attr('ctf_applied', 1)

	nima = len(data)
	# initialize data for the reference preparation function
	ref_data = [ mask3D, max(center,0), None, None ]#  for center -1 switch of centering by user function
	cs = [0.0]*3
	# do the projection matching
	for N_step in xrange(lstp):
 		for Iter in xrange(max_iter):
			print_msg("\nITERATION #%3d\n"%(N_step*max_iter+Iter+1))
			if(N_step*max_iter+Iter+1 == 2): funct = "individual"
			elif(N_step*max_iter+Iter+1 == 1): funct = "conemove"
			else:  funct = "conepsi"
			print  funct
			from utilities    import get_params_proj, set_params_proj
			from alignment    import refprojs
			from random       import random
			if(funct == "individual"):
				volft,kb = prep_vol( vol )
				refrings = prepare_refrings( volft, kb, nx, delta[N_step], ref_a, sym, numr, MPI=False)
				del volft,kb

				for im in xrange( nima ):

					if an[N_step] == -1:	
						peak, pixel_error = proj_ali_incore(data[im],refrings,numr,xrng[N_step],yrng[N_step],step[N_step])
					else:
						peak, pixel_error = proj_ali_incore_local(data[im],refrings,numr,xrng[N_step],yrng[N_step],step[N_step],an[N_step])
				if center == -1:
					cs[0], cs[1], cs[2], dummy, dummy = estimate_3D_center(data)
					msg = "Average center x = %10.3f	Center y = %10.3f	 Center z = %10.3f\n"%(cs[0], cs[1], cs[2])
					print_msg(msg)
					rotate_3D_shift(data, [-cs[0], -cs[1], -cs[2]])

				if CTF:   vol1 = recons3d_4nn_ctf(data, range(0, nima, 2), snr, 1, sym)
				else:	   vol1 = recons3d_4nn(data, range(0, nima, 2), sym)
				if CTF:   vol2 = recons3d_4nn_ctf(data, range(1, nima, 2), snr, 1, sym)
				else:	   vol2 = recons3d_4nn(data, range(1, nima, 2), sym)

				fscc = fsc_mask(vol1, vol2, mask3D, 1.0, os.path.join(outdir, "resolution%04d"%(N_step*max_iter+Iter+1)))
				del vol1
				del vol2

				# calculate new and improved 3D
				if CTF:  vol = recons3d_4nn_ctf(data, range(nima), snr, 1, sym)
				else:	 vol = recons3d_4nn(data, range(nima), sym)
				# store the reference volume
				drop_image(vol, os.path.join(outdir, "vol%04d.hdf"%(N_step*max_iter+Iter+1)))
				ref_data[2] = vol
				ref_data[3] = fscc

				#  call user-supplied function to prepare reference image, i.e., center and filter it
				vol, dummy = user_func(ref_data)

				drop_image(vol, os.path.join(outdir, "volf%04d.hdf"%(N_step*max_iter+Iter+1)))
			elif(funct == "conepsi"):
				refa = even_angles(8.06)  #   Here it should be a parameter that specifies size of the cone
				angs = []
				for im in xrange( nima ):
					phi,theta, dummy1, dummy2, dummy3 = get_params_proj(data[im])
					angs.append([phi, theta])
				ass = assign_projangles(angs, refa)
				del angs, refa
				#  repeate alignment for each cone
				for kass in xrange(len(ass)):
					nicon = len(ass[kass])
					for ibi in xrange(2):  # number of big steps, there should be a criterion here....
						orgspi = [0.0]*nicon
						for im in xrange( nicon ):
								tmp0,tmp1,tmp2,tmp3,tmp4 = get_params_proj(data[ass[kass][im]])
								orgpsi[im] = 360.0*random()
								set_oarams_proj(data[ass[kass][im]],[tmp0,tmp1,orgpsi[im],tmp3,tmp4 ])
						keepgoing = True
						while(keepgoing):  # number of times psis are corrected
							if CTF:  vol = recons3d_4nn_ctf(data, ass[kass], snr, 1, sym, npad = 2)
							else:	 vol = recons3d_4nn(data, ass[kass], sym, npad = 2)
							volft,kb = prep_vol( vol )
							del vol
							keepgoing = False
							for im in xrange( nicon ):
								ref_a = [get_params_proj(data[ass[kass][im]])]
								refrings = refprojs( volft, kb, ref_a, last_ring, mask2D, cnx, cnx, numr, "F", wr )
								refrings[0].set_attr("phi",   ref_a[0][0])
								refrings[0].set_attr("theta", ref_a[0][1])
								refrings[0].set_attr("psi",   ref_a[0][2])
								peak, pixel_error = proj_ali_incore(data[ass[kass][im]], refrings, numr, xrng[N_step], yrng[N_step], step[N_step])
								tmp0,tmp1,tmp2,tmp3,tmp4 = get_params_proj(data[ass[kass][im]])
								if(orgpsi[im] != tmp2):
									keepgoing = True
									orgpsi[im] = tmp2
								#from sys import exit
								#exit()
						del volft,kb
						# calculate new and improved 3D
						#print  N_step*max_iter+Iter+1,Util.infomask(Util.muln_img(vol,vol),mask3D,True)
						# store the reference volume
						#drop_image(vol, os.path.join(outdir, "vol%04d.hdf"%(N_step*max_iter+Iter+1)))

						#from utilities import write_text_rows
						#write_text_rows(auta, os.path.join(outdir, "angles%04d.txt"%(N_step*max_iter+Iter+1)))
						auta = [0.0]*nicon
						for im in xrange( nicon ):
							tmp = get_params_proj(data[ass[kass][im]])
							auta[im] = tmp[2]
						cb = -1.0
						dang = 180.0/(2*3.14*last_ring)
						ding = int(180.0/dang +0.5)
						for iq in xrange(ding):
							for im in xrange( nicon ):
								tmp = get_params_proj(data[ass[kass][im]])
								set_params_proj(data[ass[kass][im]],[tmp[0],tmp[1],(auta[im]+dang*iq)%360.0,tmp[3],tmp[4]])
								if CTF:  vol = recons3d_4nn_ctf(data, ass[kass], snr, 1, sym, npad = 2)
								else:	 vol = recons3d_4nn(data, ass[kass], sym, npad = 2)
							pip = Util.infomask(Util.muln_img(vol,vol),mask3D,True)
							if(pip[0]>cb):
								cb = pip[0]
								mb = iq
								volt = vol.copy()
						#print  mb,cb
						for im in xrange( nicon ):
							tmp0,tmp1,tmp2,tmp3,tmp4 = get_params_proj(data[ass[kass][im]])
							tmp2 = (auta[im]+dang*mb)%360.0
						set_params_proj(data[ass[kass][im]],[tmp0,tmp1,tmp2,tmp3,tmp4])
						vol = volt.copy()
						del volt
			elif(funct == "conemove"):
				refa = even_angles(8.06)  #   Here it should be a parameter that specifies size of the cone
				angs = []
				for im in xrange( nima ):
					phi,theta,dummy1, dummy2, dummy3 = get_params_proj(data[im])
					angs.append([phi, theta])
				ass = assign_projangles(angs, refa)
				kass=len(ass)//2
				print  kass, len(ass[kass]),refa[kass]
				del angs, refa
				if CTF:  vol = recons3d_4nn_ctf(data, ass[kass], snr, 1, sym, npad = 2)
				else:	 vol = recons3d_4nn(data, ass[kass], sym, npad = 2)
				drop_image(vol, os.path.join(outdir, "Avol%04d.hdf"%(N_step*max_iter+Iter+1)))
				if CTF:  vol = recons3d_4nn_ctf(data, range(nima), snr, 1, sym)
				else:	 vol = recons3d_4nn(data, range(nima), sym)
				drop_image(vol, os.path.join(outdir, "Bvol%04d.hdf"%(N_step*max_iter+Iter+1)))
				
				from sys import exit
				exit()
			#  END OF BIG LOOP
		#  END OF CONES

			#  here we write header info
			from utilities import write_headers
			if CTF:
				for dat in data:  dat.del_attr('ctf_applied')
			write_headers(stack, data, list_of_particles)
			if CTF:
				for dat in data:  dat.set_attr('ctf_applied', 1)

	print_end_msg("ali3d_n")'''

def nearest_ref( vecs, vec ) :

	best_s = -1.0
	best_i = -1

	for i in xrange( len(vecs) ):
		s = abs(vecs[i][0]*vec[0] + vecs[i][1]*vec[1] + vecs[i][2]*vec[2])
		if s > best_s:
			best_s = s
			best_i = i

	return best_i

def assign_torefs(vecs, refangles):
	from utilities import getvec
	refnormal = [None]*len(refangles)
	for i in xrange(len(refangles)):
		refnormal[i] = getvec( refangles[i][0], refangles[i][1] )
	assignments = []
	for i in xrange(len(refnormal)):
		mm = nearest_ref( vecs, refnormal[i] )
		assignments.append(mm)
		del vecs[mm]
	return assignments
"""
def var_mpi_new(stack, outdir, scratch, fl, aa, radccc, writelp, writestack, frepa = "default", pca=False, pcamask=None, pcanvec=None, CTF = False):
	from string       import atoi, replace, split, atof
	from utilities    import get_im, circumference, model_circle, model_blank
	from utilities    import bcast_EMData_to_all, reduce_EMData_to_root
	from filter       import filt_tanl
	from mpi          import mpi_comm_rank, mpi_comm_size, mpi_barrier, mpi_bcast, mpi_reduce, mpi_send, mpi_recv
	from mpi          import MPI_COMM_WORLD, MPI_INT, MPI_FLOAT, MPI_SUM
	from utilities    import bcast_list_to_all, bcast_number_to_all, reduce_EMData_to_root, bcast_EMData_to_all 
	from utilities    import send_attr_dict
	from utilities    import get_params_proj, file_type, print_msg, print_begin_msg, print_end_msg
	from applications import MPI_start_end
	from projection   import prep_vol
	from statistics   import im_diff
	from morphology   import square_root
	from time import time	
	import os
	from sys import exit

        myid = mpi_comm_rank( MPI_COMM_WORLD )
        ncpu = mpi_comm_size( MPI_COMM_WORLD )
	main_node = 0
	nx = 0
        if myid==0:
		if os.path.exists(outdir):
			nx = 1
			ERROR('Output directory exists, please change the name and restart the program', " var_mpi", 0)
		else:
			os.system( "mkdir " + outdir )
	nx = mpi_bcast(nx, 1, MPI_INT, 0, MPI_COMM_WORLD)
	nx = int(nx[0])
	if(nx != 0):
		from sys import exit
		exit()

	mpi_barrier( MPI_COMM_WORLD )
	
	snr = 1.0  #  PUT IT AS A PARAMETER
	CTF = False
	
	debug = False  ############################################

	if myid == main_node:
		print_begin_msg("3D PCA")
		print_msg("Input stack                 : %s\n"%(stack))
		print_msg("Output directory            : %s\n"%(outdir))
		print_msg("Scratch directory           : %s\n"%(scratch))
		print_msg("Tangent low-pass filter   fl: %f\n"%(fl))
		print_msg("Tangent low-pass filter   aa: %f\n"%(aa))
		print_msg("Mask radius                 : %f\n"%(radccc))
		print_msg("Number of CPUs              : %i\n"%(ncpu))

       	if(len(scratch) > 4):
		if(scratch[:4] == "bdb:"):
			#scratch_file = scratch + "#diffproj%05d"%myid
			scratch_volume = scratch + "#gvbs%05d"%myid
		else:
			#scratch_file = os.path.join(scratch,"diffproj%05d.hdf"%myid)
			scratch_volume = os.path.join(scratch,"gvbs%05d.hdf"%myid)
	else:
		#scratch_file = os.path.join(scratch,"diffproj%05d.hdf"%myid)
		scratch_volume = os.path.join(scratch,"gvbs%05d.hdf"%myid)

# Step 0 - read input data
	if myid == main_node:
		start_time = time()
       		if(file_type(stack) == "bdb"):
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
	total_nima = bcast_number_to_all(nima, source_node = main_node)

	if myid != main_node:
		list_of_particles = [-1]*total_nima
	list_of_particles = bcast_list_to_all(list_of_particles, source_node = main_node)

	image_start, image_end = MPI_start_end(total_nima, ncpu, myid)
	# create a list of images for each node
	list_per_node = list_of_particles[image_start: image_end]
	nima = len(list_per_node)

	if debug:
		finfo.write("image_start, image_end: %d %d\n" %(image_start, image_end))
		finfo.flush()
	data = EMData.read_images(stack, list_of_particles)
	'''
	for i in xrange(ncpu):
		if myid == i:
			data = EMData.read_images(stack, list_of_particles)
		mpi_barrier(MPI_COMM_WORLD)
	'''
	if myid== main_node:
		print_msg( "Time to read data: %d\n" % (time()-start_time) )
		start_time = time()


	t = EMUtil.get_all_attributes(stack,"xform.projection")
	nimages = len(t)

	if fl>0.0:
		for k in list_of_particles:
			data[k] = filt_tanl(data[k], fl ,aa)
	from utilities import getvec, even_angles
	#  Here only those that are on the list of particles should remain, do it later
	vecs = [None]*nimages
	for i in xrange(nimages):
		d = t[i].get_params("spider")
		d = getvec(d["phi"],d["theta"])
		vecs[i] = [d[0],d[1],d[2]]

	bunch = 18#len(even_angles(30.0))
	nset = nimages//bunch
	from reconstruction import recons3d_wbp
	for i in xrange(nset):
		# set of vectors in a given bunch
		#refa = even_angles(30.0,11.0/(nset-1)*i)
		#print len(vecs)
		bset = []
		for k in xrange(i,nimages,nset):  bset.append(k)
		#bset = assign_torefs(vecs,refa)
		recons3d_wbp(data, bset).write_image("bset.hdf", i)
		print i
	exit()
# Step 1 - compute sum_CTF^2
	nx = data[image_start].get_xsize()
	ny = data[image_start].get_ysize()
	nz = nx
	if(radccc < 1):  radcir = min(nx,ny)//2-2
	else:            radcir = radccc
	mask2D = model_circle( radcir, nx, ny)

	if CTF:
		from filter import filt_ctf
		from fundamentals import fftip, fft
		from morphology  import  ctf_img
		#  Compute sum_CTF^2+1.0/snr  in 2D
		sumCTF2 = EMData(nx,ny,1, False)
		previous_defocus = -1.0
		for k in xrange(image_start, image_end):
			ctf_params = data[k].get_attr( "ctf" )	
			if ctf_params.defocus != previous_defocus:
				previous_defocus = ctf_params.defocus
				vb = ctf_img(nx, ctf_params,1, ny)
				Util.mul_img(vb, vb)
			Util.add_img(sumCTF2, vb)
		reduce_EMData_to_root(sumCTF2, myid)
		if myid == 0:  sumCTF2 += (1.0/snr)
		bcast_EMData_to_all( sumCTF2, myid )

	#filter data
	if fl>0.0:
		for k in list_of_particles:
			data[k] = filt_tanl(data[k], fl, aa)

	if myid== main_node and CTF:
		print_msg( "Time to compute sum CTF^2 %d\n" % (time()-start_time) )
		start_time = time()

# Step 3 - compute overall 3D
	from reconstruction import recons3d_wbp	
	'''
	avg1 = recons3d_wbp(data, list_of_particles)

	if( myid == 0):
		#avg1.write_image(os.path.join(outdir, "avgb.hdf"))
		#sumCTF2.write_image(os.path.join(outdir, "s2.hdf"))
		Util.mult_scalar(avg1, 1.0/float(total_nima)).write_image(os.path.join(outdir,"avg1.hdf"),0)  # leave the volume without division!
	bcast_EMData_to_all( avg1, myid )
	if myid== main_node:
		print_msg( "Time to compute 3D: %d\n" % (time()-start_time) )
		start_time = time()
	'''
# Step 4 - resample
	'''
	#  Will be subtracting, so change the sign
	for k in xrange(total_nima):
		Util.mul_scalar(data[k], -1.0)

	lin = range(total_nima)
	from random import shuffle
	from reconstruction import one_swbp
	nvpercpu = 20000/ncpu
	remove= 100
	for b in xrange(nvpercpu):
		bv = avg1.copy()
		shuffle(lin)
		for k in xrange(remove):
			one_swbp(bv, data[lin[k]], list_of_particles[lin[k]], dm)
		Util.mul_scalar(bv, 1.0/(total_nima-remove))
		bv.write_image(scratch_volume,b)
	if myid== main_node:
		print_msg( "Time to resample projections: %d\n" % (time()-start_time) )
		start_time = time()
	exit()
	'''
# Step 4 - bootstrap

	from random import randint
	from reconstruction import one_swbp
	for b in xrange(10000//ncpu):
		lin = [0]*total_nima
		for k in xrange(total_nima):  lin[k] = list_of_particles[randint(0,total_nima-1)]
		avg1 = recons3d_wbp(data, lin)
		Util.mul_scalar(avg1, 1.0/(total_nima))
		avg1.write_image(scratch_volume,b)
		if myid== main_node:
			print_msg( "Time to bootstrap : %d\n" % (time()-start_time) )
			start_time = time()
	exit()

       	if(len(scratch) > 4):
		if(scratch[:4] == "bdb:"):
			scratch_dir = scratch[4:]
		else:
			scratch_dir = scratch
	else:
		scratch_dir = scratch
	pcaer = pcanalyzer(pcamask, outdir, pcanvec, True, scratch_dir)



	var1 = square_root(var1)

	if( myid == 0 ):
		vb = backproject_swbp(data[0], list_of_particles[k], dm)
		
		#if(fl > 0.0):
		#	vb = filt_tanl(vb, fl, aa)
		vb = Util.divn_img(vb, var1)
		refstat = Util.infomask(vb, pcamask, True)
	else:             refstat = [0.0,0.0,0.0,0.0]
	refstat = mpi_bcast(refstat, 4, MPI_FLOAT, 0, MPI_COMM_WORLD)
	refstat = map(float, refstat)

	avg2 = model_blank(nx,ny,nz)

	for k in xrange(nima):
		vb = Util.divn_img( backproject_swbp(data[k], list_of_particles[k], dm) , var1)
		#if(fl > 0.0):
		#	vb = filt_tanl(vb, fl, aa)
		pc = Util.infomask(vb, pcamask, True)
		vb -= pc[0]
		vb *= (refstat[1]/pc[1])
		Util.add_img(avg2, vb)
		pcaer.insert(vb)
	del data,vb
	reduce_EMData_to_root(avg2, myid)
	if( myid == 0):
		Util.mul_scalar(avg2, 1.0/float(total_nima))
	else:   avg2 =  model_blank(nx,ny,nz)
	bcast_EMData_to_all( avg2, myid )
	if(myid == 0):  avg2.write_image(os.path.join(outdir,"avg2.hdf"))

	assert not(avg2 is None)
	if myid== main_node:
		print_msg( "Time to divide backprojected difference projections by the standard deviation: %d\n" % (time()-start_time) )
		start_time = time()

# Step 6 - do the PCA on standardized projections
	if myid == main_node:  start_time = time()
	
	pcaer.setavg( avg2 )

	eigs = pcaer.analyze()
	from mpi import mpi_barrier
	mpi_barrier(MPI_COMM_WORLD)
	if myid==0:
		eigfile = os.path.join(outdir, "eigvol.hdf")
		for i in xrange( len(eigs) ):
			eigs[i].write_image( eigfile, i )
	if myid== main_node:
		print_msg( "Time to do PCA: %d\n" % (time()-start_time) )
		start_time = time()
	# Clean up masked files
	os.system("rm -f "+ os.path.join(scratch_dir , "maskedimg%04d.bin" % myid ) )
	if myid == main_node: print_end_msg("3D PCA")
"""
"""
# proportional resampling
def var_mpi_new(stack, outdir, scratch, fl, aa, radccc, writelp, writestack, frepa = "default", pca=False, pcamask=None, pcanvec=None, CTF = False):
	from string       import atoi, replace, split, atof
	from utilities    import get_im, circumference, model_circle, model_blank
	from utilities    import bcast_EMData_to_all, reduce_EMData_to_root
	from filter       import filt_tanl
	from mpi          import mpi_comm_rank, mpi_comm_size, mpi_barrier, mpi_bcast, mpi_reduce, mpi_send, mpi_recv
	from mpi          import MPI_COMM_WORLD, MPI_INT, MPI_FLOAT, MPI_SUM
	from utilities    import bcast_list_to_all, bcast_number_to_all, reduce_EMData_to_root, bcast_EMData_to_all 
	from utilities    import send_attr_dict
	from utilities    import get_params_proj, file_type, print_msg, print_begin_msg, print_end_msg
	from applications import MPI_start_end
	from projection   import prep_vol
	from statistics   import im_diff
	from morphology   import square_root
	from time import time	
	import os
	from sys import exit

        myid = mpi_comm_rank( MPI_COMM_WORLD )
        ncpu = mpi_comm_size( MPI_COMM_WORLD )
	main_node = 0
	nx = 0
        if myid==0:
		if os.path.exists(outdir):
			nx = 1
			ERROR('Output directory exists, please change the name and restart the program', " var_mpi", 0)
		else:
			os.system( "mkdir " + outdir )
	nx = mpi_bcast(nx, 1, MPI_INT, 0, MPI_COMM_WORLD)
	nx = int(nx[0])
	if(nx != 0):
		from sys import exit
		exit()

	mpi_barrier( MPI_COMM_WORLD )
	
	snr = 1.0  #  PUT IT AS A PARAMETER
	CTF = True
	
	debug = False  ############################################

	if myid == main_node:
		print_begin_msg("3D PCA")
		print_msg("Input stack                 : %s\n"%(stack))
		print_msg("Output directory            : %s\n"%(outdir))
		print_msg("Scratch directory           : %s\n"%(scratch))
		print_msg("Tangent low-pass filter   fl: %f\n"%(fl))
		print_msg("Tangent low-pass filter   aa: %f\n"%(aa))
		print_msg("Mask radius                 : %f\n"%(radccc))
		print_msg("Number of CPUs              : %i\n"%(ncpu))


       	if(len(scratch) > 4):
		if(scratch[:4] == "bdb:"):
			#scratch_file = scratch + "#diffproj%05d"%myid
			scratch_volume = scratch + "#gvbs%05d"%myid
		else:
			#scratch_file = os.path.join(scratch,"diffproj%05d.hdf"%myid)
			scratch_volume = os.path.join(scratch,"gvbs%05d.hdf"%myid)
	else:
		#scratch_file = os.path.join(scratch,"diffproj%05d.hdf"%myid)
		scratch_volume = os.path.join(scratch,"gvbs%05d.hdf"%myid)

# Step 0 - read input data
	if myid == main_node:
		nima = 50000
	else:
		nima = 0
	total_nima = bcast_number_to_all(nima, source_node = main_node)

	list_of_particles = range(total_nima)
	'''

	image_start, image_end = MPI_start_end(total_nima, ncpu, myid)
	# create a list of images for each node
	list_per_node = list_of_particles[image_start: image_end]
	nima = len(list_per_node)

	if debug:
		finfo.write("image_start, image_end: %d %d\n" %(image_start, image_end))
		finfo.flush()
	data = [None]*total_nima
	'''
	'''
	for i in xrange(ncpu):
		if myid == i:
			data = EMData.read_images(stack, list_of_particles)
		mpi_barrier(MPI_COMM_WORLD)
	'''
	if myid == main_node:  start_time = time()
	data = EMData.read_images(stack, list_of_particles)
	if myid== main_node:
		print_msg( "Time to read data: %d\n" % (time()-start_time) )
		start_time = time()
	nx = data[0].get_xsize()
	ny = data[0].get_ysize()
	nz = nx
	if(radccc < 1):  radcir = min(nx,ny)//2-2
	else:            radcir = radccc
	mask2D = model_circle( radcir, nx, ny)
	if CTF:
		from filter import filt_ctf
		for im in xrange(total_nima):
			st = Util.infomask(data[im], mask2D, False)
			data[im] -= st[0]
			st = data[im].get_attr_default("ctf_applied", 0)
			if(st == 0):
				ctf_params = data[im].get_attr("ctf")
				data[im] = filt_ctf(data[im], ctf_params)
				data[im].set_attr('ctf_applied', 1)
	if myid== main_node:
		print_msg( "Time to apply CTF: %d\n" % (time()-start_time) )
		start_time = time()

# Step 4 - jackknife-d proportionally drawn.
	if  CTF:
		from reconstruction import recons3d_4nn_ctf
	else:
		from reconstruction import recons3d_4nn
	from mpi import MPI_TAG_UB
	from utilities import even_angles
	refa = even_angles(10.4)#14.5)#6.51)  # yields 500;  14.5) # yields 100
	from utilities import assign_projangles
	from random import randint, seed, jumpahead
	ang = [None]*total_nima
	for k in xrange(total_nima):
	        phi, theta, dummy, dummy, dummy = get_params_proj(data[k])
	        ang[k] = [phi, theta]
	lass = assign_projangles(ang, refa)
	ivol = len(lass[0])
	for i in xrange(1,len(lass)):
	        mivol = min(len(lass[i]),ivol)
	d = 1.0/mivol
	if myid == 0:  print  "lass",len(lass),mivol,d

	seed()
	jumpahead(17*myid+123)

	ivol = 500
	##print "ivol",ivol
	nb = 0
	nlass= len(lass)
	for b in xrange(0,ivol):
		mass= [None]*nlass
		for l in xrange(nlass):
			mass[l] = lass[l][:]
		jack = range(total_nima)
		for l in xrange(nlass):
			to_go = int(d*len(lass[l])+0.5)
			for i in xrange(to_go):
				ti = randint(0,len(mass[l])-1)
				td =mass[l][ti]
				del mass[l][ti]
				del jack[jack.index(td)]

		if  CTF:
			recons3d_4nn_ctf(data, jack, snr = snr, npad=2).write_image(scratch_volume, nb)
		else:
			recons3d_4nn(data, jack, npad=2).write_image(scratch_volume, nb)
		nb += 1
		mpi_barrier(MPI_COMM_WORLD)
		if myid == 0: print "  NEW ITERATION  ",b

	if myid== main_node:
		print_msg( "Time to jacknife projections: %d\n" % (time()-start_time) )
		start_time = time()
"""
"""
#This worked on 07/25/10
def var_mpi_new(stack, outdir, scratch, fl, aa, radccc, writelp, writestack, frepa = "default", pca=False, pcamask=None, pcanvec=None, CTF = False):
	from string       import atoi, replace, split, atof
	from utilities    import get_im, circumference, model_circle, model_blank
	from utilities    import bcast_EMData_to_all, reduce_EMData_to_root
	from filter       import filt_tanl
	from mpi          import mpi_comm_rank, mpi_comm_size, mpi_barrier, mpi_bcast, mpi_reduce, mpi_send, mpi_recv
	from mpi          import MPI_COMM_WORLD, MPI_INT, MPI_FLOAT, MPI_SUM
	from utilities    import bcast_list_to_all, bcast_number_to_all, reduce_EMData_to_root, bcast_EMData_to_all 
	from utilities    import send_attr_dict
	from utilities    import get_params_proj, file_type, print_msg, print_begin_msg, print_end_msg
	from applications import MPI_start_end
	from projection   import prep_vol
	from statistics   import im_diff
	from morphology   import square_root
	from time import time	
	import os
	from sys import exit

        myid = mpi_comm_rank( MPI_COMM_WORLD )
        ncpu = mpi_comm_size( MPI_COMM_WORLD )
	main_node = 0
	nx = 0
        if myid==0:
		if os.path.exists(outdir):
			nx = 1
			ERROR('Output directory exists, please change the name and restart the program', " var_mpi", 0)
		else:
			os.system( "mkdir " + outdir )
	nx = mpi_bcast(nx, 1, MPI_INT, 0, MPI_COMM_WORLD)
	nx = int(nx[0])
	if(nx != 0):
		from sys import exit
		exit()

	mpi_barrier( MPI_COMM_WORLD )
	
	snr = 1.0  #  PUT IT AS A PARAMETER
	CTF = True
	
	debug = False  ############################################

	if myid == main_node:
		print_begin_msg("3D PCA")
		print_msg("Input stack                 : %s\n"%(stack))
		print_msg("Output directory            : %s\n"%(outdir))
		print_msg("Scratch directory           : %s\n"%(scratch))
		print_msg("Tangent low-pass filter   fl: %f\n"%(fl))
		print_msg("Tangent low-pass filter   aa: %f\n"%(aa))
		print_msg("Mask radius                 : %f\n"%(radccc))
		print_msg("Number of CPUs              : %i\n"%(ncpu))


       	if(len(scratch) > 4):
		if(scratch[:4] == "bdb:"):
			#scratch_file = scratch + "#diffproj%05d"%myid
			scratch_volume = scratch + "#gvbs%05d"%myid
		else:
			#scratch_file = os.path.join(scratch,"diffproj%05d.hdf"%myid)
			scratch_volume = os.path.join(scratch,"gvbs%05d.hdf"%myid)
	else:
		#scratch_file = os.path.join(scratch,"diffproj%05d.hdf"%myid)
		scratch_volume = os.path.join(scratch,"gvbs%05d.hdf"%myid)

# Step 0 - read input data
	if myid == main_node:
		start_time = time()
       		if(file_type(stack) == "bdb"):
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
	total_nima = bcast_number_to_all(nima, source_node = main_node)

	if myid != main_node:
		list_of_particles = [-1]*total_nima
	list_of_particles = bcast_list_to_all(list_of_particles, source_node = main_node)
	'''

	image_start, image_end = MPI_start_end(total_nima, ncpu, myid)
	# create a list of images for each node
	list_per_node = list_of_particles[image_start: image_end]
	nima = len(list_per_node)

	if debug:
		finfo.write("image_start, image_end: %d %d\n" %(image_start, image_end))
		finfo.flush()
	data = [None]*total_nima
	'''
	'''
	for i in xrange(ncpu):
		if myid == i:
			data = EMData.read_images(stack, list_of_particles)
		mpi_barrier(MPI_COMM_WORLD)
	'''
	data = EMData.read_images(stack, list_of_particles)
	if myid== main_node:
		print_msg( "Time to read data: %d\n" % (time()-start_time) )
		start_time = time()
	nx = data[0].get_xsize()
	ny = data[0].get_ysize()
	nz = nx
	if(radccc < 1):  radcir = min(nx,ny)//2-2
	else:            radcir = radccc
	mask2D = model_circle( radcir, nx, ny)
	if CTF:
		from filter import filt_ctf
		for im in xrange(total_nima):
			st = Util.infomask(data[im], mask2D, False)
			data[im] -= st[0]
			st = data[im].get_attr_default("ctf_applied", 0)
			if(st == 0):
				ctf_params = data[im].get_attr("ctf")
				data[im] = filt_ctf(data[im], ctf_params)
				data[im].set_attr('ctf_applied', 1)
	if myid== main_node:
		print_msg( "Time to apply CTF: %d\n" % (time()-start_time) )
		start_time = time()
	'''
# Step 1 - compute sum_CTF^2

	if CTF:
		from filter import filt_ctf
		from fundamentals import fftip, fft
		from morphology  import  ctf_img
		#  Compute sum_CTF^2+1.0/snr  in 2D
		sumCTF2 = EMData(nx,ny,1, False)
		previous_defocus = -1.0
		for k in xrange(image_start, image_end):
			ctf_params = data[k].get_attr( "ctf" )	
			if ctf_params.defocus != previous_defocus:
				previous_defocus = ctf_params.defocus
				vb = ctf_img(nx, ctf_params,1, ny)
				Util.mul_img(vb, vb)
			Util.add_img(sumCTF2, vb)
		reduce_EMData_to_root(sumCTF2, myid)
		if myid == 0:  sumCTF2 += (1.0/snr)
		bcast_EMData_to_all( sumCTF2, myid )

	if myid== main_node and CTF:
		print_msg( "Time to compute sum CTF^2 %d\n" % (time()-start_time) )
		start_time = time()

# Step 2 - weight projections
	nsym = 1
	nimages = total_nima
	ntripletsWnsym = nsym*nimages
	dm=[0.0]*(9*ntripletsWnsym)
	ss=[0.0]*(6*ntripletsWnsym)
	if myid == main_node:
		t = EMUtil.get_all_attributes(stack,"xform.projection")
		#  Here only those that are on the list of particles should remain, do it later
		count = 0
		for i in xrange(nimages):
			d = t[i].get_params("spider")
			DMnSS = Util.CANG(d["phi"],d["theta"],d["psi"])
			dm[(count*9) :(count+1)*9] = DMnSS["DM"]
			ss[(count*6) :(count+1)*6] = DMnSS["SS"]
			count += 1
		del t, d
	dm = bcast_list_to_all(dm, source_node = main_node)
	ss = bcast_list_to_all(ss, source_node = main_node)
	from reconstruction import weight_swbp	
	for k in xrange(image_start, image_end):
		t = data[k].get_attr("xform.projection")
		if(fl > 0.0):
			df = filt_tanl( data[k], fl, aa )
		else:
			df = data[k]
		if  CTF:
			df = Util.divn_img(fft(df),sumCTF2)
			df = fft(filt_ctf(df, data[k].get_attr( "ctf" )))
		df.set_attr("xform.projection",t)
		df.set_attr("active",1)
		data[k] = weight_swbp(df, list_of_particles[k], dm, ss, "general")
		#  Data contains weighted 2D diffprojdata
		#  If CTF, it is also multiplied by the CTF and divided by the sum_CTF^2
		data[k].set_attr("xform.projection",t)
		data[k].set_attr("active",1)
	del ss, df
	if myid== main_node:
		print_msg( "Time to weight projections: %d\n" % (time()-start_time) )
		start_time = time()
	'''
	'''
# Step 3 - compute overall 3D
	from reconstruction import backproject_swbp	
	avg1 = model_blank(nx,ny,nz)
	for k in xrange(image_start, image_end):
		Util.add_img(avg1, backproject_swbp(data[k], list_of_particles[k], dm) )
	reduce_EMData_to_root(avg1, myid)

	if( myid == 0):
		#avg1.write_image(os.path.join(outdir, "avgb.hdf"))
		#sumCTF2.write_image(os.path.join(outdir, "s2.hdf"))
		Util.mult_scalar(avg1, 1.0/float(total_nima)).write_image(os.path.join(outdir,"avg1.hdf"),0)  # leave the volume without division!
	bcast_EMData_to_all( avg1, myid )
	if myid == main_node:  start_time = time()

	if  CTF:  
		from reconstruction import recons3d_4nn_ctf
		volr = recons3d_4nn_ctf(data, list_of_particles, snr = 1.0)
	else:     
		from reconstruction import recons3d_4nn
		volr = recons3d_4nn(data, list_of_particles)

	if myid== main_node:
		volr.write_image(os.path.join(outdir, "avg.hdf"))
		del volr
		print_msg( "Time to compute 3D: %d\n" % (time()-start_time) )
		start_time = time()
# Step 4 - redistribute particles
	# Now we have to redistribute the particles
	mpi_barrier(MPI_COMM_WORLD)
	for i in xrange(ncpu):
	        ib, ie = MPI_start_end(total_nima, ncpu, i)
	        for k in xrange(ib, ie):
			if(i != myid): data[k] = model_blank(nx,ny)
			else:	t = data[k].get_attr("xform.projection")
	        	bcast_EMData_to_all( data[k], myid, i)
			if(i != myid):
				data[k].set_attr("xform.projection",t)
				data[k].set_attr("active",1)
		mpi_barrier(MPI_COMM_WORLD)
	del t

	if myid== main_node:
		print_msg( "Time to broadcast projections: %d\n" % (time()-start_time) )
		start_time = time()
	'''
	'''
	avg1 = model_blank(nx,ny,nz)
	for k in xrange(image_start, image_end):
		Util.add_img(avg1, backproject_swbp(data[k], list_of_particles[k], dm) )
	reduce_EMData_to_root(avg1, myid)

	if( myid == 0):
		#avg1.write_image(os.path.join(outdir, "avgb.hdf"))
		#sumCTF2.write_image(os.path.join(outdir, "s2.hdf"))
		Util.mult_scalar(avg1, 1.0/float(total_nima)).write_image(os.path.join(outdir,"avg2.hdf"),0)  # leave the volume without division!
	bcast_EMData_to_all( avg1, myid )
	if myid== main_node:
		print_msg( "Time to compute 3D: %d\n" % (time()-start_time) )
		start_time = time()
	#exit()
	'''

# Step 4 - resample
	'''
	#  Will be subtracting, so change the sign
	for k in xrange(total_nima):
		Util.mul_scalar(data[k], -1.0)

	lin = range(total_nima)
	from random import shuffle
	from reconstruction import one_swbp
	nvpercpu = 20000/ncpu
	remove= 100
	for b in xrange(nvpercpu):
		bv = avg1.copy()
		shuffle(lin)
		for k in xrange(remove):
			one_swbp(bv, data[lin[k]], list_of_particles[lin[k]], dm)
		Util.mul_scalar(bv, 1.0/(total_nima-remove))
		bv.write_image(scratch_volume,b)
	if myid== main_node:
		print_msg( "Time to resample projections: %d\n" % (time()-start_time) )
		start_time = time()
	exit()
	'''
	'''
# Step 4 - bootstrap

	from random import randint
	from reconstruction import one_swbp
	for b in xrange(10000//ncpu):
		lin = [0]*total_nima
		for k in xrange(total_nima):  lin[k] = randint(0,total_nima-1)
		lin.sort()
		sel = []
		q = -1
		for k in xrange(total_nima):
			if(lin[k] != q):
				sel.append([lin[k],1])
				q = lin[k]
			else:
				sel[-1][1] +=1

		bv = model_blank(nx,nx,nx)
		for k in xrange(len(sel)):
			one_swbp(bv, Util.mult_scalar(data[sel[k][0]], float(sel[k][1])), list_of_particles[sel[k][0]], dm)
		Util.mul_scalar(bv, 1.0/(total_nima))
		bv.write_image(scratch_volume,b)
	if myid== main_node:
		print_msg( "Time to bootstrap projections: %d\n" % (time()-start_time) )
		start_time = time()
	exit()

	'''
# Step 4 - jackknife-d evenly drawn.
	if  CTF:
		from reconstruction import recons3d_4nn_ctf
	else:
		from reconstruction import recons3d_4nn
	'''
	d = 0.001 #percent to be taken out
	#from utilities import write_text_file
	takeout = int(total_nima*d)
	ivol = total_nima//min(takeout, total_nima-takeout)
	nb = 0
	for b in xrange(ivol):
		if(b%ncpu == myid):
			if d<0.5:
	        		jack = [0]*total_nima
	        		for j in xrange(total_nima):  jack[j] = list_of_particles[j]
	        		for j in xrange(b,total_nima,ivol):
	        			jack[j] = -1
	        		try:
	        			while True:
	        				j = jack.index(-1)
	        				del jack[j]
	        		except:
	        			pass
			else:
				jack = []
				for j in xrange(b,total_nima,ivol):
					jack.append(list_of_particles[j])
			#write_text_file(jack,"sjk%05d.txt"%b)
			print "STAMP",myid,total_nima,takeout,ivol,b,len(jack)
	        	recons3d_4nn(data, jack, npad=2).write_image(scratch_volume,nb)
	        	nb += 1
	'''
	from mpi import MPI_TAG_UB
	if myid == main_node:
		from utilities import even_angles
		refa = even_angles(6.51)  # yields 500;  14.5) # yields 100
		from utilities import assign_projangles
		from random import randint
		ang = [None]*total_nima
		for k in xrange(total_nima):
		        phi, theta, dummy, dummy, dummy = get_params_proj(data[k])
		        ang[k] = [phi, theta]
		lass = assign_projangles(ang, refa)
		ivol = len(lass[0])
		for i in xrange(1,len(lass)):
		        ivol = min(len(lass[i]),ivol)
		lst = total_nima - len(lass)
		print  "lass",len(lass)
	else:
		ivol = 0
		lst  = 0
	ivol = bcast_number_to_all(ivol, source_node = main_node)
	lst = bcast_number_to_all(lst, source_node = main_node)
	#ivol = 1000
	print "ivol",ivol
	nb = 0
	for b in xrange(0,ivol,ncpu):
		if myid == 0:
			print  "ho  ",ivol,b,ncpu
		active_nodes = min(ivol,b+ncpu)%ncpu
		if myid == 0:  print "active nodes  A",active_nodes
		if active_nodes == 0:  active_nodes = ncpu
		if myid == 0:  print  "active nodes B",active_nodes,total_nima
	        if(myid != main_node):
	        	for i in xrange(active_nodes):
				if i == myid:
	        			jack = mpi_recv(lst, MPI_INT, main_node, MPI_TAG_UB, MPI_COMM_WORLD)
					print  " receiving  jack  ",i,len(jack)
		        		if  CTF:
						recons3d_4nn_ctf(data, jack, snr = snr, npad=2).write_image(scratch_volume, nb)
					else:
						recons3d_4nn(data, jack, npad=2).write_image(scratch_volume, nb)
		        		nb += 1
	        else:
			print "active_nodes  ",active_nodes
	        	for i in xrange(active_nodes):
	        		list_to_send = [0]*total_nima
	        		for l in xrange(total_nima): list_to_send[l] = list_of_particles[l]
	        		for l in xrange(len(lass)):
	        			ti = randint(0,len(lass[l])-1)
	        			td =lass[l][ti]
	        			del lass[l][ti]           #  HERE  -- commment out to get resampling within a cone
	        			del list_to_send[list_to_send.index(td)]
	        		if(i != main_node):
					print  " sending jack to",i
	        			mpi_send(list_to_send, lst, MPI_INT, i, MPI_TAG_UB, MPI_COMM_WORLD)
	        		else:
					print  "doing jack  ",i
	        			jack = [0]*lst
	        			for l in xrange(lst):  jack[l] = list_to_send[l]

			print  myid,len(jack)
			if  CTF:
				recons3d_4nn_ctf(data, jack, snr = snr, npad=2).write_image(scratch_volume, nb)
			else:
		        	recons3d_4nn(data, jack, npad=2).write_image(scratch_volume, nb)
	        	nb += 1
		mpi_barrier(MPI_COMM_WORLD)
		from time import sleep
		sleep(1)
		if myid == 0: print "  NEW ITERATION  ",b
			
	if myid== main_node:
		print_msg( "Time to jacknife projections: %d\n" % (time()-start_time) )
		start_time = time()
	exit()

       	if(len(scratch) > 4):
		if(scratch[:4] == "bdb:"):
			scratch_dir = scratch[4:]
		else:
			scratch_dir = scratch
	else:
		scratch_dir = scratch
	pcaer = pcanalyzer(pcamask, outdir, pcanvec, True, scratch_dir)



	var1 = square_root(var1)

	if( myid == 0 ):
		vb = backproject_swbp(data[0], list_of_particles[k], dm)
		
		#if(fl > 0.0):
		#	vb = filt_tanl(vb, fl, aa)
		vb = Util.divn_img(vb, var1)
		refstat = Util.infomask(vb, pcamask, True)
	else:             refstat = [0.0,0.0,0.0,0.0]
	refstat = mpi_bcast(refstat, 4, MPI_FLOAT, 0, MPI_COMM_WORLD)
	refstat = map(float, refstat)

	avg2 = model_blank(nx,ny,nz)

	for k in xrange(nima):
		vb = Util.divn_img( backproject_swbp(data[k], list_of_particles[k], dm) , var1)
		#if(fl > 0.0):
		#	vb = filt_tanl(vb, fl, aa)
		pc = Util.infomask(vb, pcamask, True)
		vb -= pc[0]
		vb *= (refstat[1]/pc[1])
		Util.add_img(avg2, vb)
		pcaer.insert(vb)
	del data,vb
	reduce_EMData_to_root(avg2, myid)
	if( myid == 0):
		Util.mul_scalar(avg2, 1.0/float(total_nima))
	else:   avg2 =  model_blank(nx,ny,nz)
	bcast_EMData_to_all( avg2, myid )
	if(myid == 0):  avg2.write_image(os.path.join(outdir,"avg2.hdf"))

	assert not(avg2 is None)
	if myid== main_node:
		print_msg( "Time to divide backprojected difference projections by the standard deviation: %d\n" % (time()-start_time) )
		start_time = time()

# Step 6 - do the PCA on standardized projections
	if myid == main_node:  start_time = time()
	
	pcaer.setavg( avg2 )

	eigs = pcaer.analyze()
	from mpi import mpi_barrier
	mpi_barrier(MPI_COMM_WORLD)
	if myid==0:
		eigfile = os.path.join(outdir, "eigvol.hdf")
		for i in xrange( len(eigs) ):
			eigs[i].write_image( eigfile, i )
	if myid== main_node:
		print_msg( "Time to do PCA: %d\n" % (time()-start_time) )
		start_time = time()
	# Clean up masked files
	os.system("rm -f "+ os.path.join(scratch_dir , "maskedimg%04d.bin" % myid ) )
	if myid == main_node: print_end_msg("3D PCA")

"""

# This one contains CTF and changed filtration.  It is REALLY OK. 07/16/2010
def var_mpi_new(stack, outdir, scratch, fl, aa, radccc, writelp, writestack, frepa = "default",\
       pca=False, pcamask=None, pcanvec=None, CTF = False):
	from string       import atoi, replace, split, atof
	from utilities    import get_im, circumference, model_circle, model_blank
	from utilities    import bcast_EMData_to_all, reduce_EMData_to_root
	from filter       import filt_tanl
	from mpi          import mpi_comm_rank, mpi_comm_size, mpi_barrier, mpi_bcast, mpi_reduce
	from mpi          import MPI_COMM_WORLD, MPI_INT, MPI_FLOAT, MPI_SUM
	from utilities    import bcast_list_to_all, bcast_number_to_all, reduce_EMData_to_root, bcast_EMData_to_all 
	from utilities    import send_attr_dict
	from utilities    import get_params_proj, file_type, print_msg, print_begin_msg, print_end_msg
	from applications import MPI_start_end
	from projection   import prep_vol
	from statistics   import im_diff
	from morphology   import square_root
        from reconstruction import recons3d_swbp, backproject_swbp
	from time import time	
	import os
	from sys import exit

        myid = mpi_comm_rank( MPI_COMM_WORLD )
        ncpu = mpi_comm_size( MPI_COMM_WORLD )
	main_node = 0
	nx = 0
        if myid==0:
		if os.path.exists(outdir):
			nx = 1
			ERROR('Output directory exists, please change the name and restart the program', " var_mpi", 0)
		else:
			os.system( "mkdir " + outdir )
	nx = mpi_bcast(nx, 1, MPI_INT, 0, MPI_COMM_WORLD)
	nx = int(nx[0])
	if(nx != 0):
		from sys import exit
		exit()

	mpi_barrier( MPI_COMM_WORLD )
	
	snr = 1.0  #  PUT IT AS A PARAMETER
	CTF = True
	
	debug = False  ############################################

	if myid == main_node:
		print_begin_msg("3D PCA")
		print_msg("Input stack                 : %s\n"%(stack))
		print_msg("Output directory            : %s\n"%(outdir))
		print_msg("Scratch directory           : %s\n"%(scratch))
		print_msg("Tangent low-pass filter   fl: %f\n"%(fl))
		print_msg("Tangent low-pass filter   aa: %f\n"%(aa))
		print_msg("Mask radius                 : %f\n"%(radccc))
		print_msg("Number of CPUs              : %i\n"%(ncpu))


	'''
       	if(len(scratch) > 4):
		if(scratch[:4] == "bdb:"):
			#scratch_file = scratch + "#diffproj%05d"%myid
			scratch_volume = scratch + "#gvbs%05d"%myid
		else:
			#scratch_file = os.path.join(scratch,"diffproj%05d.hdf"%myid)
			scratch_volume = os.path.join(scratch,"gvbs%05d.hdf"%myid)
	else:
		#scratch_file = os.path.join(scratch,"diffproj%05d.hdf"%myid)
		scratch_volume = os.path.join(scratch,"gvbs%05d.hdf"%myid)
	'''

# Step 0 - read input data
	if myid == main_node:
		start_time = time()
		if(file_type(stack) == "bdb"):
			from EMAN2db import db_open_dict
			dummy = db_open_dict(stack, True)
		
		# horatio active_refactoring Jy51i1EwmLD4tWZ9_00000_1
		# if(file_type(stack) == "bdb"):
		# 	from EMAN2db import db_open_dict
		# 	dummy = db_open_dict(stack, True)
		# active = EMUtil.get_all_attributes(stack, 'active')
		# list_of_particles = []
		# for im in xrange(len(active)):
		# 	if active[im]:  list_of_particles.append(im)
		# del active
		
		nima = EMUtil.get_image_count(stack)
		list_of_particles = range(nima)
		
	else:
		nima = 0
	total_nima = bcast_number_to_all(nima, source_node = main_node)

	if myid != main_node:
		list_of_particles = [-1]*total_nima
	list_of_particles = bcast_list_to_all(list_of_particles, source_node = main_node)

	image_start, image_end = MPI_start_end(total_nima, ncpu, myid)
	# create a list of images for each node
	list_of_particles = list_of_particles[image_start: image_end]
	nima = len(list_of_particles)
	if debug:
		finfo.write("image_start, image_end: %d %d\n" %(image_start, image_end))
		finfo.flush()

	data = EMData.read_images(stack, list_of_particles)
        nx = data[0].get_xsize()
        ny = data[0].get_ysize()
        nz = nx

	if myid== main_node:
		print_msg( "Time to read data: %d\n" % (time()-start_time) )
		start_time = time()

# Step 1 - Compute 3D
	if myid == main_node:  start_time = time()
	if  CTF:  
		from reconstruction import recons3d_4nn_ctf_MPI
		volr = recons3d_4nn_ctf_MPI(myid, data, snr = 1.0)
	else:     
		from reconstruction import recons3d_4nn_MPI
		volr = recons3d_4nn_MPI(myid, data)
	from sys import exit
	
	if myid == main_node:
		if(fl > 0.0):
			volr = filt_tanl( volr, fl, aa )
		volr.write_image(os.path.join(outdir,"avg1.hdf"),0)
	bcast_EMData_to_all(volr, myid, main_node)
	if myid== main_node:
		print_msg( "Time to compute 3D: %d\n" % (time()-start_time) )
		start_time = time()

# Step 2 - compute squared CTF
	if(radccc < 1):  radcir = min(nx,ny)//2-2
	else:            radcir = radccc
	mask2D = model_circle( radcir, nx, ny)

	from projection import prgs, prep_vol
	if CTF:
		from filter import filt_ctf
		from fundamentals import fftip, fft
		from morphology  import  ctf_img
		fftip(volr)
		#  Compute sum_CTF^2+1.0/snr  in 2D
		sumCTF2 = EMData(nx,ny,1, False)
		previous_defocus = -1.0
		for k in xrange(nima):
			ctf_params = data[k].get_attr( "ctf" )
			if ctf_params.defocus != previous_defocus:
				previous_defocus = ctf_params.defocus
				vb = ctf_img(nx, ctf_params,1, ny)#, nz)
				Util.mul_img(vb, vb)
			Util.add_img(sumCTF2, vb)
		reduce_EMData_to_root(sumCTF2, myid)
		if myid == 0:  sumCTF2 += (1.0/snr)
		bcast_EMData_to_all( sumCTF2, myid )
	else:
		vol1,kb = prep_vol(volr)
		del volr
	if myid== main_node:
		print_msg( "Time to compute CTF^2: %d\n" % (time()-start_time) )
		start_time = time()

# Step 3 - backproject weighted difference projections
	if myid == main_node:  start_time = time()
	nsym = 1
	nimages = total_nima
	ntripletsWnsym = nsym*nimages
	dm=[0.0]*(9*ntripletsWnsym)
	ss=[0.0]*(6*ntripletsWnsym)
	if myid == main_node:
		t = EMUtil.get_all_attributes(stack,"xform.projection")
		#  Here only those that are on the list of particles should remain, do it later
		count = 0
		for i in xrange(nimages):
			d = t[i].get_params("spider")
			DMnSS = Util.CANG(d["phi"],d["theta"],d["psi"])
			dm[(count*9) :(count+1)*9] = DMnSS["DM"]
			ss[(count*6) :(count+1)*6] = DMnSS["SS"]
			count += 1
		del t, d
	dm = bcast_list_to_all(dm, source_node = main_node)
	ss = bcast_list_to_all(ss, source_node = main_node)
	avg1 = model_blank(nx,ny,nz)
	var1 = model_blank(nx,ny,nz)
	if  CTF:  previous_defocus = -1.0
	for k in xrange(nima):
		t = data[k].get_attr("xform.projection")
		dprm = t.get_params("spider")
		if CTF:
			ctf_params = data[k].get_attr( "ctf" )
			if ctf_params.defocus != previous_defocus:
				previous_defocus = ctf_params.defocus
				#  projections will be mulitplied by the CTF, so they can be compared with the data
				vol1, kb = prep_vol(fft(filt_ctf(volr, ctf_params)))
		proj = prgs(vol1, kb, [dprm["phi"],dprm["theta"],dprm["psi"],-dprm["tx"],-dprm["ty"]])
		if(fl > 0.0):
			dat = filt_tanl( data[k], fl, aa )
			#df = dat
			df,a , b = im_diff(proj, dat, mask2D)
		else:
			df,a , b = im_diff(proj, data[k], mask2D)
		df.set_attr("xform.projection",t)
		
		# horatio active_refactoring Jy51i1EwmLD4tWZ9_00000_1
		# df.set_attr("active",1)
		if  CTF:
			df = Util.divn_img(fft(df),sumCTF2)
			df = fft(filt_ctf(df, ctf_params))
		vb, data[k] = recons3d_swbp(df, list_of_particles[k], dm, ss, "general")
		#  Data contains weighted 2D diffprojdata
		#  If CTF, it is also multiplied by the CTF and divided by the sum_CTF^2
		data[k].set_attr("xform.projection",t)
		# horatio active_refactoring Jy51i1EwmLD4tWZ9_00000_1
		# data[k].set_attr("active",1)

		#vb, jojo = recons3d_swbp(data[k], list_of_particles[k], dm, ss, "general")
		#if(fl > 0.0):
		#	vb = filt_tanl(vb, fl, aa)
		Util.add_img(avg1, vb)
		Util.add_img2(var1 , vb)
	#del ss, vb, dprm, proj, df
	if myid== main_node:
		print_msg( "Time to backproject weighted difference projections: %d\n" % (time()-start_time) )
		start_time = time()

# Step 4 - compute average and variance of difference projections
	if myid == main_node:  start_time = time()
	reduce_EMData_to_root(avg1, myid)
	reduce_EMData_to_root(var1, myid)
	
	if( myid == 0):
		#avg1.write_image(os.path.join(outdir, "avgb.hdf"))
		#sumCTF2.write_image(os.path.join(outdir, "s2.hdf"))
		Util.mul_scalar(avg1, 1.0/float(total_nima))
		avg1.write_image(os.path.join(outdir,"avg2.hdf"),0)
		avg2 = Util.muln_img(avg1, avg1)
		Util.mul_scalar(avg2, float(total_nima))
		Util.sub_img(var1, avg2)
		Util.mul_scalar(var1, 1.0/float(total_nima-1) )
		var1.write_image(os.path.join(outdir,"var1.hdf"),0)
		#p1, p2, p3, p4 = Util.infomask(var1, None, True)
		#print_msg( "First variance statistics: ave, std dev, min, max:  %15.3g  %15.3g  %15.3g  %15.3g\n" % (p1, p2, p3, p4) )
	else:
		avg1 = model_blank(nx,ny,nz)
		var1 = model_blank(nx,ny,nz)
	if myid== main_node:
		print_msg( "Time to compute the variance: %d\n" % (time()-start_time) )
		start_time = time()
	#bcast_EMData_to_all( var1, myid )
	#del avg1

# Step 5 - divide difference projections by the standard deviation (square root of the variance)
	'''
	if myid == main_node:  start_time = time()

	from statistics import pcanalyzer
	if(myid == 0):
		if(pcamask != None):
			pcamask = get_im( pcamask )
		else:  pcamask = model_blank(nx,ny,nz,1.0)
	else:           pcamask = model_blank(nx,ny,nz)
	bcast_EMData_to_all(pcamask, myid)

       	if(len(scratch) > 4):
		if(scratch[:4] == "bdb:"):
			scratch_dir = scratch[4:]
		else:
			scratch_dir = scratch
	else:
		scratch_dir = scratch

	pcaer = pcanalyzer(pcamask, outdir, pcanvec, True, scratch_dir)

	var1 = square_root(var1)

	if( myid == 0 ):
		vb = backproject_swbp(data[0], list_of_particles[k], dm)
		
		#if(fl > 0.0):
		#	vb = filt_tanl(vb, fl, aa)
		vb = Util.divn_img(vb, var1)
		refstat = Util.infomask(vb, pcamask, True)
	else:             refstat = [0.0,0.0,0.0,0.0]
	refstat = mpi_bcast(refstat, 4, MPI_FLOAT, 0, MPI_COMM_WORLD)
	refstat = map(float, refstat)

	avg2 = model_blank(nx,ny,nz)

	for k in xrange(nima):
		vb = Util.divn_img( backproject_swbp(data[k], list_of_particles[k], dm) , var1)
		#if(fl > 0.0):
		#	vb = filt_tanl(vb, fl, aa)
		pc = Util.infomask(vb, pcamask, True)
		vb -= pc[0]
		vb *= (refstat[1]/pc[1])
		Util.add_img(avg2, vb)
		pcaer.insert(vb)
	del data,vb
	reduce_EMData_to_root(avg2, myid)
	if( myid == 0):
		Util.mul_scalar(avg2, 1.0/float(total_nima))
	else:   avg2 =  model_blank(nx,ny,nz)
	bcast_EMData_to_all( avg2, myid )
	if(myid == 0):  avg2.write_image(os.path.join(outdir,"avg2.hdf"))

	assert not(avg2 is None)
	if myid== main_node:
		print_msg( "Time to divide backprojected difference projections by the standard deviation: %d\n" % (time()-start_time) )
		start_time = time()

# Step 6 - do the PCA on standardized projections
	if myid == main_node:  start_time = time()
	
	pcaer.setavg( avg2 )

	eigs = pcaer.analyze()
	from mpi import mpi_barrier
	mpi_barrier(MPI_COMM_WORLD)
	if myid==0:
		eigfile = os.path.join(outdir, "eigvol.hdf")
		for i in xrange( len(eigs) ):
			eigs[i].write_image( eigfile, i )
	if myid== main_node:
		print_msg( "Time to do PCA: %d\n" % (time()-start_time) )
		start_time = time()
	# Clean up masked files
	os.system("rm -f "+ os.path.join(scratch_dir , "maskedimg%04d.bin" % myid ) )
	if myid == main_node: print_end_msg("3D PCA")
	'''
"""
#  This one worked and produced newt9 on 07/13/10.
def OKvar_mpi_new(stack, outdir, scratch, fl, aa, radccc, writelp, writestack, frepa = "default", pca=False, pcamask=None, pcanvec=None, CTF = False):
	from string       import atoi, replace, split, atof
	from utilities    import get_im, circumference, model_circle, model_blank
	from utilities    import bcast_EMData_to_all, reduce_EMData_to_root
	from filter       import filt_tanl
	from mpi          import mpi_comm_rank, mpi_comm_size, mpi_barrier, mpi_bcast, mpi_reduce
	from mpi          import MPI_COMM_WORLD, MPI_INT, MPI_FLOAT, MPI_SUM
	from utilities    import bcast_list_to_all, bcast_number_to_all, reduce_EMData_to_root, bcast_EMData_to_all 
	from utilities    import send_attr_dict
	from utilities    import get_params_proj, file_type, print_msg, print_begin_msg, print_end_msg
	from applications import MPI_start_end
	from projection   import prep_vol
	from statistics   import im_diff
	from morphology   import square_root
	from time import time	
	import os
	from sys import exit

        myid = mpi_comm_rank( MPI_COMM_WORLD )
        ncpu = mpi_comm_size( MPI_COMM_WORLD )
	main_node = 0
	nx = 0
        if myid==0:
		if os.path.exists(outdir):
			nx = 1
			ERROR('Output directory exists, please change the name and restart the program', " var_mpi", 0)
		else:
			os.system( "mkdir " + outdir )
	nx = mpi_bcast(nx, 1, MPI_INT, 0, MPI_COMM_WORLD)
	nx = int(nx[0])
	if(nx != 0):
		from sys import exit
		exit()

	mpi_barrier( MPI_COMM_WORLD )
	debug = False  ############################################

	if myid == main_node:
		print_begin_msg("3D PCA")
		print_msg("Input stack                 : %s\n"%(stack))
		print_msg("Output directory            : %s\n"%(outdir))
		print_msg("Scratch directory           : %s\n"%(scratch))
		print_msg("Tangent low-pass filter   fl: %f\n"%(fl))
		print_msg("Tangent low-pass filter   aa: %f\n"%(aa))
		print_msg("Mask radius                 : %f\n"%(radccc))
		print_msg("Number of CPUs              : %i\n"%(ncpu))


	'''
       	if(len(scratch) > 4):
		if(scratch[:4] == "bdb:"):
			#scratch_file = scratch + "#diffproj%05d"%myid
			scratch_volume = scratch + "#gvbs%05d"%myid
		else:
			#scratch_file = os.path.join(scratch,"diffproj%05d.hdf"%myid)
			scratch_volume = os.path.join(scratch,"gvbs%05d.hdf"%myid)
	else:
		#scratch_file = os.path.join(scratch,"diffproj%05d.hdf"%myid)
		scratch_volume = os.path.join(scratch,"gvbs%05d.hdf"%myid)
	'''

# Step 0 - read input data
	if myid == main_node:
		start_time = time()
       		if(file_type(stack) == "bdb"):
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
	total_nima = bcast_number_to_all(nima, source_node = main_node)

	if myid != main_node:
		list_of_particles = [-1]*total_nima
	list_of_particles = bcast_list_to_all(list_of_particles, source_node = main_node)

	image_start, image_end = MPI_start_end(total_nima, ncpu, myid)
	# create a list of images for each node
	list_of_particles = list_of_particles[image_start: image_end]
	nima = len(list_of_particles)
	if debug:
		finfo.write("image_start, image_end: %d %d\n" %(image_start, image_end))
		finfo.flush()

	data = EMData.read_images(stack, list_of_particles)
	if myid== main_node:
		print_msg( "Time to read data: %d\n" % (time()-start_time) )
		start_time = time()

# Step 1 - Compute 3D
	if myid == main_node:  start_time = time()
	from reconstruction import recons3d_4nn_MPI, recons3d_swbp, backproject_swbp
	if  CTF:  vol1 = recons3d_4nn_MPI(myid, data)#  FIX
	else:
		vol1 = recons3d_4nn_MPI(myid, data)
	from sys import exit
	
	if myid == main_node:
		if(fl > 0.0):
			vol1 = filt_tanl( vol1, fl, aa )
		vol1.write_image(os.path.join(outdir,"avg1.hdf"),0)
	bcast_EMData_to_all(vol1, myid, main_node)
	from projection import prgs, prep_vol
	vol1,kb = prep_vol(vol1)
	if myid== main_node:
		print_msg( "Time to compute 3D: %d\n" % (time()-start_time) )
		start_time = time()

# Step 2 - compute difference projections
	nx = data[0].get_xsize()
	ny = data[0].get_ysize()
	nz = nx
	if(radccc < 1):  radcir = min(nx,ny)//2-2
	else:            radcir = radccc
	mask2D = model_circle( radcir, nx, ny)

# Step 3 - backproject weighted difference projections
	if myid == main_node:  start_time = time()
	nsym = 1
	nimages = total_nima
	ntripletsWnsym = nsym*nimages
	dm=[0.0]*(9*ntripletsWnsym)
	ss=[0.0]*(6*ntripletsWnsym)
	if myid == main_node:
		t = EMUtil.get_all_attributes(stack,"xform.projection")
		#  Here only those that are on the list of particles should remain, do it later
		count = 0
		for i in xrange(nimages):
			d = t[i].get_params("spider")
			DMnSS = Util.CANG(d["phi"],d["theta"],d["psi"])
			dm[(count*9) :(count+1)*9] = DMnSS["DM"]
			ss[(count*6) :(count+1)*6] = DMnSS["SS"]
			count += 1
		del t, d
	dm = bcast_list_to_all(dm, source_node = main_node)
	ss = bcast_list_to_all(ss, source_node = main_node)
	avg1 = model_blank(nx,ny,nz)
	var1 = model_blank(nx,ny,nz)
	for k in xrange(nima):
		t = data[k].get_attr("xform.projection")
		dprm = t.get_params("spider")
		proj = prgs(vol1, kb, [dprm["phi"],dprm["theta"],dprm["psi"],-dprm["tx"],-dprm["ty"]])
		if(fl > 0.0):
			dat = filt_tanl( data[k], fl, aa )
			df,a , b = im_diff(proj, dat, mask2D)
		else:
			df,a , b = im_diff(proj, data[k], mask2D)
		df.set_attr("xform.projection",t)
		df.set_attr("active",1)
		vb, data[k] = recons3d_swbp(df, list_of_particles[k], dm, ss, "general")
		#  Data contains weighted 2D diffprojdata
		data[k].set_attr("xform.projection",t)
		data[k].set_attr("active",1)
		#vb, jojo = recons3d_swbp(data[k], list_of_particles[k], dm, ss, "general")
		if(fl > 0.0):
			vb = filt_tanl(vb, fl, aa)
		Util.add_img(avg1, vb)
		Util.add_img2(var1 , vb)
	#del ss, vb, dprm, proj, df
	if myid== main_node:
		print_msg( "Time to backproject weighted difference projections: %d\n" % (time()-start_time) )
		start_time = time()
# Step 4 - compute average and variance of difference projections
	if myid == main_node:  start_time = time()
	reduce_EMData_to_root(avg1, myid)
	reduce_EMData_to_root(var1, myid)
	if( myid == 0):
		Util.mul_scalar(avg1, 1.0/float(total_nima))
		avg1.write_image(os.path.join(outdir,"avg2.hdf"),0)
		avg2 = Util.muln_img(avg1, avg1)
		Util.mul_scalar(avg2, float(total_nima))
		Util.sub_img(var1, avg2)
		Util.mul_scalar(var1, 1.0/float(total_nima-1) )
		p1, p2, p3, p4 = Util.infomask(var1, None, True)
		print_msg( "First variance statistics: ave, std dev, min, max:  %15.3g  %15.3g  %15.3g  %15.3g\n" % (p1, p2, p3, p4) )
	else:
		avg1 = model_blank(nx,ny,nz)
		var1 = model_blank(nx,ny,nz)

	#bcast_EMData_to_all( avg1, myid )
	bcast_EMData_to_all( var1, myid )
	del avg1

# Step 5 - divide difference projections by the standard deviation (square root of the variance)
	if myid == main_node:  start_time = time()

	if myid == main_node:  var1.write_image(os.path.join(outdir,"var1.hdf"),0)

	from statistics import pcanalyzer
	if(myid == 0):
		if(pcamask != None):
			pcamask = get_im( pcamask)
		else:  pcamask = model_blank(nx,ny,nz,1.0)
	else:           pcamask = model_blank(nx,ny,nz)
	bcast_EMData_to_all(pcamask, myid)


       	if(len(scratch) > 4):
		if(scratch[:4] == "bdb:"):
			scratch_dir = scratch[4:]
		else:
			scratch_dir = scratch
	else:
		scratch_dir = scratch

	pcaer = pcanalyzer(pcamask, outdir, pcanvec, True, scratch_dir)

	var1 = square_root(var1)

	if( myid == 0 ):
		vb = backproject_swbp(data[0], list_of_particles[k], dm)
		if(fl > 0.0):
			vb = filt_tanl(vb, fl, aa)
		vb = Util.divn_img(vb, var1)
		refstat = Util.infomask(vb, pcamask, True)
	else:             refstat = [0.0,0.0,0.0,0.0]
	refstat = mpi_bcast(refstat, 4, MPI_FLOAT, 0, MPI_COMM_WORLD)
	refstat = map(float, refstat)

	avg2 = model_blank(nx,ny,nz)

	for k in xrange(nima):
		vb = Util.divn_img(backproject_swbp(data[k], list_of_particles[k], dm), var1)
		if(fl > 0.0):
			vb = filt_tanl(vb, fl, aa)
		pc = Util.infomask(vb, pcamask, True)
		vb -= pc[0]
		vb *= (refstat[1]/pc[1])
		Util.add_img(avg2, vb)
		pcaer.insert(vb)
	del data,vb
	reduce_EMData_to_root(avg2, myid)
	if( myid == 0):
		Util.mul_scalar(avg2, 1.0/float(total_nima))
	else:   avg2 =  model_blank(nx,ny,nz)
	bcast_EMData_to_all( avg2, myid )
	if(myid == 0):  avg2.write_image(os.path.join(outdir,"avg2.hdf"))

	assert not(avg2 is None)
	if myid== main_node:
		print_msg( "Time to divide backprojected difference projections by the standard deviation: %d\n" % (time()-start_time) )
		start_time = time()

# Step 6 - do the PCA on standardized projections
	if myid == main_node:  start_time = time()
	
	pcaer.setavg( avg2 )

	eigs = pcaer.analyze()
	from mpi import mpi_barrier
	mpi_barrier(MPI_COMM_WORLD)
	if myid==0:
		eigfile = os.path.join(outdir, "eigvol.hdf")
		for i in xrange( len(eigs) ):
			eigs[i].write_image( eigfile, i )
	if myid== main_node:
		print_msg( "Time to do PCA: %d\n" % (time()-start_time) )
		start_time = time()
	# Clean up masked files
	os.system("rm -f "+ os.path.join(scratch_dir , "maskedimg%04d.bin" % myid ) )
	if myid == main_node: print_end_msg("3D PCA")

# In memory
def var_mpi_new___(stack, outdir, scratch, fl, aa, radccc, writelp, writestack, frepa = "default", pca=False, pcamask=None, pcanvec=None, CTF = False):
	from string       import atoi, replace, split, atof
	from utilities    import get_im, circumference, model_circle, model_blank
	from utilities    import bcast_EMData_to_all, reduce_EMData_to_root
	from filter       import filt_tanl
	from mpi          import mpi_comm_rank, mpi_comm_size, mpi_barrier, mpi_bcast, mpi_reduce
	from mpi          import MPI_COMM_WORLD, MPI_INT, MPI_FLOAT, MPI_SUM
	from utilities    import bcast_list_to_all, bcast_number_to_all, reduce_EMData_to_root, bcast_EMData_to_all 
	from utilities    import send_attr_dict
	from utilities    import get_params_proj, file_type, print_msg, print_begin_msg, print_end_msg
	from applications import MPI_start_end
	from projection   import prep_vol
	from statistics   import im_diff
	from morphology   import square_root
	from time import time	
	import os
	from sys import exit

        myid = mpi_comm_rank( MPI_COMM_WORLD )
        ncpu = mpi_comm_size( MPI_COMM_WORLD )
	main_node = 0
	nx = 0
        if myid==0:
		if os.path.exists(outdir):
			nx = 1
			ERROR('Output directory exists, please change the name and restart the program', " var_mpi", 0)
		else:
			os.system( "mkdir " + outdir )
	nx = mpi_bcast(nx, 1, MPI_INT, 0, MPI_COMM_WORLD)
	nx = int(nx[0])
	if(nx != 0):
		from sys import exit
		exit()

	mpi_barrier( MPI_COMM_WORLD )
	debug = False  ############################################

	if myid == main_node:
		print_begin_msg("3D PCA")
		print_msg("Input stack                 : %s\n"%(stack))
		print_msg("Output directory            : %s\n"%(outdir))
		print_msg("Scratch directory           : %s\n"%(scratch))
		print_msg("Tangent low-pass filter   fl: %f\n"%(fl))
		print_msg("Tangent low-pass filter   aa: %f\n"%(aa))
		print_msg("Mask radius                 : %f\n"%(radccc))
		print_msg("Number of CPUs              : %i\n"%(ncpu))


	'''
       	if(len(scratch) > 4):
		if(scratch[:4] == "bdb:"):
			#scratch_file = scratch + "#diffproj%05d"%myid
			scratch_volume = scratch + "#gvbs%05d"%myid
		else:
			#scratch_file = os.path.join(scratch,"diffproj%05d.hdf"%myid)
			scratch_volume = os.path.join(scratch,"gvbs%05d.hdf"%myid)
	else:
		#scratch_file = os.path.join(scratch,"diffproj%05d.hdf"%myid)
		scratch_volume = os.path.join(scratch,"gvbs%05d.hdf"%myid)
	'''

# Step 0 - read input data
	if myid == main_node:
		start_time = time()
       		if(file_type(stack) == "bdb"):
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
	total_nima = bcast_number_to_all(nima, source_node = main_node)

	if myid != main_node:
		list_of_particles = [-1]*total_nima
	list_of_particles = bcast_list_to_all(list_of_particles, source_node = main_node)

	image_start, image_end = MPI_start_end(total_nima, ncpu, myid)
	# create a list of images for each node
	list_of_particles = list_of_particles[image_start: image_end]
	nima = len(list_of_particles)
	if debug:
		finfo.write("image_start, image_end: %d %d\n" %(image_start, image_end))
		finfo.flush()

	data = EMData.read_images(stack, list_of_particles)
	if myid== main_node:
		print_msg( "Time to read data: %d\n" % (time()-start_time) )
		start_time = time()

	from reconstruction import recons3d_4nn_MPI, recons3d_swbp, backproject_swbp

# Step 2 - compute difference projections
	nx = data[0].get_xsize()
	ny = data[0].get_ysize()
	nz = nx
	if(radccc < 1):  radcir = min(nx,ny)//2-2
	else:            radcir = radccc
	mask2D = model_circle( radcir, nx, ny)

	if(myid == 0):
		if(pcamask != None):
			pcamask = get_im( pcamask)
		else:  pcamask = model_blank(nx,ny,nz,1.0)
	else:           pcamask = model_blank(nx,ny,nz)
	bcast_EMData_to_all(pcamask, myid)

# Step 3 - backproject weighted difference projections
	if myid == main_node:  start_time = time()
	nsym = 1
	nimages = total_nima
	ntripletsWnsym = nsym*nimages
	dm=[0.0]*(9*ntripletsWnsym)
	ss=[0.0]*(6*ntripletsWnsym)
	if myid == main_node:
		t = EMUtil.get_all_attributes(stack,"xform.projection")
		#  Here only those that are on the list of particles should remain, do it later
		count = 0
		for i in xrange(nimages):
			d = t[i].get_params("spider")
			DMnSS = Util.CANG(d["phi"],d["theta"],d["psi"])
			dm[(count*9) :(count+1)*9] = DMnSS["DM"]
			ss[(count*6) :(count+1)*6] = DMnSS["SS"]
			count += 1
		del t, d
	dm = bcast_list_to_all(dm, source_node = main_node)
	ss = bcast_list_to_all(ss, source_node = main_node)

	avg1 = model_blank(nx,ny,nz)
	var1 = model_blank(nx,ny,nz)
	for k in xrange(nima):
		'''
		t = data[k].get_attr("xform.projection")
		dprm = t.get_params("spider")
		proj = prgs(vol1, kb, [dprm["phi"],dprm["theta"],dprm["psi"],-dprm["tx"],-dprm["ty"]])
		if(fl > 0.0):
			dat = filt_tanl( data[k], fl, aa )
			df,a , b = im_diff(proj, dat, mask2D)
		else:
			df,a , b = im_diff(proj, data[k], mask2D)
		df.set_attr("xform.projection",t)
		df.set_attr("active",1)
		'''
		vb, data[k] = recons3d_swbp(data[k], list_of_particles[k], dm, ss, "general")
		'''
		#  Data contains weighted 2D projdata
		data[k].set_attr("xform.projection",t)
		data[k].set_attr("active",1)
		'''
		if(fl > 0.0):
			vb = filt_tanl(vb, fl, aa)
		Util.add_img(avg1, vb)
		Util.add_img2(var1 , vb)
	#del ss, vb, dprm, proj, df
	if myid == main_node:
		print_msg( "Time to backproject weighted difference projections: %d\n" % (time()-start_time) )
		start_time = time()

# Step 4 - compute average and variance of difference projections
	if myid == main_node:  start_time = time()
	reduce_EMData_to_root(avg1, myid)
	reduce_EMData_to_root(var1, myid)
	if( myid == 0):
		Util.mul_scalar(avg1, 1.0/float(total_nima))
		avg2 = Util.muln_img(avg1, avg1)
		Util.mul_scalar(avg2, float(total_nima))
		Util.sub_img(var1, avg2)
		Util.mul_scalar(var1, 1.0/float(total_nima-1) )
		p1, p2, p3, p4 = Util.infomask(var1, None, True)
		print_msg( "First variance statistics: ave, std dev, min, max:  %15.3g  %15.3g  %15.3g  %15.3g\n" % (p1, p2, p3, p4) )
	else:
		avg1 = model_blank(nx,ny,nz)
		var1 = model_blank(nx,ny,nz)

	bcast_EMData_to_all( avg1, myid )
	bcast_EMData_to_all( var1, myid )
	#del avg1

	from statistics import pcanalyzebck
	pcaer = pcanalyzebck(pcamask, pcanvec, data, list_of_particles, dm, var1, fl, aa, True)

# Step 5 - divide difference projections by the standard deviation (square root of the variance)
	if myid == main_node:  start_time = time()

	if myid == main_node:  var1.write_image(os.path.join(outdir,"var1.hdf"),0)
	'''
	var1 = square_root(var1)

	if( myid == 0 ):
		vb = backproject_swbp(data[0], list_of_particles[k], dm)
		if(fl > 0.0):
			vb = filt_tanl(vb, fl, aa)
		#vb = Util.divn_img(vb, var1)
		refstat = Util.infomask(vb, pcamask, True)
	else:             refstat = [0.0,0.0,0.0,0.0]
	refstat = mpi_bcast(refstat, 4, MPI_FLOAT, 0, MPI_COMM_WORLD)
	refstat = map(float, refstat)

	avg2 = model_blank(nx,ny,nz)

	for k in xrange(nima):
		#vb = Util.divn_img(backproject_swbp(data[k], list_of_particles[k], dm), var1)
		vb = backproject_swbp(data[k], list_of_particles[k], dm)
		if(fl > 0.0):
			vb = filt_tanl(vb, fl, aa)
		pc = Util.infomask(vb, pcamask, True)
		vb -= pc[0]
		vb *= (refstat[1]/pc[1])
		#Util.add_img(avg2, vb)
		pcaer.insert(vb)
	del data,vb
	#reduce_EMData_to_root(avg2, myid)
	#if( myid == 0):
	#	Util.mul_scalar(avg2, 1.0/float(total_nima))
	#else:   avg2 =  model_blank(nx,ny,nz)
	#bcast_EMData_to_all( avg2, myid )
	#if(myid == 0):  avg2.write_image(os.path.join(outdir,"avg2.hdf"))

	#assert not(avg2 is None)
	if myid== main_node:
		print_msg( "Time to divide backprojected difference projections by the standard deviation: %d\n" % (time()-start_time) )
		start_time = time()
	'''
# Step 6 - do the PCA on standardized projections
	if myid == main_node:  start_time = time()

	pcaer.setavg( avg1 )

	eigs = pcaer.analyze()

	if myid==0:
		eigfile = os.path.join(outdir, "eigvol.hdf")
		for i in xrange( len(eigs) ):
			eigs[i].write_image( eigfile, i )
	if myid== main_node:
		print_msg( "Time to do PCA: %d\n" % (time()-start_time) )
		start_time = time()
	if myid == main_node: print_end_msg("3D PCA")

# I do not know what this one was supposed to be....
def var_mpi_new_(stack, outdir, scratch, fl, aa, radccc, writelp, writestack, frepa = "default", pca=False, pcamask=None, pcanvec=None, CTF = False):
	from string       import atoi, replace, split, atof
	from utilities    import get_im, circumference, model_circle, model_blank
	from utilities    import bcast_EMData_to_all, reduce_EMData_to_root
	from filter       import filt_tanl
	from mpi          import mpi_comm_rank, mpi_comm_size, mpi_barrier, mpi_bcast, mpi_reduce
	from mpi          import MPI_COMM_WORLD, MPI_INT, MPI_FLOAT, MPI_SUM
	from utilities    import bcast_list_to_all, bcast_number_to_all, reduce_EMData_to_root, bcast_EMData_to_all 
	from utilities    import send_attr_dict
	from utilities    import get_params_proj, file_type, print_msg, print_begin_msg, print_end_msg
	from applications import MPI_start_end
	from projection   import prep_vol,prgs
	from statistics   import im_diff
	from morphology   import square_root
	from time import time	
	import os
	from sys import exit

        myid = mpi_comm_rank( MPI_COMM_WORLD )
        ncpu = mpi_comm_size( MPI_COMM_WORLD )
	main_node = 0
	nx = 0
        if myid==0:
		if os.path.exists(outdir):
			nx = 1
			ERROR('Output directory exists, please change the name and restart the program', " var_mpi", 0)
		else:
			os.system( "mkdir " + outdir )
	nx = mpi_bcast(nx, 1, MPI_INT, 0, MPI_COMM_WORLD)
	nx = int(nx[0])
	if(nx != 0):
		from sys import exit
		exit()

	mpi_barrier( MPI_COMM_WORLD )
	debug = False  ############################################

	if myid == main_node:
		print_begin_msg("3D PCA")
		print_msg("Input stack                 : %s\n"%(stack))
		print_msg("Output directory            : %s\n"%(outdir))
		print_msg("Scratch directory           : %s\n"%(scratch))
		print_msg("Tangent low-pass filter   fl: %f\n"%(fl))
		print_msg("Tangent low-pass filter   aa: %f\n"%(aa))
		print_msg("Mask radius                 : %f\n"%(radccc))
		print_msg("Number of CPUs              : %i\n"%(ncpu))


	'''
       	if(len(scratch) > 4):
		if(scratch[:4] == "bdb:"):
			#scratch_file = scratch + "#diffproj%05d"%myid
			scratch_volume = scratch + "#gvbs%05d"%myid
		else:
			#scratch_file = os.path.join(scratch,"diffproj%05d.hdf"%myid)
			scratch_volume = os.path.join(scratch,"gvbs%05d.hdf"%myid)
	else:
		#scratch_file = os.path.join(scratch,"diffproj%05d.hdf"%myid)
		scratch_volume = os.path.join(scratch,"gvbs%05d.hdf"%myid)
	'''
# Step 0 - read input data
	if myid == main_node:
		start_time = time()
       		if(file_type(stack) == "bdb"):
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
	total_nima = bcast_number_to_all(nima, source_node = main_node)

	if myid != main_node:
		list_of_particles = [-1]*total_nima
	list_of_particles = bcast_list_to_all(list_of_particles, source_node = main_node)

	image_start, image_end = MPI_start_end(total_nima, ncpu, myid)
	# create a list of images for each node
	list_of_particles = list_of_particles[image_start: image_end]
	nima = len(list_of_particles)
	if debug:
		finfo.write("image_start, image_end: %d %d\n" %(image_start, image_end))
		finfo.flush()

	data = EMData.read_images(stack, list_of_particles)
	if myid== main_node:
		print_msg( "Time to read data: %d\n" % (time()-start_time) )
		start_time = time()

	from reconstruction import recons3d_4nn_MPI, recons3d_swbp, backproject_swbp

# Step 2 - compute difference projections
	nx = data[0].get_xsize()
	ny = data[0].get_ysize()
	nz = nx
	if(radccc < 1):  radcir = min(nx,ny)//2-2
	else:            radcir = radccc
	mask2D = model_circle( radcir, nx, ny)



	from statistics import pcanalyzer
	if(myid == 0):
		if(pcamask != None):
			pcamask = get_im( pcamask)
		else:  pcamask = model_blank(nx,ny,nz,1.0)
	else:           pcamask = model_blank(nx,ny,nz)
	bcast_EMData_to_all(pcamask, myid)


       	if(len(scratch) > 4):
		if(scratch[:4] == "bdb:"):
			scratch_dir = scratch[4:]
		else:
			scratch_dir = scratch
	else:
		scratch_dir = scratch
	pcaer = pcanalyzer(pcamask, outdir, pcanvec, True, scratch_dir)

# Step 3 - backproject weighted difference projections
	if myid == main_node:  start_time = time()
	nsym = 1
	nimages = total_nima
	ntripletsWnsym = nsym*nimages
	dm=[0.0]*(9*ntripletsWnsym)
	ss=[0.0]*(6*ntripletsWnsym)
	if myid == main_node:
		t = EMUtil.get_all_attributes(stack,"xform.projection")
		#  Here only those that are on the list of particles should remain, do it later
		count = 0
		for i in xrange(nimages):
			d = t[i].get_params("spider")
			DMnSS = Util.CANG(d["phi"],d["theta"],d["psi"])
			dm[(count*9) :(count+1)*9] = DMnSS["DM"]
			ss[(count*6) :(count+1)*6] = DMnSS["SS"]
			count += 1
		del t, d
	dm = bcast_list_to_all(dm, source_node = main_node)
	ss = bcast_list_to_all(ss, source_node = main_node)

	avg1 = model_blank(nx,ny,nz)
	var1 = model_blank(nx,ny,nz)
	for k in xrange(nima):
		t = data[k].get_attr("xform.projection")
		dprm = t.get_params("spider")
		proj = prgs(vol1, kb, [dprm["phi"],dprm["theta"],dprm["psi"],-dprm["tx"],-dprm["ty"]])
		if(fl > 0.0):
			dat = filt_tanl( data[k], fl, aa )
			df,a , b = im_diff(proj, dat, mask2D)
		else:
			df,a , b = im_diff(proj, data[k], mask2D)
		df.set_attr("xform.projection",t)
		df.set_attr("active",1)
		vb, data[k] = recons3d_swbp(data[k], list_of_particles[k], dm, ss,"general")
		#  Data contains weighted 2D projdata
		data[k].set_attr("xform.projection",t)
		data[k].set_attr("active",1)
		if(fl > 0.0):
			vb = filt_tanl(vb, fl, aa)
		pcaer.insert(vb)
		Util.add_img(avg1, vb)
		Util.add_img2(var1 , vb)
	#del ss, vb, dprm, proj, df
	if myid== main_node:
		print_msg( "Time to backproject weighted difference projections: %d\n" % (time()-start_time) )
		start_time = time()

# Step 4 - compute average and variance of difference projections
	if myid == main_node:  start_time = time()
	reduce_EMData_to_root(avg1, myid)
	reduce_EMData_to_root(var1, myid)
	if( myid == 0):
		Util.mul_scalar(avg1, 1.0/float(total_nima))
		avg2 = Util.muln_img(avg1, avg1)
		Util.mul_scalar(avg2, float(total_nima))
		Util.sub_img(var1, avg2)
		Util.mul_scalar(var1, 1.0/float(total_nima-1) )
		p1, p2, p3, p4 = Util.infomask(var1, None, True)
		print_msg( "First variance statistics: ave, std dev, min, max:  %15.3g  %15.3g  %15.3g  %15.3g\n" % (p1, p2, p3, p4) )
	else:
		avg1 = model_blank(nx,ny,nz)
		var1 = model_blank(nx,ny,nz)

	bcast_EMData_to_all( avg1, myid )
	bcast_EMData_to_all( var1, myid )
	#del avg1

# Step 5 - divide difference projections by the standard deviation (square root of the variance)
	if myid == main_node:  start_time = time()

	if myid == main_node:  var1.write_image(os.path.join(outdir,"var1.hdf"),0)
	var1 = square_root(var1)

	if( myid == 0 ):
		vb = backproject_swbp(data[0], list_of_particles[k], dm)
		if(fl > 0.0):
			vb = filt_tanl(vb, fl, aa)
		#vb = Util.divn_img(vb, var1)
		refstat = Util.infomask(vb, pcamask, True)
	else:             refstat = [0.0,0.0,0.0,0.0]
	refstat = mpi_bcast(refstat, 4, MPI_FLOAT, 0, MPI_COMM_WORLD)
	refstat = map(float, refstat)

	avg2 = model_blank(nx,ny,nz)

	for k in xrange(nima):
		#vb = Util.divn_img(backproject_swbp(data[k], list_of_particles[k], dm), var1)
		vb = backproject_swbp(data[k], list_of_particles[k], dm)
		if(fl > 0.0):
			vb = filt_tanl(vb, fl, aa)
		pc = Util.infomask(vb, pcamask, True)
		vb -= pc[0]
		vb *= (refstat[1]/pc[1])
		#Util.add_img(avg2, vb)
		pcaer.insert(vb)
	del data,vb
	#reduce_EMData_to_root(avg2, myid)
	#if( myid == 0):
	#	Util.mul_scalar(avg2, 1.0/float(total_nima))
	#else:   avg2 =  model_blank(nx,ny,nz)
	#bcast_EMData_to_all( avg2, myid )
	#if(myid == 0):  avg2.write_image(os.path.join(outdir,"avg2.hdf"))

	#assert not(avg2 is None)
	if myid== main_node:
		print_msg( "Time to divide backprojected difference projections by the standard deviation: %d\n" % (time()-start_time) )
		start_time = time()
# Step 6 - do the PCA on standardized projections
	if myid == main_node:  start_time = time()
	
	pcaer.setavg( avg1 )

	eigs = pcaer.analyze()

	if myid==0:
		eigfile = os.path.join(outdir, "eigvol.hdf")
		for i in xrange( len(eigs) ):
			eigs[i].write_image( eigfile, i )
	if myid== main_node:
		print_msg( "Time to do PCA: %d\n" % (time()-start_time) )
		start_time = time()

	# Clean up masked files
	os.system("rm -f "+ os.path.join(scratch_dir , "maskedimg%04d.bin" % myid ) )
	if myid == main_node: print_end_msg("3D PCA")
"""

'''
###############################################################################################################################################
#Ran's test for ali3d
def MPI_start_end_rantest(nima, nproc, myid):
	image_start = int(round(float(nima)/nproc*myid))
	image_end   = int(round(float(nima)/nproc*(myid+1)))
	return image_start, image_end

'''

"""	
def ali3d_rantest_old(stack, ref_vol, outdir, maskfile = None, ir = 1, ou = -1, rs = 1, 
            xr = "4 2 2 1", yr = "-1", ts = "1 1 0.5 0.25", delta = "10 6 4 4", an = "-1",
	    center = -1, maxit = 5, CTF = False, snr = 1.0,  ref_a = "S", sym = "c1",  user_func_name = "ref_ali3d",
	    fourvar = True, debug = False):

	from alignment       import Numrinit, prepare_refrings, proj_ali_incore, proj_ali_incore_local
	from utilities       import model_circle, get_image, drop_image, get_input_from_string
	from utilities       import bcast_list_to_all, bcast_number_to_all, reduce_EMData_to_root, bcast_EMData_to_all 
	from utilities       import send_attr_dict
	from utilities       import get_params_proj, file_type
	from fundamentals    import rot_avg_image
	import os
	import types
	from utilities       import print_begin_msg, print_end_msg, print_msg
	from mpi             import mpi_bcast, mpi_comm_size, mpi_comm_rank, MPI_FLOAT, MPI_COMM_WORLD, mpi_barrier
	from mpi             import mpi_reduce, MPI_INT, MPI_SUM
	from filter          import filt_ctf
	from projection      import prep_vol, prgs
	from statistics      import hist_list, varf3d_MPI


	number_of_proc = mpi_comm_size(MPI_COMM_WORLD)
	myid           = mpi_comm_rank(MPI_COMM_WORLD)
	main_node = 0
	
	if os.path.exists(outdir):  ERROR('Output directory exists, please change the name and restart the program', "ali3d_rantest", 1, myid)
	mpi_barrier(MPI_COMM_WORLD)
	
	if myid == main_node:
		print_begin_msg("ali3d_rantest")
		os.mkdir(outdir)
	mpi_barrier(MPI_COMM_WORLD)

	if debug:
		from time import sleep
		while not os.path.exists(outdir):
			print  "Node ",myid,"  waiting..."
			sleep(5)

		info_file = os.path.join(outdir, "progress%04d"%myid)
		finfo = open(info_file, 'w')
	else:
		finfo = None

	xrng        = get_input_from_string(xr)
	if  yr == "-1":  yrng = xrng
	else          :  yrng = get_input_from_string(yr)
	step        = get_input_from_string(ts)
	delta       = get_input_from_string(delta)
	lstp = min(len(xrng), len(yrng), len(step), len(delta))
	if an == "-1":
		an = [-1] * lstp
	else:
		an = get_input_from_string(an)

	first_ring  = int(ir)
	rstep       = int(rs)
	last_ring   = int(ou)
	max_iter    = int(maxit)
	center      = int(center)

	vol     = EMData()
	vol.read_image(ref_vol)
	nx      = vol.get_xsize()
	if last_ring < 0:	last_ring = int(nx/2) - 2

	if myid == main_node:
		import user_functions
		user_func = user_functions.factory[user_func_name]

		print_msg("Input stack                 : %s\n"%(stack))
		print_msg("Reference volume            : %s\n"%(ref_vol))	
		print_msg("Output directory            : %s\n"%(outdir))
		print_msg("Maskfile                    : %s\n"%(maskfile))
		print_msg("Inner radius                : %i\n"%(first_ring))
		print_msg("Outer radius                : %i\n"%(last_ring))
		print_msg("Ring step                   : %i\n"%(rstep))
		print_msg("X search range              : %s\n"%(xrng))
		print_msg("Y search range              : %s\n"%(yrng))
		print_msg("Translational step          : %s\n"%(step))
		print_msg("Angular step                : %s\n"%(delta))
		print_msg("Angular search range        : %s\n"%(an))
		print_msg("Maximum iteration           : %i\n"%(max_iter))
		print_msg("Center type                 : %i\n"%(center))
		print_msg("CTF correction              : %s\n"%(CTF))
		print_msg("Signal-to-Noise Ratio       : %f\n"%(snr))
		print_msg("Reference projection method : %s\n"%(ref_a))
		print_msg("Symmetry group              : %s\n\n"%(sym))

	if maskfile:
		if type(maskfile) is types.StringType: mask3D = get_image(maskfile)
		else:                                  mask3D = maskfile
	else: mask3D = model_circle(last_ring, nx, nx, nx)

	numr	= Numrinit(first_ring, last_ring, rstep, "F")
	mask2D  = model_circle(last_ring,nx,nx) - model_circle(first_ring,nx,nx)

	fscmask = model_circle(last_ring,nx,nx,nx)
	if CTF:
		from reconstruction import rec3D_MPI
		from filter         import filt_ctf
	else:	 from reconstruction import rec3D_MPI_noCTF

	'''
	if myid == main_node:
       		if(file_type(stack) == "bdb"):
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
	total_nima = bcast_number_to_all(nima, source_node = main_node)

	if myid != main_node:
		list_of_particles = [-1]*total_nima
	list_of_particles = bcast_list_to_all(list_of_particles, source_node = main_node)

	image_start, image_end = MPI_start_end_rantest(total_nima, number_of_proc, myid)
	# create a list of images for each node
	list_of_particles = list_of_particles[image_start: image_end]
	nima = len(list_of_particles)
	'''
	if(file_type(stack) == "bdb"):
		from EMAN2db import db_open_dict
		dummy = db_open_dict(stack, True)
	list_of_particles = range(EMUtil.get_image_count(stack))
	if myid==main_node:
		nima=EMUtil.get_image_count(stack)
	else:
		nima=0
	total_nima = bcast_number_to_all(nima, source_node = main_node)	
	
	if myid != main_node:
		list_of_particles = [-1]*total_nima
	list_of_particles = bcast_list_to_all(list_of_particles, source_node = main_node)

	image_start, image_end = MPI_start_end_rantest(total_nima, number_of_proc, myid)
	# create a list of images for each node
	list_of_particles = list_of_particles[image_start: image_end]
	nima = len(list_of_particles)

	
	if debug:
		finfo.write("image_start, image_end: %d %d\n" %(image_start, image_end))
		finfo.flush()

	data = EMData.read_images(stack, list_of_particles)
	if fourvar:  original_data = []
	for im in xrange(nima):
		data[im].set_attr('ID', list_of_particles[im])
		if fourvar: original_data.append(data[im].copy())
		if CTF:
			ctf_params = data[im].get_attr("ctf")
			st = Util.infomask(data[im], mask2D, False)
			data[im] -= st[0]
			data[im] = filt_ctf(data[im], ctf_params)
			data[im].set_attr('ctf_applied', 1)

	if debug:
		finfo.write( '%d loaded  \n' % nima )
		finfo.flush()
	if myid == main_node:
		# initialize data for the reference preparation function
		ref_data = [ mask3D, max(center,0), None, None, None, None ]
		# for method -1, switch off centering in user function

	from time import time	

	#  this is needed for gathering of pixel errors
	disps = []
	recvcount = []
	for im in xrange(number_of_proc):
		if( im == main_node ):  disps.append(0)
		else:                  disps.append(disps[im-1] + recvcount[im-1])
		ib, ie = MPI_start_end_rantest(total_nima, number_of_proc, im)
		recvcount.append( ie - ib )

	pixer = [0.0]*nima
	cs = [0.0]*3
	total_iter = 0
	# do the projection matching
	for N_step in xrange(lstp):
		terminate = 0
		Iter = -1
 		while(Iter < max_iter-1 and terminate == 0):
			print "loop2"
			Iter += 1
			total_iter += 1
			if myid == main_node:
				start_time = time()
				print_msg("\nITERATION #%3d,  inner iteration #%3d\nDelta = %4.1f, an = %5.2f, xrange = %5.2f, yrange = %5.2f, step = %5.2f\n"%(total_iter, Iter, delta[N_step], an[N_step], xrng[N_step],yrng[N_step],step[N_step]))

			volft,kb = prep_vol( vol )
			refrings = prepare_refrings( volft, kb, nx, delta[N_step], ref_a, sym, numr, True)
			del volft,kb
			if myid== main_node:
				print_msg( "Time to prepare rings: %d\n" % (time()-start_time) )
				start_time = time()

			for im in xrange( nima ):

				if an[N_step] == -1:
					peak, pixer[im] = proj_ali_incore(data[im],refrings,numr,xrng[N_step],yrng[N_step],step[N_step],finfo)
				else:
					peak, pixer[im] = proj_ali_incore_local(data[im],refrings,numr,xrng[N_step],yrng[N_step],step[N_step],an[N_step],finfo)

			if myid == main_node:
				print_msg("Time of alignment = %d\n"%(time()-start_time))
				start_time = time()

			#output pixel errors
			from mpi import mpi_gatherv
			recvbuf = mpi_gatherv(pixer, nima, MPI_FLOAT, recvcount, disps, MPI_FLOAT, main_node, MPI_COMM_WORLD)
			mpi_barrier(MPI_COMM_WORLD)
			terminate = 0
			if(myid == main_node):
				recvbuf = map(float, recvbuf)
				from statistics import hist_list
				lhist = 20
				region, histo = hist_list(recvbuf, lhist)
				if(region[0] < 0.0):  region[0] = 0.0
				msg = "      Histogram of pixel errors\n      ERROR       number of particles\n"
				print_msg(msg)
				for lhx in xrange(lhist):
					msg = " %10.3f     %7d\n"%(region[lhx], histo[lhx])
					print_msg(msg)
				# Terminate if 95% within 1 pixel error
				im = 0
				for lhx in xrange(lhist):
					if(region[lhx] > 1.0): break
					im += histo[lhx]
				if(im/float(total_nima) > 0.95):  terminate = 1
				del region, histo
			del recvbuf
			terminate = mpi_bcast(terminate, 1, MPI_INT, 0, MPI_COMM_WORLD)
			terminate = int(terminate[0])

			if(center == -1):
				from utilities      import estimate_3D_center_MPI, rotate_3D_shift
				cs[0], cs[1], cs[2], dummy, dummy = estimate_3D_center_MPI(data, total_nima, myid, number_of_proc, main_node)				
				if myid == main_node:
					msg = " Average center x = %10.3f        Center y = %10.3f        Center z = %10.3f\n"%(cs[0], cs[1], cs[2])
					print_msg(msg)
				cs = mpi_bcast(cs, 3, MPI_FLOAT, main_node, MPI_COMM_WORLD)
				cs = [-float(cs[0]), -float(cs[1]), -float(cs[2])]
				rotate_3D_shift(data, cs)

			if CTF: vol, fscc = rec3D_MPI(data, snr, sym, fscmask, os.path.join(outdir, "resolution%04d"%(total_iter)), myid, main_node)
			else:    vol, fscc = rec3D_MPI_noCTF(data, sym, fscmask, os.path.join(outdir, "resolution%04d"%(total_iter)), myid, main_node)

			if myid == main_node:
				print_msg("3D reconstruction time = %d\n"%(time()-start_time))
				start_time = time()

			if fourvar:
			#  Compute Fourier variance
				for im in xrange(nima):
					original_data[im].set_attr( 'xform.projection', data[im].get_attr('xform.projection') )
				varf = varf3d_MPI(original_data, ssnr_text_file = os.path.join(outdir, "ssnr%04d"%(total_iter)), mask2D = None, reference_structure = vol, ou = last_ring, rw = 1.0, npad = 1, CTF = CTF, sign = 1, sym =sym, myid = myid)
				if myid == main_node:
					print_msg("Time to calculate 3D Fourier variance= %d\n"%(time()-start_time))
					start_time = time()
					varf = 1.0/varf
			else:  varf = None

			if myid == main_node:
				drop_image(vol, os.path.join(outdir, "vol%04d.hdf"%(total_iter)))
				ref_data[2] = vol
				ref_data[3] = fscc
				ref_data[4] = varf
				#  call user-supplied function to prepare reference image, i.e., center and filter it
				vol, cs = user_func(ref_data)
				drop_image(vol, os.path.join(outdir, "volf%04d.hdf"%(total_iter)))

			del varf
			bcast_EMData_to_all(vol, myid, main_node)
			# write out headers, under MPI writing has to be done sequentially
			mpi_barrier(MPI_COMM_WORLD)
			par_str = ['xform.projection', 'ID']
	        	if myid == main_node:
	        		if(file_type(stack) == "bdb"):
	        			#from utilities import recv_attr_dict_bdb
					print "rantest"
	        			recv_attr_dict_rantest(main_node, stack, data, par_str, image_start, image_end, number_of_proc, "log.txt")
					#write_attr_rantest(main_node, stack, data, par_str, image_start, image_end, number_of_proc,outdir)
	        		else:
	        			from utilities import recv_attr_dict
	        			recv_attr_dict(main_node, stack, data, par_str, image_start, image_end, number_of_proc)
	        	else:	       send_attr_dict(main_node, data, par_str, image_start, image_end)
	if myid == main_node: print_end_msg("ali3d_MPI")

def recv_attr_dict_rantest(main_node, stack, data, list_params, image_start, image_end, number_of_proc, fexport, comm = -1):
	import types
	from  utilities import  get_arb_params, set_arb_params, write_text_row
	from  mpi 	import mpi_recv
	from  mpi 	import MPI_FLOAT, MPI_INT, MPI_TAG_UB, MPI_COMM_WORLD
	#  bdb version!
	# This is done on the main node, so for images from the main node, simply write headers
	
	if comm == -1: comm = MPI_COMM_WORLD
	
	DB = db_open_dict(stack)
	fexp=open(fexport,'w')
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
			print "hello"
			header = mpi_recv(len_list*(img_end-img_begin), MPI_FLOAT, n, MPI_TAG_UB, comm)
			for im in xrange(img_begin, img_end):
				par_begin = (im-img_begin)*len_list
				#nvalue = []
				ilis = 0
				for il in xrange(len(list_params)):
					if(ink[il] == 1):
						#nvalue.append(int(header[par_begin+ilis]))
						nvalue= int(header[par_begin+ilis])
						ilis += 1
						fexp.write("%15s   \n"%str(nvalue))
					elif ink[il]==0:
						#nvalue.append(float(header[par_begin+ilis]))
						nvalue=float(header[par_begin+ilis])
						ilis += 1
						fexp.write("%15s   \n"%str(nvalue))
					else:
						assert ink[il]==2
						t = Transform()
						tmp = []
						for iii in xrange(par_begin+ilis, par_begin+ilis+12):
							tmp.append(float(header[iii]))
						t.set_matrix(tmp)
						print t, "Heloo"
						if list_params[il] == "xform.align2d":
							d = t.get_params("2D")
							fexp.write("%15.5f %15.5f %15.5f %10d %10.3f \n"%(d["alpha"],d["tx"],d["ty"],d["mirror"],d["scale"]))
						elif list_params[il] == "xform.projection":
							d = t.get_params("spider")
							fexp.write("%15.5f %15.5f %15.5f %15.5f %15.5f \n"%(d["phi"],d["theta"],d["psi"],-d["tx"],-d["ty"]))
							
						elif list_params[il] == "xform.align3d":
							d = t.get_params("spider")
							fexp.write("%15.5f %15.5f %15.5f %15.5f %15.5f %15.5f %10d %10.3f \n"%(d["phi"],d["theta"],d["psi"],d["tx"],d["ty"],d["tz"],d["mirror"],d["scale"]))

						elif list_params[il] == "ctf":
							fexp.write("%15.5f %15.5f %15.5f %15.5f %15.5f %15.5f \n"%(t.defocus, t.cs, t.voltage, t.apix, t.bfactor, t.ampcont))
						ilis += 12	
		else:
			for n in xrange(image_start, image_end):
				for param in list_params:
					t=data[n-image_start].get_attr(param)
					if param == "xform.align2d":
						d = t.get_params("2D")
						fexp.write("%15.5f %15.5f %15.5f %10d %10.3f \n"%(d["alpha"],d["tx"],d["ty"],d["mirror"],d["scale"]))
					elif param == "xform.projection":
						d = t.get_params("spider")
						fexp.write("%15.5f %15.5f %15.5f %15.5f %15.5f \n"%(d["phi"],d["theta"],d["psi"],-d["tx"],-d["ty"]))
							
					elif param == "xform.align3d":
						d = t.get_params("spider")
						fexp.write("%15.5f %15.5f %15.5f %15.5f %15.5f %15.5f %10d %10.3f \n"%(d["phi"],d["theta"],d["psi"],d["tx"],d["ty"],d["tz"],d["mirror"],d["scale"]))

					elif param == "ctf":
						fexp.write("%15.5f %15.5f %15.5f %15.5f %15.5f %15.5f \n"%(t.defocus, t.cs, t.voltage, t.apix, t.bfactor, t.ampcont))
							

					else:
						fexp.write("%15s   \n"%str(t))	
					
	fexp.close()				
	DB.close()
def ali3d_rantest(stack, ref_vol, outdir, maskfile = None, ir = 1, ou = -1, rs = 1, 
            xr = "4 2 2 1", yr = "-1", ts = "1 1 0.5 0.25", delta = "10 6 4 4", an = "-1", deltapsi = "-1", startpsi = "-1",
	    center = -1, maxit = 5, CTF = False, snr = 1.0,  ref_a = "S", sym = "c1",  user_func_name = "ref_ali3d",
	    fourvar = True, npad = 4, debug = False, termprec = 0.0):

	from alignment       import Numrinit, prepare_refrings, proj_ali_incore, proj_ali_incore_local
	from utilities       import model_circle, get_image, drop_image, get_input_from_string
	from utilities       import bcast_list_to_all, bcast_number_to_all, reduce_EMData_to_root, bcast_EMData_to_all
	from utilities       import send_attr_dict
	from utilities       import get_params_proj, file_type
	from fundamentals    import rot_avg_image
	import os
	import types
	from applications     import MPI_start_end
	from utilities       import print_begin_msg, print_end_msg, print_msg
	from mpi             import mpi_bcast, mpi_comm_size, mpi_comm_rank, MPI_FLOAT, MPI_COMM_WORLD, mpi_barrier
	from mpi             import mpi_reduce, MPI_INT, MPI_SUM, mpi_recv, mpi_send, MPI_TAG_UB
	from filter          import filt_ctf
	from projection      import prep_vol, prgs
	from statistics      import hist_list, varf3d_MPI



	number_of_proc = mpi_comm_size(MPI_COMM_WORLD)
	myid           = mpi_comm_rank(MPI_COMM_WORLD)
	main_node = 0
	
	if os.path.exists(outdir):  ERROR('Output directory exists, please change the name and restart the program', "ali3d_MPI", 1, myid)
	mpi_barrier(MPI_COMM_WORLD)
	
	if myid == main_node:
		print_begin_msg("ali3d_MPI")
		os.mkdir(outdir)
	mpi_barrier(MPI_COMM_WORLD)

	if debug:
		from time import sleep
		while not os.path.exists(outdir):
			print  "Node ",myid,"  waiting..."
			sleep(5)

		info_file = os.path.join(outdir, "progress%04d"%myid)
		finfo = open(info_file, 'w')
	else:
		finfo = None

	xrng        = get_input_from_string(xr)
	if  yr == "-1":  yrng = xrng
	else          :  yrng = get_input_from_string(yr)
	step        = get_input_from_string(ts)
	delta       = get_input_from_string(delta)
	lstp = min(len(xrng), len(yrng), len(step), len(delta))
	if an == "-1":
		an = [-1] * lstp
	else:
		an = get_input_from_string(an)

	if deltapsi == "-1":
		deltapsi = [-1] * lstp
	else:
		deltapsi = get_input_from_string(deltapsi)

	if startpsi == "-1":
		startpsi = [-1] * lstp
	else:
		startpsi = get_input_from_string(startpsi)

	first_ring  = int(ir)
	rstep       = int(rs)
	last_ring   = int(ou)
	max_iter    = int(maxit)
	center      = int(center)

	vol     = EMData()
	vol.read_image(ref_vol)
	nx      = vol.get_xsize()
	if last_ring < 0:	last_ring = int(nx/2) - 2

	if myid == main_node:
		import user_functions
		user_func = user_functions.factory[user_func_name]

		print_msg("Input stack                 : %s\n"%(stack))
		print_msg("Reference volume            : %s\n"%(ref_vol))	
		print_msg("Output directory            : %s\n"%(outdir))
		print_msg("Maskfile                    : %s\n"%(maskfile))
		print_msg("Inner radius                : %i\n"%(first_ring))
		print_msg("Outer radius                : %i\n"%(last_ring))
		print_msg("Ring step                   : %i\n"%(rstep))
		print_msg("X search range              : %s\n"%(xrng))
		print_msg("Y search range              : %s\n"%(yrng))
		print_msg("Translational step          : %s\n"%(step))
		print_msg("Angular step                : %s\n"%(delta))
		print_msg("Angular search range        : %s\n"%(an))
		print_msg("Delta psi                   : %s\n"%(deltapsi))
		print_msg("Start psi                   : %s\n"%(startpsi))
		print_msg("Maximum iteration           : %i\n"%(max_iter))
		print_msg("Percentage of change for termination: %f\n"%(termprec))
		print_msg("Center type                 : %i\n"%(center))
		print_msg("CTF correction              : %s\n"%(CTF))
		print_msg("Signal-to-Noise Ratio       : %f\n"%(snr))
		print_msg("Reference projection method : %s\n"%(ref_a))
		print_msg("Symmetry group              : %s\n\n"%(sym))

	if maskfile:
		if type(maskfile) is types.StringType: mask3D = get_image(maskfile)
		else:                                  mask3D = maskfile
	else: mask3D = model_circle(last_ring, nx, nx, nx)

	numr	= Numrinit(first_ring, last_ring, rstep, "F")
	mask2D  = model_circle(last_ring,nx,nx) - model_circle(first_ring,nx,nx)

	fscmask = model_circle(last_ring,nx,nx,nx)
	if CTF:
		from reconstruction import rec3D_MPI
		from filter         import filt_ctf
	else:	 from reconstruction import rec3D_MPI_noCTF
	
	if(file_type(stack) == "bdb"):
		from EMAN2db import db_open_dict
		dummy = db_open_dict(stack, True)
	if myid==main_node:
		nima=EMUtil.get_image_count(stack)
	else:
		nima=0
		
	total_nima = bcast_number_to_all(nima, source_node = main_node)	
	image_start, image_end = MPI_start_end_rantest(total_nima, number_of_proc, myid)
	nima = image_end-image_start
	

	if debug:
		finfo.write("image_start, image_end: %d %d\n" %(image_start, image_end))
		finfo.flush()
	
	list_of_img=range(image_start,image_end)
	
	data = EMData.read_images(stack, list_of_img)
	if fourvar:  original_data = []
	for im in xrange(nima):
		data[im].set_attr('ID', list_of_img[im])
		if fourvar: original_data.append(data[im].copy())
		if CTF:
			ctf_params = data[im].get_attr("ctf")
			st = Util.infomask(data[im], mask2D, False)
			data[im] -= st[0]
			data[im] = filt_ctf(data[im], ctf_params)
			data[im].set_attr('ctf_applied', 1)

	if debug:
		finfo.write( '%d loaded  \n' % nima )
		finfo.flush()
	if myid == main_node:
		# initialize data for the reference preparation function
		ref_data = [ mask3D, max(center,0), None, None, None, None ]
		# for method -1, switch off centering in user function

	from time import time	

	#  this is needed for gathering of pixel errors
	disps = []
	recvcount = []
	for im in xrange(number_of_proc):
		if( im == main_node ):  disps.append(0)
		else:                  disps.append(disps[im-1] + recvcount[im-1])
		ib, ie = MPI_start_end(total_nima, number_of_proc, im)
		recvcount.append( ie - ib )

	pixer = [0.0]*nima
	cs = [0.0]*3
	total_iter = 0
	# do the projection matching
	for N_step in xrange(lstp):
		terminate = 0
		Iter = -1
 		while(Iter < max_iter-1 and terminate == 0):
			Iter += 1
			total_iter += 1
			if myid == main_node:
				start_time = time()
				print_msg("\nITERATION #%3d,  inner iteration #%3d\nDelta = %4.1f, an = %5.2f, xrange = %5.2f, yrange = %5.2f, step = %5.2f\n"%(total_iter, Iter, delta[N_step], an[N_step], xrng[N_step],yrng[N_step],step[N_step]))

			volft,kb = prep_vol( vol )
			refrings = prepare_refrings( volft, kb, nx, delta[N_step], ref_a, sym, numr, True)
			del volft,kb
			if myid== main_node:
				print_msg( "Time to prepare rings: %d\n" % (time()-start_time) )
				start_time = time()

			for im in xrange( nima ):

				if deltapsi[N_step] > 0.0:
					from alignment import proj_ali_incore_delta
					peak, pixer[im] = proj_ali_incore_delta(data[im],refrings,numr,xrng[N_step],yrng[N_step],step[N_step],startpsi[N_step],deltapsi[N_step],finfo)						
				elif an[N_step] == -1:
					peak, pixer[im] = proj_ali_incore(data[im],refrings,numr,xrng[N_step],yrng[N_step],step[N_step],finfo)
				else:
					peak, pixer[im] = proj_ali_incore_local(data[im],refrings,numr,xrng[N_step],yrng[N_step],step[N_step],an[N_step],finfo)

			if myid == main_node:
				print_msg("Time of alignment = %d\n"%(time()-start_time))
				start_time = time()

			#output pixel errors
			from mpi import mpi_gatherv
			recvbuf = mpi_gatherv(pixer, nima, MPI_FLOAT, recvcount, disps, MPI_FLOAT, main_node, MPI_COMM_WORLD)
			mpi_barrier(MPI_COMM_WORLD)
			terminate = 0
			if(myid == main_node):
				recvbuf = map(float, recvbuf)
				from statistics import hist_list
				lhist = 20
				region, histo = hist_list(recvbuf, lhist)
				if(region[0] < 0.0):  region[0] = 0.0
				msg = "      Histogram of pixel errors\n      ERROR       number of particles\n"
				print_msg(msg)
				for lhx in xrange(lhist):
					msg = " %10.3f     %7d\n"%(region[lhx], histo[lhx])
					print_msg(msg)
				# Terminate if 95% within 1 pixel error
				im = 0
				for lhx in xrange(lhist):
					if(region[lhx] > 1.0): break
					im += histo[lhx]
				precn = 100*float(total_nima-im)/float(total_nima)
				msg = " Number of particles that changed orientations %7d, percentage of total: %5.1f\n"%(total_nima-im, precn)
				print_msg(msg)
				if(precn <= termprec):  terminate = 1
				del region, histo
			del recvbuf
			terminate = mpi_bcast(terminate, 1, MPI_INT, 0, MPI_COMM_WORLD)
			terminate = int(terminate[0])

			if(center == -1):
				from utilities      import estimate_3D_center_MPI, rotate_3D_shift
				cs[0], cs[1], cs[2], dummy, dummy = estimate_3D_center_MPI(data, total_nima, myid, number_of_proc, main_node)				
				if myid == main_node:
					msg = " Average center x = %10.3f        Center y = %10.3f        Center z = %10.3f\n"%(cs[0], cs[1], cs[2])
					print_msg(msg)
				cs = mpi_bcast(cs, 3, MPI_FLOAT, main_node, MPI_COMM_WORLD)
				cs = [-float(cs[0]), -float(cs[1]), -float(cs[2])]
				rotate_3D_shift(data, cs)

			if CTF: vol, fscc = rec3D_MPI(data, snr, sym, fscmask, os.path.join(outdir, "resolution%04d"%(total_iter)), myid, main_node, npad = npad)
			else:    vol, fscc = rec3D_MPI_noCTF(data, sym, fscmask, os.path.join(outdir, "resolution%04d"%(total_iter)), myid, main_node, npad = npad)

			if myid == main_node:
				print_msg("3D reconstruction time = %d\n"%(time()-start_time))
				start_time = time()

			if fourvar:
			#  Compute Fourier variance
				for im in xrange(nima):
					original_data[im].set_attr( 'xform.projection', data[im].get_attr('xform.projection') )
				varf = varf3d_MPI(original_data, ssnr_text_file = os.path.join(outdir, "ssnr%04d"%(total_iter)), mask2D = None, reference_structure = vol, ou = last_ring, rw = 1.0, npad = 1, CTF = CTF, sign = 1, sym =sym, myid = myid)
				if myid == main_node:
					print_msg("Time to calculate 3D Fourier variance= %d\n"%(time()-start_time))
					start_time = time()
					varf = 1.0/varf
			else:  varf = None

			if myid == main_node:
				drop_image(vol, os.path.join(outdir, "vol%04d.hdf"%(total_iter)))
				ref_data[2] = vol
				ref_data[3] = fscc
				ref_data[4] = varf
				#  call user-supplied function to prepare reference image, i.e., center and filter it
				vol, cs = user_func(ref_data)
				drop_image(vol, os.path.join(outdir, "volf%04d.hdf"%(total_iter)))

			del varf
			bcast_EMData_to_all(vol, myid, main_node)
			# write out headers, under MPI writing has to be done sequentially
			mpi_barrier(MPI_COMM_WORLD)
			par_str = ['xform.projection']
			m=5
			if myid == main_node:
				fexp=open(os.path.join(outdir, "Iter_record_%04d_%04d.txt"%(N_step,Iter)),"w")
				for n in xrange(number_of_proc):
					if n!=main_node:
						t=mpi_recv(recvcount[n]*m,MPI_FLOAT, n, MPI_TAG_UB,MPI_COMM_WORLD)
						for i in xrange(recvcount[n]):
							for j in xrange(m):
								fexp.write("%15.5f"%t[j+i*m])
							fexp.write("\n")
					else:
						for i in xrange(recvcount[0]):
							phi, theta, psi, s2x, s2y=get_params_proj(data[i])
							fexp.write("%15.5f %15.5f %15.5f %15.5f %15.5f "%(phi, theta, psi, s2x, s2y))
							fexp.write("\n")
				fexp.close()								
	        	else:	       
				nvalue=[0]*5*(image_end-image_start)
				for i in xrange(image_end-image_start):
					phi, theta, psi, s2x, s2y=get_params_proj(data[i])
					nvalue[5*i + 0] = phi
					nvalue[5*i + 1] = theta
					nvalue[5*i + 2] = psi
					nvalue[5*i + 3] = s2x
					nvalue[5*i + 4] = s2y
				mpi_send(nvalue, len(nvalue), MPI_FLOAT,main_node,MPI_TAG_UB, MPI_COMM_WORLD)
				del nvalue
	if myid == main_node:
	   	if(file_type(stack) == "bdb"):
	        	from utilities import recv_attr_dict_bdb
	        	recv_attr_dict_bdb(main_node, stack, data, par_str, image_start, image_end, number_of_proc)
	        else:
	        	from utilities import recv_attr_dict
	        	recv_attr_dict(main_node, stack, data, par_str, image_start, image_end, number_of_proc)
		print_msg("Time to write header information= %d\n"%(time()-start_time))
		start_time = time()
	else:	       send_attr_dict(main_node, data, par_str, image_start, image_end)

	if myid == main_node: print_end_msg("ali3d_MPI")
"""


def recons3d_n_new(prj_stack, pid_list, vol_stack, CTF=False, snr=1.0, sign=1, npad=4, sym="c1", listfile = "", group = -1, verbose=0, MPI=False,rx=1,ry=1):
       	if MPI:
		#print("Calling recons3d_n_MPI")
		recons3d_n_MPI_new(prj_stack, pid_list, vol_stack, CTF, snr, 1, npad, sym, listfile, group, verbose,rx,ry)
		return
	#if MPI:
		#recons3d_n_MPI_new(prj_stack, pid_list, vol_stack, CTF, snr, 1, npad, sym, listfile, group, verbose)
		#return
	#	exit()
	from reconstruction import recons3d_4nn_ctf, recons3d_4nn
	from utilities import drop_image
	from utilities import print_begin_msg, print_end_msg, print_msg

	print_begin_msg("recons3d_nn4_rect")
	print_msg("Input stack                 : %s\n"%(prj_stack))
	print_msg("Output volume               : %s\n"%(vol_stack))
	print_msg("Padding factor              : %i\n"%(npad))
	print_msg("CTF correction              : %s\n"%(CTF))
	print_msg("Signal-to-Noise Ratio       : %f\n"%(snr))
	print_msg("CTF sign                    : %i\n"%(sign))
	print_msg("Symmetry group              : %s\n\n"%(sym))
	if(listfile):
		from utilities import read_text_file
		pid_list = read_text_file(listfile, 0)
		pid_list = map(int, pid_list)
		print_msg("Reconstruction for images listed in file : %s\n\n"%(listfile))
	elif(group > -1):
		print_msg("Reconstruction for group             : %i\n\n"%(group))
		tmp_list = EMUtil.get_all_attributes(prj_stack, 'group')
		pid_list = []
		for i in xrange(len(tmp_list)):
			if(tmp_list[i] == group):  pid_list.append(i)
		del tmp_list
	if CTF: vol = recons3d_4nn_ctf_new(prj_stack,rx,ry, pid_list, snr, 1, sym, verbose, npad)
	else:   vol = recons3d_4nn_new(prj_stack,rx,ry, pid_list, sym, npad)
	
	if(vol_stack[-3:] == "spi"):
		drop_image(vol, vol_stack, "s")
	else:
		drop_image(vol, vol_stack)
	print_end_msg("recons3d_n")

def recons3d_4nn_new(stack_name,rx,ry,list_proj=[], symmetry="c1", npad=4, snr=None, weighting=1, varsnr=True):
	"""
	Perform a 3-D reconstruction using Pawel's FFT Back Projection algoritm.
	   
	Input:
	   stack_name - name of the file with projection data.
	   
	   list_proj -  list of projections to be used in the reconstruction

	   symmetry - Point group of the target molecule (defaults to "C1")
	   
	   npad - 

	   Angles and shifts are passed in the file header as set_attr. Keywords are phi, theta, psi, sx, sy

	   Return:  3D reconstructed volume image

	   Usage:
	     vol = recons3d_4nn(filepattern, list_proj, symmetry)
	"""
	import types
	if list_proj == []:
		if type(stack_name) == types.StringType: nima = EMUtil.get_image_count(stack_name)
		else : nima = len(stack_name)
		list_proj = xrange(nima) 
	# read first image to determine the size to use
	if type(stack_name) == types.StringType:
		proj = EMData()
		proj.read_image(stack_name, list_proj[0])
	else:    proj = stack_name[list_proj[0]].copy()

	size = proj.get_xsize()
	# sanity check -- image must be square
	if size != proj.get_ysize():
		ERROR("input data has to be square","recons3d_4nn",1)
	# reconstructor
	
	if snr is None:
		params = {"sizeprojection":size, "npad":npad, "symmetry":symmetry, "weighting":weighting,"xratio":rx,"yratio":ry}
	else:
		params = {"sizeprojection":size, "npad":npad, "symmetry":symmetry, "weighting":weighting, "snr":snr, "varsnr":int(varsnr)}
	
	r = Reconstructors.get("nn4_rect", params)
	r.setup()

	if type(stack_name) == types.StringType:
		for i in xrange(len(list_proj)):
			proj.read_image(stack_name, list_proj[i])
			# horatio active_refactoring Jy51i1EwmLD4tWZ9_00000_1
			# active = proj.get_attr_default('active', 1)
			# if(active == 1):
			# 	xform_proj = proj.get_attr( "xform.projection" )
			# 	r.insert_slice(proj, xform_proj )
		
			xform_proj = proj.get_attr( "xform.projection" )
			r.insert_slice(proj, xform_proj )
	else:
		for i in list_proj:
			# horatio active_refactoring Jy51i1EwmLD4tWZ9_00000_1
			# active = stack_name[i].get_attr_default('active', 1)
			# if(active == 1):
			# 	xform_proj = stack_name[i].get_attr( "xform.projection" )
			# 	r.insert_slice(stack_name[i], xform_proj )
			xform_proj = stack_name[i].get_attr( "xform.projection" )
			r.insert_slice(stack_name[i], xform_proj )
	return r.finish(True)

def recons3d_4nn_ctf_new(stack_name,rx,ry,list_proj = [], snr = 10.0, sign=1, symmetry="c1", verbose=0, npad=4):
	"""Perform a 3-D reconstruction using Pawel's FFT Back Projection algoritm.
	   
	   Input:
	    stack_name - name of the stack file on a disk,
	                 each image has to have the following attributes set:
			 psi, theta, phi, sx, sy, defocus, 
	    list_proj - list of images from stack_name to be included in the reconstruction
	    symmetry	 -- Point group of the target molecule (defaults to "C1")

	   Return:  3d reconstructed volume image

	   Usage:
	     
	     anglelist = getAngles("myangles.txt") # not yet written
	     vol = do_reconstruction(filepattern, start, end, anglelist, symmetry)
	"""
	print "nn4_ctf_rect reconstructor is called"
	print "rx==%f,ry==%f\t"%(rx,ry)
	import types
	# read first image to determine the size to use
	if list_proj == []:	
		if type(stack_name) == types.StringType: nima = EMUtil.get_image_count(stack_name)
		else : nima = len(stack_name)
		list_proj = xrange(nima) 
	# read first image to determine the size to use
	if type(stack_name) == types.StringType:
		proj = EMData()
		proj.read_image(stack_name, list_proj[0])
	else:    proj = stack_name[list_proj[0]].copy()
	
	# convert angles to transform (rotation) objects
	# horatio active_refactoring Jy51i1EwmLD4tWZ9_00000_1
	# active = proj.get_attr_default('active', 1)
	size   = proj.get_xsize()
	"""

	"""
	# reconstructor
	#params = {"size":size, "npad":npad, "symmetry":symmetry, "snr":snr, "sign":sign}
	params = {"sizeprojection":size, "npad":npad, "symmetry":symmetry, "snr":snr, "sign":sign,"xratio":rx,"yratio":ry}
	r = Reconstructors.get("nn4_ctf_rect", params)
	r.setup()


	if type(stack_name) == types.StringType:
		for i in xrange(len(list_proj)):
		#for i in range(1):
			#print "i=======%d\t"%(i)
			# horatio active_refactoring Jy51i1EwmLD4tWZ9_00000_1
			# proj.read_image(stack_name, list_proj[i])
			# active = proj.get_attr_default('active', 1)
			# if(active == 1):
			# 	xform_proj = proj.get_attr( "xform.projection" )
			# 	r.insert_slice(proj, xform_proj )
			proj.read_image(stack_name, list_proj[i])
			xform_proj = proj.get_attr( "xform.projection" )
			r.insert_slice(proj, xform_proj )
	else:
		for i in xrange(len(list_proj)):
			# horatio active_refactoring Jy51i1EwmLD4tWZ9_00000_1
			# active = stack_name[list_proj[i]].get_attr_default('active', 1)
			# if(active == 1):
			# 	xform_proj = stack_name[list_proj[i]].get_attr( "xform.projection" )
			# 	r.insert_slice(stack_name[list_proj[i]], xform_proj )
			xform_proj = stack_name[list_proj[i]].get_attr( "xform.projection" )
			r.insert_slice(stack_name[list_proj[i]], xform_proj )
	return r.finish(True)


def recons3d_n_MPI_new(prj_stack, pid_list, vol_stack, CTF, snr, sign, npad, sym, listfile, group, verbose,rx,ry):
	#from reconstruction import recons3d_4nn_ctf_MPI 
	#recons3d_4nn_MPI_new
	from applications import MPI_start_end
	from utilities      import get_im, drop_image, bcast_number_to_all
	from utilities      import print_begin_msg, print_end_msg, print_msg
	from string         import replace
	from time           import time
	from mpi 	    import mpi_comm_size, mpi_comm_rank, mpi_bcast, MPI_INT, MPI_COMM_WORLD

	myid  = mpi_comm_rank(MPI_COMM_WORLD)
	nproc = mpi_comm_size(MPI_COMM_WORLD)
	time_start = time()
	
	if(myid == 0):
		print_begin_msg("recons3d_n_MPI_new")
		print_msg("Input stack  	       : %s\n"%(prj_stack))
		print_msg("Output volume	       : %s\n"%(vol_stack))
		print_msg("Padding factor	       : %i\n"%(npad))
		print_msg("CTF correction	       : %s\n"%(CTF))
		print_msg("Signal-to-Noise Ratio       : %f\n"%(snr))
		print_msg("CTF sign		       : %i\n"%(sign))
		print_msg("Symmetry group	       : %s\n\n"%(sym))
		
		#print "rx==%f,ry==%f\t"%(rx,ry)
		if(listfile):
			from utilities import read_text_file
			pid_list = read_text_file(listfile, 0)
			pid_list = map(int, pid_list)
			print_msg("Reconstruction for images listed in file : %s\n\n"%(listfile))
		elif(group > -1):
			print_msg("Reconstruction for group		: %i\n\n"%(group))
			tmp_list = EMUtil.get_all_attributes(prj_stack, 'group')
			pid_list = []
			for i in xrange(len(tmp_list)):
				if(tmp_list[i] == group):  pid_list.append(i)
			del tmp_list
		nima = len(pid_list)
	else:
		nima = 0

	nima = bcast_number_to_all(nima, source_node = 0)

	if(listfile or group > -1):
		if myid != 0:
			pid_list = [-1]*nima
		pid_list = mpi_bcast(pid_list, nima, MPI_INT, 0, MPI_COMM_WORLD)
		pid_list = map(int, pid_list)
	else:
		pid_list = range(nima)

	if verbose==0:
		finfo = None
	else:
		infofile = "progress%04d.txt"%(myid+1)
		finfo = open( infofile, 'w' )

	image_start, image_end = MPI_start_end(nima, nproc, myid)

	prjlist = EMData.read_images(prj_stack, pid_list[image_start:image_end])
	del pid_list

	if CTF: vol = recons3d_4nn_ctf_MPI_new(myid, prjlist, snr, sign, sym, finfo, npad,rx,ry)
	else:	vol = recons3d_4nn_MPI_new(myid, prjlist, sym, finfo, npad,rx,ry)
	if myid == 0 :
		if(vol_stack[-3:] == "spi"):
			drop_image(vol, vol_stack, "s")
		else:
			drop_image(vol, vol_stack)
		if not(finfo is None):
			finfo.write( "result written to " + vol_stack + "\n")
			finfo.write( "Total time: %10.3f\n" % (time()-time_start) )
			finfo.flush()

def recons3d_4nn_MPI_new(myid, prjlist, symmetry="c1", info=None, npad = 4,rx=1,ry=1):
	from utilities import reduce_EMData_to_root
	
	if( len(prjlist) == 0 ):
		ERROR("empty input list","recons3d_4nn_MPI",1)

	imgsize = prjlist[0].get_xsize()
	imgsizex= prjlist[0].get_xsize()
	imgsizey= prjlist[0].get_ysize()
	if prjlist[0].get_ysize() != imgsize:
		ERROR("input data has to be square","recons3d_4nn_MPI",1)
	

	fftvol = EMData()
	weight = EMData()

	#params = {"size":imgsize, "npad":npad, "symmetry":symmetry, "fftvol":fftvol, "weight":weight}
	params = {"sizeprojection":imgsize, "npad":npad, "symmetry":symmetry, "fftvol":fftvol,"weight":weight,"xratio":rx,"yratio":ry}
	r = Reconstructors.get( "nn4_rect", params )
	r.setup()

	if( not (info is None) ): nimg = 0
	for prj in prjlist :
		if prj.get_xsize() != imgsize or prj.get_ysize() != imgsize:
			ERROR("inconsistent image size","recons3d_4nn_MPI",1)

		# horatio active_refactoring Jy51i1EwmLD4tWZ9_00000_1
		# active = prj.get_attr_default('active', 1)
		# if(active == 1):
		# 	xform_proj = prj.get_attr( "xform.projection" )
		# 	r.insert_slice(prj, xform_proj )
		xform_proj = prj.get_attr( "xform.projection" )
		r.insert_slice(prj, xform_proj )
		if( not (info is None) ):
			nimg += 1
			info.write("Image %4d inserted.\n" %(nimg) )
			info.flush()

	if not (info is None): 
		info.write( "Begin reducing ...\n" )
		info.flush()

	reduce_EMData_to_root(fftvol, myid)
	reduce_EMData_to_root(weight, myid)
	sizex=int(rx*imgsize)
	sizey=int(ry*imgsize)
		
	if(myid==0):
		print "recon 4nn mpi myid==%drx==%f,ry==%fpx=%dpy=%d\t"%(myid,rx,ry,imgsizex,imgsizey)
		
	if myid == 0 :  
		vol = r.finish(True)
		print "volume size%d %d %d \t"%(vol.get_xsize(),vol.get_ysize(),vol.get_zsize())
	else:
		from utilities import model_blank
		vol = model_blank(sizex,sizey,imgsize)
	return vol

def recons3d_4nn_ctf_MPI_new(myid, prjlist, snr, sign=1, symmetry="c1", info=None, npad = 4,rx=1,ry=1):
	"""
		recons3d_4nn_ctf - calculate CTF-corrected 3-D reconstruction from a set of projections using three Eulerian angles, two shifts, and CTF settings for each projeciton image
		Input
			stack: name of the stack file containing projection data, projections have to be squares
			list_proj: list of projections to be included in the reconstruction
			snr: Signal-to-Noise Ratio of the data 
			sign: sign of the CTF 
			symmetry: point-group symmetry to be enforced, each projection will enter the reconstruction in all symmetry-related directions.
	"""
	#if(myid==0):
		#print "recon_ctf 4nn mpi myid==%drx==%f,ry==%f\t"%(myid,rx,ry)
	#print "3d_4nn_ctf_MPI_new  myid==%drx==%f,ry==%f\t"%(myid,rx,ry)
	from utilities import reduce_EMData_to_root
	if( len(prjlist) == 0 ):
	    ERROR("empty input list","recons3d_4nn_ctf_MPI",1)

	imgsize = prjlist[0].get_xsize()
	if prjlist[0].get_ysize() != imgsize:
		ERROR("input data has to be square","recons3d_4nn_ctf_MPI_new",1)


	fftvol = EMData()
	weight = EMData()

	#params = {"size":imgsize, "npad":npad, "snr":snr, "sign":sign, "symmetry":symmetry, "fftvol":fftvol, "weight":weight}
	params = {"sizeprojection":imgsize, "npad":npad, "snr":snr, "sign":sign, "symmetry":symmetry, "fftvol":fftvol, "weight":weight,"xratio":rx,"yratio":ry}
	r = Reconstructors.get( "nn4_ctf_rect", params )
	r.setup()

	if( not (info is None) ): nimg = 0
	for prj in prjlist :
		if prj.get_xsize() != imgsize or prj.get_ysize() != imgsize:
			ERROR("inconsistent image size","recons3d_4nn_ctf_MPI",1)

		# horatio active_refactoring Jy51i1EwmLD4tWZ9_00000_1
		# active = prj.get_attr_default('active', 1)
		# if(active == 1):
		# 	xform_proj = prj.get_attr( "xform.projection" )
		# 	r.insert_slice(prj, xform_proj )
		xform_proj = prj.get_attr( "xform.projection" )
		r.insert_slice(prj, xform_proj )
		if( not (info is None) ):
			nimg += 1
			info.write(" %4d inserted\n" %(nimg) )
			info.flush()

	if( not (info is None) ): 
		info.write( "begin reduce\n" )
		info.flush()

	reduce_EMData_to_root(fftvol, myid)
	reduce_EMData_to_root(weight, myid)
	
	sizex=int(rx*imgsize)
	sizey=int(ry*imgsize)

	if( not (info is None) ): 
		info.write( "after reduce\n" )
		info.flush()

	if myid == 0 :
		vol = r.finish(True)
		print "recon_ctf 4nn mpi myid==%drx==%f,ry==%fsizex=%dsizey=%d\t"%(myid,rx,ry,sizex,sizey)
	else:
		from utilities import model_blank
		vol = model_blank(sizex,sizey,imgsize)

	return vol

def project3d_new(volume, stack, mask = None, delta = 5, method = "S", phiEqpsi = "Minus", symmetry = "c1", listagls = None , listctfs = None, noise = None):
	from projection    import   prgs
	#prep_vol
	from utilities     import   even_angles, read_text_row, set_params_proj, model_gauss_noise, info
	from string        import   split
	from filter        import   filt_ctf,filt_gaussl
	import os
	import types

	if listagls is None:
		angles = even_angles(delta, symmetry = symmetry, method = method, phiEqpsi = phiEqpsi)
	elif(type(listagls) is types.StringType):
		angles = read_text_row(listagls, "", "")
	else:
		angles = listagls

	# try to parse the CTFs list. this is either not set (None), a filename or a list of values
	if listctfs is None:
		# not set, so simply ignore it 
		pass
	elif (type(listctfs) is types.StringType):
		# a string, so assume this is a filename and try to open the file
		try:
			ctfs = read_text_row(listctfs, "", "")
		except:
			ctfs = [None for ii in xrange(len(angles))]
	else:
		# assume this a list of len(angles)
		ctfs = listctfs

	if not noise is None:
		# try to convert noise string to float. ignore noise if this fails
		try:
			noise_level = float(noise)
		except:
			noise_level = None
	# ignore noise, since it was not requested
	else:
		noise_level = None

	if(type(volume) is types.StringType):
		vol = EMData()
		vol.read_image(volume)
		if(mask):
			if(type(mask) is types.StringType):
				maski = EMData()
				maski.read_image(volume)
				Util.mul_img(vol, maski)
				del maski
			else:
				Util.mul_img(vol, mask)
		volft, kbx, kby, kbz = prep_vol_new(vol)
	else:
		if(mask):
			if(type(mask) is types.StringType):
				maski = EMData()
				maski.read_image(volume)
				Util.mul_img(vol, maski)
				del maski
			else:
				Util.mul_img(vol, mask)
		volft, kbx, kby, kbz = prep_vol_new(volume)


	if(type(stack) is types.StringType):
		Disk = True
		os.system("rm -f  "+stack)	
	else:
		out = []
		Disk = False
	
	s2x=0
	s2y=0
	for i in xrange(len(angles)):
	#for i in xrange(1):
		#print "i==%d\t"%(i)
		if(len(angles[i]) == 3):
			proj = prgs_new(volft, kbx, kby, kbz, [angles[i][0], angles[i][1], angles[i][2], 0.0, 0.0])
			set_params_proj(proj, [angles[i][0], angles[i][1], angles[i][2], 0.0, 0.0])
		else:
			proj = prgs_new(volft, kbx, kby, kbz, [angles[i][0], angles[i][1], angles[i][2], -angles[i][3], -angles[i][4]])
			set_params_proj(proj, angles[i])

		# horatio active_refactoring Jy51i1EwmLD4tWZ9_00000_1
		# proj.set_attr_dict({'active':1})

		# add noise, if noise is set. this is two-fold: application of noise before
		#    ctf filtering and after it.
		if noise is not None:
			try:
				# no mask, so call w/ false
				noise_ima = model_gauss_noise(noise_level,proj.get_xsize(),
							      proj.get_ysize())
			except:
				pass
			else:
				proj += noise_ima

		# apply ctf, if ctf option is set and if we can create a valid CTF object
		try:
			ctf = EMAN2Ctf()
			# order of values is the one applied in sxheader for import / export!
			ctf.from_dict({ "defocus":ctfs[i][0], "cs":ctfs[i][1], "voltage":ctfs[i][2], 
					"apix":ctfs[i][3], "bfactor":ctfs[i][4], "ampcont":ctfs[i][5] })
		except:
			# there are no ctf values, so ignore this and set no values
			proj.set_attr( "error",1)
		else:
			# setting of values worked, so apply ctf and set the header info correctly
			proj = filt_ctf(proj,ctf)
			proj.set_attr( "ctf",ctf)
			proj.set_attr( "ctf_applied",0)

		# add second noise level that is not affected by CTF
		if noise is not None:
			try:
				noise_ima = model_gauss_noise(noise_level,proj.get_xsize(),
							      proj.get_ysize())
			except:
				pass
			else:
				proj += noise_ima

		if(Disk):
			proj.write_image(stack, i)
		else: 
			out.append(proj)
	if(not Disk):  return out
	
def project3d_new_fast(volume, stack, mask = None, delta = 5, method = "S", phiEqpsi = "Minus", symmetry = "c1", listagls = None , listctfs = None, noise = None):
	from projection    import   prgs
	#prep_vol
	from utilities     import   even_angles, read_text_row, set_params_proj, model_gauss_noise, info
	from string        import   split
	from filter        import   filt_ctf,filt_gaussl
	import os
	import types

	if listagls is None:
		angles = even_angles(delta, symmetry = symmetry, method = method, phiEqpsi = phiEqpsi)
	elif(type(listagls) is types.StringType):
		angles = read_text_row(listagls, "", "")
	else:
		angles = listagls

	# try to parse the CTFs list. this is either not set (None), a filename or a list of values
	if listctfs is None:
		# not set, so simply ignore it 
		pass
	elif (type(listctfs) is types.StringType):
		# a string, so assume this is a filename and try to open the file
		try:
			ctfs = read_text_row(listctfs, "", "")
		except:
			ctfs = [None for ii in xrange(len(angles))]
	else:
		# assume this a list of len(angles)
		ctfs = listctfs

	if not noise is None:
		# try to convert noise string to float. ignore noise if this fails
		try:
			noise_level = float(noise)
		except:
			noise_level = None
	# ignore noise, since it was not requested
	else:
		noise_level = None

	if(type(volume) is types.StringType):
		vol = EMData()
		vol.read_image(volume)
		sizez=vol.get_zsize()
		if(mask):
			if(type(mask) is types.StringType):
				maski = EMData()
				maski.read_image(volume)
				Util.mul_img(vol, maski)
				del maski
			else:
				Util.mul_img(vol, mask)
		volft, kbx, kby, kbz = prep_vol_new(vol)
	else:
		if(mask):
			if(type(mask) is types.StringType):
				maski = EMData()
				maski.read_image(volume)
				Util.mul_img(vol, maski)
				del maski
			else:
				Util.mul_img(vol, mask)
		volft, kbx, kby, kbz = prep_vol_new(volume)
		sizez=volume.get_zsize()


	if(type(stack) is types.StringType):
		Disk = True
		os.system("rm -f  "+stack)	
	else:
		out = []
		Disk = False
	
	s2x=0
	s2y=0
	for i in xrange(len(angles)):
	#for i in xrange(1):
		#print "i==%d\t"%(i)
		if(len(angles[i]) == 3):
			proj = prgs_new_fast(volft, sizez,kbx, kby, kbz, [angles[i][0], angles[i][1], angles[i][2], 0.0, 0.0])
			set_params_proj(proj, [angles[i][0], angles[i][1], angles[i][2], 0.0, 0.0])
		else:
			proj = prgs_new_fast(volft, sizez,kbx, kby, kbz, [angles[i][0], angles[i][1], angles[i][2], -angles[i][3], -angles[i][4]])
			set_params_proj(proj, angles[i])

		# horatio active_refactoring Jy51i1EwmLD4tWZ9_00000_1
		# proj.set_attr_dict({'active':1})

		# add noise, if noise is set. this is two-fold: application of noise before
		#    ctf filtering and after it.
		if noise is not None:
			try:
				# no mask, so call w/ false
				noise_ima = model_gauss_noise(noise_level,proj.get_xsize(),
							      proj.get_ysize())
			except:
				pass
			else:
				proj += noise_ima

		# apply ctf, if ctf option is set and if we can create a valid CTF object
		try:
			ctf = EMAN2Ctf()
			# order of values is the one applied in sxheader for import / export!
			ctf.from_dict({ "defocus":ctfs[i][0], "cs":ctfs[i][1], "voltage":ctfs[i][2], 
					"apix":ctfs[i][3], "bfactor":ctfs[i][4], "ampcont":ctfs[i][5] })
		except:
			# there are no ctf values, so ignore this and set no values
			proj.set_attr( "error",1)
		else:
			# setting of values worked, so apply ctf and set the header info correctly
			proj = filt_ctf(proj,ctf)
			proj.set_attr( "ctf",ctf)
			proj.set_attr( "ctf_applied",0)

		# add second noise level that is not affected by CTF
		if noise is not None:
			try:
				noise_ima = model_gauss_noise(noise_level,proj.get_xsize(),
							      proj.get_ysize())
			except:
				pass
			else:
				proj += noise_ima

		if(Disk):
			proj.write_image(stack, i)
		else: 
			out.append(proj)
	if(not Disk):  return out
	
	
def prep_vol_new(vol):
	"""
		Name
			prep_vol - prepare the volume for calculation of gridding projections and generate the interpolants.
		Input
			vol: input volume for which projections will be calculated using prgs
			the input is rectaungalar real space image
		Output
			volft: volume prepared for gridding projections using prgs
			kb: interpolants (tabulated Kaiser-Bessel function)
	"""
	
	# prepare the volume
	
	Mx     = vol.get_xsize()
	My     = vol.get_ysize()
	Mz     = vol.get_zsize()
	# padd two times
	npad  = 2
	Nx     = Mx*npad
	Ny     = My*npad
	Nz     = Mz*npad
	# support of the window
	K     = 6
	alpha = 1.75
	rx     = Mx/2
	ry     = My/2
	rz     = Mz/2
	#v     = K/2.0/N
	kbx    = Util.KaiserBessel(alpha, K, rx, K/(2.*Nx), Nx)
	kby    = Util.KaiserBessel(alpha, K, ry, K/(2.*Ny), Ny)
	kbz    = Util.KaiserBessel(alpha, K, rz, K/(2.*Nz), Nz)
	volft = vol.copy()
	volft.divkbsinh_rect(kbx,kby,kbz)
	volft = volft.norm_pad(False, npad)
	volft.do_fft_inplace()
	volft.center_origin_fft()
	volft.fft_shuffle()
	return  volft,kbx,kby,kbz

def prgs_new(volft, kbx, kby, kbz, params):
	"""
		Name
			prg - calculate 2-D projection of a 3-D volume
		Input
			vol: input volume, all dimensions have to be the same (nx=ny=nz)
			params: input parameters given as a list [phi, theta, psi, s2x, s2y], projection in calculated using the three Eulerian angles and then shifted by sx,sy
		Output
			proj: generated 2-D projection
	"""
	#  params:  phi, theta, psi, sx, sy
	from fundamentals import fft
	from utilities import set_params_proj
	R = Transform({"type":"spider", "phi":params[0], "theta":params[1], "psi":params[2]})
	temp = volft.extract_plane_rect(R,kbx,kby,kbz)
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
	
def prgs_new_fast(volft,sizez,kbx, kby, kbz, params):
	"""
		Name
			prg - calculate 2-D projection of a 3-D volume
		Input
			vol: input volume, all dimensions have to be the same (nx=ny=nz)
			params: input parameters given as a list [phi, theta, psi, s2x, s2y], projection in calculated using the three Eulerian angles and then shifted by sx,sy
		Output
			proj: generated 2-D projection
	"""
	#  params:  phi, theta, psi, sx, sy
	#print "fast projection is called"
	from fundamentals import fft
	from utilities import set_params_proj
	R = Transform({"type":"spider", "phi":params[0], "theta":params[1], "psi":params[2]})
	temp = volft.extract_plane_rect_fast(R,kbx,kby,kbz)
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
	temp=Util.pad(temp,sizez,sizez,1)
	#print (temp.get_xsize(),temp.get_ysize(),temp.get_zsize(),sizez)
	return temp
def ihrsr_new(stack, ref_vol, outdir, maskfile, ir, ou, rs, xr, ynumber, 
          txs, delta, an, maxit, CTF, snr, dp,ndp,dp_step, dphi,ndphi,dphi_step, psi_max,
	  rmin, rmax, fract, nise, npad, sym, user_func_name, datasym,
	  fourvar, debug = False, MPI = False):
	if MPI:
		ihrsr_MPI_new(stack, ref_vol, outdir, maskfile, ir, ou, rs, xr, ynumber, 
			txs, delta, an, maxit, CTF, snr, dp,ndp,dp_step, dphi, ndphi,dphi_step,psi_max,
			rmin, rmax, fract, nise, npad, sym, user_func_name, datasym,
			fourvar, debug)
		return

	from utilities      import model_circle, drop_image
	from utilities      import get_image, get_input_from_string
	from utilities      import get_params_proj, set_params_proj
	from alignment	    import proj_ali_helical, helios, Numrinit, prepare_refrings
	from projection     import prep_vol
	from statistics     import ccc
	from fundamentals   import cyclic_shift, rot_shift3D
	import os
	import types
	from utilities      import print_begin_msg, print_end_msg, print_msg
	print_begin_msg("ihrsr")

	import user_functions
	user_func = user_functions.factory[user_func_name]

	if os.path.exists(outdir):  ERROR('Output directory exists, please change the name and restart the program', "ihrsr", 1)
	os.mkdir(outdir)
        
	xrng        = get_input_from_string(xr)
	#Guozhi Tao changed----since I cannot test this code, just modify the code to make sure no crash--begin
	if  ny == "-1":  yrng = xrng
	else          :  yrng = get_input_from_string(xr)
	#since I cannot test this code, just modify the code to make sure no crash--end
	step        = get_input_from_string(txs)
	delta       = get_input_from_string(delta)
	lstp = min( len(xrng), len(yrng), len(step), len(delta) )
	if (an == "-1"):
		an = [-1] * lstp
	else:
		an = get_input_from_string(an)
	first_ring  = int(ir)
	rstep       = int(rs)
	last_ring   = int(ou)
	max_iter    = int(maxit)

	print_msg("Input stack                               : %s\n"%(stack))
	print_msg("Reference volume                          : %s\n"%(ref_vol))	
	print_msg("Output directory                          : %s\n"%(outdir))
	print_msg("Maskfile                                  : %s\n"%(maskfile))
	print_msg("Inner radius                              : %i\n"%(first_ring))

	vol     = EMData()
	vol.read_image(ref_vol)
	nx      = vol.get_xsize()
	if (last_ring == -1):	last_ring = nx//2 - 2

	print_msg("Outer radius                              : %i\n"%(last_ring))
	print_msg("Ring step                                 : %i\n"%(rstep))
	print_msg("X search range                            : %s\n"%(xrng))
	print_msg("Y search range                            : %s\n"%(yrng))
	print_msg("Translational step                        : %s\n"%(step))
	print_msg("Angular step                              : %s\n"%(delta))
	print_msg("Angular search range                      : %s\n"%(an))
	print_msg("max radius for helical search (in Ang)    : %f\n"%(rmax))
	print_msg("fraction of volume used for helical search: %f\n"%(fract))
	print_msg("initial symmetry - angle                  : %f\n"%(dphi))
	print_msg("initial symmetry - axial rise             : %f\n"%(dp))
	print_msg("Maximum number of iterations              : %i\n"%(max_iter))
	print_msg("CTF correction                            : %s\n"%(CTF))
	print_msg("Signal-to-Noise Ratio                     : %f\n"%(snr))
	print_msg("Symmetry group                            : %s\n"%(sym))
	print_msg("symmetry doc file                         : %s\n"%(datasym))
	print_msg("npad                                      : %i\n"%(npad))
	print_msg("User function                             : %s\n"%(user_func_name))

	if (maskfile) :
		if (type(maskfile) is types.StringType): mask3D = get_image(maskfile)
		else                                  : mask3D = maskfile
	else          :   mask3D = model_circle(last_ring, nx, nx, nx)
	mask2D = model_circle(last_ring, nx, nx) - model_circle(first_ring, nx, nx)
	numr   = Numrinit(first_ring, last_ring, rstep, "F")

	if CTF:
		from reconstruction import recons3d_4nn_ctf
		from filter import filt_ctf

	else: from reconstruction import recons3d_4nn


	#drop_image(vol, os.path.join(outdir,"ref_vol00.hdf"))
	sym = "c1"
	symref = "s"
	ref_a= "P"

	# horatio active_refactoring Jy51i1EwmLD4tWZ9_00000_1
	# active = EMUtil.get_all_attributes(stack, 'active')
	# list_of_particles = []
	# for im in xrange(len(active)):
	# 	if(active[im]):  list_of_particles.append(im)
	# del active
	
	nima = EMUtil.get_image_count(stack)
	list_of_particles = range(nima)

	data = EMData.read_images(stack, list_of_particles)
        for im in xrange(len(data)):
                data[im].set_attr('ID', list_of_particles[im])
	nima = len(data)
	if(data[0].get_attr_default('ctf_applied', 2) > 0):  ctf_applied = True
	else:   ctf_applied = False
	if CTF:
		pixel_size = data[0].get_attr('ctf').apix
	else:
		pixel_size = data[0].get_attr('pixel_size')
	print_msg("Pixel size in Angstroms                   : %f\n\n"%(pixel_size))

	finfo = None#open("desperado", 'w')
	# do the projection matching
	drop_image(vol, os.path.join(outdir, "aligned0000.hdf"))
	for N_step in xrange(lstp):
		for Iter in xrange(max_iter):
			print_msg("ITERATION #%3d\n"%(N_step*max_iter + Iter+1))

			if CTF:
				if ctf_applied == False:
					previous_defocus = -1.0
				else:
					volft,kb = prep_vol( vol )
					refrings = prepare_refrings( volft, kb, nx, delta[N_step], ref_a, symref, numr, MPI=False, phiEqpsi = "Zero")
					del volft, kb
			else:
				volft,kb = prep_vol( vol )
				refrings = prepare_refrings( volft, kb, nx, delta[N_step], ref_a, symref, numr, MPI=False, phiEqpsi = "Zero")
				del volft, kb
			sx = 0.0
			sy = 0.0
			phihi = [0.0]*nima
			for im in xrange( nima ):
				if CTF and ctf_applied == False:
					ctf = data[im].get_attr( "ctf" )
					if ctf.defocus != previous_defocus:
						previous_defocus = ctf.defocus
						ctfvol = filt_ctf(vol, ctf)
						volft,kb = prep_vol( ctfvol )
						refrings = prepare_refrings( volft, kb, nx, delta[N_step], ref_a, symref, numr, MPI=False, phiEqpsi = "Zero")
						del volft, kb

				#if an[N_step] == -1:
				#	peak, pixel_error = proj_ali_helical(data[im],refrings,numr,xrng[N_step],yrng[N_step],step[N_step])
				#else:
				peak, phihi[im], sxi, syi, pixel_error = proj_ali_helical(data[im],refrings,numr,xrng[N_step],yrng[N_step],step[N_step], psi_max)

				sx += sxi
				sy += syi
				#print_msg("Image %i,  psi %9.2f,    s2x  %9.2f, s2y  %9.2f,  peak  %10.3e \n"%(im, paramali[2], paramali[3], paramali[4], peak))

			# histogram of phi's
			from statistics import hist_list
			lhist = 30
			region, histo = hist_list(phihi, lhist)
			msg = "      Histogram of phi angles\n      phi         number of particles\n"
			print_msg(msg)
			for lhx in xrange(lhist):
				msg = " %10.3f     %7d\n"%(region[lhx], histo[lhx])
				print_msg(msg)
			del phihi, region, histo, msg
			#center projections
			sx /= nima
			sy /= nima
			for im in xrange( nima ):
				paramali = get_params_proj(data[im])
				set_params_proj(data[im], [paramali[0], paramali[1], paramali[2], paramali[3] - sx, paramali[4] - sy])

			#  3D stuff
			#  calculate new and improved 3D
			if(CTF): vol = recons3d_4nn_ctf(data, range(nima), snr, npad = npad)
			else:	 vol = recons3d_4nn(data, range(nima), npad = npad)

			ref_data = [vol]
			vol = user_func(ref_data)

			# store the reference volume
			drop_image(vol, os.path.join(outdir, "unsymmetrized%04d.hdf"%(N_step*max_iter+Iter+1)))
			if(N_step*max_iter+Iter+1 > 2):
				vol, dp, dphi = helios(vol, pixel_size, dp, dphi, fract, rmax, rmin)
			else:
				#  in the first two steps the symmetry is imposed
				vol = vol.helicise(pixel_size,dp, dphi, fract, rmax, rmin)
			print_msg("New delta z and delta phi      : %s,    %s\n\n"%(dp,dphi))
			fofo = open(os.path.join(outdir,datasym),'a')
			fofo.write('  %12.4f   %12.4f\n'%(dp,dphi))
			fofo.close()
			drop_image(vol, os.path.join(outdir, "aligned%04d.hdf"%(N_step*max_iter+Iter+1)) )
			#  here we  write header info
			from utilities import write_headers
			write_headers( stack, data, list_of_particles)
	print_end_msg("ihrsr_new")

def ihrsr_MPI_new(stack, ref_vol, outdir, maskfile, ir, ou, rs, xr, ynumber, 
	txs, delta, an, maxit, CTF, snr, dp,ndp,sndp, dphi,ndphi,sndphi, psi_max,
	rmin, rmax, fract, nise, npad, sym, user_func_name, datasym,
	fourvar, debug):

	from alignment      import Numrinit, proj_ali_helical 
	from utilities      import model_circle, get_image, drop_image, get_input_from_string
	from utilities      import bcast_list_to_all, bcast_number_to_all, reduce_EMData_to_root, bcast_EMData_to_all
	from utilities      import send_attr_dict
	from utilities      import get_params_proj, set_params_proj, file_type
	from fundamentals   import rot_avg_image
	from pixel_error    import max_3D_pixel_error
	import os
	import types
	from utilities      import print_begin_msg, print_end_msg, print_msg
	from mpi            import mpi_bcast, mpi_comm_size, mpi_comm_rank, MPI_FLOAT, MPI_COMM_WORLD, mpi_barrier
	from mpi            import mpi_recv,  mpi_send, MPI_TAG_UB
	from mpi            import mpi_reduce, MPI_INT, MPI_SUM
	from filter         import filt_ctf
	#from projection     import prep_vol, prgs
	from development    import prep_vol_new,prgs_new,prepare_refrings_new,helios7_new
	from statistics     import hist_list, varf3d_MPI
	from applications   import MPI_start_end


	number_of_proc = mpi_comm_size(MPI_COMM_WORLD)
	myid           = mpi_comm_rank(MPI_COMM_WORLD)
	main_node = 0
	
	if myid == 0:
		if os.path.exists(outdir):  nx = 1
		else:  nx = 0
	else:  nx = 0
	ny = bcast_number_to_all(nx, source_node = main_node)
	
	if ny == 1:  ERROR('Output directory exists, please change the name and restart the program', "ihrsr_MPI", 1,myid)
	mpi_barrier(MPI_COMM_WORLD)

	if myid == main_node:
		os.mkdir(outdir)
	mpi_barrier(MPI_COMM_WORLD)
	
	#ndp    = 12
	#sndp   = 0.1
	#ndphi  = 12
	#sndphi = 0.1
	nlprms = (2*ndp+1)*(2*ndphi+1)
	if nlprms< number_of_proc:
		ERROR('number of cpus is larger than the number of helical search, please reduce it or at this moment modify ndp,dphi in the program', "ihrsr_MPI", 1,myid)
	


	if debug:
		from time import sleep
		while not os.path.exists(outdir):
			print  "Node ",myid,"  waiting..."
			sleep(5)

		info_file = os.path.join(outdir, "progress%04d"%myid)
		finfo = open(info_file, 'w')
	else:
		finfo = None

	sym = "c1"
	symref = "s"
	ref_a= "P"

	xrng        = get_input_from_string(xr)
	ynumber	    = get_input_from_string(ynumber)
	for i in xrange(len(ynumber)):
		if(ynumber[i]%2==1):
			ynumber[i]=ynumber[i]+1
	yrng =[]
	
	for i in xrange(len(xrng)):
		yrng.append(dp/2)
	
	stepx        = get_input_from_string(txs)
	delta       = get_input_from_string(delta)
	lstp = min(len(xrng), len(yrng), len(stepx), len(delta))
	if an == "-1":
		an = [-1] * lstp
	else:
		an = get_input_from_string(an)

	first_ring  = int(ir)
	rstep       = int(rs)
	last_ring   = int(ou)
	max_iter    = int(maxit)

	vol     = EMData()
	vol.read_image(ref_vol)
	nx      = vol.get_xsize()
	ny      = vol.get_ysize()
	nz      = vol.get_zsize()
	if last_ring < 0:	last_ring = int(nx/2) - 2

	if myid == main_node:
		import user_functions
		user_func = user_functions.factory[user_func_name]

		print_msg("Input stack                               : %s\n"%(stack))
		print_msg("Reference volume                          : %s\n"%(ref_vol))	
		print_msg("Output directory                          : %s\n"%(outdir))
		print_msg("Maskfile                                  : %s\n"%(maskfile))
		print_msg("Inner radius                              : %i\n"%(first_ring))
		print_msg("Outer radius                              : %i\n"%(last_ring))
		print_msg("Ring step                                 : %i\n"%(rstep))
		print_msg("X search range                            : %s\n"%(xrng))
		print_msg("Y number                                  : %s\n"%(ynumber))
		print_msg("Translational stepx                        : %s\n"%(stepx))
		print_msg("Angular step                              : %s\n"%(delta))
		print_msg("Angular search range                      : %s\n"%(an))
		print_msg("min radius for helical search (in pix)    : %f\n"%(rmin))
		print_msg("max radius for helical search (in pix)    : %f\n"%(rmax))
		print_msg("fraction of volume used for helical search: %f\n"%(fract))
		print_msg("initial symmetry - angle                  : %f\n"%(dphi))
		print_msg("initial symmetry - axial rise             : %f\n"%(dp))
		print_msg("Maximum iteration                         : %i\n"%(max_iter))
		print_msg("Data with CTF                             : %s\n"%(CTF))
		print_msg("Signal-to-Noise Ratio                     : %f\n"%(snr))
		print_msg("symmetry output doc file                  : %s\n"%(datasym))
		print_msg("number of times initial symmetry is imposed: %i\n"%(nise))
		print_msg("npad                                      : %i\n"%(npad))
		print_msg("User function                             : %s\n"%(user_func_name))

	if maskfile:
		if type(maskfile) is types.StringType: mask3D = get_image(maskfile)
		else:                                  mask3D = maskfile
	#else: mask3D = model_circle(last_ring, nx, nx, nx)

	numr	= Numrinit(first_ring, last_ring, rstep, "F")
	mask2D  = model_circle(last_ring,nz,nz) - model_circle(first_ring,nz,nz)

	if CTF:
		from development import recons3d_4nn_ctf_MPI_new
		from filter         import filt_ctf
	else:	 from development import recons3d_4nn_MPI_new

	if myid == main_node:
		if(file_type(stack) == "bdb"):
			from EMAN2db import db_open_dict
			dummy = db_open_dict(stack, True)
		
		# horatio active_refactoring Jy51i1EwmLD4tWZ9_00000_1		
		# active = EMUtil.get_all_attributes(stack, 'active')
		# list_of_particles = []
		# for im in xrange(len(active)):
		# 	if active[im]:  list_of_particles.append(im)
		# del active
		# nima = len(list_of_particles)

		nima = EMUtil.get_image_count(stack)
		list_of_particles = range(nima)
	else:
		nima = 0
	total_nima = bcast_number_to_all(nima, source_node = main_node)

	if myid != main_node:
		list_of_particles = [-1]*total_nima
	list_of_particles = bcast_list_to_all(list_of_particles, source_node = main_node)

	image_start, image_end = MPI_start_end(total_nima, number_of_proc, myid)
	# create a list of images for each node
	list_of_particles = list_of_particles[image_start: image_end]
	nima = len(list_of_particles)
	if debug:
		finfo.write("image_start, image_end: %d %d\n" %(image_start, image_end))
		finfo.flush()

	data = EMData.read_images(stack, list_of_particles)
	if fourvar:  original_data = []
	for im in xrange(nima):
		data[im].set_attr('ID', list_of_particles[im])
		if fourvar: original_data.append(data[im].copy())
		if CTF:
			st = Util.infomask(data[im], mask2D, False)
			data[im] -= st[0]
			st = data[im].get_attr_default("ctf_applied", 0)
			if(st == 0):
				ctf_params = data[im].get_attr("ctf")
				data[im] = filt_ctf(data[im], ctf_params)
				data[im].set_attr('ctf_applied', 1)

	if debug:
		finfo.write( '%d loaded  \n' % nima )
		finfo.flush()
	if CTF:
		pixel_size = data[0].get_attr('ctf').apix
	else:
		pixel_size = data[0].get_attr('pixel_size')
	for i in xrange(len(xrng)):
		yrng[i]=dp/(2*pixel_size)
		
	if myid == main_node:
		print_msg("Pixel size in Angstroms                   : %f\n\n"%(pixel_size))
		print_msg("Y search range (pix) initialized as       : %s\n\n"%(yrng))
		

	from time import time	

	#  this is needed for gathering of pixel errors
	disps = []
	recvcount = []
	for im in xrange(number_of_proc):
		if( im == main_node ):  disps.append(0)
		else:                  disps.append(disps[im-1] + recvcount[im-1])
		ib, ie = MPI_start_end(total_nima, number_of_proc, im)
		recvcount.append( ie - ib )

	#jeanmod
	total_iter = 0
	# do the projection matching
	for N_step in xrange(lstp):
		terminate = 0
		Iter = 0
 		while(Iter < max_iter-1 and terminate == 0):
			yrng[N_step]=float(dp)/(2*pixel_size) #will change it later according to dp
			if(ynumber[N_step]==0):
				stepy=0.0
			else:
				stepy=(2*yrng[N_step]/ynumber[N_step])
				
			pixer  = [0.0]*nima
			modphi = [0.0]*nima
			Iter += 1
			total_iter += 1
			if myid == main_node:
				start_time = time()
				print_msg("\nITERATION #%3d,  inner iteration #%3d\nDelta = %4.1f, an = %5.2f, xrange = %5.2f,stepx = %5.2f, yrange = %5.2f,  stepy = %5.2f, ynumber = %3d\n"%(total_iter, Iter, delta[N_step], an[N_step], xrng[N_step],stepx[N_step],yrng[N_step],stepy,ynumber[N_step]))

			volft,kbx,kby,kbz = prep_vol_new( vol )
			refrings = prepare_refrings_new( volft, kbx,kby,kbz, nz, delta[N_step], ref_a, symref, numr, MPI = True, phiEqpsi = "Zero")
			del volft,kbx,kby,kbz
			if myid== main_node:
				print_msg( "Time to prepare rings: %d\n" % (time()-start_time) )
				start_time = time()

			for im in xrange( nima ):

				peak, phihi, theta, psi, sxi, syi, t1 = proj_ali_helical(data[im],refrings,numr,xrng[N_step],yrng[N_step],stepx[N_step],ynumber[N_step],psi_max,finfo,)
				if(peak > -1.0e22):
					#Jeanmod: wrap y-shifts back into box within rise of one helical unit by changing phi
					dpp = (float(dp)/pixel_size)
					dyi = (float(syi)/dpp)-int(syi/dpp)
					#jdelta = 0.75/dpp # jdelta could be adjustable input parameter
					jdelta=0.0
						
					if(abs(syi)<=(0.5+jdelta)*dpp):
						synew  = syi
						phinew = phihi
					else:
						if dyi < -0.5-jdelta:  eyi = dyi+1.0
						elif dyi > 0.5+jdelta:  eyi = dyi-1.0
						else:                   eyi = dyi
	       					synew  = eyi*dpp
	        				phinew = phihi+dphi*float(synew-syi)/dpp
						phinew = phinew%360

					t2 = Transform({"type":"spider","phi":phinew,"theta":theta,"psi":psi})
					t2.set_trans(Vec2f(-sxi, -synew))
					data[im].set_attr("xform.projection", t2)
					pixer[im]  = max_3D_pixel_error(t1, t2, numr[-3])
					modphi[im] = phinew
				else:
					# peak not found, parameters not modified
					pixer[im]  = 0.0
					phihi, theta, psi, sxi, syi = get_params_proj(data[im])
					modphi[im] = phihi

			if myid == main_node:
				print_msg("Time of alignment = %d\n"%(time()-start_time))
				start_time = time()

			#output pixel errors
			from mpi import mpi_gatherv
			recvbuf = mpi_gatherv(pixer, nima, MPI_FLOAT, recvcount, disps, MPI_FLOAT, main_node, MPI_COMM_WORLD)
			mpi_barrier(MPI_COMM_WORLD)
			terminate = 0
			if(myid == main_node):
				recvbuf = map(float, recvbuf)
				from utilities import write_text_file
				write_text_file([range(len(recvbuf)), recvbuf], os.path.join(outdir, "pixer_%04d_%04d.txt"%(N_step+1,Iter)) )
				from statistics import hist_list
				lhist = 20
				region, histo = hist_list(recvbuf, lhist)
				if(region[0] < 0.0):  region[0] = 0.0
				msg = "      Histogram of pixel errors\n      ERROR       number of particles\n"
				print_msg(msg)
				for lhx in xrange(lhist):
					msg = " %10.3f     %7d\n"%(region[lhx], histo[lhx])
					print_msg(msg)
				# Terminate if 95% within 1 pixel error
				im = 0
				for lhx in xrange(lhist):
					if(region[lhx] > 1.0): break
					im += histo[lhx]
				#if(im/float(total_nima) > 0.95):  terminate = 1
				del region, histo
			#output distribution of phi
			#jeanmod
			recvbuf = mpi_gatherv(modphi, nima, MPI_FLOAT, recvcount, disps, MPI_FLOAT, main_node, MPI_COMM_WORLD)
			#end jeanmod
			mpi_barrier(MPI_COMM_WORLD)
			if(myid == main_node):
				recvbuf = map(float, recvbuf)
				lhist = 30
				region, histo = hist_list(recvbuf, lhist)
				msg = "\n      Distribution of phi\n      phi         number of particles\n"
				print_msg(msg)
				for lhx in xrange(lhist):
					msg = " %10.3f     %7d\n"%(region[lhx], histo[lhx])
					print_msg(msg)
				del region, histo
			del recvbuf
			terminate = mpi_bcast(terminate, 1, MPI_INT, 0, MPI_COMM_WORLD)
			terminate = int(terminate[0])
			del pixer, modphi
			if myid == main_node:
				print_msg("Time to compute pixer = %d\n"%(time()-start_time))
				start_time = time()
			# write out headers, under MPI writing has to be done sequentially
			mpi_barrier(MPI_COMM_WORLD)
			m = 5
			from mpi import mpi_recv, mpi_send, MPI_TAG_UB, MPI_COMM_WORLD, MPI_FLOAT
			if myid == main_node:
				
				fexp = open(os.path.join(outdir, "parameters_%04d_%04d.txt"%(N_step+1,Iter)),"w")
				for n in xrange(number_of_proc):
					if n!=main_node:
						t = mpi_recv(recvcount[n]*m,MPI_FLOAT, n, MPI_TAG_UB, MPI_COMM_WORLD)
						for i in xrange(recvcount[n]):
							for j in xrange(m):
								fexp.write(" %15.5f  "%t[j+i*m])
							fexp.write("\n")
					else:
						t = [0.0]*m
						for i in xrange(recvcount[myid]):
							t = get_params_proj(data[i])
							for j in xrange(m):
								fexp.write(" %15.5f  "%t[j])
							fexp.write("\n")
				fexp.close()
				del t
	        	else:
				nvalue = [0.0]*m*recvcount[myid]
				t = [0.0]*m
				for i in xrange(recvcount[myid]):
					t = get_params_proj(data[i])
					for j in xrange(m):
						nvalue[j + i*m] = t[j]
				mpi_send(nvalue, recvcount[myid]*m, MPI_FLOAT, main_node, MPI_TAG_UB, MPI_COMM_WORLD)
				del nvalue
			if myid == main_node:
				print_msg("Time to write parameters = %d\n"%(time()-start_time))
				start_time = time()
			
			rx=float(nx)/nz
			ry=float(ny)/nz

			if CTF: vol = recons3d_4nn_ctf_MPI_new(myid, data, snr = snr, npad = npad, rx=rx, ry=ry)
			else:    vol = recons3d_4nn_MPI_new(myid, data, npad = npad, rx=rx, ry=ry)

			if myid == main_node:
				print_msg("\n3D reconstruction time = %d\n"%(time()-start_time))
				start_time = time()

			if fourvar:
			#  Compute Fourier variance
				for im in xrange(nima):
					original_data[im].set_attr( 'xform.projection', data[im].get_attr('xform.projection') )
				varf = varf3d_MPI(original_data, ssnr_text_file = os.path.join(outdir, "ssnr%04d"%(total_iter)), mask2D = None, reference_structure = vol, ou = last_ring, rw = 1.0, npad = 1, CTF = CTF, sign = 1, sym =sym, myid = myid)
				if myid == main_node:
					print_msg("Time to calculate 3D Fourier variance= %d\n"%(time()-start_time))
					start_time = time()
					varf = 1.0/varf
					drop_image(varf, os.path.join(outdir, "varf%04d.hdf"%(total_iter)))
			else:  varf = None

			#search for helical symmetry
			if myid == main_node:
				drop_image(vol, os.path.join(outdir, "vol%04d.hdf"%(total_iter)))
							
			if(total_iter > nise):
				bcast_EMData_to_all(vol, myid, main_node)
				#from filter import filt_gaussl
				#vol = filt_gaussl(vol, 0.25)

				if myid == main_node:
					lprms = []
					for i in xrange(-ndp,ndp+1,1):
						for j in xrange(-ndphi,ndphi+1,1):
							lprms.append( dp   + i*sndp)
							lprms.append( dphi + j*sndphi)
					#print "lprms===",lprms
					recvpara = []
					for im in xrange(number_of_proc):
						helic_ib, helic_ie = MPI_start_end(nlprms, number_of_proc, im)
						recvpara.append(helic_ib )
						recvpara.append(helic_ie )

				para_start, para_end = MPI_start_end(nlprms, number_of_proc, myid)

				list_dps     = [0.0]*((para_end-para_start)*2)
				list_fvalues = [-1.0]*((para_end-para_start)*1)

				if myid == main_node:
					for n in xrange(number_of_proc):
						if n!=main_node: mpi_send(lprms[2*recvpara[2*n]:2*recvpara[2*n+1]], 2*(recvpara[2*n+1]-recvpara[2*n]), MPI_FLOAT, n, MPI_TAG_UB, MPI_COMM_WORLD)
						else:    list_dps = lprms[2*recvpara[2*0]:2*recvpara[2*0+1]]
				else:
					list_dps = mpi_recv((para_end-para_start)*2, MPI_FLOAT, main_node, MPI_TAG_UB, MPI_COMM_WORLD)

				list_dps = map(float, list_dps)

				local_pos = [0.0, 0.0, -1.0e20]
				for i in xrange(para_end-para_start):
					fvalue = helios7_new(vol, pixel_size, list_dps[i*2], list_dps[i*2+1], fract, rmax, rmin)
					if(fvalue >= local_pos[2]):
						local_pos = [list_dps[i*2], list_dps[i*2+1], fvalue ]
				if myid == main_node:
					list_return = [0.0]*(3*number_of_proc)
					for n in xrange(number_of_proc):
						if n != main_node: list_return[3*n:3*n+3]                 = mpi_recv(3,MPI_FLOAT, n, MPI_TAG_UB, MPI_COMM_WORLD)
 						else:              list_return[3*main_node:3*main_node+3]  = local_pos[:]
				else:
					mpi_send(local_pos, 3, MPI_FLOAT, main_node, MPI_TAG_UB, MPI_COMM_WORLD)

				if myid == main_node:	
					maxvalue = list_return[2]
					for i in xrange(number_of_proc):
						if( list_return[i*3+2] >= maxvalue ):
							maxvalue = list_return[i*3+2]
							dp       = list_return[i*3+0]
							dphi     = list_return[i*3+1]
					dp   = float(dp)
					dphi = float(dphi)
					#print  "  GOT dp dphi",dp,dphi

					vol  = vol.helicise_rect(pixel_size,dp, dphi, fract, rmax, rmin)
					print_msg("New delta z and delta phi      : %s,    %s\n\n"%(dp,dphi))		
			
			else:
				if myid==main_node:
					#  in the first nise steps the symmetry is imposed
					vol = vol.helicise_rect(pixel_size,dp, dphi, fract, rmax, rmin)
					print_msg("Imposed delta z and delta phi      : %s,    %s\n\n"%(dp,dphi))
			if(myid==main_node):
				fofo = open(os.path.join(outdir,datasym),'a')
				fofo.write('  %12.4f   %12.4f\n'%(dp,dphi))
				fofo.close()
				ref_data = [vol]
				if  fourvar:  ref_data.append(varf)
				vol = user_func(ref_data)

				drop_image(vol, os.path.join(outdir, "volf%04d.hdf"%(total_iter)))
				print_msg("\nSymmetry search and user function time = %d\n"%(time()-start_time))
				start_time = time()

			bcast_EMData_to_all(vol, myid, main_node)
			dp   = bcast_number_to_all(dp,   source_node = main_node)
			dphi = bcast_number_to_all(dphi, source_node = main_node)
			#
			del varf
	par_str = ["xform.projection"]
	if myid == main_node:
	   	if(file_type(stack) == "bdb"):
	        	from utilities import recv_attr_dict_bdb
	        	recv_attr_dict_bdb(main_node, stack, data, par_str, image_start, image_end, number_of_proc)
	        else:
	        	from utilities import recv_attr_dict
	        	recv_attr_dict(main_node, stack, data, par_str, image_start, image_end, number_of_proc)
		print_msg("Time to write header information= %d\n"%(time()-start_time))
		start_time = time()
	else:	       send_attr_dict(main_node, data, par_str, image_start, image_end)
	if myid == main_node: print_end_msg("ihrsr_MPI")
	
	
def prepare_refrings_new( volft, kbx,kby,kbz, nz, delta, ref_a, sym, numr, MPI=False, phiEqpsi = "Minus"):
        #from projection   import prep_vol, prgs
	from development   import prgs_new
        from math         import sin, cos, pi
	from applications import MPI_start_end
	from utilities    import even_angles
	from alignment	  import ringwe
	# generate list of Eulerian angles for reference projections
	#  phi, theta, psi
	mode = "F"
	ref_angles = even_angles(delta, symmetry=sym, method = ref_a, phiEqpsi = phiEqpsi)
	wr_four  = ringwe(numr, mode)
	cnx = nz//2 + 1
	cny = nz//2 + 1
	qv = pi/180.
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

	refrings = []     # list of (image objects) reference projections in Fourier representation

	sizex = numr[ len(numr)-2 ] + numr[ len(numr)-1 ] - 1

        for i in xrange(num_ref):
		prjref = EMData()
		prjref.set_size(sizex, 1, 1)
		refrings.append(prjref)

        for i in xrange(ref_start, ref_end):
		prjref = prgs_new(volft, kbx,kby,kbz, [ref_angles[i][0], ref_angles[i][1], ref_angles[i][2], 0.0, 0.0])
		cimage = Util.Polar2Dm(prjref, cnx, cny, numr, mode)  # currently set to quadratic....
		Util.Normalize_ring(cimage, numr)

		Util.Frngs(cimage, numr)
		Util.Applyws(cimage, numr, wr_four)
		refrings[i] = cimage

	if MPI:
		from utilities import bcast_EMData_to_all
		for i in xrange(num_ref):
			for j in xrange(ncpu):
				ref_start,ref_end = MPI_start_end(num_ref,ncpu,j)
				if i >= ref_start and i < ref_end: rootid = j

			bcast_EMData_to_all( refrings[i], myid, rootid )

	for i in xrange(len(ref_angles)):
		n1 = sin(ref_angles[i][1]*qv)*cos(ref_angles[i][0]*qv)
		n2 = sin(ref_angles[i][1]*qv)*sin(ref_angles[i][0]*qv)
		n3 = cos(ref_angles[i][1]*qv)
		refrings[i].set_attr_dict( {"n1":n1, "n2":n2, "n3":n3} )
		refrings[i].set_attr("phi", ref_angles[i][0])
		refrings[i].set_attr("theta", ref_angles[i][1])
		refrings[i].set_attr("psi", ref_angles[i][2])

	return refrings

def helios_func_new(params, data):
	sm = data[0].helicise_rect(data[2], params[0], params[1], data[3], data[4], data[5])
	#try other sim creteria
	q = sm.cmp("dot", sm, {"negative":0})
	#q = sm.cmp("dot", data[0], {"negative":0})# corelation  with the recon data
	#print  params,q
	return  q	
	
def helios7_new(vol, pixel_size, dp, dphi, section_use = 0.75, radius = 0.0, rmin = 0.0):
	from development    import helios_func_new
	nx = vol.get_xsize()
	ny = vol.get_ysize()
	nz = vol.get_zsize()
	if(radius <= 0.0):    radius = nx//2-1
	params = [dp, dphi]
	data=[vol, params, pixel_size, section_use, radius, rmin]
	q = helios_func_new([dp, dphi], data)
	return q



def mrefeq_ali2d(stack, refim, outdir, maskfile=None, ir=1, ou=-1, rs=1, xrng=0, yrng=0, step=1, center=1, maxit=0, CTF=False, snr=1.0, user_func_name="ref_ali2d", rand_seed=1000, MPI=False):
	"""
		Name
			mref_ali2d - Perform 2-D multi-reference alignment of an image series
		Input
			stack: set of 2-D images in a stack file, images have to be squares
			refim: set of initial reference 2-D images in a stack file 
			maskfile: optional maskfile to be used in the alignment
			inner_radius: inner radius for rotational correlation > 0
			outer_radius: outer radius for rotational correlation < nx/2-1
			ring_step: step between rings in rotational correlation >0
			x_range: range for translation search in x direction, search is +/xr 
			y_range: range for translation search in y direction, search is +/yr 
			translation_step: step of translation search in both directions
			center: center the average
			max_iter: maximum number of iterations the program will perform
			CTF: if this flag is set, the program will use CTF information provided in file headers
			snr: signal-to-noise ratio of the data
			rand_seed: the seed used for generating random numbers
			MPI: whether to use MPI version
		Output
			output_directory: directory name into which the output files will be written.
			header: the alignment parameters are stored in the headers of input files as 'xform.align2d'.
	"""
	import os
	result = mrefeq_ali2df(stack, refim, maskfile, ir, ou, rs, xrng, yrng, step, center, maxit, CTF, snr, user_func_name, rand_seed)
	newrefim = os.path.join(outdir,"multi_ref.hdf")
	if os.path.exists(outdir):
		from sys import exit
		exit()
	os.mkdir(outdir)
	for j in xrange(len(result)):  result[j].write_image(newrefim, j)

	return
# 2D multi-reference alignment using rotational ccf in polar coordinates and quadratic interpolation
	if MPI:
		mref_ali2d_MPI(stack, refim, outdir, maskfile, ir, ou, rs, xrng, yrng, step, center, maxit, CTF, snr, user_func_name, rand_seed)
		return

	from utilities      import   model_circle, combine_params2, inverse_transform2, drop_image, get_image
	from utilities	    import   center_2D, get_im, get_params2D, set_params2D
	from statistics     import   fsc
	from alignment      import   Numrinit, ringwe, fine_2D_refinement
	from fundamentals   import   rot_shift2D, fshift
	from morphology     import   ctf_2
	from filter         import   filt_btwl, filt_params
	from random         import   seed, randint
	import os
	import sys

	from utilities      import   print_begin_msg, print_end_msg, print_msg
	
	# create the output directory, if it does not exist
	if os.path.exists(outdir):  ERROR('Output directory exists, please change the name and restart the program', "mref_ali2d", 1)
	os.mkdir(outdir)
	
	first_ring=int(ir); last_ring=int(ou); rstep=int(rs); max_iter=int(maxit)
	if max_iter == 0:
		max_iter  = 10
		auto_stop = True
	else:
		auto_stop = False

	print_begin_msg("mrefEQ_ali2d")

	print_msg("Input stack                 : %s\n"%(stack))
	print_msg("Reference stack             : %s\n"%(refim))	
	print_msg("Output directory            : %s\n"%(outdir))
	print_msg("Maskfile                    : %s\n"%(maskfile))
	print_msg("Inner radius                : %i\n"%(first_ring))

	ima = EMData()
	ima.read_image(stack, 0)
	nx = ima.get_xsize()
	# default value for the last ring
	if last_ring == -1: last_ring = nx/2-2

	print_msg("Outer radius                : %i\n"%(last_ring))
	print_msg("Ring step                   : %i\n"%(rstep))
	print_msg("X search range              : %i\n"%(xrng))
	print_msg("Y search range              : %i\n"%(yrng))
	print_msg("Translational step          : %i\n"%(step))
	print_msg("Center type                 : %i\n"%(center))
	print_msg("Maximum iteration           : %i\n"%(max_iter))
	print_msg("CTF correction              : %s\n"%(CTF))
	print_msg("Signal-to-Noise Ratio       : %f\n"%(snr))
	print_msg("Random seed                 : %i\n\n"%(rand_seed))

	output = sys.stdout

	import user_functions
	user_func = user_functions.factory[user_func_name]

	if maskfile:
		import types
		if type(maskfile) is types.StringType:  mask = get_image(maskfile)
		else: mask = maskfile
	else: mask = model_circle(last_ring, nx, nx)
	#  references
	refi = []
	numref = EMUtil.get_image_count(refim)
	#  CTF stuff
	if CTF:
		ctf_params = ima.get_attr("ctf")
		data_had_ctf = ima.get_attr("ctf_applied")
		ctm = ctf_2(nx, ctf_params)
		lctf = len(ctm)
		ctf2 = [[[0.0]*lctf for k in xrange(2)] for j in xrange(numref)]

	# IMAGES ARE SQUARES! center is in SPIDER convention
	cnx = nx/2+1
	cny = cnx

	mode = "F"
	#precalculate rings
	numr = Numrinit(first_ring, last_ring, rstep, mode)
	wr = ringwe(numr, mode)
	# reference images
	params = []
	#read all data
	data = EMData.read_images(stack)
	nima = len(data)
	# prepare the reference
	ima.to_zero()
	for j in xrange(numref):
		temp = EMData()
		temp.read_image(refim, j)
		#  eve, odd, numer of even, number of images.  After frc, totav
		refi.append([temp, ima.copy(), 0])
	seed(rand_seed)
	a0 = -1.
	again = True
	Iter = 0

	ref_data = [mask, center, None, None]

	while Iter < max_iter and again:
		#again = False
		ringref = []
		#print "numref",numref
		for j in xrange(numref):
			refi[j][0].process_inplace("normalize.mask", {"mask":mask, "no_sigma":1})
			cimage = Util.Polar2Dm(refi[j][0], cnx, cny, numr, mode)
			Util.Frngs(cimage, numr)
			Util.Applyws(cimage, numr, wr)
			ringref.append(cimage)
			# zero refi
			refi[j][0].to_zero()
			refi[j][1].to_zero()
			refi[j][2] = 0
			if CTF:
				for i in xrange(lctf): 
					ctf2[j][0][i] = 0.0
					ctf2[j][1][i] = 0.0
		assign = [[] for i in xrange(numref)]
		sx_sum = [0.0]*numref
		sy_sum = [0.0]*numref
		peak_list = [[] for i in xrange(numref)]
		for im in xrange(nima):
			if CTF:
				ctf_params = data[im].get_attr("ctf")
				if data[im].get_attr("ctf_applied") == 0:
					st = Util.infomask(data[im], mask, False)
					data[im] -= st[0]
					from filter import filt_ctf
					data[im] = filt_ctf(data[im], ctf_params)
					data[im].set_attr('ctf_applied', 1)
			alpha, sx, sy, mirror, scale = get_params2D(data[im])
			alphai, sxi, syi, scalei = inverse_transform2(alpha, sx, sy)
			# normalize
			data[im].process_inplace("normalize.mask", {"mask":mask, "no_sigma":0})
			# align current image to the reference
			#[angt, sxst, syst, mirrort, xiref, peakt] = Util.multiref_polar_ali_2d(data[im], 
			#	ringref, xrng, yrng, step, mode, numr, cnx+sxi, cny+syi)
			#iref = int(xiref)

			temp = Util.multiref_polar_ali_2d_peaklist(data[im],
				ringref, xrng, yrng, step, mode, numr, cnx+sxi, cny+syi)
			for iref in xrange(numref):
				[alphan, sxn, syn, mn] = combine_params2(0.0, -sxi, -syi, 0, temp[iref*5+1], temp[iref*5+2], temp[iref*5+3], int(temp[iref*5+4]))
				peak_list[iref].extend([temp[iref*5], alphan, sxn, syn, mn])
		del temp


		'''
		from numpy import empty, float32
		d = empty((numref, nima), dtype = float32)
		for iref in xrange(numref):
			for im in xrange(nima):
				d[iref][im] = peak_list[iref][im*5]
		'''

		from sys import exit
		from utilities import write_text_file
		#write_text_file(peak_list,"ttt.txt")
		#exit()
		d = [0.0]*(numref*nima)
		for iref in xrange(numref):
			for im in xrange(nima):
				d[iref*nima+im] = float(peak_list[iref][im*5])
		id_list_long = Util.assign_groups(d, numref, nima)
		id_list = [[] for i in xrange(numref)]
		maxasi = nima/numref
		for i in xrange(maxasi*numref):
			id_list[i/maxasi].append(id_list_long[i])
		for i in xrange(nima%maxasi):
			id_list[id_list_long[-1]].append(id_list_long[maxasi*numref+i])
		belongsto = [0]*nima
		for im in xrange(nima):
			for iref in xrange(numref):
				try:
					i = id_list[iref].index(im)
					belongsto[im] = iref
					break
				except:
					pass

		del d

		for im in xrange(nima):
			qt = -1.0e23
			#print im,[d[iref][im] for iref in xrange(numref)],
			'''
			for iref in xrange(numref):
				if(d[iref][im] > qt):
					qt = d[iref][im]
					matchref = iref
			'''
			#print  matchref
			matchref = belongsto[im]
			# combine parameters and set them to the header, ignore previous angle and mirror
			#[alphan, sxn, syn, mn] = combine_params2(0.0, -sxi, -syi, 0, angt, sxst, syst, int(mirrort))
			alphan = peak_list[matchref][im*5+1]
			sxn    = peak_list[matchref][im*5+2]
			syn    = peak_list[matchref][im*5+3]
			mn     = peak_list[matchref][im*5+4]
			set_params2D(data[im], [alphan, sxn, syn, int(mn), scale])
			if mn == 0: sx_sum[iref] += sxn
			else: sx_sum[iref] -= sxn
			sy_sum[iref] += syn
			data[im].set_attr('assign', matchref)
			#data[im].set_attr('assign', iref)
			# apply current parameters and add to the average
			temp = rot_shift2D(data[im], alphan, sxn, syn, mn)
			it = im%2
			Util.add_img(refi[matchref][it], temp)
			if CTF:
				ctm = ctf_2(nx, ctf_params)
				for i in xrange(lctf):  ctf2[iref][it][i] += ctm[i]
			assign[matchref].append(im)
			refi[matchref][2] += 1
		del ringref
		if again:
			a1 = 0.0
			for j in xrange(numref):
				msg = "   group #%3d   number of particles = %7d\n"%(j, refi[j][2])
				print_msg(msg)
				if refi[j][2] < 4:
					#ERROR("One of the references vanished","mref_ali2d",1)
					#  if vanished, put a random image there
					assign[j] = []
					assign[j].append(randint(0, nima-1))
					refi[j][0] = data[assign[j][0]].copy()
				else:
					max_inter = 0  # switch off fine refi.
					br = 1.75
					#  the loop has to 
					for INter in xrange(max_inter+1):
						# Calculate averages at least ones, meaning even if no within group refinement was requested
						if CTF:
							for i in xrange(lctf):  ctm[i] = 1.0 / (ctf2[j][0][i] + 1.0/snr)
							from filter import filt_table
							av1 = filt_table(refi[j][0], ctm)
							for i in xrange(lctf):  ctm[i] = 1.0 / (ctf2[j][1][i] + 1.0/snr)
							av2 = filt_table(refi[j][1], ctm)
							frsc = fsc(av1, av2, 1.0, os.path.join(outdir,"drm_%03d_%04d.txt"%(Iter, j)))
							#Now the total average
							for i in xrange(lctf):  ctm[i] = 1.0 / (ctf2[j][0][i] + ctf2[j][1][i] + 1.0/snr)
							refi[j][0] = filt_table(Util.addn_img(refi[j][0], refi[j][1]), ctm)
						else:
							frsc = fsc(refi[j][0], refi[j][1], 1.0, os.path.join(outdir,"drm_%03d_%04d.txt"%(Iter, j)))
							Util.add_img(refi[j][0], refi[j][1])
							Util.mul_scalar(refi[j][0], 1.0/float(refi[j][2]))
							
						ref_data[2] = refi[j][0]
						ref_data[3] = frsc						
						refi[j][0], cs = user_func(ref_data)
						if center == -1:
							cs[0] = sx_sum[j]/len(assign[j])
							cs[1] = sy_sum[j]/len(assign[j])
							refi[j][0] = fshift(refi[j][0], -cs[0], -cs[1])
						for i in xrange(len(assign[j])):
							im = assign[j][i]
							alpha, sx, sy, mirror, scale =  get_params2D(data[im])
							alphan, sxn, syn, mirrorn = combine_params2(alpha, sx, sy, mirror, 0.0, -cs[0], -cs[1], 0)
							set_params2D(data[im], [alphan, sxn, syn, int(mirrorn), scale])
						# refine images within the group
						#  Do the refinement only if max_inter>0, but skip it for the last iteration.
						if INter < max_inter:
							fine_2D_refinement(data, br, mask, refi[j][0], j)
							#  Calculate updated average
							refi[j][0].to_zero()
							refi[j][1].to_zero()
							for i in xrange(len(assign[j])):
								im = assign[j][i]
								alpha, sx, sy, mirror, scale = get_params2D(data[im])
								# apply current parameters and add to the average
								temp = rot_shift2D(data[im], alpha, sx, sy, mn)
								it = im%2
								Util.add_img(refi[j][it], temp)

				# write the current average
				TMP = []
				for i_tmp in xrange(len(assign[j])):  TMP.append(float(assign[j][i_tmp]))
				TMP.sort()
				refi[j][0].set_attr_dict({'ave_n': refi[j][2], 'members': TMP })
				del TMP
				# replace the name of the stack with reference with the current one
				newrefim = os.path.join(outdir,"aqm%03d.hdf"%Iter)
				refi[j][0].write_image(newrefim, j)
				a1 += refi[j][0].cmp("dot", refi[j][0], {"negative":0, "mask":mask})
			Iter += 1
			msg = "ITERATION #%3d        criterion = %20.7e\n"%(Iter,a1)
			print_msg(msg)
			if a1 < a0:
				if auto_stop == True:	break
			else:	a0 = a1

	newrefim = os.path.join(outdir,"multi_ref.hdf")
	for j in xrange(numref):  refi[j][0].write_image(newrefim, j)
	if CTF:
		if data_had_ctf == 0:
			for im in xrange(nima): data[im].set_attr('ctf_applied', 0)
	from utilities import write_headers
	write_headers(stack, data, range(nima))
	print_end_msg("mrefEQ_ali2d")


def mrefeq_ali2df(stack, refim, maskfile=None, ir=1, ou=-1, rs=1, xrng=0, yrng=0, step=1, center=1, maxit=0, CTF=False, snr=1.0, user_func_name="ref_ali2d", rand_seed=1000, MPI=False):
	"""
		Name
			mref_ali2d - Perform 2-D multi-reference alignment of an image series
		Input
			stack: set of 2-D images in a stack file, images have to be squares
			refim: set of initial reference 2-D images in a stack file 
			maskfile: optional maskfile to be used in the alignment
			inner_radius: inner radius for rotational correlation > 0
			outer_radius: outer radius for rotational correlation < nx/2-1
			ring_step: step between rings in rotational correlation >0
			x_range: range for translation search in x direction, search is +/xr 
			y_range: range for translation search in y direction, search is +/yr 
			translation_step: step of translation search in both directions
			center: center the average
			max_iter: maximum number of iterations the program will perform
			CTF: if this flag is set, the program will use CTF information provided in file headers
			snr: signal-to-noise ratio of the data
			rand_seed: the seed used for generating random numbers
			MPI: whether to use MPI version
		Output
			output_directory: directory name into which the output files will be written.
			header: the alignment parameters are stored in the headers of input files as 'xform.align2d'.
	"""
	# 2D multi-reference alignment using rotational ccf in polar coordinates and quadratic interpolation
	if MPI:
		mrefeq_ali2df_MPI(stack, refim, maskfile, ir, ou, rs, xrng, yrng, step, center, maxit, CTF, snr, user_func_name, rand_seed)
		return
	from utilities      import   model_circle, combine_params2, inverse_transform2, drop_image, get_image
	from utilities	    import   center_2D, get_im, get_params2D, set_params2D, model_blank
	from statistics     import   fsc
	from alignment      import   Numrinit, ringwe, fine_2D_refinement
	from fundamentals   import   rot_shift2D, fshift
	from morphology     import   ctf_2
	from filter         import   filt_tanl, filt_params
	from random         import   seed, randint
	import os
	import sys
	from time import localtime
	from utilities      import   print_begin_msg, print_end_msg, print_msg
	
	first_ring=int(ir); last_ring=int(ou); rstep=int(rs); max_iter=int(maxit)
	if max_iter == 0:
		max_iter  = 10
		auto_stop = True
	else:
		auto_stop = False

	if(type(stack) == type("")):
		#read all data
		data = EMData.read_images(stack)
	else:
		data = stack
	nima = len(data)
	nx = data[0].get_xsize()
	# default value for the last ring
	if last_ring == -1: last_ring = nx/2-2

	#import user_functions
	#user_func = user_functions.factory[user_func_name]
	seed(rand_seed)
	if maskfile:
		import types
		if type(maskfile) is types.StringType:  mask = get_image(maskfile)
		else: mask = maskfile
	else: mask = model_circle(last_ring, nx, nx)
	#  references
	if(type(refim) == type("")):
		refi = EMData.read_images(refim)
	else:
		refi = refim
	numref = len(refi)
	from random import shuffle
	numr = range(nima)
	shuffle(numr)
	for iref in xrange(numref):
		refi[iref] = data[numr[iref]].copy()
	#  CTF stuff
	if CTF:
		ctf_params = data[0].get_attr("ctf")
		data_had_ctf = data[0].get_attr("ctf_applied")
		ctm = ctf_2(nx, ctf_params)
		lctf = len(ctm)
		ctf2 = [[[0.0]*lctf for k in xrange(2)] for j in xrange(numref)]

	# IMAGES ARE SQUARES! center is in SPIDER convention
	cnx = nx/2+1
	cny = cnx

	mode = "F"
	#precalculate rings
	numr = Numrinit(first_ring, last_ring, rstep, mode)
	wr = ringwe(numr, mode)
	# reference images
	params = []
	a0 = -1.0
	again = True
	Iter = 0

	#ref_data = [mask, center, None, None]

	while Iter < max_iter and again:
		print  "ITERATION ",Iter,max_iter,localtime()[:5]
		#again = False
		ringref = []
		#print "numref",numref
		for j in xrange(numref):
			refi[j].process_inplace("normalize.mask", {"mask":mask, "no_sigma":1})
			cimage = Util.Polar2Dm(refi[j], cnx, cny, numr, mode)
			Util.Frngs(cimage, numr)
			Util.Applyws(cimage, numr, wr)
			ringref.append(cimage)
			if CTF:
				for i in xrange(lctf): 
					ctf2[j][0][i] = 0.0
					ctf2[j][1][i] = 0.0
		peak_list = [[] for i in xrange(numref)]
		d = [0.0]*(numref*nima)
		for im in xrange(nima):
			if CTF:
				ctf_params = data[im].get_attr("ctf")
				if data[im].get_attr("ctf_applied") == 0:
					st = Util.infomask(data[im], mask, False)
					data[im] -= st[0]
					from filter import filt_ctf
					data[im] = filt_ctf(data[im], ctf_params)
					data[im].set_attr('ctf_applied', 1)
			alpha, sx, sy, mirror, scale = get_params2D(data[im])
			alphai, sxi, syi, scalei = inverse_transform2(alpha, sx, sy)
			# normalize
			data[im].process_inplace("normalize.mask", {"mask":mask, "no_sigma":0})
			# align current image to the reference

			temp = Util.multiref_polar_ali_2d_peaklist(data[im],
				ringref, xrng, yrng, step, mode, numr, cnx+sxi, cny+syi)
			for iref in xrange(numref):
				[alphan, sxn, syn, mn] = \
				   combine_params2(0.0, -sxi, -syi, 0, temp[iref*5+1], temp[iref*5+2], temp[iref*5+3], int(temp[iref*5+4]))
				peak_list[iref].extend([alphan, sxn, syn, mn])
				d[iref*nima+im] = temp[iref*5]
		del ringref
		del temp

		id_list_long = Util.assign_groups(d, numref, nima)
		id_list = [[] for i in xrange(numref)]
		maxasi = nima/numref
		for i in xrange(maxasi*numref):
			id_list[i/maxasi].append(id_list_long[i])
		for i in xrange(nima%maxasi):
			id_list[id_list_long[-1]].append(id_list_long[maxasi*numref+i])
		belongsto = [0]*nima
		for im in xrange(nima):
			for iref in xrange(numref):
				try:
					i = id_list[iref].index(im)
					belongsto[im] = iref
					break
				except:
					pass

		del id_list
		del d
		members = [0]*numref
		assign = [[] for i in xrange(numref)]
		sx_sum = [0.0]*numref
		sy_sum = [0.0]*numref

		for im in xrange(nima):
			matchref = belongsto[im]
			# combine parameters and set them to the header, ignore previous angle and mirror
			#[alphan, sxn, syn, mn] = combine_params2(0.0, -sxi, -syi, 0, angt, sxst, syst, int(mirrort))
			alphan = peak_list[matchref][im*4+0]
			sxn    = peak_list[matchref][im*4+1]
			syn    = peak_list[matchref][im*4+2]
			mn     = peak_list[matchref][im*4+3]
			set_params2D(data[im], [alphan, sxn, syn, int(mn), scale])
			if mn == 0: sx_sum[iref] += sxn
			else:       sx_sum[iref] -= sxn
			sy_sum[iref] += syn
			data[im].set_attr('assign', matchref)
			# apply current parameters and add to the average
			Util.add_img(refi[matchref], rot_shift2D(data[im], alphan, sxn, syn, mn))
			if CTF:
				ctm = ctf_2(nx, ctf_params)
				for i in xrange(lctf):  ctf2[iref][it][i] += ctm[i]
			assign[matchref].append(im)
			members[matchref] += 1
		del ringref
		if again:
			for j in xrange(numref):  refi[j].write_image("resa%03d.hdf"%Iter,j)
			a1 = 0.0
			for j in xrange(numref):
				if members[j] < 5:
					#  if vanished, put a random image there
					assign[j] = []
					assign[j].append(randint(0, nima-1))
					refi[j] = data[assign[j][0]].copy()
				else:
					#  Golden rule when to do within group refinement
					if(Iter<4 or Iter%3 != 0):
						if CTF:
							for i in xrange(lctf):  ctm[i] = 1.0 / (ctf2[j][0][i] + 1.0/snr)
							from filter import filt_table
							av1 = filt_table(refi[j][0], ctm)
							for i in xrange(lctf):  ctm[i] = 1.0 / (ctf2[j][1][i] + 1.0/snr)
							av2 = filt_table(refi[j][1], ctm)
							#frsc = fsc(av1, av2, 1.0, os.path.join(outdir,"drm_%03d_%04d.txt"%(Iter, j)))
							#Now the total average
							for i in xrange(lctf):  ctm[i] = 1.0 / (ctf2[j][0][i] + ctf2[j][1][i] + 1.0/snr)
							refi[j][0] = filt_table(Util.addn_img(refi[j][0], refi[j][1]), ctm)
						else:
							#frsc = fsc(refi[j][0], refi[j][1], 1.0, os.path.join(outdir,"drm_%03d_%04d.txt"%(Iter, j)))
							#remove odd-even?
							Util.mul_scalar(refi[j], 1.0/float(members[j]))

						#ref_data[2] = refi[j]
						#ref_data[3] = frsc						
						#refi[j], cs = user_func(ref_data)
						center = -1
						cs = [0.0, 0.0]
						if center == -1:
							cs[0] = sx_sum[j]/len(assign[j])
							cs[1] = sy_sum[j]/len(assign[j])
							refi[j] = fshift(refi[j], -cs[0], -cs[1])
						for i in xrange(len(assign[j])):
							im = assign[j][i]
							alpha, sx, sy, mirror, scale =  get_params2D(data[im])
							alphan, sxn, syn, mirrorn = combine_params2(alpha, sx, sy, mirror, 0.0, -cs[0], -cs[1], 0)
							set_params2D(data[im], [alphan, sxn, syn, int(mirrorn), scale])
					else:
						# refine images within the group
						#  Do the refinement only if max_inter>0, but skip it for the last iteration.
						print "within group", j,localtime()[:5]
						randomize = True
						refi[j] = within_group_refinement([data[im] for im in assign[j]], mask, randomize, 
						          ir, ou, rstep, [xrng], [yrng], [step], dst=0.0, maxit=15, FH = 0.3)
				refi[j] = filt_tanl( refi[j], 0.1 + (0.15/3.0)*(Iter%3), 0.2)
				# write the current average
				TMP = []
				for i_tmp in xrange(len(assign[j])):  TMP.append(float(assign[j][i_tmp]))
				TMP.sort()
				refi[j].set_attr_dict({'ave_n': refi[j], 'members': TMP, 'n_objects': members[j] })
				del TMP
				# replace the name of the stack with reference with the current one
				#newrefim = os.path.join(outdir,"aqm%03d.hdf"%Iter)
				#refi[j][0].write_image(newrefim, j)
				a1 += refi[j].cmp("dot", refi[j], {"negative":0, "mask":mask})
			#msg = "ITERATION #%3d        criterion = %20.7e\n"%(Iter,a1)
			#print_msg(msg)
			if a1 < a0:
				if auto_stop == True:	break
			else:	a0 = a1
		for j in xrange(numref):  refi[j].write_image("resi%03d.hdf"%Iter,j)
		Iter += 1

	if CTF:
		if data_had_ctf == 0:
			for im in xrange(nima): data[im].set_attr('ctf_applied', 0)
	return  [refi[j] for j in xrange(numref)]


def mrefeq_ali2df_MPI(stack, refim, maskfile = None, ir=1, ou=-1, rs=1, xrng=0, yrng=0, step=1, center=1, maxit=10, CTF=False, snr=1.0, user_func_name="ref_ali2d", rand_seed=1000, stability = False):
# 2D multi-reference alignment using rotational ccf in polar coordinates and quadratic interpolation

	from utilities      import   model_circle, combine_params2, inverse_transform2, drop_image, get_image, get_im
	from utilities      import   reduce_EMData_to_root, bcast_EMData_to_all, bcast_number_to_all
	from applications   import   MPI_start_end
	from utilities      import   send_attr_dict
	from utilities	    import   center_2D
	from statistics     import   fsc_mask
	from alignment      import   Numrinit, ringwe
	from fundamentals   import   rot_shift2D, fshift
	from utilities      import   get_params2D, set_params2D
	from random         import   seed, randint
	from morphology     import   ctf_2
	from filter         import   filt_tanl, filt_params
	from numpy          import   reshape, shape
	from utilities      import   print_msg, print_begin_msg, print_end_msg
	import os
	import sys
	from mpi 	    import mpi_comm_size, mpi_comm_rank, MPI_COMM_WORLD, mpi_barrier
	from mpi 	    import mpi_reduce, mpi_bcast, mpi_barrier, mpi_recv, mpi_send
	from mpi 	    import MPI_SUM, MPI_FLOAT, MPI_INT, MPI_TAG_UB

	if comm == -1: comm = MPI_COMM_WORLD		

	number_of_proc = mpi_comm_size(comm)
	myid = mpi_comm_rank(comm)
	main_node = 0
	

	first_ring=int(ir); last_ring=int(ou); rstep=int(rs); max_iter=int(maxit)
	if max_iter == 0:
		max_iter  = 10
		auto_stop = True
	else:
		auto_stop = False

	if(type(stack) == type("")):
		#read all data
		alldata = EMData.read_images(stack)
	else:
		alldata = stack
	nx = alldata[0].get_xsize()
	# default value for the last ring
	if last_ring == -1: last_ring = nx/2-2

	nima = len(alldata)
	
	image_start, image_end = MPI_start_end(nima, number_of_proc, myid)
	data = [None]*(image_end-image_start)
	for im in xrange(image_start, image_end):
		data[im-image_start] = alldata[im]


	first_ring=int(ir); last_ring=int(ou); rstep=int(rs); max_iter=int(maxit)
	'''
	#  delete this, assume refim is always given
	if maskfile:
		import  types
		if type(maskfile) is types.StringType:  mask = get_image(maskfile)
		else: mask = maskfile
	else : mask = model_circle(last_ring, nx, nx)
	if(type(refim) == type("")):
		refi = EMData.read_images(refim)
	else:
		refi = refim
	numref = len(refi)
	if(  myid == main_node):
		from time import localtime
		from random import shuffle
		numr = range(nima)
		shuffle(numr)
		for iref in xrange(numref):
			refi[iref] = alldata[numr[iref]].copy()
		del numr
	for j in xrange(numref):
		bcast_EMData_to_all(refi[j], myid, main_node)
	'''
	#  CTF stuff
	if CTF:
		ctf_params = ima.get_attr("ctf")
		data_had_ctf = ima.get_attr("ctf_applied")
		ctm = ctf_2(nx, ctf_params)
		lctf = len(ctm)

	# IMAGES ARE SQUARES! center is in SPIDER convention
	cnx = nx/2+1
	cny = cnx

	mode = "F"
	#precalculate rings
	numr = Numrinit(first_ring, last_ring, rstep, mode)
	wr = ringwe(numr, mode)
	# reference images
	again = True
	params = []
	#  for each node read its share of data
	data = EMData.read_images(stack, range(image_start, image_end))
	for im in xrange(image_start, image_end):
		data[im-image_start].set_attr('ID', im)
		if CTF:
			ctf_params = data[im-image_start].get_attr( "ctf" )
			if data[im-image_start].get_attr("ctf_applied") == 0:
				st = Util.infomask(data[im-image_start], mask, False)
				data[im-image_start] -= st[0]
				from filter import filt_ctf
				data[im-image_start] = filt_ctf(data[im-image_start], ctf_params)
				data[im-image_start].set_attr('ctf_applied', 1)
	if myid == main_node:  seed(rand_seed)

	a0 = -1.0
	again = True
	Iter = -1


	'''
	# MPI stuff
	if myid == main_node:
		recvcount = [ndata]
		disp = []
		for n in xrange(number_of_proc):
			if n != main_node:
				lnn = mpi_recv(1, MPI_INT, n, n*100, MPI_COMM_WORLD)
				lnn = int(lnn[0])
				recvcount.append(lnn)
			if n == 0: disp.append(0)
			else: disp.append(disp[n-1]+recvcount[n-1])
	else:
		mpi_send(ndata, 1, MPI_INT, main_node, myid*100, MPI_COMM_WORLD)
		recvcount = [0.0]*number_of_proc
		disp = [0.0]*number_of_proc

	recvcount = mpi_bcast(recvcount, number_of_proc, MPI_INT, main_node, MPI_COMM_WORLD)  
	disp = mpi_bcast(disp, number_of_proc, MPI_INT, main_node, MPI_COMM_WORLD)  
	recvcount = map(int, recvcount)
	disp = map(int, disp)
	'''

	ref_data = [mask, center, None, None]
	fl = 0.1
	FH = 0.3
	main_iter = 0
	while( main_iter < max_iter ):
		Iter += 1
		if myid == main_node: print Iter,len(data),localtime()[0:5]
		ringref = []
		for j in xrange(numref):
			refi[j].process_inplace("normalize.mask", {"mask":mask, "no_sigma":1}) # normalize reference images to N(0,1)
			cimage = Util.Polar2Dm(refi[j] , cnx, cny, numr, mode)
			Util.Frngs(cimage, numr)
			Util.Applyws(cimage, numr, wr)
			ringref.append(cimage)
		if CTF: ctf2 = [[[0.0]*lctf for k in xrange(2)] for j in xrange(numref)]
		peak_list = [[] for i in xrange(numref)]
		#  nima is the total number of images, not the one on this node, tha latter is (image_end-image_start)
		d = [0.0]*(numref*nima)
		# begin MPI section
		for im in xrange(image_start, image_end):
			alpha, sx, sy, mirror, scale = get_params2D(alldata[im])
			alphai, sxi, syi, scalei = inverse_transform2(alpha, sx, sy)
			# normalize
			alldata[im].process_inplace("normalize.mask", {"mask":mask, "no_sigma":0}) # subtract average under the mask
			# align current image to the reference

			temp = Util.multiref_polar_ali_2d_peaklist(alldata[im],
				ringref, xrng, yrng, step, mode, numr, cnx+sxi, cny+syi)
			for iref in xrange(numref):
				[alphan, sxn, syn, mn] = \
				   combine_params2(0.0, -sxi, -syi, 0, temp[iref*5+1], temp[iref*5+2], temp[iref*5+3], int(temp[iref*5+4]))
				peak_list[iref].extend([alphan, sxn, syn, mn])
				d[iref*nima+im] = temp[iref*5]
		del ringref
		del temp


		d = mpi_reduce(d, numref*nima, MPI_FLOAT, MPI_SUM, main_node, comm)
		if myid == main_node:
			print  Iter,localtime()[0:5]
			d = map(float, d)
			id_list_long = Util.assign_groups(d, numref, nima)
			id_list = [[] for i in xrange(numref)]
			maxasi = nima/numref
			for i in xrange(maxasi*numref):
				id_list[i/maxasi].append(id_list_long[i])
			for i in xrange(nima%maxasi):
				id_list[id_list_long[-1]].append(id_list_long[maxasi*numref+i])
			for iref in xrange(numref):
				id_list[iref].sort()

			belongsto = [0]*nima
			for im in xrange(nima):
				for iref in xrange(numref):
					try:
						i = id_list[iref].index(im)
						belongsto[im] = iref
						break
					except:
						pass
			del id_list
		else:
			belongsto = [0]*nima
		mpi_barrier(comm)
		belongsto = mpi_bcast(belongsto, nima, MPI_INT, main_node, comm)
		belongsto = map(int, belongsto)
		if myid == main_node: print " belongsto  ",myid, localtime()[:5]
		#  Compute partial averages
		members = [0]*numref
		sx_sum = [0.0]*numref
		sy_sum = [0.0]*numref
		for j in xrange(numref):  refi[j].to_zero()
		for im in xrange(image_start, image_end):
			matchref = belongsto[im]
			# combine parameters and set them to the header, ignore previous angle and mirror
			#[alphan, sxn, syn, mn] = combine_params2(0.0, -sxi, -syi, 0, angt, sxst, syst, int(mirrort))
			alphan = peak_list[matchref][(im-image_start)*4+0]
			sxn    = peak_list[matchref][(im-image_start)*4+1]
			syn    = peak_list[matchref][(im-image_start)*4+2]
			mn     = peak_list[matchref][(im-image_start)*4+3]
			set_params2D(alldata[im], [alphan, sxn, syn, int(mn), scale])
			if mn == 0: sx_sum[iref] += sxn
			else:       sx_sum[iref] -= sxn
			sy_sum[iref] += syn
			alldata[im].set_attr('assign', matchref)
			# apply current parameters and add to the average
			Util.add_img(refi[matchref], rot_shift2D(alldata[im], alphan, sxn, syn, mn))
			if CTF:
				ctm = ctf_2(nx, ctf_params)
				for i in xrange(lctf):  ctf2[iref][it][i] += ctm[i]
			members[matchref] += 1
		sx_sum = mpi_reduce(sx_sum, numref, MPI_FLOAT, MPI_SUM, main_node, comm)
		sx_sum = map(float, sx_sum)
		sy_sum = mpi_reduce(sy_sum, numref, MPI_FLOAT, MPI_SUM, main_node, comm)
		sy_sum = map(float, sy_sum)
		members = mpi_reduce(members, numref, MPI_INT, MPI_SUM, main_node, comm)
		members = map(int, members)
		if  myid == main_node:  print " averages done  ",localtime()[0:5]
		for j in xrange(numref):
			reduce_EMData_to_root(refi[j], myid, main_node)
			if( myid == main_node):
				if members[j] < 5:
					#  if vanished, put a random image there
					assign[j] = []
					assign[j].append(randint(0, nima-1))
					refi[j] = alldata[assign[j][0]].copy()
				else:
					#  Golden rule when to do within group refinement
					Util.mul_scalar(refi[j], 1.0/float(members[j]))
					refi[j] = filt_tanl( refi[j], fl, 0.2)
				#refi[j].write_image("avim%03d.hdf"%Iter,j)
			bcast_EMData_to_all(refi[j], myid, main_node)
		do_within_group = 0
		if(myid == main_node):
			fl += 0.05
			if(fl == FH):
				fl = 0.1
				do_within_group = 1
			
		do_within_group = mpi_bcast(do_within_group, 1, MPI_INT, main_node, MPI_COMM_WORLD)
		do_within_group = int(do_within_group)
			


		if(do_within_group == 1):
			#  Here the assumption is that numref > number_of_proc, check it.
			if stability:
				ref1 = [None]*numref
				ref2 = [None]*numref
				# and so on
			main_iter += 1
			for j in xrange(myid,numref,number_of_proc):
				assign = []
				for im in xrange(nima):
					if(j == belongsto[im]):  assign.append(im)
				if(len(assign) >4):  # this makes no sense, as all groups should have the same number of assigned elements, please try to fix it.  PAP.
						if myid == main_node:  print "within group", j,localtime()[:5]
						randomize = True
						refi[j] = within_group_refinement([alldata[im] for im in assign], mask, randomize, ir, ou, rstep,
						                                 [xrng], [yrng], [step], dst=00.0, maxit=15, FH = 0.3)
						if stability and main_iter%2==0:
							ref1[j] = within_group_refinement()
							ref2[j] = within_group_refinement()
							# and so on
				del assign
			mpi_barrier(comm)
			from utilities import recv_EMData, send_EMData 
			for sts in xrange(1,number_of_proc):
				for j in xrange(sts,numref,number_of_proc):
					if( myid == 0):
						#  receive
						#print " receive ",j,sts
						refi[j] = recv_EMData(sts, 17+3*sts)
					elif(myid == sts):
						# send
						#print  "  send  ",j,sts,myid
						send_EMData(refi[j], main_node, 17+3*sts)
					mpi_barrier(comm)
			for j in xrange(numref):
				bcast_EMData_to_all(refi[j], myid, main_node)
			if stability:
				a = a
				# and so on...
		if( myid == main_node):
			for j in xrange(numref):
				refi[j].write_image("avim%03d.hdf"%Iter,j)
	return refi

def ihrsr_n(stack, ref_vol, outdir, maskfile, ir, ou, rs, xr, ynumber, 
          txs, delta, initial_theta, delta_theta, an, maxit, CTF, snr, dp, ndp, dp_step, dphi, ndhpi, dphi_step, psi_max,
	  rmin, rmax, fract, nise, npad, sym, user_func_name, datasym,
	  fourvar, debug = False, MPI = False, chunk = 0.2):
	if MPI:
		from development import ihrsr_n_MPI
		ihrsr_n_MPI(stack, ref_vol, outdir, maskfile, ir, ou, rs, xr, ynumber, 
			txs, delta, initial_theta, delta_theta, an, maxit, CTF, snr, dp, ndp, dp_step, dphi, ndhpi, dphi_step, psi_max,
			rmin, rmax, fract, nise, npad, sym, user_func_name, datasym,
			fourvar, debug, chunk)
		return
	print_end_msg("ihrsr_n")

def ihrsr_n_MPI(stack, ref_vol, outdir, maskfile, ir, ou, rs, xr, ynumber, 
	txs, delta, initial_theta, delta_theta, an, maxit, CTF, snr, dp, ndp, dp_step, dphi, ndphi, dphi_step, psi_max,
	rmin, rmax, fract, nise, npad, sym, user_func_name, datasym,
	fourvar, debug, chunk):

	from alignment      import Numrinit, prepare_refrings, proj_ali_helical, proj_ali_helical_90, proj_ali_helical_local, proj_ali_helical_90_local, helios,helios7
	from utilities      import model_circle, get_image, drop_image, get_input_from_string
	from utilities      import bcast_list_to_all, bcast_number_to_all, reduce_EMData_to_root, bcast_EMData_to_all
	from utilities      import send_attr_dict
	from utilities      import get_params_proj, set_params_proj, file_type
	from fundamentals   import rot_avg_image
	from pixel_error    import max_3D_pixel_error
	import os
	import types
	from utilities      import print_begin_msg, print_end_msg, print_msg
	from mpi            import mpi_bcast, mpi_comm_size, mpi_comm_rank, MPI_FLOAT, MPI_COMM_WORLD, mpi_barrier
	from mpi            import mpi_recv,  mpi_send, MPI_TAG_UB
	from mpi            import mpi_reduce, MPI_INT, MPI_SUM
	from filter         import filt_ctf
	from projection     import prep_vol, prgs
	from statistics     import hist_list, varf3d_MPI
	from applications   import MPI_start_end
	from EMAN2 import Vec2f
	from string    import lower,split
	from math import cos, pi
	from random import shuffle

	number_of_proc = mpi_comm_size(MPI_COMM_WORLD)
	myid           = mpi_comm_rank(MPI_COMM_WORLD)
	main_node = 0
	
	if myid == 0:
		if os.path.exists(outdir):  nx = 1
		else:  nx = 0
	else:  nx = 0
	ny = bcast_number_to_all(nx, source_node = main_node)
	
	if ny == 1:  ERROR('Output directory exists, please change the name and restart the program', "ihrsr_MPI", 1,myid)
	mpi_barrier(MPI_COMM_WORLD)

	if myid == main_node:
		os.mkdir(outdir)
		import global_def
		global_def.LOGFILE =  os.path.join(outdir, global_def.LOGFILE)
	mpi_barrier(MPI_COMM_WORLD)
	
	'''ndp    = 12
	sndp   = 0.1
	ndphi  = 12
	sndphi = 0.1'''
	nlprms = (2*ndp+1)*(2*ndphi+1)
	if nlprms< number_of_proc:
		ERROR('number of cpus is larger than the number of helical search, please reduce it or at this moment modify ndp,dphi in the program', "ihrsr_MPI", 1,myid)
	


	if debug:
		from time import sleep
		while not os.path.exists(outdir):
			print  "Node ",myid,"  waiting..."
			sleep(5)

		info_file = os.path.join(outdir, "progress%04d"%myid)
		finfo = open(info_file, 'w')
	else:
		finfo = None


	#sym = "c1"
	symref = "s"+sym

	ref_a= "P"
	symmetryLower = sym.lower()
	symmetry_string = split(symmetryLower)[0]

	xrng        = get_input_from_string(xr)
	ynumber	    = get_input_from_string(ynumber)
	for i in xrange(len(ynumber)):
		if(ynumber[i]%2==1):
			ynumber[i]=ynumber[i]+1
	yrng =[]
	
	for i in xrange(len(xrng)):
		yrng.append(dp/2)
	
	stepx        = get_input_from_string(txs)
	delta       = get_input_from_string(delta)
	lstp = min(len(xrng), len(yrng), len(stepx), len(delta))
	if an == "-1":
		an = [-1] * lstp
	else:
		an = get_input_from_string(an)

	first_ring  = int(ir)
	rstep       = int(rs)
	last_ring   = int(ou)
	max_iter    = int(maxit)

	vol     = EMData()
	vol.read_image(ref_vol)
	nx      = vol.get_xsize()
	ny      = vol.get_ysize()
	nz      = vol.get_zsize()
	if( nx==nz & ny==nz):
		xysize = -1
	else:
		xysize=nx
	if last_ring < 0:	last_ring = int(nx/2) - 2

	if myid == main_node:
		import user_functions
		user_func = user_functions.factory[user_func_name]

		print_msg("Input stack                               : %s\n"%(stack))
		print_msg("Reference volume                          : %s\n"%(ref_vol))	
		print_msg("Output directory                          : %s\n"%(outdir))
		print_msg("Maskfile                                  : %s\n"%(maskfile))
		print_msg("Inner radius                              : %i\n"%(first_ring))
		print_msg("Outer radius                              : %i\n"%(last_ring))
		print_msg("Ring step                                 : %i\n"%(rstep))
		print_msg("X search range                            : %s\n"%(xrng))
		print_msg("Y number                                  : %s\n"%(ynumber))
		print_msg("Translational stepx                        : %s\n"%(stepx))
		print_msg("Angular step                              : %s\n"%(delta))
		print_msg("Angular search range                      : %s\n"%(an))
		print_msg("min radius for helical search (in pix)    : %5.4f\n"%(rmin))
		print_msg("max radius for helical search (in pix)    : %5.4f\n"%(rmax))
		print_msg("fraction of volume used for helical search: %5.4f\n"%(fract))
		print_msg("initial symmetry - angle                  : %5.4f\n"%(dphi))
		print_msg("initial symmetry - axial rise             : %5.4f\n"%(dp))
		print_msg("Maximum iteration                         : %i\n"%(max_iter))
		print_msg("Data with CTF                             : %s\n"%(CTF))
		print_msg("Signal-to-Noise Ratio                     : %5.4f\n"%(snr))
		print_msg("symmetry output doc file                  : %s\n"%(datasym))
		print_msg("number of times initial symmetry is imposed: %i\n"%(nise))
		print_msg("npad                                      : %i\n"%(npad))
		print_msg("User function                             : %s\n"%(user_func_name))
		print_msg("chunk                                     : %f\n"%(chunk))

	if maskfile:
		if type(maskfile) is types.StringType: mask3D = get_image(maskfile)
		else:                                  mask3D = maskfile
	#else: mask3D = model_circle(last_ring, nx, nx, nx)

	numr	= Numrinit(first_ring, last_ring, rstep, "F")
	mask2D  = model_circle(last_ring,nz,nz) - model_circle(first_ring,nz,nz)

	if CTF:
		from reconstruction import recons3d_4nn_ctf_MPI
		from filter         import filt_ctf
	else:	 from reconstruction import recons3d_4nn_MPI

	if myid == main_node:
		if(file_type(stack) == "bdb"):
			from EMAN2db import db_open_dict
			dummy = db_open_dict(stack, True)
				
		# horatio active_refactoring Jy51i1EwmLD4tWZ9_00000_1
		# active = EMUtil.get_all_attributes(stack, 'active')
		# list_of_particles = []
		# for im in xrange(len(active)):
		# 	if active[im]:  list_of_particles.append(im)
		# del active
		# nima = len(list_of_particles)
		nima = EMUtil.get_image_count(stack)
		list_of_particles = range(nima)

	else:
		nima = 0
	total_nima = bcast_number_to_all(nima, source_node = main_node)

	if myid != main_node:
		list_of_particles = [-1]*total_nima
	list_of_particles = bcast_list_to_all(list_of_particles, source_node = main_node)

	image_start, image_end = MPI_start_end(total_nima, number_of_proc, myid)
	# create a list of images for each node
	list_of_particles = list_of_particles[image_start: image_end]
	nima = len(list_of_particles)
	if debug:
		finfo.write("image_start, image_end: %d %d\n" %(image_start, image_end))
		finfo.flush()

	data = EMData.read_images(stack, list_of_particles)
	if fourvar:  original_data = []
	for im in xrange(nima):
		data[im].set_attr('ID', list_of_particles[im])
		if fourvar: original_data.append(data[im].copy())
		if CTF:
			st = Util.infomask(data[im], mask2D, False)
			data[im] -= st[0]
			st = data[im].get_attr_default("ctf_applied", 0)
			if(st == 0):
				ctf_params = data[im].get_attr("ctf")
				data[im] = filt_ctf(data[im], ctf_params)
				data[im].set_attr('ctf_applied', 1)

	if debug:
		finfo.write( '%d loaded  \n' % nima )
		finfo.flush()
	if CTF:
		pixel_size = data[0].get_attr('ctf').apix
	else:
		pixel_size = data[0].get_attr('pixel_size')
	for i in xrange(len(xrng)):
		yrng[i]=dp/(2*pixel_size)
		
	if myid == main_node:
		print_msg("Pixel size in Angstroms                   : %5.4f\n\n"%(pixel_size))
		print_msg("Y search range (pix) initialized as       : %s\n\n"%(yrng))
	#  this is needed for determining chunk, each chunk data used for aligment--anima, data used for reconstruction--rnima
	nchunk = int(1.0/chunk+0.5) 
	chunk_list=[0]*2*nchunk
	anima = [0]*nchunk
	rnima =[0]*nchunk
	for i in xrange(nchunk):
		chunk_list[2*i+0] = int(round(float(nima)/nchunk*i))
		chunk_list[2*i+1]   = int(round(float(nima)/nchunk*(i+1)))
		anima[i] = chunk_list[2*i+1]-chunk_list[2*i+0]
		rnima[i] = nima - anima[i]
	

	from time import time	

	#  this is needed for gathering of pixel errors
	disps = []
	recvcount = []
	for im in xrange(number_of_proc):
		if( im == main_node ):  disps.append(0)
		else:                  disps.append(disps[im-1] + recvcount[im-1])
		ib, ie = MPI_start_end(total_nima, number_of_proc, im)
		recvcount.append( ie - ib )

	#jeanmod
	total_iter = 0
	# do the projection matching
	
	for N_step in xrange(lstp):
		terminate = 0
		Iter = 0
 		while(Iter < max_iter and terminate == 0):
			yrng[N_step]=float(dp)/(2*pixel_size) #will change it later according to dp
			if(ynumber[N_step]==0):
				stepy=0.0
			else:
				stepy=(2*yrng[N_step]/ynumber[N_step])
				
			pixer  = [0.0]*nima
			modphi = [0.0]*nima
			Iter += 1
			total_iter += 1
			if myid == main_node:
				start_time = time()
				print_msg("\nITERATION #%3d,  inner iteration #%3d\nDelta = %4.1f, an = %5.4f, xrange = %5.4f,stepx = %5.4f, yrange = %5.4f,  stepy = %5.4f, ynumber = %3d\n"%(total_iter, Iter, delta[N_step], an[N_step], xrng[N_step],stepx[N_step],yrng[N_step],stepy,ynumber[N_step]))
			c1 = range(nima)
			shuffle(c1)
			c1_set = set(c1)
			for nch in xrange(nchunk):
				# data is all the input data of each processor, len(data) = nima
				# rdata is data used for reconstruction, len(rdata) = rnima
				# adata is the data used for aligment, len(data) = anima
				list1 = []
				adata = [None]* anima[nch]
				rdata = [None]*rnima[nch]
	
				for k in xrange( len(adata) ):
					adata[k] = data[c1[chunk_list[2*nch]+k]]
					list1.append( c1[chunk_list[2*nch]+k] )
	
				list1_set = set( list1 )
				list2 = list(  c1_set.difference(list1_set) )
				for k in xrange( len(rdata) ):   
					rdata[k]= data[ list2[k] ]
				if myid == main_node:
					start_time = time()
				if CTF: vol = recons3d_4nn_ctf_MPI(myid, rdata, symmetry=sym, snr = snr, npad = npad, xysize = xysize)
				else:    vol = recons3d_4nn_MPI(myid, rdata, symmetry=sym, npad = npad, xysize = xysize)
				del rdata, list2
				if myid==main_node:
					#  for each chunk, we directly impose helical parameters
					vol = vol.helicise(pixel_size,dp, dphi, fract, rmax, rmin)
					ref_data = [vol]
					vol = user_func(ref_data)
					vol = vol.helicise(pixel_size,dp, dphi, fract, rmax, rmin)


				bcast_EMData_to_all(vol, myid, main_node)

				if( xysize == -1 ):
					volft,kb = prep_vol( vol )
					refrings = prepare_refrings( volft, kb, nz, delta[N_step], ref_a, symref, numr, MPI = True, phiEqpsi = "Zero", initial_theta =initial_theta, delta_theta = delta_theta)
					del volft,kb
				else:
					volft, kbx, kby, kbz = prep_vol( vol )
					refrings = prepare_refrings( volft, kbz, nz, delta[N_step], ref_a, symref, numr, MPI = True, phiEqpsi = "Zero", kbx = kbx, kby = kby, initial_theta =initial_theta, delta_theta = delta_theta)
					del volft, kbx, kby, kbz
				
				#split refrings to two list: refrings1 (theta =90), and refrings2( theat not 90)
				refrings1= []
				refrings2= []
				sn = int(symmetry_string[1:])
				for i in xrange( len(refrings) ):
					if( sn%2 ==0 and abs( refrings[i].get_attr('n3') ) <1.0e-6 and (symmetry_string[0] == "c" or symmetry_string[0] =="d" ) ):
						refrings1.append( refrings[i])

					else:
						refrings2.append( refrings[i])
								
				del refrings
				for im in xrange( anima[nch] ):
					peak1 = None
					peak2 = None
					if ( len(refrings1) > 0):
						if  an[N_step] == -1:
							peak1, phihi1, theta1, psi1, sxi1, syi1, t11 = proj_ali_helical_90(adata[im],refrings1,numr,xrng[N_step],yrng[N_step],stepx[N_step],ynumber[N_step],psi_max,finfo,)
						else:
							peak1, phihi1, theta1, psi1, sxi1, syi1, t11 = proj_ali_helical_90_local(adata[im],refrings1,numr,xrng[N_step],yrng[N_step],stepx[N_step],ynumber[N_step], an[N_step], psi_max,  finfo,)
					if( len(refrings2) > 0):
						if  an[N_step] == -1:
							peak2, phihi2, theta2, psi2, sxi2, syi2, t12 = proj_ali_helical(adata[im],refrings2,numr,xrng[N_step],yrng[N_step],stepx[N_step],ynumber[N_step],psi_max,finfo,)
						else:
							peak2, phihi2, theta2, psi2, sxi2, syi2, t12 = proj_ali_helical_local(adata[im],refrings2,numr,xrng[N_step],yrng[N_step],stepx[N_step],ynumber[N_step], an[N_step], psi_max,finfo,)
					if peak1 is None: 
						peak = peak2
						phihi = phihi2
						theta = theta2
						psi = psi2
						sxi = sxi2
						syi = syi2
						t1 = t12
					elif peak2 is None:
						peak = peak1
						phihi = phihi1
						theta = theta1
						psi = psi1
						sxi = sxi1
						syi = syi1
						t1 = t11
					else:
						if(peak1 >= peak2):
							peak = peak1
							phihi = phihi1
							theta = theta1
							psi = psi1
							sxi = sxi1
							syi = syi1
							t1 = t11
						else:
							peak = peak2
							phihi = phihi2
							theta = theta2
							psi = psi2
							sxi = sxi2
							syi = syi2
							t1 = t12
					#peak, phihi, theta, psi, sxi, syi, t1 = proj_ali_helical(data[im],refrings,numr,xrng[N_step],yrng[N_step],stepx[N_step],ynumber[N_step],psi_max,finfo,)
					if(peak > -1.0e22):
						#Guozhi Tao: wrap y-shifts back into box within rise of one helical unit by changing phi
						jdelta=0.0
						dpp = (float(dp)/pixel_size)

						tp = Transform({"type":"spider","phi":phihi,"theta":theta,"psi":psi})
						tp.set_trans( Vec2f( -sxi, -syi ) )
						trans1 = tp.get_pre_trans()
						dyi = (float(trans1[2])/dpp)-int(trans1[2]/dpp)

						if dyi < -0.5-jdelta:  eyi = dyi+1.0
						elif dyi > 0.5+jdelta:  eyi = dyi-1.0
						else:                   eyi = dyi

						nperiod = float(eyi*dpp-trans1[2])/dpp
						th = Transform({"type":"spider","phi": -nperiod*dphi, "tz":nperiod*dpp})
						tfinal = tp*th
						ddd = tfinal.get_params("spider")
						sxnew = - ddd["tx"]
						synew = - ddd["ty"]
						phinew  = ddd["phi"]

						phihi = phinew
						sxi   = sxnew
						syi   = synew

						# unique range identified by [k0,k1], [k2,k3]
						tp = Transform({"type":"spider","phi":phihi,"theta":theta,"psi":psi})
						tp.set_trans( Vec2f( -sxi, -syi ) )
						k0 = 0.0
						k2 = k0+180

						if( abs( tp.at(2,2) )<1.0e-6 ):
							if (symmetry_string[0] =="c"):
								if sn%2 == 0:
									k1=360.0/sn
								else:
									k1=360.0/2/sn
							elif (symmetry_string[0] =="d"):
								if sn%2 == 0:
									k1=360.0/2/sn
								else:
									k1=360.0/4/sn
						else:
							k1=360.0/sn
						k3 = k1 +180
						from utilities import get_sym
						T = get_sym_sparx(symmetry_string[0:])

						for i in xrange( len(T) ):
							ttt = tp*Transform({"type":"spider","phi":T[i][0],"theta":T[i][1],"psi":T[i][2]})
							d1 = ttt.get_params("spider")

							if ( abs( tp.at(2,2) )<1.0e-6 ):
								if( sn%2==1 ): # theta=90 and n odd, only one of the two region match

									if( ( d1['phi'] <= float(k1) and d1['phi'] >= float(k0) ) or ( d1['phi'] < float(k3) and d1['phi'] >= float(k2) )):

										sxnew = - d1["tx"]
										synew = - d1["ty"]
										phinew = d1['phi']
										thetanew = d1["theta"]
										psinew = d1["psi"]
								else: #for theta=90 and n even, there is no mirror version during aligment, so only consider region [k0,k1]

									if( d1['phi'] <= float(k1) and d1['phi'] >= float(k0) ) :

										sxnew = - d1["tx"]
										synew = - d1["ty"]
										phinew = d1['phi']
										thetanew = d1["theta"]
										psinew = d1["psi"]

							else: #theta !=90, # if theta >90, put the projection into [k2,k3]. Otherwise put it into the region [k0,k1]
								if( sn==1):
									sxnew = sxi
									synew = syi
									phinew = phihi
									thetanew = theta
									psinew = psi
								else:

									if (tp.at(2,2) >0.0): #theta <90

										if( d1['phi'] <= float(k1) and d1['phi'] >= float(k0) ):
											if( cos( pi*float( d1['theta'] )/180.0 )>0.0 ):

												sxnew = - d1["tx"]
												synew = - d1["ty"]
												phinew = d1['phi']
												thetanew = d1["theta"]
												psinew = d1["psi"]

									else:
										if(  d1['phi'] <= float(k3) and d1['phi'] >= float(k2) ):
											if( cos( pi*float( d1['theta'] )/180.0 )<0.0 ):

												sxnew = - d1["tx"]
												synew = - d1["ty"]
												phinew = d1['phi']
												thetanew = d1["theta"]
												psinew = d1["psi"]

							del ttt,d1


						t2 = Transform({"type":"spider","phi":phinew,"theta":thetanew,"psi":psinew})
						t2.set_trans(Vec2f(-sxnew, -synew))
						adata[im].set_attr("xform.projection", t2)
						pixer[list1[im]]  = max_3D_pixel_error(t1, t2, numr[-3])
						modphi[list1[im]] = phinew
					else:
						# peak not found, parameters not modified
						pixer[list1[im]]  = 0.0
						phihi, theta, psi, sxi, syi = get_params_proj(adata[im])
						modphi[list1[im]] = phihi
				del adata, list1, list1_set
				del refrings1, refrings2
				
			#output pixel errors after all chunks are used
			from mpi import mpi_gatherv
			recvbuf = mpi_gatherv(pixer, nima, MPI_FLOAT, recvcount, disps, MPI_FLOAT, main_node, MPI_COMM_WORLD)
			del pixer		
			mpi_barrier(MPI_COMM_WORLD)
			terminate = 0
			if(myid == main_node):
				recvbuf = map(float, recvbuf)
				from utilities import write_text_file
				write_text_file([range(len(recvbuf)), recvbuf], os.path.join(outdir, "pixer_%04d_%04d.txt"%(N_step+1,Iter)) )
				from statistics import hist_list
				lhist = 20
				region, histo = hist_list(recvbuf, lhist)
				if(region[0] < 0.0):  region[0] = 0.0
				msg = "      Histogram of pixel errors\n      ERROR       number of particles\n"
				print_msg(msg)
				for lhx in xrange(lhist):
					msg = " %10.2f     %7d\n"%(region[lhx], histo[lhx])
					print_msg(msg)
				# Terminate if 95% within 1 pixel error
				im = 0
				for lhx in xrange(lhist):
					if(region[lhx] > 1.0): break
					im += histo[lhx]
				#if(im/float(total_nima) > 0.95):  terminate = 1
				del region, histo
			#output distribution of phi
			#jeanmod
			recvbuf = mpi_gatherv(modphi, nima, MPI_FLOAT, recvcount, disps, MPI_FLOAT, main_node, MPI_COMM_WORLD)
			#end jeanmod
			mpi_barrier(MPI_COMM_WORLD)
			del modphi
			if(myid == main_node):
				recvbuf = map(float, recvbuf)
				lhist = 30
				region, histo = hist_list(recvbuf, lhist)
				msg = "\n      Distribution of phi\n      phi         number of particles\n"
				print_msg(msg)
				for lhx in xrange(lhist):
					msg = " %10.2f     %7d\n"%(region[lhx], histo[lhx])
					print_msg(msg)
				del region, histo
			del recvbuf
			terminate = mpi_bcast(terminate, 1, MPI_INT, 0, MPI_COMM_WORLD)
			terminate = int(terminate[0])
			if myid == main_node:
				print_msg("Time to compute pixer = %d\n"%(time()-start_time))
				start_time = time()
			# write out headers, under MPI writing has to be done sequentially
			mpi_barrier(MPI_COMM_WORLD)
			m = 5
			from mpi import mpi_recv, mpi_send, MPI_TAG_UB, MPI_COMM_WORLD, MPI_FLOAT
			if myid == main_node:
				
				fexp = open(os.path.join(outdir, "parameters_%04d_%04d.txt"%(N_step+1,Iter)),"w")
				for n in xrange(number_of_proc):
					if n!=main_node:
						t = mpi_recv(recvcount[n]*m,MPI_FLOAT, n, MPI_TAG_UB, MPI_COMM_WORLD)
						for i in xrange(recvcount[n]):
							for j in xrange(m):
								fexp.write(" %15.5f  "%t[j+i*m])
							fexp.write("\n")
					else:
						t = [0.0]*m
						for i in xrange(recvcount[myid]):
							t = get_params_proj(data[i])
							for j in xrange(m):
								fexp.write(" %15.5f  "%t[j])
							fexp.write("\n")
				fexp.close()
				del t
	        	else:
				nvalue = [0.0]*m*recvcount[myid]
				t = [0.0]*m
				for i in xrange(recvcount[myid]):
					t = get_params_proj(data[i])
					for j in xrange(m):
						nvalue[j + i*m] = t[j]
				mpi_send(nvalue, recvcount[myid]*m, MPI_FLOAT, main_node, MPI_TAG_UB, MPI_COMM_WORLD)
				del nvalue
			if myid == main_node:
				print_msg("Time to write parameters = %d\n"%(time()-start_time))
				start_time = time()

			if CTF: vol = recons3d_4nn_ctf_MPI(myid, data, symmetry=sym, snr = snr, npad = npad, xysize = xysize)
			else:    vol = recons3d_4nn_MPI(myid, data, symmetry=sym, npad = npad, xysize = xysize)

			if myid == main_node:
				print_msg("\n3D reconstruction time = %d\n"%(time()-start_time))
				start_time = time()

			if fourvar:
			#  Compute Fourier variance
				for im in xrange(nima):
					original_data[im].set_attr( 'xform.projection', data[im].get_attr('xform.projection') )
				varf = varf3d_MPI(original_data, ssnr_text_file = os.path.join(outdir, "ssnr%04d"%(total_iter)), mask2D = None, reference_structure = vol, ou = last_ring, rw = 1.0, npad = 1, CTF = CTF, sign = 1, sym =sym, myid = myid)
				if myid == main_node:
					print_msg("Time to calculate 3D Fourier variance= %d\n"%(time()-start_time))
					start_time = time()
					varf = 1.0/varf
					drop_image(varf, os.path.join(outdir, "varf%04d.hdf"%(total_iter)))
			else:  varf = None

			#search for helical symmetry
			if myid == main_node:
				drop_image(vol, os.path.join(outdir, "vol%04d.hdf"%(total_iter)))
							
			if(total_iter > nise):
				bcast_EMData_to_all(vol, myid, main_node)
				#from filter import filt_gaussl
				#vol = filt_gaussl(vol, 0.25)

				if myid == main_node:
					lprms = []
					for i in xrange(-ndp,ndp+1,1):
						for j in xrange(-ndphi,ndphi+1,1):
							lprms.append( dp   + i*dp_step)
							lprms.append( dphi + j*dphi_step)
					#print "lprms===",lprms
					recvpara = []
					for im in xrange(number_of_proc):
						helic_ib, helic_ie = MPI_start_end(nlprms, number_of_proc, im)
						recvpara.append(helic_ib )
						recvpara.append(helic_ie )

				para_start, para_end = MPI_start_end(nlprms, number_of_proc, myid)

				list_dps     = [0.0]*((para_end-para_start)*2)
				list_fvalues = [-1.0]*((para_end-para_start)*1)

				if myid == main_node:
					for n in xrange(number_of_proc):
						if n!=main_node: mpi_send(lprms[2*recvpara[2*n]:2*recvpara[2*n+1]], 2*(recvpara[2*n+1]-recvpara[2*n]), MPI_FLOAT, n, MPI_TAG_UB, MPI_COMM_WORLD)
						else:    list_dps = lprms[2*recvpara[2*0]:2*recvpara[2*0+1]]
				else:
					list_dps = mpi_recv((para_end-para_start)*2, MPI_FLOAT, main_node, MPI_TAG_UB, MPI_COMM_WORLD)

				list_dps = map(float, list_dps)

				local_pos = [0.0, 0.0, -1.0e20]
				for i in xrange(para_end-para_start):
					fvalue = helios7(vol, pixel_size, list_dps[i*2], list_dps[i*2+1], fract, rmax, rmin)
					if(fvalue >= local_pos[2]):
						local_pos = [list_dps[i*2], list_dps[i*2+1], fvalue ]
				if myid == main_node:
					list_return = [0.0]*(3*number_of_proc)
					for n in xrange(number_of_proc):
						if n != main_node: list_return[3*n:3*n+3]                 = mpi_recv(3,MPI_FLOAT, n, MPI_TAG_UB, MPI_COMM_WORLD)
 						else:              list_return[3*main_node:3*main_node+3]  = local_pos[:]
				else:
					mpi_send(local_pos, 3, MPI_FLOAT, main_node, MPI_TAG_UB, MPI_COMM_WORLD)

				if myid == main_node:	
					maxvalue = list_return[2]
					for i in xrange(number_of_proc):
						if( list_return[i*3+2] >= maxvalue ):
							maxvalue = list_return[i*3+2]
							dp       = list_return[i*3+0]
							dphi     = list_return[i*3+1]
					dp   = float(dp)
					dphi = float(dphi)
					#print  "  GOT dp dphi",dp,dphi

					vol  = vol.helicise(pixel_size,dp, dphi, fract, rmax, rmin)
					print_msg("New delta z and delta phi      : %s,    %s\n\n"%(dp,dphi))		
			
			else:
				if myid==main_node:
					#  in the first nise steps the symmetry is imposed
					vol = vol.helicise(pixel_size,dp, dphi, fract, rmax, rmin)
					print_msg("Imposed delta z and delta phi      : %s,    %s\n\n"%(dp,dphi))
			if(myid==main_node):
				fofo = open(os.path.join(outdir,datasym),'a')
				fofo.write('  %12.4f   %12.4f\n'%(dp,dphi))
				fofo.close()
				ref_data = [vol]
				if  fourvar:  ref_data.append(varf)
				vol = user_func(ref_data)
				vol = vol.helicise(pixel_size,dp, dphi, fract, rmax, rmin)

				drop_image(vol, os.path.join(outdir, "volf%04d.hdf"%(total_iter)))
				print_msg("\nSymmetry search and user function time = %d\n"%(time()-start_time))
				start_time = time()

			bcast_EMData_to_all(vol, myid, main_node)
			dp   = bcast_number_to_all(dp,   source_node = main_node)
			dphi = bcast_number_to_all(dphi, source_node = main_node)
			#
			del varf
	par_str = ["xform.projection"]
	if myid == main_node:
	   	if(file_type(stack) == "bdb"):
	        	from utilities import recv_attr_dict_bdb
	        	recv_attr_dict_bdb(main_node, stack, data, par_str, image_start, image_end, number_of_proc)
	        else:
	        	from utilities import recv_attr_dict
	        	recv_attr_dict(main_node, stack, data, par_str, image_start, image_end, number_of_proc)
		print_msg("Time to write header information= %d\n"%(time()-start_time))
		start_time = time()
	else:	       send_attr_dict(main_node, data, par_str, image_start, image_end)
	if myid == main_node: print_end_msg("ihrsr_MPI")

def generate_sections(volstkpref,stack,miclistfile,fname_params, xysize, zsize,pixel_size,dphi,dp,rmax,rmin=0,section_use=-1,fract= 0.67,Nuse=-1,nsegthr=3):
	from utilities import read_text_row, get_im, set_params_proj, model_blank
	from reconstruction import recons3d_4nn_ctf
	from filter import filt_tanl, filt_gaussl
	from math import fmod
	from mpi            import mpi_comm_size, mpi_comm_rank, MPI_COMM_WORLD, mpi_barrier

	nproc = mpi_comm_size(MPI_COMM_WORLD)
	myid = mpi_comm_rank(MPI_COMM_WORLD)
	main_node=0
	mpi_barrier(MPI_COMM_WORLD)

	params=read_text_row(fname_params)
	miclist = read_text_row(miclistfile)
	N = len(miclist)
	volstk = volstkpref+'%i.hdf'%myid
	winxy = int(rmax)*2 + 4
	dpp = (float(dp)/pixel_size) 

	if section_use < 0:
		if fmod(dpp,1) == 0: 
			if (int(dpp))%2 == 0:
				section_use = int(2*dpp)
			else:
				section_use = 2*(1+int(dpp/2)) + int(dpp)
		else:
			section_use = 2*(1+int(dpp/2)) + int(dpp)+1
	a=get_im(stack,0)
	nx = a.get_xsize()
	ny = a.get_ysize()
	counter = 0
	if Nuse > 0:
		N = Nuse
	for ifil in xrange(N):
		if ifil%nproc != myid:
			continue

		mic = miclist[ifil][6:]
		mic = map(int, mic)
		if len(mic) < nsegthr:
			continue
		data=EMData.read_images(stack, mic)
		for i in xrange(len(mic)):
			set_params_proj(data[i],[params[mic[i]][0], params[mic[i]][1],params[mic[i]][2],params[mic[i]][3],params[mic[i]][4]])

		fvol = recons3d_4nn_ctf(data, list_proj=[], snr=1.0, sign=1, symmetry='c1', verbose=0, npad=2, xysize=xysize, zsize=zsize)
		fvol = fvol.helicise(pixel_size, dp, dphi, fract, rmax, rmin)
		fvol.set_attr('filament', ifil)
		fvol.set_attr('nfil', len(mic))
		# window out relevant part
		fvol = Util.window(fvol,winxy,winxy,section_use)
		fvol.write_image(volstk, counter)
		counter += 1

def avgvol_helical(volname, stack, fname_params,miclistfile,pixel_size,dphi,dp,nx, ny,nz ,rmax,rmin=0.0,fract= 0.67, nsegthr=3):
	from mpi            import mpi_comm_size, mpi_comm_rank, MPI_COMM_WORLD, mpi_barrier
	from utilities      import reduce_EMData_to_root,set_params_proj, read_text_row, model_blank
	from reconstruction import recons3d_4nn_ctf
	nproc = mpi_comm_size(MPI_COMM_WORLD)
	myid = mpi_comm_rank(MPI_COMM_WORLD)
	main_node=0
	mpi_barrier(MPI_COMM_WORLD)

	miclist = read_text_row(miclistfile)
	N = len(miclist)

	counter = 0
	avgvol = model_blank(nx,ny,nz)
	params=read_text_row(fname_params)
	for ifil in xrange(N):
		if ifil%nproc != myid:
			continue

		mic = miclist[ifil][6:]
		mic = map(int, mic)
		if len(mic) < nsegthr:
			continue
		data = EMData.read_images(stack,mic)
		for i in xrange(len(mic)):
			set_params_proj(data[i],[params[mic[i]][0], params[mic[i]][1],params[mic[i]][2],params[mic[i]][3],params[mic[i]][4]])

		fvol = recons3d_4nn_ctf(data, list_proj=[], snr=1.0, sign=1, symmetry='c1', verbose=0, npad=2, xysize=nx, zsize=nz)
		fvol = fvol.helicise(pixel_size, dp, dphi, fract, rmax, rmin)
		Util.add_img(avgvol, fvol)
		counter = counter + 1
		print counter
	Util.mul_scalar(avgvol, 1.0/counter)
	reduce_EMData_to_root(avgvol, myid)
	if myid == main_node:
		Util.mul_scalar(avgvol, 1.0/nproc)
		avgvol.write_image(volname)

def get_nfil_segs(fmiclist, fnewmiclist, fstack, fnewstack, nrows):
	from utilities import read_text_row, write_text_row, get_im
	miclist = read_text_row(fmiclist)
	newmiclist = miclist[0:nrows]
	write_text_row(newmiclist,fnewmiclist)
	
	segslist = []
	for ifil in xrange(nrows):
		mic = miclist[ifil][6:]
		mic = map(int, mic)
		segslist += mic
	counter = 0
	for iseg in segslist:
		a = get_im(fstack, iseg)
		a.write_image(fnewstack, counter)
		counter += 1

def get_params_by_thr(fmiclist, fstack='', fnewstack='', fparams='', fnewparams='', segthr=3):
	from utilities import read_text_row, get_im, write_text_row
	
	'''

	fstack: file name of existing image stack
	fnewstack: file name of new stack. If fstack is not empty, then extract from fstack those segments that belong to filaments with segthr or more segments and write images to fnewstack.
	
	fparams: file name of parameters
	fnewparams: if fparams is not empty, then extract from fparams those parameters of segments that belong to filaments with segthr or more segments and write images to fnewparams.
	
	'''
	miclist = read_text_row(fmiclist)
	segslist = []
	nm = len(miclist)
	for ifil in xrange(nm):
		mic = miclist[ifil][6:]
		mic = map(int, mic)
		if len(mic) >= segthr:
			segslist += mic
	
	if len(fparams) > 0:
		params = read_text_row(fparams)
		newparams=[]
		for iseg in segslist:
			newparams.append(params[iseg])

		write_text_row(newparams, fnewparams)
	
	if len(fstack) > 0:
		counter = 0
		for iseg in segslist:
			a = get_im(fstack, iseg)
			a.write_image(fnewstack, counter)
			counter += 1


def reconsvolhelical(stack, fparams, newvol, xysize,zsize, dp, dphi, pixel_size, rmax, rmin=0,fract = 0.67,sym='c1'):
	from utilities import read_text_row, set_params_proj
	from reconstruction import recons3d_4nn_ctf
	from os import system
	
	params=read_text_row(fparams)
	nima = EMUtil.get_image_count(stack)
	data = EMData.read_images(stack, range(nima))
	
	for i in xrange(nima):
		set_params_proj(data[i], [params[i][0], params[i][1],params[i][2],params[i][3],params[i][4]])
	
	fvol = recons3d_4nn_ctf(data, list_proj=[], snr=1.0, sign=1, symmetry=sym, verbose=0, npad=2, xysize=xysize, zsize=zsize)
	fvol = fvol.helicise(pixel_size, dp, dphi, fract, rmax, rmin)
	fvol.write_image(newvol)

def makefilstack(stack, newstack, nfil, miclist, thr=3):
	from utilities import read_text_row, get_im
	miclist = read_text_row(miclist)
	mic= miclist[nfil][6:]
	if len(mic) <thr:
		print "number of segments in filament less than 3"
		print "will not generate stack"
		sys.exit()
	else:
		print len(mic), "segments"
	mic = map(int,mic)
	nima = len(mic)
	for i in xrange(nima):
		ind = mic[i]
		a = get_im(stack, ind)
		a.set_attr('ID_%s'%stack, ind)
		a.write_image(newstack, i)

def helicalrecons3D_MPI(nx,ny,nz, stack, miclist, dp, pixel_size, dphi,avgvol, user_func_name, rmax, rmin=0,CTF=False, fract=0.67,outdir='.'):
	from mpi              import mpi_bcast, mpi_comm_size, mpi_comm_rank, MPI_COMM_WORLD, mpi_barrier, MPI_INT, MPI_TAG_UB, MPI_FLOAT, mpi_recv, mpi_send
	from utilities        import get_params_proj, read_text_row, model_cylinder,pad, set_params3D, get_params3D, reduce_EMData_to_root, bcast_EMData_to_all, bcast_number_to_all, drop_image, bcast_EMData_to_all, model_blank
	from utilities        import send_attr_dict, file_type
	from fundamentals     import resample, rot_shift3D
	from applications     import MPI_start_end
	from math             import fmod, atan, pi
	from utilities        import model_blank
	from filter           import filt_tanl, filt_ctf
	import os
	from statistics       import fsc_mask
	from copy             import copy
	from os               import sys
	from time import time	
	from reconstruction import recons3d_4nn_ctf
	myid  = mpi_comm_rank(MPI_COMM_WORLD)
	nproc = mpi_comm_size(MPI_COMM_WORLD)
	main_node = 0
	mpi_barrier(MPI_COMM_WORLD)
	miclist = read_text_row(miclist)
	if(myid == main_node):
		fil_list = []
		nm = len(miclist)
		for ifil in xrange(nm):
			mic = miclist[ifil][6:]
			mic = map(int, mic)
			fil_list.append(ifil)
		nvols = len(fil_list)
	else:
		nvols = 0

	total_nvols = bcast_number_to_all(nvols, source_node = main_node)
	
	if myid != main_node:
		fil_list = [-1]*total_nvols
	fil_list = mpi_bcast(fil_list, total_nvols, MPI_INT, 0, MPI_COMM_WORLD)
	fil_list = map(int, fil_list)
	
	dpp = (float(dp)/pixel_size)
	fil_start, fil_end = MPI_start_end(total_nvols, nproc, myid)
	fil_list = fil_list[fil_start:fil_end]
	nvols = len(fil_list)
	# inds[i] =[start, end] for the segments in i-th filament
	segslist = []
	inds = [[] for i in xrange(nvols)]
	seg_start = None
	seg_end = None
	for i in xrange(nvols):
		ifil = fil_list[i]
		mic = miclist[ifil][6:]
		mic = map(int, mic)
		if i == 0:
			seg_start = mic[0]
		if i == nvols - 1:
			seg_end = mic[len(mic)-1]
		segslist += mic
		start = 0
		if i > 0:
			start = inds[i-1][1]+1
		inds[i] = [start, start+len(mic) - 1]
	if seg_start != segslist[0] or seg_end != segslist[len(segslist)-1]:
		print "something wrong with segment index"
		sys.exit()
	data=EMData.read_images(stack, segslist)
	nima = len(data)
	volsum = model_blank(nx,ny,nz)
	for ivol in xrange(nvols):
		#print ivol, nvols,myid
		if CTF:
			vol = recons3d_4nn_ctf(data[inds[ivol][0]:inds[ivol][1]], list_proj=[], snr=1.0, sign=1, symmetry='c1', verbose=0, npad=2, xysize=nx, zsize=nz)
		Util.add_img(volsum, vol)
	# see var_MPI in applications.py
	reduce_EMData_to_root(volsum, myid)
	if myid == main_node:
		import user_functions
		user_func = user_functions.factory[user_func_name]
		# calculate average
		Util.mul_scalar(volsum, 1.0 /float(total_nvols) )
		volsum = volsum.helicise(pixel_size,dp, dphi, fract, rmax, rmin)
		ref_data = [volsum]
		volf = user_func(ref_data)
		volf = volf.helicise(pixel_size,dp, dphi, fract, rmax, rmin)
		drop_image(volf, os.path.join(outdir, avgvol))


def diskaliD_MPI(stack, ref_vol, outdir, maskfile, dp, dphi, pixel_size, user_func_name, zstep=1.0, fract=0.67, rmax=70, rmin=0, \
					 CTF=False, maxit=1, sym = "d1"):
	#  This version is for D symmetries, which have to be handled differently from C.
	from mpi              import mpi_bcast, mpi_comm_size, mpi_comm_rank, MPI_COMM_WORLD, mpi_barrier, MPI_INT, MPI_TAG_UB, MPI_FLOAT, mpi_recv, mpi_send
	from utilities        import get_params_proj, read_text_row, model_cylinder,pad, set_params3D, get_params3D, reduce_EMData_to_root, bcast_EMData_to_all, bcast_number_to_all, drop_image, bcast_EMData_to_all, model_blank
	from utilities        import send_attr_dict, file_type, sym_vol, get_image
	from fundamentals     import resample, rot_shift3D
	from applications     import MPI_start_end, cylindrical_trans, alihelical3, stack_disks, ordersegments
	from math             import fmod, atan, pi
	from utilities        import model_blank
	from filter           import filt_tanl, filt_ctf
	import os
	from statistics       import fsc_mask
	from copy             import copy
	from os               import sys
	from time             import time	
	from alignment        import Numrinit, ringwe
	from reconstruction   import recons3d_wbp
	from morphology       import ctf_2
	import types
	
	myid  = mpi_comm_rank(MPI_COMM_WORLD)
	nproc = mpi_comm_size(MPI_COMM_WORLD)
	main_node = 0

	mpi_barrier(MPI_COMM_WORLD)

	if myid == main_node:
		if not(os.path.exists(outdir)):
			os.mkdir(outdir)
	mpi_barrier(MPI_COMM_WORLD)

	dpp = (float(dp)/pixel_size)
	rise = int(dpp)
	winxy = int(rmax)*2 + 4
	if(float(rise) != dpp):
		print "  dpp has to be integer multiplicity of the pixel size"
		sys.exit()
	rise3 = 3*rise
	# for resampling to polar rmin>1
	rminpolar = max(1,rmin)

	# For D symmetry we will use the fact that Dn symmetry is a composition of D1 and Cn.
	#  implicit symmetrizations are done using C
	sym = "C"+sym[1:]

	from utilities import get_im
	refvol = get_im(ref_vol)
	ref_nx = refvol.get_xsize()
	ref_ny = refvol.get_ysize()
	ref_nz = refvol.get_zsize()
	if(ref_nz < rise3):
		print  "  reference volumes has to be at least 3*rise long "
		sys.exit()
	
	if maskfile:
		if type(maskfile) is types.StringType: mask3D = get_image(maskfile)
		else:                                  mask3D = maskfile
	else:
		if rmin > 0:
			mask3D = model_cylinder(rmax, ref_nx, ref_ny, ref_nz) - model_cylinder(rmin, ref_nx, ref_ny, ref_nz)
		else:
			mask3D = model_cylinder(rmax, ref_nx, ref_ny, ref_nz)

	filaments = ordersegments(stack)
	
	total_nfils = len(filaments)
	fstart, fend = MPI_start_end(total_nfils, nproc, myid)
	filaments = filaments[fstart:fend]
	nfils = len(filaments)

	list_of_particles = []
	indcs = []
	k = 0
	for i in xrange(nfils):
		list_of_particles += filaments[i]
		k1 = k+len(filaments[i])
		indcs.append([k,k1])
		k = k1

	if( myid == 0 ):
		# Build a 3D dedicated correction function
		if CTF:
			from morphology import ctf_2
			cc = EMUtil.get_all_attributes(stack, 'ctf')
			ctf2 = ctf_2(ref_nz, cc[0])
			ncc = len(ctf2)
			for i in xrange(1,len(cc)):
				temp = ctf_2(ref_nz, cc[i])
				for k in xrange(ncc): ctf2[k] += temp[k]
			del temp
		from math import sqrt
		rrc = model_blank(ref_nz, ref_nz, ref_nz)
		rc = ref_nz//2
		for i in xrange(ref_nz):
			ic = (i-rc)**2
			for j in xrange(ref_nz):
				jc = (j-rc)**2
				dc = sqrt(ic+jc)
				for k in xrange(ref_nz):
					if CTF:
						rr = sqrt((k-rc)**2 + ic + jc)
						rin = int(rr)
						drin = rr-rin
						rrc.set_value_at(i,j,k, dc/(1.0+(1.0-drin)*ctf2[rin] + drin*ctf2[rin+1]) )
					else:
						rrc.set_value_at(i,j,k, dc )

	data = EMData.read_images(stack, list_of_particles)
	nima = len(data)

	data_nx = data[0].get_xsize()
	data_ny = data[0].get_xsize()
	mask2D  = pad(model_blank(2*int(rmax), data_ny, 1, 1.0), data_nx, data_ny, 1, 0.0)
	for im in xrange(nima):
		data[im].set_attr('ID', list_of_particles[im])
		if CTF:
			ctf_params = data[im].get_attr("ctf")
			st = Util.infomask(data[im], mask2D, False)
			data[im] -= st[0]
			data[im] = filt_ctf(data[im], ctf_params)
			data[im].set_attr('ctf_applied', 1)
	del list_of_particles, mask2D

	# This is the part sufficient to get polar resampling
	refvol = Util.window(refvol, winxy, winxy, rise)
	rr = ref_nz//2-2
	# do full sized reconstruction with the projection parameters BEFORE they are modified by disk alignment step
	fullvolsum0 = model_blank(ref_nx, ref_ny, ref_nz)
	filvols = []
	data_slices = []
	start = time()
	Torg = [None]*len(data)
	for i in xrange(len(data)):  Torg[i] = data[i].get_attr("xform.projection")
	msk = model_cylinder(rmax-1, ref_nx, ref_ny, ref_nz) - model_cylinder(rmax-2, ref_nx, ref_ny, ref_nz)
	ms3 = model_cylinder(rmax, ref_nx, ref_ny, ref_nz)
	for ivol in xrange(nfils):
		fullvol0 = Util.window(recons3d_wbp(data, list_proj=range(indcs[ivol][0],indcs[ivol][1]), method = None, symmetry=sym, radius=rr), ref_nx, ref_ny, ref_nz, 0, 0, 0)
		fullvol0 = fullvol0.helicise(pixel_size, dp, dphi, fract, rmax, rmin)
		fullvol0 = sym_vol(fullvol0, symmetry=sym)
		stat = Util.infomask(fullvol0, msk, True)
		fullvol0 -= stat[0]
		fullvol0 *= ms3
		"""
		from filter import filt_tanl
		fullvol0 = filt_tanl(fullvol0, 0.3, 0.2)
		fullvol0 = fullvol0.helicise(pixel_size, dp, dphi, fract, rmax, rmin)
		fullvol0 = sym_vol(fullvol0, symmetry=sym)
		"""
		#  Reference volume is fully symmetrized
		Util.add_img(fullvolsum0, sym_vol(fullvol0, "D1"))

		# to allow phi search, 3disks are only C-symmetrized
		filvols.append(Util.window(fullvol0, winxy, winxy, rise3))

		data_slices.append(cylindrical_trans(filvols[ivol], rminpolar, rmax, rise, True))
		#print  " DONE ", myid,ivol
	#  original projection data no longer needed
	del data, fullvol0, ms3, msk
	if myid == main_node:
		tt = time()
		print " TIME to do initial reconstructions :", tt-start
		start = tt
	reduce_EMData_to_root(fullvolsum0, myid)
	if myid == main_node and False:
		fullvolsum0.write_image("verify0.hdf")
		rrc.write_image("verify1.hdf")
		fullvolsum0 = Util.window( pad(fullvolsum0, ref_nz, ref_nz, ref_nz, 0.0).filter_by_image(rrc), ref_nx, ref_ny, ref_nz, 0, 0, 0)
		fullvolsum0 = fullvolsum0.helicise(pixel_size, dp, dphi, fract, rmax, rmin)
		fullvolsum0 = sym_vol(fullvolsum0, symmetry=sym)
		stat = Util.infomask(fullvolsum0,  model_cylinder(rmax-1, ref_nx, ref_ny, ref_nz) - model_cylinder(rmax-2, ref_nx, ref_ny, ref_nz), True)
		fullvolsum0 -= stat[0]
		fullvolsum0 *= model_cylinder(rmax-1, ref_nx, ref_ny, ref_nz)
		fullvolsum0.write_image("verify2.hdf")
		#for i in xrange(len(filvols)):
		#	filvols[i].write_image("verify7.hdf",i)
	mpi_barrier(MPI_COMM_WORLD)
	#sys.exit()

	T_filament = [None]*nfils
	for iter in xrange(1,maxit+1):
		if myid == main_node: print " ITERATION:  ", iter
		refslices = cylindrical_trans(refvol, rminpolar, rmax, rise)

		for ivol in xrange(nfils):
			#if myid == main_node: print "  ALI  ",myid,ivol,Util.infomask(data_slices[ivol][0], None, True)
			T_filament[ivol] = alihelical3(data_slices[ivol], refslices, zstep, dphi, rise, rminpolar, rmax, "D"+sym[1:])

		refvol = model_blank(winxy, winxy, rise3)
		for ivol in xrange(nfils):
			d = T_filament[ivol].get_params('spider')
			#if myid == main_node: print  d["phi"], d["theta"], d["psi"], d["tz"]
			#if myid == main_node:
			#	pad(rot_shift3D(filvols[ivol],  d["phi"], d["theta"], d["psi"], sz=d["tz"]), ref_nx, ref_ny, ref_nz, 0.0).write_image("verify1.hdf")
			Util.add_img(refvol, rot_shift3D(filvols[ivol],  d["phi"], d["theta"], d["psi"], sz=d["tz"]))
		reduce_EMData_to_root(refvol, myid)

		if myid == main_node:
			#refvol.write_image("verify3.hdf")
			#Util.window(refvol, winxy, winxy, rise).write_image("verify4.hdf")
			refvol =  stack_disks(Util.window(refvol, winxy, winxy, rise), winxy, winxy, ref_nz, dphi, rise)
			#refvol.write_image("verify5.hdf")
			refvol = sym_vol(sym_vol(refvol, symmetry=sym), "D1")
			#refvol.write_image("verify6.hdf")
			if CTF:
				refvol = Util.window( pad(refvol, ref_nz, ref_nz, ref_nz, 0.0).filter_by_image(rrc), ref_nx, ref_ny, ref_nz, 0, 0, 0)
				refvol = refvol.helicise(pixel_size, dp, dphi, fract, rmax, rmin)
				refvol = sym_vol(sym_vol(refvol, symmetry=sym), "D1")
				stat = Util.infomask(refvol,  model_cylinder(rmax-1, ref_nx, ref_ny, ref_nz) - model_cylinder(rmax-2, ref_nx, ref_ny, ref_nz), True)
				refvol -= stat[0]
				refvol *= model_cylinder(rmax-1, ref_nx, ref_ny, ref_nz)

			import user_functions
			user_func = user_functions.factory[user_func_name]

			ref_data = [refvol, mask3D]
			refvol = user_func(ref_data)
			refvol = refvol.helicise( pixel_size , dp, dphi, fract, rmax, rmin)
			refvol = sym_vol(sym_vol(refvol, symmetry=sym), "D1")
			if(iter == maxit):  drop_image(refvol, os.path.join(outdir, 'fvheli.hdf'))
			refvol = Util.window(refvol, winxy, winxy, rise)

			tt = time()
			print " TIME of one iteration :", tt-start
			start = tt
		else:
			refvol = model_blank(winxy, winxy, rise)

		bcast_EMData_to_all(refvol, myid, main_node)

	del refvol, data_slices

	# update rotations of individual images
	data = [EMData() for i in xrange(nima)]
	forg = []
	helisym = Transform({"type":"spider","phi":dphi,"tz":dpp})
	ihelisym = helisym.inverse()
	from utilities import get_params_proj, set_params_proj
	permitrange = rise/2.0
	
	for ivol in xrange(nfils):
		#  This is for printout
		d = T_filament[ivol].get_params('spider')
		#print  d["phi"],d["theta"],d["psi"],d["tx"],d["ty"],d["tz"]
		if( d["tz"] < 0.0 ):  d = (helisym*T_filament[ivol]).get_params('spider')
		#print  d["phi"],d["theta"],d["psi"],d["tx"],d["ty"],d["tz"]
		#forg.append([d["phi"],d["theta"],d["psi"],d["tx"],d["ty"],d["tz"]])
		#  Use inverse transformation to modify projection directions.
		Tinv = T_filament[ivol].inverse()
		for im in xrange(indcs[ivol][0],indcs[ivol][1]):
			#d = Torg[im].get_params('spider')
			#print  d["phi"],d["theta"],d["psi"],d["tx"],d["ty"],d["tz"]
			Torg[im] = Torg[im]*Tinv
			d = Torg[im].get_params('spider')
			d["tx"] = -d["tx"]
			d["ty"] = -d["ty"]
			#if myid == 0: print  " A ",d["phi"],d["theta"],d["psi"],d["tx"],d["ty"]
			#finfo.write(" %d  %f  %f7.2  %f  %f  %f\n" %(im,d["phi"],d["theta"],d["psi"],d["tx"],d["ty"]))
			#finfo.flush()
			psi = d["psi"]
			while( d["ty"] < -permitrange ):
				if( psi > 180.0):  Torg[im] = Torg[im]*helisym
				else:              Torg[im] = Torg[im]*ihelisym
				d = Torg[im].get_params('spider')
				d["tx"] = -d["tx"]
				d["ty"] = -d["ty"]
			while( d["ty"] > permitrange ):
				if( psi > 180.0):  Torg[im] = Torg[im]*ihelisym
				else:              Torg[im] = Torg[im]*helisym
				d = Torg[im].get_params('spider')
				d["tx"] = -d["tx"]
				d["ty"] = -d["ty"]
			#finfo.write(" %d  %f  %f7.2  %f  %f  %f\n" %(im,d["phi"],d["theta"],d["psi"],d["tx"],d["ty"]))
			#finfo.flush()

			#if myid == 0: print  " B ",d["phi"],d["theta"],d["psi"],d["tx"],d["ty"]
			set_params_proj(data[im], [d["phi"],d["theta"],d["psi"],d["tx"],d["ty"]])
			data[im].set_attr('ID', filaments[ivol][im-indcs[ivol][0]])
	#from utilities import write_text_row
	#write_text_row(forg,"forg%03d.txt"%myid)
	"""
	forg = []
	for ivol in xrange(nfils):
		for im in xrange(inds[ivol][0],inds[ivol][1]):
				T = data[im].get_attr("xform.projection")
				d = T.get_params('spider')
				forg.append([d["phi"],d["theta"],d["psi"],-d["tx"],-d["ty"],data[im].get_attr('ID')])
	write_text_row(forg,"pheader%03d.txt"%myid)
	"""
	mpi_barrier(MPI_COMM_WORLD)

	# write out headers, under MPI writing has to be done sequentially
	from mpi import mpi_recv, mpi_send, MPI_TAG_UB, MPI_COMM_WORLD, MPI_FLOAT
	par_str = ['xform.projection', 'ID']
	if myid == main_node:
		if(file_type(stack) == "bdb"):
			from utilities import recv_attr_dict_bdb
			recv_attr_dict_bdb(main_node, stack, data, par_str, 0, nima, nproc)
		else:
			from utilities import recv_attr_dict
			recv_attr_dict(main_node, stack, data, par_str, 0, nima, nproc)
	else:
		send_attr_dict(main_node, data, par_str, 0, nima)


"""
#  This is a python version of what I put on 10/09/2012 in emdata_sparx.cpp PAP.
def heli(vol,  pixel_size=1.82, dp=27.6, dphi=166.5, section_use=0.67, radius=70., minrad=0.):
	from utilities import model_blank
	nx = vol.get_xsize()
	ny = vol.get_ysize()
	nz = vol.get_zsize()
	if(int(section_use*nz+0.5)>=nz-2):
		print "too large section use!"
		return  vol.copy()
	hel = model_blank(nx,ny,nz)
	nxc = nx//2
	nyc = ny//2
	nzc = nz//2
	volen = nz*pixel_size
	nzcp  = nzc*pixel_size
	sectl = nz*pixel_size*section_use
	nb = nzcp - sectl/2.0
	ne = nzcp + sectl/2.0
	numst = int( nz*pixel_size/dp )
	numri = int(sectl/dp)
	print  nz, nz*pixel_size, nzcp, sectl, nb, ne, numst, numri

	if(radius < 0.0): r2 = (nxc-1)**2
	else: r2 = radius*radius
	if(minrad < 0.0): ir = 0.0
	else: ir = minrad*minrad
	for k in xrange(nz):
		nq = 0
		for ist in xrange(numst):
			z = k*pixel_size + ist*dp
			phi = ist*dphi
			if( z >= volen ):
				z = k*pixel_size + (ist-numst)*dp
				phi = (ist-numst)*dphi 
			ca = cos(phi *pi/180.0);#  in python in radians!
			sa = sin(phi *pi/180.0);
			print k,ist,z,nb,ne
			if(z >= nb and z <= ne ):
				nq += 1
				if( nq > numri ): break
				print k,z
				zz = z/pixel_size
				for j in xrange(ny):
					jy = j - nyc;
					jj = jy*jy;
					for i in xrange(nx):
						ix = i - nxc;
						d2 = float((ix*ix + jj))
						if(d2 <= r2 and d2>=ir):
							xx =  ix*ca + jy*sa + nxc;
							yy = -ix*sa + jy*ca + nyc;
							#hel[i,j,k] += vol[int(xx+0.5), int(yy+0.5), int(zz+0.5)]

		if(nq <numri):
			print "incorrect number of repeats",numri,nq,k
			break
	#hel = binarize(hel,0.5)
	return hel
	

dp         = 27.6
pixel_size = 1.84
dphi       = 166.5
rmax       = 70
rmin       = 0
dpp = (float(dp)/pixel_size)

refvol = get_im('bigstackrecons_filt_elmar.hdf')

#refvol = Util.window(refvol,160,160,110)

o = heli(refvol, section_use=0.95)
o.write_image("heli.hdf")
exit()

"""

def make_test_stack(vol, stack, proj_stack, dphi, dp, pixel_size, fil_attr='filament', N = -1, randfile=None, phizero=False):

	'''
	INPUT
	
	vol: Name of volume from which projections are calculated
	stack: Name of data stack containing filament and ctf information
	proj_stack: Name of projection stack that will be generated from predicted parameters
	dphi: azimuthal rotation of helical symmetry
	dp: rise of helical symmetry
	pixel_size: pixel size
	fil_attr: attribute under which filament identifier is stored in header.
	randfile: name of text file to save the applied random rotations/shifts to
	
	OUTPUT
	
	Writes to disk a stack of projections named proj_stack calculated from vol using projection parameters predicted based on helical symmetry parameters.
	'''
	
	from random         import uniform, random
	from projection     import prep_vol, prgs
	from applications   import get_dist, ordersegments
	from filter			import filt_ctf
	from utilities		import get_im, write_text_row
	
	
	filaments = ordersegments(stack)
	
	ptclcoords = EMUtil.get_all_attributes(stack, 'ptcl_source_coord')

	if( N < 1 ): N = len(filaments)	
	volft, kbx, kby, kbz = prep_vol( get_im(vol) )

	THETA = 90.
	PSI   = 90.
	SX    = 0.

	rise = dp/pixel_size  # in pixels
	permitrange = rise/2.0
	ima = EMData()
	counter = 0
	randphiy=[]
	for ifil in xrange(N):
		mic = filaments[ifil]
		nsegs = len(mic)

		# This is the random rotation/shift to be applied to vol
		y0   = permitrange*random()
		if phizero:  phi0=0
		else:        phi0 = 360.0*random()
			
		rota = Transform({"type":"spider","phi":phi0,"theta":0,"psi":0,"tx":0,"ty":0,"tz":y0})

		# By multiplying projection direction of each segment by inverse of rota, the resulting reconstruction should be the 
		# refvol rotated/shifted by the random amout phi0/y0
		irota = rota.inverse()

		refcoords = ptclcoords[mic[0]]

		for iseg in xrange(nsegs):
			dA = pixel_size * get_dist(ptclcoords[mic[iseg]], refcoords) # distance in Angstroms between segments mic[iseg] and mic[0]

			iphi = (dA/dp*dphi)%360.0
			iy   = (dA%dp)/pixel_size # in pixels
			if( iy > permitrange ): iy -= rise

			ima.read_image(stack, mic[iseg], True)
			ct   = ima.get_attr('ctf')
			fifi = ima.get_attr(fil_attr)
			coord = ima.get_attr('ptcl_source_coord')
			tmp = filt_ctf(prgs(volft, kbz, [iphi, THETA,PSI,SX, iy], kbx, kby),ct)
			tmp.set_attr(fil_attr, fifi)

			t = tmp.get_attr( "xform.projection")
			torg = t.get_params('spider')
			print  ifil, iseg, mic[iseg], [torg['phi'], torg['theta'], torg['psi'], -torg['tx'], -torg['ty']], [iphi, THETA,PSI,SX, iy]
			#print  ifil, iseg, mic[iseg],  [iphi, THETA,PSI,SX, iy]

			if randfile != None:
				t = tmp.get_attr( "xform.projection")
				torg = t.get_params('spider')

				t = t*irota
				trot = t.get_params('spider')
				tmp.set_attr( "xform.projection" ,t)

				randphiy.append([torg['phi'], torg['ty'], trot['phi'], trot['ty'],ifil,(dA/pixel_size)/rise])
			tmp.set_attr( "ptcl_source_coord" ,coord)
			tmp.write_image(proj_stack, counter)
			counter += 1

	if randfile != None:
		write_text_row(randphiy, randfile)

def make_test_stack2(vol, stack, proj_stack, dphi, dp, pixel_size, fil_attr='filament', N = -1, randfile=None,phizero=True, yzero=False,xzero=True):

	'''
	INPUT
	
	vol: Name of volume from which projections are calculated
	stack: Name of data stack containing filament and ctf information
	proj_stack: Name of projection stack that will be generated from predicted parameters
	dphi: azimuthal rotation of helical symmetry
	dp: rise of helical symmetry
	pixel_size: pixel size
	fil_attr: attribute under which filament identifier is stored in header.
	randfile: name of text file to save the applied random rotations/shifts to
	
	OUTPUT
	
	Writes to disk a stack of projections named proj_stack calculated from vol using projection parameters predicted based on helical symmetry parameters.
	'''
	
	from random         import uniform, random
	from projection     import prep_vol, prgs
	from applications   import get_dist,ordersegments
	from filter			import filt_ctf
	from utilities		import get_im, write_text_row
	
	
	filaments = ordersegments(stack)
	ptclcoords = EMUtil.get_all_attributes(stack, 'ptcl_source_coord')

	if( N < 1 ): N = len(filaments)	
	volft, kbx, kby, kbz = prep_vol( get_im(vol) )

	THETA = 90.
	PSI   = 90.

	rise = dp/pixel_size  # in pixels
	permitrange = rise/2.0
	ima = EMData()
	counter = 0
	randphiy=[]
	for ifil in xrange(N):
		mic = filaments[ifil]
		nsegs = len(mic)
		
		phi0 = 0
		if not(phizero):
			phi0 = 360.0*random()
		y0 = 0
		if not(yzero):
			y0   = permitrange*random()
		s0 = 0
		if not(xzero):
			s0 = 10.0*random()
			if random > 0.5:
				s0 = s0 * -1.0
		refcoords = ptclcoords[mic[0]]

		for iseg in xrange(nsegs):
			dA = pixel_size * get_dist(ptclcoords[mic[iseg]], refcoords) # distance in Angstroms between segments mic[iseg] and mic[0]

			iphi = (phi0 + (dA/dp)*dphi)%360.0
			iy   = y0 + (dA%dp)/pixel_size # in pixels
			if( iy > permitrange ): iy -= rise

			ima.read_image(stack, mic[iseg], True)
			ct   = ima.get_attr('ctf')
			fifi = ima.get_attr(fil_attr)
			coord = ima.get_attr('ptcl_source_coord')
			tmp = filt_ctf(prgs(volft, kbz, [iphi, THETA,PSI, s0, iy], kbx, kby),ct)
			tmp.set_attr(fil_attr, fifi)

			tmp.set_attr( "ptcl_source_coord" ,coord)
			torg = (tmp.get_attr('xform.projection')).get_params('spider')
			
			randphiy.append([iseg, ifil, s0])
			tmp.write_image(proj_stack, counter)
			counter += 1
			
	if randfile != None:
		write_text_row(randphiy, randfile)

def rotline(alpha, v0):
	'''
	alpha: in degrees
	
	rotates vector with endpoints (0,0) and v0 by angle alpha around origin, and return the slope of the line. 
	'''
	from math import pi, cos, sin
	
	qv = pi/180.0
	
	# (0,0) and (x0,y0) are two points on the line l
	x0 = v0[0]*cos(qv*alpha) - v0[1]*sin(qv*alpha)
	y0 = v0[0]*sin(qv*alpha) + v0[1]*cos(qv*alpha)
	
	if x0 != 0:
		# slope of line l
		m = y0/x0
		return m
	else:
		return 0

def Gehelix_MPI(stack, ref_vol, outdir, delta, psi_max, search_rng, range, ywobble, pixel_size, dp, dphi, fract, rmax, rmin, maskfile = None, \
	    maxit = 1, CTF = False, snr = 1.0, sym = "c1",  user_func_name = "helical", npad = 2, debug = False):

	from alignment       import Numrinit, prepare_refrings, proj_ali_incore, proj_ali_incore_local, proj_ali_incore_local_psi
	from alignment       import ringwe, ang_n
	from utilities       import model_circle, get_image, drop_image, get_input_from_string, peak_search, model_cylinder, pad, model_blank
	from utilities       import bcast_list_to_all, bcast_number_to_all, reduce_EMData_to_root, bcast_EMData_to_all
	from utilities       import send_attr_dict, sym_vol, get_input_from_string
	from utilities       import get_params_proj, set_params_proj, file_type, compose_transform2
	from fundamentals    import rot_avg_image, ccf, fft, rot_shift2D
	import os
	import types
	from utilities       import print_begin_msg, print_end_msg, print_msg
	from mpi             import mpi_bcast, mpi_comm_size, mpi_comm_rank, MPI_FLOAT, MPI_COMM_WORLD, mpi_barrier
	from mpi             import mpi_reduce, MPI_INT, MPI_SUM
	from filter          import filt_ctf
	from projection      import prep_vol, prgs
	#from statistics      import hist_list, varf3d_MPI, fsc_mask
	from applications	 import MPI_start_end, ordersegments, header
	from time            import time	
	from copy 			 import copy
	from math 			 import sqrt

	nproc     = mpi_comm_size(MPI_COMM_WORLD)
	myid      = mpi_comm_rank(MPI_COMM_WORLD)
	main_node = 0

	search_rng   = int(search_rng)
	range        = int(range)
	if(pixel_size < 0.0 or dp < 0.0 ):  ERROR('Helical symmetry parameters have to be provided', "ehelix_MPI", 1, myid)

	if os.path.exists(outdir):  ERROR('Output directory %s  exists, please change the name and restart the program'%outdir, "ehelix_MPI", 1, myid)
	mpi_barrier(MPI_COMM_WORLD)

	if myid == main_node:
		print_begin_msg("ehelix_MPI")
		os.mkdir(outdir)
	mpi_barrier(MPI_COMM_WORLD)

	if debug:
		from time import sleep
		while not os.path.exists(outdir):
			print  "Node ",myid,"  waiting..."
			sleep(5)

		info_file = os.path.join(outdir, "progress%04d"%myid)
		finfo = open(info_file, 'w')
	else:
		finfo = None
	max_iter    = int(maxit)
	
	vol     = EMData()
	vol.read_image(ref_vol)
	nx      = vol.get_xsize()
	ny      = vol.get_ysize()
	nz      = vol.get_zsize()

	if(sym[0] == "d"  or sym[0] == "D"):  Dsym = True
	else:                                 Dsym = False

	xysize = nx
	zsize = -1
	#  For the time being only one delta!!!
	delta       = get_input_from_string(delta)[0]

	if myid == main_node:
		import user_functions
		user_func = user_functions.factory[user_func_name]

		print_msg("Input stack                               : %s\n"%(stack))
		print_msg("Reference volume                          : %s\n"%(ref_vol))	
		print_msg("Output directory                          : %s\n"%(outdir))
		print_msg("Maskfile                                  : %s\n"%(maskfile))
		print_msg("Angular step                              : %s\n"%(delta))
		print_msg("Maximum range for psi search              : %s\n"%(psi_max))
		print_msg("X-search range                            : %f\n"%(search_rng))
		print_msg("X-search wobble                           : %f\n"%(range))
		print_msg("Pixel size [A]                            : %f\n"%(pixel_size))
		print_msg("dp [A]                                    : %f\n"%(dp))
		print_msg("dphi                                      : %f\n"%(dphi))
		print_msg("Maximum iteration                         : %i\n"%(max_iter))
		print_msg("CTF correction                            : %s\n"%(CTF))
		print_msg("Signal-to-Noise Ratio                     : %f\n"%(snr))
		print_msg("Symmetry group                            : %s\n\n"%(sym))

	if maskfile:
		if type(maskfile) is types.StringType: mask3D = get_image(maskfile)
		else:                                  mask3D = maskfile
	else: mask3D = model_cylinder(rmax, nx, ny, nz)

	fscmask = mask3D  #model_circle(last_ring,nx,nx,nx)  For a fancy mask circle would work better  PAP 7/21/11
	if CTF:
		from reconstruction import recons3d_4nn_ctf_MPI
		from filter         import filt_ctf
	else:	 from reconstruction import recons3d_4nn_MPI

	filaments = ordersegments(stack)

	total_nfils = len(filaments)
	if myid == main_node:
		print "total number of filaments: ", total_nfils
	if total_nfils< nproc:
		ERROR('number of CPUs (%i) is larger than the number of filaments (%i), please reduce the number of CPUs used'%(nproc, total_nfils), "ehelix_MPI", 1,myid)

	fstart, fend = MPI_start_end(total_nfils, nproc, myid)
	filaments = filaments[fstart:fend]
	nfils     = len(filaments)

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

	data_nx = data[0].get_xsize()
	data_ny = data[0].get_xsize()
	mask2D  = pad(model_blank(2*int(rmax), data_ny, 1, 1.0), data_nx, data_ny, 1, 0.0)
	fdata = [None]*nima
	for im in xrange(nima):
		data[im].set_attr('ID', list_of_particles[im])
		if CTF:
			ctf_params = data[im].get_attr("ctf")
			st = Util.infomask(data[im], mask2D, False)
			data[im] -= st[0]
			data[im] = filt_ctf(data[im], ctf_params)
			data[im].set_attr('ctf_applied', 1)
		#  if FindPsi,  apply the angle to data[im], do fft and put in fdata[im]
		set_params_proj(data[im], [0.0, 90.0, 90.0, 0.0, 0.0])
		fdata[im] = fft( data[im] )
	del list_of_particles

	if debug:
		finfo.write( '%d loaded  \n' % nima )
		finfo.flush()
	if myid == main_node:
		# initialize data for the reference preparation function
		# No centering for helical reconstruction
		ref_data = [None, mask3D, None, None, None ]

	rise = int(dp/pixel_size)

	FindPsi = False  # !!!!!!!!

	phiwobble = int(float(ywobble)/rise*dphi/delta+0.5)
	print  " phi wobble ",phiwobble,phiwobble*delta


	nwx = 2*search_rng+3
	nwy = rise+2*ywobble+2
	nwxc = nwx//2
	nwyc = nwy//2
	nphi = int(360.0/delta + 0.5)
	#print  "  params  ",nwx,nwy,nwxc,nwyc,nphi
	if FindPsi:
		mode = "F"
		cnx = data_nx//2+1
		cny = cnx
		numr = Numrinit(1, data_nx//2-2, 1, mode)
		wr   = ringwe(numr, mode)
		maxrin = numr[len(numr)-1]
		crefim = [None]*nphi
	terminate = 0
	Iter = 0
 	while Iter < max_iter and terminate == 0:
		Iter += 1
		if myid == main_node:
			start_time = time()
			print_msg("\nITERATION #%3d\n"%(Iter))
		volft, kbx, kby, kbz = prep_vol( vol )

		#print  "ITERATION   ",Iter

		#del vol



		refproj = [None]*nphi
		if( not Dsym):  rotproj = [None]*nphi
		"""
		from fundamentals import rot_shift3D
		volft180x, kbx, kby, kbz = prep_vol( rot_shift3D(vol,-90,180.0,90) )
		volft180y, kbx, kby, kbz = prep_vol( rot_shift3D(vol,180.0) )
		"""
		refproj = EMData.read_images("refstraight.hdf")
		#Dsym = True  #  FOR DEBUGGING
		for iphi in xrange(nphi):
			#refproj[iphi] = Util.muln_img(prgs( volft, kbz, [delta*iphi, 90.0, 90.0, 0.0, 0.0], kbx, kby ), mask2D )
			#refproj[iphi].write_image("refstraight.hdf",iphi)
			"""
			fft(fft(refproj[iphi]).conjg()).write_image("refconj.hdf",iphi)
			rot_shift2D(refproj[iphi],180.).write_image("ref180.hdf",iphi)
			Util.muln_img(prgs( volft180y, kbz, [delta*iphi, 90.0, 90.0, 0.0, 0.0], kbx, kby ), mask2D ).write_image("ref180y.hdf",iphi)
			Util.muln_img(prgs( volft180x, kbz, [delta*iphi, 90.0, 90.0, 0.0, 0.0], kbx, kby ), mask2D ).write_image("ref180x.hdf",iphi)
			"""
			if FindPsi:
				temp = Util.Polar2Dm(refproj[iphi], cnx, cny, numr, mode)
				Util.Frngs(temp, numr)
				Util.Applyws(temp, numr, wr)
				crefim[iphi] = temp
			#  rotated in-plane by 180 are equivalent to rot_shift3D(vol,-90,180.0,90) with phi running as phi
			if(not Dsym):  rotproj[iphi] = fft( rot_shift2D(refproj[iphi],180.0) )
			refproj[iphi] = fft( refproj[iphi] )
		#exit()
		for ifil in xrange(nfils):
			if myid == main_node: print " process filament  ",ifil,time()-start_time; start_time = time()
			start = indcs[ifil][0]
			c0 = data[start].get_attr('ptcl_source_coord')
			ccfs = [[None for i in xrange(nphi)] for j in xrange(indcs[ifil][1] - start)]
			if(not Dsym):  ccfr = [[None for i in xrange(nphi)] for j in xrange(indcs[ifil][1] - start)]
			for im in xrange(start, indcs[ifil][1]):
				for iphi in xrange(nphi):
					ccfs[im-start][iphi] = Util.window(ccf( refproj[iphi], fdata[im] ), nwx, nwy)
					if(not Dsym):  ccfr[im-start][iphi] = Util.window(ccf( rotproj[iphi], fdata[im] ), nwx, nwy)
			if myid == main_node: print " ccfs computed     ",ifil,time()-start_time; start_time = time()
			xshiftlocal  = [0.0]*(indcs[ifil][1]-start)
			if(not Dsym): xrshiftlocal  = [0.0]*(indcs[ifil][1]-start)
			mxshiftlocal = [0.0]*(indcs[ifil][1]-start)
			yshiftlocal  = [0.0]*(indcs[ifil][1]-start)
			if(not Dsym): yrshiftlocal  = [0.0]*(indcs[ifil][1]-start)
			myshiftlocal = [0.0]*(indcs[ifil][1]-start)
			philocal  = [0.0]*(indcs[ifil][1]-start)
			if(not Dsym): phirlocal  = [0.0]*(indcs[ifil][1]-start)
			mphilocal = [0.0]*(indcs[ifil][1]-start)
			tmax = -1.0e23
			for ix in xrange(1,nwx-1):                                         #  X shift
				#print "im: ", len(ccfs), ix,time()-start_time
				six = ix - nwxc
				for iy in xrange(1+ywobble,nwy-ywobble-1):                     #  Y shift
					siy = iy - nwyc
					yshiftlocal[0] = float(iy-nwyc)
					for iphi in xrange(nphi):                                  #  phi search
						#qphi = iphi*delta
						philocal[0]  = iphi*delta
						phirlocal[0] = iphi*delta
						# we use the first segment as a reference, so there is no interpolation, just copy the correlation
						#  Select largest correlation within +/- range pixels from the location we explore
						mxm = -1.0e23
						for iux in xrange(max(1, ix - range), min(nwx - 1, ix+range+1), 1):    #  X wobble
							qcf = ccfs[0][iphi].get_value_at(iux,iy)
							if(qcf > mxm):
								mxm = qcf
								xshiftlocal[0] = float(iux-nwxc)
						#mxm = ccfs[0][iphi].get_value_at(ix,iy)
						if(not Dsym):
							mxr = -1.0e23
							for iux in xrange(max(1, ix - range), min(nwx - 1, ix+range+1), 1):     # Xr wobble
								qcf = ccfr[0][iphi].get_value_at(iux,iy)
								if(qcf > mxr):
									mxr = qcf
									xrshiftlocal[0] = float(iux-nwxc)
							#mxr = ccfr[0][iphi].get_value_at(ix,iy)
						#print  ix,six,iy,siy,iphi,mxm,mxr,xshiftlocal[0],xrshiftlocal[0]
						for im in xrange(1,indcs[ifil][1]-start):                                    #  predicted locations
							#print "im: ", len(ccfs), im,time()-start_time
							# dst is distance between segment 0 and current segment in pixels
							cim = data[indcs[ifil][0] + im].get_attr('ptcl_source_coord')
							dst = sqrt((c0[0] - cim[0])**2 + (c0[1] - cim[1])**2)
							#dst = 15.0
							#print im,dst,rise
							# predict for all remaining segments assuming number 0
							#  has parameters (qphi, six, siy)
							# Assume for now inter-segment distances are multiples of rise -- jia
							pphi = (philocal[0] + (dst/rise)*dphi)%360.0                          #  predicted phi with full angular accuracy, not an integer  #USED TO BE -
							pix = six # predicted x shift
							piy = siy #  predicted y shift
							xix = pix + nwxc
							yiy = piy + nwyc
							#  Local x search
							fix = int(xix)
							xdif = xix - fix
							xrem = 1.0 - xdif
							fiy = int(yiy)
							ydif = yiy - fiy
							yrem = 1.0 - ydif
							ciq = -1.0e23
							# interpolate correlation at pphi
							ttphi = int(pphi/delta + 0.5)%nphi
							for lphi in xrange(-phiwobble,phiwobble+1):                                         #  phi wobble
								tphi = (ttphi+lphi)%nphi
								for iux in xrange(max(1, fix - range), min(nwx - 1, fix+range+1), 1):           #  X wobble
									for iuy in xrange(max(1, fiy - ywobble), min(nwy - 1, fiy+ywobble+1), 1):   #  Y wobble
										qcf = xrem*yrem*ccfs[im][tphi].get_value_at(iux,iuy) + xdif*yrem*ccfs[im][tphi].get_value_at(iux+1,iuy) + xrem*ydif*ccfs[im][tphi].get_value_at(iux,iuy+1) + xdif*ydif*ccfs[im][tphi].get_value_at(iux+1,iuy+1)
										if(qcf > ciq):
											ciq = qcf
											xshiftlocal[im] = iux + xdif - nwxc
											yshiftlocal[im] = iuy + ydif - nwyc
											philocal[im] = tphi
								#print  "straight ",ix,iy,iphi, "    ", six,siy,qphi, "    ",pix,piy,xix,yiy,pphi,   "    ",fix,fiy,tphi
							#ciq = xrem*yrem*ccfs[im][tphi].get_value_at(fix,fiy) + xdif*yrem*ccfs[im][tphi].get_value_at(fix+1,fiy) + xrem*ydif*ccfs[im][tphi].get_value_at(fix,fiy+1) + xdif*ydif*ccfs[im][tphi].get_value_at(fix+1,fiy+1) 
							mxm += ciq
							# now for rotated
							if(not Dsym):
								# Assume for now inter-segment distances are multiples of rise -- jia
								pphi = (phirlocal[0] + (dst/rise)*dphi)%360.0 #  predicted phi for rotated 180 defs with full angular accuracy, not an integer
								pix = six # predicted x shift
								piy = siy #  predicted y shift
								xix = pix + nwxc
								yiy = piy + nwyc
								fix = int(xix)
								xdif = xix - fix
								xrem = 1.0 - xdif
								fiy = int(yiy)
								ydif = yiy - fiy
								yrem = 1.0 - ydif
								ciq = -1.0e23
								# interpolate correlation at pphi
								ttphi = int(pphi/delta + 0.5)%nphi
								for lphi in xrange(-phiwobble,phiwobble+1):                                           #  phi wobble
									tphi = (ttphi+lphi)%nphi
									for iux in xrange(max(1, fix - range), min(nwx - 1, fix+range+1), 1):             #  X wobble
										for iuy in xrange(max(1, fiy - ywobble), min(nwy - 1, fiy+ywobble+1), 1):     #  Y wobble
											qcf = xrem*yrem*ccfs[im][tphi].get_value_at(iux,iuy) + xdif*yrem*ccfs[im][tphi].get_value_at(iux+1,iuy) + xrem*ydif*ccfs[im][tphi].get_value_at(iux,iuy+1) + xdif*ydif*ccfs[im][tphi].get_value_at(iux+1,iuy+1)
											if(qcf > ciq):
												ciq = qcf
												xrshiftlocal[im] = iux + xdif - nwxc
												yrshiftlocal[im] = iuy + ydif - nwyc
												phirlocal[im] = tphi
								#print  "rotated ",ix,iy,iphi, "    ", six,siy,qphi, "    ",pix,piy,xix,yiy,pphi,   "    ",fix,fiy,tphi,iux + xdif - nwxc,ciq
								#ciq = xrem*yrem*ccfr[im][tphi].get_value_at(fix,fiy) + xdif*yrem*ccfr[im][tphi].get_value_at(fix+1,fiy) + xrem*ydif*ccfr[im][tphi].get_value_at(fix,fiy+1) + xdif*ydif*ccfr[im][tphi].get_value_at(fix+1,fiy+1) 
								mxr += ciq
							else:
								mxr = mxm-1.e5
						# The parameters are stored only for the first segment, the remaining one will have to be recomputed
						#if myid == main_node:  print  mxm,mxr
						#if Iter == 1:  mxm = 2*mxr
						if( mxr > mxm ):
							if(mxr > tmax):
								tmax = mxr
								mpsi = 270.0
								#mphi = iphi*delta
								for im in xrange(len(xshiftlocal)):  mxshiftlocal[im] = xrshiftlocal[im]
								for im in xrange(len(yshiftlocal)):  myshiftlocal[im] = yrshiftlocal[im]
								for im in xrange(len(phirlocal)):    mphilocal[im]    = phirlocal[im]*delta
								#msx = six
								#msy = siy
								#if myid == main_node:  print  ifil,ix, iy, iphi,pphi,tphi,mxshiftlocal
						else:
							if(mxm > tmax):
								tmax = mxm
								mpsi = 90.0
								#mphi = iphi*delta
								for im in xrange(len(xshiftlocal)):  mxshiftlocal[im] = xshiftlocal[im]
								for im in xrange(len(yshiftlocal)):  myshiftlocal[im] = yshiftlocal[im]
								for im in xrange(len(philocal)):     mphilocal[im]    = philocal[im]*delta
								#msx = six
								#msy = siy
								#if myid == main_node:  print  ifil,ix, iy, iphi,mxm,mxshiftlocal
			#print "  PARAMETERS FOR  0 ",mphi, 90.0, mpsi,mxm,mxr
			for im in xrange(start, indcs[ifil][1]):
				# Do the prediction using set (mphi, theta=90., mpsi, msx, msy)
				#cim = data[im].get_attr('ptcl_source_coord')
				#dst = sqrt((c0[0] - cim[0])**2 + (c0[1] - cim[1])**2)
				imstart = im - start
				psx  = mxshiftlocal[imstart]
				psy  = myshiftlocal[imstart]
				pphi = mphilocal[imstart]
				"""
				if mpsi == 90.0:
					psx  = mxshiftlocal[imstart] #msx
					psy  = myshiftlocal[imstart]
					pphi = mphilocal[imstart]# ( mphilocal[imstart] + (dst/rise)*dphi )%360.0
				else:
					psx  = mxshiftlocal[imstart] #msx
					psy  = myshiftlocal[imstart]
					pphi = ( mphilocal[imstart] + (dst/rise)*dphi )%360.0             #  USED TO BE - !!!!!!!!!!!!!!!!!!!!!!
				"""
				#print "  PARAMETERS FOR IM ",im,pphi, 90.0, mpsi, psx, psy
				if FindPsi:
					iphi = int(pphi/delta + 0.5)%nphi
					#print  " ref number and current parameters reduced to 2D  ",iphi,0.0, psx, psy				
					#  I should only care what the previous residual angle was
					ophi, otheta, opsi3, opx3, opy3 = get_params_proj(data[im])
					#print " old 3D params in data ",ophi, otheta, opsi3, opx3, opy3
					if(abs(opsi3 - 90.) < abs(opsi3 - 270.0)):  gamma =  90.
					else:                                       gamma = 270.0
					oalpha, osx, osy, junk = compose_transform2(0, opx3, opy3, 1.0, gamma-opsi3, 0, 0, 1.0) # reduce 3D to 2D
					#print " old 3D params, -> 2D ",oalpha, osx, osy
					# combine previous with the current in plane
					#print " current 2D combined with old 2D rotation",oalpha, csx, csy
					#  Find what the shift is without the angle
					junk, nnsx, nnsy, junk = compose_transform2(0.0, psx, psy, 1.0, -oalpha, 0., 0., 1.0)
					#print " 2D shift without angle ",nnsx, nnsy

					#rot_shift2D(data[im], 0.0, nnsx, nnsy).write_image("shifted.hdf",im)
					#fft(refproj[iphi]).write_image("reference.hdf",im)
					#from utilities import info

					cimage = Util.Polar2Dm(data[im], cnx+nnsx, cny+nnsy, numr, mode)
					Util.Frngs(cimage, numr)
					temp = Util.Crosrng_msg_s( cimage, crefim[iphi], numr)

					#from utilities import write_text_file
					#write_text_file([temp[qqq] for qqq in xrange(maxrin)],"ccf1d.txt")
					#exit()

					ipr = int(psi_max*maxrin/360. + 0.5)
					if(mpsi == 270.0):  incpsi = maxrin//2
					else:               incpsi = 0
					qn = -1.0e23
					for ips in xrange(-ipr,ipr+1,1):
						tot = (ips + incpsi + maxrin)%maxrin
						tval = temp.get_value_at(tot)
						#print  ips,incpsi,tot,tval
						if(tval > qn):
							qn = tval
							bestang = ang_n(tot+1.0, mode, maxrin)
					#print " best angle ",bestang
					bestang = (bestang - (mpsi-90.0))%360.0
					#print " angle applied ",bestang
					#rot_shift2D(data[im],-bestang).write_image("rotated.hdf",im)
					fdata[im] = fft( rot_shift2D(data[im], -bestang) )
					#print  " New composed 3D  ",mpsi,bestang, nnsx, nnsy

					epsi = (bestang+mpsi)%360.0
					psx = nnsx; psy = nnsy
					#print  " New composed 3D  ",pphi,90.0,epsi, psx, psy
					#exit()
				else:
					epsi = mpsi
				print  "  %7.2f  %7.2f  %7.2f  %7.2f  %7.2f"%(pphi, 90.0, epsi, psx, psy)
				set_params_proj(data[im], [pphi, 90.0, epsi, psx, psy])
				#print get_params_proj(data[im])
			#exit()

			if myid == main_node: print " params computed  ",ifil,time()-start_time; start_time = time()

		del ccfs
		if(not Dsym):  del ccfr

		del refproj, volft
		if(not Dsym):  del rotproj

		#print  " ======================================================================================================================="

		if myid == main_node:
			print_msg("Time of alignment = %d\n"%(time()-start_time))
			start_time = time()


		"""
		if CTF:  vol = recons3d_4nn_ctf_MPI(myid, data, symmetry=sym, snr = snr, npad = npad, xysize = xysize, zsize = zsize)
		else:    vol = recons3d_4nn_MPI(myid, data, symmetry=sym, npad = npad, xysize = xysize, zsize = zsize)

		if myid == main_node:
			print_msg("3D reconstruction time = %d\n"%(time()-start_time))
			start_time = time()

		if myid == main_node:
			vol = vol.helicise(pixel_size, dp, dphi, fract, rmax, rmin)
			vol = sym_vol(vol, symmetry=sym)
			ref_data[0] = vol
			vol = user_func(ref_data)
			vol = vol.helicise(pixel_size, dp, dphi, fract, rmax, rmin)
			vol = sym_vol(vol, symmetry=sym)
			vol.write_image(os.path.join(outdir, "volf%03d.hdf"%Iter))
			#if(Iter == max_iter-1):  drop_image(vol, os.path.join(outdir, "volfshift.hdf"))

		bcast_EMData_to_all(vol, myid, main_node)
		# write out headers, under MPI writing has to be done sequentially
		mpi_barrier(MPI_COMM_WORLD)
		par_str = ['xform.projection', 'ID']
		if myid == main_node:
			if(file_type(stack) == "bdb"):
				from utilities import recv_attr_dict_bdb
				recv_attr_dict_bdb(main_node, stack, data, par_str, 0, nima, nproc)
			else:
				from utilities import recv_attr_dict
				recv_attr_dict(main_node, stack, data, par_str, 0, nima, nproc)
			print_msg("Time to write header information= %d\n"%(time()-start_time))
			start_time = time()
		else:		send_attr_dict(main_node, data, par_str, 0, nima)
		if myid == main_node:
			# write params to text file
			header(stack, params='xform.projection', fexport=os.path.join(outdir, "parameters%04d.txt"%Iter))
		"""
	if myid == main_node: print_end_msg("ehelix_MPI")


def Iehelix_MPI(stack, ref_vol, outdir, delta, psi_max, search_rng, range, pixel_size, dp, dphi, fract, rmax, rmin, maskfile = None, \
	    maxit = 1, CTF = False, snr = 1.0, sym = "c1",  user_func_name = "helical", npad = 2, debug = False):

	from alignment       import Numrinit, prepare_refrings, proj_ali_incore, proj_ali_incore_local, proj_ali_incore_local_psi
	from alignment       import ringwe, ang_n
	from utilities       import model_circle, get_image, drop_image, get_input_from_string, peak_search, model_cylinder, pad, model_blank
	from utilities       import bcast_list_to_all, bcast_number_to_all, reduce_EMData_to_root, bcast_EMData_to_all
	from utilities       import send_attr_dict, sym_vol, get_input_from_string
	from utilities       import get_params_proj, set_params_proj, file_type, compose_transform2
	from fundamentals    import rot_avg_image, ccf, fft, rot_shift2D
	import os
	import types
	from utilities       import print_begin_msg, print_end_msg, print_msg
	from mpi             import mpi_bcast, mpi_comm_size, mpi_comm_rank, MPI_FLOAT, MPI_COMM_WORLD, mpi_barrier
	from mpi             import mpi_reduce, MPI_INT, MPI_SUM
	from filter          import filt_ctf
	from projection      import prep_vol, prgs
	#from statistics      import hist_list, varf3d_MPI, fsc_mask
	from applications	 import MPI_start_end, ordersegments, header
	from time            import time	
	from copy 			 import copy
	from math 			 import sqrt
	
	from utilities import info
	
	nproc     = mpi_comm_size(MPI_COMM_WORLD)
	myid      = mpi_comm_rank(MPI_COMM_WORLD)
	main_node = 0

	search_rng   = int(search_rng)
	range        = int(range)
	if(pixel_size < 0.0 or dp < 0.0 ):  ERROR('Helical symmetry parameters have to be provided', "ehelix_MPI", 1, myid)

	if os.path.exists(outdir):  ERROR('Output directory %s  exists, please change the name and restart the program'%outdir, "ehelix_MPI", 1, myid)
	mpi_barrier(MPI_COMM_WORLD)

	if myid == main_node:
		print_begin_msg("ehelix_MPI")
		os.mkdir(outdir)
	mpi_barrier(MPI_COMM_WORLD)

	if debug:
		from time import sleep
		while not os.path.exists(outdir):
			print  "Node ",myid,"  waiting..."
			sleep(5)

		info_file = os.path.join(outdir, "progress%04d"%myid)
		finfo = open(info_file, 'w')
	else:
		finfo = None
	max_iter    = int(maxit)
	
	vol     = EMData()
	vol.read_image(ref_vol)
	nx      = vol.get_xsize()
	ny      = vol.get_ysize()
	nz      = vol.get_zsize()

	if(sym[0] == "d"  or sym[0] == "D"):  Dsym = True
	else:                                 Dsym = False

	xysize = nx
	zsize = -1
	#  For the time being only one delta!!!
	delta       = get_input_from_string(delta)[0]

	if myid == main_node:
		import user_functions
		user_func = user_functions.factory[user_func_name]

		print_msg("Input stack                               : %s\n"%(stack))
		print_msg("Reference volume                          : %s\n"%(ref_vol))	
		print_msg("Output directory                          : %s\n"%(outdir))
		print_msg("Maskfile                                  : %s\n"%(maskfile))
		print_msg("Angular step                              : %s\n"%(delta))
		print_msg("Maximum range for psi search              : %s\n"%(psi_max))
		print_msg("X-search range                            : %f\n"%(search_rng))
		print_msg("X-search wobble                           : %f\n"%(range))
		print_msg("Pixel size [A]                            : %f\n"%(pixel_size))
		print_msg("dp [A]                                    : %f\n"%(dp))
		print_msg("dphi                                      : %f\n"%(dphi))
		print_msg("Maximum iteration                         : %i\n"%(max_iter))
		print_msg("CTF correction                            : %s\n"%(CTF))
		print_msg("Signal-to-Noise Ratio                     : %f\n"%(snr))
		print_msg("Symmetry group                            : %s\n\n"%(sym))

	if maskfile:
		if type(maskfile) is types.StringType: mask3D = get_image(maskfile)
		else:                                  mask3D = maskfile
	else: mask3D = model_cylinder(rmax, nx, ny, nz)

	fscmask = mask3D  #model_circle(last_ring,nx,nx,nx)  For a fancy mask circle would work better  PAP 7/21/11
	if CTF:
		from reconstruction import recons3d_4nn_ctf_MPI
		from filter         import filt_ctf
	else:	 from reconstruction import recons3d_4nn_MPI

	filaments = ordersegments(stack)
	print filaments
	total_nfils = len(filaments)
	if myid == main_node:
		print "total number of filaments: ", total_nfils
	if total_nfils< nproc:
		ERROR('number of CPUs (%i) is larger than the number of filaments (%i), please reduce the number of CPUs used'%(nproc, total_nfils), "ehelix_MPI", 1,myid)

	fstart, fend = MPI_start_end(total_nfils, nproc, myid)
	filaments = filaments[fstart:fend]
	nfils     = len(filaments)

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

	data_nx = data[0].get_xsize()
	data_ny = data[0].get_xsize()
	mask2D  = pad(model_blank(2*int(rmax), data_ny, 1, 1.0), data_nx, data_ny, 1, 0.0)
	fdata = [None]*nima
	for im in xrange(nima):
		data[im].set_attr('ID', list_of_particles[im])
		if CTF:
			ctf_params = data[im].get_attr("ctf")
			st = Util.infomask(data[im], mask2D, False)
			data[im] -= st[0]
			data[im] = filt_ctf(data[im], ctf_params)
			data[im].set_attr('ctf_applied', 1)
		#  if FindPsi,  apply the angle to data[im], do fft and put in fdata[im]
		set_params_proj(data[im], [0.0, 90.0, 90.0, 0.0,0.0])
		#fdata[im] = fft( data[im] )
	del list_of_particles

	if debug:
		finfo.write( '%d loaded  \n' % nima )
		finfo.flush()
	if myid == main_node:
		# initialize data for the reference preparation function
		# No centering for helical reconstruction
		ref_data = [None, mask3D, None, None, None ]

	rise = int(dp/pixel_size)

	FindPsi = True  # !!!!!!!!
	ywobble = 0
	phiwobble = int(float(ywobble)/rise*dphi/delta+0.5)
	print  " phi wobble ",phiwobble,phiwobble*delta

	mask2D  = pad(model_blank(2*int(rmax), 3*rise, 1, 1.0), data_nx, data_ny, 1, 0.0)
	for im in xrange(nima):
		data[im]*=mask2D
		fdata[im] = fft( data[im] )



	nwx = 2*search_rng+3
	nwy = rise+2*ywobble+2
	nwxc = nwx//2
	nwyc = nwy//2
	nphi = int(360.0/delta + 0.5)
	#print  "  params  ",nwx,nwy,nwxc,nwyc,nphi
	if FindPsi:
		mode = "F"
		cnx = data_nx//2+1
		cny = cnx
		numr = Numrinit(1, data_nx//2-2, 1, mode)
		wr   = ringwe(numr, mode)
		maxrin = numr[len(numr)-1]
		crefim = [None]*nphi
	terminate = 0
	Iter = 0
 	while Iter < max_iter and terminate == 0:
		Iter += 1
		if myid == main_node:
			start_time = time()
			print_msg("\nITERATION #%3d\n"%(Iter))
		volft, kbx, kby, kbz = prep_vol( vol )

		print  "ITERATION   ",Iter

		#del vol



		refproj = [None]*nphi
		if( not Dsym):  rotproj = [None]*nphi
		"""
		from fundamentals import rot_shift3D
		volft180x, kbx, kby, kbz = prep_vol( rot_shift3D(vol,-90,180.0,90) )
		volft180y, kbx, kby, kbz = prep_vol( rot_shift3D(vol,180.0) )
		"""
		refproj = EMData.read_images("refstraight.hdf")
		#Dsym = True  #  FOR DEBUGGING
		for iphi in xrange(nphi):
			#refproj[iphi] = Util.muln_img(prgs( volft, kbz, [delta*iphi, 90.0, 90.0, 0.0, 0.0], kbx, kby ), mask2D )
			#refproj[iphi].write_image("refstraight.hdf",iphi)
			"""
			fft(fft(refproj[iphi]).conjg()).write_image("refconj.hdf",iphi)
			rot_shift2D(refproj[iphi],180.).write_image("ref180.hdf",iphi)
			Util.muln_img(prgs( volft180y, kbz, [delta*iphi, 90.0, 90.0, 0.0, 0.0], kbx, kby ), mask2D ).write_image("ref180y.hdf",iphi)
			Util.muln_img(prgs( volft180x, kbz, [delta*iphi, 90.0, 90.0, 0.0, 0.0], kbx, kby ), mask2D ).write_image("ref180x.hdf",iphi)
			"""
			if FindPsi:
				temp = Util.Polar2Dm(refproj[iphi], cnx, cny, numr, mode)
				Util.Frngs(temp, numr)
				Util.Applyws(temp, numr, wr)
				crefim[iphi] = temp
			#  rotated in-plane by 180 are equivalent to rot_shift3D(vol,-90,180.0,90) with phi running as phi
			if(not Dsym):  rotproj[iphi] = fft( rot_shift2D(refproj[iphi],180.0) )
			refproj[iphi] = fft( refproj[iphi] )
		#exit()
		for ifil in xrange(nfils):
			if myid == main_node: print " process filament  ",ifil,time()-start_time; start_time = time()
			start = indcs[ifil][0]
			c0 = data[start].get_attr('ptcl_source_coord')
			ccfs = [[None for i in xrange(nphi)] for j in xrange(indcs[ifil][1] - start)]
			if(not Dsym):  ccfr = [[None for i in xrange(nphi)] for j in xrange(indcs[ifil][1] - start)]
			for im in xrange(start, indcs[ifil][1]):
				for iphi in xrange(nphi):
					ccfs[im-start][iphi] = Util.window(ccf( refproj[iphi], fdata[im] ), nwx, nwy)
					if(not Dsym):  ccfr[im-start][iphi] = Util.window(ccf( rotproj[iphi], fdata[im] ), nwx, nwy)
			if myid == main_node: print " ccfs computed     ",ifil,time()-start_time; start_time = time()
			xshiftlocal  = [0.0]*(indcs[ifil][1]-start)
			if(not Dsym): xrshiftlocal  = [0.0]*(indcs[ifil][1]-start)
			mxshiftlocal = [0.0]*(indcs[ifil][1]-start)
			yshiftlocal  = [0.0]*(indcs[ifil][1]-start)
			if(not Dsym): yrshiftlocal  = [0.0]*(indcs[ifil][1]-start)
			myshiftlocal = [0.0]*(indcs[ifil][1]-start)
			philocal  = [0.0]*(indcs[ifil][1]-start)
			if(not Dsym): phirlocal  = [0.0]*(indcs[ifil][1]-start)
			mphilocal = [0.0]*(indcs[ifil][1]-start)
			for im in xrange(indcs[ifil][1]-start):                                 #  Loop over images
				tmax = -1.0e23
				mxm = -1.0e23
				mxr = -1.0e23
				for ix in xrange(1,nwx-1):                                            #   X search
					#print "im: ", len(ccfs), ix,time()-start_time
					six = ix - nwxc
					for iy in xrange(1+ywobble,nwy-ywobble-1):                        #   Y search
						siy = iy - nwyc
						for iphi in xrange(nphi):                                      #  phi search
							qphi = iphi*delta
							# we use the first segment as a reference, so there is no interpolation, just copy the correlation
							#  Select largest correlation within +/- range pixels from the location we explore
							for iux in xrange(max(1, ix - range), min(nwx - 1, ix+range+1), 1):
								qcf = ccfs[im][iphi].get_value_at(iux,iy)
								if(qcf > mxm):
									mxm = qcf
									xshiftlocal[im] = float(iux-nwxc)
									yshiftlocal[im] = float(iy-nwyc)
									philocal[im]    = qphi
							#mxm = ccfs[0][iphi].get_value_at(ix,iy)
							if(not Dsym):
								for iux in xrange(max(1, ix - range), min(nwx - 1, ix+range+1), 1):
									qcf = ccfr[im][iphi].get_value_at(iux,iy)
									if(qcf > mxr):
										mxr = qcf
										xrshiftlocal[im] = float(iux-nwxc)
										yrshiftlocal[im] = float(iy-nwyc)
										phirlocal[im]    = qphi
								#mxr = ccfr[0][iphi].get_value_at(ix,iy)
							#print  ix,six,iy,siy,iphi,mxm,mxr,xshiftlocal[0],xrshiftlocal[0]

				#if myid == main_node:  print  mxm,mxr
				#if Iter == 1:  mxm = 2*mxr
				if( mxr > mxm ):
					tmax = mxr
					mpsi = 270.0
					#mphi = iphi*delta
					#for im in xrange(len(xshiftlocal)):  mxshiftlocal[im] = xrshiftlocal[im]
					#for im in xrange(len(yshiftlocal)):  myshiftlocal[im] = yrshiftlocal[im]
					#for im in xrange(len(phirlocal)):    mphilocal[im]    = phirlocal[im]*delta
					#msx = six
					#msy = siy
					#if myid == main_node:  print  ifil,ix, iy, iphi,pphi,tphi,mxshiftlocal
					print  " %7.2f  %7.2f  %7.2f  %7.2f  %7.2f"%(phirlocal[im]*delta,90.0,mpsi,xrshiftlocal[im],yrshiftlocal[im]),im
				else:
					tmax = mxm
					mpsi = 90.0
					#mphi = iphi*delta
					#for im in xrange(len(xshiftlocal)):  mxshiftlocal[im] = xshiftlocal[im]
					#for im in xrange(len(yshiftlocal)):  myshiftlocal[im] = yshiftlocal[im]
					#for im in xrange(len(philocal)):     mphilocal[im]    = philocal[im]*delta
					#msx = six
					#msy = siy
					#if myid == main_node:  print  ifil,ix, iy, iphi,mxm,mxshiftlocal
					print  " %7.2f  %7.2f  %7.2f  %7.2f  %7.2f"%(philocal[im]*delta,90.0,mpsi,xshiftlocal[im],yshiftlocal[im]),im
		exit()

"""
			#print "  PARAMETERS FOR  0 ",mphi, 90.0, mpsi,mxm,mxr
			for im in xrange(start, indcs[ifil][1]):
				# Do the prediction using set (mphi, theta=90., mpsi, msx, msy)
				cim = data[im].get_attr('ptcl_source_coord')
				dst = sqrt((c0[0] - cim[0])**2 + (c0[1] - cim[1])**2)
				imstart = im - start
				if mpsi == 90.0:
					psx  = mxshiftlocal[imstart] #msx
					psy  = myshiftlocal[imstart]
					pphi = ( mphilocal[imstart] + (dst/rise)*dphi )%360.0
				else:
					psx  = mxshiftlocal[imstart] #msx
					psy  = myshiftlocal[imstart]
					pphi = ( mphilocal[imstart] + (dst/rise)*dphi )%360.0             #  USED TO BE - !!!!!!!!!!!!!!!!!!!!!!
				#print "  PARAMETERS FOR IM ",im,pphi, 90.0, mpsi, psx, psy
				if FindPsi:
					iphi = int(pphi/delta + 0.5)%nphi
					#print  " ref number and current parameters reduced to 2D  ",iphi,0.0, psx, psy				
					#  I should only care what the previous residual angle was
					ophi, otheta, opsi3, opx3, opy3 = get_params_proj(data[im])
					#print " old 3D params in data ",ophi, otheta, opsi3, opx3, opy3
					if(abs(opsi3 - 90.) < abs(opsi3 - 270.0)):  gamma =  90.
					else:                                       gamma = 270.0
					oalpha, osx, osy, junk = compose_transform2(0, opx3, opy3, 1.0, gamma-opsi3, 0, 0, 1.0) # reduce 3D to 2D
					#print " old 3D params, -> 2D ",oalpha, osx, osy
					# combine previous with the current in plane
					#print " current 2D combined with old 2D rotation",oalpha, csx, csy
					#  Find what the shift is without the angle
					junk, nnsx, nnsy, junk = compose_transform2(0.0, psx, psy, 1.0, -oalpha, 0., 0., 1.0)
					#print " 2D shift without angle ",nnsx, nnsy

					#rot_shift2D(data[im], 0.0, nnsx, nnsy).write_image("shifted.hdf",im)
					#fft(refproj[iphi]).write_image("reference.hdf",im)
					#from utilities import info

					cimage = Util.Polar2Dm(data[im], cnx+nnsx, cny+nnsy, numr, mode)
					Util.Frngs(cimage, numr)
					temp = Util.Crosrng_msg_s( cimage, crefim[iphi], numr)

					#from utilities import write_text_file
					#write_text_file([temp[qqq] for qqq in xrange(maxrin)],"ccf1d.txt")
					#exit()

					ipr = int(psi_max*maxrin/360. + 0.5)
					if(mpsi == 270.0):  incpsi = maxrin//2
					else:               incpsi = 0
					qn = -1.0e20
					for ips in xrange(-ipr,ipr+1,1):
						tot = (ips + incpsi + maxrin)%maxrin
						tval = temp.get_value_at(tot)
						#print  ips,incpsi,tot,tval
						if(tval > qn):
							qn = tval
							bestang = ang_n(tot+1.0, mode, maxrin)
					#print " best angle ",bestang
					bestang = (bestang - (mpsi-90.0))%360.0
					#print " angle applied ",bestang
					#rot_shift2D(data[im],-bestang).write_image("rotated.hdf",im)
					fdata[im] = fft( rot_shift2D(data[im], -bestang) )
					#print  " New composed 3D  ",mpsi,bestang, nnsx, nnsy

					epsi = (bestang+mpsi)%360.0
					psx = nnsx; psy = nnsy
					#print  " New composed 3D  ",pphi,90.0,epsi, psx, psy
					#exit()
				else:
					epsi = mpsi
				#print  pphi, 90.0, epsi, psx, psy
				set_params_proj(data[im], [pphi, 90.0, epsi, psx, psy])
				print get_params_proj(data[im])
			#exit()

			if myid == main_node: print " params computed  ",ifil,time()-start_time; start_time = time()
		del ccfs
		if(not Dsym):  del ccfr

		del refproj, volft
		if(not Dsym):  del rotproj

		#print  " ======================================================================================================================="

		if myid == main_node:
			print_msg("Time of alignment = %d\n"%(time()-start_time))
			start_time = time()

"""

"""
		if CTF:  vol = recons3d_4nn_ctf_MPI(myid, data, symmetry=sym, snr = snr, npad = npad, xysize = xysize, zsize = zsize)
		else:    vol = recons3d_4nn_MPI(myid, data, symmetry=sym, npad = npad, xysize = xysize, zsize = zsize)

		if myid == main_node:
			print_msg("3D reconstruction time = %d\n"%(time()-start_time))
			start_time = time()

		if myid == main_node:
			vol = vol.helicise(pixel_size, dp, dphi, fract, rmax, rmin)
			vol = sym_vol(vol, symmetry=sym)
			ref_data[0] = vol
			vol = user_func(ref_data)
			vol = vol.helicise(pixel_size, dp, dphi, fract, rmax, rmin)
			vol = sym_vol(vol, symmetry=sym)
			vol.write_image(os.path.join(outdir, "volf%03d.hdf"%Iter))
			#if(Iter == max_iter-1):  drop_image(vol, os.path.join(outdir, "volfshift.hdf"))

		bcast_EMData_to_all(vol, myid, main_node)
		# write out headers, under MPI writing has to be done sequentially
		mpi_barrier(MPI_COMM_WORLD)
		par_str = ['xform.projection', 'ID']
		if myid == main_node:
			if(file_type(stack) == "bdb"):
				from utilities import recv_attr_dict_bdb
				recv_attr_dict_bdb(main_node, stack, data, par_str, 0, nima, nproc)
			else:
				from utilities import recv_attr_dict
				recv_attr_dict(main_node, stack, data, par_str, 0, nima, nproc)
			print_msg("Time to write header information= %d\n"%(time()-start_time))
			start_time = time()
		else:		send_attr_dict(main_node, data, par_str, 0, nima)
		if myid == main_node:
			# write params to text file
			header(stack, params='xform.projection', fexport=os.path.join(outdir, "parameters%04d.txt"%Iter))

#	if myid == main_node: print_end_msg("ehelix_MPI")
"""

'''
defXehelix_MPI(stack, ref_vol, outdir, delta, psi_max, search_rng, range, ywobble, pixel_size, dp, dphi, fract, rmax, rmin, maskfile = None, \
	    maxit = 1, CTF = False, snr = 1.0, sym = "c1",  user_func_name = "helical", npad = 2, debug = False):

	from alignment       import Numrinit, prepare_refrings, proj_ali_incore, proj_ali_incore_local, proj_ali_incore_local_psi
	from alignment       import ringwe, ang_n
	from utilities       import model_circle, get_image, drop_image, get_input_from_string, peak_search, model_cylinder, pad, model_blank
	from utilities       import bcast_list_to_all, bcast_number_to_all, reduce_EMData_to_root, bcast_EMData_to_all
	from utilities       import send_attr_dict, sym_vol, get_input_from_string
	from utilities       import get_params_proj, set_params_proj, file_type, compose_transform2
	from fundamentals    import rot_avg_image, ccf, fft, rot_shift2D
	import os
	import types
	from utilities       import print_begin_msg, print_end_msg, print_msg
	from mpi             import mpi_bcast, mpi_comm_size, mpi_comm_rank, MPI_FLOAT, MPI_COMM_WORLD, mpi_barrier
	from mpi             import mpi_reduce, MPI_INT, MPI_SUM
	from filter          import filt_ctf
	from projection      import prep_vol, prgs
	#from statistics      import hist_list, varf3d_MPI, fsc_mask
	from applications	 import MPI_start_end, header
	from pixel_error     import ordersegments
	from time            import time
	from copy 			 import copy
	from math 			 import sqrt

	def rot2pad(imi, alpha=0.0, sx=0.0, sy=0.0):
		from utilities    import pad
		from fundamentals import rot_shift2D
		lnx = imi.get_xsize()
		lny = imi.get_ysize()
		ln = max(lnx,lny)
		if lnx == lny: return rot_shift2D(imi,alpha,sx,sy)
		else:          return Util.window(rot_shift2D(pad(imi,ln,ln,1,"circumference"), alpha,sx,sy), lnx, lny,1, 0,0,0)


	nproc     = mpi_comm_size(MPI_COMM_WORLD)
	myid      = mpi_comm_rank(MPI_COMM_WORLD)
	main_node = 0

	search_rng   = int(search_rng)
	range        = int(range)
	if(pixel_size < 0.0 or dp < 0.0 ):  ERROR('Helical symmetry parameters have to be provided', "ehelix_MPI", 1, myid)

	if os.path.exists(outdir):  ERROR('Output directory %s  exists, please change the name and restart the program'%outdir, "ehelix_MPI", 1, myid)
	mpi_barrier(MPI_COMM_WORLD)

	if myid == main_node:
		os.mkdir(outdir)
		import global_def
		global_def.LOGFILE =  os.path.join(outdir, global_def.LOGFILE)
		print_begin_msg("ehelix_MPI")
	mpi_barrier(MPI_COMM_WORLD)

	if debug:
		from time import sleep
		while not os.path.exists(outdir):
			print  "Node ",myid,"  waiting..."
			sleep(5)

		info_file = os.path.join(outdir, "progress%04d"%myid)
		finfo = open(info_file, 'w')
	else:
		finfo = None
	max_iter = int(maxit)
	
	vol     = EMData()
	vol.read_image(ref_vol)
	nx      = vol.get_xsize()
	ny      = vol.get_ysize()
	nz      = vol.get_zsize()

	if(sym[0] == "d"  or sym[0] == "D"):  Dsym = True
	else:                                 Dsym = False

	#  For the time being only one delta!!!
	delta       = get_input_from_string(delta)[0]

	if myid == main_node:
		import user_functions
		user_func = user_functions.factory[user_func_name]

		print_msg("Input stack                               : %s\n"%(stack))
		print_msg("Reference volume                          : %s\n"%(ref_vol))	
		print_msg("Output directory                          : %s\n"%(outdir))
		print_msg("Maskfile                                  : %s\n"%(maskfile))
		print_msg("Angular step                              : %s\n"%(delta))
		print_msg("Maximum range for psi search              : %s\n"%(psi_max))
		print_msg("X-search range                            : %f\n"%(search_rng))
		print_msg("X-search wobble                           : %f\n"%(range))
		print_msg("Pixel size [A]                            : %f\n"%(pixel_size))
		print_msg("dp [A]                                    : %f\n"%(dp))
		print_msg("dphi                                      : %f\n"%(dphi))
		print_msg("Maximum iteration                         : %i\n"%(max_iter))
		print_msg("CTF correction                            : %s\n"%(CTF))
		print_msg("Signal-to-Noise Ratio                     : %f\n"%(snr))
		print_msg("Symmetry group                            : %s\n\n"%(sym))

	if maskfile:
		if type(maskfile) is types.StringType: mask3D = get_image(maskfile)
		else:                                  mask3D = maskfile
	else: mask3D = model_cylinder(rmax, nx, ny, nz)

	fscmask = mask3D
	if CTF:
		from reconstruction import recons3d_4nn_ctf_MPI
		from filter         import filt_ctf
	else:	 from reconstruction import recons3d_4nn_MPI
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

	fstart, fend = MPI_start_end(total_nfils, nproc, myid)
	filaments = filaments[fstart:fend]
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
	print  " READ IMAGES ", myid
	FindPsi = True  # !!!!!!!!

	rise = int(dp/pixel_size)

	Iter = 1
	data_nx = data[0].get_xsize()
	data_ny = data[0].get_ysize()
	data_nn = max(data_nx, data_ny)
	mask2D  = pad(model_blank(2*int(rmax), data_ny, 1, 1.0), data_nx, data_ny, 1, 0.0)
	fdata = [None]*nima
	for im in xrange(nima):
		data[im].set_attr('ID', list_of_particles[im])
		# write out headers, under MPI writing has to be done sequentially
	mpi_barrier(MPI_COMM_WORLD)
	start_time = time()
	par_str = ['xform.projection', 'ID']
	if myid == main_node:
		if(file_type(stack) == "bdb"):
			sts = time()
			print  "  writing !!! ",myid
			from utilities import recv_attr_dict_bdb
			recv_attr_dict_bdb(main_node, stack, data, par_str, 0, nima, nproc)
			print  " main wrote stuff",myid,time()-sts
		else:
			from utilities import recv_attr_dict
			recv_attr_dict(main_node, stack, data, par_str, 0, nima, nproc)
		print_msg("Time to write header information= %d\n"%(time()-start_time))
		start_time = time()
	else:
		print  "  sending header info  ",myid, time()
		send_attr_dict(main_node, data, par_str, 0, nima)
	if myid == main_node:
		# write params to text file
		header(stack, params='xform.projection', fexport=os.path.join(outdir, "parameters%04d.txt"%Iter))


	if myid == main_node: print_end_msg("ehelix_MPI")
'''

def filament_balanced_load(stack, nproc):
	from pixel_error     import ordersegments
	from utilities 		 import chunks_distribution
	infils = EMUtil.get_all_attributes(stack, "filament")
	ptlcoords = EMUtil.get_all_attributes(stack, 'ptcl_source_coord')
	filaments = ordersegments(infils, ptlcoords)
	del ptlcoords
	total_nfils = len(filaments)

	print "total number of filaments: ", total_nfils
	if total_nfils< nproc:
		ERROR('number of CPUs (%i) is larger than the number of filaments (%i), please reduce the number of CPUs used'%(nproc, total_nfils), "filament_balanced_load", 1,0)
	temp = chunks_distribution([[len(filaments[i]), i] for i in xrange(len(filaments))], nproc)
	del filaments
	qmin = 1000000
	qmax = -1
	for myid in xrange(nproc):
		q = 0
		for i in xrange(len(temp[myid])): q+= temp[myid][i][0]
		qmin = min(qmin,q)
		qmax = max(qmax,q)
		print  myid, len(temp[myid]), q
	print qmin,qmax 
'''
defXehelix_MPI(stack, ref_vol, outdir, delta, psi_max, search_rng, range, ywobble, pixel_size, dp, dphi, fract, rmax, rmin, FindPsi = True, maskfile = None, \
	    maxit = 1, CTF = False, snr = 1.0, sym = "c1",  user_func_name = "helical", npad = 2, debug = False):

	from alignment       import Numrinit, prepare_refrings, proj_ali_incore, proj_ali_incore_local, proj_ali_incore_local_psi
	from alignment       import ringwe, ang_n
	from utilities       import model_circle, get_image, drop_image, get_input_from_string, peak_search, model_cylinder, pad, model_blank
	from utilities       import bcast_list_to_all, bcast_number_to_all, reduce_EMData_to_root, bcast_EMData_to_all
	from utilities       import send_attr_dict, sym_vol, get_input_from_string
	from utilities       import get_params_proj, set_params_proj, file_type, compose_transform2
	from fundamentals    import rot_avg_image, ccf, fft, rot_shift2D
	import os
	import types
	from utilities       import print_begin_msg, print_end_msg, print_msg, chunks_distribution
	from mpi             import mpi_bcast, mpi_comm_size, mpi_comm_rank, MPI_FLOAT, MPI_COMM_WORLD, mpi_barrier
	from mpi             import mpi_reduce, MPI_INT, MPI_SUM
	from filter          import filt_ctf
	from projection      import prep_vol, prgs
	#from statistics      import hist_list, varf3d_MPI, fsc_mask
	from applications	 import MPI_start_end, header
	from pixel_error     import ordersegments
	from time            import time
	from copy 			 import copy
	from math 			 import sqrt

	def rot2pad(imi, alpha=0.0, sx=0.0, sy=0.0):
		from utilities    import pad
		from fundamentals import rot_shift2D
		lnx = imi.get_xsize()
		lny = imi.get_ysize()
		ln = max(lnx,lny)
		if lnx == lny: return rot_shift2D(imi,alpha,sx,sy)
		else:          return Util.window(rot_shift2D(pad(imi,ln,ln,1,"circumference"), alpha,sx,sy), lnx, lny,1, 0,0,0)


	nproc     = mpi_comm_size(MPI_COMM_WORLD)
	myid      = mpi_comm_rank(MPI_COMM_WORLD)
	main_node = 0

	search_rng   = int(search_rng)
	range        = int(range)
	if(pixel_size < 0.0 or dp < 0.0 ):  ERROR('Helical symmetry parameters have to be provided', "ehelix_MPI", 1, myid)

	if os.path.exists(outdir):  ERROR('Output directory %s  exists, please change the name and restart the program'%outdir, "ehelix_MPI", 1, myid)
	mpi_barrier(MPI_COMM_WORLD)

	if myid == main_node:
		os.mkdir(outdir)
		import global_def
		global_def.LOGFILE =  os.path.join(outdir, global_def.LOGFILE)
		print_begin_msg("ehelix_MPI")
	mpi_barrier(MPI_COMM_WORLD)

	if debug:
		from time import sleep
		while not os.path.exists(outdir):
			print  "Node ",myid,"  waiting..."
			sleep(5)

		info_file = os.path.join(outdir, "progress%04d"%myid)
		finfo = open(info_file, 'w')
	else:
		finfo = None
	max_iter = int(maxit)
	
	vol     = EMData()
	vol.read_image(ref_vol)
	nx      = vol.get_xsize()
	ny      = vol.get_ysize()
	nz      = vol.get_zsize()

	if(sym[0] == "d"  or sym[0] == "D"):  Dsym = True
	else:                                 Dsym = False

	#  For the time being only one delta!!!
	delta       = get_input_from_string(delta)[0]

	if myid == main_node:
		import user_functions
		user_func = user_functions.factory[user_func_name]

		print_msg("Input stack                               : %s\n"%(stack))
		print_msg("Reference volume                          : %s\n"%(ref_vol))	
		print_msg("Output directory                          : %s\n"%(outdir))
		print_msg("Maskfile                                  : %s\n"%(maskfile))
		print_msg("Angular step                              : %s\n"%(delta))
		print_msg("Search for psi                            : %s\n"%(FindPsi))
		if FindPsi:  print_msg("Maximum range for psi search              : %s\n"%(psi_max))
		print_msg("X-search range                            : %f\n"%(search_rng))
		print_msg("X-search wobble                           : %f\n"%(range))
		print_msg("Pixel size [A]                            : %f\n"%(pixel_size))
		print_msg("dp [A]                                    : %f\n"%(dp))
		print_msg("dphi                                      : %f\n"%(dphi))
		print_msg("Maximum iteration                         : %i\n"%(max_iter))
		print_msg("CTF correction                            : %s\n"%(CTF))
		print_msg("Signal-to-Noise Ratio                     : %f\n"%(snr))
		print_msg("Symmetry group                            : %s\n\n"%(sym))

	if maskfile:
		if type(maskfile) is types.StringType: mask3D = get_image(maskfile)
		else:                                  mask3D = maskfile
	else: mask3D = model_cylinder(rmax, nx, ny, nz)

	fscmask = mask3D
	if CTF:
		from reconstruction import recons3d_4nn_ctf_MPI
		from filter         import filt_ctf
	else:	 from reconstruction import recons3d_4nn_MPI

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
	print  " READ IMAGES ", myid,nima,nproc

	rise = int(dp/pixel_size)

	data_nx = data[0].get_xsize()
	data_ny = data[0].get_ysize()
	if(nx < data_nx):
		data_nx = nx
		for im in xrange(nima):  data[im]=Util.window(data[im], data_nx, data_ny, 1, 0, 0, 0)
	data_nn = max(data_nx, data_ny)
	mask2D  = pad(model_blank(2*int(rmax), data_ny, 1, 1.0), data_nx, data_ny, 1, 0.0)
	fdata = [None]*nima
	for im in xrange(nima):
		data[im].set_attr('ID', list_of_particles[im])
		st = Util.infomask(data[im], mask2D, False)
		data[im] -= st[0]
		data[im] /= st[1]
		if CTF:
			qctf = data[im].get_attr("ctf_applied")
			if qctf == 0:
				ctf_params = data[im].get_attr("ctf")
				data[im] = filt_ctf(data[im], ctf_params)
				data[im].set_attr('ctf_applied', 1)
			elif qctf != 1:
				ERROR('Incorrectly set qctf flag', "ehelix_MPI", 1,myid)
		#  if FindPsi,  apply the angle to data[im], do fft and put in fdata[im]
		if FindPsi:
			phi,theta,psi,tsx,tsy = get_params_proj(data[im])
			if( theta != 0.0):
				if(abs(psi - 90.) < abs(psi - 270.0)):  gamma =  90.0
				else:                                   gamma = 270.0
				fdata[im] = fft( rot2pad(data[im], gamma-psi) )
			else:
				set_params_proj(data[im], [0.0, 90.0, 90.0, 0.0,0.0])
				fdata[im] = fft( data[im] )				
		else:
			set_params_proj(data[im], [0.0, 90.0, 90.0, 0.0,0.0])
			fdata[im] = fft( data[im] )
	del list_of_particles

	if debug:
		finfo.write( '%d loaded  \n' % nima )
		finfo.flush()
	if myid == main_node:
		# initialize data for the reference preparation function
		# No centering for helical reconstruction
		ref_data = [None, mask3D, None, None, None ]

	phiwobble = int(float(ywobble)/rise*dphi/delta+0.5)  # phiwobble is NOT in degrees, it is in nphi units

	nwx = 2*search_rng+3
	nwy = rise+2*ywobble+2
	nwxc = nwx//2
	nwyc = nwy//2
	nphi = int(360.0/delta + 0.5)
	#print  "  params  ",nwx,nwy,nwxc,nwyc,nphi
	if FindPsi:
		mode = "F"
		cnx = data_nn//2+1
		cny = cnx
		numr = Numrinit(1, data_nn//2-2, 1, mode)
		wr   = ringwe(numr, mode)
		maxrin = numr[len(numr)-1]
		crefim = [None]*nphi
	else:
		#  have to initialize them, otherwise problem with passing the arguments
		mode = "F"
		cnx = data_nx//2+1
		cny = cnx
		numr = []
		wr   = []
		maxrin = 0
		crefim = []

	terminate = 0
	Iter = 0
 	while Iter < max_iter and terminate == 0:
		Iter += 1
		if myid == main_node:
			start_time = time()
			print_msg("\nITERATION #%3d\n"%(Iter))
		volft, kbx, kby, kbz = prep_vol( vol )
		del vol

		refproj = [None]*nphi
		if( not Dsym):  rotproj = [None]*nphi
		else:           rotproj = []

		#refproj = EMData.read_images("refstraight45.hdf")   #   Here
		for iphi in xrange(nphi):
			refproj[iphi] = Util.window(  prgs( volft, kbz, [delta*iphi, 90.0, 90.0, 0.0, 0.0], kbx, kby ), data_nx, nz,1, 0,0,0)
			#refproj[iphi].write_image("refstraight45.hdf",iphi)
			st = Util.infomask(refproj[iphi] , mask2D, True)
			refproj[iphi] -= st[0]
			refproj[iphi] /= st[1]
			refproj[iphi] = Util.muln_img(refproj[iphi], mask2D )

			if FindPsi:
				temp = Util.Polar2Dm(pad(refproj[iphi], data_nn, data_nn, 1, "circumference"), cnx, cny, numr, mode)
				Util.Frngs(temp, numr)
				Util.Applyws(temp, numr, wr)
				crefim[iphi] = temp
			#  rotated in-plane by 180 are equivalent to rot_shift3D(vol,-90,180.0,90) with phi running as phi
			if(not Dsym):  rotproj[iphi] = fft( rot2pad(refproj[iphi],180.0) )
			refproj[iphi] = fft( refproj[iphi] )
		#exit()
		#if myid == main_node:  
		astart_time = time()
		for ifil in xrange(nfils):
			if myid == main_node:  start_time = time()
			if myid == main_node:
				print_msg("Process filament %4d %d\n"%(ifil,time()-start_time));start_time = time()
			ldata = [pad(data[im], data_nn, data_nn, 1, "circumference") for im in xrange(indcs[ifil][0],indcs[ifil][1])]
			Util.constrained_helix(ldata, fdata[indcs[ifil][0]:indcs[ifil][1]], refproj, rotproj, [float(dp), float(dphi), float(rise), float(delta)], [int(nphi), int(phiwobble), int(range), int(ywobble), int(Dsym), int(nwx), int(nwy), int(nwxc), int(nwyc)], FindPsi, float(psi_max), crefim, numr, int(maxrin), mode, int(cnx), int(cny))
			#constrained_helix     (ldata, fdata[indcs[ifil][0]:indcs[ifil][1]], refproj, rotproj,  dp, dphi, rise, delta ,  nphi, phiwobble, range, ywobble, Dsym, nwx, nwy, nwxc, nwyc , FindPsi, psi_max, crefim, numr, maxrin, mode, cnx, cny, myid, main_node)
			#constrained_helix     (data[indcs[ifil][0]:indcs[ifil][1]], fdata[indcs[ifil][0]:indcs[ifil][1]], refproj, rotproj,  dp, dphi, rise, delta ,  nphi, phiwobble, range, ywobble, Dsym, nwx, nwy, nwxc, nwyc , FindPsi, psi_max, crefim, numr, maxrin, mode, cnx, cny, myid, main_node)
			for im in xrange(indcs[ifil][0], indcs[ifil][1]):
				temp = Util.get_transform_params(ldata[im-indcs[ifil][0]], "xform.projection", "spider")
				set_params_proj(data[im],[temp["phi"],temp["theta"],temp["psi"],-temp["tx"],-temp["ty"]])
			if FindPsi:
				for im in xrange(indcs[ifil][0], indcs[ifil][1]):  fdata[im] = rot2pad(data[im], ldata[im-indcs[ifil][0]].get_attr("bestang"))
			#print  "Parameters computed for filament",myid,ifil,time()-start_time;start_time = time()
			if myid == main_node:
				print_msg("Parameters computed for filament %4d %d\n"%(ifil,time()-start_time));start_time = time()
		del ldata
		del refproj, volft
		if(not Dsym):  del rotproj
		print  "Time of alignment = ",myid,time()-astart_time
		mpi_barrier(MPI_COMM_WORLD)
		#if myid == main_node:
		#	print_msg("Time of alignment = %\n"%(time()-astart_time));start_time = time()

		if CTF:  vol = recons3d_4nn_ctf_MPI(myid, data, symmetry=sym, snr = snr, npad = npad, xysize = nx, zsize = nz)
		else:    vol = recons3d_4nn_MPI(myid, data, symmetry=sym, npad = npad, xysize = nx, zsize = nz)

		if myid == main_node:
			print_msg("3D reconstruction time = %d\n"%(time()-start_time));start_time = time()

		if myid == main_node:
			#vol.write_image(os.path.join(outdir, "vol%03d.hdf"%Iter))
			vol = vol.helicise(pixel_size, dp, dphi, fract, rmax, rmin)
			vol = sym_vol(vol, symmetry=sym)
			ref_data[0] = vol
			vol = user_func(ref_data)
			vol = vol.helicise(pixel_size, dp, dphi, fract, rmax, rmin)
			vol = sym_vol(vol, symmetry=sym)
			vol.write_image(os.path.join(outdir, "volf%03d.hdf"%Iter))
			#if(Iter == max_iter-1):  drop_image(vol, os.path.join(outdir, "volfshift.hdf"))

		bcast_EMData_to_all(vol, myid, main_node)
		# write out headers, under MPI writing has to be done sequentially
		mpi_barrier(MPI_COMM_WORLD)
		par_str = ['xform.projection', 'ID']
		if myid == main_node:
			start_time = time()
			if(file_type(stack) == "bdb"):
				from utilities import recv_attr_dict_bdb
				recv_attr_dict_bdb(main_node, stack, data, par_str, 0, nima, nproc)
			else:
				from utilities import recv_attr_dict
				recv_attr_dict(main_node, stack, data, par_str, 0, nima, nproc)
			print_msg("Time to write header information= %d\n"%(time()-start_time))
			start_time = time()
		else:		send_attr_dict(main_node, data, par_str, 0, nima)
		if myid == main_node:
			# write params to text file
			header(stack, params='xform.projection', fexport=os.path.join(outdir, "parameters%04d.txt"%Iter))
	if myid == main_node: print_end_msg("ehelix_MPI")
'''


def constrained_helix(data, fdata, refproj, rotproj, dp, dphi, rise, delta, nphi, symrestrict, phiwobble, rng, ywobble, ystep,\
					 Dsym, nwx, nwy, nwxc, nwyc, FindPsi, psi_max, crefim, numr, maxrin, mode, cnx, cny, myid, main_node):
	from alignment       import ringwe, ang_n
	from fundamentals    import ccf, fft, rot_shift2D
	from utilities       import get_params_proj, set_params_proj, file_type, compose_transform2, get_dist
	from math            import sqrt
	from time import time
	#  Here rise is in PIXELS

	start_time = time()
	ndata = len(data)
	#print  " WINDOW  ",nwx,nwy,nwxc,nwyc
	coords = [None]*ndata
	ccfs = [[None for i in xrange(nphi)] for j in xrange(ndata)]
	if(not Dsym):  ccfr = [[None for i in xrange(nphi)] for j in xrange(ndata)]
	for im in xrange(ndata):
		coords[im] = data[im].get_attr('ptcl_source_coord')
		for iphi in xrange(nphi):
			ccfs[im][iphi] = Util.window(ccf( refproj[iphi], fdata[im] ), nwx, nwy)
			if(not Dsym):  ccfr[im][iphi] = Util.window(ccf( rotproj[iphi], fdata[im] ), nwx, nwy)
	if myid == main_node: print " ccfs computed     ",time()-start_time; start_time = time()

	dxshiftlocal = [0.0]*ndata
	dyshiftlocal = [0.0]*ndata
	philocal     = [0.0]*ndata
	dphilocal    = [0.0]*ndata
	dirma = -1.0e23
	for idir in xrange(-1,2,2):
		xshiftlocal  = [0.0]*ndata
		if(not Dsym): xrshiftlocal  = [0.0]*ndata
		mxshiftlocal = [0.0]*ndata
		yshiftlocal  = [0.0]*ndata
		if(not Dsym): yrshiftlocal  = [0.0]*ndata
		myshiftlocal = [0.0]*ndata
		philocal  = [0.0]*ndata
		if(not Dsym): phirlocal  = [0.0]*ndata
		mphilocal = [0.0]*ndata
		tmax = -1.0e23
		for ix in xrange(1,nwx-1):                                         #  X shift
			#print "im: ", len(ccfs), ix,time()-start_time
			six = ix - nwxc
			siy = -rise/2.0 - ystep +ywobble
			#print "will do first while",rise,siy,ystep,ywobble
			while( siy < rise/2.0 - ystep - ywobble ):
				siy += ystep
				#for iy in xrange(1+ywobble,nwy-ywobble-1):                     #  Y shift
				#siy = iy - nwyc
				yshiftlocal[0]  = siy
				yrshiftlocal[0] = siy
				yiy = siy + nwyc
				fiy = int(yiy)
				ydif = yiy - fiy
				yrem = 1.0 - ydif
				iuy = fiy
				
				xrem = 1.0
				xdif = 0.0

				for iphi in xrange(nphi/symrestrict):#(102,103):#(58,59):#102,103):#nphi):                                  #  phi search
					qphi = iphi*delta
					philocal[0]  = qphi
					phirlocal[0] = (180.0 - qphi)%360.0
					# we use the first segment as a reference, so there is no interpolation, just copy the correlation
					#  Select largest correlation within +/- range pixels from the location we explore
					mxm = -1.0e23
					for iux in xrange(max(1, ix - rng), min(nwx - 1, ix+rng+1), 1):        #  X wobble
						#qcf = ccfs[0][iphi].get_value_at(iux,iy)
						qcf = xrem*yrem*ccfs[0][iphi].get_value_at(iux,iuy)   \
							+ xdif*yrem*ccfs[0][iphi].get_value_at(iux+1,iuy) \
							+ xrem*ydif*ccfs[0][iphi].get_value_at(iux,iuy+1) \
							+ xdif*ydif*ccfs[0][iphi].get_value_at(iux+1,iuy+1)
						if(qcf > mxm):
							mxm = qcf
							xshiftlocal[0] = float(iux-nwxc)
					#mxm = ccfs[0][iphi].get_value_at(ix,iy)
					if(not Dsym):
						mxr = -1.0e23
						for iux in xrange(max(1, ix - rng), min(nwx - 1, ix+rng+1), 1):     # Xr wobble
							#qcf = ccfr[0][iphi].get_value_at(iux,iy)
							qcf = xrem*yrem*ccfr[0][iphi].get_value_at(iux,iuy)   \
								+ xdif*yrem*ccfr[0][iphi].get_value_at(iux+1,iuy) \
								+ xrem*ydif*ccfr[0][iphi].get_value_at(iux,iuy+1) \
								+ xdif*ydif*ccfr[0][iphi].get_value_at(iux+1,iuy+1)
							if(qcf > mxr):
								mxr = qcf
								xrshiftlocal[0] = float(iux-nwxc)
						#mxr = ccfr[0][iphi].get_value_at(ix,iy)
					#print  ix,six,iy,siy,iphi,philocal[0],phirlocal[0],mxm,mxr,xshiftlocal[0],xrshiftlocal[0],yshiftlocal[0],yrshiftlocal[0]
					#print "  QTS    %4d   %12.6e  %12.6e"%(iphi,mxm,mxr)
					for im in xrange(1,ndata):                                                    #  predicted locations
						#print "im: ", len(ccfs), im,time()-start_time
						# dst is distance between segment 0 and current segment in pixels
						dst = get_dist(coords[im], coords[0])
						# predict for all remaining segments assuming number 0 has parameters (qphi, six, siy)
						# 
						qd  = round((siy + dst)/rise)
						piy = siy + dst - rise*qd                                                   #  predicted y shift
						#print " PREDICTED SHIFTS " ,six,siy,dst,rise,qd,piy
																									#  predicted phi with full angular accuracy, not an integer
						pphi = ((philocal[0] + idir*dphi*qd))%360.0
						pix  = six                                                                  #  predicted x shift
						#print " PREDICTIONS ",im,six,siy,qphi,dst,pix,piy,pphi


						xix = pix + nwxc
						#  Local x search
						fix = int(xix)
						xdif = xix - fix
						xrem = 1.0 - xdif


						ciq = -1.0e23
						# interpolate correlation at pphi
						ttphi = int(pphi/delta + 0.5)%nphi
						#print " LOOP ",phiwobble,fix,rng,fiy,ywobble
						for lphi in xrange(-phiwobble,phiwobble+1):                                         #  phi wobble
							tphi = (ttphi+lphi)%nphi
							for iux in xrange(max(1, fix - rng), min(nwx - 1, fix+rng+1), 1):             #  X wobble
								#for iuy in xrange(max(1, fiy - ywobble), min(nwy - 1, fiy+ywobble+1), 1):   #  Y wobble
								tiy = -ywobble - ystep + piy + nwyc
								#print "will do second while",rise,piy,tiy,ystep,ywobble
								while( tiy <= ywobble - ystep + piy + nwyc ):
									tiy += ystep
									iuy = int(tiy)
									#print " will do tiy", iuy,tiy,rise,ystep,piy,nwyc,nwy

									if(iuy > -1 and iuy < nwy-1):
										#print " IN iuy ",iuy
										ydif = yiy - fiy
										yrem = 1.0 - ydif

										qcf = xrem*yrem*ccfs[im][tphi].get_value_at(iux,iuy)   \
											+ xdif*yrem*ccfs[im][tphi].get_value_at(iux+1,iuy) \
											+ xrem*ydif*ccfs[im][tphi].get_value_at(iux,iuy+1) \
											+ xdif*ydif*ccfs[im][tphi].get_value_at(iux+1,iuy+1)
										#print  " INTERPOL ",lphi,tphi,iux,iuy,fiy,ywobble, nwy, qcf,xrem,xdif,yrem,ydif
										if(qcf > ciq):
											ciq = qcf
											xshiftlocal[im] = iux + xdif - nwxc
											yshiftlocal[im] = tiy - nwyc
											philocal[im]    = tphi*delta
							#print  "straight   ",im,xshiftlocal[im], yshiftlocal[im], philocal[im],six,siy,qphi, "    ",pix,piy,xix,yiy,pphi,   "    ",fix,fiy,tphi
						#ciq = xrem*yrem*ccfs[im][tphi].get_value_at(fix,fiy) + xdif*yrem*ccfs[im][tphi].get_value_at(fix+1,fiy) + xrem*ydif*ccfs[im][tphi].get_value_at(fix,fiy+1) + xdif*ydif*ccfs[im][tphi].get_value_at(fix+1,fiy+1) 
						#print " PREDICTIONS ",im,six,siy,qphi,dst,pix,piy,pphi
						#print "  S    %4d   %12.6e  %12.6e  %12.6e"%(iphi,mxm,ciq,mxm+ciq)
						#print  " ttphi,phiwobble,ywobble,fiy ",im, ttphi,phiwobble,nwy,fiy,ywobble,nwx,fix,rng
						#from  sys import exit
						#exit()
						mxm += ciq
						# now for rotated
						if(not Dsym):
							# 																										#  predicted phi for rotated 180 defs with full angular accuracy, not an integer
							pphi = ((phirlocal[0] + idir*dphi*qd))%360.0
							ciq = -1.0e23
							# interpolate correlation at pphi
							for lphi in xrange(-phiwobble,phiwobble+1):                                           #  phi wobble
								ttphi = (pphi + lphi*delta)%360.0
								tphi =  int(((180.-ttphi)%360.0)/delta+0.5)%nphi
								for iux in xrange(max(1, fix - rng), min(nwx - 1, fix+rng+1), 1):               #  X wobble
									#for iuy in xrange(max(1, fiy - ywobble), min(nwy - 1, fiy+ywobble+1), 1):     #  Y wobble
									
									tiy = -ywobble - ystep + piy + nwyc 
									while( tiy <= ywobble - ystep  + piy + nwyc ):
										tiy += ystep

										iuy = int(tiy)

										if(iuy > -1 and iuy < nwy-1):

											ydif = yiy - fiy
											yrem = 1.0 - ydif

											qcf = xrem*yrem*ccfr[im][tphi].get_value_at(iux,iuy)   \
												+ xdif*yrem*ccfr[im][tphi].get_value_at(iux+1,iuy) \
												+ xrem*ydif*ccfr[im][tphi].get_value_at(iux,iuy+1) \
												+ xdif*ydif*ccfr[im][tphi].get_value_at(iux+1,iuy+1)
											if(qcf > ciq):
												ciq = qcf
												xrshiftlocal[im] = iux + xdif - nwxc
												yrshiftlocal[im] = tiy - nwyc
												phirlocal[im]    = int(ttphi/delta + 0.5)*delta
												#print  "rotated ",pphi,pix,piy,xix,yiy,  xrshiftlocal[im],yrshiftlocal[im], phirlocal[im],"    ",fix,fiy,tphi,iux + xdif - nwxc,ciq

							#ciq = xrem*yrem*ccfr[im][tphi].get_value_at(fix,fiy) + xdif*yrem*ccfr[im][tphi].get_value_at(fix+1,fiy) + xrem*ydif*ccfr[im][tphi].get_value_at(fix,fiy+1) + xdif*ydif*ccfr[im][tphi].get_value_at(fix+1,fiy+1) 
							#print "  R  ",iphi,mxr,ciq,mxr+ciq
							mxr += ciq
						else:
							mxr = mxm-1.e5
					# The parameters are stored only for the first segment, the remaining one will have to be recomputed
					#if myid == main_node:  print  "maxima",mxm,mxr,tmax
					#if Iter == 1:  mxm = 2*mxr
					#mxm = -mxr
					if( mxr > mxm ):
						if(mxr > tmax):
							tmax = mxr
							mpsi = 270.0
							for im in xrange(ndata):  mxshiftlocal[im] = xrshiftlocal[im]
							for im in xrange(ndata):  myshiftlocal[im] = yrshiftlocal[im]
							for im in xrange(ndata):  mphilocal[im]    = (180.0-phirlocal[im])%360.0
					else:
						if(mxm > tmax):
							tmax = mxm
							mpsi = 90.0
							for im in xrange(ndata):  mxshiftlocal[im] = xshiftlocal[im]
							for im in xrange(ndata):  myshiftlocal[im] = yshiftlocal[im]
							for im in xrange(ndata):  mphilocal[im]    = philocal[im]
		#print " idir  ", idir,tmax
		if(tmax>dirma):
			dirma = tmax
			dpsi = mpsi
			for im in xrange(ndata):  dxshiftlocal[im] = mxshiftlocal[im]
			for im in xrange(ndata):  dyshiftlocal[im] = myshiftlocal[im]
			for im in xrange(ndata):  dphilocal[im]    = mphilocal[im]

	#print "  PARAMETERS FOR  0 ",mphi, 90.0, mpsi,mxm,mxr
	#print  "  max found ", tmax
	for im in xrange(ndata):
		# Do the prediction using set (mphi, theta=90., mpsi, msx, msy)
		#cim = data[im].get_attr('ptcl_source_coord')
		#dst = sqrt((c0[0] - cim[0])**2 + (c0[1] - cim[1])**2)
		psx  = dxshiftlocal[im]
		psy  = dyshiftlocal[im]
		pphi = dphilocal[im]
		#print "  PARAMETERS FOR IM ",im,pphi, 90.0, mpsi, psx, psy
		if FindPsi:
			iphi = int(pphi/delta + 0.5)%nphi
			#print  " ref number and current parameters reduced to 2D  ",iphi,0.0, psx, psy				
			#  I should only care what the previous residual angle was
			ophi, otheta, opsi3, opx3, opy3 = get_params_proj(data[im])
			#print " old 3D params in data ",ophi, otheta, opsi3, opx3, opy3
			if(abs(opsi3 - 90.) < abs(opsi3 - 270.0)):  gamma =  90.0
			else:                                       gamma = 270.0
			oalpha, osx, osy, junk = compose_transform2(0, opx3, opy3, 1.0, gamma-opsi3, 0, 0, 1.0) # reduce 3D to 2D
			#print " old 3D params, -> 2D ",oalpha, osx, osy
			# combine previous with the current in plane
			#print " current 2D combined with old 2D rotation",oalpha, csx, csy
			#  Find what the shift is without the angle
			junk, nnsx, nnsy, junk = compose_transform2(0.0, psx, psy, 1.0, -oalpha, 0., 0., 1.0)
			#print " 2D shift without angle ",nnsx, nnsy

			#rot_shift2D(data[im], 0.0, nnsx, nnsy).write_image("shifted.hdf",im)
			#fft(refproj[iphi]).write_image("reference.hdf",im)
			#from utilities import info

			cimage = Util.Polar2Dm(data[im], cnx+nnsx, cny+nnsy, numr, mode)
			Util.Frngs(cimage, numr)
			temp = Util.Crosrng_msg_s( cimage, crefim[iphi], numr)

			#from utilities import write_text_file
			#write_text_file([temp[qqq] for qqq in xrange(maxrin)],"ccf1d.txt")
			#exit()

			ipr = int(psi_max*maxrin/360. + 0.5)
			if(dpsi == 270.0):  incpsi = maxrin//2
			else:               incpsi = 0
			qn = -1.0e23
			for ips in xrange(-ipr,ipr+1,1):
				tot = (ips + incpsi + maxrin)%maxrin
				tval = temp.get_value_at(tot)
				#print  ips,incpsi,tot,tval
				if(tval > qn):
					qn = tval
					bestang = ang_n(tot+1.0, mode, maxrin)
			#print " best angle ",bestang
			bestang = (bestang - (dpsi-90.0))%360.0
			#print " angle applied ",-bestang
			#rot_shift2D(data[im],-bestang).write_image("rotated.hdf",im)
			#fdata[im] = fft( rot_shift2D(data[im], -bestang) )
			#print  " New composed 3D  ",dpsi,bestang, nnsx, nnsy

			epsi = (bestang+dpsi)%360.0
			psx = nnsx; psy = nnsy
			#print  " New composed 3D  ",pphi, 90.0, epsi, psx, psy
			#exit()
		else:
			epsi = dpsi
			bestang = 0.0
		data[im].set_attr("bestang", 360.0-bestang)
		#print  "  %7.2f  %7.2f  %7.2f  %7.2f  %7.2f"%(pphi, 90.0, epsi, psx, psy)
		set_params_proj(data[im], [pphi, 90.0, epsi, psx, psy])


'''
def constrained_helix_SHC(data, fdata, refproj, rotproj, dp, dphi, rise, delta, nphi, phiwobble, rng, ywobble, Dsym, nwx, nwy, nwxc, nwyc, FindPsi, psi_max, crefim, numr, maxrin, mode, cnx, cny, myid, main_node):
	from alignment       import ringwe, ang_n
	from fundamentals    import ccf, fft, rot_shift2D
	from utilities       import get_params_proj, set_params_proj, file_type, compose_transform2
	from math            import sqrt
	from random          import shuffle
	from time import time
	
	#start_time = time()
	#ndata = len(data)

	c0 = data[0].get_attr('ptcl_source_coord')
	ccfs = [[None for i in xrange(nphi)] for j in xrange(ndata)]
	if(not Dsym):  ccfr = [[None for i in xrange(nphi)] for j in xrange(ndata)]
	"""
	for im in xrange(ndata):
		for iphi in xrange(nphi):
			ccfs[im][iphi] = Util.window(ccf( refproj[iphi], fdata[im] ), nwx, nwy)
			if(not Dsym):  ccfr[im][iphi] = Util.window(ccf( rotproj[iphi], fdata[im] ), nwx, nwy)
	if myid == main_node: print " ccfs computed     ",time()-start_time; start_time = time()
	"""

	dxshiftlocal = [0.0]*ndata
	dyshiftlocal = [0.0]*ndata
	philocal  = [0.0]*ndata
	dphilocal = [0.0]*ndata

	previousmax = data[0].get_attr("previousmax")
	#print previousmax, data[0].get_attr("filament")
	bailout = False  # this is a dreadful way of going around the problem of no goto statement in python
	tempidir = [-1,1]                                                        # SHC
	shuffle(tempidir)
	blah = 0
	for idir in tempidir:
		if bailout:  break
		xshiftlocal  = [0.0]*ndata
		if(not Dsym): xrshiftlocal  = [0.0]*ndata
		yshiftlocal  = [0.0]*ndata
		if(not Dsym): yrshiftlocal  = [0.0]*ndata
		philocal  = [0.0]*ndata
		if(not Dsym): phirlocal  = [0.0]*ndata

		tempix = range(1,nwx-1)                                              # SHC
		shuffle(tempix)
		for ix in tempix:                                                   #  X shift
			if bailout:  break
			#print "im: ", len(ccfs), ix,time()-start_time
			six = ix - nwxc
			tempiy = range(1+ywobble,nwy-ywobble-1)  #  SHC
			shuffle(tempiy)
			for iy in tempiy:                                                #  Y shift
				if bailout:  break
				siy = iy - nwyc
				yshiftlocal[0]  = float(iy-nwyc)
				yrshiftlocal[0] = float(iy-nwyc)
				tempnphi = range(nphi)                                        #  SHC
				shuffle(tempnphi)
				for iphi in tempnphi:                                         #  phi search
					if bailout:  break
					#qphi = iphi*delta
					philocal[0]  = iphi*delta
					phirlocal[0] = (180.0 - iphi*delta)%360.0
					# we use the first segment as a reference, so there is no interpolation, just copy the correlation
					#  Select largest correlation within +/- rng pixels from the location we explore
					mxm = -1.0e23
					for iux in xrange(max(1, ix - rng), min(nwx - 1, ix+rng+1), 1):                           #  X wobble
						if(ccfs[0][iphi] == None):  ccfs[0][iphi] = Util.window(ccf( refproj[iphi], fdata[0] ), nwx, nwy)
						qcf = ccfs[0][iphi].get_value_at(iux,iy)
						if(qcf > mxm):
							mxm = qcf
							xshiftlocal[0] = float(iux-nwxc)
					#mxm = ccfs[0][iphi].get_value_at(ix,iy)
					if(not Dsym):
						mxr = -1.0e23
						for iux in xrange(max(1, ix - rng), min(nwx - 1, ix+rng+1), 1):                           # X wobble
							if(ccfr[0][iphi] == None):  ccfr[0][iphi] = Util.window(ccf( rotproj[iphi], fdata[0] ), nwx, nwy)
							qcf = ccfr[0][iphi].get_value_at(iux,iy)
							if(qcf > mxr):
								mxr = qcf
								xrshiftlocal[0] = float(iux-nwxc)
						#mxr = ccfr[0][iphi].get_value_at(ix,iy)
					#print  ix,six,iy,siy,iphi,philocal[0],phirlocal[0],mxm,mxr,xshiftlocal[0],xrshiftlocal[0],yshiftlocal[0],yrshiftlocal[0]

					tempupdown = range(2)                                      #  SHC
					shuffle(tempupdown)                                        #  What follows is an awkward construct to randomize the order in which up and down searches are done
					for updown in tempupdown:
						if bailout:  break
						blah += 1
						if(updown == 0):
							#  Straight
							for im in xrange(1,ndata):                                    #  predicted locations
								#print "im: ", len(ccfs), im,time()-start_time
								# dst is distance between segment 0 and current segment in pixels
								cim = data[im].get_attr('ptcl_source_coord')
								dst = sqrt((c0[0] - cim[0])**2 + (c0[1] - cim[1])**2)
								#dst = 15.0
								#print im,dst,rise
								# predict for all remaining segments assuming number 0
								#  has parameters (qphi, six, siy)
								# Assume for now inter-segment distances are multiples of rise -- jia
								pphi = (philocal[0] + idir*(dst/rise)*dphi)%360.0                          #  predicted phi with full angular accuracy, not an integer
								pix = six # predicted x shift
								piy = siy #  predicted y shift
								xix = pix + nwxc
								yiy = piy + nwyc
								#  Local x search
								fix = int(xix)
								xdif = xix - fix
								xrem = 1.0 - xdif
								fiy = int(yiy)
								ydif = yiy - fiy
								yrem = 1.0 - ydif
								ciq = -1.0e23
								# interpolate correlation at pphi
								ttphi = int(pphi/delta + 0.5)%nphi
								for lphi in xrange(-phiwobble,phiwobble+1):                                                                  #  phi wobble
									tphi = (ttphi+lphi)%nphi
									if(ccfs[im][tphi] == None):  ccfs[im][tphi] = Util.window(ccf( refproj[tphi], fdata[im] ), nwx, nwy)
									for iux in xrange(max(1, fix - rng), min(nwx - 1, fix+rng+1), 1):                                        #  X wobble
										for iuy in xrange(max(1, fiy - ywobble), min(nwy - 1, fiy+ywobble+1), 1):                            #  Y wobble											
											qcf = xrem*yrem*ccfs[im][tphi].get_value_at(iux,iuy) + xdif*yrem*ccfs[im][tphi].get_value_at(iux+1,iuy) + xrem*ydif*ccfs[im][tphi].get_value_at(iux,iuy+1) + xdif*ydif*ccfs[im][tphi].get_value_at(iux+1,iuy+1)
											if(qcf > ciq):
												ciq = qcf
												xshiftlocal[im] = iux + xdif - nwxc
												yshiftlocal[im] = iuy + ydif - nwyc
												philocal[im] = tphi*delta
									#print  "straight ",ix,iy,iphi, "    ", six,siy,qphi, "    ",pix,piy,xix,yiy,pphi,   "    ",fix,fiy,tphi
								#ciq = xrem*yrem*ccfs[im][tphi].get_value_at(fix,fiy) + xdif*yrem*ccfs[im][tphi].get_value_at(fix+1,fiy) + xrem*ydif*ccfs[im][tphi].get_value_at(fix,fiy+1) + xdif*ydif*ccfs[im][tphi].get_value_at(fix+1,fiy+1) 
								#print "  S  ",iphi,mxm,ciq,mxm+ciq
								mxm += ciq

							# The parameters are stored only for the first segment, the remaining ones will have to be recomputed
							if(mxm > previousmax):
								previousmax = mxm
								dpsi = 90.0
								for im in xrange(ndata):  dxshiftlocal[im] = xshiftlocal[im]
								for im in xrange(ndata):  dyshiftlocal[im] = yshiftlocal[im]
								for im in xrange(ndata):  dphilocal[im]    = philocal[im]
								#  Bail out
								bailout = True


						else:
							# now for rotated
							if(not Dsym):
								for im in xrange(1,ndata):                                    #  predicted locations
									#print "im: ", len(ccfs), im,time()-start_time
									# dst is distance between segment 0 and current segment in pixels
									cim = data[im].get_attr('ptcl_source_coord')
									dst = sqrt((c0[0] - cim[0])**2 + (c0[1] - cim[1])**2)
									#dst = 15.0
									#print im,dst,rise
									# predict for all remaining segments assuming number 0
									#  has parameters (qphi, six, siy)
									# Assume for now inter-segment distances are multiples of rise -- jia
									pphi = (phirlocal[0] + idir*(dst/rise)*dphi)%360.0 #  predicted phi for rotated 180 defs with full angular accuracy, not an integer
									pix = six # predicted x shift
									piy = siy #  predicted y shift
									xix = pix + nwxc
									yiy = piy + nwyc
									fix = int(xix)
									xdif = xix - fix
									xrem = 1.0 - xdif
									fiy = int(yiy)
									ydif = yiy - fiy
									yrem = 1.0 - ydif
									ciq = -1.0e23
									# interpolate correlation at pphi
									for lphi in xrange(-phiwobble,phiwobble+1):                                           #  phi wobble
										ttphi = (pphi + lphi*delta)%360.0
										tphi =  int(((180.-ttphi)%360.0)/delta+0.5)%nphi
										if(ccfr[im][tphi] == None):  ccfr[im][tphi] = Util.window(ccf( rotproj[tphi], fdata[im] ), nwx, nwy)
										for iux in xrange(max(1, fix - rng), min(nwx - 1, fix+rng+1), 1):             #  X wobble
											for iuy in xrange(max(1, fiy - ywobble), min(nwy - 1, fiy+ywobble+1), 1):     #  Y wobble
												qcf = xrem*yrem*ccfr[im][tphi].get_value_at(iux,iuy) + xdif*yrem*ccfr[im][tphi].get_value_at(iux+1,iuy) + xrem*ydif*ccfr[im][tphi].get_value_at(iux,iuy+1) + xdif*ydif*ccfr[im][tphi].get_value_at(iux+1,iuy+1)
												if(qcf > ciq):
													ciq = qcf
													xrshiftlocal[im] = iux + xdif - nwxc
													yrshiftlocal[im] = iuy + ydif - nwyc
													phirlocal[im]    = int(ttphi/delta + 0.5)*delta
													#print  "rotated ",pphi,pix,piy,xix,yiy,  xrshiftlocal[im],yrshiftlocal[im], phirlocal[im],"    ",fix,fiy,tphi,iux + xdif - nwxc,ciq
									#ciq = xrem*yrem*ccfr[im][tphi].get_value_at(fix,fiy) + xdif*yrem*ccfr[im][tphi].get_value_at(fix+1,fiy) + xrem*ydif*ccfr[im][tphi].get_value_at(fix,fiy+1) + xdif*ydif*ccfr[im][tphi].get_value_at(fix+1,fiy+1) 
									#print "  R  ",iphi,mxr,ciq,mxr+ciq
									mxr += ciq
								# The parameters are stored only for the first segment, the remaining ones will have to be recomputed
								if(mxr > previousmax):
									previousmax = mxr
									dpsi = 270.0
									for im in xrange(ndata):  dxshiftlocal[im] = xrshiftlocal[im]
									for im in xrange(ndata):  dyshiftlocal[im] = yrshiftlocal[im]
									for im in xrange(ndata):  dphilocal[im]    = (180.0-phirlocal[im])%360.0
									#  Bail out
									bailout = True

	#print "blah blah: ", blah
	# if got here, it did not find anything better and should return
	if not bailout: return   -1
	
	# if found better should end up here
	# print "found better ", idir,previousmax, data[0].get_attr('filament')


	#print "  PARAMETERS FOR  0 ",mphi, 90.0, mpsi,mxm,mxr
	#print  "  max found ", tmax
	for im in xrange(ndata):
		data[im].set_attr("previousmax", previousmax)
		# Do the prediction using set (mphi, theta=90., mpsi, msx, msy)
		#cim = data[im].get_attr('ptcl_source_coord')
		#dst = sqrt((c0[0] - cim[0])**2 + (c0[1] - cim[1])**2)
		psx  = dxshiftlocal[im]
		psy  = dyshiftlocal[im]
		pphi = dphilocal[im]
		#print "  PARAMETERS FOR IM ",im,pphi, 90.0, mpsi, psx, psy
		if FindPsi:
			iphi = int(pphi/delta + 0.5)%nphi
			#print  " ref number and current parameters reduced to 2D  ",iphi,0.0, psx, psy				
			#  I should only care what the previous residual angle was
			ophi, otheta, opsi3, opx3, opy3 = get_params_proj(data[im])
			#print " old 3D params in data ",ophi, otheta, opsi3, opx3, opy3
			if(abs(opsi3 - 90.) < abs(opsi3 - 270.0)):  gamma =  90.0
			else:                                       gamma = 270.0
			oalpha, osx, osy, junk = compose_transform2(0, opx3, opy3, 1.0, gamma-opsi3, 0, 0, 1.0) # reduce 3D to 2D
			#print " old 3D params, -> 2D ",oalpha, osx, osy
			# combine previous with the current in plane
			#print " current 2D combined with old 2D rotation",oalpha, csx, csy
			#  Find what the shift is without the angle
			junk, nnsx, nnsy, junk = compose_transform2(0.0, psx, psy, 1.0, -oalpha, 0., 0., 1.0)
			#print " 2D shift without angle ",nnsx, nnsy

			#rot_shift2D(data[im], 0.0, nnsx, nnsy).write_image("shifted.hdf",im)
			#fft(refproj[iphi]).write_image("reference.hdf",im)
			#from utilities import info

			cimage = Util.Polar2Dm(data[im], cnx+nnsx, cny+nnsy, numr, mode)
			Util.Frngs(cimage, numr)
			temp = Util.Crosrng_msg_s( cimage, crefim[iphi], numr)

			#from utilities import write_text_file
			#write_text_file([temp[qqq] for qqq in xrange(maxrin)],"ccf1d.txt")
			#exit()

			ipr = int(psi_max*maxrin/360. + 0.5)
			if(dpsi == 270.0):  incpsi = maxrin//2
			else:               incpsi = 0
			qn = -1.0e23
			for ips in xrange(-ipr,ipr+1,1):
				tot = (ips + incpsi + maxrin)%maxrin
				tval = temp.get_value_at(tot)
				#print  ips,incpsi,tot,tval
				if(tval > qn):
					qn = tval
					bestang = ang_n(tot+1.0, mode, maxrin)
			#print " best angle ",bestang
			bestang = (bestang - (dpsi-90.0))%360.0
			#print " angle applied ",-bestang
			#rot_shift2D(data[im],-bestang).write_image("rotated.hdf",im)
			#fdata[im] = fft( rot_shift2D(data[im], -bestang) )
			#print  " New composed 3D  ",dpsi,bestang, nnsx, nnsy

			epsi = (bestang+dpsi)%360.0
			psx = nnsx; psy = nnsy
			#print  " New composed 3D  ",pphi, 90.0, epsi, psx, psy
			#exit()
		else:
			epsi = dpsi
			bestang = 0.0
		data[im].set_attr("bestang", 360.0-bestang)
		#print  "  %7.2f  %7.2f  %7.2f  %7.2f  %7.2f"%(pphi, 90.0, epsi, psx, psy)
		set_params_proj(data[im], [pphi, 90.0, epsi, psx, psy])
		#print get_params_proj(data[im])
	return  blah

def yshift_to_phi(yshift, psi, pixel_size, dp, dphi):
	
	# If psi=90: -1 rise of yshift corresponds to rotation by -1 dphi
	#			 1 rise of yshift corresponds to rotation by 1 dphi
	# IF psi=270, -1 rise of yshift corresponds to rotation by 1 dphi
	#			  1 rise of yshift corresponds to rotation by -1 dphi
	dpp = dp/pixel_size
	sgn_yshift = 1
	if yshift < 0:  sgn_yshift = -1
	mdp, yshift = divmod(abs(yshift), dpp)
	yshift = yshift * sgn_yshift

	adphi = mdp * dphi
	if abs(psi - 90.0) >= 90.0:  adphi *= -sgn_yshift
	else:                        adphi *= sgn_yshift

	# check if yshift is larger than 0.5*dpp
	if abs(yshift) > 0.5*dpp:
		sgn_yshift=1
		if yshift < 0:
			sgn_yshift=-1
			yshift = dpp - abs(yshift)
		else:
			yshift = -(dpp - abs(yshift))
		
		if yshift > 0.0:
			if abs(psi - 90.0) < 90.0:
				adphi += dphi
			else:
				adphi -= dphi
		if yshift < 0:
			if abs(psi - 90.0) < 90.0:
				adphi -= dphi
			else:
				adphi += dphi		
	return yshift, adphi
'''

"""
def generate_projections_helical(vol, stack, proj_stack, consparams, trueparams, dphi, dp, pixel_size, fil_attr='filament', phirand=0.0, psirand=0.0, yrand=0.0, xrand=0.0, N = -1, phizero=True, yzero=True, xzero=True, CTF=True, xpermitrange=0.0):

	'''
	
	The program writes to disk a stack of projections named proj_stack calculated from vol 
	using projection parameters predicted based on helical symmetry parameters.

	The projections are generated using predicted value plus a random perturbation based 
	on phirand/psirand/yrand/xrand
		
	The phi/psi/x/y of the projection orientation parameters of segment is set to the 
	predicted value.

	
	INPUT
	
	vol: 			Name of volume from which projections are calculated
	stack: 			Name of data stack containing filament and ctf information
	proj_stack: 	Name of projection stack that will be generated from predicted parameters
	consparams:		Name of text file to which the predicted UNPERTURBED parameters are saved to.
	dphi: 			azimuthal rotation of helical symmetry
	dp: 			rise of helical symmetry
	pixel_size: 	pixel size
	fil_attr: 		attribute under which filament identifier is stored in header.
	N:				number of filaments to generate segment projections for. 
					Default is N = -1, in which case all filaments are generated.
	
	phirand:		A random float in the range [-phirand, phirand] is added to the predicted
					phi of each segment.
	psirand:		A random float in the range [-psirand, psirand] is added to the predicted
					psi of each segment.
	xrand:			A random float in the range [-xrand, xrand] is added to the predicted
					x of ALL segments belonging to same filament. (PIXELS, need to change to angstroms at some point)
	yrand:			A random float in the range [-yrand, yrand] is added to the predicted
					y of ALL segments belonging to same filament. (PIXELS, need to change to angstroms at some point)
									
	phizero:		if True, then the phi of the segment used to predict orientation parameters of the remaining segments (based on helical symmetry) is set to zero.
			 		if False, then it's set to a random float between 0 and 360.
	yzero:   		if True, then the y shift of the segment used to predict orientation parameters of the remaining segments (based on helical symmetry) is set to zero.
			 		if False, then it's set to a random float in the range [-dpp/2, dpp/2].
	xpermitrage: 	the amount in PIXELS (need to change to Angstroms at some point) of the x-shift of the segment used to predict orientation parameters of the remaining segments.
	xzero:   		if True, then the x shift of the segment used to predict orientation parameters of the remaining segments (based on helical symmetry) is set to zero.
			 		if False, then it's set to a random float in the range [-xpermitrange, xpermitrange].
	
	Note: If phizero=0, yzero=0 and xzero=0, then the reconstructions of the filaments 
		  should be aligned, i.e., the reconstruction using ALL of the segments should 
		  yield vol (assuming there are enough segments generated in the first place).
	
	'''
	
	from random         import uniform, random
	from projection     import prep_vol, prgs
	from applications   import get_dist,ordersegments, header
	from filter			import filt_ctf
	from utilities		import get_im, write_text_row
	
	
	filaments = ordersegments(stack)
	ptclcoords = EMUtil.get_all_attributes(stack, 'ptcl_source_coord')

	if( N < 1 ): N = len(filaments)	
	if N > len(filaments):
		N = len(filaments)
	
	volft, kbx, kby, kbz = prep_vol( get_im(vol) )

	THETA = 90.
	PSI   = 90.

	rise = dp/pixel_size  # in pixels
	permitrange = rise/2.0
	ima = EMData()
	counter = 0
	projparams=[]
	for ifil in xrange(N):
		#print "ifil: ", ifil
		mic = filaments[ifil]
		nsegs = len(mic)
		
		phi0 = 0
		if not(phizero):
			phi0 = 360.0*random()
		y0 = 0
		if not(yzero):
			y0 = permitrange*(random()-0.5)
		x0 = 0
		if not(xzero):
			x0 = xpermitrange*(random()-0.5)
		
		refcoords = ptclcoords[mic[0]]
		
		# calculate perturbations for x and y shifts, which for now we assume is same for 
		# segments belonging to same filament
		pty = yrand*(random() - 0.5)
		ptx = xrand*(random() - 0.5)
		
		for iseg in xrange(nsegs):
			dA = pixel_size * get_dist(ptclcoords[mic[iseg]], refcoords) # distance in Angstroms between segments mic[iseg] and mic[0]

			iphi = (phi0 + (dA/dp)*dphi)%360.0
			iy   = y0 + ((dA%dp)/pixel_size) # in pixels
			if( iy > permitrange ): iy -= rise
			if( iy < -permitrange ): iy += rise
			
			ima.read_image(stack, mic[iseg], True)
			fifi = ima.get_attr(fil_attr)
			coord = ima.get_attr('ptcl_source_coord')
			
			t = Transform({'type':'spider','phi':iphi,'theta':THETA,'psi':PSI,'tx':x0,'ty':iy})
			torg = Transform(t)
			
			# add perturbations
			sp = t.get_params('spider')
			newphi = sp['phi'] + phirand *(random()-0.5)
			newpsi = sp['psi'] + psirand *(random()-0.5)
			newtx  = sp['tx']  + ptx
			newty  = sp['ty']  + pty
			
			t = Transform({'type':'spider','phi':newphi,'theta':sp['theta'],'psi':newpsi,'tx':newtx,'ty':newty})
			sp = t.get_params('spider')
			tmp = prgs(volft, kbz, [sp['phi'], sp['theta'],sp['psi'], sp['tx'], sp['ty']], kbx, kby)
			projparams.append([sp['phi'], sp['theta'],sp['psi'], -sp['tx'], -sp['ty']])
			
			if CTF:
				ct   = ima.get_attr('ctf')
				tmp = filt_ctf(tmp,ct)
			
			tmp.set_attr( "xform.projection" ,torg)
			tmp.set_attr(fil_attr, fifi)
			tmp.set_attr( "ptcl_source_coord" ,coord)
			tmp.write_image(proj_stack, counter)
			counter += 1
	
	header(proj_stack, params='xform.projection', fexport=consparams)
	write_text_row(projparams, trueparams)
	
def generate_projections_helical_xwobble(vol, stack, proj_stack, consparams, trueparams, set_to_proj_ori, dphi, dp, pixel_size, fil_attr='filament', phirand=0.0, psirand=0.0, NFIL = -1, NSEGS = -1, phizero=True, yzero=True, xzero=True, CTF=True, xpermitrange=0.0, xwobble_max=0.0, ywobble_max=0.0, PSI=90.0, startphi=0.0, startx = 0.0, starty=0.0):

	'''
	
	The program writes to disk a stack of projections named proj_stack calculated from vol 
	using projection parameters predicted based on helical symmetry parameters.

	The projections are generated using predicted value plus a random perturbation based 
	on phirand/psirand/yrand/xrand
		
	The phi/psi/x/y of the projection orientation parameters of segment is set to the 
	predicted value.

	
	INPUT
	
	vol: 			Name of volume from which projections are calculated
	stack: 			Name of data stack containing filament and ctf information
	proj_stack: 	Name of projection stack that will be generated from predicted parameters
	consparams:		Name of text file to which the predicted UNPERTURBED parameters are saved to.
	set_to_proj_ori: If True, then set xform.projection attribute of each projection to the 
					projection parameters used to calculate the projection.
					Otherwise, set it to predicted parameters.
	dphi: 			azimuthal rotation of helical symmetry
	dp: 			rise of helical symmetry
	pixel_size: 	pixel size
	fil_attr: 		attribute under which filament identifier is stored in header.
	NFIL:				number of filaments to generate segment projections for. 
					Default is N = -1, in which case all filaments are generated.
	NSEGS:			Generate projections from just enough filaments so the resulting number
					of projections is not less than NSEGS.
					Default is NSEGS=-1, in which case this parameters is set to the 
					number of segments in input stack.
	phirand:		phi pertubation: a random float in the range [-phirand, phirand] is added to the predicted
					phi of each segment. If the predicted phi is < 180 (> 180), then the perturbed phi is also
					guaranteed to be < 180 (>180). This is because in local search, if the initial
					parameters is mirrored, then only mirrored are checked, and similarly for not mirrored.
					So if predicted is mirrored while perturbed is not (or vice versa), then the local search
					will never find exact answers even if there is no noise.
	psirand:		psi perturbation: a random float in the range [-psirand, psirand] is added to the predicted
					psi of each segment.
									
	phizero:		if True, then the phi of the segment used to predict orientation parameters of the remaining segments (based on helical symmetry) is set to zero.
			 		if False, then it's set to a random float between 0 and 360.
	yzero:   		if True, then the y shift of the segment used to predict orientation 
					parameters of the remaining segments (based on helical symmetry) is set to zero.
			 		if False, then it's set to a random float in the range [-dpp/2, dpp/2].
	xpermitrange: 	when xzero is False, the x-shift (in Angstroms) of the 
					segment used to predict orientation parameters of the remaining 
					segments is set to a random float in the range [-xpermitrange, xpermitrange].
	xzero:   		if True, then the x shift of the segment used to predict orientation 
					parameters of the remaining segments (based on helical symmetry) is 
					set to zero.
					if False, then it's set to a random float in the range 
					[-xpermitrange, xpermitrange].
	xwobble_max:  	The maximum amount in Angstroms each segment can "wobble" along
					x-axis (is not the same for all segments windowed from same filament).
	ywobble_max:  	The maximum amount in Angstroms each segment can "wobble" along
					y-axis (is not the same for all segments windowed from same filament).
					
	Note: If phizero=0, yzero=0 and xzero=0, then the reconstructions of the filaments 
		  should be aligned, i.e., the reconstruction using ALL of the segments should 
		  yield vol (assuming there are enough segments generated in the first place).
	
	'''
	
	from random         import uniform, random
	from projection     import prep_vol, prgs
	from applications   import get_dist,ordersegments, header
	from filter			import filt_ctf
	from utilities		import get_im, write_text_row
	import sys
	
	xpermitrange = xpermitrange/pixel_size
	xwobble_max = xwobble_max/pixel_size
	ywobble_max = ywobble_max/pixel_size
	
	filaments = ordersegments(stack)
	ptclcoords = EMUtil.get_all_attributes(stack, 'ptcl_source_coord')

	if( NFIL < 1 ): NFIL = len(filaments)	
	if NFIL > len(filaments):
		NFIL = len(filaments)
	if NSEGS < 0:
		NSEGS = EMUtil.get_image_count(stack)
	volft, kbx, kby, kbz = prep_vol( get_im(vol) )

	THETA = 90.

	rise = dp/pixel_size  # in pixels
	permitrange = rise/2.0
	ima = EMData()
	counter = 0
	projparams=[]
	totsegs = 0
	for ifil in xrange(NFIL):
		if totsegs > NSEGS:
			break
		mic = filaments[ifil]
		nsegs = len(mic)
		
		totsegs += nsegs
		
		phi0 = startphi
		if not(phizero):
			phi0 = 360.0*random()
		y0 = starty
		if not(yzero):
			y0 = permitrange*2.0*(random()-0.5)
		x0 = startx
		if not(xzero):
			x0 = xpermitrange*2.0*(random()-0.5)
		
		refcoords = ptclcoords[mic[0]]
		
		# x and y perturbations
		xwob_prev = int( (random()-0.5) * 2.0 * xwobble_max)
		ywob_prev = int( (random()-0.5) * 2.0 * ywobble_max)
		
		for iseg in xrange(nsegs):
			dA = pixel_size * get_dist(ptclcoords[mic[iseg]], refcoords) # distance in Angstroms between segments mic[iseg] and mic[0]
			#print "segment ", iseg, " : ", (dA/dp)
			iphi = (phi0 + (dA/dp)*dphi)%360.0
			iy   = y0 + ((dA%dp)/pixel_size) # in pixels
			if( iy > permitrange ): iy -= rise
			if( iy < -permitrange ): iy += rise
			
			ima.read_image(stack, mic[iseg], True)
			fifi = ima.get_attr(fil_attr)
			coord = ima.get_attr('ptcl_source_coord')
			
			t = Transform({'type':'spider','phi':iphi,'theta':THETA,'psi':PSI,'tx':x0,'ty':iy})
			torg = Transform(t)
			
			# add perturbations
			sp = t.get_params('spider')
			phiptr = phirand * 2.0 * (random()-0.5)
			phi_perturb = (sp['phi'] + phiptr)%360.
			if (sp['phi'] < 180 and phi_perturb >= 180) or (sp['phi'] >= 180 and phi_perturb < 180):
				phi_perturb  = (sp['phi'] - phiptr)%360.
			if (sp['phi'] < 180 and phi_perturb >= 180) or (sp['phi'] >= 180 and phi_perturb < 180):
				print "phi perturbation too large, decrease it! ", sp['phi'], 
				sys.exit()
			newphi = phi_perturb
			newpsi = sp['psi'] + psirand * 2.0 * (random()-0.5)
			newtx  = sp['tx']
			newty  = sp['ty'] 
			
			# add x-wobble
			# current wobble should differ from previous wobble by at most one Angstrom
			xwob = (random()-0.5)*2.0 + xwob_prev
			
			if xwob < 0:
				xwob = int(xwob - 0.5)
			else:
				xwob = int(xwob + 0.5)
				
			if abs(xwob) > xwobble_max:
				if xwob < 0:
					xwob = -xwobble_max
				else:
					xwob = xwobble_max
			xwob_prev = xwob
			
			newtx += xwob
			
			# add y-wobble
			# current wobble should differ from previous wobble by at most one Angstrom
			ywob = (random()-0.5)*2.0 + ywob_prev
			
			if ywob < 0:
				ywob = int(ywob - 0.5)
			else:
				ywob = int(ywob + 0.5)
				
			if abs(ywob) > ywobble_max:
				if ywob < 0:
					ywob = -ywobble_max
				else:
					ywob = ywobble_max
			ywob_prev = ywob
			
			newty += ywob
			
			t = Transform({'type':'spider','phi':newphi,'theta':sp['theta'],'psi':newpsi,'tx':newtx,'ty':newty})
			sp = t.get_params('spider')
			tmp = prgs(volft, kbz, [sp['phi'], sp['theta'],sp['psi'], sp['tx'], sp['ty']], kbx, kby)
			projparams.append([sp['phi'], sp['theta'],sp['psi'], -sp['tx'], -sp['ty']])
			
			if CTF:
				ct   = ima.get_attr('ctf')
				tmp = filt_ctf(tmp,ct)
			
			if not(set_to_proj_ori):
				tmp.set_attr( "xform.projection" ,torg)
				
			tmp.set_attr(fil_attr, fifi)
			tmp.set_attr( "ptcl_source_coord" ,coord)
			tmp.write_image(proj_stack, counter)
			counter += 1
	
	header(proj_stack, params='xform.projection', fexport=consparams)
	write_text_row(projparams, trueparams)
"""

def newparams_3D_to_2D(stack):
	from utilities import print_begin_msg, print_end_msg, print_msg, set_params2D, write_header
	from development import newparams_3D_2D
	
	print_begin_msg("newparams_3D_to_2D")
	print_msg("Input stack                 : %s\n\n"%(stack))

	nima = EMUtil.get_image_count(stack)
	ima = EMData()
	for im in xrange(nima):
		ima.read_image(stack, im, True)
		from utilities import set_params_proj, get_params_proj
		phi,theta,psi,s2x,s2y = get_params_proj( ima )
		alpha, sx, sy, mirror = newparams_3D_2D(phi, theta, psi, s2x, s2y)
		set_params2D(ima, [alpha, sx, sy, mirror, 1.0])
		write_header(stack, ima, im)

def newparams_3D_2D(phi, theta, psi, s2x, s2y):
	"""
		Convert 3D alignment parameters (phi, theta, psi, s2x, s2y)  # there is no mirror in 3D! 
		into 2D alignment parameters (alpha, sx, sy, mirror)
	"""
	from utilities import compose_transform2, print_msg
	
	if theta > 90.0:
		mirror = 1
		print_msg(" theta  > 90 ... dealing with it same as theta=90 for now...\n")
		if abs(psi - 90) < 90: # psi ~ 90
			alpha, sx, sy, scalen = compose_transform2(0, s2x, s2y, 1.0, 90.0-psi, 0, 0, 1.0)
		else:
			alpha, sx, sy, scalen = compose_transform2(0, s2x, s2y, 1.0, 270.0-psi, 0, 0, 1.0)
	else:
		mirror = 0
		if abs(psi - 90) < 90: # psi ~ 90
			alpha, sx, sy, scalen = compose_transform2(0, s2x, s2y, 1.0, 90.0-psi, 0, 0, 1.0)
		else:
			alpha, sx, sy, scalen = compose_transform2(0, s2x, s2y, 1.0, 270.0-psi, 0, 0, 1.0)
	return  alpha, sx, sy, mirror

def applyxshift(stack, newstack):
	from utilities import get_im, get_params2D, set_params2D
	from fundamentals import cyclic_shift
	
	nima	=EMUtil.get_image_count(stack)
	for im in xrange(nima):
		prj = get_im(stack,im)
		alpha, sx, sy, mirror, scale = get_params2D(prj)
		prj = cyclic_shift(prj, int(sx))
		set_params2D(prj, [0.0,0.,0.0,0,1])
		prj.write_image(newstack, im)
		print "   ",im,
	print " "


def ali3d_shc(stack, ref_vol, outdir, maskfile = None, ir = 1, ou = -1, rs = 1, 
            xr = "4 2 2 1", yr = "-1", ts = "1 1 0.5 0.25", delta = "10 6 4 4", an = "-1", apsi = "-1", deltapsi = "-1", startpsi = "-1",
	    center = -1, maxit = 5, CTF = False, snr = 1.0,  ref_a = "S", sym = "c1",
	    user_func_name = "ref_ali3d", fourvar = True, npad = 4, debug = False, MPI = False, termprec = 0.0):
	"""
		Name
			ali3d - Perform 3-D projection matching given initial reference volume and image series
		Input
			stack: set of 2-D images in a stack file, images have to be squares
			ref_vol: initial reference volume
			outdir: directory name into which the results will be written
			maskfile: filename of the file containing 3D mask.
			ir: inner radius for rotational correlation > 0 
			ou: outer radius for rotational correlation <int(nx/2)-1 
			rs: steps between rings in rotational correlation >0
			xr: range for translation search in x direction in each iteration, search is +/xr
			yr: range for translation search in y direction in each iteration, search is +/yr
			ts: step size of the translation search in both directions, search is -xr, -xr+ts, 0, xr-ts, xr, can be fractional.
			delta: angular step for the reference projections in respective iterations
			an: angular neighborhood for local searches
			center: average center method
			max_iter: maximum iterations at each angle step
			CTF: if the flag is present, program will use the CTF information stored in file headers
			snr: signal noise ratio used in the 3D reconstruction
			ref_a: method for creating quasi-uniform distribution of the projection directions of reference projections: "S" - spiral
			sym: symmetry of the refined structure
			function: name of the user-supplied-function
			MPI: if presetm use MPI version
		Output
			output_directory: directory name into which the output files will be written.
			header: the alignment parameters are stored in the headers of input files as Transform Object xform.proj
	"""
	if MPI:
		ali3d_shcMPI(stack, ref_vol, outdir, maskfile, ir, ou, rs, xr, yr, ts,
	        	delta, an, apsi, deltapsi, startpsi, center, maxit, CTF, snr, ref_a, sym, user_func_name,
				fourvar, npad, debug, termprec)



# parameters: list of (all) projections | reference volume | ...
def ali3d_shc_simple_mpi(stack, ref_vol, maskfile = None, ir = 1, ou = -1, rs = 1, 
        xr = "4 2 2 1", yr = "-1", ts = "1 1 0.5 0.25", delta = "10 6 4 4", an = "-1", apsi = "-1", deltapsi = "-1", startpsi = "-1",
	    center = -1, maxit = 5, CTF = False, snr = 1.0,  ref_a = "S", sym = "c1",  user_func_name = "ref_ali3d",
	    fourvar = False, npad = 4, debug = False, termprec = 0.0, mpi_comm = None, log = None, proc_of_checked_refs=2 ):

	from alignment       import Numrinit, prepare_refrings, proj_ali_incore, proj_ali_incore_local, proj_ali_incore_local_psi
	from utilities       import model_circle, get_image, drop_image, get_input_from_string
	from utilities       import reduce_EMData_to_root, bcast_EMData_to_all
	from utilities       import send_attr_dict
	from utilities       import get_params_proj, file_type, set_params_proj
	from fundamentals    import rot_avg_image
	import os
	import types
	from utilities       import print_begin_msg, print_end_msg, print_msg
	from mpi             import mpi_bcast, mpi_comm_size, mpi_comm_rank, MPI_FLOAT, MPI_COMM_WORLD, mpi_barrier
	from mpi             import mpi_reduce, MPI_INT, MPI_SUM
	from filter          import filt_ctf
	from projection      import prep_vol, prgs
	from statistics      import hist_list, varf3d_MPI
	from applications    import MPI_start_end
	from math            import sqrt
	from random import random
	from utilities       import wrap_mpi_gatherv, wrap_mpi_bcast
	from reconstruction  import recons3d_4nn_MPI, recons3d_4nn_ctf_MPI
	from development     import shc

	if mpi_comm == None: mpi_comm = MPI_COMM_WORLD

	if log == None:
		from logger import Logger
		log = Logger()

	number_of_proc = mpi_comm_size(mpi_comm)
	myid           = mpi_comm_rank(mpi_comm)
	main_node = 0
	
	if myid == main_node:
		log.add("Start")
	
	xrng        = get_input_from_string(xr)
	if  yr == "-1":  yrng = xrng
	else          :  yrng = get_input_from_string(yr)
	step        = get_input_from_string(ts)
	delta       = get_input_from_string(delta)
	lstp = min(len(xrng), len(yrng), len(step), len(delta))
	if an == "-1":
		an = [-1] * lstp
	else:
		an = get_input_from_string(an)

	if apsi == "-1":
		apsi = [-1] * lstp
	else:
		apsi = get_input_from_string(apsi)

	if deltapsi == "-1":
		deltapsi = [-1] * lstp
	else:
		deltapsi = get_input_from_string(deltapsi)

	if startpsi == "-1":
		startpsi = [-1] * lstp
	else:
		startpsi = get_input_from_string(startpsi)

	first_ring  = int(ir)
	rstep       = int(rs)
	last_ring   = int(ou)
	max_iter    = int(maxit)
	center      = int(center)

	vol = ref_vol
	nx      = vol.get_xsize()
	if last_ring < 0:	last_ring = int(nx/2) - 2

	if myid == main_node:
		import user_functions
		user_func = user_functions.factory[user_func_name]

	if maskfile:
		mask3D = maskfile
	else: 
		mask3D = model_circle(last_ring, nx, nx, nx)

	numr	= Numrinit(first_ring, last_ring, rstep, "F")
	mask2D  = model_circle(last_ring,nx,nx) - model_circle(first_ring,nx,nx)

	fscmask = mask3D  #model_circle(last_ring,nx,nx,nx)  For a fancy mask circle would work better  PAP 7/21/11
	if CTF:
		from reconstruction import rec3D_MPI
		from filter         import filt_ctf
	else:	 from reconstruction import rec3D_MPI_noCTF

	if myid == main_node:
		
		# horatio active_refactoring Jy51i1EwmLD4tWZ9_00000_1
		# active = [ stack[im].get_attr('active') for im in xrange(len(stack)) ]
		# list_of_particles = []
		# for im in xrange(len(active)):
		# 	if active[im]:  list_of_particles.append(im)
		# del active
		# nima = len(list_of_particles)
		
		nima = len(stack)
		list_of_particles = range(nima)
			
		
		
		
	else:
		nima = 0
	total_nima = wrap_mpi_bcast(nima, main_node, mpi_comm)

	if myid != main_node:
		list_of_particles = [-1]*total_nima
	list_of_particles = wrap_mpi_bcast(list_of_particles, main_node, mpi_comm)

	image_start, image_end = MPI_start_end(total_nima, number_of_proc, myid)
	# create a list of images for each node
	list_of_particles = list_of_particles[image_start: image_end]
	nima = len(list_of_particles)

	data = [ stack[im] for im in list_of_particles ]
	if fourvar:  original_data = []
	for im in xrange(nima):
		data[im].set_attr('ID', list_of_particles[im])
		if fourvar: original_data.append(data[im].copy())
		ctf_applied = data[im].get_attr_default('ctf_applied', 0)
		if CTF and ctf_applied == 0:
			ctf_params = data[im].get_attr("ctf")
			st = Util.infomask(data[im], mask2D, False)
			data[im] -= st[0]
			data[im] = filt_ctf(data[im], ctf_params)
			data[im].set_attr('ctf_applied', 1)

	if myid == main_node:
		# initialize data for the reference preparation function
		ref_data = [ mask3D, max(center,0), None, None, None, None ]
		# for method -1, switch off centering in user function

	from time import time	

	#  this is needed for gathering of pixel errors
	disps = []
	recvcount = []
	for im in xrange(number_of_proc):
		if im == main_node :  disps.append(0)
		else:                 disps.append(disps[im-1] + recvcount[im-1])
		ib, ie = MPI_start_end(total_nima, number_of_proc, im)
		recvcount.append(ie-ib)

	pixer = [0.0]*nima
	cs = [0.0]*3
	total_iter = 0
	final_params = None
	final_volume = None
	final_volume_filtered = None
	# do the projection matching
	for N_step in xrange(lstp):
		
		terminate = 0
		Iter = 0
 		while Iter < max_iter and terminate == 0:

			Iter += 1
			total_iter += 1

			mpi_barrier(mpi_comm)
			if myid == main_node:
				log.add("ITERATION #%3d,  inner iteration #%3d\nDelta = %4.1f, an = %5.2f, xrange = %5.2f, yrange = %5.2f, step = %5.2f, delta psi = %5.2f, start psi = %5.2f\n"%(total_iter, Iter, delta[N_step], an[N_step], xrng[N_step],yrng[N_step],step[N_step],deltapsi[N_step],startpsi[N_step]))
				start_time = time()

			#=========================================================================
			# build references
			volft, kb = prep_vol(vol)
			refrings = prepare_refrings(volft, kb, nx, delta[N_step], ref_a, sym, numr, MPI=mpi_comm)
			del volft, kb
			#=========================================================================
			
 			if myid == main_node:
 				log.add("Time to prepare rings: %f\n" % (time()-start_time))
 				start_time = time()
			
			#=========================================================================
			if total_iter == 1:
				# adjust params to references, calculate psi+shifts, calculate previousmax
				for im in xrange(nima):
					stable = data[im].get_attr_default("stable", 0)
					if stable == 0:
						data[im].set_attr("previousmax", -1.0e23)
						data[im].set_attr("stable", 1)
					else:
						peak, pixer[im] = proj_ali_incore_local(data[im],refrings,numr,xrng[N_step],yrng[N_step],step[N_step],1.0)
						data[im].set_attr("previousmax", peak)
				if myid == main_node:
					log.add("Time to calculate first psi+shifts+previousmax: %f\n" % (time()-start_time))
					start_time = time()
			#=========================================================================

			mpi_barrier(mpi_comm)
			if myid == main_node:
				start_time = time()
			#=========================================================================
			# alignment
			number_of_checked_refs = 0
			for im in xrange(nima):
				peak, pixer[im], checked_refs = shc(data[im],refrings,numr,xrng[N_step],yrng[N_step],step[N_step],an[N_step])
				number_of_checked_refs += checked_refs
			#=========================================================================
			mpi_barrier(mpi_comm)
			if myid == main_node:
				log.add("Time of alignment = %f\n"%(time()-start_time))
				start_time = time()

			#=========================================================================
			#output pixel errors, check stop criterion
			from mpi import mpi_gatherv
			recvbuf = mpi_gatherv(pixer, nima, MPI_FLOAT, recvcount, disps, MPI_FLOAT, main_node, mpi_comm)
			total_checked_refs = wrap_mpi_gatherv([number_of_checked_refs], main_node, mpi_comm)
			terminate = 0
			if myid == main_node:
				recvbuf = map(float, recvbuf)
				total_checked_refs = sum(total_checked_refs)
				from statistics import hist_list
				lhist = 20
				region, histo = hist_list(recvbuf, lhist)
				log.add("=========================")
				for lhx in xrange(lhist):
					msg = " %10.3f     %7d"%(region[lhx], histo[lhx])
					log.add(msg)
				if (max(recvbuf) < 0.5) and (sum(recvbuf) < 1.0):
					terminate = 1
				if total_checked_refs >= proc_of_checked_refs * (len(refrings) * total_nima):
					terminate = 1
			terminate = wrap_mpi_bcast(terminate, main_node, mpi_comm)
			del recvbuf
			#=========================================================================

			#=========================================================================
			# centering
			if center == -1 and sym[0] == 'c':
				from utilities      import estimate_3D_center_MPI, rotate_3D_shift
				cs[0], cs[1], cs[2], dummy, dummy = estimate_3D_center_MPI(data, total_nima, myid, number_of_proc, main_node, mpi_comm=mpi_comm)
# 				if myid == main_node:
# 					msg = " Average center x = %10.3f        Center y = %10.3f        Center z = %10.3f\n"%(cs[0], cs[1], cs[2])
# 					print_msg(msg)
				if int(sym[1]) > 1:
					cs[0] = cs[1] = 0.0
# 					if myid == main_node:
# 						print_msg("For symmetry group cn (n>1), we only center the volume in z-direction\n")
				cs = mpi_bcast(cs, 3, MPI_FLOAT, main_node, mpi_comm)
				cs = [-float(cs[0]), -float(cs[1]), -float(cs[2])]
				rotate_3D_shift(data, cs)
			#=========================================================================

			#=========================================================================
			# volume reconstruction
			vol_previous = vol
			mpi_barrier(mpi_comm)
			if myid == main_node:
				start_time = time()
			if CTF: vol = recons3d_4nn_ctf_MPI(myid, data, snr, symmetry=sym, npad=npad, mpi_comm=mpi_comm)
			else:   vol = recons3d_4nn_MPI    (myid, data,      symmetry=sym, npad=npad, mpi_comm=mpi_comm)
	
			# log
			if myid == main_node:
				log.add("3D reconstruction time = %f\n"%(time()-start_time))
				start_time = time()
			if fourvar:
			#  Compute Fourier variance
				for im in xrange(nima):
					original_data[im].set_attr( 'xform.projection', data[im].get_attr('xform.projection') )
				varf = varf3d_MPI(original_data, ssnr_text_file = None, mask2D = None, reference_structure = vol, ou = last_ring, rw = 1.0, npad = 1, CTF = CTF, sign = 1, sym =sym, myid = myid, mpi_comm=mpi_comm)
				if myid == main_node:
					log.add("Time to calculate 3D Fourier variance= %f\n"%(time()-start_time))
					start_time = time()
					varf = 1.0/varf
			else:
				varf = None
			# user functions + save volume
			if myid == main_node:
				ref_data[2] = vol
				ref_data[3] = None #fscc
				ref_data[4] = varf
				#  call user-supplied function to prepare reference image, i.e., center and filter it
				vol, cs = user_func(ref_data)
			del varf
			# broadcast volume
			bcast_EMData_to_all(vol, myid, main_node, comm=mpi_comm)
			#=========================================================================
			# log
			if myid == main_node:
				log.add("user function + Bcast time = %d\n"%(time()-start_time))
				start_time = time()

	#=========================================================================
	# gather parameters
	params = []
	previousmax = []
	for im in data:
		t = get_params_proj(im)
		p = im.get_attr("previousmax")
		params.append( [t[0], t[1], t[2], t[3], t[4]] )
		previousmax.append(p)
	params = wrap_mpi_gatherv(params, main_node, mpi_comm)
	previousmax = wrap_mpi_gatherv(previousmax, main_node, mpi_comm)

	if myid == main_node: 
		log.add("Finish")
		return vol, params, previousmax
	else:
		return None, None, None  # results for the other processes


def shc(data, refrings, numr, xrng, yrng, step, an, finfo=None):
	from utilities    import compose_transform2
	from math         import cos, sin, pi
	from EMAN2 import Vec2f

	ID = data.get_attr("ID")

	number_of_checked_refs = 0

	mode = "F"
	nx   = data.get_xsize()
	ny   = data.get_ysize()
	#  center is in SPIDER convention
	cnx  = nx//2 + 1
	cny  = ny//2 + 1

	ant = cos(an*pi/180.0)
	#phi, theta, psi, sxo, syo = get_params_proj(data)
	t1 = data.get_attr("xform.projection")
	dp = t1.get_params("spider")
	if finfo:
		finfo.write("Image id: %6d\n"%(ID))
		#finfo.write("Old parameters: %9.4f %9.4f %9.4f %9.4f %9.4f\n"%(phi, theta, psi, sxo, syo))
		finfo.write("Old parameters: %9.4f %9.4f %9.4f %9.4f %9.4f\n"%(dp["phi"], dp["theta"], dp["psi"], -dp["tx"], -dp["ty"]))
		finfo.flush()

	[ang, sxs, sys, mirror, iref, peak, checked_refs] = Util.shc(data, refrings, xrng, yrng, step, ant, mode, numr, cnx+dp["tx"], cny+dp["ty"])
	iref=int(iref)
	number_of_checked_refs += int(checked_refs)
	#[ang,sxs,sys,mirror,peak,numref] = apmq_local(projdata[imn], ref_proj_rings, xrng, yrng, step, ant, mode, numr, cnx-sxo, cny-syo)
	#ang = (ang+360.0)%360.0
	if iref > -1:
		# The ormqip returns parameters such that the transformation is applied first, the mirror operation second.
		# What that means is that one has to change the the Eulerian angles so they point into mirrored direction: phi+180, 180-theta, 180-psi
		angb, sxb, syb, ct = compose_transform2(0.0, sxs, sys, 1, -ang, 0.0, 0.0, 1)
		if  mirror:
			phi   = (refrings[iref].get_attr("phi")+540.0)%360.0
			theta = 180.0-refrings[iref].get_attr("theta")
			psi   = (540.0-refrings[iref].get_attr("psi")+angb)%360.0
			s2x   = sxb - dp["tx"]
			s2y   = syb - dp["ty"]
		else:
			phi   = refrings[iref].get_attr("phi")
			theta = refrings[iref].get_attr("theta")
			psi   = (refrings[iref].get_attr("psi")+angb+360.0)%360.0
			s2x   = sxb - dp["tx"]
			s2y   = syb - dp["ty"]

		#set_params_proj(data, [phi, theta, psi, s2x, s2y])
		t2 = Transform({"type":"spider","phi":phi,"theta":theta,"psi":psi})
		t2.set_trans(Vec2f(-s2x, -s2y))
		data.set_attr("xform.projection", t2)
		from pixel_error import max_3D_pixel_error
		pixel_error = max_3D_pixel_error(t1, t2, numr[-3])
		if finfo:
			finfo.write( "New parameters: %9.4f %9.4f %9.4f %9.4f %9.4f %10.5f  %11.3e\n\n" %(phi, theta, psi, s2x, s2y, peak, pixel_error))
			finfo.flush()
		return peak, pixel_error, number_of_checked_refs
	else:
		return -1.0e23, 0.0, number_of_checked_refs


def ali3d_shc2(stack, ref_vol, outdir, maskfile = None, ir = 1, ou = -1, rs = 1, 
            xr = "4 2 2 1", yr = "-1", ts = "1 1 0.5 0.25", delta = "10 6 4 4", an = "-1", apsi = "-1", deltapsi = "-1", startpsi = "-1",
	    center = -1, maxit = 5, CTF = False, snr = 1.0,  ref_a = "S", sym = "c1",
	    user_func_name = "ref_ali3d", fourvar = True, npad = 4, debug = False, MPI = False, termprec = 0.0):
	"""
		Name
			ali3d - Perform 3-D projection matching given initial reference volume and image series
		Input
			stack: set of 2-D images in a stack file, images have to be squares
			ref_vol: initial reference volume
			outdir: directory name into which the results will be written
			maskfile: filename of the file containing 3D mask.
			ir: inner radius for rotational correlation > 0 
			ou: outer radius for rotational correlation <int(nx/2)-1 
			rs: steps between rings in rotational correlation >0
			xr: range for translation search in x direction in each iteration, search is +/xr
			yr: range for translation search in y direction in each iteration, search is +/yr
			ts: step size of the translation search in both directions, search is -xr, -xr+ts, 0, xr-ts, xr, can be fractional.
			delta: angular step for the reference projections in respective iterations
			an: angular neighborhood for local searches
			center: average center method
			max_iter: maximum iterations at each angle step
			CTF: if the flag is present, program will use the CTF information stored in file headers
			snr: signal noise ratio used in the 3D reconstruction
			ref_a: method for creating quasi-uniform distribution of the projection directions of reference projections: "S" - spiral
			sym: symmetry of the refined structure
			function: name of the user-supplied-function
			MPI: if presetm use MPI version
		Output
			output_directory: directory name into which the output files will be written.
			header: the alignment parameters are stored in the headers of input files as Transform Object xform.proj
	"""
	if MPI:
		ali3d_shcMPI2(stack, ref_vol, outdir, maskfile, ir, ou, rs, xr, yr, ts,
	        	delta, an, apsi, deltapsi, startpsi, center, maxit, CTF, snr, ref_a, sym, user_func_name,
				fourvar, npad, debug, termprec)

def ali3d_shcMPI2(stack, ref_vol, outdir, maskfile = None, ir = 1, ou = -1, rs = 1, 
            xr = "4 2 2 1", yr = "-1", ts = "1 1 0.5 0.25", delta = "10 6 4 4", an = "-1", apsi = "-1", deltapsi = "-1", startpsi = "-1",
	    center = -1, maxit = 5, CTF = False, snr = 1.0,  ref_a = "S", sym = "c1",  user_func_name = "ref_ali3d",
	    fourvar = True, npad = 4, debug = False, termprec = 0.0 ):

	from alignment       import Numrinit, prepare_refrings, proj_ali_incore, proj_ali_incore_local, proj_ali_incore_local_psi
	from utilities       import model_circle, get_image, drop_image, get_input_from_string
	from utilities       import bcast_list_to_all, bcast_number_to_all, reduce_EMData_to_root, bcast_EMData_to_all
	from utilities       import send_attr_dict, angle_between_projections_directions
	from utilities       import get_params_proj, file_type, set_params_proj
	from fundamentals    import rot_avg_image
	import os
	import types
	from utilities       import print_begin_msg, print_end_msg, print_msg
	from mpi             import mpi_bcast, mpi_comm_size, mpi_comm_rank, MPI_FLOAT, MPI_COMM_WORLD, mpi_barrier
	from mpi             import mpi_reduce, MPI_INT, MPI_SUM
	from filter          import filt_ctf
	from projection      import prep_vol, prgs
	from statistics      import hist_list, varf3d_MPI
	from applications    import MPI_start_end
	from development     import shc
	from math            import sqrt
	from random import random

	number_of_proc = mpi_comm_size(MPI_COMM_WORLD)
	myid           = mpi_comm_rank(MPI_COMM_WORLD)
	main_node = 0
	
	if os.path.exists(outdir):  ERROR('Output directory exists, please change the name and restart the program', "ali3d_MPI", 1, myid)
	mpi_barrier(MPI_COMM_WORLD)

	if myid == main_node:
		os.mkdir(outdir)
		import global_def
		global_def.LOGFILE =  os.path.join(outdir, global_def.LOGFILE)
		print_begin_msg("ali3d_shcMPI")
	mpi_barrier(MPI_COMM_WORLD)

	if debug:
		from time import sleep
		while not os.path.exists(outdir):
			print  "Node ",myid,"  waiting..."
			sleep(5)

		info_file = os.path.join(outdir, "progress%04d"%myid)
		finfo = open(info_file, 'w')
	else:
		finfo = None

	xrng        = get_input_from_string(xr)
	if  yr == "-1":  yrng = xrng
	else          :  yrng = get_input_from_string(yr)
	step        = get_input_from_string(ts)
	delta       = get_input_from_string(delta)
	lstp = min(len(xrng), len(yrng), len(step), len(delta))
	if an == "-1":
		an = [-1] * lstp
	else:
		an = get_input_from_string(an)

	if apsi == "-1":
		apsi = [-1] * lstp
	else:
		apsi = get_input_from_string(apsi)

	if deltapsi == "-1":
		deltapsi = [-1] * lstp
	else:
		deltapsi = get_input_from_string(deltapsi)

	if startpsi == "-1":
		startpsi = [-1] * lstp
	else:
		startpsi = get_input_from_string(startpsi)

	first_ring  = int(ir)
	rstep       = int(rs)
	last_ring   = int(ou)
	max_iter    = int(maxit)
	center      = int(center)

	vol     = EMData()
	vol.read_image(ref_vol)
	nx      = vol.get_xsize()
	if last_ring < 0:	last_ring = int(nx/2) - 2

	if myid == main_node:
		import user_functions
		user_func = user_functions.factory[user_func_name]

		print_msg("Input stack                 : %s\n"%(stack))
		print_msg("Reference volume            : %s\n"%(ref_vol))	
		print_msg("Output directory            : %s\n"%(outdir))
		print_msg("Maskfile                    : %s\n"%(maskfile))
		print_msg("Inner radius                : %i\n"%(first_ring))
		print_msg("Outer radius                : %i\n"%(last_ring))
		print_msg("Ring step                   : %i\n"%(rstep))
		print_msg("X search range              : %s\n"%(xrng))
		print_msg("Y search range              : %s\n"%(yrng))
		print_msg("Translational step          : %s\n"%(step))
		print_msg("Angular step                : %s\n"%(delta))
		print_msg("Angular search range (phi and theta)       : %s\n"%(an))
		print_msg("Angular search range (psi)                 : %s\n"%(apsi))
		print_msg("Delta psi                   : %s\n"%(deltapsi))
		print_msg("Start psi                   : %s\n"%(startpsi))
		print_msg("Maximum iteration           : %i\n"%(max_iter))
		print_msg("Percentage of change for termination: %f\n"%(termprec))
		print_msg("Center type                 : %i\n"%(center))
		print_msg("CTF correction              : %s\n"%(CTF))
		print_msg("Signal-to-Noise Ratio       : %f\n"%(snr))
		print_msg("Reference projection method : %s\n"%(ref_a))
		print_msg("Symmetry group              : %s\n\n"%(sym))

	if maskfile:
		if type(maskfile) is types.StringType: mask3D = get_image(maskfile)
		else:                                  mask3D = maskfile
	else: mask3D = model_circle(last_ring, nx, nx, nx)

	numr	= Numrinit(first_ring, last_ring, rstep, "F")
	mask2D  = model_circle(last_ring,nx,nx) - model_circle(first_ring,nx,nx)

	fscmask = mask3D  #model_circle(last_ring,nx,nx,nx)  For a fancy mask circle would work better  PAP 7/21/11
	if CTF:
		from reconstruction import rec3D_MPI
		from filter         import filt_ctf
	else:	 from reconstruction import rec3D_MPI_noCTF

	if myid == main_node:
		if file_type(stack) == "bdb":
			from EMAN2db import db_open_dict
			dummy = db_open_dict(stack, True)

		# horatio active_refactoring Jy51i1EwmLD4tWZ9_00000_1
		# active = EMUtil.get_all_attributes(stack, 'active')
		# list_of_particles = []
		# for im in xrange(len(active)):
		# 	if active[im]:  list_of_particles.append(im)
		# del active
		# nima = len(list_of_particles)
		
		nima = EMUtil.get_image_count(stack)
		list_of_particles = range(nima)



	else:
		nima = 0
	total_nima = bcast_number_to_all(nima, source_node = main_node)

	if myid != main_node:
		list_of_particles = [-1]*total_nima
	list_of_particles = bcast_list_to_all(list_of_particles, source_node = main_node)

	image_start, image_end = MPI_start_end(total_nima, number_of_proc, myid)
	# create a list of images for each node
	list_of_particles = list_of_particles[image_start: image_end]
	nima = len(list_of_particles)
	if debug:
		finfo.write("image_start, image_end: %d %d\n" %(image_start, image_end))
		finfo.flush()

	data = EMData.read_images(stack, list_of_particles)
	if fourvar:  original_data = []
	for im in xrange(nima):
		data[im].set_attr('ID', list_of_particles[im])
		if fourvar: original_data.append(data[im].copy())
		if CTF:
			ctf_params = data[im].get_attr("ctf")
			st = Util.infomask(data[im], mask2D, False)
			data[im] -= st[0]
			data[im] = filt_ctf(data[im], ctf_params)
			data[im].set_attr('ctf_applied', 1)

	if debug:
		finfo.write( '%d loaded  \n' % nima )
		finfo.flush()
	if myid == main_node:
		# initialize data for the reference preparation function
		ref_data = [ mask3D, max(center,0), None, None, None, None ]
		# for method -1, switch off centering in user function

	from time import time	

	#  this is needed for gathering of pixel errors
	disps = []
	recvcount = []
	for im in xrange(number_of_proc):
		if im == main_node :  disps.append(0)
		else:                 disps.append(disps[im-1] + recvcount[im-1])
		ib, ie = MPI_start_end(total_nima, number_of_proc, im)
		recvcount.append(ie-ib)

# 	# save stable set
# 	stable_set = set()
# 	for i in xrange(nima):
# 		img = data[i]
# 		if img.get_attr("stable") == 1:
# 			stable_set.add(i)

	pixer = [0.0]*nima
	cs = [0.0]*3
	total_iter = 0
	final_params = None
	final_volume = None
	final_volume_filtered = None
	# do the projection matching
	for N_step in xrange(lstp):
		
		terminate = 0
		Iter = 0
 		while Iter < max_iter and terminate == 0:

			Iter += 1
			total_iter += 1

			if myid == main_node:
				start_time = time()
				print_msg("\nITERATION #%3d,  inner iteration #%3d\nDelta = %4.1f, an = %5.2f, xrange = %5.2f, yrange = %5.2f, step = %5.2f, delta psi = %5.2f, start psi = %5.2f\n"%(total_iter, Iter, delta[N_step], an[N_step], xrng[N_step],yrng[N_step],step[N_step],deltapsi[N_step],startpsi[N_step]))

			#=========================================================================
			# build references
			volft, kb = prep_vol(vol)
			refrings = prepare_refrings(volft, kb, nx, delta[N_step], ref_a, sym, numr, True)
			del volft, kb
			#=========================================================================
			
			if myid == main_node:
				print_msg("Time to prepare rings: %d\n" % (time()-start_time))
				start_time = time()
			
			#=========================================================================
			# adjust params to references, calculate psi+shifts, calculate previousmax
			if total_iter <= 1:
				for im in xrange(nima):
					data[im].set_attr("previousmax", -1.0e23)
# 			for im in xrange(nima):
# 				stable = data[im].get_attr_default("stable", 0)
# 				if stable == 0:
# 					data[im].set_attr("previousmax", -1.0e23)
# 					if total_iter > 1:
# 						data[im].set_attr("stable", 1)
# 				else:
# 					peak, pixer[im] = proj_ali_incore_local(data[im],refrings,numr,xrng[N_step],yrng[N_step],step[N_step],1.0,finfo)
# 					data[im].set_attr("previousmax", peak)
			#=========================================================================

			if myid == main_node:
				print_msg("Time to calculate psi+shifts+previousmax: %d\n" % (time()-start_time))
				start_time = time()

			#=========================================================================
			# alignment
			for im in xrange(nima):
				#=============================================================================
				# find reference rings to check for each projection
				if data[im].get_attr_default("stable", 0) == 1 and an[N_step] > 0.0:
					phi, theta, psi, sx, sy = get_params_proj(data[im])
					refrings_to_check = []
					for refr in refrings:
						ref_phi   = refr.get_attr("phi")
						ref_theta = refr.get_attr("theta")
						if angle_between_projections_directions([phi, theta], [ref_phi, ref_theta]) <= an[N_step]:
							refrings_to_check.append(refr)
				else:
					refrings_to_check = refrings
				#=============================================================================
				#print im, get_params_proj(data[im])
				peak, pixer[im], checked_refs = shc(data[im],refrings_to_check,numr,xrng[N_step],yrng[N_step],step[N_step],an[N_step],finfo)
				#print im, get_params_proj(data[im])
				"""
				if deltapsi[N_step] > 0.0:
					from alignment import proj_ali_incore_delta
					peak, pixer[im] = proj_ali_incore_delta(data[im],refrings,numr,xrng[N_step],yrng[N_step],step[N_step],startpsi[N_step],deltapsi[N_step],finfo)						
				elif an[N_step] == -1:
					peak, pixer[im] = proj_ali_incore(data[im],refrings,numr,xrng[N_step],yrng[N_step],step[N_step],finfo)
				else:
					if apsi[N_step] == -1:
						# it requires attribute previousmax to be set
						peak, pixer[im] = shc(data[im],refrings,numr,xrng[N_step],yrng[N_step],step[N_step],an[N_step],finfo)
					else:
						peak, pixer[im] = proj_ali_incore_local_psi(data[im],refrings,numr,xrng[N_step],yrng[N_step],step[N_step],an[N_step],apsi[N_step],finfo)
				"""
			#=========================================================================

			if myid == main_node:
				print_msg("Time of alignment = %d\n"%(time()-start_time))
				start_time = time()

			#=========================================================================
			#output pixel errors, check stop criterion
			from mpi import mpi_gatherv
			recvbuf = mpi_gatherv(pixer, nima, MPI_FLOAT, recvcount, disps, MPI_FLOAT, main_node, MPI_COMM_WORLD)
			mpi_barrier(MPI_COMM_WORLD)
			terminate = 0
			if myid == main_node:
				recvbuf = map(float, recvbuf)
				from statistics import hist_list
				lhist = 20
				region, histo = hist_list(recvbuf, lhist)
				if region[0] < 0.0:  region[0] = 0.0
				msg = "      Histogram of pixel errors\n      ERROR       number of particles\n"
				print_msg(msg)
				for lhx in xrange(lhist):
					msg = " %10.3f     %7d\n"%(region[lhx], histo[lhx])
					print_msg(msg)
				# Terminate if 95% within 1 pixel error
				#im = 0
				#for lhx in xrange(lhist):
				#	if region[lhx] > 1.0: break
				#	im += histo[lhx]
				#precn = 100*float(total_nima-im)/float(total_nima)
				#msg = " Number of particles that changed orientations %7d, percentage of total: %5.1f\n"%(total_nima-im, precn)
				#print_msg(msg)
				#if precn <= termprec:  terminate = 1
				if max(recvbuf) < 0.001: terminate = 1
				del region, histo
			del recvbuf
			terminate = mpi_bcast(terminate, 1, MPI_INT, 0, MPI_COMM_WORLD)
			terminate = int(terminate[0])
			#=========================================================================

			#=========================================================================
			# centering
			if center == -1 and sym[0] == 'c':
				from utilities      import estimate_3D_center_MPI, rotate_3D_shift
				cs[0], cs[1], cs[2], dummy, dummy = estimate_3D_center_MPI(data, total_nima, myid, number_of_proc, main_node)
				if myid == main_node:
					msg = " Average center x = %10.3f        Center y = %10.3f        Center z = %10.3f\n"%(cs[0], cs[1], cs[2])
					print_msg(msg)
				if int(sym[1]) > 1:
					cs[0] = cs[1] = 0.0
					if myid == main_node:
						print_msg("For symmetry group cn (n>1), we only center the volume in z-direction\n")
				cs = mpi_bcast(cs, 3, MPI_FLOAT, main_node, MPI_COMM_WORLD)
				cs = [-float(cs[0]), -float(cs[1]), -float(cs[2])]
				rotate_3D_shift(data, cs)
			#=========================================================================

			#=========================================================================
			# write out headers, under MPI writing has to be done sequentially
			mpi_barrier(MPI_COMM_WORLD)
			par_str = ['xform.projection', 'previousmax', 'ID']
			if myid == main_node:
				if(file_type(stack) == "bdb"):
					from utilities import recv_attr_dict_bdb
					recv_attr_dict_bdb(main_node, stack, data, par_str, image_start, image_end, number_of_proc)
				else:
					from utilities import recv_attr_dict
					recv_attr_dict(main_node, stack, data, par_str, image_start, image_end, number_of_proc)
				# save parameters to file
				paro = [None]*total_nima
				projs_headers = EMData.read_images(stack, range(total_nima), True)
				for im in xrange(total_nima):
					a1,a2,a3,a4,a5 = get_params_proj(projs_headers[im])
					previousmax = projs_headers[im].get_attr("previousmax")
					paro[im] = [a1,a2,a3,a4,a5,previousmax]
				from utilities import write_text_row
				write_text_row(paro,os.path.join(outdir, "params%04d.txt"%(total_iter)))
				final_params = paro
				del projs_headers
				del paro
				# ------- end of saving parameters to file
				print_msg("Time to write header information= %d\n"%(time()-start_time))
				start_time = time()
			else:
				send_attr_dict(main_node, data, par_str, image_start, image_end)
			#=========================================================================

			#=========================================================================
			# volume reconstruction
			vol_previous = vol
			if CTF: vol, fscc = rec3D_MPI(data, snr, sym, fscmask, os.path.join(outdir, "resolution%04d"%(total_iter)), myid, main_node, npad = npad)
			else:   vol, fscc = rec3D_MPI_noCTF(data, sym, fscmask, os.path.join(outdir, "resolution%04d"%(total_iter)), myid, main_node, npad = npad)
			# log
			if myid == main_node:
				print_msg("3D reconstruction time = %d\n"%(time()-start_time))
				start_time = time()
			if fourvar:
			#  Compute Fourier variance
				for im in xrange(nima):
					original_data[im].set_attr( 'xform.projection', data[im].get_attr('xform.projection') )
				varf = varf3d_MPI(original_data, ssnr_text_file = os.path.join(outdir, "ssnr%04d"%(total_iter)), mask2D = None, reference_structure = vol, ou = last_ring, rw = 1.0, npad = 1, CTF = CTF, sign = 1, sym =sym, myid = myid)
				if myid == main_node:
					print_msg("Time to calculate 3D Fourier variance= %d\n"%(time()-start_time))
					start_time = time()
					varf = 1.0/varf
			else:
				varf = None
			# user functions + save volume
			if myid == main_node:
				drop_image(vol, os.path.join(outdir, "vol%04d.hdf"%(total_iter)))
				final_volume = vol
				ref_data[2] = vol
				ref_data[3] = fscc
				ref_data[4] = varf
				#  call user-supplied function to prepare reference image, i.e., center and filter it
				vol, cs = user_func(ref_data)
				drop_image(vol, os.path.join(outdir, "volf%04d.hdf"%(total_iter)))
				final_volume_filtered = vol
				print_msg("L2 between this and previous volume: " + str(sqrt(vol.cmp("SqEuclidean",vol_previous,{"mask":mask3D,"zeromask":0,"normto":0}))) + "\n")
				print_msg("Dot product of the volume: " + str(vol.cmp("dot", vol, {"negative":0, "mask":mask3D})) + "\n")
			del varf
			# broadcast volume
			bcast_EMData_to_all(vol, myid, main_node)
			#=========================================================================

	if myid == main_node: 
		write_text_row(final_params, os.path.join(outdir, "params.txt"))
		drop_image(final_volume_filtered, os.path.join(outdir, "volume_filt.hdf"))
		drop_image(final_volume         , os.path.join(outdir, "volume.hdf"     ))
		print_end_msg("ali3d_shcMPI")



def assign_projs_to_refs_hungarian(projs, refrings, numr, xrng, yrng, step):
	from utilities  import compose_transform2
	from global_def import Util
	from pixel_error import max_3D_pixel_error
	nx = projs[0].get_xsize()
	ny = projs[0].get_ysize()
	mode = "F"
	cnx = nx//2+1
	cny = ny//2+1
	peaks  = [0]*len(projs)   # matrix of cross-correlations (first index - images, second - references)
	params = [0]*len(projs)   # matrix of transformations    (first index - images, second - references)
	pixel_errors = [0]*len(projs) 
	max_peak = 0
	for im in xrange(len(projs)):
		t1 = projs[im].get_attr("xform.projection")
		dp = t1.get_params("spider")
		results = Util.multiref_polar_ali_2d_peaklist(projs[im], refrings, xrng, yrng, step, mode, numr, cnx+dp["tx"], cny+dp["ty"])
		row_peaks  = [0]*len(refrings)
		row_params = [0]*len(refrings)
		row_pix_err = [0]*len(refrings)
		for iref in xrange(len(refrings)):
			peak   = results[iref*5 + 0]
			ang    = results[iref*5 + 1]
			sxs    = results[iref*5 + 2]
			sys    = results[iref*5 + 3]
			mirror = results[iref*5 + 4]
			# The ormqip returns parameters such that the transformation is applied first, the mirror operation second.
			#  What that means is that one has to change the the Eulerian angles so they point into mirrored direction: phi+180, 180-theta, 180-psi
			angb, sxb, syb, ct = compose_transform2(0.0, sxs, sys, 1, -ang, 0.0, 0.0, 1)
			if mirror:
				phi   = (refrings[iref].get_attr("phi")+540.0)%360.0
				theta = 180.0-refrings[iref].get_attr("theta")
				psi   = (540.0-refrings[iref].get_attr("psi")+angb)%360.0
				s2x   = sxb - dp["tx"]
				s2y   = syb - dp["ty"]
			else:
				phi   = refrings[iref].get_attr("phi")
				theta = refrings[iref].get_attr("theta")
				psi   = (refrings[iref].get_attr("psi")+angb+360.0)%360.0
				s2x   = sxb - dp["tx"]
				s2y   = syb - dp["ty"]
			#set_params_proj(data, [phi, theta, psi, s2x, s2y])
			t2 = Transform({"type":"spider","phi":phi,"theta":theta,"psi":psi})
			t2.set_trans(Vec2f(-s2x, -s2y))
			row_peaks[iref] = peak
			row_params[iref] = t2
			row_pix_err[iref] = max_3D_pixel_error(t1, t2, numr[-3])
			if peak > max_peak:
				max_peak = peak
		peaks[im] = row_peaks
		params[im] = row_params
		pixel_errors[im] = row_pix_err
	
	for row in peaks:
		for i in xrange(len(row)):
			row[i] = max_peak + 10 - row[i]
	
	return peaks, params, pixel_errors


def ali3d_saturn2_MPI(stack, ref_vol, outdir, maskfile = None, ir = 1, ou = -1, rs = 1, 
            xr = "4 2 2 1", yr = "-1", ts = "1 1 0.5 0.25", delta = "10 6 4 4", an = "-1", apsi = "-1", deltapsi = "-1", startpsi = "-1",
	    center = -1, maxit = 5, CTF = False, snr = 1.0,  ref_a = "S", sym = "c1",
	    user_func_name = "ref_ali3d", fourvar = True, npad = 4, debug = False, termprec = 0.0):

	from alignment      import proj_ali_incore_local
	from utilities      import model_circle, drop_image, get_image, get_input_from_string
	from utilities      import get_params_proj, file_type, bcast_EMData_to_all
	from utilities      import estimate_3D_center_MPI, rotate_3D_shift
	from filter         import filt_params, fit_tanh, filt_tanl, filt_ctf
	from statistics     import fsc_mask, hist_list
	import os
	import types
	from utilities      import print_begin_msg, print_end_msg, print_msg
	from alignment      import Numrinit, prepare_refrings
	from projection     import prep_vol
	from math           import sqrt
	from mpi            import mpi_bcast, mpi_comm_size, mpi_comm_rank, MPI_FLOAT, mpi_barrier, MPI_INT, MPI_COMM_WORLD
	from utilities      import wrap_mpi_bcast, wrap_mpi_gatherv
	from applications   import MPI_start_end
	from random         import shuffle

	import user_functions
	user_func = user_functions.factory[user_func_name]

	mpi_comm = MPI_COMM_WORLD
	main_node = 0

	number_of_proc = mpi_comm_size(mpi_comm)
	myid           = mpi_comm_rank(mpi_comm)

	if myid == 0:
		if os.path.exists(outdir):  ERROR('Output directory exists, please change the name and restart the program', "ali3d", 1)
		os.mkdir(outdir)
		import global_def
		global_def.LOGFILE =  os.path.join(outdir, global_def.LOGFILE)
		print_begin_msg("ali3d_saturn")

	xrng        = get_input_from_string(xr)
	if  yr == "-1":  yrng = xrng
	else          :  yrng = get_input_from_string(yr)
	step        = get_input_from_string(ts)
	delta       = get_input_from_string(delta)
	lstp = min(len(xrng), len(yrng), len(step), len(delta))
	if an == "-1":
		an = [-1] * lstp
	else:
		an = get_input_from_string(an)
	first_ring  = int(ir)
	rstep       = int(rs)
	last_ring   = int(ou)
	max_iter    = int(maxit)
	center      = int(center)

	vol     = EMData()
	vol.read_image(ref_vol)
	nx      = vol.get_xsize()
	if last_ring == -1:
		last_ring = nx/2 - 2

	if myid == main_node:
		print_msg("Input stack                 : %s\n"%(stack))
		print_msg("Reference volume            : %s\n"%(ref_vol))	
		print_msg("Output directory            : %s\n"%(outdir))
		print_msg("Maskfile                    : %s\n"%(maskfile))
		print_msg("Inner radius                : %i\n"%(first_ring))
		print_msg("Outer radius                : %i\n"%(last_ring))
		print_msg("Ring step                   : %i\n"%(rstep))
		print_msg("X search range              : %s\n"%(xrng))
		print_msg("Y search range              : %s\n"%(yrng))
		print_msg("Translational step          : %s\n"%(step))
		print_msg("Angular step                : %s\n"%(delta))
		print_msg("Angular search range        : %s\n"%(an))
		print_msg("Maximum iteration           : %i\n"%(max_iter))
		print_msg("Center type                 : %i\n"%(center))
		print_msg("CTF correction              : %s\n"%(CTF))
		print_msg("Signal-to-Noise Ratio       : %f\n"%(snr))
		print_msg("Reference projection method : %s\n"%(ref_a))
		print_msg("Symmetry group              : %s\n\n"%(sym))

	if maskfile :
		if type(maskfile) is types.StringType: mask3D = get_image(maskfile)
		else                                  : mask3D = maskfile
	else          :   mask3D = model_circle(last_ring, nx, nx, nx)
	mask2D = model_circle(last_ring, nx, nx) - model_circle(first_ring, nx, nx)
	numr   = Numrinit(first_ring, last_ring, rstep, "F")

	if CTF:
		from reconstruction import recons3d_4nn_ctf_MPI
		from filter         import filt_ctf
	else: 
		from reconstruction import recons3d_4nn_MPI

	#if debug:  outf = file(os.path.join(outdir, "progress"), "w")
	#else:      outf = None

	# ========================================================
#	active = EMUtil.get_all_attributes(stack, 'active')
#	list_of_particles = []
#	for im in xrange(len(active)):
#		if(active[im]):  list_of_particles.append(im)
#	del active
#	data = EMData.read_images(stack, list_of_particles)
#	for im in xrange(len(data)):
#		data[im].set_attr('ID', list_of_particles[im])
#		if CTF:
#			ctf_params = data[im].get_attr("ctf")
#			st = Util.infomask(data[im], mask2D, False)
#			data[im] -= st[0]
#			data[im] = filt_ctf(data[im], ctf_params)
#			data[im].set_attr('ctf_applied', 1)
	# =====================================================
	if myid == 0:
		if file_type(stack) == "bdb":
			from EMAN2db import db_open_dict
			dummy = db_open_dict(stack, True)

		# horatio active_refactoring Jy51i1EwmLD4tWZ9_00000_1
		# active = EMUtil.get_all_attributes(stack, 'active')
		# list_of_particles = []
		# for im in xrange(len(active)):
		# 	if active[im]:  list_of_particles.append(im)
		# del active
		# nima = len(list_of_particles)
		
		nima = EMUtil.get_image_count(stack)
		list_of_particles = range(nima)
	
	else:
		nima = 0
	total_nima = wrap_mpi_bcast(nima, main_node)

	if myid != 0:
		list_of_particles = [-1]*total_nima
	list_of_particles = wrap_mpi_bcast(list_of_particles, 0)

	image_start, image_end = MPI_start_end(total_nima, number_of_proc, myid)
	max_nima = (total_nima - 1) / number_of_proc + 1
	# create a list of images for each node
	list_of_particles = list_of_particles[image_start: image_end]
	nima = len(list_of_particles)
	if debug:
		finfo.write("image_start, image_end: %d %d\n" %(image_start, image_end))
		finfo.flush()

	data = EMData.read_images(stack, list_of_particles)
	
	nima = len(data)
	
	if fourvar:
		original_data = []
	for im in xrange(nima):
		data[im].set_attr('ID', list_of_particles[im])
		if fourvar:
			original_data.append(data[im].copy())
		if CTF:
			ctf_params = data[im].get_attr("ctf")
			st = Util.infomask(data[im], mask2D, False)
			data[im] -= st[0]
			data[im] = filt_ctf(data[im], ctf_params)
			data[im].set_attr('ctf_applied', 1)
	# =====================================================
	
	# initialize data for the reference preparation function
	ref_data = [ mask3D, max(center,0), None, None ]#  for center -1 switch of centering by user function
	
	cs = [0.0]*3
	# do the projection matching
	for N_step in xrange(lstp):
		for Iter in xrange(max_iter):
			if myid == 0:
				print_msg("\nITERATION #%3d\n"%(N_step*max_iter+Iter+1))

			volft, kb = prep_vol(vol)
			refrings = prepare_refrings( volft, kb, nx, delta[N_step], ref_a, sym, numr, MPI=True)
			del volft, kb
			rorder = range(nima)
			shuffle(rorder)
			pixel_error = [0.0]*nima
			
			if an[N_step] == -1:
				peaks, transforms, pixel_errors = assign_projs_to_refs_hungarian(data, refrings, numr, xrng[N_step],yrng[N_step],step[N_step])
				
				peaks  = wrap_mpi_gatherv(peaks , 0, mpi_comm)
				
				# =============== Hungarian algorithm
				if myid == 0:
					from statistics import Munkres
					m = Munkres()
					assignment_raw = m.compute(peaks)
					assignment = [0] * total_nima
					for row, col in assignment_raw:
						assignment[row] = col
				else:
					assignment = None
				assignment = wrap_mpi_bcast(assignment, 0, mpi_comm)
				
				for im in xrange(nima):
					iref = assignment[image_start + im]
					pixel_error[im] = pixel_errors[im][iref]
					data[im].set_attr("xform.projection", transforms[im][iref])
			else:
				ERROR('NOT IMPLEMENTED', "ali3d", 1)
			
			pixel_error = wrap_mpi_gatherv(pixel_error, 0, mpi_comm)
			if myid == 0:
				lhist = 20
				region, histo = hist_list(pixel_error, lhist)
				if region[0] < 0.0:  region[0] = 0.0
				msg = "      Histogram of pixel errors\n      ERROR       number of particles\n"
				print_msg(msg)
				for lhx in xrange(lhist):
					msg = " %10.3f     %7d\n"%(region[lhx], histo[lhx])
					print_msg(msg)
				# Terminate if 95% within 1 pixel error
				im = 0
				for lhx in xrange(lhist):
					if region[lhx] > 1.0: break
					im += histo[lhx]
				precn = 100*float(nima-im)/float(nima)
				msg = " Number of particles that changed orientations %7d, percentage of total: %5.1f\n"%(nima-im, precn)
				print_msg(msg)
			#if precn <= termprec:  terminate = 1   from shc
			
			if center == -1 and sym[0] == 'c':
				cs[0], cs[1], cs[2], dummy, dummy = estimate_3D_center_MPI(data, total_nima, myid, number_of_proc, main_node)
				if myid == 0:
					msg = "Average center x = %10.3f        Center y = %10.3f        Center z = %10.3f\n"%(cs[0], cs[1], cs[2])
					print_msg(msg)
				if int(sym[1]) > 1:
					cs[0] = cs[1] = 0.0
					if myid == 0:
						print_msg("For symmetry group cn (n>1), we only center the volume in z-direction\n")
				rotate_3D_shift(data, [-cs[0], -cs[1], -cs[2]])
			"""
			if CTF:    vol1 = recons3d_4nn_ctf(data, range(0, nima, 2), snr, 1, sym)
			else:	   vol1 = recons3d_4nn(data, range(0, nima, 2), sym)
			if CTF:    vol2 = recons3d_4nn_ctf(data, range(1, nima, 2), snr, 1, sym)
			else:	   vol2 = recons3d_4nn(data, range(1, nima, 2), sym)

			fscc = fsc_mask(vol1, vol2, mask3D, 1.0, os.path.join(outdir, "resolution%04d"%(N_step*max_iter+Iter+1)))
			del vol1
			del vol2
			"""
			# calculate new and improved 3D
			vol_previous = vol
			if CTF:
				vol = recons3d_4nn_ctf_MPI(myid, data, snr, 1, sym)
			else:
				vol = recons3d_4nn_MPI(myid, data, sym)
			# store the reference volume
			if myid == 0:
				drop_image(vol, os.path.join(outdir, "vol%04d.hdf"%(N_step*max_iter+Iter+1)))
				ref_data[2] = vol
				ref_data[3] = None#fscc
				current_result = vol
				#  call user-supplied function to prepare reference image, i.e., center and filter it
				vol, dummy = user_func(ref_data)
				print_msg("L2 between this and previous volume: " + str(sqrt(vol.cmp("SqEuclidean",vol_previous,{"mask":mask3D,"zeromask":0,"normto":0}))) + "\n")
				print_msg("Dot product of the volume: " + str(vol.cmp("dot", vol, {"negative":0, "mask":mask3D})) + "\n")

				drop_image(vol, os.path.join(outdir, "volf%04d.hdf"%(N_step*max_iter+Iter+1)))
				
			paro = [None]*nima
			for im in xrange(nima):
				a1,a2,a3,a4,a5 = get_params_proj(data[im])
				paro[im] = [a1,a2,a3,a4,a5]
			paro = wrap_mpi_gatherv(paro, 0)
			if myid == 0:
				from utilities import write_text_row
				write_text_row(paro,os.path.join(outdir, "params%04d.txt"%(N_step*max_iter+Iter+1)))
				del paro
			"""
			#  here we write header info
			from utilities import write_headers
			#from utilities import write_select_headers
			if CTF:
				for dat in data:  dat.set_attr('ctf_applied',0)
			for dat in data:  dat.del_attr('referencenumber')
			write_headers(stack, data, list_of_particles)
			#list_params= ['ID','xform.projection']
			#write_select_headers(stack, data, list_params)
			if CTF:
				for dat in data:  dat.set_attr('ctf_applied', 1)
			"""
			bcast_EMData_to_all(vol, myid, 0)

	if myid == 0:
		print_end_msg("ali3d")
	#return current_result


def ali3d_saturn2(stack, ref_vol, outdir, maskfile = None, ir = 1, ou = -1, rs = 1, 
            xr = "4 2 2 1", yr = "-1", ts = "1 1 0.5 0.25", delta = "10 6 4 4", an = "-1", apsi = "-1", deltapsi = "-1", startpsi = "-1",
	    center = -1, maxit = 5, CTF = False, snr = 1.0,  ref_a = "S", sym = "c1",
	    user_func_name = "ref_ali3d", fourvar = True, npad = 4, debug = False, MPI = False, termprec = 0.0):
	"""
		Name
			ali3d - Perform 3-D projection matching given initial reference volume and image series
		Input
			stack: set of 2-D images in a stack file, images have to be squares
			ref_vol: initial reference volume
			outdir: directory name into which the results will be written
			maskfile: filename of the file containing 3D mask.
			ir: inner radius for rotational correlation > 0 
			ou: outer radius for rotational correlation <int(nx/2)-1 
			rs: steps between rings in rotational correlation >0
			xr: range for translation search in x direction in each iteration, search is +/xr
			yr: range for translation search in y direction in each iteration, search is +/yr
			ts: step size of the translation search in both directions, search is -xr, -xr+ts, 0, xr-ts, xr, can be fractional.
			delta: angular step for the reference projections in respective iterations
			an: angular neighborhood for local searches
			center: average center method
			max_iter: maximum iterations at each angle step
			CTF: if the flag is present, program will use the CTF information stored in file headers
			snr: signal noise ratio used in the 3D reconstruction
			ref_a: method for creating quasi-uniform distribution of the projection directions of reference projections: "S" - spiral
			sym: symmetry of the refined structure
			function: name of the user-supplied-function
			MPI: if presetm use MPI version
		Output
			output_directory: directory name into which the output files will be written.
			header: the alignment parameters are stored in the headers of input files as Transform Object xform.proj
	"""
	if MPI:
		ali3d_saturn2_MPI(stack, ref_vol, outdir, maskfile, ir, ou, rs, xr, yr, ts, delta, an, apsi, deltapsi, startpsi, center, maxit, CTF
					, snr, ref_a, sym, user_func_name, fourvar, npad, debug, termprec)
		return
	else:
		ERROR('NOT IMPLEMENTED', "ali3d", 1)


# returns [ref_id, peak, pixel_error, Transform3D]*number_of_assigned_refs sorted by peak (desc)
def proj_ali_incore_peaklist(data, refrings, numr, xrng, yrng, step, finfo=None, number_of_assigned_refs=1):
	from utilities   import compose_transform2
	from EMAN2       import Vec2f
	from pixel_error import max_3D_pixel_error

	ID = data.get_attr("ID")
	if finfo:
		from utilities    import get_params_proj
		phi, theta, psi, s2x, s2y = get_params_proj(data)
		finfo.write("Image id: %6d\n"%(ID))
		finfo.write("Old parameters: %9.4f %9.4f %9.4f %9.4f %9.4f\n"%(phi, theta, psi, s2x, s2y))
		finfo.flush()

	mode = "F"
	#  center is in SPIDER convention
	nx   = data.get_xsize()
	ny   = data.get_ysize()
	cnx  = nx//2 + 1
	cny  = ny//2 + 1

	#phi, theta, psi, sxo, syo = get_params_proj(data)
	t1 = data.get_attr("xform.projection")
	dp = t1.get_params("spider")
	#[ang, sxs, sys, mirror, iref, peak] = Util.multiref_polar_ali_2d(data, refrings, xrng, yrng, step, mode, numr, cnx-sxo, cny-syo)
	# [ peak0, ang0, sx0, sy0, mirror0, peak1, ang1, ... ]
	results = Util.multiref_polar_ali_2d_peaklist(data, refrings, xrng, yrng, step, mode, numr, cnx+dp["tx"], cny+dp["ty"])
	
	# sorting by peak
	temp_list = []
	for i in xrange(len(refrings)):
		temp_list.append([results[i*5], i])
	temp_list.sort(reverse=True)
	temp_list = temp_list[:number_of_assigned_refs]
	
	# building output list
	output = []
	for t in temp_list:
		iref   = t[1]
		peak   = results[iref*5 + 0]
		ang    = results[iref*5 + 1]
		sxs    = results[iref*5 + 2]
		sys    = results[iref*5 + 3]
		mirror = results[iref*5 + 4]
		#[ang,sxs,sys,mirror,peak,numref] = apmq(projdata[imn], ref_proj_rings, xrng, yrng, step, mode, numr, cnx-sxo, cny-syo)
		#ang = (ang+360.0)%360.0
		# The ormqip returns parameters such that the transformation is applied first, the mirror operation second.
		#  What that means is that one has to change the the Eulerian angles so they point into mirrored direction: phi+180, 180-theta, 180-psi
		angb, sxb, syb, ct = compose_transform2(0.0, sxs, sys, 1, -ang, 0.0, 0.0, 1)
		if mirror:
			phi   = (refrings[iref].get_attr("phi")+540.0)%360.0
			theta = 180.0-refrings[iref].get_attr("theta")
			psi   = (540.0-refrings[iref].get_attr("psi")+angb)%360.0
			s2x   = sxb - dp["tx"]
			s2y   = syb - dp["ty"]
		else:
			phi   = refrings[iref].get_attr("phi")
			theta = refrings[iref].get_attr("theta")
			psi   = (refrings[iref].get_attr("psi")+angb+360.0)%360.0
			s2x   = sxb - dp["tx"]
			s2y   = syb - dp["ty"]
		#set_params_proj(data, [phi, theta, psi, s2x, s2y])
		t2 = Transform({"type":"spider","phi":phi,"theta":theta,"psi":psi})
		t2.set_trans(Vec2f(-s2x, -s2y))
		pixel_error = max_3D_pixel_error(t1, t2, numr[-3])
		output.append([ iref, peak, pixel_error, t2 ])

	return output


def ali3d_saturnMPI(stack, ref_vol, outdir, maskfile = None, ir = 1, ou = -1, rs = 1, 
            xr = "4 2 2 1", yr = "-1", ts = "1 1 0.5 0.25", delta = "10 6 4 4", an = "-1", apsi = "-1", deltapsi = "-1", startpsi = "-1",
	    center = -1, maxit = 5, CTF = False, snr = 1.0,  ref_a = "S", sym = "c1",
	    user_func_name = "ref_ali3d", fourvar = True, npad = 4, debug = False, termprec = 0.0):

	from alignment      import proj_ali_incore_local
	from utilities      import model_circle, drop_image, get_image, get_input_from_string
	from utilities      import get_params_proj, file_type, bcast_EMData_to_all
	from utilities      import estimate_3D_center_MPI, rotate_3D_shift
	from filter         import filt_params, fit_tanh, filt_tanl, filt_ctf
	from statistics     import fsc_mask, hist_list
	import os
	import types
	from utilities      import print_begin_msg, print_end_msg, print_msg
	from alignment      import Numrinit, prepare_refrings
	from projection     import prep_vol
	from math           import sqrt
	from mpi            import mpi_bcast, mpi_comm_size, mpi_comm_rank, MPI_FLOAT, mpi_barrier, MPI_INT, MPI_COMM_WORLD
	from utilities      import wrap_mpi_bcast, wrap_mpi_gatherv
	from applications   import MPI_start_end
	from random         import shuffle

	import user_functions
	user_func = user_functions.factory[user_func_name]

	mpi_comm = MPI_COMM_WORLD
	main_node = 0

	number_of_proc = mpi_comm_size(mpi_comm)
	myid           = mpi_comm_rank(mpi_comm)

	if myid == 0:
		if os.path.exists(outdir):  ERROR('Output directory exists, please change the name and restart the program', "ali3d", 1)
		os.mkdir(outdir)
		import global_def
		global_def.LOGFILE =  os.path.join(outdir, global_def.LOGFILE)
		print_begin_msg("ali3d_saturn")

	xrng        = get_input_from_string(xr)
	if  yr == "-1":  yrng = xrng
	else          :  yrng = get_input_from_string(yr)
	step        = get_input_from_string(ts)
	delta       = get_input_from_string(delta)
	lstp = min(len(xrng), len(yrng), len(step), len(delta))
	if an == "-1":
		an = [-1] * lstp
	else:
		an = get_input_from_string(an)
	first_ring  = int(ir)
	rstep       = int(rs)
	last_ring   = int(ou)
	max_iter    = int(maxit)
	center      = int(center)

	vol     = EMData()
	vol.read_image(ref_vol)
	nx      = vol.get_xsize()
	if last_ring == -1:
		last_ring = nx/2 - 2

	if myid == main_node:
		print_msg("Input stack                 : %s\n"%(stack))
		print_msg("Reference volume            : %s\n"%(ref_vol))	
		print_msg("Output directory            : %s\n"%(outdir))
		print_msg("Maskfile                    : %s\n"%(maskfile))
		print_msg("Inner radius                : %i\n"%(first_ring))
		print_msg("Outer radius                : %i\n"%(last_ring))
		print_msg("Ring step                   : %i\n"%(rstep))
		print_msg("X search range              : %s\n"%(xrng))
		print_msg("Y search range              : %s\n"%(yrng))
		print_msg("Translational step          : %s\n"%(step))
		print_msg("Angular step                : %s\n"%(delta))
		print_msg("Angular search range        : %s\n"%(an))
		print_msg("Maximum iteration           : %i\n"%(max_iter))
		print_msg("Center type                 : %i\n"%(center))
		print_msg("CTF correction              : %s\n"%(CTF))
		print_msg("Signal-to-Noise Ratio       : %f\n"%(snr))
		print_msg("Reference projection method : %s\n"%(ref_a))
		print_msg("Symmetry group              : %s\n\n"%(sym))

	if maskfile :
		if type(maskfile) is types.StringType: mask3D = get_image(maskfile)
		else                                  : mask3D = maskfile
	else          :   mask3D = model_circle(last_ring, nx, nx, nx)
	mask2D = model_circle(last_ring, nx, nx) - model_circle(first_ring, nx, nx)
	numr   = Numrinit(first_ring, last_ring, rstep, "F")

	if CTF:
		from reconstruction import recons3d_4nn_ctf_MPI
		from filter         import filt_ctf
	else: 
		from reconstruction import recons3d_4nn_MPI

	#if debug:  outf = file(os.path.join(outdir, "progress"), "w")
	#else:      outf = None

	# ========================================================
#	active = EMUtil.get_all_attributes(stack, 'active')
#	list_of_particles = []
#	for im in xrange(len(active)):
#		if(active[im]):  list_of_particles.append(im)
#	del active
#	data = EMData.read_images(stack, list_of_particles)
#	for im in xrange(len(data)):
#		data[im].set_attr('ID', list_of_particles[im])
#		if CTF:
#			ctf_params = data[im].get_attr("ctf")
#			st = Util.infomask(data[im], mask2D, False)
#			data[im] -= st[0]
#			data[im] = filt_ctf(data[im], ctf_params)
#			data[im].set_attr('ctf_applied', 1)
	# =====================================================
	if myid == 0:
		if file_type(stack) == "bdb":
			from EMAN2db import db_open_dict
			dummy = db_open_dict(stack, True)

		# horatio active_refactoring Jy51i1EwmLD4tWZ9_00000_1
		# active = EMUtil.get_all_attributes(stack, 'active')
		# list_of_particles = []
		# for im in xrange(len(active)):
		# 	if active[im]:  list_of_particles.append(im)
		# del active
		# nima = len(list_of_particles)
		nima = EMUtil.get_image_count(stack)
		list_of_particles = range(nima)
	
		
	else:
		nima = 0
	total_nima = wrap_mpi_bcast(nima, main_node)

	if myid != 0:
		list_of_particles = [-1]*total_nima
	list_of_particles = wrap_mpi_bcast(list_of_particles, 0)

	image_start, image_end = MPI_start_end(total_nima, number_of_proc, myid)
	max_nima = (total_nima - 1) / number_of_proc + 1
	# create a list of images for each node
	list_of_particles = list_of_particles[image_start: image_end]
	nima = len(list_of_particles)
	if debug:
		finfo.write("image_start, image_end: %d %d\n" %(image_start, image_end))
		finfo.flush()

	data = EMData.read_images(stack, list_of_particles)
	
	nima = len(data)
	
	if fourvar:
		original_data = []
	for im in xrange(nima):
		data[im].set_attr('ID', list_of_particles[im])
		if fourvar:
			original_data.append(data[im].copy())
		if CTF:
			ctf_params = data[im].get_attr("ctf")
			st = Util.infomask(data[im], mask2D, False)
			data[im] -= st[0]
			data[im] = filt_ctf(data[im], ctf_params)
			data[im].set_attr('ctf_applied', 1)
	# =====================================================
	
	# initialize data for the reference preparation function
	ref_data = [ mask3D, max(center,0), None, None ]#  for center -1 switch of centering by user function
	
	cs = [0.0]*3
	# do the projection matching
	for N_step in xrange(lstp):
		for Iter in xrange(max_iter):
			if myid == 0:
				print_msg("\nITERATION #%3d\n"%(N_step*max_iter+Iter+1))

			volft, kb = prep_vol(vol)
			refrings = prepare_refrings( volft, kb, nx, delta[N_step], ref_a, sym, numr, MPI=True)
			del volft, kb
			rorder = range(nima)
			shuffle(rorder)
			pixel_error = [0.0]*nima
			
			for im in xrange(nima):
				if an[N_step] == -1:	
					# ref params - list of lists sorted by peak in descending order ( [ref_id, peak, pixel_error, Transform3D] )
					ref_params = proj_ali_incore_peaklist(data[rorder[im]],refrings,numr,xrng[N_step],yrng[N_step],step[N_step], number_of_assigned_refs=number_of_proc)
					#print ref_params
					ref_index = []
					for i in ref_params:
						ref_index.append(i[0])
					
					all_ref_index = wrap_mpi_gatherv([ref_index], 0)
					
					assigned_ref = None
					if myid == 0:
						assigned_ref = [-1] * number_of_proc
						temp_set = set()
						for i in xrange(number_of_proc):
							iref = -1
							for iref in all_ref_index[i]:
								if iref not in temp_set:
									temp_set.add(iref)
									break
							assigned_ref[i] = iref
					
					assigned_ref = wrap_mpi_bcast(assigned_ref, 0)
					iref = assigned_ref[myid]
					for r in ref_params:
						if r[0] == iref:
							break
					ref_params = r
					peak = ref_params[1]
					pixel_error[im] = ref_params[2]
					data[rorder[im]].set_attr("xform.projection", ref_params[3])
					
					# deleting occupied references
					assigned_ref.sort(reverse=True)
					for i in assigned_ref:
						if i >= 0:
							del refrings[i]
				else:
					ERROR('NOT IMPLEMENTED', "ali3d", 1)
					# NOT implemented
					#peak, pixel_error[im] = proj_ali_incore_local(data[im],refrings,numr,xrng[N_step],yrng[N_step],step[N_step],an[N_step])
			
			for im in xrange(nima,max_nima):
				wrap_mpi_gatherv([[]], 0)
				wrap_mpi_bcast(None, 0)
			
			pixel_error = wrap_mpi_gatherv(pixel_error, 0)
			if myid == 0:
				lhist = 20
				region, histo = hist_list(pixel_error, lhist)
				if region[0] < 0.0:  region[0] = 0.0
				msg = "      Histogram of pixel errors\n      ERROR       number of particles\n"
				print_msg(msg)
				for lhx in xrange(lhist):
					msg = " %10.3f     %7d\n"%(region[lhx], histo[lhx])
					print_msg(msg)
				# Terminate if 95% within 1 pixel error
				im = 0
				for lhx in xrange(lhist):
					if region[lhx] > 1.0: break
					im += histo[lhx]
				precn = 100*float(nima-im)/float(nima)
				msg = " Number of particles that changed orientations %7d, percentage of total: %5.1f\n"%(nima-im, precn)
				print_msg(msg)
			#if precn <= termprec:  terminate = 1   from shc
			
			if center == -1 and sym[0] == 'c':
				cs[0], cs[1], cs[2], dummy, dummy = estimate_3D_center_MPI(data, total_nima, myid, number_of_proc, main_node)
				if myid == 0:
					msg = "Average center x = %10.3f        Center y = %10.3f        Center z = %10.3f\n"%(cs[0], cs[1], cs[2])
					print_msg(msg)
				if int(sym[1]) > 1:
					cs[0] = cs[1] = 0.0
					if myid == 0:
						print_msg("For symmetry group cn (n>1), we only center the volume in z-direction\n")
				rotate_3D_shift(data, [-cs[0], -cs[1], -cs[2]])
			"""
			if CTF:    vol1 = recons3d_4nn_ctf(data, range(0, nima, 2), snr, 1, sym)
			else:	   vol1 = recons3d_4nn(data, range(0, nima, 2), sym)
			if CTF:    vol2 = recons3d_4nn_ctf(data, range(1, nima, 2), snr, 1, sym)
			else:	   vol2 = recons3d_4nn(data, range(1, nima, 2), sym)

			fscc = fsc_mask(vol1, vol2, mask3D, 1.0, os.path.join(outdir, "resolution%04d"%(N_step*max_iter+Iter+1)))
			del vol1
			del vol2
			"""
			# calculate new and improved 3D
			vol_previous = vol
			if CTF:
				vol = recons3d_4nn_ctf_MPI(myid, data, snr, 1, sym)
			else:
				vol = recons3d_4nn_MPI(myid, data, sym)
			# store the reference volume
			if myid == 0:
				drop_image(vol, os.path.join(outdir, "vol%04d.hdf"%(N_step*max_iter+Iter+1)))
				ref_data[2] = vol
				ref_data[3] = None#fscc
				current_result = vol
				#  call user-supplied function to prepare reference image, i.e., center and filter it
				vol, dummy = user_func(ref_data)
				print_msg("L2 between this and previous volume: " + str(sqrt(vol.cmp("SqEuclidean",vol_previous,{"mask":mask3D,"zeromask":0,"normto":0}))) + "\n")
				print_msg("Dot product of the volume: " + str(vol.cmp("dot", vol, {"negative":0, "mask":mask3D})) + "\n")

				drop_image(vol, os.path.join(outdir, "volf%04d.hdf"%(N_step*max_iter+Iter+1)))
				
			paro = [None]*nima
			for im in xrange(nima):
				a1,a2,a3,a4,a5 = get_params_proj(data[im])
				paro[im] = [a1,a2,a3,a4,a5]
			paro = wrap_mpi_gatherv(paro, 0)
			if myid == 0:
				from utilities import write_text_row
				write_text_row(paro,os.path.join(outdir, "params%04d.txt"%(N_step*max_iter+Iter+1)))
				del paro
			"""
			#  here we write header info
			from utilities import write_headers
			#from utilities import write_select_headers
			if CTF:
				for dat in data:  dat.set_attr('ctf_applied',0)
			for dat in data:  dat.del_attr('referencenumber')
			write_headers(stack, data, list_of_particles)
			#list_params= ['ID','xform.projection']
			#write_select_headers(stack, data, list_params)
			if CTF:
				for dat in data:  dat.set_attr('ctf_applied', 1)
			"""
			bcast_EMData_to_all(vol, myid, 0)

	if myid == 0:
		print_end_msg("ali3d")
	#return current_result


'''
def ali3d_a(stack, ref_vol, outdir, maskfile = None, ir = 1, ou = -1, rs = 1, 
            xr = "4 2 2 1", yr = "-1", ts = "1 1 0.5 0.25", delta="10 6 4 4", an="-1", 
	    center = 1.0, maxit = 5, CTF = False, ref_a = "S", sym="c1", user_func_name="ref_ali3d"):

	from utilities      import model_circle, drop_image
	from utilities      import get_image, get_input_from_string
	from utilities      import get_arb_params, set_arb_params
	from filter         import filt_params, filt_btwl, filt_from_fsc, filt_table, fit_tanh, filt_tanl
	from alignment	    import proj_ali_incore, proj_ali_incore_local
	from statistics     import fsc_mask
	from fundamentals   import fft
	from reconstruction import recons3d_nn_SSNR
	import os
	import types
	from utilities      import print_begin_msg, print_end_msg, print_msg

	import user_functions
	user_func = user_functions.factory[user_func_name]

	if CTF : from reconstruction import recons3d_4nn_ctf
	else   : from reconstruction import recons3d_4nn

	if os.path.exists(outdir):  ERROR('Output directory exists, please change the name and restart the program', " ali3d_a", 1)
	os.mkdir(outdir)
	import global_def
	global_def.LOGFILE =  os.path.join(outdir, global_def.LOGFILE)
	print_begin_msg("ali3d_a")

	xrng        = get_input_from_string(xr)
	if  yr == "-1":  yrng = xrng
	else          :  yrng = get_input_from_string(yr)
	step        = get_input_from_string(ts)
	delta       = get_input_from_string(delta)
	lstp = min( len(xrng), len(yrng), len(step), len(delta) )
	if (an == "-1"):
		an = []
		for i in xrange(lstp):   an.append(-1)
	else:
		an = get_input_from_string(an)
	first_ring  = int(ir)
	rstep       = int(rs)
	last_ring   = int(ou)
	max_iter    = int(maxit)
	center      = int(center)

	print_msg("Input stack                 : %s\n"%(stack))
	print_msg("Reference volume            : %s\n"%(ref_vol))	
	print_msg("Output directory            : %s\n"%(outdir))
	print_msg("Maskfile                    : %s\n"%(maskfile))
	print_msg("Inner radius                : %i\n"%(first_ring))

	vol     = EMData()
	vol.read_image(ref_vol)
	nx      = vol.get_xsize()
	if (last_ring == -1):	last_ring = nx//2 - 2

	print_msg("Outer radius                : %i\n"%(last_ring))
	print_msg("Ring step                   : %i\n"%(rstep))
	print_msg("X search range              : %s\n"%(xrng))
	print_msg("Y search range              : %s\n"%(yrng))
	print_msg("Translational step          : %s\n"%(step))
	print_msg("Angular step                : %s\n"%(delta))
	print_msg("Angular search range        : %s\n"%(an))
	print_msg("Maximum iteration           : %i\n"%(max_iter))
	print_msg("Center type                 : %i\n"%(center))
	print_msg("CTF correction              : %s\n"%(CTF))
	print_msg("Reference projection method : %s\n"%(ref_a))
	print_msg("Symmetry group              : %s\n\n"%(sym))

	if (maskfile) :
		if (type(maskfile) is types.StringType): mask3D = get_image(maskfile)
		else                                  : mask3D = maskfile
	else          :   mask3D = model_circle(last_ring, nx, nx, nx)
	mask = model_circle(last_ring, nx, nx)
	#drop_image(vol, os.path.join(outdir,"ref_vol00.hdf"))

	del ima
	data = EMData.read_images(stack)
	nima = len(data)
	for im in xrange(nima):
		data[im].set_attr('ID', im)
		if(CTF):
			if(data[im].get_attr_default('ctf_applied', 2) > 0):
				ERROR("data cannot be ctf-applied","ali2d_a",1)
			st = Util.infomask(data[im], mask, False)
			data[im] -= st[0]

	# initialize data for the reference preparation function
	ref_data = [mask3D, ceter, None, None]

	# do the projection matching
	for N_step in xrange(lstp):
 		for Iter in xrange(max_iter):
			print_msg("ITERATION #%3d\n"%(N_step*max_iter + Iter+1))

			if(an[N_step] == -1):	peak, pixel_error = proj_ali_incore(vol, mask3D, data, first_ring, last_ring, rstep, xrng[N_step], yrng[N_step], step[N_step], delta[N_step], ref_a, sym, CTF, MPI=False)
			else:	                peak, pixel_error = proj_ali_incore_local(vol, mask3D, data, first_ring, last_ring, rstep, xrng[N_step], yrng[N_step], step[N_step], delta[N_step], an[N_step], ref_a, sym, CTF, MPI=False)
			#  3D stuff
			if(CTF): vol1 = recons3d_4nn_ctf(data, range(0,nima,2), 1.0e6, 1, sym)
			else:	 vol1 = recons3d_4nn(data, range(0,nima,2), sym)
			if(CTF): vol2 = recons3d_4nn_ctf(data, range(1,nima,2), 1.0e6, 1, sym)
			else:	 vol2 = recons3d_4nn(data, range(1,nima,2), sym)

			fscc = fsc_mask(vol1, vol2, mask3D, 1.0, os.path.join(outdir, "resolution%04d"%(N_step*max_iter+Iter+1)))
			del vol1
			del vol2

			# calculate new and improved 3D
			if(CTF): vol = recons3d_4nn_ctf(data, range(nima), 1.0e6, 1, sym)
			else:	 vol = recons3d_4nn(data, range(nima), sym)
			# store the reference volume
			drop_image(vol, os.path.join(outdir, "vol%04d.hdf"%(N_step*max_iter+Iter+1)))
			did, vav = recons3d_nn_SSNR(data,  mask2D = None, ring_width=1, npad =1, sign=1, symmetry = sym, CTF =CTF)
			vol = fft(Util.divn_filter(fft(vol), vav))
			ref_data[2] = vol
			ref_data[3] = fscc
			#  call user-supplied function to prepare reference image, i.e., center and filter it
			vol, cs = user_func( ref_data )
			if center == 1:
				from utilities import rotate_3D_shift
				rotate_3D_shift(data, [-cs[0], -cs[1], -cs[2]])
			drop_image(vol, os.path.join(outdir, "volf%04d.hdf"%(N_step*max_iter+Iter+1)))
	#  here we  write header info
	from utilities import write_headers
	write_headers( stack, data, range(nima))
	print_end_msg("ali3d_a")

'''

'''
def Xehelix_MPI_test(stack, ref_vol, outdir, delta, psi_max, search_rng, range, ywobble, pixel_size, dp, dphi, fract, rmax, rmin, FindPsi = True, maskfile = None, \
	    maxit = 1, CTF = False, snr = 1.0, sym = "c1",  user_func_name = "helical", npad = 2, debug = False):

	from alignment       import Numrinit, prepare_refrings, proj_ali_incore, proj_ali_incore_local, proj_ali_incore_local_psi
	from alignment       import ringwe, ang_n
	from utilities       import model_circle, get_image, drop_image, get_input_from_string, peak_search, model_cylinder, pad, model_blank
	from utilities       import bcast_list_to_all, bcast_number_to_all, reduce_EMData_to_root, bcast_EMData_to_all
	from utilities       import send_attr_dict, sym_vol, get_input_from_string
	from utilities       import get_params_proj, set_params_proj, file_type, compose_transform2
	from fundamentals    import rot_avg_image, ccf, fft, rot_shift2D
	import os
	import types
	from utilities       import print_begin_msg, print_end_msg, print_msg, chunks_distribution
	from mpi             import mpi_bcast, mpi_comm_size, mpi_comm_rank, MPI_FLOAT, MPI_COMM_WORLD, mpi_barrier
	from mpi             import mpi_reduce, MPI_INT, MPI_SUM
	from filter          import filt_ctf
	from projection      import prep_vol, prgs
	#from statistics      import hist_list, varf3d_MPI, fsc_mask
	from applications	 import MPI_start_end, header
	from pixel_error     import ordersegments
	from time            import time
	from copy 			 import copy
	from math 			 import sqrt

	def rot2pad(imi, alpha=0.0, sx=0.0, sy=0.0):
		from utilities    import pad
		from fundamentals import rot_shift2D
		lnx = imi.get_xsize()
		lny = imi.get_ysize()
		ln = max(lnx,lny)
		if lnx == lny: return rot_shift2D(imi,alpha,sx,sy)
		else:          return Util.window(rot_shift2D(pad(imi,ln,ln,1,"circumference"), alpha,sx,sy), lnx, lny,1, 0,0,0)


	nproc     = mpi_comm_size(MPI_COMM_WORLD)
	myid      = mpi_comm_rank(MPI_COMM_WORLD)
	main_node = 0

	search_rng   = int(search_rng)
	range        = int(range)
	if(pixel_size < 0.0 or dp < 0.0 ):  ERROR('Helical symmetry parameters have to be provided', "ehelix_MPI", 1, myid)

	if os.path.exists(outdir):  ERROR('Output directory %s  exists, please change the name and restart the program'%outdir, "ehelix_MPI", 1, myid)
	mpi_barrier(MPI_COMM_WORLD)

	if myid == main_node:
		os.mkdir(outdir)
		import global_def
		global_def.LOGFILE =  os.path.join(outdir, global_def.LOGFILE)
		print_begin_msg("ehelix_MPI")
	mpi_barrier(MPI_COMM_WORLD)

	if debug:
		from time import sleep
		while not os.path.exists(outdir):
			print  "Node ",myid,"  waiting..."
			sleep(5)

		info_file = os.path.join(outdir, "progress%04d"%myid)
		finfo = open(info_file, 'w')
	else:
		finfo = None
	max_iter = int(maxit)
	
	vol     = EMData()
	vol.read_image(ref_vol)
	nx      = vol.get_xsize()
	ny      = vol.get_ysize()
	nz      = vol.get_zsize()

	if(sym[0] == "d"  or sym[0] == "D"):  Dsym = True
	else:                                 Dsym = False

	#  For the time being only one delta!!!
	delta       = get_input_from_string(delta)[0]

	if myid == main_node:
		import user_functions
		user_func = user_functions.factory[user_func_name]

		print_msg("Input stack                               : %s\n"%(stack))
		print_msg("Reference volume                          : %s\n"%(ref_vol))	
		print_msg("Output directory                          : %s\n"%(outdir))
		print_msg("Maskfile                                  : %s\n"%(maskfile))
		print_msg("Angular step                              : %s\n"%(delta))
		print_msg("Search for psi                            : %s\n"%(FindPsi))
		if FindPsi:  print_msg("Maximum range for psi search              : %s\n"%(psi_max))
		print_msg("X-search range                            : %f\n"%(search_rng))
		print_msg("X-search wobble                           : %f\n"%(range))
		print_msg("Pixel size [A]                            : %f\n"%(pixel_size))
		print_msg("dp [A]                                    : %f\n"%(dp))
		print_msg("dphi                                      : %f\n"%(dphi))
		print_msg("Maximum iteration                         : %i\n"%(max_iter))
		print_msg("CTF correction                            : %s\n"%(CTF))
		print_msg("Signal-to-Noise Ratio                     : %f\n"%(snr))
		print_msg("Symmetry group                            : %s\n\n"%(sym))

	if maskfile:
		if type(maskfile) is types.StringType: mask3D = get_image(maskfile)
		else:                                  mask3D = maskfile
	else: mask3D = model_cylinder(rmax, nx, ny, nz)

	fscmask = mask3D
	if CTF:
		from reconstruction import recons3d_4nn_ctf_MPI
		from filter         import filt_ctf
	else:	 from reconstruction import recons3d_4nn_MPI

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
	print  " READ IMAGES ", myid,nima,nproc

	rise = int(dp/pixel_size)

	data_nx = data[0].get_xsize()
	data_ny = data[0].get_ysize()
	if(nx < data_nx):
		data_nx = nx
		for im in xrange(nima):  data[im]=Util.window(data[im], data_nx, data_ny, 1, 0, 0, 0)
	data_nn = max(data_nx, data_ny)
	mask2D  = pad(model_blank(2*int(rmax), data_ny, 1, 1.0), data_nx, data_ny, 1, 0.0)
	fdata = [None]*nima
	for im in xrange(nima):
		data[im].set_attr('ID', list_of_particles[im])
		st = Util.infomask(data[im], mask2D, False)
		data[im] -= st[0]
		data[im] /= st[1]
		if CTF:
			qctf = data[im].get_attr("ctf_applied")
			if qctf == 0:
				ctf_params = data[im].get_attr("ctf")
				data[im] = filt_ctf(data[im], ctf_params)
				data[im].set_attr('ctf_applied', 1)
			elif qctf != 1:
				ERROR('Incorrectly set qctf flag', "ehelix_MPI", 1,myid)
		#  if FindPsi,  apply the angle to data[im], do fft and put in fdata[im]
		if FindPsi:
			phi,theta,psi,tsx,tsy = get_params_proj(data[im])
			if( theta != 0.0):
				if(abs(psi - 90.) < abs(psi - 270.0)):  gamma =  90.0
				else:                                   gamma = 270.0
				fdata[im] = fft( rot2pad(data[im], gamma-psi) )
			else:
				set_params_proj(data[im], [0.0, 90.0, 90.0, 0.0,0.0])
				fdata[im] = fft( data[im] )				
		else:
			set_params_proj(data[im], [0.0, 90.0, 90.0, 0.0,0.0])
			fdata[im] = fft( data[im] )
	del list_of_particles

	if debug:
		finfo.write( '%d loaded  \n' % nima )
		finfo.flush()
	if myid == main_node:
		# initialize data for the reference preparation function
		# No centering for helical reconstruction
		ref_data = [None, mask3D, None, None, None ]

	phiwobble = int(float(ywobble)/rise*dphi/delta+0.5)  # phiwobble is NOT in degrees, it is in nphi units

	nwx = 2*search_rng+3
	nwy = rise+2*ywobble+2
	nwxc = nwx//2
	nwyc = nwy//2
	nphi = int(360.0/delta + 0.5)
	#print  "  params  ",nwx,nwy,nwxc,nwyc,nphi
	if FindPsi:
		mode = "F"
		cnx = data_nn//2+1
		cny = cnx
		numr = Numrinit(1, data_nn//2-2, 1, mode)
		wr   = ringwe(numr, mode)
		maxrin = numr[len(numr)-1]
		crefim = [None]*nphi
	else:
		#  have to initialize them, otherwise problem with passing the arguments
		mode = "F"
		cnx = data_nx//2+1
		cny = cnx
		numr = []
		wr   = []
		maxrin = 0
		crefim = []

	terminate = 0
	Iter = 0
 	while Iter < max_iter and terminate == 0:
		Iter += 1
		if myid == main_node:
			start_time = time()
			print_msg("\nITERATION #%3d\n"%(Iter))
		volft, kbx, kby, kbz = prep_vol( vol )
		del vol

		refproj = [None]*nphi
		if( not Dsym):  rotproj = [None]*nphi
		else:           rotproj = []

		#refproj = EMData.read_images("refstraight45.hdf")   #   Here
		for iphi in xrange(nphi):
			refproj[iphi] = Util.window(  prgs( volft, kbz, [delta*iphi, 90.0, 90.0, 0.0, 0.0], kbx, kby ), data_nx, nz,1, 0,0,0)
			#refproj[iphi].write_image("refstraight45.hdf",iphi)
			st = Util.infomask(refproj[iphi] , mask2D, True)
			refproj[iphi] -= st[0]
			refproj[iphi] /= st[1]
			refproj[iphi] = Util.muln_img(refproj[iphi], mask2D )

			if FindPsi:
				temp = Util.Polar2Dm(pad(refproj[iphi], data_nn, data_nn, 1, "circumference"), cnx, cny, numr, mode)
				Util.Frngs(temp, numr)
				Util.Applyws(temp, numr, wr)
				crefim[iphi] = temp
			#  rotated in-plane by 180 are equivalent to rot_shift3D(vol,-90,180.0,90) with phi running as phi
			if(not Dsym):  rotproj[iphi] = fft( rot2pad(refproj[iphi],180.0) )
			refproj[iphi] = fft( refproj[iphi] )
		#exit()
		#if myid == main_node:  
		astart_time = time()
		for ifil in xrange(nfils):
			if myid == main_node:  start_time = time()
			if myid == main_node:
				print_msg("Process filament %4d %d\n"%(ifil,time()-start_time));start_time = time()
			ldata = [pad(data[im], data_nn, data_nn, 1, "circumference") for im in xrange(indcs[ifil][0],indcs[ifil][1])]
			Util.constrained_helix_test(ldata, fdata[indcs[ifil][0]:indcs[ifil][1]], refproj, rotproj, [float(dp), float(dphi), float(rise), float(delta)], [int(nphi), int(phiwobble), int(range), int(ywobble), int(Dsym), int(nwx), int(nwy), int(nwxc), int(nwyc)], FindPsi, float(psi_max), crefim, numr, int(maxrin), mode, int(cnx), int(cny))
			#constrained_helix     (ldata, fdata[indcs[ifil][0]:indcs[ifil][1]], refproj, rotproj,  dp, dphi, rise, delta ,  nphi, phiwobble, range, ywobble, Dsym, nwx, nwy, nwxc, nwyc , FindPsi, psi_max, crefim, numr, maxrin, mode, cnx, cny, myid, main_node)
			#constrained_helix     (data[indcs[ifil][0]:indcs[ifil][1]], fdata[indcs[ifil][0]:indcs[ifil][1]], refproj, rotproj,  dp, dphi, rise, delta ,  nphi, phiwobble, range, ywobble, Dsym, nwx, nwy, nwxc, nwyc , FindPsi, psi_max, crefim, numr, maxrin, mode, cnx, cny, myid, main_node)
			for im in xrange(indcs[ifil][0], indcs[ifil][1]):
				temp = Util.get_transform_params(ldata[im-indcs[ifil][0]], "xform.projection", "spider")
				set_params_proj(data[im],[temp["phi"],temp["theta"],temp["psi"],-temp["tx"],-temp["ty"]])
			if FindPsi:
				for im in xrange(indcs[ifil][0], indcs[ifil][1]):  fdata[im] = rot2pad(data[im], ldata[im-indcs[ifil][0]].get_attr("bestang"))
			#print  "Parameters computed for filament",myid,ifil,time()-start_time;start_time = time()
			if myid == main_node:
				print_msg("Parameters computed for filament %4d %d\n"%(ifil,time()-start_time));start_time = time()
		del ldata
		del refproj, volft
		if(not Dsym):  del rotproj
		print  "Time of alignment = ",myid,time()-astart_time
		mpi_barrier(MPI_COMM_WORLD)
		#if myid == main_node:
		#	print_msg("Time of alignment = %\n"%(time()-astart_time));start_time = time()

		if CTF:  vol = recons3d_4nn_ctf_MPI(myid, data, symmetry=sym, snr = snr, npad = npad, xysize = nx, zsize = nz)
		else:    vol = recons3d_4nn_MPI(myid, data, symmetry=sym, npad = npad, xysize = nx, zsize = nz)

		if myid == main_node:
			print_msg("3D reconstruction time = %d\n"%(time()-start_time));start_time = time()

		if myid == main_node:
			#vol.write_image(os.path.join(outdir, "vol%03d.hdf"%Iter))
			vol = vol.helicise(pixel_size, dp, dphi, fract, rmax, rmin)
			vol = sym_vol(vol, symmetry=sym)
			ref_data[0] = vol
			vol = user_func(ref_data)
			vol = vol.helicise(pixel_size, dp, dphi, fract, rmax, rmin)
			vol = sym_vol(vol, symmetry=sym)
			vol.write_image(os.path.join(outdir, "volf%03d.hdf"%Iter))
			#if(Iter == max_iter-1):  drop_image(vol, os.path.join(outdir, "volfshift.hdf"))

		bcast_EMData_to_all(vol, myid, main_node)
		# write out headers, under MPI writing has to be done sequentially
		mpi_barrier(MPI_COMM_WORLD)
		par_str = ['xform.projection', 'ID']
		if myid == main_node:
			start_time = time()
			if(file_type(stack) == "bdb"):
				from utilities import recv_attr_dict_bdb
				recv_attr_dict_bdb(main_node, stack, data, par_str, 0, nima, nproc)
			else:
				from utilities import recv_attr_dict
				recv_attr_dict(main_node, stack, data, par_str, 0, nima, nproc)
			print_msg("Time to write header information= %d\n"%(time()-start_time))
			start_time = time()
		else:		send_attr_dict(main_node, data, par_str, 0, nima)
		if myid == main_node:
			# write params to text file
			header(stack, params='xform.projection', fexport=os.path.join(outdir, "parameters%04d.txt"%Iter))
	if myid == main_node: print_end_msg("ehelix_MPI")
'''
def consistency(params, coords,pixel_size, dp, dphi):
	''' 
	determine consistency of parameters of segments belonging to ONE filament
	and dst is distance between adjacent segments in pixels
	
	params: parameters for segments in this filament
	coords: coordinates for segments in this filament
	'''
	from development import predict
	from pixel_error   import angle_diff
	from copy import copy

	def get_dist(c1, c2):
		from math import sqrt
		d = sqrt((c1[0] - c2[0])**2 + (c1[1] - c2[1])**2)
		return d

	dpp = dp/pixel_size
	ns = len(params)
	phig = [0.0]*ns
	yg = [0.0]*ns
	for j in xrange(ns): 
		phig[j] = params[j][0]
		yg[j] = params[j][4]

	distances = [0.0]*ns
	for i in xrange(1,ns):  distances[i] = get_dist( coords[i], coords[0] )

	terr = 1.e23
	consphi=[0.0]*ns
	for idir in xrange(-1,2,2):
		phierr = []
		#  get phi's
		ddphi = pixel_size/dp*idir*dphi
		phis = [0.0]*ns
		for i in xrange(ns):
			yy = distances[i] #+ yg[i]
			phis[i] = (yy*ddphi)%360.0
		# find the overall angle
		angdif = angle_diff(phis,phig)
		#print " angdif ",angdif
		lerr = 0.0
		for i in xrange(ns):
			anger = (phis[i]+angdif - phig[i] + 360.0)%360.0
			if( anger > 180.0 ): anger -= 360.0
			lerr += abs(anger)
		if(lerr < terr):
			terr = lerr
			for i in xrange(ns):
				cphi = (phis[i]+angdif)%360.0
				consphi[i] = cphi

	consy = sum(yg)/ns
	err = 0.0
	for i in xrange(ns):
		yerr  = (abs(yg[i] - consy)/dpp)
		err += yerr
		dv = (abs(consphi[i] - phig[i]))%360.
		dv = min(dv, 360.0-dv)
		err += dv/dphi
	
	avgerr = err/ns
	return avgerr

# Calculate averages per filament to a given stack (wrap for ave_var in statistics)
def ave_ali_filament(name_stack, name_out = None, ali = False, active = False, param_to_save_size = None, set_as_member_id = None):
	from statistics 	import ave_var, add_ave_varf, k_means_list_active
	from utilities  	import file_type
	from pixel_error  	import ordersegments
	"""
	   This function is called by sxave_ali.py
	"""
	N = EMUtil.get_image_count(name_stack)
	if ali:
		mode = 'a'
	else:
		mode = ''
	
	infils = EMUtil.get_all_attributes(name_stack, "filament")
	ptlcoords = EMUtil.get_all_attributes(name_stack, 'ptcl_source_coord')
	filaments = ordersegments(infils, ptlcoords)
	total_nfils = len(filaments)
	
	for i in xrange(total_nfils):
		fildata = EMData.read_images(name_stack, filaments[i])
		ave, var = ave_var(fildata, mode)
		ave.set_attr('members', filaments[i])
		ext = file_type(name_stack)
		if name_out is None:
			if ext == 'bdb': name = name_stack.split(':')[1] + '.hdf'
			else:            name = name_stack
			ave.write_image('ave_' + name, i)
		else:
			ave.write_image(name_out, i)

def sum_lines(img):
	from utilities import model_blank
	nx  = img.get_xsize()
	ny  = img.get_ysize()

	b = model_blank(nx)

	for j in xrange(ny):  Util.add_img(b, Util.window(img,nx,1,1,0,j-ny//2,0))

	b/=ny
	return b

def symsearch_MPI(ref_vol, outdir, maskfile, dp, ndp, dp_step, dphi, ndphi, dphi_step,\
	rmin, rmax, fract, sym, user_func_name, datasym,\
	pixel_size, debug):

	from alignment      import Numrinit, prepare_refrings, proj_ali_helical, proj_ali_helical_90, proj_ali_helical_local, proj_ali_helical_90_local, helios,helios7
	from utilities      import model_circle, get_image, drop_image, get_input_from_string, pad, model_blank, sym_vol
	from utilities      import bcast_list_to_all, bcast_number_to_all, reduce_EMData_to_root, bcast_EMData_to_all
	from utilities      import send_attr_dict
	from utilities      import get_params_proj, set_params_proj, file_type
	from fundamentals   import rot_avg_image
	from pixel_error    import max_3D_pixel_error
	import os
	import types
	from utilities      import print_begin_msg, print_end_msg, print_msg
	from mpi            import mpi_bcast, mpi_comm_size, mpi_comm_rank, MPI_FLOAT, MPI_COMM_WORLD, mpi_barrier
	from mpi            import mpi_recv,  mpi_send, MPI_TAG_UB
	from mpi            import mpi_reduce, MPI_INT, MPI_SUM
	from filter         import filt_ctf
	from projection     import prep_vol, prgs
	from statistics     import hist_list, varf3d_MPI
	from applications   import MPI_start_end
	from EMAN2 import Vec2f
	from string    import lower,split
	from math import cos, pi

	number_of_proc = mpi_comm_size(MPI_COMM_WORLD)
	myid           = mpi_comm_rank(MPI_COMM_WORLD)
	main_node = 0

	if myid == 0:
		if os.path.exists(outdir):  nx = 1
		else:  nx = 0
	else:  nx = 0
	ny = bcast_number_to_all(nx, source_node = main_node)

	if ny == 1:  ERROR('Output directory exists, please change the name and restart the program', "symsearch_MPI", 1,myid)
	mpi_barrier(MPI_COMM_WORLD)

	if myid == main_node:
		os.mkdir(outdir)
		import global_def
		global_def.LOGFILE =  os.path.join(outdir, global_def.LOGFILE)
	mpi_barrier(MPI_COMM_WORLD)

	nlprms = (2*ndp+1)*(2*ndphi+1)
	if nlprms< number_of_proc:
		ERROR('number of CPUs is larger than the number of helical search, please reduce it or at this moment modify ndp,dphi in the program', "ihrsr_MPI", 1,myid)

	if debug:
		from time import sleep
		while not os.path.exists(outdir):
			print  "Node ",myid,"  waiting..."
			sleep(5)

		info_file = os.path.join(outdir, "progress%04d"%myid)
		finfo = open(info_file, 'w')
	else:
		finfo = None
		
	vol     = EMData()
	vol.read_image(ref_vol)
	
	if myid == main_node:
		import user_functions
		user_func = user_functions.factory[user_func_name]

		print_msg("Reference volume                          : %s\n"%(ref_vol))	
		print_msg("Output directory                          : %s\n"%(outdir))
		print_msg("Maskfile                                  : %s\n"%(maskfile))
		print_msg("min radius for helical search (in pix)    : %5.4f\n"%(rmin))
		print_msg("max radius for helical search (in pix)    : %5.4f\n"%(rmax))
		print_msg("fraction of volume used for helical search: %5.4f\n"%(fract))
		print_msg("initial symmetry - angle                  : %5.4f\n"%(dphi))
		print_msg("initial symmetry - axial rise             : %5.4f\n"%(dp))
		print_msg("symmetry output doc file                  : %s\n"%(datasym))
		print_msg("User function                             : %s\n"%(user_func_name))

	if maskfile:
		if type(maskfile) is types.StringType: mask3D = get_image(maskfile)
		else:                                  mask3D = maskfile
	else: mask3D = None
	#else: mask3D = model_circle(last_ring, nx, nx, nx)

	from time import time
	#from filter import filt_gaussl
	#vol = filt_gaussl(vol, 0.25)
	start_time = time()
	if myid == main_node:
		lprms = []
		for i in xrange(-ndp,ndp+1,1):
			for j in xrange(-ndphi,ndphi+1,1):
				lprms.append( dp   + i*dp_step)
				lprms.append( dphi + j*dphi_step)
		#print "lprms===",lprms
		recvpara = []
		for im in xrange(number_of_proc):
			helic_ib, helic_ie = MPI_start_end(nlprms, number_of_proc, im)
			recvpara.append(helic_ib )
			recvpara.append(helic_ie )

	para_start, para_end = MPI_start_end(nlprms, number_of_proc, myid)

	list_dps     = [0.0]*((para_end-para_start)*2)
	list_fvalues = [-1.0]*((para_end-para_start)*1)

	if myid == main_node:
		for n in xrange(number_of_proc):
			if n!=main_node: mpi_send(lprms[2*recvpara[2*n]:2*recvpara[2*n+1]], 2*(recvpara[2*n+1]-recvpara[2*n]), MPI_FLOAT, n, MPI_TAG_UB, MPI_COMM_WORLD)
			else:    list_dps = lprms[2*recvpara[2*0]:2*recvpara[2*0+1]]
	else:
		list_dps = mpi_recv((para_end-para_start)*2, MPI_FLOAT, main_node, MPI_TAG_UB, MPI_COMM_WORLD)

	list_dps = map(float, list_dps)

	local_pos = [0.0, 0.0, -1.0e20]
	for i in xrange(para_end-para_start):
		fvalue = helios7(vol, pixel_size, list_dps[i*2], list_dps[i*2+1], fract, rmax, rmin)
		if(fvalue >= local_pos[2]):
			local_pos = [list_dps[i*2], list_dps[i*2+1], fvalue ]
	if myid == main_node:
		list_return = [0.0]*(3*number_of_proc)
		for n in xrange(number_of_proc):
			if n != main_node: list_return[3*n:3*n+3]                 = mpi_recv(3,MPI_FLOAT, n, MPI_TAG_UB, MPI_COMM_WORLD)
			else:              list_return[3*main_node:3*main_node+3]  = local_pos[:]
	else:
		mpi_send(local_pos, 3, MPI_FLOAT, main_node, MPI_TAG_UB, MPI_COMM_WORLD)

	if myid == main_node:	
		maxvalue = list_return[2]
		for i in xrange(number_of_proc):
			if( list_return[i*3+2] >= maxvalue ):
				maxvalue = list_return[i*3+2]
				dp       = list_return[i*3+0]
				dphi     = list_return[i*3+1]
		dp   = float(dp)
		dphi = float(dphi)
		#print  "  GOT dp dphi",dp,dphi

		vol  = vol.helicise(pixel_size, dp, dphi, fract, rmax, rmin)
		drop_image(vol, os.path.join(outdir, "vol.hdf"))

		print_msg("New delta z and delta phi      : %s,    %s\n\n"%(dp,dphi))		
		
		
	if(myid==main_node):
		fofo = open(os.path.join(outdir,datasym),'a')
		fofo.write('  %12.4f   %12.4f\n'%(dp,dphi))
		fofo.close()
		vol = sym_vol(vol, symmetry=sym)
		ref_data = [vol, mask3D]
		#if  fourvar:  ref_data.append(varf)
		vol = user_func(ref_data)
		vol = vol.helicise(pixel_size, dp, dphi, fract, rmax, rmin)
		vol = sym_vol(vol, symmetry=sym)
		drop_image(vol, os.path.join(outdir, "volf.hdf"))
		print_msg("\nSymmetry search and user function time = %d\n"%(time()-start_time))
		start_time = time()
	
	# del varf
	if myid == main_node: print_end_msg("symsearch_MPI")


def nlocal_ali3d_MPI(stack, outdir, maskfile, ou = -1,  delta = 2, ts=0.25, center = -1, maxit = 10, 
                CTF = False, snr = 1.0, sym = "c1", chunk = -1.0, user_func_name = "ref_ali3d",
		    	fourvar = True, npad = 4, debug = False):
	fourvar = False
	"""
		
	"""
	from alignment        import eqproj_cascaded_ccc
	from filter           import filt_ctf
	from projection       import prep_vol
	from utilities        import bcast_number_to_all, model_circle, get_params_proj, set_params_proj
	from utilities        import bcast_EMData_to_all, bcast_list_to_all, send_attr_dict
	from utilities        import get_image, drop_image, file_type
	from utilities        import amoeba_multi_level, rotate_3D_shift, estimate_3D_center_MPI
	from utilities        import print_begin_msg, print_end_msg, print_msg
	from reconstruction   import rec3D_MPI, rec3D_MPI_noCTF
	from statistics       import varf3d_MPI
	from math             import pi
	from mpi              import mpi_bcast, mpi_comm_size, mpi_comm_rank, MPI_FLOAT, MPI_COMM_WORLD, mpi_barrier
	from mpi              import mpi_reduce, MPI_INT, MPI_SUM
	from applications     import MPI_start_end
	from EMAN2 import Processor
	import os
	import sys
	from EMAN2 import Vec2f


	number_of_proc = mpi_comm_size(MPI_COMM_WORLD)
	myid = mpi_comm_rank(MPI_COMM_WORLD)

	if CTF:
		from filter import filt_ctf

	main_node = 0
	
	if os.path.exists(outdir):  ERROR('Output directory exists, please change the name and restart the program', "local_ali3d_MPI ", 1,myid)
	mpi_barrier(MPI_COMM_WORLD)

	if myid == main_node:
		os.mkdir(outdir)
		import global_def
		global_def.LOGFILE =  os.path.join(outdir, global_def.LOGFILE)
		print_begin_msg("local_ali3d_MPI")
		import user_functions
		user_func = user_functions.factory[user_func_name]
		if CTF:
			ima = EMData()
			ima.read_image(stack, 0)
			ctf_applied = ima.get_attr("ctf_applied")
			del ima
			if ctf_applied == 1:  ERROR("local_ali3d does not work for CTF-applied data", "local_ali3d_MPI", 1,myid)
	mpi_barrier(MPI_COMM_WORLD)

	if debug:
		from time import sleep
		while not os.path.exists(outdir):
			print  "Node ",myid,"  waiting..."
			sleep(5)

		finfo = open(os.path.join(outdir, "progress%04d"%myid), 'w')
	else:
		finfo = None

	last_ring   = int(ou)
	max_iter    = int(maxit)
	center      = int(center)

	if myid == main_node:
		if(file_type(stack) == "bdb"):
			from EMAN2db import db_open_dict
			dummy = db_open_dict(stack, True)
		# horatio active_refactoring Jy51i1EwmLD4tWZ9_00000_1
		# active = EMUtil.get_all_attributes(stack, 'active')
		# list_of_particles = []
		# for im in xrange(len(active)):
		# 	if(active[im]):  list_of_particles.append(im)
		# del active
		# nima = len(list_of_particles)
		
		nima = EMUtil.get_image_count(stack)
		list_of_particles = range(nima)
	
		ima     = EMData()
		ima.read_image(stack, 0)
		nx      = ima.get_xsize()
		del ima
	else:
		nima = 0
		nx = 0
	total_nima = bcast_number_to_all(nima, source_node = main_node)
	nx = bcast_number_to_all(nx, source_node = main_node)

	if myid != main_node:
		list_of_particles = [-1]*total_nima
	list_of_particles = bcast_list_to_all(list_of_particles, source_node = main_node)

	image_start, image_end = MPI_start_end(total_nima, number_of_proc, myid)
	# create a list of images for each node
	list_of_particles = list_of_particles[image_start: image_end]
	nima = len(list_of_particles)

	if last_ring < 0:	last_ring = int(nx/2) - 2

	if chunk <= 0.0:  chunk = 1.0
	n_of_chunks = int(1.0/chunk)

	if myid == main_node:
		import user_functions
		user_func = user_functions.factory[user_func_name]

		print_msg("Input stack                 : %s\n"%(stack))
		print_msg("Output directory            : %s\n"%(outdir))
		print_msg("Maskfile                    : %s\n"%(maskfile))
		print_msg("Outer radius                : %i\n"%(last_ring))
		print_msg("Angular search range        : %s\n"%(delta))
		print_msg("Shift search range          : %f\n"%(ts))
		print_msg("Maximum iteration           : %i\n"%(max_iter))
		print_msg("Center type                 : %i\n"%(center))
		print_msg("CTF correction              : %s\n"%(CTF))
		print_msg("Signal-to-Noise Ratio       : %f\n"%(snr))
		print_msg("Symmetry group              : %s\n"%(sym))
		print_msg("Chunk size                  : %f\n\n"%(chunk))

	if maskfile:
		import  types
		if type(maskfile) is types.StringType:  mask3D = get_image(maskfile)
		else:                                  mask3D = maskfile
	else:
		mask3D = model_circle(last_ring, nx, nx, nx)
	mask2D = model_circle(last_ring, nx, nx)

	if debug:
		finfo.write( "image_start, image_end: %d %d\n" %(image_start, image_end) )
		finfo.flush()

	if debug:
		finfo.write("  chunk = "+str(chunk)+"   ")
		finfo.write("\n")
		finfo.flush()
		finfo.write("  Number of chunks = "+str(n_of_chunks)+"   ")
		finfo.write("\n")
		finfo.flush()

	dataim = EMData.read_images(stack, list_of_particles)
	for im in xrange(len(dataim)):
		dataim[im].set_attr('ID', list_of_particles[im])
	del list_of_particles

	if debug:
		finfo.write("  First image on this processor: "+str(image_start)+"   ")
		finfo.write("  Last  image on this processor: "+str(image_end)+"   ")
		finfo.write("\n")
		finfo.flush()

	if myid == main_node:
		# initialize data for the reference preparation function
		ref_data = [ mask3D, max(center,0), None, None, None ]
		# for method -1, switch off centering in user function
		ref_data.append( None )

	from time import time	
		


	disps = []
	recvcount = []
	for im in xrange(number_of_proc):
		if( im == main_node ):  disps.append(0)
		else:                  disps.append(disps[im-1] + recvcount[im-1])
		ib, ie = MPI_start_end(total_nima, number_of_proc, im)
		recvcount.append( ie - ib )

	#  this is needed for gathering of pixel errors
	pixer = [0.0]*nima
	#  data: 0 - volf , 1 - kb, 2 - image, 3 - mask2D, 4 - , 5 - params, 6 - 
	data = [None]*7
	data[3] = mask2D
	cs = [0.0]*3

	for iteration in xrange(maxit+1):
		if myid == main_node:
			start_time = time()
			print_msg("ITERATION #%3d\n"%(iteration+1))
		if debug:
			finfo.write("  iteration = "+str(iteration)+"   ")
			finfo.write("\n")
			finfo.flush()
		for ic in xrange(n_of_chunks):
			if(center == -1):
				if debug:
					finfo.write("  begin centering \n")
					finfo.flush()
				cs[0], cs[1], cs[2], dummy, dummy = estimate_3D_center_MPI(dataim, total_nima, myid, number_of_proc, main_node)
				cs = mpi_bcast(cs, 3, MPI_FLOAT, main_node, MPI_COMM_WORLD)
				cs = [-float(cs[0]), -float(cs[1]), -float(cs[2])]
				rotate_3D_shift(dataim, cs)
				if myid == main_node:
					msg = "Average center x = %10.3f        Center y = %10.3f        Center z = %10.3f\n"%(cs[0], cs[1], cs[2])
					print_msg(msg)
					print_msg("Time to center = %d\n"%(time()-start_time))
					start_time = time()
			# compute updated 3D before each chunk
 	    		# resolution
			if debug:
				finfo.write("  begin reconstruction = "+str(image_start))
				finfo.write("\n")
				finfo.flush()

			if CTF: vol, fscc = rec3D_MPI(dataim, snr, sym, mask3D, os.path.join(outdir, "resolution%03d_%03d"%(iteration, ic)), myid, main_node, npad = npad)
			else:   vol, fscc = rec3D_MPI_noCTF(dataim, sym, mask3D, os.path.join(outdir, "resolution%03d_%03d"%(iteration, ic)), myid, main_node, npad = npad)

			if myid == main_node:
				drop_image(vol, os.path.join(outdir, "vol%03d_%03d.hdf"%(iteration, ic) ))
				print_msg("3D reconstruction time = %d\n"%(time()-start_time))
				start_time = time()
			if debug:
				finfo.write("  done reconstruction = "+str(image_start))
				finfo.write("\n")
				finfo.flush()

			if fourvar:
			#  Compute Fourier variance
				varf = varf3d_MPI(dataim, ssnr_text_file = os.path.join(outdir, "ssnr%03d_%03d"%(iteration, ic)), mask2D = None, reference_structure = vol, ou = ou, rw = 1.0, npad = 1, CTF = CTF, sign = 1, sym =sym, myid = myid)
				if myid == main_node:
					varf = 1.0/varf
					print_msg("Time to calculate 3D Fourier variance= %d\n"%(time()-start_time))
					start_time = time()

			else:  varf = None
			if myid == main_node:
				ref_data[2] = vol
				ref_data[3] = fscc
				ref_data[4] = varf
				#  call user-supplied function to prepare reference image, i.e., filter it
				# When center = -1, which is by default, we use the average center method
				ref_data[1] = 0
				vol, dummy = user_func(ref_data)
				drop_image(vol, os.path.join(outdir, "volf%03d_%03d.hdf"%(iteration, ic)))
			del varf

			# in last iteration return here
			if(iteration == maxit):
				if myid == main_node: print_end_msg("local_ali3d_MPI")
				return
			bcast_EMData_to_all(vol, myid, main_node)

			if not CTF:
				data[0], data[1] = prep_vol(vol)

			image_start_in_chunk = image_start + ic*nima/n_of_chunks
			image_end_in_chunk   = image_start + (ic+1)*nima/n_of_chunks
			if debug:
				finfo.write("Chunk "+str(ic)+"   Number of images in this chunk: "+str(image_end_in_chunk-image_start_in_chunk)+"\n")
				finfo.write("First image in this chunk: "+str(image_start_in_chunk)+"   Last image in this chunk: "+str(image_end_in_chunk-1)+"\n")
				finfo.flush()
			if CTF:  previous_defocus = -1.0
			for imn in xrange(image_start_in_chunk, image_end_in_chunk):
				if CTF:
					ctf_params = dataim[imn-image_start].get_attr( "ctf" )
					if ctf_params.defocus != previous_defocus:
						previous_defocus = ctf_params.defocus
						data[0], data[1] = prep_vol(filt_ctf(vol, ctf_params))

				data[2] = dataim[imn-image_start]

				#phi, theta, psi, tx, ty = get_params_proj(dataim[imn-image_start])
				t1 = dataim[imn-image_start].get_attr("xform.projection")
				dp = t1.get_params("spider")
				data[5] = [dp["phi"], dp["theta"], dp["psi"], -dp["tx"], -dp["ty"]]
				from development import EVAL
				if debug:
					initial = EVAL(data[5], data) # this is if we need initial discrepancy
					finfo.write("Image "+str(imn)+"\n")
					finfo.write('Old  %6.1f  %6.1f  %6.1f   %5.2f  %5.2f  %11.4e\n'%(data[5][0], data[5][1], data[5][2], data[5][3], data[5][4], initial))
					#print "init  ",data[5],initial
				#from random import random
				#data[5] = [(random()-0.5)*2,(random()-0.5)*2]  #  HERE !!!!!!!!!!!
				from development import cshcca
				initialPoint = [dp["phi"], dp["theta"], dp["psi"], -dp["tx"], -dp["ty"]]
				initialStepSizes = [delta/1.2,delta/1.2, delta/1.2, ts/1.2, ts/1.2 ]
				optm_params, score = cshcca( initialPoint, initialStepSizes,  EVAL, data)
				#print "after   ",optm_params, score

				if debug:
					finfo.write('New  %6.1f  %6.1f  %6.1f   %5.2f  %5.2f  %11.4e\n'%(optm_params[0], optm_params[1], optm_params[2], optm_params[3], optm_params[4], score))
					finfo.flush()

				#from sys import exit
				#exit()
				t2 = Transform({"type":"spider","phi":optm_params[0],"theta":optm_params[1],"psi":optm_params[2]})
				t2.set_trans(Vec2f(-optm_params[3], -optm_params[4]))
				dataim[imn-image_start].set_attr("xform.projection", t2)
				from pixel_error import max_3D_pixel_error
				pixer[imn-image_start] = max_3D_pixel_error(t1, t2, last_ring)
				#set_params_proj(dataim[imn-image_start], optm_params[0])
				if( myid == main_node ):
					print_msg( "Time to process %6d particles : %d\n" % (image_end_in_chunk-image_start_in_chunk, time()-start_time) )
					start_time = time()

			# release memory of volft
			data[0] = None

			# write out headers, under MPI writing has to be done sequentially
			mpi_barrier(MPI_COMM_WORLD)
			par_str = ['xform.projection', 'ID']
			if myid == main_node:
				from utilities import file_type
				if(file_type(stack) == "bdb"):
					from utilities import recv_attr_dict_bdb
					recv_attr_dict_bdb(main_node, stack, dataim, par_str, image_start, image_end, number_of_proc)
				else:
					from utilities import recv_attr_dict
					recv_attr_dict(main_node, stack, dataim, par_str, image_start, image_end, number_of_proc)
			else:	        send_attr_dict(main_node, dataim, par_str, image_start, image_end)
			if myid == main_node:
				print_msg("Time to write header information= %d\n"%(time()-start_time))
				start_time = time()

		#output pixel errors after all headers were processed
		from mpi import mpi_gatherv
		recvbuf = mpi_gatherv(pixer, nima, MPI_FLOAT, recvcount, disps, MPI_FLOAT, main_node, MPI_COMM_WORLD)
		mpi_barrier(MPI_COMM_WORLD)
		terminate = 0
		if(myid == main_node):
			recvbuf = map(float, recvbuf)
			from statistics import hist_list
			lhist = 20
			region, histo = hist_list(recvbuf, lhist)
			if(region[0] < 0.0):  region[0] = 0.0
			msg = "      Histogram of pixel errors\n      ERROR       number of particles\n"
			print_msg(msg)
			for lhx in xrange(lhist):
				msg = " %10.3f     %7d\n"%(region[lhx], histo[lhx])
				print_msg(msg)
			# Terminate if 95% within 1 pixel error
			im = 0
			for lhx in xrange(lhist):
				if(region[lhx] > 1.0): break
				im += histo[lhx]
			if(im/float(total_nima) > 0.95):  terminate = 1
			del region, histo
		del recvbuf
		terminate = mpi_bcast(terminate, 1, MPI_INT, 0, MPI_COMM_WORLD)
		terminate = int(terminate[0])

def cshcca( initialPoint, initialStepSizes,  EVAL, data, someAcceleration=1.2, epsilon = 1.0e-3):
	npar = len(initialPoint)
	currentPoint = [initialPoint[i] for i in xrange(npar)] 
	stepSize     = [initialStepSizes[i] for i in xrange(npar)]
	acceleration = someAcceleration
	candidate = [0.0]*5
	candidate[0] = -acceleration;
	candidate[1] = -1.0 / acceleration;
	candidate[2] = 0.0;
	candidate[3] = 1.0 / acceleration;
	candidate[4] = acceleration;
	before = EVAL(data[5], data);  #  This checks value at previously best position
	#print  " inside init",data[5],before
	ita=0
	iopa=2.0
	ropa = range(npar)
	rcan = range(5)
	from random import shuffle
	for l in xrange(1000000):
		iopa /= 2.0
		stepSize     = [initialStepSizes[i]*iopa for i in xrange(npar)]
		ct = True
		while  ct:
			shuffle(ropa)
			for  i in ropa:  # randomize order
				best = 2;
				#bestScore = -1.0e23;
				shuffle(rcan)
				for j in rcan:  # randomize order
					currentPoint[i] = currentPoint[i] + stepSize[i] * candidate[j];
					temp = EVAL(currentPoint, data);
					#print " evaluated ",temp
					#if(temp > bestScore):
					if(temp > before):
						#bestScore = temp;
						before = temp;
						best = j;
						#print  "will return  ",currentPoint, before
						return currentPoint, before
					currentPoint[i] = currentPoint[i] - stepSize[i] * candidate[j];
				if candidate[best] != 0.0:
					currentPoint[i] = currentPoint[i] + stepSize[i] * candidate[best];
					stepSize[i] = stepSize[i] * candidate[best];
			qv = EVAL(currentPoint, data)
			ita+=1
			"""
			print " >>  %3d"%ita,
			for m in xrange(npar):  print  "%8.2e"%currentPoint[m],
			for m in xrange(npar):  print  "%8.2e"%stepSize[m],
			print "  %8.3e  %8.3e"%(qv,before)
			"""
			if(qv > before): before = qv
			else:  ct = False
		if( iopa < epsilon ):
			return currentPoint, qv;

def EVAL(x, data):
	from statistics import ccc
	from projection import prgs
	return ccc(prgs(data[0], data[1], x), data[2], data[3])
	


def ehelix_MPI_lastringtest(stack, ref_vol, outdir, seg_ny, delta, psi_max, search_rng, rng, ywobble, pixel_size, dp, dphi, fract, rmax, rmin, FindPsi = True, maskfile = None, \
	    maxit = 1, CTF = False, snr = 1.0, sym = "c1",  user_func_name = "helical", npad = 2, debug = False, doExhaustive=False, termprec=5.0):

	from alignment       import Numrinit, prepare_refrings, proj_ali_incore, proj_ali_incore_local, proj_ali_incore_local_psi
	from alignment       import ringwe, ang_n
	from utilities       import model_circle, get_image, drop_image, get_input_from_string, peak_search, model_cylinder, pad, model_blank
	from utilities       import bcast_list_to_all, bcast_number_to_all, reduce_EMData_to_root, bcast_EMData_to_all, write_text_file
	from utilities       import send_attr_dict, sym_vol, get_input_from_string
	from utilities       import get_params_proj, set_params_proj, file_type, compose_transform2
	from fundamentals    import rot_avg_image, ccf, fft, rot_shift2D
	import os
	import types
	from utilities       import print_begin_msg, print_end_msg, print_msg, chunks_distribution
	from mpi             import mpi_bcast, mpi_comm_size, mpi_comm_rank, MPI_FLOAT, MPI_COMM_WORLD, mpi_barrier
	from mpi             import mpi_reduce, MPI_INT, MPI_SUM
	from filter          import filt_ctf
	from projection      import prep_vol, prgs
	#from statistics      import hist_list, varf3d_MPI, fsc_mask
	from applications	 import MPI_start_end, header
	from pixel_error     import ordersegments, max_3D_pixel_error
	from time            import time
	from copy 			 import copy
	from math 			 import sqrt
	'''
	def rot2pad(imi, alpha=0.0, sx=0.0, sy=0.0):
		from utilities    import pad
		from fundamentals import rot_shift2D
		lnx = imi.get_xsize()
		lny = imi.get_ysize()
		ln = max(lnx,lny)
		if lnx == lny: return rot_shift2D(imi,alpha,sx,sy)
		else:          return Util.window(rot_shift2D(pad(imi,ln,ln,1,"circumference"), alpha,sx,sy), lnx, lny,1, 0,0,0)
	'''

	nproc     = mpi_comm_size(MPI_COMM_WORLD)
	myid      = mpi_comm_rank(MPI_COMM_WORLD)
	main_node = 0

	search_rng   = int(search_rng)
	rng        = int(rng)
	if(pixel_size < 0.0 or dp < 0.0 ):  ERROR('Helical symmetry parameters have to be provided', "ehelix_MPI", 1, myid)

	if os.path.exists(outdir):  ERROR('Output directory %s  exists, please change the name and restart the program'%outdir, "ehelix_MPI", 1, myid)
	mpi_barrier(MPI_COMM_WORLD)

	if myid == main_node:
		os.mkdir(outdir)
		import global_def
		global_def.LOGFILE =  os.path.join(outdir, global_def.LOGFILE)
		print_begin_msg("ehelix_MPI")
	mpi_barrier(MPI_COMM_WORLD)

	if debug:
		from time import sleep
		while not os.path.exists(outdir):
			print  "Node ",myid,"  waiting..."
			sleep(5)

		info_file = os.path.join(outdir, "progress%04d"%myid)
		finfo = open(info_file, 'w')
	else:
		finfo = None
	max_iter = int(maxit)
	
	vol     = EMData()
	vol.read_image(ref_vol)
	nx      = vol.get_xsize()
	ny      = vol.get_ysize()
	nz      = vol.get_zsize()
	
	# Only handle square volumes now!
	if nz < nx:
		ERROR('Do not handle square volumes .... nz cannot be less than nx', "ehelix_MPI", 1, myid)
	
	# Pad to square
	if nz > nx:
		nx = nz
		ny = nz	
		vol = pad(vol, nx, ny,nz,background=0.0)
	
	
	if(sym[0] == "d"  or sym[0] == "D"):  Dsym = True
	else:                                 Dsym = False

	#  For the time being only one delta!!!
	delta       = get_input_from_string(delta)[0]

	if myid == main_node:
		import user_functions
		user_func = user_functions.factory[user_func_name]

		print_msg("Input stack                               : %s\n"%(stack))
		print_msg("Reference volume                          : %s\n"%(ref_vol))	
		print_msg("Output directory                          : %s\n"%(outdir))
		print_msg("Maskfile                                  : %s\n"%(maskfile))
		print_msg("Angular step                              : %s\n"%(delta))
		print_msg("Search for psi                            : %s\n"%(FindPsi))
		if FindPsi:  print_msg("Maximum range for psi search              : %s\n"%(psi_max))
		print_msg("X-search range                            : %f\n"%(search_rng))
		print_msg("X-search wobble                           : %f\n"%(rng))
		print_msg("Pixel size [A]                            : %f\n"%(pixel_size))
		print_msg("dp [A]                                    : %f\n"%(dp))
		print_msg("dphi                                      : %f\n"%(dphi))
		print_msg("Maximum iteration                         : %i\n"%(max_iter))
		print_msg("CTF correction                            : %s\n"%(CTF))
		print_msg("Signal-to-Noise Ratio                     : %f\n"%(snr))
		print_msg("Symmetry group                            : %s\n"%(sym))
		print_msg("seg_ny                        		     : %i\n\n"%(seg_ny))
		print_msg("termprec                       		     : %f\n\n"%(termprec))

	if maskfile:
		if type(maskfile) is types.StringType: mask3D = get_image(maskfile)
		else:                                  mask3D = maskfile
	else: mask3D = model_cylinder(rmax, nx, ny, nz)

	fscmask = mask3D
	if CTF:
		from reconstruction import recons3d_4nn_ctf_MPI
		from filter         import filt_ctf
	else:	 from reconstruction import recons3d_4nn_MPI

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

	#  this is needed for gathering of pixel errors
	disps = []
	recvcount = [] #recvcount[i] is number of images on proc i
	total_nima = EMUtil.get_image_count(stack)
	if myid == main_node:
		maptoimgID = [0]*total_nima # maptoimgID[i] is the element in recvbuf corresponding to img ID i
		counter = 0
	for im in xrange(nproc):
		if im == main_node :  disps.append(0)
		else:                 disps.append(disps[im-1] + recvcount[im-1])
		
		temp = chunks_distribution([[len(filaments[i]), i] for i in xrange(len(filaments))], nproc)[im:im+1][0]
		filaments_im = [filaments[temp[i][1]] for i in xrange(len(temp))]
		nfils_im = len(filaments_im)
		nimgs = sum([len(filaments_im[i]) for i in xrange(nfils_im)])
		recvcount.append(nimgs)
		
		if myid == main_node:
			for ii in xrange(nfils_im):
				numsegs = len(filaments_im[ii])
				for jj in xrange(numsegs):
					maptoimgID[int(filaments_im[ii][jj])] = counter
					counter += 1	

	if sum(recvcount) != total_nima:
		ERROR("recvcount is not calculated correctly!", "ehelix_MPI", 1, myid)
	
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
	print  " READ IMAGES ", myid,nima,nproc

	rise = int(dp/pixel_size)

	data_nx = data[0].get_xsize()
	data_ny = data[0].get_ysize()
	
	if data_nx != data_ny:
		ERROR('Input projections must be square.', "ehelix_streak_MPI", 1,myid)
	
	if(nx != data_ny):
		ERROR('Height of reference volume must be same as dimension of input projections', "ehelix_streak_MPI", 1,myid)
		
	data_nn = max(data_nx, data_ny)
	segmask = pad(model_blank(2*int(rmax), seg_ny, 1, 1.0), data_nx, data_ny, 1, 0.0)
	fdata = [None]*nima
	resetatone = False
	for im in xrange(nima):
		data[im].set_attr('ID', list_of_particles[im])
		st = Util.infomask(data[im], segmask, False)
		data[im] -= st[0]
		data[im] /= st[1]
		if CTF:
			qctf = data[im].get_attr("ctf_applied")
			if qctf == 0:
				ctf_params = data[im].get_attr("ctf")
				data[im] = filt_ctf(data[im], ctf_params)
				data[im].set_attr('ctf_applied', 1)
			elif qctf != 1:
				ERROR('Incorrectly set qctf flag', "ehelix_MPI", 1,myid)
		#  if FindPsi,  apply the angle to data[im], do fft and put in fdata[im]
		if FindPsi:
			phi,theta,psi,tsx,tsy = get_params_proj(data[im])
			if( theta != 0.0):
				if(abs(psi - 90.) < abs(psi - 270.0)):  gamma =  90.0
				else:                                   gamma = 270.0
				fdata[im] = fft( segmask*rot_shift2D(data[im], gamma-psi) )
			else:
				set_params_proj(data[im], [0.0, 90.0, 90.0, 0.0,0.0])
				fdata[im] = fft( segmask*data[im] )
		else:
			set_params_proj(data[im], [0.0, 90.0, 90.0, 0.0,0.0])
			fdata[im] = fft( segmask*data[im] )
		'''
		# check previous max and if does not exist set it to -1.e23
		p = data[im].get_attr_default('previousmax',-1.0e23)
		if( p == -1.0e23 ):
			resetatone = True
			data[im].set_attr('previousmax', p)
		'''
	del list_of_particles

	if debug:
		finfo.write( '%d loaded  \n' % nima )
		finfo.flush()
	if myid == main_node:
		# initialize data for the reference preparation function
		# No centering for helical reconstruction
		ref_data = [None, mask3D, None, None, None ]

	phiwobble = int(float(ywobble)/rise*dphi/delta+0.5)  # phiwobble is NOT in degrees, it is in nphi units

	nwx = 2*search_rng+3
	nwy = rise+2*ywobble+2
	nwxc = nwx//2
	nwyc = nwy//2
	nphi = int(360.0/delta + 0.5)
	#print  "  params  ",nwx,nwy,nwxc,nwyc,nphi
	if FindPsi:
		mode = "F"
		cnx = data_nn//2+1
		cny = cnx
		numr = Numrinit(1, data_nn//2-2, 1, mode)
		wr   = ringwe(numr, mode)
		maxrin = numr[len(numr)-1]
		crefim = [None]*nphi
	else:
		#  have to initialize them, otherwise problem with passing the arguments
		mode = "F"
		cnx = data_nx//2+1
		cny = cnx
		numr = []
		wr   = []
		maxrin = 0
		crefim = []

	terminate = 0
	Iter = 0

	pixer = [0.0]*nima
	
 	while Iter < max_iter:
		Iter += 1
		if myid == main_node:
			start_time = time()
			print_msg("\nITERATION #%3d\n"%(Iter))

		volft, kbz = prep_vol( vol )
		del vol

		refproj = [None]*nphi
		if( not Dsym):  rotproj = [None]*nphi
		else:           rotproj = []

		for iphi in xrange(nphi):
			refproj[iphi] = prgs( volft, kbz, [delta*iphi, 90.0, 90.0, 0.0, 0.0])
			st = Util.infomask(refproj[iphi] , segmask, True)
			refproj[iphi] -= st[0]
			refproj[iphi] /= st[1]
			refproj[iphi] = Util.muln_img(refproj[iphi], segmask )

			if FindPsi:
				temp = Util.Polar2Dm(refproj[iphi], cnx, cny, numr, mode)
				Util.Frngs(temp, numr)
				Util.Applyws(temp, numr, wr)
				crefim[iphi] = temp
			#  rotated in-plane by 180 are equivalent to rot_shift3D(vol,-90,180.0,90) with phi running as phi
			if(not Dsym):  rotproj[iphi] = fft( segmask * (rot_shift2D(refproj[iphi],180.0)) )
			refproj[iphi] = fft( segmask*(refproj[iphi]) )
		del volft
		#exit()
		#if myid == main_node:  
		astart_time = time()
		terminate = 0
		for ifil in xrange(nfils):
			if myid == main_node:  start_time = time()
			if myid == main_node:
				start_time = time()
			t1 = []
			for im in xrange(indcs[ifil][0], indcs[ifil][1]):
				t1.append(data[im].get_attr("xform.projection"))
			ldata = [data[im] for im in xrange(indcs[ifil][0],indcs[ifil][1])]
			#for im in xrange(len(ldata)):  ldata[im].set_attr("bestang", 10000.0)
			Util.constrained_helix_exhaustive(ldata, fdata[indcs[ifil][0]:indcs[ifil][1]], refproj, rotproj, [float(dp), float(dphi), float(rise), float(delta)], [int(nphi), int(phiwobble), int(rng), int(ywobble), int(Dsym), int(nwx), int(nwy), int(nwxc), int(nwyc)], FindPsi, float(psi_max), crefim, numr, int(maxrin), mode, int(cnx), int(cny))

			'''
			if doExhaustive:
				Util.constrained_helix_exhaustive(ldata, fdata[indcs[ifil][0]:indcs[ifil][1]], refproj, rotproj, [float(dp), float(dphi), float(rise), float(delta)], [int(nphi), int(phiwobble), int(rng), int(ywobble), int(Dsym), int(nwx), int(nwy), int(nwxc), int(nwyc)], FindPsi, float(psi_max), crefim, numr, int(maxrin), mode, int(cnx), int(cny))
				terminate = 0
			else:
				tempch = Util.constrained_helix(ldata, fdata[indcs[ifil][0]:indcs[ifil][1]], refproj, rotproj, [float(dp), float(dphi), float(rise), float(delta)], [int(nphi), int(phiwobble), int(rng), int(ywobble), int(Dsym), int(nwx), int(nwy), int(nwxc), int(nwyc)], FindPsi, float(psi_max), crefim, numr, int(maxrin), mode, int(cnx), int(cny))
				#print "tempch, Iter, myid: ", tempch, Iter, myid
				#tempch = constrained_helix_SHC(ldata, fdata[indcs[ifil][0]:indcs[ifil][1]], refproj, rotproj,  dp, dphi, rise, delta ,  nphi, phiwobble, rng, ywobble, Dsym, nwx, nwy, nwxc, nwyc , FindPsi, psi_max, crefim, numr, maxrin, mode, cnx, cny, myid, main_node)
				if tempch > -1:
					#if myid == main_node:
					#	print_msg("tempch %d\n"%tempch)
					terminate = 0
			'''
			for im in xrange(indcs[ifil][0], indcs[ifil][1]):
				temp = Util.get_transform_params(ldata[im-indcs[ifil][0]], "xform.projection", "spider")
				set_params_proj(data[im],[temp["phi"],temp["theta"],temp["psi"],-temp["tx"],-temp["ty"]])
				t2 = data[im].get_attr("xform.projection")
				pixer[im]  = max_3D_pixel_error(t1[im-indcs[ifil][0]], t2, data_nn//2-2)
				#print data[im].get_attr('ID'),pixer[im], numr[-3]
				#if not(doExhaustive):
				#	if Iter == 1 and resetatone:  data[im].set_attr('previousmax',-1.0e23)

			if FindPsi:
				for im in xrange(indcs[ifil][0], indcs[ifil][1]):
					fdata[im] = fft( segmask*rot_shift2D(data[im], ldata[im-indcs[ifil][0]].get_attr("bestang") ) )
					#bestang = ldata[im-indcs[ifil][0]].get_attr("bestang")
					#if( bestang < 10000.0): fdata[im] = fft( segmask*rot_shift2D(data[im], bestang ) )
				
			#print  "Parameters computed for filament",myid,ifil,time()-start_time;start_time = time()
			if myid == main_node:
				start_time = time()
		del ldata
		del refproj
		if(not Dsym):  del rotproj
		#print  "Time of alignment = ",myid,time()-astart_time
		mpi_barrier(MPI_COMM_WORLD)
		
		#output pixel errors
		from mpi import mpi_gatherv
		recvbuf = mpi_gatherv(pixer, nima, MPI_FLOAT, recvcount, disps, MPI_FLOAT, main_node, MPI_COMM_WORLD)
		mpi_barrier(MPI_COMM_WORLD)
		terminate = 0
		if myid == main_node:
			recvbuf = map(float, recvbuf)
			recvbuford = []
			for ii in xrange(total_nima):
				recvbuford.append( recvbuf[maptoimgID[ii]])
			write_text_file([range(len(recvbuford)), recvbuford], os.path.join(outdir, "pixer_%04d.txt"%(Iter)) )
			from statistics import hist_list
			lhist = 20
			region, histo = hist_list(recvbuford, lhist)
			if region[0] < 0.0:  region[0] = 0.0
			msg = "      Histogram of pixel errors\n      ERROR       number of particles\n"
			print_msg(msg)
			for lhx in xrange(lhist):
				msg = " %10.3f     %7d\n"%(region[lhx], histo[lhx])
				print_msg(msg)
			# Terminate if (100-termprec)% within 1 pixel error
			im_same = 0
			for lhx in xrange(lhist):
				if region[lhx] > 1.0: break
				im_same += histo[lhx]
			precn = 100*float(total_nima-im_same)/float(total_nima)
			msg = " Number of particles that changed orientations %7d, percentage of total: %5.1f\n"%(total_nima-im_same, precn)
			print_msg(msg)
			if precn <= termprec:  terminate = 1
			#print "main_node iter, terminate: ", Iter, terminate
			del region, histo, recvbuford
		del recvbuf
		terminate = mpi_bcast(terminate, 1, MPI_INT, 0, MPI_COMM_WORLD)
		terminate = int(terminate[0])
		#print "my id, Iter, terminate: ", myid, Iter, terminate
		#if myid == main_node:
		#	print_msg("Time of alignment = %\n"%(time()-astart_time));start_time = time()
		if CTF:  vol = recons3d_4nn_ctf_MPI(myid, data, symmetry=sym, snr = snr, npad = npad)
		else:    vol = recons3d_4nn_MPI(myid, data, symmetry=sym, npad = npad)
		if myid == main_node:
			print_msg("3D reconstruction time = %d\n"%(time()-start_time));start_time = time()

		if myid == main_node:
			#vol.write_image(os.path.join(outdir, "vol%03d.hdf"%Iter))
			vol = vol.helicise(pixel_size, dp, dphi, fract, rmax, rmin)
			vol = sym_vol(vol, symmetry=sym)
			ref_data[0] = vol
			vol = user_func(ref_data)
			vol = vol.helicise(pixel_size, dp, dphi, fract, rmax, rmin)
			vol = sym_vol(vol, symmetry=sym)
			vol.write_image(os.path.join(outdir, "volf%03d.hdf"%Iter))
			#if(Iter == max_iter-1):  drop_image(vol, os.path.join(outdir, "volfshift.hdf"))

		bcast_EMData_to_all(vol, myid, main_node)
		# write out headers, under MPI writing has to be done sequentially
		mpi_barrier(MPI_COMM_WORLD)
		par_str = ['xform.projection', 'ID']   #, 'previousmax']
		if myid == main_node:
			start_time = time()
			if(file_type(stack) == "bdb"):
				from utilities import recv_attr_dict_bdb
				recv_attr_dict_bdb(main_node, stack, data, par_str, 0, nima, nproc)
			else:
				from utilities import recv_attr_dict
				recv_attr_dict(main_node, stack, data, par_str, 0, nima, nproc)
			print_msg("Time to write header information= %d\n"%(time()-start_time))
			start_time = time()
		else:		send_attr_dict(main_node, data, par_str, 0, nima)
		if myid == main_node:
			# write params to text file
			header(stack, params='xform.projection', fexport=os.path.join(outdir, "parameters%04d.txt"%Iter))
			#header(stack, params='previousmax', fexport=os.path.join(outdir, "previousmax%04d.txt"%Iter))
		mpi_barrier(MPI_COMM_WORLD)
		
		if terminate > 0:
			print_end_msg("ehelix_MPI")
			return



def ali3d_saturn(stack, ref_vol, outdir, maskfile = None, ir = 1, ou = -1, rs = 1, 
            xr = "4 2 2 1", yr = "-1", ts = "1 1 0.5 0.25", delta = "10 6 4 4", an = "-1", apsi = "-1", deltapsi = "-1", startpsi = "-1",
	    center = -1, maxit = 5, CTF = False, snr = 1.0,  ref_a = "S", sym = "c1",
	    user_func_name = "ref_ali3d", fourvar = True, npad = 4, debug = False, MPI = False, termprec = 0.0, gamma=0.0):
	"""
		Name
			ali3d - Perform 3-D projection matching given initial reference volume and image series
		Input
			stack: set of 2-D images in a stack file, images have to be squares
			ref_vol: initial reference volume
			outdir: directory name into which the results will be written
			maskfile: filename of the file containing 3D mask.
			ir: inner radius for rotational correlation > 0 
			ou: outer radius for rotational correlation <int(nx/2)-1 
			rs: steps between rings in rotational correlation >0
			xr: range for translation search in x direction in each iteration, search is +/xr
			yr: range for translation search in y direction in each iteration, search is +/yr
			ts: step size of the translation search in both directions, search is -xr, -xr+ts, 0, xr-ts, xr, can be fractional.
			delta: angular step for the reference projections in respective iterations
			an: angular neighborhood for local searches
			center: average center method
			max_iter: maximum iterations at each angle step
			CTF: if the flag is present, program will use the CTF information stored in file headers
			snr: signal noise ratio used in the 3D reconstruction
			ref_a: method for creating quasi-uniform distribution of the projection directions of reference projections: "S" - spiral
			sym: symmetry of the refined structure
			function: name of the user-supplied-function
			MPI: if presetm use MPI version
		Output
			output_directory: directory name into which the output files will be written.
			header: the alignment parameters are stored in the headers of input files as Transform Object xform.proj
	"""
	if MPI:
		ali3d_saturnMPI(stack, ref_vol, outdir, maskfile, ir, ou, rs, xr, yr, ts, delta, an, apsi, deltapsi, startpsi, center, maxit, CTF
					, snr, ref_a, sym, user_func_name, fourvar, npad, debug, termprec, gamma)
		return

	from alignment      import proj_ali_incore, proj_ali_incore_local
	from utilities      import model_circle, drop_image, get_image, get_input_from_string
	from utilities      import get_params_proj
	from utilities      import estimate_3D_center, rotate_3D_shift
	from filter         import filt_params, fit_tanh, filt_tanl, filt_ctf
	from statistics     import fsc_mask
	import os
	import types
	from utilities      import print_begin_msg, print_end_msg, print_msg, angle_between_projections_directions
	from alignment      import Numrinit, prepare_refrings
	from projection     import prep_vol
	from math           import sqrt

	import user_functions
	user_func = user_functions.factory[user_func_name]

	if os.path.exists(outdir):  ERROR('Output directory exists, please change the name and restart the program', "ali3d", 1)
	os.mkdir(outdir)
	import global_def
	global_def.LOGFILE =  os.path.join(outdir, global_def.LOGFILE)
	print_begin_msg("ali3d_saturn")

	xrng        = get_input_from_string(xr)
	if  yr == "-1":  yrng = xrng
	else          :  yrng = get_input_from_string(yr)
	step        = get_input_from_string(ts)
	delta       = get_input_from_string(delta)
	lstp = min(len(xrng), len(yrng), len(step), len(delta))
	if an == "-1":
		an = [-1] * lstp
	else:
		an = get_input_from_string(an)
	first_ring  = int(ir)
	rstep       = int(rs)
	last_ring   = int(ou)
	max_iter    = int(maxit)
	center      = int(center)

	print_msg("Input stack                 : %s\n"%(stack))
	print_msg("Reference volume            : %s\n"%(ref_vol))	
	print_msg("Output directory            : %s\n"%(outdir))
	print_msg("Maskfile                    : %s\n"%(maskfile))
	print_msg("Inner radius                : %i\n"%(first_ring))

	vol     = EMData()
	vol.read_image(ref_vol)
	nx      = vol.get_xsize()
	if last_ring == -1:	last_ring = nx/2 - 2

	print_msg("Outer radius                : %i\n"%(last_ring))
	print_msg("Ring step                   : %i\n"%(rstep))
	print_msg("X search range              : %s\n"%(xrng))
	print_msg("Y search range              : %s\n"%(yrng))
	print_msg("Translational step          : %s\n"%(step))
	print_msg("Angular step                : %s\n"%(delta))
	print_msg("Angular search range        : %s\n"%(an))
	print_msg("Maximum iteration           : %i\n"%(max_iter))
	print_msg("Center type                 : %i\n"%(center))
	print_msg("CTF correction              : %s\n"%(CTF))
	print_msg("Signal-to-Noise Ratio       : %f\n"%(snr))
	print_msg("Reference projection method : %s\n"%(ref_a))
	print_msg("Symmetry group              : %s\n\n"%(sym))

	if maskfile :
		if type(maskfile) is types.StringType: mask3D = get_image(maskfile)
		else                                  : mask3D = maskfile
	else          :   mask3D = model_circle(last_ring, nx, nx, nx)
	mask2D = model_circle(last_ring, nx, nx) - model_circle(first_ring, nx, nx)
	numr   = Numrinit(first_ring, last_ring, rstep, "F")

	if CTF:
		from reconstruction import recons3d_4nn_ctf
		from filter         import filt_ctf
	else: from reconstruction import recons3d_4nn

	if debug:  outf = file(os.path.join(outdir, "progress"), "w")
	else:      outf = None

	# horatio active_refactoring Jy51i1EwmLD4tWZ9_00000_1
	# active = EMUtil.get_all_attributes(stack, 'active')
	# list_of_particles = []
	# for im in xrange(len(active)):
	# 	if(active[im]):  list_of_particles.append(im)
	# del active

	nima = EMUtil.get_image_count(stack)
	list_of_particles = range(nima)
	
	data = EMData.read_images(stack, list_of_particles)
	for im in xrange(len(data)):
		data[im].set_attr('ID', list_of_particles[im])
	if CTF:
		ctf_params = data[im].get_attr("ctf")
		st = Util.infomask(data[im], mask2D, False)
		data[im] -= st[0]
		data[im] = filt_ctf(data[im], ctf_params)
		data[im].set_attr('ctf_applied', 1)

	nima = len(data)
	# initialize data for the reference preparation function
	ref_data = [ mask3D, max(center,0), None, None ]#  for center -1 switch of centering by user function
	from random import shuffle
	cs = [0.0]*3
	# do the projection matching
	for N_step in xrange(lstp):
 		for Iter in xrange(max_iter):
			print_msg("\nITERATION #%3d\n"%(N_step*max_iter+Iter+1))

			volft, kb = prep_vol(vol)
			refrings = prepare_refrings( volft, kb, nx, delta[N_step], ref_a, sym, numr, MPI=False)
			del volft, kb
			rorder = range(nima)
			shuffle(rorder)
			pixel_error = [0.0]*nima
			for im in xrange(nima):
				if an[N_step] == -1:	
					peak, pixel_error[im] = proj_ali_incore(data[rorder[im]],refrings,numr,xrng[N_step],yrng[N_step],step[N_step])
					iref = data[rorder[im]].get_attr("referencenumber")
					temp_range = range(len(refrings))
					temp_range.reverse()
					phi   = refrings[iref].get_attr("phi")
					theta = refrings[iref].get_attr("theta")
					psi   = refrings[iref].get_attr("psi")
					paramA = [phi, theta, psi]
					for irr in temp_range:
						phi   = refrings[irr].get_attr("phi")
						theta = refrings[irr].get_attr("theta")
						psi   = refrings[irr].get_attr("psi")
						paramB = [phi, theta, psi]
						if angle_between_projections_directions(paramA, paramB) < gamma:
							del refrings[irr]
				else:
					peak, pixel_error[im] = proj_ali_incore_local(data[im],refrings,numr,xrng[N_step],yrng[N_step],step[N_step],an[N_step])

			from statistics import hist_list
			lhist = 20
			region, histo = hist_list(pixel_error, lhist)
			if region[0] < 0.0:  region[0] = 0.0
			msg = "      Histogram of pixel errors\n      ERROR       number of particles\n"
			print_msg(msg)
			for lhx in xrange(lhist):
				msg = " %10.3f     %7d\n"%(region[lhx], histo[lhx])
				print_msg(msg)
			# Terminate if 95% within 1 pixel error
			im = 0
			for lhx in xrange(lhist):
				if region[lhx] > 1.0: break
				im += histo[lhx]
			precn = 100*float(nima-im)/float(nima)
			msg = " Number of particles that changed orientations %7d, percentage of total: %5.1f\n"%(nima-im, precn)
			print_msg(msg)
			#if precn <= termprec:  terminate = 1   from shc

			if center == -1 and sym[0] == 'c':
				cs[0], cs[1], cs[2], dummy, dummy = estimate_3D_center(data)
				msg = "Average center x = %10.3f        Center y = %10.3f        Center z = %10.3f\n"%(cs[0], cs[1], cs[2])
				print_msg(msg)
				if int(sym[1]) > 1:
					cs[0] = cs[1] = 0.0
					print_msg("For symmetry group cn (n>1), we only center the volume in z-direction\n")
				rotate_3D_shift(data, [-cs[0], -cs[1], -cs[2]])
			"""
			if CTF:    vol1 = recons3d_4nn_ctf(data, range(0, nima, 2), snr, 1, sym)
			else:	   vol1 = recons3d_4nn(data, range(0, nima, 2), sym)
			if CTF:    vol2 = recons3d_4nn_ctf(data, range(1, nima, 2), snr, 1, sym)
			else:	   vol2 = recons3d_4nn(data, range(1, nima, 2), sym)

			fscc = fsc_mask(vol1, vol2, mask3D, 1.0, os.path.join(outdir, "resolution%04d"%(N_step*max_iter+Iter+1)))
			del vol1
			del vol2
			"""
			# calculate new and improved 3D
			vol_previous = vol
			if CTF:  vol = recons3d_4nn_ctf(data, range(nima), snr, 1, sym)
			else:	 vol = recons3d_4nn(data, range(nima), sym)
			# store the reference volume
			drop_image(vol, os.path.join(outdir, "vol%04d.hdf"%(N_step*max_iter+Iter+1)))
			ref_data[2] = vol
			ref_data[3] = None#fscc
			current_result = vol
			#  call user-supplied function to prepare reference image, i.e., center and filter it
			vol, dummy = user_func(ref_data)
			print_msg("L2 between this and previous volume: " + str(sqrt(vol.cmp("SqEuclidean",vol_previous,{"mask":mask3D,"zeromask":0,"normto":0}))) + "\n")
			print_msg("Dot product of the volume: " + str(vol.cmp("dot", vol, {"negative":0, "mask":mask3D})) + "\n")

			drop_image(vol, os.path.join(outdir, "volf%04d.hdf"%(N_step*max_iter+Iter+1)))
			paro = [None]*nima
			for im in xrange(nima):
				a1,a2,a3,a4,a5 = get_params_proj(data[im])
				paro[im] = [a1,a2,a3,a4,a5]
			from utilities import write_text_row
			write_text_row(paro,os.path.join(outdir, "params%04d.txt"%(N_step*max_iter+Iter+1)))
			del paro
			"""
			#  here we write header info
			from utilities import write_headers
			#from utilities import write_select_headers
			if CTF:
				for dat in data:  dat.set_attr('ctf_applied',0)
			for dat in data:  dat.del_attr('referencenumber')
			write_headers(stack, data, list_of_particles)
			#list_params= ['ID','xform.projection']
			#write_select_headers(stack, data, list_params)
			if CTF:
				for dat in data:  dat.set_attr('ctf_applied', 1)
			"""
	print_end_msg("ali3d")
	return current_result

def ali3d_MPI_chunks(stack, ref_vol, outdir, maskfile = None, ir = 1, ou = -1, rs = 1, 
            xr = "4 2 2 1", yr = "-1", ts = "1 1 0.5 0.25", delta = "10 6 4 4", an = "-1", apsi = "-1", deltapsi = "-1", startpsi = "-1",
	    center = -1, maxit = 5, CTF = False, snr = 1.0,  ref_a = "S", sym = "c1",  user_func_name = "ref_ali3d",
	    fourvar = True, npad = 4, debug = False, termprec = 0.0):

	from alignment       import Numrinit, prepare_refrings, proj_ali_incore, proj_ali_incore_local, proj_ali_incore_local_psi
	from utilities       import model_circle, get_image, drop_image, get_input_from_string
	from utilities       import bcast_list_to_all, bcast_number_to_all, reduce_EMData_to_root, bcast_EMData_to_all
	from utilities       import send_attr_dict
	from utilities       import get_params_proj, file_type
	from fundamentals    import rot_avg_image
	import os
	import types
	from utilities       import print_begin_msg, print_end_msg, print_msg
	from mpi             import mpi_bcast, mpi_comm_size, mpi_comm_rank, MPI_FLOAT, MPI_COMM_WORLD, mpi_barrier
	from mpi             import mpi_reduce, MPI_INT, MPI_SUM
	from filter          import filt_ctf
	from projection      import prep_vol, prgs
	from statistics      import hist_list, varf3d_MPI
	from applications    import MPI_start_end


	number_of_proc = mpi_comm_size(MPI_COMM_WORLD)
	myid           = mpi_comm_rank(MPI_COMM_WORLD)
	main_node = 0
	
	if os.path.exists(outdir):  ERROR('Output directory exists, please change the name and restart the program', "ali3d_MPI", 1, myid)
	mpi_barrier(MPI_COMM_WORLD)

	if myid == main_node:
		os.mkdir(outdir)
		import global_def
		global_def.LOGFILE =  os.path.join(outdir, global_def.LOGFILE)
		print_begin_msg("ali3d_MPI")
	mpi_barrier(MPI_COMM_WORLD)

	if debug:
		from time import sleep
		while not os.path.exists(outdir):
			print  "Node ",myid,"  waiting..."
			sleep(5)

		info_file = os.path.join(outdir, "progress%04d"%myid)
		finfo = open(info_file, 'w')
	else:
		finfo = None

	xrng        = get_input_from_string(xr)
	if  yr == "-1":  yrng = xrng
	else          :  yrng = get_input_from_string(yr)
	step        = get_input_from_string(ts)
	delta       = get_input_from_string(delta)
	lstp = min(len(xrng), len(yrng), len(step), len(delta))
	if an == "-1":
		an = [-1] * lstp
	else:
		an = get_input_from_string(an)

	if apsi == "-1":
		apsi = [-1] * lstp
	else:
		apsi = get_input_from_string(apsi)

	if deltapsi == "-1":
		deltapsi = [-1] * lstp
	else:
		deltapsi = get_input_from_string(deltapsi)

	if startpsi == "-1":
		startpsi = [-1] * lstp
	else:
		startpsi = get_input_from_string(startpsi)

	first_ring  = int(ir)
	rstep       = int(rs)
	last_ring   = int(ou)
	max_iter    = int(maxit)
	center      = int(center)

	vol     = EMData()
	vol.read_image(ref_vol)
	nx      = vol.get_xsize()
	if last_ring < 0:	last_ring = int(nx/2) - 2

	if myid == main_node:
		import user_functions
		user_func = user_functions.factory[user_func_name]

		print_msg("Input stack                 : %s\n"%(stack))
		print_msg("Reference volume            : %s\n"%(ref_vol))	
		print_msg("Output directory            : %s\n"%(outdir))
		print_msg("Maskfile                    : %s\n"%(maskfile))
		print_msg("Inner radius                : %i\n"%(first_ring))
		print_msg("Outer radius                : %i\n"%(last_ring))
		print_msg("Ring step                   : %i\n"%(rstep))
		print_msg("X search range              : %s\n"%(xrng))
		print_msg("Y search range              : %s\n"%(yrng))
		print_msg("Translational step          : %s\n"%(step))
		print_msg("Angular step                : %s\n"%(delta))
		print_msg("Angular search range (phi and theta)       : %s\n"%(an))
		print_msg("Angular search range (psi)                 : %s\n"%(apsi))
		print_msg("Delta psi                   : %s\n"%(deltapsi))
		print_msg("Start psi                   : %s\n"%(startpsi))
		print_msg("Maximum iteration           : %i\n"%(max_iter))
		print_msg("Percentage of change for termination: %f\n"%(termprec))
		print_msg("Center type                 : %i\n"%(center))
		print_msg("CTF correction              : %s\n"%(CTF))
		print_msg("Signal-to-Noise Ratio       : %f\n"%(snr))
		print_msg("Reference projection method : %s\n"%(ref_a))
		print_msg("Symmetry group              : %s\n\n"%(sym))
		print_msg("User function               : %s\n"%(user_func_name))

	if maskfile:
		if type(maskfile) is types.StringType: mask3D = get_image(maskfile)
		else:                                  mask3D = maskfile
	else: mask3D = model_circle(last_ring, nx, nx, nx)

	numr	= Numrinit(first_ring, last_ring, rstep, "F")
	mask2D  = model_circle(last_ring,nx,nx) - model_circle(first_ring,nx,nx)

	fscmask = mask3D  #model_circle(last_ring,nx,nx,nx)  For a fancy mask circle would work better  PAP 7/21/11
	if CTF:
		from reconstruction import rec3D_MPI
		from filter         import filt_ctf
	else:	 from reconstruction import rec3D_MPI_noCTF

	if myid == main_node:
		if file_type(stack) == "bdb":
			from EMAN2db import db_open_dict
			dummy = db_open_dict(stack, True)
		# horatio active_refactoring Jy51i1EwmLD4tWZ9_00000_1
		# active = EMUtil.get_all_attributes(stack, 'active')
		# list_of_particles = []
		# for im in xrange(len(active)):
		# 	if active[im]:  list_of_particles.append(im)
		# del active
		# nima = len(list_of_particles)

		nima = EMUtil.get_image_count(stack)
		list_of_particles = range(nima)
		
	else:
		nima = 0
	total_nima = bcast_number_to_all(nima, source_node = main_node)

	if myid != main_node:
		list_of_particles = [-1]*total_nima
	list_of_particles = bcast_list_to_all(list_of_particles, source_node = main_node)

	image_start, image_end = MPI_start_end(total_nima, number_of_proc, myid)
	# create a list of images for each node
	list_of_particles = list_of_particles[image_start: image_end]
	nima = len(list_of_particles)
	if debug:
		finfo.write("image_start, image_end: %d %d\n" %(image_start, image_end))
		finfo.flush()

	data = EMData.read_images(stack, list_of_particles)
	if fourvar:  original_data = []
	for im in xrange(nima):
		data[im].set_attr('ID', list_of_particles[im])
		if fourvar: original_data.append(data[im].copy())
		if CTF:
			ctf_params = data[im].get_attr("ctf")
			st = Util.infomask(data[im], mask2D, False)
			data[im] -= st[0]
			data[im] = filt_ctf(data[im], ctf_params)
			data[im].set_attr('ctf_applied', 1)

	if debug:
		finfo.write( '%d loaded  \n' % nima )
		finfo.flush()
	if myid == main_node:
		# initialize data for the reference preparation function
		ref_data = [ mask3D, max(center,0), None, None, None, None ]
		# for method -1, switch off centering in user function

	from time import time	

	#  this is needed for gathering of pixel errors
	disps = []
	recvcount = []
	for im in xrange(number_of_proc):
		if im == main_node :  disps.append(0)
		else:                 disps.append(disps[im-1] + recvcount[im-1])
		ib, ie = MPI_start_end(total_nima, number_of_proc, im)
		recvcount.append(ie-ib)

	pixer = [0.0]*nima
	peaks = [0.0]*nima
	cs = [0.0]*3
	total_iter = 0

# Consider the following subtractions before chunk calculations
#  1. input data
#  2. vol
#  3. possible 3d mask
#  4. prep_vol
#  5. Numrinit: points on rings
	mem_vol = vol.get_size()
	#  We know how much is the output volume from prep_vol, there is no need to do it experimentally.
	# How? TD
	volft, kb = prep_vol(vol)
	mem_prep_vol = volft.get_size()
	del volft, kb

	mem_mask3d = mask3D.get_size()	
	mem_inp = data[0].get_size() * nima

	mem_occupied = mem_vol + mem_prep_vol + mem_mask3d + mem_inp
	mem_refring = sum([numr[3*i+2] for i in range(len(numr/3))])
	
	gb_to_byte = float(1024**3)
	float_bytes = sys.getsizeof(float())
	#  It appears it confused words with bytes.  PAP
	#  Please put he explanation what the equation is.  Where 1024**3 is coming from??  
	#  In addition, it is unclear whether the equation is meant to be int or fload, now both re confused.
	# Revized TD
	schunk = int((global_def.MEM_PER_CPU_GB*gb_to_byte - mem_occupied*float_bytes) / mem_refring*float_bytes)


# 	phiEqpsi is an argument of prepare_refrings(). In ali3d_MPI when prepare_refrings() is called no
#   no argument for phiEqpsi is given and the default for it is "Minus"
#   Since even_angles() which was moved out of prepare_refrings(), needs this variable,
#   it is assigned here before even_angles() is called TD
	phiEqpsi = "Minus"  #  Where this is coming from ?? PAP
	# do the projection matching
	for N_step in xrange(lstp):
		terminate = 0
		Iter = -1

		if initial_theta is None:
			ref_angles = even_angles(delta[N_step], symmetry=sym, method = ref_a, phiEqpsi = phiEqpsi)
		else:
			if delta_theta is None: delta_theta = 1.0
			ref_angles = even_angles(delta[N_step], theta1 = initial_theta, theta2 = delta_theta, symmetry=sym, method = ref_a, phiEqpsi = phiEqpsi)
		nrefs = len(refrings)
		nchunk = nrefs/schunk + 1

		while Iter < max_iter-1 and terminate == 0:
			Iter += 1
			total_iter += 1
			if myid == main_node:
				start_time = time()
				print_msg("\nITERATION #%3d,  inner iteration #%3d\nDelta = %4.1f, an = %5.2f, xrange = %5.2f, yrange = %5.2f, step = %5.2f, delta psi = %5.2f, start psi = %5.2f\n"%(total_iter, Iter, delta[N_step], an[N_step], xrng[N_step],yrng[N_step],step[N_step],deltapsi[N_step],startpsi[N_step]))

			volft, kb = prep_vol(vol)
			for ch in range(nchunk):
				#  Generally, you need some kind of comments to clarify what you're doing  PAP
				# index ranges for the chunks with balanced load distribution TD
				ibeg = int(round(floor(nrefs)/nchunk*ch))
				iend = int(round(floor(nrefs)/nchunk*(ch+1)))
				#  Why do you need sym ??  PAP
				# sym is an argument passed to ali_3d_MPI() and prepare_refrings()
				#  If not specified the default will be used in prepare_refrings() 
				#  and it wouldn't matter what the caller of ali3d_MPI() passed in sym TD
				refrings = prepare_refrings_chunks(volft, kb, nx, ref_angles[ibeg:iend], ref_a, sym, numr, True, ant = max(an[N_step],0.0)*1.1)  # 1.1 is to have extra safety

				if myid == main_node:
					print_msg("Time to prepare rings: %d\n" % (time()-start_time))
					start_time = time()
	
				for im in xrange(nima):
					if deltapsi[N_step] > 0.0:
						from alignment import proj_ali_incore_delta
						peak, pxr = proj_ali_incore_delta(data[im],refrings,numr,xrng[N_step],yrng[N_step],step[N_step],startpsi[N_step],deltapsi[N_step],finfo)						
					elif an[N_step] == -1:
						peak, pxr = proj_ali_incore(data[im],refrings,numr,xrng[N_step],yrng[N_step],step[N_step],finfo)
					else:
						if apsi[N_step] == -1:
							peak, pxr = proj_ali_incore_local(data[im],refrings,numr,xrng[N_step],yrng[N_step],step[N_step],an[N_step],finfo, sym = sym)
						else:
							peak, pxr = proj_ali_incore_local_psi(data[im],refrings,numr,xrng[N_step],yrng[N_step],step[N_step],an[N_step],apsi[N_step],finfo)
					if( peak > peaks[im]):
						peaks[im] = peak
						pixer[im] = pxr
						data[im].set_attr("previousmax", peak)
				

			del volft, kb
			if myid == main_node:
				print_msg("Time of alignment = %d\n"%(time()-start_time))
				start_time = time()

			#output pixel errors
			from mpi import mpi_gatherv
			recvbuf = mpi_gatherv(pixer, nima, MPI_FLOAT, recvcount, disps, MPI_FLOAT, main_node, MPI_COMM_WORLD)
			mpi_barrier(MPI_COMM_WORLD)
			terminate = 0
			if myid == main_node:
				recvbuf = map(float, recvbuf)
				from statistics import hist_list
				lhist = 20
				region, histo = hist_list(recvbuf, lhist)
				if region[0] < 0.0:  region[0] = 0.0
				msg = "      Histogram of pixel errors\n      ERROR       number of particles\n"
				print_msg(msg)
				for lhx in xrange(lhist):
					msg = " %10.3f     %7d\n"%(region[lhx], histo[lhx])
					print_msg(msg)
				# Terminate if 95% within 1 pixel error
				im = 0
				for lhx in xrange(lhist):
					if region[lhx] > 1.0: break
					im += histo[lhx]
				precn = 100*float(total_nima-im)/float(total_nima)
				msg = " Number of particles that changed orientations %7d, percentage of total: %5.1f\n"%(total_nima-im, precn)
				print_msg(msg)
				if precn <= termprec:  terminate = 1
				del region, histo
			del recvbuf
			terminate = mpi_bcast(terminate, 1, MPI_INT, 0, MPI_COMM_WORLD)
			terminate = int(terminate[0])

			if center == -1 and sym[0] == 'c':
				from utilities      import estimate_3D_center_MPI, rotate_3D_shift
				cs[0], cs[1], cs[2], dummy, dummy = estimate_3D_center_MPI(data, total_nima, myid, number_of_proc, main_node)
				if myid == main_node:
					msg = " Average center x = %10.3f        Center y = %10.3f        Center z = %10.3f\n"%(cs[0], cs[1], cs[2])
					print_msg(msg)
				if int(sym[1]) > 1:
					cs[0] = cs[1] = 0.0
					if myid == main_node:
						print_msg("For symmetry group cn (n>1), we only center the volume in z-direction\n")
				cs = mpi_bcast(cs, 3, MPI_FLOAT, main_node, MPI_COMM_WORLD)
				cs = [-float(cs[0]), -float(cs[1]), -float(cs[2])]
				rotate_3D_shift(data, cs)

			# write out headers, under MPI writing has to be done sequentially
			mpi_barrier(MPI_COMM_WORLD)
			par_str = ['xform.projection', 'previousmax', 'ID']
			if myid == main_node:
	   			if(file_type(stack) == "bdb"):
	        			from utilities import recv_attr_dict_bdb
	        			recv_attr_dict_bdb(main_node, stack, data, par_str, image_start, image_end, number_of_proc)
	        		else:
	        			from utilities import recv_attr_dict
	        			recv_attr_dict(main_node, stack, data, par_str, image_start, image_end, number_of_proc)
				print_msg("Time to write header information= %d\n"%(time()-start_time))
				start_time = time()
	        	else:	       send_attr_dict(main_node, data, par_str, image_start, image_end)

			if CTF: vol, fscc = rec3D_MPI(data, snr, sym, fscmask, os.path.join(outdir, "resolution%04d"%(total_iter)), myid, main_node, npad = npad)
			else:   vol, fscc = rec3D_MPI_noCTF(data, sym, fscmask, os.path.join(outdir, "resolution%04d"%(total_iter)), myid, main_node, npad = npad)

			if myid == main_node:
				print_msg("3D reconstruction time = %d\n"%(time()-start_time))
				start_time = time()

			if fourvar:
			#  Compute Fourier variance
				for im in xrange(nima):
					original_data[im].set_attr( 'xform.projection', data[im].get_attr('xform.projection') )
				varf = varf3d_MPI(original_data, ssnr_text_file = os.path.join(outdir, "ssnr%04d"%(total_iter)), mask2D = None, reference_structure = vol, ou = last_ring, rw = 1.0, npad = 1, CTF = CTF, sign = 1, sym =sym, myid = myid)
				if myid == main_node:
					print_msg("Time to calculate 3D Fourier variance= %d\n"%(time()-start_time))
					start_time = time()
					varf = 1.0/varf
			else:  varf = None

			if myid == main_node:
				drop_image(vol, os.path.join(outdir, "vol%04d.hdf"%(total_iter)))
				ref_data[2] = vol
				ref_data[3] = fscc
				ref_data[4] = varf
				#  call user-supplied function to prepare reference image, i.e., center and filter it
				vol, cs = user_func(ref_data)
				drop_image(vol, os.path.join(outdir, "volf%04d.hdf"%(total_iter)))

			del varf
			bcast_EMData_to_all(vol, myid, main_node)


	if myid == main_node: print_end_msg("ali3d_MPI")

def prepare_refrings_chunks( volft, kb, nz, ref_angles, ref_a, sym = "c1", numr = None, MPI=False, \
						phiEqpsi = "Minus", kbx = None, kby = None, initial_theta = None, \
						delta_theta = None, ant = -1.0):
	"""
		Generate quasi-evenly distributed reference projections converted to rings
		ant - neighborhood for local searches.  I believe the fastest way to deal with it in case of point-group symmetry
		        it to generate cushion projections within an/2 of the border of unique zone
	"""
	from projection   import prep_vol, prgs
	from applications import MPI_start_end
	from utilities    import even_angles, getfvec
	from types import BooleanType

	# mpi communicator can be sent by the MPI parameter
	if type(MPI) is BooleanType:
		if MPI:
			from mpi import MPI_COMM_WORLD
			mpi_comm = MPI_COMM_WORLD
	else:
		mpi_comm = MPI
		MPI = True

	# generate list of Eulerian angles for reference projections
	#  phi, theta, psi
	mode = "F"
# 	if initial_theta is None:
# 		ref_angles = even_angles(delta, symmetry=sym, method = ref_a, phiEqpsi = phiEqpsi)
# 	else:
# 		if delta_theta is None: delta_theta = 1.0
# 		ref_angles = even_angles(delta, theta1 = initial_theta, theta2 = delta_theta, symmetry=sym, method = ref_a, phiEqpsi = phiEqpsi)
	wr_four  = ringwe(numr, mode)
	cnx = nz//2 + 1
	cny = nz//2 + 1
	num_ref = len(ref_angles)

	if MPI:
		from mpi import mpi_comm_rank, mpi_comm_size
		myid = mpi_comm_rank( mpi_comm )
		ncpu = mpi_comm_size( mpi_comm )
	else:
		ncpu = 1
		myid = 0
	
	ref_start, ref_end = MPI_start_end(num_ref, ncpu, myid)

	refrings = []     # list of (image objects) reference projections in Fourier representation

	sizex = numr[len(numr)-2] + numr[len(numr)-1]-1

	for i in xrange(num_ref):
		prjref = EMData()
		prjref.set_size(sizex, 1, 1)
		refrings.append(prjref)

	if kbx is None:
		for i in xrange(ref_start, ref_end):
			prjref = prgs(volft, kb, [ref_angles[i][0], ref_angles[i][1], ref_angles[i][2], 0.0, 0.0])
			cimage = Util.Polar2Dm(prjref, cnx, cny, numr, mode)  # currently set to quadratic....
			Util.Normalize_ring(cimage, numr)
			Util.Frngs(cimage, numr)
			Util.Applyws(cimage, numr, wr_four)
			refrings[i] = cimage
	else:
		for i in xrange(ref_start, ref_end):
			prjref = prgs(volft, kb, [ref_angles[i][0], ref_angles[i][1], ref_angles[i][2], 0.0, 0.0], kbx, kby)
			cimage = Util.Polar2Dm(prjref, cnx, cny, numr, mode)  # currently set to quadratic....
			Util.Normalize_ring(cimage, numr)
			Util.Frngs(cimage, numr)
			Util.Applyws(cimage, numr, wr_four)
			refrings[i] = cimage

	if MPI:
		from utilities import bcast_EMData_to_all
		for i in xrange(num_ref):
			for j in xrange(ncpu):
				ref_start, ref_end = MPI_start_end(num_ref, ncpu, j)
				if i >= ref_start and i < ref_end: rootid = j
			bcast_EMData_to_all(refrings[i], myid, rootid, comm=mpi_comm)

	for i in xrange(len(ref_angles)):
		n1,n2,n3 = getfvec(ref_angles[i][0], ref_angles[i][1])
		refrings[i].set_attr_dict( {"phi":ref_angles[i][0], "theta":ref_angles[i][1], "psi":ref_angles[i][2], "n1":n1, "n2":n2, "n3":n3} )

	return refrings




def ali3d_MPI(stack, ref_vol, outdir, maskfile = None, ir = 1, ou = -1, rs = 1,
            xr = "4 2 2 1", yr = "-1", ts = "1 1 0.5 0.25", delta = "10 6 4 4", an = "-1", apsi = "-1", deltapsi = "-1", startpsi = "-1",
	    center = -1, maxit = 5, CTF = False, snr = 1.0,  ref_a = "S", sym = "c1",  user_func_name = "ref_ali3d",
	    fourvar = True, npad = 2, debug = False, termprec = 0.0):

	from alignment       import Numrinit, prepare_refrings, proj_ali_incore, proj_ali_incore_local, proj_ali_incore_local_psi
	from utilities       import model_circle, get_image, drop_image, get_input_from_string
	from utilities       import bcast_list_to_all, bcast_number_to_all, reduce_EMData_to_root, bcast_EMData_to_all
	from utilities       import send_attr_dict
	from utilities       import get_params_proj, file_type
	from fundamentals    import rot_avg_image
	import os
	import types
	from utilities       import print_begin_msg, print_end_msg, print_msg
	from mpi             import mpi_bcast, mpi_comm_size, mpi_comm_rank, MPI_FLOAT, MPI_COMM_WORLD, mpi_barrier
	from mpi             import mpi_reduce, MPI_INT, MPI_SUM
	from filter          import filt_ctf
	from projection      import prep_vol, prgs
	from statistics      import hist_list, varf3d_MPI
	from applications    import MPI_start_end


	number_of_proc = mpi_comm_size(MPI_COMM_WORLD)
	myid           = mpi_comm_rank(MPI_COMM_WORLD)
	main_node = 0

	if os.path.exists(outdir):  ERROR('Output directory exists, please change the name and restart the program', "ali3d_MPI", 1, myid)
	mpi_barrier(MPI_COMM_WORLD)

	if myid == main_node:
		os.mkdir(outdir)
		import global_def
		global_def.LOGFILE =  os.path.join(outdir, global_def.LOGFILE)
		print_begin_msg("ali3d_MPI")
	mpi_barrier(MPI_COMM_WORLD)

	if debug:
		from time import sleep
		while not os.path.exists(outdir):
			print  "Node ",myid,"  waiting..."
			sleep(5)

		info_file = os.path.join(outdir, "progress%04d"%myid)
		finfo = open(info_file, 'w')
	else:
		finfo = None

	xrng        = get_input_from_string(xr)
	if  yr == "-1":  yrng = xrng
	else          :  yrng = get_input_from_string(yr)
	step        = get_input_from_string(ts)
	delta       = get_input_from_string(delta)
	lstp = min(len(xrng), len(yrng), len(step), len(delta))
	if an == "-1":
		an = [-1] * lstp
	else:
		an = get_input_from_string(an)

	if apsi == "-1":
		apsi = [-1] * lstp
	else:
		apsi = get_input_from_string(apsi)

	if deltapsi == "-1":
		deltapsi = [-1] * lstp
	else:
		deltapsi = get_input_from_string(deltapsi)

	if startpsi == "-1":
		startpsi = [-1] * lstp
	else:
		startpsi = get_input_from_string(startpsi)

	first_ring  = int(ir)
	rstep       = int(rs)
	last_ring   = int(ou)
	max_iter    = int(maxit)
	center      = int(center)

	vol     = EMData()
	vol.read_image(ref_vol)
	nx      = vol.get_xsize()
	if last_ring < 0:	last_ring = int(nx/2) - 2

	if myid == main_node:
		import user_functions
		user_func = user_functions.factory[user_func_name]

		print_msg("Input stack                 : %s\n"%(stack))
		print_msg("Reference volume            : %s\n"%(ref_vol))
		print_msg("Output directory            : %s\n"%(outdir))
		print_msg("Maskfile                    : %s\n"%(maskfile))
		print_msg("Inner radius                : %i\n"%(first_ring))
		print_msg("Outer radius                : %i\n"%(last_ring))
		print_msg("Ring step                   : %i\n"%(rstep))
		print_msg("X search range              : %s\n"%(xrng))
		print_msg("Y search range              : %s\n"%(yrng))
		print_msg("Translational step          : %s\n"%(step))
		print_msg("Angular step                : %s\n"%(delta))
		print_msg("Angular search range (phi and theta)       : %s\n"%(an))
		print_msg("Angular search range (psi)                 : %s\n"%(apsi))
		print_msg("Delta psi                   : %s\n"%(deltapsi))
		print_msg("Start psi                   : %s\n"%(startpsi))
		print_msg("Maximum iteration           : %i\n"%(max_iter))
		print_msg("Percentage of change for termination: %f\n"%(termprec))
		print_msg("Center type                 : %i\n"%(center))
		print_msg("CTF correction              : %s\n"%(CTF))
		print_msg("Signal-to-Noise Ratio       : %f\n"%(snr))
		print_msg("Reference projection method : %s\n"%(ref_a))
		print_msg("Symmetry group              : %s\n\n"%(sym))
		print_msg("User function               : %s\n"%(user_func_name))

	if maskfile:
		if type(maskfile) is types.StringType: mask3D = get_image(maskfile)
		else:                                  mask3D = maskfile
	else: mask3D = model_circle(last_ring, nx, nx, nx)

	numr	= Numrinit(first_ring, last_ring, rstep, "F")
	mask2D  = model_circle(last_ring,nx,nx) - model_circle(first_ring,nx,nx)

	fscmask = mask3D  #model_circle(last_ring,nx,nx,nx)  For a fancy mask circle would work better  PAP 7/21/11
	if CTF:
		from reconstruction import rec3D_MPI
		from filter         import filt_ctf
	else:	 from reconstruction import rec3D_MPI_noCTF

	if myid == main_node:
		if file_type(stack) == "bdb":
			from EMAN2db import db_open_dict
			dummy = db_open_dict(stack, True)
		# horatio active_refactoring Jy51i1EwmLD4tWZ9_00000_1
		# active = EMUtil.get_all_attributes(stack, 'active')
		# list_of_particles = []
		# for im in xrange(len(active)):
		# 	if active[im]:  list_of_particles.append(im)
		# del active
		# nima = len(list_of_particles)
		
		nima = EMUtil.get_image_count(stack)
		list_of_particles = range(nima)
	
	else:
		nima = 0
	total_nima = bcast_number_to_all(nima, source_node = main_node)

	if myid != main_node:
		list_of_particles = [-1]*total_nima
	list_of_particles = bcast_list_to_all(list_of_particles, source_node = main_node)

	image_start, image_end = MPI_start_end(total_nima, number_of_proc, myid)
	# create a list of images for each node
	list_of_particles = list_of_particles[image_start: image_end]
	nima = len(list_of_particles)
	if debug:
		finfo.write("image_start, image_end: %d %d\n" %(image_start, image_end))
		finfo.flush()

	data = EMData.read_images(stack, list_of_particles)
	if fourvar:  original_data = []
	for im in xrange(nima):
		data[im].set_attr('ID', list_of_particles[im])
		if fourvar: original_data.append(data[im].copy())
		if CTF:
			ctf_params = data[im].get_attr("ctf")
			st = Util.infomask(data[im], mask2D, False)
			data[im] -= st[0]
			data[im] = filt_ctf(data[im], ctf_params)
			data[im].set_attr('ctf_applied', 1)

	if debug:
		finfo.write( '%d loaded  \n' % nima )
		finfo.flush()
	if myid == main_node:
		# initialize data for the reference preparation function
		ref_data = [ mask3D, max(center,0), None, None, None, None ]
		# for method -1, switch off centering in user function

	from time import time

	#  this is needed for gathering of pixel errors
	disps = []
	recvcount = []
	for im in xrange(number_of_proc):
		if im == main_node :  disps.append(0)
		else:                 disps.append(disps[im-1] + recvcount[im-1])
		ib, ie = MPI_start_end(total_nima, number_of_proc, im)
		recvcount.append(ie-ib)

	pixer = [0.0]*nima
	cs = [0.0]*3
	total_iter = 0
	# do the projection matching
	for N_step in xrange(lstp):
		terminate = 0
		Iter = -1
 		while Iter < max_iter-1 and terminate == 0:
			Iter += 1
			total_iter += 1
			if myid == main_node:
				start_time = time()
				print_msg("\nITERATION #%3d,  inner iteration #%3d\nDelta = %4.1f, an = %5.2f, xrange = %5.2f, yrange = %5.2f, step = %5.2f, delta psi = %5.2f, start psi = %5.2f\n"%(total_iter, Iter, delta[N_step], an[N_step], xrng[N_step],yrng[N_step],step[N_step],deltapsi[N_step],startpsi[N_step]))


			volft, kb = prep_vol(vol)
			refrings = prepare_refrings(volft, kb, nx, delta[N_step], ref_a, sym, numr, True)


			del volft, kb
			if myid == main_node:
				print_msg("Time to prepare rings: %d\n" % (time()-start_time))
				start_time = time()

			for im in xrange(nima):
				if deltapsi[N_step] > 0.0:
					print "A1"
					from alignment import proj_ali_incore_delta
					peak, pixer[im] = proj_ali_incore_delta(data[im],refrings,numr,xrng[N_step],yrng[N_step],step[N_step],startpsi[N_step],deltapsi[N_step],finfo)
				elif an[N_step] == -1:
					print "A2"
					peak, pixer[im] = proj_ali_incore(data[im],refrings,numr,xrng[N_step],yrng[N_step],step[N_step],finfo)
				else:

					if apsi[N_step] == -1:
						print numr,xrng[N_step],yrng[N_step],step[N_step],an[N_step],finfo, sym
						peak, pixer[im] = proj_ali_incore_local(data[im],refrings,numr,xrng[N_step],yrng[N_step],step[N_step],an[N_step],finfo, sym = sym)
						print "A3"
					else:
						peak, pixer[im] = proj_ali_incore_local_psi(data[im],refrings,numr,xrng[N_step],yrng[N_step],step[N_step],an[N_step],apsi[N_step],finfo)
						print "A4"
				data[im].set_attr("previousmax", peak)


			import sys
			print "I'm here"
			sys.exit(0)


			if myid == main_node:
				print_msg("Time of alignment = %d\n"%(time()-start_time))
				start_time = time()

			#output pixel errors
			from mpi import mpi_gatherv
			recvbuf = mpi_gatherv(pixer, nima, MPI_FLOAT, recvcount, disps, MPI_FLOAT, main_node, MPI_COMM_WORLD)
			mpi_barrier(MPI_COMM_WORLD)
			terminate = 0
			if myid == main_node:
				recvbuf = map(float, recvbuf)
				from statistics import hist_list
				lhist = 20
				region, histo = hist_list(recvbuf, lhist)
				if region[0] < 0.0:  region[0] = 0.0
				msg = "      Histogram of pixel errors\n      ERROR       number of particles\n"
				print_msg(msg)
				for lhx in xrange(lhist):
					msg = " %10.3f     %7d\n"%(region[lhx], histo[lhx])
					print_msg(msg)
				# Terminate if 95% within 1 pixel error
				im = 0
				for lhx in xrange(lhist):
					if region[lhx] > 1.0: break
					im += histo[lhx]
				precn = 100*float(total_nima-im)/float(total_nima)
				msg = " Number of particles that changed orientations %7d, percentage of total: %5.1f\n"%(total_nima-im, precn)
				print_msg(msg)
				if precn <= termprec:  terminate = 1
				del region, histo
			del recvbuf
			terminate = mpi_bcast(terminate, 1, MPI_INT, 0, MPI_COMM_WORLD)
			terminate = int(terminate[0])

			if center == -1 and sym[0] == 'c':
				from utilities      import estimate_3D_center_MPI, rotate_3D_shift
				cs[0], cs[1], cs[2], dummy, dummy = estimate_3D_center_MPI(data, total_nima, myid, number_of_proc, main_node)
				if myid == main_node:
					msg = " Average center x = %10.3f        Center y = %10.3f        Center z = %10.3f\n"%(cs[0], cs[1], cs[2])
					print_msg(msg)
				if int(sym[1]) > 1:
					cs[0] = cs[1] = 0.0
					if myid == main_node:
						print_msg("For symmetry group cn (n>1), we only center the volume in z-direction\n")
				cs = mpi_bcast(cs, 3, MPI_FLOAT, main_node, MPI_COMM_WORLD)
				cs = [-float(cs[0]), -float(cs[1]), -float(cs[2])]
				rotate_3D_shift(data, cs)

			# write out headers, under MPI writing has to be done sequentially
			mpi_barrier(MPI_COMM_WORLD)
			par_str = ['xform.projection', 'previousmax', 'ID']
			if myid == main_node:
	   			if(file_type(stack) == "bdb"):
	        			from utilities import recv_attr_dict_bdb
	        			recv_attr_dict_bdb(main_node, stack, data, par_str, image_start, image_end, number_of_proc)
	        		else:
	        			from utilities import recv_attr_dict
	        			recv_attr_dict(main_node, stack, data, par_str, image_start, image_end, number_of_proc)
				print_msg("Time to write header information= %d\n"%(time()-start_time))
				start_time = time()
	        	else:	       send_attr_dict(main_node, data, par_str, image_start, image_end)

			if CTF: vol, fscc = rec3D_MPI(data, snr, sym, fscmask, os.path.join(outdir, "resolution%04d"%(total_iter)), myid, main_node, npad = npad)
			else:   vol, fscc = rec3D_MPI_noCTF(data, sym, fscmask, os.path.join(outdir, "resolution%04d"%(total_iter)), myid, main_node, npad = npad)

			if myid == main_node:
				print_msg("3D reconstruction time = %d\n"%(time()-start_time))
				start_time = time()

			if fourvar:
			#  Compute Fourier variance
				for im in xrange(nima):
					original_data[im].set_attr( 'xform.projection', data[im].get_attr('xform.projection') )
				varf = varf3d_MPI(original_data, ssnr_text_file = os.path.join(outdir, "ssnr%04d"%(total_iter)), mask2D = None, reference_structure = vol, ou = last_ring, rw = 1.0, npad = 1, CTF = CTF, sign = 1, sym =sym, myid = myid)
				if myid == main_node:
					print_msg("Time to calculate 3D Fourier variance= %d\n"%(time()-start_time))
					start_time = time()
					varf = 1.0/varf
			else:  varf = None

			if myid == main_node:
				drop_image(vol, os.path.join(outdir, "vol%04d.hdf"%(total_iter)))
				ref_data[2] = vol
				ref_data[3] = fscc
				ref_data[4] = varf
				#  call user-supplied function to prepare reference image, i.e., center and filter it
				vol, cs = user_func(ref_data)
				drop_image(vol, os.path.join(outdir, "volf%04d.hdf"%(total_iter)))

			del varf
			bcast_EMData_to_all(vol, myid, main_node)


	if myid == main_node: print_end_msg("ali3d_MPI")



def doCones_Horatio01():

	import os
	os.chdir("/Users/hvoicu/EMAN2/src/test_scripts")

	import sys

	redirect_to_stdout = True
	# redirect_to_stdout = False

	# How my volume was generated
	# os.system("sxpdb2em.py tteftu_with_tRNA.pdb tmp.hdf --apix=12.2 --box=64")
	# os.system("e2proc3d.py tmp.hdf my_volume.hdf --process=filter.lowpass.tanh:cutoff_abs=0.25:fall_off=0.35")


	from utilities import get_symt, even_angles

	from math import sin, cos, radians
	from sys import exit
	from morphology import  bracket_def
	from utilities import assign_projangles
	from projection import prep_vol


	def get_latest_output_file_increment_value(file_prefix):
		import os
		import glob
		list_of_files = glob.glob(file_prefix + "_????_*")

		count_values = []
		for fn in list_of_files:
			try:
				count_values.append(int(fn[(len(file_prefix)+1):(len(file_prefix)+1+4)]))
			except:
				pass
		if count_values == []:
			return 1
		else:
			return max(count_values) + 1


	delta = 4
	step, an = (0.5,20)
	numberofcones = 2
	cone_neighborhood = 5

	for delta in [4, 5]:
		for step in [0.5, 1]:
			for numberofcones in range(5,21,5):
				for cone_neighborhood in [1.5, 2, 3, 4, 5]:

					# def do_cones_horatio():


					ref_a="S"
					my_vol = EMData("my_volume.hdf")

					def computenumberofrefs(x, dat):
						#  dat = [sym, desired number of refs]
						return (len(even_angles(x,symmetry = dat[0])) - dat[1])**2


					# an = angular neighborhood
					xrng, yrng, = (0,0)

					volft,kb = prep_vol(my_vol)
					nx = my_vol.get_xsize()
					nz = my_vol.get_zsize()
					cnx = cny = nz//2 + 1

					numr = Numrinit(1,15)
					mode = "F"
					wr_four = ringwe(numr, mode)
					finfo = None



					projangles = even_angles(delta,theta2=180.0)

					# proj_data = prgs(volft, kb, projangles)
					proj_data = project3d(my_vol, listagls=projangles)


					import copy

					original_proj_data = []
					for im in proj_data:
						original_proj_data.append(copy.deepcopy(im))

					import random
					random.seed(123)


					# outside_an = [None]*len(proj_data)
					ran_phi_theta = [None]*len(proj_data)
					modification_parameters = [None]*len(proj_data)
					for idx, im in enumerate(proj_data):
						im.set_attr('previousmax', -1.0e23)
						im.set_attr("ID", idx)
						# horatio active_refactoring Jy51i1EwmLD4tWZ9_00000_1
						# im.set_attr("active", 1)
						random_phi = random.randint(1,an) - an/2
						random_theta = random.randint(1,an) - an/2
						ran_phi_theta[idx] = (random_phi, random_theta)
						# outside_an[idx] = int((abs(random_phi)>an) or abs(random_theta)>an)
						modification_parameters[idx]= [(projangles[idx][0]+random_phi)%360, (projangles[idx][1] + random_theta)%180, projangles[idx][2], 0, 0]
						set_params_proj(im, modification_parameters[idx])
						# set_params_proj(im, [projangles[idx][0], projangles[idx][1] + 10, projangles[idx][2], 0, 0], "xform.anchor")

					original_proj_data_modified = []
					for im in proj_data:
						original_proj_data_modified.append(copy.deepcopy(im))


					sym = 'c1'
					#  Number of reference image that fit into the memory.  Will have to be hardwired.



					if redirect_to_stdout:
						old_stdout = sys.stdout
						# sys.stdout = open('output_cones=%02d_delta=%.2f_step=%.2f_an=%.2f_cone_neighborhood_times_an=%.2f.txt'%(numberofcones,delta, step, an, cone_neighborhood), 'w')
						file_prefix = "output_regular"
						sys.stdout = open('%s_%04d_cones=%02d_delta=%.2f_step=%.2f_an=%.2f_cone_neighborhood_times_an=%.2f.txt'%(file_prefix, get_latest_output_file_increment_value(file_prefix), numberofcones,delta, step, an, cone_neighborhood), 'w')

						# sys.stdout = open('output_.txt', 'w')
						sys.stdout.flush()

					print " number of cones ",numberofcones,'   ',len(projangles)

					if( numberofcones == 1):
						print "  One cone, i.e., standard code"
						refrings = prepare_refrings(volft, kb, nx, delta, ref_a, sym, numr, MPI = False)
						for im in xrange(len(projangles)):
							peak, pixer = proj_ali_incore_local(proj_data[im], refrings, numr, xrng, yrng, step, an, finfo, sym = sym)
							proj_data[im].set_attr("previousmax", peak)

					else:
						rs = delta
						h = 1.0
						dat = [sym, numberofcones]
						def1, def2 = bracket_def(computenumberofrefs,dat, rs, h)
						def1, val  = goldsearch_astigmatism(computenumberofrefs, dat, def1, def2, tol=1.0)
						print '  computed delta ',def1
						print 'computed number of cones  ',  len(even_angles(def1,symmetry = sym))
						coneangles = even_angles(def1, symmetry = sym)
						#  assign exclusively angl to coneangles
						assignments = assign_projangles(projangles, coneangles)
						print 'assignments  ',len(assignments)
						# for k in xrange(len(assignments)):
						# 	# print k,projangles[assignments[k][0]],projangles[assignments[k][-1]]
						# 	for j in range(10):
						# 		print k,projangles[assignments[k][j]]
						# 	print
						# 	for j in reversed(range(10)):
						# 		print k,projangles[assignments[k][-j]]
						# 	print
						# 	print

							# projangles[assignments[k][-2]],projangles[assignments[k][-1]]
						#  for each coneangle compute refimages within 1.5*an degrees coneangle
						print "Finished 1"

						for k in xrange(len(coneangles)):
							if(len(assignments[k]) > 0):
								refsincone = even_angles(delta,symmetry = sym)
								ant = cone_neighborhood*an
								refsincone = cone_ang( refsincone, coneangles[k][0], coneangles[k][1], ant )
								refrings = refprojs( volft, kb, refsincone, cnx, cny, numr, mode, wr_four )
								#    match projections to its cone using an as a distance.
								for im in assignments[k]:
									# peak, pixer = proj_ali_incore_local(data[im], refrings, numr, xrng[N_step], yrng[N_step], step[N_step], an[N_step], finfo, sym = sym)
									# print numr,xrng,yrng,step,an, sym
									# sys.stdout.flush()
									peak, pixer = proj_ali_incore_local(proj_data[im], refrings, numr, xrng, yrng, step, an, finfo, sym = sym)
									proj_data[im].set_attr("previousmax", peak)

						print "Finished 2"

					dif_proj_angles = [None]*len(proj_data)
					dif_proj_angles_processed = [None]*len(proj_data)
					for im in range(len(dif_proj_angles)):
						a = get_params_proj(original_proj_data[im])[0]
						b = get_params_proj(original_proj_data_modified[im])[0]
						phi_diff = 360 - abs(a-b)%360 if (abs(a-b)%360 > 180) else abs(a-b)%360
						a = get_params_proj(original_proj_data[im])[1]
						b = get_params_proj(original_proj_data_modified[im])[1]
						theta_diff = 180 - abs(a-b)%180 if (abs(a-b)%180 > 90) else abs(a-b)%180
						dif_proj_angles[im] = phi_diff + theta_diff

						a = get_params_proj(original_proj_data[im])[0]
						b = get_params_proj(proj_data[im])[0]
						phi_diff = 360 - abs(a-b)%360 if (abs(a-b)%360 > 180) else abs(a-b)%360
						a = get_params_proj(original_proj_data[im])[1]
						b = get_params_proj(proj_data[im])[1]
						theta_diff = 180 - abs(a-b)%180 if (abs(a-b)%180 > 90) else abs(a-b)%180
						dif_proj_angles_processed[im] = phi_diff + theta_diff

					import numpy as np

					print
					print "Original and modified differences, Avg:", sum(dif_proj_angles)/len(dif_proj_angles)
					print "Original and modified differences, Max:", max(dif_proj_angles)
					print "Original and modified differences, Min:", min(dif_proj_angles)
					print "Original and modified differences, Percentiles:"
					for pp in range(10,100,10) + range(91,100,1) + [99.1 + x * 0.1 for x in range(0, 10)] : print pp,np.percentile(dif_proj_angles, pp)
					print
					print
					print "Original and processed differences, Avg:", sum(dif_proj_angles_processed)/len(dif_proj_angles_processed)
					print "Original and processed differences, Max:", max(dif_proj_angles_processed)
					print "Original and processed differences, Min:", min(dif_proj_angles_processed)
					print "Original and processed differences, Percentiles:"
					for pp in range(10,100,10) + range(91,100,1) + [99.1 + x * 0.1 for x in range(0, 10)] : print pp,np.percentile(dif_proj_angles_processed, pp)

					print

					for im in range(len(dif_proj_angles)):
						# if (abs(ran_phi_theta[im][0])+ abs(ran_phi_theta[im][1]))*1.3 < dif_proj_angles[im]:
						# if dif_proj_angles_processed[im] > 21:
							print "random phi theta displacement:", ran_phi_theta[im]
							print "Difference between original and modified: ",dif_proj_angles[im]
							print "Difference between original and processed: ",dif_proj_angles_processed[im]
							print "Original data:  ", get_params_proj(original_proj_data[im])
							print "Modification parameters:  ", modification_parameters[im]
							print "Modified data:  ", get_params_proj(original_proj_data_modified[im])
							print "Processed data: ", get_params_proj(proj_data[im])
							print


					print "============="
					for im in range(len(dif_proj_angles)):
						# if (abs(ran_phi_theta[im][0])+ abs(ran_phi_theta[im][1]))*1.3 < dif_proj_angles[im]:
						if dif_proj_angles_processed[im] > 21:
							print "random phi theta displacement:", ran_phi_theta[im]
							print "Difference between original and modified: ",dif_proj_angles[im]
							print "Difference between original and processed: ",dif_proj_angles_processed[im]
							print "Original data:  ", get_params_proj(original_proj_data[im])
							print "Modification parameters:  ", modification_parameters[im]
							print "Modified data:  ", get_params_proj(original_proj_data_modified[im])
							print "Processed data: ", get_params_proj(proj_data[im])
							print

					if redirect_to_stdout:
						sys.stdout.flush()
						sys.stdout.close()
						sys.stdout = old_stdout




def doCones_Horatio():

	import os
	os.chdir("/Users/hvoicu/EMAN2/src/test_scripts")



	# How my volume was generated
	# os.system("sxpdb2em.py tteftu_with_tRNA.pdb tmp.hdf --apix=12.2 --box=64")
	# os.system("e2proc3d.py tmp.hdf my_volume.hdf --process=filter.lowpass.tanh:cutoff_abs=0.25:fall_off=0.35")


	from utilities import get_symt, even_angles

	from math import sin, cos, radians
	from sys import exit
	from morphology import  bracket_def
	from utilities import assign_projangles
	from projection import prep_vol

	ref_a="S"
	# cnx=0
	# cny=0

	my_vol = EMData("my_volume.hdf")

	def computenumberofrefs(x, dat):
		#  dat = [sym, desired number of refs]
		return (len(even_angles(x,symmetry = dat[0])) - dat[1])**2

	# an = angular neighborhood
	xrng, yrng, step, an = (8,8,1,20)
	xrng, yrng, step, an = (0,0,0.5,20)

	volft,kb = prep_vol(my_vol)
	nx = my_vol.get_xsize()
	nz = my_vol.get_zsize()
	cnx = cny = nz//2 + 1

	numr = Numrinit(1,15)
	mode = "F"
	wr_four = ringwe(numr, mode)
	finfo = None


	delta = 0.5
	delta = 2
	delta = 4

	projangles = even_angles(delta,theta2=180.0)

	# proj_data = prgs(volft, kb, projangles)
	proj_data = project3d(my_vol, listagls=projangles)


	import copy

	original_proj_data = []
	for im in proj_data:
		original_proj_data.append(copy.deepcopy(im))

	import random

	# outside_an = [None]*len(proj_data)
	for idx, im in enumerate(proj_data):
		im.set_attr('previousmax', -1.0e23)
		im.set_attr("ID", idx)
		# horatio active_refactoring Jy51i1EwmLD4tWZ9_00000_1
		# im.set_attr("active", 1)
		random_phi = random.randint(3,4*an) - 2* an
		random_theta = random.randint(3,4*an) - 2*an
		# outside_an[idx] = int((abs(random_phi)>an) or abs(random_theta)>an)
		set_params_proj(im, [projangles[idx][0]+random.randint(3,30), projangles[idx][1] + random.randint(3,30), projangles[idx][2], 0, 0])
		# set_params_proj(im, [projangles[idx][0], projangles[idx][1] + 10, projangles[idx][2], 0, 0], "xform.anchor")

	original_proj_data_modified = []
	for im in proj_data:
		original_proj_data_modified.append(copy.deepcopy(im))

	sym = 'c1'
	#  Number of reference image that fit into the memory.  Will have to be hardwired.
	numberofrefs = 500

	totrefs = len(even_angles(delta,symmetry = sym))
	numberofcones = max(totrefs/numberofrefs,1)
	print " number of cones ",numberofcones,'   ',len(projangles)

	# if( numberofcones == 1):
	if 0:
		print "  One cone, i.e., standard code"
		refrings = prepare_refrings(volft, kb, nx, delta, ref_a, sym, numr, MPI = True)
		for im in xrange(len(projangles)):
			peak, pixer = proj_ali_incore_local(proj_data[im], refrings, numr, xrng, yrng, step, an, finfo, sym = sym)
			proj_data[im].set_attr("previousmax", peak)
	else:
		rs = delta
		h = 1.0
		dat = [sym, numberofcones]
		def1, def2 = bracket_def(computenumberofrefs,dat, rs, h)
		def1, val  = goldsearch_astigmatism(computenumberofrefs, dat, def1, def2, tol=1.0)
		print '  computed delta ',def1
		print 'computed number of cones  ',  len(even_angles(def1,symmetry = sym))
		coneangles = even_angles(def1, symmetry = sym)
		#  assign exclusively angl to coneangles
		assignments = assign_projangles(projangles, coneangles)
		print 'assignments  ',len(assignments)
		for k in xrange(len(assignments)):
			# print k,projangles[assignments[k][0]],projangles[assignments[k][-1]]
			for j in range(10):
				print k,projangles[assignments[k][j]]
			print
			for j in reversed(range(10)):
				print k,projangles[assignments[k][-j]]
			print
			print


			# projangles[assignments[k][-2]],projangles[assignments[k][-1]]
		#  for each coneangle compute refimages within 1.5*an degrees coneangle
	print "Finished 1"

	for k in xrange(len(coneangles)):
		if(len(assignments[k]) > 0):
			refsincone = even_angles(delta,symmetry = sym)
			ant = 1.5*an
			refsincone = cone_ang( refsincone, coneangles[k][0], coneangles[k][1], ant )
			refrings = refprojs( volft, kb, refsincone, cnx, cny, numr, mode, wr_four )
			#    match projections to its cone using an as a distance.
			for im in assignments[k]:
				# peak, pixer = proj_ali_incore_local(data[im], refrings, numr, xrng[N_step], yrng[N_step], step[N_step], an[N_step], finfo, sym = sym)
				# print numr,xrng,yrng,step,an, sym
				# sys.stdout.flush()
				peak, pixer = proj_ali_incore_local(proj_data[im], refrings, numr, xrng, yrng, step, an, finfo, sym = sym)
				proj_data[im].set_attr("previousmax", peak)

	print "Finished 2"

	dif_proj_angles = [None]*len(proj_data)
	dif_proj_angles = [None]*len(proj_data)
	for im in range(len(dif_proj_angles)):
		dif_proj_angles[im] = sum([abs(a - b) for a, b in zip(get_params_proj(original_proj_data[im]), get_params_proj(proj_data[im]))])
		dif_proj_angles[im] = sum([abs(a - b) for a, b in zip(get_params_proj(original_proj_data[im]), get_params_proj(proj_data[im]))])

	for im in range(200):
		print dif_proj_angles[im]
		print get_params_proj(original_proj_data[im])
		print get_params_proj(original_proj_data_modified[im])
		print get_params_proj(proj_data[im])
		print



def do_volume_mrk01(data, Tracker, iter, mpi_comm = None):
	"""
		data - projections (scattered between cpus) or the volume.  If volume, just do the volume processing
		options - the same for all cpus
		return - volume the same for all cpus
	"""
	from EMAN2          import Util
	from mpi            import mpi_comm_rank, MPI_COMM_WORLD
	from filter         import filt_table
	from reconstruction import recons3d_4nn_MPI, recons3d_4nn_ctf_MPI
	from utilities      import bcast_EMData_to_all
	import types

	if(mpi_comm == None):  mpi_comm = MPI_COMM_WORLD
	myid = mpi_comm_rank(mpi_comm)

	#=========================================================================
	# volume reconstruction
	if( type(data) == types.ListType ):
		if Tracker["constants"]["CTF"]:
			vol = recons3d_4nn_ctf_MPI(myid, data, Tracker["constants"]["snr"], \
					symmetry=Tracker["constants"]["sym"], npad=Tracker["constants"]["npad"], mpi_comm=mpi_comm)
		else:
			vol = recons3d_4nn_MPI    (myid, data,\
					symmetry=Tracker["constants"]["sym"], npad=Tracker["constants"]["npad"], mpi_comm=mpi_comm)
	else:
		vol = data

	if myid == 0:
		from morphology import threshold
		from filter     import filt_tanl, filt_btwl
		from utilities  import model_circle, get_im
		import types
		nx = vol.get_xsize()
		if(Tracker["constants"]["mask3D"] == None):
			mask3D = model_circle(int(Tracker["constants"]["radius"]*float(nx)/float(Tracker["constants"]["nnxo"])+0.5), nx, nx, nx)
		elif(Tracker["constants"]["mask3D"] == "auto"):
			from utilities import adaptive_mask
			mask3D = adaptive_mask(vol)
		else:
			if( type(Tracker["constants"]["mask3D"]) == types.StringType ):  mask3D = get_im(Tracker["constants"]["mask3D"])
			else:  mask3D = (Tracker["constants"]["mask3D"]).copy()
			nxm = mask3D.get_xsize()
			if( nx != nxm):
				from fundamentals import rot_shift3D
				mask3D = Util.window(rot_shift3D(mask3D,scale=float(nx)/float(nxm)),nx,nx,nx)
				nxm = mask3D.get_xsize()
				assert(nx == nxm)

		stat = Util.infomask(vol, mask3D, False)
		vol -= stat[0]
		Util.mul_scalar(vol, 1.0/stat[1])
		vol = threshold(vol)
		#Util.mul_img(vol, mask3D)
		if( Tracker["PWadjustment"] ):
			from utilities    import read_text_file
			from fundamentals import rops_table, fftip, fft
			rt = read_text_file( Tracker["PWadjustment"] )
			fftip(vol)
			ro = rops_table(vol)
			#  Here unless I am mistaken it is enough to take the beginning of the reference pw.
			for i in xrange(1,len(ro)):  ro[i] = (rt[i]/ro[i])**Tracker["upscale"]
			if( type(Tracker["lowpass"]) == types.ListType ):
				vol = fft( filt_table( filt_table(vol, Tracker["lowpass"]), ro) )
			else:
				vol = fft( filt_table( filt_tanl(vol, Tracker["lowpass"], Tracker["falloff"]), ro) )
		else:
			if( type(Tracker["lowpass"]) == types.ListType ):
				vol = filt_table(vol, Tracker["lowpass"])
			else:
				vol = filt_tanl(vol, Tracker["lowpass"], Tracker["falloff"])
		stat = Util.infomask(vol, mask3D, False)
		vol -= stat[0]
		Util.mul_scalar(vol, 1.0/stat[1])
		vol = threshold(vol)
		vol = filt_btwl(vol, 0.38, 0.5)
		Util.mul_img(vol, mask3D)
		del mask3D
		# vol.write_image('toto%03d.hdf'%iter)
	# broadcast volume
	bcast_EMData_to_all(vol, myid, 0, comm=mpi_comm)
	#=========================================================================
	return vol



#***************************************************************************************************
#*************************** horatio code start ****************************************************
#**************************  2016-01-05--01-34-37-250  *********************************************
#***************************************************************************************************






####################################################################################################
############################# start tNqKI3B53v1WXxq ################################################
####################################################################################################
#########  this block of functions is only used for     sali3d_base_h_01 ###########################


def save_object (obj, filename):
	import cPickle as pickle
	with open(filename, 'wb') as output:
		pickle.dump(obj, output, -1)

def load_object(filename):
	import cPickle as pickle
	with open(filename, 'rb') as output:
		obj = pickle.load(output)
	return obj

def individual_process(file_name_of_pickled_object_for_which_we_want_to_know_the_increase_in_process_memory_size):
	import gc, psutil, sys, os
	gc.disable()
	mem1 = psutil.Process(os.getpid()).get_memory_info()[0]
	my_object = load_object(file_name_of_pickled_object_for_which_we_want_to_know_the_increase_in_process_memory_size)
	mem2 = psutil.Process(os.getpid()).get_memory_info()[0]
	# print "mem2={:,}".format(mem2)
	# print "mem1={:,}".format(mem1)
	# print "mem2-mem1={:,}".format(mem2-mem1)
	sys.stdout.write("%ld"%(mem2-mem1))
	sys.stdout.flush()
	gc.enable()

def random_string(length_of_randomstring = 16):
	import random
	chars=map(chr, range(97, 123)) # a..z
	chars.extend(map(chr, range(65, 91))) # A..Z
	chars.extend(map(chr, range(48, 58))) # 0..9
	random_string = ""
	for i in xrange(length_of_randomstring):
		random_string += chars[random.randint(0,len(chars)-1)]
	return random_string

def total_size_of_object_in_memory(my_object):
	import inspect, os, subprocess

	file_name_my_object = random_string()
	while os.path.exists(file_name_my_object):
		file_name_my_object = random_string()

	save_object(my_object, file_name_my_object)

	file_name_my_python_code = random_string() + ".py"
	while os.path.exists(file_name_my_python_code):
		file_name_my_python_code = random_string() + ".py"

	fp = open(file_name_my_python_code, "w")
	fp.write("#!/usr/bin/env python\n\n")
	fp.write("from EMAN2 import *\n")
	fp.write("from sparx import *\n")

	for line in inspect.getsourcelines(load_object)[0]: fp.write(line)
	for line in inspect.getsourcelines(individual_process)[0]: fp.write(line)
	fp.write("individual_process('%s')"%file_name_my_object)
	fp.close()
	os.system("chmod +x ./%s"%file_name_my_python_code)

	import sys
	current_env = os.environ.copy()
	current_env['PYTHONPATH'] = ':'.join(sys.path)
	
	output = 0
	for i in xrange(10):
		# output += 0.1*int(subprocess.Popen(["./%s"%file_name_my_python_code], stdout = subprocess.PIPE, stderr = subprocess.STDOUT).communicate()[0])
		output += 0.1*int(subprocess.Popen(["./%s"%file_name_my_python_code], stdout = subprocess.PIPE, stderr = subprocess.STDOUT, env = current_env).communicate()[0])
	os.system("rm ./%s"%file_name_my_python_code)
	os.system("rm ./%s"%file_name_my_object)
	# print "output=", output
	return int(output) + 1


def determine_maximum_number_of_processes_per_node_from_all_nodes_that_belong_to_the_same_mpi_run():
	import os, socket
	from mpi import mpi_barrier, MPI_COMM_WORLD

	hostname = socket.gethostname()
	file_prefix = "WKDkSGYtLDTW9Nb2Vcu1SpsptFpEIod_mpi_process_count_"
	os.system("touch %s%s_%d"%(file_prefix, hostname, os.getpid()))
	mpi_barrier(MPI_COMM_WORLD)
	import glob
	list_of_files = glob.glob(file_prefix + "*")
	mpi_barrier(MPI_COMM_WORLD)
	hostname_list=[]
	for fn in list_of_files:
		hostname_list.append(fn[(len(file_prefix)):(len(file_prefix)+len(hostname))])
	from collections import Counter
	counter = Counter(hostname_list)
	os.system("rm %s%s_%d"%(file_prefix, hostname, os.getpid()))
	return max(counter.values())

def calculate_number_of_cones(volft, kb, delta, sym, cnx, cny, numr, mode, wr_four):

	import sys
	from alignment import prepare_refrings, refprojs, Numrinit, ringwe
	from morphology import bracket_def, goldsearch_astigmatism
	from applications import computenumberofrefs
	from utilities import even_angles, assign_projangles, cone_ang, print_from_process
	
	
	LOW_LIMIT_FOR_NUMBER_OF_REFERENCES_THAT_FIT_MEMORY = 100
	FRACTION_OF_MEMORY_THAT_CAN_BE_ALLOCATED = 0.9 # do not allocate all available memory
	FRACTION_OF_MEMORY_THAT_CAN_BE_ALLOCATED = 0.000125 # yields about 21 cones
	FRACTION_OF_MEMORY_THAT_CAN_BE_ALLOCATED = 0.000125/4 # yields about 103 cones
	LEAVE_THIS_FRACTION_OF_TOTAL_MEMORY_UNALLOCATED = 0.05  # for 64GB this represents about 3.2GB

	refsincone= even_angles(delta, symmetry = sym)

	total_number_of_references = len(refsincone)

	try:
		refrings = refprojs(volft, kb, refsincone[:LOW_LIMIT_FOR_NUMBER_OF_REFERENCES_THAT_FIT_MEMORY], cnx, cny, numr, mode, wr_four )
	except Exception:
		print "Not enough memory for allocating LOW_LIMIT_FOR_NUMBER_OF_REFERENCES_THAT_FIT_MEMORY. Exit."
		sys.exit()


	# from total_size_of_object_in_memory import total_size_of_object_in_memory
	refrings_memory_increase = total_size_of_object_in_memory(refrings)
	
	memory_for_one_item = refrings_memory_increase/LOW_LIMIT_FOR_NUMBER_OF_REFERENCES_THAT_FIT_MEMORY + 1

	import psutil
	machine_memory_that_can_be_allocated = psutil.avail_phymem() - (psutil.TOTAL_PHYMEM*LEAVE_THIS_FRACTION_OF_TOTAL_MEMORY_UNALLOCATED)
	machine_memory_that_can_be_allocated *= FRACTION_OF_MEMORY_THAT_CAN_BE_ALLOCATED

	error_status = [0]
	if machine_memory_that_can_be_allocated <= 0:
		print "Not enough memory for allocating refrings. Exit."
		error_status = [1]
		
	from mpi import mpi_reduce, MPI_INT, MPI_SUM, MPI_COMM_WORLD, mpi_comm_rank, mpi_comm_size
	from utilities import if_error_all_processes_quit_program
	error_status = mpi_reduce(error_status, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD)
	if_error_all_processes_quit_program(error_status)	

	number_of_concurrent_processes_per_node = determine_maximum_number_of_processes_per_node_from_all_nodes_that_belong_to_the_same_mpi_run()
	number_of_references_that_fit_in_memory = (machine_memory_that_can_be_allocated/number_of_concurrent_processes_per_node)/memory_for_one_item

	myid = mpi_comm_rank(MPI_COMM_WORLD)
	number_of_processes = mpi_comm_size(MPI_COMM_WORLD)
	
	all_cones_estimates = [0]*number_of_processes
	from math import ceil
	all_cones_estimates[myid] = max(int(ceil(total_number_of_references/number_of_references_that_fit_in_memory)),1)
	
	mpi_reduce(all_cones_estimates, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD)
	if myid == 0:
		number_of_cones_to_return = max(all_cones_estimates)
	else:
		number_of_cones_to_return = 0
		
	from mpi import mpi_bcast
	number_of_cones_to_return = mpi_bcast(number_of_cones_to_return, 1, MPI_INT, 0, MPI_COMM_WORLD)[0]
	return number_of_cones_to_return


def generate_indices_and_refrings(number_of_cones, nima, projangles, volft, kb, nx, delta, an, rangle, ref_a, sym, numr, MPI, phiEqpsi = "Zero"):
	
	from alignment import prepare_refrings, refprojs, Numrinit, ringwe, generate_list_of_reference_angles_for_search
	from alignment import reduce_indices_so_that_angles_map_only_to_asymmetrix_unit_and_keep_mirror_info
	from morphology import bracket_def, goldsearch_astigmatism
	from applications import computenumberofrefs
	from utilities import even_angles, assign_projangles_f, assign_projangles, cone_ang_f
	from utilities import cone_ang_with_index
	import sys
	from projection import prep_vol

	cnx = cny = nx//2 + 1
	# numr = Numrinit(1,15)
	mode = "F"
	wr_four = ringwe(numr, mode)

	if an <= 0.0:
		#=========================================================================
		# prepare reference angles
		ref_angles = even_angles(delta, symmetry=sym, method = ref_a, phiEqpsi = "Zero")
		#  Modify 0,0,0 s it can be properly inverted
		if( ref_angles[0][0] == 0.0  and ref_angles[0][1] == 0.0 ):
			ref_angles[0][0] = 0.01
			ref_angles[0][1] = 0.01
		if( rangle > 0.0 ):
			# shake
			from utilities import rotate_shift_params
			ref_angles = rotate_shift_params(anglelist, [ delta*rangle, delta*rangle, delta*rangle ])
		
		#=========================================================================
		# build references
		# volft, kb = prep_vol(vol)
		refrings = prepare_refrings(volft, kb, nx, delta, ref_angles, sym, numr, MPI=MPI, phiEqpsi = "Zero")
		del volft, kb
		#=========================================================================		


		# refrings = prepare_refrings(volft, kb, nx, delta, ref_a, sym, numr, MPI = False)
		list_of_reference_angles = \
			generate_list_of_reference_angles_for_search([[refrings[lr].get_attr("phi"), refrings[lr].get_attr("theta")] for lr in xrange(len(refrings))], sym=sym)
		
		
		# print "\nexiting through NO CONEs (an < 0.0) generate_indices_and_refrings\n"
		sys.stdout.flush()
		yield range(nima), refrings, list_of_reference_angles

	else:	
		number_of_cones = calculate_number_of_cones(volft, kb, delta, sym, cnx, cny, numr, mode, wr_four)
		
		# number_of_cones = 30
		# number_of_cones = 4
		# number_of_cones = 1
		if( number_of_cones == 1):
			print "  One cone, i.e., standard code"
			sys.stdout.flush()			
			
			ref_angles = even_angles(delta, symmetry=sym, method = ref_a, phiEqpsi = "Zero")
			#  Modify 0,0,0 s it can be properly inverted
			if( ref_angles[0][0] == 0.0  and ref_angles[0][1] == 0.0 ):
				ref_angles[0][0] = 0.01
				ref_angles[0][1] = 0.01
			if( rangle > 0.0 ):
				# shake
				from utilities import rotate_shift_params
				ref_angles = rotate_shift_params(anglelist, [ delta*rangle, delta*rangle, delta*rangle ])
			
			#=========================================================================
			# build references
			# volft, kb = prep_vol(vol)
			refrings = prepare_refrings(volft, kb, nx, delta, ref_angles, sym, numr, MPI, phiEqpsi = "Zero")
			del volft, kb
			#=========================================================================		
	
	
			# refrings = prepare_refrings(volft, kb, nx, delta, ref_a, sym, numr, MPI = False)
			list_of_reference_angles = \
				generate_list_of_reference_angles_for_search([[refrings[lr].get_attr("phi"), refrings[lr].get_attr("theta")] for lr in xrange(len(refrings))], sym=sym)
			
			yield range(nima), refrings, list_of_reference_angles

		else:
			# delta = 1.29
			# sym = "c5"
			# rs = delta; h = 1.0
			rs = delta; h = 0.1
			dat = [sym, number_of_cones, "P"]
			def computenumberofrefs(x, dat):
				return (len(even_angles(x,symmetry = dat[0])) - dat[1])**2
	
			def1, def2 = bracket_def(computenumberofrefs, dat, rs, h)
			if def1 == None:
				delta_cone = def2
			else:
				delta_cone, val  = goldsearch_astigmatism(computenumberofrefs, dat, def1, def2, tol=1.0)
			# coneangles = even_angles(delta_cone, theta2=180.0, symmetry = sym, method='P')
			coneangles = even_angles(delta_cone, symmetry = sym, method='P')
			# assignments = assign_projangles_f(projangles, coneangles)

			mapped_projangles = [[0.0, 0.0, 0.0] for i in xrange(len(projangles))]
			
			for i in xrange(len(projangles)):
				mapped_projangles[i][1] = projangles[i][1]
				if projangles[i][1] < 90:
					mapped_projangles[i][0] = projangles[i][0]%(360/int(sym[1]))
				else:
					mapped_projangles[i][0] = ((projangles[i][0]+180)%(360/int(sym[1])) + 180)%360
					# mapped_projangles[i][0] = (projangles[i][0] + 180)%(360/int(sym[1])) + 180

			#active
			assignments = assign_projangles(mapped_projangles, coneangles)
			largest_angles_in_cones = Util.get_largest_angles_in_cones(mapped_projangles, coneangles)
			
			number_of_cones = len(coneangles)
			
			# I0xDS5gejz3yqarg
			print "number_of_cones999:", number_of_cones
			
			all_refs_angles_within_asymmetric_unit = even_angles(delta, symmetry=sym, method = "S", phiEqpsi = "Zero")
			len_of_all_refs_angles_within_asymmetric_unit = len(all_refs_angles_within_asymmetric_unit)
			
			all_refs_angles_within_asymmetric_unit_plus_mirror_and_symmetries = generate_list_of_reference_angles_for_search(all_refs_angles_within_asymmetric_unit, sym)
			
			for k in xrange(len(coneangles)):
				if(len(assignments[k]) > 0):
					filtered_refsincone_plus_mirror_and_symmetries_with_original_index, original_index = \
					cone_ang_with_index(all_refs_angles_within_asymmetric_unit_plus_mirror_and_symmetries, coneangles[k][0], coneangles[k][1], min(largest_angles_in_cones[k] + an/2 + 1.5*delta, 180))

					reduced_original_index = [i % len_of_all_refs_angles_within_asymmetric_unit for i in original_index]
					set_of_reduced_original_index = sorted(list(set(reduced_original_index)))
					for i in xrange(len(filtered_refsincone_plus_mirror_and_symmetries_with_original_index)):
						filtered_refsincone_plus_mirror_and_symmetries_with_original_index[i] += \
						[set_of_reduced_original_index.index(reduced_original_index[i])]
					filtered_refsincone_plus_mirror_and_symmetries_with_original_index_and_refrings_index = \
						filtered_refsincone_plus_mirror_and_symmetries_with_original_index
					
					from mpi import MPI_COMM_WORLD, mpi_comm_rank 
					myid = mpi_comm_rank(MPI_COMM_WORLD)
					
					filtered_refsincone_plus_mirror_and_symmetries_with_original_index_and_refrings_index[0].append(len_of_all_refs_angles_within_asymmetric_unit)
						
					angles_used_to_generate_refrings = 	[all_refs_angles_within_asymmetric_unit[i] for i in set_of_reduced_original_index]
					refrings = prepare_refrings(volft, kb, nx, delta, angles_used_to_generate_refrings, sym, numr, MPI = False, phiEqpsi = "Zero")
					
					sys.stdout.flush()

					yield assignments[k], refrings, filtered_refsincone_plus_mirror_and_symmetries_with_original_index_and_refrings_index
				
				else:
					yield [],[],[]

	



# def generate_indices_and_refrings____PREVIOUS_VERSION(number_of_cones, nima, projangles, volft, kb, nx, delta, an, ref_a, sym, numr, MPI, phiEqpsi = "Zero"):
# 	
# 	from alignment import prepare_refrings, refprojs, Numrinit, ringwe, generate_list_of_reference_angles_for_search
# 	from alignment import reduce_indices_so_that_angles_map_only_to_asymmetrix_unit_and_keep_mirror_info
# 	from morphology import bracket_def, goldsearch_astigmatism
# 	from applications import computenumberofrefs
# 	from utilities import even_angles, assign_projangles_f, assign_projangles, cone_ang_f
# 	from utilities import cone_ang_with_index
# 	import sys
# 
# 	# print "\nentering generate_indices_and_refrings\n"
# 	sys.stdout.flush()
# 
# 
# 	cnx = cny = nx//2 + 1
# 	# numr = Numrinit(1,15)
# 	mode = "F"
# 	wr_four = ringwe(numr, mode)
# 
# 	if an <= 0.0:
# 		# refrings = prepare_refrings(volft, kb, nx, delta, ref_a, sym, numr, MPI = False)
# 		# # print "\nexiting through NO CONEs (an < 0.0) generate_indices_and_refrings\n"
# 		# sys.stdout.flush()
# 		# yield range(nima), refrings
# 
# 		all_refs_angles_within_asymmetric_unit = even_angles(delta, symmetry = sym)
# 		refrings = refprojs( volft, kb, all_refs_angles_within_asymmetric_unit, cnx, cny, numr, mode, wr_four )
# 
# 		from alignment import generate_list_of_reference_angles_for_search
# 		list_of_reference_angles = \
# 			generate_list_of_reference_angles_for_search(all_refs_angles_within_asymmetric_unit, sym=sym)
# 		
# 		yield range(nima), refrings, list_of_reference_angles
# 
# 
# 	else:	
# 		# number_of_cones = calculate_number_of_cones(volft, kb, delta, sym, cnx, cny, numr, mode, wr_four)
# 		number_of_cones = 20
# 		if( number_of_cones == 1):
# 			print "  One cone, i.e., standard code"
# 			refrings = prepare_refrings(volft, kb, nx, delta, ref_a, sym, numr, MPI = False)
# 			print "\nexiting through NO CONEs generate_indices_and_refrings\n"
# 			sys.stdout.flush()
# 			
# 			# all_refs_angles_within_asymmetric_unit = even_angles(delta, symmetry = sym)
# 			# refrings = refprojs( volft, kb, all_refs_angles_within_asymmetric_unit, cnx, cny, numr, mode, wr_four )
# 			# 
# 			# from alignment import generate_list_of_reference_angles_for_search
# 			# list_of_reference_angles = \
# 			# 	generate_list_of_reference_angles_for_search(all_refs_angles_within_asymmetric_unit, sym=sym)
# 			list_of_reference_angles = \
# 				generate_list_of_reference_angles_for_search([[refrings[lr].get_attr("phi"), refrings[lr].get_attr("theta")] for lr in xrange(len(refrings))], sym=sym)
# 			
# 			yield range(nima), refrings, list_of_reference_angles
# 		else:
# 			rs = delta; h = 1.0
# 			dat = [sym, number_of_cones, "P"]
# 			def computenumberofrefs(x, dat):
# 				return (len(even_angles(x,symmetry = dat[0])) - dat[1])**2
# 	
# 			def1, def2 = bracket_def(computenumberofrefs, dat, rs, h)
# 			delta_cone, val  = goldsearch_astigmatism(computenumberofrefs, dat, def1, def2, tol=1.0)
# 			# coneangles = even_angles(delta_cone, theta2=180.0, symmetry = sym, method='P')
# 			coneangles = even_angles(delta_cone, symmetry = sym, method='P')
# 			# assignments = assign_projangles_f(projangles, coneangles)
# 			assignments = assign_projangles(projangles, coneangles)
# 			largest_angles_in_cones = Util.get_largest_angles_in_cones(projangles, coneangles)
# 			# assert(len(coneangles) ==  number_of_cones)
# 			number_of_cones = len(coneangles)
# 			print "number_of_cones:", number_of_cones
# 
# 			all_refs_angles_within_asymmetric_unit = even_angles(delta, symmetry = sym)
# 			len_of_all_refs_angles_within_asymmetric_unit = len(all_refs_angles_within_asymmetric_unit)
# 			all_refs_angles_within_asymmetric_unit_plus_mirror_and_symmetries = generate_list_of_reference_angles_for_search(all_refs_angles_within_asymmetric_unit, sym)
# 			
# 			for k in xrange(len(coneangles)):
# 				if(len(assignments[k]) > 0):
# 					filtered_refsincone_plus_mirror_and_symmetries_with_original_index, original_index = \
# 					cone_ang_with_index(all_refs_angles_within_asymmetric_unit_plus_mirror_and_symmetries, coneangles[k][0], coneangles[k][1], largest_angles_in_cones[k] + an/2 + 1.5*delta)
# 					# refsincone = even_angles(delta, theta2=180.0, symmetry = sym)
# 					# print k, "largest_angles_in_cones[k]", largest_angles_in_cones[k]
# 					# aaa = cone_ang_f( all_refs_angles_within_asymmetric_unit_plus_mirror_and_symmetries, coneangles[k][0], coneangles[k][1], largest_angles_in_cones[k] + an)
# 					# filtered_refsincone_plus_mirror_and_symmetries_with_original_index, original_index = \
# 					# 	cone_ang_with_index(all_refs_angles_within_asymmetric_unit_plus_mirror_and_symmetries, coneangles[k][0], coneangles[k][1], min(180, 30 + largest_angles_in_cones[k] + an/2 + 1.5*delta))
# 
# 						# cone_ang_with_index(all_refs_angles_within_asymmetric_unit_plus_mirror_and_symmetries, coneangles[k][0], coneangles[k][1], 180)
# 						# cone_ang_with_index(all_refs_angles_within_asymmetric_unit_plus_mirror_and_symmetries, coneangles[k][0], coneangles[k][1], 80)
# 						# cone_ang_with_index(all_refs_angles_within_asymmetric_unit_plus_mirror_and_symmetries, coneangles[k][0], coneangles[k][1], largest_angles_in_cones[k] + an/2 + 1.5*delta)
# 					
# 					
# 					reduced_original_index = [i % len_of_all_refs_angles_within_asymmetric_unit for i in original_index]
# 					set_of_reduced_original_index = list(set(reduced_original_index))
# 					for i in xrange(len(filtered_refsincone_plus_mirror_and_symmetries_with_original_index)):
# 						filtered_refsincone_plus_mirror_and_symmetries_with_original_index[i] += \
# 						[set_of_reduced_original_index.index(reduced_original_index[i])]
# 					filtered_refsincone_plus_mirror_and_symmetries_with_original_index_and_refrings_index = \
# 						filtered_refsincone_plus_mirror_and_symmetries_with_original_index	
# 						
# 					filtered_refsincone_plus_mirror_and_symmetries_with_original_index_and_refrings_index[0].append(len_of_all_refs_angles_within_asymmetric_unit)	
# 						
# 					angles_used_to_generate_refrings = 	[all_refs_angles_within_asymmetric_unit[i] for i in set_of_reduced_original_index]
# 					
# 					# filtered_all_refs_angles_reduced_to_asymmetrix_unit_with_mirror_info =\
# 					# reduce_indices_so_that_angles_map_only_to_asymmetrix_unit_and_keep_mirror_info(all_refs_angles_within_asymmetric_unit, )
# 					
# 					# the mirror info does not affect the function
# 					refrings = refprojs( volft, kb, angles_used_to_generate_refrings, cnx, cny, numr, mode, wr_four )
# 					
# 					# print "\nexiting through CONEs generate_indices_and_refrings\n"
# 					yield assignments[k], refrings, filtered_refsincone_plus_mirror_and_symmetries_with_original_index_and_refrings_index
# 					# yield assignments[k], refrings, filtered_refrings_index_angles_with_mirror_info
# 			
# # ####################################################################################################################################
# # 			cone_neighborhood = 1.5
# # 			rs = delta
# # 			h = 1.0
# # 			# dat = [sym, number_of_cones]
# # 			dat = [sym, number_of_cones, "S"]
# # 			# print dat
# # 			
# # 			def computenumberofrefs(x, dat):
# # 				#  dat = [sym, desired number of refs]
# # 				return (len(even_angles(x,symmetry = dat[0])) - dat[1])**2
# # 	
# # 			def1, def2 = bracket_def(computenumberofrefs, dat, rs, h)
# # 			# def1 could be NoneType
# # 		
# # 			def1, val  = goldsearch_astigmatism(computenumberofrefs, dat, def1, def2, tol=1.0)
# # 			# print 'computed delta ',def1
# # 			# print 'computed number of cones  ',  len(even_angles(def1,symmetry = sym))
# # 			coneangles = even_angles(def1, symmetry = sym)
# # 			#  assign exclusively angl to coneangles
# # 			assignments = assign_projangles(projangles, coneangles)
# # 			# print 'assignments  ',len(assignments)
# # 			for k in xrange(len(coneangles)):
# # 				if(len(assignments[k]) > 0):
# # 					refsincone = even_angles(delta/8,symmetry = sym)
# # 					ant = cone_neighborhood*an
# # 					refsincone = cone_ang( refsincone, coneangles[k][0], coneangles[k][1], ant )
# # 					refrings = refprojs( volft, kb, refsincone, cnx, cny, numr, mode, wr_four )
# # 					print "\nexiting through CONEs generate_indices_and_refrings\n"
# # 					yield assignments[k], refrings 
# # ####################################################################################################################################


# updated from pawels version 2016-01-07--14-18-45-568
# def sali3d_base_h_01(stack, ref_vol = None, Tracker = None, mpi_comm = None, log = None):
def sali3d_base_h_01(stack, ref_vol = None, Tracker = None, rangle = 0.0, rshift = 0.0, mpi_comm = None, log = None):	
	"""
		parameters: list of (all) projections | reference volume is optional, the data is shrank, 
		  the program does not know anything about shrinking| ...
		Data is assumed to be CTF multiplied and the ctf_applied flag to be set.
		The alignment done depends on nsoft:
					 nsoft = 0 & an = -1: exhaustive deterministic
					 nsoft = 0 & an > 0 : local deterministic
					 nsoft = 1 shc
					 nsoft >1  shc_multi
		
	"""

	from alignment       import Numrinit, prepare_refrings
	from alignment       import proj_ali_incore,  proj_ali_incore_zoom,  proj_ali_incore_local, proj_ali_incore_local_zoom
	from alignment       import shc, center_projections_3D, ringwe
	from utilities       import bcast_number_to_all, bcast_EMData_to_all, 	wrap_mpi_gatherv, wrap_mpi_bcast, model_blank, print_from_process
	from utilities       import get_im, file_type, model_circle, get_input_from_string, get_params_proj, set_params_proj, pad
	from utilities       import even_angles
	from mpi             import mpi_bcast, mpi_comm_size, mpi_comm_rank, MPI_FLOAT, MPI_COMM_WORLD, mpi_barrier, mpi_reduce, MPI_INT, MPI_SUM
	from projection      import prep_vol
	from statistics      import hist_list
	from applications    import MPI_start_end
	from filter          import filt_ctf, filt_table
	from global_def      import Util
	from fundamentals    import resample, fshift
	from multi_shc       import shc_multi
	#from development     import do_volume_mrk01
	import user_functions
	from EMAN2           import EMUtil, EMData
	import types
	from time            import time

	nsoft            = Tracker["nsoft"]
	saturatecrit     = Tracker["saturatecrit"]
	pixercutoff      = Tracker["pixercutoff"]
	zoom             = Tracker["zoom"]
	center           = Tracker["constants"]["center"]
	CTF              = Tracker["constants"]["CTF"]
	ref_a            = Tracker["constants"]["ref_a"]
	rstep            = Tracker["constants"]["rs"]
	sym              = Tracker["constants"]["sym"]
	first_ring       = 1
	last_ring        = Tracker["radius"]
	xr               = Tracker["xr"]
	yr               = Tracker["yr"]
	ts               = Tracker["ts"]
	an               = Tracker["an"]
	delta            = Tracker["delta"]
	max_iter         = Tracker["maxit"]

	if mpi_comm == None: mpi_comm = MPI_COMM_WORLD

	if log == None:
		from logger import Logger
		log = Logger()

	number_of_proc = mpi_comm_size(mpi_comm)
	myid           = mpi_comm_rank(mpi_comm)
	main_node      = 0

	if myid == main_node:
		log.add("Start sali3d_base_h_01, nsoft = %1d"%nsoft)

	xrng        = get_input_from_string(xr)
	if  yr == "-1":  yrng = xrng
	else          :  yrng = get_input_from_string(yr)
	step        = get_input_from_string(ts)
	delta       = get_input_from_string(delta)
	lstp = min(len(xrng), len(yrng), len(step), len(delta))
	if an == "-1":
		an = [-1] * lstp
	else:
		an = get_input_from_string(an)

	if( type(stack) is types.StringType ):
		if myid == main_node:
			total_nima = EMUtil.get_image_count( stack )
		else:
			total_nima = 0
		total_nima = wrap_mpi_bcast(total_nima, main_node, mpi_comm)
		list_of_particles = range(total_nima)
		image_start, image_end = MPI_start_end(total_nima, number_of_proc, myid)
		# create a list of images for each node
		list_of_particles = list_of_particles[image_start: image_end]
		nima = len(list_of_particles)

	else:
		list_of_particles = range(len(stack))
		nima = len(list_of_particles)
		total_nima = len(list_of_particles)
		total_nima = mpi_reduce(total_nima, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD)
		total_nima = mpi_bcast(total_nima, 1, MPI_INT, 0, MPI_COMM_WORLD)
		total_nima = int(total_nima[0])


	if myid == 0 or myid == 5 or myid == 14:
		# finfo = None
		import os
		outdir = "./"
		# info_file = os.path.join(outdir, "progress%04d"%myid)
		
		# finfo = open(info_file, 'w')
		if "filename_first_time" not in sali3d_base_h_01.__dict__:
			import datetime
			sali3d_base_h_01.filename_first_time = os.path.join(outdir, "progress%04d_%s"%(myid, datetime.datetime.now().strftime('%Y-%m-%d--%H-%M-%S-%f')[:-3]))
			finfo = open(sali3d_base_h_01.filename_first_time, 'w')
		else: 
			finfo = open(sali3d_base_h_01.filename_first_time, 'a')
			finfo.write("=======================================================================================\n")
		# finfo = open(info_file, 'a')
	else:
		finfo = None

	if( myid == main_node):
		if( type(stack) is types.StringType ):  mask2D = get_im(stack, list_of_particles[0])
		else:                                   mask2D = stack[list_of_particles[0]]
		nx = mask2D.get_xsize()
	else:  nx = 0
	nx  = bcast_number_to_all(nx, source_node = main_node)
	mx = 2*nx
	if last_ring < 0:	last_ring = int(nx/2) - 2

	numr	= Numrinit(first_ring, last_ring, rstep, "F")

	data = [None]*nima
	for im in xrange(nima):
		if( type(stack) is types.StringType ):  data[im] = get_im(stack, list_of_particles[im])
		else:                                   data[im] = stack[list_of_particles[im]]
	mpi_barrier(mpi_comm)


	if myid == main_node:
		start_time = time()

	#  Read	template volume if provided or reconstruct it
	"""
	user_func = Tracker["constants"]["user_func"]
	if ref_vol:
		#vol = do_volume_mrk01(ref_vol, Tracker, 0, mpi_comm)
		ref_data = [ref_vol, Tracker, 0, mpi_comm]
		vol = user_func(ref_data)
	else:
		#vol = do_volume_mrk01(data, Tracker, 0, mpi_comm)
		ref_data = [data, Tracker, 0, mpi_comm]
		vol = user_func(ref_data)
	"""
	vol = ref_vol
	# log
	if myid == main_node:
		log.add("Setting of reference 3D reconstruction time = %10.1f\n"%(time()-start_time))
		start_time = time()


	pixer = [0.0]*nima
	#historyofchanges = [0.0, 0.5, 1.0]
	#par_r = [[] for im in list_of_particles ]
	cs = [0.0]*3
	total_iter = 0
	# do the projection matching,  it has a loop over iterations here, 
	#  but it can only do one iteration as many settings are done in meridien.  Perturbations are a good example, there is one per each iteration.
	if zoom: lstp = 1
	
	
	# import os
	# all_dirs = [d for d in os.listdir(".") if not os.path.isdir(d)]
	# import re; r = re.compile("^list_of_reference_angles_iter.*$")
	# all_dirs = filter(r.match, all_dirs)
	# if len(all_dirs) > 0:
	# 	from mpi import mpi_finalize
	# 	mpi_finalize()
	# 	import sys
	# 	sys.exit()
		
	# if "count_function_calls" not in sali3d_base_h_01.__dict__:
	# 	sali3d_base_h_01.count_function_calls = 0
	
	# sali3d_base_h_01.count_function_calls += 1
	
	for N_step in xrange(lstp):
		# calculate_number_of_cones(volft, kb, delta, sym, cnx, cny, numr, mode, wr_four)
		cnx = cny = nx//2 + 1
		# numr = Numrinit(1,15)
		mode = "F"
		wr_four = ringwe(numr, mode)
		volft, kb = prep_vol(vol)
		# number_of_cones = calculate_number_of_cones(volft, kb, delta[N_step], sym, cnx, cny, numr, mode, wr_four)
		number_of_cones = 4

		terminate = 0
		Iter = 0
		while Iter < max_iter and terminate == 0:

			Iter += 1
			total_iter += 1


			mpi_barrier(mpi_comm)
			if myid == main_node:
				log.add("ITERATION #%3d,  inner iteration #%3d"%(total_iter, Iter))
				log.add("Delta = %7.4f, an = %7.4f, xrange = %7.4f, yrange = %7.4f, step = %7.4f   %7.4f  %7.4f\n"%\
							(delta[N_step], an[N_step], xrng[N_step], yrng[N_step], step[N_step], rshift, rangle))
				start_time = time()



			#=========================================================================
			# build references
			volft, kb = prep_vol(vol)
			projangles = [[] for i in xrange(nima)]
			for im in xrange(nima):
				projangles[im] = get_params_proj(data[im])[:3]
					
					
			# print "XXXXXXXXXXXXXXXYYYYYYYYYYYYYZZZZZZZZZ", an[N_step]
			# sys.stdout.flush()
			# from mpi import mpi_finalize
			# mpi_finalize()
			# import sys
			# sys.exit()

			# start loop here ZCcw2oL8ZbcGbU		
			# print "nima", nima
			# from utilities import mpi_exit
			# mpi_exit()
			# for image_indices, refrings, filtered_all_refs_angles_reduced_to_asymmetrix_unit_with_mirror_info in generate_indices_and_refrings(number_of_cones, nima, projangles, volft, kb, nx, delta[N_step], an[N_step],
			cone_count = 0
			# an[N_step] = 1.5 * delta[N_step]
			for image_indices, refrings, list_of_reference_angles in generate_indices_and_refrings(number_of_cones, nima, projangles, volft, kb, nx, delta[N_step], an[N_step],	rangle, ref_a, sym, numr, MPI=mpi_comm, phiEqpsi = "Zero"):
				
				# if myid in [0,1,10]:
				# 	import json; f = open("list_of_reference_angles_func_call%03d_iter%02d_cone%d_myid%03d.json"%(sali3d_base_h_01.count_function_calls, Iter, cone_count, myid), 'w')
				# 	json.dump(list_of_reference_angles,f); f.close()
				cone_count += 1
				print "cone_count", cone_count

				# #=========================================================================
				# # build references
				# volft, kb = prep_vol(vol)
				# # refrings = prepare_refrings(volft, kb, nx, delta[N_step], ref_a, sym, numr, MPI=mpi_comm, phiEqpsi = "Zero")
				# refrings = prepare_refrings(volft, kb, nx, delta[N_step], ref_a, sym, numr, MPI=mpi_comm, phiEqpsi = "Zero")
				# del volft, kb
				# #=========================================================================
				
				# if myid == main_node:
				# 	from utilities import write_text_row
				# 	write_text_row(list_of_reference_angles, "list_of_reference_angles8.txt")
				# 	
				# 	import sys
				# 	sys.stdout.flush()
				# 
				# mpi_barrier(MPI_COMM_WORLD)
				# 
				# from mpi import mpi_finalize
				# mpi_finalize()
				# import sys
				# sys.exit()
				
				
				if myid == main_node:
					log.add("Time to prepare rings: %10.1f\n" % (time()-start_time))
					start_time = time()
	
				#=========================================================================
				#  there is no need for previousmax for deterministic searches
				if total_iter == 1 and nsoft > 0:
					if(an[N_step] < 0.0):
						# adjust params to references, calculate psi+shifts, calculate previousmax
						# for im in xrange(nima):
						for im in image_indices:
							previousmax = data[im].get_attr_default("previousmax", -1.0e23)
							if(previousmax == -1.0e23):
								peak, pixer[im] = proj_ali_incore_local(data[im], refrings, list_of_reference_angles, numr, \
										xrng[N_step], yrng[N_step], step[N_step], delta[N_step]*2.5, sym = sym)
								data[im].set_attr("previousmax", peak)
					else:
						#  Here it is supposed to be shake and bake for local SHC, but it would have to be signaled somehow
						for im in xrange(nima):
							data[im].set_attr("previousmax", -1.0e23)
					if myid == main_node:
						log.add("Time to calculate first psi+shifts+previousmax: %10.1f\n" % (time()-start_time))
						start_time = time()
				#=========================================================================
	
				# cannot have barriers in this loop because the number of cones might vary from process to process. Even if they are the same cones might be empty!
				# mpi_barrier(mpi_comm)
				if myid == main_node:  start_time = time()
				#=========================================================================
				# alignment
				#number_of_checked_refs = 0
				par_r = [0]*max(2,(nsoft+1))
				if(an[N_step] > 0):
					pass
					# these are already calculated by the generator at the top of the loop
					# generate list of angles
					# from alignment import generate_list_of_reference_angles_for_search
					# list_of_reference_angles = \
					# generate_list_of_reference_angles_for_search([[refrings[lr].get_attr("phi"), refrings[lr].get_attr("theta")] for lr in xrange(len(refrings))], sym=sym)			
				else:  list_of_reference_angles = [[1.0,1.0]]
				error_status = 0
				from utilities import if_error_all_processes_quit_program

				# for im in xrange(nima):
				for im in image_indices:
					if Tracker["constants"]["pwsharpening"] :
						#  High-pass filtration of data[im]
						try:
							stmp = data[im].get_attr("ptcl_source_image")
						except:
							try:
								stmp = data[im].get_attr("ctf")
								stmp = round(stmp.defocus,4)
							except:
								ERROR("Either ptcl_source_image or ctf has to be present in the header.","meridien",1, myid)
						try:
							indx = Tracker["bckgnoise"][1].index(stmp)
						except:
							ERROR("Problem with indexing ptcl_source_image.","meridien",1, myid)
	
						tempdata = Util.window(pad(filt_table(data[im],[Tracker["bckgnoise"][0][i,indx] for i in xrange(nx)]), mx, mx,1,0.0), nx, nx)
					else:  tempdata = data[im].copy()
					if(nsoft == 0):
						if(an[N_step] == -1):
							#  In zoom option each projection goes through shift zoom alignment
							if  zoom: peak, pixer[im] = proj_ali_incore_zoom(tempdata, refrings, numr, \
															xrng, yrng, step, finfo = finfo, sym=sym)
							else:  peak, pixer[im] = proj_ali_incore(tempdata, refrings, numr, \
													xrng[N_step], yrng[N_step], step[N_step], finfo = finfo, sym=sym, delta_psi = delta[N_step], rshift = rshift*xrng[N_step])
						else:
							if  zoom: peak, pixer[im] = proj_ali_incore_local_zoom(tempdata, refrings, list_of_reference_angles, numr, \
										xrng, yrng, step, an, finfo = finfo, sym=sym)
							else:  
								
								peak, pixer[im] = proj_ali_incore_local(tempdata, refrings, list_of_reference_angles, numr, \
										xrng[N_step], yrng[N_step], step[N_step], an[N_step], finfo = finfo, sym=sym, delta_psi = delta[N_step], rshift = rshift)
						if(pixer[im] == 0.0):  par_r[0] += 1
					elif(nsoft == 1):
						tempdata.set_attr("previousmax", data[im].get_attr("previousmax"))
						peak, pixer[im], number_of_checked_refs, iref = \
							shc(tempdata, refrings, list_of_reference_angles, numr, xrng[N_step], yrng[N_step], step[N_step], an[N_step], sym, finfo = finfo)
						if(pixer[im] == 0.0):  par_r[0] += 1
						data[im].set_attr("previousmax", tempdata.get_attr("previousmax"))
					elif(nsoft > 1):
						#  This is not functional
						peak, pixer[im], checked_refs, number_of_peaks = shc_multi(data[im], refrings, numr, \
													xrng[N_step], yrng[N_step], step[N_step], an[N_step], nsoft, sym, finfo = finfo)
						par_r[number_of_peaks] += 1
						#number_of_checked_refs += checked_refs
					data[im].set_attr("xform.projection", tempdata.get_attr("xform.projection"))
				if len(image_indices)>0: del tempdata
				if(an[N_step] > 0):  del list_of_reference_angles
				#=========================================================================
				# if_error_all_processes_quit_program(error_status)
				# cannot have barriers in this loop because the number of cones might vary from process to process. Even if they are the same cones might be empty!
				# mpi_barrier(mpi_comm)
				if myid == main_node:
					#print  data[0].get_attr_dict()
					log.add("Time of alignment = %10.1f\n"%(time()-start_time))
					start_time = time()
				
				# end loop here ZCcw2oL8ZbcGbU
				
			del volft, kb
			mpi_barrier(mpi_comm)
			print_from_process(0, "passed")
			#=========================================================================
			#  Pixer errors available here are useless as they are done for shifts on the reduced image size.
			#output pixel errors, check stop criterion
			all_pixer = wrap_mpi_gatherv(pixer, 0, mpi_comm)
			par_r = mpi_reduce(par_r, len(par_r), MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD)
			#total_checked_refs = wrap_mpi_gatherv([number_of_checked_refs], main_node, mpi_comm)
			terminate = 0
			if myid == main_node:
				#total_checked_refs = sum(total_checked_refs)
				if(nsoft < 2):  par_r[1] = total_nima - par_r[0]
				log.add("=========== Number of better orientations found ==============")
				for lhx in xrange(len(par_r)):
					msg = "            %5d     %7d"%(lhx, par_r[lhx])
					log.add(msg)
				log.add("_______________________________________________________")
				changes = par_r[0]/float(total_nima)
				"""  Have to think about it PAP
				if(  changes > saturatecrit ):
					if( Iter == 1 ):
						log.add("Will continue even though %4.2f images did not find better orientations"%saturatecrit)
					else:
						terminate = 1
						log.add("...............")
						log.add(">>>>>>>>>>>>>>>   Will terminate as %4.2f images did not find better orientations"%saturatecrit)
				if( terminate == 0 ):
					historyofchanges.append(changes)
					historyofchanges = historyofchanges[:3]
					historyofchanges.sort()
					if( (historyofchanges[-1]-historyofchanges[0])/2/(historyofchanges[-1]+historyofchanges[0]) <0.05 ):
						terminate = 1
						log.add("...............")
						log.add(">>>>>>>>>>>>>>>   Will terminate as orientations do not improve anymore")
				"""

				lhist = 20
				region, histo = hist_list(all_pixer, lhist)
				log.add("=========== Histogram of pixel errors ==============")
				for lhx in xrange(lhist):
					msg = "          %10.3f     %7d"%(region[lhx], histo[lhx])
					log.add(msg)
				log.add("____________________________________________________")
				if(nsoft<2 and terminate == 0):
					lhx = 0
					for msg in all_pixer:
						if(msg < pixercutoff): lhx += 1
					lhx = float(lhx)/float(total_nima)
					log.add(">>> %4.2f images had pixel error <%5.2f"%(lhx,pixercutoff))
					"""
					if( lhx > saturatecrit):
						if( Iter == 1 ):
							log.add("Will continue even though %4.2f images had pixel error < %5.2f"%(saturatecrit,pixercutoff))
						else:
							terminate = 1
							log.add("...............")
							log.add(">>>>>>>>>>>>>>>   Will terminate as %4.2f images had pixel error < %5.2f"%(saturatecrit,pixercutoff))
					"""
			terminate = True  #wrap_mpi_bcast(terminate, main_node, mpi_comm)
			#=========================================================================
			mpi_barrier(mpi_comm)
			if myid == main_node:
				#print  data[0].get_attr_dict()
				log.add("Time to compute histograms = %10.1f\n"%(time()-start_time))
				start_time = time()


			#=========================================================================
			mpi_barrier(mpi_comm)
			if( terminate or (Iter == max_iter) ):
				# gather parameters
				params = []
				for im in xrange(nima):
					t = get_params_proj(data[im])
					params.append( [t[0], t[1], t[2], t[3], t[4]] )
				params = wrap_mpi_gatherv(params, main_node, mpi_comm)
			# centering and volume reconstruction if not terminating
			else:
				#=========================================================================
				# centering
				if center == -1 and sym[0] == 'c':
					from utilities      import estimate_3D_center_MPI, rotate_3D_shift
					cs[0], cs[1], cs[2], dummy, dummy = estimate_3D_center_MPI(data, total_nima, myid, number_of_proc, main_node, mpi_comm=mpi_comm)
					if myid == main_node:
						msg = " Average center x = %10.3f        Center y = %10.3f        Center z = %10.3f\n"%(cs[0], cs[1], cs[2])
						log.add(msg)
					if int(sym[1]) > 1:
						cs[0] = cs[1] = 0.0
						if myid == main_node:
							log.add("For symmetry group cn (n>1), we only center the volume in z-direction\n")
					cs = mpi_bcast(cs, 3, MPI_FLOAT, main_node, mpi_comm)
					cs = [-float(cs[0]), -float(cs[1]), -float(cs[2])]
					rotate_3D_shift(data, cs)
				#=========================================================================
				if myid == main_node:
					start_time = time()
				#vol = do_volume_mrk01(data, Tracker, total_iter, mpi_comm)
				ref_data = [data, Tracker, total_iter, mpi_comm]
				user_func = Tracker["constants"] ["user_func"]
				vol = user_func(ref_data)
				#if myid == main_node:  vol.write_image('soft/smvol%04d.hdf'%total_iter)
				# log
				if myid == main_node:
					log.add("3D reconstruction time = %10.1f\n"%(time()-start_time))
					start_time = time()
			#=========================================================================

			"""
			#=========================================================================
			if(False):  #total_iter%1 == 5 or terminate):
				# gather parameters
				params = []
				previousmax = []
				for im in data:
					t = get_params_proj(im)
					params.append( [t[0], t[1], t[2], t[3]/shrinkage, t[4]/shrinkage] )
					#if(t[3] >0.0 or t[4]>0.0):  print  "  ERRROR  ",t
					previousmax.append(im.get_attr("previousmax"))
				assert(nima == len(params))
				params = wrap_mpi_gatherv(params, 0, mpi_comm)
				if myid == 0:
					assert(total_nima == len(params))
				previousmax = wrap_mpi_gatherv(previousmax, 0, mpi_comm)
				if myid == main_node:
					from utilities import write_text_row, write_text_file
					write_text_row(params, "soft/params%04d.txt"%total_iter)
					write_text_file(previousmax, "soft/previousmax%04d.txt"%total_iter)


				del previousmax, params
				i = 1
				while data[0].has_attr("xform.projection" + str(i)):
					params = []
					previousmax = []
					for im in data:

						try:
							#print  im.get_attr("xform.projection" + str(i))
							t = get_params_proj(im,"xform.projection" + str(i))
						except:
							print " NO XFORM  ",myid, i,im.get_attr('ID')
							from sys import exit
							exit()

						params.append( [t[0], t[1], t[2], t[3]/shrinkage, t[4]/shrinkage] )
					assert(nima == len(params))
					params = wrap_mpi_gatherv(params, 0, mpi_comm)
					if myid == 0:
						assert(total_nima == len(params))
					if myid == main_node:
						write_text_row(params, "soft/params-%04d-%04d.txt"%(i,total_iter))
					del previousmax, params
					i+=1


			if( ( terminate or (Iter == max_iter) ) and (myid == main_node) ):
				if( type(stack) is types.StringType ):
					from EMAN2 import Vec2f, Transform
					from EMAN2db import db_open_dict
					DB = db_open_dict(stack)
					for im in xrange(len(params)):
						t = Transform({"type":"spider","phi":params[im][0],"theta":params[im][1],"psi":params[im][2]})
						t.set_trans(Vec2f(-params[im][3], -params[im][4]))
						DB.set_attr(particle_ids[im], "xform.projection", t)
					DB.close()
				else:
					for im in xrange(len(params)): set_params_proj(stack[particle_ids[im]], params[im])
			"""


	if myid == main_node:
		log.add("Finish sali3d_base, nsoft = %1d"%nsoft)
	return params


# original 2016-01-07--14-17-51-354
# def sali3d_base_h_01(stack, ref_vol = None, Tracker = None, mpi_comm = None, log = None):
# 	"""
# 		parameters: list of (all) projections | reference volume is optional, the data is shrank, 
# 		  the program does not know anything about shrinking| ...
# 		Data is assumed to be CTF multiplied and the ctf_applied flag to be set.
# 		The alignment done depends on nsoft:
# 					 nsoft = 0 & an = -1: exhaustive deterministic
# 					 nsoft = 0 & an > 0 : local deterministic
# 					 nsoft = 1 shc
# 					 nsoft >1  shc_multi
# 		
# 	"""
# 
# 	from alignment       import Numrinit, prepare_refrings, proj_ali_incore,  proj_ali_incore_zoom,  proj_ali_incore_local, proj_ali_incore_local_zoom
# 	from alignment       import shc, center_projections_3D, ringwe
# 	from utilities       import bcast_number_to_all, bcast_EMData_to_all, 	wrap_mpi_gatherv, wrap_mpi_bcast, model_blank, print_from_process
# 	from utilities       import get_im, file_type, model_circle, get_input_from_string, get_params_proj, set_params_proj
# 	from mpi             import mpi_bcast, mpi_comm_size, mpi_comm_rank, MPI_FLOAT, MPI_COMM_WORLD, mpi_barrier, mpi_reduce, MPI_INT, MPI_SUM
# 	from projection      import prep_vol
# 	from statistics      import hist_list
# 	from applications    import MPI_start_end
# 	from filter          import filt_ctf
# 	from global_def      import Util
# 	from fundamentals    import resample, fshift
# 	from multi_shc       import shc_multi
# 	# from development     import do_volume_mrk01
# 	import user_functions
# 	from EMAN2           import EMUtil, EMData
# 	import types
# 	from time            import time
# 
# 	nsoft            = Tracker["nsoft"]
# 	saturatecrit     = Tracker["saturatecrit"]
# 	pixercutoff      = Tracker["pixercutoff"]
# 	zoom             = Tracker["zoom"]
# 	center           = Tracker["constants"]["center"]
# 	CTF              = Tracker["constants"]["CTF"]
# 	ref_a            = Tracker["constants"]["ref_a"]
# 	rstep            = Tracker["constants"]["rs"]
# 	sym              = Tracker["constants"]["sym"]
# 	first_ring       = 1
# 	last_ring        = Tracker["radius"]
# 	xr               = Tracker["xr"]
# 	yr               = Tracker["yr"]
# 	ts               = Tracker["ts"]
# 	an               = Tracker["an"]
# 	delta            = Tracker["delta"]
# 	max_iter         = Tracker["maxit"]
# 
# 	if mpi_comm == None: mpi_comm = MPI_COMM_WORLD
# 
# 	if log == None:
# 		from logger import Logger
# 		log = Logger()
# 
# 	number_of_proc = mpi_comm_size(mpi_comm)
# 	myid           = mpi_comm_rank(mpi_comm)
# 	main_node      = 0
# 
# 	if myid == main_node:
# 		log.add("Start sali3d_base_h_01, nsoft = %1d"%nsoft)
# 
# 	xrng        = get_input_from_string(xr)
# 	if  yr == "-1":  yrng = xrng
# 	else          :  yrng = get_input_from_string(yr)
# 	step        = get_input_from_string(ts)
# 	delta       = get_input_from_string(delta)
# 	lstp = min(len(xrng), len(yrng), len(step), len(delta))
# 	if an == "-1":
# 		an = [-1] * lstp
# 	else:
# 		an = get_input_from_string(an)
# 
# 	if( type(stack) is types.StringType ):
# 		if myid == main_node:
# 			total_nima = EMUtil.get_image_count( stack )
# 		else:
# 			total_nima = 0
# 		total_nima = wrap_mpi_bcast(total_nima, main_node, mpi_comm)
# 		list_of_particles = range(total_nima)
# 		image_start, image_end = MPI_start_end(total_nima, number_of_proc, myid)
# 		# create a list of images for each node
# 		list_of_particles = list_of_particles[image_start: image_end]
# 		nima = len(list_of_particles)
# 
# 	else:
# 		list_of_particles = range(len(stack))
# 		nima = len(list_of_particles)
# 		total_nima = len(list_of_particles)
# 		total_nima = mpi_reduce(total_nima, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD)
# 		total_nima = mpi_bcast(total_nima, 1, MPI_INT, 0, MPI_COMM_WORLD)
# 		total_nima = int(total_nima[0])
# 
# 
# 	if myid == 0:
# 		# finfo = None
# 		import os
# 		outdir = "./"
# 		info_file = os.path.join(outdir, "progress%04d"%myid)
# 		finfo = open(info_file, 'w')
# 	else:
# 		finfo = None
# 
# 	if( myid == main_node):
# 		if( type(stack) is types.StringType ):  mask2D = get_im(stack, list_of_particles[0])
# 		else:                                   mask2D = stack[list_of_particles[0]]
# 		nx = mask2D.get_xsize()
# 	else:  nx = 0
# 	nx  = bcast_number_to_all(nx, source_node = main_node)
# 	if last_ring < 0:	last_ring = int(nx/2) - 2
# 
# 	numr	= Numrinit(first_ring, last_ring, rstep, "F")
# 
# 	data = [None]*nima
# 	for im in xrange(nima):
# 		if( type(stack) is types.StringType ):  data[im] = get_im(stack, list_of_particles[im])
# 		else:                                   data[im] = stack[list_of_particles[im]]
# 	mpi_barrier(mpi_comm)
# 
# 	if myid == main_node:
# 		start_time = time()
# 
# 	#  Read	template volume if provided or reconstruct it
# 	#  Apply initfl first, meaning true fl has to be preserved
# 	#fl = Tracker["lowpass"]
# 	#Tracker["lowpass"] = Tracker["initialfl"]
# 	user_func = Tracker["constants"] ["user_func"]
# 	if ref_vol:
# 		# vol = do_volume_mrk01(ref_vol, Tracker, 0, mpi_comm)
# 		ref_data = [ref_vol, Tracker, 0, mpi_comm]
# 		vol = user_func(ref_data)
# 	else:
# 		# vol = do_volume_mrk01(data, Tracker, 0, mpi_comm)
# 		ref_data = [data, Tracker, 0, mpi_comm]
# 		vol = user_func(ref_data)
# 	#  Restore desired fl
# 	#Tracker["lowpass"] = fl
# 
# 	# log
# 	if myid == main_node:
# 		log.add("Setting of reference 3D reconstruction time = %10.1f\n"%(time()-start_time))
# 		start_time = time()
# 
# 
# 	pixer = [0.0]*nima
# 	historyofchanges = [0.0, 0.5, 1.0]
# 	#par_r = [[] for im in list_of_particles ]
# 	cs = [0.0]*3
# 	total_iter = 0
# 	# do the projection matching
# 	if zoom: lstp = 1
# 	
# 	
# 	# import os
# 	# all_dirs = [d for d in os.listdir(".") if not os.path.isdir(d)]
# 	# import re; r = re.compile("^list_of_reference_angles_iter.*$")
# 	# all_dirs = filter(r.match, all_dirs)
# 	# if len(all_dirs) > 0:
# 	# 	from mpi import mpi_finalize
# 	# 	mpi_finalize()
# 	# 	import sys
# 	# 	sys.exit()
# 		
# 	# if "count_function_calls" not in sali3d_base_h_01.__dict__:
# 	# 	sali3d_base_h_01.count_function_calls = 0
# 	
# 	# sali3d_base_h_01.count_function_calls += 1
# 	
# 	for N_step in xrange(lstp):
# 		# calculate_number_of_cones(volft, kb, delta, sym, cnx, cny, numr, mode, wr_four)
# 		cnx = cny = nx//2 + 1
# 		# numr = Numrinit(1,15)
# 		mode = "F"
# 		wr_four = ringwe(numr, mode)
# 		volft, kb = prep_vol(vol)
# 		# number_of_cones = calculate_number_of_cones(volft, kb, delta[N_step], sym, cnx, cny, numr, mode, wr_four)
# 		number_of_cones = 4
# 
# 		terminate = 0
# 		Iter = 0
# 		while Iter < max_iter and terminate == 0:
# 
# 			Iter += 1
# 			total_iter += 1
# 
# 			mpi_barrier(mpi_comm)
# 			if myid == main_node:
# 				log.add("ITERATION #%3d,  inner iteration #%3d"%(total_iter, Iter))
# 				log.add("Delta = %5.2f, an = %5.2f, xrange = %5d, yrange = %5d, step = %5.2f\n"%\
# 							(delta[N_step], an[N_step], xrng[N_step], yrng[N_step], step[N_step]))
# 				start_time = time()
# 
# 
# 
# 			#=========================================================================
# 			# build references
# 			volft, kb = prep_vol(vol)
# 			projangles = [[] for i in xrange(nima)]
# 			for im in xrange(nima):
# 				projangles[im] = get_params_proj(data[im])[:3]
# 					
# 					
# 			# print "XXXXXXXXXXXXXXXYYYYYYYYYYYYYZZZZZZZZZ", an[N_step]
# 			# sys.stdout.flush()
# 			# from mpi import mpi_finalize
# 			# mpi_finalize()
# 			# import sys
# 			# sys.exit()
# 
# 							
# 			# start loop here ZCcw2oL8ZbcGbU		
# 			# print "nima", nima
# 			# from utilities import mpi_exit
# 			# mpi_exit()
# 			# for image_indices, refrings, filtered_all_refs_angles_reduced_to_asymmetrix_unit_with_mirror_info in generate_indices_and_refrings(number_of_cones, nima, projangles, volft, kb, nx, delta[N_step], an[N_step],
# 			cone_count = 0
# 			for image_indices, refrings, list_of_reference_angles in generate_indices_and_refrings(number_of_cones, nima, projangles, volft, kb, nx, delta[N_step], an[N_step],	ref_a, sym, numr, MPI=mpi_comm, phiEqpsi = "Zero"):
# 				
# 				# if myid in [0,1,10]:
# 				# 	import json; f = open("list_of_reference_angles_func_call%03d_iter%02d_cone%d_myid%03d.json"%(sali3d_base_h_01.count_function_calls, Iter, cone_count, myid), 'w')
# 				# 	json.dump(list_of_reference_angles,f); f.close()
# 				cone_count += 1 
# 
# 				
# 				# #=========================================================================
# 				# # build references
# 				# volft, kb = prep_vol(vol)
# 				# # refrings = prepare_refrings(volft, kb, nx, delta[N_step], ref_a, sym, numr, MPI=mpi_comm, phiEqpsi = "Zero")
# 				# refrings = prepare_refrings(volft, kb, nx, delta[N_step], ref_a, sym, numr, MPI=mpi_comm, phiEqpsi = "Zero")
# 				# del volft, kb
# 				# #=========================================================================
# 				
# 				# if myid == main_node:
# 				# 	from utilities import write_text_row
# 				# 	write_text_row(list_of_reference_angles, "list_of_reference_angles8.txt")
# 				# 	
# 				# 	import sys
# 				# 	sys.stdout.flush()
# 				# 
# 				# mpi_barrier(MPI_COMM_WORLD)
# 				# 
# 				# from mpi import mpi_finalize
# 				# mpi_finalize()
# 				# import sys
# 				# sys.exit()
# 				
# 				
# 				if myid == main_node:
# 					log.add("Time to prepare rings: %10.1f\n" % (time()-start_time))
# 					start_time = time()
# 	
# 				#=========================================================================
# 				#  there is no need for previousmax for deterministic searches
# 				if total_iter == 1 and nsoft > 0:
# 					if(an[N_step] < 0.0):
# 						# adjust params to references, calculate psi+shifts, calculate previousmax
# 						# for im in xrange(nima):
# 						for im in image_indices:
# 							previousmax = data[im].get_attr_default("previousmax", -1.0e23)
# 							if(previousmax == -1.0e23):
# 								peak, pixer[im] = proj_ali_incore_local(data[im], refrings, list_of_reference_angles, numr, \
# 										xrng[N_step], yrng[N_step], step[N_step], delta[N_step]*2.5, sym = sym)
# 								data[im].set_attr("previousmax", peak)
# 					else:
# 						#  Here it is supposed to be shake and bake for local SHC, but it would have to be signaled somehow
# 						for im in xrange(nima):
# 							data[im].set_attr("previousmax", -1.0e23)
# 					if myid == main_node:
# 						log.add("Time to calculate first psi+shifts+previousmax: %10.1f\n" % (time()-start_time))
# 						start_time = time()
# 				#=========================================================================
# 	
# 				# mpi_barrier(mpi_comm)
# 				if myid == main_node:  start_time = time()
# 				#=========================================================================
# 				# alignment
# 				#number_of_checked_refs = 0
# 				par_r = [0]*max(2,(nsoft+1))
# 				# for im in xrange(nima):
# 				for im in image_indices:
# 					# from mpi import MPI_COMM_WORLD, mpi_comm_rank 
# 					# myid = mpi_comm_rank(MPI_COMM_WORLD)
# 					# if myid==0:
# 					# 	print "image_indices", image_indices
# 					# 	import sys
# 					# 	sys.stdout.flush()
# 					# mpi_barrier(MPI_COMM_WORLD)
# 					# from mpi import mpi_finalize
# 					# mpi_finalize()
# 					# import sys
# 					# sys.exit()
# 					
# 					if(nsoft == 0):
# 						if(an[N_step] == -1):
# 							#  In zoom option each projection goes through shift zoom alignment
# 							if  zoom: 
# 								
# 								# if cone_count == 3:
# 								# 	if myid == main_node:
# 								# 		peak, pixer[im] = proj_ali_incore_zoom(data[im], refrings, numr, \
# 								# 								xrng, yrng, step, sym=sym)
# 								# 	mpi_barrier(MPI_COMM_WORLD)
# 								# 	from mpi import mpi_finalize
# 								# 	mpi_finalize()
# 								# 	import sys
# 								# 	sys.exit()
# 								peak, pixer[im] = proj_ali_incore_zoom(data[im], refrings, numr, \
# 														xrng, yrng, step, sym=sym)
# 
# 								
# 							else:  peak, pixer[im] = proj_ali_incore(data[im], refrings, numr, \
# 															xrng[N_step], yrng[N_step], step[N_step], sym=sym)
# 						else:
# 							if  zoom: peak, pixer[im] = proj_ali_incore_local_zoom(data[im], refrings, list_of_reference_angles, numr, \
# 										xrng, yrng, step, an, finfo = finfo, sym=sym)
# 							else:  
# 								peak, pixer[im] = proj_ali_incore_local(data[im], refrings, list_of_reference_angles, numr, \
# 										xrng[N_step], yrng[N_step], step[N_step], an[N_step], finfo = finfo, sym=sym)
# 						if(pixer[im] == 0.0):  par_r[0] += 1
# 					elif(nsoft == 1):
# 						peak, pixer[im], number_of_checked_refs, iref = \
# 							shc(data[im], refrings, list_of_reference_angles, numr, xrng[N_step], yrng[N_step], step[N_step], an[N_step], sym, finfo = finfo)
# 						if(pixer[im] == 0.0):  par_r[0] += 1
# 					elif(nsoft > 1):
# 						peak, pixer[im], checked_refs, number_of_peaks = shc_multi(data[im], refrings, numr, \
# 													xrng[N_step], yrng[N_step], step[N_step], an[N_step], nsoft, sym, finfo = finfo)
# 						par_r[number_of_peaks] += 1
# 						#number_of_checked_refs += checked_refs
# 	
# 				#=========================================================================
# 				# mpi_barrier(mpi_comm)
# 				if myid == main_node:
# 					#print  data[0].get_attr_dict()
# 					log.add("Time of alignment = %10.1f\n"%(time()-start_time))
# 					start_time = time()
# 				
# 				# end loop here ZCcw2oL8ZbcGbU
# 				
# 			del volft, kb
# 			mpi_barrier(mpi_comm)
# 			print_from_process(0, "passed")
# 			#=========================================================================
# 			#output pixel errors, check stop criterion
# 			all_pixer = wrap_mpi_gatherv(pixer, 0, mpi_comm)
# 			par_r = mpi_reduce(par_r, len(par_r), MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD)
# 			#total_checked_refs = wrap_mpi_gatherv([number_of_checked_refs], main_node, mpi_comm)
# 			terminate = 0
# 			if myid == main_node:
# 				#total_checked_refs = sum(total_checked_refs)
# 				if(nsoft < 2):  par_r[1] = total_nima - par_r[0]
# 				log.add("=========== Number of better orientations found ==============")
# 				for lhx in xrange(len(par_r)):
# 					msg = "            %5d     %7d"%(lhx, par_r[lhx])
# 					log.add(msg)
# 				log.add("_______________________________________________________")
# 				changes = par_r[0]/float(total_nima)
# 				if(  changes > saturatecrit ):
# 					if( Iter == 1 ):
# 						log.add("Will continue even though %4.2f images did not find better orientations"%saturatecrit)
# 					else:
# 						terminate = 1
# 						log.add("...............")
# 						log.add(">>>>>>>>>>>>>>>   Will terminate as %4.2f images did not find better orientations"%saturatecrit)
# 				if( terminate == 0 ):
# 					historyofchanges.append(changes)
# 					historyofchanges = historyofchanges[:3]
# 					historyofchanges.sort()
# 					"""  Have to think about it PAP
# 					if( (historyofchanges[-1]-historyofchanges[0])/2/(historyofchanges[-1]+historyofchanges[0]) <0.05 ):
# 						terminate = 1
# 						log.add("...............")
# 						log.add(">>>>>>>>>>>>>>>   Will terminate as orientations do not improve anymore")
# 					"""
# 
# 				lhist = 20
# 				region, histo = hist_list(all_pixer, lhist)
# 				log.add("=========== Histogram of pixel errors ==============")
# 				for lhx in xrange(lhist):
# 					msg = "          %10.3f     %7d"%(region[lhx], histo[lhx])
# 					log.add(msg)
# 				log.add("____________________________________________________")
# 				if(nsoft<2 and terminate == 0):
# 					lhx = 0
# 					for msg in all_pixer:
# 						if(msg < pixercutoff): lhx += 1
# 					lhx = float(lhx)/float(total_nima)
# 					log.add(">>> %4.2f images had pixel error <%5.2f"%(lhx,pixercutoff))
# 					if( lhx > saturatecrit):
# 						if( Iter == 1 ):
# 							log.add("Will continue even though %4.2f images had pixel error < %5.2f"%(saturatecrit,pixercutoff))
# 						else:
# 							terminate = 1
# 							log.add("...............")
# 							log.add(">>>>>>>>>>>>>>>   Will terminate as %4.2f images had pixel error < %5.2f"%(saturatecrit,pixercutoff))
# 			terminate = wrap_mpi_bcast(terminate, main_node, mpi_comm)
# 			#=========================================================================
# 			mpi_barrier(mpi_comm)
# 			if myid == main_node:
# 				#print  data[0].get_attr_dict()
# 				log.add("Time to compute histograms = %10.1f\n"%(time()-start_time))
# 				start_time = time()
# 
# 
# 			#=========================================================================
# 			mpi_barrier(mpi_comm)
# 			if( terminate or (Iter == max_iter) ):
# 				# gather parameters
# 				params = []
# 				for im in xrange(nima):
# 					t = get_params_proj(data[im])
# 					params.append( [t[0], t[1], t[2], t[3], t[4]] )
# 				params = wrap_mpi_gatherv(params, main_node, mpi_comm)
# 			# centering and volume reconstruction if not terminating
# 			else:
# 				#=========================================================================
# 				# centering
# 				if center == -1 and sym[0] == 'c':
# 					from utilities      import estimate_3D_center_MPI, rotate_3D_shift
# 					cs[0], cs[1], cs[2], dummy, dummy = estimate_3D_center_MPI(data, total_nima, myid, number_of_proc, main_node, mpi_comm=mpi_comm)
# 					if myid == main_node:
# 						msg = " Average center x = %10.3f        Center y = %10.3f        Center z = %10.3f\n"%(cs[0], cs[1], cs[2])
# 						log.add(msg)
# 					if int(sym[1]) > 1:
# 						cs[0] = cs[1] = 0.0
# 						if myid == main_node:
# 							log.add("For symmetry group cn (n>1), we only center the volume in z-direction\n")
# 					cs = mpi_bcast(cs, 3, MPI_FLOAT, main_node, mpi_comm)
# 					cs = [-float(cs[0]), -float(cs[1]), -float(cs[2])]
# 					rotate_3D_shift(data, cs)
# 				#=========================================================================
# 				if myid == main_node:
# 					start_time = time()
# 				# vol = do_volume_mrk01(data, Tracker, total_iter, mpi_comm)
# 				ref_data = [data, Tracker, total_iter, mpi_comm]
# 				user_func = Tracker["constants"] ["user_func"]
# 				vol = user_func(ref_data)
# 				#if myid == main_node:  vol.write_image('soft/smvol%04d.hdf'%total_iter)
# 				# log
# 				if myid == main_node:
# 					log.add("3D reconstruction time = %10.1f\n"%(time()-start_time))
# 					start_time = time()
# 			#=========================================================================
# 
# 			"""
# 			#=========================================================================
# 			if(False):  #total_iter%1 == 5 or terminate):
# 				# gather parameters
# 				params = []
# 				previousmax = []
# 				for im in data:
# 					t = get_params_proj(im)
# 					params.append( [t[0], t[1], t[2], t[3]/shrinkage, t[4]/shrinkage] )
# 					#if(t[3] >0.0 or t[4]>0.0):  print  "  ERRROR  ",t
# 					previousmax.append(im.get_attr("previousmax"))
# 				assert(nima == len(params))
# 				params = wrap_mpi_gatherv(params, 0, mpi_comm)
# 				if myid == 0:
# 					assert(total_nima == len(params))
# 				previousmax = wrap_mpi_gatherv(previousmax, 0, mpi_comm)
# 				if myid == main_node:
# 					from utilities import write_text_row, write_text_file
# 					write_text_row(params, "soft/params%04d.txt"%total_iter)
# 					write_text_file(previousmax, "soft/previousmax%04d.txt"%total_iter)
# 
# 
# 				del previousmax, params
# 				i = 1
# 				while data[0].has_attr("xform.projection" + str(i)):
# 					params = []
# 					previousmax = []
# 					for im in data:
# 
# 						try:
# 							#print  im.get_attr("xform.projection" + str(i))
# 							t = get_params_proj(im,"xform.projection" + str(i))
# 						except:
# 							print " NO XFORM  ",myid, i,im.get_attr('ID')
# 							from sys import exit
# 							exit()
# 
# 						params.append( [t[0], t[1], t[2], t[3]/shrinkage, t[4]/shrinkage] )
# 					assert(nima == len(params))
# 					params = wrap_mpi_gatherv(params, 0, mpi_comm)
# 					if myid == 0:
# 						assert(total_nima == len(params))
# 					if myid == main_node:
# 						write_text_row(params, "soft/params-%04d-%04d.txt"%(i,total_iter))
# 					del previousmax, params
# 					i+=1
# 
# 
# 			if( ( terminate or (Iter == max_iter) ) and (myid == main_node) ):
# 				if( type(stack) is types.StringType ):
# 					from EMAN2 import Vec2f, Transform
# 					from EMAN2db import db_open_dict
# 					DB = db_open_dict(stack)
# 					for im in xrange(len(params)):
# 						t = Transform({"type":"spider","phi":params[im][0],"theta":params[im][1],"psi":params[im][2]})
# 						t.set_trans(Vec2f(-params[im][3], -params[im][4]))
# 						DB.set_attr(particle_ids[im], "xform.projection", t)
# 					DB.close()
# 				else:
# 					for im in xrange(len(params)): set_params_proj(stack[particle_ids[im]], params[im])
# 			"""
# 
# 
# 	if myid == main_node:
# 		log.add("Finish sali3d_base, nsoft = %1d"%nsoft)
# 	return params




#########  this block of functions is only used for     sali3d_base_h_01 ###########################
####################################################################################################
############################# end tNqKI3B53v1WXxq #################################################
####################################################################################################



#***************************************************************************************************
#*************************** horatio code end ****************************************************
#**************************  2016-01-05--01-34-37-250  *********************************************
#***************************************************************************************************


####################################################################################################
############################# begin obsolete code #################################################
####################################################################################################

def local_ali3d(stack, outdir, maskfile = None, ou = -1,  delta = 2, ts=0.25, center = -1, maxit = 10, 
           CTF = False, snr = 1.0, sym = "c1", chunk = -1.0, user_func_name = "ref_ali3d",
	     fourvar = True, npad = 4, debug = False, MPI = False):
	"""
		
	"""

	if MPI:
		local_ali3d_MPI(stack, outdir, maskfile, ou, delta, ts, center, maxit,
				CTF, snr, sym, chunk, user_func_name, 
				fourvar, npad, debug)
		return

	from alignment      import eqproj_cascaded_ccc
	from projection     import prep_vol
	from utilities      import model_circle, get_params_proj, set_params_proj
	from utilities      import get_image, drop_image
	from utilities      import amoeba_multi_level, rotate_3D_shift, estimate_3D_center
	from math           import pi
	from statistics     import fsc_mask
	from utilities      import print_begin_msg, print_end_msg, print_msg
	from EMAN2 import Processor
	import os 
	import sys

	if os.path.exists(outdir): ERROR('Output directory exists, please change the name and restart the program', "local_ali3d", 1)
	os.mkdir(outdir)
	import global_def
	global_def.LOGFILE =  os.path.join(outdir, global_def.LOGFILE)
	
	print_begin_msg('local_ali3d')

	import user_functions
	user_func = user_functions.factory[user_func_name]

	if CTF:
		ima = EMData()
		ima.read_image(stack, 0)
		ctf_applied = ima.get_attr("ctf_applied")
		del ima
		if ctf_applied == 1:  ERROR("local_ali3d does not work for CTF-applied data", "local_ali3d", 1)
		from reconstruction import recons3d_4nn_ctf
		from filter         import filt_ctf
	else   : from reconstruction import recons3d_4nn

	last_ring   = int(ou)
	max_iter    = int(maxit)
	center      = int(center)

	ima     = EMData()
	ima.read_image(stack, 0)
	nx      = ima.get_xsize()
	del ima
	if last_ring == -1:	last_ring = nx//2 - 2

	print_msg("Input stack                 : %s\n"%(stack))
	print_msg("Output directory            : %s\n"%(outdir))
	print_msg("Maskfile                    : %s\n"%(maskfile))
	print_msg("Outer radius                : %i\n"%(last_ring))
	print_msg("Angular search range        : %s\n"%(delta))
	print_msg("Translation search range    : %f\n"%(ts))
	print_msg("Maximum iteration           : %i\n"%(max_iter))
	print_msg("Center type                 : %i\n"%(center))
	print_msg("CTF correction              : %s\n"%(CTF))
	print_msg("Signal-to-Noise Ratio       : %f\n"%(snr))
	print_msg("Symmetry group              : %s\n"%(sym))
	if chunk <= 0.0:  chunk = 1.0
	print_msg("Chunk size                  : %f\n\n"%(chunk))
	print_msg("User function               : %s\n"%(user_func_name))
	
	if maskfile:
		import  types
		if type(maskfile) is types.StringType:  mask3D = get_image(maskfile)
		else:                                  mask3D = maskfile
	else:
		mask3D = model_circle(last_ring, nx, nx, nx)
	mask2D = model_circle(last_ring, nx, nx)


	if debug:  outf = file(os.path.join(outdir, "progress"), "w")
	else:      outf = None

	# horatio active_refactoring Jy51i1EwmLD4tWZ9_00000_1
	# active = EMUtil.get_all_attributes(stack, 'active')
	# list_of_particles = []
	# for im in xrange(len(active)):
	# 	if(active[im]):  list_of_particles.append(im)
	# del active
	
	nima = EMUtil.get_image_count(stack)
	list_of_particles = range(nima)
	
	dataim = EMData.read_images(stack, list_of_particles)
	nima = len(dataim)

	if debug:
		outf.write("  data read")
		outf.write("\n")
		outf.flush()

	n_of_chunks = int(1.0/chunk)
	
	if debug:
		outf = file(os.path.join(outdir, "progress"), "w")
		outf.write("  chunk = "+str(chunk)+"   ")
		outf.write("\n")
		outf.flush()
		outf.write("  chunk = "+str(n_of_chunks)+"   ")
		outf.write("\n")
		outf.flush()

	# initialize data for the reference preparation function
	ref_data = [mask3D, center, None, None]

	M = nx
	npad = 2
	N = M*npad
	K = 6
	alpha = 1.75
	r = M/2
	v = K/2.0/N
	params = {"filter_type": Processor.fourier_filter_types.KAISER_SINH_INVERSE, "alpha":alpha, "K":K, "r":r, "v":v, "N":N}

	data = [None]*7
	data[3] = mask2D
	cs = [0.0]*3

	for iteration in xrange(maxit+1):
		print_msg("ITERATION #%3d\n"%(iteration+1))
		for ic in xrange(n_of_chunks):
			if(center == -1):
				cs[0], cs[1], cs[2], dummy, dummy = estimate_3D_center(dataim)				
				rotate_3D_shift(dataim, [-cs[0], -cs[1], -cs[2]])
				msg = "Average center x = %10.3f        Center y = %10.3f        Center z = %10.3f\n"%(cs[0], cs[1], cs[2])
				print_msg(msg)				
			# compute updated 3D at the beginning of each chunk
			#  3D stuff
			if CTF: vol1 = recons3d_4nn_ctf(dataim, range(0, nima, 2), snr, 1, sym)
			else:	   vol1 = recons3d_4nn(dataim, range(0, nima, 2), sym, snr = snr)

			if CTF: vol2 = recons3d_4nn_ctf(dataim, range(1, nima, 2), snr, 1, sym)
			else:	   vol2 = recons3d_4nn(dataim, range(1, nima, 2), sym, snr = snr)

			# resolution
			fscc = fsc_mask(vol1, vol2, mask3D, 1.0, os.path.join(outdir, "resolution%04d"%(iteration*n_of_chunks+ic+1)))
			del vol1
			del vol2

			# calculate new and improved 3D
			if CTF: vol = recons3d_4nn_ctf(dataim, range(nima), snr, 1, sym)
			else:	   vol = recons3d_4nn(dataim, range(nima), sym, snr = snr)

			# store the reference volume
			drop_image(vol, os.path.join(outdir, "vol%04d.hdf"%(iteration*n_of_chunks+ic+1)))
			ref_data[2] = vol
			ref_data[3] = fscc

			#  call user-supplied function to prepare reference image, i.e., center and filter it
			vol, dummy = user_func(ref_data)

			drop_image(vol, os.path.join(outdir, "volf%04d.hdf"%(iteration*n_of_chunks+ic+1)))
			if(iteration == maxit):
				#  in last iteration quit here
				print_end_msg("local_ali3d")
				return

			if not CTF:
				data[0], data[1] = prep_vol(vol)

			image_start_in_chunk = ic*nima/n_of_chunks
			image_end_in_chunk   = (ic+1)*nima/n_of_chunks
			if debug:
				outf.write("image_start_in_chunk "+str(image_start_in_chunk)+"\n")
				outf.write("\n")
				outf.write("image_end_in_chunk "+str(image_end_in_chunk)+"\n")
				outf.write("\n")
				outf.flush()
			if CTF:  previous_defocus = -1.0
			for imn in xrange(image_start_in_chunk, image_end_in_chunk):
				if CTF:
					ctf_params = dataim[imn].get_attr( "ctf" )
					if ctf_params.defocus != previous_defocus:
						previous_defocus = ctf_params.defocus
						data[0], data[1] = prep_vol(filt_ctf(vol, ctf_params))

				data[2] = dataim[imn]

				if ts > 0.0:
					refi = dataim[imn].FourInterpol(nx*2, nx*2, 1, False)
					data[4] = Processor.EMFourierFilter(refi, params)
				
				phi, theta, psi, tx, ty = get_params_proj(dataim[imn])
				atparams = [phi, theta, psi]
				data[5] = [tx, ty]
				data[6] = ts
				data[5][0] *= -1
				data[5][1] *= -1

				if debug:
					initial, dummy  = eqproj_cascaded_ccc(atparams, data)  # this is if we need initial discrepancy
					outf.write("Image "+str(imn)+"\n")
					outf.write('Old  %8.3f  %8.3f  %8.3f  %8.3f  %8.3f  %11.4f'%(phi,theta,psi,tx,ty,initial))
					outf.write("\n")
					
				# change signs of shifts for projections
			
				weight_phi = max(delta, delta*abs((atparams[1]-90.0)/180.0*pi))

				optm_params = amoeba_multi_level(atparams, [weight_phi, delta, weight_phi], eqproj_cascaded_ccc, 1.e-4, 1.e-4, 500, data)
				optm_params[0].append(optm_params[3][0])
				optm_params[0].append(optm_params[3][1])
				optm_params[0][3] *= -1
				optm_params[0][4] *= -1
				
				if debug:
					outf.write('New  %8.3f  %8.3f  %8.3f  %8.3f  %8.3f  %11.4f  %4d'%(optm_params[0][0], optm_params[0][1], optm_params[0][2], optm_params[0][3], optm_params[0][4], optm_params[1], optm_params[2]))
					outf.write("\n")
					outf.flush()

				set_params_proj(dataim[imn], optm_params[0])

			#  here we write header infomation
			from utilities import write_headers
			#write_headers(stack, dataim, list_of_particles)

def local_ali3d_base_MPI(stack, templatevol, ali3d_options, shrinkage = 1.0,
		    	mpi_comm = None, log= None, chunk = -1.0, saturatecrit = 0.95, pixercutoff = 1.0, debug = False ):
	"""
		
	"""
	from alignment        import eqproj_cascaded_ccc
	from filter           import filt_ctf
	from projection       import prep_vol
	from fundamentals     import resample
	from utilities        import bcast_number_to_all, model_circle, get_params_proj, set_params_proj
	from utilities        import bcast_EMData_to_all, bcast_list_to_all, send_attr_dict, wrap_mpi_bcast, wrap_mpi_gatherv
	from utilities        import get_image, drop_image, file_type, get_im, get_input_from_string, model_blank
	from utilities        import amoeba_multi_level, rotate_3D_shift, estimate_3D_center_MPI
	from utilities        import print_begin_msg, print_end_msg, print_msg
	from multi_shc        import do_volume
	from statistics       import varf3d_MPI
	from math             import pi
	from mpi              import mpi_bcast, mpi_comm_size, mpi_comm_rank, MPI_FLOAT, MPI_COMM_WORLD, mpi_barrier
	from mpi              import mpi_reduce, MPI_INT, MPI_SUM
	from EMAN2 import Processor
	from EMAN2 import Vec2f, Transform
	import os
	import sys
	import types

	maxit  = ali3d_options.maxit
	ou     = ali3d_options.ou
	ts     = get_input_from_string(ali3d_options.ts)[0]
	delta  = get_input_from_string(ali3d_options.delta)[0]
	sym    = ali3d_options.sym
	sym    = sym[0].lower() + sym[1:]
	center = ali3d_options.center
	CTF    = ali3d_options.CTF
	fourvar = False



	if log == None:
		from logger import Logger
		log = Logger()


	if mpi_comm == None:
		mpi_comm = mpi_comm

	number_of_proc = mpi_comm_size(mpi_comm)
	myid           = mpi_comm_rank(mpi_comm)
	main_node = 0

	if myid == main_node:
		log.add("Start local_ali3d_base")

	#if os.path.exists(outdir):  ERROR('Output directory exists, please change the name and restart the program', "local_ali3d_MPI ", 1,myid)
	#mpi_barrier(mpi_comm)

	"""
	if myid == main_node:
		os.mkdir(outdir)
		import global_def
		global_def.LOGFILE =  os.path.join(outdir, global_def.LOGFILE)
		print_begin_msg("local_ali3d_MPI")
		import user_functions
		user_func = user_functions.factory[user_func_name]
		if CTF:
			ima = EMData()
			ima.read_image(stack, 0)
			ctf_applied = ima.get_attr_default("ctf_applied", 0)
			del ima
			if ctf_applied == 1:  ERROR("local_ali3d does not work for CTF-applied data", "local_ali3d_MPI", 1,myid)
	mpi_barrier(mpi_comm)
	"""
	if debug:
		from time import sleep
		while not os.path.exists(outdir):
			print  "Node ",myid,"  waiting..."
			sleep(5)

		finfo = open(os.path.join(outdir, "progress%04d"%myid), 'w')
	else:
		finfo = None

	last_ring   = int(ou)
	center      = int(center)

	if( type(stack) is types.StringType ):
		if myid == main_node:
			if(file_type(stack) == "bdb"):
				from EMAN2db import db_open_dict
				dummy = db_open_dict(stack, True)
			# horatio active_refactoring Jy51i1EwmLD4tWZ9_00000_1
			# active = EMUtil.get_all_attributes(stack, 'active')
			# list_of_particles = []
			# for im in xrange(len(active)):
			# 	if(active[im]):  list_of_particles.append(im)
			# del active
		
			nima = EMUtil.get_image_count(stack)
			list_of_particles = range(nima)
	
			total_nima = len(list_of_particles)
		else:
			list_of_particles = None
			total_nima = 0
		total_nima = wrap_mpi_bcast(total_nima, main_node, mpi_comm)
		total_nima = int(total_nima[0])
		list_of_particles = wrap_mpi_bcast(list_of_particles, main_node, mpi_comm)
		if myid == main_node:
			particle_ids = [0]*total_nima
			for i in xrange(total_nima):  particle_ids[i] = list_of_particles[i]
		image_start, image_end = MPI_start_end(total_nima, number_of_proc, myid)
		# create a list of images for each node
		list_of_particles = list_of_particles[image_start: image_end]
		nima = len(list_of_particles)

	else:
		list_of_particles = range(len(stack))
		nima = len(list_of_particles)
		total_nima = len(list_of_particles)
		total_nima = mpi_reduce(total_nima, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD)
		total_nima = mpi_bcast(total_nima, 1, MPI_INT, 0, MPI_COMM_WORLD)
		total_nima = int(total_nima[0])
		image_start = 0
		image_and   = nima

	if(myid == main_node):
		if( type(stack) is types.StringType ):  dataim = get_im(stack, list_of_particles[0])
		else:                                   dataim = stack[list_of_particles[0]]
		onx      = dataim.get_xsize()
		if(shrinkage == 1.0):  nx = onx
		else:		
			st = resample(dataim, shrinkage)
			nx = st.get_xsize()
	else:
		nx = 0
		onx = 0
	nx  = bcast_number_to_all(nx, source_node = main_node)
	onx = bcast_number_to_all(onx, source_node = main_node)


	if last_ring < 0:	last_ring = int(onx/2) - 2
	mask2D  = model_circle(last_ring, onx, onx)
	if(shrinkage < 1.0):
		last_ring  = int(last_ring*shrinkage)
		ali3d_options.ou = last_ring


	dataim = [None]*nima
	for im in xrange(nima):
		if( type(stack) is types.StringType ):  dataim[im] = get_im(stack, list_of_particles[im])
		else:                                   dataim[im] = stack[list_of_particles[im]]
		dataim[im].set_attr('ID', list_of_particles[im])
		ctf_applied = dataim[im].get_attr_default('ctf_applied', 0)
		if CTF :
			ctf_params = dataim[im].get_attr("ctf")
			if ctf_applied == 0:
				st = Util.infomask(dataim[im], mask2D, False)
				dataim[im] -= st[0]
			else:
				ERROR("Projection data cannot be CTF-applied","local_ali3d_base",1,myid)
		if(shrinkage != 1.0):
			phi,theta,psi,sx,sy = get_params_proj(dataim[im])
			dataim[im] = resample(dataim[im], shrinkage)
			sx *= shrinkage
			sy *= shrinkage
			set_params_proj(dataim[im], [phi,theta,psi,sx,sy])
			if CTF :
				ctf_params.apix /= shrinkage
				dataim[im].set_attr('ctf', ctf_params)

	mask2D  = model_circle(last_ring, nx, nx)


	if chunk <= 0.0:  chunk = 1.0
	n_of_chunks = int(1.0/chunk)

	"""
	if myid == main_node:
		import user_functions
		user_func = user_functions.factory[user_func_name]

		print_msg("Input stack                 : %s\n"%(stack))
		print_msg("Output directory            : %s\n"%(outdir))
		print_msg("Maskfile                    : %s\n"%(ali3d_options.mask3D))
		print_msg("Outer radius                : %i\n"%(last_ring))
		print_msg("Angular search range        : %s\n"%(delta))
		print_msg("Shift search range          : %f\n"%(ts))
		print_msg("Maximum iteration           : %i\n"%(max_iter))
		print_msg("Center type                 : %i\n"%(center))
		print_msg("CTF correction              : %s\n"%(CTF))
		print_msg("Signal-to-Noise Ratio       : %f\n"%(snr))
		print_msg("Symmetry group              : %s\n"%(sym))
		print_msg("Chunk size                  : %f\n\n"%(chunk))
		print_msg("User function               : %s\n"%(user_func_name))
	"""

	import  types
	if ali3d_options.mask3D:
		if type(ali3d_options.mask3D) is types.StringType:
			if myid == main_node:
				mask3D = get_im(ali3d_options.mask3D)
			else:
				mask3D = model_blank(nx, nx, nx)
		else:
			mask3D = ali3d_options.mask3D.copy()
		if myid == main_node:
			i = mask3D.get_xsize()
			if( shrinkage != 1.0 ):
				if( i != nx ):
					mask3D = resample(mask3D, shrinkage)
		bcast_EMData_to_all(mask3D, myid, main_node)
	else:
		mask3D = model_circle(last_ring, nx, nx, nx)

	#  Read	template volume if provided
	if templatevol:
		if type(templatevol) is types.StringType:
			if myid == main_node:
				vol = get_im(templatevol)
				i = vol.get_xsize()
				if( shrinkage != 1.0 ):
					if( i != nx ):
						vol = resample(vol, shrinkage)
			else:
				vol = model_blank(nx, nx, nx)
		else:
			if myid == main_node:
				i = templatevol.get_xsize()
				if( shrinkage != 1.0 ):
					if( i != nx ):
						vol = resample(templatevol, shrinkage)
				else:
					vol = templatevol.copy()
			else:
				vol = model_blank(nx, nx, nx)
		bcast_EMData_to_all(vol, myid, main_node)
		del templatevol
	else:
		vol = None

	if debug:
		finfo.write( "image_start, image_end: %d %d\n" %(image_start, image_end) )
		finfo.flush()

	if debug:
		finfo.write("  chunk = "+str(chunk)+"   ")
		finfo.write("\n")
		finfo.flush()
		finfo.write("  Number of chunks = "+str(n_of_chunks)+"   ")
		finfo.write("\n")
		finfo.flush()

	del list_of_particles

	if debug:
		finfo.write("  First image on this processor: "+str(image_start)+"   ")
		finfo.write("  Last  image on this processor: "+str(image_end)+"   ")
		finfo.write("\n")
		finfo.flush()

	if myid == main_node:
		# initialize data for the reference preparation function
		ref_data = [ mask3D, max(center,0), None, None, None ]
		# for method -1, switch off centering in user function
		ref_data.append( None )

	from time import time	
	if myid == main_node:
		log.add("Dimensions used (nx, onx, last_ring, shrinkage)  %5d    %5d     %5d     %6.3f\n"%(nx, onx, last_ring, shrinkage))
		start_time = time()


		
	M = nx
	npad = 2
	N = M*npad
	K = 6
	alpha = 1.75
	r = M/2
	v = K/2.0/N
	params = {"filter_type": Processor.fourier_filter_types.KAISER_SINH_INVERSE, "alpha":alpha, "K":K, "r":r, "v":v, "N":N}

	data = [None]*7
	data[3] = mask2D
	cs = [0.0]*3

	for iteration in xrange(maxit):
		if myid == main_node:
			start_time = time()
			log.add("ITERATION #%3d\n"%(iteration+1))
		if debug:
			finfo.write("  iteration = "+str(iteration)+"   ")
			finfo.write("\n")
			finfo.flush()
		pixer = [0.0]*nima
		for ic in xrange(n_of_chunks):
			# In the very first step the volume has to be computed if it was not provided by the user
			if( ((iteration > 0) and (ic > 0)) or vol == None):
				if(center == -1 and sym[0] == 'c'):
					if debug:
						finfo.write("  begin centering \n")
						finfo.flush()
					cs[0], cs[1], cs[2], dummy, dummy = estimate_3D_center_MPI(dataim, total_nima, myid, number_of_proc, main_node)
					cs = mpi_bcast(cs, 3, MPI_FLOAT, main_node, mpi_comm)
					cs = [-float(cs[0]), -float(cs[1]), -float(cs[2])]
					if int(sym[1]) > 1:
						cs[0] = cs[1] = 0.0
					rotate_3D_shift(dataim, cs)
					if myid == main_node:
						msg = "Average center x = %10.3f        Center y = %10.3f        Center z = %10.3f\n"%(cs[0], cs[1], cs[2])
						log.add(msg)
						log.add("Time to center = %d\n"%(time()-start_time))
						start_time = time()
				# compute updated 3D before each chunk
					# resolution
				if debug:
					finfo.write("  begin reconstruction = "+str(image_start))
					finfo.write("\n")
					finfo.flush()

				#  Do the 3D
				vol = do_volume(dataim, ali3d_options, iteration, mpi_comm)

				if myid == main_node:
					#drop_image(vol, os.path.join(outdir, "vol%03d_%03d.hdf"%(iteration, ic) ))
					log.add("3D reconstruction time = %d"%(time()-start_time))
					start_time = time()
				if debug:
					finfo.write("  done reconstruction = "+str(image_start))
					finfo.write("\n")
					finfo.flush()

				if fourvar:
				#  Compute Fourier variance
					varf = varf3d_MPI(dataim, ssnr_text_file = os.path.join(outdir, "ssnr%03d_%03d"%(iteration, ic)), mask2D = None, reference_structure = vol, ou = ou, rw = 1.0, npad = 1, CTF = CTF, sign = 1, sym =sym, myid = myid)
					if myid == main_node:
						varf = 1.0/varf
						print_msg("Time to calculate 3D Fourier variance= %d\n"%(time()-start_time))
						start_time = time()
				else:  varf = None

			if CTF:
				previous_defocus = -1.0
				#vol = fft(pad(vol, N, N, N))
			else:
				data[0], data[1] = prep_vol(vol)

			image_start_in_chunk = ic*nima/n_of_chunks
			image_end_in_chunk   = (ic+1)*nima/n_of_chunks
			if debug:
				finfo.write("Chunk "+str(ic)+"   Number of images in this chunk: "+str(image_end_in_chunk-image_start_in_chunk)+"\n")
				finfo.write("First image in this chunk: "+str(image_start_in_chunk)+"   Last image in this chunk: "+str(image_end_in_chunk-1)+"\n")
				finfo.flush()
			for imn in xrange(image_start_in_chunk, image_end_in_chunk):
				if CTF:
					ctf_params = dataim[imn].get_attr( "ctf" )
					if ctf_params.defocus != previous_defocus:
						previous_defocus = ctf_params.defocus
						data[0], data[1] = prep_vol(filt_ctf(vol, ctf_params))

				data[2] = dataim[imn]
				if ts > 0.0:
					refi = dataim[imn].FourInterpol(nx*2, nx*2, 1, True)
					data[4] = Processor.EMFourierFilter(refi, params)

				#phi, theta, psi, tx, ty = get_params_proj(dataim[imn])
				t1 = dataim[imn].get_attr("xform.projection")
				dp = t1.get_params("spider")
				atparams = [dp["phi"], dp["theta"], dp["psi"]]
				data[5]  = [dp["tx"], dp["ty"]]
				if debug:
					# we have to distiguish between no shift situation, which is done through ccc, and shift, which is done using gridding in 2D
					if(ts == 0.0):  data[6] = 0.0
					else:           data[6] = -1.0#ts#-1.0
					initial, dummy = eqproj_cascaded_ccc(atparams, data)  # this is if we need initial discrepancy
					finfo.write("Image "+str(imn)+"\n")
					finfo.write('Old  %6.1f  %6.1f  %6.1f   %5.2f  %5.2f  %11.4e\n'%(atparams[0],atparams[1],atparams[2], -dummy[0], -dummy[1], initial))
				# change signs of shifts for projections
				data[6] = ts
				#from random import random
				#data[5] = [(random()-0.5)*2,(random()-0.5)*2]  #  HERE !!!!!!!!!!!

				weight_phi = max(delta, delta*abs((atparams[1]-90.0)/180.0*pi))
				optm_params = amoeba_multi_level(atparams, [weight_phi, delta, weight_phi], eqproj_cascaded_ccc, 1.0, 1.e-2, 500, data)
				optm_params[0].append(optm_params[3][0])
				optm_params[0].append(optm_params[3][1])
				optm_params[0][3] *= -1
				optm_params[0][4] *= -1

				if debug:
					finfo.write('New  %6.1f  %6.1f  %6.1f   %5.2f  %5.2f  %11.4e  %4d\n'%(optm_params[0][0], optm_params[0][1], optm_params[0][2], optm_params[0][3], optm_params[0][4], optm_params[1], optm_params[2]))
					finfo.flush()

				#from sys import exit
				#exit()
				t2 = Transform({"type":"spider","phi":optm_params[0][0],"theta":optm_params[0][1],"psi":optm_params[0][2]})
				t2.set_trans(Vec2f(-optm_params[0][3], -optm_params[0][4]))
				dataim[imn].set_attr("xform.projection", t2)
				from pixel_error import max_3D_pixel_error
				pixer[imn] = max_3D_pixel_error(t1, t2, last_ring)
				#set_params_proj(dataim[imn], optm_params[0])
				#if( myid == main_node and imn%4 == 0):
				#	log.add( "Time to process %6d particles : %d\n" % (imn, time()-start_time) )
				#	start_time = time()
			if( myid == main_node ):
				log.add( "Time to process %6d particles : %d" % (image_end_in_chunk-image_start_in_chunk, time()-start_time) )
				start_time = time()

			# release memory
			data[0] = None


		#output pixel errors after all headers were processed
		from mpi import mpi_gatherv
		pixer = wrap_mpi_gatherv(pixer, main_node, mpi_comm)
		mpi_barrier(mpi_comm)
		terminate = 0
		if(myid == main_node):
			pixer = map(float, pixer)
			from statistics import hist_list
			lhist = 20
			region, histo = hist_list(pixer, lhist)
			log.add(" ")
			log.add("=========== Histogram of pixel errors ==============")
			for lhx in xrange(lhist):
				msg = "          %10.3f     %7d"%(region[lhx], histo[lhx])
				log.add(msg)
			log.add("____________________________________________________\n")


			# Terminate if saturatecrit% within pixercutoff pixel error
			im = 0
			for lhx in xrange(lhist):
				if(region[lhx] > pixercutoff): break
				im += histo[lhx]
			lhx = im/float(total_nima)
			if( lhx > saturatecrit):
				if( iteration == 1 ):
					log.add("First iteration, will continue even though %4.2f images did not find better orientations"%saturatecrit)
				else:
					terminate = 1
					log.add("..............................................................")
					log.add(">>>>>>>>>>>>>>>   Will terminate as %4.2f images did not find better orientations"%saturatecrit)
			del region, histo
		del pixer
		terminate = mpi_bcast(terminate, 1, MPI_INT, 0, mpi_comm)
		terminate = int(terminate[0])
		if terminate:  break


	del vol
	# gather parameters
	params = []
	for im in dataim:
		t = get_params_proj(im)
		params.append( [t[0], t[1], t[2], t[3]/shrinkage, t[4]/shrinkage] )
	params = wrap_mpi_gatherv(params, main_node, mpi_comm)

	if( myid == main_node ):
		"""
		if( type(stack) is types.StringType ):
			from EMAN2 import Vec2f, Transform
			from EMAN2db import db_open_dict
			DB = db_open_dict(stack)
			for im in xrange(len(params)):
				t = Transform({"type":"spider","phi":params[im][0],"theta":params[im][1],"psi":params[im][2]})
				t.set_trans(Vec2f(-params[im][3], -params[im][4]))
				DB.set_attr(particle_ids[im], "xform.projection", t)
			DB.close()
		else:
			for im in xrange(len(params)): set_params_proj(stack[particle_ids[im]], params[im])

		log.add("Time to write header information= %d\n"%(time()-start_time))
		"""
		log.add("Finish local_ ali3d_base")

	return  params


####################################################################################################
############################# end obsolete code #################################################
####################################################################################################
