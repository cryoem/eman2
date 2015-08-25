#
# Author: Pawel A.Penczek, 09/09/2006 (Pawel.A.Penczek@uth.tmc.edu)
# Copyright (c) 2000-2006 The University of Texas - Houston Medical School
#
# This software is issued under a joint BSD/GNU license. You may use the
# source code in this file under either license. However, note that the
# complete EMAN2 and SPARX software packages have some GPL dependencies,
# so you are responsible for compliance with the licenses of these packages
# if you opt to use BSD licensing. The warranty disclaimer below holfds
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

def iter_isac(stack, ir, ou, rs, xr, yr, ts, maxit, CTF, snr, dst, FL, FH, FF, init_iter, main_iter, iter_reali, \
			  match_first, max_round, match_second, stab_ali, thld_err, indep_run, thld_grp, img_per_grp, \
			  generation, candidatesexist = False, random_seed=None, new = False):
	from global_def   import ERROR, EMData, Transform
	from pixel_error  import multi_align_stability
	from utilities    import model_blank, write_text_file, get_params2D
	from utilities    import gather_EMData, bcast_EMData_to_all, send_EMData, recv_EMData
	from mpi          import mpi_comm_size, mpi_comm_rank, MPI_COMM_WORLD, MPI_FLOAT, MPI_INT
	from mpi          import mpi_bcast, mpi_barrier, mpi_send, mpi_recv, mpi_comm_split
	from random       import randint, seed
	from time         import localtime, strftime
	from applications import within_group_refinement
	import os

	number_of_proc = mpi_comm_size(MPI_COMM_WORLD)
	myid = mpi_comm_rank(MPI_COMM_WORLD)
	main_node = 0

	seed(myid)
	rand1 = randint(1,1000111222)
	seed(random_seed)
	rand2 = randint(1,1000111222)
	seed(rand1 + rand2)

	if main_iter%iter_reali != 0:
		ERROR("main_iter should be a multiple of iter_reali, please reset them and restart the program", "iter_isac", 1, myid)
	mpi_barrier(MPI_COMM_WORLD)

	if generation == 0:
		ERROR("Generation should begin from 1, please reset it and restart the program", "iter_isac", 1, myid)
	mpi_barrier(MPI_COMM_WORLD)

	if indep_run < 2 or indep_run > 4:
		ERROR("indep_run must equal 2, 3 or 4, please reset it and restart the program", "iter_isac", 1, myid)
	mpi_barrier(MPI_COMM_WORLD)

	if number_of_proc % indep_run != 0:
		ERROR("Number of MPI processes must be a multiplicity of indep_run, please reset it and restart the program", "iter_isac", 1, myid)
	mpi_barrier(MPI_COMM_WORLD)

	ali_params_dir = "ali_params_generation_%d"%generation
	if os.path.exists(ali_params_dir):  
		ERROR('Output directory %s for alignment parameters exists, please either change its name or delete it and restart the program'%ali_params_dir, "iter_isac", 1, myid)
	mpi_barrier(MPI_COMM_WORLD)

	
	if new: alimethod = "SHC"
	else:   alimethod = ""

	if myid == main_node:
		print "****************************************************************************************************"
		print "*                                                                                                  *"
		print "*                 Beginning of the ISAC program                "+strftime("%a, %d %b %Y %H:%M:%S", localtime())+"           *"
		print "*                                                                                                  *"
		print "* Iterative Stable Alignment and Clustering                                                        *"
		print "* By Zhengfan Yang, Jia Fang, Francisco Asturias and Pawel A. Penczek                              *"
		print "*                                                                                                  *"
		print '* REFERENCE: Z. Yang, J. Fang, J. Chittuluru, F. J. Asturias and P. A. Penczek, "Iterative Stable  *'
		print '*            Alignment and Clustering of 2D Transmission Electron Microscope Images",              *' 
		print '*            Structure 20, 237-247, February 8, 2012.                                              *'
		print "*                                                                                                  *"
		print "* Last updated: 07/23/2015 PAP                                                                     *"
		print "****************************************************************************************************"
		print "*                                       Generation %3d                                             *"%(generation)
		#print " alignment method  ",alimethod
		print "****************************************************************************************************"

	color = myid%indep_run
	key = myid/indep_run
	group_comm = mpi_comm_split(MPI_COMM_WORLD, color, key)
	group_main_node = 0

	# Read data on each processor, there are two ways, one is read on main_node and send them to all other nodes
	# The other way is all nodes reading it one by one, we have to test to determine which way is better.
	# The test shows that way 1 (18s) is way faster then way 2 (197s) on the test on 16 nodes.
	# The drawback of way 1 is it cannot have all attibutes, but I assume this is not important.

	# Method 1:
	if myid == main_node:
		alldata = EMData.read_images(stack)
		ndata = len(alldata)
		# alldata_n stores the original index of the particle (i.e., the index before running Generation 1)  
		alldata_n = [0]*ndata
		if generation > 1:
			for i in xrange(ndata): alldata_n[i] = alldata[i].get_attr('data_n')
		else:
			for i in xrange(ndata): alldata_n[i] = i
		nx = alldata[0].get_xsize()
	else:
		ndata = 0
		nx = 0
	ndata = mpi_bcast(ndata, 1, MPI_INT, main_node, MPI_COMM_WORLD)
	ndata = int(ndata[0])
	nx = mpi_bcast(nx, 1, MPI_INT, main_node, MPI_COMM_WORLD)
	nx = int(nx[0])

	if myid != main_node:
		alldata = [model_blank(nx, nx) for i in xrange(ndata)]
	mpi_barrier(MPI_COMM_WORLD)
	data = [None]*ndata
	tdummy = Transform({"type":"2D"})
	for im in xrange(ndata):
		bcast_EMData_to_all(alldata[im], myid, main_node)
		mpi_barrier(MPI_COMM_WORLD)  # has to be here, otherwise it chokes on our cluster.  PAP
		# This is the absolute ID, the only time we use it is
		# when setting the members of 4-way output. All other times, the id in 'members' is 
		# the relative ID.
		alldata[im].set_attr_dict({"xform.align2d": tdummy, "ID": im})
		data[im] = alldata[im]
	mpi_barrier(MPI_COMM_WORLD)
	'''
	# Method 2:
	alldata = EMData.read_images(stack)
	ndata = len(alldata)	
	# alldata_n stores the original index of the particle (i.e., the index before running Generation 1)  
	alldata_n = [0]*ndata
	if generation > 1:
		for i in xrange(ndata): alldata_n[i] = alldata[i].get_attr('data_n')
	else:
		for i in xrange(ndata): alldata_n[i] = i
	nx = alldata[0].get_xsize()
	data = [None]*ndata
	tdummy = Transform({"type":"2D"})
	for im in xrange(ndata):
		# This is the absolute ID, the only time we use it is
		# when setting the members of 4-way output. All other times, the id in 'members' is 
		# the relative ID.
		alldata[im].set_attr_dict({"xform.align2d": tdummy, "ID": im})
		data[im] = alldata[im]
	mpi_barrier(MPI_COMM_WORLD)
	'''

	ali_params_filename = "ali_params_%d"%color

		
	avg_num = 0
	Iter = 1
	match_initialization = False
	avg_first_stage = "class_averages_candidate_generation_%d.hdf"%generation

	if  not candidatesexist:
		if myid == main_node:
			print "******************************************************************************************"
			print "*            Beginning of the first phase           "+strftime("%a, %d %b %Y %H:%M:%S", localtime())+"            *"
			print "*                                                                                        *"
			print "* The first phase is an exploratory phase. In this phase, we set the criteria very       *"
			print "* loose and try to find as many candidate class averages as possible. This phase         *"
			print "* typically should have 10 to 20 rounds (default = 20). The candidate class averages are *"
			print "* stored in class_averages_candidate_generation_n.hdf.                                   *"
			print "******************************************************************************************"

		# I am adding here Artificial Intelligence for stopping 
		#  The program should stop if
		#	(a)  three times in a row it could not find new stable groups
		couldnt_find_stable = 0
		#	(b)  if number of groups to process is less than three
		K = ndata/img_per_grp

		while Iter <= max_round and couldnt_find_stable < 3 and K > 3:
			if myid == main_node: 
				print "################################################################################"
				print "#           Beginning of Round %2d           "%Iter+strftime("%a, %d %b %Y %H:%M:%S", localtime())+"          #"
				print "################################################################################"
				print "     Initialization of averages using EQ-mref"
				print "********************************************************************************"
				print "     We will process:  %d current images divided equally between %d groups"%(ndata, K)

			# Generate random averages for each group
			if key == group_main_node:
				refi = generate_random_averages(data, K, 9023)
				#refi = generate_random_averages(data, K, Iter)
				#refi = generate_random_averages(data, K, -1)
				###for j in xrange(len(refi)):  refi[j].write_image("refim_%d.hdf"%color, j)
			else:
				refi = [model_blank(nx, nx) for i in xrange(K)]

			for i in xrange(K):
				bcast_EMData_to_all(refi[i], key, group_main_node, group_comm)

			# Generate inital averages
			###if myid == main_node: print "	 Generating initial averages ",color,myid,localtime()[:5]
			refi = isac_MPI(data, refi, maskfile=None, outname=None, ir=ir, ou=ou, rs=rs, xrng=xr, yrng=yr, step=ts, 
					maxit=maxit, isac_iter=init_iter, CTF=CTF, snr=snr, rand_seed=-1, color=color, comm=group_comm, 
					stability=False, FL=FL, FH=FH, FF=FF, dst=dst, method = alimethod)

			# gather the data on main node
			if match_initialization:                #  This is not executed at all.  It was always this way, at least since version 1.1 by Piotr
				if key == group_main_node:          # as all refims are initialized the same way and also the flag is set to False!
					###print "Begin gathering ...", myid, len(refi)  #  It will append data
					refi = gather_EMData(refi, indep_run, myid, main_node)
				if myid == main_node:
					# Match all averages in the initialization and select good ones
					#print "before matching, len = ", len(refi)
					current_refim = get_unique_averages(refi, indep_run)
					# If data_good is too few, add some random ones, otherwise, cut to K
					###print " found data good = ", len(current_refim)
					if len(current_refim) > K:
						current_refim = current_refim[:K]
					elif len(current_refim) < K:
						defi = K - len(current_refim)
						for i in xrange(defi):
							current_refim.append(refi[randint(0, indep_run*K-1)].copy())
				else:
					current_refim = [model_blank(nx, nx) for i in xrange(K)]

				mpi_barrier(MPI_COMM_WORLD)	
			else:
				current_refim = refi

			# broadcast current_refim to all nodes
			for i in xrange(K):
				bcast_EMData_to_all(current_refim[i], myid, main_node)
			mpi_barrier(MPI_COMM_WORLD)

			###if key == group_main_node:
			###	for i in xrange(K):
			###		#  Each color has the same set of refim
			###		current_refim[i].write_image("init_group%d_round%d.hdf"%(color, Iter), i)

			# Run ISAC
			if myid == main_node:
				print "**********************************************************************"
				print "     Processing of candidate averages   "+strftime("%a, %d %b %Y %H:%M:%S", localtime())
				print "**********************************************************************"
	
			for mloop in xrange(1, match_first+1):
				if myid == main_node:
					print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
					print "     Loop %3d for 2-way matching   "%mloop+strftime("%a, %d %b %Y %H:%M:%S", localtime())
					print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
				refi = isac_MPI(data, current_refim, maskfile=None, outname=None, ir=ir, ou=ou, rs=rs, xrng=xr, yrng=yr, step=ts,
						maxit=maxit, isac_iter=main_iter, CTF=CTF, snr=snr, rand_seed=-1, color=color, comm=group_comm,
						stability=True, stab_ali=stab_ali, iter_reali=iter_reali, thld_err=thld_err, FL=FL, FH=FH, FF=FF, dst=dst, method = alimethod)

				all_ali_params = [[] for i in xrange(4)]
				for im in data:
					alpha, sx, sy, mirror, scale = get_params2D(im)
					all_ali_params[0].append(alpha)
					all_ali_params[1].append(sx)
					all_ali_params[2].append(sy)
					all_ali_params[3].append(mirror)
					#all_ali_params[4].append(scale)
				if key == group_main_node:
					final_ali_params_filename = ali_params_filename + "_" + str(mloop)
					#if os.path.exists(final_ali_params_filename):
					#	os.remove(final_ali_params_filename)
					write_text_file(all_ali_params, final_ali_params_filename)
				del all_ali_params

				# gather the data from the group main node to the main node
				if key == group_main_node:
					refi = gather_EMData(refi, indep_run, myid, main_node)

					###for i in xrange(len(refi)):
					###	#  Each color has the same set of refim
					###	refi[i].write_image("refi%d_round%d.hdf"%(color, Iter), i)

				if mloop != match_first:
					if myid == main_node:
						current_refim = match_2_way(data, refi, indep_run, thld_grp, FH, FF, suffix="_"+str(mloop) )
					else:
						current_refim = [model_blank(nx, nx) for i in xrange(K)]
					for k in xrange(K):
						bcast_EMData_to_all(current_refim[k], myid, main_node)
				mpi_barrier(MPI_COMM_WORLD)
			del current_refim

			# Run Matching
			if myid == main_node:
				###print " Before matching ...  ", color, myid,localtime()[:5] #len(data), len(refi), indep_run
				matched_data = match_2_way(data, refi, indep_run, thld_grp, FH, FF, suffix="_"+str(mloop) )
				members = []
				for im in matched_data:
					im.write_image(avg_first_stage, avg_num)
					avg_num += 1
					members.extend(im.get_attr('members'))

				# Because it's 2-way matching, it is possible some members are accounted for twice, we must delete the duplicate ones.   Yang 03/28/11
				members.sort()
				for i in xrange(len(members)-1, 0, -1):
					if members[i] == members[i-1]: del members[i]
				for i in xrange(len(members)-1): assert members[i]!=members[i+1]
				mem_len = len(members)
				print "In Round #%d, we found %d stable and reproducible averages, accounted for %d particles.  "%(Iter, len(matched_data), mem_len)
			else:
				mem_len = 0
				members = []
			mem_len = mpi_bcast(mem_len, 1, MPI_INT, main_node, MPI_COMM_WORLD)
			mem_len = int(mem_len[0])

			if mem_len > 0:
				# In members we have absolute ID
				members = mpi_bcast(members, mem_len, MPI_INT, main_node, MPI_COMM_WORLD)
				members = map(int, members)

				# Take out the good ones and use the remaining ones for initialization again
				nndata = ndata-len(members)
				newdata = [-1]*nndata
				ll = 0
				for i in xrange(ndata):
					abs_id = data[i].get_attr("ID")
					if abs_id not in members:
						newdata[ll] = abs_id
						ll += 1
				data = [alldata[im] for im in newdata]
				del newdata
				for im in data:
					im.set_attr("xform.align2d", tdummy)
				ndata = nndata

				couldnt_find_stable = 0
				K = ndata/img_per_grp
			else:
				couldnt_find_stable += 1
			Iter += 1
			mpi_barrier(MPI_COMM_WORLD)
		
		if myid == main_node:
			#  We will return after candidate averages are prepared so their calculation can be independently
			print "******************************************************************************************"
			print "*              End of the first phase             "+strftime("%a, %d %b %Y %H:%M:%S", localtime())+"              *"
			print "******************************************************************************************"
		return
	#  If candidates exist start from here
	refim_stack = avg_first_stage


	if myid == main_node:
		print ""
		print "******************************************************************************************"
		print "*           Beginning of the second phase         "+strftime("%a, %d %b %Y %H:%M:%S", localtime())+"             *"
		print "*                                                                                        *"
		print "* The second phase is where the actual class averages are generated, it typically has    *"
		print "* 3~9 iterations (default = 5) of matching. The first half of iterations are 2-way       *"
		print "* matchings, the second half of iterations are 3-way matchings, and the last iteration is*"
		print "* 4-way matching. In the second phase, three files will be generated:                    *"
		print "* class_averages_generation_n.hdf : class averages generated in this generation          *"
		print "* generation_n_accounted.txt      : IDs of accounted particles in this generation        *"
		print "* generation_n_unaccounted.txt    : IDs of unaccounted particles in this generation      *"
		print "******************************************************************************************"
		try:
			refim = EMData.read_images(refim_stack)
			print "* Using existing %4d candidate class averages                                            *"%len(refim)
		except:
			refim = []
		nrefim = len(refim)
	else:
		nrefim = 0
	nrefim = mpi_bcast(nrefim, 1, MPI_INT, main_node, MPI_COMM_WORLD)			# number of ref
	nrefim = int(nrefim[0])
	if(nrefim == 0):  ERROR("sxisac","Candidate averages do not exist",1,myid)

	if myid != main_node:
		refim = [model_blank(nx, nx) for i in xrange(nrefim)]
	mpi_barrier(MPI_COMM_WORLD)

	nn = [0]*nrefim
	for i in xrange(nrefim):
		bcast_EMData_to_all(refim[i], myid, main_node)						   # ref + n_objects
		if myid == main_node: n_objects = refim[i].get_attr('n_objects')
		else: n_objects = 0
		n_objects = mpi_bcast(n_objects, 1, MPI_INT, main_node, MPI_COMM_WORLD)
		nn[i] = int(n_objects[0])
		if myid != main_node:  refim[i].set_attr('n_objects', int(n_objects[0]))
		mpi_barrier(MPI_COMM_WORLD)
	nn.sort()

	if len(nn) > 0: img_per_grp = nn[-1]
	refim_all = refim

	two_way_loop = match_second/2
	ndata = len(alldata)
	K = ndata/img_per_grp
	for mloop in xrange(1, match_second+1):
		if mloop <= two_way_loop:
			wayness = 2
		elif mloop != match_second:
			if indep_run >= 3:
				wayness = 3
			else:
				wayness = 2
		else:
			if indep_run >= 4:
				wayness = 4
			else:
				wayness = indep_run
		if myid == main_node:		
			print "################################################################################"
			print "#       Iteration %2d for %d-way matching       "%(mloop, wayness)+strftime("%a, %d %b %Y %H:%M:%S", localtime())+"        #"
			print "################################################################################"

			members = []
			for im in refim_all:
				if im.get_attr('n_objects') > 1:
					members.extend(im.get_attr('members'))
			members.sort()
			for i in xrange(len(members)-1, 0, -1):
				if members[i] == members[i-1]: del members[i]
			n_members = len(members)
		else:
			n_members = 0
		n_members = mpi_bcast(n_members, 1, MPI_INT, main_node, MPI_COMM_WORLD)		   # n_members
		n_members = int(n_members[0])
		if myid != main_node:
			members = [0]*n_members
		members = mpi_bcast(members, n_members, MPI_INT, main_node, MPI_COMM_WORLD)	   # members
		members = map(int, members)

		ndata = len(alldata)
		nleft = ndata-n_members
		data_left = [None]*nleft
		c = 0
		for i in xrange(ndata):
			if i not in members:
				data_left[c] = alldata[i]
				c += 1

		K_left = nleft/img_per_grp
		if K_left > 0:
			if myid == main_node: 
				print "**********************************************************************"
				print "        Generating initial averages for unaccounted for images        "
				print "**********************************************************************"
				print "   Number of images unaccounted for = %d     Number of groups = %d"%(nleft, K_left)

			# Generate random averages for each group
			if key == group_main_node:
				refim_left = generate_random_averages(data_left, K_left)
				#for j in xrange(K_left):  refim_left[j].write_image("refim_left_%d.hdf"%color, j)
			else:
				refim_left = [model_blank(nx, nx) for i in xrange(K_left)]

			for i in xrange(K_left):
				bcast_EMData_to_all(refim_left[i], key, group_main_node, group_comm)		  # Within one SAC

			# Generate initial averages for the unaccounted images
			refim_left = isac_MPI(data_left, refim_left, maskfile=None, outname=None, ir=ir, ou=ou, rs=rs, xrng=xr, yrng=yr, step=ts, 
					maxit=maxit, isac_iter=init_iter, CTF=CTF, snr=snr, rand_seed=-1, color=color, comm=group_comm, stability=False, 
					FL=FL, FH=FH, FF=FF, dst=dst, method = alimethod)

			if len(refim) < K:
				# This will only happen in the first iteration, if applicable
				for k in xrange(K_left):
					refim.append(refim_left[k])
				refim = refim[:K]
			else:
				refim = refim[:K]
				ileft = 0
				for k in xrange(K):
					if refim[k].get_attr('n_objects') == 1:
						refim[k] = refim_left[ileft]
						ileft += 1
						if ileft >= K_left:  break
			mpi_barrier(MPI_COMM_WORLD)

#			if key == group_main_node:
#				for i in xrange(K):
#					refim[i].write_image("init_group%d_2nd_phase_round%d.hdf"%(color, mloop), i)

		# Run ISAC
		if myid == main_node:
			print "**********************************************************************"
			print "     Run the main part of ISAC program   "+strftime("%a, %d %b %Y %H:%M:%S", localtime())
			print "**********************************************************************"
			print "    Number of images = %d               Number of groups = %d"%(ndata, K)

		refim = isac_MPI(alldata, refim, maskfile=None, outname=None, ir=ir, ou=ou, rs=rs, xrng=xr, yrng=yr, step=ts, 
				maxit=maxit, isac_iter=main_iter, CTF=CTF, snr=snr, rand_seed=-1, color=color, comm=group_comm, 
				stability=True, stab_ali=stab_ali, iter_reali=iter_reali, thld_err=thld_err, FL=FL, FH=FH, FF=FF, dst=dst, method = alimethod)

		all_ali_params = [[] for i in xrange(4)]
		for im in alldata:
			alpha, sx, sy, mirror, scale = get_params2D(im)
			all_ali_params[0].append(alpha)
			all_ali_params[1].append(sx)			
			all_ali_params[2].append(sy)
			all_ali_params[3].append(mirror)
			#all_ali_params[4].append(scale)
		if key == group_main_node:
			final_ali_params_filename = ali_params_filename + "_" + str(mloop)
			#if os.path.exists(final_ali_params_filename):
			#	os.remove(final_ali_params_filename)
			write_text_file(all_ali_params, final_ali_params_filename)

		# gather refim to the main node
		if key == group_main_node:
			refim = gather_EMData(refim, indep_run, myid, main_node)
#			for i in xrange(len(refim)):
#				refim[i].write_image("log_mainPart_" + str(color) + "_" + str(mloop) + ".hdf", i)

		if mloop != match_second:
			if myid == main_node:
				print "**********************************************************************"
				print "     Run the %d-way matching algorithm  "%wayness+strftime("%a, %d %b %Y %H:%M:%S", localtime())
				print "**********************************************************************"
				# In this last two-way loop, we find all unique 2-way matches and use it as the starting
				# point of three-way match
				if mloop == two_way_loop:
					refim_all = match_2_way(alldata, refim, indep_run, thld_grp, FH, FF, suffix="_"+str(mloop) )
					# If they are enough, good; otherwise, add some random images into it.
					if len(refim_all) > K:
						print "Since the number of unique 2-way matches is larger than the number of groups (%d), we only use the first %d of them."%(K, K)
						refim_all = refim_all[:K]
					elif len(refim_all) < K:
						defi = K - len(refim_all)
						print "Since the number of unique 2-way matches is smaller than the number of groups (%d), we have to append %d random images."%(K, defi)
						for i in xrange(defi):
							# put some dummy avgs here
							temp_id = randint(0, ndata-1)
							ave = alldata[temp_id].copy()
							ave.set_attr_dict({"members": [temp_id], "n_objects": 1})
							refim_all.append(ave)
					for i in xrange(K*(indep_run-1)):
						refim_all.append(refim_all[i%K])
				else:
					refim_all = match_2_way(alldata, refim, indep_run, thld_grp, FH, FF, find_unique = False, wayness = wayness, suffix="_"+str(mloop) )
			else:
				refim_all = [model_blank(nx, nx) for i in xrange(K*indep_run)]
			for k in xrange(K*indep_run):
				bcast_EMData_to_all(refim_all[k], myid, main_node)
				if myid == main_node: n_objects = refim_all[k].get_attr('n_objects')
				else: n_objects = 0
				n_objects = mpi_bcast(n_objects, 1, MPI_INT, main_node, MPI_COMM_WORLD)
				if myid != main_node:  refim_all[k].set_attr('n_objects', int(n_objects[0]))
			refim = refim_all[color*K:(color+1)*K]
#			if key == group_main_node:
#				for k in xrange(K): refim[k].write_image("%d_way_match_%02d_%02d.hdf"%(wayness, mloop, color), k)
		mpi_barrier(MPI_COMM_WORLD)
#		if key == group_main_node:
#			for i in xrange(len(refim)):
#				refim[i].write_image("log_afterMatching_" + str(color) + "_" + str(mloop) + ".hdf", i)

	if key == group_main_node:
		final_ali_params_filename = ali_params_filename + "_" + str(mloop)
		if os.path.exists(final_ali_params_filename):
			os.remove(final_ali_params_filename)

	if myid == main_node:
		print "**********************************************************************"
		print "       Run the final %d-way matching algorithm  "%indep_run+strftime("%a, %d %b %Y %H:%M:%S", localtime())
		print "**********************************************************************"

	# Run 4-way Matching
	# Comment by Zhengfan Yang on 6/20/11
	# The original design was way too slow, we have to send all stable sets to each node and let each node to do the realignment
	# and send them back, even though the code will be much more complicated.
	# I have decided that main node should not do realignment, otherwise it could clog the whole operation if it happened to have
	# a very large group.  The main node is used to send and collect information.

	if myid == main_node:
		STB_PART = match_independent_runs(alldata, refim, indep_run, thld_grp)
		l_STB = len(STB_PART)
		os.mkdir(ali_params_dir)
		print  "  l_STB   ",l_STB
	else:
		l_STB = 0
		from time import sleep
		while not os.path.exists(ali_params_dir):   sleep(5)
	mpi_barrier(MPI_COMM_WORLD)
	l_STB = mpi_bcast(l_STB, 1, MPI_INT, main_node, MPI_COMM_WORLD)
	l_STB = int(l_STB[0])

	if myid == main_node:
		for i in xrange(l_STB):
			node_to_run = i%(number_of_proc-1)+1
			mpi_send(len(STB_PART[i]), 1, MPI_INT, node_to_run, i+10000, MPI_COMM_WORLD)
			mpi_send(STB_PART[i], len(STB_PART[i]), MPI_INT, node_to_run, i+20000, MPI_COMM_WORLD)

		members_acc = []
		ave_num = 0
		for i in xrange(l_STB):
			node_to_run = i%(number_of_proc-1)+1
			l_stable_members = mpi_recv(1, MPI_INT, node_to_run, i+30000, MPI_COMM_WORLD)
			l_stable_members = int(l_stable_members[0])
			stable_members = mpi_recv(l_stable_members, MPI_INT, node_to_run, i+40000, MPI_COMM_WORLD)
			stable_members = map(int, stable_members)
			mirror_consistent_rate = mpi_recv(1, MPI_FLOAT, node_to_run, i+50000, MPI_COMM_WORLD)
			mirror_consistent_rate = float(mirror_consistent_rate[0])
			pix_err = mpi_recv(1, MPI_FLOAT, node_to_run, i+60000, MPI_COMM_WORLD)
			pix_err = float(pix_err[0])

			print  "Group %d ...... Mirror consistent rate = %f"%(i, mirror_consistent_rate)
			print  "Group %d ...... Average pixel error = %f"%(i, pix_err)
			print  "Group %d ...... Size of stable subset = %d"%(i, l_stable_members)
			print  "Group %d ......"%i,

			if l_stable_members <= thld_grp:
				print "Size of stable subset smaller than the threshold, discarded\n"
				continue
			print "Size of stable subset larger than the threshold, kept\n"

			ave = recv_EMData(node_to_run, i+70000)
			stable_members_ori = [0]*l_stable_members
			for j in xrange(l_stable_members): stable_members_ori[j] = alldata_n[stable_members[j]]
			ave.set_attr_dict({"members": stable_members_ori, "n_objects": l_stable_members})
			ave.write_image("class_averages_generation_%d.hdf"%generation, ave_num)
			mpi_send(ave_num, 1, MPI_INT, node_to_run, i+80000, MPI_COMM_WORLD)
			ave_num += 1
			members_acc.extend(stable_members_ori)

		members_acc.sort()
		for i in xrange(len(members_acc)-1): assert members_acc[i] != members_acc[i+1]

		members_unacc = [0]*(ndata-len(members_acc))
		c = 0
		for i in xrange(ndata):
			if alldata_n[i] in members_acc: continue
			members_unacc[c] = alldata_n[i]
			c += 1

		print "In the second phase, we found %d stable and reproducible averages that account for %d particles.  "%(ave_num, len(members_acc))
		#  The following will write a zero-length file if the list is empty
		write_text_file(members_acc, "generation_%d_accounted.txt"%generation)
		write_text_file(members_unacc, "generation_%d_unaccounted.txt"%generation)
		print "******************************************************************************************"
		print "*     End of the second phase             "+strftime("%a, %d %b %Y %H:%M:%S", localtime())+"            *"
		print "******************************************************************************************"
	else:
		for i in xrange(l_STB):
			node_to_run = i%(number_of_proc-1)+1
			if myid != node_to_run: continue
			l_STB_PART = mpi_recv(1, MPI_INT, main_node, i+10000, MPI_COMM_WORLD)
			l_STB_PART = int(l_STB_PART[0])
			STB_PART = mpi_recv(l_STB_PART, MPI_INT, main_node, i+20000, MPI_COMM_WORLD)
			STB_PART = map(int, STB_PART)
			STB_PART.sort()

			class_data = [None]*l_STB_PART
			members_id = [0]*l_STB_PART
			for im in xrange(l_STB_PART):
				class_data[im] = alldata[STB_PART[im]]
				members_id[im] = alldata[STB_PART[im]].get_attr('ID')
			for im in xrange(l_STB_PART-1):
				assert members_id[im] != members_id[im+1]

			ali_params = [[] for j in xrange(stab_ali)]
			for ii in xrange(stab_ali):
				ave = within_group_refinement(class_data, None, True, ir, ou, rs, [xr], [yr], [ts], \
												dst, maxit, FH, FF, method = alimethod)
				for im in xrange(l_STB_PART):
					alpha, sx, sy, mirror, scale = get_params2D(class_data[im])
					ali_params[ii].extend([alpha, sx, sy, mirror])
			if ou == -1:  ou = nx/2-2
			stable_set, mirror_consistent_rate, pix_err = multi_align_stability(ali_params, 0.0, 10000.0, thld_err, False, ou*2)

			l_stable_set = len(stable_set)
			stable_set_id = [0]*l_stable_set
			all_alpha     = [0]*l_stable_set
			all_sx        = [0]*l_stable_set
			all_sy        = [0]*l_stable_set
			all_mirror    = [0]*l_stable_set
			#all_scale     = [1.0]*l_stable_set
			for j in xrange(l_stable_set): 
				stable_set_id[j] = members_id[stable_set[j][1]]
				all_alpha[j]     = stable_set[j][2][0]
				all_sx[j]        = stable_set[j][2][1]
				all_sy[j]        = stable_set[j][2][2]
				all_mirror[j]    = stable_set[j][2][3]

			mpi_send(l_stable_set, 1, MPI_INT, main_node, i+30000, MPI_COMM_WORLD)
			mpi_send(stable_set_id, l_stable_set, MPI_INT, main_node, i+40000, MPI_COMM_WORLD)
			mpi_send(mirror_consistent_rate, 1, MPI_FLOAT, main_node, i+50000, MPI_COMM_WORLD)
			mpi_send(pix_err, 1, MPI_FLOAT, main_node, i+60000, MPI_COMM_WORLD)

			if l_stable_set > thld_grp:
				send_EMData(ave, main_node, i+70000)		
				ave_num = mpi_recv(1, MPI_INT, main_node, i+80000, MPI_COMM_WORLD)
				ave_num = int(ave_num[0])
				write_text_file([all_alpha, all_sx, all_sy, all_mirror], "%s/ali_params_%03d"%(ali_params_dir, ave_num))

	mpi_barrier(MPI_COMM_WORLD)
	
	if myid == main_node:
		print "****************************************************************************************************"
		print "*                                                                                                  *"
		print "*                   End of the ISAC program                 "+strftime("%a, %d %b %Y %H:%M:%S", localtime())+"            *"
		print "*                                                                                                  *"
		print "****************************************************************************************************"
	return

		
# stack - list of images (filename also accepted)
# refim - list of reference images (filename also accepted)
# maskfile - image with mask (filename also accepted)
# CTF - not supported
# snr - not supported
# stability - when True stability checking is performed
# stab_ali - used only when stability=True, 
# iter_reali - used only when stability=True - for each iteration with index holds (index of iteration % iter_reali == 0) stability checking is performed
def isac_MPI(stack, refim, maskfile = None, outname = "avim", ir=1, ou=-1, rs=1, xrng=0, yrng=0, step=1, 
			 maxit=30, isac_iter=10, CTF=False, snr=1.0, rand_seed=-1, color=0, comm=-1, 
			 stability=False, stab_ali=5, iter_reali=1, thld_err=1.732, FL=0.1, FH=0.3, FF=0.2, dst=90.0, method = ""):
	
	from global_def   import EMData, Util
	from alignment	  import Numrinit, ringwe, search_range
	from applications import MPI_start_end, within_group_refinement
	from filter	      import filt_tanl
	from fundamentals import rot_shift2D, fshift, fft
	from pixel_error  import multi_align_stability
	from statistics   import ave_series
	from utilities	  import model_circle, model_blank, combine_params2, inverse_transform2, get_image
	from utilities	  import reduce_EMData_to_root, bcast_EMData_to_all
	from utilities	  import get_params2D, set_params2D
	from random	      import seed, randint, jumpahead
	from mpi		  import mpi_comm_size, mpi_comm_rank, MPI_COMM_WORLD
	from mpi		  import mpi_reduce, mpi_bcast, mpi_barrier, mpi_recv, mpi_send
	from mpi		  import MPI_SUM, MPI_FLOAT, MPI_INT, MPI_TAG_UB
	from numpy        import zeros, float32
	from time         import localtime, strftime
	import os
	import sys

	if comm == -1: comm = MPI_COMM_WORLD		

	number_of_proc = mpi_comm_size(comm)
	myid = mpi_comm_rank(comm)
	my_abs_id = mpi_comm_rank(MPI_COMM_WORLD)
	main_node = 0

	first_ring=int(ir); last_ring=int(ou); rstep=int(rs); max_iter=int(isac_iter)

	if type(stack) == type(""):
		#read all data
		alldata = EMData.read_images(stack)
	else:
		alldata = stack
	nx = alldata[0].get_xsize()

	nima = len(alldata)
	#  Explicitly force all parameters to be zero on input
	for im in xrange(nima):  set_params2D(alldata[im], [0.,0.,0.,0, 1.0])
		
	
	image_start, image_end = MPI_start_end(nima, number_of_proc, myid)

	if maskfile:
		import  types
		if type(maskfile) is types.StringType:  mask = get_image(maskfile)
		else: mask = maskfile
	else : mask = model_circle(last_ring, nx, nx)
	if type(refim) == type(""):
		refi = EMData.read_images(refim)
	else:
		# It's safer to make a hard copy here. Although I am not sure, I believe a shallow copy
		# has messed up the program.
		#   This is really strange.  It takes much memory without any need.  PAP 01/17/2015
		#      However, later I made changes so refi is deleted from time to time.  All to be checked.
		# refi = refim
		refi = [None for i in xrange(len(refim))]
		for i in xrange(len(refim)):  refi[i] = refim[i].copy()
	numref = len(refi)

	#  CTF stuff
#	if CTF:
#		ctf_params = ima.get_attr("ctf")
#		data_had_ctf = ima.get_attr("ctf_applied")
#		ctm = ctf_2(nx, ctf_params)
#		lctf = len(ctm)

	# IMAGES ARE SQUARES! center is in SPIDER convention
	cnx = nx/2+1
	cny = cnx

	mode = "F"
	#precalculate rings
	numr = Numrinit(first_ring, last_ring, rstep, mode)
	wr = ringwe(numr, mode)
	# reference images
	#  for each node read its share of data
	#data = EMData.read_images(stack, range(image_start, image_end))
#	for im in xrange(image_start, image_end):
#		#data[im-image_start].set_attr('ID', im)
#		if CTF:
#			ctf_params = alldata[im].get_attr( "ctf" )
#			if alldata[im].get_attr("ctf_applied") == 0:
#				st = Util.infomask(alldata[im], mask, False)
#				alldata[im] -= st[0]
#				from filter import filt_ctf
#				alldata[im] = filt_ctf(alldata[im], ctf_params)
#				alldata[im].set_attr('ctf_applied', 1)

	if rand_seed > -1:      seed(rand_seed)
	else:                   seed(randint(1,2000111222))
	if myid != main_node:   jumpahead(17*myid + 12345)

	fl = FL
	Iter = -1
	main_iter = 0

	while main_iter < max_iter:
		Iter += 1
		###if my_abs_id == main_node: print "Iteration within isac_MPI = ", Iter, "	main_iter = ", main_iter, "	len data = ", image_end-image_start, localtime()[0:5], myid
		mashi = cnx-ou-2
		for j in xrange(numref):
			refi[j].process_inplace("normalize.mask", {"mask":mask, "no_sigma":1}) # normalize reference images to N(0,1)
			###if myid == main_node:
			###	refi[j].write_image("refincoming%02d_round%02d.hdf"%(color, Iter), j)
			cimage = Util.Polar2Dm(refi[j] , cnx, cny, numr, mode)
			Util.Frngs(cimage, numr)
			Util.Applyws(cimage, numr, wr)
			refi[j] = cimage.copy()


#		if CTF: ctf2 = [[[0.0]*lctf for k in xrange(2)] for j in xrange(numref)]
		peak_list = [zeros(4*(image_end-image_start), dtype=float32) for i in xrange(numref)]
		#  nima is the total number of images, not the one on this node, the latter is (image_end-image_start)
		#    d matrix required by EQ-Kmeans can be huge!!  PAP 01/17/2015
		d = zeros(numref*nima, dtype=float32)
		# begin MPI section
		for im in xrange(image_start, image_end):
			alpha, sx, sy, mirror, scale = get_params2D(alldata[im])
			##  TEST WHETHER PARAMETERS ARE WITHIN RANGE
			alphai, sxi, syi, scalei = inverse_transform2(alpha, sx, sy)
			# If shifts are outside of the permissible range, reset them
			if(abs(sxi)>mashi or abs(syi)>mashi):
				sxi = 0.0
				syi = 0.0
				set_params2D(alldata[im],[0.0,0.0,0.0,0,1.0])
			# normalize
			alldata[im].process_inplace("normalize.mask", {"mask":mask, "no_sigma":0}) # subtract average under the mask
			ny = nx
			txrng = search_range(nx, ou, sxi, xrng, "ISAC")
			txrng = [txrng[1],txrng[0]]
			tyrng = search_range(ny, ou, syi, yrng, "ISAC")
			tyrng = [tyrng[1],tyrng[0]]

			# align current image to the reference
			temp = Util.multiref_polar_ali_2d_peaklist(alldata[im], refi, txrng, tyrng, step, mode, numr, cnx+sxi, cny+syi)
			for iref in xrange(numref):
				[alphan, sxn, syn, mn] = \
				   combine_params2(0.0, -sxi, -syi, 0, temp[iref*5+1], temp[iref*5+2], temp[iref*5+3], int(temp[iref*5+4]))
				peak_list[iref][(im-image_start)*4+0] = alphan
				peak_list[iref][(im-image_start)*4+1] = sxn
				peak_list[iref][(im-image_start)*4+2] = syn
				peak_list[iref][(im-image_start)*4+3] = mn
				d[iref*nima+im] = temp[iref*5]

		# ???  This is attempt to do mref with restricted searches.  It does not work out as some classes may require
		#      much larger shifts to center averages than other.

		"""
		mashi = cnx-ou-2  # needed for maximum shift
		for im in xrange(image_start, image_end):
			alpha, sx, sy, mirror, scale = get_params2D(alldata[im])
			alphai, sxi, syi, scalei     = inverse_transform2(alpha, sx, sy)
			# normalize
			alldata[im].process_inplace("normalize.mask", {"mask":mask, "no_sigma":0}) # subtract average under the mask

			ny = nx
			#  The search range procedure was adjusted for 3D searches, so since in 2D the order of operations is inverted, we have to invert ranges
			txrng = search_range(nx, ou, sxi, xrng, "ISAC")
			txrng = [txrng[1],txrng[0]]
			tyrng = search_range(ny, ou, syi, yrng, "ISAC")
			tyrng = [tyrng[1],tyrng[0]]

			# align current image to all references - THIS IS REALLY TIME CONSUMING PAP 01/17/2015
			temp = Util.multiref_polar_ali_2d_peaklist(alldata[im], refi, txrng, tyrng, step, mode, numr, cnx+sxi, cny+syi)
			for iref in xrange(numref):
				from utilities import inverse_transform2
				[alphan, sxn, syn, mn] = \
				   combine_params2(0.0, -sxi, -syi, 0, temp[iref*5+1], temp[iref*5+2], temp[iref*5+3], int(temp[iref*5+4]))
				alphan, sxn, syn, mn = inverse_transform2(alphan, sxn, syn, mn)
				sxn = min(max(round(sxn,2),-mashi),mashi)
				syn = min(max(round(syn,2),-mashi),mashi)
				alphan, sxn, syn, mn = inverse_transform2(alphan, sxn, syn, mn)
				peak_list[iref][(im-image_start)*4+0] = alphan
				peak_list[iref][(im-image_start)*4+1] = sxn
				peak_list[iref][(im-image_start)*4+2] = syn
				peak_list[iref][(im-image_start)*4+3] = mn
				qd0,qd1,qd2,qd3 = inverse_transform2(alphan, sxn, syn, mn)
				if(abs(qd1)>mashi or abs(qd2)>mashi):  print  " multiref2 ",sxi,syi,temp[iref*5+1], temp[iref*5+2], temp[iref*5+3], int(temp[iref*5+4]),alphan, sxn, syn, mn,qd0,qd1,qd2,qd3
				d[iref*nima+im] = temp[iref*5]
		"""
		"""
		#  This version does cyclic shifts of images to center them prior to multiref 
		#      to keep them within permissible range of translations.
		for im in xrange(image_start, image_end):
			alphai, sxi, syi, mirrori, scale = get_params2D(alldata[im])
			lx = int(round(sxi,0))
			ly = int(round(syi,0))
			sxi -= lx
			syi -= ly
			tempdata = alldata[im].copy()
			Util.cyclicshift(tempdata, {"dx":-lx,"dy":-ly})

			# normalize
			tempdata.process_inplace("normalize.mask", {"mask":mask, "no_sigma":0}) # subtract average under the mask

			#  Since shifts are now a fraction of pixel, we do not have to worry about checking the ranges
			# align current image to all references
			temp = Util.multiref_polar_ali_2d_peaklist(tempdata, refi, [xrng,xrng], [yrng,yrng], step, mode, numr, cnx+sxi, cny+syi)
			for iref in xrange(numref):
				alphan, sxn, syn, mn = inverse_transform2(-temp[iref*5+1], -temp[iref*5+2]+sxi+lx, -temp[iref*5+3]+syi+ly, 0)
				mn = int(temp[iref*5+4])
				peak_list[iref][(im-image_start)*4+0] = alphan
				peak_list[iref][(im-image_start)*4+1] = sxn
				peak_list[iref][(im-image_start)*4+2] = syn
				peak_list[iref][(im-image_start)*4+3] = mn
				d[iref*nima+im] = temp[iref*5]
		"""

		del refi, temp

		d = mpi_reduce(d, numref*nima, MPI_FLOAT, MPI_SUM, main_node, comm)  #  RETURNS numpy array
		if myid != main_node:
		 	del d
		mpi_barrier(comm) # to make sure that slaves freed the matrix d

		if myid == main_node:
			#  PAP 03/20/2015  added cleaning of long lists...
			id_list_long = Util.assign_groups(str(d.__array_interface__['data'][0]), numref, nima) # string with memory address is passed as parameters
			del d
			id_list = [[] for i in xrange(numref)]
			maxasi = nima/numref
			for i in xrange(maxasi*numref):
				id_list[i/maxasi].append(id_list_long[i])
			for i in xrange(nima%maxasi):
				id_list[id_list_long[-1]].append(id_list_long[maxasi*numref+i])
			for iref in xrange(numref):
				id_list[iref].sort()
			del id_list_long

			belongsto = [0]*nima
			for iref in xrange(numref):
				for im in id_list[iref]: belongsto[im] = iref
		else:
			belongsto = [0]*nima
		mpi_barrier(comm)
		belongsto = mpi_bcast(belongsto, nima, MPI_INT, main_node, comm)
		belongsto = map(int, belongsto)
		###if my_abs_id == main_node: print "Completed EQ-mref within isac_MPI = ", Iter, "	main_iter = ", main_iter , localtime()[0:5], color, myid


		#  Compute partial averages
		members = [0]*numref
		sx_sum = [0.0]*numref
		sy_sum = [0.0]*numref
		refi = [model_blank(nx,ny) for j in xrange(numref)]
		for im in xrange(image_start, image_end):
			matchref = belongsto[im]
			alphan = float(peak_list[matchref][(im-image_start)*4+0])
			sxn = float(peak_list[matchref][(im-image_start)*4+1])
			syn = float(peak_list[matchref][(im-image_start)*4+2])
			mn = int(peak_list[matchref][(im-image_start)*4+3])
			if mn == 0: sx_sum[matchref] += sxn
			else:	   sx_sum[matchref] -= sxn
			sy_sum[matchref] += syn
			# apply current parameters and add to the average
			Util.add_img(refi[matchref], rot_shift2D(alldata[im], alphan, sxn, syn, mn))
#			if CTF:
#				ctm = ctf_2(nx, ctf_params)
#				for i in xrange(lctf):  ctf2[matchref][it][i] += ctm[i]
			members[matchref] += 1
		sx_sum = mpi_reduce(sx_sum, numref, MPI_FLOAT, MPI_SUM, main_node, comm)
		sy_sum = mpi_reduce(sy_sum, numref, MPI_FLOAT, MPI_SUM, main_node, comm)
		members = mpi_reduce(members, numref, MPI_INT, MPI_SUM, main_node, comm)
		if myid != main_node:
			sx_sum = [0.0]*numref
			sy_sum = [0.0]*numref
			members = [0.0]*numref
		sx_sum = mpi_bcast(sx_sum, numref, MPI_FLOAT, main_node, comm)
		sy_sum = mpi_bcast(sy_sum, numref, MPI_FLOAT, main_node, comm)
		members = mpi_bcast(members, numref, MPI_INT, main_node, comm)
		sx_sum = map(float, sx_sum)
		sy_sum = map(float, sy_sum)
		members = map(int, members)

		for j in xrange(numref):
			sx_sum[j] /= float(members[j])
			sy_sum[j] /= float(members[j])

		for im in xrange(image_start, image_end):
			matchref = belongsto[im]
			alphan = float(peak_list[matchref][(im-image_start)*4+0])
			sxn = float(peak_list[matchref][(im-image_start)*4+1])
			syn = float(peak_list[matchref][(im-image_start)*4+2])
			mn = int(peak_list[matchref][(im-image_start)*4+3])
			if mn == 0:
				set_params2D(alldata[im], [alphan, sxn-sx_sum[matchref], syn-sy_sum[matchref], mn, scale])
			else:
				set_params2D(alldata[im], [alphan, sxn+sx_sum[matchref], syn-sy_sum[matchref], mn, scale])

		del peak_list

		for j in xrange(numref):
			reduce_EMData_to_root(refi[j], myid, main_node, comm)
			if myid == main_node:
				# Golden rule when to do within group refinement
				Util.mul_scalar(refi[j], 1.0/float(members[j]))
				refi[j] = filt_tanl(refi[j], fl, FF)
				refi[j] = fshift(refi[j], -sx_sum[j], -sy_sum[j])
				set_params2D(refi[j], [0.0, 0.0, 0.0, 0, 1.0])

		if myid == main_node:
			#  this is most likely meant to center them, if so, it works poorly, 
			#      it has to be checked and probably a better method used PAP 01/17/2015
			dummy = within_group_refinement(refi, mask, True, first_ring, last_ring, rstep, [xrng], [yrng], [step], dst, maxit, FH, FF)
			ref_ali_params = []
			for j in xrange(numref):
				alpha, sx, sy, mirror, scale = get_params2D(refi[j])
				refi[j] = rot_shift2D(refi[j], alpha, sx, sy, mirror)
				ref_ali_params.extend([alpha, sx, sy, mirror])
		else:
			ref_ali_params = [0.0]*(numref*4)
		ref_ali_params = mpi_bcast(ref_ali_params, numref*4, MPI_FLOAT, main_node, comm)
		ref_ali_params = map(float, ref_ali_params)

		for j in xrange(numref):
			bcast_EMData_to_all(refi[j], myid, main_node, comm)

		for j in xrange(numref):
			bcast_EMData_to_all(refi[j], myid, main_node, comm)

		###if myid == main_node:
		###	print  "  WRITING refaligned  for color:",color
		###	for j in xrange(numref):
		###		refi[j].write_image("refaligned%02d_round%02d.hdf"%(color, Iter), j)

		# Compensate the centering to averages
		for im in xrange(image_start, image_end):
			matchref = belongsto[im]
			alpha, sx, sy, mirror, scale = get_params2D(alldata[im])
			alphan, sxn, syn, mirrorn = combine_params2(alpha, sx, sy, mirror, ref_ali_params[matchref*4], ref_ali_params[matchref*4+1], \
				ref_ali_params[matchref*4+2], int(ref_ali_params[matchref*4+3]))
			set_params2D(alldata[im], [alphan, sxn, syn, int(mirrorn), 1.0])

		do_within_group = 0
		fl += 0.05
		if fl >= FH:
			fl = FL
			do_within_group = 1

		# Here stability does not need to be checked for each main iteration, it only needs to
		# be done for every 'iter_reali' iterations. If one really wants it to be checked each time
		# simple set iter_reali to 1, which is the default value right now.
		check_stability = (stability and (main_iter%iter_reali==0))

		if do_within_group == 1:
			###if my_abs_id == main_node: print "Doing within group alignment .......", localtime()[0:5]

			# Broadcast the alignment parameters to all nodes
			for i in xrange(number_of_proc):
				im_start, im_end = MPI_start_end(nima, number_of_proc, i)
				if myid == i:
					ali_params = []
					for im in xrange(image_start, image_end):
						alpha, sx, sy, mirror, scale = get_params2D(alldata[im])
						ali_params.extend([alpha, sx, sy, mirror])
				else:
					ali_params = [0.0]*((im_end-im_start)*4)
				ali_params = mpi_bcast(ali_params, len(ali_params), MPI_FLOAT, i, comm)
				ali_params = map(float, ali_params)
				for im in xrange(im_start, im_end):
					alpha = ali_params[(im-im_start)*4]
					sx = ali_params[(im-im_start)*4+1]
					sy = ali_params[(im-im_start)*4+2]
					mirror = int(ali_params[(im-im_start)*4+3])
					set_params2D(alldata[im], [alpha, sx, sy, mirror, 1.0])

			main_iter += 1

			# There are two approaches to scatter calculations among MPI processes during stability checking.
			# The first one is the original one. I added the second method.
			# Here we try to estimate the calculation time for both approaches.
			stab_calc_time_method_1 = stab_ali * ((numref-1) // number_of_proc + 1)
			stab_calc_time_method_2 = (numref * stab_ali - 1) // number_of_proc + 1
			#if my_abs_id == main_node: print "Times estimation: ", stab_calc_time_method_1, stab_calc_time_method_2

			# When there is no stability checking or estimated calculation time of new method is greater than 80% of estimated calculation time of original method 
			# then the original method is used. In other case. the second (new) method is used.
			#if (not check_stability) or (stab_calc_time_method_2 > 0.80 * stab_calc_time_method_1):
			#  For the time being only use this method as the other one is not worked out as far as parameter ranges go.
			if True :
				###if my_abs_id == main_node: print "Within group refinement and checking within group stability, original approach .......", check_stability, "  ",localtime()[0:5]
				# ====================================== standard approach is used, calculations are parallelized by scatter groups (averages) among MPI processes
				if check_stability and main_iter == max_iter-1: gpixer = []
				for j in xrange(myid, numref, number_of_proc):
					assign = []
					for im in xrange(nima):
						if j == belongsto[im]:  assign.append(im)

					randomize = True  # I think there is no reason not to be True
					class_data = [alldata[im] for im in assign]
					refi[j] = within_group_refinement(class_data, mask, randomize, first_ring, last_ring, rstep, \
													[xrng], [yrng], [step], dst, maxit, FH, FF, method = method)

					if check_stability:
						###if my_abs_id == main_node: print "Checking within group stability, original approach .......", check_stability, "  ",localtime()[0:5]
						ali_params = [[] for qq in xrange(stab_ali)]
						for ii in xrange(stab_ali):
							if ii > 0:  # The first one does not have to be repeated
								dummy = within_group_refinement(class_data, mask, randomize, first_ring, last_ring, rstep, [xrng], [yrng], [step], \
																dst, maxit, FH, FF, method = method)
							for im in xrange(len(class_data)):
								alpha, sx, sy, mirror, scale = get_params2D(class_data[im])
								ali_params[ii].extend([alpha, sx, sy, mirror])

						stable_set, mirror_consistent_rate, err = multi_align_stability(ali_params, 0.0, 10000.0, thld_err, False, last_ring*2)
						if main_iter == max_iter-1 :  gpixer.append(err)

						###print  "Color %1d, class %4d ...... Size of the group = %4d and of the stable subset = %4d  Mirror consistent rate = %5.3f  Average pixel error prior to class pruning = %10.2f"\
						###				%(color, j, len(class_data), len(stable_set), mirror_consistent_rate, err)

						# If the size of stable subset is too small (say 1, 2), it will cause many problems, so we manually increase it to 5
						while len(stable_set) < 5:
							duplicate = True
							while duplicate:
								duplicate = False
								p = randint(0, len(class_data)-1)
								for ss in stable_set:
									if p == ss[1]: duplicate = True
							stable_set.append([100.0, p, [0.0, 0.0, 0.0, 0]])
						stable_data = []
						stable_members = []
						for err in stable_set:
							im = err[1]
							stable_members.append(assign[im])
							stable_data.append(class_data[im])
							set_params2D( class_data[im], [err[2][0], err[2][1], err[2][2], int(err[2][3]), 1.0] )
						stable_members.sort()

						refi[j] = filt_tanl(ave_series(stable_data), FH, FF)
						refi[j].set_attr('members', stable_members)
						refi[j].set_attr('n_objects', len(stable_members))
						del stable_members
					# end of stability
					del assign
				if check_stability and main_iter == max_iter-1:
					#  gather all pixers and print a histogram
					from utilities import wrap_mpi_gatherv
					gpixer = wrap_mpi_gatherv(gpixer, main_node, comm)
					if my_abs_id == main_node and color == 0:
						from statistics   import hist_list
						lhist = 12
						region, histo = hist_list(gpixer, lhist)
						print  "\n=== Histogram of average within-class pixel errors prior to class pruning ==="
						for lhx in xrange(lhist):  print   "     %10.3f     %7d"%(region[lhx], histo[lhx])
						print  "=============================================================================\n"
					del gpixer
				mpi_barrier(comm)

				for im in xrange(nima):
					done_on_node = belongsto[im]%number_of_proc
					if myid == done_on_node:
						alpha, sx, sy, mirror, scale = get_params2D(alldata[im])
						ali_params = [alpha, sx, sy, mirror]
					else:
						ali_params = [0.0]*4
					ali_params = mpi_bcast(ali_params, 4, MPI_FLOAT, done_on_node, comm)
					ali_params = map(float, ali_params)
					set_params2D(alldata[im], [ali_params[0], ali_params[1], ali_params[2], int(ali_params[3]), 1.0])

			else:
				###if my_abs_id == main_node: print "Checking within group stability, new approach .......", localtime()[0:5]
				# ================================================ more complicated approach is used - runs of within_group_refinement are scattered among MPI processes
				refi = isac_stability_check_mpi(alldata, numref, belongsto, stab_ali, thld_err, mask, first_ring, last_ring, rstep, xrng, yrng, step, \
												dst, maxit, FH, FF, method, comm)

			for j in xrange(numref):
				bcast_EMData_to_all(refi[j], myid, j%number_of_proc, comm)

			if check_stability:
				# In this case, we need to set the 'members' attr using stable members from the stability test
				for j in xrange(numref):
					done_on_node = j%number_of_proc
					if done_on_node != main_node:
						if myid == main_node:
							mem_len = mpi_recv(1, MPI_INT, done_on_node, MPI_TAG_UB, comm)
							mem_len = int(mem_len[0])
							members = mpi_recv(mem_len, MPI_INT, done_on_node, MPI_TAG_UB, comm)
							members = map(int, members)
							refi[j].set_attr_dict({'members': members,'n_objects': mem_len})
						elif myid == done_on_node:
							members = refi[j].get_attr('members')
							mpi_send(len(members), 1, MPI_INT, main_node, MPI_TAG_UB, comm)
							mpi_send(members, len(members), MPI_INT, main_node, MPI_TAG_UB, comm)
			###if myid == main_node:  print "within group alignment done. ", localtime()[0:5]
			###if myid == main_node:
			###	print  "  WRITING refrealigned  for color:",color
			###	for j in xrange(numref):
			###		refi[j].write_image("refrealigned%02d_round%02d.hdf"%(color, Iter), j)

		# end of do_within_group
		mpi_barrier(comm)

		if myid == main_node:
			#  I added a switch here.  I do not think we will need those in the future.  PAP 03/26
			if outname != None:
				final_outname = outname+'%02d_%03d.hdf'%(color, Iter)
				if os.path.exists(final_outname):
					os.remove(final_outname)
			if check_stability:
				# In this case, the attr 'members' is defined as the stable members, its setting is done
				# in the code before
				if outname != None:
					for j in xrange(numref):
						refi[j].write_image(final_outname, j)
			else:
				for j in xrange(numref):
					refi[j].set_attr_dict({'members': id_list[j], 'n_objects': len(id_list[j])})
					if outname != None:
						refi[j].write_image(final_outname, j)
			del id_list
		mpi_barrier(comm)

	return refi


#  I blocked usage of this procedure.  Yang called it a "new" approach, but the standard one seems to be timewise good enough. 07/09/2015 PAP
# This routine runs within_group_refinement several times in parallel for multiple groups (parameter randomize must be set to True)
# All MPI processes must have the same values of all parameters
# This function returns list of references images with numref elements, elements with index holds (index % mpi_comm_size(comm) == mpi_comm_rank(comm)) contains corresponding reference images, 
# rest of elements contains blank images
'''
def isac_stability_check_mpi(alldata, numref, belongsto, stab_ali, thld_err, mask, first_ring, last_ring, rstep, xrng, yrng, step, \
								dst, maxit, FH, FF, alimethod, comm):
	from applications import within_group_refinement
	from mpi		  import mpi_comm_size, mpi_comm_rank, mpi_barrier, mpi_bcast, mpi_send, mpi_recv, MPI_FLOAT, MPI_TAG_UB
	from pixel_error  import multi_align_stability
	from utilities	import get_params2D, set_params2D, model_blank
	from filter	   import filt_tanl
	from random	   import randint
	from statistics   import ave_series
	from time         import localtime, strftime
	
	myid = mpi_comm_rank(comm)
	number_of_proc = mpi_comm_size(comm)
	nima = len(alldata)
	refi = []
	for i in xrange(numref):
		refi.append(model_blank(alldata[0].get_xsize(), alldata[0].get_ysize()))

	# --- divide images into groups
	grp_images = []
	grp_indexes = []
	for i in xrange(numref):
		grp_images.append([])
		grp_indexes.append([])
	for i in xrange(nima):
		grp_images[belongsto[i]].append(alldata[i])
		grp_indexes[belongsto[i]].append(i)
	
	# --- prepare mapping (group id, run index) -> process id
	grp_run_to_mpi_id = []
	for ig in xrange(numref):
		grp_run_to_mpi_id.append([])
		for ii in xrange(stab_ali):
			grp_run_to_mpi_id[ig].append( (ig * stab_ali + ii) % number_of_proc )
	
	# --- prepare structure for alignment parameters: (group id, run index) -> []
	grp_run_to_ali_params = []
	for ig in xrange(numref):
		grp_run_to_ali_params.append([])
		for ii in xrange(stab_ali):
			grp_run_to_ali_params[ig].append([])
	
	# --- run within_group_refinement
#	print "run within_group_refinement"
	for ig in xrange(numref):
		for ii in xrange(stab_ali):
			if grp_run_to_mpi_id[ig][ii] == myid:
				within_group_refinement(grp_images[ig], mask, True, first_ring, last_ring, rstep, [xrng], [yrng], [step], \
											dst, maxit, FH, FF, method = alimethod)
				for im in (grp_images[ig]):
					alpha, sx, sy, mirror, scale = get_params2D(im)
					grp_run_to_ali_params[ig][ii].extend([alpha, sx, sy, mirror])

	# --- send obtained alignment parameters to target MPI process, alignment parameters from last runs are broadcasted and copied to images' headers
#	print "send obtained alignment parameters to target MPI process"
	for ig in xrange(numref):
		grp_size = len(grp_images[ig])
		for ii in xrange(stab_ali-1):
			src_proc_id = grp_run_to_mpi_id[ig][ii]
			trg_proc_id = ig % number_of_proc
			if src_proc_id == trg_proc_id:
				continue
			if myid == src_proc_id:
				mpi_send(grp_run_to_ali_params[ig][ii], 4*grp_size, MPI_FLOAT, trg_proc_id, MPI_TAG_UB, comm)
			if myid == trg_proc_id:
				grp_run_to_ali_params[ig][ii] = mpi_recv(4*grp_size, MPI_FLOAT, src_proc_id, MPI_TAG_UB, comm)
				grp_run_to_ali_params[ig][ii] = map(float, grp_run_to_ali_params[ig][ii])

	for ig in xrange(numref):
		grp_size = len(grp_images[ig])
		ii = stab_ali - 1
		src_proc_id = grp_run_to_mpi_id[ig][ii]
		if myid != src_proc_id:
			grp_run_to_ali_params[ig][ii] = [0.0] * (4 * grp_size)
		grp_run_to_ali_params[ig][ii] = mpi_bcast(grp_run_to_ali_params[ig][ii], 4*grp_size, MPI_FLOAT, src_proc_id, comm)
		grp_run_to_ali_params[ig][ii] = map(float, grp_run_to_ali_params[ig][ii])
		for i in xrange(len(grp_images[ig])):
			set_params2D(grp_images[ig][i], [grp_run_to_ali_params[ig][ii][4*i+0], grp_run_to_ali_params[ig][ii][4*i+1], grp_run_to_ali_params[ig][ii][4*i+2], int(grp_run_to_ali_params[ig][ii][4*i+3]), 1.0])

	# ======================= rest of stability checking is analogical to the old approach
#	print "rest of code..."
	for j in xrange(myid, numref, number_of_proc):
		
		stable_set, mirror_consistent_rate, err = multi_align_stability(grp_run_to_ali_params[j], 0.0, 10000.0, thld_err, False, last_ring*2)
		#print  "Stability check, class %d ...... Size of the group = %d and of the stable subset = %d, Pixer threshold = %f, Mirror consistent rate = %f,  Average pixel error = %f"\
		#			%(j, len(class_data), len(stable_set),thld_err, mirror_consistent_rate, err)

		# If the size of stable subset is too small (say 1, 2), it will cause many problems, so we manually increase it to 5
		while len(stable_set) < 5:
			duplicate = True
			while duplicate:
				duplicate = False
				p = randint(0, len(grp_images[j])-1)
				for ss in stable_set:
					if p == ss[1]: duplicate = True
			stable_set.append([100.0, p, [0.0, 0.0, 0.0, 0]])

		stable_data = []
		stable_members = []
		for err in stable_set:
			im = err[1]
			stable_members.append(grp_indexes[j][im])
			stable_data.append(grp_images[j][im])
			set_params2D( grp_images[j][im], [err[2][0], err[2][1], err[2][2], int(err[2][3]), 1.0] )
		stable_members.sort()

		refi[j] = filt_tanl(ave_series(stable_data), FH, FF)
		refi[j].set_attr('members', stable_members)
		refi[j].set_attr('n_objects', len(stable_members))
		del stable_members
		# end of stability
	
	mpi_barrier(comm)
	return refi
'''


def match_independent_runs(data, refi, n_group, T):

	from numpy	     import array
	from statistics  import k_means_stab_bbenum

	K = len(refi)/n_group
	Parts = []
	for i in xrange(n_group):
		part = []
		for k in xrange(K):
			lid = refi[i*K+k].get_attr('members')
			lid = array(lid, 'int32')
			lid.sort()
			part.append(lid)
		Parts.append(part)

	#print Parts
	#print "Before matching = ", localtime()[:5]
	MATCH, STB_PART, CT_s, CT_t, ST, st = k_means_stab_bbenum(Parts, T=T, J=50, max_branching=40, stmult=0.1, branchfunc=2)
	#print "After matching = ", localtime()[:5]

	# I commented out next three, not much use printing them,  PAP.
	#print MATCH
	#print STB_PART
	#print CT_s
	#print CT_t
	#print ST
	#print st

	cost_by_match_thresh = []
	for i in xrange(len(CT_s)):
		if CT_s[i] > T:
			cost_by_match_thresh.append(CT_s[i])

	print "%d-way match: total cost of matches over threshold: "%(len(Parts)), sum(cost_by_match_thresh)
	print "%d-way match: total number of matches over threshold: "%(len(Parts)), len(cost_by_match_thresh)
	print "%d-way match: cost by match over threshold: "%(len(Parts)), cost_by_match_thresh
	print " "

	STB_PART_cleaned = []
	for i in xrange(len(STB_PART)):
		if len(STB_PART[i]) > T:
			STB_PART_cleaned.append(STB_PART[i])

	return STB_PART_cleaned



def match_2_way(data, refi, indep_run, thld_grp, FH, FF, find_unique=True, wayness=2, suffix=""):

	from utilities  import read_text_row, set_params2D
	from statistics import ave_series, k_means_stab_bbenum
	from random	    import randint, shuffle
	from filter	    import filt_tanl
	from numpy	    import array

	K = len(refi)/indep_run
	run = range(indep_run)
	shuffle(run)

	#print run

	reproducible_avgs = []
		
	for irun in xrange(indep_run):
		filename = "ali_params_%d"%run[irun] + suffix
		all_ali_params = read_text_row(filename)
	
		Parts = []
		part = [] 
		for k in xrange(K): 
			lid = refi[run[irun]*K+k].get_attr('members') 
			lid = array(lid, 'int32') 
			lid.sort() 
			part.append(lid)
		Parts.append(part)

		part = [] 
		for k in xrange(K): 
			lid = refi[run[(irun+1)%indep_run]*K+k].get_attr('members') 
			lid = array(lid, 'int32') 
			lid.sort() 
			part.append(lid)
		Parts.append(part)

		if wayness == 3:
			part = [] 
			for k in xrange(K): 
				lid = refi[run[(irun+2)%indep_run]*K+k].get_attr('members') 
				lid = array(lid, 'int32') 
				lid.sort() 
				part.append(lid)
			Parts.append(part)

		#print Parts
		MATCH, STB_PART, CT_s, CT_t, ST, st = k_means_stab_bbenum(Parts, T=thld_grp, J=50, max_branching=40, stmult=0.1, branchfunc=2)

		cost_by_match_thresh = []
		for i in xrange(len(CT_s)):
			if CT_s[i] > thld_grp:
				cost_by_match_thresh.append(CT_s[i])

		"""  I have no idea why would anybody want to see it  03/25/2014  PAP
		if wayness == 3:     print "%d-way match  %d-%d-%d:"%(len(Parts),run[irun], run[(irun+1)%indep_run], run[(irun+2)%indep_run])
		else:                print "%d-way match  %d-%d:"%(len(Parts),run[irun], run[(irun+1)%indep_run])
		print "  total cost of matches over threshold: ", sum(cost_by_match_thresh)
		print "  total number of matches over threshold: ", len(cost_by_match_thresh)
		print "  cost by match over threshold: ", cost_by_match_thresh
		print " "
		"""
		for i in xrange(len(STB_PART)):
			if len(STB_PART[i]) > 0:
				class_data = []
				members_id = []
				for im in STB_PART[i]:
					set_params2D(data[im], [all_ali_params[im][0], all_ali_params[im][1], all_ali_params[im][2], int(all_ali_params[im][3]), 1.0])
					class_data.append(data[im])
					members_id.append(data[im].get_attr('ID'))
				ave = ave_series(class_data)
				ave.set_attr_dict({"members": members_id, "n_objects": len(members_id)})
				ave = filt_tanl(ave, FH, FF)
				reproducible_avgs.append(ave)
			else:
				# put some dummy avgs here
				temp_id = randint(0, len(data)-1)
				ave = data[temp_id].copy()
				ave.set_attr_dict({"members": [temp_id], "n_objects": 1})
				reproducible_avgs.append(ave)

	if find_unique:
		# Here the idea is like this: in all reproducible averages, find the ones that are unique
		reproducible_avgs_unique = get_unique_averages(reproducible_avgs, indep_run)
		print "Found %d unique class averages through %d-way matching"%(len(reproducible_avgs_unique), wayness)
		return reproducible_avgs_unique
	else:
		print "Found %d class averages through %d-way matching"%(len(reproducible_avgs), wayness)
		return reproducible_avgs



def generate_random_averages(data, K, rand_seed = -1):
	from random import shuffle, seed, randint	
	#  I prefer to take random images....  PAP
	if rand_seed == -1:  seed(randint(1,2000111222))
	else:                seed(rand_seed)
	ndata = len(data)
	ll = range(ndata)
	shuffle(ll)
	return [data[ll[i]].copy() for i in xrange(K)]

	'''
	from alignment import align2d
	from fundamentals import rot_shift2D
	for im in xrange(K,ndata):
		wnmr = im%K
		alpha,sx,sy,mirror,peak = align2d(data[ll[im]], avgs[wnmr], 1,1,0.5,1,30)
		Util.add_img(avgs[wnmr], rot_shift2D(data[ll[im]], alpha, sx, sy, mirror))

	l = [[] for i in xrange(K)]
	for im in xrange(ndata):
		l[randint(0, K-1)].append(im)

	avgs = []
	for k in xrange(K):
		if len(l[k]) > 1:
			temp_stack = [None]*len(l[k])
			for i in xrange(len(l[k])):
				temp_stack[i] = data[l[k][i]]
			ave = aveq(temp_stack, mode = "")
		elif len(l[k]) == 1: ave = data[l[k][0]]
		else: ave = data[randint(0, ndata)]
		avgs.append(ave)
	'''
	return avgs



def get_unique_averages(data, indep_run, m_th=0.45):
	
	size_all = len(data)
	size = size_all/indep_run
	assert size_all%indep_run == 0	

	# Meaning of flag
	# 0 - not matched yet
	# 1 - matched, kept
	# 2 - matched, discard, as we only want to keep a unique copy
	# 3 - not matched, kept
	flag = [0]*size_all
	
	for i in xrange(size_all-size):
		m1 = map(int, data[i].get_attr('members'))
		if len(m1) < 5: continue
		for j in xrange((i/size+1)*size, size_all):
			m2 = map(int, data[j].get_attr('members'))
			if len(m2) < 5: continue
			m = set(m1).intersection(set(m2))
			if float(len(m)) > m_th*(len(m1)+len(m2))/2:
				if flag[i] == 0 and flag[j] == 0:
					flag[i] = 1
					flag[j] = 2
					break
				elif flag[i] == 2 or flag[j] == 2:
					flag[i] = 2
					flag[j] = 2
				else:
					print "Impossible: Something is wrong!"
		if flag[i] == 0: flag[i] = 3

	data_good = []
	# Give priority to the averages that appeared more than once
	for im in xrange(size_all):
		if flag[im] == 1:	 data_good.append(data[im])
	for im in xrange(size_all):
		if flag[im] == 3:	 data_good.append(data[im])
	
	return data_good

"""
#  This program removes from candidate averages numbers of accounted for images
#  It seems to work but it would have to be tested should we decide to go with recycling of candidates.
from EMAN2 import *
from sparx import *
la = map(int, read_text_file('generation_1_accounted.txt'))

lu = map(int, read_text_file('generation_1_unaccounted.txt'))

nn = max(max(la), max(lu))

na = range(nn)
for i in xrange(len(la)):
	na[la[i]] = -1

l = 0
for i in xrange(nn):
	if(na[i] > -1):
		na[i] = l
		l += 1

d = EMData.read_images('class_averages_candidate_generation_1.hdf')

l=0
for i in xrange(len(d)):
	li = d[i].get_attr('members')
	ou = []
	for k in xrange(len(li)):
		m = na[li[k]]
		if(m>-1):  ou.append(m)
	if(len(ou)>0):
		d[i].set_attr('members',ou)
		d[i].write_image('class_averages_candidate_generation_2.hdf',l)
		l += 1
	else:
		print ' Group  ',i,'  skipped'
"""
