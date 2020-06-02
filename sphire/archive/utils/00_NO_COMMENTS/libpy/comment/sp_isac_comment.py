





















































































































































'''0
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





































































































































































































































































































































































































































































































































































































































































































































"""1
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
				from sp_utilities import inverse_transform2
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
"""2
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






















































































































































































































































































































'''3
def isac_stability_check_mpi(alldata, numref, belongsto, stab_ali, thld_err, mask, first_ring, last_ring, rstep, xrng, yrng, step, \
								dst, maxit, FH, FF, alimethod, comm):
	from sp_applications import within_group_refinement
	from mpi		  import mpi_comm_size, mpi_comm_rank, mpi_barrier, mpi_bcast, mpi_send, mpi_recv, MPI_FLOAT
	from sp_pixel_error  import multi_align_stability
	from sp_utilities	import get_params2D, set_params2D, model_blank
	from sp_filter	   import filt_tanl
	from random	   import randint
	from sp_statistics   import ave_series
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
				mpi_send(grp_run_to_ali_params[ig][ii], 4*grp_size, MPI_FLOAT, trg_proc_id, SPARX_MPI_TAG_UNIVERSAL, comm)
			if myid == trg_proc_id:
				grp_run_to_ali_params[ig][ii] = mpi_recv(4*grp_size, MPI_FLOAT, src_proc_id, SPARX_MPI_TAG_UNIVERSAL, comm)
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








































































































"""  I have no idea why would anybody want to see it  03/25/2014  PAP4
		if wayness == 3:     print "%d-way match  %d-%d-%d:"%(len(Parts),run[irun], run[(irun+1)%indep_run], run[(irun+2)%indep_run])
		else:                print "%d-way match  %d-%d:"%(len(Parts),run[irun], run[(irun+1)%indep_run])
		print "  total cost of matches over threshold: ", sum(cost_by_match_thresh)
		print "  total number of matches over threshold: ", len(cost_by_match_thresh)
		print "  cost by match over threshold: ", cost_by_match_thresh
		print " "
		"""








































'''5
	from sp_alignment import align2d
	from sp_fundamentals import rot_shift2D
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













































"""6
#  This program removes from candidate averages numbers of accounted for images
#  It seems to work but it would have to be tested should we decide to go with recycling of candidates.
from EMAN2 import *
from sp_sparx import *
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




