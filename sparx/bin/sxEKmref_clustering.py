#!/usr/bin/env python
#
#  08/13/2015
#  New version.  
#
#
import os
import global_def
from   global_def import *
from   optparse  import OptionParser
import sys
from   numpy     import array
import types
from   logger    import Logger, BaseLogger_Files

def mref_ali3d_MPI_beta(stack, ref_vol, outdir, maskfile=None, focus = None, maxit=1, ir=1, ou=-1, rs=1, \
            xr ="4 2  2  1", yr="-1", ts="1 1 0.5 0.25",   delta="10  6  4  4", an="-1", center = -1, \
            nassign = 3, nrefine= 1, CTF = False, snr = 1.0,  ref_a="S", sym="c1",
			user_func_name="ref_ali3d", npad = 2, debug = False, fourvar=False, termprec = 0.0,\
			mpi_comm = None, log = None,frequency_low_pass=.4):
	from utilities      import model_circle, reduce_EMData_to_root, bcast_EMData_to_all, bcast_number_to_all, drop_image
	from utilities      import bcast_string_to_all, bcast_list_to_all, get_image, get_input_from_string, get_im
	from utilities      import get_arb_params, set_arb_params, drop_spider_doc, send_attr_dict
	from utilities      import get_params_proj, set_params_proj, model_blank, wrap_mpi_bcast
	from filter         import filt_params, filt_btwl, filt_ctf, filt_table, fit_tanh, filt_tanl
	from utilities      import rotate_3D_shift,estimate_3D_center_MPI
	from alignment      import Numrinit, prepare_refrings, proj_ali_incore
	from random         import randint, random
	from filter         import filt_ctf
	from utilities      import print_begin_msg, print_end_msg, print_msg, read_text_file
	from projection     import prep_vol, prgs, project, prgq, gen_rings_ctf
	from morphology     import binarize
	import os
	import types
	from mpi            import mpi_bcast, mpi_comm_size, mpi_comm_rank, MPI_FLOAT, MPI_COMM_WORLD, mpi_barrier
	from mpi            import mpi_reduce, mpi_gatherv, mpi_scatterv, MPI_INT, MPI_SUM
        from applications import MPI_start_end

	if mpi_comm == None: mpi_comm = MPI_COMM_WORLD

	if log == None:
		from logger import Logger
		log = Logger()

	number_of_proc = mpi_comm_size(mpi_comm)
	myid           = mpi_comm_rank(mpi_comm)
	main_node = 0

	if myid == 0:
		if os.path.exists(outdir):  nx = 1
		else:  nx = 0
	else:  nx = 0
	ny = bcast_number_to_all(nx, source_node = main_node)
	
	if ny == 1:  ERROR('Output directory exists, please change the name and restart the program', "mref_ali3d_MPI", 1,myid)
	mpi_barrier(MPI_COMM_WORLD)

	if myid == main_node:	
		os.mkdir(outdir)
		import global_def
		global_def.LOGFILE =  os.path.join(outdir, global_def.LOGFILE)
		log.add("Equal Kmeans-modified K-means  ")
	mpi_barrier(MPI_COMM_WORLD)

	from time import time	

	if debug:
		from time import sleep
		while not os.path.exists(outdir):
			print  "Node ",myid,"  waiting..."
			sleep(5)

		finfo = open(os.path.join(outdir, "progress%04d"%myid), 'w')
		frec  = open( os.path.join(outdir, "recons%04d"%myid), "w" )
	else:
		finfo = None
		frec  = None

	xrng        = get_input_from_string(xr)
	if  yr == "-1":  yrng = xrng
	else          :  yrng = get_input_from_string(yr)
	step        = get_input_from_string(ts)
	delta       = get_input_from_string(delta)
	lstp = min( len(xrng), len(yrng), len(step), len(delta) )
	if (an == "-1"):
		an = []
		for i in xrange(len(xrng)):   an.append(-1)
	else:
		from  alignment	    import proj_ali_incore_local
		an      = get_input_from_string(an)

	first_ring  = int(ir)
	rstep       = int(rs)
	last_ring   = int(ou)
	center      = int(center)

	numref = EMUtil.get_image_count(ref_vol)
	volref     = EMData()
	volref.read_image(stack, 0)
	nx      = volref.get_xsize()
	if last_ring < 0:	last_ring = nx//2 - 2

	if (myid == main_node):
		import user_functions
		user_func = user_functions.factory[user_func_name]
		log.add("mref_ali3d_MPI")
		log.add("Input stack                               : %s"%(stack))
		log.add("Reference volumes                         : %s"%(ref_vol))	
		log.add("Number of reference volumes               : %i"%(numref))
		log.add("Output directory                          : %s"%(outdir))
		log.add("User function                             : %s"%(user_func_name))
		if(focus != None):  \
		log.add("Maskfile 3D for focused clustering        : %s"%(focus))
		log.add("Overall 3D mask applied in user function  : %s"%(maskfile))
		log.add("Inner radius                              : %i"%(first_ring))
		log.add("Outer radius                              : %i"%(last_ring))
		log.add("Ring step                                 : %i"%(rstep))
		log.add("X search range                            : %s"%(xrng))
		log.add("Y search range                            : %s"%(yrng))
		log.add("Translational step                        : %s"%(step))
		log.add("Angular step                              : %s"%(delta))
		log.add("Angular search range                      : %s"%(an))
		log.add("Number of assignments in each iteration   : %i"%(nassign))
		log.add("Number of alignments in each iteration    : %i"%(nrefine))
		log.add("Number of iterations                      : %i"%(lstp*maxit) )
		log.add("Center type                               : %i"%(center))
		log.add("CTF correction                            : %s"%(CTF))
		log.add("Signal-to-Noise Ratio                     : %f"%(snr))
		log.add("Reference projection method               : %s"%(ref_a))
		log.add("Symmetry group                            : %s"%(sym))
		log.add("Percentage of change for termination      : %f"%(termprec))
		log.add("User function                             : %s"%(user_func_name))

	if(maskfile):
		if(type(maskfile) is types.StringType): mask3D = get_image(maskfile)
		else: 	                                mask3D = maskfile
	else        :  mask3D = model_circle(last_ring, nx, nx, nx)

	numr     = Numrinit(first_ring, last_ring, rstep, "F")
	mask2D   = model_circle(last_ring, nx, nx)
	if(first_ring > 1):  mask2D -= model_circle(first_ring, nx, nx)


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
	'''
	if(myid == main_node):	
		total_nima = EMUtil.get_image_count(stack)
		list_of_particles = range(total_nima)
	
	else:
		total_nima =0

	total_nima = bcast_number_to_all(total_nima, source_node = main_node)

	if(myid != main_node):
		list_of_particles = [-1]*total_nima

	list_of_particles = bcast_list_to_all(list_of_particles, myid,  source_node = main_node)

	image_start, image_end = MPI_start_end(total_nima, number_of_proc, myid)
	# create a list of images for each node
	list_of_particles = list_of_particles[image_start: image_end]
	nima = len(list_of_particles)
	'''
	if debug:
		finfo.write( "Image_start, image_end: %d %d\n" %(image_start, image_end) )
		finfo.flush()

	start_time = time()
	data = [None]*nima
	#  Here the assumption is that input are always volumes.  It should be most likely be changed so optionally these are group assignments.
	#  Initialize Particle ID and set group number to non-existant -1
	for im in xrange(nima):
		if( type(stack) is types.StringType ):
			data[im] = get_im(stack, list_of_particles[im])
			data[im].set_attr_dict({'ID':list_of_particles[im], 'group':-1})
		else:
			data[im] = stack[list_of_particles[im]]
			#  NOTE: in case data comes in, it would have to have ID set as there is no way to tell here what was the original ordering.
			data[im].set_attr_dict({ 'group':-1})
	if(myid == 0):
		log.add( "Time to read data: %d" % (time()-start_time) );start_time = time()

	if fourvar:
		from reconstruction import rec3D_MPI
		from statistics     import varf3d_MPI
		#  Compute Fourier variance
		vol, fscc = rec3D_MPI(data, snr, sym, model_circle(last_ring, nx, nx, nx), os.path.join(outdir, "resolution0000"), myid, main_node, finfo=frec, npad=npad)
		varf = varf3d_MPI(data, os.path.join(outdir, "ssnr0000"), None, vol, last_ring, 1.0, 1, CTF, 1, sym, myid)
		if myid == main_node:   
			varf = 1.0/varf
			varf.write_image( os.path.join(outdir,"varf0000.hdf") )
	else:
		varf = None

	if myid == main_node:
		refdata = [None]*7
		for  iref in xrange(numref):
			vol = get_im(ref_vol, iref).write_image(os.path.join(outdir, "vol0000.hdf"), iref)
		refdata[0] = numref
		refdata[1] = outdir
		refdata[2] = frequency_low_pass
		refdata[3] = 0
		#refdata[4] = varf
		refdata[5] = mask3D
		refdata[6] = False # whether to align on 50S, this only happens at refinement step
		user_func(refdata)
		#vol.write_image(os.path.join(outdir, "volf0000.hdf"), iref)
	mpi_barrier( MPI_COMM_WORLD )

	if CTF:
		if(data[0].get_attr_default("ctf_applied",0) > 0):  ERROR("mref_ali3d_MPI does not work for CTF-applied data", "mref_ali3d_MPI", 1, myid)
		from reconstruction import rec3D_MPI
	else:
		from reconstruction import rec3D_MPI_noCTF

	if debug:
		finfo.write( '%d loaded  \n' % len(data) )
		finfo.flush()

	#  this is needed for gathering of pixel errors
	disps = []
	recvcount = []
	for im in xrange(number_of_proc):
		if( im == main_node ):  disps.append(0)
		else:                   disps.append(disps[im-1] + recvcount[im-1])
		ib, ie = MPI_start_end(total_nima, number_of_proc, im)
		recvcount.append( ie - ib )

	total_iter = 0
	tr_dummy = Transform({"type":"spider"})

	if(focus != None):
		if(myid == main_node):
			vol = get_im(focus)
		else:
			vol =  model_blank(nx, nx, nx)
		bcast_EMData_to_all(vol, myid, main_node)
		focus, kb = prep_vol(vol)

	Niter = int(lstp*maxit*(nassign + nrefine) )
	for Iter in xrange(Niter):
		N_step = (Iter%(lstp*(nassign+nrefine)))/(nassign+nrefine)
		if Iter%(nassign+nrefine) < nassign:
			runtype = "ASSIGNMENT"
		else:
			runtype = "REFINEMENT"

		total_iter += 1
		if(myid == main_node):
			log.add("\n%s ITERATION #%3d,  inner iteration #%3d\nDelta = %4.1f, an = %5.2f, xrange = %5.2f, yrange = %5.2f, step = %5.2f  \
			"%(runtype, total_iter, Iter, delta[N_step], an[N_step], xrng[N_step],yrng[N_step],step[N_step]))
			start_ime = time()
	
		peaks =  [ [ -1.0e23 for im in xrange(nima) ] for iref in xrange(numref) ]
		if runtype=="REFINEMENT":
 			trans = [ [ tr_dummy for im in xrange(nima) ] for iref in xrange(numref) ]
			pixer = [ [  0.0     for im in xrange(nima) ] for iref in xrange(numref) ]
			if(an[N_step] > 0):
				from utilities    import even_angles
				ref_angles = even_angles(delta[N_step], symmetry=sym, method = ref_a, phiEqpsi = "Zero")
				# generate list of angles
				from alignment import generate_list_of_reference_angles_for_search
				list_of_reference_angles = \
				generate_list_of_reference_angles_for_search(ref_angles, sym=sym)
				del ref_angles
			else:  list_of_reference_angles = [[1.0,1.0]]

		cs = [0.0]*3
		for iref in xrange(numref):
			if(myid == main_node):
				volft = get_im(os.path.join(outdir, "volf%04d.hdf"%(total_iter-1)), iref)
			else:
				volft =  model_blank(nx, nx, nx)
			bcast_EMData_to_all(volft, myid, main_node)

			volft, kb = prep_vol(volft)
			if CTF:
				previous_defocus = -1.0
				if runtype=="REFINEMENT":
					start_time = time()
					prjref = prgq( volft, kb, nx, delta[N_step], ref_a, sym, MPI=True)
					if(myid == 0):
						log.add( "Calculation of projections: %d" % (time()-start_time) );start_time = time()
					del volft, kb

			else:
				if runtype=="REFINEMENT":
					start_time = time()
					refrings = prepare_refrings( volft, kb, nx, delta[N_step], ref_a, sym, numr)
					if(myid == 0):
						log.add( "Initial time to prepare rings: %d" % (time()-start_time) );start_time = time()
					del volft, kb


			start_time = time()
			for im in xrange(nima):
				if(CTF):
					ctf = data[im].get_attr( "ctf" )
					if runtype=="REFINEMENT":
						if(ctf.defocus != previous_defocus):
							previous_defocus = ctf.defocus
							rstart_time = time()
							refrings = gen_rings_ctf( prjref, nx, ctf, numr)
							if(myid == 0):
								log.add( "Repeated time to prepare rings: %d" % (time()-rstart_time) );rstart_time = time()

				if runtype=="ASSIGNMENT":
					phi,tht,psi,s2x,s2y = get_params_proj(data[im])
					ref = prgs( volft, kb, [phi,tht,psi,-s2x,-s2y])
					if CTF:  ref = filt_ctf( ref, ctf )
					if(focus != None):  mask2D = binarize( prgs( focus, kb, [phi,tht,psi,-s2x,-s2y]) )  #  Should be precalculated!!
					peak = ref.cmp("ccc",data[im],{"mask":mask2D, "negative":0})
					if not(finfo is None):
						finfo.write( "ID, iref, peak: %6d %d %8.5f\n" % (list_of_particles[im],iref,peak) )
				else:
					if(an[N_step] == -1):
						peak, pixel_error = proj_ali_incore(data[im], refrings, numr, xrng[N_step], yrng[N_step], step[N_step])
					else:
						peak, pixel_error = proj_ali_incore_local(data[im], refrings, list_of_reference_angles, numr, \
																	xrng[N_step], yrng[N_step], step[N_step], an[N_step],sym=sym)
					if not(finfo is None):
						phi,tht,psi,s2x,s2y = get_params_proj(data[im])
						finfo.write( "ID, iref, peak,t rans: %6d %d %f %f %f %f %f %f\n"%(list_of_particles[im],iref,peak,phi,tht,psi,s2x,s2y) )
						finfo.flush()

				peaks[iref][im] = peak
				if runtype=="REFINEMENT":
					pixer[iref][im] = pixel_error
					trans[iref][im] = data[im].get_attr( "xform.projection" )

			if(myid == 0):
				log.add( "Time to process particles for reference %3d: %d" % (iref, time()-start_time) );start_time = time()


		if runtype=="ASSIGNMENT":  del volft, kb, ref
		else:
			if CTF: del prjref
			del refrings
			if(an[N_step] > 0): del list_of_reference_angles


		#  send peak values to the main node, do the assignments, and bring them back
		from numpy import float32, empty, inner, abs
		if( myid == 0 ):
			dtot = empty( (numref, total_nima), dtype = float32)
		for  iref in xrange(numref):
			recvbuf = mpi_gatherv(peaks[iref], nima, MPI_FLOAT, recvcount, disps, MPI_FLOAT, main_node, MPI_COMM_WORLD)
			if( myid == 0 ): dtot[iref] = recvbuf
		del recvbuf


		#  The while loop over even angles delta should start here.
		#  prepare reference directions
		from utilities import even_angles, getvec
		refa = even_angles(60.0)
		numrefang = len(refa)
		refanorm = empty( (numrefang, 3), dtype = float32)
		for i in xrange(numrefang):
			tmp = getvec(refa[i][0], refa[i][1])
			for j in xrange(3):
				refanorm[i][j] = tmp[j]
		del  refa, tmp

		transv = empty( (nima, 3), dtype = float32)
		if runtype=="ASSIGNMENT":
			for im in xrange(nima):
				trns = data[im].get_attr( "xform.projection" )
				for j in xrange(3):
					transv[im][j] = trns.at(2,j)
		else:
			# For REFINEMENT we have a problem, as the exact angle is known only after the next step of assigning projections.
			# So, we will assume it is the one with max peak
			for im in xrange(nima):
				qt = -1.0e23
				it = -1
				for iref in xrange(numref):
					pt = peaks[iref][im]
					if(pt > qt):
						qt = pt
						it = iref
				for j in xrange(3):
					transv[im][j] = trans[it][im].at(2,j)
		#  We have all vectors, now create a list of assignments of images to references
		refassign = [-1]*nima
		for im in xrange(nima):
			refassign[im] = abs(inner(refanorm,transv[im])).argmax()
		assigntorefa = mpi_gatherv(refassign, nima, MPI_INT, recvcount, disps, MPI_INT, main_node, MPI_COMM_WORLD)
		assigntorefa = map(int, assigntorefa)

		del refassign, refanorm, transv


		"""
		#  Trying to use ISAC code for EQ-Kmeans  PAP 03/21/2015
		if myid == main_node:

			for imrefa in xrange(numrefang):
				from utilities import findall
				N = findall(imrefa, assigntorefa)
				current_nima = len(N)
				if( current_nima >= numref and report_error == 0):
					tasi = [[] for iref in xrange(numref)]
					maxasi = current_nima//numref
					nt = current_nima
					kt = numref
					K = range(numref)

					d = empty( (numref, current_nima), dtype = float32)
					for ima in xrange(current_nima):
						for iref in xrange(numref):  d[iref][ima] = dtot[iref][N[ima]]

			d = empty( (numref, total_nima), dtype = float32)
			for ima in xrange(total_nima):
				for iref in xrange(numref):  d[iref][ima] = dtot[iref][N[ima]]
			id_list_long = Util.assign_groups(str(d.__array_interface__['data'][0]), numref, nima) # string with memory address is passed as parameters
			del d
			id_list = [[] for i in xrange(numref)]
			maxasi = total_nima/numref
			for i in xrange(maxasi*numref):
				id_list[i/maxasi].append(id_list_long[i])
			for i in xrange(total_nima%maxasi):
				id_list[id_list_long[-1]].append(id_list_long[maxasi*numref+i])
			for iref in xrange(numref):
				id_list[iref].sort()

			assignment = [0]*total_nima
			for iref in xrange(numref):
				for im in id_list[iref]: assignment[im] = iref
		else:
			assignment = [0]*total_nima
		mpi_barrier(MPI_COMM_WORLD)
		#belongsto = mpi_bcast(belongsto, nima, MPI_INT, main_node, MPI_COMM_WORLD)
		#belongsto = map(int, belongsto)
		"""


		if myid == main_node:
			SA = False
			asi = [[] for iref in xrange(numref)]
			report_error = 0
			for imrefa in xrange(numrefang):
				from utilities import findall
				N = findall(imrefa, assigntorefa)
				current_nima = len(N)
				if( current_nima >= numref and report_error == 0):
					tasi = [[] for iref in xrange(numref)]
					maxasi = current_nima//numref
					nt = current_nima
					kt = numref
					K = range(numref)

					d = empty( (numref, current_nima), dtype = float32)
					for ima in xrange(current_nima):
						for iref in xrange(numref):  d[iref][ima] = dtot[iref][N[ima]]

					while nt > 0 and kt > 0:
						l = d.argmax()
						group = l//current_nima
						ima   = l-current_nima*group
						if SA:
							J = [0.0]*numref
							sJ = 0
							Jc = [0.0]*numref
							for iref in xrange(numref):
								J[iref] = exp(d[iref][ima]/T)
								sJ += J[iref]
							for iref in xrange(numref):
								J[iref] /= sJ
							Jc[0] = J[0]
							for iref in xrange(1, numref):
								Jc[iref] = Jc[iref-1]+J[iref]
							sss = random()
							for group in xrange(numref):
								if( sss <= Jc[group]): break
						tasi[group].append(N[ima])
						N[ima] = -1
						for iref in xrange(numref):  d[iref][ima] = -1.e10
						nt -= 1
						masi = len(tasi[group])
						if masi == maxasi:
							for im in xrange(current_nima):  d[group][im] = -1.e10
							kt -= 1
					else:
						for ima in xrange(current_nima):
							if N[ima] > -1:
								qm = -1.e10
								for iref in xrange(numref):
									qt = dtot[iref][N[ima]]
									if( qt > qm ):
										qm = qt
										group = iref
								tasi[group].append(N[ima])

					del d, N, K
					if  SA:  del J, Jc
					for iref in xrange(numref):
						asi[iref] += tasi[iref]
					del tasi
				else:
					report_error = 1
			#  This should be deleted only once we know that the number of images is sufficiently large, see below.
			del dtot

		else:
			assignment = []
			report_error = 0

		report_error = bcast_number_to_all(report_error, source_node = main_node)
		if report_error == 1:  ERROR('Number of images within a group too small', "mref_ali3d_MPI", 1, myid)
		if myid == main_node:
			assignment = [0]*total_nima
			for iref in xrange(numref):
				for im in xrange(len(asi[iref])):
					assignment[asi[iref][im]] = iref
			del asi
		
		"""
		if myid == main_node:
			assignment = [0]*total_nima
			for iref in xrange(numref):
				for im in xrange(len(asi[iref])):
					assignment[asi[iref][im]] = iref
			del asi
		"""

		'''
		if myid == main_node:
			SA = False
			maxasi = total_nima//numref
			asi = [[] for iref in xrange(numref)]
			nt = total_nima
			kt = numref
			K = range(numref)
			N = range(total_nima)

			while nt > 0 and kt > 1:
				l = d.argmax()
				group = l//total_nima
				ima   = l-total_nima*group
				if SA:
					J = [0.0]*numref
					sJ = 0
					Jc = [0.0]*numref
					for iref in xrange(numref):
						J[iref] = exp(d[iref][ima]/T)
						sJ += J[iref]
					for iref in xrange(numref):
						J[iref] /= sJ
					Jc[0] = J[0]
					for iref in xrange(1, numref):
						Jc[iref] = Jc[iref-1]+J[iref]
					sss = random()
					for group in xrange(numref):
						if( sss <= Jc[group]): break
				asi[group].append(N[ima])
				for iref in xrange(numref):  d[iref][ima] = -1.e10
				nt -= 1
				masi = len(asi[group])
				if masi == maxasi:
					for im in xrange(total_nima):  d[group][im] = -1.e10
					kt -= 1
			else:
				mas = [len(asi[iref]) for iref in xrange(numref)]
				group = mas.index(min(mas))
				del mas
				for im in xrange(total_nima):
					kt = 0
					go = True
					while(go and kt < numref):
						if d[kt][im] > -1.e10:
							asi[group].append(im)
							go = False
						kt += 1

			assignment = [0]*total_nima
			for iref in xrange(numref):
				for im in xrange(len(asi[iref])):
					assignment[asi[iref][im]] = iref

			del asi, d, N, K
			if  SA:  del J, Jc


		else:
			assignment = []
		'''

		assignment = mpi_scatterv(assignment, recvcount, disps, MPI_INT, recvcount[myid], MPI_INT, main_node, MPI_COMM_WORLD)
		assignment = map(int, assignment)


		#  compute number of particles that changed assignment and how many are in which group
		nchng = 0
		npergroup = [0]*numref
		for im in xrange(nima):
			iref = data[im].get_attr('group')
			npergroup[assignment[im]] += 1
			if( iref != assignment[im]): nchng += 1
			data[im].set_attr('group', assignment[im])
		nchng = mpi_reduce(nchng, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD)
		npergroup = mpi_reduce(npergroup, numref, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD)
		npergroup = map(int, npergroup)
		terminate = 0
		if( myid == 0 ):
			nchng = int(nchng[0])
			precn = 100*float(nchng)/float(total_nima)
			msg = " Number of particles that changed assignments %7d, percentage of total: %5.1f"%(nchng, precn)
			log.add(msg)
			msg = " Group       number of particles"
			log.add(msg)
			for iref in xrange(numref):
				msg = " %5d       %7d"%(iref+1, npergroup[iref])
				log.add(msg)
			if(precn <= termprec):  terminate = 1
		terminate = mpi_bcast(terminate, 1, MPI_INT, 0, MPI_COMM_WORLD)
		terminate = int(terminate[0])

		if runtype=="REFINEMENT":
			for im in xrange(nima):
				data[im].set_attr('xform.projection', trans[assignment[im]][im])
				pixer[0][im] = pixer[assignment[im]][im]
			pixer = pixer[0]

			if(center == -1):
				cs[0], cs[1], cs[2], dummy, dummy = estimate_3D_center_MPI(data, total_nima, myid, number_of_proc, main_node)				
				if myid == main_node:
					msg = " Average center x = %10.3f        Center y = %10.3f        Center z = %10.3f"%(cs[0], cs[1], cs[2])
					log.add(msg)
				cs = mpi_bcast(cs, 3, MPI_FLOAT, main_node, MPI_COMM_WORLD)
				cs = [-float(cs[0]), -float(cs[1]), -float(cs[2])]
				rotate_3D_shift(data, cs)
			#output pixel errors
			recvbuf = mpi_gatherv(pixer, nima, MPI_FLOAT, recvcount, disps, MPI_FLOAT, main_node, MPI_COMM_WORLD)
			mpi_barrier(MPI_COMM_WORLD)
			if(myid == main_node):
				recvbuf = map(float, recvbuf)
				from statistics import hist_list
				lhist = 20
				region, histo = hist_list(recvbuf, lhist)
				if(region[0] < 0.0):  region[0] = 0.0
				msg = "      Histogram of pixel errors\n      ERROR       number of particles"
				log.add(msg)
				for lhx in xrange(lhist):
					msg = " %10.3f      %7d"%(region[lhx], histo[lhx])
					log.add(msg)
				del region, histo
			del recvbuf

		fscc = [None]*numref

		if fourvar and runtype=="REFINEMENT":
			sumvol = model_blank(nx, nx, nx)

		start_time = time()
		for iref in xrange(numref):
			#  3D stuff
			from time import localtime, strftime
			if(CTF): volref, fscc[iref] = rec3D_MPI(data, snr, sym, model_circle(last_ring, nx, nx, nx),\
			 os.path.join(outdir, "resolution_%02d_%04d"%(iref, total_iter)), myid, main_node, index = iref, npad = npad, finfo=frec)
			else:    volref, fscc[iref] = rec3D_MPI_noCTF(data, sym, model_circle(last_ring, nx, nx, nx),\
			 os.path.join(outdir, "resolution_%02d_%04d"%(iref, total_iter)), myid, main_node, index = iref, npad = npad, finfo=frec)
			if(myid == 0):
				log.add( "Time to compute 3D: %d" % (time()-start_time) );start_time = time()

			if(myid == main_node):
				volref.write_image(os.path.join(outdir, "vol%04d.hdf"%( total_iter)), iref)
				if fourvar and runtype=="REFINEMENT":
					sumvol += volref
			del volref

		if runtype=="REFINEMENT":
			if fourvar:
				varf = varf3d_MPI(data, os.path.join(outdir, "ssnr%04d"%total_iter), None,sumvol,last_ring, 1.0, 1, CTF, 1, sym, myid)
				if myid == main_node:   
					varf = 1.0/varf
					varf.write_image( os.path.join(outdir,"varf%04d.hdf"%total_iter) )

		if(myid == main_node):
			frcs={}
			for iref in xrange(numref):
				frc=read_text_file(os.path.join(outdir, "resolution_%02d_%04d"%(iref, total_iter)),-1)
				frcs[iref]=frc
			refdata = [None]*7
			refdata[0] = numref
			refdata[1] = outdir
			refdata[2] = frequency_low_pass
			refdata[3] = total_iter
			refdata[4] = varf
			refdata[5] = mask3D
			refdata[6] = (runtype=="REFINEMENT") # whether to align on 50S, this only happens at refinement step
			user_func(refdata)

		mpi_barrier(MPI_COMM_WORLD)
		if terminate ==0: # headers are only updated when the program is going to terminate
			start_time = time()
			if nrefine!=0:
				par_str = ['xform.projection', 'ID', 'group']
			else:
				par_str = ['group', 'ID' ]
	        	if myid == main_node:
				from utilities import file_type
	        		if(file_type(stack) == "bdb"):
	        			from utilities import recv_attr_dict_bdb
	        			recv_attr_dict_bdb(main_node, stack, data, par_str, image_start, image_end, number_of_proc)
	        		else:
	        			from utilities import recv_attr_dict
	        			recv_attr_dict(main_node, stack, data, par_str, image_start, image_end, number_of_proc)
	        	else:		send_attr_dict(main_node, data, par_str, image_start, image_end)
			if(myid == 0):
				log.add( "Time to write headers: %d\n" % (time()-start_time) );start_time = time()
			mpi_barrier(MPI_COMM_WORLD)
			if myid==main_node:
				log.add("mref_ali3d_MPI terminated due to small number of objects changing assignments")
				#from utilities import cmdexecute
                        	#cmd = "{} {} {} {}".format("sxheader.py",stack,"--params=xform.projection", "--export="+os.path.join(outdir, "ali3d_params_%03d.txt"%total_iter))
                        	#cmdexecute(cmd)
			mpi_barrier(MPI_COMM_WORLD)		
		if terminate==1:
			break
	if myid==main_node:
		log.add("mref_ali3d_MPI finishes")
		from utilities import cmdexecute
		cmd = "{} {} {} {}".format("sxheader.py",stack,"--params=xform.projection", "--export="+os.path.join(outdir,"ali3d_params.txt"))
		cmdexecute(cmd)
def print_upper_triangular_matrix(data_table_dict,N_indep,log_main):
        msg =""
        for i in xrange(N_indep):
                msg +="%7d"%i
        log_main.add(msg)
        for i in xrange(N_indep):
                msg ="%5d "%i
                for j in xrange(N_indep):
                        if i<j:
                                msg +="%5.2f "%data_table_dict[(i,j)]
                        else:
                                msg +="      "
                log_main.add(msg)
			
def print_a_line_with_timestamp(string_to_be_printed ):                 
	line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
 	print(line,string_to_be_printed)
	return string_to_be_printed

def convertasi(asig,K):
        p = []
        for k in xrange(K):
                l = []
                for i in xrange(len(asig)):
                        if(asig[i] == k): l.append(i)
                l = array(l, 'int32')
                l.sort()
                p.append(l)
        return p

def prepare_ptp(data_list,K):
	num_of_pt=len(data_list)
	ptp=[]
	for ipt in xrange(num_of_pt):
		ptp.append([])
	for ipt in xrange(num_of_pt):
		nc =len(data_list[ipt])
		asig=[-1]*nc
		for i in xrange(nc):
        		asig[i]=data_list[ipt][i]
		ptp[ipt] = convertasi(asig,K)
	return ptp

def print_dict(dict,theme):
        line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
        print(line,theme)
        spaces = "                    "
        for key, value in sorted( dict.items() ):
                if(key != "constants"):  print("                    => ", key+spaces[len(key):],":  ",value)

def checkstep(item, keepchecking, myid, main_node):
	from utilities import bcast_number_to_all
        if(myid == main_node):
                if keepchecking:
                        if(os.path.exists(item)):
                                doit = 0
                        else:
                                doit = 1
                                keepchecking = False
                else:
                        doit = 1
        else:
                doit = 1
        doit = bcast_number_to_all(doit, source_node = main_node)
        return doit, keepchecking

def get_resolution_mrk01(vol, radi, nnxo, fscoutputdir, mask_option):
        # this function is single processor
        #  Get updated FSC curves, user can also provide a mask using radi variable
        import types
	from statistics import fsc
	from utilities import model_circle, get_im
	from filter import fit_tanh1
        if(type(radi) == int):
                if(mask_option is None):  mask = model_circle(radi,nnxo,nnxo,nnxo)
                else:                           mask = get_im(mask_option)
        else:  mask = radi
        nfsc = fsc(vol[0]*mask,vol[1]*mask, 1.0,os.path.join(fscoutputdir,"fsc.txt") )
        currentres = -1.0
        ns = len(nfsc[1])
        #  This is actual resolution, as computed by 2*f/(1+f)
        for i in xrange(1,ns-1):
                if ( nfsc[1][i] < 0.333333333333333333333333):
                        currentres = nfsc[0][i-1]
                        break
        if(currentres < 0.0):
                print("  Something wrong with the resolution, cannot continue")
                mpi_finalize()
                exit()
        """
        lowpass = 0.5
        ns = len(nfsc[1])
        #  This is resolution used to filter half-volumes
        for i in xrange(1,ns-1):
                if ( nfsc[1][i] < 0.5 ):
                        lowpass = nfsc[0][i-1]
                        break
        """
        lowpass, falloff = fit_tanh1(nfsc, 0.01)

        return  round(lowpass,4), round(falloff,4), round(currentres,2)

def get_K_equal_assignment(input_list,K):
	import types
	from utilities import read_text_file
	if type(K) != IntType:
		return None
	if type(input_list) == StringType:
		 ll=read_text_file(input_list)
	elif type(input_list) == types.ListType:
		ll=input_list
	else:
		print "Error"
		return None
def partition_to_groups(alist,K):
	res =[]
	for igroup in xrange(K):
		this_group =[]
		for imeb in xrange(len(alist)):
			if alist[imeb] ==igroup:
				this_group.append(imeb)
		this_group.sort()
		res.append(this_group)
	return res

def partition_independent_runs(run_list,K):
	indep_runs_groups={}
	for indep in xrange(len(run_list)):
		this_run = run_list[indep]
		groups = partition_to_groups(this_run,K)
		indep_runs_groups[indep]=groups 
	return indep_runs_groups

def get_outliers(total_number,plist):
        tlist={}
        for i in xrange(total_number):
                tlist[i]=i
        for a in plist:
                del tlist[a]
        out =[]
        for a in tlist:
                out.append(a)
        return out

def merge_groups(stable_members_list):
	alist=[]
	for i in xrange(len(stable_members_list)):
		for j in xrange(len(stable_members_list[i])):
			alist.append(stable_members_list[i][j])
	return alist

def save_alist(Tracker,name_of_the_text_file,alist):
	from utilities import write_text_file
	import os
	log       =Tracker["constants"]["log_main"]
	myid      =Tracker["constants"]["myid"]
	main_node =Tracker["constants"]["main_node"]
	dir_to_save_list =Tracker["this_dir"]
	if myid==main_node:
		file_name=os.path.join(dir_to_save_list,name_of_the_text_file)
		write_text_file(alist, file_name)

def reconstruct_grouped_vols(Tracker,name_of_grouped_vols,grouped_plist):
	from mpi import  MPI_COMM_WORLD,mpi_barrier
	from applications import recons3d_n_MPI
	from utilities import cmdexecute
	import os
	myid		 = Tracker["constants"]["myid"]
	main_node	 = Tracker["constants"]["main_node"]
	main_dir         = Tracker["this_dir"]
	number_of_groups = Tracker["number_of_groups"]
	data_stack       = Tracker["this_data_stack"]
	log              = Tracker["constants"]["log_main"]
	####
	if myid ==main_node:
		log.add("Reconstruct a group of volumes")
	for igrp in xrange(number_of_groups):
		#pid_list=Tracker["two_way_stable_member"][igrp]
		pid_list =grouped_plist[igrp]
                vol_stack=os.path.join(main_dir,"TMP_init%03d.hdf"%igrp)
                tlist=[]
                for b in pid_list:
                	tlist.append(int(b))
                recons3d_n_MPI(data_stack,tlist,vol_stack,Tracker["constants"]["CTF"],Tracker["constants"]["snr"],\
                Tracker["constants"]["sign"],Tracker["constants"]["npad"],Tracker["constants"]["sym"], \
		Tracker["constants"]["listfile"],Tracker["constants"]["group"],\
                Tracker["constants"]["verbose"],Tracker["constants"]["xysize"],Tracker["constants"]["zsize"])
                newvols=os.path.join(main_dir,"TMP_init*.hdf")
        if myid==main_node:
       		#cmd ="{} {} {}".format("sxcpy.py",newvols,Tracker["EKKMREF"][inew]["refvols"]) # Make stacks
		cmd ="{} {} {}".format("sxcpy.py",newvols,name_of_grouped_vols)  
        	cmdexecute(cmd)
                cmd ="{} {}".format("rm",newvols)
                cmdexecute(cmd)
        mpi_barrier(MPI_COMM_WORLD)

def N_independent_reconstructions(Tracker):
	from applications import recons3d_n_MPI
	from utilities import write_text_file, read_text_file, cmdexecute
	import os
	from mpi import mpi_barrier,MPI_COMM_WORLD
	from random import shuffle
	myid      =Tracker["constants"]["myid"]
	main_node =Tracker["constants"]["main_node"]
	initdir   =Tracker["this_dir"]
	data_stack=Tracker["this_data_stack"]
	total_stack =Tracker["this_total_stack"] 
	log_main =Tracker["constants"]["log_main"]
	if myid ==main_node:
		log_main.add("-----------Independent reconstructions---------------")
        for irandom in xrange(Tracker["constants"]["indep_runs"]):
        	ll=range(total_stack)
                shuffle(ll)
                if myid ==main_node:
                	log_main.add("Initial random assignments     "+os.path.join(initdir,"random_list%d.txt"%irandom))
                        write_text_file(ll,os.path.join(initdir,"random_list%d.txt"%irandom))
                        log_main.add("preset ali3d_parameters ")
	mpi_barrier(MPI_COMM_WORLD)
	if myid ==main_node:
         	if Tracker["importali3d"]!="":
         		cmd = "{} {} {} {}".format("sxheader.py",data_stack,"--params=xform.projection", "--import="+ \
				Tracker["importali3d"])
                	cmdexecute(cmd)
			log_main.add("ali3d_parameters is preset to "+Tracker["importali3d"])
		else:
			log_main.add("ali3d parameters in this run are not altered !")
	for iter_indep in xrange(Tracker["constants"]["indep_runs"]):
       		ll=read_text_file(os.path.join(initdir,"random_list%d.txt"%iter_indep))
                linit_list =[]
                for igrp in xrange(Tracker["number_of_groups"]):
                       	alist=ll[(total_stack*igrp)//Tracker["number_of_groups"]:(total_stack*(igrp+1))//Tracker["number_of_groups"]]
                       	alist.sort()
                       	linit_list.append(alist)
                	linit_list_file_name=os.path.join(initdir,"list1_grp%0d_%0d.txt"%(igrp,iter_indep))
			if myid ==main_node:
                		write_text_file(alist,linit_list_file_name)
                       	volstack =os.path.join(initdir,"TMP_vol%03d.hdf"%igrp)
                       	recons3d_n_MPI(data_stack,alist,volstack,Tracker["constants"]["CTF"],Tracker["constants"]["snr"],Tracker["constants"]["sign"],\
                       	Tracker["constants"]["npad"],Tracker["constants"]["sym"],Tracker["constants"]["listfile"],Tracker["constants"]["group"],\
			Tracker["constants"]["verbose"],Tracker["constants"]["xysize"],Tracker["constants"]["zsize"])
                newvols=os.path.join(initdir,"TMP_vol*.hdf")
                computed_refvol_name=os.path.join(initdir,"vol_run_%03d.hdf"%iter_indep)
                filtered_volumes =os.path.join(initdir,"volf_run_%03d.hdf"%iter_indep)
                filter_options ="--process=filter.lowpass.tanh:cutoff_abs="+str(Tracker["low_pass"])+":fall_off=.1"
                if myid==main_node:
                       cmd ="{} {} {}".format("sxcpy.py",newvols,computed_refvol_name) # Make stacks 
                       cmdexecute(cmd)
                       cmd ="{} {}".format("rm",newvols)
                       cmdexecute(cmd)
                       cmd ="{}  {}  {} {}".format("e2proc3d.py",computed_refvol_name,filtered_volumes,filter_options)
                       cmdexecute(cmd)
		mpi_barrier(MPI_COMM_WORLD)

def N_independent_Kmref(mpi_comm,Tracker):
	from mpi import mpi_barrier,MPI_COMM_WORLD
	from applications import Kmref_ali3d_MPI
	from utilities import cmdexecute
	from logger import Logger,BaseLogger_Files
	myid            = Tracker["constants"]["myid"]
	main_node       = Tracker["constants"]["main_node"]
	log_main        = Tracker["constants"]["log_main"]
	data_stack      = Tracker["this_data_stack"]
	if myid ==main_node:
		log_main.add("-----------%5d independent runs of Kmref -----------------"%Tracker["constants"]["indep_runs"])
	for iter_indep in xrange(Tracker["constants"]["indep_runs"]):
        	log_Kmref       =Logger(BaseLogger_Files())
                log_Kmref.prefix=Tracker["KMREF"][iter_indep]["output_dir"]+"/"
                empty_group =Kmref_ali3d_MPI(data_stack,Tracker["KMREF"][iter_indep]["refvols"],Tracker["KMREF"][iter_indep]["output_dir"], \
		Tracker["constants"]["mask3D"],Tracker["constants"]["focus3Dmask"],Tracker["constants"]["maxit"],\
		Tracker["constants"]["ir"],Tracker["constants"]["radius"],Tracker["constants"]["rs"],\
                Tracker["constants"]["xr"],Tracker["constants"]["yr"],Tracker["constants"]["ts"],Tracker["constants"]["delta"],\
		Tracker["constants"]["nassign"],Tracker["constants"]["nrefine"],Tracker["constants"]["an"],Tracker["constants"]["center"],\
                Tracker["constants"]["CTF"],Tracker["constants"]["snr"],Tracker["constants"]["ref_a"],Tracker["constants"]["sym"],\
                Tracker["constants"]["user_func"],Tracker["constants"]["npad"],Tracker["constants"]["debug"],\
		Tracker["constants"]["fourvar"],Tracker["constants"]["stoprnct"],mpi_comm,log_Kmref)
                if myid==main_node:
                	cmd = "{} {} {} {}".format("sxheader.py",data_stack,"--params=group", \
					"--export="+Tracker["KMREF"][iter_indep]["partition"])
                	cmdexecute(cmd)
                mpi_barrier(MPI_COMM_WORLD)
	
def do_independent_EKmref(Tracker):
	from mpi import mpi_barrier,MPI_COMM_WORLD
	from applications import Kmref_ali3d_MPI
	from utilities import cmdexecute
	from logger import Logger,BaseLogger_Files
	myid            = Tracker["constants"]["myid"]
        main_node       = Tracker["constants"]["main_node"]
        log_main        = Tracker["constants"]["log_main"]
        data_stack      = Tracker["this_data_stack"]
	iruns           = Tracker["this_iruns"]
	for iter_indep in xrange(iruns):
		log_Kmref       =Logger(BaseLogger_Files())
		log_Kmref.prefix=Tracker["EKMREF"][iter_indep]["output_dir"]+"/"
		empty_group =Kmref_ali3d_MPI(data_stack,Tracker["EKMREF"][iter_indep]["refvols"],Tracker["EKMREF"][iter_indep]["output_dir"], \
                Tracker["constants"]["mask3D"],Tracker["constants"]["focus3Dmask"],Tracker["constants"]["maxit"],\
                Tracker["constants"]["ir"],Tracker["constants"]["radius"],Tracker["constants"]["rs"],\
                Tracker["constants"]["xr"],Tracker["constants"]["yr"],Tracker["constants"]["ts"],Tracker["constants"]["delta"],\
		Tracker["constants"]["an"],Tracker["constants"]["center"],\
                Tracker["constants"]["nassign"],Tracker["constants"]["nrefine"],\
                Tracker["constants"]["CTF"],Tracker["constants"]["snr"],Tracker["constants"]["ref_a"],Tracker["constants"]["sym"],\
                Tracker["constants"]["user_func"],Tracker["constants"]["npad"],Tracker["constants"]["debug"],\
                Tracker["constants"]["fourvar"],Tracker["constants"]["stoprnct"],Tracker["constants"]["mpi_comm"],log_Kmref)
                if myid==main_node:
                        cmd = "{} {} {} {}".format("sxheader.py",data_stack,"--params=group",\
					 "--export="+Tracker["EKMREF"][iter_indep]["partition"])
                        cmdexecute(cmd)
                mpi_barrier(MPI_COMM_WORLD)

def N_independent_mref(mpi_comm,Tracker):
	from mpi import mpi_barrier,MPI_COMM_WORLD
	from applications import mref_ali3d_MPI
	from utilities import cmdexecute
	from logger import Logger,BaseLogger_Files
	myid      = Tracker["constants"]["myid"]
	main_node = Tracker["constants"]["main_node"]
	log_main  = Tracker["constants"]["log_main"]
	data_stack= Tracker["this_data_stack"]
	if myid ==main_node:
        	log_main.add("-----------%5d independent runs of Equal Kmeans -----------------"%Tracker["constants"]["indep_runs"])
	for iter_indep in xrange(Tracker["constants"]["indep_runs"]):
		if myid==main_node:
			log_main.add(" %d th run  starts !  "%iter_indep)
        	outdir = Tracker["EMREF"][iter_indep]["output_dir"]
                #doit, keepchecking = checkstep(outdir, keepchecking, myid, main_node)
                log_Emref=Logger(BaseLogger_Files())
                log_Emref.prefix=outdir+"/"
                if Tracker["importali3d"]!="" and myid ==main_node:
                	cmd = "{} {} {} {}".format("sxheader.py",data_stack,"--params=xform.projection", \
				 "--import="+Tracker["importali3d"])
                        cmdexecute(cmd)
               	mpi_barrier(MPI_COMM_WORLD)
               	mref_ali3d_MPI_beta(data_stack,Tracker["EMREF"][iter_indep]["refvols"],outdir,Tracker["constants"]["mask3D"] ,\
               	Tracker["constants"]["focus3Dmask"],Tracker["constants"]["maxit"],Tracker["constants"]["ir"],\
		Tracker["constants"]["radius"],Tracker["constants"]["rs"],Tracker["constants"]["xr"],\
		Tracker["constants"]["yr"],Tracker["constants"]["ts"],\
               	Tracker["constants"]["delta"],Tracker["constants"]["an"],Tracker["constants"]["center"],\
		Tracker["constants"]["nassign"],Tracker["constants"]["nrefine"],\
               	Tracker["constants"]["CTF"],Tracker["constants"]["snr"],Tracker["constants"]["ref_a"],Tracker["constants"]["sym"],\
               	Tracker["constants"]["user_func"],Tracker["constants"]["npad"],Tracker["constants"]["debug"],\
		Tracker["constants"]["fourvar"],Tracker["constants"]["stoprnct"], mpi_comm,log_Emref,Tracker["frequency_low_pass"])
               	if myid==main_node:
             		cmd = "{} {} {} {}".format("sxheader.py",Tracker["this_data_stack"],"--params=group",\
			 "--export="+Tracker["EMREF"][iter_indep]["partition"])
                        cmdexecute(cmd)
                        if Tracker["constants"]["nrefine"] !=0:
                        	cmd = "{} {} {} {}".format("sxheader.py",Tracker["this_data_stack"],"--params=xform.projection",\
				 "--export="+Tracker["EMREF"][iter_indep]["ali3d"])
                        	cmdexecute(cmd)
		mpi_barrier(MPI_COMM_WORLD)

def prepare_EMREF_dict(Tracker):
	main_dir =Tracker["this_dir"]
	EMREF_iteration  ={}
        for iter_indep in xrange(Tracker["constants"]["indep_runs"]):
        	output_for_this_run ={}
        	output_for_this_run["output_dir"]= os.path.join(main_dir,"EMREF_independent_%03d"%iter_indep)
        	output_for_this_run["partition"] = os.path.join(output_for_this_run["output_dir"],"list2.txt")
        	output_for_this_run["refvols"]   = os.path.join(main_dir,"volf_run_%03d.hdf"%iter_indep)
        	output_for_this_run["ali3d"]     = os.path.join(output_for_this_run["output_dir"],"ali3d_params.txt")
        	EMREF_iteration[iter_indep]=output_for_this_run
        Tracker["EMREF"]=EMREF_iteration

def prepare_KMREF_dict(Tracker):
	main_dir =Tracker["this_dir"]
	KMREF_iteration  ={}
        for iter_indep in xrange(Tracker["constants"]["indep_runs"]):
        	output_for_this_run ={}
                output_for_this_run["output_dir"]= os.path.join(main_dir,"KMREF_independent_%03d"%iter_indep)
                output_for_this_run["partition"] = os.path.join(output_for_this_run["output_dir"],"list2.txt")
                output_for_this_run["refvols"]   = os.path.join(main_dir,"volf_run_%03d.hdf"%iter_indep)
                output_for_this_run["ali3d"]     = os.path.join(output_for_this_run["output_dir"],"ali3d_params.txt")
                KMREF_iteration[iter_indep]=output_for_this_run
        Tracker["KMREF"]=KMREF_iteration

def set_refvols_table(initdir,vol_pattern,Tracker):
	# vol_pattern="volf_run_%03d.hdf"; volf_stable_%3d.hdf
	import os
	initial_random_vol={}
        for iter_indep in xrange(Tracker["constants"]["indep_runs"]):
        	filtered_volumes =os.path.join(initdir,vol_pattern%iter_indep)
                initial_random_vol[iter_indep]=filtered_volumes
        Tracker["refvols"]=initial_random_vol

def prepare_EKMREF_dict(Tracker):
	main_dir = Tracker["this_dir"]
	KMREF_iteration  ={}
        for iter_indep in xrange(Tracker["constants"]["indep_runs"]):
        	output_for_this_run ={}
                output_for_this_run["output_dir"]=os.path.join(main_dir,"EKMREF%03d"%iter_indep)
                output_for_this_run["partition"] =os.path.join(output_for_this_run["output_dir"],"list2.txt")
                output_for_this_run["ali3d"]     = os.path.join(output_for_this_run["output_dir"],"ali3d_params.txt")
                KMREF_iteration[iter_indep]=output_for_this_run
        Tracker["EKMREF"]=KMREF_iteration

def partition_ali3d_params_of_outliers(Tracker):
	from utilities import read_text_file,write_text_file
	from mpi import mpi_barrier, MPI_COMM_WORLD
	myid 	   = Tracker["constants"]["myid"]
	main_node  = Tracker["constants"]["main_node"]
	outlier_list=read_text_file(Tracker["this_unaccounted"])
        accounted_list=read_text_file(Tracker["this_accounted"])
	ali3d_params_of_outliers = []
	ali3d_params_of_accounted  = []
	if Tracker["constants"]["importali3d"] !="":
     		ali3d_params_of_outliers = []
        	ali3d_params_of_accounted  = []
		ali3d_params=read_text_file(Tracker["this_ali3d"])
		for i in xrange(len(outlier_list)):
			ali3d_params_of_outliers.append(ali3d_params[outlier_list[i]])
		if myid ==main_node:
			write_text_file(ali3d_params_of_outliers,Tracker["ali3d_of_outliers"])
		for i in xrange(len(accounted_list)):
			ali3d_params_of_accounted.append(ali3d_params[accounted_list[i]])
		if myid ==main_node:
			write_text_file(ali3d_params_of_accounted,Tracker["ali3d_of_accounted"])
	Tracker["number_of_unaccounted"]=len(outlier_list)
	Tracker["average_members_in_a_group"]= len(accounted_list)/float(Tracker["number_of_groups"])
	mpi_barrier(MPI_COMM_WORLD)

def do_two_way_comparison(Tracker):
	from mpi import mpi_barrier, MPI_COMM_WORLD
	from utilities import read_text_file
	from statistics import k_means_match_clusters_asg_new
	######
	myid              =Tracker["constants"]["myid"]
	main_node         =Tracker["constants"]["main_node"]
	log_main          =Tracker["constants"]["log_main"]
	total_stack       =Tracker["this_total_stack"]
	main_dir          =Tracker["this_dir"]
	number_of_groups  =Tracker["number_of_groups"]
	######
	if myid ==main_node:
		msg="-------Two_way comparisons analysis of %3d independent runs of equal Kmeans-------"%Tracker["constants"]["indep_runs"]
		log_main.add(msg)
        total_partition=[]
	if Tracker["constants"]["indep_runs"]<2:
		if myid ==main_node:
			log_main.add(" Error! One single run cannot make two-way comparison")
		from mip import mpi_finalize
		mpi_finalize()
                exit()
        for iter_indep in xrange(Tracker["constants"]["indep_runs"]):
		if Tracker["constants"]["mode"][0:1]!="K":
			partition_list=read_text_file(Tracker["EMREF"][iter_indep]["partition"])
			if myid ==main_node: 
				log_main.add("Equal Kmeans  with %d         %d"%(total_stack,number_of_groups)+"\n")
		else: 
			partition_list=read_text_file(Tracker["KMREF"][iter_indep]["partition"])
			if myid ==main_node: 
				log_main.add("Kmeans  with %d         %d"%(total_stack,number_of_groups)+"\n")
       		total_partition.append(partition_list)
       	### Two-way comparision is carried out on all nodes 
       	ptp=prepare_ptp(total_partition,number_of_groups)
	indep_runs_to_groups =partition_independent_runs(total_partition,number_of_groups)
	total_pop=0
	two_ways_stable_member_list={}
	avg_two_ways               =0.0
	avg_two_ways_square        =0.
	scores                     ={}
       	for iptp in xrange(len(ptp)):
        	for jptp in xrange(len(ptp)):
               		newindeces, list_stable, nb_tot_objs = k_means_match_clusters_asg_new(ptp[iptp],ptp[jptp])
                       	tt =0.
			if myid ==main_node and iptp<jptp:
				aline=print_a_line_with_timestamp("Two-way comparison between independent run %3d and %3d"%(iptp,jptp))
				log_main.add(aline)
                        for m in xrange(len(list_stable)):
                        	a=list_stable[m]
                               	tt +=len(a)
				if myid==main_node and iptp<jptp:
					aline=print_a_line_with_timestamp("Group %d  number of stable members %10d  "%(m,len(a)))
					log_main.add(aline)
					aline=print_a_line_with_timestamp("The comparison is between %3d th group of %3d th run and %3d th group of %3d th run" \
					 									%(newindeces[m][0],iptp,newindeces[m][1],jptp))
					log_main.add(aline)
					aline=print_a_line_with_timestamp("The  %3d th group of %3d th run contains %6d members" \
						        %(iptp,newindeces[m][0],len(indep_runs_to_groups[iptp][newindeces[m][0]])))
					log_main.add(aline)
					aline=print_a_line_with_timestamp("The  %3d th group of %3d th run contains %6d members"%(jptp,newindeces[m][1],\
								len(indep_runs_to_groups[jptp][newindeces[m][1]]))) 
			if myid==main_node and iptp<jptp:
				unaccounted=total_stack-tt
				ratio_unaccounted = 100.-tt/total_stack*100.
				ratio_accounted   = tt/total_stack*100
				aline           = print_a_line_with_timestamp("Accounted data is %6d, %5.2f "%(int(tt),ratio_accounted))
				log_main.add(aline)
				aline=print_a_line_with_timestamp("Unaccounted data is %6d, %5.2f"%(int(unaccounted),ratio_unaccounted))
				log_main.add(aline)
                       	rate=tt/total_stack*100.0
			scores[(iptp,jptp)] =rate
			if iptp<jptp:
				avg_two_ways 	    +=rate
				avg_two_ways_square +=rate**2
				total_pop +=1
                       	if myid ==main_node and iptp<jptp:
                               	aline=print_a_line_with_timestamp("The two-way comparison stable member total ratio %3d %3d %5.3f  "%(iptp,jptp,rate))
				log_main.add(aline)
                       	new_list=[]
                       	for a in list_stable:
                       		a.tolist()
                                new_list.append(a)
                        two_ways_stable_member_list[(iptp,jptp)]=new_list
	if myid ==main_node:
		log_main.add("two_way comparison is done!")
	#### Score each independent run by pairwise summation
	summed_scores =[]
	two_way_dict  ={}
	for ipp in xrange(len(ptp)):
		avg_scores =0.0
		for jpp in xrange(len(ptp)):
			if ipp!=jpp:
				avg_scores +=scores[(ipp,jpp)]
		avg_rate =avg_scores/(len(ptp)-1)
		summed_scores.append(avg_rate)
		two_way_dict[avg_rate] =ipp
	#### Select two independent runs that have the first two highest scores
	if Tracker["constants"]["indep_runs"]>2:
		rate1=max(summed_scores)
		run1 =two_way_dict[rate1]
		summed_scores.remove(rate1)
		rate2 =max(summed_scores)
		run2 = two_way_dict[rate2]
	else:
		run1 =0
		run2 =1
		rate1  =max(summed_scores)
		rate2  =rate1
	Tracker["two_way_stable_member"]      =two_ways_stable_member_list[(run1,run2)]
	Tracker["pop_size_of_stable_members"] =1
	if myid ==main_node:
                log_main.add("Get outliers of the selected comparison")
	####  Save accounted ones and unaccounted ones
	counted_list = merge_groups(two_ways_stable_member_list[(run1,run2)])
	outliers     = get_outliers(total_stack,counted_list)
	if myid ==main_node:
                log_main.add("Save outliers")
	save_alist(Tracker,"Unaccounted.txt",outliers)
	mpi_barrier(MPI_COMM_WORLD)
	save_alist(Tracker,"Accounted.txt",counted_list)
	mpi_barrier(MPI_COMM_WORLD)
	Tracker["this_unaccounted_dir"]=main_dir
	Tracker["this_unaccounted"]    =os.path.join(main_dir,"Unaccounted.txt")
	Tracker["this_accounted"]      =os.path.join(main_dir,"Accounted.txt")
	Tracker["ali3d_of_outliers"] =os.path.join(main_dir,"ali3d_params_of_outliers.txt")
	Tracker["ali3d_of_accounted"]  =os.path.join(main_dir, "ali3d_params_of_accounted.txt")
	if myid==main_node:
		log_main.add(" Selected indepedent runs      %5d and  %5d"%(run1,run2))
		log_main.add(" Their pair-wise averaged rates are %5.2f  and %5.2f "%(rate1,rate2))		
	from math import sqrt
	avg_two_ways = avg_two_ways/total_pop
	two_ways_std = sqrt(avg_two_ways_square/total_pop-avg_two_ways**2)
	net_rate     = avg_two_ways-1./number_of_groups*100.
	Tracker["net_rate"]=net_rate
	if myid ==main_node: 
		msg="average of two-way comparison  %5.3f"%avg_two_ways
		log_main.add(msg)
		msg="net rate of two-way comparison  %5.3f"%net_rate
		log_main.add(msg)
		msg="std of two-way comparison %5.3f"%two_ways_std
		log_main.add(msg)
		msg ="Score table of two_way comparison when Kgroup =  %5d"%number_of_groups
		log_main.add(msg)
		print_upper_triangular_matrix(scores,Tracker["constants"]["indep_runs"],log_main)
	del two_ways_stable_member_list
	Tracker["score_of_this_comparison"]=(avg_two_ways,two_ways_std,net_rate)
	mpi_barrier(MPI_COMM_WORLD)
	
def do_EKmref(Tracker):
	from mpi import mpi_barrier, MPI_COMM_WORLD
	from utilities import cmdexecute
	## Use stable members of only one two-way comparison
	main_dir  =Tracker["this_dir"]
	pop_size  =Tracker["pop_size_of_stable_members"]# Currently is 1, can be changed
	log_main  =Tracker["constants"]["log_main"]
	myid      =Tracker["constants"]["myid"]
	main_node =Tracker["constants"]["main_node"]
	import os
	###
	EKMREF_iteration ={}
	if myid ==main_node:
		log_main.add("-----------------EKmref---------------------------")
	for inew in xrange(1):# new independent run
        	output_for_this_run ={}
                output_for_this_run["output_dir"]=os.path.join(main_dir,"EKMREF%03d"%inew)
                output_for_this_run["partition"] =os.path.join(output_for_this_run["output_dir"],"list2.txt")
                output_for_this_run["refvols"]=os.path.join(main_dir,"vol_stable_member_%03d.hdf"%inew)
                EKMREF_iteration[inew]=output_for_this_run
        Tracker["EKMREF"]=EKMREF_iteration
        for inew in xrange(1):
        	if myid ==main_node:
                	msg="Create %dth two-way reference  model"%inew
                	log_main.add(msg)
                name_of_grouped_vols=Tracker["EKMREF"][inew]["refvols"]
                #Tracker["this_data_stack"] =Tracker["constants"]["stack"]
                grouped_plist=Tracker["two_way_stable_member"]
                reconstruct_grouped_vols(Tracker,name_of_grouped_vols,grouped_plist)
	########
	Tracker["this_iruns"]=1
	Tracker["data_stack_of_accounted"]="bdb:"+os.path.join(Tracker["this_dir"],"accounted")
	partition_ali3d_params_of_outliers(Tracker)
	if myid==main_node:
		cmd = "{} {} {} {} {}".format("e2bdb.py",Tracker["this_data_stack"],"--makevstack",\
			Tracker["data_stack_of_accounted"],"--list="+Tracker["this_accounted"])
       		cmdexecute(cmd)
		#cmd = "{} {} {} {}  {}".format("sxheader.py", Tracker["data_stack_of_counted"],"--params=xform.projection",\
		#		 "--import="+Tracker["ali3d_of_counted"], "--consecutive ")
		#cmdexecute(cmd)
	mpi_barrier(MPI_COMM_WORLD)
	Tracker["this_data_stack"]=Tracker["data_stack_of_accounted"]
	do_independent_EKmref(Tracker)

def Kgroup_guess(Tracker,data_stack):
	myid      = Tracker["constants"]["myid"]
	main_node = Tracker["constants"]["main_node"]
	log_main  = Tracker["constants"]["log_main"]
	if myid==main_node:
		msg="Program Kgroup_guess"
		log_main.add(msg)
		msg="Estimate how many groups the dataset can be partitioned into"
		log_main.add(msg)
	number_of_groups=2
	Tracker["this_data_stack"]=data_stack
	maindir                   =Tracker["this_dir"]
	terminate                 =0
	while terminate !=1:
		do_EKTC_for_a_given_number_of_groups(Tracker,number_of_groups)
		(score_this, score_std_this, score_net)=Tracker["score_of_this_comparison"]
		number_of_groups +=1
		if number_of_groups >10:
			terminate  =1
def do_N_groups(Tracker,data_stack):
	myid      = Tracker["constants"]["myid"]
        main_node = Tracker["constants"]["main_node"]
        log_main  = Tracker["constants"]["log_main"]
        if myid==main_node:
                msg="Program do_N_groups"
                log_main.add(msg)
                msg="Estimate how many groups the dataset can be partitioned into"
                log_main.add(msg)
        Tracker["this_data_stack"]=data_stack
        maindir                   =Tracker["this_dir"]
	number_of_groups          =Tracker["number_of_groups"]
	for Ngroup in xrange(2,number_of_groups+1):
		Tracker["this_data_stack"]=data_stack
		Tracker["this_dir"]       =maindir
		Tracker["number_of_groups"]=Ngroup
		do_EKTC_for_a_given_number_of_groups(Tracker,Ngroup)
                (score_this, score_std_this, score_net)=Tracker["score_of_this_comparison"]
		Tracker["this_dir"]       =os.path.join(maindir,"GRP%03d"%Tracker["number_of_groups"])
		do_EKmref(Tracker)

def do_EKTC_for_a_given_number_of_groups(Tracker,given_number_of_groups):
	import os
	from utilities import cmdexecute
	Tracker["number_of_groups"]= given_number_of_groups
	from mpi import mpi_barrier,MPI_COMM_WORLD
	log_main  = Tracker["constants"]["log_main"]
	myid      = Tracker["constants"]["myid"]
	main_node = Tracker["constants"]["main_node"]
	maindir   = Tracker["this_dir"]
	workdir   = os.path.join(maindir,"GRP%03d"%Tracker["number_of_groups"])
	if myid ==main_node:
		msg ="Number of group is %5d"%given_number_of_groups
		log_main.add(msg)
 		cmd="{} {}".format("mkdir",workdir)
        	cmdexecute(cmd)
	mpi_barrier(MPI_COMM_WORLD)
	Tracker["this_dir"]=workdir
	N_independent_reconstructions(Tracker)
	prepare_EMREF_dict(Tracker)
	N_independent_mref(Tracker["constants"]["mpi_comm"],Tracker)
	do_two_way_comparison(Tracker)
	Tracker["this_dir"]=maindir # always reset maindir

def do_EKTC_for_a_large_number_RUN(Tracker):
        import os
        from utilities import cmdexecute
        Tracker["number_of_groups"]=Tracker["constants"]["Kgroup"] 
        from mpi import mpi_barrier,MPI_COMM_WORLD
        log_main  = Tracker["constants"]["log_main"]
        myid      = Tracker["constants"]["myid"]
        main_node = Tracker["constants"]["main_node"]
        N_independent_reconstructions(Tracker)
        prepare_EMREF_dict(Tracker)
        N_independent_mref(Tracker["constants"]["mpi_comm"],Tracker)
        do_two_way_comparison(Tracker)

def low_pass_filter_search(Tracker):
        import os
        from utilities import cmdexecute,write_text_file
        from mpi import mpi_barrier,MPI_COMM_WORLD
	from logger import Logger, BaseLogger_Files
        log_main                   = Tracker["constants"]["log_main"]
        myid                       = Tracker["constants"]["myid"]
        main_node                  = Tracker["constants"]["main_node"]
        maindir                    = Tracker["this_dir"]
	number_of_groups           = Tracker["number_of_groups"]
	frequency_start_search     = Tracker["constants"]["frequency_start_search"]
	frequency_stop_search      = Tracker["constants"]["frequency_stop_search"]
	frequency_search_step      = Tracker["constants"]["frequency_search_step"]
	N_steps                    = int((frequency_stop_search-frequency_start_search)/frequency_search_step)
	search_result = []
	log_lpfsearch       =Logger(BaseLogger_Files())
        log_lpfsearch.prefix=maindir+"/"
	if myid ==main_node:
		log_lpfsearch.add(" program low_pass_filter_search  ")
		log_lpfsearch.add(" start frequency is %5.3f"%frequency_start_search)
		log_lpfsearch.add(" stop frequency is %5.3f"%frequency_stop_search)
		log_lpfsearch.add(" search step is %5.3f"%frequency_search_step)
		log_lpfsearch.add(" total iteration is %5d"%N_steps)
	for iter_lpf in xrange(N_steps):
		frequency_low_pass           =frequency_start_search + iter_lpf*frequency_search_step
		Tracker["frequency_low_pass"]=frequency_low_pass 
		workdir                      =os.path.join(maindir,"lpf"+str(frequency_low_pass)) 
        	if myid ==main_node:
                	msg ="the low-pass filter applied frequency is %5.3f"%frequency_low_pass
                	log_lpfsearch.add(msg)
                	cmd="{} {}".format("mkdir",workdir)
                	cmdexecute(cmd)
        	mpi_barrier(MPI_COMM_WORLD)
        	Tracker["this_dir"]=workdir
        	N_independent_reconstructions(Tracker)
        	prepare_EMREF_dict(Tracker)
        	N_independent_mref(Tracker["constants"]["mpi_comm"],Tracker)
        	do_two_way_comparison(Tracker)
		net_rate = Tracker["net_rate"]
		search_result.append([net_rate,frequency_low_pass])
		if myid ==main_node:
                	msg =" frequency  %5.3f   net rate  %5.3f "%(frequency_low_pass,net_rate)
                        log_lpfsearch.add(msg)
		mpi_barrier(MPI_COMM_WORLD)
	if myid ==main_node:
		write_text_file(search_result,os.path.join(maindir, "lpf_group%03d.txt"%number_of_groups))
def main():
	from logger import Logger, BaseLogger_Files
        arglist = []
        i = 0
        while( i < len(sys.argv) ):
            if sys.argv[i]=='-p4pg':
                i = i+2
            elif sys.argv[i]=='-p4wd':
                i = i+2
            else:
                arglist.append( sys.argv[i] )
                i = i+1
	progname = os.path.basename(arglist[0])
	usage = progname + " stack  outdir refvols  <mask> --focus=3Dmask --ir=inner_radius --radius=outer_radius --rs=ring_step --xr=x_range --yr=y_range  --ts=translational_searching_step " +\
	" --delta=angular_step --an=angular_neighborhood --center=1 --nassign=reassignment_number --nrefine=alignment_number --maxit=max_iter --stoprnct=percentage_to_stop " + \
	" --debug --fourvar=fourier_variance --CTF --snr=1.0 --ref_a=S --sym=c1 --function=user_function --independent=indenpendent_runs_for_equalmref  --Kgroup=number_of_groups  --resolution --mode=EK --import=params.txt"
	parser = OptionParser(usage,version=SPARXVERSION)
	parser.add_option("--focus",    type="string",       default=None,             help="3D mask for focused clustering ")
	parser.add_option("--ir",       type= "int",         default=1, 	       help="inner radius for rotational correlation > 0 (set to 1)")
	parser.add_option("--radius",   type= "int",         default="-1",	       help="outer radius for rotational correlation <nx-1 (set to the radius of the particle)")
	parser.add_option("--maxit",	type= "int",         default=5, 	       help="maximum number of iteration")
	parser.add_option("--rs",       type= "int",         default="1",	       help="step between rings in rotational correlation >0 (set to 1)" ) 
	parser.add_option("--xr",       type="string",       default="4 2 1 1 1",      help="range for translation search in x direction, search is +/-xr ")
	parser.add_option("--yr",       type="string",       default="-1",	       help="range for translation search in y direction, search is +/-yr (default = same as xr)")
	parser.add_option("--ts",       type="string",       default="0.25",           help="step size of the translation search in both directions direction, search is -xr, -xr+ts, 0, xr-ts, xr ")
	parser.add_option("--delta",    type="string",       default="10 6 4  3   2",  help="angular step of reference projections")
	parser.add_option("--an",       type="string",       default="-1",	       help="angular neighborhood for local searches")
	parser.add_option("--center",   type="int",          default=0,	               help="0 - if you do not want the volume to be centered, 1 - center the volume using cog (default=0)")
	parser.add_option("--nassign",  type="int",          default=0, 	       help="number of reassignment iterations performed for each angular step (set to 3) ")
	parser.add_option("--nrefine",  type="int",          default=1, 	       help="number of alignment iterations performed for each angular step (set to 1) ")
	parser.add_option("--CTF",      action="store_true", default=False,            help="Consider CTF correction during the alignment ")
	parser.add_option("--snr",      type="float",        default=1.0,              help="Signal-to-Noise Ratio of the data")   
	parser.add_option("--stoprnct", type="float",        default=0.0,              help="Minimum percentage of assignment change to stop the program")   
	parser.add_option("--ref_a",    type="string",       default="S",              help="method for generating the quasi-uniformly distributed projection directions (default S) ")
	parser.add_option("--sym",      type="string",       default="c1",             help="symmetry of the structure ")
	parser.add_option("--function", type="string",       default="ref_ali3dm",     help="name of the reference preparation function")
	parser.add_option("--MPI",      action="store_true", default=False,            help="Use MPI version ")
	parser.add_option("--npad",     type="int",          default= 2,               help="padding size for 3D reconstruction")
	parser.add_option("--debug",    action="store_true", default=False,            help="debug ")
	parser.add_option("--fourvar",  action="store_true", default=False,            help="compute and use fourier variance")
	parser.add_option("--independent", type="int",       default= 3,               help="number of independent run")
	parser.add_option("--Kgroup",      type="int",       default= 5,               help="number of groups")
	parser.add_option("--resolution",  type="float",     default= .40,             help="structure is low-pass-filtered to this resolution for clustering" )
	parser.add_option("--mode",        type="string",    default="EK_only",        help="mode options: EK_only, Kgroup_guess, auto_search, lpf_search,large_number_run" )
	parser.add_option("--importali3d", type="string",    default="",               help="import the xform.projection parameters as the initial configuration for 3-D reconstruction" )
	parser.add_option("--do_unaccounted",  action="store_true",default=False,        help="continue clustering on unaccounted images" )
	parser.add_option("--Kgroup_guess",  action="store_true",default=False,        help="Guess the possible number of groups existing in one dataset" )
	parser.add_option("--frequency_start_search",  type="float",default=.10,       help="start frequency for low pass filter search")
	parser.add_option("--frequency_stop_search",   type="float",default=.40,       help="stop frequency for low pass filter search")
	parser.add_option("--frequency_search_step",   type="float",default=.02,       help="frequency step for low pass filter search")
	parser.add_option("--members_per_group", type="int",default=1000,             help="minimum number of images per group")
	(options, args) = parser.parse_args(arglist[1:])
	if len(args) < 1  or len(args) > 4:
    		print "usage: " + usage
    		print "Please run '" + progname + " -h' for detailed options"
	else:

		if global_def.CACHE_DISABLE:
			from utilities import disable_bdb_cache
			disable_bdb_cache()
            	if len(args)>2:
            		mask_file = args[2]
            	else:
                	mask_file = None

		orgstack                        =args[0]
		masterdir                       =args[1]
		global_def.BATCH = True
		#---initialize MPI related variables
		from mpi import mpi_init, mpi_comm_size, MPI_COMM_WORLD, mpi_comm_rank,mpi_barrier,mpi_bcast, mpi_bcast, MPI_INT
		sys.argv = mpi_init(len(sys.argv),sys.argv)
		nproc     = mpi_comm_size(MPI_COMM_WORLD)
		myid      = mpi_comm_rank(MPI_COMM_WORLD)
		mpi_comm = MPI_COMM_WORLD
		main_node = 0
		# import some utilites
		from utilities import get_im,bcast_number_to_all,cmdexecute,write_text_file,read_text_file  #get_shrink_data
		from applications import recons3d_n_MPI, mref_ali3d_MPI, Kmref_ali3d_MPI
		from statistics import k_means_match_clusters_asg_new,k_means_stab_bbenum 
		# Create the main log file
	        from logger import Logger,BaseLogger_Files
             	if myid ==main_node:
             		log_main=Logger(BaseLogger_Files())
                     	log_main.prefix=masterdir+"/"
                else:
                        log_main =None
		#--- fill input parameters into dictionary named after Constants
		Constants		         ={}
		Constants["stack"]               = args[0]
        	Constants["masterdir"]           = masterdir
		Constants["mask3D"]              = mask_file
		Constants["focus3Dmask"]         = options.focus
		Constants["indep_runs"]          = options.independent
		Constants["stoprnct"]            = options.stoprnct
		Constants["Kgroup"]              = options.Kgroup
		Constants["CTF"]                 = options.CTF
		Constants["npad"]                = options.npad
		Constants["maxit"]               = options.maxit
		Constants["ir"]                  = options.ir 
		Constants["radius"]              = options.radius 
		Constants["nassign"]             = options.nassign
		Constants["snr"]                 = options.snr
		Constants["rs"]                  = options.rs 
		Constants["xr"]                  = options.xr
		Constants["yr"]                  = options.yr
		Constants["ts"]                  = options.ts
		Constants["ref_a"]               = options.ref_a
		Constants["delta"]               = options.delta
		Constants["an"]                  = options.an
		Constants["sym"]                 = options.sym
		Constants["center"]              = options.center
		Constants["nrefine"]             = options.nrefine
		Constants["fourvar"]             = options.fourvar 
		Constants["user_func"]           = options.function
		Constants["resolution"]          = options.resolution
		Constants["debug"]               = options.debug
		Constants["sign"]                = 1
		Constants["listfile"]            = ""
		Constants["xysize"]              =-1
		Constants["zsize"]               =-1
		Constants["group"]               =-1
		Constants["verbose"]             = 0
		Constants["mode"]                = options.mode
		Constants["mpi_comm"]            = mpi_comm
		Constants["main_log_prefix"]     =args[1]
		Constants["importali3d"]         =options.importali3d
		Constants["myid"]	         =myid
		Constants["main_node"]           =main_node
		Constants["log_main"]            =log_main
		Constants["do_unaccounted"]        =options.do_unaccounted
		Constants["frequency_search_step"] = options.frequency_search_step
		Constants["frequency_start_search"] = options.frequency_start_search
		Constants["frequency_stop_search"] = options.frequency_stop_search
		Constants["members_per_group"]    = options.members_per_group
		# -----------------------------------------------------
		#
		# Create and initialize Tracker dictionary with input options
		Tracker = 			    		{}
		Tracker["constants"]=				Constants
		Tracker["maxit"]          = Tracker["constants"]["maxit"]
        	Tracker["radius"]         = Tracker["constants"]["radius"]
        	Tracker["xr"]             = ""
        	Tracker["yr"]             = "-1"  # Do not change!
        	Tracker["ts"]             = 1
        	Tracker["an"]             = "-1"
        	Tracker["delta"]          = "2.0"
        	Tracker["zoom"]           = True
        	Tracker["nsoft"]          = 0
        	Tracker["local"]          = False
        	Tracker["PWadjustment"]   = ""
        	Tracker["upscale"]        = 0.5
        	Tracker["applyctf"]       = True  #  Should the data be premultiplied by the CTF.  Set to False for local continuous.
        	#Tracker["refvol"]         = None
        	Tracker["nxinit"]         = 64
        	Tracker["nxstep"]         = 32
        	Tracker["icurrentres"]    = -1
        	Tracker["ireachedres"]    = -1
        	Tracker["lowpass"]        = 0.4
        	Tracker["falloff"]        = 0.2
        	#Tracker["inires"]         = options.inires  # Now in A, convert to absolute before using
        	Tracker["fuse_freq"]      = 50  # Now in A, convert to absolute before using
        	Tracker["delpreviousmax"] = False
        	Tracker["anger"]          = -1.0
        	Tracker["shifter"]        = -1.0
        	Tracker["saturatecrit"]   = 0.95
        	Tracker["pixercutoff"]    = 2.0
        	Tracker["directory"]      = ""
        	Tracker["previousoutputdir"] = ""
        	#Tracker["eliminated-outliers"] = False
        	Tracker["mainiteration"]  = 0
        	Tracker["movedback"]      = False
        	#Tracker["state"]          = Tracker["constants"]["states"][0] 
		Tracker["global_resolution"] =0.0
		#--------------------------------------------------------------------
		#
		# Get the pixel size; if none, set to 1.0, and the original image size
        	if(myid == main_node):
                	line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
                	print(line,"INITIALIZATION OF Equal-Kmeans & Kmeans Clustering")

                	a = get_im(orgstack)
                	nnxo = a.get_xsize()
                	if( Tracker["nxinit"] > nnxo ):
                        	ERROR("Image size less than minimum permitted $d"%Tracker["nxinit"],"sxmeridien",1)
                        	nnxo = -1
                	else:
                        	if Tracker["constants"]["CTF"]:
                                	i = a.get_attr('ctf')
                                	pixel_size = i.apix
                                	fq = pixel_size/Tracker["fuse_freq"]
                        	else:
                                	pixel_size = 1.0
                                	#  No pixel size, fusing computed as 5 Fourier pixels
                                	fq = 5.0/nnxo
                        	del a
        	else:
                	nnxo = 0
                	fq = 0.0
                	pixel_size = 1.0
        	nnxo = bcast_number_to_all(nnxo, source_node = main_node)
        	if( nnxo < 0 ):
                	mpi_finalize()
                	exit()
        	pixel_size = bcast_number_to_all(pixel_size, source_node = main_node)
        	fq         = bcast_number_to_all(fq, source_node = main_node)
        	Tracker["constants"]["nnxo"]         = nnxo
        	Tracker["constants"]["pixel_size"]   = pixel_size
        	Tracker["fuse_freq"]    = fq
        	del fq, nnxo, pixel_size

        	if(Tracker["constants"]["radius"]  < 1):
                	Tracker["constants"]["radius"]  = Tracker["constants"]["nnxo"]//2-2
        	elif((2*Tracker["constants"]["radius"] +2) > Tracker["constants"]["nnxo"]):
                	ERROR("Particle radius set too large!","sxEKmref_clustering.py",1,myid)

####-----------------------------------------------------------------------------------------
		# Master directory
		if myid == main_node:
			if masterdir =="":
				timestring = strftime("_%d_%b_%Y_%H_%M_%S", localtime())
				masterdir ="master_EKMREF"+timestring
				li =len(masterdir)
				cmd="{} {}".format("mkdir", masterdir)
				cmdexecute(cmd)
			else:
				li = 0
				keepchecking =1			
		else:
			li=0
			keepchecking =1
	    	li = mpi_bcast(li,1,MPI_INT,main_node,MPI_COMM_WORLD)[0]
		if li>0:
			masterdir = mpi_bcast(masterdir,li,MPI_CHAR,main_node,MPI_COMM_WORLD)
			masterdir = string.join(masterdir,"")
		if myid ==main_node:
			print_dict(Tracker["constants"], "Permanent settings of Equal-Kmeans & Kmeans Clustering")
		######### create a vstack from input stack to the local stack in masterdir
		# stack name set to default
		Tracker["constants"]["stack"] = "bdb:"+masterdir+"/rdata"
	   	if(myid == main_node):
                	if keepchecking:
                        	if(os.path.exists(os.path.join(masterdir,"EMAN2DB/rdata.bdb"))):  doit = False
                        	else:  doit = True
                	else:  doit = True
                	if  doit:
                        	if(orgstack[:4] == "bdb:"):     cmd = "{} {} {}".format("e2bdb.py", orgstack,"--makevstack="+Tracker["constants"]["stack"])
                        	else:  cmd = "{} {} {}".format("sxcpy.py", orgstack, Tracker["constants"]["stack"])
                        	cmdexecute(cmd)
                        	cmd = "{} {}".format("sxheader.py  --consecutive  --params=originalid", Tracker["constants"]["stack"])
                        	cmdexecute(cmd)
                        	keepchecking = False
                	total_stack = EMUtil.get_image_count(Tracker["constants"]["stack"])
        	else:
                	total_stack = 0
		mpi_barrier(MPI_COMM_WORLD)
        	total_stack = bcast_number_to_all(total_stack, source_node = main_node)
		
		###----------------------------------------------------------------------------------
		# Initial data analysis 
		from random import shuffle
		# Compute the resolution 
		partids =[None]*2
		for procid in xrange(2):partids[procid] =os.path.join(masterdir, "chunk%01d.txt"%procid)
		inivol =[None]*2
		for procid in xrange(2):inivol[procid]  =os.path.join(masterdir,"vol%01d.hdf"%procid)
					
		if myid ==main_node:
			ll=range(total_stack)
			shuffle(ll)
			l1=ll[0:total_stack//2]
			l2=ll[total_stack//2:]
			del ll
			l1.sort()
			l2.sort()
			write_text_file(l1,partids[0])
			write_text_file(l2,partids[1])
		mpi_barrier(MPI_COMM_WORLD)
		ll =[]
		for procid in xrange(2):
			l=read_text_file(partids[procid])
			ll.append(l)
		### create two volumes to estimate resolution
		vols=[]
		for procid in xrange(2):
			recons3d_n_MPI(Tracker["constants"]["stack"],ll[procid],inivol[procid],Tracker["constants"]["CTF"],\
			Tracker["constants"]["snr"],Tracker["constants"]["sign"],Tracker["constants"]["npad"],\
			Tracker["constants"]["sym"],Tracker["constants"]["listfile"],Tracker["constants"]["group"],\
			Tracker["constants"]["verbose"],Tracker["constants"]["xysize"],Tracker["constants"]["zsize"])
		mpi_barrier(MPI_COMM_WORLD)
		for procid in xrange(2):
			if myid==main_node:
				vols.append(get_im(inivol[procid]))
			else:
				vols.append(None)
		if myid ==main_node:
			low_pass, falloff,currentres =get_resolution_mrk01(vols,Tracker["constants"]["radius"],\
			Tracker["constants"]["nnxo"],masterdir,Tracker["constants"]["mask3D"])
			if low_pass >Tracker["constants"]["resolution"]:
				low_pass= Tracker["constants"]["resolution"]
		else:
			low_pass    =0.0
			falloff     =0.0
			currentres  =0.0
		bcast_number_to_all(currentres,source_node = main_node)
		bcast_number_to_all(low_pass,source_node = main_node)
		bcast_number_to_all(falloff,source_node = main_node)
		Tracker["currentres"]         = currentres
		Tracker["low_pass"]           = low_pass
		Tracker["falloff"]            = falloff
		Tracker["frequency_low_pass"] =.4
		if myid ==main_node:
			log_main.add("The command-line inputs are as following:")
			log_main.add("**********************************************************")
		for a in sys.argv:
			if myid ==main_node:
				log_main.add(a)
			else:
				continue
		if myid ==main_node:
			log_main.add("**********************************************************")
		### START EK Clustering
		if myid ==main_node:
                        log_main.add("----------EK Clustering------- ")
		if Tracker["constants"]["mode"]=="Kgroup_guess":
			if myid ==main_node:
                        	log_main.add("------Current mode is Kgroup_guess------")
			Tracker["this_dir"] =os.path.join(masterdir,"Kgroup_guess")
			if myid ==main_node:
                        	cmd="{} {}".format("mkdir",Tracker["this_dir"])
                        	cmdexecute(cmd)
                	Tracker["this_data_stack"]  ="bdb:"+os.path.join(Tracker["this_dir"],"data")
                	Tracker["this_total_stack"] =total_stack
                	Tracker["number_of_groups"] =Tracker["constants"]["Kgroup"]
                	Tracker["this_ali3d"]       =Tracker["constants"]["importali3d"]
                	Tracker["importali3d"]      =Tracker["constants"]["importali3d"]
	                ### Create data stack
         		if myid ==main_node:
                        	log_main.add("----create data stack ")
                        	cmd = "{} {} {} {} ".format("e2bdb.py",orgstack,"--makevstack",\
                        		Tracker["this_data_stack"])
                        	cmdexecute(cmd)
                        	cmd = "{} {}".format("sxheader.py  --consecutive  --params=originalid", Tracker["this_data_stack"])
                        	cmdexecute(cmd)
                        	if Tracker["importali3d"] !="":
                                	cmd = "{} {} {} {} ".format("sxheader.py", Tracker["this_data_stack"],"--params=xform.projection",\
                                 	"--import="+Tracker["importali3d"])
                        	cmdexecute(cmd)
                	mpi_barrier(MPI_COMM_WORLD)
			Kgroup_guess(Tracker,Tracker["this_data_stack"])
			if myid ==main_node:
				msg ="Kgroup guessing is done!"
				log_main.add(msg)
		elif Tracker["constants"]["mode"]=="lpf_search":
			if myid==main_node:
                                msg ="--------Current mode is lpf_search: search for stringent low-pass filter--------"
                                log_main.add(msg)
			Tracker["this_dir"] =os.path.join(masterdir,"lpf_search")
			if myid ==main_node:
					cmd="{} {}".format("mkdir",Tracker["this_dir"])
                                	cmdexecute(cmd)
		        Tracker["this_data_stack"]  ="bdb:"+os.path.join(Tracker["this_dir"],"data")
                        Tracker["this_total_stack"] =total_stack
                        Tracker["number_of_groups"] =Tracker["constants"]["Kgroup"]
                        Tracker["this_ali3d"]       =Tracker["constants"]["importali3d"]
                        Tracker["importali3d"]      =Tracker["constants"]["importali3d"]
			if myid ==main_node:
                        	log_main.add("----create data stack ")
                                cmd = "{} {} {} {} ".format("e2bdb.py",orgstack,"--makevstack",\
                                        Tracker["this_data_stack"])
                                cmdexecute(cmd)
                                cmd = "{} {}".format("sxheader.py  --consecutive  --params=originalid", Tracker["this_data_stack"])
                                cmdexecute(cmd)
                                if Tracker["importali3d"] !="":
                                        cmd = "{} {} {} {} ".format("sxheader.py", Tracker["this_data_stack"],"--params=xform.projection",\
                                        "--import="+Tracker["importali3d"])
                                cmdexecute(cmd)
                        mpi_barrier(MPI_COMM_WORLD)
                        low_pass_filter_search(Tracker)
                        if myid ==main_node:
                                msg ="lpf_search is done!"
                                log_main.add(msg)
		elif Tracker["constants"]["mode"]=="large_number_run":
                        if myid==main_node:
                                msg ="--------Current mode is large number run--------"
                                log_main.add(msg)
			Tracker["this_dir"] =os.path.join(masterdir,"large_number_run")
	             	if myid ==main_node:
                                cmd="{} {}".format("mkdir",Tracker["this_dir"])
                                cmdexecute(cmd)
			Tracker["this_data_stack"]  ="bdb:"+os.path.join(Tracker["this_dir"],"data")
                        Tracker["this_total_stack"] =total_stack
                        Tracker["number_of_groups"] =Tracker["constants"]["Kgroup"]
                        Tracker["this_ali3d"]       =Tracker["constants"]["importali3d"]
                        Tracker["importali3d"]      =Tracker["constants"]["importali3d"]
			### Create data stack
			if myid ==main_node:
                        	log_main.add("----Create data stack ")
                                cmd = "{} {} {} {} ".format("e2bdb.py",orgstack,"--makevstack",\
                                Tracker["this_data_stack"])
                                cmdexecute(cmd)
                                cmd = "{} {}".format("sxheader.py  --consecutive  --params=originalid", Tracker["this_data_stack"])
                                cmdexecute(cmd)
                       		if Tracker["importali3d"] !="":
                                	cmd = "{} {} {} {} ".format("sxheader.py", Tracker["this_data_stack"],"--params=xform.projection",\
                                	 "--import="+Tracker["importali3d"])
                        		cmdexecute(cmd)
                        mpi_barrier(MPI_COMM_WORLD)
			do_EKTC_for_a_large_number_RUN(Tracker)
		elif Tracker["constants"]["mode"]=="EK_only":
			if myid==main_node:
				msg ="--------Current mode is EK_only--------"
                                log_main.add(msg)
			Tracker["this_dir"] =os.path.join(masterdir,"EK_only")
	             	if myid ==main_node:
                                cmd="{} {}".format("mkdir",Tracker["this_dir"])
                                cmdexecute(cmd)
			Tracker["this_data_stack"]  ="bdb:"+os.path.join(Tracker["this_dir"],"data")
                        Tracker["this_total_stack"] =total_stack
                        Tracker["number_of_groups"] =Tracker["constants"]["Kgroup"]
                        Tracker["this_ali3d"]       =Tracker["constants"]["importali3d"]
                        Tracker["importali3d"]      =Tracker["constants"]["importali3d"]
			### Create data stack
			if myid ==main_node:
                        	log_main.add("----Create data stack ")
                                cmd = "{} {} {} {} ".format("e2bdb.py",orgstack,"--makevstack",\
                                Tracker["this_data_stack"])
                                cmdexecute(cmd)
                                cmd = "{} {}".format("sxheader.py  --consecutive  --params=originalid", Tracker["this_data_stack"])
                                cmdexecute(cmd)
                       		if Tracker["importali3d"] !="":
                                	cmd = "{} {} {} {} ".format("sxheader.py", Tracker["this_data_stack"],"--params=xform.projection",\
                                	 "--import="+Tracker["importali3d"])
                        		cmdexecute(cmd)
                        mpi_barrier(MPI_COMM_WORLD)
			N_independent_reconstructions(Tracker)
			prepare_EMREF_dict(Tracker)
                	N_independent_mref(mpi_comm,Tracker)
                	do_two_way_comparison(Tracker)
                	mpi_barrier(MPI_COMM_WORLD)
			do_EKmref(Tracker)
			if myid ==main_node:
                                msg ="Equal-kmeans and K-means is done!"
                                log_main.add(msg)
			if Tracker["constants"]["do_uncounted"]:
				if myid ==main_node:
					msg ="----Continue EK clustering on those uncounted members-----"
					log_main.add(msg)
                                Tracker["this_dir"]   =os.path.join(Tracker["this_unaccounted_dir"],"Unaccounted")
				Tracker["last_data_stack"] =Tracker["constants"]["stack"]
				Tracker["this_data_stack"] ="bdb:"+os.path.join(Tracker["this_dir"],"data")
				Tracker["number_of_groups"]=Tracker["constants"]["Kgroup"]
				if myid ==main_node:
                        		cmd="{} {}".format("mkdir",Tracker["this_dir"])
                        		cmdexecute(cmd)
					cmd = "{} {} {} {} {}".format("e2bdb.py",Tracker["last_data_stack"],"--makevstack",\
                                        Tracker["this_data_stack"],"--list="+Tracker["this_unaccounted"])
                                        cmdexecute(cmd)
					#cmd = "{} {} {} {}  {}".format("sxheader.py", Tracker["this_data_stack"],"--params=xform.projection",\
                                        #      "--import="+Tracker["ali3d_of_outliers"], "--consecutive")
                                        #cmdexecute(cmd)
				mpi_barrier(MPI_COMM_WORLD)
				Tracker["this_total_stack"] =Tracker["number_of_unaccounted"]
				do_N_groups(Tracker,Tracker["this_data_stack"])
                        	mpi_barrier(MPI_COMM_WORLD)
		elif Tracker["constants"]["mode"]=="auto_search":
			generation =0
			if myid ==main_node:
				log_main.add("Current mode is auto_search")
                        	log_main.add("clustering generation %3d"%generation)
                	Tracker["this_dir"]   =os.path.join(masterdir,"generation%03d"%generation)
                	if myid ==main_node:
                        	cmd="{} {}".format("mkdir",Tracker["this_dir"])
                        	cmdexecute(cmd)
                	Tracker["this_data_stack"]  ="bdb:"+os.path.join(Tracker["this_dir"],"data")
                	Tracker["this_total_stack"] =total_stack
                	Tracker["number_of_groups"] =Tracker["constants"]["Kgroup"]
                	Tracker["this_ali3d"]       =Tracker["constants"]["importali3d"]
                	Tracker["importali3d"]      =Tracker["constants"]["importali3d"]
			### Create data stack
			if myid ==main_node:
                        	log_main.add("----Create data stack ")
                                cmd = "{} {} {} {} ".format("e2bdb.py",orgstack,"--makevstack",\
                                        Tracker["this_data_stack"])
                                cmdexecute(cmd)
                                cmd = "{} {}".format("sxheader.py  --consecutive  --params=originalid", Tracker["this_data_stack"])
                                cmdexecute(cmd)
                        	if Tracker["importali3d"] !="":
                                	cmd = "{} {} {} {} ".format("sxheader.py", Tracker["this_data_stack"],"--params=xform.projection",\
                                 	"--import="+Tracker["importali3d"])
                        		cmdexecute(cmd)
                        mpi_barrier(MPI_COMM_WORLD)
			if myid==main_node:
				log_main.add("Create referece volumes for %3d independent runs"%Tracker["constants"]["indep_runs"])
			N_independent_reconstructions(Tracker)
			if myid ==main_node:
				log_main.add("----------Equal-Kmeans---------")
			prepare_EMREF_dict(Tracker)
			N_independent_mref(mpi_comm,Tracker)
			do_two_way_comparison(Tracker)
			mpi_barrier(MPI_COMM_WORLD)
			if myid ==main_node:
				log_main.add("--------EKmref------- ")
			do_EKmref(Tracker)
			if myid ==main_node:
                        	log_main.add("----------EKmref---------")
			number_of_groups=min(Tracker["number_of_groups"],int(Tracker["number_of_unaccounted"]/float(max(Tracker["constants"]["members_per_group"],Tracker["average_members_in_a_group"]))))
			Tracker["last_dir"]   =os.path.join(masterdir,"generation%03d"%generation)
			Tracker["last_data_stack"]="bdb:"+os.path.join(Tracker["last_dir"],"data")
			if myid ==main_node:
				log_main.add("Kgroup for the next generation is %d"%number_of_groups)
				log_main.add("the number of outliers is  %5d"%Tracker["number_of_unaccounted"]) 
				log_main.add("average_members_in_a_group  %5d"%int(Tracker["average_members_in_a_group"]))
			mpi_barrier(MPI_COMM_WORLD)
			while number_of_groups >=2 and Tracker["constants"]["do_unaccounted"]:
				generation            +=1
				if myid ==main_node:
                        		log_main.add("clustering generation %03d"%generation)
                		Tracker["this_dir"]   =os.path.join(masterdir,"generation%03d"%generation)
                		if myid ==main_node:
                        		cmd="{} {}".format("mkdir",Tracker["this_dir"])
                        		cmdexecute(cmd)
                		Tracker["this_data_stack"]  ="bdb:"+os.path.join(Tracker["this_dir"],"data")
                		Tracker["this_total_stack"] =Tracker["number_of_unaccounted"]
                		Tracker["number_of_groups"] =number_of_groups
				Tracker["this_ali3d"]       =Tracker["ali3d_of_outliers"]
				if myid ==main_node:
	                        	log_main.add("Kgroup for the next generation is %d"%number_of_groups)
                                	log_main.add("the number of outliers is  %5d"%Tracker["number_of_unaccounted"])
                                	log_main.add("average_members_in_a_group  %5d"%int(Tracker["average_members_in_a_group"]))
				### Create data stack
				if myid==main_node:
					log_main.add("----Create data stack unaccounted in the last generation")
                			cmd = "{} {} {} {} {}".format("e2bdb.py",Tracker["last_data_stack"],"--makevstack",\
                        		Tracker["this_data_stack"],"--list="+Tracker["this_unaccounted"])
                			cmdexecute(cmd)
					cmd = "{} {}".format("sxheader.py  --consecutive", Tracker["this_data_stack"])
                        		cmdexecute(cmd)
                			#cmd = "{} {} {} {}".format("sxheader.py", Tracker["this_data_stack"],"--params=xform.projection",\
                                 	#"--import="+Tracker["this_ali3d"])
					#cmdexecute(cmd)
				mpi_barrier(MPI_COMM_WORLD)
				if myid ==main_node:
                        		log_main.add("-----Create referece volumes for %3d independent runs"%Tracker["constants"]["indep_runs"])
                		N_independent_reconstructions(Tracker)
				if myid ==main_node:
                        		log_main.add("---------Equal-Kmeans----------")
                		prepare_EMREF_dict(Tracker)
                		N_independent_mref(mpi_comm,Tracker)
				do_two_way_comparison(Tracker)
                		mpi_barrier(MPI_COMM_WORLD)
                		if myid ==main_node:
                        		log_main.add("---------EKmref---------")
                		do_EKmref(Tracker)
                		if myid ==main_node:
                        		log_main.add("------- end of EKmref-------")
                		number_of_groups=min(Tracker["constants"]["Kgroup"],int(Tracker["number_of_unaccounted"]/float(max(Tracker["constants"]["members_per_group"],Tracker["average_members_in_a_group"]))))
				Tracker["last_dir"]   =os.path.join(masterdir,"generation%03d"%generation)
                		Tracker["last_data_stack"]="bdb:"+os.path.join(Tracker["last_dir"],"data")
		# Finish program
               	mpi_barrier(MPI_COMM_WORLD)
               	from mpi import mpi_finalize
               	mpi_finalize()
		exit()
if __name__ == "__main__":
	main()
