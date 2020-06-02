




































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































from __future__ import print_function
def volume_reconstruction(data, options, mpi_comm):
	pass#IMPORTIMPORTIMPORT from mpi import mpi_comm_rank
	pass#IMPORTIMPORTIMPORT from sp_reconstruction import recons3d_4nn_MPI, recons3d_4nn_ctf_MPI
	pass#IMPORTIMPORTIMPORT from sp_utilities import bcast_EMData_to_all, model_circle
	
	myid = mpi.mpi_comm_rank(mpi_comm)
	sym  = options.sym
	sym = sym[0].lower() + sym[1:]
	npad      = options.npad
	user_func = options.user_func
	CTF       = options.CTF
	snr       = options.snr
	center    = options.center
	#=========================================================================
	# volume reconstruction
	if CTF: vol = sp_reconstruction.recons3d_4nn_ctf_MPI(myid, data, snr, symmetry=sym, npad=npad, mpi_comm=mpi_comm)
	else:   vol = sp_reconstruction.recons3d_4nn_MPI    (myid, data,      symmetry=sym, snr=snr, npad=npad, mpi_comm=mpi_comm)

	if myid == 0:
		nx = data[0].get_xsize()
		last_ring   = int(options.ou)
		mask3D = sp_utilities.model_circle(last_ring, nx, nx, nx)
		ref_data = [ mask3D, max(center,0), None, None, None, None ]
		ref_data[2] = vol
		ref_data[3] = None #fscc
		ref_data[4] = None#varf
		#  call user-supplied function to prepare reference image, i.e., center and filter it
		vol, cs = user_func(ref_data)


	# broadcast volume
	sp_utilities.bcast_EMData_to_all(vol, myid, 0, comm=mpi_comm)
	#=========================================================================
	return vol


def volume_recsp(data, options):

	pass#IMPORTIMPORTIMPORT from sp_reconstruction import recons3d_4nn, recons3d_4nn_ctf
	pass#IMPORTIMPORTIMPORT from sp_utilities import bcast_EMData_to_all, model_circle
	
	sym  = options.sym
	sym = sym[0].lower() + sym[1:]
	npad      = options.npad
	user_func = options.user_func
	CTF       = options.CTF
	snr       = options.snr
	center    = options.center
	#=========================================================================
	# volume reconstruction
	if CTF: vol = sp_reconstruction.recons3d_4nn_ctf(data, snr, symmetry=sym, npad=npad)
	else:   vol = sp_reconstruction.recons3d_4nn(data,      symmetry=sym, npad=npad)

	nx = data[0].get_xsize()
	last_ring   = int(options.ou)
	mask3D = sp_utilities.model_circle(last_ring, nx, nx, nx)
	ref_data = [ mask3D, max(center,0), None, None, None, None ]
	ref_data[2] = vol
	ref_data[3] = None #fscc
	ref_data[4] = None#varf
	#  call user-supplied function to prepare reference image, i.e., center and filter it
	vol, cs = user_func(ref_data)

	#=========================================================================
	return vol





















































































































































































































































































































































def proj_ali_incore_multi(data, refrings, numr, xrng = 0.0, yrng = 0.0, step=1.0, an = 1.0, nsoft = -1, finfo=None, sym="c1"):
	pass#IMPORTIMPORTIMPORT from sp_utilities    import compose_transform2
	pass#IMPORTIMPORTIMPORT from math         import cos, pi, radians, degrees
	pass#IMPORTIMPORTIMPORT from EMAN2 import Vec2f, Transform
	pass#IMPORTIMPORTIMPORT from sp_global_def import Util
	pass#IMPORTIMPORTIMPORT from sp_global_def import ERROR

	mode = "F"
	nx   = data.get_xsize()
	ny   = data.get_ysize()
	#  center is in SPIDER convention
	cnx  = nx//2 + 1
	cny  = ny//2 + 1
	ant = math.cos(math.radians(an))

	#phi, theta, psi, sxo, syo = get_params_proj(data)
	t1 = data.get_attr("xform.projection")
	dp = t1.get_params("spider")
	if finfo:
		ID = data.get_attr("ID")
		finfo.write("Image id: %6d\n"%(ID))
		finfo.write("Old parameters: %9.4f %9.4f %9.4f %9.4f %9.4f\n"%(dp["phi"], dp["theta"], dp["psi"], -dp["tx"], -dp["ty"]))
		finfo.flush()
		
	ou = numr[-3]
	sxi = dp["tx"]
	syi = dp["ty"]
	txrng = [0.0]*2 
	tyrng = [0.0]*2
	sp_global_def.ERROR("proj_ali_incore_multi","Needs corrections",1)
	txrng[0] = max(0,min(cnx+sxi-ou, xrng+sxi))
	txrng[1] = max(0, min(nx-cnx-sxi-ou, xrng-sxi))
	tyrng[0] = max(0,min(cny+syi-ou, yrng+syi))
	tyrng[1] = max(0, min(ny-cny-syi-ou, yrng-syi))
		
	#print "Old parameters: %9.4f %9.4f %9.4f %9.4f %9.4f\n"%(dp["phi"], dp["theta"], dp["psi"], -dp["tx"], -dp["ty"])
	#[ang, sxs, sys, mirror, iref, peak, checked_refs] = Util.shc(data, refrings, xrng, yrng, step, ant, mode, numr, cnx+dp["tx"], cny+dp["ty"])
	peaks = EMAN2_cppwrap.Util.multiref_polar_ali_2d_peaklist_local(data, refrings, txrng, tyrng, step, ant, mode, numr, cnx+sxi, cny+syi)
	peaks_count = len(peaks) / 5
	#pixel_error = 0.0
	peak = 0.0
	if( peaks_count > 0 ):
		if( nsoft == -1 ):  nsoft = peaks_count
		params = [None]*peaks_count
		#                                              peak         iref      ang  sxs  sys 
		for i in range(peaks_count):  params[i] = [ peaks[i*5+0], int(peaks[i*5+4]), peaks[i*5+1], peaks[i*5+2], peaks[i*5+3]]
		params.sort(reverse=True)
		if(peaks_count < nsoft ):
			for i in range(peaks_count,nsoft,1): params.insert(0,params[0])
			peaks_count = nsoft
		elif( peaks_count > nsoft ):  peaks_count = min(peaks_count, nsoft)
		ws = sum([params[i][0] for i in range(peaks_count)])
		for i in range(peaks_count):
			iref   = params[i][1]
			ang    = params[i][2]
			sxs    = params[i][3]
			sys    = params[i][4]
			#mirror = 0
			peak   = params[i][0]/ws
			# The ormqip returns parameters such that the transformation is applied first, the mirror operation second.
			# What that means is that one has to change the Eulerian angles so they point into mirrored direction: phi+180, 180-theta, 180-psi
			angb, sxb, syb, ct = sp_utilities.compose_transform2(0.0, sxs, sys, 1, -ang, 0.0, 0.0, 1)
			"""Multiline Comment17"""
			#MULTILINEMULTILINEMULTILINE 17
				#MULTILINEMULTILINEMULTILINE 17
				#MULTILINEMULTILINEMULTILINE 17
				#MULTILINEMULTILINEMULTILINE 17
				#MULTILINEMULTILINEMULTILINE 17
				#MULTILINEMULTILINEMULTILINE 17
			#MULTILINEMULTILINEMULTILINE 17
			#MULTILINEMULTILINEMULTILINE 17
			phi   = refrings[iref].get_attr("phi")
			theta = refrings[iref].get_attr("theta")
			psi   = (refrings[iref].get_attr("psi")+angb+360.0)%360.0
			s2x   = sxb - sxi
			s2y   = syb - syi

			t2 = EMAN2_cppwrap.Transform({"type":"spider","phi":phi,"theta":theta,"psi":psi})
			t2.set_trans(EMAN2_cppwrap.Vec2f(-s2x, -s2y))
			if i == 0:
				data.set_attr("xform.projection", t2)
			else:
				data.set_attr("xform.projection" + str(i), t2)
			if i == 0:
				data.set_attr("weight", peak)
			else:
				data.set_attr("weight" + str(i), peak)
			#from pixel_error import max_3D_pixel_error
			#pixel_error += max_3D_pixel_error(t1, t2, numr[-3])
			if finfo:
				finfo.write( "New parameters: %9.4f %9.4f %9.4f %9.4f %9.4f %10.5f  %11.3e\n\n" %(phi, theta, psi, s2x, s2y, peak, pixel_error))
				finfo.flush()
			#print  "New parameters: %9.4f %9.4f %9.4f %9.4f %9.4f %10.5f  %11.3e\n\n" %(phi, theta, psi, s2x, s2y, peak, pixel_error)

		# remove old xform.projection
		i = max(peaks_count, 1)
		while data.has_attr("xform.projection" + str(i)):
			data.del_attr("xform.projection" + str(i))
			i += 1
		i = max(peaks_count, 1)
		while data.has_attr("weight" + str(i)):
			data.del_attr("weight" + str(i))
			i += 1
		#pixel_error /= peaks_count
	return ws
	"""Multiline Comment18"""
		#MULTILINEMULTILINEMULTILINE 18

	#MULTILINEMULTILINEMULTILINE 18
	#MULTILINEMULTILINEMULTILINE 18

def shc_multi(data, refrings, numr, xrng, yrng, step, an, nsoft, sym, finfo=None):
	pass#IMPORTIMPORTIMPORT from sp_utilities    import compose_transform2
	pass#IMPORTIMPORTIMPORT from sp_fundamentals import mirror
	pass#IMPORTIMPORTIMPORT from math         import cos, pi
	pass#IMPORTIMPORTIMPORT from EMAN2 import Vec2f, Transform
	pass#IMPORTIMPORTIMPORT from sp_global_def import Util

	ID = data.get_attr("ID")

	mode = "F"
	nx   = data.get_xsize()
	ny   = data.get_ysize()
	#  center is in SPIDER convention
	cnx  = nx//2 + 1
	cny  = ny//2 + 1

	ant = math.cos(an*math.pi/180.0)

	if finfo:
		t1 = data.get_attr("xform.projection")
		dp = t1.get_params("spider")
		finfo.write("Image id: %6d\n"%(ID))
		finfo.write("Old parameters: %9.4f %9.4f %9.4f %5.2f %5.2f   %10.3f\n"%(dp["phi"], dp["theta"], dp["psi"], -dp["tx"], -dp["ty"], \
					data.get_attr("previousmax")))
		
		i = 1
		while data.has_attr("xform.projection" + str(i)):
			t1 = data.get_attr("xform.projection" + str(i))
			dp = t1.get_params("spider")
			finfo.write("Add parameters: %9.4f %9.4f %9.4f %5.2f %5.2f   %10.3f\n"%(dp["phi"], dp["theta"], dp["psi"], -dp["tx"], -dp["ty"]))
			i += 1

	t1 = data.get_attr("xform.projection")

	#[ang, sxs, sys, mir, iref, peak, checked_refs] = Util.shc(data, refrings, xrng, yrng, step, ant, mode, numr, cnx+dp["tx"], cny+dp["ty"])
	#peaks = Util.shc_multipeaks(data, refrings, xrng, yrng, step, ant, mode, numr, cnx+dp["tx"], cny+dp["ty"], nsoft)
	#  Do not shift the image to prevent sliding away

	ou = numr[-3]
	txrng = [0.0]*2 
	tyrng = [0.0]*2
	txrng[0] = max(0,min(cnx-ou, xrng))
	txrng[1] = max(0, min(nx-cnx-ou, xrng))
	tyrng[0] = max(0,min(cny-ou, yrng))
	tyrng[1] = max(0, min(ny-cny-ou, yrng))
		
	peaks = EMAN2_cppwrap.Util.shc_multipeaks(data, refrings, txrng, tyrng, step, ant, mode, numr, cnx, cny, nsoft)
	peaks_count = len(peaks) / 7
	pixel_error = 0.0
	number_of_checked_refs = 0
	peak = 0.0
	if( peaks_count > 0 ):
		params = [None]*peaks_count
		#                                              peak         iref                  ang        sxs           sys           mir           checked references
		for i in range(peaks_count):  params[i] = [ peaks[i*7+5], int(peaks[i*7+4]), peaks[i*7+0], peaks[i*7+1], peaks[i*7+2], int(peaks[i*7+3]), int(peaks[i*7+6])]
		#  Make sure nothing is repeated
		if(peaks_count>1):
			taken = [params[k][1] for k in range(peaks_count)]
			pass#IMPORTIMPORTIMPORT from sp_utilities import findall
			i = 0
			while(i<peaks_count):
				ll = sp_utilities.findall(taken[i], taken)
				if(len(ll) > 1):
					sp_global_def.sxprint("  PROBLEM, found the same orientation more than once !  ")
					for k in range(len(params)):  sp_global_def.sxprint(params[k])
					ll.sort(reverse=True)
					for k in range(0,len(ll)-1):
						del params[k]
						peaks_count -= 1
					taken = [params[k][1] for k in range(peaks_count)]
				i+=1
		params.sort(reverse=True)
		ws = sum([params[i][0] for i in range(peaks_count)])  # peaks could be stretched
		for i in range(peaks_count):
			ang    = params[i][2]
			sxs    = params[i][3]
			sys    = params[i][4]
			mir    = params[i][5]
			iref   = params[i][1]
			#peak   = peaks[i*7+5]
			#checked_refs = int(peaks[i*7+6])
			#number_of_checked_refs += checked_refs
			#if(sxs>0.0 or sys >0.0):  print  "  SERROR in shc_multi  ",i,params[i]

			# The ormqip returns parameters such that the transformation is applied first, the mir operation second.
			# What that means is that one has to change the the Eulerian angles so they point into mired direction: phi+180, 180-theta, 180-psi
			angb, sxb, syb, ct = sp_utilities.compose_transform2(0.0, sxs, sys, 1, -ang, 0.0, 0.0, 1)
			if  mir:
				phi   = (refrings[iref].get_attr("phi")+540.0)%360.0
				theta = 180.0-refrings[iref].get_attr("theta")
				psi   = (540.0-refrings[iref].get_attr("psi")+angb)%360.0
				s2x   = sxb #- dp["tx"]
				s2y   = syb #- dp["ty"]
			else:
				phi   = refrings[iref].get_attr("phi")
				theta = refrings[iref].get_attr("theta")
				psi   = (refrings[iref].get_attr("psi")+angb+360.0)%360.0
				s2x   = sxb #- dp["tx"]
				s2y   = syb #- dp["ty"]
			#if(sxs>0.0 or sys >0.0):  print  "  SERROR2 in shc_multi  ",i,phi,theta,psi,s2x,s2y

			t2 = EMAN2_cppwrap.Transform({"type":"spider","phi":phi,"theta":theta,"psi":psi})
			t2.set_trans(EMAN2_cppwrap.Vec2f(-s2x, -s2y))
			#print i,phi,theta,psi
			if i == 0:
				data.set_attr("xform.projection", t2)
				data.set_attr("weight", params[i][0]/ws)
			else:
				data.set_attr("xform.projection" + str(i), t2)
				data.set_attr("weight" + str(i), params[i][0]/ws)
			pass#IMPORTIMPORTIMPORT from sp_pixel_error import max_3D_pixel_error
			pixel_error += sp_pixel_error.max_3D_pixel_error(t1, t2, numr[-3])
			#  preserve params, they might be needed if peaks_count<nsoft
			params[i] = [params[i][0], phi, theta, psi, s2x, s2y, iref]

		# Now set previousmax to a value halfway through
		data.set_attr("previousmax", params[peaks_count//2][0])

		#  Add orientations around the main peak with exclusion of those already taken
		#  Allow matching to bsoft>nsoft, but still keep nsoft.  This should allow elongated neighborhood
		bsoft = 2*nsoft
		if(peaks_count<nsoft):
			tempref = [refrings[i] for i in range(len(refrings))]
			taken   = [params[i][6] for i in range(peaks_count)]
			taken.sort(reverse=True)
			if(len(taken) > 1):
				for k in range(1,len(taken)):
					dod = []
					if( taken[k] == taken[k-1] ):
						sp_global_def.sxprint("  PROBLEM 2, entries duplicated  ",taken)
						dod.append(k)
				if(len(dod) >0):
					for k in dod:  del taken[k]
			#  delete taken
			try:
				for i in range(peaks_count):  del  tempref[taken[i]]
			except:
				sp_global_def.sxprint("  failed deleting tempref ")
				sp_global_def.sxprint(i,peaks_count,nsoft)
				sp_global_def.sxprint(" taken ",taken)
				sp_global_def.sxprint(len(tempref), len(refrings))
				pass#IMPORTIMPORTIMPORT from sys import exit
				exit()

			pass#IMPORTIMPORTIMPORT from sp_utilities import getfvec
			t1 = data.get_attr("xform.projection")
			dp = t1.get_params("spider")
			n1,n2,n3 = sp_utilities.getfvec(dp["phi"],dp["theta"])
			datanvec = [n1,n2,n3]
			if(int(sym[1:]) >1):
				iq = len(tempref)
				iq3 = 3*iq
				iq6 = 6*iq
				tempref += (tempref+tempref)
				refvecs = [None]*3*iq3
				dphi = 360.0/int(sym[1:])
				for i in range(iq):
					phi   = tempref[i].get_attr("phi")
					theta = tempref[i].get_attr("theta")
					n1,n2,n3 = sp_utilities.getfvec(phi-dphi,theta)
					refvecs[3*i+0] = n1
					refvecs[3*i+1] = n2
					refvecs[3*i+2] = n3
					n1,n2,n3 = sp_utilities.getfvec(phi+dphi,theta)
					refvecs[3*i+0+iq6] = n1
					refvecs[3*i+1+iq6] = n2
					refvecs[3*i+2+iq6] = n3
				for i in range(iq,2*iq):
					n1 = tempref[i].get_attr("n1")
					n2 = tempref[i].get_attr("n2")
					n3 = tempref[i].get_attr("n3")
					refvecs[3*i+0] = n1
					refvecs[3*i+1] = n2
					refvecs[3*i+2] = n3
			else: 
				refvecs = [None]*3*len(tempref)
				for i in range(len(tempref)):
					n1 = tempref[i].get_attr("n1")
					n2 = tempref[i].get_attr("n2")
					n3 = tempref[i].get_attr("n3")
					refvecs[3*i+0] = n1
					refvecs[3*i+1] = n2
					refvecs[3*i+2] = n3
			pass#IMPORTIMPORTIMPORT from sp_utilities import nearestk_to_refdir
			nrst = sp_utilities.nearestk_to_refdir(refvecs, datanvec, howmany = bsoft-peaks_count)
			del refvecs
			#  it does not use mir, do it by hand
			if( dp["theta"] > 90.0 ):  tdata = sp_fundamentals.mirror(data)
			else:                      tdata = data.copy()
			#  delete from tdata higher xform and weight and keep only base one as it will be used to do orientation search.
			#  In addition, zero shifts as here we always search around the origin to prevent sliding away.
			i = 1
			while tdata.has_attr("xform.projection" + str(i)):
				tdata.del_attr("xform.projection" + str(i))
				i += 1
			i = 1
			while tdata.has_attr("weight" + str(i)):
				tdata.del_attr("weight" + str(i))
				i += 1
			#  Search
			#if( dp["theta"] > 90.0 ): 
			#	print "  IS MIRRORED  ",dp["phi"], dp["theta"], dp["psi"], -dp["tx"], -dp["ty"]
			pws = proj_ali_incore_multi(tdata, [tempref[k] for k in nrst], numr, xrng, yrng, step, 180.0, bsoft-peaks_count, sym=sym)
			#  Can there be a problem with (0,0) direction??  PAP  05/25/2014
			for i in range(bsoft-peaks_count):
				if i == 0:    t1 = tdata.get_attr("xform.projection")
				else:         t1 = tdata.get_attr("xform.projection" + str(i))
				d = t1.get_params("spider")
				phi   = d["phi"]
				theta = d["theta"]
				psi   = d["psi"]
				s2x   = d["tx"]
				s2y   = d["ty"]
				if( dp["theta"] > 90.0 ):
					#  Change parameters if mirrored
					#print "  BEFORE MIRRORED  ",i,phi, theta, psi, s2x, s2y
					phi   = (phi+540.0)%360.0
					theta = 180.0-theta
					psi   = (540.0-psi)%360.0
					#print "  AFTER MIRRORED  ",i,phi, theta, psi, s2x, s2y
				if i == 0 :   w = tdata.get_attr("weight")
				else:         w = tdata.get_attr("weight" + str(i))
				w *= pws  # remove normalization
				params.append([w,  phi, theta, psi, s2x, s2y])
			#  From now on process nsfot largest
			params.sort(reverse=True)
			ws = sum([params[i][0] for i in range(nsoft)])  # peaks could be stretched
			for i in range(nsoft):
				#print  "  ADDITIONAL SOFT ASSIGNMENT  ",i,peaks_count,params[i][1],params[i][1],params[i][2],params[i][0]/ws
				t2 = EMAN2_cppwrap.Transform({"type":"spider","phi":params[i][1],"theta":params[i][2],"psi":params[i][3]})
				t2.set_trans(EMAN2_cppwrap.Vec2f(-params[i][4], -params[i][5]))
				#print i,phi,theta,psi
				if i == 0:
					data.set_attr("xform.projection", t2)
					data.set_attr("weight", params[i][0]/ws)
				else:
					data.set_attr("xform.projection" + str(i), t2)
					data.set_attr("weight" + str(i), params[i][0]/ws)

		if finfo:
			t1 = data.get_attr("xform.projection")
			dp = t1.get_params("spider")
			#finfo.write("Image id: %6d\n"%(ID))
			finfo.write("New parameters: %9.4f %9.4f %9.4f %5.2f %5.2f  %10.3f  %5.3f\n"%(dp["phi"], dp["theta"], dp["psi"], -dp["tx"], -dp["ty"], \
						data.get_attr("previousmax"), data.get_attr("weight")))
			i = 1
			while data.has_attr("xform.projection" + str(i)):
				t1 = data.get_attr("xform.projection" + str(i))
				dp = t1.get_params("spider")
				finfo.write("Add parameters: %9.4f %9.4f %9.4f %5.2f %5.2f         %5.3f\n"%(dp["phi"], dp["theta"], dp["psi"], -dp["tx"], -dp["ty"], \
						data.get_attr("weight" + str(i)) ) )
				finfo.flush()
				i += 1
		"""Multiline Comment19"""
		#MULTILINEMULTILINEMULTILINE 19
		#MULTILINEMULTILINEMULTILINE 19
		#MULTILINEMULTILINEMULTILINE 19
			#MULTILINEMULTILINEMULTILINE 19
			#MULTILINEMULTILINEMULTILINE 19
		#MULTILINEMULTILINEMULTILINE 19
		#MULTILINEMULTILINEMULTILINE 19
			#MULTILINEMULTILINEMULTILINE 19
			#MULTILINEMULTILINEMULTILINE 19
		#MULTILINEMULTILINEMULTILINE 19
		pixel_error /= peaks_count
		peak = params[0][0]  # It is not used anywhere, but set it to the maximum.
	
	#  if it did not find any higher peaks would do nothing and return peaks_count=0
	return peak, pixel_error, number_of_checked_refs, peaks_count


# parameters: list of (all) projections | reference volume | ...
def ali3d_multishc_soft(stack, ref_vol, ali3d_options, mpi_comm = None, log = None, nsoft=2 ):

	pass#IMPORTIMPORTIMPORT from sp_alignment       import Numrinit, prepare_refrings, proj_ali_incore_local
	pass#IMPORTIMPORTIMPORT from sp_utilities       import get_im, file_type, model_circle, get_input_from_string, get_params_proj, wrap_mpi_gatherv, wrap_mpi_bcast
	pass#IMPORTIMPORTIMPORT from mpi             import mpi_bcast, mpi_comm_size, mpi_comm_rank, MPI_FLOAT, MPI_COMM_WORLD, mpi_barrier, mpi_reduce, MPI_INT, MPI_SUM
	pass#IMPORTIMPORTIMPORT from sp_projection      import prep_vol
	pass#IMPORTIMPORTIMPORT from sp_statistics      import hist_list
	pass#IMPORTIMPORTIMPORT from sp_applications    import MPI_start_end
	pass#IMPORTIMPORTIMPORT from sp_filter          import filt_ctf
	pass#IMPORTIMPORTIMPORT from sp_global_def      import Util
	pass#IMPORTIMPORTIMPORT from EMAN2           import EMUtil, EMData
	pass#IMPORTIMPORTIMPORT import types
	pass#IMPORTIMPORTIMPORT from time            import time

	ir     = ali3d_options.ir
	rs     = ali3d_options.rs
	ou     = ali3d_options.ou
	xr     = ali3d_options.xr
	yr     = ali3d_options.yr
	ts     = ali3d_options.ts
	an     = ali3d_options.an
	sym    = ali3d_options.sym
	sym    = sym[0].lower() + sym[1:]
	delta  = ali3d_options.delta
	center = ali3d_options.center
	CTF    = ali3d_options.CTF
	ref_a  = ali3d_options.ref_a

	if mpi_comm == None:
		mpi_comm = mpi.MPI_COMM_WORLD

	if log == None:
		pass#IMPORTIMPORTIMPORT from sp_logger import Logger
		log = sp_logger.Logger()

	number_of_proc = mpi.mpi_comm_size(mpi_comm)
	myid           = mpi.mpi_comm_rank(mpi_comm)
	main_node = 0

	if myid == main_node:
		log.add("Start ali3d_multishc_soft")

	xrng        = sp_utilities.get_input_from_string(xr)
	if  yr == "-1":  yrng = xrng
	else          :  yrng = sp_utilities.get_input_from_string(yr)
	step        = sp_utilities.get_input_from_string(ts)
	delta       = sp_utilities.get_input_from_string(delta)
	lstp = min(len(xrng), len(yrng), len(step), len(delta))
	if an == "-1":
		an = [-1] * lstp
	else:
		an = sp_utilities.get_input_from_string(an)

	first_ring  = int(ir)
	rstep       = int(rs)
	last_ring   = int(ou)
	max_iter    = int(ali3d_options.maxit)
	center      = int(center)

	if( type(ref_vol) is bytes ):  vol = sp_utilities.get_im(ref_vol)
	else:	vol = ref_vol
	nx      = vol.get_xsize()
	if last_ring < 0:	last_ring = int(nx/2) - 2

	numr	= sp_alignment.Numrinit(first_ring, last_ring, rstep, "F")
	mask2D  = sp_utilities.model_circle(last_ring,nx,nx) - sp_utilities.model_circle(first_ring,nx,nx)

	if( type(stack) is bytes ):
		if myid == main_node:
			if sp_utilities.file_type(stack) == "bdb":
				pass#IMPORTIMPORTIMPORT from EMAN2db import db_open_dict
				dummy = EMAN2db.db_open_dict(stack, True)
			# horatio active_refactoring Jy51i1EwmLD4tWZ9_00000_1
			# active = EMUtil.get_all_attributes(stack, 'active')
			# list_of_particles = []
			# for im in xrange(len(active)):
			# 	if active[im]:  list_of_particles.append(im)
			# del active
			nima = EMAN2_cppwrap.EMUtil.get_image_count(stack)
			list_of_particles = list(range(nima))
			
			total_nima = len(list_of_particles)
		else:
			list_of_particles = None
			total_nima = 0

	else:
		if myid == main_node:
			list_of_particles = list(range(len(stack)))
			total_nima = len(list_of_particles)
		else:
			list_of_particles = None
			total_nima = None
	total_nima = sp_utilities.wrap_mpi_bcast(total_nima, main_node, mpi_comm)
	list_of_particles = sp_utilities.wrap_mpi_bcast(list_of_particles, main_node, mpi_comm)

	image_start, image_end = sp_applications.MPI_start_end(total_nima, number_of_proc, myid)
	# create a list of images for each node
	list_of_particles = list_of_particles[image_start: image_end]
	nima = len(list_of_particles)

	if( type(stack) is bytes ):  data = EMAN2_cppwrap.EMData.read_images(stack, list_of_particles)
	else:                                   data = [ stack[im] for im in list_of_particles ]
	for im in range(nima):
		data[im].set_attr('ID', list_of_particles[im])
		ctf_applied = data[im].get_attr_default('ctf_applied', 0)
		if CTF and ctf_applied == 0:
			ctf_params = data[im].get_attr("ctf")
			st = EMAN2_cppwrap.Util.infomask(data[im], mask2D, False)
			data[im] -= st[0]
			data[im] = sp_filter.filt_ctf(data[im], ctf_params)
			data[im].set_attr('ctf_applied', 1)

	pixer = [0.0]*nima
	#par_r = [[] for im in list_of_particles ]
	cs = [0.0]*3
	total_iter = 0
	# do the projection matching
	for N_step in range(lstp):

		terminate = 0
		Iter = 0
		while Iter < max_iter and terminate == 0:

			Iter += 1
			total_iter += 1

			mpi.mpi_barrier(mpi_comm)
			if myid == main_node:
				log.add("ITERATION #%3d,  inner iteration #%3d"%(total_iter, Iter))
				log.add("Delta = %4.1f, an = %5.2f, xrange = %5.2f, yrange = %5.2f, step = %5.2f\n"%(delta[N_step], an[N_step], xrng[N_step],yrng[N_step],step[N_step]))
				start_time = time.time()

			#=========================================================================
			# build references
			volft, kb = sp_projection.prep_vol(vol)
			refrings = sp_alignment.prepare_refrings(volft, kb, nx, delta[N_step], ref_a, sym, numr, MPI=mpi_comm)
			del volft, kb
			#=========================================================================

			if myid == main_node:
				log.add("Time to prepare rings: %f\n" % (time.time()-start_time))
				start_time = time.time()
			
			#=========================================================================
			if total_iter == 1:
				# adjust params to references, calculate psi+shifts, calculate previousmax
				for im in range(nima):
					previousmax = data[im].get_attr_default("previousmax", -1.0e23)
					if(previousmax == -1.0e23):
						peak, pixer[im] = sp_alignment.proj_ali_incore_local(data[im],refrings,numr,xrng[N_step],yrng[N_step],step[N_step],10.0, sym=sym)
						data[im].set_attr("previousmax", peak*0.9)
				if myid == main_node:
					log.add("Time to calculate first psi+shifts+previousmax: %f\n" % (time.time()-start_time))
					start_time = time.time()
			#=========================================================================

			mpi.mpi_barrier(mpi_comm)
			if myid == main_node:
				start_time = time.time()
			#=========================================================================
			# alignment
			#number_of_checked_refs = 0
			par_r = [0]*max(2,(nsoft+1))
			for im in range(nima):
				sp_global_def.ERROR("shc_multi","Needs corrections")
				peak, pixer[im], checked_refs, number_of_peaks = shc_multi(data[im], refrings, numr, xrng[N_step], yrng[N_step], step[N_step],\
																			an[N_step], nsoft, sym)
				#number_of_checked_refs += checked_refs
				par_r[number_of_peaks] += 1
				#print  myid,im,number_of_peaks
				#t = get_params_proj(data[im])
				#if(t[3] >0.0 or t[4]>0.0):  print  "  MERRROR  ",t
				
			#=========================================================================
			mpi.mpi_barrier(mpi_comm)
			if myid == main_node:
				#print  data[0].get_attr_dict()
				log.add("Time of alignment = %f\n"%(time.time()-start_time))
				start_time = time.time()

			#=========================================================================
			#output pixel errors, check stop criterion
			all_pixer = sp_utilities.wrap_mpi_gatherv(pixer, 0, mpi_comm)
			par_r = mpi.mpi_reduce(par_r, len(par_r), mpi.MPI_INT, mpi.MPI_SUM, 0, mpi.MPI_COMM_WORLD)
			#total_checked_refs = wrap_mpi_gatherv([number_of_checked_refs], main_node, mpi_comm)
			terminate = 0
			if myid == main_node:
				#total_checked_refs = sum(total_checked_refs)
				log.add("=========== Number of better peaks found ==============")
				for lhx in range(nsoft+1):
					msg = "            %5d     %7d"%(lhx, par_r[lhx])
					log.add(msg)
				log.add("_______________________________________________________")

				lhist = 20
				region, histo = sp_statistics.hist_list(all_pixer, lhist)
				log.add("=========== Histogram of pixel errors ==============")
				for lhx in range(lhist):
					msg = "          %10.3f     %7d"%(region[lhx], histo[lhx])
					log.add(msg)
				log.add("____________________________________________________")
				if (max(all_pixer) < 0.5) and (sum(all_pixer)/total_nima < 0.05):
					terminate = 1
					log.add("...............")
					log.add(">>>>>>>>>>>>>>>   Will terminate due to small pixel errors")
			terminate = sp_utilities.wrap_mpi_bcast(terminate, main_node, mpi_comm)
			#=========================================================================

			#=========================================================================
			# centering
			if center == -1 and sym[0] == 'c':
				pass#IMPORTIMPORTIMPORT from sp_utilities      import estimate_3D_center_MPI, rotate_3D_shift
				cs[0], cs[1], cs[2], dummy, dummy = sp_utilities.estimate_3D_center_MPI(data, total_nima, myid, number_of_proc, main_node, mpi_comm=mpi_comm)
				if myid == main_node:
					msg = " Average center x = %10.3f        Center y = %10.3f        Center z = %10.3f\n"%(cs[0], cs[1], cs[2])
					log.add(msg)
				if int(sym[1]) > 1:
					cs[0] = cs[1] = 0.0
					if myid == main_node:
						log.add("For symmetry group cn (n>1), we only center the volume in z-direction\n")
				cs = mpi.mpi_bcast(cs, 3, mpi.MPI_FLOAT, main_node, mpi_comm)
				cs = [-float(cs[0]), -float(cs[1]), -float(cs[2])]
				sp_utilities.rotate_3D_shift(data, cs)
			#=========================================================================

			#=========================================================================
			# volume reconstruction
			mpi.mpi_barrier(mpi_comm)
			if myid == main_node:
				start_time = time.time()
			vol = volume_reconstruction(data, ali3d_options, mpi_comm)
			if myid == main_node:  vol.write_image('soft/smvol%04d.hdf'%total_iter)
			# log
			if myid == main_node:
				log.add("3D reconstruction time = %f\n"%(time.time()-start_time))
				start_time = time.time()
			#=========================================================================

			#=========================================================================
			if(total_iter%1 == 0 or terminate):
				# gather parameters
				params = []
				previousmax = []
				for im in data:
					t = sp_utilities.get_params_proj(im)
					params.append( [t[0], t[1], t[2], t[3], t[4]] )
					#if(t[3] >0.0 or t[4]>0.0):  print  "  ERRROR  ",t
					previousmax.append(im.get_attr("previousmax"))
				assert(nima == len(params))
				params = sp_utilities.wrap_mpi_gatherv(params, 0, mpi_comm)
				if myid == 0:
					assert(total_nima == len(params))
				previousmax = sp_utilities.wrap_mpi_gatherv(previousmax, 0, mpi_comm)
				if myid == main_node:
					pass#IMPORTIMPORTIMPORT from sp_utilities import write_text_row, write_text_file
					sp_utilities.write_text_row(params, "soft/params%04d.txt"%total_iter)
					sp_utilities.write_text_file(previousmax, "soft/previousmax%04d.txt"%total_iter)
				del previousmax, params
				i = 1
				while data[0].has_attr("xform.projection" + str(i)):
					params = []
					previousmax = []
					for im in data:

						try:
							#print  im.get_attr("xform.projection" + str(i))
							t = sp_utilities.get_params_proj(im,"xform.projection" + str(i))
						except:
							sp_global_def.sxprint(" NO XFORM  ",myid, i,im.get_attr('ID'))
							pass#IMPORTIMPORTIMPORT from sys import exit
							exit()

						params.append( [t[0], t[1], t[2], t[3], t[4]] )
						#if(t[3] >0.0 or t[4]>0.0):  print  "  ERRROR  ",i,t
					assert(nima == len(params))
					params = sp_utilities.wrap_mpi_gatherv(params, 0, mpi_comm)
					if myid == 0:
						assert(total_nima == len(params))
					if myid == main_node:
						sp_utilities.write_text_row(params, "soft/params-%04d-%04d.txt"%(i,total_iter))
					del previousmax, params
					i+=1


	if myid == main_node:
		log.add("Finished ali3d_multishc_soft")
		return #params, vol, previousmax, par_r
	else:
		return #None, None, None, None  # results for the other processes


# data - projections (scattered between cpus) or the volume.  If volume, just do the volume processing
# options - the same for all cpus
# return - volume the same for all cpus




















































































def no_of_processors_restricted_by_data__do_volume(projections, ali3d_options, iter, mpi_comm):
	pass#IMPORTIMPORTIMPORT from mpi import mpi_comm_rank, mpi_comm_size, mpi_finalize, mpi_comm_split, mpi_barrier, MPI_COMM_WORLD
	pass#IMPORTIMPORTIMPORT from sp_utilities      import bcast_EMData_to_all
	pass#IMPORTIMPORTIMPORT from sp_applications import MPI_start_end

	mpi_size = mpi.mpi_comm_size(mpi_comm)
	n_projs = len(projections)
	mpi_rank = mpi.mpi_comm_rank(mpi_comm)

	if (mpi_size > n_projs):
		working = int(not(mpi_rank < n_projs))
		mpi_subcomm = mpi.mpi_comm_split(mpi_comm, working,  mpi_rank - working*n_projs)
		mpi_subsize = mpi.mpi_comm_size(mpi_subcomm)
		mpi_subrank = mpi.mpi_comm_rank(mpi_subcomm)
		if (mpi_rank < n_projs):
			proj_begin, proj_end = sp_applications.MPI_start_end(n_projs, mpi_subsize, mpi_subrank)
			ref_vol = do_volume(projections[proj_begin:proj_end], ali3d_options, 0, mpi_comm=mpi_subcomm)
		else:
			pass#IMPORTIMPORTIMPORT from sp_utilities import model_blank
			nx = projections[0].get_xsize()
			ref_vol = sp_utilities.model_blank(nx,nx,nx)
		sp_utilities.bcast_EMData_to_all(ref_vol, mpi_rank, 0, comm=mpi_comm)
	else:
		proj_begin, proj_end = sp_applications.MPI_start_end(n_projs, mpi_size, mpi_rank)
		ref_vol = do_volume(projections[proj_begin:proj_end], ali3d_options, 0, mpi_comm=mpi_comm)

	return ref_vol

























































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































