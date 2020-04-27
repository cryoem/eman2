






































from __future__ import print_function
def rec2D(  lines, idrange=None, snr=None ):
	""" Perform a 2D reconstruction on a set of 1D lines using nearest neighbouring reverse FFT algorithm.
		Input: a set of 1D lines
		Output: a 2D image
	"""
	pass#IMPORTIMPORTIMPORT from EMAN2 import Reconstructors

	assert len(lines) > 0


	size = lines[0].get_xsize();

	if snr is None:
		params = {"size":size, "npad":4, "ndim":2}
	else: 
		params = {"size":size, "npad":4, "ndim":2, "snr":snr}

	r = EMAN2_cppwrap.Reconstructors.get("nn4", params)
	r.setup()

	if idrange is None:
		idrange = list(range( len(lines)))

	t = EMAN2_cppwrap.Transform()
	for i in idrange:
		r.insert_slice( lines[i], t )

	return r.finish(True)


















def insert_slices_pdf(reconstructor, proj):
	xforms =   proj.get_attr("xform.projection") 
	weights =  proj.get_attr_default("weight", 1.0) 
	reconstructor.insert_slice( proj, xforms, weights )
	ixform = 0
	while True:
		ixform += 1
		xform_proj = proj.get_attr_default("xform.projection" + str(ixform), None)
		if xform_proj == None:
			return 
		weights = proj.get_attr_default("weight" + str(ixform), 1.0)
		reconstructor.insert_slice( proj, xforms, weights)











































































































































































































































































































































































































































def recons3d_4nnw_MPI(myid, prjlist, bckgdata, snr = 1.0, sign=1, symmetry="c1", finfo=None, npad=2, xysize=-1, zsize=-1, mpi_comm=None, smearstep = 0.0, fsc = None):
	"""
		recons3d_4nn_ctf - calculate CTF-corrected 3-D reconstruction from a set of projections using three Eulerian angles, two shifts, and CTF settings for each projeciton image
		Input
			stack: name of the stack file containing projection data, projections have to be squares
			prjlist: list of projections to be included in the reconstruction or image iterator
			bckgdata = [get_im("tsd.hdf"),read_text_file("data_stamp.txt")]
			snr: Signal-to-Noise Ratio of the data 
			sign: sign of the CTF 
			symmetry: point-group symmetry to be enforced, each projection will enter the reconstruction in all symmetry-related directions.
	"""
	pass#IMPORTIMPORTIMPORT from sp_utilities  import reduce_EMData_to_root, pad
	pass#IMPORTIMPORTIMPORT from EMAN2      import Reconstructors
	pass#IMPORTIMPORTIMPORT from sp_utilities  import iterImagesList, set_params_proj, model_blank
	pass#IMPORTIMPORTIMPORT from mpi        import MPI_COMM_WORLD
	pass#IMPORTIMPORTIMPORT import types

	if mpi_comm == None:
		mpi_comm = mpi.MPI_COMM_WORLD

	if type(prjlist) == list:
		prjlist = sp_utilities.iterImagesList(prjlist)
	if not prjlist.goToNext():
		sp_global_def.ERROR("empty input list","recons3d_4nnw_MPI",1)
	imgsize = prjlist.image().get_xsize()
	if prjlist.image().get_ysize() != imgsize:
		imgsize = max(imgsize, prjlist.image().get_ysize())
		dopad = True
	else:
		dopad = False
	prjlist.goToPrev()

	#  Do the FSC shtick.
	bnx     = imgsize*npad//2+1
	if  fsc:
		pass#IMPORTIMPORTIMPORT from math import sqrt
		pass#IMPORTIMPORTIMPORT from sp_utilities import reshape_1d
		t = [0.0]*len(fsc)
		for i in range(len(fsc)):
			t[i] = min(max(fsc[i],0.0), 0.999)
		t = sp_utilities.reshape_1d(t,len(t),npad*len(t))
		refvol = sp_utilities.model_blank(bnx,1,1,0.0)
		for i in range(len(fsc)):  refvol.set_value_at(i,t[i])
	else:
		refvol = sp_utilities.model_blank(bnx,1,1,1.0)
	refvol.set_attr("fudge", 1.0)

	fftvol = EMAN2_cppwrap.EMData()
	weight = EMAN2_cppwrap.EMData()

	if( smearstep > 0.0 ):
		#if myid == 0:  print "  Setting smear in prepare_recons_ctf"
		ns = 1
		smear = []
		for j in range(-ns,ns+1):
			if( j != 0):
				for i in range(-ns,ns+1):
					for k in range(-ns,ns+1):
						smear += [i*smearstep,j*smearstep,k*smearstep,1.0]
		# Deal with theta = 0.0 cases
		prj = []
		for i in range(-ns,ns+1):
			for k in range(-ns,ns+1):
				prj.append(i+k)
		for i in range(-2*ns,2*ns+1,1):
			smear += [i*smearstep,0.0,0.0,float(prj.count(i))]
		#if myid == 0:  print "  Smear  ",smear
		fftvol.set_attr("smear", smear)

	if (xysize == -1 and zsize == -1 ):
		params = {"size":imgsize, "npad":npad, "snr":snr, "sign":sign, "symmetry":symmetry, "refvol":refvol, "fftvol":fftvol, "weight":weight}
		r = EMAN2_cppwrap.Reconstructors.get( "nn4_ctfw", params )
	else:
		if ( xysize != -1 and zsize != -1):
			rx = float(xysize)/imgsize
			ry = float(xysize)/imgsize
			rz = float(zsize)/imgsize
		elif( xysize != -1):
			rx = float(xysize)/imgsize
			ry = float(xysize)/imgsize
			rz = 1.0
		else:
			rx = 1.0
			ry = 1.0
			rz = float(zsize)/imgsize
		#  There is an error here with sizeprojection  PAP 10/22/2014
		params = {"size":sizeprojection, "npad":npad, "snr":snr, "sign":sign, "symmetry":symmetry, "fftvol":fftvol, "weight":weight,"xratio":rx,"yratio":ry,"zratio":rz}
		r = EMAN2_cppwrap.Reconstructors.get( "nn4_ctf_rect", params )
	r.setup()

	
	#from utilities import model_blank, get_im, read_text_file
	#bckgdata = [get_im("tsd.hdf"),read_text_file("data_stamp.txt")]

	nnx = bckgdata[0].get_xsize()
	nny = bckgdata[0].get_ysize()
	bckgnoise = []
	for i in range(nny):
		prj = sp_utilities.model_blank(nnx)
		for k in range(nnx):  prj[k] = bckgdata[0].get_value_at(k,i)
		bckgnoise.append(prj)

	datastamp = bckgdata[1]
	if not (finfo is None): nimg = 0
	while prjlist.goToNext():
		prj = prjlist.image()
		try:
			stmp = nnx/0
			stmp = prj.get_attr("ptcl_source_image")
		except:
			try:
				stmp = prj.get_attr("ctf")
				stmp = round(stmp.defocus,4)
			except:
				sp_global_def.ERROR("Either ptcl_source_image or ctf has to be present in the header.","recons3d_4nnw_MPI",1, myid)
		try:
			indx = datastamp.index(stmp)
		except:
			sp_global_def.ERROR("Problem with indexing ptcl_source_image.","recons3d_4nnw_MPI",1, myid)

		if dopad:
			prj = sp_utilities.pad(prj, imgsize, imgsize, 1, "circumference")

		prj.set_attr("bckgnoise", bckgnoise[indx])
		insert_slices(r, prj)
		if not (finfo is None):
			nimg += 1
			finfo.write(" %4d inserted\n" %(nimg) )
			finfo.flush()
	del sp_utilities.pad
	if not (finfo is None): 
		finfo.write( "begin reduce\n" )
		finfo.flush()
		
	sp_utilities.reduce_EMData_to_root(fftvol, myid, comm=mpi_comm)
	sp_utilities.reduce_EMData_to_root(weight, myid, comm=mpi_comm)

	if not (finfo is None): 
		finfo.write( "after reduce\n" )
		finfo.flush()

	if myid == 0 :
		dummy = r.finish(True)
	else:
		pass#IMPORTIMPORTIMPORT from sp_utilities import model_blank
		if ( xysize == -1 and zsize == -1 ):
			fftvol = sp_utilities.model_blank(imgsize, imgsize, imgsize)
		else:
			if zsize == -1:
				fftvol = sp_utilities.model_blank(xysize, xysize, imgsize)
			elif xysize == -1:
				fftvol = sp_utilities.model_blank(imgsize, imgsize, zsize)
			else:
				fftvol = sp_utilities.model_blank(xysize, xysize, zsize)
	return fftvol






















































































































































































































































































































































































































def recons3d_4nnf_MPI(myid, list_of_prjlist, bckgdata, snr = 1.0, sign=1, symmetry="c1", finfo=None, npad=2, mpi_comm=None, smearstep = 0.0):
	"""
		recons3d_4nn_ctf - calculate CTF-corrected 3-D reconstruction from a set of projections using three Eulerian angles, two shifts, and CTF settings for each projeciton image
		Input
			list_of_prjlist: list of lists of projections to be included in the reconstruction
			bckgdata = [get_im("tsd.hdf"),read_text_file("data_stamp.txt")]
			snr: Signal-to-Noise Ratio of the data 
			sign: sign of the CTF
			symmetry: point-group symmetry to be enforced, each projection will enter the reconstruction in all symmetry-related directions.
	"""
	pass#IMPORTIMPORTIMPORT from sp_utilities  import reduce_EMData_to_root, random_string, get_im
	pass#IMPORTIMPORTIMPORT from EMAN2      import Reconstructors
	pass#IMPORTIMPORTIMPORT from sp_utilities  import model_blank
	pass#IMPORTIMPORTIMPORT from mpi        import MPI_COMM_WORLD, mpi_barrier
	pass#IMPORTIMPORTIMPORT import types
	pass#IMPORTIMPORTIMPORT from sp_statistics import fsc
	pass#IMPORTIMPORTIMPORT import datetime
	
	if mpi_comm == None:
		mpi_comm = mpi.MPI_COMM_WORLD

	imgsize = list_of_prjlist[0].get_xsize()
	
	if( smearstep > 0.0 ):
		#if myid == 0:  print "  Setting smear in prepare_recons_ctf"
		ns = 1
		smear = []
		for j in range(-ns,ns+1):
			if( j != 0):
				for i in range(-ns,ns+1):
					for k in range(-ns,ns+1):
						smear += [i*smearstep,j*smearstep,k*smearstep,1.0]
		# Deal with theta = 0.0 cases
		prj = []
		for i in range(-ns,ns+1):
			for k in range(-ns,ns+1):
				prj.append(i+k)
		for i in range(-2*ns,2*ns+1,1):
			smear += [i*smearstep,0.0,0.0,float(prj.count(i))]
		#if myid == 0:  print "  Smear  ",smear

	#from utilities import model_blank, get_im, read_text_file
	#bckgdata = [get_im("tsd.hdf"),read_text_file("data_stamp.txt")]

	nnx = bckgdata[0].get_xsize()
	nny = bckgdata[0].get_ysize()
	bckgnoise = []
	for i in range(nny):
		prj = sp_utilities.model_blank(nnx)
		for k in range(nnx):  prj[k] = bckgdata[0].get_value_at(k,i)
		bckgnoise.append(prj)

	datastamp = bckgdata[1]

	#  Do the FSC shtick.
	bnx     = imgsize*npad//2+1
	refvol = sp_utilities.model_blank(bnx)  # fill fsc with zeroes so the first reconstruction is done using simple Wiener filter.
	refvol.set_attr("fudge", 1.0)

	results_list = []
	fftvol_file =[]
	weight_file = []

	for iset in range(2):
		if not (finfo is None): nimg = 0

		fftvol = EMAN2_cppwrap.EMData()
		weight = EMAN2_cppwrap.EMData()
		if( smearstep > 0.0 ):  fftvol.set_attr("smear", smear)
	
		params = {"size":imgsize, "npad":npad, "snr":snr, "sign":sign, "symmetry":symmetry, "refvol":refvol, "fftvol":fftvol, "weight":weight}
		r = EMAN2_cppwrap.Reconstructors.get( "nn4_ctfw", params )
		r.setup()

		for image in list_of_prjlist[iset]:
			try:
				#raise ValueError('A very specific thing happened')
				stmp = image.get_attr("ptcl_source_image")
			except:
				try:
					stmp = image.get_attr("ctf")
					stmp = round(stmp.defocus,4)
				except:
					sp_global_def.ERROR("Either ptcl_source_image or ctf has to be present in the header.","recons3d_4nnw_MPI",1, myid)
			try:
				indx = datastamp.index(stmp)
			except:
				sp_global_def.ERROR("Problem with indexing ptcl_source_image.","recons3d_4nnf_MPI",1, myid)
	
			image.set_attr("bckgnoise", bckgnoise[indx])
			insert_slices(r, image)
			if not (finfo is None):
				nimg += 1
				finfo.write(" %4d inserted\n" %(nimg) )
				finfo.flush()

		if not (finfo is None): 
			finfo.write( "begin reduce\n" )
			finfo.flush()
	
		sp_utilities.reduce_EMData_to_root(fftvol, myid, main_node, comm=mpi_comm)
		sp_utilities.reduce_EMData_to_root(weight, myid, main_node, comm=mpi_comm)
		
		if not (finfo is None): 
			finfo.write( "after reduce\n" )
			finfo.flush()

		if myid == 0:
			"""Multiline Comment4"""
			#MULTILINEMULTILINEMULTILINE 4
			#MULTILINEMULTILINEMULTILINE 4
			#MULTILINEMULTILINEMULTILINE 4
			#MULTILINEMULTILINEMULTILINE 4
			#MULTILINEMULTILINEMULTILINE 4
			#MULTILINEMULTILINEMULTILINE 4
			dummy = r.finish(True)
			results_list.append(fftvol)
			#if(iset == 0):  fftvol.write_image(results_list[-1])

		mpi.mpi_barrier(mpi_comm)

	if myid == 0:
		fourier_shell_correlation = sp_statistics.fsc(results_list[0], results_list[1], 1.0)[1]
		"""Multiline Comment5"""
		#MULTILINEMULTILINEMULTILINE 5
		#MULTILINEMULTILINEMULTILINE 5
		#MULTILINEMULTILINEMULTILINE 5
		#MULTILINEMULTILINEMULTILINE 5
		#MULTILINEMULTILINEMULTILINE 5
			#MULTILINEMULTILINEMULTILINE 5

		#MULTILINEMULTILINEMULTILINE 5
		#MULTILINEMULTILINEMULTILINE 5
			#MULTILINEMULTILINEMULTILINE 5
			#MULTILINEMULTILINEMULTILINE 5
			#MULTILINEMULTILINEMULTILINE 5
			#MULTILINEMULTILINEMULTILINE 5
				#MULTILINEMULTILINEMULTILINE 5
			#MULTILINEMULTILINEMULTILINE 5

			#MULTILINEMULTILINEMULTILINE 5
			#MULTILINEMULTILINEMULTILINE 5
			#MULTILINEMULTILINEMULTILINE 5

			#MULTILINEMULTILINEMULTILINE 5
			#MULTILINEMULTILINEMULTILINE 5


		#MULTILINEMULTILINEMULTILINE 5
		#MULTILINEMULTILINEMULTILINE 5
		#MULTILINEMULTILINEMULTILINE 5
		#MULTILINEMULTILINEMULTILINE 5
	mpi.mpi_barrier(mpi_comm)
	if myid == 0:
		return results_list[0], results_list[1], fourier_shell_correlation
	else:
		return None, None, None

def recons3d_4nnfs_MPI(myid, main_node, prjlist, upweighted = True, finfo=None, mpi_comm=None, smearstep = 0.0, CTF = True, compensate = False, target_size=-1):
	"""
		recons3d_4nn_ctf - calculate CTF-corrected 3-D reconstruction from a set of projections using three Eulerian angles, two shifts, and CTF settings for each projeciton image
		Input
			list_of_prjlist: list of lists of projections to be included in the reconstruction
	"""
	pass#IMPORTIMPORTIMPORT from sp_utilities  import reduce_EMData_to_root, random_string, get_im
	pass#IMPORTIMPORTIMPORT from EMAN2      import Reconstructors
	pass#IMPORTIMPORTIMPORT from sp_utilities  import model_blank
	pass#IMPORTIMPORTIMPORT from sp_filter		import filt_table
	pass#IMPORTIMPORTIMPORT from mpi        import MPI_COMM_WORLD, mpi_barrier
	pass#IMPORTIMPORTIMPORT import types
	pass#IMPORTIMPORTIMPORT from sp_statistics import fsc
	pass#IMPORTIMPORTIMPORT import datetime
	pass#IMPORTIMPORTIMPORT from sp_reconstruction import insert_slices_pdf
	
	if mpi_comm == None:
		mpi_comm = mpi.MPI_COMM_WORLD
	imgsize = prjlist[0].get_ysize()  # It can be Fourier, so take y-size
	"""Multiline Comment6"""
	#MULTILINEMULTILINEMULTILINE 6
		#MULTILINEMULTILINEMULTILINE 6
		#MULTILINEMULTILINEMULTILINE 6
		#MULTILINEMULTILINEMULTILINE 6
		#MULTILINEMULTILINEMULTILINE 6
			#MULTILINEMULTILINEMULTILINE 6
				#MULTILINEMULTILINEMULTILINE 6
					#MULTILINEMULTILINEMULTILINE 6
						#MULTILINEMULTILINEMULTILINE 6
		#MULTILINEMULTILINEMULTILINE 6
		#MULTILINEMULTILINEMULTILINE 6
		#MULTILINEMULTILINEMULTILINE 6
			#MULTILINEMULTILINEMULTILINE 6
				#MULTILINEMULTILINEMULTILINE 6
		#MULTILINEMULTILINEMULTILINE 6
			#MULTILINEMULTILINEMULTILINE 6
		#MULTILINEMULTILINEMULTILINE 6
	#MULTILINEMULTILINEMULTILINE 6
	#from utilities import model_blank, get_im, read_text_file
	#bckgdata = [get_im("tsd.hdf"),read_text_file("data_stamp.txt")]

	
	"""Multiline Comment7"""
	#MULTILINEMULTILINEMULTILINE 7
	#MULTILINEMULTILINEMULTILINE 7
	#MULTILINEMULTILINEMULTILINE 7
	"""Multiline Comment8"""
	#MULTILINEMULTILINEMULTILINE 8
	#MULTILINEMULTILINEMULTILINE 8
		#MULTILINEMULTILINEMULTILINE 8
		#MULTILINEMULTILINEMULTILINE 8
		#MULTILINEMULTILINEMULTILINE 8
	#MULTILINEMULTILINEMULTILINE 8
	#datastamp = bckgdata[1]
	"""Multiline Comment9"""
	#MULTILINEMULTILINEMULTILINE 9
	#MULTILINEMULTILINEMULTILINE 9
		#MULTILINEMULTILINEMULTILINE 9
		#MULTILINEMULTILINEMULTILINE 9
		#MULTILINEMULTILINEMULTILINE 9
			#MULTILINEMULTILINEMULTILINE 9
			#MULTILINEMULTILINEMULTILINE 9
			#MULTILINEMULTILINEMULTILINE 9
			#MULTILINEMULTILINEMULTILINE 9
		#MULTILINEMULTILINEMULTILINE 9
	#MULTILINEMULTILINEMULTILINE 9
		#MULTILINEMULTILINEMULTILINE 9
		#MULTILINEMULTILINEMULTILINE 9
	#MULTILINEMULTILINEMULTILINE 9

	refvol = sp_utilities.model_blank(target_size)
	refvol.set_attr("fudge", 1.0)


	if CTF: do_ctf = 1
	else:   do_ctf = 0
	if not (finfo is None): nimg = 0

	fftvol = EMAN2_cppwrap.EMData()
	weight = EMAN2_cppwrap.EMData()
	#if( smearstep > 0.0 ):  fftvol.set_attr("smear", smear)


	pass#IMPORTIMPORTIMPORT from sp_utilities import info
	params = {"size":target_size, "npad":2, "snr":1.0, "sign":1, "symmetry":"c1", "refvol":refvol, "fftvol":fftvol, "weight":weight, "do_ctf": do_ctf}
	r = EMAN2_cppwrap.Reconstructors.get( "nn4_ctfw", params )
	r.setup()
	for image in prjlist:
		if not upweighted: insert_slices_pdf(r, sp_filter.filt_table(image, image.get_attr("bckgnoise")) )
		else:              insert_slices_pdf(r, image)

	if not (finfo is None): 
		finfo.write( "begin reduce\n" )
		finfo.flush()

	sp_utilities.reduce_EMData_to_root(fftvol, myid, main_node, comm=mpi_comm)
	sp_utilities.reduce_EMData_to_root(weight, myid, main_node, comm=mpi_comm)

	if not (finfo is None): 
		finfo.write( "after reduce\n" )
		finfo.flush()


	if myid == main_node:
		dummy = r.finish(compensate)
	mpi.mpi_barrier(mpi_comm)

	if myid == main_node: return fftvol, weight, refvol
	else: return None, None, None

def recons3d_4nnstruct_MPI(myid, main_node, prjlist, paramstructure, refang, delta, upweighted = True, mpi_comm=None, CTF = True, target_size=-1, avgnorm = 1.0, norm_per_particle = None):
	"""
		recons3d_4nn_ctf - calculate CTF-corrected 3-D reconstruction from a set of projections using three Eulerian angles, two shifts, and CTF settings for each projeciton image
		Input
			list_of_prjlist: list of lists of projections to be included in the reconstruction
	"""
	pass#IMPORTIMPORTIMPORT from sp_utilities  import reduce_EMData_to_root, random_string, get_im, findall
	pass#IMPORTIMPORTIMPORT from EMAN2      import Reconstructors
	pass#IMPORTIMPORTIMPORT from sp_utilities  import model_blank
	pass#IMPORTIMPORTIMPORT from sp_filter		import filt_table
	pass#IMPORTIMPORTIMPORT from mpi        import MPI_COMM_WORLD, mpi_barrier
	pass#IMPORTIMPORTIMPORT import types
	pass#IMPORTIMPORTIMPORT import datetime
	
	if mpi_comm == None: mpi_comm = mpi.MPI_COMM_WORLD

	refvol = sp_utilities.model_blank(target_size)
	refvol.set_attr("fudge", 1.0)


	if CTF: do_ctf = 1
	else:   do_ctf = 0

	fftvol = EMAN2_cppwrap.EMData()
	weight = EMAN2_cppwrap.EMData()

	pass#IMPORTIMPORTIMPORT from sp_utilities import info
	params = {"size":target_size, "npad":2, "snr":1.0, "sign":1, "symmetry":"c1", "refvol":refvol, "fftvol":fftvol, "weight":weight, "do_ctf": do_ctf}
	r = EMAN2_cppwrap.Reconstructors.get( "nn4_ctfw", params )
	r.setup()
	
	if norm_per_particle == None: norm_per_particle = len(prjlist)*[1.0]

	for im in range(len(prjlist)):
		#  parse projection structure, generate three lists:
		#  [ipsi+iang], [ishift], [probability]
		#  Number of orientations for a given image
		numbor = len(paramstructure[im][2])
		ipsiandiang = [ paramstructure[im][2][i][0]/1000  for i in range(numbor) ]
		allshifts   = [ paramstructure[im][2][i][0]%1000  for i in range(numbor) ]
		probs       = [ paramstructure[im][2][i][1] for i in range(numbor) ]
		#  Find unique projection directions
		tdir = list(set(ipsiandiang))
		bckgn = prjlist[im][0].get_attr("bckgnoise")
		#  For each unique projection direction:
		for ii in range(len(tdir)):
			#  Find the number of times given projection direction appears on the list, it is the number of different shifts associated with it.
			lshifts = sp_utilities.findall(tdir[ii], ipsiandiang)
			toprab  = 0.0
			for ki in range(len(lshifts)):  toprab += probs[lshifts[ki]]
			recdata = EMAN2_cppwrap.Util.mult_scalar(prjlist[im][allshifts[lshifts[0]]], probs[lshifts[0]]/toprab)
			recdata.set_attr_dict({"padffted":1, "is_complex":0})
			for ki in range(1,len(lshifts)):
				EMAN2_cppwrap.Util.add_img(recdata, EMAN2_cppwrap.Util.mult_scalar(prjlist[im][allshifts[lshifts[ki]]], probs[lshifts[ki]]/toprab))
			recdata.set_attr_dict({"padffted":1, "is_complex":1})
			if not upweighted:  recdata = sp_filter.filt_table(recdata, bckgn )
			recdata.set_attr("bckgnoise", bckgn )
			ipsi = tdir[ii]%100000
			iang = tdir[ii]/100000
			r.insert_slice( recdata, EMAN2_cppwrap.Transform({"type":"spider","phi":refang[iang][0],"theta":refang[iang][1],"psi":refang[iang][2]+ipsi*delta}), toprab*avgnorm/norm_per_particle[im])
	#  clean stuff
	del bckgn, recdata, tdir, ipsiandiang, allshifts, probs


	sp_utilities.reduce_EMData_to_root(fftvol, myid, main_node, comm=mpi_comm)
	sp_utilities.reduce_EMData_to_root(weight, myid, main_node, comm=mpi_comm)

	if myid == main_node:
		dummy = r.finish(True)
	mpi.mpi_barrier(mpi_comm)

	if myid == main_node: return fftvol, weight, refvol
	else: return None, None, None






















































































def recons3d_4nnstruct_MPI_test(myid, main_node, prjlist, paramstructure, refang, parameters, delta, upweighted = True, mpi_comm=None, CTF = True, target_size=-1, avgnorm = 1.0, norm_per_particle = None):
	"""
		recons3d_4nn_ctf - calculate CTF-corrected 3-D reconstruction from a set of projections using three Eulerian angles, two shifts, and CTF settings for each projeciton image
		Input
			list_of_prjlist: list of lists of projections to be included in the reconstruction
	"""
	pass#IMPORTIMPORTIMPORT from sp_utilities  import reduce_EMData_to_root, random_string, get_im, findall
	pass#IMPORTIMPORTIMPORT from EMAN2      import Reconstructors
	pass#IMPORTIMPORTIMPORT from sp_utilities  import model_blank
	pass#IMPORTIMPORTIMPORT from sp_filter		import filt_table
	pass#IMPORTIMPORTIMPORT from mpi        import MPI_COMM_WORLD, mpi_barrier
	pass#IMPORTIMPORTIMPORT import types
	pass#IMPORTIMPORTIMPORT from sp_statistics import fsc
	pass#IMPORTIMPORTIMPORT import datetime
	pass#IMPORTIMPORTIMPORT from sp_reconstruction import insert_slices_pdf
	
	if mpi_comm == None: mpi_comm = mpi.MPI_COMM_WORLD

	imgsize = prjlist[0][0].get_ysize()  # It can be Fourier, so take y-size

	refvol = sp_utilities.model_blank(target_size)
	refvol.set_attr("fudge", 1.0)


	if CTF: do_ctf = 1
	else:   do_ctf = 0

	fftvol = EMAN2_cppwrap.EMData()
	weight = EMAN2_cppwrap.EMData()

	pass#IMPORTIMPORTIMPORT from sp_utilities import info
	params = {"size":target_size, "npad":2, "snr":1.0, "sign":1, "symmetry":"c1", "refvol":refvol, "fftvol":fftvol, "weight":weight, "do_ctf": do_ctf}
	r = EMAN2_cppwrap.Reconstructors.get( "nn4_ctfw", params )
	r.setup()
	
	if norm_per_particle == None: norm_per_particle = len(prjlist)*[1.0]

	for im in range(len(prjlist)):
		#  parse projection structure, generate three lists:
		#  [ipsi+iang], [ishift], [probability]
		#  Number of orientations for a given image
		bckgn = prjlist[im][0].get_attr("bckgnoise")
		recdata = EMAN2_cppwrap.Util.mult_scalar(prjlist[im][0], parameters[im][5])
		recdata.set_attr_dict({"padffted":1, "is_complex":1})
		pass#IMPORTIMPORTIMPORT from sp_fundamentals import fshift
		recdata = sp_fundamentals.fshift(recdata, parameters[im][3], parameters[im][4])
		if not upweighted:  recdata = sp_filter.filt_table(recdata, bckgn )
		r.insert_slice( recdata, EMAN2_cppwrap.Transform({"type":"spider","phi":parameters[im][0],"theta":parameters[im][1],"psi":parameters[im][2]}), 1.0)
		"""Multiline Comment10"""
		#MULTILINEMULTILINEMULTILINE 10
		#MULTILINEMULTILINEMULTILINE 10
		#MULTILINEMULTILINEMULTILINE 10
		#MULTILINEMULTILINEMULTILINE 10
		#MULTILINEMULTILINEMULTILINE 10
		#MULTILINEMULTILINEMULTILINE 10
		#MULTILINEMULTILINEMULTILINE 10
		#MULTILINEMULTILINEMULTILINE 10
		#MULTILINEMULTILINEMULTILINE 10
			#MULTILINEMULTILINEMULTILINE 10
			#MULTILINEMULTILINEMULTILINE 10
			#MULTILINEMULTILINEMULTILINE 10
			#MULTILINEMULTILINEMULTILINE 10
			#MULTILINEMULTILINEMULTILINE 10
			#MULTILINEMULTILINEMULTILINE 10
			#MULTILINEMULTILINEMULTILINE 10
				#MULTILINEMULTILINEMULTILINE 10
			#MULTILINEMULTILINEMULTILINE 10
			#MULTILINEMULTILINEMULTILINE 10
			#MULTILINEMULTILINEMULTILINE 10
			#MULTILINEMULTILINEMULTILINE 10
			#MULTILINEMULTILINEMULTILINE 10
			#MULTILINEMULTILINEMULTILINE 10
		#MULTILINEMULTILINEMULTILINE 10
	#  clean stuff
	del bckgn, recdata


	sp_utilities.reduce_EMData_to_root(fftvol, myid, main_node, comm=mpi_comm)
	sp_utilities.reduce_EMData_to_root(weight, myid, main_node, comm=mpi_comm)

	if myid == main_node:
		dummy = r.finish(True)
	mpi.mpi_barrier(mpi_comm)

	if myid == main_node: return fftvol, weight, refvol
	else: return None, None, None







































































































































































































def recons3d_nn_SSNR(stack_name,  mask2D = None, ring_width=1, npad =1, sign=1, symmetry="c1", CTF = False, random_angles = 0):

	"""
	Perform a 3-D reconstruction using nearest neighbor interpolation and 
	calculate 3D spectral signal-to-noise ratio (SSNR)	   
	Input : stack_name - Name of the file with projection data.
		CTF        - 
	        symmetry   - Point group of the target molecule (default "c1")
		npad       - Times of padding applied, default is 1
		sign       - Currently not used, may be used in the future
		w          - The thickness of the shell, default is 1
		filename   - The filename in which you can save the SSNR results 
	Return: reconstructed 3D SSNR volume
        Usage : vol = recons3d_nn_SSNR(stack_name, CTF, symmetry, npad, snr, sign, w, filename])
	CTF true:
	variance at one voxel  = Gamma^2d->3d [ |F_k^2D|^2   +  ctf^2*|P^2D->3D(F^3D)|^2 -
	          -2*Real(conj(F_k^2D)*ctf*P^2D->3D(F^3D))]	
	signal  at one voxel   = Gamma^2d->3d [ |F_k^2D|^2  ]
	SSNR =  sum_rot [ wght*signal/Kn ]/sum_rot[ wght*variance /(Kn(Kn-1))] -1
	Notice: wght is always turned on during SSNR calculation.
	"""
	pass#IMPORTIMPORTIMPORT import types
	pass#IMPORTIMPORTIMPORT from EMAN2 import Reconstructors

	# Yang add a safety on 05/22/07
	if type(stack_name) == bytes: nima = EMAN2_cppwrap.EMUtil.get_image_count(stack_name)
	else :                                   nima = len(stack_name)
	# read first image to determine the size to use
	if type(stack_name) == bytes:
		proj = EMAN2_cppwrap.EMData()
		proj.read_image(stack_name, 0)
	else:    
		proj = stack_name[0].copy()
	#active = proj.get_attr('active')
	size   = proj.get_xsize()
	# sanity check -- image must be square
	if size != proj.get_ysize(): sp_global_def.ERROR("input data has to be square","recons3d_nn_SSNR",1)
	# reconstructor
	SSNR = EMAN2_cppwrap.EMData()
	fftvol = EMAN2_cppwrap.EMData()
	weight = EMAN2_cppwrap.EMData()
	weight2 = EMAN2_cppwrap.EMData()
	vol_ssnr = EMAN2_cppwrap.EMData()
	params = {"size":size, "npad":npad, "symmetry":symmetry, "SSNR":SSNR, "w":ring_width, "fftvol":fftvol, "weight":weight, "weight2":weight2, "vol_ssnr":vol_ssnr}
	if CTF:
		weight3 = EMAN2_cppwrap.EMData()
		params["sign"] = sign
		params["weight3"] = weight3
		r = EMAN2_cppwrap.Reconstructors.get("nnSSNR_ctf", params)
	else:
		r = EMAN2_cppwrap.Reconstructors.get("nnSSNR", params)
	r.setup()

	for i in range(nima):
		if type(stack_name) == bytes:
			proj.read_image(stack_name, i)
		else:
			proj = stack_name[i]
		# horatio active_refactoring Jy51i1EwmLD4tWZ9_00000_1
		# active = proj.get_attr_default('active', 1)
		# if(active == 1):
		if(random_angles  == 2):
			pass#IMPORTIMPORTIMPORT from  random import  random
			phi    = 360.0*random.random()
			theta  = 180.0*random.random()
			psi    = 360.0*random.random()
			xform_proj = EMAN2_cppwrap.Transform( {"type":"spider", "phi":phi, "theta":theta, "psi":psi} )
		elif(random_angles  == 3):
			pass#IMPORTIMPORTIMPORT from  random import  random
			phi    = 360.0*random.random()
			theta  = 180.0*random.random()
			psi    = 360.0*random.random()
			tx     = 6.0*(random.random() - 0.5)
			ty     = 6.0*(random.random() - 0.5)
			xform_proj = EMAN2_cppwrap.Transform( {"type":"spider", "phi":phi, "theta":theta, "psi":psi, "tx":tx, "ty":ty} )
		elif(random_angles  == 1):
			pass#IMPORTIMPORTIMPORT from  random import  random
			old_xform_proj = proj.get_attr( "xform.projection" )
			dict = old_xform_proj.get_rotation( "spider" )
			dict["psi"] = 360.0*random.random()
			xform_proj = EMAN2_cppwrap.Transform( dict )
		else:
			xform_proj = proj.get_attr( "xform.projection" )

		if mask2D:
			stats = EMAN2_cppwrap.Util.infomask(proj, mask2D, True)
			proj -= stats[0]
			proj *= mask2D
		r.insert_slice(proj, xform_proj)
		# horatio active_refactoring Jy51i1EwmLD4tWZ9_00000_1  END

	dummy = r.finish(True)
	outlist = [[] for i in range(6)]
	nn = SSNR.get_xsize()
	for i in range(1,nn): outlist[0].append((float(i)-0.5)/(float(nn-1)*2))
	for i in range(1,nn):
		if(SSNR(i,1,0) > 0.0):
			outlist[1].append(max(0.0,(SSNR(i,0,0)/SSNR(i,1,0)-1.)))      # SSNR
		else:
			outlist[1].append(0.0)
	for i in range(1,nn): outlist[2].append(SSNR(i,1,0)/SSNR(i,2,0))	          # variance
	for i in range(1,nn): outlist[3].append(SSNR(i,2,0))				  # number of points in the shell
	for i in range(1,nn): outlist[4].append(SSNR(i,3,0))				  # number of added Fourier points
	for i in range(1,nn): outlist[5].append(SSNR(i,0,0))				  # square of signal
	return [outlist, vol_ssnr]

def recons3d_nn_SSNR_MPI(myid, prjlist, mask2D, ring_width=1, npad =1, sign=1, symmetry="c1", CTF = False, random_angles = 0, mpi_comm = None):
	pass#IMPORTIMPORTIMPORT from sp_utilities import reduce_EMData_to_root
	pass#IMPORTIMPORTIMPORT from EMAN2 import Reconstructors
	pass#IMPORTIMPORTIMPORT from mpi import MPI_COMM_WORLD

	if mpi_comm == None:
		mpi_comm = mpi.MPI_COMM_WORLD

	if( len(prjlist) == 0 ):    sp_global_def.ERROR("empty input list","recons3d_nn_SSNR_MPI",1)
	imgsize = prjlist[0].get_xsize()
	if prjlist[0].get_ysize() != imgsize:  sp_global_def.ERROR("input data has to be square","recons3d_nn_SSNR_MPI",1)
	fftvol   = EMAN2_cppwrap.EMData()
	weight   = EMAN2_cppwrap.EMData()
	weight2  = EMAN2_cppwrap.EMData()
	SSNR     = EMAN2_cppwrap.EMData()
	vol_ssnr = EMAN2_cppwrap.EMData()
	params = {"size":imgsize, "npad":npad, "symmetry":symmetry, "SSNR":SSNR, "fftvol":fftvol, "weight":weight, "weight2":weight2, "vol_ssnr":vol_ssnr, "w":ring_width }
	if CTF:
		weight3  = EMAN2_cppwrap.EMData()
		params["sign"] = sign
		params["weight3"] = weight3
		r = EMAN2_cppwrap.Reconstructors.get("nnSSNR_ctf", params)
	else:
		r = EMAN2_cppwrap.Reconstructors.get("nnSSNR", params)
	r.setup()

	if prjlist[0].get_xsize() != imgsize or prjlist[0].get_ysize() != imgsize: sp_global_def.ERROR("inconsistent image size","recons3d_nn_SSNR_MPI",1)
	for prj in prjlist:
		# horatio active_refactoring Jy51i1EwmLD4tWZ9_00000_1
		# active = prj.get_attr_default('active', 1)
		# if active == 1:
		if random_angles  == 2:
			pass#IMPORTIMPORTIMPORT from  random import  random
			phi	 = 360.0*random.random()
			theta    = 180.0*random.random()
			psi	 = 360.0*random.random()
			xform_proj = EMAN2_cppwrap.Transform( {"type":"spider", "phi":phi, "theta":theta, "psi":psi} )
		elif random_angles  == 3:
			pass#IMPORTIMPORTIMPORT from  random import  random
			phi    = 360.0*random.random()
			theta  = 180.0*random.random()
			psi    = 360.0*random.random()
			tx     = 6.0*(random.random() - 0.5)
			ty     = 6.0*(random.random() - 0.5)
			xform_proj = EMAN2_cppwrap.Transform( {"type":"spider", "phi":phi, "theta":theta, "psi":psi, "tx":tx, "ty":ty} )
		elif random_angles  == 1:
			pass#IMPORTIMPORTIMPORT from  random import  random
			old_xform_proj = prj.get_attr( "xform.projection" )
			dict = old_xform_proj.get_rotation( "spider" )
			dict["psi"] = 360.0*random.random()
			xform_proj = EMAN2_cppwrap.Transform( dict )
		else:
			xform_proj = prj.get_attr( "xform.projection" )
		if mask2D:
			stats = EMAN2_cppwrap.Util.infomask(prj, mask2D, True)
			prj -= stats[0]
			prj *= mask2D
		r.insert_slice(prj, xform_proj )
		# horatio active_refactoring Jy51i1EwmLD4tWZ9_00000_1 END

	#from utilities import info
	sp_utilities.reduce_EMData_to_root(weight,  myid, 0, comm=mpi_comm)
	sp_utilities.reduce_EMData_to_root(fftvol,  myid, 0, comm=mpi_comm)
	sp_utilities.reduce_EMData_to_root(weight2, myid, 0, comm=mpi_comm)
	if CTF:
		sp_utilities.reduce_EMData_to_root(weight3, myid, 0, comm=mpi_comm)
	if myid == 0 :
		dummy = r.finish(True)		
		outlist = [[] for i in range(6)]
		nn = SSNR.get_xsize()
		for i in range(1,nn): outlist[0].append((float(i)-0.5)/(float(nn-1)*2))
		for i in range(1,nn):
			if SSNR(i,1,0) > 0.0:
				outlist[1].append(max(0.0,(SSNR(i,0,0)/SSNR(i,1,0)-1.)))     # SSNR
			else:
				outlist[1].append(0.0)
		for i in range(1,nn): 
			if SSNR(i,2,0) > 0.0:
				outlist[2].append(SSNR(i,1,0)/SSNR(i,2,0))	          # variance
			else:
				outlist[2].append(0.0)
		for i in range(1,nn): outlist[3].append(SSNR(i,2,0))				  # number of points in the shell
		for i in range(1,nn): outlist[4].append(SSNR(i,3,0))				  # number of added Fourier points
		for i in range(1,nn): outlist[5].append(SSNR(i,0,0))				  # square of signal
		return [outlist, vol_ssnr]


class memory_store(object):
	def __init__(self, npad):
		self.m_npad = npad
		self.m_imgs = []

	def add_image(self, img):
		self.m_imgs.append(img)

	def get_image(self, id):
		return self.m_imgs[id]

def bootstrap_nn(proj_stack, volume_stack, list_proj, niter, media="memory", npad=4, symmetry="c1", output=-1, CTF=False, snr=1.0, sign=1, myseed=None ):
	pass#IMPORTIMPORTIMPORT from random import seed
	pass#IMPORTIMPORTIMPORT from random import randint
	pass#IMPORTIMPORTIMPORT from time   import time
	pass#IMPORTIMPORTIMPORT from sys    import stdout
	pass#IMPORTIMPORTIMPORT from sp_utilities import set_ctf
	pass#IMPORTIMPORTIMPORT from EMAN2 import Reconstructors

	if(output == -1):
		pass#IMPORTIMPORTIMPORT import sys
		output=os.sys.stdout

	if not(myseed is None):
		random.seed(myseed) 

	nimages = len(list_proj)
	if nimages == 0 :
		sp_global_def.sxprint("empty list of projections input!")
		return None


	if media=="memory" :
		store = memory_store(npad)
	else :
		store = EMAN2_cppwrap.file_store(media,npad, 0, CTF)
		if not(output is None):
			output.flush()

	proj = EMAN2_cppwrap.EMData()
	proj.read_image(proj_stack,list_proj[0])

	size = proj.get_xsize()
	if size != proj.get_ysize():
		sp_global_def.sxprint("Image projections must be square!")
		return None

	overall_start = time.time()
	for i in range(niter):
		iter_start = time.time()
		mults = nimages*[0]
		for j in range(nimages):
			imgid = random.randint(0,nimages-1)
			mults[imgid]=mults[imgid]+1

		if CTF:
			params = {"size":size, "npad":npad, "symmetry":symmetry, "snr":snr, "sign":sign}
			r = EMAN2_cppwrap.Reconstructors.get("nn4_ctf", params);
		else:
			params = {"size":size, "npad":npad, "symmetry":symmetry, "snr":snr}
			r = EMAN2_cppwrap.Reconstructors.get("nn4", params);

		r.setup()
		if not(output is None):
			output.write( "Bootstrap volume %8d " % i )
			output.flush()

		store.restart()

		if not(output is None):
			output.write( "Inserting images " )
			output.flush()

		for j in range(nimages):
			if mults[j] > 0 :
				img_j = EMAN2_cppwrap.EMData()
				store.get_image( j, img_j );
				phi_j = img_j.get_attr( "phi" )
				tht_j = img_j.get_attr( "theta" )
				psi_j = img_j.get_attr( "psi" )
				tra_j = EMAN2_cppwrap.Transform( {"type":"spider", "phi":phi_j, "theta":tht_j, "psi":psi_j} )

				if CTF:
					cs = img_j.get_attr( "Cs" )
					pixel   = img_j.get_attr( "Pixel_size" )
					defocus = img_j.get_attr( "defocus" )
					voltage = img_j.get_attr( "voltage" )
					ampcont = img_j.get_attr( "amp_contrast" )
					bfactor = 0.0
					sp_utilities.set_ctf( img_j, [defocus, cs, voltage, pixel, bfactor, ampcont] )

				r.insert_slice(img_j, tra_j, mults[j])

				#[mean,sigma,min,max]= Util.infomask(img_j, None, False)
				#output.write( "img %4d %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f\n" % (j, mean, sigma, min, max, phi, theta, psi) )
				#output.flush()

		if not(output is None):
			output.write( "Finishing... " )
			output.flush( )

		vol = r.finish(True)

		if not(output is None):
			output.write( "Writing... " )
			output.flush()

		vol.write_image(volume_stack,i)

		if not(output is None):
			output.write( " done!" )
			output.write( " time %15.3f %15.3f \n" % (time.time()-iter_start,time.time()-overall_start) )
			output.flush()


def recons3d_em(projections_stack, max_iterations_count = 100, radius = -1, min_avg_abs_voxel_change = 0.01, use_weights = False, symmetry = "c1"):
	"""
	Reconstruction algorithm basing on the Expectation Maximization method
		projections_stack            -- file or list with projections
		max_iterations_count         -- stop criterion 
		min_avg_abs_voxel_change     -- stop criterion 
		use_weights                  -- true == multiply projections by extra weights
		symmetry                     -- type of symmetry
	#
	"""
	pass#IMPORTIMPORTIMPORT from time import clock
	pass#IMPORTIMPORTIMPORT from sp_utilities import model_blank, model_circle, model_square
	pass#IMPORTIMPORTIMPORT from sp_morphology import threshold_to_minval
	pass#IMPORTIMPORTIMPORT import types
	min_allowed_divisor = 0.0001

	if type(projections_stack) is bytes:
		projections = EMAN2_cppwrap.EMData.read_images(projections_stack)
	else:
		projections = projections_stack

	if len(projections) == 0:
		sp_global_def.ERROR("Stack of projections cannot be empty", "recons3d_em")

	nx = projections[0].get_xsize()
	if (projections[0].get_ysize() != nx) or (projections[0].get_zsize() != 1):
		sp_global_def.ERROR("This procedure works only for square images", "recons3d_em")

	if radius < 0:  radius = nx // 2 - 1
	sphere2D = sp_utilities.model_circle(radius, nx, nx)   
	sphere3D = sp_utilities.model_circle(radius, nx, nx, nx)
	solution = sp_utilities.model_blank(nx, nx, nx)
	a = sp_utilities.model_blank(nx, nx, nx) # normalization volume
	e2D = sp_utilities.model_square(nx, nx, nx)
	sphere3D_volume = sp_utilities.model_blank(nx,nx,nx).cmp("lod",sphere3D,{"negative":0,"normalize":0})
	#print "Parameters:  size=%d  radius=%d  projections_count=%d  max_iterations_count=%d min_avg_abs_voxel_change=%f" % (
	#					nx, radius, len(projections), max_iterations_count, min_avg_abs_voxel_change )

	# ----- create initial solution, calculate weights and normalization image (a)
	projections_angles = []  # list of lists of angles
	projections_data   = []  # list of lists of projections' images with weights
	for proj in projections:
		angles = [] # list of angles
		data = []   # list of projections' images with weights
		RA = proj.get_attr( "xform.projection" )
		EMAN2_cppwrap.Util.mul_img( proj, sphere2D )
		for j in range(RA.get_nsym(symmetry)):
			angdict = RA.get_sym_sparx(symmetry,j).get_rotation("spider") 
			angles.append( [angdict["phi"], angdict["theta"], angdict["psi"]] )
			chao_params = {"anglelist":angles[j],"radius":radius}
			EMAN2_cppwrap.Util.add_img( solution, proj.backproject("chao", chao_params) )
			EMAN2_cppwrap.Util.add_img( a, e2D.backproject("chao", chao_params) )
			if use_weights:
				proj3Dsphere = sphere3D.project("chao", chao_params)
				EMAN2_cppwrap.Util.mul_scalar( proj3Dsphere, 1.0 / EMAN2_cppwrap.Util.infomask(proj3Dsphere, None, True)[3] )
				EMAN2_cppwrap.Util.mul_img( proj, proj3Dsphere )
			data.append(proj)
		projections_angles.append(angles)
		projections_data.append(data)
	a = sp_morphology.threshold_to_minval(a, min_allowed_divisor)  # make sure that voxels' values are not too small (image a is divisior)
	EMAN2_cppwrap.Util.mul_img( solution, sphere3D )
	EMAN2_cppwrap.Util.div_img( solution, a )
	#print "Projections loading COMPLETED"
	# ----- iterations
	prev_avg_absolute_voxel_change = 999999999.0
	time_projection = 0.0
	time_backprojection = 0.0
	time_iterations = time.clock()
	for iter_no in range(max_iterations_count):
		q = sp_utilities.model_blank(nx, nx, nx)
		for i in range(len(projections_angles)):
			for j in range(len(projections_angles[i])):
				chao_params = {"anglelist":projections_angles[i][j],"radius":radius}
				time_start = time.clock()
				w = solution.project("chao", chao_params)
				time_projection += time.clock() - time_start
				p = projections_data[i][j] / sp_morphology.threshold_to_minval(w, min_allowed_divisor)
				time_start = time.clock()
				q += p.backproject("chao", chao_params)
				time_backprojection += time.clock() - time_start
		EMAN2_cppwrap.Util.div_img( q, a )
		EMAN2_cppwrap.Util.mul_img( q, solution ) # q <- new solution  
		avg_absolute_voxel_change = q.cmp("lod",solution,{"mask":sphere3D,"negative":0,"normalize":0}) / sphere3D_volume
		if avg_absolute_voxel_change > prev_avg_absolute_voxel_change:
			#print "Finish and return last good solution"
			break
		prev_avg_absolute_voxel_change = avg_absolute_voxel_change
		solution = q
		#print "Iteration ", iter_no, ",  avg_abs_voxel_change=", avg_absolute_voxel_change 
		if min_avg_abs_voxel_change > avg_absolute_voxel_change:
			break
	time_iterations = time.clock() - time_iterations
	# ----- return solution and exit
	#print "Times: iterations=", time_iterations, "  project=", time_projection, "  backproject=", time_backprojection
	return solution


def recons3d_em_MPI(projections_stack, output_file, max_iterations_count = 100, radius = -1, min_norm_absolute_voxel_change = 0.01, use_weights = False, symmetry = "c1", min_norm_squared_voxel_change = 0.0001):
	"""
	Reconstruction algorithm basing on the Expectation Maximization method.
		projections_stack              -- file or list with projections
		max_iterations_count           -- stop criterion 
		min_norm_absolute_voxel_change -- stop criterion (set -1 to switch off) 
		min_norm_squared_voxel_change  -- stop criterion (set -1 to switch off)
		use_weights                    -- true == multiply projections by extra weights
		symmetry                       -- type of symmetry
	#
	"""
	pass#IMPORTIMPORTIMPORT from time import clock
	pass#IMPORTIMPORTIMPORT from sp_utilities import model_blank, model_circle, model_square, circumference
	pass#IMPORTIMPORTIMPORT from sp_morphology import threshold_to_minval
	pass#IMPORTIMPORTIMPORT import types
	pass#IMPORTIMPORTIMPORT from string import replace
	pass#IMPORTIMPORTIMPORT from mpi import mpi_comm_size, mpi_comm_rank, MPI_COMM_WORLD
	pass#IMPORTIMPORTIMPORT from sp_utilities import reduce_EMData_to_root, bcast_EMData_to_all, bcast_number_to_all, send_EMData, recv_EMData
	min_allowed_divisor = 0.0001

	mpi_n = mpi.mpi_comm_size(mpi.MPI_COMM_WORLD)
	mpi_r = mpi.mpi_comm_rank(mpi.MPI_COMM_WORLD)

	# ----- read projections 
	if type(projections_stack) is bytes:
		all_projs_count = EMAN2_cppwrap.EMUtil.get_image_count(projections_stack)
	else:
		all_projs_count = len(projections_stack)

	if all_projs_count < mpi_n:
		sp_global_def.ERROR("Number of projections cannot be less than number of MPI processes", "recons3d_em")

	projs_begin = (mpi_r * all_projs_count) // mpi_n
	projs_end = ((mpi_r+1) * all_projs_count) // mpi_n

	if type(projections_stack) is bytes:
		projections = EMAN2_cppwrap.EMData.read_images(projections_stack, list(range(projs_begin,projs_end)))
	else:
		#projections = projections_stack[projs_begin:projs_end]
		projections = projections_stack
	# ----------------------------------------------

	nx = projections[0].get_xsize()
	if (projections[0].get_ysize() != nx) or (projections[0].get_zsize() != 1):
		sp_global_def.ERROR("This procedure works only for square images", "recons3d_em")

	if radius < 0: radius = nx // 2 - 1
	sphere2D = sp_utilities.model_circle(radius, nx, nx)   
	sphere3D = sp_utilities.model_circle(radius, nx, nx, nx)
	solution = sp_utilities.model_blank(nx, nx, nx)
	a = sp_utilities.model_blank(nx, nx, nx) # normalization volume
	e2D = sp_utilities.model_square(nx, nx, nx)
	if mpi_r == 0:
		sp_global_def.sxprint("MPI processes: ", mpi_n)
		sp_global_def.sxprint("Parameters:  size=%d  radius=%d  projections_count=%d  max_iterations_count=%d min_norm_absolute_voxel_change=%f" % (
						nx, radius, all_projs_count, max_iterations_count, min_norm_absolute_voxel_change ))	

	# ----- create initial solution, calculate weights and normalization image (a)
	projections_angles = []  # list of lists of angles
	projections_data   = []  # list of lists of projections' images with weights
	for proj in projections:
		angles = [] # list of angles
		data = []   # list of projections' images with weights
		RA = proj.get_attr( "xform.projection" )
		EMAN2_cppwrap.Util.mul_img( proj, sphere2D )
		for j in range(RA.get_nsym(symmetry)):
			angdict = RA.get_sym_sparx(symmetry,j).get_rotation("spider") 
			angles.append( [angdict["phi"], angdict["theta"], angdict["psi"]] )
			chao_params = {"anglelist":angles[j],"radius":radius}
			EMAN2_cppwrap.Util.add_img( solution, proj.backproject("chao", chao_params) )
			EMAN2_cppwrap.Util.add_img( a, e2D.backproject("chao", chao_params) )
			if use_weights:
				proj3Dsphere = sphere3D.project("chao", chao_params)
				EMAN2_cppwrap.Util.mul_scalar( proj3Dsphere, 1.0 / EMAN2_cppwrap.Util.infomask(proj3Dsphere, None, True)[3] )
				EMAN2_cppwrap.Util.mul_img( proj, proj3Dsphere )
			data.append(proj)
		projections_angles.append(angles)
		projections_data.append(data)
	# reduce_scatter(solution)
	sp_utilities.reduce_EMData_to_root(solution, mpi_r)
	sp_utilities.bcast_EMData_to_all  (solution, mpi_r)
	# reduce_scatter(a)
	sp_utilities.reduce_EMData_to_root(a, mpi_r)
	sp_utilities.bcast_EMData_to_all  (a, mpi_r)
	# ------------------------
	a = sp_morphology.threshold_to_minval(a, min_allowed_divisor)  # make sure that voxels' values are not too small (image a is divisior)
	EMAN2_cppwrap.Util.mul_img( solution, sphere3D )
	EMAN2_cppwrap.Util.div_img( solution, a )
	if mpi_r == 0: sp_global_def.sxprint("Projections loading COMPLETED")
	# ----- iterations
	prev_avg_absolute_voxel_change = 999999999.0
	time_projection = 0.0
	time_backprojection = 0.0
	time_iterations = time.clock()
	for iter_no in range(max_iterations_count):
		q = sp_utilities.model_blank(nx, nx, nx)
		for i in range(len(projections_angles)):
			for j in range(len(projections_angles[i])):
				chao_params = {"anglelist":projections_angles[i][j],"radius":radius}
				time_start = time.clock()
				w = solution.project("chao", chao_params)
				time_projection += time.clock() - time_start
				p = projections_data[i][j] / sp_morphology.threshold_to_minval(w, min_allowed_divisor)
				time_start = time.clock()
				q += p.backproject("chao", chao_params)
				time_backprojection += time.clock() - time_start
		# reduce_scatter(q)
		sp_utilities.reduce_EMData_to_root(q, mpi_r)
		sp_utilities.bcast_EMData_to_all  (q, mpi_r)
		# ----------------------
		EMAN2_cppwrap.Util.div_img( q, a )
		EMAN2_cppwrap.Util.mul_img( q, solution ) # q <- new solution  
		norm_absolute_voxel_change = q.cmp("lod",solution,{"mask":sphere3D,"negative":0,"normalize":0}) / q.cmp("lod",sp_utilities.model_blank(nx,nx,nx),{"mask":sphere3D,"negative":0,"normalize":0})
		norm_squared_voxel_change  = q.cmp("sqEuclidean",solution,{"mask":sphere3D}) / q.cmp("sqEuclidean",sp_utilities.model_blank(nx,nx,nx),{"mask":sphere3D})
		if norm_absolute_voxel_change > prev_avg_absolute_voxel_change:
			if mpi_r == 0: sp_global_def.sxprint("Finish and return last good solution")
			break
		prev_avg_absolute_voxel_change = norm_absolute_voxel_change
		solution = q
		solution = sp_utilities.circumference(solution, radius-2, radius)
		if (iter_no+1)%5 == 0 and mpi_r == 0:
			solution.write_image(string.replace(output_file, ".hdf", "_%03d.hdf"%(iter_no+1)))
		if mpi_r == 0: sp_global_def.sxprint("Iteration ", iter_no+1, ",  norm_abs_voxel_change=", norm_absolute_voxel_change, ",  norm_squared_voxel_change=", norm_squared_voxel_change) 
		if min_norm_absolute_voxel_change > norm_absolute_voxel_change or min_norm_squared_voxel_change > norm_squared_voxel_change:
			break
	time_iterations = time.clock() - time_iterations
	# ----- return solution and exit
	if mpi_r == 0: sp_global_def.sxprint("Times: iterations=", time_iterations, "  project=", time_projection, "  backproject=", time_backprojection)
	return solution


def recons3d_sirt(stack_name, list_proj, radius, lam=1.0e-4, maxit=100, symmetry="c1", tol=0.001):
	"""
	SIRT
		
		tol   -- convergence tolerance
		lam   -- damping parameter
		maxit -- maximum number of iterations
	#
	"""
	pass#IMPORTIMPORTIMPORT from math import sqrt
	pass#IMPORTIMPORTIMPORT from sp_utilities import model_circle
	#  analyze the symmetries Phil's code has all symmetries ready...
	nsym=1

	#  get image size from the first image
	data = EMAN2_cppwrap.EMData()
	data.read_image(stack_name,list_proj[0])
	nx = data.get_xsize()
	mask2d=sp_utilities.model_circle(radius,nx,nx)  # SIRT works for squares only!
	mask2d = 1.0 - mask2d  # invert the mask to get average in corners
	nangles = len(list_proj)
	#
	mask3d=sp_utilities.model_circle(radius,nx,nx,nx) # a 3D mask for error calculation
	#
	# create a volume to hold the reconstruction 
	#
	xvol = EMAN2_cppwrap.EMData()
	xvol.set_size(nx,nx,nx)
	xvol.to_zero()
	#
	# create a volume to hold trans(P)*P*xvol
	#
	pxvol = xvol.copy()
	#  array of symmetrized angles
	symangles=3*[0.0]
	angles = []

	# start iterating
	iter  = 1
	while iter <= maxit:
		if (iter == 1):
			#
			# backproject 2D images first to create the right-hand side of the
			# the normal equation
			#
			bvol = EMAN2_cppwrap.EMData()
			bvol.set_size(nx,nx,nx)
			bvol.to_zero()
			for i in range(nangles):
				# read projections and do initial  backprojection
				data.read_image(stack_name,list_proj[i])
				stat = EMAN2_cppwrap.Util.infomask(data, mask2d, False)
				data = data-stat[0]   # subtract the background average in the corners
				
				RA = data.get_attr( "xform.projection" )

				angles.append(RA)
				#ATTENTION
				#for transform in Symmetry3D.get_symmetries(symmetry):
					#Tf = transform*RA
					# though why do you go from 1 to nysm? why not 0 to nysm-1 ? It should be
					# equivalent unless I am missing something
					#angdict = Tf.get_params("spider")
					# then continue as before
				for ns in range(1,nsym+1):
					# multiply myangles by symmetry using Phil's Transform class
					Tf=RA.get_sym_sparx(symmetry,ns) #
					angdict = Tf.get_rotation("spider")
					#   Chao - please check the order of phi, theta, psi
					symangles[0] = angdict["phi"]
					symangles[1] = angdict["theta"]
					symangles[2] = angdict["psi"]
					myparams = {"anglelist":symangles, "radius":radius}
					bvol += data.backproject("chao", myparams)
			old_rnorm = bnorm = math.sqrt(bvol.cmp("dot",bvol,{"mask":mask3d,"negative":0}))
			grad  = bvol
		else:
			#  Insert your favorite MPI here
			pxvol.to_zero() 
			for i in range(nangles):
				# just a single slice of phi, theta, psi
				RA = angles[i]
				for ns in range(1,nsym+1):
					# multiply myangles by symmetry using Phil's Transform class
					Tf = RA.get_sym_sparx(symmetry,ns)#Tf.get_rotation()
					angdict = Tf.get_rotation("spider")
					#				    Chao - please check the order of phi, theta, psi
					symangles[0] = angdict["phi"]
					symangles[1] = angdict["theta"]
					symangles[2] = angdict["psi"]
					myparams = {"anglelist":symangles, "radius":radius}
					data  = xvol.project("chao", myparams) 
					pxvol += data.backproject("chao",myparams)
			grad  = bvol - pxvol

		rnorm = math.sqrt(grad.cmp("dot",grad,{"mask":mask3d,"negative":0}))
		sp_global_def.sxprint('iter = %3d,  rnorm = %6.3f,  rnorm/bnorm = %6.3f' % (iter,rnorm,rnorm/bnorm))
		if (rnorm < tol or rnorm > old_rnorm): break
		old_rnorm = rnorm
		xvol = xvol + lam*grad
		iter = iter + 1

	return  xvol

def recons3d_wbp(stack_name, list_proj, method = "general", const=1.0E4, symmetry="c1", radius=None): 
	"""
		Weigthed back-projection algorithm.
		stack_name - disk stack with projections or in-core list of images
		list_proj  - list of projections to be used in the reconstruction
		method - "general" Radermacher's Gaussian, "exact" MvHs triangle
		const  - for "general" 1.0e4 works well, for "exact" it should be the diameter of the object
		symmetry - point group symmetry of the object
	""" 
	pass#IMPORTIMPORTIMPORT import types
	pass#IMPORTIMPORTIMPORT from sp_utilities import get_im

	if type(stack_name) == bytes:
		B = EMAN2_cppwrap.EMData()
		B.read_image(stack_name,list_proj[0])
	else : B = stack_name[list_proj[0]].copy()

	ny = B.get_ysize()  # have to take ysize, because xsize is different for real and fft images

	if radius == None: radius = (ny - 1) // 2

	CUBE = EMAN2_cppwrap.EMData()
	CUBE.set_size(ny, ny, ny)
	CUBE.to_zero()

	nsym = EMAN2_cppwrap.Transform.get_nsym(symmetry)
	nimages = len(list_proj)

	ss = [0.0]*(6*nsym*nimages)
	symmetry_transforms = [ [None for i in range(nsym)] for j in range(nimages) ] # list of nimages lists of nsym elements
	for iProj in range(nimages):
		if type(stack_name) == bytes:
			B.read_image(stack_name,list_proj[iProj], True)
		else:
			B = stack_name[list_proj[iProj]]
		transform = B.get_attr("xform.projection")
		for iSym in range(nsym):
			symmetry_transforms[iProj][iSym] = transform.get_sym_sparx(symmetry, iSym)
			d = symmetry_transforms[iProj][iSym].get_params("spider")
			DMnSS = EMAN2_cppwrap.Util.CANG(d["phi"], d["theta"], d["psi"])
			ss[ (iProj*nsym+iSym)*6 : (iProj*nsym+iSym+1)*6 ] = DMnSS["SS"]

	if method=="exact":    
		const = int(const)

	for iProj in range(nimages):
		proj = sp_utilities.get_im(stack_name, list_proj[iProj])
		for iSym in range(nsym):
			B = proj.copy()
			B.set_attr("xform.projection", symmetry_transforms[iProj][iSym])
			if   method=="general":  EMAN2_cppwrap.Util.WTF(B, ss, const, iProj*nsym+iSym+1)  # counting in WTF start from 1!
			elif method=="exact"  :  EMAN2_cppwrap.Util.WTM(B, ss, const, iProj*nsym+iSym+1)  # counting in WTM start from 1!
			EMAN2_cppwrap.Util.BPCQ(B, CUBE, radius)

	return CUBE


def recons3d_vwbp(stack_name, list_proj, method = "general", const=1.0E4, symmetry="c1",outstack="bdb:temp"): 
	"""
		Weigthed back-projection algorithm.
		stack_name - disk stack with projections or in-core list of images
		list_proj  - list of projections to be used in the reconstruction
		method - "general" Radermacher's Gaussian, "exact" MvHs triangle
		const  - for "general" 1.0e4 works well, for "exact" it should be the diameter of the object
		symmetry - point group symmetry of the object

		WARNING - symmetries not implemented!!!!!!!!!
	""" 
	pass#IMPORTIMPORTIMPORT import types
	pass#IMPORTIMPORTIMPORT from sp_utilities import get_im

	if type(stack_name) == bytes:
		B = EMAN2_cppwrap.EMData()
		B.read_image(stack_name,list_proj[0])
	else : B = stack_name[list_proj[0]].copy()

	ny = B.get_ysize()  # have to take ysize, because xsize is different for real and fft images

	nsym = EMAN2_cppwrap.Transform.get_nsym(symmetry)
	nimages = len(list_proj)

	ss = [0.0]*(6*nsym*nimages)
	symmetry_transforms = [ [None for i in range(nsym)] for j in range(nimages) ] # list of nimages lists of nsym elements
	for iProj in range(nimages):
		if type(stack_name) == bytes:
			B.read_image(stack_name,list_proj[iProj], True)
		else:
			B = stack_name[list_proj[iProj]]
		transform = B.get_attr("xform.projection")
		for iSym in range(nsym):
			symmetry_transforms[iProj][iSym] = transform.get_sym_sparx(symmetry, iSym)
			d = symmetry_transforms[iProj][iSym].get_params("spider")
			DMnSS = EMAN2_cppwrap.Util.CANG(d["phi"], d["theta"], d["psi"])
			ss[ (iProj*nsym+iSym)*6 : (iProj*nsym+iSym+1)*6 ] = DMnSS["SS"]

	if method=="exact":    
		const = int(const)

	for iProj in range(nimages):
		if(iProj%100 == 0):  sp_global_def.sxprint("BPCQ  ",iProj)
		CUBE = EMAN2_cppwrap.EMData()
		CUBE.set_size(ny, ny, ny)
		CUBE.to_zero()
		proj = sp_utilities.get_im(stack_name, list_proj[iProj])
		for iSym in range(nsym):
			B = proj.copy()
			B.set_attr("xform.projection", symmetry_transforms[iProj][iSym])
			if   method=="general":  EMAN2_cppwrap.Util.WTF(B, ss, const, iProj*nsym+iSym+1)  # counting in WTF start from 1!
			elif method=="exact"  :  EMAN2_cppwrap.Util.WTM(B, ss, const, iProj*nsym+iSym+1)  # counting in WTM start from 1!
			pass#IMPORTIMPORTIMPORT from sp_filter import filt_tanl
			B = sp_filter.filt_tanl(B, 0.3, 0.1)
			EMAN2_cppwrap.Util.BPCQ(B, CUBE, (B.get_ysize()-1)//2)
		CUBE.write_image(outstack, iProj)
	return CUBE


def prepare_wbp(stack_name, list_proj, method = "general", const=1.0E4, symmetry="c1"):
	"""
		Prepare auxiliary arrays dm and ss.
		Weigthed back-projection algorithm.
		stack_name - disk stack with projections or in-core list of images
		list_proj  - list of projections to be used in the reconstruction
		method - "general" Radermacher's Gaussian, "exact" MvHs triangle
		const  - for "general" 1.0e4 works well, for "exact" it should be the diameter of the object
		symmetry - point group symmetry of the object
	"""
	pass#IMPORTIMPORTIMPORT import types

	if type(stack_name) == bytes:
		B = EMAN2_cppwrap.EMData()
		B.read_image(stack_name,list_proj[0])
	else : B = stack_name[list_proj[0]].copy()

	nx = B.get_xsize()

	RA = EMAN2_cppwrap.Transform()
	nsym = RA.get_nsym(symmetry)

	nimages = len(list_proj)
	ntripletsWnsym = nsym*nimages
	dm=[0.0]*(9*ntripletsWnsym)
	ss=[0.0]*(6*ntripletsWnsym)
	count = 0
	pass#IMPORTIMPORTIMPORT from sp_utilities import get_params_proj
	for i in range(nimages):
		if type(stack_name) == bytes:
			B.read_image(stack_name,list_proj[i], True)
			PHI, THETA, PSI, s2x, s2y = sp_utilities.get_params_proj( B )
		else:  
			PHI, THETA, PSI, s2x, s2y = sp_utilities.get_params_proj( stack_name[list_proj[i]] )
		DMnSS = EMAN2_cppwrap.Util.CANG(PHI,THETA,PSI)
		dm[(count*9) :(count+1)*9] = DMnSS["DM"]
		ss[(count*6) :(count+1)*6] = DMnSS["SS"]
		count += 1
	return dm,ss


def recons3d_swbp(A, transform, L, ss, method = "general", const=1.0E4, symmetry="c1"):
	"""
	    Take one projection, but angles form the entire set.  Build the weighting function for the given projection taking into account all,
		apply it, and backproject.
		The projection number is L counting from zero.
		Weigthed back-projection algorithm.
		stack_name - disk stack with projections or in-core list of images
		method - "general" Radermacher's Gaussian, "exact" MvHs triangle
		const  - for "general" 1.0e4 works well, for "exact" it should be the diameter of the object
		symmetry - point group symmetry of the object
	"""
	B = A.copy()
	nx = B.get_xsize()
	if(method=="exact"  ):    const = int(const)
	nsym = 1
	CUBE = EMAN2_cppwrap.EMData()
	CUBE.set_size(nx, nx, nx)
	CUBE.to_zero()

	org_transform = B.get_attr("xform.projection")

	count = 0
	for j in range(1):
		count += 1   # It has to be there as counting in WTF and WTM start from 1!
		if   (method=="general"):    EMAN2_cppwrap.Util.WTF(B, ss, const, L+1)
		elif (method=="exact"  ):    EMAN2_cppwrap.Util.WTM(B, ss, const, L+1)

		B.set_attr("xform.projection", transform)
		EMAN2_cppwrap.Util.BPCQ(B, CUBE, (B.get_ysize()-1)//2)
		
	B.set_attr("xform.projection", org_transform)
	return CUBE, B

def weight_swbp(A, L, ss, method = "general", const=1.0E4, symmetry="c1"):
	"""
	    Take one projection, but angles form the entire set.  Build the weighting function for the given projection taking into account all,
		apply it, return weighted projection.
		The projection number is L counting from zero.
		Weigthed back-projection algorithm.
		stack_name - disk stack with projections or in-core list of images
		method - "general" Radermacher's Gaussian, "exact" MvHs triangle
		const  - for "general" 1.0e4 works well, for "exact" it should be the diameter of the object
		symmetry - point group symmetry of the object
	"""

	if(method=="exact"  ):    const = int(const)
	nsym = 1
	B = A.copy()
	count = 0
	for j in range(1):
		count += 1   # It has to be there as counting in WTF and WTM start from 1!
		if   (method=="general"):    EMAN2_cppwrap.Util.WTF(B, ss, const, L+1)
		elif (method=="exact"  ):    EMAN2_cppwrap.Util.WTM(B, ss, const, L+1)

	return B

def backproject_swbp(B, transform = None, symmetry="c1"): 
	"""
	    Take one projection, but angles form the entire set.  Build the weighting function for the given projection taking into account all,
		apply it, and backproject.
		The projection number is L counting from zero.
		Weigthed back-projection algorithm.
		stack_name - disk stack with projections or in-core list of images
		method - "general" Radermacher's Gaussian, "exact" MvHs triangle
		const  - for "general" 1.0e4 works well, for "exact" it should be the diameter of the object
		symmetry - point group symmetry of the object
	""" 

	ny = B.get_ysize()
	CUBE = EMAN2_cppwrap.EMData()
	CUBE.set_size(ny, ny, ny)
	CUBE.to_zero()

	org_transform = B.get_attr("xform.projection")
	if transform != None:
		B.set_attr("xform.projection", transform)
	EMAN2_cppwrap.Util.BPCQ(B, CUBE, (B.get_ysize()-1)//2)
	B.set_attr("xform.projection", org_transform)

	return CUBE

def one_swbp(CUBE, B, transform = None, symmetry="c1"): 
	"""
	        Take one projection, but angles form the entire set.  Build the weighting function for the given projection taking into account all,
		apply it, and backproject.
		The projection number is L counting from zero.
		method - "general" Rademacher's Gaussian, "exact" MvHs triangle
		const  - for "general" 1.0e4 works well, for "exact" it should be the diameter of the object
		symmetry - point group symmetry of the object
	""" 
	org_transform = B.get_attr("xform.projection")
	if transform != None:
		B.set_attr("xform.projection", transform)
	EMAN2_cppwrap.Util.BPCQ(B, CUBE, (B.get_ysize()-1)//2)  
	B.set_attr("xform.projection", org_transform)

def prepare_recons(data, symmetry, myid, main_node_half, half_start, step, index, finfo=None, npad = 2, mpi_comm=None):
	pass#IMPORTIMPORTIMPORT from random     import randint
	pass#IMPORTIMPORTIMPORT from sp_utilities  import reduce_EMData_to_root
	pass#IMPORTIMPORTIMPORT from mpi        import mpi_barrier, MPI_COMM_WORLD
	pass#IMPORTIMPORTIMPORT from EMAN2 import Reconstructors

	if mpi_comm == None:
		mpi_comm = mpi.MPI_COMM_WORLD

	nx = data[0].get_xsize()

	fftvol_half = EMAN2_cppwrap.EMData()
	weight_half = EMAN2_cppwrap.EMData()
	half_params = {"size":nx, "npad":npad, "symmetry":symmetry, "fftvol":fftvol_half, "weight":weight_half}
	half = EMAN2_cppwrap.Reconstructors.get( "nn4", half_params )
	half.setup()

	group = -1
	for i in range(half_start, len(data), step):
		if(index >-1 ):  group = data[i].get_attr('group')
		if(group == index):
			# horatio active_refactoring Jy51i1EwmLD4tWZ9_00000_1
			# if( data[i].get_attr_default('active',1) == 1):
			# 	xform_proj = data[i].get_attr( "xform.projection" )
			# 	half.insert_slice(data[i], xform_proj )
			xform_proj = data[i].get_attr( "xform.projection" )
			half.insert_slice(data[i], xform_proj )

	if not(finfo is None):
		finfo.write( "begin reduce half\n" )
		finfo.flush()

	sp_utilities.reduce_EMData_to_root(fftvol_half, myid, main_node_half, mpi_comm)
	sp_utilities.reduce_EMData_to_root(weight_half, myid, main_node_half, mpi_comm)

	if not(finfo is None):
		finfo.write( "after reduce half\n" )
		finfo.flush()

	if myid == main_node_half:
		tmpid = random.randint(0, 1000000)
		fftvol_half_file = ("fftvol_half%d.hdf" % tmpid)
		weight_half_file = ("weight_half%d.hdf" % tmpid)
		fftvol_half.write_image(fftvol_half_file)
		weight_half.write_image(weight_half_file)
	mpi.mpi_barrier(mpi_comm)

	fftvol_half = None
	weight_half = None

	if myid == main_node_half:  return fftvol_half_file, weight_half_file

	return None, None
































def prepare_recons_ctf(nx, data, snr, symmetry, myid, main_node_half, half_start, step, finfo=None, npad = 2, mpi_comm=None, smearstep = 0.0):
	pass#IMPORTIMPORTIMPORT from random     import randint
	pass#IMPORTIMPORTIMPORT from sp_utilities  import reduce_EMData_to_root
	pass#IMPORTIMPORTIMPORT from mpi        import mpi_barrier, MPI_COMM_WORLD
	pass#IMPORTIMPORTIMPORT from EMAN2 import Reconstructors

	if mpi_comm == None:
		mpi_comm = mpi.MPI_COMM_WORLD

	fftvol_half = EMAN2_cppwrap.EMData()

	if( smearstep > 0.0 ):
		#if myid == 0:  print "  Setting smear in prepare_recons_ctf"
		ns = 1
		smear = []
		for j in range(-ns,ns+1):
			if( j != 0):
				for i in range(-ns,ns+1):
					for k in range(-ns,ns+1):
						smear += [i*smearstep,j*smearstep,k*smearstep,1.0]
		# Deal with theta = 0.0 cases
		prj = []
		for i in range(-ns,ns+1):
			for k in range(-ns,ns+1):
				prj.append(i+k)
		for i in range(-2*ns,2*ns+1,1):
			smear += [i*smearstep,0.0,0.0,float(prj.count(i))]
		#if myid == 0:  print "  Smear  ",smear
		fftvol_half.set_attr("smear", smear)

	weight_half = EMAN2_cppwrap.EMData()
	half_params = {"size":nx, "npad":npad, "snr":snr, "sign":1, "symmetry":symmetry, "fftvol":fftvol_half, "weight":weight_half}
	half = EMAN2_cppwrap.Reconstructors.get( "nn4_ctf", half_params )
	half.setup()

	for i in range(half_start, len(data), step):
		xform_proj = data[i].get_attr( "xform.projection" )
		half.insert_slice(data[i], xform_proj )

	if not(finfo is None):
		finfo.write( "begin reduce half\n" )
		finfo.flush()

	sp_utilities.reduce_EMData_to_root(fftvol_half, myid, main_node_half, mpi_comm)
	sp_utilities.reduce_EMData_to_root(weight_half, myid, main_node_half, mpi_comm)

	if not(finfo is None):
		finfo.write( "after reduce half\n" )
		finfo.flush()

	if myid == main_node_half:
		tmpid = random.randint(0, 1000000) 
		fftvol_half_file = ("fftvol_half%d.hdf" % tmpid)
		weight_half_file = ("weight_half%d.hdf" % tmpid)
		fftvol_half.write_image(fftvol_half_file)
		weight_half.write_image(weight_half_file)
	mpi.mpi_barrier(mpi_comm)

	fftvol_half = None
	weight_half = None

	if myid == main_node_half:
		return fftvol_half_file, weight_half_file

	return None,None


def recons_from_fftvol(size, fftvol, weight, symmetry, npad = 2):
	pass#IMPORTIMPORTIMPORT from EMAN2 import Reconstructors

	params = {"size":size, "npad":npad, "symmetry":symmetry, "fftvol":fftvol, "weight":weight}
	r = EMAN2_cppwrap.Reconstructors.get("nn4", params)
	r.setup()
	dummy = r.finish(True)
	return fftvol


def recons_ctf_from_fftvol(size, fftvol, weight, snr, symmetry, weighting=1, npad = 2):
	pass#IMPORTIMPORTIMPORT from EMAN2 import Reconstructors

	params = {"size":size, "npad":npad, "snr":snr, "sign":1, "symmetry":symmetry, "fftvol":fftvol, "weight":weight, "weighting":weighting}
	r = EMAN2_cppwrap.Reconstructors.get("nn4_ctf", params)
	r.setup()
	dummy = r.finish(True)
	return fftvol

def recons_ctf_from_fftvol_using_nn4_ctfw(size, fftvol, weight, snr, symmetry, weighting=1, npad = 2):
	pass#IMPORTIMPORTIMPORT from EMAN2 import Reconstructors

	params = {"size":size, "npad":npad, "snr":snr, "sign":1, "symmetry":symmetry, "fftvol":fftvol, "weight":weight, "weighting":weighting}
	# r = Reconstructors.get("nn4_ctf", params)
	r = EMAN2_cppwrap.Reconstructors.get("nn4_ctfw", params)
	r.setup()
	dummy = r.finish(True)
	return fftvol


def get_image_size( imgdata, myid ):
	pass#IMPORTIMPORTIMPORT from mpi import mpi_gather, mpi_bcast, MPI_COMM_WORLD, MPI_INT
	nimg = len(imgdata)

	nimgs = mpi.mpi_gather( nimg, 1, mpi.MPI_INT, 1, mpi.MPI_INT, 0, mpi.MPI_COMM_WORLD )

	if myid==0:
		src = -1
		for i in range( len(nimgs) ):
			if int(nimgs[i]) > 0 :
				src = i
				break
		if src==-1:
			return 0
	else:
		src = -1

	size_src = mpi.mpi_bcast( src, 1, mpi.MPI_INT, 0, mpi.MPI_COMM_WORLD )

	if myid==int(size_src[0]):
		assert nimg > 0
		size = imgdata[0].get_xsize()
	else:
		size = -1

	nx = mpi.mpi_bcast( size, 1, mpi.MPI_INT, size_src[0], mpi.MPI_COMM_WORLD )
	return int(nx[0])


def rec3D_MPI(data, snr = 1.0, symmetry = "c1", mask3D = None, fsc_curve = None, \
		myid = 0, main_node = 0, rstep = 1.0, odd_start=0, eve_start=1, finfo=None, \
		index=-1, npad = 2, mpi_comm=None, smearstep = 0.0):
	'''
	  This function is to be called within an MPI program to do a reconstruction on a dataset kept 
	  in the memory, computes reconstruction and through odd-even, in order to get the resolution
	'''
	pass#IMPORTIMPORTIMPORT import os
	pass#IMPORTIMPORTIMPORT from sp_statistics import fsc_mask
	pass#IMPORTIMPORTIMPORT from sp_utilities  import model_blank, model_circle, get_image, send_EMData, recv_EMData
	pass#IMPORTIMPORTIMPORT from mpi        import mpi_comm_size, MPI_COMM_WORLD
	
	if mpi_comm == None:
		mpi_comm = mpi.MPI_COMM_WORLD
	
	nproc = mpi.mpi_comm_size(mpi_comm)

	if nproc==1:
		assert main_node==0
		main_node_odd = main_node
		main_node_eve = main_node
		main_node_all = main_node
	elif nproc==2:
		main_node_odd = main_node
		main_node_eve = (main_node+1)%2
		main_node_all = main_node

		tag_voleve     = 1000
		tag_fftvol_eve = 1001
		tag_weight_eve = 1002
	else:
		#spread CPUs between different nodes to save memory
		main_node_odd = main_node
		main_node_eve = (int(main_node)+nproc-1)%int(nproc)
		main_node_all = (int(main_node)+nproc//2)%int(nproc)

		tag_voleve     = 1000
		tag_fftvol_eve = 1001
		tag_weight_eve = 1002

		tag_fftvol_odd = 1003
		tag_weight_odd = 1004
		tag_volall     = 1005


	if index != -1 :
		grpdata = []
		for i in range(len(data)):
			if data[i].get_attr('group') == index:
				grpdata.append(data[i])
		imgdata = grpdata
	else:
		imgdata = data

	nx = get_image_size(imgdata, myid)
	if nx == 0:
		sp_global_def.ERROR("Warning: no images were given for reconstruction, this usually means there is an empty group, returning empty volume", "rec3D", 0)
		return sp_utilities.model_blank( 2, 2, 2 ), None

	fftvol_odd_file, weight_odd_file = prepare_recons_ctf(nx, imgdata, snr, symmetry, myid, main_node_odd, odd_start, 2, finfo, npad, mpi_comm=mpi_comm, smearstep = smearstep)
	fftvol_eve_file, weight_eve_file = prepare_recons_ctf(nx, imgdata, snr, symmetry, myid, main_node_eve, eve_start, 2, finfo, npad, mpi_comm=mpi_comm, smearstep = smearstep)
	del imgdata

	if nproc == 1:
		fftvol = sp_utilities.get_image(fftvol_odd_file)
		weight = sp_utilities.get_image(weight_odd_file)
		volodd = recons_ctf_from_fftvol(nx, fftvol, weight, snr, symmetry, npad = npad)

		fftvol = sp_utilities.get_image(fftvol_eve_file)
		weight = sp_utilities.get_image(weight_eve_file)
		voleve = recons_ctf_from_fftvol(nx, fftvol, weight, snr, symmetry, npad = npad)

		if( not mask3D ):
			nx = volodd.get_xsize()
			ny = volodd.get_ysize()
			nz = volodd.get_zsize()
			mask3D = sp_utilities.model_circle(min(nx,ny,nz)//2 - 2, nx,ny,nz)
		fscdat = sp_statistics.fsc_mask( volodd, voleve, mask3D, rstep, fsc_curve)
		del  volodd, voleve, mask3d

		fftvol = sp_utilities.get_image( fftvol_odd_file )
		fftvol_tmp = sp_utilities.get_image(fftvol_eve_file)
		fftvol += fftvol_tmp
		fftvol_tmp = None

		weight = sp_utilities.get_image( weight_odd_file )
		weight_tmp = sp_utilities.get_image(weight_eve_file)
		weight += weight_tmp
		weight_tmp = None

		volall = recons_ctf_from_fftvol(nx, fftvol, weight, snr, symmetry, npad = npad)
		os.system( "rm -f " + fftvol_odd_file + " " + weight_odd_file )
		os.system( "rm -f " + fftvol_eve_file + " " + weight_eve_file )

		return volall,fscdat

	if nproc == 2:
		if myid == main_node_odd:
			fftvol = sp_utilities.get_image( fftvol_odd_file )
			weight = sp_utilities.get_image( weight_odd_file )
			volodd = recons_ctf_from_fftvol(nx, fftvol, weight, snr, symmetry, npad = npad)
			voleve = sp_utilities.recv_EMData(main_node_eve, tag_voleve, mpi_comm)
			
			if( not mask3D ):
				nx = volodd.get_xsize()
				ny = volodd.get_ysize()
				nz = volodd.get_zsize()
				mask3D = sp_utilities.model_circle(min(nx,ny,nz)//2 - 2, nx,ny,nz)
			fscdat = sp_statistics.fsc_mask( volodd, voleve, mask3D, rstep, fsc_curve)
			del  volodd, voleve, mask3D
		else:
			assert myid == main_node_eve
			fftvol = sp_utilities.get_image( fftvol_eve_file )
			weight = sp_utilities.get_image( weight_eve_file )
			voleve = recons_ctf_from_fftvol(nx, fftvol, weight, snr, symmetry, npad = npad)
			sp_utilities.send_EMData(voleve, main_node_odd, tag_voleve, mpi_comm)

		if myid == main_node_odd:
			fftvol = sp_utilities.get_image( fftvol_odd_file )
			fftvol_tmp = sp_utilities.recv_EMData( main_node_eve, tag_fftvol_eve, mpi_comm)
			fftvol += fftvol_tmp
			fftvol_tmp = None

			weight = sp_utilities.get_image( weight_odd_file )
			weight_tmp = sp_utilities.recv_EMData( main_node_eve, tag_weight_eve, mpi_comm)
			weight += weight_tmp
			weight_tmp = None

			volall = recons_ctf_from_fftvol(nx, fftvol, weight, snr, symmetry, npad = npad)
			os.system( "rm -f " + fftvol_odd_file + " " + weight_odd_file )

			return volall,fscdat
		else:
			assert myid == main_node_eve
			fftvol = sp_utilities.get_image( fftvol_eve_file )
			weight = sp_utilities.get_image( weight_eve_file )
			sp_utilities.send_EMData(fftvol, main_node_odd, tag_fftvol_eve, mpi_comm)
			sp_utilities.send_EMData(weight, main_node_odd, tag_weight_eve, mpi_comm)
			os.system( "rm -f " + fftvol_eve_file + " " + weight_eve_file )
			return sp_utilities.model_blank(nx,nx,nx), None

	# cases from all other number of processors situations
	if myid == main_node_odd:
		fftvol = sp_utilities.get_image( fftvol_odd_file )
		sp_utilities.send_EMData(fftvol, main_node_eve, tag_fftvol_odd, mpi_comm)

		if not(finfo is None):
			finfo.write("fftvol odd sent\n")
			finfo.flush()

		weight = sp_utilities.get_image( weight_odd_file )
		sp_utilities.send_EMData(weight, main_node_all, tag_weight_odd, mpi_comm)

		if not(finfo is None):
			finfo.write("weight odd sent\n")
			finfo.flush()

		volodd = recons_ctf_from_fftvol(nx, fftvol, weight, snr, symmetry, npad = npad)
		del fftvol, weight
		voleve = sp_utilities.recv_EMData(main_node_eve, tag_voleve, mpi_comm)

		if( not mask3D ):
			nx = volodd.get_xsize()
			ny = volodd.get_ysize()
			nz = volodd.get_zsize()
			mask3D = sp_utilities.model_circle(min(nx,ny,nz)//2 - 2, nx,ny,nz)

		fscdat = sp_statistics.fsc_mask(volodd, voleve, mask3D, rstep, fsc_curve)
		del  volodd, voleve, mask3D
		volall = sp_utilities.recv_EMData(main_node_all, tag_volall, mpi_comm)
		os.system( "rm -f " + fftvol_odd_file + " " + weight_odd_file )
		return volall, fscdat

	if myid == main_node_eve:
		ftmp = sp_utilities.recv_EMData(main_node_odd, tag_fftvol_odd, mpi_comm)
		fftvol = sp_utilities.get_image( fftvol_eve_file )
		EMAN2_cppwrap.Util.add_img( ftmp, fftvol )
		sp_utilities.send_EMData(ftmp, main_node_all, tag_fftvol_eve, mpi_comm)
		del ftmp

		weight = sp_utilities.get_image( weight_eve_file )
		sp_utilities.send_EMData(weight, main_node_all, tag_weight_eve, mpi_comm)

		voleve = recons_ctf_from_fftvol(nx, fftvol, weight, snr, symmetry, npad = npad)
		sp_utilities.send_EMData(voleve, main_node_odd, tag_voleve, mpi_comm)
		os.system( "rm -f " + fftvol_eve_file + " " + weight_eve_file );

		return sp_utilities.model_blank(nx,nx,nx), None


	if myid == main_node_all:
		fftvol = sp_utilities.recv_EMData(main_node_eve, tag_fftvol_eve, mpi_comm)
		if not(finfo is None):
			finfo.write( "fftvol odd received\n" )
			finfo.flush()

		weight = sp_utilities.recv_EMData(main_node_odd, tag_weight_odd, mpi_comm)
		weight_tmp = sp_utilities.recv_EMData(main_node_eve, tag_weight_eve, mpi_comm)
		EMAN2_cppwrap.Util.add_img( weight, weight_tmp )
		weight_tmp = None

		volall = recons_ctf_from_fftvol(nx, fftvol, weight, snr, symmetry, npad = npad)
		sp_utilities.send_EMData(volall, main_node_odd, tag_volall, mpi_comm)

		return sp_utilities.model_blank(nx,nx,nx),None

	return sp_utilities.model_blank(nx,nx,nx),None


def rec3D_MPI_with_getting_odd_even_volumes_from_files(fftvol_files, weight_files, reconstructed_vol_files,\
		nx, snr = 1.0, symmetry = "c1", mask3D = None, fsc_curve = None, \
		myid = 0, main_node = 0, rstep = 1.0, finfo=None, \
		npad = 2, mpi_comm=None):
	'''
	  This function is to be called within an MPI program to do a reconstruction on a dataset kept 
	  in the memory, computes reconstruction and through odd-even, in order to get the resolution
	'''
	
	fftvol_odd_file, fftvol_eve_file = fftvol_files 	
	weight_odd_file, weight_eve_file = weight_files
	reconstructed_odd_vol_files, reconstructed_eve_vol_files = reconstructed_vol_files
		
	pass#IMPORTIMPORTIMPORT import os
	pass#IMPORTIMPORTIMPORT from sp_statistics import fsc_mask
	pass#IMPORTIMPORTIMPORT from sp_utilities  import model_blank, model_circle, get_image, send_EMData, recv_EMData
	pass#IMPORTIMPORTIMPORT from mpi        import mpi_comm_size, MPI_COMM_WORLD
	
	if mpi_comm == None:
		mpi_comm = mpi.MPI_COMM_WORLD
	
	nproc = mpi.mpi_comm_size(mpi_comm)

	if nproc==1:
		assert main_node==0
		main_node_odd = main_node
		main_node_eve = main_node
		main_node_all = main_node
	elif nproc==2:
		main_node_odd = main_node
		main_node_eve = (main_node+1)%2
		main_node_all = main_node

		tag_voleve     = 1000
		tag_fftvol_eve = 1001
		tag_weight_eve = 1002
	else:
		#spread CPUs between different nodes to save memory
		main_node_odd = main_node
		main_node_eve = (int(main_node)+nproc-1)%int(nproc)
		main_node_all = (int(main_node)+nproc//2)%int(nproc)

		tag_voleve     = 1000
		tag_fftvol_eve = 1001
		tag_weight_eve = 1002

		tag_fftvol_odd = 1003
		tag_weight_odd = 1004
		tag_volall     = 1005

	if nproc == 1:
		fftvol = sp_utilities.get_image(fftvol_odd_file)
		weight = sp_utilities.get_image(weight_odd_file)
		# volodd = recons_ctf_from_fftvol(nx, fftvol, weight, snr, symmetry, npad = npad)
		volodd = sp_utilities.get_image(reconstructed_odd_vol_files)

		fftvol = sp_utilities.get_image(fftvol_eve_file)
		weight = sp_utilities.get_image(weight_eve_file)
		# voleve = recons_ctf_from_fftvol(nx, fftvol, weight, snr, symmetry, npad = npad)
		voleve = sp_utilities.get_image(reconstructed_eve_vol_files)

		if( not mask3D ):
			nx = volodd.get_xsize()
			ny = volodd.get_ysize()
			nz = volodd.get_zsize()
			mask3D = sp_utilities.model_circle(min(nx,ny,nz)//2 - 2, nx,ny,nz)
		fscdat = sp_statistics.fsc_mask( volodd, voleve, mask3D, rstep, fsc_curve)
		del  volodd, voleve, mask3d

		fftvol = sp_utilities.get_image( fftvol_odd_file )
		fftvol_tmp = sp_utilities.get_image(fftvol_eve_file)
		fftvol += fftvol_tmp
		fftvol_tmp = None

		weight = sp_utilities.get_image( weight_odd_file )
		weight_tmp = sp_utilities.get_image(weight_eve_file)
		weight += weight_tmp
		weight_tmp = None

		volall = recons_ctf_from_fftvol_using_nn4_ctfw(nx, fftvol, weight, snr, symmetry, npad = npad)
		os.system( "rm -f " + fftvol_odd_file + " " + weight_odd_file )
		os.system( "rm -f " + fftvol_eve_file + " " + weight_eve_file )

		return volall,fscdat

	if nproc == 2:
		if myid == main_node_odd:
			fftvol = sp_utilities.get_image( fftvol_odd_file )
			weight = sp_utilities.get_image( weight_odd_file )
			# volodd = recons_ctf_from_fftvol(nx, fftvol, weight, snr, symmetry, npad = npad)
			volodd = sp_utilities.get_image(reconstructed_odd_vol_files)
			voleve = sp_utilities.recv_EMData(main_node_eve, tag_voleve, mpi_comm)
			
			if( not mask3D ):
				nx = volodd.get_xsize()
				ny = volodd.get_ysize()
				nz = volodd.get_zsize()
				mask3D = sp_utilities.model_circle(min(nx,ny,nz)//2 - 2, nx,ny,nz)
			fscdat = sp_statistics.fsc_mask( volodd, voleve, mask3D, rstep, fsc_curve)
			del  volodd, voleve, mask3D
		else:
			assert myid == main_node_eve
			fftvol = sp_utilities.get_image( fftvol_eve_file )
			weight = sp_utilities.get_image( weight_eve_file )
			# voleve = recons_ctf_from_fftvol(nx, fftvol, weight, snr, symmetry, npad = npad)
			voleve = sp_utilities.get_image(reconstructed_eve_vol_files)
			sp_utilities.send_EMData(voleve, main_node_odd, tag_voleve, mpi_comm)

		if myid == main_node_odd:
			fftvol = sp_utilities.get_image( fftvol_odd_file )
			fftvol_tmp = sp_utilities.recv_EMData( main_node_eve, tag_fftvol_eve, mpi_comm)
			fftvol += fftvol_tmp
			fftvol_tmp = None

			weight = sp_utilities.get_image( weight_odd_file )
			weight_tmp = sp_utilities.recv_EMData( main_node_eve, tag_weight_eve, mpi_comm)
			weight += weight_tmp
			weight_tmp = None

			volall = recons_ctf_from_fftvol_using_nn4_ctfw(nx, fftvol, weight, snr, symmetry, npad = npad)
			os.system( "rm -f " + fftvol_odd_file + " " + weight_odd_file )

			return volall,fscdat
		else:
			assert myid == main_node_eve
			fftvol = sp_utilities.get_image( fftvol_eve_file )
			weight = sp_utilities.get_image( weight_eve_file )
			sp_utilities.send_EMData(fftvol, main_node_odd, tag_fftvol_eve, mpi_comm)
			sp_utilities.send_EMData(weight, main_node_odd, tag_weight_eve, mpi_comm)
			os.system( "rm -f " + fftvol_eve_file + " " + weight_eve_file )
			return sp_utilities.model_blank(nx,nx,nx), None

	# cases from all other number of processors situations
	if myid == main_node_odd:
		fftvol = sp_utilities.get_image( fftvol_odd_file )
		sp_utilities.send_EMData(fftvol, main_node_eve, tag_fftvol_odd, mpi_comm)

		if not(finfo is None):
			finfo.write("fftvol odd sent\n")
			finfo.flush()

		weight = sp_utilities.get_image( weight_odd_file )
		sp_utilities.send_EMData(weight, main_node_all, tag_weight_odd, mpi_comm)

		if not(finfo is None):
			finfo.write("weight odd sent\n")
			finfo.flush()

		# volodd = recons_ctf_from_fftvol(nx, fftvol, weight, snr, symmetry, npad = npad)
		volodd = sp_utilities.get_image(reconstructed_odd_vol_files)
		
		del fftvol, weight
		voleve = sp_utilities.recv_EMData(main_node_eve, tag_voleve, mpi_comm)

		if( not mask3D ):
			nx = volodd.get_xsize()
			ny = volodd.get_ysize()
			nz = volodd.get_zsize()
			mask3D = sp_utilities.model_circle(min(nx,ny,nz)//2 - 2, nx,ny,nz)

		fscdat = sp_statistics.fsc_mask(volodd, voleve, mask3D, rstep, fsc_curve)
		del  volodd, voleve, mask3D
		volall = sp_utilities.recv_EMData(main_node_all, tag_volall, mpi_comm)
		os.system( "rm -f " + fftvol_odd_file + " " + weight_odd_file )
		return volall, fscdat

	if myid == main_node_eve:
		ftmp = sp_utilities.recv_EMData(main_node_odd, tag_fftvol_odd, mpi_comm)
		fftvol = sp_utilities.get_image( fftvol_eve_file )
		EMAN2_cppwrap.Util.add_img( ftmp, fftvol )
		sp_utilities.send_EMData(ftmp, main_node_all, tag_fftvol_eve, mpi_comm)
		del ftmp

		weight = sp_utilities.get_image( weight_eve_file )
		sp_utilities.send_EMData(weight, main_node_all, tag_weight_eve, mpi_comm)

		# voleve = recons_ctf_from_fftvol(nx, fftvol, weight, snr, symmetry, npad = npad)
		voleve = sp_utilities.get_image(reconstructed_eve_vol_files)
		
		sp_utilities.send_EMData(voleve, main_node_odd, tag_voleve, mpi_comm)
		os.system( "rm -f " + fftvol_eve_file + " " + weight_eve_file )

		return sp_utilities.model_blank(nx,nx,nx), None


	if myid == main_node_all:
		fftvol = sp_utilities.recv_EMData(main_node_eve, tag_fftvol_eve, mpi_comm)
		if not(finfo is None):
			finfo.write( "fftvol odd received\n" )
			finfo.flush()

		weight = sp_utilities.recv_EMData(main_node_odd, tag_weight_odd, mpi_comm)
		weight_tmp = sp_utilities.recv_EMData(main_node_eve, tag_weight_eve, mpi_comm)
		EMAN2_cppwrap.Util.add_img( weight, weight_tmp )
		weight_tmp = None

		volall = recons_ctf_from_fftvol_using_nn4_ctfw(nx, fftvol, weight, snr, symmetry, npad = npad)
		sp_utilities.send_EMData(volall, main_node_odd, tag_volall, mpi_comm)

		return sp_utilities.model_blank(nx,nx,nx),None

	return sp_utilities.model_blank(nx,nx,nx),None


def rec3D_MPI_noCTF(data, symmetry = "c1", mask3D = None, fsc_curve = None, myid = 2, main_node = 0, \
		rstep = 1.0, odd_start=0, eve_start=1, finfo=None, index = -1, npad = 2, mpi_comm=None):
	'''
	  This function is to be called within an MPI program to do a reconstruction on a dataset kept in the memory 
	  Computes reconstruction and through odd-even, in order to get the resolution
	  if index > -1, projections should have attribute group set and only those whose group matches index will be used in the reconstruction
	    this is for multireference alignment
	'''
	pass#IMPORTIMPORTIMPORT import os
	pass#IMPORTIMPORTIMPORT from sp_statistics import fsc_mask
	pass#IMPORTIMPORTIMPORT from sp_utilities  import model_blank, get_image,send_EMData, recv_EMData
	pass#IMPORTIMPORTIMPORT from mpi        import mpi_comm_size, MPI_COMM_WORLD
	
	if mpi_comm == None:
		mpi_comm = mpi.MPI_COMM_WORLD
	
	nproc = mpi.mpi_comm_size(mpi_comm)

	if nproc==1:
		assert main_node==0
		main_node_odd = main_node
		main_node_eve = main_node
		main_node_all = main_node
	elif nproc==2:
		main_node_odd = main_node
		main_node_eve = (main_node+1)%2
		main_node_all = main_node

		tag_voleve     = 1000
		tag_fftvol_eve = 1001
		tag_weight_eve = 1002
	else:
		#spread CPUs between different nodes to save memory
		main_node_odd = main_node
		main_node_eve = (int(main_node)+nproc-1)%int(nproc)
		main_node_all = (int(main_node)+nproc//2)%int(nproc)

		tag_voleve     = 1000
		tag_fftvol_eve = 1001
		tag_weight_eve = 1002

		tag_fftvol_odd = 1003
		tag_weight_odd = 1004
		tag_volall     = 1005

	nx = data[0].get_xsize()

	fftvol_odd_file,weight_odd_file = prepare_recons(data, symmetry, myid, main_node_odd, odd_start, 2, index, finfo, npad, mpi_comm=mpi_comm)
	fftvol_eve_file,weight_eve_file = prepare_recons(data, symmetry, myid, main_node_eve, eve_start, 2, index, finfo, npad, mpi_comm=mpi_comm) 

	if nproc == 1:
		fftvol = sp_utilities.get_image( fftvol_odd_file )
		weight = sp_utilities.get_image( weight_odd_file )
		volodd = recons_from_fftvol(nx, fftvol, weight, symmetry, npad)

		fftvol = sp_utilities.get_image( fftvol_eve_file )
		weight = sp_utilities.get_image( weight_eve_file )
		voleve = recons_from_fftvol(nx, fftvol, weight, symmetry, npad)

		fscdat = sp_statistics.fsc_mask( volodd, voleve, mask3D, rstep, fsc_curve)
		del  volodd, voleve

		fftvol = sp_utilities.get_image( fftvol_odd_file )
		EMAN2_cppwrap.Util.add_img( fftvol, sp_utilities.get_image(fftvol_eve_file) )

		weight = sp_utilities.get_image( weight_odd_file )
		EMAN2_cppwrap.Util.add_img( weight, sp_utilities.get_image(weight_eve_file) )

		volall = recons_from_fftvol(nx, fftvol, weight, symmetry, npad)
		os.system( "rm -f " + fftvol_odd_file + " " + weight_odd_file );
		os.system( "rm -f " + fftvol_eve_file + " " + weight_eve_file );
		return volall,fscdat

	if nproc == 2:
		if myid == main_node_odd:
			fftvol = sp_utilities.get_image( fftvol_odd_file )
			weight = sp_utilities.get_image( weight_odd_file )
			volodd = recons_from_fftvol(nx, fftvol, weight, symmetry, npad)
			voleve = sp_utilities.recv_EMData(main_node_eve, tag_voleve, mpi_comm)
			fscdat = sp_statistics.fsc_mask( volodd, voleve, mask3D, rstep, fsc_curve)
			del  volodd, voleve
		else:
			assert myid == main_node_eve
			fftvol = sp_utilities.get_image( fftvol_eve_file )
			weight = sp_utilities.get_image( weight_eve_file )
			voleve = recons_from_fftvol(nx, fftvol, weight, symmetry, npad)
			sp_utilities.send_EMData(voleve, main_node_odd, tag_voleve, mpi_comm)

		if myid == main_node_odd:
			fftvol = sp_utilities.get_image( fftvol_odd_file )
			fftvol_tmp = sp_utilities.recv_EMData( main_node_eve, tag_fftvol_eve, mpi_comm)
			EMAN2_cppwrap.Util.add_img( fftvol, fftvol_tmp )
			fftvol_tmp = None

			weight = sp_utilities.get_image( weight_odd_file )
			weight_tmp = sp_utilities.recv_EMData( main_node_eve, tag_weight_eve, mpi_comm)
			EMAN2_cppwrap.Util.add_img( weight, weight_tmp )
			weight_tmp = None
			volall = recons_from_fftvol(nx, fftvol, weight, symmetry, npad)
			os.system( "rm -f " + fftvol_odd_file + " " + weight_odd_file );
			return volall,fscdat
		else:
			assert myid == main_node_eve
			fftvol = sp_utilities.get_image( fftvol_eve_file )
			sp_utilities.send_EMData(fftvol, main_node_odd, tag_fftvol_eve, mpi_comm)

			weight = sp_utilities.get_image( weight_eve_file )
			sp_utilities.send_EMData(weight, main_node_odd, tag_weight_eve, mpi_comm)
			os.system( "rm -f " + fftvol_eve_file + " " + weight_eve_file );
			return sp_utilities.model_blank(nx,nx,nx), None
	# cases from all other number of processors situations
	if myid == main_node_odd:
		fftvol = sp_utilities.get_image( fftvol_odd_file )
		sp_utilities.send_EMData(fftvol, main_node_eve, tag_fftvol_odd, mpi_comm)

		if not(finfo is None):
			finfo.write("fftvol odd sent\n")
			finfo.flush()

		weight = sp_utilities.get_image( weight_odd_file )
		sp_utilities.send_EMData(weight, main_node_all, tag_weight_odd, mpi_comm)

		if not(finfo is None):
			finfo.write("weight odd sent\n")
			finfo.flush()

		volodd = recons_from_fftvol(nx, fftvol, weight, symmetry, npad)
		del fftvol, weight
		voleve = sp_utilities.recv_EMData(main_node_eve, tag_voleve, mpi_comm)
		fscdat = sp_statistics.fsc_mask(volodd, voleve, mask3D, rstep, fsc_curve)
		del  volodd, voleve
		volall = sp_utilities.recv_EMData(main_node_all, tag_volall, mpi_comm)
		os.system( "rm -f " + fftvol_odd_file + " " + weight_odd_file );
		return volall,fscdat

	if myid == main_node_eve:
		ftmp = sp_utilities.recv_EMData(main_node_odd, tag_fftvol_odd, mpi_comm)
		fftvol = sp_utilities.get_image( fftvol_eve_file )
		EMAN2_cppwrap.Util.add_img( ftmp, fftvol )
		sp_utilities.send_EMData(ftmp, main_node_all, tag_fftvol_eve, mpi_comm)
		del ftmp

		weight = sp_utilities.get_image( weight_eve_file )
		sp_utilities.send_EMData(weight, main_node_all, tag_weight_eve, mpi_comm)

		voleve = recons_from_fftvol(nx, fftvol, weight, symmetry, npad)
		sp_utilities.send_EMData(voleve, main_node_odd, tag_voleve, mpi_comm)
		os.system( "rm -f " + fftvol_eve_file + " " + weight_eve_file );

		return sp_utilities.model_blank(nx,nx,nx), None


	if myid == main_node_all:
		fftvol = sp_utilities.recv_EMData(main_node_eve, tag_fftvol_eve, mpi_comm)
		if not(finfo is None):
			finfo.write( "fftvol odd received\n" )
			finfo.flush()

		weight = sp_utilities.recv_EMData(main_node_odd, tag_weight_odd, mpi_comm)
		weight_tmp = sp_utilities.recv_EMData(main_node_eve, tag_weight_eve, mpi_comm)
		EMAN2_cppwrap.Util.add_img( weight, weight_tmp )
		weight_tmp = None

		volall = recons_from_fftvol(nx, fftvol, weight, symmetry, npad)
		sp_utilities.send_EMData(volall, main_node_odd, tag_volall, mpi_comm)

		return sp_utilities.model_blank(nx,nx,nx),None


	return sp_utilities.model_blank(nx,nx,nx),None
	
def prepare_recons_ctf_two_chunks(nx,data,snr,symmetry,myid,main_node_half,chunk_ID,finfo=None,npad=2,mpi_comm=None,smearstep = 0.0):
	pass#IMPORTIMPORTIMPORT from random     import randint
	pass#IMPORTIMPORTIMPORT from sp_utilities  import reduce_EMData_to_root
	pass#IMPORTIMPORTIMPORT from mpi        import mpi_barrier, MPI_COMM_WORLD
	pass#IMPORTIMPORTIMPORT from EMAN2 import Reconstructors

	if mpi_comm == None:
		mpi_comm = mpi.MPI_COMM_WORLD

	fftvol_half = EMAN2_cppwrap.EMData()

	if( smearstep > 0.0 ):
		#if myid == 0:  print "  Setting smear in prepare_recons_ctf"
		ns = 1
		smear = []
		for j in range(-ns,ns+1):
			if( j != 0):
				for i in range(-ns,ns+1):
					for k in range(-ns,ns+1):
						smear += [i*smearstep,j*smearstep,k*smearstep,1.0]
		# Deal with theta = 0.0 cases
		prj = []
		for i in range(-ns,ns+1):
			for k in range(-ns,ns+1):
				prj.append(i+k)
		for i in range(-2*ns,2*ns+1,1):
			 smear += [i*smearstep,0.0,0.0,float(prj.count(i))]
		#if myid == 0:  print "  Smear  ",smear
		fftvol_half.set_attr("smear", smear)

	weight_half = EMAN2_cppwrap.EMData()
	half_params = {"size":nx, "npad":npad, "snr":snr, "sign":1, "symmetry":symmetry, "fftvol":fftvol_half, "weight":weight_half}
	half = EMAN2_cppwrap.Reconstructors.get( "nn4_ctf", half_params )
	half.setup()
	for i in range(len(data)):
		if data[i].get_attr("chunk_id") == chunk_ID:
			xform_proj = data[i].get_attr( "xform.projection" )
			half.insert_slice(data[i], xform_proj )
	if not(finfo is None):
		finfo.write( "begin reduce half\n" )
		finfo.flush()

	sp_utilities.reduce_EMData_to_root(fftvol_half, myid, main_node_half, mpi_comm)
	sp_utilities.reduce_EMData_to_root(weight_half, myid, main_node_half, mpi_comm)

	if not(finfo is None):
		finfo.write( "after reduce half\n" )
		finfo.flush()

	if myid == main_node_half:
		tmpid = random.randint(0, 1000000) 
		fftvol_half_file = ("fftvol_half%d.hdf" % tmpid)
		weight_half_file = ("weight_half%d.hdf" % tmpid)
		fftvol_half.write_image(fftvol_half_file)
		weight_half.write_image(weight_half_file)
	mpi.mpi_barrier(mpi_comm)

	fftvol_half = None
	weight_half = None

	if myid == main_node_half:
		return fftvol_half_file, weight_half_file

	return None,None
	
def rec3D_two_chunks_MPI(data, snr = 1.0, symmetry = "c1", mask3D = None, fsc_curve = None, \
		myid = 0, main_node = 0, rstep = 1.0, finfo=None, \
		index=-1, npad = 2, mpi_comm=None, smearstep = 0.0):
	'''
	  This function is to be called within an MPI program to do a reconstruction on a dataset kept 
	  in the memory, computes reconstruction and through odd-even, in order to get the resolution
	'''
	pass#IMPORTIMPORTIMPORT import os
	pass#IMPORTIMPORTIMPORT from sp_statistics import fsc_mask
	pass#IMPORTIMPORTIMPORT from sp_utilities  import model_blank, model_circle, get_image, send_EMData, recv_EMData
	pass#IMPORTIMPORTIMPORT from mpi        import mpi_comm_size, MPI_COMM_WORLD
	
	if mpi_comm == None:
		mpi_comm = mpi.MPI_COMM_WORLD
	
	nproc = mpi.mpi_comm_size(mpi_comm)

	if nproc==1:
		assert main_node==0
		main_node_odd = main_node
		main_node_eve = main_node
		main_node_all = main_node
	elif nproc==2:
		main_node_odd = main_node
		main_node_eve = (main_node+1)%2
		main_node_all = main_node

		tag_voleve     = 1000
		tag_fftvol_eve = 1001
		tag_weight_eve = 1002
	else:
		#spread CPUs between different nodes to save memory
		main_node_odd = main_node
		main_node_eve = (int(main_node)+nproc-1)%int(nproc)
		main_node_all = (int(main_node)+nproc//2)%int(nproc)

		tag_voleve     = 1000
		tag_fftvol_eve = 1001
		tag_weight_eve = 1002

		tag_fftvol_odd = 1003
		tag_weight_odd = 1004
		tag_volall     = 1005


	if index != -1 :
		grpdata = []
		for i in range(len(data)):
			if data[i].get_attr('group') == index:
				grpdata.append(data[i])
		imgdata = grpdata
	else:
		imgdata = data

	nx = get_image_size(imgdata, myid)
	if nx == 0:
		sp_global_def.ERROR("Warning: no images were given for reconstruction, this usually means there is an empty group, returning empty volume", "rec3D", 0)
		return sp_utilities.model_blank( 2, 2, 2 ), None

	fftvol_odd_file,weight_odd_file = prepare_recons_ctf_two_chunks(nx, imgdata, snr, symmetry, myid, main_node_odd, 0, finfo, npad, mpi_comm=mpi_comm, smearstep = smearstep)
	fftvol_eve_file,weight_eve_file = prepare_recons_ctf_two_chunks(nx, imgdata, snr, symmetry, myid, main_node_eve, 1, finfo, npad, mpi_comm=mpi_comm, smearstep = smearstep)
	del imgdata

	if nproc == 1:
		fftvol = sp_utilities.get_image(fftvol_odd_file)
		weight = sp_utilities.get_image(weight_odd_file)
		volodd = recons_ctf_from_fftvol(nx, fftvol, weight, snr, symmetry, npad = npad)

		fftvol = sp_utilities.get_image(fftvol_eve_file)
		weight = sp_utilities.get_image(weight_eve_file)
		voleve = recons_ctf_from_fftvol(nx, fftvol, weight, snr, symmetry, npad = npad)

		if( not mask3D ):
			nx = volodd.get_xsize()
			ny = volodd.get_ysize()
			nz = volodd.get_zsize()
			mask3D = sp_utilities.model_circle(min(nx,ny,nz)//2 - 2, nx,ny,nz)
		fscdat = sp_statistics.fsc_mask( volodd, voleve, mask3D, rstep, fsc_curve)
		del  volodd, voleve, mask3d

		fftvol = sp_utilities.get_image( fftvol_odd_file )
		fftvol_tmp = sp_utilities.get_image(fftvol_eve_file)
		fftvol += fftvol_tmp
		fftvol_tmp = None

		weight = sp_utilities.get_image( weight_odd_file )
		weight_tmp = sp_utilities.get_image(weight_eve_file)
		weight += weight_tmp
		weight_tmp = None

		volall = recons_ctf_from_fftvol(nx, fftvol, weight, snr, symmetry, npad = npad)
		os.system( "rm -f " + fftvol_odd_file + " " + weight_odd_file )
		os.system( "rm -f " + fftvol_eve_file + " " + weight_eve_file )

		return volall,fscdat

	if nproc == 2:
		if myid == main_node_odd:
			fftvol = sp_utilities.get_image( fftvol_odd_file )
			weight = sp_utilities.get_image( weight_odd_file )
			volodd = recons_ctf_from_fftvol(nx, fftvol, weight, snr, symmetry, npad = npad)
			voleve = sp_utilities.recv_EMData(main_node_eve, tag_voleve, mpi_comm)
			
			if( not mask3D ):
				nx = volodd.get_xsize()
				ny = volodd.get_ysize()
				nz = volodd.get_zsize()
				mask3D = sp_utilities.model_circle(min(nx,ny,nz)//2 - 2, nx,ny,nz)
			fscdat = sp_statistics.fsc_mask( volodd, voleve, mask3D, rstep, fsc_curve)
			del  volodd, voleve, mask3D
		else:
			assert myid == main_node_eve
			fftvol = sp_utilities.get_image( fftvol_eve_file )
			weight = sp_utilities.get_image( weight_eve_file )
			voleve = recons_ctf_from_fftvol(nx, fftvol, weight, snr, symmetry, npad = npad)
			sp_utilities.send_EMData(voleve, main_node_odd, tag_voleve, mpi_comm)

		if myid == main_node_odd:
			fftvol = sp_utilities.get_image( fftvol_odd_file )
			fftvol_tmp = sp_utilities.recv_EMData( main_node_eve, tag_fftvol_eve, mpi_comm)
			fftvol += fftvol_tmp
			fftvol_tmp = None

			weight = sp_utilities.get_image( weight_odd_file )
			weight_tmp = sp_utilities.recv_EMData( main_node_eve, tag_weight_eve, mpi_comm)
			weight += weight_tmp
			weight_tmp = None

			volall = recons_ctf_from_fftvol(nx, fftvol, weight, snr, symmetry, npad = npad)
			os.system( "rm -f " + fftvol_odd_file + " " + weight_odd_file )

			return volall,fscdat
		else:
			assert myid == main_node_eve
			fftvol = sp_utilities.get_image( fftvol_eve_file )
			weight = sp_utilities.get_image( weight_eve_file )
			sp_utilities.send_EMData(fftvol, main_node_odd, tag_fftvol_eve, mpi_comm)
			sp_utilities.send_EMData(weight, main_node_odd, tag_weight_eve, mpi_comm)
			os.system( "rm -f " + fftvol_eve_file + " " + weight_eve_file )
			return sp_utilities.model_blank(nx,nx,nx), None

	# cases from all other number of processors situations
	if myid == main_node_odd:
		fftvol = sp_utilities.get_image( fftvol_odd_file )
		sp_utilities.send_EMData(fftvol, main_node_eve, tag_fftvol_odd, mpi_comm)

		if not(finfo is None):
			finfo.write("fftvol odd sent\n")
			finfo.flush()

		weight = sp_utilities.get_image( weight_odd_file )
		sp_utilities.send_EMData(weight, main_node_all, tag_weight_odd, mpi_comm)

		if not(finfo is None):
			finfo.write("weight odd sent\n")
			finfo.flush()

		volodd = recons_ctf_from_fftvol(nx, fftvol, weight, snr, symmetry, npad = npad)
		del fftvol, weight
		voleve = sp_utilities.recv_EMData(main_node_eve, tag_voleve, mpi_comm)

		if( not mask3D ):
			nx = volodd.get_xsize()
			ny = volodd.get_ysize()
			nz = volodd.get_zsize()
			mask3D = sp_utilities.model_circle(min(nx,ny,nz)//2 - 2, nx,ny,nz)

		fscdat = sp_statistics.fsc_mask(volodd, voleve, mask3D, rstep, fsc_curve)
		del  volodd, voleve, mask3D
		volall = sp_utilities.recv_EMData(main_node_all, tag_volall, mpi_comm)
		os.system( "rm -f " + fftvol_odd_file + " " + weight_odd_file )
		return volall, fscdat

	if myid == main_node_eve:
		ftmp = sp_utilities.recv_EMData(main_node_odd, tag_fftvol_odd, mpi_comm)
		fftvol = sp_utilities.get_image( fftvol_eve_file )
		EMAN2_cppwrap.Util.add_img( ftmp, fftvol )
		sp_utilities.send_EMData(ftmp, main_node_all, tag_fftvol_eve, mpi_comm)
		del ftmp

		weight = sp_utilities.get_image( weight_eve_file )
		sp_utilities.send_EMData(weight, main_node_all, tag_weight_eve, mpi_comm)

		voleve = recons_ctf_from_fftvol(nx, fftvol, weight, snr, symmetry, npad = npad)
		sp_utilities.send_EMData(voleve, main_node_odd, tag_voleve, mpi_comm)
		os.system( "rm -f " + fftvol_eve_file + " " + weight_eve_file );

		return sp_utilities.model_blank(nx,nx,nx), None


	if myid == main_node_all:
		fftvol = sp_utilities.recv_EMData(main_node_eve, tag_fftvol_eve, mpi_comm)
		if not(finfo is None):
			finfo.write( "fftvol odd received\n" )
			finfo.flush()

		weight = sp_utilities.recv_EMData(main_node_odd, tag_weight_odd, mpi_comm)
		weight_tmp = sp_utilities.recv_EMData(main_node_eve, tag_weight_eve, mpi_comm)
		EMAN2_cppwrap.Util.add_img( weight, weight_tmp )
		weight_tmp = None

		volall = recons_ctf_from_fftvol(nx, fftvol, weight, snr, symmetry, npad = npad)
		sp_utilities.send_EMData(volall, main_node_odd, tag_volall, mpi_comm)

		return sp_utilities.model_blank(nx,nx,nx),None

	return sp_utilities.model_blank(nx,nx,nx),None





