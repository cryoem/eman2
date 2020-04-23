

















































































































































































































































































'''0
# secondrun
def secondrunrecons3d_4nnw_MPI(myid, prjlist, prevol, symmetry="c1", finfo=None, npad=2, mpi_comm=None):
	from sp_utilities     import reduce_EMData_to_root, pad, get_params_proj
	from EMAN2         import Reconstructors
	from sp_utilities     import iterImagesList, model_blank, model_circle, reshape_1d, read_text_file
	from sp_fundamentals  import fft, rops
	from mpi           import MPI_COMM_WORLD
	import types

	if mpi_comm == None:
		mpi_comm = MPI_COMM_WORLD

	if type(prjlist) == types.ListType:
		prjlist = iterImagesList(prjlist)

	if not prjlist.goToNext():
		ERROR("empty input list","recons3d_4nn_MPI",1)

	imgsize = prjlist.image().get_xsize()
	bigsize = imgsize*npad
	bnx     = bigsize//2+1

	prjlist.goToPrev()

	fftvol = EMData()
	weight = EMData()
	print "   NEW"
	#t = read_text_file('fromrun8model.txt',4)
	#  GET FSC
	t = read_text_file('data_model.txt',4)
	from math import sqrt
	for i in xrange(len(t)):
		t[i] = max(t[i],0.0)
		#  This is what is used to get the SSNR
		t[i] = sqrt(2*t[i]/(1.0+t[i]))
	t = reshape_1d(t,len(t),npad*len(t))
	refvol = model_blank(2*bnx,1,1,0.5)
	for i in xrange(len(t)):  refvol.set_value_at(i,t[i])
	"""
	from math import tanh,pi
	fl = 0.15
	aa = 0.15
	for i in xrange(bnx):
		r = float(i)/bigsize
		refvol.set_value_at(i, 0.5*( tanh(pi*(r+fl)/2.0/fl/aa) - tanh(pi*(r-fl)/2.0/2.0/fl/aa) ) )
		print "  FILTER  ",i,refvol.get_value_at(i)
	"""
	#print " DONE refvol"
	params = {"size":imgsize, "npad":npad, "symmetry":symmetry, "fftvol":fftvol, "refvol":refvol, "weight":weight, "weighting":0, "snr":1.0}
	r = Reconstructors.get( "nn4_ctfw", params )
	r.setup()

	from sp_projection import prep_vol, prgs
	from sp_filter import filt_ctf
	#volft,kb = prep_vol(prevol)

	#mask2d = model_circle(imgsize//2-2, imgsize,imgsize)
	#maskbi = model_circle(imgsize//2-2, bigsize,bigsize)
	#  noise model of 2D data.
	models = [None]
	for ml in xrange(len(models)):
		temp = read_text_file('sigma2.txt',2)
		temp = reshape_1d(temp, len(temp), 2*len(temp))
		models[ml] = model_blank(len(temp)+10)
		for lm in xrange(len(temp)):  models[ml].set_value_at(lm,1.0/(temp[lm]*4*imgsize**2/npad))
		#from sys import exit
		#print "metadata/model-%04d.txt"%groupkeys[1][ml]
		#for lm in xrange(len(temp)):  print  lm,models[ml].get_value_at(lm)
		#exit()


	if not (finfo is None): nimg = 0
	ll = 0
	while prjlist.goToNext():
		prj = prjlist.image()

		# horatio active_refactoring Jy51i1EwmLD4tWZ9_00000_1
		# active = prj.get_attr_default('active', 1)
		# if(active == 1):
		if ll%100 == 0:  print "  moved  ",ll
		ll +=1
		prj.set_attr("sigmasq2", models[0])
		#if ll == 0:
		#	write_text_file([range(bigsize),[pqdif[i] for i in xrange(bigsize)] ],"pqdif.txt")
		#	ll+=1
		insert_slices(r, prj)
		if( not (finfo is None) ):
			nimg += 1
			info.write("Image %4d inserted.\n" %(nimg) )
			info.flush()

		if ll%100 == 0:  print "  moved  ",ll
		ll +=1
		prj.set_attr("sigmasq2", models[0])
		#if ll == 0:
		#	write_text_file([range(bigsize),[pqdif[i] for i in xrange(bigsize)] ],"pqdif.txt")
		#	ll+=1
		insert_slices(r, prj)
		if( not (finfo is None) ):
			nimg += 1
			info.write("Image %4d inserted.\n" %(nimg) )
			info.flush()

	if not (finfo is None):
		info.write( "Begin reducing ...\n" )
		info.flush()

	#del qdif, pqdif, mask2d, maskbi

	reduce_EMData_to_root(fftvol, myid, comm=mpi_comm)
	reduce_EMData_to_root(weight, myid, comm=mpi_comm)

	if myid == 0:
		print  "  STARTING FINISH"
		dummy = r.finish(True)
	else:
		from sp_utilities import model_blank
		fftvol = model_blank(imgsize, imgsize, imgsize)
	return fftvol
'''
'''1
#chc5
def recons3d_4nnw_MPI(myid, prjlist, prevol, symmetry="c1", finfo=None, npad=2, mpi_comm=None):
	from sp_utilities     import reduce_EMData_to_root, pad, get_params_proj
	from EMAN2         import Reconstructors
	from sp_utilities     import iterImagesList, model_blank, model_circle, reshape_1d, read_text_file
	from sp_fundamentals  import fft, rops
	from mpi           import MPI_COMM_WORLD
	import types

	if mpi_comm == None:
		mpi_comm = MPI_COMM_WORLD

	if type(prjlist) == types.ListType:
		prjlist = iterImagesList(prjlist)

	if not prjlist.goToNext():
		ERROR("empty input list","recons3d_4nn_MPI",1)

	imgsize = prjlist.image().get_xsize()
	bigsize = imgsize*npad
	bnx     = bigsize//2+1

	prjlist.goToPrev()

	fftvol = EMData()
	weight = EMData()
	"""
	if myid == 0:
		model_blank(bnx, bigsize, bigsize)
		temp = fft(pad(prevol,bigsize,bigsize,bigsize,0.0))
		temp.set_attr("is_complex",0)
		st = 0.5/(bigsize*bigsize)
		for kk in xrange(bigsize):
			for jj in xrange(bigsize):
				for ii in xrange(0,bnx,2):
					#print ii,jj,kk,temp.get_value_at(ii,jj,kk), temp.get_value_at(ii,jj,kk+1)
					refvol.set_value_at_fast(ii//2,jj,kk,st/((temp.get_value_at(ii,jj,kk))**2+(temp.get_value_at(ii+1,jj,kk))**2) )
					#refvol.set_value_at_fast(ii//2,jj,kk,1.0 )
		refvol.set_value_at_fast(0,0,0,0.0)
		del temp

		st = rops(pad(prevol,bigsize,bigsize,bigsize,0.0))*(bigsize**6)/4.
		from sp_utilities import info
		from sp_utilities import write_text_file
		#zizi = [st.get_value_at(i) for i in xrange(st.get_xsize())]
		#for i in xrange(st.get_xsize()):  st.set_value_at(i,1.0)#/st.get_value_at(i))
		#info(st,None,"refvol")
		refvol = model_blank(bigsize,1,1,1.0)
		for i in xrange(st.get_xsize()):  refvol.set_value_at(i,1.0/(211*st.get_value_at(i)))
	else:  refvol = EMData()
	"""
	print "   NEW"
	#t = read_text_file('fromrun8model.txt',4)
	t = read_text_file('../for-pawel/fsc-relion.txt',1)
	from math import sqrt
	for i in xrange(len(t)):
		t[i] = max(t[i],0.0)
		t[i] = sqrt(2*t[i]/(1.0+t[i]))
	t = reshape_1d(t,len(t),npad*len(t))
	refvol = model_blank(2*bnx,1,1,0.5)
	#for i in xrange(len(t)):  refvol.set_value_at(i,t[i])
	"""
	from math import tanh,pi
	fl = 0.15
	aa = 0.15
	for i in xrange(bnx):
		r = float(i)/bigsize
		refvol.set_value_at(i, 0.5*( tanh(pi*(r+fl)/2.0/fl/aa) - tanh(pi*(r-fl)/2.0/2.0/fl/aa) ) )
		print "  FILTER  ",i,refvol.get_value_at(i)
	"""
	#print " DONE refvol"
	params = {"size":imgsize, "npad":npad, "symmetry":symmetry, "fftvol":fftvol, "refvol":refvol, "weight":weight, "weighting":0, "snr":1.0}
	r = Reconstructors.get( "nn4_ctfw", params )
	r.setup()

	from sp_projection import prep_vol, prgs
	from sp_filter import filt_ctf
	#volft,kb = prep_vol(prevol)

	#mask2d = model_circle(imgsize//2-2, imgsize,imgsize)
	#maskbi = model_circle(imgsize//2-2, bigsize,bigsize)

	groupkeys = read_text_file("groupkeys.txt",-1)
	for ml in xrange(3):  groupkeys[ml] = map(int, groupkeys[ml])
	models = [None]*len(groupkeys[0])
	for ml in xrange(len(models)):
		temp = read_text_file("metadata/model-%04d.txt"%groupkeys[1][ml],-1)
		temp = reshape_1d(temp[2], len(temp[0]), 2*len(temp[0]))
		models[ml] = model_blank(len(temp)+10)
		for lm in xrange(len(temp)):  models[ml].set_value_at(lm,1.0/(temp[lm]*4*imgsize**2/npad))
		#from sys import exit
		#print "metadata/model-%04d.txt"%groupkeys[1][ml]
		#for lm in xrange(len(temp)):  print  lm,models[ml].get_value_at(lm)
		#exit()


	if not (finfo is None): nimg = 0
	ll = 0
	while prjlist.goToNext():
		prj = prjlist.image()
		if ll%100 == 0:  print "  moved  ",ll
		ll +=1
		ml = prj.get_attr('groupindex')#int(prj.get_attr('data_path')[4:8])
		prj.set_attr("sigmasq2", models[groupkeys[1].index(ml)])
		insert_slices(r, prj)
		if( not (finfo is None) ):
			nimg += 1
			info.write("Image %4d inserted.\n" %(nimg) )
			info.flush()

	if not (finfo is None):
		info.write( "Begin reducing ...\n" )
		info.flush()

	#del qdif, pqdif, mask2d, maskbi

	reduce_EMData_to_root(fftvol, myid, comm=mpi_comm)
	reduce_EMData_to_root(weight, myid, comm=mpi_comm)

	if myid == 0:
		print  "  STARTING FINISH"
		dummy = r.finish(True)
	else:
		from sp_utilities import model_blank
		fftvol = model_blank(imgsize, imgsize, imgsize)
	return fftvol
'''































































































































































'''2
def recons3d_4nnf_MPI(myid, list_of_prjlist, bckgdata, snr = 1.0, sign=1, symmetry="c1", finfo=None, npad=2, mpi_comm=None, smearstep = 0.0, main_node = 0):
	"""
		recons3d_4nn_ctf - calculate CTF-corrected 3-D reconstruction from a set of projections using three Eulerian angles, two shifts, and CTF settings for each projeciton image
		Input
			list_of_prjlist: list of lists of projections to be included in the reconstruction
			bckgdata = [get_im("tsd.hdf"),read_text_file("data_stamp.txt")]
			snr: Signal-to-Noise Ratio of the data 
			sign: sign of the CTF
			symmetry: point-group symmetry to be enforced, each projection will enter the reconstruction in all symmetry-related directions.
	"""
	from sp_utilities  import reduce_EMData_to_root, get_im, send_string_to_all
	from sp_statistics import fsc
	from sp_utilities  import get_image, send_EMData, recv_EMData

	from EMAN2      import Reconstructors
	from sp_utilities  import model_blank, cmdexecute
	from mpi        import MPI_COMM_WORLD, mpi_barrier, mpi_comm_size 
	import types
	from sp_statistics import fsc
	import datetime
	
	if mpi_comm == None:
		mpi_comm = MPI_COMM_WORLD
	main_node = 0
	imgsize = list_of_prjlist[0][0].get_xsize()
	
	if( smearstep > 0.0 ):
		#if myid == 0:  print "  Setting smear in prepare_recons_ctf"
		ns = 1
		smear = []
		for j in xrange(-ns,ns+1):
			if( j != 0):
				for i in xrange(-ns,ns+1):
					for k in xrange(-ns,ns+1):
						smear += [i*smearstep,j*smearstep,k*smearstep,1.0]
		# Deal with theta = 0.0 cases
		prj = []
		for i in xrange(-ns,ns+1):
			for k in xrange(-ns,ns+1):
				prj.append(i+k)
		for i in xrange(-2*ns,2*ns+1,1):
			smear += [i*smearstep,0.0,0.0,float(prj.count(i))]
		#if myid == 0:  print "  Smear  ",smear

	#from utilities import model_blank, get_im, read_text_file
	#bckgdata = [get_im("tsd.hdf"),read_text_file("data_stamp.txt")]

	nnx = bckgdata[0].get_xsize()
	nny = bckgdata[0].get_ysize()
	bckgnoise = []
	for i in xrange(nny):
		prj = model_blank(nnx)
		for k in xrange(nnx):  prj[k] = bckgdata[0].get_value_at(k,i)
		bckgnoise.append(prj)

	datastamp = bckgdata[1]

	#  Do the FSC shtick.
	bnx     = imgsize*npad//2+1
	refvol = model_blank(bnx)  # fill fsc with zeroes so the first reconstruction is done using simple Wiener filter.
	refvol.set_attr("fudge", 1.0)

	reconstructed_vols_list = []
	fftvol_file =[]
	weight_file = []
	

	volall_files = []
	volall = []
	fscdat_list = []
	no_of_splits = 2
	for iset in xrange(2):
		reconstructed_vol_files = []
		for iter_no_of_splits in range(no_of_splits):
			if not (finfo is None): nimg = 0
	
			fftvol = EMData()
			weight = EMData()
			if( smearstep > 0.0 ):  fftvol.set_attr("smear", smear)
		
			params = {"size":imgsize, "npad":npad, "snr":snr, "sign":sign, "symmetry":symmetry, "refvol":refvol, "fftvol":fftvol, "weight":weight}
			r = Reconstructors.get( "nn4_ctfw", params )
			r.setup()

			for image in list_of_prjlist[iset][slice(iter_no_of_splits,len(list_of_prjlist[iset]), no_of_splits)]:	
			# for image in list_of_prjlist[iset]:
				try:
					# lll = iset/0
					#raise ValueError('A very specific thing happened')
					stmp = image.get_attr("ptcl_source_image")
				except:
					try:
						stmp = image.get_attr("ctf")
						stmp = round(stmp.defocus,4)
					except:
						ERROR("Either ptcl_source_image or ctf has to be present in the header.","recons3d_4nnw_MPI    %f"%stmp,1, myid)
				try:
					indx = datastamp.index(stmp)
				except:
					ERROR("Problem with indexing ptcl_source_image.","recons3d_4nnf_MPI   %s"%stmp,1, myid)

			if not (finfo is None): 
				finfo.write( "begin reduce\n" )
				finfo.flush()
		
			reduce_EMData_to_root(fftvol, myid, main_node, comm=mpi_comm)
			reduce_EMData_to_root(weight, myid, main_node, comm=mpi_comm)
			
			if not (finfo is None): 
				finfo.write( "after reduce\n" )
				finfo.flush()
	
			if myid == 0:
				tmpid = datetime.datetime.now().strftime('%Y-%m-%d--%I-%M-%f')[:-3]
				fftvol_file.append("fftvol__%s__idx%d_split%d.hdf"%(tmpid, iset, iter_no_of_splits))
				weight_file.append("weight__%s__idx%d_split%d.hdf"%(tmpid, iset, iter_no_of_splits))
				fftvol.write_image(fftvol_file[-1])
				weight.write_image(weight_file[-1])
				del weight
	
				dummy = r.finish(True)
				### finish returns fftvol, will rename to 

				if(iter_no_of_splits == 0):
					reconstructed_vol_files.append("rvol__%s__idx%d_split%d.hdf"%(tmpid, iset, iter_no_of_splits))
					fftvol.write_image(reconstructed_vol_files[-1])
					del fftvol
				else:
					reconstructed_vol_files.append(fftvol)
	
		mpi_barrier(mpi_comm)
			
		### at this point we calculate V[iset]
		### fftvol_file = [FE, FO]
		### weight_file = [WE, WO]
		### reconstructed_vol_files = [VE(disk), VO(memory)]
	
		########################################	
		### code not tested 2015-12-21--13-01-49 
		########################################
		
		nproc = mpi_comm_size(mpi_comm)
		assert nproc > 2
		main_node_odd = main_node
		main_node_eve = (int(main_node)+nproc-1)%int(nproc)
		main_node_all = (int(main_node)+nproc//2)%int(nproc)

		tag_voleve     = 1000
		tag_fftvol_eve = 1001
		tag_weight_eve = 1002

		tag_fftvol_odd = 1003
		tag_weight_odd = 1004
		tag_volall     = 1005
		
		fftvol_file = eval(send_string_to_all(str(fftvol_file), source_node = main_node))
		weight_file = eval(send_string_to_all(str(weight_file), source_node = main_node))
		
		fftvol_odd_file, fftvol_eve_file = fftvol_file[-2:]
		weight_odd_file, weight_eve_file = weight_file[-2:]
		
		if myid == main_node_odd:
			fftvol = get_image( fftvol_odd_file )
			send_EMData(fftvol, main_node_eve, tag_fftvol_odd, mpi_comm)

			weight = get_image( weight_odd_file )
			send_EMData(weight, main_node_all, tag_weight_odd, mpi_comm)
			volodd = reconstructed_vol_files[1]
			del fftvol, weight
			voleve = get_im(reconstructed_vol_files[0])
			fscdat_list.append(fsc(volodd, voleve, 1.0)[1])
			del  volodd, voleve

			volall = recv_EMData(main_node_all, tag_volall, mpi_comm)
			if iset == 0:
				volall_files.append("volall__%s__idx%d.hdf"%(tmpid, iset))
				volall.write_image(volall_files[-1])
			else:
				volall_files.append(volall)

		if myid == main_node_eve:
			ftmp = recv_EMData(main_node_odd, tag_fftvol_odd, mpi_comm)
			fftvol = get_image( fftvol_eve_file )
			Util.add_img( ftmp, fftvol )

			send_EMData(ftmp, main_node_all, tag_fftvol_eve, mpi_comm)
			del ftmp
			weight = get_image( weight_eve_file )
			send_EMData(weight, main_node_all, tag_weight_eve, mpi_comm)

		if myid == main_node_all:


			fftvol = recv_EMData(main_node_eve, tag_fftvol_eve, mpi_comm)
			weight = recv_EMData(main_node_odd, tag_weight_odd, mpi_comm)
			weight_tmp = recv_EMData(main_node_eve, tag_weight_eve, mpi_comm)
			Util.add_img( weight, weight_tmp )

			weight_tmp = None
			volall = recons_ctf_from_fftvol_using_nn4_ctfw(imgsize, fftvol, weight, snr, symmetry, npad = npad)
			send_EMData(volall, main_node_odd, tag_volall, mpi_comm)

	### at this point we calculate fourier_shell_correlation for V[0], V[1]

	if myid == 0:
		fourier_shell_correlation = fsc(get_im(volall_files[0]), volall_files[1], 1.0)[1]
		fourier_shell_correlation[0] = 1.0

		from math import sqrt
		from sp_utilities import reshape_1d
		t = [0.0]*len(fourier_shell_correlation)
		t = reshape_1d(t,len(t),npad*len(t))
		for i in xrange(len(t):
			t[i] = min(max(t[i], 0.0), 0.999)


		ovol = []
		for idx in range(2):
			fftvol = get_im(fftvol_file[idx])
			weight = get_im(weight_file[idx])
			refvol = model_blank(bnx,1,1,0.0)
			for i in xrange(len(t)):  
				refvol.set_value_at(i, t[i])
			refvol.set_attr("fudge", 1.0)
			
			params = {"size":imgsize, "npad":npad, "snr":snr, "sign":sign, "symmetry":symmetry, "refvol":refvol, "fftvol":fftvol, "weight":weight}
			r = Reconstructors.get("nn4_ctfw", params)
			r.setup()
			
			dummy = r.finish(True)
			ovol.append(fftvol)


		cmd = "{} {} {} {} {} {}".format("rm -f", fftvol_file[0], fftvol_file[1], weight_file[0], weight_file[1], results_list[0] )
		junk = cmdexecute(cmd)

	mpi_barrier(mpi_comm)
	if myid == 0:
		return ovol[0], ovol[1], fourier_shell_correlation, fscdat_list[0], fscdat_list[1]
	else:
		return None, None, None, None, None

'''
'''3
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
	from sp_utilities  import reduce_EMData_to_root, random_string, get_im
	from EMAN2      import Reconstructors
	from sp_utilities  import model_blank
	from mpi        import MPI_COMM_WORLD, mpi_barrier
	import types
	from sp_statistics import fsc
	import datetime
	
	if mpi_comm == None:
		mpi_comm = MPI_COMM_WORLD
	main_node = 0
	imgsize = list_of_prjlist[0][0].get_xsize()
	
	if( smearstep > 0.0 ):
		#if myid == 0:  print "  Setting smear in prepare_recons_ctf"
		ns = 1
		smear = []
		for j in xrange(-ns,ns+1):
			if( j != 0):
				for i in xrange(-ns,ns+1):
					for k in xrange(-ns,ns+1):
						smear += [i*smearstep,j*smearstep,k*smearstep,1.0]
		# Deal with theta = 0.0 cases
		prj = []
		for i in xrange(-ns,ns+1):
			for k in xrange(-ns,ns+1):
				prj.append(i+k)
		for i in xrange(-2*ns,2*ns+1,1):
			smear += [i*smearstep,0.0,0.0,float(prj.count(i))]
		#if myid == 0:  print "  Smear  ",smear

	#from utilities import model_blank, get_im, read_text_file
	#bckgdata = [get_im("tsd.hdf"),read_text_file("data_stamp.txt")]

	nnx = bckgdata[0].get_xsize()
	nny = bckgdata[0].get_ysize()
	bckgnoise = []
	for i in xrange(nny):
		prj = model_blank(nnx)
		for k in xrange(nnx):  prj[k] = bckgdata[0].get_value_at(k,i)
		bckgnoise.append(prj)

	datastamp = bckgdata[1]

	#  Do the FSC shtick.
	bnx     = imgsize*npad//2+1
	refvol = model_blank(bnx)  # fill fsc with zeroes so the first reconstruction is done using simple Wiener filter.
	refvol.set_attr("fudge", 1.0)

	results_list = []
	fftvol_file =[]
	weight_file = []

	for iset in xrange(2):
		if not (finfo is None): nimg = 0

		fftvol = EMData()
		weight = EMData()
		if( smearstep > 0.0 ):  fftvol.set_attr("smear", smear)
	
		params = {"size":imgsize, "npad":npad, "snr":snr, "sign":sign, "symmetry":symmetry, "refvol":refvol, "fftvol":fftvol, "weight":weight}
		r = Reconstructors.get( "nn4_ctfw", params )
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
					ERROR("Either ptcl_source_image or ctf has to be present in the header.","recons3d_4nnw_MPI",1, myid)
			try:
				indx = datastamp.index(stmp)
			except:
				ERROR("Problem with indexing ptcl_source_image.","recons3d_4nnf_MPI",1, myid)
	
			image.set_attr("bckgnoise", bckgnoise[indx])
			insert_slices(r, image)
			if not (finfo is None):
				nimg += 1
				finfo.write(" %4d inserted\n" %(nimg) )
				finfo.flush()

		if not (finfo is None): 
			finfo.write( "begin reduce\n" )
			finfo.flush()
	
		reduce_EMData_to_root(fftvol, myid, main_node, comm=mpi_comm)
		reduce_EMData_to_root(weight, myid, main_node, comm=mpi_comm)
		
		if not (finfo is None): 
			finfo.write( "after reduce\n" )
			finfo.flush()

		if myid == 0:
			tmpid = datetime.datetime.now().strftime('%Y-%m-%d--%I-%M-%f')[:-3]
			fftvol_file.append("fftvol__%s__idx%d.hdf"%(tmpid, iset))
			weight_file.append("weight__%s__idx%d.hdf"%(tmpid, iset))
			fftvol.write_image(fftvol_file[-1])
			weight.write_image(weight_file[-1])

			dummy = r.finish(True)
			results_list.append("rvol__%s__idx%d.hdf"%(tmpid, iset))
			if(iset == 0):  fftvol.write_image(results_list[-1])

		mpi_barrier(mpi_comm)

	if myid == 0:
		fourier_shell_correlation = fsc(get_im(results_list[0]), fftvol, 1.0)[1]

		from math import sqrt
		from sp_utilities import reshape_1d
		t = [0.0]*len(fourier_shell_correlation)
		t = reshape_1d(t,len(t),npad*len(t))
		for i in xrange(len(t)):
			t[i] = min(max(t[i], 0.0), 0.999)

		ovol = []
		for idx in range(2):
			fftvol = get_im(fftvol_file[idx])
			weight = get_im(weight_file[idx])
			refvol = model_blank(bnx,1,1,0.0)
			for i in xrange(min(bnx,len(t))):  
				refvol.set_value_at(i, t[i])
			refvol.set_attr("fudge", 1.0)
			
			params = {"size":imgsize, "npad":npad, "snr":snr, "sign":sign, "symmetry":symmetry, "refvol":refvol, "fftvol":fftvol, "weight":weight}
			r = Reconstructors.get("nn4_ctfw", params)
			r.setup()
			
			dummy = r.finish(True)
			ovol.append(fftvol)


		cmd = "{} {} {} {} {} {}".format("rm -f", fftvol_file[0], fftvol_file[1], weight_file[0], weight_file[1], results_list[0] )
		import subprocess
		outcome = subprocess.call(cmd, shell=True)

	mpi_barrier(mpi_comm)
	if myid == 0:
		return ovol[0], ovol[1], fourier_shell_correlation
	else:
		return None, None, None
'''













































































































"""4
			tmpid = datetime.datetime.now().strftime('%Y-%m-%d--%I-%M-%f')[:-3]
			fftvol_file.append("fftvol__%s__idx%d.hdf"%(tmpid, iset))
			weight_file.append("weight__%s__idx%d.hdf"%(tmpid, iset))
			fftvol.write_image(fftvol_file[-1])
			weight.write_image(weight_file[-1])
			"""








"""5
		from math import sqrt
		from sp_utilities import reshape_1d
		t = [0.0]*len(fourier_shell_correlation)
		t = reshape_1d(t,len(t),npad*len(t))
		for i in xrange(len(t)):
			t[i] = min(max(t[i], 0.0), 0.999)

		ovol = []
		for idx in range(2):
			fftvol = get_im(fftvol_file[idx])
			weight = get_im(weight_file[idx])
			refvol = model_blank(bnx,1,1,0.0)
			for i in xrange(min(bnx,len(t))):  
				refvol.set_value_at(i, t[i])
			refvol.set_attr("fudge", 1.0)
			
			params = {"size":imgsize, "npad":npad, "snr":snr, "sign":sign, "symmetry":symmetry, "refvol":refvol, "fftvol":fftvol, "weight":weight}
			r = Reconstructors.get("nn4_ctfw", params)
			r.setup()
			
			dummy = r.finish(True)
			ovol.append(fftvol)


		cmd = "{} {} {} {} {} {}".format("rm -f", fftvol_file[0], fftvol_file[1], weight_file[0], weight_file[1], results_list[0] )
		import subprocess
		outcome = subprocess.call(cmd, shell=True)
		"""

























'''6
	if( smearstep > 0.0 ):
		#if myid == 0:  print "  Setting smear in prepare_recons_ctf"
		ns = 1
		smear = []
		for j in xrange(-ns,ns+1):
			if( j != 0):
				for i in xrange(-ns,ns+1):
					for k in xrange(-ns,ns+1):
						smear += [i*smearstep,j*smearstep,k*smearstep,1.0]
		# Deal with theta = 0.0 cases
		prj = []
		for i in xrange(-ns,ns+1):
			for k in xrange(-ns,ns+1):
				prj.append(i+k)
		for i in xrange(-2*ns,2*ns+1,1):
			smear += [i*smearstep,0.0,0.0,float(prj.count(i))]
		#if myid == 0:  print "  Smear  ",smear
	'''




"""7
	nnx = bckgdata[0].get_xsize()
	nny = bckgdata[0].get_ysize()
	"""
'''8
	bckgnoise = []
	for i in xrange(1):
		prj = model_blank(600,1,1,1)
		#for k in xrange(nnx):  prj[k] = bckgdata[i].get_value_at(k,i)
		bckgnoise.append(prj)
	'''

"""9
	#  Do the FSC shtick.
	if cfsc:
		bnx     = len(cfsc)*npad*2
		refvol  = model_blank(bnx)
		if(npad > 1):
			from sp_utilities import reshape_1d
			bfsc = reshape_1d(cfsc, len(cfsc), bnx)
			for i in xrange(bnx):  refvol[i] = bfsc[i]
			del bfsc
		else:  refvol[i] = cfsc[i]
	else:
		#  Set refvol to longer array so in finish it can be used to return regularization part
		refvol = model_blank(target_size)  # fill fsc with zeroes so the first reconstruction is done using simple Wiener filter.
	"""
























































































































































































































































'''10
		numbor = len(paramstructure[im][2])
		ipsiandiang = [ paramstructure[im][2][i][0]/1000  for i in xrange(numbor) ]
		allshifts   = [ paramstructure[im][2][i][0]%1000  for i in xrange(numbor) ]
		probs       = [ paramstructure[im][2][i][1] for i in xrange(numbor) ]
		#  Find unique projection directions
		tdir = list(set(ipsiandiang))
		bckgn = prjlist[im][0].get_attr("bckgnoise")
		#  For each unique projection direction:
		for ii in xrange(len(tdir)):
			#  Find the number of times given projection direction appears on the list, it is the number of different shifts associated with it.
			lshifts = findall(tdir[ii], ipsiandiang)
			toprab  = 0.0
			for ki in xrange(len(lshifts)):  toprab += probs[lshifts[ki]]
			recdata = Util.mult_scalar(prjlist[im][allshifts[lshifts[0]]], probs[lshifts[0]]/toprab)
			recdata.set_attr_dict({"padffted":1, "is_complex":0})
			for ki in xrange(1,len(lshifts)):
				Util.add_img(recdata, Util.mult_scalar(prjlist[im][allshifts[lshifts[ki]]], probs[lshifts[ki]]/toprab))
			recdata.set_attr_dict({"padffted":1, "is_complex":1})
			if not upweighted:  recdata = filt_table(recdata, bckgn )
			recdata.set_attr("bckgnoise", bckgn )
			ipsi = tdir[ii]%100000
			iang = tdir[ii]/100000
			r.insert_slice( recdata, Transform({"type":"spider","phi":refang[iang][0],"theta":refang[iang][1],"psi":refang[iang][2]+ipsi*delta}), toprab*avgnorm/norm_per_particle[im])
		'''







































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































'''  Not used anywhere?  07/29/2015  PAP11
def prepare_recons_ctf_fftvol(data, snr, symmetry, myid, main_node_half, pidlist, finfo=None, npad = 2, mpi_comm=None):
	from sp_utilities import reduce_EMData_to_root
	from EMAN2 import Reconstructors
	from mpi import MPI_COMM_WORLD

	if mpi_comm == None:
		mpi_comm = MPI_COMM_WORLD

	nx = data[0].get_xsize()

	fftvol_half = EMData()
	weight_half = EMData()
	half_params = {"size":nx, "npad":npad, "snr":snr, "sign":1, "symmetry":symmetry, "fftvol":fftvol_half, "weight":weight_half}
	half = Reconstructors.get( "nn4_ctf", half_params )
	half.setup()

	for i in pidlist:
		xform_proj = data[i].get_attr( "xform.projection" )
		half.insert_slice(data[i], xform_proj )

	if not(finfo is None):
		finfo.write( "begin reduce half\n" )
		finfo.flush()

	reduce_EMData_to_root(fftvol_half, myid, main_node_half, mpi_comm)
	reduce_EMData_to_root(weight_half, myid, main_node_half, mpi_comm)

	return fftvol_half, weight_half
'''




























































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































