











































































































"""0
			peaks = ormq_peaks(ima, cimage, xrng, yrng, step, mode, numr, cnx+sxi, cny+syi)
			[angt, sxst, syst, mirrort, peakt, select] = sim_anneal(peaks, T, step, mode, maxrin)
			[alphan, sxn, syn, mn] = combine_params2(0.0, -sxi, -syi, 0, angt, sxst, syst, mirrort)
			data[im].set_attr_dict({"select":select, "peak":peakt})
			set_params2D(data[im], [alphan, sxn, syn, mn, 1.0], ali_params)
			"""


















































'''1
#  commented out as it does not seem to be used anywhere PAP 03/02/2015
def ali2d_single_iter_fast(data, dimage, params, numr, wr, cs, tavg, cnx, cny, \
							xrng, yrng, step, maxrange, nomirror = False, mode="F", \
							random_method="", T=1.0, ali_params="xform.align2d", delta = 0.0):
	"""
		single iteration of 2D alignment using ormq
		if CTF = True, apply CTF to data (not to reference!)
	"""
	from sp_utilities import combine_params2, inverse_transform2, get_params2D, set_params2D
	from sp_alignment import ormq, ornq

	# 2D alignment using rotational ccf in polar coords and quadratic interpolation
	cimage = Util.Polar2Dm(tavg, cnx, cny, numr, mode)
	Util.Frngs(cimage, numr)
	Util.Applyws(cimage, numr, wr)

	maxrin = numr[-1]
	#sx_sum = 0
	#sy_sum = 0
	for im in xrange(len(data)):
		#alpha, sx, sy, mirror, dummy = get_params2D(data[im], ali_params)
		#alpha, sx, sy, mirror        = combine_params2(alpha, sx, sy, mirror, 0.0, -cs[0], -cs[1], 0)
		#alphai, sxi, syi, scalei     = inverse_transform2(alpha, sx, sy)
		"""
		# align current image to the reference
		if random_method == "SA":
			peaks = ormq_peaks(ima, cimage, xrng, yrng, step, mode, numr, cnx+sxi, cny+syi)
			[angt, sxst, syst, mirrort, peakt, select] = sim_anneal(peaks, T, step, mode, maxrin)
			[alphan, sxn, syn, mn] = combine_params2(0.0, -sxi, -syi, 0, angt, sxst, syst, mirrort)
			data[im].set_attr_dict({"select":select, "peak":peakt})
			set_params2D(data[im], [alphan, sxn, syn, mn, 1.0], ali_params)
		else:
			if nomirror:  [angt, sxst, syst, mirrort, peakt] = ornq(ima, cimage, xrng, yrng, step, mode, numr, cnx+sxi, cny+syi)
			else:
		"""
		[angt, sxst, syst, mirrort, peakt] = ormq_fast(dimage[im], cimage, xrng, yrng, step, params[im][1], params[im][2], maxrange, numr, mode, delta)
		# combine parameters and set them to the header, ignore previous angle and mirror
		#[alphan, sxn, syn, mn] = combine_params2(0.0, -sxi, -syi, 0, angt, sxst, syst, mirrort)
		#set_params2D(data[im], [alphan, sxn, syn, mn, 1.0], ali_params)
		params[im] = [angt, sxst, syst, mirrort]

		#if mirrort == 0: sx_sum += sxst
		#else:            sx_sum -= sxst
		#sy_sum += syst

#	return sx_sum, sy_sum
'''










'''2
def Applyws(circ, numr, wr):
	"""
	  Apply weights to FTs of rings
	"""
	nring = len(numr)/3
	maxrin = numr[len(numr)-1]
	for i in xrange(nring):
		numr3i = numr[2+i*3]
		numr2i = numr[1+i*3]-1
		w = wr[i]
		circ[numr2i] *= w
		if (numr3i == maxrin): circ[numr2i+1] *= w
		else: circ[numr2i+1] *= 0.5*w
		for j in xrange(2,numr3i):
			jc = j+numr2i
			circ[jc] *= w
'''





























































































































































































































































































































































































































































































































































































































'''3
def Numrinit(first_ring, last_ring, skip=1, mode="F"):
	#  This is to test equal length rings
	"""This function calculates the necessary information for the 2D 
	   polar interpolation. For each ring, three elements are recorded:
	   numr[i*3]:  Radius of this ring
	   numr[i*3+1]: Total number of samples of all inner rings+1
	   		(Or, the beginning point of this ring)
	   numr[i*3+2]: Number of samples of this ring. This number is an 
	   		FFT-friendly power of the 2.
			
	   "F" means a full circle interpolation
	   "H" means a half circle interpolation
	"""
	MAXFFT = 32768
	from math import pi

	if (mode == 'f' or mode == 'F'): dpi = 2*pi
	else:                            dpi = pi
	numr = []
	lcirc = 1
	#  This is for testing equal length rings
	ip = 128
	for k in xrange(first_ring, last_ring+1, skip):
		numr.append(k)
		numr.append(lcirc)
		numr.append(ip)
		lcirc += ip		
	return  numr
'''


































































































'''4
			# The following code is used when mirror is not considered
			retvals = Util.Crosrng_e(crefim, cimage, numr, 0, 0.0)
			qn = retvals["qn"]
			if qn >= peak:
				sx = -ix
				sy = -iy
				ang = ang_n(retvals["tot"], mode, numr[-1])
				peak = qn
				mirror = 0
			'''

















"""5
	lkx = int(xrng[0]*istep)
	rkx = int(xrng[-1]*istep)

	lky = int(yrng[0]*istep)
	rky = int(yrng[-1]*istep)
	"""






















"""6
	co =  cos(radians(ang))
	so = -sin(radians(ang))
	sxs = sx*co - sy*so
	sys = sx*so + sy*co
	"""



































































'''7
def process_peak_1d_pad(peaks, step, mode, numr, nx):
	from math import pi, cos, sin
	
	peak_num = len(peaks)
	maxrin = numr[-1]
	for i in xrange(peak_num):
		peaks[i][1]= ang_n(peaks[i][1]+1-nx/2, mode, maxrin)
	peaks.sort(reverse=True)
	
	return peaks

def find_position(list_a, t):
	"""
	The function determines how many elements of list_a is larger or equal than t.
	Here we assume list_a is sorted reversely.
	"""
	if list_a[0] < t:
		return 0
	elif list_a[-1] >= t:
		return len(list_a)
	else:
		K = len(list_a)
		k_min = 1
		k_max = K-1
		while k_min != k_max:
			k = (k_min+k_max)/2
			if list_a[k] < t:
				if list_a[k-1] >= t:
					k_min = k
					k_max = k
				else:
					k_max = k-1
			else:
				k_min = k+1
		return k_min


def select_major_peaks(g, max_major_peaks, min_height, dim):

	from sp_filter import filt_gaussl
	from sp_fundamentals import fft
	from sp_utilities import peak_search
	
	G = fft(g)
	
	found = False
	min_fl = 0.001
	max_fl = 0.5
	
	while found == False:
		fl = (min_fl+max_fl)/2
		peakg = peak_search(fft(filt_gaussl(G, fl)), 1000)
		K = len(peakg)
		list_a = [0.0]*K
		for i in xrange(K):  list_a[i] = peakg[i][dim+1]
		k = find_position(list_a, min_height)
		if k > max_major_peaks: 
			max_fl = fl
		elif k < max_major_peaks:
			min_fl = fl
		else:
			found = True
		if max_fl - min_fl < 0.001: found = True
		 
	return peakg[0:k] 


def select_major_peaks_Gaussian_fitting(peak):

	# Generate the histogram of the angle
	ang_bin = [0]*30
	for i in xrange(len(angle)):
		#angle.append(peak[i][1])
		bin_num = int(angle[i]/12)
		ang_bin[bin_num] += 1
	ang_bin_index = []
	for i in xrange(30):
		ang_bin_index.append([ang_bin[i], i])
	ang_bin_index.sort(reverse=True)
	print ang_bin
	print ang_bin_index
	
	K = 5
	A = [0.0]*K
	mu = [0.0]*K
	sigma = [0.0]*K
	
	for k in xrange(K):
		A[k] = ang_bin_index[k][0]
		mu[k] = ang_bin_index[k][1]*12+6
		sigma[k] = 5.0
	
	print A, mu, sigma 
	
	
	return []


def ormq_peaks_major(image, crefim, xrng, yrng, step, mode, numr, cnx, cny):
	"""
	Determine shift and rotation between image and reference image (crefim)
	crefim should be as FT of polar coords with applied weights
	consider mirror
	quadratic interpolation
	cnx, cny in FORTRAN convention
	"""
	from sp_utilities import peak_search, pad
	
	ccfs = EMData()
	ccfm = EMData()
	ccfs_compress = EMData()
	ccfm_compress = EMData()

	Util.multiref_peaks_compress_ali2d(image, crefim, xrng, yrng, step, mode, numr, cnx, cny, ccfs, ccfm, ccfs_compress, ccfm_compress)
	
	nx = ccfs.get_xsize()
	ny = ccfs.get_ysize()
	nz = ccfs.get_zsize()

	peaks = peak_search(ccfs, 1000)
	peakm = peak_search(ccfm, 1000)

	peaks = process_peak(peaks, step, mode, numr)
	peakm = process_peak(peakm, step, mode, numr)

	max_major_peaks = 5
	min_height = 0.7
	
	peaks_major = select_major_peaks(pad(ccfs_compress, nx*2), max_major_peaks, min_height, 1)
	peakm_major = select_major_peaks(pad(ccfm_compress, nx*2), max_major_peaks, min_height, 1)

	peaks_major = process_peak_1d_pad(peaks_major, step, mode, numr, nx)
	peakm_major = process_peak_1d_pad(peakm_major, step, mode, numr, nx)	

	"""
	ccfs_compress = EMData(nx, 1, 1, True)
	ccfm_compress = EMData(nx, 1, 1, True)

	for x in xrange(nx):
		slices = Util.window(ccfs, 1, ny-2, nz-2, x-nx/2, 0, 0)
		slicem = Util.window(ccfm, 1, ny-2, nz-2, x-nx/2, 0, 0)
		
		[means, dummy, dummy, dummy] = Util.infomask(slices, None, True)
		[meanm, dummy, dummy, dummy] = Util.infomask(slicem, None, True)
		
		ccfs_compress.set_value_at(x, 0, 0, means)
		ccfm_compress.set_value_at(x, 0, 0, meanm)

	peaks_major = select_major_peaks_Gaussian_fitting(peaks)
	peakm_major = select_major_peaks_Gaussian_fitting(peakm)
	
	fs = Util.window(ccfs, nx, ny-2, nz-2, 0, 0, 0)
	fm = Util.window(ccfm, nx, ny-2, nz-2, 0, 0, 0)

	dummy1, dummy2, mins, dummy3 = Util.infomask(fs, None, True)
	dummy1, dummy2, minm, dummy3 = Util.infomask(fm, None, True)

	gs = threshold_to_minval(ccfs, mins)
	gm = threshold_to_minval(ccfm, minm)

	peaks_major = select_major_peaks(gs, max_major_peaks, min_height, 3)
	peakm_major = select_major_peaks(gm, max_major_peaks, min_height, 3)
	
	peaks_major = select_major_peaks(pad(ccfs_compress, nx*2), max_major_peaks, min_height, 1)
	peakm_major = select_major_peaks(pad(ccfm_compress, nx*2), max_major_peaks, min_height, 1)
	time4 = time()
	peaks = peak_search(ccfs, 1000)
	peakm = peak_search(ccfm, 1000)
	time5 = time()
	peaks = process_peak(peaks, step, mode, numr)
	peakm = process_peak(peakm, step, mode, numr)
	peaks_major = process_peak_1d_pad(peaks_major, step, mode, numr)
	peakm_major = process_peak_1d_pad(peakm_major, step, mode, numr)	
	time6 = time()
	"""
	
	return peaks, peakm, peaks_major, peakm_major
'''












































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































"""8
	co =  cos(radians(ang))
	so = -sin(radians(ang))
	sxs = sx*co - sy*so
	sys = sx*so + sy*co
	"""






































'''9
def prepare_refproj(volprep, refang, delta_psi = 1.0, mempercpu = 1.e9, kb3D = None):
	from sp_projection import prgs,prgl
	from sp_fundamentals import fft
	from math import sqrt
	ny = volprep.get_ysize()
	if kb3D:  ny /= 2
	npsi = int(360./delta_psi)
	nang = len(refang)
	
	if( (4.0*(ny+2)*ny)*nang*npsi < mempercpu ):
		refproj = []
		nrefproj = 0
		for i in xrange(nang):
			###if myid == main_node:  print "  Angle :",i,time()-at
			for j in xrange(npsi):
				psi = j*delta_psi
				if kb3D:  temp = prgs(volprep, kb3D, [refang[i][0],refang[i][1],psi, 0.0,0.0])
				else:     temp = prgl(volprep,[ refang[i][0],refang[i][1],psi, 0.0,0.0], 1, False)
				temp.set_attr("is_complex",0)
				nrmref = sqrt(temp.cmp("dot", temp, dict(negative = 0)))
				refproj.append([temp, nrmref])
				nrefproj += 1
	else:
		refproj = None
	return refproj
'''


















































'''10
					toto = -1
					for k in xrange(lentop):
						if(peak > newpar[kl][-1][k][1]):
							toto = k
							break
					if( toto == 0 ):  newpar[kl][-1] = [[im + iangpsi, peak]] + newpar[kl][-1][:lentop-1]
					elif(toto > 0 ):  newpar[kl][-1] = newpar[kl][-1][:toto-1] + [[im + iangpsi, peak]] + newpar[kl][-1][toto:lentop-1]
					'''



















































































'''11
						toto = -1
						for k in xrange(lentop):
							if(peak > newpar[kl][-1][k][1]):
								toto = k
								break
						if( toto == 0 ):  newpar[kl][-1] = [[im + iangpsi, peak]] + newpar[kl][-1][:lentop-1]
						elif(toto > 0 ):  newpar[kl][-1] = newpar[kl][-1][:toto-1] + [[im + iangpsi, peak]] + newpar[kl][-1][toto:lentop-1]
						'''





















































































'''12
							toto = -1
							for k in xrange(lentop):
								if(peak > newpar[kl][-1][k][1]):
									toto = k
									break
							if( toto == 0 ):  newpar[kl][-1] = [[im + iangpsi, peak]] + newpar[kl][-1][:lentop-1]
							elif(toto > 0 ):  newpar[kl][-1] = newpar[kl][-1][:toto-1] + [[im + iangpsi, peak]] + newpar[kl][-1][toto:lentop-1]
							'''



















































































































































































































































































































































































"""13
	[ang, sxs, sys, mirror, iref, peak] = \
		Util.multiref_polar_ali_helicon_90_local(data, refrings, xrng, yrng, stepx, ant, psi_max, mode, numr, cnx-tx, cny-ty, int(ynumber), yrnglocal)
	"""


























































































































































































































































































































































































































'''14
def align2dshc(image, refim, xrng=0, yrng=0, step=1, first_ring=1, last_ring=0, rstep=1, mode = "F"):
	"""  Determine shift and rotation between image and reference image
	     quadratic interpolation
	     Output: ang, sxs, sys, mirror, peak
	"""
	#from utilities import print_col
	from sp_alignment import Numrinit, ringwe
	step = float(step)
	nx = refim.get_xsize()
	ny = refim.get_ysize()
	MAX_XRNG = nx/2
	MAX_YRNG=ny/2
	if xrng >= MAX_XRNG:
		ERROR('Translation search range in x is at most %d'%MAX_XRNG, "align2d ", 1)
	if yrng >= MAX_YRNG:
		ERROR('Translation search range in y is at most %d'%MAX_YRNG, "align2d ", 1)
	if(last_ring == 0):  last_ring = nx/2-2-int(max(xrng,yrng))
	# center in SPIDER convention
	cnx = nx//2+1
	cny = ny//2+1
	#precalculate rings
	numr = Numrinit(first_ring, last_ring, rstep, mode)
	wr   = ringwe(numr, mode)
	#cimage=Util.Polar2Dmi(refim, cnx, cny, numr, mode, kb)
	crefim = Util.Polar2Dm(refim, cnx, cny, numr, mode)
	#crefim = Util.Polar2D(refim, numr, mode)
	#print_col(crefim)
	Util.Frngs(crefim, numr)
	Util.Applyws(crefim, numr, wr)
	#return ormq(image, crefim, xrng, yrng, step, mode, numr, cnx, cny)
	return   Util.shc(image, [crefim], xrng, yrng, step, -1.0, mode, numr, cnx, cny, "c1")
'''















































































































































































































































































'''15
def parabl(Z):
	#  Original with Fortran indexing

	C1 = (26.*Z[1,1] - Z[1,2] + 2*Z[1,3] - Z[2,1] - 19.*Z[2,2] - 7.*Z[2,3] + 2.*Z[3,1] - 7.*Z[3,2] + 14.*Z[3,3])/9.

	C2 = (8.* Z[1,1] - 8.*Z[1,2] + 5.*Z[2,1] - 8.*Z[2,2] + 3.*Z[2,3] +2.*Z[3,1] - 8.*Z[3,2] + 6.*Z[3,3])/(-6.)

	C3 = (Z[1,1] - 2.*Z[1,2] + Z[1,3] + Z[2,1] -2.*Z[2,2] + Z[2,3] + Z[3,1] - 2.*Z[3,2] + Z[3,3])/6.

	C4 = (8.*Z[1,1] + 5.*Z[1,2] + 2.*Z[1,3] -8.*Z[2,1] -8.*Z[2,2] - 8.*Z[2,3] + 3.*Z[3,2] + 6.*Z[3,3])/(-6.)

	C5 = (Z[1,1] - Z[1,3] - Z[3,1] + Z[3,3])/4.

	C6 = (Z[1,1] + Z[1,2] + Z[1,3] - 2.*Z[2,1] - 2.*Z[2,2] -2.*Z[2,3] + Z[3,1] + Z[3,2] + Z[3,3])/6.

	DENOM = 4. * C3 * C6 - C5 * C5
	if(DENOM == 0.):
		return 0.0, 0.0, 0.0

	YSH   = (C4*C5 - 2.*C2*C6) / DENOM - 2.
	XSH   = (C2*C5 - 2.*C4*C3) / DENOM - 2.

	PEAKV = 4.*C1*C3*C6 - C1*C5*C5 - C2*C2*C6 + C2*C4*C5 - C4*C4*C3
	PEAKV = PEAKV / DENOM
	#print  "  in PARABL  ",XSH,YSH,Z[2,2],PEAKV

	XSH = min(max( XSH, -1.0), 1.0)
	YSH = min(max( YSH, -1.0), 1.0)

	return XSH, YSH, PEAKV
'''




























































































































































































































































































































"""16
	fft(ima).write_image('ima.hdf')
	for i in xrange(nr):  fft(ref[i]).write_image('ref.hdf',i)
	from sys import exit
	exit()
	"""
























"""17
				if(pp[0]>ma1):
					ma1 = pp[0]
					oma1 = pp+[XSH, YSH,int(pp[4])+XSH, int(pp[5])+YSH, PEAKV,(i-nc)*psistep]
				"""
























"""18
				if(pp[0]>ma3):
					ma3 = pp[0]
					oma3 = pp+[XSH, YSH,int(pp[4])+XSH, int(pp[5])+YSH, PEAKV,(i-nc)*psistep]
				"""








"""19
		print oma1
		print oma2
		print  "        %6.2f %6.2f  %6.2f"%(oma2[-1],oma2[-4],oma2[-3])
		"""




"""20
		print oma3
		print oma4
		"""












































































































"""21
		N = refs[0].get_ysize()  # assumed square image
		# prepare 
		#npad = 2
		#N = nx*npad
		K = 6
		alpha = 1.75
		r = nx/2
		v = K/2.0/N
		kb = Util.KaiserBessel(alpha, K, r, v, N)
		"""









"""22
			ref[i] = rot_shift2D(refs,(i-nc)*psistep, interpolation_method = 'gridding')
			ref[i] = ref[i].FourInterpol(N, N, 1,0)
			#fft(ref[i]).write_image('refprep.hdf')
			"""

































"""23
	fft(ima).write_image('ima.hdf')
	for i in xrange(nr):  fft(ref[i]).write_image('ref.hdf',i)
	from sys import exit
	exit()
	"""






















"""24
				loc = w.calc_max_location()
				PEAKV = w.get_value_at(loc[0],loc[1])
				print "  Did not find a peak  :",i,loc[0]-wxc, loc[1]-wyc, PEAKV
				if(PEAKV>ma2):
						ma2  = PEAKV
						oma2 = pp+[loc[0]-wxc, loc[1]-wyc, loc[0]-wxc, loc[1]-wyc, PEAKV,(i-nc)*psistep]
				"""









"""25
				if(pp[0]>ma1):
					ma1 = pp[0]
					oma1 = pp+[XSH, YSH,int(pp[4])+XSH, int(pp[5])+YSH, PEAKV,(i-nc)*psistep]
				"""














"""26
				loc = w.calc_max_location()
				PEAKV = w.get_value_at(loc[0],loc[1])
				if(PEAKV>ma4):
					ma4  = PEAKV
					oma4 = pp+[loc[0], loc[1], loc[0], loc[1], PEAKV,(i-nc)*psistep]
				"""









"""27
				if(pp[0]>ma3):
					ma3 = pp[0]
					oma3 = pp+[XSH, YSH,int(pp[4])+XSH, int(pp[5])+YSH, PEAKV,(i-nc)*psistep]
				"""








"""28
		print oma1
		print oma2
		print  "        %6.2f %6.2f  %6.2f"%(oma2[-1],oma2[-4],oma2[-3])
		"""




"""29
		print oma3
		print oma4
		"""






















"""30
	M = inima.get_xsize()
	alpha = 1.75
	K = 6
	N = M*2  # npad*image size
	r = M/2
	v = K/2.0/N
	params = {"filter_type" : Processor.fourier_filter_types.KAISER_SINH_INVERSE,
	          "alpha" : alpha, "K":K,"r":r,"v":v,"N":N}
	kb = Util.KaiserBessel(alpha, K, r, v, N)
	"""






































"""31
	fft(ima).write_image('ima.hdf')
	for i in xrange(nr):  fft(ref[i]).write_image('ref.hdf',i)
	from sys import exit
	exit()
	"""






















"""32
				loc = w.calc_max_location()
				PEAKV = w.get_value_at(loc[0],loc[1])
				print "  Did not find a peak  :",i,loc[0]-wxc, loc[1]-wyc, PEAKV
				if(PEAKV>ma2):
						ma2  = PEAKV
						oma2 = pp+[loc[0]-wxc, loc[1]-wyc, loc[0]-wxc, loc[1]-wyc, PEAKV,(i-nc)*psistep]
				"""









"""33
				if(pp[0]>ma1):
					ma1 = pp[0]
					oma1 = pp+[XSH, YSH,int(pp[4])+XSH, int(pp[5])+YSH, PEAKV,(i-nc)*psistep]
				"""














"""34
				loc = w.calc_max_location()
				PEAKV = w.get_value_at(loc[0],loc[1])
				if(PEAKV>ma4):
					ma4  = PEAKV
					oma4 = pp+[loc[0], loc[1], loc[0], loc[1], PEAKV,(i-nc)*psistep]
				"""









"""35
				if(pp[0]>ma3):
					ma3 = pp[0]
					oma3 = pp+[XSH, YSH,int(pp[4])+XSH, int(pp[5])+YSH, PEAKV,(i-nc)*psistep]
				"""








"""36
		print oma1
		print oma2
		print  "        %6.2f %6.2f  %6.2f"%(oma2[-1],oma2[-4],oma2[-3])
		"""












































"""37
	M = inima.get_xsize()
	alpha = 1.75
	K = 6
	N = M*2  # npad*image size
	r = M/2
	v = K/2.0/N
	params = {"filter_type" : Processor.fourier_filter_types.KAISER_SINH_INVERSE,
	          "alpha" : alpha, "K":K,"r":r,"v":v,"N":N}
	kb = Util.KaiserBessel(alpha, K, r, v, N)
	"""


























































"""38
	fft(ima).write_image('ima.hdf')
	for i in xrange(nr):  fft(ref[i]).write_image('ref.hdf',i)
	from sys import exit
	exit()
	"""



























"""39
				loc = w.calc_max_location()
				PEAKV = w.get_value_at(loc[0],loc[1])
				#print "  Did not find a peak  :",i,loc[0]-wxc, loc[1]-wyc, PEAKV
				if(PEAKV>ma2):
						ma2  = PEAKV
						oma2 = pp+[loc[0]-wxc, loc[1]-wyc, loc[0]-wxc, loc[1]-wyc, PEAKV,(i-nc)*psistep]
				"""









"""40
				if(pp[0]>ma1):
					ma1 = pp[0]
					oma1 = pp+[XSH, YSH,int(pp[4])+XSH, int(pp[5])+YSH, PEAKV,(i-nc)*psistep]
				"""

















"""41
				loc = w.calc_max_location()
				PEAKV = w.get_value_at(loc[0],loc[1])
				if(PEAKV>ma4):
					ma4  = PEAKV
					oma4 = pp+[loc[0], loc[1], loc[0], loc[1], PEAKV,(i-nc)*psistep]
				"""









"""42
				if(pp[0]>ma3):
					ma3 = pp[0]
					oma3 = pp+[XSH, YSH,int(pp[4])+XSH, int(pp[5])+YSH, PEAKV,(i-nc)*psistep]
				"""











"""43
		print oma1
		print oma2
		print  "        %6.2f %6.2f  %6.2f"%(oma2[-1],oma2[-4],oma2[-3])
		"""












































"""44
	M = inima.get_xsize()
	alpha = 1.75
	K = 6
	N = M*2  # npad*image size
	r = M/2
	v = K/2.0/N
	params = {"filter_type" : Processor.fourier_filter_types.KAISER_SINH_INVERSE,
	          "alpha" : alpha, "K":K,"r":r,"v":v,"N":N}
	kb = Util.KaiserBessel(alpha, K, r, v, N)
	"""

















































"""45
	fft(ima).write_image('ima.hdf')
	for i in xrange(nr):  fft(ref[i]).write_image('ref.hdf',i)
	from sys import exit
	exit()
	"""

















































"""46
				if(pp[0]>ma1):
					ma1 = pp[0]
					oma1 = pp+[XSH, YSH,int(pp[4])+XSH, int(pp[5])+YSH, PEAKV,(i-nc)*psistep]
				"""







































"""47
				if(pp[0]>ma3):
					ma3 = pp[0]
					oma3 = pp+[XSH, YSH,int(pp[4])+XSH, int(pp[5])+YSH, PEAKV,(i-nc)*psistep]
				"""









"""48
		print oma1
		print oma2
		print  "        %6.2f %6.2f  %6.2f"%(oma2[-1],oma2[-4],oma2[-3])
		"""


















































































































































































































































































































































'''49
def Xshc0(data, cimages, refrings, numr, xrng, yrng, step, an = -1.0, sym = "c1", finfo=None):
	from sp_utilities    import compose_transform2
	from math         import cos, sin, degrees, radians
	from EMAN2 import Vec2f

	number_of_checked_refs = 0

	mode = "F"
	nx   = data.get_xsize()
	ny   = data.get_ysize()
	#  center is in SPIDER convention
	cnx  = nx//2 + 1
	cny  = ny//2 + 1

	if( an>= 0.0):  ant = cos(radians(an))
	else:           ant = -1.0
	#phi, theta, psi, sxo, syo = get_params_proj(data)
	t1 = data.get_attr("xform.projection")
	#dp = t1.get_params("spider")
	if finfo:
		#finfo.write("Old parameters: %9.4f %9.4f %9.4f %9.4f %9.4f\n"%(phi, theta, psi, sxo, syo))
		finfo.write("Old parameters: %9.4f %9.4f %9.4f %9.4f %9.4f\n"%(dp["phi"], dp["theta"], dp["psi"], -dp["tx"], -dp["ty"]))
		finfo.flush()

	previousmax = data.get_attr("previousmax")
	cimages[0].set_attr("xform.projection", t1)
	cimages[0].set_attr("previousmax", previousmax)
	#  The code for shc does not work for local searches!  PAP 01/27/2015
	#  Do not use previous shifts so the image does not slide away
	[ang, sxs, sys, mirror, iref, peak, checked_refs] = Util.shc0(cimages, refrings, xrng, yrng, step, ant, mode, numr, cnx, cny, sym )  #+dp["tx"], cny+dp["ty"])
	iref=int(iref)
	number_of_checked_refs += int(checked_refs)
	#[ang,sxs,sys,mirror,peak,numref] = apmq_local(projdata[imn], ref_proj_rings, xrng, yrng, step, ant, mode, numr, cnx-sxo, cny-syo)
	#ang = (ang+360.0)%360.0

	if peak <= previousmax:
		return -1.0e23, 0.0, number_of_checked_refs, -1
		"""
		# there is no better solutions - if the current position is free, we don't change anything
		last_phi = dp["phi"]
		last_theta = dp["theta"]
		found_current_location = False
		for ir in xrange(len(refrings)):
			r = refrings[ir]
			if abs(last_phi - r.get_attr("phi")) < 0.1 and abs(last_theta - r.get_attr("theta")) < 0.1:
				found_current_location = True
				break
		if found_current_location:
			return -1.0e23, 0.0, number_of_checked_refs, ir
		"""
	else:
		# The ormqip returns parameters such that the transformation is applied first, the mirror operation second.
		# What that means is that one has to change the the Eulerian angles so they point into mirrored direction: phi+180, 180-theta, 180-psi
		angb, sxb, syb, ct = compose_transform2(0.0, sxs, sys, 1, -ang, 0.0, 0.0, 1)
		if  mirror:
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

		#set_params_proj(data, [phi, theta, psi, s2x, s2y])
		t2 = Transform({"type":"spider","phi":phi,"theta":theta,"psi":psi})
		t2.set_trans(Vec2f(-s2x, -s2y))
		data.set_attr("xform.projection", t2)
		data.set_attr("previousmax", peak)
		#  Find the pixel error that is minimum over symmetry transformations
		from sp_pixel_error import max_3D_pixel_error
		ts = t2.get_sym_proj(sym)
		if(len(ts) > 1):
			# only do it if it is not c1
			pixel_error = +1.0e23
			for ut in ts:
				# we do not care which position minimizes the error
				pixel_error = min(max_3D_pixel_error(t1, ut, numr[-3]), pixel_error)
		else:
			pixel_error = max_3D_pixel_error(t1, t2, numr[-3])
		if finfo:
			finfo.write( "New parameters: %9.4f %9.4f %9.4f %9.4f %9.4f %10.5f  %11.3e\n\n" %(phi, theta, psi, s2x, s2y, peak, pixel_error))
			finfo.flush()
		return peak, pixel_error, number_of_checked_refs, iref
'''








































"""50
		# there is no better solutions - if the current position is free, we don't change anything
		last_phi = dp["phi"]
		last_theta = dp["theta"]
		found_current_location = False
		for ir in xrange(len(refrings)):
			r = refrings[ir]
			if abs(last_phi - r.get_attr("phi")) < 0.1 and abs(last_theta - r.get_attr("theta")) < 0.1:
				found_current_location = True
				break
		if found_current_location:
			return -1.0e23, 0.0, number_of_checked_refs, ir
		"""




























































































"""51
	if an == "-1":
		an = [-1] * lstp
	else:
		an = get_input_from_string(an)
	"""






"""52
		import os
		outdir = "./"
		info_file = os.path.join(outdir, "progress%04d"%myid)
		finfo = open(info_file, 'w')
		"""


































































































































































































'''53
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

def total_size_of_object_in_memory(my_object):
	import inspect, os, subprocess
	from sp_utilities import random_string

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
	fp.write("from sp_sparx import *\n")

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
		output += 0.1*int(subprocess.Popen(["./%s"%file_name_my_python_code], stdout = subprocess.PIPE, stderr = subprocess.STDOUT, env = current_env).communicate()[0])
	os.system("rm ./%s"%file_name_my_python_code)
	os.system("rm ./%s"%file_name_my_object)
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
	from sp_alignment import prepare_refrings, refprojs, Numrinit, ringwe
	from sp_morphology import bracket_def, goldsearch_astigmatism
	from sp_applications import computenumberofrefs
	from sp_utilities import even_angles, assign_projangles, cone_ang, print_from_process
	
	
	LOW_LIMIT_FOR_NUMBER_OF_REFERENCES_THAT_FIT_MEMORY = 100
	# FRACTION_OF_MEMORY_THAT_CAN_BE_ALLOCATED = 0.9 # do not allocate all available memory
	# FRACTION_OF_MEMORY_THAT_CAN_BE_ALLOCATED = 0.000125 # yields about 21 cones
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
	from sp_utilities import if_error_then_all_processes_exit_program
	error_status = mpi_reduce(error_status, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD)
	if_error_then_all_processes_exit_program(error_status)	

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


def generate_indices_and_refrings(nima, projangles, volft, kb, nx, delta, an, rangle, ref_a, sym, numr, MPI, phiEqpsi = "Zero"):
	
	from sp_alignment import prepare_refrings, refprojs, Numrinit, ringwe, generate_list_of_reference_angles_for_search
	from sp_alignment import reduce_indices_so_that_angles_map_only_to_asymmetrix_unit_and_keep_mirror_info
	from sp_morphology import bracket_def, goldsearch_astigmatism
	from sp_applications import computenumberofrefs
	from sp_utilities import even_angles, assign_projangles_f, assign_projangles
	from sp_utilities import cone_ang_with_index
	import sys
	from sp_projection import prep_vol

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
			from sp_utilities import rotate_shift_params
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
		
		# use at least 10 cones
		if number_of_cones > 1 and number_of_cones < 10:
			number_of_cones = 10
		
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
				from sp_utilities import rotate_shift_params
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




### end: code that supports cone implementation
########################################################################################################################

def frame_alignment(movie_stack, particle_radius, templates, x_half_size, psi_half_size, y_half_size = None, apply_alignment_in_place = False):
	
	from sp_utilities import model_circle, list_prod, calculate_space_size
	from sp_statistics import ccc
	import numpy as np
	from sp_fundamentals import rot_shift2D
	
	if y_half_size == None:
		y_half_size = x_half_size
	
	NUMBER_OF_FRAMES = len(movie_stack)
	# x_half_size, y_half_size, psi_half_size = 5, 5, 5
	image_movement_space_size = x_length, y_length, psi_length = calculate_space_size(x_half_size, y_half_size, psi_half_size)
	
	nx = movie_stack[0].get_xsize()
	mask = model_circle(particle_radius, nx, nx)
	
	space_size = [NUMBER_OF_FRAMES] + image_movement_space_size
	new_ccEMData = EMData(list_prod(space_size), 1)

	for k in range(NUMBER_OF_FRAMES):
		# print "k", k
		for psi_i in range(psi_length):
			for y_i in range(y_length):
				for x_i in range(x_length):
					Util.write_nd_array(new_ccEMData, space_size, [k, x_i, y_i, psi_i], ccc(movie_stack[k], templates[x_i][y_i][psi_i], mask))
	
	result = Util.max_sum_along_line_in_nd_array(new_ccEMData, space_size, NUMBER_OF_FRAMES)
	
	if apply_alignment_in_place:
		[x_i, y_i, psi_i, x_j, y_j, psi_j] = result
		xxx = np.linspace(-x_half_size, x_half_size, x_length)
		yyy = np.linspace(-y_half_size, y_half_size, y_length)
		ppp = np.linspace(-psi_half_size, psi_half_size, psi_length)
		
		x = np.linspace(xxx[x_i], xxx[x_j], NUMBER_OF_FRAMES)
		y = np.linspace(yyy[y_i], yyy[y_j], NUMBER_OF_FRAMES)
		psi = np.linspace(ppp[psi_i], ppp[psi_j], NUMBER_OF_FRAMES)
		
		
		for i in range(NUMBER_OF_FRAMES):
			movie_stack[i] = rot_shift2D(movie_stack[i], -psi[i], -x[i], -y[i], 0, interpolation_method="linear")
			return None
			
	else:
		return result
'''

'''54
def XXali2d_single_iter(data, numr, wr, cs, tavg, cnx, cny, \
						xrng, yrng, step, nomirror = False, mode="F", CTF=False, \
						random_method="", T=1.0, ali_params="xform.align2d", delta = 0.0):
	"""
		single iteration of 2D alignment using ormq
		if CTF = True, apply CTF to data (not to reference!)
	"""
	from sp_utilities import combine_params2, inverse_transform2, get_params2D, set_params2D
	from sp_alignment import ormq, ornq

	if CTF:
		from sp_filter  import filt_ctf

	# 2D alignment using rotational ccf in polar coords and quadratic interpolation
	cimage = Util.Polar2Dm(tavg, cnx, cny, numr, mode)
	Util.Frngs(cimage, numr)
	Util.Applyws(cimage, numr, wr)

	maxrin = numr[-1]
	sx_sum = 0.0
	sy_sum = 0.0
	for im in xrange(len(data)):
		if CTF:
			#Apply CTF to image
			ctf_params = data[im].get_attr("ctf")
			ima = filt_ctf(data[im], ctf_params, True)
		else:
			ima = data[im]
		alpha, sx, sy, mirror, dummy = get_params2D(data[im], ali_params)
		alpha, sx, sy, mirror        = combine_params2(alpha, sx, sy, mirror, 0.0, -cs[0], -cs[1], 0)
		alphai, sxi, syi, scalei     = inverse_transform2(alpha, sx, sy)

		# align current image to the reference
		if random_method == "SA":
			peaks = ormq_peaks(ima, cimage, xrng, yrng, step, mode, numr, cnx+sxi, cny+syi)
			[angt, sxst, syst, mirrort, peakt, select] = sim_anneal(peaks, T, step, mode, maxrin)
			[alphan, sxn, syn, mn] = combine_params2(0.0, -sxi, -syi, 0, angt, sxst, syst, mirrort)
			data[im].set_attr_dict({"select":select, "peak":peakt})
			set_params2D(data[im], [alphan, sxn, syn, mn, 1.0], ali_params)
		else:
			if nomirror:  [angt, sxst, syst, mirrort, peakt] = ornq(ima, cimage, xrng, yrng, step, mode, numr, cnx+sxi, cny+syi)
			else:	      [angt, sxst, syst, mirrort, peakt] = ormq(ima, cimage, [xrng,xrng], [yrng,yrng], step, mode, numr, cnx+sxi, cny+syi, delta)
			# combine parameters and set them to the header, ignore previous angle and mirror
			[alphan, sxn, syn, mn] = combine_params2(0.0, -sxi, -syi, 0, angt, sxst, syst, mirrort)
			set_params2D(data[im], [alphan, sxn, syn, mn, 1.0], ali_params)

		if mn == 0: sx_sum += sxn
		else:       sx_sum -= sxn
		sy_sum += syn

	return sx_sum, sy_sum
'''




