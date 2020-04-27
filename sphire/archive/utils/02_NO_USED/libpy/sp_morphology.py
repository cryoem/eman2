































































































































































































































































from __future__ import print_function
def invert(im):
	"""
	 Invert contrast of an image (while keeping the average)
	"""
	p = EMAN2_cppwrap.Util.infomask(im, None, True)[0]
	return ((-1.0*im) + 2*p)

#def compress(img, value = 0.0, frange=1.0):
#	return img.process( "threshold.compress", {"value": value, "range": frange } )

def expn(img, a = 1.0, b=0.0):
	"""
		Name
			expn - generate image whose pixels are generated of raising to a given power pixels of the input image
		Input
			image: input real image
		Output
			the output image whose pixels are given by o=ir
			r: exponent
	"""
	return img.process( "math.exp", {"low": 1.0/a, "high":b})













def alog10(img):
	return img.process( "math.log")





































def threshold_to_zero(img, minval = 0.0):
	"""
		Name
			threshold_to_zero - replace values below given threshold by zero and values above by (value-threshold)
		Input
			img: input image
			minval: value below which image pixels will be set to zero
		Output
			thresholded image.
	"""
	return img.process( "threshold.belowtozero_cut", {"minval": minval } )

def threshold_to_minval(img, minval = 0.0):
	"""
		Name
			threshold_to_minval - replace values below given threshold by the threshold value
		Input
			img: input image
			minval: value below which image pixels will be set to this value
		Output
			thresholded image.
	"""
	return img.process( "threshold.belowtominval", {"minval": minval } )












def threshold_inside(img, minval, maxval):
	"""
		Name
			threshold_inside - replace values inside given thresholds by zeroes
		Input
			img: input image
			minval, maxval: image pixels that have values between these thresholds will be set to zero.
		WARNING: This function is written in python and thus very slow for large images/volumes
	"""
	im = img.copy()
	nx = im.get_xsize()
	ny = im.get_ysize()
	nz = im.get_zsize()
	for z in range(nz):
		for y in range(ny):
			for x in range(nx):
				q = im.get_value_at(x,y,z)
				if( q>minval and q < maxval): im.set_value_at(x,y,z,0.0)
	return im

def threshold_maxval(img, maxval = 0.0):
	"""
		Name
			threshold_maxval - replace values above given threshold by the threshold value
		Input
			img: input image
			maxval: value above which image pixels will be set to this value
		Output
			thresholded image.
	"""
	st = EMAN2_cppwrap.Util.infomask(img, None, True)
	return img.process( "threshold.clampminmax", {"minval": st[2], "maxval": maxval } )

def notzero(img):
	"""
		Name
			notzero - replace values that are not zero by 1.0
		Input
			img: input image
		Output
			binary image.
	"""
	return img.process( "threshold.notzero" )

def linchange(a, fct):
	"""
	reinterpolate a line given as a list by a factor of fct.
	Useful for adjusting 1D power spectra, uses linear interplation
	"""
	fctf = float(fct)
	n = len(a)
	m = int(n*fctf+0.5)
	o = [0.0]*m
	for i in range(m):
		x = i/fctf
		j = min(int(x), n-2)
		dx = x-j
		o[i] = (1.0-dx)*a[j] + dx*a[j+1]
	return o

## Fitting of "ideal" 3D FSC, as elaborated in Resolution, Meth Enz, 2010, Eq.3.25













































































































































































































































































































def compare_ctfs(nx, ctf1, ctf2):
	sign = 1.0
	dict       = ctf1.to_dict()
	dz         = dict["defocus"]
	cs         = dict["cs"]
	voltage    = dict["voltage"]
	pixel_size = dict["apix"]
	b_factor   = dict["bfactor"]
	ampcont    = dict["ampcont"]
	dza        = dict["dfdiff"]
	azz        = dict["dfang"]
	cim1 = EMAN2_cppwrap.Util.ctf_rimg(nx, 1, 1, dz, pixel_size, voltage, cs, ampcont, b_factor, dza, azz, sign)
	dict       = ctf2.to_dict()
	dz         = dict["defocus"]
	cs         = dict["cs"]
	voltage    = dict["voltage"]
	pixel_size2= dict["apix"]
	b_factor   = dict["bfactor"]
	ampcont    = dict["ampcont"]
	dza        = dict["dfdiff"]
	azz        = dict["dfang"]
	if(pixel_size != pixel_size2):
		sp_global_def.ERROR("CTFs have different pixel sizes, pixel size from the first one used", "compare_ctfs", 0)
	cim2 = EMAN2_cppwrap.Util.ctf_rimg(nx, 1, 1, dz, pixel_size, voltage, cs, ampcont, b_factor, dza, azz, sign)
	for i in range(nx//2,nx):
		if(cim1.get_value_at(i)*cim2.get_value_at(i) < 0.0):
			limi = i-nx//2
			break
	return limi, pixel_size*nx/limi
	

###----D-----------------------		










def defocus_get(fnam_roo, volt=300, Pixel_size=1, Cs=2, wgh=.1, f_start=0, f_stop=-1, docf="a", skip="#", round_off=1, nr1=3, nr2=6):
	"""
	
		1.Estimating envelope function and baseline noise using constrained simplex method
		  so as to extract CTF imprints from 1D power spectrum
		2. Based one extracted ctf imprints, perform exhaustive defocus searching to get 
		   defocus which matches the extracted CTF imprints 
	"""

	pass#IMPORTIMPORTIMPORT from math 	import sqrt, atan
	pass#IMPORTIMPORTIMPORT from sp_utilities 	import read_text_row
	roo     = []
	res     = []
	if(docf == "a"):
		TMP_roo = sp_utilities.read_text_row(fnam_roo, "a", skip)
		for i in range(len(TMP_roo)):	roo.append(TMP_roo[i][1])
	else:
		TMP_roo=sp_utilities.read_text_row(fnam_roo,"s",";")
		for i in range(len(TMP_roo)):	roo.append(TMP_roo[i][2])
	Res_roo = []
	Res_TE  = []	
	if f_start == 0 : 	i_start=0
	else: 			i_start=int(Pixel_size*2.*len(roo)/f_start)
	if f_stop <= i_start : 	i_stop=len(roo)
	else: 			i_stop=int(Pixel_size*2.*len(roo)/f_stop)

	TE  = defocus_env_baseline_fit(roo, i_start, i_stop, int(nr1), 4)
	Pn1 = defocus_env_baseline_fit(roo, i_start, i_stop, int(nr2), 3)
	Res_roo = []
	Res_TE  = []	
	for i in range(len(roo)):
		Res_roo.append(roo[i] - Pn1[i])
		Res_TE.append( TE[i]  - Pn1[i])
	#
	defocus=defocus_guess(Res_roo, Res_TE, volt, Cs, Pixel_size, wgh, i_start, i_stop, 2, round_off)
	del    roo
	del    TE
	del    Pn1
	del    Res_TE
	del    Res_roo	
	return defocus

def defocus_gett(roo, voltage=300.0, Pixel_size=1.0, Cs=2.0, wgh=0.1, f_start=0.0, f_stop=-1.0, round_off=1.0, nr1=3, nr2=6, parent=None):
	"""
	
		1. Estimate envelope function and baseline noise using constrained simplex method
		   so as to extract CTF imprints from 1D power spectrum
		2. Based one extracted ctf imprints, perform exhaustive defocus searching to get 
		   defocus which matches the extracted CTF imprints 
	"""
	pass#IMPORTIMPORTIMPORT from sp_utilities  import generate_ctf
	#print "CTF params:", voltage, Pixel_size, Cs, wgh, f_start, f_stop, round_off, nr1, nr2, parent
	if f_start == 0 : 	i_start = 0
	else: 			    i_start = int(Pixel_size*2.*len(roo)*f_start)
	if f_stop <= f_start : 	i_stop  = len(roo)
	else:                   i_stop  = min(len(roo), int(Pixel_size*2.*len(roo)*f_stop))

	#print "f_start, i_start, f_stop, i_stop:", f_start, i_start, f_stop, i_stop
	TE  = defocus_env_baseline_fit(roo, i_start, i_stop, int(nr1), 4)
	Pn1 = defocus_env_baseline_fit(roo, i_start, i_stop, int(nr2), 3)
	Res_roo = []
	Res_TE  = []
	for i in range(len(roo)):
		Res_roo.append(roo[i] - Pn1[i])
		Res_TE.append( TE[i]  - Pn1[i])

	defocus = defocus_guess(Res_roo, Res_TE, voltage, Cs, Pixel_size, wgh, i_start, i_stop, 2, round_off)

	nx  = int(len(Res_roo)*2)
	ctf = ctf_2(nx, sp_utilities.generate_ctf([defocus, Cs, voltage, Pixel_size, 0.0, wgh]))
	if (parent is not None):
		parent.ctf_data=[roo, Res_roo, Res_TE]
		parent.i_start = i_start
		parent.i_stop = i_stop
		pass#IMPORTIMPORTIMPORT from sp_utilities import write_text_file
		sp_utilities.write_text_file([roo, Res_roo, Res_TE, ctf], "procpw.txt")
	else:
		pass#IMPORTIMPORTIMPORT from sp_utilities import write_text_file
		sp_utilities.write_text_file([roo, Res_roo, Res_TE, ctf], "procpw.txt")
	return defocus

def defocus_get_Eudis(fnam_roo, volt=300, Pixel_size=1, Cs=2, wgh=.1, f_start=0, f_stop=-1, docf="a" ,skip="#", round_off=1, nr1=3, nr2=6):
	"""
		1. Estimating envelope function and baseline noise using constrained simplex method
		   so as to extract CTF imprints from 1D power spectrum
		2. Based one extracted ctf imprints, perform exhaustive defocus searching to get 
		   defocus which matches the extracted CTF imprints
		3. It returns Euclidean distance for defocus selection 
	"""
	pass#IMPORTIMPORTIMPORT from math 	import sqrt, atan
	pass#IMPORTIMPORTIMPORT from sp_utilities 	import read_text_row, generate_ctf
	roo     = []
	res     = []
	if docf == "a":
		TMP_roo = sp_utilities.read_text_row(fnam_roo, "a", skip)
		for i in range(len(TMP_roo)): # remove first record
			roo.append(TMP_roo[i][1])
	else:
		skip = ";"
		TMP_roo = sp_utilities.read_text_row(fnam_roo, "s", skip)
		for i in range(len(TMP_roo)): # remove first record
			roo.append(TMP_roo[i][2])
	Res_roo = []
	Res_TE  = []	
	if f_start == 0 : 	i_start=0
	else: 			i_start=int(Pixel_size*2.*len(roo)/f_start)
	if f_stop <= i_start :	i_stop=len(roo)
	else: 			i_stop=int(Pixel_size*2.*len(roo)/f_stop)	
	TE  = defocus_env_baseline_fit(roo, i_start, i_stop, int(nr1), 4)
	Pn1 = defocus_env_baseline_fit(roo, i_start, i_stop, int(nr2), 3)
	Res_roo = []
	Res_TE  = []
	for i in range(len(roo)):
		Res_roo.append( roo[i] - Pn1[i] )
		Res_TE.append(  TE[i]  - Pn1[i] )
	#
	defocus=defocus_guess(Res_roo, Res_TE, volt, Cs, Pixel_size, wgh, i_start, i_stop, 2, round_off)
	nx  = int(len(roo)*2)
	ctf = ctf_2(nx, sp_utilities.generate_ctf([defocus,Cs,voltage,Pixel_size, 0.0, wgh]))
	for i in range(len(Res_TE)):
		ctf[i]=ctf[i]*Res_TE[i]
	dis = defocus_L2_euc(ctf, Res_roo, i_start, i_stop)
	return [defocus, dis]

def defocus_L2_euc(v1,v2, ist,istp):
	pass#IMPORTIMPORTIMPORT from math import sqrt
	dis    = 0.0
	pw_sum = 0.0
	if ist == istp :	sp_global_def.ERROR("No pw2 curve is included  ", "defocus_L2_euc", 0)
	else:			tfeq = istp-ist
	for i in range(ist,istp,1):
		dis+=    (v1[i]-v2[2])**2
		pw_sum+= (v1[i])**2
	if pw_sum <= 0:		sp_global_def.ERROR("negative or zero power ", "defocus_L2_euc", 1)
	if dis    <= 0:		sp_global_def.ERROR("bad fitting, change options settings and try again  ", "defocus_L2_euc", 0)
	else:
		res = numpy.sqrt(dis)/numpy.sqrt(pw_sum)/tfeq	
		return res



































































































































def defocus_get_fast(indir, writetodoc="w", Pixel_size=1, volt=120, Cs=2, wgh=.1, round_off=100, dz_max0=50000, f_l0=30, f_h0=5, nr_1=5, nr_2=5, prefix="roo", docf="a",skip="#", micdir="no", print_screen="p"):
	"""
		Estimate defocus using user supplied 1D power spectrum area
		writetodoc="a" return the estimated defoci in a list, and write them down also in a text file
		writetodoc="l" output estimated defocus in a list
		writetodoc="w" output estimated defocus in a text file
	"""
	pass#IMPORTIMPORTIMPORT import os
	pass#IMPORTIMPORTIMPORT import types
	pass#IMPORTIMPORTIMPORT from sp_utilities import set_arb_params, get_image
	if writetodoc[0]   != "a" and writetodoc[0]   != "l" and writetodoc[0] != "a": 	writetodoc= "a"
	if print_screen[0] != "p" and print_screen[0] != "n"			     : 	print_screen = "n"
	if os.path.exists(indir) == False: 	sp_global_def.ERROR("roodir doesn't exist", "defocus_get_fast",1)
	ctf_dicts = ["defocus", "Pixel_size", "voltage", "Cs", "amp_contrast", "B_factor", "sign"] 
	flist = os.listdir(indir)
	res   = []
	f_l   = f_l0
	f_h   = f_h0
	if(f_l <= 1 and f_l> 0)	:
		 f_l = 1./f_l
		 f_h = 1./f_h
	if(f_h > f_l or f_l <= 0 or f_h <= 0): 
		f_h  = 8
		f_l  = 30
	if nr_1       <=  1 	:	nr_1      =  5.
	if nr_2       <=  1 	:	nr_2      =  5.
	if round_off  <=  0	: 	round_off =  100.
	if dz_max0    <=  1	: 	dz_max0   =  100000.
	dz_max=dz_max0
	if writetodoc[0] == "w" or writetodoc == "a":
		fdefo_nam = "defocus.txt"
		out       =  open(fdefo_nam, "w")
		out.write("#defocus: %s\n")
	defocus = 0
	ncount  = 0	
	for i, v in enumerate(flist):
		(fnam, fext) = os.path.splitext(v)
		if(fnam[0:len(prefix)] == prefix):
			ncount   += 1
			fnam_root = fnam[len(prefix):]
			nr1       = int(nr_1)
			nr2       = int(nr_2)
			istart    = int(f_l)
			istop     = int(f_h)
			fnam_roo  = os.path.join(indir, v)
			defocus = defocus_get(fnam_roo, volt, Pixel_size, Cs, wgh, istart, istop, docf, skip, round_off, nr1, nr2)
			if(defocus > dz_max):
				while(nr1 <= 7 or nr2 <= 7):
					nr1 += 1
					nr2 += 1
					defocus = defocus_get(fnam_roo, volt, Pixel_size, Cs, wgh,istart, istop, docf,skip, round_off, nr1, nr2)
					#if(print_screen[0] == "p" or print_screen[0] == "P" ): print "defocus",defocus,"Euclidean distance",dis,"starting feq",istart,"stop freq",istop,"P R E", nr1,"P R B", nr2
					if(defocus<dz_max): break
			if(defocus > dz_max):
				while(nr1 >= 2 and nr2 >= 2):
					nr1 -= 1
					nr2 -= 1
					defocus = defocus_get(fnam_roo, volt,Pixel_size, Cs, wgh, istart, istop, docf, skip, round_off, nr1, nr2)
					#if(print_sreen[0] == "p" or print_screen=="P"): print "defocus",defocus,"Euclidean distance",dis,"starting feq",istart,"stop freq",istop,"P R E", nr1,"P R B", nr2
					if(defocus < dz_max): break
			if(defocus > dz_max):
				while(istart > istop):
					nr1    =  5
					nr2    =  5
					istart -=.5
					defocus = defocus_get(fnam_roo, volt, Pixel_size, Cs, wgh, istart, istop, docf,skip, round_off, nr1, nr2)
					#if(print_screen[0] == "p" or print_screen == "P"): print "defocus",defocus,"Euclidean distance",dis,"starting feq",istart,"stop freq",istop,"P R E", nr1,"P R B", nr2
					if(defocus < dz_max): break
			if(defocus > dz_max):
				while(istart > istop):
					nr1     = 5										    	
					nr2     = 5
					istop  += 0.5
					defocus = defocus_get(fnam_roo, volt, Pixel_size, Cs, wgh, istart, istop, docf, skip, round_off, nr1, nr2)
					if(print_screen == "p" or print_screen == "P"): sp_global_def.sxprint("defocus",defocus,"Euclidean distance", inspect.dis, "starting feq", istart, "stop freq", istop,"P R E", nr1,"P R B", nr2)
					if(defocus < dz_max): 				break
			if(defocus >= dz_max): 					sp_global_def.ERROR("defocus_get_fast fails at estimating defocus", fnam, action = 0)
			sp_global_def.sxprint("", flist[i], '%5d'%(defocus)) 	# screen output, give the user a general impression about estimated defoci
			if(writetodoc[0] == "w" or writetodoc[0] != "l"):	out.write("%d\t%5d\t%s\n" % (ncount,defocus,flist[i]))
			if(writetodoc[0] == "l"):				res.append(defocus)
			if type(micdir) is bytes : 
				ctf_param = [defocus, Pixel_size, volt, Cs, wgh, 0, 1]
				mic_name  = os.path.join(micdir,""+ fnam_root+ ".hdf")
				if os.path.exists(mic_name) :
					e = sp_utilities.get_image (mic_name)
					U______set_arb_params(e, ctf_param, ctf_dicts)  # THIS IS INCORRECT< PLEASE CHANGE
					e.write_image(mic_name,0, EMAN2_cppwrap.EMUtil.ImageType.IMAGE_HDF, True)
					sp_global_def.sxprint("ctf parameters is written back into headers of ", mic_name)
				#else :  print  mic_name, " Not found"
	if(len(res) == 0 and  writetodoc == "l" ):				sp_global_def.ERROR("No input file is found, check the input directory of file prefix", indir, 1)
	else:
		if(writetodoc[0] == "a"):
			out.close()
			return res
		if(writetodoc[0] == "l"): 	return res
		if(writetodoc[0] == "w"): 	out.close()

def defocus_get_fast_MPI(indir, writetodoc="w", Pixel_size=1, volt=120, Cs=2, wgh=.1, round_off=100, dz_max0=50000, f_l0=30, f_h0=5, nr_1=5, nr_2=5, prefix_of_="roo", docf="a",skip="#",print_screen="no"):
	"""
		Estimate defocus using user defined 1D power spectrum area
		writetodoc="a" return the estimated defoci in a list, and write them down also in a text file
		writetodoc="l" output estimated defocus in a list
		writetodoc="w" output estimated defocus in a text file
	"""
	pass#IMPORTIMPORTIMPORT import os
	pass#IMPORTIMPORTIMPORT import sys
	if os.path.exists(indir) == False: 	sp_global_def.ERROR("roodir doesn't exist", "defocus_get_fast",1)
	flist = os.listdir(indir)
	for i, v in enumerate(flist):
		micname                  = os.path.join(indir,v)
		(filename, filextension) = os.path.splitext(v)
		if(filename[0:len(prefix_of_)] == prefix_of_):
			mic_name_list.append(micname)
			nima += 1
	if nima < 1: 	sp_global_def.ERROR("No  is found, check either directory or prefix of s is correctly given","pw2sp",1)
	
	sys.argv       = mpi.mpi_init(len(sys.argv),sys.argv)
	number_of_proc = mpi.mpi_comm_size(mpi.MPI_COMM_WORLD)
	myid           = mpi.mpi_comm_rank(mpi.MPI_COMM_WORLD)
	#  chose a random node as a main one...
	main_node      = 0
	if(myid == 0): main_node = random.randint(0,number_of_proc-1)
	main_node      = mpi.mpi_bcast(main_node, 1, mpi.MPI_INT, 0, mpi.MPI_COMM_WORLD)

	if(myid == main_node):
		if os.path.exists(outdir):  os.system('rm -rf '+outdir)
		os.mkdir(outdir)
	if(number_of_proc <= nima ):	nimage_per_node = nima/number_of_proc
	else: 				nimage_per_node = 1 
	image_start    = myid * nimage_per_node
	if(myid == number_of_proc-1):  image_end = nima
	else:                          image_end = image_start + nimage_per_node
	
	if writetodoc[0]   != "a" and writetodoc[0]   != "l" and writetodoc[0] != "a": 	writetodoc= "a"
	if print_screen[0] != "p" and print_screen[0] != "n"			     : 	print_screen = "n"
	res   = []
	f_l   = f_l0
	f_h   = f_h0
	if(f_l <= 1 and f_l> 0)	:
		 f_l = 1./f_l
		 f_h = 1./f_h
	if(f_h > f_l or f_l <= 0 or f_h <= 0): 
		f_h  = 8
		f_l  = 30
	if nr_1       <=  1 	:	nr_1      =  5.
	if nr_2       <=  1 	:	nr_2      =  5.
	if round_off  <=  0	: 	round_off =  100.
	if dz_max0    <=  1	: 	dz_max0   =  100000.
	dz_max = dz_max0
	if writetodoc[0] == "w" or writetodoc == "a":
		fdefo_nam = "defocus.txt"
		out       =  open(fdefo_nam, "w")
		out.write("#defocus: %s\n")
	defocus = 0
	ncount  = 0
	nr1	= int(nr_1)
	nr2	= int(nr_2)
	istart	= int(f_l )
	istop	= int(f_h )
	for i in range(image_start,image_end):
		filename=mic_name_list[i] 
		sp_global_def.sxprint('%-15s%-30s'%("s # ",filename))
		(f_nam, filextension) = os.path.splitext(filename)
		fnam_roo     = "particle_"+f_nam[len(prefix_of_)+len(indir)+2:]+filextension	
#	for i, v in enumerate(flist):
		ncount   += 1
		defocus = defocus_get(fnam_roo, volt, Pixel_size, Cs, wgh, istart, istop, docf, skip, round_off, nr1, nr2)
		if(defocus > dz_max):
			while(nr1 <= 7 or nr2 <= 7):
				nr1 += 1
				nr2 += 1
				defocus = defocus_get(fnam_roo, volt, Pixel_size, Cs, wgh,istart, istop, docf,skip, round_off, nr1, nr2)
				if(print_screen[0] == "p" or print_screen[0] == "P" ): sp_global_def.sxprint("defocus",defocus,"Euclidean distance",inspect.dis,"starting feq",istart,"stop freq",istop,"P R E", nr1,"P R B", nr2)
				if(defocus<dz_max): break
		if(defocus > dz_max):
			while(nr1 >= 2 and nr2 >= 2):
				nr1 -= 1
				nr2 -= 1
				defocus = defocus_get(fnam_roo, volt,Pixel_size, Cs, wgh, istart, istop, docf, skip, round_off, nr1, nr2)
				if(print_sreen[0] == "p" or print_screen=="P"): sp_global_def.sxprint("defocus",defocus,"Euclidean distance",inspect.dis,"starting feq",istart,"stop freq",istop,"P R E", nr1,"P R B", nr2)
				if(defocus < dz_max): break
		if(defocus > dz_max):
			while(istart > istop):
				nr1    =  5
				nr2    =  5
				istart -=.5
				defocus = defocus_get(fnam_roo, volt, Pixel_size, Cs, wgh, istart, istop, docf,skip, round_off, nr1, nr2)
				if(print_screen[0] == "p" or print_screen == "P"): sp_global_def.sxprint("defocus",defocus,"Euclidean distance",inspect.dis,"starting feq",istart,"stop freq",istop,"P R E", nr1,"P R B", nr2)
				if(defocus < dz_max): break
		if(defocus > dz_max):
			while(istart > istop):
				nr1     = 5										    	
				nr2     = 5
				istop  += 0.5
				defocus = defocus_get(fnam_roo, volt, Pixel_size, Cs, wgh, istart, istop, docf, skip, round_off, nr1, nr2)
				if(print_screen == "p" or print_screen == "P"): sp_global_def.sxprint("defocus",defocus,"Euclidean distance", inspect.dis, "starting feq", istart, "stop freq", istop,"P R E", nr1,"P R B", nr2)
				if(defocus < dz_max): 				break
		if(defocus >= dz_max): 					sp_global_def.ERROR("defocus_get_fast fails at estimating defocus", fnam, action = 0)
		sp_global_def.sxprint("", flist[i], '%10.3g'(defocus)) 	# screen output, give the user a general impression about estimated defoci
		if(writetodoc[0] == "w" or writetodoc[0] != "l"):	out.write("%d\t%f\t%s\n" % (ncount,defocus,flist[i]))
		if(writetodoc[0] == "l"):				res.append(defocus)
	if(len(res) == 0 and  writetodoc == "l" ):				sp_global_def.ERROR("No input file is found, check the input directory of file prefix", indir, 1)
	else:
		if writetodoc[0] == "a":
			out.close()
			return res
	if(writetodoc[0] == "l"): 	return res
	if(writetodoc[0] == "w"): 	out.close()

def defocus_get_slow(indir, writetodoc="w", Pixel_size=1, volt=120, Cs=2, wgh=.1, round_off=100, dz_max0=50000, f_l0=30, f_h0=5, prefix="roo", docf="s", skip=";",micdir="", print_screen="p"):
	"""
		Estimate defocus using user provided 1D power spectrum
		mode=1 return the estimated defoci in a list, and writes them down also in a text file
		mode=2 output estimated defocus in a list
		mode=3 output estimated defocus in a text file
		This is a slow version, more accurate than no s version
	"""
	pass#IMPORTIMPORTIMPORT from sp_morphology import defocus_get_Eudis
	pass#IMPORTIMPORTIMPORT import os
	if writetodoc[0]   != "a" and writetodoc[0]   != "l" and writetodoc[0] != "a" : writetodoc   = "a"
	if print_screen[0] != "p" and print_screen[0] != "n": 				print_screen = "n" 
	if os.path.exists(indir) == False: 	sp_global_def.ERROR("roodir doesn't exist", "defocus_get_slow",1)
	flist=os.listdir(indir)
	res  = []
	f_l  = f_l0
	f_h  = f_h0
	if f_l <= 1 and f_l > 0:
		 f_l = 1./f_l
		 f_h = 1./f_h
	if f_h > f_l or f_l <= 0 or f_h <= 0: 
		f_h=8.  # angstrom
		f_l=30. # angstrom 
	if round_off <= 0: 	round_off = 100.
	if dz_max0   <= 1: 	dz_max    = 100000.
	dz_max = dz_max0
	if( writetodoc[0] == "w" or writetodoc == "a" ):
		fdefo_nam = "defocus.txt"
		out = open(fdefo_nam, "w")
		out.write("#Coordinates: %s\n")
	ncount = 0	
	for i, v in enumerate(flist):
		(fnam, fext) = os.path.splitext(v)
		if(fnam[0:len(prefix)] == prefix):
			istart   = int(f_l)
			istop    = int(f_h)
			fnam_roo = os.path.join(indir,v)
			Mdis     = 1.e22
			defo     = 0.0
			for nr1 in range(2,7,1):
				for nr2 in range(2,7,1):
					[defocus, dis]     = defocus_get_Eudis(fnam_roo, volt, Pixel_size, Cs, wgh, istart, istop, docf, skip, round_off, nr1, nr2)
					if(print_screen[0]=="p"): sp_global_def.sxprint("defocus",defocus,"Euclidean distance",dis,"starting feq",istart,"stop freq",istop,"P R E", nr1,"P R B", nr2)
					if(Mdis > dis):
						defo = defocus
						Mdis = dis
			if(defo > dz_max):
				istart-= 1.
				for nr1 in range(3,5,1):
					for nr2 in range(2,4,1):
						[defocus, dis] = defocus_get_Eudis(fnam_roo, volt, Pixel_size, Cs, wgh, istart, istop, docf, skip, round_off, nr1, nr2)
						if(Mdis>dis):
							defo = defocus
							Mdis = dis
			if(defo >= dz_max): 	sp_global_def.ERROR("defo_get_s fails at estimating defocus from ", fnam, 0)
			else:				sp_global_def.sxprint("", flist[i], defo) # screen output, give the user a general impression about estimated defoci		
			if writetodoc    == "w" or writetodoc[0] == "a":out.write("%d\t%f\t%s\n" % (ncount, defo, fdefo_nam))
			if writetodoc[0] == "l" : 	res.append(defo)
	if  len(res) == 0 and writetodoc == "l" :  sp_global_def.ERROR("No input file, check the input directory", indir, 1)
	else:
		if writetodoc[0] == "a":
			out.close()
			return res
		if writetodoc[0] == "l": 	return res
		if writetodoc[0] == "w":	out.close()

def flcc(t, e):
	"""
		Fast local cross correlation function 
		See Alan Roseman's paper in Ultramicroscopy
	"""
	pass#IMPORTIMPORTIMPORT from sp_utilities import model_blank
	pass#IMPORTIMPORTIMPORT from sp_fundamentals import ccf
	tmp        = EMAN2_cppwrap.EMData()
	mic_avg_sq = EMAN2_cppwrap.EMData()
	mic_sq     = EMAN2_cppwrap.EMData()
	mask       = sp_utilities.model_blank(t.get_xsize(), t.get_ysize(), 1)	
	mask       +=1. 
	[mean_t, sigma_t, imin_t, imax_t] = EMAN2_cppwrap.Util.infomask(t,None,False)
	nx         = e.get_xsize()
	ny         = e.get_ysize()		
	n_pixelt   = t.get_xsize()*t.get_ysize()  # get total pixels in template   
	n_pixele   = nx*ny  # get total pixels in mic
	t          = (t-mean_t)/sigma_t # normalize the template such that the average of template is zero.
	t_pad      = EMAN2_cppwrap.Util.pad(t,    nx, ny, 1, {"background":0}, 0, 0, 0)
	m_pad      = EMAN2_cppwrap.Util.pad(mask, nx, ny, 1, {"background":0}, 0, 0, 0) # create a mask (blank, value=1 )file and pad to size of mic   	 	
	tmp        = sp_fundamentals.ccf(e, m_pad)/n_pixele # calculate the local average
	mic_avg_sq = tmp*tmp    # calculate average square
	tmp        = e*e
	mic_sq     = sp_fundamentals.ccf(tmp,m_pad)/n_pixelt 	  # calculate the average of squared mic	       
	tmp        = mic_sq-mic_avg_sq*n_pixelt   #  
	mic_var    = tmp.get_pow(.5)              # Calculate the local variance of the image 
	cc_map     = sp_fundamentals.ccf(e,t_pad)
	cc_map    /= (mic_var*n_pixelt) # Normalize the cross correlation map 
	return cc_map

##-----------------------------img formation parameters related functions---------------------------------



























def imf_B_factor_get(res_N, x, ctf_params):
	pass#IMPORTIMPORTIMPORT from scipy.optimize import fmin
	nx    = len(res_N)*2
	ctf   = ctf_1d(nx, ctf_params)
	p     = [1,1]
	xopt  = scipy.optimize.fmin(residuals_B1, p, (res_N,x))
	p     = xopt
	xopt1 = scipy.optimize.fmin(residuals_B2, p, (res_N,ctf[1][0:nx-1], x))
	sp_global_def.sxprint(xopt)
	return xopt

def imf_residuals_B1(p,y,x):
	"""
		Give the initial guess of B-factor
	"""
	pass#IMPORTIMPORTIMPORT from numpy import exp
	C,B = p
	err = 0.0
	for i in range(len(y)):
		err+= abs(y[i] - C*numpy.exp(-B*x[i]*x[i]))  # should be 4*B
	return err

def imf_residuals_B2(p,y,ctf,x):
	"""
		fit B-factor in case of considering CTF effect
	""" 
	pass#IMPORTIMPORTIMPORT from numpy import exp
	C,B = p
	err = 0.0
	for i in range(len(y)):
		err+= abs(y[i] - ctf[i]*C*numpy.exp(-B*x[i]*x[i]))  # should be 4*B
	return err

def imf_params_get(fstrN, fstrP, ctf_params, pu, nrank, q, lowf=0.01):
	"""
		Extract image formation parameters using optimization method
		Output params: 1. freq; 2.Pn1; 3.B factor.4. C; 5. C*Pu; 6. Pn2
	"""
	params = []
	w      = []
	pw_N   = get_1dpw_list(fstrN)
	pw_P   = get_1dpw_list(fstrP)
	t_N    = imf_params_cl1(pw_N,nrank,3,ctf_params[0])
	t_P    = imf_params_cl1(pw_P,nrank,3,ctf_params[0])
	res_N  = []
	res_P  = []
	for i in range(len(t_N[0])):
		res_N.append(t_N[2][i] - t_N[1][i])
		res_P.append(t_P[2][i] - t_N[1][i])
	params.append(t_N[0]) # freq
	params.append(t_N[1]) # baseline
#	params.append(t_N[1])
	parm1  = imf_B_factor_get(res_N,t_N[0],ctf_params)
	params.append(parm1[1])
	n_lowf = lowf*ctf_params[0]*len(res_P)*2
	n_lowf = int(n_lowf)
	for i in range(len(res_P)):
		if(i <= n_lowf): w.append(0.)
		else:            w.append(1.)
	parm2 = imf_fit_pu(res_P,t_N[0],ctf_params,pu,parm1[0],parm1[1],q,w)
	params.append(parm2[1])
	params.append(parm2[0])
	for i in range(len(res_N)):
		res_N[i] *= q
	params.append(res_N)
	return params

def imf_fit_pu(res_P, x, ctf_params, pu, C, B, q, w):
	pass#IMPORTIMPORTIMPORT from scipy.optimize import fmin
	res   = []
	nx    = len(res_P)*2
	ctf   = ctf_1d(nx, ctf_params)
	for i in range(len(pu)):
		res_P[i] = res_P[i]-q*C*ctf[1][i]*w[i]
		pu[i]   *= ctf[1][i]
	p     = [1]
	xopt  = scipy.optimize.fmin(residuals_pu,p,(res_P,pu,x))
	res.append(pu)
	res.append(xopt[0])
	return res

def imf_residuals_pu(p,y,pu,x):
	"""
		fit B-factor in case of considering CTF effect
	""" 
	pass#IMPORTIMPORTIMPORT from numpy import exp
	C   = p
	err = 0.0
	for i in range(len(y)):
		err+= abs(y[i] - C*pu[i])
	return err

def residuals_simplex(args, data):
	err      = 0.0
	for i in range(len(data[0])):  err -= (data[0][i] - (args[0] + (args[1]/(data[1][i]/args[2]+1.0)**2)))**2
	return err

def residuals_lsq(p,y,x):
	c1,c2,c3 = p
	err	 = []
	for i in range(len(y)):
		err.append(abs(y[i] - c1-c2/(x[i]+c3)**2))
	return err

def residuals_lsq_peak(p,y,x,c):
	pass#IMPORTIMPORTIMPORT from numpy import exp
	d1,d2,d3 = p
	c1,c2,c3 = c
	err	 = []
	for i in range(len(y)):
		tmp1 = numpy.exp(-(x[i] - d2)**2/d3)
		tmp2 = numpy.exp(c1)*numpy.exp(c2/(x[i] + c3)**2)
		err.append(abs(y[i] - tmp2 - d1*tmp1))
	return err

def residual_1dpw2(list_1dpw2, polynomial_rankB = 2, Pixel_size = 1, cut_off = 0):
	"""
		calculate signal residual from 1D rotationally averaged power spectra 
	"""
	background = []
	freq       = []
	out = EMAN2_cppwrap.Util.pw_extract(list_1dpw2[0:cut_off + 1], polynomial_rankB, 3, Pixel_size )
	for i in range(len(list_1dpw2)):
		j = i*2
		k = i*2+1
		if i <= cut_off:
			res.append(list_1dpw2[i]-background[i])
			freq.append(i/(2*Pixel_size*len(list_1dpw2)))
		else : 
			res.append(0.0)
			freq.append(i/(2*Pixel_size*len(list_1dpw2)))
	return res, freq






















































































































































































































def get_shrink_3dmask(nxinit, mask_file_name):
	pass#IMPORTIMPORTIMPORT from sp_utilities import get_im
	pass#IMPORTIMPORTIMPORT from sp_fundamentals import resample
	pass#IMPORTIMPORTIMPORT from sp_morphology   import binarize
	mask3d = sp_utilities.get_im(mask_file_name)
	nx2 = nxinit
	nx1 = mask3d.get_xsize()
	if nx1 == nx2:
		return mask3d
	else:
		shrinkage = float(nx2)/nx1
		mask3d    = binarize(sp_fundamentals.resample(mask3d,shrinkage),0.5) #added 0.5 here to fix binarization problem
		return mask3d









































































def adaptive_mask_mass(vol, mass=2000, Pixel_size=3.6):
	pass#IMPORTIMPORTIMPORT from sp_utilities  import gauss_edge, model_blank
	pass#IMPORTIMPORTIMPORT from sp_morphology import binarize, threshold, dilation
	pass#IMPORTIMPORTIMPORT from sp_filter     import filt_gaussl
	nx = vol.get_xsize()
	a = sp_filter.filt_gaussl(vol, 0.15, True)
	TH = a.find_3d_threshold(mass, Pixel_size)
	a = binarize(a,TH)
	d = a.delete_disconnected_regions(0,0,0)

	d = dilation(d, sp_utilities.model_blank(3,3,3,1.0), "BINARY")
	#d = filt_dilation(d, model_blank(3,3,3,1.0), "BINARY")
	d = sp_utilities.gauss_edge(d)
	return d
	#Util.mul_img(vol, d)
	#return threshold(vol, 0.0)






































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































def bracket_original(f, x1, h):
	c = 1.618033989 
	f1 = f(x1)
	x2 = x1 + h; f2 = f(x2)
	# Determine downhill direction and change sign of h if needed
	if f2 > f1:
		h = -h
		x2 = x1 + h; f2 = f(x2)
		# Check if minimum between x1 - h and x1 + h
		if f2 > f1: return x2,x1 - h 
	# Search loop
	for i in range (100):    
		h = c*h
		x3 = x2 + h; f3 = f(x3)
		if f3 > f2: return x1,x3
		x1 = x2; x2 = x3
		f1 = f2; f2 = f3
	sp_global_def.sxprint("Bracket did not find a mimimum")


 








































































































































def simpw2d(defocus, data2d):
	pass#IMPORTIMPORTIMPORT from sp_utilities import generate_ctf
	pass#IMPORTIMPORTIMPORT from sp_morphology import ctf_rimg
	pass#IMPORTIMPORTIMPORT from math import sqrt
	
	#             0        1     2      3     4         5             6                      7           
	#           [defocus, cs, voltage, apix, bfactor, ampcont, astigmatism_amplitude, astigmatism_angle]
	#
	#             0        1             2      3    4         5           6        7            8                     9            10
	#  data2d = [nx, experimental_pw, defocus, Cs, voltage, Pixel_size, bfactor, ampcont, astigmatism_amplitude, astigmatism_angle, mask]
	
	defocust = max(min(defocus, 6.0), 0.01)
	data2d[7] = max(min(data2d[7],99.0), 1.0)
	ct = ctf_rimg(data2d[0], sp_utilities.generate_ctf([defocust, data2d[3], data2d[4], data2d[5], data2d[6], data2d[7], data2d[8], data2d[9]]), sign=0, ny=data2d[0])
	q2 = ct.cmp("dot", ct, dict(negative = 0, mask = data2d[10], normalize = 0))#Util.infomask(ct*ct, data2d[10], True)[0]
	q1 = ct.cmp("dot", data2d[1], dict(negative = 0, mask = data2d[10], normalize = 0))
	"""Multiline Comment18"""
	#MULTILINEMULTILINEMULTILINE 18
	#MULTILINEMULTILINEMULTILINE 18
	#MULTILINEMULTILINEMULTILINE 18
	#MULTILINEMULTILINEMULTILINE 18
	#MULTILINEMULTILINEMULTILINE 18
	return  -q1/q2


def simpw1dc(defocus, data):
	pass#IMPORTIMPORTIMPORT import numpy as np
	pass#IMPORTIMPORTIMPORT from sp_morphology import ctf_2
	pass#IMPORTIMPORTIMPORT from sp_utilities import generate_ctf
	
	#[defocus, cs, voltage, apix, bfactor, ampcont, astigmatism_amplitude, astigmatism_angle]
	#  data = [subpw[i_start:i_stop], envelope[i_start:i_stop], nx, defocus, Cs, voltage, Pixel_size, ampcont, i_start, i_stop]
	# data[1] - envelope
	ct = data[1]*numpy.array( ctf_2(data[2], sp_utilities.generate_ctf([defocus, data[4], data[5], data[6], 0.0, data[7], 0.0, 0.0]))[data[8]:data[9]], numpy.float32)
	sp_global_def.sxprint(" 1d  ",sum(data[0]*ct),numpy.linalg.norm(ct,2))
	return  2.0-sum(data[0]*ct)/numpy.linalg.norm(ct,2),ctf_2(data[2], sp_utilities.generate_ctf([defocus, data[4], data[5], data[6], 0.0, data[7], 0.0, 0.0]))

def simpw2dc(defocus, data2d):
	pass#IMPORTIMPORTIMPORT from sp_utilities import generate_ctf
	pass#IMPORTIMPORTIMPORT from sp_morphology import ctf2_rimg
	pass#IMPORTIMPORTIMPORT from math import sqrt
	
	#             0        1     2      3     4         5             6                      7           
	#           [defocus, cs, voltage, apix, bfactor, ampcont, astigmatism_amplitude, astigmatism_angle]
	#
	#             0        1             2      3    4         5           6        7            8                     9            10
	#  data2d = [nx, experimental_pw, defocus, Cs, voltage, Pixel_size, bfactor, ampcont, astigmatism_amplitude, astigmatism_angle, mask]
	

	ct = ctf2_rimg(data2d[0], sp_utilities.generate_ctf([defocus, data2d[3], data2d[4], data2d[5], data2d[6], data2d[7], data2d[8], data2d[9]]), ny=data2d[0])
	pass#IMPORTIMPORTIMPORT from sp_utilities import info
	q1 = ct.cmp("dot", data2d[1], dict(negative = 0, mask = data2d[10], normalize = 0))
	q2 = numpy.sqrt(ct.cmp("dot", ct, dict(negative = 0, mask = data2d[10], normalize = 0)))
	"""Multiline Comment19"""
	#MULTILINEMULTILINEMULTILINE 19
	#MULTILINEMULTILINEMULTILINE 19
	#MULTILINEMULTILINEMULTILINE 19
	sp_global_def.sxprint(" 2d  ",q1,q2)
	return  2.0-q1/q2,ct









































































































































































































































































































































































































































































































































































































def simpw1d_crf(defocus, data):
	pass#IMPORTIMPORTIMPORT from sp_morphology import ctf_1d
	pass#IMPORTIMPORTIMPORT import numpy as np
	pass#IMPORTIMPORTIMPORT from sp_utilities import generate_ctf
	
	#[defocus, Cs, volt, Pixel_size, 0.0, ampcont]
	# data[1] - envelope
	ct = data[1]*numpy.array( ctf_1d(data[2], sp_utilities.generate_ctf([defocus, data[4], data[5], data[6], 0.0, data[7], 0.0, 0.0]))[data[8]:data[9]], numpy.float32)
	return  2.0-sum(data[0]*ct)/numpy.linalg.norm(ct,2)

def linregnp(y):
	pass#IMPORTIMPORTIMPORT import numpy as np
	ny = len(y)
	ff = type(y[0])
	x = numpy.array(list(range(ny)), ff)
	return  numpy.linalg.lstsq(numpy.vstack([x, numpy.ones(ny,ff)]).T, y)

def defocusgett_crf(roo, nx, voltage=300.0, Pixel_size=1.0, Cs=2.0, ampcont=0.1, f_start=-1.0, f_stop=-1.0, round_off=1.0, nr1=3, nr2=6, parent=None, DEBug=False):
	#(roo, nx, voltage=300.0, Pixel_size=1.0, Cs=2.0, ampcont=0.1, f_start=-1.0, f_stop=-1.0, round_off=1.0, nr1=3, nr2=6, parent=None):
	"""
	
		1. Estimate envelope function and baseline noise using constrained simplex method
		   so as to extract CTF imprints from 1D power spectrum
		2. Based one extracted ctf imprints, perform exhaustive defocus searching to get 
		   defocus which matches the extracted CTF imprints 
	"""
	pass#IMPORTIMPORTIMPORT from sp_utilities  import generate_ctf
	pass#IMPORTIMPORTIMPORT from sp_morphology import bracket_def, goldsearch_astigmatism, ctflimit, simpw1d_crf, ctf_1d
	pass#IMPORTIMPORTIMPORT import numpy as np


	#print "CTF params:", voltage, Pixel_size, Cs, wgh, f_start, f_stop, round_off, nr1, nr2, parent

	if f_start == 0 : 	    i_start = 0
	else: 			        i_start = int(Pixel_size*nx*f_start+0.5)
	if f_stop <= f_start :
		i_stop  = len(roo)
		adjust_fstop = True
	else:
		i_stop  = min(len(roo), int(Pixel_size*nx*f_stop+0.5))
		adjust_fstop = False

	nroo = len(roo)
	"""Multiline Comment26"""
	#MULTILINEMULTILINEMULTILINE 26
	#MULTILINEMULTILINEMULTILINE 26
	#MULTILINEMULTILINEMULTILINE 26
	#MULTILINEMULTILINEMULTILINE 26
	#MULTILINEMULTILINEMULTILINE 26
	#print "f_start, i_start, f_stop, i_stop:", f_start, i_start, f_stop, i_stop-1
	#  THERE IS NO NEED FOR BASELINE!!!
	#write_text_file([roo,baseline,subpw],"dbg.txt")
	#print "IN defocusgett  ",np.min(subpw),np.max(subpw)
	envelope = numpy.array([1.0]*i_stop*2, numpy.float32) # movingaverage(  abs( np.array(roo, np.float32) )   , nroo//4, 3)
	#write_text_file([roo,baseline,subpw],"dbgt.txt")

	#print "IN defocusgett  ",np.min(subpw),np.max(subpw),np.min(envelope)
	#envelope = np.ones(nroo, np.float32)
	defocus = 0.0
	data = [roo[i_start:i_stop], envelope[i_start:i_stop], nx, defocus, Cs, voltage, Pixel_size, ampcont, i_start, i_stop]
	#for i in xrange(nroo):
	#	print  i,"   ",roo[i],"   ",baseline[i],"   ",subpw[i],"   ",envelope[i]
	h = 0.1
	#def1, def2 = bracket(simpw1d, data, h)
	#if DEBug:  print "first bracket ",def1, def2,simpw1d(def1, data),simpw1d(def2, data)
	#def1=0.1
	ndefs = 18
	defound = []
	for  idef in range(ndefs):
		def1 = (idef+1)*0.5
		def1, def2 = bracket_def(simpw1d_crf, data, def1, h)
		#if DEBug:  print "second bracket ",idef,def1, def2,simpw1d_crf(def1, data),simpw1d(def2, data),h
		def1, val2 = goldsearch_astigmatism(simpw1d_crf, data, def1, def2, tol=1.0e-3)
		#if DEBug:  print "golden ",idef,def1, val2,simpw1d_crf(def1, data)
		if def1>0.0:  defound.append([val2,def1])
	defound.sort()
	del defound[3:]
	if DEBug:  sp_global_def.sxprint(" BEST DEF CANDIDATES",defound)
	if adjust_fstop:
		pass#IMPORTIMPORTIMPORT from sp_morphology import ctflimit
		newstop,fnewstop = ctflimit(nx, defound[0][1], Cs, voltage, Pixel_size)
		if DEBug:  
			sp_global_def.sxprint("newstop  ",int(newstop),fnewstop,i_stop,newstop,nx, defound[0][1])
		if( newstop != i_stop):
			i_stop = newstop
			data = [roo[i_start:i_stop], envelope[i_start:i_stop], nx, defocus, Cs, voltage, Pixel_size, ampcont, i_start, i_stop]
			h = 0.05
			for idef in range(3):
				def1, def2 = bracket_def(simpw1d_crf, data, defound[idef][1], h)
				if DEBug:  sp_global_def.sxprint(" adjusted def ",def1,def2)
				def1, val2 = goldsearch_astigmatism(simpw1d_crf, data, def1, def2, tol=1.0e-3)
				if DEBug:  sp_global_def.sxprint("adjusted golden ",def1, val2,simpw1d_crf(def1, data))
				if def1>0.0:  defound[idef] = [val2,def1]
			defound.sort()
	def1 = defound[0][1]
	if DEBug: sp_global_def.sxprint(" ultimate defocus",def1,defound)

	#defocus = defocus_guessn(Res_roo, voltage, Cs, Pixel_size, ampcont, i_start, i_stop, 2, round_off)
	#print simpw1d(def1, data),simpw1d(4.372, data)
	"""Multiline Comment27"""
	#MULTILINEMULTILINEMULTILINE 27
	#MULTILINEMULTILINEMULTILINE 27
	#MULTILINEMULTILINEMULTILINE 27
	#MULTILINEMULTILINEMULTILINE 27
	#MULTILINEMULTILINEMULTILINE 27
	if DEBug and False:
		qm = 1.e23
		toto = []
		for i in range(1000,100000,5):
			dc = float(i)/10000.0
			qt = simpw1d_crf(dc, data)
			toto.append([dc,qt])
			if(qt<qm):
				qm=qt
				defi = dc
		sp_utilities.write_text_row(toto,"toto1.txt")
		sp_global_def.sxprint(" >>>>>>>>>  ",defi,simpw1d_crf(defi, data))#,generate_ctf([defi, Cs, voltage, Pixel_size, 0.0, ampcont])
		#def1 = defi
	#exit()
	ctf1d = ctf_1d(nx, sp_utilities.generate_ctf([def1, Cs, voltage, Pixel_size, 0.0, ampcont]))

	return def1, ctf1d, None, envelope, i_start, i_stop

def envelopegett_crf(defold, roo, nx, voltage=300.0, Pixel_size=1.0, Cs=2.0, ampcont=0.1, f_start=-1.0, f_stop=-1.0, round_off=1.0, nr1=3, nr2=6, parent=None, DEBug=False):
	"""
	
		
	"""
	pass#IMPORTIMPORTIMPORT from sp_utilities  import generate_ctf
	pass#IMPORTIMPORTIMPORT from sp_morphology import ctf_1d
	pass#IMPORTIMPORTIMPORT import numpy as np


	#print "CTF params:", voltage, Pixel_size, Cs, wgh, f_start, f_stop, round_off, nr1, nr2, parent

	if f_start == 0 : 	    i_start = 0
	else: 			        i_start = int(Pixel_size*nx*f_start+0.5)
	if f_stop <= f_start :
		i_stop  = len(roo)
		adjust_fstop = True
	else:
		i_stop  = min(len(roo), int(Pixel_size*nx*f_stop+0.5))
		adjust_fstop = False

	nroo = len(roo)
	"""Multiline Comment28"""
	#MULTILINEMULTILINEMULTILINE 28
	#MULTILINEMULTILINEMULTILINE 28
	#MULTILINEMULTILINEMULTILINE 28
	#MULTILINEMULTILINEMULTILINE 28
	#MULTILINEMULTILINEMULTILINE 28
	#print "f_start, i_start, f_stop, i_stop:", f_start, i_start, f_stop, i_stop-1
	#  THERE IS NO NEED FOR BASELINE!!!
	#write_text_file([roo,baseline,subpw],"dbg.txt")
	#print "IN defocusgett  ",np.min(subpw),np.max(subpw)
	envelope = []#movingaverage(  abs( np.array(roo, np.float32) )   , nroo//4, 3)

	if adjust_fstop:
		pass#IMPORTIMPORTIMPORT from sp_morphology import ctflimit
		newstop,fnewstop = ctflimit(nx, defold, Cs, voltage, Pixel_size)
		if DEBug:  sp_global_def.sxprint("newstop  ",int(newstop*0.7),fnewstop*0.7,i_stop)
	
	ctf1d = ctf_1d(nx, sp_utilities.generate_ctf([defold, Cs, voltage, Pixel_size, 0.0, ampcont]))

	return envelope, i_start, i_stop

def fufu(args,data):
	pass#IMPORTIMPORTIMPORT from sp_morphology import fastigmatism2
	return -fastigmatism2(args[1],[data[0], data[1], data[2], args[0], data[4], data[5], data[6], data[7], data[8]])

#  
# NOTE: 2016/03/21 Toshio Moriya
# getastcrfNOE() function does not work with the new output format of cter_mrk()
# 
def getastcrfNOE(refvol, datfilesroot, voltage=300.0, Pixel_size= 1.264, Cs = 2.0, wgh = 7.0, kboot=16, DEBug = False):
	"""

	#####################   Estimation from crossresolution   ############################


	#  mpirun -np 1 python /Volumes/pawel//ED/getastcrf.py <reference volume>  <name string of the data file>


	"""

	pass#IMPORTIMPORTIMPORT from sp_applications import MPI_start_end
	pass#IMPORTIMPORTIMPORT from sp_utilities import read_text_file, write_text_file, get_im, model_blank, model_circle, amoeba, generate_ctf
	pass#IMPORTIMPORTIMPORT from sys import exit
	pass#IMPORTIMPORTIMPORT import numpy as np
	pass#IMPORTIMPORTIMPORT import os
	pass#IMPORTIMPORTIMPORT from mpi  import mpi_comm_size, mpi_comm_rank, MPI_COMM_WORLD, mpi_barrier
	pass#IMPORTIMPORTIMPORT from sp_fundamentals import tilemic, rot_avg_table
	pass#IMPORTIMPORTIMPORT from sp_morphology import threshold, bracket_def, bracket, goldsearch_astigmatism, defocus_baseline_fit, simpw1d, movingaverage, localvariance, defocusgett, defocus_guessn, defocusget_from_crf, make_real, fastigmatism, fastigmatism1, fastigmatism2, fastigmatism3, simctf, simctf2, simctf2out, fupw,ctf2_rimg
	pass#IMPORTIMPORTIMPORT from sp_alignment import Numrinit, ringwe
	pass#IMPORTIMPORTIMPORT from sp_statistics import table_stat
	pass#IMPORTIMPORTIMPORT from sp_pixel_error import angle_ave

	myid = mpi.mpi_comm_rank(mpi.MPI_COMM_WORLD)
	ncpu = mpi.mpi_comm_size(mpi.MPI_COMM_WORLD)
	main_node = 0
	
	UseOldDef = True

	#f_start = 0.022
	#f_stop  = 0.24
	
	#lenroot = len(nameroot)


	#ll = read_text_file("lookup",-1)
	#ll = [[367], [12031], [25]]
	ll = sp_utilities.read_text_file("lookup")

	#set_start, set_end = MPI_start_end(len(ll[0]), ncpu, myid)
	#for k in xrange(3):
	#	ll[k] = map(int,ll[k])[set_start:set_end]

	volft,kb = sp_projection.prep_vol(sp_utilities.get_im(refvol))


	totresi = []
	for ifi in range(len(ll)):
		#  There should be something here that excludes sets with too few images
		#namics = datfilesroot+"%05d_%06d"%(ll[0][ifi],ll[1][ifi])
		namics = datfilesroot+"%05d"%ll[ifi]
		d = EMAN2_cppwrap.EMData.read_images( namics )
		nimi = len(d)
		nx = d[0].get_xsize()
		ny2 = nx//2
		if UseOldDef:
			defold, csi, volti, apixi, bfcti, ampconti, astampi, astangi = sp_utilities.get_ctf(d[0])
			if DEBug:  sp_global_def.sxprint(" USING OLD CTF  ",defold, csi, volti, apixi, bfcti, ampconti, astampi, astangi)

		fa  = [None]*nimi
		fb  = [None]*nimi
		fbc = [None]*nimi

		for imi in range(nimi):
			phi,theta,psi,tx,ty = sp_utilities.get_params_proj(d[imi])
			# next is test
			#psi = 0.0
			#d[imi] = prgs(volft, kb, [phi,theta,psi,-tx,-ty])
			####
			fa[imi] = sp_fundamentals.fft( d[imi] )
			#prgs(volft, kb, [phi,theta,psi,-tx,-ty]).write_image("bdb:projs",imi)
			#fb[imi] = filt_ctf(fa[imi] , generate_ctf([1.9, Cs, voltage, Pixel_size, 0.0, wgh, 0.9, 177.]),False) + fft(model_gauss_noise(2., nx,nx)) #!!!!!!!!!!fft(get_im("bdb:projs",imi))  #
	
			#  next modified for test
			fb[imi] = sp_fundamentals.fft( sp_projection.prgs(volft, kb, [phi,theta,psi,-tx,-ty]) )   #fa[imi].copy()#
			#fbc[imi] = fb[imi].conjg()  #  THIS IS WRONG PAP 06/19/2018
			# next is test
			#fa[imi] = filt_ctf(fa[imi] , generate_ctf([defold, Cs, voltage, Pixel_size, 0.0, wgh, 0.9, 77.]),False) + fft(model_gauss_noise(2., nx,nx)) 

		del d  # I do not think projections are needed anymore
		adefocus = [0.0]*kboot
		aamplitu = [0.0]*kboot
		aangle   = [0.0]*kboot
	
		if True:  #try:
			for nboot in range(kboot):
				if(nboot == 0): boot = list(range(nimi))
				else:
					pass#IMPORTIMPORTIMPORT from random import randint
					for imi in range(nimi): boot[imi] = random.randint(0,nimi-1)
	
				qs = sp_utilities.model_blank(nx,nx)
				qa = sp_utilities.model_blank(nx,nx)
				qb = sp_utilities.model_blank(nx,nx)
				crf1d = []
				for imboot in range(nimi):
					imi = boot[imboot]
		
					temp = sp_statistics.fsc(fa[imi],fb[imi])[1]
					if( len(crf1d) == 0 ): crf1d = [0.0]*len(temp)
					for k in range(len(temp)):  crf1d[k] += temp[k]
					t  = make_real( EMAN2_cppwrap.Util.muln_img(fa[imi], fbc[imi]) )
					EMAN2_cppwrap.Util.mul_scalar(t, 1.0/(float(nx)**4))
					EMAN2_cppwrap.Util.add_img(qs , t)

					EMAN2_cppwrap.Util.add_img(qa, EMAN2_cppwrap.periodogram(fa[imi]))
					EMAN2_cppwrap.Util.add_img(qb, EMAN2_cppwrap.periodogram(fb[imi]))

				for k in range(len(temp)):  crf1d[k] /= nimi
				"""Multiline Comment29"""
				#MULTILINEMULTILINEMULTILINE 29
				#MULTILINEMULTILINEMULTILINE 29
				#MULTILINEMULTILINEMULTILINE 29
				#MULTILINEMULTILINEMULTILINE 29
				#MULTILINEMULTILINEMULTILINE 29
				#MULTILINEMULTILINEMULTILINE 29
				#MULTILINEMULTILINEMULTILINE 29
	
				pass#IMPORTIMPORTIMPORT from math import sqrt
				nc = nx//2
				tqa = [0.0]*(nc+1)
				for i in range(nx):
					for j in range(nx):
						r = numpy.sqrt((i-nc)**2 + (j-nc)**2)
						ir = int(r)
						if(ir<nc):
							dr = r - ir
							qqa = qa.get_value_at(i,j)
							tqa[ir]   += (1.0-dr)*qqa
							tqa[ir+1] +=       dr*qqa
				for i in range(nc+1): tqa[i] = numpy.sqrt(max(tqa[i],0.0))

				divs = sp_utilities.model_blank(nx, nx, 1, 1.0)
				for i in range(nx):
					for j in range(nx):
						r = numpy.sqrt((i-nc)**2 + (j-nc)**2)
						ir = int(r)
						if(ir<nc):
							dr = r - ir
							divs.set_value_at(i,j,  (1.-dr)*tqa[ir] + dr*tqa[ir+1] )
				#if(nboot == 0): qs.write_image("rs1.hdf")
				EMAN2_cppwrap.Util.div_img(qs, divs)
				qs.set_value_at(ny2,ny2,1.0)
				#if(nboot == 0): write_text_file(crf1d,"crf1d.txt")
				#if(nboot == 0): qs.write_image("rs2.hdf")


				sroo = sp_fundamentals.rot_avg_table(qs)
				lenroo = len(sroo)
				#  Find a break point
				bp = 1.e23
				for i in range(1,max(3,lenroo//4)):
					if( sroo[i] <0.0 ):
						istart = max(3,i//2)
						break

				#istart = 25
				#print istart
	
				f_start = istart/(Pixel_size*nx)
				#print namics[ifi],istart,f_start


				if UseOldDef:
					envelope, istart, istop = envelopegett_crf(defold, crf1d, nx, voltage=voltage, Pixel_size=Pixel_size, Cs=Cs, ampcont=wgh, f_start=f_start, f_stop=-1.0, round_off=1.0, nr1=3, nr2=6, parent=None, DEBug=DEBug)
					defc = defold
				else:
					defc, ctf1d, baseline, envelope, istart, istop = defocusgett_crf(crf1d, nx, voltage=voltage, Pixel_size=Pixel_size, Cs=Cs, ampcont=wgh, f_start=f_start, f_stop=-1.0, round_off=1.0, nr1=3, nr2=6, parent=None, DEBug=DEBug)
					if DEBug:  sp_global_def.sxprint("  RESULT ",namics,defc, istart, istop)
					if DEBug:
						freq = list(range(len(crf1d)))
						for i in range(len(crf1d)):  freq[i] = float(i)/nx/Pixel_size
						sp_utilities.write_text_file([freq, crf1d, ctf1d, envelope.tolist()],"ravg%05d.txt"%ifi)
				#mpi_barrier(MPI_COMM_WORLD)
				#   NOT USING ENVELOPE!
				"""Multiline Comment30"""
				#MULTILINEMULTILINEMULTILINE 30
				#MULTILINEMULTILINEMULTILINE 30
				#MULTILINEMULTILINEMULTILINE 30
				#MULTILINEMULTILINEMULTILINE 30
					#MULTILINEMULTILINEMULTILINE 30
						#MULTILINEMULTILINEMULTILINE 30
						#MULTILINEMULTILINEMULTILINE 30
						#MULTILINEMULTILINEMULTILINE 30
							#MULTILINEMULTILINEMULTILINE 30
							#MULTILINEMULTILINEMULTILINE 30
				#MULTILINEMULTILINEMULTILINE 30

				#exit()
				#istop = nx//4
				#mask = model_circle(istop-1,nx,nx)*(model_blank(nx,nx,1,1.0)-model_circle(istart,nx,nx))
				#qse = qa*envl
				#(qse*mask).write_image("rs2.hdf")
				#qse.write_image("rs3.hdf")
				##  SIMULATION
				#bang = 0.7
				#qse = ctf2_rimg(nx, generate_ctf([defc,Cs,voltage,Pixel_size,0.0,wgh, bang, 37.0]) )
				#qse.write_image("rs3.hdf")

				mask = sp_utilities.model_circle(istop-1,nx,nx)*(sp_utilities.model_blank(nx,nx,1,1.0)-sp_utilities.model_circle(istart,nx,nx))
				qse = qs #* envl
				#if(nboot == 0): (qs*mask).write_image("rs5.hdf")


				cnx = nx//2+1
				cny = cnx
				mode = "H"
				numr = sp_alignment.Numrinit(istart, istop, 1, mode)
				wr   = sp_alignment.ringwe(numr, mode)

				crefim = EMAN2_cppwrap.Util.Polar2Dm(qse, cnx, cny, numr, mode)
				EMAN2_cppwrap.Util.Frngs(crefim, numr)
				EMAN2_cppwrap.Util.Applyws(crefim, numr, wr)


				"""Multiline Comment31"""
				#MULTILINEMULTILINEMULTILINE 31
				#MULTILINEMULTILINEMULTILINE 31
				#MULTILINEMULTILINEMULTILINE 31
				#MULTILINEMULTILINEMULTILINE 31
				#MULTILINEMULTILINEMULTILINE 31
				#MULTILINEMULTILINEMULTILINE 31

				#pc = ctf2_rimg(nx,generate_ctf([defc,Cs,voltage,Pixel_size,0.0,wgh]))
				#print ccc(pc*envl, subpw, mask)

				bang = 0.0
				bamp = 0.0
				bdef = defc
				bold = 1.e23
				dstep = 0.1
				while( True):
					data = [qse, mask, nx, bamp, Cs, voltage, Pixel_size, wgh, bang]
					h = 0.05*bdef
					#print "  bdef  at the beginning of while loop   ",nboot,bdef
					amp1, amp2 = bracket_def(simctf, data, bdef*0.9, h)
					#print "bracketing of the defocus  ",nboot,amp1, amp2
					amp1, val2 = goldsearch_astigmatism(simctf, data, amp1, amp2, tol=1.0e-3)
					#print "golden defocus ",amp1, val2,simctf(amp1, data)
					#,simctf(1.1, [qse, mask, nx, 0.5, Cs, voltage, Pixel_size, wgh, 77.])
					#print " ttt ",time()-srtt
					#bdef = 1.1
					#bamp = 0.2
					#write_text_file(ctf_1d(data[2], generate_ctf([amp1, data[4], data[5], data[6], 0.0, data[7], 0.0, 0.0])), "gctf1d.txt")
					#exit()
	
					astdata = [crefim, numr, nx, bdef, Cs, voltage, Pixel_size, wgh, bang]
	
					h = 0.01
					amp1,amp2 = bracket(fastigmatism2, astdata, h)
					#print "  astigmatism bracket  ",nboot,amp1,amp2,astdata[-1]
					#print " ttt ",time()-srtt
	
					bamp, bcc = goldsearch_astigmatism(fastigmatism2, astdata, amp1, amp2, tol=1.0e-3)
	
					#junk = fastigmatism2(bamp,astdata)
					#bang = astdata[-1]
					#print " ang within the loop   ",bdef, bamp,astdata[-1],junk, fastigmatism2(0.5,[crefim, numr, nx, 1.1, Cs, voltage, Pixel_size, wgh, bang])
					#print  fastigmatism2(0.0,astdata)
					#print astdata[-1]
					#print "  golden search ",bamp,data[-1], fastigmatism2(bamp,data), fastigmatism2(0.0,data)
					#print " ttt ",time()-srtt
					#bamp = 0.5
					#bang = 277
					dama = sp_utilities.amoeba([bdef, bamp],[0.2,0.2], fufu, 1.e-4,1.e-4,500, astdata)

					if DEBug:  sp_global_def.sxprint("AMOEBA    ",nboot,dama)
					bdef = dama[0][0]
					bamp = dama[0][1]
					astdata = [crefim, numr, nx, bdef, Cs, voltage, Pixel_size, wgh, bang]
					junk = fastigmatism2(bamp, astdata)
					bang = astdata[-1]
					if DEBug:  sp_global_def.sxprint(" after amoeba ", nboot,bdef, bamp, bang)
					#  The looping here is blocked as one shot at amoeba is good enough.  To unlock it, remove - from bold.
					if(bcc < -bold): bold = bcc
					else:           break

				adefocus[nboot] = bdef
				aamplitu[nboot] = bamp
				aangle[nboot]   = bang
				#from sys import exit
				if DEBug:  sp_global_def.sxprint("this is what I found  ",nboot,bdef,bamp,bang)
				#exit()

			#print " ttt ",time()-srtt
			#from sys import exit
			#exit()
			ad1,ad2,ad3,ad4 = sp_statistics.table_stat(adefocus)
			reject = []
			thr = 3*numpy.sqrt(ad3)
			for i in range(len(adefocus)):
				if(abs(adefocus[i]-ad1)>thr):
					if DEBug:  sp_global_def.sxprint(adefocus[i],ad1,thr)
					reject.append(i)
			if(len(reject)>0):
				if DEBug:  sp_global_def.sxprint("  Number of rejects  ",namics,len(reject))
				for i in range(len(reject)-1,-1,-1):
					del adefocus[i]
					del aamplitu[i]
					del aangle[i]
			if(len(adefocus)<2):
				sp_global_def.sxprint("  After rejection of outliers too few estimated defocus values for :",namics)
			else:
				#print "adefocus",adefocus
				#print  "aamplitu",aamplitu
				#print "aangle",aangle
				ad1,ad2,ad3,ad4 = sp_statistics.table_stat(adefocus)
				bd1,bd2,bd3,bd4 = sp_statistics.table_stat(aamplitu)
				cd1,cd2 = sp_pixel_error.angle_ave([2*q for q in aangle])  # Have to use this trick as the function works for range [0,360]
				cd1/=2
				cd2/=2
				temp = 0.0
				sp_global_def.sxprint(namics,ad1, Cs, voltage, Pixel_size, temp, wgh, bd1, cd1, numpy.sqrt(max(0.0,ad2)),numpy.sqrt(max(0.0,bd2)),cd2) 
				totresi.append( [ namics, ad1, Cs, voltage, Pixel_size, temp, wgh, bd1, cd1, numpy.sqrt(max(0.0,ad2)),numpy.sqrt(max(0.0,bd2)),cd2 ])
				#if ifi == 4 : break
				"""Multiline Comment32"""
				#MULTILINEMULTILINEMULTILINE 32
					#MULTILINEMULTILINEMULTILINE 32
					#MULTILINEMULTILINEMULTILINE 32
					#MULTILINEMULTILINEMULTILINE 32
				#MULTILINEMULTILINEMULTILINE 32
				#ctf_rimg(nx,generate_ctf([ad1, Cs, voltage, Pixel_size, temp, wgh, bd1, cd1])).write_image("ctf1.hdf")
				lnsb = len(crf1d)
	
				try:		crot2 = rotavg_ctf(ctf_rimg(nx,sp_utilities.generate_ctf([ad1, Cs, voltage, Pixel_size, temp, wgh, bd1, cd1])), ad1, Cs, voltage, Pixel_size, bd1, cd1)[:lnsb]
				except:     crot2 = [0.0]*lnsb
				try:		pwrot2 = rotavg_ctf(qs, ad1, Cs, voltage, Pixel_size, bd1, cd1)[:lnsb]
				except:     pwrot2 = [0.0]*lnsb
				try:		crot1 = rotavg_ctf(ctf_rimg(nx,sp_utilities.generate_ctf([ad1, Cs, voltage, Pixel_size, temp, wgh, bd1, cd1])), ad1, Cs, voltage, Pixel_size, 0.0, 0.0)[:lnsb]
				except:     crot1 = [0.0]*lnsb
				try:		pwrot1 = rotavg_ctf(qs, ad1, Cs, voltage, Pixel_size, 0.0, 0.0)[:lnsb]
				except:     pwrot1 = [0.0]*lnsb
				freq = list(range(lnsb))
				for i in range(len(freq)):  freq[i] = float(i)/nx/Pixel_size
				#fou = "crfrot/rotinf%05d_%06d.txt"%(ll[0][ifi],ll[1][ifi])
				fou = "crfrot/rotinf%05d.txt"%(ll[ifi])
				#  #1 - rotational averages without astigmatism, #2 - with astigmatism
				sp_utilities.write_text_file([list(range(len(crot1))), freq, pwrot1, crot1, pwrot2, crot2],fou)
				cmd = "echo "+"    "+namics+"  >>  "+fou
				os.system(cmd)
		else:  #except:
			sp_global_def.sxprint(namics,"     FAILED")
	#from utilities import write_text_row
	outf = open( "partcrf/partcrf_%05d"%myid, "w")
	for i in range(len(totresi)):
		for k in range(1,len(totresi[i])): outf.write("  %12.5g"%totresi[i][k])
		outf.write("  %s\n"%totresi[i][0])
	outf.close()		

########################################################################
# end of code used for estimation of cross resolution
########################################################################


################
#
#  CTER code (07/10/2017)
#
################
# 
# 











































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































def Xdefocusgett_vpp2(qse, roo, nx, xdefc, xampcont, voltage=300.0, Pixel_size=1.0, Cs=2.0, f_start=-1.0, f_stop=-1.0, round_off=1.0, nr1=3, nr2=6, parent=None, DEBug=False):
	"""
		1. Estimate envelope function and baseline noise using constrained simplex method
		   so as to extract CTF imprints from 1D power spectrum
		2. Based one extracted ctf imprints, perform exhaustive defocus searching to get
		   defocus which matches the extracted CTF imprints
	"""
	pass#IMPORTIMPORTIMPORT from sp_utilities  import generate_ctf
	pass#IMPORTIMPORTIMPORT import numpy as np
	pass#IMPORTIMPORTIMPORT from sp_morphology import ctf_2, bracket_def, defocus_baseline_fit, ctflimit, simpw1d, goldsearch_astigmatism

	#print "CTF params:", voltage, Pixel_size, Cs, wgh, f_start, f_stop, round_off, nr1, nr2, parent

	if f_start == 0 : 	    i_start = 0
	else: 			        i_start = int(Pixel_size*nx*f_start+0.5)
	if f_stop <= f_start :
		i_stop  = len(roo)
		adjust_fstop = True
	else:
		i_stop  = min(len(roo), int(Pixel_size*nx*f_stop+0.5))
		adjust_fstop = False

	nroo = len(roo)

	if DEBug:  sp_global_def.sxprint("f_start, i_start, f_stop, i_stop:", f_start, i_start, f_stop, i_stop-1)
	#TE  = defocus_env_baseline_fit(roo, i_start, i_stop, int(nr1), 4)
	#baseline = defocus_baseline_fit(roo, i_start, i_stop, int(nr2), 3)

	baseline = defocus_baseline_fit(roo, i_start,nroo, int(nr2), 3)
	subpw = numpy.array(roo, numpy.float32) - baseline
	subpw[0] = subpw[1]
	#write_text_file([roo,baseline,subpw],"dbg.txt")
	#print "IN defocusgett  ",np.min(subpw),np.max(subpw)
	for i in range(len(subpw)):  subpw[i] = max(subpw[i],0.0)
	#print "IN defocusgett  ",np.min(subpw),np.max(subpw)
	#envelope = movingaverage(  subpw   , nroo//4, 3)
	envelope = numpy.array([1.0]*len(subpw), numpy.float32)
	#write_text_file([roo,baseline,subpw],"dbgt.txt")

	#print "IN defocusgett  ",np.min(subpw),np.max(subpw),np.min(envelope)
	#envelope = np.ones(nroo, np.float32)
	defocus = 0.0
	ampcont = 0.0
	data = [subpw[i_start:i_stop], envelope[i_start:i_stop], nx, defocus, Cs, voltage, Pixel_size, ampcont, i_start, i_stop]
	wn = 512
	pass#IMPORTIMPORTIMPORT from sp_utilities import model_circle, model_blank, amoeba
	pass#IMPORTIMPORTIMPORT from sp_alignment import Numrinit, ringwe
	mask = sp_utilities.model_circle(i_stop - 1, wn, wn) * (sp_utilities.model_blank(wn, wn, 1, 1.0) - sp_utilities.model_circle(i_start, wn, wn))
	pass#IMPORTIMPORTIMPORT from sp_fundamentals import rot_avg_table
	zizi = sp_fundamentals.rot_avg_table(qse)[i_start:i_stop]
	pass#IMPORTIMPORTIMPORT from sp_utilities import write_text_file
	dudi = subpw[i_start:i_stop]
	#print dudi.tolist()
	#print zizi

	cnx = wn // 2 + 1
	cny = cnx
	mode = "H"
	numr = sp_alignment.Numrinit(i_start, i_stop, 1, mode)
	wr = sp_alignment.ringwe(numr, mode)
	
	crefim = EMAN2_cppwrap.Util.Polar2Dm(qse*mask, cnx, cny, numr, mode)
	sp_global_def.sxprint("  CREFIM    ",EMAN2_cppwrap.Util.infomask(qse*mask,None,True),EMAN2_cppwrap.Util.infomask(crefim,None,True))
	EMAN2_cppwrap.Util.Frngs(crefim, numr)
	EMAN2_cppwrap.Util.Applyws(crefim, numr, wr)
	bdef = 0.
	baco = 0.0  #  amplitude contrast
	bamp = 0.0      #  initial astigmatism amplitude
	bang = 0.0      #  initial astigmatism angle
	astdata = [crefim, numr, wn, bdef, Cs, voltage, Pixel_size, baco, bamp, bang, mask]
	data2d = [nx, qse, bdef, Cs, voltage, Pixel_size, 0.0, baco, bamp, bang, mask]

	sp_global_def.sxprint(" i_start:i_stop",i_start,i_stop)

	qm = 1.e23
	dm = 1.e23
	dp = 1.0e23
	toto = []
	"""Multiline Comment42"""
	#MULTILINEMULTILINEMULTILINE 42
	#MULTILINEMULTILINEMULTILINE 42
		#MULTILINEMULTILINEMULTILINE 42
		#MULTILINEMULTILINEMULTILINE 42
		#MULTILINEMULTILINEMULTILINE 42
		#MULTILINEMULTILINEMULTILINE 42
			#MULTILINEMULTILINEMULTILINE 42
			#MULTILINEMULTILINEMULTILINE 42
			#MULTILINEMULTILINEMULTILINE 42
			#MULTILINEMULTILINEMULTILINE 42
			#MULTILINEMULTILINEMULTILINE 42
			#MULTILINEMULTILINEMULTILINE 42
			#MULTILINEMULTILINEMULTILINE 42
			#MULTILINEMULTILINEMULTILINE 42
			#MULTILINEMULTILINEMULTILINE 42
			#MULTILINEMULTILINEMULTILINE 42
			#MULTILINEMULTILINEMULTILINE 42
			#MULTILINEMULTILINEMULTILINE 42
			#MULTILINEMULTILINEMULTILINE 42
			#MULTILINEMULTILINEMULTILINE 42
	#MULTILINEMULTILINEMULTILINE 42
	for aa in range(0,20,4):
		a = xampcont + aa
		data[7] = float(a)
		sp_global_def.sxprint("  fdasfdsfa  ",a)
		for i in range(0,2000,200):
			dc = xdefc + float(i-1000)/10000.0
			qt = simpw1d(dc, data)
			ju1 = dc # defocus
			ju2 = float(a) # amp contrast
			ju3 = 0.0  # astigma amp
			dama = sp_utilities.amoeba([ju1,ju2,ju3], [0.005, 2.0, 0.002], fupw_vpp, 1.e-4, 1.e-4, 200, astdata)
			data2d[7] = float(a)
			zigi = simpw2d(dc, data2d)
			qma = -dama[-2]
			sp_global_def.sxprint(" amoeba  %7.2f  %7.2f  %12.6g  %12.6g  %12.6g  %7.2f  %7.2f  %7.2f "%(dc,data[7],qma,zigi,qt,dama[0][0],dama[0][1],dama[0][2]), dama)
			toto.append([dc,data[7],qt,zigi,qma])
			if(qma<dp):
				dp = qma
				dpefi = dama[0][0]
				dpmpcont = dama[0][1]
			if(zigi<dm):
				dm = zigi
				ddefi = dc
				dampcont = data[7]
			if(qt<qm):
				qm = qt
				defi = dc
				ampcont = data[7]
	if DEBug:
		pass#IMPORTIMPORTIMPORT from sp_utilities import write_text_row
		sp_utilities.write_text_row(toto,"toto1.txt")
		sp_global_def.sxprint(" repi3  ",dp,dpefi,dpmpcont)
		sp_global_def.sxprint(" resi2  ",qm,defi,ampcont)
		sp_global_def.sxprint(" resi1  ",dm,ddefi,dampcont)
		
		#print " >>>>>>>>>  ",defi,simpw1d(defi, data)#,generate_ctf([defi, Cs, voltage, Pixel_size, 0.0, ampcont])
		#def1 = defi
	#exit()
	pass#IMPORTIMPORTIMPORT from sp_morphology import ctf2_rimg, ctf_rimg, square_root
	ctf2 = ctf_rimg(nx, sp_utilities.generate_ctf([defi, Cs, voltage, Pixel_size, 0.0, ampcont]), sign=0)
	cq = ctf_1d(nx, sp_utilities.generate_ctf([defi, Cs, voltage, Pixel_size, 0.0, ampcont]), doabs = True)[20:150]
	qse.write_image("qse.hdf")
	ctf2.write_image("c1.hdf")
	ctf22 = ctf_rimg(nx, sp_utilities.generate_ctf([ddefi, Cs, voltage, Pixel_size, 0.0, dampcont]), sign=0)
	ci = ctf_1d(nx, sp_utilities.generate_ctf([ddefi, Cs, voltage, Pixel_size, 0.0, dampcont]), doabs = True)[20:150]
	dq = ctf_1d(nx, sp_utilities.generate_ctf([dpefi, Cs, voltage, Pixel_size, 0.0, dpmpcont]), doabs = True)[20:150]
	sp_utilities.write_text_file([dudi.tolist(),zizi,cq,ci,dq],"pwds.txt")
	ctf22.write_image("c2.hdf")
	"""Multiline Comment43"""
	#MULTILINEMULTILINEMULTILINE 43
	#MULTILINEMULTILINEMULTILINE 43
	#MULTILINEMULTILINEMULTILINE 43
	#MULTILINEMULTILINEMULTILINE 43
	return defi, ampcont, subpw, ctf2, baseline, envelope, i_start, i_stop


def Xdefocusgett_vpp22(qse, roo, nx, voltage=300.0, Pixel_size=1.0, Cs=2.0, f_start=-1.0, f_stop=-1.0, round_off=1.0, nr1=3, nr2=6, parent=None, DEBug=False):
	"""
		1. Estimate envelope function and baseline noise using constrained simplex method
		   so as to extract CTF imprints from 1D power spectrum
		2. Based one extracted ctf imprints, perform exhaustive defocus searching to get
		   defocus which matches the extracted CTF imprints
	"""
	pass#IMPORTIMPORTIMPORT from sp_utilities  import generate_ctf
	pass#IMPORTIMPORTIMPORT import numpy as np
	pass#IMPORTIMPORTIMPORT from sp_morphology import ctf_2, bracket_def, defocus_baseline_fit, ctflimit, simpw1d, goldsearch_astigmatism

	#print "CTF params:", voltage, Pixel_size, Cs, wgh, f_start, f_stop, round_off, nr1, nr2, parent

	if f_start == 0 : 	    i_start = 0
	else: 			        i_start = int(Pixel_size*nx*f_start+0.5)
	if f_stop <= f_start :
		i_stop  = len(roo)
		adjust_fstop = True
	else:
		i_stop  = min(len(roo), int(Pixel_size*nx*f_stop+0.5))
		adjust_fstop = False

	nroo = len(roo)

	if DEBug:  sp_global_def.sxprint("f_start, i_start, f_stop, i_stop:", f_start, i_start, f_stop, i_stop-1)
	#TE  = defocus_env_baseline_fit(roo, i_start, i_stop, int(nr1), 4)
	#baseline = defocus_baseline_fit(roo, i_start, i_stop, int(nr2), 3)

	baseline = defocus_baseline_fit(roo, i_start,nroo, int(nr2), 3)
	subpw = numpy.array(roo, numpy.float32) - baseline
	subpw[0] = subpw[1]
	#write_text_file([roo,baseline,subpw],"dbg.txt")
	#print "IN defocusgett  ",np.min(subpw),np.max(subpw)
	for i in range(len(subpw)):  subpw[i] = max(subpw[i],0.0)
	#print "IN defocusgett  ",np.min(subpw),np.max(subpw)
	#envelope = movingaverage(  subpw   , nroo//4, 3)
	envelope = numpy.array([1.0]*len(subpw), numpy.float32)
	#write_text_file([roo,baseline,subpw],"dbgt.txt")

	#print "IN defocusgett  ",np.min(subpw),np.max(subpw),np.min(envelope)
	#envelope = np.ones(nroo, np.float32)
	defocus = 0.0
	ampcont = 0.0
	data = [subpw[i_start:i_stop], envelope[i_start:i_stop], nx, defocus, Cs, voltage, Pixel_size, ampcont, i_start, i_stop]
	wn = 512
	pass#IMPORTIMPORTIMPORT from sp_utilities import model_circle, model_blank, amoeba
	pass#IMPORTIMPORTIMPORT from sp_alignment import Numrinit, ringwe
	mask = sp_utilities.model_circle(i_stop - 1, wn, wn) * (sp_utilities.model_blank(wn, wn, 1, 1.0) - sp_utilities.model_circle(i_start, wn, wn))
	pass#IMPORTIMPORTIMPORT from sp_fundamentals import rot_avg_table
	zizi = sp_fundamentals.rot_avg_table(qse)[i_start:i_stop]
	pass#IMPORTIMPORTIMPORT from sp_utilities import write_text_file
	dudi = subpw[i_start:i_stop]
	#print dudi.tolist()
	#print zizi

	cnx = wn // 2 + 1
	cny = cnx
	mode = "H"
	numr = sp_alignment.Numrinit(i_start, i_stop, 1, mode)
	wr = sp_alignment.ringwe(numr, mode)
	
	crefim = EMAN2_cppwrap.Util.Polar2Dm(qse*mask, cnx, cny, numr, mode)
	EMAN2_cppwrap.Util.Frngs(crefim, numr)
	EMAN2_cppwrap.Util.Applyws(crefim, numr, wr)
	bdef = 0.
	baco = 0.0  #  amplitude contrast
	bamp = 0.0      #  initial astigmatism amplitude
	bang = 0.0      #  initial astigmatism angle
	astdata = [crefim, numr, wn, bdef, Cs, voltage, Pixel_size, baco, bamp, bang, mask]
	data2d = [nx, qse, bdef, Cs, voltage, Pixel_size, 0.0, baco, bamp, bang, mask]

	sp_global_def.sxprint(" i_start:i_stop",i_start,i_stop)

	qm = 1.e23
	dm = 1.e23
	dp = 1.0e23
	toto = []
	"""Multiline Comment44"""
	#MULTILINEMULTILINEMULTILINE 44
	#MULTILINEMULTILINEMULTILINE 44
		#MULTILINEMULTILINEMULTILINE 44
		#MULTILINEMULTILINEMULTILINE 44
		#MULTILINEMULTILINEMULTILINE 44
		#MULTILINEMULTILINEMULTILINE 44
			#MULTILINEMULTILINEMULTILINE 44
			#MULTILINEMULTILINEMULTILINE 44
			#MULTILINEMULTILINEMULTILINE 44
			#MULTILINEMULTILINEMULTILINE 44
			#MULTILINEMULTILINEMULTILINE 44
			#MULTILINEMULTILINEMULTILINE 44
			#MULTILINEMULTILINEMULTILINE 44
			#MULTILINEMULTILINEMULTILINE 44
			#MULTILINEMULTILINEMULTILINE 44
			#MULTILINEMULTILINEMULTILINE 44
			#MULTILINEMULTILINEMULTILINE 44
			#MULTILINEMULTILINEMULTILINE 44
			#MULTILINEMULTILINEMULTILINE 44
			#MULTILINEMULTILINEMULTILINE 44
	#MULTILINEMULTILINEMULTILINE 44
	for a in range(5,96,10):
		data[7] = float(a)
		sp_global_def.sxprint("  fdasfdsfa  ",a)
		for i in range(1000,100000,5000):
			dc = float(i)/10000.0
			qt = simpw1d(dc, data)
			ju1 = dc # defocus
			ju2 = float(a) # amp contrast
			ju3 = 0.0  # astigma amp
			dama = sp_utilities.amoeba([ju1,ju2,ju3], [0.002, 0.001, 0.002], fupw_vpp, 1.e-4, 1.e-4, 1, astdata)
			data2d[7] = float(a)
			zigi = simpw2d(dc, data2d)
			qma = dama[-2]/42.
			sp_global_def.sxprint(" amoeba  %7.2f  %7.2f  %12.6g  %12.6g  %12.6g  %7.2f  %7.2f  %7.2f "%(dc,data[7],qma,zigi,qt,dama[0][0],dama[0][1],dama[0][2]), dama)
			toto.append([dc,data[7],qt,zigi,qma])
			if(qma<dp):
				dp = qma
				dpefi = dc
				dpmpcont = data[7]
			if(zigi<dm):
				dm = zigi
				ddefi = dc
				dampcont = data[7]
			if(qt<qm):
				qm = qt
				defi = dc
				ampcont = data[7]
	if DEBug:
		pass#IMPORTIMPORTIMPORT from sp_utilities import write_text_row
		sp_utilities.write_text_row(toto,"toto1.txt")
		sp_global_def.sxprint(" repi3  ",dp,dpefi,dpmpcont)
		sp_global_def.sxprint(" resi2  ",qm,defi,ampcont)
		sp_global_def.sxprint(" resi1  ",dm,ddefi,dampcont)
		
		#print " >>>>>>>>>  ",defi,simpw1d(defi, data)#,generate_ctf([defi, Cs, voltage, Pixel_size, 0.0, ampcont])
		#def1 = defi
	#exit()
	pass#IMPORTIMPORTIMPORT from sp_morphology import ctf2_rimg, ctf_rimg, square_root
	ctf2 = ctf_rimg(nx, sp_utilities.generate_ctf([defi, Cs, voltage, Pixel_size, 0.0, ampcont]), sign=0)
	cq = ctf_1d(nx, sp_utilities.generate_ctf([defi, Cs, voltage, Pixel_size, 0.0, ampcont]), doabs = True)[20:150]
	qse.write_image("qse.hdf")
	ctf2.write_image("c1.hdf")
	ctf22 = ctf_rimg(nx, sp_utilities.generate_ctf([ddefi, Cs, voltage, Pixel_size, 0.0, dampcont]), sign=0)
	ci = ctf_1d(nx, sp_utilities.generate_ctf([ddefi, Cs, voltage, Pixel_size, 0.0, dampcont]), doabs = True)[20:150]
	dq = ctf_1d(nx, sp_utilities.generate_ctf([dpefi, Cs, voltage, Pixel_size, 0.0, dpmpcont]), doabs = True)[20:150]
	sp_utilities.write_text_file([dudi.tolist(),zizi,cq,ci,dq],"pwds.txt")
	ctf22.write_image("c2.hdf")
	"""Multiline Comment45"""
	#MULTILINEMULTILINEMULTILINE 45
	#MULTILINEMULTILINEMULTILINE 45
	#MULTILINEMULTILINEMULTILINE 45
	#MULTILINEMULTILINEMULTILINE 45
	return defi, ampcont, subpw, ctf2, baseline, envelope, i_start, i_stop










