#!/usr/bin/env python

# Author: James Michael Bell, jmbell@bcm.edu, 10/18/2015

from EMAN2 import *
import numpy as np
from IPython import embed
from copy import deepcopy

def main(args):

	progname = os.path.basename(sys.argv[0])
	usage = """e2moviealigner.py [options] <ddd_movie_stack> <ddd_movie_stack> ... <ddd_movie_stack>
	
	Determines the optimal whole-frame alignment of a DDD movie. It can be
	used to determine the proper x and y translations for each movie frame and
	can perform the actual alignment according to the translations.
	
	Example: e2moviealigner.py --dark dark.hdf --gain gain.hdf movie.hdf -v 9
	"""
	
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	
	parser.add_argument("--dark",type=str,default=None,help="Perform dark image correction using the specified image file")
	parser.add_argument("--gain",type=str,default=None,help="Perform gain image normalization using the specified image file")
	parser.add_argument("--gaink2",type=str,default=None,help="Perform gain image normalization. Gatan K2 gain images are the reciprocal of DDD gain images.")
	parser.add_argument("--boxsize", type=int, help="Set the boxsize used to compute power spectra across movie frames",default=512)
	parser.add_argument("--maxshift", type=int, help="Set the maximum frame translation distance in pixels.",default=24)
	parser.add_argument("--step", type=float, help="Set step size for cross coherence calculation.",default=1.0)
	parser.add_argument("--nnorm", type=float, help="Set the norm to be used for fixing axes. Default is sqrt(2)",default=np.sqrt(2))
	parser.add_argument("--fixaxes",action="store_true",default=False,help="Tries to identify bad pixels and fill them in with sane values instead")
	parser.add_argument("--fixbadlines",action="store_true",default=False,help="If you wish to remove detector-specific bad lines, you must specify this flag and --xybadlines.")
	parser.add_argument('--xybadlines', help="Specify the list of bad pixel coordinates for your detector. Will only be used if --fixbadlines is also specified.", nargs=2, default=['3106,3093','3621,3142','4719,3494'])
	parser.add_argument("--rotavg", action="store_true",default=False,help="Use a rotationally averaged coherent power spectrum.")
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-2)
	
	(options, args) = parser.parse_args()
	
	if len(args) < 1:
		print("You must specify at least one ddd movie stack.")
		exit(1)
	
	pid=E2init(sys.argv)
	
	for fname in args:
		
		if not os.path.isfile(fname):
			print("Sorry, {} does not exist".format(fname))
			continue
		else: print("Processing {}".format(fname))
		
		hdr = EMData(fname,0,True)
		if (hdr['nx'] / options.boxsize - 1) < 2 or (hdr['ny'] / options.boxsize - 1) < 2:
			print("You will need to use a smaller box size with your data.")
			sys.exit(1)
		
		if options.gain or options.dark or options.gaink2:
			if options.verbose: print("Correcting frames")
			fname = DirectDetectorUtil.correct_frames(options,fname)
		
		incfile = "{}_incoherent_pws.hdf".format(fname[:-4])
		aligned = "{}_aligned.hdf".format(fname[:-4])
		average = aligned[:-4]+"_average.hdf"
		
		try: os.remove('{}_pairwise_coherent_pws.hdf'.format(fname[:-4]))
		except OSError: pass
		try: os.remove('{}_ccf.hdf'.format(fname[:-4]))
		except OSError: pass
		try: os.remove('{}_best_coherent_pws.hdf'.format(fname[:-4]))
		except OSError: pass
		try: os.remove('{}_pwcc.hdf'.format(fname[:-4]))
		except OSError: pass
		
		n=EMUtil.get_image_count(fname)
		pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
		npairs = len(pairs)
		
		if options.verbose: print("Loading frames from {}".format(fname))
		all_frames = []
		for i in range(n):
			print2line("\t{}/{}".format(i,n))
			f = EMData(fname,i)
			if options.fixbadlines:
				for line in options.xybadlines:
					coords = map(int,line.split(','))
					f.process_inplace('math.xybadline',{'xloc':coords[0],'yloc':coords[1]})
			if options.fixaxes:
				f.process_inplace('filter.xyaxes0',{'neighbor':1,'neighbornorm':options.nnorm,'x':1,'y':1})
			all_frames.append(f)
		
		try: ok = EMUtil.is_same_size(EMData(options.boxsize,options.boxsize),EMData(incfile))
		except: ok = False
		if ok:
			if options.verbose: print("Loading incoherent power spectrum") 
			if options.rotavg: ips = EMData(incfile,1) # 1: ROT-AVG
			else: ips = EMData(incfile,0) # 0: NOT ROT-AVG
		else:
			if options.verbose: print("Generating incoherent power spectrum")
			ips = incoherent_pws(all_frames,bs=options.boxsize)
			ips.do_ift().write_image(incfile,0)
			ips.process_inplace('math.rotationalaverage')
			ips.do_ift().write_image(incfile,1)
		#ips.process_inplace('filter.highpass.gauss',{'cutoff_pixels':1})
		
		ips_fit = fit_ctf(ips)
		
		if options.verbose: print("Generating coefficient matrix")
		A = gen_coeff_matrix(n)
		#print("A: ({} X {})".format(A.shape[0],A.shape[1]))
		
		if options.verbose: print("Computing ordinate vector components")
		bx = []
		by = []
		ccf_xs = []
		ccf_ys = []
		ctr = 0
		for i,(a,b) in enumerate(pairs):
			
			# UCSF
			print2line("{}/{}".format(str(i).rjust(3),npairs))
			ccf_x,ccf_y,ccf = ccf_ordinate(all_frames[a],all_frames[b],i,fname)
			ccf_xs.append(ccf_x)
			ccf_ys.append(ccf_y)
			ccf.write_image('{}_ccf.hdf'.format(fname[:-4]),ctr)
			
			# ALTERNATE
			pair = [all_frames[a],all_frames[b]]
			pc = PairwiseCoherence(fname,pair,ips_fit,bs=options.boxsize)
			x,y,pwcc,cps = pc.maximize(options.maxshift,options.step)
			pwcc.write_image('{}_pwcc.hdf'.format(fname[:-4]),ctr)
			
			bx.append(x) # NEED TO CHECK...is x,y from center, allowing for negative translations?
			by.append(y) # OR is x,y from bottom left so all translations are positive (or negative)?
			ctr+=1
		
		b = np.concatenate([bx,by])
		#print("b: ({})".format(b.shape[0]))
		
		if options.verbose: print("\nOptimizing alignment")
		#results = np.linalg.lstsq(A,b)
		#r=results[0]
		#print("r: ({})".format(r.shape[0]))
		results_x = np.linalg.lstsq(A,bx)
		results_y = np.linalg.lstsq(A,by)
		rx = results_y[0]
		ry = results_y[0]
		
		print("\nFrame\tX\tY")
		avg = Averagers.get('mean')
		for i,(tx,ty) in enumerate(zip(rx,ry)):
			#tx = round(t[0],2)
			#ty = round(t[1],2)
			print("{}/{}\t{}\t{}".format(str(i+1).rjust(2),n,tx,ty))
			f = EMData(fname,i)
			f.process_inplace('xform',{'transform':Transform({'type':'eman','tx':tx,'ty':ty})})
			avg.add_image(f)
			f.write_image(aligned,i)
		avg = avg.finish()
		avg.write_image(average,0)
	
	return

def fit_ctf(ips):
	# fit defocus, 
	# subtract background, 
	#return 2D CTF fit.
	return ips

def gen_coeff_matrix(n): # N frames
	m = n*(n-1)/2
	#A0 = np.zeros([m,n])
	A = np.zeros([m,n])
	for i in range(n):
		a = A[n*i:n*i+n,i:]
		for row in range(a.shape[0]): 
			col = row + 1 # include diagonal
			a[row][:col] = 1
	#Ax0 = np.concatenate([A,A0])
	#Ay0 = np.concatenate([A0,A])
	return A #np.hstack([Ax0,Ay0])

def incoherent_pws(frames,bs=512):
	mx = np.arange(bs+50,frames[0]['nx']-bs+50,bs)
	my = np.arange(bs+50,frames[1]['ny']-bs+50,bs)
	regions = {}
	for i in xrange(len(frames)): 
		regions[i] = [[x,y] for y in my for x in mx]
	ips = Averagers.get('mean')
	for i in xrange(len(frames)):
		print2line("\t{}/{}".format(i+1,len(frames)))
		img = frames[i]
		frame_avg = Averagers.get('mean')
		for r in regions[i]:
			reg = frames[i].get_clip(Region(r[0],r[1],bs,bs))
			reg.process_inplace("normalize.unitlen")
			reg.do_fft_inplace()
			reg.ri2inten()
			frame_avg.add_image(reg)
		ips.add_image(frame_avg.finish())
	ips = ips.finish()
	ips.process_inplace("math.sqrt")
	ips.process_inplace('normalize.edgemean')
	return ips

def ccf_ordinate(f1,f2,i,fname): # Li (2012)
	#This is essentially the UCSF aligner. Need to impose B-factor?
	ccf = f1.calc_ccf(f2)
	#ccf.process_inplace('xform.phaseorigin.tocenter')
	#ccf.write_image('{}_ccf.hdf'.format(fname[:-4]),i)
	#mloc = ccf.calc_max_location()
	#x = ccf.get_xsize()/2-mloc[0]
	#y = ccf.get_ysize()/2-mloc[1]
	ali = f1.align("translational",f2) # align to ccf peak
	#x,y,z = ali.calc_max_location_wrap(-1,-1,-1)
	x,y = ali["xform.align2d"].get_trans_2d()
	print2line("\t\t\t({},{})".format(x,y))
	return x,y,ccf

class PairwiseCoherence:
	
	def __init__(self,fname,frames,ips,bs=512):
		self.fname = fname
		self.frames = frames
		self.ips_r = ips
		self.ips_c = ips.do_fft()
		self.bs = bs
		self.n = len(frames)
		mx = np.arange(self.bs+50,self.frames[0]['nx']-self.bs+50,self.bs)
		my = np.arange(self.bs+50,self.frames[0]['ny']-self.bs+50,self.bs)
		self._tiles = {}
		for i in xrange(self.n):
			self._tiles[i] = np.array([[x,y] for y in my for x in mx])
		self._stacks = {}
		for ir in xrange(len(self._tiles[0])):
			self._stacks[ir] = np.array([self._tiles[i][ir] for i in xrange(self.n)])
		self._orig_tiles = deepcopy(self._tiles)
		self._orig_stacks = deepcopy(self._stacks)
	
	def _calc_coherent_power_spectrum(self):
		cps_avg = Averagers.get('mean')
		for s in xrange(len(self._stacks)):
			stack_avg = Averagers.get('mean')
			for i,r in enumerate(self._stacks[s]):
				reg = Region(r[0],r[1],self.bs,self.bs)
				stack_avg.add_image(self.frames[i].get_clip(reg))
			avg = stack_avg.finish()
			avg.process_inplace('normalize.unitlen')
			avg.do_fft_inplace()
			avg.ri2inten()
			cps_avg.add_image(avg)
		cps = cps_avg.finish()
		cps.process_inplace('math.sqrt')
		cps.process_inplace('normalize.edgemean')
		#cps.process_inplace('filter.highpass.gauss',{'cutoff_pixels':1})
		return cps

	def maximize(self,maxshift=10,incr=1.0):
		shifts = np.arange(-maxshift/2,maxshift/2+incr,incr)
		ls = len(shifts)
		coherence = np.zeros([ls,ls])
		max_coh = 0.0
		for i,tx in enumerate(shifts):
			print2line("\t{}/{}".format(str(i+1).rjust(len(str(ls))),ls))
			for j,ty in enumerate(shifts):
				print2line("\t\t{}/{}".format(str(j+1).rjust(len(str(ls))),ls))
				tr = np.array([tx,ty])
				for n in range(len(self._tiles[1])):
					self._tiles[1][n] = np.add(self._orig_tiles[1][n],tr)
					self._stacks[n][1] = np.add(self._orig_stacks[n][1],tr)
				cps = self._calc_coherent_power_spectrum()
				cps.process_inplace('xform.phaseorigin.tocenter')
				cps.process_inplace('filter.lowpass.gauss',{'cutoff_abs':0.4})
				cps.do_ift().write_image('{}_pairwise_coherent_pws.hdf'.format(self.fname[:-4]),-1)
				#cps.process_inplace('filter.highpass.gauss',{'cutoff_pixels':1})
				coh = -1*cps.cmp('frc',self.ips_c) 
				#np.power(cps.cmp('dot',self.ips_c),2)
				#np.sqrt(-1*cps.cmp('dot',cps))
				if coh > max_coh:
					max_coh = coh
					best_cps = cps
				coherence[i,j] = coh
		pwcc = from_numpy(coherence)
		x,y,z = pwcc.calc_max_location_wrap(-1,-1,0)
		print2line("\t\t\t({},{})".format(x,y))
		best_cps.do_ift().write_image('{}_best_coherent_pws.hdf'.format(self.fname[:-4]),-1)
		return x,y,pwcc,best_cps

class DirectDetectorUtil:
	
	@classmethod
	def correct_frames(cls,options,fname,outfile=None):
		options.path = fname
		if options.dark:
			if options.verbose > 8: print("Generating dark reference")
			dark = cls.dark_correct(options)
		else: dark = None
		if options.gain:
			if options.verbose > 8: print("Generating gain reference")
			gain = cls.gain_correct(options,dark)
		elif options.gaink2:
			if options.verbose > 8: print("Generating K2 gain reference")
			gain = EMData(options.gaink2)
		else: gain = None
		hdr = EMData(options.path,0,True)
		if options.path[-4:].lower() in (".mrc"): nd = hdr['nz']
		else: nd = EMUtil.get_image_count(options.path)
		if not outfile: outfile = options.path[:-4] + "_corrected.hdf"
		for i in xrange(nd):
			if options.verbose:
				print "Correcting frame: {}/{}	\r".format(i+1,nd),
				sys.stdout.flush()
			if options.path[-4:].lower() in (".mrc"):
				r = Region(0,0,i,nx,ny,1)
				im=EMData(path,0,False,r)
			else: im=EMData(options.path,i)
			if dark: im.sub(dark)
			if gain: im.mult(gain)
			im.process_inplace("threshold.clampminmax",{"minval":0,"maxval":im["mean"]+im["sigma"]*3.5,"tozero":1})
			# fixes clear outliers as well as values which were exactly zero
			im.process_inplace("threshold.outlier.localmean",{"sigma":3.5,"fix_zero":1}) 
			im.process_inplace("normalize.edgemean")
			im.write_image(outfile,i)
		else: outfile = options.path
		return outfile
	
	@classmethod
	def dark_correct(cls,options):
		hdr = EMData(options.dark,0,True)
		if options.path[-4:].lower() in (".mrc"): nd = hdr['nz']
		else: nd = EMUtil.get_image_count(options.dark)
		dark=EMData(options.dark,0)
		if nd>1:
			sigd=dark.copy()
			sigd.to_zero()
			a=Averagers.get("mean",{"sigma":sigd,"ignore0":1})
			for i in xrange(0,nd):
				if options.verbose:
					print "Summing dark: {}/{}	\r".format(i+1,nd),
					sys.stdout.flush()
				t=EMData(options.dark,i)
				t.process_inplace("threshold.clampminmax",{"minval":0,"maxval":t["mean"]+t["sigma"]*3.5,"tozero":1})
				a.add_image(t)
			dark=a.finish()
			sigd.write_image(options.dark.rsplit(".",1)[0]+"_sig.hdf")
			if options.fixbadlines:
				# Theoretically a "perfect" pixel would have zero sigma, but in reality, the opposite is true
				sigd.process_inplace("threshold.binary",{"value":sigd["sigma"]/10.0})  
				dark.mult(sigd)
			dark.write_image(options.dark.rsplit(".",1)[0]+"_sum.hdf")
		dark.process_inplace("threshold.clampminmax.nsigma",{"nsigma":3.0})
		dark2=dark.process("normalize.unitlen")
		return dark
	
	@classmethod
	def gain_correct(cls,options,dark):
		hdr = EMData(options.gain,0,True)
		if options.path[-4:].lower() in (".mrc"): nd = hdr['nz']
		else: nd = EMUtil.get_image_count(options.gain)
		gain=EMData(options.gain,0)
		if nd>1:
			sigg=gain.copy()
			sigg.to_zero()
			a=Averagers.get("mean",{"sigma":sigg,"ignore0":1})
			for i in xrange(0,nd):
				if options.verbose:
					print "Summing gain: {}/{}	\r".format(i+1,nd),
					sys.stdout.flush()
				t=EMData(options.gain,i)
				t.process_inplace("threshold.clampminmax",{"minval":0,"maxval":t["mean"]+t["sigma"]*3.5,"tozero":1})
				a.add_image(t)
			gain=a.finish()
			sigg.write_image(options.gain.rsplit(".",1)[0]+"_sig.hdf")
			# Theoretically a "perfect" pixel would have zero sigma, but in reality, the opposite is true
			if options.fixbadlines: sigg.process_inplace("threshold.binary",{"value":sigg["sigma"]/10.0})
			if dark!=None:
				sigd=EMData(options.dark.rsplit(".",1)[0]+"_sig.hdf",0,False)
				sigg.mult(sigd)
			gain.mult(sigg)
			gain.write_image(options.gain.rsplit(".",1)[0]+"_sum.hdf")
		if dark!=None : gain.sub(dark)	# dark correct the gain-reference
		# normalize so gain reference on average multiplies by 1.0
		gain.mult(1.0/gain["mean"])
		# setting zero values to zero helps identify bad pixels
		gain.process_inplace("math.reciprocal",{"zero_to":0.0})	 
		return gain

def print2line(text):
	sys.stdout.write(text+"\r")
	sys.stdout.flush()

# def coherent_pws(frames,bs=512):
# 	mx = np.arange(bs+50,frames[0]['nx']-bs+50,bs)
# 	my = np.arange(bs+50,frames[1]['ny']-bs+50,bs)
# 	regions = {}
# 	for i in xrange(len(frames)):
# 		regions[i] = [[x,y] for y in my for x in mx]
# 	stacks = {}
# 	for ir in xrange(len(regions[0])):
# 		stacks[ir] = [regions[i][ir] for i in xrange(len(frames))]
# 	cps = Averagers.get('mean')
# 	for s in xrange(len(stacks)):
# 		print2line("\t\t{}/{}".format(str(s+1).rjust(3),len(stacks)))
# 		stack_avg = Averagers.get('mean')
# 		for i,r in enumerate(stacks[s]):
# 			f = frames[i]
# 			stack_avg.add_image(f.get_clip(Region(r[0],r[1],bs,bs)))
# 		avg = stack_avg.finish()
# 		avg.process_inplace('normalize.unitlen')
# 		avg.do_fft_inplace()
# 		avg.ri2inten()
# 		cps.add_image(avg)
# 	cps = cps.finish()
# 	cps.process_inplace('math.sqrt')
# 	cps.process_inplace('normalize.edgemean')
# 	return cps

# def cps_ordinate(f1,f2,ips): #f1 is stationary
# 	# compute "cross coherence" between frames
# 	c = PairwiseCoherence([f1,f2],ips)
# 	# To determine x/y translations that maximize power:
# 		# Use integrated coherent power spectrum (may require less signal)
# 		# OR use "COHERENCE" in x and y between CPS/IPS as attempted before
# 	x,y = c.maximize()
# 	# return relative x and y shifts
# 	return x,y

# in MAIN
		# cohfile = "{}_coherent.hdf".format(fname[:-4])
		# pairfile = '{}_pairwise_coherent_pws.hdf'.format(fname[:-4])
		# if not os.path.isfile(cohfile):
		# 	print("\tCoherent (ALL FRAMES)")
		# 	cps_all = coherent_pws(all_frames)
		# 	cps_all.do_ift().write_image(cohfile,0)
		# if not os.path.isfile(pairfile):
		# 	print("\tCoherent (PAIRWISE)")
		# 	for i,(a,b) in enumerate(pairs):
		# 		print2line("\t{}/{}".format(i+1,npairs))
		# 		cps = coherent_pws([all_frames[a],all_frames[b]],bs=options.boxsize)
		# 		cps.do_ift().write_image(pairfile,i)

			#cps = EMData(pairfile,i)
			#x,y = ccf_ordinate(all_frames[pair[0]],all_frames[pair[1]])
			#x,y = cps_ordinate(f1,f2,ips)

if __name__ == "__main__":
	main(sys.argv[1:])