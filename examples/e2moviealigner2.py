#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')
from Simplex import Simplex
from EMAN2 import *
import numpy as np
import matplotlib.pyplot as plt

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """e2moviealigner.py [options] <ddd_movie_stack> <ddd_movie_stack> ... <ddd_movie_stack>

	Determines the optimal whole-frame alignment of a DDD movie. It can be
	used to determine the proper x and y translations for each movie frame and
	can perform the actual alignment according to the translations.
	
	Example: e2moviealigner.py --dark dark.hdf --gain gain.hdf movie.hdf -v 9
	"""
	
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_argument("--dark",type=str,default=None,help="Perform dark image correction using the specified image file")
	parser.add_argument("--gain",type=str,default=None,help="Perform gain image correction using the specified image file")
	parser.add_argument("--gaink2",type=str,default=None,help="Perform gain image correction. Gatan K2 gain images are the reciprocal of DDD gain images.")
	parser.add_argument("--boxsize", type=int, help="Set the boxsize used to compute power spectra across movie frames",default=512)
	parser.add_argument("--maxshift", type=int, help="Set the maximum radial frame translation distance (in pixels) from the initial frame alignment.",default=1.5)
	parser.add_argument("--step",type=str,default="0,1",help="Specify <first>,<step>,[last]. Processes only a subset of the input data. ie- 0,2 would process all even particles. Same step used for all input files. [last] is exclusive. Default= 0,1 (first image skipped)")
	parser.add_argument("--fixbadpixels",action="store_true",default=False,help="Tries to identify bad pixels in the dark/gain reference, and fills images in with sane values instead")
	parser.add_argument('--xybadlines', help="Specify the list of bad pixel coordinates for your detector. Will only be used if --fixbadpixels is also specified.", nargs=2, default=['3106,3093','3621,3142','4719,3494'])
	parser.add_argument("--maxiters", type=int, help="Set the maximum number of iterations for the simplex minimizer to run before stopping. Default is 250.",default=250)
	parser.add_argument("--epsilon", type=float, help="Set the learning rate for the simplex minimizer. Smaller is better, but takes longer. Default is 0.001.",default=0.001)
	parser.add_argument("--show", action="store_true",help="Show average of movie frames before and after alignment.",default=False)
	parser.add_argument("--onlycorrect", action="store_true",help="Only correct the data. Do not align.",default=False)
	parser.add_argument("--nocorrection", action="store_true",help="Do not correct data. Only align.",default=False)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-2)

	(options, args) = parser.parse_args()

	if len(args) < 1:
		print usage
		parser.error("You must specify at least one input DDD stack.")
		sys.exit(1)
	
	if options.boxsize: bs = options.boxsize
	hdr = EMData(args[0],0,True)
	if (hdr['nx'] / bs - 1) < 2 or (hdr['ny'] / bs - 1) < 2:
		print("You will need to use a smaller box size with your data.")
		sys.exit(1)

	if options.maxshift < 0:
		print("The parameter '--maxshift' must be greater than 0")
		sys.exit(1)
	
	if options.xybadlines: options.xybadlines = [map(int,s.split(',')) for s in options.xybadlines]
	
	pid=E2init(sys.argv)

	for fname in args:
		
		print("Processing {}".format(fname))
		
		if not options.nocorrection:
			if options.gain or options.dark or options.gaink2 or options.fixbadpixels:
				if options.verbose: print("Correcting movie frames...")
				corrected = MovieAligner.correct_frames(options,fname)
			elif options.onlycorrect:
				if not options.gain and not options.gaink2: print("No gain reference to apply.")
				if not options.dark: print("No dark reference to apply")
				if not options.fixbadpixels: print("Not fixing bad pixels")
				sys.exit(1)
		else: corrected = fname
		
		if not options.onlycorrect:
			if options.verbose: print("Computing alignment parameters...")
			alignment = MovieAligner(corrected,bs=options.boxsize)
			if options.verbose: print("Optimizing alignment...")
			alignment.optimize(options)
			if options.show: display([alignment.before,alignment.after])
	
	E2end(pid)

class MovieAligner:
	
	def __init__(self, fname, bs):
		self.path = fname
		self.hdr = EMData(fname,0,True)
		self.hdr['nimg'] = EMUtil.get_image_count(fname)
		self.base = fname.split('.')[0]
		self.frames = [EMData(fname,i) for i in xrange(self.hdr['nimg'])]
		self.translations = np.zeros((self.hdr['nimg'],2))
		self.before = self.save(descriptor="unaligned",save_frames=False)
		self.fullbox=good_boxsize(min(self.hdr['nx']-10,self.hdr['ny']-10),False)
		self.bs = bs
		if bs >= good_boxsize: bs = self.fullbox
		mx = np.arange(50,self.hdr['nx']-50,bs)
		my = np.arange(50,self.hdr['ny']-50,bs)
		self._regions = {}
		for i in xrange(self.hdr['nimg']):
			self._regions[i] = np.array([[x,y] for y in my for x in mx])
		self._stacks = {}
		for ir in xrange(len(self._regions[0])):
			self._stacks[ir] = [self._regions[i][ir] for i in xrange(self.hdr['nimg'])]
		self.iteration = 0
		self.update_energy()
		self.translations = np.random.randint(-4,4,size=(self.hdr['nimg'],2))
	
	def update_energy(self):
		if self.iteration == 0:
			self.calc_incoherent_pws()
			self.write_ips("incoherent_pws.hdf")
		self.calc_coherent_pws()
		self.write_cps("coherent_pws.hdf",n=self.iteration)
		self.energy = self.calc_energy()
		self.iteration += 1
	
	def calc_energy(self):
		inc_pws = np.asarray(self.ips_ctf_fit)
		coh_pws = np.asarray(self.oned_cps)
		return -np.dot(inc_pws,coh_pws)
	
	def calc_incoherent_pws(self):
		ips = Averagers.get('mean')
		for i in xrange(self.hdr['nimg']):
			img = self.frames[i]
			img_ips = Averagers.get('mean')
			for r in self._regions[i]:
				reg = self.frames[i].get_clip(Region(r[0],r[1],self.bs,self.bs))
				reg.process_inplace("normalize.unitlen")
				reg.do_fft_inplace()
				reg.ri2inten()
				img_ips.add_image(reg)
			ips.add_image(img_ips.finish())
		self.ips = ips.finish()
		self.ips.process_inplace("math.sqrt") 
		self.ips["is_intensity"]=0 
		self.ips.process_inplace('math.rotationalaverage')
		self.real_ips = self.ips.do_ift()
		self.oned_ips, self.ips_ctf, self.ips_ctf_fit  = self.fit_defocus(self.ips)
	
	def calc_coherent_pws(self):
		cps = Averagers.get('mean')
		for s in xrange(len(self._stacks)):
			stack_cps = Averagers.get('mean')
			for i,r in enumerate(self._stacks[s]):
				stack_cps.add_image(self.frames[i].get_clip(Region(r[0],r[1],self.bs,self.bs)))
			avg = stack_cps.finish()
			avg.process_inplace('normalize.unitlen') 
			avg.do_fft_inplace()
			avg.ri2inten()
			cps.add_image(avg)
		self.cps = cps.finish()
		self.cps.process_inplace('math.sqrt')
		self.cps["is_intensity"]=0 
		self.real_cps = self.cps.do_ift()
		self.oned_cps, self.cps_ctf, self.cps_ctf_fit = self.fit_defocus(self.cps)
	
	def write_cps(self,name,n=0):
		self.real_cps['is_intensity'] = 0
		self.real_cps.write_image(name,n)
	
	def write_ips(self,name,n=0):
		self.real_ips['is_intensity'] = 0
		self.real_ips.write_image(name,n)
	
	def fit_defocus(self,img):
		ds=1.0/(img["apix_x"]*self.bs)
		ns=min(int(floor(.25/ds)),img["ny"]/2)
		oned=np.array(img.calc_radial_dist(ns,0,1.0,1)[1:])
		oned=np.log10(oned)
		oned-=min(oned)
		oned/=max(oned)
		ctf=EMAN2Ctf()
		ctf.voltage=300.0
		ctf.cs=4.7
		ctf.ampcont=10
		ctf.bfactor=0
		ctf.dfdiff=0
		s=[ds*(i+1) for i,j in enumerate(oned)]
		dfl=[]
		ql=[]
		for df in np.arange(0.25,5,.01):
			ctf.defocus=df
			curve=np.array(ctf.compute_1d(ns*2,ds,Ctf.CtfType.CTF_AMP)[1:])
			curve*=curve
			#curve-=0.5	# duh
			zeros=[int(ctf.zero(i)/ds) for i in xrange(15)]
			zeros=[i for i in zeros if i<len(curve) and i>0]
			onedbg=self.bgsub(oned,zeros)
			qual=curve.dot(onedbg)
			dfl.append(df)
			ql.append(qual)
		a=np.argmax(ql)
		df=dfl[a]
		ctf.defocus = df
		curve=np.array(ctf.compute_1d(ns*2,ds,Ctf.CtfType.CTF_AMP)[1:])
		curve*=curve
		#curve-=0.5
		zeros=[int(ctf.zero(i)/ds) for i in xrange(15)]
		zeros=[i for i in zeros if i<len(curve) and i>0]
		onedbg=self.bgsub(oned,zeros)
		qual=curve.dot(onedbg)
		ctf.qual = qual
		return onedbg, ctf, curve
	
	@staticmethod
	def bgsub(curve,zeros):
		floc=min(zeros[0]/2,8)
		itpx=[curve[:floc].argmin()]+list(zeros)+[len(curve)-1]
		itpy=[min(curve[i-1:i+2]) for i in itpx]
		itpy[0]=curve[:floc].min()
		itpx=np.array(itpx)
		itpy=np.array(itpy)
		ret=curve-np.interp(range(len(curve)),itpx,itpy)
		ret[:floc]=0
		return ret
	
	def optimize(self, options):
		if options.verbose: print("Performing simplex minimization")
		state = self.translations
		sm = Simplex(self._compares,state,[1]*len(state),data=self)
		result, error, iters = sm.minimize(options.epsilon,options.maxiters,monitor=1)
		self.translations = result
		print("\nFrame\tTranslation")
		with open(self.path[:-4]+"_results.txt",'w') as results:
			title = "\nError: {}\nIters: {}".format(error,iters)
			results.write(title+"\n")
			for i,t in enumerate(self.translations):
				info = "{}\t( {}, {} )".format(i+1,t[0],t[1])
				print(info)
				results.write(info+"\n")
		print("")
		self.after = self.save(descriptor="aligned",save_frames=True)
	
	def update_frame_params(self,fnum,trans):
		self.translations[fnum] = trans # update translations
		for rnum,reg in enumerate(self.regions[fnum]): # update regions first
			self._regions[fnum][rnum] = np.add(reg,trans) # add transform vector from "origin"
		for rnum,reg in enumerate(self._regions[fnum]): # update stacks after regions are updated
			self._stacks[rnum][fnum] = reg
	
	@staticmethod
	def _compares(translations,aligner):
		middle = aligner.hdr['nimg']/2
		for i,t in enumerate(translations):
			ct = aligner.translations[i]
			if i != middle and t[0] != ct[0] and t[1] != ct[1]:
				aligner.update_frame_params(i,t)
		aligner.update_energy()
		return aligner.energy
	
	def save(self,descriptor,save_frames=False):
		if save_frames:
			tofile = self.path[:-4]+"_{}_frames.hdf".format(descriptor)
			for i,t in enumerate(self.translations):
				tform = Transform({'type':'eman','tx':int(t[0]),'ty':int(t[1])})
				f = self.frames[i].copy()
				f.transform(tform)
				if save_frames: f.write_image(tofile,i)
			fromfile = tofile
		else: fromfile = self.path
		avg = Averagers.get('mean')
		for i in range(len(self.frames)):
			f = EMData(fromfile,i)
			avg.add_image(f)
		average = avg.finish()
		average.write_image(self.path[:-4]+"_{}_average.hdf".format(descriptor))
		return average
	
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
			if options.fixbadpixels:
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
			if options.fixbadpixels: sigg.process_inplace("threshold.binary",{"value":sigg["sigma"]/10.0})
			if dark!=None:
				sigd=EMData(options.dark.rsplit(".",1)[0]+"_sig.hdf",0,False)
				sigg.mult(sigd)
			gain.mult(sigg)
			gain.write_image(options.gain.rsplit(".",1)[0]+"_sum.hdf")
		if dark!=None : gain.sub(dark)	# dark correct the gain-reference
		gain.mult(1.0/gain["mean"])	 # normalize so gain reference on average multiplies by 1.0
		gain.process_inplace("math.reciprocal",{"zero_to":0.0})	 # setting zero values to zero helps identify bad pixels
		return gain
	
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
		if dark != None or gain != None and options.verbose >= 8: print("Applying references")
		hdr = EMData(options.path,0,True)
		if options.path[-4:].lower() in (".mrc"): nd = hdr['nz']
		else: nd = EMUtil.get_image_count(options.path)
		step = options.step.split(",")
		if len(step) == 3: last = int(step[2])
		else: last = nd
		first = int(step[0])
		step  = int(step[1])
		if options.dark or options.gain or options.gaink2 or options.fixbadpixels:
			if not outfile: outfile = options.path[:-4] + "_cleaned.hdf"
			for i in xrange(first,last,step):
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
				if options.fixbadpixels:
					# fixes clear outliers as well as values which were exactly zero
					im.process_inplace("threshold.outlier.localmean",{"sigma":3.5,"fix_zero":1})
					im.process_inplace('filter.xyaxes0')
					for bl in options.xybadlines:
						# remove detector-specific bad lines
						im.process_inplace('math.xybadline',{'xloc':bl[0],'yloc':bl[1]})
				im.process_inplace("normalize.edgemean")
				im.write_image(outfile,i-first)
		else: outfile = options.path
		return outfile

if __name__ == "__main__":
	main()
