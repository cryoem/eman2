#!/usr/bin/env python

# Author: James Michael Bell, jmbell@bcm.edu, 10/18/2015

import matplotlib
#matplotlib.use('Agg')
from Simplex import Simplex
from EMAN2 import *
import numpy as np
import matplotlib.pyplot as plt
try:
	from scipy import optimize
except:
	print("You must install scipy to use this program.")
	print("To do this, you could use either 'sudo pip install scipy' or 'sudo easy_install scipy'.")
	sys.exit(1)

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
	parser.add_argument("--maxshift", type=int, help="Set the maximum frame translation distance in pixels.",default=10)
	parser.add_argument("--maxiter", type=int, help="Set the maximum iterations for optimization (this is NOT equivalent to function evaluations).",default=25)
	#parser.add_argument("--stepsize", type=int, help="Set the step size for optimization.",default=1.0)
	#parser.add_argument("--temperature", type=float, help="Set the temperature (rejection criterion) for optimization.",default=0.5)
	parser.add_argument("--neighbornorm", type=int, help="Set the norm to be used for fixing axes. Default is 2",default=2)
	#parser.add_argument("--popsize", type=int, help="Population size to be used in differential evolution.",default=25)
	parser.add_argument("--fixaxes",action="store_true",default=False,help="Tries to identify bad pixels and fill them in with sane values instead")
	parser.add_argument("--fixbadlines",action="store_true",default=False,help="If you wish to remove detector-specific bad lines, you must specify this flag and --xybadlines.")
	parser.add_argument('--xybadlines', help="Specify the list of bad pixel coordinates for your detector. Will only be used if --fixbadlines is also specified.", nargs=2, default=['3106,3093','3621,3142','4719,3494'])
	parser.add_argument("--show", action="store_true",help="Show average of movie frames before and after alignment.",default=True)
	parser.add_argument("--plot",action="store_true",default=False,help="Plot the 1D power spectra and exit.")
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-2)
	
	(options, args) = parser.parse_args()
	
	if len(args) < 1:
		print usage
		parser.error("You must specify at least one input DDD movie stack.")
		sys.exit(1)
	
	if options.boxsize: bs = options.boxsize
	hdr = EMData(args[0],0,True)
	if (hdr['nx'] / bs - 1) < 2 or (hdr['ny'] / bs - 1) < 2:
		print("You will need to use a smaller box size with your data.")
		sys.exit(1)
	
	pid=E2init(sys.argv)
	
	for fname in args:
		
		if options.verbose and len(args) > 1: print("Processing {}".format(fname))
		
		if options.gain or options.dark or options.gaink2:
			if options.verbose: print("Correcting frames")
			fname = MovieAligner.correct_frames(options,fname)
		
		if options.verbose: print("Acquiring parameters")
		alignment = MovieAligner(options,fname)
		
		if options.verbose: print("Optimizing alignment")
		alignment.optimize(options)
		
		if options.show: display([alignment.before,alignment.after])
	
	E2end(pid)

class MovieAligner:
	
	def __init__(self, options, fname):
		self.path = fname
		self.bs = options.boxsize
		self.hdr = EMData(fname,0,True)
		self.hdr['nimg'] = EMUtil.get_image_count(fname)
		self.base = fname.split('.')[0]
		self.nnorm = options.neighbornorm
		self.frames = []
		for i in xrange(self.hdr['nimg']):
			f = EMData(fname,i)
			if options.fixbadlines:
				for line in options.xybadlines:
					coords = map(int,line.split(','))
					f.process_inplace('math.xybadline',{'xloc':coords[0],'yloc':coords[1]})
			self.frames.append(f)
		self.translations = np.zeros((self.hdr['nimg'],2))
		self.best_translations = np.zeros((self.hdr['nimg'],2))
		self.fixaxes = options.fixaxes
		self.before = self.save(descriptor="unaligned",save_frames=False) 
		mx = np.arange(self.bs+50,self.hdr['nx']-self.bs+50,self.bs)
		my = np.arange(self.bs+50,self.hdr['ny']-self.bs+50,self.bs)
		self._regions = {}
		for i in xrange(self.hdr['nimg']):
			self._regions[i] = [[x,y] for y in my for x in mx]
		self._stacks = {}
		for ir in xrange(len(self._regions[0])):
			self._stacks[ir] = [self._regions[i][ir] for i in xrange(self.hdr['nimg'])]
		self.iter = 0
		self.static = int(self.hdr['nimg']/3)
		self.energies = [np.inf]
		self.verbose = options.verbose
		self.plot = options.plot
		self.cps_counter = 0
		#self.low_cutoff = 0.07
		self.calc_energy()
	
	def calc_incoherent_pws(self):
		ips = Averagers.get('mean')
		for i in xrange(self.hdr['nimg']):
			img = self.frames[i]
			frame_avg = Averagers.get('mean')
			for r in self._regions[i]:
				reg = self.frames[i].get_clip(Region(r[0],r[1],self.bs,self.bs))
				reg.process_inplace("normalize.unitlen")
				reg.do_fft_inplace()
				reg.ri2inten()
				frame_avg.add_image(reg)
			ips.add_image(frame_avg.finish())
		self.ips = ips.finish()
		self.ips.process_inplace("math.sqrt")
		self.ips.process_inplace('normalize.edgemean')
		#self.ips.process_inplace('math.rotationalaverage')
		self.oned_ips, self.ips_ctf_fit  = self.fit_defocus(self.ips)
	
	def calc_coherent_pws(self):
		cps = Averagers.get('mean')
		for s in xrange(len(self._stacks)):
			stack_avg = Averagers.get('mean')
			for i,r in enumerate(self._stacks[s]):
				stack_avg.add_image(self.frames[i].get_clip(Region(r[0],r[1],self.bs,self.bs)))
			avg = stack_avg.finish()
			avg.process_inplace('normalize.unitlen')
			avg.do_fft_inplace()
			avg.ri2inten()
			cps.add_image(avg)
		self.cps = cps.finish()
		self.cps.process_inplace('math.sqrt')
		self.cps.process_inplace('normalize.edgemean')
		self.oned_cps, self.cps_ctf_fit = self.fit_defocus(self.cps)
	
	def calc_energy(self):
		if self.iter == 0:
			self.calc_incoherent_pws()
			self.write_ips()
		self.calc_coherent_pws()
		if self.iter == 0: self.write_cps()
		# get (windowed) 1D power spectra
		start = int(len(self.oned_cps)*(1./8.))
		stop = int(len(self.oned_cps)*(7./8.))
		#filter out low resolution
		ips = self.ips_ctf_fit[start:stop]
		cps = self.oned_cps[start:stop]
		# normalize ips/cps
		#ips/=np.sqrt(ips.dot(ips))
		#cps/=np.sqrt(cps.dot(cps))
		#ramp=np.linspace(0.5,1.25,len(cps))
		#ramp=np.exp(ramp)
		#iwt = 1.0
		#cwt = 1.0
		if self.plot:
			if self.iter % 50 == 0:
				plt.figure()
				plt.plot(ips,label='ips')
				plt.plot(cps,label='cps')
				#plt.plot(ips*ramp*iwt,label='ips * {} * ramp'.format(iwt))
				#plt.plot(cps*ramp*cwt,label='cps * {} * ramp'.format(cwt))
				plt.legend()
				plt.tight_layout()
				plt.show()
		# compute dot product between ips/cps
		#print(np.std(cps))
		#print(np.mean(cps))
		energy = -np.dot(ips,cps)
		# if we're subtracting the background, the minima should be ~0
		if energy < min(self.energies) and min(cps) < 0.25: 
			self.write_cps()
			self.best_translations = self.translations
		self.energies.append(energy)
		self.iter += 1
		return energy
	
	def write_cps(self,fname="_pws_coherent.hdf",n=None):
		if self.cps_counter == 0:
			print("Removing old coherent power spectrum")
			try: os.remove(self.path[:-4]+fname)
			except OSError: pass
		self.cps[0,0] = 0
		self.real_cps = self.cps.do_ift()
		self.real_cps.process_inplace('filter.xyaxes0',{'neighbor':1,'neighbornorm':self.nnorm,'x':1,'y':1})
		self.real_cps['is_intensity'] = 0
		if n == None:
			self.real_cps.write_image(self.path[:-4]+fname,self.cps_counter)
			self.cps_counter += 1
		else: self.real_cps.write_image(self.path[:-4]+fname,n)
	
	def write_ips(self,fname="_pws_incoherent.hdf"):
		self.ips[0,0] = 0
		self.real_ips = self.ips.do_ift()
		self.real_ips.process_inplace('filter.xyaxes0',{'neighbor':1,'neighbornorm':self.nnorm,'x':1,'y':1})
		self.real_ips['is_intensity'] = 0
		self.real_ips.write_image(self.path[:-4]+fname,0)
	
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
		bgl=[]
		for df in np.arange(0.6,5,.01):
			ctf.defocus=df
			curve=np.array(ctf.compute_1d(ns*2,ds,Ctf.CtfType.CTF_AMP)[1:])
			curve*=curve
			zeros=[int(ctf.zero(i)/ds) for i in xrange(15)]
			zeros=[i for i in zeros if i<len(curve) and i>0]
			onedbg,bg=self.bgsub(oned,zeros)
			qual=curve.dot(onedbg)
			dfl.append(df)
			ql.append(qual)
			bgl.append(bg)
		a=np.argmax(ql)
		ctf.defocus = dfl[a]
		curve=np.array(ctf.compute_1d(ns*2,ds,Ctf.CtfType.CTF_AMP)[1:])
		curve*=curve
		if self.fixaxes: 
			img.process_inplace('filter.xyaxes0',{'neighbor':1,'neighbornorm':self.nnorm,'x':1,'y':1})
			oned=np.array(img.calc_radial_dist(ns,0,1.0,1)[1:])
			oned=np.log10(oned)
			oned-=min(oned)
			oned/=max(oned)
		onedbg=oned-bgl[a]
		qual=curve.dot(onedbg)
		return np.asarray(onedbg), np.asarray(curve)
	
	@staticmethod
	def bgsub(curve,zeros):
		floc=min(zeros[0]/2,8)
		itpx=[curve[:floc].argmin()]+list(zeros)+[len(curve)-1]
		itpy=[min(curve[i-1:i+2]) for i in itpx]
		itpy[0]=curve[:floc].min()
		itpx=np.array(itpx)
		itpy=np.array(itpy)
		bg = np.interp(range(len(curve)),itpx,itpy)
		ret=curve-bg
		ret[:floc]=0
		return ret,bg
	
	def optimize(self, options, bounds=None):
		if options.verbose: disp=True
		else: disp = False
		
		if bounds == None:
			ms = options.maxshift
			c1 = self.static
			c2 = self.hdr['nimg'] - self.static
			ub1 = [round(ms*((c1+1)-(i+1))/(c1+1),1) for i in range(c1)]
			ub2 = [round(ms*(i+1)/(c2+1),1)+1 for i in range(c2)]
			bounds = []
			for bound in ub1+ub2:
				b = int(round(bound))
				bounds.append((-b,b)) # X
				bounds.append((-b,b)) # Y
		self.bounds = bounds
		
		# if options.verbose > 8:
		# 	print("Initializing optimization with the following bounds:")
		# 	print("Frame\tMaximum Shift Magnitude")
		# 	bds = np.array(bounds).reshape((self.hdr['nimg'],2,2)).astype(int)
		# 	for i,bd in enumerate(bds):
		# 		print("{}\t{}".format(i+1,round(np.sqrt(bd[0,1]**2+bd[1,1]**2),0)))
		# 	print("")
		
		init_state = np.random.randint(-1,1,size=(self.hdr['nimg'],2)).flatten()
		self.init_state = init_state
		
		scb = SimplexCallBack(self,options)
		cb = None #CallBack(self,options)
		ts = None #TakeStep(self,options)
		at = None #AcceptTest(self,options)
		
		#res = optimize.basinhopping(self._compares, init_state, minimizer_kwargs={'method':'Nelder-Mead','args':self,'callback':scb}, T=options.temperature, niter=options.maxiter, disp=disp, callback=cb, take_step=ts, accept_test=at)
		#res = optimize.differential_evolution(self._compares,bounds,args=(self,),polish=False, maxiter=options.maxiter,popsize=15,disp=disp)
		res = optimize.minimize(self._compares, init_state, method='Nelder-Mead', callback=scb, options={'maxiter':options.maxiter,'disp': True}, args=self)
		
		print("\nFrame\tTranslation")
		with open(self.path[:-4]+"_results.txt",'w') as results:
			info = "\nEnergy: {}\nIters: {}\nFunc Evals: {}\n".format(res.fun,res.nit,res.nfev)
			results.write(info)
			for i,t in enumerate(self.best_translations):
				info = "{}\t( {}, {} )".format(i+1,int(t[0]),int(t[1]))
				print(info); results.write(info+"\n")
		self.after = self.save(descriptor="aligned",save_frames=True)
	
	@staticmethod
	def _compares(ts,aligner):
		translations = np.asarray(ts).reshape((aligner.hdr['nimg'],2))
		translations = np.round(translations).astype(int) # only allow integral shifts
		for frame_num,trans in enumerate(translations):
			if frame_num != aligner.static:
				if not bool(np.all(trans == aligner.translations[frame_num])):
					aligner.translations[frame_num] = trans
					for reg_num in range(len(aligner._regions[frame_num])):
						aligner._regions[frame_num][reg_num] = trans
						aligner._stacks[reg_num][frame_num] = trans
		energy = aligner.calc_energy()
		return energy
	
	def save(self,descriptor,save_frames=False):
		frame_file = self.path[:-4]+"_{}_frames.hdf".format(descriptor)
		average_file = self.path[:-4]+"_{}_average.hdf".format(descriptor)
		avg = Averagers.get('mean')
		for i,t in enumerate(self.translations):
			f = self.frames[i].copy()
			if self.fixaxes:
				f.process_inplace('filter.xyaxes0',{'neighbor':1,'neighbornorm':self.nnorm,'x':1,'y':1})
			f.translate(int(round(t[0])),int(round(t[1])),0)
			if save_frames: f.write_image(frame_file,i)
			avg.add_image(f)
		average = avg.finish()
		average.write_image(average_file)
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

class TakeStep(object): # Impose soft bounds on the optimization
	
	def __init__(self, aligner, options):
		self.maxshift = options.maxshift
		self.static = aligner.static
	
	def __call__(self, x):
		step = np.random.uniform(-self.maxshift/2, self.maxshift/2, x.shape).astype(int)
		step[2*self.static] = 0
		step[2*self.static+1] = 0
		return x+step

class AcceptTest(object): # Impose hard bounds on the optimiation
	
	def __init__(self, aligner, options):
		self.xmin = [min(i) for i in aligner.bounds]
		self.xmax = [max(i) for i in aligner.bounds]
	
	def __call__(self, **kwargs):
		x = kwargs["x_new"]
		tmax = bool(np.all(x <= self.xmax))
		tmin = bool(np.all(x >= self.xmin))
		return tmax and tmin

class CallBack(object):
	
	def __init__(self, aligner, options):
		self.aligner = aligner
		self.counter = 0
	
	def __call__(self,x,fx,accepted):
		i = str(self.aligner.iter)
		n = str(self.aligner.cps_counter)
		b = str(min(self.aligner.energies))
		w = str(max(self.aligner.energies[1:]))
		c = str(self.aligner.energies[-1])
		if self.aligner.verbose: 
			print("fEvals: {} CPS: {} Energy: {} Accepted: {}".format(i,n,fx,accepted))
		self.aligner.write_cps(fname="_pws_coherent_callback.hdf",n=self.counter)
		self.counter += 1
		#self.aligner.best_translations = self.aligner.translations

class SimplexCallBack(object):
	
	def __init__(self, aligner, options):
		self.aligner = aligner
		self.counter = 0
	
	def __call__(self,x):
		i = str(self.aligner.iter)
		n = str(self.aligner.cps_counter)
		b = str(min(self.aligner.energies))
		w = str(max(self.aligner.energies[1:]))
		c = str(self.aligner.energies[-1])
		if self.aligner.verbose: 
			sys.stdout.write("fEvals: {}  CPS: {}  Energy: {}  Best: {}  Worst: {}\r".format(i,n,c,b,w))
			sys.stdout.flush()
		self.aligner.write_cps(fname="_pws_coherent_simplex_callback.hdf",n=self.counter)
		self.counter += 1

if __name__ == "__main__":
	main()
