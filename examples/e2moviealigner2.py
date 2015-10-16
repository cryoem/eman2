#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')
from Simplex import Simplex
from EMAN2 import *
import numpy as np
import matplotlib.pyplot as plt
try:
	from scipy.optimize import minimize
	#from scipy.optimize import differential_evolution
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
	parser.add_argument("--maxshift", type=int, help="Set the maximum frame translation distance in pixels.",default=5)
	parser.add_argument("--maxiter", type=int, help="Set the maximum iterations for optimization.",default=500)
	parser.add_argument("--step",type=str,default="0,1",help="Specify <first>,<step>,[last]. Processes only a subset of the input data. ie- 0,2 would process all even particles. Same step used for all input files. [last] is exclusive. Default= 0,1 (first image skipped)")
	parser.add_argument("--fixaxes",action="store_true",default=True,help="Tries to identify bad pixels and fill them in with sane values instead")
	parser.add_argument("--fixbadlines",action="store_true",default=False,help="If you wish to remove detector-specific bad lines, you must specify this flag and --xybadlines.")
	parser.add_argument('--xybadlines', help="Specify the list of bad pixel coordinates for your detector. Will only be used if --fixbadlines is also specified.", nargs=2, default=['3106,3093','3621,3142','4719,3494'])
	parser.add_argument("--show", action="store_true",help="Show average of movie frames before and after alignment.",default=False)
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
		
		if options.gain or options.dark or options.gaink2 or options.fixbadlines:
			if options.verbose: print("Correcting movie frames")
			fname = MovieAligner.correct_frames(options,fname)
		
		if options.verbose: print("Acquiring alignment parameters")
		alignment = MovieAligner(fname,bs=options.boxsize,fixaxes=options.fixaxes)
		if options.verbose: print("Optimizing frame alignment")
		alignment.optimize(options)
		
		if options.show:
			br = alignment.before
			bc = alignment.before.do_fft()
			ar = alignment.after
			ac = alignment.after.do_fft()
			display([br,bc,ar,ac])
	
	E2end(pid)

class MovieAligner:
	
	def __init__(self, fname, bs,fixaxes=False):
		self.path = fname
		self.bs = bs
		self.fixaxes = fixaxes
		self.hdr = EMData(fname,0,True)
		self.hdr['nimg'] = EMUtil.get_image_count(fname)
		self.base = fname.split('.')[0]
		self.frames = []
		for i in xrange(self.hdr['nimg']):
			f = EMData(fname,i)
			if fixaxes: f.process_inplace('filter.xyaxes0',{'x':0,'y':0})
			self.frames.append(f)
		self.translations = np.zeros((self.hdr['nimg'],2))
		self.best_translations = np.zeros((self.hdr['nimg'],2))
		self.before = self.save(descriptor="unaligned",save_frames=False) 
		border = bs+50
		mx = np.arange(border,self.hdr['nx']-border,bs)
		my = np.arange(border,self.hdr['ny']-border,bs)
		self._regions = {}
		for i in xrange(self.hdr['nimg']): 
			self._regions[i] = np.array([[x,y] for y in my for x in mx])
		self._stacks = {}
		for ir in xrange(len(self._regions[0])):
			self._stacks[ir] = [self._regions[i][ir] for i in xrange(self.hdr['nimg'])]
		self.iter = 0
		self.static_fnum = int(self.hdr['nimg']/3)
		self.energies = [np.inf]
		self.verbose = False
		self.cps_counter = 0
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
		self.ips["is_intensity"]=0
		#self.ips.process_inplace('filter.highpass.gauss',{'cutoff_freq':0.0005})
		self.ips.process_inplace('math.rotationalaverage')
		self.real_ips = self.ips.do_ift()
		self.real_ips['is_intensity'] = 0
		self.oned_ips, self.ips_ctf, self.ips_ctf_fit  = self.fit_defocus(self.ips)
	
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
		self.cps["is_intensity"]=0
		#self.cps.process_inplace('filter.highpass.gauss',{'cutoff_freq':0.0005})
		self.real_cps = self.cps.do_ift()
		self.real_cps['is_intensity'] = 0
		self.oned_cps, self.cps_ctf, self.cps_ctf_fit = self.fit_defocus(self.cps)
	
	def calc_energy(self):
		if self.iter == 0:
			self.calc_incoherent_pws()
			self.write_ips()
		self.calc_coherent_pws()
		ips_ctf_fit = np.asarray(self.ips_ctf_fit)
		oned_cps = np.asarray(self.oned_cps)
		# normalize
		ips_ctf_fit/=np.sqrt(ips_ctf_fit.dot(ips_ctf_fit))
		oned_cps/=np.sqrt(oned_cps.dot(oned_cps))
		# regularize
		ips_ctf_fit -= 0.1
		# compare
		compared = np.dot(ips_ctf_fit,oned_cps)
		# scale
		energy = compared #np.log(1-compared)
		#c = self.cps.process('normalize.unitlen')
		#i = self.ips.process('normalize.unitlen')
		#energy = -EMData.cmp(c,'dot',i)
		if energy < min(self.energies):
			self.write_cps()
			self.best_translations = self.translations
		self.energies.append(energy)
		self.iter += 1
		return energy
	
	def write_cps(self):
		self.real_cps.write_image(self.path[:-4]+"_pws_coherent.hdf",self.cps_counter)
		self.cps_counter += 1
	
	def write_ips(self):
		self.real_ips.write_image(self.path[:-4]+"_pws_incoherent.hdf",0)
	
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
	
	def optimize(self, options,bounds=None):
		if options.verbose: self.verbose = True
		state = np.random.randint(-3,3,size=(self.hdr['nimg'],2)).flatten()
		#sm = Simplex(self._compares,state,[1]*len(state),data=self)
		#minimum, error, iters = sm.minimize(0.01,options.maxiter,monitor=1)
		res = minimize(self._compares, state, method='Nelder-Mead', options={'maxiter':options.maxiters,'disp': True,'xtol':0.25}, args=self)
		# if bounds == None:
# 			ms = options.maxshift
# 			c1 = self.static_fnum
# 			c2 = self.hdr['nimg'] - self.static_fnum
# 			ub1 = [round(ms*((c1+1)-(i+1))/(c1+1),1) for i in range(c1)]
# 			ub2 = [round(ms*(i+1)/(c2+1),1) for i in range(c2)]
# 			bounds = []
# 			for b in ub1+ub2:
# 				bounds.append((-int(round(b)),int(round(b))))
# 				bounds.append((-int(round(b)),int(round(b))))
# 		#if options.verbose > 6:
# 		#	print("Initializing optimization with the following bounds:")
# 		#	print("Frame\tLower\t\tUpper")
# 		#	bds = np.array(bounds).reshape((self.hdr['nimg'],2,2)).astype(int)
# 		#	for i,bd in enumerate(bds):
# 		#		print("{}\t( {}, {} )\t( {}, {} )".format(i+1,bd[0,0],bd[1,0],bd[0,1],bd[1,1]))
# 		res = differential_evolution(self._compares, bounds, args=(self,), polish=True, maxiter=options.maxiter, popsize=25, strategy='best2bin', mutation=(0.5,1), recombination=0.7, disp=True, tol=0.01)
		if options.verbose > 6: print(res.message)
		info = "\nEnergy: {}\nIters: {}\nFunc Evals: {}".format(res.fun,res.nit,res.nfev)
		print(info)
		print("\nFrame\tTranslation")
		with open(self.path[:-4]+"_results.txt",'w') as results:
			results.write(info+"\n")
			for i,t in enumerate(self.best_translations):
				info = "{}\t( {}, {} )".format(i+1,int(t[0]),int(t[1]))
				print(info); results.write(info+"\n")
		self.after = self.save(descriptor="aligned",save_frames=True)
	
	@staticmethod
	def _compares(ts,aligner):
		translations = np.asarray(ts).reshape((aligner.hdr['nimg'],2))
		translations = np.round(translations).astype(int)
		for frame_num,trans in enumerate(translations):
			if frame_num != aligner.static_fnum: # keep one frame in place
				aligner.translations[frame_num] = trans.astype(int) # only allow integral shifts
				for reg_num in range(len(aligner._regions[frame_num])): 
					#new_coords = trans #np.add(current_coords,trans)
					aligner._regions[frame_num][reg_num] = trans
					aligner._stacks[reg_num][frame_num] = trans
		energy = aligner.calc_energy()
		if aligner.verbose:
			i = str(aligner.iter)
			b = str(min(aligner.energies))
			c = str(aligner.energies[-1])
			if i > 1: w = str(max(aligner.energies[1:])).ljust(4)
			else: w = l
			out = "{}    Energy: {}    Best: {}    Worst: {}\r".format(i,c,b,w)
			sys.stdout.write(out)
			sys.stdout.flush()
		return energy
	
	def save(self,descriptor,save_frames=False):
		frame_file = self.path[:-4]+"_{}_frames.hdf".format(descriptor)
		average_file = self.path[:-4]+"_{}_average.hdf".format(descriptor)
		avg = Averagers.get('mean')
		for i,t in enumerate(self.translations):
			f = self.frames[i].copy()
			#tf = Transform({'type':'eman','tx':int(round(t[0])),'ty':int(round(t[1]))})
			#f.transform(tf)
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
		step = options.step.split(",")
		if len(step) == 3: last = int(step[2])
		else: last = nd
		first = int(step[0])
		step  = int(step[1])
		if options.xybadlines:
			options.xybadlines = [map(int,s.split(',')) for s in options.xybadlines]
		if not outfile: outfile = options.path[:-4] + "_corrected.hdf"
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
			if not "corrected" in fname:
				im.process_inplace("threshold.clampminmax",{"minval":0,"maxval":im["mean"]+im["sigma"]*3.5,"tozero":1})
				# fixes clear outliers as well as values which were exactly zero
				im.process_inplace("threshold.outlier.localmean",{"sigma":3.5,"fix_zero":1}) 
			#if options.fixaxes: im.process_inplace('filter.xyaxes0')
			if options.fixbadlines:
				for bl in options.xybadlines:
					# remove detector-specific bad lines
					im.process_inplace('math.xybadline',{'xloc':bl[0],'yloc':bl[1]}) 
			im.process_inplace("normalize.edgemean")
			im.write_image(outfile,i-first)
		else: outfile = options.path
		return outfile

if __name__ == "__main__":
	main()
