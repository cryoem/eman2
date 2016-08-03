#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import os
from EMAN2 import *
import subprocess
import shutil
#import pickle
import time
import multiprocessing

colors = ['b', 'g', 'r', 'c', 'm', 'y']

def which(prog):
    cmd = "which {}".format(prog)
    process = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE)
    return process.communicate()[0].replace("\n","")

global pkgs
pkgs = {"EMAN2":"movie_ccf.py",
        "UCSF":"dosefgpu_driftcorr",
        "UNBLUR":"unblur_openmp_7_17_15.exe",
        "DE":"DE_process_frames-2.8.1.py",
        "IMOD":"alignframes",
        "LMBGFS":"alignframes_lmbfgs.exe"}

for pkg in pkgs.keys():
    path=which(pkgs[pkg])
    if path == "": print("Could not find {}.".format(pkgs[pkg]))
    pkgs[pkg] = path

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options]

	Program to characterize DDD motion correction algorithms.
	"""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_argument("--box", default=2048, type=int, help="Size of box to use for sub-region selection.")
	parser.add_argument("--apix", default=None, type=float, help="Apix of input ddd frames. Will search the header by default.")
	parser.add_argument("--skipalign",action="store_true",default=False,help="If you wish to skip running alignments, specify this option.")
	parser.add_argument("--hcreg", type=str, help="Center coordinates of high contrast region. If not specified, region will be computed automatically. Format should be: x,y", default=None)
	parser.add_argument("--lcreg", type=str, help="Center coordinates of low contrast region. If not specified, region will be computed automatically. Format should be: x,y", default=None)
	parser.add_argument("--show", action="store_true",help="Show average of movie frames before and after alignment.",default=True)
	parser.add_argument("--plot",action="store_true",default=False,help="Plot the 1D power spectra and exit.")
	parser.add_argument("--include", type=str, help="Comma separated list of packages to include during comparison.", default="DE,EMAN2,IMOD,UCSF,UNBLUR")
	parser.add_argument("--exclude", type=str, help="Comma separated list of packages to exclude during comparison.", default="")
	# DDD frame correction
	parser.add_argument("--dark",type=str,default=None,help="Perform dark image correction using the specified image file")
	parser.add_argument("--gain",type=str,default=None,help="Perform gain image correction using the specified image file")
	parser.add_argument("--gaink2",type=str,default=None,help="Perform gain image correction. Gatan K2 gain images are the reciprocal of DDD gain images.")
	parser.add_argument("--fixaxes",action="store_true",default=False,help="Tries to identify bad pixels and fill them in with sane values instead")
	parser.add_argument("--neighbornorm", type=int, help="Set the norm to be used for fixing axes. Default is 2",default=2)
	parser.add_argument("--fixbadlines",action="store_true",default=False,help="If you wish to remove detector-specific bad lines, you must specify this flag and --xybadlines.")
	parser.add_argument('--xybadlines', help="Specify the list of bad pixel coordinates for your detector. Will only be used if --fixbadlines is also specified.", nargs=2, default=['3106,3093','3621,3142','4719,3494'])
	# program options
	parser.add_argument("--threads", default=multiprocessing.cpu_count()-2, type=int, help="Number of threads to use with each aligner. By default will use all but 2 cpus.")
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-2)

	global options
	(options, args) = parser.parse_args()

	if len(args) < 1:
		print usage
		parser.error("You must specify a single DDD movie stack in HDF or MRC format.")
		sys.exit(1)

	options.include = [o.upper() for o in options.include.split(",") if o not in options.exclude.split(",") if pkgs[o] != ""]

	included = [pkg for pkg in pkgs if pkg in options.include]
	if len(included) < 1:
		print("You must specify at least one package to run this program (DE, EMAN2, IMOD, UCSF, UNBLUR)")
		sys.exit(1)

	if options.apix:
		apix = options.apix
	else:
		try: apix = hdr["apix_x"]
		except: apix = 1.0

	bs = options.box
	logid=E2init(sys.argv,options.ppid)

	for arg in args:

		fname = arg
		hdr = EMData(fname,0,True)

		if options.gain or options.dark or options.gaink2:
			if options.verbose: print("Correcting frames before processing")
			fname = FrameCorrector.correct_frames(options,fname)

		if options.verbose: print("Processing {}".format(fname))

		# PART 0: Setup alignment data
		(basen,ext) = os.path.splitext(fname)
		bdir = "{}".format(basen)

		try: os.makedirs(bdir)
		except: pass

		hcname = "{}/hictrst.hdf".format(bdir)
		lcname = "{}/loctrst.hdf".format(bdir)

		frames = load_frames(fname)
		avg = average_frames(frames)
		avg.write_image("{}/{}_avg.hdf".format(bdir,basen))

		if options.hcreg == None or options.lcreg == None:
			hcreg, lcreg = get_hclc_regions(avg,n=64,bs=bs,edge=10)
			hcx,hcy,_= hcreg.get_origin()
			lcx,lcy,_= lcreg.get_origin()
		else:
			hcreg = get_region(options.hcreg,bs)
			lcreg = get_region(options.lcreg,bs)
			hcx,hcy,_ = hcreg.get_origin()
			lcx,lxy,_ = lcreg.get_origin()

		with open("{}/params.txt".format(bdir),"w") as p:
			p.write("# boxsize\thigh contrast coordinates (x,y)\tlow contrast coordinates (x,y)\n")
			p.write("{}\t{},{}\t{},{}\n".format(bs,hcx,hcy,lcx,lcy))

		# write regions to disk (add large enough edge to region so we can clip down to desired box size after)
		hcrs = []
		lcrs = []
		for i,f in enumerate(frames):
			hcf = f.get_clip(hcreg)
			hcf.write_image(hcname,i)
			hcrs.append(hcf)
			lcf = f.get_clip(lcreg)
			lcf.write_image(lcname,i)
			lcrs.append(lcf)
		hcrs_avg = average_frames(hcrs)
		lcrs_avg = average_frames(lcrs)
		hcrs_avg.write_image("{}/hictrst_noali_avg.hdf".format(bdir),0)
		lcrs_avg.write_image("{}/loctrst_noali_avg.hdf".format(bdir),0)

		trans_wmg = {} # whole micrograph
		trans_lo = {} # low contrast region
		trans_hi = {} # high contrast region

		# process low and high contrast regions
		cwd = os.getcwd()
		localfn = fname.split("/")[-1]
		localhc = hcname.split("/")[-1]
		locallc = lcname.split("/")[-1]

		runtimes = {key:[] for key in pkgs.keys()}

		for pkg in sorted(pkgs.keys()):
			if pkg in options.include:
				prog = pkgs[pkg].split("/")[-1]
				if not options.skipalign: print("Running {} on {}".format(prog,fname))
				else: print("Analyzing {} alignment".format(prog,fname))

				pdir = "{}/{}".format(bdir,pkg)
				try: os.makedirs(pdir)
				except: pass

				if pkg == "DE":

					if not options.skipalign:
						src = "/".join(pkgs[pkg].split("/")[:-1]) + "/DE_process_frames.cfg"
						dst = cwd+"/DE_process_frames.cfg"
						try: os.symlink(src,dst)
						except: pass

						for fn in [fname,hcname,lcname]: # whole micrograph, high contrast, low contrast
							o,e,rt=run("{} {} {} --run_cores {}".format(pkgs[pkg],pdir,fn,options.threads))
							runtimes[pkg].append(rt)

					for f in os.listdir(pdir):
						if "translations_x" in f:
							if "hictrst" in f: hi_x = os.path.join(pdir,f)
							elif "loctrst" in f: lo_x = os.path.join(pdir,f)
							else: wmg_x = os.path.join(pdir,f)
						elif "translations_y" in f:
							if "hictrst" in f: hi_y = os.path.join(pdir,f)
							elif "loctrst" in f: lo_y = os.path.join(pdir,f)
							else: wmg_y = os.path.join(pdir,f)

					try:
						trans_wmg[pkg] = parse_de(wmg_x,wmg_y)
						trans_hi[pkg] = parse_de(hi_x,hi_y)
						trans_lo[pkg] = parse_de(lo_x,lo_y)
					except: failed(pkg)

				if pkg == "EMAN2":

					if not options.skipalign:
						for fn in [fname,hcname,lcname]: # whole micrograph, high contrast, low contrast
							lfn = fn.split("/")[-1] # local file name
							if not os.path.isfile("{}/{}".format(pdir,lfn)): shutil.copy2(fn,"{}/{}".format(pdir,lfn))

							o,e,rt=run("{} {} {}".format(pkgs[pkg],lfn,options.threads),cwd=pdir)
							runtimes[pkg].append(rt)

					for f in os.listdir(pdir):
						if "hictrst_info.txt" in f: hi = os.path.join(pdir,f)
						elif "loctrst_info.txt" in f: lo = os.path.join(pdir,f)
						elif "info.txt" in f: wmg = os.path.join(pdir,f)

					try:
						trans_wmg[pkg] = parse_eman2(wmg)
						trans_hi[pkg] = parse_eman2(hi)
						trans_lo[pkg] = parse_eman2(lo)
					except: failed(pkg)

				if pkg == "IMOD":

					if not options.skipalign:
						for fn in [fname,hcname,lcname]: # whole micrograph, high contrast, low contrast
							lfn = fn.split("/")[-1] # local file name
							if not os.path.isfile("{}/{}".format(pdir,lfn)): shutil.copy2(fn,"{}/{}".format(pdir,lfn))

							o,e,rt=run("{} -gpu 0 -n -xfext xf -input {}".format(pkgs[pkg],lfn),cwd=pdir)
							runtimes[pkg].append(rt)

					for f in os.listdir(pdir):
						if "hictrst.xf" in f: hi = os.path.join(pdir,f)
						elif "loctrst.xf" in f: lo = os.path.join(pdir,f)
						elif ".xf" in f: wmg = os.path.join(pdir,f)

					try:
						trans_wmg[pkg] = parse_imod(wmg)
						trans_hi[pkg] = parse_imod(hi)
						trans_lo[pkg] = parse_imod(lo)
					except: failed(pkg)

				if pkg == "UCSF":

					if not options.skipalign:
						for fn in [fname,hcname,lcname]: # whole micrograph, high contrast, low contrast
							lfn = fn.split("/")[-1]
							lext = lfn.split(".")[-1]
							if lext == "hdf": # HDF format not supported by UCSF
								shutil.copy2(fn,"{}/{}".format(pdir,lfn))
								lfn_new = lfn.replace(".hdf",".mrcs")
								if not os.path.isfile("{}/{}".format(pdir,lfn_new)):
									o,e,rt=run("e2proc2d.py {} {}".format(lfn,lfn_new),cwd=pdir)
									os.remove("{}/{}".format(pdir,lfn))
								lfn = lfn_new
								lext = "mrcs"
							else:
								if not os.path.isfile("{}/{}".format(pdir,lfn)): shutil.copy2(fn,"{}/{}".format(pdir,lfn))

							o,e,rt=run("{} {} -srs 1 -ssc 1 -atm 1".format(pkgs[pkg],lfn),cwd=pdir)
							runtimes[pkg].append(rt)

					for f in os.listdir(pdir):
						if "hictrst_Log.txt" in f: hi = os.path.join(pdir,f)
						elif "loctrst_Log.txt" in f: lo = os.path.join(pdir,f)
						elif "_Log.txt" in f: wmg = os.path.join(pdir,f)

					try:
						trans_wmg[pkg] = parse_ucsf(wmg)
						trans_hi[pkg] = parse_ucsf(hi)
						trans_lo[pkg] = parse_ucsf(lo)
					except: failed(pkg)

				if pkg == "UNBLUR":

					if not options.skipalign:
						template = """{prog} <<EOF
{fname}
{nframes}
{alignsum}
{shiftfile}
{apix}
{dosefilt}
{saveali}
{aliname}
{advopts}
{frcname}
{minsrch}
{maxsrch}
{bfact}
{vfmwidth}
{hfmwidth}
{thresh}
{maxiter}
{verbose}
EOF
"""

						for fn in [fname,hcname,lcname]:
							lfn = fn.split("/")[-1]
							lext = lfn.split(".")[-1]

							if lext == "hdf": # HDF format not supported by unblur
								shutil.copy2(fn,"{}/{}".format(pdir,lfn))
								lfn_new = lfn.replace(".hdf",".mrcs")
								if not os.path.isfile("{}/{}".format(pdir,lfn_new)):
									o,e,rt=run("e2proc2d.py {} {}".format(lfn,lfn_new),cwd=pdir)
									os.remove("{}/{}".format(pdir,lfn))
								lfn = lfn_new
								lext = "mrcs"
							else:
								if not os.path.isfile("{}/{}".format(pdir,lfn)): shutil.copy2(fn,"{}/{}".format(pdir,lfn))

							if lext == "mrcs": # MRCS format not supported by unblur
								lfn_new = lfn.replace(".mrcs",".mrc")
								try: shutil.move("{}/{}".format(pdir,lfn),"{}/{}".format(pdir,lfn_new))
								except: pass # file likely already exists
								lfn = lfn_new
								lext = "mrc"

							nfs = len(frames)
							alis =lfn.replace(".mrc","_ali_sum.mrc")
							shft =lfn.replace(".mrc","_shifts.txt")
							ali=lfn.replace(".mrc","_ali.mrc")
							frc=lfn.replace(".mrc","_frc.txt")

							cmd = template.format(prog=pkgs[pkg],fname=lfn, nframes=nfs, alignsum=alis, shiftfile=shft, apix=apix, dosefilt="NO", saveali="YES", aliname=ali, advopts="NO", frcname=frc, minsrch=2.0, maxsrch=200.0, bfact=1500, vfmwidth=1, hfmwidth=1, thresh=0.1, maxiter=10, verbose="NO")

							o,e,rt=run(cmd,shell=True,cwd=pdir)
							runtimes[pkg].append(rt)

					for f in os.listdir(pdir):
						if "hictrst_shifts.txt" in f: hi = os.path.join(pdir,f)
						elif "loctrst_shifts.txt" in f: lo = os.path.join(pdir,f)
						elif "_shifts.txt" in f: wmg = os.path.join(pdir,f)

					try:
						trans_wmg[pkg] = parse_unblur(wmg)
						trans_hi[pkg] = parse_unblur(hi)
						trans_lo[pkg] = parse_unblur(lo)
					except: failed(pkg)

				if pkg == "LMBFGS":

					if not options.skipalign:
						template = """{prog} << EOF
{movielist}
{boxsize},{psize},{nsigma},{rmax1},{rmax2},{smooth}
{bfactor}
{framefirst},{framelast},{zeroframe}
{factr}
{inpath}
{outpath}
shf
vec
EOF
"""

						for fn in [fname,hcname,lcname]:
							lfn = fn.split("/")[-1]
							lext = lfn.split(".")[-1]

							mlfn = "{}/movielist.txt".format(pdir)
							with open(mlfn,"w") as movielist:
								movielist.write("{}".format(lfn))

							if lext == "hdf": # HDF format not supported by unblur
								shutil.copy2(fn,"{}/{}".format(pdir,lfn))
								lfn_new = lfn.replace(".hdf",".mrcs")
								if not os.path.isfile("{}/{}".format(pdir,lfn_new)):
									o,e,rt=run("e2proc2d.py {} {}".format(lfn,lfn_new),cwd=pdir)
									os.remove("{}/{}".format(pdir,lfn))
								lfn = lfn_new
								lext = "mrcs"
							else:
								if not os.path.isfile("{}/{}".format(pdir,lfn)): shutil.copy2(fn,"{}/{}".format(pdir,lfn))

							#bs = 3600 # defined earlier
							ps = 1.45
							nsig = 5
							rm1 = 500
							rm2 =100
							smooth = "0d0"
							bfact = 2000
							ffirst = 1
							flast = 0
							fmiddle = int(len(frames)/2)
							factr = "1d1"
							inpath = pdir
							outpath = pdir

							cmd = template.format(movielist=mlfn, boxsize=bs, psize=ps, nsigma=nsig, rmax1=rm1, rmax2=rm2, smooth=smooth, bfactor=bfact, framefirst=ffirst, framelast=flast, zeroframe=fmiddle, factr=factr, inpath=inpath, outpath=outpath)

							o,e,rt=run(cmd,shell=True,cwd=pdir,exe="/bin/csh")
							runtimes[pkg].append(rt)

					for f in os.listdir(pdir):
						if "hictrst.shf" in f: hi = os.path.join(pdir,f)
						elif "loctrst.shf" in f: lo = os.path.join(pdir,f)
						elif ".shf" in f: wmg = os.path.join(pdir,f)

					try:
						trans_wmg[pkg] = parse_unblur(wmg)
						trans_hi[pkg] = parse_unblur(hi)
						trans_lo[pkg] = parse_unblur(lo)
					except: failed(pkg)

		# PART 1: How quickly do these alorithms run
		if not options.skipalign:
			with open("{}/runtimes.txt".format(bdir),"w") as f:
				f.write("PKG\tRUNTIME\n")
				for pkg in options.include:
					q,r,s = runtimes[pkg]
					#print(pkg,runtimes)
					f.write("{}\t{}\n".format(pkg,q,r,s))

		# PART 2: How does calculated motion differ between the most commonly used alignment algorithms?
		trans_wmg, mags_wmg, fig1 = plot_trans(trans_wmg,bdir)
		if options.plot: plt.show()

		#pkl = "{}/trans_wmg.p".format(bdir)
		#pickle.dump(trans_wmg,open(pkl,"wb"))
		#pkl = "{}/mags_wmg.p".format(bdir)
		#pickle.dump(mags_wmg,open(pkl,"wb"))

		# PART 3: How different are predicted frame translations with and without high contrast features (real space comparison)
		trans_hi, trans_lo, ssdfs, fig2 = plot_differences(trans_hi,trans_lo,bdir)
		if options.plot: plt.show()

		#pkl = "{}/trans_hi.p".format(bdir)
		#pickle.dump(trans_hi,open(pkl,"wb"))
		#pkl = "{}/trans_lo.p".format(bdir)
		#pickle.dump(trans_lo,open(pkl,"wb"))

		# Part 4: Compare agreement in fourier space (coherent/incoherent power spectra)
		frames_hi = load_frames(hcname)
		frames_lo = load_frames(lcname)

		ftypes = {} # need to get regions from original micrographs to avoid edge effects.
		ftypes["hi_wmg"] = shift_by(frames_hi,trans_wmg)
		ftypes["lo_wmg"] = shift_by(frames_lo,trans_wmg)
		ftypes["lo_lo"] = shift_by(frames_lo,trans_lo)
		ftypes["lo_hi"] = shift_by(frames_lo,trans_hi)
		ftypes["hi_lo"] = shift_by(frames_hi,trans_lo)
		ftypes["hi_hi"] = shift_by(frames_hi,trans_hi)
		scores = calc_cips_scores(ftypes)

		#pkl = "{}/scores.p".format(bdir)
		#pickle.dump(scores,open(pkl,"wb"))

		fig3 = plot_scores(scores,bdir)
		if options.plot: plt.show()

	print("DONE")

def failed(pkg):
	print("Could not find frame shifts for {} aligner. Perhaps you should use the --align option.".format(pkg))
	sys.exit(1)

def get_region(c,bs):
	x,y=map(int,c.split(","))
	return Region(x-bs/2,y-bs/2,bs,bs)

def get_hclc_regions(img,n=100,bs=2048,edge=128):
    rbs = int(np.sqrt(bs))
    idx = np.argmin([(gbs-rbs)**2 for gbs in good_box_sizes])
    nbs = good_box_sizes[idx]
    xmin = bs/2 + edge
    xmax = img["nx"] - xmin
    ymin = bs/2 + edge
    ymax = img["ny"] - ymin
    coords = np.random.randint(low=max(xmin,ymin),high=min(xmax,ymax),size=(n,2))
    crds = []
    regs = []
    energy = []
    for i,(x,y) in enumerate(coords):
        r = Region(x-bs/2,y-bs/2,bs,bs)
        reg = img.get_clip(r)
        reg.process_inplace("math.fft.resample",{"n":(bs/nbs)})
        reg.process_inplace("filter.highpass.gauss",{"cutoff_pixels":5})
        reg.process_inplace("filter.lowpass.gauss",{"cutoff_abs":0.1})
        reg.process_inplace("mask.gaussian",{"inner_radius":nbs/6,"outer_radius":nbs/3})
        regs.append(reg)
        energy.append(reg["sigma"]*(reg["maximum"]-reg["minimum"])**2)
        crds.append("{},{}".format(x,y))
    regs.sort(key=dict(zip(regs, energy)).get)
    crds.sort(key=dict(zip(crds, energy)).get)
    hc = map(int,crds[-1].split(","))
    lc = map(int,crds[0].split(","))
    hcr = Region(hc[0]-bs/2,lc[1]-bs/2,bs,bs)
    lcr = Region(lc[0]-bs/2,lc[1]-bs/2,bs,bs)
    return hcr,lcr

def run(cmd,shell=False,cwd=None,exe="/bin/sh"):
	if options.verbose: print(cmd.replace("\n"," "))
	if cwd == None:
		cwd = os.getcwd()
	if shell == False:
		cmd = cmd.split()
	process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE, shell=shell, cwd=cwd, executable=exe)
	start = time.time()
	out, err = process.communicate()
	runtime = time.time() - start
	if options.verbose: print("Runtime: {}".format(runtime))
	return out, err, runtime

def shift_by(frames,trans):
	d = {}
	for key in trans.keys():
		if key != "ALL":
			d[key] = []
			for i, frame in enumerate(frames):
				f = frame.copy()
				f.translate(trans[key][i][0],trans[key][i][1],0)
				d[key].append(f)
	return d

def calc_cips_scores(ftypes,bs=512):
	scores = {}
	for pkg in [pkg for pkg in pkgs.keys() if pkg in options.include]:
		scores[pkg] = []
		for ftype in sorted(ftypes.keys()):
			cps = calc_coherent_pws(ftypes[ftype][pkg],bs)
			ips = calc_incoherent_pws(ftypes[ftype][pkg],bs)
			frc = -1 * EMData.cmp(cps,'frc',ips)
			dot = -1 * EMData.cmp(cps,'dot',ips)
			osb = EMData.cmp(cps,'optsub',ips)
			sse = EMData.cmp(cps,'sqeuclidean',ips)
			scores[pkg].append([frc,dot,osb,sse])
		scores[pkg] = np.asarray(scores[pkg])
	return scores

def plot_scores(scores,bdir): # 1 axis per metric
	metrics = np.asarray(["frc","dot","optsub","sqeuclidean"]).reshape(2,2)
	fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(16,16))
	for (i,j),metric in np.ndenumerate(metrics):
		ax = axs[i][j]
		for key in scores.keys():
			if metric == "frc":
				ax.plot(scores[key][:,0],label=key)
			elif metric == "dot":
				ax.plot(scores[key][:,1],label=key)
			elif metric == "optsub":
				ax.plot(scores[key][:,2],label=key)
			elif metric == "sqeuclidean":
				ax.plot(scores[key][:,3],label=key)
		ax.legend(loc="best")
		ax.set_xticks([0,1,2,3,4,5])
		ax.set_xticklabels(["hi_hi", "hi_lo", "hi_wmg", "lo_hi","lo_lo","lo_wmg"])
		ax.set_xlabel("Aligned To")
		ax.set_title(metric)
	fig.tight_layout()
	plt.savefig("{}/scores.png".format(bdir))
	return fig

def plot_trans(trans,bdir,nsig=1):
	ats = np.hstack([trans[k] for k in sorted(trans.keys())])
	np.savetxt("{}/all_trans.txt".format(bdir),ats)

	trans["ALL"] = np.dstack([trans[key] for key in trans.keys()])
	trans["MEAN"] = np.mean(trans["ALL"],axis=2)
	trans["VAR"] = np.var(trans["ALL"],axis=2)
	trans["STD"] = np.std(trans["ALL"],axis=2)
	err = nsig * trans["STD"]
	if nsig == 1: siglabel = "STD"
	else: siglabel = "{}x STD".format(nsig)
	exclude = ["STD","VAR","ALL"]
	xx = np.arange(1,len(trans["MEAN"][:,0])+1,1)

	mags = {"MEAN":np.sqrt(trans["MEAN"][:,0]**2+trans["MEAN"][:,1]**2)}
	for key in trans.keys():
		if key not in exclude + ["MEAN"]:
			mags[key] = np.sqrt(trans[key][:,0]**2+trans[key][:,1]**2) - mags["MEAN"]

	fig = plt.figure()

	ax1 = fig.add_subplot(211)

	cctr = 0
	for key in trans.keys():
		if key not in exclude:
			if key != "MEAN":
				ax1.plot(trans[key][:,0],trans[key][:,1],colors[cctr],label=key,alpha=0.8)
			else:
				ax1.plot(trans[key][:,0],trans[key][:,1],'k-',linewidth=1.5,label=key,alpha=1)
				ax1.errorbar(trans[key][:,0],trans[key][:,1],xerr=err[:,0],yerr=err[:,1],color='k',alpha=0.6,label=siglabel)
			cctr += 1
	ax1.legend(loc="best")#,bbox_to_anchor=(0.95, 0.45))
	ax1.set_title("Measured Drift")
	ax1.set_xlabel("X Frame Translation (pixels)")
	ax1.set_ylabel("Y Frame Translations (pixels)")

	ax2 = fig.add_subplot(212)

	cctr = 0
	for key in mags.keys():
		if key != "MEAN":
			ax2.plot(xx,mags[key],colors[cctr],label=key,alpha=1)
		else:
			ax2.plot(xx,np.zeros_like(xx),"k--",label=key,alpha=1)
			ax2.errorbar(xx,np.zeros_like(xx),xerr=err[:,0],yerr=err[:,1],color='k',alpha=0.6,label=siglabel)
		cctr+=1
	ax2.set_xlabel("Frame Number (#)")
	ax2.set_ylabel("Relative shift magnitude (pixels)")
	ax2.set_title("Frame shift magnitude (relative to MEAN)")
	ax2.legend(loc="best")

	fig.tight_layout()

	plt.savefig("{}/trans.png".format(bdir))
	return trans, mags, fig

def plot_differences(trans_hi,trans_lo,bdir,nsig=1):
	hts = np.hstack([trans_hi[k] for k in sorted(trans_hi.keys())])
	lts = np.hstack([trans_lo[k] for k in sorted(trans_lo.keys())])
	np.savetxt("{}/hi_trans.txt".format(bdir),hts)
	np.savetxt("{}/lo_trans.txt".format(bdir),lts)

	trans_hi["ALL"] = np.dstack([trans_hi[key] for key in trans_hi.keys()])
	trans_hi["MEAN"] = np.mean(trans_hi["ALL"],axis=2)
	trans_hi["VAR"] = np.var(trans_hi["ALL"],axis=2)
	trans_hi["STD"] = np.std(trans_hi["ALL"],axis=2)
	err_hi = nsig * trans_hi["STD"]
	trans_lo["ALL"] = np.dstack([trans_lo[key] for key in trans_lo.keys()])
	trans_lo["MEAN"] = np.mean(trans_lo["ALL"],axis=2)
	trans_lo["VAR"] = np.var(trans_lo["ALL"],axis=2)
	trans_lo["STD"] = np.std(trans_lo["ALL"],axis=2)
	err_lo = nsig * trans_lo["STD"]
	if nsig == 1: siglabel = "STD"
	else: siglabel = "{}x STD".format(nsig)
	exclude = ["STD","VAR","ALL"]
	xx = np.arange(1,len(trans_hi["MEAN"][:,0])+1,1)
	mags_hi = {"MEAN":np.sqrt(trans_hi["MEAN"][:,0]**2+trans_hi["MEAN"][:,1]**2)}
	mags_lo = {"MEAN":np.sqrt(trans_lo["MEAN"][:,0]**2+trans_lo["MEAN"][:,1]**2)}
	for key in trans_hi.keys():
		if key not in exclude + ["MEAN"]:
			mags_hi[key] = np.sqrt(trans_hi[key][:,0]**2+trans_hi[key][:,1]**2) - mags_hi["MEAN"]
			mags_lo[key] = np.sqrt(trans_lo[key][:,0]**2+trans_lo[key][:,1]**2) - mags_lo["MEAN"]
	ssdfs = {}

	grid_size = (4,2)
	ax1 = plt.subplot2grid(grid_size,(0,0),rowspan=2,colspan=2)
	cctr = 0
	with open("{}/sse.txt".format(bdir),"w") as sse:
		for key in trans_hi.keys():
			if key not in exclude:
				if key != "MEAN":
					ssdfs[key] = np.linalg.norm(np.sum(np.square((trans_hi[key]-trans_lo[key])),axis=0))
					print("{}\t{}".format(key,ssdfs[key]))
					ax1.plot(trans_hi[key][:,0],trans_hi[key][:,1],"{}-".format(colors[cctr]),label="{} high".format(key))
					ax1.plot(trans_lo[key][:,0],trans_lo[key][:,1],"{}.".format(colors[cctr]),label="{} low".format(key))
				else:
					ssdfs[key] = np.linalg.norm(np.sum(np.square((trans_hi[key]-trans_lo[key])),axis=0))
					print("{}\t{}".format(key,ssdfs[key]))
					ax1.plot(trans_hi[key][:,0],trans_hi[key][:,1],"k--".format(colors[cctr]),linewidth=2,label="{} high".format(key))
					ax1.plot(trans_lo[key][:,0],trans_lo[key][:,1],"k.".format(colors[cctr]),linewidth=2,label="{} low".format(key))
				sse.write("{}\t{}\n".format(key,ssdfs[key]))
				cctr+=1
	ax1.set_xlabel("X-shift (pixels)")
	ax1.set_ylabel("Y-shift (pixels)")
	ax1.set_title("Movie trajectory")
	ax1.legend(loc="best")
	ax2 = plt.subplot2grid(grid_size,(2,0),rowspan=2,colspan=2)
	cctr = 0
	for key in mags_hi.keys():
		if key != "MEAN":
			ax2.plot(xx,mags_hi[key],label="{} high".format(key),alpha=1)
		else:
			#ax2.plot(xx,np.zeros_like(xx),"k--",label=key,alpha=1)
			ax2.errorbar(xx,np.zeros_like(xx),xerr=err_hi[:,0],yerr=err_hi[:,1],color='k',alpha=0.6,label="{} high".format(siglabel),linewidth=2)
	ax2.set_xlabel("Frame Number (#)")
	ax2.set_ylabel("Relative shift magnitude (pixels)")
	ax2.set_title("Frame shift magnitude (relative to MEAN)")
	ax2.legend(loc="best")
	cctr = 0
	for key in mags_lo.keys():
		if key != "MEAN":
			ax2.plot(xx,mags_lo[key],"{}--".format(colors[cctr]),label="{} low".format(key),alpha=1)
		else:
			ax2.plot(xx,np.zeros_like(xx),"k--",label=key,alpha=1)
			ax2.errorbar(xx,np.zeros_like(xx),xerr=err_lo[:,0],yerr=err_lo[:,1],color='k',alpha=0.6,label="{} low".format(siglabel))
		cctr+=1
	ax2.set_xlabel("Frame Number (#)")
	ax2.set_ylabel("Relative shift magnitude (pixels)")
	ax2.set_title("Frame shift magnitude (relative to MEAN)")
	ax2.legend(loc="best")
	fig = plt.gcf()
	fig.tight_layout()
	plt.savefig("{}/differences.png".format(bdir))
	return trans_hi, trans_lo, ssdfs, fig

def parse_eman2(fn):
	with open(fn) as f:
		trans = [txt.split("\t")[1:3] for txt in f.readlines()[1:]]
	trans = np.asarray(trans).astype(float)
	return np.cumsum(trans - np.mean(trans,axis=0),axis=1)

def parse_de(fnx,fny):
	with open(fnx) as xf:
		xs = np.asarray(xf.read().replace("\n","").split("\t"))[1:].astype(float)
	with open(fny) as yf:
		ys = np.asarray(yf.read().replace("\n","").split("\t"))[1:].astype(float)
	trans = np.vstack([xs,ys]).T
	return np.cumsum(trans,axis=1)# - np.mean(trans,axis=0),axis=1)

def parse_imod(fn):
	with open(fn) as f:
		trans = np.asarray([np.asarray(txt.split()[4:6]).astype(float) for txt in f.readlines()])
	return np.cumsum(trans - np.mean(trans,axis=0),axis=1)

def parse_ucsf(fn):
	trans = []
	with open(fn) as f:
		for txt in f.readlines():
			if "......Shift of Frame #" in txt:
				trans.append(np.asarray(txt.replace("......Add Frame #","").split(":")[-1].split())[:2])
	trans = np.asarray(trans).astype(float)
	return np.cumsum(trans - np.mean(trans,axis=0),axis=1)

def parse_unblur(fn):
	lines = [] #[xs,ys]
	with open(fn) as f:
		for txt in f.readlines():
			if "#" not in txt: lines.append(txt.split())
			elif "Pixel size (A):" in txt: apix = float(txt.split()[-1])
	trans = np.asarray([[x,y] for (x,y) in zip(lines[0],lines[1])]).astype(float)
	trans /= apix
	return np.cumsum(trans - np.mean(trans,axis=0),axis=1)

def calc_incoherent_pws(frames,bs=512):
    mx = np.arange(bs+50,frames[0]['nx']-bs+50,bs)
    my = np.arange(bs+50,frames[1]['ny']-bs+50,bs)
    regions = {}
    for i in xrange(len(frames)): regions[i] = [[x,y] for y in my for x in mx]
    ips = Averagers.get('mean')
    for i in xrange(len(frames)):
        img = frames[i].copy()
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
    ips.process_inplace('math.rotationalaverage')
    return ips

def calc_coherent_pws(frames,bs=512):
    mx = np.arange(bs+50,frames[0]['nx']-bs+50,bs)
    my = np.arange(bs+50,frames[1]['ny']-bs+50,bs)
    regions = {}
    for i in xrange(len(frames)):
        regions[i] = [[x,y] for y in my for x in mx]
    stacks = {}
    for ir in xrange(len(regions[0])):
        stacks[ir] = [regions[i][ir] for i in xrange(len(frames))]
    cps = Averagers.get('mean')
    for s in xrange(len(stacks)):
        stack_avg = Averagers.get('mean')
        for i,r in enumerate(stacks[s]):
            stack_avg.add_image(frames[i].copy().get_clip(Region(r[0],r[1],bs,bs)))
        avg = stack_avg.finish()
        avg.process_inplace('normalize.unitlen')
        avg.do_fft_inplace()
        avg.ri2inten()
        cps.add_image(avg)
    cps = cps.finish()
    cps.process_inplace('math.sqrt')
    cps.process_inplace('normalize.edgemean')
    return cps

def load_frames(fn): return [EMData(fn,i) for i in range(EMUtil.get_image_count(fn))]

def shift_frames(frames,trans):
	shifted = []
	for frame,(x,y) in zip(frames,trans):
		f = frame.copy()
		f.translate(x,y,0)
		shifted.append(f)
	return shifted

def average_frames(frames):
	avgr = Averagers.get("mean")
	avgr.add_image_list(frames)
	return avgr.finish()


class FrameCorrector:

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


if __name__ == "__main__":
	main()
