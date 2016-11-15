#!/usr/bin/env python

# Author: Michael Bell 09/2016

import os

if os.getenv("DISPLAY") == None:
	import matplotlib as mpl
	mpl.use('Agg')

import matplotlib.pyplot as plt

from EMAN2 import *
from scipy import ndimage
import numpy as np
import os
import shutil
import time
import multiprocessing
import subprocess
import gc
import psutil

colors = ['b', 'g', 'r', 'c', 'm', 'y']

def which(prog):
	cmd = "which {}".format(prog)
	process = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE)
	return process.communicate()[0].replace("\n","")

global pkgs
pkgs = {"EMAN2":"e2ddd_movie.py",
		"UCSF":"dosefgpu_driftcorr", # motioncor 2.1
		"UNBLUR":"unblur_openmp_7_17_15.exe",
		"DE":"DE_process_frames-2.8.1.py",
		"IMOD":"alignframes",
		"LMBFGS":"alignframes_lmbfgs.exe"}

found = []
for pkg in pkgs.keys():
	path = which(pkgs[pkg])
	pkgs[pkg] = path # replace program name with the output of "which"
	if path != "": found.append(pkg)
print("Found the following alignment packages in your current environment: {}\n".format(", ".join(found)))

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options]

	Program to characterize DDD motion correction algorithms.
	"""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	# Program arguments
	parser.add_argument("--apix", default=None, type=float, help="Apix of input ddd frames. Will search the header by default.")
	parser.add_argument("--skipalign",action="store_true",default=False,help="If you wish to skip running alignments, specify this option.")
	parser.add_argument("--plot",action="store_true",default=False,help="Plot the 1D power spectra and exit.")
	parser.add_argument("--include", type=str, help="Comma separated list of packages to include during comparison.", default="DE,EMAN2,IMOD,UCSF,UNBLUR,LMBFGS")
	parser.add_argument("--exclude", type=str, help="Comma separated list of packages to exclude during comparison.", default="")
	# DDD frame correction
	parser.add_argument("--dark",type=str,default=None,help="Perform dark image correction using the specified image file")
	parser.add_argument("--gain",type=str,default=None,help="Perform gain image correction using the specified image file")
	parser.add_argument("--gaink2",type=str,default=None,help="Perform gain image correction. Gatan K2 gain images are the reciprocal of DDD gain images.")
	parser.add_argument("--step",type=str,default="1,1,31",help="Specify <first>,<step>,[last]. Processes only a subset of the input data. ie- 0,2 would process all even frames. Same step used for all input files. [last] is exclusive. Default= 1,1,31.")
	# Gold subtraction
	parser.add_argument("--goldsize", default=70, type=float, help="Gold diameter in pixels. Default is 70.")
	parser.add_argument("--boxsize", default=128, type=float, help="Box size to use when tiling image to compute noise. Default is 128.")
	parser.add_argument("--oversample", default=6, type=float, help="Oversample by this much when tiling image to compute noise. Higher values increase smoothness of transition between noise from one region to its neighbor, but require significantly more time to process. Default is 6.")
	# Other options
	parser.add_argument("--threads", default=1, type=int, help="Number of threads to use with each aligner. Default is 1.")
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-2)

	global options
	(options, args) = parser.parse_args()

	if len(args) < 1:
		print usage
		parser.error("You must specify a single DDD movie stack in HDF or MRC format.")
		sys.exit(1)

	options.include = [o.upper() for o in options.include.split(",") if o.upper() not in options.exclude.split(",") if pkgs[o] != ""]

	included = [pkg for pkg in pkgs if pkg in options.include]
	if len(included) < 1:
		print("You must specify at least one package to run this program (DE, EMAN2, IMOD, UCSF, UNBLUR)")
		sys.exit(1)

	if options.apix:
		apix = options.apix
	else:
		try: apix = hdr["apix_x"]
		except: apix = 1.0

	if not options.skipalign: print("Performing alignments with the following packages: {}\n".format(", ".join(options.include)))

	logid=E2init(sys.argv,options.ppid)

	for arg in args:

		fname = arg
		hdr = EMData(fname,0,True)

		if options.verbose: print("Processing {}".format(fname))

		# PART 1: Setup alignment data
		(basen,ext) = os.path.splitext(fname)
		bdir = "{}".format(basen)

		try: os.makedirs(bdir)
		except: pass

		hcname = "{}/hictrst.hdf".format(bdir)
		lcname = "{}/loctrst.hdf".format(bdir)

		localhc = hcname.split("/")[-1]
		locallc = lcname.split("/")[-1]

		# Correct Frames

		if options.gain or options.dark or options.gaink2:
			if os.path.isfile("{}/{}_proc_corr.hdf".format(bdir,basen)):
				print("Dark/gain corrections were already performed")
			else:
				if options.verbose: print("Generating reference images for dark/gain correction(s)")
				if options.dark: dark = "--dark {}".format(options.dark)
				else: dark = ""
				if options.gain: gain = "--gain {}".format(options.gain)
				else: gain = ""
				if options.gaink2: gaink2 = "--gaink2 {}".format(options.gaink2)
				else: gaink2 = ""
				cmd = "e2ddd_movie.py {} {} {} {} --fixbadpixels --frames -v9".format(fname,dark,gain,gaink2)
				times = run(cmd,shell=True)
				if options.dark:
					dfn = options.dark.split(".")[0]
					os.unlink("{}_sig.hdf".format(dfn))
					os.unlink("{}_sum.hdf".format(dfn))
				if options.gain:
					gfn = options.gain.split(".")[0]
					os.unlink("{}_sig.hdf".format(gfn))
					os.unlink("{}_sum.hdf".format(gfn))
				shutil.move("{}_proc_corr.hdf".format(basen),bdir)
			fname = "{}/{}_proc_corr.hdf".format(bdir,basen)

		try: first,step,last = [int(i) for i in options.step.split(",")]
		except:
			print("Could not parse options.step. Please revise.")
			sys.exit(1)
		frames_hictrst = load_frames(fname,first,last,step)
		write_frames(frames_hictrst,hcname)
		middle = int(len(frames_hictrst)/2)

		hictrst_avg = average_frames(frames_hictrst)
		hictrst_avg.write_image("{}/hictrst_avg_noali.hdf".format(bdir))

		# Subtract gold

		if os.path.isfile(lcname):
			print("Found gold subtracted frames from previous analysis.")
		else:
			try: os.unlink("{}/hictrst_proc.hdf".format(bdir))
			except: pass
			if options.verbose: print("Erasing gold and other very high contrast features")
			cmd = "e2erasefiducials.py hictrst.hdf --goldsize={} --boxsize={} --parallel=thread:{} --lowpass --oversample={} --verbose={} --apix={}"
			cmd = cmd.format(options.goldsize,options.boxsize,options.threads,options.oversample,options.verbose,apix)
			times = run(cmd,cwd=bdir,shell=True)
			print(locallc)
			times = run("e2proc2d.py hictrst_proc.hdf {}".format(locallc),shell=True,cwd=bdir)

		frames_loctrst = load_frames(lcname)
		loctrst_avg = average_frames(frames_loctrst)
		loctrst_avg.write_image("{}/loctrst_avg_noali.hdf".format(bdir))

		cwd = os.getcwd()

		trans_lo = {}
		trans_hi = {}

		# PART 2: Run alignment programs
		#if os.path.isfile("{}/runtimes.txt".format(bdir)):
		#	os.unlink("{}/runtimes.txt".format(bdir))

		for pkg in sorted(pkgs.keys()):
			if pkg in options.include:
				prog = pkgs[pkg].split("/")[-1]

				if not options.skipalign: print("Running {} on {}".format(prog,fname))
				else: print("Parsing previous {} alignment results".format(prog,fname))

				pdir = "{}/{}".format(bdir,pkg)
				try: os.makedirs(pdir)
				except: pass

				if pkg == "DE":
					dev = "CPU"
					if not options.skipalign:
						src = "/".join(pkgs[pkg].split("/")[:-1]) + "/DE_process_frames.cfg"
						dst = cwd+"/DE_process_frames.cfg"
						try: os.symlink(src,dst)
						except: pass
						for fname in [hcname,lcname]:
							times=run("{} {} {} --run_cores {}".format(pkgs[pkg],pdir,fname,options.threads),shell=True,clear=True)
							write_runtime(bdir,pkg,dev,times,options.threads)
						try: os.unlink(dst)
						except: pass
					for f in os.listdir(pdir):
						if "translations_x" in f:
							if "hictrst" in f: hi_x = os.path.join(pdir,f)
							elif "loctrst" in f: lo_x = os.path.join(pdir,f)
						elif "translations_y" in f:
							if "hictrst" in f: hi_y = os.path.join(pdir,f)
							elif "loctrst" in f: lo_y = os.path.join(pdir,f)
					try:
						trans_hi[pkg] = parse_de(hi_x,hi_y)
						trans_lo[pkg] = parse_de(lo_x,lo_y)
					except: failed(pkg)

				if pkg == "EMAN2":
					dev = "CPU"
					if not options.skipalign:
						for fn in [hcname,lcname]: # high contrast, low contrast
							lfn = fn.split("/")[-1] # local file name
							if not os.path.isfile("{}/{}".format(pdir,lfn)):
								shutil.copy2(fn,"{}/{}".format(pdir,lfn))
							times=run("{} {} --align_frames --threads={} --allali --normalize".format(pkgs[pkg],lfn,options.threads),cwd=pdir,shell=True,clear=True)
							write_runtime(bdir,pkg,dev,times,options.threads)
					for f in os.listdir(pdir):
						if "hictrst_proc_info.txt" in f: hi = os.path.join(pdir,f)
						elif "loctrst_proc_info.txt" in f: lo = os.path.join(pdir,f)
					try:
						trans_hi[pkg] = parse_eman2(hi,middle)
						trans_lo[pkg] = parse_eman2(lo,middle)
					except: failed(pkg)

				if pkg == "IMOD":
					dev = "GPU"
					if not options.skipalign:
						for fn in [hcname,lcname]: # high contrast, low contrast
							lfn = fn.split("/")[-1] # local file name
							if not os.path.isfile("{}/{}".format(pdir,lfn)):
								shutil.copy2(fn,"{}/{}".format(pdir,lfn))
							alif = lfn.split(".")[0]+"_ali.mrc"
							times=run("{} -gpu 0 -xfext xf -input {} -output {}".format(pkgs[pkg],lfn,alif),cwd=pdir,shell=True,clear=True)
							write_runtime(bdir,pkg,dev,times,options.threads)
					for f in os.listdir(pdir):
						if "hictrst.xf" in f: hi = os.path.join(pdir,f)
						elif "loctrst.xf" in f: lo = os.path.join(pdir,f)
					try:
						trans_hi[pkg] = parse_imod(hi,middle)
						trans_lo[pkg] = parse_imod(lo,middle)
					except: failed(pkg)

				if pkg == "UCSF":
					dev = "GPU"
					if not options.skipalign:
						for fn in [hcname,lcname]: # high contrast, low contrast
							lfn = fn.split("/")[-1]
							lext = lfn.split(".")[-1]
							if lext == "hdf": # HDF format not supported by UCSF
								shutil.copy2(fn,"{}/{}".format(pdir,lfn))
								lfn_new = lfn.replace(".hdf",".mrcs")
								if not os.path.isfile("{}/{}".format(pdir,lfn_new)):
									times=run("e2proc2d.py {} {}".format(lfn,lfn_new),cwd=pdir,shell=True)
									os.remove("{}/{}".format(pdir,lfn))
								lfn = lfn_new
								lext = "mrcs"
							else:
								if not os.path.isfile("{}/{}".format(pdir,lfn)):
									shutil.copy2(fn,"{}/{}".format(pdir,lfn))
							times=run("{} {} -srs 1 -ssc 1 -atm 1".format(pkgs[pkg],lfn),cwd=pdir,shell=True,clear=True)
							write_runtime(bdir,pkg,dev,times,options.threads)
					for f in os.listdir(pdir):
						if "hictrst_Log.txt" in f: hi = os.path.join(pdir,f)
						elif "loctrst_Log.txt" in f: lo = os.path.join(pdir,f)
					try:
						trans_hi[pkg] = parse_ucsf(hi,middle)
						trans_lo[pkg] = parse_ucsf(lo,middle)
					except: failed(pkg)

				if pkg == "UNBLUR":
					dev = "CPU"
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

						for fn in [hcname,lcname]:
							lfn = fn.split("/")[-1]
							lext = lfn.split(".")[-1]
							if lext == "hdf": # HDF format not supported by unblur
								shutil.copy2(fn,"{}/{}".format(pdir,lfn))
								lfn_new = lfn.replace(".hdf",".mrcs")
								if not os.path.isfile("{}/{}".format(pdir,lfn_new)):
									times=run("e2proc2d.py {} {}".format(lfn,lfn_new),cwd=pdir,shell=True)
									os.remove("{}/{}".format(pdir,lfn))
								lfn = lfn_new
								lext = "mrcs"
							else:
								if not os.path.isfile("{}/{}".format(pdir,lfn)):
									shutil.copy2(fn,"{}/{}".format(pdir,lfn))
							if lext == "mrcs": # MRCS format not supported by unblur
								lfn_new = lfn.replace(".mrcs",".mrc")
								try:
									shutil.move("{}/{}".format(pdir,lfn),"{}/{}".format(pdir,lfn_new))
								except: pass # file likely already exists
								lfn = lfn_new
								lext = "mrc"
							nfs = len(frames_hictrst)
							alis =lfn.replace(".mrc","_ali_sum.mrc")
							shft =lfn.replace(".mrc","_shifts.txt")
							ali=lfn.replace(".mrc","_ali.mrc")
							frc=lfn.replace(".mrc","_frc.txt")
							cmd = template.format(prog=pkgs[pkg],fname=lfn, nframes=nfs, alignsum=alis, shiftfile=shft, apix=apix, dosefilt="NO", saveali="YES", aliname=ali, advopts="NO", frcname=frc, minsrch=2.0, maxsrch=200.0, bfact=1500, vfmwidth=1, hfmwidth=1, thresh=0.1, maxiter=10, verbose="NO")
							times=run(cmd,shell=True,cwd=pdir,clear=True)
							write_runtime(bdir,pkg,dev,times,options.threads)
					for f in os.listdir(pdir):
						if "hictrst_shifts.txt" in f: hi = os.path.join(pdir,f)
						elif "loctrst_shifts.txt" in f: lo = os.path.join(pdir,f)
					try:
						trans_hi[pkg] = parse_unblur(hi,middle)
						trans_lo[pkg] = parse_unblur(lo,middle)
					except: failed(pkg)

				if pkg == "LMBFGS":
					dev = "CPU"
					if not options.skipalign:
						template = """time {prog} << eot
{movielist}
{boxsize},{psize},{nsigma},{rmax1},{rmax2},{smooth}
{bfactor}
{framefirst},{framelast},{zeroframe}
{factr}
{inpath}
{outpath}
{shfext}
{vecext}
eot
"""

						for fn in [hcname,lcname]:
							lfn = fn.split("/")[-1]
							lext = lfn.split(".")[-1]
							if lext == "hdf": # HDF format not supported by alignframes_lmbfgs
								shutil.copy2(fn,"{}/{}".format(pdir,lfn))
								lfn_new = lfn.replace(".hdf",".mrcs")
								if not os.path.isfile("{}/{}".format(pdir,lfn_new)):
									times=run("e2proc2d.py {} {}".format(lfn,lfn_new),cwd=pdir,shell=True)
									try: os.unlink("{}/{}".format(pdir,lfn))
									except: pass
								lfn = lfn_new
								lext = "mrcs"
							else:
								if not os.path.isfile("{}/{}".format(pdir,lfn)):
									shutil.copy2(fn,"{}/{}".format(pdir,lfn))
							mlfn = "{}/movielist.txt".format(pdir)
							with open(mlfn,"w") as movielist:
								movielist.write("{}".format(lfn))
							bs = int(max(3600,min(frames_hictrst[0]["nx"]*0.75,frames_hictrst[0]["ny"]*0.75)))
							ps = apix
							nsig = 5
							rm1 = 500
							rm2 =100
							smooth = "0d0"
							bfact = 2000
							ffirst = 1
							flast = 0
							fmiddle = middle
							factr = "1d1"
							inpath = "./"
							outpath = "./"
							cmd = template.format(prog=pkgs[pkg],movielist="movielist.txt", boxsize=bs, psize=ps,
								nsigma=nsig, rmax1=rm1, rmax2=rm2, smooth=smooth, bfactor=bfact, framefirst=ffirst,
								framelast=flast, zeroframe=fmiddle, factr=factr, inpath=inpath, outpath=outpath,
								shfext="shf", vecext="vec")
							times=run(cmd,shell=True,exe="/bin/csh",cwd=pdir,clear=True)
							write_runtime(bdir,pkg,dev,times,options.threads)
					for f in os.listdir(pdir):
						if f == "hictrst.shf":
							hi = os.path.join(pdir,f)
						elif f == "loctrst.shf":
							lo = os.path.join(pdir,f)
					try:
						trans_hi[pkg] = parse_lmbfgs(hi,middle)
						trans_lo[pkg] = parse_lmbfgs(lo,middle)
					except: failed(pkg)

		print("")

		# STEP 3: Analyze results

		# Question 1: How quickly do frame alignment alorithms run
		print("RUNTIMES")
		with open("{}/runtimes.txt".format(bdir),"r") as f:
			for l in f: print(l.strip())

		print("")

		# Question 2: How does calculated motion differ between the most commonly used alignment algorithms?
		print("PREDICTED TRANSLATIONS")
		trans_orig, mags_orig = plot_traj(trans_hi,bdir,plot=options.plot)

		print("")

		# Question 3: How different are predicted frame translations with and without high contrast features?
		print("CONTRAST-BASED TRANSLATION COMPARISON")
		trans_hi, trans_lo, ssdfs = plot_differences(trans_hi,trans_lo,bdir,plot=options.plot)

		print("")

		# Question 4: How do alignment algorithms influence power spectral coherence?
		print("COHERENCE-BASED POWER SPECTRAL COMPARISON")
		# no algorithm seems to produce measurably better coherence (using calc_cips_scores).
		hi_ips,lo_ips,hi_cps,lo_cps = calc_coherence(frames_hictrst,frames_loctrst,trans_hi,trans_lo,plot=options.plot):

		print("")

		# Question 5: How do alignment algorithms influence real-space contrast?
		print("CONTRAST-BASED IMAGE COMPARISON")
		#measure_contrast(frames_hictrst)
		
#############################
# Parse Output Translations #
#############################

def parse_imod(fn,middle):
    with open(fn) as f:
        trans = [np.asarray(txt.split()[4:6]) for txt in f.readlines()]
    xf = np.asarray(trans).astype(float)
    return xf-xf[middle]

def parse_de(fnx,fny):
    with open(fnx) as xf:
        xs = np.asarray(xf.read().replace("\n","").split("\t"))[1:].astype(float)
    with open(fny) as yf:
        ys = np.asarray(yf.read().replace("\n","").split("\t"))[1:].astype(float)
    xf = np.vstack([xs,ys]).T
    return xf

def parse_ucsf(fn,middle):
    trans = []
    with open(fn) as f:
        for txt in f.readlines():
            if "......Shift of Frame #" in txt:
                trans.append(np.asarray(txt.replace("......Add Frame #","").split(":")[-1].split())[:2])
    xf = np.asarray(trans).astype(float)
    return xf-xf[middle]

def parse_unblur(fn,middle):
    lines = []
    with open(fn) as f:
        for txt in f.readlines():
            if "#" not in txt:
                lines.append(txt.split())
            elif "Pixel size (A):" in txt:
                apix = float(txt.split()[-1])
    trans = -1*np.asarray([[x,y] for (x,y) in zip(lines[0],lines[1])]).astype(float)
    xf = trans / apix
    return -xf+xf[middle]

def parse_lmbfgs(fn,middle):
	lines = []
	with open(fn) as f:
		for i,l in enumerate(f):
			if i != 0: lines.append(l.strip().split()) # skip header on line 0
	xf = -1*np.asarray(lines)[:,1:].astype(float) # skip frame numbers
	return xf - xf[middle]

def parse_eman2(fn,middle):
    with open(fn) as f:
        trans = [txt.split("\t")[1:3] for txt in f.readlines()[1:]]
    xf = np.asarray(trans).astype(float)
    return xf - xf[middle]

####################
# Compare and plot #
####################

def plot_traj(trans,bdir,nsig=1,plot=False):
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
	fig = plt.figure(figsize=(20,10))
	ax1 = fig.add_subplot(211)
	cctr = 0
	for key in trans.keys():
		if key not in exclude:
			if key != "MEAN":
				ax1.plot(trans[key][:,0],trans[key][:,1],colors[cctr],label=key,alpha=0.8)
				ax1.scatter(trans[key][:,0],trans[key][:,1],color=colors[cctr],alpha=0.8)
			else:
				ax1.plot(trans[key][:,0],trans[key][:,1],'k-',linewidth=1.5,label=key,alpha=1)
				ax1.scatter(trans[key][:,0],trans[key][:,1],color='k',alpha=1)
				ax1.errorbar(trans[key][:,0],trans[key][:,1],xerr=err[:,0],yerr=err[:,1],color='k',alpha=0.25,label=siglabel)
			cctr += 1
	ax1.legend(loc="best")#,bbox_to_anchor=(0.95, 0.45))
	ax1.set_title("Frame Trajectory")
	ax1.set_xlabel("X Frame Translation (pixels)")
	ax1.set_ylabel("Y Frame Translations (pixels)")
	ax2 = fig.add_subplot(212)
	cctr = 0
	print("PKG\tABS(DIST)")
	for key in mags.keys():
		if key != "MEAN":
			ax2.plot(xx,mags[key],colors[cctr],label=key,alpha=1)
		else:
			ax2.plot(xx,np.zeros_like(xx),"k--",label=key,alpha=1)
			ax2.errorbar(xx,np.zeros_like(xx),xerr=err[:,0],yerr=err[:,1],color='k',alpha=0.6,label=siglabel)
		cctr+=1
		print("{}\t{}".format(key,np.sum(np.abs(mags[key]))))
	ax2.set_xlabel("Frame Number (#)")
	ax2.set_ylabel("Relative shift magnitude (pixels)")
	ax2.set_title("Frame shift magnitude (relative to MEAN)")
	ax2.legend(loc="best")
	fig.tight_layout()
	plt.savefig("{}/trans.png".format(bdir))
	if plot: plt.show()
	return trans, mags

def plot_differences(trans_hi,trans_lo,bdir,nsig=1,plot=False):
	hts = []
	for k in sorted(trans_hi.keys()):
		if k not in ["ALL","MEAN","STD","VAR"]:
			hts.append(trans_hi[k])
	hts = np.hstack(hts)
	np.savetxt("{}/hi_trans.txt".format(bdir),hts)

	lts = []
	for k in sorted(trans_lo.keys()):
		if k not in ["ALL","MEAN","STD","VAR"]:
			lts.append(trans_lo[k])
	lts = np.hstack(lts)
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

	fig2 = plt.figure(figsize=(20,12))
	ax1 = fig2.add_subplot(211)
	cctr = 0
	print("PKG\tRMSD(HI,LO)")
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
					ax1.plot(trans_hi[key][:,0],trans_hi[key][:,1],"k--",linewidth=2,label="{} high".format(key))
					ax1.plot(trans_lo[key][:,0],trans_lo[key][:,1],"k.",linewidth=2,label="{} low".format(key))
				sse.write("{}\t{}\n".format(key,ssdfs[key]))
				cctr+=1
	ax1.set_xlabel("X-shift (pixels)")
	ax1.set_ylabel("Y-shift (pixels)")
	ax1.set_title("Frame trajectory")
	ax1.legend(loc="best")
	ax2 = fig2.add_subplot(212)
	cctr = 0
	for key in mags_hi.keys():
		if key != "MEAN":
			ax2.plot(xx,mags_hi[key],label="{} high".format(key),alpha=1)
		else:
			#ax2.plot(xx,np.zeros_like(xx),"k--",label=key,alpha=1)
			ax2.errorbar(xx,np.zeros_like(xx),yerr=err_hi[:,1],color='k',alpha=0.6,label="{} high".format(siglabel),linewidth=2)
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
	fig2.tight_layout()
	plt.savefig("{}/differences.png".format(bdir))
	if plot: plt.show()
	return trans_hi, trans_lo, ssdfs

#########################
# Compute Power Spectra #
#########################

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

def calc_coherence(frames_hi,frames_lo,all_trans_hi,all_trans_lo,plot=False):
	# use incoherent power spectrum as "control"
	hi_twod_ipws = calc_incoherent_pws(frames_hi)
	hi_oned_ipws = make_oned_pws(hi_twod_ipws)

	lo_twod_ipws = calc_incoherent_pws(frames_lo)
	lo_oned_ipws = make_oned_pws(lo_twod_ipws)

	hi_oned_cpws = {}
	lo_oned_cpws = {}
	
	for key in all_trans_hi.keys():
		shifted_hi = shift_frames(frames_hi,all_trans_hi[key])
		hi_twod_cpws = calc_coherent_pws(shifted_hi)
		hi_oned_cpws[key] = make_oned_pws(hi_twod_cpws)

		shifted_lo = shift_frames(frames_lo,all_trans_lo[key])
		lo_twod_cpws = calc_coherent_pws(shifted_lo)
		lo_oned_cpws[key] = make_oned_pws(lo_twod_cpws)

	if plot: plot_oned_spectra(hi_ipws,lo_ipws,hi_cpws,lo_cpws)
	return hi_oned_ipws,lo_oned_ipws,hi_oned_cpws,lo_oned_cpws

def make_oned_pws(pws):
	return pws.calc_radial_dist(pws["nx"]/2.,0,1.0,1)[1:]

def plot_oned_spectra(hi_ipws,lo_ipws,hi_cpws,lo_cpws):
	fig = plt.figure()

	ax1 = fig.add_subplot(211)
	ax1.plot(hi_ipws,label="Incoherent")
	for alg in hi_cpws.keys():
		ax1.plot(hi_cpws[alg],label=alg)
	ax1.legend()
	ax1.set_title("High Contrast Power Spectra")

	ax2 = fig.add_subplot(212)
	ax2.plot(lo_ipws,label="Incoherent")
	for alg in lo_cpws.keys():
		ax2.plot(lo_cpws[alg],label=alg)
	ax2.legend()
	ax2.set_title("Low Contrast Power Spectra")

	plt.show()

def calc_cips_scores(ftypes,bs=1024): #cips = coherent & incoherent power spectra
	# for comparison of coherent and incoherent power spectra...not such a good metric.
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

#####################
# Utility functions #
#####################

def load_frames(fn,first=0,last=-1,step=1):
	nimgs = EMUtil.get_image_count(fn)

	if first < 0: first = 0
	elif first > nimgs: first = 0

	if last == -1: last = nimgs
	elif last > nimgs: last = nimgs
	elif last < 0: last = nimgs

	if step < 1: step = 1
	elif step > nimgs: step = nimgs

	frames = []
	for i in range(first,last,step):
		frames.append(EMData(fn,i))
	return frames

def write_frames(frames,fn):
	for i,f in enumerate(frames):
		f.write_image(fn,i)

def average_frames(frames):
	avgr = Averagers.get("mean")
	avgr.add_image_list(frames)
	return avgr.finish()

def failed(pkg):
	print("ERROR: Could not find frame shifts for {} aligner.".format(pkg))
	sys.exit(1)

def shift_frames(frames,trans):
	shifted = []
	for frame,(x,y) in zip(frames,trans):
		f = frame.copy()
		f.translate(x,y,0)
		shifted.append(f)
	return shifted

def write_runtime(bdir,pkg,device,times,threads):
	walltime = times[0]
	usrtime = times[1]
	systime = times[2]
	with open("{}/runtimes.txt".format(bdir),"a") as rt:
		if device == "GPU":
			rt.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(pkg,device,walltime,usrtime,systime,"NA"))
		else: # CPU
			rt.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(pkg,device,walltime,usrtime,systime,threads))

def run(cmd,shell=False,cwd=None,exe="/bin/sh",clear=False):
	if options.verbose: print(cmd.replace("\n"," "))
	if cwd == None:
		cwd = os.getcwd()
	if shell == False:
		cmd = cmd.split()
	if clear:
		try:
			cc = subprocess.Popen("/home/jmbell/src/utils/clearcache") # "sync; echo 3 > /proc/sys/vm/drop_caches"
			cc.communicate()
			print("Memory cache cleared.")
		except:
			print("Memory cache not cleared. Do not trust runtime results.")
	p = psutil.Popen("/usr/bin/time -p "+cmd, shell=shell, cwd=cwd, executable=exe,stderr=subprocess.PIPE) # ,stdout=subprocess.PIPE, stdin=subprocess.PIPE,

	_,err=p.communicate() # write stdout to terminal. only keep timing info.

	wt,ut,st = err.strip().split("\n")[-3:] # anything before this is exit status 1 or an actual error.

	times = [float(wt.split(" ")[-1]),float(ut.split(" ")[-1]),float(st.split(" ")[-1])]

	if options.verbose: print(err.strip())

	return times


def fit_defocus(img):
	ds=1.0/(img["apix_x"]*img["nx"])
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

if __name__ == "__main__":
	main()
