#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import os
from EMAN2 import *

def main():

	e2_fn = "eman2/platelet-wt-murine-03_012_-0.9_info.txt"
	e2_trans = parse_eman2(e2_fn)

	de_fnx = "de/platelet-wt-murine-03_012_-0.9/platelet-wt-murine-03_012_-0.9.mrcs.translations_x.txt"
	de_fny = "de/platelet-wt-murine-03_012_-0.9/platelet-wt-murine-03_012_-0.9.mrcs.translations_y.txt"
	de_trans = parse_de(de_fnx,de_fny)

	imod_fn = "imod/platelet-wt-murine-03_012_-0.9.xf"
	imod_trans = parse_imod(imod_fn)

	ucsf_fn = "ucsf/platelet-wt-murine-03_012_-0.9_Log.txt"
	ucsf_trans = parse_ucsf(ucsf_fn)

	unblur_fn = "unblur/platelet-wt-murine-03_012_-0.9_shifts.txt"
	unblur_trans = parse_unblur(unblur_fn)

	all_trans = np.dstack([e2_trans,de_trans,imod_trans,ucsf_trans,unblur_trans])
	all_mean = np.mean(all_trans,axis=2)
	all_var = np.var(all_trans,axis=2)
	all_std = np.std(all_trans,axis=2)

	x = all_mean[:,0][::skip]
	y = all_mean[:,1][::skip]
	nsig = 1
	err = nsig * all_std

	fig = plt.figure(figsize=(16,16))
	ax = fig.add_subplot(221)
	ax.plot(e2_trans[:,0],e2_trans[:,1],label="EMAN",alpha=0.9) # EMAN2
	ax.plot(de_trans[:,0],de_trans[:,1],label="DE Script",alpha=0.9) # DE
	ax.plot(imod_trans[:,0],imod_trans[:,1],label="IMOD",alpha=0.9) # IMOD
	ax.plot(ucsf_trans[:,0],ucsf_trans[:,1],label="UCSF",alpha=0.9) # UCSF
	ax.plot(unblur_trans[:,0],unblur_trans[:,1],label="UNBLUR",alpha=0.9) # UNBLUR
	ax.plot(x,y,'k-',linewidth=1.5,label="$\mu$",alpha=1)
	if nsig == 1: ax.errorbar(x,y,xerr=err[:,0],yerr=err[:,1],color='k',alpha=0.6,label="$\sigma$")
	else: ax.errorbar(x,y,xerr=err[:,0],yerr=err[:,1],color='k',alpha=0.6,label="{} $\sigma$".format(nsig))
	ax.legend(loc="best")#loc="upper left",bbox_to_anchor=(1.1, 1.05))
	ax.set_title("Measured Drift")
	ax.set_xlabel("X Frame Translation (pixels)")
	ax.set_ylabel("Y Frame Translations (pixels)")

	xx = np.arange(1,len(x)+1,1)

	mean_mag = np.sqrt(x**2+y**2)
	e2_mag = np.sqrt(e2_trans[:,0]**2+e2_trans[:,1]**2) - mean_mag
	de_mag = np.sqrt(de_trans[:,0]**2+de_trans[:,1]**2) - mean_mag
	imod_mag = np.sqrt(imod_trans[:,0]**2+imod_trans[:,1]**2) - mean_mag
	ucsf_mag = np.sqrt(ucsf_trans[:,0]**2+ucsf_trans[:,1]**2) - mean_mag
	unblur_mag = np.sqrt(unblur_trans[:,0]**2+unblur_trans[:,1]**2) - mean_mag

	#fig = plt.figure(figsize=(16,8))
	ax = fig.add_subplot(222)
	ax.plot(xx,e2_mag,color="b",label="EMAN",alpha=1) # EMAN2
	ax.plot(xx,de_mag,color="g",label="DE Script",alpha=1) # DE
	ax.plot(xx,imod_mag,color="r",label="IMOD",alpha=1) # IMOD
	ax.plot(xx,ucsf_mag,color="c",label="UCSF",alpha=1) # UCSF
	ax.plot(xx,unblur_mag,color="m",label="UNBLUR",alpha=1) # UNBLUR
	ax.plot(xx,np.zeros_like(xx),"k--",label="$\mu$",alpha=1)
	ax.set_xlabel("Frame Number")
	ax.set_ylabel("Shift Magnitude (relative mean shift magnitude)")
	ax.set_title("Frame Translations relative to $\mu$")
	ax.legend(loc="upper center")
	plt.show()



pkgs = {"EMAN2":"movie_ccf.py",
		"UCSF":"dosefgpu_driftcorr",
		"UNBLUR":"unblur_openmp_7_17_15.exe",
		"DE":"DE_process_frames-2.8.1.py",
		"IMOD":"alignframes"}

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
	return np.cumsum(trans - np.mean(trans,axis=0),axis=1)

def parse_imod(fn):
	with open(fn) as f:
		trans = np.asarray([np.asarray(txt.split()[4:6]).astype(float) for txt in f.readlines()])
	imod_trans = np.cumsum(trans - np.mean(trans,axis=0),axis=1)

def parse_ucsf(fn):
	fn = "ucsf/platelet-wt-murine-03_012_-0.9_Log.txt"
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

if __name__ == "__main__":
	main()
