#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import os
from EMAN2 import *

def main():

	trans = {}

	e2_fn = "eman2/platelet-wt-murine-03_012_-0.9_info.txt"
	trans["EMAN2"] = parse_eman2(e2_fn)

	de_fnx = "de/platelet-wt-murine-03_012_-0.9/platelet-wt-murine-03_012_-0.9.mrcs.translations_x.txt"
	de_fny = "de/platelet-wt-murine-03_012_-0.9/platelet-wt-murine-03_012_-0.9.mrcs.translations_y.txt"
	trans["DE"] = parse_de(de_fnx,de_fny)

	imod_fn = "imod/platelet-wt-murine-03_012_-0.9.xf"
	trans["IMOD"] = parse_imod(imod_fn)

	ucsf_fn = "ucsf/platelet-wt-murine-03_012_-0.9_Log.txt"
	trans["UCSF"] = parse_ucsf(ucsf_fn)

	unblur_fn = "unblur/platelet-wt-murine-03_012_-0.9_shifts.txt"
	trans["UNBLUR"] = parse_unblur(unblur_fn)

	for key in trans.keys():
		print(key,trans[key].shape)

	plot_trans(trans)

def plot_trans(trans,nsig=1):

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

	fig = plt.figure(figsize=(16,16))

	ax1 = fig.add_subplot(211)

	for key in trans.keys():
		if key not in exclude:
			if key == "MEAN":
				ax1.plot(trans[key][:,0],trans[key][:,1],'k-',linewidth=1.5,label=key,alpha=1)
				ax1.errorbar(trans[key][:,0],trans[key][:,1],xerr=err[:,0],yerr=err[:,1],color='k',alpha=0.6,label=siglabel)
			else:
				ax1.plot(trans[key][:,0],trans[key][:,1],label=key,alpha=0.8)
	ax1.legend(loc="upper left",bbox_to_anchor=(0.95, 0.45))
	ax1.set_title("Measured Drift")
	ax1.set_xlabel("X Frame Translation (pixels)")
	ax1.set_ylabel("Y Frame Translations (pixels)")

	ax2 = fig.add_subplot(212)

	for key in mags.keys():
		if key != "MEAN":
			ax2.plot(xx,mags[key],label=key,alpha=1)
		else:
			ax2.plot(xx,np.zeros_like(xx),"k--",label=key,alpha=1)
			ax2.errorbar(xx,np.zeros_like(xx),xerr=err[:,0],yerr=err[:,1],color='k',alpha=0.6,label=siglabel)
	ax2.set_xlabel("Frame Number (#)")
	ax2.set_ylabel("Relative shift magnitude (pixels)")
	ax2.set_title("Frame shift magnitude (relative to MEAN)")
	ax2.legend(loc="upper left",bbox_to_anchor=(0.95, 0.95))

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
	return np.cumsum(trans - np.mean(trans,axis=0),axis=1)

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
