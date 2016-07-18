#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import os
from EMAN2 import *

colors = ['b', 'g', 'r', 'c', 'm', 'y']

def main():
	show = True

	# PART 1: How does calculated motion differ between the most commonly used alignment algorithms?

	trans_wmg = {}

	de_fnx = "de/platelet-wt-murine-03_012_-0.9/platelet-wt-murine-03_012_-0.9.mrcs.translations_x.txt"
	de_fny = "de/platelet-wt-murine-03_012_-0.9/platelet-wt-murine-03_012_-0.9.mrcs.translations_y.txt"
	trans_wmg["DE"] = parse_de(de_fnx,de_fny)

	e2_fn = "eman2/platelet-wt-murine-03_012_-0.9_info.txt"
	trans_wmg["EMAN2"] = parse_eman2(e2_fn)

	imod_fn = "imod/platelet-wt-murine-03_012_-0.9.xf"
	trans_wmg["IMOD"] = parse_imod(imod_fn)

	ucsf_fn = "ucsf/platelet-wt-murine-03_012_-0.9_Log.txt"
	trans_wmg["UCSF"] = parse_ucsf(ucsf_fn)

	unblur_fn = "unblur/platelet-wt-murine-03_012_-0.9_shifts.txt"
	trans_wmg["UNBLUR"] = parse_unblur(unblur_fn)

	trans_wmg, mags_wmg = plot_trans(trans_wmg,show=show)

	# PART 2: How different are predicted frame translations with and without high contrast features?

	trans_lo = {} # low contrast
	trans_hi = {} # high contrast

	de_fnx = "de/cell_region/cell_region.mrcs.translations_x.txt"
	de_fny = "de/cell_region/cell_region.mrcs.translations_y.txt"
	trans_lo["DE"] = parse_de(de_fnx,de_fny)
	de_fnx = "de/feducial_region/feducial_region.mrcs.translations_x.txt"
	de_fny = "de/feducial_region/feducial_region.mrcs.translations_y.txt"
	trans_hi["DE"] = parse_de(de_fnx,de_fny)

	e2_fn = "eman2/cell_region_info.txt"
	trans_lo["EMAN2"] = parse_eman2(e2_fn)
	e2_fn = "eman2/feducial_region_info.txt"
	trans_hi["EMAN2"] = parse_eman2(e2_fn)

	imod_fn = "imod/cell_region.xf"
	trans_lo["IMOD"] = parse_imod(imod_fn)
	imod_fn = "imod/feducial_region.xf"
	trans_hi["IMOD"] = parse_imod(imod_fn)

	ucsf_fn = "ucsf/cell_region_Log.txt"
	trans_lo["UCSF"] = parse_ucsf(ucsf_fn)
	ucsf_fn = "ucsf/feducial_region_Log.txt"
	trans_hi["UCSF"] = parse_ucsf(ucsf_fn)

	unblur_fn = "unblur/cell_region_shifts.txt"
	trans_lo["UNBLUR"] = parse_unblur(unblur_fn)
	unblur_fn = "unblur/feducial_region_shifts.txt"
	trans_hi["UNBLUR"] = parse_unblur(unblur_fn)

	ssdfs = plot_differences(trans_hi,trans_lo,show=show)

	# Part 3?

	#frames = load_frames("platelet-wt-murine-03_012_-0.9.mrcs")
	frames_hi = load_frames("feducial_region.mrcs")
	frames_lo = load_frames("cell_region.mrcs")

	ftypes = {}
	ftypes["hi_wmg"] = shift_by(frames_hi,trans_wmg)
	ftypes["lo_wmg"] = shift_by(frames_lo,trans_wmg)
	ftypes["lo_lo"] = shift_by(frames_lo,trans_lo)
	ftypes["lo_hi"] = shift_by(frames_lo,trans_hi)
	ftypes["hi_lo"] = shift_by(frames_hi,trans_lo)
	ftypes["hi_hi"] = shift_by(frames_hi,trans_hi)

	scores = calc_cips_scores(ftypes)
	plot_scores(scores,show=show)

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

def calc_cips_scores(ftypes,bs=256):
	scores = {}
	for pkg in pkgs.keys():
		scores[pkg] = []
		for ftype in sorted(ftypes.keys()):
			cps = calc_coherent_pws(ftypes[ftype][pkg],bs)
			ips = calc_incoherent_pws(ftypes[ftype][pkg],bs)
			frc = -1 * EMData.cmp(cps,'frc',ips)
			dot = -1 * EMData.cmp(cps,'dot',ips)
			osb = -1 * EMData.cmp(cps,'optsub',ips)
			sse = -1 * EMData.cmp(cps,'sqeuclidean',ips)
			scores[pkg].append([frc,dot,osb,sse])
		scores[pkg] = np.asarray(scores[pkg])
	return scores

def plot_scores(scores,show=False): # 1 axis per metric
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
	if show: plt.show()
	plt.savefig("scores.jpg")

def plot_trans(trans,nsig=1,show=False):
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
	ax2.legend(loc="best")#,bbox_to_anchor=(0.95, 0.95))

	fig.tight_layout()

	if show: plt.show()
	plt.savefig("trans.jpg")
	return trans, mags

def plot_differences(trans_hi,trans_lo,nsig=1,show=False):
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
	plt.savefig("differences.jpg")
	if show: plt.show()
	return trans_hi, trans_lo, ssdfs

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

def calc_incoherent_pws(frames,bs=256):
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

def calc_coherent_pws(frames,bs=256):
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

pkgs = {"EMAN2":"movie_ccf.py",
		"UCSF":"dosefgpu_driftcorr",
		"UNBLUR":"unblur_openmp_7_17_15.exe",
		"DE":"DE_process_frames-2.8.1.py",
		"IMOD":"alignframes"}

if __name__ == "__main__":
	main()
