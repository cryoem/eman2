#!/usr/bin/env python
# === Stage 2: e2tomogram_stage2.py using JSON metadata ===
import sys
import json
import numpy as np
from types import SimpleNamespace
from EMAN2 import *

from datetime import datetime
from concurrent.futures import ThreadPoolExecutor, as_completed
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")



#     e2tomogram_stage2.py FullFrames1495FN.json SatMorningJun7th


#### make tomogram by tiles
#### this is faster and has less artifacts. but takes a lot of memory (~4x the tomogram)

def make_tomogram_tile(imgs, tltpm, options, basenameForOutputs, errtlt=[], clipz=-1):
	time0=time.time()
	ntilts=len(imgs)
	img0 = EMData(imgs[0][0], imgs[0][1])
	img0.write_image('test000.hdf')# PRB
	imgFN = imgs[0][0]
	scale=img0["apix_x"]/options.apix_init
	imgsz=min(img0["nx"],img0["ny"])
	if imgsz<=1024*1.1:
		b=1
	elif imgsz<=2048*1.1:
		b=2
	else:
		b=4
		#print("tiling only support for 1k and 2k tomograms...")
		#return make_tomogram(imgs, tltpm, options, errtlt=errtlt)

	print("Making bin{:d} tomogram by tiling...".format(int(np.round(scale))))
	tpm=np.array(tltpm.copy())
	print(scale)
	print(type(tpm))
	print('len(tpm) = ' + str(len(tpm)))
	#print(tpm)
	tpm[:,:2]/=scale; #   temp PRB

	print('errtlt')
	print(type(errtlt))

	if len(errtlt)==0:
		errtlt=np.zeros(ntilts)
		nrange=list(range(ntilts))
	else:
		nrange=np.argsort(errtlt)[:int(ntilts*options.tltkeep)]
		et=errtlt[nrange]
		nrange=nrange[et<500]
		#print(et)

	print('nrange ',len(nrange), len(np.unique(np.array(nrange))),np.unique(np.array(nrange)) ); #PRB
	print("Using {} out of {} tilts..".format(len(nrange), ntilts))
	nx, ny=img0["nx"], img0["ny"]

	print('nx,ny',nx,ny)#PRB

	if options.autoclipxy:
		bds=[]
		for t in tpm:
			rot=Transform({"type":"2d","tx":t[0], "ty":t[1],"alpha":t[2]})
			p=np.array([rot.transform([nx/2, ny/2, 0]),rot.transform([nx/2, -ny/2, 0])])
			p=np.max(abs(p), axis=0)

			bds.append(p)

		bds=np.array(bds)
		bds=np.median(abs(bds), axis=0)*2 ## so we clip a rectangle area that covers half of the tilt images

		outx=good_size(bds[0])
		outy=good_size(bds[1])
		print("Final tomogram shape: {} x {}".format(outx, outy))
	else:
		outx=outy=1024*b

	#outxy=1024*b
	sz=good_size(max(int(clipz*0.8), 256*b)) #### this is the output 3D size
	step=sz//2 #### distance between each tile
	if options.extrapad:
		pad=good_size(sz*2) #### this is the padded size in fourier space
	else:
		pad=good_size(sz*1.4) #### this is the padded size in fourier space

	if clipz>0:
		outz=clipz
	else:
		outz=sz#good_boxsize(sz*1.2)

	print('b , size , step, pad')
	print(b , sz , step, pad)
	#### we make 2 tomograms with half a box shift and average them together to compensate for boundary artifacts.


	#options.moretile=True
	if options.moretile:
		full3d=EMData(outx, outy, outz)
		#mem=(outx*outy*outz*4+pad*pad*pad*options.threads*4)
		#print("This will take {}x{}x{}x4 + {}x{}x{}x{}x4 = {:.1f} GB of memory...".format(outx, outy, outz, pad, pad, pad,options.threads, mem/1024**3))
		wtcon=1
	else:
		full3d=[EMData(outx, outy, outz), EMData(outx, outy, outz)]
		#mem=(outx*outy*outz*2*4+pad*pad*pad*options.threads*4)
		#print("This will take {}x{}x{}x2x4 + {}x{}x{}x{}x4 = {:.1f} GB of memory...".format(outx, outy, outz, pad, pad, pad,options.threads, mem/1024**3))
		wtcon=2.5


	# jsd=queue.Queue(0)
	# jobs=[]
	nstepx=int(outx/step/2)
	nstepy=int(outy/step/2)

	print(f"nstepx = {nstepx}, nstepy = {nstepy}")
	tilecount_per_tilt = sum([
		len(range(-nstepy, nstepy+1)) if options.moretile
		else len(range(-nstepy + stepx % 2, nstepy+1, 2))
		for stepx in range(-nstepx, nstepx+1)
	])
	print(f"Estimated tiles per tilt: {tilecount_per_tilt}")
	print(f"Expected total tiles: {tilecount_per_tilt * ntilts}")

	if options.ctf!=None:
		ctf=EMAN2Ctf()
		ctf.from_dict({
			"defocus":1.0, "voltage":options.ctf["voltage"], "bfactor":0., "cs":options.ctf["cs"],"ampcont":0, "apix":img0["apix_x"]})
		dfs=[]
	else:
		ctf=None

	# define number of tilts that we will read at same time:   PRB
	nImgChunk = 32

	ntiles = (2*nstepy+1)*(2*nstepx+1);#look below if not moretile  PRB


	print('nImgChunk,ntiles,ntilts')
	print(nImgChunk,ntiles,ntilts)



	if 0:  return

	if options.ctf is not None:
		print("Doing Ctf correction. Average defocus {:.2f}".format(np.mean(dfs)))

	#############
	if options.normslice:
		options.addnoise = True
		stepx, stepy, threed = make_tile(jobs[len(jobs)//2])

		img0 = threed.numpy().copy()
		std0 = np.std(img0, axis=(1, 2))
		std0[std0 == 0] = 1
		std0 = 1.0 / std0
		std3 = np.zeros_like(img0)
		std3 += std0[:, None, None]
		maskz = from_numpy(std3)

	options.addnoise = False

	# Non-round falloff mask
	x, y = np.indices((sz, sz), dtype=float) / sz - 0.5
	f = wtcon + np.exp(-(x**2 + y**2) / 0.1) - np.exp(-((np.abs(x) - 0.5)**2 + (np.abs(y) - 0.5)**2) / 0.1)
	f3 = np.repeat(f[None, :, :], outz, axis=0)
	msk = from_numpy(f3).copy()

	jobs = []
	for stepx in range(-nstepx, nstepx + 1):
		if options.moretile:
			yrange = range(-nstepy, nstepy + 1)
		else:
			yrange = range(-nstepy + stepx % 2, nstepy + 1, 2)
			ntiles = len(yrange)*(2*nstepx+1)

		for stepy in yrange:
			jobs.append((imgFN, nrange, tpm, sz, pad, step, stepx, stepy, outz, ctf, options))

	# Main tile loop, single-threaded

	# for ijob, job in enumerate(jobs):

	total = len(jobs)
	max_workers = min(total, os.cpu_count() or 3)
	max_workers = 5
	print(f"Using {max_workers} threads for {total} tiles")

	# Submit all tiles to worker threads; each worker runs make_tile(job)
	with ThreadPoolExecutor(max_workers=max_workers) as ex:
		futures = [ex.submit(make_tile, job) for job in jobs]

		# tileData = m[i2].get_clip(Region(tx - pad // 2, ty - pad // 2, pad, pad), fill=0)
		completed = 0
		for fut in as_completed(futures):
			# Each future returns (stepx, stepy, threed) from make_tile(args)
			stepx, stepy, threed = fut.result()

			# Post-processing & insertion remain on the main thread (safest with EMAN2)
			threed.mult(msk)
			if options.normslice:
				threed.mult(maskz)

			if options.moretile:
				full3d.insert_scaled_sum(
					threed,
					(int(stepx * step + outx // 2),
					int(stepy * step + outy // 2),
					outz // 2)
				)
			else:
				full3d[stepx % 2].insert_clip(
					threed,
					(int(stepx * step + outx // 2 - sz // 2),
					int(stepy * step + outy // 2 - sz // 2),
					outz // 2 - threed["nz"] // 2)
				)

			del threed
			completed += 1
			timestamp = datetime.now().strftime("%H:%M:%S")
			print(f"[{timestamp}] Completed {completed}/{total} (tile {stepx},{stepy})")

	# Finalizing
	if not options.moretile:
		full3d = full3d[0] + full3d[1]

	a = full3d.numpy()
	std = np.std(a[[0, -1]])
	full3d.mult(1.0 / std)

	full3d["zshift"] = 0
	apix = img0["apix_x"]

	full3d["apix_x"] = full3d["apix_y"] = full3d["apix_z"] = apix

	print("Reconstruction done ({:.1f} s). Now writing tomogram to disk...".format(time.time() - time0))
	return full3d


def make_tile(args):
	# Unpack inputs (no queue used)
	tiltsFN, tiltsIndices, tpm, sz, pad, step, stepx, stepy, outz, ctf, options = args

	recon = Reconstructors.get("fourier", {
		"sym": 'c1',
		"size": [pad, pad, pad]
#		"mode": options.reconmode
	})
	recon.setup()

	nImg = EMUtil.get_image_count(tiltsFN)
	blockSize = getattr(options, "blocksize", 25)  # smaller block to cut peak RSS

	# Loop over frames in blocks; read each block once, then clip ROIs for this tile
	for i in range(0, nImg, blockSize):
		if i%400==0: print(i)

		# Skip I/O if none of this block's frames are requested for this tile
		block_range = range(i, min(i + blockSize, nImg))
		if not any((ii in tiltsIndices) for ii in block_range):
			continue

		m = EMData.read_images(tiltsFN, list(block_range))  # sequential I/O for the block

		for i2 in range(len(m)):
			ii = i + i2
			if ii not in tiltsIndices:
				continue

			t = tpm[ii]
			pos = [stepx * step, stepy * step, 0]
			pxf = get_xf_pos(t, pos)
			tx = int(m[i2]["nx"] // 2 + pxf[0])
			ty = int(m[i2]["ny"] // 2 + pxf[1])

			# ROI from the already-loaded frame
			tileData = m[i2].get_clip(Region(tx - pad // 2, ty - pad // 2, pad, pad), fill=0)

			# Thread-safe CTF: build fresh per-slice object (don't mutate shared 'ctf')
			if ctf is not None and getattr(options, "ctf", None) is not None:
				rot = Transform({"type": "xyz", "xtilt": float(t[4]), "ytilt": float(t[3])})
				p1 = rot.transform(pos)
				pz = p1[2] * m[i2]["apix_x"] / 10000.0

				ctf_local = EMAN2Ctf()
				ctf_local.from_dict({
					"defocus": options.ctf["defocus"][ii] - pz,
					"voltage": options.ctf["voltage"],
					"bfactor": 0.0,
					"cs": options.ctf["cs"],
					"ampcont": 0,
					"apix": m[i2]["apix_x"]
				})
				ctf_local.set_phase(options.ctf["phase"][ii] * np.pi / 180.0)
				tileData["ctf"] = ctf_local

				fft1 = tileData.do_fft()
				flipim = fft1.copy()
				ctf_local.compute_2d_complex(flipim, Ctf.CtfType.CTF_SIGN)
				fft1.mult(flipim)
				tileData = fft1.do_ift()
				del flipim, fft1, ctf_local  # free temps early

			tileData.process_inplace("filter.ramp")
			tileData.process_inplace("xform", {"alpha": -t[2]})

			if getattr(options, "addnoise", False):
				tileData.process_inplace("math.addnoise", {"noise": 100})
				tileData.process_inplace("normalize")

			xf = Transform({"type": "xyz", "ytilt": t[3], "xtilt": t[4]})
			dy = (pad // 2) - np.cos(t[3] * np.pi / 180.0) * pad / 2

			msk = EMData(pad, pad)
			msk.to_one()
			edge = (sz // 10)
			msk.process_inplace("mask.zeroedge2d", {
				"x0": dy + edge,
				"x1": dy + edge,
				"y0": edge,
				"y1": edge
			})
			msk.process_inplace("mask.addshells.gauss", {
				"val1": 0,
				"val2": edge
			})

			tileData.mult(msk)
			del msk

			mp = recon.preprocess_slice(tileData, xf)
			recon.insert_slice(mp, xf, 1)
			del tileData, mp  # free per-slice buffers promptly

		del m  # free the whole block before reading the next one

	threed = recon.finish(True)

	#if options.reconmode == "gauss_2":  #PRB 06/04/2025;   should be just fourier
	if 1:
		threed.process_inplace("math.gausskernelfix", {"gauss_width": 4.0})

	threed.clip_inplace(Region((pad - sz) // 2, (pad - sz) // 2, (pad - outz) // 2, sz, sz, outz))
	if 0:  # PRB  06/04/2025
		threed.process_inplace("filter.lowpass.gauss", {"cutoff_abs": options.filterto})

	# Return results directly instead of using jsd
	return stepx, stepy, threed

def make_tile_Oct7th(args):
    tiltsFN, tiltsIndices, tpm, sz, pad, step, stepx, stepy, outz, _ctf_ignored, options = args

    # Reconstructor at the (smaller) pad
    recon = Reconstructors.get("fourier", {"sym": "c1", "size": [pad, pad, pad]})
    recon.setup()

    # Read header once to get dims & apix (header-only read is cheap)
    hdr = EMData(tiltsFN, 0, True)  # True => header only
    nx, ny = int(hdr["nx"]), int(hdr["ny"])
    apix = float(hdr["apix_x"])

    blockSize = 8  # smaller block to keep RSS down

    # Work set for ROI reading
    indices = sorted(tiltsIndices)
    nImg = EMUtil.get_image_count(tiltsFN)

    for i0 in range(0, nImg, blockSize):
        i1 = min(i0 + blockSize, nImg)
        for i in range(i0, i1):
            if i not in tiltsIndices:
                continue

            t = tpm[i]
            pos = [stepx * step, stepy * step, 0]
            pxf = get_xf_pos(t, pos)
            tx = int(nx // 2 + pxf[0])
            ty = int(ny // 2 + pxf[1])

            # --- ROI READ: read just the tile (no full-frame in RAM)
            reg = Region(tx - pad // 2, ty - pad // 2, pad, pad)
            try:
                tileData = EMData(tiltsFN, i, False, reg)  # ROI from file
            except Exception:
                # Fallback if file driver can’t ROI-read this format
                full = EMData.read_images(tiltsFN, [i])[0]
                tileData = full.get_clip(reg, fill=0)
                del full

            # CTF: build a fresh local CTF (never mutate a shared object)
            if options.ctf is not None:
                rot = Transform({"type": "xyz", "xtilt": float(t[4]), "ytilt": float(t[3])})
                p1 = rot.transform(pos)
                pz = p1[2] * apix / 10000.0

                ctf_local = EMAN2Ctf()
                ctf_local.from_dict({
                    "defocus": options.ctf["defocus"][i] - pz,
                    "voltage": options.ctf["voltage"],
                    "bfactor": 0.,
                    "cs": options.ctf["cs"],
                    "ampcont": 0,
                    "apix": apix,
                })
                ctf_local.set_phase(options.ctf["phase"][i] * np.pi / 180.0)
                tileData["ctf"] = ctf_local

                # Apply CTF in Fourier domain with minimal temporaries
                fft1 = tileData.do_fft()
                filt = fft1.copy()                         # filter buffer
                ctf_local.compute_2d_complex(filt, Ctf.CtfType.CTF_SIGN)
                fft1.mult(filt)
                tileData = fft1.do_ift()
                del filt, fft1, ctf_local

            # Small, in-place steps only
            tileData.process_inplace("filter.ramp")
            tileData.process_inplace("xform", {"alpha": -t[2]})

            if getattr(options, "addnoise", False):
                tileData.process_inplace("math.addnoise", {"noise": 100})
                tileData.process_inplace("normalize")

            xf = Transform({"type": "xyz", "ytilt": t[3], "xtilt": t[4]})
            dy = (pad // 2) - np.cos(t[3] * np.pi / 180.0) * pad / 2

            # Build mask per-slice (cheap) and apply
            msk = EMData(pad, pad); msk.to_one()
            edge = (sz // 10)
            msk.process_inplace("mask.zeroedge2d", {
                "x0": dy + edge, "x1": dy + edge, "y0": edge, "y1": edge
            })
            msk.process_inplace("mask.addshells.gauss", {"val1": 0, "val2": edge})
            tileData.mult(msk)
            del msk

            mp = recon.preprocess_slice(tileData, xf)
            recon.insert_slice(mp, xf, 1)
            del tileData, mp

        # encourage early free of temps referenced by C-level
        import gc; gc.collect()

    threed = recon.finish(True)  # True frees recon’s internal buffers
    # Light, final touches
    threed.process_inplace("math.gausskernelfix", {"gauss_width": 4.0})
    threed.clip_inplace(Region((pad - sz)//2, (pad - sz)//2, (pad - outz)//2, sz, sz, outz))
    return stepx, stepy, threed




def make_tile_Orig(args):
	# Unpack inputs (no queue used)
	tiltsFN, tiltsIndices, tpm, sz, pad, step, stepx, stepy, outz, ctf, options = args

	recon = Reconstructors.get("fourier", {
		"sym": 'c1',
		"size": [pad, pad, pad]
#		"mode": options.reconmode
	})
	recon.setup()

	nImg = EMUtil.get_image_count(tiltsFN)
	blockSize =25;# Static, consistent with other chunking

	for i in range(0,nImg,blockSize):
		if i%200 ==0: print(i)
		m = EMData.read_images(tiltsFN,range(i, min(i+blockSize,nImg)))
		for i2  in range(len(m)):
			if i+i2 not in tiltsIndices: continue
			t = tpm[i+i2]
			pos=[stepx*step,stepy*step,0]
			pxf=get_xf_pos(t, pos)
			tx = int( m[i2]["nx"]//2 +pxf[0])
			ty = int( m[i2]["ny"]//2 +pxf[1])
			tileData = m[i2].get_clip(Region(tx-pad//2,ty-pad//2, pad, pad), fill=0)
			if ctf!=None:
				rot=Transform({"type":"xyz","xtilt":float(t[4]),"ytilt":float(t[3])})
				p1=rot.transform(pos)
				pz=p1[2]*m[i2]["apix_x"]/10000.
				ctf.defocus=options.ctf["defocus"][i]-pz
				ctf.set_phase(options.ctf["phase"][i]*np.pi/180.)
				dfs.append(ctf.defocus)
				tileData["ctf"]=ctf
				fft1 = tileData.do_fft()
				flipim = fft1.copy()
				ctf.compute_2d_complex(flipim, Ctf.CtfType.CTF_SIGN)
				fft1.mult(flipim)
				tileData = fft1.do_ift()

			tileData.process_inplace("filter.ramp")
			tileData.process_inplace("xform", {"alpha": -t[2]})

			if options.addnoise:
				tileData.process_inplace("math.addnoise", {"noise": 100})
				tileData.process_inplace("normalize")

			xf = Transform({"type": "xyz", "ytilt": t[3], "xtilt": t[4]})
			dy = (pad // 2) - np.cos(t[3] * np.pi / 180.) * pad / 2

			msk = EMData(pad, pad)
			msk.to_one()
			edge = (sz // 10)
			msk.process_inplace("mask.zeroedge2d", {
				"x0": dy + edge,
				"x1": dy + edge,
				"y0": edge,
				"y1": edge
			})
			msk.process_inplace("mask.addshells.gauss", {
				"val1": 0,
				"val2": edge
			})

			tileData.mult(msk)
			mp = recon.preprocess_slice(tileData, xf)
			recon.insert_slice(mp, xf, 1)

	threed = recon.finish(True)

	#if options.reconmode == "gauss_2":  #PRB 06/04/2025;   should be just fourier
	if 1:
		threed.process_inplace("math.gausskernelfix", {"gauss_width": 4.0})

	threed.clip_inplace(Region((pad - sz) // 2, (pad - sz) // 2, (pad - outz) // 2, sz, sz, outz))
	if 0:#PRB  06/04/2025
		threed.process_inplace("filter.lowpass.gauss", {"cutoff_abs": options.filterto})

	# Return results directly instead of using jsd
	return stepx, stepy, threed


#### reconstruct tomogram...
def make_tomogram(imgs, tltpm, options, outname=None, padr=1.2,  errtlt=[], clipz=-1, doclip=True):
	num=len(imgs)
	scale=imgs[0]["apix_x"]/options.apix_init
	print("Making bin{:d} tomogram...".format(int(np.round(scale))))
	tlt_params=tltpm.copy()
	tlt_params[:,:2]/=scale

	#### sort tilt by loss to exclude the worst ones
	if len(errtlt)==0:
		errtlt=np.zeros(num)
		nrange=list(range(num))
	else:
		for nid in range(num):
			ytlt=tlt_params[nid][3]
			if ytlt < options.tltrange[0] or ytlt > options.tltrange[1]:
				errtlt[nid] = np.inf # ensure this tilt is excluded if outside desired tilt range
		nrange=np.argsort(errtlt)[:int(num*options.tltkeep)]
		et=errtlt[nrange]
		nrange=nrange[et<500]

	nx=imgs[0]["nx"]
	ny=imgs[0]["ny"]
	outxy=good_size(min(nx, ny))

	pad=good_size(outxy*padr)
	#############
	#clipz=options.clipz
	zthick=good_size(max(clipz*1.2, pad//2))
	if options.verbose:
		print("\t Image size: {:d} x {:d}".format(nx, ny))
		print("\tPadded volume to: {:d} x {:d} x {:d}".format(pad, pad, zthick))

	#recon=Reconstructors.get("fourier", {"sym":'c1',"size":[pad,pad,zthick], "mode":options.reconmode})# PRB 06/04/2025
	recon=Reconstructors.get("fourier", {"sym":'c1',"size":[pad,pad,zthick], "mode":"std"})
	recon.setup()

	#### prepare jobs
	jobs=[]
	for nid in range(num):
		exclude= nid not in nrange
		tpm=tlt_params[nid]
		pxf=get_xf_pos(tlt_params[nid], [0,0,0])
		xform={"type":"xyz","ztilt":tpm[2],"ytilt":tpm[3], "xtilt":tpm[4], "tx":pxf[0], "ty":pxf[1]}
		jobs.append([nid,imgs[nid],  recon, pad, xform, exclude, options])

	#### starting threads
	thrds=[threading.Thread(target=reconstruct,args=(i)) for i in jobs]
	thrtolaunch=0
	tsleep=threading.active_count()
	while thrtolaunch<len(thrds) or threading.active_count()>tsleep:
		if thrtolaunch<len(thrds) :
			while (threading.active_count()==options.threads ) : time.sleep(.1)
			if options.verbose : print("Inserting slice {}/{}".format(thrtolaunch,len(thrds)))
			thrds[thrtolaunch].start()
			thrtolaunch+=1
		else: time.sleep(1)
	for t in thrds: t.join()

	threed=recon.finish(True)

	threed.process_inplace("math.gausskernelfix",{"gauss_width":4.0})
	threed.process_inplace("normalize")
	if 0: #PRB 06/04/2025
		threed.process_inplace("filter.lowpass.gauss",{"cutoff_abs":options.filterto})
	#print(threed["nx"], threed["ny"], threed["nz"])
	if doclip:
		if clipz<0:
			threed.clip_inplace(Region(((pad-outxy)//2), ((pad-outxy)//2), 0, outxy, outxy, zthick))
		else:
			threed.clip_inplace(Region(((pad-outxy)//2), ((pad-outxy)//2), (zthick-clipz)//2, outxy, outxy, clipz))
	threed["zshift"]=0

	apix=imgs[0]["apix_x"]
	threed["apix_x"]=threed["apix_y"]=threed["apix_z"]=apix

	#if options.normslice and outxy>1000:
		#threed.process_inplace("normalize.rows")
	if outname:
		if options.compressbits<0: threed.write_image(outname)
		else: threed.write_compressed(outname,0,options.compressbits,nooutliers=True)
		if options.verbose: print("Map written to {}.".format(outname))

	return threed

#### reconstruction function for the subprocesses
def reconstruct(nid, img, recon, pad, xform,  exclude, options):
	m=img.copy()
	#### the ramp filter and decay edge helps soften the edge artifacts
	m.process_inplace("filter.ramp")
	m.process_inplace("normalize")
	m.process_inplace("mask.decayedge2d", {"width":int(pad//40)})
	p2=m.get_clip(Region(m["nx"]//2-pad//2,m["ny"]//2-pad//2, pad, pad), fill=0)
	#### give up on the subpixel accuracy since it does not really matter here..
	p2.translate(-int(xform["tx"]), -int(xform["ty"]), 0)
	p2.rotate(-xform["ztilt"],0,0)
	xf=Transform({"type":"xyz","ytilt":xform["ytilt"],"xtilt":xform["xtilt"]})

	#### mask out the extra information on the edge of high tilt
	dy=p2["nx"]//2-np.cos(xform["ytilt"]*np.pi/180.)*m["nx"]//2
	msk=p2.copy()
	msk.to_one()
	edge=int(old_div(pad,20))
	msk.process_inplace("mask.zeroedge2d",{"x0":dy+edge, "x1":dy+edge, "y0":edge, "y1":edge})
	msk.process_inplace("mask.addshells.gauss",{"val1":0, "val2":edge})
	p2.mult(msk)

	if not exclude:
		p3=recon.preprocess_slice(p2, xf)
		recon.insert_slice(p3,xf,1)

#### get 2D position on a tilt given 3D location
def get_xf_pos(tpm, pk):
	xf0=Transform({"type":"xyz","xtilt":tpm[4],"ytilt":tpm[3],"ztilt":tpm[2],"tx":tpm[0], "ty":tpm[1]})
	p1=xf0.transform([pk[0], pk[1], pk[2]])

	return [p1[0], p1[1]]



def main():
    usage = "e2tomogram_stage2.py  <full_options.json> <baseNameForOutputs>"
    parser = EMArgumentParser(usage=usage)
    #parser.add_pos_argument("--json", type=str, help="Path to full options JSON file")
    #parser.add_argument("--ppid", type=int, help="PID of the parent process")

    (options,args) =parser.parse_args()
    logid = E2init(sys.argv)

    baseNameForOutputs = args[1]

    # Load full options
    with open(args[0]) as f:
        optdict = json.load(f)

    # Convert to SimpleNamespace for attribute-style access
    fullopts = SimpleNamespace(**optdict)

    # Load data

    #tlt_params_FN = fullopts.basename +"_tlt_params.npy"
    #tlt_params = np.load(tlt_params_FN, allow_pickle=True)
    tlt_params = getattr(fullopts,"tlt_params")

    #pks = np.load("pks.npy", allow_pickle=True).tolist() if os.path.exists("pks.npy") else None

	#imgout = EMData.read_images(fullopts.inputname)

    stackFN = fullopts.inputname
    print(stackFN)

    nimg    = EMUtil.get_image_count(stackFN);


    imgout =[]
    for jImgNow in range(nimg):
        imgout.append([stackFN,jImgNow])


    #if pks is not None:
    #    fullopts.pks = pks

    # Provide safe defaults for missing fields
    errtlt = getattr(fullopts, "ali_loss", [0.0] * len(imgout))
    errtlt = np.array(errtlt)
    clipz = getattr(fullopts, "clipz", -1)




    # Run reconstruction
    threed = make_tomogram_tile(imgout, tlt_params, fullopts,baseNameForOutputs, errtlt=errtlt, clipz=clipz)

    # Save output
    outname = getattr(fullopts, "basename", "reconstruction.hdf")
    #if getattr(fullopts, "compressbits", -1) >= 0:
    #    threed.write_compressed(outname, 0, fullopts.compressbits, nooutliers=True)
    #else:
    #    threed.write_image(outname)
    os.makedirs("tomograms", exist_ok=True)
    outfile = os.path.join("tomograms", baseNameForOutputs + ".hdf")
    threed.write_image(outfile)

    print(f"[Stage 2] Reconstruction complete. Output: {outfile}")
    E2end(logid)

if __name__ == '__main__':
    main()

#     e2tomogram_stage2.py FullFrames1495FN.json SatMorningJun7th
