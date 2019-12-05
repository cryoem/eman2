#!/usr/bin/env python

from EMAN2 import *
import numpy as np
import os
import sys
from scipy.optimize import minimize
from scipy.optimize import fmin
from joblib import Parallel, delayed 
from EMAN2_utils import *

def main():
	parser = EMArgumentParser()
	parser.add_argument("--tomogram", required=True, help="Path to tomogram.")
	parser.add_argument("--raw_tiltseries", required=True, help="Path to tiltseries.")
	parser.add_argument("--filterto",default=0.125,type=float,help="cutoff_abs for lowpass filter applied to tilts and reconstruction projections. Useful to exclude high-resolution information we can't currently correct.")
	parser.add_argument("--niters", default=3,type=int,help="Number of iterations to run.")
	parser.add_argument("--threads", default=8,type=int,help="Number of threads to run in parallel on a single computer.")
	parser.add_argument("--outpath", default="ampcor", help="Name of output directory. Numbers will be appended to ensure each instance is written to a unique folder. Default is ampcor_XX.")
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	(options, args) = parser.parse_args()

	niters = options.niters

	info_file = info_name(options.raw_tiltseries)
	if os.path.isfile(info_file):
		jd = js_open_dict(info_file)
		try: 
			tiltparams = np.asarray(jd['tlt_params'])
			jd.close()
		except:
			jd.close()
			print("Could not locate tlt_params key in {}. Try running this through e2tomogram.py before proceeding.".format(info_file))
			sys.exit(1)
	else:
		print("Could not locate metadata file ({}) for {}. Exiting.".format(info_file,options.raw_tiltseries))
		sys.exit(1)

	if options.verbose: print("Loading input data...")

	if options.verbose: print("Reading in tomogram: {}".format(options.raw_tiltseries))
	r_phase=EMData(options.tomogram)
	rnx = r_phase["nx"]
	rny = r_phase["ny"]
	rnz = r_phase["nz"]
	apix = r_phase["apix_x"]

	## Not sure if this is correct, but it seems like it might help
	## amplitude contrast should result in img["minimum"]>=0, 
	## but phase effects could make that more complicated
	# if r_phase["minimum"]<0: 
	# 	r_phase -= r_phase["minimum"]

	if options.verbose: print("Reading in raw tiltseries: {}".format(options.raw_tiltseries))
	binfac = 1
	for f in [2,4,8,16]:
		if "_bin{}".format(f) in options.tomogram:
			binfac = f
	if binfac == 1: 
		tilts = EMData.read_images(options.raw_tiltseries)
	else:
		# not sure if this is read/write limited, but the fft resample process takes quite a bit of time on large images.
		ntlt = EMUtil.get_image_count(options.raw_tiltseries)
		tilts = Parallel(n_jobs=options.threads, verbose=0, backend="threading")(delayed(get_binned_tilt) (i,binfac,options) for i in range(ntlt))

	outdir = make_path(options.outpath)

	# initial parameters for optimization. Likely not the best values.
	incid = []
	for t in tilts:
		inc = t+t["minimum"]
		incid.append(inc["maximum"])
	incid = np.max(incid) # incident intensity (approximated)
	alpha = 1.0
	beta = 0.1 # ...just something to start with
	x0 = [alpha/beta,incid] # initial guess for optimization

	if options.verbose: print("Generating initial projections of reconstructed tomogram")

	# Generate initial projections of reconstructed tomogram

	projs = Parallel(n_jobs=options.threads, verbose=0, backend="threading")(delayed(make_projection) (tp,r_phase,binfac) for tp in tiltparams)
	for i,p in enumerate(projs):
		p.process_inplace("filter.lowpass.gauss",{"cutoff_abs":options.filterto})
		p.write_image("{}/projs_00.hdf".format(outdir),i)

	r_phase = make_tomogram(projs,tiltparams,options,outname="{}/threed_00.hdf".format(outdir), padr=1.2, clipz=rnz)
	
	for itr in range(niters):

		thisi = str(itr).zfill(2)
		nexti = str(itr+1).zfill(2)
		
		if options.verbose: print("Generating amplitude contrast representation")

		# Get differences between reconstruction projections and original tilt series
		
		logdiffs = []
		for i in range(len(tiltparams)):
			m = tilts[i].copy()
			if itr == 0: m.write_image("{}/tiltseries_{}.hdf".format(outdir,thisi),i)

			p = EMData("{}/projs_{}.hdf".format(outdir,thisi),i)
			p.clip_inplace(Region((rnx-m["nx"])//2,(rny-m["ny"])//2,m["nx"],m["ny"]))

			dmsk = p.process("threshold.notzero")
			p.process_inplace("normalize.toimage",{"to":m,"fourieramp":False,"ignore_zero":True})

			pzero = p*(1-dmsk)
			d = p - m*dmsk + pzero["mean_nonzero"]

			d.process_inplace("math.log")
			d.process_inplace("threshold.belowtozero",{"minval":0.0})

			d.write_image("{}/logdiffs_{}.hdf".format(outdir,nexti),i)
			logdiffs.append(d)
		
		# Reconstruct amplitude contrast representation

		r_amp = make_tomogram(logdiffs,tiltparams,options,outname="{}/ramplitude_{}.hdf".format(outdir,nexti), padr=1.2, clipz=rnz)
		# cmd = "e2make3dpar.py --input {inp} --output {out} --keep {keep} --apix {apix} --pad {pad} --sym {sym} --mode {mode} --altedgemask --threads {ncores} --outsize {nx}".format(
		# 	inp="{}/logdiffs_{}.hdf".format(outdir,nexti), 
		# 	out="{}/ramplitude_{}.hdf".format(outdir,nexti), 
		# 	keep=1.0, apix=apix, pad=-1, sym="c1", 
		# 	mode="gauss_2", ncores=options.threads, nx=rnx)
		# run(cmd,options)
		# r_amp = EMData("{}/ramplitude_{}.hdf".format(outdir,nexti))
		#r_amp.clip_inplace(Region(0,0,(r_amp["nz"]-rnz)//2,rnx,rny,rnz))

		if options.verbose: print("Optimizing parameters relating phase and amplitude contrast representations")

		# Optimize parameters relating phase and amplitude contrast representations	

		x_opt = fmin(cost, x0, disp=True,args=(r_phase,r_amp))
		print(x_opt)
		
		# construct phi volume from phase and amplitude contrast representations

		new_phi = (r_phase + x_opt[0]*r_amp/np.log(x_opt[1]))/2.
		new_phi.write_image("{}/phivol_{}.hdf".format(outdir,nexti),0)
		
		if options.verbose: print("Generating corrected tiltseries")

		# generate corrected tiltseries

		projs = Parallel(n_jobs=options.threads, verbose=0, backend="threading")(delayed(make_projection) (tp,new_phi,binfac) for tp in tiltparams)
		for i,p in enumerate(projs):
			p.process_inplace("filter.lowpass.gauss",{"cutoff_abs":options.filterto})
			p.write_image("{}/projs_{}.hdf".format(outdir,nexti),i)
		
		if options.verbose: print("Reconstructing corrected tiltseries")

		r_phase = make_tomogram(projs,tiltparams,options,outname="{}/threed_{}.hdf".format(outdir,nexti), padr=1.2, clipz=rnz)
		# Reconstruct new phase contrast representation
		#cmd = "e2make3dpar.py --input {inp} --output {out} --keep {keep} --apix {apix} --pad {pad} --sym {sym} --altedgemask --mode {mode} --threads {ncores} --outsize {nx}".format(
		#	inp="{}/projs_{}.hdf".format(outdir,nexti),
		#	out="{}/threed_{}.hdf".format(outdir,nexti),
		#	keep=1.0,apix=apix,pad=-1,sym='c1',
		#	mode='gauss_2',ncores=options.threads,nx=rnx)
		#run(cmd,options)
		
		# define new phase contrast representation and tiltseries

		#r_phase = EMData("{}/threed_{}.hdf".format(outdir,nexti))

		# Should we use the raw tiltseries or projections from the last iteration's reconstruction?

		#tilts = EMData.read_images("{}/projs_{}.hdf".format(outdir,nexti)) 

	print("DONE")

# Optimize parameters by minimizing ||R_phase+(a/b)*R_amp/log(I0)||^2
def cost(x,r_phase,r_amp):
	if x[1] <= 0.0 : return np.inf
	energy = r_phase - x[0] * ( r_amp / np.log(x[1]) )
	return np.sqrt(energy["square_sum"])/(energy["nx"]*energy["ny"]*energy["nz"])

def get_binned_tilt(i,binfac,options):
	t = EMData(options.raw_tiltseries,i)
	t.process_inplace("math.fft.resample",{"n":binfac})
	t.process_inplace("filter.lowpass.gauss",{"cutoff_abs":options.filterto})
	return t

def make_projection(param,threed,binfac):
	xf=Transform({"type":"xyz","tx":param[0]/binfac,"ty":param[1]/binfac,"ytilt":param[3],"xtilt":param[4],"ztilt":param[2]})
	p = threed.project("standard",xf)
	return p

def run(cmd,options):
	if options.verbose: print("{}: {}".format(time.ctime(time.time()),cmd))
	ret=launch_childprocess(cmd)
	if ret != 0: print("Error running: {}".format(cmd))
	return

def make_tomogram(imgs, tltpm, options, outname=None, padr=1.2, clipz=-1):
	num=len(imgs)
	scale=1.#old_div(imgs[0]["apix_x"],options.apix)
	print("Making bin{:d} tomogram...".format(int(np.round(scale))))
	ttparams=tltpm.copy()
	ttparams[:,:2]/=scale

	errtlt=np.zeros(num)
	nrange=list(range(num))

	nx=imgs[0]["nx"]
	ny=imgs[0]["ny"]
	outxy=good_size(max(nx, ny))

	pad=good_size(outxy*padr)
	zthick=good_size(pad//2)
	if options.verbose:
		print("\t Image size: {:d} x {:d}".format(nx, ny))
		print("\tPadded volume to: {:d} x {:d} x {:d}".format(pad, pad, zthick))
	
	recon=Reconstructors.get("fourier", {"sym":'c1',"size":[pad,pad,zthick], "mode":"gauss_2"})
	recon.setup()
	
	#### prepare jobes
	jobs=[]
	for nid in range(num):
		exclude= nid not in nrange
		tpm=ttparams[nid]
		pxf=get_xf_pos(ttparams[nid], [0,0,0])
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
	threed.process_inplace("normalize")
	threed.process_inplace("filter.lowpass.gauss",{"cutoff_abs":options.filterto})
	
	if clipz<0:
		threed.clip_inplace(Region(((pad-outxy)//2), ((pad-outxy)//2), 0, outxy, outxy, zthick))
	else:
		threed.clip_inplace(Region(((pad-outxy)//2), ((pad-outxy)//2), (zthick-clipz)//2, outxy, outxy, clipz))
	threed["zshift"]=0

	apix=imgs[0]["apix_x"]
	threed["apix_x"]=threed["apix_y"]=threed["apix_z"]=apix

	if outname:
		threed.write_image(outname)
		if options.verbose: print("Map written to {}.".format(outname))

	return threed

#### reconstruction function for the subprocesses
def reconstruct(nid, img, recon, pad, xform,  exclude, options):
	m=img.copy()
	#### the ramp filter and decay edge helps soften the edge artifacts
	m.process_inplace("filter.ramp")
	m.process_inplace("normalize")
	m.process_inplace("mask.decayedge2d", {"width":int(pad//20)})
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
	### first project the point to xy plane
	xf0=Transform({"type":"xyz","xtilt":float(tpm[4]),"ytilt":float(tpm[3])})
	p0=[pk[0], pk[1], pk[2]]
	p1=xf0.transform(p0)#).astype(int)

	### undo the 2D alignment
	xf1=Transform({"type":"2d","tx":tpm[0], "ty":tpm[1],"alpha":tpm[2]})
	p2=xf1.transform([p1[0], p1[1]])

	return [p2[0], p2[1]]

if __name__ == "__main__":
	main()