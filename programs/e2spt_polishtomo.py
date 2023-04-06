#!/usr/bin/env python
# Muyuan Chen 2023-03
from EMAN2 import *
import numpy as np

emdir=e2getinstalldir()
sys.path.insert(0,os.path.join(emdir,'bin'))
from e2tomogram import *

def main():
	
	usage="""
	Polish a tomogram given the subtilt refinement of particles inside. 
	e2spt_polishtomo.py --fname tomograms/xxxx.hdf --path spt_xx 
	"""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--path", type=str,help="path of spt refinement", default=None)
	parser.add_argument("--fname", type=str,help="name of tomogram", default=None)
	parser.add_argument("--nneighbor", type=int,help="number of neighbors", default=5)
	parser.add_argument("--res", type=float,help="lowpass filter the output to the target resolution.", default=50)
	parser.add_argument("--makeraw", action="store_true", default=False ,help="skip polish for testing")

	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)
		
	path=options.path

	print("Gathering metadata...")
	info3d=load_lst_params(f"{path}/particle_info_3d.lst")
	info2d=load_lst_params(f"{path}/particle_info_2d.lst")
	for i in range(99,0,-1):
		fm=f"{path}/aliptcls3d_{i:02d}.lst"
		if os.path.isfile(fm):
			options.loadali3d=fm
			break
	print("using 3d alignment from {}".format(options.loadali3d))

	for i in range(99,0,-1):
		fm=f"{path}/aliptcls2d_{i:02d}.lst"
		if os.path.isfile(fm):
			lst=load_lst_params(fm, range(10))
			options.loadali2d=fm
			break
	print("using 2d alignment from {}".format(options.loadali2d))

	alipm=load_lst_params(options.loadali2d)
	for i,a in zip(info2d, alipm):
		i["pastxf"]=a["xform.projection"]
		i["score"]=a["score"]
		if "defocus" in a:
			i["defocus"]=a["defocus"]
		else:
			i["defocus"]=0

	alipm=load_lst_params(options.loadali3d)
	for i,a in zip(info3d, alipm):
		i["xform.align3d"]=a["xform.align3d"]
		i["score"]=a["score"]

	#info2d=read_aliptcls(aliptcls, info2d)
	filenames, ptclcount=np.unique([d["src"] for d in info3d], return_counts=True)

	print("load {} particles from {} tomograms".format(len(info3d), len(filenames)))
	
	
	fname=[f for f in filenames if base_name(f)==base_name(options.fname)]
	if len(fname)==0:
		print("cannot find file")
		print(options.fname)
		print(filenames)
		
	fname=fname[0]	
	print(f"Polishing tomogram {fname}")
	d3d=[d for d in info3d if d["src"]==fname]
	f=fname.replace("particles3d", "particles")
	tid=np.array([d["tilt_id"] for d in info2d if d["src"]==f])
	tid=np.sort(np.unique(tid))
	print("Loading {} particles from {} tilts...".format(len(d3d), len(tid)))

	sel_coord=[]
	sel_score=[]
	sel_defocus=[]
	sel_dxy=[]
	tltang=[]
	for td in tid:

		d3ds=[d for d in info3d if d["src"]==fname]
		d2d=[]
		d3d=[]
		for d3 in d3ds:
			d2=[info2d[d] for d in d3["idx2d"]]
			d2=[d for d in d2 if d["tilt_id"]==td]
			if len(d2)==0: continue
			d2d.append(d2[0])
			d3d.append(d3)

		xfali=[d["pastxf"] for d in d2d]

		coord=np.array([d["coord"] for d in d3d])
		txfs=[d["xform.align3d"].inverse() for d in d3d]
		coord-=np.array([t.get_trans() for t in txfs])

		xfpj=[d["xform.projection"] for d in d2d]
		tltang.append(np.mean([x.get_params("xyz")["ytilt"] for x in xfpj]))
		xfraw=[a*b for a,b in zip(xfpj, txfs)]


		pastxf=([b*a.inverse()for a,b in zip(xfraw, xfali)])
		dxy=np.array([a.get_trans() for a in pastxf])
		score=[d["score"] for d in d2d]
		defocus=[d["defocus"] for d in d2d]

		sel_score.append(score)
		sel_dxy.append(dxy)
		sel_coord.append(coord)
		sel_defocus.append(defocus)

	plt_scr=[np.mean(s) for s in sel_score]
	plt_def=[np.mean(abs(np.array(s))) for s in sel_defocus]
	plt_dxy=[np.mean(np.linalg.norm(d, axis=1)) for d in sel_dxy]
	sel_tid=int(np.mean(tid))
	tltang=np.array(tltang)
		
	info=js_open_dict(info_name(fname))
	tfile=info["tlt_file"]
	if "defocus" in info:
		print("Loading CTF information. will do phase flipping for tomograms")
		options.ctf={	"defocus":info["defocus"], "phase":info["phase"], 
				"cs":info["cs"], "voltage":info["voltage"]}
	else:
		options.ctf=None
		
	imgs_raw=EMData.read_images(tfile)

	imgs_1k=[]
	for m in imgs_raw:
		e=m.process("math.fft.resample",{"n":4})
		e.process_inplace("normalize.edgemean")
		imgs_1k.append(e)
		
	tltpm=np.array(info["tlt_params"]).copy()
	
	options.threads=12
	imgs=imgs_1k
	e=EMData(options.loadali3d,0, True)
	options.xfscale=imgs[0]["apix_x"]/e["apix_x"]
	print(f"apix of tomogram: {imgs[0]['apix_x']:.2f}, apix of particles: {e['apix_x']:.2f}, scale xform by {options.xfscale:.2f}")
	options.filterto=imgs[0]["apix_x"]/options.res

	num=len(imgs)
	scale=4
	imgsz=min(imgs[0]["nx"],imgs[0]["ny"])

	print("Making bin{:d} tomogram by tiling...".format(int(np.round(scale))))
	tpm=tltpm.copy()
	tpm[:,:2]/=scale

	nrange=list(range(num))
	nx, ny=imgs[0]["nx"], imgs[0]["ny"]
	
	
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
	print("tomogram shape: {} x {}".format(outx, outy))

	clipz=256
	sz=256 #### this is the output 3D size 
	pad=good_size(sz*1.4) #### this is the padded size in fourier space
	options.outz=outz=clipz
	options.step=step=sz//2

	full3d=EMData(outx, outy, outz)
	mem=(outx*outy*outz*4+pad*pad*pad*options.threads*4)
	print("This will take {}x{}x{}x4 + {}x{}x{}x{}x4 = {:.1f} GB of memory...".format(outx, outy, outz, pad, pad, pad,options.threads, mem/1024**3))
	wtcon=1
	jsd=queue.Queue(0)
	
	jobs=[]
	nstepx=int(outx/step/2)
	nstepy=int(outy/step/2)
	coord=sel_coord[0]/options.xfscale

	tpos=[]
	
	if options.ctf!=None:
		ctf=EMAN2Ctf()
		ctf.from_dict({
			"defocus":1.0, "voltage":options.ctf["voltage"], "bfactor":0., "cs":options.ctf["cs"],"ampcont":0, "apix":imgs[0]["apix_x"]})
		dfs=[]
		
	for stepx in range(-nstepx,nstepx+1):
		#### shift y by half a tile
		yrange=range(-nstepy,nstepy+1)
		
		for stepy in yrange:
			tiles=[]
			for i in range(num):
				if i in nrange:
					t=tpm[i]
					pos=[stepx*step,stepy*step,0]
					pxf=get_xf_pos(t, pos)
					img=imgs[i]
					m=img.get_clip(Region(img["nx"]//2-pad//2+pxf[0],img["ny"]//2-pad//2+pxf[1], pad, pad), fill=0)
					
					if options.ctf!=None:
						rot=Transform({"type":"xyz","xtilt":float(t[4]),"ytilt":float(t[3])})
						p1=rot.transform(pos)
						pz=p1[2]*img["apix_x"]/10000.
						ctf.defocus=options.ctf["defocus"][i]-pz
						ctf.set_phase(options.ctf["phase"][i]*np.pi/180.)
						dfs.append(ctf.defocus)
						m["ctf"]=ctf
					
					tiles.append(m)

			jobs.append((jsd, tiles, tpm, sz, pad, stepx, stepy, coord, sel_dxy, options))

	x,y=np.indices((sz,sz),dtype=float)/sz-.5
	#f=.25-(x**2+y**2)/2 + ((abs(x)-0.5)**2+(abs(y)-0.5)**2)/2
	f=wtcon+np.exp(-(x**2+y**2)/0.1) - np.exp(-((abs(x)-0.5)**2+(abs(y)-0.5)**2)/0.1)
	f3=np.repeat(f[None, :,:], outz, axis=0)
	msk=from_numpy(f3).copy()
	
	thrds=[threading.Thread(target=reconstruct_tile,args=([i])) for i in jobs]
	print("now start threads...")
	thrtolaunch=0
	tsleep=threading.active_count()
	
	while thrtolaunch<len(thrds) or threading.active_count()>tsleep or not jsd.empty():
		if thrtolaunch<len(thrds) :
			while (threading.active_count()==options.threads ) : time.sleep(.1)
			thrds[thrtolaunch].start()
			thrtolaunch+=1
		else: time.sleep(.1)
		
		if not jsd.empty():
			stepx, stepy, threed=jsd.get()
			threed.mult(msk)
			#### insert the cubes to corresponding tomograms
			full3d.insert_scaled_sum(
					threed,(int(stepx*step+outx//2),int(stepy*step+outy//2), outz//2))
			
	for t in thrds: t.join()
	
	apix=imgs[0]["apix_x"]
	full3d["apix_x"]=full3d["apix_y"]=full3d["apix_z"]=apix
	
	full3d.process_inplace("normalize")
	if options.makeraw:
		outname="tomograms/{}__bin4_raw.hdf".format(base_name(fname))
	else:
		outname="tomograms/{}__bin4_polish.hdf".format(base_name(fname))
	full3d.write_compressed(outname,0,8,nooutliers=True)
	
	print(f"Output written to {outname}")
	E2end(logid)

def reconstruct_tile(job):
	jsd, tiles, tpm, sz, pad, stepx, stepy, coord, sel_dxy, options=job
	recon=Reconstructors.get("fourier", {"sym":'c1',"size":[pad,pad,pad], "mode":"trilinear"})
	recon.setup()
	step=options.step
	outz=options.outz
	pos=np.array([(stepx*step),int(stepy*step), 0])
	d=np.linalg.norm(coord-pos, axis=1)
	nnb=options.nneighbor
	di=np.argsort(d)[:nnb]
	d=d[di]

	wt=np.exp(-d*.1)
	wt/=np.sum(wt)

	for tid in range(len(tiles)):
		t=tpm[tid]
		m=tiles[tid].copy()

		if m.has_attr('ctf'):
			ctf=m["ctf"]
			fft1=m.do_fft()
			flipim=fft1.copy()
			ctf.compute_2d_complex(flipim,Ctf.CtfType.CTF_SIGN)
			fft1.mult(flipim)
			m=fft1.do_ift()
			
		dxy=sel_dxy[tid].copy()[:,:2]
		dx=dxy[di]/options.xfscale
		dx=np.sum(dx*wt[:,None], axis=0)
		# print(tid, dx)
		dt=Transform()
		dt.set_trans(dx.tolist())

		m.process_inplace("filter.ramp")
		# m.process_inplace("xform",{"transform":dt})
		xf=Transform({"type":"xyz","ytilt":t[3],"xtilt":t[4],"ztilt":t[2]})
		if not options.makeraw:
			xf=dt.inverse()*xf

		dy=(pad//2)-np.cos(t[3]*np.pi/180.)*pad/2
		msk=EMData(pad, pad)
		msk.to_one()
		edge=(sz//10)
		msk.process_inplace("mask.zeroedge2d",{"x0":dy+edge, "x1":dy+edge, "y0":edge, "y1":edge})
		msk.process_inplace("mask.addshells.gauss",{"val1":0, "val2":edge})

		m.mult(msk)
		mp=recon.preprocess_slice(m, xf)
		recon.insert_slice(mp,xf,1)


	threed=recon.finish(True)

	threed.clip_inplace(Region((pad-sz)//2, (pad-sz)//2, (pad-outz)//2, sz, sz, outz))
	threed.process_inplace("filter.lowpass.gauss",{"cutoff_abs":options.filterto})
	jsd.put( [stepx, stepy, threed])
	return 


if __name__ == '__main__':
	main()
	
