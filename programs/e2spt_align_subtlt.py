#!/usr/bin/env python

from past.utils import old_div
from future import standard_library
standard_library.install_aliases()
from builtins import range
from EMAN2 import *
import time
import os
import threading
import queue
from sys import argv,exit
from EMAN2jsondb import JSTask
import numpy as np
from scipy.optimize import minimize

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """Usage: e2spt_align.py [options] <subvolume_stack> <reference>
This program is part of the 'new' hierarchy of e2spt_ programs. It performs one iteration of a classical subtomogram refinement, ie -  aligning particles with missing wedge to a reference in 3-D

The reference may be <volume> or <volume>,<n>
"""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_argument("--threads", default=4,type=int,help="Number of alignment threads to run in parallel on a single computer. This is the only parallelism supported by e2spt_align at present.", guitype='intbox', row=24, col=2, rowspan=1, colspan=1, mode="refinement")
	parser.add_argument("--iter",type=int,help="Iteration number within path. Default = start a new iteration",default=0)
	parser.add_argument("--goldcontinue",action="store_true",help="Will use even/odd refs corresponding to specified reference to continue refining without phase randomizing again",default=False)
	parser.add_argument("--path",type=str,default=None,help="Path to a folder where results should be stored, following standard naming conventions (default = spt_XX)")
	parser.add_argument("--sym",type=str,default="c1",help="Symmetry of the input. Must be aligned in standard orientation to work properly.")
	parser.add_argument("--maxres",type=float,help="Maximum resolution to consider in alignment (in A, not 1/A)",default=0)
	parser.add_argument("--minres",type=float,help="Minimum resolution to consider in alignment (in A, not 1/A)",default=0)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higher number means higher level of verboseness")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--maxshift", type=int, help="maximum shift for subtilt refine",default=8)
	parser.add_argument("--parallel", type=str,help="Thread/mpi parallelism to use", default=None)
	parser.add_argument("--fromscratch",action="store_true",help=".",default=False)
	parser.add_argument("--skipali3d",action="store_true",help="",default=False)
	parser.add_argument("--skipali2d",action="store_true",help="",default=False)
	parser.add_argument("--debug",action="store_true",help=".",default=False)
	parser.add_argument("--plst",type=str,default=None,help="list of 2d particle with alignment parameters. will reconstruct before alignment.")
	parser.add_argument("--smooth",type=float,help="smooth local motion by this factor. smoother local motion with larger numbers",default=-1)
	parser.add_argument("--smoothN",type=float,help="number of neighboring particles used for smoothing. default 15",default=15)

	(options, args) = parser.parse_args()
	
	if options.path == None:
		options.path=num_path_new("spt_")

	if options.iter<=0 :
		fls=[int(i[15:17]) for i in os.listdir(options.path) if i[:15]=="particle_parms_" and str.isdigit(i[15:17])]
		if len(fls)==0 : options.iter=1
		else: options.iter=max(fls)+1
		
	if options.parallel==None:
		options.parallel="thread:{}".format(options.threads)
		

	# file may be "name" or "name,#"
	reffile=args[1].split(",")[0]
	try: refn=int(args[1].split(",")[1])
	except: refn=0
	
	NTHREADS=max(options.threads+1,2)		# we have one thread just writing results

	logid=E2init(sys.argv, options.ppid)
	refnames=[reffile, reffile]

	if options.goldcontinue:
		ref=[]
		try:
			refnames=[reffile[:-4]+"_even.hdf", reffile[:-4]+"_odd.hdf"]
			ref.append(EMData(refnames[0],0))
			ref.append(EMData(refnames[1],0))
			
		except:
			print("Error: cannot find one of reference files, eg: ",EMData(reffile[:-4]+"_even.hdf",0))
			
	n=-1
	tasks=[]
	readjson=False
	if args[0].endswith(".lst") or args[0].endswith(".hdf"):
		nptcl=EMUtil.get_image_count(args[0])
		tasks.extend([(args[0],i,refnames, i%2) for i in range(nptcl)])
	
	elif args[0].endswith(".json"):
		js=js_open_dict(args[0])
		readjson=True
		jsinput=dict(js)
		keys=sorted(js.keys())
		nptcl=len(keys)
		for k in keys:
			src, ii=eval(k)
			dic=js[k]
			xf=dic["xform.align3d"]
			tasks.append([src, ii, refnames, ii%2, xf])
	
	
	if options.plst:
		xinfo=[[] for i in range(nptcl)]
		xfkey=["type","alt","az","phi","tx","ty","tz","alpha","scale"]
		lst=LSXFile(options.plst)
		for i in range(lst.n):
			l=lst.read(i)
			dc=eval(l[2])
			pid=dc["pid"]
			dxf=Transform({k:dc[k] for k in dc.keys() if k in xfkey})
			xinfo[pid].append([dc["tid"],i, dxf])
		
		lst=None
		for i,xn in enumerate(xinfo):
			sid=np.argsort([x[0] for x in xn])
			xinfo[i]=[[xn[s][1],xn[s][2]] for s in sid]
			
		for i in range(nptcl):
			ii=tasks[i][1]
			tasks[i].append(xinfo[ii])

	from EMAN2PAR import EMTaskCustomer
	etc=EMTaskCustomer(options.parallel, module="e2spt_align_subtlt.SptAlignTask")
	
	num_cpus = etc.cpu_est()
	options.nowtime=time.time()
	if options.debug:
		tasks=tasks[:num_cpus*4]
	print("{} jobs on {} CPUs".format(len(tasks), num_cpus))
	njob=num_cpus
	
	tids=[]
	for i in range(njob):
		t=tasks[i::njob]
		task=SptAlignTask(t, options)
		if options.debug:
			ret=task.execute(print)
			return 
		tid=etc.send_task(task)
		tids.append(tid)

	while 1:
		st_vals = etc.check_task(tids)
		if -100 in st_vals:
			print("Error occurs in parallelism. Exit")
			return
		E2progress(logid, np.mean(st_vals)/100.)
		
		if np.min(st_vals) == 100: break
		time.sleep(5)
	
	output=[]
	for i in tids:
		rets=etc.get_results(i)[1]
		for ret in rets:
			fsp,n,dic=ret
			output.append([fsp, n, dic])
		
	
	angs={}
	data=[]
	
	for out in output:
		fsp,n,dc=out
		angs[(fsp, n)]={
			"xform.align3d":dc["xform.align3d"], 
			"score":dc["score"]}
		
		ts=dc.pop("transscr")
		for i, ii in enumerate(dc["imgidx"]):
			dxf=dc["imgxfs"][i].get_params("eman")
			#dxf["score"]=float(dc["imgscore"][i])
			dxf["pid"]=n; dxf["tid"]=i; dxf["class"]=n%2
			tsmult=dc["tsmult"]
			data.append({
				"src":dc["imgsrc"], "srci":ii, "pid":n, "tid":i,
				"xf":dxf, "coord":dc["coord"], "tsscr":ts[i]
				})
	
	m=options.maxshift
	trans=np.indices((m*2+1, m*2+1)).reshape((2,-1)).T-m
	alldtx=np.zeros((len(data), 3))
	fnames=np.unique([d["src"] for d in data])
	maxtid=np.max([d["tid"] for d in data])+1
	print("Smoothing trajectories...")
	for fname in fnames:
		print(fname)

		for tid in range(maxtid):
			idx=[i for i,d in enumerate(data) if d["src"]==fname and d["tid"]==tid]
			if len(idx)==0:
				continue
			
			coord=np.array([data[i]["coord"] for i in idx])
			scrs=np.array([data[i]["tsscr"] for i in idx])
			
			if options.smooth<=0:
				s=np.argmin(scrs, 1)
				sv=np.min(scrs, axis=1)
				newdtx=np.hstack([trans[s], sv[:,None]])
				if tid==0:print("  Skip smoothing")
			else:	
				newdtx=[]
				
				for i, crd in enumerate(coord):
					dst=np.linalg.norm(coord-crd, axis=1)
					srtid=np.argsort(dst)[:options.smoothN]
					dst=dst[srtid]
					wt=np.exp(-dst/options.smooth)

					scr=scrs[srtid].copy()
					scr=np.sum(scr*wt[:,None], axis=0)
					pos=trans[np.argmin(scr)].copy()
					
					pos=pos.astype(float)*tsmult


					newdtx.append([pos[0], pos[1], np.min(scr)])

				sys.stdout.write("\r   {}: {:.4f}".format(tid, np.mean(scr)))
				sys.stdout.flush()
				newdtx=np.array(newdtx)
			
			alldtx[idx]=newdtx
		print()
		 
	for i,d in enumerate(data):
		dxf=d["xf"]
		dxf["tx"]+=alldtx[i][0]
		dxf["ty"]+=alldtx[i][1]
		dxf["score"]=alldtx[i][2]
		d["xf"]=dxf
	
	fm="{}/aliptcls_{:02d}.lst".format(options.path, options.iter)
	if os.path.isfile(fm): os.remove(fm)
	lout=LSXFile(fm, False)
	for i,d in enumerate(data):
		lout.write(-1, d["srci"], d["src"], d["xf"])
	lout=None
	
	f="{}/aliptcls_ts_{:02d}.hdf".format(options.path,options.iter)
	if os.path.isfile(f): os.remove(f)
	t=np.array([d["tsscr"] for d in data])
	m=from_numpy(t).copy()
	m.write_image(f)
	#print(f, t.shape)
	
	out="{}/particle_parms_{:02d}.json".format(options.path,options.iter)
	if os.path.isfile(out):
		os.remove(out)
	js=js_open_dict(out)
	js.update(angs)
	js.close()

	del etc

	E2end(logid)

class SptAlignTask(JSTask):
	
	
	def __init__(self, data, options):
		
		JSTask.__init__(self,"SptAlign",data,{},"")
		self.options=options
	
	
	def execute(self, callback):
		
		def testxf(x, return_xf=False):
			
			txf=Transform({"type":"eman", "tx":x[0], "ty":x[1], "tz":x[2],"alt":x[3], "az":x[4], "phi":x[5]})
			xfs=[p*txf for p in pjxfsmall]
			pjs=[refsmall.project('gauss_fft',{"transform":x, "returnfft":1}) for x in xfs]

			fscs=np.array([m.calc_fourier_shell_correlation(p) for m,p in zip(imgsmallpjs, pjs)])
			fscs=fscs.reshape((len(fscs), 3, -1))[:,1,:]
			scr=-np.mean(fscs[:, minp:maxp])
			
			if return_xf:
				ss=-np.mean(fscs[:, minp:maxp], axis=1)
				return scr, txf, xfs, ss
			else:
				return scr
		
		
		callback(0)
		options=self.options
		rets=[]
		
		refnames=self.data[0][2]
		ref=EMData(refnames[0],0,True)
		by=ref["ny"]
		ny=good_size(by*1.4)
		apix=ref["apix_x"]
		
		refs=[]
		for r in refnames:
			ref=EMData(r,0)
			ref=ref.get_clip(Region((by-ny)/2, (by-ny)/2,(by-ny)/2, ny, ny,ny))
			refft=ref.do_fft()
			refft.process_inplace("xform.phaseorigin.tocenter")
			refft.process_inplace("xform.fourierorigin.tocenter")
			refs.append(refft)
		
		#### for initial coarse alignment
		astep=7.5
		sym=Symmetries.get(options.sym)
		xfcrs=sym.gen_orientations("saff",{"delta":astep,"phitoo":astep,"inc_mirror":1})
		
		#### for later iterative alignment
		ssrg=2**np.arange(4,12, dtype=int)
		ssrg[:2]=ssrg[:2]*3/2
		ssrg=ssrg.tolist()
		ssrg=[24]+ssrg
		
		#### resolution range
		if options.minres>0:
			minp=ceil(ny*ref["apix_x"]/options.minres)
		else:
			minp=2
			
		if options.maxres>0:
			maxp=ceil(ny*ref["apix_x"]/options.maxres)
			maxy=good_size(maxp*3)
			maxy=int(min(maxy, ny))
		else:
			maxy=ny
			maxp=int(ny*.4)
			
			
		if options.debug: print("max res: {:.2f}, max box size {}".format(options.maxres, maxy))
		for di,data in enumerate(self.data):
			
			#### prepare inputs
			fsp=data[0]
			fid=data[1]
			ref=refs[data[3]]
			if len(data)>4:
				lastxf=data[4].inverse()
			else:
				lastxf=None
			
			## 3D particle
			img=EMData(fsp, fid)
			img.process_inplace("filter.highpass.gauss",{"cutoff_pixels":4})
			img.process_inplace("normalize.edgemean")
			img.process_inplace("mask.soft",{"outer_radius":-1})
			img=img.get_clip(Region((by-ny)/2, (by-ny)/2,(by-ny)/2, ny, ny,ny))
			
			img=img.do_fft()
			img.process_inplace("xform.phaseorigin.tocenter")
			img.process_inplace("xform.fourierorigin.tocenter")
			
			imgsrc=img["class_ptcl_src"]
			imgidx=img["class_ptcl_idxs"]
			imgcoord=img["ptcl_source_coord"]
			
			## projection transforms for the 2d particles
			pjxfs=[]
			if options.plst:
				#### a 2d particle list with transforms is provided
				## reconstruct 3d particle from 2d paticles
				dts=[]
				pjxfs00=[]
				xinfo=data[5]
				parms = {"size":[ny]*3,"mode":"trilinear","sym":"c1"}
				recon=Reconstructors.get("fourier", parms)
				recon.setup()
				for i, txf in xinfo:
					e=EMData(options.plst, i)
					e-=e.get_edge_mean()
					e.clip_inplace(Region((e["nx"]-ny)//2, (e["ny"]-ny)//2, ny, ny))
					xf=txf*(lastxf.inverse())
					pxf=e["xform.projection"]
					dt=xf.get_trans()-pxf.get_trans()
					dts.append(dt)
					ts=Transform(xf)
					ts.set_rotation({"type":"eman"})
					ts.invert() 
					e=recon.preprocess_slice(e,ts)
					recon.insert_slice(e,xf,1)
					pjxfs.append(xf)
					pjxfs00.append(pxf)
				
				img=recon.finish(True)
				img.process_inplace("normalize.circlemean", {"radius":by*.5,"width":4})
				img.process_inplace("mask.soft",{"outer_radius":by*.5,"width":8})
				#img["xform.align3d"]=lastxf.inverse()
				#img=img.get_clip(Region((ny-by)/2, (ny-by)/2,(ny-by)/2, by, by,by))
				#img.write_image("tmpptcls.hdf",fid)
				#continue
			
				img=img.do_fft()
				img.process_inplace("xform.phaseorigin.tocenter")
				img.process_inplace("xform.fourierorigin.tocenter")
				dts=np.array(dts)[:,:2]
			else:
				#### otherwise, read projection transform from header
				dts=np.zeros((len(imgidx), 2))
				for i in imgidx: 
					m=EMData(imgsrc, i, True)
					pjxfs.append(m["xform.projection"])
					
				pjxfs00=pjxfs

			#### make 2d particle projections from 3d particles
			imgpjs=[]
			for pxf in pjxfs:
				m=img.project('gauss_fft',{"transform":pxf, "returnfft":1})
				m.process_inplace("xform.fourierorigin.tocenter")
				imgpjs.append(m)
			
			#### start from coarse alignment / do refine search around previous solution
			if options.fromscratch:
				curxfs=[]
				npos=32
				ifirst=0
			else:
				curxfs=[lastxf]
				npos=1
				ifirst=len(ssrg)-1
				
			#### optionally skip 3d alignment. need to prepare results for subsequent 2d alignment
			if options.skipali3d:
				ilast=-1
				ss=maxy
				refsmall=ref.get_clip(Region(0,(ny-ss)//2, (ny-ss)//2, ss+2, ss, ss))
				newscore=[-1]
				imgsmallpjs=[]
				for pj in imgpjs:
					m=pj.get_clip(Region(0,(ny-ss)//2, ss+2, ss))
					m.process_inplace("xform.fourierorigin.tocenter")
					imgsmallpjs.append(m)
					
			else:
				ilast=len(ssrg)
			
			#### 3d alignment loop. increase fourier box size and reduce solutions every iteration
			if options.debug: print("Align particle ( {}, {} )".format(fsp, fid))
			for si in range(ifirst, ilast):
				ss=ssrg[si]
				if ss>=maxy: 
					ss=maxy
			
				mxsft=ss//6
				astep=89.999/floor((np.pi/(3*np.arctan(2./ss)))) ### copied from spt tree aligner
				newxfs=[]
				score=[]
				if options.debug: print("size {}, solution {}, astep {:.1f}, xmult {:.1f}".format(ss, npos, astep, ny/ss))

				refsmall=ref.get_clip(Region(0,(ny-ss)//2, (ny-ss)//2, ss+2, ss, ss))	
				
				if si==0:
					## coarse alignment using 3d particles
					## need to shift origin to rotate 3d volume
					imgsmall=img.get_clip(Region(0,(ny-ss)//2, (ny-ss)//2, ss+2, ss, ss))
					refsmall.process_inplace("xform.fourierorigin.tocenter")
					imgsmall.process_inplace("xform.fourierorigin.tocenter")
					
					for xf0 in xfcrs:
						xf=Transform(xf0)
						refrot=refsmall.process("xform", {"transform":xf})
						ccf=imgsmall.calc_ccf(refrot)
						pos=ccf.calc_max_location_wrap(mxsft, mxsft, mxsft)
						#xf.translate(pos[0], pos[1],pos[2])
				
						#refrot=refsmall.process("xform", {"transform":xf})
						refrot=refrot.process("xform", {"tx":pos[0], "ty":pos[1], "tz":pos[2]})
						scr=imgsmall.cmp("ccc.tomo.thresh", refrot)
						score.append(scr)
						
						newxfs.append(xf)
						if options.debug: 
							sys.stdout.write("\r {}/{}    {:.3f}  ".format(len(newxfs), len(xfcrs), score[-1]).ljust(40))
							sys.stdout.flush()
					
				else:
					## refine orientation by projections
					## make small projections
					imgsmallpjs=[]
					for pj in imgpjs:
						m=pj.get_clip(Region(0,(ny-ss)//2, ss+2, ss))
						m.process_inplace("xform.fourierorigin.tocenter")
						imgsmallpjs.append(m)
					
					pjxfsmall=[Transform(p) for p in pjxfs]
					for p in pjxfsmall:
						p.set_trans(p.get_trans()*ss/ny)
						
					for xf in curxfs:
						x=xf.get_params("eman")
						x0=[x["tx"]*ss/ny, x["ty"]*ss/ny, x["tz"]*ss/ny, x["alt"], x["az"], x["phi"]]
						
						simplex=np.vstack([[0,0,0,0,0,0], np.eye(6)])
						simplex[4:]*=astep
						
						res=minimize(testxf, x0,  method='Nelder-Mead', options={'ftol': 1e-3, 'disp': False, "maxiter":50,"initial_simplex":simplex+x0})
						
						x=res.x
						txf=Transform({"type":"eman", "tx":x[0], "ty":x[1], "tz":x[2],"alt":x[3], "az":x[4], "phi":x[5]})
						newxfs.append(txf)
						score.append(float(res.fun))
						
						if options.debug: 
							scr0 = testxf(x0)
							sys.stdout.write("\r {}/{}    {:.3f}    {:.3f}     ".format(len(newxfs), len(curxfs), scr0, score[-1]))
							sys.stdout.flush()
				
				if options.debug: print()
					
				curxfs=newxfs

				newxfs=[]
				newscore=[]
				idx=np.argsort(score)
				for i in idx:
					dt=[(x.inverse()*curxfs[i]).get_params("spin")["omega"] for x in newxfs]
					if len(dt)==0 or np.min(dt)>astep*4:
						newxfs.append(curxfs[i])
						newscore.append(score[i])
					if len(newxfs)>=npos:
						break
					
				for xi,xf in enumerate(newxfs):
					## translation between -ny//2 and ny//2
					x=np.array(xf.get_trans()*ny/ss)
					x=x%ny
					x=x-ny*(x>ny//2)
					
					xf.set_trans(x.tolist())
						
					if options.debug: print(newscore[xi], xf)
					
				npos=max(1, npos//2)
				curxfs=newxfs
				if ss>=maxy:
					break
				
			
			##############
			
			xfout=Transform(curxfs[0])
			score=newscore[0]
			
			thisxf=Transform(curxfs[0])
			thisxf.set_trans(thisxf.get_trans()*maxy/ny)
			
			pjxfsmall=[Transform(p) for p in pjxfs]
			for p in pjxfsmall:
				p.set_trans(p.get_trans()*maxy/ny)
				
			imgxfs=[p*thisxf for p in pjxfsmall]
			pjs=[refsmall.project('gauss_fft',{"transform":x, "returnfft":1}) for x in imgxfs]
			m=options.maxshift
			trans=np.indices((m*2+1, m*2+1)).reshape((2,-1)).T-m
			if options.skipali2d:
				scr=np.zeros((len(trans), len(pjxfs)))+1
				fscs=np.array([m.calc_fourier_shell_correlation(p) for m,p in zip(imgsmallpjs, pjs)])
				fscs=fscs.reshape((len(fscs), 3, -1))[:,1,minp:maxp]
				scr[len(trans)//2]=-np.mean(fscs, axis=1)
			else:
				scr=[]

				for t in trans.tolist():
					pjts=[p.process("xform",{"tx":t[0],"ty":t[1]}) for p in pjs]
					
					fscs=np.array([m.calc_fourier_shell_correlation(p) for m,p in zip(imgsmallpjs, pjts)])
					fscs=fscs.reshape((len(fscs), 3, -1))[:,1,minp:maxp]
					scr.append(-np.mean(fscs, axis=1))
					
				scr=np.array(scr)
			#print(scr[len(trans)//2])
				
			dts=dts*maxy/ny
			s0=scr.reshape((m*2+1, m*2+1, -1))
			s1=np.zeros_like(scr)
			dt=np.round(dts).astype(int)
			for i in range(scr.shape[1]):
				s=do_trans(s0[:,:,i], dt[i])
				s=s.flatten()
				s1[:,i]=s
			
			scr=s1.copy()
			
			imgxfs=[p*xfout for p in pjxfs00]
			
			c={"xform.align3d":xfout.inverse(), "score":score, "imgsrc":imgsrc, "imgidx":imgidx, "imgxfs":imgxfs, "transscr":scr.T, "coord":imgcoord, "tsmult":float(ny/maxy)}
						
			rets.append((fsp,fid,c))
			
			if options.debug:
				
				if lastxf:
					x=Transform(lastxf)
					x.set_trans(x.get_trans()*maxy/ny)
					x=x.get_params("eman")
					c=[x["tx"], x["ty"], x["tz"], x["alt"], x["az"], x["phi"]]
					s0 = testxf(c)
				else:
					s0=0
				
				x=xfout.get_params("eman")
				c=[x["tx"], x["ty"], x["tz"], x["alt"], x["az"], x["phi"]]
				
				ti=thisxf.inverse()
				pjxfs=[p*ti for p in imgxfs]
				x=thisxf.get_params("eman")
				x0=[x["tx"], x["ty"], x["tz"], x["alt"], x["az"], x["phi"]]
				s1 = testxf(x0)
				print(fid,  np.round(c,2),s0, score,s1)
				
				print('#############')
				exit()
				
			else:
				callback(len(rets)*100//len(self.data))

		
		return rets
		

def do_trans(a, x):
	n=len(a)
	b=np.zeros_like(a)
	x=np.array(x)
	i=np.array([[0,n],[0,n]])+x[:,None]
	i=np.clip(i, 0, n).astype(int)
	j=np.array([[0,n],[0,n]])-x[:,None]
	j=np.clip(j, 0, n).astype(int)
	b[i[0,0]:i[0,1],i[1,0]:i[1,1]]=a[j[0,0]:j[0,1], j[1,0]:j[1,1]]
	return b
		

if __name__ == "__main__":
	main()

