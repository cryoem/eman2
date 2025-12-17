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
	parser.add_argument("--scale", type=float,help="shrink factor for particles", default=-1)
	parser.add_argument("--maxstd", type=float,help="skip particles with 2d translation beyond N std", default=3)	
	parser.add_argument("--local", action="store_true",help="also fit local motion.", default=False)


	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)
		
	path=options.path

	print("Gathering metadata...")
	info3d=load_lst_params(f"{path}/particle_info_3d.lst")

	for i in range(99,0,-1):
		fm=f"{path}/aliptcls2d_{i:02d}.lst"
		if os.path.isfile(fm):
			lst=load_lst_params(fm, range(10))
			options.loadali2d=fm
			break
	print("using 2d alignment from {}".format(options.loadali2d))

	lst2d=load_lst_params(options.loadali2d)
	for l in lst2d:
		l["src"]=l["src"].split('/')[1].split('__')[0]
		
	fnames=np.unique([base_name(l["src"]) for l in info3d])
	print("load {} particles from {} tomograms".format(len(info3d), len(fnames)))
	
	if options.scale>0:
		scale=options.scale
	else:
		e=EMData(info3d,0,True)
		apix=e["apix_x"]
		js=dict(js_open_dict(info_name(fnames[0])))
		scale=apix/js["apix_unbin"]
		
	######################
	######################
	
	#alldata={}
	for fname in fnames:
		info=info_name(fname)
		js=js_open_dict(info)
		tltpm=np.array(js["tlt_params"])[:,:5]
		
		lst_sel=[l for l in lst2d if fname == l["src"]]
		tlt=[l["tilt_id"] for l in lst_sel]
		tlt=np.sort(np.unique(tlt))
		
		pids=[l["ptcl3d_id"] for l in lst_sel]
		pids=np.sort(np.unique(pids))
		
		trans=[]
		for tid in range(len(tltpm)):
			lst_tlt={l["ptcl3d_id"]:l for l in lst_sel if l["tilt_id"]==tid}
			ts=[lst_tlt[i]["dxf"].get_trans() if i in lst_tlt else [0,0,-100] for i in pids]
			
			trans.append(ts)
		
		trans=np.array(trans)
		
		info3d_sel=[info3d[i] for i in pids]
		coord=np.array([i["coord"] for i in info3d_sel])
		print(fname, len(lst_sel), trans.shape, coord.shape, tltpm.shape)
		
		if options.local:
			cfs=[]
			for tid in range(len(tltpm)):
				ts=trans[tid,:,:2].copy()*scale
				d=np.linalg.norm(ts, axis=1)
				#print(ts.shape, d.shape)
				t=np.mean(d)+np.std(d)*options.maxstd
				goodi=d<t
				ts=ts[goodi]
				cd=coord[goodi]
				A = np.vstack((cd[:,0], cd[:,1], np.ones(len(cd)))).T
				coeffs, _,_,_=np.linalg.lstsq(A, ts, rcond=None)
				ts_fit=np.dot(A, coeffs)
				cfs.append(coeffs)
				
			cfs=np.array(cfs).reshape((len(cfs), -1))
			xf2=tltpm.copy()
			xf2=np.hstack([xf2, cfs])
			#print(xf2.shape)
				
			js["tlt_params"]=xf2.copy()
			js.close()
			
		else:
			ts=trans[...,:2].copy()
			d=np.linalg.norm(ts, axis=2)
			t=np.mean(d, axis=1)+np.std(d, axis=1)*options.maxstd
			ts=np.array([np.mean(ts[i, d[i]<t[i]],axis=0) for i in range(len(ts))])
			ts[np.isnan(ts)]=0
			
			xf2=tltpm.copy()
			xf2[:,[0,1]]+=ts*scale
			#xf2=np.hstack([np.arange(nimg)[:,None], xf2])
			js["tlt_params"]=xf2.copy()
			js.close()
	
	E2end(logid)


if __name__ == '__main__':
	main()
	
