#!/usr/bin/env python
from EMAN2 import *
from EMAN2jsondb import JSTask
import numpy as np
from scipy.optimize import minimize
from scipy.ndimage import map_coordinates

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """
"""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	
	parser.add_argument("--path",type=str,help="Path to a folder where results should be stored, following standard naming conventions",default="subtlt_00")
	parser.add_argument("--ref",type=str,help="reference map",default=None)
	parser.add_argument("--iter",type=int,help="Iteration number within path. Default = start a new iteration",default=0)
	parser.add_argument("--maxres",type=float,help="Maximum resolution to consider in alignment (in A, not 1/A)",default=0)
	parser.add_argument("--minres",type=float,help="Minimum resolution to consider in alignment (in A, not 1/A)",default=0)
	parser.add_argument("--maxshift",type=float,help="max shift in pixel. default default box size/6",default=-1)
	parser.add_argument("--order",type=int,help="order for polynomial fit. default is 2",default=2)
	parser.add_argument("--l2",type=float,help="l2 regularization for polynomial fit. default is 1e-4",default=1e-4)
	
	#parser.add_argument("--aliptcls2d",type=str,help="optional aliptcls input. the program can start search from the position from last run.",default="")
	parser.add_argument("--aliptcls3d",type=str,help="optional aliptcls input.",default="")

	parser.add_argument("--parallel", type=str,help="Thread/mpi parallelism to use", default="thread:4")
	parser.add_argument("--flatten",action="store_true",help="flatten power spetrum before alignment.",default=False)
	parser.add_argument("--debug",action="store_true",help="for testing.",default=False)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higher number means higher level of verboseness")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)

	(options, args) = parser.parse_args()
	logid=E2init(sys.argv, options.ppid)
	
	options.info2dname="{}/particle_info_2d.lst".format(options.path)
	options.info3dname="{}/particle_info_3d.lst".format(options.path)
	
	lst=load_lst_params(options.info2dname)
	n=len(lst)
	fnames=np.unique([l["src"] for l in lst])
	
	tasks=[]
	for fm in fnames:		
		ii=[i for i,l in enumerate(lst) if l["src"]==fm]
		tid=np.unique([lst[i]["tilt_id"] for i in ii])
		for t in tid:
			it=[i for i in ii if lst[i]["tilt_id"]==t]
			tasks.append(it)
			 
	print("{:d} tomograms, {:d} tasks total.".format(len(fnames), len(tasks)))
	
	from EMAN2PAR import EMTaskCustomer
	etc=EMTaskCustomer(options.parallel, module="e2spt_subtlt_global.SptAlignTask")
	
	num_cpus = etc.cpu_est()
	if options.debug:
		tasks=tasks[:num_cpus*2]
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
	
	output=[None]*n
	for i in tids:
		rets=etc.get_results(i)[1]
		for r in rets:
			output[r[0]]=r[1]
		
	del etc
	
	fm="{}/aliptcls2d_{:02d}.lst".format(options.path, options.iter)
	save_lst_params(output, fm)
	
	E2end(logid)

class SptAlignTask(JSTask):
	
	
	def __init__(self, data, options):
		
		JSTask.__init__(self,"SptAlign",data,{},"")
		self.options=options
	
	
	def execute(self, callback):
		
		def test_trans(tc):
			tc=tc.reshape((-1,2))
			p_new=np.dot(xy1, tc)+maxshift/2
			pp=np.vstack([np.arange(len(p_new)), p_new[:,1], p_new[:,0]])
			c=map_coordinates(cf, pp, order=1)
			c=-np.mean(c)
			c+=np.sum(tc[:-2]**2)*options.l2
			return c
		
		callback(0)
		options=self.options
		
			
		### load references
		refnames=[options.ref[:-4]+f"_{eo}.hdf" for eo in ["even", "odd"]]
			
		ref=EMData(refnames[0],0,True)
		ny=ref["ny"]
		apix=ref["apix_x"]
		
		refs=[]
		for r in refnames:
			ref=EMData(r,0)
			refft=ref.do_fft()
			refft.process_inplace("xform.phaseorigin.tocenter")
			refft.process_inplace("xform.fourierorigin.tocenter")
			refs.append(refft)
				
		#### resolution range
		if options.minres>0:
			minp=ceil(ny*ref["apix_x"]/options.minres)
		else:
			minp=2
			
		if options.maxres>0:
			maxp=ceil(ny*ref["apix_x"]/options.maxres)
		else:
			maxp=ny//2			
		
		if options.maxshift<0:
			options.maxshift=ny//6
		maxshift=options.maxshift
		
		sf=XYData()
		sf.set_xy_list([0.,1.], [1.,1.])
		
		if options.debug:
			print(f"minp: {minp}, maxp : {maxp}, maxshift: {maxshift}")
		
		rets=[]
		for iid, di in enumerate(self.data):
			info2d=load_lst_params(options.info2dname, di)
			assert info2d[0]["src"]==info2d[-1]["src"]
			assert info2d[0]["tilt_id"]==info2d[-1]["tilt_id"]
			
			if options.debug:
				print("tomo: {}, tilt {}, {} particles.".format(base_name(info2d[0]["src"]), info2d[0]["tilt_id"], len(di)))
			
			i3=[i["idx3d"] for i in info2d]
			info3d=load_lst_params(options.info3dname, i3)
			ali3d=load_lst_params(options.aliptcls3d, i3)
			
			xf3d=[a["xform.align3d"].inverse() for a in ali3d]
			
			if options.debug:  print("reading images...")
			imgs=[]
			for l in info2d:
				m=EMData(l["src"], l["idx"])
				pad=m['ny']
				m=m.get_clip(Region( (pad-ny)//2,(pad-ny)//2,ny,ny))
				if options.flatten:
					m.process_inplace("filter.setstrucfac", {"strucfac":sf})
					
				m.process_inplace("filter.lowpass.gauss",{"cutoff_pixels":maxp})
				m.process_inplace("filter.highpass.gauss",{"cutoff_pixels":minp})
				m.process_inplace("normalize.edgemean")
				
				a=m.do_fft()
				a.process_inplace("xform.phaseorigin.tocenter")
				a.process_inplace("xform.fourierorigin.tocenter")
				a.process_inplace("xform.fourierorigin.tocenter")
				imgs.append(a)
				
				
			if options.debug:  print("projecting references...")
			xfs=[m["xform.projection"]*x for m,x in zip(imgs, xf3d)]
			cls=[l["class"] for l in info3d]
			pjs=[refs[c].project('gauss_fft',{"transform":x, "returnfft":1}) for x,c in zip (xfs, cls)]
			#print(xfs[0])
			#### cross correlation
			
			ccf=[m.calc_ccf(p) for p,m in zip(pjs, imgs)]
			ccf=[c.process("xform.phaseorigin.tocenter") for c in ccf]
			ccf=[c.get_clip(Region( (ny-maxshift)//2,(ny-maxshift)//2,maxshift,maxshift)) for c in ccf]

			cf=sum(ccf, ccf[0]*0)
			p_global=np.array(cf.calc_max_location())[:2]-maxshift//2
			if options.debug:  print("global shift : ",p_global)
			
			
			cf=np.array([c.numpy().copy() for c in ccf])
			cf/=np.std(cf, axis=(1,2))[:,None,None]
			
			#### optimization

			trans_coef=np.zeros((3*options.order+1,2))
			trans_coef[-1]=p_global.copy()
			
			crd=np.array([i["coord"] for i in info3d])
			if options.debug:  
				print("particle location : {} to {}".format(np.min(crd, axis=0), np.max(crd, axis=0)))
			
			crd/=1000
			
			xy1= np.ones((len(crd),1))
			for i in range(options.order):
				xy1=np.hstack([crd**(i+1), xy1])
			
			res=minimize(test_trans, trans_coef.flatten(),
						method='Powell',options={'ftol': 1e-4, 'disp': False, "maxiter":10})
			
			if options.debug: print(res)
			
			tc_new=res.x.reshape((-1,2))
			p_new=np.dot(xy1, tc_new)
			
			for i, l in enumerate(info2d):
				dxf=Transform({"type":"eman", "tx":float(p_new[i,0]), "ty":float(p_new[i,1])})
				xf=dxf*xfs[i]
				l["xform.projection"]=xf
				l["class"]=cls[i]
				l["ptcl3d_id"]=l.pop("idx3d")
				l["dxf"]=dxf
				
				pj=refs[cls[i]].project('gauss_fft',{"transform":xf, "returnfft":1})
				frc=pj.cmp("frc", imgs[i], {"ampweight":0, "sweight":0, "snrweight":0, "zeromask":0, "pmin":minp, "pmax":maxp})
				l["score"]=frc
				
				rets.append((di[i],l))

			callback((iid+1)*100//len(self.data))
		
		return rets

def make_projs(refs, rid, xfs):
	pjs=[]
	for ci,xf in enumerate(xfs):
		ref=refs[rid[ci]]
		pj=ref.project('gauss_fft',{"transform":xf, "returnfft":1})
		pjs.append(pj)
		
	return pjs


if __name__ == "__main__":
	main()

