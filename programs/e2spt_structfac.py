#!/usr/bin/env python
# Muyuan Chen 2019-06
from future import standard_library
standard_library.install_aliases()
from builtins import range
from EMAN2 import *
import numpy as np
from scipy.optimize import minimize

def main():
	
	usage=" "
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--even", type=str,help="even threed file", default=None)
	parser.add_argument("--label", type=str,help="particles label", default=None)
	parser.add_argument("--sfout", type=str,help="output", default="sf.txt")
	parser.add_argument("--cutoff", type=float,help="cutoff", default=20)
	parser.add_argument("--res", type=float,help="lowpass resolution. default 15", default=15)
	parser.add_argument("--sqrt", action="store_true", default=False ,help="sqrt on structure factor curve. maybe better for high res maps.")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)

	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)
	
	if len(args)==0:
		fname=options.even
		e=EMData(fname)
		o=EMData(fname.replace("even","odd"))
		data=(e+o)/2
		fname=fname.replace('_even','')
		data.write_image(fname)
	else:
		fname=args[0]
		data=e=EMData(fname)

	bxsz=e["nx"]
	apix=e["apix_x"]

	
	#data.process_inplace("mask.soft",{"outer_radius":72,"width":16})
	
	dataf = data.do_fft()
	curve = dataf.calc_radial_dist((data["ny"]//2), 0, 1.0, False)
	curve=np.array([i/dataf["nx"]*dataf["ny"]*dataf["nz"] for i in curve])
	if options.sqrt: curve=np.sqrt(curve)
	
	if options.label:
		print("weighting by fsc...")
		pname="ptcl2d.lst"
		run("e2proclst.py particles/*__{}.hdf --create {}".format(options.label, pname))
		n=EMUtil.get_image_count(pname)
		idx=np.arange(n, dtype=int)
		np.random.shuffle(idx)
		ds=1.0/(apix*bxsz)

		c1d=[]
		for i in idx[:1000].tolist():
			e=EMData(pname, i, True)
			ctf=e["ctf"]
			ctf.bfactor=0
			c=ctf.compute_1d(bxsz,ds,Ctf.CtfType.CTF_ABS)
			c1d.append(c)
		
		c=np.mean(c1d, axis=0)
		if os.path.isfile(pname):
			os.remove(pname)
			
		c[0]=1
		curve=curve/c
	
	curve[0]=0
	curve/=np.max(curve)
	
	#### classic structural factor
	s=np.arange(len(curve))*1./(apix*bxsz)
	cv=[2.80906, 1.97088, 0.71626, 0.44646, 0.17816, -0.03370, -0.27300, -0.43296, -0.47462, -0.49176, -0.51401, -0.50851, -0.50085, -0.51879, -0.50726, -0.44237, -0.46572, -0.41184, -0.37315, -0.36693, -0.38623, -0.36812, -0.38944, -0.44176, -0.49944, -0.59203, -0.67172, -0.70637, -0.75822, -0.82767, -0.80866, -0.79560, -0.79147, -0.77391, -0.75435, -0.74013, -0.71295, -0.67304, -0.63188, -0.59686, -0.56459, -0.53561, -0.52926, -0.51478, -0.52703, -0.54996, -0.56983, -0.59393, -0.61916, -0.64065, -0.65594, -0.66507, -0.67619, -0.69587, -0.72263, -0.74979, -0.77228, -0.79427, -0.81728, -0.84210, -0.86782, -0.88952, -0.90666, -0.92398, -0.93935, -0.95353, -0.96825, -0.98245, -0.99630, -1.00828, -1.01905, -1.02951,-1.04015, -1.04975, -1.05807, -1.06691, -1.07601, -1.08674, -1.09222, -1.09494, -1.09815, -1.10561, -1.11427, -1.11832, -1.11867, -1.11744, -1.12003, -1.12583, -1.13025, -1.13495, -1.13707, -1.13804, -1.14301, -1.14933, -1.14846, -1.14018, -1.12828, -1.11983, -1.12223]
	cvy=np.power(10,np.array(cv))
	cvx=np.arange(len(cv), dtype=float)/200.+0.002

	cvyp=np.interp(s, cvx, cvy)
	cvyp=cvyp/np.max(cvyp[1:])
	
	
	cutoff=1./options.cutoff
	def fit_curve(pm, fitting=True):
		cc,b1,b2=pm
		b1=abs(b1); b2=abs(b2)
		cut=1./cc
		mult=np.exp((s)*b1)
		mult[s>cut]*=np.exp((-s[s>cut])*b2)/np.exp(-cut*b2)
		co=curve*mult
		err=np.log(1+cvyp)-np.log(1+co)
		err=np.mean(abs(err[s>cutoff]))
		if fitting:
			return err
		else:
			return co

	init=[10,10,10]
	co=fit_curve(init, False)
	from scipy.optimize import minimize
	res=minimize(fit_curve, init, 
				method='Nelder-Mead',options={'ftol': 1e-4, 'disp': False, "maxiter":30})
	co=fit_curve(res.x, False)

	np.savetxt(options.sfout,np.vstack([s,co]).T)
	
	sf=XYData()
	sf.set_xy_list(s.tolist(), co.tolist())
	# datasf=data.copy()
	datasf=data.process("filter.setstrucfac",{"apix":data["apix_x"],"strucfac":sf})
	datasf.process_inplace("filter.lowpass.gauss", {"cutoff_freq":1./options.res})
	datasf.write_image(fname[:-4]+"sf.hdf")
	
	print("Output written to {}".format(fname[:-4]+"sf.hdf"))
	
	E2end(logid)
	
def run(cmd):
	print(cmd)
	launch_childprocess(cmd)
	
	
if __name__ == '__main__':
	main()
	
