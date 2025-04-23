#!/usr/bin/env python
# Muyuan Chen 2020-07
from EMAN2 import *
import numpy as np
from multiprocessing import Pool
#import threading
#import queue

def local_ccc(c):
	ef=em.process("filter.bandpass.gauss", {"cutoff_abs":.05, "center":c})
	of=om.process("filter.bandpass.gauss", {"cutoff_abs":.05, "center":c})
	#ef=em.process("filter.bandpass.tophat", {"apix":1, "high_cutoff_frequency":c+.05, "low_cutoff_frequency":c-.05})
	#of=om.process("filter.bandpass.tophat", {"apix":1, "high_cutoff_frequency":c+.05, "low_cutoff_frequency":c-.05})
	ef.mult(mask)
	of.mult(mask)
	scr=[]
	l=msk["ny"]
	#print(c)
	# vout=EMData(ny,ny,ny)
	for ii in ind.tolist():
		x,y,z=ii
		efc=ef.get_clip(Region(x-l//2,y-l//2,z-l//2,l,l,l))
		ofc=of.get_clip(Region(x-l//2,y-l//2,z-l//2,l,l,l))
		s=-efc.cmp("ccc", ofc,{"mask":msk})
		scr.append(s)
		
	#jsd.put(c, np.array(scr))
	return np.array(scr)


def main():
	
	usage=" "
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--even", type=str,help="", default="")
	parser.add_argument("--odd", type=str,help="", default="")
	parser.add_argument("--output", type=str,help="", default="")
	parser.add_argument("--setsf", type=str,help="", default="")
	parser.add_argument("--mask", type=str,help="", default="")
	parser.add_argument("--step", type=int,help="", default=4)
	parser.add_argument("--winsize", type=int,help="", default=17)
	parser.add_argument("--cut", type=float,help="", default=.2)
	parser.add_argument("--sym", type=str,help="", default="c1")	
	parser.add_argument("--overwrite", action="store_true", default=False ,help="overwrite even/odd input")
	parser.add_argument("--gauss", action="store_true", default=False ,help="gauss instead of tophat")
	parser.add_argument("--ppid", type=int,help="", default=-1)
	parser.add_argument("--fscvol", type=str,help="use volume as fsc input", default=None)	

	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)
	
	global em,om,mask,ind,msk
	if options.odd=="":
		options.odd=options.even.replace("_even", "_odd")
	em=EMData(options.even)
	om=EMData(options.odd)
	if options.setsf:
		sf=XYData()
		sf.read_file(options.setsf)
		em.process_inplace("filter.setstrucfac", {"strucfac":sf})
		om.process_inplace("filter.setstrucfac", {"strucfac":sf})
	if options.mask:
		mask=EMData(options.mask)
	else:
		mask=em.copy()
		mask.to_one()
		
	ny=em["ny"]
	step=options.step
	ns=ny//step-1
	ind=np.indices((ns,ns,ns)).reshape((3,-1)).T+1
	ind=ind*step
	# print(ind)
	lnx=options.winsize
	msk=EMData(lnx, lnx, lnx)
	msk.to_one()
	msk.process_inplace("mask.gaussian",{"inner_radius":lnx//6,"outer_radius":lnx//6})
	# vout=EMData(ny,ny,ny)
	cs=(np.arange(20)*.025+.0125)#.tolist()
	if options.fscvol==None:
		pl=Pool(20)
		scrs=pl.map(local_ccc, cs)
		pl.close()
	
		
		scrs=np.array(scrs)
		cut=options.cut
		res=np.sum(np.cumsum(scrs<cut, axis=0)==0, axis=0)#-1
		res[res==len(cs)]=len(cs)-1
		res=cs[res]
		ii=np.indices((ns,ns,ns)).reshape((3,-1))
		fvol=np.zeros((ns,ns,ns))
		fvol[ii[0], ii[1], ii[2]]=np.array(res)
		f=from_numpy(fvol.T.copy())
		f.clip_inplace(Region((ns-ny)//2,(ns-ny)//2,(ns-ny)//2, ny, ny, ny))
		f.scale(step)
		#tx=lnx//4
		#f.translate(tx,tx,tx)
		
	else:
		f=EMData(options.fscvol)
		
	f.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.3})
	f.process_inplace("xform.applysym",{"sym":options.sym})
	f.mult(mask)
	f.set_attr_dict(em.get_attr_dict())
	
	mps=[]
	for mp in [em, om]:
		eout=EMData(ny,ny,ny)
		wt=eout.copy()

		for c in cs:
			m=f.process("threshold.binaryrange", {"low":c-.05, "high":c+.05})
			m.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.3})
			if options.gauss:
				ef=mp.process("filter.lowpass.gauss",{"cutoff_abs":c+.05})
			else:
				ef=mp.process("filter.lowpass.tophat",{"cutoff_abs":c+.05})
			#ef.process_inplace("normalize.circlemean",{"radius":-8,"width":5})
			ef.mult(m)
			eout.add(ef)
			wt.add(m)


		wt.process("math.reciprocal")
		eout.mult(wt)
		eout.process_inplace("xform.applysym",{"sym":options.sym})
		eout.process_inplace("normalize.circlemean",{"radius":-8,"width":5})
		
		eout.mult(mask)
		eout.set_attr_dict(mp.get_attr_dict())
		
		mps.append(eout)
		
	if options.overwrite:
		mps[0].write_image(options.even)
		mps[1].write_image(options.odd)
	else:
		mps[0].write_image(options.even[:-4]+"_out.hdf")
		mps[1].write_image(options.odd[:-4]+"_out.hdf")
	mp=mps[0]+mps[1]
	mp.mult(.5)
	if options.output=="":
		options.output=options.even[:-4].replace("_even","")+".hdf"
		
	mp.write_image(options.output)
	
	f.mult(1./mp["apix_x"])
	f.write_image(options.output.replace("threed","fscvol"))
	E2end(logid)
	if options.fscvol==None: pl.join()
	
def run(cmd):
	print(cmd)
	launch_childprocess(cmd)
	
	
if __name__ == '__main__':
	main()
	
