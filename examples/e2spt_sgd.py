#!/usr/bin/env python
# Muyuan Chen 2017-04
from __future__ import print_function
from EMAN2 import *
import numpy as np
import threading
import Queue
#from e2spt_align import alifn

def alifn(jsd,fsp,i,a,options):
	t=time.time()
	b=EMData(fsp,i).do_fft()
	b.process_inplace("xform.phaseorigin.tocorner")

	# we align backwards due to symmetry
	c=a.xform_align_nbest("rotate_translate_3d_tree",b,{"verbose":0,"sym":options.sym,"sigmathis":0.1,"sigmato":1.0},1)
	for cc in c : cc["xform.align3d"]=cc["xform.align3d"].inverse()

	jsd.put((fsp,i,c[0]))
	if options.verbose>1 : print("{}\t{}\t{}\t{}".format(fsp,i,time.time()-t,c[0]["score"]))


def main():
	
	usage=" "
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--path", type=str,help="path", default=None)
	parser.add_argument("--ref", type=str,help="ref", default=None)
	parser.add_argument("--sym", type=str,help="symmetry", default="c1")
	parser.add_argument("--batchsize", type=int,help="batch size", default=12)
	parser.add_argument("--niter", type=int,help="iterations", default=50)
	parser.add_argument("--learnrate", type=float,help="learnrate", default=.1)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	
	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)
	
	from EMAN2_utils import makepath
	options = makepath(options,'sptsgd')
	
	#if options.path==None:
	#	for i in range(100):
	#		pname="sptsgd_{:02d}".format(i)
	#		if not os.path.isdir(pname):
	#			os.mkdir(pname)
	#			options.path=pname
	#			break
	#	else:
	#		print("something is wrong...")
	#		exit()
			
	path=options.path
	print("Writing in {}..".format(path))
	fname=args[0]
	num=EMUtil.get_image_count(fname)
	batchsize=options.batchsize
	
	if not options.ref:
		
		tt=parsesym("c1")
		xfs=tt.gen_orientations("rand",{"n":batchsize})
		idx=np.arange(num)
		np.random.shuffle(idx)
		avgr=Averagers.get("mean.tomo")
		for i in range(batchsize):
			p=EMData(fname, idx[i])
			p.transform(xfs[i])
			avgr.add_image(p)
		ref=avgr.finish()
		ref.process_inplace('filter.lowpass.gauss', {"cutoff_freq":.01})
		ref.process_inplace('filter.lowpass.randomphase', {"cutoff_freq":.01})
		ref.process_inplace("xform.applysym",{"sym":options.sym})
		ref.write_image(os.path.join(path,"ref.hdf"))
	else:
		ref=EMData(options.ref)
		
	learnrate=options.learnrate
	lrmult=.98
	tmpout=os.path.join(path,"tmpout.hdf")
	try: os.remove(tmpout)
	except: pass
	ref.write_image(tmpout,-1)
	print("iteration, learning rate, mean gradient")
	for it in range(options.niter):
		idx=np.arange(num)
		np.random.shuffle(idx)
		nbatch=num/batchsize
		cc=[]
		for ib in range(nbatch):
			jsd=Queue.Queue(0)
			thrds=[threading.Thread(target=alifn,args=(jsd,fname,i,ref,options)) for i in idx[ib*batchsize:(ib+1)*batchsize]]
			for t in thrds:
				t.start()
			angs={}
			while threading.active_count()>1:
				time.sleep(1)
				while not jsd.empty():
					fsp,n,d=jsd.get()
					angs[(fsp,n)]=d
			avgr=Averagers.get("mean.tomo")
			#print angs
			for ks in angs.keys():
				d=angs[ks]
				p=EMData(ks[0], ks[1])
				p.transform(d["xform.align3d"])
				avgr.add_image(p)
			avg=avgr.finish()
			avg.process_inplace("xform.applysym",{"sym":options.sym})
			dmap=avg-ref
			ddm=dmap*dmap
			#print "{:d}\t{:.3f}\t{:.3f}".format(it, ddm["mean_nonzero"], np.mean(scr))
			cc.append(ddm["mean_nonzero"])
			ref=ref+learnrate*dmap
		ref.process_inplace("xform.centerofmass")
		ref.write_image(tmpout,-1)
		
		
		#ref.write_image(tmpout,-1)
		
		print("\t{:d}, {:.3f}, {:.5f}".format(it, learnrate, np.mean(cc)))
		learnrate*=lrmult	
	
	ref.write_image(os.path.join(path,"output.hdf"))
	print("Done")
	E2end(logid)
	
def run(cmd):
	print(cmd)
	launch_childprocess(cmd)
	
	
if __name__ == '__main__':
	main()
	