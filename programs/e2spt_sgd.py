#!/usr/bin/env python
# Muyuan Chen 2017-04
from past.utils import old_div
from future import standard_library
standard_library.install_aliases()
from builtins import range
from EMAN2 import *
import sys
import numpy as np
import threading
import queue
#from e2spt_align import alifn

def alifn(jsd,fsp,i,a,options):
	t=time.time()
	b=EMData(fsp,int(i))#.do_fft()
	b.process_inplace("filter.lowpass.gauss", {"cutoff_freq":options.filterto})
	if options.shrink>1:
		b.process_inplace("math.fft.resample",{"n":options.shrink})
	b=b.do_fft()
	b.process_inplace("xform.phaseorigin.tocorner")
	
	aligndic={"maxshift":10,"verbose":options.verbose,"sym":"c1","sigmathis":0.1,"sigmato":1.0}
	if options.refine and b.has_attr("xform.align3d"):
		ntry=8
		initxf=b["xform.align3d"]
		xfs=[]
		
		for ii in range(len(xfs), ntry):
			ixf=initxf.get_params("eman")
			ixf["phi"]=np.random.rand()*360.
			ixf["alt"]=ixf["alt"]+180*(np.random.rand()>.5)
			ixf=Transform(ixf)
			
			v=np.random.rand(3)-0.5
			nrm=np.linalg.norm(v)
			if nrm>0:
				v=v/nrm
			else:
				v=(0,0,1)
			xf=Transform({"type":"spin", "n1":v[0], "n2":v[1], "n3":v[2],
					"omega":30.0*np.random.randn()/3.0})
			xfs.append(xf*ixf)
		
			
		aligndic["initxform"]=xfs
		aligndic["maxang"]=30.0
		aligndic["randphi"]=True
		aligndic["rand180"]=True

	# we align backwards due to symmetry
	c=a.xform_align_nbest("rotate_translate_3d_tree",b,aligndic,1)
	for cc in c : cc["xform.align3d"]=cc["xform.align3d"].inverse()

	jsd.put((fsp,i,c[0]))
	if options.verbose>1 : print("{}\t{}\t{}\t{}".format(fsp,i,time.time()-t,c[0]["score"]))


def main():

	usage="""prog <particle stack>
	Generate initial model from subtomogram particles using stochastic gradient descent.
	"""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_pos_argument(name="particles",help="Specify particles to use to generate an initial model.", default="", guitype='filebox', browser="EMSetsTable(withmodal=True,multiselect=False)", row=0, col=0,rowspan=1, colspan=2, mode="model")

	parser.add_header(name="orblock1", help='Just a visual separation', title="Options", row=1, col=1, rowspan=1, colspan=1, mode="model")

	parser.add_argument('--reference','--ref', type=str, default="", help="""3D reference for initial model generation. No reference is used by default.""", guitype='filebox', browser="EMBrowserWidget(withmodal=True,multiselect=False)", row=2, col=0,rowspan=1, colspan=2, mode="model")

	parser.add_argument("--mask", type=str,help="Mask file to be applied to initial model", default=None, guitype='filebox', browser="EMBrowserWidget(withmodal=True,multiselect=False)", row=3, col=0,rowspan=1, colspan=2, mode="model")

	parser.add_argument("--sym", type=str,help="symmetry", default="c1", guitype='strbox',row=4, col=0,rowspan=1, colspan=1, mode="model")
	#parser.add_argument("--gaussz", type=float,help="Extra gauss filter at z direction. seems useless...", default=-1)

	parser.add_argument("--filterto", type=float,help="Fiter map to frequency after each iteration. Default is 0.02", default=.02, guitype='floatbox',row=6, col=0,rowspan=1, colspan=1, mode="model")

	parser.add_argument("--fourier", action="store_true", default=False ,help="gradient descent in fourier space", guitype='boolbox',row=6, col=1,rowspan=1, colspan=1, mode="model[True]")

	parser.add_argument("--batchsize", type=int,help="SGD batch size", default=12,guitype='intbox',row=9, col=0,rowspan=1, colspan=1, mode="model")
	parser.add_argument("--learnrate", type=float,help="Learning rate. Default is 0.1", default=.1,guitype='floatbox',row=9, col=1,rowspan=1, colspan=1, mode="model")

	parser.add_argument("--niter", type=int,help="Number of iterations", default=5, guitype='intbox',row=10, col=0,rowspan=1, colspan=1, mode="model")
	parser.add_argument("--nbatch", type=int,help="Number of batches per iteration", default=10,guitype='intbox',row=10, col=1,rowspan=1, colspan=1, mode="model")

	parser.add_argument("--applysym", action="store_true", default=False ,help="apply symmetry", guitype='boolbox',row=11, col=0,rowspan=1, colspan=1, mode="model")
	
	parser.add_argument("--writemovie", action="store_true", default=False ,help="write all temporary files as a stack")
	parser.add_argument("--shrink", type=int,help="Shrink factor for particles", default=1,guitype='intbox',row=11, col=1,rowspan=1, colspan=1, mode="model")

	parser.add_argument("--path", type=str,help="path of output", default=None)
	
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")

	parser.add_argument("--ppid", type=int,help="ppid", default=-2)
	
	parser.add_argument("--refine", action="store_true", default=False ,help="start from xform.align3d in header", guitype='boolbox',row=12, col=0,rowspan=1, colspan=1, mode="model")


	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)

	if options.path==None:
		for i in range(100):
			pname="sptsgd_{:02d}".format(i)
			if not os.path.isdir(pname):
				os.mkdir(pname)
				options.path=pname
				break
		else:
			print("something is wrong...")
			exit()
	else:
		if not os.path.isdir(options.path):
			os.mkdir(options.path)
		else:
			print("Overwritting {}...".format(options.path))
			os.system("rm {}/*".format(options.path))

	path=options.path
	print("Writing in {}..".format(path))
	fname=args[0]

	ref=make_ref(fname, options)
	
	options.input_ptcls=fname
	options.cmd=' '.join(sys.argv)
	js=js_open_dict("{}/0_spt_params.json".format(options.path))
	js.update(vars(options))
	js.close()

	num=EMUtil.get_image_count(fname)
	batchsize=options.batchsize
	learnrate=options.learnrate

	options.threads=options.batchsize
	#nbatch=len(idxs[0])/batchsize

	jspast=None
	filterto=options.filterto
	for itr in range(1, options.niter+1):

		nbatch=options.nbatch
		if itr==1: nbatch*=2
		print('#'*nbatch)

		jspm=js_open_dict("{}/particle_parms_{:02d}.json".format(options.path,itr))
		if itr>1:
			jspm.update(jspast)

		print("Iteration {}:".format(itr))

		ref0=ref.copy()

		#tmpout=os.path.join(path,"tmpout_{:02d}_{}.hdf".format(itr, eo))
		#ref.write_image(tmpout,-1)
		cc=[]
		for ib in range(nbatch):

			jsd=queue.Queue(0)
			idx=np.arange(num)
			np.random.shuffle(idx)
			thrds=[threading.Thread(target=alifn,args=(jsd,fname,i,ref0,options)) for i in idx[:batchsize]]
			#thrds=[threading.Thread(target=alifn,args=(jsd,fname,i,ref,options)) for i in idx[ib*batchsize:(ib+1)*batchsize]]
			for t in thrds:
				t.start()
			angs={}
			while threading.active_count()>1 or not jsd.empty():
				time.sleep(1)
				while not jsd.empty():
					fsp,n,d=jsd.get()
					angs[(fsp,n)]=d
			avgr=Averagers.get("mean.tomo")
			#print(angs)
			for ks in list(angs.keys()):
				d=angs[ks]
				jspm[ks]=d
				p=EMData(str(ks[0]), int(ks[1]))
				if options.shrink>1:
					p.process_inplace("math.fft.resample",{"n":options.shrink})
				p.transform(d["xform.align3d"])
				avgr.add_image(p)
			avg=avgr.finish()
			avg.process_inplace('filter.lowpass.gauss', {"cutoff_freq":filterto})
			#if options.gaussz>0:
				#avg.process_inplace('filter.lowpass.gauss', {"cutoff_freq":options.gaussz})
			avg.process_inplace('normalize.edgemean')
			#if options.applysym:
				#avg.process_inplace("xform.applysym",{"sym":options.sym,"averager":"mean.tomo"})
			
			if options.fourier:
				avgft=avg.do_fft()
				refft=ref.do_fft()
				avgft.process_inplace("mask.wedgefill",{"fillsource":refft, "thresh_sigma":1})
				avg=avgft.do_ift()

				#dmap=avgft-refft
				#refft=refft+learnrate*dmap
				#refnew=refft.do_ift()
				#refnew.process_inplace('normalize.edgemean')
				#dmap=refnew-ref
				#ref=refnew.copy()
			#else:

			dmap=avg-ref
			ref=ref+learnrate*dmap
			ref.process_inplace('normalize.edgemean')

			ddm=dmap*dmap
			cc.append(ddm["mean_nonzero"])
			
			ref0=ref.copy()
			if options.mask:
				if os.path.isfile(options.mask):
				
					ref0.process_inplace("mask.fromfile", {"filename": options.mask})
				else:
					pdict=parsemodopt(options.mask)
					ref0.process_inplace(pdict[0], pdict[1])
					
				ref.write_image(os.path.join(path,"output_nomask.hdf"))

			#ref.write_image(tmpout,-1)
			ref0.write_image(os.path.join(path,"output.hdf"))
			
			
			if options.writemovie:
				#ref0.write_image(os.path.join(path,"output_all.hdf"), -1)
				ref.write_image(os.path.join(path,"allref.hdf"), -1)
				avg.write_image(os.path.join(path,"allavg.hdf"), -1)
				
			sys.stdout.write('#')
			sys.stdout.flush()

		print("  mean gradient {:.3f}".format(np.mean(cc)))

		jspast=jspm.data.copy()
		jspm=None
		#### if symmetry exist, first align to symmetry axis
		if options.sym!="c1" and len(options.reference)==0:
			ref=sym_search(ref, options)
		if options.applysym:
			ref.process_inplace("xform.applysym", {"averager":"mean.tomo", "sym":options.sym})
			
		if len(options.reference)==0:
			ref.process_inplace("xform.centerofmass")
		ref.write_image("{}/threed_{:02d}.hdf".format(path, itr), 0)

	#ref.write_image(os.path.join(path,"output.hdf"))
	print("Done")
	E2end(logid)


def sym_search(e, options):
	print("Align to symmetry axis...")
	ntry=20
	sym=options.sym
	s=parsesym(sym)
	if options.refine:
		oris=[Transform()]
	else:
		oris=s.gen_orientations("rand",{"n":ntry, "phitoo":True})
		
	jsd=queue.Queue(0)
	thrds=[threading.Thread(target=sym_ali,args=([e, o, sym, jsd])) for o in oris]
	for t in thrds: t.start()
	for t in thrds: t.join()

	alis=[]
	while not jsd.empty():
		alis.append(jsd.get())
	#print alis
	scr=[a["score"] for a in alis]
	im=np.argmin(scr)
	a=alis[im]
	a.process_inplace("normalize.edgemean")
	#a.process_inplace("xform.centerofmass")
	
	#outname=fname[:-4]+"_sym.hdf"
	#a.write_image(outname)
	return a

def sym_ali(e,o,sym,jsd):
	s=old_div(e["nx"],4)
	a=e.align("symalignquat", e, {"sym":sym, "xform.align3d":o,"maxshift":s}, "ccc.tomo.thresh")
	jsd.put(a)


def make_ref(fname, options):

	num=EMUtil.get_image_count(fname)
	rfile="{}/ref.hdf".format(options.path)
	if len(options.reference)==0:
		print("Making random references...")
		
		tt=parsesym("c1")
		xfs=tt.gen_orientations("rand",{"n":options.batchsize})
		idx=np.arange(num)
		np.random.shuffle(idx)
		avgr=Averagers.get("mean.tomo")
		for i in range(options.batchsize):
			p=EMData(fname, int(idx[i]))
			p.process_inplace('normalize.edgemean')
			p.process_inplace('xform.centerofmass')
			p.transform(xfs[i])
			avgr.add_image(p)
		ref=avgr.finish()
		ref.process_inplace('filter.lowpass.gauss', {"cutoff_freq":.01})
		ref.process_inplace('filter.lowpass.randomphase', {"cutoff_freq":.01})
		if options.shrink>1:
			ref.process_inplace("math.fft.resample",{"n":options.shrink})
		#ref.process_inplace("xform.applysym",{"sym":options.sym})
		ref.process_inplace('normalize.edgemean')
		ref.process_inplace('xform.centerofmass')
		ref.process_inplace("mask.soft",{"outer_radius":-4})
		ref.write_image(rfile)
	else:
		er=EMData(options.reference)
		ep=EMData(fname,0)
		pp=" --process filter.lowpass.gauss:cutoff_freq={:.3f} --process filter.lowpass.randomphase:cutoff_freq={:.3f} --process normalize.edgemean".format(options.filterto, options.filterto)
		if options.shrink>1:
			pp+=" --process math.fft.resample:n={}".format(options.shrink)


		if ep["apix_x"]!=er["apix_x"] or ep["nx"]!=er["nx"]:
			print("apix mismatch {:.2f} vs {:.2f}".format(ep["apix_x"], er["apix_x"]))
			rs=er["apix_x"]/ep["apix_x"]

			if rs>1.:
				run("e2proc3d.py {} {}/ref.hdf --clip {} --scale {} {}".format(options.reference, options.path, ep["nx"], rs, pp))
			else:
				run("e2proc3d.py {} {}/ref.hdf --scale {} --clip {} {}".format(options.reference, options.path,  rs, ep["nx"], pp))
		else:
			run("e2proc3d.py {} {}/ref.hdf {}".format(options.reference, options.path, pp))
		
		ref=EMData(rfile)

	return ref
def run(cmd):
	print(cmd)
	launch_childprocess(cmd)


if __name__ == '__main__':
	main()

