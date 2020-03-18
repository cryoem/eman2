#!/usr/bin/env python
# Muyuan Chen 2017-03
from builtins import range
from EMAN2 import *
import numpy as np
from EMAN2_utils import *

def main():
	
	usage="""prog <particle stack> --ref <reference> [options]
	Iterative subtomogram refinement.  
	"""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_pos_argument(name="particles",help="Specify particles to use to generate an initial model.", default="", guitype='filebox', browser="EMSetsTable(withmodal=True,multiselect=False)", row=0, col=0,rowspan=1, colspan=3, mode="model")
	parser.add_argument("--reference","--ref",help="""3D reference for iterative alignment/averaging. No reference is used by default. <name> or <name>,#""", default="", guitype='filebox', browser="EMBrowserWidget(withmodal=True,multiselect=False)", row=1, col=0,rowspan=1, colspan=3, mode="model")

	parser.add_header(name="orblock1", help='Just a visual separation', title="Options", row=2, col=1, rowspan=1, colspan=1, mode="model")

	parser.add_argument("--mask", type=str,help="Mask file to be applied to initial model", default="", guitype='filebox', browser="EMBrowserWidget(withmodal=True,multiselect=False)", row=3, col=0,rowspan=1, colspan=3, mode="model")

	parser.add_argument("--niter", type=int,help="Number of iterations", default=5, guitype='intbox',row=4, col=0,rowspan=1, colspan=1, mode="model")
	parser.add_argument("--sym", type=str,help="symmetry", default="c1", guitype='strbox',row=4, col=1,rowspan=1, colspan=1, mode="model")
	
	parser.add_argument("--mass", type=float,help="mass. default -1 will skip by mass normalization", default=-1, guitype='floatbox',row=5, col=0,rowspan=1, colspan=1, mode="model")
	parser.add_argument("--localfilter", action="store_true", default=False ,help="use tophat local", guitype='boolbox',row=5, col=2,rowspan=1, colspan=1, mode="model")

	parser.add_argument("--goldstandard", type=int,help="initial resolution for gold standard refinement", default=-1, guitype='intbox',row=6, col=0,rowspan=1, colspan=1, mode="model")
	parser.add_argument("--goldcontinue", action="store_true", default=False ,help="continue from an existing gold standard refinement", guitype='boolbox',row=6, col=1,rowspan=1, colspan=1, mode="model")

	parser.add_argument("--setsf", type=str,help="structure factor", default=None, guitype='filebox', browser="EMBrowserWidget(withmodal=True,multiselect=False)", row=7, col=0,rowspan=1, colspan=3, mode="model")

	parser.add_argument("--pkeep", type=float,help="fraction of particles to keep", default=0.8, guitype='floatbox',row=8, col=0,rowspan=1, colspan=1, mode="model")

	parser.add_argument("--maxtilt",type=float,help="Explicitly zeroes data beyond specified tilt angle. Assumes tilt axis exactly on Y and zero tilt in X-Y plane. Default 90 (no limit).",default=90.0, guitype='floatbox',row=8, col=2,rowspan=1, colspan=1, mode="model")

	parser.add_argument("--threads", type=int,help="threads", default=12, guitype='intbox',row=9, col=0,rowspan=1, colspan=1, mode="model")

	parser.add_argument("--path", type=str,help="Specify name of refinement folder. Default is spt_XX.", default=None)#, guitype='strbox', row=10, col=0,rowspan=1, colspan=3, mode="model")
	parser.add_argument("--maxang",type=float,help="maximum anglular difference in refine mode.",default=30)
	parser.add_argument("--maxshift",type=float,help="maximum shift in pixel.",default=-1)

	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-2)
	parser.add_argument("--parallel", type=str,help="Thread/mpi parallelism to use", default="")
	parser.add_argument("--refine",action="store_true",help="local refinement from xform.init in header.",default=False)
	parser.add_argument("--randphi",action="store_true",help="randomize phi for refine search",default=False)
	parser.add_argument("--resume",action="store_true",help="resume from previous run",default=False)
	
	#parser.add_argument("--masktight", type=str,help="Mask_tight file", default="")

	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)
	ptcls=args[0]
	ref=options.reference.split(",")[0]
	try: refn=int(options.reference.split(",")[1])
	except: refn=0

	if options.resume:
		try:
			itr=int(ref[-6:-4])
		except:
			print("Can only resume using previous threed_xx.hdf as reference")
			return
		
		options.path=os.path.dirname(ref)
		print("Work in path {}".format(options.path))
		startitr=itr+1
		try:
			fsc=np.loadtxt(os.path.join(options.path, "fsc_masked_{:02d}.txt".format(itr)))
			rs=1./fsc[fsc[:,1]<0.3, 0][0]
		except:
			rs=20
		curres=rs
		if curres>20.:
			curres=20
		print("Start resolution {:.1f}".format(curres))
		
	else:
		curres=20
		startitr=1
		
		
	if options.path==None: options.path = make_path("spt") 
	if options.parallel=="":
		options.parallel="thread:{}".format(options.threads)
	
	#### make a list file if the particles are not in a lst
	if not ptcls.endswith(".lst"):
		ptcllst="{}/input_ptcls.lst".format(options.path)
		run("e2proclst.py {} --create {}".format(ptcls, ptcllst))
		ptcls=ptcllst

	options.input_ptcls=ptcls
	options.input_ref=ref
	options.cmd=' '.join(sys.argv)
	for i in range(10):
		fm="{}/{}_spt_params.json".format(options.path, i)
		if not os.path.isfile(fm):
			js=js_open_dict(fm)
			js.update(vars(options))
			js.close()
			break
		
	
	
	
	ep=EMData(ptcls,0)
	
	if options.goldcontinue==False:
		er=EMData(ref,refn,True)
		if er["apix_x"]==1.0 : print("Warning, A/pix exactly 1.0. You may wish to double-check that this is correct!")
		if abs(1-ep["apix_x"]/er["apix_x"])>0.01 or ep["nx"]!=er["nx"]:
			print("apix mismatch {:.2f} vs {:.2f}".format(ep["apix_x"], er["apix_x"]))
			rs=er["apix_x"]/ep["apix_x"]
			if rs>1.:
				run("e2proc3d.py {} {}/model_input.hdf --clip {} --scale {} --process mask.soft:outer_radius=-1 --first {} --last {}".format(ref, options.path, ep["nx"], rs,refn,refn))
			else:
				run("e2proc3d.py {} {}/model_input.hdf --scale {} --clip {} --process mask.soft:outer_radius=-1 --first {} --last {}".format(ref, options.path, rs, ep["nx"],refn,refn))
		else:	
				run("e2proc3d.py {} {}/model_input.hdf --process mask.soft:outer_radius=-1 --first {} --last {}".format(ref, options.path, refn,refn))
		ref="{}/model_input.hdf".format(options.path)
		
	
	for itr in range(startitr,options.niter+startitr):
		
		#### generate alignment command first
		gd=""
		if options.goldstandard>0 and itr==1:
			curres=options.goldstandard
			gd=" --goldstandard {}".format(options.goldstandard)
		if options.goldcontinue or (options.goldstandard>0 and itr>1):
			gd=" --goldcontinue".format(options.goldstandard)
		
		if options.refine:
			gd+=" --refine --maxang {:.1f}".format(options.maxang)
			if options.randphi:
				gd+=" --randphi"
		#curres=0
		
		if options.maxshift>0:
			gd+=" --maxshift {:.1f}".format(options.maxshift)

		cmd="e2spt_align.py {} {} --parallel {} --path {} --iter {} --sym {} --nsoln 1 {}".format(ptcls, ref,  options.parallel, options.path, itr, options.sym, gd)
		
		#### in case e2spt_align get segfault....
		ret=1
		while ret>0:
			try: os.remove(os.path.join(options.path, "particle_parms_{:02d}.json".format(itr)))
			except:pass

			ret=run(cmd)
		
		js=js_open_dict(os.path.join(options.path, "particle_parms_{:02d}.json".format(itr)))
		score=[]
		for k in list(js.keys()):
			score.append(float(js[k]["score"]))
		
		s=""
		if options.pkeep<1:
			simthr=np.sort(score)[int(len(score)*options.pkeep)]
			s+=" --simthr {:f}".format(simthr)
		if options.maxtilt<90.:
			s+=" --maxtilt {:.1f}".format(options.maxtilt)
			
		
		#run("e2spt_average.py --threads {} --path {} --sym {} --skippostp {}".format(options.threads, options.path, options.sym, s))
		run("e2spt_average.py --parallel {} --path {} --sym {} --skippostp {}".format(options.parallel, options.path, options.sym, s))
		
		even=os.path.join(options.path, "threed_{:02d}_even.hdf".format(itr))
		odd=os.path.join(options.path, "threed_{:02d}_odd.hdf".format(itr))
		combine=os.path.join(options.path, "threed_{:02d}.hdf".format(itr))
		
		if options.setsf==None:
			#### do a simple amplitute correction when no sf provided
			data=EMData(combine)
			bxsz=data["nx"]
			apix=data["apix_x"]
			dataf = data.do_fft()
			curve = dataf.calc_radial_dist((data["ny"]//2), 0, 1.0, False)
			curve=np.array([i/dataf["nx"]*dataf["ny"]*dataf["nz"] for i in curve])
			s=np.arange(len(curve))*1./(apix*bxsz)
			sf=XYData()
			sf.set_xy_list(s.tolist(), curve.tolist())
			
			
		for f in [even, odd]:
			e=EMData(f)
			e.del_attr("xform.align3d")
			e.write_image(f.replace("threed_{:02d}_".format(itr), "threed_raw_"))
			if options.setsf==None:
				e.process_inplace("filter.setstrucfac",{"apix":e["apix_x"],"strucfac":sf})
			e.write_image(f)
			
		
		msk=options.mask
		if len(msk)>0:
			if os.path.isfile(msk):
				msk=" --automask3d mask.fromfile:filename={}".format(msk)
			else:
				msk=" --automask3d {}".format(msk)

		#if len(options.masktight)>0:
			#if os.path.isfile(options.masktight):
				#msk+=" --automask3dtight mask.fromfile:filename={}".format(options.masktight)
			#else:
				#msk+=" --automask3dtight {}".format(options.masktight)

		s=""
		if options.goldstandard>0:
			s+=" --align"
		
		if options.setsf:
			s+=" --setsf {}".format(options.setsf)
			
		if options.localfilter:
			s+=" --tophat local "
			
		ppcmd="e2refine_postprocess.py --even {} --odd {} --output {} --iter {:d} --mass {} --restarget {} --threads {} --sym {} {} {} ".format(even, odd, os.path.join(options.path, "threed_{:02d}.hdf".format(itr)), itr, options.mass, curres, options.threads, options.sym, msk, s)
		run(ppcmd)
		
		#if options.localnorm:
			#for f in [even, odd]:
				#run("e2proc3d.py {} {} --process normalize --process normalize.local:threshold=1:radius=16".format(f,f))
			
		ref=os.path.join(options.path, "threed_{:02d}.hdf".format(itr))
		fsc=np.loadtxt(os.path.join(options.path, "fsc_masked_{:02d}.txt".format(itr)))
		
		fi=fsc[:,1]<0.3
		if np.sum(fi)==0:
			print("something wrong with the FSC curve. Cannot estimate resolution. Please check.")
		else:
			rs=1./fsc[fi, 0][0]
			print("Resolution (FSC<0.3) is ~{:.1f} A".format(rs))
		curres=rs
		if curres>40.:
			curres=40
		

	E2end(logid)
	
def run(cmd):
	print(cmd)
	ret=launch_childprocess(cmd)
	return ret
	
	
if __name__ == '__main__':
	main()

