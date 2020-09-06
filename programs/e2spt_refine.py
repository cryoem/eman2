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
	parser.add_argument("--maskalign", type=str,help="Mask file applied to 3D alignment reference in each iteration. Not applied to the average, which will follow normal masking routine.", default="", guitype='filebox', browser="EMBrowserWidget(withmodal=True,multiselect=False)", row=4, col=0,rowspan=1, colspan=3, mode="model")

	parser.add_argument("--niter", type=int,help="Number of iterations", default=5, guitype='intbox',row=5, col=0,rowspan=1, colspan=1, mode="model")
	parser.add_argument("--sym", type=str,help="symmetry", default="c1", guitype='strbox',row=5, col=1,rowspan=1, colspan=1, mode="model")
	
	parser.add_argument("--mass", type=float,help="mass. default -1 will skip by mass normalization", default=-1, guitype='floatbox',row=5, col=2,rowspan=1, colspan=1, mode="model")
	parser.add_argument("--localfilter", action="store_true", default=False ,help="use tophat local", guitype='boolbox',row=6, col=2,rowspan=1, colspan=1, mode="model")

	parser.add_argument("--goldstandard", type=int,help="initial resolution for gold standard refinement", default=-1, guitype='intbox',row=6, col=0,rowspan=1, colspan=1, mode="model")
	parser.add_argument("--goldcontinue", action="store_true", default=False ,help="continue from an existing gold standard refinement", guitype='boolbox',row=6, col=1,rowspan=1, colspan=1, mode="model")

	parser.add_argument("--setsf", type=str,help="structure factor", default=None, guitype='filebox', browser="EMBrowserWidget(withmodal=True,multiselect=False)", row=7, col=0,rowspan=1, colspan=3, mode="model")

	parser.add_argument("--pkeep", type=float,help="fraction of particles to keep", default=0.8, guitype='floatbox',row=8, col=0,rowspan=1, colspan=1, mode="model")

	parser.add_argument("--maxtilt",type=float,help="Explicitly zeroes data beyond specified tilt angle. Assumes tilt axis exactly on Y and zero tilt in X-Y plane. Default 90 (no limit).",default=90.0, guitype='floatbox',row=8, col=2,rowspan=1, colspan=1, mode="model")


	parser.add_argument("--path", type=str,help="Specify name of refinement folder. Default is spt_XX.", default=None)#, guitype='strbox', row=10, col=0,rowspan=1, colspan=3, mode="model")
	parser.add_argument("--maxang",type=float,help="maximum anglular difference in refine mode.",default=30)
	parser.add_argument("--maxshift",type=float,help="maximum shift in pixel.",default=-1)

	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-2)
	parser.add_argument("--threads", type=int,help="threads", default=12, guitype='intbox',row=9, col=0,rowspan=1, colspan=1, mode="model")
	parser.add_argument("--parallel", type=str,help="Thread/mpi parallelism to use", default="")
	parser.add_argument("--transonly",action="store_true",help="translational alignment only",default=False)
	parser.add_argument("--refine",action="store_true",help="local refinement from xform in header.",default=False)
	parser.add_argument("--randphi",action="store_true",help="randomize phi for refine search",default=False)
	parser.add_argument("--rand180",action="store_true",help="include 180 degree rotation for refine search",default=False)
	parser.add_argument("--resume",action="store_true",help="resume from previous run",default=False)
	parser.add_argument("--scipy",action="store_true",help="test scipy refinement",default=False)
	parser.add_argument("--breaksym",action="store_true",help="break symmetry",default=False)
	parser.add_argument("--breaksymsym", type=str,help="Specify a different symmetry for breaksym.", default=None)
	parser.add_argument("--symalimask",type=str,default=None,help="This will translationally realign each asymmetric unit to the previous map masked by the specified mask. While this invokes symalimasked in e2spt_average, this isn't the same, it is a mask, not a masked reference. ")
	
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
	
	msk=options.mask
	if len(msk)>0:
		if os.path.isfile(msk):
			msk=" --automask3d mask.fromfile:filename={}".format(msk)
		else:
			msk=" --automask3d {}".format(msk)

	#### make a list file if the particles are not in a lst
	if ptcls.endswith(".json"):
		jstmp=js_open_dict(ptcls)
		ky=jstmp.keys()[0]
		pt,ii=eval(ky)
		ep=EMData(pt, ii)
		
	elif ptcls.endswith(".lst"):
		ep=EMData(ptcls,0)
		
	else:
		ptcllst="{}/input_ptcls.lst".format(options.path)
		run("e2proclst.py {} --create {}".format(ptcls, ptcllst))
		ptcls=ptcllst
		ep=EMData(ptcls,0)

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

		# the alignment ref may be masked using a different file, or just copied
		ar=EMData(ref,0)
		if (len(options.maskalign)>0):
			m=EMData(options.maskalign,0)
			ar.mult(m)
		ar.write_image(f"{options.path}/alignref.hdf",0)
		
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
			if options.rand180:
				gd+=" --rand180"
			if itr>startitr:
				ptcls=os.path.join(options.path, "particle_parms_{:02d}.json".format(itr-1))

		if options.transonly: gd+=" --transonly"

		
		if options.breaksym:
			gd+=" --breaksym"
		if options.breaksymsym:
			gd+=" --breaksymsym {}".format(options.breaksymsym)
		
		if options.maxshift>0:
			gd+=" --maxshift {:.1f}".format(options.maxshift)
		if options.scipy:
			gd+=" --scipy"

		cmd="e2spt_align.py {} {}/alignref.hdf --parallel {} --path {} --iter {} --sym {} --maxres {} {}".format(ptcls, options.path,  options.parallel, options.path, itr, options.sym, curres, gd)
		
		ret=run(cmd)
		
		
		s=""
		if options.maxtilt<90.:
			s+=" --maxtilt {:.1f}".format(options.maxtilt)

		# we apply the symmetric subunit mask provided to the current reference and send it to e2spt_average to do a final translational alignment
		if options.symalimask!=None: 
			cmd=f"e2proc3d.py {ref} {options.path}/ref_mono.hdf --multfile {options.symalimask}"
			run(cmd)
			s+=f" --symalimasked={options.path}/ref_mono.hdf"
			
		
		run("e2spt_average.py --threads {} --path {} --sym {} --skippostp {}".format(options.threads, options.path, options.sym, s))
		#run("e2spt_average.py --parallel {} --path {} --sym {} --keep {:.3f} --iter {} --skippostp {}".format(options.parallel, options.path, options.sym, options.pkeep, itr, s))
		
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
			
		
		s=""
		if options.goldstandard>0:
			s+=" --align"
		
		if options.setsf:
			s+=" --setsf {}".format(options.setsf)
			
		if options.localfilter:
			s+=" --tophat local "
			
		# if we are doing local symmetry refinement or breaking the symmetry
		# it's a bit counterproductive if we then apply symmetry here (as was happening before 8/22/20)
		syms=f"--sym {options.sym}"
		if options.symalimask!=None or options.breaksym: syms=""
		run(f"e2refine_postprocess.py --even {even} --odd {odd} --output {options.path}threed_{itr:02d}.hdf --iter {itr:d} --tomo --mass {options.mass} --threads {options.threads} {syms} {msk} {s}")

		try: symn=int(options.sym[1:])
		except: symn=0
		if options.symalimask!=None and not options.breaksym and symn>0:
			os.rename(f"{options.path}threed_{itr:02d}.hdf",f"{options.path}threed_{itr:02d}_nosym.hdf")
			phir=360.0/(symn*2.0)
			if   options.sym[0].lower()=="c":
				run(f"e2proc3d.py {options.path}threed_{itr:02d}_nosym.hdf {options.path}threed_{itr:02d}.hdf --process mask.cylinder:phicen=0:phirange={phir-5.0}:phitriangle=1:phitrirange=10.0 --sym {options.sym}")
			elif options.sym[0].lower()=="d":
				run(f"e2proc3d.py {options.path}threed_{itr:02d}_nosym.hdf {options.path}threed_{itr:02d}.hdf --process mask.cylinder:phicen=0:phirange={phir-5.0}:phitriangle=1:phitrirange=10.0:zmin={data['nz']/2}:zmax={data['nz']} --sym {options.sym}")

		ref=os.path.join(options.path, "threed_{:02d}.hdf".format(itr))
		fsc=np.loadtxt(os.path.join(options.path, "fsc_masked_{:02d}.txt".format(itr)))
		
		fi=fsc[:,1]<0.2
		if np.sum(fi)==0:
			print("something wrong with the FSC curve. Cannot estimate resolution. Please check.")
		else:
			rs=1./fsc[fi, 0][0]
			print("Resolution (FSC<0.2) is ~{:.1f} A".format(rs))
		curres=rs*.8
		if curres>40.:
			curres=40
		

	E2end(logid)
	
def run(cmd):
	print(cmd)
	ret=launch_childprocess(cmd)
	return ret
	
	
if __name__ == '__main__':
	main()

