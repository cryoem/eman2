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
	parser.add_pos_argument(name="particles",help="Specify particles to use for refinement. May be a sets/*.lst file or a .json file from a previous run.", default="", guitype='filebox', browser="EMSetsTable(withmodal=True,multiselect=False)", row=0, col=0,rowspan=1, colspan=3, mode="model")
	parser.add_argument("--reference","--ref",help="""3D reference for iterative alignment/averaging. No reference is used by default. <name> or <name>,#""", default="", guitype='filebox', browser="EMBrowserWidget(withmodal=True,multiselect=False)", row=1, col=0,rowspan=1, colspan=3, mode="model")

	parser.add_header(name="orblock1", help='Just a visual separation', title="Options", row=2, col=1, rowspan=1, colspan=1, mode="model")

	parser.add_argument("--mask", type=str,help="Mask file to be applied to initial model", default="", guitype='filebox', browser="EMBrowserWidget(withmodal=True,multiselect=False)", row=3, col=0,rowspan=1, colspan=3, mode="model")
	parser.add_argument("--maskalign", type=str,help="Mask file applied to 3D alignment reference in each iteration. Not applied to the average, which will follow normal masking routine.", default=None, guitype='filebox', browser="EMBrowserWidget(withmodal=True,multiselect=False)", row=4, col=0,rowspan=1, colspan=3, mode="model")

	parser.add_argument("--niter", type=int,help="Number of iterations", default=5, guitype='intbox',row=5, col=0,rowspan=1, colspan=1, mode="model")
	parser.add_argument("--sym", type=str,help="symmetry", default="c1", guitype='strbox',row=5, col=1,rowspan=1, colspan=1, mode="model")
	parser.add_argument("--maxres",type=float,help="Maximum resolution (the smaller number) to consider in alignment (in A, not 1/A)",default=0)
	parser.add_argument("--minres",type=float,help="Minimum resolution (the larger number) to consider in alignment (in A, not 1/A, default=200)",default=200)
	
	parser.add_argument("--mass", type=float,help="mass. default -1 will skip by mass normalization", default=-1, guitype='floatbox',row=5, col=2,rowspan=1, colspan=1, mode="model")
	parser.add_argument("--localfilter", action="store_true", default=False ,help="Deprecated. Please use --tophat")
	parser.add_argument("--tophat", type=str, default=None,help = "'local', 'localwiener' or 'global'. Instead of imposing a uniform Wiener filter, use a tophat filter (global similar to Relion). local is a local tophat filter, localwiener is a localized Wiener filter", guitype='strbox', row=6, col=2,rowspan=1, colspan=1, mode="model['local']")

	parser.add_argument("--goldstandard", type=int,help="Phase randomization resolution for gold standard refinement in A. Not equivalent to restarget in e2refine_easy.", default=-1, guitype='intbox',row=6, col=0,rowspan=1, colspan=1, mode="model")
	parser.add_argument("--goldcontinue", action="store_true", default=False ,help="continue from an existing gold standard refinement", guitype='boolbox',row=6, col=1,rowspan=1, colspan=1, mode="model")
	parser.add_argument("--restarget", default=0, type=float,help="The resolution you reasonably expect to achieve in the current refinement run in A.")

	parser.add_argument("--setsf", type=str,help="structure factor", default=None, guitype='filebox', browser="EMBrowserWidget(withmodal=True,multiselect=False)", row=7, col=0,rowspan=1, colspan=3, mode="model")

	parser.add_argument("--pkeep", type=float,help="fraction of particles to keep", default=0.8, guitype='floatbox',row=8, col=0,rowspan=1, colspan=1, mode="model")

	parser.add_argument("--maxtilt",type=float,help="Explicitly zeroes data beyond specified tilt angle. Assumes tilt axis exactly on Y and zero tilt in X-Y plane. Default 90 (no limit).",default=90.0, guitype='floatbox',row=8, col=2,rowspan=1, colspan=1, mode="model")


	parser.add_argument("--path", type=str,help="Specify name of refinement folder. Default is spt_XX.", default=None)#, guitype='strbox', row=10, col=0,rowspan=1, colspan=3, mode="model")
	parser.add_argument("--smooth", type=float,help="smoothing factor for subtlt.", default=40)
	parser.add_argument("--maxang",type=float,help="maximum anglular difference in refine mode.",default=-1)
	parser.add_argument("--maxshift",type=float,help="maximum shift in pixel.",default=-1)

	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-2)
	parser.add_argument("--threads", type=int,help="threads", default=12, guitype='intbox',row=9, col=0,rowspan=1, colspan=1, mode="model")
	parser.add_argument("--parallel", type=str,help="Thread/mpi parallelism to use", default="")
	parser.add_argument("--transonly",action="store_true",help="translational alignment only",default=False)
	parser.add_argument("--refine",action="store_true",help="local refinement from xform in header.",default=False)
	parser.add_argument("--refinetry",type=int,help="Number of local perturbations to initialize local alignment with (default=8)",default=8)
	parser.add_argument("--realign",action="store_true",help="realigns the average to the initial reference to prevent drift in C1 refinements",default=False)
	parser.add_argument("--randphi",action="store_true",help="randomize phi for refine search",default=False)
	parser.add_argument("--rand180",action="store_true",help="include 180 degree rotation for refine search",default=False)
	parser.add_argument("--test180",action="store_true",help="Test for improved alignment with 180 degree rotations even during refine alignment",default=False)
	parser.add_argument("--resume",action="store_true",help="resume from previous run",default=False)
	parser.add_argument("--scipy",action="store_true",help="test scipy refinement",default=False)
	parser.add_argument("--breaksym",action="store_true",help="break symmetry",default=False)
	parser.add_argument("--breaksymsym", type=str,help="Specify a different symmetry for breaksym.", default=None)
	parser.add_argument("--symalimask",type=str,default=None,help="This will translationally realign each asymmetric unit to the previous map masked by the specified mask. While this invokes symalimasked in e2spt_average, this isn't the same, it is a mask, not a masked reference. ")
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higher number means higher level of verboseness")
	
	#parser.add_argument("--masktight", type=str,help="Mask_tight file", default="")

	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)
	ptcls=args[0]
	ref=options.reference.split(",")[0]
	try: refn=int(options.reference.split(",")[1])
	except: refn=0
	
	if (options.maxang>0 and not options.refine) or (options.maxang<=0 and options.refine):
		print("Error: --maxang and --refine must be specified together")
		sys.exit(1)

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
			rs=60
		curres=rs
		if curres>75.:
			curres=75
		print("Start resolution {:.1f}".format(curres))
		
	else:
		if options.maxres!=0 : 
			curres=options.maxres
		elif options.goldstandard>0:
			curres=options.goldstandard
		else: 
			curres=60
		startitr=1
		
	if options.localfilter:
		if options.tophat==None: options.tophat="local"
		
	if options.path==None: options.path = num_path_new("spt") 
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

	boxsize=ep["ny"]
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
		
	#if options.maskalign!=None:
		#maskalign=EMData(options.maskalign,0)
	
	for itr in range(startitr,options.niter+startitr):

		# the alignment ref may be masked using a different file, or just copied
		if options.goldstandard>0 or options.goldcontinue :
			if itr==1: 
				if options.verbose>1: print(f"{ref} -> alignref_even.hdf & odd with random phase past {options.goldstandard} A")
				ar=EMData(ref,0)
				ar.process_inplace("filter.lowpass.randomphase",{"cutoff_freq":1.0/options.goldstandard})
				ar.write_compressed(f"{options.path}/alignref_even.hdf",0,12,erase=True)
				
				# no harm in repeating
				ar.process_inplace("filter.lowpass.randomphase",{"cutoff_freq":1.0/options.goldstandard})
				ar.write_compressed(f"{options.path}/alignref_odd.hdf",0,12,erase=True)
			else:
				refe=ref[:-4]+"_even.hdf"
				refo=ref[:-4]+"_odd.hdf"
				if options.verbose>1: print(f"{refe} -> alignref_even.hdf    {refo} -> alignref_odd.hdf")
				
				ar=EMData(refe,0)
				ar.write_compressed(f"{options.path}/alignref_even.hdf",0,12,erase=True)

				ar=EMData(refo,0)
				ar.write_compressed(f"{options.path}/alignref_odd.hdf",0,12,erase=True)

		else:
			ar=EMData(ref,0)
			#if (len(options.maskalign)>0):
				#m=EMData(options.maskalign,0)
				#ar.mult(m)
			ar.write_compressed(f"{options.path}/alignref_even.hdf",0,12,erase=True)
			ar.write_compressed(f"{options.path}/alignref_odd.hdf",0,12,erase=True)
		
		refsize=ar["ny"]
		if refsize!=boxsize:
			print("reference and particle have different box size. something must be wrong...")
			exit()
		
		#### generate alignment command first
		### just deal with gold standard references in spt_refine instead of spt_align..
		gd=" --goldcontinue"		
		
		if options.refine:
			gd+=" --refine --maxang {:.1f} --refinetry {:d} ".format(options.maxang,options.refinetry)
			if options.randphi:
				gd+=" --randphi"
			if options.rand180:
				gd+=" --rand180"
			if options.test180:
				gd+=" --test180"
			if itr>startitr:
				ptcls=os.path.join(options.path, "particle_parms_{:02d}.json".format(itr-1))

		if options.transonly: gd+=" --transonly"

		
		if options.breaksym:
			gd+=" --breaksym"
		if options.breaksymsym:
			gd+=" --breaksymsym {}".format(options.breaksymsym)
		
		if options.maxshift>0:
			gd+=" --maxshift {:.1f}".format(options.maxshift)
		
		if options.maskalign!=None:
			gd+=f" --mask {options.maskalign}" 

		cmd="e2spt_align.py {} {}/alignref.hdf --parallel {} --path {} --iter {} --sym {} --minres {} --maxres {} {}".format(ptcls, options.path,  options.parallel, options.path, itr, options.sym, options.minres, curres*.75, gd)
		
		if options.scipy:
			cmd="e2spt_align_subtlt.py {} {}/alignref.hdf --parallel {} --path {} --iter {} --sym {} --minres {} --maxres {} --smooth {} --fromscratch --goldcontinue".format(ptcls, options.path,  options.parallel, options.path, itr, options.sym, options.minres, curres*.75, options.smooth)
		
		ret=run(cmd)
		
		
		s=""
		if options.maxtilt<90.:
			s+=" --maxtilt {:.1f}".format(options.maxtilt)

		# we apply the symmetric subunit mask provided to the current reference and send it to e2spt_average to do a final translational alignment
		if options.symalimask!=None: 
			cmd=f"e2proc3d.py {ref} {options.path}/ref_mono.hdf --multfile {options.symalimask}"
			run(cmd)
			s+=f" --symalimasked={options.path}/ref_mono.hdf"
			
		even=os.path.join(options.path, "threed_{:02d}_even.hdf".format(itr))
		odd=os.path.join(options.path, "threed_{:02d}_odd.hdf".format(itr))
		combine=os.path.join(options.path, "threed_{:02d}.hdf".format(itr))
		
		if options.scipy:
			run(f"e2spa_make3d.py --parallel {options.parallel} --input {options.path}/aliptcls_{itr:02d}.lst --output {even} --keep {options.pkeep:.3f} --clsid 0 --outsize {boxsize}")
			run(f"e2spa_make3d.py --parallel {options.parallel} --input {options.path}/aliptcls_{itr:02d}.lst --output {odd} --keep {options.pkeep:.3f} --clsid 1 --outsize {boxsize}")
		
		else:
			run("e2spt_average.py --parallel {} --path {} --sym {} --keep {:.3f} --iter {} --skippostp {}".format(options.parallel, options.path, options.sym, options.pkeep, itr, s))
		
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
			if options.verbose>1: print("No structure factor so computing simple SF curve for correction")
			
			
		for f in [even, odd]:
			e=EMData(f)
			e.del_attr("xform.align3d")
			e.write_compressed(f.replace("threed_{:02d}_".format(itr), "threed_raw_"),0,12,erase=True)
			if options.setsf==None:
				e.process_inplace("filter.setstrucfac",{"apix":e["apix_x"],"strucfac":sf})
			e.write_compressed(f,0,12,erase=True)
			
		
		s=""
		if options.goldstandard>0:
			s+=" --align"
		
		if options.setsf:
			s+=" --setsf {}".format(options.setsf)
			
		if options.tophat!=None:
			s+=f" --tophat {options.tophat} "
		
		# if we are doing local symmetry refinement or breaking the symmetry
		# it's a bit counterproductive if we then apply symmetry here (as was happening before 8/22/20)
		syms=f"--sym {options.sym}"
		if options.symalimask!=None or options.breaksym: syms=""
		if options.restarget>0: restarget=options.restarget
		else: restarget=curres
		run(f"e2refine_postprocess.py --even {even} --odd {odd} --output {options.path}/threed_{itr:02d}.hdf --iter {itr:d} --tomo --mass {options.mass} --threads {options.threads} --restarget {restarget} {syms} {msk} {s}")

		try: symn=int(options.sym[1:])
		except: symn=0
		if options.symalimask!=None and not options.breaksym and symn>0:
			os.rename(f"{options.path}threed_{itr:02d}.hdf",f"{options.path}/threed_{itr:02d}_nosym.hdf")
			phir=360.0/(symn*2.0)
			if   options.sym[0].lower()=="c":
				run(f"e2proc3d.py {options.path}threed_{itr:02d}_nosym.hdf {options.path}/threed_{itr:02d}.hdf --process mask.cylinder:phicen=0:phirange={phir-5.0}:phitriangle=1:phitrirange=10.0 --sym {options.sym}")
			elif options.sym[0].lower()=="d":
				run(f"e2proc3d.py {options.path}threed_{itr:02d}_nosym.hdf {options.path}/threed_{itr:02d}.hdf --process mask.cylinder:phicen=0:phirange={phir-5.0}:phitriangle=1:phitrirange=10.0:zmin={data['nz']/2}:zmax={data['nz']} --sym {options.sym}")

		ref=os.path.join(options.path, "threed_{:02d}.hdf".format(itr))
		fsc=np.loadtxt(os.path.join(options.path, "fsc_masked_{:02d}.txt".format(itr)))
		
		fi=fsc[:,1]<0.2
		if np.sum(fi)==0:
			print("something wrong with the FSC curve. Cannot estimate resolution. Please check.")
		else:
			rs=1./fsc[fi, 0][0]
			print("Resolution (FSC<0.2) is ~{:.1f} A".format(rs))
		curres=rs
		if curres>50.:
			curres=50
			
		# realign to the initial reference to prevent drifting. Particularly important with alignment masks
		if options.realign:
			print("Realigning to initial reference")
			os.rename(f"{combine}","tmp.hdf")
			run(f"e2proc3d.py tmp.hdf {combine} --alignref {options.path}/model_input.hdf --align rotate_translate_3d_tree --compressbits 12")	# align
			xform=EMData(combine,0,True)["xform.align3d"]		# recover alignment orientation
			
			# this is critical for the next round, but also makes sense for self-consistency
			a=EMData(even)
			a.process_inplace("xform",{"transform":xform})
			a.write_compressed(even,0,12,erase=True)

			a=EMData(odd)
			a.process_inplace("xform",{"transform":xform})
			a.write_compressed(odd,0,12,erase=True)
			
			print("Updating particle orientations from alignment")
			angs=js_open_dict("{}/particle_parms_{:02d}.json".format(options.path,itr))		# now we want to update the particle orientations as well for the next round
			for k in angs.keys():
				parm=angs.get(k,True)
				parm["xform.align3d"]=parm["xform.align3d"]*xform
				angs.setval(k,parm,True)
			angs.sync()

	E2end(logid)
	
def run(cmd):
	print(cmd)
	ret=launch_childprocess(cmd)
	return ret
	
	
if __name__ == '__main__':
	main()

