#!/usr/bin/env python
# Muyuan Chen 2017-03
from __future__ import print_function
from __future__ import division
from past.utils import old_div
from builtins import range
from EMAN2 import *
import numpy as np
from EMAN2_utils import *

def main():
	
	usage=" "
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_pos_argument(name="particles",help="Specify particles to use to generate an initial model.", default="", guitype='filebox', browser="EMSetsTable(withmodal=True,multiselect=False)", row=0, col=0,rowspan=1, colspan=3, mode="model")
	parser.add_argument("--reference","--ref",help="""3D reference for initial model generation. No reference is used by default.""", default="", guitype='filebox', browser="EMBrowserWidget(withmodal=True,multiselect=False)", row=1, col=0,rowspan=1, colspan=3, mode="model")

	parser.add_header(name="orblock1", help='Just a visual separation', title="Options", row=2, col=1, rowspan=1, colspan=1, mode="model")

	parser.add_argument("--mask", type=str,help="Mask file to be applied to initial model", default="", guitype='filebox', browser="EMBrowserWidget(withmodal=True,multiselect=False)", row=3, col=0,rowspan=1, colspan=3, mode="model")

	parser.add_argument("--niter", type=int,help="Number of iterations", default=5, guitype='intbox',row=4, col=0,rowspan=1, colspan=1, mode="model")
	parser.add_argument("--sym", type=str,help="symmetry", default="c1", guitype='strbox',row=4, col=1,rowspan=1, colspan=1, mode="model")
	#parser.add_argument("--gaussz", type=float,help="Extra gauss filter at z direction", default=-1, guitype='floatbox',row=4, col=2,rowspan=1, colspan=1, mode="model")

	parser.add_argument("--mass", type=float,help="mass", default=500.0, guitype='floatbox',row=5, col=0,rowspan=1, colspan=1, mode="model")
	#parser.add_argument("--tarres", type=float,help="target resolution", default=10, guitype='floatbox',row=5, col=1,rowspan=1, colspan=1, mode="model")
	parser.add_argument("--localfilter", action="store_true", default=False ,help="use tophat local", guitype='boolbox',row=5, col=2,rowspan=1, colspan=1, mode="model")

	parser.add_argument("--goldstandard", type=int,help="initial resolution for gold standard refinement", default=-1, guitype='intbox',row=6, col=0,rowspan=1, colspan=1, mode="model")
	parser.add_argument("--goldcontinue", action="store_true", default=False ,help="continue from an existing gold standard refinement", guitype='boolbox',row=6, col=1,rowspan=1, colspan=1, mode="model")

	parser.add_argument("--setsf", type=str,help="structure factor", default=None, guitype='filebox', browser="EMBrowserWidget(withmodal=True,multiselect=False)", row=7, col=0,rowspan=1, colspan=3, mode="model")

	parser.add_argument("--pkeep", type=float,help="fraction of particles to keep", default=0.8, guitype='floatbox',row=8, col=0,rowspan=1, colspan=1, mode="model")

	#parser.add_argument("--tkeep", type=float,help="fraction of tilt to keep. only used when tltrefine is specified. default is 0.5", default=0.5, guitype='floatbox',row=8, col=1,rowspan=1, colspan=1, mode="model")

	parser.add_argument("--maxtilt",type=float,help="Explicitly zeroes data beyond specified tilt angle. Assumes tilt axis exactly on Y and zero tilt in X-Y plane. Default 90 (no limit).",default=90.0, guitype='floatbox',row=8, col=2,rowspan=1, colspan=1, mode="model")

	#parser.add_argument("--tltrefine", action="store_true", default=False ,help="refine individual subtilts", guitype='boolbox',row=9, col=0,rowspan=1, colspan=1, mode="model")
	#parser.add_argument("--tltrefine_unmask", action="store_true", default=False ,help="use unmasked maps for tilt refinement", guitype='boolbox',row=9, col=1,rowspan=1, colspan=1, mode="model")

	parser.add_argument("--threads", type=int,help="threads", default=12, guitype='intbox',row=9, col=0,rowspan=1, colspan=1, mode="model")

	parser.add_argument("--path", type=str,help="Specify name of refinement folder. Default is spt_XX.", default=None)#, guitype='strbox', row=10, col=0,rowspan=1, colspan=3, mode="model")
	parser.add_argument("--wtori",type=float,help="Weight for using the prior orientation in the particle header. default is -1, i.e. not used.",default=-1)

	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-2)

	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)
	ptcls=args[0]
	ref=options.reference

	if options.path==None: options.path = make_path("spt") 
	
	#### make a list file if the particles are not in a lst
	if not ptcls.endswith(".lst"):
		ptcllst="{}/input_ptcls.lst".format(options.path)
		run("e2proclst.py {} --create {}".format(ptcls, ptcllst))
		ptcls=ptcllst

	options.input_ptcls=ptcls
	options.input_ref=ref
	options.cmd=' '.join(sys.argv)
	js=js_open_dict("{}/0_spt_params.json".format(options.path))
	js.update(vars(options))
	js.close()
	
	ep=EMData(ptcls,0)
	
	if options.goldcontinue==False:
		er=EMData(ref,0)
		if abs(1-ep["apix_x"]/er["apix_x"])>0.01 or ep["nx"]!=er["nx"]:
			print("apix mismatch {:.2f} vs {:.2f}".format(ep["apix_x"], er["apix_x"]))
			rs=er["apix_x"]/ep["apix_x"]
			if rs>1.:
				run("e2proc3d.py {} {}/model_input.hdf --clip {} --scale {} --process mask.soft:outer_radius=-1".format(ref, options.path, ep["nx"], rs))
			else:
				run("e2proc3d.py {} {}/model_input.hdf --scale {} --clip {} --process mask.soft:outer_radius=-1".format(ref, options.path, rs, ep["nx"]))
			
			ref="{}/model_input.hdf".format(options.path)
		
	
	curres=40.
	for itr in range(1,options.niter+1):
		
		#### generate alignment command first
		gd=""
		if options.goldstandard>0 and itr==1:
			gd=" --goldstandard {}".format(options.goldstandard)
		if options.goldcontinue or (options.goldstandard>0 and itr>1):
			gd=" --goldcontinue".format(options.goldstandard)
			
		#curres=0
		cmd="e2spt_align.py {} {} --threads {} --path {} --iter {} --sym {} --wtori {} {}".format(ptcls, ref,  options.threads, options.path, itr, options.sym, options.wtori, gd)
		
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
			
		
		run("e2spt_average.py --threads {} --path {} --sym {} --skippostp {}".format(options.threads, options.path, options.sym, s))
		
		msk=options.mask
		if len(msk)>0:
			if os.path.isfile(msk):
				msk=" --automask3d mask.fromfile:filename={}".format(msk)
			else:
				msk=" --automask3d {}".format(msk)
		#msk=" --automask3d mask.fromfile:filename=Subtomograms/atpase_mask.hdf"
		#msk=" --automask3d mask.cylinder:outer_radius=14:phirange=360:zmax=50.0:zmin=14.0 --automask3d2 mask.addshells.gauss:val1=0:val2=8"
		s=""
		if options.goldstandard>0:
			s+=" --align"
		
		if options.setsf:
			s+=" --setsf {}".format(options.setsf)
			
		if options.localfilter:
			s+=" --tophat local "
			
		#if options.gaussz>0:
			#s+=" --m3dpostprocess filter.lowpass.gaussz:cutoff_freq={}".format(options.gaussz)
			
		os.system("rm {}/mask*.hdf {}/*unmasked.hdf".format(options.path, options.path))
		ppcmd="e2refine_postprocess.py --even {} --odd {} --output {} --iter {:d} --mass {} --restarget {} --threads {} --sym {} {} {} ".format(os.path.join(options.path, "threed_{:02d}_even.hdf".format(itr)), os.path.join(options.path, "threed_{:02d}_odd.hdf".format(itr)), os.path.join(options.path, "threed_{:02d}.hdf".format(itr)), itr, options.mass, curres, options.threads, options.sym, msk, s)
		run(ppcmd)
		
		ref=os.path.join(options.path, "threed_{:02d}.hdf".format(itr))
		fsc=np.loadtxt(os.path.join(options.path, "fsc_masked_{:02d}.txt".format(itr)))
		rs=1./fsc[fsc[:,1]<0.3, 0][0]
		print("Resolution (FSC<0.3) is ~{:.1f} A".format(rs))
		curres=rs
		if curres>40.:
			curres=40
		
		# if options.tltrefine:# and itr%2==0:
			
		# 	os.system("rm {}/mask*.hdf {}/*unmasked.hdf".format(options.path, options.path))
		# 	cmd="e2spt_tltrefine.py --path {} --iter {} --threads {} --sym {} --maxres {} --keep {}".format(options.path, itr, options.threads, options.sym, rs, options.tkeep)
		# 	if options.tltrefine_unmask:
		# 		cmd+=" --unmask "
		# 	run(cmd)
		# 	#ppcmd=ppcmd.replace("threed_{:02d}".format(itr), "threed_{:02d}_ali".format(itr))
		# 	run(ppcmd)
		# 	for eo in ["", "_even", "_odd"]:
		# 		os.rename("{}/threed_{:02d}{}.hdf".format(options.path, itr, eo), 
		# 				"{}/threed_{:02d}_ali{}.hdf".format(options.path, itr, eo))
		# 	refnew=os.path.join(options.path, "threed_{:02d}_ali.hdf".format(itr))
		# 	if os.path.isfile(refnew):
		# 		ref=refnew
		# 		fscs=[os.path.join(options.path,f) for f in os.listdir(options.path) if f.startswith("fsc") and f.endswith("{:02d}.txt".format(itr))]
		# 		for f in fscs:
		# 			os.rename(f, f[:-4]+"_ali.txt")
		# 	else:
		# 		print("tilt refinement failed for some reason...")
		# 		return
			
	

	E2end(logid)
	
def run(cmd):
	print(cmd)
	ret=launch_childprocess(cmd)
	return ret
	
	
if __name__ == '__main__':
	main()

