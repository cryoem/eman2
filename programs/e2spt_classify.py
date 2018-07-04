#!/usr/bin/env python
# Muyuan Chen 2016-10

from EMAN2 import *
import numpy as np
from shutil import copyfile

def main():
	
	usage=" "
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_pos_argument(name="particles",help="Specify particles to use to generate an initial model.", default="", guitype='filebox', browser="EMSPTParticleTable(withmodal=True,multiselect=False)", row=0, col=0, rowspan=1, colspan=3,mode="multi")
	parser.add_argument("--references","--refs", type=str,help="3D reference volumes", default=None, guitype="filebox", browser="EMBrowserWidget(withmodal=True,multiselect=False)", row=1, col=0,rowspan=1, colspan=3,mode="multi")

	parser.add_header(name="orblock1", help='Just a visual separation', title="Options", row=2, col=1, rowspan=1, colspan=1,mode="multi")

	parser.add_argument("--mask", type=str,help="mask", default='', guitype="filebox", browser="EMBrowserWidget(withmodal=True,multiselect=False)", row=3, col=0, rowspan=1, colspan=3, mode="multi")
	parser.add_argument("--niter", type=int,help="iterations", default=3, guitype="intbox", row=4, col=0,rowspan=1, colspan=1,mode="multi")
	parser.add_argument("--sym", type=str,help="sym", default="c1", guitype="strbox", row=4, col=1,rowspan=1, colspan=1,mode="multi")
	parser.add_argument("--tarres", type=float,help="target resolution", default=20.0, guitype="floatbox", row=4, col=2,rowspan=1, colspan=1,mode="multi")
	parser.add_argument("--mass", type=float,help="mass", default=500, guitype="floatbox", row=5, col=0,rowspan=1, colspan=1,mode="multi")
	parser.add_argument("--localfilter", action="store_true", default=False ,help="use tophat local", guitype="boolbox", row=5, col=1,rowspan=1, colspan=1,mode="multi")
	parser.add_argument("--threads", type=int,help="threads", default=12, guitype="intbox", row=5, col=2,rowspan=1, colspan=1,mode="multi")
	parser.add_argument("--path", type=str,help="path", default=None)
	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)
	
		
	ptcls=args[0]
	refs=options.refs.split(',')
	
	
	if options.path==None:
		for i in range(100):
			pname="spt_{:02d}".format(i)
			if not os.path.isdir(pname):
				os.mkdir(pname)
				options.path=pname
				break
		else:
			print "something is wrong..."
			exit()
	else:
		try: 
			os.mkdir(options.path)
		except:
			pass
	
	path=options.path
	print "Working in path {}...".format(path)
	
	for i,r in enumerate(refs):
		rf="{}/init_ref_{:02d}.hdf".format(path, i)
		run("e2proc3d.py {} {} --process filter.lowpass.randomphase:cutoff_freq={} --process filter.lowpass.gauss:cutoff_freq={} --process normalize".format(
			r, rf, 0.02, 0.02))
		refs[i]=rf
	
	for itr in range(1,options.niter+1):
		jsali="{}/particle_parms_{:02d}.json".format(path, itr)
		jscls=["{}/particle_parms_{:02d}_cls{:02d}.json".format(path, itr, i) for i in range(len(refs))]
		for ir, rf in enumerate(refs):
			cmd="e2spt_align.py {} {} --threads {} --path {} --iter {} --sym {}".format(ptcls, rf,  options.threads, path, itr, options.sym)
			
			
			### in case e2spt_align get segfault....
			ret=1
			while ret>0:
				try: os.remove(jsali)
				except:pass

				ret=run(cmd)
				
			os.rename(jsali, jscls[ir])
		
		
		jss=[js_open_dict(j) for j in jscls]
		ks=[np.asarray(sorted(js.keys())) for js in jss]
		

		if np.sum(ks[0]!=ks[1])>0:
			print "Something wrong..."
			print [len(k) for k in ks]
		else:
			keys=ks[0]
		
		score=[]
		for i in range(len(refs)):
			score.append([jss[i][k]["score"] for k in keys])
		score=np.array(score)
		
		#print "score shape:", score.shape
		clsmx=np.argmin(score, axis=0)
		u,c=np.unique(clsmx, return_counts=True)
		print "Pariticles in each class: ", ', '.join(c.astype(str).tolist())
		
		np.savetxt("{}/classmx_{:02d}.txt".format(path, itr), clsmx)
		
		for ir, rf in enumerate(refs):
			try: os.remove(jsali)
			except:pass
		
			clsid=(score[ir,:]==np.min(score,axis=0))
			ks=keys[clsid]
			dic={}
			for k in ks:
				dic[k]=jss[ir][k]
			#print dic
			js=js_open_dict(jsali)
			js.update(dic)
			js=None
			
			run("e2spt_average.py --threads {} --path {} --sym {} --iter {}".format(options.threads, options.path, options.sym, itr))
			
			fnms=["{}/threed_{:02d}.hdf".format(path, itr),
				"{}/threed_{:02d}_even.hdf".format(path, itr),
				"{}/threed_{:02d}_odd.hdf".format(path, itr),
				"{}/fsc_masked_{:02d}.txt".format(path, itr),
				"{}/fsc_unmasked_{:02d}.txt".format(path, itr),
				"{}/fsc_maskedtight_{:02d}.txt".format(path, itr),
				"{}/fscvol_{:02d}.hdf".format(path, itr)]
			
			othercmd=""
			if options.localfilter:
				othercmd+=" --tophat local "
			msk=options.mask
			if len(msk)>0:
				if os.path.isfile(msk):
					msk=" --automask3d mask.fromfile:filename={}".format(msk)
				else:
					msk=" --automask3d {}".format(msk)
			othercmd+=msk
				
			ppcmd="e2refine_postprocess.py --even {} --odd {} --output {} --iter {:d} --mass {} --restarget {} --threads {} --sym {} {} ".format(fnms[1], fnms[2], fnms[0], itr, options.mass, options.tarres, options.threads, options.sym, othercmd)
			run(ppcmd)
		
			
			
			for fm in fnms:
				if os.path.isfile(fm):
					os.rename(fm, fm.replace(".", "_cls{:02d}.".format(ir)))
			
		refs=["{}/threed_{:02d}_cls{:02d}.hdf".format(path, itr, ir) for ir in range(len(refs))]
		
		

	E2end(logid)
	
def run(cmd):
	print cmd
	launch_childprocess(cmd)
	
	
if __name__ == '__main__':
	main()
	