#!/usr/bin/env python
# Muyuan Chen 2021-01
from EMAN2 import *
import numpy as np

def main():
	req={
		"name":"_rlnImageName", 
		"mname":"_rlnMicrographName", 
		"dfu":"_rlnDefocusU", 
		"dfv":"_rlnDefocusV",
		"dfang":"_rlnDefocusAngle",
		"psi":"_rlnAnglePsi", 
		"tilt":"_rlnAngleTilt",
		"rot":"_rlnAngleRot",
		"txa":"_rlnOriginXAngst",
		"tya":"_rlnOriginYAngst",
		"tx":"_rlnOriginX",
		"ty":"_rlnOriginY",
		"class":"_rlnRandomSubset",
		}
	
	usage=" read a relion refinement star file and convert to a EMAN list format. currently only read {}".format(','.join(req.values()))
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--skipheader", type=int,help="skip the first N lines of the star file (before the parameter names starting with _). default is auto", default=-1)
	parser.add_argument("--voltage", type=int,help="voltage", default=-1)
	parser.add_argument("--cs", type=float,help="cs", default=-1)
	parser.add_argument("--amp", type=float,help="amplitude contrast (0-100)", default=-1)
	parser.add_argument("--apix", type=float,help="apix", default=-1)
	parser.add_argument("--output", type=str,help="output lst file", default="sets/all_relion.lst")
	parser.add_argument("--shrink", type=str,help="shrink factor", default=None)
	parser.add_argument("--head", type=int,help="use only first N particles. for testing", default=-1)
	parser.add_argument("--clip", type=int,help="clip to size after ctf before saving", default=-1)
	parser.add_argument("--skipwrite", action="store_true", default=False ,help="skip file writing")
	parser.add_argument("--invert", action="store_true", default=False ,help="invert contrast")
	parser.add_argument("--onestack", type=str, default=None ,help="save all particles in one stack")
	parser.add_argument("--lst2star", type=str, default=None ,help="convert the xform from a lst file back to the star file in the same order. ")
	parser.add_argument("--lst2star_subset", type=str, default=None ,help="convert only a subset to star. still need to specify the original lst output in --lst2star")
	parser.add_argument("--make3d", action="store_true", default=False ,help="do 3d reconstruction after import")
	parser.add_argument("--sym", type=str, default="c1" ,help="symmetry for make3d. ")

	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)
	
	if options.lst2star:
		options.skipwrite=True
		options.output=None
		
	if options.shrink:
		shrink=options.shrink.split(',')
		
	f=open(args[0],'r')
	alllines=f.readlines()
	f.close()
	
	############
	has_optics=False
	opt_keys=["_rlnVoltage", "_rlnSphericalAberration", "_rlnAmplitudeContrast", "_rlnImagePixelSize"]
	opt_idx=[-1 for k in opt_keys]
	for i in range(100):
		m=alllines[i]
		if m.startswith("data_optics"):
			has_optics=True
			continue
			
		if has_optics:
			for j,k in enumerate(opt_keys):
				if m.startswith(k):
					opt_idx[j]=int(m.split()[1][1:])-1
					continue
				
			if max(opt_idx)>=0 and not m.startswith("_"):
				break
					
	if min(opt_idx)<0:
		print(opt_idx)
		print("Cannot find all optics parameters. please specify through options")
		print("Required keys are:")
		print('\n  '.join(opt_keys))
		
	else:
		opt=alllines[i].split()
		if options.voltage<0: 
			options.voltage=float(opt[opt_idx[0]])
		print(f"Voltage: {options.voltage}")
		
		if options.cs<0: 
			options.cs=float(opt[opt_idx[1]])
		print(f"CS: {options.cs}")
		
		if options.amp<0: 
			options.amp=float(opt[opt_idx[2]])*100
		print(f"Amp: {options.amp}")
		
		if options.apix<0: 
			options.apix=float(opt[opt_idx[3]])
		print(f"Apix: {options.apix}")
		
		if options.skipheader<0:
			options.skipheader=i
		print(f"Start reading data from line {i}")
	
	############
	keys=[]
	starti=options.skipheader
	for i in range(starti,1000):
		if alllines[i].startswith('_'):
			keys.append(alllines[i][:-1])
		elif len(keys)>0 and i>starti+5:
			starti=i
			break
			
	lines=alllines[starti:]
	l=lines[0].split()
	print(keys)
	rid={}
	for r in req.keys():
		
		for i, k in enumerate(keys):
			if k.startswith(req[r]+' '):
				rid[r]=i
				break
		else:
			print("key {} does not exist".format(req[r]))

	print("##################")
	print("found keys:")
	for r in rid.keys():
		print(r, req[r], rid[r])
		
	if "tx" in rid:
		txa=False
		print("translation in pixel")
	else:
		txa=True
		print("translation in angstrom")

	apix=options.apix
	fnames=[]
	
	try:os.mkdir("particles")
	except: pass
	
	if options.head>0:
		nptcl=options.head
	else:
		nptcl=len(lines)
		
	if options.skipwrite==False and options.onestack:
		try: os.remove(options.onestack)
		except: pass
					  
	src0=""
	mcount={}
	for i in range(nptcl):
		ln=lines[i]
		if len(ln)<5: 
			continue
		l=ln.split()
		if "name" in rid:
			ii,src=l[rid["name"]].split('@')
			
		elif "mname" in rid:
			ky=l[rid["mname"]]
			if ky not in mcount:
				mcount[ky]=0
				
			mcount[ky]+=1 ## since star index start from 1
			
			src=ky
			ii=mcount[ky]
		
		else:
			print("error: no micrograph or particle file name")
			exit()
			
		ii=int(ii)-1
		if options.onestack:
			fm=options.onestack
			io=i
		else:
			io=ii
			if '/' in src:
				fm="particles"+src[src.rfind('/'):src.rfind('.')]+"__flip.hdf"
			else:
				fm="particles/"+src[:src.rfind('.')]+"__flip.hdf"
			
		if not os.path.isfile(src):
			print("{} does not exist".format(src))
			continue

		if not options.skipwrite:
			if src!=src0:
				e=EMData(src, 0, True)
				nx=e["nx"]
				ny=e["ny"]
				nz=e["nz"]
				
			src0=src
			#try:	
			if src.endswith('.mrc') and nz>1:
				e=EMData(src, 0, False, Region(0, 0, ii, nx, ny, 1))
			else:
				e=EMData(src, ii)
			#except:
				#print("something wrong with {} image {}".format(src, ii))
				#exit()
			
			du, dv, ang=float(l[rid["dfu"]]), float(l[rid["dfv"]]), float(l[rid["dfang"]])

			ctf = EMAN2Ctf()
			ctf.from_dict({"defocus":(du+dv)/20000.,"dfang":ang,
					"dfdiff":abs(du-dv)/10000,"voltage":options.voltage,
					"cs":options.cs,"ampcont":options.amp,"apix":options.apix})
			
			fft1=e.do_fft()
			flipim=fft1.copy()
			ctf.compute_2d_complex(flipim,Ctf.CtfType.CTF_SIGN)
			fft1.mult(flipim)
			e=fft1.do_ift()
			e["ctf"]=ctf
			e.process_inplace("normalize.edgemean")
			if options.invert: e.mult(-1)
			e.process_inplace("mask.soft",{"outer_radius":-4,"width":4})
			e["src"]=src
			e["srci"]=ii
			e["apix_x"]=options.apix
			e["apix_y"]=options.apix
			
			if options.clip>0:
				e.clip_inplace(Region((nx-options.clip)//2, (ny-options.clip)//2, options.clip, options.clip))
			e.write_compressed(fm, io,12,nooutliers=True)
			if options.shrink:
				for s in shrink:
					fms="{}_bin{}.hdf".format(fm[:-4],s.replace('.','_'))
					ex=e.copy()
					ex.process_inplace("math.fft.resample",{"n":float(s)})
					ex.process_inplace("normalize.edgemean")
					ex.write_compressed(fms, io,12,nooutliers=True)
			
		fnames.append([i, fm ,io])
		if i%100==0 or i==len(lines)-1 : sys.stdout.write("\r{}/{}	  ".format(i, len(lines)))
		sys.stdout.flush()
	
	print()

	try:os.mkdir("sets")
	except: pass
	#lst=LSXFile(options.output)
	towrite=[]
	for fm in fnames:
		i=fm[0]
		ln=lines[i]
		if len(ln)<5: 
			continue
		
		dic={"src":fm[1], "idx":fm[2]}
		l=ln.split()
		
		if "psi" in rid:
			psi,tlt,rot=float(l[rid["psi"]]), float(l[rid["tilt"]]), float(l[rid["rot"]])
			if txa:
				try: tx,ty=float(l[rid["txa"]])/apix, float(l[rid["tya"]])/apix
				except: tx,ty=0,0		# for symmetry replicated files, this may be right?
			else:
				tx,ty=float(l[rid["tx"]]), float(l[rid["ty"]])
				
			d={"type":"spider", "psi":psi, "theta":tlt, "phi":rot,"tx":-tx, "ty":-ty}
			d=Transform(d)
			dic["xform.projection"]=d
		
		#fm=fnames[i]
		if "class" in rid:
			dic["class"]=int(l[rid["class"]])-1
		towrite.append(dic)
		#lst.write(i, fm[2], fm[1], str(d))
		
	if options.output:
		save_lst_params(towrite, options.output)
		pout=options.output
		if options.shrink:
			for s in shrink:
				ss=float(s)
				if options.onestack:
					tg="bin"+s.replace('.','_')
					run(f"e2proclst.py {options.output} --create {options.output[:-4]}_{tg}.lst")
					run(f"e2proclst.py {options.output[:-4]}_{tg}.lst --refile {options.onestack[:-4]}_{tg}.hdf --scale {1./ss}")
					
				else:
					tg="flip_bin"+s.replace('.','_')
					run(f"e2proclst.py {options.output} --create {options.output[:-4]}_{tg}.lst")
					run(f"e2proclst.py {options.output[:-4]}_{tg}.lst --retype {tg} --scale {1./ss}")
					
				pout=f"{options.output[:-4]}_{tg}.lst"
				
		print("aligned list written to {}".format(options.output))
			
		if options.make3d:
			path=num_path_new("r3d_")
			run(f"e2proclst.py {pout} --create {path}/ptcls_00.lst")
			for eo in ["even", "odd"]:
				run(f"e2spa_make3d.py --input {path}/ptcls_00.lst --output {path}/threed_00_{eo}.hdf --keep 1 --parallel thread:32 --clsid {eo} --sym {options.sym}")
				run(f"e2proc3d.py {path}/threed_00_{eo}.hdf {path}/threed_raw_{eo}.hdf")
			
			if not os.path.isfile("sf.txt"):
				run(f"e2spt_structfac.py --even {path}/threed_raw_even.hdf --res 5")
			
			run(f"e2refine_postprocess.py --even {path}/threed_00_even.hdf --res 5 --setsf sf.txt --tophat localwiener --sym {options.sym} --thread 32")
			
			
	if options.lst2star:
		lst=load_lst_params(options.lst2star)
		print("{} particles from star, {} particles from lst".format(len(fnames), len(lst)))
		if options.lst2star_subset:
			lst_sub=load_lst_params(options.lst2star_subset)
			print("selecting {} particles in subset".format(len(lst_sub)))
			ldic={"{}_{}".format(base_name(l["src"]), l["idx"]):i for i,l in enumerate(lst_sub)}
			
		oname=args[0][:-5]+"_from_eman.star"
		f=open(oname,'w')
		for ln in alllines[:starti]:
			f.write(ln)
		for fm in fnames:
			i=fm[0]
			ln=lines[i]
			l=ln.split()
			
			ls=lst[i]
			if options.lst2star_subset:
				ky="{}_{}".format(base_name(ls["src"]), ls["idx"])
				if ky in ldic:
					k=ldic[ky]
					ls=lst_sub[k]
				else:
					continue
			
			d=Transform(ls["xform.projection"])
			d=d.get_params("spider")
			
			
			c=[f"{x:.6f}" for x in (d["psi"], d["theta"], d["phi"])]
			l[rid["psi"]], l[rid["tilt"]], l[rid["rot"]]=c[0], c[1], c[2]
			if txa:
				c=[f"{x:.6f}" for x in (-d["tx"]*apix, -d["ty"]*apix)]
				l[rid["txa"]],l[rid["tya"]]=c[0], c[1]
				
			else:
				c=[f"{x:.6f}" for x in (-d["tx"], -d["ty"])]
				l[rid["tx"]],l[rid["ty"]]=c[0], c[1]
			
			lo=' '.join(l)+'\n'
			#print(lo)
									
			f.write(lo)
			
		f.close()
		print(oname)
		
		
	E2end(logid)
	
def run(cmd):
	print(cmd)
	launch_childprocess(cmd)
	
	
if __name__ == '__main__':
	main()
	
