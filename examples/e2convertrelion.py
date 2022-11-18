#!/usr/bin/env python
# Muyuan Chen 2021-01
from EMAN2 import *
import numpy as np

def main():
	req={
		"name":"_rlnImageName", 
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
		}
	
	usage=" read a relion refinement star file and convert to a EMAN list format. currently only read {}".format(','.join(req.values()))
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--skipheader", type=int,help="skip the first N lines of the star file (before the parameter names starting with _). default is 1", default=1)
	parser.add_argument("--voltage", type=int,help="voltage", default=300)
	parser.add_argument("--cs", type=float,help="cs", default=2.7)
	parser.add_argument("--amp", type=float,help="amplitude contrast. default 0 (0-100)", default=10)
	parser.add_argument("--apix", type=float,help="apix", default=1.0)
	parser.add_argument("--output", type=str,help="output lst file", default="sets/all_relion.lst")
	parser.add_argument("--head", type=int,help="use only first N particles. for testing", default=-1)
	parser.add_argument("--skipwrite", action="store_true", default=False ,help="skip file writing")
	parser.add_argument("--onestack", type=str, default=None ,help="save all particles in one stack")

	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)
	
	
	f=open(args[0],'r')
	lines=f.readlines()
	f.close()
	
	keys=[]
	starti=options.skipheader
	for i in range(starti,1000):
		if lines[i].startswith('_'):
			keys.append(lines[i][:-1])
		elif i>starti+5:
			starti=i
			break
			
	lines=lines[starti:]
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
					  
	for i in range(nptcl):
		ln=lines[i]
		if len(ln)<5: 
			continue
		l=ln.split()
		ii,src=l[rid["name"]].split('@')
		ii=int(ii)-1
		if options.onestack:
			fm=options.onestack
			io=i
		else:
			io=ii
			if '/' in src:
				fm="particles"+src[src.rfind('/'):src.rfind('.')]+".hdf"
			else:
				fm="particles/"+src[:src.rfind('.')]+".hdf"
			
		if not os.path.isfile(src):
			print("{} does exist".format(src))
			continue

		if not options.skipwrite:
			#print(src,ii)
			try:
				e=EMData(src, ii)
			except:
				print("something wrong with {} image {}".format(src, ii))
				exit()
			
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
			e.process_inplace("mask.soft",{"outer_radius":-4,"width":4})
			e["src"]=src
			e["srci"]=ii
			e["apix_x"]=options.apix
			e["apix_y"]=options.apix
			
			e.write_compressed(fm, io,12,nooutliers=True)
			
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
		l=ln.split()
		psi,tlt,rot=float(l[rid["psi"]]), float(l[rid["tilt"]]), float(l[rid["rot"]])
		if txa:
			try: tx,ty=float(l[rid["txa"]])/apix, float(l[rid["tya"]])/apix
			except: tx,ty=0,0		# for symmetry replicated files, this may be right?
		else:
			tx,ty=float(l[rid["tx"]]), float(l[rid["ty"]])
			
		d={"type":"spider", "psi":psi, "theta":tlt, "phi":rot,"tx":-tx, "ty":-ty}
		d=Transform(d)
		
		#fm=fnames[i]
		dic={"src":fm[1], "idx":fm[2], "xform.projection":d}
		towrite.append(dic)
		#lst.write(i, fm[2], fm[1], str(d))
		
	
	save_lst_params(towrite, options.output)
	print("aligned list written to {}".format(options.output))
	E2end(logid)
	
def run(cmd):
	print(cmd)
	launch_childprocess(cmd)
	
	
if __name__ == '__main__':
	main()
	
