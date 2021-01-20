#!/usr/bin/env python
# Muyuan Chen 2021-01
from EMAN2 import *
import numpy as np

def main():
	req=[
		"_rlnImageName", "_rlnDefocusU", "_rlnDefocusV",
		"_rlnDefocusAngle","_rlnAnglePsi", "_rlnAngleTilt",
		"_rlnAngleRot","_rlnOriginXAngst","_rlnOriginYAngst",
		]
	
	usage=" read a relion refinement star file and convert to a EMAN list format. currently only read {}".format(','.join(req))
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--skipheader", type=int,help="skip the first N lines of the star file (before the parameter names starting with _). default is 5", default=5)
	parser.add_argument("--voltage", type=int,help="voltage", default=300)
	parser.add_argument("--cs", type=float,help="cs", default=2.7)
	parser.add_argument("--amp", type=float,help="amplitude contrast. default 0", default=0)
	parser.add_argument("--apix", type=float,help="apix", default=1.0)
	parser.add_argument("--output", type=str,help="output lst file", default="sets/all_relion.lst")
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
	
	rid=[-1]*len(req)
	for i, k in enumerate(keys):
		li=l[i]
		for j,r in enumerate(req):
			if k.startswith(r):
				rid[j]=i

	for i,ri in enumerate(rid):
		print(i, ri, keys[ri], l[ri])

	apix=options.apix
	fnames=[]
	
	try:os.mkdir("particles")
	except: pass

	for i in range(len(lines)):
		ln=lines[i]
		if len(ln)<5: 
			continue
		l=ln.split()
		ii,src=l[rid[0]].split('@')
		ii=int(ii)-1
	#	 src="data/"+src
		e=EMData(src, ii)
		
		du, dv, ang=float(l[rid[1]]), float(l[rid[2]]), float(l[rid[3]])

		ctf = EMAN2Ctf()
		ctf.from_dict({"defocus":(du+dv)/20000.,"dfang":ang,
					"dfdiff":abs(du-dv)/10000,"voltage":options.voltage,"cs":options.cs,"ampcont":options.amp,"apix":options.apix})
		
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
		
		fm="particles"+src[src.rfind('/'):src.rfind('.')]+".hdf"
		
		e.write_image(fm,ii)
		fnames.append([i, fm ,ii])
		sys.stdout.write("\r{}/{}	  ".format(i, len(lines)))
		sys.stdout.flush()
	
	print()

	lst=LSXFile(options.output)
	for i in range(len(lines)):
		ln=lines[i]
		if len(ln)<5: 
			continue
		l=ln.split()
		psi,tlt,rot=float(l[rid[4]]), float(l[rid[5]]), float(l[rid[6]])
		tx,ty=float(l[rid[7]])/apix, float(l[rid[8]])/apix
		d={"type":"spider", "psi":psi, "theta":tlt, "phi":rot,"tx":-tx, "ty":-ty}
		d=Transform(d).get_params("eman")
		
		fm=fnames[i]
		lst.write(i, fm[2], fm[1], str(d))
		
	lst=None
	print("aligned list written to {}".format(options.output))
	E2end(logid)
	
def run(cmd):
	print(cmd)
	launch_childprocess(cmd)
	
	
if __name__ == '__main__':
	main()
	