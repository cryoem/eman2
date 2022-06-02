#!/usr/bin/env python
# Muyuan Chen 2017-10
from future import standard_library
standard_library.install_aliases()
from builtins import range
from EMAN2 import *
import numpy as np

def main():
	
	usage="""Make 3D particles from subtilt refinement results
	spt_remake3dptcl.py --path <subtlt_xx> --iter <final iteration #> --label <label of new 3D particles>"""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--path", type=str,help="path to subtilt refinement", default="subtlt_00")
	parser.add_argument("--iter", type=int,help="iteration number", default=1)
	parser.add_argument("--label", type=str,help="new particle label", default="")
	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)
	
	threed=EMData("{}/threed_{:02d}.hdf".format(options.path, options.iter))
	params=[]
	
	for eo in ["even", "odd"]:
		lst0=LSXFile("{}/ali_ptcls_00_{}.lst".format(options.path, eo), True)
		lname="{}/ali_ptcls_{:02d}_{}.lst".format(options.path, options.iter, eo)
		lst1=LSXFile(lname, True)
		bx=threed["nx"]
		n=lst0.n
		print("Loading {} particles from {}...".format(n, lname))
		for i in range(n):
			l0=lst0.read(i)
			l1=lst1.read(i)
			dic0=eval(l0[2])
			dic1=eval(l1[2])
			# scr=dic0.pop("score")
			scr=dic1.pop("score")
			xf0=Transform(dic0)
			xf1=Transform(dic1)

			si,src=l0[:2]
			e=EMData(src, si, True)
			xf=xf1*xf0.inverse()
			xpj=e["xform.projection"]
			xali=xpj.inverse()*xf0
			xali.invert()
			pm=[src, si, e["file_threed"], e["model_id"], xf*xpj, xali]
			params.append(pm)
			sys.stdout.write("\r{}/{}".format(i+1, n))
			sys.stdout.flush()
		print()


	if options.label=="":
		p=params[0][0]
		label=p[p.find('__')+2:-4]+"remake"
		print("using label {}".format(label))
	else:
		label=options.label
		
	srcs=sorted(np.unique([t[2] for t in params]))
	p3d=e["nx"]
	apix=e["apix_x"]
	print("Making 3d particles...")
	for src in srcs:
		### per tomogram
		pms=[t for t in params if t[2]==src]
		pids=sorted(np.unique([t[3] for t in pms]))
		print("{}: {} particles".format(src, len(pids)))
		fname="particles3d/{}__{}.hdf".format(base_name(src), label)
		print("  writing to {}".format(fname))
		if os.path.isfile(fname):
			os.remove(fname)
		for ii,ip in enumerate(pids):
			### per 3d particle
			tids=[t for t in pms if t[3]==ip]
			recon=Reconstructors.get("fourier", {"sym":'c1', "size":[p3d, p3d, p3d], "mode":"gauss_2"})
			recon.setup()
			for ti in tids:
				### per tilt
				e0=EMData(ti[0],ti[1])
				xform=ti[4]
	#			 xform=e0["xform.projection"]
				x=xform.get_trans()
				trans=Transform({"type":"2d", "tx":-x[0], "ty":-x[1]})
				e1=recon.preprocess_slice(e0, trans)
				recon.insert_slice(e1,xform,1)
			
			
			threed=recon.finish(True)
			
			threed.process_inplace("math.gausskernelfix",{"gauss_width":4.0})
			threed=threed.get_clip(Region((p3d-bx)//2,(p3d-bx)//2,(p3d-bx)//2,bx,bx,bx))
			threed["apix_x"]=threed["apix_y"]=threed["apix_z"]=apix
			threed["xform.align3d"]=tids[0][5]
			threed["model_id"]=ip
			threed["ptcl_source_coord"]=e0["ptcl_source_coord_3d"]
			threed.write_image(fname, -1)
			
			sys.stdout.write("\r    {}/{}".format(ii+1, len(pids)))
			sys.stdout.flush()
		print()
			
	E2end(logid)
	
def run(cmd):
	print(cmd)
	launch_childprocess(cmd)
	
	
if __name__ == '__main__':
	main()
	