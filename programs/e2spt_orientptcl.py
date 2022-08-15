#!/usr/bin/env python
# Muyuan Chen 2022-08
from EMAN2 import *
import numpy as np
import scipy.spatial.distance as scidist

def main():
	
	usage="""
	"""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--ptcllabel", type=str,help="partcle label", default=None)
	parser.add_argument("--cntlabel", type=str,help="label of center coordinate particles", default=None)
	parser.add_argument("--output", type=str,help=".lst output file", default=None)
	parser.add_argument("--ppid", type=int,help="ppid", default=-2)

	
	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)
	
	pnames=[f for f in os.listdir("particles3d/") if f.endswith(options.ptcllabel+".hdf")]
	pnames=[f"particles3d/{f}" for f in pnames]
	allplst=[]

	for pid, pname in enumerate(pnames):
	
		n=EMUtil.get_image_count(pname)
		cname=pname.replace(options.ptcllabel+".hdf",options.cntlabel+".hdf")
		hdrs=[EMData(pname, i, True) for i in range(n)]
		if not os.path.isfile(cname):
			print("error: no center particle in {}".format(cname))
			return
		
		chdr=EMData.read_images(cname)
		print("{} : {} particles, {} centers".format(pname, len(hdrs), len(chdr)))

		coord=[h["ptcl_source_coord"] for h in hdrs]
		coord=np.array(coord)
		cnt=[h["ptcl_source_coord"] for h in chdr]
		
		cnt=np.array(cnt)

		dst=scidist.cdist(coord, cnt)
		di=np.argmin(dst, axis=1)
		dr=coord-cnt[di]
		dr=dr/np.linalg.norm(dr, axis=1)[:,None]

		plst=[]
		for i,r in enumerate(dr):
				tf_dir=Transform()
				tf_dir.set_rotation(r.tolist())
				tf_dir.invert()
				l={"src":pname, "idx":i, "xform.align3d":tf_dir}
				plst.append(l)
			
		allplst.extend(plst)
		
	if options.output==None:
		options.output="sets/"+options.ptcllabel+"_xf.lst"
		
	save_lst_params(allplst, options.output)
	print("Output saved to {}".format(options.output))
	E2end(logid)
	
	
	
if __name__ == '__main__':
	main()
	
