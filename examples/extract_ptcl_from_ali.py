#!/usr/bin/env python
# Muyuan Chen 2015-03
from EMAN2 import *
import numpy as np

def main():
	
	usage=" "
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--tomo", type=str,help="reconstructed tomo name", default=None)
	parser.add_argument("--ali", type=str,help="imod ali file name", default=None)
	parser.add_argument("--tlt", type=str,help="imod tlt file name", default=None)
	parser.add_argument("--xtilt", type=float,help="imod xtilt value (from tomopitch.log)", default=0)
	parser.add_argument("--ptclout", type=str,help="particle output", default=None)
	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)

	tomoname=options.tomo
	js=js_open_dict(info_name(tomoname))
	box=np.array([[b[0],b[1],b[3]] for b in js["boxes"]])
	js=None
	
	tlts=np.loadtxt(options.tlt)
	ptclout=options.ptclout
	try: os.remove(ptclout)
	except: pass

	alifile=options.ali

	sz=32
	e0=EMData(tomoname,0,True)
	boxes=box.copy()
	boxes[:,2]-=e0["nz"]/2
	boxes[:,1]-=e0["ny"]/2
	boxes[:,0]-=e0["nx"]/2
	xtlt=options.xtilt
	allb=[]
	
	for bi,b in enumerate(boxes):
		for i,t in enumerate(tlts):
	
			tr=Transform()
			tr=Transform({"type":"xyz","xtilt":xtlt, "ytilt":-t})
			
			p=tr.transform(Vec3f(b.astype(int).tolist()))
			p[0]+=e0["nx"]/2
			p[1]+=e0["ny"]/2
			allb.append([p[0],p[1], "manual", i])
			e=EMData(alifile, 0, False, Region(p[0]-sz/2,p[1]-sz/2,i,sz,sz,1))
			e.process_inplace("normalize.edgemean")
			e["box"]=[int(p[0]),int(p[1])]
			e["model_id"]=bi
			e.write_image(ptclout,-1)
	js0=js_open_dict(info_name(alifile))
	js0["boxes"]=allb
	js0=None
	
	E2end(logid)
	
def run(cmd):
	print cmd
	launch_childprocess(cmd)
	
	
if __name__ == '__main__':
	main()
	