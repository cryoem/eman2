#!/usr/bin/env python
from past.utils import old_div
from builtins import range
# Muyuan Chen 2018-08
from EMAN2 import *
import numpy as np
from EMAN2_utils import *

def main():
	
	usage="""prog --path <spt_xx> --iter <X> --tomo <tomogram>
	map aligned particles back to tomograms """
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--path", type=str,help="spt_xx path", default="",guitype='strbox',row=0, col=0,rowspan=1, colspan=1)
	parser.add_argument("--iter", type=int,help="iteration number", default=1,guitype='intbox',row=0, col=1,rowspan=1, colspan=1)
	parser.add_argument("--tomo", type=str,help="tomogram file name", default="", guitype='filebox', browser="EMBrowserWidget(withmodal=True, startpath='tomograms')", row=2, col=0,rowspan=1, colspan=2,)
	parser.add_argument("--avg", type=str,help="3D average. will use spt_xx/threed_xx by default", default="")
	parser.add_argument("--ptcl", type=str,help="particle input. will read from 0_spt_params by default", default="")
	parser.add_argument("--keep", type=float,help="propotion to keep. will exclude bad particles if this is smaller than 1.0", default=1.0)
	parser.add_argument("--ppid", type=int,help="ppid...", default=-1)
	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)
	
	path=options.path
	itr=options.iter
	
	if len(options.ptcl)>0:
		ptclin=options.ptcl
	else:
		js=js_open_dict("{}/0_spt_params.json".format(path))
		ptclin=str(js["input_ptcls"])
		js=None
	
	js=js_open_dict("{}/particle_parms_{:02d}.json".format(path, itr))
	
	lst=LSXFile(ptclin, True)
	
	tomo=EMData(options.tomo)
	
	if len(options.avg)==0:
		options.avg="{}/threed_{:02d}.hdf".format(path, itr)
	avg=EMData(options.avg, 0)
	print("Using averaged map {}".format(options.avg))
		
	
	apix_tomo=tomo["apix_x"]
	apix_ptcl=avg["apix_x"]
	shrink=apix_tomo/apix_ptcl
	
	print("Apix from tomogram {:.2f}, apix from average {:.2f}. Shrinking factor {:.1f}".format(
		apix_tomo, apix_ptcl, shrink))
	
	avg.process_inplace("math.fft.resample",{"n":shrink})
	avg.process_inplace("normalize.edgemean")
	
	tomo.to_zero()
	
	bname=base_name(options.tomo)
	ptclid=[]
	scr=[]
	for i in range(lst.n):
		l=lst.read(i)
		if bname==base_name(l[1]):
			ky=str((ptclin,i))
			if js.has_key(ky):
				ptclid.append(i)
				val=js[ky]
				scr.append(val["score"])
			
	nptcl=int(len(scr)*options.keep)
	if options.keep<1.0:
		sthr=np.sort(scr)[nptcl]
	else:
		sthr=100
	
	pts=[]
	print("{:d} particles total.".format(int(nptcl)))
	for i in ptclid:
		l=lst.read(i)
		ky=str((ptclin,i))
		val=js[ky]
		s=val["score"]
		if s>sthr:
			continue
		a=EMData(l[1], l[0],True)
		xf=Transform(val["xform.align3d"])
		
		if tomo.has_attr("zshift"):
			zs=tomo["zshift"]
		else:
			zs=0
		crd=np.array(a["ptcl_source_coord"])/shrink + [
			tomo["nx"]//2, tomo["ny"]//2, tomo["nz"]//2 + zs/shrink]
		ts=np.array(xf.get_trans())
		xf.set_trans((ts/shrink).tolist())
		xf=xf.inverse()
		t=avg.process("xform", {"transform":xf})
		pts.append(crd-ts/4)
		tomo.insert_scaled_sum(t,crd.tolist())
		print("\t{}/{} finished.".format(len(pts), nptcl), end='\r')
		sys.stdout.flush()
			

	pts=np.array(pts)
	pfile="{}/ptcls_pos_{:02d}.pdb".format(path, itr)
	tfile="{}/ptcls_in_tomo_{}_{:02d}.hdf".format(path, bname, itr)
	tomo.write_image(tfile)
	numpy2pdb(fname=pfile, data=pts)
	print("Map {} particles to the tomogram".format(len(pts)))
	print("Particle coordinates written to {}".format(pfile))
	print("Map with particles written to {}".format(tfile))
	print("Done.")
	js=None
	E2end(logid)
	
	
	
if __name__ == '__main__':
	main()
	