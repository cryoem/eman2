#!/usr/bin/env python
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
	parser.add_argument("--avg", type=str,help="3D volume to insert. spt_xx/threed_xx if unspecified", default="",guitype='filebox', browser="EMBrowserWidget(withmodal=True, startpath='.')", row=3, col=0,rowspan=1, colspan=2)
	parser.add_argument("--postxf", type=str,help="extra shift after alignment", default="")
	parser.add_argument("--pdb", type=str,help="pdb input", default=None)
	parser.add_argument("--keep", type=float,help="propotion to keep. will exclude bad particles if this is smaller than 1.0", default=1.0)
	parser.add_argument("--gui",action="store_true",help="open the resulting map and tomogram in a GUI display",default=False,guitype="boolbox",row=4, col=0,rowspan=1, colspan=1)
	parser.add_argument("--new",action="store_true",help="Results from e2spt_refine_new",default=False,guitype="boolbox",row=4, col=1,rowspan=1, colspan=1)
	parser.add_argument("--ppid", type=int,help="ppid...", default=-1)
	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)
	
	path=options.path
	itr=options.iter
	try:
		postxf=options.postxf.split(',')
		postxf=[float(i) for i in postxf]
		print("post shift : ", postxf)
	except:
		postxf=[0,0,0]
	
		
	tomo=EMData(options.tomo)
	if options.gui: tomo_orig=tomo.copy()
	
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
	
	if options.pdb:
		pdb=pdb2numpy(options.pdb)
		print(f"load pdb of {len(pdb)} atoms")
		pdb/=avg["apix_x"]
		print(pdb)
		pdb-=np.mean(pdb, axis=0)
		
		#pdb-=avg["nx"]//2
		
	tomo.to_zero()
	
	ptcl=[]
	llfile=""
	bname=base_name(options.tomo)
	if options.new:
		alipm=load_lst_params("{}/aliptcls3d_{:02d}.lst".format(path, itr))
		
		for p in alipm:
			if base_name(p["src"])==bname:
				ptcl.append((p["score"], p["src"], p["idx"], p["xform.align3d"]))
	else:
		js=js_open_dict("{}/particle_parms_{:02d}.json".format(path, itr))
		for k in js.keys():
			fsp,i=eval(k)
			if fsp[-4:]==".lst" :
				if llfile!=fsp : 
					llfile=fsp
					lsx=LSXFile(fsp,True)
				i,fsp,x=lsx.read(i)
			if base_name(fsp)==bname:
				ptcl.append((js[k]["score"],fsp,i,js[k]["xform.align3d"]))
	
	ptcl.sort()
	#print(ptcl)
			
	nptcl=int(len(ptcl)*options.keep)
	if options.keep<1.0:
		sthr=ptcl[nptcl][0]
	else:
		sthr=100
	
	pts=[]
	print("{:d} particles total.".format(int(nptcl)))
	if nptcl==0:
		print("No particles. Exiting")
		sys.exit(0)
		
	nc=0
	for s,fsp,i,xf in ptcl:
		if s>sthr:
			continue
		a=EMData(fsp, i,True)
		
		if tomo.has_attr("zshift"):
			zs=tomo["zshift"]
		else:
			zs=0
		crd=np.array(a["ptcl_source_coord"])/shrink + [
			tomo["nx"]//2, tomo["ny"]//2, tomo["nz"]//2 + zs/shrink]
		ts=np.array(xf.get_trans())
		ts+=postxf
		xf.set_trans((ts/shrink).tolist())
		xf=xf.inverse()
		t=avg.process("xform", {"transform":xf})
		#if a["ptcl_source_coord"][0]<0 or a["ptcl_source_coord"][1]<0:
			#continue
		if options.pdb:
			p1=np.array([xf.transform(p.tolist()) for p in pdb])
			pts.extend(p1+crd)
		else:
			pts.append(crd-ts/4)
			tomo.insert_scaled_sum(t,crd.tolist())
			
		nc+=1
		print("\t{}/{} finished.".format(nc, nptcl), end='\r')
		sys.stdout.flush()
			

	pts=np.array(pts)
	print(pts.shape)
	tfile="{}/ptcls_in_tomo_{}_{:02d}.hdf".format(path, bname, itr)
	if options.pdb:
		pass
	else:
		
		tomo.write_compressed(tfile,0,8)
	
	pfile="{}/ptcls_pos_{:02d}.pdb".format(path, itr)
	numpy2pdb(fname=pfile, data=pts)
	np.savetxt("{}/ptcls_pos_{:02d}.txt".format(path, itr), pts)
	print("Map {} particles to the tomogram".format(nc))
	print("Particle coordinates written to {}".format(pfile))
	print("Map with particles written to {}".format(tfile))
	js=None
	E2end(logid)
	
	if options.gui:
		print("opening GUI")
		from eman2_gui.emapplication import get_application, EMApp
		from eman2_gui.emscene3d import EMScene3D
		from eman2_gui.emdataitem3d import EMDataItem3D, EMIsosurface, EMSliceItem3D
		from eman2_gui.emshape import EMShape

		app = EMApp()
		view=EMScene3D()
		
		# slice through the original tomogram
		tomo_orig_di = EMDataItem3D(tomo_orig, transform=Transform())
		view.insertNewNode('Data', tomo_orig_di, parentnode=view)
		volslice = EMSliceItem3D(tomo_orig_di, transform=Transform())
		view.insertNewNode("Slice", volslice, parentnode=tomo_orig_di)
		
		# isosurface of reconstituted particles
		tomo_di = EMDataItem3D(tomo, transform=Transform())
		view.insertNewNode('Data', tomo_di, parentnode=view)
		isosurface = EMIsosurface(tomo_di, transform=Transform())
		view.insertNewNode("Iso", isosurface, parentnode=tomo_di)

		view.show()
		try: view.raise_()
		except: pass

		app.execute()
			
	print("Done.")

	
	
if __name__ == '__main__':
	main()
	
