#!/usr/bin/env python
# Muyuan Chen 2022-08
from EMAN2 import *
import numpy as np
import scipy.spatial.distance as scipydist
from scipy.optimize import minimize

def gather_metadata(ptcls):
	info3d=[]
	info2d=[]
	data=load_lst_params(ptcls)
	for ii,dt in enumerate(data):
		img=EMData(dt["src"], dt["idx"], True)
		imgsrc=img["class_ptcl_src"]
		imgidx=img["class_ptcl_idxs"]

		try: rhdrs=EMData.read_images(imgsrc,imgidx,IMAGE_UNKNOWN,True)
		except:
			print(f"couldn't read {imgidx} from {imgsrc}")
			return

		idx2d=[]
		for k,i in enumerate(imgidx): 
			e=rhdrs[k]
			dc={"src":imgsrc,"idx":i,
				"idx3d":ii, "xform.projection":e["xform.projection"], "tilt_id":e["tilt_id"]}
			idx2d.append(len(info2d))
			info2d.append(dc)

		dc={"src":dt["src"], "idx":dt["idx"],
			"coord":img["ptcl_source_coord"], "idx2d":idx2d,
			"class":img["orig_class"], "xform.align3d":img["xform.align3d"]
		}

		info3d.append(dc)

		sys.stdout.write("\r {}/{}".format(ii+1, len(data)))
		sys.stdout.flush()
	print()
	return info2d, info3d

def main():
	
	usage=" "
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--ptcls", type=str,help="input 3d particles from sets/", default=None)
	parser.add_argument("--ali2d", type=str,help="existing aligned 2d particles from spt_xx, used to interpolate local subtilt translation", default=None)
	parser.add_argument("--ali3d", type=str,help="use align3d from a list in place of the info from --ptcls header. require same particles as --ptcls", default=None)
	parser.add_argument("--path", type=str,help="new refinement path", default=None)
	#parser.add_argument("--sym", type=str,help="", default="c1")
	parser.add_argument("--userot", action="store_true", default=False, help="use rotational subtilt alignment as well. very slow and may not be as useful..")

	(options, args) = parser.parse_args()
	
	logid=E2init(sys.argv)
	
	if options.path==None: options.path=num_path_new("spt_")
	print("Writing in {}".format(options.path))
	
	print("Gather metadata from {}".format(options.ptcls))
	info2d1, info3d1 = gather_metadata(options.ptcls)
	print("Total {} 3d particles and {} 2d particles".format(len(info3d1), len(info2d1)))
	if options.ali3d:
		a3=load_lst_params(options.ali3d)
		if len(a3)==len(info3d1):
			for i,a in zip(info3d1, a3):
				i["xform.align3d"]=a["xform.align3d"]
		else:
			print("error. --ali3d and --ptcls do not match")
	
	fali2d=options.ali2d
	path0=os.path.dirname(fali2d)
	p=fali2d.rfind('.')
	itr=int(fali2d[p-2:p])
	fali3d=[f"{path0}/aliptcls3d_{i:02d}.lst" for i in range(itr)]
	fali3d=[f for f in fali3d if os.path.isfile(f)]
	fali3d=fali3d[-1]
	finfo3=f"{path0}/particle_info_3d.lst"
	print("Reading previous alignment")
	print(f"Loading particle info from {finfo3}")
	print(f"        3d alignment from {fali3d}")
	print(f"        2d alignment from {fali2d}")
	info3d0=load_lst_params(finfo3)
	ali2d0=load_lst_params(fali2d)
	ali3d0=load_lst_params(fali3d)
	e0=EMData(info3d0[0]["src"], info3d0[0]["idx"], True)
	e1=EMData(info3d1[0]["src"], info3d1[0]["idx"], True)
	scale=e1["apix_x"]/e0["apix_x"]
	print(f"Scale from old to new alignment {scale:.3f}")

	print("Sanity check on one tomogram")
	fname=base_name(info3d0[0]["src"])
	i3d0=[i for i in info3d0 if base_name(i["src"])==fname]
	i3d1=[i for i in info3d1 if base_name(i["src"])==fname]
	pos0=np.array([i["coord"] for i in i3d0])
	pos1=np.array([i["coord"] for i in i3d1])
	pos1*=scale
	tid=0
	print(f"tomogram {fname}, tilt id {tid}")
	a2d0=[i for i in ali2d0 if base_name(i["src"])==fname]
	v0=np.array([a["dxf"].get_trans() for a in a2d0 if a["tilt_id"]==tid])
	dst=scipydist.cdist(pos0,pos1)
	print(f"{len(pos0)} particles from old refinement, {len(pos1)} particles from new refinement")
	print("mean neighbor distance between particles from two refinement: {:.4f}".format(np.mean(np.min(dst, axis=0))))
	dst=np.exp(-dst/10)
	dst=dst/np.sum(dst,axis=0)
	v1=np.matmul(v0.T,dst).T
	print("mean subtilt translation {:.4f} from old refinement".format(np.mean(np.linalg.norm(v0, axis=1))))
	print("mean subtilt translation {:.4f} from new refinement".format(np.mean(np.linalg.norm(v1, axis=1))))
	
	ali3d1=[]
	for a in info3d1:
		dc={"src":a["src"],"idx":a["idx"],"xform.align3d":a["xform.align3d"]}
		if "class" in a:
			dc["class"]=a["class"]
			
		ali3d1.append(dc)
			
	ali2d1=[]
	for a in info2d1:
		
		dc2={"src":a["src"],"idx":a["idx"],
			"ptcl3d_id":a["idx3d"],
			"tilt_id":a["tilt_id"]}
		
		i3=info3d1[a["idx3d"]]
		if "class" in i3:
			dc2["class"]=i3["class"]
		else:
			dc2["class"]=a["idx3d"]%2
			
		ali2d1.append(dc2)
		
		
	nnb=15
	def test_rot(x):
		xfr=Transform({"type":"xyz","xtilt":x[0],"ytilt":x[1],"ztilt":x[2]}).inverse()
		drspin=[(a*xfr).get_params("spin")["omega"] for a in dr]
		scr=np.dot(drspin, ds)
		return scr
	
	print("Compute subtilt translation for new refinement")
	fnames=np.unique([i["src"] for i in info3d1])
	fnames=[base_name(f) for f in fnames]
	for fname in fnames:
		# i2d1=[i for i in info2d1 if base_name(i["src"])==fname]
		i3d1=[i for i in info3d1 if base_name(i["src"])==fname]
		i3d0=[i for i in info3d0 if base_name(i["src"])==fname]
		
		a2d0=[[ali2d0[i] for i in i3["idx2d"]] for i3 in i3d0]
		i2d1=[[info2d1[i] for i in i3["idx2d"]] for i3 in i3d1]
		a2d1=[[ali2d1[i] for i in i3["idx2d"]] for i3 in i3d1]
		print(fname, len(i3d0), len(i3d1))
		
		txfs=[d["xform.align3d"].inverse() for d in i3d1]
		
		pos0=np.array([i["coord"] for i in i3d0])
		pos1=np.array([i["coord"] for i in i3d1])*scale
		tids=[a["tilt_id"] for a1 in i2d1 for a in a1]
		for tid in sorted(np.unique(tids)):
			ad0=[a for a2 in a2d0 for a in a2 if a["tilt_id"]==tid]
			if1=[a for a1 in i2d1 for a in a1 if a["tilt_id"]==tid]
			ad1=[a for a1 in a2d1 for a in a1 if a["tilt_id"]==tid]
			v0=np.array([a["dxf"].get_trans() for a in ad0])
			i0=[tid in [a["tilt_id"] for a in a2] for a2 in a2d0 ]
			p0=pos0[i0]
			i1=[tid in [a["tilt_id"] for a in a1] for a1 in a2d1 ]
			p1=pos1[i1]
			
			dst=scipydist.cdist(p0,p1)
			dst=np.exp(-dst/10)
			dst=dst/np.sum(dst,axis=0)
			v1=np.matmul(v0.T,dst).T
			v1/=scale
			
			scr=np.array([a["score"] for a in ad0])
			scr=np.dot(scr, dst)
        
			if options.userot:
				rts=[]
				for ip in range(len(ad1)):
					nid=np.argsort(-dst[:,ip])[:nnb]
					ds=dst[nid, ip].copy()
					ds/=np.sum(ds)
					dr=[Transform(ad0[i]["dxf"].get_rotation()) for i in nid]
					x0=[0,0,0]
					res=minimize(test_rot, x0,
							method='Powell',options={'ftol': 1e-3, 'disp': False, "maxiter":5})
					rts.append(res.x)
				rts=np.array(rts)
			else:
				rts=np.zeros((len(ad1),  3))
        
			xfpj=[a["xform.projection"] for a in if1]
			xfraw=[a*b for a,b in zip(xfpj, txfs)]

			dxf=[Transform({"type":"xyz","xtilt":x[0],"ytilt":x[1],"ztilt":x[2],"tx":v[0], "ty":v[1]}) for x,v in zip(rts,v1)]
			xfali=[d*x for x, d in zip(xfraw, dxf)]
			
			for i,a in enumerate(ad1):
				a["dxf"]=dxf[i]
				a["xform.projection"]=xfali[i]
				a["score"]=float(scr[i])
				
				
	save_lst_params(ali3d1, "{}/aliptcls3d_01.lst".format(options.path))
	save_lst_params(ali2d1, "{}/aliptcls2d_01.lst".format(options.path))
	save_lst_params(info2d1, "{}/particle_info_2d.lst".format(options.path))
	save_lst_params(info3d1, "{}/particle_info_3d.lst".format(options.path))
	
	#if options.lstin:
		#lst0=load_lst_params(options.lstin)
		#key0=["{}_{:04d}".format(l["src"], l["idx"]) for l in lst0]
		#dic0={k:i for i,k in enumerate(key0)}
	#lst1=load_lst_params(options.ptcls)
	
	#lout3d=[]
	#lout2d=[]
	#info3d=[]
	#info2d=[]
	"""
	ndxf=0
	for li,l in enumerate(lst1):
		e=EMData(l["src"], l["idx"], True)
		ky="{}_{:04d}".format(e["orig_ptcl"], e["orig_idx"])
		
		if options.lstin:
			ii=dic0[ky]
			l0=lst0[ii]
			scr=l0["score"]
			if "orig_idx" in l0:
				ii=l0["orig_idx"]
		else:
			scr=1
			ii=li
			
		if "xform.align3d" in l:
			xf3d=l["xform.align3d"]
		else:
			xf3d=e["xform.align3d"]

		dc3={"src":l["src"], "idx":l["idx"], 
			"score":scr, "xform.align3d":xf3d,
			"orig_idx":ii}
		lout3d.append(dc3)

		imgsrc=e["class_ptcl_src"]
		imgidx=e["class_ptcl_idxs"]
		coord=e["ptcl_source_coord"]
		boxsz=e["nx"]
		
		rhdrs=EMData.read_images(imgsrc,imgidx,IMAGE_UNKNOWN,True)
		idx2d=[]
		for k,i in enumerate(imgidx): 
			h=rhdrs[k]
			xf=Transform(h["xform.projection"])
			xf=xf*(xf3d.inverse())
			dc2={"src":imgsrc,"idx":i, "class":ii%2,
				"ptcl3d_id":li, "xform.projection":xf,
				"tilt_id":h["tilt_id"]}
			if h.has_attr("orig_dxf"):
				dc2["dxf"]=h["orig_dxf"]
				ndxf+=1
			else:
				dc2["dxf"]=Transform()
				
			dc={"src":imgsrc,"idx":i,
				"idx3d":li,
				"xform.projection":Transform(h["xform.projection"]), "tilt_id":h["tilt_id"]}
			
			idx2d.append(len(lout2d))
			lout2d.append(dc2)
			info2d.append(dc)
				
		dc={"src":l["src"], "idx":l["idx"],
			"coord":coord, "idx2d":idx2d, "xform.align3d":xf3d}
		info3d.append(dc)
		
		sys.stdout.write("\r {}/{}".format(li+1, len(lst1)))
		sys.stdout.flush()
		
	if options.make3d:
		for eo in ["even", "odd"]:
			run(f"e2spa_make3d.py --input {options.path}/aliptcls2d_01.lst --output {options.path}/threed_00_{eo}.hdf --clsid {eo} --parallel thread:32")
		
	"""
	
	E2end(logid)

	
if __name__ == '__main__':
	main()
	
