#!/usr/bin/env python
# Muyuan Chen 2023-12
from EMAN2 import *
import numpy as np

def main():
	
	usage="""prog --ptcls spt_xx/aliptcls3d_xx.lst --tomo <tomogram>
	Map aligned particles back to tomograms. 
	This will generate a script to place the objects in unreal engine 5."""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--ptcls", type=str,help="particle list file", default="")
	parser.add_argument("--tomo", type=str,help="tomogram file name", default="")
	parser.add_argument("--pyout", type=str,help="output py file for ue5", default="from_tomo.py")
	parser.add_argument("--ue_objs", type=str,help="objects to place in ue5. multiple objects separated by comma", default="")
	parser.add_argument("--ue_dir", type=str,help="folder that contains objects to place in ue5. default is /Game/", default="/Game/")
	# parser.add_argument("--postxf", type=str,help="extra shift after alignment", default="")
	parser.add_argument("--ppid", type=int,help="ppid...", default=-1)
	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)
	
	ptcl=[]
	bname=base_name(options.tomo)
	alipm=load_lst_params(options.ptcls)
	path=os.path.dirname(options.ptcls)
	js=js_open_dict(info_name(options.tomo))
	apix_unbin=js["apix_unbin"]
	info3d=load_lst_params("{}/particle_info_3d.lst".format(path))
	for i,p in enumerate(alipm):
		if base_name(p["src"])==bname:
			ptcl.append((p["score"], p["src"], p["idx"], p["xform.align3d"], info3d[i]["coord"]))
	

	nptcl=len(ptcl)	
	print("{:d} particles total.".format(int(nptcl)))
	if nptcl==0:
		print("No particles. Exiting")
		return

	tomo=EMData(options.tomo, 0, True)
	
	pos=np.array([p[-1] for p in ptcl])
	posstr="["
	
	for ii, p in enumerate(pos):
		
		xf=Transform(ptcl[ii][-2])
		crd=p.copy()*apix_unbin #+ [tomo["nx"]//2, tomo["ny"]//2, tomo["nz"]//2]
		ts=np.array(xf.get_trans())*apix_unbin
		xf.set_trans(ts.tolist())
		xf=xf.inverse()
		# xp=xf.transform(0,0,-50)
		
		t2=xf.get_params("xyz")
		t2["ztilt"]-=90.
		t2=Transform(t2)
		x=t2.get_params("quaternion")
		s=crd+[x["tx"], x["ty"], x["tz"]]
		x=[x["e0"], x["e1"], x["e2"], x["e3"]]
		
		
		posstr+='['
		posstr+=", ".join([f"{i:.6f}" for i in s])
		posstr+=', '
		posstr+=", ".join([f"{i:.6f}" for i in x])
		posstr+='],\n'
		
	posstr+=']'
	objs=options.ue_objs.split(',')
	objs=','.join([f'"{o}"' for o in objs])
	
	
	pyfile="""
import unreal
import numpy as np

pos={p}
objlst=[{o}]
uedir="{d}"
""".format(p=posstr, o=objs, d=options.ue_dir)
	
	pyfile+="""
scale=0.1
pos=np.array(pos)
pos[:,:3]=pos[:,[2,0,1]]/scale

pp=np.mean(pos, axis=0)

root0=unreal.EditorLevelLibrary.spawn_actor_from_class(unreal.StaticMeshActor,(pp[0],pp[1],pp[2]))
root0.set_actor_label(f"root_00")

nsym=1
for ip,pp in enumerate(pos):
	v=unreal.Vector(pp[0],pp[1],pp[2])
	q=unreal.Quat(pp[3], pp[4], pp[5], pp[6])
	r1=q.rotator()

	r90=unreal.Rotator(0,0,-90)
	r1=r90.combine(r1)
	root=unreal.EditorLevelLibrary.spawn_actor_from_class(unreal.StaticMeshActor,(pp[0],pp[1],pp[2]), q.rotator())
	root.set_actor_label(f"object_{ip:03d}")
	root.attach_to_actor(root0,"none",unreal.AttachmentRule.KEEP_WORLD,unreal.AttachmentRule.KEEP_WORLD,unreal.AttachmentRule.KEEP_WORLD)
	
	for o in objlst:
		obj=unreal.EditorAssetLibrary.load_asset(uedir+f"{o}.{o}")

		# act=unreal.EditorLevelLibrary.spawn_actor_from_object(obj,v, q.rotator())
		act=unreal.EditorLevelLibrary.spawn_actor_from_object(obj,v, r1)
		act.set_actor_relative_scale3d((scale, scale, scale))
		act.attach_to_actor(root,"none",unreal.AttachmentRule.KEEP_WORLD,unreal.AttachmentRule.KEEP_WORLD,unreal.AttachmentRule.KEEP_WORLD)
		act.set_actor_label(o)

p=-np.mean(pos, axis=0)
c=(0,0,5000)
center=unreal.Vector(p[0]+c[0],p[1]+c[1],p[2]+c[2])
root0.add_actor_world_offset(center, False, True)
root0.add_actor_world_rotation(unreal.Rotator(0,90,0), False, True)
"""
	
	f=open(options.pyout, 'w')
	f.write(pyfile)
	f.close()
	print(f"output written to {options.pyout}")

	
	
if __name__ == '__main__':
	main()
	
