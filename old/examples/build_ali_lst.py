#!/usr/bin/env python
# Muyuan Chen 2016-05
from EMAN2 import *
import numpy as np
import json

def main():
	
	usage="prog [lst name] --path  "
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--path", type=str,help="path", default="refine_01")
	parser.add_argument("--replace", type=str,help="replace particle set", default=None)
	parser.add_argument("--szmult", type=float,help="size mult", default=None)
	parser.add_argument("--flatten", action="store_true",help="flatten the euler angle distribution", default=False)
	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)
	
	with open('/'.join([options.path,"0_refine_parms.json"])) as json_file:
		js = json.load(json_file)
	
	if not options.szmult:
		if not options.replace:
			options.szmult=1.
		else:
			inp=str(js["input"][0])
			e=EMData(inp, 0, True)
			srcname= e["data_source"]
			e=EMData(srcname, 0, True)
			apix0=e["apix_x"]
			if options.replace:
				tail=srcname.find("__")
				if tail<0:
					print "cannot find particles to replace"
					exit()
				repname= srcname[:tail+2]+options.replace+srcname[-4:]
			e1=EMData(repname, 0, True)
			apix1=e1["apix_x"]
			options.szmult=apix0/apix1
			
	
	eo=["even","odd"]
	for eoid in range(2):
		last=str(js["last_{}".format(eo[eoid])])
		clsmxfile=last.replace("threed","classmx")
		ptcls=str(js["input"][eoid])
		
		lstname=args[0][:-4]+'_'+eo[eoid]+".lst"
		print lstname,clsmxfile, ptcls
		try: os.remove(lstname)
		except: pass
		lst=LSXFile(lstname, False)
		clsfile=clsmxfile.replace("classmx","projections")
		
		clsmx_0=EMData(clsmxfile,0)
		clsmx_ang=EMData(clsmxfile,4)
		clsmx_tx=EMData(clsmxfile,2)
		clsmx_ty=EMData(clsmxfile,3)
		clsmx_mr=EMData(clsmxfile,5)
		clsimgs=EMData.read_images(clsfile,[],True)
		#print clsmxfile
		#print clsmx_0["nx"],clsmx_0["ny"]
		#exit()
		num=clsmx_0["ny"]
		for i in range(num):
			clsid=int(clsmx_0[0,i])
			clsang=clsmx_ang[0,i]
			clstx=clsmx_tx[0,i]
			clsty=clsmx_ty[0,i]
			#print clsid,clsang,clstx
			cls=clsimgs[clsid]
			tr=Transform()
			tr.set_rotation({"type":"2d", "alpha":clsang})
			tr.translate(clstx,clsty)
			tr.set_mirror(clsmx_mr[0,i]>0)

			rot=cls["xform.projection"]
			rr=rot.get_params("eman")
			tr_inv=tr.inverse()
			rt=tr_inv.get_params("eman")
			rr["phi"]=rt["phi"]
			rr["mirror"]=rt["mirror"]
			
			rr["tx"]=rt["tx"]*options.szmult
			rr["ty"]=rt["ty"]*options.szmult
			
			#print Transform(eval(ss))
			e=EMData(ptcls,i,True)
			srcname= e["data_source"]
			if options.replace:
				tail=srcname.find("__")
				if tail<0:
					print "cannot find particles to replace"
					exit()
					
				srcname= srcname[:tail+2]+options.replace+srcname[-4:]
			#exit()
			lst.write(-1,e["data_n"], srcname,str(rr))
			
		
		
		
		lst=None
		
	E2end(logid)
	
def run(cmd):
	print cmd
	launch_childprocess(cmd)
	
	
if __name__ == '__main__':
	main()
	