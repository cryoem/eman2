#!/usr/bin/env python
# Muyuan Chen 2021-03
from EMAN2 import *
import numpy as np

def main():
	
	usage="convert the json and lst file produced by two versions of e2spt pipeline back and forth. run [prog] <input> <output>"
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)
	
	if args[0].endswith(".json"):
		jsname=args[0]
		js=dict(js_open_dict(jsname)).copy()
		keys=list(js.keys())
		src=eval(keys[0])[0]
		info3d=[]
		n=EMUtil.get_image_count(src)

		lst=LSXFile(src)
		for ii in range(n):
			l=lst.read(ii)
				
			k=str((src, ii))
			j=js[k]
			dc={"idx":l[0],"src":l[1], "class":ii%2,
				"xform.align3d":j["xform.align3d"], "score":j["score"]}
			
			info3d.append(dc)

			sys.stdout.write("\r {}/{}".format(ii, n))
			sys.stdout.flush()


		fm3d=args[1]
		save_lst_params(info3d, fm3d)
	
	elif args[0].endswith(".lst"):
		lstname=args[0]
		info3d=load_lst_params(lstname)
		dic={}
		for i,a in enumerate(info3d):
			k=str((lstname, i))
			d={"xform.align3d":a["xform.align3d"], "score":a["score"]}
			dic[k]=d
		js=js_open_dict(args[1])
		js.update(dic)
		js.close()

	E2end(logid)
	
def run(cmd):
	print(cmd)
	launch_childprocess(cmd)
	
	
if __name__ == '__main__':
	main()
	