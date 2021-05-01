#!/usr/bin/env python
# Muyuan Chen 2021-03
from EMAN2 import *
import numpy as np

def main():
	
	usage=" "
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)
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

	E2end(logid)
	
def run(cmd):
	print(cmd)
	launch_childprocess(cmd)
	
	
if __name__ == '__main__':
	main()
	