#!/usr/bin/env python
# Muyuan Chen 2019-07
from __future__ import print_function
from __future__ import division
from past.utils import old_div
from future import standard_library
standard_library.install_aliases()
from builtins import range
from EMAN2 import *
import numpy as np

def main():
	
	usage="read alignment info from particle header and generate alignment json file. \n [prog] <ptcl input> <json output> "
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--randphi", action="store_true",help="randomize phi", default=False)
	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)
	
	inp=args[0]
	out=args[1]
	if os.path.isfile(out):
		os.remove(out)
	dic={}
	n=EMUtil.get_image_count(inp)
	for i in range(n):
		e=EMData(inp, i, True)
		xf=e["xform.align3d"]
		if options.randphi:
			pm=xf.get_params("eman")
			pm["phi"]=np.random.rand()*360
			xf=Transform(pm)

		d={"xform.align3d":xf, "score":-1}
		dic[(inp, i)]=d
		
	js=js_open_dict(out)
	js.update(dic)
	js.close()
	

	E2end(logid)
	
def run(cmd):
	print(cmd)
	launch_childprocess(cmd)
	
	
if __name__ == '__main__':
	main()
	