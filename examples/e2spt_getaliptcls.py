#!/usr/bin/env python
from __future__ import division

# Muyuan Chen 2019-05
from EMAN2 import *
import numpy as np

def main():
	
	usage="""get a stack of aligned particles from a spt alignment.  
	e2spt_getaliptcls.py spt_xx/particle_parms_xx.json --output <particle stack>"
	"""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--output", type=str,help="output", default="")
	parser.add_argument("--rand", type=int,help="take N random particles", default=-1)
	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)
	
	jsfile=args[0]
	if options.output=="":
		output=jsfile.replace("particle_parms", "aliptcls")[:-4]+"hdf"
	else:
		output=options.output
	
	print("Writing output to {}...".format(output))
	try: os.remove(output)
	except: pass
	js=js_open_dict(jsfile)
	keys=sorted(list(js.keys()))
	
	if options.rand>0:
		np.random.shuffle(keys)
		keys=keys[:options.rand]
	
	for i,ky in enumerate(keys):
		src, ii=eval(ky)
		e=EMData(src, ii)
		xf=js[ky]["xform.align3d"]
		e.transform(xf)
		e.write_image(output, -1)
		sys.stdout.write("\r  {}/{}".format(i, len(keys)))
		sys.stdout.flush()
	
	js.close()
	print("\nDone")
	E2end(logid)
	
def run(cmd):
	print(cmd)
	launch_childprocess(cmd)
	
	
if __name__ == '__main__':
	main()
	