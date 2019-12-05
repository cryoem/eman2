#!/usr/bin/env python
# Muyuan Chen 2016-10
from EMAN2 import *
import numpy as np

def main():
	
	usage=" "
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--output", type=str,help="output file name", default="cat.hdf")
	parser.add_argument("--boxsize", type=int,help="box size", default=64)
	parser.add_argument("--apix", type=float,help="a/pix", default=1)
	parser.add_argument("--pos", type=float,help="pos of motion", default=0)
	parser.add_argument("--size", type=float,help="relative size", default=.75)
	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)
	
	sz0=options.boxsize
	e=EMData(sz0,sz0,sz0)
	e.to_zero()
	sz=sz0*options.size
	#### body
	body_a=sz*.5*.45
	body_b=sz*.5*.3
	body_c=sz*.5*.3
	a=e.copy()
	a.process_inplace("testimage.ellipsoid",{"a":body_a, "b":body_b, "c":body_c})
	a.process_inplace("mask.addshells.gauss",{"val1":0, "val2":2})
	e.add(a)
	
	#### head
	head_a=sz*.5*.3
	head_b=sz*.5*.3
	head_c=sz*.5*.3
	a=e.copy()
	a.to_zero()
	a.process_inplace("testimage.ellipsoid",{"a":head_a, "b":head_b, "c":head_c})
	a.process_inplace("mask.addshells.gauss",{"val1":0, "val2":2})
	a.translate(body_a,0,body_c+head_c/2.)
	e.add(a)
	
	
	#### ears
	ear_a=sz*.5*.12
	ear_b=sz*.5*.13
	ear_c=sz*.5*.13
	a=e.copy()
	a.to_zero()
	a.process_inplace("testimage.ellipsoid",{"a":ear_a, "b":ear_b, "c":ear_c})
	a.process_inplace("mask.addshells.gauss",{"val1":0, "val2":1})
	a.process_inplace("mask.zeroedge3d",{"x1":sz0*.5-1})
	b=a.copy()
	a.translate(body_a+head_a*.4,head_b,body_c+head_c+ear_c/2.)
	b.translate(body_a+head_a*.4,-head_b,body_c+head_c+ear_c/2.)
	e.add(a)
	e.add(b)
	
	#### legs
	leg_a=sz*.5*.05
	leg_b=sz*.5*.05
	leg_c=sz*.5*.3
	a=e.copy()
	a.to_zero()
	a.process_inplace("testimage.ellipsoid",{"a":leg_a, "b":leg_b, "c":leg_c})
	a.process_inplace("mask.addshells.gauss",{"val1":0, "val2":1})
	for i in [-1,1]:
		for j in [-1,1]:
			b=a.copy()
			b.translate(i*body_a*.5, j*body_b*.5, -body_c)
			e.add(b)
	
	#### tail
	tailin=sz*.5*.2
	tailout=sz*.5*.3
	tailthk=sz0*.5*(1-.05)
	tailarc=sz0*.5*1
	a=e.copy()
	a.to_one()
	a.process_inplace("mask.sharp",{"inner_radius":tailin, "outer_radius":tailout})
	a.process_inplace("mask.zeroedge3d",{"y1":0,"y0":tailarc,"z0":tailthk, "z1":tailthk})
	a.process_inplace("mask.addshells.gauss",{"val1":0, "val2":1})
	a.rotate(0,180.*options.pos,0)
	a.translate(-body_a-tailin,0,0)
	e.add(a)
	
	
	####
	e.process_inplace("threshold.clampminmax",{"maxval":1, "minval":0})
	e.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.3})
	e.process_inplace("filter.highpass.gauss",{"cutoff_abs":.05})
	e["apix_x"]=e["apix_y"]=e["apix_z"]=options.apix
	e.process_inplace("normalize")
	#e.process_inplace("normalize.bymass",{"mass":500})
	e.write_image(options.output)
	
	#display(e)
	E2end(logid)
	
def run(cmd):
	print(cmd)
	launch_childprocess(cmd)
	
	
if __name__ == '__main__':
	main()
	