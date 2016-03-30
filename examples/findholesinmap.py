#!/usr/bin/env python
# Muyuan Chen 2016-03
# Look for holes in a 3d density map
from EMAN2 import *
import numpy as np
import scipy.ndimage as ndimage

def main():
	
	usage="Look for holes in a 3d density map.\n\tfindholesinmap.py map.mrc\n**********scipy is required!*********"
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--thr", type=float,help="Threshold for the isosurface", default=1)
	parser.add_argument("--closeiter", type=int,help="Number of iterations for the closing operation", default=10)
	parser.add_argument("--filter_res", type=float,help="Resolution for the final filter", default=10)
	parser.add_argument("--output", type=str,help="output file name", default=None)
	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)
	
	e=EMData(args[0])
	img=e.numpy()
	if options.output==None:
		options.output=args[0][:-4]+"_holes.hdf"
	
	apix=e["apix_x"]
	img_open=ndimage.binary_closing(img>options.thr,iterations=options.closeiter)
	m=img.copy()
	m[m<0]=0
	m/=np.max(m)
	hole=img_open-m
	a=from_numpy(hole)

	
	a["apix_x"]=apix
	a["apix_y"]=apix
	a["apix_z"]=apix
	a.process_inplace("filter.lowpass.gauss",{"cutoff_freq":1./options.filter_res})

	a.write_image(options.output)
	
	E2end(logid)
	
def run(cmd):
	print cmd
	launch_childprocess(cmd)
	
	
if __name__ == '__main__':
	main()
	