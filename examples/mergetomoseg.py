#!/usr/bin/env python
# Muyuan Chen 2017-06
from EMAN2 import *
import numpy as np

def main():
	
	usage=""" This program merges multiple tomogram segmentation outputs into a multilevel mask file, where each type of feature is labeled as an interger value. Note that this program reads all maps into memory so it can be quite ineffcient when the map is large.
	
	[prog] seg1.hdf seg2.hdf seg3.hdf --output mask.hdf
	
	"""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--output", type=str,help="Output mask file.", default=None)
	parser.add_argument("--thr", type=float,help="Threshold to define nothingness, assuming the manual annotation  in the given training set to be 1. Default is 0.8", default=.8)
	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)
	
	if not options.output:
		print "Output file name is required. Exit."
		return
	if len(args)<2:
		print "Less than 2 inputs. Nothing to merge. Exit."
		return
	
	imgs=[]
	shp=[]
	lbs=[]
	apix=1
	for ii,a in enumerate(args):
		e=EMData(a)
		
		m=e.numpy().copy()
		if len(shp)==0:
			shp=m.shape
			apix=e["apix_x"]
		else:
			if shp!=m.shape:
				print "Error: Shape of the input files does not match ({} and {}). Exit.".format(shp, m.shape)
				return
		imgs.append(m)
		try:
			lbs.append(e["nnet_src"])
			print "{} : {}".format(ii+1, e["nnet_src"])
		except: 
			lbs.append(a)
		
	
	m0=np.zeros_like(imgs[0])+options.thr
	imgs=[m0]+imgs
	
	mm=np.argmax(np.array(imgs), 0)
	
	e=from_numpy(mm.copy())
	e["apix_x"]=e["apix_y"]=e["apix_z"]=apix
	
	e["labels"]=lbs
	
	e.write_image(options.output)

	E2end(logid)
	
def run(cmd):
	print cmd
	launch_childprocess(cmd)
	
	
if __name__ == '__main__':
	main()
	