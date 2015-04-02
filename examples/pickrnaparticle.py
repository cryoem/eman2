#!/usr/bin/env python
# Muyuan Chen 2015-03-24
# pick rna particles
from EMAN2 import *
import json


def main():
	usage=""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--ub", type=float,help="upper bound", default=.5)
	parser.add_argument("--lb", type=float,help="lower bound", default=.5)
	(options, args) = parser.parse_args()
	
	
	img=EMData(args[0])
	
	print "Generating mask.."
	img.mult(-1)
	img.process_inplace("filter.lowpass.gauss",{"cutoff_abs":0.01})
	img.process_inplace("normalize.edgemean")
	img.process_inplace("filter.ramp")
	img.process_inplace("threshold.belowtozero",{"minval":img["mean_nonzero"]+img["sigma_nonzero"]})
	img.process_inplace("mask.addshells",{"nshells":20})
	img.mult(-1)
	img.process_inplace("normalize.edgemean")
	img.process_inplace("threshold.binaryrange",{"low":1e-10,"high":100})
	img.write_image("masktmp.hdf")
	
	print "Balancing image.."
	img=EMData(args[0])
	img.mult(-1)
	img.sub(img["minimum"])
	#img.process_inplace("filter.ramp")
	#img.process_inplace("normalize.edgemean")
	#img.process_inplace("threshold.belowtozero",{"minval":img["mean"]})
	#img.process_inplace("mask.fromfile",{"filename":"masktmp.hdf"})
	
	img.process_inplace("filter.lowpass.gauss",{"cutoff_abs":0.2})
	img.process_inplace("mask.fromfile",{"filename":"masktmp.hdf"})
	img.sub(img["mean_nonzero"])
	img.process_inplace("mask.fromfile",{"filename":"masktmp.hdf"})
	#img.process_inplace("filter.ramp")
	msk=img.process("filter.lowpass.gauss",{"cutoff_abs":0.005})
	img.sub(msk)
	
	print "Binarizing image.."
	img.process_inplace("threshold.belowtozero",{"minval":img["mean"]})
	img.process_inplace("mask.dust3d",{"voxels":100,"threshold":img["mean_nonzero"]})
	img.process_inplace("threshold.belowtozero",{"minval":img["mean_nonzero"]})
	#img.process_inplace("threshold.belowtozero",{"minval":img["mean_nonzero"]})
	img.process_inplace("mask.dust3d",{"voxels":400,"threshold":0})
	img.process_inplace("morph.object.density",{"thresh":img["mean"]})
	#img.process_inplace("threshold.belowtozero",{"minval":1e-3})
	#img.process_inplace("filter.lowpass.gauss",{"cutoff_abs":0.2})
	#img.process_inplace("threshold.binaryrange",{"low":img["mean"],"high":100})
	#exit()
	img.process_inplace("mask.fromfile",{"filename":"masktmp.hdf"})
	img.process_inplace("threshold.abovetozero",{"maxval":img["mean_nonzero"]+img["sigma_nonzero"]*options.ub})
	img.process_inplace("threshold.belowtozero",{"minval":img["mean_nonzero"]+img["sigma_nonzero"]*options.lb})
	
	
	print "Find centers.."
	imlabel=img.process("morph.object.label",{"thresh":img["mean"],"write_centers":True})
	cnts=imlabel["obj_centers"]
	ct=[]
	for c in cnts:
		ct.append(c)
	print ct
	print cnts
	print imlabel["obj_centers"]
	box=[]
	for i in range(0,len(cnts),3):
		box.append([cnts[i]/1,cnts[i+1]/1,"manual"])
	jn={}
	jn["boxes"]=box
	#print jn
	mn=args[0]
	img.write_image("{n}_tmp.hdf".format(n=mn[:-4]))
	f = open('info/{n}_info.json'.format(n=mn[:-4]), 'w')
	json.dump(jn, f, indent=0)
	f.close()
	
	#img.write_image(sys.argv[2])


if __name__ == '__main__':
	main()
	