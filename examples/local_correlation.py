#!/usr/bin/env python
# This program will compute the 'local correlation' between two maps
# maps are hardcoded as a,b below, as is the local region
from EMAN2 import *

# read the images
a=EMData("bdb:refine_03#threed_filt_00")
b=EMData("bdb:refine_03#threed_filt_01")

# Used for normalization
asq=a.process("math.squared")
bsq=b.process("math.squared")
a.mult(b)

# These define the integration region as a Gaussian weight around each point with a 'resolution' of ~33 A
asq.process_inplace("filter.lowpass.gauss",{"cutoff_freq":.03})
bsq.process_inplace("filter.lowpass.gauss",{"cutoff_freq":.03})
a.process_inplace("filter.lowpass.gauss",{"cutoff_freq":.03})

# normalize
asq.process_inplace("math.sqrt")
bsq.process_inplace("math.sqrt")
asq.process_inplace("math.reciprocal")
bsq.process_inplace("math.reciprocal")
a.mult(asq)
a.mult(bsq)

# masking, to eliminate spurious artifacts at the edges where there is almost no density
b.process_inplace("threshold.binary",{"value":b["mean"]+b["sigma"]/4})
a.mult(b)	# this is a mask of sorts

a.write_image("loc_cor.hdf")
