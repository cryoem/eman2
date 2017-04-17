#!/usr/bin/env python

# Simple script to extract the same pixel from a stack of images. Normally this would be used, eg, with aligned particles from a class-average to look at the distribution of values contributing to the average

from EMAN2 import *

import sys

a=EMData.read_images(sys.argv[1])

x=int(sys.argv[2])
y=int(sys.argv[3])
#out=file("pixel_{:03d}_{:03d}.txt".format(x,y),"w")
outl=file("ccp.txt","w")

rsum=sum(a)
rsum.mult(1.0/len(a))

#	i.process_inplace("filter.highpass.gauss",{"cutoff_freq":0.1})
#	out.write("{}\n".format(i[x,y]))
	
for n,i in enumerate(a):
										   
	ic=i.copy()
	ic.mult(rsum)
	ic.process_inplace("filter.lowpass.gauss",{"cutoff_freq":.02})
	ic.process_inplace("threshold.belowtozero",{"minval":0})
#	ic.process_inplace("math.sqrt")
#	ic.write_image("cfim.hdf",n)
	
	if n<50 :
		icm=i.copy()
		icm.mult(ic)
#		rsumc=rsum.copy()
#		rsumc.mult(ic)
		frc=icm.calc_fourier_shell_correlation(rsum,16.0)
		fout=open("frc_{:03d}.txt".format(n),"w")
		for j in xrange(len(frc)/3):
			fout.write("{}\t{}\n".format(frc[j],frc[j+len(frc)//3]))
			
		icm.write_image("cmcmp.hdf",n*2)
		rsum.write_image("cmcmp.hdf",n*2+1)

print "local"

#for x in xrange(-100,150,50):
	#for y in xrange(-100,150,50):
		#for n,i in enumerate(a):
			#im=i.process("mask.gaussian",{"outer_radius":25,"dx":x,"dy":y})
			#lsum=rsum.process("mask.gaussian",{"inner_radius":10,"outer_radius":25,"dx":x,"dy":y})
			#c=lsum.calc_ccf(im)
			#c.process_inplace("xform.phaseorigin.tocenter")
			#c.write_image("ccl.hdf",-1)
			#l=c.calc_max_location()
			#outl.write("{}\t{}\n".format(l[0],l[1]))

#	c=rsum.calc_ccf(i)
#	c.process_inplace("xform.phaseorigin.tocenter")
#	c.write_image("cc.hdf",n)
#	try: csum.add(c)
#	except: csum=c

#csum.write_image("csum.hdf")
