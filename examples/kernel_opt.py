# a program to generate a good kernel for gridding-based reconstruction
# This is configured to optimize the kernel with a 256^3 box size.
# the space is 3x oversampled with a 768^3 box, which gives us some interpolated values between points,
# which helps smooth antialiasing. That is, we oversample, then force the data outside the 256^3 volume
# to be as nearly zero as possible, while retaining significant amplitudes within as much of the 256^3 volume as possible
# The process iterates between real and Fourier spaces, forcing appropriate values to zero in each round
from __future__ import print_function
from __future__ import division
from past.utils import old_div
from builtins import range
from EMAN2 import *

targsize=7	# goal is to have an nxn kernel with 3x oversampling
ksize=old_div((targsize*3-1),2)

kstore=EMData(ksize+1,ksize*2+1,ksize*2+1)

map=EMData(768,768,768)
map.to_one()
map.process_inplace("mask.gaussian",{"inner_radius":72,"outer_radius":72})

# we iterate 10 times making sure the things we want to be zero are zero
for i in range(64):
#	map.process_inplace("mask.gaussian",{"inner_radius":104,"outer_radius":8})
#	map.process_inplace("mask.gaussian",{"inner_radius":104,"outer_radius":2})
	map.process_inplace("mask.zeroedge3d",{"x0":266,"x1":266,"y0":266,"y1":266,"z0":266,"z1":266})
	map.process_inplace("xform.phaseorigin.tocorner")
	fft=map.do_fft()

	# extract the kernel we just constructed, also display it
	for x in range(ksize+1):
		for y in range(-ksize,ksize+1):
			for z in range(-ksize,ksize+1):
				kstore[x,y+ksize,z+ksize]=fft.get_complex_at(x,y,z).real
				print("{:7.4f}".format(old_div(kstore[x,y+ksize,z+ksize],(256*256*256))),end="")
			if (abs(y)<=ksize and x<=ksize) : print("")
		if (x<=ksize) : print("\n")
	print("------")

	# clear the existing FFT
	fft.to_zero()
	
	# replace the real kernel elements
	for x in range(ksize+1):
		for y in range(-ksize,ksize+1):
			for z in range(-ksize,ksize+1):
				fft.set_complex_at(x,y,z,kstore[x,y+ksize,z+ksize])


	map=fft.do_ift()
	map.process_inplace("xform.phaseorigin.tocenter")
#	display(map,True)

map.write_image("final_real.hdf")
kstore.write_image("final_kernel.hdf")

# note that this extends 1 value beyond the edge of the array, this last value is a 0 to deal with roundoff errors when using the kernel
print("const float FourierInserter3DMode7::kernel[9][9][9] = {")
for z in range(ksize+2):
	for y in range(ksize+2):
		for x in range(ksize+2):
			print("{:7.7f},".format(old_div(kstore[x,y+ksize,z+ksize],(256*256*256))),end="")
		print("")
		
print("}")
#display(map,True)

# make an oversampled map to look for artifacts from neighboring unit cells
#fft2=EMData(256*3,256*3,256*3).do_fft()
#fft2.to_zero()
#for x in range(3):
#	for y in range(-2,3):
#		for z in range(-2,3):
#			v=fft.get_complex_at(x,y,z)
#			for xx in range(x*3-1,x*3+2):
#				for yy in range(y*3-1,y*3+2):
#					for zz in range(z*3-1,z*3+2):
#						fft2.set_complex_at(xx,yy,zz,v)
#
#
#map2=fft2.do_ift()
#map2.process_inplace("xform.phaseorigin.tocenter")
#display(map2,True)
