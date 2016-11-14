#!/bin/env python

from EMAN2  import *
from sparx  import *

n = 128
a = threshold(test_image(size=(n,n))) # a= model_circle(32,n,n)
n*=2
a = pad(a,n,n, background=0.0)
info(a)
ac = fft(a)
piw = periodogram(ac)
mask = binarize(filt_gaussl(a,0.2),0.9)#mask = model_circle(40,n,n)
av= model_blank(n,n)
fsc(a,fshift(a,-0.5,0.5), 1.0, "fsc0")
#exit()


trials = 10
for trt in xrange(trials):

	im = model_gauss_noise(10.0,n,n)

	im = im.replace_amplitudes(ac)

	ci = -1.0
	for i in xrange(2000):
		imn = im.replace_amplitudes(ac)
		imn = threshold(imn*mask)

		#if(i%500 == 0):  
		cc = ccc(piw,periodogram(imn))
		if(abs(cc-ci)<1.0e-7):  break
		ci = cc
		#print  trt,i,"   ",ccc(imn,im),"   ",ccc(imn,a),"   ",cc
		im = imn.copy()
	print  i,trt,"   ",ccc(imn,im),"   ",ccc(im,a),"   ",cc
	fsc(piw,periodogram(imn),1.0,"dfsc%05d"%trt)
	if(trt == 0): refi = im.copy()
	else:
		p = peak_search(ccf(im,refi))
		print  trt,p
		im = fshift(im,-p[0][4],-p[0][5])
		print  trt," after shift  ",ccc(im,a)
	av += im
	im.write_image("o4.hdf",trt)
	imd = im-a
	info(imd)
	imd.write_image("o5.hdf",trt)


av /= trials
print  "   ",ccc(av,a)
info(av)
im = av.copy()
p = peak_search(ccf(im,a))
print  p
im = fshift(im,-p[0][4],-p[0][5])
print  " after shift  ",ccc(im,a)
f = fsc(im,a,1.0,"fsc1")

im.write_image("o6.hdf",0)
imd = im-a
info(imd)
imd.write_image("o6.hdf",1)

for i in xrange(5000):
	imn = im.replace_amplitudes(ac)
	imn = threshold(imn*mask)

	if(i%1000 == 0):  print  i,"   ",ccc(imn,im),"   ",ccc(imn,a)
	im = imn.copy()
print  " before shift  ",ccc(im,a)
p = peak_search(ccf(im,a))
print  p
im = fshift(im,-p[0][4],-p[0][5])
print  " after shift  ",ccc(im,a)
f = fsc(im,a,1.0,"fsc2")
im.write_image("o6.hdf",2)
imd = im-a
info(imd)
imd.write_image("o6.hdf",3)
a.write_image("o6.hdf",4)
mask.write_image("o6.hdf",5)
