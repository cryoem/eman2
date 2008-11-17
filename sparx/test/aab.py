#!/usr/bin/env python

from EMAN2  import *
from sparx  import *
from math import pi, sin, cos
from sys import exit
fc = 0.12
sb = Util.sincBlackman(20, fc)
#exit()
"""
M = 20
fc = 0.12
M2 = M/2
#d = model_blank(256)
#info(d)
q = 0.0
ff = []
for x in xrange(0,M2+1):
	if(x == 0):
		xx = 1.e-7
		f = sin(2*pi*fc*(xx))/(xx)*(0.52-0.5*cos(2*pi*(xx-M2)/M)+0.08*cos(4*pi*(xx-M2)/M))
		#2*pi*fc
	else:
		f = sin(2*pi*fc*(x))/(x)*(0.52-0.5*cos(2*pi*(x-M2)/M)+0.08*cos(4*pi*(x-M2)/M))
	print x,f
	ff.append(f)
for x in xrange(0,M+1):
	if(x == M2):
		xx = x+1.e-7
		f = sin(2*pi*fc*(xx-M2))/(xx-M2)*(0.52-0.5*cos(2*pi*xx/M)+0.08*cos(4*pi*xx/M))
		#2*pi*fc
	else:
		f = sin(2*pi*fc*(x-M2))/(x-M2)*(0.52-0.5*cos(2*pi*x/M)+0.08*cos(4*pi*x/M))
	ff.append(f)
	print x,f
	#d.set_value_at(x,f)
exit()
#print_image(d)
d = power(periodogram(d),0.5)
#print_image(d)

nx = d.get_xsize()
for i in xrange(nx):
	print i,"   ",d.get_value_at(i)

exit()
kernel = model_blank(M+1,M+1)
for y in xrange(0,M+1):
	for x in xrange(0,M+1):
		kernel.set_value_at(x,y,ff[x]*ff[y])

"""
mx = 4096
my = 4096
#mx = 8964
#my = 8964
#mx = 8964
#mx = 10120
#my = 13356
#Wed Aug 13 13:43:25 2008
#Wed Aug 13 13:43:36 2008
#Wed Aug 13 13:44:05 2008

"""
#  assume we want to subsample q times
q = 3.7
# cutoff freq is 0.5/q
cf = 0.5/q
# the gaussian halfwidth in real space is r
r = 1.0/(2*pi*cf)
print  " Gaussian halfwidth :",cf,r
#  It is a good question how wide gaussian kernel should be...
n = 7
kernel = model_gauss(r, n, n)
"""
"""
n = M+1
info(kernel)
print_col(kernel,n//2)
ss = Util.infomask(kernel, None, True)
kernel /=(ss[0]*n*n)
info(kernel)
print_col(kernel,n//2)
"""
#a = model_gauss_noise(1.0,mx,my)
a = model_circle(0.5,mx,my)
 
scale = fc/0.5
print ttime()
#b = rsconvolution(a,kernel)
b = a.downsample(sb, scale)
print ttime()
info(a)
info(b)
mmx = b.get_xsize()
mmy = b.get_ysize()
drop_image(b, "ito.hdf")
p = power(periodogram(b),0.5)*mmx*mmy
drop_image(p, "oto.hdf")
d = []
for i in xrange(mmx):
	d.append(p.get_value_at(i,mmy//2))
del p
write_text_file(d,'toto.txt')
print ttime()

b = filt_tanl(a, fc, 0.05)

print ttime()
exit()
mmx = smallprime(mx,3)
mmy = smallprime(my,3)
a = Util.window(a,mmx,mmy,1,0,0,0)
b = filt_gaussl(a, 0.2)
print ttime()

