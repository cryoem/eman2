#!/bin/env python
from EMAN2  import *
from sparx  import *

from sys import exit
from random import seed,gauss
from math import pi


	
def prepij(image):
	M=image.get_xsize()
	npad = 2
	N = M*npad
	K = 6
	alpha = 1.75
	r = M/2
	v = K/2.0/N
	kb = Util.KaiserBessel(alpha, K, r, K/(2.*N), N)
	#out = rotshift2dg(image, angle*pi/180., sx, sy, kb,alpha)
	o = image.FourInterpol(2*M, 2*M, 1, 0)
	params = {"filter_type": Processor.fourier_filter_types.KAISER_SINH_INVERSE,
	          "alpha":alpha, "K":K, "r":r, "v":v, "N":N}
	q = Processor.EMFourierFilter(o, params)
	return  q,kb



"""
nxn=65
nyn=65
for iq in xrange(3,111):
	nx = 2*iq+1
	ny = nx
	print  "   nx =",nx
	e = model_circle(0.5,nx,ny)
	#e = model_gauss(3,15,15)
	[ave,d1,dm,d,nx,ny,nz]=info(e)
	#e /= d
	peak_search(e,print_screen = True)
	o = e.FourInterpol(2*nx, 2*ny, 1, 0)
	#o = e.FourInterpol(nxn, nyn, 1, 0)
	u = fft(o)
	peak_search(u,print_screen = True)
	#drop_image(u,"xxx1.hdf")
	#print_image(o)
	b = pad(cyclic_shift(e,-2,+2),2*nx, 2*ny)
	peak_search(b,print_screen = True)
	cf = ccf(b,o)
	peak_search(cf,print_screen = True)
#print_image(cf)
exit()
print  "  EVEN  "
#e = EMData()
#e.read_image("tf2d0001.tfc")
e = model_circle(0.5,6,6)
[ave,d,d,d,nx,ny,nz]=info(e)
#e -= ave
o = e.FourInterpol(2*nx, 2*ny, 1, 0)
o = fft(o)
peak_search(o,print_screen = True)
print_image(o)
b = fft(pad(cyclic_shift(e,-2,+2),2*nx, 2*ny))
cf = ccf(b,o)
peak_search(cf,print_screen = True)
print_image(cf)

exit()





"""





# padded input
nx = ny = 127
sx = 0.314
sy = -0.273
#sx = 1
#sy = -1
#e = model_circle(0.5,nx,ny)
e = model_gauss(1.,nx,ny)
[ave,d1,dm,d,nx,ny,nz]=info(e)
e /= d
#info(e)
peak_search(e,print_screen = True)

#u,kb=prepi(e)
#es = u.rot_scale_conv(0., 0.114,-0.573, kb)
#print_image(fft(e))
es = e.FourInterpol(2*nx, 2*ny, 1, 0)
#print_image(es)
qe = fshift(e, sx, sy)
qes = qe.FourInterpol(2*nx, 2*ny, 1, 0)
#print_image(qes)
#info(qes)
#es = fshift(e, 0.0,0.0)

u = ccf(qes,es)  # large ccf
info(u, None, " large ccf ")
#print_image(fft(u))
peak_search(u,print_screen = True)
drop_image(u,"pi.hdf")
eg,kb = prepij(e)
eg = fshift(eg,1,1)
u = ccf(qes,eg)
info(u, None, " grid ccf ")
peak_search(u,print_screen = True)
drop_image(u,"pig.hdf")

ns = 200
#cfc = model_blank(2*ns+1,2*ns+1)
xma = -1.0e23
xold = float(nx//2)
yold = float(ny//2)
for i in xrange(-ns,ns+1):
	for j in xrange(-ns,ns+1):
		value = u.get_pixel_conv(xold+i/float(ns),yold+j/float(ns),1.0,kb)
		#cfc.set_value_at(i+ns,j+ns,value)
		if(value > xma):
			xma = value
			im = xold+i/float(ns) - float(nx//2)
			jm = yold+j/float(ns) - float(ny//2)
print  "  PEAK  ",im,jm,xma



exit()





# padded ccf
nx = ny = 7
sx = 0.314
sy = -0.273
#sx = 1
#sy = -1
#e = model_circle(0.5,nx,ny)
e = model_gauss(1.,nx,ny)
[ave,d1,dm,d,nx,ny,nz]=info(e)
e /= d
info(e)
peak_search(e,print_screen = True)

#u,kb=prepi(e)
#es = u.rot_scale_conv(0., sx, sy, kb)
#print_image(fft(e))
#es = e.FourInterpol(2*nx, 2*ny, 1, 0)
#print_image(es)
qes = fshift(e,sx, sy)
#qes = fshift(qes, 0.5,0.5)
#info(qes)
#print_image(fft(qes))
#es = fshift(e, 0.0,0.0)

u = ccf(qes,e)   #small
qu = u.FourInterpol(2*nx, 2*ny, 1, 0)  # large ccf
#print_image(qu)
peak_search(fft(qu),print_screen = True)
qc = fshift(qu,1,1)
#info(qc)
#print_image(qc)
qcf = fft(qc)
#info(qcf)
peak_search(qcf,print_screen = True)
drop_image(qcf,"pc.hdf")







#c = ccf(es,e)
u,kb=prepi(u)  #small ccf
peak_search(u,print_screen = True)

#qe,kb=prepij(e)
#u = ccf(qes,qe)
drop_image(u,"pic.hdf")

ns = 200
#cfc = model_blank(2*ns+1,2*ns+1)
xma = -1.0e23
xold = float(nx//2)
yold = float(ny//2)
for i in xrange(-ns,ns+1):
	for j in xrange(-ns,ns+1):
		value = u.get_pixel_conv(xold+i/float(ns),yold+j/float(ns),1.0,kb)
		#cfc.set_value_at(i+ns,j+ns,value)
		if(value > xma):
			xma = value
			im = xold+i/float(ns) - float(nx//2)
			jm = yold+j/float(ns) - float(ny//2)
print  "  PEAK  ",im,jm,xma
#drop_image(cfc,"cfc.hdf")



exit()



u,kb=prepi(e)
mask=model_circle(30,nx,ny)

sx = -7.34782
sy = 3.78236
#sx = -7.0
#sy = 3.0
angle = 0.0

gs = u.rot_scale_conv(angle*pi/180., sx, sy, kb)
peak_search(e,print_screen = True)
peak_search(gs,print_screen = True)
cf = ccfnp(gs,e)
drop_image(cf,"cf.hdf")
peak_search(cf,print_screen = True)

u,kb=prepi(gs)
res = u.rot_scale_conv(-angle*pi/180., -sx, -sy, kb)
print ccc(e,res,mask)

fs = fshift(e,sx,sy)
ff = ccfnp(fs,e)
drop_image(ff,"ff.hdf")
peak_search(ff,print_screen = True)

print ccc(e,fshift(fs,-sx,-sy),mask)




refi,kb  = prepij(e)
#refi = e.FourInterpol(2*nx, 2*ny, 1, 0)
o = fs.FourInterpol(2*nx, 2*ny, 1, 0)
#info(refi)
print  "  padded "
qt = fft(o)
drop_image(qt,"qt.hdf")
a=peak_search(qt,print_screen = True)



# Calculate cross-correlation fucntion in Fourier space
product = ccf(o,refi)
info(product)
drop_image(product,"product.hdf")

ccfg = product.rot_scale_conv(0., 0., 0., kb)
info(ccfg)

a=peak_search(ccfg,print_screen = True)
print " on gridding ",a
drop_image(ccfg,"ccfg.hdf")

xma = -1.0e23
xold = float(nx//2)+int(a[0][4])
yold = float(ny//2)+int(a[0][5])
ns = 3
for i in xrange(-ns,ns+1):
	for j in xrange(-ns,ns+1):
		value = product.get_pixel_conv(xold+i,yold+j,1.0,kb)
		#print  i,j,"   ",value
		if(value > xma):
			xma = value
			im = xold+i - float(nx//2)
			jm = yold+j - float(ny//2)
print  "  PEAK  ",im,jm,xma
ns = 100
cfc = model_blank(2*ns+1,2*ns+1)
xma = -1.0e23
for i in xrange(-ns,ns+1):
	for j in xrange(-ns,ns+1):
		value = product.get_pixel_conv(xold+i/float(ns),yold+j/float(ns),1.0,kb)
		cfc.set_value_at(i+ns,j+ns,value)
		#print  i,j,"   ",value
		if(value > xma):
			xma = value
			im = xold+i/float(ns) - float(nx//2)
			jm = yold+j/float(ns) - float(ny//2)
print  "  PEAK  ",im,jm,xma
drop_image(cfc,"cfc.hdf")
