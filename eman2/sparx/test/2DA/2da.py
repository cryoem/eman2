#!/usr/bin/env python
from EMAN2  import *
from sparx  import *
from sys import exit
alpha = 40.0
n = 512
#sx = 7.11
#sy =-11.98
a = test_image(size=(n,n))
a = filt_tophatl(a,0.49)
#a = filt_gaussh(a,0.11)
#st = Util.infomask(a, None, True)
#a -=st[0]
#a = model_gauss_noise(1.0,n,n)
#a = get_im("/Users/pawel/ng12_2.hdf")
n = a.get_xsize()
#a = filt_tophatb(a,0.03,0.49)
m = model_circle(n//2-6,n,n)
for alpha in xrange(0,91):
	h = rot_shift2D(rot_shift2D(a,alpha,interpolation_method = "ftgridding"), -alpha, interpolation_method="ftgridding")
	#h = filt_tophatb(h,0.03,0.49)
	q = rot_shift2D(rot_shift2D(a,alpha,interpolation_method = "quadratic"), -alpha, interpolation_method="quadratic")
	g = rot_shift2D(rot_shift2D(a,alpha,interpolation_method = "gridding"), -alpha, interpolation_method="gridding")
	l = rot_shift2D(rot_shift2D(a,alpha,interpolation_method = "linear"), -alpha, interpolation_method="linear")
	f=fsc(a*m,h*m)
	fq = fsc(a*m,q*m)
	fg = fsc(a*m,g*m)
	fl = fsc(a*m,l*m)
	print ccc(a,l,m), ccc(a,q,m), ccc(a,g,m), ccc(a,h,m)
exit()
"""
info(a)
info(g)
"""
d, ap, bp = im_diff(h,a,m)

#info(h)
h = ap*h-bp
#info(h)

dropImage(m*periodogram((a-h)*m),"fg.hdf")
dropImage(m*periodogram((a-g)*m),"gg.hdf")
write_text_file([fl[1], fq[1], fg[1], f[1]], "rfl_40.txt")
write_text_file([rot_avg_table(periodogram((a-g)*m)), rot_avg_table(periodogram((a-h)*m))], "ror.txt")

for i in xrange(10):
	fl = i*0.05
	fh = fl + 0.05
	dag = square(m*filt_tophatb((a-g)*m,fl,fh))
	sg = Util.infomask(dag, None, True)
	dah = square(m*filt_tophatb((a-h)*m,fl,fh))
	sh = Util.infomask(dah, None, True)
	print  i,sg[0], sh[0]
