#!/usr/bin/env python

from EMAN2 import *
from numpy import *

def wrt(ary,s,fsp):
	f=file(fsp,"w")
	for j,i in enumerate(ary): f.write("%f\t%g\n"%(s[j],i))

ctf=EMAN2Ctf()
ctf.voltage=300.0
ctf.cs=4.1
ctf.ampcont=10.0
ctf.apix=2.0
ctf.bfactor=100.0

ds=1.0/(ctf.apix*2*256)
s=arange(0,ds*256,ds)
ctf.dsbg=ds
ctf.background=[.1]*364

ctf.defocus=1.0
curve=array(ctf.compute_1d(len(s)*2,ds,Ctf.CtfType.CTF_AMP))
wrt(curve,s,"c1.txt")

ctf.defocus=1.0
curve=array(ctf.compute_1d(len(s)*2,ds,Ctf.CtfType.CTF_SIGN))
wrt(curve,s,"c2.txt")

ctf.defocus=1.0
curve=array(ctf.compute_1d(len(s)*2,ds,Ctf.CtfType.CTF_TOTAL))
wrt(curve,s,"c3.txt")

img=EMData(258,256)
img.to_zero()
img.set_complex(1)

ctf.compute_2d_complex(img,Ctf.CtfType.CTF_AMP,None)
img.write_image("c1.hdf")

ctf.compute_2d_complex(img,Ctf.CtfType.CTF_SIGN,None)
img.write_image("c2.hdf")

ctf.compute_2d_complex(img,Ctf.CtfType.CTF_TOTAL,None)
img.write_image("c3.hdf")

