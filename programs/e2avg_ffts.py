#!/usr/bin/env python

#
# Author: Mike Schmid, 02/2006  
# Copyright (c) 2000-2006 Baylor College of Medicine
#
# This software is issued under a joint BSD/GNU license. You may use the
# source code in this file under either license. However, note that the
# complete EMAN2 and SPARX software packages have some GPL dependencies,
# so you are responsible for compliance with the licenses of these packages
# if you opt to use BSD licensing. The warranty disclaimer below holds
# in either instance.
#
# This complete copyright notice must be included in any revised version of the
# source code. Additional authorship citations may be added, but existing
# author citations must be preserved.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  2111-1307 USA
#
#

# gives straight, masked weight, and fractional weighted averages.
from EMAN2 import *
import sys
import math

print "WARNING: This program does not currently conform to EMAN2 usage standards"

args=sys.argv

def avg_this_pixel (filename, imagedict, thr, i,j,k):
#    print i
    sumrst=0.
    sumist=0.
    sumr2=0.
    sumi2=0.
    sumrth=0.
    sumith=0.
    suma=0.
#    suma2=0.
    num=0
    numth=0
    for filename in imagedict.keys():
#        print filename, i,j,k
        re= imagedict[filename].get_value_at(i,j,k)
        im= imagedict[filename].get_value_at(i+1,j,k)
#        print i,j,k, re,im
        num=num+1
        amp2=re*re + im*im
        amp=math.sqrt(amp2)
        suma=suma+amp
#        suma2=suma2+amp2
        sumrst=sumrst+re
        sumist=sumist+im
        sumr2=sumr2 + (re*amp)
        sumi2=sumi2 + (im*amp)
        if (amp>float(thr)):
#            print i,j,k,amp
            numth=numth+1
            sumrth=sumrth+re
            sumith=sumith+im
    rst=sumrst/num
    ist=sumist/num
    rth=0.
    ith=0.
    if (numth>0):
        rth=sumrth/numth
        ith=sumith/numth
    rwt=0.
    iwt=0.
    if (suma>0):
        rwt=sumr2/suma
        iwt=sumi2/suma
#        process (filename, imagedict[filename] )
#        print filename, thr, i,j,k
    return (rst,ist,rth,ith,rwt,iwt)

imagedict={}

if (len(args)<4) :
    print "Usage:\navg_ffts.py <out> thr <in1> <in2>...\n"
    print "All files same size (not checked), at least 2 files to average"
    print "Produces 3 files, "
    print "st-out, straight average, should be comparable to real space avg; "
    print "ma-out, amps < thr ignored in number to avg; "
    print "wt-out, finalamp= amp1*(amp1/ampsum) + amp2*(amp2/ampsum) +..."
    sys.exit()

st_fft=EMData()
ma_fft=EMData()
wt_fft=EMData()

st_fft.set_complex(1)
ma_fft.set_complex(1)
wt_fft.set_complex(1)

st_fft.set_ri(1)
ma_fft.set_ri(1)
wt_fft.set_ri(1)

st_fft.to_zero()
ma_fft.to_zero()
wt_fft.to_zero()

outst="st-%s"%(args[1])
outma="ma-%s"%(args[1])
outwt="wt-%s"%(args[1])

thr=args[2]
#print args[3:]
#got1=0
for filename in args[3:] :
    next_image=EMData()
    next_image.read_image(filename)
    nx=next_image.get_xsize()
    ny=next_image.get_ysize()
    nz=next_image.get_zsize()
    next_fft=next_image.do_fft()
    key = filename
    value = next_fft
    imagedict.update({key:value})
st_fft.set_size(nx+2,ny,nz)
ma_fft.set_size(nx+2,ny,nz)
wt_fft.set_size(nx+2,ny,nz)
i=0
while i<nx:
    j=0
    print i
    while j<ny:
        k=0
#        print i,j
        while k<nz:
            a=avg_this_pixel (filename, imagedict, thr,i,j,k )
#            print a
            st_fft.set_value_at(i,j,k,a[0])
            st_fft.set_value_at(i+1,j,k,a[1])
            ma_fft.set_value_at(i,j,k,a[2])
            ma_fft.set_value_at(i+1,j,k,a[3])
            wt_fft.set_value_at(i,j,k,a[4])
            wt_fft.set_value_at(i+1,j,k,a[5])
#            print a
            k=k+1
        j=j+1
    i=i+2
print "here"
st_image=st_fft.do_ift()
ma_image=ma_fft.do_ift()
wt_image=wt_fft.do_ift()
print "and here"
st_image.write_image(outst)
ma_image.write_image(outma)
wt_image.write_image(outwt)
print "here too"
#print "coming out %f %f %f"%(a[0],a[1],a[2])
