#!/bin/env python
from EMAN2 import *
from sparx  import *

mode = "F"

a = test_image()

b = rot_shift2D(a,12.6,3.1,-2.2)

numr = Numrinit(1,51)
#print "START ",ttime()
cimage=Util.Polar2Dm(a, 65, 65, numr, mode)
#print "Polar ",ttime()

Util.Frngs(cimage, numr)
wr=ringwe(numr ,mode)
Applyws(cimage, numr, wr)
#print "Frngs",ttime()

print  ormq(b, cimage,4,4,1,mode,numr,65,65)
print  ornq(b, cimage,4,4,1,mode,numr,65,65)

a=mirror(a)
cimage=Util.Polar2Dm(a, 65, 65, numr, mode)
Util.Frngs(cimage, numr)
Applyws(cimage, numr, wr)

print  ormq(b, cimage,4,4,1,mode,numr,65,65)
print  ornq(b, cimage,4,4,1,mode,numr,65,65)
