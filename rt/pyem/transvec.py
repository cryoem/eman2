#!/bin/env python

from EMAN2 import *

m1 = Matrix4f([2,4,9,10, 1,9,3,12, 6,4,3,9, 12,3,9,1])
t = Transform(m1)

v1 = [1,0,0]
v2 = [1,1,1]
v3 = [1,2,3]

vs = [v1,v2,v3]

for v in vs:
	newv = v * t
	print newv.get_as_list()
	
