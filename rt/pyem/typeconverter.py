#!/bin/env python

from EMAN2 import *

def cmp(l1, t1):
	assert(l1[0] == t1[0] and l1[1] == t1[1] and l1[2] == t1[2])

nlist = [10, 20, 30]

vec3i = TestUtil.from_Vec3i()
v1 = vec3i.as_list()
TestUtil.to_Vec3i(nlist)
cmp(v1, nlist)

ip1 = TestUtil.from_IntPoint()
TestUtil.to_IntPoint(nlist)
cmp(ip1, nlist)

