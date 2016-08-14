#!/usr/bin/env python

from EMAN2 import *

e = EMData()

nx = 3
ny = 2
e.set_size(nx,ny)
e.set_value_at(0,0,0.0)
e.set_value_at(1,0,1.0)
e.set_value_at(2,0,2.0)
e.set_value_at(0,1,3.0)
e.set_value_at(1,1,4.0)
e.set_value_at(2,1,5.0)

d = e.get_2dview()
print 'data in boost array, get_2dview():'
print d

print
print 'data in numpy array:'
a = EMNumPy.em2numpy(e)
print "shape is: ", a.shape
print a

print
print 'data by get_value_at(x,y):' 
for j in range(2):
    for i in range(3):
        print e.get_value_at(i, j),
    print
