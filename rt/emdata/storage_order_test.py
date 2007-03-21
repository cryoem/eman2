#!/usr/bin/env python

from EMAN2 import *

e = EMData()
e.set_size(2,3)
e.set_value_at(0,0,1)
e.set_value_at(0,1,2)
e.set_value_at(0,2,3)
e.set_value_at(1,0,4)
e.set_value_at(1,1,5)
e.set_value_at(1,2,6)

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
for j in range(3):
    for i in range(2):
        print e.get_value_at(i, j),
    print
