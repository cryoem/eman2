#!/usr/bin/env python

from EMAN2 import *

try:
    stack = sys.argv[1]
#    if stack[-3:] != "hdf":
#        print("Stack must be in HDF format")
#        exit(1)

    a = EMData(stack,0)

#    if a['nz'] <= 1:
#        print("You must specify a 3D stack")
#        exit(1)

    b = EMData(stack,1)
    a.write_image(stack.replace('.hdf','_0.hdf'))
    b.write_image(stack.replace('.hdf','_1.hdf'))

except:
    print("You must specify a 3D HDF stack.")
    print("Example: e2split3dstack.py 3dstack.hdf")
    print("Result will be stored in the directory containing the 3D stack.")

