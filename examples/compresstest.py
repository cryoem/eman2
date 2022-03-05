#!/usr/bin/env python

from EMAN2 import *

# unsigned int
a=EMData(4,4,1)
for i in range(16): a[i]=i

a.write_compressed("ctest.hdf",0,5)
a.write_compressed("ctest.hdf",1,4)
a.write_compressed("ctest.hdf",2,3)
a.write_compressed("ctest.hdf",3,2)

b=EMData("ctest.hdf",0)
c=EMData("ctest.hdf",1)
d=EMData("ctest.hdf",2)
e=EMData("ctest.hdf",3)

for i in range(16):
	print(f"{a[i]:1.0f}\t{b[i]:2.3f}\t{c[i]:2.3f}\t{d[i]:2.3f}\t{e[i]:2.3f}")

os.unlink("ctest.hdf")
