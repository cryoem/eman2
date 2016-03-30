from EMAN2 import *
from time import time

a=test_image(size=(64,64))

t1=time()
print "HDF"
for j in range(50):
	for i in range(10000): a.write_image("xx.hdf",-1)
	t2=time()
	os.system("sync")
	t3=time()
	print j*10000,t3-t1

t1=time()
print "BDB"
for j in range(50):
	for i in range(10000): a.write_image("bdb:test",-1)
	t2=time()
	os.system("sync")
	t3=time()
	print j*10000,t3-t1

print ""
t1=time()
for i in range(20000): a.write_image("x.hed",-1)
t2=time()
os.system("sync")
t3=time()
print "hed",t2-t1,t3-t1

t1=time()
for i in range(20000): a.write_image("x.hdf",-1)
t2=time()
os.system("sync")
t3=time()
print "hdf",t2-t1,t3-t1

t1=time()
for i in range(20000): a.write_image("x.spi",-1)
t2=time()
os.system("sync")
t3=time()
print "spi",t2-t1,t3-t1

t1=time()
for i in range(20000): a.write_image("bdb:x",-1)
t2=time()
os.system("sync")
t3=time()
print "bdb",t2-t1,t3-t1

