#!/usr/bin/env python

import os
import sys

if len(sys.argv) != 2:
   print "Please use e2fa.py <number of iterations to run FreAlign>"
   exit(-1)
iter = int(sys.argv[1])
os.system('clear')

os.system('cp 3DMapInOut.mrc 3DMapInOut.mrc.old')
os.system('cp ptcl_meta_data ptcl_meta_data.old')
for i in range(iter):
   os.system('frealign_v8.exe < card.txt')
   if i < 10:
      k = '0' + str(i) 
   else:
      k = str(i)
   s = "cp 3DMapInOut.mrc OutputMap_" + k +".mrc"
   os.system(s)
   s = "cp OutParam OutParam_" + k
   os.system(s)
   os.system('mv OutParam ptcl_meta_data')
   s = "cp OutParamShift OutParamShift_" + k
   os.system(s)

os.system('mv 3DMapInOut.mrc.old 3DMapInOut.mrc')
os.system('mv ptcl_meta_data.old ptcl_meta_data')

print "e2runfrealign.py finished"

exit(0)
