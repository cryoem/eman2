#!/usr/bin/env python
# Author: Stephen Murray (scmurray@bcm.edu), 2/14/11
# Copyright (c) 2000-2011 Baylor Colelge of Medicine

# Official copyright notice. EMAN2 is distributed under a joint GPL/BSD license. Please copy
# this statement from one of the other programs. You must agree to use this license if your
# code is distributed with EMAN2. While you may use your own institution for the copyright notice
# the terms of the GPL/BSD license permit us to redistribute it.

# import block
import os
import sys

if len(sys.argv) != 1:
   print "Please use e2runfrealign.py"
   exit(-1)

E2n=E2init(sys.argv)

os.system('clear')
os.system('cp 3DMapInOut.mrc 3DMapInOut.mrc.old')
os.system('cp ptcl_meta_data ptcl_meta_data.old')
dir_list = os.listdir('.')
dir_list.sort()
high = 0
for item in dir_list:
   if len(item) > 6:
      if item[:4] == 'card':
         item = item.replace(".txt",'')
         item = item.replace("card",'')
         if int(item) > high:
            high = int(item)

if high == 0:
   high = high + 1
for i in range(high):
   os.system('frealign_v8.exe < card' + str(i) + '.txt')
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
   os.system('cp 3DMapInOut.mrc.old 3DMapInOut.mrc')

os.system('mv 3DMapInOut.mrc.old 3DMapInOut.mrc')
os.system('mv ptcl_meta_data.old ptcl_meta_data')

print "e2runfrealign.py finished"

E2end(E2n)
