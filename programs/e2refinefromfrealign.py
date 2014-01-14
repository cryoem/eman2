#!/usr/bin/env python
# Author: Stephen Murray (scmurray@bcm.edu), 2/14/11
# Copyright (c) 2000-2011 Baylor Colelge of Medicine

# Official copyright notice. EMAN2 is distributed under a joint GPL/BSD license. Please copy
# this statement from one of the other programs. You must agree to use this license if your
# code is distributed with EMAN2. While you may use your own institution for the copyright notice
# the terms of the GPL/BSD license permit us to redistribute it.

# import block
from EMAN2 import *
from EMAN2db import db_open_dict
import pyemtbx.options
import os
import sys
from subprocess import *


IN_META_FILE = "ptcl_meta_data"
OUTFILE = "diff.txt"
progname = os.path.basename(sys.argv[0])
usage = """ [options] <name of frealign directory>
Help info needs to go here.....
"""

parser = EMArgumentParser(usage,version=EMANVERSION)

parser.add_pos_argument(name="frealigndir",help="The Frealign directory to use.", default="", guitype='dirbox', dirbasename='frealign',  row=0, col=0,rowspan=1, colspan=3)
parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
parser.add_argument("--icosahedral_symmetry", action="store_true", help="Does it have icosahedral symmetry?", default=False, guitype='boolbox', row=1, col=0, rowspan=1, colspan=1)
#row=2, col=0, rowspan=1, colspan=1)

optionList = pyemtbx.options.get_optionlist(sys.argv[1:])

(options, args) = parser.parse_args()

if len(args) != 1:
   print "usage:" + usage
   print "Please run'" + progname + " -h' for detailed options"
   sys.exit(1)

E2n = E2init(args,options.ppid)

dir = args[0]

in_img = EMData(dir + "/3DMapInOut.mrc")
fa_out_img = EMData(dir + "/OutputMap.mrc")
a = fa_out_img.process("normalize.toimage", {"to": in_img, "ignore_zero":1})
a.write_image(dir + "/OutputMap_Normalized.mrc")

for option1 in optionList:
	if option1 == "icosahedral_symmetry":
		s1 = "e2proc3d.py " + dir + "/OutputMap_Normalized.mrc " + dir + "/OutputMap_Normalized --icos2to5"
		call(s,shell=True)
		

OUT_META_FILE = "OutParam"
FINAL_META_FILE = OUT_META_FILE
#fsc_dict = db_open_dict("bdb:" + dir + "#convergence.results")
#a = EMData(dir + "/3DMapInOut.mrc")
#b = EMData(dir + "/OutputMap_Normalized.mrc")
#s =  "init_00_fsc"
#fsc = a.calc_fourier_shell_correlation(b)
#fsc_len = len(fsc)
#for i in range(0,fsc_len/3-1):
   #fsc[i] = fsc[i]/b['apix_x']
#fsc_dict[s] = [fsc[:fsc_len/3-1],fsc[fsc_len/3:fsc_len/3*2-1]]  

#f=open(dir + "/" + FINAL_META_FILE, 'r')
#l3=f.readlines()
#f.close()
#nyquist = []
#fsc_val = []
#start = 0

#for i in range(len(l3)):
   #lst=l3[i].split()
   #if l3[i][0] == 'C' and len(lst) > 1:
      #if lst[1] =='Average' and lst[2] != 'phase':
         #break
   #if start == 1:
      #nyquist.append(1/float(lst[2]))
      #fsc_val.append(float(lst[5]))
   #if l3[i][0] == 'C' and len(lst) > 1:
      #if lst[1] =='NO.':
         #start = 1
#s = 'even_odd_fsc'
#fsc_dict[s] = [nyquist,fsc_val]

   
f = open(dir + "/" + IN_META_FILE, 'r')
l2 = f.readlines()
f.close()
ptcl_dict_in = {}
k = 0
for i in range(len(l2)):
   ptcl = {}
   lst = l2[i].split()
   if l2[i][0] == 'C':
      k = k+1
   else:
      ptcl[0] = lst[1] # Psi
      ptcl[1] = lst[2] # Theta
      ptcl[2] = lst[3] # Phi
      ptcl[3] = lst[4] # Shift in X
      ptcl[4] = lst[5] # Shift in Y
      ptcl_dict_in[i-k] = ptcl
f = open(dir + "/" + OUT_META_FILE, 'r')
l = f.readlines()
f.close()
ptcl_dict_out = {}
k = 0
for i in range(len(l)):
   ptcl = {}
   lst = l[i].split()
   if l[i][0] == 'C':
      k = k+1
   else:
      ptcl[0] = lst[1] # Psi
      ptcl[1] = lst[2] # Theta
      ptcl[2] = lst[3] # Phi
      ptcl[3] = lst[4] # Shift in X
      ptcl[4] = lst[5] # Shift in Y
      ptcl_dict_out[i-k] = ptcl


if len(ptcl_dict_out) != len(ptcl_dict_in):
   print "3D map has been fixed but FreAlign was not run in mode 1 so no fsc can be done"
   exit(-1)
f = open(dir + "/" + OUTFILE,'w')
f.write("Particle    Psi-Diff   Theta-Diff   Phi-Diff   X-Shift-Diff  Y-Shift-Diff\n")
phi_sum = theta_sum = psi_sum = x_sum = y_sum = 0
for i in range(len(ptcl_dict_in)):
   s='{0:8d}{1:12.3f}{2:13.3f}{3:11.3f}{4:15.3f}{5:14.3f}\n'.format(i+1, float(ptcl_dict_in[i][0])-float(ptcl_dict_out[i][0]), float(ptcl_dict_in[i][1])-float(ptcl_dict_out[i][1]), float(ptcl_dict_in[i][2])-float(ptcl_dict_out[i][2]), float(ptcl_dict_in[i][3])-float(ptcl_dict_out[i][3]), float(ptcl_dict_in[i][4])-float(ptcl_dict_out[i][4]))
   if abs(float(ptcl_dict_in[i][0])-float(ptcl_dict_out[i][0])) > 180:
      phi_sum = phi_sum + 360 - abs(float(ptcl_dict_in[i][0])-float(ptcl_dict_out[i][0]))
   else:       
      phi_sum = phi_sum + abs(float(ptcl_dict_in[i][0])-float(ptcl_dict_out[i][0]))

   if abs(float(ptcl_dict_in[i][1])-float(ptcl_dict_out[i][1])) > 180:
      theta_sum = theta_sum + 360 - abs(float(ptcl_dict_in[i][1])-float(ptcl_dict_out[i][1]))
   else:
      theta_sum = theta_sum + abs(float(ptcl_dict_in[i][1])-float(ptcl_dict_out[i][1]))

   if abs(float(ptcl_dict_in[i][2])-float(ptcl_dict_out[i][2])) > 180:
      psi_sum = psi_sum + 360 - abs(float(ptcl_dict_in[i][2])-float(ptcl_dict_out[i][2]))
   else:
      psi_sum = psi_sum + abs(float(ptcl_dict_in[i][2])-float(ptcl_dict_out[i][2]))
 
      x_sum = x_sum + abs(float(ptcl_dict_in[i][3])-float(ptcl_dict_out[i][3]))

      y_sum = y_sum + abs(float(ptcl_dict_in[i][4])-float(ptcl_dict_out[i][4]))            
   f.write(s)    
s=' Average{0:12.3f}{1:13.3f}{2:11.3f}{3:15.3f}{4:14.3f}\n'.format(phi_sum/len(ptcl_dict_in), theta_sum/len(ptcl_dict_in), psi_sum/len(ptcl_dict_in), x_sum/len(ptcl_dict_in), y_sum/len(ptcl_dict_in))
f.write(s)
f.close()

print "e2refinefromfrealign.py finished"

E2end(E2n)
