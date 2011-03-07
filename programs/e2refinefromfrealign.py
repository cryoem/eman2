#!/usr/bin/env python

########################################################################################################
# This program is designed to parse the output files from FreAlign and to convert the data back to     #
# Eman2.                                                                                               #
# Developed in the Ludtke Lab by SCM                                                                   #
# Version 1.0                                                                                          #
######################################################################################################## 

from EMAN2 import *
from optparse import OptionParser
import pyemtbx.options
import os
import sys
IN_META_FILE = "ptcl_meta_data"
OUTFILE = "diff.txt"

progname = os.path.basename(sys.argv[0])
usage = progname + """ [options] <name of frealign directory>

"""


if len(sys.argv) != 2:
   print "usage:" + usage
   print "Please run'" + progname + " -h' for detailed options"
   sys.exit(1)

high = 0
dir = sys.argv[1]

for item in os.listdir(dir):
   if len(item) == 11:
      if item[0:8] == 'OutParam':
         if int(item[9:11]) > high:
            high = int(item[9:11])


for l_iter in range(high+1):
   if l_iter < 10:
      iter = '0' + str(l_iter)
   else:
      iter = str(l_iter)

   in_img = EMData(dir + "/3DMapInOut.mrc")
   fa_out_img = EMData(dir + "/OutputMap_" + iter + ".mrc")
   a = fa_out_img.process("normalize.toimage", {"to": in_img, "ignore_zero":1})
   a.write_image(dir + "/OutputMap_Normalized_" + iter + ".mrc")


OUT_META_FILE = "OutParam"
FINAL_META_FILE = OUT_META_FILE + "_" + iter
fsc_dict = db_open_dict("bdb:" + dir + "#convergence.results")
a = EMData(dir + "/3DMapInOut.mrc")
b = EMData(dir + "/OutputMap_Normalized_00.mrc")
s =  "00_init_fsc"
fsc = a.calc_fourier_shell_correlation(b)
fsc_len = len(fsc)
fsc_dict[s] = [fsc[:fsc_len/3-1],fsc[fsc_len/3:fsc_len/3*2-1]]




for k in range(1,high):
   if k < 10:
      iter = '0' + str(k)    
   else:
      iter = str(k)   
   if k-1 < 10:
      iter_minus = '0' + str(k-1)
   else:
      iter_minus = str(k-1)


   a = EMData(dir + "/OutputMap_Normalized_" + iter + ".mrc")
   b = EMData(dir + "/OutputMap_Normalized_" + iter_minus + ".mrc")
   s = iter_minus + "_" + iter + "_fsc"
   fsc = a.calc_fourier_shell_correlation(b)
   fsc_len = len(fsc)
   fsc_dict[s] = [fsc[:fsc_len/3-1],fsc[fsc_len/3:fsc_len/3*2-1]]

f=open(dir + "/" + FINAL_META_FILE, 'r')
l3=f.readlines()
f.close()
ring_rad = []
fsc_val = []
start = 0

for i in range(len(l3)):
   lst=l3[i].split()
   if l3[i][0] == 'C' and len(lst) > 1:
      if lst[1] =='Average':
         break
   if start == 1:
      ring_rad.append(int(lst[3]))
      fsc_val.append(int(lst[5]))
   if l3[i][0] == 'C' and len(lst) > 1:
      if lst[1] =='NO.':
         start = 1
s = 'even_odd_' + iter + '_fsc'
fsc_dict[s] = [ring_rad,fsc_val]

   
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
f = open(dir + "/" + OUT_META_FILE + "_" + iter, 'r')
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
