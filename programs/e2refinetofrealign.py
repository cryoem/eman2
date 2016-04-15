#!/usr/bin/env python
# Author: Stephen Murray (scmurray@bcm.edu), 2/14/11
# Copyright (c) 2000-2013 Baylor College of Medicine

# Official copyright notice. EMAN2 is distributed under a joint GPL/BSD license. Please copy
# this statement from one of the other programs. You must agree to use this license if your
# code is distributed with EMAN2. While you may use your own institution for the copyright notice
# the terms of the GPL/BSD license permit us to redistribute it.


#import block
from EMAN2 import *
from EMAN2db import db_open_dict
import pyemtbx.options
import os
import sys
from subprocess import *

SEP = 1                                # First is always best
CLS = 7                               # Number of images in the cls_result_## bdb
CLASS = 0
DX = 2
DY = 3
DALPHA = 4                             # index of the dalpha image
FLIP = 5
APERPIX = 100000                       # set as 10 microns
ANGAST = 0                             # eman2 throws out astygmatic particles so the angle is always 0
PRESA = 0                              # default
SPACE = ' '                            # to clean up the string manipulation later
FBEAUT = 'F'
FFILT = 'F'
FBFACT = 'F'
IFSC = '0'
FSTAT = 'F'
IMEM = '2'
RREC = '10.0'
RMAX1 = '200.0'
RMAX2 = '25.0'
IFLAG = '1'
progname = os.path.basename(sys.argv[0])
usage = """ prog [options] <name of refinement directory> <iteration number>

This program will extract the necessary parameters from an EMAN2 refinement and create a card file (card.txt)
and a particle meta data file in the correct formats for FreAlign. 

Examples:

The Eman2 refinement directory is refine_04 and the iteration you want to use within it is the 7th iteration:
e2refinetofrealign.py refine_04 7 

"""

parser = EMArgumentParser(usage,version=EMANVERSION)

parser.add_header(name="frealignheader", help='Options below this label are specific to e2refinetofrealign', title="### e2refinetofrealign options: Works with only FreAlign v9.07 and newer ###", row=0, col=0, rowspan=1, colspan=3)
parser.add_pos_argument(name="dir",help="The refinement directory to use for FreAlign.", default="", guitype='dirbox', dirbasename='refine',  row=1, col=0,rowspan=1, colspan=1)
parser.add_pos_argument(name="refineiter",help="The refinement iteration to use.", default="", guitype='intbox',  row=1, col=1,rowspan=1, colspan=1)
parser.add_argument("--mode", type=str, help="Mode to run FreAlign in: Mode 1 - Refinement and Reconstruction, Mode 3 - Simple Search and Refinement", default='1', guitype='combobox', choicelist="""'1','3'""", row=2, col=0, rowspan=1, colspan=1)
parser.add_argument("--fbeaut", action="store_true", help="(T/F)Apply extra real space symmetry averaging and masking to beautify final map prior to output", default=False, guitype='boolbox', row=4, col=0, rowspan=1, colspan=1)
parser.add_argument("--ffilt", action="store_true", help="(T/F)Apply Single Particle Wiener filter to final reconstruction", default=False, guitype='boolbox', row=4, col=1, rowspan=1, colspan=1)
parser.add_argument("--fstat", action="store_true", help="(T/F)Calculate additional statistics in resolution table at end (QFACT, SSNR, CC, etc.). T Uses more than 50 percent more memory.", default=False, guitype='boolbox', row=4, col=2, rowspan=1, colspan=1)
parser.add_argument("--reslow", type=float, help="Resolution of the data included in the alignment. This is the low resolution value. ex:200", default=200.0, guitype='floatbox', row=5, col=0, rowspan=1, colspan=2)
parser.add_argument("--reshigh", type=float, help="Resolution of the data included in the alignment. This is the high resolution value. ex:25", default=25.0, guitype='floatbox', row=5, col=2, rowspan=1, colspan=1)
parser.add_argument("--rrec", type=float, help="Resolution of reconstruction in angstroms. It is the resolution to which the reconstruction is calculated.", default = 10.0, guitype='floatbox', row=6, col=0, rowspan=1, colspan=1)
parser.add_argument("--rclas", type=float, help="High resloution limit used for classification", default = 10.0, guitype='floatbox', row=6, col=1, rowspan=1, colspan=1)
parser.add_argument("--thresh", type=float, help="Phase Residual cutoff. Particles with a higher phase residual will not be included in the refinement ", default=0.0, guitype='floatbox', row=6, col=2, rowspan=1, colspan=1)
parser.add_argument("--mass", default=0, type=float,help="The ~mass of the particle in kilodaltons", guitype='floatbox', row=7, col=0, rowspan=1, colspan=1)
parser.add_argument("--interp", type=str, help="Type of interpolation: 0 - Nearest Neighbor, 1 - Trilinear Interpolation (More Time-Consuming)", default='0', guitype='combobox', choicelist="""'0','1'""", row=7, col=1, rowspan=1, colspan=1)
parser.add_argument("--randomizemodel", help="Optionally randomize the phases of the initial model to this resolution (in Angstroms)", default=0,  guitype='floatbox', row=8, col=0, rowspan=1, colspan=2)
parser.add_argument("--imem", type=str, help="Memory Usage: 0 - Least Memory, 3 - Most memory", default='1', guitype='combobox', choicelist="""'0','1','2','3'""",row=8, col=2, rowspan=1, colspan=1)
parser.add_argument("--verbose",type=int, help="Level of verbose; how much information do you want the program to output?(0-9)", default=1, guitype='intbox',row=9, col=0, rowspan=1, colspan=1)
#parser.add_argument("--mass", default=0, type=float,help="The ~mass of the particle in kilodaltons, used to run normalize.bymass. Due to resolution effects, not always the true mass.", guitype='floatbox', row=12, col=0, rowspan=1, colspan=1, mode="refinement['self.pm().getMass()']")


parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
optionList = pyemtbx.options.get_optionlist(sys.argv[1:])

(options, args) = parser.parse_args()

if len(args) != 2:
   print "usage:" + usage
   print "Please run'" + progname + " -h' for detailed options"
   sys.exit(1)


E2n = E2init(sys.argv,options.ppid)
# Create the E2FA directory structure if it does not exist
dir = args[0]

i = 1
found = 1
while found == 1:
   if i < 10:
      fa_run = '0' + str(i)
   else:
      fa_run = str(i)
   found = os.path.exists("frealign_" + fa_run)
   i = i+1

E2FA = "frealign_" + fa_run
os.mkdir(E2FA)
meta_data_file = "ptcl_meta_data"
OUTFILE1 = E2FA + "/" + meta_data_file # meta data required by FA for each particle

ptcl_data = {}
classes = {}
mcrscp = {}

high = args[1]
if int(high) < 10:
   high = "0" + high

s = "e2proc2d.py " + dir + "/classes_" + high +"_even.hdf class_stack_temp.hdf --interlv=" + dir + "/classes_" + high +"_odd.hdf --verbose=" + str(options.verbose)
call(s,shell=True)
class_list = {}
class_list = EMData.read_images("class_stack_temp.hdf")
num_classes = len(class_list)   # Number of classes
az_list = {}
alt_list = {}
exc = 0
command_dict = js_open_dict(dir+"/0_refine_parms.json")
classes_even = EMData.read_images(dir + "/classes_" + high +"_even.hdf")
classes_odd = EMData.read_images(dir + "/classes_" + high +"_odd.hdf")


total_len = len(az_list)
even_set_name = class_list[0].get_attr_dict()['class_ptcl_src']
odd_set_name = class_list[1].get_attr_dict()['class_ptcl_src']
all_set_name = "sets/"+base_name(even_set_name)+"_ptcls.lst"
all_set_data = EMData.read_images(all_set_name)
print all_set_name

s = "e2proc2d.py " + all_set_name + " " + E2FA + "/particlestack.mrc --twod2threed --process=normalize.edgemean --mult=-1 --verbose=" + str(options.verbose)
call(s, shell=True)
s = "e2proc3d.py " + dir + "/threed_" + high + ".hdf " + E2FA + "/3DMapInOut.mrc --process=normalize.edgemean --verbose=" + str(options.verbose)
call(s, shell=True)
s = "e2proc2d.py " + dir + "/cls_result_" + high + "_even.hdf cls_stack_temp.hdf --interlv=" +  dir + "/cls_result_" + high + "_odd.hdf --writejunk --verbose=" + str(options.verbose)
call(s, shell=True)
#brute force way to do it but need to make sure its right...
cls_even = EMData.read_images(dir + "/cls_result_" + high + "_even.hdf")
cls_odd = EMData.read_images(dir + "/cls_result_" + high + "_odd.hdf")
cls_class_list_even = cls_even[CLASS]
cls_class_list_odd = cls_odd[CLASS]
dx_list_even = cls_even[DX]
dx_list_odd = cls_odd[DX]
dy_list_even = cls_even[DY]
dy_list_odd = cls_odd[DY]
dalpha_list_even = cls_even[DALPHA]
dalpha_list_odd = cls_odd[DALPHA]
flip_list_even = cls_even[FLIP]
flip_list_odd = cls_odd[FLIP]

dx_list = {}
dy_list = {}
dalpha_list = {}
flip_list = {}
cls_class_list = {}

for i in range(cls_class_list_even.get_attr_dict()['ny']):
	cls_class_list[2*i] = cls_class_list_even[0,i]
	dx_list[2*i] = dx_list_even[0,i]
	dy_list[2*i] = dy_list_even[0,i]
	dalpha_list[2*i] = dalpha_list_even[0,i]
	flip_list[2*i] = flip_list_even[0,i]

for i in range(cls_class_list_odd.get_attr_dict()['ny']):
	cls_class_list[2*i+1] = cls_class_list_odd[0,i]
	dx_list[2*i+1] = dx_list_odd[0,i]
	dy_list[2*i+1] = dy_list_odd[0,i]
	dalpha_list[2*i+1] = dalpha_list_odd[0,i]
	flip_list[2*i+1] = flip_list_odd[0,i]

ny = len(dx_list)
   
## Write output file of particle meta data for FreAlign into the E2FA subdirectory created above
f = open(OUTFILE1,'w')
film = 0
bool_found = 0
film_dict = {}
ctf_dict = {}

for i in range(len(all_set_data)):
   ptcl = all_set_data[i]
   ctf_dict = ptcl['ctf'].to_dict()
   defocus = ctf_dict['defocus'] * 10000
   for item in film_dict:
      if ctf_dict == film_dict[item]:
         film = item
         bool_found = 1
   if not bool_found:
      film = len(film_dict)
      film_dict[film] = ctf_dict
   bool_found = 0
   mag = APERPIX / ctf_dict['apix']
   apix_shift = ctf_dict['apix']
   class_num = cls_class_list[i]
   if i%2==0:
      az=classes_even[int(class_num)].get_attr_dict()['xform.projection'].get_rotation("eman")['az']
      alt=classes_even[int(class_num)].get_attr_dict()['xform.projection'].get_rotation("eman")['alt']
   else:
      az=classes_odd[int(class_num)].get_attr_dict()['xform.projection'].get_rotation("eman")['az']
      alt=classes_odd[int(class_num)].get_attr_dict()['xform.projection'].get_rotation("eman")['alt']
   if flip_list[i] == 1:
      t = Transform({"type":"eman","az":az+180,"alt":180-alt,"phi":dalpha_list[i]*-1,"tx":dx_list[i],"ty":dy_list[i]})
   else:
      t = Transform({"type":"eman","az":az,"alt":alt,"phi":dalpha_list[i],"tx":dx_list[i],"ty":dy_list[i]})
   t = t.inverse()
   s = '{0:7d}{1:8.2f}{2:8.2f}{3:8.2f}{4:8.2f}{5:8.2f}{6:7.0f}.{7:6d}{8:9.1f}{9:9.1f}{10:8.2f}{11:7.2f}{12:6.2f}\n'.format(i+1, t.get_rotation("spider")['phi'], t.get_rotation("spider")['theta'], t.get_rotation("spider")['psi'], t.get_trans()[0]*apix_shift, t.get_trans()[1]*apix_shift, mag, film + 1, defocus, defocus, ANGAST, PRESA, PRESA)
   #s = '{0:7d}{1:8.2f}{2:8.2f}{3:8.2f}{4:8.2f}{5:8.2f}{6:7.0f}.{7:6d}{8:9.1f}{9:9.1f}{10:8.2f}{11:7.2f}{12:6.2f}\n'.format(i+1, t.get_rotation("eman")['az'], t.get_rotation("eman")['alt'], t.get_rotation("eman")['phi'], t.get_trans()[0], t.get_trans()[1], mag, film + 1, defocus, defocus, ANGAST, PRESA, PRESA)
   f.write(s)
s = "rm class_stack_temp.hdf cls_stack_temp.hdf "
call(s,shell=True)
f.close()

#####################################################################################################################
##														   #
## Card Creation For FreAlign. Please see the FreAlign README.TXT file for more information about each Card or Flag #
## 											 			   #
#####################################################################################################################


RO = str(command_dict['apix']*.375*(all_set_data[0]['nx']))
for option1 in optionList:
	if option1 == "fbeaut":
		FBEAUT = 'T'
	elif option1 == "ffilt":
		FFILT = 'T'
  	elif option1 == "rrec":
		RREC = str(options.rrec)
	elif option1 == "reslow":
		RMAX1 = str(options.reslow)
	elif option1 == "reshigh":
		RMAX2 = str(options.reshigh)
	elif option1 == "rclas":
		RCLAS = str(options.rclas)
	elif option1 == "fstat":
		FSTAT = 'T'
	elif option1 == "mass":
		MASS = str(options.mass)
	elif option1 == "interp":
		INTERP = str(options.interp)
	elif option1 == "mode":
		IFLAG = str(options.mode)
	elif option1 == "randomizemodel":
		if float(options.randomizemodel) != 0.0:
			s1 = "e2proc3d.py " + E2FA + "/3DMapInOut.mrc " + E2FA + "/3DMapInOut.mrc --process=filter.lowpass.randomphase:apix=" + str(command_dict['apix']) + ":cutoff_freq=" + str(1/float(options.randomizemodel))
			call(s1, shell=True)

OUTFILE2 = E2FA + "/card.txt"          # Cards required by FA
f = open(OUTFILE2, 'w')      # card.txt to be placed in the E2FA subdirectory created above

 
# Card 1
CFORM = 'M'

FMAG = FDEF = FASTIG = FPART = FMATCH = 'F'
IEWALD = '0'
s = CFORM + SPACE + IFLAG + SPACE + FMAG + SPACE + FDEF + SPACE + FASTIG + SPACE + FPART + SPACE + IEWALD + SPACE + FBEAUT + SPACE + FFILT + SPACE + FBFACT + SPACE + FMATCH + SPACE + IFSC + SPACE + FSTAT + SPACE + IMEM + SPACE + INTERP + '\n'
f.write(s)

# Card 2
RI = DANG = ITMAX = '0'
XSTD = '0'
PBC = '20'
BOFF = '30'
IPMAX = '10' 
PSIZE = str(command_dict['apix']) 
WGH = str(ctf_dict['ampcont']/100)
MASS = str(options.mass)
s = RO + SPACE + RI + SPACE + PSIZE + SPACE + MASS + SPACE + WGH + SPACE + XSTD + SPACE + PBC + SPACE + BOFF + SPACE + DANG + SPACE + ITMAX + SPACE + IPMAX + '\n'
f.write(s)
 
# Card 3
MASK = '1 1 1 1 1'
f.write(MASK + '\n')

# Card 4
IFIRST = '1'
ILAST = str(ny)
f.write(IFIRST + SPACE + ILAST + '\n')

# Card 5
sym = command_dict['sym']
if sym[0] == 'c':      # Convert from EMAN2 symmetry conventions to FreAlign symmetry conventions
   ASYM = 'C' + sym[1]
elif sym[0] == 'd':
   ASYM = 'D' + sym[1]
elif sym[0] == 'i':
   ASYM = 'I'
   s1 = "e2proc3d.py " + E2FA + "/3DMapInOut.mrc " + E2FA + "/3DMapInOut.mrc --icos5to2"
   call(s1,shell=True)
elif sym[0] == 'o':
   ASYM = 'O'
elif sym[0] == 't':
   ASYM = 'T'
else:
   ASYM = sym
f.write(ASYM + '\n')

# Card 6
RELMAG = 1.0
DSTEP = 10 
TARGET = 15
THRESH = 0.0
CS = ctf_dict['cs']
AKV = ctf_dict['voltage'] 
TX = TY = 0.0
s = str(RELMAG) + SPACE + str(DSTEP) + SPACE + str(TARGET) + SPACE + str(THRESH) + SPACE + str(CS) + SPACE + str(AKV) + SPACE + str(TX) + SPACE + str(TY) + "\n"
f.write(s)

# Card 7
DFSTD = '200.0'
RBFACT = '0.0'
s = RREC + SPACE + RMAX1 + SPACE + RMAX2 + SPACE + RCLAS + SPACE + DFSTD + SPACE + RBFACT + '\n' 
f.write(s)

# Card 8
FINPAT1 = 'particlestack.mrc'
f.write(FINPAT1 + '\n')

# Card 9
FINPAT2 = '/dev/null'
f.write(FINPAT2 + '\n')

# Card 10
FINPAR = meta_data_file
f.write(FINPAR + '\n')

# Card 11
FOUTPAR = 'OutParam'
f.write(FOUTPAR + '\n')

# Card 12
FOUTSH = 'OutParamShift'
f.write(FOUTSH + '\n')

# Card 6 again because it is required to deliniate that there is only one data set
f.write("0 0 0 0 0 0 0 0\n")

# Card 13
F3D = '3DMapInOut.mrc'
f.write(F3D + '\n')

# Card 14
FWEIGH = '3DWeights.mrc'
f.write(FWEIGH + '\n')


MAP1 = '3DMap_Odd.mrc'
MAP2 = '3DMap_Even.mrc'

# Card 15
f.write(MAP1 + '\n')

# Card 16
f.write(MAP2 + '\n')

# Card 17
FPHA = '3DPhaseResidual.mrc'
f.write(FPHA + '\n')

# Card 18
FPOI = '3DPointSpreadFunc.mrc'
f.write(FPOI + '\n')

f.close()

print "e2refinetofrealign.py finished"

E2end(E2n)
