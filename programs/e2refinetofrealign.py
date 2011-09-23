#!/usr/bin/env python
# Author: Stephen Murray (scmurray@bcm.edu), 2/14/11
# Copyright (c) 2000-2011 Baylor Colelge of Medicine

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
CLS = 6                                # Number of images in the cls_result_## bdb
DALPHA = 4                             # index of the dalpha image
APERPIX = 100000                       # set as 10 microns
ANGAST = 0                             # eman2 throws out astygmatic particles so the angle is always 0
PRESA = 0                              # default
SPACE = ' '                            # to clean up the string manipulation later
FBEAUT = 'F'
FCREF = 'F'
IFSC = '0'
FSTAT = 'F'
IBLOW = '2'
RREC = '10.0'
RMAX1 = '200.0'
RMAX2 = '25.0'

progname = os.path.basename(sys.argv[0])
usage = """ prog [options] <name of refinement directory> <iteration number>

This program will extract the necessary parameters from an EMAN2 refinement and create a card file (card.txt)
and a particle meta data file in the correct formats for FreAlign. 

Examples:

The Eman2 refinement directory is refine_04 and the iteration you want to use within it is the 7th iteration:
e2refinetofrealign.py refine_04 7 

"""

parser = EMArgumentParser(usage,version=EMANVERSION)

parser.add_pos_argument(name="dir",help="The refinement directory to use for FreAlign.", default="", guitype='combobox', choicelist='glob.glob("refine*")', positional=True, row=0, col=0,rowspan=1, colspan=2)
parser.add_pos_argument(name="refineiter",help="The refinement iteration to use.", default="", guitype='intbox', positional=True, row=0, col=2,rowspan=1, colspan=1)
parser.add_header(name="frealignheader", help='Options below this label are specific to e2refinetofrealign', title="### e2refinetofrealign options ###", row=1, col=0, rowspan=1, colspan=3)
parser.add_argument("--fbeaut", action="store_true", help="(T/F)Apply extra real space symmetry averaging and masking to beautify final map prior to output", default=False, guitype='boolbox', row=2, col=0, rowspan=1, colspan=1)
parser.add_argument("--fcref", action="store_true", help="(T/F)Apply FOM filter to final reconstruction using function SQRT(2.0*FSC/(1.0+FSC))", default=False, guitype='boolbox', row=2, col=1, rowspan=1, colspan=1)
parser.add_argument("--fstat", action="store_true", help="(T/F)Calculate additional statistics in resolution table at end (QFACT, SSNR, CC, etc.). T Uses more than 50% more memory.", default=False, guitype='boolbox', row=2, col=2, rowspan=1, colspan=1)
parser.add_argument("--rrec", type=float, help="Resolution of reconstruction in angstroms. It is the resolution to which the reconstruction is calculated.", default = 10.0, guitype='floatbox', row=4, col=0, rowspan=1, colspan=2)
parser.add_argument("--reslow", type=float, help="Resolution of the data included in the alignment. This is the low resolution value. ex:200", default=200.0, guitype='floatbox', row=3, col=0, rowspan=1, colspan=2)
parser.add_argument("--reshigh", type=float, help="Resolution of the data included in the alignment. This is the high resolution value. ex:25", default=25.0, guitype='floatbox', row=3, col=2, rowspan=1, colspan=1)
parser.add_argument("--thresh", type=float, help="Phase Residual cutoff. Particles with a higher phase residual will not be included in the refinement ", default=90.0, guitype='floatbox', row=4, col=2, rowspan=1, colspan=1)

optionList = pyemtbx.options.get_optionlist(sys.argv[1:])

(options, args) = parser.parse_args()

if len(args) != 2:
   print "usage:" + usage
   print "Please run'" + progname + " -h' for detailed options"
   sys.exit(1)


E2n = E2init(sys.argv)
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
dir = "bdb:" + dir

db = db_open_dict(dir + "#projections_" + high)                     # Projections

num_images = EMUtil.get_image_count(dir + "#projections_" + high)   # Number of classes

cmd = db_open_dict(dir + "#register")                               # Register containing parameters used in the refinement
cmd = cmd['cmd_dict']
inpt1 = cmd['input']						    # Input filename for refinement
inpt = inpt1.replace('_phase_flipped-hp','')
inpt = inpt.replace('_phase_flipped','')
inpt = inpt.replace('_wiener_filtered','')
inpt = inpt+'_original_data'
part2 = db_open_dict(inpt1)
if not 'ctf' in part2[0].get_attr_dict():
   print "This command cannot be run with data that has not been CTF corrected"
   exit(-1)
part = db_open_dict(inpt)
tmp = part2[0]['ctf'].to_dict()                                      # CTF has things like Defocus, apix, etc

for i in range(num_images):
   classes[i]= db[i]['xform.projection'].get_rotation("eman")        # Has the angles needed to get the Euler Angles

nx = db[0]['nx']
cls_dir = dir + "#cls_result_" + high                               # CLS results. 

s = "e2proc2d.py " + inpt + " " + E2FA + "/particlestack.mrc --twod2threed --process=normalize.edgemean"
call(s, shell=True)
s = "e2proc2d.py " + dir + "#threed_filt_" + high + " " + E2FA + "/3DMapInOut.mrc --threed2threed --process=normalize.edgemean"
call(s, shell=True)

# Retrieve the pixel values of each image in the CLS stack.
for i in range(CLS):
   img = EMData(cls_dir, i)
   ny = img['ny']
   ptcl = {}
   for j in range(ny):
       ptcl[j] = img[SEP-1,j]        
   ptcl_data[i] = ptcl


# Write output file of particle meta data for FreAlign into the E2FA subdirectory created above
f = open(OUTFILE1,'w')
film = 0
bool_found = 0
film_dict = {}
for i in range(ny):
   ctf_dict = part2[i]['ctf'].to_dict()
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
   class_num = int(ptcl_data[0][i])
   if ptcl_data[5][i]:
      t = Transform({"type":"eman","az":classes[class_num]['az'],"alt":classes[class_num]['alt']*-1,"phi":ptcl_data[DALPHA][i]+180,"tx":ptcl_data[2][i],"ty":ptcl_data[3][i]})
   else:
      t = Transform({"type":"eman","az":classes[class_num]['az'],"alt":classes[class_num]['alt'],"phi":ptcl_data[DALPHA][i],"tx":ptcl_data[2][i],"ty":ptcl_data[3][i]})
   t = t.inverse()
   s = '{0:7d}{1:8.2f}{2:8.2f}{3:8.2f}{4:8.2f}{5:8.2f}{6:7.0f}.{7:6d}{8:9.1f}{9:9.1f}{10:8.2f}{11:7.2f}{12:6.2f}\n'.format(i+1, t.get_rotation("mrc")['phi'], t.get_rotation("mrc")['theta'], t.get_rotation("mrc")['omega'], t.get_trans()[0], t.get_trans()[1], mag, film + 1, defocus, defocus, ANGAST, PRESA, PRESA)
   f.write(s)
f.close()



####################################################################################################################
#														   #
# Card Creation For FreAlign. Please see the FreAlign README.TXT file for more information about each Card or Flag #
# 											 			   #
####################################################################################################################


RO = str(tmp['apix']*.375*(part[0]['nx']))

for option1 in optionList:
   if option1 == "fbeaut":
      FBEAUT = 'T'

   elif option1 == "fcref":
      FCREF = 'T'
  
   elif option1 == "rrec":
      RREC = str(options.rrec)

   elif option1 == "reslow":
      RMAX1 = str(options.reslow)

   elif option1 == "reshigh":
      RMAX2 = str(options.reshigh)

   elif option1 == "fstat":
         FSTAT = 'T'

OUTFILE2 = E2FA + "/card.txt"          # Cards required by FA
f = open(OUTFILE2, 'w')      # card.txt to be placed in the E2FA subdirectory created above

 
# Card 1
CFORM = 'M'
IFLAG = '1'
FMAG = FDEF = FASTIG = FPART = FMATCH = 'F'
IEWALD = '0'
s = CFORM + SPACE + IFLAG + SPACE + FMAG + SPACE + FDEF + SPACE + FASTIG + SPACE + FPART + SPACE + IEWALD + SPACE + FBEAUT + SPACE + FCREF + SPACE + FMATCH + SPACE + IFSC + SPACE + FSTAT + SPACE + IBLOW + '\n'
f.write(s)

# Card 2
RI = DANG = ITMAX = '0'
XSTD = '1'
PBC = '100'
BOFF = '60'
IPMAX = '10' 
PSIZE = str(tmp['apix']) 
WGH = str(tmp['ampcont']/100)
s = RO + SPACE + RI + SPACE + PSIZE + SPACE + WGH + SPACE + XSTD + SPACE + PBC + SPACE + BOFF + SPACE + DANG + SPACE + ITMAX + SPACE + IPMAX + '\n'
f.write(s)
 
# Card 3
MASK = '1 1 1 1 1'
f.write(MASK + '\n')

# Card 4
IFIRST = '1'
ILAST = str(ny)
f.write(IFIRST + SPACE + ILAST + '\n')

# Card 5
sym = cmd['sym']
if sym[0] == 'c':      # Convert from EMAN2 symmetry conventions to FreAlign symmetry conventions
   ASYM = 'C' + sym[1]
elif sym[0] == 'd':
   ASYM = 'D' + sym[1]
elif sym[0] == 'i':
   ASYM = 'I'
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
THRESH = 90
CS = tmp['cs']
AKV = tmp['voltage'] 
TX = TY = 0.0
s = str(RELMAG) + SPACE + str(DSTEP) + SPACE + str(TARGET) + SPACE + str(THRESH) + SPACE + str(CS) + SPACE + str(AKV) + SPACE + str(TX) + SPACE + str(TY) + "\n"
f.write(s)

# Card 7
DFSTD = '200.0'
RBFACT = '0.0'
s = RREC + SPACE + RMAX1 + SPACE + RMAX2 + SPACE + DFSTD + SPACE + RBFACT + '\n' 
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
