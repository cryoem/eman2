#!/usr/bin/env python
# Author: Stephen Murray (scmurray@bcm.edu), 12/12/11
# Copyright (c) 2000-2011 Baylor College of Medicine

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
import shutil
from subprocess import *

DEFANGLE = 	'0'
amplitude_contrast = 	'0.1'
RELION_RUN = 	"relion_refine_mpi "

progname = os.path.basename(sys.argv[0])
usage = """ prog [options] <set name> <reference map name>

This program will extract the necessary parameters from an EMAN2 refinement and 

Examples:


"""

print "Running e2refinetorelion3d.py"
# Program Options
parser = EMArgumentParser(usage, version=EMANVERSION)
parser.add_header(name="relion3dheader", help='Options below this label are specific to e2refinetorelion3d', title="### e2refinetorelion3d options ###", row=0, col=0, rowspan=1, colspan=3)
parser.add_pos_argument(name="set_name",help="The set name of the set of particles.", default="", guitype='filebox', browser="EMBrowserWidget(withmodal=True,multiselect=False)", row=1, col=0, rowspan=1, colspan=2)
parser.add_pos_argument(name="refmap", type=str, help="Reference Map", default="", guitype='filebox', browser="EMBrowserWidget(withmodal=True,multiselect=False)", row=2, col=0, rowspan=1, colspan=2)
parser.add_argument("--greyscale", action="store_true", help="Is the reference map in greyscale?", default=False, guitype='boolbox', row=3, col=0, rowspan=1, colspan=1)
parser.add_argument("--refctfcorrected", action="store_true", help="Has the reference map been ctf corrected?", guitype='boolbox', default=False,row=3, col=1, rowspan=1, colspan=1)
parser.add_argument("--solventmask",type=str, help="Location of the mask to be used", guitype='filebox',default="", browser="EMBrowserWidget(withmodal=True,multiselect=False)", filecheck=False,row=4, col=0, rowspan=1, colspan=2)
parser.add_argument("--symmgroup", help="Symmentry group", guitype='combobox', default='C',choicelist="""'C','D','T','O','I'""", row=5, col=0, rowspan=1, colspan=1)
parser.add_argument("--symmnumber",type=int,help="Symmetry number",default=1, guitype='intbox', row=5, col=1, rowspan=1, colspan=1)
parser.add_argument("--voltage",type=int, default=None, help="Voltage of the Microscope (kV)", guitype='intbox', row=6, col=0, rowspan=1, colspan=1)
parser.add_argument("--cs", type=float, default=None, help="Spherical Aberration", guitype='floatbox', row=6, col=1, rowspan=1, colspan=1)
parser.add_argument("--apix", type=float, default=None, help="Angstrom per pixel", guitype='floatbox', row=6, col=2, rowspan=1, colspan=1)
parser.add_argument("--imagemaskd", type=float, help="Diameter of the image mask", default=-1, guitype='floatbox', row=7, col=0, rowspan=1, colspan=1)
parser.add_argument("--defocus", type=float, help="Defocus in A", default=10000, guitype='floatbox', row=7, col=1, rowspan=1, colspan=1)
parser.add_argument("--amplitudecontrast", type=float, help="Amplitude Contrast value for the micrographs", default=0.1, guitype='floatbox', row=7, col=2, rowspan=1, colspan=2)
parser.add_argument("--numiter", type=int, help="# of iterations", default=25, guitype='intbox', row=8, col=0, rowspan=1, colspan=1)
parser.add_argument("--numclasses", type=int, help="# of classes", default=8, guitype='intbox', row=8, col=1, rowspan=1, colspan=1)
parser.add_argument("--ctfcorrect", action="store_true", help="(T/F)Do CTF Correction?", default=False, guitype='boolbox', row=9, col=0, rowspan=1, colspan=1)
parser.add_argument("--intensitycorrection", action="store_true", help="(T/F)Do intensity correction?", default=False, guitype='boolbox', row=9, col=1, rowspan=1, colspan=2)
parser.add_argument("--onlyflipphase", action="store_true", help="(T/F)Only flip phases?", default=False, guitype='boolbox', row=10, col=0, rowspan=1, colspan=1)
parser.add_argument("--dataphaseflipped", action="store_true", help="(T/F)Has the data been phase flipped already?", default=False, guitype='boolbox', row=10, col=1, rowspan=1, colspan=2)
parser.add_argument("--ignoretofirstpeak", action="store_true", help="(T/F)Ignore CTF's until the first peak?", default=False, guitype='boolbox', row=11, col=0, rowspan=1, colspan=2)
parser.add_argument("--maskrefstructure", action="store_true", help="(T/F)Mask reference structures?", default=False, guitype='boolbox', row=11, col=1, rowspan=1, colspan=2)
parser.add_argument("--splithalves", action="store_true", help="(T/F)Split data into random halves?", default=False, guitype='boolbox', row=12, col=0, rowspan=1, colspan=1)
parser.add_argument("--joinhalves", action="store_true", help="(T/F)Join random halves?", default=False, guitype='boolbox', row=12, col=1, rowspan=1, colspan=1)
parser.add_argument("--pad", type=int,help="Padding factor",default=2, guitype='combobox', choicelist='1,2,3', row=13 ,col=0, rowspan=1, colspan=1)
parser.add_argument("--lowpass",type=float,help="Initial low-pass filter (Ang)",default=60.0, guitype='floatbox', row=13, col=1, rowspan=1, colspan=1)
parser.add_argument("--regparam", type=float, help="Regularization Parameter T (weights experimental data vs. prior", default=1.0, guitype='floatbox', row=14, col=0, rowspan=1, colspan=1)
parser.add_argument("--healpix", type=str, default=30, help="Angular Sampling Interval (Degrees)", guitype='combobox', choicelist='30,15,7.5,3.7,1.8,0.9,0.5,0.2,0.1', row=14, col=1, rowspan=1, colspan=1)
parser.add_argument("--local", type=float,help="Perform local angular search at what range?",default=None, guitype='floatbox', row=15, col=0, rowspan=1, colspan=1)
parser.add_argument("--inplaneang", type=float, help="In-plane angular sampling", default=None, guitype='floatbox', row=15, col=1, rowspan=1, colspan=1)
parser.add_argument("--offsetrange", type=float, help="Offset search range (pix)", default=10.0, guitype='floatbox', row=16, col=0, rowspan=1, colspan=1) 
parser.add_argument("--offsetstep", type=float, help="Offset search step (pix)", default=2.0, guitype='floatbox', row=16, col=1, rowspan=1, colspan=1)
parser.add_argument("--oversampling",help="Oversampling order",default='1', guitype='combobox', choicelist='0,1,2', row=17, col=0, rowspan=1, colspan=1)
parser.add_argument("--threads", type=int, help="# of threads", default=1, guitype='intbox', row=18, col=0, rowspan=1, colspan=1)
parser.add_argument("--maxmemory", type=float, help="Maximum memory available for each node", guitype='floatbox', row=18, col=1, rowspan=1, colspan=1)
parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
optionList = pyemtbx.options.get_optionlist(sys.argv[1:])
(options, args) = parser.parse_args()

#Check for basic usage
if len(args) != 2:
   print "usage:" + usage
   print "Please run'" + progname + " -h' for detailed options"
   sys.exit(1)

current_dir = os.getcwd()
# Create the E2RLN directory structure if it does not exist
i = 1
found = 1
while found == 1:
	if i < 10:
		rln_run = '0' + str(i)
	else:
		rln_run = str(i)
	found = os.path.exists("relion3d_" + rln_run)
	i = i+1
E2RLN = "relion3d_" + rln_run
os.mkdir(E2RLN)

#E2n = E2init(sys.argv,options.ppid)
set_name, refmap = args[0], args[1]
s1 = "e2proc2d.py " + refmap + " " + E2RLN + "/3DRefMap.mrc --threed2threed --process=normalize.edgemean --verbose=0"
call(s1, shell=True)
header = EMData(set_name,0,True)
h_keys = header.get_attr_dict().keys()
nx,ny,nz=header['nx'],header['ny'],header['nz']
num_images = EMUtil.get_image_count(set_name)
mrc = False
if num_images > nz:
	num_ptcl = num_images
else:
	num_ptcl = nz
	mrc = True
project_db = db_open_dict("bdb:.#project")

if not header.get_attr_dict().__contains__('data_source'):
	print "The input stack/database of particles is invalid - it does not contain information tying each particle to a micrograph. Please provide another input file"
	print "Exiting e2refinetorelion3d"
	shutil.rmtree(E2RLN)
	exit(-1)

if options.apix == None:
	if header.get_attr_dict().__contains__('apix_x'):
		apix = header['apix_x']
	elif project_db.__contains__('global.apix'):
		apix = project_db['global.apix']
	else:
		print "An Angstrom per pixel was not found in the project database, the images themselves, and was not provided via a command line option. Please provide one"
		print "Exiting e2refinetorelion3d"
		shutil.rmtree(E2RLN)
		exit(-1)
else:
	apix = options.apix

if options.cs == None:
	if header.get_attr_dict().__contains__('cs'):
		cs = header['cs']
	elif project_db.__contains__('global.microscope_cs'):
		cs = project_db['global.microscope_cs']
	else:
		print "A spherical aberration value was not found in the project database, the images themselves, and was not provided via a command line option. Please provide one"
		print "Exiting e2refinetorelion3d"
		shutil.rmtree(E2RLN)
		exit(-1)
else:
	cs = options.cs
	
if options.voltage == None:
	if header.get_attr_dict().__contains__('voltage'):
		voltage = header['voltage']
	elif project_db.__contains__('global.microscope_voltage'):
		voltage = project_db['global.microscope_voltage']
	else:
		print "A microscope voltage was not found in the project database, the images themselves, and was not provided via a command line option. Please provide one (in kV)"
		print "Exiting e2refinetorelion3d"
		shutil.rmtree(E2RLN)
		exit(-1)
else:
	voltage = options.voltage	


## Create the particle stack files (.mrcs files) and the .star file that RELION needs as an inputs.
set_orig = set_name
#set_db = db_open_dict(set_orig)
i = 0
old_src = EMData(set_name,0)['data_source']

if mrc:
	s =  "e2proc2d.py " + set_orig + " " + E2RLN + "/ptcl_stack.hdf --threed2twod --verbose=0"
else:
	s =  "e2proc2d.py " + set_orig + " " + E2RLN + "/ptcl_stack.hdf --verbose=0"
call(s,shell=True)
ctf_corr = 0
for option1 in optionList:
	if option1 == "ctf":
		ctf_corr = 1
	if option1 == "amplitudecontrast":
		amplitude_contrast = str(options.amplitudecontrast)		
if ctf_corr == 1:
	s = "relion_star_loopheader rlnImageName rlnMicrographName rlnDefocusU rlnDefocusV rlnDefocusAngle rlnVoltage rlnSphericalAberration rlnAmplitudeContrast > " + E2RLN + "/all_images.star"
	DEF1 = DEF2 = str(options.defocus)
else:
	s = "relion_star_loopheader rlnImageName rlnMicrographName rlnVoltage rlnAmplitudeContrast > " + E2RLN + "/all_images.star"
call(s,shell=True)



print "Converting EMAN2 Files to Formats Compatible with RELION"
for k in range(num_ptcl):
	src = EMData(set_name,k)['data_source']
	if (src != old_src):
		s = "e2proc2d.py " + E2RLN + "/ptcl_stack.hdf" + " " + E2RLN + "/" + old_src.split('?')[0].replace('bdb:particles#','') + ".hdf --verbose=0 --first=" + str(i) + " --last=" + str(k-1)
		call(s, shell=True)
		if (k-i-1) == 0:
			s = "e2proc2d.py " + E2RLN + "/" + old_src.split('?')[0].replace('bdb:particles#','') + ".hdf " + E2RLN + "/" + old_src.split('?')[0].replace('bdb:particles#','') + ".mrc --verbose=0"
                        call(s, shell=True)
		else:
			s = "e2proc2d.py " + E2RLN + "/" + old_src.split('?')[0].replace('bdb:particles#','') + ".hdf " + E2RLN + "/" + old_src.split('?')[0].replace('bdb:particles#','') + ".mrc --verbose=0 --twod2threed" 
			call(s, shell=True)
		s1 = E2RLN + "/" + old_src.split('?')[0].replace('bdb:particles#','') + ".mrc"
		s2 = s1 + "s"
		shutil.move(s1, s2)
		s = "relion_star_datablock_stack " +  str(k-i) + " " +  E2RLN + "/" + old_src.split('?')[0].replace('bdb:particles#','') + ".mrcs " + E2RLN + "/" + old_src.split('?')[0].replace('bdb:particles#','') + ".mrcs " + str(voltage) + " " + str(amplitude_contrast) + "  >> " + E2RLN + "/all_images.star" 
		call(s,shell=True)
		s = "rm " + E2RLN + "/" + old_src.split('?')[0].replace('bdb:particles#','') + ".hdf" 
		call(s,shell=True)
		i = k
		old_src = src
	if (k+1) == num_ptcl:
		s = "e2proc2d.py " + E2RLN + "/ptcl_stack.hdf" + " " + E2RLN + "/" + src.split('?')[0].replace('bdb:particles#','') + ".hdf --verbose=0 --first=" + str(i) + " --last=" + str(k)
		call(s, shell=True)
		s = "e2proc2d.py " + E2RLN + "/" + src.split('?')[0].replace('bdb:particles#','') + ".hdf " + E2RLN + "/" + src.split('?')[0].replace('bdb:particles#','') + ".mrc --verbose=0"
                call(s, shell=True)

		s1 = E2RLN + "/" + src.split('?')[0].replace('bdb:particles#','') + ".mrc"
		s2 = s1 + "s"
		shutil.move(s1, s2)
		s = "relion_star_datablock_stack 1 "+  E2RLN + "/" + src.split('?')[0].replace('bdb:particles#','') + ".mrcs " + E2RLN + "/" + src.split('?')[0].replace('bdb:particles#','') + ".mrcs " + str(voltage) + " " + amplitude_contrast + "  >> " + E2RLN + "/all_images.star" 
		call(s,shell=True)
		s = "rm " + E2RLN + "/" + src.split('?')[0].replace('bdb:particles#','') + ".hdf" 
		call(s,shell=True)
		i = k
		old_src = src
		break


s = "rm " + E2RLN + "/ptcl_stack.hdf"
call(s,shell=True)
print "File Conversion Complete"

# Create the run directory structure if it does not exist
i = 1
found = 1
while found == 1:
	if i < 10:
		run_dir = '0' + str(i)
	else:
		run_dir = str(i)
	found = os.path.exists(E2RLN + "/run" + run_dir)
	i = i+1
RUNDIR = E2RLN + "/run" + run_dir
os.mkdir(RUNDIR)

#Parse the options and create the command to run Relion
s = RELION_RUN + "--i " + E2RLN + "/all_images.star --o " + RUNDIR + "/" + E2RLN + " --angpix " + str(apix)
grey = 0

for option1 in optionList:
	if option1 == "imagemaskd":
		s = s + " --particle_diameter " + str(options.imagemaskd)
	elif option1 == "numiter":
		s = s + " --iter " + str(options.numiter)
	elif option1 == "regparam":
		s = s + " --tau2_fudge " + str(options.regparam)
	elif option1 == "maskrefstructure":
		s = s + " --flatten_solvent"
	elif option1 == "intensitycorrection":
		s = s + " --scale"
	elif option1 == "splithalves":
		s = s + " --split_random_halves"
	elif option1 == "ctfcorrect":
		s = s + " --ctf"
	elif option1 == "onlyflipphase":
		s = s + " --only_flip_phases" 
	elif option1 == "dataphaseflipped":
		s = s + " --ctf_phase_flipped"
	elif option1 == "ignoretofirstpeak":
		s = s + " --ctf_intact_first_peak"
	elif option1 == "refctfcorrected":
		s = s + " --ctf_corrected_ref"
	elif option1 == "solventmask":
		s1 = "e2proc3d.py " + str(options.solventmask) + " " + E2RLN + "/mask.mrc"
		call(s1,shell=True)
		s = s + " --solvent_mask " + E2RLN + "/mask.mrc"
	elif option1 == "numclasses":
		s = s + " --K " + str(options.numclasses)
	elif option1 == "oversampling":
		s = s + " --oversampling " + str(options.oversampling)
	elif option1 == "inplaneang":
		s = s + " --psi_step " + str((2 ** int(options.oversampling)) * options.inplaneang)
	elif option1 == "offsetrange":
		s = s + " --offset_range " + str(options.offsetrange)
	elif option1 == "offsetstep":
		s = s + " --offset_step " + str(2 * options.offsetstep)
	elif option1 == "threads":
		s = s + " --j " + str(options.threads)
	elif option1 == "maxmemory":
		s = s + " --max_memory " + str(options.maxmemory)
	elif option1 == "greyscale":
		grey = 1
	elif option1 == "pad":
		s = s + " --pad " + str(options.pad)
	elif option1 == "lowpass":
		s = s + " --ini_high " + str(options.lowpass)
	elif option1 == "local":
		s = s + " --sigma_ang " + str(options.local / 3)
	elif option1 == "queue":
		submit_queue = str(options.queue)
	elif option1 == "healpix":
		if options.healpix == '30':
			s = s + " --healpix_order 0"
		elif options.healpix =='15':
			s = s + " --healpix_order 1"
		elif options.healpix =='7.5':
			s = s + " --healpix_order 2"
		elif options.healpix =='3.7':
			s = s + " --healpix_order 3"
		elif options.healpix =='1.8':
			s = s + " --healpix_order 4"
		elif options.healpix =='0.9':
			s = s + " --healpix_order 5"
		elif options.healpix =='0.5':
			s = s + " --healpix_order 6"
		elif options.healpix =='0.2':
			s = s + " --healpix_order 7"
		elif options.healpix =='0.1':
			s = s + " --healpix_order 8"
		else:
			print "Invalid angular sampling interval (--healpix). Defaulting to 30 degrees"
			s = s + " --healpix_order 0"
			
s = s + " --sym " + str(options.symmgroup) + str(options.symmnumber) + " --ref " + E2RLN + "/3DRefMap.mrc"

if grey != 1:
	s = s + " --firstiter_cc"
	

print "Generating qsub file"
###### create files to make the runs ######
qsub_file = current_dir + "/" + RUNDIR + "/qsub.sh"
f = open(qsub_file, 'w')
f.write("""#!/bin/sh

#PBS -V
#PBS -N """ + E2RLN + """
#PBS -l walltime=40:00:00
#PBS -l nodes=###NumberofNodes###:ppn=###ProcessorsPerNode###
#PBS -q ###QueueName###

cd ###DirectoryDataIsIn###
mpiexec -bynode -n ###NumberofNodes### """ + s + """ 
echo "done" """)

f.close()
print "************** Relion Command *******************"
print s
print "*************************************************"
print "e2refinetorelion3d.py has completed successfully"

