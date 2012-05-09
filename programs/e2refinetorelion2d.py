#!/usr/bin/env python
# Author: Stephen Murray (scmurray@bcm.edu), 12/05/11
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
import shutil
from subprocess import *
from collections import defaultdict
import pprint

DEFANGLE = 	'0'
amplitude_contrast = '0.1'
RELION_RUN = 	"relion_refine "

progname = os.path.basename(sys.argv[0])
usage = """ prog [options] <set name>

This program will extract the necessary parameters from an EMAN2 refinement and determine class averages using Relion 2D ML

Examples:


"""

# Program Options
parser = EMArgumentParser(usage, version=EMANVERSION)
parser.add_header(name="relion2dheader", help='Options below this label are specific to e2refinetorelion2d', title="### e2refinetorelion options ###", row=0, col=0, rowspan=1, colspan=3)
parser.add_pos_argument(name="set_name",help="The set name of the set of particles.", default="", guitype='filebox', browser="EMBrowserWidget(withmodal=True,multiselect=False)", row=1, col=0, rowspan=1, colspan=2)
parser.add_argument("--numiter", type=int, help="# of iterations to refine", default=25, guitype='intbox',row=2, col=0, rowspan=1, colspan=1 )
parser.add_argument("--numclasses", type=int, help="# of classes", default=8, guitype='intbox',row=2, col=1, rowspan=1, colspan=1 )
parser.add_argument("--ctfcorrect", action="store_true", help="(T/F)Do CTF Correction?", default=False, guitype='boolbox', row=3, col=0, rowspan=1, colspan=1)
parser.add_argument("--onlyflipphase", action="store_true", help="(T/F)Only flip phases?", default=False, guitype='boolbox',row=5, col=0, rowspan=1, colspan=1 )
parser.add_argument("--dataphaseflipped", action="store_true", help="(T/F)Has the data been phase flipped already?", default=False, guitype='boolbox',row=5, col=1, rowspan=1, colspan=2 )
parser.add_argument("--ignoretofirstpeak", action="store_true", help="(T/F)Ignore CTF's until the first peak?", default=False, guitype='boolbox',row=4, col=0, rowspan=1, colspan=2 )
parser.add_argument("--intensitycorrection", action="store_true", help="(T/F)Do intensity correction?", default=False, guitype='boolbox',row=3, col=1, rowspan=1, colspan=2 )
parser.add_argument("--maskrefstructure", action="store_true", help="(T/F)Mask reference structures?", default=False, guitype='boolbox',row=4, col=1, rowspan=1, colspan=1 )
parser.add_argument("--splithalves", action="store_true", help="(T/F)Split data into random halves?", default=False, guitype='boolbox',row=6, col=0, rowspan=1, colspan=1 )
parser.add_argument("--joinhalves", action="store_true", help="(T/F)Join random halves?", default=False, guitype='boolbox',row=6, col=1, rowspan=1, colspan=2 )
parser.add_argument("--regparam", type=float, help="Regularization Parameter T (weights experimental data vs. prior", default=1.0, guitype='floatbox',row=7, col=0, rowspan=1, colspan=1 )
parser.add_argument("--imagemaskd", type=float, help="Diameter of the image mask", default=-1, guitype='floatbox',row=7, col=1, rowspan=1, colspan=1 )
parser.add_argument("--inplaneang", type=float, help="In-plane angular sampling", default=None, guitype='floatbox',row=8, col=0, rowspan=1, colspan=1 )
parser.add_argument("--oversampling",help="Oversampling order",default='0', guitype='combobox', choicelist='0,1,2', row=9, col=0, rowspan=1, colspan=1 )
parser.add_argument("--offsetrange", type=float, help="Offset search range (pix)", default=10.0, guitype='floatbox',row=8, col=1, rowspan=1, colspan=1 ) 
parser.add_argument("--offsetstep", type=float, help="Offset search step (pix)", default=2.0, guitype='floatbox',row=8, col=2, rowspan=1, colspan=1 )
parser.add_argument("--threads", type=int, help="# of threads", default=1, guitype='intbox',row=11, col=0, rowspan=1, colspan=1 )
parser.add_argument("--echo",action="store_true", default=False, help="Echo Relion Command to terminal only")
parser.add_argument("--defocus", type=float, help="Defocus in A", default=10000, guitype='floatbox',row=9, col=1, rowspan=1, colspan=1 )
parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
parser.add_argument("--voltage",type=int, default=None, help="Voltage of the Microscope (kV)", guitype='intbox', row=10, col=0, rowspan=1, colspan=1)
parser.add_argument("--cs", type=float, default=None, help="Spherical Aberration", guitype='floatbox', row=10, col=1, rowspan=1, colspan=1)
parser.add_argument("--apix", type=float, default=None, help="Angstrom per pixel", guitype='floatbox', row=10, col=2, rowspan=1, colspan=1)
optionList = pyemtbx.options.get_optionlist(sys.argv[1:])
(options, args) = parser.parse_args()


#Check for basic usage
if len(args) != 1:
   print "usage:" + usage
   print "Please run'" + progname + " -h' for detailed options"
   sys.exit(1)


#E2n = E2init(sys.argv,options.ppid)
set_name = args[0]

E2n = E2init(sys.argv,options.ppid)

# Create the E2RLN directory structure if it does not exist
i = 1
found = 1
while found == 1:
	if i < 10:
		rln_run = '0' + str(i)
	else:
		rln_run = str(i)
	found = os.path.exists("relion2d_" + rln_run)
	i = i+1
E2RLN = "relion2d_" + rln_run
os.mkdir(E2RLN)

verbosity = '0'
for option1 in optionList:
	if option1 == "v":
		verbosity=str(options.v)

header= EMData(set_name,0,True)
h_keys=header.get_attr_dict().keys()
nx,ny,nz=header['nx'],header['ny'],header['nz']
num_images = EMUtil.get_image_count(set_name)
mrc = False
if num_images > nz:
        num_ptcl = num_images
else:
        num_ptcl = nz
        mrc = True
project_db = db_open_dict("bdb:.#project")

if options.apix == None:
        if header.get_attr_dict().__contains__('apix_x'):
                apix = header['apix_x']
        elif project_db.__contains__('global.apix'):
                apix = project_db['global.apix']
        else:
                print "An Angstrom per pixel was not found in the project database, the images themselves, and was not provided via a command line option. Please provide another input file"
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
                print "A spherical aberration value was not found in the project database, the images themselves, and was not provided via a command line option. Please provide another input file"
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
                print "A microscope voltage was not found in the project database, the images themselves, and was not provided via a command line option. Please provide another input file"
                print "Exiting e2refinetorelion3d"
                shutil.rmtree(E2RLN)
                exit(-1)
else:
        voltage = options.voltage

project_db = db_open_dict("bdb:.#project")

# Create the particle stack files (.mrcs files) and the .star file that RELION needs as an inputs.
set_orig = set_name
i = 0
old_src = EMData(set_name,0)['data_source']
if mrc:
        s =  "e2proc2d.py " + set_orig + " " + E2RLN + "/ptcl_stack.hdf --threed2twod --verbose=0"
else:
        s =  "e2proc2d.py " + set_orig + " " + E2RLN + "/ptcl_stack.hdf --verbose=0"
call(s,shell=True)

ctf_corr = 0
for option1 in optionList:
	if option1 == "ctfcorrect":
		ctf_corr = 1
if ctf_corr == 1:
	s = "relion_star_loopheader rlnImageName rlnMicrographName rlnDefocusU rlnDefocusV rlnDefocusAngle rlnVoltage rlnSphericalAberration rlnAmplitudeContrast > " + E2RLN + "/all_images.star"
	DEF1 = DEF2 = str(options.defocus)
else:
	s = "relion_star_loopheader rlnImageName rlnMicrographName rlnVoltage rlnAmplitudeContrast > " + E2RLN + "/all_images.star"

call(s,shell=True)
print "Converting EMAN2 Files to Formats Compatible with RELION"

for k in range(num_images):
	src = EMData(set_name,k)['data_source']
	if src != old_src:
		s = "e2proc2d.py " + E2RLN + "/ptcl_stack.hdf" + " " + E2RLN + "/" + old_src.replace('bdb:particles#','') + ".hdf --first=" + str(i) + " --last=" + str(k-1) + " --verbose=" + verbosity
		call(s, shell=True)
		if (k-i-1) == 0:
			s = "e2proc2d.py " + E2RLN + "/" + old_src.replace('bdb:particles#','') + ".hdf " + E2RLN + "/" + old_src.replace('bdb:particles#','') + ".mrc --verbose=" + verbosity
			call(s,shell=True)
		else:
			s = "e2proc2d.py " + E2RLN + "/" + old_src.replace('bdb:particles#','') + ".hdf " + E2RLN + "/" + old_src.replace('bdb:particles#','') + ".mrc --twod2threed --verbose=" + verbosity 
			call(s, shell=True)
		s1 = E2RLN + "/" + old_src.replace('bdb:particles#','') + ".mrc"
		s2 = s1 + "s"
		shutil.move(s1, s2)
		if ctf_corr == 1:
			s = "relion_star_datablock_stack " +  str(k-i) + " " +  E2RLN + "/" + old_src.replace('bdb:particles#','') + ".mrcs " + E2RLN + "/" + old_src.replace('bdb:particles#','') + ".mrcs " + str(voltage) + " " + amplitude_contrast + "  >> " + E2RLN + "/all_images.star" 
		else:
			s = "relion_star_datablock_stack " +  str(k-i) + " " +  E2RLN + "/" + old_src.replace('bdb:particles#','') + ".mrcs " + E2RLN + "/" + old_src.replace('bdb:particles#','') + ".mrcs " + str(options.defocus) + " " + str(options.defocus) + " 0 " +str(voltage) + " " + str(options.cs) + " " + amplitude_contrast + "  >> " + E2RLN + "/all_images.star" 
		call(s,shell=True)
		s = "rm " + E2RLN + "/" + old_src.replace('bdb:particles#','') + ".hdf" 
		call(s,shell=True)
		i = k
		old_src = src
	if (k+1) == num_images:
		s = "e2proc2d.py " + E2RLN + "/ptcl_stack.hdf" + " " + E2RLN + "/" + src.replace('bdb:particles#','') + ".hdf --first=" + str(k) + " --last=" + str(k) + " --verbose=" + verbosity
		call(s, shell=True)
		s = "e2proc2d.py " + E2RLN + "/" + src.replace('bdb:particles#','') + ".hdf " + E2RLN + "/" + src.replace('bdb:particles#','') + ".mrc --verbose=" + verbosity
		call(s,shell=True)
		s1 = E2RLN + "/" + src.replace('bdb:particles#','') + ".mrc"
		s2 = s1 + "s"
		shutil.move(s1, s2)
		if ctf_corr == 1:
			s = "relion_star_datablock_stack 1 "+  E2RLN + "/" + src.replace('bdb:particles#','') + ".mrcs " + E2RLN + "/" + src.replace('bdb:particles#','') + ".mrcs "  + str(options.defocus) + " " + str(options.defocus) + " 0 " + str(voltage) + " " + str(options.cs) + " " + amplitude_contrast + "  >> " + E2RLN + "/all_images.star" 
		else:
			s = "relion_star_datablock_stack 1 "+  E2RLN + "/" + src.replace('bdb:particles#','') + ".mrcs " + E2RLN + "/" + src.replace('bdb:particles#','') + ".mrcs " + str(voltage) + " " + amplitude_contrast + "  >> " + E2RLN + "/all_images.star" 
		call(s,shell=True)
		s = "rm " + E2RLN + "/" + src.replace('bdb:particles#','') + ".hdf" 
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


echo=False
#Parse the options and create the command to run Relion
s = RELION_RUN + "--i " + E2RLN + "/all_images.star --o " + RUNDIR + "/" + E2RLN + " --angpix " + str(apix)
for option1 in optionList:
	if option1 == "echo":
		echo=True
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
	elif option1 == "numclasses":
		s = s + " --K " + str(options.numclasses)
	elif option1 == "oversampling":
		s = s + " --oversampling " + str(options.oversampling)
	elif option1 == "inplaneang":
		s = s + " --psi_step " + str((2 ** int(options.oversampling)) * options.inplaneang)
	elif option1 == "offsetrange":
		s = s + " --offset_range " + str(options.offsetrange)
	elif option1 == "offsetstep":
		s = s + " --offset_step " + str(options.offsetstep)
	elif option1 == "threads":
		s = s + " --j " + str(options.threads)
print "Running Relion 2D Class Averaging"
print "************** Relion Command *******************"
print s
print "*************************************************"
if not(echo):
	call(s,shell=True)
	print "Relion 2D Class Averaging Complete"


	# Clean the Directory up
	print "Cleaning up Directory Structure"
	os.mkdir(E2RLN + "/input_files")
	s = "mv " + E2RLN + "/*.mrcs " + E2RLN + "/input_files"
	call(s,shell=True)
	s1 = E2RLN + "/all_images.star"
	s2 = E2RLN + "/input_files/all_images.star" 
	shutil.move(s1,s2)
	print "Directory Structure Cleaned"




	#Move the classes created by Relion into an eman2 style format
	os.mkdir(E2RLN + "/tmp")
	if (options.numiter) < 10:
		s = "e2proc2d.py " + RUNDIR + "/" + E2RLN + "_it00" + str(options.numiter) + "_classes.mrcs " + E2RLN + "/tmp/classes.mrc --unstacking --verbose=" + verbosity
		data_file =  RUNDIR + "/" + E2RLN + "_it00" + str(options.numiter) + "_data.star"
		call(s,shell=True)
	elif (options.numiter) < 100:
		s = "e2proc2d.py " + RUNDIR + "/" + E2RLN + "_it0" + str(options.numiter) + "_classes.mrcs " + E2RLN + "/tmp/classes.mrc --unstacking --verbose=" + verbosity
		data_file =  RUNDIR + "/" + E2RLN + "_it0" + str(options.numiter) + "_data.star"
		call(s,shell=True)
	elif (options.numiter) < 1000:
		s = "e2proc2d.py " + RUNDIR + "/" + E2RLN + "_it" + str(options.numiter) + "_classes.mrcs " + E2RLN + "/tmp/classes.mrc --unstacking --verbose=" + verbosity
		data_file =  RUNDIR + "/" + E2RLN + "_it" + str(options.numiter) + "_data.star"
		call(s,shell=True)
	db_classes = db_open_dict("bdb:" + E2RLN + "#classes_" + rln_run)
	
	for i in range(options.k):
		a = EMData(E2RLN + "/tmp/classes-" + str(i+1) + ".mrc")
		a['apix_x'] = a['apix_y'] = a['apix_z'] = apix
		a['origin_z'] = 0
		a['class_ptcl_src'] = set_orig
		db_classes[i] = a


	d = defaultdict(list)
	f = open(data_file, 'r')
	col = 0
	row = 0
	for line in f.readlines():
		col = col + 1
		row = row + 1
		splits = line.split()
		if len(splits) > 0:
			if splits[0] == "loop_":
				col = 0
			if splits[0] == "_rlnClassNumber":
				class_col = col
		if line[0] == '_':
			row = 0
		if len(splits) > 1:
			d[int(splits[class_col-1])].append(row)

	f.close()
	db_classes = db_open_dict("bdb:" + E2RLN + "#classes_" + rln_run)
	for key in d.keys():
		a = db_classes[int(key)-1]
		a['class_ptcl_idxs']=d[key]
		a['ptcl_repr']=len(d[key])
		db_classes[int(key)-1] = a

