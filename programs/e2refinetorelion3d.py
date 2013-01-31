#!/usr/bin/env python
# Author: Stephen Murray (scmurray@bcm.edu), 12/12/11
# Copyright (c) 2000-2011 Baylor College of Medicine

# Official copyright notice. EMAN2 is distributed under a joint GPL/BSD license. Please copy
# this statement from one of the other programs. You must agree to use this license if your
# code is distributed with EMAN2. While you may use your own institution for the copyright notice
# the terms of the GPL/BSD license permit us to redistribute it.

#****************UPDATED for RELION 1.1 RELEASE on 12/4/12******************

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
# Required Program Options and Parameters (GUI and Command Line)
parser = EMArgumentParser(usage, version=EMANVERSION)
parser.add_header(name="relion3dheader", help="Options below this label are specific to e2refinetorelion3d", title="   ### e2refinetorelion3d.py Options ###", row=0, col=0, rowspan=1, colspan=3)
#I/O Options
parser.add_header(name="io", help="Options in this section pertain to I/O", title="---I/O Options---", row=1, col=0, rowspan=1, colspan=1) 
parser.add_pos_argument(name="set_name",help="The set name of the set of particles.", default="", guitype='filebox', browser="EMBrowserWidget(withmodal=True,multiselect=False)", row=2, col=0, rowspan=2, colspan=2)
parser.add_pos_argument(name="refmap", type=str, help="Reference Map", default="", guitype='filebox', browser="EMBrowserWidget(withmodal=True,multiselect=False)", row=4, col=0, rowspan=2, colspan=2)
parser.add_argument("--randomizemodel", help="Optionally randomize the phases of the initial model to this resolution (in Angstroms)", default=0,  guitype='floatbox', row=6, col=0, rowspan=1, colspan=1)
parser.add_argument("--greyscale", action="store_true", help="Is the reference map in greyscale?", default=False, guitype='boolbox', row=6, col=1, rowspan=1, colspan=1)
parser.add_argument("--symmgroup", help="Symmetry group", guitype='combobox', default='C',choicelist="""'C','D','T','O','I'""", row=7, col=0, rowspan=1, colspan=1)
parser.add_argument("--symmnumber",type=int,help="Symmetry number",default=1, guitype='intbox', row=7, col=1, rowspan=1, colspan=1)
#CTF Options
parser.add_header(name="ctf", help="Options in this section pertain to CTF parameters", title="---CTF Options---",row=8, col=0, rowspan=1, colspan=3)
parser.add_argument("--ctfcorrect", action="store_true", help="(T/F)Do CTF Correction?", default=False, guitype='boolbox', row=9, col=0, rowspan=1, colspan=1)
parser.add_argument("--dataphaseflipped", action="store_true", help="(T/F)Has the data been phase flipped already?", default=False, guitype='boolbox', row=10, col=0, rowspan=1, colspan=2)
parser.add_argument("--ignoretofirstpeak", action="store_true", help="(T/F)Ignore CTF's until the first peak?", default=False, guitype='boolbox', row=11, col=0, rowspan=1, colspan=2)
#Optimisation Options
parser.add_header(name="optimisation", help="Options in this section pertain to Optimisation parameters", title="---Optimisation Options---", row=12, col=0, rowspan=1, colspan=3)
parser.add_argument("--lowpass",type=float,help="Initial low-pass filter (Ang)",default=60.0, guitype='floatbox', row=13, col=0, rowspan=1, colspan=1)
parser.add_argument("--imagemaskd", type=float, help="Diameter of the image mask", default=-1, guitype='floatbox', row=13, col=1, rowspan=1, colspan=1)
parser.add_argument("--solventmask",type=str, help="Location of the mask to be used", guitype='filebox',default="", browser="EMBrowserWidget(withmodal=True,multiselect=False)", filecheck=False,row=14, col=0, rowspan=2, colspan=2)
#Sampling Options
parser.add_header(name="sampling", help="Options in this section pertain to Sampling parameters", title="---Sampling Options---", row=16, col=0, rowspan=1, colspan=3)
parser.add_argument("--healpix", type=str, default=7.5, help="Angular Sampling Interval (Degrees)", guitype='combobox', choicelist='30,15,7.5,3.7,1.8,0.9,0.5,0.2,0.1', row=17, col=0, rowspan=1, colspan=1)
parser.add_argument("--auto_healpix", type=str, default=1.8, help="Local angular search value", guitype='combobox', choicelist='30,15,7.5,3.7,1.8,0.9,0.5,0.2,0.1', row=17, col=1, rowspan=1, colspan=1)
parser.add_argument("--offsetrange", type=float, help="Offset search range (pix)", default=10.0, guitype='floatbox', row=18, col=0, rowspan=1, colspan=1) 
parser.add_argument("--offsetstep", type=float, help="Offset search step (pix)", default=2.0, guitype='floatbox', row=18, col=1, rowspan=1, colspan=1)
#Run Options
parser.add_header(name="running", help="Options in this section pertain to running RELION itself", title="---Run Options---", row=19, col=0, rowspan=1, colspan=3)
parser.add_argument("--threads", type=int, help="# of threads", default=1, guitype='intbox', row=20, col=0, rowspan=1, colspan=1)
parser.add_argument("--maxmemory", type=float, help="Maximum memory (in GB) available for each thread", guitype='floatbox', row=20, col=1, rowspan=1, colspan=1)


# Command line only parameters (or expert parameters in hte GUI)
parser.add_header(name="Expert Options", help="Options in this section are for expert users", title="---Expert Options---", expert=True, row=21, col=0, rowspan=1, colspan=3)
parser.add_argument("--amplitudecontrast", type=float, help="Amplitude Contrast value for the micrographs", default=0.07, guitype='floatbox',row=22, col=0, rowspan=1, colspan=1, expert=True)
parser.add_argument("--intensitycorrection", action="store_true", help="(T/F)Do intensity correction?", default=False, guitype='boolbox', row=22, col=1, rowspan=1, colspan=1, expert=True)
parser.add_argument("--print_symmetry", action="store_true",help="Print all symmetry transformation matrices, and exit", default=False, guitype='boolbox', row=23, col=0, rowspan=1, colspan=1, expert=True)
parser.add_argument("--nearest_neighbor", action="store_true", help="Perform nearest-neighbor instead of linear Fourier-space interpolation", default=False, guitype='boolbox', row=23, col=1, rowspan=1, colspan=1, expert=True)
parser.add_argument("--pad", type=int,help="Padding factor",  default=1, guitype='combobox', choicelist='1,2,3', row=24, col=0, rowspan=1, colspan=1, expert=True)
parser.add_argument("--oversampling",help="Oversampling order", default=0, guitype='combobox', choicelist='0,1,2,3,4', row=24, col=1, rowspan=1, colspan=1, expert=True)
parser.add_argument("--limit_tilt",type=int, help="Limited tilt angle: positive for keeping side views, negative for keeping top views", default=-91, guitype='intbox', row=25, col=0, rowspan=1, colspan=1, expert=True)
parser.add_argument("--verbosity", type=int, help="Set the level of verbosity for the code", default=1, guitype='combobox', choicelist='0,1,2,3,4,5,6,7,8,9', row=25, col=1, rowspan=1, colspan=1, expert=True)
parser.add_argument("--onlyflipphase", action="store_true", help="(T/F)Only flip phases?", default=False, guitype='boolbox', row=26, col=0, rowspan=1, colspan=1, expert=True)

parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
#parser.add_argument("--inplaneang", type=float, help="In-plane angular sampling", default=row=26, col=0, rowspan=1, colspan=2, expert=True)
#parser.add_argument("--refctfcorrected", action="store_true", help="Has the reference map been ctf corrected?")
#parser.add_argument("--maskrefstructure", action="store_true", help="(T/F)Mask reference structures?", default=False, guitype='boolbox', row=5, col=0, rowspan=1, colspan=2)
#parser.add_argument("--splithalves", action="store_true", help="(T/F)Split data into random halves?")
#parser.add_argument("--joinhalves", action="store_true", help="(T/F)Join random halves?")
#parser.add_argument("--regparam", type=float, help="Regularization Parameter T (weights experimental data vs. prior", default=1.0, guitype='floatbox', row=14, col=0, rowspan=1, colspan=1)
#parser.add_argument("--local", type=float,help="Perform local angular search at what range?",default=None, guitype='floatbox', row=15, col=0, rowspan=1, colspan=1)
#parser.add_argument("--nmpi", type=int, help="# of MPI procs to use", default=1, guitype='intbox', row=17, col=1, rowspan=1, colspan=1)
#parser.add_argument("--autosample", action="store_true", help="Perform automated orientational sampling?", default=False, guitype='boolbox', row=20, col=0, rowspan=1, colspan=1)
#parser.add_argument("--numiter", type=int, help="# of iterations", default=25, guitype='intbox', row=8, col=0, rowspan=1, colspan=1)
#parser.add_argument("--numclasses", type=int, help="# of classes", default=8, guitype='intbox', row=8, col=1, rowspan=1, colspan=1)

optionList = pyemtbx.options.get_optionlist(sys.argv[1:])
(options, args) = parser.parse_args()

#Check for basic usage
if len(args) != 2:
   print "usage:" + usage
   print "Please run'" + progname + " -h' for detailed options"
   sys.exit(1)

for option1 in optionList:
	if option1 == "ctfcorrect":
		ctf_corr = 1
	if option1 == "amplitudecontrast":
		amplitude_contrast = str(options.amplitudecontrast)
	#if option1 == "auto_healpix":
	#	if options.autosample == False:
	#		print "The --auto_healpix option can only be used when the --autosample option is checked. Please rectify this and re-run if you wish to use the auto sampling"
			

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
os.mkdir(E2RLN + "/stacks")
#E2n = E2init(sys.argv,options.ppid)
set_name, refmap = args[0], args[1]
num_images = EMUtil.get_image_count(set_name)
s1 = "e2proc3d.py " + refmap + " " + E2RLN + "/3DRefMap.mrc --process=normalize.edgemean --verbose="+str(options.verbosity)
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


if header.get_attr_dict().__contains__('apix') and header.get_attr_dict().__contains__('apix_x'):
	if header['apix_x'] == 1:
		apix=header['apix']
		print """***An "apix" value was found in the header. This is an old style and apix values should be correctly stored in apix_x, apix_y, apix_z. The value in apix is """ + str(apix) + """. Be aware this may not be the value you intended***"""
	else:
		apix=['apix_x']
elif header.get_attr_dict().__contains__('apix_x'):
	print """***An "apix" value was found in the header. This is an old style and apix values should be correctly stored in apix_x, apix_y, apix_z. The value in apix is """ + str(apix) + """. Be aware this may not be the value you intended***"""
	apix = header['apix_x']
elif header.get_attr_dict().__contains__('apix'):
	apix = header['apix']
	print """***An "apix" value was found in the header. This is an old style and apix values should be correctly stored in apix_x, apix_y, apix_z. The value in apix is """ + str(apix) + """. Be aware this may not be the value you intended***"""
elif project_db.__contains__('global.apix'):
	apix = project_db['global.apix']
else:
	print "An Angstrom per pixel was not found in the project database nor in the images themselves. Please ensure it exists in the data header"
	print "Exiting e2refinetorelion3d"
	shutil.rmtree(E2RLN)
	exit(-1)

if header.get_attr_dict().__contains__('cs'):
	cs = header['cs']
elif project_db.__contains__('global.microscope_cs'):
	cs = project_db['global.microscope_cs']
else:
	print "A spherical aberration value was not found in the project database nor in the images themselves. Please ensure it exists in the data header"
	print "Exiting e2refinetorelion3d"
	shutil.rmtree(E2RLN)
	exit(-1)
	

if header.get_attr_dict().__contains__('voltage'):
	voltage = header['voltage']
elif project_db.__contains__('global.microscope_voltage'):
	voltage = project_db['global.microscope_voltage']
else:
	print "A microscope voltage was not found in the project database nor in the images themselves. Please ensure it exists in the data header"
	print "Exiting e2refinetorelion3d"
	shutil.rmtree(E2RLN)
	exit(-1)


## Create the particle stack files (.mrcs files) and the .star file that RELION needs as an inputs.
set_orig = set_name
#set_db = db_open_dict(set_orig)
i = 0
old_src = EMData(set_name,0)['data_source']
s =  "e2proc2d.py " + set_orig + " " + E2RLN + "/ptcl_stack.hdf --verbose="+str(options.verbosity)
call(s,shell=True)
ctf_corr = 0
for option1 in optionList:
	if option1 == "ctfcorrect":
		ctf_corr = 1
	if option1 == "amplitudecontrast":
		amplitude_contrast = str(options.amplitudecontrast)
	if option1 == "autosample":
		asample = True
if ctf_corr == 1:

	dblist = db_list_dicts("bdb:./sets")	
	for db in dblist:
		db_src=set_name.replace("bdb:./sets#",'').replace("_original_data",'')
		if not db.find(db_src):
			db_set=EMData("bdb:./sets#"+db,0,True)
			if db_set.get_attr_dict().__contains__('ctf') and (EMUtil.get_image_count("bdb:sets#"+db) == num_images):
				ctf_value=True
#				defocus = db_set['ctf'].to_dict()['defocus']*1000
				break
	print "CTF information being pulled from: " + db
	s = "relion_star_loopheader rlnImageName rlnMicrographName rlnDefocusU rlnDefocusV rlnDefocusAngle rlnVoltage rlnSphericalAberration rlnAmplitudeContrast > " + E2RLN + "/all_images.star"
	defocus1 = defocus2 = 0
else:
	s = "relion_star_loopheader rlnImageName rlnMicrographName rlnVoltage rlnAmplitudeContrast > " + E2RLN + "/all_images.star"
call(s,shell=True)

print "Converting EMAN2 Files to Formats Compatible with RELION"
temp = EMData(set_name,0)
#print '# \tMicrograph \t\tDefocus \tVoltage CS \tAPIX'
for k in range(num_ptcl):
#	print k, '\t',temp['data_source'].split('?')[0].split('#')[1],'\t', temp['defocus'],'\t', temp['voltage'],'\t', temp['cs'],'\t', '1.8'	
	src = EMData(set_name,k)['data_source']
	if (src != old_src):
		temp=EMData("bdb:./sets#"+db,k-1)
		s = "e2proc2d.py " + E2RLN + "/ptcl_stack.hdf" + " " + E2RLN + "/" + old_src.split('?')[0].replace('bdb:particles#','') + ".hdf --verbose="+str(options.verbosity)+" --step=" + str(i) + ",1 --last=" + str(k-1)
		call(s, shell=True)
		if (k-i-1) == 0:
			s = "e2proc2d.py " + E2RLN + "/" + old_src.split('?')[0].replace('bdb:particles#','') + ".hdf " + E2RLN + "/" + old_src.split('?')[0].replace('bdb:particles#','') + ".mrc --verbose="+str(options.verbosity) + " --process=normalize.edgemean"
                        call(s, shell=True)
		else:
			s = "e2proc2d.py " + E2RLN + "/" + old_src.split('?')[0].replace('bdb:particles#','') + ".hdf " + E2RLN + "/" + old_src.split('?')[0].replace('bdb:particles#','') + ".mrc --verbose=" + str(options.verbosity) + " --process=normalize.edgemean --twod2threed" 
			call(s, shell=True)
		s1 = E2RLN + "/" + old_src.split('?')[0].replace('bdb:particles#','') + ".mrc"
		stemp="e2proc3d.py " + s1 + " " + s1 + " --process=normalize"
		call(stemp, shell=True)
		s2 = E2RLN + "/stacks/" + old_src.split('?')[0].replace('bdb:particles#','') + ".mrcs"
		shutil.move(s1, s2)
		if ctf_corr == 1:
			defocus1 = defocus2 = str(temp['ctf'].to_dict()['defocus']*1000)
			s = "relion_star_datablock_stack " + str(k-i) + " " + E2RLN + "/stacks/" + old_src.split('?')[0].replace('bdb:particles#','') + ".mrcs " + E2RLN + "/stacks/" + old_src.split('?')[0].replace('bdb:particles#','') + ".mrcs " + str(defocus1) + " " + str(defocus2) + " 0 " + str(voltage) + " " + str(cs) + " " + str(amplitude_contrast) + " >> " + E2RLN + "/all_images.star"
		else:
			s = "relion_star_datablock_stack " + str(k-i) + " " + E2RLN + "/stacks/" + old_src.split('?')[0].replace('bdb:particles#','') + ".mrcs " + E2RLN + "/stacks/" + old_src.split('?')[0].replace('bdb:particles#','') + ".mrcs " + str(voltage) + " " + str(amplitude_contrast) + "  >> " + E2RLN + "/all_images.star" 
		call(s,shell=True)
		s = "rm " + E2RLN + "/" + old_src.split('?')[0].replace('bdb:particles#','') + ".hdf" 
		call(s,shell=True)
		i = k
		old_src = src
	elif (k+1) == num_ptcl:
		diff = k-i
		temp=EMData("bdb:./sets#"+db,k)
		s = "e2proc2d.py " + E2RLN + "/ptcl_stack.hdf" + " " + E2RLN + "/" + src.split('?')[0].replace('bdb:particles#','') + ".hdf --verbose="+str(options.verbosity)+" --step=" + str(i) + ",1 --last=" + str(k)
		call(s, shell=True)
		if (k-i-1) == 0:
			s = "e2proc2d.py " + E2RLN + "/" + src.split('?')[0].replace('bdb:particles#','') + ".hdf " + E2RLN + "/" + src.split('?')[0].replace('bdb:particles#','') + ".mrc --verbose="+str(options.verbosity) + " --process=normalize.edgemean"
		else:
			s = "e2proc2d.py " + E2RLN + "/" + src.split('?')[0].replace('bdb:particles#','') + ".hdf " + E2RLN + "/" + src.split('?')[0].replace('bdb:particles#','') + ".mrc --verbose="+str(options.verbosity) + " --process=normalize.edgemean --twod2threed"
		call(s, shell=True)
		s1 = E2RLN + "/" + src.split('?')[0].replace('bdb:particles#','') + ".mrc"
		stemp="e2proc3d.py " + s1 + " " + s1 + " --process=normalize"
		call(stemp, shell=True)
		s2 = E2RLN + "/stacks/" + src.split('?')[0].replace('bdb:particles#','') + ".mrcs"
		shutil.move(s1, s2)
		if ctf_corr == 1:
			defocus1 = defocus2 = str(temp['ctf'].to_dict()['defocus']*1000)
			s = "relion_star_datablock_stack " + str(k-i+1) + " " + E2RLN + "/stacks/" + old_src.split('?')[0].replace('bdb:particles#','') + ".mrcs " + E2RLN + "/stacks/" + old_src.split('?')[0].replace('bdb:particles#','') + ".mrcs " + str(defocus1) + " " + str(defocus2) + " 0 " + str(voltage) + " " + str(cs) + " " + str(amplitude_contrast) + " >> " + E2RLN + "/all_images.star"
		else:
			s = "relion_star_datablock_stack " + str(k-i+1) + " " +  E2RLN + "/stacks/" + old_src.split('?')[0].replace('bdb:particles#','') + ".mrcs " + E2RLN + "/stacks/" + old_src.split('?')[0].replace('bdb:particles#','') + ".mrcs " + str(voltage) + " " + amplitude_contrast + "  >> " + E2RLN + "/all_images.star" 
		call(s,shell=True)
		s = "rm " + E2RLN + "/" + src.split('?')[0].replace('bdb:particles#','') + ".hdf" 
		call(s,shell=True)
		i = k
		old_src = src
		

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
s = RELION_RUN + "--i " + E2RLN + "/all_images.star --o " + RUNDIR + "/" + E2RLN + " --angpix " + str(apix) + " --auto_sampling --split_random_halves --norm"
grey = oversample = 0

for option1 in optionList:
	if option1 == "imagemaskd":
		s = s + " --particle_diameter " + str(options.imagemaskd)
	#elif option1 == "autosample":
	#	s = s + " --auto_sampling"
	elif option1 == "auto_healpix":
		if options.auto_healpix == '30':
			s = s + " --auto_local_healpix_order 0"
		elif options.auto_healpix =='15':
			s = s + " --auto_local_healpix_order 1"
		elif options.auto_healpix =='7.5':
			s = s + " --auto_local_healpix_order 2"
		elif options.auto_healpix =='3.7':
			s = s + " --auto_local_healpix_order 3"
		elif options.auto_healpix =='1.8':
			s = s + " --auto_local_healpix_order 4"
		elif options.auto_healpix =='0.9':
			s = s + " --auto_local_healpix_order 5"
		elif options.auto_healpix =='0.5':
			s = s + " --auto_local_healpix_order 6"
		elif options.auto_healpix =='0.2':
			s = s + " --auto_local_healpix_order 7"
		elif options.auto_healpix =='0.1':
			s = s + " --auto_local_healpix_order 8"
		else:
			print "Invalid angular sampling interval (--auto_healpix). Defaulting to 1.8 degrees"
			s = s + " --auto_local_healpix_order 4"
#	elif option1 == "numiter":
#		s = s + " --iter " + str(options.numiter)
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
	elif option1 == "verbosity":
		if options.verbosity == 0:
			s = s + " --verb 0"
		else:
			s = s + " --verb 1"
	elif option1 == "limit_tilt":
		s = s + " --limit_tilt " + str(options.limit_tilt)
	elif option1 == "print_symmetry":
		s = s + " --print_symmetry_ops"
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
#	elif option1 == "numclasses":
#		s = s + " --K " + str(options.numclasses)
	elif option1 == "oversampling":
		s = s + " --oversampling " + str(options.oversampling)
		oversample = 1
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
	elif option1 == "randomizemodel":
		if float(options.randomizemodel) != 0:
			s1 = "e2proc3d.py " + E2RLN + "/3DRefMap.mrc " + E2RLN + "/3DRefMap.mrc --process=filter.lowpass.randomphase:apix=" + str(apix) + ":cutoff_freq=" + str(1/float(options.randomizemodel))
			call(s1,shell=True)
	elif option1 == "pad":
		s = s + " --pad " + str(options.pad)
	elif option1 == "lowpass":
		s = s + " --ini_high " + str(options.lowpass)
#	elif option1 == "local":
#		s = s + " --sigma_ang " + str(options.local / 3)
	elif option1 == "nearest_neighbor":
		s = s + " --NN"
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
			print "Invalid angular sampling interval (--healpix). Defaulting to 7.5 degrees"
			s = s + " --healpix_order 2"

s = s + " --sym " + str(options.symmgroup) + str(options.symmnumber) + " --ref " + E2RLN + "/3DRefMap.mrc"

if grey != 1:
	s = s + " --firstiter_cc"
if oversample:
	s = s + " --oversampling 1"

print "Generating qsub file"
###### create files to make the runs ######
qsub_file = current_dir + "/" + RUNDIR + "/qsub.pbs"
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

