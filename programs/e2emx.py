#!/usr/bin/env python
#********************************************************************************
# Author: Stephen Murray (scmurray@bcm.edu), 6/12/13
# Copyright (c) 2000-2013 Baylor College of Medicine

# Official copyright notice. EMAN2 is distributed under a joint GPL/BSD license. Please copy
# this statement from one of the other programs. You must agree to use this license if your
# code is distributed with EMAN2. While you may use your own institution for the copyright notice
# the terms of the GPL/BSD license permit us to redistribute it.
#********************************************************************************


#import block
from EMAN2 import *
import pyemtbx.options
import os
import sys
import shutil
from subprocess import *
import xml.etree.ElementTree

cwd = os.getcwd()
progname = os.path.basename(sys.argv[0])
usage = """ prog [options]

This program will extract the required information from an EMAN2 project and output it in EMX form. 

"""

print "Running e2emx.py"
# Required Program Options and Parameters (GUI and Command Line)
parser = EMArgumentParser(usage, version=EMANVERSION)
parser.add_argument("--export_whole_project", action="store_true", help="This option will create an emx directory, where it will export the eman2 project into EMX format", default=False)
parser.add_argument("--import_box_coordinates", type=str, help="Import box coordinates and corresponding micrographs")
parser.add_argument("--import_ctf", type=str, help="Import ctf information and corresponding micrographs")
parser.add_argument("--import_2d_alignment", type=str, help="Import particles and corresponding transformation information")
parser.add_argument("--refinedefocus",  action="store_true", help="Will use EMAN2 CTF fitting to refine the defocus by SNR optimization (+-0.1 micron from the current values, no astigmatism adjustment)")
parser.add_argument("--refitdefocus",  action="store_true", help="Will use EMAN2 CTF fitting to refit the defocus values (+-0.1 micron, astigmatism unchanged)")
parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higher number means higher level of verboseness")
parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
optionList = pyemtbx.options.get_optionlist(sys.argv[1:])
(options, args) = parser.parse_args()

#Check for basic usage
#if len(args) != 2:
	#print "usage:" + usage
	#print "Please run'" + progname + " -h' for detailed options"
	#sys.exit(1)


for option1 in optionList:
	if option1 == "export_whole_project":
		if not os.path.exists("./emx"):
			print "Creating EMX directory"
			os.mkdir("emx")
		else:
			print "EMX directory already exists"

		f = open("./emx/particles.emx",'w')

		f.write("""<?xml version='1.0' encoding='utf-8'?>
<EMX version="1.0">
<!--
##########################################################################
#               EMX Exchange file 
#               Produced by e2emx.py (EMAN2, version 2.1)
# 
#  This is a EMX file.
#
#  Information on this file format is available at 
#  http://i2pc.cnb.csic.es/emx
##########################################################################
-->\n""")

		#Check to see if the micrographs folder exists
		dir_list = os.listdir(cwd)
		if "micrographs" in dir_list:
			print "-----Writing Micrograph Information-----"
			for image in os.listdir(cwd + "/micrographs"):
				if not image[0] == '.':
					s = "e2proc2d.py " + cwd + "/micrographs/" + image + " " + cwd + "/emx/" + base_name(image) + ".mrc --verbose=0"
					call(s,shell=True)
					f.write("<micrograph fileName=\"" + base_name(image) + ".mrc\">\n")
					micrograph = EMData()
					micrograph.read_image(cwd + "/micrographs/" + image)
					f.write("  <pixelSpacing>\n")
					f.write("    <X unit =\"A/px\">" + str(micrograph['apix_x']) + "</X>\n")
					f.write("    <Y unit =\"A/px\">" + str(micrograph['apix_y']) + "</Y>\n")
					f.write("  </pixelSpacing>\n")
					f.write("</micrograph>\n")
			#f.write("  <acceleratingVoltage unit=\"kV\">" + str(particle['ctf'].to_dict()['voltage']) + "</acceleratingVoltage>\n")
			#f.write("  <amplitudeContrast>" + str(particle['ctf'].to_dict()['ampcont']) + "</amplitudeContrast>\n")
							
			f.write("\n")
		if "particles" in dir_list:
			print "-----Writing Particle Information-----"
			for ptcl_by_micrograph in os.listdir(cwd + "/particles"):
				if ptcl_by_micrograph.find("__") == -1 and ptcl_by_micrograph.find(".") != 0:
					if os.path.exists(cwd + "/particles/" + ptcl_by_micrograph.replace(".hdf",'') + "__ctf_flip.hdf"):
						particle_stack = EMData().read_images(cwd + "/particles/" + ptcl_by_micrograph.replace(".hdf",'') + "__ctf_flip.hdf")
						s1 = "e2proc2d.py " + cwd + "/particles/" + ptcl_by_micrograph.replace(".hdf",'') + "__ctf_flip.hdf emx/" + ptcl_by_micrograph.replace(".hdf.",".mrcs")
						call(s1,shell=True)
					else:
						particle_stack = EMData().read_images(cwd + "/particles/" + ptcl_by_micrograph)
						s1 = "e2proc2d.py " + cwd + "/particles/" + ptcl_by_micrograph + " emx/" + ptcl_by_micrograph.replace(".hdf",".mrcs")
						call(s1,shell=True)
					#num_images = len(particle_stack)
					
					#print num_images
					index = 1
					for particle in particle_stack:
						f.write("<particle fileName=\"" + ptcl_by_micrograph.replace(".hdf",".mrc") +"\" index=\"" + str(index) + "\">\n")
						f.write("  <micrograph fileName=\"" + base_name(ptcl_by_micrograph).split("_")[0] + ".mrc\"/>\n")
						f.write("  <centerCoord>\n")
						f.write("    <X unit=\"px\">" + str(particle['ptcl_source_coord'][0]) + "</X>\n")
						f.write("    <Y unit=\"px\">" + str(particle['ptcl_source_coord'][1]) + "</Y>\n")
						f.write("  </centerCoord>\n")
						f.write("  <boxSize>\n")
						f.write("    <X unit=\"px\">" + str(particle['nx']) + "</X>\n")
						f.write("    <Y unit=\"px\">" + str(particle['ny']) + "</Y>\n")
						f.write("  </boxSize>\n")
						if particle.get_attr_dict().__contains__("ctf"):
							f.write("  <defocusU unit=\"nm\">" + str(particle['ctf'].to_dict()['defocus']*1000) + "</defocusU>\n")
							f.write("  <defocusV unit=\"nm\">" + str(particle['ctf'].to_dict()['defocus']*1000) + "</defocusV>\n")
							f.write("  <defocusUAngle unit=\"deg\">0.0</defocusUAngle>\n")
						f.write("</particle>\n")
						index = index + 1
		f.write("</EMX>")
		f.close()
	elif option1 == "import_ctf":
		found_per_particle = False
		et = xml.etree.ElementTree.parse(options.import_ctf)
		emx = et.getroot()
		micro_dict = {}
		part_dict = {}
		part_list = []
		first_index = 0
		last_part_filename=""
		for item in emx:
			if item.tag == "micrograph":
				micrograph_filename = item.attrib['fileName']
				temp_dict = {}
				for micrograph_attrib in item.attrib:
					if micrograph_attrib == "index":
						micrograph_index = item.attrib['index']
						temp_dict['index']= micrograph_index
					elif micrograph_attrib == "fileName":
						pass
					else:
						print "Unknown tag: " + micrograph_attrib
				for item2 in item:
					if item2.tag == "acceleratingVoltage":
						voltage = item2.text #in kilovolts
						temp_dict['voltage']= voltage
					elif item2.tag == "amplitudeContrast":
						ampcont = item2.text #0->1
						temp_dict['ampcont']= ampcont
					elif item2.tag == "cs":
						cs = item2.text #in mm
						temp_dict['cs']=cs
					elif item2.tag == "defocusU":
						defocus1 = item2.text # in nm
						temp_dict['defocusU']=float (defocus1) /1000.0
					elif item2.tag == "defocusV":
						defocus2 = item2.text # in nm
						temp_dict['defocusV']=float(defocus2) /1000.0
					elif item2.tag == "defocusUAngle":
						defocus_angle = item2.text # in degrees (0->180)
						temp_dict['dfang']=defocus_angle
					elif item2.tag == "fom":
						fom = item2.text # figure-of-merit (0->1)
						temp_dict['fom']=fom
					elif item2.tag == "pixelSpacing":
						for item3 in item2:
							if item3.tag == "X":
								apix_x = float(item3.text)
								temp_dict['apix_x']=apix_x
							elif item3.tag == "Y":
								apix_y = float(item3.text)
								temp_dict['apix_y']=apix_y
							elif item3.tag == "Z":
								apix_z = float(item3.text)
								temp_dict['apix_z']=apix_z
					else:
						print "Unknown tag: " + item2.tag
				micro_dict[micrograph_filename] = temp_dict
				ctf=EMAN2Ctf()
				ctf.from_dict({"defocus":(float(defocus1)+float(defocus2))/2000,"dfang":float(defocus_angle),"dfdiff":abs(float(defocus1)-float(defocus2))/1000,"voltage":float(voltage),"cs":float(cs),"ampcont":float(ampcont),"apix":float(apix_x)})
				jdb = js_open_dict(info_name(micrograph_filename))
				jdb['ctf_frame']=[512,ctf,(256,256),tuple(),5,1]
				jdb.setval("ctf_frame",jdb['ctf_frame'],deferupdate=True)
				micro_dict[micrograph_filename]['stack']=micrograph_filename
			elif item.tag == "particle":
				print "particle"
				temp_dict={}
				foundU=foundV=foundAng=foundapix = False
				for particle_attrib in item.attrib:
					if particle_attrib == "fileName":
						particle_filename = item.attrib['fileName']
						temp_dict['particle_filename'] = particle_filename
					elif particle_attrib == "index":
						particle = item.attrib['index']
						temp_dict['index']=particle
					else:
						print "Unknown tag: " + particle_attrib
				for item2 in item:
					#temp_dict={}
					if item2.tag == "defocusU":
						temp_dict['defocusU'] = float(item2.text) / 1000 # in nm
						foundU = True
					elif item2.tag == "defocusV":
						temp_dict['defocusV'] = float(item2.text) /1000 # in nm
						foundV = True
					elif item2.tag == "defocusUAngle":
						temp_dict['defocusUAngle'] = item2.text # in degrees (0->180)
						foundAng = True
					elif item2.tag == "micrograph":
						for micrograph_attrib in item2.attrib:
							if micrograph_attrib == "fileName":
								particle_micrograph_filename = item2.attrib['fileName']
								temp_dict['particle_micrograph_filename']=particle_micrograph_filename
							elif micrograph_attrib == "index":
								particle_micrograph_index = item2.attrib['index']
								temp_dict['index']=particle_micrograph_index
							else:
								print "Unknown tag: " + micrograph_attrib
					elif item2.tag == "pixelSpacing":
						for item3 in item2:
							if item3.tag == "X":
								apix_x = float(item3.text)
								foundapix = True
							elif item3.tag == "Y":
								foundapix = True
								apix_y = float(item3.text)
							elif item3.tag == "Z":
								foundapix = True
								apix_z = float(item3.text)
							else:
								print "Unknown Tag: " + item3.tag
				part_list.append(temp_dict)
				if particle_micrograph_filename != last_part_filename:
					micro_dict[particle_micrograph_filename]['first_index'] = int(particle) - 1
					last_part_filename = particle_micrograph_filename
					if 'last_index' not in micro_dict[particle_micrograph_filename].keys():
						micro_dict[particle_micrograph_filename]['last_index'] = int(particle)-1
				else:
					micro_dict[particle_micrograph_filename]['last_index'] = int(particle)-1
					found_per_particle=True
				if foundapix:
					micro_dict[particle_micrograph_filename]["apix_x"]=apix_x
				else:
					apix_x=micro_dict[particle_micrograph_filename]['apix_x']
				ctf=EMAN2Ctf()
				ctf.from_dict({"defocus":(float(defocus1)+float(defocus2))/2000,"dfang":float(defocus_angle),"dfdiff":abs(float(defocus1)-float(defocus2))/1000,"voltage":float(voltage),"cs":float(cs),"ampcont":float(ampcont),"apix":float(apix_x)})
				jdb = js_open_dict(info_name(particle_micrograph_filename))
				jdb['ctf']=[512,ctf,(256,256),tuple(),5,1]
				jdb.setval("ctf",jdb['ctf'],deferupdate=True)
				micro_dict[particle_micrograph_filename]['stack']=particle_filename
		if options.refinedefocus : 
			dfopt="--curdefocushint --refinebysnr"
			if options.verbose>0 : print "CTF Refinement"
		elif options.refitdefocus : 
			dfopt="--curdefocushint"
			if options.verbose>0 : print "CTF Refit"
		else: 
			dfopt="--curdefocusfix"
			if options.verbose>0 : print "Computing particle SNRs"
		if not os.path.exists("particles"):
			os.mkdir("particles")
		for item in micro_dict.keys():
			print "e2proc2d.py {} particles/{}_ptcls.hdf --threed2twod --first {} --last {}".format(micro_dict[item]['stack'],base_name(item),micro_dict[item]['first_index'],micro_dict[item]['last_index']) 
			launch_childprocess("e2proc2d.py {} particles/{}_ptcls.hdf --threed2twod --first {} --last {}".format(micro_dict[item]['stack'],base_name(item),micro_dict[item]['first_index'],micro_dict[item]['last_index']))
			launch_childprocess("e2ctf.py particles/{}_ptcls.hdf --voltage {} --cs {} --ac {} --apix {} --autofit --zerook --storeparm --astigmatism {} -v {}".format(base_name(item),micro_dict[item]['voltage'],micro_dict[item]['cs'],micro_dict[item]['ampcont'],micro_dict[item]['apix_x'],dfopt,options.verbose-1))
		if found_per_particle:
			print "Per-particle defocus values or angles found. Please note that we do not support import of this information at the moment. Using the per-micrograph information provided"
			for part in part_list:
				pdb = EMData("particles/"+base_name(part['particle_micrograph_filename'])+"_ptcls.hdf",int(part['index'])-1,True)
				pdbctf = pdb['ctf'].to_dict()
				pdbctf['dfang'] = part['defocusUAngle']
				pdbctf['defocus'] = (float(part['defocusU'])+float(part['defocusV']))/2
				pdbctf['dfdiff'] = abs(float(part['defocusU'])-float(part['defocusV']))
				ctf = EMAN2Ctf()
				ctf.from_dict(pdbctf)
				pdb['ctf'] = ctf
				pdb.write_image("particles/"+base_name(part['particle_micrograph_filename'])+"_ptcls.hdf",int(part['index'])-1)
	elif option1 == "import_box_coordinates":
		current_micrograph = ""
		et = xml.etree.ElementTree.parse(options.import_box_coordinates)
		emx = et.getroot()
		for item in emx:
			if item.tag == "micrograph":
				for micrograph_attrib in item.attrib:
					if micrograph_attrib == "fileName":
						micrograph_filename = item.attrib['fileName']
					elif micrograph_attrib == "index":
						micrograph_index = item.attrib['index']
					else:
						print "Unknown tag: " + micrograph_attrib
				for item2 in item:
					if item2.tag == "acceleratingVoltage":
						voltage = item2.text #in kilovolts
					elif item2.tag == "amplitudeContrast":
						ampcont = item2.text #0->1
					elif item2.tag == "cs":
						cs = item2.text #in mm
					elif item2.tag == "defocusU":
						defocus1 = item2.text # in nm
					elif item2.tag == "defocusV":
						defocus2 = item2.text # in nm
					elif item2.tag == "defocusUAngle":
						defocus_angle = item2.text # in degrees (0->180)
					elif item2.tag == "fom":
						fom = item2.text # figure-of-merit (0->1)
					elif item2.tag == "pixelSpacing":
						for item3 in item2:
							if item3.tag == "X":
								apix_x = item3.text
							elif item3.tag == "Y":
								apix_y = item3.text
							elif item3.tag == "Z":
								apix_z = item3.text
					else:
						print "Unknown tag: " + item2.tag
			elif item.tag == "particle":
				for particle_attrib in item.attrib:
					if particle_attrib == "fileName":
						particle_filename = item.attrib['fileName']
					elif particle_attrib == "index":
						particle = item.attrib['index']
					else:
						print "Unknown tag: " + particle_attrib
				for item2 in item:
					if item2.tag == "boxSize":
						for item3 in item2:
							if item3.tag == "X":
								nx = item3.text #in pixels
							elif item3.tag == "Y":
								ny = item3.text #in pixels
							elif item3.tag == "Z":
								nz = item3.text #in pixels
					elif item2.tag == "centerCoord":
						for item3 in item2:
							if item3.tag == "X":
								center_x = item3.text #in pixels
							elif item3.tag == "Y":
								center_y = item3.text #in pixels
							elif item3.tag == "Z":
								center_z = item3.text #in pixels
					elif item2.tag == "defocusU":
						defocus_particle_1 = item2.text # in nm
					elif item2.tag == "defocusV":
						defocus_particle_2 = item2.text # in nm
					elif item2.tag == "defocusUAngle":
						defocus_particle_angle = item2.text # in degrees (0->180)
					elif item2.tag == "fom":
						particle_fom = item2.text # figure-of-merit (0->1)
					elif item2.tag == "micrograph":
						for micrograph_attrib in item2.attrib:
							if micrograph_attrib == "fileName":
								particle_micrograph_filename = item2.attrib['fileName']
							elif micrograph_attrib == "index":
								particle_micrograph_index = item2.attrib['index']
							else:
								print "Unknown tag: " + micrograph_attrib
					elif item2.tag == "pixelSpacing":
						for item3 in item2:
							if item3.tag == "X":
								apix_x = item3.text
							elif item3.tag == "Y":
								apix_y = item3.text
							elif item3.tag == "Z":
								apix_z = item3.text
					#elif item2.tag == "transformationMatrix":
						#do something...
					else:
						print "Unknown Tag: " + item2.tag
				db = js_open_dict(info_name(particle_micrograph_filename))
				tup = (float(center_x),float(center_y),"manual")
				if 'boxes' not in db:
					db['boxes'] = [tup]
				else:
					db['boxes'].append(tup)
					db.setval("boxes",db['boxes'],deferupdate=True)


	elif option1 == "import_2d_alignment":
		et = xml.etree.ElementTree.parse(options.import_2d_alignment)
		emx = et.getroot()
		for item in emx:
			if item.tag == "particle":
				for particle_attrib in item.attrib:
					if particle_attrib == "fileName":
						particle_filename = item.attrib['fileName']
					elif particle_attrib == "index":
						particle_index = item.attrib['index']
					else:
						print "Unknown tag: " + particle_attrib
			for item2 in item:
				if item2.tag == "transformationMatrix":
					t = Transform([float(item2.find('t11').text),float(item2.find('t12').text),float(item2.find('t13').text),float(item2.find('t14').text),float(item2.find('t21').text),float(item2.find('t22').text),float(item2.find('t23').text),float(item2.find('t24').text),float(item2.find('t31').text),float(item2.find('t32').text),float(item2.find('t33').text),float(item2.find('t34').text)])
					print t















				##if item2.tag == "boxSize":
					##for item3 in item2:
						##if item3.tag == "X":
							##nx = item3.text #in pixels
						##elif item3.tag == "Y":
							##ny = item3.text #in pixels
						##elif item3.tag == "Z":
							##nz = item3.text #in pixels
				##elif item2.tag == "centerCoord":
					##for item3 in item2:
						##if item3.tag == "X":
							##center_x = item3.text #in pixels
						##elif item3.tag == "Y":
							##center_y = item3.text #in pixels
						##elif item3.tag == "Z":
							##center_z = item3.text #in pixels
				##elif item2.tag == "defocusU":
					##defocus_particle_1 = item2.text # in nm
				##elif item2.tag == "defocusV":
					##defocus_particle_2 = item2.text # in nm
				##elif item2.tag == "defocusUAngle":
					##defocus_particle_angle = item2.text # in degrees (0->180)
				##elif item2.tag == "fom":
					##particle_fom = item2.text # figure-of-merit (0->1)
				##elif item2.tag == "micrograph":
					##for micrograph_attrib in item2.attrib:
						##if micrograph_attrib == "fileName":
							##particle_micrograph_filename = item2.attrib['fileName']
						##elif micrograph_attrib == "index":
							##particle_micrograph_index = item2.attrib['index']
						##else:
							##print "Unknown tag: " + micrograph_attrib
				##elif item2.tag == "pixelSpacing":
					##for item3 in item2:
						##if item3.tag == "X":
							##apix_x = item3.text
						##elif item3.tag == "Y":
							##apix_y = item3.text
						##elif item3.tag == "Z":
							##apix_z = item3.text
				#if item2.tag == "transformationMatrix":
					##do something...
				#else:
					#print "Unknown Tag: " + item2.tag
			#db = js_open_dict(info_name(particle_micrograph_filename))
			#tup = (float(center_x),float(center_y),"manual")
			#if 'boxes' not in db:
				#db['boxes'] = [tup]
			#else:
				#db['boxes'].append(tup)
				#db.setval("boxes",db['boxes'],deferupdate=True)
#js_close_dict(info_name(particle_micrograph_filename))
print "e2emx.py finished!"

















