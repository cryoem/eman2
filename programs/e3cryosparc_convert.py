#!/usr/bin/env python
#
# Author: Steve Ludtke 04/16/14 (sludtke@bcm.edu)
# Author: Anya Porter 05/21/24 (anastasia.porter@bcm.edu)
# Copyright (c) 2014- Baylor College of Medicine
#
# This software is issued under a joint BSD/GNU license. You may use the
# source code in this file under either license. However, note that the
# complete EMAN2 and SPARX software packages have some GPL dependencies,
# so you are responsible for compliance with the licenses of these packages
# if you opt to use BSD licensing. The warranty disclaimer below holds
# in either instance.
#
# This complete copyright notice must be included in any revised version of the
# source code. Additional authorship citations may be added, but existing
# author citations must be preserved.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  2111-1307 USA
#
#


from EMAN3 import *
from math import *
import numpy as np
import os
import traceback
import re
import line_profiler

@line_profiler.profile
def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options] <CryoSparc .cs file>

Import a CryoSparc .cs file and accompanying images to an EMAN3 style .lst file in a new project.


"""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	#options associated with e2refine.py
	parser.add_header(name="text1", help='Important instructions', title="Use this to create an EMAN3 project from a CryoSparc .cs File:", row=0, col=0, rowspan=1, colspan=3)
	parser.add_header(name="text2", help='Important instructions', title="* cd <folder with CryoSparc .cs file>", row=1, col=0, rowspan=1, colspan=3)
	parser.add_header(name="text3", help='Important instructions', title="* run e3projectmanager, and use this tool", row=2, col=0, rowspan=1, colspan=3)
	parser.add_header(name="text4", help='Important instructions', title="* exit PM, cd eman3, run PM from new eman3 folder", row=3, col=0, rowspan=1, colspan=3)
	parser.add_pos_argument(name="cs_file",help="Select .cs file", default="", guitype='filebox', browser="EMParticlesEditTable(withmodal=True,multiselect=False)",  row=6, col=0,rowspan=1, colspan=3)
	parser.add_argument("--phaseflip",action="store_true",help="If set, will also generate a set of phase-flipped particles. Phase flipped particles will not permit defocus refinement.", guitype='boolbox', row=8, col=0, rowspan=1, colspan=1)
	parser.add_argument("--apix", default=0, type=float,help="The angstrom per pixel of the input particles, if not found in the file.", guitype='floatbox', row=8, col=1, rowspan=1, colspan=1)
	parser.add_argument("--voltage", default=0, type=float,help="Microscope voltage in kV, if not found in the file", guitype='floatbox', row=9, col=0, rowspan=1, colspan=1)
	parser.add_argument("--cs", default=0, type=float,help="Spherical aberration in mm, if not found in the file", guitype='floatbox', row=9, col=1, rowspan=1, colspan=1)
	parser.add_argument("--ac", default=0, type=float,help="Amplitude contrast as a percentage, eg - 10, not 0.1, if not found in the file", guitype='floatbox', row=9, col=2, rowspan=1, colspan=1)
	parser.add_argument("--particlebits", default=6, type=float,help="Significant bits to retain in HDF files for raw particles, 6 is usually more than sufficient (default 6)", guitype='intbox', row=10, col=0, rowspan=1, colspan=1)
	parser.add_argument("--clip", type=int, help="Adjust box size to specified --clip value AFTER CTF phase flipping (if requested). Note that this makes phase un-flipping impossible",default=-1)
	parser.add_argument("--dftolerance",default=0.001, type=float,help="If defocus has to be used to group the particles by micrograph, and the defocus varies per particle, this is the amount of variation in microns to permit within a single file")
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higher number means higher level of verboseness")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)

	(options, args) = parser.parse_args()


	logid=E3init(sys.argv,options.ppid)

	try: os.mkdir("eman3profile2")
	except: pass
	try: os.mkdir("eman3profile2/particles")
	except: pass
	try: os.mkdir("eman3profile2/sets")
	except: pass
	os.chdir("eman3profile2")	# many things need to happen with the project directory as a base

	if options.verbose>0 : print("Parsing .cs file")
	print("filename: ", args[0])
	csfile = np.load("../"+args[0])

	voltage=options.voltage
	cs=options.cs
	apix=options.apix
	ac=options.ac

	try:
		if np.min(csfile["ctf/accel_kv"]) == np.max(csfile["ctf/accel_kv"]):
			voltage=float(csfile["ctf/accel_kv"][0]) # Check to see if it is self consistent for voltage/Cs? (min == max is probably least expensive)
		else:
			print("ctf/accel_kv differs between particles and/or micrographs")
	except: print("Missing ctf/accel_kv, specify --voltage")
	try:
		if np.min(csfile["ctf/cs_mm"]) == np.max(csfile["ctf/cs_mm"]):
			cs=float(csfile["ctf/cs_mm"][0])
		else:
			print("ctf/cs_mm differs between particles and/or micrographs")
	except: print("Missing ctf/cs_mm, specify --cs")
	try: ac=float(csfile["ctf/amp_contrast"][0])
	except: print("Missing ctf/amp_contrast, specify --ac") # Find relion/cryosparc paper to see how they handle phase plate and ctf/amp_contrast and ctf/phase_shift_rad for relion ctffind3 determines ctf for relion
	try: apix=float(csfile["blob/psize_A"][0])
	except: print("Missing blob/psize_A, specify --apix")

	nptcl=len(csfile["ctf/df1_A"])
	ctferrprt=True
	# Need to add something that inverts if blob/sign is not -1
	# Create "ugnums" as a per-particle micrograph identifier
	# Try using the micrograph id first
	try:
		ugnames=set(csfile["location/micrograph_uid"])
		if options.verbose>1: print(f"Identified {len(ugnames)} micrograph uids")
		if len(ugnames)>1 and nptcl/len(ugnames)>3:
			nmtonum={}
			for i,n in enumerate(ugnames): nmtonum[n]=i # Map names to numbers
			ugnums=[nmtonum[n] for n in csfile["location/micrograph_uid"]]
			if options.verbose>0: print("Using location/micrograph_uid to group particles")
	except:
		# If that fails, try using blob/path
		try:
			ugnames=set(csfile["blob/path"])
			# If blob/path seems to be valid and document multiple microgrpahs, we use it
			if len(ugnames)>1 and nptcl/len(ugnames)>3:
				nmtonum={}
				for i,n in enumerate(ugnames): nmtonum[n]=i
				ugnums=[nmtonum[n] for n in csfile["blob/path"]]
				if options.verbose>0: print("Using blob/path to group particles")
			elif options.verbose>0: print("Unable to group particles using blob/path")
		except:
			# final try. If the same defocus was assigned to all images in a micrograph, we can use that
			dfrng=[int(0.0002*df/options.dftolerance) for df in csfile["ctf/df1_A"]] # defocus converted to integers covering the acceptable range of values, 0.002 microns&range
			dfs=set(dfrng)
			if len(dfs)<nptcl/5:
				if options.verbose>0: print(f"Using ctf/df1_A to group particles into {len(dfs)} groups")
				n=0
				ugnums=[0]
				df=dfrng
				for i in range(1,nptcl):
					if df[i]!=df[i-1]: n+=1
					ugnums.append(n)
			if len(dfs)>=nptcl/5 or ugnums[-1]>nptcl/5:
				print(df[:1000])
				print(f'WARNING: Unable to group particles usefully by micrograph ({len(dfs)} defoci, {len(ugnums)} migrographs), collapsing to a single file. Consider rerunning with sufficiently large --dftolerance')
				ugnums=np.zeros(nptcl,"int32")

	# copy particles by micrograph
	# Currently doing this one same as e3relion_convert.py. Optimize once it's working?
	ugnum=-1
	try: os.unlink("sets/fromcs.lst")
	except: pass
	lst=LSXFile("sets/fromcs.lst")
	if options.phaseflip:
		try: os.unlink("sets/fromcs__ctf_flip.lst")
		except: pass
		lstf=LSXFile("sets/fromcs__ctf_flip.lst")
		try: os.unlink("sets/fromcs__ctf_flip_ds5.lst")
		except: pass
		lstft=LSXFile("sets/fromcs__ctf_flip_ds5.lst")
		nbx=good_size(int(options.clip*apix/1.8))
	lastctf=[]
	t0=time.time()
#	for i in range(nptcl):
	for i in range(50000):
		if options.verbose>0 and time.time()-t0>1.0:
			print(f"  {i+1}/{nptcl}\r",end="")
			t0=time.time()
		# Convert CryoSparc info to EMAN2 conventions
		try:
			df1=float(csfile["ctf/df1_A"][i]) # df1_A directly matches rlnDefocusU
			df2=float(csfile["ctf/df2_A"][i]) # df2_A directly matches rlnDefocusV
			dfang=float(csfile["ctf/df_angle_rad"][i])
			defocus=(df1+df2)/20000.0
			dfdiff=(df1-df2)/10000.0
		except:
			if ctferrprt: print("ERROR: could not determine defocus from cs file")
			ctferrprt=False
			defocus,dfiff,dfang=1.0, 0.0, 0.0  # 1 micron default if we can't find a good defocus
		try: bfactor=float(csfile["ctf/bfactor"][i])
		except: bfactor=50.0
		try:
			pose = csfile["alignments3D/pose"][i]
			twist = np.linalg.norm(pose)
			pose = pose/twist
			twist = twist*180/pi
			shift = csfile["alignments3D/shift"][i]
			xform=Transform({"type":"spin","n1":float(pose[0]),"n2":float(pose[1]),"n3":float(pose[2]),"omega":float(twist),"tx":float(-shift[0]),"ty":float(-shift[1])})
		except Exception:
			if i==0:
				traceback.print_exc()
				print("No usable particle orientations found in .cs file!")
			xform=Transform()

		# expensive to compute, but since CryoSparc has per-particle CTF, we preserve this even when "grouping by micrograph"
		ctf=EMAN2Ctf()
		ctf.from_dict({"defocus":defocus,"dfang":dfang,"dfdiff":dfdiff,"voltage":voltage,"cs":cs,"ampcont":ac,"apix":apix,"bfactor":bfactor})

		# begin a new micrograph file
		if ugnums[i]!=ugnum:
			ugnum=ugnums[i] # Add micrograph names to info so could relate (need to update wiki too)
			ptclname=f"particles/micro_{ugnum:05d}.hdf:{int(options.particlebits)}"
			ptclfname=f"particles/micro_{ugnum:05d}__ctf_flip.hdf:{int(options.particlebits)}"
			ptclfsname=f"particles/micro_{ugnum:05d}__ctf_flip_ds5.hdf:{int(options.particlebits)}"
			ptclno=0
			jdb=js_open_dict(info_name(ptclname))

			# Make a "micrograph" CTF entry for each set of different defocuses to use when fitting
			jdb["ctf_frame"]=[1024,ctf,(256,256),tuple(),5,1]

		imfsp=csfile["blob/path"][i].decode('UTF-8')[1:] # blob/path is numpy.bytes and has a leading > to be removed
		imfsp=re.sub(r'/\d+_', '/', imfsp) #They also have the micrograph uid added to the filepath at the start of the filename, after the relative path
		imn=csfile["blob/idx"][i]

		img=EMData("../"+imfsp,imn) # CryoSparc and Eman start with 0
		img["apix_x"]=apix
		img["apix_y"]=apix
		img["apix_z"]=apix
		img["ctf"]=ctf                                          # CTF with good general parameters in the header, can be overridden in lst files
		img["ctf_phase_flipped"]=0

		# write the particle to an image stack and add an entry to the .lst file (-1 appends)
		if options.clip>0:
			cp=img.get_clip(Region((img["nx"]-options.clip)//2,(img["ny"]-options.clip)//2,options.clip,options.clip))
			cp.write_image(ptclname,ptclno)
		else: img.write_image(ptclname,ptclno)
		lst[-1]=(ptclno,ptclname.split(":")[0],{"xform.projection":xform})

		# write the phase flipped image/lst if requested
		if options.phaseflip:
			imgf=img.do_fft()
			if  lastctf!=ctf.to_vector():
				flipim=imgf.copy()
				ctf.compute_2d_complex(flipim,Ctf.CtfType.CTF_SIGN)
			imgf.mult(flipim)
			imgc=imgf.do_ift()
			imgc["apix_x"]=apix
			imgc["apix_y"]=apix
			imgc["apix_z"]=apix
			imgc["ctf"]=ctf
			imgc["ctf_phase_flipped"]=1
			if options.clip>0:
				cp=imgc.get_clip(Region((img["nx"]-options.clip)//2,(img["ny"]-options.clip)//2,options.clip,options.clip))
				cp.write_image(ptclfname,ptclno)
				shr=cp.process("math.fft.resample",{"n":options.clip/nbx})
				shr.write_image(ptclfsname,ptclno)
			else:
				imgc.write_image(ptclfname,ptclno)
				shr=imgc.process("math.fft.resample",{"n":imgc["nx"]/nbx})
				shr.write_image(ptclfsname,ptclno)
			lstf[-1]=(ptclno,ptclfname.split(":")[0],{"xform.projection":xform})
			lstft[-1]=(ptclno,ptclfsname.split(":")[0],{"xform.projection":xform})

			lastctf=ctf.to_vector()
		ptclno+=1
	E3end(logid)




if __name__ == "__main__":
    main()
