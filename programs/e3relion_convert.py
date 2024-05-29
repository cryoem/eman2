#!/usr/bin/env python
#
# Author: Steve Ludtke 04/16/14 (sludtke@bcm.edu)
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
from EMAN3star import StarFile
from math import *
import numpy as np
import os
import traceback


def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options] <relion STAR file>

Import a Relion star file and accompanying images to an EMAN3 style .lst file in a new project.


"""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	#options associated with e2refine.py
	parser.add_header(name="text1", help='Important instructions', title="Use this to create an EMAN3 project from a Relion Star File:", row=0, col=0, rowspan=1, colspan=3)
	parser.add_header(name="text2", help='Important instructions', title="* cd <folder with Relion STAR file>", row=1, col=0, rowspan=1, colspan=3)
	parser.add_header(name="text3", help='Important instructions', title="* run e3projectmanager, and use this tool", row=2, col=0, rowspan=1, colspan=3)
	parser.add_header(name="text4", help='Important instructions', title="* exit PM, cd eman3, run PM from new eman2 folder", row=3, col=0, rowspan=1, colspan=3)
	parser.add_pos_argument(name="star_file",help="Select STAR file", default="", guitype='filebox', browser="EMParticlesEditTable(withmodal=True,multiselect=False)",  row=6, col=0,rowspan=1, colspan=3)
	parser.add_argument("--phaseflip", action="store_true",help="If set, will also generate a set of phase-flipped particles. Phase flipped particles will not permit defocus refinement.", guitype='floatbox', row=8, col=0, rowspan=1, colspan=1)
	parser.add_argument("--apix", default=0, type=float,help="The angstrom per pixel of the input particles, if not found in the file.", guitype='floatbox', row=8, col=0, rowspan=1, colspan=1)
	parser.add_argument("--voltage", default=0, type=float,help="Microscope voltage in kV, if not found in STAR file", guitype='floatbox', row=8, col=0, rowspan=1, colspan=1)
	parser.add_argument("--cs", default=0, type=float,help="Spherical aberration in mm, if not found in STAR file", guitype='floatbox', row=8, col=0, rowspan=1, colspan=1)
	parser.add_argument("--ac", default=0, type=float,help="Amplitude contrast as a percentage, eg - 10, not 0.1, if not found in STAR file", guitype='floatbox', row=8, col=0, rowspan=1, colspan=1)
	parser.add_argument("--particlebits", default=6, type=float,help="Significant bits to retain in HDF files for raw particles, 6 is usually more than sufficient (default 6)", guitype='intbox', row=9, col=0, rowspan=1, colspan=1)
	parser.add_argument("--clip", type=int, help="Adjust box size to specified --clip value AFTER CTF phase flipping (if requested). Note that this makes phase un-flipping impossible",default=-1)
	parser.add_argument("--dftolerance",default=0.001, type=float,help="If defocus has to be used to group the particles by micrograph, and the defocus varies per particle, this is the amount of variation in microns to permit within a single file")
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higher number means higher level of verboseness")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)

	(options, args) = parser.parse_args()


	logid=E3init(sys.argv,options.ppid)

	try: os.mkdir("eman3")
	except: pass
	try: os.mkdir("eman3/particles")
	except: pass
	try: os.mkdir("eman3/sets")
	except: pass
	os.chdir("eman3")	# many things need to happen with the project directory as a base

	if options.verbose>0 : print("Parsing STAR file")
	star=StarFile("../"+args[0])

	voltage=options.voltage
	cs=options.cs
	apix=options.apix
	ac=options.ac

	# newer relion file
	if "optics" in star:
		if len(star["optics"]["rlnVoltage"])>1 :
			print("WARNING: Multiple optics groups detected in star file. EMAN3 does not support this at present. The first optics group will be used for all data, which may be incorrect.")
		try: voltage=float(star["optics"]["rlnVoltage"][0])
		except: print("data_optics missing rlnVoltage")
		try: cs=float(star["optics"]["rlnSphericalAberration"][0])
		except: print("data_optics missing rlnSphericalAberration")
		try: ac=float(star["optics"]["rlnAmplitudeContrast"][0])
		except: print("data_optics missing rlnAmplitudeContrast")
		try: apix=float(star["optics"]["rlnImagePixelSize"][0])
		except: print("data_optics missing rlnImagePixelSize")
		rkey="particles"	# data_particles in new files
		if "rlnOriginXAngst" in star[rkey]:			# ran into a file like this, just in case, catch it
			star[rkey]["rlnOriginXAngstrom"]=star[rkey]["rlnOriginXAngst"]
			star[rkey]["rlnOriginYAngstrom"]=star[rkey]["rlnOriginYAngst"]
	#older relion file
	else:
		if len(star.keys())>1 : print(f"WARNING: Star file has multiple data_ blocks ({star.keys()}), but no data_optics block. Not sure how to support this.")
		rkey=list(star.keys())[0]			# generally just one block in old files usually "data_"
		try: voltage=float(star[rkey]["rlnVoltage"][0])
		except: print("missing rlnVoltage, specify --voltage")
		try: cs=float(star[rkey]["rlnSphericalAberration"][0])
		except: print("missing rlnSphericalAberration, specify --cs")
		try: ac=float(star[rkey]["rlnAmplitudeContrast"][0])
		except: print("missing rlnAmplitudeContrast, specify --ac")
		try: apix=float(star[rkey]["rlnDetectorPixelSize"][0])*10000.0/float(star[rkey]["rlnMagnification"][0])
		except: print("missing rlnDetectorPixelSize and/or rlnMagnification, specify --apix")

	# at this point rkey should contain the star key for the particle metadata

	nptcl=len(star[rkey]["rlnAnglePsi"])
	ctferrprt=True

	### Here we create "ugnums" a per-particle micrograph identifier
	# first we try using the micrograph name field
	try:
		ugnames=set(star[rkey]["rlnMicrographName"])
		if options.verbose>1: print(f"Identified {len(ugnames)} rlnMicrographNames")

		# if the rlnMicrographName seems to be valid, and document multiple micrographs, we use it
		if len(ugnames)>1 and nptcl/len(ugnames)>3:
			nmtonum={}
			for i,n in enumerate(ugnames): nmtonum[n]=i		# map micrograph names to numbers
			ugnums=[nmtonum[n] for n in star[rkey]["rlnMicrographName"]]
			if options.verbose>0: print("Using rlnMicrographName to group particles")
	except:
		# if that fails, we see if imagename will be useful. Since that relies on separate mrcs stacks, in most cases it will also fail
		try:
			ugnames=set([i.split("@")[1] for i in star[rkey]["rlnImageName"]])

			# if the rlnMicrographName seems to be valid, and document multiple micrographs, we use it
			if len(ugnames)>1 and nptcl/len(ugnames)>3:
				nmtonum={}
				for i,n in enumerate(ugnames): nmtonum[n]=i		# map micrograph names to numbers
				ugnums=[nmtonum[n] for n in star[rkey]["rlnImageName"].split("@")[1]]
				if options.verbose>0: print("Using rlnImageName to group particles")
			if options.verbose>0: print("Unable to group particles using rlnImageName")
		except:
			# final try. If the same defocus was assigned to all images in a micrograph, we can use that
			dfrng=[int(.0002*df/options.dftolerance) for df in star[rkey]["rlnDefocusU"]]  # defocus converted to integers covering the acceptable range of values, .002 microns&range
			dfs=set(dfrng)
			if len(dfs)<nptcl/5:
				if options.verbose>0: print(f"Using rlnDefocusU to group particles into {len(dfs)} groups")
				n=0
				ugnums=[0]
				df=dfrng
				for i in range(1,nptcl):
					if df[i]!=df[i-1]: n+=1
					ugnums.append(n)
			if len(dfs)>=nptcl/5 or ugnums[-1]>nptcl/5:
				print(df[:1000])
				print(f'WARNING: Unable to group particles usefully by micrograph ({len(dfs)} defoci, {len(ugnums)} "micrographs"), collapsing to a single file. Consider rerunning with sufficiently large --dftolerance')
				ugnums=np.zeros(nptcl,"int32")

	# copy particles by micrograph
	# doing this one at a time in the first version. May go back and optimize once it's working
	ugnum=-1
	try: os.unlink("sets/fromstar.lst")
	except: pass
	lst=LSXFile("sets/fromstar.lst")
	if options.phaseflip: 
		try: os.unlink("sets/fromstar__ctf_flip.lst")
		except: pass
		lstf=LSXFile("sets/fromstar__ctf_flip.lst")
		try: os.unlink("sets/fromstar__ctf_flip_ds5.lst")
		except: pass
		lstft=LSXFile("sets/fromstar__ctf_flip_ds5.lst")
	lastctf=[]
	t0=time.time()
	for i in range(nptcl):
		if options.verbose>0 and time.time()-t0>1.0:
			print(f"  {i+1}/{nptcl}\r",end="")
			t0=time.time()

		# Convert Relion info to EMAN2 conventions
		try:
			dfu=star[rkey]["rlnDefocusU"][i]
			dfv=star[rkey]["rlnDefocusV"][i]
			dfang=star[rkey]["rlnDefocusAngle"][i]
			defocus=(dfu+dfv)/20000.0
			dfdiff=(dfu-dfv)/10000.0
		except:
			if ctferrprt: print("ERROR: could not determine defocus from STAR file")
			ctferrprt=False
			defocus,dfdiff,dfang=1.0,0.0,0.0		# 1 micron default if we can't find a good defocus
		try: bfactor=star[rkey]["rlnCtfBfactor"][i]
		except: bfactor=50.0
		try: xform=Transform({"type":"spider","phi":star[rkey]["rlnAngleRot"][i],"theta":star[rkey]["rlnAngleTilt"][i],"psi":star[rkey]["rlnAnglePsi"][i],"tx":-star[rkey]["rlnOriginX"][i],"ty":-star[rkey]["rlnOriginY"][i]})
		except:
			try: xform=Transform({"type":"spider","phi":star[rkey]["rlnAngleRot"][i],"theta":star[rkey]["rlnAngleTilt"][i],"psi":star[rkey]["rlnAnglePsi"][i],"tx":-star[rkey]["rlnOriginXAngstrom"][i]/apix,"ty":-star[rkey]["rlnOriginYAngstrom"][i]/apix})
			except:
				traceback.print_exc()
				if i==0: print("No usable particle orientations found in STAR file!")
				xform=Transform()

		# expensive to compute, but since so many relion files have per-particle CTF, we preserve this even when "grouping by micrograph"
		ctf=EMAN2Ctf()
		ctf.from_dict({"defocus":defocus,"dfang":dfang,"dfdiff":dfdiff,"voltage":voltage,"cs":cs,"ampcont":ac,"apix":apix,"bfactor":bfactor})

		# begin a new micrograph file
		if ugnums[i]!=ugnum:
			ugnum=ugnums[i]
			ptclname=f"particles/micro_{ugnum:05d}.hdf:{options.particlebits}"
			ptclfname=f"particles/micro_{ugnum:05d}__ctf_flip.hdf:{options.particlebits}"
			ptclfsname=f"particles/micro_{ugnum:05d}__ctf_flip_ds5.hdf:{options.particlebits}"
			ptclno=0
			jdb=js_open_dict(info_name(ptclname))

			# Make a "micrograph" CTF entry for each set of different defocuses to use when fitting
			jdb["ctf_frame"]=[1024,ctf,(256,256),tuple(),5,1]

		imn,imfsp=star[rkey]["rlnImageName"][i].split("@")
		imn=int(imn)

		img=EMData("../"+imfsp,imn-1)		# relion starts with #1, EMAN #0
		img["apix_x"]=apix
		img["apix_y"]=apix
		img["apix_z"]=apix
		img["ctf"]=ctf						# CTF with good general parameters in the header, can be overridden in lst files
		img["ctf_phase_flipped"]=0

		# write the particle to an image stack, and add an entry to the .lst file (-1 appends)
		if options.clip>0:
			cp=img.get_clip(Region((img["nx"]-options.clip)//2,(img["ny"]-options.clip)//2,options.clip,options.clip))
			cp.write_image(ptclname,ptclno)
		else: img.write_image(ptclname,ptclno)
		lst[-1]=(ptclno,ptclname.split(":")[0],{"xform.projection":xform})

		# write the phase flipped image/lst if requested
		if options.phaseflip:
			imgf=img.do_fft()
			if lastctf!=ctf.to_vector():
				flipim=imgf.copy()
				ctf.compute_2d_complex(flipim,Ctf.CtfType.CTF_SIGN)
			imgf.mult(flipim)
			imgc=imgf.do_ift()
			imgc["apix_x"]=apix
			imgc["apix_y"]=apix
			imgc["apix_z"]=apix
			imgc["ctf"]=ctf						# CTF with good general parameters in the header, can be overridden in lst files
			imgc["ctf_phase_flipped"]=1
			if options.clip>0:
				cp=imgc.get_clip(Region((img["nx"]-options.clip)//2,(img["ny"]-options.clip)//2,options.clip,options.clip))
				cp.write_image(ptclfname,ptclno)
				nbx=good_size(int(options.clip*apix/1.8))
				shr=cp.process("math.fft.resample",{"n":options.clip/nbx})
				shr.write_image(ptclfsname,ptclno)
			else:
				imgc.write_image(ptclfname,ptclno)
				nbx=good_size(int(imgc["nx"]*apix/1.8))
				shr=imgc.process("math.fft.resample",{"n":imgc["nx"]/nbx})
				shr.write_image(ptclfsname,ptclno)
			lstf[-1]=(ptclno,ptclfname.split(":")[0],{"xform.projection":xform})
			lstft[-1]=(ptclno,ptclfsname.split(":")[0],{"xform.projection":xform})

		lastctf=ctf.to_vector()

		ptclno+=1


	E3end(logid)

if __name__ == "__main__":
    main()
