#!/usr/bin/env python
# This program performs simple processing of .json files

# Author: Steven Ludtke, 07/26/2017 (sludtke@bcm.edu), modified: May 15, 2017 (Jesus GalazMontoya)
# Copyright (c) 2017- Baylor College of Medicine
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

from EMAN2 import *
from math import *
import os
import sys

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """Usage:\nprocjson.py [options] <json 1> <json 2> ... 
	
	Provides simple utility functions for examining and modifying JSON files. Note that JSON (JavaScript Object Notation) 
	is a text format, so files can also be read/processed with standard text editing and processing tools."""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	####################
	
	parser.add_argument("--allinfo", action="store_true", default=False, help="Uses all of the .json files in info/ rather than specifying a list on the command line")
	parser.add_argument("--listkeys",action="store_true", default=False, help="Lists all of the keys in all of the specified info files")
	parser.add_argument("--remaplstkeys",action="store_true", default=False, help="For JSON files where the keys are image name,# pairs referencing a .lst file, will replace each key with the original image")
	parser.add_argument("--dump",action="store_true", default=False, help="Nicely print the entire contents of the JSON file to the screen. More readable than simply looking at the file.")
	parser.add_argument("--retype",type=str, default=None, help="For JSON files where the keys are image name,# pairs, will change the __type value in the image name in all keys")
	parser.add_argument("--extractkey", type=str, default=None, help="This will extract a single named value from each specified file. Output will be multicolumn if the referenced label is an object, such as CTF.")
	parser.add_argument("--extractspt", action="store_true", default=False, help="This will extract the parameters from a particle_parms JSON file in SPT projects as a multicolumn text file.")
	parser.add_argument("--extractboxes", action="store_true", default=False, help="This will save 2-D particle box locations to one or more text files")
	parser.add_argument("--extractboxes3d", action="store_true", default=False, help="This will save 3-D particle box locations to one or more text files")
	parser.add_argument("--tilebox3d", type=str, default=None, help="Generate a tiling of 3-D box locations for a tomogram on the z=0 plane: <x0>,<y0>,<dxy>[,name]")
	parser.add_argument("--autoloadboxes", action="store_true", default=False, help="Provide a list of .txt files using the same naming convention produced by extractboxes[3d]. Boxes from files will be added to any existing boxes in corresponding info/*.json files")
	parser.add_argument("--removekey", type=str, default=None, help="DANGER! This will remove all data associated with the named key from all listed .json files.")
	parser.add_argument("--addkey", type=str, default=None, help="add a simple key in the format of key:value. Will try to conver value to float if possible.")
	parser.add_argument("--removeptcl", type=str, default=None, help="Remove tomo particle box locations for the specified label.")

	parser.add_argument("--output", type=str, default="jsoninfo.txt", help="Output for text operations (not JSON) filename. default = jsoninfo.txt")
	parser.add_argument("--setoption",type=str, default=None, help="Set a single option in application preferences, eg - display2d.autocontrast:true")
	parser.add_argument("--listoptions",action="store_true", default=False, help="List all currently set user application preferences")
	parser.add_argument("--infofile", action="store_true", default=False, help="Find and process the info file corresponding to the input file.")

	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, help="verbose level [0-9], higher number means higher level of verboseness",default=1)


	(options, args) = parser.parse_args()

	if options.setoption!=None:
		try:
			k,v=options.setoption.split(":")
			com,opt=k.split(".")
			if isinstance(v,str) and v.lower()=="true" : v=True
			elif isinstance(v,str) and v.lower()=="false" : v=False
			E2setappval(com,opt,v)
		except:
#			import traceback
#			traceback.print_exc()
			print("ERROR: could not write preferences. Must be of form 'program.option:value'")
			sys.exit(1)
		sys.exit(0)
	
	if options.infofile:
		newargs=[info_name(a) for a in args]
		for a,n in zip(args, newargs):
			print("{} -> {}".format(a, n))
		args=newargs

	if options.tilebox3d is not None:
		for ifsp in args:
			js=js_open_dict(ifsp)
			b3d=js["boxes_3d"]
			name="tiled"
			try:
				x0,y0,dxy=[int(i) for i in options.tilebox3d.split(",")]
			except:
				try: 
					x0,y0,dxy,name=options.tilebox3d.split(",")
					x0,y0,dxy=int(x0),int(y0),int(dxy)
				except:	
					print("ERROR: --tilebox3d <x0>,<y0>,<dxy> in raw tilt-series coordinates")
					sys.exit(1)
			try:
				cls=js["class_list"]
				clsno=max([int(k) for k in cls.keys()])+1
				cls[str(clsno)]={"boxsize":dxy,"name":name}
			except:
				clsno=0
				cls={str(clsno):{"boxsize":dxy,"name":tiled}}
			js["class_list"]=cls

			for y in range(y0,-y0+dxy//2,dxy):
				for x in range(x0,-x0+dxy//2,dxy):
					b3d.append([x,y,0,"manual",0.0,clsno])

			js["boxes_3d"]=b3d
		
	if options.remaplstkeys:
		try: os.remove("tmp_tmp.json")
		except: pass
		for ifsp in args:
			ls=None
			js=js_open_dict(ifsp)
			jsout=js_open_dict("tmp_tmp.json")
			for k in js.keys():
				fsp,n=eval(k)
				if ls!=fsp:
					ls=fsp
					lsx=LSXFile(fsp)
				n2,fsp2,tmp=lsx.read(n)
				jsout.setval(str((fsp2,n2)),js[k],True)
				
			jsout.sync()
			jsout=None
			js=None
			os.remove(ifsp)
			os.rename("tmp_tmp.json",ifsp)

	if options.retype!=None:
		try: os.remove("tmp_tmp.json")
		except: pass
		for ifsp in args:
			ls=None
			js=js_open_dict(ifsp)
			jsout=js_open_dict("tmp_tmp.json")
			for k in js.keys():
				fsp,n=eval(k)
				fsp2=fsp.rsplit("__",1)[0]+"__"+options.retype+".hdf"
				jsout.setval(str((fsp2,n)),js[k],True)
				
			jsout.sync()
			jsout=None
			js=None
			os.remove(ifsp)
			os.rename("tmp_tmp.json",ifsp)
		
	if options.listoptions:
		prefs=E2getappvals()
		if len(prefs)==0 : 
			print("No preferences have been set. Please see 'http://eman2.org/ApplicationPreferences' for a list of available preferences.")
			sys.exit(0)
		for a,s,v,gl in prefs: print("{:>30s} {:25s} {}".format(".".join((a,s)),str(v),gl))
		sys.exit(0)

	if options.allinfo:
		args=["info/{}".format(i) for i in os.listdir("info") if ".json" in i]

	if options.verbose>1: print(len(args)," json files to process")

	if len(args)<1 :
		parser.error("At least one lst file required")
		sys.exit(1)

	logid=E2init(sys.argv,options.ppid)

	if options.listkeys:
		allkeys=set()
		for fsp in args:
			js=js_open_dict(fsp)
			allkeys.update(js.keys())
			js.close()
		
		for k in sorted(list(allkeys)): print(k)
		
	if options.dump:
		from pprint import pprint
		allkeys=set()
		for fsp in args:
			if len(args)>1: print("\n",fsp,":")
			js=js_open_dict(fsp)
			pprint(dict(js.items()))
			js.close()

	if options.autoloadboxes:
		for fsp in args:
			fb=os.path.basename(fsp)
			if fb.startswith("box_"):
				try:
					iname=fb[4:].split("__")[0]				# base name for json file
					tp=fb.split("__")[1].rsplit(".",1)[0]	# box type
				except:
					error_exit("ERROR: --autoloadboxes expects a list of .txt files of the form box[3d]_<micrograph>__<boxtype>.txt, as produced by --extractboxes or --extractboxes3d")
				js=js_open_dict(info_name(iname))
				try: bxs=js["boxes"]
				except: bxs=[]
				for l in open(fsp,"r"):
					if l[0]=="#": continue
					bxs.append([float(l.split()[0]),float(l.split()[1]),tp])
				js["boxes"]=bxs
			elif fb.startswith("box3d_"):
				try:
					iname=fb[6:].split("__")[0]				# base name for json file
					tp=fb.split("__")[1].rsplit(".",1)[0]	# box type
				except:
					error_exit("ERROR: --autoloadboxes expects a list of .txt files of the form box[3d]_<micrograph>__<boxtype>.txt, as produced by --extractboxes or --extractboxes3d")
				js=js_open_dict(info_name(iname))
				try: bxs=js["boxes_3d"]
				except: bxs=[]
				for l in open(fsp,"r"):
					if l[0]=="#": continue
					ls=l.split()
					bx=[float(ls[0]),float(ls[1]),float(ls[2]),tp]
					if len(ls)>3: bx.append(float(ls[3]))
					if len(ls)>4: bx.append(int(ls[4]))
					bxs.append(bx)
				js["boxes_3d"]=bxs

			else:
				error_exit("ERROR: --autoloadboxes expects a list of .txt files of the form box[3d]_<micrograph>__<boxtype>.txt, as produced by --extractboxes or --extractboxes3d")

	if options.extractboxes:
		allkeys=set()
		for fsp in args:
			js=js_open_dict(fsp)
			if "boxes" in js:
				bxs=js["boxes"]			# array of all boxes for this file
				tps={b[2]:open(f"box_{os.path.basename(fsp)[:-5]}__{b[2].replace('.','-')}.txt","w") for b in bxs}	# dict of open files for each box type
				for b in bxs:
					tps[b[2]].write(f"{b[0]:1.1f}\t{b[1]:1.1f}\n")
			js.close()
			tps=None		# implicitly close the files

	if options.extractboxes3d:
		allkeys=set()
		for fsp in args:
			js=js_open_dict(fsp)
			if "boxes_3d" in js:
				bxs=js["boxes_3d"]			# array of all boxes for this file
				tps={b[3]:open(f"box3d_{os.path.basename(fsp)[:-5]}__{b[3].replace('.','-')}.txt","w") for b in bxs}	# dict of open files for each box type
				for b in bxs:
					tps[b[3]].write(f"{b[0]:1.1f}\t{b[1]:1.1f}\t{b[2]:1.1f}")
					if len(b)>5 : tps[b[3]].write(f"\t{b[4]:1.1f}\t{b[5]:d}\n")
					elif len(b)>4 : tps[b[3]].write(f"\t{b[4]:1.1f}\n")
					else : tps[b[3]].write("\n")
			js.close()
			tps=None		# implicitly close the files


	if options.extractspt :
		out=open(options.output,"w")
		out.write("# file_ptcl,score,trans_x,trans_y,trans_z,az,alt,phi,rel_trans_x,rel_trans_y,rel_trans_z,rel_alt,rel_az,rel_phi\n")
		nf=0
		for fsp in args:
			js=js_open_dict(fsp)
			ks=[(eval(k)[1],eval(k)[0],k) for k in js.keys()]
			for k in sorted(ks):
				xf=js[k[2]]["xform.align3d"]
				if "xform.start" in js[k[2]] : xfd=xf*js[k[2]]["xform.start"][0].inverse()
				else: xfd=None
				tr=xf.get_trans()
				rt=xf.get_rotation()
				out.write(f"{k[0]}\t{js[k[2]]['score']}\t{tr[0]}\t{tr[1]}\t{tr[2]}\t{rt['az']}\t{rt['alt']}\t{rt['phi']}")
				if xfd!=None:
					tr=xfd.get_trans()
					rt=xfd.get_rotation()
					out.write(f"\t{tr[0]}\t{tr[1]}\t{tr[2]}\t{rt['az']}\t{rt['alt']}\t{rt['phi']}")
				out.write(f"\t# {k[2]}\n")

	if options.extractkey :
		out=open(options.output,"w")
		nf=0
		for fsp in args:
			js=js_open_dict(fsp)
			if options.extractkey in js:
				v=js[options.extractkey]
				nf+=1
				if isinstance(v,list) and len(v)>0 and isinstance(v[0],EMAN2Ctf): v=v[0]
				
				if isinstance(v,EMAN2Ctf) :
					out.write("{:1.5f}\t{:1.1f}\t{:1.4f}\t{:1.5f}\t{:1.2f}\t{:1.1f}\t{:1.3f}\t{:1.2f}\t# {}\n".format(v.defocus,v.bfactor,v.apix,v.dfdiff,v.dfang,v.voltage,v.cs,v.get_phase(),fsp[:-5]))
				elif isinstance(v,EMData) :
					out.write("{}\t{}\t{}\t{}\t{}\t# {}\n".format(v["nx"],v["ny"],v["nz"],v["mean"],v["sigma"],fsp[:-5]))
				elif isinstance(v,Transform) :
					dct=v.get_params("eman")
					out.write("{}\t{}\t{}\t{}\t{}\t{}\t# {}\n".format(dct["az"],dct["alt"],dct["phi"],dct["tx"],dct["ty"],dct["tz"],fsp[:-5]))
				else:
					out.write("{}\t# {}\n".format(v,fsp[:-5]))
			
			js.close()
		if options.verbose: print("{} found in {} JSON files".format(options.extractkey,nf))
					
	if options.removekey:
		jsb={}
		
		nf=0
		for ii, fsp in enumerate(args):
			
			if options.verbose>5:
				sys.stdout.write("\r{}/{} finished.".format(ii, len(args)))
				sys.stdout.flush()
			js=js_open_dict(fsp)
			if options.removekey in js:
				nf+=1
				v=js[options.removekey]
				jsb[fsp]=v
				del js[options.removekey]
			js.close()

		js=js_open_dict("backup_removed.json")
		js.update(jsb)
		js.close()
		print("Removed {} from {} files. Backup stored in backup_removed.json".format(options.removekey,nf))
			
	if options.addkey:
		key=options.addkey.split(':')[0]
		val=options.addkey.split(':')[1]
		try: val=float(val)
		except: pass
	
		print(key, val)
		for ii, fsp in enumerate(args):
			js=js_open_dict(fsp)
			js[key]=val
			js.close()
			
			
	if options.removeptcl:
		
		for ii, fsp in enumerate(args):
				
			js=dict(js_open_dict(fsp)).copy()
			if "class_list" in js:
				ky=list(js["class_list"].keys())
				torm=[k for k in ky if js["class_list"][k]["name"]==options.removeptcl]
				if len(torm)==0:
					continue
				boxnew=js["boxes_3d"]
				nb=len(boxnew)
				for ir in torm:
					js["class_list"].pop(ir)
					boxnew=[b for b in boxnew if b[-1]!=int(ir)]
					
				print("{} -> remove {} classes, {} particles".format(base_name(fsp), len(torm), nb-len(boxnew)))
				js["boxes_3d"]=boxnew
				
				f=js_open_dict(fsp)
				f.update(js)
				f.close()
				
			else:
				continue
			

	E2end(logid)


if __name__ == "__main__":
	main()
