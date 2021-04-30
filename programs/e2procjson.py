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
	parser.add_argument("--retype",type=str, default=None, help="For JSON files where the keys are image name,# pairs, will change the __type value in the image name in all keys")
	parser.add_argument("--extractkey", type=str, default=None, help="This will extract a single named value from each specified file. Output will be multicolumn if the referenced label is an object, such as CTF.")
	parser.add_argument("--extractspt", action="store_true", default=False, help="This will extract the parameters from a particle_parms JSON file in SPT projects as a multicolumn text file.")
	parser.add_argument("--removekey", type=str, default=None, help="DANGER! This will remove all data associated with the named key from all listed .json files.")
	parser.add_argument("--output", type=str, default="jsoninfo.txt", help="Output for text operations (not JSON) filename. default = jsoninfo.txt")
	parser.add_argument("--setoption",type=str, default=None, help="Set a single option in application preferences, eg - display2d.autocontrast:true")
	parser.add_argument("--listoptions",action="store_true", default=False, help="List all currently set user application preferences")
	
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
					out.write("{}\t{}\t{}\t{}\t{}\t{}\t# {}\n".format(v["az"],v["alt"],v["phi"],v["tx"],v["ty"],v["tz"],fsp[:-5]))
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
			


	E2end(logid)


if __name__ == "__main__":
	main()
