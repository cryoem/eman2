#!/usr/bin/env python
#
# Author: Steven Ludtke, 04/10/2003 (sludtke@bcm.edu), Last update: 11/27/201 (Jesus Galaz)
# Copyright (c) 2000-2006 Baylor College of Medicine
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

from builtins import range
import pprint
from EMAN2 import *
import sys

def get_data_type_string(datatype):
	dtstring = {
		0 : "UNKNOWN",
		1 : "CHAR",
		2 : "UNSIGNED CHAR",
		3 : "SHORT",
		4 : "UNSIGNED SHORT",
		5 : "INT",
		6 : "UNSIGNED INT",
		7 : "FLOAT",
		8 : "DOUBLE",
		9 : "SHORT_COMPLEX",
		10 : "USHORT_COMPLEX",
		11 : "FLOAT_COMPLEX"
	}
	return dtstring.get(datatype, 'UNKNOWN')

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options] imagefile ...

	This program can be used to extract various metadata/header information from images of any file format,
	including BDB databases (though e2bdb.py has additional database-specific functionality).
	"""
	
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	
	parser.add_argument("-H", "--header", action="store_true",help="Show all header information",default=False)
	parser.add_argument("-N", "--number", type=int, help="Image number for single image info",default=-1)
	parser.add_argument("-Q", "--quality", type=int, help="Include only images with a single quality value (integer 0-9)",default=-1)
	parser.add_argument("--dfmin", type=float, help="Include only images with defocus >= the specified value",default=None)
	parser.add_argument("--dfmax", type=float, help="Include only images with defocus <= the specified value",default=None)
	parser.add_argument("--nameonly",action="store_true",help="Only display the matching filenames. No other info.",default=False)
	parser.add_argument("-s", "--stat", action="store_true",help="Show statistical information about the image(s).",default=False)
	parser.add_argument("-O", "--outliers", action="store_true",help="Statistics on number of outlier values, must be used with --stat",default=False)
	parser.add_argument("-E", "--euler", action="store_true",help="Show Euler angles from header",default=False)
	parser.add_argument("-a", "--all", action="store_true",help="Show info for all images in file",default=False)
	parser.add_argument("-C", "--check", action="store_true",help="Checks to make sure all image numbers are populated with images, and that all images have valid CTF parameters",default=False)
	parser.add_argument("-c", "--count", action="store_true",help="Just show a count of the number of particles in each file",default=False)
	parser.add_argument("--extractkey", action="append",dest="extractkey", default=None, help="This will extract a single named header value from each image as a text file. Output will be multicolumn if the referenced label is an object, such as CTF.")

	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-2)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higher number means higher level of verboseness")
	
	(options, args) = parser.parse_args()

	if len(args)<1:
		print(usage)
		parser.error("Specify image file")

	if options.header and options.stat : 
		print("Only one of --header and --stat may be specified")
		sys.exit(1)


	if options.extractkey is not None:
		nf=0
		for fsp in args:
			out=open(f"{fsp.rsplit('.',1)[0]}_{options.extractkey[0]}.txt","w")
			if fsp.endswith(".lst"):
				lsx=LSXFile(fsp)
				hdrs=[e[2] for e in lsx]
			else: hdrs=[h.get_attr_dict() for h in EMData.read_images(fsp,None,True)]

			for i,h in enumerate(hdrs):
				for k in options.extractkey:
					if k in hdrs[0].keys():
						v=h[k]

						if isinstance(v,EMAN2Ctf) :
							out.write("{:1.5f}\t{:1.1f}\t{:1.4f}\t{:1.5f}\t{:1.2f}\t{:1.1f}\t{:1.3f}\t{:1.2f}".format(v.defocus,v.bfactor,v.apix,v.dfdiff,v.dfang,v.voltage,v.cs,v.get_phase()))
						elif isinstance(v,Transform) :
							dct=v.get_params("eman")
							out.write("{}\t{}\t{}\t{}\t{}\t{}".format(dct["az"],dct["alt"],dct["phi"],dct["tx"],dct["ty"],dct["tz"]))
						elif isinstance(v,list) or isinstance(v,tuple):
							out.write([str(i) for i in v].join("\t"))
						else:
							out.write("{}".format(v))
					else: print(options.extractkey," not present in ",fsp)
					if k == options.extractkey[-1] : out.write("\n")
					else: out.write("\t")

	nimgs=0
	for imagefile in args:
		if options.quality!=-1 or options.dfmin!=None or options.dfmax!=None:
			try: 
				db=js_open_dict(info_name(imagefile))
				if options.quality!=-1 : 
					q=db["quality"]
				db.close()
				if q!=options.quality : continue
			except:
				sys.stderr.write("No quality for {}. Including\n".format(imagefile))
				
		if options.nameonly :
			print(imagefile)
			continue
		
		if imagefile.lower()[:4]!="bdb:" :
			try: nimg = EMUtil.get_image_count(imagefile)
			except:
				print(imagefile," is not a recognized image format")
				continue
			nimgs+=nimg
			imgtype = EMUtil.get_image_type(imagefile)
			imgtypename = EMUtil.get_imagetype_name(imgtype)
			if imgtype == EMUtil.ImageType.IMAGE_SPIDER and not options.stat: 
				image_index = -1
				
			d=EMData()
			try: d.read_image(imagefile, max(options.number,0), True)
			except :
				print("Image read error (%s)"%imagefile)
				continue
			if d["HostEndian"]!=d["ImageEndian"]: swapped="(swap)"
			else: swapped = ""
			if options.count : print("%d\t"%(nimg), end=' ')
			if d["nz"]==1 : print("%s\t %d images in %s format %s\t%d x %d"%(imagefile,nimg,imgtypename,swapped,d["nx"],d["ny"]))
			else : print("%s\t %d images in %s format %s\t%d x %d x %d"%(imagefile,nimg,imgtypename,swapped,d["nx"],d["ny"],d["nz"]))
			
		else: 
			print("ERROR: BDB no longer supported")
			sys.exit(1)
			
#		else :
#			dct=db_open_dict(imagefile,ro=True)
#			d=dct.get_header(0)
#			nimg=len(dct)
#			nimgs+=nimg
#			
#			d=EMData(imagefile, max(options.number,0), True)
#			if d["nz"]==1 : print("%s\t %d images in BDB format\t%d x %d"%(imagefile,len(dct),d["nx"],d["ny"]))
#			else : print("%s\t %d images in BDB format\t%d x %d x %d"%(imagefile,len(dct),d["nx"],d["ny"],d["nz"]))

		if options.all or options.check:
			imgn = list(range(nimg))
		elif options.number >= 0:
			imgn = [options.number]
		elif options.header or options.stat or options.euler:
			# If we request header information, without specifying an image #, assume the first image.
			imgn = [0]
		else:
			imgn = []

		nptcl=0
		d=EMData()
		for i in imgn:
			if options.stat : d.read_image(imagefile,i)
			else : 
				try: d.read_image(imagefile,i,True) #Jesus
				except: print("Error reading image ",imagefile,i)
			if options.check:
				try: 
					ctf=d["ctf"]	
					snr=ctf.snr[0]
				except:
					print("Error with CTF on ",imagefile,i)
			else:
				print("%d. "%i, end='')
				try:
					print("%d ptcl\t"%d["ptcl_repr"], end=' ')
					nptcl+=d["ptcl_repr"]
				except:
					print("\t",end='')
			
			if options.stat :
				print("apix=%-5.2f min=%-10.4g max=%-10.4g mean=%-10.4g sigma=%-9.4g skewness=%-9.4g kurtosis=%-9.4g moment_inertia=%9.4g radius_gyration=%9.4g"%(d["apix_x"],d["minimum"],d["maximum"],d["mean"],d["sigma"],d["skewness"],d["kurtosis"],d["moment_inertia"],d["radius_gyration"]), end=' ')
				if d.has_attr("stored_rendermin"):
					print(f"\n\tstored_rendermin={d['stored_rendermin']:-10.4g} stored_rendermax={d['stored_rendermax']:-10.4g} stored_renderbits={d['stored_renderbits']:-10.4g} stored_truncated={int(d['stored_truncated']):d} ({100.0*d['stored_truncated']/(d['nx']*d['ny']*d['nz']):1.6f}%) ",end=' ')
				if options.outliers:
					print("\nOutliers: {:^5.0f} {:^5.0f} {:^5.0f} {:^5.0f} {:^5.0f} {:^5.0f} {:^5.0f} {:^5.0f} {:^5.0f} {:^5.0f} {:^5.0f}".format(*d.calc_hist(11,d["mean"]-5.5*d["sigma"],d["mean"]+5.5*d["sigma"])))
				try:
					c=d["ctf"]
					print("   defocus=%-6.2f B=%-1.0f"%(c.defocus,c.bfactor))
				except:
					print(" ")
						
			if options.euler:
#				d=EMData(imagefile,i,True) #Jesus
				try: print("%s"%(str(d["xform.projection"])))
				except : print("No transform information", end=' ')

			if options.header :
#				d=EMData(imagefile,i, True) #Jesus
				print("")
				keys=list(d.get_attr_dict().keys())
				keys.sort()
				for k in keys:
					print("\t%s: %s"%(k,str(d[k])))
				print("======================")
			
	
	if nimgs>1 : print("%d total images"%nimgs)
	try : 
		if nptcl>0 : print("representing %d particles"%nptcl)
	except: pass
	if options.outliers : print("Outlier centers are mean+n*sigma: -5 -4 -3 -2 -1 0 1 2 3 4 5") 
	
if __name__ == "__main__":
	main()
