#!/usr/bin/env python
#
# Author: John Flanagan (jfflanag@bcm.edu)
# Copyright (c) 2000-2011 Baylor College of Medicine


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
import os, re
from EMAN2 import *
import numpy as np
from EMAN2_utils import natural_sort

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options] --stackname myfile.hdf <image1> <image2> <image3> ...
	This program will combine many image files into a single output file. 
	
	If the output name has a ".lst" extension:
	the output is a formatted text file, one line per image, describing the file containing the actual
	image data in a searchable form. .lst files can be used as if they contained actual images in any
	EMAN2 programs.
	
	If the output is a normal image file (.hdf, .spi, etc.) then the images will be copied into the
	output file sequentially in the order provided on the command-line. Some file formats will not
	support multiple images, or multiple volumes. Appropriate errors will be raised in these cases.
	HDF is the only format supporting full metadata retention for stacks of images or volumes.
	
	The output file will be emptied and overwritten!
	"""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_pos_argument(name="stack_files",help="List of images to be stacked. Selecting a folder to use all images inside.", default="micrographs", guitype='filebox', browser="EMParticlesEditTable(withmodal=True,multiselect=True)",  row=0, col=0,rowspan=1, colspan=3, nosharedb=True,mode="default")
	parser.add_pos_argument(name="tilt_images",help="List of images to be stacked. Input order will determine the order of images in output tiltseries.", default="", guitype='filebox', browser="EMBrowserWidget(withmodal=True)",  row=0, col=0,rowspan=1, colspan=3, nosharedb=True,mode="tomo")
	parser.add_argument("--output",type=str,help="Name of the output stack to build (Extension will be .hdf unless specified). Note, all tiltseries will be stored in the 'tiltseries' directory.", default=None, guitype='strbox',row=2, col=0, rowspan=1, colspan=1, mode="default,tomo")
	parser.add_argument("--tilts",action="store_true",default=False,help="Write results to 'tiltseries' directory in current project.", guitype='boolbox',row=4, col=0, rowspan=1, colspan=1,mode="tomo[True]")
	parser.add_argument("--guess",action="store_true",default=False,help="Guess how to split micrographs into tilt series and the order of images in each tilt series from file names. Tilt angles must be incuded in file names. May and may not work depending on the file name format...", guitype='boolbox',row=4, col=1, rowspan=1, colspan=1,mode="tomo[False]")	
	parser.add_argument("--guesscol", type=int, help="column to separate tilt series if the guess mode fails",default=-1)
	parser.add_argument("--mdoc", type=str, help="Read info from mdoc files",default=None)

	#parser.add_argument("--rawtlt",type=str,help="Name of tilt angles text file.\nNote, angles must correspond to stack file names in alphabetical/numerical order.", default="", guitype='filebox', browser="EMBrowserWidget(withmodal=True,multiselect=True)",row=3, col=0, rowspan=1, colspan=1, mode="tomo")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, help="verbose level [0-9], higher number means higher level of verboseness",default=1)
	parser.add_argument("--tltang",type=str,help="a text file for tilt angle. will sort the images in each stack accrodingly", default=None)

	(options, args) = parser.parse_args()

	if options.output==None:
		if options.guess or options.mdoc:
			pass
		else:
			print("--output is required (output file)")
			sys.exit(1)

	if len(args)==1 and os.path.isdir(args[0]):
		print("input is a directory. reading all mrc/mrcs/hdf files in it...")
		path=args[0]
		args=[]
		ext=["mrc", "mrcs", "hdf"]
		for f in os.listdir(path):
			for e in ext:
				if f.endswith(e):
					args.append(os.path.join(path, f))
		print("found {} files".format(len(args)))
		
	if options.tilts:
			
		try:
			os.mkdir("tiltseries")
		except:
			pass
		
		if options.guess:	

			lst=[]
			lstpos=[]
			for ag in args:
				l=[]
				s0=""
				p=[0]
				a=ag.replace(".", "_")
				for i,c in enumerate(a[:-1]):
					if c.isdigit() or (s0=="" and c=='-' and a[i+1].isdigit()):
						s0+=c
					elif len(s0)>0:
						l.append(int(s0))
						p.append(i)
						s0=""
				lst.append(l)
				lstpos.append(p)
			
			l0=len(lst[0])
			for i,l in enumerate(lst):
				if len(l)!=l0:
					print("File name mismatch detected")
					print("{}\n\t-> {}".format(args[0], lst[0]))
					print("{}\n\t-> {}".format(args[i], l))
					break
			try:
				lst=np.array(lst, dtype=float)
			except:
				print("Something is wrong in filename formatting. exit.")
				return
				
			print("File name of the first input:")
			print("\t{}".format(args[0]))
			print("{} Columns".format(len(lst[0])))
			print("\t"+',  '.join(lst[0].astype(int).astype(str)))
			dt=[]
			for i in range(len(lst[0])):
				mn,mx=np.min(lst[:,i]), np.max(lst[:,i])
				#print("{}: range from {} to {}".format(i, mn, mx))
				dt.append(abs(mn+60)+abs(mx-60))
			
			ic=np.argmin(dt)
			print("Guess column {} is for tilt angles,\n\tranging from {:.1f} to {:.1f}.".format(ic, np.min(lst[:,ic]), np.max(lst[:,ic])))
			
			
			c=lst[:, ic]
			if len(c)>len(np.unique(c)):
				print("Multiple tilt series exist...")
				if options.guesscol<0:
					it=np.where(np.std(lst,axis=0)>0)[0][0]
					print("Guess column {} separates different tilt series".format(it))
				else:
					it=options.guesscol
					print("Separates different tilt series using column {} ".format(it))
			
				fid=sorted(np.unique(lst[:,it]))
				print("\t{} files, from {:.0f} to {:.0f}.".format(len(fid), np.min(lst[:,it]), np.max(lst[:,it])))
				
				tlts=[np.where(lst[:,it]==t)[0] for t in fid]
			else:
				tlts=[np.arange(len(lst), dtype=int)]
				it=0
			
			for tid in tlts:
				l=lst[tid]
				#print(l[:,ic])
				if len(l)>len(np.unique(l[:,ic])):
					print("    duplicated tilt images exist...")
					aid=np.argsort(l[:,ic-1])
					tid2=[]
					aid2=[]
					for ii in aid:
						if l[ii, ic] not in aid2:
							aid2.append(l[ii, ic])
							tid2.append(tid[ii])
						else:
							aid2[-1]=l[ii,ic]
							tid2[-1]=tid[ii]
					tid=np.array(tid2, dtype=int)
					print("    keeping {} out of {} images".format(len(tid), len(l)))
					l=lst[tid]
					
				aid=np.argsort(l[:,ic])
				fnames=[args[i] for i in tid[aid]]
				p=lstpos[tid[0]][it+1]
				prefix=fnames[0][fnames[0].rfind('/')+1:p]
				prefix=prefix.replace("__", "_")
				prefix=prefix.replace(".", "_")
				
				lstname=os.path.join("tiltseries", prefix+'.lst')
				
				print("{} : {} images -> {}".format(prefix, len(fnames), lstname))
				if options.tltang:
					if os.path.isdir(options.tltang):
						tnames=os.listdir(options.tltang)
						ts=[t[:t.rfind('.')] for t in tnames]
						ts=[k for k,t in enumerate(ts) if t in prefix]
						if len(ts)==0:
							print("cannot fine tilt file")
							break
						elif len(ts)>1:
							print("multiple matching files exist")
							print([tnames[k] for k in ts])
							break
						
						tname=tnames[ts[0]]
						tname="{}/{}".format(options.tltang, tname)
						print("using tlt file",tname)
					else:
						tname=options.tltang
						
					ang=np.loadtxt(tname)
					if len(ang)==len(fnames):
						print("Sorting by tilt angle file")
						srt=np.argsort(ang)
						fnames=[fnames[i] for i in srt]
						ang=np.sort(ang)
					else:
						print("tilt file length mismatch",len(ang), len(fnames))
						break
				
				lout=[{"src":fm,"idx":0} for fm in fnames]
				if options.tltang:
					for li in range(len(lout)):
						lout[li]["tilt_angle"]=ang[li]
				
				save_lst_params(lout, lstname)
				

		elif options.mdoc:
			print("Parsing mdoc files")
			if os.path.isdir(options.mdoc):
				print("input is a directory. reading all mdoc files in it...")
				mdoc=[os.path.join(options.mdoc, f) for f in os.listdir(options.mdoc) if f.endswith(".mdoc")]
			else:
				mdoc=[options.mdoc]
				
			for md in mdoc:
				print(md)
				f=open(md, 'r')
				lines=f.readlines()
				fnames=[l for l in lines if l.startswith("SubFramePath")]
				ang=[l for l in lines if l.startswith("TiltAngle")]
				ang=[float(l.split('=')[-1]) for l in ang]
				srt=np.argsort(ang)
				ang=np.sort(ang)
				lst=[]
				for l in fnames:
					p0=max(l.rfind('/'), l.rfind('\\'))+1
					p1=l.rfind('.')
					l=l[p0:p1]
					l1=l.replace('.', '_')
					match=[a for a in args if l in a]
					match+=[a for a in args if l1 in a]
					if len(match)==0:
						print("error: image file for {} does not exist".format(l))
						break
					if len(match)>1:
						print("error: multiple images for {} exist".format(l))
						for a in match:
							print('\t',a)
						break
					
					
					lst.append({"src":match[0],"idx":0})
				
				else:
					lst=[lst[i] for i in srt]
					for i,l in enumerate(lst):
						l["tilt_angle"]=ang[i]
					tname=md[md.rfind('/')+1:]
					tname=tname[:tname.find('.')]+".lst"
					tname=os.path.join("tiltseries", tname)
					print(f"{len(lst)} images -> {tname}")
					save_lst_params(lst, tname)
				
				f.close()
			
		else:
			stdir = os.path.join(".","tiltseries")
			options.output = "{}/{}".format(stdir,options.output)

			if options.output.split(".")[-1] not in ["hdf","mrc","mrcs"]:
				options.output = options.output + ".hdf"

			# remove existing output file
			if os.path.exists(options.output) :
				print("The file {} already exists.".format(options.output))
				print("Please move, rename, or remove this file to generate an alternate version with this program.")
				sys.exit(1)

			
			for n,arg in enumerate(args):
				img = EMData(arg)
				img.write_image(options.output,n)
			
			
		# if options.rawtlt:
		# 	try:
		# 		angles = np.loadtxt(options.rawtlt)
		# 	except:
		# 		print("Error: Could not read tilt angles from {}".format(options.rawtlt))
		# 		sys.exit(1)
		# 	if len(angles) != len(args):
		# 		print("Error: There are not enough tilt angles in this tilt angles file.")
		# 		sys.exit(1)

		# tlt_assoc = {}
		# for i,arg in enumerate(args):
		# 	if options.rawtlt: tlt_assoc[angles[i]] = arg
		# 	else:
		# 		db=js_open_dict(info_name(arg,nodir=True))
		# 		ang = float(db["tilt_angle"])
		# 		tlt_assoc[ang] = arg
		# 		db.close()

		#ordered_angles = sorted([float(a) for a in tlt_assoc.keys()])
		#sorted_args = [tlt_assoc[a] for a in ordered_angles] # order args according to tilt angle parameter

		#series_db=js_open_dict(info_name(options.output,nodir=True))

		#series_db["tilt_angles"] = ordered_angles

		#for n,(angle,arg) in enumerate(zip(ordered_angles,sorted_args)):

			#series_db[angle] = arg

			#nimg = EMUtil.get_image_count(arg) # number of images in each input file as it is processed

			# if options.verbose:
			# 	if nimg==1: print(arg)
			# 	else: print(arg,nimg)

			#for i in xrange(nimg):

			#img=EMData(arg,0)
			#img["tilt_angle"] = angle

			# if os.path.isfile(info_name(arg,nodir=True)):
			# 	db=js_open_dict(info_name(arg,nodir=True))
			# 	try: # this data may already be present
			# 		img["SerialEM.tilt_angle"] = db["tilt_angle"]
			# 		img["SerialEM.intensity"] = db["intensity"]
			# 		img["SerialEM.exposure_time"] = db["exposure_time"]
			# 		img["SerialEM.exposure_dose"] = db["exposure_dose"]
			# 		img["SerialEM.sub_frame_count"] = db["sub_frame_count"]
			# 		img["SerialEM.prior_record_dose"] = db["prior_record_dose"]
			# 		img["SerialEM.frames_per_second"] = db["frames_per_second"]
			# 	except: pass
			# 	db.close()

			#img.write_image(options.output,n)

		#series_db.close()

	else:

		# remove existing output file
		if os.path.exists(options.output) :
			try: os.unlink(options.output)
			except:
				print("ERROR: Unable to remove ",options.output,". Cannot proceed")
				sys.exit(1)

		# if output is LSX format, we handle it differently, with a specific object for these files
		if options.output[-4:].lower()==".lst" :
			outfile=LSXFile(options.output)
		else: outfile=None

		n=0		# number of images in output file
		for infile in args:
			nimg = EMUtil.get_image_count(infile)		# number of images in each input file as it is processed

			if options.verbose :
				if nimg==1 : print(infile)
				else : print(infile,nimg)

			for i in range(nimg):
				if outfile!=None:
					outfile.write(n,i,infile)
				else:
					img=EMData(infile,i)
					img.write_image(options.output,n)
				n+=1

		if options.verbose : print(n," total images written to ",options.output)


if __name__ == "__main__":
	main()
