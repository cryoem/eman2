#!/usr/bin/env python

#
# Author: Steven Ludtke, 4/2/2010 
# Copyright (c) 2010 Baylor College of Medicine
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


from EMAN2 import *
from EMAN2db import db_open_dict
from math import *
import os
import sys
import datetime


def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options] 
	This program is still in its early stages. Eventually will provide a variety of tools for 
	evaluating a single particle reconstruction refinement run. Currently only provides a 
	single option to compare the parameters used for different refinement runs in a single project."""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	
	parser.add_argument("--timing", default=False, action="store_true", help="report on how long each step of the refinement process took during the first iteration of each run")
	parser.add_argument("--parmcmp",  default=False, action="store_true",help="Compare parameters used in different refinement rounds")
	parser.add_argument("--parmpair",default=None,type=str,help="Specify iter,iter to compare the parameters used between 2 itertions.")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	#options associated with e2refine.py
	#parser.add_argument("--iter", dest = "iter", type = int, default=0, help = "The total number of refinement iterations to perform")
	#parser.add_argument("--check", "-c", dest="check", default=False, action="store_true",help="Checks the contents of the current directory to verify that e2refine.py command will work - checks for the existence of the necessary starting files and checks their dimensions. Performs no work ")
	#parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	#parser.add_argument("--input", dest="input", default=None,type=str, help="The name of the image containing the particle data")

	(options, args) = parser.parse_args()

	if options.parmpair :
		try: 
			na,nb=options.parmpair.split(",")
			na=int(na)
			nb=int(nb)
		except:
			print "Please specify <iter1>,<iter2>"
			sys.exit(1)
		db=db_open_dict("bdb:refine_%02d#register"%na,ro=True)
		pa=db["cmd_dict"]
		db=db_open_dict("bdb:refine_%02d#register"%nb,ro=True)
		pb=db["cmd_dict"]

		ks=set(pa.keys())
		ks|=set(pb.keys())
		ks=list(ks)
		ks.sort()

		for k in ks:
			try :
				if pa[k]!=pb[k] : print "%s: %s -> %s"%(k,pa[k],pb[k])
			except:
				try: print "%s: %s -> None"%(k,pa[k])
				except: print "%s: None -> %s"%(k,pb[k])

	if options.timing:
		#dl=[i for i in os.listdir(".") if "refine_" in i]		# list of all refine_ directories
		#dl.sort()

		hist=[]
		fin=open(".eman2log.txt","r")
		while 1:
			line=fin.readline()
			if len(line)==0 : break
			
			spl=line.split("\t")
			try: com=spl[4].split()[0].split("/")[-1]
			except : continue

			if com in ("e2refine.py","e2refine_easy.py","e2refinemulti.py","e2project3d.py","e2simmx.py","e2simmx2stage.py","e2classify.py","e2classaverage.py","e2make3d.py") : hist.append((com,spl))
			
		n=0
		while n<len(hist):
			com=hist[n][1][4]
			ttime=timestamp_diff(hist[n][1][0],hist[n][1][1])
			if hist[n][0] in ("e2refine.py","e2refine_easy.py"): 
				pl=com.find("--path=")
				parl=com.find("--parallel=")
				print "%s\t%1.2f hours\te2refine %s"%(difftime(ttime),ttime/3600.0,com[pl+7:].split()[0]),
				if parl>0: print com[parl+11:].split()[0]
				else: print " "

			else:
				print "\t%s\t%1.2f hours\t%s"%(difftime(ttime),ttime/3600.0,hist[n][0])
			
			n+=1
			

	if options.parmcmp:
		dl=[i for i in os.listdir(".") if "refine_" in i]		# list of all refine_ directories
		dl.sort()
		
		# extract an array of various parameters
		parmlist=[(" ","it","shr","sep","2s","simcmp","simr","simrcmp","clsit","clscmp","clsr","clsrcmp","csf","3sf","3post","initmdl","sym")]
		for d in dl:
			if not os.path.isdir(d) : continue
			try :
				db=db_open_dict("bdb:%s#register"%d,ro=True)
				parms=db["cmd_dict"]
			except :
				print d," missing refinement log"
				continue
			try:
				p=[d,parms["iter"],parms.get("shrink",1),parms["sep"],parms.get("twostage",0),parms["simcmp"],str(parms["simralign"])[:1],parms["simraligncmp"],parms["classiter"],
					parms["classcmp"],str(parms["classralign"])[:1],parms["classraligncmp"],parms["classrefsf"][-2:-1],parms["m3dsetsf"],parms["m3dpostprocess"],parms["model"],parms["sym"]]
			except:
				print d, "incomplete parameters"
				continue
		
			if p[2]==None or p[2]=="None": p[2]=1
			
			parmlist.append(tuple(p))
		# find the max length of each element in the array
		nparm=len(parmlist[0])
		lens=[0]*nparm
		for i in xrange(nparm):
			for p in parmlist:
				try : lens[i]=max(lens[i],len(str(p[i])))
				except : print "Error with len of ",i,p[i]

		# build a format string
		fmt="%%-%0ds "%lens[0]
		totlen=lens[0]+1
		for i in xrange(1,nparm): 
			if totlen+lens[i]>120 and lens[i]+lens[0]+1<120 :
				fmt+="\n          "
				totlen=10
			fmt+="%%-%0ds "%lens[i]
			totlen+=lens[i]+1
			
		# print results
		for p in parmlist:
			print fmt%p

		
if __name__ == "__main__":
    main()
