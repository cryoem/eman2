#!/usr/bin/env python

#
# Author: Steven Ludtke, 11/13/2008 (sludtke@bcm.edu)
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

# e2bdb.py  11/13/2008 Steven Ludtke
# This program allows manipulation and querying of the local database

from math import *
import time
import os
import sys
import re
import traceback

from EMAN2 import EMAN2DB, EMUtil, EMANVERSION
from EMAN2db import db_open_dict, db_list_dicts, db_remove_dict, e2gethome
from EMAN2 import *

def main():
	global debug
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options] <path or db> ...
	
Various utilities related to BDB databases.

examples :
e2bdb.py -c   Is perhaps the most critical function, as it cleans up the database cache. See the Wiki for more.
e2bdb.py <path> -s    will list the contents of the database in a directory in bdb: notation
e2bdb.py <path> -l    Will give useful summary info about stacks in a directory
e2bdb.py <database> --dump    Gives a mechanism to dump all of the metadata in a database, even if the database contains no images
"""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_argument("--cleanup","-c",action="store_true",default=False,help="This option will clean up the database cache so files can safely be moved or accessed on another computer via NFS.")
	parser.add_argument("--force","-F",action="store_true",default=False,help="This will force an action that would normally fail due to failed checks.")
	parser.add_argument("--delete",action="store_true",default=False,help="This will delete (or at least empty) the named database(s)")
	parser.add_argument("--all","-a",action="store_true",help="List per-particle info",default=False)
	parser.add_argument("--long","-l",action="store_true",help="Long listing",default=False)
	parser.add_argument("--short","-s",action="store_true",help="Dense listing of names only",default=False)
	parser.add_argument("--filt",type=str,help="Only include dictionary names containing the specified string",default=None)
	parser.add_argument("--filtexclude",type=str,help="Exclude dictionary names containing the specified string",default=None)
	parser.add_argument("--match",type=str,help="Only include dictionaries matching the provided Python regular expression",default=None)
	parser.add_argument("--exclude",type=str,help="The name of a database containing a list of exclusion keys",default=None)
	parser.add_argument("--dump","-D",action="store_true",help="List contents of an entire database, eg 'e2bdb.py -D refine_01#register",default=False)
	parser.add_argument("--smalldump",action="store_true",help="Lists contents of an entire database, but only list 2 items per dictionary to better see headers",default=False)
	parser.add_argument("--extractplots",action="store_true",help="If a database contains sets of plots, such as bdb:refine_xx#convergence.results, this will extract the plots as text files.")
	parser.add_argument("--check",action="store_true",help="Check for self-consistency and errors in the structure of specified databases",default=False)
	parser.add_argument("--nocache",action="store_true",help="Don't use the database cache for this operation",default=False)
	parser.add_argument("--merge",action="store_true",help="This will merge the contents of BDB 2-N into BDB 1 (including BDB 1's contents)",default=False)


	parser.add_argument("--makevstack",type=str,help="Creates a 'virtual' BDB stack with its own metadata, but the binary data taken from the (filtered) list of stacks",default=None)
	parser.add_argument("--appendvstack",type=str,help="Appends to/creates a 'virtual' BDB stack with its own metadata, but the binary data taken from the (filtered) list of stacks",default=None)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higher number means higher level of verboseness")
	parser.add_argument("--list",type=str,help="Specify the name of a file with a list of images to use in creation of virtual stacks. Please see source for details.",default=None)
	parser.add_argument("--exlist",type=str,help="Specify the name of a file with a list of images to exclude in creation of virtual stacks. Please see source for details.",default=None)
	parser.add_argument("--restore",type=str,help="Write changes in the derived virtual stack back to the original stack",default=None)
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--checkctf",action="store_true",help="Verfies that all images in the file contain CTF information, and gives some basic statistics",default=False)

	parser.add_argument("--step",type=str,default="0,1",help="Specify <init>,<step>[,<max>]. Processes only a subset of the input data. For example, 0,2 would process only the even numbered particles")
	(options, args) = parser.parse_args()

	if options.nocache : EMAN2db.BDB_CACHE_DISABLE=True

	if options.cleanup : 
		db_cleanup(options.force)
		sys.exit(0)

	try : options.step=int(options.step.split(",")[0]),int(options.step.split(",")[1]),int(options.step.split(",")[2])		# convert strings to tuple
	except:
		try:
			options.step=int(options.step.split(",")[0]),int(options.step.split(",")[1])
		except:
			print "Invalid --step specification"
			sys.exit(1)

	if options.all : options.long=1
	if len(args)==0 : args.append("bdb:.")
	
	logid=0
	if options.makevstack : 
		logid=E2init(sys.argv)
		vstack=db_open_dict(options.makevstack)
		vstackn=0
	elif options.appendvstack :
		logid=E2init(sys.argv)
		vstack=db_open_dict(options.appendvstack)
		vstackn=len(vstack)
	else : vstack=None
		
	if options.merge :
		print "WARNING: Merge mode\nCombining contents of: ",", ".join(args[1:])
		print "into ",args[0]
		
		if raw_input("Proceed (y/n) :").lower() != "y" :
			print "Aborting"
			sys.exit(1)
		
		
		for i,path in enumerate(args):
			if path.lower()[:4]=="bdb:" and not "#" in path : path="bdb:.#"+path[4:]
			if path.lower()[:4]!="bdb:" : path="bdb:"+path
			
			if i==0 :
				outdb=db_open_dict(path)
				continue
			
			indb=db_open_dict(path,True)
			for k in indb.keys():
				outdb[k]=indb[k]
				
		print "Merging complete"
		sys.exit(0)
	
	for path in args:
		if path.lower()[:4]=="bdb:" and not "#" in path :
			uu = os.path.split(path)
			if(uu[0] == ''):    path="bdb:.#"+path[4:]
			else:               path=uu[0]+"#"+uu[1]
		if path.lower()[:4]!="bdb:" : path="bdb:"+path
		if '#' in path :
			if len(args)>1 : print "\n",path,":"
			path,dbs=path.rsplit("#",1)
			path+="#"
			dbs=[dbs]
		else:
			if not '#' in path and path[-1]!='/' : path+='#'			
			if len(args)>1 : print "\n",path[:-1],":"
			dbs=db_list_dicts(path)
			

		dbs.sort()
		if options.filt:
			dbs=[db for db in dbs if options.filt in db]
			
		if options.filtexclude:
			dbs=[db for db in dbs if options.filtexclude not in db]

		if options.match!=None:
			dbs=[db for db in dbs if re.match(options.match,db)]
		
		
		if options.list :
			if options.makevstack==None and options.appendvstack==None :
				print "ERROR, this option is used for virtual stack creation, please add makevstack or appendvstack options, and restart"
				sys.exit(1)
			vdata=open(options.list,'r').readlines()
			n=len(vdata[0].split())
			slist=[]
			for line in vdata:
				line=line.split()
				for i in xrange(n):
					val=int(line[i])
					slist.append(val)     		
			del n,val,vdata
		
		if options.exlist :
			if options.makevstack==None:
				print "ERROR, this option is used for virtual stack creation, please add makevstack or appendvstack options, and restart"
				sys.exit(1)
			vdata=open(options.exlist,'r').readlines()
			n=len(vdata[0].split())
			slist=[]
			for line in vdata:
				line=line.split()
				for i in xrange(n):
					val=int(line[i])
					slist.append(val)     
			n = EMUtil.get_image_count(args[0])
			good = set(range(n)) - set(slist)
			slist = [i for i in good]
			slist.sort()
			del n,val,vdata,good
			
		if options.makevstack!=None or options.appendvstack!=None :
			
			vspath=os.path.realpath(vstack.path)+"/"
			if options.verbose>2 : print "vspath: ",vspath
			for db in dbs:
				dct,keys=db_open_dict(path+db,ro=True,with_keys=True)
				if dct==vstack : continue
				if len(options.step)==2 :
					if keys == None: vals = xrange(options.step[0],len(dct),options.step[1])
					else: vals = keys[options.step[0]::options.step[1]]		# we apply --step even if we have a list of keys
				else:
					if keys == None: vals = xrange(options.step[0],options.step[2],options.step[1])
					else: vals = keys[options.step[0]:options.step[2]:options.step[1]]		# we apply --step even if we have a list of keys

				if options.list !=None or options.exlist != None: vals=slist
				for n in vals:
					try: d=dct.get(n,nodata=1).get_attr_dict()
					except:
						traceback.print_exc()
						print "---\nerror reading ",db,n 
						continue
					# This block converts an absolute path to the actual data to a relative path
					try:
						dpath=os.path.realpath(dct.get_data_path(n))
						if options.verbose>2 : print "dpath: ",dpath
						if os.name == 'nt':
							vspath=vspath.replace("\\", '/')
							dpath=dpath.replace('\\', '/')
						rpath=makerelpath(vspath,dpath)
						if options.verbose>2 : print "rpath: ",rpath
					except:
						print "error with data_path ",db,n
						continue
					d["data_path"]=rpath
					d["data_n"]=n
					d["data_source"]= path+db
					if d["data_path"]==None :
						print "error with data_path ",db,n
						continue
					vstack[vstackn]=d
					vstackn+=1
					if vstackn%100==0:
						try:
							print "\r  ",vstackn,"     ",
							sys.stdout.flush()
						except: pass	
				print "\r  ",vstackn,"     "
				dct.close()

		try: maxname=max([len(s) for s in dbs])
		except:
			print "Error reading ",path

		if options.restore :
			nima = EMUtil.get_image_count(options.restore)
			IB = db_open_dict(options.restore)
			source_old = None
			if len(options.step)==3 : nima=min(options.step[2],nima)
			for i in xrange(options.step[0],nima,options.step[1]):
				source = IB.get_header(i)
				source_path = source["source_path"]
				ID = source["source_n"]
				if( source_old != source_path):
					if( source_old != None):  DB.close()
					DB = db_open_dict(source_path,ro=True)
					source_old = source_path
				target = DB.get_header( ID )
				try:
					source["data_path"] = target["data_path"]
					source["data_n"]    = target["data_n"]
					source["source_path"] = target["source_path"]
					source["source_n"]    = target["source_n"]
				except:
					#  top level does not have data_path
					del source['data_path']
					del source['data_n']
					source["source_path"] = target["source_path"]
					source["source_n"]    = target["source_n"]
				DB.set_header(ID, source)
			DB.close()

		if options.extractplots :
			for db in dbs:
				print "####  Extracting plots from ",db
				dct=db_open_dict(path+db,ro=True)
				
				#### Dump
				keys=dct.keys()
				keys.sort()
				for k in keys:
					v=dct[k]
					try :
						ns=[len(i) for i in v]
						fsp=db+"-"+k+".txt"
						print "%s  (%d columns)"%(fsp,len(ns))
						out=file(fsp,"w")
						for i in range(ns[0]):
							for j in range(len(ns)):
								out.write(str(v[j][i]))
								if j<len(ns)-1 : out.write("\t")
							out.write("\n")
						out.close()
					except: continue
				dct.close()
			
		if options.smalldump :
			for db in dbs:
				print "##### ",db
				dct=db_open_dict(path+db,ro=True)
				
				#### Dump
				keys=dct.keys()
				keys.sort()
				if len(options.step)==3 : keys=keys[:options.step[2]]
				for k in keys[options.step[0]::options.step[1]]:
					v=dct[k]
					print "%s : "%k,
					if isinstance (v,list) or isinstance(v,tuple)  :
						for i in v: print "\n\t%s"%str(i),
						print ""
					elif isinstance(v,dict) :
						ks2=v.keys()
						ks2.sort()
						kc=0
						for i in ks2:
							if kc>=2 :
								print "..."
								break
							print "\n\t%s : %s"%(i,v[i]),
							kc+=1
						print ""
					else : print str(v)
				dct.close()
		if options.checkctf:
			for db in dbs:
				print "##### CTF -> ",db
				dct=db_open_dict(path+db,ro=True)
				keys=dct.keys()
				if len(options.step)==3 : keys=keys[:options.step[2]]
				defocus=set()
				for k in keys[options.step[0]::options.step[1]]:
					v=dct.get_header(k)
					try:
						ctf=v["ctf"]
					except:
						if k!="maxrec" : print "CTF missing on image %s"%k
						continue
					
					defocus.add(ctf.defocus)

				defocus=list(defocus)
				print "Defocuses found: ",
				for i in defocus: print "%1.3f, "%i,
				print "\n\nRange: %1.3f - %1.3f  (%d unique values)"%(min(defocus),max(defocus),len(defocus))

					

		if options.dump :
			for db in dbs:
				print "##### ",db
				dct=db_open_dict(path+db,ro=True)
				
				#### Dump
				keys=dct.keys()
				if len(options.step)==3 : keys=keys[:options.step[2]]
				keys.sort()
				for k in keys[options.step[0]::options.step[1]]:
					v=dct[k]
					print "%s : "%k,
					if isinstance (v,list) or isinstance(v,tuple)  :
						for i in v: print "\n\t%s"%str(i),
						print ""
					elif isinstance(v,dict) :
						ks2=v.keys()
						ks2.sort()
						for i in ks2:
							print "\n\t%s : %s"%(i,v[i]),
						print ""
					else : print str(v)
				dct.close()
			
		# long listing, one db per line
		elif options.long :
			width=maxname+3
			fmt="%%-%ds %%-07d %%14s  %%s"%width
			fmt2="%%-%ds (not an image stack)"%width
			total=[0,0]
			for db in dbs:
				dct=db_open_dict(path+db,True)
								
				### Info on all particles
				if options.all :
					mx=len(dct)
					if len(options.step)==3 : mx=min(mx,options.step[2])
					for i in range(options.step[0],mx,options.step[1]):
						try: 
							im=dct[i]
							if im==None : raise Exception
						except: continue
						print "%d. %d x %d x %d\tA/pix=%1.2f\tM=%1.4f\tS=%1.4f\tSk=%1.4f"%(i,im["nx"],im["ny"],im["nz"],im["apix_x"],im["mean"],im["sigma"],im["skewness"]),
						try: print "\t%s"%str(im["model_id"])
						except: pass
						try:
							print "\tdf=%1.3f\tB=%1.1f"%(im["ctf"].defocus,im["ctf"].bfactor)
						except: print " "

				first=EMData()
				try: 
					first.read_image(path+db,0,True)
					size=first.get_xsize()*first.get_ysize()*first.get_zsize()*len(dct)*4;
					total[0]+=len(dct)
					total[1]+=size
					print fmt%(db,len(dct),"%dx%dx%d   apix: %1.2f"%(first.get_xsize(),first.get_ysize(),first.get_zsize(),first["apix_x"]),human_size(size)),
				except:
					print fmt2%db
				try: print "\tdf: %1.3f\tB: %1.0f"%(first["ctf"].defocus,first["ctf"].bfactor)
				except: print ""
				dct.close()
			print fmt%("TOTAL",total[0],"",human_size(total[1]))
		elif options.check :
			from cPickle import loads
			for db in dbs:
				dct=db_open_dict(path+db,ro=True)
				dct.realopen()
				keys=dct.bdb.keys()
				allkvp={}
				for k in keys:
					s1,s2=k.split("\x80",1)		# start of a pickled string. 
					s2=loads("\x80"+s2)			# the pickled part
					if len(s1)>0 :				# If anything unpickled, then it is an axbxc prefix identifying the location of a binary
						st=allkvp.setdefault(s1,set()) # set of all positions seen so far
						v=loads(dct.bdb.get(k))	# position in binary file
						if v in st : print "Error: value %d seen multiple times in %s (%s,%s)"%(v,db,s1,s2)
						st.add(v)
				print "%s : "%db,
				for i in allkvp.keys(): 
					if options.verbose>0 : print "%s %d/%d\t"%(i,len(allkvp[i]),int(max(allkvp[i]))+1),
					if len(allkvp[i])!=int(max(allkvp[i])+1) : print "\nMismatch found in %s. Could be normal if file has been rewritten multiple times, but is unusual"%db
				if options.verbose>0 : print ""
				else : print " done"
				dct.close()

		elif options.short :
			for db in dbs:
				print path+db,
			print " "

		elif not options.makevstack and not options.appendvstack :
			# Nicely formatted 'ls' style display
			cols=int(floor(80.0/(maxname+3)))
			width=80/cols
			rows=int(ceil(float(len(dbs))/cols))
			
			fmt="%%-%ds"%width
			for r in range(rows):
				for c in range(cols):
					try: print fmt%dbs[r+c*rows],
					except: pass
				print " "

		if options.delete :
			if not options.force :
				print "You are requesting to delete the following databases:"
				for db in dbs:
					print db," ",
				if raw_input("\nAre you sure (y/n) ? ")[0].lower()!='y' :
					print "Aborted"
					sys.exit(1)
			
			for db in dbs: db_remove_dict(path+db)
			

	if logid : E2end(logid)

def makerelpath(p1,p2):
	"""Takes a pair of paths /a/b/c/d and /a/b/e/f/g and returns a relative path to b from a, ../../e/f/g"""
	
	p1s=[i for i in p1.split("/") if len(i)>0]
	p2s=[i for i in p2.split("/") if len(i)>0]

	for dv in range(min(len(p1s),len(p2s))):
		if p1s[dv]!=p2s[dv] : break
	else: dv+=1

	p1s=p1s[dv:]
	p2s=p2s[dv:]
	
	return "../"*len(p1s)+"/".join(p2s)

def db_cleanup(force=False):
	"""This is an important utility function to clean up the database environment so databases can safely be moved or used remotely
	from other machines. If working on a cluster, this routine should be called on any machine which has opened a database before that
	database is written to on another machine"""
	
	if(sys.platform == 'win32'):
		force = True
		path="eman2db-%s"%os.getenv('USERNAME')
		print "Database cleanup is in force mode on windows machines"
	else:		
		path="eman2db-%s"%os.getenv("USER","anyone")
	
	if not force :
		try:
			pipe=os.popen("lsof","r")
			op=[l.split() for l in pipe if path in l]
			ret=pipe.close()
			if ret!=None :
				pipe=os.popen("/usr/sbin/lsof","r")
				op=[l.split() for l in pipe if path in l]
				ret=pipe.close()
				if ret!=None: raise Exception
		except:
			print "Error : could not check for running EMAN2 jobs, please make sure the 'lsof' command is installed and functioning, or insure no EMAN2 commands are running and run e2bdb.py -cF"
			sys.exit(1)
		
		# someone is still using the cache
		if len(op)>0 :
			s=set()
			for i in op:
				s.add(i[1])
		
			print "These processes are actively using the cache. Please exit them and try again :"
			for i in s: 
				try: print os.popen("ps %s"%i,"r").readlines()[-1]
				except: print i
			
			reply=raw_input("Would you like me to kill all of these jobs (YES/NO) : ")
			if reply != "YES" : 
				print "Not killing jobs. Please exit them manually then retry."
				return

			for i in s: os.kill(int(i),15)
			print "Signal sent to kill job(s). Please retry e2bdb.py -c"
			return

	# ok, properly close the cache and delete it
	try:
		d=EMAN2DB()
		d.close()		# Properly 'close' the environment before we delete it
	except:
		traceback.print_exc()
		print """
*******************

A serious error occured in the database cache. This normally happens if you try to access a corrupt database file. Please follow the following steps to minimize the chance of data loss:
1. run "db_recover -h %s"   (note, this may be called db4.8_recover or something similar if not using the EMAN2 binary distribution)
2. run e2bdb.py -c again
3. If you are aware which image file caused this error to occur in the first place, you can try accessing it again. If it triggers this same failure, repeat steps 1 and 2 then manually delete the offending image database inside the EMAN2DB directory"""%path
		sys.exit(1)
	
	if(sys.platform == 'win32'):
		import shutil
		shutil.rmtree('C:'+e2gethome()+'/.eman2/EMAN2DB')
		shutil.rmtree('C:/tmp/'+path)
	else:
		os.system("rm -rf /tmp/%s"%path)
	print "Database cache removed. Now safe to access databases from another machine or delete existing databases"
	
	

def human_size(size):
	if size>1000000000: return "%1.2f gb"%(size/1000000000)
	elif size>1000000: return "%1.2f mb"%(size/1000000)
	else: return "%1.2f kb"%(size/1000)
	return str(size)
			
if __name__ == "__main__":
	main()
