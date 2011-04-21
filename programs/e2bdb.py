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

from optparse import OptionParser
from math import *
import time
import os
import sys
import re
import traceback

from EMAN2 import EMAN2DB, EMUtil, EMANVERSION
from EMAN2db import db_open_dict, db_list_dicts, db_remove_dict, e2gethome
from EMAN2 import E2init, E2end

def main():
	global debug
	progname = os.path.basename(sys.argv[0])
	usage = """%prog [options] <path or db> ...
	
Various utilities related to BDB databases.

examples :
e2bdb.py -c   Is perhaps the most critical function, as it cleans up the database cache. See the Wiki for more.
e2bdb.py <path> -s    will list the contents of the database in a directory in bdb: notation
e2bdb.py <path> -l    Will give useful summary info about stacks in a directory
e2bdb.py <database> --dump    Gives a mechanism to dump all of the metadata in a database, even if the database contains no images
"""

	parser = OptionParser(usage=usage,version=EMANVERSION)

	parser.add_option("--cleanup","-c",action="store_true",default=False,help="This option will clean up the database cache so files can safely be moved or accessed on another computer via NFS.")
	parser.add_option("--force","-F",action="store_true",default=False,help="This will force an action that would normally fail due to failed checks.")
	parser.add_option("--delete",action="store_true",default=False,help="This will delete (or at least empty) the named database(s)")
	parser.add_option("--all","-a",action="store_true",help="List per-particle info",default=False)
	parser.add_option("--long","-l",action="store_true",help="Long listing",default=False)
	parser.add_option("--short","-s",action="store_true",help="Dense listing of names only",default=False)
	parser.add_option("--filt",type="string",help="Only include dictionary names containing the specified string",default=None)
	parser.add_option("--filtexclude",type="string",help="Exclude dictionary names containing the specified string",default=None)
	parser.add_option("--match",type="string",help="Only include dictionaries matching the provided Python regular expression",default=None)
	parser.add_option("--exclude",type="string",help="The name of a database containing a list of exclusion keys",default=None)
	parser.add_option("--dump","-D",action="store_true",help="List contents of an entire database, eg 'e2bdb.py -D refine_01#register",default=False)
	parser.add_option("--check",action="store_true",help="Check for self-consistency and errors in the structure of specified databases",default=False)
	parser.add_option("--nocache",action="store_true",help="Don't use the database cache fof this operation",default=False)
	
	parser.add_option("--makevstack",type="string",help="Creates a 'virtual' BDB stack with its own metadata, but the binary data taken from the (filtered) list of stacks",default=None)
	parser.add_option("--appendvstack",type="string",help="Appends to/creates a 'virtual' BDB stack with its own metadata, but the binary data taken from the (filtered) list of stacks",default=None)
	parser.add_option("--verbose", "-v", dest="verbose", action="store", metavar="n", type="int", default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	parser.add_option("--list",type="string",help="Specify the name of a file with a list of images to use in creation of virtual stacks. Please see source for details.",default=None)
	parser.add_option("--restore",type="string",help="Write changes in the derived virtual stack back to the original stack",default=None)


	(options, args) = parser.parse_args()

	if options.nocache : BDB_CACHE_DISABLE=True

	if options.cleanup : 
		db_cleanup(options.force)
		sys.exit(0)

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
		
	
	for path in args:
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
		
		
		if options.makevstack!=None or options.appendvstack!=None :
			
			vspath=os.path.realpath(vstack.path)+"/"
			if options.verbose>2 : print "vspath: ",vspath
			for db in dbs:
				dct,keys=db_open_dict(path+db,with_keys=True)
				if dct==vstack : continue
				vals = keys
				if keys == None: vals = range(len(dct))
				if options.list !=None: vals=slist
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
			continue

		if options.restore :
			nima = EMUtil.get_image_count(options.restore)
			IB = db_open_dict(options.restore)
			data_old = None
			for i in xrange(nima):
				source = IB.get_header(i)
				data_source = source["data_source"]
				ID = source["data_n"]
				if( data_old != data_source):
					if( data_old != None):  DB.close()
					DB = db_open_dict(data_source)
					data_old = data_source
				target = DB.get_header( ID )
				try:
					source["data_source"] = target["data_source"]
					source["data_n"]      = target["data_n"]
				except:
					del source['data_source']
					del source['data_n']
					del source['data_path']
				
				DB.set_header(ID, source)
			DB.close()
			continue

			
		if options.dump :
			for db in dbs:
				print "##### ",db
				dct=db_open_dict(path+db)
				
				#### Dump
				keys=dct.keys()
				keys.sort()
				for k in keys:
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
				dct=db_open_dict(path+db)
								
				### Info on all particles
				if options.all :
					for i in range(len(dct)):
						im=dct[i]
						print "%d. %d x %d x %d\tA/pix=%1.2f\tM=%1.4f\tS=%1.4f\tSk=%1.4f"%(i,im["nx"],im["ny"],im["nz"],im["apix_x"],im["mean"],im["sigma"],im["skewness"]),
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
				dct=db_open_dict(path+db)
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
