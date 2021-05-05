#!/usr/bin/env python
#
# Author: Steven Ludtke, 11/25/2011 (sludtke@bcm.edu)
# Copyright (c) 2011- Baylor College of Medicine
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

from future import standard_library
standard_library.install_aliases()
from builtins import range
from builtins import object
from EMAN2 import *
from EMAN2db import db_open_dict, db_list_dicts
from math import *
import time
import os
import sys
import re
from pickle import dumps,loads,dump,load
from zlib import compress,decompress
from subprocess import Popen,PIPE
import traceback
import atexit
import random
import string

def main():
	global debug
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options]

	This program supports remote access of images and metadata from other machines. It is launched as a daemon
	and maintains open connections to databases and other computers for transfers.
	"""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	#parser.add_argument("--cleanup",action="store_true",default=False,help="This option will clean up the database cache so files can safely be moved or accessed on another computer via NFS.")
	parser.add_argument("--kill",action="store_true",default=False,help="This will kill an existing server instance before any other operations.")
	parser.add_argument("--ssh",type=str,help="Connect to a remote ssh server",default=None)
	parser.add_argument("--emen2",type=str,help="Connect to a remote EMEN2 server",default=None)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higher number means higher level of verboseness")

	(options, args) = parser.parse_args()

	run_daemon(options,args)
	
def run_daemon(options,args):
	import sys, os 

	# UNIX double-fork magic
	try: 
		pid = os.fork() 
		if pid > 0: sys.exit(0) # exit parent once
	except OSError as e: 
		print("Daemon error") 
		sys.exit(1)

	try: 
		os.setsid() 
		os.umask(0) 	# do we really want this
		os.chdir("/") 	# or this ?
	except:
		pass

	# do second fork
	try: 
		pid = os.fork() 
		if pid > 0: sys.exit(0) 	# exit the second time
	except OSError as e: 
		print("Daemon error 2") 
		sys.exit(1) 

	# ok, we got here, so we should be running in a parentless daemon now
    
class daemon(object):
	
	def __init__(self,options,args):
		
		self.listen_sock=None
		
		# Find a good random port and bind to it
		while self.listen_sock==None:
			self.listen_port=random.randint(10000,15000)
			self.listen_sock=socket()
			try: self.listen_sock.bind(("",self.listen_port))
			except: continue
			
		# This is a magic string for basic security
		self.magic="".join([random.choice(string.letters) for i in range(20)])
		
		# This file should be readable by the user only, and contains
		# "magic" string for security, port number, and PID
		out=open(e2gethome()+"/.eman2/remoted.txt","w")
		out.write("%s\n%d\n%d\n"%(self.magic,self.listen_port,os.getpid()))
		out.close()
    
    
def client_error(stdout,msg):
	"Return an error to the server in the form of !!!!message\n"
	msg=msg.replace("\n","")
	stdout.write("!!!!%s\n"%msg)

def write_obj(stdout,obj):
	"Writes an object to stdout as a length\\n followed by a compressed pickled string"
	xmit=compress(dumps(obj,-1),3)
	stdout.write("%d\n"%len(xmit))
	stdout.write(xmit)

def read_obj(stdin):
	"Reads a length followed by a compressed pickled string and returns the encapsulated object"
	size=stdin.readline()
	try: size=int(size)
	except:
		if size[:4]=="!!!!" : raise Exception(size[4:])
		raise Exception("Unknown error : "+size)
	return loads(decompress(stdin.read(size)))

def write_chunk(stdout,obj):
	"Writes an object to stdout as a length\\n followed by a compressed pickled string"
	if obj==None or len(obj)==0 : 
		stdout.write("0\n")
		return
	xmit=compress(obj,3)
	stdout.write("%d\n"%len(xmit))
	stdout.write(xmit)

def read_chunk(stdin):
	"Reads a length followed by a compressed string"
	size=stdin.readline()
	try: size=int(size)
	except:
		if size[:4]=="!!!!" : raise Exception(size[4:])
		raise Exception("Unknown error : "+size)
	if size==0 : return ""
	return decompress(stdin.read(size))

def send_file(stdout,path):
	"Sends a file to stdout as a set of chunks terminated with a 0 length chunk"
	fin=open(path,"rb")
	while 1:
		data=fin.read(1000000)
		if len(data)==0 :break
		write_chunk(stdout,data)
		if len(data)<1000000 : break
	write_chunk(stdout,None)

def recv_file(stdin,path):
	"Receives a file into path. Reads a set of chunks terminated with a zero-length chunk"
	try :os.makedirs(os.path.dirname(path))
	except: pass
	out=open(path,"w")
	while 1:
		chunk=read_chunk(stdin)
		if len(chunk)==0 or chunk==None : break
		out.write(chunk)

def send_bdb(stdout,path):
	"Sends a BDB to stdout as a set of compressed pickled key/value pairs, terminated by a None key"
	db=db_open_dict(path)
	keys=list(db.keys())
	for k in keys:
		write_obj(stdout,k)
		write_obj(stdout,db[k])
	write_obj(stdout,None)
	db.close()

def recv_bdb(stdin,path):
	"Receives a BDB from stdin as a set of None terminated compressed pickled key/value pairs"
	try :os.makedirs(path[4:].split("#")[0])
	except: pass
#	sys.stderr.write("open %s\n"%path)
	db=db_open_dict(path)
	db.realopen()
	db.bdb.truncate()			# empty the existing database
	while (1):
		k=read_obj(stdin)
		if k==None or (isinstance(k,str) and len(k)==0) : break
		db[k]=read_obj(stdin)
	db.close()

def get_bdb_list(path):
	"Given a bdb:* path, returns a list of items in the path"
	
	dicts=db_list_dicts(path)
	if '#' in path :
		if path.split("#")[-1] in dicts :		# check for existance of single specified dictionary
			return [path]
		else :
			return []
			
	return ["%s#%s"%(path,d) for d in dicts]

def get_dir_list_recurse(path,basepath=None):
	"Recursively lists the contents of a directory, including BDB contents"
	
	if ("EMAN2DB") in path:
		print("ERROR : EMAN2DB may not be specified as a path to copy. Use bdb: specifier instead.")
		return []
	
	if path[:4].lower()=="bdb:" :
		if basepath!=None : return get_bdb_list("bdb:"+basepath+"/"+path[4:])
		return get_bdb_list(path)
		
	if basepath!=None : path=basepath+"/"+path
	
	try :
		flist=os.listdir(path)
	except: return []
	
	rlist=[]
	for f in flist:
		fpath="%s/%s"%(path,f)
		if f=="EMAN2DB" : rlist.extend(get_bdb_list("bdb:%s"%path))
		elif os.path.isdir(fpath) : rlist.extend(get_dir_list_recurse(fpath))
		elif os.path.isfile(fpath) : rlist.append(fpath)
		
	if basepath!=None : return [i.replace(basepath+"/","") for i in rlist]
	return rlist
	
	
def scp_client():
	"""This acts as a client for interacting with pyscp.py over stdin/out"""
	
	stdout=sys.stdout
	stdin=sys.stdin
	
	stdout.write("HELO\n")
	
	while 1:
		stdout.flush()
		try : com=stdin.readline().strip()
		except : break
		if com=="exit" :
			break
		# mkdir command. Specify path to create. Returns 'OK' on success
		if com=="mkdir" :
			path=stdin.readline().strip()
			if path[:4].lower()=="bdb:" : path=path[4:].split("#")[0]
			if not os.path.exists(path) : 
				try : os.makedirs(path)
				except:
					client_error(stdout,'Failed to makedirs %s'%path)
					continue
			stdout.write("OK\n")
			continue
		
		# List a path. Returns #\npath\npath... where # is the number of returned lines
		if com=="listrecurse" :
			path=stdin.readline().strip()
			basepath=stdin.readline().strip()
			if len(basepath)==0 : basepath=None
			if path[:4].lower()=="bdb:" :
				plist=get_bdb_list(path)
			else :
				plist=get_dir_list_recurse(path,basepath)
			
			stdout.write("%d\n"%len(plist))
			for p in plist : stdout.write(p+"\n")
			continue
		
		# Return image header. Specify path/n#\n, return client_write return dictionary
		if com=="getheader":
			path=stdin.readline().strip()
			n=int(stdin.readline().strip())
			try: hdr=EMData(path,n,True)
			except: client_error(stdout,"Could not read header %s (%d)"%(path,n))
			write_obj(stdout,hdr.get_attr_dict())
			continue
		
		# Return an entire image as EMData object. Specify path\n#\n, client_write return EMData
		if com=="getimage" :
			path=stdin.readline().strip()
			n=int(stdin.readline().strip())
			try: img=EMData(path,n)
			except: client_error(stdout,"Could not read header %s (%d)"%(path,n))
			write_obj(stdout,img)
			continue
		
		# Writes an image on the client
		if com=="putimage" :
			path=stdin.readline().strip()
			n=int(stdin.readline().strip())
			img=read_obj(stdin)
			img.write_image(path,n)
			continue
		
		# Return a binary file object in chunks with a final 0 chunk
		if com=="getfile" :
			path=stdin.readline().strip()
			send_file(stdout,path)
			continue
		
		# Writes a binary file object, path\nseries of chunks ending with empty chunk
		if com=="putfile" :
			path=stdin.readline().strip()
			recv_file(stdin,path)
			continue
		
		# return a bdb dict as a set of n pickled key/value pairs
		if com=="getbdb" :
			path=stdin.readline().strip()
			send_bdb(stdout,path)
			continue
		
		# puts a bdb dict as a set of None terminated key/value pairs
		if com=="putbdb" :
			path=stdin.readline().strip()
			recv_bdb(stdin,path)
			continue
		
class scp_proxy(object):
	def __init__(self,userhost,verbose=0):
		"""Opens a connection to the remote host and establishes the remote client. userhost should be of the form "user@host"""
		self.verbose=verbose
		# Open the remote SSH process
		try :
			self.ssh=Popen(("ssh",userhost,"e2scp.py --client"),stdin=PIPE,stdout=PIPE)
			self.stdout=self.ssh.stdout		# read from this
			self.stdin=self.ssh.stdin		# write to this
		except:
			print("ssh to remote machine failed : ",("ssh",host,"e2ssh.py --client"))
			traceback.print_exc()
			sys.exit(2)
		
		while 1:
			ln=self.stdout.readline().strip()
			if len(ln)==0 : 
				print("Error running e2scp.py on the remote machine. EMAN2 installed ?")
				sys.exit(3)
			if ln=="HELO" : 
				if self.verbose : print("Connection established")
				break
			if self.verbose >1 : print("*** ",ln)
		
		atexit.register(self.close)
		
	def close(self):
		"""Close the connection"""
		if self.ssh!=None :
			try:
				self.stdin.write("exit\n")
#			self.stdin.flush()
				self.ssh.kill()
			except: pass
			self.ssh=None

	def mkdir(self,path):
		"""Create a path on the remote host, at all necessary levels"""
		
		self.stdin.write("mkdir\n%s\n"%path)
		self.stdin.flush()
		r=self.stdout.readline().strip()
		if r!="OK" : raise Exception("Error in creating remote path (%s)"%(r))

	def listrecurse(self,path,basepath=""):
		"""Recursively list the contents of a remote path, may be a directory or a BDB specifier. If specified
		will reference paths with respect to basepath."""
		self.stdin.write("listrecurse\n%s\n%s\n"%(path,basepath))
		r=int(self.stdout.readline().strip())
		ret=[]
		for i in range(r):
			ret.append(self.stdout.readline().strip())
			
		return ret
	
	def getheader(self,path,n):
		"""Return the header of a remote image as a dictionary"""
		self.stdin.write("getheader\n%s\n%d\n"%(path,n))
		self.stdin.flush()
		return read_obj(self.stdout,img)
		
	def getimage(self,path,n):
		"""Return a single reomote EMData object from a file"""
		self.stdin.write("getimage\n%s\n%d\n"%(path,n))
		self.stdin.flush()
		return read_obj(self.stdout,img)
		
	def putimage(self,path,n,img):
		"""Write a single EMData object into a remote file""" 
		self.stdin.write("putimage\n%s\n%d\n"%(path,n))
		self.write_obj(self.stdin,img)
		self.stdin.flush()
		
	def getfile(self,remotepath,localpath):
		"""Retrieve a single binary remote file and write to path. Streams for minimal memory usage."""
		self.stdin.write("getfile\n%s\n"%(remotepath))
		self.stdin.flush()
		recv_file(self.stdout,localpath)
		
	def putfile(self,localpath,remotepath):
		self.stdin.write("putfile\n%s\n"%(remotepath))
		self.stdin.flush()
		send_file(self.stdin,localpath)
		
	def getbdb(self,remotepath,localpath):
		"""Copyies a remote BDB to a local machine"""
		self.stdin.write("getbdb\n%s\n"%(remotepath))
		self.stdin.flush()
		recv_bdb(self.stdout,localpath)
		
	def putbdb(self,localpath,remotepath):
		"""Copies a local BDB to a remote machine"""
		self.stdin.write("putbdb\n%s\n"%(remotepath))
		self.stdin.flush()
		send_bdb(self.stdin,localpath)


		
if __name__ == "__main__":
	main()
