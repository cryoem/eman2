#!/usr/bin/env python

#
# Author: Steven Ludtke, 3/1/2010 (sludtke@bcm.edu)
# Copyright (c) 2010- Baylor College of Medicine
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
from optparse import OptionParser
from math import *
import time
import os
import sys
import re
from cPickle import dumps,loads,dump,load
from zlib import compress,decompress
from subprocess import Popen,PIPE
import traceback

def main():
	global debug
	progname = os.path.basename(sys.argv[0])
	usage = """%prog [options] <path or db> ... <target>
	
This program is used to copy directories or files including BDB databases between machines. Requires a
properly configured SSH client on the local machine, and SSH server running on target machine. 

- Sources may be of the form: path, bdb:path, user@host:path, user@host:bdb:path
- Target may not include bdb: specifiers, but must be a directory (local or remote)
- Only sources or target may include user@host specification, not both.
- user is not optional in remote specification
"""

	parser = OptionParser(usage=usage,version=EMANVERSION)

	#parser.add_option("--cleanup","-c",action="store_true",default=False,help="This option will clean up the database cache so files can safely be moved or accessed on another computer via NFS.")
	#parser.add_option("--filt",type="string",help="Only include dictionary names containing the specified string",default=None)
	parser.add_option("--client","-c",action="store_true",default=False,help="This option is for internal use only. Do NOT specify.")
	parser.add_option("--verbose", "-v", dest="verbose", action="store", metavar="n", type="int", default=0, help="verbose level [0-9], higner number means higher level of verboseness")

	(options, args) = parser.parse_args()

	if options.client : 
		scp_client()
		sys.exit(0)

	# Make sure at most one of source or dest is remote
	if '@' in args[-1] and '@' in "".join(args[:-1]) :
		print "Remote specification may not be in both source and target"
		sys.exit(1)
	
	# target is remote
	if '@' in args[-1] :
		# Open the remote SSH process
		host=args[-1].split(":")[0]
		try :
			ssh=Popen(("ssh",host,"e2ssh.py --client"),stdin=PIPE,stdout=PIPE)
		except:
			print "ssh to remote machine failed : ",("ssh",host,"e2ssh.py --client")
			traceback.print_exc()
			sys.exit(2)
		
		while 1:
			ln=ssh.stdout.readline().strip()
			if len(ln)==0 : continue
			if ln=="HELO" : 
				if options.verbose : print "Connection established"
				break
			if options.verbose >1 : print "*** ",ln
		
		
		
		

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
		if size[:4]=="!!!!" : raise Exception,size[4:]
		raise Exception,"Unknown error : "+size
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
		if size[:4]=="!!!!" : raise Exception,size[4:]
		raise Exception,"Unknown error : "+size
	if size==0 : return ""
	return decompress(stdin.read(size))

def send_file(stdout,path):
	"Sends a file to stdout as a set of chunks terminated with a 0 length chunk"
	fin=file(path,"rb")
	while 1:
		data=fin.read(1000000)
		if len(data)==0 :break
		write_chunk(stdout,data)
		if len(data)<1000000 : break

def recv_file(stdin,path):
	"Receives a file into path. Reads a set of chunks terminated with a zero-length chunk"
	out=file(path,"w")
	while 1:
		chunk=read_chunk(stdin)
		if len(chunk)==0 or chunk==None : break
		out.write(chunk)

def send_bdb(stdout,path):
	"Sends a BDB to stdout as a set of compressed pickled key/value pairs, terminated by a None key"
	db=db_open_dict(path)
	keys=db.keys()
	for k in keys:
		write_obj(stdout,k)
		write_obj(stdout,db[k])
	write_obj(stdout,None)
	db.close()

def recv_bdb(stdin,path):
	"Receives a BDB from stdin as a set of None terminated compressed pickled key/value pairs"
	db=db_open_dict(path)
	db.bdb.truncate()			# empty the existing database
	while (1):
		k=read_obj(stdin)
		if k==None or len(k)==0 : break
		db[k]=read_obj(stdin)
	db.close()

def get_bdb_list(ath):
	"Given a bdb:* path, returns a list of items in the path"
	
	dicts=db_list_dicts(path)
	if '#' in path :
		if path.split("#")[-1] in dicts :		# check for existance of single specified dictionary
			return [path]
		else :
			return []
			
	return ["%s#%s"%(path,d) for d in dicts]

def get_dir_list_recurse(path):
	"Recursively lists the contents of a directory, including BDB contents"
	
	try : flist=os.listdir(path)
	except: return []
	
	rlist=[]
	for f in flist:
		fpath="%s/%s"%(path,f)
		if f=="EMAN2DB" : rlist.extend(get_bdb_list("bdb:%s"%path))
		elif os.path.isdir(fpath) : rlist.extend(get_dir_list_recurse(fpath))
		elif os.path.isfile(fpath) : rlist.append(fpath)
		
	return rlist
	
	
def scp_client():
	"""This acts as a client for interacting with pyscp.py over stdin/out"""
	
	stdout=sys.stdout
	stdin=sys.stdin
	
	stdout.write("HELO/n")
	
	while 1:
		stdout.flush()
		try : com=stdin.readline().strip()
		except : break
		if com=="bye" :
			break
		# mkdir command. Specify path to create. Returns 'OK' on success
		if com=="mkdir" :
			path=stdin.readline().strip()
			if path[:4].lower()=="bdb:" : path=path[4:].split("#")[0]
			try : os.makedirs(path)
			except:
				client_error('Failed to makedirs %s'%path)
				continue
			stdout.write("OK\n")
			continue
		
		# List a path. Returns #\npath\npath... where # is the number of returned lines
		if com=="listrecurse" :
			path=stdin.readline().strip()
			if path[:4].lower()=="bdb:" :
				plist=get_bdb_list(path)
			else :
				plist=get_dir_list_recurse(path)
			
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
			path=stdin.readline.strip()
			recv_bdb(stdin,path)
			continue
		
		
		
if __name__ == "__main__":
	main()
