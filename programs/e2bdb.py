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

from EMAN2 import *
from optparse import OptionParser
from math import *
import time
import os
import sys

def main():
	global debug
	progname = os.path.basename(sys.argv[0])
	usage = """%prog [options] <path or db> ...
	
Various utilities related to BDB databases."""

	parser = OptionParser(usage=usage,version=EMANVERSION)

	parser.add_option("--long","-l",action="store_true",help="Long listing",default=False)
	parser.add_option("--short","-s",action="store_true",help="Dense listing of names only",default=False)
	parser.add_option("--filt",type="string",help="Only include dictionaries containing the specified string (since wildcards like * cannot be used)",default=None)

	(options, args) = parser.parse_args()

	if len(args)==0 : args.append("bdb:.")
		
	for path in args:
		if path.lower()[:4]!="bdb:" : path="bdb:"+path
		if not '#' in path and path[-1]!='/' : path+='/'
		if len(args)>1 : print path,":"
		
		dbs=db_list_dicts(path)
		
		dbs.sort()
		if options.filt:
			dbs=[db for db in dbs if options.filt in db]
		
		try: maxname=max([len(s) for s in dbs])
		except: 
			print "Error reading ",path
			continue
			
		# long listing, one db per line
		if options.long :
			width=maxname+3
			fmt="%%-%ds %%-07d %%dx%%dx%%d  %%s"%width
			fmt2="%%-%ds (not an image stack)"%width
			for db in dbs:
				dct=db_open_dict(path+"#"+db)
				first=EMData()
				try: 
					first.read_image(path+"#"+db,0,True)
					size=first.get_xsize()*first.get_ysize()*first.get_zsize()*len(dct)*4;
					if size>1000000000: size="%1.2f gb"%(size/1000000000)
					elif size>1000000: size="%1.2f mb"%(size/1000000)
					else: size="%1.2f kb"%(size/1000)
					print fmt%(db,len(dct),first.get_xsize(),first.get_ysize(),first.get_zsize(),size)
				except:
					print fmt2%db
		elif options.short :
			for db in dbs:
				print "bdb:"+db,
			print " "

		else :
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
			
			
if __name__ == "__main__":
	main()
