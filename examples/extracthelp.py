#!/usr/bin/env python

# extracthelp.py -  Steven Ludtke  09/14/2010

# This program will extract the command-line help options from another EMAN2 program and output
# them in a wiki-friendly format

import sys
from optparse import OptionParser

parser = OptionParser(usage="doesn't matter",version="1.0")

lines=[i for i in file(sys.argv[1],"r")]


for j,i in enumerate(lines): 
 	if "parser.add_option" in i:
		try : exec i.strip()
		except:
			try: exec (i+lines[j+1]).strip()
			except: print (i+lines[j+1]).strip()
	
#print parser.__dict__.keys()
for i in parser.option_list:
	try: o_long=i._long_opts[0]
	except: o_long=" "
	
	try: o_short=i._short_opts[0]
	except: o_short=" "
	
	try: 
		o_help=i.help
		if o_help==None : raise Exception
	except: o_help=" "
	
	try: 
		o_type=i.type
		if o_type==None : raise Exception
	except: o_type="bool"
	
	print "||%s||%s||%s||%s||"%(o_short,o_long,o_type,o_help)
