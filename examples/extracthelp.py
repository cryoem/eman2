#!/usr/bin/env python

# extracthelp.py -  Steven Ludtke  09/14/2010
# redone for EMAN2.1 5/15/14

# This program will extract the command-line help options from another EMAN2 program and output
# them in a wiki-friendly format

import sys
#from optparse import OptionParser
#from EMAN2 import *

options=[]
def parse(s):
	global options
	
	sub="myfn"+s[s.find("("):]
	return eval(sub)
	

def myfn(*args,**kargs):
	return (args,kargs)

lines=[i for i in file(sys.argv[1],"r")]

for j,i in enumerate(lines): 
 	if "parser.add_argument" in i:
		try : op=parse(i.strip())
		except:
			try: op=parse((i+lines[j+1]).strip())
			except: print (i+lines[j+1]).strip()
			
		com={"name":op[0][0][2:]}
		com.update(op[1])
		if len(op[0])>1 : com["short"]=op[0][1]

		options.append(com)

#print options

#print parser.__dict__.keys()
#print parser.optionslist
#print "0000000000000000000000000"
typelookup={ type(""):"string",type(1):"int",type(1.0):"float","bool":"bool" }

for i in options:
	try: o_long=i["name"]
	except: o_long=" "
	
	try: o_short=i["short"]
	except: o_short=" "
	
	try: 
		o_help=i["help"]
		if o_help==None : raise Exception
	except: o_help=" "
	
	try: 
		o_type=typelookup[i["type"]]
		if o_type==None : raise Exception
	except: o_type="bool"
	
	print "||%s||%s||%s||%s||"%(o_short,o_long,o_type,o_help)
