#!/bin/env/python

# extracthelp.py -  Steven Ludtke  09/14/2010

# This program will extract the command-line help options from another EMAN2 program and output
# them in a wiki-friendly format

from sys import argv
from optparse import OptionParser

parser = OptionParser(usage="doesn't matter",version="1.0")

lines=[i for i in file(sys.argv[1],"r") if "parser.add_option" in i]

for i in lines: 

