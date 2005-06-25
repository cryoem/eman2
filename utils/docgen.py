#!/usr/bin/env python

# Purpose: to generate Processor manuals automatically
#
# Usage:
#    1. ./docgen.py
#    2. copy processors.html, processor_groups.html to ../doc
#

from EMAN2 import *
import sys
import time

	
def write_header(output):
	output.write("<head><title>EMAN2 Processor Manual</title></head><body>\n")
	output.write("<h1> <center> <font color=\"blue\">EMAN2 Processor Manual </font></center></h1>\n")
	output.write("<br>")
	output.write("<i>Last modified on " + time.strftime('%a, %d %b %Y %H:%M:%S %Z'))
	output.write("<br>")
	output.write("<i>This document is automatically generated. Please don't edit it.</i>\n")

	output.write("<br><br>")
	output.write("<table border=1 cellspacing=4 cellpadding=4>")
	output.write("<tr>\n")
	output.write("  <td align=center> <font size='+1'><b>Processor Name </b></font></td>\n")
	output.write("  <td align=center> <font size='+1'><b>Parameters </b></font></td>\n")
	output.write("  <td align=center> <font size='+1'><b>Description </b></font> </td>\n")
	output.write("</tr>\n")

def write_tail(output):
	output.write("</table>\n</body>\n")

def write_processor(processorname, output, bgcolor):
	output.write("<tr bgcolor=" + bgcolor + ">\n")
	output.write("  <td> <font color='0000a0'> <b>"  + processorname + " </b></font> </td>\n")

	processor = Processors.get(processorname)

	typedict = processor.get_param_types()
	output.write("  <td>")
	typekeys = typedict.keys()

	for typekey in typekeys:
		output.write(typedict.get_type(typekey).lower() + " ")
		output.write("<font color=green> <b>" + typekey + "</b></font>")
		output.write(": " + typedict.get_desc(typekey))
		output.write("<br>")

	output.write("</td>\n")
	output.write("  <td>" + processor.get_desc() + "</td>\n")
	output.write("</tr>\n")

def write_group(groupname, output):
	groupclass = None
	processor_names = Processors.get_list()
	
	if groupname == "Processor":
		groupclass = Processor
	if groupname == "RealPixelFilter":
		groupclass = RealPixelFilter
	elif groupname == "BoxStatFilter":
		groupclass = BoxStatFilter
	elif groupname == "ComplexPixelFilter":
		groupclass = ComplexPixelFilter
	elif groupname == "CoordinateFilter":
		groupclass = CoordinateFilter
	elif groupname == "FourierFilter":
		groupclass = FourierFilter
	elif groupname == "NormalizeFilter":
		groupclass = NormalizeFilter

	groupbg = "a0a0f0"
	output.write("<tr bgcolor=" + groupbg + ">\n")
	output.write("  <td> <font color='a000a0' size=+1> <b>"  + groupname + " </b></font> </td>\n")
	output.write("  <td> </td>")
	output.write("  <td><b>");
	if groupclass:
		output.write(groupclass.get_group_desc())
		
	output.write("</b></td>\n")
	output.write("</tr>\n")
		


def write_single_processors():
	out = open("processors.html", "wb")
	write_header(out)
	processor_names = Processors.get_list()
	bgcolor1 = 'f0f0fa'
	bgcolor2 = 'fafaf0'
	bgcolor = bgcolor1
	
	for processorname in processor_names:
		if bgcolor == bgcolor1:
			bgcolor = bgcolor2
		else:
			bgcolor = bgcolor1

		write_processor(processorname, out, bgcolor)
	write_tail(out)
	
	

def write_group_processors():
	gout = open("processor_groups.html", "wb")
	write_header(gout)
	processorgroups = group_processors()
	groupnames = processorgroups.keys()
	bgcolor1 = 'f0f0fa'
	bgcolor2 = 'fafaf0'
	bgcolor = bgcolor1

	write_group("Processor", gout)

	sorted_groupnames = []
	for groupname in groupnames:
		if groupname != "Others":
			sorted_groupnames.append(groupname)

	sorted_groupnames.sort()
	sorted_groupnames.append("Others")

	for groupname in sorted_groupnames:
		groupitems = processorgroups[groupname]
		write_group(groupname, gout)
		
		for processorname in groupitems:
			if bgcolor == bgcolor1:
				bgcolor = bgcolor2
			else:
				bgcolor = bgcolor1

			write_processor(processorname, gout, bgcolor)
			
		
def main():	
	write_single_processors()
	write_group_processors()
	
if __name__ == '__main__':
    main()
