#!/usr/bin/env python

from EMAN2 import *
import sys
import time


bgcolor1 = 'f0f0fa'
bgcolor2 = 'fafaf0'
bgcolor = bgcolor1
	
def write_header(output):
	output.write("<head><title>EMAN2 Filter Manual</title></head><body>\n")
	output.write("<h1> <center> <font color=\"blue\">EMAN2 Filter Manual </font></center></h1>\n")
	output.write("<br>")
	output.write("<i>Last modified on " + time.strftime('%a, %d %b %Y %H:%M:%S %Z'))
	output.write("<br>")
	output.write("<i>This document is automatically generated. Please don't edit it.</i>\n")

	output.write("<br><br>")
	output.write("<table border=1 cellspacing=4 cellpadding=4>")
	output.write("<tr>\n")
	output.write("  <td align=center> <font size='+1'><b>Filter Name </b></font></td>\n")
	output.write("  <td align=center> <font size='+1'><b>Parameters </b></font></td>\n")
	output.write("  <td align=center> <font size='+1'><b>Description </b></font> </td>\n")
	output.write("</tr>\n")

def write_tail(output):
	output.write("</table>\n</body>\n")

def write_filter(filtername, output):
	if bgcolor == bgcolor1:
			bgcolor = bgcolor2
		else:
			bgcolor = bgcolor1

	output.write("<tr bgcolor=" + bgcolor + ">\n")
	output.write("  <td> <font color='0000a0'> <b>"  + filtername + " </b></font> </td>\n")
	filter = Filters.get(filtername)

	typedict = filter.get_param_types()
	output.write("  <td>")
	typekeys = typedict.keys()

	for typekey in typekeys:
		output.write(typedict.get_type(typekey).lower() + " ")
		output.write("<font color=green> <b>" + typekey + "</b></font>")
		output.write(": " + typedict.get_desc(typekey))
		output.write("<br>")

	output.write("</td>\n")
	output.write("  <td>" + filter.get_desc() + "</td>\n")
	output.write("</tr>\n")

def write_group(groupname, output):
	groupclass = None

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

	


		


def write_single_filters():
	out = open("filters.html", "wb")
	write_header(out)
	filter_names = Filters.get_list()
	
	for filtername in filter_names:
		write_filter(filtername, out)
	write_tail(out)
	
	

def write_group_filters():
	gout = open("filter_groups.html", "wb")
	write_header(gout)

	filtergroups = group_filters()

	groupnames = filtergroups.keys()
	
	for groupname in groupnames:
		groupitems = filtergroups[groupname]

		if groupname == "RealPixelFilter":
			

		for filtername in groupitems:
			#write_filter(filtername, gout)
			print "  ", filtername

		
def main():
	write_single_filters()
	
if __name__ == '__main__':
    main()
