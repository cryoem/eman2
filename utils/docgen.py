#!/usr/bin/env python

from EMAN2 import *
import sys
import time

out = open("filters.html", "wb")
out.write("<head><title>EMAN2 Filter Manual</title></head><body>\n")
out.write("<h1> <center> <font color=\"blue\">EMAN2 Filter Manual </font></center></h1>\n")
out.write("<br>")
out.write("<i>Last modified on " + time.strftime('%a, %d %b %Y %H:%M:%S %Z'))
out.write("<br>")
out.write("<i>This document is automatically generated. Please don't edit it.</i>\n")


out.write("<br><br>")
out.write("<table border=1 cellspacing=4 cellpadding=4>")
out.write("<tr>\n")
out.write("  <td align=center> <font size='+1'><b>Filter Name </b></font></td>\n")
out.write("  <td align=center> <font size='+1'><b>Parameters </b></font></td>\n")
out.write("  <td align=center> <font size='+1'><b>Description </b></font> </td>\n")
out.write("</tr>\n")


filter_names = Filters.get_list()
bgcolor1 = 'f0f0fa'
bgcolor2 = 'fafaf0'
bgcolor = bgcolor1

for filtername in filter_names:
	if bgcolor == bgcolor1:
		bgcolor = bgcolor2
	else:
		bgcolor = bgcolor1
	
	out.write("<tr bgcolor=" + bgcolor + ">\n")
	out.write("  <td> <font color='0000a0'> <b>"  + filtername + " </b></font> </td>\n")
	filter = Filters.get(filtername)
	
	typedict = filter.get_param_types()
	out.write("  <td>")
	typekeys = typedict.keys()
	
	for typekey in typekeys:
		out.write(typedict.get_type(typekey).lower() + " ")
		out.write("<font color=green> <b>" + typekey + "</b></font>")
		out.write(": " + typedict.get_desc(typekey))
		out.write("<br>")

	out.write("</td>\n")
	out.write("  <td>" + filter.get_desc() + "</td>\n")
	out.write("</tr>\n")
	
out.write("</table>\n</body>\n")

	
