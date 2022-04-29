#!/usr/bin/env python
#
# Author: Steven Ludtke, 04/10/2003 (sludtke@bcm.edu)
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

# Purpose: to generate Processor manuals automatically
#
# Usage:
#    1. ./docgen.py
#    2. copy processors.html, processor_groups.html to ../doc
#

from EMAN2 import *
import sys
import time
import os
import shutil
import glob

    
def write_header(output, name):
    output.write(b"<head><title>EMAN2 " + name.encode() + b" Manual</title></head><body>\n")
    output.write(b"<h1> <center> <font color=\"blue\">EMAN2 " + name.encode() + b" Manual </font></center></h1>\n")
    output.write(b"<br>")
    output.write(b"<i>Last modified on " + time.strftime('%a, %d %b %Y %H:%M:%S %Z').encode())
    output.write(b"<br>")
    output.write(b"<i>This document is automatically generated. Please don't edit it.</i>\n")

    output.write(b"<br><br>")
    output.write(b"<table border=1 cellspacing=4 cellpadding=4>")
    output.write(b"<tr>\n")
    output.write(b"  <td align=center> <font size='+1'><b>" + name.encode() + b" Name </b></font></td>\n")
    output.write(b"  <td align=center> <font size='+1'><b>Parameters </b></font></td>\n")
    output.write(b"  <td align=center> <font size='+1'><b>Description </b></font> </td>\n")
    output.write(b"</tr>\n")


def write_tail(output):
    output.write(b"</table>\n</body>\n")


def write_processor(processorname, output, bgcolor):
    output.write(b"<tr bgcolor=" + bgcolor.encode() + b">\n")
    output.write(b"  <td> <font color='0000a0'> <b>" + processorname.encode() + b" </b></font> </td>\n")

    processor = Processors.get(processorname)

    typedict = processor.get_param_types()
    output.write(b"  <td>")
    typekeys = list(typedict.keys())

    for typekey in typekeys:
        output.write(typedict.get_type(typekey).lower().encode() + b" ")
        output.write(b"<font color=green> <b>" + typekey.encode() + b"</b></font>")
        output.write(b": " + typedict.get_desc(typekey).encode())
        output.write(b"<br>")

    output.write(b"</td>\n")
    output.write(b"  <td>" + processor.get_desc().encode() + b"</td>\n")
    output.write(b"</tr>\n")


def write_cmp(cmp_name, output, bgcolor):
    output.write(b"<tr bgcolor=" + bgcolor.encode() + b">\n")
    output.write(b"  <td> <font color='0000a0'> <b>" + cmp_name.encode() + b" </b></font> </td>\n")
    
    cmp = Cmps.get(cmp_name)
    
    typedict = cmp.get_param_types()
    output.write(b"  <td>")
    typekeys = list(typedict.keys())
    
    for typekey in typekeys:
        output.write(typedict.get_type(typekey).lower().encode() + b" ")
        output.write(b"<font color=green> <b>" + typekey.encode() + b"</b></font>")
        output.write(b": " + typedict.get_desc(typekey).encode())
        output.write(b"<br>")
    
    output.write(b"</td>\n")
    output.write(b"  <td>" + cmp.get_desc().encode() + b"</td>\n")
    output.write(b"</tr>\n")


def write_aligner(aligner_name, output, bgcolor):
    output.write(b"<tr bgcolor=" + bgcolor.encode() + b">\n")
    output.write(b"  <td> <font color='0000a0'> <b>" + aligner_name.encode() + b" </b></font> </td>\n")
    
    align = Aligners.get(aligner_name)
    
    typedict = align.get_param_types()
    output.write(b"  <td>")
    typekeys = list(typedict.keys())
    
    for typekey in typekeys:
        output.write(typedict.get_type(typekey).lower().encode() + b" ")
        output.write(b"<font color=green> <b>" + typekey.encode() + b"</b></font>")
        output.write(b": " + typedict.get_desc(typekey).encode())
        output.write(b"<br>")
    
    output.write(b"</td>\n")
    output.write(b"  <td>" + align.get_desc().encode() + b"</td>\n")
    output.write(b"</tr>\n")


def write_projector(projector_name, output, bgcolor):
    output.write(b"<tr bgcolor=" + bgcolor.encode() + b">\n")
    output.write(b"  <td> <font color='0000a0'> <b>" + projector_name.encode() + b" </b></font> </td>\n")
    
    project = Projectors.get(projector_name)
    
    typedict = project.get_param_types()
    output.write(b"  <td>")
    typekeys = list(typedict.keys())
    
    for typekey in typekeys:
        output.write(typedict.get_type(typekey).lower().encode() + b" ")
        output.write(b"<font color=green> <b>" + typekey.encode() + b"</b></font>")
        output.write(b": " + typedict.get_desc(typekey).encode())
        output.write(b"<br>")
    
    output.write(b"</td>\n")
    output.write(b"  <td>" + project.get_desc().encode() + b"</td>\n")
    output.write(b"</tr>\n")


def write_reconstructor(reconstructor_name, output, bgcolor):
    output.write(b"<tr bgcolor=" + bgcolor.encode() + b">\n")
    output.write(b"  <td> <font color='0000a0'> <b>" + reconstructor_name.encode() + b" </b></font> </td>\n")
    
    reconstruct = Reconstructors.get(reconstructor_name)
    
    typedict = reconstruct.get_param_types()
    output.write(b"  <td>")
    typekeys = list(typedict.keys())
    
    for typekey in typekeys:
        output.write(typedict.get_type(typekey).lower().encode() + b" ")
        output.write(b"<font color=green> <b>" + typekey.encode() + b"</b></font>")
        output.write(b": " + typedict.get_desc(typekey).encode())
        output.write(b"<br>")

    output.write(b"</td>\n")
    output.write(b"  <td>" + reconstruct.get_desc().encode() + b"</td>\n")
    output.write(b"</tr>\n")


def write_averager(averager_name, output, bgcolor):
    output.write(b"<tr bgcolor=" + bgcolor.encode() + b">\n")
    output.write(b"  <td> <font color='0000a0'> <b>" + averager_name.encode() + b" </b></font> </td>\n")
    
    average = Averagers.get(averager_name)

    typedict = average.get_param_types()
    output.write(b"  <td>")
    typekeys = list(typedict.keys())
    
    for typekey in typekeys:
        output.write(typedict.get_type(typekey).lower().encode() + b" ")
        output.write(b"<font color=green> <b>" + typekey.encode() + b"</b></font>")
        output.write(b": " + typedict.get_desc(typekey).encode())
        output.write(b"<br>")
    
    output.write(b"</td>\n")
    output.write(b"  <td>" + average.get_desc().encode() + b"</td>\n")
    output.write(b"</tr>\n")


def write_analyzer(analyzer_name, output, bgcolor):
    output.write(b"<tr bgcolor=" + bgcolor.encode() + b">\n")
    output.write(b"  <td> <font color='0000a0'> <b>" + analyzer_name.encode() + b" </b></font> </td>\n")
    
    analyze = Analyzers.get(analyzer_name)
    
    typedict = analyze.get_param_types()
    output.write(b"  <td>")
    typekeys = list(typedict.keys())
    
    for typekey in typekeys:
        output.write(typedict.get_type(typekey).lower().encode() + b" ")
        output.write(b"<font color=green> <b>" + typekey.encode() + b"</b></font>")
        output.write(b": " + typedict.get_desc(typekey).encode())
        output.write(b"<br>")
    
    output.write(b"</td>\n")
    output.write(b"  <td>" + analyze.get_desc().encode() + b"</td>\n")
    output.write(b"</tr>\n")


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
    output.write(b"<tr bgcolor=" + groupbg.encode() + b">\n")
    output.write(b"  <td> <font color='a000a0' size=+1> <b>" + groupname.encode() + b" </b></font> </td>\n")
    output.write(b"  <td> </td>")
    output.write(b"  <td><b>")
    if groupclass:
        output.write(groupclass.get_group_desc().encode())
        
    output.write(b"</b></td>\n")
    output.write(b"</tr>\n")


def write_single_processors():
    out = open("processors.html", "wb")
    write_header(out, 'Processor')
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


def write_single_comparators():
    out = open("cmps.html", "wb")
    write_header(out, 'Cmp')
    cmp_names = Cmps.get_list()
    bgcolor1 = 'f0f0fa'
    bgcolor2 = 'fafaf0'
    bgcolor = bgcolor1
    
    for cmp_name in cmp_names:
        if bgcolor == bgcolor1:
            bgcolor = bgcolor2
        else:
            bgcolor = bgcolor1
        
        write_cmp(cmp_name, out, bgcolor)
    write_tail(out)


def write_single_aligners():
    out = open("aligners.html", "wb")
    write_header(out, 'Aligner')
    aligner_names = Aligners.get_list()
    bgcolor1 = 'f0f0fa'
    bgcolor2 = 'fafaf0'
    bgcolor = bgcolor1
    
    for aligner_name in aligner_names:
        if bgcolor == bgcolor1:
            bgcolor = bgcolor2
        else:
            bgcolor = bgcolor1
    
        write_aligner(aligner_name, out, bgcolor)
    write_tail(out)


def write_single_projectors():
    out = open("projectors.html", "wb")
    write_header(out, 'Projector')
    projector_names = Projectors.get_list()
    bgcolor1 = 'f0f0fa'
    bgcolor2 = 'fafaf0'
    bgcolor = bgcolor1
    
    for projector_name in projector_names:
        if bgcolor == bgcolor1:
            bgcolor = bgcolor2
        else:
            bgcolor = bgcolor1
    
        write_projector(projector_name, out, bgcolor)
    write_tail(out)


def write_single_reconstructors():
    out = open("reconstructors.html", "wb")
    write_header(out, 'Reconstructor')
    reconstructor_names = Reconstructors.get_list()
    bgcolor1 = 'f0f0fa'
    bgcolor2 = 'fafaf0'
    bgcolor = bgcolor1
    
    for reconstructor_name in reconstructor_names:
        if bgcolor == bgcolor1:
            bgcolor = bgcolor2
        else:
            bgcolor = bgcolor1
    
        write_reconstructor(reconstructor_name, out, bgcolor)
    write_tail(out)


def write_single_averagers():
    out = open("averagers.html", "wb")
    write_header(out, 'Averager')
    averager_names = Averagers.get_list()
    bgcolor1 = 'f0f0fa'
    bgcolor2 = 'fafaf0'
    bgcolor = bgcolor1
    
    for averager_name in averager_names:
        if bgcolor == bgcolor1:
            bgcolor = bgcolor2
        else:
            bgcolor = bgcolor1
    
        write_averager(averager_name, out, bgcolor)
    write_tail(out)


def write_single_analyzers():
    out = open("analyzers.html", "wb")
    write_header(out, 'Analyzer')
    analyzer_names = Analyzers.get_list()
    bgcolor1 = 'f0f0fa'
    bgcolor2 = 'fafaf0'
    bgcolor = bgcolor1
    
    for analyzer_name in analyzer_names:
        if bgcolor == bgcolor1:
            bgcolor = bgcolor2
        else:
            bgcolor = bgcolor1
    
        write_analyzer(analyzer_name, out, bgcolor)
    write_tail(out)


def write_group_processors():
    gout = open("processor_groups.html", "wb")
    write_header(gout, 'Processor Group')
    processorgroups = group_processors()
    groupnames = list(processorgroups.keys())
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
    write_single_comparators()
    write_single_aligners()
    write_single_projectors()
    write_single_reconstructors()
    write_single_averagers()
    write_single_analyzers()


    outdir = 'modular_class_html'
    os.makedirs(outdir, exist_ok=True)

    for f in glob.glob('*.html'):
        fnew = os.path.join(outdir, f)
        print(f"Moving: {f} -> {fnew}")
        shutil.move(f, fnew)


if __name__ == '__main__':
    main()
