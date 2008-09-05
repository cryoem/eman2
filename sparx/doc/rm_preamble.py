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

# This program attempts to remove all preamble tags from an input LaTeX file and outputs
# the result.  The output file name ends in .tex while the input file name is appended
# with a .old suffix.

from os import rename

#-----------------------------------------------
def process(x):
    try:
        istream=open(x,"r");
    except IOError:
        print "Cannot find file "+x
        return

    ostream=open(x+".tex","w");

    line=istream.readline();    #take out the preamble
    while((line.find("\\begin{document}")==-1)and(line!="")):
        ostream.write("%"+line);
        line=istream.readline();

    if (line==""):  #no \begin{doc} found, so reset files
        try:            
            istream.close();
            ostream.close();
            istream=open(x,"r");
            ostream=open(x+".tex","w");
            line=istream.readline();
        except IOError:
            pass

    remove=["\\title", "\maketitle", "\documentclass", "\\begin{document}",
            "\end{document}", "\\tableofcontents"]; 
        
    while(line!=""):   #remove tags we don't want
        ok=1;
        for word in remove:
            if (line.find(word)!=-1):
                ok=0;
        if (ok==1):
            ostream.write(line);
        else:
            ostream.write("%"+line);
        line=istream.readline();

    try:
        istream.close();
        ostream.close();
        rename(x,x+".old")
        rename(x+".tex",x)
        print "Input File '%s' renamed to : %s.old\t\tOutput File: %s" % (x,x,x)
#        print "Input File "+x+"\t\tOutput file: "+x+".tex"
    except IOError:
        print "Had trouble closing the file "+x+" or "+x+".tex"            

#------------------------------------------------


import sys

if (len(sys.argv)<2):
    sys.stderr.write("Need a source LaTeX file");
    exit(1);

sys.argv.pop(0);
for x in sys.argv:
    process(x)
