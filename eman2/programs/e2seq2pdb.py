#!/usr/bin/env python

# Author: Matthew Baker
# Copyright (c) 2000-2012 Baylor College of Medicine

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


import sys
from EMAN2 import *
from optparse import OptionParser
import os
import commands
from sys import argv

#           1         2         3         4         5         6         7
# 01234567890123456789012345678901234567890123456789012345678901234567890123456789
# ATOM      1  N   ASP L   1      57.169  62.936   7.756  1.00 30.65      1ACY 129
# HEADER    COMPLEX(ANTIBODY/HIV-1 FRAGMENT)        10-FEB-94   1ACY      1ACY   2


def getThreeLetter(aa):
	if aa=="A": threeletter="ALA"
	elif aa=="C": threeletter="CYS"
	elif aa=="D": threeletter="ASP"
	elif aa=="E": threeletter="GLU"
	elif aa=="F": threeletter="PHE"
	elif aa=="G": threeletter="GLY"
	elif aa=="H": threeletter="HIS"
	elif aa=="I": threeletter="ILE"
	elif aa=="K": threeletter="LYS"
	elif aa=="L": threeletter="LEU"
	elif aa=="M": threeletter="MET"
	elif aa=="N": threeletter="ASN"
	elif aa=="P": threeletter="PRO"
	elif aa=="Q": threeletter="GLN"
	elif aa=="R": threeletter="ARG"
	elif aa=="S": threeletter="SER"
	elif aa=="T": threeletter="THR"
	elif aa=="V": threeletter="VAL"
	elif aa=="W": threeletter="TRP"
	elif aa=="Y": threeletter="TYR"
	else: threeletter="UNK"
	
	return threeletter
		

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """%prog [options] <PDBmodel> <sequence> <outputPDB> 
	
	Builds all possible chain traces
	"""
	
	parser = OptionParser(usage=usage,version=EMANVERSION)
	(options, args) = parser.parse_args()

	seqfile=open(args[1], "r")
	seqlines=seqfile.readlines()
	seqfile.close()	
	aaa=[]
	for sline in seqlines:
		for a in sline:
			aaa.append(getThreeLetter(a))

	##################################get PDB
	pdbfile=open(args[0], "r")
	lines=pdbfile.readlines()
	pdbfile.close()


	outfile=args[2]
	out=open(outfile,"w")
	i=0
	for line in lines:
		isatom=str(line[0:6].strip())
		if (isatom=="ATOM"):
			i=i+1
			out.write("%s%5d%s%s%s%4d%s"%(line[:6],i,line[11:17],aaa[i-1],line[20:22],i,line[26:]))
				

	revout=outfile.split('.')[0]+"-r.pdb"
	rout=open(revout,"w")			
	lines.reverse()
	i=0
	for line in lines:
		isatom=str(line[0:6].strip())
		if (isatom=="ATOM"):
			i=i+1
			rout.write("%s%5d%s%s%s%4d%s"%(line[:6],i,line[11:17],aaa[i-1],line[20:22],i,line[26:]))				
	out.close()				
	rout.close()
	

	
# If executed as a program
if __name__ == '__main__':
	main() 