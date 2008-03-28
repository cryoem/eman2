#!/usr/bin/python

from EMAN2 import *
from optparse import OptionParser
import os
import sys

#constants
HEADER_ONLY=True
HEADER_AND_DATA=False


def main():
	progname = os.path.basename(sys.argv[0])
	usage = """%prog [options] <input file> <tlt file> <output file>"""

	parser = OptionParser(usage=usage,version=EMANVERSION)
	(options, args) = parser.parse_args()
	if len(args)<3 : 
		print usage
		parser.error("Specify all input and output files")

	angles_filename=args[1]
	ali_filename=args[0]
	output_stack=args[2]

	logid=E2init(sys.argv)
	f=file(angles_filename,'r')
	lines=f.readlines()
	angles=[]
	for line in lines:
		print str.strip(line)
		angles.append(float(line))
		
	print angles

	input_image=EMData()
	input_image.read_image(str(ali_filename),0,HEADER_ONLY)
	nx = input_image.get_attr('nx')
	ny = input_image.get_attr('ny')
	nz = input_image.get_attr('nz')

	print nx,ny,nz
	for z_index in range(0,nz):
		roi=Region(0,0,z_index,nx,ny,1)
		input_image=EMData()
		input_image.read_image(ali_filename,0, HEADER_AND_DATA, roi)
		input_image.set_rotation(90,angles[z_index],90)
		input_image.set_attr('ptcl_repr',1)
		input_image.write_image(output_stack,-1)
		
	E2end(logid)
# If executed as a program
if __name__ == '__main__':
	main()	

