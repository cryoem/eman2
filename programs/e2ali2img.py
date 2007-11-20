#!/usr/bin/python

from EMAN2 import EMData,Region

#constants
HEADER_ONLY=True
HEADER_AND_DATA=False


angles_filename='input.tlt'
ali_filename='input.mrc'
output_stack='outpu.img'


f=file(angles_filename,'r')
lines=f.readlines()
angles=[]
for line in lines:
	print str.strip(line)
	angles.append(float(line))
print angles

input_image=EMData()
input_image.read_image(ali_filename,0,HEADER_ONLY)
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

