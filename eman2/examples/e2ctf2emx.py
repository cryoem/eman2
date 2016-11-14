#!/usr/bin/env python

from EMAN2 import *
from sys import argv

print '<EMX version="1.0">'
for f in sorted(argv[1:]):
	base=base_name(f,nodir=True)[:-5]
	db=js_open_dict(f)
	try:
		ctf=db["ctf"][0]
		#ctf=db["ctf_frame"][1]
	except:
		continue

	print '<micrograph fileName="{}">'.format(base)
	print '<defocusU unit="nm">{:1.3f}</defocusU>'.format(1000.0*(ctf.defocus+ctf.dfdiff/2.0))
	print '<defocusV unit="nm">{:1.3f}</defocusV>'.format(1000.0*(ctf.defocus-ctf.dfdiff/2.0))
	print '<defocusUAngle unit="deg">{:1.3f}</defocusUAngle>'.format(ctf.dfang)
	print '</micrograph>'

print '</EMX>'
