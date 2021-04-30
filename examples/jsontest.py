#!/usr/bin/env python
from EMAN2 import *
import random

# Steve Ludtke - 04/29/21
# Simple program to test the required precision in string representations of transformation matrices
# and potentially other string representations


# az,alt,phi
orts=[
[0,0,0],[0.1,0,0],[0,0.1,0],[0.1,0.1,0.1],
[90,0,0],[90,0.1,0],[90,0.1,0.1],
[90,90,0],[90,90.1,0],[90,90,0.1],[90,90.1,0.1],[90.1,90.1,0],
[90,10,-90],[90.1,10,-90.1],[90,10.1,-90],[90,0.1,-90.1]]

for i in range(20):
	orts.append([random.uniform(0,90),random.uniform(0,180),random.uniform(0,360)])

js=js_open_dict("/tmp/test.json")
for x in orts:
	t=Transform({"type":"eman","az":x[0],"alt":x[1],"phi":x[2]})
	js[str(x)]=t

js.close()	# note that this just forces a re-open/reread
for x in orts:
	t=Transform({"type":"eman","az":x[0],"alt":x[1],"phi":x[2]})
	r=js[str(x)]
	#(t*r.inverse()).printme()
	print(x,(t*r.inverse()).get_rotation("spin")["omega"],r.get_rotation("eman"))
