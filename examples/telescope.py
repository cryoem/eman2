from EMAN2 import *

def rawframes(noframes,outfile=None):
	ret=[]
	for i in range(noframes):
		a=EMData()
		a.read_image("/dev/video0",0)
		ret.append(a)
	if outfile: write_set(ret,outfile)
	return ret
		
def rawframesdisk(noframes,outfile):
	for i in range(noframes):
		a=EMData()
		a.read_image("/dev/video0",0)
		a.write_image(outfile,-1)

def write_set(lst,outfile):
	for i in lst:
		i.write_image(outfile,-1)

a=rawframes(1)
display(a)
