from EMAN2 import *
from emimage import *
import time

def live():
	b=EMData()
	b.read_image("/dev/video0",0)
	c=EMImage(b)
	for i in range(10):
		b.read_image("/dev/video0",0)
		c.setdata(b)
		time.sleep(1)

def rawavg(noframes):
	b=EMData()
	b.read_image("/dev/video0",0)
	for i in range(noframes):
		a=EMData()
		a.read_image("/dev/video0",0)
		b+=a
		if i%10==0 : print i
	return b

def aliavg(noframes):
	b=EMData()
	b.read_image("/dev/video0",0)
	for i in range(noframes):
		a=EMData()
		a.read_image("/dev/video0",0)
#		ba=a.align("translational",b,{},"optvariance",{"matchfilt":1})
		ba=a.align("translational",b,{},"dot",{})
		b+=ba
		print i
	return b


def rawframes(noframes,outfile=None):
	ret=[]
	b=EMData()
	b.read_image("/dev/video0",0)
	for i in range(noframes):
		a=EMData()
		a.read_image("/dev/video0",0)
		ret.append(a)
		b+=a
	if outfile: write_set(ret,outfile)
	return (ret,b)
		
def rawframesdisk(noframes,outfile):
	for i in range(noframes):
		a=EMData()
		a.read_image("/dev/video0",0)
		a.write_image(outfile,-1)
		if i%10==0: print i

def avgframesdisk(noframes,noavg,outfile):
	for i in range(noframes):
	  a=EMData()
	  a.read_image("/dev/video0",0)
	  for j in range(noavg-1):
		b=EMData()
		b.read_image("/dev/video0",0)
		a+=b
	  a.write_image(outfile,-1)
	  if i%10==0: print i

def write_set(lst,outfile):
	for i in lst:
		i.write_image(outfile,-1)

a=rawframes(10)[0][9]
display(a)
