from EMAN3 import *
import numpy as np
import sys

imgs=EMData.read_images(sys.argv[1])
inp=[img.numpy() for img in imgs]

for i in range(1,len(imgs)):
	try: c=inp[i]-inp[i-1]
	except:	continue	# change of stage
	print(f"{i}\t{np.fabs(c).mean()}\t{c.std()}\t{c.max()}")

