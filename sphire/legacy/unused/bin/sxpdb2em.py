"""1
	while(bigb*fcbig <= 512):   # maximum boxsize = 512
		fcbig +=1
	fcbig -= 1
	if(fcbig == 0):
		print  "Pixel size too small resulting in a box size >512"
		sys.exit()
	"""
"""2
	fftip(outmap)
	# Atom in Fourier space has sigma = 0.41 [1/A]
	outmap = filt_gaussl(outmap,pixelbig*0.41)
	# prepare for Fourier truncation
	if(fcbig > 1):
		outmap = filt_tanl(outmap,0.5/fcbig, 0.02)
		outmap = outmap.FourTruncate(box[0],box[1],box[2], False)
	outmap = fft(outmap)
	"""