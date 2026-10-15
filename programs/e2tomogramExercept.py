

for job in jobs:
		stepx, stepy, threed = make_tile(job)  # directly call the function
		threed.mult(msk)
		if options.normslice:
			threed.mult(maskz)

		if options.moretile:
			full3d.insert_scaled_sum(
				threed, (int(stepx*step + outx//2), int(stepy*step + outy//2), outz//2)
			)
		else:
			full3d[stepx % 2].insert_clip(
				threed,
				(int(stepx*step + outx//2 - sz//2),
				 int(stepy*step + outy//2 - sz//2),
				 outz//2 - threed["nz"]//2)
			)




    thrds=[threading.Thread(target=make_tile,args=([i])) for i in jobs]
	print("now start threads...")
	thrtolaunch=0
	tsleep=threading.active_count()

	#### non-round fall off. this is mathematically correct but seem to have grid artifacts
	#f=np.zeros((sz,sz))
	x,y=np.indices((sz,sz),dtype=float)/sz-.5
	#f=.25-(x**2+y**2)/2 + ((abs(x)-0.5)**2+(abs(y)-0.5)**2)/2
	f=wtcon+np.exp(-(x**2+y**2)/0.1) - np.exp(-((abs(x)-0.5)**2+(abs(y)-0.5)**2)/0.1)
	f3=np.repeat(f[None, :,:], outz, axis=0)
	msk=from_numpy(f3).copy()
	#####


	while thrtolaunch<len(thrds) or threading.active_count()>tsleep or not jsd.empty():
		if thrtolaunch<len(thrds) :
			while (threading.active_count()==options.threads ) : time.sleep(.1)
			thrds[thrtolaunch].start()
			thrtolaunch+=1
		else: time.sleep(.1)

		if not jsd.empty():
			stepx, stepy, threed=jsd.get()
			#threed["pos"]=[stepx, stepy]
			threed.mult(msk)
			if options.normslice:
				threed.mult(maskz)
			#### insert the cubes to corresponding tomograms
			if options.moretile:
				full3d.insert_scaled_sum(
					threed,(int(stepx*step+outx//2),int(stepy*step+outy//2), outz//2))
			else:
				full3d[stepx%2].insert_clip(
				threed,
				(int(stepx*step+outx//2-sz//2),
				int(stepy*step+outy//2-sz//2),
				outz//2-threed["nz"]//2))


	for t in thrds: t.join()

	if not options.moretile:
		full3d=full3d[0]+full3d[1]
	#full3d.process_inplace("normalize")
	a=full3d.numpy()
	std=np.std(a[[0,-1]])
	full3d.mult(1./std)

	#### skip the tomogram positioning step because there is some contrast difference at the boundary that sometimes breaks the algorithm...
	full3d["zshift"]=0

	apix=imgs[0]["apix_x"]
	full3d["apix_x"]=full3d["apix_y"]=full3d["apix_z"]=apix
	#if options.normslice:
		#full3d.process_inplace("normalize.rows")

	print("Reconstruction done ({:.1f} s). Now writting tomogram to disk...".format(time.time()-time0))
	return full3d
