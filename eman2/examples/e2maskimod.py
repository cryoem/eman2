#!/usr/bin/env python
# Muyuan Chen 2015-03
from EMAN2 import *

def main():
	
	usage="""maskimod.py <density map> <imod model> [options]
	
	maskimod.py tomograph.hdf imodmodel.imod --output output.hdf
	Read imod ascii file and mask out the corresponding density map."""
	
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--stack", type=int,help="number of stack (temp)", default=1)
	parser.add_argument("--mskfile", type=str,help="mask file", default=None)
	parser.add_argument("--output", type=str,help="output file", default="output.hdf")
	parser.add_argument("--s1", type=int,help="width of the binary mask", default=0)
	parser.add_argument("--s2", type=int,help="width of the gaussian mask", default=10)
	parser.add_argument("--byobj", action="store_true",help="Generate a density map for each object", default=False )
	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)
	
	densitymap=args[0]
	
	
	data=EMData(densitymap)
	data.mult(-1)
	data.process_inplace("normalize.edgemean")
	mask=data.copy()
	mask.to_zero()
	
	isscatter=False
	
	for stk in range(options.stack):
		
		if options.stack==1:
			imodfile=args[1]
		else:
			imodfile=args[1]+str(stk)+".imod"
		print "Processing ",imodfile,"....."
		
		#imodinfo -a aaaa > ascii2.txt
		
		### check if it is binary file
		textchars = bytearray([7,8,9,10,12,13,27]) + bytearray(range(0x20, 0x100))
		is_binary_string = lambda bytes: bool(bytes.translate(None, textchars))
		if is_binary_string(open(imodfile, 'rb').read(1024)):
			print "Binary file detected. Converting to ascii file..."
			os.system("imodinfo -a {} > {}".format(imodfile,imodfile+'.txt'))
			imodfile=imodfile+'.txt'
		imodin=open(imodfile)
		lines=imodin.readlines()
		imodin.close()
		
		ci=0
		reading=0
		
		for ln in lines:
			if ln.startswith("object"):
				if ci>0 and options.byobj:
					print "Generating map for object ",ci,"......"
					isscatter=False
					outnm=options.output[:-4]+str(ci)+'.hdf'
					
					mask.process_inplace("mask.addshells.gauss",{"val1":options.s1, "val2":options.s2})
					if options.mskfile:
						msknm=options.mskfile[:-4]+str(ci)+'.hdf'
						if(os.path.isfile(msknm)): os.remove(msknm)
						mask.write_image(msknm)
					
					if(os.path.isfile(outnm)): os.remove(outnm)
					mask.mult(data)
					mask.write_image(outnm)
					mask.to_zero()
					#print "done"
				ci+=1
				reading=0
				#print
				
			if ln.startswith("scattered"):
				isscatter=True
				
			if ln.startswith("contour"):
				lt=ln.split()
				reading=int(lt[3])
				lastpt=[-1,-1,-1]
				#print reading
				continue
				
			if reading>0:
				
				reading-=1
				if isscatter:
					point=[int(round(float(i))) for i in ln.split()]
					mask.set_value_at(point[0],point[1],point[2],1)
				
				else:
					point=[float(i) for i in ln.split()]
					if lastpt[0]>=0:
						lst=interpolate(lastpt,point)
						for p in lst:
							mask.set_value_at(p[0],p[1],p[2],1)
						#print lst
					lastpt=point
			
	
	mask.process_inplace("mask.addshells.gauss",{"val1":options.s1, "val2":options.s2})
	
	if options.byobj:
		if options.mskfile: msknm=options.mskfile[:-4]+str(ci)+'.hdf'
		outnm=options.output[:-4]+str(ci)+'.hdf'
		print "Generating map for object ",ci,"......"
	else:
		if options.mskfile: msknm=options.mskfile
		outnm=options.output
		print "Writing output...."
	#display(mask)
	if options.mskfile:
		if(os.path.isfile(msknm)): os.remove(msknm)
		mask.write_image(msknm)
	
	if(os.path.isfile(outnm)): os.remove(outnm)
	
	
	data.mult(mask)
	data.write_image(outnm)
	#launch_childprocess("e2display.py "+options.output)
	E2end(logid)

def interpolate(pt1,pt2):
	
	df=[i[1]-i[0] for i in zip(pt1,pt2)]
	d=Util.hypot3(df[0],df[1],df[2])
	if d==0:
		return []
	#step=[float(i)/d for i in df]
	#print pt1,pt2
	ret=[]
	for u in range(int(round(d))+1):
		t=u/d
		ret.append([int(round((1-t)*i[0]+t*i[1])) for i in zip(pt1,pt2)])
	return ret

if __name__ == '__main__':
	main()
	