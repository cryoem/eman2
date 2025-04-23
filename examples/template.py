from math import *
from EMAN3 import *
from EMAN3tensor import *

target=EMData("lamella_3_1_clip_evenmore.hdf",0)
target.mult(-1.0)
nx=target["nx"]
ny=target["ny"]
nz=target["nz"]
targetf=target.do_fft()

template=EMData("long_tmplt.hdf",0)
templatesca=template.process("math.fft.resample",{"n":target["apix_x"]/template["apix_x"]})
nxt1=templatesca["nx"]
clp=good_size_small(templatesca["nx"])
templatescaclip=templatesca.get_clip(Region((nxt1-clp)//2,(nxt1-clp)//2,(nxt1-clp)//2,clp,clp,clp))

owner=EMData(target["nx"],target["ny"],target["nz"])
avg=Averagers.get("minmax",{"max":True,"owner":owner})
orts=[]
for az in range(0,360,5):
	for alt in range(0,180,5):
		trot=templatescaclip.process("xform",{"transform":Transform({"type":"eman","az":az,"alt":alt})}).get_clip(Region((clp-nx)//2,(clp-ny)//2,(clp-nz)//2,nx,ny,nz))
		trotf=trot.do_fft()
		ccf=trotf.calc_ccf(targetf,fp_flag.CIRCULANT,True)
		ccf["ortid"]=len(orts)
		orts.append((az,alt))
		avg.add_image(ccf)
		print(az,alt)

peaks=avg.finish()
peaks.write_image("out.hdf:-1",0)
owner.write_image("outo.hdf:-1",0)

