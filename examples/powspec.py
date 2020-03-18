#!/usr/bin/env python

from past.utils import old_div
from builtins import range
from EMAN2 import *
import sys

# because this is in examples but we need access to a function within a program located in bin
sys.path.append(sys.executable.replace("/python",""))

import e2ctf
from math import *
import numpy as np

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options] micrographs/<micrograph_name>

	This program will allow you to serially compute and write the 1D background subtracted power spectra for multiple micrographs.
	
	For the best background subtraction results, it is highly recommended to run this program within a project where CTF correction has been performed (i.e. running e2rawdata.py and/or e2ctf.py).
	"""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_pos_argument(name="images",help="Calculate the power spectra of these images.", default="")#, guitype='filebox', browser="EMBrowserWidget(withmodal=True,multiselect=True)",  row=0, col=0,rowspan=1, colspan=2, mode="eval")
	parser.add_argument("--apix",type=float,help="Angstroms per pixel for all images",default=None)#, guitype='floatbox', row=3, col=0, rowspan=1, colspan=1, mode="eval['self.pm().getAPIX()']")
	parser.add_argument("--constbfactor",type=float,help="Set B-factor to fixed specified value, negative value autofits. Default is 200.0.",default=-1.0)#, guitype='floatbox', row=8, col=0, rowspan=1, colspan=1, mode='eval[-1.0]')
	parser.add_argument("--voltage",type=float,help="Microscope voltage in KV. Default is 200.0.",default=None)# guitype='floatbox', row=3, col=1, rowspan=1, colspan=1, mode="eval['self.pm().getVoltage()']")
	parser.add_argument("--cs",type=float,help="Microscope Cs (spherical aberation). Default is 4.1.",default=None)#, guitype='floatbox', row=4, col=0, rowspan=1, colspan=1, mode="eval['self.pm().getCS()']")
	parser.add_argument("--ac",type=float,help="Amplitude contrast (percentage). Default is 10.0",default=None)#, guitype='floatbox', row=4, col=1, rowspan=1, colspan=1, mode="eval")
	parser.add_argument("--box",type=int,help="Forced box size in grid mode.. ",default=512)#, guitype='intbox', row=5, col=0, rowspan=1, colspan=1, mode="eval")
	parser.add_argument("--oversamp",type=float,help="Oversample power spectrum. Default 1.0.",default=1.0)#, guitype='intbox', row=5, col=0, rowspan=1, colspan=1, mode="eval")
	parser.add_argument("--nobgsub",action="store_true",help="Skip background subtraction",default=False)#, guitype='floatbox', row=3, col=0, rowspan=1, colspan=1, mode="eval['self.pm().getAPIX()']")

	parser.add_argument("--pad",type=int,help="Exclude this many pixels around the edge of input micrograph(s) when computing power spectra.",default=64)#, guitype='intbox', row=5, col=0, rowspan=1, colspan=1, mode="eval")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	(options, args) = parser.parse_args()

	logid=E2init(sys.argv,options.ppid)

	try: os.mkdir('micrographs')
	except: pass

	if options.box<64 :
		print("Box size too small. Using larger box size of 512 pixels).")
		options.box=512
	box=options.box

	voltage=options.voltage
	cs=options.cs
	ac=options.ac
	constbfactor=options.constbfactor

	pad = options.pad

	if options.oversamp < 1: print("Oversamp must be greater than or equal to 1. Using --oversamp=1 instead.")
	oversamp = max(options.oversamp,1)

	for curset,arg in enumerate(args):
		sys.stdout.write("\r{}/{}: {}".format(curset+1,len(args),arg))
		sys.stdout.flush()
		data = EMData(arg,0)

		if options.apix == None:
			try: apix = data["apix_x"]
			except:
				print("Could not find apix. Please specify using --apix or assign the corresponding header parameter, i.e. 'apix_x'.")
				sys.exit(1)
		else: apix = options.apix

		noctfflag = False
		try: 
			ctf=js_open_dict(info_name(arg,nodir=True))["ctf"][0]
			print("")
			print("\tDefocus: {}".format(ctf.defocus))
			print("\tVoltage: {}".format(ctf.voltage))
			print("\tApix: {}".format(ctf.apix))
			print("\tCs: {}".format(ctf.cs))
			print("\tAC: {}".format(ctf.ampcont))
		except:
			print("Could not find CTF parameters. Skipping background subtraction.")
			#ctf = EMAN2Ctf()
			#ctf.from_dict({'defocus':0.0,'dfdiff':0.0,'dfang':0.0,'bfactor':200.0,'ampcont':10.0,'voltage':200.0,'cs':4.1,'apix':apix,'dsbg':-1})
			#if options.voltage!=None : ctf.voltage=options.voltage
			#if options.ac != None: ctf.ampcont = options.ac
			#if options.cs!=None : ctf.cs=options.cs
			#if options.constbfactor>0 : ctf.bfactor=options.constbfactor
			noctfflag = True

		ds=1.0/(apix*box*oversamp)
		nx=old_div(data["nx"],box)-1
		cumulfft=EMData(box,box)
		cumulfft.do_fft_inplace()
		nbx=0
		for ix,x in enumerate(range(nx)):
			for iy,y in enumerate(range(old_div(data["ny"],box)-1)):
				# read the data and make the FFT
				clip=data.get_clip(Region(x*box+old_div(box,2),y*box+old_div(box,2),box,box))
				clip.process_inplace("normalize.edgemean")
				if oversamp>1 : clip=clip.get_clip(Region(0,0,box*oversamp,box*oversamp))		# since we aren't using phases, doesn't matter if we center it or not
				fft=clip.do_fft()
				fft.ri2inten()
				cumulfft+=fft
				nbx+=1

		cumulfft.mult(1.0/(nbx*box**2))
		cumulfft.process_inplace("math.sqrt")
		cumulfft["is_intensity"]=0				# These 2 steps are done so the 2-D display of the FFT looks better. Things would still work properly in 1-D without it

		fftbg=cumulfft.process("math.nonconvex")
		fft1d=cumulfft.calc_radial_dist(old_div(cumulfft.get_ysize(),2),0.0,1.0,1)	# note that this handles the ri2inten averages properly

		if options.nobgsub or noctfflag:
			s=np.arange(0,ds*len(fft1d),ds)

			pwsfn = "micrographs/{}-pws.txt".format(base_name(arg))
			try: os.remove(pwsfn)
			except: pass

			with open(pwsfn,"w") as pwsf:
				pwsf.write("# s(1/A); PWS(s)\n")
				for i in range(len(fft1d)):
					pwsf.write("{}\t{}\n".format(s[i],fft1d[i]))

		else:
			# Compute 1-D curve and background
			bg_1d=e2ctf.low_bg_curve(fft1d,ds)
			ctf.background=bg_1d
			ctf.dsbg=ds

			fft1d=np.asarray(fft1d)
			ds=ctf.dsbg
			bg_1d=list(ctf.background)

			xyd=XYData()

			# Find the minimum value near the origin, which we'll use as a zero (though it likely should not be)
			mv=(fft1d[1],1)
			fz=int(old_div(ctf.zero(0),(ds*2)))
			for lz in range(1,fz):
				mv=min(mv,(fft1d[lz],lz))

			xyd.insort(mv[1],mv[0])

			# now we add all of the zero locations to our XYData object
			for i in range(100):
				z=int(old_div(ctf.zero(i),ds))
				if z>=len(bg_1d)-1: break
				if fft1d[z-1]<fft1d[z] and fft1d[z-1]<fft1d[z+1]: mv=(z-1,fft1d[z-1])
				elif fft1d[z]<fft1d[z+1] : mv=(z,fft1d[z])
				else : mv=(z+1,fft1d[z+1])
				xyd.insort(mv[0],mv[1])

			# new background is interpolated XYData
			ctf.background=[xyd.get_yatx_smooth(i,1) for i in range(len(bg_1d))]

			# if our first point (between the origin and the first 0) is too high, we readjust it once
			bs=[fft1d[i]-ctf.background[i] for i in range(fz)]
			if min(bs)<0 :
				mv=(bs[0],fft1d[0],0)
				for i in range(1,fz): mv=min(mv,(bs[i],fft1d[i],i))
				xyd.set_x(0,mv[2])
				xyd.set_y(0,mv[1])
				
				ctf.background=[xyd.get_yatx_smooth(i,1) for i in range(len(bg_1d))]

			bg1d=np.array(ctf.background)
			r=len(ctf.background)
			s=np.arange(0,ds*r,ds)

			try: bgsub=fft1d-bg1d
			except:
				print("Error computing bgsub on this image")
				continue

			fit=np.array(ctf.compute_1d(len(s)*2,ds,Ctf.CtfType.CTF_AMP))		# The fit curve
			fit=fit*fit			# squared

			# auto-amplitude for b-factor adjustment
			rto,nrto=0,0
			for i in range(int(.04/ds)+1,min(int(0.15/ds),len(s)-1)):
				if bgsub[i]>0 :
					rto+=fit[i]
					nrto+=fabs(bgsub[i])
			if nrto==0 : rto=1.0
			else : rto/=nrto
			fit=[old_div(fit[i],rto) for i in range(len(s))]

			pwsfn = "micrographs/{}-pws.txt".format(base_name(arg))

			try: os.remove(pwsfn)
			except: pass

			with open(pwsfn,"w") as pwsf:
				pwsf.write("# s(1/A); PWS-BG(s); CTF Fit(s); PWS(s); BG(s)\n")
				for i in range(len(fft1d)):
					pwsf.write("{}\t{}\t{}\t{}\t{}\n".format(s[i],bgsub[i],fit[i],fft1d[i],bg1d[i]))

	print("\nPower spectra files saved within the 'micrographs' directory.")

if __name__ == "__main__":
	main()

