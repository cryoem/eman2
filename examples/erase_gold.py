#!/usr/bin/env python

from EMAN2 import *

def main():
  progname = os.path.basename(sys.argv[0])
  usage = """prog [options] stack1.hdf stack2.mrcs ...

  Program to erase gold from DDD movies.
  """

  parser = EMArgumentParser(usage=usage,version=EMANVERSION)

  parser.add_argument("--average", default=False, action="store_true", help="Erase gold from average of input stack(s).")
  parser.add_argument("--downsamp", default=1, type=int, help="Downsample the input stack(s). Default is 1, i.e. no downsampling.")
  parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
  parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-2)

  (options, args) = parser.parse_args()

  for arg in args:
    if options.verbose: print(arg)
    frames = load(arg,ds=options.downsample,inv=True)

    if options.average:
      avg = average(frames)
      avg.process_inplace("normalize")
      proc = avg.copy()
      proc.process_inplace("filter.highpass.gauss",{"cutoff_pixels":25})
      proc.process_inplace("filter.lowpass.gauss",{"cutoff_abs":0.1})
      proc.process_inplace("threshold.clampminmax",{"maxval":proc["maximum"],"minval":proc["mean"]+3*proc["sigma"],"tozero":True})
      proc.process_inplace("threshold.binary",{"value":proc["mean"]})
      proc.process_inplace("mask.addshells.gauss",{"val1":0,"val2":8})
      masked = avg-(proc*avg)
      noise = local_noise(masked)
      result = masked + noise * proc
      result *= -1
      result.write_image("{}_proc.hdf")
    
    else:
      for i,f in enumerate(frames):
        f.process_inplace("normalize")
        proc = f.copy()
        proc.process_inplace("filter.highpass.gauss",{"cutoff_pixels":25})
        proc.process_inplace("filter.lowpass.gauss",{"cutoff_abs":0.1})
        proc.process_inplace("threshold.clampminmax",{"maxval":proc["maximum"],"minval":proc["mean"]+3*proc["sigma"],"tozero":True})
        proc.process_inplace("threshold.binary",{"value":proc["mean"]})
        proc.process_inplace("mask.addshells.gauss",{"val1":0,"val2":8})
        masked = f-(proc*f)
        noise = local_noise(masked)
        result = masked + noise * proc
        result *= -1
        result.write_image("{}_proc.hdf",i)


def load(fn,inv=False,ds=1):
  frames = []
  for i in range(EMUtil.get_image_count(fn)):
    f = EMData(fn,i)
    if ds != 1: f.process_inplace("math.fft.resample",{"n":ds})
    if inv: f*=-1
    frames.append(f)
  return frames


def average(frames):
    avgr = Averagers.get("mean")
    avgr.add_image_list(frames)
    return avgr.finish()


def local_noise(img,bs=None):
  if bs == None: 
    bs = np.sqrt(img["nx"]).astype(int)
  localnoise = EMData(img["nx"],img["ny"])
  localnoise.to_zero()
  sz = (bs,bs)
  mx = np.arange(0,img['nx']+bs,bs)
  my = np.arange(0,img['ny']+bs,bs)
  for x in mx:
    for y in my:
      r = img.get_clip(Region(x-bs/2,y-bs/2,bs,bs))
      n = EMData()
      n.to_zero()
      n.process_inplace("math.addnoise",{"noise":r["mean_nonzero"]})
      n.process_inplace("math.addsignoise",{"noise":r["sigma_nonzero"]})
      localnoise.insert_clip(n,(x-bs/2,y-bs/2))
  return localnoise


if __name__ == "__main__":
  main()
