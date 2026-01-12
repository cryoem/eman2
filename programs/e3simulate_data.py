from EMAN3 import *
from EMAN3jax import *
import gemmi

try: os.mkdir(".jaxcache")
except: pass

# We cache the JIT compilation results to speed up future runs
jax.config.update("jax_compilation_cache_dir", "./.jaxcache")
jax.config.update("jax_persistent_cache_min_entry_size_bytes", -1)
jax.config.update("jax_persistent_cache_min_compile_time_secs", 2)
jax.config.update("jax_persistent_cache_enable_xla_caches", "xla_gpu_per_fusion_autotune_cache_dir")

jax.config.update("jax_default_matmul_precision", "float32")

def main():

	usage="""e3simulate_data.py <3D map/gauss/model>

	This program will take an input PDB/MMCIF/Gauss model/map and produce a set of simulated projection images. At present this uses a
	basic weak phase-approximation model of the CTF, white Gaussian noise and fixed B-factor.

	"""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--output", type=str,help="Output file for simulated projections", default="proj.hdf:8")
	parser.add_argument("--lstout", type=str,help="LST output containing projection metadata", default="proj.lst")
	parser.add_argument("--nptcl",type=int,help="number of simulated particles to generate", default=1000)
	parser.add_argument("--noise", type=float, help="relative noise level to add, default=0. typical ~1", default=0)
	parser.add_argument("--defocus", type=str, help="Defocus range for optional CTF, eg (1,2).Underfocus positive. Default=1,2", default="1,2")
	parser.add_argument("--ctfamp", action="store_true",help="Multiply CTF amplitude (including phase flip)",default=False)
	parser.add_argument("--ctfphaseflip", action="store_true",help="CTF phase flipping",default=None)
	parser.add_argument("--lowpass", type=float, help="Apply a lowpass filter to the images before noise addition, Nyquist=0.5, default=0", default=0)
	parser.add_argument("--ac", type=float, help="Fixed amplitude contrast to use in the CTF as a fractional value, default=0.1  (10%)", default=0.1)
	parser.add_argument("--apix", type=float, help="A/pix override for output data. Assumes input data has correct A/pix when relevant.", default=None)
	parser.add_argument("--boxsize",type=int,help="particle box size in pixels, default=auto", default=None)
	parser.add_argument("--sym", dest = "sym", help = "Specify symmetry - choices are: c<n>, d<n>, h<n>, tet, oct, icos. This does not enforce symmetry on the model/map, it only restricts orientation generation to a single asymmetric unit",default="c1")

	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higher number means higher level of verbosity")

	(options, args) = parser.parse_args()
	jax_set_device(dev=0,maxmem=options.gpuram)

	options.defocus=[float(i) for i in options.defocus.split(",")]
	if len(options.defocus)==1: options.defocus.append(options.defocus[0])
	if len(options.defocus)>2: error_exit("Defocus should be provided as <min>,<max>")
	if options.ctfamp and (options.ctfphaseflip is not None and not options.ctfphaseflip): error_exit("No support for CTF amplitude without phase-flip")

	### Load the map/model once at the beginning
	# Gaussian model
	if args[0].endswith(".txt"):
		inmodel=Gaussians(args[0])
		inmap=None
		modsca=1.0
		if options.apix is None or options.boxsize is None: error_exit("--apix and --boxsize required with Gaussian model inputs")
	# PDB or MMCIF
	elif args[0].endswith(".pdb") or args[0].endswith(".ent") or args[0].endswith(".cif"):
		struc=gemmi.read_structure(args[0])
		inmodel=np.array([[i.atom.pos.x, i.atom.pos.y, i.atom.pos.z, i.atom.element.atomic_number ] for i in struc[0].all()])
		if options.apix is None or options.boxsize is None: error_exit("--apix and --boxsize required with PDB/MMCIF inputs")
		modsca=1.0/(options.apix*options.boxsize)
		inmodel[:][:3]*=modsca
		inmap=None
	# assume some sort of volume
	else:
		inmodel=None
		inmap=EMData(args[0])
		if (options.apix is not None and options.apix!=inmap["apix_x"]) or (options.boxsize is not None and options.boxsize!=inmap["nz"]):
			inmap=inmap.process("xform.scale",{"clip":options.boxsize,"scale":inmap["apix_x"]/options.apix})

	if options.apix<=0:
		if inmodel is None: options.apix=inmap["apix_x"]
		else: error_exit("--apix is required when the input is a model")

	if options.boxsize is None:
		if inmodel is None: options.boxsize=inmap["nz"]
		else: error_exit("--boxsize is required when the input is a model")

	llo=E3init(sys.argv,options.ppid)

	# Number of particles to generate at a time is memory constrained
	step=max(min(1000000000//(options.boxsize**2),options.nptcl//10),2)	# ~4G of RAM
	ngrp=((options.nptcl-1)//step)+1

	# We make a set of CTF objects to use
	ctf=EMAN3CTF(options.defocus[0],0,0,300,2.7,0.1,options.apix)
	dfstep=(options.defocus[1]-options.defocus[0])/(ngrp-1)
	ctfstack=ctf.compute_2d_stack_complex(options.boxsize, "amplitude" if options.ctfamp else "sign", options.defocus, "defocus", apix=options.apix, use_astig=False, ewald_sphere=False, beam_tilt=False, defocus_step=dfstep)

	# if lowpass filter selected apply to CTF images
	if options.lowpass>0:
		lp=jax_gaussfilt_2d(options.boxsize,options.lowpass)
		ctfstack.mult(lp)

	lst=LSXFile(options.lstout)
	for i in range(0,options.nptcl,step):
		eulers = sym_object.gen_orientations("rand",{"n":step,"phitoo":1})
		orts=Orientations.init_from_transforms(eulers)
		if inmap is None:
			projs=inmodel.project_simple(orts,options.boxsize)
		else:
			projs=EMStack2D([inmap.project("standard",{"transform":xform}) for xform in eulers])

		if options.ctfamp or options.ctfphaseflip:
			projsf=projs.do_fft()
			projsf.mult(ctfstack[i//step])
			projs=projsf.do_ift()

		projs.write_images(options.output,n0=i)
		for j in range(i,min(options.nptcl,i+step)):
			ctf.defocus=options.defocus[0]+dfstep*(j-i)
			dct={"xform.projection":eulers[j-i],"ctf":ctf}
			lst[j]=options.output,j,dct

	E3end(llo)
