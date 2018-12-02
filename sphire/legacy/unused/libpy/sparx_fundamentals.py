"""1
		new_nx = int(nx*sub_rate+0.5)
		new_ny = int(ny*sub_rate+0.5)
		if nz==1:
			new_nz = 1
		else:
			new_nz = int(ny*sub_rate+0.5)
		if ( nx!=ny and nz==1 ):
			nn = max(new_nx, new_ny)
			e = Util.pad(img, nn, nn,  1, 0, 0, 0, "circumference")
			e, kb = prepi(e)
			e = Util.window( e.rot_scale_conv_new(0.0, 0.0, 0.0, kb, sub_rate), new_nx, new_ny, 1, 0,0,0)
		 
		elif ((nx!=ny or nx!=nz or ny!=nz) and nz>1):
			nn = max(new_nx, new_ny,new_nz)
			e = Util.pad(img, nn, nn,  nn, 0, 0, 0, "circumference")
			e, kb = prepi3D(e)
			e = Util.window( e.rot_scale_conv_new_3D(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, kb, sub_rate), new_nx, new_ny, new_nz, 0,0,0)
		else:
			if nz==1:
				e, kb = prepi(Util.pad(img, new_nx, new_ny, 1, 0, 0, 0, "circumference"))
				e = e.rot_scale_conv_new(0.0, 0.0, 0.0, kb, sub_rate)
			else:
				e, kb = prepi3D(Util.pad(img, new_nx, new_ny, new_nz, 0, 0, 0, "circumference"))
				e = e.rot_scale_conv_new_3D(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, kb, sub_rate)
		"""
'''2
def gridrot_shift2D(image, ang = 0.0, sx = 0.0, sy = 0.0, scale = 1.0):
	"""
		Rotate and shift an image using gridding in Fourier space.
	"""
	pass#IMPORTIMPORTIMPORT from EMAN2 import Processor
	pass#IMPORTIMPORTIMPORT from fundamentals import fftip, fft

	nx = image.get_xsize()
	# split shift into integer and fractional parts
	isx = int(sx)
	fsx = sx - isx
	isy = int(sy)
	fsy = sy - isy
	# prepare 
	npad = 2
	N = nx*npad
	K = 6
	alpha = 1.75
	r = nx/2
	v = K/2.0/N
	kb = Util.KaiserBessel(alpha, K, r, v, N)

	image1 = image.copy()  # This step is needed, otherwise image will be changed outside the function
	# divide out gridding weights
	image1.divkbsinh(kb)
	# pad and center image, then FFT
	image1 = image1.norm_pad(False, npad)
	fftip(image1)
	# Put the origin of the (real-space) image at the center
	image1.center_origin_fft()
	# gridding rotation
	image1 = image1.fouriergridrot2d(ang, scale, kb)
	if(fsx != 0.0 or fsy != 0.0):
		params = {"filter_type" : Processor.fourier_filter_types.SHIFT,	"x_shift" : float(fsx), "y_shift" : float(fsy), "z_shift" : 0.0 }
		image1 = Processor.EMFourierFilter(image1, params)
	# put the origin back in the corner
	image1.center_origin_fft()
	# undo FFT and remove padding (window)
	image1 = fft(image1)
	image1 = image1.window_center(nx)
	Util.cyclicshift(image1,{"dx":isx,"dy":isy,"dz":0})
	return image1
'''
"""3
def rotmatrix(phi,theta,psi):
	pass#IMPORTIMPORTIMPORT from math import sin,cos,radians
	rphi   = radians(phi)
	rtheta = radians(theta)
	rpsi   = radians(psi)
	mat = [[0.0]*3,[0.0]*3,[0.0]*3]

	mat[0][0] =  cos(rpsi)*cos(rtheta)*cos(rphi) - sin(rpsi)*sin(rphi)
	mat[1][0] = -sin(rpsi)*cos(rtheta)*cos(rphi) - cos(rpsi)*sin(rphi)
	mat[2][0] =            sin(rtheta)*cos(rphi)


	mat[0][1] =  cos(rpsi)*cos(rtheta)*sin(rphi) + sin(rpsi)*cos(rphi)
	mat[1][1] = -sin(rpsi)*cos(rtheta)*sin(rphi) + cos(rpsi)*cos(rphi)
	mat[2][1] =            sin(rtheta)*sin(rphi)


	mat[0][2] = -cos(rpsi)*sin(rtheta)
	mat[1][2] =  sin(rpsi)*sin(rtheta)
	mat[2][2] =            cos(rtheta)
	return mat
"""
"""4

def mulmat_np(m1,m2):
	pass#IMPORTIMPORTIMPORT import numpy as np
	mat1 = np.matrix(m1,dtype="f8")
	mat2 = np.matrix(m2,dtype="f8")
	mat1 = np.array(mat1*mat2)
	return [list(q) for q in mat1]

def rotmatrix_np(phi,theta,psi):
	pass#IMPORTIMPORTIMPORT import numpy as np
	mat = np.matrix(((0.,0.,0.),(0.,0.,0.),(0.,0.,0.)), dtype = "f8")
	rphi   = np.radians(np.float64(phi))
	rtheta = np.radians(np.float64(theta))
	rpsi   = np.radians(np.float64(psi))
	cosphi = np.cos(rphi)
	sinphi = np.sin(rphi)
	costheta = np.cos(rtheta)
	sintheta = np.sin(rtheta)
	cospsi = np.cos(rpsi)
	sinpsi = np.sin(rpsi)

	mat[0,0] =  cospsi*costheta*cosphi - sinpsi*sinphi
	mat[1,0] = -sinpsi*costheta*cosphi - cospsi*sinphi
	mat[2,0] =            sintheta*cosphi


	mat[0,1] =  cospsi*costheta*sinphi + sinpsi*cosphi
	mat[1,1] = -sinpsi*costheta*sinphi + cospsi*cosphi
	mat[2,1] =            sintheta*sinphi


	mat[0,2] = -cospsi*sintheta
	mat[1,2] =  sinpsi*sintheta
	mat[2,2] =            costheta
	return mat


def recmat_np(mat):
	#from math import np.arccos,np.np.arcsin,np.arctan2,degrees,pi
	pass#IMPORTIMPORTIMPORT import numpy as np
	'''
	def sign(x):
		if( x >= 0.0 ): return 1
		else:  return -1
	mat = [[0.0]*3,[0.0]*3,[0.0]*3]
	# limit precision
	for i in range(3):
		for j in range(3):
			mat[i,j] = inmat[i,j]
			#if(abs(inmat[i,j])<1.0e-8):  mat[i,j] = 0.0
			#else: mat[i,j] = inmat[i,j]
	for i in range(3):
		for j in range(3):  print  "     %14.8f"%mat[i,j],
		print ""
	'''
	if(mat[2,2] == 1.0):
		theta = 0.0
		psi = 0.0
		if( mat[0,0] == 0.0 ):
			phi = np.np.arcsin(mat[0,1])
		else:
			phi = np.arctan2(mat[0,1],mat[0,0])
	elif(mat[2,2] == -1.0):
		theta = pi
		psi = 0.0
		if(mat[0,0] == 0.0):
			phi = np.np.arcsin(-mat[0,1])
		else:
			phi = np.arctan2(-mat[0,1],-mat[0,0])
	else:
		theta = np.arccos(mat[2,2])
		st = np.sign(theta)
		#print theta,st,mat[2,0]
		if(mat[2,0] == 0.0):
			if( st != np.sign(mat[2,1]) ):
				phi = 1.5*pi
			else:
				phi = 0.5*pi
		else:
			phi = np.arctan2(st*mat[2,1], st*mat[2,0])

		#print theta,st,mat[0,2],mat[1,2]
		if(mat[0,2] == 0.0):
			if( st != np.sign(mat[1,2]) ):
				psi = 1.5*pi
			else:
				psi = 0.5*pi
		else:
			psi = np.arctan2(st*mat[1,2], -st*mat[0,2])
	pi2 = 2*np.pi
	#return  degrees(round(phi,10)%pi2),degrees(round(theta,10)%pi2),degrees(round(psi,10)%pi2)
	#return  degrees(round(phi,10)%pi2)%360.0,degrees(round(theta,10)%pi2)%360.0,degrees(round(psi,10)%pi2)%360.0
	return  np.degrees(np.mod(phi,pi2)),np.degrees(np.mod(theta,pi2)),np.degrees(np.mod(psi,pi2))
"""
"""5
			#  These angles were translated from eman to spider, but the do not agree with definitions of subunit above
			for l1 in range(30,271,120):
				for l2 in range(30,271,120):
					self.symangles.append([float(l1),lvl1,float(l2)])
			"""
"""6
	def reduce_normal(self, phi, theta, psi, inc_mirror):
		pass#IMPORTIMPORTIMPORT from math import degrees, radians, sin, cos, tan, atan, acos, sqrt
		return False
	"""
"""Create a list of Euler angles suitable for projections.7
			   method is either 'S' - for Saff algorithm
							  or   'P' - for Penczek '94 algorithm
					  'S' assumes phi1<phi2 and phi2-phi1>> delta ;
			   phiEqpsi  - set this to 'Minus', if you want psi=-phi;
			"""
"""8
		#  Helical symmetry should not be here
		elif(self.sym[0]  == "s"):

			#if symetry is "s", deltphi=delta, theata intial=theta1, theta end=90, delttheta=theta2
			# for helical, theta1 cannot be 0.0
			if theta1 > 90.0:
				ERROR('theta1 must be less than 90.0 for helical symmetry', 'even_angles', 1)
			if theta1 == 0.0: theta1 =90.0
			theta_number = int((90.0 - theta1)/theta2)
			#for helical, symmetry = s or scn
			cn = int(self.sym[2:])
			for j in range(theta_number,-1, -1):

				if( j == 0):
					if (self.sym[1] =="c"):
						if cn%2 == 0:
							k=int(359.99/cn/delta)
						else:
							k=int(359.99/2/cn/delta)
					elif (self.sym[1] =="d"):
						if cn%2 == 0:
							k=int(359.99/2/cn/delta)
						else:
							k=int(359.99/4/cn/delta)
					else:
						ERROR("For helical strucutre, we only support scn and sdn symmetry","even_angles",1)

				else:
					if (self.sym[1] =="c"):
						k=int(359.99/cn/delta)
					elif (self.sym[1] =="d"):
						k=int(359.99/2/cn/delta)

				for i in range(k+1):
						angles.append([i*delta,90.0-j*theta2,90.0])
			"""
def absi(e):
	if e.is_complex():
		e.set_attr_dict({"is_complex_ri":1})
	return e.absi()

# Autocorrelation functions
def acf(e, center=True):
	"""
		Name
			acf - calculate the circulant autocorrelation function of an image
		Input
			e: input image, can be either real or Fourier
			center: if set to True (default), the origin of the result is at the center; if set to False, the origin is at (0,0).
		Output
			circulant autocorrelation function of the input image. Real.
	"""
	pass#IMPORTIMPORTIMPORT from EMAN2 import autocorrelation, fp_flag
	return EMAN2_cppwrap.autocorrelation(e, EMAN2_cppwrap.fp_flag.CIRCULANT, center)

def acfn(e, center=True):
	"""
		Name
			acfn - calculate the normalized circulant autocorrelation function of an image
		Input
			e: input image (real)
			center: if set to True (default), the origin of the result is at the center; if set to False, the origin is at (0,0).
		Output
			normalized circulant autocorrelation function of an input image. Real.
	"""
	pass#IMPORTIMPORTIMPORT from EMAN2 import autocorrelation, fp_flag
	return EMAN2_cppwrap.autocorrelation(e, EMAN2_cppwrap.fp_flag.CIRCULANT_NORMALIZED, center)

def acfp(e, center=True):
	pass#IMPORTIMPORTIMPORT from EMAN2 import autocorrelation, fp_flag
	return EMAN2_cppwrap.autocorrelation(e, EMAN2_cppwrap.fp_flag.PADDED, center)

def acfnp(e, center=True):
	"""
		Name
			acfnp - calculate the normalized autocorrelation function of an image
		Input
			e: input image (real)
			center: if set to True (default), the origin of the result is at the center; if set to False, the origin is at (0,0).
		Output
			normalized autocorrelation function of the input image. Real. 
	"""
	pass#IMPORTIMPORTIMPORT from EMAN2 import autocorrelation, fp_flag
	return EMAN2_cppwrap.autocorrelation(e, EMAN2_cppwrap.fp_flag.PADDED_NORMALIZED, center)

def acfpl(e, center=True):
	"""
		Name
			acfpl - calculate the normalized autocorrelation function of an image
		Input
			e: input image (real)
			center: if set to True (default), the origin of the result is at the center; if set to False, the origin is at (0,0).
		Output
			normalized autocorrelation function of the input image. Real. 
			
	"""
	pass#IMPORTIMPORTIMPORT from EMAN2 import autocorrelation, fp_flag
	return EMAN2_cppwrap.autocorrelation(e, EMAN2_cppwrap.fp_flag.PADDED_LAG, center)

def acfnpl(e, center=True):
	"""
		Name
			acfnpl - calculate the normalized autocorrelation function of an image
		Input
			e: input image (real)
			center: if set to True (default), the origin of the result is at the center; if set to False, the origin is at (0,0).
		Output
			autocorrelation function of the input image. Real. 
	"""
	pass#IMPORTIMPORTIMPORT from EMAN2 import autocorrelation, fp_flag
	return EMAN2_cppwrap.autocorrelation(e, EMAN2_cppwrap.fp_flag.PADDED_NORMALIZED_LAG, center)

def __buildweights(m, kb):
	weights = EMAN2_cppwrap.EMData()
	weights.set_size(m,m,1)
	for iy in range(m):
		wy = kb.sinhwin(iy-m//2)
		for ix in range(m):
			wx = kb.sinhwin(ix-m//2)
			weights.set_value_at(ix,iy,wx*wy)
	return weights

# shortcuts to Fourier product functions
# Correlation functions
def ccf(e, f, center=True):
	"""
	Return the circulant cross-correlation function of images e and f.
	Input images may be real or complex.  Output image is real.
	1-D, 2-D, or 3-D images supported.
	"""
	pass#IMPORTIMPORTIMPORT from EMAN2 import correlation, fp_flag
	return EMAN2_cppwrap.correlation(e,f,EMAN2_cppwrap.fp_flag.CIRCULANT, center)

def ccfn(e, f, center=True):
	"""
		Name
			ccfn - calculate the normalized circulant cross-correlation function between two images.
		Input
			e: input image (real)
			ref: second input image (real) (in the alignment problems, it should be the reference image).
			center: if set to True (default), the origin of the result is at the center
		Output
			normalized circulant cross-correlation function between image and ref. Real.
	"""
	pass#IMPORTIMPORTIMPORT from EMAN2 import correlation, fp_flag
	return EMAN2_cppwrap.correlation(e,f,EMAN2_cppwrap.fp_flag.CIRCULANT_NORMALIZED, center)

def ccfp(e, f, center=True):
	"""
		Name
			ccfp - calculate the cross-correlation function between two images
		Input
			e: input image (real)
			ref: second input image (real)
			center: if set to True (default), the origin of the result is at the center
		Output
			cross-correlation function between image and ref. Real.
	"""
	pass#IMPORTIMPORTIMPORT from EMAN2 import correlation, fp_flag
	return EMAN2_cppwrap.correlation(e,f,EMAN2_cppwrap.fp_flag.PADDED, center)

def ccfnp(e, f, center=True):
	"""
		Name
			ccfnp - calculate the normalized cross-correlation function between two images.
		Input
			e: input image (real)
			ref: second input image (real).
			center: if set to True (default), the origin of the result is at the center
		Output
			normalized cross-correlation function between image and ref. Real.
	"""
	pass#IMPORTIMPORTIMPORT from EMAN2 import correlation, fp_flag
	return EMAN2_cppwrap.correlation(e,f,EMAN2_cppwrap.fp_flag.PADDED_NORMALIZED, center)

def ccfpl(e, f, center=True):
	"""
		Name
			ccfpl - calculate the cross-correlation function between two images	
		Input
			e: input image (real)
			ref: second input image (real) 
			center: if set to True (default), the origin of the result is at the center
		Output
			cross-correlation function between image and ref. Real. 
	"""
	pass#IMPORTIMPORTIMPORT from EMAN2 import correlation, fp_flag
	return EMAN2_cppwrap.correlation(e,f, EMAN2_cppwrap.fp_flag.PADDED_LAG, center)

def ccfnpl(e, f, center=True):
	"""
		Name
			ccfnpl - calculate the normalized cross-correlation function between two images
		Input
			e: input image (real)
			ref: second input image (real) 
			center: if set to True (default), the origin of the result is at the center
	"""
	pass#IMPORTIMPORTIMPORT from EMAN2 import correlation, fp_flag
	return EMAN2_cppwrap.correlation(e,f,EMAN2_cppwrap.fp_flag.PADDED_NORMALIZED_LAG, center)
    
# Convolution functions
def cnv(e, f, center=True):
	"""
		Name
			cnv - calculate the circulant convolution function between two images
		Input
			e: input image, can be either real or Fourier
			ref: second input image, can be either real or Fourier.
			center: if set to True (default), the origin of the result is at the center
		Output
			circulant convolution function between image and ref. Real.
	"""
	pass#IMPORTIMPORTIMPORT from EMAN2 import convolution, fp_flag
	return EMAN2_cppwrap.convolution(e,f,EMAN2_cppwrap.fp_flag.CIRCULANT, center)

def cnvn(e, f, center=True):
	"""
		Name
			cnvn - calculate the normalized circulant convolution function between two images
		Input
			e: input image (real)
			ref: second input image (real).
			center: if set to True (default), the origin of the result is at the center	
		Output
			normalized circulant convolution function between image and ref. Real. 
	"""
	pass#IMPORTIMPORTIMPORT from EMAN2 import convolution, fp_flag
	return EMAN2_cppwrap.convolution(e,f,EMAN2_cppwrap.fp_flag.CIRCULANT_NORMALIZED, center)

def cnvp(e, f, center=True):
	"""
		Name
			cnvp - calculate the convolution function between two images 
		Input
			e: input image (real)
			ref: second input image (real).
			center: if set to True (default), the origin of the result is at the center			
		Output
			convolution function between image and ref. Real.
	"""
	pass#IMPORTIMPORTIMPORT from EMAN2 import convolution, fp_flag
	return EMAN2_cppwrap.convolution(e,f,EMAN2_cppwrap.fp_flag.PADDED, center)

def cnvnp(e, f, center=True):
	"""
		Name
			cnvnp - calculate the normalized convolution function between two images
		Input
			e: input image (real)
			ref: second input image (real) 
			center: if set to True (default), the origin of the result is at the center
		Output
			normalized convolution function between image and ref. Real. 
	"""
	pass#IMPORTIMPORTIMPORT from EMAN2 import convolution, fp_flag
	return EMAN2_cppwrap.convolution(e,f,EMAN2_cppwrap.fp_flag.PADDED_NORMALIZED, center)

def cnvpl(e, f, center=True):
	"""
		Name
			cnvpl - calculate the convolution function between two images
		Input
			e: input image (real)
			ref: second input image (real) 
			center: if set to True (default), the origin of the result is at the center
		Output
			convolution function between image and ref. Real. 
	"""
	pass#IMPORTIMPORTIMPORT from EMAN2 import convolution, fp_flag
	return EMAN2_cppwrap.convolution(e,f,EMAN2_cppwrap.fp_flag.PADDED_LAG, center)

def cnvnpl(e, f, center=True):
	"""
		Name
			cnvnpl - calculate the normalized convolution function between two images
		Input
			e: input image (real)
			ref:second input image (real) 
			center: if set to True (default), the origin of the result is at the center
		Output
			convolution function between image and ref. Real.
	"""
	pass#IMPORTIMPORTIMPORT from EMAN2 import convolution, fp_flag
	return EMAN2_cppwrap.convolution(e,f,EMAN2_cppwrap.fp_flag.PADDED_NORMALIZED_LAG, center)
    
    
# Selfcorrelation functions
def scf(e, center=True):
	"""
		Name
			scf - calculate the circulant self-correlation function of an image
		Input
			e: input image, can be either real or Fourier
			center: if set to True (default), the origin of the result is at the center
		Output
			circulant self-correlation function of the input image. Real.
	"""
	pass#IMPORTIMPORTIMPORT from EMAN2 import self_correlation, fp_flag
	return EMAN2_cppwrap.self_correlation(e, EMAN2_cppwrap.fp_flag.CIRCULANT, center)

def scfn(e, center=True):
	"""
		Name
			scfn - calculate the normalized circulant self-correlation function
		Input
			e: input image, can be either real or Fourier
			center: if set to True (default), the origin of the result is at the center
		Output
			normalized circulant self-correlation function of an input image. Real.
	"""
	pass#IMPORTIMPORTIMPORT from EMAN2 import self_correlation, fp_flag
	return EMAN2_cppwrap.self_correlation(e, EMAN2_cppwrap.fp_flag.CIRCULANT_NORMALIZED, center)

def scfp(e, center=True):
	"""
		Name
			scfp - calculate the self-correlation function of an image
		Input
			e: input image (real)
			center: if set to True (default), the origin of the result is at the center
		Output
			self-correlation function of the input image. Real. 
	"""
	pass#IMPORTIMPORTIMPORT from EMAN2 import self_correlation, fp_flag
	return EMAN2_cppwrap.self_correlation(e, EMAN2_cppwrap.fp_flag.PADDED, center)

def scfnp(e, center=True):
	"""
		Name
			scfnp - calculate the normalized self-correlation function of an image
		Input
			e: input image (real)
			center: if set to True (default), the origin of the result is at the center
		Output
			normalized self-correlation function of the input image. Real.
	"""
	pass#IMPORTIMPORTIMPORT from EMAN2 import self_correlation, fp_flag
	return EMAN2_cppwrap.self_correlation(e, EMAN2_cppwrap.fp_flag.PADDED_NORMALIZED, center)

def scfpl(e, center=True):
	"""
		Name
			scfpl - calculate the self-correlation function of an image
		Input
			e: input image (real)
			center: if set to True (default), the origin of the result is at the center
		Output
			self-correlation function of the input image. Real.
	"""
	pass#IMPORTIMPORTIMPORT from EMAN2 import self_correlation, fp_flag
	return EMAN2_cppwrap.self_correlation(e, EMAN2_cppwrap.fp_flag.PADDED_LAG, center)

def scfnpl(e, center=True):
	"""
		Name
			scfnpl - calculate the normalized self-correlation function of an image
		Input
			e: input image (real)
			center: if set to True (default), the origin of the result is at the center
		Output
			self-correlation function of the input image. Real.
	"""
	pass#IMPORTIMPORTIMPORT from EMAN2 import self_correlation, fp_flag
	return EMAN2_cppwrap.self_correlation(e, EMAN2_cppwrap.fp_flag.PADDED_NORMALIZED_LAG, center)
 
def fft_(e, npad=1):
	"""Out-of-place fft / ift
		No padding performed, and fft-extension along x removed after ift. Zhong added in July,5,06
	"""
	if (e.is_complex()):
		return e.do_ift()
	else:
		if npad > 1:
			f = e.norm_pad(False, npad)
			f.do_fft_inplace()
			return f
		else:
			f = e.norm_pad(False, 1)
			f.do_fft_inplace()
			return f

def image_decimate(img, decimation=2, fit_to_fft = True, frequency_low=0, frequency_high=0):
	"""
		Window 2D image to FFT-friendly size, apply Butterworth low pass filter,
		and decimate image by integer factor
	"""
	pass#IMPORTIMPORTIMPORT from filter       import filt_btwl
	pass#IMPORTIMPORTIMPORT from fundamentals import smallprime
	pass#IMPORTIMPORTIMPORT from utilities    import get_image
	if type(img)     == str:	img=utilities.get_image(img)
	nz       = img.get_zsize()
	if( nz > 1):                    global_def.ERROR("This command works only for 2-D images", "image_decimate", 1)
	if decimation    <= 1  : 	global_def.ERROR("Improper decimation ratio", "image_decimate", 1)
	if(decimation    == 1.0): 	return  img.copy()
	if frequency_low <= 0  :	
		frequency_low     = 0.5/decimation-0.02
		if frequency_low <= 0 : global_def.ERROR("Butterworth pass-band frequency is too low","image_decimation",1)			
		frequency_high    = min(0.5/decimation + 0.02, 0.499)
	if fit_to_fft:
		nx       = img.get_xsize()
		ny       = img.get_ysize()
		nx_fft_m = smallprime(nx)
		ny_fft_m = smallprime(ny)
		e        = EMAN2_cppwrap.Util.window(img, nx_fft_m, ny_fft_m, 1, 0,0,0)
		e        = filter.filt_btwl(e, frequency_low, frequency_high)
	else:
		e        = filter.filt_btwl(img, frequency_low, frequency_high)
	return EMAN2_cppwrap.Util.decimate(e, int(decimation), int(decimation), 1)

def fdownsample(img, sub_rate=0.5, RetReal = True):
	"""
		resample image based on the value of sub_rate.
		the input image can be either 2D image or 3D volume.
		sub_rate < 1.0, subsampling the image.
		sub_rate > 1.0, upsampling the image using new gridding interpolation.
		fit_to_fft will change the ouput image size to an fft_friendly size
	"""

	pass#IMPORTIMPORTIMPORT from fundamentals import fdecimate
	pass#IMPORTIMPORTIMPORT from utilities    import get_pixel_size, set_pixel_size

	if type(img) == str:
		pass#IMPORTIMPORTIMPORT from utilities    import get_image
		img = utilities.get_image(img)
	nx = img.get_xsize()
	if img.is_complex():
		nx -= (2-nx%2)
	ny = img.get_ysize()
	nz = img.get_zsize()
	if( ny == 1):  global_def.ERROR("Only 2D or 3D images allowed","resample",1)
	if sub_rate == 1.0: return  img.copy()
	elif sub_rate < 1.0:
		nnx = int(nx*sub_rate+0.5)
		nny = int(ny*sub_rate+0.5)
		nnz = int(nz*sub_rate+0.5)
		e = fdecimate(img, nnx, nny, nnz, RetReal = RetReal)
	else:  #  sub_rate>1
		global_def.ERROR("fdownsample","upscaling not implemented",1)
		"""Multiline Comment0"""

	# Automatically adjust pixel size for ctf parameters
	pass#IMPORTIMPORTIMPORT from utilities import get_pixel_size, set_pixel_size
	apix = utilities.get_pixel_size(e)
	apix /= sub_rate
	utilities.set_pixel_size(e, apix)
	cc = e.get_attr_default("xform.projection", None)
	if cc:
		cp = cc.get_params("spider")
		cp["tx"] *= sub_rate
		cp["ty"] *= sub_rate
		pass#IMPORTIMPORTIMPORT from utilities import set_params_proj
		utilities.set_params_proj(e, [cp["phi"], cp["theta"], cp["psi"], -cp["tx"], -cp["ty"]]) # have to invert as set inverts them again
	cc = e.get_attr_default("xform.align2d", None)
	if cc:
		cp = cc.get_params("2D")
		cp["tx"] *= sub_rate
		cp["ty"] *= sub_rate
		pass#IMPORTIMPORTIMPORT from utilities import set_params2D
		utilities.set_params2D(e, [cp["alpha"], cp["tx"], cp["ty"], cp["mirror"], cp["scale"]])

	return 	e


def prepf(image, npad = 2):
	"""
		Name
			prepf - prepare 2-D image for Fourier interpolation rotation/shift
		Input
			image: input image that is going to be rotated and shifted using fourier_rotate_shift2d
		Output
			imageft: Fourier space image prepared for Fourier interpolation rotation/shift
	"""
	
	cimage = EMAN2_cppwrap.Util.pad(image,2*(image.get_xsize()),2*(image.get_ysize()),1,0,0,0,"0.0")
	cimage.set_attr("npad",npad)
	cimage.div_sinc(1)
	cimage = cimage.norm_pad(False, 1)
	cimage.do_fft_inplace()
	cimage.center_origin_fft()
	cimage.fft_shuffle()
	cimage.set_attr("npad",npad)
	return cimage

def prep_refim_gridding(refim, wr, numr, mode = "F"):
	pass#IMPORTIMPORTIMPORT from fundamentals import prepi
	nx = refim.get_xsize()
	ny = refim.get_ysize()
	cnx = nx//2+1
	cny = ny//2+1
	#precalculate rings
	temp,kb = prepi(refim)
	crefim = EMAN2_cppwrap.Util.Polar2Dmi(temp, cnx, cny, numr, mode, kb)
	EMAN2_cppwrap.Util.Normalize_ring(cimage, numr, 0 )
	EMAN2_cppwrap.Util.Frngs(crefim, numr)
	EMAN2_cppwrap.Util.Applyws(crefim, numr, wr)
	return  crefim,kb

def prepg(image, kb):
	"""
		Name
			prepg - prepare 2-D image for rotation/shift using gridding method.
		Input
			image: input image that is going to be rotated and shifted using rtshgkb
			kb: interpolants (tabulated Kaiser-Bessel function)
		Output
			imageft: image prepared for gridding rotation and shift
	"""
	pass#IMPORTIMPORTIMPORT from EMAN2 import Processor

	M = image.get_xsize()
	# padd two times
	npad = 2
	N = M*npad
	# support of the window
	K = 6
	alpha = 1.75
	r = M/2
	v = K/2.0/N
	# first pad it with zeros in Fourier space
	o = image.FourInterpol(2*M,2*M,1,0)
	params = {"filter_type" : EMAN2_cppwrap.Processor.fourier_filter_types.KAISER_SINH_INVERSE,
		  "alpha" : alpha, "K":K,"r":r,"v":v,"N":N}
	q = EMAN2_cppwrap.Processor.EMFourierFilter(o,params)
	return  fft(q)
	
def rot_avg(e):
	"""Rotational average.
	   Returns a 1-D image containing a rotational average of image e.
	"""
	return e.rotavg()

def rot_avg_image(image_to_be_averaged):
	"""
	Rotational average
	Returns a 2-D or 3-D image containing a rotational average of image e
	"""
	pass#IMPORTIMPORTIMPORT import types
	pass#IMPORTIMPORTIMPORT from utilities import get_im
	if type(image_to_be_averaged) is bytes: image_to_be_averaged = utilities.get_im(image_to_be_averaged)
	return image_to_be_averaged.rotavg_i()

def ro_textfile(e, filename, helpful_string=""):
	"""Rotational average stored as a text file.
	   Saves a text file (suitable for gnuplot) of the rotational average of e.
	"""
	out = open(filename, "w")
	out.write("#Rotational average: %s\n" % (helpful_string));
	f = e.rotavg()
	nr = f.get_xsize()
	for ir in range(nr):
		out.write("%d\t%12.5g\n" % (ir, f.get_value_at(ir)))
	out.close()

def rops(e):
	"""Rotational average of the power spectrum.
	   Returns a 1-D image containing a rotational average
	   of the periodogram of image e.
	"""
	pass#IMPORTIMPORTIMPORT from EMAN2 import periodogram
	ps = EMAN2_cppwrap.periodogram(e)
	return ps.rotavg()

def rops_textfile(e, filename, helpful_string="", lng = False):
	"""Rotational average of the periodogram stored as a text file.
	   Saves a text file (suitable for gnuplot) of the rotational average 
	   of the periodogram of image e.
	"""
	pass#IMPORTIMPORTIMPORT from EMAN2 import periodogram
	out = open(filename, "w")
	if helpful_string != "": out.write("#Rotational average: %s\n" % (helpful_string))
	ps = EMAN2_cppwrap.periodogram(e)
	f = ps.rotavg()
	nr = f.get_xsize()
	table = [0.0]*nr
	for ir in range(nr): table[ir] = f.get_value_at(ir)
	if lng:
		pass#IMPORTIMPORTIMPORT from math import log
		for ir in range(1,nr): table[ir] = numpy.log(table[ir])
		table[0] = table[1]
	for ir in range(nr): out.write("%d\t%12.5g\n" % (ir, table[ir]))
	out.close()
	
def rops_dir(indir, output_dir = "1dpw2_dir"):
	"""
		Calculate 1D rotationally averaged power spectra from
		image stack listed in a directory
	"""
	pass#IMPORTIMPORTIMPORT from EMAN2 import periodogram
	pass#IMPORTIMPORTIMPORT import os
	flist = os.listdir(indir)
	print(flist)
	if os.path.exists(output_dir) is False: os.mkdir(output_dir)
	for i, v in enumerate(flist):
		(filename, filextension) = os.path.splitext(v)
		nima = EMAN2_cppwrap.EMUtil.get_image_count(os.path.join(indir,v))
		print(nima)
		for im in range(nima):
			e = EMAN2_cppwrap.EMData()
			file_name = os.path.join(indir,v)
			e.read_image(file_name, im)
			tmp1 = EMAN2_cppwrap.periodogram(e)
			tmp  = tmp1.rotavg()
			if im == 0:
				sum_ima  = utilities.model_blank(tmp.get_xsize())
				sum_ima += tmp
			else :  sum_ima += tmp
		table = []
		nr = sum_ima.get_xsize()
		for ir in range(nr):  table.append([sum_ima.get_value_at(ir)])
		utilities.drop_spider_doc(os.path.join(output_dir, "1dpw2_"+filename+".txt"), table)


def rotshift2dg(image, ang, dx, dy, kb, scale = 1.0):
	"""Rotate and shift an image using gridding
	"""
	pass#IMPORTIMPORTIMPORT from math import radians
	pass#IMPORTIMPORTIMPORT from EMAN2 import Processor

	M = image.get_xsize()
	alpha = 1.75
	K = 6
	N = M*2  # npad*image size
	r = M/2
	v = K/2.0/N
	# first pad it with zeros in Fourier space
	o = image.FourInterpol(N,N,1,0)
	params = {"filter_type" : EMAN2_cppwrap.Processor.fourier_filter_types.KAISER_SINH_INVERSE,
	          "alpha" : alpha, "K":K,"r":r,"v":v,"N":N}
	q = EMAN2_cppwrap.Processor.EMFourierFilter(o,params)
	o = fft(q)

	# gridding rotation
	return o.rot_scale_conv(numpy.radians(ang), dx, dy, kb, scale)

"""Multiline Comment1"""

def ft2polargrid(image, ring_length, nb, ne):
	"""
		resample to polar coordinates using gridding in Fourier space.
	"""
	pass#IMPORTIMPORTIMPORT from EMAN2 import Processor
	pass#IMPORTIMPORTIMPORT from fundamentals import fftip, fft

	nx = image.get_xsize()
	# prepare 
	npad = 2
	N = nx*npad
	K = 6
	alpha = 1.75
	r = nx/2
	v = K/2.0/N
	kb = EMAN2_cppwrap.Util.KaiserBessel(alpha, K, r, v, N)

	image1 = image.copy()  # This step is needed, otherwise image will be changed outside the function
	# divide out gridding weights
	image1.divkbsinh(kb)
	# pad and center image, then FFT
	image1 = image1.norm_pad(False, npad)
	fftip(image1)
	# Put the origin of the (real-space) image at the center
	image1.center_origin_fft()
	return image1.ft2polargrid(ring_length, nb, ne, kb)


def rot_shift3D_grid(img, phi=0.0, theta=0.0, psi=0.0, sx=0.0, sy=0.0, sz=0.0, scale=1.0, kb=None, mode="background", wrap=False):
	"""
		rotate/shift/scale image using the gridding method.
		if kb = None, the image is prepared in this function before the rot/shift operation.
		If kb is NOT None, then the supplied (img,kb) must be the output from prepi3D.
		'mode' specifies what will be put in corners, should they stick out:
			background - leave original values
			cyclic - use pixels from the image using wrap-around transformation
		'wrap': option for using wraparound pixels during translations
	"""

	if scale == 0.0 :  global_def.ERROR("scale=0 not allowed", "rot_shift3D_grid", 1)

	if mode == "cyclic":
		pass#IMPORTIMPORTIMPORT from math import radians
		if kb == None:
			o, kb = prepi3D(img)
		else:
			o = img
		# gridding rotation/shift:
		#if  mirror: o.process_inplace("xform.mirror", {"axis":'x'})
		return o.rot_scale_conv_new_3D(numpy.radians(phi), numpy.radians(theta), numpy.radians(psi), sx, sy, sz, kb, scale, wrap)
	elif mode == "background":
		pass#IMPORTIMPORTIMPORT from math import radians
		if kb == None:
			o, kb = prepi3D(img)
		else:
			o = img
		# gridding rotation
		#if  mirror: o.process_inplace("xform.mirror", {"axis":'x'})
		return o.rot_scale_conv_new_background_3D(numpy.radians(phi), numpy.radians(theta), numpy.radians(psi), sx, sy, sz, kb, scale, wrap)	
	else: global_def.ERROR("rot_shift3D_grid mode not valid", "rot_shift3D_grid", 1)


def sinc2inv(nx):
	pass#IMPORTIMPORTIMPORT from math import sqrt
	s = sincinv(nx)
	return [i*i for i in s]

def sincinv(nx):
	pass#IMPORTIMPORTIMPORT from math import pi,sin
	cdf =numpy.pi/nx
	npad = 1
	nxb = nx/2/npad
	nxe = nxb + (nx/npad)%2
	s = [1.0]*nx
	for i in range( -nxb, nxe):
		if( i != 0 ):
			rrr=abs(i)
			s[i+nxb] = (rrr*cdf)/numpy.sin(rrr*cdf)
	return s

def welch_pw2(img, win_size=512, overlp_x=50, overlp_y=50, edge_x=0, edge_y=0):
	""" 
		Calculate the power spectrum using Welch periodograms (overlapped periodogram)
	"""
	pass#IMPORTIMPORTIMPORT from fundamentals import window2d, ramp
	pass#IMPORTIMPORTIMPORT from EMAN2 import periodogram
	nx = img.get_xsize()
	ny = img.get_ysize()
	nx_fft = smallprime(nx)
	ny_fft = smallprime(ny)
	x_gaussian_hi = 1./win_size
	pass#IMPORTIMPORTIMPORT from filter    import filt_gaussh
	e_fil = filter.filt_gaussh(window2d(img,nx_fft,ny_fft,"l"), x_gaussian_hi)
	x38 = 100/(100-overlp_x) # normalization of % of the overlap in x 
	x39 = 100/(100-overlp_y) # normalization of % of the overlap in y
	x26 = int(x38*((nx-2*edge_x)/win_size-1)+1)  # number of pieces horizontal dim.(X)
	x29 = int(x39*((ny-2*edge_y)/win_size-1)+1)  # number of pieces vertical dim.(Y)
	iz = 0	
	pw2 = EMAN2_cppwrap.EMData()
	for iy in range(1, x29+1):	
		x21 = (win_size/x39)*(iy-1) + edge_y  #  y-direction it should start from 0 if edge_y=0	      
		for ix in  range(1, x26+1):			 
			x22 = (win_size/x38)*(ix-1) + edge_x  # x-direction it should start from 0 if edge_x =0
			wi  = window2d(e_fil, win_size, win_size, "l", x22, x21)
			iz  = iz+1
			if (iz == 1): pw2  = EMAN2_cppwrap.periodogram(ramp(wi))
			else:         pw2 += EMAN2_cppwrap.periodogram(ramp(wi))
	return  pw2/float(iz)

def welch_pw2_tilt_band(img,theta,num_bnd=-1,overlp_y=50,edge_x=0,edge_y=0,win_s=256):
	""" 
		1. Calculate the power spectra of tilt bands
		2. The tilt micrograph is rotated such that the tilt axis is vertical (along Y axis)
		3. edge_x and edge_y are removed from the micrograph
	""" 
	pass#IMPORTIMPORTIMPORT from EMAN2 import periodogram
	nx = img.get_xsize()
	ny = img.get_ysize()
	num1 = int(nx-2*edge_x)
	num2 = int(ny-2*edge_y)
	nx_fft = smallprime(num1)
	ny_fft = smallprime(num2)
	img1 = window2d(img,nx_fft,ny_fft,"l",edge_x,edge_y)
	if(num_bnd == -1):
		num_bnd = int(nx_fft/win_s)
		win_x   = int(win_s)
	else:
		win_x = int(nx_fft/num_bnd)
		win_x = int(smallprime(win_x))
	win_y = win_x
	x_gaussian_hi = 1./win_x
	del img
	pass#IMPORTIMPORTIMPORT from filter import filt_gaussh
	pass#IMPORTIMPORTIMPORT from utilities import drop_image, rot_image
	# The input img is rotated such that tilt axis is vertical
	img2  = rot_image(img1,theta, 0, 0, 1.0,1.0)	
	e_fil = filter.filt_gaussh(img2, x_gaussian_hi)
	del img1
	del img2
	x39 = 100/(100-overlp_y) # normalization of % of the overlap in y
	x29 = int(x39*((ny)/win_y-1)+1)  # number of pieces vertical dim.(Y)
	pw2 = EMAN2_cppwrap.EMData()
	pw2_band = []
	for ix in  range(1, num_bnd+1):
		x22 = (win_x)*(ix-1)# x-direction it should start from 0 if edge_x =0
		iz=0
		for iy in range(1, x29+1):	
			x21 = (win_y/x39)*(iy-1) #  y-direction it should start from 0 if edge_y=0	      			 
			wi = window2d(e_fil,win_x, win_y,"l",x22, x21)
			iz = iz+1
			if (iz == 1): pw2  = EMAN2_cppwrap.periodogram(ramp(wi))
			else:         pw2 += EMAN2_cppwrap.periodogram(ramp(wi))
		pw2/=float(iz)
		# drop_image(pw2,"band%03d"%(ix))
		pw2_band.append(pw2)	
	return 	pw2_band


def tilemic(img, win_size=512, overlp_x=50, overlp_y=50, edge_x=0, edge_y=0):
	""" 
		Calculate set of periodograms for tiles.  Returns a list.
	"""
	pass#IMPORTIMPORTIMPORT from fundamentals import window2d, ramp
	pass#IMPORTIMPORTIMPORT from EMAN2 import periodogram
	nx = img.get_xsize()
	ny = img.get_ysize()
	nx_fft = smallprime(nx)
	ny_fft = smallprime(ny)
	x_gaussian_hi = 1./win_size
	pass#IMPORTIMPORTIMPORT from filter    import filt_gaussh
	e_fil = filter.filt_gaussh(window2d(img,nx_fft,ny_fft,"l"), x_gaussian_hi)
	x38 = 100/(100-overlp_x) # normalization of % of the overlap in x 
	x39 = 100/(100-overlp_y) # normalization of % of the overlap in y
	x26 = int(x38*((nx-2*edge_x)/win_size-1)+1)  # number of pieces horizontal dim.(X)
	x29 = int(x39*((ny-2*edge_y)/win_size-1)+1)  # number of pieces vertical dim.(Y)
	pw2 = []
	for iy in range(1, x29+1):	
		x21 = (win_size/x39)*(iy-1) + edge_y  #  y-direction it should start from 0 if edge_y=0	      
		for ix in  range(1, x26+1):			 
			x22 = (win_size/x38)*(ix-1) + edge_x  # x-direction it should start from 0 if edge_x =0
			wi  = ramp( window2d(e_fil, win_size, win_size, "l", x22, x21) )
			st = EMAN2_cppwrap.Util.infomask(wi, None, True)
			wi = (wi - st[0])/st[1]*win_size
			pw2.append(EMAN2_cppwrap.periodogram(wi))
	return  pw2


def bracket(f,x1,h):
	c = 1.618033989 
	f1 = f(x1)
	x2 = x1 + h; f2 = f(x2)
	# Determine downhill direction and change sign of h if needed
	if f2 > f1:
		h = -h
		x2 = x1 + h; f2 = f(x2)
		# Check if minimum between x1 - h and x1 + h
		if f2 > f1: return x2,x1 - h 
	# Search loop
	for i in range (100):    
		h = c*h
		x3 = x2 + h; f3 = f(x3)
		if f3 > f2: return x1,x3
		x1 = x2; x2 = x3
		f1 = f2; f2 = f3
	print("Bracket did not find a mimimum")        
 
