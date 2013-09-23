from EMAN2_cppwrap import *
from global_def import *


iteration_count = -1


def ali3d_reference1( ref_data ):
	from utilities      import print_msg
	from filter         import fit_tanh1, filt_tanl
	from fundamentals   import fshift
	from morphology     import threshold
	#  Prepare the reference in 3D alignment, i.e., low-pass filter and center.
	#  Input: list ref_data
	#   0 - mask
	#   1 - center flag
	#   2 - raw average
	#   3 - fsc result
	#  Output: filtered, centered, and masked reference image
	#  apply filtration (FSC) to reference image:

	print_msg("reference1\n")
	cs = [0.0]*3

	stat = Util.infomask(ref_data[2], ref_data[0], False)
	volf = ref_data[2] - stat[0]
	Util.mul_scalar(volf, 1.0/stat[1])
	volf = threshold(volf)
	Util.mul_img(volf, ref_data[0])
	#fl, aa = fit_tanh1(ref_data[3], 0.1)
	fl = 0.15
	aa = 0.2
	msg = "Tangent filter:  cut-off frequency = %10.3f        fall-off = %10.3f\n"%(fl, aa)
	print_msg(msg)
	volf = filt_tanl(volf, fl, aa)
	if ref_data[1] == 1:
		cs = volf.phase_cog()
		msg = "Center x = %10.3f        Center y = %10.3f        Center z = %10.3f\n"%(cs[0], cs[1], cs[2])
		print_msg(msg)
		volf  = fshift(volf, -cs[0], -cs[1], -cs[2])
	return  volf, cs



def ali3d_reference2( ref_data ):
	from utilities      import print_msg, get_im
	from filter         import fit_tanh1, filt_tanl, filt_table
	from fundamentals   import fshift, rops_table, fft
	from morphology     import threshold
	#  Prepare the reference in 3D alignment, i.e., low-pass filter and center.
	#  Input: list ref_data
	#   0 - mask
	#   1 - center flag
	#   2 - raw average
	#   3 - fsc result
	#  Output: filtered, centered, and masked reference image
	#  apply filtration (FSC) to reference image:

	print_msg("reference2\n")
	cs = [0.0]*3

	stat = Util.infomask(ref_data[2], ref_data[0], False)
	volf = ref_data[2] - stat[0]
	Util.mul_scalar(volf, 1.0/stat[1])
	volf = threshold(volf)
	Util.mul_img(volf, ref_data[0])
	o = rops_table(get_im("../model_structure.hdf"))
	n = rops_table(volf)
	#from math import sqrt
	ff = [0.0]*len(n)
	for i in xrange(len(o)):  ff[i] = o[i]/n[i]
	#fl, aa = fit_tanh1(ref_data[3], 0.1)
	fl = 0.20
	aa = 0.2
	msg = "Tangent filter:  cut-off frequency = %10.3f        fall-off = %10.3f\n"%(fl, aa)
	print_msg(msg)
	volf = fft( filt_tanl( filt_table(fft(volf),ff), fl, aa) )
	if ref_data[1] == 1:
		cs = volf.phase_cog()
		msg = "Center x = %10.3f        Center y = %10.3f        Center z = %10.3f\n"%(cs[0], cs[1], cs[2])
		print_msg(msg)
		volf  = fshift(volf, -cs[0], -cs[1], -cs[2])
	return  volf, cs



def ali3d_reference3( ref_data ):
	from utilities      import print_msg
	from filter         import fit_tanh1, filt_tanl
	from fundamentals   import fshift
	from morphology     import threshold
	#  Prepare the reference in 3D alignment, i.e., low-pass filter and center.
	#  Input: list ref_data
	#   0 - mask
	#   1 - center flag
	#   2 - raw average
	#   3 - fsc result
	#  Output: filtered, centered, and masked reference image
	#  apply filtration (FSC) to reference image (commented out).
	global  iteration_count
	iteration_count += 1
	print_msg("User function reference3, iteration #%3d\n"%(iteration_count+1) )
	cs = [0.0]*3

	stat = Util.infomask(ref_data[2], ref_data[0], False)
	volf = ref_data[2] - stat[0]
	Util.mul_scalar(volf, 1.0/stat[1])
	volf = threshold(volf)
	Util.mul_img(volf, ref_data[0])
	#fl, aa = fit_tanh1(ref_data[3], 0.1)
	fl = min(0.35, 0.15 + 0.03*(iteration_count//2))
	aa = 0.2
	msg = "Tangent filter:  cut-off frequency = %10.3f        fall-off = %10.3f\n"%(fl, aa)
	print_msg(msg)
	volf = filt_tanl(volf, fl, aa)
	if ref_data[1] == 1:
		cs = volf.phase_cog()
		msg = "Center x = %10.3f        Center y = %10.3f        Center z = %10.3f\n"%(cs[0], cs[1], cs[2])
		print_msg(msg)
		volf  = fshift(volf, -cs[0], -cs[1], -cs[2])
	return  volf, cs

def localali3d_reference3( ref_data ):
	from utilities      import print_msg
	from filter         import fit_tanh1, filt_tanl
	from fundamentals   import fshift
	from morphology     import threshold
	#  Prepare the reference in 3D alignment, i.e., low-pass filter and center.
	#  Input: list ref_data
	#   0 - mask
	#   1 - center flag
	#   2 - raw average
	#   3 - fsc result
	#  Output: filtered, centered, and masked reference image
	#  apply filtration (FSC) to reference image:

	print_msg("User function localali3d_reference3\n")
	cs = [0.0]*3

	stat = Util.infomask(ref_data[2], ref_data[0], False)
	volf = ref_data[2] - stat[0]
	Util.mul_scalar(volf, 1.0/stat[1])
	volf = threshold(volf)
	Util.mul_img(volf, ref_data[0])
	#fl, aa = fit_tanh1(ref_data[3], 0.1)
	fl = 0.4
	aa = 0.2
	msg = "Tangent filter:  cut-off frequency = %10.3f        fall-off = %10.3f\n"%(fl, aa)
	print_msg(msg)
	volf = filt_tanl(volf, fl, aa)
	if ref_data[1] == 1:
		cs = volf.phase_cog()
		msg = "Center x = %10.3f        Center y = %10.3f        Center z = %10.3f\n"%(cs[0], cs[1], cs[2])
		print_msg(msg)
		volf  = fshift(volf, -cs[0], -cs[1], -cs[2])
	return  volf, cs	
	

def localali3d_reference4( ref_data ):
	from utilities      import print_msg
	from filter         import fit_tanh, filt_tanl
	from fundamentals   import fshift
	from morphology     import threshold
	#  Prepare the reference in 3D alignment, i.e., low-pass filter and center.
	#  Input: list ref_data
	#   0 - mask
	#   1 - center flag
	#   2 - raw average
	#   3 - fsc result
	#  Output: filtered, centered, and masked reference image
	#  apply filtration (FSC) to reference image:

	print_msg("reference4\n")
	cs = [0.0]*3

	stat = Util.infomask(ref_data[2], ref_data[0], False)
	volf = ref_data[2] - stat[0]
	Util.mul_scalar(volf, 1.0/stat[1])
	#volf = threshold(volf)
	#Util.mul_img(volf, ref_data[0])
	#fl, aa = fit_tanh(ref_data[3])
	fl = 0.40
	aa = 0.1
	msg = "Tangent filter:  cut-off frequency = %10.3f        fall-off = %10.3f\n"%(fl, aa)
	print_msg(msg)
	volf = filt_tanl(volf, fl, aa)
	if ref_data[1] == 1:
		cs = volf.phase_cog()
		msg = "Center x = %10.3f        Center y = %10.3f        Center z = %10.3f\n"%(cs[0], cs[1], cs[2])
		print_msg(msg)
		volf  = fshift(volf, -cs[0], -cs[1], -cs[2])
	return  volf, cs
	
