/*
 * Author: Steven Ludtke, 04/10/2003 (sludtke@bcm.edu)
 * Copyright (c) 2000-2006 Baylor College of Medicine
 *
 * This software is issued under a joint BSD/GNU license. You may use the
 * source code in this file under either license. However, note that the
 * complete EMAN2 and SPARX software packages have some GPL dependencies,
 * so you are responsible for compliance with the licenses of these packages
 * if you opt to use BSD licensing. The warranty disclaimer below holds
 * in either instance.
 *
 * This complete copyright notice must be included in any revised version of the
 * source code. Additional authorship citations may be added, but existing
 * author citations must be preserved.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 * */

#include "cmp.h"
#include "emdata.h"
#include "ctf.h"
#include "plugins/cmp_template.h"
#undef max
#include <climits>

#ifdef EMAN2_USING_CUDA
// Only CCC, DOT  and CCC.TOMO are cuda enabled
#include "cuda/cuda_processor.h"
#include "cuda/cuda_cmp.h"
#endif // EMAN2_USING_CUDA

using namespace EMAN;

const string CccCmp::NAME = "ccc";
const string LodCmp::NAME = "lod";
const string SqEuclideanCmp::NAME = "sqeuclidean";
const string DotCmp::NAME = "dot";
const string TomoCccCmp::NAME = "ccc.tomo";
const string TomoWedgeCccCmp::NAME = "ccc.tomo.thresh";
const string TomoWedgeFscCmp::NAME = "fsc.tomo.auto";
const string TomoFscCmp::NAME = "fsc.tomo";
const string QuadMinDotCmp::NAME = "quadmindot";
const string OptVarianceCmp::NAME = "optvariance";
const string OptSubCmp::NAME = "optsub";
const string PhaseCmp::NAME = "phase";
const string FRCCmp::NAME = "frc";

template <> Factory < Cmp >::Factory()
{
	force_add<CccCmp>();
	force_add<LodCmp>();
	force_add<SqEuclideanCmp>();
	force_add<DotCmp>();
	force_add<TomoCccCmp>();
	force_add<TomoWedgeCccCmp>();
	force_add<TomoWedgeFscCmp>();
	force_add<TomoFscCmp>();
	force_add<QuadMinDotCmp>();
	force_add<OptVarianceCmp>();
	force_add<OptSubCmp>();
	force_add<PhaseCmp>();
	force_add<FRCCmp>();
//	force_add<XYZCmp>();
}

void Cmp::validate_input_args(const EMData * image, const EMData *with) const
{
	
#ifdef EMAN2_USING_CUDA  
	if (image->getcudarwdata() && with->getcudarwdata()) {
		//no need for futher checking, which will induce an expensive copy from device to host
		return;
	}
#endif
	if (!image) {
		throw NullPointerException("compared image");
	}
	if (!with) {
		throw NullPointerException("compare-with image");
	}

	if (!EMUtil::is_same_size(image, with)) {
		throw ImageFormatException( "images not same size");
	}

	float *d1 = image->get_data();
	if (!d1) {
		throw NullPointerException("image contains no data");
	}

	float *d2 = with->get_data();
	if (!d2) {
		throw NullPointerException("compare-with image data");
	}
}

//  It would be good to add code for complex images!  PAP
float CccCmp::cmp(EMData * image, EMData *with) const
{
	ENTERFUNC;
	if (image->is_complex() || with->is_complex())
		throw ImageFormatException( "Complex images not supported by CMP::CccCmp");
	validate_input_args(image, with);

	const float *const d1 = image->get_const_data();
	const float *const d2 = with->get_const_data();

	float negative = (float)params.set_default("negative", 1);
	if (negative) negative=-1.0; else negative=1.0;

	double avg1 = 0.0, var1 = 0.0, avg2 = 0.0, var2 = 0.0, ccc = 0.0;
	long n = 0;
	size_t totsize = image->get_xsize()*image->get_ysize()*image->get_zsize();

	bool has_mask = false;
	EMData* mask = 0;
	if (params.has_key("mask")) {
		mask = params["mask"];
		if(mask!=0) {has_mask=true;}
	}
#ifdef EMAN2_USING_CUDA
	if (image->getcudarwdata() && with->getcudarwdata()) {
		//cout << "CUDA ccc cmp" << endl;
		float* maskdata = 0;
		if(has_mask && !mask->getcudarwdata()){
			mask->copy_to_cuda();
			maskdata = mask->getcudarwdata();
		}
		float ccc = ccc_cmp_cuda(image->getcudarwdata(), with->getcudarwdata(), maskdata, image->get_xsize(), image->get_ysize(), image->get_zsize());
		ccc *= negative;
		//cout << "CUDA CCC is: " << ccc << endl;
		return ccc;
	}
#endif
	if (has_mask) {
		const float *const dm = mask->get_const_data();
		for (size_t i = 0; i < totsize; ++i) {
			if (dm[i] > 0.5) {
				avg1 += double(d1[i]);
				var1 += d1[i]*double(d1[i]);
				avg2 += double(d2[i]);
				var2 += d2[i]*double(d2[i]);
				ccc += d1[i]*double(d2[i]);
				n++;
			}
		}
	} else {
		for (size_t i = 0; i < totsize; ++i) {
			avg1 += double(d1[i]);
			var1 += d1[i]*double(d1[i]);
			avg2 += double(d2[i]);
			var2 += d2[i]*double(d2[i]);
			ccc += d1[i]*double(d2[i]);
		}
		n = totsize;
	}

	avg1 /= double(n);
	var1 = var1/double(n) - avg1*avg1;
	avg2 /= double(n);
	var2 = var2/double(n) - avg2*avg2;
	ccc = ccc/double(n) - avg1*avg2;
	ccc /= sqrt(var1*var2);
	if (!Util::goodf(&ccc)) ccc=-2.0;		// Steve - if one image was 0, this returned nan, which messes certain things up. -2.0 is out of range, and should serve as a proxy
	ccc *= negative;
	return static_cast<float>(ccc);
	EXITFUNC;
}


// Added by JAK 11/12/10
// L^1-norm difference of two maps, after normalization.
float LodCmp::cmp(EMData * image, EMData *with) const
{
	ENTERFUNC;
	if (image->is_complex() || with->is_complex())
		throw ImageFormatException( "Complex images not yet supported by CMP::LodCmp");
	validate_input_args(image, with);

	const float *const d1 = image->get_const_data();
	const float *const d2 = with->get_const_data();

	float negative = (float)params.set_default("negative", 1);
	if (negative) negative=-1.0; else negative=1.0;

	double avg1 = 0.0, sig1 = 0.0, avg2 = 0.0, sig2 = 0.0, lod = 0.0;
	long n = 0;
	size_t totsize = image->get_xsize()*image->get_ysize()*image->get_zsize();
	size_t i;
	
	bool has_mask = false;
	EMData* mask = 0;
	if (params.has_key("mask")) {
		mask = params["mask"];
		if(mask!=0) {has_mask=true;}
	}

	int normalize = 0;
	if (params.has_key("normalize")) {
		normalize = params["normalize"];
	}

	if (normalize) {
		if (has_mask) {
			const float *const dm = mask->get_const_data();
			for (i = 0; i < totsize; ++i) {
				if (dm[i] > 0.5) {
					avg1 += double(d1[i]);
					avg2 += double(d2[i]);
					n++;
				}
			}
		} else {
			for (i = 0; i < totsize; ++i) {
				avg1 += double(d1[i]);
				avg2 += double(d2[i]);
			}
			n = totsize;
		}

		avg1 /= double(n);
		avg2 /= double(n);

		if (has_mask) {
			const float *const dm = mask->get_const_data();
			for (i = 0; i < totsize; ++i) {
				if (dm[i] > 0.5) {
					sig1 += fabs(double(d1[i])-avg1);
					sig2 += fabs(double(d2[i])-avg2);
				}
			}
		} else {
			for (i = 0; i < totsize; ++i) {
				sig1 += fabs(double(d1[i])-avg1);
				sig2 += fabs(double(d2[i])-avg2);
			}
		}
	} else {
		avg1 = 0.0; avg2 = 0.0;
		sig1 = 1.0; sig2 = 1.0;
	}

	if (has_mask) {
		const float *const dm = mask->get_const_data();
		for (i = 0; i < totsize; ++i) {
			if (dm[i] > 0.5) {
				lod += fabs((double(d1[i])-avg1)/sig1 - (double(d2[i])-avg2)/sig2);
			}
		}
	} else {
		for (i = 0; i < totsize; ++i) {
			lod += fabs((double(d1[i])-avg1)/sig1 - (double(d2[i])-avg2)/sig2);
		}
	}
	
	lod *= (-0.5);
	lod *= negative;
	return static_cast<float>(lod);
	EXITFUNC;
}


//float SqEuclideanCmp::cmp(EMData * image, EMData *withorig) const
float SqEuclideanCmp::cmp(EMData *image,EMData * withorig ) const
{
	ENTERFUNC;
	EMData *with = withorig;
	validate_input_args(image, with);

	int zeromask = params.set_default("zeromask",0);
	int normto = params.set_default("normto",0);

	if (normto) {
		if (zeromask) with = withorig->process("normalize.toimage",Dict("to",image));
		else with = withorig->process("normalize.toimage",Dict("to",image,"ignore_zero",0));
		with->set_attr("deleteme",1);
		if ((float)(with->get_attr("norm_mult"))<=0) {		// This means the normalization inverted the image, a clear probablity of noise bias, so we undo the normalization
			delete with;
			with=withorig;
		}
	}

	const float *const y_data = with->get_const_data();
	const float *const x_data = image->get_const_data();
	double result = 0.;
	float n = 0.0f;
	if(image->is_complex() && with->is_complex()) {
	// Implemented by PAP  01/09/06 - please do not change.  If in doubts, write/call me.
		int nx  = with->get_xsize();
		int ny  = with->get_ysize();
		int nz  = with->get_zsize();
		nx = (nx - 2 + with->is_fftodd()); // nx is the real-space size of the input image
		int lsd2 = (nx + 2 - nx%2) ; // Extended x-dimension of the complex image

		int ixb = 2*((nx+1)%2);
		int iyb = ny%2;
		//
		if(nz == 1) {
		//  it looks like it could work in 3D, but it does not, really.
		for ( int iz = 0; iz <= nz-1; iz++) {
			double part = 0.;
			for ( int iy = 0; iy <= ny-1; iy++) {
				for ( int ix = 2; ix <= lsd2 - 1 - ixb; ix++) {
						size_t ii = ix + (iy  + iz * ny)* lsd2;
						part += (x_data[ii] - y_data[ii])*double(x_data[ii] - y_data[ii]);
				}
			}
			for ( int iy = 1; iy <= ny/2-1 + iyb; iy++) {
				size_t ii = (iy  + iz * ny)* lsd2;
				part += (x_data[ii] - y_data[ii])*double(x_data[ii] - y_data[ii]);
				part += (x_data[ii+1] - y_data[ii+1])*double(x_data[ii+1] - y_data[ii+1]);
			}
			if(nx%2 == 0) {
				for ( int iy = 1; iy <= ny/2-1 + iyb; iy++) {
					size_t ii = lsd2 - 2 + (iy  + iz * ny)* lsd2;
					part += (x_data[ii] - y_data[ii])*double(x_data[ii] - y_data[ii]);
					part += (x_data[ii+1] - y_data[ii+1])*double(x_data[ii+1] - y_data[ii+1]);
				}

			}
			part *= 2;
			part += (x_data[0] - y_data[0])*double(x_data[0] - y_data[0]);
			if(ny%2 == 0) {
				int ii = (ny/2  + iz * ny)* lsd2;
				part += (x_data[ii] - y_data[ii])*double(x_data[ii] - y_data[ii]);
			}
			if(nx%2 == 0) {
				int ii = lsd2 - 2 + (0  + iz * ny)* lsd2;
				part += (x_data[ii] - y_data[ii])*double(x_data[ii] - y_data[ii]);
				if(ny%2 == 0) {
					int ii = lsd2 - 2 +(ny/2  + iz * ny)* lsd2;
					part += (x_data[ii] - y_data[ii])*double(x_data[ii] - y_data[ii]);
				}
			}
			result += part;
		}
		n = (float)nx*(float)ny*(float)nz*(float)nx*(float)ny*(float)nz;

		} else { //This 3D code is incorrect, but it is the best I can do now 01/09/06 PAP
		int ky, kz;
		int ny2 = ny/2; int nz2 = nz/2;
		for ( int iz = 0; iz <= nz-1; iz++) {
			if(iz>nz2) kz=iz-nz; else kz=iz;
			for ( int iy = 0; iy <= ny-1; iy++) {
				if(iy>ny2) ky=iy-ny; else ky=iy;
				for ( int ix = 0; ix <= lsd2-1; ix++) {
				// Skip Friedel related values
				if(ix>0 || (kz>=0 && (ky>=0 || kz!=0))) {
						size_t ii = ix + (iy  + iz * ny)* (size_t)lsd2;
						result += (x_data[ii] - y_data[ii])*double(x_data[ii] - y_data[ii]);
					}
				}
			}
		}
		n = ((float)nx*(float)ny*(float)nz*(float)nx*(float)ny*(float)nz)/2.0f;
		}
	} else {		// real space
		size_t totsize = (size_t)image->get_xsize()*image->get_ysize()*image->get_zsize();
		if (params.has_key("mask")) {
		  EMData* mask;
		  mask = params["mask"];
  		  const float *const dm = mask->get_const_data();
		  for (size_t i = 0; i < totsize; i++) {
			   if (dm[i] > 0.5) {
				double temp = x_data[i]- y_data[i];
				result += temp*temp;
				n++;
			   }
		  }
		} 
		else if (zeromask) {
			n=0;
			for (size_t i = 0; i < totsize; i++) {
				if (x_data[i]==0 || y_data[i]==0) continue;
				double temp = x_data[i]- y_data[i];
				result += temp*temp;
				n++;
			}
			
		}
		else {
		  for (size_t i = 0; i < totsize; i++) {
				double temp = x_data[i]- y_data[i];
				result += temp*temp;
		   }
		   n = (float)totsize;
		}
	}
	result/=n;

	EXITFUNC;
	if (with->has_attr("deleteme")) delete with;
	float ret = (float)result;
	if (!Util::goodf(&ret)) return FLT_MAX;
	return ret;
}


// Even though this uses doubles, it might be wise to recode it row-wise
// to avoid numerical errors on large images
float DotCmp::cmp(EMData* image, EMData* with) const
{
	ENTERFUNC;
	
	validate_input_args(image, with);

	int normalize = params.set_default("normalize", 0);
	float negative = (float)params.set_default("negative", 1);
	if (negative) negative=-1.0; else negative=1.0;
#ifdef EMAN2_USING_CUDA // SO far only works for real images I put CUDA first to avoid running non CUDA overhead (calls to getdata are expensive!!!!)
	if(image->is_complex() && with->is_complex()) {
	} else {
		if (image->getcudarwdata() && with->getcudarwdata()) {
			//cout << "CUDA dot cmp" << endl;
			float* maskdata = 0;
			bool has_mask = false;
			EMData* mask = 0;
			if (params.has_key("mask")) {
				mask = params["mask"];
				if(mask!=0) {has_mask=true;}
			}
			if(has_mask && !mask->getcudarwdata()){
				mask->copy_to_cuda();
				maskdata = mask->getcudarwdata();
			}

			float result = dot_cmp_cuda(image->getcudarwdata(), with->getcudarwdata(), maskdata, image->get_xsize(), image->get_ysize(), image->get_zsize());
			result *= negative;

			return result;
			
		}
	}
#endif
	const float *const x_data = image->get_const_data();
	const float *const y_data = with->get_const_data();

	double result = 0.;
	long n = 0;
	if(image->is_complex() && with->is_complex()) {
	// Implemented by PAP  01/09/06 - please do not change.  If in doubts, write/call me.
		int nx  = with->get_xsize();
		int ny  = with->get_ysize();
		int nz  = with->get_zsize();
		nx = (nx - 2 + with->is_fftodd()); // nx is the real-space size of the input image
		int lsd2 = (nx + 2 - nx%2) ; // Extended x-dimension of the complex image

		int ixb = 2*((nx+1)%2);
		int iyb = ny%2;
		//
		if(nz == 1) {
		//  it looks like it could work in 3D, but does not
		for ( int iz = 0; iz <= nz-1; ++iz) {
			double part = 0.;
			for ( int iy = 0; iy <= ny-1; ++iy) {
				for ( int ix = 2; ix <= lsd2 - 1 - ixb; ++ix) {
					size_t ii = ix + (iy  + iz * ny)* lsd2;
					part += x_data[ii] * double(y_data[ii]);
				}
			}
			for ( int iy = 1; iy <= ny/2-1 + iyb; ++iy) {
				size_t ii = (iy  + iz * ny)* lsd2;
				part += x_data[ii] * double(y_data[ii]);
				part += x_data[ii+1] * double(y_data[ii+1]);
			}
			if(nx%2 == 0) {
				for ( int iy = 1; iy <= ny/2-1 + iyb; ++iy) {
					size_t ii = lsd2 - 2 + (iy  + iz * ny)* lsd2;
					part += x_data[ii] * double(y_data[ii]);
					part += x_data[ii+1] * double(y_data[ii+1]);
				}

			}
			part *= 2;
			part += x_data[0] * double(y_data[0]);
			if(ny%2 == 0) {
				size_t ii = (ny/2  + iz * ny)* lsd2;
				part += x_data[ii] * double(y_data[ii]);
			}
			if(nx%2 == 0) {
				size_t ii = lsd2 - 2 + (0  + iz * ny)* lsd2;
				part += x_data[ii] * double(y_data[ii]);
				if(ny%2 == 0) {
					int ii = lsd2 - 2 +(ny/2  + iz * ny)* lsd2;
					part += x_data[ii] * double(y_data[ii]);
				}
			}
			result += part;
		}
		if( normalize ) {
		//  it looks like it could work in 3D, but does not
		double square_sum1 = 0., square_sum2 = 0.;
		for ( int iz = 0; iz <= nz-1; ++iz) {
			for ( int iy = 0; iy <= ny-1; ++iy) {
				for ( int ix = 2; ix <= lsd2 - 1 - ixb; ++ix) {
					size_t ii = ix + (iy  + iz * ny)* lsd2;
					square_sum1 += x_data[ii] * double(x_data[ii]);
					square_sum2 += y_data[ii] * double(y_data[ii]);
				}
			}
			for ( int iy = 1; iy <= ny/2-1 + iyb; ++iy) {
				size_t ii = (iy  + iz * ny)* lsd2;
				square_sum1 += x_data[ii] * double(x_data[ii]);
				square_sum1 += x_data[ii+1] * double(x_data[ii+1]);
				square_sum2 += y_data[ii] * double(y_data[ii]);
				square_sum2 += y_data[ii+1] * double(y_data[ii+1]);
			}
			if(nx%2 == 0) {
				for ( int iy = 1; iy <= ny/2-1 + iyb; ++iy) {
					size_t ii = lsd2 - 2 + (iy  + iz * ny)* lsd2;
					square_sum1 += x_data[ii] * double(x_data[ii]);
					square_sum1 += x_data[ii+1] * double(x_data[ii+1]);
					square_sum2 += y_data[ii] * double(y_data[ii]);
					square_sum2 += y_data[ii+1] * double(y_data[ii+1]);
				}

			}
			square_sum1 *= 2;
			square_sum1 += x_data[0] * double(x_data[0]);
			square_sum2 *= 2;
			square_sum2 += y_data[0] * double(y_data[0]);
			if(ny%2 == 0) {
				int ii = (ny/2  + iz * ny)* lsd2;
				square_sum1 += x_data[ii] * double(x_data[ii]);
				square_sum2 += y_data[ii] * double(y_data[ii]);
			}
			if(nx%2 == 0) {
				int ii = lsd2 - 2 + (0  + iz * ny)* lsd2;
				square_sum1 += x_data[ii] * double(x_data[ii]);
				square_sum2 += y_data[ii] * double(y_data[ii]);
				if(ny%2 == 0) {
					int ii = lsd2 - 2 +(ny/2  + iz * ny)* lsd2;
					square_sum1 += x_data[ii] * double(x_data[ii]);
					square_sum2 += y_data[ii] * double(y_data[ii]);
				}
			}
		}
		result /= sqrt(square_sum1*square_sum2);
		} 
			else  result /= ((float)nx*(float)ny*(float)nz*(float)nx*(float)ny*(float)nz);
			// This concludes the 2D part.
		} else {
		// The 3D code
		int ky, kz;
		int ny2 = ny/2; int nz2 = nz/2;
		for ( int iz = 0; iz <= nz-1; ++iz) {
			if(iz>nz2) kz=iz-nz; else kz=iz;
			for ( int iy = 0; iy <= ny-1; ++iy) {
				if(iy>ny2) ky=iy-ny; else ky=iy;
				for (int ix = 0; ix <= lsd2-1; ++ix) {
						int p = 2;
						if  (ix <= 1 || ix >= lsd2-2 && nx%2 == 0)   p = 1;
						size_t ii = ix + (iy  + iz * ny)* (size_t)lsd2;
						result += p * x_data[ii] * double(y_data[ii]);
				}
			}
		}
		if( normalize ) {
		double square_sum1 = 0., square_sum2 = 0.;
		int ky, kz;
		int ny2 = ny/2; int nz2 = nz/2;
		for ( int iz = 0; iz <= nz-1; ++iz) {
			if(iz>nz2) kz=iz-nz; else kz=iz;
			for ( int iy = 0; iy <= ny-1; ++iy) {
				if(iy>ny2) ky=iy-ny; else ky=iy;
				for (int ix = 0; ix <= lsd2-1; ++ix) {
						int p = 2;
						if  (ix <= 1 || ix >= lsd2-2 && nx%2 == 0)   p = 1;
						size_t ii = ix + (iy  + iz * ny)* (size_t)lsd2;
						square_sum1 += p * x_data[ii] * double(x_data[ii]);
						square_sum2 += p * y_data[ii] * double(y_data[ii]);
				}
			}
		}
		result /= sqrt(square_sum1*square_sum2);
		} 
			else result /= ((float)nx*(float)ny*(float)nz*(float)nx*(float)ny*(float)nz);
		}
	} else {
		// This part is for when two images are real, which is much easier because 2-D or 3-D
		// is the same code.
		size_t totsize = (size_t)image->get_xsize() * image->get_ysize() * image->get_zsize();

		double square_sum1 = 0., square_sum2 = 0.;

		if (params.has_key("mask")) {
			EMData* mask;
			mask = params["mask"];
			const float *const dm = mask->get_const_data();
			if (normalize) {
				for (size_t i = 0; i < totsize; i++) {
					if (dm[i] > 0.5) {
						square_sum1 += x_data[i]*double(x_data[i]);
						square_sum2 += y_data[i]*double(y_data[i]);
						result += x_data[i]*double(y_data[i]);
					}
				}
			} else {
				for (size_t i = 0; i < totsize; i++) {
					if (dm[i] > 0.5) {
						result += x_data[i]*double(y_data[i]);
						n++;
					}
				}
			}
		} else {
			// this little bit of manual loop unrolling makes the dot product as fast as sqeuclidean with -O2
			for (size_t i=0; i<totsize; i++) result+=x_data[i]*y_data[i];

			if (normalize) {
				square_sum1 = image->get_attr("square_sum");
				square_sum2 = with->get_attr("square_sum");
			} else n = totsize;
		}
		if (normalize) result /= (sqrt(square_sum1*square_sum2)); else result /= n;
	}


	EXITFUNC;
	return (float) (negative*result);
}

// This implements the technique of Mike Schmid where by the cross correlation is normalized
// in an effort to remove the effects of the missing wedge. Somewhat of a heuristic solution, but it seems
// to work. Basically it relies on the observation that 'good' matchs will conentrate the correlation
// signal in the peak, wheras 'bad' correlations will distribute the signal.
// John Flanagan 18-10-2010

float TomoCccCmp::cmp(EMData * image, EMData *with) const
{
	ENTERFUNC;
	EMData* ccf = params.set_default("ccf",(EMData*) NULL);
	bool ccf_ownership = false;
	bool norm = params.set_default("norm",true);
	float negative = (float)params.set_default("negative", 1);
	if (negative) negative=-1.0; else negative=1.0;
	
#ifdef EMAN2_USING_CUDA	
	if(image->getcudarwdata() && with->getcudarwdata()){
		if (!ccf) {
			ccf = image->calc_ccf(with);
			ccf_ownership = true;
		}
		//cout << "using CUDA" << endl;
		float2 stats = get_stats_cuda(ccf->getcudarwdata(), ccf->get_xsize(), ccf->get_ysize(), ccf->get_zsize());
		float best_score = get_value_at_wrap_cuda(ccf->getcudarwdata(), 0, 0, 0, ccf->get_xsize(), ccf->get_ysize(), ccf->get_zsize());
		if(norm) {
			best_score = negative*(best_score - stats.x)/sqrt(stats.y);
		} else {
			best_score = negative*best_score;
		}
		
		if (ccf_ownership) delete ccf; ccf = 0;
		
		EXITFUNC;
		return best_score;
		
	}
#endif

	if (!ccf) {
		ccf = image->calc_ccf(with);
		ccf_ownership = true;
	}
	if (norm) ccf->process_inplace("normalize");
	
	float best_score = ccf->get_value_at_wrap(0,0,0);
        if (ccf_ownership) delete ccf; ccf = 0;

	return negative*best_score;

}

// CCC in Fourier space with small valued pixels presumed to be in missing wedge, and excluded from similarity
// Fractional overlap returned in "image"

float TomoWedgeCccCmp::cmp(EMData * image, EMData *with) const
{
	ENTERFUNC;
	if (!image->is_complex() || !with->is_complex())  throw InvalidCallException("Error: TomoWedgeCccCmp requires complex images");
	if (image->get_xsize()!=with->get_xsize() || image->get_ysize()!=with->get_ysize() || image->get_zsize()!=with->get_zsize())  throw InvalidCallException("Error: TomoWedgeCccCmp requires complex images");

	float sigmaimg = params.set_default("sigmaimg",0.5f);
	float sigmawith = params.set_default("sigmawith",0.5f);
	int negative = params.set_default("negative",1);

	// Note 'sigma' is the STD of real and imag values treated independently
	// s1 and s2 are threshold values squared (for speed)
	float s1=pow((float)image->get_attr("sigma")*sigmaimg,(float)2.0);		
	float s2=pow((float)with->get_attr("sigma")*sigmawith,(float)2.0);
	
	
	double sum=0;
	double sumsq1=0;
	double sumsq2=0;
	double norm=0;
	for (int z=0; z<image->get_zsize(); z++) {
		for (int y=0; y<image->get_ysize(); y++) {
			for (int x=0; x<image->get_xsize(); x+=2) {
				float v1r=image->get_value_at(x,y,z);
				float v1i=image->get_value_at(x+1,y,z);
				float v1=Util::square_sum(v1r,v1i);
				if (v1<s1) continue;
				
				float v2r=with->get_value_at(x,y,z);
				float v2i=with->get_value_at(x+1,y,z);
				float v2=Util::square_sum(v2r,v2i);
				if (v2<s2) continue;
				
				sum+=v1r*v2r+v1i*v2i;
				sumsq1+=v1;
				if (Util::is_nan(sumsq1)) { printf("TomoWedgeCccCmp: NaN encountered: %d %d %d %f %f %f\n",x,y,z,v1r,v1i,v1); }
				sumsq2+=v2;
				norm+=1.0;
			}
		}
	}
	image->set_attr("fft_overlap",(float)(2.0*norm/(image->get_xsize()*image->get_ysize()*image->get_zsize())));
//	printf("%f\t%f\t%f\t%f\t%f\n",s1,s2,sumsq1,sumsq2,norm);
	
	if (negative) sum*=-1.0;
	return float(sum/(sqrt(sumsq1)*sqrt(sumsq2)));
}

// CCC in Fourier space with small valued pixels presumed to be in missing wedge, and excluded from similarity
// Fractional overlap returned in "image"

float TomoWedgeFscCmp::cmp(EMData * image, EMData *with) const
{
	ENTERFUNC;
	if (!image->is_complex() || !with->is_complex())  throw InvalidCallException("Error: TomoWedgeFscCmp requires complex images");
	if (image->get_xsize()!=with->get_xsize() || image->get_ysize()!=with->get_ysize() || image->get_zsize()!=with->get_zsize())  throw InvalidCallException("Error: TomoWedgeFscCmp requires 2 images the same size");

	int nx=image->get_xsize();
	int ny=image->get_ysize();
	int nz=image->get_zsize();
	
	// User can pass in a sigma vector if they like, otherwise we call it 1/10 of the standard deviation in each shell
	// This has proven to be a good cutoff for identifying the missing wedge voxels, without throwing away too many
	// voxels if the image is complete in Fourier space (no wedge)
	vector<float> sigmaimg; 
	if (params.has_key("sigmaimg")) sigmaimg=params["sigmaimg"];
	else {
		sigmaimg=image->calc_radial_dist(nx/2,0,1,4);
		for (int i=0; i<nx/2; i++) sigmaimg[i]*=.1;
	}
	
	vector<float> sigmawith;
	if (params.has_key("sigmawith")) sigmawith = params["sigmawith"];
	else {
		sigmawith=image->calc_radial_dist(nx/2,0,1,4);
		for (int i=0; i<nx/2; i++) sigmawith[i]*=.1;
	}

	int negative = params.set_default("negative",1);
	
	double sum=0;
	double sumsq1=0;
	double sumsq2=0;
	double norm=0;
	for (int z=0; z<nz; z++) {
		for (int y=0; y<ny; y++) {
			for (int x=0; x<nx; x+=2) {
				int r=int(Util::hypot3(x/2,y<ny/2?y:ny-y,z<nz/2?z:nz-z));	// origin at 0,0; periodic
				
				float v1r=image->get_value_at(x,y,z);
				float v1i=image->get_value_at(x+1,y,z);
				float v1=Util::square_sum(v1r,v1i);
				if (v1<sigmaimg[r]) continue;
				
				float v2r=with->get_value_at(x,y,z);
				float v2i=with->get_value_at(x+1,y,z);
				float v2=Util::square_sum(v2r,v2i);
				if (v2<sigmawith[r]) continue;
				
				sum+=v1r*v2r+v1i*v2i;
				sumsq1+=v1;
				if (Util::is_nan(sumsq1)) { printf("TomoWedgeCccCmp: NaN encountered: %d %d %d %f %f %f\n",x,y,z,v1r,v1i,v1); }
				sumsq2+=v2;
				norm+=1.0;
			}
		}
	}
	image->set_attr("fft_overlap",(float)(2.0*norm/(image->get_xsize()*image->get_ysize()*image->get_zsize())));
//	printf("%f\t%f\t%f\t%f\t%f\n",s1,s2,sumsq1,sumsq2,norm);
	
	if (negative) sum*=-1.0;
	return float(sum/(sqrt(sumsq1)*sqrt(sumsq2)));
}

float TomoFscCmp::cmp(EMData * image, EMData *with) const
{
	ENTERFUNC;
	bool usecpu = 1;
	bool del_imagefft = 0;
	bool del_withfft = 0;
	float score = 1.0f;
	
	//get parameters
	if (!image->has_attr("spt_wedge_mean") || !image->has_attr("spt_wedge_sigma"))  throw InvalidCallException("Rubbish!!! Image Subtomogram does not have mena and/or sigma amps metadata");
	// BAD DESIGN!!!! The fact that I have to load attrs into a variable before I can do operations on them is silly
	
	//get threshold information
	float image_meanwedgeamp = image->get_attr("spt_wedge_mean");
	float image_sigmawedgeamp = image->get_attr("spt_wedge_sigma");
	
	// Use with amp stats if avaliable otherwise effectivly ignore
	float with_meanwedgeamp = image_meanwedgeamp;
	float with_sigmawedgeamp = image_sigmawedgeamp;
	if (with->has_attr("spt_wedge_mean") && with->has_attr("spt_wedge_sigma"))
	{
		with_meanwedgeamp = with->get_attr("spt_wedge_mean");
		with_sigmawedgeamp = with->get_attr("spt_wedge_sigma");
	}
	
	// Find threshold
	float sigmas = params.set_default("sigmas",5.0f);
	float img_amp_thres = pow(image_meanwedgeamp + sigmas*image_sigmawedgeamp, 2.0f);
	float with_amp_thres = pow(with_meanwedgeamp + sigmas*with_sigmawedgeamp, 2.0f);
	
	// take negative of score
	float negative = (float)params.set_default("negative", 1.0f);
	if (negative) negative=-1.0; else negative=1.0;
	//get apix, use param apix, if not specified use apix_x, if this is not specified then apix=1.0
	float apix = params.set_default("apix",image->get_attr_default("apix_x", 1.0f));
	//get min and max res
	float minres = params.set_default("minres",std::numeric_limits<float>::max());
	float maxres = params.set_default("maxres", 0.0f);
	
	//Check to ensure that images are complex
	EMData* image_fft = image;
	EMData* with_fft = with;
	
#ifdef EMAN2_USING_CUDA	
	//do CUDA FFT, does not allow minres, maxres yet
	if(EMData::usecuda == 1 && image->getcudarwdata() && with->getcudarwdata()) {
		if(!image->is_complex()){
			del_imagefft = 1;
			image_fft = image->do_fft_cuda();
		}
		if(!with->is_complex()){
			del_withfft = 1;
			with_fft = with->do_fft_cuda();
		}
		score = fsc_tomo_cmp_cuda(image_fft->getcudarwdata(), with_fft->getcudarwdata(), img_amp_thres, with_amp_thres, 0.0, 0.0, image_fft->get_xsize(), image_fft->get_ysize(), image_fft->get_zsize());
		usecpu = 0;
	}
#endif	
	if(usecpu){
		if(!image->is_complex()){
			del_imagefft = 1;
			image_fft = image->do_fft();
		}
		if(!with->is_complex()){
			del_withfft = 1;
			with_fft = with->do_fft();
		}
		
		//loop over all voxels
		int count = 0;
		double sum_imgamp_sq = 0.0;
		double sum_withamp_sq = 0.0;
		double cong = 0.0;
		float* img_data = image_fft->get_data();
		float* with_data = with_fft->get_data();
	
		int nx  = image_fft->get_xsize();
		int ny  = image_fft->get_ysize();
		int nz  = image_fft->get_zsize();
		int ny2 = ny/2;
		int nz2 = nz/2;
	
		//compute FSC
		int ii, kz, ky;
		for (int iz = 0; iz < nz; iz++) {
			if(iz > nz2) kz = nz-iz; else kz=iz;
			for (int iy = 0; iy < ny; iy++) {
				if(iy > ny2) ky = ny-iy; else ky=iy;
				for (int ix = 0; ix < nx; ix+=2) {
					//compute spatial frequency and convert to resolution stat
					float freq = std::sqrt(kz*kz + ky*ky + ix*ix/4.0f)/float(nz);
					float reso = apix/freq;
					
					//only look within a resolution domain
					if(reso < minres && reso > maxres){
						ii = ix + (iy  + iz * ny)* nx;
						float img_r = img_data[ii];
						float img_i = img_data[ii+1];
						float with_r = with_data[ii];
						float with_i = with_data[ii+1];
						double img_amp_sq = img_r*img_r + img_i*img_i;
						double with_amp_sq = with_r*with_r + with_i*with_i;

						if((img_amp_sq >  img_amp_thres) && (with_amp_sq >  with_amp_thres)){
							count ++;
							sum_imgamp_sq += img_amp_sq;
							sum_withamp_sq += with_amp_sq;
							cong += img_r*with_r + img_i*with_i;
						}
					}
				}
			}
		}
		
		if(count > 0){ 
			score = (float)(cong/sqrt(sum_imgamp_sq*sum_withamp_sq));
		}else{
			score = 1.0f;
		}
	}
	
	//avoid mem leaks
	if(del_imagefft) delete image_fft;
	if(del_withfft) delete with_fft;
	
	return negative*score;
	EXITFUNC;
}

// Even though this uses doubles, it might be wise to recode it row-wise
// to avoid numerical errors on large images
float QuadMinDotCmp::cmp(EMData * image, EMData *with) const
{
	ENTERFUNC;
	validate_input_args(image, with);

	if (image->get_zsize()!=1) throw InvalidValueException(0, "QuadMinDotCmp supports 2D only");

	int nx=image->get_xsize();
	int ny=image->get_ysize();

	int normalize = params.set_default("normalize", 0);
	float negative = (float)params.set_default("negative", 1);

	if (negative) negative=-1.0; else negative=1.0;

	double result[4] = { 0,0,0,0 }, sq1[4] = { 0,0,0,0 }, sq2[4] = { 0,0,0,0 } ;

	vector<int> image_saved_offsets = image->get_array_offsets();
	vector<int> with_saved_offsets = with->get_array_offsets();
	image->set_array_offsets(-nx/2,-ny/2);
	with->set_array_offsets(-nx/2,-ny/2);
	int i,x,y;
	for (y=-ny/2; y<ny/2; y++) {
		for (x=-nx/2; x<nx/2; x++) {
			int quad=(x<0?0:1) + (y<0?0:2);
			result[quad]+=(*image)(x,y)*(*with)(x,y);
			if (normalize) {
				sq1[quad]+=(*image)(x,y)*(*image)(x,y);
				sq2[quad]+=(*with)(x,y)*(*with)(x,y);
			}
		}
	}
	image->set_array_offsets(image_saved_offsets);
	with->set_array_offsets(with_saved_offsets);

	if (normalize) {
		for (i=0; i<4; i++) result[i]/=sqrt(sq1[i]*sq2[i]);
	} else {
		for (i=0; i<4; i++) result[i]/=nx*ny/4;
	}

	float worst=static_cast<float>(result[0]);
	for (i=1; i<4; i++) if (static_cast<float>(result[i])<worst) worst=static_cast<float>(result[i]);

	EXITFUNC;
	return (float) (negative*worst);
}

float OptVarianceCmp::cmp(EMData * image, EMData *with) const
{
	ENTERFUNC;
	validate_input_args(image, with);

	int keepzero = params.set_default("keepzero", 1);
	int invert = params.set_default("invert",0);
	int matchfilt = params.set_default("matchfilt",1);
	int matchamp = params.set_default("matchamp",0);
	int radweight = params.set_default("radweight",0);
	int dbug = params.set_default("debug",0);

	size_t size = (size_t)image->get_xsize() * image->get_ysize() * image->get_zsize();


	EMData *with2=NULL;
	if (matchfilt) {
		EMData *a = image->do_fft();
		EMData *b = with->do_fft();

		vector <float> rfa=a->calc_radial_dist(a->get_ysize()/2,0.0f,1.0f,1);
		vector <float> rfb=b->calc_radial_dist(b->get_ysize()/2,0.0f,1.0f,1);

		float avg=0;
		for (size_t i=0; i<a->get_ysize()/2.0f; i++) {
			rfa[i]=(rfb[i]==0?0.0f:(rfa[i]/rfb[i]));
			avg+=rfa[i];
		}

		avg/=a->get_ysize()/2.0f;
		for (size_t i=0; i<a->get_ysize()/2.0f; i++) {
			if (rfa[i]>avg*10.0) rfa[i]=10.0;			// If some particular location has a small but non-zero value, we don't want to overcorrect it
		}
		rfa[0]=0.0;

		if (dbug) b->write_image("a.hdf",-1);

		b->apply_radial_func(0.0f,1.0f/a->get_ysize(),rfa);
		with2=b->do_ift();

		if (dbug) b->write_image("a.hdf",-1);
		if (dbug) a->write_image("a.hdf",-1);

/*		if (dbug) {
			FILE *out=fopen("a.txt","w");
			for (int i=0; i<a->get_ysize()/2.0; i++) fprintf(out,"%d\t%f\n",i,rfa[i]);
			fclose(out);

			out=fopen("b.txt","w");
			for (int i=0; i<a->get_ysize()/2.0; i++) fprintf(out,"%d\t%f\n",i,rfb[i]);
			fclose(out);
		}*/


		delete a;
		delete b;

		if (dbug) {
			with2->write_image("a.hdf",-1);
			image->write_image("a.hdf",-1);
		}

//		with2->process_inplace("matchfilt",Dict("to",this));
//		x_data = with2->get_data();
	}

	// This applies the individual Fourier amplitudes from 'image' and
	// applies them to 'with'
	if (matchamp) {
		EMData *a = image->do_fft();
		EMData *b = with->do_fft();
		size_t size2 = (size_t)a->get_xsize() * a->get_ysize() * a->get_zsize();

		a->ri2ap();
		b->ri2ap();

		const float *const ad=a->get_const_data();
		float * bd=b->get_data();

		for (size_t i=0; i<size2; i+=2) bd[i]=ad[i];
		b->update();

		b->ap2ri();
		with2=b->do_ift();
//with2->write_image("a.hdf",-1);
		delete a;
		delete b;
	}

	const float * x_data;
	if (with2) x_data=with2->get_const_data();
	else x_data = with->get_const_data();
	const float *const y_data = image->get_const_data();

	size_t nx = image->get_xsize();
	float m = 0;
	float b = 0;

	// This will write the x vs y file used to calculate the density
	// optimization. This behavior may change in the future
	if (dbug) {
		FILE *out=fopen("dbug.optvar.txt","w");
		if (out) {
			for (size_t i=0; i<size; i++) {
				if ( !keepzero || (x_data[i] && y_data[i])) fprintf(out,"%g\t%g\n",x_data[i],y_data[i]);
			}
			fclose(out);
		}
	}


	Util::calc_least_square_fit(size, x_data, y_data, &m, &b, keepzero);
	if (m == 0) {
		m = FLT_MIN;
	}
	b = -b / m;
	m = 1.0f / m;

	// While negative slopes are really not a valid comparison in most cases, we
	// still want to detect these instances, so this if is removed
/*	if (m < 0) {
		b = 0;
		m = 1000.0;
	}*/

	double  result = 0;
	int count = 0;

	if (radweight) {
		if (image->get_zsize()!=1) throw ImageDimensionException("radweight option is 2D only");
		if (keepzero) {
			for (size_t i = 0,y=0; i < size; y++) {
				for (size_t x=0; x<nx; i++,x++) {
					if (y_data[i] && x_data[i]) {
#ifdef	_WIN32
						if (invert) result += Util::square(x_data[i] - (y_data[i]-b)/m)*(_hypot((float)x,(float)y)+nx/4.0);
						else result += Util::square((x_data[i] * m) + b - y_data[i])*(_hypot((float)x,(float)y)+nx/4.0);
#else
						if (invert) result += Util::square(x_data[i] - (y_data[i]-b)/m)*(hypot((float)x,(float)y)+nx/4.0);
						else result += Util::square((x_data[i] * m) + b - y_data[i])*(hypot((float)x,(float)y)+nx/4.0);
#endif
						count++;
					}
				}
			}
			result/=count;
		}
		else {
			for (size_t i = 0,y=0; i < size; y++) {
				for (size_t x=0; x<nx; i++,x++) {
#ifdef	_WIN32
					if (invert) result += Util::square(x_data[i] - (y_data[i]-b)/m)*(_hypot((float)x,(float)y)+nx/4.0);
					else result += Util::square((x_data[i] * m) + b - y_data[i])*(_hypot((float)x,(float)y)+nx/4.0);
#else
					if (invert) result += Util::square(x_data[i] - (y_data[i]-b)/m)*(hypot((float)x,(float)y)+nx/4.0);
					else result += Util::square((x_data[i] * m) + b - y_data[i])*(hypot((float)x,(float)y)+nx/4.0);
#endif
				}
			}
			result = result / size;
		}
	}
	else {
		if (keepzero) {
			for (size_t i = 0; i < size; i++) {
				if (y_data[i] && x_data[i]) {
					if (invert) result += Util::square(x_data[i] - (y_data[i]-b)/m);
					else result += Util::square((x_data[i] * m) + b - y_data[i]);
					count++;
				}
			}
			result/=count;
		}
		else {
			for (size_t i = 0; i < size; i++) {
				if (invert) result += Util::square(x_data[i] - (y_data[i]-b)/m);
				else result += Util::square((x_data[i] * m) + b - y_data[i]);
			}
			result = result / size;
		}
	}
	scale = m;
	shift = b;

	image->set_attr("ovcmp_m",m);
	image->set_attr("ovcmp_b",b);
	if (with2) delete with2;
	EXITFUNC;

#if 0
	return (1 - result);
#endif

	return static_cast<float>(result);
}

float PhaseCmp::cmp(EMData * image, EMData *with) const
{
	ENTERFUNC;

	int snrweight = params.set_default("snrweight", 0);
	int snrfn = params.set_default("snrfn",0);
	int ampweight = params.set_default("ampweight",0);
	int zeromask = params.set_default("zeromask",0);
	float minres = params.set_default("minres",500.0f);
	float maxres = params.set_default("maxres",10.0f);

	if (snrweight && snrfn) throw InvalidCallException("SNR weight and SNRfn cannot both be set in the phase comparator");

	EMData *image_fft = NULL;
	EMData *with_fft = NULL;

	int ny = image->get_ysize();
//	int np = (int) ceil(Ctf::CTFOS * sqrt(2.0f) * ny / 2) + 2;
	int np = 0;
	vector<float> snr;

	// weighting based on SNR estimate from CTF
	if (snrweight) {
		Ctf *ctf = NULL;
		if (!image->has_attr("ctf")) {
			if (!with->has_attr("ctf")) throw InvalidCallException("SNR weight with no CTF parameters");
			ctf=with->get_attr("ctf");
		}
		else ctf=image->get_attr("ctf");

		float ds=1.0f/(ctf->apix*ny);
		snr=ctf->compute_1d(ny,ds,Ctf::CTF_SNR); // note that this returns ny/2 values
		for (int i=0; i<snr.size(); i++) {
			if (snr[i]<=0) snr[i]=0.001;		// make sure that points don't get completely excluded due to SNR estimation issues, or worse, contribute with a negative weight
		}
		if(ctf) {delete ctf; ctf=0;}
		np=snr.size();
	}
	// weighting based on empirical SNR function (not really good)
	else if (snrfn==1) {
		np=ny/2;
		snr.resize(np);

		for (int i = 0; i < np; i++) {
			float x2 = 10.0f*i/np;
			snr[i] = x2 * exp(-x2);
		}
	}
	// no weighting
	else {
		np=ny/2;
		snr.resize(np);

		for (int i = 0; i < np; i++) snr[i]=1.0;
	}

	// Min/max modifications to weighting
	float pmin,pmax;
	if (minres>0) pmin=((float)image->get_attr("apix_x")*image->get_ysize())/minres;		//cutoff in pixels, assume square
	else pmin=0;
	if (maxres>0) pmax=((float)image->get_attr("apix_x")*image->get_ysize())/maxres;
	else pmax=0;

//	printf("%f\t%f\n",pmin,pmax);

	// We use 'soft' edges for the Fourier cutoffs to minimize issues with pixel interpolation
	for (int i = 0; i < np; i++) {
		if (pmin>0) snr[i]*=(tanh(5.0f*(i-pmin)/pmin)+1.0f)/2.0f;
		if (pmax>0) snr[i]*=(1.0f-tanh(i-pmax))/2.0f;
//		printf("%d\t%f\n",i,snr[i]);
	}

	if (zeromask) {
		image_fft=image->copy();
		with_fft=with->copy();
		
		if (image_fft->is_complex()) image_fft->do_ift_inplace();
		if (with_fft->is_complex()) with_fft->do_ift_inplace();
		
		int sz=image_fft->get_xsize()*image_fft->get_ysize()*image_fft->get_zsize();
		float *d1=image_fft->get_data();
		float *d2=with_fft->get_data();
		
		for (int i=0; i<sz; i++) {
			if (d1[i]==0.0 || d2[i]==0.0) { d1[i]=0.0; d2[i]=0.0; }
		}
		
		image_fft->update();
		with_fft->update();
		image_fft->do_fft_inplace();
		with_fft->do_fft_inplace();
		image_fft->set_attr("free_me",1); 
		with_fft->set_attr("free_me",1); 
	}
	else {
		if (image->is_complex()) image_fft=image;
		else {
			image_fft=image->do_fft();
			image_fft->set_attr("free_me",1);
		}
		
		if (with->is_complex()) with_fft=with;
		else {
			with_fft=with->do_fft();
			with_fft->set_attr("free_me",1);
		}
	}
// 	image_fft->ri2ap();
// 	with_fft->ri2ap();

	const float *const image_fft_data = image_fft->get_const_data();
	const float *const with_fft_data = with_fft->get_const_data();
	double sum = 0;
	double norm = FLT_MIN;
	size_t i = 0;
//	int nx=image_fft->get_xsize();
	    ny=image_fft->get_ysize();
	int nz=image_fft->get_zsize();
	int nx2=image_fft->get_xsize()/2;
	int ny2=image_fft->get_ysize()/2;
//	int nz2=image_fft->get_zsize()/2;

	// This can never happen any more...
	if (np==0) {
		for (int z = 0; z < nz; z++){
			for (int y = 0; y < ny; y++) {
				for (int x = 0; x < nx2; x ++) {
					float a;
					if (ampweight) a= (float)hypot(with_fft_data[i],with_fft_data[i+1]);
					else a=1.0f;
					sum += Util::angle_err_ri(image_fft_data[i],image_fft_data[i+1],with_fft_data[i],with_fft_data[i+1]) * a;
					norm += a;
					i += 2;
				}
			}
		}
		
	}
	else if (nz==1) {
		for (int y = 0; y < ny; y++) {
			for (int x = 0; x < nx2; x ++) {
				int r=Util::hypot_fast_int(x,y>ny/2?ny-y:y);
				if (r>=ny2) { i+=2; continue; }		// we only have snr values to the box radius
				
				float a;
				if (ampweight) a= (float)hypot(with_fft_data[i],with_fft_data[i+1]);
				else a=1.0f;
				a*=snr[r];
				sum += Util::angle_err_ri(image_fft_data[i],image_fft_data[i+1],with_fft_data[i],with_fft_data[i+1]) * a;
				norm += a;
				i += 2;
			}
		}
	}
	else {
		for (int z = 0; z < nz; z++){
			for (int y = 0; y < ny; y++) {
				for (int x = 0; x < nx2; x ++) {
					int r=(int)Util::hypot3(x,y>ny/2?ny-y:y,z>nz/2?nz-z:z);
					if (r>=ny2) { i+=2; continue; }		// we only have snr values to the box radius
					
					float a;
					if (ampweight) a= (float)hypot(with_fft_data[i],with_fft_data[i+1]);
					else a=1.0f;
					a*=snr[r];
					sum += Util::angle_err_ri(image_fft_data[i],image_fft_data[i+1],with_fft_data[i],with_fft_data[i+1]) * a;
					norm += a;
					i += 2;
				} 
			}
		}
		
	}

	EXITFUNC;

	if( image_fft->has_attr("free_me") )
	{
		delete image_fft;
		image_fft = 0;
	}
	if( with_fft->has_attr("free_me") )
	{
		delete with_fft;
		with_fft = 0;
	}
#if 0
	return (1.0f - sum / norm);
#endif
	return (float)(sum / norm);
}

float FRCCmp::cmp(EMData * image, EMData * with) const
{
	ENTERFUNC;
	validate_input_args(image, with);

	int snrweight = params.set_default("snrweight", 0);
	int ampweight = params.set_default("ampweight", 0);
	int sweight = params.set_default("sweight", 1);
	int nweight = params.set_default("nweight", 0);
	int zeromask = params.set_default("zeromask",0);
	float minres = params.set_default("minres",500.0f);
	float maxres = params.set_default("maxres",10.0f);

	vector < float >fsc;
	bool use_cpu = true;

	if (use_cpu) {
		if (zeromask) {
			image=image->copy();
			with=with->copy();
		
			int sz=image->get_xsize()*image->get_ysize()*image->get_zsize();
			float *d1=image->get_data();
			float *d2=with->get_data();
		
			for (int i=0; i<sz; i++) {
				if (d1[i]==0.0 || d2[i]==0.0) { d1[i]=0.0; d2[i]=0.0; }
			}
		
			image->update();
			with->update();
			image->do_fft_inplace();
			with->do_fft_inplace();
			image->set_attr("free_me",1); 
			with->set_attr("free_me",1); 
		}


		if (!image->is_complex()) {
			image=image->do_fft(); 
			image->set_attr("free_me",1); 
		}
		if (!with->is_complex()) { 
			with=with->do_fft(); 
			with->set_attr("free_me",1); 
		}

		fsc = image->calc_fourier_shell_correlation(with,1);
	}
	
	int ny = image->get_ysize();
	int ny2=ny/2+1;

	// The fast hypot here was supposed to speed things up. Little effect
// 	if (image->get_zsize()>1) fsc = image->calc_fourier_shell_correlation(with,1);
// 	else {
// 		double *sxy = (double *)malloc(ny2*sizeof(double)*4);
// 		double *sxx = sxy+ny2;
// 		double *syy = sxy+2*ny2;
// 		double *norm= sxy+3*ny2;
//
// 		float *df1=image->get_data();
// 		float *df2=with->get_data();
// 		int nx2=image->get_xsize();
//
// 		for (int y=-ny/2; y<ny/2; y++) {
// 			for (int x=0; x<nx2/2; x++) {
// 				if (x==0 && y<0) continue;	// skip Friedel pair
// 				short r=Util::hypot_fast_int(x,y);
// 				if (r>ny2-1) continue;
// 				int l=x*2+(y<0?ny+y:y)*nx2;
// 				sxy[r]+=df1[l]*df2[l]+df1[l+1]*df2[l+1];
// 				sxx[r]+=df1[l]*df1[l];
// 				syy[r]+=df2[l]*df2[l];
// 				norm[r]+=1.0;
// 			}
// 		}
// 		fsc.resize(ny2*3);
// 		for (int r=0; r<ny2; r++) {
// 			fsc[r]=r*0.5/ny2;
// 			fsc[ny2+r]=sxy[r]/(sqrt(sxx[r])*sqrt(syy[r]));
// 			fsc[ny2*2+r]=norm[r];
// 		}
// 		free(sxy);
// 	}

	vector<float> snr;
	if (snrweight) {
		Ctf *ctf = NULL;
		if (!image->has_attr("ctf")) {
			if (!with->has_attr("ctf")) throw InvalidCallException("SNR weight with no CTF parameters");
			ctf=with->get_attr("ctf");
		}
		else ctf=image->get_attr("ctf");

		float ds=1.0f/(ctf->apix*ny);
		snr=ctf->compute_1d(ny,ds,Ctf::CTF_SNR);
		for (int i=0; i<snr.size(); i++) {
			if (snr[i]<=0) snr[i]=0.001;		// make sure that points don't get completely excluded due to SNR estimation issues, or worse, contribute with a negative weight
		}
		if(ctf) {delete ctf; ctf=0;}
	}

	vector<float> amp;
	if (ampweight) amp=image->calc_radial_dist(ny/2,0,1,0);

	// Min/max modifications to weighting
	float pmin,pmax;
	if (minres>0) pmin=((float)image->get_attr("apix_x")*image->get_ysize())/minres;		//cutoff in pixels, assume square
	else pmin=0;
	if (maxres>0) pmax=((float)image->get_attr("apix_x")*image->get_ysize())/maxres;
	else pmax=0;

	double sum=0.0, norm=0.0;

	for (int i=0; i<ny/2; i++) {
		double weight=1.0;
		if (sweight) weight*=fsc[(ny2)*2+i];
		if (ampweight) weight*=amp[i];
		if (snrweight) weight*=snr[i];
//		if (snrweight)  {
//			if (snr[i]>0) weight*=sqrt(snr[i]);
//			else weight=0;
//		}
//if(snr[i]<0) printf("snr[%d] = %1.5g\n",i,snr[i]);
		if (pmin>0) weight*=(tanh(5.0*(i-pmin)/pmin)+1.0)/2.0;
		if (pmax>0) weight*=(1.0-tanh(i-pmax))/2.0;
		
		sum+=weight*fsc[ny2+i];
		norm+=weight;
//		printf("%d\t%f\t%f\n",i,weight,fsc[ny/2+1+i]);
	}

	// This performs a weighting that tries to normalize FRC by correcting from the number of particles represented by the average
	sum/=norm;
	if (nweight && with->get_attr_default("ptcl_repr",0) && sum>=0 && sum<1.0) {
		sum=sum/(1.0-sum);							// convert to SNR
		sum/=(float)with->get_attr_default("ptcl_repr",0);	// divide by ptcl represented
		sum=sum/(1.0+sum);							// convert back to correlation
	}

	if (image->has_attr("free_me")) delete image;
	if (with->has_attr("free_me")) delete with;

	EXITFUNC;

	if (!Util::goodf(&sum)) sum=-2.0;	// normally should be >-1.0

	//.Note the negative! This is because EMAN2 follows the convention that
	// smaller return values from comparitors indicate higher similarity -
	// this enables comparitors to be used in a generic fashion.
	return (float)-sum;
}

float OptSubCmp::cmp(EMData * image, EMData * with) const
{
	ENTERFUNC;
	validate_input_args(image, with);

// 	int snrweight = params.set_default("snrweight", 0);
// 	int ampweight = params.set_default("ampweight", 0);
// 	int sweight = params.set_default("sweight", 1);
// 	int nweight = params.set_default("nweight", 0);
	int ctfweight = params.set_default("ctfweight",0);
	int zeromask = params.set_default("zeromask",0);
//	int swap = params.set_default("swap",0);
	float minres = params.set_default("minres",200.0f);
	float maxres = params.set_default("maxres",10.0f);
	EMData *mask = params.set_default("mask",(EMData *)NULL);
	
//	float ds=1.0f/((float)image->get_attr("apix_x")*(int)image->get_ysize());
	float apix=(float)image->get_attr("apix_x");

	// Sometimes we will get the "raw" image with CTF as image and sometimes as with, we always want to subtract the less noisy
	// reference, so if one has CTF parameters (even if ctfweight isn't used) we always pass it in as the primary
	EMData *diff;
	if (image->has_attr("ctf")) diff=image->process("math.sub.optimal",Dict("ref",with,"return_presigma",1,"low_cutoff_frequency",apix/minres ,"high_cutoff_frequency",apix/maxres,"ctfweight",ctfweight));
	else diff=with->process("math.sub.optimal",Dict("ref",image,"return_presigma",1,"low_cutoff_frequency",apix/minres ,"high_cutoff_frequency",apix/maxres,"ctfweight",ctfweight));

	if (mask!=NULL) diff->mult(*mask);
	if (zeromask) {
		EMData *tmp=with->process("threshold.notzero");
		diff->mult(*tmp);
		delete tmp;
	}
	
// 	diff->process_inplace("filter.highpass.tophat",Dict("cutoff_freq",(float)1.0/minres));
// 	diff->process_inplace("filter.lowpass.tophat",Dict("cutoff_freq",(float)1.0/maxres));
	
	float ret=(float)diff->get_attr("sigma")/(float)diff->get_attr("sigma_presub");
	delete diff;
	return ret;
	
	
// 	EMData *diff=image->process("math.sub.optimal",Dict("ref",with,"return_fft",1));
// 	
// 	// Very expensive to use a mask, since it requires another ift/fft pair
// 	if (mask!=NULL) {
// 		EMData *tmp = diff->do_ift();
// 		tmp->mult(*mask);
// 		delete diff;
// 		diff=tmp->do_fft();
// 		delete tmp;
// 	}
// 	
// 	// This gives us basically the 1-D power spectrum of what's left after subtraction (and masking)
// 	vector<float> dist=diff->calc_radial_dist(diff->get_ysize()/2,0.0f,1.0f,1);
// 	int s0=int(floor(1.0/(minres*ds)));
// 	int s1=int(ceil(1.0/(maxres*ds)));
// 	if (s0<2) s0=2;
// 	
// 	if (s1<=s0) throw InvalidCallException("OptSubCmp error. minres must be greater than maxres.");
// //	printf("%d %d\n",s0,s1);
// 	
// 	double sum=0.0f,sum2=0.0f;
// 	for (int i=s0; i<s1; i++) sum+=dist[i]*i;
// 	sum/=image->get_size()*s1;
// 	
// 	delete diff;
// 	return sum;
}


void EMAN::dump_cmps()
{
	dump_factory < Cmp > ();
}

map<string, vector<string> > EMAN::dump_cmps_list()
{
	return dump_factory_list < Cmp > ();
}

/* vim: set ts=4 noet: */
