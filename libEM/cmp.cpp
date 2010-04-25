/**
 * $Id$
 */

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

using namespace EMAN;

const string CccCmp::NAME = "ccc";
const string SqEuclideanCmp::NAME = "sqeuclidean";
const string DotCmp::NAME = "dot";
const string TomoDotCmp::NAME = "dot.tomo";
const string QuadMinDotCmp::NAME = "quadmindot";
const string OptVarianceCmp::NAME = "optvariance";
const string PhaseCmp::NAME = "phase";
const string FRCCmp::NAME = "frc";

template <> Factory < Cmp >::Factory()
{
	force_add<CccCmp>();
	force_add<SqEuclideanCmp>();
	force_add<DotCmp>();
	force_add<TomoDotCmp>();
	force_add<QuadMinDotCmp>();
	force_add<OptVarianceCmp>();
	force_add<PhaseCmp>();
	force_add<FRCCmp>();
//	force_add<XYZCmp>();
}

void Cmp::validate_input_args(const EMData * image, const EMData *with) const
{
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
	var2 = var1/double(n) - avg2*avg2;
	ccc = ccc/double(n) - avg1*avg2;
	ccc /= sqrt(var1*var2);
	ccc *= negative;
	return static_cast<float>(ccc);
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
		with = withorig->process("normalize.toimage",Dict("to",image));
		with->set_attr("deleteme",1);
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
		//  it looks like it could work in 3D, but it is not, really.
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

		}else{ //This 3D code is incorrect, but it is the best I can do now 01/09/06 PAP
		int ky, kz;
		int ny2 = ny/2; int nz2 = nz/2;
		for ( int iz = 0; iz <= nz-1; iz++) {
			if(iz>nz2) kz=iz-nz; else kz=iz;
			for ( int iy = 0; iy <= ny-1; iy++) {
				if(iy>ny2) ky=iy-ny; else ky=iy;
				for ( int ix = 0; ix <= lsd2-1; ix++) {
				// Skip Friedel related values
				if(ix>0 || (kz>=0 && (ky>=0 || kz!=0))) {
						size_t ii = ix + (iy  + iz * ny)* lsd2;
						result += (x_data[ii] - y_data[ii])*double(x_data[ii] - y_data[ii]);
					}
				}
			}
		}
		n = ((float)nx*(float)ny*(float)nz*(float)nx*(float)ny*(float)nz)/2.0f;
		}
	} else {		// real space
		size_t totsize = image->get_xsize()*image->get_ysize()*image->get_zsize();
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
	return static_cast<float>(result);
}


// Even though this uses doubles, it might be wise to recode it row-wise
// to avoid numerical errors on large images
float DotCmp::cmp(EMData* image, EMData* with) const
{
	ENTERFUNC;
	validate_input_args(image, with);

	const float *const x_data = image->get_const_data();
	const float *const y_data = with->get_const_data();

	int normalize = params.set_default("normalize", 0);
	float negative = (float)params.set_default("negative", 1);

	if (negative) negative=-1.0; else negative=1.0;
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
		} else  result /= ((float)nx*(float)ny*(float)nz*(float)nx*(float)ny*(float)nz);

		} else { //This 3D code is incorrect, but it is the best I can do now 01/09/06 PAP
		int ky, kz;
		int ny2 = ny/2; int nz2 = nz/2;
		for ( int iz = 0; iz <= nz-1; ++iz) {
			if(iz>nz2) kz=iz-nz; else kz=iz;
			for ( int iy = 0; iy <= ny-1; ++iy) {
				if(iy>ny2) ky=iy-ny; else ky=iy;
				for ( int ix = 0; ix <= lsd2-1; ++ix) {
					// Skip Friedel related values
					if(ix>0 || (kz>=0 && (ky>=0 || kz!=0))) {
						size_t ii = ix + (iy  + iz * ny)* lsd2;
						result += x_data[ii] * double(y_data[ii]);
					}
				}
			}
		}
		if( normalize ) {
		//  still incorrect
		double square_sum1 = 0., square_sum2 = 0.;
		int ky, kz;
		int ny2 = ny/2; int nz2 = nz/2;
		for ( int iz = 0; iz <= nz-1; ++iz) {
			if(iz>nz2) kz=iz-nz; else kz=iz;
			for ( int iy = 0; iy <= ny-1; ++iy) {
				if(iy>ny2) ky=iy-ny; else ky=iy;
				for ( int ix = 0; ix <= lsd2-1; ++ix) {
					// Skip Friedel related values
					if(ix>0 || (kz>=0 && (ky>=0 || kz!=0))) {
						size_t ii = ix + (iy  + iz * ny)* lsd2;
						square_sum1 += x_data[ii] * double(x_data[ii]);
						square_sum2 += y_data[ii] * double(y_data[ii]);
					}
				}
			}
		}
		result /= sqrt(square_sum1*square_sum2);
		} else result /= ((float)nx*(float)ny*(float)nz*(float)nx*(float)ny*(float)nz/2);
		}
	} else {
		size_t totsize = image->get_xsize() * image->get_ysize() * image->get_zsize();

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


float TomoDotCmp::cmp(EMData * image, EMData *with) const
{
	ENTERFUNC;
	float threshold = params.set_default("threshold",0.f);
	if (threshold < 0.0f) throw InvalidParameterException("The threshold parameter must be greater than or equal to zero");

	if ( threshold > 0) {
		EMData* ccf = params.set_default("ccf",(EMData*) NULL);
		bool ccf_ownership = false;
		if (!ccf) {
			ccf = image->calc_ccf(with);
			ccf_ownership = true;
		}
		bool norm = params.set_default("norm",false);
		if (norm) ccf->process_inplace("normalize");
		int tx = params.set_default("tx",0); int ty = params.set_default("ty",0); int tz = params.set_default("tz",0);
		float best_score = ccf->get_value_at_wrap(tx,ty,tz)/static_cast<float>(image->get_size());
		EMData* ccf_fft = ccf->do_fft();// so cuda works, or else we could do an fft_inplace - though honestly doing an fft inplace is less efficient anyhow
		if (ccf_ownership) delete ccf; ccf = 0;
		ccf_fft->process_inplace("threshold.binary.fourier",Dict("value",threshold));
		float map_sum =  ccf_fft->get_attr("mean");
		if (map_sum == 0.0f) throw UnexpectedBehaviorException("The number of voxels in the Fourier image with an amplitude above your threshold is zero. Please adjust your parameters");
		best_score /= map_sum;
		delete ccf_fft; ccf_fft = 0;
		return -best_score;
	} else {
		return -image->dot(with);
	}


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

	size_t size = image->get_xsize() * image->get_ysize() * image->get_zsize();


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
		size_t size2 = a->get_xsize() * a->get_ysize() * a->get_zsize();

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

#ifdef EMAN2_USING_CUDA
 	if (image->gpu_operation_preferred()) {
// 		cout << "Cuda cmp" << endl;
 		EXITFUNC;
 		return cuda_cmp(image,with);
 	}
#endif

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
	int nx=image_fft->get_xsize();
	    ny=image_fft->get_ysize();
	int nz=image_fft->get_zsize();
	int nx2=image_fft->get_xsize()/2;
	int ny2=image_fft->get_ysize()/2;
	int nz2=image_fft->get_zsize()/2;

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

#ifdef EMAN2_USING_CUDA
#include "cuda/cuda_cmp.h"
float PhaseCmp::cuda_cmp(EMData * image, EMData *with) const
{
	ENTERFUNC;
	validate_input_args(image, with);

	typedef vector<EMData*> EMDatas;
	static EMDatas hist_pyramid;
	static EMDatas norm_pyramid;
	static EMData weighting;
	static int image_size = 0;

	int size;
	EMData::CudaDataLock imagelock(image);
	EMData::CudaDataLock withlock(with);

	if (image->is_complex()) {
		size = image->get_xsize();
	} else {
		int nx = image->get_xsize()+2;
		nx -= nx%2;
		size = nx*image->get_ysize()*image->get_zsize();
	}
	if (size != image_size) {
		for(unsigned int i =0; i < hist_pyramid.size(); ++i) {
			delete hist_pyramid[i];
			delete norm_pyramid[i];
		}
		hist_pyramid.clear();
		norm_pyramid.clear();
		int s = size;
		if (s < 1) throw UnexpectedBehaviorException("The image is 0 size");
		int p2 = 1;
		while ( s != 1 ) {
			s /= 2;
			p2 *= 2;
		}
		if ( p2 != size ) {
			p2 *= 2;
			s = p2;
		}
		if (s != 1) s /= 2;
		while (true) {
			EMData* h = new EMData();
			h->set_size_cuda(s); h->to_value(0.0);
			hist_pyramid.push_back(h);
			EMData* n = new EMData();
			n->set_size_cuda(s); n->to_value(0.0);
			norm_pyramid.push_back(n);
			if ( s == 1) break;
			s /= 2;
		}
		int nx = image->get_xsize()+2;
		nx -= nx%2; // for Fourier stuff
		int ny = image->get_ysize();
		int nz = image->get_zsize();
		weighting.set_size_cuda(nx,ny,nz);
		// Size of weighting need only be half this, but does that translate into faster code?
		weighting.set_size_cuda(nx/2,ny,nz);
		float np = (int) ceil(Ctf::CTFOS * sqrt(2.0f) * ny / 2) + 2;
		EMDataForCuda tmp = weighting.get_data_struct_for_cuda();
		calc_phase_weights_cuda(&tmp,np);
		//weighting.write_image("phase_wieghts.hdf");
		image_size = size;
	}

	EMDataForCuda hist[hist_pyramid.size()];
	EMDataForCuda norm[hist_pyramid.size()];

	EMDataForCuda wt = weighting.get_data_struct_for_cuda();
	EMData::CudaDataLock lock1(&weighting);
	for(unsigned int i = 0; i < hist_pyramid.size(); ++i ) {
		hist[i] = hist_pyramid[i]->get_data_struct_for_cuda();
		hist_pyramid[i]->cuda_lock();
		norm[i] = norm_pyramid[i]->get_data_struct_for_cuda();
		norm_pyramid[i]->cuda_lock();
	}

	EMData *image_fft = image->do_fft_cuda();
	EMDataForCuda left = image_fft->get_data_struct_for_cuda();
	EMData::CudaDataLock lock2(image_fft);
	EMData *with_fft = with->do_fft_cuda();
	EMDataForCuda right = with_fft->get_data_struct_for_cuda();
	EMData::CudaDataLock lock3(image_fft);

	mean_phase_error_cuda(&left,&right,&wt,hist,norm,hist_pyramid.size());
	float result;
	float* gpu_result = hist_pyramid[hist_pyramid.size()-1]->get_cuda_data();
	cudaError_t error = cudaMemcpy(&result,gpu_result,sizeof(float),cudaMemcpyDeviceToHost);
	if ( error != cudaSuccess) throw UnexpectedBehaviorException( "CudaMemcpy (host to device) in the phase comparator failed:" + string(cudaGetErrorString(error)));

	delete image_fft; image_fft=0;
	delete with_fft; with_fft=0;

	for(unsigned int i = 0; i < hist_pyramid.size(); ++i ) {
// 		hist_pyramid[i]->write_image("hist.hdf",-1); // debug
// 		norm_pyramid[i]->write_image("norm.hdf",-1); // debug
		hist_pyramid[i]->cuda_unlock();
		norm_pyramid[i]->cuda_unlock();
	}

	EXITFUNC;
	return result;

}

#endif // EMAN2_USING_CUDA


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

	static vector < float >default_snr;

// 	if (image->get_zsize() > 1) {
// 		throw ImageDimensionException("2D only");
// 	}

//	int nx = image->get_xsize();
	int ny = image->get_ysize();
	int ny2=ny/2+1;

	vector < float >fsc;

		

	fsc = image->calc_fourier_shell_correlation(with,1);

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


	//.Note the negative! This is because EMAN2 follows the convention that
	// smaller return values from comparitors indicate higher similarity -
	// this enables comparitors to be used in a generic fashion.
	return (float)-sum;
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
