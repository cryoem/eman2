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

#include "processor.h"
#include "sparx/processor_sparx.h"
#include "ctf.h"
#include "xydata.h"
#include "emdata.h"
#include "emassert.h"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_wavelet.h>
#include <gsl/gsl_wavelet2d.h>
#include <algorithm>
#include <ctime>

using namespace EMAN;
using std::reverse;

template <> Factory < Processor >::Factory()
{
	force_add(&LowpassSharpCutoffProcessor::NEW);
	force_add(&HighpassSharpCutoffProcessor::NEW);
	force_add(&LowpassGaussProcessor::NEW);
	force_add(&HighpassGaussProcessor::NEW);

	force_add(&LowpassTanhProcessor::NEW);
	force_add(&HighpassTanhProcessor::NEW);
	force_add(&HighpassButterworthProcessor::NEW);
	force_add(&AmpweightFourierProcessor::NEW);
	force_add(&Wiener2DFourierProcessor::NEW);

	force_add(&LinearPyramidProcessor::NEW);
	force_add(&LinearRampProcessor::NEW);
	force_add(&AbsoluateValueProcessor::NEW);
	force_add(&BooleanProcessor::NEW);
	force_add(&ValuePowProcessor::NEW);
	force_add(&ValueSquaredProcessor::NEW);
	force_add(&ValueSqrtProcessor::NEW);

	force_add(&ClampingProcessor::NEW);
	force_add(&NSigmaClampingProcessor::NEW);
	
	force_add(&ToZeroProcessor::NEW);
	force_add(&ToMinvalProcessor::NEW);
	force_add(&CutToZeroProcessor::NEW);
	force_add(&BinarizeProcessor::NEW);
	force_add(&CollapseProcessor::NEW);
	force_add(&LinearXformProcessor::NEW);

	force_add(&ExpProcessor::NEW);
	force_add(&RangeThresholdProcessor::NEW);
	force_add(&SigmaProcessor::NEW);
	force_add(&LogProcessor::NEW);

	force_add(&PaintProcessor::NEW);
	force_add(&MaskSharpProcessor::NEW);
	force_add(&MaskEdgeMeanProcessor::NEW);
	force_add(&MaskNoiseProcessor::NEW);
	force_add(&MaskGaussProcessor::NEW);
	force_add(&MaskGaussInvProcessor::NEW);
	force_add(&XrayPixelProcessor::NEW);
	
	force_add(&MaxShrinkProcessor::NEW);
	force_add(&MinShrinkProcessor::NEW);
	
	force_add(&MakeRadiusSquaredProcessor::NEW);
	force_add(&MakeRadiusProcessor::NEW);
	
	force_add(&ComplexNormPixel::NEW);

	force_add(&LaplacianProcessor::NEW);
	force_add(&ZeroConstantProcessor::NEW);
	
	force_add(&BoxMedianProcessor::NEW);
	force_add(&BoxSigmaProcessor::NEW);
	force_add(&BoxMaxProcessor::NEW);
	
	force_add(&MinusPeakProcessor::NEW);
	force_add(&PeakOnlyProcessor::NEW);
	force_add(&DiffBlockProcessor::NEW);

	force_add(&CutoffBlockProcessor::NEW);
	force_add(&GradientRemoverProcessor::NEW);
	force_add(&GradientPlaneRemoverProcessor::NEW);
	force_add(&VerticalStripeProcessor::NEW);
	force_add(&RealToFFTProcessor::NEW);
	force_add(&SigmaZeroEdgeProcessor::NEW);
	force_add(&RampProcessor::NEW);

	force_add(&BeamstopProcessor::NEW);
	force_add(&MeanZeroEdgeProcessor::NEW);
	force_add(&AverageXProcessor::NEW);
	force_add(&ZeroEdgeRowProcessor::NEW);
	force_add(&ZeroEdgePlaneProcessor::NEW);

	force_add(&BilateralProcessor::NEW);
	
	force_add(&ConvolutionProcessor::NEW);
	
	force_add(&NormalizeStdProcessor::NEW);
	force_add(&NormalizeUnitProcessor::NEW);
	force_add(&NormalizeUnitSumProcessor::NEW);
	force_add(&NormalizeMaskProcessor::NEW);
	force_add(&NormalizeEdgeMeanProcessor::NEW);
	force_add(&NormalizeCircleMeanProcessor::NEW);
	force_add(&NormalizeLREdgeMeanProcessor::NEW);
	force_add(&NormalizeMaxMinProcessor::NEW);
	force_add(&NormalizeRowProcessor::NEW);
	
	force_add(&HistogramBin::NEW);
	
	force_add(&NormalizeToStdProcessor::NEW);
	force_add(&NormalizeToFileProcessor::NEW);
	force_add(&NormalizeToLeastSquareProcessor::NEW);

	force_add(&RadialAverageProcessor::NEW);
	force_add(&RadialSubstractProcessor::NEW);
	force_add(&FlipProcessor::NEW);
	force_add(&MirrorProcessor::NEW);

	force_add(&AddNoiseProcessor::NEW);
	force_add(&AddSigmaNoiseProcessor::NEW);
	force_add(&AddRandomNoiseProcessor::NEW);

// 	force_add(&Phase180Processor::NEW);
	force_add(&PhaseToCenterProcessor::NEW);
	force_add(&PhaseToCornerProcessor::NEW);
// 	force_add(&FourierOriginShiftProcessor::NEW);
	force_add(&FourierToCenterProcessor::NEW);
	force_add(&FourierToCornerProcessor::NEW);
	force_add(&AutoMask2DProcessor::NEW);
	force_add(&AutoMask3DProcessor::NEW);
	force_add(&AutoMask3D2Processor::NEW);
	force_add(&AddMaskShellProcessor::NEW);

	force_add(&ToMassCenterProcessor::NEW);
	force_add(&ACFCenterProcessor::NEW);
	force_add(&SNRProcessor::NEW);

	force_add(&XGradientProcessor::NEW);
	force_add(&YGradientProcessor::NEW);
	force_add(&ZGradientProcessor::NEW);
	
	force_add(&FileFourierProcessor::NEW);

	force_add(&SymSearchProcessor::NEW);
	force_add(&LocalNormProcessor::NEW);

	force_add(&IndexMaskFileProcessor::NEW);
	force_add(&CoordinateMaskFileProcessor::NEW);
	force_add(&SetSFProcessor::NEW);

	force_add(&SmartMaskProcessor::NEW);
	force_add(&IterBinMaskProcessor::NEW);
	
	force_add(&TestImageGaussian::NEW);
	force_add(&TestImagePureGaussian::NEW);
	force_add(&TestImageSinewave::NEW);
	force_add(&TestImageSinewaveCircular::NEW);
	force_add(&TestImageSquarecube::NEW);
	force_add(&TestImageCirclesphere::NEW);
	force_add(&TestImageAxes::NEW);
	force_add(&TestImageNoiseUniformRand::NEW);
	force_add(&TestImageNoiseGauss::NEW);
	force_add(&TestImageScurve::NEW);
	force_add(&TestImageCylinder::NEW);
	force_add(&TestImageGradient::NEW);
	force_add(&TestTomoImage::NEW);
	
	force_add(&TomoTiltEdgeMaskProcessor::NEW);
	force_add(&TomoTiltAngleWeightProcessor::NEW);
	
	force_add(&NewLowpassTopHatProcessor::NEW);
	force_add(&NewHighpassTopHatProcessor::NEW);
	force_add(&NewBandpassTopHatProcessor::NEW);
	force_add(&NewHomomorphicTopHatProcessor::NEW);
	force_add(&NewLowpassGaussProcessor::NEW);
	force_add(&NewHighpassGaussProcessor::NEW);
	force_add(&NewBandpassGaussProcessor::NEW);
	force_add(&NewHomomorphicGaussProcessor::NEW);
	force_add(&NewInverseGaussProcessor::NEW);
	force_add(&NewLowpassButterworthProcessor::NEW);
	force_add(&NewHighpassButterworthProcessor::NEW);
	force_add(&NewHomomorphicButterworthProcessor::NEW);
	force_add(&NewLowpassTanhProcessor::NEW);
	force_add(&NewHighpassTanhProcessor::NEW);
	force_add(&NewBandpassTanhProcessor::NEW);
	force_add(&NewHomomorphicTanhProcessor::NEW);
	force_add(&NewRadialTableProcessor::NEW);
	force_add(&InverseKaiserI0Processor::NEW);
	force_add(&InverseKaiserSinhProcessor::NEW);	 
	force_add(&CCDNormProcessor::NEW);
	force_add(&CTF_Processor::NEW);
	force_add(&SHIFTProcessor::NEW);

	force_add(&WaveletProcessor::NEW);
	force_add(&FFTProcessor::NEW);
	force_add(&RadialProcessor::NEW);
}

EMData* Processor::process(const EMData * const image)
{
	EMData * result = image->copy();
	process_inplace(result);
	return result;
}

void ImageProcessor::process_inplace(EMData * image)
{
	if (!image) {
		LOGWARN("NULL image");
		return;
	}

	EMData *processor_image = create_processor_image();

	if (image->is_complex()) {
		(*image) *= *processor_image;
	}
	else {
		EMData *fft = image->do_fft();
		(*fft) *= (*processor_image);
		EMData *ift = fft->do_ift();

		float *data = image->get_data();
		float *t = data;
		float *ift_data = ift->get_data();

		data = ift_data;
		ift_data = t;

		ift->update();

		if( fft )
		{
			delete fft;
			fft = 0;
		}

		if( ift )
		{
			delete ift;
			ift = 0;
		}
	}

	image->update();
}

#define FFTRADIALOVERSAMPLE 4
void FourierProcessor::process_inplace(EMData * image)
{
	if (!image) {
		LOGWARN("NULL Image");
		return;
	}

	int array_size = FFTRADIALOVERSAMPLE * image->get_ysize();
	float step=0.5f/array_size;
	
	vector < float >yarray(array_size);

	create_radial_func(yarray);

	if (image->is_complex()) {
		image->apply_radial_func(0, step, yarray);
	}
	else {
		EMData *fft = image->do_fft();
		fft->apply_radial_func(0, step, yarray);
		EMData *ift = fft->do_ift();

		memcpy(image->get_data(),ift->get_data(),ift->get_xsize()*ift->get_ysize()*ift->get_zsize()*sizeof(float));

		ift->update();

		if( fft )
		{
			delete fft;
			fft = 0;
		}

		if( ift )
		{
			delete ift;
			ift = 0;
		}
	}

	image->update();
}

void AmpweightFourierProcessor::process_inplace(EMData * image)
{
	EMData *fft;
	float *fftd;
	int i,f=0;
//	static float sum1=0,sum1a=0;
//	static double sum2=0,sum2a=0;

	if (!image) {
		LOGWARN("NULL Image");
		return;
	}

	if (!image->is_complex()) {
		fft = image->do_fft();
		fftd = fft->get_data();
		f=1;
	}
	else {
		fft=image;
		fftd=image->get_data();
	}
	float *sumd = NULL;
	if (sum) sumd=sum->get_data();
//printf("%d %d    %d %d\n",fft->get_xsize(),fft->get_ysize(),sum->get_xsize(),sum->get_ysize());
	int n = fft->get_xsize()*fft->get_ysize()*fft->get_zsize();
	for (i=0; i<n; i+=2) {
		float c;
		if (dosqrt) c=pow(fftd[i]*fftd[i]+fftd[i+1]*fftd[i+1],0.25f);
		else c = static_cast<float>(hypot(fftd[i],fftd[i+1]));
		if (c==0) c=1.0e-30f;	// prevents divide by zero in normalization
		fftd[i]*=c;
		fftd[i+1]*=c;
		if (sumd) { sumd[i]+=c; sumd[i+1]+=0; }
		
		// debugging
/*		if (i==290*1+12) {
			sum1+=fftd[i];
			sum2+=fftd[i];
			printf("%f\t%f\t%f\t%f\t%f\t%f\n",sum1,sum2,fftd[i],fftd[i+1],sumd[i],c);
		}
		if (i==290*50+60) {
			sum1a+=fftd[i];
			sum2a+=fftd[i];
			printf("%f\t%f\t%f\t%f\t%f\t%f\n",sum1a,sum2a,fftd[i],fftd[i+1],sumd[i],c);
	}*/
	}
	
	if (f) {
		fft->update();
		EMData *ift=fft->do_ift();
		memcpy(image->get_data(),ift->get_data(),n*sizeof(float));
		delete fft;
		delete ift;
	}

	sum->update();
	image->update();

}

void LinearPyramidProcessor::process_inplace(EMData *image) {

	if (image->get_zsize()!=1) { throw ImageDimensionException("Only 2-D images supported"); }

	float *d=image->get_data();
	int nx=image->get_xsize();
	int ny=image->get_ysize();
	
	for (int y=0; y<ny; y++) {
		for (int x=0; x<nx; x++) {
			int l=x+y*nx;
			d[l]*=1.0-abs(x-nx/2)*abs(y-ny/2)*4.0/(nx*ny);
		}
	}
	image->update();
}

EMData * Wiener2DAutoAreaProcessor::process(const EMData * image)
{
// TODO NOT IMPLEMENTED YET !!!
	EMData *ret = 0;
	const EMData *fft;
	float *fftd;
	int f=0;

	if (!image) {
		LOGWARN("NULL Image");
		return ret;
	}
	throw NullPointerException("Processor not yet implemented");

	if (!image->is_complex()) {
		fft = image->do_fft();
		fftd = fft->get_data();
		f=1;
	}
	else {
		fft=image;
		fftd=image->get_data();
	}
	
	return ret;
}

void Wiener2DAutoAreaProcessor::process_inplace(EMData *image) {
	EMData *tmp=process(image);
	memcpy(image->get_data(),tmp->get_data(),image->get_xsize()*image->get_ysize()*image->get_zsize()*sizeof(float));
	delete tmp;
	image->update();
	return;
}


EMData * Wiener2DFourierProcessor::process(const EMData * image)
{
	EMData *ret = NULL;
/*	const EMData *fft;
	float *fftd;
	int f=0;

	if (!image) {
		LOGWARN("NULL Image");
		return ret;
	}

	if (!image->is_complex()) {
		fft = image->do_fft();
		fftd = fft->get_data();
		f=1;
	}
	else {
		fft=image;
		fftd=image->get_data();
	}
	powd=image->get_data();

	int bad=0;
	for (int i=0; i<image->get_xsize()*image->get_ysize(); i+=2) {
		snr=(fftd[i]*fftd[i]+fftd[i+1]*fftd[i+1]-powd[i])/powd[i];
		if (snr<0) { bad++; snr=0; }
		
	}

	print("%d bad pixels\n",snr);
*/	return ret;

}

void Wiener2DFourierProcessor::process_inplace(EMData *image) {
	EMData *tmp=process(image);
	memcpy(image->get_data(),tmp->get_data(),image->get_xsize()*image->get_ysize()*image->get_zsize()*sizeof(float));
	delete tmp;
	image->update();
	return;
}

void LowpassSharpCutoffProcessor::create_radial_func(vector < float >&radial_mask) const
{
	Assert(radial_mask.size() > 0);
	float x = 0.0f , step = 0.5f/radial_mask.size();
	for (size_t i = 0; i < radial_mask.size(); i++) {
		if (x <= lowpass) {
			radial_mask[i] = 1.0f;
		}
		else {
			radial_mask[i] = 0;
		}
		x += step;
	}
}


void HighpassSharpCutoffProcessor::create_radial_func(vector < float >&radial_mask) const
{
	Assert(radial_mask.size() > 0);
	float x = 0.0f , step = 0.5f/radial_mask.size();
	for (size_t i = 0; i < radial_mask.size(); i++) {
		if (x >= highpass) {
			radial_mask[i] = 1.0f;
		}
		else {
			radial_mask[i] = 0;
		}
		x += step;
	}
}


void LowpassGaussProcessor::create_radial_func(vector < float >&radial_mask) const
{
//	printf("rms = %d  lp = %f\n",radial_mask.size(),lowpass);
//	Assert(radial_mask.size() > 0);		// not true, negative numbers do inverse filter processing
	float x = 0.0f , step = 0.5f/radial_mask.size();
	float sig = 1;
	if (lowpass > 0) {
		sig = -1;
	}

	for (size_t i = 0; i < radial_mask.size(); i++) {
		radial_mask[i] = exp(sig * x * x / (lowpass * lowpass));
		x += step;
	}

}

void HighpassGaussProcessor::create_radial_func(vector < float >&radial_mask) const
{
	Assert(radial_mask.size() > 0);
	float x = 0.0f , step = 0.5f/radial_mask.size();
	for (size_t i = 0; i < radial_mask.size(); i++) {
		radial_mask[i] = 1.0f - exp(-x * x / (highpass * highpass));
		x += step;
	}
}


void LowpassTanhProcessor::create_radial_func(vector < float >&radial_mask) const
{
	Assert(radial_mask.size() > 0);
	float x = 0.0f , step = 0.5f/radial_mask.size();
	for (size_t i = 0; i < radial_mask.size(); i++) {
		radial_mask[i] = tanh((lowpass - x)*60.0f) / 2.0f + 0.5f;
		x += step;
	}
}

void HighpassTanhProcessor::create_radial_func(vector < float >&radial_mask) const
{
	Assert(radial_mask.size() > 0);
	float x = 0.0f , step = 0.5f/radial_mask.size();
	for (size_t i = 0; i < radial_mask.size(); i++) {
		radial_mask[i] = tanh(x - highpass) / 2.0f + 0.5f;
		x += step;
	}
}


void HighpassButterworthProcessor::create_radial_func(vector < float >&radial_mask) const
{
	Assert(radial_mask.size() > 0);
	float x = 0.0f , step = 0.5f/radial_mask.size();
	for (size_t i = 0; i < radial_mask.size(); i++) {
		float t = highpass / 1.5f / (x + 0.001f);
		radial_mask[i] = 1.0f / (1.0f + t * t);
		x += step;
	}
}


void LinearRampProcessor::create_radial_func(vector < float >&radial_mask) const
{
	Assert(radial_mask.size() > 0);
	float x = 0.0f , step = 0.5f/radial_mask.size();
	size_t size=radial_mask.size();
	for (size_t i = 0; i < size; i++) {
		radial_mask[i] = intercept + ((slope - intercept) * i) / (size - 1.0f);
		x += step;
	}
}


void RealPixelProcessor::process_inplace(EMData * image)
{
	if (!image) {
		LOGWARN("NULL Image");
		return;
	}

	maxval = image->get_attr("maximum");
	mean = image->get_attr("mean");
	sigma = image->get_attr("sigma");

	calc_locals(image);

	size_t size = (size_t)image->get_xsize() *
		          (size_t)image->get_ysize() *
		          (size_t)image->get_zsize();
	float *data = image->get_data();

	for (size_t i = 0; i < size; i++) {
		process_pixel(&data[i]);
	}
	image->update();
}

void CoordinateProcessor::process_inplace(EMData * image)
{
	if (!image) {
		LOGWARN("NULL Image");
		return;
	}

	maxval = image->get_attr("maximum");
	mean = image->get_attr("mean");
	sigma = image->get_attr("sigma");
	nx = image->get_xsize();
	ny = image->get_ysize();
	nz = image->get_zsize();
	is_complex = image->is_complex();

	calc_locals(image);


	if (!is_valid()) {
		return;
	}

	float *data = image->get_data();
	int i = 0;

	for (int z = 0; z < nz; z++) {
		for (int y = 0; y < ny; y++) {
			for (int x = 0; x < nx; x++) {
				process_pixel(&data[i], x, y, z);
				i++;
			}
		}
	}
	image->update();
}

void PaintProcessor::process_inplace(EMData *image) {
	int nx=image->get_xsize();
	int ny=image->get_ysize();
	int nz=image->get_zsize();

	if (nz==1) {
		for (int j=(y<r2?0:y-r2); j<(y+r2>ny?ny:y+r2); j++) {
			for (int i=(x<r2?0:x-r2); i<(x+r2>nx?nx:x+r2); i++) {
				float r=sqrt(float(Util::square(i-x)+Util::square(j-y)));
				if (r>r2 && r>r1) continue;
				if (r>r1) image->set_value_at(i,j,0,v2*(r-r1)/(r2-r1)+v1*(r2-r)/(r2-r1));
				else image->set_value_at(i,j,0,v1);
			}
		}
	}
	else {
		for (int k=(z<r2?0:z-r2); k<(z+r2>nz?nz:z+r2); k++) {
			for (int j=(y<r2?0:y-r2); j<(y+r2>ny?ny:y+r2); j++) {
				for (int i=(x<r2?0:x-r2); i<(x+r2>nx?nx:x+r2); i++) {
				float r=sqrt(float(Util::square(i-x)+Util::square(j-y)+Util::square(k-z)));
				if (r>r2 && r>r1) continue;
				if (r>r1) image->set_value_at(i,j,k,v2*(r-r1)/(r2-r1)+v1*(r2-r)/(r2-r1));
				else image->set_value_at(i,j,k,v1);
				}
			}
		}
	}
	image->update();
}

void CircularMaskProcessor::calc_locals(EMData *)
{
	xc = nx / 2.0f + dx;
	yc = ny / 2.0f + dy;
	zc = nz / 2.0f + dz;


	if (outer_radius <= 0) {
		outer_radius = nx / 2;
		outer_radius_square = outer_radius * outer_radius;
	}

	if (inner_radius <= 0) {
		inner_radius_square = 0;
	}
}


void MaskEdgeMeanProcessor::calc_locals(EMData * image)
{
	if (!image) {
		throw NullPointerException("NULL image");
	}
	int nitems = 0;
	float sum = 0;
	float *data = image->get_data();
	int i = 0;

	for (int z = 0; z < nz; z++) {
		for (int y = 0; y < ny; y++) {
			for (int x = 0; x < nx; x++) {
				float x1 = sqrt((x - xc) * (x - xc) + (y - yc) * (y - yc) + (z - zc) * (z - zc));
				if (x1 <= outer_radius + ring_width && x1 >= outer_radius - ring_width) {
					sum += data[i];
					nitems++;
				}
				i++;
			}
		}
	}

	ring_avg = sum / nitems;
}

void ComplexPixelProcessor::process_inplace(EMData * image)
{
	if (!image) {
		LOGWARN("NULL image");
		return;
	}
	if (!image->is_complex()) {
		LOGWARN("cannot apply complex processor on a real image. Nothing is done.");
		return;
	}

	size_t size = (size_t)image->get_xsize() *
		          (size_t)image->get_ysize() *
		          (size_t)image->get_zsize();
	float *data = image->get_data();

	image->ri2ap();

	for (size_t i = 0; i < size; i += 2) {
		process_pixel(data);
		data += 2;
	}

	image->update();
	image->ap2ri();
}



void AreaProcessor::process_inplace(EMData * image)
{
	if (!image) {
		LOGWARN("NULL Image");
		return;
	}

	float *data = image->get_data();

	nx = image->get_xsize();
	ny = image->get_ysize();
	nz = image->get_zsize();

	int n = (areasize - 1) / 2;
	matrix_size = areasize * areasize;

	if (nz > 1) {
		matrix_size *= areasize;
	}

	float *matrix = new float[matrix_size];
	kernel = new float[matrix_size];

	int cpysize = areasize * (int) sizeof(float);
	int start = (nx * ny + nx + 1) * n;

	int xend = nx - n;
	int yend = ny - n;

	int zstart = n;
	int zend = nz - n;

	int zbox_start = 0;
	int zbox_end = areasize;

	if (nz == 1) {
		zstart = 0;
		zend = 1;
		zbox_end = 1;
	}

	size_t nsec = (size_t)nx * (size_t)ny;
	int box_nsec = areasize * areasize;

	create_kernel();

	size_t total_size = (size_t)nx * (size_t)ny * (size_t)nz;
	float *data2 = new float[total_size];
	memcpy(data2, data, total_size * sizeof(float));


	for (int z = zstart; z < zend; z++) {
		for (int y = n; y < yend; y++) {
			for (int x = n; x < xend; x++) {

				size_t k = z * nsec + y * nx + x;

				for (int bz = zbox_start; bz < zbox_end; bz++) {
					for (int by = 0; by < areasize; by++) {
						memcpy(&matrix[bz * box_nsec + by * areasize],
							   &data2[k - start + bz * nsec + by * nx], cpysize);
					}
				}

				process_pixel(&data[k], (float) x, (float) y, (float) z, matrix);
			}
		}
	}

	if( matrix )
	{
		delete[]matrix;
		matrix = 0;
	}
	
	if( kernel )
	{
		delete[]kernel;
		kernel = 0;
	}
	image->update();
}


void LaplacianProcessor::create_kernel() const
{
	if (nz == 1) {
		memset(kernel, 0, areasize * areasize);
		kernel[1] = -0.25f;
		kernel[3] = -0.25f;
		kernel[5] = -0.25f;
		kernel[7] = -0.25f;
		kernel[4] = 1;
	}
	else {
		memset(kernel, 0, areasize * areasize * areasize);
		kernel[4] = -1.0f / 6.0f;
		kernel[10] = -1.0f / 6.0f;
		kernel[12] = -1.0f / 6.0f;
		kernel[14] = -1.0f / 6.0f;
		kernel[16] = -1.0f / 6.0f;
		kernel[22] = -1.0f / 6.0f;
		kernel[13] = 1;
	}
}

void BoxStatProcessor::process_inplace(EMData * image)
{
	if (!image) {
		LOGWARN("NULL Image");
		return;
	}

	int nx = image->get_xsize();
	int ny = image->get_ysize();
	int nz = image->get_zsize();

	int n = params.set_default("radius",1);
	int areasize = 2 * n + 1;

	int matrix_size = areasize * areasize;
	if (nz > 1) {
		matrix_size *= areasize;
	}

	float *array = new float[matrix_size];
//	image->process_inplace("normalize");

	float *data = image->get_data();
	size_t total_size = (size_t)nx * (size_t)ny * (size_t)nz;
	float *data2 = new float[total_size];
	memcpy(data2, data, total_size * sizeof(float));

	int z_begin = 0;
	int z_end = 1;
	int nzz=0;
	if (nz > 1) {
		z_begin = n;
		z_end = nz - n;
		nzz=n;
	}

	int nxy = nx * ny;

	for (int k = z_begin; k < z_end; k++) {
		int knxy = k * nxy;

		for (int j = n; j < ny - n; j++) {
			int jnx = j * nx;

			for (int i = n; i < nx - n; i++) {
				int s = 0;

				for (int i2 = i - n; i2 <= i + n; i2++) {
					for (int j2 = j - n; j2 <= j + n; j2++) {
						for (int k2 = k - nzz; k2 <= k + nzz; k2++) {
							array[s] = data2[i2 + j2 * nx + k2 * nxy];
							s++;
						}
					}
				}

				process_pixel(&data[i + jnx + knxy], array, matrix_size);
			}
		}
	}

	image->update();

	if( data2 )
	{
		delete[]data2;
		data2 = 0;
	}
}

void DiffBlockProcessor::process_inplace(EMData * image)
{
	if (!image) {
		LOGWARN("NULL Image");
		return;
	}

	int nz = image->get_zsize();

	if (nz > 1) {
		LOGERR("%s Processor doesn't support 3D", get_name().c_str());
		throw ImageDimensionException("3D model not supported");
	}

	int nx = image->get_xsize();
	int ny = image->get_ysize();

	int v1 = params["cal_half_width"];
	int v2 = params["fill_half_width"];

	int v0 = v1 > v2 ? v1 : v2;

	if (v2 <= 0) {
		v2 = v1;
	}

	float *data = image->get_data();

	for (int y = v0; y <= ny - v0 - 1; y += v2) {
		for (int x = v0; x <= nx - v0 - 1; x += v2) {

			float sum = 0;
			for (int y1 = y - v1; y1 <= y + v1; y1++) {
				for (int x1 = x - v1; x1 <= x + v1; x1++) {
					sum += data[x1 + y1 * nx];
				}
			}
			float mean = sum / ((v1 * 2 + 1) * (v1 * 2 + 1));

			for (int j = y - v2; j <= y + v2; j++) {
				for (int i = x - v2; i <= x + v2; i++) {
					data[i + j * nx] = mean;
				}
			}
		}
	}

	image->update();
}


void CutoffBlockProcessor::process_inplace(EMData * image)
{
	if (!image) {
		LOGWARN("NULL Image");
		return;
	}
	int nz = image->get_zsize();

	if (nz > 1) {
		LOGERR("%s Processor doesn't support 3D", get_name().c_str());
		throw ImageDimensionException("3D model not supported");
	}

	int nx = image->get_xsize();
	int ny = image->get_ysize();

	float value1 = params["value1"];
	float value2 = params["value2"];

	int v1 = (int) value1;
	int v2 = (int) value2;
	if (v2 > v1 / 2) {
		LOGERR("invalid value2 '%f' in CutoffBlockProcessor", value2);
		return;
	}

	if (v2 <= 0) {
		v2 = v1;
	}

	float *data = image->get_data();
	int y = 0, x = 0;
	for (y = 0; y <= ny - v1; y += v1) {
		for (x = 0; x <= nx - v1; x += v1) {

			EMData *clip = image->get_clip(Region(x, y, v1, v1));
			EMData *fft = clip->do_fft();

			float *fft_data = fft->get_data();
			float sum = 0;
			int nitems = 0;

			for (int i = -v2; i < v2; i++) {
				for (int j = 0; j < v2; j++) {
					if (j == 0 && i == 0) {
						continue;
					}

					if (hypot(j, i) < value2) {
						int t = j * 2 + (i + v1 / 2) * (v1 + 2);
						sum += (fft_data[t] * fft_data[t] + fft_data[t + 1] * fft_data[t + 1]);
						nitems++;
					}
				}
			}

			if( clip )
			{
				delete clip;
				clip = 0;
			}

			float mean = sum / nitems;

			for (int i = y; i < y + v1; i++) {
				for (int j = x; j < x + v1; j++) {
					data[i * nx + j] = mean;
				}
			}
		}
	}

	memset(&data[y * nx], 0, (ny - y) * nx * sizeof(float));

	for (int i = 0; i < ny; i++) {
		memset(&data[i * nx + x], 0, (nx - x) * sizeof(float));
	}

	image->update();
}


template<class LogicOp>
EMData* BooleanShrinkProcessor::process(const EMData *const image)
{
	// The basic idea of this code is to iterate through each pixel in the output image
	// determining its value by investigation a region of the input image
	
	if (!image) throw NullPointerException("Attempt to max shrink a null image");

	if (image->is_complex() ) throw ImageFormatException("Can not max shrink a complex image");
	
	
	int shrink = params.set_default("shrink",2);
	int search = params.set_default("search",2);
	
	if ( shrink < 0 ) throw InvalidValueException(shrink, "Can not shrink by a value less than 0");
	

	int nz = image->get_zsize();
	int ny = image->get_ysize();
	int nx = image->get_xsize();
	
	if (nx == 1 && ny == 1 && nz == 1 ) return image->copy();
	
	LogicOp op;
	EMData* return_image = new EMData();
	
	int shrinkx = shrink;
	int shrinky = shrink;
	int shrinkz = shrink;
	
	int searchx = search;
	int searchy = search;
	int searchz = search;
		
	// Clamping the shrink values to the dimension lengths
	// ensures that the return image has non zero dimensions
	if ( shrinkx > nx ) shrinkx = nx;
	if ( shrinky > ny ) shrinky = ny;
	if ( shrinkz > nz ) shrinkz = nz;
	
	if ( nz == 1 && ny == 1 )
	{
		return_image->set_size(nx/shrinkx);
		for(int i = 0; i < nx/shrinkx; ++i)
		{
			float tmp = op.get_start_val();
			for(int s=0; s < searchx; ++s)
			{
				int idx = shrinkx*i+s;
				// Don't ask for memory beyond limits
				if ( idx > nx ) break;
				else
				{
					float val = image->get_value_at(idx);
					if ( op( val,tmp) ) tmp = val;
				}	
			}
			return_image->set_value_at(i,tmp);
		}
	}
	else if ( nz == 1 )
	{
		int ty = ny/shrinky;
		int tx = nx/shrinkx;
		return_image->set_size(tx,ty);
		for(int y = 0; y < ty; ++y) {
			for(int x = 0; x < tx; ++x) {
				float tmp = op.get_start_val();
				for(int sy=0; sy < searchy; ++sy) {
					int yidx = shrinky*y+sy;
					if ( yidx >= ny) break;
					for(int sx=0; sx < searchx; ++sx) {
						int xidx = shrinkx*x+sx;
						if ( xidx >= nx) break;
					
						float val = image->get_value_at(xidx,yidx);
						if ( op( val,tmp) ) tmp = val;
					}
				}
				return_image->set_value_at(x,y,tmp);
			}
		}
	}
	else
	{
		int tz = nz/shrinkz;
		int ty = ny/shrinky;
		int tx = nx/shrinkx;
		
		return_image->set_size(tx,ty,tz);
		for(int z = 0; z < tz; ++z) {
			for(int y = 0; y < ty; ++y) {
				for(int x = 0; x < tx; ++x) {
					float tmp = op.get_start_val();
					
					for(int sz=0; sz < searchz; ++sz) {
						int zidx = shrinkz*z+sz;
						if ( zidx >= nz) break;
						
						for(int sy=0; sy < searchy; ++sy) {
							int yidx = shrinky*y+sy;
							if ( yidx >= ny) break;
							
							for(int sx=0; sx < searchx; ++sx) {
								int xidx = shrinkx*x+sx;
								if ( xidx >= nx) break;
								float val = image->get_value_at(xidx,yidx,zidx);
								if ( op( val,tmp) ) tmp = val;
							}
						}
					}
					return_image->set_value_at(x,y,z,tmp);
				}
			}
		}
	}
	return_image->update();
	
	return return_image;
}

template<class LogicOp>
void BooleanShrinkProcessor::process_inplace(EMData * image)
{
	// The basic idea of this code is to iterate through each pixel in the output image
	// determining its value by investigation a region of the input image
	if (!image) throw NullPointerException("Attempt to max shrink a null image");

	if (image->is_complex() ) throw ImageFormatException("Can not max shrink a complex image");
	
	
	int shrink = params.set_default("shrink",2);
	int search = params.set_default("search",2);
	
	if ( shrink < 0 ) throw InvalidValueException(shrink, "Can not shrink by a value less than 0");
	

	int nz = image->get_zsize();
	int ny = image->get_ysize();
	int nx = image->get_xsize();
	
	LogicOp op;
	
	int shrinkx = shrink;
	int shrinky = shrink;
	int shrinkz = shrink;
	
	int searchx = search;
	int searchy = search;
	int searchz = search;
		
	// Clamping the shrink values to the dimension lengths
	// ensures that the return image has non zero dimensions
	if ( shrinkx > nx ) shrinkx = nx;
	if ( shrinky > ny ) shrinkx = ny;
	if ( shrinkz > nz ) shrinkx = nz;
	
	if (nx == 1 && ny == 1 && nz == 1 ) return;
	
	if ( nz == 1 && ny == 1 )
	{
		for(int i = 0; i < nx/shrink; ++i)
		{
			float tmp = op.get_start_val();
			for(int s=0; s < searchx; ++s)
			{
				int idx = shrinkx*i+s;
				if ( idx > nx ) break;
				else
				{
					float val = image->get_value_at(idx);
					if ( op( val,tmp) ) tmp = val;
				}	
			}
			image->set_value_at(i,tmp);
		}
		
		image->set_size(nx/shrinkx);
	}
	else if ( nz == 1 )
	{
		int ty = ny/shrinky;
		int tx = nx/shrinkx;
		for(int y = 0; y < ty; ++y) {
			for(int x = 0; x < tx; ++x) {
				float tmp = op.get_start_val();
				for(int sy=0; sy < searchy; ++sy) {
					int yidx = shrinky*y+sy;
					if ( yidx >= ny) break;
					for(int sx=0; sx < searchx; ++sx) {
						int xidx = shrinkx*x+sx;
						if ( xidx >= nx) break;
						
						float val = image->get_value_at(xidx,yidx);
						if ( op( val,tmp) ) tmp = val;
					}
				}
				(*image)(x+tx*y) = tmp;
			}
		}
		image->set_size(tx,ty);
	}
	else
	{	
		int tnxy = nx/shrinkx*ny/shrinky;
		int tz = nz/shrinkz;
		int ty = ny/shrinky;
		int tx = nx/shrinkx;
		
		for(int z = 0; z < tz; ++z) {
			for(int y = 0; y < ty; ++y) {
				for(int x = 0; x < tx; ++x) {
					float tmp = op.get_start_val();
					for(int sz=0; sz < searchz; ++sz) {
						int zidx = shrinkz*z+sz;
						if ( zidx >= nz) break;
						for(int sy=0; sy < searchy; ++sy) {
							int yidx = shrinky*y+sy;
							if ( yidx >= ny) break;
							for(int sx=0; sx < shrinkx; ++sx) {
								int xidx = shrinkx*x+sx;
								if ( xidx >= nx) break;
								
								float val = image->get_value_at(xidx,yidx,zidx);
								if ( op( val,tmp) ) tmp = val;
							}
						}
					}
					(*image)(x+tx*y+tnxy*z) = tmp;
				}
			}
		}
		image->set_size(tx,ty,tz);
	}
	
	image->update();
}



void GradientRemoverProcessor::process_inplace(EMData * image)
{
	if (!image) {
		LOGWARN("NULL Image");
		return;
	}

	int nz = image->get_zsize();
	if (nz > 1) {
		LOGERR("%s Processor doesn't support 3D model", get_name().c_str());
		throw ImageDimensionException("3D model not supported");
	}

	int nx = image->get_xsize();
	int ny = image->get_ysize();
	float *dy = new float[ny];
	float m = 0;
	float b = 0;
	float sum_y = 0;
	float *data = image->get_data();

	for (int i = 0; i < ny; i++) {
		Util::calc_least_square_fit(nx, 0, data + i * nx, &m, &b, false);
		dy[i] = b;
		sum_y += m;
	}

	float mean_y = sum_y / ny;
	float sum_x = 0;
	Util::calc_least_square_fit(ny, 0, dy, &sum_x, &b, false);

	for (int j = 0; j < ny; j++) {
		for (int i = 0; i < nx; i++) {
			data[i + j * nx] -= i * sum_x + j * mean_y + b;
		}
	}

	image->update();
}

#include <gsl/gsl_linalg.h>
void GradientPlaneRemoverProcessor::process_inplace(EMData * image)
{
	if (!image) {
		LOGWARN("NULL Image");
		return;
	}

	int nz = image->get_zsize();
	if (nz > 1) {
		LOGERR("%s Processor doesn't support 3D model", get_name().c_str());
		throw ImageDimensionException("3D map not supported");
	}

	int nx = image->get_xsize();
	int ny = image->get_ysize();
	float *d = image->get_data();
	EMData *mask = 0;
	float *dm = 0;
	if (params.has_key("mask")) {
		mask = params["mask"];
		if (nx!=mask->get_xsize() || ny!=mask->get_ysize()) {
			LOGERR("%s Processor requires same size mask image", get_name().c_str());
			throw ImageDimensionException("wrong size mask image");
		}
		dm = mask->get_data();
	}
	int count = 0;
	if (dm) {
		for(int i=0; i<nx*ny; i++) {
			if(dm[i]) count++;
		}
	}
	else {
		count = nx * ny;
	}
	if(count<3) {
		LOGERR("%s Processor requires at least 3 pixels to fit a plane", get_name().c_str());
		throw ImageDimensionException("too few usable pixels to fit a plane");
	}
	// Allocate the working space
	gsl_vector *S=gsl_vector_calloc(3);
	gsl_matrix *A=gsl_matrix_calloc(count,3);
	gsl_matrix *V=gsl_matrix_calloc(3,3);
	
	double m[3] = {0, 0, 0};
	int index=0;
	if (dm) {
		for(int j=0; j<ny; j++){
			for(int i=0; i<nx; i++){
				int ij=j*nx+i;
				if(dm[ij]) {
					m[0]+=i;	// x
					m[1]+=j;	// y
					m[2]+=d[ij];	// z
					/*printf("index=%d/%d\ti,j=%d,%d\tval=%g\txm,ym,zm=%g,%g,%g\n", \
					        index,count,i,j,d[ij],m[0]/(index+1),m[1]/(index+1),m[2]/(index+1));*/
					index++;
				}
			}
		}
	}
	else {
		for(int j=0; j<ny; j++){
			for(int i=0; i<nx; i++){
				int ij=j*nx+i;
					m[0]+=i;	// x
					m[1]+=j;	// y
					m[2]+=d[ij];	// z
					/*printf("index=%d/%d\ti,j=%d,%d\tval=%g\txm,ym,zm=%g,%g,%g\n", \
					        index,count,i,j,d[ij],m[0]/(index+1),m[1]/(index+1),m[2]/(index+1));*/
					index++;
			}
		}
	}
	
	for(int i=0; i<3; i++) m[i]/=count;	// compute center of the plane
	
	index=0;
	if (dm) {
		for(int j=0; j<ny; j++){
			for(int i=0; i<nx; i++){
				int ij=j*nx+i;
				if(dm[ij]) {
					//printf("index=%d/%d\ti,j=%d,%d\tval=%g\n",index,count,i,j,d[index]);
					gsl_matrix_set(A,index,0,i-m[0]);
					gsl_matrix_set(A,index,1,j-m[1]);
					gsl_matrix_set(A,index,2,d[ij]-m[2]);
					index++;
				}
			}
		}
		mask->update(); 
	}
	else {
		for(int j=0; j<ny; j++){
			for(int i=0; i<nx; i++){
				int ij=j*nx+i;
					//printf("index=%d/%d\ti,j=%d,%d\tval=%g\n",index,count,i,j,d[index]);
					gsl_matrix_set(A,index,0,i-m[0]);
					gsl_matrix_set(A,index,1,j-m[1]);
					gsl_matrix_set(A,index,2,d[ij]-m[2]);
					index++;
			}
		}
	}
	
	// SVD decomposition and use the V vector associated with smallest singular value as the plan normal
	gsl_linalg_SV_decomp_jacobi(A, V, S);
	
	double n[3];
	for(int i=0; i<3; i++) n[i] = gsl_matrix_get(V, i, 2);
	
	#ifdef DEBUG
	printf("S=%g,%g,%g\n",gsl_vector_get(S,0), gsl_vector_get(S,1), gsl_vector_get(S,2));
	printf("V[0,:]=%g,%g,%g\n",gsl_matrix_get(V,0,0), gsl_matrix_get(V,0,1),gsl_matrix_get(V,0,2));
	printf("V[1,:]=%g,%g,%g\n",gsl_matrix_get(V,1,0), gsl_matrix_get(V,1,1),gsl_matrix_get(V,1,2));
	printf("V[2,:]=%g,%g,%g\n",gsl_matrix_get(V,2,0), gsl_matrix_get(V,2,1),gsl_matrix_get(V,2,2));
	printf("Fitted plane: p0=%g,%g,%g\tn=%g,%g,%g\n",m[0],m[1],m[2],n[0],n[1],n[2]);
	#endif
	
	int changeZero = 0;
	if (params.has_key("changeZero")) changeZero = params["changeZero"];
	if (changeZero) {
		for(int j=0; j<nx; j++){
			for(int i=0; i<ny; i++){
				int ij = j*nx+i;
				d[ij]-=static_cast<float>(-((i-m[0])*n[0]+(j-m[1])*n[1])/n[2]+m[2]);
			}
		}
	}
	else {
		for(int j=0; j<nx; j++){
			for(int i=0; i<ny; i++){
				int ij = j*nx+i;
				if(d[ij]) d[ij]-=static_cast<float>(-((i-m[0])*n[0]+(j-m[1])*n[1])/n[2]+m[2]);
			}
		}
	}
	image->update();
	// set return plane parameters
	vector< float > planeParam;
	planeParam.resize(6);
	for(int i=0; i<3; i++) planeParam[i] = static_cast<float>(n[i]);
	for(int i=0; i<3; i++) planeParam[i+3] = static_cast<float>(m[i]);
	params["planeParam"]=EMObject(planeParam);
}

void VerticalStripeProcessor::process_inplace(EMData * image)
{
	if (!image) {
		LOGWARN("NULL Image");
		return;
	}

	int nx = image->get_xsize();
	int ny = image->get_ysize();
	int nz = image->get_zsize();

	float *data = image->get_data();
	float sigma = image->get_attr("sigma");

	for (int k = 0; k < nz; k++) {
		for (int i = 0; i < nx; i++) {
			double sum = 0;
			for (int j = ny / 4; j < 3 * ny / 4; j++) {
				sum += data[i + j * nx];
			}

			float mean = (float)sum / (ny / 2);
			for (int j = 0; j < ny; j++) {
				data[i + j * nx] = (data[i + j * nx] - mean) / sigma;
			}
		}
	}

	image->update();
}

void RealToFFTProcessor::process_inplace(EMData *image)
{
	if (!image) {
		LOGWARN("NULL Image");
		return;
	}
	
	//Note : real image only!
	if(image->is_complex()) {
		LOGERR("%s Processor only operates on real images", get_name().c_str());
		throw ImageFormatException("apply to real image only");
	}

	// Note : 2D only!	
	int nz = image->get_zsize();
	if (nz > 1) {
		LOGERR("%s Processor doesn't support 3D models", get_name().c_str());
		throw ImageDimensionException("3D model not supported");
	}
	
	EMData *ff=image->do_fft();
	ff->ri2ap();
	
	int nx=image->get_xsize();
	int ny=image->get_ysize();
	
	int x,y;
	float norm=static_cast<float>(nx*ny);

	for (y=0; y<ny; y++) image->set_value_at(0,y,0);
	
	for (x=1; x<nx/2; x++) {
		for (y=0; y<ny; y++) {
			int y2;
			if (y<ny/2) y2=y+ny/2;
			else if (y==ny/2) y2=ny;
			else y2=y-ny/2;
			image->set_value_at(x,y,ff->get_value_at(nx-x*2,ny-y2)/norm);
		}
	}
	
	for (x=nx/2; x<nx; x++) {
		for (y=0; y<ny; y++) {
			int y2;
			if (y<ny/2) y2=y+ny/2;
			else y2=y-ny/2;
			image->set_value_at(x,y,ff->get_value_at(x*2-nx,y2)/norm);
		}
	}
	
	image->update();
	if( ff )
	{
		delete ff;
		ff = 0;
	}
}

void SigmaZeroEdgeProcessor::process_inplace(EMData * image)
{
	if (!image) {
		LOGWARN("NULL Image");
		return;
	}

	if (image->get_zsize() > 1) {
		LOGERR("%s Processor doesn't support 3D model", get_name().c_str());
		throw ImageDimensionException("3D model not supported");
	}
	float *d = image->get_data();
	int i = 0;
	int j = 0;

	int nx = image->get_xsize();
	int ny = image->get_ysize();

	for (j = 0; j < ny; j++) {
		for (i = 0; i < nx - 1; i++) {
			if (d[i + j * nx] != 0) {
				break;
			}
		}

		float v = d[i + j * nx];
		while (i >= 0) {
			d[i + j * nx] = v;
			i--;
		}

		for (i = nx - 1; i > 0; i--) {
			if (d[i + j * nx] != 0)
				break;
		}
		v = d[i + j * nx];
		while (i < nx) {
			d[i + j * nx] = v;
			i++;
		}
	}

	for (i = 0; i < nx; i++) {
		for (j = 0; j < ny; j++) {
			if (d[i + j * nx] != 0)
				break;
		}

		float v = d[i + j * nx];
		while (j >= 0) {
			d[i + j * nx] = v;
			j--;
		}

		for (j = ny - 1; j > 0; j--) {
			if (d[i + j * nx] != 0)
				break;
		}
		v = d[i + j * nx];
		while (j < ny) {
			d[i + j * nx] = v;
			j++;
		}
	}


	image->update();
}



void BeamstopProcessor::process_inplace(EMData * image)
{
	if (!image) {
		LOGWARN("NULL Image");
		return;
	}
	if (image->get_zsize() > 1) {
		LOGERR("BeamstopProcessor doesn't support 3D model");
		throw ImageDimensionException("3D model not supported");
	}

	float value1 = params["value1"];
	float value2 = params["value2"];
	float value3 = params["value3"];

	float thr = fabs(value1);
	float *data = image->get_data();
	int cenx = (int) value2;
	int ceny = (int) value3;

	int nx = image->get_xsize();
	int ny = image->get_ysize();

	if (cenx <= 0) {
		cenx = nx / 2;
	}

	if (ceny <= 0) {
		ceny = ny / 2;
	}

	int mxr = (int) floor(sqrt(2.0f) * nx / 2);

	float *mean_values = new float[mxr];
	float *sigma_values = new float[mxr];
	double sum = 0;
	int count = 0;
	double square_sum = 0;

	for (int i = 0; i < mxr; i++) {
		sum = 0;
		count = 0;
		square_sum = 0;
		int nitems = 6 * i + 2;

		for (int j = 0; j < nitems; j++) {
			float ang = j * 2 * M_PI / nitems;
			int x0 = (int) floor(cos(ang) * i + cenx);
			int y0 = (int) floor(sin(ang) * i + ceny);

			if (x0 < 0 || y0 < 0 || x0 >= nx || y0 >= ny) {
				continue;
			}

			float f = data[x0 + y0 * nx];
			sum += f;
			square_sum += f * f;
			count++;
		}

		mean_values[i] = (float)sum / count;
		sigma_values[i] = (float) sqrt(square_sum / count - mean_values[i] * mean_values[i]);
	}


	for (int k = 0; k < 5; k++) {
		for (int i = 0; i < mxr; i++) {
			sum = 0;
			count = 0;
			square_sum = 0;
			int nitems = 6 * i + 2;
			double thr1 = mean_values[i] - sigma_values[i] * thr;
			double thr2 = mean_values[i] + sigma_values[i];

			for (int j = 0; j < nitems; j++) {
				float ang = j * 2 * M_PI / nitems;
				int x0 = (int) floor(cos(ang) * i + cenx);
				int y0 = (int) floor(sin(ang) * i + ceny);

				if (x0 < 0 || y0 < 0 || x0 >= nx || y0 >= ny ||
					data[x0 + y0 * nx] < thr1 || data[x0 + y0 * nx] > thr2) {
					continue;
				}

				sum += data[x0 + y0 * nx];
				square_sum += data[x0 + y0 * nx] * data[x0 + y0 * nx];
				count++;
			}

			mean_values[i] = (float) sum / count;
			sigma_values[i] = (float) sqrt(square_sum / count - mean_values[i] * mean_values[i]);
		}
	}

	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny; j++) {

			int r = Util::round(hypot((float) i - cenx, (float) j - ceny));

			if (value1 < 0) {
				if (data[i + j * nx] < (mean_values[r] - sigma_values[r] * thr)) {
					data[i + j * nx] = 0;
				}
				else {
					data[i + j * nx] -= mean_values[r];
				}
				continue;
			}
			if (data[i + j * nx] > (mean_values[r] - sigma_values[r] * thr)) {
				continue;
			}
			data[i + j * nx] = mean_values[r];
		}
	}

	if( mean_values )
	{
		delete[]mean_values;
		mean_values = 0;
	}

	if( sigma_values )
	{
		delete[]sigma_values;
		sigma_values = 0;
	}

	image->update();
}



void MeanZeroEdgeProcessor::process_inplace(EMData * image)
{
	if (!image) {
		LOGWARN("NULL Image");
		return;
	}
	if (image->get_zsize() > 1) {
		LOGERR("MeanZeroEdgeProcessor doesn't support 3D model");
		throw ImageDimensionException("3D model not supported");
	}

	int nx = image->get_xsize();
	int ny = image->get_ysize();
	Dict dict = image->get_attr_dict();
	float mean_nonzero = dict.get("mean_nonzero");

	float *d = image->get_data();
	int i = 0;
	int j = 0;

	for (j = 0; j < ny; j++) {
		for (i = 0; i < nx - 1; i++) {
			if (d[i + j * nx] != 0) {
				break;
			}
		}

		if (i == nx - 1) {
			i = -1;
		}

		float v = d[i + j * nx] - mean_nonzero;

		while (i >= 0) {
			v *= 0.9f;
			d[i + j * nx] = v + mean_nonzero;
			i--;
		}


		for (i = nx - 1; i > 0; i--) {
			if (d[i + j * nx] != 0) {
				break;
			}
		}

		if (i == 0) {
			i = nx;
		}

		v = d[i + j * nx] - mean_nonzero;

		while (i < nx) {
			v *= .9f;
			d[i + j * nx] = v + mean_nonzero;
			i++;
		}
	}


	for (i = 0; i < nx; i++) {
		for (j = 0; j < ny; j++) {
			if (d[i + j * nx] != 0)
				break;
		}

		float v = d[i + j * nx] - mean_nonzero;

		while (j >= 0) {
			v *= .9f;
			d[i + j * nx] = v + mean_nonzero;
			j--;
		}

		for (j = ny - 1; j > 0; j--) {
			if (d[i + j * nx] != 0)
				break;
		}

		v = d[i + j * nx] - mean_nonzero;

		while (j < ny) {
			v *= .9f;
			d[i + j * nx] = v + mean_nonzero;
			j++;
		}
	}

	image->update();
}



void AverageXProcessor::process_inplace(EMData * image)
{
	if (!image) {
		LOGWARN("NULL Image");
		return;
	}

	float *data = image->get_data();
	int nx = image->get_xsize();
	int ny = image->get_ysize();
	int nz = image->get_zsize();
	size_t nxy = (size_t)nx * ny;

	for (int z = 0; z < nz; z++) {
		for (int x = 0; x < nx; x++) {
			double sum = 0;
			for (int y = 0; y < ny; y++) {
				sum += data[x + y * nx + z * nxy];
			}
			float mean = (float) sum / ny;

			for (int y = 0; y < ny; y++) {
				data[x + y * nx + z * nxy] = mean;
			}
		}
	}

	image->update();
}


void ZeroEdgeRowProcessor::process_inplace(EMData * image)
{
	if (!image) {
		LOGWARN("NULL Image");
		return;
	}

	if (image->get_zsize() > 1) {
		LOGERR("ZeroEdgeRowProcessor is not supported in 3D models");
		throw ImageDimensionException("3D model not supported");
	}

	int nx = image->get_xsize();
	int ny = image->get_ysize();

	float *d = image->get_data();
	int top_nrows = params["y0"];
	int bottom_nrows = params["y1"];

	int left_ncols = params["x0"];
	int right_ncols = params["x1"];

	size_t row_size = nx * sizeof(float);

	memset(d, 0, top_nrows * row_size);
	memset(d + (ny - bottom_nrows) * nx, 0, bottom_nrows * row_size);

	for (int i = top_nrows; i < ny - bottom_nrows; i++) {
		memset(d + i * nx, 0, left_ncols * sizeof(float));
		memset(d + i * nx + nx - right_ncols, 0, right_ncols * sizeof(float));
	}
	image->update();
}

void ZeroEdgePlaneProcessor::process_inplace(EMData * image)
{
	if (!image) {
		LOGWARN("NULL Image");
		return;
	}

	if (image->get_zsize() <= 1) {
		LOGERR("ZeroEdgePlaneProcessor only support 3D models");
		throw ImageDimensionException("3D model only");
	}

	int nx = image->get_xsize();
	int ny = image->get_ysize();
	int nz = image->get_zsize();

	float *d = image->get_data();
	
	int x0=params["x0"];
	int x1=params["x1"];
	int y0=params["y0"];
	int y1=params["y1"];
	int z0=params["z0"];
	int z1=params["z1"];
	
	size_t row_size = nx * sizeof(float);
	size_t nxy = nx * ny;
	size_t sec_size = nxy * sizeof(float);
	size_t y0row = y0 * row_size;
	size_t y1row = y1 * row_size;
	int max_y = ny-y1;
	size_t x0size = x0*sizeof(float);
	size_t x1size = x1*sizeof(float);
	
	memset(d,0,z0*sec_size);					// zero -z
	memset(d+(nxy*(nz-z1)),0,sec_size*z1);	    // zero +z
	
	for (int z=z0; z<nz-z1; z++) {
		memset(d+z*nxy,0,y0row);			// zero -y
		memset(d+z*nxy+(ny-y1)*nx,0,y1row);	// zero +y

		int znxy = z * nxy;
		int znxy2 = znxy + nx - x1;
		
		for (int y=y0; y<max_y; y++) {
			memset(d+znxy+y*nx,0,x0size);	// zero -x
			memset(d+znxy2+y*nx,0,x1size);	// zero +x
		}
	}
	
	image->update();
}


float NormalizeProcessor::calc_sigma(EMData * image) const
{
	return image->get_attr("sigma");
}

void NormalizeProcessor::process_inplace(EMData * image)
{
	if (!image) {
		LOGWARN("cannot do normalization on NULL image");
		return;
	}

	if (image->is_complex()) {
		LOGWARN("cannot do normalization on complex image");
		return;
	}

	float sigma = calc_sigma(image);
	if (sigma == 0 || !Util::goodf(&sigma)) {
		LOGWARN("cannot do normalization on image with sigma = 0");
		return;
	}

	float mean = calc_mean(image);

	size_t size = (size_t)image->get_xsize() * (size_t)image->get_ysize() *
		          (size_t)image->get_zsize();
	float *data = image->get_data();

	for (size_t i = 0; i < size; i++) {
		data[i] = (data[i] - mean) / sigma;
	}

	image->update();
}

float NormalizeUnitProcessor::calc_sigma(EMData * image) const
{
	if (!image) {
		LOGWARN("NULL Image");
		return 0;
	}
	float ret=sqrt((float)image->get_attr("square_sum"));
	return ret==0.0f?1.0f:ret;
}

float NormalizeUnitSumProcessor::calc_sigma(EMData * image) const
{
	if (!image) {
		LOGWARN("NULL Image");
		return 0;
	}
	float ret=(float)image->get_attr("mean")*image->get_xsize()*image->get_ysize()*image->get_zsize();
	return ret==0.0f?1.0f:ret;
}

float NormalizeMaskProcessor::calc_sigma(EMData * image) const
{
	if (!image) {
		LOGWARN("NULL Image");
		return 0;
	}
	EMData *mask = params["mask"];
	int no_sigma = params["no_sigma"];
	
	if(no_sigma == 0) {
		return 1;
	}
	else {
		if (!EMUtil::is_same_size(mask, image)) {
			LOGERR("normalize.maskProcessor: mask and image must be the same size");
			throw ImageDimensionException("mask and image must be the same size");
		}
	
		float *data = image->get_data();
		float *mask_data = mask->get_data();
		size_t size = image->get_xsize() * image->get_ysize() * image->get_zsize();
		double sum = 0;
		double sq2 = 0;
		size_t n_norm = 0;
	
		for (size_t i = 0; i < size; i++) {
			if (mask_data[i] > 0.5f) {
				sum += data[i];
				sq2 += data[i]*double (data[i]);
				n_norm++;
			}
		}
		return sqrt(static_cast<float>((sq2 - sum * sum /n_norm)/(n_norm -1))) ;
	}
}

float NormalizeMaskProcessor::calc_mean(EMData * image) const
{
	if (!image) {
		LOGWARN("NULL Image");
		return 0;
	}
	EMData *mask = params["mask"];

	if (!EMUtil::is_same_size(mask, image)) {
		LOGERR("normalize.maskProcessor: mask and image must be the same size");
		throw ImageDimensionException("mask and image must be the same size");
	}

	float *data = image->get_data();
	float *mask_data = mask->get_data();
	size_t size = image->get_xsize() * image->get_ysize() * image->get_zsize();
	double sum = 0;
	size_t n_norm = 0;

	for (size_t i = 0; i < size; i++) {
		if (mask_data[i] > 0.5f) {
			sum += data[i];
			n_norm++;
		}
	}

	float mean = 0;
	if (n_norm == 0) {
		mean = image->get_edge_mean();
	}
	else {
		mean = (float) sum / n_norm;
	}

	return mean;
}

float NormalizeEdgeMeanProcessor::calc_mean(EMData * image) const
{
	if (!image) {
		LOGWARN("NULL Image");
		return 0;
	}
	return image->get_edge_mean();
}

float NormalizeCircleMeanProcessor::calc_mean(EMData * image) const
{
	if (!image) {
		LOGWARN("NULL Image");
		return 0;
	}
	return image->get_circle_mean();
}


float NormalizeMaxMinProcessor::calc_sigma(EMData * image) const
{
	if (!image) {
		LOGWARN("NULL Image");
		return 0;
	}
	float maxval = image->get_attr("maximum");
	float minval = image->get_attr("minimum");
	return (maxval + minval) / 2;
}

float NormalizeMaxMinProcessor::calc_mean(EMData * image) const
{
	if (!image) {
		LOGWARN("NULL Image");
		return 0;
	}
	float maxval = image->get_attr("maximum");
	float minval = image->get_attr("minimum");
	return (maxval - minval) / 2;
}

float NormalizeLREdgeMeanProcessor::calc_mean(EMData * image) const
{
	if (!image) {
		LOGWARN("NULL Image");
		return 0;
	}
	double sum = 0;
	int nx = image->get_xsize();
	int ny = image->get_ysize();
	int nz = image->get_zsize();
	float *d = image->get_data();
	int nyz = ny * nz;

	for (int i = 0; i < nyz; i++) {
		int l = i * nx;
		int r = l + nx - 2;
		sum += d[l] + d[l + 1] + d[r] + d[r + 1];
	}
	float mean = (float) sum / (4 * nyz);
	return mean;
}

void NormalizeRowProcessor::process_inplace(EMData * image)
{
	if (!image) {
		LOGWARN("NULL Image");
		return;
	}

	if (image->get_zsize() > 1) {
		LOGERR("row normalize only works for 2D image");
		return;
	}

	float *rdata = image->get_data();
	int nx = image->get_xsize();
	int ny = image->get_ysize();

	for (int y = 0; y < ny; y++) {
		double row_sum = 0;
		for (int x = 0; x < nx; x++) {
			row_sum += rdata[x + y * nx];
		}

		double row_mean = row_sum / nx;
		if (row_mean <= 0) {
			row_mean = 1;
		}

		for (int x = 0; x < nx; x++) {
			rdata[x + y * nx] /= (float)row_mean;
		}
	}

	image->update();
}

float NormalizeStdProcessor::calc_mean(EMData * image) const
{
	if (!image) {
		LOGWARN("NULL Image");
		return 0;
	}
	return image->get_attr("mean");
}

void NormalizeToStdProcessor::process_inplace(EMData * image)
{
	if (!image) {
		LOGWARN("NULL Image");
		return;
	}
	EMData *noisy = 0;
	if (params.has_key("noisy")) {
		noisy = params["noisy"];
	}
	else {
		noisy = new EMData();
		noisy->read_image((const char *) params["noisyfile"]);
	}

	int nx = image->get_xsize();
	int ny = image->get_ysize();
	int nz = image->get_zsize();

	int size = nx * ny * nz;

	if (!EMUtil::is_same_size(image, noisy)) {
		LOGERR("normalize to on different sized images");
		return;
	}

	int keepzero = params["keepzero"];
	int invert = params["invert"];
	float mult_factor = params["mult"];
	float add_factor = params["add"];
	float sigmax = params["sigmax"];

	float *this_data = image->get_data();
	float *noisy_data = noisy->get_data();
	float m = 0;
	float b = 0;

	if (!invert) {
		Util::calc_least_square_fit(size, this_data, noisy_data, &m, &b, keepzero,sigmax*(float)noisy->get_attr("sigma"));
	}

	if (invert || m < 0) {
		Util::calc_least_square_fit(size, noisy_data, this_data, &m, &b, keepzero,sigmax*(float)image->get_attr("sigma"));

		if (m < 0) {
			b = 0;
			m = 1;
		}
		else if (m > 0) {
			b = -b / m;
			m = 1.0f / m;
		}
	}

	(*image) *= m;

	if (keepzero) {
		float *rdata = image->get_data();
		for (int i = 0; i < size; i++) {
			if (rdata[i]) {
				rdata[i] += b;
			}
		}
		image->update();
	}
	else {
		(*image) += b;
	}

	mult_factor *= m;
	add_factor *= b;

	params["mult"] = mult_factor;
	params["add"] = add_factor;

	if (!params.has_key("noisy")) {
		if( noisy )
		{
			delete noisy;
			noisy = 0;
		}
	}
}


void NormalizeToLeastSquareProcessor::process_inplace(EMData * image)
{
	if (!image) {
		LOGWARN("NULL Image");
		return;
	}

	EMData *to = params["to"];

	float low_threshold = FLT_MIN;
	string low_thr_name = "low_threshold";
	if (params.has_key(low_thr_name)) {
		low_threshold = params[low_thr_name];
	}

	float high_threshold = FLT_MAX;
	string high_thr_name = "high_threshold";
	if (params.has_key(high_thr_name)) {
		high_threshold = params[high_thr_name];
	}

	float *rawp = image->get_data();
	float *refp = to->get_data();

	int nx = image->get_xsize();
	int ny = image->get_ysize();
	int nz = image->get_zsize();
	int size = nx * ny * nz;

	float sum_x = 0;
	float sum_y = 0;
	int count = 0;

	for (int i = 0; i < size; i++) {
		if (refp[i] >= low_threshold && refp[i] <= high_threshold && refp[i] != 0) {
			count++;
			sum_x += refp[i];
			sum_y += rawp[i];
		}
	}

	float sum_x_mean = sum_x / count;
	float sum_tt = 0;
	float b = 0;

	for (int i = 0; i < size; i++) {
		if (refp[i] >= low_threshold && refp[i] <= high_threshold && refp[i] != 0) {
			float t = refp[i] - sum_x_mean;
			sum_tt += t * t;
			b += t * rawp[i];
		}
	}

	b /= sum_tt;

	float a = (sum_y - sum_x * b) / count;
	float scale = 1 / b;
	float shift = -a / b;

	for (int i = 0; i < size; i++) {
		rawp[i] = (rawp[i] - a) / b;
	}

	image->update();

	params["scale"] = scale;
	params["shift"] = shift;
}



void BilateralProcessor::process_inplace(EMData * image)
{
	if (!image) {
		LOGWARN("NULL Image");
		return;
	}
	
	float distance_sigma = params["distance_sigma"];
	float value_sigma = params["value_sigma"];
	int max_iter = params["niter"];
	int half_width = params["half_width"];

	if (half_width < distance_sigma) {
		LOGWARN("localwidth(=%d) should be larger than distance_sigma=(%f)\n",
							half_width, distance_sigma);
	}

	distance_sigma *= distance_sigma;

	float image_sigma = image->get_attr("sigma");
	if (image_sigma > value_sigma) {
		LOGWARN("image sigma(=%f) should be smaller than value_sigma=(%f)\n",
							image_sigma, value_sigma);
	}
	value_sigma *= value_sigma;

	int nx = image->get_xsize();
	int ny = image->get_ysize();
	int nz = image->get_zsize();

	if(nz==1) {	//for 2D image
		int width=nx, height=ny; 

		int i,j,m,n;
	
		float tempfloat1,tempfloat2,tempfloat3;
		int   index1,index2,index;
		int   Iter;
		int   tempint1,tempint3;
	
		tempint1=width;
		tempint3=width+2*half_width;
	
		float* mask=(float*)calloc((2*half_width+1)*(2*half_width+1),sizeof(float));
		float* OrgImg=(float*)calloc((2*half_width+width)*(2*half_width+height),sizeof(float));
		float* NewImg=image->get_data();
	
		for(m=-(half_width);m<=half_width;m++)
			for(n=-(half_width);n<=half_width;n++) {
	       	   index=(m+half_width)*(2*half_width+1)+(n+half_width);
	       	   mask[index]=exp((float)(-(m*m+n*n)/distance_sigma/2.0));
	  	}
	
		//printf("entering bilateral filtering process \n");
	
		Iter=0;
		while(Iter<max_iter) {
			for(i=0;i<height;i++)
	    		for(j=0;j<width;j++) {
		    		index1=(i+half_width)*tempint3+(j+half_width);
					index2=i*tempint1+j;
	        		OrgImg[index1]=NewImg[index2];      
	   		}
	
			// Mirror Padding
			for(i=0;i<height;i++){	
				for(j=0;j<half_width;j++) OrgImg[(i+half_width)*tempint3+(j)]=OrgImg[(i+half_width)*tempint3+(2*half_width-j)];
				for(j=0;j<half_width;j++) OrgImg[(i+half_width)*tempint3+(j+width+half_width)]=OrgImg[(i+half_width)*tempint3+(width+half_width-j-2)];
	   		}
			for(i=0;i<half_width;i++){	
				for(j=0;j<(width+2*half_width);j++) OrgImg[i*tempint3+j]=OrgImg[(2*half_width-i)*tempint3+j];
				for(j=0;j<(width+2*half_width);j++) OrgImg[(i+height+half_width)*tempint3+j]=OrgImg[(height+half_width-2-i)*tempint3+j];
			}
	
			//printf("finish mirror padding process \n");
			//now mirror padding have been done
	
			for(i=0;i<height;i++){
				//printf("now processing the %d th row \n",i);
				for(j=0;j<width;j++){
					tempfloat1=0.0; tempfloat2=0.0;
					for(m=-(half_width);m<=half_width;m++)
						for(n=-(half_width);n<=half_width;n++){
							index =(m+half_width)*(2*half_width+1)+(n+half_width);
							index1=(i+half_width)*tempint3+(j+half_width);
							index2=(i+half_width+m)*tempint3+(j+half_width+n);
							tempfloat3=(OrgImg[index1]-OrgImg[index2])*(OrgImg[index1]-OrgImg[index2]);
	
							tempfloat3=mask[index]*(1.0f/(1+tempfloat3/value_sigma));	// Lorentz kernel
							//tempfloat3=mask[index]*exp(tempfloat3/Sigma2/(-2.0));	// Guassian kernel
							tempfloat1+=tempfloat3;
	
							tempfloat2+=tempfloat3*OrgImg[(i+half_width+m)*tempint3+(j+half_width+n)];     
	        			}
					NewImg[i*width+j]=tempfloat2/tempfloat1;
				}
			}
	   		Iter++;
	    }
	
	    //printf("have finished %d  th iteration\n ",Iter);
//		doneData();
		free(mask); 
		free(OrgImg);
		// end of BilaFilter routine
		
	}
	else {	//3D case
		int width = nx;
		int height = ny;
		int slicenum = nz;
	
		int slice_size = width * height;
		int new_width = width + 2 * half_width;
		int new_slice_size = (width + 2 * half_width) * (height + 2 * half_width);
	
		int width1 = 2 * half_width + 1;
		int mask_size = width1 * width1;
		int old_img_size = (2 * half_width + width) * (2 * half_width + height);
	
		int zstart = -half_width;
		int zend = -half_width;
		int is_3d = 0;
		if (nz > 1) {
			mask_size *= width1;
			old_img_size *= (2 * half_width + slicenum);
			zend = half_width;
			is_3d = 1;
		}
	
		float *mask = (float *) calloc(mask_size, sizeof(float));
		float *old_img = (float *) calloc(old_img_size, sizeof(float));
	
		float *new_img = image->get_data();
	
		for (int p = zstart; p <= zend; p++) {
			int cur_p = (p + half_width) * (2 * half_width + 1) * (2 * half_width + 1);
	
			for (int m = -half_width; m <= half_width; m++) {
				int cur_m = (m + half_width) * (2 * half_width + 1) + half_width;
	
				for (int n = -half_width; n <= half_width; n++) {
					int l = cur_p + cur_m + n;
					mask[l] = exp((float) (-(m * m + n * n + p * p * is_3d) / distance_sigma / 2.0f));
				}
			}
		}
	
		int iter = 0;
		while (iter < max_iter) {
			for (int k = 0; k < slicenum; k++) {
				int cur_k1 = (k + half_width) * new_slice_size * is_3d;
				int cur_k2 = k * slice_size;
	
				for (int i = 0; i < height; i++) {
					int cur_i1 = (i + half_width) * new_width;
					int cur_i2 = i * width;
	
					for (int j = 0; j < width; j++) {
						int k1 = cur_k1 + cur_i1 + (j + half_width);
						int k2 = cur_k2 + cur_i2 + j;
						old_img[k1] = new_img[k2];
					}
				}
			}
	
			for (int k = 0; k < slicenum; k++) {
				int cur_k = (k + half_width) * new_slice_size * is_3d;
	
				for (int i = 0; i < height; i++) {
					int cur_i = (i + half_width) * new_width;
	
					for (int j = 0; j < half_width; j++) {
						int k1 = cur_k + cur_i + j;
						int k2 = cur_k + cur_i + (2 * half_width - j);
						old_img[k1] = old_img[k2];
					}
	
					for (int j = 0; j < half_width; j++) {
						int k1 = cur_k + cur_i + (width + half_width + j);
						int k2 = cur_k + cur_i + (width + half_width - j - 2);
						old_img[k1] = old_img[k2];
					}
				}
	
	
				for (int i = 0; i < half_width; i++) {
					int i2 = i * new_width;
					int i3 = (2 * half_width - i) * new_width;
					for (int j = 0; j < (width + 2 * half_width); j++) {
						int k1 = cur_k + i2 + j;
						int k2 = cur_k + i3 + j;
						old_img[k1] = old_img[k2];
					}
	
					i2 = (height + half_width + i) * new_width;
					i3 = (height + half_width - 2 - i) * new_width;
					for (int j = 0; j < (width + 2 * half_width); j++) {
						int k1 = cur_k + i2 + j;
						int k2 = cur_k + i3 + j;
						old_img[k1] = old_img[k2];
					}
				}
			}
	
			for (int k = 0; k < slicenum; k++) {
				int cur_k = (k + half_width) * new_slice_size;
	
				for (int i = 0; i < height; i++) {
					int cur_i = (i + half_width) * new_width;
	
					for (int j = 0; j < width; j++) {
						float f1 = 0;
						float f2 = 0;
						int k1 = cur_k + cur_i + (j + half_width);
	
						for (int p = zstart; p <= zend; p++) {
							int cur_p1 = (p + half_width) * (2 * half_width + 1) * (2 * half_width + 1);
							int cur_p2 = (k + half_width + p) * new_slice_size;
	
							for (int m = -half_width; m <= half_width; m++) {
								int cur_m1 = (m + half_width) * (2 * half_width + 1);
								int cur_m2 = cur_p2 + cur_i + m * new_width + j + half_width;
	
								for (int n = -half_width; n <= half_width; n++) {
									int k = cur_p1 + cur_m1 + (n + half_width);
									int k2 = cur_m2 + n;
									float f3 = Util::square(old_img[k1] - old_img[k2]);
	
									f3 = mask[k] * (1.0f / (1 + f3 / value_sigma));
									f1 += f3;
									int l1 = cur_m2 + n;
									f2 += f3 * old_img[l1];
								}
	
								new_img[k * height * width + i * width + j] = f2 / f1;
							}
						}
					}
				}
			}
			iter++;
		}
		if( mask ) {
			free(mask);
			mask = 0;
		}
	
		if( old_img ) {
			free(old_img);
			old_img = 0;
		}
	}
	
	image->update();
}

void RadialAverageProcessor::process_inplace(EMData * image)
{
	if (!image || image->is_complex()) {
		LOGWARN("only works on real image. do nothing.");
		return;
	}

	float *rdata = image->get_data();
	int nx = image->get_xsize();
	int ny = image->get_ysize();

	vector < float >dist = image->calc_radial_dist(nx / 2, 0, 1,0);

	int c = 0;
	for (int y = 0; y < ny; y++) {
		for (int x = 0; x < nx; x++, c++) {
			float r = (float) hypot(x - nx / 2.0f, y - ny / 2.0f);

			int i = (int) floor(r);
			r -= i;
			if (i >= 0 && i < nx / 2 - 1) {
				rdata[c] = dist[i] * (1.0f - r) + dist[i + 1] * r;
			}
			else if (i < 0) {
				rdata[c] = dist[0];
			}
			else {
				rdata[c] = 0;
			}
		}
	}

	image->update();
}



void RadialSubstractProcessor::process_inplace(EMData * image)
{
	if (!image || image->is_complex()) {
		LOGWARN("only works on real image. do nothing.");
		return;
	}

	float *rdata = image->get_data();
	int nx = image->get_xsize();
	int ny = image->get_ysize();

	vector < float >dist = image->calc_radial_dist(nx / 2, 0, 1,0);

	int c = 0;
	for (int y = 0; y < ny; y++) {
		for (int x = 0; x < nx; x++, c++) {
			float r = (float) hypot(x - nx / 2, y - ny / 2);
			int i = (int) floor(r);
			r -= i;
			if (i >= 0 && i < nx / 2 - 1) {
				rdata[c] -= dist[i] * (1.0f - r) + dist[i + 1] * r;
			}
			else {
				rdata[c] = 0;
			}
		}
	}

	image->update();
}



void FlipProcessor::process_inplace(EMData * image)
{
	if (!image) {
		LOGWARN("NULL Image");
		return;
	}
	
	string axis = (const char*)params["axis"];

	float *d = image->get_data();
	int nx = image->get_xsize();
	int ny = image->get_ysize();
	int nz = image->get_zsize();

	int nxy = nx * ny;

	if (axis == "x" || axis == "X") {		// horizontal flip
		for(int z = 0; z < nz; ++z) {
			for(int y = 0; y < ny; ++y) {
				for(int x = 0; x < nx / 2; ++x) {
					std::swap(d[z*nxy + y*nx + x], d[z*nxy + y*nx + (nx-x-1)]);
				}
			}
		}
	}
	else if (axis == "y" || axis == "Y") {		// vertical flip
		for(int z=0; z<nz; ++z) {
			for(int y=0; y<ny/2; ++y) {
				for(int x=0; x<nx; ++x) {
					std::swap(d[z*nxy + y*nx +x], d[z*nxy + (ny -y -1)*nx +x]);
				}
			}
		} 
	}
	else if (axis == "z" || axis == "Z") {		//z axis flip
		for(int z=0; z<nz/2; ++z) {
			for(int y=0; y<ny; ++y) {
				for(int x=0; x<nx; ++x) {
					std::swap(d[z*nxy + y*nx + x], d[(nz-z-1)*nxy + y*nx + x]);
				}
			}
		}
	}

	image->update();
}
 
void AddNoiseProcessor::process_inplace(EMData * image)
{
	if (!image) {
		LOGWARN("NULL Image");
		return;
	}

	if(params.has_key("seed")) {
#ifdef _WIN32
		srand((int)params["seed"]);
#else 
		srandom((int)params["seed"]);
#endif
	}
	else {
#ifdef _WIN32
		srand(time(0));
#else
		srandom(time(0));
#endif
	}

	float addnoise = params["noise"];
	addnoise *= get_sigma(image);
	float *dat = image->get_data();
	size_t size = image->get_xsize() * image->get_ysize() * image->get_zsize();

	for (size_t j = 0; j < size; j++) {
		dat[j] += Util::get_gauss_rand(addnoise, addnoise / 2);
	}

	image->update();
}

float AddSigmaNoiseProcessor::get_sigma(EMData * image)
{
	if (!image) {
		LOGWARN("NULL Image");
		return 0;
	}
	return image->get_attr("sigma");
}

void FourierToCornerProcessor::process_inplace(EMData * image)
{
	if ( !image->is_complex() ) throw ImageFormatException("Can not Fourier origin shift an image that is not complex");
	
	int nx=image->get_xsize();
	int ny=image->get_ysize();
	int nz=image->get_zsize();
	
	int nxy = nx*ny;
	
	if ( ny == 1 && nz == 1 ){
		cout << "Warning- attempted	Fourier origin shift a 1D image - no action taken" << endl;
		return;
	}
	int yodd = (ny%2==1);
	int zodd = (nz%2==1);
	
	float* rdata = image->get_data();
	
	float tmp[2];
	float* p1;
	float* p2;
		
	if (yodd){
		// Swap the middle slice (with respect to the y direction) with the bottom slice
		// shifting all slices above the middles slice upwards by one pixel, stopping
		// at the middle slice, not if nz = 1 we are not talking about slices, we are
		// talking about rows
		float prev[2];
		for( int s = 0; s < nz; s++ ) {
			for( int c =0; c < nx; c += 2 ) {
				prev[0] = rdata[s*nxy+ny/2*nx+c];
				prev[1] = rdata[s*nxy+ny/2*nx+c+1];	
				for( int r = 0; r <= ny/2; ++r ) {
					float* p1 = &rdata[s*nxy+r*nx+c];
					tmp[0] = p1[0];
					tmp[1] = p1[1];
						
					p1[0] = prev[0];
					p1[1] = prev[1];
						
					prev[0] = tmp[0];
					prev[1] = tmp[1];
				}
			}
		}
	}
		
	// Shift slices (3D) or rows (2D) correctly in the y direction
	for( int s = 0; s < nz; ++s ) {
		for( int r = 0 + yodd; r < ny/2+yodd; ++r ) {
			for( int c =0; c < nx; c += 2 ) {
				p1 = &rdata[s*nxy+r*nx+c];
				p2 = &rdata[s*nxy+(r+ny/2)*nx+c];
					
				tmp[0] = p1[0];
				tmp[1] = p1[1];
					
				p1[0] = p2[0];
				p1[1] = p2[1];
					
				p2[0] = tmp[0];
				p2[1] = tmp[1];
			}
		}
	}
	
	if ( nz != 1 )
	{
		
		if (zodd){
			// Swap the middle slice (with respect to the z direction) and the front slice
			// shifting all behind the front slice towards the middle a distance of 1 voxel,
			// stopping at the middle slice.
			float prev[2];
			for( int r = 0; r < ny; ++r ) {
				for( int c =0; c < nx; c += 2 ) {
					prev[0] = rdata[nz/2*nxy+r*nx+c];
					prev[1] = rdata[nz/2*nxy+r*nx+c+1];	
					for( int s = 0; s <= nz/2; ++s ) {
						float* p1 = &rdata[s*nxy+r*nx+c];
						tmp[0] = p1[0];
						tmp[1] = p1[1];
					
						p1[0] = prev[0];
						p1[1] = prev[1];
					
						prev[0] = tmp[0];
						prev[1] = tmp[1];
					}
				}
			}
		}
		
		// Shift slices correctly in the z direction
		for( int s = 0+zodd; s < nz/2 + zodd; ++s ) {
			for( int r = 0; r < ny; ++r ) {
				for( int c =0; c < nx; c += 2 ) {
					p1 = &rdata[s*nxy+r*nx+c];
					p2 = &rdata[(s+nz/2)*nxy+r*nx+c];
					
					tmp[0] = p1[0];
					tmp[1] = p1[1];
					
					p1[0] = p2[0];
					p1[1] = p2[1];
					
					p2[0] = tmp[0];
					p2[1] = tmp[1];
				}
			}
		}
		
	}
}

void FourierToCenterProcessor::process_inplace(EMData * image)
{
	if ( !image->is_complex() ) throw ImageFormatException("Can not Fourier origin shift an image that is not complex");
	
	int nx=image->get_xsize();
	int ny=image->get_ysize();
	int nz=image->get_zsize();
	
	int nxy = nx*ny;
	
	if ( ny == 1 && nz == 1 ){
		cout << "Warning- attempted	Fourier origin shift a 1D image - no action taken" << endl;
		return;
	}
	
	int yodd = (ny%2==1);
	int zodd = (nz%2==1);
	
	float* rdata = image->get_data();
	
	float tmp[2];
	float* p1;
	float* p2;
		
	if (yodd){
		// In 3D this is swapping the bottom slice (with respect to the y direction) and the middle slice,
		// shifting all slices below the middle slice down one. In 2D it is equivalent, but in terms of rows.
		float prev[2];
		for( int s = 0; s < nz; s++ ) {
			for( int c =0; c < nx; c += 2 ) {
				prev[0] = rdata[s*nxy+c];
				prev[1] = rdata[s*nxy+c+1];	
				for( int r = ny/2; r >= 0; --r ) {
					float* p1 = &rdata[s*nxy+r*nx+c];
					tmp[0] = p1[0];
					tmp[1] = p1[1];
						
					p1[0] = prev[0];
					p1[1] = prev[1];
						
					prev[0] = tmp[0];
					prev[1] = tmp[1];
				}
			}
		}
	}
			
	// 3D - Shift slices correctly in the y direction, 2D - shift rows
	for( int s = 0; s < nz; ++s ) {
		for( int r = 0; r < ny/2; ++r ) {
			for( int c =0; c < nx; c += 2 ) {
				p1 = &rdata[s*nxy+r*nx+c];
				p2 = &rdata[s*nxy+(r+ny/2+yodd)*nx+c];
					
				tmp[0] = p1[0];
				tmp[1] = p1[1];
					
				p1[0] = p2[0];
				p1[1] = p2[1];
					
				p2[0] = tmp[0];
				p2[1] = tmp[1];
			}
		}
	}
	
	if ( nz != 1 )  {
		if (zodd){
			// Swap the front slice (with respect to the z direction) and the middle slice
			// shifting all slices behind the middles slice towards the front slice 1 voxel.
			float prev[2];
			for( int r = 0; r < ny; ++r ) {
				for( int c =0; c < nx; c += 2 ) {
					prev[0] = rdata[r*nx+c];
					prev[1] = rdata[r*nx+c+1];	
					for( int s = nz/2; s >= 0; --s ) {
						float* p1 = &rdata[s*nxy+r*nx+c];
						tmp[0] = p1[0];
						tmp[1] = p1[1];
					
						p1[0] = prev[0];
						p1[1] = prev[1];
					
						prev[0] = tmp[0];
						prev[1] = tmp[1];
					}
				}
			}
		}
		
		// Shift slices correctly in the y direction
		for( int s = 0; s < nz/2; ++s ) {
			for( int r = 0; r < ny; ++r ) {
				for( int c =0; c < nx; c += 2 ) {
					p1 = &rdata[s*nxy+r*nx+c];
					p2 = &rdata[(s+nz/2+zodd)*nxy+r*nx+c];
					
					tmp[0] = p1[0];
					tmp[1] = p1[1];
					
					p1[0] = p2[0];
					p1[1] = p2[1];
					
					p2[0] = tmp[0];
					p2[1] = tmp[1];
				}
			}
		}
	}
}

void Phase180Processor::fourier_phaseshift180(EMData * image)
{
	if ( !image->is_complex() ) throw ImageFormatException("Can not handle images that are not complex in fourier phase shift 180");
	
	int nx = image->get_xsize();
	int ny = image->get_ysize();
	int nz = image->get_zsize();

	int nxy = nx * ny;
	
	float *rdata = image->get_data();
	
	int of=0;
	if (((ny/2)%2)+((nz/2)%2)==1) of=1;
	
	for (int k = 0; k < nz; k++) {
		int k2 = k * nxy;
		
		for (int j = 0; j < ny; j++) {
			int i = ((k+j)%2==of?2:0);
			int j2 = j * nx + k2;
			
			for (; i < nx; i += 4) {
				rdata[i + j2] *= -1.0f;
				rdata[i + j2 + 1] *= -1.0f;
			}
		}
	}
}

void Phase180Processor::swap_corners_180(EMData * image)
{
	int nx = image->get_xsize();
	int ny = image->get_ysize();
	int nz = image->get_zsize();
	
	int xodd = (nx % 2) == 1;
	int yodd = (ny % 2) == 1;
	int zodd = (nz % 2) == 1;
	
	int nxy = nx * ny;
	
	float *rdata = image->get_data();
	
	if ( ny == 1 && nz == 1 ){
		throw ImageDimensionException("Error, cannot handle 1D images. This function should not have been called");
	}
	else if ( nz == 1 ) {
		
		// Swap the bottom left and top right corners
		for ( int r = 0; r < ny/2; ++r ) {
			for ( int c = 0; c < nx/2; ++c) {
				int idx1 = r*nx + c;
				int idx2 = (r+ny/2+yodd)*nx + c + nx/2+xodd;
				float tmp = rdata[idx1];
				rdata[idx1] = rdata[idx2];
				rdata[idx2] = tmp;
			} 
		}
		
		// Swap the top left and bottom right corners
		for ( int r = ny-1; r >= (ny/2+yodd); --r ) {
			for ( int c = 0; c < nx/2; ++c) {
				int idx1 = r*nx + c;
				int idx2 = (r-ny/2-yodd)*nx + c + nx/2+xodd;
				float tmp = rdata[idx1];
				rdata[idx1] = rdata[idx2];
				rdata[idx2] = tmp;
			}
		}
	}
	else // nx && ny && nz are greater than 1
	{
		float tmp;
		// Swap the bottom left front and back right top quadrants
		for ( int s = 0; s < nz/2; ++s ) {
			for ( int r = 0; r < ny/2; ++r ) {
				for ( int c = 0; c < nx/2; ++ c) {
					int idx1 = s*nxy+r*nx+c;
					int idx2 = (s+nz/2+zodd)*nxy+(r+ny/2+yodd)*nx+c+nx/2+xodd;
					tmp = rdata[idx1];
					rdata[idx1] = rdata[idx2];
					rdata[idx2] = tmp;
				}
			}
		}
		// Swap the bottom right front and back left top quadrants
		for ( int s = 0; s < nz/2; ++s ) {
			for ( int r = 0; r < ny/2; ++r ) {
				for ( int c = nx-1; c >= (nx/2+xodd); --c) {
					int idx1 = s*nxy+r*nx+c;
					int idx2 = (s+nz/2+zodd)*nxy+(r+ny/2+yodd)*nx+c-nx/2-xodd;
					tmp = rdata[idx1];
					rdata[idx1] = rdata[idx2];
					rdata[idx2] = tmp;
				}
			}
		}
		// Swap the top right front and back left bottom quadrants
		for ( int s = 0; s < nz/2; ++s ) {
			for ( int r = ny-1; r >= (ny/2+yodd); --r ) {
				for ( int c = nx-1; c >= (nx/2+xodd); --c) {
					int idx1 = s*nxy+r*nx+c;
					int idx2 = (s+nz/2+zodd)*nxy+(r-ny/2-yodd)*nx+c-nx/2-xodd;
					tmp = rdata[idx1];
					rdata[idx1] = rdata[idx2];
					rdata[idx2] = tmp;
				}
			}
		}
		// Swap the top left front and back right bottom quadrants
		for ( int s = 0; s < nz/2; ++s ) {
			for ( int r = ny-1; r >= (ny/2+yodd); --r ) {
				for ( int c = 0; c < nx/2; ++ c) {
					int idx1 = s*nxy+r*nx+c;
					int idx2 = (s+nz/2+zodd)*nxy+(r-ny/2-yodd)*nx+c+nx/2+xodd;
					tmp = rdata[idx1];
					rdata[idx1] = rdata[idx2];
					rdata[idx2] = tmp;
				}
			}
		}
	}
}

void Phase180Processor::swap_central_slices_180(EMData * image)
{
	int nx = image->get_xsize();
	int ny = image->get_ysize();
	int nz = image->get_zsize();
	
	int xodd = (nx % 2) == 1;
	int yodd = (ny % 2) == 1;
	int zodd = (nz % 2) == 1;
	
	int nxy = nx * ny;
	
	float *rdata = image->get_data();
	
	if ( ny == 1 && nz == 1 ){
		throw ImageDimensionException("Error, cannot handle 1D images. This function should not have been called");
	}
	else if ( nz == 1 ) {
		float tmp;
		if ( yodd ) {
			// Iterate along middle row, swapping values where appropriate
			int r = ny/2;
			for ( int c = 0; c < nx/2; ++c ) {
				int idx1 = r*nx + c;
				int idx2 = r*nx + c + nx/2+ xodd;
				tmp = rdata[idx1];
				rdata[idx1] = rdata[idx2];
				rdata[idx2] = tmp;
			}
		}
		
		if ( xodd )	{
			// Iterate along the central column, swapping values where appropriate
			int c = nx/2;
			for (  int r = 0; r < ny/2; ++r ) {
				int idx1 = r*nx + c;
				int idx2 = (r+ny/2+yodd)*nx + c;
				tmp = rdata[idx1];
				rdata[idx1] = rdata[idx2];
				rdata[idx2] = tmp;
			}
		}
	}
	else // nx && ny && nz are greater than 1
	{
		float tmp;
		if ( xodd ) {
			// Iterate along the x = nx/2 slice, swapping values where appropriate
			int c = nx/2;
			for( int s = 0; s < nz/2; ++s ) {
				for ( int r = 0; r < ny/2; ++r ) {
					int idx1 = s*nxy+r*nx+c;
					int idx2 = (s+nz/2+zodd)*nxy+(r+ny/2+yodd)*nx+c;
					tmp = rdata[idx1];
					rdata[idx1] = rdata[idx2];
					rdata[idx2] = tmp;
				}
			}
			
			for( int s = nz-1; s >= (nz/2+zodd); --s ) {
				for ( int r = 0; r < ny/2; ++r ) {
					int idx1 = s*nxy+r*nx+c;
					int idx2 = (s-nz/2-zodd)*nxy+(r+ny/2+yodd)*nx+c;
					tmp = rdata[idx1];
					rdata[idx1] = rdata[idx2];
					rdata[idx2] = tmp;
				}
			}
		}
		if ( yodd ) {
			// Iterate along the y = ny/2 slice, swapping values where appropriate
			int r = ny/2;
			for( int s = 0; s < nz/2; ++s ) {
				for ( int c = 0; c < nx/2; ++c ) {
					int idx1 = s*nxy+r*nx+c;
					int idx2 =(s+nz/2+zodd)*nxy+r*nx+c+nx/2+xodd;
					tmp = rdata[idx1];
					rdata[idx1] = rdata[idx2];
					rdata[idx2] = tmp;
				}
			}
			
			for( int s = nz-1; s >= (nz/2+zodd); --s ) {
				for ( int c = 0; c < nx/2; ++c ) {
					int idx1 = s*nxy+r*nx+c;
					int idx2 = (s-nz/2-zodd)*nxy+r*nx+c+nx/2+xodd;
					tmp = rdata[idx1];
					rdata[idx1] = rdata[idx2];
					rdata[idx2] = tmp;
				}
			}
		}
		if ( zodd ) {
			// Iterate along the z = nz/2 slice, swapping values where appropriate
			int s = nz/2;
			for( int r = 0; r < ny/2; ++r ) {
				for ( int c = 0; c < nx/2; ++c ) {
					int idx1 = s*nxy+r*nx+c;
					int idx2 = s*nxy+(r+ny/2+yodd)*nx+c+nx/2+xodd;
					tmp = rdata[idx1];
					rdata[idx1] = rdata[idx2];
					rdata[idx2] = tmp;
				}
			}
			
			for( int r = ny-1; r >= (ny/2+yodd); --r ) {
				for ( int c = 0; c < nx/2; ++c ) {
					int idx1 = s*nxy+r*nx+c;
					int idx2 = s*nxy+(r-ny/2-yodd)*nx+c+nx/2+xodd;
					tmp = rdata[idx1];
					rdata[idx1] = rdata[idx2];
					rdata[idx2] = tmp;
				}
			}
		}
	}
}

void PhaseToCornerProcessor::process_inplace(EMData * image)
{
	if (!image)	throw NullPointerException("Error: attempt to phase shift a null image");

	if (image->is_complex()) {
		fourier_phaseshift180(image);
		return;
	}

	int nx = image->get_xsize();
	int ny = image->get_ysize();
	int nz = image->get_zsize();

	if ( ny == 1 && nz == 1 && nx == 1) return;
	
	int nxy = nx * ny;
	
	float *rdata = image->get_data();

	bool xodd = (nx % 2) == 1;
	bool yodd = (ny % 2) == 1;
	bool zodd = (nz % 2) == 1;
	
	if ( ny == 1 && nz == 1 ){
		if (xodd){
			// Put the last pixel in the center, shifting the contents 
			// to right of the center one step to the right
			float in_x = rdata[nx-1];
			float tmp;
			for ( int i = nx/2; i < nx; ++i ) {
				tmp = rdata[i];
				rdata[i] = in_x;
				in_x = tmp;
			}
		}
		// now the operation is straight forward
		for ( int i = 0; i < nx/2; ++i ) {
			int idx = i+nx/2+xodd;
			float tmp = rdata[i];
			rdata[i] = rdata[idx];
			rdata[idx] = tmp;
		}
		
	}
	else if ( nz == 1 ) {
		if (yodd) {
			// Tranfer the top row into the middle row,
			// shifting all pixels above and including the current middle up one.
			for ( int c = 0; c < nx; ++c ) {
				// Get the value in the top row
				float last_val = rdata[(ny-1)*nx + c];
				float tmp;
				for ( int r = ny/2; r < ny; ++r ){
					int idx =r*nx+c;
					tmp = rdata[idx];
					rdata[idx] = last_val;
					last_val = tmp;
				}
			}
		}
		
		if (xodd) {
			// Transfer the right most column into the center column
			// Shift all columns right of and including center to the right one pixel
			for ( int r  = 0; r < ny; ++r ) {
				float last_val = rdata[(r+1)*nx -1];
				float tmp;
				for ( int c = nx/2; c < nx; ++c ){
					int idx =r*nx+c;
					tmp = rdata[idx];
					rdata[idx] = last_val;
					last_val = tmp;
				}
			}
		}
		// It is important central slice shifting come after the previous two operations
		swap_central_slices_180(image);
		// Now the corners of the image can be shifted...
		swap_corners_180(image);
		
	}
	else
	{
		float tmp;
		if (zodd) {
			// Tranfer the back slice into the middle slice,
			// shifting all pixels beyond and including the middle slice back one.
			for (int r = 0; r < ny; ++r){
				for (int c = 0; c < nx; ++c) {
					float last_val = rdata[(nz-1)*nxy+r*nx+c];
					for (int s = nz/2; s < nz; ++s) {
						int idx = s*nxy+r*nx+c;
						tmp = rdata[idx];
						rdata[idx] = last_val;
						last_val = tmp;
					}
				}
			}
		}
		if (yodd) {
			// Tranfer the top slice into the middle slice,
			// shifting all pixels above and including the middle slice up one.
			for (int s = 0; s < nz; ++s) {
				for (int c = 0; c < nx; ++c) {
				float last_val = rdata[s*nxy+(ny-1)*nx+c];
					for (int r = ny/2; r < ny; ++r){
						int idx = s*nxy+r*nx+c;
						tmp = rdata[idx];
						rdata[idx] = last_val;
						last_val = tmp;
					}
				}
			}
		}
		if (xodd) {
			// Transfer the right most slice into the central slice
			// Shift all pixels to right of and including center slice to the right one pixel
			for (int s = 0; s < nz; ++s) {
				for (int r = 0; r < ny; ++r) {
					float last_val = rdata[s*nxy+r*nx+nx-1];
					for (int c = nx/2; c < nx; ++c){
						int idx = s*nxy+r*nx+c;
						tmp = rdata[idx];
						rdata[idx] = last_val;
						last_val = tmp;
					}
				}
			}
		}
		// Now swap the various parts in the central slices
		swap_central_slices_180(image);
		// Now shift the corners
		swap_corners_180(image);
	}
}


void PhaseToCenterProcessor::process_inplace(EMData * image)
{
	if (!image)	throw NullPointerException("Error: attempt to phase shift a null image");

	if (image->is_complex()) {
		fourier_phaseshift180(image);
		return;
	}
	
	int nx = image->get_xsize();
	int ny = image->get_ysize();
	int nz = image->get_zsize();

	if ( ny == 1 && nz == 1 && nx == 1) return;
	
	int nxy = nx * ny;
	
	float *rdata = image->get_data();

	bool xodd = (nx % 2) == 1;
	bool yodd = (ny % 2) == 1;
	bool zodd = (nz % 2) == 1;
	
	if ( ny == 1 && nz == 1 ){
		if (xodd) {
			// Put the center pixel at the end, shifting the contents 
			// to right of the center one step to the left
			float in_x = rdata[nx/2];
			float tmp;
			for ( int i = nx-1; i >= nx/2; --i ) {
				tmp = rdata[i];
				rdata[i] = in_x;
				in_x = tmp;
			}
		}
		// now the operation is straight forward
		for ( int i = 0; i < nx/2; ++i ) {
			int idx = i + nx/2;
			float tmp = rdata[i];
			rdata[i] = rdata[idx];
			rdata[idx] = tmp;
		}	
	}
	else if ( nz == 1 ){
		// The order in which these operations occur literally undoes what the
		// PhaseToCornerProcessor did to the image.
		// First, the corners sections of the image are swapped appropriately
		swap_corners_180(image);
		// Second, central pixel lines are swapped
		swap_central_slices_180(image);
		
		float tmp;
		// Third, appropriate sections of the image are cyclically shifted by one pixel
		if (xodd) {
			// Transfer the middle column to the far right 
			// Shift all from the far right to (but not including the) middle one to the left
			for ( int r  = 0; r < ny; ++r ) {
				float last_val = rdata[r*nx+nx/2];
				for ( int c = nx-1; c >=  nx/2; --c ){
					int idx = r*nx+c;
					tmp = rdata[idx];
					rdata[idx] = last_val;
					last_val = tmp;
				}
			}
		}	
		if (yodd) {
			// Tranfer the middle row to the top,
			// shifting all pixels from the top row down one, until  but not including the) middle
			for ( int c = 0; c < nx; ++c ) {
				// Get the value in the top row
				float last_val = rdata[ny/2*nx + c];
				for ( int r = ny-1; r >= ny/2; --r ){
					int idx = r*nx+c;
					tmp = rdata[idx];
					rdata[idx] = last_val;
					last_val = tmp;
				}
			}
		}
	}
	else
	{
		// The order in which these operations occur literally undoes the
		// PhaseToCornerProcessor operation - in 3D.
		// First, the corner quadrants of the voxel volume are swapped
		swap_corners_180(image);
		// Second, appropriate parts of the central slices are swapped
		swap_central_slices_180(image);

		float tmp;
		// Third, appropriate sections of the image are cyclically shifted by one voxel
		if (xodd) {
			// Transfer the central slice in the x direction to the far right
			// moving all slices on the far right toward the center one pixel, until
			// the center x slice is ecountered
			for (int s = 0; s < nz; ++s) {
				for (int r = 0; r < ny; ++r) {
					float last_val = rdata[s*nxy+r*nx+nx/2];
					for (int c = nx-1; c >= nx/2; --c){
						int idx = s*nxy+r*nx+c;
						tmp = rdata[idx];
						rdata[idx] = last_val;
						last_val = tmp;
					}
				}
			}
		}
		if (yodd) {
			// Tranfer the central slice in the y direction to the top
			// shifting all pixels below it down on, until the center y slice is encountered.
			for (int s = 0; s < nz; ++s) {
				for (int c = 0; c < nx; ++c) {
					float last_val = rdata[s*nxy+ny/2*nx+c];
					for (int r = ny-1; r >= ny/2; --r){
						int idx = s*nxy+r*nx+c;
						tmp = rdata[idx];
						rdata[idx] = last_val;
						last_val = tmp;
					}
				}
			}
		}
		if (zodd) {
			// Tranfer the central slice in the z direction to the back
			// shifting all pixels beyond and including the middle slice back one.
			for (int r = 0; r < ny; ++r){
				for (int c = 0; c < nx; ++c) {
					float last_val = rdata[nz/2*nxy+r*nx+c];
					for (int s = nz-1; s >= nz/2; --s) {
						int idx = s*nxy+r*nx+c;
						tmp = rdata[idx];
						rdata[idx] = last_val;
						last_val = tmp;
					}
				}
			}
		}
		
		
	}
}

void AutoMask2DProcessor::process_inplace(EMData * image)
{
	if (!image) {
		LOGWARN("NULL Image");
		return;
	}
	
	int nz = image->get_zsize();
	if (nz > 1) {
		LOGERR("%s Processor doesn't support 3D model", get_name().c_str());
		throw ImageDimensionException("3D model not supported");
	}
	
	float threshold = params["threshold"];
	float filter = 0.1f;
	if (params.has_key("filter")) {
		filter = params["filter"];
	}

	EMData *d = image->copy();

	if (threshold > 0) {
		threshold = -threshold;
		d->mult(-1);
	}

	int ny = image->get_ysize();

	d->process_inplace("eman1.filter.lowpass.gaussian", Dict("lowpass", (filter * ny / 2)));
	d->process_inplace("eman1.filter.highpass.gaussian", Dict("highpass", 0));

	d->process_inplace("normalize");

	int d_nx = d->get_xsize();
	int d_ny = d->get_ysize();
	int d_size = d_nx * d_ny;

	float *dn = (float *) calloc(d_size, sizeof(float));
	float *dc = (float *) calloc(d_size, sizeof(float));
	float *dd = d->get_data();

	int k = 0;
	const float dn_const = 10;
	float d_sigma = d->get_attr_dict().get("sigma");
	float threshold_sigma = threshold * d_sigma;

	for (int l = 0; l < d_ny; l++) {
		for (int j = 0; j < d_nx; j++) {
			if (hypot(l - d_ny / 2, j - d_nx / 2) >= d_ny / 2) {
				dn[k] = -dn_const;
			}
			else if (dd[k] < threshold_sigma) {
				dn[k] = dn_const;
			}
			else {
				dn[k] = 0;
			}
			k++;
		}
	}

	int ch = 1;
	while (ch) {
		ch = 0;
		memcpy(dc, dn, d_size * sizeof(float));

		for (int l = 1; l < d_ny - 1; l++) {
			for (int j = 1; j < d_nx - 1; j++) {
				k = j + l * d_nx;
				if (dn[k] == dn_const) {
					if (dn[k - 1] == -dn_const || dn[k + 1] == -dn_const ||
						dn[k - d_nx] == -dn_const || dn[k + d_nx] == -dn_const) {
						dn[k] = -dn_const;
						ch = 1;
					}
				}
			}
		}
	}

	k = 0;
	float threshold_sigma2 = threshold * d_sigma * 2;

	for (int l = 0; l < d_ny; l++) {
		for (int j = 0; j < d_nx; j++) {
			if (l == 0 || l == d_ny - 1 || j == 0 || j == d_nx - 1) {
				dn[k] = -dn_const;
			}
			else if (dd[k] > threshold_sigma2 && dn[k] == dn_const) {
				dn[k] = 0;
			}
			k++;
		}
	}

	ch = 1;
	const float v_const = 0.25f;

	while (ch) {
		ch = 0;
		memcpy(dc, dn, d_size * sizeof(float));

		for (int l = 1; l < d_ny - 1; l++) {
			for (int j = 1; j < d_nx - 1; j++) {
				int k = j + l * d_nx;

				if (dn[k] != -dn_const && dn[k] != dn_const) {
					float v = Util::get_max(dc[k - 1], dc[k + 1], dc[k - d_nx], dc[k + d_nx]);
					float v2 = Util::get_min(dc[k - 1], dc[k + 1], dc[k - d_nx], dc[k + d_nx]);

					if (v2 >= 0 && v > v_const) {
						dn[k] = v - v_const;
					}
					else if (v <= 0 && v2 < -v_const) {
						dn[k] = v2 + v_const;
					}
					else if (v > .25f && v2 < -v_const) {
						dn[k] = (v + v2) / 2;
					}
					if (dn[k] != dc[k]) {
						ch++;
					}
				}
			}
		}
	}

	for (int k = 0; k < d_size; k++) {
		if (dn[k] >= -5) {
			dn[k] = 1;
		}
		else {
			dn[k] = 0;
		}
	}

	dn[d_nx * d_ny / 2 + d_nx / 2] = 2;
	ch = 1;
	while (ch) {
		ch = 0;
		for (int l = 1; l < d_ny - 1; l++) {
			for (int j = 1; j < d_nx - 1; j++) {
				int k = j + l * d_nx;
				if (dn[k] == 1) {
					float v = Util::get_max(dn[k - 1], dn[k + 1], dn[k - d_nx], dn[k + d_nx]);
					if (v == 2.0f) {
						dn[k] = 2.0f;
						ch = 1;
					}
				}
			}
		}
	}

	for (ch = 0; ch < 4; ch++) {
		memcpy(dc, dn, d_nx * d_ny * sizeof(float));

		for (int l = 1; l < d_ny - 1; l++) {
			for (int j = 1; j < d_nx - 1; j++) {
				int k = j + l * d_nx;

				if (dc[k] != 2.0f) {
					float v = Util::get_max(dc[k - 1], dc[k + 1], dc[k - d_nx], dc[k + d_nx]);
					if (v == 2.0f) {
						dn[k] = 2.0f;
					}
				}
			}
		}
	}

	for (int k = 0; k < d_size; k++) {
		if (dn[k] == 2.0f) {
			dn[k] = 1;
		}
		else {
			dn[k] = 0;
		}
	}

	memcpy(dd, dn, d_size * sizeof(float));
	if( dn )
	{
		free(dn);
		dn = 0;
	}

	if( dc )
	{
		free(dc);
		dc = 0;
	}

	d->update();

	image->mult(*d);
	if( d )
	{
		delete d;
		d = 0;
	}
}


void AddRandomNoiseProcessor::process_inplace(EMData * image)
{
	if (!image) {
		LOGWARN("NULL Image");
		return;
	}
	
	if (!image->is_complex()) {
		LOGERR("AddRandomNoise Processor only works for complex image");
		throw ImageFormatException("only work for complex image");
	}

	int n = params["n"];
	float x0 = params["x0"];
	float dx = params["dx"];
	vector < float >y = params["y"];

	int interpolation = 1;
	if (params.has_key("interpolation")) {
		interpolation = params["interpolation"];
	}

	if(params.has_key("seed")) {
#ifdef _WIN32
		srand((int)params["seed"]);
#else
		srandom((int)params["seed"]);
#endif 
	}
	else {
#ifdef _WIN32
		srand(time(0));
#else
		srandom(time(0));
#endif 
	}

	int nx = image->get_xsize();
	int ny = image->get_ysize();
	int nz = image->get_zsize();

	image->ap2ri();
	float *rdata = image->get_data();

	int k = 0;
	float half_nz = 0;
	if (nz > 1) {
		half_nz = nz / 2.0f;
	}

	const float sqrt_2 = sqrt((float) 2);

	for (int h = 0; h < nz; h++) {
		for (int j = 0; j < ny; j++) {
			for (int i = 0; i < nx; i += 2, k += 2) {
				float r = sqrt(Util::hypot3(i / 2.0f, j - ny / 2.0f, h - half_nz));
				r = (r - x0) / dx;
				int l = 0;
				if (interpolation) {
					l = (int) floor(r);
				}
				else {
					l = (int) floor(r + 0.5f);
				}
				r -= l;
				float f = 0;
				if (l >= n - 2) {
					f = y[n - 1];
				}
				else if (l < 0) {
					l = 0;
				}
				else {
					if (interpolation) {
						f = (y[l] * (1 - r) + y[l + 1] * r);
					}
					else {
						f = y[l];
					}
				}
				f = Util::get_gauss_rand(sqrt(f), sqrt(f) / 3);
				float a = Util::get_frand(0.0f, (float)(2 * M_PI));
				if (i == 0) {
					f *= sqrt_2;
				}
				rdata[k] += f * cos(a);
				rdata[k + 1] += f * sin(a);
			}
		}
	}

	image->update();
}

void AddMaskShellProcessor::process_inplace(EMData * image)
{
	if (!image) {
		LOGWARN("NULL Image");
		return;
	}

	int nx = image->get_xsize();
	int ny = image->get_ysize();
	int nz = image->get_zsize();

	if (ny == 1) {
		LOGERR("Tried to add mask shell to 1d image");
		return;
	}

	int num_shells = params["nshells"];

	float *d = image->get_data();
	float k = 0.99999f;
	int nxy = nx * ny;

	if (nz == 1) {
		for (int i = 0; i < num_shells; i++) {
			for (int y = 1; y < ny - 1; y++) {
				int cur_y = y * nx;

				for (int x = 1; x < nx - 1; x++) {
					int j = x + cur_y;
					if (!d[j] && (d[j - 1] > k || d[j + 1] > k || d[j + nx] > k || d[j - nx] > k)) {
						d[j] = k;
					}
				}
			}
			k -= 0.00001f;
		}
	}
	else {
		for (int i = 0; i < num_shells; i++) {
			for (int z = 1; z < nz - 1; z++) {
				int cur_z = z * nx * ny;

				for (int y = 1; y < ny - 1; y++) {
					int cur_y = y * nx + cur_z;

					for (int x = 1; x < nx - 1; x++) {
						int j = x + cur_y;

						if (!d[j] && (d[j - 1] > k || d[j + 1] > k || d[j + nx] > k ||
									  d[j - nx] > k || d[j - nxy] > k || d[j + nxy] > k)) {
							d[j] = k;
						}
					}
				}
			}

			k -= 0.00001f;
		}
	}

	size_t size = nx * ny * nz;
	for (size_t i = 0; i < size; i++) {
		if (d[i]) {
			d[i] = 1;
		}
		else {
			d[i] = 0;
		}
	}

	image->update();
}

void ToMassCenterProcessor::process_inplace(EMData * image)
{
	if (!image) {
		LOGWARN("NULL Image");
		return;
	}

	int int_shift_only = params["int_shift_only"];
	image->process_inplace("normalize");

	float *rdata = image->get_data();
	int nx = image->get_xsize();
	int ny = image->get_ysize();
	int nz = image->get_zsize();

	float sigma = image->get_attr("sigma");
	float mean = image->get_attr("mean");
	float xm = 0;
	float ym = 0;
	float zm = 0;
	float m = 0;
	int nxy = nx * ny;

	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny; j++) {
			int j2 = nx * j;
			for (int k = 0; k < nz; k++) {
				int l = i + j2 + k * nxy;
				if (rdata[l] >= sigma * .75 + mean) {
					xm += i * rdata[l];
					ym += j * rdata[l];
					zm += k * rdata[l];
					m += rdata[l];
				}
			}
		}
	}

	image->update();

	xm /= m;
	ym /= m;
	zm /= m;

	float dx = 0;
	float dy = 0;
	float dz = 0;

	if (int_shift_only) {
		dx = -(floor(xm + 0.5f) - nx / 2);
		dy = -(floor(ym + 0.5f) - ny / 2);
		dz = 0;
		if (nz > 1) {
			dz = -(floor(zm + 0.5f) - nz / 2);
		}

	}
	else {
		dx = -(xm - nx / 2);
		dy = -(ym - ny / 2);
		if (nz > 1) {
			dz = -(zm - nz / 2);
		}
	}

	image->translate(dx, dy, dz);
}

void ACFCenterProcessor::process_inplace(EMData * image)
{
	if (!image) {
		LOGWARN("NULL Image");
		return;
	}
	
	Dict params1;
	params1["intonly"] = 1;
	params1["maxshift"] = image->get_xsize() / 4;
	image->align("translational", 0, params1);

}

void SNRProcessor::process_inplace(EMData * image)
{
	if (!image) {
		return;
	}

	int wiener = params["wiener"];
	const char *snrfile = params["snrfile"];

	XYData sf;
	int err = sf.read_file(snrfile);
	if (err) {
		LOGERR("couldn't read structure factor file!");
		return;
	}


	for (size_t i = 0; i < sf.get_size(); i++) {
		if (sf.get_y(i) <= 0) {
			sf.set_y(i, -4.0f);
		}
		else {
			sf.set_y(i, log10(sf.get_y(i)));
		}
	}
	sf.update();

	Ctf *image_ctf = image->get_ctf();

	vector < float >ctf;
	if (wiener) {
		ctf = image_ctf->compute_1d(image->get_ysize(), Ctf::CTF_WIENER_CTF_CORRECTION1, &sf);
	}
	else {
		ctf = image_ctf->compute_1d(image->get_ysize(), Ctf::CTF_ABS_SNR, &sf);
	}

	image->process_inplace("normalize.circlemean");

	int nx = image->get_xsize();
	int ny = image->get_ysize();

	Region clip_r(-nx / 2, -ny / 2, nx * 2, ny * 2);
	EMData *d3 = image->get_clip(clip_r);
	EMData *d2 = d3->do_fft();

	d2->apply_radial_func(0, 2.0f / Ctf::CTFOS, ctf, 0);

	if( d3 )
	{
		delete d3;
		d3 = 0;
	}

	if( image )
	{
		delete image;
		image = 0;
	}

	EMData *d1 = d2->do_ift();
	int d1_nx = d1->get_xsize();
	int d1_ny = d1->get_ysize();
	Region d1_r(d1_nx / 4, d1_ny / 4, d1_nx / 2, d1_ny / 2);

	image = d1->get_clip(d1_r);

	if( d1 )
	{
		delete d1;
		d1 = 0;
	}

	if( d2 )
	{
		delete d2;
		d2 = 0;
	}
}

void FileFourierProcessor::process_inplace(EMData * image)
{
	if (!image) {
		LOGWARN("NULL Image");
		return;
	}
	const char *filename = params["filename"];
	float apix = params["apix"];

	FILE *in = fopen(filename, "rb");
	if (!in) {
		LOGERR("FileFourierProcessor: cannot open file '%s'", filename);
		return;
	}

	float f = 0;
	int n = 0;
	while (fscanf(in, " %f %f", &f, &f) == 2) {
		n++;
	}
	rewind(in);

	vector < float >xd(n);
	vector < float >yd(n);

	float sf = apix * image->get_xsize();

	for (int i = 0; fscanf(in, " %f %f", &xd[i], &yd[i]) == 2; i++) {
		xd[i] *= sf;
	}

	if (xd[2] - xd[1] != xd[1] - xd[0]) {
		LOGWARN("Warning, x spacing appears nonuniform %g!=%g\n",
				xd[2] - xd[1], xd[1] - xd[0]);
	}

	EMData *d2 = image->do_fft();
	if( image )
	{
		delete image;
		image = 0;
	}

	d2->apply_radial_func(xd[0], xd[1] - xd[0], yd, 1);
	image = d2->do_ift();
}

void LocalNormProcessor::process_inplace(EMData * image)
{
	if (!image) {
		LOGWARN("NULL Image");
		return;
	}
	float apix = params["apix"];
	float threshold = params["threshold"];
	float radius = params["radius"];

	if (apix > 0) {
		int ny = image->get_ysize();
		radius = ny * apix / radius;
		//printf("Norm filter radius=%1.1f\n", radius);
	}

	EMData *blur = image->copy();
	EMData *maskblur = image->copy();

	maskblur->process_inplace("threshold.binary", Dict("value", threshold));
	maskblur->process_inplace("eman1.filter.lowpass.gaussian", Dict("lowpass", radius));
	maskblur->process_inplace("eman1.filter.highpass.tanh", Dict("highpass", -10.0f));
	maskblur->process_inplace("threshold.belowtozero", Dict("minval", 0.001f));
	maskblur->process_inplace("threshold.belowtozero", Dict("minval", 0.001f));


	blur->process_inplace("threshold.belowtozero", Dict("minval", threshold));
	blur->process_inplace("eman1.filter.lowpass.gaussian", Dict("lowpass", radius));
	blur->process_inplace("eman1.filter.highpass.tanh", Dict("highpass", -10.0f));

	maskblur->div(*blur);
	image->mult(*maskblur);
	maskblur->write_image("norm.mrc", 0, EMUtil::IMAGE_MRC);

	if( maskblur )
	{
		delete maskblur;
		maskblur = 0;
	}

	if( blur )
	{
		delete blur;
		blur = 0;
	}
}


void SymSearchProcessor::process_inplace(EMData * image)
{
	if (!image) {
		LOGWARN("NULL Image");
		return;
	}
	float thresh = params["thresh"];
	int output_symlabel = params["output_symlabel"];

	// set up all the symmetry transforms for all the searched symmetries
	const vector<string> sym_list = params["sym"];
	int sym_num = sym_list.size();
	vector< vector< Transform3D > > transforms(sym_num);
	vector< float* > symvals(sym_num);
	for (int i =0; i < sym_num; i++) {
		Transform3D r;
		int nsym = r.get_nsym(sym_list[i]);
		vector<Transform3D> sym_transform(nsym);
		for (int s=0; s<nsym; s++) {
			sym_transform[s] = r.get_sym(sym_list[i], s);
		}
		transforms[i] = sym_transform;
		symvals[i] = new float[nsym]; // new float(nsym);
	}
	
	EMData *orig = image->copy();

	image->to_zero();
	
	int nx= image->get_xsize();
	int ny= image->get_ysize();
	int nz= image->get_zsize();
	int xy = nx * ny;
	float * data = image->get_data();
	float * sdata = orig->get_data();

	EMData *symlabel = 0;
	float * ldata = symlabel->get_data();
	if (output_symlabel) {
		symlabel = image->copy();
		symlabel->to_zero();
		ldata = symlabel->get_data();
	}
	
	for (int k = 0; k < nz; k++) {
		for (int j = 0; j < ny; j++) {
			for(int i = 0; i < nx; i++) {
				long index = k * nx * ny + j * nx + i;
				float val = sdata[ index ];
				float bestmean = val, bestsymlevel = FLT_MAX;
				int bestsym = 0;
				for( int sym = 0; sym< sym_num; sym++) {
					int cur_sym_num = transforms[sym].size();
					float *symval = symvals[sym];
					// first find out all the symmetry related location values
					for( int s = 0; s < cur_sym_num; s++){
						Transform3D r = transforms[sym][s];
						float x2 = (float)(r[0][0] * (i-nx/2) + r[0][1] * (j-ny/2) + r[0][2] * (k-nz/2) + nx / 2);
						float y2 = (float)(r[1][0] * (i-nx/2) + r[1][1] * (j-ny/2) + r[1][2] * (k-nz/2) + ny / 2);
						float z2 = (float)(r[2][0] * (i-nx/2) + r[2][1] * (j-ny/2) + r[2][2] * (k-nz/2) + nz / 2);
		
						if (x2 >= 0 && y2 >= 0 && z2 >= 0 && x2 < (nx - 1) && y2 < (ny - 1)
							&& z2 < (nz - 1)) {
							float x = (float)Util::fast_floor(x2);
							float y = (float)Util::fast_floor(y2);
							float z = (float)Util::fast_floor(z2);
		
							float t = x2 - x;
							float u = y2 - y;
							float v = z2 - z;
		
							int ii = (int) (x + y * nx + z * xy);
		
							symval[s]=
								Util::trilinear_interpolate(sdata[ii], sdata[ii + 1], sdata[ii + nx],
															sdata[ii + nx + 1], sdata[ii + nx * ny],
															sdata[ii + xy + 1], sdata[ii + xy + nx],
															sdata[ii + xy + nx + 1], t, u, v);
						}
						else {
							symval[s] = 0.0 ;
						}
					}
					float tmean=0, tsigma=0;
					for( int s = 0; s < cur_sym_num; s++) {
						tmean += symval[s];
						tsigma += symval[s] * symval[s];
					}
					tmean /= cur_sym_num;
					tsigma = tsigma/cur_sym_num - tmean*tmean;
					if (tsigma < bestsymlevel ) {
						bestsymlevel = tsigma;
						bestmean = tmean;
						bestsym = sym;
					}
				}
				if ( bestsymlevel > thresh) {
					if (output_symlabel) ldata[index] = (float)bestsym;
					data[index] = bestmean;
				}
				else {
					if (output_symlabel) ldata[index] = -1;
					data[index] = val;
				}
			}
		}
	}
	if( orig )
	{
		delete orig;
		orig = 0;
	}
	for (int i =0; i < sym_num; i++) {
		if( symvals[i] )
		{
			delete symvals[i];
			symvals[i] = 0;
		}
	}
	if (symlabel) params.put("symlabel_map", EMObject(symlabel));
}


void IndexMaskFileProcessor::process_inplace(EMData * image)
{
	if (!image) {
		LOGWARN("NULL Image");
		return;
	}

	const char *filename = params["filename"];
	EMData *msk = new EMData();
	msk->read_image(filename);
	if (!EMUtil::is_same_size(image, msk)) {
		LOGERR("IndexMaskFileProcessor: Mask size different than image");
		return;
	}

	if ((int) params["ismaskset"] != 0) {
		msk->process_inplace("threshold.binaryrange", Dict("low", 0.5f, "high", 1.5f));
	}

	image->mult(*msk);
	if( msk )
	{
		delete msk;
		msk = 0;
	}
}


void CoordinateMaskFileProcessor::process_inplace(EMData * image)
{
	if (!image) {
		LOGWARN("NULL Image");
		return;
	}

	const char *filename = params["filename"];
	EMData *msk = new EMData();
	msk->read_image(filename);

	int nx = image->get_xsize();
	int ny = image->get_ysize();
	int nz = image->get_zsize();

	int xm = msk->get_xsize();
	int ym = msk->get_ysize();
	int zm = msk->get_zsize();

	float apix = image->get_attr("apix_x");
	float apixm = msk->get_attr("apix_x");

	float xo = image->get_attr("origin_row");
	float yo = image->get_attr("origin_col");
	float zo = image->get_attr("origin_sec");

	float xom = msk->get_attr("origin_row");
	float yom = msk->get_attr("origin_col");
	float zom = msk->get_attr("origin_sec");

	float *dp = image->get_data();
	float *dpm = msk->get_data();
	int nxy = nx * ny;

	for (int k = 0; k < nz; k++) {
		float zc = zo + k * apix;
		if (zc <= zom || zc >= zom + zm * apixm) {
			memset(&(dp[k * nxy]), 0, sizeof(float) * nxy);
		}
		else {
			int km = (int) ((zc - zom) / apixm);

			for (int j = 0; j < ny; j++) {
				float yc = yo + j * apix;
				if (yc <= yom || yc >= yom + ym * apixm) {
					memset(&(dp[k * nxy + j * nx]), 0, sizeof(float) * nx);
				}
				else {
					int jm = (int) ((yc - yom) / apixm);

					for (int i = 0; i < nx; i++) {
						float xc = xo + i * apix;
						if (xc <= xom || xc >= xom + xm * apixm) {
							dp[k * nxy + j * nx + i] = 0;
						}
						else {
							int im = (int) ((xc - xom) / apixm);
							if (dpm[km * xm * ym + jm * xm + im] <= 0) {
								dp[k * nxy + j * nx + i] = 0;
							}
						}
					}
				}
			}
		}
	}

	image->update();
	msk->update();
	if( msk )
	{
		delete msk;
		msk = 0;
	}
}


void SetSFProcessor::process_inplace(EMData * image)
{
	if (!image) {
		LOGWARN("NULL Image");
		return;
	}

	const char *sffile = params["filename"];
	float apix = params["apix"];

	if (apix <= 0) {
		LOGERR("Must specify apix with setsf");
		return;
	}

	int nx = image->get_xsize();
	int ny = image->get_ysize();
	int nz = image->get_zsize();

	if (nx != ny || ny != nz) {
		LOGERR("SetSF processor works only on even cubic images");
		return;
	}

	XYData sf;
	int err = sf.read_file(sffile);
	if (err) {
		LOGERR("couldn't read structure factor file!");
		return;
	}

	EMData *dataf = image->do_fft();

	vector < float >curve = dataf->calc_radial_dist(nx, 0, 0.5f,1);
	float step = 1.0f / (apix * 2.0f * nx);
	Util::save_data(0, step, curve, "proc3d.raddist");

	double sum = 0;
	for (int i = 0; i < nx; i++) {
		if (curve[i] > 0) {
			curve[i] = sqrt(sf.get_yatx(i * step) / curve[i]);
		}
		else {
			curve[i] = 0;
		}
		sum += curve[i];
	}

	float avg = (float) sum / nx;
	for (int i = 0; i < nx; i++) {
		curve[i] /= avg;
	}

	float lv2 = curve[0];

	for (int i = 1; i < nx - 1; i++) {
		float lv = curve[i];
		curve[i] = 0.25f * lv2 + 0.5f * curve[i] + 0.25f * curve[i + 1];
		lv2 = lv;
	}
	Util::save_data(0, step, curve, "proc3d.filt");

	dataf->apply_radial_func(0, 0.5f, curve);
	if( image )
	{
		delete image;
		image = 0;
	}

	image = dataf->do_ift();
}




void SmartMaskProcessor::process_inplace(EMData * image)
{
	if (!image) {
		LOGWARN("NULL Image");
		return;
	}

	float mask = params["mask"];

	int nx = image->get_xsize();
	int ny = image->get_ysize();
	int nz = image->get_zsize();

	float *dat = image->get_data();
	double sma = 0;
	int smn = 0;
	for (int k = 0; k < nz; k++) {
		for (int j = 0; j < ny; j++) {
			for (int i = 0; i < nx; i++, dat++) {
				float r =
					sqrt((float) Util::square(i - nx / 2) + Util::square(j - ny / 2) +
						 Util::square(k - nz / 2));
				if (r > mask - 1.5f && r < mask - 0.5f) {
					sma += *dat;
					smn++;
				}
			}
		}
	}

	float smask = (float) sma / smn;
	image->update();

	dat = image->get_data();
	for (int k = 0; k < nz; k++) {
		for (int j = 0; j < ny; j++) {
			for (int i = 0; i < nx; i++, dat++) {
				float r =
					sqrt((float) Util::square(i - nx / 2) + Util::square(j - ny / 2) +
						 Util::square(k - nz / 2));
				if (r > mask - .5) {
					*dat = 0;
				}
				else {
					*dat -= smask;
				}
			}
		}
	}

	image->update();
}

void AutoMask3DProcessor::search_nearby(float *dat, float *dat2, int nx, int ny, int nz, float threshold)
{
	Assert(dat != 0);
	Assert(dat2 != 0);
	Assert(nx > 0);
	Assert(ny > 0);
	
	bool done = false;
	int nxy = nx * ny;

	while (!done) {
		done = true;
		for (int k = 1; k < nz - 1; k++) {
			int k2 = k * nxy;
			for (int j = 1; j < ny - 1; j++) {
				int l = j * nx + k2 + 1;

				for (int i = 1; i < nx - 1; i++) {
					if (dat[l] >= threshold && dat2[l]) {
						if (dat2[l - 1] || dat2[l + 1] ||
							dat2[l - nx] || dat2[l + nx] || dat2[l - nxy] || dat2[l + nxy]) {
							dat2[l] = 1.0f;
							done = false;
						}
					}
					l++;
				}
			}
		}
	}
}

void AutoMask3DProcessor::fill_nearby(float *dat2, int nx, int ny, int nz)
{
	Assert(dat2 != 0);
	Assert(nx > 0);
	Assert(ny > 0);
	Assert(nz >= 0);
	
	int nxy = nx * ny;

	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny; j++) {
			int j2 = j * nx + i;
			int k0 = 0;
			for (int k = 0; k < nz; k++) {
				if (dat2[j2 + k * nxy]) {
					k0 = k;
					break;
				}
			}

			if (k0 != nz) {
				int k1 = nz - 1;
				for (int k = nz - 1; k >= 0; k--) {
					if (dat2[j2 + k * nxy]) {
						k1 = k;
						break;
					}
				}

				for (int k = k0 + 1; k < k1; k++) {
					dat2[j2 + k * nxy] = 1.0f;
				}
			}
		}
	}

	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < nz; j++) {
			int j2 = j * nxy + i;
			int k0 = 0;
			for (int k = 0; k < ny; k++) {
				if (dat2[k * nx + j2]) {
					k0 = k;
					break;
				}
			}

			if (k0 != ny) {
				int k1 = ny - 1;
				for (int k = ny - 1; k >= 0; k--) {
					if (dat2[k * nx + j2]) {
						k1 = k;
						break;
					}
				}

				for (int k = k0 + 1; k < k1; k++) {
					dat2[k * nx + j2] = 1.0f;
				}
			}
		}
	}

	for (int i = 0; i < ny; i++) {
		for (int j = 0; j < nz; j++) {
			int j2 = i * nx + j * nxy;
			int k0 = 0;
			for (int k = 0; k < nx; k++) {
				if (dat2[k + j2]) {
					k0 = k;
					break;
				}
			}
			if (k0 != nx) {
				int k1 = nx - 1;
				for (int k = nx - 1; k >= 0; k--) {
					if (dat2[k + j2]) {
						k1 = k;
						break;
					}
				}

				for (int k = k0 + 1; k < k1; k++) {
					dat2[k + j2] = 1.0f;
				}
			}
		}
	}

}

void AutoMask3DProcessor::process_inplace(EMData * image)
{
	if (!image) {
		LOGWARN("NULL Image");
		return;
	}

	int nx = image->get_xsize();
	int ny = image->get_ysize();
	int nz = image->get_zsize();

	EMData *amask = new EMData();
	amask->set_size(nx, ny, nz);

	float sig = 0;
	float mean = 0;

	if (params.has_key("threshold1") && params.has_key("threshold2")) {
		sig = image->get_attr("sigma");
		mean = image->get_attr("mean");
	}

	float *dat = image->get_data();
	float *dat2 = amask->get_data();

	float t = 0;
	if (params.has_key("threshold1")) {
		t = params["threshold1"];
	}
	else {
		t = mean + sig * 2.5f;
	}

	int l = 0;
	for (int k = 0; k < nz; k++) {
		for (int j = 0; j < ny; j++) {
			for (int i = 0; i < nx; i++) {
				if (dat[l] > t) {
					dat2[l] = 1.0f;
				}
				l++;
			}
		}
	}


	if (params.has_key("threshold2")) {
		t = params["threshold2"];
	}
	else {
		t = mean + sig * 0.5f;
	}

	search_nearby(dat, dat2, nx, ny, nz, t);

	int nxy = nx * ny;

	for (int k = 1; k < nz - 1; k++) {
		for (int j = 1; j < ny - 1; j++) {
			int l = j * nx + k * nxy + 1;
			for (int i = 1; i < nx - 1; i++, l++) {
				if (dat2[l - 1] == 1.0f || dat2[l + 1] == 1.0f ||
					dat2[l - nx] == 1.0f || dat2[l + nx] == 1.0f ||
					dat2[l - nxy] == 1.0f || dat2[l + nxy] == 1.0f) {
					dat2[l] = 2.0f;
				}
			}
		}
	}

	size_t size = nx * ny * nz;
	for (size_t i = 0; i < size; i++) {
		if (dat2[i] == 2.0f) {
			dat2[i] = 1.0f;
		}
	}

	fill_nearby(dat2, nx, ny, nz);

	image->update();
	amask->update();

	image->mult(*amask);
	amask->write_image("mask.mrc", 0, EMUtil::IMAGE_MRC);
	if( amask )
	{
		delete amask;
		amask = 0;
	}
}


void AutoMask3D2Processor::process_inplace(EMData * image)
{
	if (!image) {
		LOGWARN("NULL Image");
		return;
	}

	int radius = params["radius"];
	float threshold = params["threshold"];
	int nshells = params["nshells"];

	int nx = image->get_xsize();
	int ny = image->get_ysize();
	int nz = image->get_zsize();

	EMData *amask = new EMData();
	amask->set_size(nx, ny, nz);

	float *dat = image->get_data();
	float *dat2 = amask->get_data();

	// start with an initial sphere
	int l = 0;
	for (int k = -nz / 2; k < nz / 2; k++) {
		for (int j = -ny / 2; j < ny / 2; j++) {
			for (int i = -nx / 2; i < nx / 2; i++) {
				if (abs(k) > radius || abs(j) > radius || abs(i) > radius) {
					continue;
				}
				if ( (k * k + j * j + i * i) > (radius*radius) || dat[l] < threshold) {
					continue;
				}
				dat2[l] = 1.0f;
				l++;
			}
		}
	}

	AutoMask3DProcessor::search_nearby(dat, dat2, nx, ny, nz, threshold);

	amask->update();
	float val1 = fabs((float)nshells);
	float val2 = val1 > 2.0f ? 2.0f : 0.0f;

	amask->process_inplace("mask.addshells.gauss", Dict("val1", val1, "val2", val2));

	dat2 = amask->get_data();

	if (nshells < 0) {
		AutoMask3DProcessor::fill_nearby(dat2, nx, ny, nz);
	}

	image->update();
	amask->update();

	EMData *norm = amask->copy();
	norm->process_inplace("threshold.binary");

	EMData *norm2 = norm->copy();
	norm->process_inplace("mask.addshells", Dict("nshells", 1));
	norm->sub(*norm2);
	if( norm2 )
	{
		delete norm2;
		norm2 = 0;
	}

	norm->process_inplace("mask.addshells", Dict("nshells", 2));
	image->process_inplace("normalize.mask", Dict("mask", norm, "no_sigma", 0));
	
	if( norm )
	{
		delete norm;
		norm = 0;
	}

	image->mult(*amask);
	amask->write_image("mask.mrc", 0, EMUtil::IMAGE_MRC);
}

void IterBinMaskProcessor::process_inplace(EMData * image)
{
	if (!image) {
		LOGWARN("NULL Image");
		return;
	}

	float val1 = params["val1"];
	float val2 = params["val2"];

	float *d = image->get_data();
	int nx = image->get_xsize();
	int ny = image->get_ysize();
	int nz = image->get_zsize();

	const int nxy = nx * ny;
	size_t size = nx * ny * nz;
	const float default_val = 100000.0f;


	for (size_t i = 0; i < size; i++) {
		if (d[i] == 0) {
			d[i] = default_val;
		}
	}

	if (nz != 1) {
		for (int l = 1; l <= (int) val1; l++) {
			for (int k = 1; k < nz - 1; k++) {
				int k2 = k * nxy;
				for (int j = 1; j < ny - 1; j++) {
					int j2 = j * nx + k2;
					for (int i = 1; i < nx - 1; i++) {
						int t = i + j2;

						if (d[t - 1] <= l || d[t + 1] <= l || d[t + nx] <= l ||
							d[t - nx] <= l || d[t + nxy] <= l || d[t - nxy] <= l) {
							d[t] = (float) l + 1;
						}
					}
				}
			}
		}
	}
	else {
		for (int l = 1; l <= (int) val1; l++) {
			for (int j = 1; j < ny - 1; j++) {
				for (int i = 1; i < nx - 1; i++) {
					int t = i + j * nx;
					if (d[t - 1] <= l || d[t + 1] <= l || d[t + nx] <= l || d[t - nx] <= l)
						d[t] = (float) l + 1;
				}
			}
		}
	}

	for (size_t i = 0; i < size; i++) {
		if (d[i] == default_val) {
			d[i] = 0;
		}
		else if (d[i] < val2) {
			float f1 = -Util::square(d[i] - val2);
			float f2 = Util::square(2.0f * (val1 - val2));
			d[i] = exp(f1 / f2);
		}
		else {
			d[i] = 1.0f;
		}
	}

	image->update();
}

void TestImageProcessor::preprocess(const EMData * const image)
{
	if (!image) {
		LOGWARN("NULL Image");
		return;
	}
	
	nx = image->get_xsize();
	ny = image->get_ysize();
	nz = image->get_zsize();
}

void TestImageGaussian::process_inplace(EMData * image)
{
	preprocess(image);
	
	float sigma = params["sigma"];
	string axis = (const char*)params["axis"];
	float c = params["c"];
	
	float *dat = image->get_data();
	float r; //this is the distance of pixel from the image center(nx/2, ny/2, nz/2)
	float x2, y2, z2; //this is the coordinates of this pixel from image center
	for (int k = 0; k < nz; k++) {
		for (int j = 0; j < ny; j++) {
			for (int i = 0; i < nx; i++, dat++) {
				x2 = (float)( i - nx/2 );
				y2 = (float)( j - ny/2 );
				z2 = (float)( k - nz/2 );
				
				if(axis==""){
					r = (float)sqrt(x2*x2+y2*y2+z2*z2);
				}
				else if(axis == "x"){
					float lc = -c;
					float rc = c;
					r = ( (float)sqrt((x2-lc)*(x2-lc)+y2*y2+z2*z2) + 
						  (float)sqrt((x2-rc)*(x2-rc)+y2*y2+z2*z2) ) /2.0f - c;
				}
				else if(axis == "y"){
					float lc = -c;
					float rc = c;
					r = ( (float)sqrt(x2*x2+(y2-lc)*(y2-lc)+z2*z2) + 
						  (float)sqrt(x2*x2+(y2-rc)*(y2-rc)+z2*z2) ) /2.0f - c;
				}
				else if(axis == "z"){
					if( nz == 1 ){
						throw InvalidValueException(0, "This is a 2D image, no asymmetric feature for z axis");
					}
					float lc = -c;
					float rc = c;
					r = ( (float)sqrt(x2*x2+y2*y2+(z2-lc)*(z2-lc)) + 
						  (float)sqrt(x2*x2+y2*y2+(z2-rc)*(z2-rc)) ) /2.0f - c;
				}
				else{
					throw InvalidValueException(0, "please specify a valid axis for asymmetric features");
				}
				//the amplitude of the pixel is proportional to the distance of this pixel from the center 
				*dat = (float)gsl_ran_gaussian_pdf((double)r,(double)sigma);
			}
		}
	}
	
	image->update();
}

void TestImageGradient::process_inplace(EMData * image)
{
	string axis = params.set_default("axis", "x");
	
	float m = params.set_default("m", 1.0);
	float b = params.set_default("b", 0.0);
	
	if ( axis != "z" && axis != "y" && axis != "x") throw InvalidParameterException("Axis must be x,y or z");
	
	preprocess(image);
	
	if ( axis == "x")
	{
		for(int k=0; k<nz;++k) {
			for(int j=0; j<ny; ++j) {
				for(int i=0; i <nx; ++i) {
					image->set_value_at(i,j,k,m*i+b);
				}
			}
		}
	}
	else if ( axis == "y")
	{
		for(int k=0; k<nz;++k) {
			for(int j=0; j<ny; ++j) {
				for(int i=0; i <nx; ++i) {
					image->set_value_at(i,j,k,m*j+b);
				}
			}
		}
	}
	else if ( axis == "z")
	{
		for(int k=0; k<nz;++k) {
			for(int j=0; j<ny; ++j) {
				for(int i=0; i <nx; ++i) {
					image->set_value_at(i,j,k,m*k+b);
				}
			}
		}
	}
	image->update();
}

void TestImageAxes::process_inplace(EMData * image)
{
	preprocess(image);
	
	float fill = params.set_default("fill", 1.0f);
	// get the central coordinates
	int cx = nx/2;
	int cy = ny/2;
	int cz = nz/2;
	
	// Offsets are used to detect when "the extra pixel" needs to be filled in
	// They are implemented on the assumption that for odd dimensions
	// the "center pixel" is the center pixel, but for even dimensions the "center
	// pixel" is displaced in the positive direction by 1
	int xoffset = (nx % 2 == 0? 1:0);
	int yoffset = (ny % 2 == 0? 1:0);
	int zoffset = (nz % 2 == 0? 1:0);
	
	// This should never occur - but if indeed it did occur, the code in this function
	// would break - the function would proceed into the final "else" and seg fault
	// It is commented out but left for clarity
// 	if ( nx < 1 || ny < 1 || nz < 1 ) throw ImageDimensionException("Error: one of the image dimensions was less than zero");
	
	if ( nx == 1 && ny == 1 && nz == 1 )
	{
		(*image)(0) = fill;
	}
	else if ( ny == 1 && nz == 1 )
	{
		int radius = params.set_default("radius", cx );
		if ( radius > cx ) radius = cx;
		
		(*image)(cx) = fill;
		for ( int i = 1; i <= radius-xoffset; ++i ) (*image)(cx+i) = fill;
		for ( int i = 1; i <= radius; ++i ) (*image)(cx-i) = fill;
	}
	else if ( nz == 1 )
	{
		int min = ( nx < ny ? nx : ny );
		min /= 2;
		
		int radius = params.set_default("radius", min );
		if ( radius > min ) radius = min;
		
		(*image)(cx,cy) = fill;
		
		for ( int i = 1; i <= radius-xoffset; ++i ) (*image)(cx+i,cy) = fill;
		for ( int i = 1; i <= radius-yoffset; ++i )(*image)(cx,cy+i) = fill;
				
		for ( int i = 1; i <= radius; ++i )
		{
			(*image)(cx-i,cy) = fill;
			(*image)(cx,cy-i) = fill;
		}
		
	}
	else
	{
		// nx > 1 && ny > 1 && nz > 1
		int min = ( nx < ny ? nx : ny );
		if (nz < min ) min = nz;
		min /= 2;
		
		int radius = params.set_default("radius", min);
		if ( radius > min ) radius = min;

		
		(*image)(cx,cy,cz) = fill;
		for ( int i = 1; i <=radius-xoffset; ++i ) (*image)(cx+i,cy,cz) = fill;
		for ( int i = 1; i <=radius-yoffset; ++i ) (*image)(cx,cy+i,cz) = fill;
		for ( int i = 1; i <=radius-zoffset; ++i ) (*image)(cx,cy,cz+i) = fill;
		for ( int i = 1; i <= radius; ++i )
		{
			(*image)(cx-i,cy,cz) = fill;
			(*image)(cx,cy-i,cz) = fill;
			(*image)(cx,cy,cz-i) = fill;
		}
	}
	
	image->update();
}

void TestImageScurve::process_inplace(EMData * image)
{
	preprocess(image);
	
	int dim_size = image->get_ndim();
	if( 2 != dim_size ) {
		throw ImageDimensionException("work for 2D image only");
	} 
	
	int nx = image->get_xsize();
	int ny = image->get_ysize();
	image->to_zero();
	
	for (int i=0; i<100; i++) {
		int x=static_cast<int>( nx/2+nx/6.0*sin(i*2.0*3.14159/100.0) );
		int y=ny/4+i*ny/200;
		for (int xx=x-nx/10; xx<x+nx/10; xx++) {
			for (int yy=y-ny/10; yy<y+ny/10; yy++) {
				(*image)(xx,yy)+=exp(-pow(static_cast<float>(hypot(xx-x,yy-y))*30.0f/nx,2.0f))*(sin(static_cast<float>((xx-x)*(yy-y)))+.5f);
			}
		}
	}
	
	image->update();
}

void TestImagePureGaussian::process_inplace(EMData * image)
{
	preprocess(image);
	
	float x_sigma = params["x_sigma"];
	float y_sigma = params["y_sigma"];
	float z_sigma = params["z_sigma"];

	float x_center = params["x_center"];
	float y_center = params["y_center"];
	float z_center = params["z_center"];

	int nx = image->get_xsize();
	int ny = image->get_ysize();
	int nz = image->get_zsize();

	float x_twosig2 = 2*x_sigma*x_sigma;
	float y_twosig2 = 2*y_sigma*y_sigma;
	float z_twosig2 = 2*z_sigma*z_sigma;

	float sr2pi = sqrt( 2.0f*(float)pi );
	float norm  = 1.0f/ ( x_sigma*y_sigma*z_sigma*sr2pi*sr2pi*sr2pi );

	for (int iz=0; iz < nz; iz++) {
		float z = static_cast<float>(iz) - z_center;
		for (int iy=0; iy < ny; iy++) {
			float y = static_cast<float>(iy) - y_center;
			for (int ix=0; ix < nx; ix++) {
				float x = static_cast<float>(ix) - x_center;
				float sum = x*x/x_twosig2 + y*y/y_twosig2 + z*z/z_twosig2;
				float val = norm*exp(-sum);
				(*image)(ix,iy,iz) = val;
			}
		}
	}
	image->update();
}

void TestImageSinewave::process_inplace(EMData * image)
{
	preprocess(image);
	
	if(!params.has_key("wave_length")) {
		LOGERR("%s wave_length is required parameter", get_name().c_str());
		throw InvalidParameterException("wave_length parameter is required.");
	}
	float wave_length = params["wave_length"];
	
	string axis = "";
	if(params.has_key("axis")) {
		axis = (const char*)params["axis"];
	}
	
	float phase = 0;
	if(params.has_key("phase")) {
		phase = params["phase"];
	}
	
	int ndim = image->get_ndim();
	float * dat = image->get_data();
	
	if(ndim==1) {	//1D
		for(int i=0; i<nx; ++i, ++dat) {
			*dat = sin(i*(2.0f*M_PI/wave_length) - phase*180/M_PI);
		}
	}
	else if(ndim==2) {	//2D
		float alpha = 0;
		if(params.has_key("az")) {
			alpha = params["az"];
		}
		for(int j=0; j<ny; ++j) {
			for(int i=0; i<nx; ++i, ++dat) {
				if(alpha != 0) {
					*dat = sin((i*sin((180-alpha)*M_PI/180)+j*cos((180-alpha)*M_PI/180))*(2.0f*M_PI/wave_length) - phase*M_PI/180); 
				}
				else if(axis.compare("y")==0 || axis.compare("Y")==0) {
					*dat = sin(j*(2.0f*M_PI/wave_length) - phase*M_PI/180);
				}
				else {
					*dat = sin(i*(2.0f*M_PI/wave_length) - phase*M_PI/180);
				}
			} 
		}
	}
	else {	//3D 
		float az = 0;
		if(params.has_key("az")) {
			az = params["az"];
		}
		float alt = 0;
		if(params.has_key("alt")) {
			alt = params["alt"];
		}
		float phi = 0;
		if(params.has_key("phi")) {
			phi = params["phi"];
		}
		
		for(int k=0; k<nz; ++k) {
			for(int j=0; j<ny; ++j) {
				for(int i=0; i<nx; ++i, ++dat) {
					if(axis.compare("z")==0 || axis.compare("Z")==0) {
						*dat = sin(k*(2.0f*M_PI/wave_length) - phase*M_PI/180);
					}
					else if(axis.compare("y")==0 || axis.compare("Y")==0) {
						*dat = sin(j*(2.0f*M_PI/wave_length) - phase*M_PI/180);
					}
					else {
						*dat = sin(i*(2.0f*M_PI/wave_length) - phase*M_PI/180);
					}
				}
			}
		}
		
		if(az != 0 || alt != 0 || phi != 0) {
			image->rotate(az, alt, phi);
		}
	}
	
	image->update();	
}

void TestImageSinewaveCircular::process_inplace(EMData * image)
{
	preprocess(image);
	
	float wave_length = params["wave_length"];
	string axis = (const char*)params["axis"];
	float c = params["c"];
	float phase = params["phase"];
	
	float *dat = image->get_data();
	float r; //this is the distance of pixel from the image center(nx/2, ny/2, nz/2)
	float x2, y2, z2; //this is the coordinates of this pixel from image center
	for (int k = 0; k < nz; k++) {
		for (int j = 0; j < ny; j++) {
			for (int i = 0; i < nx; i++, dat++) {
				x2 = (float)( i - nx/2 );
				y2 = (float)( j - ny/2 );
				z2 = (float)( k - nz/2 );
				if(axis == ""){
					r = (float)sqrt(x2*x2+y2*y2+z2*z2);
				}
				else if(axis == "x"){
					float lc = -c;
					float rc = c;
					r = ( (float)sqrt((x2-lc)*(x2-lc)+y2*y2+z2*z2) + 
						  (float)sqrt((x2-rc)*(x2-rc)+y2*y2+z2*z2) ) /2.0f - c;
				}
				else if(axis == "y"){
					float lc = -c;
					float rc = c;
					r = ( (float)sqrt(x2*x2+(y2-lc)*(y2-lc)+z2*z2) + 
						  (float)sqrt(x2*x2+(y2-rc)*(y2-rc)+z2*z2) ) /2.0f - c;
				}
				else if(axis == "z"){
					if( nz == 1 ){
						throw InvalidValueException(0, "This is a 2D image, no asymmetric feature for z axis");
					}
					float lc = -c;
					float rc = c;
					r = ( (float)sqrt(x2*x2+y2*y2+(z2-lc)*(z2-lc)) + 
						  (float)sqrt(x2*x2+y2*y2+(z2-rc)*(z2-rc)) ) /2.0f - c;
				}
				else{
					throw InvalidValueException(0, "please specify a valid axis for asymmetric features");
				}
				*dat = sin( r * (2.0f*M_PI/wave_length) - phase*180/M_PI);
			}
		}
	}
	
	image->update();	
}

void TestImageSquarecube::process_inplace(EMData * image)
{
	preprocess(image);
	
	float edge_length = params["edge_length"];
	string axis = (const char*)params["axis"];
	float odd_edge = params["odd_edge"];
	int fill = (int)params["fill"];
	
	float *dat = image->get_data();
	float x2, y2, z2; //this coordinates of this pixel from center
	float xEdge, yEdge, zEdge; //half of edge length for this cube
	if(axis == ""){
		xEdge = edge_length/2.0f;
		yEdge = edge_length/2.0f;
		zEdge = edge_length/2.0f;
	}
	else if(axis == "x"){
		xEdge = odd_edge/2.0f;
		yEdge = edge_length/2.0f;
		zEdge = edge_length/2.0f;
	}
	else if(axis == "y"){
		xEdge = edge_length/2.0f;
		yEdge = odd_edge/2.0f;
		zEdge = edge_length/2.0f;
	}
	else if(axis == "z"){
		if( nz == 1 ){
			throw InvalidValueException(0, "This is a 2D image, no asymmetric feature for z axis");
		}
		xEdge = edge_length/2.0f;
		yEdge = edge_length/2.0f;
		zEdge = odd_edge/2.0f;
	}
	else{
		throw InvalidValueException(0, "please specify a valid axis for asymmetric features");
	}
	for (int k = 0; k < nz; k++) {
		for (int j = 0; j < ny; j++) {
			for (int i = 0; i < nx; i++, dat++) {
				x2 = (float)fabs((float)i - nx/2);
				y2 = (float)fabs((float)j - ny/2);
				z2 = (float)fabs((float)k - nz/2);
				if( x2<=xEdge && y2<=yEdge && z2<=zEdge ) {
					if( !fill) {
						*dat = 0;
					} 
					else {
						*dat = 1;
					}
				}
				else {
					if( !fill ) {
						*dat = 1;
					}
					else {
						*dat = 0;
					}
				}
			}
		}
	}
	
	image->update();
}

void TestImageCirclesphere::process_inplace(EMData * image)
{
	preprocess(image);
	
	float radius = params.set_default("radius",nx/2.0f);
	string axis = (const char*)params["axis"];
	float c =  params.set_default("c",nx/2.0f);
	int fill = params.set_default("fill",1);
	
	float *dat = image->get_data();
	float x2, y2, z2; //this is coordinates of this pixel from center
	float r = 0.0f;
	float asy = 0.0f;
	if(axis == ""){
		asy = radius;
	}
	else if(axis == "x" || axis == "y"){
		asy = c;
	}
	else if(axis=="z"){
		if( nz == 1 ){
			throw InvalidValueException(0, "This is a 2D image, no asymmetric feature for z axis");
		}			
		asy = c;
	}
	else{
		throw InvalidValueException(0, "please specify a valid axis for asymmetric features");
	}
	
	for (int k = 0; k < nz; k++) {
		for (int j = 0; j < ny; j++) {
			for (int i = 0; i < nx; i++, dat++) {				
				x2 = fabs((float)i - nx/2);
				y2 = fabs((float)j - ny/2);
				z2 = fabs((float)k - nz/2);
				if( axis == "" ){
					r = (x2*x2)/(radius*radius) + (y2*y2)/(radius*radius) + (z2*z2)/(radius*radius);
				}
				else if (axis == "x"){
					r = (x2*x2)/(asy*asy) + (y2*y2)/(radius*radius) + (z2*z2)/(radius*radius);
				}
				else if(axis == "y"){
					r = (x2*x2)/(radius*radius) + (y2*y2)/(asy*asy) + (z2*z2)/(radius*radius);
				}
				else if(axis=="z"){
					r = (x2*x2)/(radius*radius) + (y2*y2)/(radius*radius) + (z2*z2)/(asy*asy);
				}
				if( r<=1 ) {
					if( !fill) {
						*dat = 0;
					} 
					else {
						*dat = 1;
					}
				}
				else {
					if( !fill ) {
						*dat = 1;
					}
					else {
						*dat = 0;
					}
				}
			}
		}
	}
	
	image->update();
}

void TestImageNoiseUniformRand::process_inplace(EMData * image)
{
	preprocess(image);
	
	if(params.has_key("seed")) {
#ifdef _WIN32
		srand((int)params["seed"]);
#else
		srandom((int)params["seed"]);
#endif 
	}
	else {
#ifdef _WIN32
		srand(time(0));
#else 
		srandom(time(0)); //generate a seed by current time
#endif 
	}
	
	float *dat = image->get_data();
	for (int i=0; i<nx*ny*nz; i++) {
		dat[i] = (float)rand()/(float)RAND_MAX;
	}
	
	image->update();
}

void TestImageNoiseGauss::process_inplace(EMData * image)
{
	preprocess(image);
	
	float sigma = params["sigma"];
	if (sigma<=0) sigma = 1.0;

	float mean = params["mean"];
	
	if (params.has_key("seed")) {
#ifdef _WIN32
		srand((int)params["seed"]);
#else
		srandom((int)params["seed"]);
#endif
	}
	
	float *dat = image->get_data();
	for (int i=0; i<nx*ny*nz; i++) {
		dat[i] = Util::get_gauss_rand(mean, sigma);
	}
	
	image->update();
}

void TestImageCylinder::process_inplace(EMData * image)
{
	preprocess(image);
	
	int nx = image->get_xsize();
	int ny = image->get_ysize();
	int nz = image->get_zsize();
	
	if(nz == 1) {
		throw ImageDimensionException("This processor only apply to 3D image");
	} 
	
	float radius = params["radius"];
#ifdef _WIN32
	if(radius > _cpp_min(nx, ny)/2.0) {
#else
	if(radius > std::min(nx, ny)/2.0) {
#endif
		throw InvalidValueException(radius, "radius must be <= min(nx, ny)/2");
	}
	
	float height;
	if(params.has_key("height")) {
		height = params["height"];
		if(height > nz) {
			throw InvalidValueException(radius, "height must be <= nz");
		}
	}
	else {
		height = static_cast<float>(nz);
	}
	
	float *dat = image->get_data();
	float x2, y2; //this is coordinates of this pixel from center axle
	float r = 0.0f;
	for (int k = 0; k < nz; k++) {
		for (int j = 0; j < ny; j++) {
			for (int i = 0; i < nx; i++, dat++) {				
				x2 = fabs((float)i - nx/2);
				y2 = fabs((float)j - ny/2);
				r = (x2*x2)/(radius*radius) + (y2*y2)/(radius*radius);
				
				if(r<=1 && k>=(nz-height)/2 && k<=(nz+height)/2) {
					*dat = 1;
				}
				else {
					*dat = 0;
				}
			}
		}
	}
	
	image->update();
}

void RampProcessor::process_inplace(EMData * image)
{
	if (!image) {
		return;
	}

	int nz = image->get_zsize();
	if (nz > 1) {
		LOGERR("%s Processor doesn't support 3D model", get_name().c_str());
		throw ImageDimensionException("3D model not supported");
	}

	int nsam = image->get_xsize();
	int nrow = image->get_ysize();
	int n1 = nsam / 2;
	double sx1 = double(n1)*double(nsam+1);
	if ( nsam % 2 == 1 ) 
		sx1 += 1 + n1;
	sx1 *= nrow;
	int n2 = nrow / 2;
	double sx2 = double(n2)*double(nrow+1);
	if ( nrow % 2 == 1 ) 
		sx2 += 1 + n2;
	sx2 *= nsam;
	float *data = image->get_data();
	float *row = NULL; // handy pointer for values in a specific row of the data
	// statistical sums
	double syx1 = 0, syx2 = 0, sy = 0, sx1q = 0, sx2q = 0, syq = 0;
	for (int j=1; j <= nrow; j++) {
		row = data + (j-1)*nsam - 1; // "-1" so that we can start counting at 1
		for (int i=1; i<=nsam; i++) {
			syx1 += row[i]*i;
			syx2 += row[i]*j;
			sy += row[i];
			sx1q += i*i;
			sx2q += j*j;
			syq += row[i]*double(row[i]);
		}
	}
	// least-squares
	float dn = float(nsam)*float(nrow);
	double qyx1 = syx1 - sx1*sy / dn;
	double qyx2 = syx2 - sx2*sy / dn;
	double qx1x2 = 0.0;
	double qx1 = sx1q - sx1*sx1 / dn;
	double qx2 = sx2q - sx2*sx2 / dn;
	double qy = syq - sy*sy / dn;
	double c = qx1*qx2 - qx1x2*qx1x2;
	if ( c > FLT_EPSILON ) {
		double b1 = (qyx1*qx2 - qyx2*qx1x2) / c;
		double b2 = (qyx2*qx1 - qyx1*qx1x2) / c;
		double a = (sy - b1*sx1 - b2*sx2) / dn;
		double d = a + b1 + b2;
		for (int i=1; i<=nrow; i++) {
			qy = d;
			row = data + (i-1)*nsam - 1;
			for (int k=1; k<=nsam; k++) {
				row[k] -= static_cast<float>(qy);
				qy += b1;
			}
			d += b2;
		}
	} // image not altered if c is zero

	image->update();
}


void CCDNormProcessor::process_inplace(EMData * image)
{
	if (!image) {
	  Log::logger()->set_level(Log::ERROR_LOG);
	  Log::logger()->error("Null image during call to CCDNorm\n");
	  return;
	}
	if (image->get_zsize() > 1) {
	  Log::logger()->set_level(Log::ERROR_LOG);
	  Log::logger()->error("CCDNorm does not support 3d images\n");	  
	  return;
	}

	int xs = image->get_xsize();
	int ys = image->get_ysize();

	// width of sample area on either side of the seams
	int width = params["width"];

	width%=(xs > ys ? xs/2 : ys/2);  // make sure width is a valid value
	if (width==0) {
	  width=1;
	}

	// get the 4 "seams" of the image
	float *left, *right, *top, *bottom;

	double *temp;
	temp= (double*)malloc((xs > ys ? xs*width : ys*width)*sizeof(double));
	if (temp==NULL) {
	  Log::logger()->set_level(Log::ERROR_LOG);
	  Log::logger()->error("Could not allocate enough memory during call to CCDNorm\n");
	  return;
	}

	int x, y, z;

	// the mean values of each seam and the average
	double mL,mR,mT,mB;

	// how much to shift each seam
	double nl,nr,nt,nb;

	// quad. shifting amount
	double q1,q2,q3,q4;

	// calc. the mean for each quadrant
	for (z=0; z<width; z++) {
	  left = image->get_col(xs/2 -1-z)->get_data();
	  for (y=0; y<ys; y++)
	    temp[z*ys+y]=left[y];
	}
	mL=gsl_stats_mean(temp,1,ys*width);

	for (z=0; z<width; z++) {
	  right = image->get_col(xs/2 +z)->get_data();
	  for (y=0; y<ys; y++)
	    temp[z*ys+y]=right[y];
	}
	mR=gsl_stats_mean(temp,1,ys*width);

	for (z=0; z<width; z++) {
	  top = image->get_row(ys/2 -1-z)->get_data();
	  for (x=0; x<xs; x++)
	    temp[z*xs+x]=top[x];
	}
	mT=gsl_stats_mean(temp,1,xs*width);

	for (z=0; z<width; z++) {
	  bottom = image->get_row(ys/2 +z)->get_data();
	  for (x=0; x<xs; x++)
	    temp[z*xs+x]=bottom[x];
	}
	mB=gsl_stats_mean(temp,1,xs*width);

	free(temp);

	nl=(mL+mR)/2-mL;
	nr=(mL+mR)/2-mR;
	nt=(mT+mB)/2-mT;
	nb=(mT+mB)/2-mB;

	q1=nl+nt;
	q2=nr+nt;
	q3=nr+nb;
	q4=nl+nb;

	// change the pixel values
	for (x = 0; x < xs / 2; x++)
	  for (y = 0; y < ys / 2; y++) {
	    image->set_value_at_fast(x, y, image->get_value_at(x, y) + static_cast<float>(q1));
	  }
	for (x = xs / 2; x < xs; x++)
	  for (y = 0; y < ys / 2; y++) {
	    image->set_value_at_fast(x, y, image->get_value_at(x, y) + static_cast<float>(q2));
	  }
	for (x = xs / 2; x < xs; x++)
	  for (y = ys / 2; y < ys; y++) {
	    image->set_value_at_fast(x, y, image->get_value_at(x, y) + static_cast<float>(q3));
	  }
	for (x = 0; x < xs / 2; x++)
	  for (y = ys / 2; y < ys; y++) {
	    image->set_value_at_fast(x, y, image->get_value_at(x, y) + static_cast<float>(q4));
	  }

}

void WaveletProcessor::process_inplace(EMData *image) 
{
	if (image->get_zsize() != 1) {
			LOGERR("%s Processor doesn't support 3D", get_name().c_str());
			throw ImageDimensionException("3D model not supported");
	}	

	int i,nx,ny;
	const gsl_wavelet_type * T;
	nx=image->get_xsize();
	ny=image->get_ysize();
	
	if (nx != ny && ny!=1) throw ImageDimensionException("Wavelet transform only supports square images");
	float l=log((float)nx)/log(2.0f);
	if (l!=floor(l)) throw ImageDimensionException("Wavelet transform size must be power of 2");

	// Unfortunately GSL works only on double() arrays
	// eventually we should put our own wavelet code in here
	// but this will work for now
	double *cpy = (double *)malloc(nx*ny*sizeof(double));
	
	for (i=0; i<nx*ny; i++) cpy[i]=image->get_value_at(i,0,0);

	string tp = (const char*)params["type"];
	if (tp=="daub") T=gsl_wavelet_daubechies;
	else if (tp=="harr") T=gsl_wavelet_haar;
	else if (tp=="bspl") T=gsl_wavelet_bspline;
	else throw InvalidStringException(tp,"Invalid wavelet name, 'daub', 'harr' or 'bspl'");
	
	int K=(int)params["ord"];
	gsl_wavelet_direction dir;
	if ((int)params["dir"]==1) dir=forward;
	else dir=backward;

	gsl_wavelet *w = gsl_wavelet_alloc(T, K);
	gsl_wavelet_workspace *work = gsl_wavelet_workspace_alloc(nx);

	if (ny==1) gsl_wavelet_transform (w,cpy, 1, nx, dir, work);
	else gsl_wavelet2d_transform (w, cpy, nx,nx,ny, dir, work);

	gsl_wavelet_workspace_free (work);
	gsl_wavelet_free (w);

	for (i=0; i<nx*ny; i++) image->set_value_at_fast(i,0,0,static_cast<float>(cpy[i]));
	
	free(cpy);
}

void FFTProcessor::process_inplace(EMData* image)
{
	if( params.has_key("dir") ) {
		if ((int)params["dir"]==-1) {
			image->do_ift_inplace();
		}
		else {
			image->do_fft_inplace();
		}
	}
}

void RadialProcessor::process_inplace(EMData * image)
{
	if (!image) {
		LOGWARN("NULL Image");
		return;
	}
	
	//Note : real image only!
	if(image->is_complex()) {
		LOGERR("%s Processor only operates on real images", get_name().c_str());
		throw ImageFormatException("apply to real image only");
	}
	
	vector<float> table = params["table"];
	vector<float>::size_type tsize = table.size();
	
	int nx = image->get_xsize();
	int ny = image->get_ysize();
	int nz = image->get_zsize();
	
	int nx2 = nx / 2;
	int ny2 = ny / 2;
	int nz2 = nz / 2;
	float sz[3];
	sz[0] = static_cast<float>(nx2);
	sz[1] = static_cast<float>(ny2);
	sz[2] = static_cast<float>(nz2);
	float szmax = *std::max_element(&sz[0], &sz[3]);
	float maxsize;
	if(nz>1) {
		maxsize = (float)(1.8f * szmax);
	}
	else{
		maxsize = (float)(1.5f * szmax);
	}
	for(int i=tsize+1; i<maxsize; i++) {
		table.push_back(0.0f);
	}
	
	float dx = 1.0f / (float)nx;
	float dy = 1.0f / (float)ny;
	float dz = 1.0f / (float)nz;
	float dx2 = dx * dx;
	float dy2 = dy * dy;
	float dz2 = dz * dz;
	int iz, iy, ix, jz, jy, jx;
	float argx, argy, argz;
	float rf, df, f;
	int ir;
	for(iz=1; iz<=nz; iz++) {
		jz = iz - 1;
		if(jz > nz2) {
			jz -= nz; 
		}
		argz = float(jz*jz) * dz2;
		
		for(iy=1; iy<=ny; iy++) {
			jy = iy - 1;
			if(jy > ny2) {
				jy -= ny;
			}
			argy = argz + float(jy*jy) * dy2;
			
			for(ix=1; ix<=nx; ix++) {
				jx = ix -1;
				argx = argy + float(jx*jx)*dx2;
				
				rf = sqrt(argx)*2.0f*nx2;
				ir = int(rf);
				df = rf - float(ir);
				f = table[ir] + df*(table[ir+1]-table[ir]);
				
				(*image)(ix-1,iy-1,iz-1) *= f;
			}
		}
	}
	
	image->update();
}

void MirrorProcessor::process_inplace(EMData *image) 
{
	if (!image) {
		LOGWARN("NULL Image");
		return;
	}
	
	string axis = (const char*)params["axis"];

	float* data = image->EMData::get_data();
	
	int nx = image->EMData::get_xsize();
	int ny = image->EMData::get_ysize();
	int nz = image->EMData::get_zsize();	 
	int nxy = nx*ny;
	
	int x_start = 1-nx%2;
	int y_start = 1-ny%2;
	int z_start = 1-nz%2;
	
	if (axis == "x" || axis == "X") {
 	        for (int iz = 0; iz < nz; iz++)
			for (int iy = 0; iy < ny; iy++) {
				int offset = nx*iy + nxy*iz;
                         	 reverse(&data[offset+x_start],&data[offset+nx]);
		         }
         }
	 
         if (axis == "y" || axis == "Y") {
                float *tmp = new float[nx];
	         int nhalf   = ny/2;
 		 int beg     = 0;
                 for (int iz = 0; iz < nz; iz++) {
		      beg = iz*nxy;
		      for (int iy = y_start; iy < nhalf; iy++) {
		      		memcpy(tmp, &data[beg+iy*nx], nx*sizeof(float));
				memcpy(&data[beg+iy*nx], &data[beg+(ny-iy-1+y_start)*nx], nx*sizeof(float));
				memcpy(&data[beg+(ny-iy-1+y_start)*nx], tmp, nx*sizeof(float));
			}
		}
		delete[] tmp;
	}			 
	 
         if (axis == "z" || axis == "Z") {
	 	if(1-z_start){
			int nhalf = nz/2;
			float *tmp = new float[nxy];
			for(int iz = 0;iz<nhalf;iz++){
				memcpy(tmp,&data[iz*nxy],nxy*sizeof(float));
				memcpy(&data[iz*nxy],&data[(nz-iz-1)*nxy],nxy*sizeof(float));
				memcpy(&data[(nz-iz-1)*nxy],tmp,nxy*sizeof(float));
				}
			delete[] tmp;
		}
		else{
			float *tmp = new float[nx];
			int nhalf   = nz/2;
			int beg     = 0;
			for (int iy = 0; iy < ny; iy++) {
				beg = iy*nx;
				for (int iz = z_start; iz < nhalf; iz++) {
					memcpy(tmp, &data[beg+ iz*nxy], nx*sizeof(float));
					memcpy(&data[beg+iz*nxy], &data[beg+(nz-iz-1+z_start)*nxy], nx*sizeof(float));
					memcpy(&data[beg+(nz-iz-1+z_start)*nxy], tmp, nx*sizeof(float));
				}
			}
			delete[] tmp;
		}
         }
	 
	 image->update();
}

 
int EMAN::multi_processors(EMData * image, vector < string > processornames)
{
	Assert(image != 0);
	Assert(processornames.size() > 0);
	
	for (size_t i = 0; i < processornames.size(); i++) {
		image->process_inplace(processornames[i]);
	}
	return 0;
}


void ClampingProcessor::process_inplace( EMData* image )
{
	
	if ( image->is_complex() ) throw ImageFormatException("Error: clamping processor does not work on complex images");
	
	float min = params.set_default("minval",default_min);
	float max = params.set_default("maxval",default_max);
	
	// Okay, throwing such an error is probably overkill - but atleast the user will get a loud message
	// saying what went wrong.
	if ( max < min ) throw InvalidParameterException("Error: minval was greater than maxval, aborting");
	
	int size = image->get_xsize()*image->get_ysize()*image->get_zsize();
	for(int i = 0; i < size; ++i )
	{
		float * data = &image->get_data()[i];
		if ( *data < min ) *data = min;
		else if ( *data > max ) *data = max;
	}
}

#define deg2rad 0.017453292519943295
void TestTomoImage::insert_rectangle( EMData* image, const Region& region, const float& value, const Transform3D* const t3d )
{
	int startx = (int)region.origin[0] - (int)region.size[0]/2;
	int starty = (int)region.origin[1] - (int)region.size[1]/2;
	int startz = (int)region.origin[2] - (int)region.size[2]/2;
	
	int endx  = (int)region.origin[0] + (int)region.size[0]/2;
	int endy  = (int)region.origin[1] + (int)region.size[1]/2;
	int endz  = (int)region.origin[2] + (int)region.size[2]/2;
	
	if ( t3d != NULL ) {
		for ( float z = startz; z < endz; z += 0.25 ) {
			for ( float y = starty; y < endy; y += 0.25 ) {
				for ( float x = startx; x < endx; x += 0.25 ) {
					float xt = (float) x - region.origin[0];
					float yt = (float) y - region.origin[1];
					float zt = (float) z - region.origin[2];
					Vec3f v((float)xt,(float)yt,(float)zt);
					v = (*t3d)*v;
					image->set_value_at((int)(v[0]+region.origin[0]),(int)(v[1]+region.origin[1]),(int)(v[2]+region.origin[2]), value);
				}
			}
		}
	} else {	
		for ( int z = startz; z < endz; ++z ) {
			for ( int y = starty; y < endy; ++y ) {
				for ( int x = startx; x < endx; ++x ) {
					image->set_value_at(x,y,z, value);
				}
			}
		}
	}
}

void TestTomoImage::insert_solid_ellipse( EMData* image, const Region& region, const float& value, const Transform3D* const t3d )
{
	int originx = (int)region.origin[0];
	int originy = (int)region.origin[1];
	int originz = (int)region.origin[2];
	
	int radiusx  = (int)region.size[0];
	int radiusy  = (int)region.size[1];
	int radiusz = (int)region.size[2];
	
	int maxrad = ( radiusx > radiusy ? radiusx : radiusy );
	maxrad = ( maxrad > radiusz ? maxrad : radiusz );
	
	float length = radiusx*radiusx + radiusy*radiusy + radiusz*radiusz;
	length = sqrtf(length);
	
	float xinc = radiusx/length/2.0;
	float yinc = radiusy/length/2.0;
	float zinc = radiusz/length/2.0;
	
	for ( float rad = 0; rad <= 2.0*length; rad += 1.0 ) {
		float rx = xinc + rad*xinc;
		float ry = yinc + rad*yinc;
		float rz = zinc + rad*zinc;
		for ( float beta = -90.0; beta <= 90.0; beta += 0.2 ) {
			for ( float lambda = -180.0; lambda <= 180.0; lambda += 0.2 ) {
				float b = deg2rad * beta;
				float l = deg2rad * lambda;
				if ( t3d != NULL ) {
					float z = rz * sin(b);
					float y = ry * cos(b) * sin(l);
					float x = rx * cos(b) * cos(l);
					Vec3f v(x,y,z);
					v = (*t3d)*v;
					image->set_value_at((int)(v[0]+originx),(int)(v[1]+originy),(int)(v[2]+originz), value);
				}
				else {
					int z = (int) (originz + rz * sin(b));
					int y = (int) (originy + ry * cos(b) * sin(l));
					int x = (int) (originx + rx * cos(b) * cos(l));
					image->set_value_at(x,y,z, value);
				}
			}
		}
	}
}

void TestTomoImage::insert_hollow_ellipse( EMData* image, const Region& region, const float& value, const int& radius, const Transform3D* const t3d )
{
	int originx = (int)region.origin[0];
	int originy = (int)region.origin[1];
	int originz = (int)region.origin[2];
	
	int radiusx  = (int)region.size[0];
	int radiusy  = (int)region.size[1];
	int radiusz = (int)region.size[2];
	
	for ( int rad = 0; rad <= radius; ++rad ) {
		int tmprad = rad - radius/2;
		int rx = radiusx + tmprad;
		int ry = radiusy + tmprad;
		int rz = radiusz + tmprad;
		for ( float beta = -90.0; beta <= 90.0; beta += 0.2 ) {
			for ( float lambda = -180.0; lambda <= 180.0; lambda += 0.2 ) {
				float b = deg2rad * beta;
				float l = deg2rad * lambda;
				
				if ( t3d != NULL ) {
					float z = rz * sin(b);
					float y = ry * cos(b) * sin(l);
					float x = rx * cos(b) * cos(l);
					Vec3f v(x,y,z);
					v = (*t3d)*v;
					image->set_value_at((int)(v[0]+originx),(int)(v[1]+originy),(int)(v[2]+originz), value);
				}
				else {
					int z = (int) (originz + rz * sin(b));
					int y = (int) (originy + ry * cos(b) * sin(l));
					int x = (int) (originx + rx * cos(b) * cos(l));
					image->set_value_at(x,y,z, value);
				}
			}
		}
	}
}

void TestTomoImage::process_inplace( EMData* image )
{
	float nx = 240;
	float ny = 240;
	float nz = 60;
	
	// This increment is used to simplified positioning
	// It's an incremental factor that matches the grid size of the paper
	// that I drew this design on before implementing it in code
	float inc = 1.0/22.0;
	
	image->set_size((int)nx,(int)ny,(int)nz);
	
	Region region(nx/2.0,ny/2.,nz/2.,.4*nx,0.4*ny,0.4*nz);
	insert_solid_ellipse(image, region, 0.1);
	insert_hollow_ellipse(image, region, 0.2, 3);	
	
	// Center x, center z, bottom y ellipsoids that grow progessively smaller
	{
		Region region(nx/2.,ny*4.*inc,nz/2.,2.*inc*nx,0.5*inc*ny,1.*inc*nz);
		insert_solid_ellipse(image, region, 0.3);
	}
	
	{
		Region region(nx/2.,ny*5.5*inc,nz/2.,1.5*inc*nx,.5*inc*ny,1.*inc*nz);
		insert_solid_ellipse(image, region, 0.0);
	}
	
	{
		Region region(nx/2.,ny*7.*inc,nz/2.,1.*inc*nx,0.5*inc*ny,1.*inc*nz);
		insert_solid_ellipse(image, region, 0.3);
	}
	
	{
		Region region(nx/2.,ny*8.5*inc,nz/2.,0.75*inc*nx,0.5*inc*ny,1.*inc*nz);
		insert_solid_ellipse(image, region, 0.0);
	}
	
	// Center x, center z, bottom y ellipsoids that grow progessively smaller
	{
		Region region(nx/2.,ny*18.*inc,nz/2.,2.*inc*nx,0.5*inc*ny,1.*inc*nz);
		insert_solid_ellipse(image, region, 0.3);
	}
	
	{
		Region region(nx/2.,ny*16.5*inc,nz/2.,1.5*inc*nx,.5*inc*ny,1.*inc*nz);
		insert_solid_ellipse(image, region, 0.3);
	}
	
	{
		Region region(nx/2.,ny*15.*inc,nz/2.,1.*inc*nx,0.5*inc*ny,1.*inc*nz);
		insert_solid_ellipse(image, region, 0.3);
	}
	
	{
		Region region(nx/2.,ny*13.5*inc,nz/2.,0.75*inc*nx,0.5*inc*ny,1.*inc*nz);
		insert_solid_ellipse(image, region, 0.3);
	}
	
	// Left ellipsoids from the bottom up
	{
		Region region(nx*7.*inc,ny*5.*inc,nz/2.,1.*inc*nx,0.75*inc*ny,1.5*inc*nz);
		insert_solid_ellipse(image, region, 0.25);
	}
	
	{
		Region region(nx*7.*inc,ny*7.*inc,nz/2.,1.5*inc*nx,0.75*inc*ny,1.5*inc*nz);
		insert_solid_ellipse(image, region, 0.25);
	}
	
	{
		Region region(nx*7.*inc,ny*9.*inc,nz/2.,2.*inc*nx,0.75*inc*ny,1.*inc*nz);
		insert_solid_ellipse(image, region, 0.25);
	}
	
	{
		Region region(nx*7.*inc,ny*11.*inc,nz/2.,2.5*inc*nx,0.75*inc*ny,1.*inc*nz);
		insert_solid_ellipse(image, region, 0.25);
	}
	
	{
		Region region(nx*7.*inc,ny*13.*inc,nz/2.,3.*inc*nx,0.75*inc*ny,1.*inc*nz);
		insert_solid_ellipse(image, region, 0.25);
	}
	
	// Right rectangle from the top down
	{
		Region region(nx*15.*inc,ny*17.*inc,nz/2.,1.*inc*nx,1.5*inc*ny,1.5*inc*nz);
		insert_rectangle(image, region, 0.25);
	}
	{
		Region region(nx*15.*inc,ny*15.*inc,nz/2.,1.5*inc*nx,1.5*inc*ny,1.5*inc*nz);
		insert_rectangle(image, region, 0.25);
	}
	{
		Region region(nx*15.*inc,ny*13.*inc,nz/2.,2.*inc*nx,1.5*inc*ny,1.5*inc*nz);
		insert_rectangle(image, region, 0.25);
	}
	{
		Region region(nx*15.*inc,ny*11.*inc,nz/2.,2.5*inc*nx,1.5*inc*ny,1.5*inc*nz);
		insert_rectangle(image, region, 0.25);
	}
	{
		Region region(nx*15.*inc,ny*9.*inc,nz/2.,3.*inc*nx,1.5*inc*ny,1.5*inc*nz);
		insert_rectangle(image, region, 0.25);
	}
	
	// Center rotated rectangle
	{
		Region region(nx/2.,ny/2.,nz/2.,2.*inc*nx,2.5*inc*ny,1.*inc*nz);
		Transform3D t3d(-25.0,0.,0.);
		insert_rectangle(image, region, 0.4, &t3d);
	}
	
	// Rotated ellipsoids
	{
		Region region(nx*6.8*inc,ny*16.*inc,nz/2.,1.5*inc*nx,0.5*inc*ny,0.5*inc*nz);
		Transform3D t3d(45.0,0.,0.);
		insert_solid_ellipse(image, region, 0.2, &t3d);
	}
	{
		Region region(nx*7.2*inc,ny*16.*inc,nz/2.,1.5*inc*nx,0.5*inc*ny,0.5*inc*nz);
		Transform3D t3d(135.0,0.,0.);
		insert_solid_ellipse(image, region, 0.3, &t3d);
	}
	
	// Dense small ellipsoids
	{
		Region region(nx*4.*inc,ny*8.*inc,nz/2.,0.5*inc*nx,0.5*inc*ny,0.5*inc*nz);
		insert_solid_ellipse(image, region, 2.05);
	}
	{
		Region region(nx*8.*inc,ny*18.*inc,nz/2.,0.5*inc*nx,0.5*inc*ny,0.5*inc*nz);
		insert_solid_ellipse(image, region, 2.05);
	}
	{
		Region region(nx*14.*inc,ny*18.2*inc,nz/2.,0.5*inc*nx,0.5*inc*ny,0.5*inc*nz);
		insert_solid_ellipse(image, region, 2.05);
	}
	
	{
		Region region(nx*18.*inc,ny*14.*inc,nz/2.,0.5*inc*nx,0.5*inc*ny,0.5*inc*nz);
		insert_solid_ellipse(image, region, 2.05);
	}
	
	{
		Region region(nx*17.*inc,ny*7.5*inc,nz/2.,0.5*inc*nx,0.5*inc*ny,0.5*inc*nz);
		insert_solid_ellipse(image, region, 2.05);
	}
	
	// Dense small rectangles
	{
		Region region(nx*18.*inc,ny*11.5*inc,nz/2.,1.*inc*nx,1.*inc*ny,1.*inc*nz);
		Transform3D t3d(45.0,0.,0.);
		insert_rectangle(image, region, 1.45, &t3d);
	}
	{
		Region region(nx*3.5*inc,ny*10.5*inc,nz/2.,1.*inc*nx,1.*inc*ny,1.*inc*nz);
		Transform3D t3d(45.0,0.,0.);
		insert_rectangle(image, region, 1.45, &t3d);
	}
	
	// Insert small cluster of spheres
	{
		Region region(nx*14.*inc,ny*7.5*inc,nz/2.,0.5*inc*nx,0.5*inc*ny,0.5*inc*nz);
		insert_solid_ellipse(image, region, .35);
	}
	{
		Region region(nx*15.*inc,ny*7.5*inc,nz/2.,0.25*inc*nx,0.25*inc*ny,0.25*inc*nz);
		insert_solid_ellipse(image, region, .35);
	}
	{
		Region region(nx*13.5*inc,ny*6.5*inc,nz/2.,0.25*inc*nx,0.25*inc*ny,0.25*inc*nz);
		insert_solid_ellipse(image, region, .35);
	}
	{
		Region region(nx*14.5*inc,ny*6.5*inc,nz/2.,0.25*inc*nx,0.25*inc*ny,0.25*inc*nz);
		insert_solid_ellipse(image, region, .35);
	}
	{
		Region region(nx*15.5*inc,ny*6.5*inc,nz/2.,0.25*inc*nx,0.25*inc*ny,0.25*inc*nz);
		insert_solid_ellipse(image, region, .35);
	}
	
	{
		Region region(nx*14.*inc,ny*5.5*inc,nz/2.,0.25*inc*nx,0.25*inc*ny,0.25*inc*nz);
		insert_solid_ellipse(image, region, .35);
	}
	{
		Region region(nx*15.*inc,ny*5.5*inc,nz/2.,0.25*inc*nx,0.25*inc*ny,0.25*inc*nz);
		insert_solid_ellipse(image, region, .35);
	}
	{
		Region region(nx*16.*inc,ny*5.5*inc,nz/2.,0.25*inc*nx,0.25*inc*ny,0.25*inc*nz);
		insert_solid_ellipse(image, region, .35);
	}

	{
		Region region(nx*14.5*inc,ny*4.5*inc,nz/2.,0.25*inc*nx,0.25*inc*ny,0.25*inc*nz);
		insert_solid_ellipse(image, region, .35);
	}
	{
		Region region(nx*15.5*inc,ny*4.5*inc,nz/2.,0.25*inc*nx,0.25*inc*ny,0.25*inc*nz);
		insert_solid_ellipse(image, region, .35);
	}
	
	// Insert feducials around the outside of the "cell"
// 	for ( float i = 0.; i < 3.; i += 1. ) {
// 		for ( float j = 0.; j < 3.; j += 1. ) {
// 			Region region(nx*2.+i*inc,ny*2.+j*inc,nz/2.,0.05*inc*nx,0.05*inc*ny,0.05*inc*nz);
// 			insert_solid_ellipse(image, region, 2.0);
// 		}
// 	}
	
}

void NSigmaClampingProcessor::process_inplace(EMData *image)
{
	float nsigma = params.set_default("nsigma",default_sigma);
	float sigma = image->get_attr("sigma");
	float mean = image->get_attr("mean");
	params.set_default("minval",mean - nsigma*sigma);
	params.set_default("maxval",mean + nsigma*sigma);
				
	ClampingProcessor::process_inplace(image);
}

void HistogramBin::process_inplace(EMData *image)
{
	float min = image->get_attr("minimum");
	float max = image->get_attr("maximum");
	float nbins = params.set_default("nbins",default_bins);
	bool debug = params.set_default("debug",false);
	
	vector<int> debugscores;
	if ( debug ) {
		debugscores = vector<int>((int)nbins, 0);
	}
	
	if ( nbins < 0 ) throw InvalidParameterException("nbins must be greater than 0");
	
	float bin_width = (max-min)/nbins;
	float bin_val_offset = bin_width/2.0;
	
	int size = image->get_xsize()*image->get_ysize()*image->get_zsize();
	float* dat = image->get_data();
	
	for(int i = 0; i < size; ++i ) {
		float val = dat[i];
		val -= min;
		int bin = (int) (val/bin_width);
		
		// This makes the last interval [] and not [)
		if (bin == nbins) bin -= 1;
		
		dat[i] = min + bin*bin_width + bin_val_offset;
		if ( debug ) {
			debugscores[bin]++;
		}
	}
	
	if ( debug ) {
		int i = 0;
		for( vector<int>::const_iterator it = debugscores.begin(); it != debugscores.end(); ++it, ++i)
			cout << "Bin " << i << " has " << *it << " pixels in it" << endl;
	}
	
}
void ConvolutionProcessor::process_inplace(EMData* image)
{
	bool complexflag = false;
	EMData* null = 0;
	EMData* with = params.set_default("with", null);
	if ( with == NULL ) throw InvalidParameterException("Error - the image required for the convolution is null");
	
	if ( image->is_complex() != with->is_complex() ) throw ImageFormatException("Error - one of the images is complex and one is real. Can not proceed");
	
	if ( image->is_complex() ) complexflag = true;
	else {
		image->process_inplace("xform.phaseorigin.tocorner");
		image->do_fft_inplace();
		with->process_inplace("xform.phaseorigin.tocorner");
		with->do_fft_inplace();
	}
	
	int nx = image->get_xsize();
	int ny = image->get_ysize();
	int nz = image->get_zsize();
	
	for (int iz = 0; iz < nz; iz++) {
		for (int iy = 0; iy < ny; iy++) {
			for (int ix = 0; ix < nx/2; ix++) {
				image->cmplx(ix,iy,iz) *= with->cmplx(ix,iy,iz);
			}
		}
	}
	
	if ( complexflag == false ) {
		image->do_ift_inplace();
		image->postift_depad_corner_inplace();
		image->process_inplace("xform.phaseorigin.tocenter");
	}
	
}

void XGradientProcessor::process_inplace( EMData* image )
{
	if (image->is_complex()) throw ImageFormatException("Cannot edge detect a complex image");
	
	EMData* e = new EMData();
	int nx = image->get_xsize();
	int ny = image->get_ysize();
	int nz = image->get_zsize();
	
	if ( nz == 1 && ny == 1 ) {
		if ( nx < 3 ) throw ImageDimensionException("Error - cannot edge detect an image with less than three pixels");
		
		e->set_size(3,1,1);
		e->set_value_at(0,-1);
		e->set_value_at(2, 1);
		
		Region r = Region(-nx/2+1,nx);
		e->clip_inplace(r);
	} else if ( nz == 1 ) {
		if ( nx < 3 || ny < 3 ) throw ImageDimensionException("Error - cannot edge detect an image with less than three pixels");
		e->set_size(3,3,1);
		e->set_value_at(0,0,-1);
		e->set_value_at(0,1,-2);
		e->set_value_at(0,2,-1);
		
		e->set_value_at(2,0,1);
		e->set_value_at(2,1,2);
		e->set_value_at(2,2,1);
		Region r = Region(-nx/2+1,-ny/2+1,nx,ny);
		e->clip_inplace(r);
	} else {
		if ( nx < 3 || ny < 3 || nz < 3) throw ImageDimensionException("Error - cannot edge detect an image with less than three pixels");
		e->set_size(3,3,3);
		e->set_value_at(0,0,0,-1);
		e->set_value_at(0,1,0,-1);
		e->set_value_at(0,2,0,-1);
		e->set_value_at(0,0,1,-1);
		e->set_value_at(0,1,1,-8);
		e->set_value_at(0,2,1,-1);
		e->set_value_at(0,0,2,-1);
		e->set_value_at(0,1,2,-1);
		e->set_value_at(0,2,2,-1);
		
		e->set_value_at(2,0,0,1);
		e->set_value_at(2,1,0,1);
		e->set_value_at(2,2,0,1);
		e->set_value_at(2,0,1,1);
		e->set_value_at(2,1,1,8);
		e->set_value_at(2,2,1,1);
		e->set_value_at(2,0,2,1);
		e->set_value_at(2,1,2,1);
		e->set_value_at(2,2,2,1);
		
		Region r = Region(-nx/2+1,-ny/2+1,-nz/2+1,nx,ny,nz);
		e->clip_inplace(r);
	}
	
	Dict conv_parms;
	conv_parms["with"] = e;
	image->process_inplace("convolution", conv_parms);
	
	delete e;
}

class GaussianFunctoid
{
	public:
		GaussianFunctoid(const float sigma, const float mean = 0.0) : m_mean(mean), m_sigma_squared(sigma*sigma) {}
		~GaussianFunctoid() {}
		
		float operator()(const float distance ) 
		{
			return exp( -(distance-m_mean)*(distance-m_mean)/ (m_sigma_squared ));
		}
	private:
		float m_mean, m_sigma_squared;
};

void TomoTiltAngleWeightProcessor::process_inplace( EMData* image )
{
	bool fim = params.set_default("angle_fim", false);
	float alt;
	if ( fim ) {
		alt = image->get_attr("euler_alt");
	}
	else alt = params.set_default("angle", 0.0);
	
	float cosine = cos(alt*M_PI/180.0f);
	float mult_fac =  1.0f/(cosine);
	image->mult( mult_fac );
}

void TomoTiltEdgeMaskProcessor::process_inplace( EMData* image )
{
	bool biedgemean = params.set_default("biedgemean", false);
	bool edgemean = params.set_default("edgemean", false);
	// You can only do one of these - so if someone specifies them both the code complains loudly
	if (biedgemean && edgemean) throw InvalidParameterException("The edgemean and biedgemean options are mutually exclusive");
	
	bool fim = params.set_default("angle_fim", false);
	float alt;
	if ( fim ) {
		alt = image->get_attr("euler_alt");
	}
	else alt = params.set_default("angle", 0.0);
	
	
	float cosine = cos(alt*M_PI/180.0f);
		
	// Zero the edges
	int nx = image->get_xsize();
	int ny = image->get_ysize();
	int x_clip = static_cast<int>( (float) nx * ( 1.0 - cosine ) / 2.0);
	
	float x1_edge_mean = 0.0;
	float x2_edge_mean = 0.0;

	if ( biedgemean )
	{
		float edge_mean = 0.0;
			
		// Accrue the pixel densities on the side strips
		for ( int i = 0; i < ny; ++i ) {
			edge_mean += image->get_value_at(x_clip, i );
			edge_mean += image->get_value_at(nx - x_clip-1, i );
		}
		// Now make it so the mean is stored
		edge_mean /= 2*ny;
		
		// Now shift pixel values accordingly
		for ( int i = 0; i < ny; ++i ) {
			for ( int j = nx-1; j >= nx - x_clip; --j) {
				image->set_value_at(j,i,edge_mean);
			}
			for ( int j = 0; j < x_clip; ++j) {
				image->set_value_at(j,i,edge_mean);
			}
		}
		x1_edge_mean = edge_mean;
		x2_edge_mean = edge_mean;
	}
	else if (edgemean)
	{
		for ( int i = 0; i < ny; ++i ) {
			x1_edge_mean += image->get_value_at(x_clip, i );
			x2_edge_mean += image->get_value_at(nx - x_clip-1, i );
		}
		x1_edge_mean /= ny;
		x2_edge_mean /= ny;
		
		for ( int i = 0; i < ny; ++i ) {
			for ( int j = 0; j < x_clip; ++j) {
				image->set_value_at(j,i,x1_edge_mean);
			}
			for ( int j = nx-1; j >= nx - x_clip; --j) {
				image->set_value_at(j,i,x2_edge_mean);
			}
		}
	}
	else
	{
		// The edges are just zeroed - 
		Dict zero_dict;
		zero_dict["x0"] = x_clip;
		zero_dict["x1"] = x_clip;
		zero_dict["y0"] = 0;
		zero_dict["y1"] = 0;
		image->process_inplace( "mask.zeroedge2d", zero_dict );
	}

	int gauss_rad = params.set_default("gauss_falloff", 0);
	if ( gauss_rad != 0)
	{
		// If the gaussian falloff distance is greater than x_clip, it will technically
		// go beyond the image boundaries. Thus we clamp gauss_rad so this cannot happen.
		// Therefore, there is potential here for (benevolent) unexpected behavior.
		if ( gauss_rad > x_clip ) gauss_rad = x_clip;
		
		float gauss_sigma = params.set_default("gauss_sigma", 3.0f);
		if ( gauss_sigma < 0 ) throw InvalidParameterException("Error - you must specify a positive, non-zero gauss_sigma");
		float sigma = (float) gauss_rad/gauss_sigma;
		
		GaussianFunctoid gf(sigma);
		
		for ( int i = 0; i < ny; ++i ) {
			
			float left_value = image->get_value_at(x_clip, i );
			float scale1 = left_value-x1_edge_mean;
			
			float right_value = image->get_value_at(nx - x_clip - 1, i );
			float scale2 = right_value-x2_edge_mean;
			
			for ( int j = 1; j < gauss_rad; ++j )
			{
				image->set_value_at(x_clip-j, i, scale1*gf((float)j)+x1_edge_mean );
				image->set_value_at(nx - x_clip + j-1, i, scale2*gf((float)j)+x2_edge_mean);
			}
		}
	}
	
	image->update();
}

void YGradientProcessor::process_inplace( EMData* image )
{
	if (image->is_complex()) throw ImageFormatException("Cannot edge detect a complex image");
	
	EMData* e = new EMData();
	int nx = image->get_xsize();
	int ny = image->get_ysize();
	int nz = image->get_zsize();
	
	if ( nz == 1 && ny == 1 ) {
		throw ImageDimensionException("Error - cannot detect Y edges for an image that that is 1D!");
	} else if ( nz == 1 ) {
		if ( nx < 3 || ny < 3 ) throw ImageDimensionException("Error - cannot edge detect an image with less than three pixels");
		e->set_size(3,3,1);
		e->set_value_at(0,0,-1);
		e->set_value_at(1,0,-2);
		e->set_value_at(2,0,-1);
		
		e->set_value_at(0,2,1);
		e->set_value_at(1,2,2);
		e->set_value_at(2,2,1);
		Region r = Region(-nx/2+1,-ny/2+1,nx,ny);
		e->clip_inplace(r);
	} else {
		if ( nx < 3 || ny < 3 || nz < 3) throw ImageDimensionException("Error - cannot edge detect an image with less than three pixels");
		e->set_size(3,3,3);
		e->set_value_at(0,0,0,-1);
		e->set_value_at(1,0,0,-1);
		e->set_value_at(2,0,0,-1);
		e->set_value_at(0,0,1,-1);
		e->set_value_at(1,0,1,-8);
		e->set_value_at(2,0,1,-1);
		e->set_value_at(0,0,2,-1);
		e->set_value_at(1,0,2,-1);
		e->set_value_at(2,0,2,-1);
		
		e->set_value_at(0,2,0,1);
		e->set_value_at(1,2,0,1);
		e->set_value_at(2,2,0,1);
		e->set_value_at(0,2,1,1);
		e->set_value_at(1,2,1,8);
		e->set_value_at(2,2,1,1);
		e->set_value_at(0,2,2,1);
		e->set_value_at(1,2,2,1);
		e->set_value_at(2,2,2,1);
		
		Region r = Region(-nx/2+1,-ny/2+1,-nz/2+1,nx,ny,nz);
		e->clip_inplace(r);
	}
	
	Dict conv_parms;
	conv_parms["with"] = e;
	image->process_inplace("convolution", conv_parms);
	
	delete e;
}


void ZGradientProcessor::process_inplace( EMData* image )
{
	if (image->is_complex()) throw ImageFormatException("Cannot edge detect a complex image");
	
	EMData* e = new EMData();
	int nx = image->get_xsize();
	int ny = image->get_ysize();
	int nz = image->get_zsize();
	
	if ( nx < 3 || ny < 3 || nz < 3) throw ImageDimensionException("Error - cannot edge detect in the z direction with any dimension being less than three pixels");
	
	e->set_size(3,3,3);
	e->set_value_at(0,0,0,-1);
	e->set_value_at(1,0,0,-1);
	e->set_value_at(2,0,0,-1);
	e->set_value_at(0,1,0,-1);
	e->set_value_at(1,1,0,-8);
	e->set_value_at(2,1,0,-1);
	e->set_value_at(0,2,0,-1);
	e->set_value_at(1,2,0,-1);
	e->set_value_at(2,2,0,-1);
	
	e->set_value_at(0,0,2,1);
	e->set_value_at(1,0,2,1);
	e->set_value_at(2,0,2,1);
	e->set_value_at(0,1,2,1);
	e->set_value_at(1,1,2,8);
	e->set_value_at(2,1,2,1);
	e->set_value_at(0,2,2,1);
	e->set_value_at(1,2,2,1);
	e->set_value_at(2,2,2,1);
	
	Region r = Region(-nx/2+1,-ny/2+1,-nz/2+1,nx,ny,nz);
	e->clip_inplace(r);
	
	Dict conv_parms;
	conv_parms["with"] = e;
	image->process_inplace("convolution", conv_parms);
	
	delete e;
}

void EMAN::dump_processors()
{
	dump_factory < Processor > ();
}

map<string, vector<string> > EMAN::dump_processors_list()
{
	return dump_factory_list < Processor > ();
}

map<string, vector<string> > EMAN::group_processors()
{
	map<string, vector<string> > processor_groups;

	vector <string> processornames = Factory<Processor>::get_list();

	for (size_t i = 0; i < processornames.size(); i++) {
		Processor * f = Factory<Processor>::get(processornames[i]);
		if (dynamic_cast<RealPixelProcessor*>(f) != 0) {
			processor_groups["RealPixelProcessor"].push_back(f->get_name());				
		}
		else if (dynamic_cast<BoxStatProcessor*>(f)  != 0) {
			processor_groups["BoxStatProcessor"].push_back(f->get_name());				
		}
		else if (dynamic_cast<ComplexPixelProcessor*>(f)  != 0) {
			processor_groups["ComplexPixelProcessor"].push_back(f->get_name());				
		}
		else if (dynamic_cast<CoordinateProcessor*>(f)  != 0) {
			processor_groups["CoordinateProcessor"].push_back(f->get_name());				
		}
		else if (dynamic_cast<FourierProcessor*>(f)  != 0) {
			processor_groups["FourierProcessor"].push_back(f->get_name());				
		}
		else if (dynamic_cast<NewFourierProcessor*>(f)  != 0) {
			processor_groups["FourierProcessor"].push_back(f->get_name());				
		}
		else if (dynamic_cast<NormalizeProcessor*>(f)  != 0) {
			processor_groups["NormalizeProcessor"].push_back(f->get_name());				
		}
		else {
			processor_groups["Others"].push_back(f->get_name());				
		}
	}
		
	return processor_groups;
}

/* vim: set ts=4 noet: */
