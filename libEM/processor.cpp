/**
 * $Id$
 */

#include "emconstants.h"
#include "processor.h"
#include "ctf.h"
#include "xydata.h"
#include "emdata.h"
#include "Assert.h"
#include <time.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

using namespace EMAN;

template <> Factory < Processor >::Factory()
{
	force_add(&LowpassSharpCutoffProcessor::NEW);
	force_add(&HighpassSharpCutoffProcessor::NEW);
	force_add(&LowpassGaussProcessor::NEW);
	force_add(&HighpassGaussProcessor::NEW);

	force_add(&LowpassTanhProcessor::NEW);
	force_add(&HighpassTanhProcessor::NEW);
	force_add(&HighpassButterworthProcessor::NEW);
	
	force_add(&LinearRampProcessor::NEW);
	force_add(&AbsoluateValueProcessor::NEW);
	force_add(&BooleanProcessor::NEW);
	force_add(&ValueSquaredProcessor::NEW);
	force_add(&ValueSqrtProcessor::NEW);

	force_add(&ToZeroProcessor::NEW);
	force_add(&BinarizeProcessor::NEW);
	force_add(&CollapseProcessor::NEW);
	force_add(&LinearXformProcessor::NEW);

	force_add(&ExpProcessor::NEW);
	force_add(&RangeThresholdProcessor::NEW);
	force_add(&SigmaProcessor::NEW);
	force_add(&LogProcessor::NEW);

	force_add(&MaskSharpProcessor::NEW);
	force_add(&MaskEdgeMeanProcessor::NEW);
	force_add(&MaskNoiseProcessor::NEW);
	force_add(&MaskGaussProcessor::NEW);
	force_add(&MaskGaussInvProcessor::NEW);
	
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
	
	force_add(&NormalizeStdProcessor::NEW);
	force_add(&NormalizeUnitProcessor::NEW);
	force_add(&NormalizeMaskProcessor::NEW);
	force_add(&NormalizeEdgeMeanProcessor::NEW);
	force_add(&NormalizeCircleMeanProcessor::NEW);
	force_add(&NormalizeLREdgeMeanProcessor::NEW);
	force_add(&NormalizeMaxMinProcessor::NEW);
	force_add(&NormalizeRowProcessor::NEW);
	
	force_add(&NormalizeToStdProcessor::NEW);
	force_add(&NormalizeToFileProcessor::NEW);
	force_add(&NormalizeToLeastSquareProcessor::NEW);

	force_add(&RadialAverageProcessor::NEW);
	force_add(&RadialSubstractProcessor::NEW);
	force_add(&FlipProcessor::NEW);

	force_add(&AddNoiseProcessor::NEW);
	force_add(&AddRandomNoiseProcessor::NEW);

	force_add(&Phase180Processor::NEW);
	force_add(&FourierOriginShiftProcessor::NEW);
	force_add(&AutoMask2DProcessor::NEW);
	force_add(&AutoMask3DProcessor::NEW);
	force_add(&AddMaskShellProcessor::NEW);

	force_add(&ToMassCenterProcessor::NEW);
	force_add(&ACFCenterProcessor::NEW);
	force_add(&SNRProcessor::NEW);

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
	force_add(&TestImageSquarecube::NEW);
	force_add(&TestImageCirclesphere::NEW);
	force_add(&TestImageNoiseUniformRand::NEW);
	force_add(&TestImageNoiseGauss::NEW);
	force_add(&TestImageScurve::NEW);

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

}


void ImageProcessor::process(EMData * image)
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

		ift->done_data();

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

	image->done_data();
}

#define FFTRADIALOVERSAMPLE 4
void FourierProcessor::process(EMData * image)
{
	if (!image) {
		LOGWARN("NULL Image");
		return;
	}

	int array_size = FFTRADIALOVERSAMPLE * image->get_ysize();
	float step=0.5/array_size;
	
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

		ift->done_data();

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

	image->done_data();
}

void LowpassSharpCutoffProcessor::create_radial_func(vector < float >&radial_mask) const
{
	Assert(radial_mask.size() > 0);
	float x = 0 , step = 0.5/radial_mask.size();
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
	float x = 0 , step = 0.5/radial_mask.size();
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
	float x = 0 , step = 0.5/radial_mask.size();
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
	float x = 0 , step = 0.5/radial_mask.size();
	for (size_t i = 0; i < radial_mask.size(); i++) {
		radial_mask[i] = 1.0f - exp(-x * x / (highpass * highpass));
		x += step;
	}
}


void LowpassTanhProcessor::create_radial_func(vector < float >&radial_mask) const
{
	Assert(radial_mask.size() > 0);
	float x = 0 , step = 0.5/radial_mask.size();
	for (size_t i = 0; i < radial_mask.size(); i++) {
		radial_mask[i] = tanh(lowpass - x) / 2.0f + 0.5f;
		x += step;
	}
}

void HighpassTanhProcessor::create_radial_func(vector < float >&radial_mask) const
{
	Assert(radial_mask.size() > 0);
	float x = 0 , step = 0.5/radial_mask.size();
	for (size_t i = 0; i < radial_mask.size(); i++) {
		radial_mask[i] = tanh(x - highpass) / 2.0f + 0.5f;
		x += step;
	}
}


void HighpassButterworthProcessor::create_radial_func(vector < float >&radial_mask) const
{
	Assert(radial_mask.size() > 0);
	float x = 0 , step = 0.5/radial_mask.size();
	for (size_t i = 0; i < radial_mask.size(); i++) {
		float t = highpass / 1.5f / (x + 0.001f);
		radial_mask[i] = 1.0f / (1.0f + t * t);
		x += step;
	}
}


void LinearRampProcessor::create_radial_func(vector < float >&radial_mask) const
{
	Assert(radial_mask.size() > 0);
	float x = 0 , step = 0.5/radial_mask.size();
	float size=radial_mask.size();
	for (size_t i = 0; i < size; i++) {
		radial_mask[i] = intercept + ((slope - intercept) * i) / (size - 1.0f);
		x += step;
	}
}


void RealPixelProcessor::process(EMData * image)
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
	image->done_data();
}

void CoordinateProcessor::process(EMData * image)
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

void ComplexPixelProcessor::process(EMData * image)
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

	image->done_data();
	image->ap2ri();
}



void AreaProcessor::process(EMData * image)
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

	image->done_data();

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

void BoxStatProcessor::process(EMData * image)
{
	if (!image) {
		LOGWARN("NULL Image");
		return;
	}

	int nx = image->get_xsize();
	int ny = image->get_ysize();
	int nz = image->get_zsize();

	int n = 1;
	int areasize = 2 * n + 1;

	int matrix_size = areasize * areasize;
	if (nz > 1) {
		matrix_size *= areasize;
	}

	float *array = new float[matrix_size];
//	image->process("eman1.normalize");

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

	image->done_data();

	if( data2 )
	{
		delete[]data2;
		data2 = 0;
	}
}

void DiffBlockProcessor::process(EMData * image)
{
	if (!image) {
		LOGWARN("NULL Image");
		return;
	}

	int nz = image->get_zsize();

	if (nz > 1) {
		LOGERR("%s Processor doesn't support 3D", get_name().c_str());
		return;
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

	image->done_data();
}


void CutoffBlockProcessor::process(EMData * image)
{
	if (!image) {
		LOGWARN("NULL Image");
		return;
	}
	int nz = image->get_zsize();

	if (nz > 1) {
		LOGERR("%s Processor doesn't support 3D", get_name().c_str());
		return;
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

	image->done_data();
}


void GradientRemoverProcessor::process(EMData * image)
{
	if (!image) {
		LOGWARN("NULL Image");
		return;
	}

	int nz = image->get_zsize();
	if (nz > 1) {
		LOGERR("%s Processor doesn't support 3D model", get_name().c_str());
		return;
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

	image->done_data();
}



void VerticalStripeProcessor::process(EMData * image)
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

	image->done_data();
}

void RealToFFTProcessor::process(EMData *image)
{
// Note : 2D only!
if (!image || image->is_complex() || image->get_zsize()>1) return;

EMData *ff=image->do_fft();
ff->ri2ap();

int nx=image->get_xsize();
int ny=image->get_ysize();

int x,y;

for (y=0; y<ny; y++) image->set_value_at(0,y,0);

for (x=1; x<nx/2; x++) {
	for (y=0; y<ny; y++) {
		float y2;
		if (y<ny/2) y2=y+ny/2;
		else y2=y-ny/2;
		image->set_value_at(x,y,ff->get_value_at(nx-x*2,static_cast<int>(ny-y2)));
	}
}

for (x=nx/2; x<nx; x++) {
	for (y=0; y<ny; y++) {
		float y2;
		if (y<ny/2) y2=y+ny/2;
		else y2=y-ny/2;
		image->set_value_at(x,y,ff->get_value_at(x*2-nx,static_cast<int>(y2)));
	}
}

image->update();
if( ff )
{
	delete ff;
	ff = 0;
}
}

void SigmaZeroEdgeProcessor::process(EMData * image)
{
	if (!image) {
		LOGWARN("NULL Image");
		return;
	}

	if (image->get_zsize() > 1) {
		LOGERR("%s Processor doesn't support 3D model", get_name().c_str());
		return;
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


	image->done_data();
}



void BeamstopProcessor::process(EMData * image)
{
	if (!image) {
		LOGWARN("NULL Image");
		return;
	}
	if (image->get_zsize() > 1) {
		LOGERR("BeamstopProcessor doesn't support 3D model");
		return;
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

	image->done_data();
}



void MeanZeroEdgeProcessor::process(EMData * image)
{
	if (!image) {
		LOGWARN("NULL Image");
		return;
	}
	if (image->get_zsize() > 1) {
		LOGERR("MeanZeroEdgeProcessor doesn't support 3D model");
		return;
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

	image->done_data();
}



void AverageXProcessor::process(EMData * image)
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

	image->done_data();
}


void ZeroEdgeRowProcessor::process(EMData * image)
{
	if (!image) {
		LOGWARN("NULL Image");
		return;
	}

	if (image->get_zsize() > 1) {
		LOGERR("ZeroEdgeRowProcessor is not supported in 3D models");
		return;
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
	image->done_data();
}

void ZeroEdgePlaneProcessor::process(EMData * image)
{
	if (!image) {
		LOGWARN("NULL Image");
		return;
	}

	if (image->get_zsize() <= 1) {
		LOGERR("Use ZeroEdgeRowProcessor for 2D images");
		return;
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
	
	image->done_data();
}


float NormalizeProcessor::calc_sigma(EMData * image) const
{
	return image->get_attr("sigma");
}

void NormalizeProcessor::process(EMData * image)
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

	image->done_data();
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
	float sigma = 1;
	if ((int) params["no_sigma"] == 0) {
		sigma = image->get_attr("sigma");
	}
	return sigma;
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
		return 0;
	}

	float *data = image->get_data();
	float *mask_data = mask->get_data();
	size_t size = image->get_xsize() * image->get_ysize() * image->get_zsize();
	double sum = 0;
	size_t n_norm = 0;

	for (size_t i = 0; i < size; i++) {
		if (mask_data[i] && data[i]) {
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

void NormalizeRowProcessor::process(EMData * image)
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

	image->done_data();
}

float NormalizeStdProcessor::calc_mean(EMData * image) const
{
	if (!image) {
		LOGWARN("NULL Image");
		return 0;
	}
	return image->get_attr("mean");
}

void NormalizeToStdProcessor::process(EMData * image)
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

	float *this_data = image->get_data();
	float *noisy_data = noisy->get_data();
	float m = 0;
	float b = 0;

	if (!invert) {
		Util::calc_least_square_fit(size, this_data, noisy_data, &m, &b, 1);
	}

	if (invert || m < 0) {
		Util::calc_least_square_fit(size, noisy_data, this_data, &m, &b, 1);

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
		image->done_data();
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


void NormalizeToLeastSquareProcessor::process(EMData * image)
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

	image->done_data();
	image->update();

	params["scale"] = scale;
	params["shift"] = shift;
}



void BilateralProcessor::process(EMData * image)
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

	image->done_data();
	
	if( mask )
	{
		free(mask);
		mask = 0;
	}

	if( old_img )
	{
		free(old_img);
		old_img = 0;
	}

	image->update();
}

void RadialAverageProcessor::process(EMData * image)
{
	if (!image || image->is_complex()) {
		LOGWARN("only works on real image. do nothing.");
		return;
	}

	float *rdata = image->get_data();
	int nx = image->get_xsize();
	int ny = image->get_ysize();

	vector < float >dist = image->calc_radial_dist(nx / 2, 0, 1);

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

	image->done_data();
	image->update();
}



void RadialSubstractProcessor::process(EMData * image)
{
	if (!image || image->is_complex()) {
		LOGWARN("only works on real image. do nothing.");
		return;
	}

	float *rdata = image->get_data();
	int nx = image->get_xsize();
	int ny = image->get_ysize();

	vector < float >dist = image->calc_radial_dist(nx / 2, 0, 1);

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

	image->done_data();
	image->update();
}



void FlipProcessor::process(EMData * image)
{
	if (!image) {
		LOGWARN("NULL Image");
		return;
	}

	const char *axis = params["axis"];

	float *d = image->get_data();
	int nx = image->get_xsize();
	int ny = image->get_ysize();
	int nz = image->get_zsize();

	int nxy = nx * ny;

	if (axis == "x") {			// horizonal flip
		for (int z = 0; z < nz; z++) {
			for (int y = 0; y < ny; y++) {
				int i = y * nx + z * nxy;
				int j = i + nx;

				for (int x = 1; x < nx / 2; x++) {
					float t = d[i + x];
					d[i + x] = d[j - x];
					d[j - x] = t;
				}
			}
		}
	}
	else if (axis == "y") {		// vertical flip
		if (nz == 1) {
			for (int i = 1; i < ny / 2; i++) {
				for (int j = 0; j < nx; j++) {
					int l1 = i * nx + j;
					int l2 = (ny - i) * nx + j;
					float t = d[l1];
					d[l1] = d[l2];
					d[l2] = t;
				}
			}
			for (int j = 0; j < nx; j++) {
				d[j] = d[nx + j];
			}
		}
		else {
			for (int i = 1; i < nz / 2; i++) {
				for (int j = 0; j < ny; j++) {
					for (int k = 0; k < nx; k++) {
						int l1 = i * nx * ny + j * nx + k;
						int l2 = (nz - i) * nxy + j * nx + k;
						float t = d[l1];
						d[l1] = d[l2];
						d[l2] = t;
					}
				}
			}

			for (int j = 0; j < ny; j++) {
				for (int k = 0; k < nx; k++) {
					d[j * nx + k] = d[nxy + j * nx + k];
				}
			}
		}
	}
	else if (axis == "z") {
		LOGWARN("flip around z axis is not implemented yet. do nothing");
	}

	image->done_data();
}



void AddNoiseProcessor::process(EMData * image)
{
	if (!image) {
		LOGWARN("NULL Image");
		return;
	}

	float addnoise = params["noise"];
	addnoise *= get_sigma(image);
	float *dat = image->get_data();
	size_t size = image->get_xsize() * image->get_ysize() * image->get_zsize();

	for (size_t j = 0; j < size; j++) {
		dat[j] += Util::get_gauss_rand(addnoise, addnoise / 2);
	}

	image->done_data();
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

void FourierOriginShiftProcessor::process(EMData * image)
{

// To quote Pawel "funky reordering"
	float *d = image->get_data();
	int nx2=image->get_xsize();
	int ny=image->get_ysize();
	int nz=image->get_zsize();
	
	
	if (nz == 1) {
		int l = ny / 2 * nx2;

		for (int i = 0; i < ny / 2; i++) {
			int inx2 = i * nx2;
			for (int j = 0; j < nx2; j++) {
				int k = j + inx2;
				float f = d[k];
				d[k] = d[k + l];
				d[k + l] = f;
			}
		}
	}
	else if (ny != 1) {
		char *t = new char[nx2 * sizeof(float)];

		int k = nx2 * ny * (nz + 1) / 2;
		int l = nx2 * ny * (nz - 1) / 2;
		size_t jj = nx2 * sizeof(float);
		int ii = 0;

		for (int j = 0; j < nz / 2; j++) {
			for (int i = 0; i < ny; i++) {
				memcpy(t, d + ii, jj);

				if (i < ny / 2) {
					memcpy(d + ii, d + ii + k, jj);
					memcpy(d + ii + k, t, jj);
				}
				else {
					memcpy(d + ii, d + ii + l, jj);
					memcpy(d + ii + l, t, jj);
				}
				ii += nx2;
			}
		}
		if( t )
		{
			delete[]t;
			t = 0;
		}
	}
}

void Phase180Processor::process(EMData * image)
{
	if (!image) {
		LOGWARN("NULL Image");
		return;
	}

	int nx = image->get_xsize();
	int ny = image->get_ysize();
	int nz = image->get_zsize();

	int nxy = nx * ny;
	
	
	float *rdata = image->get_data();

	if (image->is_complex()) {
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
	else {
		int nxyz = nxy * nz;
		int half_nx = nx / 2;
		int half_nxy = nxy / 2;
		int half_nxyz = nxyz / 2;
		
		int half_nx_nxy = half_nx + half_nxy;
		int half_nx_nxyz = half_nx + half_nxyz;
		int half_nxy_nxyz = half_nxy + half_nxyz;
		int half_nx_nxy_nxyz = half_nx + half_nxy + half_nxyz;
		
		int half_ny = ny / 2;
		int half_nz = nz / 2;

		
		if (nz == 1) {
			for (int j = 0; j < half_ny; j++) {
				int j2 = j * nx;
				int j3 = j2 + half_nx_nxy;
				
				for (int i = 0; i < half_nx; i++) {
					int l = i + j2;
					int m = i + j3;
					
					float t = rdata[l];
					rdata[l] = rdata[m];
					rdata[m] = t;
				}
			}
			for (int j = 0; j < half_ny; j++) {
				int j1 = j * nx;
				int j2 = j1 + half_nx;
				int j3 = j1 + half_nxy;
				
				for (int i = 0; i < half_nx; i++) {
					int l = i + j2;
					int m = i + j3;
					
					float t = rdata[l];
					rdata[l] = rdata[m];
					rdata[m] = t;
				}
			}
		}
		else {
			for (int k = 0; k < half_nz; k++) {
				int k2 = k * nxy;
				
				for (int j = 0; j < half_ny; j++) {
					int j2 = j * nx + k2;
					int j3 = j2 + half_nx_nxy_nxyz;
					
					for (int i = 0; i < half_nx; i++) {
						int l = i + j2;
						int m = i + j3;
						
						float t = rdata[l];
						rdata[l] = rdata[m];
						rdata[m] = t;
					}
				}
			}

			for (int k = 0; k < half_nz; k++) {
				int k2 = k * nxy;
				
				for (int j = 0; j < half_ny; j++) {
					int j2 = j * nx + k2 + half_nx;
					int j3 = j * nx + k2 + half_nxy_nxyz;
					
					for (int i = 0; i < half_nx; i++) {
						int l = i + j2;
						int m = i + j3;
						
						float t = rdata[l];
						rdata[l] = rdata[m];
						rdata[m] = t;
					}
				}
			}

			for (int k = 0; k < half_nz; k++) {
				int k2 = k * nxy;
				
				for (int j = 0; j < half_ny; j++) {
					int j2 = j * nx + k2 + half_nx_nxy;
					int j3 = j * nx + k2 + half_nxyz;
					
					for (int i = 0; i < half_nx; i++) {
						int l = i + j2;
						int m = i + j3;
						
						float t = rdata[l];
						rdata[l] = rdata[m];
						rdata[m] = t;
					}
				}
			}

			for (int k = 0; k < half_nz; k++) {
				int k2 = k * nxy;
				
				for (int j = 0; j < half_ny; j++) {
					int j2 = j * nx + k2 + half_nxy;
					int j3 = j * nx + k2 + half_nx_nxyz;
					
					for (int i = 0; i < half_nx; i++) {
						int l = i + j2;
						int m = i + j3;
						
						float t = rdata[l];
						rdata[l] = rdata[m];
						rdata[m] = t;
					}
				}
			}

		}
	}

	image->done_data();
}

void AutoMask2DProcessor::process(EMData * image)
{
	if (!image) {
		LOGWARN("NULL Image");
		return;
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

	d->process("eman1.filter.lowpass.gaussian", Dict("lowpass", (filter * ny / 2)));
	d->process("eman1.filter.highpass.gaussian", Dict("highpass", 0));

	d->process("eman1.normalize");

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

	d->done_data();

	image->mult(*d);
	if( d )
	{
		delete d;
		d = 0;
	}
}



void AddRandomNoiseProcessor::process(EMData * image)
{
	if (!image) {
		LOGWARN("NULL Image");
		return;
	}
	
	if (!image->is_complex()) {
		LOGERR("AddRandomNoise Processor only works for complex image");
		return;
	}

	int n = params["n"];
	float x0 = params["x0"];
	float dx = params["dx"];
	vector < float >y = params["y"];

	int interpolation = 1;
	if (params.has_key("interpolation")) {
		interpolation = params["interpolation"];
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

	image->done_data();
}






void AddMaskShellProcessor::process(EMData * image)
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

	image->done_data();
}

void ToMassCenterProcessor::process(EMData * image)
{
	if (!image) {
		LOGWARN("NULL Image");
		return;
	}

	int int_shift_only = params["int_shift_only"];
	image->process("eman1.normalize");

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

	image->done_data();

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

void ACFCenterProcessor::process(EMData * image)
{
	if (!image) {
		LOGWARN("NULL Image");
		return;
	}

	int is3d = params["is3d"];

	Dict params1;
	params1["intonly"] = 1;
	if (!is3d) {
		params1["maxshift"] = image->get_xsize() / 4;
		image->align("translational", 0, params1);
	}
	else {
		image->align("translational3d", 0, params1);
	}
}


void SNRProcessor::process(EMData * image)
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

	image->process("eman1.normalize.circlemean");

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

void FileFourierProcessor::process(EMData * image)
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

void LocalNormProcessor::process(EMData * image)
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
		printf("Norm filter radius=%1.1f\n", radius);
	}

	EMData *blur = image->copy();
	EMData *maskblur = image->copy();

	maskblur->process("eman1.threshold.binary", Dict("value", threshold));
	maskblur->process("eman1.filter.lowpass.gaussian", Dict("lowpass", radius));
	maskblur->process("eman1.filter.highpass.tanh", Dict("highpass", -10.0f));
	maskblur->process("eman1.threshold.belowtozero", Dict("minval", 0.001f));
	maskblur->process("eman1.threshold.belowtozero", Dict("minval", 0.001f));


	blur->process("eman1.threshold.belowtozero", Dict("minval", threshold));
	blur->process("eman1.filter.lowpass.gaussian", Dict("lowpass", radius));
	blur->process("eman1.filter.highpass.tanh", Dict("highpass", -10.0f));

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


void SymSearchProcessor::process(EMData * image)
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


void IndexMaskFileProcessor::process(EMData * image)
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
		msk->process("eman1.threshold.binaryrange", Dict("low", 0.5f, "high", 1.5f));
	}

	image->mult(*msk);
	if( msk )
	{
		delete msk;
		msk = 0;
	}
}


void CoordinateMaskFileProcessor::process(EMData * image)
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

	image->done_data();
	image->update();
	msk->done_data();
	if( msk )
	{
		delete msk;
		msk = 0;
	}
}


void SetSFProcessor::process(EMData * image)
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

	vector < float >curve = dataf->calc_radial_dist(nx, 0, 0.5f);
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




void SmartMaskProcessor::process(EMData * image)
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
	image->done_data();

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

	image->done_data();
}

void AutoMask3DProcessor::search_nearby(float *dat, float *dat2,
									 int nx, int ny, int nz, float threshold)
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

void AutoMask3DProcessor::process(EMData * image)
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

	image->done_data();
	amask->done_data();

	image->mult(*amask);
	amask->write_image("mask.mrc", 0, EMUtil::IMAGE_MRC);
	if( amask )
	{
		delete amask;
		amask = 0;
	}
}


void AutoMask3D2Processor::process(EMData * image)
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
				if (sqrt((float) (k * k + j * j + i * i)) > radius || dat[l] < threshold) {
					continue;
				}
				dat2[l] = 1.0f;
				l++;
			}
		}
	}

	AutoMask3DProcessor::search_nearby(dat, dat2, nx, ny, nz, threshold);

	amask->done_data();
	float val1 = fabs((float)nshells);
	float val2 = val1 > 2.0f ? 2.0f : 0.0f;

	amask->process("eman1.mask.addshells.gauss", Dict("val1", val1, "val2", val2));

	dat2 = amask->get_data();

	if (nshells < 0) {
		AutoMask3DProcessor::fill_nearby(dat2, nx, ny, nz);
	}

	image->done_data();
	amask->done_data();

	EMData *norm = amask->copy();
	norm->process("eman1.threshold.binary");

	EMData *norm2 = norm->copy();
	norm->process("eman1.mask.addshells", Dict("nshells", 1));
	norm->sub(*norm2);
	if( norm2 )
	{
		delete norm2;
		norm2 = 0;
	}

	norm->process("eman1.mask.addshells", Dict("nshells", 2));
	image->process("eman1.normalize.mask", Dict("mask", norm, "no_sigma", 0));
	
	if( norm )
	{
		delete norm;
		norm = 0;
	}

	image->mult(*amask);
	amask->write_image("mask.mrc", 0, EMUtil::IMAGE_MRC);
}

void IterBinMaskProcessor::process(EMData * image)
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

	image->done_data();
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

void TestImageGaussian::process(EMData * image)
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
	
	image->done_data();
	image->update();
}

void TestImageScurve::process(EMData * image)
{
	preprocess(image);
	
	int nx = image->get_xsize();
	int ny = image->get_ysize();
	image->to_zero();
	MArray2D imdat = image->get_2dview();
	
	for (int i=0; i<100; i++) {
		int x=static_cast<int>( nx/2+nx/6.0*sin(i*2.0*3.14159/100.0) );
		int y=ny/4+i*ny/200;
		for (int xx=x-nx/10; xx<x+nx/10; xx++) {
			for (int yy=y-ny/10; yy<y+ny/10; yy++) {
				imdat[xx][yy]+=exp(-pow(hypot(xx-x,yy-y)*30.0/nx,2.0))*(sin((xx-x)*(yy-y))+.5);
			}
		}
	}
	
	image->update();
}

void TestImagePureGaussian::process(EMData * image)
{
	preprocess(image);
	
	float sigma = params["sigma"];
	int nx = image->get_xsize();
	int ny = image->get_ysize();
	int nz = image->get_zsize();
	int xc = nx/2, yc = ny/2, zc = nz/2;
	float twosig2 = 2*sigma*sigma;
	int d = image->get_ndim();
	float norm = pow(twosig2*pi,-float(d)/2);
	MArray3D imdat = image->get_3dview();
	for (int iz=0; iz < nz; iz++) {
		float z = iz - zc;
		for (int iy=0; iy < ny; iy++) {
			float y = iy - yc;
			for (int ix=0; ix < nx; ix++) {
				float x = ix - xc;
				float r2 = x*x + y*y + z*z;
				float val = norm*exp(-r2/twosig2);
				imdat[ix][iy][iz] = val;
			}
		}
	}
	image->done_data();
	image->update();
}

void TestImageSinewave::process(EMData * image)
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
				*dat = fabs( sin( r * (1.0f/wave_length) + phase) );
			}
		}
	}
	
	image->done_data();
	image->update();	
}

void TestImageSquarecube::process(EMData * image)
{
	preprocess(image);
	
	float edge_length = params["edge_length"];
	string axis = (const char*)params["axis"];
	float odd_edge = params["odd_edge"];
	string fill = (const char*)params["fill"];
	
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
					if( fill == "no" ) {
						*dat = 0;
					} 
					else {
						*dat = 1;
					}
				}
				else {
					if( fill == "no" ) {
						*dat = 1;
					}
					else {
						*dat = 0;
					}
				}
			}
		}
	}
	
	image->done_data();
	image->update();
}

void TestImageCirclesphere::process(EMData * image)
{
	preprocess(image);
	
	float radius = params["radius"];
	string axis = (const char*)params["axis"];
	float c = params["c"];
	string fill = (const char*)params["fill"];
	
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
					if( fill == "no" ) {
						*dat = 0;
					} 
					else {
						*dat = 1;
					}
				}
				else {
					if( fill == "no" ) {
						*dat = 1;
					}
					else {
						*dat = 0;
					}
				}
			}
		}
	}
	
	image->done_data();
	image->update();
}

void TestImageNoiseUniformRand::process(EMData * image)
{
	preprocess(image);
	
	srand(time(0)); //generate a seed by current time
	float *dat = image->get_data();
	for (int i=0; i<nx*ny*nz; i++) {
		dat[i] = (float)rand()/(float)RAND_MAX;
	}
	
	image->done_data();
	image->update();
}

void TestImageNoiseGauss::process(EMData * image)
{
	preprocess(image);
	
	float sigma = params["noise_level"];
	if( sigma == 0 ) {
		sigma = 0.5f;
	}
	
	//gsl_rng_random_glibc2 is the glibc version random number generator
	const gsl_rng * gen = gsl_rng_alloc(gsl_rng_random_glibc2);
	float *dat = image->get_data();
	for (int i=0; i<nx*ny*nz; i++) {
		dat[i] = (float)gsl_ran_gaussian( gen, sigma );
	}
	
	image->done_data();
	image->update();
}

void RampProcessor::process(EMData * image)
{
	if (!image) {
		return;
	}

	int nz = image->get_zsize();
	if (nz > 1) {
		LOGERR("%s Processor doesn't support 3D model", get_name().c_str());
		return;
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
				row[k] -= qy;
				qy += b1;
			}
			d += b2;
		}
	} // image not altered if c is zero

	image->done_data();
}


int EMAN::multi_processors(EMData * image, vector < string > processornames)
{
	Assert(image != 0);
	Assert(processornames.size() > 0);
	
	for (size_t i = 0; i < processornames.size(); i++) {
		image->process(processornames[i]);
	}
	return 0;
}


void EMAN::dump_processors()
{
	dump_factory < Processor > ();
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
