/**
 * $Id$
 */
#include "cmp.h"
#include "emdata.h"
#include "log.h"
#include "util.h"
#include "ctf.h"

using namespace EMAN;

template <> Factory < Cmp >::Factory()
{
	force_add(&DotCmp::NEW);
	force_add(&VarianceCmp::NEW);
	force_add(&PhaseCmp::NEW);
	force_add(&FRCCmp::NEW);
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
		throw NullPointerException("compared image data");
	}
	
	float *d2 = with->get_data();
	if (!d2) {
		throw NullPointerException("compare-with image data");
	}
}


float DotCmp::cmp(EMData * image, EMData *with) const
{
	ENTERFUNC;
	validate_input_args(image, with);

	float *d1 = image->get_data();
	float *d2 = with->get_data();

	int evenonly = params.set_default("evenonly", 0);

	double result = 0;
	size_t size = image->get_xsize() * image->get_ysize() * image->get_zsize();
	int step = 1;

	if (evenonly) {
		step = 2;
	}

	for (size_t i = 0; i < size; i += step) {
		result += (*d1) * (*d2);
		d1 += step;
		d2 += step;
	}
#if 0
	double square_sum1 = image->get_attr_dict().get("square_sum");
	double square_sum2 = with->get_attr_dict().get("square_sum");

	result = 2 * result / (square_sum1 + square_sum2);
#endif
	
	EXITFUNC;
	return (float) result;
}


float OptVarianceCmp::cmp(EMData * image, EMData *with) const
{
	ENTERFUNC;
	validate_input_args(image, with);

	int matchfilt = params.set_default("matchfilt",0);
	
	float *y_data = with->get_data();
	float *x_data = image->get_data();
	
	EMData *with2=NULL;
	if (matchfilt) {
		with2=with->copy();
//		with2->process("matchfilt",Dict("to",this));
		y_data = with2->get_data();
	}
		

	size_t size = image->get_xsize() * image->get_ysize() * image->get_zsize();
	float m = 0;
	float b = 0;

	Util::calc_least_square_fit(size, x_data, y_data, &m, &b, 1);
	if (m == 0) {
		m = FLT_MIN;
	}
	b = -b / m;
	m = 1.0f / m;
	if (m < 0) {
		b = 0;
		m = 1.0f;
	}

	int keepzero = params.set_default("keepzero", 0);
	double  result = 0;
	int count = 0;


	if (keepzero) {
		for (size_t i = 0; i < size; i++) {
			if (y_data[i] && x_data[i]) {
				result += Util::square((x_data[i] * m) + b - y_data[i]);
				count++;
			}
		}
		result/=count;
	}
	else {
		for (size_t i = 0; i < size; i++) {
			result += Util::square((x_data[i] * m) + b - y_data[i]);
		}
		result = result / size;
	}
	scale = m;
	shift = b;
	
	EXITFUNC;
	
#if 0
	return (1 - result);
#endif
	
	return result;
}

float VarianceCmp::cmp(EMData * image, EMData *with) const
{
	ENTERFUNC;
	validate_input_args(image, with);

	float *y_data = with->get_data();
	float *x_data = image->get_data();
	double result = 0;
	
	size_t size = image->get_xsize() * image->get_ysize() * image->get_zsize();
	
	for (size_t i = 0; i < size; i++) {
		result += Util::square(x_data[i]- y_data[i]);
	}
	result/=size;
	
	EXITFUNC;
	
#if 0
	return (1 - result);
#endif
	
	return result;
}

float PhaseCmp::cmp(EMData * image, EMData *with) const
{
	ENTERFUNC;
	validate_input_args(image, with);

	static float *dfsnr = 0;
	static int nsnr = 0;

	if (image->get_zsize() > 1) {
		throw ImageDimensionException("2D only");
	}

	int nx = image->get_xsize();
	int ny = image->get_ysize();

	int np = (int) ceil(Ctf::CTFOS * sqrt(2.0f) * ny / 2) + 2;

	if (nsnr != np) {
		nsnr = np;
		dfsnr = (float *) realloc(dfsnr, np * sizeof(float));

		float w = Util::square(nx / 4.0f);

		for (int i = 0; i < np; i++) {
			float x2 = Util::square(i / (float) Ctf::CTFOS);
			dfsnr[i] = (1.0f - exp(-x2 / 4.0f)) * exp(-x2 / w);
		}

		Util::save_data(0, 1.0f / Ctf::CTFOS, dfsnr, np, "filt.txt");
	}

	EMData *image_fft = image->do_fft();
	image_fft->ri2ap();
	EMData *with_fft = with->do_fft();
	with_fft->ri2ap();

	float *image_fft_data = image_fft->get_data();
	float *with_fft_data = with_fft->get_data();
	double sum = 0;
	double norm = FLT_MIN;
	int i = 0;

	for (float y = -ny / 2.0f; y < ny / 2.0f; y++) {
		for (int x = 0; x < nx + 2; x += 2) {
			int r = Util::round(hypot(x / 2, y) * Ctf::CTFOS);
			float a = dfsnr[r] * with_fft_data[i];

			sum += Util::angle_sub_2pi(image_fft_data[i + 1], with_fft_data[i + 1]) * a;
			norm += a;
			i += 2;
		}
	}
	EXITFUNC;
	
	if( image_fft )
	{
		delete image_fft;
		image_fft = 0;
	}
	if( with_fft )
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
	
	static vector < float >default_snr;

	if (image->get_zsize() > 1) {
		throw ImageDimensionException("2D only");
	}

	int nx = image->get_xsize();
	int ny = image->get_ysize();

	vector < float >snr = params["snr"];
	vector < float >fsc_array;

	if (snr.size() == 0) {
		int np = (int) ceil(Ctf::CTFOS * sqrt(2.0f) * ny / 2) + 2;

		fsc_array = image->calc_fourier_shell_correlation(with);

		if (default_snr.size() != (unsigned int) np) {
			default_snr = vector < float >(np);
			float w = Util::square(nx / 8.0f);

			for (int i = 0; i < np; i++) {
				float x2 = Util::square(i / (float) Ctf::CTFOS);
				default_snr[i] = (1.0f - exp(-x2 / 4.0f)) * exp(-x2 / w);
			}
		}
	}

	double sum = 0;
	double norm = 0;

	const int n = Ctf::CTFOS;
	for (int i = 0; i < ny / 2; i++) {
		sum += fsc_array[i] * i * default_snr[i * n + n / 2];
		norm += i * default_snr[i * n + n / 2];
	}

	EXITFUNC;
	
	return (float)(sum / norm);
}

void EMAN::dump_cmps()
{
	dump_factory < Cmp > ();
}
