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


float DotCmp::cmp(EMData * em, Transform *) const
{
	EMData *with = params["with"];
	if (!with || !EMUtil::is_same_size(em, with)) {
		return 0;
	}

	float *d1 = em->get_data();
	float *d2 = with->get_data();

	int evenonly = params.set_default("evenonly", 0);

	double result = 0;
	size_t size = em->get_xsize() * em->get_ysize() * em->get_zsize();
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
	double square_sum1 = em->get_attr_dict().get("square_sum");
	double square_sum2 = with->get_attr_dict().get("square_sum");

	result = 2 * result / (square_sum1 + square_sum2);
#endif
	return (float) result;
}

// scale and shift cannot be returned
float VarianceCmp::cmp(EMData * em, Transform *) const
{
	EMData *with = params["with"];
	if (!with || !EMUtil::is_same_size(em, with)) {
		return 0;
	}

	float *x_data = with->get_data();
	float *y_data = em->get_data();

	size_t size = em->get_xsize() * em->get_ysize() * em->get_zsize();
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
	float result = 0;
	int count = 0;

	for (size_t i = 0; i < size; i++) {
		if (y_data[i] && x_data[i]) {
			result += Util::square((x_data[i] * m) + b - y_data[i]);
			count++;
		}
	}

	if (keepzero) {
		result = result / count;
	}
	else {
		result = result / size;
	}
	params["scale"] = m;
	params["shift"] = b;
#if 0
	return (1 - result);
#endif
	return result;
}

float PhaseCmp::cmp(EMData * em, Transform *) const
{
	static float *dfsnr = 0;
	static int nsnr = 0;

	EMData *with = params["with"];
	if (!with || !EMUtil::is_same_size(em, with) || em->get_zsize() > 1) {
		return 0;
	}

	int nx = em->get_xsize();
	int ny = em->get_ysize();

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

	EMData *em_fft = em->do_fft();
	em_fft->ri2ap();
	EMData *with_fft = with->do_fft();
	with_fft->ri2ap();

	float *em_fft_data = em_fft->get_data();
	float *with_fft_data = with_fft->get_data();
	double sum = 0;
	double norm = FLT_MIN;
	int i = 0;

	for (float y = -ny / 2.0f; y < ny / 2.0f; y++) {
		for (int x = 0; x < nx + 2; x += 2) {
			int r = Util::round(hypot(x / 2, y) * Ctf::CTFOS);
			float a = dfsnr[r] * with_fft_data[i];

			sum += Util::angle_sub_2pi(em_fft_data[i + 1], with_fft_data[i + 1]) * a;
			norm += a;
			i += 2;
		}
	}
#if 0
	return (1.0f - sum / norm);
#endif
	return (sum / norm);
}


float FRCCmp::cmp(EMData * em, Transform *) const
{
	static vector < float >default_snr;

	EMData *with = params["with"];
	if (!with || !EMUtil::is_same_size(em, with) || em->get_zsize() > 1) {
		return 0;
	}

	int nx = em->get_xsize();
	int ny = em->get_ysize();

	vector < float >snr = params["snr"].get_farray();
	vector < float >fsc_array;

	if (snr.size() == 0) {
		int np = (int) ceil(Ctf::CTFOS * sqrt(2.0f) * ny / 2) + 2;

		fsc_array = em->calc_fourier_shell_correlation(with);

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

	return (sum / norm);
}

void dump_cmps()
{
	dump_factory < Cmp > ();
}
