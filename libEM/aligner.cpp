/**
 * $Id$
 */


#include "aligner.h"
#include "log.h"
#include "emdata.h"
#include "filter.h"
#include "util.h"
#include <float.h>
#include <gsl/gsl_multimin.h>

using namespace EMAN;

template <> Factory < Aligner >::Factory()
{
	force_add(&TranslationalAligner::NEW);
	force_add(&Translational3DAligner::NEW);
	force_add(&RotationalAligner::NEW);
	force_add(&RotatePrecenterAligner::NEW);
	force_add(&RotateCHAligner::NEW);
	force_add(&RotateTranslateAligner::NEW);
	force_add(&RotateTranslateBestAligner::NEW);
	force_add(&RotateTranslateRadonAligner::NEW);
	force_add(&RotateFlipAligner::NEW);
	force_add(&RotateTranslateFlipAligner::NEW);
	force_add(&RTFSlowAligner::NEW);
	force_add(&RTFSlowestAligner::NEW);
	force_add(&RTFBestAligner::NEW);
	force_add(&RTFRadonAligner::NEW);
	force_add(&RefineAligner::NEW);
}


EMData *TranslationalAligner::align(EMData * this_img, const string&) const
{
	if (!this_img) {
		return 0;
	}

	int nz = this_img->get_zsize();
	if (nz > 1) {
		LOGERR("%s doesn't support 3D alignment", get_name().c_str());
		return 0;
	}

	EMData *to = params["to"];
	if (to && !EMUtil::is_same_size(this_img, to)) {
		LOGERR("%s: images must be the same size", get_name().c_str());
		return 0;
	}

	int useparent = params["useparent"];

	EMData *cf = 0;
	EMData *parent = this_img->get_parent();

	if (useparent && parent != 0) {
		cf = parent->calc_ccf(to, true);
	}
	else {
		cf = this_img->calc_ccf(to, true);
	}

	int nx = this_img->get_xsize();
	int ny = this_img->get_ysize();
	int maxshift = params["maxshift"];

	if (maxshift <= 0) {
		maxshift = ny / 8;
	}
	if (maxshift > nx / 2 - 1) {
		maxshift = nx / 2 - 1;
	}

	if (maxshift > ny / 2 - 1) {
		maxshift = ny / 2 - 1;
	}


	float *cf_data = cf->get_data();

	float neg = (float)cf->get_attr("mean") - (float)cf->get_attr("minimum");
	float pos = (float)cf->get_attr("maximum") - (float)cf->get_attr("mean");

	int flag = 1;
	if (neg > pos) {
		flag = -1;
	}

	int peak_x = nx / 2;
	int peak_y = ny / 2;

	float max_value = -FLT_MAX;
	float min_value = FLT_MAX;

	for (int j = ny / 2 - maxshift; j < ny / 2 + maxshift; j++) {
		for (int i = nx / 2 - maxshift; i < nx / 2 + maxshift; i++) {
			int l = i + j * nx;

			if (cf_data[l] * flag > max_value) {
				max_value = cf_data[l] * flag;
				peak_x = i;
				peak_y = j;
			}
			if (cf_data[l] < min_value) {
				min_value = cf_data[l];
			}
		}
	}

	Vec3f pre_trans = this_img->get_translation();
	Vec3f cur_trans = Vec3f ((float)(nx / 2 - peak_x), (float)(ny / 2 - peak_y), 0);

	Vec3f result;
	if (useparent && parent) {
		result = cur_trans - pre_trans;
	}
	else {
		result = cur_trans;
	}

	if (!to) {
		cur_trans /= 2.0f;
	}

	int intonly = params["intonly"];

	if (intonly) {
		cur_trans[0] = floor(cur_trans[0] + 0.5f);
		cur_trans[1] = floor(cur_trans[1] + 0.5f);
	}

	this_img->translate(cur_trans);

	float score = (float)hypot(result[0], result[1]);
	cf->set_attr("align_score", max_value);
	cf->set_attr("translational.dx",result[0]); 
	cf->set_attr("translational.dy",result[1]); 
	cf->set_attr("translational.dz",result[2]); 
	cf->done_data();

	return cf;
}



EMData *Translational3DAligner::align(EMData * this_img, const string&) const
{
	if (!this_img) {
		return 0;
	}

	EMData *to = params["to"];
	int useparent = params.set_default("useparent", 0);
	params.set_default("intonly", 0);

	if (to && !EMUtil::is_same_size(this_img, to)) {
		LOGERR("%s: images must be same size", get_name().c_str());
		return 0;
	}

	if (!to) {
		LOGWARN("%s: ACF", get_name().c_str());
	}

	EMData *cf = 0;
	EMData *parent = this_img->get_parent();

	if (useparent && parent != 0) {
		cf = parent->calc_ccf(to, true);
	}
	else {
		cf = this_img->calc_ccf(to, true);
	}

	float *cf_data = cf->get_data();

	float neg = (float)cf->get_attr("mean") - (float)cf->get_attr("minimum");
	float pos = (float)cf->get_attr("maximum") - (float)cf->get_attr("mean");

	int flag = 1;
	if (neg > pos) {
		flag = -1;
	}

	int nx = this_img->get_xsize();
	int ny = this_img->get_ysize();
	int nz = this_img->get_zsize();

	int peak_x = nx / 2;
	int peak_y = ny / 2;
	int peak_z = nz / 2;

	float max_value = -FLT_MAX;
	float min_value = FLT_MAX;

	int nsec = nx * ny;

	for (int k = nz / 4; k < 3 * nz / 4; k++) {
		for (int j = ny / 4; j < 3 * ny / 4; j++) {
			for (int i = nx / 4; i < 3 * nx / 4; i++) {
				if (cf_data[i + j * nx + k * nsec] * flag > max_value) {
					max_value = cf_data[i + j * nx + k * nsec] * flag;
					peak_x = i;
					peak_y = j;
					peak_z = k;
				}

				if (cf_data[i + j * nx + k * nsec] < min_value) {
					min_value = cf_data[i + j * nx + k * nsec];
				}
			}
		}
	}

	float tx = (float)(nx / 2 - peak_x);
	float ty = (float)(ny / 2 - peak_y);
	float tz = (float)(nz / 2 - peak_z);

	float score = 0;
	if (useparent && parent) {
		Vec3f trans_v = this_img->get_translation();
		score = Util::hypot3(tx - trans_v[0], ty - trans_v[1], tz - trans_v[2]);
	}
	else {
		score = Util::hypot3(tx, ty, tz);
	}
	cf->set_attr("align_score", score);

	if (!to) {
		tx /= 2;
		ty /= 2;
		tz /= 2;
	}

	this_img->translate(tx, ty, tz);

	cf->done_data();

	return cf;
}


EMData *RotationalAligner::align(EMData * this_img, const string&) const
{
	EMData *to = params["to"];

	if (!to) {
		return 0;
	}

	bool premasked = true;
	EMData *this_img2 = this_img->make_rotational_footprint(premasked);
	EMData *to2 = to->make_rotational_footprint(premasked);

	int this_img2_nx = this_img2->get_xsize();

	EMData *cf = this_img2->calc_ccfx(to2, 0, this_img->get_ysize());

	float *data = cf->get_data();
	float peak = 0;
	int peak_index = 0;

	Util::find_max(data, this_img2_nx, &peak, &peak_index);
	this_img->rotate((float)(-peak_index * M_PI / this_img2_nx), 0, 0);
	cf->set_attr("align_score", peak);
	cf->set_attr("rotational",-peak_index * M_PI / this_img2_nx);

	cf->done_data();

	return cf;;
}


EMData *RotatePrecenterAligner::align(EMData * this_img, const string&) const
{
	EMData *to = params["to"];
	if (!to) {
		return 0;
	}

	int ny = this_img->get_ysize();
	int size = Util::calc_best_fft_size((int) (M_PI * ny * 1.5));
	EMData *e1 = this_img->unwrap(4, ny * 7 / 16, size, 0, 0, 1);
	EMData *e2 = to->unwrap(4, ny * 7 / 16, size, 0, 0, 1);
	EMData *cf = e1->calc_ccfx(e2, 0, ny);

	float *data = cf->get_data();

	float peak = 0;
	int peak_index = 0;
	Util::find_max(data, size, &peak, &peak_index);
	float a = (float) ((1.0f - 1.0f * peak_index / size) * M_PI * 2);
	this_img->rotate(a, 0, 0);

	cf->set_attr("align_score", peak);
	cf->set_attr("rotational",a);
	cf->done_data();

	delete e1;
	e1 = 0;

	delete e2;
	e2 = 0;

	return cf;
}


EMData *RotateCHAligner::align(EMData * this_img, const string&) const
{
	static vector < EMData * >ralfp;
	static int rali = 0;
	static int ralo = 0;
	static int rals = 0;

	if (!this_img) {
		return 0;
	}

	EMData *to = params["to"];
	int irad = params.set_default("irad", 8);
	int orad = params.set_default("orad", 0);

	int nx = this_img->get_xsize();
	int ny = this_img->get_ysize();

	const int ralrad = 5;
	const int ralang = 8;

	if (nx != ny) {
		return 0;
	}

	if (irad <= 0) {
		irad = 6;
	}

	if (orad <= 2 || orad > nx / 2) {
		orad = nx / 2 - 2;
	}

	float center = 0;
	int max_num = 0;

	if (nx != rals || irad != rali || orad != ralo) {
		for (size_t i = 0; i < ralfp.size(); i++) {
			delete ralfp[i];
			ralfp[i] = 0;
		}
		ralfp.clear();

		rals = nx;
		ralo = orad;
		rali = irad;

		float wid = (float)((ralo - rali) / (ralrad + 1));

		for (int i = 0; i < ralang; i++) {
			for (int j = 0; j < ralrad; j++) {
				float cen = (j + i) * wid + rali;
				int na = 1 << (i + 1);

				if (cen * M_PI >= na) {
					EMData *d1 = new EMData();
					d1->set_size(nx, ny, 1);
					EMData *d2 = new EMData();
					d2->set_size(nx, ny, 1);

					center = cen;
					max_num = na;

					float *dat1 = d1->get_data();
					float *dat2 = d2->get_data();

					for (int x = 0; x < nx; x++) {
						for (int y = 0; y < ny; y++) {
							float a = 0;
							if (y != ny / 2 || x != nx / 2) {
								a = atan2((float) (y - ny / 2), (float) (x - nx / 2));
							}

							float r = (float)hypot((y - ny / 2), (x - nx / 2)) - cen;
							dat1[x + y * nx] = sin(a * na) * exp(-r * r / (wid * wid / 2.4f));
							dat2[x + y * nx] = cos(a * na) * exp(-r * r / (wid * wid / 2.4f));
						}
					}

					d1->done_data();
					d2->done_data();

					ralfp.push_back(d1);
					ralfp.push_back(d2);
				}
			}
		}
	}

	unsigned int i = 0;
	float ndot = 0;
	float aa = 0;

	for (i = 0; i < ralfp.size(); i += 2) {
		EMData *d1 = ralfp[i];
		EMData *d2 = ralfp[i + 1];

		float ta = this_img->dot(d1);
		float tb = this_img->dot(d2);
		float p1 = 0;
		if (ta != 0 || tb != 0) {
			p1 = atan2(ta, tb);
		}

		ta = to->dot(d1);
		tb = to->dot(d2);

		float a2 = (float)hypot(ta, tb);
		float p2 = 0;
		if (ta != 0 || tb != 0) {
			p2 = atan2(ta, tb);
		}

		float ca = p2 - p1;
		ca = (float)(ca - floor(ca / (2 * M_PI)) * 2 * M_PI);

		if (ca > M_PI) {
			ca -= (float) (2.0 * M_PI);
		}

		if (ndot > 0) {
			float dl = (float) (2.0 * M_PI / max_num);
			float ep = aa / ndot;
			ep = (float)((ep - floor(ep / dl) * dl) * 2.0 * M_PI / dl);
			ca -= ep;
			ca = (float) (ca - floor(ca / (2.0 * M_PI)) * 2.0 * M_PI);

			if (ca > M_PI) {
				ca -= (float) (2.0 * M_PI);
			}

			ca = aa / ndot + ca / max_num;
		}

		aa += ca * center * a2;
		ndot += center * a2;

		printf("%f\t%d\n", ca * 180.0 / M_PI, i);
	}

	printf("%f\t%d\n", aa / ndot * 180.0 / M_PI, i + 5);
	this_img->rotate(aa / ndot, 0, 0);
	this_img->set_attr("align_score", aa / ndot);
	return 0;
}


EMData *RotateTranslateAligner::align(EMData * this_img, const string & cmp_name) const
{
	params.set_default("maxshift", -1);
#if 0
	int usedot = params.set_default("usedot", 0);
	if (usedot) {
		cmp_name = "Dot";
	}
#endif
	EMData *this_copy = this_img->copy();
	this_img->align("Rotational", params);

	EMData *this_copy2 = this_copy->copy();
	this_copy2->set_parent(this_copy->get_parent());

	Dict trans_params;
	trans_params["to"] = params["to"];
	trans_params["useparent"] = 0;
	trans_params["intonly"] = 1;
	trans_params["maxshift"] = params["maxshift"];

	this_copy->align("Translational", trans_params);
	this_copy2->rotate_180();

	this_copy2->align("Translational", trans_params);

	EMData *to = params["to"];

	float dot1 = 0;
	float dot2 = 0;
	Dict cmp_params;
	dot1 = this_copy->cmp(cmp_name, to, cmp_params);
	dot2 = this_copy2->cmp(cmp_name, to, cmp_params);

	EMData *result = 0;
	if (dot1 > dot2) {
		this_copy->set_attr("align_score", dot1);
		delete this_copy2;
		this_copy2 = 0;
		result = this_copy;
	}
	else {
		this_copy2->set_attr("align_score", dot2);
		delete this_copy;
		this_copy = 0;
		result = this_copy2;
	}

	return result;
}


EMData *RotateTranslateBestAligner::align(EMData * this_img, const string & cmp_name) const
{
	params.set_default("maxshift", -1);

	EMData *this_copy = this_img->copy();
	this_img->align("Rotational", params);
	
	Dict rotation = this_img->get_transform().get_rotation(Transform::EMAN);
	float cda = rotation["alt"];

	EMData *this_copy2 = this_copy->copy();
	this_copy2->set_parent(this_copy->get_parent());

	Dict trans_params;
	trans_params["to"] = params["to"];
	trans_params["intonly"] = 0;
	trans_params["useparent"] = 1;
	trans_params["maxshift"] = params["maxshift"];
	this_copy->align("Translational", trans_params);

	Vec3f trans_v = this_copy->get_translation();
	float cdx = trans_v[0] * cos(cda) + trans_v[1] * sin(cda);
	float cdy = trans_v[0] * sin(cda) + trans_v[1] * cos(cda);

	params["alt"] = cda;
	params["dx"] = cdx;
	params["dy"] = cdy;

	this_copy->align("Refine", params);

	float cda2 = cda + (float)M_PI;
	this_copy2->rotate_180();

	this_copy2->align("Translational", trans_params);
	Vec3f trans_v2 = this_copy2->get_translation();

	cdx = trans_v2[0] * cos(cda2) + trans_v2[1] * sin(cda2);
	cdy = -trans_v2[0] * sin(cda2) + trans_v2[1] * cos(cda2);

	params["alt"] = cda2;
	params["dx"] = cdx;
	params["dy"] = cdy;

	this_copy2->align("Refine", params);
	EMData * with = params["to"];
	float dot1 = this_copy->cmp(cmp_name, with, params);
	float dot2 = this_copy2->cmp(cmp_name, with, params);

	EMData *result = 0;
	if (dot1 > dot2) {
		this_copy->set_attr("align_score", dot1);
		delete this_copy2;
		this_copy2 = 0;
		result = this_copy;
	}
	else {
		this_copy2->set_attr("align_score", dot2);
		delete this_copy;
		this_copy = 0;
		result = this_copy2;
	}

	return result;
}



EMData *RotateTranslateRadonAligner::align(EMData * this_img, const string&) const
{
	EMData *to = params["to"];
	int maxshift = params.set_default("maxshift", -1);
	EMData *radonwith = params.set_default("radonwith", (EMData *) 0);
	EMData *radonthis = params.set_default("radonthis", (EMData *) 0);

	int nx = this_img->get_xsize();
	int ny = this_img->get_ysize();
	int size = nx;

	if (nx != ny) {
		LOGERR("%s: images must be square", get_name().c_str());
		return 0;
	}

	if (to && EMUtil::is_same_size(this_img, to)) {
		LOGERR("%s: images must be same size", get_name().c_str());
		return 0;
	}

	float *vert = new float[size];
	int drw = 0;
	int drt = 0;

	if (!radonwith) {
		drw = 1;
		radonwith = to->do_radon();
	}
	if (!radonthis) {
		drt = 1;
		radonthis = this_img->do_radon();
	}
	if (maxshift <= 0) {
		maxshift = size / 8;
	}

	EMData *t1 = radonthis->copy(false);
	radonthis->write_image("radon.hed", 0, EMUtil::IMAGE_IMAGIC);

	float si = 0;
	float co = 0;
	float max = 0;
	float ta = 0;
	int lda = 0;

	for (int j = 0; j < 3; j++) {
		if (j) {
			float *d = radonthis->get_data();
			float *d2 = t1->get_data();

			for (int x = 0; x < size; x++) {
				for (int y = maxshift; y < size - maxshift + 1; y++) {
					d2[x + y * size] =
						d[x + (y + (int) floor(max * sin(x * 2.0 * M_PI / size + ta))) * size];
				}
			}
		}

		t1->write_image("radon.hed", j + 1, EMUtil::IMAGE_IMAGIC);

		EMData *r1 = EMUtil::vertical_acf(t1, maxshift);
		EMData *r2 = EMUtil::vertical_acf(radonwith, maxshift);
		EMData *ccf = r1->calc_ccfx(r2, 0, maxshift);
		r1->write_image("racf.hed", 0, EMUtil::IMAGE_IMAGIC);
		r2->write_image("racf.hed", 1, EMUtil::IMAGE_IMAGIC);

		delete r1;
		r1 = 0;
		delete r2;
		r2 = 0;

		float *d = ccf->get_data();
		float peak_value = 0;
		int peak_x = 0;

		for (int i = 0; i < size; i++) {
			if (d[i] > peak_value) {
				peak_value = d[i];
				peak_x = i % size;
			}
		}

		delete ccf;
		ccf = 0;

		lda = peak_x;
		if (peak_x > size / 2) {
			lda = size - peak_x;
		}

		printf("R Peak %d\t%g\t%1.2f\n", lda, peak_value, lda * 360.0 / size);

		d = radonthis->get_data();
		float *d2 = radonwith->get_data();
		int x2 = lda % size;

		for (int x = 0; x < size / 2; x++) {
			float best = 0;
			int besti = 0;
			for (int i = -maxshift; i < maxshift; i++) {
				float dot = 0;
				for (int y = maxshift; y < size - maxshift; y++) {
					dot += d2[x2 + y * size] * d[x + (y + i) * size];
				}

				if (dot > best) {
					best = dot;
					besti = i;
				}
			}

			vert[x] = (float)besti;
			printf("%d\t%d\n", x, besti);
			x2 = (x2 + 1) % size;
		}

		si = co = max = ta = 0;
		for (int x = 0; x < size / 2; x++) {
			si += (float) sin(x * 2.0 * M_PI / size) * vert[x];
			co += (float) cos(x * 2 * M_PI / size) * vert[x];
			if (fabs(vert[x]) > max) {
				max = fabs(vert[x]);
			}
		}

		float inten = (si * si + co * co) / (size * (float)M_PI);
		ta = atan2(co, si);
		printf("x, y = %g, %g\ta, p=%g / %g, %g\n", co, si, sqrt(inten), max, 180.0f / (float)M_PI * ta);
		max = floor(sqrt(inten) + 0.5f);

		t1->done_data();
		radonwith->done_data();
	}

	delete t1;
	t1 = 0;

	t1 = this_img->copy();

	t1->rotate_translate(-lda * (float)M_PI * 2.0f / size, 0, 0, -max * cos(ta), -max * sin(ta), 0);

	if (drt) {
		delete radonthis;
		radonthis = 0;
	}

	if (drw) {
		delete radonwith;
		radonwith = 0;
	}

	delete[]vert;
	vert = 0;

	return t1;
}



EMData *RotateFlipAligner::align(EMData * this_img, const string&) const
{
	EMData *to = params["to"];
	EMData *flip = params.set_default("to", (EMData *) 0);
	params.set_default("imask", 0);

	EMData *this_copy = this_img->copy();
	this_copy->align("Rotational", params);

	float dot1 = this_copy->dot(to);
	float dot2 = 0;

	EMData *this_copy2 = this_img->copy();

	if (!flip) {
		this_copy2->filter("xform.flip", Dict("axis", "y"));
	}

	this_copy2->align("Rotational", params);
	dot2 = this_copy2->dot(to);

	EMData *result = 0;

	if (!this_copy) {
		result = this_copy2;
	}
	else if (!this_copy2) {
		result = this_copy;
	}
	else {
		if (dot1 > dot2) {
			this_copy->set_flipped(0);
			delete this_copy2;
			this_copy2 = 0;
			result = this_copy;
		}
		else {
			this_copy2->set_flipped(1);
			delete this_copy;
			this_copy = 0;
			result = this_copy2;
		}
	}

	return result;
}

EMData *RotateTranslateFlipAligner::align(EMData * this_img,
										  const string & given_cmp_name) const
{
	EMData *with = params.set_default("to", (EMData *) 0);
	EMData *flip = params.set_default("flip", (EMData *) 0);
	
	params.set_default("maxshift", -1);
	string cmp_name = given_cmp_name;
	int usedot = params.set_default("usedot", 1);
	if (usedot) {
		cmp_name = "Dot";
	}

	EMData *this_copy = this_img->align("RotateTranslate", params);
	EMData *this_copy2 = 0;

	if (flip) {
		this_copy2 = flip->align("RotateTranslate", params);
	}
	else {
		this_img->filter("xform.flip", Dict("axis", "x"));
		this_copy2 = this_img->align("RotateTranslate", params);
	}

	if (!this_copy) {
		LOGERR("%s failed", get_name().c_str());
		return this_copy2;
	}
	if (!this_copy2) {
		LOGERR("%s flip failed", get_name().c_str());
		return this_copy;
	}

	float dot1 = 0;
	float dot2 = 0;

	if (usedot) {
		Dict cmp_params;
		dot1 = this_copy->cmp(cmp_name, with, cmp_params);
		dot2 = this_copy2->cmp(cmp_name, with, cmp_params);

		if (usedot == 2) {
			Vec3f trans = this_copy->get_translation();
			Dict rot = this_copy->get_transform().get_rotation(Transform::EMAN);

			Vec3f trans2 = this_copy2->get_translation();
			Dict rot2 = this_copy2->get_transform().get_rotation(Transform::EMAN);

			printf("%f vs %f  (%1.1f, %1.1f  %1.2f) (%1.1f, %1.1f  %1.2f)\n",
				   dot1, dot2, trans[0], trans[1], (float)rot["alt"] * 180. / M_PI,
				   trans2[0], trans2[1], (float)rot2["alt"] * 180. / M_PI);
		}
	}
	else {
		Dict cmp_params;
		cmp_params["keepzero"] = 1;

		dot1 = this_copy->cmp(cmp_name, with, cmp_params);
		dot2 = this_copy2->cmp(cmp_name, with, cmp_params);
	}

	EMData *result = 0;

	if (dot1 > dot2) {
		this_copy2->set_attr("flipped",0);

		if (!flip) {
			this_img->filter("xform.flip", Dict("axis", "x"));
		}

		delete this_copy2;
		this_copy2 = 0;
		result = this_copy;
	}
	else {
		this_copy2->set_attr("flipped",1);
		delete this_copy;
		this_copy = 0;
		result = this_copy2;
	}

	return result;
}



EMData *RTFSlowAligner::align(EMData * this_img, const string & cmp_name) const
{
	EMData *to = params["to"];
	EMData *flip = params.set_default("flip", (EMData *) 0);
	int maxshift = params.set_default("maxshift", -1);

	EMData *to_copy = to->copy(false);
	int ny = this_img->get_ysize();

	int xst = (int) floor(2 * M_PI * ny);
	xst = Util::calc_best_fft_size(xst);

	to_copy->median_shrink(2);

	int to_copy_r2 = to_copy->get_ysize() / 2 - 2 - maxshift / 2;
	EMData *tmp = to_copy->unwrap(4, to_copy_r2, xst / 2, 0, 0, true);
	delete to_copy;
	to_copy = 0;
	to_copy = tmp;

	EMData *wsc = to_copy->copy(false);
	to = to->unwrap(4, to->get_ysize() / 2 - 2 - maxshift, xst, 0, 0, true);
	EMData *to_copy2 = to->copy(false);

	EMData *df = 0;
	if (flip) {
		df = flip->copy();
	}
	else {
		df = this_img->copy(false);
		df->filter("xform.flip", Dict("axis", "x"));
	}

	EMData *dns = this_img->copy(false);
	EMData *dfs = df->copy(false);
	dns->median_shrink(2);
	dfs->median_shrink(2);

	int dflip = 0;
	float bestval = FLT_MAX;
	float bestang = 0;
	int bestflip = 0;
	int bestdx = 0;
	int bestdy = 0;

	int half_maxshift = maxshift / 2;
	for (dflip = 0; dflip < 2; dflip++) {
		EMData *u = 0;

		if (dflip) {
			u = dfs;
		}
		else {
			u = dns;
		}

		int ur2 = u->get_ysize() / 2 - 2 - half_maxshift;
		
		Dict cmp_params;
		cmp_params["keepzero"] = 1;
		
		for (int dy = -half_maxshift; dy <= half_maxshift; dy++) {
			for (int dx = -half_maxshift; dx <= half_maxshift; dx++) {
				if (hypot(dx, dy) <= half_maxshift) {
					EMData *uw = u->unwrap(4, ur2, xst / 2, dx, dy, true);
					EMData *uwc = uw->copy(false);
					EMData *a = uw->calc_ccfx(to_copy);

					uwc->rotate_x(a->calc_max_index());

					float cm = uwc->cmp(cmp_name, wsc, cmp_params);

					if (cm < bestval) {
						bestval = cm;
						bestang = (float) (2.0 * M_PI * a->calc_max_index() / a->get_xsize());
						bestdx = dx;
						bestdy = dy;
						bestflip = dflip;
					}

					delete a;
					a = 0;
					delete uw;
					uw = 0;
					delete uwc;
					uwc = 0;
				}
			}
		}
	}
	delete dns;
	dns = 0;
	delete dfs;
	dfs = 0;
	delete to_copy;
	to_copy = 0;
	delete wsc;
	wsc = 0;

	bestdx *= 2;
	bestdy *= 2;
	bestval = FLT_MAX;

	int bestdx2 = bestdx;
	int bestdy2 = bestdy;


	for (dflip = 0; dflip < 2; dflip++) {
		EMData *u = 0;

		if (dflip) {
			u = df;
			bestdx2 = -bestdx2;
			bestang = -bestang;
		}
		else {
			u = this_img;
		}


		for (int dy = bestdy2 - 3; dy <= bestdy2 + 3; dy++) {
			for (int dx = bestdx2 - 3; dx <= bestdx2 + 3; dx++) {
				if (hypot(dx, dy) <= maxshift) {
					EMData *uw = u->unwrap(4, u->get_ysize() / 2 - 2 - maxshift, xst, dx, dy, true);
					EMData *uwc = uw->copy(false);
					EMData *a = uw->calc_ccfx(to);

					uwc->rotate_x(a->calc_max_index());
					Dict cmp_params;
					cmp_params["keepzero"] = 1;
					float cm = uwc->cmp(cmp_name, to_copy2, cmp_params);

					if (cm < bestval) {
						bestval = cm;
						bestang = (float)(2.0 * M_PI * a->calc_max_index() / a->get_xsize());
						bestdx = dx;
						bestdy = dy;
						bestflip = dflip;
					}
					delete a;
					a = 0;
					delete uw;
					uw = 0;
					delete uwc;
					uwc = 0;
				}
			}
		}
	}

	delete to;
	to = 0;
	delete to_copy;
	to_copy = 0;

	if (bestflip) {
		df->rotate_translate((float)bestang, 0.0f, 0.0f, (float)-bestdx, (float)-bestdy, 0.0f);
		df->set_flipped(1);
		return df;
	}

	delete df;
	df = 0;

	EMData *dn = 0;
	if (dflip) {
		dn = this_img->copy();
	}
	else {
		dn = this_img->copy(false);
	}

	dn->rotate_translate((float)bestang, 0.0f, 0.0f, (float)-bestdx, (float)-bestdy, 0.0f);

	return dn;
}


EMData *RTFSlowestAligner::align(EMData * this_img, const string & cmp_name) const
{
	EMData *to = params["to"];
	EMData *flip = params.set_default("flip", (EMData *) 0);
	int maxshift = params.set_default("maxshift", -1);

	EMData *dn = this_img->copy();
	EMData *df = 0;

	if (flip) {
		df = flip->copy();
	}
	else {
		df = this_img->copy();
		df->filter("xform.flip", Dict("axis", "x"));
		df = df->copy();
	}

	int nx = this_img->get_xsize();

	if (maxshift < 0) {
		maxshift = nx / 10;
	}

	float astep = atan2(2.0f, (float)nx);

	EMData *dns = dn->copy(false);
	EMData *dfs = df->copy(false);
	EMData *to_copy = to->copy(false);

	dns->median_shrink(2);
	dfs->median_shrink(2);
	to_copy->median_shrink(2);
	dns = dns->copy();
	dfs = dfs->copy();

	int bestflip = 0;
	int bestdx = 0;
	int bestdy = 0;

	float bestang = 0;
	float bestval = FLT_MAX;

	int dflip = 0;
	int half_maxshift = maxshift / 2;

	for (dflip = 0; dflip < 2; dflip++) {
		EMData *u = 0;

		if (dflip) {
			u = dfs;
		}
		else {
			u = dns;
		}

		for (int dy = -half_maxshift; dy <= half_maxshift; dy++) {
			for (int dx = -half_maxshift; dx <= half_maxshift; dx++) {
				if (hypot(dx, dy) <= maxshift) {
					for (float ang = -astep * 2.0f; ang <= (float)2 * M_PI; ang += astep * 4.0f) {
						u->rotate_translate(ang, 0.0f, 0.0f, (float)dx, (float)dy, 0.0f);

						Dict cmp_params;
						float lc = u->cmp(cmp_name, to_copy, cmp_params);

						if (lc < bestval) {
							bestval = lc;
							bestang = ang;
							bestdx = dx;
							bestdy = dy;
							bestflip = dflip;
						}
					}
				}
			}
		}
	}

	delete dns->get_parent();
	delete dfs->get_parent();
	delete dns;
	dns = 0;
	delete dfs;
	dfs = 0;
	delete to_copy;
	to_copy = 0;

	bestdx *= 2;
	bestdy *= 2;
	bestval = FLT_MAX;

	int bestdx2 = bestdx;
	int bestdy2 = bestdy;
	float bestang2 = bestang;

	for (dflip = 0; dflip < 2; dflip++) {
		EMData *u = 0;
		if (dflip) {
			u = df;
		}
		else {
			u = dn;
		}

		if (dflip) {
			bestdx2 = -bestdx2;
		}

		for (int dy = bestdy2 - 3; dy <= bestdy2 + 3; dy++) {
			for (int dx = bestdx2 - 3; dx <= bestdx2 + 3; dx++) {
				if (hypot(dx, dy) <= maxshift) {
					for (float ang = bestang2 - astep * 6.0f; ang <= bestang2 + astep * 6.0f;
						 ang += astep) {
						u->rotate_translate(ang, 0.0f, 0.0f, (float)dx, (float)dy, 0.0f);

						Dict cmp_params;
						cmp_params["keepzero"] = 1;
						float lc = u->cmp(cmp_name, to_copy, cmp_params);

						if (lc < bestval) {
							bestval = lc;
							bestang = ang;
							bestdx = dx;
							bestdy = dy;
							bestflip = dflip;
						}
					}
				}
			}
		}
	}

	if (bestflip) {
		delete dn;
		dn = 0;

		df->rotate_translate(bestang, 0.0f, 0.0f, (float)bestdx, (float)bestdy, 0.0f);

		if (!dflip) {
			delete df->get_parent();
			df->set_parent(0);
		}
		df->set_flipped(1);
		return df;
	}

	dn->rotate_translate(bestang, 0.0f, 0.0f, (float)bestdx, (float)bestdy, 0.0f);

	if (!dflip) {
		delete df->get_parent();
	}

	delete df;
	df = 0;

	return dn;
}

EMData *RTFBestAligner::align(EMData * this_img, const string & cmp_name) const
{
	EMData *flip = params.set_default("flip", (EMData *) 0);
	params.set_default("maxshift", -1);


	EMData *this_copy = this_img->align("RotateTranslateBest", params);
	EMData *flip_copy = 0;

	if (!flip) {
		this_img->filter("xform.flip", Dict("axis", "x"));
		flip_copy = this_img->align("RotateTranslateBest", params);
	}
	else {
		flip_copy = flip->align("RotateTranslateBest", params);
	}

	if (!this_copy) {
		LOGERR("%s align failed", get_name().c_str());
		return flip_copy;
	}
	else if (!flip_copy) {
		LOGERR("%s align failed", get_name().c_str());
		return this_copy;
	}

	EMData * with = params["to"];
	float this_cmp = this_copy->cmp(cmp_name, with, params);
	float flip_cmp = flip_copy->cmp(cmp_name, with, params);

	EMData *result = 0;

	if (this_cmp > flip_cmp) {
		this_copy->set_flipped(0);
		if (!flip) {
			this_img->filter("xform.flip", Dict("axis", "x"));
		}
		delete flip_copy;
		flip_copy = 0;
		result = this_copy;
	}
	else {
		flip_copy->set_flipped(1);
		delete this_copy;
		this_copy = 0;
		result = flip_copy;
	}

	return result;
}


EMData *RTFRadonAligner::align(EMData * this_img, const string&) const
{
	EMData *to = params["to"];
	params.set_default("maxshift", -1);
	EMData *thisf = params.set_default("thisf", (EMData *) 0);
	EMData *radonwith = params.set_default("radonwith", (EMData *) 0);
	params.set_default("radonthis", (EMData *) 0);
	params.set_default("radonthisf", (EMData *) 0);

	int drw = 0;
	if (!radonwith) {
		drw = 1;
		radonwith = to->do_radon();
		params["radonwith"] = radonwith;
	}

	EMData *r1 = this_img->align("RotateTranslateRadon", params);
	EMData *r2 = 0;

	if (!thisf) {
		this_img->filter("xform.flip", Dict("axis", "x"));
		r2 = this_img->align("RTFRadon", params);
		this_img->filter("xform.flip", Dict("axis", "x"));
	}
	else {
		r2 = thisf->align("RTFRadon", params);
	}

	float r1_score = r1->dot(to);
	float r2_score = r2->dot(to);

	if (drw) {
		delete radonwith;
		radonwith = 0;
	}

	EMData *result = 0;
	if (r1_score < r2_score) {
		delete r1;
		r1 = 0;

		if (!thisf) {
			this_img->filter("xform.flip", Dict("axis", "x"));
		}
		result = r2;
	}
	else {
		delete r2;
		r2 = 0;
		result = r1;
	}

	return result;
}


static double refalifn(const gsl_vector * v, void *params)
{
	Dict *dict = (Dict *) params;

	double x = gsl_vector_get(v, 0);
	double y = gsl_vector_get(v, 1);
	double a = gsl_vector_get(v, 2);

	EMData *this_img = (*dict)["this"];
	EMData * with = (*dict)["with"];
	this_img->rotate_translate((float)a, 0.0f, 0.0f, (float)x, (float)y, 0.0f);

	return this_img->cmp("FRC", with, *dict);
}

static double refalifnfast(const gsl_vector * v, void *params)
{
	Dict *dict = (Dict *) params;
	EMData *this_img = (*dict)["this"];
	EMData *img_to = (*dict)["to"];

	double x = gsl_vector_get(v, 0);
	double y = gsl_vector_get(v, 1);
	double a = gsl_vector_get(v, 2);

	double r = this_img->dot_rotate_translate(img_to, (float)x, (float)y, (float)a);
	int nsec = this_img->get_xsize() * this_img->get_ysize();
	double result = 1.0 - r / nsec;

	return result;
}


EMData *RefineAligner::align(EMData * this_img, const string & cmp_name) const
{
	EMData *to = params["to"];
	if (!to) {
		return 0;
	}

	EMData *result = 0;

	int ny = this_img->get_ysize();

	float salt = params["alt"];
	float sdx = params["dx"];
	float sdy = params["dy"];
	float sdz = params["dz"];

	float dda = atan(2.0f / ny);

	int mode = params.set_default("mode", 0);

	if (this_img->get_parent() == 0) {
		LOGWARN("%s: no parent", get_name().c_str());
	}

	if (mode == 1 || mode == 2) {
		int np = 3;
		Dict gsl_params;
		gsl_params["this"] = this_img;
		gsl_params["with"] = to;
		gsl_params["snr"] = params["snr"];

		const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex;

		gsl_vector *ss = gsl_vector_alloc(np);
		gsl_vector_set(ss, 0, 1.0f);
		gsl_vector_set(ss, 1, 1.0f);
		gsl_vector_set(ss, 2, 0.1f);

		gsl_vector *x = gsl_vector_alloc(np);
		gsl_vector_set(x, 0, sdx);
		gsl_vector_set(x, 1, sdy);
		gsl_vector_set(x, 2, salt);

		gsl_multimin_function minex_func;

		if (mode == 2) {
			minex_func.f = &refalifnfast;
		}
		else {
			minex_func.f = &refalifn;
		}
		minex_func.n = np;
		minex_func.params = (void *) &gsl_params;

		gsl_multimin_fminimizer *s = gsl_multimin_fminimizer_alloc(T, np);
		gsl_multimin_fminimizer_set(s, &minex_func, x, ss);


		int rval = GSL_CONTINUE;
		int status = GSL_SUCCESS;
		int iter = 1;

		while (rval == GSL_CONTINUE && iter < 28) {
			iter++;
			status = gsl_multimin_fminimizer_iterate(s);
			if (status) {
				break;
			}

			rval = gsl_multimin_test_size(gsl_multimin_fminimizer_size(s), 0.04f);
		}

		this_img->rotate_translate((float)gsl_vector_get(s->x, 2), 0, 0,
								   (float)gsl_vector_get(s->x, 0), (float)gsl_vector_get(s->x, 1), 0);
		gsl_vector_free(x);
		gsl_vector_free(ss);
		gsl_multimin_fminimizer_free(s);
	}
	else {
		float best = 0;
		Dict cmp_params;
		
		if (mode == 0) {
			best = this_img->cmp(cmp_name, to, cmp_params);
		}
		else if (mode < 0) {
			printf("Start %f  %f,%f\n", salt * 180.0 / M_PI, sdx, sdy);
		}

		int j = 0;
		float f = 0;

		for (int i = 0; i < 5 && j; i++) {
			j = 0;
			float last_alt = salt;
			for (float daz = salt - dda; daz <= salt + dda; daz += dda / 2.0f) {
				if (daz != salt) {
					last_alt = daz;
					this_img->rotate_translate(daz, 0, 0, sdx, sdy, sdz);
					f = this_img->cmp(cmp_name, to, cmp_params);

					if (f > best) {
						best = f;
						salt = daz;
						j = 1;
						break;
					}
				}
			}

			for (float dx = sdx - 1.0f; dx <= sdx + 1.0; dx += 0.25f) {
				for (float dy = sdy - 1.0f; dy <= sdy + 1.0; dy += 0.25f) {
					this_img->rotate_translate(last_alt, 0, 0, dx, dy, 0);
					f = this_img->cmp(cmp_name, to, cmp_params);

					if (f > best) {
						best = f;
						sdx = dx;
						sdy = dy;
						j = 1;
						break;
					}
				}
			}

			if (mode < 0) {
				printf("-- %f  %f,%f (%1.3g) \n", salt * 180.0 / M_PI, sdx, sdy, best);
			}
		}

		this_img->rotate_translate(salt, 0, 0, sdx, sdy, 0);

		result = this_img->copy();
		result->set_attr("align_score", best);
	}

	return result;
}

void EMAN::dump_aligners()
{
	dump_factory < Aligner > ();
}
