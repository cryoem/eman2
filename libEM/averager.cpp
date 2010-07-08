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

#include "averager.h"
#include "emdata.h"
#include "xydata.h"
#include "ctf.h"
#include <cstring>
#include "plugins/averager_template.h"

using namespace EMAN;

const string ImageAverager::NAME = "mean";
const string MinMaxAverager::NAME = "minmax";
const string IterationAverager::NAME = "iteration";
const string WeightingAverager::NAME = "snr_weight";
const string CtfCAverager::NAME = "ctfc";
const string CtfCWAverager::NAME = "ctfcw";
const string CtfCWautoAverager::NAME = "ctfw.auto";
const string CtfCAutoAverager::NAME = "ctf.auto";

template <> Factory < Averager >::Factory()
{
	force_add<ImageAverager>();
	force_add<MinMaxAverager>();
	force_add<IterationAverager>();
//	force_add<WeightingAverager>();
	// These commented out until we're happy they're working. (d.woolford, Feb 3rd 2009)
//	force_add(&CtfCAverager::NEW);
//	force_add(&CtfCWAverager::NEW);
	force_add<CtfCWautoAverager>();
	force_add<CtfCAutoAverager>();
//	force_add<XYZAverager>();
}

void Averager::mult(const float& s)
{
	if ( result != 0 )
	{
		result->mult(s);
	}
	else throw NullPointerException("Error, attempted to multiply the result image, which is NULL");
}


void Averager::add_image_list(const vector<EMData*> & image_list)
{
	for (size_t i = 0; i < image_list.size(); i++) {
		add_image(image_list[i]);
	}
}

ImageAverager::ImageAverager()
	: sigma_image(0), nimg_n0(0), ignore0(0), nimg(0)
{

}

void ImageAverager::add_image(EMData * image)
{
	if (!image) {
		return;
	}

	if (nimg >= 1 && !EMUtil::is_same_size(image, result)) {
		LOGERR("%sAverager can only process same-size Image",
			   get_name().c_str());
		return;
	}

	nimg++;

	int nx = image->get_xsize();
	int ny = image->get_ysize();
	int nz = image->get_zsize();
	size_t image_size = nx * ny * nz;

	if (nimg == 1) {
		result = new EMData();
		result->set_size(nx, ny, nz);
		sigma_image = params.set_default("sigma", (EMData*)0);
		ignore0 = params["ignore0"];

		nimg_n0 = new int[image_size];
		for (size_t i = 0; i < image_size; i++) {
			nimg_n0[i] = 0;
		}
	}

	float *result_data = result->get_data();
	float *sigma_image_data = 0;
	if (sigma_image) {
		sigma_image->set_size(nx, ny, nz);
		sigma_image_data = sigma_image->get_data();
	}

	float * image_data = image->get_data();

	if (!ignore0) {
		for (size_t j = 0; j < image_size; j++) {
			float f = image_data[j];
			result_data[j] += f;
			if (sigma_image_data) {
				sigma_image_data[j] += f * f;
			}
		}
	}
	else {
		for (size_t j = 0; j < image_size; j++) {
			float f = image_data[j];
			if (f) {
				result_data[j] += f;
				if (sigma_image_data) {
					sigma_image_data[j] += f * f;
				}
				nimg_n0[j]++;
			}
		}
	}
}

EMData * ImageAverager::finish()
{
	if (result && nimg > 1) {
		size_t image_size = result->get_xsize() * result->get_ysize() * result->get_zsize();
		float * result_data = result->get_data();

		if (!ignore0) {
			for (size_t j = 0; j < image_size; j++) {
				result_data[j] /= nimg;
			}

			if (sigma_image) {
				float * sigma_image_data = sigma_image->get_data();

				for (size_t j = 0; j < image_size; j++) {
					float f1 = sigma_image_data[j] / nimg;
					float f2 = result_data[j];
					sigma_image_data[j] = sqrt(f1 - f2 * f2);
				}

				sigma_image->update();
			}
		}
		else {
			for (size_t j = 0; j < image_size; j++) {
				result_data[j] /= nimg_n0[j];
			}
			if (sigma_image) {
				float * sigma_image_data = sigma_image->get_data();

				for (size_t j = 0; j < image_size; j++) {
					float f1 = sigma_image_data[j] / nimg_n0[j];
					float f2 = result_data[j];
					sigma_image_data[j] = sqrt(f1 - f2 * f2);
				}

				sigma_image->update();
			}
		}

		result->update();
	}

	if( nimg_n0 )
	{
		delete [] nimg_n0;
		nimg_n0 = 0;
	}

	return result;
}

#if 0
EMData *ImageAverager::average(const vector < EMData * >&image_list) const
{
	if (image_list.size() == 0) {
		return 0;
	}

	EMData *sigma_image = params["sigma"];
	int ignore0 = params["ignore0"];

	EMData *image0 = image_list[0];

	int nx = image0->get_xsize();
	int ny = image0->get_ysize();
	int nz = image0->get_zsize();
	size_t image_size = nx * ny * nz;

	EMData *result = new EMData();
	result->set_size(nx, ny, nz);

	float *result_data = result->get_data();
	float *sigma_image_data = 0;

	if (sigma_image) {
		sigma_image->set_size(nx, ny, nz);
		sigma_image_data = sigma_image->get_data();
	}

	int c = 1;
	if (ignore0) {
		for (size_t j = 0; j < image_size; j++) {
			int g = 0;
			for (size_t i = 0; i < image_list.size(); i++) {
				float f = image_list[i]->get_value_at(j);
				if (f) {
					g++;
					result_data[j] += f;
					if (sigma_image_data) {
						sigma_image_data[j] += f * f;
					}
				}
			}
			if (g) {
				result_data[j] /= g;
			}
		}
	}
	else {
		float *image0_data = image0->get_data();
		if (sigma_image_data) {
			memcpy(sigma_image_data, image0_data, image_size * sizeof(float));

			for (size_t j = 0; j < image_size; j++) {
				sigma_image_data[j] *= sigma_image_data[j];
			}
		}

		image0->update();
		memcpy(result_data, image0_data, image_size * sizeof(float));

		for (size_t i = 1; i < image_list.size(); i++) {
			EMData *image = image_list[i];

			if (EMUtil::is_same_size(image, result)) {
				float *image_data = image->get_data();

				for (size_t j = 0; j < image_size; j++) {
					result_data[j] += image_data[j];
				}

				if (sigma_image_data) {
					for (size_t j = 0; j < image_size; j++) {
						sigma_image_data[j] += image_data[j] * image_data[j];
					}
				}

				image->update();
				c++;
			}
		}

		for (size_t j = 0; j < image_size; j++) {
			result_data[j] /= static_cast < float >(c);
		}
	}

	if (sigma_image_data) {
		for (size_t j = 0; j < image_size; j++) {
			float f1 = sigma_image_data[j] / c;
			float f2 = result_data[j];
			sigma_image_data[j] = sqrt(f1 - f2 * f2);
		}
	}

	if (sigma_image_data) {
		sigma_image->update();
	}

	result->update();
	return result;
}
#endif

MinMaxAverager::MinMaxAverager()
	: nimg(0)
{
	/*move max out of initializer list, since this max(0) is considered as a macro
	 * in Visual Studio, which we defined somewhere else*/
	max = 0;
}

void MinMaxAverager::add_image(EMData * image)
{
	if (!image) {
		return;
	}

	if (nimg >= 1 && !EMUtil::is_same_size(image, result)) {
		LOGERR("%sAverager can only process same-size Image",
			   get_name().c_str());
		return;
	}

	nimg++;

	int nx = image->get_xsize();
	int ny = image->get_ysize();
	int nz = image->get_zsize();

	if (nimg == 1) {
		result = image->copy();
		max = params["max"];
		return;
	}

	for (int z=0; z<nz; z++) {
		for (int y=0; y<ny; y++) {
			for (int x=0; x<nx; x++) {
				if (result->get_value_at(x,y,z)>image->get_value_at(x,y,z))
					{ if (!max) result->set_value_at(x,y,z,image->get_value_at(x,y,z)); }
				else { if (max) result->set_value_at(x,y,z,image->get_value_at(x,y,z)); }
			}
		}
	}

}

EMData *MinMaxAverager::finish()
{
	result->update();
	if (result && nimg > 1) return result;

	return NULL;
}


IterationAverager::IterationAverager() : nimg(0)
{

}

void IterationAverager::add_image( EMData * image)
{
	if (!image) {
		return;
	}

	if (nimg >= 1 && !EMUtil::is_same_size(image, result)) {
		LOGERR("%sAverager can only process same-size Image",
							 get_name().c_str());
		return;
	}

	nimg++;

	int nx = image->get_xsize();
	int ny = image->get_ysize();
	int nz = image->get_zsize();
	size_t image_size = nx * ny * nz;

	if (nimg == 1) {
		result = new EMData();
		result->set_size(nx, ny, nz);
		sigma_image = new EMData();
		sigma_image->set_size(nx, ny, nz);
	}

	float *image_data = image->get_data();
	float *result_data = result->get_data();
	float *sigma_image_data = sigma_image->get_data();

	for (size_t j = 0; j < image_size; j++) {
		float f = image_data[j];
		result_data[j] += f;
		sigma_image_data[j] += f * f;
	}


}

EMData * IterationAverager::finish()
{
	if (nimg < 1) {
		return result;
	}

	int nx = result->get_xsize();
	int ny = result->get_ysize();
	int nz = result->get_zsize();
	size_t image_size = nx * ny * nz;

	float *result_data = result->get_data();
	float *sigma_image_data = sigma_image->get_data();

	for (size_t j = 0; j < image_size; j++) {
		result_data[j] /= nimg;
		float f1 = sigma_image_data[j] / nimg;
		float f2 = result_data[j];
		sigma_image_data[j] = sqrt(f1 - f2 * f2) / sqrt((float)nimg);
	}

	result->update();
	sigma_image->update();

	result->append_image("iter.hed");
	float sigma = sigma_image->get_attr("sigma");
	float *sigma_image_data2 = sigma_image->get_data();
	float *result_data2 = result->get_data();
	float *d2 = new float[nx * ny];
	size_t sec_size = nx * ny * sizeof(float);

	memcpy(d2, result_data2, sec_size);
	memcpy(sigma_image_data2, result_data2, sec_size);

	printf("Iter sigma=%f\n", sigma);

	for (int k = 0; k < 1000; k++) {
		for (int i = 1; i < nx - 1; i++) {
			for (int j = 1; j < ny - 1; j++) {
				int l = i + j * nx;
				float c1 = (d2[l - 1] + d2[l + 1] + d2[l - nx] + d2[l + nx]) / 4.0f - d2[l];
				float c2 = fabs(result_data2[l] - sigma_image_data2[l]) / sigma;
				result_data2[l] += c1 * Util::eman_erfc(c2) / 100.0f;
			}
		}

		memcpy(d2, result_data2, sec_size);
	}

	if( d2 )
	{
		delete[]d2;
		d2 = 0;
	}

	sigma_image->update();
	if( sigma_image )
	{
		delete sigma_image;
		sigma_image = 0;
	}

	result->update();
	result->append_image("iter.hed");


	return result;
}

CtfCWautoAverager::CtfCWautoAverager()
	: nimg(0)
{

}


void CtfCWautoAverager::add_image(EMData * image)
{
	if (!image) {
		return;
	}



	EMData *fft=image->do_fft();

	if (nimg >= 1 && !EMUtil::is_same_size(fft, result)) {
		LOGERR("%s Averager can only process images of the same size", get_name().c_str());
		return;
	}

	nimg++;
	if (nimg == 1) {
		result = fft->copy_head();
		result->to_zero();
	}

	Ctf *ctf = (Ctf *)image->get_attr("ctf");
	float b=ctf->bfactor;
	ctf->bfactor=100.0;			// FIXME - this is a temporary fixed B-factor which does a (very) little sharpening

	EMData *snr = result -> copy();
	ctf->compute_2d_complex(snr,Ctf::CTF_SNR);
	EMData *ctfi = result-> copy();
	ctf->compute_2d_complex(ctfi,Ctf::CTF_AMP);

	ctf->bfactor=b;	// return to its original value

	float *outd = result->get_data();
	float *ind = fft->get_data();
	float *snrd = snr->get_data();
	float *ctfd = ctfi->get_data();

	size_t sz=snr->get_xsize()*snr->get_ysize();
	for (size_t i = 0; i < sz; i+=2) {
		if (snrd[i]<0) snrd[i]=0;
		ctfd[i]=fabs(ctfd[i]);
		if (ctfd[i]<.05) {
			if (snrd[i]<=0) ctfd[i]=.05f;
			else ctfd[i]=snrd[i]*10.0f;
		}
		outd[i]+=ind[i]*snrd[i]/ctfd[i];
		outd[i+1]+=ind[i+1]*snrd[i]/ctfd[i];
	}

	if (nimg==1) {
		snrsum=snr->copy_head();
		float *ssnrd=snrsum->get_data();
		// we're only using the real component, and we need to start with 1.0
		for (size_t i = 0; i < sz; i+=2) { ssnrd[i]=1.0; ssnrd[i+1]=0.0; }
	}
	snrsum->add(*snr);

	delete ctf;
	delete fft;
	delete snr;
	delete ctfi;
}

EMData * CtfCWautoAverager::finish()
{
/*	EMData *tmp=result->do_ift();
	tmp->write_image("ctfcw.hdf",0);
	delete tmp;

	tmp=snrsum->do_ift();
	tmp->write_image("ctfcw.hdf",1);
	delete tmp;*/

//	snrsum->write_image("snrsum.hdf",-1);
	size_t sz=result->get_xsize()*result->get_ysize();
	float *snrsd=snrsum->get_data();
	float *outd=result->get_data();

	for (size_t i=0; i<sz; i+=2) {
		outd[i]/=snrsd[i];		// snrsd contains total SNR+1
		outd[i+1]/=snrsd[i];
	}
	result->update();
	result->set_attr("ptcl_repr",nimg);
	result->set_attr("ctf_snr_total",snrsum->calc_radial_dist(snrsum->get_ysize()/2,0,1,false));
	result->set_attr("ctf_wiener_filtered",true);
	
	delete snrsum;
	EMData *ret=result->do_ift();
	delete result;
	result=NULL;
	return ret;
}

CtfCAutoAverager::CtfCAutoAverager()
	: nimg(0)
{

}


void CtfCAutoAverager::add_image(EMData * image)
{
	if (!image) {
		return;
	}



	EMData *fft=image->do_fft();

	if (nimg >= 1 && !EMUtil::is_same_size(fft, result)) {
		LOGERR("%s Averager can only process images of the same size", get_name().c_str());
		return;
	}

	nimg++;
	if (nimg == 1) {
		result = fft->copy_head();
		result->to_zero();
	}

	Ctf *ctf = (Ctf *)image->get_attr("ctf");
	float b=ctf->bfactor;
	ctf->bfactor=0;			// NO B-FACTOR CORRECTION !

	EMData *snr = result -> copy();
	ctf->compute_2d_complex(snr,Ctf::CTF_SNR);
	EMData *ctfi = result-> copy();
	ctf->compute_2d_complex(ctfi,Ctf::CTF_AMP);

	ctf->bfactor=b;	// return to its original value

	float *outd = result->get_data();
	float *ind = fft->get_data();
	float *snrd = snr->get_data();
	float *ctfd = ctfi->get_data();

	size_t sz=snr->get_xsize()*snr->get_ysize();
	for (size_t i = 0; i < sz; i+=2) {
		if (snrd[i]<0) snrd[i]=0;
		ctfd[i]=fabs(ctfd[i]);
		
		// This limits the maximum possible amplification in CTF correction to 10x
		if (ctfd[i]<.05) {
			if (snrd[i]<=0) ctfd[i]=.05f;
			else ctfd[i]=snrd[i]*10.0f;
		}
		
		// SNR weight and CTF correction
		outd[i]+=ind[i]*snrd[i]/ctfd[i];
		outd[i+1]+=ind[i+1]*snrd[i]/ctfd[i];
	}

	if (nimg==1) {
		snrsum=snr->copy_head();
		float *ssnrd=snrsum->get_data();
		// we're only using the real component, for Wiener filter we put 1.0 in R, but for just SNR weight we use 0
		for (size_t i = 0; i < sz; i+=2) { ssnrd[i]=0.0; ssnrd[i+1]=0.0; }
	}
	snrsum->add(*snr);

	delete ctf;
	delete fft;
	delete snr;
	delete ctfi;
}

EMData * CtfCAutoAverager::finish()
{
/*	EMData *tmp=result->do_ift();
	tmp->write_image("ctfcw.hdf",0);
	delete tmp;

	tmp=snrsum->do_ift();
	tmp->write_image("ctfcw.hdf",1);
	delete tmp;*/

//	snrsum->write_image("snrsum.hdf",-1);
	size_t sz=result->get_xsize()*result->get_ysize();
	float *snrsd=snrsum->get_data();
	float *outd=result->get_data();

	for (size_t i=0; i<sz; i+=2) {
		outd[i]/=snrsd[i];		// snrsd contains total SNR
		outd[i+1]/=snrsd[i];
	}
	result->update();
	result->set_attr("ptcl_repr",nimg);
	result->set_attr("ctf_snr_total",snrsum->calc_radial_dist(snrsum->get_ysize()/2,0,1,false));
	result->set_attr("ctf_wiener_filtered",false);
	
	delete snrsum;
	EMData *ret=result->do_ift();
	delete result;
	result=NULL;
	return ret;
}

#if 0
EMData *IterationAverager::average(const vector < EMData * >&image_list) const
{
	if (image_list.size() == 0) {
		return 0;
	}

	EMData *image0 = image_list[0];

	int nx = image0->get_xsize();
	int ny = image0->get_ysize();
	int nz = image0->get_zsize();
	size_t image_size = nx * ny * nz;

	EMData *result = new EMData();
	result->set_size(nx, ny, nz);

	EMData *sigma_image = new EMData();
	sigma_image->set_size(nx, ny, nz);

	float *image0_data = image0->get_data();
	float *result_data = result->get_data();
	float *sigma_image_data = sigma_image->get_data();

	memcpy(result_data, image0_data, image_size * sizeof(float));
	memcpy(sigma_image_data, image0_data, image_size * sizeof(float));

	for (size_t j = 0; j < image_size; j++) {
		sigma_image_data[j] *= sigma_image_data[j];
	}

	image0->update();

	int nc = 1;
	for (size_t i = 1; i < image_list.size(); i++) {
		EMData *image = image_list[i];

		if (EMUtil::is_same_size(image, result)) {
			float *image_data = image->get_data();

			for (size_t j = 0; j < image_size; j++) {
				result_data[j] += image_data[j];
				sigma_image_data[j] += image_data[j] * image_data[j];
			}

			image->update();
			nc++;
		}
	}

	float c = static_cast < float >(nc);
	for (size_t j = 0; j < image_size; j++) {
		float f1 = sigma_image_data[j] / c;
		float f2 = result_data[j] / c;
		sigma_image_data[j] = sqrt(f1 - f2 * f2) / sqrt(c);
	}


	for (size_t j = 0; j < image_size; j++) {
		result_data[j] /= c;
	}

	result->update();
	sigma_image->update();

	result->append_image("iter.hed");

	float sigma = sigma_image->get_attr("sigma");
	float *sigma_image_data2 = sigma_image->get_data();
	float *result_data2 = result->get_data();
	float *d2 = new float[nx * ny];
	size_t sec_size = nx * ny * sizeof(float);

	memcpy(d2, result_data2, sec_size);
	memcpy(sigma_image_data2, result_data2, sec_size);

	printf("Iter sigma=%f\n", sigma);

	for (int k = 0; k < 1000; k++) {
		for (int i = 1; i < nx - 1; i++) {
			for (int j = 1; j < ny - 1; j++) {
				int l = i + j * nx;
				float c1 = (d2[l - 1] + d2[l + 1] + d2[l - nx] + d2[l + nx]) / 4.0f - d2[l];
				float c2 = fabs(result_data2[l] - sigma_image_data2[l]) / sigma;
				result_data2[l] += c1 * Util::eman_erfc(c2) / 100.0f;
			}
		}

		memcpy(d2, result_data2, sec_size);
	}

	if( d2 )
	{
		delete[]d2;
		d2 = 0;
	}

	sigma_image->update();
	if( sigma_image )
	{
		delete sigma_image;
		sigma_image = 0;
	}

	result->update();
	result->append_image("iter.hed");

	return result;
}
#endif


CtfAverager::CtfAverager() :
	sf(0), curves(0), need_snr(false), outfile(0),
	image0_fft(0), image0_copy(0), snri(0), snrn(0),
	tdr(0), tdi(0), tn(0),
	filter(0), nimg(0), nx(0), ny(0), nz(0)
{

}

void CtfAverager::add_image(EMData * image)
{
	if (!image) {
		return;
	}

	if (nimg >= 1 && !EMUtil::is_same_size(image, result)) {
		LOGERR("%sAverager can only process same-size Image",
							 get_name().c_str());
		return;
	}

	if (image->get_zsize() != 1) {
		LOGERR("%sAverager: Only 2D images are currently supported",
							 get_name().c_str());
	}

	string alg_name = get_name();

	if (alg_name == "CtfCW" || alg_name == "CtfCWauto") {
		if (image->get_ctf() != 0 && !image->has_ctff()) {
			LOGERR("%sAverager: Attempted CTF Correction with no ctf parameters",
								 get_name().c_str());
		}
	}
	else {
		if (image->get_ctf() != 0) {
			LOGERR("%sAverager: Attempted CTF Correction with no ctf parameters",
								 get_name().c_str());
		}
	}

	nimg++;


	if (nimg == 1) {
		image0_fft = image->do_fft();

		nx = image0_fft->get_xsize();
		ny = image0_fft->get_ysize();
		nz = image0_fft->get_zsize();

		result = new EMData();
		result->set_size(nx - 2, ny, nz);


		if (alg_name == "Weighting" && curves) {
			if (!sf) {
				LOGWARN("CTF curve in file will contain relative, not absolute SNR!");
			}
			curves->set_size(Ctf::CTFOS * ny / 2, 3, 1);
			curves->to_zero();
		}


		if (alg_name == "CtfC") {
			filter = params["filter"];
			if (filter == 0) {
				filter = 22.0f;
			}
			float apix_y = image->get_attr_dict().get("apix_y");
			float ds = 1.0f / (apix_y * ny * Ctf::CTFOS);
			filter = 1.0f / (filter * ds);
		}

		if (alg_name == "CtfCWauto") {
			int nxy2 = nx * ny/2;

			snri = new float[ny / 2];
			snrn = new float[ny / 2];
			tdr = new float[nxy2];
			tdi = new float[nxy2];
			tn = new float[nxy2];

			for (int i = 0; i < ny / 2; i++) {
				snri[i] = 0;
				snrn[i] = 0;
			}

			for (int i = 0; i < nxy2; i++) {
				tdr[i] = 1;
				tdi[i] = 1;
				tn[i] = 1;
			}
		}

		image0_copy = image0_fft->copy_head();
		image0_copy->ap2ri();
		image0_copy->to_zero();
	}

	Ctf::CtfType curve_type = Ctf::CTF_AMP;
	if (alg_name == "CtfCWauto") {
		curve_type = Ctf::CTF_AMP;
	}

	float *src = image->get_data();
	image->ap2ri();
	Ctf *image_ctf = image->get_ctf();
	int ny2 = image->get_ysize();

	vector<float> ctf1 = image_ctf->compute_1d(ny2,1.0f/(image_ctf->apix*image->get_ysize()), curve_type);

	if (ctf1.size() == 0) {
		LOGERR("Unexpected CTF correction problem");
	}

	ctf.push_back(ctf1);

	vector<float> ctfn1;
	if (sf) {
		ctfn1 = image_ctf->compute_1d(ny2,1.0f/(image_ctf->apix*image->get_ysize()), Ctf::CTF_SNR, sf);
	}
	else {
		ctfn1 = image_ctf->compute_1d(ny2,1.0f/(image_ctf->apix*image->get_ysize()), Ctf::CTF_SNR);
	}

	ctfn.push_back(ctfn1);

	if (alg_name == "CtfCWauto") {
		int j = 0;
		for (int y = 0; y < ny; y++) {
			for (int x = 0; x < nx / 2; x++, j += 2) {
#ifdef	_WIN32
				float r = (float)_hypot((float)x, (float)(y - ny / 2.0f));
#else
				float r = (float)hypot((float)x, (float)(y - ny / 2.0f));
#endif	//_WIN32
				int l = static_cast < int >(Util::fast_floor(r));

				if (l >= 0 && l < ny / 2) {
					int k = y*nx/2 + x;
					tdr[k] *= src[j];
					tdi[k] *= src[j + 1];
#ifdef	_WIN32
					tn[k] *= (float)_hypot(src[j], src[j + 1]);
#else
					tn[k] *= (float)hypot(src[j], src[j + 1]);
#endif	//_WIN32
				}
			}
		}
	}


	float *tmp_data = image0_copy->get_data();

	int j = 0;
	for (int y = 0; y < ny; y++) {
		for (int x = 0; x < nx / 2; x++, j += 2) {
			float r = Ctf::CTFOS * sqrt(x * x + (y - ny / 2.0f) * (y - ny / 2.0f));
			int l = static_cast < int >(Util::fast_floor(r));
			r -= l;

			float f = 0;
			if (l <= Ctf::CTFOS * ny / 2 - 2) {
				f = (ctf1[l] * (1 - r) + ctf1[l + 1] * r);
			}
			tmp_data[j] += src[j] * f;
			tmp_data[j + 1] += src[j + 1] * f;
		}
	}

	EMData *image_fft = image->do_fft();
	image_fft->update();
	if(image_ctf) {delete image_ctf; image_ctf=0;}
}

EMData * CtfAverager::finish()
{
	int j = 0;
	for (int y = 0; y < ny; y++) {
		for (int x = 0; x < nx / 2; x++, j += 2) {
#ifdef	_WIN32
			float r = (float) _hypot(x, y - ny / 2.0f);
#else
			float r = (float) hypot(x, y - ny / 2.0f);
#endif
			int l = static_cast < int >(Util::fast_floor(r));
			if (l >= 0 && l < ny / 2) {
				int k = y*nx/2 + x;
				snri[l] += (tdr[k] + tdi[k]/tn[k]);
				snrn[l] += 1;
			}
		}
	}

	for (int i = 0; i < ny / 2; i++) {
		snri[i] *= nimg / snrn[i];
	}

	if(strcmp(outfile, "") != 0) {
		Util::save_data(0, 1, snri, ny / 2, outfile);
	}


	float *cd = 0;
	if (curves) {
		cd = curves->get_data();
	}

	for (int i = 0; i < Ctf::CTFOS * ny / 2; i++) {
		float ctf0 = 0;
		for (int j = 0; j < nimg; j++) {
			ctf0 += ctfn[j][i];
			if (ctf[j][i] == 0) {
				ctf[j][i] = 1.0e-12f;
			}

			if (curves) {
				cd[i] += ctf[j][i] * ctfn[j][i];
				cd[i + Ctf::CTFOS * ny / 2] += ctfn[j][i];
				cd[i + 2 * Ctf::CTFOS * ny / 2] += ctfn[j][i];
			}
		}

		string alg_name = get_name();

		if (alg_name == "CtfCW" && need_snr) {
			snr[i] = ctf0;
		}

		float ctf1 = ctf0;
		if (alg_name == "CtfCWauto") {
			ctf1 = snri[i / Ctf::CTFOS];
		}

		if (ctf1 <= 0.0001f) {
			ctf1 = 0.1f;
		}

		if (alg_name == "CtfC") {
			for (int j = 0; j < nimg; j++) {
				ctf[j][i] = exp(-i * i / (filter * filter)) * ctfn[j][i] / (fabs(ctf[j][i]) * ctf1);
			}
		}
		else if (alg_name == "Weighting") {
			for (int j = 0; j < nimg; j++) {
				ctf[j][i] = ctfn[j][i] / ctf1;
			}
		}
		else if (alg_name == "CtfCW") {
			for (int j = 0; j < nimg; j++) {
				ctf[j][i] = (ctf1 / (ctf1 + 1)) * ctfn[j][i] / (ctf[j][i] * ctf1);
			}
		}
		else if (alg_name == "CtfCWauto") {
			for (int j = 0; j < nimg; j++) {
				ctf[j][i] = ctf1 * ctfn[j][i] / (fabs(ctf[j][i]) * ctf0);
			}
		}
	}


	if (curves) {
		for (int i = 0; i < Ctf::CTFOS * ny / 2; i++) {
			cd[i] /= cd[i + Ctf::CTFOS * ny / 2];
		}
		curves->update();
	}

	image0_copy->update();

	float *result_data = result->get_data();
	EMData *tmp_ift = image0_copy->do_ift();
	float *tmp_ift_data = tmp_ift->get_data();
	memcpy(result_data, tmp_ift_data, (nx - 2) * ny * sizeof(float));

	tmp_ift->update();
	result->update();

	if( image0_copy )
	{
		delete image0_copy;
		image0_copy = 0;
	}

	if (snri) {
		delete[]snri;
		snri = 0;
	}

	if (snrn) {
		delete[]snrn;
		snrn = 0;
	}

	if( snri )
	{
		delete [] snri;
		snri = 0;
	}
	if( snrn )
	{
		delete [] snrn;
		snrn = 0;
	}
	if( tdr )
	{
		delete [] tdr;
		tdr = 0;
	}
	if( tdi )
	{
		delete [] tdi;
		tdi = 0;
	}
	if( tn )
	{
		delete [] tn;
		tn = 0;
	}

	return result;
}

#if 0
EMData *CtfAverager::average(const vector < EMData * >&image_list) const
{
	if (image_list.size() == 0) {
		return 0;
	}

	EMData *image0 = image_list[0];
	if (image0->get_zsize() != 1) {
		LOGERR("Only 2D images are currently supported");
		return 0;
	}

	string alg_name = get_name();

	if (alg_name == "CtfCW" || alg_name == "CtfCWauto") {
		if (image0->get_ctf() != 0 && !image0->has_ctff()) {
			LOGERR("Attempted CTF Correction with no ctf parameters");
			return 0;
		}
	}
	else {
		if (image0->get_ctf() != 0) {
			LOGERR("Attempted CTF Correction with no ctf parameters");
			return 0;
		}
	}

	size_t num_images = image_list.size();
	vector < float >*ctf = new vector < float >[num_images];
	vector < float >*ctfn = new vector < float >[num_images];
	float **src = new float *[num_images];

	Ctf::CtfType curve_type = Ctf::CTF_ABS_AMP;
	if (alg_name == "CtfCWauto") {
		curve_type = Ctf::CTF_AMP;
	}

	for (size_t i = 0; i < num_images; i++) {
		EMData *image = image_list[i]->do_fft();
		image->ap2ri();
		src[i] = image->get_data();
		Ctf *image_ctf = image->get_ctf();
		int ny = image->get_ysize();
		ctf[i] = image_ctf->compute_1d(ny, curve_type);

		if (ctf[i].size() == 0) {
			LOGERR("Unexpected CTF correction problem");
			return 0;
		}

		if (sf) {
			ctfn[i] = image_ctf->compute_1d(ny, Ctf::CTF_ABS_SNR, sf);
		}
		else {
			ctfn[i] = image_ctf->compute_1d(ny, Ctf::CTF_RELATIVE_SNR);
		}

		if(image_ctf) {delete image_ctf; image_ctf=0;}
	}

	EMData *image0_fft = image0->do_fft();

	int nx = image0_fft->get_xsize();
	int ny = image0_fft->get_ysize();
	int nz = image0_fft->get_zsize();

	EMData *result = new EMData();
	result->set_size(nx - 2, ny, nz);

	float *cd = 0;
	if (alg_name == "Weighting" && curves) {
		if (!sf) {
			LOGWARN("CTF curve in file will contain relative, not absolute SNR!");
		}
		curves->set_size(Ctf::CTFOS * ny / 2, 3, 1);
		curves->to_zero();
		cd = curves->get_data();
	}

	float filter = 0;
	if (alg_name == "CtfC") {
		filter = params["filter"];
		if (filter == 0) {
			filter = 22.0f;
		}
		float apix_y = image0->get_attr_dict().get("apix_y");
		float ds = 1.0f / (apix_y * ny * Ctf::CTFOS);
		filter = 1.0f / (filter * ds);
	}

	float *snri = 0;
	float *snrn = 0;

	if (alg_name == "CtfCWauto") {
		snri = new float[ny / 2];
		snrn = new float[ny / 2];

		for (int i = 0; i < ny / 2; i++) {
			snri[i] = 0;
			snrn[i] = 0;
		}

		int j = 0;
		for (int y = 0; y < ny; y++) {
			for (int x = 0; x < nx / 2; x++, j += 2) {
				float r = hypot(x, y - ny / 2.0f);
				int l = static_cast < int >(Util::fast_floor(r));

				if (l >= 0 && l < ny / 2) {
					float tdr = 1;
					float tdi = 1;
					float tn = 1;

					for (size_t i = 0; i < num_images; i++) {
						tdr *= src[i][j];
						tdi *= src[i][j + 1];
						tn *= hypot(src[i][j], src[i][j + 1]);
					}

					tdr += tdi / tn;
					snri[l] += tdr;
					snrn[l] += 1;
				}
			}
		}

		for (int i = 0; i < ny / 2; i++) {
			snri[i] *= num_images / snrn[i];
		}
		if (outfile != "") {
			Util::save_data(0, 1, snri, ny / 2, outfile);
		}
	}

	for (int i = 0; i < Ctf::CTFOS * ny / 2; i++) {
		float ctf0 = 0;
		for (size_t j = 0; j < num_images; j++) {
			ctf0 += ctfn[j][i];
			if (ctf[j][i] == 0) {
				ctf[j][i] = 1.0e-12;
			}

			if (curves) {
				cd[i] += ctf[j][i] * ctfn[j][i];
				cd[i + Ctf::CTFOS * ny / 2] += ctfn[j][i];
				cd[i + 2 * Ctf::CTFOS * ny / 2] += ctfn[j][i];
			}
		}

		if (alg_name == "CtfCW" && need_snr) {
			snr[i] = ctf0;
		}

		float ctf1 = ctf0;
		if (alg_name == "CtfCWauto") {
			ctf1 = snri[i / Ctf::CTFOS];
		}

		if (ctf1 <= 0.0001f) {
			ctf1 = 0.1f;
		}

		if (alg_name == "CtfC") {
			for (size_t j = 0; j < num_images; j++) {
				ctf[j][i] = exp(-i * i / (filter * filter)) * ctfn[j][i] / (fabs(ctf[j][i]) * ctf1);
			}
		}
		else if (alg_name == "Weighting") {
			for (size_t j = 0; j < num_images; j++) {
				ctf[j][i] = ctfn[j][i] / ctf1;
			}
		}
		else if (alg_name == "CtfCW") {
			for (size_t j = 0; j < num_images; j++) {
				ctf[j][i] = (ctf1 / (ctf1 + 1)) * ctfn[j][i] / (ctf[j][i] * ctf1);
			}
		}
		else if (alg_name == "CtfCWauto") {
			for (size_t j = 0; j < num_images; j++) {
				ctf[j][i] = ctf1 * ctfn[j][i] / (fabs(ctf[j][i]) * ctf0);
			}
		}
	}


	if (curves) {
		for (int i = 0; i < Ctf::CTFOS * ny / 2; i++) {
			cd[i] /= cd[i + Ctf::CTFOS * ny / 2];
		}
		curves->update();
	}

	EMData *image0_copy = image0_fft->copy_head();
	image0_copy->ap2ri();

	float *tmp_data = image0_copy->get_data();

	int j = 0;
	for (int y = 0; y < ny; y++) {
		for (int x = 0; x < nx / 2; x++, j += 2) {
			float r = Ctf::CTFOS * sqrt(x * x + (y - ny / 2.0f) * (y - ny / 2.0f));
			int l = static_cast < int >(Util::fast_floor(r));
			r -= l;

			tmp_data[j] = 0;
			tmp_data[j + 1] = 0;

			for (size_t i = 0; i < num_images; i++) {
				float f = 0;
				if (l <= Ctf::CTFOS * ny / 2 - 2) {
					f = (ctf[i][l] * (1 - r) + ctf[i][l + 1] * r);
				}
				tmp_data[j] += src[i][j] * f;
				tmp_data[j + 1] += src[i][j + 1] * f;
			}
		}
	}

	image0_copy->update();

	float *result_data = result->get_data();
	EMData *tmp_ift = image0_copy->do_ift();
	float *tmp_ift_data = tmp_ift->get_data();
	memcpy(result_data, tmp_ift_data, (nx - 2) * ny * sizeof(float));

	tmp_ift->update();

	if( image0_copy )
	{
		delete image0_copy;
		image0_copy = 0;
	}

	for (size_t i = 0; i < num_images; i++) {
		EMData *img = image_list[i]->do_fft();
		img->update();
	}

	if( src )
	{
		delete[]src;
		src = 0;
	}

	if( ctf )
	{
		delete[]ctf;
		ctf = 0;
	}

	if( ctfn )
	{
		delete[]ctfn;
		ctfn = 0;
	}

	if (snri) {
		delete[]snri;
		snri = 0;
	}

	if (snrn) {
		delete[]snrn;
		snrn = 0;
	}

	result->update();
	return result;
}
#endif


void EMAN::dump_averagers()
{
	dump_factory < Averager > ();
}

map<string, vector<string> > EMAN::dump_averagers_list()
{
	return dump_factory_list < Averager > ();
}
